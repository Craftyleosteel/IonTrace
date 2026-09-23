/**
 * Trajectory integration - Newton's second law for an ion in a static field.
 *
 * The equation of motion, in full, for what IonTrace currently models:
 *
 *     m d2r/dt2 = q E(r)
 *
 * That is the whole of it. There is no magnetic term, no gas collision term,
 * no ion-ion Coulomb term and no relativistic correction. Every one of those
 * omissions is deliberate for this build and is listed in docs/PHYSICS.md;
 * the point of keeping the right-hand side this small is that each term can
 * be validated on its own before another is added.
 *
 * Two integrators are provided deliberately.
 *
 *   RK4      Fourth-order Runge-Kutta. Highest accuracy per step, and the
 *            method SIMION uses. Not symplectic, so its energy error grows
 *            slowly and secularly over very long flights.
 *
 *   Verlet   Velocity Verlet. Only second order, but symplectic: its energy
 *            error oscillates about zero instead of drifting, which is the
 *            behaviour that matters once periodic (trap) fields arrive.
 *
 * Running the same case through both and comparing is a genuine check on the
 * result rather than on one method's self-consistency, so the test suite does
 * exactly that.
 */

import { PLANAR } from './grid.js';
import {
  currentShares,
  spaceChargeField,
  lineChargeDensity,
  coulombField,
} from './spacecharge.js';

/**
 * @typedef {object} IonState
 * @property {number} mass    kg
 * @property {number} charge  C (signed)
 * @property {number} x       signed transverse position, m
 * @property {number} z       axial position, m
 * @property {number} vx      transverse velocity, m/s
 * @property {number} vz      axial velocity, m/s
 * @property {number} t       elapsed time, s
 */

/**
 * Acceleration a = (q/m) E at the ion's current position.
 *
 * `extra` adds a field that does not come from the electrode solve - at
 * present the beam's own space charge. It is supplied by the caller as a
 * constant vector for the whole step rather than re-evaluated at each RK4
 * stage, because the self-field depends on where every OTHER ion is, and
 * those have no defined position at an intermediate stage. Freezing the
 * self-force across a step is the standard particle-in-cell treatment. It
 * costs accuracy: the space-charge part of the motion is effectively second
 * order even though the electrode part stays fourth.
 */
export function accelerationAt(field, qOverM, x, y, z, t, extra) {
  const E = field.fieldAt3D(x, y, z, t);
  if (extra) {
    return {
      ax: qOverM * (E.Ex + extra.Ex),
      ay: qOverM * (E.Ey + (extra.Ey ?? 0)),
      az: qOverM * (E.Ez + extra.Ez),
    };
  }
  return { ax: qOverM * E.Ex, ay: qOverM * E.Ey, az: qOverM * E.Ez };
}

/**
 * One classical fourth-order Runge-Kutta step.
 *
 * The state is (x, y, z, vx, vy, vz) and the derivative is
 * (vx, vy, vz, q Ex / m, q Ey / m, q Ez / m). Each stage needs exactly one
 * field evaluation.
 *
 * Stages are evaluated at t, t + dt/2, t + dt/2 and t + dt. For a static field
 * that is redundant; for a radio-frequency one it is the difference between
 * fourth order and second.
 */
export function stepRK4(field, ion, dt, extra) {
  const qm = ion.charge / ion.mass;
  const { x, z, vx, vz, t } = ion;
  const y = ion.y ?? 0;
  const vy = ion.vy ?? 0;
  const h = 0.5 * dt;

  const a1 = accelerationAt(field, qm, x, y, z, t, extra);
  const k1 = { dx: vx, dy: vy, dz: vz, dvx: a1.ax, dvy: a1.ay, dvz: a1.az };

  const a2 = accelerationAt(
    field, qm, x + h * k1.dx, y + h * k1.dy, z + h * k1.dz, t + h, extra
  );
  const k2 = {
    dx: vx + h * k1.dvx,
    dy: vy + h * k1.dvy,
    dz: vz + h * k1.dvz,
    dvx: a2.ax, dvy: a2.ay, dvz: a2.az,
  };

  const a3 = accelerationAt(
    field, qm, x + h * k2.dx, y + h * k2.dy, z + h * k2.dz, t + h, extra
  );
  const k3 = {
    dx: vx + h * k2.dvx,
    dy: vy + h * k2.dvy,
    dz: vz + h * k2.dvz,
    dvx: a3.ax, dvy: a3.ay, dvz: a3.az,
  };

  const a4 = accelerationAt(
    field, qm, x + dt * k3.dx, y + dt * k3.dy, z + dt * k3.dz, t + dt, extra
  );
  const k4 = {
    dx: vx + dt * k3.dvx,
    dy: vy + dt * k3.dvy,
    dz: vz + dt * k3.dvz,
    dvx: a4.ax, dvy: a4.ay, dvz: a4.az,
  };

  const sixth = dt / 6;
  const w = (p, q, r, s) => sixth * (p + 2 * q + 2 * r + s);

  return {
    ...ion,
    x: x + w(k1.dx, k2.dx, k3.dx, k4.dx),
    y: y + w(k1.dy, k2.dy, k3.dy, k4.dy),
    z: z + w(k1.dz, k2.dz, k3.dz, k4.dz),
    vx: vx + w(k1.dvx, k2.dvx, k3.dvx, k4.dvx),
    vy: vy + w(k1.dvy, k2.dvy, k3.dvy, k4.dvy),
    vz: vz + w(k1.dvz, k2.dvz, k3.dvz, k4.dvz),
    t: t + dt,
  };
}

/**
 * One velocity-Verlet step.
 *
 *     r(t+dt) = r(t) + v(t) dt + 1/2 a(t) dt^2
 *     v(t+dt) = v(t) + 1/2 [ a(t) + a(t+dt) ] dt
 *
 * Costs one field evaluation per step (the new acceleration is reused as the
 * next step's old one only if the caller loops; here it is recomputed, which
 * is simpler and still cheap).
 */
export function stepVerlet(field, ion, dt, extra) {
  const qm = ion.charge / ion.mass;
  const y0 = ion.y ?? 0;
  const vy0 = ion.vy ?? 0;
  const half = 0.5 * dt * dt;

  const a0 = accelerationAt(field, qm, ion.x, y0, ion.z, ion.t, extra);

  const x = ion.x + ion.vx * dt + half * a0.ax;
  const y = y0 + vy0 * dt + half * a0.ay;
  const z = ion.z + ion.vz * dt + half * a0.az;

  const a1 = accelerationAt(field, qm, x, y, z, ion.t + dt, extra);

  return {
    ...ion,
    x,
    y,
    z,
    vx: ion.vx + 0.5 * (a0.ax + a1.ax) * dt,
    vy: vy0 + 0.5 * (a0.ay + a1.ay) * dt,
    vz: ion.vz + 0.5 * (a0.az + a1.az) * dt,
    t: ion.t + dt,
  };
}

export const INTEGRATORS = { rk4: stepRK4, verlet: stepVerlet };

/**
 * Default step-size aggressiveness: the fraction of a grid step an ion is
 * allowed to advance per step. Defined once so the library default, the JSDoc
 * and docs/PHYSICS.md cannot drift apart from each other.
 */
export const DEFAULT_CFL = 0.05;

/** Kinetic energy in joules. */
export function kineticEnergy(ion) {
  const vy = ion.vy ?? 0;
  return 0.5 * ion.mass * (ion.vx * ion.vx + vy * vy + ion.vz * ion.vz);
}

/**
 * Total energy in joules: kinetic plus electrostatic potential energy qV.
 *
 * In a static field this is a constant of the motion, so tracking it is the
 * single most informative diagnostic available - it catches integrator error,
 * a mis-scaled field, a unit slip and a bad interpolation all at once.
 */
export function totalEnergy(field, ion) {
  return (
    kineticEnergy(ion) +
    ion.charge * field.potentialAt3D(ion.x, ion.y ?? 0, ion.z, ion.t)
  );
}

/**
 * Choose a step size for the ion's current state.
 *
 * Two limits are imposed, and the smaller wins:
 *
 *   travel        |v| dt  <= cfl * h, so the ion cannot skip over grid cells
 *                 and miss field structure between them.
 *   acceleration  1/2 |a| dt^2 <= cfl * h, which takes over near rest, where
 *                 the travel limit alone would permit an unbounded step.
 *
 * `cfl` below 1 keeps the ion sampling the field several times per cell.
 */
export function suggestTimeStep(field, ion, cfl = DEFAULT_CFL) {
  const h = field.lengthScale;
  const qm = ion.charge / ion.mass;
  const y = ion.y ?? 0;
  const vy = ion.vy ?? 0;
  const { ax, ay, az } = accelerationAt(field, qm, ion.x, y, ion.z, ion.t);

  const speed = Math.hypot(ion.vx, vy, ion.vz);
  const accel = Math.hypot(ax, ay, az);

  const byTravel = speed > 0 ? (cfl * h) / speed : Infinity;
  const byAccel = accel > 0 ? Math.sqrt((2 * cfl * h) / accel) : Infinity;
  // A time-dependent field imposes its own limit: the step must resolve the
  // fastest oscillation present, or the ion samples an aliased field rather
  // than the real one. Without this an RF quadrupole would look stable at any
  // amplitude simply because the step skipped over the field's reversals.
  //
  // At the default cfl this is twenty steps per RF period. Resolving the
  // oscillation is not the same as resolving the grid: an ion can be moving
  // slowly enough that the travel limit permits a step spanning several
  // reversals, and that step would integrate a field the ion never saw.
  const byPeriod = field.shortestPeriod ? cfl * field.shortestPeriod : Infinity;
  const dt = Math.min(byTravel, byAccel, byPeriod);

  // A completely motionless ion in a null field has no natural scale; fall
  // back to something finite so the caller's loop still terminates.
  return Number.isFinite(dt) ? dt : 1e-9;
}

/**
 * Fly one ion until it leaves the domain, strikes an electrode, or runs out
 * of budget.
 *
 * @param {import('./field.js').Field} field
 * @param {IonState} ion  Initial state, in SI units.
 * @param {object} [opts]
 * @param {'rk4'|'verlet'} [opts.method]
 * @param {number} [opts.cfl]        Step-size aggressiveness (default DEFAULT_CFL).
 * @param {number} [opts.maxSteps]   Iteration cap (default 200000).
 * @param {number} [opts.maxTime]    Flight-time cap in seconds.
 * @param {number} [opts.recordEvery] Keep every Nth point (default 1).
 * @returns {{points: IonState[], stop: StopReason, energyDrift: number}}
 *
 * `stop` is one of:
 *   'exited'     left forwards through the far face - transmitted
 *   'reflected'  came back out of the entrance face - NOT transmitted
 *   'electrode'  struck metal, including the outer wall
 *   'time-limit' / 'step-limit'  ran out of budget
 */
export function flyIon(field, ion, opts = {}) {
  const step = INTEGRATORS[opts.method ?? 'rk4'];
  if (!step) throw new Error(`Unknown integrator "${opts.method}"`);

  const cfl = opts.cfl ?? DEFAULT_CFL;
  const maxSteps = opts.maxSteps ?? 200000;
  const maxTime = opts.maxTime ?? Infinity;
  const recordEvery = opts.recordEvery ?? 1;

  let current = { ...ion };
  const points = [current];

  const E0 = totalEnergy(field, current);
  const scale = Math.abs(E0) > 0 ? Math.abs(E0) : kineticEnergy(current) || 1;
  let worstDrift = 0;

  let stop = 'step-limit';

  for (let n = 1; n <= maxSteps; n++) {
    const dt = suggestTimeStep(field, current, cfl);
    const next = step(field, current, dt);

    // Metal is tested first. The outer wall is real hardware, so a step that
    // overshoots it must be a strike and not an escape - and since `strikes`
    // already ignores the open end faces, testing it first cannot steal a
    // legitimate exit.
    //
    // Resolution is one grid node, which biases every aperture *inward* by
    // h/2 rather than merely blurring it: transmission is systematically
    // pessimistic. See docs/PHYSICS.md section 3.6.
    if (field.strikes(next.x, next.y ?? 0, next.z)) {
      current = next;
      stop = 'electrode';
      points.push(current);
      break;
    }

    // Leaving ends the flight, and the direction matters. An ion that comes
    // back out of the entrance has been REFLECTED, not transmitted; reporting
    // both as "exited" turns a working ion mirror into a lens with a negative
    // focal length. A third case exists once elements can be misaligned: an
    // ion can escape through a gap without reaching either end, and calling
    // that "exited" would count a lost ion as transmitted.
    const where = field.classify(next.x, next.y ?? 0, next.z);
    if (where !== 'inside') {
      current = next;
      stop = where;
      points.push(current);
      break;
    }

    current = next;

    const drift = Math.abs(totalEnergy(field, current) - E0) / scale;
    if (drift > worstDrift) worstDrift = drift;

    if (n % recordEvery === 0) points.push(current);

    if (current.t >= maxTime) {
      stop = 'time-limit';
      if (points[points.length - 1] !== current) points.push(current);
      break;
    }
  }

  if (points[points.length - 1] !== current) points.push(current);

  return { points, stop, energyDrift: worstDrift };
}

/**
 * Create a beam flight that can be advanced a chunk at a time.
 *
 * Unlike `flyIon`, which takes one ion to its end before starting the next,
 * this advances every ion together on a SHARED time step. It has to: the
 * force on any ion depends on where all the others are at that instant, so a
 * beam whose members sit at different times has no defined mutual force. Ions
 * that terminate stop contributing, which is correct - an ion that has struck
 * metal or left the modelled region is no longer part of the beam.
 *
 * Two repulsion models, because they answer different questions:
 *
 *   'beam'     Each trajectory is a RING of charge and the field follows from
 *              Gauss's law on the enclosed current. Correct for a continuous
 *              beam; driven by `beamCurrent`.
 *   'coulomb'  Each trajectory is one point charge feeling every other one
 *              directly. Correct for a countable bunch or cloud; driven by
 *              `ionsPerParticle`, the number of real ions each simulated
 *              particle stands for.
 *   'none'     No self-interaction.
 *
 * With repulsion off, trajectories are bit-identical to `flyIon` - verified by
 * test, since a self-field model that perturbs the answer when switched off is
 * worse than none.
 *
 * @param {import('./field.js').Field} field
 * @param {IonState[]} ions
 * @param {object} [opts]
 * @param {'none'|'beam'|'coulomb'} [opts.repulsion]
 * @param {number} [opts.beamCurrent]     Beam current in amperes ('beam').
 * @param {number} [opts.ionsPerParticle] Macro-weight ('coulomb').
 * @param {number} [opts.softening]       Plummer softening length, metres.
 * @param {'rk4'|'verlet'} [opts.method]
 * @param {number} [opts.cfl]
 * @param {number} [opts.maxSteps]
 * @param {number} [opts.recordEvery]
 * @returns {{tracks: object[], steps: number, done: boolean,
 *            advance: (budget?: number) => object}}
 */
export function createFlight(field, ions, opts = {}) {
  const step = INTEGRATORS[opts.method ?? 'rk4'];
  if (!step) throw new Error(`Unknown integrator "${opts.method}"`);

  const cfl = opts.cfl ?? DEFAULT_CFL;
  const maxSteps = opts.maxSteps ?? 200000;
  const recordEvery = opts.recordEvery ?? 1;
  const repulsion = opts.repulsion ?? (opts.beamCurrent ? 'beam' : 'none');
  const beamCurrent = opts.beamCurrent ?? 0;
  const ionsPerParticle = opts.ionsPerParticle ?? 1;

  const planar = field.grid?.symmetry === PLANAR;

  // Ring model only: each ion's share of the beam current is fixed at launch
  // from its starting radius and conserved thereafter.
  const shares = currentShares(ions.map((ion) => Math.abs(ion.x)));

  const tracks = ions.map((ion, index) => {
    const state = { ...ion };
    const E0 = totalEnergy(field, state);
    return {
      state,
      points: [state],
      stop: 'step-limit',
      active: true,
      index,
      E0,
      scale: Math.abs(E0) > 0 ? Math.abs(E0) : kineticEnergy(state) || 1,
      energyDrift: 0,
    };
  });

  // Ring model: inside one grid step of the axis the ring sampling has no
  // resolution anyway, and 1/r would otherwise be dominated by discretisation
  // noise.
  const ringSoftening = field.lengthScale / 2;
  // Discrete model: bounds the 1/r^2 force during a close pass, which a finite
  // time step would otherwise turn into energy from nowhere.
  const coulombSoftening = opts.softening ?? field.lengthScale;

  const flight = {
    tracks,
    steps: 0,
    done: false,
    advance,
  };

  /**
   * Advance the whole beam by at most `budget` steps.
   *
   * Splitting the flight into chunks is what lets the UI draw it live. It
   * changes nothing about the result: the time step is chosen from each ion's
   * own state, never from wall-clock time, so the trajectory is identical
   * whether it is run in one call or a hundred, on a fast machine or a slow
   * one.
   */
  function advance(budget = Infinity) {
    let taken = 0;
    while (taken < budget && flight.steps < maxSteps) {
      const live = tracks.filter((t) => t.active);
      if (live.length === 0) break;

      // One step for everyone, so the beam stays synchronised. It has to be:
      // the force between ions depends on where they all are at the same
      // instant, and a beam whose members sit at different times has no
      // defined mutual force at all.
      let dt = Infinity;
      for (const t of live) {
        dt = Math.min(dt, suggestTimeStep(field, t.state, cfl));
      }
      if (!Number.isFinite(dt) || dt <= 0) break;

      const extras = selfFields(live);

      for (let k = 0; k < live.length; k++) {
        const t = live[k];
        const next = step(field, t.state, dt, extras[k]);

        if (field.strikes(next.x, next.y ?? 0, next.z)) {
          t.state = next;
          t.stop = 'electrode';
          t.active = false;
        } else {
          const where = field.classify(next.x, next.y ?? 0, next.z);
          if (where !== 'inside') {
            t.state = next;
            t.stop = where;
            t.active = false;
          } else {
            t.state = next;
            const drift = Math.abs(totalEnergy(field, next) - t.E0) / t.scale;
            if (drift > t.energyDrift) t.energyDrift = drift;
          }
        }

        if (!t.active || flight.steps % recordEvery === 0) t.points.push(t.state);
      }

      flight.steps++;
      taken++;
    }

    if (flight.steps >= maxSteps || tracks.every((t) => !t.active)) finish();
    return flight;
  }

  /** Self-field on each live ion, as a constant vector for this step. */
  function selfFields(live) {
    if (repulsion === 'coulomb' && ionsPerParticle !== 0 && live.length > 1) {
      // Every particle is a point charge and feels every other one directly.
      const { Ex, Ey, Ez } = coulombField(
        live.map((t) => t.state.x),
        live.map((t) => t.state.y ?? 0),
        live.map((t) => t.state.z),
        live.map((t) => t.state.charge),
        ionsPerParticle,
        coulombSoftening
      );
      return live.map((_, k) =>
        Ex[k] === 0 && Ey[k] === 0 && Ez[k] === 0
          ? null
          : { Ex: Ex[k], Ey: Ey[k], Ez: Ez[k] }
      );
    }

    if (repulsion === 'beam' && beamCurrent !== 0) {
      const radii = live.map((t) => Math.hypot(t.state.x, t.state.y ?? 0));
      const liveShares = live.map((t) => shares[t.index]);
      // The density is set by how fast the beam is moving HERE: a decelerated
      // beam is a denser one, which is why space charge bites hardest in the
      // slow part of an optic.
      const meanVz =
        live.reduce((s, t) => s + Math.abs(t.state.vz), 0) / live.length;
      const lambda = lineChargeDensity(beamCurrent, meanVz);
      const sign = Math.sign(live[0].state.charge) || 1;
      const Er = spaceChargeField(radii, liveShares, lambda * sign, ringSoftening);
      return live.map((t, k) => {
        if (Er[k] === 0) return null;
        const x = t.state.x;
        const y = t.state.y ?? 0;
        if (planar) return { Ex: Er[k], Ey: 0, Ez: 0 };
        // Resolve the radial self-field onto x and y by direction cosines.
        const r = Math.hypot(x, y);
        if (r === 0) return null;
        return { Ex: (Er[k] * x) / r, Ey: (Er[k] * y) / r, Ez: 0 };
      });
    }

    return live.map(() => null);
  }

  function finish() {
    if (flight.done) return;
    for (const t of tracks) {
      if (t.points[t.points.length - 1] !== t.state) t.points.push(t.state);
      delete t.active;
      delete t.index;
      delete t.E0;
      delete t.scale;
    }
    flight.done = true;
  }

  return flight;
}

/** Run a beam flight to completion. Thin wrapper over `createFlight`. */
export function flyBeam(field, ions, opts = {}) {
  const flight = createFlight(field, ions, opts);
  while (!flight.done) flight.advance(2000);
  return { tracks: flight.tracks, steps: flight.steps };
}

/**
 * True if the grid node nearest (z, r) is metal the ion has struck.
 *
 * Nodes on an open face are excluded. Those nodes carry a fixed potential
 * because the Laplace problem needs the domain closed, but they represent the
 * entrance and exit apertures of the modelled region, not a surface. Counting
 * them as strikes would report every successfully transmitted ion as having
 * crashed into the far wall.
 */
export function electrodeHit(grid, z, r) {
  const i = Math.round((z - grid.z0) / grid.step);
  const j = Math.round(r / grid.step);
  if (i < 0 || j < 0 || i >= grid.nz || j >= grid.nr) return false;
  if (i === 0 && grid.openFaces?.zMin) return false;
  if (i === grid.nz - 1 && grid.openFaces?.zMax) return false;
  return grid.isElectrode(i, j);
}
