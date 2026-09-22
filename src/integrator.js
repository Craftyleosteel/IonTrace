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

/** Acceleration a = (q/m) E at the ion's current position. */
export function accelerationAt(field, qOverM, x, z) {
  const { Ex, Ez } = field.fieldAtCartesian(x, z);
  return { ax: qOverM * Ex, az: qOverM * Ez };
}

/**
 * One classical fourth-order Runge-Kutta step.
 *
 * The state is (x, z, vx, vz) and the derivative is
 * (vx, vz, q Ex / m, q Ez / m). Because the acceleration depends only on
 * position, each stage needs exactly one field evaluation.
 */
export function stepRK4(field, ion, dt) {
  const qm = ion.charge / ion.mass;
  const { x, z, vx, vz } = ion;

  const a1 = accelerationAt(field, qm, x, z);
  const k1 = { dx: vx, dz: vz, dvx: a1.ax, dvz: a1.az };

  const a2 = accelerationAt(field, qm, x + 0.5 * dt * k1.dx, z + 0.5 * dt * k1.dz);
  const k2 = {
    dx: vx + 0.5 * dt * k1.dvx,
    dz: vz + 0.5 * dt * k1.dvz,
    dvx: a2.ax,
    dvz: a2.az,
  };

  const a3 = accelerationAt(field, qm, x + 0.5 * dt * k2.dx, z + 0.5 * dt * k2.dz);
  const k3 = {
    dx: vx + 0.5 * dt * k2.dvx,
    dz: vz + 0.5 * dt * k2.dvz,
    dvx: a3.ax,
    dvz: a3.az,
  };

  const a4 = accelerationAt(field, qm, x + dt * k3.dx, z + dt * k3.dz);
  const k4 = {
    dx: vx + dt * k3.dvx,
    dz: vz + dt * k3.dvz,
    dvx: a4.ax,
    dvz: a4.az,
  };

  const sixth = dt / 6;
  return {
    ...ion,
    x: x + sixth * (k1.dx + 2 * k2.dx + 2 * k3.dx + k4.dx),
    z: z + sixth * (k1.dz + 2 * k2.dz + 2 * k3.dz + k4.dz),
    vx: vx + sixth * (k1.dvx + 2 * k2.dvx + 2 * k3.dvx + k4.dvx),
    vz: vz + sixth * (k1.dvz + 2 * k2.dvz + 2 * k3.dvz + k4.dvz),
    t: ion.t + dt,
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
export function stepVerlet(field, ion, dt) {
  const qm = ion.charge / ion.mass;
  const a0 = accelerationAt(field, qm, ion.x, ion.z);

  const x = ion.x + ion.vx * dt + 0.5 * a0.ax * dt * dt;
  const z = ion.z + ion.vz * dt + 0.5 * a0.az * dt * dt;

  const a1 = accelerationAt(field, qm, x, z);

  return {
    ...ion,
    x,
    z,
    vx: ion.vx + 0.5 * (a0.ax + a1.ax) * dt,
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
  return 0.5 * ion.mass * (ion.vx * ion.vx + ion.vz * ion.vz);
}

/**
 * Total energy in joules: kinetic plus electrostatic potential energy qV.
 *
 * In a static field this is a constant of the motion, so tracking it is the
 * single most informative diagnostic available - it catches integrator error,
 * a mis-scaled field, a unit slip and a bad interpolation all at once.
 */
export function totalEnergy(field, ion) {
  return kineticEnergy(ion) + ion.charge * field.potentialAtCartesian(ion.x, ion.z);
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
  const h = field.grid.step;
  const qm = ion.charge / ion.mass;
  const { ax, az } = accelerationAt(field, qm, ion.x, ion.z);

  const speed = Math.hypot(ion.vx, ion.vz);
  const accel = Math.hypot(ax, az);

  const byTravel = speed > 0 ? (cfl * h) / speed : Infinity;
  const byAccel = accel > 0 ? Math.sqrt((2 * cfl * h) / accel) : Infinity;
  const dt = Math.min(byTravel, byAccel);

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

  const { grid } = field;
  const zMin = grid.z0;
  const zMax = grid.z0 + grid.zLength;
  const rMax = grid.rLength;

  // Cylindrical geometry folds the signed transverse coordinate onto a
  // radius; planar geometry does not, because there its transverse axis runs
  // from 0 to rLength with no symmetry about zero. Folding it anyway would
  // mirror the electrode map about y = 0 and let an ion below the floor keep
  // flying through an extrapolated field.
  const planar = grid.symmetry === PLANAR;
  const transverse = (x) => (planar ? x : Math.abs(x));

  let current = { ...ion };
  const points = [current];

  const E0 = totalEnergy(field, current);
  const scale = Math.abs(E0) > 0 ? Math.abs(E0) : kineticEnergy(current) || 1;
  let worstDrift = 0;

  let stop = 'step-limit';

  for (let n = 1; n <= maxSteps; n++) {
    const dt = suggestTimeStep(field, current, cfl);
    const next = step(field, current, dt);

    const r = transverse(next.x);

    // Metal is tested first. The outer radial wall is real hardware, so a
    // step that overshoots it must be a strike and not an escape - and since
    // electrodeHit already ignores the open end faces, testing it first
    // cannot steal a legitimate exit.
    //
    // Resolution is one grid node, which biases every aperture *inward* by
    // h/2 rather than merely blurring it: transmission is systematically
    // pessimistic. See docs/PHYSICS.md section 3.6.
    if (r > rMax || r < 0 || electrodeHit(grid, next.z, r)) {
      current = next;
      stop = 'electrode';
      points.push(current);
      break;
    }

    // Leaving through an open face ends the flight, and the direction matters.
    // An ion that comes back out of the entrance has been REFLECTED, not
    // transmitted; reporting both as "exited" turns a working ion mirror into
    // a lens with a negative focal length.
    if (next.z > zMax) {
      current = next;
      stop = 'exited';
      points.push(current);
      break;
    }
    if (next.z < zMin) {
      current = next;
      stop = 'reflected';
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
