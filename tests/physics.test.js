/**
 * Physics validation suite.
 *
 * The organising principle: validate each layer against something external to
 * it, never against itself.
 *
 *   constants     against CODATA values and textbook numbers
 *   Laplace       against closed-form harmonic functions
 *   integrators   against an analytically solvable mock field, so integrator
 *                 error is isolated from field-solver error
 *   einzel lens   against conservation laws and symmetries that hold
 *                 regardless of how well the numerics perform
 *
 * The last group is the important one for a user: energy conservation and
 * mirror symmetry are properties the *real physics* has, so a violation is a
 * bug no matter how plausible the trajectory looks on screen.
 */

import { describe, it, assert, assertClose, assertRelClose } from './harness.js';

import {
  ELEMENTARY_CHARGE,
  ATOMIC_MASS_UNIT,
  ELECTRON_MASS,
  speedFromKineticEnergy,
  relativisticError,
  eVToJoules,
  amuToKg,
  mmToM,
  joulesToEV,
} from '../src/constants.js';

import { PotentialArray, CYLINDRICAL, PLANAR, NO_ELECTRODE } from '../src/grid.js';
import { relax, maxResidual, solveBasis, optimalOmega } from '../src/laplace.js';
import { Field } from '../src/field.js';
import {
  stepRK4,
  stepVerlet,
  flyIon,
  totalEnergy,
  kineticEnergy,
} from '../src/integrator.js';
import { makeIon, parallelBeam, axialCrossing } from '../src/ion.js';
import { buildEinzelLens } from '../src/geometries/einzel.js';

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

/**
 * Impose a known analytic potential on the domain boundary, relax the
 * interior, and return both the numerical solution and the worst interior
 * error.
 *
 * This is the sharpest test available for a Laplace solver: if `exact` is
 * harmonic, it is the unique solution of the boundary value problem, so any
 * interior discrepancy is solver error and nothing else.
 */
function solveAgainstExact(grid, exact, opts = {}) {
  const id = grid.addElectrode('boundary');
  grid.paintEnclosure(id);

  const phi = new Float64Array(grid.nz * grid.nr);
  for (let j = 0; j < grid.nr; j++) {
    for (let i = 0; i < grid.nz; i++) {
      if (grid.isElectrode(i, j)) {
        phi[grid.idx(i, j)] = exact(grid.zAt(i), grid.rAt(j));
      }
    }
  }

  const report = relax(grid, phi, { tolerance: 1e-12, maxSweeps: 200000, ...opts });

  let worst = 0;
  for (let j = 0; j < grid.nr; j++) {
    for (let i = 0; i < grid.nz; i++) {
      if (grid.isElectrode(i, j)) continue;
      const err = Math.abs(phi[grid.idx(i, j)] - exact(grid.zAt(i), grid.rAt(j)));
      if (err > worst) worst = err;
    }
  }

  return { phi, worst, report };
}

/**
 * A mock Field with an analytically known solution: a linear restoring force
 * along x, which makes the motion simple harmonic.
 *
 *     phi(x, z) = 1/2 k x^2      =>      Ex = -k x,   Ez = 0
 *
 * For charge q and mass m this gives omega = sqrt(q k / m) and, for an ion
 * released from rest at x = A,
 *
 *     x(t) = A cos(omega t),   vx(t) = -A omega sin(omega t)
 *
 * Using this instead of a solved grid lets integrator error be measured
 * without any contamination from interpolation or discretisation, which is
 * the only way to see an integrator's true order of accuracy.
 */
function harmonicField(k) {
  return {
    grid: { step: 1e-3, symmetry: PLANAR, z0: 0, zLength: 1, rLength: 1, nz: 3, nr: 3 },
    fieldAtCartesian(x) {
      return { Ex: -k * x, Ez: 0 };
    },
    potentialAtCartesian(x) {
      return 0.5 * k * x * x;
    },
  };
}

/** Integrate `steps` fixed-size steps and return the final state. */
function integrateFixed(stepFn, field, ion, dt, steps) {
  let s = ion;
  for (let n = 0; n < steps; n++) s = stepFn(field, s, dt);
  return s;
}

/** Empirical order of accuracy from errors at step h and h/2. */
function observedOrder(errCoarse, errFine) {
  return Math.log2(errCoarse / errFine);
}

/* ------------------------------------------------------------------ */
/* constants and units                                                 */
/* ------------------------------------------------------------------ */

describe('Constants and unit conversion', () => {
  it('uses the exact SI-2019 elementary charge', () => {
    assertClose(ELEMENTARY_CHARGE, 1.602176634e-19, 0, 'e is exact by definition');
  });

  it('uses the CODATA atomic mass unit', () => {
    assertRelClose(ATOMIC_MASS_UNIT, 1.66053906892e-27, 1e-12, 'u');
  });

  it('reproduces the textbook speed of a 1 eV electron', () => {
    // v = sqrt(2E/m) = sqrt(2 x 1.602177e-19 / 9.109384e-31) ~ 5.931e5 m/s
    const v = speedFromKineticEnergy(eVToJoules(1), ELECTRON_MASS);
    assertRelClose(v, 5.9309e5, 1e-4, '1 eV electron speed');
  });

  it('reproduces the speed of a 10 eV, 100 u singly charged ion', () => {
    // v = sqrt(2 x 10 x 1.602177e-19 / (100 x 1.660539e-27)) ~ 4392.7 m/s
    const v = speedFromKineticEnergy(eVToJoules(10), amuToKg(100));
    assertRelClose(v, 4392.7, 1e-4, '10 eV, 100 u ion speed');
  });

  it('scales speed as the square root of energy', () => {
    const m = amuToKg(40);
    const v1 = speedFromKineticEnergy(eVToJoules(10), m);
    const v4 = speedFromKineticEnergy(eVToJoules(40), m);
    assertRelClose(v4 / v1, 2, 1e-12, 'quadrupling energy doubles speed');
  });

  it('round-trips energy through the eV conversion', () => {
    assertRelClose(joulesToEV(eVToJoules(137.036)), 137.036, 1e-12, 'eV round trip');
  });

  it('confirms the non-relativistic assumption holds for ions', () => {
    // A 1000 eV, 100 u ion moves at 4.4e4 m/s, so beta ~ 1.5e-4 and the
    // leading Newtonian error is (3/4) beta^2 ~ 1.6e-8. Checking against that
    // estimate rather than an arbitrary small number makes the test fail if
    // the error model itself is wrong, not just if it is large.
    const m = amuToKg(100);
    const v = speedFromKineticEnergy(eVToJoules(1000), m);
    const beta = v / 299792458;
    const err = relativisticError(v);
    assertRelClose(err, 0.75 * beta * beta, 1e-6, 'leading-order relativistic error');
    assert(err < 1e-6, `keV-scale ions must be safely Newtonian, got ${err}`);
  });

  it('flags a case where the non-relativistic assumption fails', () => {
    // A 1 MeV electron is emphatically not Newtonian; the helper must say so.
    const v = speedFromKineticEnergy(eVToJoules(1e6), ELECTRON_MASS);
    assert(
      relativisticError(v) > 0.1,
      'MeV electrons must be flagged as relativistic'
    );
  });

  it('rejects a neutral particle', () => {
    let threw = false;
    try {
      makeIon({ mass: 100, charge: 0, energy: 10 });
    } catch {
      threw = true;
    }
    assert(threw, 'a charge of zero is not a physically meaningful ion here');
  });
});

/* ------------------------------------------------------------------ */
/* Laplace solver                                                      */
/* ------------------------------------------------------------------ */

describe('Laplace solver', () => {
  it('reproduces a harmonic function exactly in planar geometry', () => {
    // phi = z^2 - y^2 satisfies d2/dz2 + d2/dy2 = 2 - 2 = 0, and the
    // five-point stencil is exact for quadratics, so the only error allowed
    // is the relaxation tolerance.
    const grid = new PotentialArray({ nz: 41, nr: 31, step: 1e-3, symmetry: PLANAR });
    const exact = (z, y) => z * z - y * y;
    const { worst } = solveAgainstExact(grid, exact);
    assert(worst < 1e-9, `planar harmonic error ${worst.toExponential(3)} too large`);
  });

  it('reproduces a harmonic function exactly in cylindrical geometry', () => {
    // phi = z^2 - r^2/2 is axisymmetric-harmonic:
    //   d2/dr2 = -1,  (1/r) d/dr = -1,  d2/dz2 = +2  ->  sum 0.
    // Crucially this exercises the on-axis stencil, where the (1/r) dphi/dr
    // term is singular and is replaced by its limit. A mis-derived axis
    // stencil fails here and nowhere else.
    const grid = new PotentialArray({ nz: 41, nr: 31, step: 1e-3, symmetry: CYLINDRICAL });
    const exact = (z, r) => z * z - 0.5 * r * r;
    const { worst } = solveAgainstExact(grid, exact);
    assert(
      worst < 1e-9,
      `cylindrical harmonic error ${worst.toExponential(3)} too large`
    );
  });

  it('keeps the on-axis solution consistent with the off-axis solution', () => {
    // Independent of the test above: the axis row must not be a special case
    // that happens to satisfy its own stencil while disagreeing with the
    // interior. Compare phi(r=0) against the exact value directly.
    const grid = new PotentialArray({ nz: 41, nr: 31, step: 1e-3, symmetry: CYLINDRICAL });
    const exact = (z, r) => z * z - 0.5 * r * r;
    const { phi } = solveAgainstExact(grid, exact);
    for (let i = 1; i < grid.nz - 1; i++) {
      assertClose(phi[grid.idx(i, 0)], exact(grid.zAt(i), 0), 1e-9, `axis node i=${i}`);
    }
  });

  it('converges at second order in the grid step', () => {
    // 1/sqrt(r^2 + z^2) is the axisymmetric point-charge potential, harmonic
    // everywhere except the origin. Keeping the origin outside the domain
    // gives a smooth non-polynomial solution that the stencil cannot
    // reproduce exactly, so a genuine discretisation error appears and its
    // scaling with h can be measured.
    const exact = (z, r) => 1 / Math.sqrt(r * r + z * z);

    const coarse = new PotentialArray({
      nz: 41, nr: 21, step: 0.5e-3, symmetry: CYLINDRICAL, z0: 5e-3,
    });
    const fine = new PotentialArray({
      nz: 81, nr: 41, step: 0.25e-3, symmetry: CYLINDRICAL, z0: 5e-3,
    });

    const ec = solveAgainstExact(coarse, exact).worst;
    const ef = solveAgainstExact(fine, exact).worst;
    const order = observedOrder(ec, ef);

    assert(
      order > 1.7 && order < 2.3,
      `expected 2nd-order convergence, observed order ${order.toFixed(2)} ` +
        `(errors ${ec.toExponential(2)} -> ${ef.toExponential(2)})`
    );
  });

  it('leaves a small Laplace residual everywhere', () => {
    const { grid, field } = buildEinzelLens({ gridStep: 1 });
    field.setVoltages({ housing: 0, entrance: 0, centre: -500, exit: 0 });
    const res = maxResidual(grid, field.phi);
    // Residual is h^2 * div grad phi in volts, against a 500 V scale.
    assert(res < 1e-5, `max Laplace residual ${res.toExponential(3)} V too large`);
  });

  it('obeys the maximum principle', () => {
    // A harmonic function attains its extrema only on the boundary. An
    // interior node exceeding every electrode potential would mean the solver
    // has invented charge that is not there.
    const { grid, field } = buildEinzelLens({ gridStep: 1 });
    const volts = { housing: 0, entrance: 0, centre: -500, exit: 0 };
    field.setVoltages(volts);
    const values = Object.values(volts);
    const lo = Math.min(...values);
    const hi = Math.max(...values);

    for (let k = 0; k < field.phi.length; k++) {
      if (grid.electrodeId[k] !== NO_ELECTRODE) continue;
      assert(
        field.phi[k] >= lo - 1e-9 && field.phi[k] <= hi + 1e-9,
        `free node potential ${field.phi[k]} escaped the electrode range [${lo}, ${hi}]`
      );
    }
  });

  it('rejects an over-relaxation factor outside the stable range', () => {
    const grid = new PotentialArray({ nz: 11, nr: 11, step: 1e-3 });
    grid.addElectrode('b');
    grid.paintEnclosure(0);
    let threw = false;
    try {
      relax(grid, new Float64Array(121), { omega: 2.5 });
    } catch {
      threw = true;
    }
    assert(threw, 'SOR diverges for omega >= 2 and must be refused');
  });

  it('picks an over-relaxation factor in the convergent range', () => {
    const omega = optimalOmega(247, 33);
    assert(omega > 1 && omega < 2, `omega ${omega} outside (1, 2)`);
  });
});

/* ------------------------------------------------------------------ */
/* field evaluation                                                    */
/* ------------------------------------------------------------------ */

describe('Field evaluation', () => {
  it('produces exactly zero radial field on the cylindrical axis', () => {
    // Rotational symmetry leaves a radial field at r = 0 with no direction to
    // point in, so any non-zero value is unphysical and would deflect an
    // on-axis ion out of the axis.
    const { grid, field } = buildEinzelLens({ gridStep: 1 });
    field.setVoltages({ housing: 0, entrance: 0, centre: -500, exit: 0 });
    for (let i = 0; i < grid.nz; i++) {
      assertClose(field.Er[grid.idx(i, 0)], 0, 0, `Er on axis at i=${i}`);
    }
  });

  it('recovers a uniform field from a linear potential', () => {
    // phi = -E0 y is harmonic, so the solver reproduces it exactly, and its
    // gradient must come back as a constant E0 everywhere - including at the
    // rim, where one-sided differences are used.
    const E0 = 1000; // V/m
    const grid = new PotentialArray({ nz: 31, nr: 21, step: 0.5e-3, symmetry: PLANAR });
    const { phi } = solveAgainstExact(grid, (z, y) => -E0 * y);

    // The solved map is used directly as the single basis solution, so a unit
    // volt on the lone "boundary" electrode reproduces it exactly.
    const field = new Field(grid, [phi]);
    field.setVoltages([1]);

    for (const [z, y] of [[5e-3, 3e-3], [10e-3, 7e-3], [2e-3, 9e-3]]) {
      const { Ez, Er } = field.fieldAt(z, y);
      assertRelClose(Er, E0, 1e-9, `transverse field at (${z}, ${y})`);
      assertClose(Ez, 0, 1e-6, `axial field at (${z}, ${y})`);
    }
  });

  it('superposes basis solutions to the same answer as a direct solve', () => {
    // Fast adjust is only legitimate because Laplace's equation is linear and
    // every basis solution shares the same homogeneous outer boundary. This
    // checks that claim numerically rather than trusting it.
    const { grid, field } = buildEinzelLens({ gridStep: 1 });
    const volts = { housing: 0, entrance: 120, centre: -500, exit: -30 };
    field.setVoltages(volts);

    const direct = new Float64Array(grid.nz * grid.nr);
    for (let k = 0; k < direct.length; k++) {
      const id = grid.electrodeId[k];
      if (id !== NO_ELECTRODE) direct[k] = volts[grid.electrodeNames[id]];
    }
    relax(grid, direct, { tolerance: 1e-11, maxSweeps: 200000 });

    let worst = 0;
    for (let k = 0; k < direct.length; k++) {
      worst = Math.max(worst, Math.abs(direct[k] - field.phi[k]));
    }
    assert(
      worst < 1e-6,
      `fast adjust differs from direct solve by ${worst.toExponential(3)} V`
    );
  });

  it('scales linearly with applied voltage', () => {
    const { field } = buildEinzelLens({ gridStep: 1 });
    field.setVoltages({ housing: 0, entrance: 0, centre: -100, exit: 0 });
    const a = field.potentialAt(mmToM(45), mmToM(2));
    field.setVoltages({ housing: 0, entrance: 0, centre: -300, exit: 0 });
    const b = field.potentialAt(mmToM(45), mmToM(2));
    assertRelClose(b, 3 * a, 1e-9, 'tripling the electrode voltage triples phi');
  });

  it('rejects an unknown electrode name', () => {
    const { field } = buildEinzelLens({ gridStep: 1 });
    let threw = false;
    try {
      field.setVoltages({ middle: -500 });
    } catch {
      threw = true;
    }
    assert(threw, 'a typo in an electrode name must not silently ground it');
  });
});

/* ------------------------------------------------------------------ */
/* integrators, isolated from the field solver                         */
/* ------------------------------------------------------------------ */

describe('Integrators against an analytic solution', () => {
  const k = 1e10; // V/m^2, gives a convenient oscillation period
  const field = harmonicField(k);
  const mass = amuToKg(100);
  const charge = ELEMENTARY_CHARGE;
  const omega = Math.sqrt((charge * k) / mass);
  const amplitude = 1e-3; // m

  const start = { mass, charge, x: amplitude, z: 0, vx: 0, vz: 0, t: 0 };
  const exactX = (t) => amplitude * Math.cos(omega * t);
  const exactVx = (t) => -amplitude * omega * Math.sin(omega * t);

  it('RK4 tracks simple harmonic motion', () => {
    const period = (2 * Math.PI) / omega;
    const steps = 2000;
    const dt = period / steps;
    const end = integrateFixed(stepRK4, field, start, dt, steps);
    assertRelClose(end.x, exactX(steps * dt), 1e-8, 'x after one full period');
    assertClose(end.vx, exactVx(steps * dt), amplitude * omega * 1e-8, 'vx after one period');
  });

  // Order is measured at a deliberately generic phase, not at a whole number
  // of periods. At t = T the oscillator sits at a turning point where
  // dx/dt = 0, so a phase error delta shifts x only by A delta^2 / 2 and the
  // measurement sees the wrong thing entirely: RK4's amplitude decay (5th
  // order) and the square of Verlet's phase error (4th order). Sampling
  // mid-swing, where dx/dt is large, exposes the phase error linearly and
  // recovers the methods' true orders.
  const SAMPLE_PHASE = 0.3; // fraction of a period

  it('RK4 converges at fourth order', () => {
    const period = (2 * Math.PI) / omega;
    const tEnd = SAMPLE_PHASE * period;
    const errAt = (steps) =>
      Math.abs(integrateFixed(stepRK4, field, start, tEnd / steps, steps).x - exactX(tEnd));
    const order = observedOrder(errAt(60), errAt(120));
    assert(
      order > 3.6 && order < 4.4,
      `expected 4th-order RK4, observed ${order.toFixed(2)}`
    );
  });

  it('velocity Verlet converges at second order', () => {
    const period = (2 * Math.PI) / omega;
    const tEnd = SAMPLE_PHASE * period;
    const errAt = (steps) =>
      Math.abs(integrateFixed(stepVerlet, field, start, tEnd / steps, steps).x - exactX(tEnd));
    const order = observedOrder(errAt(200), errAt(400));
    assert(
      order > 1.7 && order < 2.3,
      `expected 2nd-order Verlet, observed ${order.toFixed(2)}`
    );
  });

  it('keeps Verlet energy bounded over a long flight', () => {
    // The symplectic property: over hundreds of oscillations Verlet's energy
    // error must oscillate rather than accumulate. This is the reason the
    // method is kept alongside RK4 for future trapping work.
    const period = (2 * Math.PI) / omega;
    const dt = period / 200;
    let s = start;
    const E0 = totalEnergy(field, s);
    let worst = 0;
    for (let n = 0; n < 200 * 300; n++) {
      s = stepVerlet(field, s, dt);
      worst = Math.max(worst, Math.abs(totalEnergy(field, s) - E0) / Math.abs(E0));
    }
    assert(worst < 1e-3, `Verlet energy drift ${worst.toExponential(3)} over 300 periods`);
  });

  it('is exact for a uniform field', () => {
    // With constant acceleration the true solution is quadratic in t, which
    // both integrators reproduce identically. Any discrepancy here is a
    // coding error, not a truncation error.
    const a = 9.6486e8;
    const uniform = {
      grid: { step: 1e-3 },
      fieldAtCartesian: () => ({ Ex: (a * mass) / charge, Ez: 0 }),
      potentialAtCartesian: (x) => (-a * mass * x) / charge,
    };
    const ion = { mass, charge, x: 0, z: 0, vx: 0, vz: 4392.7, t: 0 };
    const dt = 1e-8;
    const steps = 200;
    const T = dt * steps;

    const rk4 = integrateFixed(stepRK4, uniform, ion, dt, steps);
    const verlet = integrateFixed(stepVerlet, uniform, ion, dt, steps);
    const expected = 0.5 * a * T * T;

    assertRelClose(rk4.x, expected, 1e-12, 'RK4 parabola');
    assertRelClose(verlet.x, expected, 1e-12, 'Verlet parabola');
    assertRelClose(rk4.z, 4392.7 * T, 1e-12, 'RK4 uniform drift');
  });
});

/* ------------------------------------------------------------------ */
/* the einzel lens as a whole                                          */
/* ------------------------------------------------------------------ */

describe('Einzel lens physics', () => {
  // One build shared across the tests below; solving is the expensive part.
  const built = buildEinzelLens({ gridStep: 0.5 });
  const { field, geometry } = built;
  const VOLTS = { housing: 0, entrance: 0, centre: -2000, exit: 0 };
  field.setVoltages(VOLTS);

  const SPEC = { mass: 100, charge: 1, energy: 1000, z: 0.5 };

  it('does no net work on a transmitted ion', () => {
    // The defining property of an einzel lens: entrance and exit sit at the
    // same potential, so an ion leaves with the energy it arrived with. This
    // holds exactly in the real device, independently of how well the lens
    // focuses, which makes it a strong end-to-end check on the field, the
    // units and the integrator at once.
    const ion = makeIon({ ...SPEC, x: 2 });
    const { points, stop } = flyIon(field, ion, { cfl: 0.05 });
    assert(stop === 'exited', `ion did not transmit (stop = ${stop})`);

    const last = points[points.length - 1];
    assert(last.z > mmToM(geometry.bounds.z6), 'ion must leave downstream of the lens');

    const kIn = joulesToEV(kineticEnergy(points[0]));
    const kOut = joulesToEV(kineticEnergy(last));
    assertRelClose(kOut, kIn, 2e-3, `KE in ${kIn.toFixed(3)} eV, out ${kOut.toFixed(3)} eV`);
  });

  it('conserves total energy along the flight', () => {
    const ion = makeIon({ ...SPEC, x: 2 });
    const { energyDrift } = flyIon(field, ion, { cfl: 0.05 });
    // The bound is loose on purpose. Force comes from E interpolated after
    // nodal differencing, while this diagnostic uses phi interpolated
    // directly; the two are consistent only to O(h^2), so a small residual
    // drift is expected discretisation error rather than an integrator fault.
    // See docs/PHYSICS.md.
    assert(
      energyDrift < 5e-3,
      `total energy drifted by ${energyDrift.toExponential(3)} of its initial value`
    );
  });

  it('reduces energy drift as the grid is refined', () => {
    // This pins the *cause* of the drift above. If it is discretisation error
    // in the field representation, refining the grid must shrink it. If it
    // were an integrator or units fault it would not care about the grid at
    // all, so this test distinguishes the two explanations.
    const driftAt = (gridStep) => {
      const { field: f } = buildEinzelLens({ gridStep });
      f.setVoltages(VOLTS);
      return flyIon(f, makeIon({ ...SPEC, x: 2 }), { cfl: 0.05 }).energyDrift;
    };
    const coarse = driftAt(1.0);
    const fine = driftAt(0.5);
    assert(
      fine < coarse,
      `halving the grid step did not reduce energy drift ` +
        `(${coarse.toExponential(2)} -> ${fine.toExponential(2)}), so the drift ` +
        'is not discretisation error and needs a different explanation'
    );
  });

  it('keeps an on-axis ion on the axis', () => {
    const ion = makeIon({ ...SPEC, x: 0, angle: 0 });
    const { points } = flyIon(field, ion, { cfl: 0.05 });
    for (const p of points) {
      assertClose(p.x, 0, 1e-15, `on-axis ion strayed to x = ${p.x} m at z = ${p.z}`);
    }
  });

  it('is mirror symmetric about the axis', () => {
    // The geometry is rotationally symmetric, so trajectories launched at +x
    // and -x must be exact reflections. An asymmetry would point at a sign
    // error in the radial field or in the axis folding.
    const up = flyIon(field, makeIon({ ...SPEC, x: 3 }), { cfl: 0.05 });
    const down = flyIon(field, makeIon({ ...SPEC, x: -3 }), { cfl: 0.05 });

    assert(up.points.length === down.points.length, 'mirrored flights differ in length');
    for (let n = 0; n < up.points.length; n++) {
      assertClose(up.points[n].x, -down.points[n].x, 1e-15, `x mirror at step ${n}`);
      assertClose(up.points[n].z, down.points[n].z, 1e-15, `z mirror at step ${n}`);
    }
  });

  it('focuses a parallel beam to a point downstream of the lens', () => {
    const rays = parallelBeam({ ...SPEC, count: 5, maxOffset: 3 });
    const foci = [];
    for (const ray of rays) {
      if (ray.x === 0) continue; // the axial ray never crosses
      const { points } = flyIon(field, ray, { cfl: 0.05 });
      const f = axialCrossing(points);
      assert(f !== null, `ray at x = ${ray.x} m never crossed the axis`);
      foci.push({ x0: Math.abs(ray.x), f });
    }

    assert(foci.length >= 4, 'expected at least four non-axial rays');
    for (const { f } of foci) {
      assert(
        f > mmToM(geometry.bounds.z4),
        `focus at z = ${(f * 1e3).toFixed(1)} mm is not downstream of the lens centre`
      );
    }
  });

  it('shows positive spherical aberration', () => {
    // A real electrostatic lens is stronger for rays further off axis, so the
    // outermost ray must cross the axis *before* the innermost one. An ideal
    // thin lens would put them at the same point; reproducing the aberration
    // is evidence the fringe fields are being integrated, not idealised away.
    const inner = flyIon(field, makeIon({ ...SPEC, x: 1 }), { cfl: 0.05 });
    const outer = flyIon(field, makeIon({ ...SPEC, x: 4 }), { cfl: 0.05 });

    const fInner = axialCrossing(inner.points);
    const fOuter = axialCrossing(outer.points);
    assert(fInner !== null && fOuter !== null, 'both rays must cross the axis');
    assert(
      fOuter < fInner,
      `outer ray focused at ${(fOuter * 1e3).toFixed(2)} mm, inner at ` +
        `${(fInner * 1e3).toFixed(2)} mm; expected the outer ray to focus earlier`
    );
  });

  it('agrees between RK4 and velocity Verlet', () => {
    // Two methods of different order and different structure landing on the
    // same trajectory is evidence about the trajectory, not about either
    // method's internal consistency.
    const spec = { ...SPEC, x: 3 };
    const a = flyIon(field, makeIon(spec), { method: 'rk4', cfl: 0.02 });
    const b = flyIon(field, makeIon(spec), { method: 'verlet', cfl: 0.02 });

    const fa = axialCrossing(a.points);
    const fb = axialCrossing(b.points);
    assert(fa !== null && fb !== null, 'both integrators must produce a crossing');
    assertRelClose(fb, fa, 1e-4, 'focal position from RK4 vs Verlet');
  });

  it('strengthens focusing as the centre voltage rises', () => {
    // Monotonicity is a weaker claim than a focal-length formula but it is
    // one an ion optician would immediately check, and it fails loudly if the
    // sign of the radial force is wrong.
    const focusAt = (Vc) => {
      field.setVoltages({ housing: 0, entrance: 0, centre: Vc, exit: 0 });
      const { points } = flyIon(field, makeIon({ ...SPEC, x: 3 }), { cfl: 0.05 });
      return axialCrossing(points);
    };

    // Both biases must focus *inside* the modelled region for the comparison
    // to mean anything. A weak lens has a long focal length: at -1000 V this
    // geometry focuses past the end of the domain, so axialCrossing correctly
    // reports no crossing and the test would be comparing against nothing.
    const weak = focusAt(-1500);
    const strong = focusAt(-3000);
    field.setVoltages(VOLTS); // restore for any later test

    assert(
      weak !== null && strong !== null,
      'both biases must focus within the domain for this comparison to apply'
    );
    assert(
      strong < weak,
      `stronger bias focused at ${(strong * 1e3).toFixed(1)} mm, weaker at ` +
        `${(weak * 1e3).toFixed(1)} mm; expected the stronger lens to focus sooner`
    );
  });

  it('focuses for either polarity of the centre electrode', () => {
    // An einzel lens is convergent whether the centre electrode accelerates
    // or decelerates the ion. The decelerating case must still transmit, so
    // the bias is kept well below the beam energy.
    for (const Vc of [-1500, 600]) {
      field.setVoltages({ housing: 0, entrance: 0, centre: Vc, exit: 0 });
      const { points, stop } = flyIon(field, makeIon({ ...SPEC, x: 3 }), { cfl: 0.05 });
      assert(stop === 'exited', `ion failed to transmit at Vc = ${Vc} V (${stop})`);
      const f = axialCrossing(points);
      assert(f !== null, `no focus at Vc = ${Vc} V`);
    }
    field.setVoltages(VOLTS);
  });

  it('reflects an ion that cannot climb the centre electrode barrier', () => {
    // A decelerating bias above the beam energy is a potential barrier. The
    // ion must turn around rather than tunnel through it, which is both
    // correct physics and a check that the axial force reverses sign.
    field.setVoltages({ housing: 0, entrance: 0, centre: 3000, exit: 0 });
    const { points } = flyIon(field, makeIon({ ...SPEC, x: 0 }), { cfl: 0.05 });
    const last = points[points.length - 1];
    assert(last.vz < 0, `expected reflection, but the ion left with vz = ${last.vz} m/s`);
    field.setVoltages(VOLTS);
  });
});
