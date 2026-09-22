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

import {
  describe,
  it,
  assert,
  assertClose,
  assertRelClose,
  assertFinite,
} from './harness.js';

import {
  ELEMENTARY_CHARGE,
  ATOMIC_MASS_UNIT,
  ELECTRON_MASS,
  VACUUM_PERMITTIVITY,
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
  flyBeam,
  totalEnergy,
  kineticEnergy,
} from '../src/integrator.js';
import {
  currentShares,
  spaceChargeField,
  lineChargeDensity,
} from '../src/spacecharge.js';
import { makeIon, parallelBeam, axialCrossing, focalCrossing } from '../src/ion.js';
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

  // Math.max, never `if (err > worst)`. A comparison silently drops NaN, so a
  // solution that has blown up entirely would report zero error and pass.
  let worst = 0;
  for (let j = 0; j < grid.nr; j++) {
    for (let i = 0; i < grid.nz; i++) {
      if (grid.isElectrode(i, j)) continue;
      const value = phi[grid.idx(i, j)];
      assertFinite(value, `potential at node (${i}, ${j})`);
      worst = Math.max(worst, Math.abs(value - exact(grid.zAt(i), grid.rAt(j))));
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
/* the metal/vacuum interface                                          */
/* ------------------------------------------------------------------ */

describe('Field at a conductor surface', () => {
  /**
   * A grounded slab filling y <= 2 mm, with a linear potential above it.
   *
   * The potential map is written analytically rather than solved, so the only
   * thing under test is the differencing. Everything is registered as one
   * electrode because the electrode MAP is what the gradient code consults;
   * the values come from the supplied basis array.
   */
  function slabField(step = 0.5e-3, slabTop = 2e-3, top = 10e-3, V = 100) {
    const nr = Math.round(top / step) + 1;
    const grid = new PotentialArray({ nz: 21, nr, step, symmetry: PLANAR });
    const id = grid.addElectrode('metal');
    grid.paintEnclosure(id);
    grid.paint(id, (z, y) => y <= slabTop + 1e-12);

    const phi = new Float64Array(grid.nz * nr);
    for (let j = 0; j < nr; j++) {
      const y = grid.rAt(j);
      const value = y <= slabTop ? 0 : (V * (y - slabTop)) / (top - slabTop);
      for (let i = 0; i < grid.nz; i++) phi[grid.idx(i, j)] = value;
    }

    const field = new Field(grid, [phi]);
    field.setVoltages([1]);
    return { grid, field, expected: -V / (top - slabTop) };
  }

  it('returns the full surface field, not half of it', () => {
    // A central difference taken AT a surface node reaches one node into the
    // conductor, where the potential is pinned at the electrode value. It
    // therefore returns (phi_vacuum - V) / 2h where the true surface
    // derivative is (phi_vacuum - V) / h - exactly half, and being a factor
    // rather than a truncation term it does not shrink with refinement.
    for (const step of [1e-3, 0.5e-3, 0.25e-3]) {
      const { grid, field, expected } = slabField(step);
      const j = Math.round(2e-3 / step); // the surface node
      const Ey = field.Er[grid.idx(10, j)];
      assertRelClose(
        Ey,
        expected,
        1e-9,
        `surface field at h = ${step * 1e3} mm (half would be ${expected / 2})`
      );
    }
  });

  it('gives zero field inside a conductor', () => {
    const { grid, field } = slabField();
    const j = Math.round(1e-3 / 0.5e-3); // well inside the slab
    assertClose(field.Er[grid.idx(10, j)], 0, 1e-9, 'field inside metal');
  });

  it('does not extrapolate the field outside the domain', () => {
    // Clamping the cell index but not the interpolation fraction lets the
    // bilinear form run away outside the grid: the default lens reported
    // 807 V at r = 66 mm inside a grounded 16 mm housing.
    const { field } = buildEinzelLens({ gridStep: 1 });
    field.setVoltages({ housing: 0, entrance: 0, centre: -500, exit: 0 });

    let inside = 0;
    for (let k = 0; k < field.phi.length; k++) {
      inside = Math.max(inside, Math.abs(field.phi[k]));
    }
    for (const [z, r] of [[0.06, 0.021], [0.06, 0.066], [0.06, -0.004], [0.4, 0.008]]) {
      const value = Math.abs(field.potentialAt(z, r));
      assert(
        value <= inside + 1e-9,
        `potential at (${z}, ${r}) m is ${value} V, outside the solved range of ` +
          `${inside} V - the interpolator is extrapolating`
      );
    }
  });
});

/* ------------------------------------------------------------------ */
/* absolute scale and geometry                                         */
/* ------------------------------------------------------------------ */

describe('Absolute scale', () => {
  it('pins the millimetre conversion to an absolute length', () => {
    // Every lens assertion in this suite is a ratio or an ordering, and most
    // express their bounds through the same mmToM they rely on. A global
    // rescale - mmToM returning centimetres, say - would leave all of them
    // green while making the modelled device ten times too big. This is the
    // one assertion that cannot be satisfied that way.
    assertRelClose(mmToM(1), 1e-3, 0, 'one millimetre in metres');
    assertRelClose(mmToM(123), 0.123, 1e-15, '123 mm in metres');

    const { grid, geometry } = buildEinzelLens({ gridStep: 0.5 });
    assertRelClose(geometry.totalLength, 123, 1e-12, 'lens length in mm');
    assertRelClose(grid.zLength, 0.123, 1e-12, 'lens length in metres');
    assertRelClose(grid.step, 5e-4, 1e-15, 'grid step in metres');
    assertRelClose(grid.rLength, 0.016, 1e-12, 'housing radius in metres');
  });

  it('recovers a known field strength in volts per metre', () => {
    // An absolute V/m, with the gap written as a literal in metres so it does
    // not route through the converter under test.
    const gap = 0.008; // m
    const V = 100;
    const { grid, field } = (() => {
      const g = new PotentialArray({ nz: 21, nr: 21, step: 0.5e-3, symmetry: PLANAR });
      const { phi } = solveAgainstExact(g, (z, y) => (V * y) / 0.01);
      const f = new Field(g, [phi]);
      f.setVoltages([1]);
      return { grid: g, field: f };
    })();
    assertRelClose(
      field.fieldAt(5e-3, 5e-3).Er,
      -V / 0.01,
      1e-9,
      'transverse field in V/m'
    );
    assert(grid.step === 0.5e-3, 'grid step literal');
    assert(gap > 0, 'gap defined');
  });

  it('paints the geometry that was requested', () => {
    // Nothing else checks that the painted metal matches the specification,
    // so an electrode silently painted at half its length would pass the
    // whole suite.
    const { grid, geometry } = buildEinzelLens({ gridStep: 0.5 });
    const extent = (name) => {
      const id = grid.electrodeNames.indexOf(name);
      let zMin = Infinity;
      let zMax = -Infinity;
      let rMin = Infinity;
      let rMax = -Infinity;
      for (let j = 0; j < grid.nr; j++) {
        for (let i = 0; i < grid.nz; i++) {
          if (grid.electrodeId[grid.idx(i, j)] !== id) continue;
          zMin = Math.min(zMin, grid.zAt(i));
          zMax = Math.max(zMax, grid.zAt(i));
          rMin = Math.min(rMin, grid.rAt(j));
          rMax = Math.max(rMax, grid.rAt(j));
        }
      }
      return { zMin, zMax, rMin, rMax };
    };

    const b = geometry.bounds;
    const h = grid.step;
    for (const [name, z0, z1] of [
      ['entrance', b.z1, b.z2],
      ['centre', b.z3, b.z4],
      ['exit', b.z5, b.z6],
    ]) {
      const e = extent(name);
      assertClose(e.zMin, mmToM(z0), h, `${name} starts at z = ${z0} mm`);
      assertClose(e.zMax, mmToM(z1), h, `${name} ends at z = ${z1} mm`);
      assertClose(e.rMin, mmToM(geometry.boreRadius), h, `${name} bore radius`);
      assertClose(e.rMax, mmToM(geometry.outerRadius), h, `${name} outer radius`);
    }

    // The device is axially symmetric by construction.
    const counts = grid.electrodeNodeCounts();
    const idOf = (n) => grid.electrodeNames.indexOf(n);
    assert(
      counts[idOf('entrance')] === counts[idOf('exit')],
      'entrance and exit cylinders must be painted identically'
    );
  });

  it('refuses a domain rim that no electrode owns', () => {
    // The relaxation never updates rim nodes, so an unpainted rim is frozen
    // at 0 V - an invisible grounded box no voltage setting can override.
    const grid = new PotentialArray({ nz: 11, nr: 11, step: 1e-3, symmetry: PLANAR });
    grid.addElectrode('partial');
    grid.paint(0, (z, y) => y <= 0); // only the floor
    let threw = false;
    try {
      solveBasis(grid);
    } catch {
      threw = true;
    }
    assert(threw, 'an unowned rim must be refused, not silently grounded');
  });
});

/* ------------------------------------------------------------------ */
/* focus detection                                                     */
/* ------------------------------------------------------------------ */

describe('Axial crossing', () => {
  const pt = (x, z, vz = 1) => ({ x, z, vz, vx: 0 });

  it('returns the first forward crossing, not the last', () => {
    // The back focal point is where the ray FIRST reaches the axis. A ray
    // that crosses, diverges and is turned again by a later element focuses
    // at the first crossing; reporting the last would name the wrong plane.
    const path = [pt(1, 0), pt(-1, 2), pt(1, 4), pt(-1, 6)];
    assertClose(axialCrossing(path), 1, 1e-12, 'first crossing');
  });

  it('finds a crossing that lands exactly on the axis', () => {
    // A strict sign product misses this: a.x * 0 is never negative.
    assertClose(axialCrossing([pt(1, 0), pt(0.5, 1), pt(0, 2)]), 2, 1e-12, 'exact landing');
  });

  it('reports no crossing for a ray that never leaves the axis', () => {
    assert(
      axialCrossing([pt(0, 0), pt(0, 1), pt(0, 2)]) === null,
      'an on-axis ray has no crossing, and certainly not one at the far wall'
    );
  });

  it('reports no crossing for a ray diverging from the axis', () => {
    assert(axialCrossing([pt(0, 0), pt(1, 1), pt(2, 2)]) === null, 'diverging ray');
  });

  it('ignores crossings made while travelling backwards', () => {
    const returning = [pt(1, 4, -1), pt(-1, 2, -1), pt(-2, 0, -1)];
    assert(axialCrossing(returning) === null, 'a reflected ray has no focus');
  });

  it('extrapolates a focus that lies beyond the modelled region', () => {
    // Beyond the last electrode the ray is straight, so the crossing follows
    // from the exit position and slope. Without this, a weak lens whose focus
    // lies past the end of the grid reports "no focus" even though it works.
    const exit = { x: 2e-3, z: 0.1, vx: -1e3, vz: 1e5, t: 0 };
    const f = focalCrossing([{ x: 3e-3, z: 0.09, vx: -1e3, vz: 1e5, t: 0 }, exit]);
    assert(f !== null && f.extrapolated, 'expected an extrapolated focus');
    // z = 0.1 - 2e-3 * (1e5 / -1e3) = 0.1 + 0.2 = 0.3
    assertRelClose(f.z, 0.3, 1e-12, 'extrapolated crossing');
  });

  it('does not extrapolate a focus for a diverging ray', () => {
    const exit = { x: 2e-3, z: 0.1, vx: +1e3, vz: 1e5, t: 0 };
    assert(focalCrossing([exit, exit]) === null, 'diverging ray has no focus ahead');
  });
});

/* ------------------------------------------------------------------ */
/* space charge                                                        */
/* ------------------------------------------------------------------ */

describe('Space charge', () => {
  /** Rays evenly sampling a uniform-density beam of radius R. */
  const sampleBeam = (R, n) =>
    Array.from({ length: n + 1 }, (_, i) => (i * R) / n);

  it('carries current shares that sum to one', () => {
    const shares = currentShares(sampleBeam(4e-3, 12));
    const total = shares.reduce((a, b) => a + b, 0);
    assertRelClose(total, 1, 1e-12, 'current shares must partition the beam');
  });

  it('weights shares by annulus area, not by ray count', () => {
    // A ray at radius r stands for the annulus halfway to each neighbour, so
    // its share is (r_out^2 - r_in^2) / R^2. Treating every ray as an equal
    // slice would put far too much charge near the axis and overstate the
    // defocusing. Checked against the closed form rather than an arbitrary
    // ratio.
    const R = 4e-3;
    const n = 10;
    const radii = sampleBeam(R, n);
    const shares = currentShares(radii);
    const d = R / n;

    for (let i = 0; i <= n; i++) {
      const inner = i === 0 ? 0 : radii[i] - d / 2;
      const outer = i === n ? R : radii[i] + d / 2;
      const expected = (outer * outer - inner * inner) / (R * R);
      assertRelClose(shares[i], expected, 1e-12, `share of ray ${i}`);
    }
  });

  it('treats a ray and its mirror image as one ring', () => {
    // parallelBeam launches signed offsets, so every non-zero radius appears
    // twice - the two halves of a single ring in the meridional plane. They
    // must split one annulus between them, not claim one each, or the beam
    // carries twice the current it was given.
    const signed = [-3e-3, -1.5e-3, 0, 1.5e-3, 3e-3];
    const shares = currentShares(signed.map(Math.abs));
    const unique = currentShares([0, 1.5e-3, 3e-3]);

    assertRelClose(shares[0] + shares[4], unique[2], 1e-12, 'outer ring total');
    assertRelClose(shares[1] + shares[3], unique[1], 1e-12, 'middle ring total');
    assertRelClose(shares[0], shares[4], 1e-12, 'mirrored halves share equally');

    // And they must feel an identical field, or the beam is not symmetric.
    const E = spaceChargeField(signed.map(Math.abs), shares, 1e-9);
    assertRelClose(E[0], E[4], 1e-12, 'mirrored halves feel the same field');
    assertRelClose(E[1], E[3], 1e-12, 'mirrored halves feel the same field');
  });

  it('reproduces the analytic field of a uniform cylindrical beam', () => {
    // Gauss's law on a uniform beam of radius R and line density lambda gives
    //
    //     E_r(r) = lambda r / (2 pi eps0 R^2)   for r <= R
    //
    // i.e. the field rises LINEARLY from the axis to the beam edge. This is
    // the single most important check on the model: it is a closed-form
    // result the ring sampling has to reproduce.
    const R = 4e-3;
    const n = 40;
    const radii = sampleBeam(R, n);
    const shares = currentShares(radii);
    const lambda = 1e-9; // C/m

    const E = spaceChargeField(radii, shares, lambda);
    const k = 1 / (2 * Math.PI * VACUUM_PERMITTIVITY);

    // Skip the innermost rays: within a couple of ring spacings of the axis
    // the discrete annuli cannot resolve the 1/r, and the model is documented
    // as unreliable there.
    for (let i = 5; i <= n; i++) {
      const expected = (lambda * radii[i] * k) / (R * R);
      assertRelClose(E[i], expected, 0.05, `E_r at r = ${radii[i] * 1e3} mm`);
    }
  });

  it('converges towards the analytic field as rays are added', () => {
    const R = 4e-3;
    const lambda = 1e-9;
    const k = 1 / (2 * Math.PI * VACUUM_PERMITTIVITY);

    const errorAt = (n) => {
      const radii = sampleBeam(R, n);
      const E = spaceChargeField(radii, currentShares(radii), lambda);
      // Compare at the beam edge, where the sampling is best resolved.
      const expected = (lambda * R * k) / (R * R);
      return Math.abs(E[n] - expected) / expected;
    };

    assert(
      errorAt(40) < errorAt(10),
      `refining the ray sampling must improve the field ` +
        `(${errorAt(10).toExponential(2)} -> ${errorAt(40).toExponential(2)})`
    );
  });

  it('gives exactly zero field on the axis', () => {
    // By symmetry there is no direction for a radial field to point in at
    // r = 0, exactly as for the electrode field.
    const radii = sampleBeam(4e-3, 10);
    const E = spaceChargeField(radii, currentShares(radii), 1e-9);
    assertClose(E[0], 0, 0, 'space-charge field on the axis');
  });

  it('shields the interior from charge outside it', () => {
    // The cylindrical shell theorem: a uniform shell of charge exerts no net
    // force on anything inside it. Adding current at large radius must not
    // change the field felt at small radius.
    const inner = [1e-3, 2e-3];
    const withOuter = [1e-3, 2e-3, 8e-3];
    const lambda = 1e-9;

    // Give the inner two rays identical absolute shares in both cases, so the
    // only difference is the presence of the outer ring.
    const sharesA = [0.25, 0.75];
    const sharesB = [0.25, 0.75, 4.0]; // outer ring carries far more

    const a = spaceChargeField(inner, sharesA, lambda);
    const b = spaceChargeField(withOuter, sharesB, lambda);

    assertRelClose(b[0], a[0], 1e-12, 'inner ray must not feel the outer shell');
    assertRelClose(b[1], a[1], 1e-12, 'middle ray must not feel the outer shell');
  });

  it('scales linearly with beam current', () => {
    const radii = sampleBeam(4e-3, 10);
    const shares = currentShares(radii);
    const a = spaceChargeField(radii, shares, 1e-9);
    const b = spaceChargeField(radii, shares, 3e-9);
    for (let i = 1; i < radii.length; i++) {
      assertRelClose(b[i], 3 * a[i], 1e-12, `tripling lambda at ray ${i}`);
    }
  });

  it('makes a slower beam denser', () => {
    // lambda = I / v: the same current at half the speed is twice the charge
    // per unit length. This is why space charge bites hardest where an optic
    // decelerates the beam.
    assertRelClose(
      lineChargeDensity(1e-6, 1000),
      2 * lineChargeDensity(1e-6, 2000),
      1e-12,
      'halving the speed doubles the line density'
    );
  });
});

describe('Space charge in flight', () => {
  const built = buildEinzelLens({ gridStep: 0.5 });
  const { field } = built;
  field.setVoltages({ housing: 0, entrance: 0, centre: -2000, exit: 0 });
  const SPEC = { mass: 100, charge: 1, energy: 1000, z: 0.5 };
  const beam = () => parallelBeam({ ...SPEC, count: 9, maxOffset: 3 });

  it('changes nothing at all when the current is zero', () => {
    // A self-field model that perturbs the answer when switched off is worse
    // than no model, so this is checked against flyIon to machine precision
    // rather than to a tolerance.
    const single = flyIon(field, makeIon({ ...SPEC, x: 3 }), { cfl: 0.05 });
    const { tracks } = flyBeam(field, [makeIon({ ...SPEC, x: 3 })], {
      cfl: 0.05,
      beamCurrent: 0,
    });
    const a = single.points[single.points.length - 1];
    const b = tracks[0].points[tracks[0].points.length - 1];
    assertClose(b.x, a.x, 1e-18, 'final x with space charge disabled');
    assertClose(b.z, a.z, 1e-18, 'final z with space charge disabled');
    assert(tracks[0].stop === single.stop, 'stop reason must match flyIon');
  });

  it('pushes the focus downstream as the current rises', () => {
    // Space charge is repulsive for a single-species beam, so it always
    // opposes the lens and the crossover must move DOWNSTREAM. If any current
    // pulled it upstream, the sign of the self-field would be wrong.
    //
    // Kept in the weak-space-charge regime where a crossover still exists;
    // the strong regime is the next test.
    const focusAt = (beamCurrent) => {
      const { tracks } = flyBeam(field, beam(), { cfl: 0.05, beamCurrent });
      return focalCrossing(tracks[tracks.length - 1].points);
    };

    let previous = -Infinity;
    for (const I of [0, 1e-6, 3e-6]) {
      const f = focusAt(I);
      assert(f !== null, `no crossover at ${I} A, where one is still expected`);
      assert(
        f.z > previous,
        `focus moved upstream at ${I} A: ` +
          `${(previous * 1e3).toFixed(1)} mm -> ${(f.z * 1e3).toFixed(1)} mm`
      );
      previous = f.z;
    }
  });

  it('replaces the point focus with a finite waist at high current', () => {
    // This is the defining behaviour of a space-charge-dominated beam and it
    // is not a numerical artefact. The self-field goes as 1/r, so as the beam
    // converges the repulsion diverges: it cannot be brought to a point. The
    // beam reaches a minimum radius - a waist - and expands again, and that
    // waist grows with current.
    //
    // A test that insisted on a crossover at every current would be asserting
    // physics that is simply false.
    const waistOf = (beamCurrent) => {
      const { tracks } = flyBeam(field, beam(), { cfl: 0.05, beamCurrent });
      const outer = tracks[tracks.length - 1];
      let waist = Infinity;
      for (const p of outer.points) {
        // Downstream of the last electrode, so this is the beam's own waist
        // and not something happening inside the lens.
        if (p.z < mmToM(83)) continue;
        waist = Math.min(waist, Math.abs(p.x));
      }
      return { waist, crossover: focalCrossing(outer.points) };
    };

    let previous = -Infinity;
    for (const I of [1e-5, 2e-5, 5e-5]) {
      const { waist, crossover } = waistOf(I);
      assert(
        crossover === null,
        `expected no point focus at ${I} A, but one was reported`
      );
      assert(
        waist > 1e-5,
        `expected a finite waist at ${I} A, got ${(waist * 1e3).toFixed(4)} mm`
      );
      assert(
        waist > previous,
        `waist must grow with current: ${(previous * 1e3).toFixed(3)} mm -> ` +
          `${(waist * 1e3).toFixed(3)} mm at ${I} A`
      );
      previous = waist;
    }
  });

  it('keeps an on-axis ion on the axis', () => {
    // The self-field vanishes on the axis by symmetry, so space charge must
    // not deflect the central ray.
    const { tracks } = flyBeam(field, beam(), { cfl: 0.05, beamCurrent: 2e-5 });
    const centre = tracks[(tracks.length - 1) / 2];
    assertClose(centre.start?.x ?? centre.points[0].x, 0, 1e-18, 'centre ray starts on axis');
    for (const p of centre.points) {
      assertClose(p.x, 0, 1e-15, `centre ray strayed to x = ${p.x}`);
    }
  });

  it('stays mirror symmetric under space charge', () => {
    const { tracks } = flyBeam(field, beam(), { cfl: 0.05, beamCurrent: 2e-5 });
    const n = tracks.length;
    for (let i = 0; i < (n - 1) / 2; i++) {
      const lo = tracks[i];
      const hi = tracks[n - 1 - i];
      assert(lo.points.length === hi.points.length, `pair ${i} differs in length`);
      const a = lo.points[lo.points.length - 1];
      const b = hi.points[hi.points.length - 1];
      assertClose(a.x, -b.x, 1e-15, `pair ${i} transverse mirror`);
      assertClose(a.z, b.z, 1e-15, `pair ${i} axial mirror`);
    }
  });

  it('flies every ion on a common time base', () => {
    // The self-field is only defined for a beam whose members are at the same
    // instant, so lockstep is a correctness requirement, not an optimisation.
    const { tracks } = flyBeam(field, beam(), { cfl: 0.05, beamCurrent: 1e-5 });
    const reference = tracks[0];
    for (const t of tracks) {
      const steps = Math.min(t.points.length, reference.points.length);
      for (let n = 0; n < steps; n++) {
        assertClose(
          t.points[n].t,
          reference.points[n].t,
          1e-18,
          `ion times diverged at sample ${n}`
        );
      }
    }
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

  it('reports a reflected ion as reflected, not transmitted', () => {
    // A decelerating bias above the beam energy turns the lens into an ion
    // mirror. The ion leaves through the ENTRANCE face, and calling that
    // "exited" would count it as transmitted and let its backwards axis
    // crossing masquerade as a focal length.
    // Launched on-axis, so the return path is clean. An off-axis ion above
    // the barrier is usually driven into metal on the way back out, which is
    // itself correct behaviour but tests a different branch.
    field.setVoltages({ housing: 0, entrance: 0, centre: 3000, exit: 0 });
    const { points, stop } = flyIon(field, makeIon({ ...SPEC, x: 0 }), { cfl: 0.05 });
    const last = points[points.length - 1];
    field.setVoltages(VOLTS);

    assert(stop === 'reflected', `expected 'reflected', got '${stop}'`);
    assert(last.vz < 0, 'a reflected ion must be travelling backwards');
    assert(last.z < points[0].z, 'a reflected ion must end upstream of its launch');
    assert(
      focalCrossing(points) === null,
      'a reflected ion has no focus; its return path must not be reported as one'
    );
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

  it('gives the same focus for every ion mass at fixed energy', () => {
    // Electrostatic optics depends only on E/q: the force is qE and the
    // trajectory shape is set by the ratio of kinetic to potential energy, in
    // which mass cancels entirely. Only the flight TIME scales, as sqrt(m).
    //
    // This is real device physics and is external to the code, so it cannot
    // be satisfied by a self-consistent-but-wrong implementation. It is also
    // the sharpest available check on the unit chain: a stray factor of u, of
    // e, or of 1000 anywhere between makeIon and the force would break the
    // mass independence immediately.
    const foci = [];
    const times = [];
    for (const mass of [1, 100, 10000]) {
      const { points } = flyIon(field, makeIon({ ...SPEC, mass, x: 3 }), { cfl: 0.05 });
      const f = axialCrossing(points);
      assert(f !== null, `no focus for mass ${mass} u`);
      foci.push(f);
      times.push(points[points.length - 1].t);
    }
    assertRelClose(foci[1], foci[0], 1e-6, 'focus at 100 u vs 1 u');
    assertRelClose(foci[2], foci[0], 1e-6, 'focus at 10000 u vs 1 u');

    // Flight time must scale as sqrt(m): 100x the mass is 10x the time.
    assertRelClose(times[1] / times[0], 10, 1e-4, 'flight time scaling with mass');
    assertRelClose(times[2] / times[0], 100, 1e-4, 'flight time scaling with mass');
  });

  it('transmits a negative ion at mirrored polarity identically', () => {
    // Reversing the sign of both the charge and every electrode gives an
    // identical force at every point, so the trajectory must be identical.
    // Nothing else in the suite flies anything but a singly-charged cation.
    const a = flyIon(field, makeIon({ ...SPEC, x: 3 }), { cfl: 0.05 });

    field.setVoltages({ housing: 0, entrance: 0, centre: 2000, exit: 0 });
    const b = flyIon(field, makeIon({ ...SPEC, charge: -1, x: 3 }), { cfl: 0.05 });
    field.setVoltages(VOLTS);

    assert(b.stop === 'exited', `mirrored anion did not transmit (${b.stop})`);
    assertRelClose(
      axialCrossing(b.points),
      axialCrossing(a.points),
      1e-9,
      'anion at mirrored polarity must follow the cation trajectory exactly'
    );
  });

  it('reports a ray that hits the housing as an electrode strike', () => {
    // The outer wall is real metal, unlike the end faces. An ion driven into
    // it must be a strike, not an escape.
    const steep = makeIon({ ...SPEC, x: 0, angle: 45 });
    const { points, stop } = flyIon(field, steep, { cfl: 0.05 });
    const last = points[points.length - 1];
    assert(stop === 'electrode', `expected 'electrode', got '${stop}'`);
    assert(
      Math.abs(last.x) <= built.grid.rLength + 1e-12,
      'a strike must be recorded at or inside the wall, not beyond it'
    );
  });

  it('never produces a non-finite potential or field', () => {
    // Comparison-based error accumulators silently pass on NaN, so a blown-up
    // solve can look perfect. Check the arrays directly.
    for (let k = 0; k < field.phi.length; k++) {
      assertFinite(field.phi[k], `phi[${k}]`);
      assertFinite(field.Ez[k], `Ez[${k}]`);
      assertFinite(field.Er[k], `Er[${k}]`);
    }
  });

  it('converges every basis solution', () => {
    for (let e = 0; e < built.solverReports.length; e++) {
      assert(
        built.solverReports[e].converged,
        `basis solution for "${built.grid.electrodeNames[e]}" did not converge`
      );
    }
  });

  it('reduces energy drift at second order, not merely below a threshold', () => {
    // A bare bound is a regression pin: it is flat in the time step, so it
    // cannot detect a broken integrator, and it fails for non-physics reasons
    // as soon as the default grid changes. The documented claim is that the
    // drift is the O(h^2) force/potential interpolation mismatch, so test the
    // scaling instead of the value.
    const driftAt = (gridStep) => {
      const { field: f } = buildEinzelLens({ gridStep });
      f.setVoltages(VOLTS);
      return flyIon(f, makeIon({ ...SPEC, x: 2 }), { cfl: 0.05 }).energyDrift;
    };
    const coarse = driftAt(1.0);
    const mid = driftAt(0.5);
    const fine = driftAt(0.25);

    const r1 = coarse / mid;
    const r2 = mid / fine;
    assert(
      r1 > 3 && r1 < 5.5 && r2 > 3 && r2 < 5.5,
      `expected ~4x drift reduction per halving, got ${r1.toFixed(2)} and ${r2.toFixed(2)}`
    );
  });

  it('keeps energy drift independent of the time step', () => {
    // The counterpart to the test above. If the drift really is a grid
    // artefact, shrinking the time step by 20x must not change it - and if it
    // DOES change, the drift is integrator error after all and the O(h^2)
    // story in docs/PHYSICS.md is wrong. This is the test that a degraded
    // integrator cannot pass.
    const at = (cfl) =>
      flyIon(field, makeIon({ ...SPEC, x: 2 }), { cfl }).energyDrift;
    const coarse = at(0.2);
    const fine = at(0.01);
    assertRelClose(fine, coarse, 0.1, 'energy drift must not depend on the time step');
  });
});
