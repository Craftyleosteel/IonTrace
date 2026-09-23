/**
 * Beamline and element validation.
 *
 * The quadrupole is the interesting one: it is the first element whose field
 * is not axisymmetric, so it is checked against the closed-form ideal
 * quadrupole potential, against Mathieu stability theory, and against the
 * symmetry the rod geometry must have. The composition layer is checked for
 * the things that go wrong when elements are chained - coordinate translation,
 * element lookup at the joins, and the step size adapting to the most
 * demanding element rather than the first.
 */

import { describe, it, assert, assertClose, assertRelClose } from './harness.js';

import { Beamline } from '../src/beamline.js';
import { createElement, ELEMENT_TYPES, needsRebuild } from '../src/elements/index.js';
import { createQuadrupole, MATHIEU_Q_LIMIT } from '../src/elements/quadrupole.js';
import { createDrift } from '../src/elements/drift.js';
import { makeIon } from '../src/ion.js';
import { flyIon, flyBeam, kineticEnergy } from '../src/integrator.js';
import {
  mmToM,
  mToMm,
  joulesToEV,
  ELEMENTARY_CHARGE,
  ATOMIC_MASS_UNIT,
} from '../src/constants.js';

/* ------------------------------------------------------------------ */
/* quadrupole field                                                    */
/* ------------------------------------------------------------------ */

describe('Quadrupole field', () => {
  // Coarse enough to keep the suite quick, fine enough to resolve the rods.
  const quad = createQuadrupole({ gridStep: 0.25, fieldRadius: 4, rodRadius: 4.6 });
  const r0 = mmToM(4);

  it('has the sign structure of a quadrupole', () => {
    // Positive pole voltage means the potential is positive along x and
    // negative along y. An ion of positive charge is therefore pushed AWAY
    // from the axis along x and TOWARDS it along y - converging in one plane
    // and diverging in the other, at the same instant. This is the property
    // that makes a DC quadrupole a singlet and an RF one a filter.
    const onX = quad.potentialAt(r0 * 0.5, 0, quad.length / 2, 0);
    const onY = quad.potentialAt(0, r0 * 0.5, quad.length / 2, 0);
    assert(onX > 0, `potential on the x axis should be positive, got ${onX}`);
    assert(onY < 0, `potential on the y axis should be negative, got ${onY}`);
    assertRelClose(onX, -onY, 0.02, 'the two planes must be mirror images');
  });

  it('is zero on the axis to within solver residue', () => {
    // phi = (x^2 - y^2)/r0^2 vanishes at the origin and, by the reflection
    // symmetries of the rod set, so does the field.
    //
    // Unlike the cylindrical axis, this is not IMPOSED: there the radial
    // field is set to exactly zero because rotational symmetry leaves it no
    // direction to point in. Here the axis is an ordinary interior point of a
    // planar grid, so what remains is the relaxation's own residue. It must
    // therefore be judged against the field scale rather than against zero.
    const z = quad.length / 2;
    const scale = Math.abs(quad.fieldAt(r0 * 0.3, 0, z, 0).Ex);
    assert(scale > 1e3, 'the reference field should be substantial');

    assertClose(quad.potentialAt(0, 0, z, 0), 0, 1e-6, 'potential on axis');
    const E = quad.fieldAt(0, 0, z, 0);
    assert(
      Math.abs(E.Ex) / scale < 1e-6,
      `Ex on axis is ${(Math.abs(E.Ex) / scale).toExponential(2)} of the field scale`
    );
    assert(
      Math.abs(E.Ey) / scale < 1e-6,
      `Ey on axis is ${(Math.abs(E.Ey) / scale).toExponential(2)} of the field scale`
    );
  });

  it('matches the ideal quadrupole potential near the axis', () => {
    // phi = W (x^2 - y^2) / r0^2 is the closed form. Real round rods produce
    // higher multipoles, so agreement is only expected well inside r0 - which
    // is exactly the claim being tested.
    const z = quad.length / 2;
    const W = quad.params.dcVoltage + quad.params.rfAmplitude; // t = 0, phase 0

    for (const f of [0.2, 0.3, 0.4]) {
      for (const [x, y] of [[f, 0], [0, f], [f * 0.7, f * 0.7]]) {
        const X = x * r0;
        const Y = y * r0;
        const expected = (W * (X * X - Y * Y)) / (r0 * r0);
        const actual = quad.potentialAt(X, Y, z, 0);
        assertRelClose(
          actual,
          expected,
          0.06,
          `phi at (${(x).toFixed(2)}, ${(y).toFixed(2)}) r0`
        );
      }
    }
  });

  it('gives a field that rises linearly with distance from the axis', () => {
    // E = -grad phi = -2W(x, -y)/r0^2, so |E| is linear in r. This is the
    // defining property of a quadrupole and the reason its restoring force is
    // harmonic.
    const z = quad.length / 2;
    const at = (f) => quad.fieldAt(f * r0, 0, z, 0).Ex;
    const e1 = at(0.2);
    const e2 = at(0.4);
    assertRelClose(e2, 2 * e1, 0.05, 'field must double when x doubles');
  });

  it('is four-fold symmetric', () => {
    // Rotating by 90 degrees must flip the sign; rotating by 180 must leave
    // it unchanged. A painting error in one rod shows up here immediately.
    const z = quad.length / 2;
    const x = r0 * 0.35;
    const a = quad.potentialAt(x, 0, z, 0);
    const b = quad.potentialAt(0, x, z, 0);
    const c = quad.potentialAt(-x, 0, z, 0);
    const d = quad.potentialAt(0, -x, z, 0);
    assertRelClose(c, a, 1e-9, '180 degree rotation');
    assertRelClose(d, b, 1e-9, '180 degree rotation');
    assertRelClose(b, -a, 0.02, '90 degree rotation flips the sign');
  });

  it('oscillates at the drive frequency', () => {
    // W(t) = U + V cos(Omega t). At half a period the RF has reversed, so a
    // converging plane has become a diverging one.
    const z = quad.length / 2;
    const x = r0 * 0.3;
    const period = 1 / (quad.params.frequency * 1e6);
    const atZero = quad.potentialAt(x, 0, z, 0);
    const atHalf = quad.potentialAt(x, 0, z, period / 2);
    const atFull = quad.potentialAt(x, 0, z, period);
    assertRelClose(atHalf, -atZero, 1e-9, 'half a period reverses the field');
    assertRelClose(atFull, atZero, 1e-9, 'a full period returns it');
  });

  it('reports a shortest period only when the RF is on', () => {
    assertRelClose(quad.shortestPeriod, 1e-6, 1e-12, '1 MHz period');
    const dc = createQuadrupole({ gridStep: 0.5, rfAmplitude: 0 });
    assert(dc.shortestPeriod === null, 'a DC quadrupole has no period to resolve');
  });

  it('computes Mathieu a and q correctly', () => {
    // a = 8 z e U / (m r0^2 Omega^2),  q = 4 z e V / (m r0^2 Omega^2)
    const m = 100;
    const z = 1;
    const { a, q } = quad.mathieu(m, z);

    const omega = 2 * Math.PI * quad.params.frequency * 1e6;
    const denom = m * ATOMIC_MASS_UNIT * r0 * r0 * omega * omega;
    assertRelClose(
      q,
      (4 * z * ELEMENTARY_CHARGE * quad.params.rfAmplitude) / denom,
      1e-12,
      'q parameter'
    );
    assertRelClose(
      a,
      (8 * z * ELEMENTARY_CHARGE * quad.params.dcVoltage) / denom,
      1e-12,
      'a parameter'
    );
  });

  it('scales q inversely with mass and with the square of frequency', () => {
    // These scalings are what make the device a MASS filter: at fixed RF,
    // each mass sits at a different q.
    const base = quad.mathieu(100, 1).q;
    assertRelClose(quad.mathieu(200, 1).q, base / 2, 1e-12, 'q ~ 1/m');
    assertRelClose(quad.mathieu(100, 2).q, base * 2, 1e-12, 'q ~ charge');

    const faster = createQuadrupole({ gridStep: 0.5, frequency: 2 });
    assertRelClose(
      faster.mathieu(100, 1).q,
      quad.mathieu(100, 1).q / 4,
      1e-12,
      'q ~ 1/Omega^2'
    );
  });
});

/* ------------------------------------------------------------------ */
/* quadrupole dynamics                                                 */
/* ------------------------------------------------------------------ */

describe('Quadrupole dynamics', () => {
  /** A short line: drift, quadrupole, drift. */
  function line(quadParams = {}) {
    return new Beamline([
      createElement('drift', { length: 5, bore: 10 }),
      createElement('quadrupole', { gridStep: 0.3, length: 60, ...quadParams }),
      createElement('drift', { length: 5, bore: 10 }),
    ]);
  }

  const SPEC = { mass: 100, charge: 1, energy: 10 };

  it('keeps an on-axis ion on the axis', () => {
    // The field vanishes on the axis at all times, so there is nothing to
    // deflect the ion in either transverse direction.
    //
    // The bound is not machine zero, and deliberately so. The quadrupole's
    // axis field is zero by symmetry but is not IMPOSED to be, unlike the
    // cylindrical axis, so a residual of order 1e-6 V/m survives the
    // relaxation and nudges the ion by about a nanometre over the flight.
    // Demanding exact zero would be demanding that SOR converge perfectly.
    //
    // The bound is stated against the aperture, which is what makes it a
    // physical claim - the stray must be negligible compared with the
    // dimension that decides whether the ion survives - rather than an
    // arbitrary small number that would need renegotiating whenever the
    // solver tolerance moved.
    //
    // The ion must also be a STABLE one. On the axis of a quadrupole running
    // outside its stability region, an ion sits at an unstable equilibrium:
    // the restoring force is proportional to displacement, so a perturbation
    // of any size grows exponentially, and round-off is a perturbation. An
    // unstable ion leaving the axis is then correct behaviour, and testing
    // "stays on axis" with one would be testing the wrong thing.
    const bl = line();
    const tol = bl.elements[1].bore * 1e-4;
    const stable = { mass: 400, charge: 1, energy: 10 };
    assert(bl.elements[1].mathieu(stable.mass, 1).q < 0.6, 'probe ion must be stable');

    const { points } = flyIon(bl, makeIon({ ...stable, x: 0, y: 0 }), { cfl: 0.05 });
    for (const p of points) {
      assertClose(p.x, 0, tol, `strayed in x at z = ${p.z}`);
      assertClose(p.y, 0, tol, `strayed in y at z = ${p.z}`);
    }
  });

  it('confines a stable ion and loses an unstable one', () => {
    // The first stability region at a = 0 ends at q = 0.90803. An ion below
    // it should be transmitted; one well above it should not. This is the
    // single most important behavioural claim a quadrupole can make, and it
    // comes from Mathieu theory, not from this code.
    const bl = line({ rfAmplitude: 300, dcVoltage: 0, frequency: 1 });
    const quad = bl.elements[1];

    const qOf = (mass) => quad.mathieu(mass, 1).q;

    // Heavy ion: small q, comfortably inside the stable region. Kept well
    // clear of the 0.908 boundary rather than just inside it, so the test is
    // about stability rather than about how precisely the edge is resolved.
    const heavy = 400;
    // Light ion: large q, well outside it.
    const light = 40;
    assert(qOf(heavy) < 0.6, `expected a stable q, got ${qOf(heavy).toFixed(3)}`);
    assert(qOf(light) > 1.3, `expected an unstable q, got ${qOf(light).toFixed(3)}`);

    const fly = (mass) =>
      flyIon(bl, makeIon({ mass, charge: 1, energy: 10, x: 0.4, y: 0.4 }), {
        cfl: 0.05,
      });

    const stable = fly(heavy);
    const unstable = fly(light);

    assert(
      stable.stop === 'exited',
      `ion at q = ${qOf(heavy).toFixed(3)} should transmit, got ${stable.stop}`
    );
    assert(
      unstable.stop === 'electrode',
      `ion at q = ${qOf(light).toFixed(3)} should be lost, got ${unstable.stop}`
    );
  });

  it('keeps a stable ion bounded inside the field radius', () => {
    // Being transmitted is not enough: a stable ion's excursion must stay
    // bounded rather than growing along the rods.
    const bl = line({ rfAmplitude: 300, dcVoltage: 0, frequency: 1 });
    const { points } = flyIon(
      bl,
      makeIon({ mass: 400, charge: 1, energy: 10, x: 0.4, y: 0.4 }),
      { cfl: 0.05 }
    );
    const r0 = bl.elements[1].bore;
    let worst = 0;
    for (const p of points) worst = Math.max(worst, Math.hypot(p.x, p.y));
    assert(
      worst < r0,
      `stable ion reached ${mToMm(worst).toFixed(2)} mm, outside r0 = ${mToMm(r0)} mm`
    );
  });

  it('takes the ion out of the launch plane only when it starts off both axes', () => {
    // An ion launched on the x axis stays in the x-z plane, because Ey is
    // proportional to y and vanishes there. This is the quadrupole's own
    // symmetry, and it is worth pinning: if it failed, the 3D machinery would
    // be introducing motion that the physics does not have.
    const bl = line();
    const onAxis = flyIon(bl, makeIon({ mass: 200, charge: 1, energy: 10, x: 0.5, y: 0 }), {
      cfl: 0.05,
    });
    for (const p of onAxis.points) {
      // A nanometre, for the same reason as the on-axis test above: the
      // symmetry is exact but the relaxation's residue is not.
      assertClose(p.y, 0, 1e-9, 'an ion launched on the x axis must stay in that plane');
    }

    const offAxis = flyIon(
      bl,
      makeIon({ mass: 200, charge: 1, energy: 10, x: 0.5, y: 0.5 }),
      { cfl: 0.05 }
    );
    const movedInY = offAxis.points.some((p) => Math.abs(p.y) > 1e-9);
    assert(movedInY, 'an ion off both axes must move in y as well as x');
  });

  it('drives the two transverse planes in antiphase', () => {
    // At any instant the force converges in one plane and diverges in the
    // other, so an ion launched symmetrically at x = y must move oppositely
    // in the two coordinates at first.
    const bl = line({ rfAmplitude: 300, phase: 0 });
    const { points } = flyIon(
      bl,
      makeIon({ mass: 200, charge: 1, energy: 10, x: 0.5, y: 0.5 }),
      { cfl: 0.02 }
    );
    const inside = points.filter(
      (p) => p.z > bl.elements[1].zStart + 1e-3 && p.z < bl.elements[1].zEnd
    );
    assert(inside.length > 10, 'expected samples inside the rods');
    const p = inside[Math.floor(inside.length * 0.02)];
    assert(
      (Math.abs(p.x) - 0.5e-3) * (Math.abs(p.y) - 0.5e-3) < 0,
      'the two planes must respond in opposite senses'
    );
  });

  it('resolves the RF rather than aliasing it', () => {
    // If the step were free to grow past the RF period, an ion would sample
    // one phase over and over and an unstable trajectory could look stable.
    // The step must stay a small fraction of the period.
    const bl = line({ frequency: 2 });
    const { points } = flyIon(
      bl,
      makeIon({ mass: 200, charge: 1, energy: 10, x: 0.5, y: 0.2 }),
      { cfl: 0.05 }
    );
    const period = bl.shortestPeriod;
    let worst = 0;
    for (let i = 1; i < points.length; i++) {
      worst = Math.max(worst, points[i].t - points[i - 1].t);
    }
    assert(
      worst < period / 4,
      `largest step ${worst.toExponential(2)} s is not small against the ` +
        `${period.toExponential(2)} s RF period`
    );
  });
});

/* ------------------------------------------------------------------ */
/* composition                                                         */
/* ------------------------------------------------------------------ */

describe('Beamline composition', () => {
  it('lays elements end to end without gaps or overlaps', () => {
    const bl = new Beamline([
      createElement('drift', { length: 10 }),
      createElement('drift', { length: 25 }),
      createElement('drift', { length: 5 }),
    ]);
    assertRelClose(bl.length, mmToM(40), 1e-12, 'total length');
    assertRelClose(bl.elements[1].zStart, mmToM(10), 1e-12, 'second starts where first ends');
    assertRelClose(bl.elements[2].zStart, mmToM(35), 1e-12, 'third starts where second ends');
  });

  it('finds the right element at every point, including the joins', () => {
    const bl = new Beamline([
      createElement('drift', { length: 10 }),
      createElement('drift', { length: 10 }),
    ]);
    assert(bl.elementAt(mmToM(5)) === bl.elements[0], 'inside the first');
    assert(bl.elementAt(mmToM(10)) === bl.elements[1], 'the join belongs to the second');
    assert(bl.elementAt(mmToM(15)) === bl.elements[1], 'inside the second');
    assert(bl.elementAt(mmToM(20)) === bl.elements[1], 'the far end belongs to the last');
    assert(bl.elementAt(mmToM(25)) === null, 'beyond the column');
  });

  it('translates coordinates into each element', () => {
    // The lens must act where it was placed. Putting a drift in front of it
    // and asking for the field at the same OFFSET inside the lens must give
    // the same answer both times.
    const alone = new Beamline([createElement('einzel', { gridStep: 1 })]);
    const shifted = new Beamline([
      createElement('drift', { length: 30 }),
      createElement('einzel', { gridStep: 1 }),
    ]);

    const offset = mmToM(40);
    const a = alone.potentialAt3D(mmToM(2), 0, offset, 0);
    const b = shifted.potentialAt3D(mmToM(2), 0, mmToM(30) + offset, 0);
    assertRelClose(b, a, 1e-12, 'the same point inside the lens');
    assert(Math.abs(a) > 1, 'the test point should be somewhere the field is non-trivial');
  });

  it('adds, removes and reorders elements', () => {
    const bl = new Beamline([
      createElement('drift', { length: 10 }),
      createElement('aperture', { gridStep: 1, margin: 6 }),
    ]);
    const before = bl.length;

    bl.add(createElement('drift', { length: 7 }));
    assertRelClose(bl.length, before + mmToM(7), 1e-12, 'adding extends the column');

    assert(bl.move(0, 1), 'move should succeed');
    assert(bl.elements[0].type === 'aperture', 'the aperture moved to the front');
    assertRelClose(bl.elements[0].zStart, 0, 1e-15, 'and was re-laid at zero');

    assert(!bl.move(0, -1), 'cannot move past the start');

    bl.remove(0);
    assertRelClose(bl.length, mmToM(17), 1e-12, 'removing shortens the column');
  });

  it('takes its step size from the most demanding element', () => {
    // A coarse drift must not let an ion skip through a finely resolved
    // quadrupole beside it.
    const bl = new Beamline([
      createDrift({ length: 20 }),
      createElement('quadrupole', { gridStep: 0.2 }),
    ]);
    assertRelClose(bl.lengthScale, mmToM(0.2), 1e-12, 'finest element wins');
    assertRelClose(bl.shortestPeriod, 1e-6, 1e-12, 'RF period is exposed to the integrator');
  });

  it('reports no field outside the column', () => {
    const bl = new Beamline([createElement('einzel', { gridStep: 1 })]);
    const E = bl.fieldAt3D(0, 0, bl.length * 2, 0);
    assertClose(E.Ex, 0, 0, 'no field beyond the end');
    assertClose(E.Ez, 0, 0, 'no field beyond the end');
    assert(!bl.strikes(0, 0, bl.length * 2), 'nothing to strike beyond the end');
  });

  it('warns when two live elements are butted together', () => {
    // Each element is solved in isolation with grounded end faces, so the
    // field where two live ones meet is not a true solution for the pair.
    // Silently producing it would be the dishonest option.
    const bad = new Beamline([
      createElement('einzel', { gridStep: 1 }),
      createElement('aperture', { gridStep: 1, margin: 6 }),
    ]);
    assert(
      bad.warnings.some((w) => w.includes('adjacent')),
      'expected an adjacency warning'
    );

    const good = new Beamline([
      createElement('einzel', { gridStep: 1 }),
      createElement('drift', { length: 15 }),
      createElement('aperture', { gridStep: 1, margin: 6 }),
    ]);
    assert(
      !good.warnings.some((w) => w.includes('adjacent')),
      'a drift between them should clear the warning'
    );
  });

  it('stops an ion on a drift tube wall', () => {
    const bl = new Beamline([createDrift({ length: 40, bore: 3 })]);
    const { stop } = flyIon(
      bl,
      makeIon({ mass: 100, charge: 1, energy: 10, x: 0, angle: 20 }),
      { cfl: 0.05 }
    );
    assert(stop === 'electrode', `expected a wall strike, got ${stop}`);
  });

  it('carries an ion through a field-free drift in a straight line', () => {
    // The simplest possible end-to-end check on the composition layer.
    const bl = new Beamline([createDrift({ length: 50, bore: 10 })]);
    const ion = makeIon({ mass: 100, charge: 1, energy: 100, x: 1, angle: 2 });
    const { points, stop } = flyIon(bl, ion, { cfl: 0.2 });
    assert(stop === 'exited', `expected transmission, got ${stop}`);

    const last = points[points.length - 1];
    const expectedX = ion.x + (ion.vx / ion.vz) * (last.z - ion.z);
    assertRelClose(last.x, expectedX, 1e-9, 'a field-free ion travels straight');
    assertRelClose(
      joulesToEV(kineticEnergy(last)),
      100,
      1e-9,
      'and does not change energy'
    );
  });
});

/* ------------------------------------------------------------------ */
/* element registry                                                    */
/* ------------------------------------------------------------------ */

describe('Element registry', () => {
  it('builds every registered type from its own defaults', () => {
    for (const [type, spec] of Object.entries(ELEMENT_TYPES)) {
      // Coarse grids keep this quick; the point is that each type builds.
      const coarse = type === 'quadrupole' ? { gridStep: 0.5 } : { gridStep: 1 };
      const el = createElement(type, type === 'drift' ? {} : coarse);
      assert(el.length > 0, `${type} must have a positive length`);
      assert(el.bore > 0, `${type} must have a positive bore`);
      assert(typeof el.fieldAt === 'function', `${type} must expose fieldAt`);
      assert(el.label === spec.label, `${type} label should match the registry`);
    }
  });

  it('knows which parameters need a re-solve and which do not', () => {
    // This is the distinction the whole fast-adjust architecture rests on.
    assert(!needsRebuild('einzel', 'voltage'), 'voltages are a weighted sum, never a re-solve');
    assert(needsRebuild('einzel', 'boreRadius'), 'geometry changes the boundary');
    assert(!needsRebuild('quadrupole', 'rfAmplitude'), 'RF amplitude only scales the solution');
    assert(needsRebuild('quadrupole', 'fieldRadius'), 'moving the rods is geometry');
    assert(needsRebuild('einzel', 'somethingUnknown'), 'unknown keys must err towards re-solving');
  });

  it('changes an einzel voltage without re-solving', () => {
    // Fast adjust: the field must scale exactly, because it is a weighted sum
    // of solutions already computed.
    const el = createElement('einzel', { gridStep: 1, voltage: -100 });
    const a = el.potentialAt(mmToM(2), 0, el.length / 2);
    el.setVoltage(-300);
    const b = el.potentialAt(mmToM(2), 0, el.length / 2);
    assertRelClose(b, 3 * a, 1e-9, 'tripling the electrode triples the potential');
  });

  it('rejects an unknown element type', () => {
    let threw = false;
    try {
      createElement('sextupole');
    } catch {
      threw = true;
    }
    assert(threw, 'unknown types must be refused, not silently ignored');
  });

  it('exposes the Mathieu limit used to judge stability', () => {
    assertRelClose(MATHIEU_Q_LIMIT, 0.90803, 1e-5, 'first stability region cut-off');
  });
});
