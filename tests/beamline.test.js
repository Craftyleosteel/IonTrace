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

import { Beamline, exitsOf } from '../src/beamline.js';
import {
  createElement,
  ELEMENT_TYPES,
  needsRebuild,
  startingParams,
  fieldRange,
} from '../src/elements/index.js';
import {
  createQuadrupole,
  MATHIEU_Q_LIMIT,
  MATHIEU_Q_WORKING,
  amplitudeForQ,
} from '../src/elements/quadrupole.js';
import { createDrift } from '../src/elements/drift.js';
import {
  createBender,
  matchedVoltage,
  DEFLECTOR_CONSTANT,
  DEFLECTOR_DESIGN_PHASE,
} from '../src/elements/bender.js';
import {
  compose,
  translation,
  yawFrame,
  toGlobal,
  toLocal,
  vectorToGlobal,
  forwardOf,
} from '../src/frames.js';
import {
  tunableKnobs,
  scoreBeamline,
  optimizeVoltages,
  applyKnob,
  readKnobs,
  TUNABLE,
} from '../src/optimize.js';
import { symmetricEigen, beamQuality, refineNullSpace } from '../src/reduced.js';
import { axisymmetricRuns, canShareGrid, decayLength } from '../src/column.js';
import { serialise, restore } from '../src/scene.js';
import { makeIon, discBeam } from '../src/ion.js';
import { flyIon, flyBeam, kineticEnergy } from '../src/integrator.js';
import {
  mmToM,
  mToMm,
  joulesToEV,
  ELEMENTARY_CHARGE,
  ATOMIC_MASS_UNIT,
  VACUUM_PERMITTIVITY,
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

  it('solves on a grid that is exactly symmetric about the axis', () => {
    // The enclosure must share the rods' four-fold symmetry. Deriving the node
    // count as round(2R/h) + 1 and then spanning -R upward leaves the domain
    // lopsided by up to half a grid step, which is enough to make E_y non-zero
    // on the y = 0 plane - and an ion launched there is then pushed out of it.
    // Measured with the lopsided grid, a planar ion drifted 197 um, five per
    // cent of the aperture, from a 0.15 mm asymmetry in the box.
    for (const gridStep of [0.15, 0.25, 0.35, 0.4]) {
      const q = createQuadrupole({ gridStep, housingRadius: 12 });
      const g = q.grid;
      assertClose(g.zAt(0), -g.zAt(g.nz - 1), 1e-15, `x extent at h = ${gridStep}`);
      assertClose(g.rAt(0), -g.rAt(g.nr - 1), 1e-15, `y extent at h = ${gridStep}`);
      assert(g.nz % 2 === 1, `node count must be odd so the axis is a node (h = ${gridStep})`);
      // And therefore the axis really is a node, sitting at exactly zero.
      assertClose(g.zAt((g.nz - 1) / 2), 0, 1e-15, 'axis node position');
    }
  });

  it('has no transverse field on either symmetry plane', () => {
    // E_y = 0 everywhere on y = 0, and E_x = 0 everywhere on x = 0, by the
    // reflection symmetry of the rod set. This is the property that keeps a
    // planar beam planar, and it is the one the grid asymmetry destroyed.
    const z = quad.length / 2;
    const scale = Math.abs(quad.fieldAt(r0 * 0.5, 0, z, 0).Ex);
    for (const f of [0.1, 0.25, 0.5, 0.75]) {
      const onX = quad.fieldAt(f * r0, 0, z, 0);
      const onY = quad.fieldAt(0, f * r0, z, 0);
      assert(
        Math.abs(onX.Ey) / scale < 1e-4,
        `E_y on the y = 0 plane at x = ${f} r0 is ${(Math.abs(onX.Ey) / scale).toExponential(2)} of scale`
      );
      assert(
        Math.abs(onY.Ex) / scale < 1e-4,
        `E_x on the x = 0 plane at y = ${f} r0 is ${(Math.abs(onY.Ex) / scale).toExponential(2)} of scale`
      );
    }
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
    // Stated from the element's own frequency rather than a literal, so this
    // keeps testing the relation if the default frequency ever moves again.
    assertRelClose(
      quad.shortestPeriod,
      1 / (quad.params.frequency * 1e6),
      1e-12,
      'one RF period'
    );
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

    // Twice the frequency of whatever this one runs at, so the factor of four
    // is the physics and not an assumption about the default.
    const faster = createQuadrupole({ gridStep: 0.5, frequency: quad.params.frequency * 2 });
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
/* the beam is three-dimensional                                       */
/* ------------------------------------------------------------------ */

describe('Beam in three dimensions', () => {
  const line = () =>
    new Beamline([
      createElement('drift', { length: 6, bore: 6 }),
      createElement('quadrupole', { gridStep: 0.35, length: 120, rfAmplitude: 250, frequency: 2 }),
      createElement('drift', { length: 10, bore: 6 }),
    ]);
  const SPEC = { mass: 100, charge: 1, energy: 50 };

  it('spreads a disc beam over both transverse axes', () => {
    // parallelBeam puts every ion on the x axis, which is a line of ions and
    // not a beam. discBeam must actually fill the disc.
    const ions = discBeam({ ...SPEC, count: 24, radius: 2 });
    const offAxisY = ions.filter((i) => Math.abs(i.y) > 1e-9);
    assert(
      offAxisY.length > ions.length / 2,
      `expected most ions off the y = 0 plane, got ${offAxisY.length} of ${ions.length}`
    );

    // Uniform areal density: half the ions inside r/sqrt(2), half outside.
    const radii = ions.map((i) => Math.hypot(i.x, i.y)).sort((a, b) => a - b);
    const median = radii[Math.floor(radii.length / 2)];
    assertRelClose(median, mmToM(2) / Math.SQRT2, 0.12, 'median radius of a uniform disc');
    assert(radii[radii.length - 1] <= mmToM(2) + 1e-12, 'no ion outside the beam radius');
  });

  it('keeps a planar beam exactly planar in an axisymmetric column', () => {
    // Here the confinement to a plane is EXACT, not approximate. The radial
    // field is resolved onto y as E_r * (y/r), which is identically zero when
    // y is zero - no solver residue can leak in, because the multiplication
    // by y does it. So this is asserted at machine zero.
    const bl = new Beamline([
      createElement('drift', { length: 6, bore: 8 }),
      createElement('einzel', { gridStep: 0.6, voltage: -200, boreRadius: 6, entryDrift: 18, exitDrift: 18 }),
      createElement('drift', { length: 10, bore: 8 }),
    ]);
    const { points } = flyIon(bl, makeIon({ ...SPEC, x: 1.2, y: 0 }), { cfl: 0.05 });
    for (const p of points) {
      assertClose(p.y, 0, 0, `planar ion left its plane at z = ${mToMm(p.z).toFixed(1)} mm`);
    }
    assert(points.some((p) => Math.abs(p.x) < mmToM(1.1)), 'the lens should still act in x');
  });

  it('keeps a planar beam planar in a quadrupole to within solver residue', () => {
    // Same statement, weaker guarantee, and the difference is worth knowing.
    // A quadrupole's E_y comes from interpolating a solved transverse grid,
    // where it is zero on the y = 0 plane only to the relaxation's residue
    // rather than by construction. The ion therefore drifts sub-nanometre
    // amounts out of its plane - physically nothing against a 4 mm aperture,
    // but not machine zero, and it would be wrong to claim otherwise.
    // A thousandth of the aperture. The residual E_y does not simply displace
    // the ion, it DRIVES it: the quadrupole restores in y, so the residue
    // acts as a forcing term on an oscillator and the ion rings at a small
    // amplitude rather than settling at a small offset. A micron on a 4 mm
    // aperture is still nothing, but the bound has to leave room for the
    // oscillation rather than for a static offset.
    const bl = line();
    const tol = bl.elements[1].bore * 1e-3;
    const { points } = flyIon(bl, makeIon({ ...SPEC, x: 1.2, y: 0 }), { cfl: 0.05 });
    for (const p of points) {
      assertClose(p.y, 0, tol, `planar ion left its plane at z = ${mToMm(p.z).toFixed(1)} mm`);
    }
    assert(
      points.some((p) => Math.abs(p.x) > mmToM(1.3)),
      'the ion should still be moving in x'
    );
  });

  it('moves an ion in both transverse axes when it starts off both', () => {
    // The same line, the same element, one ion moved off the symmetry plane.
    const bl = line();
    const { points } = flyIon(bl, makeIon({ ...SPEC, x: 1.2, y: 0.8 }), { cfl: 0.05 });
    let maxX = 0;
    let maxY = 0;
    for (const p of points) {
      maxX = Math.max(maxX, Math.abs(p.x));
      maxY = Math.max(maxY, Math.abs(p.y));
    }
    assert(maxY > mmToM(0.5), `expected real motion in y, got ${mToMm(maxY).toFixed(3)} mm`);
    assert(maxX > mmToM(0.5), `expected real motion in x, got ${mToMm(maxX).toFixed(3)} mm`);
  });

  it('drives the two transverse planes out of phase in a quadrupole', () => {
    // The signature of a quadrupole: while one plane is being squeezed the
    // other is being stretched. Both coordinates must oscillate, and their
    // extremes must not coincide.
    const bl = line();
    const { points } = flyIon(bl, makeIon({ ...SPEC, x: 1.2, y: 1.2 }), { cfl: 0.02 });
    const q = bl.elements[1];
    const inside = points.filter((p) => p.z > q.zStart && p.z < q.zEnd);
    assert(inside.length > 50, 'need a decent sample inside the rods');

    // Count sign changes: a confined ion oscillates in both planes.
    const crossings = (key) => {
      let n = 0;
      for (let i = 1; i < inside.length; i++) {
        if (inside[i - 1][key] * inside[i][key] < 0) n++;
      }
      return n;
    };
    assert(crossings('x') > 2, `expected several x oscillations, got ${crossings('x')}`);
    assert(crossings('y') > 2, `expected several y oscillations, got ${crossings('y')}`);
  });

  it('sends a divergent disc beam outward, not sideways', () => {
    // Divergence is radial: every ion's transverse velocity points away from
    // the axis along its own azimuth, so the beam expands as a cone rather
    // than fanning out in one plane.
    for (const ion of discBeam({ ...SPEC, count: 16, radius: 2, divergence: 3 })) {
      const r = Math.hypot(ion.x, ion.y);
      if (r < 1e-9) continue;
      const vr = (ion.vx * ion.x + ion.vy * ion.y) / r;
      const vTransverse = Math.hypot(ion.vx, ion.vy);
      assertRelClose(vr, vTransverse, 1e-9, 'transverse velocity must be purely radial');
      // And the angle must scale with radius, reaching the full value at the
      // beam edge.
      const angle = (Math.atan2(vTransverse, ion.vz) * 180) / Math.PI;
      assertRelClose(angle, 3 * (r / mmToM(2)), 1e-6, 'divergence scales with radius');
    }
  });

  it('expands a divergent beam and leaves a collimated one alone', () => {
    const bl = new Beamline([createElement('drift', { length: 60, bore: 10 })]);
    // Compared against the beam's OWN launch radius, not against the nominal
    // one: a uniform-density disc places its outermost ion a little inside
    // the nominal edge, so the nominal radius is a property of the
    // distribution rather than the position of any particular ion.
    const widthAfter = (divergence) => {
      const ions = discBeam({ ...SPEC, count: 12, radius: 1, divergence });
      let launched = 0;
      let worst = 0;
      for (const ion of ions) {
        launched = Math.max(launched, Math.hypot(ion.x, ion.y));
        const { points } = flyIon(bl, ion, { cfl: 0.2 });
        const p = points[points.length - 1];
        worst = Math.max(worst, Math.hypot(p.x, p.y));
      }
      return { launched, worst };
    };
    const collimated = widthAfter(0);
    const diverging = widthAfter(3);
    assertRelClose(
      collimated.worst,
      collimated.launched,
      1e-9,
      'a collimated beam keeps the radius it was launched with'
    );
    assert(
      diverging.worst > collimated.worst * 2,
      `a 3 degree beam should widen over 60 mm: ${mToMm(collimated.worst).toFixed(2)} -> ` +
        `${mToMm(diverging.worst).toFixed(2)} mm`
    );
  });
});

/* ------------------------------------------------------------------ */
/* placements                                                          */
/* ------------------------------------------------------------------ */

describe('Frames', () => {
  it('round-trips a point through a placement', () => {
    const f = compose(translation(0.1, -0.02, 0.3), yawFrame(0.7));
    const p = [0.004, -0.003, 0.05];
    const back = toLocal(f, toGlobal(f, p));
    for (let i = 0; i < 3; i++) {
      assertClose(back[i], p[i], 1e-15, `component ${i} must survive the round trip`);
    }
  });

  it('rotates vectors without displacing them', () => {
    // A field is a vector: it turns with the frame but must not pick up the
    // frame's origin. Treating it as a point would be a category error that
    // happens to compile.
    const f = compose(translation(1, 2, 3), yawFrame(Math.PI / 2));
    const v = vectorToGlobal(f, [0, 0, 1]);
    assertClose(Math.hypot(v[0], v[1], v[2]), 1, 1e-15, 'rotation preserves length');
    assertClose(v[0], -1, 1e-15, 'a quarter turn sends +z to -x');
    assertClose(v[2], 0, 1e-15, 'and leaves nothing along z');
  });

  it('composes placements in order', () => {
    // Two quarter turns make a half turn, and the origins accumulate along
    // the rotated axes rather than the global ones.
    const quarter = compose(translation(0, 0, 1), yawFrame(Math.PI / 2));
    const half = compose(quarter, quarter);
    const dir = forwardOf(half);
    assertClose(dir[2], -1, 1e-15, 'two quarter turns reverse the direction');
    assertClose(half.o[0], -1, 1e-15, 'the second leg runs along the turned axis');
    assertClose(half.o[2], 1, 1e-15, 'the first leg ran along the original one');
  });

  it('keeps a misalignment local to its own element', () => {
    // A misaligned lens does not move the ones downstream of it: each is
    // mounted independently. Propagating the error would model a bent bench.
    const bl = new Beamline([
      createElement('drift', { length: 10, bore: 6 }),
      createElement('drift', { length: 10, bore: 6 }),
      createElement('drift', { length: 10, bore: 6 }),
    ]);
    const before = bl.elements[2].frame.o.slice();

    bl.elements[1].align = { dx: mmToM(2), dy: 0, tiltX: 0, tiltY: 0 };
    bl.layout();

    assertClose(bl.elements[1].frame.o[0], mmToM(2), 1e-15, 'the element moves');
    for (let i = 0; i < 3; i++) {
      assertClose(
        bl.elements[2].frame.o[i],
        before[i],
        1e-15,
        'the element after it must not move'
      );
    }
    assert(bl.misaligned, 'the column should report itself misaligned');

    bl.autoAlign();
    assertClose(bl.elements[1].frame.o[0], 0, 1e-15, 'auto-align restores the mount');
    assert(!bl.misaligned, 'and the column is aligned again');
  });
});

/* ------------------------------------------------------------------ */
/* bender                                                              */
/* ------------------------------------------------------------------ */

describe('Quadrupole deflector', () => {
  // r0/a = 0.95, the proportions of the reference instrument.
  //
  // No corner posts. Every test in this block is about the IDEAL quadrupole -
  // the closed-form matched voltage, the purity of phi = C X Z, the right
  // angle it produces - and a grounded post on the diagonal is a deliberate
  // departure from that ideal. Its effect belongs in its own tests, below,
  // not smeared through these as an unexplained few per cent. The electrode
  // here is 0.5 mm thick, so a post of any useful size would take most of it.
  // channelWidth is pinned for the same reason every other dimension here is:
  // this block measures a closed form against a specific instrument, so it
  // should not move when a default chosen for a different one moves. The
  // electrode is 0.5 mm thick, so its channel is a 1 mm slot and none of the
  // entrance-deflection behaviour that sets the shipped default applies.
  const PARAMS = {
    apertureRadius: 19,
    electrodeThickness: 0.5,
    boxClearance: 0.5,
    channelWidth: 10.5,
    cornerSize: 0,
    height: 30,
    gridStep: 0.5,
  };
  const ENERGY = 1000;
  const A = mmToM(20); // r0 + thickness + clearance

  /** A column with the deflector in the middle, at `voltage`. */
  const column = (extra = {}) =>
    new Beamline([
      createElement('drift', { length: 8, bore: 6 }),
      createElement('bender', { ...PARAMS, ...extra }),
      createElement('drift', { length: 25, bore: 6 }),
    ]);

  it('matches the derived design constant', () => {
    // V0 = k (T/q)(r0/a)^2 with k = coth^2(s) and cot s = tanh s. Checking the
    // root rather than the constant, because the constant is where a typo
    // would hide.
    const s = DEFLECTOR_DESIGN_PHASE;
    assertClose(Math.cos(s) / Math.sin(s), Math.tanh(s), 1e-6, 'cot s = tanh s');
    assertRelClose(DEFLECTOR_CONSTANT, Math.tanh(s) ** -2, 1e-12, 'k = coth^2 s');

    const V = matchedVoltage(PARAMS, ENERGY, 1);
    assertRelClose(V, DEFLECTOR_CONSTANT * ENERGY * 0.95 ** 2, 1e-12, 'closed form');
  });

  it('scales the matched voltage with energy, charge and aperture ratio', () => {
    const base = matchedVoltage(PARAMS, 1000, 1);
    assertRelClose(matchedVoltage(PARAMS, 2000, 1), 2 * base, 1e-12, 'V ~ T');
    assertRelClose(matchedVoltage(PARAMS, 1000, 2), base / 2, 1e-12, 'V ~ 1/q');

    // The ratio enters squared, which is the part that is easy to get wrong:
    // `a` is the half-width of the field region, not the electrode radius.
    const wider = matchedVoltage(
      { ...PARAMS, electrodeThickness: 4.5, boxClearance: 0.5 },
      1000,
      1
    );
    assertRelClose(wider, DEFLECTOR_CONSTANT * 1000 * (19 / 24) ** 2, 1e-12, 'V ~ (r0/a)^2');
  });

  it('solves a quadrupole potential in the bend plane', () => {
    // phi = C X Z: antisymmetric across each axis, and identically zero ON
    // them. This is the whole reason the device couples the two axes, so it is
    // worth checking directly rather than inferring it from a trajectory.
    const b = createBender({ ...PARAMS, voltage: 100 });
    const r = mmToM(6);
    const at = (X, Z) => b.potentialAt(X, 0, A + Z);

    assertClose(at(0, r), 0, 0.05, 'zero on the entrance axis');
    assertClose(at(r, 0), 0, 0.05, 'zero on the exit axis');
    assertRelClose(at(r, r), -at(-r, r), 1e-6, 'antisymmetric across X');
    assertRelClose(at(r, r), -at(r, -r), 1e-6, 'antisymmetric across Z');
    assertRelClose(at(r, r), at(-r, -r), 1e-6, 'symmetric under a half turn');

    // And it really is bilinear near the centre: halving both coordinates
    // should quarter the potential.
    assertRelClose(at(r / 2, r / 2), at(r, r) / 4, 0.02, 'phi ~ X Z');
  });

  it('solves on a planar grid that is exactly symmetric about the centre', () => {
    // An even node count would put the box half a step further out on one side
    // than the other, and the four electrodes would not see the same
    // enclosure - the deflector would bend by different amounts either way.
    const b = createBender(PARAMS);
    assert(b.grid.nz % 2 === 1 && b.grid.nr % 2 === 1, 'node counts must be odd');
    assertClose(b.grid.z0, -b.grid.step * (b.grid.nz - 1) / 2, 1e-15, 'centred in Z');
    assertClose(b.grid.r0, -b.grid.step * (b.grid.nr - 1) / 2, 1e-15, 'centred in X');
    assert(!b.grid.includesAxis, 'a planar grid has no symmetry axis to claim');
  });

  it('leaves the entrance and exit apertures open', () => {
    // The Laplace problem needs the box closed, but the holes in it are real.
    // Without this the beam is destroyed on the entrance plane.
    const b = createBender(PARAMS);
    assert(!b.strikes(0, 0, 0), 'the entrance is a hole, not a wall');
    assert(!b.strikes(-A, 0, A), 'and so is the exit');
    // The rest of that face is metal.
    assert(b.strikes(mmToM(15), 0, 0), 'the entrance face is otherwise solid');
    assert(b.strikes(-A, 0, A + mmToM(15)), 'and so is the exit face');
  });

  it('turns the column through a right angle', () => {
    // The whole point of the element: its exit faces somewhere else, so
    // everything downstream turns with it.
    const dir = forwardOf(column().exitFrame);
    assertClose(dir[0], -1, 1e-12, 'the column leaves along -x');
    assertClose(dir[2], 0, 1e-12, 'with nothing left along z');
  });

  it('places the element after the bend on the turned axis', () => {
    const bl = new Beamline([
      createElement('bender', PARAMS),
      createElement('drift', { length: 20, bore: 4 }),
    ]);
    const after = bl.elements[1].frame;
    // Entrance at the centre of one face, exit at the centre of the next, so
    // the exit sits one half-width across and one half-width along.
    assertClose(after.o[0], -A, 1e-12, 'exit sits one half-width to the side');
    assertClose(after.o[2], A, 1e-12, 'and one half-width downstream');
    assertClose(forwardOf(after)[0], -1, 1e-12, 'the drift runs across the original axis');
  });

  it('carries a matched ion through and out at ninety degrees', () => {
    const V = matchedVoltage(PARAMS, ENERGY, 1);
    const { points, stop } = flyIon(
      column({ voltage: V }), makeIon({ mass: 100, charge: 1, energy: ENERGY }), { cfl: 0.05 }
    );
    assert(stop === 'exited', `a matched ion should get through, got ${stop}`);

    const last = points[points.length - 1];
    const turned = (Math.atan2(-last.vx, last.vz) * 180) / Math.PI;
    // Two degrees of the nominal right angle. The residual is real and is
    // documented in the element: the electrodes subtend finite arcs and the
    // grounded box pulls the potential down near the apertures, so the solved
    // field is slightly weaker than the ideal form the constant assumes.
    assertClose(turned, 90, 2.5, 'a matched ion turns through a right angle');
  });

  it('has a transmitting window of about a tenth either side', () => {
    // Measured, and relied on by the readout in the UI. Inside the window the
    // ion gets through; well outside it, it does not.
    const V = matchedVoltage(PARAMS, ENERGY, 1);
    const ion = () => makeIon({ mass: 100, charge: 1, energy: ENERGY });
    for (const f of [0.92, 1.0, 1.08]) {
      const { stop } = flyIon(column({ voltage: V * f }), ion(), { cfl: 0.05 });
      assert(stop === 'exited', `${f.toFixed(2)} x matched should transmit, got ${stop}`);
    }
    for (const f of [0.7, 1.4]) {
      const { stop } = flyIon(column({ voltage: V * f }), ion(), { cfl: 0.05 });
      assert(stop !== 'exited', `${f.toFixed(2)} x matched should not transmit, got ${stop}`);
    }
  });

  it('disperses an off-energy ion', () => {
    // A deflector at a fixed voltage selects an energy: this is what makes it
    // an energy filter rather than merely a corner.
    const V = matchedVoltage(PARAMS, ENERGY, 1);
    const bl = column({ voltage: V });
    for (const energy of [ENERGY * 0.7, ENERGY * 1.4]) {
      const { stop } = flyIon(bl, makeIon({ mass: 100, charge: 1, energy }), { cfl: 0.05 });
      assert(
        stop !== 'exited',
        `an ion ${((energy / ENERGY - 1) * 100).toFixed(0)} % off energy should be ` +
          `dispersed, got ${stop}`
      );
    }
  });

  it('keeps an on-orbit ion in the bend plane', () => {
    // The electrodes are uniform perpendicular to the bend plane, so there is
    // no field in that direction and an ion launched in the plane stays in it.
    const V = matchedVoltage(PARAMS, ENERGY, 1);
    const { points } = flyIon(
      column({ voltage: V }), makeIon({ mass: 100, charge: 1, energy: ENERGY }), { cfl: 0.05 }
    );
    for (const p of points) {
      assertClose(p.y, 0, mmToM(1e-3), `left the bend plane at z = ${mToMm(p.z).toFixed(1)} mm`);
    }
  });

  it('bends in whichever plane it is rolled into', () => {
    // One element, one solve, any plane. A 90 degree roll turns the same
    // deflector from a horizontal corner into a vertical one.
    const h = forwardOf(new Beamline([createElement('bender', PARAMS)]).exitFrame);
    const v = forwardOf(
      new Beamline([createElement('bender', { ...PARAMS, bendPlane: 90 })]).exitFrame
    );

    assertClose(h[0], -1, 1e-12, 'an unrolled deflector turns in x');
    assertClose(h[1], 0, 1e-12, 'and not at all in y');
    assertClose(v[1], -1, 1e-12, 'a 90 degree roll turns it in y instead');
    assertClose(v[0], 0, 1e-12, 'and not at all in x');
  });

  it('does not roll the beam it passes on', () => {
    // The conjugation by the roll is what makes this true. Without it a
    // vertical deflector would also rotate "up" into "sideways" for every
    // element downstream, which is not what a bender does.
    const bl = new Beamline([
      createElement('bender', { ...PARAMS, bendPlane: 90 }),
      createElement('drift', { length: 20, bore: 4 }),
    ]);
    const localX = vectorToGlobal(bl.elements[1].frame, [1, 0, 0]);
    assertClose(localX[0], 1, 1e-12, 'transverse x survives a vertical bend');
    assertClose(localX[1], 0, 1e-12, 'with nothing leaking into y');
  });

  it('carries a matched ion round a vertical bend', () => {
    // The same physics, in the other plane. The ion should leave travelling
    // in -y with essentially no x motion at all.
    const V = matchedVoltage(PARAMS, ENERGY, 1);
    const { points, stop } = flyIon(
      column({ bendPlane: 90, voltage: V }),
      makeIon({ mass: 100, charge: 1, energy: ENERGY }),
      { cfl: 0.05 }
    );
    assert(stop === 'exited', `a matched ion should get through, got ${stop}`);

    const last = points[points.length - 1];
    const turned = (Math.atan2(-last.vy, last.vz) * 180) / Math.PI;
    assertClose(turned, 90, 2.5, 'a vertical deflector turns the beam downward');
    for (const p of points) {
      assertClose(p.x, 0, mmToM(1e-3), 'and leaves the horizontal plane alone');
    }
  });

  it('makes the column occupy the vertical plane', () => {
    // Which is what tells the view a side elevation is worth drawing.
    const flat = new Beamline([
      createElement('drift', { length: 20, bore: 4 }),
      createElement('bender', PARAMS),
    ]);
    const tall = new Beamline([
      createElement('drift', { length: 20, bore: 4 }),
      createElement('bender', { ...PARAMS, bendPlane: 90 }),
    ]);
    assert(!flat.usesVerticalPlane, 'a horizontal column stays in its plane');
    assert(tall.usesVerticalPlane, 'a vertical bend leaves it');
  });

  it('bends into any plane, not just the four right angles', () => {
    // The roll is a rotation about the beam applied when the field is read,
    // so an arbitrary angle is no harder than a quarter turn. The plane the
    // beam ends up in must be exactly the plane asked for.
    for (const deg of [0, 30, 45, 90, 135, 180, 270, -45]) {
      const bl = new Beamline([
        createElement('drift', { length: 8, bore: 6 }),
        createElement('bender', { ...PARAMS, bendPlane: deg }),
      ]);
      const d = forwardOf(bl.exitFrame);
      // Azimuth of the exit direction in the transverse plane, measured the
      // same way the roll is.
      const azimuth = (Math.atan2(-d[1], -d[0]) * 180) / Math.PI;
      // Compared modulo a full turn: +180 and -180 are the same plane.
      const off = ((((azimuth - deg) % 360) + 540) % 360) - 180;
      assertClose(off, 0, 1e-9, `a ${deg} degree roll bends into that plane`);
      assertClose(d[2], 0, 1e-9, 'and still turns through a right angle');
    }
  });

  it('transmits a matched ion at a roll no multiple of ninety', () => {
    // The plane being right is not enough; the beam has to get through it.
    const V = matchedVoltage(PARAMS, ENERGY, 1);
    for (const deg of [30, 45, 137]) {
      const { stop } = flyIon(
        column({ bendPlane: deg, voltage: V }),
        makeIon({ mass: 100, charge: 1, energy: ENERGY }),
        { cfl: 0.05 }
      );
      assert(stop === 'exited', `a ${deg} degree roll should transmit, got ${stop}`);
    }
  });

  it('refuses a geometry with no electrode left', () => {
    const refuses = (extra, why) => {
      let threw = false;
      try {
        createBender({ ...PARAMS, ...extra });
      } catch {
        threw = true;
      }
      assert(threw, why);
    };
    // The electrode is what is LEFT between the four channels and inside the
    // outer face, so either can eat it entirely.
    refuses({ channelWidth: 40 }, 'channels wider than the block leave nothing to hold a voltage');
    refuses({ channelWidth: 0 }, 'no channel leaves no way in or out');
    // The electrode here is 0.5 mm thick in a 20 mm box, so the post has to be
    // most of the box before it meets the channel.
    refuses({ cornerSize: 16 }, 'a post that reaches the channels leaves nothing to hold a voltage');
  });

  /*
    The corner posts, and the shielding that makes them free.

    They are grounded structure - tie rods, in a real instrument - standing on
    the diagonals, which is the one direction that crosses no beam channel.

    The first version of these tests asserted that they must WEAKEN the field,
    on the reasoning that they put ground on the diagonal where phi = C X Z is
    largest. That reasoning is wrong, and the suite caught it: the measured
    potential was 37.28 V with posts and 37.28 V without, identical to every
    digit printed. The posts stand BEHIND the electrode blocks, an electrode is
    a conductor, and nothing behind a conductor's surface reaches the field in
    front of it. The aperture is bounded by the arcs and the channel mouths and
    can see nothing else.

    So these tests pin the shielding instead. That is the more valuable thing
    to hold: it is what says the corner geometry is free to match the drawing,
    and it is what identifies the CHANNEL WIDTH as the only part of this shape
    that sets the voltage.
  */
  describe('corner posts', () => {
    // The shipped proportions, where a 5 mm post is a modest bite out of a
    // 20 mm wide electrode rather than most of a thin one. Channel width left
    // at its default, since that is the one dimension under test here that
    // would change the field.
    const THICK = {
      apertureRadius: 19,
      electrodeThickness: 9,
      boxClearance: 1,
      height: 30,
      gridStep: 0.5,
    };
    const CENTRE = mmToM(29); // r0 + thickness + clearance

    /** Potential at a point given in centred bend coordinates, in mm. */
    const phiAt = (b, Xmm, Zmm) => b.potentialAt(mmToM(Xmm), 0, CENTRE + mmToM(Zmm));

    it('puts ground on the diagonal where the electrode would be', () => {
      const bare = createBender({ ...THICK, cornerSize: 0, voltage: 100 });
      const posted = createBender({ ...THICK, cornerSize: 5, voltage: 100 });

      // (26, 26) mm: on the diagonal, outside the 24 mm post boundary, and
      // well inside the 28 mm electrode block. Metal either way - the question
      // is whose.
      assert(
        Math.abs(phiAt(bare, 26, 26)) > 90,
        `without a post the diagonal corner is driven metal, got ${phiAt(bare, 26, 26).toFixed(1)} V`
      );
      assertClose(
        phiAt(posted, 26, 26),
        0,
        1e-9,
        'with a post the diagonal corner is grounded'
      );
    });

    it('leaves the field in the aperture completely unchanged', () => {
      const bare = createBender({ ...THICK, cornerSize: 0, voltage: 100 });
      const posted = createBender({ ...THICK, cornerSize: 5, voltage: 100 });

      // Inside the aperture, off both axes, so the quadrupole term is what is
      // being read - and on the diagonal, the direction the posts stand in and
      // so the least favourable place to look for no effect.
      const a = Math.abs(phiAt(bare, 8, 8));
      const b = Math.abs(phiAt(posted, 8, 8));
      assert(a > 0.5, `expected a usable quadrupole potential, got ${a.toFixed(2)} V`);

      // A thousandth is far below any real effect and far above the noise two
      // SOR solves of different geometries leave behind.
      assertRelClose(
        b,
        a,
        1e-3,
        `the blocks shield the posts, so the aperture field should not move: ` +
          `${b.toFixed(4)} V with, ${a.toFixed(4)} V without`
      );
    });

    it('shields the aperture from the corner whatever stands there', () => {
      /*
        The general statement, of which the posts are one case: the field
        inside is set by the arcs and the channel mouths, so changing the
        corner CANNOT move it. Sweeping the post size is the sharpest way to
        say that - if any of this leaked through, a post four times the size
        would show it.
      */
      const ref = Math.abs(phiAt(createBender({ ...THICK, cornerSize: 0, voltage: 100 }), 8, 8));
      for (const cs of [2, 5, 10]) {
        const v = Math.abs(phiAt(createBender({ ...THICK, cornerSize: cs, voltage: 100 }), 8, 8));
        assertRelClose(v, ref, 1e-3, `a ${cs} mm post moved the aperture field`);
      }
    });

    it('is the channel width, not the corner, that sets the field', () => {
      /*
        The channel is the one part of this shape that reaches the beam, since
        it is the only one that changes the arcs themselves. It weakens the
        field, and WEAKLY - which is the measurement worth holding, because the
        obvious model says otherwise.

        F = (4/pi) cos(2 psi0), psi0 = asin(w/r0), predicts this ratio at 0.762.
        Measured: 0.915. The model treats the gaps as sitting at zero, as a gap
        in a driven boundary would, but a channel is not a gap: its two walls
        belong to neighbouring blocks at OPPOSITE polarity, so the slot goes on
        driving field into the aperture after the arc has ended. Moving the
        walls apart does not remove them, so the field barely notices.

        The bound below is deliberately loose about the value and strict about
        the shape of the answer: it must fall, and it must fall far less than
        the arc-coverage model claims. If a change ever makes 0.762 come true,
        that is the driven walls having stopped working and is worth knowing.
      */
      const narrow = createBender({ ...THICK, channelWidth: 10.5, voltage: 100 });
      const wide = createBender({ ...THICK, channelWidth: 16, voltage: 100 });
      const a = Math.abs(phiAt(narrow, 8, 8));
      const b = Math.abs(phiAt(wide, 8, 8));
      const ratio = b / a;

      assert(
        ratio < 1,
        `a wider channel should weaken the field: ${b.toFixed(2)} V at 16 mm, ` +
          `${a.toFixed(2)} V at 10.5 mm`
      );
      assert(
        ratio > 0.85,
        `the driven channel walls should keep the loss small; measured ${ratio.toFixed(3)}, ` +
          'against 0.762 for a boundary whose gaps really were grounded'
      );
    });

    it('still turns the beam through a right angle', () => {
      /*
        The posts change the voltage needed, not what the device does. The
        tuner is what finds the new voltage; this only asks that a right angle
        is still reachable.

        The sweep is deliberately wide. Neither the posts nor the filled
        corners move the aperture field at all - the blocks shield both - so
        the expectation is that a right angle lands near the ideal, and the
        width is there to catch the case where that reasoning is wrong again
        rather than to encode confidence in it. The failure message reports
        where the right angle actually landed, so a run of this suite MEASURES
        the calibration instead of assuming it.
      */
      const V0 = matchedVoltage({ ...THICK, cornerSize: 5 }, 1000, 1);
      let best = null;
      for (let f = 0.3; f <= 1.9; f += 0.05) {
        const bl = new Beamline([
          createElement('drift', { length: 8, bore: 6 }),
          createElement('bender', { ...THICK, cornerSize: 5, voltage: f * V0 }),
          createElement('drift', { length: 25, bore: 6 }),
        ]);
        const { stop, points } = flyIon(
          bl,
          makeIon({ mass: 100, charge: 1, energy: 1000 }),
          { cfl: 0.05 }
        );
        if (stop !== 'exited' || points.length < 6) continue;
        const p1 = points[points.length - 1];
        const p0 = points[points.length - 6];
        const deg = (Math.atan2(-(p1.x - p0.x), p1.z - p0.z) * 180) / Math.PI;
        const err = Math.abs(deg - 90);
        if (!best || err < best.err) best = { err, f, deg };
      }
      assert(best, 'no voltage between 0.3 and 1.9 V0 transmitted the ion at all');
      assert(
        best.err < 5,
        `best turn was ${best.deg.toFixed(1)} deg at V/V0 = ${best.f.toFixed(2)} ` +
          `(V0 = ${V0.toFixed(0)} V) - if that ratio is far from 1, suspect the arc ` +
          'coverage: F = (4/pi) cos(2 asin(w/r0)) says the channel width sets it'
      );
    });
  });
});

/* ------------------------------------------------------------------ */
/* starting values                                                     */
/* ------------------------------------------------------------------ */

describe('Elements as placed from the toolbar', () => {
  // The ions a user is plausibly simulating, including a negative one and a
  // doubly charged one - the two cases a voltage written for "an ion" gets
  // wrong.
  const IONS = [
    { mass: 4, charge: 1, energy: 10 },
    { mass: 100, charge: 1, energy: 50 },
    { mass: 100, charge: 2, energy: 500 },
    { mass: 1000, charge: 1, energy: 2000 },
    { mass: 100, charge: -1, energy: 50 },
  ];

  const describeIon = (i) => `${i.mass} u, ${i.charge > 0 ? '+' : ''}${i.charge}, ${i.energy} eV`;

  /**
   * How many ions an element passes on, placing it exactly as the toolbar
   * does.
   *
   * A detector is scored on what it COLLECTS rather than on what it passes on,
   * because passing ions on is precisely what it is not for. Counting it the
   * same way as the rest would demand that a working detector let the beam
   * through it.
   */
  function alone(type, ion) {
    const bl = new Beamline([
      createElement('drift', { length: 12, bore: 5 }),
      createElement(type, startingParams(type, ion)),
      createElement('drift', { length: 25, bore: 5 }),
    ]);
    const { tracks } = flyBeam(bl, discBeam({ ...ion, count: 9, radius: 1.0 }), {
      cfl: 0.05,
      maxSteps: 600000,
    });
    if (type === 'detector') {
      return tracks.filter((t) => {
        const p = t.points[t.points.length - 1];
        return t.stop === 'electrode' && bl.detected(p.x, p.y ?? 0, p.z);
      }).length;
    }
    return tracks.filter((t) => t.stop === 'exited').length;
  }

  it('transmits the beam for every element type', () => {
    // The bug this guards against is not subtle and was shipped: a deflector
    // placed from the toolbar arrived at zero volts, so it did not deflect at
    // all - the beam flew straight on into the far wall of the box. An einzel
    // arrived at a voltage forty times too strong for the default beam, and a
    // quadrupole at a Mathieu q of 1.83, outside the first stability region
    // entirely. Three of the five elements did nothing useful when placed.
    const ion = { mass: 100, charge: 1, energy: 50 };
    for (const type of Object.keys(ELEMENT_TYPES)) {
      const n = alone(type, ion);
      assert(n >= 7, `a ${type} placed from the toolbar transmitted only ${n} of 9`);
    }
  });

  it('follows the ion in the source, not the one the defaults were written for', () => {
    // Every starting value that is set at all is set from the beam. A fixed
    // number cannot be right for both a 4 u ion at 10 eV and a 1000 u ion at
    // 2 keV, and the scalings that make it right are physics: a lens depends
    // only on T/q, a deflector on T/q and its own proportions, a filter on
    // mass through the Mathieu q.
    for (const ion of IONS) {
      for (const type of ['einzel', 'bender', 'quadrupole']) {
        const n = alone(type, ion);
        assert(n === 9, `${type} lost ${9 - n} of 9 for ${describeIon(ion)}`);
      }
    }
  });

  it('flips polarity for a negative ion', () => {
    // The case a magnitude gets wrong. A negative ion meeting the centre
    // electrode of a lens set for a positive one is decelerated rather than
    // accelerated, and cannot climb the barrier; a deflector at the wrong
    // polarity steers it into the wall instead of round the corner.
    const pos = { mass: 100, charge: 1, energy: 50 };
    const neg = { mass: 100, charge: -1, energy: 50 };
    for (const type of ['einzel', 'bender']) {
      const a = startingParams(type, pos).voltage;
      const b = startingParams(type, neg).voltage;
      assertRelClose(b, -a, 1e-12, `${type} reverses for a negative ion`);
    }
    // The RF filter does not: it confines either sign, because the drive
    // reverses every half cycle regardless.
    assertRelClose(
      startingParams('quadrupole', neg).rfAmplitude,
      startingParams('quadrupole', pos).rfAmplitude,
      1e-12,
      'an RF filter has no polarity to get wrong'
    );
  });

  it('puts a placed quadrupole at a usable Mathieu q', () => {
    for (const ion of IONS) {
      const q = createElement('quadrupole', startingParams('quadrupole', ion));
      const m = q.mathieu(ion.mass, Math.abs(ion.charge));
      // One per cent, not exact: the starting amplitude is snapped to the
      // step its own control will show, so that the box does not read one
      // number while the element holds another. Half a step is a few tenths
      // of a per cent of the value, and q follows it.
      assertRelClose(m.q, MATHIEU_Q_WORKING, 0.01, `q for ${describeIon(ion)}`);
      assert(m.q < MATHIEU_Q_LIMIT, 'and inside the first stability region');
    }
  });

  it('inverts the Mathieu relation exactly', () => {
    // amplitudeForQ is the inverse of mathieu, and a round trip is the only
    // check that cannot be fooled by both being wrong the same way... so it
    // is checked against the closed form too.
    const params = { gridStep: 0.5, fieldRadius: 4, frequency: 2 };
    const V = amplitudeForQ(params, 100, 1, 0.5);
    const q = createElement('quadrupole', { ...params, rfAmplitude: V });
    assertRelClose(q.mathieu(100, 1).q, 0.5, 1e-12, 'round trip');

    // q = 4 z e V / (m r0^2 Omega^2), written out independently.
    const omega = 2 * Math.PI * 2e6;
    const expected =
      (0.5 * 100 * ATOMIC_MASS_UNIT * mmToM(4) ** 2 * omega ** 2) / (4 * ELEMENTARY_CHARGE);
    assertRelClose(V, expected, 1e-12, 'closed form');
  });

  it('keeps a scaled control inside the parameter’s hard limits', () => {
    // The slider may be narrowed to where the answer is, but it can never
    // offer a setting the parameter does not allow - or the UI would propose
    // voltages the element rejects.
    const field = ELEMENT_TYPES.bender.fields.find((f) => f.key === 'voltage');
    for (const energy of [1, 50, 1000, 500000]) {
      const ion = { mass: 100, charge: 1, energy };
      const r = fieldRange(field, { ...ELEMENT_TYPES.bender.defaults }, ion);
      assert(r.min >= field.min && r.max <= field.max, `${energy} eV stays within the limits`);
      assert(r.step > 0, 'and has a usable step');
    }
  });

  it('widens a scaled control to contain the value already set', () => {
    // Otherwise changing the beam energy would re-render the panel with a
    // range that excludes the current voltage, and the browser would clamp
    // it - silently moving a setting the user chose deliberately.
    const field = ELEMENT_TYPES.bender.fields.find((f) => f.key === 'voltage');
    const ion = { mass: 100, charge: 1, energy: 50 }; // matched near 40 V
    const r = fieldRange(field, { ...ELEMENT_TYPES.bender.defaults, voltage: 5000 }, ion);
    assert(r.max >= 5000, `range should reach the 5000 V already set, stops at ${r.max}`);
  });
});

/* ------------------------------------------------------------------ */
/* space charge in a drift                                             */
/* ------------------------------------------------------------------ */

describe('A beam expanding under its own charge', () => {
  const SPEC = { mass: 100, charge: 1, energy: 50 };
  const R0 = 2.0; // mm
  const LEN = 100; // mm
  const speed = Math.sqrt(
    (2 * SPEC.energy * ELEMENTARY_CHARGE) / (SPEC.mass * ATOMIC_MASS_UNIT)
  );

  /**
   * The envelope equation for a uniform cylindrical beam,
   *
   *     r'' = K / r,   K = q I / (2 pi eps0 m v^3),
   *
   * marched from a waist. Derived here rather than taken from the simulator,
   * so it is an independent prediction of how far the beam should have spread
   * and not a restatement of what the code already does.
   */
  function envelope(I, lenMm) {
    const K =
      (SPEC.charge * ELEMENTARY_CHARGE * I) /
      (2 * Math.PI * VACUUM_PERMITTIVITY * SPEC.mass * ATOMIC_MASS_UNIT * speed ** 3);
    let r = mmToM(R0);
    let rp = 0;
    const dz = 1e-6;
    for (let z = 0; z < mmToM(lenMm); z += dz) {
      rp += (K / r) * dz;
      r += rp * dz;
    }
    return r;
  }

  /** Widest ion at the end of a plain drift, bore wide enough to touch nothing. */
  function widthAfter(opts) {
    const bl = new Beamline([createElement('drift', { length: LEN, bore: 60 })]);
    const { tracks } = flyBeam(bl, discBeam({ ...SPEC, count: 41, radius: R0 }), {
      cfl: 0.02,
      maxSteps: 800000,
      ...opts,
    });
    let r = 0;
    for (const t of tracks) {
      const p = t.points[t.points.length - 1];
      r = Math.max(r, Math.hypot(p.x, p.y ?? 0));
    }
    return r;
  }

  it('spreads by the amount the envelope equation predicts', () => {
    // The commonest surprise in the simulator is a beam that sails down a
    // drift without widening, and the answer is that repulsion defaults to
    // off. This pins down that it is off by CHOICE and not by omission: with
    // it on, the spreading is right to a fraction of a per cent.
    for (const uA of [0.2, 0.5, 1, 2]) {
      const got = widthAfter({ repulsion: 'beam', beamCurrent: uA * 1e-6 });
      assertRelClose(
        got,
        envelope(uA * 1e-6, LEN),
        0.02,
        `${uA} uA over ${LEN} mm of drift`
      );
      assert(got > mmToM(R0) * 2, 'and it really has expanded');
    }
  });

  it('does nothing at all when it is off', () => {
    // Off has to mean off. A collimated beam in a field-free tube must arrive
    // exactly as wide as it started.
    assertRelClose(widthAfter({}), widthAfter({ repulsion: 'none' }), 1e-12, 'same either way');
    assertRelClose(widthAfter({}), mmToM(R0), 0.01, 'and unchanged from launch');
  });

  it('spreads more with more current', () => {
    const a = widthAfter({ repulsion: 'beam', beamCurrent: 0.5e-6 });
    const b = widthAfter({ repulsion: 'beam', beamCurrent: 2e-6 });
    assert(b > a, `2 uA should spread further than 0.5: ${mToMm(a)} vs ${mToMm(b)} mm`);
  });
});

/* ------------------------------------------------------------------ */
/* multipole guides                                                    */
/* ------------------------------------------------------------------ */

describe('Multipole guide', () => {
  const SPEC = { mass: 100, charge: 1, energy: 5 };

  /** The largest rod that fits between 2n of them around r0, at 60 %. */
  const fitRod = (poles, r0 = 5) => {
    const s = Math.sin(Math.PI / poles);
    return Math.min(Math.round(((r0 * s) / (1 - s)) * 0.6 * 10) / 10, 3);
  };

  const build = (poles, over = {}) =>
    createElement('multipole', {
      poles,
      fieldRadius: 5,
      rodRadius: fitRod(poles),
      length: 120,
      rfAmplitude: 300,
      frequency: 2,
      gridStep: 0.25,
      ...over,
    });

  /** Field magnitude averaged round a circle of radius r. */
  const ring = (e, r) => {
    let sum = 0;
    const n = 32;
    for (let k = 0; k < n; k++) {
      const a = (2 * Math.PI * k) / n;
      const f = e.field.fieldAt(r * Math.cos(a), r * Math.sin(a));
      sum += Math.hypot(f.Ez, f.Er);
    }
    return sum / n;
  };

  it('produces a field going as r^(n-1), the ideal multipole form', () => {
    // The defining property. An ideal 2n-pole has phi ~ r^n cos(n theta), so
    // |E| ~ r^(n-1): linear in r for a quadrupole, cubic for an octopole.
    // Fitted from the SOLVED field of real round rods, near the rim - for a
    // 12-pole the multipole term at a third of the field radius is a quarter
    // of a per cent of the rim field, and anything else present swamps it.
    for (const poles of [4, 6, 8, 12]) {
      const e = build(poles);
      const k =
        Math.log(ring(e, mmToM(4.5)) / ring(e, mmToM(3))) / Math.log(4.5 / 3);
      assertClose(k, poles / 2 - 1, 0.1, `${poles} rods should give |E| ~ r^${poles / 2 - 1}`);
    }
  });

  it('deepens its well as the pole count rises', () => {
    // At fixed voltage and radius, a flatter-bottomed well is also a deeper
    // one, which is why octopoles are used to hold ions and quadrupoles to
    // select them.
    let last = 0;
    for (const poles of [4, 6, 8]) {
      const { depth } = build(poles).trapping(SPEC.mass, SPEC.charge);
      assert(depth > last, `${poles} rods should trap harder than fewer`);
      last = depth;
    }
  });

  it('matches the closed-form well depth for an ideal quadrupole', () => {
    // Dehmelt: U* = q^2 E0^2 / (4 m Omega^2), with |E0| = 2V/r0 at the field
    // radius of an ideal quadrupole. Round rods at the classic 1.1487 r0 are
    // not hyperbolic, so agreement is close rather than exact - and which way
    // it misses is the point: the solved field is the weaker one.
    const e = build(4, { rodRadius: 5.74, housingRadius: 22, gridStep: 0.25 });
    const { depth } = e.trapping(100, 1);
    const r0 = mmToM(5);
    const w = 2 * Math.PI * 2e6;
    const E0 = (2 * 300) / r0;
    const ideal =
      (ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * E0 * E0) /
      (4 * 100 * ATOMIC_MASS_UNIT * w * w) /
      ELEMENTARY_CHARGE;
    assertRelClose(depth, ideal, 0.15, 'within fifteen per cent of the ideal form');
    assert(depth < ideal, 'and round rods give a weaker field than hyperbolic ones');
  });

  it('guides better with more rods', () => {
    // The reason the element exists. Same voltage, same aperture, same beam.
    const through = (poles) => {
      const bl = new Beamline([
        createElement('drift', { length: 10, bore: 4 }),
        build(poles),
        createElement('drift', { length: 15, bore: 4 }),
      ]);
      const { tracks } = flyBeam(bl, discBeam({ ...SPEC, count: 9, radius: 1.5 }), {
        cfl: 0.05,
        maxSteps: 600000,
      });
      return tracks.filter((t) => t.stop === 'exited').length;
    };
    const four = through(4);
    const eight = through(8);
    assert(eight >= four, `eight rods should guide at least as well: ${four} vs ${eight}`);
    assert(eight >= 8, `an octopole should pass nearly everything, got ${eight}/9`);
  });

  it('refuses rods that will not fit', () => {
    // Adjacent rods touching would short the two phases together, which is not
    // a multipole at all.
    let threw = false;
    try {
      createElement('multipole', { poles: 12, fieldRadius: 5, rodRadius: 3 });
    } catch {
      threw = true;
    }
    assert(threw, 'twelve fat rods around a small aperture is not a geometry');
  });

  it('flags a drive too slow for an effective potential to mean anything', () => {
    // The whole idea rests on the drive being fast compared with the ion's
    // motion. Past a Mathieu q of about 0.3 it is not, and a depth quoted from
    // it would be a number about a model that no longer applies.
    const slow = build(8, { frequency: 0.3 });
    assert(!slow.trapping(SPEC.mass, SPEC.charge).valid, 'a slow drive is flagged');
    assert(build(8, { frequency: 4 }).trapping(SPEC.mass, SPEC.charge).valid, 'a fast one is not');
  });
});

/* ------------------------------------------------------------------ */
/* ion funnels                                                         */
/* ------------------------------------------------------------------ */

describe('Ion funnel', () => {
  const SPEC = { mass: 100, charge: 1, energy: 5 };

  it('narrows its aperture down the stack', () => {
    const f = createElement('funnel', { entryRadius: 12, exitRadius: 2, rings: 14 });
    assertRelClose(f.apertureOf(0), mmToM(12), 1e-9, 'first ring');
    assertRelClose(f.apertureOf(13), mmToM(2), 1e-9, 'last ring');
    for (let k = 1; k < 14; k++) {
      assert(f.apertureOf(k) < f.apertureOf(k - 1), `ring ${k} is narrower than the one before`);
    }
  });

  it('carries a monotonic DC ramp along its axis', () => {
    /*
      The failure this guards against was total and did not look like a bug.
      With the domain's end faces grounded, a DC ramp painted on the rings
      becomes a HILL: zero at the entrance, up to 13 V inside, back to zero at
      the exit. A 5 eV beam is reflected in its entirety - 0 of 9 through at
      every width - and nothing about that says "boundary condition". The end
      faces belong to the ring beside them instead, which is what a funnel
      embedded in a longer line actually looks like.
    */
    const f = createElement('funnel', { dcEntry: 20, dcExit: 0 });
    const zs = [0, 0.25, 0.5, 0.75, 1].map((t) => t * f.length);
    const phi = zs.map((z) => f.dcField.potentialAt3D(0, 0, z));
    assertClose(phi[0], 20, 1.5, 'starts at the entry potential');
    assertClose(phi[phi.length - 1], 0, 1.5, 'ends at the exit potential');
    for (let i = 1; i < phi.length; i++) {
      assert(phi[i] < phi[i - 1], `the ramp must fall all the way: ${phi.join(' -> ')}`);
    }
  });

  it('keeps the DC and RF patterns independent over one solve', () => {
    // Two weighted readings of the same basis. Changing the amplitude must not
    // touch the ramp, and it must not re-solve.
    const f = createElement('funnel', {});
    const grid = f.grid;
    const mid = f.length / 2;
    const dcBefore = f.dcField.potentialAt3D(0, 0, mid);
    const rfBefore = f.field.potentialAt3D(0, 0, mmToM(1));

    f.setVoltage(400);
    assert(f.grid === grid, 'the same grid, not a fresh solve');
    assertClose(f.dcField.potentialAt3D(0, 0, mid), dcBefore, 1e-9, 'the ramp is untouched');
    // The RF map is held at unit amplitude and scaled by the drive, so it too
    // is unchanged - the amplitude lives in the drive, not in the map.
    assertClose(f.field.potentialAt3D(0, 0, mmToM(1)), rfBefore, 1e-9, 'the RF map is untouched');
  });

  it('drives the rings in alternating phase', () => {
    // Adjacent rings at the same potential would make the stack a tube, not a
    // funnel: there would be no RF wall at all.
    const f = createElement('funnel', { rings: 8 });
    for (let k = 0; k < 8; k++) {
      assertClose(f.field.voltageOf(`ring${k}`), k % 2 === 0 ? 1 : -1, 1e-12, `ring ${k}`);
    }

    /*
      And it reaches the beam - measured just inside each aperture, where the
      wall the ions are held off actually is.

      Not on the axis. Near the wide end the rings are far closer together
      than the aperture is across, so their opposite phases very nearly cancel
      before they get there and the axial RF is almost nothing. That is not a
      defect, it is the reason a funnel has to narrow: the wall only closes in
      on the beam once the aperture becomes comparable to the pitch.
    */
    const rim = (k) => {
      const r = f.rects[k];
      return f.field.fieldAt((r.z0 + r.z1) / 2, Math.max(r.r0 - f.grid.step * 2, f.grid.step));
    };
    for (let k = 1; k < 8; k++) {
      const a = rim(k - 1);
      const b = rim(k);
      assert(
        a.Er * b.Er < 0,
        `the RF wall must reverse from ring ${k - 1} to ring ${k}`
      );
    }
  });

  it('transmits a cold beam at its defaults', () => {
    // Which the defaults were chosen to do. They are NOT the settings a real
    // funnel runs at - see the note in the element - because without gas to
    // carry the RF heating away, a deeper wall transmits worse.
    const bl = new Beamline([
      createElement('drift', { length: 8, bore: 12 }),
      createElement('funnel', {}),
      createElement('drift', { length: 20, bore: 3 }),
    ]);
    const { tracks } = flyBeam(bl, discBeam({ ...SPEC, count: 9, radius: 1.5 }), {
      cfl: 0.05,
      maxSteps: 1200000,
    });
    const n = tracks.filter((t) => t.stop === 'exited').length;
    assert(n >= 8, `expected a cold beam through, got ${n}/9`);
  });

  it('refuses rings that would touch', () => {
    let threw = false;
    try {
      createElement('funnel', { ringThickness: 3, pitch: 2 });
    } catch {
      threw = true;
    }
    assert(threw, 'rings thicker than their pitch are not a stack');
  });
});

/* ------------------------------------------------------------------ */
/* detector                                                            */
/* ------------------------------------------------------------------ */

describe('Ion detector', () => {
  const SPEC = { mass: 100, charge: 1, energy: 50 };

  it('pulls ions in rather than waiting to be hit', () => {
    // The bias reaches out through the aperture, which is where a detector's
    // collection efficiency comes from. An ion aimed to MISS the collecting
    // surface should be bent onto it.
    const aimed = (volts) => {
      const bl = new Beamline([
        createElement('drift', { length: 20, bore: 10 }),
        createElement('detector', { voltage: volts, entranceRadius: 8, activeRadius: 8 }),
      ]);
      const { tracks } = flyBeam(bl, discBeam({ ...SPEC, count: 9, radius: 5 }), {
        cfl: 0.05,
        maxSteps: 400000,
      });
      return tracks.filter((t) => {
        const p = t.points[t.points.length - 1];
        return t.stop === 'electrode' && bl.detected(p.x, p.y ?? 0, p.z);
      }).length;
    };
    const off = aimed(0);
    const on = aimed(-3000);
    assert(on >= off, `bias should not collect less: ${off} -> ${on}`);
    assert(on >= 8, `a -3 kV detector should collect nearly everything, got ${on}/9`);
  });

  it('tells its collecting surface from its housing', () => {
    // Both are strikes; only one is a count, and the difference is the whole
    // reason the element exists.
    const d = createElement('detector', {});
    const back = d.length;
    assert(d.detects(0, 0, back), 'the middle of the surface counts');
    assert(!d.detects(mmToM(13), 0, back * 0.2), 'the housing wall does not');
    assert(!d.detects(0, 0, back * 0.3), 'and nor does empty space in front of it');
  });

  it('accelerates what it collects', () => {
    // Which matters for a real detector: secondary-electron yield depends on
    // impact energy, so a collector at a few kilovolts is not just a target.
    const bl = new Beamline([
      createElement('drift', { length: 20, bore: 10 }),
      createElement('detector', { voltage: -3000 }),
    ]);
    const { tracks } = flyBeam(bl, [makeIon({ ...SPEC, x: 0 })], { cfl: 0.05, maxSteps: 400000 });
    const t = tracks[0];
    const gained = joulesToEV(kineticEnergy(t.points[t.points.length - 1])) - SPEC.energy;
    assertRelClose(gained, 3000, 0.05, 'an ion lands with the bias added to its energy');
  });
});

/* ------------------------------------------------------------------ */
/* saving and loading                                                  */
/* ------------------------------------------------------------------ */

describe('Saving a column', () => {
  const SPEC = { mass: 100, charge: 1, energy: 50 };

  /** A branching column with something set on every branch. */
  function built() {
    const bl = new Beamline([
      createElement('drift', { length: 12, bore: 5 }),
      createElement('bender', { voltage: startingParams('bender', SPEC).voltage }),
    ]);
    const sw = bl.elements[1];
    const a = bl.add(createElement('drift', { length: 30, bore: 5 }), 2, {
      parent: sw, port: 'bend',
    });
    bl.add(createElement('einzel', { gridStep: 1, voltage: -275 }), 3, {
      parent: a, port: 'out',
    });
    bl.add(createElement('detector', { voltage: -2500, gridStep: 0.6 }), 4, {
      parent: sw, port: 'straight',
    });
    bl.elements[0].align = { dx: mmToM(0.35), dy: mmToM(-0.2), tiltX: 0, tiltY: 0 };
    return bl;
  }

  it('round-trips the tree, the settings and the alignment', () => {
    const before = built();
    const data = JSON.parse(JSON.stringify(serialise(before, { source: SPEC })));
    const { elements, problems } = restore(data);
    assert(problems.length === 0, `expected a clean load, got: ${problems.join(' ')}`);

    const after = new Beamline().adopt(elements);
    assert(after.elements.length === before.elements.length, 'same number of elements');

    for (let i = 0; i < before.elements.length; i++) {
      const a = before.elements[i];
      const b = after.elements[i];
      assert(a.typeKey === b.typeKey, `element ${i} is the same kind`);
      // The port it hangs from, and which element that is, by position.
      assert(a.from.port === b.from.port, `element ${i} keeps its port`);
      const pa = a.from.parent ? before.elements.indexOf(a.from.parent) : null;
      const pb = b.from.parent ? after.elements.indexOf(b.from.parent) : null;
      assert(pa === pb, `element ${i} keeps its parent`);
      assertClose(a.align.dx, b.align.dx, 1e-15, `element ${i} keeps its offset`);
      assertClose(a.align.dy, b.align.dy, 1e-15, `element ${i} keeps its offset`);
    }
  });

  it('rebuilds the same fields, not just the same numbers', () => {
    // A file describes a column; loading it re-solves. So the reloaded column
    // must fly the beam to the same place, not merely list the same settings.
    const before = built();
    const after = new Beamline().adopt(
      restore(JSON.parse(JSON.stringify(serialise(before, {})))).elements
    );
    const ions = () => discBeam({ ...SPEC, count: 9, radius: 1 });
    const fly = (bl) =>
      flyBeam(bl, ions(), { cfl: 0.05, maxSteps: 400000 }).tracks.map((t) => ({
        stop: t.stop,
        end: t.points[t.points.length - 1],
      }));

    const a = fly(before);
    const b = fly(after);
    for (let i = 0; i < a.length; i++) {
      assert(a[i].stop === b[i].stop, `ion ${i} meets the same fate`);
      for (const k of ['x', 'y', 'z']) {
        assertClose(a[i].end[k], b[i].end[k], 1e-12, `ion ${i} lands in the same place`);
      }
    }
  });

  it('survives a parameter that did not exist when it was written', () => {
    // Settings are merged over today's defaults, so an older file loads with
    // the new parameter at its default rather than as undefined.
    const bl = built();
    const data = serialise(bl, {});
    for (const e of data.elements) delete e.params.gridStep;
    const { elements, problems } = restore(data);
    assert(problems.length === 0, `should still load: ${problems.join(' ')}`);
    // Only the elements that have one - a drift has no grid to step over.
    const solved = elements.filter((e) => 'gridStep' in ELEMENT_TYPES[e.typeKey].defaults);
    assert(solved.length > 0, 'the column has elements with a grid step');
    assert(
      solved.every((e) => Number.isFinite(e.params.gridStep)),
      'the missing parameter takes its default'
    );
  });

  it('refuses a file that is not a column, and says why', () => {
    for (const bad of [null, {}, { format: 'something-else', elements: [] }]) {
      const { elements, problems } = restore(bad);
      assert(elements.length === 0, 'nothing is built');
      assert(problems.length > 0, 'and there is a reason given');
    }
    // A parent that is not in the file, rather than a silently broken tree.
    const { problems } = restore({
      format: 'iontrace.column',
      version: 1,
      elements: [{ type: 'drift', params: {}, parent: 7, port: 'out', align: {} }],
    });
    assert(problems.length > 0, 'a dangling parent is reported');
  });

  it('moves an element whose exit no longer exists rather than dropping it', () => {
    const data = serialise(built(), {});
    const onStraight = data.elements.findIndex((e) => e.port === 'straight');
    data.elements[onStraight].port = 'sideways'; // a port from some other version
    const { elements, problems } = restore(data);
    assert(elements.length === data.elements.length, 'nothing is lost');
    assert(
      problems.some((p) => /no longer exists/.test(p)),
      `the change is reported, got: ${problems.join(' ')}`
    );
    assert(elements[onStraight].from.port === 'bend', 'and it lands on a real exit');
  });
});

/* ------------------------------------------------------------------ */
/* branching                                                           */
/* ------------------------------------------------------------------ */

describe('Branching columns', () => {
  const SPEC = { mass: 100, charge: 1, energy: 50 };
  const V = () => startingParams('bender', SPEC).voltage;

  /** Source, deflector, and a line on each of its exits. */
  function switched() {
    const bl = new Beamline([
      createElement('drift', { length: 12, bore: 5 }),
      createElement('bender', { voltage: V() }),
    ]);
    const sw = bl.elements[1];
    bl.add(createElement('drift', { length: 30, bore: 5 }), 2, { parent: sw, port: 'bend' });
    bl.add(createElement('drift', { length: 40, bore: 5 }), 3, { parent: sw, port: 'straight' });
    return { bl, sw };
  }

  it('gives a deflector three exits and a plain element one', () => {
    // The box has an aperture on each of its four faces, because the
    // electrodes stop short of both axes in both directions. One is the
    // entrance; the other three are all reachable, and the voltage alone picks
    // which - positive bends one way, negative the other, zero goes straight.
    const b = createElement('bender', { voltage: V() });
    assert(exitsOf(b).length === 3, `a deflector has three ways out, got ${exitsOf(b).length}`);
    assert(
      exitsOf(b).map((x) => x.port).join(',') === 'bend,straight,counter',
      'bent first, because that is the reason it is in the line'
    );
    assert(exitsOf(createElement('drift', {})).length === 1, 'a drift has one');
  });

  it('points the two bends opposite ways in the same plane', () => {
    const bl = new Beamline([
      createElement('drift', { length: 12, bore: 5 }),
      createElement('bender', { voltage: V() }),
    ]);
    const sw = bl.elements[1];
    const dirOf = (port) => {
      const exit = exitsOf(sw).find((x) => x.port === port);
      return forwardOf(compose(sw.nominalFrame, exit.transform));
    };
    const a = dirOf('bend');
    const c = dirOf('counter');
    const s = dirOf('straight');

    assertClose(a[0], -1, 1e-12, 'one bend goes -x');
    assertClose(c[0], 1, 1e-12, 'the other goes +x');
    assertClose(s[2], 1, 1e-12, 'and straight carries on down the axis');
    // Exactly opposed, which is what "the other way" has to mean.
    const dot = a[0] * c[0] + a[1] * c[1] + a[2] * c[2];
    assertClose(dot, -1, 1e-12, 'the two bends are exactly opposed');
  });

  it('rolls both bends together', () => {
    // The roll sets the PLANE; the polarity picks the direction within it. A
    // quarter turn should take the pair from left/right to down/up.
    const bl = new Beamline([
      createElement('drift', { length: 12, bore: 5 }),
      createElement('bender', { voltage: V(), bendPlane: 90 }),
    ]);
    const sw = bl.elements[1];
    const dirOf = (port) =>
      forwardOf(compose(sw.nominalFrame, exitsOf(sw).find((x) => x.port === port).transform));
    assertClose(dirOf('bend')[1], -1, 1e-12, 'one bend goes -y');
    assertClose(dirOf('counter')[1], 1, 1e-12, 'the other goes +y');
    assertClose(dirOf('straight')[2], 1, 1e-12, 'straight is unaffected by the roll');
  });

  it('sends the beam to whichever port the polarity chooses', () => {
    // The point of a three-way switch: one knob, three destinations, no
    // hardware moved. Reversing the polarity mirrors the whole problem in x,
    // so the counter-bend needs no second matched voltage to be found.
    const bl = new Beamline([
      createElement('drift', { length: 12, bore: 5 }),
      createElement('bender', { voltage: V() }),
    ]);
    const sw = bl.elements[1];
    const ends = {};
    for (const x of exitsOf(sw)) {
      ends[x.port] = bl.add(createElement('drift', { length: 35, bore: 5 }), bl.elements.length, {
        parent: sw,
        port: x.port,
      });
    }

    const destination = (volts) => {
      sw.setVoltage(volts);
      const { tracks } = flyBeam(bl, discBeam({ ...SPEC, count: 9, radius: 1 }), {
        cfl: 0.05,
        maxSteps: 400000,
      });
      const tally = {};
      for (const t of tracks) {
        if (t.stop !== 'exited') continue;
        const p = t.points[t.points.length - 1];
        const e = bl.endNearest(p.x, p.y ?? 0, p.z);
        const port = Object.keys(ends).find((k) => ends[k] === e?.element) ?? '?';
        tally[port] = (tally[port] ?? 0) + 1;
      }
      return tally;
    };

    assert(destination(V())?.bend === 9, 'a positive matched voltage bends one way');
    assert(destination(0)?.straight === 9, 'zero goes straight through');
    assert(destination(-V())?.counter === 9, 'a negative one bends the other way');
  });

  it('places each branch where its own exit points', () => {
    const { bl, sw } = switched();
    const bent = bl.childAt(sw, 'bend');
    const straight = bl.childAt(sw, 'straight');

    // The bend leaves sideways; the straight path carries on down the axis.
    assertClose(forwardOf(bent.frame)[0], -1, 1e-12, 'the bent branch turns');
    assertClose(forwardOf(straight.frame)[2], 1, 1e-12, 'the straight one does not');
    assertClose(straight.frame.o[0], 0, 1e-12, 'and stays on the original axis');
    // Straight through the box is 2a; round the bend is a quarter arc of a.
    const a = mmToM(19 + 9 + 1);
    assertRelClose(straight.zStart - sw.zStart, 2 * a, 1e-9, 'straight path length');
    assertRelClose(bent.zStart - sw.zStart, (Math.PI / 2) * a, 1e-9, 'bent path length');
  });

  it('reports an open end for every unused exit', () => {
    const bl = new Beamline([
      createElement('drift', { length: 12, bore: 5 }),
      createElement('bender', { voltage: V() }),
    ]);
    const ends = bl.openEnds();
    assert(ends.length === 3, `a bare deflector has three open ends, got ${ends.length}`);
    // `switched()` builds lines on two of the three, so the third stays open.
    const { bl: full } = switched();
    assert(full.openEnds().length === 3, 'two lines and one bare port is still three ends');
    assert(
      full.openEnds().filter((e) => e.element.typeKey === 'drift').length === 2,
      'two of them now sit at the ends of those lines'
    );
  });

  it('switches the beam with the voltage, not with the topology', () => {
    // The point of the whole thing. The hardware does not move; the field
    // decides which way the ions go.
    const { bl, sw } = switched();
    const where = () => {
      const { tracks } = flyBeam(bl, discBeam({ ...SPEC, count: 9, radius: 1 }), {
        cfl: 0.05,
        maxSteps: 400000,
      });
      const out = tracks.filter((t) => t.stop === 'exited');
      const tally = {};
      for (const t of out) {
        const p = t.points.at(-1);
        const e = bl.endNearest(p.x, p.y ?? 0, p.z);
        const key = e ? bl.elements.indexOf(e.element) : -1;
        tally[key] = (tally[key] ?? 0) + 1;
      }
      return tally;
    };

    const bentIdx = bl.elements.indexOf(bl.childAt(sw, 'bend'));
    const straightIdx = bl.elements.indexOf(bl.childAt(sw, 'straight'));

    sw.setVoltage(0);
    const off = where();
    assert(off[straightIdx] === 9, `unpowered should go straight, got ${JSON.stringify(off)}`);

    sw.setVoltage(V());
    const on = where();
    assert(on[bentIdx] === 9, `matched should bend, got ${JSON.stringify(on)}`);
  });

  it('counts only ions that reach the end being aimed at', () => {
    // Otherwise a search told to maximise transmission would discover that a
    // deflector transmits everything at zero volts - straight out the back -
    // and switch it off. That is not tuning a beamline.
    const { bl, sw } = switched();
    const ions = () => discBeam({ ...SPEC, count: 9, radius: 1 });

    assert(bl.mainEnd, 'the column names the end it is aimed at');
    assert(bl.mainEnd.element === bl.childAt(sw, 'bend'), 'which is down the bend');

    sw.setVoltage(0);
    assert(scoreBeamline(bl, ions).transmitted === 0, 'straight through scores nothing');
    sw.setVoltage(V());
    assert(scoreBeamline(bl, ions).transmitted === 9, 'round the bend scores everything');
  });

  it('takes an orphaned branch with the junction it hung from', () => {
    const { bl, sw } = switched();
    assert(bl.elements.length === 4, 'four to begin with');
    bl.remove(bl.elements.indexOf(sw));
    // The first branch inherits the deflector's place; the second had nothing
    // left to hang from, so it goes too rather than floating free.
    assert(bl.elements.length === 2, `expected 2 left, got ${bl.elements.length}`);
    assert(bl.roots().length === 1, 'and the column is still one tree');
  });

  it('refuses to move an element below itself', () => {
    const { bl, sw } = switched();
    const bent = bl.childAt(sw, 'bend');
    assert(!bl.reparent(sw, bent, 'out'), 'that would cut the tree loose from the source');
    assert(bl.childAt(sw, 'bend') === bent, 'and nothing moved');
  });

  it('gives every exit a place of its own to be dropped on', () => {
    /*
      What a drop has to be scored against.

      Picking the nearest ELEMENT and then its first free port was fine while
      everything had one exit. A deflector has three, pointing three different
      ways, so the nearest element does not say where a dropped item goes: it
      landed on whichever port happened to be free first, which from the
      outside looks like the drop going somewhere at random. Every exit has a
      position - the frame a child hanging there would start at - and those
      are far enough apart to tell one from another.
    */
    const { bl, sw } = switched();
    const places = exitsOf(sw).map((x) => ({
      port: x.port,
      o: compose(sw.nominalFrame, x.transform).o,
    }));
    assert(places.length === 3, 'three exits to aim at');
    for (let i = 0; i < places.length; i++) {
      for (let j = i + 1; j < places.length; j++) {
        const d = Math.hypot(
          places[i].o[0] - places[j].o[0],
          places[i].o[1] - places[j].o[1],
          places[i].o[2] - places[j].o[2]
        );
        assert(
          d > mmToM(20),
          `${places[i].port} and ${places[j].port} are only ${mToMm(d).toFixed(1)} mm apart`
        );
      }
    }
    void bl;
  });

  it('keeps the element list in the order the beam visits it', () => {
    const { bl } = switched();
    for (const e of bl.elements) {
      const p = e.from?.parent;
      if (!p) continue;
      assert(
        bl.elements.indexOf(p) < bl.elements.indexOf(e),
        `${e.label} is listed before the element it hangs from`
      );
    }
  });
});

/* ------------------------------------------------------------------ */
/* fringe fields                                                       */
/* ------------------------------------------------------------------ */

describe('Column solves and fringe fields', () => {
  const PLATE = { voltage: -500, bore: 4, margin: 10, housingRadius: 14 };

  /** A charged plate, then a tube of the given bore. */
  function plateThenTube(bore, fringe, extra = []) {
    const bl = new Beamline([
      createElement('drift', { length: 20, bore: 12 }),
      createElement('aperture', PLATE),
      createElement('drift', { length: 60, bore }),
      ...extra,
    ]);
    bl.setFringe(fringe);
    return bl;
  }

  /** On-axis potential `d` mm downstream of the plate's middle. */
  function past(bl, d) {
    const p = bl.elements[1];
    const z = p.zStart + mmToM(p.params.margin + p.params.thickness / 2);
    return bl.potentialAt3D(0, 0, z + mmToM(d));
  }

  it('groups only the elements that can share an r-z grid', () => {
    // Through a Beamline, because a run is a chain along a BRANCH: the
    // grouping is read off the topology, not off positions in an array.
    const bl = new Beamline([
      createElement('drift', { length: 20, bore: 6 }),
      createElement('einzel', { gridStep: 0.8 }),
      createElement('drift', { length: 20, bore: 6 }),
      createElement('bender', { voltage: 40 }),
      createElement('drift', { length: 20, bore: 6 }),
      createElement('aperture', PLATE),
    ]);
    const runs = axisymmetricRuns(bl.elements);
    assert(
      JSON.stringify(runs) === JSON.stringify([[0, 1, 2], [4, 5]]),
      `runs should split at the deflector, got ${JSON.stringify(runs)}`
    );
    // A deflector is in a grounded box of its own, so ending a run there is
    // not an approximation - its field really does stop.
    assert(!canShareGrid(bl.elements[3]), 'a deflector cannot share an r-z grid');
    assert(canShareGrid(bl.elements[1]), 'a lens can');
  });

  it('does not join two branches leaving the same junction', () => {
    // They sit next to each other in the array and point in different
    // directions. Painting both onto one r-z grid would put one branch's
    // hardware on top of the other's.
    const bl = new Beamline([
      createElement('drift', { length: 12, bore: 6 }),
      createElement('bender', { voltage: 40 }),
    ]);
    const sw = bl.elements[1];
    const a = bl.add(createElement('aperture', PLATE), 2, { parent: sw, port: 'bend' });
    const b = bl.add(createElement('aperture', PLATE), 3, { parent: sw, port: 'straight' });

    const runs = axisymmetricRuns(bl.elements);
    const ia = bl.elements.indexOf(a);
    const ib = bl.elements.indexOf(b);
    assert(
      !runs.some((r) => r.includes(ia) && r.includes(ib)),
      `the two branches must not share a grid, got ${JSON.stringify(runs)}`
    );
  });

  it('excludes a misaligned element, which is no longer a body of revolution', () => {
    const e = createElement('drift', { length: 20, bore: 6 });
    assert(canShareGrid(e), 'aligned');
    e.align = { dx: mmToM(1), dy: 0, tiltX: 0, tiltY: 0 };
    assert(!canShareGrid(e), 'and not once it is nudged off the axis');
  });

  it('stops the field dead at the element boundary when off', () => {
    // The thing being fixed. An isolated element's grounded end faces are a
    // numerical device, but they behave exactly like a grounded plate.
    const off = plateThenTube(12, false);
    assert(Math.abs(past(off, 15)) < 1e-9, `field should be gone at 15 mm, got ${past(off, 15)}`);
    assert(Math.abs(past(off, 8)) > 1, 'while still inside the element it is not');
  });

  it('lets the field reach into the next element when on', () => {
    const on = plateThenTube(12, true);
    assert(Math.abs(past(on, 15)) > 10, `expected a real fringe, got ${past(on, 15)} V`);
    // And it falls off, rather than simply not being clipped.
    assert(Math.abs(past(on, 35)) < Math.abs(past(on, 15)) / 10, 'and it decays');
  });

  it('decays inside a grounded pipe at the rate the pipe sets', () => {
    // exp(-j01 z / R) and nothing else: not the bore of the element that made
    // the field, not its voltage. This is the whole basis of shielding by
    // geometry, so it is checked against the closed form rather than against
    // a previous run.
    for (const bore of [12, 8]) {
      const bl = plateThenTube(bore, true);
      const a = past(bl, 18);
      const b = past(bl, 28);
      const fitted = mmToM(10) / Math.log(Math.abs(a / b));
      assertRelClose(
        fitted,
        decayLength(mmToM(bore)),
        0.12,
        `a ${bore} mm pipe should give a ${mToMm(decayLength(mmToM(bore))).toFixed(2)} mm decay`
      );
    }
  });

  it('shields with a grounded plate, which superposition could never do', () => {
    // The reason this needs one grid rather than a sum. A grounded electrode
    // contributes nothing to a superposition - it is at zero everywhere - yet
    // it changes the field completely, because it changes the boundary of the
    // problem rather than adding to its solution.
    const open = plateThenTube(12, true);
    const shielded = plateThenTube(12, true, []);
    // Rebuild with a grounded plate close behind the live one.
    const near = new Beamline([
      createElement('drift', { length: 20, bore: 12 }),
      createElement('aperture', PLATE),
      createElement('drift', { length: 4, bore: 12 }),
      createElement('aperture', { ...PLATE, voltage: 0, margin: 3 }),
      createElement('drift', { length: 60, bore: 12 }),
    ]);
    near.setFringe(true);

    const openAt = Math.abs(past(open, 25));
    const shieldedAt = Math.abs(past(near, 25));
    assert(openAt > 5, `an unshielded plate should still be felt at 25 mm, got ${openAt} V`);
    assert(
      shieldedAt < openAt / 5,
      `a grounded plate should cut it hard: ${openAt.toFixed(2)} -> ${shieldedAt.toFixed(2)} V`
    );
    void shielded;
  });

  it('barely changes an einzel, which is its own Faraday cage', () => {
    // Its outer cylinders are grounded, so they shield the centre electrode
    // almost completely and the isolated solve was already right. Worth
    // asserting: it says the two models agree where they should, which is the
    // only reason to trust them where they differ.
    const build = (fringe) => {
      const bl = new Beamline([
        createElement('drift', { length: 20, bore: 12 }),
        createElement('einzel', {
          gridStep: 0.5, voltage: -300, boreRadius: 5,
          housingRadius: 14, entryDrift: 15, exitDrift: 15,
        }),
        createElement('drift', { length: 40, bore: 12 }),
      ]);
      bl.setFringe(fringe);
      return bl;
    };
    const a = build(false);
    const b = build(true);
    const e = a.elements[1];
    const z = e.zStart + e.length / 2;
    assertRelClose(
      b.potentialAt3D(0, 0, z),
      a.potentialAt3D(0, 0, z),
      0.01,
      'a lens is unchanged by solving it with its neighbours'
    );
  });

  it('keeps voltages a fast adjust over the shared grid', () => {
    // A column solve holds its own copy of the voltages, so changing one has
    // to reach it - but must not re-solve, or the tuner would be unusable.
    const bl = plateThenTube(12, true);
    const before = past(bl, 15);
    const grid = bl.runs[0].grid;

    bl.elements[1].setVoltage(-1000);
    bl.syncRuns();
    const after = past(bl, 15);

    assert(bl.runs[0].grid === grid, 'the same grid, not a fresh solve');
    assertRelClose(after, before * 2, 1e-9, 'and the field scales with the voltage');
  });

  it('turns off cleanly', () => {
    const bl = plateThenTube(12, true);
    assert(bl.runs.length === 1, 'one run while on');
    bl.setFringe(false);
    assert(bl.runs.length === 0, 'none while off');
    assert(Math.abs(past(bl, 15)) < 1e-9, 'and the field is contained again');
  });
});

/* ------------------------------------------------------------------ */
/* voltage tuning                                                      */
/* ------------------------------------------------------------------ */

describe('Voltage optimiser', () => {
  const DEFLECTOR = {
    apertureRadius: 19,
    electrodeThickness: 0.5,
    boxClearance: 0.5,
    // Fully specified, as above: these tests are about the tuner, and should
    // not move when a shipped default does.
    channelWidth: 10.5,
    cornerSize: 0,
    gridStep: 0.5,
  };
  const SPEC = { mass: 100, charge: 1, energy: 1000 };
  const beam = (count = 5) => () => discBeam({ ...SPEC, count, radius: 1.5 });

  const bentColumn = (voltage = 0) =>
    new Beamline([
      createElement('drift', { length: 8, bore: 6 }),
      createElement('bender', { ...DEFLECTOR, voltage }),
      createElement('drift', { length: 25, bore: 6 }),
    ]);

  it('offers only parameters that do not need a re-solve', () => {
    // The optimiser's whole premise: a voltage is a multiplier on a stored
    // solution, so hundreds of trials cost no solver time. If a tunable key
    // were ever marked `rebuild`, `tunableKnobs` must refuse rather than
    // silently re-solve a few hundred times.
    for (const [type, keys] of Object.entries(TUNABLE)) {
      for (const key of keys) {
        const field = ELEMENT_TYPES[type].fields.find((f) => f.key === key);
        assert(field, `${type}.${key} must exist in the registry`);
        assert(!needsRebuild(type, key), `${type}.${key} must not need a re-solve`);
      }
    }
  });

  it('finds every tunable voltage in a column, and nothing else', () => {
    const bl = new Beamline([
      createElement('drift', { length: 8, bore: 6 }),
      createElement('einzel', { boreRadius: 6, centreLength: 20, gap: 6 }),
      createElement('bender', DEFLECTOR),
    ]);
    const knobs = tunableKnobs(bl);
    assert(knobs.length === 2, `expected 2 knobs, got ${knobs.length}`);
    assert(
      knobs.every((k) => k.index > 0),
      'a drift has no voltage to tune'
    );
    // Limits come from the registry, so the optimiser can never propose a
    // setting the user could not have typed.
    for (const k of knobs) {
      const field = ELEMENT_TYPES[k.type].fields.find((f) => f.key === k.key);
      assert(k.min === field.min && k.max === field.max, `${k.key} limits from the registry`);
      assert(k.lo === field.min && k.hi === field.max, `${k.key} sweeps it all without an ion`);
    }
  });

  it('leaves the RF quadrupole alone', () => {
    // Not an oversight. A mass filter's voltages do change transmission, and
    // the setting that transmits most is the one that filters nothing - so a
    // search told to maximise transmission simply turns the RF off. That is
    // the optimiser working correctly on the wrong objective, and the fix is
    // to keep the filter out of its hands.
    const bl = new Beamline([
      createElement('quadrupole', { gridStep: 0.5, length: 60, rfAmplitude: 250, frequency: 2 }),
    ]);
    assert(tunableKnobs(bl).length === 0, 'a mass filter offers no transmission knobs');
    assert(!('quadrupole' in TUNABLE), 'and is not listed as tunable at all');
  });

  it('narrows the sweep to where an element says the answer is', () => {
    // The deflector's range has to reach tens of kilovolts, because real ones
    // run there. A 1 keV beam is matched near 1.7 kV. A coarse sweep of the
    // full range would step straight over that, so when the ion is known the
    // sweep is centred on the element's own closed-form estimate instead.
    const bl = bentColumn(0);
    const wide = tunableKnobs(bl)[0];
    const aimed = tunableKnobs(bl, SPEC)[0];

    assert(aimed.seed !== null, 'the deflector can suggest a voltage');
    assertRelClose(
      aimed.seed,
      matchedVoltage(DEFLECTOR, SPEC.energy, 1),
      1e-12,
      'and the suggestion is the matched voltage'
    );
    assert(aimed.hi - aimed.lo < (wide.hi - wide.lo) / 2, 'the aimed sweep is narrower');
    assert(aimed.min === wide.min && aimed.max === wide.max, 'but the hard limits are unchanged');
    // Symmetric about zero, so a negative ion's opposite polarity is still
    // reachable without anyone having to say so.
    assertClose(aimed.lo, -aimed.hi, 1e-12, 'the sweep spans both polarities');
  });

  it('always tries the suggested voltage, however coarse the sweep', async () => {
    // Five samples across the aimed range is a spacing far wider than the
    // deflector's transmitting window, so an evenly spaced scan alone cannot
    // find it. It works only because the element's own estimate is tried
    // explicitly - which is the point: whether a sweep happens to land on the
    // right value must not depend on the sample count.
    const bl = bentColumn(0);
    const knobs = tunableKnobs(bl, SPEC);
    const result = await optimizeVoltages(bl, beam(), knobs, {
      passes: 1,
      coarse: 5,
      levels: 0,
    });
    assert(
      result.transmitted === result.count,
      `a five-sample sweep should still find it, got ${result.transmitted}/${result.count}`
    );
  });

  it('scores transmission, and rewards partial progress below it', () => {
    const V = matchedVoltage(DEFLECTOR, SPEC.energy, 1);
    const good = scoreBeamline(bentColumn(V), beam());
    const dead = scoreBeamline(bentColumn(0), beam());

    assert(good.transmitted === good.count, `matched should transmit all, got ${good.transmitted}`);
    // An unpowered deflector still passes the beam - straight through, out of
    // the far side of its box. It scores nothing because it delivered nothing
    // to the end being aimed at, which is what "transmitted" has to mean once
    // a column can branch.
    assert(dead.transmitted === 0, 'an unpowered deflector delivers nothing to the bend');
    assert(good.score > dead.score, 'transmitting must outscore not transmitting');

    // The partial-credit and beam-size terms together must never be able to
    // outweigh one transmitted ion, or the search would trade beam away for
    // tidiness.
    assert(
      dead.score < 1 / dead.count,
      `a fully lost beam scored ${dead.score.toFixed(4)}, which is more than one ion is worth`
    );
  });

  it('prefers the tighter beam among settings that all transmit', () => {
    // The tie-break, measured transverse to the EXIT axis. Measured from the
    // origin instead it would be dominated by the bend offset and would say
    // nothing about the beam at all.
    const V = matchedVoltage(DEFLECTOR, SPEC.energy, 1);
    const a = scoreBeamline(bentColumn(V), beam());
    const b = scoreBeamline(bentColumn(V * 1.05), beam());
    assert(a.transmitted === b.transmitted, 'both settings transmit everything');
    assert(a.exitRadius !== null && b.exitRadius !== null, 'both report an exit radius');
    assert(
      a.exitRadius < mmToM(10) && b.exitRadius < mmToM(10),
      'the radius is transverse to the exit axis, not the distance from the origin'
    );
    // Whichever is tighter must be the one that scores higher.
    const tighter = a.exitRadius < b.exitRadius ? a : b;
    const looser = tighter === a ? b : a;
    assert(tighter.score > looser.score, 'the tighter beam wins the tie-break');
  });

  it('recovers a working deflector voltage from nothing', async () => {
    const bl = bentColumn(0);
    const knobs = tunableKnobs(bl, SPEC);
    const before = scoreBeamline(bl, beam());
    assert(before.transmitted === 0, 'starts with nothing reaching the bend');

    const result = await optimizeVoltages(bl, beam(), knobs, { passes: 1, coarse: 11, levels: 2 });

    assert(result.transmitted === result.count, `expected full transmission, got ${result.transmitted}`);
    assert(result.improved, 'and it should report that it moved something');

    // The voltage it lands on must be the physical one, not an artefact of
    // the search: within the measured window around the derived value.
    const V = matchedVoltage(DEFLECTOR, SPEC.energy, 1);
    const found = bl.elements[1].params.voltage;
    assertRelClose(found, V, 0.15, 'the tuner finds the matched voltage');
  });

  it('leaves a beamline exactly as it found it when nothing helps', async () => {
    // A drift-only column has no knobs at all; the optimiser must not invent
    // any, and must not disturb what is there.
    const bl = new Beamline([createElement('drift', { length: 20, bore: 6 })]);
    const knobs = tunableKnobs(bl);
    assert(knobs.length === 0, 'a drift column has nothing to tune');
    const result = await optimizeVoltages(bl, beam(), knobs, { passes: 1 });
    assert(!result.improved, 'nothing to improve');
    assert(result.transmitted === result.count, 'and the beam still gets through');
  });

  it('stops when asked and keeps the best setting found so far', async () => {
    const bl = bentColumn(0);
    const knobs = tunableKnobs(bl, SPEC);
    let seen = 0;
    const result = await optimizeVoltages(bl, beam(), knobs, {
      shouldStop: () => ++seen > 20,
    });
    assert(result.cancelled, 'should report that it was cancelled');
    assert(result.evaluations < 40, `should stop early, ran ${result.evaluations} trials`);
    // Whatever it kept must actually be what is set on the element.
    assert(
      bl.elements[1].params.voltage === result.values[0],
      'the element carries the value the optimiser reports'
    );
  });

  it('reaches the same answer whether the scan flies coarse or fine', async () => {
    // The search ranks settings at a coarser time step than the answer is
    // judged at, which is most of why it is fast. That is only legitimate if
    // the ranking survives it, so this compares the two directly rather than
    // trusting that it does.
    const fine = bentColumn(0);
    const coarse = bentColumn(0);
    const a = await optimizeVoltages(fine, beam(), tunableKnobs(fine, SPEC), {
      scanSpeed: 1,
      polish: false,
    });
    const b = await optimizeVoltages(coarse, beam(), tunableKnobs(coarse, SPEC), {
      scanSpeed: 4,
      polish: false,
    });
    assert(
      a.transmitted === b.transmitted,
      `same transmission: ${a.transmitted} vs ${b.transmitted}`
    );
    assertRelClose(b.values[0], a.values[0], 0.1, 'and lands on the same voltage');
  });

  it('does not fly the same settings twice', async () => {
    // Coordinate descent revisits points constantly - every narrowing level
    // re-samples its own bracket centre - and the flight is deterministic, so
    // re-running it is pure waste.
    const bl = bentColumn(0);
    const knobs = tunableKnobs(bl, SPEC);
    const withMemo = await optimizeVoltages(bl, beam(), knobs, { polish: false });
    // The budget is an upper bound that assumes nothing repeats; coming in
    // meaningfully under it is the memo doing its job.
    const perKnob = 13 + 3 + 7 * 3;
    const ceiling = 1 + 2 * knobs.length * perKnob;
    assert(
      withMemo.evaluations < ceiling * 0.9,
      `expected repeats to be skipped: ${withMemo.evaluations} of a possible ${ceiling}`
    );
  });

  it('writes a knob through to the field, not just to the parameter', async () => {
    // `applyKnob` has to call the element's setter: writing `params.voltage`
    // alone would change the readout and leave the field untouched, so the
    // optimiser would score every trial against the same trajectory.
    const bl = bentColumn(0);
    const [knob] = tunableKnobs(bl, SPEC);
    const before = bl.elements[1].potentialAt(mmToM(6), 0, mmToM(26));
    applyKnob(bl, knob, 500);
    const after = bl.elements[1].potentialAt(mmToM(6), 0, mmToM(26));
    assertClose(before, 0, 1e-12, 'no voltage, no potential');
    assert(Math.abs(after) > 1, `the field should follow the knob, got ${after}`);
  });
});

/* ------------------------------------------------------------------ */
/* reduced-Hessian refinement                                          */
/* ------------------------------------------------------------------ */

describe('Symmetric eigensolver', () => {
  /** Assert A v = lambda v and that the vectors are orthonormal. */
  function check(A, label) {
    const { values, vectors } = symmetricEigen(A);
    const n = A.length;
    assert(values.length === n, `${label}: one eigenvalue per dimension`);

    for (let k = 0; k < n; k++) {
      const v = vectors[k];
      for (let i = 0; i < n; i++) {
        const Av = A[i].reduce((s, a, j) => s + a * v[j], 0);
        assertClose(Av, values[k] * v[i], 1e-10, `${label}: A v = lambda v`);
      }
    }
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const dot = vectors[i].reduce((s, c, k) => s + c * vectors[j][k], 0);
        assertClose(dot, i === j ? 1 : 0, 1e-10, `${label}: orthonormal`);
      }
    }
    for (let k = 1; k < n; k++) {
      assert(values[k - 1] >= values[k], `${label}: sorted, stiffest first`);
    }
    return values;
  }

  it('diagonalises a matrix with known eigenvalues', () => {
    const v = check([[2, 1], [1, 2]], '[[2,1],[1,2]]');
    assertClose(v[0], 3, 1e-12, 'larger eigenvalue');
    assertClose(v[1], 1, 1e-12, 'smaller eigenvalue');
  });

  it('handles a repeated eigenvalue', () => {
    // The case that matters: two equally sloppy directions is exactly what a
    // degenerate beamline produces, and a solver that cannot give orthogonal
    // vectors for it would report a null space that is not one.
    const v = check([[4, 0], [0, 4]], 'degenerate');
    assertClose(v[0], 4, 1e-12, 'both equal');
    assertClose(v[1], 4, 1e-12, 'both equal');
  });

  it('separates a stiff direction from a flat one', () => {
    const v = check([[1, 0], [0, 1e-9]], 'nearly singular');
    assertClose(v[0], 1, 1e-12, 'the stiff one');
    assert(Math.abs(v[1]) < 1e-8, `the flat one, got ${v[1]}`);
  });

  it('handles an indefinite matrix', () => {
    // A saddle, which is where a beamline that is not yet tuned normally sits.
    const v = check([[0, 1], [1, 0]], 'saddle');
    assertClose(v[0], 1, 1e-12, 'one uphill');
    assertClose(v[1], -1, 1e-12, 'one downhill');
  });

  it('diagonalises a larger matrix', () => {
    check([[6, 2, 1], [2, 5, 1], [1, 1, 4]], '3x3');
    check([[2, -1, 0, 0], [-1, 2, -1, 0], [0, -1, 2, -1], [0, 0, -1, 2]], '4x4 tridiagonal');
  });
});

describe('Refining a tuned beamline', () => {
  const SPEC = { mass: 100, charge: 1, energy: 50 };
  const ions = () => discBeam({ ...SPEC, count: 9, radius: 1.0 });

  const columnAt = (voltage) =>
    new Beamline([
      createElement('drift', { length: 12, bore: 6 }),
      createElement('einzel', {
        gridStep: 0.8, voltage: -300, boreRadius: 6,
        housingRadius: 15, entryDrift: 14, exitDrift: 14,
      }),
      createElement('drift', { length: 14, bore: 6 }),
      createElement('bender', { voltage }),
      createElement('drift', { length: 40, bore: 6 }),
    ]);

  /*
    The deflector voltage this fixture runs at, FOUND rather than assumed.

    These tests need a column that delivers every ion, because `beamQuality`
    reports no merit at all unless the beam arrives intact - that is the next
    test down. The obvious way to get one is to set the deflector to its
    closed-form V0, and that is what this did.

    But V0 is an estimate of an ideal quadrupole, deliberately: the element
    ships with a tuner precisely because the solved field wants a few per cent
    either way, and §9.3 of the physics notes refuses to quote a calibration
    table for exactly this reason. So a fixture that pins V0 is asserting
    something the rest of the codebase says is not true, and it breaks whenever
    the electrode geometry moves - which it just did, taking this test from
    9 of 9 to 8 of 9 on a change that was nothing to do with what it tests.

    Searching costs one extra Laplace solve for the whole block, because a
    deflector voltage is only a multiplier on the solved field: the scan reuses
    one column and calls setVoltage. V0 is tried first, so a well-calibrated
    element pays a single flight for this.
  */
  let found = null;
  function deliveringVoltage() {
    if (found !== null) return found;
    const V0 = startingParams('bender', SPEC).voltage;
    const probe = columnAt(V0);
    const bender = probe.elements[3];
    let best = { arrived: -1, v: V0 };
    for (const f of [1, 0.95, 1.05, 0.9, 1.1, 0.85, 1.15]) {
      bender.setVoltage(f * V0);
      const arrived = beamQuality(probe, ions).arrived;
      if (arrived > best.arrived) best = { arrived, v: f * V0 };
      if (arrived === 9) break;
    }
    found = best.v;
    return found;
  }

  const column = () => columnAt(deliveringVoltage());

  it('measures the beam at the target plane, not where the ion stopped', () => {
    // The difference between a merit that can be differentiated and one that
    // cannot. The last recorded point sits wherever the adaptive step left
    // off, which moves as the voltages change - so a finite difference reads
    // that jitter rather than the physics.
    const bl = column();
    const target = bl.mainEnd;
    const frame = bl.endFrame(target);
    const q = beamQuality(bl, ions);
    assert(
      q.arrived === 9,
      `all nine should arrive, got ${q.arrived} at the best of seven voltages ` +
        `(${deliveringVoltage().toFixed(1)} V) - if no voltage delivers the beam, the ` +
        'loss is geometric rather than a calibration, so look at the deflector ' +
        'channel width rather than at the tuning'
    );

    // Every contributing ion is measured exactly on the plane, so the distance
    // along the exit axis is zero rather than a step's worth of overshoot.
    const { tracks } = flyBeam(bl, ions(), { cfl: 0.05, maxSteps: 200000 });
    const f = forwardOf(frame);
    let worst = 0;
    for (const t of tracks) {
      const p = t.points[t.points.length - 1];
      const along =
        (p.x - frame.o[0]) * f[0] + ((p.y ?? 0) - frame.o[1]) * f[1] + (p.z - frame.o[2]) * f[2];
      worst = Math.max(worst, along);
    }
    assert(worst > 0, 'the recorded endpoint really does overshoot the plane');
  });

  it('reports no merit when the beam is not fully delivered', () => {
    // A merit that counted a lost ion as a tighter beam would reward throwing
    // the edge of the beam away.
    const bl = column();
    bl.elements[3].setVoltage(0); // deflector off: nothing reaches the bend
    const q = beamQuality(bl, ions);
    assert(q.merit === null, `expected no merit, got ${q.merit}`);
  });

  it('tightens the beam without losing any of it', async () => {
    const bl = column();
    const knobs = tunableKnobs(bl, SPEC);
    const before = beamQuality(bl, ions);
    const r = await refineNullSpace(bl, ions, knobs, {});

    assert(r.rmsAfter <= r.rmsBefore, 'never worse than it started');
    const after = beamQuality(bl, ions);
    assert(
      after.arrived === before.arrived,
      `transmission must hold: ${before.arrived} -> ${after.arrived}`
    );
    if (r.improved) assert(r.steps > 0, 'an improvement means it actually stepped');
  });

  it('puts the voltages back if it cannot improve on them', async () => {
    // A refinement that leaves a column worse than it found it is worse than
    // no refinement.
    const bl = column();
    const knobs = tunableKnobs(bl, SPEC);
    const before = readKnobs(bl, knobs);
    const r = await refineNullSpace(bl, ions, knobs, { iterations: 1 });
    if (!r.improved) {
      const now = readKnobs(bl, knobs);
      assert(
        now.every((v, i) => v === before[i]),
        'unchanged when nothing was gained'
      );
    }
  });

  it('gives every knob a natural size from the element itself', async () => {
    // Without one, a knob's size falls back to a sixth of its full range: for
    // a lens allowed +/-40 kV that is 13 kV, so the finite-difference probe is
    // in kilovolts around an answer of a few hundred volts and the Hessian is
    // meaningless. Every tunable field declares a scale for exactly this.
    const bl = column();
    for (const k of tunableKnobs(bl, SPEC)) {
      assert(k.seed != null, `${k.label} has no natural size`);
      assert(k.seed > 0 && Number.isFinite(k.seed), `${k.label} seed ${k.seed}`);
      assert(
        k.seed < Math.abs(k.max - k.min) / 4,
        `${k.label}: a seed of ${k.seed} is a large fraction of its whole range`
      );
    }
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

  it('does not let one branch of a folded column claim another', () => {
    // Elements answer `contains` on their axial extent alone, so that an ion
    // inside an element's length but outside its bore is still that element's
    // business and is reported as hitting its wall. In a straight column that
    // is unambiguous. Fold the column through two right angles and the last
    // drift runs back alongside the first, inside its axial range but eighty
    // millimetres off its axis - and every ion entering the last drift was
    // being reported as striking the wall of the first.
    const V = matchedVoltage({ apertureRadius: 19 }, 50, 1);
    const bl = new Beamline([
      createElement('drift', { length: 12, bore: 5 }),
      createElement('bender', { voltage: V }),
      createElement('drift', { length: 30, bore: 5 }),
      createElement('bender', { voltage: V }),
      createElement('drift', { length: 30, bore: 5 }),
    ]);
    assertClose(forwardOf(bl.exitFrame)[2], -1, 1e-9, 'the column folds right back');

    // A point on the axis of the LAST drift, well inside it.
    const last = bl.elements[4];
    const g = toGlobal(last.frame, [0, 0, mmToM(15)]);
    const hit = bl.locate(g);
    assert(hit, 'a point on the last drift axis belongs to some element');
    assert(hit.index === 4, `expected the last drift, got element ${hit.index}`);
    assert(!bl.strikes(g[0], g[1], g[2]), 'and it is free space, not a wall');

    // And the beam actually survives the second bend.
    const { tracks } = flyBeam(bl, discBeam({ mass: 100, charge: 1, energy: 50, count: 9, radius: 1 }), {
      cfl: 0.05,
      maxSteps: 600000,
    });
    const n = tracks.filter((t) => t.stop === 'exited').length;
    assert(n >= 5, `a double bend should transmit most of the beam, got ${n} of 9`);
  });

  it('gives the same answer with the lookup hint as without it', () => {
    // Element lookup remembers where the last one landed, because a flight of
    // nine ions asks six hundred thousand times and almost every question is
    // about a point a fraction of a millimetre from the previous one. It is a
    // hint and nothing more: every use re-checks the point against that
    // element. This pins that it changes no answer at all.
    const V = matchedVoltage({ apertureRadius: 19 }, 50, 1);
    const make = () =>
      new Beamline([
        createElement('drift', { length: 12, bore: 5 }),
        createElement('bender', { voltage: V }),
        createElement('drift', { length: 30, bore: 5 }),
        createElement('bender', { voltage: V }),
        createElement('drift', { length: 30, bore: 5 }),
      ]);
    const ions = () => discBeam({ mass: 100, charge: 1, energy: 50, count: 9, radius: 1 });

    const hinted = flyBeam(make(), ions(), { cfl: 0.05, maxSteps: 600000 }).tracks;

    // The same column, with the hint cleared before every single lookup.
    const bl = make();
    const base = Object.getPrototypeOf(bl).locate;
    bl.locate = function (g) {
      this.layout();
      return base.call(this, g);
    };
    const plain = flyBeam(bl, ions(), { cfl: 0.05, maxSteps: 600000 }).tracks;

    assert(hinted.length === plain.length, 'same number of tracks');
    for (let i = 0; i < hinted.length; i++) {
      assert(hinted[i].stop === plain[i].stop, `ion ${i}: same fate`);
      assert(hinted[i].points.length === plain[i].points.length, `ion ${i}: same step count`);
      const a = hinted[i].points[hinted[i].points.length - 1];
      const b = plain[i].points[plain[i].points.length - 1];
      for (const k of ['x', 'y', 'z', 'vx', 'vy', 'vz']) {
        assert(a[k] === b[k], `ion ${i}: ${k} identical to the last bit`);
      }
    }
  });

  it('still blames the right element for a genuine wall strike', () => {
    // The transverse bound on that claim must not go so far as to stop an ion
    // that really does hit a wall being attributed to it.
    const bl = new Beamline([
      createElement('drift', { length: 20, bore: 5 }),
      createElement('drift', { length: 20, bore: 5 }),
    ]);
    const g = [mmToM(6), 0, mmToM(10)]; // inside the first, outside its bore
    const hit = bl.locate(g);
    assert(hit?.index === 0, 'the wall it hit is the first drift');
    assert(bl.strikes(g[0], g[1], g[2]), 'and it counts as a strike');
    // Far outside every element, though, belongs to nothing.
    assert(bl.locate([mmToM(500), 0, mmToM(10)]) === null, 'half a metre out is not a strike');
  });

  it('measures the beam transversely to the axis it is on', () => {
    // What the cross-section profile depends on. After a right-hand bend the
    // beam runs along -x, so global x is the direction of travel and global z
    // is transverse; an ion on the axis is far from the origin in x but zero
    // from its own axis. Reading the global x and y instead would draw a
    // profile that smears with distance flown.
    const V = matchedVoltage({ apertureRadius: 19 }, 50, 1);
    const bl = new Beamline([
      createElement('drift', { length: 12, bore: 5 }),
      createElement('bender', { voltage: V }),
      createElement('drift', { length: 30, bore: 5 }),
    ]);
    const after = bl.elements[2];
    const g = toGlobal(after.frame, [mmToM(1), mmToM(0.5), mmToM(20)]);

    const local = toLocal(bl.locate(g).element.frame, g);
    assertClose(mToMm(Math.hypot(local[0], local[1])), Math.hypot(1, 0.5), 1e-9, 'offset from its own axis');
    // Measured from the global axis instead it would be tens of millimetres.
    assert(
      mToMm(Math.hypot(g[0], g[1])) > 20,
      'while the global offset is large, which is the thing not to plot'
    );
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

  it('paints an element identically wherever it sits in the grid', () => {
    // Node coordinates come from `z0 + i * step`, geometry from `mm * 1e-3`.
    // Those two routes to the same number differ in the last bit, and which
    // way depends on the element's absolute position. Without a tolerance an
    // electrode whose edge lands on a node is painted one step shorter or
    // longer purely because of where it was placed, so the same element
    // solves to two different fields in two different beamlines.
    //
    // Measured before the fix: a 15 mm centre electrode came out 14.8 mm in
    // one position and 15.0 mm in another, moving the on-axis potential by
    // 7.8 V out of 300 because that node sits where the field changes at
    // 20 V/mm.
    const COMMON = {
      gridStep: 0.4, voltage: -300, boreRadius: 5, wallThickness: 2,
      outerLength: 15, centreLength: 15, gap: 4, housingRadius: 14,
    };

    const extentOf = (element, name, offsetM) => {
      const g = element.grid;
      const id = g.electrodeNames.indexOf(name);
      let lo = Infinity;
      let hi = -Infinity;
      for (let j = 0; j < g.nr; j++) {
        for (let i = 0; i < g.nz; i++) {
          if (g.electrodeId[g.idx(i, j)] !== id) continue;
          lo = Math.min(lo, g.zAt(i) + offsetM);
          hi = Math.max(hi, g.zAt(i) + offsetM);
        }
      }
      return { lo, hi };
    };

    // The same lens, reached two ways: shifted by a 10 mm drift, or with that
    // 10 mm folded into its own margins. The metal is in the same place.
    const shifted = createElement('einzel', { ...COMMON, entryDrift: 16, exitDrift: 16 });
    const direct = createElement('einzel', { ...COMMON, entryDrift: 26, exitDrift: 26 });

    for (const name of ['entrance', 'centre', 'exit']) {
      const a = extentOf(shifted, name, mmToM(10));
      const b = extentOf(direct, name, 0);
      assertClose(a.lo, b.lo, 1e-9, `${name} start must not depend on placement`);
      assertClose(a.hi, b.hi, 1e-9, `${name} end must not depend on placement`);
    }

    // And therefore the fields must agree, including in the steep region at
    // the downstream edge of the centre electrode where the bug showed up.
    let worst = 0;
    for (let zmm = 40; zmm <= 80; zmm += 0.2) {
      const a = shifted.potentialAt(0, 0, mmToM(zmm - 10));
      const b = direct.potentialAt(0, 0, mmToM(zmm));
      worst = Math.max(worst, Math.abs(a - b));
    }
    assert(
      worst < 0.5,
      `the same lens placed two ways differs by ${worst.toFixed(2)} V on axis`
    );
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
    const quad = createElement('quadrupole', { gridStep: 0.2 });
    const bl = new Beamline([createDrift({ length: 20 }), quad]);
    assertRelClose(bl.lengthScale, mmToM(0.2), 1e-12, 'finest element wins');
    assertRelClose(
      bl.shortestPeriod,
      1 / (quad.params.frequency * 1e6),
      1e-12,
      'RF period is exposed to the integrator'
    );
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
