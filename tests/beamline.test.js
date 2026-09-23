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
import { createBender, matchedVoltage } from '../src/elements/bender.js';
import {
  compose,
  translation,
  yawFrame,
  toGlobal,
  toLocal,
  vectorToGlobal,
  forwardOf,
} from '../src/frames.js';
import { makeIon, discBeam } from '../src/ion.js';
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

describe('Bender', () => {
  const PARAMS = { bendRadius: 40, bendAngle: 90, gap: 8, height: 16, gridStep: 0.3 };
  const ENERGY = 50;

  it('matches the closed-form plate voltage', () => {
    // For the central orbit, the electric force supplies the centripetal one:
    //   qE = 2T/R,  with  E = V / (R ln(r2/r1))   =>   V = (2T/q) ln(r2/r1)
    // which for a narrow gap approaches the familiar 2Td/(qR).
    const V = matchedVoltage(PARAMS, ENERGY, 1);
    const r1 = mmToM(PARAMS.bendRadius - PARAMS.gap / 2);
    const r2 = mmToM(PARAMS.bendRadius + PARAMS.gap / 2);
    assertRelClose(V, 2 * ENERGY * Math.log(r2 / r1), 1e-12, 'closed form');
    // And the narrow-gap limit, to the accuracy that limit deserves.
    assertRelClose(V, (2 * ENERGY * PARAMS.gap) / PARAMS.bendRadius, 0.01, 'narrow-gap estimate');
  });

  it('scales the matched voltage with energy and charge', () => {
    const base = matchedVoltage(PARAMS, 50, 1);
    assertRelClose(matchedVoltage(PARAMS, 100, 1), 2 * base, 1e-12, 'V ~ T');
    assertRelClose(matchedVoltage(PARAMS, 50, 2), base / 2, 1e-12, 'V ~ 1/q');
  });

  it('turns the column by its bend angle', () => {
    // The whole point of the element: its exit faces somewhere else, so
    // everything downstream turns with it.
    for (const deg of [30, 90, 127]) {
      const bl = new Beamline([
        createElement('drift', { length: 5, bore: 4 }),
        createElement('bender', { ...PARAMS, bendAngle: deg }),
      ]);
      const dir = forwardOf(bl.exitFrame);
      const turned = (Math.atan2(-dir[0], dir[2]) * 180) / Math.PI;
      assertRelClose(turned, deg, 1e-9, `a ${deg} degree bend`);
    }
  });

  it('places the element after a bend on the turned axis', () => {
    const bl = new Beamline([
      createElement('bender', { ...PARAMS, bendAngle: 90 }),
      createElement('drift', { length: 20, bore: 4 }),
    ]);
    const after = bl.elements[1].frame;
    // A 90 degree bend of radius R puts the exit at (-R, 0, R) and pointing
    // along -x, so the drift after it runs in -x from there.
    const R = mmToM(40);
    assertClose(after.o[0], -R, 1e-12, 'exit sits one radius to the side');
    assertClose(after.o[2], R, 1e-12, 'and one radius downstream');
    const dir = forwardOf(after);
    assertClose(dir[0], -1, 1e-12, 'and the drift runs across the original axis');
  });

  it('carries a matched ion round the arc and out', () => {
    const V = matchedVoltage(PARAMS, ENERGY, 1);
    const bl = new Beamline([
      createElement('drift', { length: 10, bore: 4 }),
      createElement('bender', { ...PARAMS, voltage: V }),
      createElement('drift', { length: 20, bore: 4 }),
    ]);
    const { points, stop } = flyIon(
      bl, makeIon({ mass: 100, charge: 1, energy: ENERGY }), { cfl: 0.05 }
    );
    assert(stop === 'exited', `a matched ion should get through, got ${stop}`);

    const last = points[points.length - 1];
    const turned = (Math.atan2(-last.vx, last.vz) * 180) / Math.PI;
    // Within a degree or so of nominal. The remaining error is real: the
    // grounded lids sit close enough to perturb the field from the ideal
    // cylindrical-capacitor form the matched voltage assumes.
    assertClose(turned, 90, 2, 'a matched ion follows the nominal bend');
  });

  it('puts a badly mismatched ion into a plate', () => {
    // A bender at the wrong voltage does not bend the beam slightly wrongly,
    // it disperses it onto an electrode. That is also what makes it an energy
    // filter rather than merely a corner.
    const V = matchedVoltage(PARAMS, ENERGY, 1);
    const bl = new Beamline([
      createElement('drift', { length: 10, bore: 4 }),
      createElement('bender', { ...PARAMS, voltage: V }),
      createElement('drift', { length: 20, bore: 4 }),
    ]);
    for (const energy of [ENERGY * 0.8, ENERGY * 1.2]) {
      const { stop } = flyIon(
        bl, makeIon({ mass: 100, charge: 1, energy }), { cfl: 0.05 }
      );
      assert(
        stop === 'electrode',
        `an ion ${((energy / ENERGY - 1) * 100).toFixed(0)} % off energy should be ` +
          `dispersed into a plate, got ${stop}`
      );
    }
  });

  it('keeps an on-orbit ion in the bend plane', () => {
    // A cylindrical sector bends horizontally and does essentially nothing
    // vertically - it has no vertical focusing at all in the ideal geometry,
    // which is exactly why spherical deflectors exist. An ion launched in the
    // plane must stay in it.
    const V = matchedVoltage(PARAMS, ENERGY, 1);
    const bl = new Beamline([
      createElement('drift', { length: 10, bore: 4 }),
      createElement('bender', { ...PARAMS, voltage: V }),
      createElement('drift', { length: 20, bore: 4 }),
    ]);
    const { points } = flyIon(
      bl, makeIon({ mass: 100, charge: 1, energy: ENERGY }), { cfl: 0.05 }
    );
    for (const p of points) {
      assertClose(p.y, 0, mmToM(1e-3), `left the bend plane at z = ${mToMm(p.z).toFixed(1)} mm`);
    }
  });

  it('solves on a cylindrical band that excludes the axis', () => {
    // The bender reuses the axisymmetric stencil with the roles of the axes
    // reinterpreted, but its radial band sits about the bend radius rather
    // than starting at zero - so the singular on-axis stencil must NOT be
    // applied, and row j = 0 is an ordinary wall.
    const b = createBender(PARAMS);
    assert(b.grid.r0 > 0, 'the band should not reach the axis');
    assert(!b.grid.includesAxis, 'and must not claim to contain it');
    // Which means j = 0 has to be part of the enclosure, or the problem is
    // unbounded on that side.
    assert(b.grid.isElectrode(Math.floor(b.grid.nz / 2), 0), 'inner wall must be fixed');
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
