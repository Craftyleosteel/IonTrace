/**
 * Linear quadrupole - four parallel rods, optionally driven at radio
 * frequency.
 *
 * Why this element forced the rest of IonTrace into three dimensions
 * ------------------------------------------------------------------
 * Every other element here is rotationally symmetric, so its potential does
 * not depend on the azimuth and a two-dimensional (z, r) solve is exact. A
 * quadrupole is not. Its ideal potential is
 *
 *     phi(x, y, t) = (U + V cos(Omega t)) (x^2 - y^2) / r0^2
 *
 * which is positive along x and negative along y: the azimuth is the whole
 * point. It cannot be represented on an axisymmetric grid at all, and an ion
 * in it does not stay in a plane containing the axis. Hence 3D ions.
 *
 * Note the sign structure. At any instant the field CONVERGES in one
 * transverse plane and DIVERGES in the other. A DC quadrupole is therefore
 * always unstable in one plane - it is a singlet, and needs a partner of
 * opposite sign to make a net-focusing doublet. What makes an RF quadrupole
 * different is that the sign alternates faster than the ion can escape, so
 * the time-averaged force is inward in BOTH planes for the right parameters.
 * That is the whole basis of the mass filter.
 *
 * Solve
 * -----
 * The rods are uniform along z, so the field is a two-dimensional problem in
 * the TRANSVERSE plane rather than the meridional one. The same planar
 * five-point stencil solves it, on a grid whose two axes are x and y. Real
 * round rods are painted rather than ideal hyperbolae, so the solution carries
 * the higher multipoles a real rod set has - the 12-pole in particular, which
 * is what makes rod radius against field radius an engineering choice rather
 * than a detail.
 *
 * The two rod pairs are always driven antisymmetrically, at +W(t) and -W(t),
 * so the field is exactly linear in W and one unit solution suffices:
 *
 *     E(x, y, t) = W(t) . E_unit(x, y),   W(t) = U + V cos(Omega t + phase)
 *
 * No re-solve is needed as the RF swings, and nothing is approximated by the
 * factorisation - it is the same fast-adjust superposition used everywhere
 * else, with a time-dependent coefficient.
 *
 * Hard edge
 * ---------
 * The field is taken as uniform along the rods and zero outside them. Real
 * fringe fields at the rod ends are NOT modelled, and they matter: in a real
 * mass filter the entrance fringe is a well-known cause of transmission loss,
 * because an ion crosses a region where the RF is still ramping up. Two
 * consequences to be aware of:
 *
 *   - transmission through this element is optimistic;
 *   - an off-axis ion's potential energy jumps at the rod ends, because the
 *     potential is discontinuous there. The transverse dynamics inside the
 *     rods, which is what sets stability, are unaffected.
 *
 * Softening the edge with a longitudinal envelope g(z) would be worse, not
 * better: g(z)(x^2-y^2) does not satisfy Laplace's equation unless g'' = 0, so
 * it would trade a known approximation for a field that is not a field.
 */

import { PotentialArray, PLANAR } from '../grid.js';
import { solveBasis } from '../laplace.js';
import { Field } from '../field.js';
import { mmToM, ELEMENTARY_CHARGE, ATOMIC_MASS_UNIT } from '../constants.js';

export const QUADRUPOLE_DEFAULTS = {
  fieldRadius: 4, // mm, r0 - axis to rod surface
  rodRadius: 4.6, // mm, 1.1487 x r0 is the classic round-rod optimum
  // 150 mm at 2 MHz gives a 100 u, 50 eV ion about thirty RF cycles in the
  // rods. That number matters more than it looks: Mathieu stability is an
  // asymptotic property, and an ion crossing in a handful of cycles can be
  // thrown out whatever its (a, q) says. At the previous 60 mm and 1 MHz it
  // saw six, and nothing was transmitted.
  length: 150, // mm, rod length
  dcVoltage: 0, // V, U
  rfAmplitude: 250, // V, V (zero-to-peak) - see amplitudeForQ
  frequency: 2.0, // MHz
  phase: 0, // degrees
  housingRadius: 12, // mm
  gridStep: 0.15, // mm, in the transverse plane
};

export function createQuadrupole(params = {}, solverOpts = {}) {
  const p = { ...QUADRUPOLE_DEFAULTS, ...params };
  const warnings = [];

  const r0 = mmToM(p.fieldRadius);
  const rod = mmToM(p.rodRadius);
  const housing = mmToM(p.housingRadius);

  if (p.fieldRadius + 2 * p.rodRadius >= p.housingRadius * 2) {
    warnings.push(
      'Rods nearly fill the housing; the enclosure is distorting the field.'
    );
  }
  const ratio = p.rodRadius / p.fieldRadius;
  if (Math.abs(ratio - 1.1487) > 0.15) {
    warnings.push(
      `Rod/field radius ratio is ${ratio.toFixed(3)}. The value that cancels ` +
        'the 12-pole for round rods is 1.1487; away from it the field departs ' +
        'measurably from an ideal quadrupole.'
    );
  }

  // Transverse grid: both axes are spatial coordinates in the x-y plane, so
  // it straddles the origin in each. The PotentialArray's "z" axis carries x
  // and its "r" axis carries y - the planar stencil does not care which
  // Cartesian pair it is solving.
  //
  // The node count is forced ODD and the extent derived from it, rather than
  // the other way round. Taking n = round(2R/h) + 1 and spanning -R..-R+(n-1)h
  // leaves the domain lopsided by up to half a grid step, which destroys the
  // one property this element depends on: the enclosure must share the rods'
  // four-fold symmetry. Without that, E_y is not zero on the y = 0 plane, and
  // an ion launched there is pushed out of it. Measured before this fix, a
  // planar ion drifted 197 um - five per cent of the aperture - purely
  // because the grid was 0.15 mm wider on one side than the other.
  //
  // An odd count also puts a node exactly on the axis, which is where the
  // potential is zero by symmetry and where ions most often sit.
  const half = Math.max(2, Math.round(p.housingRadius / p.gridStep));
  const step = mmToM(p.gridStep);
  const extent = half * step;
  const grid = new PotentialArray({
    nz: 2 * half + 1,
    nr: 2 * half + 1,
    step,
    symmetry: PLANAR,
    z0: -extent,
    r0: -extent,
  });

  // This grid's two axes are x and y; the beam runs perpendicular to both. So
  // the z end faces are the grounded housing on the +x and -x sides, real
  // metal, not the apertures the default assumes.
  grid.openFaces = { zMin: false, zMax: false };

  const encl = grid.addElectrode('housing');
  const poleA = grid.addElectrode('poleA'); // on the x axis
  const poleB = grid.addElectrode('poleB'); // on the y axis
  grid.paintEnclosure(encl);

  // Rod centres sit at r0 + rodRadius from the axis, so the rod SURFACE is
  // tangent to the field radius r0. That is what r0 means.
  const centre = r0 + rod;
  // The tolerance matters here for the same reason as elsewhere: a rod
  // surface passing exactly through a node must be painted the same way
  // whichever of the four rods it belongs to, or the set is not symmetric.
  const disc = (cx, cy) => (x, y) =>
    Math.hypot(x - cx, y - cy) <= rod + grid.step * 1e-6;

  const paintedA =
    grid.paint(poleA, disc(+centre, 0)) + grid.paint(poleA, disc(-centre, 0));
  const paintedB =
    grid.paint(poleB, disc(0, +centre)) + grid.paint(poleB, disc(0, -centre));
  if (paintedA === 0 || paintedB === 0) {
    throw new Error('Quadrupole rods covered no grid nodes; check the geometry');
  }

  const { basis, reports } = solveBasis(grid, solverOpts);
  if (reports.some((r) => !r.converged)) {
    warnings.push('Laplace solve did not converge; this field is not trustworthy.');
  }

  // One unit solution: +1 V on the x rods, -1 V on the y rods. Everything the
  // element ever does is this map scaled by W(t).
  const unit = new Field(grid, basis);
  unit.setVoltages({ housing: 0, poleA: 1, poleB: -1 });

  const length = mmToM(p.length);
  const omega = () => 2 * Math.PI * p.frequency * 1e6;
  const phaseRad = () => (p.phase * Math.PI) / 180;

  /** Instantaneous pole voltage W(t) = U + V cos(Omega t + phase). */
  const drive = (t) =>
    p.dcVoltage + p.rfAmplitude * Math.cos(omega() * t + phaseRad());

  return {
    type: 'quadrupole',
    label: 'Quadrupole',
    params: p,
    length,
    bore: r0,
    // The rod surfaces are tangent to the field radius, so inside it the
    // aperture is clear. See `isFree` in beamline.js.
    clearBore: r0,
    outerRadius: extent,
    lengthScale: grid.step,
    warnings,
    grid,
    field: unit,

    /**
     * Shortest timescale the field varies on. The integrator uses this to
     * keep its step short enough to resolve the RF rather than alias it -
     * without it, a large enough step would sample the same phase every time
     * and an unstable ion would look perfectly confined.
     */
    get shortestPeriod() {
      return p.rfAmplitude === 0 || p.frequency === 0
        ? null
        : 1 / (p.frequency * 1e6);
    },

    contains(x, y, zl) {
      return zl >= 0 && zl <= length;
    },

    fieldAt(x, y, zl, t) {
      if (zl < 0 || zl > length) return { Ex: 0, Ey: 0, Ez: 0 };
      // The grid's two axes are x and y, so its "axial" component is Ex and
      // its "radial" component is Ey. There is no Ez: the rods are uniform
      // along the beam, and the hard-edge model has no longitudinal field.
      const { Ez: ex, Er: ey } = unit.fieldAt(x, y);
      const w = drive(t);
      return { Ex: w * ex, Ey: w * ey, Ez: 0 };
    },

    potentialAt(x, y, zl, t) {
      if (zl < 0 || zl > length) return 0;
      return drive(t) * unit.potentialAt(x, y);
    },

    strikes(x, y, zl) {
      if (Math.hypot(x, y) > extent) return true;
      if (zl < 0 || zl > length) return false;
      // Field.strikes takes (transverse, height, axial), which is the reverse
      // of fieldAt(axial, transverse) above - so with the grid's axial axis
      // carrying x and its transverse axis y, the arguments go (y, 0, x).
      return unit.strikes(y, 0, x);
    },

    /**
     * Mathieu stability parameters for a given ion.
     *
     *     a = 8 z e U / (m r0^2 Omega^2)
     *     q = 4 z e V / (m r0^2 Omega^2)
     *
     * These are what actually decide whether an ion is transmitted. For an
     * RF-only filter (a = 0) the first stability region ends at q = 0.90803,
     * which is the high-mass cut-off; scanning U and V along a line just
     * inside the region's apex is how a quadrupole becomes a mass filter.
     */
    mathieu(massAmu, chargeStates) {
      const m = massAmu * ATOMIC_MASS_UNIT;
      const q = chargeStates * ELEMENTARY_CHARGE;
      const w = omega();
      const denom = m * r0 * r0 * w * w;
      if (denom === 0) return { a: 0, q: 0 };
      return {
        a: (8 * q * p.dcVoltage) / denom,
        q: (4 * q * p.rfAmplitude) / denom,
      };
    },

    // Drawn in the z-x plane. The x rods lie in that plane and are solid; the
    // y rods are perpendicular to it and are drawn faintly, because they are
    // there and act on the ion but are not in the slice being shown.
    rects: [
      { z0: 0, z1: length, r0: centre - rod, r1: centre + rod },
      { z0: 0, z1: length, r0: 0, r1: rod, ghost: true },
    ],
  };
}

/** High-mass cut-off of the first Mathieu stability region at a = 0. */
export const MATHIEU_Q_LIMIT = 0.90803;

/**
 * A working point well inside the first stability region.
 *
 * Not the apex. The apex is where a filter is operated when resolution is the
 * goal, and it is precisely where transmission is most fragile; 0.38 is
 * comfortably confining and comfortably clear of the q = 0.908 cut-off.
 * Measured on a 150 mm, 2 MHz filter with a 100 u, 50 eV beam: 9 of 9 at
 * q = 0.38 and at 0.23, 3 of 9 at 0.46, nothing at 1.83.
 */
export const MATHIEU_Q_WORKING = 0.38;

/**
 * The RF amplitude that puts a given ion at a given Mathieu q.
 *
 *     q = 4 z e V / (m r0^2 Omega^2)   =>   V = q m r0^2 Omega^2 / (4 z e)
 *
 * The inverse of `mathieu`, and the reason a quadrupole placed from the
 * toolbar confines the ion in the source rather than whichever ion the
 * defaults were written for. A filter's stability depends on mass, so a fixed
 * default amplitude is stable only by coincidence: at the previous defaults a
 * 100 u ion sat at q = 1.83, outside the first region entirely, and every ion
 * was lost.
 */
export function amplitudeForQ(params, massAmu, chargeStates = 1, q = MATHIEU_Q_WORKING) {
  const p = { ...QUADRUPOLE_DEFAULTS, ...params };
  const m = massAmu * ATOMIC_MASS_UNIT;
  const z = Math.abs(chargeStates) * ELEMENTARY_CHARGE;
  const r0 = mmToM(p.fieldRadius);
  const w = 2 * Math.PI * p.frequency * 1e6;
  if (z === 0) return 0;
  return (q * m * r0 * r0 * w * w) / (4 * z);
}
