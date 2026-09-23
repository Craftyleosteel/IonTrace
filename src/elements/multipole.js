/**
 * Multipole ion guide or trap - 2n rods driven in alternating RF phase.
 *
 * The quadrupole is the n = 2 member of this family and has its own element,
 * because it is used for something different: as a MASS FILTER, where what
 * matters is the Mathieu stability boundary and the fact that only a narrow
 * band of m/z survives. A hexapole or octopole is used for the opposite
 * purpose - to guide or trap ions of every mass with as little discrimination
 * as possible - and the quantity that governs that is not a stability
 * boundary but the depth of an effective potential well.
 *
 * Why higher multipoles guide better
 * ----------------------------------
 * The ideal 2n-pole potential is
 *
 *     phi = V (r/r0)^n cos(n theta),
 *
 * so the field amplitude grows as r^(n-1). For a quadrupole that is linear in
 * r; for an octopole it goes as r^3. The effective potential an ion feels when
 * the drive is fast compared with its motion is the Dehmelt result
 *
 *     U* = q^2 E0^2 / (4 m Omega^2),
 *
 * which therefore goes as r^(2n-2): flat across the middle of the guide and
 * rising steeply near the rods. That flat bottom is the point. An ion spends
 * most of its time where the RF field is weak, so it is heated less, its
 * transverse energy is less coupled to the drive, and the guide accepts a much
 * wider range of masses than a quadrupole can.
 *
 * The cost is that the confining force is weak near the axis, so a multipole
 * holds ions loosely - which is why real ones are usually operated with buffer
 * gas, and why they make good traps and poor filters.
 *
 * What is and is not modelled
 * ---------------------------
 * The field is solved, not assumed: real round rods are painted at the
 * requested radius and the Laplace solution includes whatever higher-order
 * terms that geometry produces. The effective potential is then computed from
 * that solved field rather than from the ideal form, so the number reported is
 * for the guide as drawn.
 *
 * The rods have a hard edge, like the quadrupole's, and for the same reason
 * (docs/PHYSICS.md section 8.6). There is no buffer gas, which for a device
 * usually operated with one is the omission that matters most - see the note
 * in the interface.
 */

import { PotentialArray, PLANAR } from '../grid.js';
import { solveBasis } from '../laplace.js';
import { Field } from '../field.js';
import { mmToM, ELEMENTARY_CHARGE, ATOMIC_MASS_UNIT } from '../constants.js';

export const MULTIPOLE_DEFAULTS = {
  poles: 8, // total rods; 6 = hexapole, 8 = octopole, ...
  fieldRadius: 5, // mm, r0 - axis to rod surface
  rodRadius: 2, // mm - eight of these clear each other by 1.4 mm at r0 = 5
  length: 120, // mm
  rfAmplitude: 300, // V, zero-to-peak on each rod
  frequency: 2.0, // MHz
  phase: 0, // degrees
  housingRadius: 14, // mm
  gridStep: 0.2, // mm, in the transverse plane
};

/**
 * Depth of the effective potential well, in electron-volts.
 *
 * Dehmelt's result for a fast drive: an ion in an oscillating field of
 * amplitude E0 feels a time-averaged potential energy
 *
 *     U* = q^2 E0^2 / (4 m Omega^2).
 *
 * Evaluated here at the field radius from the SOLVED field, so it accounts for
 * the real rod shape rather than assuming the ideal multipole.
 *
 * The approximation behind it is that the drive is fast compared with the
 * ion's own motion - formally that the Mathieu q is small. Above q of about
 * 0.3 the effective potential stops being a good description and a trap
 * designed on it will not behave as its depth suggests, so the element warns.
 */
export function wellDepth(fieldAmplitude, massAmu, chargeStates, frequencyMHz) {
  const m = massAmu * ATOMIC_MASS_UNIT;
  const q = Math.abs(chargeStates) * ELEMENTARY_CHARGE;
  const omega = 2 * Math.PI * frequencyMHz * 1e6;
  if (m === 0 || omega === 0) return 0;
  // Joules, then expressed in eV.
  return (q * q * fieldAmplitude * fieldAmplitude) / (4 * m * omega * omega) / ELEMENTARY_CHARGE;
}

export function createMultipole(params = {}, solverOpts = {}) {
  const p = { ...MULTIPOLE_DEFAULTS, ...params };
  const warnings = [];

  const poles = Math.max(4, Math.round(p.poles / 2) * 2);
  const order = poles / 2; // n: 2 = quadrupole, 3 = hexapole, 4 = octopole
  const r0 = mmToM(p.fieldRadius);
  const rod = mmToM(p.rodRadius);
  const housing = mmToM(p.housingRadius);

  if (r0 + 2 * rod >= housing) {
    throw new Error('Rods do not fit inside the housing; reduce the rod radius');
  }
  // Adjacent rods must not touch, or the two phases are shorted together.
  const gapBetween = 2 * (r0 + rod) * Math.sin(Math.PI / poles) - 2 * rod;
  if (gapBetween <= 0) {
    throw new Error(
      `${poles} rods of ${p.rodRadius} mm will not fit around a ${p.fieldRadius} mm ` +
        'aperture without touching; use thinner rods or fewer of them'
    );
  }
  if (gapBetween < rod * 0.4) {
    warnings.push(
      `Only ${(gapBetween * 1e3).toFixed(2)} mm between adjacent rods. Real ones need ` +
        'clearance for the RF, and the field between them is poorly resolved here.'
    );
  }

  // Transverse grid, odd node count, for the same reason as the quadrupole's:
  // the enclosure has to share the rods' symmetry or the guide is lopsided.
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
  // Both ends of this grid are housing, not beam apertures: the beam runs
  // perpendicular to the plane being solved.
  grid.openFaces = { zMin: false, zMax: false };

  const encl = grid.addElectrode('housing');
  const poleA = grid.addElectrode('poleA');
  const poleB = grid.addElectrode('poleB');
  grid.paintEnclosure(encl);

  // Rod centres sit one rod radius outside the field radius, so the rod
  // SURFACE is tangent to r0 - that is what r0 means.
  const centre = r0 + rod;
  const disc = (cx, cy) => (x, y) => Math.hypot(x - cx, y - cy) <= rod + grid.step * 1e-6;

  let paintedA = 0;
  let paintedB = 0;
  for (let k = 0; k < poles; k++) {
    const a = (k * Math.PI) / order;
    const cx = centre * Math.cos(a);
    const cy = centre * Math.sin(a);
    // Alternating phase around the ring, which is what makes it a multipole
    // rather than 2n rods at the same potential.
    if (k % 2 === 0) paintedA += grid.paint(poleA, disc(cx, cy));
    else paintedB += grid.paint(poleB, disc(cx, cy));
  }
  if (paintedA === 0 || paintedB === 0) {
    throw new Error('Multipole rods covered no grid nodes; check the geometry');
  }

  const { basis, reports } = solveBasis(grid, solverOpts);
  if (reports.some((r) => !r.converged)) {
    warnings.push('Laplace solve did not converge; this field is not trustworthy.');
  }

  // One unit solution, +1 V on one phase and -1 V on the other. Everything
  // this element does is that map scaled by the drive, so changing the
  // amplitude never re-solves.
  const unit = new Field(grid, basis);
  unit.setVoltages({ housing: 0, poleA: 1, poleB: -1 });

  const length = mmToM(p.length);
  const omega = () => 2 * Math.PI * p.frequency * 1e6;
  const drive = (t) => p.rfAmplitude * Math.cos(omega() * t + (p.phase * Math.PI) / 180);

  /** Field amplitude at the field radius, averaged around the aperture. */
  const rimField = () => {
    let sum = 0;
    const n = 64;
    for (let k = 0; k < n; k++) {
      const a = (2 * Math.PI * k) / n;
      const { Ez, Er } = unit.fieldAt(r0 * Math.cos(a), r0 * Math.sin(a));
      sum += Math.hypot(Ez, Er);
    }
    return (sum / n) * Math.abs(p.rfAmplitude);
  };

  return {
    type: 'multipole',
    // The kind of thing, not this one's pole count - that goes in the summary
    // line, where it changes without the element appearing to become a
    // different sort of element.
    label: 'Multipole guide',
    poles,
    params: p,
    length,
    bore: r0,
    clearBore: r0,
    outerRadius: extent,
    lengthScale: grid.step,
    warnings,
    grid,
    field: unit,
    order,

    get shortestPeriod() {
      return p.rfAmplitude === 0 || p.frequency === 0 ? null : 1 / (p.frequency * 1e6);
    },

    contains(x, y, zl) {
      return zl >= 0 && zl <= length;
    },

    fieldAt(x, y, zl, t) {
      if (zl < 0 || zl > length) return { Ex: 0, Ey: 0, Ez: 0 };
      // The grid's axes are x and y, so its "axial" component is Ex and its
      // "radial" one Ey. No Ez: the rods are uniform along the beam.
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
      // Field.strikes takes (transverse, height, axial), the reverse of
      // fieldAt(axial, transverse).
      return unit.strikes(y, 0, x);
    },

    /**
     * How deep the effective well is for a given ion, in electron-volts, and
     * whether the approximation it rests on is valid.
     */
    trapping(massAmu, chargeStates = 1) {
      const depth = wellDepth(rimField(), massAmu, chargeStates, p.frequency);
      // The Mathieu q of the equivalent quadrupole, as a validity check: the
      // effective potential is a fast-drive approximation and stops describing
      // the motion once q is no longer small.
      const m = massAmu * ATOMIC_MASS_UNIT;
      const z = Math.abs(chargeStates) * ELEMENTARY_CHARGE;
      const w = omega();
      const q = m && w ? (4 * z * Math.abs(p.rfAmplitude)) / (m * r0 * r0 * w * w) : 0;
      return { depth, q, valid: q < 0.3 };
    },

    /**
     * Rods, drawn in the plane the beam travels in.
     *
     * Only the rods that actually cross that plane are solid; the rest are
     * ghosted, because they are there and act on the ion but are not in the
     * slice being shown. With eight or twelve rods most of them are ghosts,
     * which is honest - a cross-section of a multipole in the beam plane is
     * mostly empty space.
     */
    get rects() {
      const out = [];
      for (let k = 0; k < poles; k++) {
        const a = (k * Math.PI) / order;
        const cy = centre * Math.cos(a); // distance from the axis in the drawn plane
        const inPlane = Math.abs(Math.sin(a)) < 1e-9;
        out.push({
          z0: 0,
          z1: length,
          r0: Math.abs(cy) - rod,
          r1: Math.abs(cy) + rod,
          ghost: !inPlane,
        });
      }
      return out;
    },
  };
}
