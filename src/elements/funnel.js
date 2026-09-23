/**
 * Ion funnel - a stack of ring electrodes with a shrinking aperture.
 *
 * Adjacent rings are driven in opposite RF phase, which builds a repulsive
 * effective potential along the inside face of the stack, and a DC gradient
 * down the stack pushes ions along it. Because the apertures narrow, the
 * confining wall closes in as the ions travel, so a wide diffuse cloud enters
 * and a narrow beam leaves. That is what a funnel is for: collecting ions out
 * of a gas-filled region and handing them to the vacuum optics downstream.
 *
 * Axisymmetric, so it solves in the r-z plane like a lens, and the ring stack
 * is just a lot of annuli painted at decreasing inner radius.
 *
 * Two fields, one basis
 * ---------------------
 * A funnel needs two independent voltage patterns at once: a DC ramp along the
 * stack and an RF alternation across it. Both are weighted sums of the same
 * per-ring unit solutions, so the Laplace problem is solved once per ring and
 * then read twice -
 *
 *     phi(r, z, t) = phi_DC(r, z) + cos(Omega t) phi_RF(r, z)
 *
 * with each term a fast adjust over the shared basis. Changing the gradient or
 * the amplitude re-weights; it never re-solves. That is why every ring gets
 * its own electrode rather than being lumped into two groups: lumping them
 * would make the RF cheap and the DC ramp impossible.
 *
 * THE OMISSION THAT MATTERS
 * -------------------------
 * A real ion funnel works at one to thirty millibar, and the gas is not
 * incidental - it is half the mechanism. Collisions damp the ions' transverse
 * energy so they settle into the effective potential minimum instead of
 * rattling across it, and they are what lets a funnel collect a warm, diffuse
 * cloud at all. IonTrace has no buffer gas (docs/PHYSICS.md section 5).
 *
 * So what this element reproduces is the FIELD of a funnel and the guiding a
 * cold, near-axis beam would get from it. What it does not reproduce is the
 * collection of a hot one, which is the job funnels exist to do. Transmission
 * here is optimistic for a cold beam and meaningless for a warm one, and the
 * element says so rather than letting a number stand unqualified.
 */

import { PotentialArray, CYLINDRICAL } from '../grid.js';
import { solveBasis } from '../laplace.js';
import { Field } from '../field.js';
import { mmToM, ELEMENTARY_CHARGE, ATOMIC_MASS_UNIT } from '../constants.js';

export const FUNNEL_DEFAULTS = {
  rings: 14, // number of ring electrodes
  entryRadius: 12, // mm, aperture of the first ring
  exitRadius: 2, // mm, aperture of the last
  ringThickness: 0.5, // mm, along the beam
  pitch: 2.5, // mm, ring centre to ring centre
  /*
    A real funnel runs near 0.5-1 MHz at 100-200 V peak-to-peak, and gets away
    with it because the gas removes the energy the RF puts in. In vacuum those
    settings do not work here and should not: at 0.7 MHz a 100 u ion crosses
    this stack in under six RF cycles, far too few for the effective potential
    to mean anything, and it is thrown into a ring. Measured, 0 of 9 through.

    Worse - and this is the honest core of it - once the drive IS fast enough,
    a DEEPER well transmits WORSE. At 4 MHz: a 1.1 eV well passes 9 of 9, a
    4.3 eV well 3 of 9, a 17 eV well none at all. The well confines, but the
    heating that comes with it has nowhere to go, and the ion climbs out. In a
    real funnel collisions carry that energy away continuously; that is not a
    detail of the device, it is the mechanism.

    So these defaults are chosen to work in vacuum, which is a different regime
    from the one the hardware was designed for. The readout says so.
  */
  rfAmplitude: 100, // V, zero-to-peak, alternating ring to ring
  frequency: 4.0, // MHz
  phase: 0, // degrees
  dcEntry: 20, // V on the first ring
  dcExit: 0, // V on the last
  housingRadius: 18, // mm
  gridStep: 0.25, // mm
};

export function createFunnel(params = {}, solverOpts = {}) {
  const p = { ...FUNNEL_DEFAULTS, ...params };
  const warnings = [];

  const rings = Math.max(3, Math.round(p.rings));
  const pitch = mmToM(p.pitch);
  const thick = mmToM(p.ringThickness);
  const housing = mmToM(p.housingRadius);
  const entryR = mmToM(p.entryRadius);
  const exitR = mmToM(p.exitRadius);

  if (p.entryRadius >= p.housingRadius) {
    throw new Error('The entry aperture must be smaller than the housing');
  }
  if (p.exitRadius <= 0) throw new Error('The exit aperture must be open');
  if (p.ringThickness >= p.pitch) {
    throw new Error('Rings of that thickness would touch; increase the pitch');
  }
  if (p.gridStep > p.ringThickness) {
    warnings.push(
      `Grid step (${p.gridStep} mm) does not resolve a ${p.ringThickness} mm ring. ` +
        'The stack is being approximated by something coarser than it is.'
    );
  }

  // A margin at each end so the fringe of the stack is inside the domain
  // rather than pressed against a grounded face.
  const margin = Math.max(pitch * 2, entryR);
  const stack = pitch * (rings - 1) + thick;
  const totalLength = stack + 2 * margin;

  const step = mmToM(p.gridStep);
  const nz = Math.round(totalLength / step) + 1;
  const nr = Math.round(housing / step) + 1;
  const grid = new PotentialArray({ nz, nr, step, symmetry: CYLINDRICAL });

  const encl = grid.addElectrode('housing');
  grid.paintEnclosure(encl);

  /** Inner radius of ring k, shrinking linearly down the stack. */
  const apertureOf = (k) => entryR + ((exitR - entryR) * k) / (rings - 1);

  const ids = [];
  for (let k = 0; k < rings; k++) {
    const id = grid.addElectrode(`ring${k}`);
    // The two end faces belong to the nearest ring rather than to the
    // grounded housing.
    //
    // This is not a detail. A funnel works by a DC gradient down the stack,
    // and grounded end caps clamp the potential to zero at BOTH ends - so any
    // ramp painted on the rings becomes a hill that rises into the stack and
    // falls back out of it. Measured with the caps grounded: a 20 V ramp put a
    // 13 V barrier in front of the entrance, and a 5 eV beam was reflected in
    // its entirety, 0 of 9 through at every beam width.
    //
    // A real funnel is not a stack sitting between two earthed plates; it is
    // embedded in a line that continues the gradient. Giving each end face the
    // potential of the ring beside it says exactly that, and keeps the
    // Dirichlet problem closed.
    if (k === 0) grid.paint(id, (z, r, i) => i === 0);
    if (k === rings - 1) grid.paint(id, (z, r, i) => i === grid.nz - 1);
    const zA = margin + k * pitch;
    const zB = zA + thick;
    const inner = apertureOf(k);
    const painted = grid.paint(
      id,
      (z, r) => grid.spans(z, zA, zB) && r >= inner - grid.step * 1e-6
    );
    if (painted === 0) {
      throw new Error(`Ring ${k} covered no grid nodes; check the pitch and step`);
    }
    ids.push({ id, k, name: `ring${k}`, inner, zA, zB });
  }

  const { basis, reports } = solveBasis(grid, solverOpts);
  if (reports.some((r) => !r.converged)) {
    warnings.push('Laplace solve did not converge; this field is not trustworthy.');
  }

  /*
    Two readings of the one basis. `dc` carries the ramp down the stack; `rf`
    carries the alternation across it, at unit amplitude so the drive can scale
    it without re-weighting.
  */
  const dc = new Field(grid, basis);
  const rf = new Field(grid, basis);

  const applyVoltages = () => {
    const dcMap = { housing: 0 };
    const rfMap = { housing: 0 };
    for (const r of ids) {
      const f = rings === 1 ? 0 : r.k / (rings - 1);
      dcMap[r.name] = p.dcEntry + (p.dcExit - p.dcEntry) * f;
      rfMap[r.name] = r.k % 2 === 0 ? 1 : -1;
    }
    dc.setVoltages(dcMap);
    rf.setVoltages(rfMap);
  };
  applyVoltages();

  const length = grid.zLength;
  const omega = () => 2 * Math.PI * p.frequency * 1e6;
  const drive = (t) => p.rfAmplitude * Math.cos(omega() * t + (p.phase * Math.PI) / 180);

  return {
    type: 'funnel',
    label: 'Ion funnel',
    params: p,
    length,
    bore: exitR,
    // The clear aperture narrows down the stack, so the smallest one is the
    // only radius that is free along the WHOLE element - which is what this
    // promise has to mean. See `isFree` in beamline.js.
    clearBore: exitR,
    outerRadius: housing,
    lengthScale: grid.step,
    warnings,
    grid,
    field: rf,
    dcField: dc,
    rings,
    apertureOf,

    get shortestPeriod() {
      return p.rfAmplitude === 0 || p.frequency === 0 ? null : 1 / (p.frequency * 1e6);
    },

    /** Re-weight both maps. Fast adjust: no relaxation happens here. */
    setVoltage(v) {
      p.rfAmplitude = v;
      applyVoltages();
    },
    retune() {
      applyVoltages();
    },

    contains(x, y, zl) {
      return zl >= 0 && zl <= length;
    },

    fieldAt(x, y, zl, t) {
      const a = dc.fieldAt3D(x, y, zl);
      const b = rf.fieldAt3D(x, y, zl);
      const w = drive(t);
      return { Ex: a.Ex + w * b.Ex, Ey: a.Ey + w * b.Ey, Ez: a.Ez + w * b.Ez };
    },

    potentialAt(x, y, zl, t) {
      return dc.potentialAt3D(x, y, zl) + drive(t) * rf.potentialAt3D(x, y, zl);
    },

    strikes(x, y, zl) {
      return dc.strikes(x, y, zl);
    },

    /**
     * Depth of the RF effective potential at the rim of ring k, in eV.
     *
     * Dehmelt's fast-drive result U* = q^2 E0^2 / (4 m Omega^2), evaluated on
     * the solved RF field just inside a ring's aperture - which is where the
     * wall the ions are being held off actually is.
     */
    wellAt(k, massAmu, chargeStates = 1) {
      const r = ids[Math.max(0, Math.min(ids.length - 1, k))];
      const probe = Math.max(r.inner - grid.step * 2, grid.step);
      const { Ez, Er } = rf.fieldAt((r.zA + r.zB) / 2, probe);
      const amp = Math.hypot(Ez, Er) * Math.abs(p.rfAmplitude);
      const m = massAmu * ATOMIC_MASS_UNIT;
      const q = Math.abs(chargeStates) * ELEMENTARY_CHARGE;
      const w = omega();
      if (!m || !w) return 0;
      return (q * q * amp * amp) / (4 * m * w * w) / ELEMENTARY_CHARGE;
    },

    /** The rings, as annuli in the r-z plane. */
    get rects() {
      return ids.map((r) => ({ z0: r.zA, z1: r.zB, r0: r.inner, r1: housing }));
    },
  };
}
