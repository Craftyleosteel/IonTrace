/**
 * Quadrupole deflector - four curved electrodes that turn the beam through a
 * right angle.
 *
 * This is the device a real beamline uses to fold its path: four electrodes
 * filling the corners of a grounded box, with concave faces forming a
 * near-circular aperture, held at +V, -V, +V, -V around the square. The beam
 * enters through the gap on one axis and leaves through the gap on the
 * perpendicular one.
 *
 * Why it is not a sector bender
 * -----------------------------
 * A cylindrical sector deflects by pushing the beam radially inward with a
 * field that is everywhere perpendicular to the orbit. A quadrupole deflector
 * does something else entirely: it sets up a two-dimensional QUADRUPOLE field
 * IN the bend plane and lets it couple the two axes. Writing the potential as
 *
 *     phi = C x z
 *
 * with the beam entering along z and leaving along x, the equations of motion
 * are coupled rather than central:
 *
 *     m x'' = -qC z        m z'' = -qC x
 *
 * so the axial velocity is traded for transverse velocity as the ion crosses.
 * The trajectory is NOT a circular arc.
 *
 * Diagonalising with u = x + z and w = x - z separates them:
 *
 *     u'' = -b u     (oscillatory)
 *     w'' = +b w     (exponential),       b = qC/m
 *
 * The exponential branch is the reason this device is touchy: half of the
 * motion grows rather than oscillates, so an ion that is off-energy or
 * off-axis diverges from the design orbit instead of merely lagging it.
 *
 * The matched voltage
 * -------------------
 * Requiring the ion to arrive at the exit aperture with BOTH the right
 * position and the right direction fixes the design point. Writing s for half
 * the transit phase, the oscillatory branch demands v/(a.omega) = cot(s) and
 * the exponential branch demands v/(a.omega) = tanh(s), so
 *
 *     cot(s) = tanh(s)      =>      s = 0.937552,  theta = 1.875104
 *
 * whence (a.omega/v)^2 = coth^2(s) = 1.85565 and
 *
 *     V0 = 1.8556 (T/q) (r0/a)^2
 *
 * with T the kinetic energy, r0 the aperture radius and a the half-width of
 * the field region.
 *
 * There is a second root, and it is a trap
 * ----------------------------------------
 * cot(s) = tanh(s) has further roots, the next near s = 2.347 (k = 0.9641),
 * and direct integration of the coupled equations confirms it is also a clean
 * ninety degrees - the ion simply takes a longer way round, bending the other
 * way. It is not usable here. On that branch the oscillatory amplitude is
 * 1.427 a rather than 1.241 a, which carries the orbit out to rho = 1.01 a:
 * past the electrode faces at r0, which is necessarily less than a. The ion
 * lands on an electrode instead of reaching the exit. The short branch peaks
 * at rho = 0.88 a and clears the aperture with room to spare, so that is the
 * one this element is designed around.
 *
 * What the formula assumes, and therefore what to expect
 * ------------------------------------------------------
 * It assumes the ideal quadrupole potential everywhere inside the box and
 * nothing outside it. The solved field is neither: real electrodes subtend
 * finite arcs, the grounded box shapes the field near the apertures, and the
 * field does not stop abruptly at the entrance plane. Measured against the
 * solved field, with V/V0 the voltage as a multiple of the ideal one and a
 * nine-ion beam of 1.5 mm radius:
 *
 *     r0/a    turns exactly 90 deg at   transmits 9/9 over
 *     0.655           V/V0 ~ 1.00         0.75 .. 1.05
 *     0.905           V/V0 ~ 1.02         0.85 .. 1.05
 *     0.950           V/V0 ~ 1.07         0.90 .. 1.10
 *
 * So the formula is a good predictor of the right angle at every proportion,
 * drifting a few per cent high as the electrodes thin and the box moves in.
 * The window is at least a tenth either side throughout, which is what the
 * tolerance shown in the interface is based on. It is asymmetric, and which
 * way it leans depends on the geometry - thick electrodes tolerate too little
 * voltage, thin ones tolerate too much - which is the kind of thing worth
 * finding with the tuner rather than reasoning about.
 *
 * So the matched voltage is a starting point of the right size, not the final
 * answer - which is exactly why a real deflector of this kind is followed by a
 * correction lens, and why this element ships with a tuner that searches the
 * solved field for the voltage that actually maximises transmission.
 */

import { PotentialArray, PLANAR } from '../grid.js';
import { solveBasis } from '../laplace.js';
import { Field } from '../field.js';
import { mmToM } from '../constants.js';
import { yawFrame, compose, translation, rollFrame, inverse } from '../frames.js';

export const BENDER_DEFAULTS = {
  apertureRadius: 19, // mm, r0 - centre to the concave electrode faces
  electrodeThickness: 9, // mm, radial thickness of each electrode
  boxClearance: 1, // mm, gap between electrode backs and the grounded box
  gapAngle: 16, // degrees of clear aperture either side of each axis
  height: 30, // mm, aperture perpendicular to the bend plane
  voltage: 0, // V on each electrode (+V and -V on the diagonals)
  bendPlane: 0, // degrees of roll: 0 bends horizontally, 90 vertically
  gridStep: 0.5, // mm, in the bend plane
};

/** The first root of cot(s) = tanh(s), and the constant it gives. */
export const DEFLECTOR_DESIGN_PHASE = 0.937552;
export const DEFLECTOR_CONSTANT = Math.tanh(DEFLECTOR_DESIGN_PHASE) ** -2; // 1.85565

/**
 * Ideal matched voltage: V0 = k (T/q) (r0/a)^2.
 *
 * `a` is the half-width of the field region - the distance from the centre to
 * the entrance aperture - not the electrode radius. The two differ by the
 * electrode thickness and the box clearance, and the ratio enters squared, so
 * the distinction is worth a factor of two in practice.
 */
export function matchedVoltage(params, energyEV, chargeStates = 1) {
  const p = { ...BENDER_DEFAULTS, ...params };
  const r0 = p.apertureRadius;
  const a = r0 + p.electrodeThickness + p.boxClearance;
  const q = Math.abs(chargeStates);
  if (q === 0 || a === 0) return 0;
  return DEFLECTOR_CONSTANT * (energyEV / q) * (r0 / a) ** 2;
}

export function createBender(params = {}, solverOpts = {}) {
  const p = { ...BENDER_DEFAULTS, ...params };
  const warnings = [];

  const r0 = mmToM(p.apertureRadius);
  const outer = mmToM(p.apertureRadius + p.electrodeThickness);
  const a = mmToM(p.apertureRadius + p.electrodeThickness + p.boxClearance);
  const halfHeight = mmToM(p.height) / 2;
  const gap = (p.gapAngle * Math.PI) / 180;

  if (p.gapAngle >= 45) {
    throw new Error('Aperture gaps of 45 degrees or more leave no electrode');
  }
  if (p.gapAngle < 5) {
    warnings.push(
      `Aperture gaps of ${p.gapAngle} degrees are narrow; the beam has little ` +
        'clearance entering and leaving.'
    );
  }

  // Bend plane, as a roll about the beam direction. One solve serves every
  // plane: the field is the same, seen from a rotated frame.
  const roll = (p.bendPlane * Math.PI) / 180;
  const cosR = Math.cos(roll);
  const sinR = Math.sin(roll);
  const intoBend = (x, y) => [cosR * x + sinR * y, -sinR * x + cosR * y];
  const outOfBend = (u, v) => [cosR * u - sinR * v, sinR * u + cosR * v];

  /**
   * The solve lives in the BEND PLANE, and the beam travels inside it.
   *
   * This is the only element here where that is so. A lens is solved in the
   * r-z plane it is symmetric about; a mass filter in the plane transverse to
   * the beam. Here the two-dimensional field and the trajectory occupy the
   * SAME plane, which is precisely why the axes couple.
   *
   * Node count forced odd so the box is exactly symmetric about the centre -
   * without that the four electrodes do not see the same enclosure and the
   * deflector is lopsided.
   */
  const step = mmToM(p.gridStep);
  const half = Math.max(6, Math.round(a / step));
  const extent = half * step;
  const grid = new PotentialArray({
    nz: 2 * half + 1,
    nr: 2 * half + 1,
    step,
    symmetry: PLANAR,
    z0: -extent,
    r0: -extent,
  });

  // The grid's z axis is the transverse bend coordinate, not the beam axis,
  // so its end faces are box walls rather than apertures. Leaving them open -
  // the default, which suits every element solved in the r-z plane - would
  // make the entrance wall invisible to the collision test. The real holes in
  // the box are carved out by `inAperture` below.
  grid.openFaces = { zMin: false, zMax: false };

  const box = grid.addElectrode('box');
  const poleA = grid.addElectrode('poleA'); // the x.z > 0 diagonal
  const poleB = grid.addElectrode('poleB'); // the x.z < 0 diagonal
  grid.paintEnclosure(box);

  /**
   * One electrode: a quadrant of the annulus between r0 and `outer`, stopping
   * `gapAngle` short of each axis so the beam has somewhere to enter and
   * leave. The concave face is what makes the field quadrupolar near the
   * centre.
   *
   * The grid's "z" axis carries the axial coordinate Z and its "r" axis the
   * transverse X, both measured from the deflector centre.
   */
  const quadrant = (sx, sz) => (Z, X) => {
    const rho = Math.hypot(X, Z);
    if (rho < r0 - step * 1e-6 || rho > outer + step * 1e-6) return false;
    if (Math.sign(X) !== sx || Math.sign(Z) !== sz) return false;
    // Angle away from the nearer axis, so the gap is symmetric about both.
    const psi = Math.atan2(Math.abs(X), Math.abs(Z));
    return psi > gap && psi < Math.PI / 2 - gap;
  };

  // Sign convention: a POSITIVE voltage bends a POSITIVE ion towards -x. That
  // choice is what ties the electrode polarity to the exit face, and the two
  // are not independent - with the diagonals swapped the same ion turns the
  // other way and leaves through the +x wall instead, which this element has
  // no aperture in. A device wired the other way round, as the published one
  // in the reference drawing is, is the mirror image of this and identical
  // physics.
  const paintedA =
    grid.paint(poleA, quadrant(1, 1)) + grid.paint(poleA, quadrant(-1, -1));
  const paintedB =
    grid.paint(poleB, quadrant(-1, 1)) + grid.paint(poleB, quadrant(1, -1));
  if (paintedA === 0 || paintedB === 0) {
    throw new Error('Deflector electrodes covered no grid nodes; check the geometry');
  }

  const { basis, reports } = solveBasis(grid, solverOpts);
  if (reports.some((r) => !r.converged)) {
    warnings.push('Laplace solve did not converge; this field is not trustworthy.');
  }

  // Unit solution, per volt on each electrode. Everything the element does is
  // this map scaled, so changing the voltage never re-solves.
  const unit = new Field(grid, basis);
  unit.setVoltages({ box: 0, poleA: -1, poleB: 1 });

  // Nominal path: a quarter turn of radius `a` from the entrance aperture to
  // the exit one. The REAL trajectory is not this curve - it is not a
  // circular arc at all - but the nominal orbit is what the element's
  // placement and the drawn reference line are built on.
  const length = (Math.PI / 2) * a;

  /** Local point -> deflector-centred bend coordinates. */
  function centred(x, y, z) {
    const [bx, by] = intoBend(x, y);
    return { X: bx, Z: z - a, h: by };
  }

  /**
   * The beam holes in the grounded box.
   *
   * The solve paints the box all the way round, because a grounded plane with
   * a hole in it is very nearly a grounded plane and solving it that way keeps
   * the Dirichlet boundary closed. But the holes are real, and an ion passing
   * through one must not be counted as hitting metal. Without this the failure
   * is total rather than subtle: every ion is destroyed the instant it reaches
   * the entrance plane, at every voltage.
   *
   * The hole is as wide as the channel it feeds: the electrodes stop `gapAngle`
   * short of each axis, so the clear width where the channel meets the box is
   * `outer sin(gap)` either side. The beam enters through the -Z face and
   * leaves through the -X one.
   */
  const holeHalf = outer * Math.sin(gap);
  const rim = step * 1.5;
  const inAperture = (X, Z) =>
    (Z < -extent + rim && Math.abs(X) < holeHalf) ||
    (X < -extent + rim && Math.abs(Z) < holeHalf);

  return {
    type: 'bender',
    label: 'Quadrupole deflector',
    params: p,
    length,
    bore: r0,
    outerRadius: Math.max(extent, halfHeight),
    lengthScale: step,
    shortestPeriod: null,
    warnings,
    grid,
    field: unit,
    curved: true,

    /** A quarter turn, conjugated by the roll so the beam is not twisted. */
    exitTransform: compose(
      compose(rollFrame(roll), compose(translation(-a, 0, a), yawFrame(Math.PI / 2))),
      inverse(rollFrame(roll))
    ),

    pathPoint(f) {
      const t = (f * Math.PI) / 2;
      const [x, y] = outOfBend(-a + a * Math.cos(t), 0);
      return [x, y, a * Math.sin(t)];
    },

    setVoltage(v) {
      // Nothing to recompute: fieldAt scales the stored unit map.
      p.voltage = v;
    },

    contains(x, y, z) {
      const { X, Z, h } = centred(x, y, z);
      // A tolerance of one grid step: the entrance plane sits exactly on the
      // box edge, and without it an ion arriving there would belong to no
      // element at all and be reported as having wandered out of the column.
      const edge = extent + step;
      return Math.abs(X) <= edge && Math.abs(Z) <= edge && Math.abs(h) <= halfHeight * 1.5;
    },

    fieldAt(x, y, z) {
      const { X, Z } = centred(x, y, z);
      if (Math.abs(X) > extent || Math.abs(Z) > extent) return { Ex: 0, Ey: 0, Ez: 0 };
      // The grid's axes are (Z, X); there is no field along the height for an
      // ideal deflector, because the electrodes are uniform in that direction.
      const { Ez: eZ, Er: eX } = unit.fieldAt(Z, X);
      const [ex, ey] = outOfBend(p.voltage * eX, 0);
      return { Ex: ex, Ey: ey, Ez: p.voltage * eZ };
    },

    potentialAt(x, y, z) {
      const { X, Z } = centred(x, y, z);
      if (Math.abs(X) > extent || Math.abs(Z) > extent) return 0;
      return p.voltage * unit.potentialAt(Z, X);
    },

    strikes(x, y, z) {
      const { X, Z, h } = centred(x, y, z);
      if (Math.abs(h) > halfHeight) return true;
      if (inAperture(X, Z)) return false;
      if (Math.abs(X) > extent || Math.abs(Z) > extent) return true;
      // Field.strikes takes (transverse, height, axial). The grid's axial axis
      // carries Z and its transverse axis carries X, so they go in that order -
      // the reverse of fieldAt(Z, X), which takes (axial, transverse).
      return unit.strikes(X, 0, Z);
    },

    /** Ideal matched voltage for this ion, in volts on each electrode. */
    matchedVoltage(energyEV, chargeStates) {
      return matchedVoltage(p, energyEV, chargeStates);
    },

    /**
     * Four curved electrodes plus the box, as polygons in the element's own
     * frame. The concave faces are drawn as real arcs, because their shape is
     * what makes the field quadrupolar rather than a detail of the picture.
     */
    shapes() {
      // A point in the deflector's own centred coordinates, put back into the
      // element's frame - rolled into the bend plane, and shifted so the
      // centre sits one half-width along the beam.
      const corner = (X, Z) => {
        const [x, y] = outOfBend(X, 0);
        return [x, y, a + Z];
      };

      const span = Math.PI / 2 - 2 * gap;
      const arcSteps = 10;
      const out = [];

      for (const [sx, sz] of [[1, 1], [-1, 1], [1, -1], [-1, -1]]) {
        const arc = (rho, reverse) => {
          const pts = [];
          for (let k = 0; k <= arcSteps; k++) {
            const f = reverse ? 1 - k / arcSteps : k / arcSteps;
            const psi = gap + f * span;
            pts.push(corner(sx * rho * Math.sin(psi), sz * rho * Math.cos(psi)));
          }
          return pts;
        };
        // Down the concave face, back along the outside.
        out.push({ points: [...arc(r0, false), ...arc(outer, true)] });
      }

      // The grounded box, as four thin walls.
      const t = Math.max(step * 2, a * 0.03);
      const wall = (x0, x1, z0, z1) => ({
        points: [corner(x0, z0), corner(x1, z0), corner(x1, z1), corner(x0, z1)],
        wall: true,
      });
      out.push(wall(-extent, extent, extent - t, extent)); // far side
      out.push(wall(extent - t, extent, -extent, extent)); // far side
      // The entrance (-Z) and exit (-X) faces are drawn in two pieces, with the
      // beam hole between them, so the picture shows the same apertures the
      // collision test uses.
      out.push(wall(-extent, -holeHalf, -extent, -extent + t));
      out.push(wall(holeHalf, extent, -extent, -extent + t));
      out.push(wall(-extent, -extent + t, -extent, -holeHalf));
      out.push(wall(-extent, -extent + t, holeHalf, extent));

      return out;
    },
  };
}
