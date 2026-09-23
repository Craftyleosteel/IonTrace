/**
 * Cylindrical sector bender - two curved plates that turn the beam.
 *
 * This is the element that makes a beamline stop being a line. Its exit faces
 * a different direction from its entrance, so everything downstream turns with
 * it, and the column becomes a path through space.
 *
 * Geometry and the reused solver
 * ------------------------------
 * Two coaxial cylindrical electrodes share a vertical axis, and the beam
 * travels between them on a circular arc of radius R. The potential does not
 * depend on the angle round that axis, so in the coordinates (rho, h) -
 * distance from the bend axis, and height above the bend plane - Laplace's
 * equation is
 *
 *     d2phi/drho2 + (1/rho) dphi/drho + d2phi/dh2 = 0
 *
 * which is EXACTLY the axisymmetric stencil the einzel lens uses. The same
 * solver applies unchanged; only the interpretation of the two axes differs.
 * For a lens, "z" is the beam axis and "r" is distance from it. Here "z" is
 * height and "r" is distance from a vertical bend axis, and the beam travels
 * in the ignorable third coordinate. That reuse is not a trick - it is the
 * same symmetry.
 *
 * The one thing that had to change is that the grid no longer contains the
 * axis: rho runs over a narrow band about R rather than starting at zero, so
 * the singular on-axis stencil does not apply and `r0 > 0` is now permitted
 * in cylindrical mode.
 *
 * Matching the beam
 * -----------------
 * For an ion of kinetic energy T and charge q to travel the arc of radius R,
 * the electric force must supply the centripetal one:
 *
 *     qE = m v^2 / R = 2T / R
 *
 * Between cylinders at rho1 and rho2 held at +/- V/2 the field on the central
 * orbit is E = V / (R ln(rho2/rho1)), so the plates must be at
 *
 *     V = (2T / q) ln(rho2 / rho1)
 *
 * which for a narrow gap d is close to the familiar 2Td/(qR). `matchedVoltage`
 * returns it, and the UI shows it beside the control, because a bender at the
 * wrong voltage does not bend the beam slightly wrongly - it puts it into a
 * plate.
 *
 * This is an ENERGY filter as much as a bender: the matched voltage depends on
 * T, so ions of the wrong energy are dispersed onto the electrodes. That is
 * what a hemispherical or 127-degree analyser is for.
 *
 * What is not modelled
 * --------------------
 *  - Fringe fields at the entrance and exit of the sector. As with the
 *    quadrupole, the field switches on abruptly, so transmission is
 *    optimistic and an off-orbit ion's potential energy jumps at the faces.
 *  - Vertical focusing from a real cylindrical bender comes only from the
 *    fringe fields, so this element has essentially none: a cylindrical
 *    sector focuses in the bend plane and not at all out of it. That is
 *    correct for the ideal geometry, and it is the reason spherical
 *    deflectors exist.
 */

import { PotentialArray, CYLINDRICAL } from '../grid.js';
import { solveBasis } from '../laplace.js';
import { Field } from '../field.js';
import { mmToM, ELEMENTARY_CHARGE, eVToJoules } from '../constants.js';
import { yawFrame, compose, translation, rollFrame, inverse } from '../frames.js';

export const BENDER_DEFAULTS = {
  bendRadius: 40, // mm, radius of the central orbit
  bendAngle: 90, // degrees
  // Which plane the bend happens in, as a roll about the beam direction.
  // 0 turns the beam horizontally, 90 turns it vertically. One element and
  // one solve serve every plane; see the note on conjugation below.
  bendPlane: 0, // degrees
  gap: 8, // mm, between the plates
  height: 16, // mm, aperture across the plates, perpendicular to the bend
  voltage: 0, // V across the plates; 0 means "use the matched value"
  gridStep: 0.3, // mm
};

/**
 * Plate voltage that puts an ion of `energyEV` and charge `chargeStates` on
 * the central orbit, in volts across the pair.
 */
export function matchedVoltage({ bendRadius, gap }, energyEV, chargeStates = 1) {
  const r1 = mmToM(bendRadius - gap / 2);
  const r2 = mmToM(bendRadius + gap / 2);
  const T = eVToJoules(energyEV);
  const q = Math.abs(chargeStates) * ELEMENTARY_CHARGE;
  if (q === 0) return 0;
  return ((2 * T) / q) * Math.log(r2 / r1);
}

export function createBender(params = {}, solverOpts = {}) {
  const p = { ...BENDER_DEFAULTS, ...params };
  const warnings = [];

  const R = mmToM(p.bendRadius);
  const gap = mmToM(p.gap);
  const alpha = (p.bendAngle * Math.PI) / 180;

  if (p.gap >= p.bendRadius) {
    throw new Error('Bender gap must be smaller than its bend radius');
  }
  if (p.gap / p.bendRadius > 0.5) {
    warnings.push(
      `Gap is ${(100 * p.gap / p.bendRadius).toFixed(0)} % of the bend radius. ` +
        'The field departs noticeably from uniform across a gap this wide, and ' +
        'the matched voltage only puts the central orbit right.'
    );
  }

  // Radial band: the plates plus a little clearance either side so the
  // enclosure is not sitting on them.
  const rInner = R - gap / 2;
  const rOuter = R + gap / 2;
  const pad = Math.max(gap * 0.25, mmToM(p.gridStep) * 2);
  const rMin = rInner - pad;
  const rMax = rOuter + pad;

  const step = mmToM(p.gridStep);
  const nr = Math.max(7, Math.round((rMax - rMin) / step) + 1);
  // Vertical: an odd count so the bend plane is a node and the geometry is
  // exactly symmetric about it, for the same reason the quadrupole needs it.
  const halfH = Math.max(3, Math.round(mmToM(p.height) / 2 / step));
  const nz = 2 * halfH + 1;

  const grid = new PotentialArray({
    nz,
    nr,
    step,
    symmetry: CYLINDRICAL,
    z0: -halfH * step,
    r0: rMin,
  });

  const lids = grid.addElectrode('lids');
  const inner = grid.addElectrode('inner');
  const outer = grid.addElectrode('outer');
  grid.paintEnclosure(lids);

  // The plates span the full height, one grid step thick in radius.
  const plate = (target) => (z, r) => grid.spans(r, target - step / 2, target + step / 2);
  const paintedIn = grid.paint(inner, plate(rInner));
  const paintedOut = grid.paint(outer, plate(rOuter));
  if (paintedIn === 0 || paintedOut === 0) {
    throw new Error('Bender plates covered no grid nodes; check the geometry');
  }

  const { basis, reports } = solveBasis(grid, solverOpts);
  if (reports.some((r) => !r.converged)) {
    warnings.push('Laplace solve did not converge; this field is not trustworthy.');
  }

  // One unit solution, so the stored map is per volt ACROSS the pair and
  // scaling it is a single multiply.
  //
  // The polarity is not arbitrary. A positive ion needs a centripetal force,
  // pointing INWARD towards the bend centre, so it needs E pointing inward -
  // and since E = -grad phi, the potential must INCREASE outward. The outer
  // plate is therefore the positive one. Getting this backwards does not bend
  // the beam slightly wrongly; it drives it straight into the outer plate,
  // which is exactly what the first version did at every voltage.
  const unit = new Field(grid, basis);
  unit.setVoltages({ lids: 0, inner: -0.5, outer: 0.5 });

  // Arc length of the central orbit. This is the element's "length" for the
  // purposes of path bookkeeping, though the beam does not travel in a
  // straight line through it.
  const length = R * alpha;

  /**
   * Local geometry.
   *
   * The beam enters at the local origin travelling along +z, and curves
   * towards -x. The centre of curvature therefore sits at local (-R, 0, 0),
   * and a point at arc angle theta lies at
   *
   *     (-R + R cos theta, y, R sin theta)
   *
   * so theta = 0 is the entrance and theta = alpha the exit.
   */
  const centre = [-R, 0, 0];

  /**
   * The bend plane, as a roll about the beam direction.
   *
   * The element is built to curve in its own x-z plane. Rolling it makes the
   * same element - the same solve, the same geometry - bend in any other
   * plane, so a vertical bender needs no separate implementation. Everything
   * below works in the ROLLED frame, and the two helpers convert.
   */
  const roll = (p.bendPlane * Math.PI) / 180;
  const cosR = Math.cos(roll);
  const sinR = Math.sin(roll);

  /** A local vector expressed in the rolled (bend) frame. */
  const intoBend = (x, y) => [cosR * x + sinR * y, -sinR * x + cosR * y];
  /** A vector in the bend frame expressed back in local coordinates. */
  const outOfBend = (u, v) => [cosR * u - sinR * v, sinR * u + cosR * v];

  /** Bend-frame coordinates: distance from the bend axis, angle, height. */
  function bendCoords(x, y, z) {
    const [bx, by] = intoBend(x, y);
    const dx = bx - centre[0];
    const rho = Math.hypot(dx, z);
    const theta = Math.atan2(z, dx);
    return { rho, theta, h: by, dx };
  }

  return {
    type: 'bender',
    label: 'Bender',
    params: p,
    length,
    bore: gap / 2,
    // The element's own TRANSVERSE half-extent - how far the metal reaches
    // across the orbit - not its bend radius. Using R + gap here would claim
    // the element was as wide as the corner it turns, which inflates the
    // view's margins and makes every aperture comparison meaningless.
    outerRadius: Math.max(gap, mmToM(p.height) / 2),
    lengthScale: step,
    shortestPeriod: null,
    warnings,
    grid,
    field: unit,
    curved: true,
    bendAngle: alpha,
    bendRadius: R,

    /**
     * The exit placement: round the arc, then facing along the new direction.
     *
     * This single value is what turns the rest of the column. Everything
     * downstream is composed onto it, so a bender needs no cooperation from
     * any other element to redirect the whole beamline.
     */
    /**
     * Conjugated by the roll: R then the bend then R inverse.
     *
     * The conjugation matters. Composing the roll and the bend without
     * undoing the roll would leave the downstream beam rotated about its own
     * axis, so a vertical bender would also turn "up" into "sideways" for
     * every element after it. What a bender should change is where the beam
     * goes, not which way is up.
     */
    exitTransform: compose(
      compose(
        rollFrame(roll),
        compose(
          translation(R * (Math.cos(alpha) - 1), 0, R * Math.sin(alpha)),
          yawFrame(alpha)
        )
      ),
      inverse(rollFrame(roll))
    ),

    /** A point a fraction of the way along the curved reference orbit. */
    pathPoint(f) {
      const th = f * alpha;
      const [x, y] = outOfBend(R * (Math.cos(th) - 1), 0);
      return [x, y, R * Math.sin(th)];
    },

    setVoltage(v) {
      p.voltage = v;
      // Nothing to recompute: the stored map is per volt across the pair and
      // `fieldAt` scales it, so this only records the setting.
    },

    /** Inside the angular sector, and within the solved radial band. */
    contains(x, y, z) {
      const { rho, theta } = bendCoords(x, y, z);
      if (theta < -1e-12 || theta > alpha + 1e-12) return false;
      return rho >= rMin - pad && rho <= rMax + pad;
    },

    fieldAt(x, y, z) {
      const { rho, theta, h, dx } = bendCoords(x, y, z);
      if (theta < 0 || theta > alpha) return { Ex: 0, Ey: 0, Ez: 0 };

      // The solved map's "z" axis is height across the plates and its "r"
      // axis is distance from the bend axis.
      const { Ez: Eh, Er: Erho } = unit.fieldAt(h, rho);
      const scale = p.voltage;
      if (rho === 0) {
        const [ex, ey] = outOfBend(0, scale * Eh);
        return { Ex: ex, Ey: ey, Ez: 0 };
      }

      // Resolve the radial field in the bend frame, then roll it back into
      // the element's own axes.
      const [ex, ey] = outOfBend((scale * Erho * dx) / rho, scale * Eh);
      return { Ex: ex, Ey: ey, Ez: (scale * Erho * z) / rho };
    },

    potentialAt(x, y, z) {
      const { rho, theta, h } = bendCoords(x, y, z);
      if (theta < 0 || theta > alpha) return 0;
      return p.voltage * unit.potentialAt(h, rho);
    },

    strikes(x, y, z) {
      const { rho, theta, h } = bendCoords(x, y, z);
      if (theta < 0 || theta > alpha) return false;
      if (Math.abs(h) > (nz - 1) * step * 0.5) return true;
      return rho <= rInner || rho >= rOuter;
    },

    /** Plate voltage that puts this ion on the central orbit. */
    matchedVoltage(energyEV, chargeStates) {
      return matchedVoltage(p, energyEV, chargeStates);
    },

    /**
     * Drawn as two curved plates, each a polygon following the arc.
     *
     * Supplied directly rather than as mirrored boxes because a bender is
     * neither axisymmetric nor straight: its plates sit at two different
     * radii on the same side of the orbit, and both follow a curve. Points
     * are in the element's own (transverse, axial) plane, and the beamline
     * puts them into global space.
     */
    shapes() {
      const steps = Math.max(8, Math.round((p.bendAngle / 90) * 18));
      const thick = Math.max(step * 2, gap * 0.12);

      // A point on the arc at angle theta and radius rho, in local (x, z).
      // Points come back in the element's own transverse plane, so a rolled
      // bender draws in the plane it actually bends in.
      const at = (theta, rho) => {
        const [x, y] = outOfBend(centre[0] + rho * Math.cos(theta), 0);
        return [x, y, rho * Math.sin(theta)];
      };

      const plate = (rho) => {
        const pts = [];
        for (let k = 0; k <= steps; k++) pts.push(at((k / steps) * alpha, rho));
        for (let k = steps; k >= 0; k--) {
          pts.push(at((k / steps) * alpha, rho + thick));
        }
        return { points: pts };
      };

      // Inner plate grows inward, outer plate outward, so neither intrudes on
      // the gap the beam travels through.
      return [plate(rInner - thick), plate(rOuter)];
    },
  };
}
