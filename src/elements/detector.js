/**
 * Ion detector - a biased collector that pulls ions in and counts them.
 *
 * The front end of a channeltron, a conversion dynode or a Faraday cup with a
 * bias supply: a grounded entrance aperture with an active surface a few
 * millimetres behind it, held several kilovolts below ground. The field
 * through that aperture is what gives a detector its collection efficiency -
 * ions that would otherwise drift past are pulled in and accelerated onto the
 * surface, which for a secondary-electron detector also matters because the
 * yield depends on impact energy.
 *
 * So this is not a target painted on the end of the beamline. It is an
 * electrode like any other, it bends trajectories before they reach it, and at
 * -3 kV it bends them a long way.
 *
 * What makes it a detector rather than a wall
 * -------------------------------------------
 * Every element already stops an ion that hits metal. What this adds is
 * knowing WHICH metal: `detects` is true only on the active surface, so an ion
 * that lands on the housing is a loss and one that lands on the face is a
 * count. Both are strikes and neither is transmission, and reporting them as
 * one number would hide the difference that the whole element exists to draw.
 *
 * Its reach depends on the field model
 * ------------------------------------
 * Solved alone, the detector's field stops at its own boundary, so it only
 * collects ions that were already heading into its aperture. With fringe
 * fields on (docs/PHYSICS.md section 10) it is solved together with whatever
 * is in front of it, and the -3 kV reaches back up the column - which is what
 * "sucks the ions in" actually means and is worth seeing both ways.
 */

import { PotentialArray, CYLINDRICAL } from '../grid.js';
import { solveBasis } from '../laplace.js';
import { Field } from '../field.js';
import { mmToM } from '../constants.js';

export const DETECTOR_DEFAULTS = {
  voltage: -3000, // V on the active surface
  entranceRadius: 6, // mm, the grounded aperture ions come through
  activeRadius: 8, // mm, the collecting surface behind it
  depth: 10, // mm, aperture to surface
  margin: 12, // mm of drift in front, so the fringe is inside the domain
  housingRadius: 14, // mm
  gridStep: 0.3, // mm
};

export function createDetector(params = {}, solverOpts = {}) {
  const p = { ...DETECTOR_DEFAULTS, ...params };
  const warnings = [];

  const entrance = mmToM(p.entranceRadius);
  const active = mmToM(p.activeRadius);
  const housing = mmToM(p.housingRadius);
  const depth = mmToM(p.depth);

  if (p.entranceRadius >= p.housingRadius) {
    throw new Error('The entrance aperture must be smaller than the housing');
  }
  if (p.activeRadius > p.housingRadius) {
    throw new Error('The active surface must fit inside the housing');
  }
  if (p.depth <= 0) throw new Error('The detector needs depth behind its aperture');

  const step = mmToM(p.gridStep);
  const front = mmToM(p.margin);
  const totalLength = front + depth;
  const nz = Math.round(totalLength / step) + 1;
  const nr = Math.round(housing / step) + 1;
  const grid = new PotentialArray({ nz, nr, step, symmetry: CYLINDRICAL });

  /*
    Ions come in and do not come out. The entrance face is an aperture, as
    everywhere else; the BACK face is the detector, so an ion reaching it has
    been collected and must be recorded as a strike rather than as having left
    the modelled region.
  */
  grid.openFaces = { zMin: true, zMax: false };

  const encl = grid.addElectrode('housing');
  const face = grid.addElectrode('face');
  grid.paintEnclosure(encl);

  // The grounded entrance plate: the front face, everything outside the
  // aperture. paintEnclosure already covered it; this is the aperture being
  // cut back out of it, which openFaces handles for the ion and which the
  // solve keeps as a Dirichlet boundary.
  const zFace = front + depth;
  const painted = grid.paint(
    face,
    (z, r) => z >= zFace - step * 1.5 && r <= active + step * 1e-6
  );
  if (painted === 0) throw new Error('Detector face covered no grid nodes');

  const { basis, reports } = solveBasis(grid, solverOpts);
  if (reports.some((r) => !r.converged)) {
    warnings.push('Laplace solve did not converge; this field is not trustworthy.');
  }

  const field = new Field(grid, basis);
  const setAll = () => field.setVoltages({ housing: 0, face: p.voltage });
  setAll();

  const length = grid.zLength;

  return {
    type: 'detector',
    label: 'Ion detector',
    params: p,
    length,
    bore: entrance,
    /*
      Deliberately no `clearBore`.

      That promise is "no metal at this radius anywhere along my length", and a
      detector is the one element that cannot make it: it is CLOSED at the
      back, and the closed part is the collecting surface - dead centre, where
      the clear bore would be widest. Declaring one short-circuits the strike
      test on the axis, so an ion flying straight into the middle of the
      collector was reported as being in free space and sailed out the back.
    */
    outerRadius: housing,
    lengthScale: grid.step,
    shortestPeriod: null,
    warnings,
    grid,
    field,

    setVoltage(v) {
      p.voltage = v;
      setAll(); // fast adjust; no relaxation here
    },

    /*
      A tolerance past the back, for the same reason the deflector has one at
      its entrance: an ion arriving at the collecting surface lands ON the
      element's last plane, and one step later it is outside a range that stops
      exactly there. With no element claiming it, it was recorded as having
      left the column - so the detector caught nothing and reported everything
      as transmitted, which is the most misleading answer available.
    */
    contains(x, y, zl) {
      return zl >= 0 && zl <= length + grid.step * 2;
    },

    fieldAt(x, y, zl) {
      return field.fieldAt3D(x, y, zl);
    },
    potentialAt(x, y, zl) {
      return field.potentialAt3D(x, y, zl);
    },
    strikes(x, y, zl) {
      // The back is solid: surface in the middle, housing around it. Anything
      // that has reached it has hit something.
      if (zl >= length - grid.step * 0.5) return true;
      return field.strikes(x, y, zl);
    },

    /**
     * Did this land on the active surface rather than on the housing?
     *
     * The distinction the element exists to draw: one is a count, the other is
     * a loss, and both are strikes.
     */
    detects(x, y, zl) {
      return zl >= zFace - grid.step * 2 && Math.hypot(x, y) <= active + grid.step;
    },

    // The collecting surface, then the grounded plate its aperture is cut in.
    rects: [
      { z0: zFace - mmToM(1), z1: zFace, r0: 0, r1: active },
      { z0: front - mmToM(1), z1: front, r0: entrance, r1: housing },
    ],
  };
}
