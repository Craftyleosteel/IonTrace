/**
 * Einzel lens as a beamline element.
 *
 * A thin wrapper over the validated geometry in ../geometries/einzel.js: the
 * physics, the tests and the documentation all still refer to that module, and
 * this only adapts it to the element interface the beamline expects.
 */

import { buildEinzelLens } from '../geometries/einzel.js';
import { mmToM } from '../constants.js';

export const EINZEL_ELEMENT_DEFAULTS = {
  boreRadius: 6, // mm
  wallThickness: 2, // mm
  outerLength: 20, // mm
  centreLength: 20, // mm
  gap: 4, // mm
  // About three bore radii. Shorter margins let this element's own grounded
  // end faces clip its fringe field - see the warning in geometries/einzel.js.
  entryDrift: 18, // mm, this element's own grounded margin
  exitDrift: 18, // mm
  housingRadius: 16, // mm
  gridStep: 0.5, // mm
  voltage: -2000, // V on the centre electrode
};

export function createEinzel(params = {}, solverOpts = {}) {
  const p = { ...EINZEL_ELEMENT_DEFAULTS, ...params };
  const built = buildEinzelLens(p, solverOpts);
  const { grid, field, geometry, warnings } = built;

  field.setVoltages({ housing: 0, entrance: 0, centre: p.voltage, exit: 0 });

  const b = geometry.bounds;
  const bore = mmToM(p.boreRadius);
  const outer = mmToM(geometry.outerRadius);
  // The grid's own extent, not the requested length. Node counts are rounded,
  // so a requested 85 mm at 0.4 mm resolution actually spans 85.2 mm, and the
  // beamline must place the next element where this one's solved domain ends.
  const length = grid.zLength;

  return {
    type: 'einzel',
    label: 'Einzel lens',
    params: p,
    length,
    bore,
    // The three cylinders all start at the bore, so inside it the lens is
    // clear along its whole length. See `isFree` in beamline.js.
    clearBore: bore,
    outerRadius: mmToM(p.housingRadius),
    lengthScale: grid.step,
    shortestPeriod: null,
    warnings,
    grid,
    field,
    geometry,

    setVoltage(v) {
      p.voltage = v;
      // Fast adjust; no relaxation happens here.
      field.setVoltages({ housing: 0, entrance: 0, centre: v, exit: 0 });
    },

    // Axial only: an ion inside the length but outside the bore is still
    // this element's business, and strikes decides it has hit metal.
    contains(x, y, zl) {
      return zl >= 0 && zl <= length;
    },

    fieldAt(x, y, zl) {
      return field.fieldAt3D(x, y, zl);
    },
    potentialAt(x, y, zl) {
      return field.potentialAt3D(x, y, zl);
    },
    strikes(x, y, zl) {
      return field.strikes(x, y, zl);
    },

    /**
     * The metal, for a solve that spans the column.
     *
     * The three cylinders and the housing pipe, and nothing else - in
     * particular not the grounded caps closing this element's own domain.
     *
     * An einzel turns out to need this least of any element here. Its outer
     * two cylinders are grounded, so they shield its centre electrode almost
     * completely: measured on the shipped geometry, 0.285 V escapes past the
     * exit cylinder from a lens sitting at 1946 V on axis, a part in seven
     * thousand. A lens is very nearly its own Faraday cage. The entry and
     * exit drifts in its own domain are therefore mostly wasted space, and
     * what the column solve buys here is not a fringe field but the ability
     * to put something else inside that space.
     */
    parts: [
      {
        name: 'housing',
        z0: 0,
        z1: length,
        r0: mmToM(p.housingRadius),
        r1: Infinity,
        voltage: () => 0,
      },
      { name: 'entrance', z0: mmToM(b.z1), z1: mmToM(b.z2), r0: bore, r1: outer, voltage: () => 0 },
      { name: 'centre', z0: mmToM(b.z3), z1: mmToM(b.z4), r0: bore, r1: outer, voltage: () => p.voltage },
      { name: 'exit', z0: mmToM(b.z5), z1: mmToM(b.z6), r0: bore, r1: outer, voltage: () => 0 },
    ],

    rects: [
      { z0: mmToM(b.z1), z1: mmToM(b.z2), r0: bore, r1: outer },
      { z0: mmToM(b.z3), z1: mmToM(b.z4), r0: bore, r1: outer },
      { z0: mmToM(b.z5), z1: mmToM(b.z6), r0: bore, r1: outer },
    ],
  };
}
