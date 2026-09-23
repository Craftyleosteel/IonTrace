/**
 * Drift - a field-free tube.
 *
 * Not decoration. Every element in IonTrace is solved in isolation with
 * grounded end faces, so its fringe field is contained inside its own
 * footprint. That is only a fair approximation if neighbouring elements are
 * far enough apart not to see each other, and a drift is how you buy that
 * separation. Butting two active elements together is a real modelling error,
 * and the beamline warns about it.
 *
 * A drift still has a bore, so an ion that has been given too much transverse
 * velocity is stopped by the tube wall rather than flying through it.
 */

import { mmToM } from '../constants.js';

export const DRIFT_DEFAULTS = {
  length: 20, // mm
  bore: 8, // mm, tube inner radius
};

export function createDrift(params = {}) {
  const p = { ...DRIFT_DEFAULTS, ...params };
  const length = mmToM(p.length);
  const bore = mmToM(p.bore);

  return {
    type: 'drift',
    label: 'Drift',
    params: p,
    length,
    bore,
    outerRadius: bore,
    // No field structure to resolve, so a drift never constrains the time
    // step. Something in the line has to, and something always does.
    lengthScale: Infinity,
    shortestPeriod: null,
    warnings: [],

    /**
     * Whether a local point falls inside this element's span.
     *
     * Deliberately axial only: an ion that is within the element's length but
     * outside its bore is still THIS element's business, and `strikes` is
     * what decides it has hit the wall. Excluding it here would leave the
     * beamline unable to find any element for it, and a lost ion would be
     * reported as having wandered out of the column rather than as a strike.
     */
    contains(x, y, zl) {
      return zl >= 0 && zl <= length;
    },

    fieldAt() {
      return { Ex: 0, Ey: 0, Ez: 0 };
    },

    potentialAt() {
      return 0;
    },

    strikes(x, y) {
      return Math.hypot(x, y) > bore;
    },

    // Drawn as the tube wall only.
    rects: [{ z0: 0, z1: length, r0: bore, r1: bore * 1.06, wall: true }],
  };
}
