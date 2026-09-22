/**
 * Ion construction - the single conversion point between the units a chemist
 * types and the SI units the integrator requires.
 *
 * Nothing downstream of makeIon() should ever see a dalton, an electronvolt,
 * a millimetre or a degree.
 */

import {
  amuToKg,
  chargesToCoulombs,
  eVToJoules,
  mmToM,
  speedFromKineticEnergy,
  relativisticError,
} from './constants.js';

/**
 * Build an ion state in SI units from practical units.
 *
 * @param {object} spec
 * @param {number} spec.mass      Mass in daltons (u).
 * @param {number} spec.charge    Charge state in elementary charges (signed).
 * @param {number} spec.energy    Initial kinetic energy in eV.
 * @param {number} [spec.x]       Initial transverse offset in mm (signed).
 * @param {number} [spec.z]       Initial axial position in mm.
 * @param {number} [spec.angle]   Launch angle from the z axis, in degrees.
 * @returns {import('./integrator.js').IonState}
 */
export function makeIon({ mass, charge, energy, x = 0, z = 0, angle = 0 }) {
  if (!(mass > 0)) throw new Error('Ion mass must be positive');
  if (charge === 0) throw new Error('A neutral particle feels no electric force');
  if (energy < 0) throw new Error('Kinetic energy must be non-negative');

  const massKg = amuToKg(mass);
  const speed = speedFromKineticEnergy(eVToJoules(energy), massKg);
  const theta = (angle * Math.PI) / 180;

  const err = relativisticError(speed);
  if (err > 1e-3) {
    console.warn(
      `IonTrace: ion at ${energy} eV is ${(err * 100).toFixed(2)}% off the ` +
        'relativistic kinetic energy; this build is non-relativistic.'
    );
  }

  return {
    mass: massKg,
    charge: chargesToCoulombs(charge),
    x: mmToM(x),
    z: mmToM(z),
    vx: speed * Math.sin(theta),
    vz: speed * Math.cos(theta),
    t: 0,
  };
}

/**
 * A fan of ions launched parallel to the axis at evenly spaced radii.
 *
 * This is the classic lens-characterisation beam: parallel rays in, and the
 * axial crossing point of the outgoing rays is the back focal point. Rays are
 * spaced from `-maxOffset` to `+maxOffset` so both halves of the lens are
 * exercised and any spurious up-down asymmetry shows up immediately.
 */
export function parallelBeam({ count = 7, maxOffset = 4, ...spec }) {
  const ions = [];
  if (count === 1) return [makeIon({ ...spec, x: 0, angle: 0 })];
  for (let n = 0; n < count; n++) {
    const frac = (2 * n) / (count - 1) - 1; // -1 .. +1
    ions.push(makeIon({ ...spec, x: frac * maxOffset, angle: 0 }));
  }
  return ions;
}

/**
 * A fan of ions launched from one point over a spread of angles - a point
 * source, used to find the image plane rather than the focal plane.
 */
export function divergentBeam({ count = 7, maxAngle = 3, ...spec }) {
  const ions = [];
  if (count === 1) return [makeIon({ ...spec, angle: 0 })];
  for (let n = 0; n < count; n++) {
    const frac = (2 * n) / (count - 1) - 1;
    ions.push(makeIon({ ...spec, angle: frac * maxAngle }));
  }
  return ions;
}

/**
 * Axial crossing of a trajectory: the z at which the path last changes the
 * sign of x, found by linear interpolation between the bracketing points.
 *
 * For a parallel input ray this is the lens's focal point. Returns null if the
 * ray never crosses, which is what a diverging (or untouched) ray does.
 */
export function axialCrossing(points) {
  for (let n = points.length - 1; n > 0; n--) {
    const a = points[n - 1];
    const b = points[n];
    if (a.x === 0) return a.z;
    if (a.x * b.x < 0) {
      const f = a.x / (a.x - b.x); // fraction of the way from a to b
      return a.z + f * (b.z - a.z);
    }
  }
  return null;
}
