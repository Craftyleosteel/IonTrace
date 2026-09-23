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
 * The state is three-dimensional. `y` and `azimuth` default to zero, which
 * launches the ion in the x-z plane; for an axisymmetric element it then
 * stays there for ever, because the field has no azimuthal component to take
 * it out. That is why the 2D picture remains exact for a lens, and why adding
 * the third dimension changes nothing for the cases that never needed it.
 *
 * @param {object} spec
 * @param {number} spec.mass      Mass in daltons (u).
 * @param {number} spec.charge    Charge state in elementary charges (signed).
 * @param {number} spec.energy    Initial kinetic energy in eV.
 * @param {number} [spec.x]       Initial transverse offset in mm (signed).
 * @param {number} [spec.y]       Initial second transverse offset in mm.
 * @param {number} [spec.z]       Initial axial position in mm.
 * @param {number} [spec.angle]   Polar launch angle from the z axis, degrees.
 * @param {number} [spec.azimuth] Direction of the transverse velocity in the
 *                                x-y plane, degrees. Only meaningful when
 *                                `angle` is non-zero.
 * @returns {import('./integrator.js').IonState}
 */
export function makeIon({
  mass,
  charge,
  energy,
  x = 0,
  y = 0,
  z = 0,
  angle = 0,
  azimuth = 0,
}) {
  if (!(mass > 0)) throw new Error('Ion mass must be positive');
  if (charge === 0) throw new Error('A neutral particle feels no electric force');
  if (energy < 0) throw new Error('Kinetic energy must be non-negative');

  const massKg = amuToKg(mass);
  const speed = speedFromKineticEnergy(eVToJoules(energy), massKg);
  const theta = (angle * Math.PI) / 180;
  const phi = (azimuth * Math.PI) / 180;

  const err = relativisticError(speed);
  if (err > 1e-3) {
    console.warn(
      `IonTrace: ion at ${energy} eV is ${(err * 100).toFixed(2)}% off the ` +
        'relativistic kinetic energy; this build is non-relativistic.'
    );
  }

  const transverse = speed * Math.sin(theta);

  return {
    mass: massKg,
    charge: chargesToCoulombs(charge),
    x: mmToM(x),
    y: mmToM(y),
    z: mmToM(z),
    vx: transverse * Math.cos(phi),
    vy: transverse * Math.sin(phi),
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
 * A round beam: ions spread over a DISC in the transverse plane.
 *
 * This is what a real beam looks like, and it is what `parallelBeam` is not.
 * That one places every ion on the x axis with y = 0, which is a line of ions
 * rather than a beam. In any element here y = 0 is a symmetry plane - for an
 * axisymmetric lens because the radial field has no azimuthal component, and
 * for a quadrupole because E_y is proportional to y - so such a beam stays in
 * that plane for ever. The motion looks two-dimensional because the SOURCE is
 * two-dimensional, not because the integrator is.
 *
 * A meridional fan is still the right tool for characterising a lens, where
 * seeing every ray cross the axis in one plane is the point. Use a disc when
 * the question is what a beam actually does - and always for a quadrupole,
 * whose whole behaviour is that the two transverse planes differ.
 *
 * Ions are placed on a Fermat (sunflower) spiral: radius proportional to
 * sqrt(k) so the areal density is uniform, with the golden angle between
 * successive ions so no two line up. It is deterministic - no random number
 * generator - so the same beam is reproducible run to run, which matters when
 * a trajectory is being compared against another.
 *
 * @param {object} spec
 * @param {number} [spec.count]      Number of ions.
 * @param {number} [spec.radius]     Beam radius in mm.
 * @param {number} [spec.divergence] Half-angle in degrees at the beam edge.
 */
export function discBeam({ count = 9, radius = 1.5, divergence = 0, ...spec }) {
  if (count <= 1) return [makeIon({ ...spec, x: 0, y: 0 })];

  // The angle between successive seeds in a sunflower head.
  const GOLDEN_ANGLE = Math.PI * (3 - Math.sqrt(5));
  const ions = [];

  for (let k = 0; k < count; k++) {
    // +0.5 keeps the innermost ion off the exact centre, so a beam of any
    // count samples the disc rather than stacking one ion on the axis.
    const frac = Math.sqrt((k + 0.5) / count);
    const r = radius * frac;
    const theta = k * GOLDEN_ANGLE;

    ions.push(
      makeIon({
        ...spec,
        x: r * Math.cos(theta),
        y: r * Math.sin(theta),
        // Divergence grows linearly with radius and points radially outward,
        // so the beam expands from a waist here rather than from a point
        // source at some arbitrary distance upstream. The ion on the axis
        // travels straight; the one at the edge gets the full half-angle.
        angle: divergence * frac,
        azimuth: (theta * 180) / Math.PI,
      })
    );
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
 * First axial crossing of a trajectory while the ion is travelling forwards,
 * found by linear interpolation between the bracketing points.
 *
 * For a parallel input ray this is the lens's focal point, which is why it
 * must be the FIRST crossing and not the last: a ray that crosses, diverges
 * and is turned again by a later element focuses at the first crossing. It
 * must also require forward motion, because a reflected ion re-crosses the
 * axis on its way back out and that crossing is not a focus of anything.
 *
 * Returns null when there is no crossing: a diverging ray, a reflected ray,
 * or a ray that lies on the axis for its whole flight - the last of which has
 * no crossing at all rather than one everywhere.
 */
export function axialCrossing(points) {
  for (let n = 1; n < points.length; n++) {
    const a = points[n - 1];
    const b = points[n];

    // A returning ion's crossing is not a focus.
    if (b.vz <= 0) continue;

    // Identically on the axis over this interval: no crossing here.
    if (a.x === 0 && b.x === 0) continue;

    // Landing exactly on the axis is a crossing; a strict sign product would
    // miss it, since a.x * 0 is never negative.
    if (b.x === 0) return b.z;

    if (a.x * b.x < 0) {
      const f = a.x / (a.x - b.x); // fraction of the way from a to b
      return a.z + f * (b.z - a.z);
    }
  }
  return null;
}

/**
 * Where a ray crosses the axis, extrapolating beyond the modelled region if
 * it has not crossed by the time it leaves.
 *
 * Without this, the reported focus is truncated by the domain: a weak lens
 * whose focus lies past the end of the grid reports "no crossing", and one
 * whose focus lies just past the end reports only the most aberrated outer
 * rays, biasing the number low. Beyond the last electrode the field is
 * negligible and the ray is straight, so the crossing follows exactly from
 * the exit position and slope:
 *
 *     z_cross = z - x * (vz / vx)
 *
 * @returns {{z: number, extrapolated: boolean} | null}
 */
export function focalCrossing(points) {
  const inside = axialCrossing(points);
  if (inside !== null) return { z: inside, extrapolated: false };

  const p = points[points.length - 1];
  // A reflected, parallel, or diverging ray has no focus ahead of it.
  if (p.vz <= 0 || p.vx === 0) return null;
  if (p.x * p.vx >= 0) return null; // moving away from the axis

  const z = p.z - p.x * (p.vz / p.vx);
  if (!Number.isFinite(z) || z <= p.z) return null;
  return { z, extrapolated: true };
}
