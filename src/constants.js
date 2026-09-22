/**
 * Physical constants and unit conversions.
 *
 * IonTrace works in strict SI internally: metres, seconds, kilograms,
 * coulombs, volts. Every value that enters or leaves the solver is
 * converted at the boundary, never inside the physics loop. Unit slips are
 * the most common source of silently-wrong trajectories, so the conversion
 * helpers below are the only sanctioned way in and out.
 *
 * Values are CODATA 2022 / SI-2019 exact where applicable.
 */

/** Elementary charge, C. Exact by the 2019 SI redefinition. */
export const ELEMENTARY_CHARGE = 1.602176634e-19;

/** Unified atomic mass unit (dalton), kg. CODATA 2022. */
export const ATOMIC_MASS_UNIT = 1.66053906892e-27;

/** Electron mass, kg. CODATA 2022. */
export const ELECTRON_MASS = 9.1093837139e-31;

/** Vacuum permittivity, F/m. CODATA 2022. */
export const VACUUM_PERMITTIVITY = 8.8541878188e-12;

/** Speed of light in vacuum, m/s. Exact. */
export const SPEED_OF_LIGHT = 299792458;

/**
 * One electronvolt in joules. Numerically equal to the elementary charge
 * because 1 eV is defined as e x 1 V, but named separately so that energy
 * expressions read as energy rather than as charge.
 */
export const ELECTRONVOLT = ELEMENTARY_CHARGE;

/** Millimetres -> metres. Geometry is authored in mm; physics runs in m. */
export const mmToM = (mm) => mm * 1e-3;

/** Metres -> millimetres, for display only. */
export const mToMm = (m) => m * 1e3;

/** Mass in daltons (u) -> kg. */
export const amuToKg = (amu) => amu * ATOMIC_MASS_UNIT;

/** Charge in elementary charges -> coulombs. A charge state may be negative. */
export const chargesToCoulombs = (z) => z * ELEMENTARY_CHARGE;

/** Kinetic energy in eV -> joules. */
export const eVToJoules = (eV) => eV * ELECTRONVOLT;

/** Kinetic energy in joules -> eV, for display. */
export const joulesToEV = (J) => J / ELECTRONVOLT;

/**
 * Non-relativistic speed of a particle of mass `massKg` with kinetic energy
 * `energyJ`:  KE = 1/2 m v^2  =>  v = sqrt(2 KE / m).
 *
 * IonTrace is explicitly non-relativistic (see docs/PHYSICS.md). This helper
 * is the single place that assumption is expressed, so `relativisticError`
 * below can be used to check it is still safe for a given case.
 */
export function speedFromKineticEnergy(energyJ, massKg) {
  if (energyJ < 0) throw new Error('Kinetic energy must be non-negative');
  if (massKg <= 0) throw new Error('Mass must be positive');
  return Math.sqrt((2 * energyJ) / massKg);
}

/**
 * Fractional error incurred by treating a particle at speed `v` as Newtonian.
 *
 * The exact relativistic kinetic energy is (gamma - 1)mc^2 and the Newtonian
 * value 1/2 mv^2 underestimates it. The quantity wanted is
 *
 *     [ (gamma - 1)mc^2 - 1/2 mv^2 ] / (gamma - 1)mc^2
 *
 * which must not be evaluated as written. For an ion, beta is of order 1e-4,
 * so gamma - 1 is around 1e-8: forming it as gamma minus one throws away
 * eight of the sixteen available digits, and the subtraction in the numerator
 * then cancels almost everything that survives. The result is wrong in its
 * first significant figure.
 *
 * Writing s = sqrt(1 - beta^2), the same expression reduces exactly to
 *
 *     beta^2 (2 + s) / (2 (1 + s))
 *
 * with no subtraction of nearly-equal quantities anywhere. It tends to the
 * familiar (3/4) beta^2 as beta -> 0 and stays exact all the way to beta -> 1.
 *
 * For singly-charged keV ions this returns about 1e-8; the non-relativistic
 * assumption is safe until it approaches a part in a thousand.
 *
 * Note the result depends only on speed. Mass cancels, which is why it is not
 * a parameter: a proton and an electron at the same speed are equally (non-)
 * relativistic, even though reaching that speed costs them very different
 * energies.
 */
export function relativisticError(speedMs) {
  const beta = speedMs / SPEED_OF_LIGHT;
  if (beta >= 1) return Infinity;
  const s = Math.sqrt(1 - beta * beta);
  return (beta * beta * (2 + s)) / (2 * (1 + s));
}
