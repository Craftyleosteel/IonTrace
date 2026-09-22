/**
 * Einzel lens - three coaxial cylinders, the standard electrostatic focusing
 * element and the geometry SIMION's introductory demo is built around.
 *
 *            entrance            centre             exit
 *          |===========|      |==========|      |===========|
 *   r      |           |      |          |      |           |
 *   ^      -------------      ------------      -------------
 *   |   ion ........................................................> z
 *   |      -------------      ------------      -------------
 *   +--> z |===========|      |==========|      |===========|
 *             V = 0             V = Vc             V = 0
 *
 * The outer two cylinders sit at the same potential (normally ground) and the
 * centre one is biased. An ion therefore leaves the lens with the kinetic
 * energy it arrived with - the lens does no net work on it, which is exactly
 * what makes it an *einzel* ("single") lens and what distinguishes it from an
 * accelerating immersion lens. That property is a free physics check, and
 * tests/physics.test.js asserts it.
 *
 * Focusing is a second-order effect: the radial field kicks inward in one
 * half of each gap and outward in the other, but the ion is slower where the
 * inward kick acts (for a decelerating centre electrode), so the inward
 * impulse wins. The net result is convergent for either sign of Vc.
 *
 * All dimensions in, and only in, millimetres. Conversion to SI happens once,
 * at the boundary, via mmToM.
 */

import { PotentialArray, CYLINDRICAL } from '../grid.js';
import { solveBasis } from '../laplace.js';
import { Field } from '../field.js';
import { mmToM } from '../constants.js';

/** Default geometry, loosely following SIMION's introductory einzel demo. */
export const EINZEL_DEFAULTS = {
  boreRadius: 6, // mm, inner radius of every cylinder
  wallThickness: 2, // mm, radial thickness of the cylinder walls
  outerLength: 20, // mm, axial length of the two grounded cylinders
  centreLength: 20, // mm, axial length of the biased cylinder
  gap: 4, // mm, axial gap between adjacent cylinders
  entryDrift: 15, // mm, field-free run-up before the first cylinder
  exitDrift: 40, // mm, drift after the lens, where the focus forms
  housingRadius: 16, // mm, radius of the grounded enclosure
  gridStep: 0.5, // mm, node spacing
};

/**
 * Build the potential array, solve its basis solutions and wrap them in a
 * Field.
 *
 * The grounded enclosure is a real modelling choice, not a formality. An
 * unbounded Dirichlet problem has no unique solution on a finite grid, so the
 * domain is closed by a conducting box held at its own electrode potential
 * ("housing"). Placing that box too close to the bore distorts the field, so
 * housingRadius should stay comfortably outside the cylinders; the returned
 * `warnings` array says so if it does not.
 *
 * @param {Partial<typeof EINZEL_DEFAULTS>} [options]
 * @param {object} [solverOpts] Passed through to the Laplace relaxation.
 */
export function buildEinzelLens(options = {}, solverOpts = {}) {
  const g = { ...EINZEL_DEFAULTS, ...options };
  const warnings = [];

  const outerRadius = g.boreRadius + g.wallThickness;
  if (outerRadius >= g.housingRadius) {
    throw new Error(
      `housingRadius (${g.housingRadius} mm) must exceed the cylinder outer ` +
        `radius (${outerRadius} mm)`
    );
  }
  if (g.housingRadius - outerRadius < g.boreRadius) {
    warnings.push(
      'Enclosure sits within one bore radius of the cylinders; the field near ' +
        'the bore may be measurably distorted by the boundary.'
    );
  }
  if (g.gridStep > g.wallThickness) {
    warnings.push(
      `Grid step (${g.gridStep} mm) is coarser than the cylinder wall ` +
        `(${g.wallThickness} mm); electrodes may be under-resolved.`
    );
  }
  if (g.gridStep > g.gap / 2) {
    warnings.push(
      `Grid step (${g.gridStep} mm) resolves the inter-electrode gap ` +
        `(${g.gap} mm) with fewer than two cells; fringe fields, which are ` +
        'what actually focus the beam, will be poorly represented.'
    );
  }

  // Axial layout, in mm from the domain's left face.
  const z1 = g.entryDrift;
  const z2 = z1 + g.outerLength;
  const z3 = z2 + g.gap;
  const z4 = z3 + g.centreLength;
  const z5 = z4 + g.gap;
  const z6 = z5 + g.outerLength;
  const totalLength = z6 + g.exitDrift;

  // Node counts. The +1 turns a cell count into a node count; rounding keeps
  // the requested extents to within half a grid step.
  const nz = Math.round(totalLength / g.gridStep) + 1;
  const nr = Math.round(g.housingRadius / g.gridStep) + 1;

  const grid = new PotentialArray({
    nz,
    nr,
    step: mmToM(g.gridStep),
    symmetry: CYLINDRICAL,
    z0: 0,
  });

  const housing = grid.addElectrode('housing');
  const entrance = grid.addElectrode('entrance');
  const centre = grid.addElectrode('centre');
  const exit = grid.addElectrode('exit');

  grid.paintEnclosure(housing);

  // Painted after the enclosure so a cylinder touching the boundary wins the
  // node, which is the physically right precedence: it is the nearer metal.
  const cylinder = (zStart, zEnd) => (z, r) =>
    z >= mmToM(zStart) &&
    z <= mmToM(zEnd) &&
    r >= mmToM(g.boreRadius) &&
    r <= mmToM(outerRadius);

  const painted = {
    entrance: grid.paint(entrance, cylinder(z1, z2)),
    centre: grid.paint(centre, cylinder(z3, z4)),
    exit: grid.paint(exit, cylinder(z5, z6)),
  };

  for (const [name, count] of Object.entries(painted)) {
    if (count === 0) {
      throw new Error(`Electrode "${name}" covered no grid nodes; check geometry`);
    }
  }

  const { basis, reports } = solveBasis(grid, solverOpts);
  const field = new Field(grid, basis);

  return {
    grid,
    field,
    geometry: { ...g, outerRadius, totalLength, bounds: { z1, z2, z3, z4, z5, z6 } },
    solverReports: reports,
    warnings,
  };
}
