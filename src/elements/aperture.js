/**
 * Aperture plate - a disc of metal with a hole in it, held at a potential.
 *
 * This is the simplest real ion-optical element there is, and it is not a
 * passive hole: a charged plate with an aperture is a LENS. The equipotential
 * surfaces bulge through the hole, so an ion crossing the plane feels a radial
 * kick. For a thin aperture between two regions of axial field the focal
 * length is the Davisson-Calbick result
 *
 *     f = 4V / (E_2 - E_1)
 *
 * with V the beam energy in volts and E the axial fields either side. IonTrace
 * does not use that formula - the field is solved - but it is the reason an
 * aperture cannot be treated as a piece of empty space with a boundary.
 *
 * Geometry is axisymmetric, so this reuses the cylindrical solver directly.
 */

import { PotentialArray, CYLINDRICAL } from '../grid.js';
import { solveBasis } from '../laplace.js';
import { Field } from '../field.js';
import { mmToM } from '../constants.js';

export const APERTURE_DEFAULTS = {
  bore: 4, // mm, hole radius
  thickness: 2, // mm, plate thickness along z
  voltage: -200, // V
  housingRadius: 14, // mm
  margin: 14, // mm of grounded drift each side, inside this element
  gridStep: 0.4, // mm
};

/**
 * @param {Partial<typeof APERTURE_DEFAULTS>} [params]
 * @param {object} [solverOpts]
 */
export function createAperture(params = {}, solverOpts = {}) {
  const p = { ...APERTURE_DEFAULTS, ...params };
  const warnings = [];

  if (p.bore >= p.housingRadius) {
    throw new Error('Aperture bore must be smaller than the housing radius');
  }
  if (p.margin < p.bore) {
    warnings.push(
      `Only ${p.margin} mm of grounded drift either side of a ${p.bore} mm ` +
        'bore. The plate’s fringe field is being clipped by this element’s own ' +
        'end faces; widen the margin for a faithful lens.'
    );
  }
  if (p.gridStep > p.thickness / 2) {
    warnings.push(
      `Grid step (${p.gridStep} mm) resolves the ${p.thickness} mm plate with ` +
        'fewer than two cells.'
    );
  }

  const totalMm = 2 * p.margin + p.thickness;
  const nz = Math.round(totalMm / p.gridStep) + 1;
  const nr = Math.round(p.housingRadius / p.gridStep) + 1;

  const grid = new PotentialArray({
    nz,
    nr,
    step: mmToM(p.gridStep),
    symmetry: CYLINDRICAL,
  });

  const housing = grid.addElectrode('housing');
  const plate = grid.addElectrode('plate');
  grid.paintEnclosure(housing);

  const zA = mmToM(p.margin);
  const zB = mmToM(p.margin + p.thickness);
  const painted = grid.paint(
    plate,
    (z, r) => z >= zA && z <= zB && r >= mmToM(p.bore)
  );
  if (painted === 0) throw new Error('Aperture plate covered no grid nodes');

  const { basis, reports } = solveBasis(grid, solverOpts);
  if (reports.some((r) => !r.converged)) {
    warnings.push('Laplace solve did not converge; this field is not trustworthy.');
  }

  const field = new Field(grid, basis);
  field.setVoltages({ housing: 0, plate: p.voltage });

  const length = mmToM(totalMm);
  const bore = mmToM(p.bore);

  return {
    type: 'aperture',
    label: 'Aperture plate',
    params: p,
    length,
    bore,
    outerRadius: mmToM(p.housingRadius),
    lengthScale: grid.step,
    shortestPeriod: null,
    warnings,
    grid,
    field,

    setVoltage(v) {
      p.voltage = v;
      // Fast adjust: no re-solve, one pass over the grid.
      field.setVoltages({ housing: 0, plate: v });
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

    rects: [
      { z0: zA, z1: zB, r0: bore, r1: mmToM(p.housingRadius) },
    ],
  };
}
