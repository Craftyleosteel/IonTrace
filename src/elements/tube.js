/**
 * A tube at a DC voltage - the single-element lens.
 *
 * One cylinder, biased, sitting inside a grounded housing with open vacuum at
 * each end. The simplest thing you can build that is not a drift and not an
 * aperture, and a real device: it is what you get when you take an einzel lens
 * and remove its two grounded guard cylinders.
 *
 * What it does to the beam, and what it does not
 * ----------------------------------------------
 * Inside the tube the ion is at a different kinetic energy - that is the whole
 * point of biasing it. Across the whole element it is not. The domain is closed
 * by grounded faces at both ends, so the ion starts and finishes at the same
 * potential, and whatever energy the entrance fringe took from it the exit
 * fringe gives back.
 *
 * So this focuses, and it does not accelerate. It is a LENS, and the tuner can
 * treat its voltage like any other lens voltage.
 *
 * Why that is not a limitation of the model but a statement about the hardware:
 * a biased tube only changes a beam's energy if the potential where the ion
 * ends differs from where it started, which means the thing after it must be at
 * a different potential too. A tube in a grounded beamline cannot do it. To
 * build an accelerating stage - a lift, or an elevator - you need this tube AND
 * a downstream element at another potential, solved together with fringe fields
 * on so the field between them is real rather than two clipped halves.
 *
 * Against the einzel
 * ------------------
 * An einzel's grounded guard cylinders shield its live electrode almost
 * completely: measured on the shipped geometry, a part in seven thousand
 * escapes past the exit cylinder. This has no such guards, so its fringe field
 * reaches out to whatever the housing allows - a decay length of
 * housingRadius/2.405, several millimetres rather than a fraction of one.
 *
 * That makes it the better element for showing what a fringe field IS, and the
 * worse one for putting next to something sensitive.
 */

import { PotentialArray, CYLINDRICAL } from '../grid.js';
import { solveBasis } from '../laplace.js';
import { Field } from '../field.js';
import { mmToM } from '../constants.js';
import { decayLength } from '../column.js';

export const TUBE_DEFAULTS = {
  bore: 6, // mm, clear radius inside the tube
  tubeLength: 30, // mm of biased cylinder
  wallThickness: 2, // mm, bore to the outside of the tube
  voltage: -500, // V on the cylinder
  margin: 18, // mm of vacuum each side, inside the grounded housing
  housingRadius: 16, // mm
  gridStep: 0.5, // mm
};

export function createTube(params = {}, solverOpts = {}) {
  const p = { ...TUBE_DEFAULTS, ...params };
  const warnings = [];

  if (p.bore <= 0) throw new Error('A tube needs a bore');
  if (p.tubeLength <= 0) throw new Error('A tube needs some length');
  if (p.bore + p.wallThickness > p.housingRadius) {
    throw new Error('The tube must fit inside the housing');
  }

  /*
    The margin is what keeps this element's own end caps out of its fringe
    field. Inside a grounded pipe of radius R the fringe dies as
    exp(-2.405 z/R), so the requirement is a number of decay lengths rather
    than a number of millimetres, and a wide housing needs a longer margin than
    a narrow one for the same honesty.
  */
  const decay = decayLength(p.housingRadius);
  const lengths = p.margin / decay;
  if (lengths < 2.5) {
    warnings.push(
      `The ${p.margin} mm margin is only ${lengths.toFixed(1)} decay lengths for a ` +
        `${p.housingRadius} mm housing, so this element's own end faces are clipping ` +
        `about ${(Math.exp(-lengths) * 100).toFixed(0)} % of its fringe field. Lengthen ` +
        'the margin, narrow the housing, or turn on fringe fields.'
    );
  }

  const step = mmToM(p.gridStep);
  const bore = mmToM(p.bore);
  const outer = mmToM(p.bore + p.wallThickness);
  const housing = mmToM(p.housingRadius);

  const totalMm = 2 * p.margin + p.tubeLength;
  const nz = Math.round(totalMm / p.gridStep) + 1;
  const nr = Math.round(p.housingRadius / p.gridStep) + 1;
  const grid = new PotentialArray({ nz, nr, step, symmetry: CYLINDRICAL });

  const housingId = grid.addElectrode('housing');
  const tubeId = grid.addElectrode('tube');
  grid.paintEnclosure(housingId);

  const zA = mmToM(p.margin);
  const zB = mmToM(p.margin + p.tubeLength);
  const painted = grid.paint(
    tubeId,
    (z, r) =>
      grid.spans(z, zA, zB) && r >= bore - step * 1e-6 && r <= outer + step * 1e-6
  );
  if (painted === 0) throw new Error('The tube covered no grid nodes; check the geometry');

  const { basis, reports } = solveBasis(grid, solverOpts);
  if (reports.some((r) => !r.converged)) {
    warnings.push('Laplace solve did not converge; this field is not trustworthy.');
  }

  const field = new Field(grid, basis);
  field.setVoltages({ housing: 0, tube: p.voltage });

  const length = grid.zLength;

  return {
    type: 'tube',
    label: 'Biased tube',
    params: p,
    length,
    bore,
    // The tube starts at the bore and the housing is further out still, so the
    // axis is clear from end to end. See `isFree` in beamline.js.
    clearBore: bore,
    outerRadius: housing,
    lengthScale: grid.step,
    shortestPeriod: null,
    warnings,
    grid,
    field,

    setVoltage(v) {
      p.voltage = v;
      field.setVoltages({ housing: 0, tube: v }); // fast adjust; no re-solve
    },

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
     * The housing pipe and the cylinder - not the grounded caps, which are the
     * numerical device that stops this element's fringe at its own boundary.
     * This element has more to gain from a column solve than any other here,
     * because it has no guard cylinders of its own to hide behind.
     */
    parts: [
      {
        name: 'housing',
        z0: 0,
        z1: length,
        r0: housing,
        r1: Infinity,
        voltage: () => 0,
      },
      { name: 'tube', z0: zA, z1: zB, r0: bore, r1: outer, voltage: () => p.voltage },
    ],

    rects: [{ z0: zA, z1: zB, r0: bore, r1: outer }],
  };
}
