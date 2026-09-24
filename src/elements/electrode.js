/**
 * A bare conductor - whatever shape of revolution you ask for, at a DC voltage.
 *
 * Every other element here is a named device with a geometry chosen for it: an
 * einzel is three cylinders, an aperture is a plate with a hole. This is the
 * one that is not. It is a single rectangle in (z, r), spun about the axis,
 * held at a voltage, inside a grounded housing - and that one primitive covers
 * most of electrostatic optics:
 *
 *     r0 = bore, r1 = housing, thin in z   an aperture plate
 *     r0 = bore, r1 = bore + wall, long    a tube
 *     r0 = 0, thin in z                    a solid disc: a stop, or a cup face
 *     a narrow band of r                   a ring, as a funnel is built from
 *     r0 large, r1 larger                  a guard ring
 *
 * Why one conductor and not a list of them
 * ----------------------------------------
 * Because the beamline already composes things, and composing them there is
 * better than composing them inside an element. Put three of these in a row
 * with fringe fields on and they are solved together on one grid, with the
 * field flowing between them - which is the same answer a multi-electrode
 * element would give, except that each piece can be selected, tuned, dragged,
 * saved and explained on its own.
 *
 * It also means this element needs nothing special from the rest of the
 * program. One conductor is one voltage, so it is an ordinary tunable knob, an
 * ordinary row in the inspector and an ordinary entry in a saved file. An
 * element carrying an array of electrodes would have to teach the optimiser,
 * the inspector and the scene format about arrays first, and would still not
 * let you tune one ring without tuning its neighbours.
 *
 * What it is NOT
 * --------------
 * Axisymmetric. This is a body of revolution, so it can be a ring but never a
 * pair of rods, a slit, or anything with corners in the transverse plane. The
 * solver would happily accept such a thing on a transverse grid - that is how
 * the quadrupole and the deflector work - but then it could not share an r-z
 * column solve with its neighbours, which is most of the point of being able
 * to stack these.
 */

import { PotentialArray, CYLINDRICAL } from '../grid.js';
import { solveBasis } from '../laplace.js';
import { Field } from '../field.js';
import { mmToM } from '../constants.js';
import { decayLength } from '../column.js';

export const ELECTRODE_DEFAULTS = {
  /*
    Long enough not to trip its own margin warning. The housing is 16 mm, so
    the fringe decay length is 16/2.405 = 6.65 mm and an honest margin is about
    2.5 of those - 17 mm either side of the metal. An element that warns the
    moment it is placed has taught the user to ignore its warnings.
  */
  length: 40, // mm, the whole element including vacuum either side
  z0: 18, // mm, where the metal starts along the element
  z1: 22, // mm, where it ends
  r0: 4, // mm, inner radius - 0 for a solid disc
  r1: 14, // mm, outer radius
  voltage: -200, // V
  housingRadius: 16, // mm, the grounded pipe around it
  gridStep: 0.4, // mm
};

export function createElectrode(params = {}, solverOpts = {}) {
  const p = { ...ELECTRODE_DEFAULTS, ...params };
  const warnings = [];

  if (!(p.z1 > p.z0)) throw new Error('The conductor needs z1 greater than z0');
  if (!(p.r1 > p.r0)) throw new Error('The conductor needs r1 greater than r0');
  if (p.z0 < 0 || p.z1 > p.length) {
    throw new Error('The conductor must sit inside the element length');
  }
  if (p.r0 < 0 || p.r1 > p.housingRadius) {
    throw new Error('The conductor must fit inside the housing');
  }

  /*
    Vacuum between the metal and this element's own end caps. Those caps are
    grounded, so too little margin means the element is shielding its own
    fringe field - and how much is "too little" is set by the housing, because
    inside a grounded pipe of radius R a disturbance dies as exp(-2.405 z/R).
  */
  const decay = decayLength(p.housingRadius);
  const gap = Math.min(p.z0, p.length - p.z1);
  const lengths = gap / decay;
  if (lengths < 2.5) {
    warnings.push(
      `Only ${gap.toFixed(1)} mm of vacuum between the conductor and this element's ` +
        `end face - ${lengths.toFixed(1)} decay lengths for a ${p.housingRadius} mm ` +
        `housing, so roughly ${(Math.exp(-lengths) * 100).toFixed(0)} % of its fringe ` +
        'field is being clipped. Lengthen the element or turn on fringe fields.'
    );
  }

  const touching = p.r1 >= p.housingRadius - 1e-9;
  if (touching) {
    warnings.push(
      'The conductor reaches the housing, so it interrupts the grounded wall rather ' +
        'than sitting inside it. That is what a plate mounted across a tube does, but ' +
        'if it was meant to be insulated from the wall, pull r1 in.'
    );
  }

  const step = mmToM(p.gridStep);
  const nz = Math.round(p.length / p.gridStep) + 1;
  const nr = Math.round(p.housingRadius / p.gridStep) + 1;
  const grid = new PotentialArray({ nz, nr, step, symmetry: CYLINDRICAL });

  const housingId = grid.addElectrode('housing');
  const metalId = grid.addElectrode('metal');

  // Enclosure first, conductor over it, so the conductor owns any node the two
  // share. Same order as the aperture plate, and for the same reason: where
  // they meet, it is the driven surface the ion sees.
  grid.paintEnclosure(housingId);

  const zA = mmToM(p.z0);
  const zB = mmToM(p.z1);
  const rA = mmToM(p.r0);
  const rB = mmToM(p.r1);
  const painted = grid.paint(
    metalId,
    (z, r) =>
      grid.spans(z, zA, zB) && r >= rA - step * 1e-6 && r <= rB + step * 1e-6
  );
  if (painted === 0) {
    throw new Error(
      'The conductor covered no grid nodes - it is thinner than the grid step'
    );
  }

  const { basis, reports } = solveBasis(grid, solverOpts);
  if (reports.some((r) => !r.converged)) {
    warnings.push('Laplace solve did not converge; this field is not trustworthy.');
  }

  const field = new Field(grid, basis);
  field.setVoltages({ housing: 0, metal: p.voltage });

  const length = grid.zLength;

  return {
    type: 'electrode',
    label: 'Conductor',
    params: p,
    length,
    bore: rA,
    /*
      A clear bore only when the metal actually leaves the axis clear. With
      r0 = 0 this conductor is a solid disc blocking the whole aperture, and
      promising a clear radius there would tell the strike test that the middle
      of a beam stop is free space - so an ion would fly straight through the
      one element whose entire purpose is to stop it.
    */
    ...(p.r0 > 0 ? { clearBore: rA } : {}),
    outerRadius: mmToM(p.housingRadius),
    lengthScale: grid.step,
    shortestPeriod: null,
    warnings,
    grid,
    field,

    setVoltage(v) {
      p.voltage = v;
      field.setVoltages({ housing: 0, metal: v }); // fast adjust; no re-solve
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
     * The housing pipe is left off when the conductor reaches it: the two
     * would then be the same surface at two different potentials, and in a
     * shared grid whichever was painted second would silently win along the
     * whole run rather than just inside this element.
     */
    parts: [
      ...(touching
        ? []
        : [
            {
              name: 'housing',
              z0: 0,
              z1: length,
              r0: mmToM(p.housingRadius),
              r1: Infinity,
              voltage: () => 0,
            },
          ]),
      { name: 'metal', z0: zA, z1: zB, r0: rA, r1: rB, voltage: () => p.voltage },
    ],

    rects: [{ z0: zA, z1: zB, r0: rA, r1: rB }],
  };
}
