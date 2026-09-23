/**
 * Column solves - fringe fields, and why they need more than one element.
 *
 * Every element in IonTrace is normally solved on its own grid, closed by
 * grounded faces at each end. Those faces are a numerical device: Laplace's
 * equation needs a closed boundary, and the ends of an element are where the
 * domain has to stop. But they are also, physically, a grounded plate a few
 * millimetres from the hardware - and a grounded plate is exactly what stops a
 * fringe field. So the isolated model does not merely omit the fringe: it
 * quietly shields it.
 *
 * Solving neighbouring elements together removes those internal faces. The
 * field then flows out of one element and into whatever is actually next to
 * it, which is the thing being asked for, and it decays at a rate set by that
 * neighbour rather than by an artefact.
 *
 * Why this cannot be done by superposition
 * ----------------------------------------
 * The obvious cheap alternative - solve each element alone, then add the
 * fields where they overlap - cannot show shielding, and shielding is the
 * whole point. A grounded plate contributes nothing to a sum: its potential is
 * zero everywhere on it, so adding its solution adds zero. Yet putting one
 * next to a charged lens changes that lens's field completely, because it
 * changes the BOUNDARY of the lens's own problem. Charge rearranges on both.
 *
 * Superposition is exact over ELECTRODE VOLTAGES on a fixed set of conductors
 * - that is the fast adjust this whole simulator is built on (docs/PHYSICS.md
 * section 3.2). It is not valid over geometry. Adding a conductor is a
 * different problem, not a bigger sum.
 *
 * So: one grid, all the metal on it, one solve.
 *
 * What decays how fast
 * --------------------
 * Inside a grounded pipe of radius R, a disturbance at one end dies along the
 * axis as the lowest Bessel mode,
 *
 *     phi(z) ~ exp(-j01 z / R),     j01 = 2.40483,
 *
 * so the decay length is R/2.405 and depends on NOTHING else - not the bore of
 * the element that made the field, not its length, not its voltage. Measured
 * against this solver over housing radii of 10 to 20 mm, fitted decay lengths
 * come out within 3 % of R/2.405.
 *
 * Two consequences worth knowing before reaching for a shield:
 *
 *   1. A narrow grounded drift tube is already an excellent shield. Five
 *      millimetres of bore gives a decay length of 2.1 mm, so a fringe is
 *      down by e^-7 within 15 mm.
 *   2. A wide housing is a poor one. Sixteen millimetres gives 6.7 mm, and the
 *      same 15 mm only knocks the field down by e^-2.2.
 *
 * Which is why the answer to "how do I stop this field reaching my detector"
 * is usually an aperture plate at ground, not more distance.
 *
 * What is left out
 * ----------------
 * Only axisymmetric elements can share an r-z grid, so a run stops at a
 * quadrupole or a deflector. That is not much of a loss: both of those sit
 * inside grounded housings of their own, which genuinely do terminate the
 * field rather than merely appearing to. A run also stops at a misaligned
 * element, which is not axisymmetric about the column any more.
 */

import { PotentialArray, CYLINDRICAL } from './grid.js';
import { solveBasis } from './laplace.js';
import { Field } from './field.js';

/** Element types whose geometry is a body of revolution about the beam. */
export const AXISYMMETRIC = new Set(['drift', 'aperture', 'einzel']);

/** First zero of J0: sets how fast a field dies inside a grounded pipe. */
export const BESSEL_J01 = 2.404825557695773;

/** Decay length of a fringe field inside a grounded pipe of radius `r`. */
export function decayLength(r) {
  return r / BESSEL_J01;
}

/** True if this element can take part in a shared r-z solve. */
export function canShareGrid(e) {
  if (!AXISYMMETRIC.has(e.typeKey) || !Array.isArray(e.parts)) return false;
  // A misaligned element is no longer a body of revolution about the column's
  // axis, so it cannot be painted on an r-z grid with its neighbours.
  const a = e.align;
  return !a || (a.dx === 0 && a.dy === 0);
}

/**
 * Maximal stretches of consecutive elements that can share a grid.
 *
 * Returned as [from, to] index pairs, inclusive. Single elements are included:
 * a lone lens between two deflectors still benefits, because solving it with
 * its own end caps removed is not the same as solving it with them.
 */
export function axisymmetricRuns(elements) {
  const runs = [];
  let from = -1;
  for (let i = 0; i <= elements.length; i++) {
    const ok = i < elements.length && canShareGrid(elements[i]);
    if (ok && from === -1) from = i;
    if (!ok && from !== -1) {
      runs.push([from, i - 1]);
      from = -1;
    }
  }
  // A stretch of nothing but drift has no electrode that can ever hold a
  // voltage, so its solution is zero everywhere and solving it is pure cost.
  // The elements' own (empty) fields say the same thing for free.
  return runs.filter(([a, b]) =>
    elements.slice(a, b + 1).some((e) => e.typeKey !== 'drift')
  );
}

/**
 * Solve one run of elements on a single grid.
 *
 * @param {object[]} elements the whole element list
 * @param {[number, number]} run inclusive index range
 * @param {object} [solverOpts]
 * @returns {{field: Field, grid: PotentialArray, from: number, to: number,
 *            z0: number, z1: number, sync: () => void, warnings: string[]}}
 */
export function buildRunField(elements, [from, to], solverOpts = {}) {
  const warnings = [];
  const run = elements.slice(from, to + 1);

  // Where the run sits along the reference orbit. Elements are laid end to
  // end, so the run spans one contiguous interval of path length.
  const z0 = run[0].zStart;
  const z1 = run[run.length - 1].zEnd;

  // The finest resolution anyone in the run asked for. Coarsening an element
  // below its own setting would change its solved field, which is not what
  // turning fringe fields on is supposed to do.
  let step = Infinity;
  for (const e of run) {
    if (Number.isFinite(e.lengthScale)) step = Math.min(step, e.lengthScale);
  }
  if (!Number.isFinite(step)) step = (z1 - z0) / 200; // a run of pure drift

  // Wide enough for the widest piece of hardware in the run.
  const rMax = Math.max(...run.map((e) => e.outerRadius));

  // Rounded UP, not to nearest. A drift's wall sits exactly at its bore, and
  // its bore may be the widest thing in the run; rounding down would put the
  // grid's last node a fraction of a step inside that wall, so the wall would
  // paint no nodes at all and the tube would have no metal in it.
  const nz = Math.max(3, Math.ceil((z1 - z0) / step) + 1);
  const nr = Math.max(3, Math.ceil(rMax / step) + 1);

  const grid = new PotentialArray({ nz, nr, step, symmetry: CYLINDRICAL });

  /*
    The run's own enclosure. Its two end faces are the same numerical device
    as before - they have to be, the problem still needs closing - but they are
    now at the ends of the RUN rather than of every element, so nothing inside
    is shielded by an artefact. Marked open so an ion may fly through them.

    Whether they are far enough from the nearest live electrode to be harmless
    is checked below, and warned about if not.
  */
  const ground = grid.addElectrode('ground');
  grid.paintEnclosure(ground);

  const ids = [];
  for (let k = 0; k < run.length; k++) {
    const e = run[k];
    const offset = e.zStart - z0;
    for (const part of e.parts) {
      // Namespaced, so two lenses in a run keep their own centre electrodes.
      const name = `${from + k}:${part.name}`;
      const id = grid.addElectrode(name);
      ids.push({ name, id, part, element: e });

      const a = offset + part.z0;
      const b = offset + part.z1;
      const r1 = Number.isFinite(part.r1) ? part.r1 : Infinity;
      const painted = grid.paint(
        id,
        (z, r) =>
          grid.spans(z, a, b) &&
          r >= part.r0 - grid.step * 1e-6 &&
          (r1 === Infinity || r <= r1 + grid.step * 1e-6)
      );
      if (painted === 0) {
        warnings.push(
          `${e.label}: "${part.name}" covered no nodes in the column solve.`
        );
      }
    }
  }

  const { basis, reports } = solveBasis(grid, solverOpts);
  if (reports.some((r) => !r.converged)) {
    warnings.push('Column solve did not converge; this field is not trustworthy.');
  }

  const field = new Field(grid, basis);

  /** Re-apply every element's voltages. Fast adjust: no re-solve. */
  const sync = () => {
    const v = {};
    for (const name of grid.electrodeNames) v[name] = 0;
    for (const { name, part } of ids) v[name] = part.voltage();
    field.setVoltages(v);
  };
  sync();

  /*
    How much field is still alive where the run's domain has to stop. This is
    the one artefact the column solve cannot remove - it can only push it
    somewhere harmless - so it is measured rather than assumed.
  */
  const peak = Math.max(
    ...ids.map(({ part }) => Math.abs(part.voltage())),
    1e-12
  );
  const atEnds = Math.max(
    Math.abs(field.potentialAt3D(0, 0, grid.step * 0.5)),
    Math.abs(field.potentialAt3D(0, 0, (nz - 1.5) * grid.step))
  );
  if (atEnds / peak > 0.01) {
    warnings.push(
      `The column's own end faces still carry ${((atEnds / peak) * 100).toFixed(1)} % ` +
        'of the largest electrode voltage, so they are clipping a real fringe ' +
        'field. Add a drift at that end of the line.'
    );
  }

  return { field, grid, from, to, z0, z1, sync, warnings };
}
