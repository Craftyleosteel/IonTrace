/**
 * The assumption the field drawing rests on.
 *
 * With column solves on, the potential map, the equipotentials and the field
 * lines are all painted from the RUN's grid rather than each element's, so
 * that a line can cross from one element into the next instead of stopping at
 * a boundary the ions do not see.
 *
 * Painting it in the right PLACE depends on one thing being true: a run's grid
 * starts exactly where its first element starts. If that ever drifts, the map
 * slides along the column and the failure is silent - a field drawn in the
 * wrong place still looks like a field.
 */

import { describe, it, assert, assertClose } from './harness.js';

import { Beamline } from '../src/beamline.js';
import { createElement } from '../src/elements/index.js';
import { axisymmetricRuns } from '../src/column.js';

/** A lens stack with no deflector, so the whole thing is one run. */
function stack() {
  return new Beamline([
    createElement('drift', { length: 10, bore: 6 }),
    createElement('aperture', { voltage: -300 }),
    createElement('drift', { length: 8, bore: 6 }),
    createElement('einzel', {
      voltage: -300, boreRadius: 6, housingRadius: 15, entryDrift: 4, exitDrift: 4,
    }),
    createElement('drift', { length: 10, bore: 6 }),
  ]);
}

describe('Field drawing across a run', () => {
  it('starts a run grid where its first element starts', () => {
    const bl = stack();
    bl.setFringe(true);
    assert(bl.runs.length > 0, 'the stack should form at least one run');

    for (const run of bl.runs) {
      const first = bl.elements[run.indices[0]];
      assertClose(run.z0, first.zStart, 1e-12, 'the run begins at its first element');
      // The grid itself is laid out from zero, so grid coordinate 0 and the
      // first element's own origin are the same point. That is what lets the
      // drawing put the run in the first element's frame.
      assertClose(run.grid.z0, 0, 1e-12, 'and its grid is measured from there');
    }
  });

  it('spans the whole run, not one element', () => {
    const bl = stack();
    bl.setFringe(true);
    const run = bl.runs[0];
    const covered = run.z1 - run.z0;
    assertClose(run.grid.zLength, covered, run.grid.step, 'the grid covers the run');

    const longest = Math.max(...run.indices.map((i) => bl.elements[i].length));
    assert(
      run.grid.zLength > longest * 1.5,
      'and is much longer than any single element in it'
    );
  });

  it('falls back to one field per element when column solves are off', () => {
    // The lines SHOULD stop at each boundary then, because the field does.
    const bl = stack();
    bl.setFringe(false);
    assert(bl.runs.length === 0, 'no shared solves with fringe fields off');
    for (let i = 0; i < bl.elements.length; i++) {
      assert(bl.runFor(i) === null, `element ${i} draws its own field`);
    }
  });

  it('breaks the run at a deflector', () => {
    /*
      A deflector is not a body of revolution, so it cannot share an r-z grid.
      The field lines stopping either side of it is therefore correct rather
      than a limitation of the drawing, and this pins that the run really does
      stop there rather than quietly painting across it.
    */
    const bl = new Beamline([
      createElement('drift', { length: 10, bore: 6 }),
      createElement('einzel', { voltage: -300, boreRadius: 6, housingRadius: 15 }),
      createElement('bender', {}),
      createElement('drift', { length: 20, bore: 6 }),
      createElement('aperture', { voltage: -200 }),
    ]);
    const runs = axisymmetricRuns(bl.elements);
    for (const indices of runs) {
      assert(
        !indices.some((i) => bl.elements[i].typeKey === 'bender'),
        'no run may contain a deflector'
      );
    }
    assert(runs.length >= 1, 'the parts either side can still form runs');
  });
});
