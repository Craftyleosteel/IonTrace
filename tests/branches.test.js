/**
 * Tuning a branching column.
 *
 * The thing worth testing is not that the optimiser runs several times - it is
 * that the answers genuinely differ. A deflector is a switch, so if two
 * branches come back with the same voltages then either the tuning is not
 * targeting anything or the column does not actually branch, and both of those
 * would pass a test that only checked it produced some numbers.
 */

import { describe, it, assert } from './harness.js';

import { Beamline } from '../src/beamline.js';
import { createElement } from '../src/elements/index.js';
import { discBeam } from '../src/ion.js';
import {
  tunableKnobs,
  optimizeBranches,
  branchLabel,
  readKnobs,
  writeKnobs,
  scoreBeamline,
} from '../src/optimize.js';

const SPEC = { mass: 100, charge: 1, energy: 50 };
const beam = (count = 5) => () => discBeam({ ...SPEC, count, radius: 1.0 });

/** A deflector with a line on its bent port and another straight through. */
function switched() {
  const bl = new Beamline([
    createElement('drift', { length: 12, bore: 5 }),
    createElement('bender', { voltage: 0 }),
  ]);
  const bend = bl.elements[1];
  bl.add(createElement('drift', { length: 30, bore: 5 }), undefined, {
    parent: bend,
    port: 'bend',
  });
  bl.add(createElement('drift', { length: 30, bore: 5 }), undefined, {
    parent: bend,
    port: 'straight',
  });
  return bl;
}

const FAST = { passes: 1, coarse: 9, levels: 1, polish: false };

describe('Tuning each branch', () => {
  it('finds one setting per open end', async () => {
    const bl = switched();
    const ends = bl.openEnds();
    assert(ends.length >= 2, `a switched column has several ends, got ${ends.length}`);

    const knobs = tunableKnobs(bl, SPEC);
    const { branches } = await optimizeBranches(bl, beam(), knobs, FAST);
    assert(
      branches.length === ends.length,
      `expected one result per end, got ${branches.length} of ${ends.length}`
    );
    for (const b of branches) {
      assert(b.label && b.label.length > 0, 'every branch is named');
      assert(Array.isArray(b.settings) && b.settings.length === knobs.length,
        'and carries a full set of voltages');
    }
  });

  it('gives the bent and the straight branch different voltages', async () => {
    /*
      The whole point. A deflector at the voltage that turns the beam does not
      pass it straight through, so these two answers must disagree - and on the
      deflector in particular, not merely somewhere in the vector.
    */
    const bl = switched();
    const knobs = tunableKnobs(bl, SPEC);
    const iBend = knobs.findIndex((k) => bl.elements[k.index].typeKey === 'bender');
    assert(iBend >= 0, 'the deflector is tunable');

    const { branches } = await optimizeBranches(bl, beam(), knobs, FAST);
    const delivered = branches.filter((b) => b.transmitted > 0);
    assert(delivered.length >= 2, `at least two branches should be reachable, got ${delivered.length}`);

    const volts = delivered.map((b) => b.settings[iBend]);
    const spread = Math.max(...volts) - Math.min(...volts);
    assert(
      spread > 1,
      `branches should want different deflector voltages, got ${volts.map((v) => v.toFixed(1)).join(', ')}`
    );
  });

  it('leaves the column exactly as it found it', async () => {
    // Tuning several destinations and leaving the column on whichever was last
    // would silently change what the user was looking at.
    const bl = switched();
    const knobs = tunableKnobs(bl, SPEC);
    const before = readKnobs(bl, knobs);
    await optimizeBranches(bl, beam(), knobs, FAST);
    const after = readKnobs(bl, knobs);
    for (let i = 0; i < before.length; i++) {
      assert(before[i] === after[i], `knob ${i} moved: ${before[i]} -> ${after[i]}`);
    }
  });

  it('starts every branch from the same voltages', async () => {
    /*
      Otherwise each tuning inherits the previous one's answer, and the result
      for a branch depends on the order the branches happened to be visited.
      Running the list twice must therefore agree.
    */
    const bl = switched();
    const knobs = tunableKnobs(bl, SPEC);
    const first = await optimizeBranches(bl, beam(), knobs, FAST);
    const second = await optimizeBranches(bl, beam(), knobs, FAST);
    for (let i = 0; i < first.branches.length; i++) {
      const a = first.branches[i].settings;
      const b = second.branches[i].settings;
      for (let k = 0; k < a.length; k++) {
        assert(a[k] === b[k], `branch ${i} knob ${k} differed between runs: ${a[k]} vs ${b[k]}`);
      }
    }
  });

  it('delivers the beam where it said it would', async () => {
    // A branch's settings are only worth keeping if applying them reproduces
    // the transmission that was reported for them.
    const bl = switched();
    const knobs = tunableKnobs(bl, SPEC);
    const { branches } = await optimizeBranches(bl, beam(), knobs, FAST);
    const best = branches.reduce((a, b) => (b.transmitted > a.transmitted ? b : a));
    assert(best.transmitted > 0, 'at least one branch should receive the beam');

    writeKnobs(bl, knobs, best.settings);
    const check = scoreBeamline(bl, beam(), { target: best.end });
    assert(
      check.transmitted === best.transmitted,
      `applying "${best.label}" gave ${check.transmitted}, not the ${best.transmitted} reported`
    );
  });

  it('reports an unreachable branch rather than dropping it', async () => {
    // "No voltage I tried puts beam out of this port" is a finding about the
    // column, and usually the one worth knowing.
    const bl = switched();
    const knobs = tunableKnobs(bl, SPEC);
    const { branches } = await optimizeBranches(bl, beam(), knobs, FAST);
    assert(
      branches.length === bl.openEnds().length,
      'every end appears in the results, reachable or not'
    );
    for (const b of branches) {
      assert(Number.isFinite(b.transmitted), `${b.label} reports a number`);
      assert(b.count > 0, `${b.label} reports how many were tried`);
    }
  });

  it('names a deflector port by which way it goes', () => {
    const bl = switched();
    const ends = bl.openEnds();
    for (const e of ends) {
      const label = branchLabel(e);
      assert(typeof label === 'string' && label.length > 0, 'every end has a label');
    }
  });

  it('knows which ends lie below an element', () => {
    /*
      What "this branch" means when something is selected. A drift partway down
      one line leads to exactly one end, so tuning it is unambiguous; the
      deflector above it leads to all of them, so it is not, and the interface
      has to offer per-branch tuning there instead of guessing.
    */
    const bl = switched();
    const [entry, bend, bentLine, straightLine] = bl.elements;

    assert(bl.endsBelow(bentLine).length === 1, 'a leaf drift has one end below it');
    assert(bl.endsBelow(straightLine).length === 1, 'so does the other line');
    assert(
      bl.endsBelow(bentLine)[0].element === bentLine,
      'and that end is its own open exit'
    );

    const below = bl.endsBelow(bend);
    assert(below.length === bl.openEnds().length, `a deflector leads everywhere, got ${below.length}`);
    assert(
      bl.endsBelow(entry).length === below.length,
      'and so does everything above it'
    );
    assert(bl.endsBelow(null).length === 0, 'nothing selected leads nowhere');
  });

  it('tunes to one named end without touching the others', async () => {
    // What the per-branch Tune button does: the same search, scored against
    // one destination instead of the main line.
    const bl = switched();
    const knobs = tunableKnobs(bl, SPEC);
    const ends = bl.openEnds();

    for (const end of ends) {
      writeKnobs(bl, knobs, knobs.map(() => 0));
      const { optimizeVoltages } = await import('../src/optimize.js');
      const r = await optimizeVoltages(bl, beam(), knobs, {
        ...FAST,
        flight: { target: end },
      });
      // Whatever it found, the score it reports must be the score of the
      // setting it left behind, measured against the end it was aiming at.
      const check = scoreBeamline(bl, beam(), { target: end });
      assert(
        check.transmitted === r.transmitted,
        `tuning to ${branchLabel(end)} reported ${r.transmitted} but delivers ${check.transmitted}`
      );
    }
  });

  it('stops when asked', async () => {
    const bl = switched();
    const knobs = tunableKnobs(bl, SPEC);
    const { branches, cancelled } = await optimizeBranches(bl, beam(), knobs, {
      ...FAST,
      shouldStop: () => true,
    });
    assert(cancelled, 'it should report having been stopped');
    assert(branches.length === 0, `nothing should have been tuned, got ${branches.length}`);
  });
});
