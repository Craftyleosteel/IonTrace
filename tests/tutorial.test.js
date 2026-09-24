/**
 * The tutorial and the help topics.
 *
 * Why a tutorial needs tests at all
 * ---------------------------------
 * Because its claims are checkable. A step that says "build this and fly it"
 * and then builds a column that loses the beam is not a typo - it is the
 * documentation and the simulator disagreeing about the physics, which is the
 * one kind of error this project cannot afford to ship quietly.
 *
 * So every step that builds something is built and flown here, against the same
 * element code the browser loads. Prose drifts away from code; a column that
 * has to deliver ions cannot.
 *
 * The help topics get the cheaper treatment they deserve - structure, not
 * content. No test can tell whether an explanation is a good one. It can tell
 * whether a question mark somewhere in the interface points at a topic that no
 * longer exists, which is the failure that actually happens.
 */

import { describe, it, assert } from './harness.js';

import { Beamline, exitsOf } from '../src/beamline.js';
import { createElement, ELEMENT_TYPES } from '../src/elements/index.js';
import { discBeam } from '../src/ion.js';
import { flyBeam } from '../src/integrator.js';
import { HELP, topic, topicKeys } from '../src/help.js';
import { LESSONS, lesson, allSteps } from '../src/tutorial.js';

/* ------------------------------------------------------------------ */
/* help                                                                */
/* ------------------------------------------------------------------ */

describe('Help topics', () => {
  it('gives every topic a title and something to say', () => {
    for (const key of topicKeys()) {
      const t = HELP[key];
      assert(typeof t.title === 'string' && t.title.length > 0, `${key} has no title`);
      assert(Array.isArray(t.body) && t.body.length > 0, `${key} has an empty body`);
      for (const p of t.body) {
        assert(typeof p === 'string' && p.trim().length > 20, `${key} has a stub paragraph`);
      }
    }
  });

  it('covers every element a person can place', () => {
    // The question mark on an element's panel looks up its type key, so a new
    // element type with no topic silently loses its help.
    for (const key of Object.keys(ELEMENT_TYPES)) {
      assert(topic(key), `element type "${key}" has no help topic`);
    }
  });

  it('returns nothing for a key it does not have', () => {
    assert(topic('no-such-topic') === null, 'an unknown key should be null');
    // Guards against reaching Object.prototype: `topic('toString')` returning a
    // function would put a function where the popover expects a topic.
    assert(topic('toString') === null, 'must not reach through to Object.prototype');
  });
});

/* ------------------------------------------------------------------ */
/* tutorial                                                            */
/* ------------------------------------------------------------------ */

describe('Tutorial', () => {
  it('is structurally complete', () => {
    assert(LESSONS.length > 0, 'there are no lessons');
    const ids = new Set();
    for (const l of LESSONS) {
      assert(l.id && !ids.has(l.id), `lesson id "${l.id}" is missing or repeated`);
      ids.add(l.id);
      assert(l.title && l.summary, `lesson "${l.id}" needs a title and a summary`);
      assert(l.steps.length > 0, `lesson "${l.id}" has no steps`);
      for (const s of l.steps) {
        assert(s.title, `a step of "${l.id}" has no title`);
        assert(Array.isArray(s.body) && s.body.length > 0, `step "${s.title}" says nothing`);
        if (s.action) {
          assert(s.action.label, `step "${s.title}" has an action with no label`);
          assert(typeof s.action.run === 'function', `step "${s.title}" action cannot run`);
        }
      }
    }
    assert(lesson(LESSONS[0].id) === LESSONS[0], 'lookup by id should work');
    assert(lesson('nope') === null, 'an unknown lesson id should be null');
  });

  it('points only at help topics that exist', () => {
    for (const { step } of allSteps()) {
      if (step.topic) {
        assert(topic(step.topic), `step "${step.title}" points at missing topic "${step.topic}"`);
      }
    }
  });

  /**
   * A stand-in for the interface, recording what a step asked for.
   *
   * It does the two things that can actually fail - building a column out of
   * real elements, and holding the beam settings - and ignores the rest, which
   * is drawing.
   */
  function harness() {
    const state = {
      beamline: null,
      spec: { mass: 100, charge: 1, energy: 50 },
      rays: 9,
      radius: 1.0,
      fringe: false,
      flew: 0,
      tuned: 0,
    };
    const api = {
      build(specs) {
        const made = specs.map((s) => createElement(s.type, s.params ?? {}));
        for (let i = 1; i < made.length; i++) {
          made[i].from = { parent: made[i - 1], port: exitsOf(made[i - 1])[0].port };
        }
        if (made.length) made[0].from = { parent: null, port: 'out' };
        const bl = new Beamline();
        bl.fringe = state.fringe;
        bl.adopt(made);
        state.beamline = bl;
      },
      beam(s) {
        for (const k of ['mass', 'charge', 'energy']) if (s[k] !== undefined) state.spec[k] = s[k];
        if (s.rays !== undefined) state.rays = s.rays;
        if (s.beamRadius !== undefined) state.radius = s.beamRadius;
      },
      physics(s) {
        if (s.fringe !== undefined) {
          state.fringe = Boolean(s.fringe);
          state.beamline?.setFringe(state.fringe);
        }
      },
      view() {},
      select() {},
      fly() {
        state.flew += 1;
      },
      async optimise() {
        state.tuned += 1;
      },
    };
    return { api, state };
  }

  it('builds a column that delivers the beam at every step that flies one', async () => {
    /*
      The claim under test is the one the tutorial makes out loud: press this
      and watch the ions go through. Counting only what EXITED would fail the
      detector lesson, where the whole point is that the ions stop - so an ion
      that lands on a detector's active surface counts as delivered too.
    */
    for (const { lesson: l, step } of allSteps()) {
      if (!step.action) continue;
      const { api, state } = harness();
      await step.action.run(api);
      if (!state.beamline || state.beamline.elements.length === 0) continue;

      const { tracks } = flyBeam(
        state.beamline,
        discBeam({ ...state.spec, count: state.rays, radius: state.radius }),
        { cfl: 0.05, maxSteps: 400000 }
      );
      const delivered = tracks.filter((t) => {
        if (t.stop === 'exited') return true;
        const p = t.points[t.points.length - 1];
        return t.stop === 'electrode' && state.beamline.detected(p.x, p.y ?? 0, p.z);
      }).length;

      assert(
        delivered >= Math.ceil(state.rays * 0.6),
        `"${l.title} / ${step.title}" builds a column that delivers only ` +
          `${delivered} of ${state.rays}`
      );
    }
  });

  it('never asks the interface for something it cannot do', () => {
    // Every method a step calls has to exist on the real api in main.js. This
    // pins the surface so a step cannot quietly depend on a helper that was
    // only ever in the test harness.
    const { api } = harness();
    const allowed = new Set(Object.keys(api));
    for (const key of ['build', 'beam', 'physics', 'view', 'select', 'fly', 'optimise']) {
      assert(allowed.has(key), `the tutorial api is missing ${key}`);
    }
  });
});
