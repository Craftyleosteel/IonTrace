/**
 * Beamline - a chain of elements that ions fly through.
 *
 * The column is a PATH, not a line. Each element carries a frame saying where
 * it sits and which way it faces, and the chain is built by starting the next
 * element wherever the last one left the beam. That is what makes elements
 * snap together, and it is also what lets a bender exist: a bender's exit
 * faces a different direction from its entrance, so everything after it turns
 * with it, without any element having to know it has been bent.
 *
 * A field query finds the element containing the point, transforms the point
 * into that element's local coordinates, asks it, and rotates the answer back.
 * Elements therefore never know where they are - which is precisely why they
 * are interchangeable.
 *
 * Misalignment
 * ------------
 * Each element may carry a small offset and tilt. These are applied to that
 * element's own placement and deliberately NOT propagated: a misaligned lens
 * does not move the ones downstream of it, because each is mounted
 * independently. Propagating them would model a bent optical bench rather
 * than a misaligned element.
 *
 * The modelling assumption, stated plainly
 * ----------------------------------------
 * **Elements are solved in isolation, not as one system.** Each is a separate
 * Dirichlet problem with grounded end faces, so its fringe field is confined
 * inside its own footprint. Measured against a single solve of a whole column:
 * the grounded faces BETWEEN elements cost almost nothing (0.004 V out of
 * 300), but each element's own margin clipping its OWN fringe costs about
 * 19 % at 0.8 bore radii, 12 % at 1.6, 6.6 % at 2.4 and 3.5 % at 3.2. So the
 * rule is not "leave a drift between elements" - it is "give each element
 * about three bore radii of its own margin". Elements warn when they do not.
 */

import {
  identityFrame,
  compose,
  toLocal,
  toGlobal,
  vectorToGlobal,
  misalignment,
  forwardOf,
} from './frames.js';
import { axisymmetricRuns, buildRunField } from './column.js';

/**
 * How far past its own outer radius an element may still claim a strike.
 *
 * Enough to cover a step that overshoots the wall, and enough for the corners
 * of a box whose diagonal exceeds its stated radius - a deflector's outer
 * radius is its half-width, and its corners sit a factor of root two further
 * out. Not enough to reach a neighbouring branch of a folded column, which is
 * the whole point.
 */
const CLAIM_MARGIN = 1.6;

/**
 * Is this local point in the element's free space?
 *
 * The answer is `e.strikes`, which for anything with a solved field means a
 * grid lookup - a hypotenuse, two roundings and an array index. That runs
 * about 170,000 times per flight of nine ions, and three quarters of those
 * calls are about a point comfortably inside the clear aperture, where the
 * answer is a foregone conclusion.
 *
 * So an element may declare `clearBore`: a transverse radius within which it
 * promises there is no metal anywhere along its length. Inside it, one
 * comparison settles the question. Outside it, or for an element that cannot
 * make that promise, the real test runs.
 *
 * A deflector cannot make it. Its electrodes wrap around the beam rather than
 * lying outside a cylinder, so a point can be well within r0 of the element's
 * axis and still be inside metal - which makes distance from that axis no
 * guide at all.
 */
function isFree(e, l) {
  const clear = e.clearBore;
  if (clear !== undefined && l[0] * l[0] + l[1] * l[1] <= clear * clear) return true;
  return !e.strikes(l[0], l[1], l[2]);
}

export class Beamline {
  constructor(elements = []) {
    this.elements = [];
    /**
     * Whether neighbouring elements are solved together.
     *
     * Off, every element is solved alone behind grounded end faces, which
     * both omits its fringe field and silently shields it. On, stretches of
     * axisymmetric elements share one grid and the field flows between them.
     * See src/column.js.
     *
     * Not a rendering option: it changes the field the ions fly through.
     */
    this.fringe = false;
    this.runs = [];
    this.runWarnings = [];
    for (const e of elements) this.add(e);
  }

  /**
   * The element the previous lookup landed in, as a hint for the next.
   *
   * Cleared whenever the column is re-laid, since an element's frame may have
   * moved underneath it. It is only ever an optimisation: every path that uses
   * it re-checks the point against that element before trusting it.
   */
  #hint = null;

  /**
   * Turn column solves on or off, re-solving as needed.
   *
   * Costs a Laplace solve per run, so it is called on a structural change or
   * when the setting is flipped - never from `layout`, which runs on every
   * drag.
   */
  setFringe(on) {
    this.fringe = Boolean(on);
    this.rebuildRuns();
  }

  rebuildRuns() {
    this.runs = [];
    this.runWarnings = [];
    if (!this.fringe) return;
    for (const run of axisymmetricRuns(this.elements)) {
      const built = buildRunField(this.elements, run);
      this.runs.push(built);
      this.runWarnings.push(...built.warnings);
    }
  }

  /** Re-apply element voltages to the column solves. No re-solve. */
  syncRuns() {
    for (const r of this.runs) r.sync();
  }

  /** The column solve covering an element, if there is one. */
  runFor(index) {
    if (!this.fringe) return null;
    return this.runs.find((r) => r.indices.includes(index)) ?? null;
  }

  /*
    The four structural changes. Each re-lays the column and, when column
    solves are on, re-solves them - which is why `rebuildRuns` lives here and
    not in `layout`: layout runs on every drag, and a Laplace solve must not.
  */

  /**
   * Put an element into the column.
   *
   * Topology is carried by `element.from = {parent, port}` - an object
   * reference, not an index, so splicing the array cannot silently re-wire the
   * tree. With no `attach` given the element is appended to the end of the
   * line that runs through whatever is currently last, which is what a plain
   * "add this" means when there is only one line.
   *
   * Anything already hanging from that port is pushed downstream and becomes
   * the new element's own child, so inserting in the middle of a branch splices
   * rather than truncates.
   */
  add(element, index = this.elements.length, attach = null) {
    element.align ??= { dx: 0, dy: 0, tiltX: 0, tiltY: 0 };

    const parent = attach ? attach.parent : this.elements[index - 1] ?? null;
    const port = attach ? attach.port : firstFreePort(this, parent);
    element.from = { parent, port };

    // Whatever was on that port now hangs off the new element instead.
    const displaced = parent ? this.childAt(parent, port) : this.roots()[0] ?? null;
    this.elements.splice(index, 0, element);
    if (displaced && displaced !== element) {
      displaced.from = { parent: element, port: exitsOf(element)[0].port };
    }

    this.layout();
    this.rebuildRuns();
    return element;
  }

  /**
   * Take an element out, joining what was above it to what was below.
   *
   * Its first child inherits its place on the parent's port. Any OTHER
   * children - a second branch off a deflector - lose their mount, so they are
   * removed with it. Silently orphaning a whole branch would be worse.
   */
  remove(index) {
    const removed = this.elements[index];
    if (!removed) return null;

    const kids = this.childrenOf(removed);
    const heir = kids[0] ?? null;
    if (heir) heir.from = { ...removed.from };
    for (const other of kids.slice(1)) this.#detach(other);

    this.elements.splice(this.elements.indexOf(removed), 1);
    this.layout();
    this.rebuildRuns();
    return removed;
  }

  /** Remove an element and everything hanging below it. */
  #detach(element) {
    for (const child of this.childrenOf(element)) this.#detach(child);
    const i = this.elements.indexOf(element);
    if (i >= 0) this.elements.splice(i, 1);
  }

  /**
   * Swap an element with the one before or after it on its own branch.
   *
   * Only along a branch: "upstream" of a junction is ambiguous, and swapping
   * across one would mean re-parenting a whole subtree into a different line.
   */
  move(index, delta) {
    const e = this.elements[index];
    if (!e) return false;

    const parent = e.from?.parent ?? null;
    if (delta < 0) {
      if (!parent) return false;
      const grand = parent.from?.parent ?? null;
      const kids = this.childrenOf(e);
      if (this.childrenOf(parent).length > 1) return false;
      // e takes the parent's mount; the parent hangs off e.
      e.from = { ...parent.from };
      parent.from = { parent: e, port: exitsOf(e)[0].port };
      for (const k of kids) k.from = { parent, port: exitsOf(parent)[0].port };
      void grand;
    } else {
      const next = this.childrenOf(e)[0];
      if (!next) return false;
      if (this.childrenOf(e).length > 1) return false;
      const kids = this.childrenOf(next);
      next.from = { ...e.from };
      e.from = { parent: next, port: exitsOf(next)[0].port };
      for (const k of kids) k.from = { parent: e, port: exitsOf(e)[0].port };
    }

    this.#reorder();
    this.layout();
    this.rebuildRuns();
    return true;
  }

  replace(index, element) {
    const old = this.elements[index];
    element.align ??= old?.align ?? { dx: 0, dy: 0, tiltX: 0, tiltY: 0 };
    element.from = old?.from ?? { parent: null, port: 'out' };
    // Children were mounted on the old object, so point them at the new one -
    // and at a port it actually has, since a rebuild can change the geometry
    // but never the kind of element.
    for (const child of this.childrenOf(old)) child.from.parent = element;
    this.elements[index] = element;
    this.layout();
    this.rebuildRuns();
    return element;
  }

  /**
   * Hang an element off a different exit, taking its line with it.
   *
   * Dropping a piece of hardware somewhere else is the one editing gesture a
   * tree needs that a list did not. Two rules keep it from producing something
   * that is not a column:
   *
   *   - nothing may be moved below itself, which would cut a loop out of the
   *     tree and leave it unreachable from the source;
   *   - a junction with both branches built is not moved, because there is no
   *     single answer to what happens to the branch that is left behind.
   *
   * Otherwise the gap it leaves heals - its own continuation takes its place -
   * and whatever was on the destination port is pushed down to hang off it.
   */
  reparent(element, parent, port) {
    if (!element || element === parent) return false;
    for (let p = parent; p; p = p.from?.parent) if (p === element) return false;
    if (this.childrenOf(element).length > 1) return false;

    const heir = this.childrenOf(element)[0] ?? null;
    if (heir) heir.from = { ...element.from };

    const displaced = parent
      ? this.childAt(parent, port)
      : this.roots().find((r) => r !== element) ?? null;

    element.from = { parent, port };
    if (displaced && displaced !== element) {
      displaced.from = { parent: element, port: exitsOf(element)[0].port };
    }

    this.#reorder();
    this.layout();
    this.rebuildRuns();
    return true;
  }

  /**
   * Take a set of elements that already know how they hang together.
   *
   * `add` decides an element's parent from where it is being inserted, which
   * is right for building a column by hand and wrong for loading one that has
   * already been built - it would overwrite the very structure being restored.
   * This takes the wiring as given, checks it forms one tree rooted somewhere,
   * and lays it out.
   */
  adopt(elements) {
    this.elements = elements.slice();
    for (const e of this.elements) {
      e.align ??= { dx: 0, dy: 0, tiltX: 0, tiltY: 0 };
      e.from ??= { parent: null, port: 'out' };
      // A parent outside this set would leave the element unreachable from the
      // source, which is not a column.
      if (e.from.parent && !this.elements.includes(e.from.parent)) {
        e.from = { parent: null, port: 'out' };
      }
    }
    this.#reorder();
    this.layout();
    this.rebuildRuns();
    return this;
  }

  /* ---------------------------------------------------------------- */
  /* topology                                                          */
  /* ---------------------------------------------------------------- */

  /** Elements with nothing upstream of them. Normally exactly one. */
  roots() {
    return this.elements.filter((e) => !e.from?.parent);
  }

  childrenOf(element) {
    if (!element) return this.roots();
    return this.elements.filter((e) => e.from?.parent === element);
  }

  childAt(element, port) {
    return this.elements.find((e) => e.from?.parent === element && e.from.port === port) ?? null;
  }

  /**
   * Was this point collected by a detector, rather than merely lost on metal?
   *
   * Both are strikes and neither is transmission, so the flight itself does
   * not need to tell them apart - but they mean opposite things to whoever is
   * reading the result, and reporting them as one number would hide the
   * distinction a detector exists to draw.
   */
  detected(x, y, z) {
    const hit = this.locate([x, y, z]);
    if (!hit || typeof hit.element.detects !== 'function') return false;
    return hit.element.detects(hit.local[0], hit.local[1], hit.local[2]);
  }

  /** The frame of a named open end, or the main exit if it is not one. */
  endFrame(target) {
    if (!target) return this.exitFrame;
    const found = this.openEnds().find(
      (o) => o.element === target.element && o.port === target.port
    );
    return found ? found.frame : this.exitFrame;
  }

  /** Every exit with nothing bolted to it: the open ends of the column. */
  openEnds() {
    const out = [];
    for (const e of this.elements) {
      for (const exit of exitsOf(e)) {
        if (!this.childAt(e, exit.port)) {
          out.push({ element: e, port: exit.port, frame: compose(e.nominalFrame, exit.transform) });
        }
      }
    }
    return out;
  }

  /** Put `elements` back into depth-first order, so display follows the beam. */
  #reorder() {
    const seen = new Set();
    const order = [];
    const walk = (e) => {
      if (!e || seen.has(e)) return;
      seen.add(e);
      order.push(e);
      for (const exit of exitsOf(e)) walk(this.childAt(e, exit.port));
    };
    for (const r of this.roots()) walk(r);
    // Anything unreachable keeps its place rather than vanishing.
    for (const e of this.elements) if (!seen.has(e)) order.push(e);
    this.elements = order;
  }

  /**
   * Walk the chain, giving every element its place.
   *
   * `nominalFrame` is where the element's mount puts it; `frame` is where it
   * actually is, once its own misalignment is applied. The chain continues
   * from the NOMINAL exit, so one element's error does not displace the rest
   * of the column.
   */
  layout() {
    // Depth first, because a branch continues from where its junction left
    // off. The recursion is over the tree, not the array: array order is for
    // display and is kept in step by `#reorder`.
    this.#hint = null;
    let longest = 0;
    let main = identityFrame();
    let mainEnd = null;

    const place = (e, cursor, path, isMain) => {
      e.nominalFrame = cursor;
      e.frame = compose(cursor, misalignment(e.align));
      // Distance along the REFERENCE PATH, not along z. For a straight column
      // the two coincide, which is why the name survives; once a bender is in
      // the line, z stops being meaningful and arc length is what orders the
      // elements. On a branching column it is distance from the source along
      // this element's OWN branch.
      e.zStart = path;
      e.zEnd = path + e.length;

      const exits = exitsOf(e);
      for (let k = 0; k < exits.length; k++) {
        const exit = exits[k];
        const next = compose(cursor, exit.transform);
        const child = this.childAt(e, exit.port);
        // The main line is whatever you reach by always taking the first exit,
        // which for a deflector is the bend - the reason it is there.
        const onMain = isMain && k === 0;
        if (child) {
          place(child, next, path + exit.length, onMain);
        } else if (onMain) {
          main = next;
          mainEnd = { element: e, port: exit.port };
          longest = Math.max(longest, path + exit.length);
        } else {
          longest = Math.max(longest, path + exit.length);
        }
      }
      if (exits.length === 0) longest = Math.max(longest, e.zEnd);
    };

    for (const r of this.roots()) place(r, identityFrame(), 0, r === this.roots()[0]);

    this.exitFrame = main;
    /**
     * The open end the column is "aimed at".
     *
     * Reached by always taking an element's first exit, which for a deflector
     * is the bend - the reason it is in the line at all. A branching column
     * has several open ends, so anything that has to mean "out" rather than
     * "out of something" needs one of them named, and this is it.
     */
    this.mainEnd = mainEnd;
    this.length = longest;
    return this;
  }

  /**
   * Which open end an ion left through, by proximity.
   *
   * Not by "is it past the plane": on a branching column an ion can be past
   * more than one exit plane at once - the two ends of a switch face different
   * ways, and a point far downstream of one may be nominally ahead of the
   * other as well. Where it actually came out is the end it is nearest.
   */
  endNearest(x, y, z) {
    let best = null;
    let bestD = Infinity;
    for (const end of this.openEnds()) {
      const o = end.frame.o;
      const d = (x - o[0]) ** 2 + ((y ?? 0) - o[1]) ** 2 + (z - o[2]) ** 2;
      if (d < bestD) {
        bestD = d;
        best = end;
      }
    }
    return best;
  }

  /**
   * The element containing a global point, with the point already expressed
   * in that element's local coordinates.
   *
   * O(N) rather than a sorted lookup, because once the path can bend and
   * elements can be nudged off their mounts there is no single coordinate to
   * sort on. For the handful of elements a column has, this is nothing.
   */
  /**
   * The element at a given distance along the reference path.
   *
   * Ordering by arc length rather than by z, because z stops being monotonic
   * as soon as the column bends back on itself.
   */
  elementAt(s) {
    for (const e of this.elements) {
      if (s >= e.zStart && s < e.zEnd) return e;
    }
    // Ambiguous once the column branches - two elements on different branches
    // sit at the same distance from the source - so this returns the first
    // match and is only used where a single line is assumed.
    const last = this.elements[this.elements.length - 1];
    if (last && s === last.zEnd) return last;
    return null;
  }

  /**
   * Which element a global point belongs to.
   *
   * Elements answer `contains` on their AXIAL extent alone: a point inside an
   * element's length but outside its bore is still that element's business,
   * because it is about to be a strike on that element's wall rather than an
   * ion that has wandered out of the column. That is right, and it is why
   * this cannot simply take the first element that says yes.
   *
   * In a straight column it could, because axial position identifies an
   * element uniquely. Once the path bends it does not. Fold a column through
   * two right angles and its last drift runs back alongside its first, inside
   * the first one's axial range but eighty millimetres off its axis - and
   * every ion entering the last drift was being reported as striking the wall
   * of the first. Nothing about that looks like a lookup failure from the
   * outside: the beam simply stops, at a plausible place, for a plausible
   * reason.
   *
   * So an element claims a point only if the point is in its free space.
   * Whichever element holds it in vacuum wins. Failing that it may still be
   * claimed as a strike, but only within reach of the element's own hardware:
   * an ion eighty-seven millimetres off the axis of a lens whose housing is
   * fourteen millimetres across did not hit that lens, and saying so killed
   * beams that had already flown the whole column and reached the exit.
   * Beyond every element's envelope the point belongs to no element, which is
   * what lets `classify` call it an exit rather than a crash.
   */
  locate(g) {
    /*
      Try the element the last lookup found, before scanning for one.

      This is the hottest path in the simulator by a wide margin. A single
      flight of nine ions runs about eight thousand integrator steps, and each
      step asks for the field four times per ion and for a strike test once -
      six hundred thousand lookups, every one of which was scanning the whole
      column and calling both `contains` and `strikes` on each element it
      passed.

      Almost all of those ask about a point a fraction of a millimetre from the
      previous one, or about another ion in the same beam a millimetre away, so
      the answer is nearly always the same element as last time. The hint is
      only ever a guess: it is accepted only when the point is in that
      element's free space, which is the same condition the scan below would
      have returned it for, and otherwise the scan runs as usual.
    */
    const hint = this.#hint;
    if (hint) {
      const e = hint.element;
      const l = toLocal(e.frame, g);
      if (e.contains(l[0], l[1], l[2]) && isFree(e, l)) {
        return { element: e, local: l, index: hint.index, free: true };
      }
    }

    let fallback = null;
    for (let i = 0; i < this.elements.length; i++) {
      const e = this.elements[i];
      const l = toLocal(e.frame, g);
      if (!e.contains(l[0], l[1], l[2])) continue;
      if (isFree(e, l)) {
        this.#hint = { element: e, index: i };
        return { element: e, local: l, index: i, free: true };
      }
      const reach = Math.hypot(l[0], l[1]) <= e.outerRadius * CLAIM_MARGIN;
      if (reach && !fallback) fallback = { element: e, local: l, index: i, free: false };
    }
    return fallback;
  }

  /* ---------------------------------------------------------------- */
  /* the interface the integrator flies through                       */
  /* ---------------------------------------------------------------- */

  /**
   * Local coordinates within the column solve covering `hit`, if any.
   *
   * A run spans several elements, so the point has to be expressed relative
   * to the run's own start rather than the element's. Every element in a run
   * is axisymmetric and unmisaligned, so their frames differ only by a shift
   * along the axis, and that shift is the difference of their path positions.
   */
  #runLocal(hit) {
    const run = this.runFor(hit.index);
    if (!run) return null;
    const zl = hit.local[2] + (hit.element.zStart - run.z0);
    return { run, x: hit.local[0], y: hit.local[1], z: zl };
  }

  fieldAt3D(x, y, z, t) {
    const hit = this.locate([x, y, z]);
    if (!hit) return { Ex: 0, Ey: 0, Ez: 0 };
    const { element, local } = hit;

    const r = this.#runLocal(hit);
    const e = r
      ? r.run.field.fieldAt3D(r.x, r.y, r.z)
      : element.fieldAt(local[0], local[1], local[2], t);
    const g = vectorToGlobal(element.frame, [e.Ex, e.Ey, e.Ez]);
    return { Ex: g[0], Ey: g[1], Ez: g[2] };
  }

  potentialAt3D(x, y, z, t) {
    const hit = this.locate([x, y, z]);
    if (!hit) return 0;
    const r = this.#runLocal(hit);
    if (r) return r.run.field.potentialAt3D(r.x, r.y, r.z);
    return hit.element.potentialAt(hit.local[0], hit.local[1], hit.local[2], t);
  }

  strikes(x, y, z) {
    const hit = this.locate([x, y, z]);
    if (!hit) return false;
    // `locate` already had to decide this to choose the element at all, so
    // asking again would be the same grid lookup twice for every step.
    if (hit.free) return false;
    return hit.element.strikes(hit.local[0], hit.local[1], hit.local[2]);
  }

  /**
   * Where an ion is relative to the column: inside it, before it, past it, or
   * wandered out of it altogether.
   *
   * With a straight contiguous line this was just a z comparison. Once the
   * path bends and elements can be misaligned, "past the end" means past the
   * exit PLANE - on the far side of the last element's exit face - and there
   * is a fourth possibility that did not exist before: an ion can leave
   * through a gap opened up by a misalignment without ever reaching either
   * end. Calling that "exited" would quietly count a lost ion as transmitted.
   */
  classify(x, y, z) {
    if (this.locate([x, y, z])) return 'inside';
    const g = [x, y, z];

    // Past ANY open end counts as having left. A branching column has several,
    // and a deflector used as a switch has one on whichever port nothing is
    // bolted to - an ion going that way has genuinely left the instrument, not
    // wandered out of it.
    for (const end of this.openEnds()) {
      const o = end.frame.o;
      const f = forwardOf(end.frame);
      const ahead = (g[0] - o[0]) * f[0] + (g[1] - o[1]) * f[1] + (g[2] - o[2]) * f[2];
      if (ahead >= 0) return 'exited';
    }

    const first = this.roots()[0];
    if (first) {
      const entryO = first.frame.o;
      const f = forwardOf(first.frame);
      const behind =
        (g[0] - entryO[0]) * f[0] +
        (g[1] - entryO[1]) * f[1] +
        (g[2] - entryO[2]) * f[2];
      if (behind <= 0) return 'reflected';
    }

    return 'lost';
  }

  get radiusLimit() {
    return Math.max(1e-3, ...this.elements.map((e) => e.outerRadius));
  }

  get lengthScale() {
    const scales = this.elements
      .map((e) => e.lengthScale)
      .filter((s) => Number.isFinite(s) && s > 0);
    return scales.length ? Math.min(...scales) : 1e-3;
  }

  get shortestPeriod() {
    const periods = this.elements
      .map((e) => e.shortestPeriod)
      .filter((p) => p != null && p > 0);
    return periods.length ? Math.min(...periods) : null;
  }

  /* ---------------------------------------------------------------- */
  /* diagnostics                                                      */
  /* ---------------------------------------------------------------- */

  get warnings() {
    const out = [];
    for (const e of this.elements) {
      for (const w of e.warnings ?? []) out.push(`${e.label}: ${w}`);
    }

    out.push(...this.runWarnings);

    for (let i = 1; i < this.elements.length; i++) {
      const a = this.elements[i - 1];
      const b = this.elements[i];
      if (a.type === 'drift' || b.type === 'drift') continue;
      // With column solves on, a shared grid IS a true solution for the pair -
      // that is the whole point of it - so this only applies to neighbours the
      // run did not cover.
      if (this.runFor(i) && this.runFor(i) === this.runFor(i - 1)) continue;
      out.push(
        `${a.label} and ${b.label} are adjacent with no drift between them. ` +
          'Each was solved in isolation with grounded end faces, so the field ' +
          'at their junction is not a true solution for the pair.'
      );
    }

    for (let i = 1; i < this.elements.length; i++) {
      const a = this.elements[i - 1];
      const b = this.elements[i];
      if (b.bore < a.bore * 0.5) {
        out.push(
          `${b.label} has less than half the bore of ${a.label} before it; ` +
            'the aperture step will clip the beam.'
        );
      }
    }

    // A misalignment big enough to open a gap the beam can escape through.
    for (const e of this.elements) {
      const off = Math.hypot(e.align.dx, e.align.dy);
      if (off > e.bore * 0.5) {
        out.push(
          `${e.label} is offset by ${(off * 1e3).toFixed(2)} mm, more than half ` +
            'its own bore. The beam will clip its aperture.'
        );
      }
    }

    return out;
  }

  /** Reset every element onto its nominal mount. */
  autoAlign() {
    for (const e of this.elements) {
      e.align = { dx: 0, dy: 0, tiltX: 0, tiltY: 0 };
    }
    return this.layout();
  }

  /** True if anything is off its nominal placement. */
  get misaligned() {
    return this.elements.some(
      (e) =>
        e.align.dx !== 0 ||
        e.align.dy !== 0 ||
        e.align.tiltX !== 0 ||
        e.align.tiltY !== 0
    );
  }

  /**
   * Drawable outline in GLOBAL coordinates.
   *
   * Each element's rectangles are expressed in its own frame as
   * (along, transverse) pairs; here they become quadrilaterals in space, so a
   * rotated or bent element draws correctly without the renderer knowing
   * anything about bends.
   */
  outline() {
    const out = [];
    for (const e of this.elements) {
      for (const shape of shapesOf(e)) {
        out.push({
          corners: shape.points.map((pt) => toGlobal(e.frame, pt)),
          // The same metal seen edge-on in the other plane. Only meaningful
          // for axisymmetric elements, where it really is the same shape.
          cornersRolled: shape.axisymmetric
            ? shape.points.map((pt) => toGlobal(e.frame, [0, pt[0], pt[2]]))
            : null,
          ghost: !!shape.ghost,
          wall: !!shape.wall,
          element: e,
        });
      }
    }
    return out;
  }

  /**
   * The reference path, as one polyline per route from the source to an open
   * end.
   *
   * Separate polylines, because a branching column has no single path and
   * drawing one flat list of points joins the end of one branch to the start
   * of the next with a line that is not there. Routes that share a trunk draw
   * it twice, which costs nothing and keeps this simple.
   *
   * Which curve an element contributes depends on the exit being taken: a
   * deflector's bent path is a quarter arc, its straight path a chord.
   */
  centreLines(perElement = 12) {
    const lines = [];
    const walk = (e, acc) => {
      for (const exit of exitsOf(e)) {
        const n = exit.path ? perElement : 1;
        const seg = [];
        for (let k = 0; k <= n; k++) {
          const f = k / n;
          seg.push(toGlobal(e.frame, exit.path ? exit.path(f) : [0, 0, f * exit.length]));
        }
        const line = acc.concat(seg);
        const child = this.childAt(e, exit.port);
        if (child) walk(child, line);
        else lines.push(line);
      }
    };
    for (const r of this.roots()) walk(r, []);
    return lines;
  }

  /** Every reference-path point, flattened. For bounds, which needs no order. */
  centreLine(perElement = 12) {
    return this.centreLines(perElement).flat();
  }

  /** Bounds of everything drawn, in all three global axes. */
  bounds() {
    const lo = [Infinity, Infinity, Infinity];
    const hi = [-Infinity, -Infinity, -Infinity];
    const see = (p) => {
      for (let i = 0; i < 3; i++) {
        lo[i] = Math.min(lo[i], p[i]);
        hi[i] = Math.max(hi[i], p[i]);
      }
    };
    for (const r of this.outline()) {
      for (const c of r.corners) see(c);
      if (r.cornersRolled) for (const c of r.cornersRolled) see(c);
    }
    for (const p of this.centreLine()) see(p);
    if (!Number.isFinite(lo[0])) {
      return { minX: -0.01, maxX: 0.01, minY: -0.01, maxY: 0.01, minZ: 0, maxZ: 0.1 };
    }
    return {
      minX: lo[0], maxX: hi[0],
      minY: lo[1], maxY: hi[1],
      minZ: lo[2], maxZ: hi[2],
    };
  }

  /**
   * Does the column leave the horizontal plane?
   *
   * The view uses this to decide whether a side elevation is worth showing.
   * A straight or horizontally-bent line is fully described by the top view,
   * and a second empty pane would be wasted space.
   */
  get usesVerticalPlane() {
    // Measured on the PATH, not on the bounding box. A column's vertical
    // extent is partly just how tall its hardware is, and for a quadrupole
    // deflector the two are the same number - it is exactly as wide across
    // the bend as the bend displaces the beam - so a size-based test cannot
    // tell a vertical bend from a horizontal one at all.
    const frames = [...this.elements.map((e) => e.frame), this.exitFrame];

    // A bend out of the plane shows up exactly in the direction the path
    // points, so this needs no tolerance beyond round-off.
    if (frames.some((f) => Math.abs(forwardOf(f)[1]) > 1e-9)) return true;

    // Failing that, only a displacement large enough to be worth drawing
    // counts - otherwise every fraction of a millimetre of misalignment would
    // open a second pane to show it in.
    const spread = Math.max(0, ...frames.map((f) => Math.abs(f.o[1])));
    return spread > this.radiusLimit * 0.5;
  }
}

/** The exit placement of a straight element of the given length. */
export function straightExit(length) {
  return { o: [0, 0, length], m: [1, 0, 0, 0, 1, 0, 0, 0, 1] };
}

/**
 * The ways out of an element.
 *
 * Almost everything has one, and says so by not mentioning it: a lens passes
 * the beam straight on, displaced by its own length. A deflector declares two,
 * because it is a junction - see `exits` in elements/bender.js.
 */
export function exitsOf(element) {
  if (Array.isArray(element?.exits) && element.exits.length) return element.exits;
  return [
    {
      port: 'out',
      label: 'Out',
      length: element.length,
      transform: element.exitTransform ?? straightExit(element.length),
    },
  ];
}

/** The first exit of `parent` with nothing already on it. */
function firstFreePort(beamline, parent) {
  if (!parent) return 'out';
  const exits = exitsOf(parent);
  const free = exits.find((x) => !beamline.childAt(parent, x.port));
  return (free ?? exits[0]).port;
}

/**
 * An element's drawable shapes, as polygons of 3D points in its own frame.
 *
 * Most elements are axisymmetric and describe themselves as a few boxes in
 * the x-z plane, which are mirrored about the axis to give the familiar
 * two-sided cross-section. A bender is neither axisymmetric nor straight, and
 * once it can be rolled it does not even lie in the x-z plane - so it
 * supplies full 3D polygons directly. Going through polygons rather than
 * boxes is what lets the renderer stay ignorant of which kind it is drawing,
 * and what lets the same shapes be projected into more than one view.
 */
function shapesOf(e) {
  if (typeof e.shapes === 'function') return e.shapes();
  const out = [];
  for (const r of e.rects ?? []) {
    for (const sign of [1, -1]) {
      out.push({
        points: [
          [sign * r.r0, 0, r.z0],
          [sign * r.r1, 0, r.z0],
          [sign * r.r1, 0, r.z1],
          [sign * r.r0, 0, r.z1],
        ],
        ghost: r.ghost,
        wall: r.wall,
        // Axisymmetric metal looks the same in any plane containing the axis,
        // so the side view draws it by rotating the same box a quarter turn.
        axisymmetric: true,
      });
    }
  }
  return out;
}
