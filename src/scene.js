/**
 * Saving and reloading a column.
 *
 * Why JSON rather than CSV
 * ------------------------
 * A CSV is a table, and a beamline is a tree. You can force one into the
 * other - a parent column and a port column, with every element's settings
 * flattened into a row whose meaning depends on its type, so a deflector's
 * "voltage" and a funnel's "rings" share a column or sprawl into thirty mostly
 * empty ones. Either way the file stops being readable by the spreadsheet that
 * was the only reason to pick CSV, and the parser has to reconstruct the
 * hierarchy anyway.
 *
 * JSON holds the shape directly, every browser parses it without a library,
 * and a person can open it and see what it says. So: JSON.
 *
 * What is saved, and what is not
 * ------------------------------
 * The saved file is the DESCRIPTION of a column, not its solution. Element
 * types and their parameters, how they hang off one another, each one's
 * misalignment, the ion in the source and the two physics switches. Nothing
 * that was computed from those: no fields, no basis solutions, no
 * trajectories. Reloading re-solves, which takes a moment and means a file
 * written by an older version gets today's physics rather than yesterday's
 * cached numbers.
 *
 * Parameters are stored by name and merged over the element's current
 * defaults, so a file that predates a new parameter still loads - the new one
 * takes its default rather than becoming undefined.
 *
 * That merge is deliberately one-sided, and the other side is reported rather
 * than swallowed. A parameter the FILE has and this version does not means it
 * was renamed or removed, and the merge cannot tell that from an ordinary
 * version bump: the old key rides along unread, the new key takes its default,
 * and the element loads looking fine while describing different hardware. So
 * `restore` lists those keys as problems. Likewise a file naming an element
 * type that no longer exists, or hanging an element on a port that has since
 * gone. The rule throughout is that a file may load imperfectly, but it may
 * not load imperfectly and quietly.
 */

import { ELEMENT_TYPES, createElement } from './elements/index.js';
import { exitsOf } from './beamline.js';

export const SCENE_FORMAT = 'iontrace.column';
export const SCENE_VERSION = 1;

/**
 * A column as a plain object.
 *
 * Parents are referenced by index into `elements`, which is in the order the
 * beam visits them, so a parent always appears before its children and a
 * reader can build the tree in one pass.
 */
export function serialise(beamline, extras = {}) {
  const index = new Map(beamline.elements.map((e, i) => [e, i]));
  return {
    format: SCENE_FORMAT,
    version: SCENE_VERSION,
    saved: new Date().toISOString(),
    source: extras.source ?? null,
    physics: {
      repulsion: extras.repulsion ?? 'none',
      beamCurrent: extras.beamCurrent ?? null,
      ionsPerParticle: extras.ionsPerParticle ?? null,
      fringe: Boolean(beamline.fringe),
    },
    elements: beamline.elements.map((e) => ({
      type: e.typeKey,
      // A copy, not a reference: `params` is live and would otherwise be
      // mutated under a caller holding the serialised object.
      params: { ...e.params },
      parent: e.from?.parent ? index.get(e.from.parent) : null,
      port: e.from?.port ?? null,
      align: { ...e.align },
    })),
  };
}

/** Everything wrong with a file, as sentences rather than a thrown error. */
function validate(data) {
  const problems = [];
  if (!data || typeof data !== 'object') return ['That file does not contain a column.'];
  if (data.format !== SCENE_FORMAT) {
    problems.push(`This is not an IonTrace column file (it says "${data.format}").`);
  }
  if (!Array.isArray(data.elements)) {
    problems.push('It lists no elements.');
    return problems;
  }
  data.elements.forEach((e, i) => {
    if (!ELEMENT_TYPES[e?.type]) {
      problems.push(`Element ${i + 1} is a "${e?.type}", which this version does not have.`);
    }
    if (e?.parent != null && (e.parent < 0 || e.parent >= data.elements.length)) {
      problems.push(`Element ${i + 1} hangs off an element that is not in the file.`);
    }
    if (e?.parent != null && e.parent >= i) {
      problems.push(`Element ${i + 1} hangs off one that comes after it.`);
    }
  });
  return problems;
}

/**
 * Rebuild the elements of a saved column.
 *
 * Returns them already wired to one another, ready to be handed to a Beamline;
 * the caller owns everything outside the column, since the ion source and the
 * physics switches live in the interface rather than in the beamline object.
 *
 * @returns {{elements: object[], source: object|null, physics: object,
 *            problems: string[]}}
 */
export function restore(data) {
  const problems = validate(data);
  if (problems.length) return { elements: [], source: null, physics: {}, problems };

  const built = data.elements.map((saved, i) => {
    // Merged over the current defaults, so a file written before a parameter
    // existed still loads and that parameter takes its default.
    const spec = ELEMENT_TYPES[saved.type];
    const e = createElement(saved.type, { ...spec.defaults, ...saved.params });
    e.align = { dx: 0, dy: 0, tiltX: 0, tiltY: 0, ...saved.align };

    /*
      A parameter the file has and this version does not.

      Merging over defaults quietly handles the case a version bump usually
      produces - a NEW parameter, absent from the file, taking its default. It
      does not handle a parameter being RENAMED or dropped, and that case looks
      identical from here: the old key rides along and is ignored, the new key
      takes its default, the file loads without complaint and the element comes
      back with different geometry from the one that was saved.

      Silently different is the worst of the available outcomes, so say it. The
      element is still built - the file is not wrong, it is just older than the
      code - but the reload is not faithful and the person deserves to know
      which setting stopped being honoured.
    */
    const known = Object.keys(spec.defaults ?? {});
    const stale = Object.keys(saved.params ?? {}).filter((k) => !known.includes(k));
    if (stale.length) {
      problems.push(
        `Element ${i + 1} (${saved.type}) was saved with ${stale.join(', ')}, which ` +
          'this version no longer has; it now uses the current default instead, so ' +
          'its geometry may differ from the saved one.'
      );
    }
    return e;
  });

  built.forEach((e, i) => {
    const saved = data.elements[i];
    const parent = saved.parent == null ? null : built[saved.parent];
    let port = saved.port;
    if (parent) {
      const ports = exitsOf(parent).map((x) => x.port);
      if (!ports.includes(port)) {
        // The element it hangs off has different exits from when this was
        // saved - a deflector gained a third, say. Fall back to its first
        // rather than leaving the element mounted on nothing.
        problems.push(
          `Element ${i + 1} was on a "${port}" exit that no longer exists; ` +
            `moved to "${ports[0]}".`
        );
        port = ports[0];
      }
    }
    e.from = { parent, port: parent ? port : 'out' };
  });

  return {
    elements: built,
    source: data.source ?? null,
    physics: data.physics ?? {},
    problems,
  };
}
