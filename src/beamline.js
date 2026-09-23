/**
 * Beamline - an ordered column of elements that ions fly through.
 *
 * Each element carries its own local field solve and its own local
 * coordinates. The beamline lays them end to end, and at any point in space it
 * finds the element containing that z and asks it, translating the coordinate.
 * This is what makes the elements composable: adding a new kind of optic means
 * describing its metal and its field, not touching anything here.
 *
 * The modelling assumption, stated plainly
 * ----------------------------------------
 * **Elements are solved in isolation, not as one system.** Each is a separate
 * Dirichlet problem with grounded end faces, so its fringe field is confined
 * inside its own footprint and stops abruptly at the boundary. The true
 * solution for the whole column would let neighbouring electrodes see each
 * other.
 *
 * How good is that? Measured, by building the same column both ways and
 * comparing the on-axis potential:
 *
 *   - The grounded end faces BETWEEN elements cost almost nothing, provided
 *     each element's own margins are adequate: splitting a lens into
 *     drift + lens + drift rather than one wider lens changed the on-axis
 *     potential by 0.004 V out of 300 V.
 *   - What actually costs accuracy is each element's own margin clipping its
 *     OWN fringe field. The error decays roughly exponentially in
 *     margin / bore: about 19 % at 0.8 bore radii, 12 % at 1.6, 6.6 % at 2.4
 *     and 3.5 % at 3.2.
 *
 * So the rule is not "leave a drift between elements" - it is "give each
 * element about three bore radii of its own margin". Elements warn when they
 * do not have it. A drift between two live elements still helps, because a
 * real grounded drift tube genuinely shields, which is exactly what the
 * isolation approximation is pretending.
 *
 * The alternative, a single solve over the whole line, is what a 3D code
 * does; it is far more expensive and is not what this build does.
 *
 * Solving the whole column at once would also destroy the property that makes
 * this interactive: each element re-solves only when ITS geometry changes, and
 * voltages never re-solve at all.
 */

export class Beamline {
  constructor(elements = []) {
    this.elements = [];
    for (const e of elements) this.add(e);
  }

  /** Append an element and re-lay the column. */
  add(element, index = this.elements.length) {
    this.elements.splice(index, 0, element);
    this.layout();
    return element;
  }

  /** Remove the element at `index`. */
  remove(index) {
    const [removed] = this.elements.splice(index, 1);
    this.layout();
    return removed;
  }

  /** Move an element one place along the column. */
  move(index, delta) {
    const target = index + delta;
    if (target < 0 || target >= this.elements.length) return false;
    const [e] = this.elements.splice(index, 1);
    this.elements.splice(target, 0, e);
    this.layout();
    return true;
  }

  replace(index, element) {
    this.elements[index] = element;
    this.layout();
    return element;
  }

  /** Assign each element its axial start, in metres. */
  layout() {
    let z = 0;
    for (const e of this.elements) {
      e.zStart = z;
      e.zEnd = z + e.length;
      z = e.zEnd;
    }
    this.length = z;
    return this;
  }

  /** The element containing this z, or null beyond the ends. */
  elementAt(z) {
    for (const e of this.elements) {
      if (z >= e.zStart && z < e.zEnd) return e;
    }
    // The very end of the column belongs to its last element rather than to
    // nothing, so an ion exactly on the exit plane is still inside the optic.
    const last = this.elements[this.elements.length - 1];
    if (last && z === last.zEnd) return last;
    return null;
  }

  /* ---------------------------------------------------------------- */
  /* the interface the integrator flies through                       */
  /* ---------------------------------------------------------------- */

  fieldAt3D(x, y, z, t) {
    const e = this.elementAt(z);
    if (!e) return { Ex: 0, Ey: 0, Ez: 0 };
    return e.fieldAt(x, y, z - e.zStart, t);
  }

  potentialAt3D(x, y, z, t) {
    const e = this.elementAt(z);
    if (!e) return 0;
    return e.potentialAt(x, y, z - e.zStart, t);
  }

  strikes(x, y, z) {
    const e = this.elementAt(z);
    if (!e) return false;
    return e.strikes(x, y, z - e.zStart);
  }

  get zRange() {
    return [0, this.length];
  }

  /** Widest the column gets, for framing the view. */
  get radiusLimit() {
    return Math.max(1e-3, ...this.elements.map((e) => e.outerRadius));
  }

  /**
   * Finest spatial detail anywhere in the column.
   *
   * The whole line is stepped at the pace of its most demanding element,
   * because an ion that resolves a coarse drift well is not thereby resolving
   * a fine quadrupole.
   */
  get lengthScale() {
    const scales = this.elements
      .map((e) => e.lengthScale)
      .filter((s) => Number.isFinite(s) && s > 0);
    return scales.length ? Math.min(...scales) : 1e-3;
  }

  /** Shortest period of any time-dependent element, or null if all static. */
  get shortestPeriod() {
    const periods = this.elements
      .map((e) => e.shortestPeriod)
      .filter((p) => p != null && p > 0);
    return periods.length ? Math.min(...periods) : null;
  }

  /* ---------------------------------------------------------------- */
  /* diagnostics                                                      */
  /* ---------------------------------------------------------------- */

  /**
   * Problems with the column as assembled, as opposed to with any one
   * element.
   */
  get warnings() {
    const out = [];
    for (const e of this.elements) {
      for (const w of e.warnings ?? []) out.push(`${e.label}: ${w}`);
    }

    // Two live elements with no grounded drift between them: each was solved
    // as though the other were not there.
    for (let i = 1; i < this.elements.length; i++) {
      const a = this.elements[i - 1];
      const b = this.elements[i];
      if (a.type === 'drift' || b.type === 'drift') continue;
      out.push(
        `${a.label} and ${b.label} are adjacent with no drift between them. ` +
          'Each was solved in isolation with grounded end faces, so the field ' +
          'at their junction is not a true solution for the pair. Insert a ' +
          'drift of at least one bore radius.'
      );
    }

    // A bore that steps outward then back in is a real aperture, but a bore
    // that steps DOWN abruptly will scrape the beam.
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

    return out;
  }

  /** Flat list of drawable electrode rectangles in beamline coordinates. */
  outline() {
    const out = [];
    for (const e of this.elements) {
      for (const r of e.rects ?? []) {
        out.push({
          z0: e.zStart + r.z0,
          z1: e.zStart + r.z1,
          r0: r.r0,
          r1: r.r1,
          ghost: !!r.ghost,
          wall: !!r.wall,
          element: e,
        });
      }
    }
    return out;
  }
}
