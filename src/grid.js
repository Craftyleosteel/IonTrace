/**
 * PotentialArray - the discretised domain that electrodes are painted onto.
 *
 * This mirrors SIMION's potential array: a regular grid of nodes, each node
 * either an *electrode* node (a fixed-potential Dirichlet boundary, owned by
 * one electrode) or a *free* node whose potential the Laplace solver must
 * find.
 *
 * Coordinate convention
 * ---------------------
 * Index (i, j) with i along z (the optic axis) and j along the transverse
 * coordinate. Node (i, j) sits at the physical point
 *
 *     z = z0 + i * h        i in [0, nz)
 *     r = j * h             j in [0, nr)
 *
 * where h is the grid step in metres. Note r starts at exactly 0, so j = 0 is
 * the symmetry axis in cylindrical mode. Data is stored row-major in i, that
 * is index = j * nz + i, so marching along z is contiguous.
 *
 * Symmetry
 * --------
 * 'cylindrical' treats the transverse coordinate as a radius with rotational
 * symmetry about r = 0. 'planar' treats it as a Cartesian y with no symmetry
 * assumption. The two differ only in the Laplace stencil (see laplace.js);
 * everything in this file is shared.
 */

export const CYLINDRICAL = 'cylindrical';
export const PLANAR = 'planar';

/** Sentinel stored in `electrodeId` for a node that is not an electrode. */
export const NO_ELECTRODE = -1;

export class PotentialArray {
  /**
   * @param {object} opts
   * @param {number} opts.nz     Node count along z (>= 3).
   * @param {number} opts.nr     Node count along r (>= 3).
   * @param {number} opts.step   Grid step h in metres (uniform in both axes).
   * @param {string} [opts.symmetry] CYLINDRICAL (default) or PLANAR.
   * @param {number} [opts.z0]   Physical z of node i = 0, in metres.
   * @param {number} [opts.r0]   Physical transverse coordinate of node j = 0.
   *                             Must be 0 in cylindrical mode, where j = 0 is
   *                             the symmetry axis by definition. In planar
   *                             mode it may be negative, which is what lets a
   *                             transverse (x, y) plane straddle the origin -
   *                             needed for a quadrupole, whose four rods sit
   *                             on both sides of it.
   */
  constructor({ nz, nr, step, symmetry = CYLINDRICAL, z0 = 0, r0 = 0 }) {
    if (!Number.isInteger(nz) || !Number.isInteger(nr)) {
      throw new Error('Grid dimensions must be integers');
    }
    if (nz < 3 || nr < 3) {
      throw new Error('Grid must be at least 3x3 for a central-difference stencil');
    }
    if (!(step > 0)) throw new Error('Grid step must be positive');
    if (symmetry !== CYLINDRICAL && symmetry !== PLANAR) {
      throw new Error(`Unknown symmetry: ${symmetry}`);
    }

    if (symmetry === CYLINDRICAL && r0 < 0) {
      throw new Error(
        'Cylindrical symmetry needs r0 >= 0; a radius cannot be negative'
      );
    }

    this.nz = nz;
    this.nr = nr;
    this.step = step;
    this.symmetry = symmetry;
    this.z0 = z0;
    this.r0 = r0;

    /** Per-node owning electrode index, or NO_ELECTRODE for a free node. */
    this.electrodeId = new Int32Array(nz * nr).fill(NO_ELECTRODE);

    /** Registered electrode names, indexed by electrode id. */
    this.electrodeNames = [];

    /**
     * Which domain faces an ion may fly out through.
     *
     * Closing the domain with a conducting box is what makes the Dirichlet
     * problem well posed, so the end faces must carry a fixed potential and
     * therefore must be electrode nodes. But they are a numerical device, not
     * hardware: in a real beamline the ion enters and leaves through
     * apertures at those planes. Marking them open keeps the field solve
     * correct while letting the trajectory code report such an ion as having
     * left the modelled region rather than as having struck an electrode.
     *
     * The outer radial wall is deliberately not listed. In the einzel
     * geometry that surface is the grounded housing, which is real metal an
     * ion genuinely can hit.
     *
     * This default assumes the grid's z axis IS the beam axis, which is true
     * of every element solved in the r-z plane. Elements that solve in a plane
     * of their own - the mass filter, in x and y, and the quadrupole
     * deflector, in the bend plane - must close both faces, because for them
     * the z ends are ordinary walls. They do that explicitly rather than by
     * omission.
     */
    this.openFaces = { zMin: true, zMax: true };
  }

  /** Flat array index of node (i, j). */
  idx(i, j) {
    return j * this.nz + i;
  }

  /** Physical z of column i, metres. */
  zAt(i) {
    return this.z0 + i * this.step;
  }

  /** Physical r (or y) of row j, metres. */
  rAt(j) {
    return this.r0 + j * this.step;
  }

  /** Total extent along z, metres. */
  get zLength() {
    return (this.nz - 1) * this.step;
  }

  /** Total extent along r, metres. */
  get rLength() {
    return (this.nr - 1) * this.step;
  }

  /** Lowest transverse coordinate the domain represents, metres. */
  get rMin() {
    return this.r0;
  }

  /** Highest transverse coordinate the domain represents, metres. */
  get rMax() {
    return this.r0 + this.rLength;
  }

  /** True if node (i, j) is a fixed-potential electrode node. */
  isElectrode(i, j) {
    return this.electrodeId[this.idx(i, j)] !== NO_ELECTRODE;
  }

  /**
   * Inclusive range test for geometry painting, tolerant of round-off.
   *
   * Use this rather than bare `>=` and `<=` when deciding whether a node lies
   * inside a piece of metal. Node coordinates are accumulated as
   * `z0 + i * step`, while the geometry is written as `mm * 1e-3`; the two
   * routes to what should be the same number differ in the last bit, and
   * which way they differ depends on where the element happens to sit in the
   * grid.
   *
   * Without a tolerance the consequence is not subtle. An electrode whose
   * edge falls exactly on a node is painted one step shorter or longer purely
   * because of its absolute position, so the *same* element placed at two
   * different points in a beamline solves to two different fields. Measured
   * on a 15 mm centre electrode at 0.4 mm resolution, that single node moved
   * the on-axis potential by 7.8 V out of 300 - because it sits where the
   * potential is changing at 20 V/mm.
   *
   * The tolerance is a millionth of a grid step: far too small to capture a
   * node that is genuinely outside, and far larger than the round-off it
   * exists to absorb. Ties resolve inclusively, and consistently.
   */
  spans(value, lo, hi) {
    const tol = this.step * 1e-6;
    return value >= lo - tol && value <= hi + tol;
  }

  /**
   * True if (i, j) lies on a physical rim of the domain.
   *
   * In cylindrical mode the axis row j = 0 is a symmetry line, not a rim, so
   * it is excluded except where the two end faces cross it. Reporting the
   * axis as boundary would contradict paintEnclosure, which deliberately
   * leaves it free so the solver can find the on-axis potential.
   */
  isBoundary(i, j) {
    if (i === 0 || i === this.nz - 1 || j === this.nr - 1) return true;
    return j === 0 && this.symmetry === PLANAR;
  }

  /**
   * Register a named electrode and return its integer id. Ids are assigned in
   * registration order and are the index used by the fast-adjust
   * superposition in field.js, so registration order is part of the contract.
   */
  addElectrode(name) {
    this.electrodeNames.push(name);
    return this.electrodeNames.length - 1;
  }

  get electrodeCount() {
    return this.electrodeNames.length;
  }

  /**
   * Mark every node satisfying `predicate(z, r, i, j)` as belonging to
   * electrode `id`. Coordinates passed to the predicate are in metres, so
   * geometry can be written in physical terms rather than node counts.
   *
   * Later calls overwrite earlier ones, so electrodes painted in sequence
   * behave like layers. Overlap is a geometry error the caller should avoid,
   * but it is not fatal - the last writer owns the node.
   */
  paint(id, predicate) {
    if (id < 0 || id >= this.electrodeNames.length) {
      throw new Error(`Unknown electrode id ${id}`);
    }
    let painted = 0;
    for (let j = 0; j < this.nr; j++) {
      const r = this.rAt(j);
      for (let i = 0; i < this.nz; i++) {
        if (predicate(this.zAt(i), r, i, j)) {
          this.electrodeId[this.idx(i, j)] = id;
          painted++;
        }
      }
    }
    return painted;
  }

  /**
   * Paint the outer rim of the domain as a single electrode. Used to impose
   * the grounded enclosure that makes the Dirichlet problem well posed.
   *
   * In cylindrical mode the axis row j = 0 is deliberately left free along the
   * interior: r = 0 is a symmetry line, not a physical wall, and clamping it
   * would wrongly force the on-axis potential to the enclosure value. The two
   * end faces i = 0 and i = nz-1 are painted over their full height, axis
   * included, because those really are walls and they close the domain.
   */
  paintEnclosure(id) {
    return this.paint(id, (z, r, i, j) => {
      if (i === 0 || i === this.nz - 1) return true; // end walls, full height
      if (j === this.nr - 1) return true; // outer wall
      // Row j = 0 is only a symmetry axis when the grid actually reaches
      // r = 0. A cylindrical band that starts further out - a bender's, which
      // spans a narrow range about its bend radius - has an ordinary wall
      // there, and leaving it free would leave the problem unbounded.
      if (j === 0) return this.symmetry === PLANAR || this.includesAxis === false;
      return false;
    });
  }

  /** True when row j = 0 really is the symmetry axis rather than a wall. */
  get includesAxis() {
    return this.symmetry === CYLINDRICAL && this.r0 === 0;
  }

  /** Count of nodes owned by each electrode, for diagnostics. */
  electrodeNodeCounts() {
    const counts = new Array(this.electrodeNames.length).fill(0);
    for (let k = 0; k < this.electrodeId.length; k++) {
      const id = this.electrodeId[k];
      if (id !== NO_ELECTRODE) counts[id]++;
    }
    return counts;
  }
}
