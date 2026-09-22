/**
 * Field - turns solved basis potentials into the electric field an ion feels.
 *
 * Two jobs:
 *
 *  1. Fast adjust. Combine the per-electrode unit solutions from
 *     laplace.solveBasis() into a potential map for the applied voltages,
 *     phi = SUM_i V_i phi_i. Changing a voltage costs one pass over the grid,
 *     not a fresh relaxation.
 *
 *  2. Differentiate and interpolate. The force on an ion is qE with
 *     E = -grad phi. The gradient is taken by central differences at the grid
 *     nodes, then bilinearly interpolated to the ion's actual position.
 *
 * Why differentiate first and interpolate second
 * ----------------------------------------------
 * The alternative - interpolate phi, then differentiate the interpolant - is
 * cheaper but worse. Bilinear phi has a piecewise-constant derivative normal
 * to each cell edge, so E would jump discontinuously as an ion crossed a cell
 * boundary, injecting spurious impulses that a trajectory integrator happily
 * accumulates into drifting energy. Differencing at the nodes and then
 * interpolating gives an E that is continuous across cell boundaries, which
 * is what the energy-conservation test in tests/physics.test.js relies on.
 */

import { CYLINDRICAL, PLANAR, NO_ELECTRODE } from './grid.js';

export class Field {
  /**
   * @param {import('./grid.js').PotentialArray} grid
   * @param {Float64Array[]} basis  One unit solution per electrode, in
   *                                electrode-id order.
   */
  constructor(grid, basis) {
    if (basis.length !== grid.electrodeCount) {
      throw new Error(
        `Expected ${grid.electrodeCount} basis solutions, got ${basis.length}`
      );
    }
    this.grid = grid;
    this.basis = basis;

    const n = grid.nz * grid.nr;
    this.phi = new Float64Array(n);
    this.Ez = new Float64Array(n);
    this.Er = new Float64Array(n);
    this.voltages = new Float64Array(grid.electrodeCount);

    this.setVoltages(this.voltages);
  }

  /**
   * Apply electrode voltages (volts) and rebuild the potential and field maps.
   * Accepts an array in electrode-id order, or an object keyed by electrode
   * name. Any electrode omitted from an object form is taken as 0 V - stating
   * it explicitly is better practice, since a floating electrode is a real
   * thing this model cannot represent.
   */
  setVoltages(voltages) {
    const { grid } = this;
    const v = this.voltages;

    if (Array.isArray(voltages) || ArrayBuffer.isView(voltages)) {
      if (voltages.length !== grid.electrodeCount) {
        throw new Error(
          `Expected ${grid.electrodeCount} voltages, got ${voltages.length}`
        );
      }
      for (let e = 0; e < v.length; e++) v[e] = voltages[e];
    } else {
      for (const name of Object.keys(voltages)) {
        const id = grid.electrodeNames.indexOf(name);
        if (id === -1) throw new Error(`Unknown electrode "${name}"`);
      }
      for (let e = 0; e < v.length; e++) {
        v[e] = voltages[grid.electrodeNames[e]] ?? 0;
      }
    }

    const { phi, basis } = this;
    phi.fill(0);
    for (let e = 0; e < v.length; e++) {
      const Ve = v[e];
      if (Ve === 0) continue; // a zeroed electrode contributes nothing
      const b = basis[e];
      for (let k = 0; k < phi.length; k++) phi[k] += Ve * b[k];
    }

    this.#computeGradient();
    return this;
  }

  /** Voltage currently applied to a named electrode. */
  voltageOf(name) {
    const id = this.grid.electrodeNames.indexOf(name);
    if (id === -1) throw new Error(`Unknown electrode "${name}"`);
    return this.voltages[id];
  }

  /**
   * Nodal E = -grad phi by central differences, one-sided at the rim.
   *
   * On the cylindrical axis Er is set to exactly zero rather than differenced.
   * That is not a numerical convenience: rotational symmetry means a non-zero
   * radial field at r = 0 would have no direction to point in, so Er(0, z) = 0
   * is exact, and imposing it keeps ions launched on-axis from acquiring a
   * spurious radial kick from round-off.
   *
   * A central difference must not straddle a conductor surface. At a node
   * lying ON an electrode, the neighbour on the metal side sits inside the
   * conductor, where the potential is constant at the electrode value. The
   * central difference then returns
   *
   *     (phi_vacuum - V) / 2h   instead of   (phi_vacuum - V) / h
   *
   * exactly half the true surface field - and being a factor, not a
   * truncation term, it does not shrink as the grid is refined. Bilinear
   * interpolation spreads that halved value a full cell into the vacuum, so
   * an ion grazing an aperture feels a systematically weak radial field. Such
   * nodes therefore use the one-sided difference on the vacuum side, which is
   * the correct surface derivative. Nodes buried inside metal keep the
   * central form and give E = 0, as a conductor's interior should.
   */
  #computeGradient() {
    const { nz, nr, step, symmetry, electrodeId } = this.grid;
    const { phi, Ez, Er } = this;
    const inv2h = 1 / (2 * step);
    const invh = 1 / step;
    const cylindrical = symmetry === CYLINDRICAL;
    const metal = (k) => electrodeId[k] !== NO_ELECTRODE;

    for (let j = 0; j < nr; j++) {
      const row = j * nz;
      for (let i = 0; i < nz; i++) {
        const k = row + i;
        const onMetal = metal(k);

        // dphi/dz
        if (i === 0) Ez[k] = -(phi[k + 1] - phi[k]) * invh;
        else if (i === nz - 1) Ez[k] = -(phi[k] - phi[k - 1]) * invh;
        else if (onMetal && metal(k - 1) && !metal(k + 1))
          Ez[k] = -(phi[k + 1] - phi[k]) * invh;
        else if (onMetal && metal(k + 1) && !metal(k - 1))
          Ez[k] = -(phi[k] - phi[k - 1]) * invh;
        else Ez[k] = -(phi[k + 1] - phi[k - 1]) * inv2h;

        // dphi/dr
        if (j === 0) {
          Er[k] = cylindrical ? 0 : -(phi[k + nz] - phi[k]) * invh;
        } else if (j === nr - 1) {
          Er[k] = -(phi[k] - phi[k - nz]) * invh;
        } else if (onMetal && metal(k - nz) && !metal(k + nz)) {
          Er[k] = -(phi[k + nz] - phi[k]) * invh;
        } else if (onMetal && metal(k + nz) && !metal(k - nz)) {
          Er[k] = -(phi[k] - phi[k - nz]) * invh;
        } else {
          Er[k] = -(phi[k + nz] - phi[k - nz]) * inv2h;
        }
      }
    }
  }

  /** True if (z, r) in metres lies inside the solved domain. */
  contains(z, r) {
    const { grid } = this;
    return (
      z >= grid.z0 &&
      z <= grid.z0 + grid.zLength &&
      r >= 0 &&
      r <= grid.rLength
    );
  }

  /**
   * Locate (z, r) within the grid and return the enclosing cell corner plus
   * the fractional offsets used by every bilinear lookup here.
   */
  #locate(z, r) {
    const { grid } = this;
    const gz = (z - grid.z0) / grid.step;
    const gr = r / grid.step;

    // Clamp so a particle exactly on the far edge still lands in a valid cell.
    let i = Math.floor(gz);
    let j = Math.floor(gr);
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i > grid.nz - 2) i = grid.nz - 2;
    if (j > grid.nr - 2) j = grid.nr - 2;

    // Clamp the fractional offsets too. Clamping only the cell index leaves
    // fz or fr outside [0, 1] for a point beyond the domain, and the bilinear
    // form then EXTRAPOLATES rather than saturating: the grounded 16 mm
    // housing of the default lens reports 807 V at r = 66 mm, growing without
    // limit. RK4's intermediate stages can land out there even when the step
    // endpoint does not, so the fabricated field would feed back into a
    // trajectory that never visibly left the grid.
    const fz = Math.min(1, Math.max(0, gz - i));
    const fr = Math.min(1, Math.max(0, gr - j));

    return { i, j, fz, fr };
  }

  /** Bilinear sample of `arr` at (z, r). */
  #sample(arr, z, r) {
    const { nz } = this.grid;
    const { i, j, fz, fr } = this.#locate(z, r);
    const k = j * nz + i;
    const a = arr[k];
    const b = arr[k + 1];
    const c = arr[k + nz];
    const d = arr[k + nz + 1];
    return (
      a * (1 - fz) * (1 - fr) +
      b * fz * (1 - fr) +
      c * (1 - fz) * fr +
      d * fz * fr
    );
  }

  /** Electrostatic potential in volts at (z, r) in metres. */
  potentialAt(z, r) {
    return this.#sample(this.phi, z, r);
  }

  /** Field in V/m at (z, r) in metres, as axial and radial components. */
  fieldAt(z, r) {
    return {
      Ez: this.#sample(this.Ez, z, r),
      Er: this.#sample(this.Er, z, r),
    };
  }

  /**
   * Field in the meridional plane, in Cartesian components.
   *
   * Trajectories are integrated in a plane containing the axis, with x the
   * signed transverse coordinate so an ion can cross the axis instead of
   * bouncing off r = 0. The radial field is resolved onto x by
   * Ex = Er * (x / |x|), which is exact for motion with zero angular
   * momentum about the axis - the case IonTrace currently supports, since no
   * azimuthal field component exists to create any.
   *
   * In planar mode there is no axis to fold about: the transverse coordinate
   * is an ordinary Cartesian y running from 0 to the domain height, so x is
   * used as given and the transverse field passes straight through.
   */
  fieldAtCartesian(x, z) {
    if (this.grid.symmetry === PLANAR) {
      const { Ez, Er } = this.fieldAt(z, x);
      return { Ex: Er, Ez };
    }
    const r = Math.abs(x);
    const { Ez, Er } = this.fieldAt(z, r);
    // At x = 0, Er is already zero by symmetry, so the sign is immaterial.
    const Ex = r > 0 ? Er * (x / r) : 0;
    return { Ex, Ez };
  }

  /** Potential in the meridional plane at signed transverse position x. */
  potentialAtCartesian(x, z) {
    return this.potentialAt(z, this.grid.symmetry === PLANAR ? x : Math.abs(x));
  }
}
