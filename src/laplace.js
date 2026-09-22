/**
 * Laplace solver - finds the electrostatic potential in the charge-free
 * region between electrodes.
 *
 * In a region with no free charge, Gauss's law reduces to Laplace's equation
 *
 *     div grad phi = 0
 *
 * with phi fixed on every electrode surface (a Dirichlet problem). IonTrace
 * solves it by finite differences with successive over-relaxation, which is
 * the same class of method SIMION uses to refine a potential array.
 *
 * Discretisation
 * --------------
 * Planar (2D Cartesian, coordinates z and y):
 *
 *     d2phi/dz2 + d2phi/dy2 = 0
 *
 * On a uniform grid of step h the five-point stencil gives the update
 *
 *     phi[i][j] = ( phi[i+1][j] + phi[i-1][j]
 *                 + phi[i][j+1] + phi[i][j-1] ) / 4
 *
 * Cylindrical (rotationally symmetric, coordinates z and r):
 *
 *     d2phi/dr2 + (1/r) dphi/dr + d2phi/dz2 = 0
 *
 * Off the axis, at r = j*h with j > 0, central differences give
 *
 *     phi[i][j] = [ phi[i][j+1] (1 + 1/(2j))
 *                 + phi[i][j-1] (1 - 1/(2j))
 *                 + phi[i+1][j] + phi[i-1][j] ] / 4
 *
 * On the axis the (1/r) dphi/dr term is singular and must be handled by its
 * limit. Symmetry forces dphi/dr = 0 at r = 0, so by L'Hopital
 * (1/r) dphi/dr -> d2phi/dr2 and the equation becomes
 *
 *     2 d2phi/dr2 + d2phi/dz2 = 0
 *
 * Mirroring the ghost node phi[i][-1] = phi[i][1] then yields
 *
 *     phi[i][0] = ( 4 phi[i][1] + phi[i+1][0] + phi[i-1][0] ) / 6
 *
 * Both stencils are second-order accurate, so halving h should cut the
 * discretisation error by about four. tests/physics.test.js checks that.
 */

import { CYLINDRICAL, NO_ELECTRODE } from './grid.js';

/**
 * Theoretically optimal over-relaxation factor for a Dirichlet Laplace
 * problem on an nz x nr grid.
 *
 * The Jacobi iteration's spectral radius on such a grid is
 * rho = (cos(pi/nz) + cos(pi/nr)) / 2, and the SOR optimum follows as
 * omega = 2 / (1 + sqrt(1 - rho^2)). Real geometries have interior electrodes
 * that this estimate ignores, so it is an approximation - but a good enough
 * one to turn an O(N^2) relaxation into roughly O(N^1.5).
 */
export function optimalOmega(nz, nr) {
  const rho = (Math.cos(Math.PI / nz) + Math.cos(Math.PI / nr)) / 2;
  return 2 / (1 + Math.sqrt(1 - rho * rho));
}

/**
 * Relax `phi` in place until the largest node change in a sweep falls below
 * `tolerance`, holding every electrode node fixed at its current value.
 *
 * @param {import('./grid.js').PotentialArray} grid
 * @param {Float64Array} phi  Potentials, pre-loaded with electrode values.
 * @param {object} [opts]
 * @param {number} [opts.tolerance]   Max nodal change to declare convergence.
 * @param {number} [opts.maxSweeps]   Hard iteration cap.
 * @param {number} [opts.omega]       Over-relaxation factor; defaults to the
 *                                    theoretical optimum for this grid.
 * @returns {{sweeps: number, change: number, converged: boolean}}
 */
export function relax(grid, phi, opts = {}) {
  const { nz, nr, symmetry, electrodeId } = grid;
  const tolerance = opts.tolerance ?? 1e-9;
  const maxSweeps = opts.maxSweeps ?? 20000;
  const omega = opts.omega ?? optimalOmega(nz, nr);

  if (omega <= 0 || omega >= 2) {
    // SOR diverges outside this interval; catching it here beats silently
    // returning a garbage field to the trajectory integrator.
    throw new Error(`Over-relaxation factor must lie in (0, 2), got ${omega}`);
  }

  const cylindrical = symmetry === CYLINDRICAL;
  let change = Infinity;
  let sweeps = 0;

  while (sweeps < maxSweeps && change > tolerance) {
    change = 0;

    // The axis row is only a free row in cylindrical mode; in planar mode the
    // enclosure has already pinned it, so the generic loop below skips it.
    if (cylindrical) {
      for (let i = 1; i < nz - 1; i++) {
        const k = i; // j = 0, so index == i
        if (electrodeId[k] !== NO_ELECTRODE) continue;
        const target = (4 * phi[nz + i] + phi[k + 1] + phi[k - 1]) / 6;
        const delta = omega * (target - phi[k]);
        phi[k] += delta;
        const mag = Math.abs(delta);
        if (mag > change) change = mag;
      }
    }

    for (let j = 1; j < nr - 1; j++) {
      const row = j * nz;
      // Radial weights are constant along a row, so hoist them out of the
      // inner loop. In planar mode both are exactly 1 and the stencil
      // collapses to the standard five-point form.
      const wUp = cylindrical ? 1 + 1 / (2 * j) : 1;
      const wDown = cylindrical ? 1 - 1 / (2 * j) : 1;

      for (let i = 1; i < nz - 1; i++) {
        const k = row + i;
        if (electrodeId[k] !== NO_ELECTRODE) continue;
        const target =
          (wUp * phi[k + nz] + wDown * phi[k - nz] + phi[k + 1] + phi[k - 1]) / 4;
        const delta = omega * (target - phi[k]);
        phi[k] += delta;
        const mag = Math.abs(delta);
        if (mag > change) change = mag;
      }
    }

    sweeps++;
  }

  return { sweeps, change, converged: change <= tolerance };
}

/**
 * Largest absolute residual of Laplace's equation over the free nodes,
 * expressed as h^2 * div grad phi in volts.
 *
 * This is an independent check on the answer rather than on the iteration: a
 * converged relaxation can still be wrong if a stencil is mis-derived, but it
 * cannot have a small residual if it is. Used by the test suite.
 */
export function maxResidual(grid, phi) {
  const { nz, nr, symmetry, electrodeId } = grid;
  const cylindrical = symmetry === CYLINDRICAL;
  let worst = 0;

  if (cylindrical) {
    for (let i = 1; i < nz - 1; i++) {
      if (electrodeId[i] !== NO_ELECTRODE) continue;
      // 2 d2phi/dr2 + d2phi/dz2, mirrored across the axis.
      const res =
        4 * phi[nz + i] - 6 * phi[i] + phi[i + 1] + phi[i - 1];
      worst = Math.max(worst, Math.abs(res));
    }
  }

  for (let j = 1; j < nr - 1; j++) {
    const row = j * nz;
    const wUp = cylindrical ? 1 + 1 / (2 * j) : 1;
    const wDown = cylindrical ? 1 - 1 / (2 * j) : 1;
    for (let i = 1; i < nz - 1; i++) {
      const k = row + i;
      if (electrodeId[k] !== NO_ELECTRODE) continue;
      const res =
        wUp * phi[k + nz] + wDown * phi[k - nz] + phi[k + 1] + phi[k - 1] - 4 * phi[k];
      worst = Math.max(worst, Math.abs(res));
    }
  }

  return worst;
}

/**
 * Solve the unit basis solutions that make fast adjustment possible.
 *
 * For each electrode in turn, Laplace's equation is solved with that
 * electrode held at 1 V and every other electrode at 0 V. Because Laplace's
 * equation is linear and the boundary conditions are purely Dirichlet, the
 * field for an arbitrary set of applied voltages is the weighted sum
 *
 *     phi(x) = SUM_i  V_i * phi_i(x)
 *
 * of these basis solutions. That is SIMION's "fast adjust": the expensive
 * relaxation is paid once per electrode, after which changing a voltage costs
 * one multiply-add per node instead of a fresh solve.
 *
 * The superposition is only valid because every basis solution satisfies the
 * *same* homogeneous outer boundary condition (the grounded enclosure sits at
 * 0 V in all of them). field.js documents the consequence for the caller.
 *
 * @returns {{basis: Float64Array[], reports: object[]}}
 */
export function solveBasis(grid, opts = {}) {
  const n = grid.nz * grid.nr;
  const basis = [];
  const reports = [];

  for (let e = 0; e < grid.electrodeCount; e++) {
    const phi = new Float64Array(n);
    for (let k = 0; k < n; k++) {
      if (grid.electrodeId[k] === e) phi[k] = 1;
    }
    const report = relax(grid, phi, opts);
    if (!report.converged) {
      console.warn(
        `IonTrace: electrode "${grid.electrodeNames[e]}" basis solution did not ` +
          `converge (change ${report.change.toExponential(2)} after ${report.sweeps} sweeps)`
      );
    }
    basis.push(phi);
    reports.push(report);
  }

  return { basis, reports };
}
