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
 * The Jacobi iteration's spectral radius on a grid with nz by nr TOTAL nodes,
 * hence nz-2 interior ones, is rho = (cos(pi/(nz-1)) + cos(pi/(nr-1))) / 2,
 * and the SOR optimum follows as omega = 2 / (1 + sqrt(1 - rho^2)). The
 * denominators are node spacings, not node counts; using nz directly biases
 * omega upward, which is the worse side to err on because the SOR spectral
 * radius degrades linearly in omega - 1 above the optimum.
 *
 * This remains an estimate. The formula is derived for an empty rectangle,
 * and interior electrodes change the spectrum; measured against a sweep of
 * omega on the shipped einzel geometry it costs roughly 30 % more sweeps than
 * the empirical optimum. That is a speed matter, not an accuracy one - the
 * converged solution is the same - and it still turns an O(N^2) relaxation
 * into roughly O(N^1.5).
 */
export function optimalOmega(nz, nr) {
  const rho =
    (Math.cos(Math.PI / (nz - 1)) + Math.cos(Math.PI / (nr - 1))) / 2;
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
  // The singular on-axis stencil applies only when the grid actually reaches
  // r = 0. A cylindrical band that starts further out - a bender's, which
  // spans a narrow range about its bend radius - has an ordinary wall at
  // j = 0, and applying the axis limit there would be solving a different
  // problem entirely.
  const hasAxis = cylindrical && grid.r0 === 0;
  let change = Infinity;
  let sweeps = 0;

  while (sweeps < maxSweeps && change > tolerance) {
    change = 0;

    if (hasAxis) {
      for (let i = 1; i < nz - 1; i++) {
        const k = i; // j = 0, so index == i
        if (electrodeId[k] !== NO_ELECTRODE) continue;
        const target = (4 * phi[nz + i] + phi[k + 1] + phi[k - 1]) / 6;
        const delta = omega * (target - phi[k]);
        phi[k] += delta;
        // Math.max, not a comparison. `NaN > change` is false, so a
        // comparison silently discards non-finite updates and leaves `change`
        // at zero - the iteration would then report convergence on a grid
        // that had blown up to NaN. Math.max propagates NaN, which falls out
        // of the loop and is reported as not converged.
        change = Math.max(change, Math.abs(delta));
      }
    }

    for (let j = 1; j < nr - 1; j++) {
      const row = j * nz;
      // Radial weights are constant along a row, so hoist them out of the
      // inner loop. In planar mode both are exactly 1 and the stencil
      // collapses to the standard five-point form.
      //
      // The weights depend on the PHYSICAL radius, not the row index, so a
      // band that starts away from the axis gets 1 +/- h/(2r) with r = r0 + jh
      // rather than 1 +/- 1/(2j). They coincide only when r0 is zero.
      const r = grid.rAt(j);
      const wUp = cylindrical ? 1 + grid.step / (2 * r) : 1;
      const wDown = cylindrical ? 1 - grid.step / (2 * r) : 1;

      for (let i = 1; i < nz - 1; i++) {
        const k = row + i;
        if (electrodeId[k] !== NO_ELECTRODE) continue;
        const target =
          (wUp * phi[k + nz] + wDown * phi[k - nz] + phi[k + 1] + phi[k - 1]) / 4;
        const delta = omega * (target - phi[k]);
        phi[k] += delta;
        // Math.max, not a comparison. `NaN > change` is false, so a
        // comparison silently discards non-finite updates and leaves `change`
        // at zero - the iteration would then report convergence on a grid
        // that had blown up to NaN. Math.max propagates NaN, which falls out
        // of the loop and is reported as not converged.
        change = Math.max(change, Math.abs(delta));
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
 * This measures how completely the ITERATION has converged, and nothing more.
 * It deliberately does not claim to validate the stencil: it evaluates the
 * same expressions `relax` iterates, so a mis-derived stencil reproduces its
 * own error here and the residual still falls to round-off. Relaxing with
 * deliberately wrong radial weights of 1 +/- 1/j gives a max error of 9.4e-2
 * against the analytic solution while this function reports 3.4e-14.
 *
 * The stencil itself is validated instead by solving a closed-form harmonic
 * function and comparing - see tests/physics.test.js.
 */
export function maxResidual(grid, phi) {
  const { nz, nr, symmetry, electrodeId } = grid;
  const cylindrical = symmetry === CYLINDRICAL;
  const hasAxis = cylindrical && grid.r0 === 0;
  let worst = 0;

  if (hasAxis) {
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
    const r = grid.rAt(j);
    const wUp = cylindrical ? 1 + grid.step / (2 * r) : 1;
    const wDown = cylindrical ? 1 - grid.step / (2 * r) : 1;
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
 * Verify that every node the relaxation never updates belongs to an electrode.
 *
 * `relax` sweeps only the interior (plus the axis row in cylindrical mode),
 * so the domain rim is frozen at whatever it was initialised to regardless of
 * whether any electrode was painted there. An unpainted rim is therefore held
 * at 0 V invisibly - a grounded box the caller never asked for and cannot
 * change by setting voltages. Catching it here turns a silently wrong field
 * into an error at build time.
 */
export function assertFrozenNodesAreOwned(grid) {
  const { nz, nr, symmetry, electrodeId } = grid;
  const axisIsFree = symmetry === CYLINDRICAL;

  for (let j = 0; j < nr; j++) {
    for (let i = 0; i < nz; i++) {
      const rim =
        i === 0 ||
        i === nz - 1 ||
        j === nr - 1 ||
        (j === 0 && !axisIsFree);
      if (!rim) continue;
      if (electrodeId[j * nz + i] === NO_ELECTRODE) {
        throw new Error(
          `Grid node (i=${i}, j=${j}) lies on the domain rim but belongs to no ` +
            'electrode. The relaxation never updates rim nodes, so it would be ' +
            'frozen at 0 V and impose an invisible grounded boundary. Paint an ' +
            'enclosure electrode over the whole rim (see paintEnclosure).'
        );
      }
    }
  }
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
 * What makes the superposition valid is NOT that some electrode sits at 0 V.
 * In this codebase the enclosure is itself a registered electrode, so its own
 * basis solution has the enclosure at 1 V and nothing is grounded at all -
 * and superposition still holds exactly. The real conditions are:
 *
 *   1. the set of Dirichlet-frozen nodes is identical in every basis solve
 *      (the geometry must not change between them), and
 *   2. those frozen nodes are PARTITIONED by the electrodes - every frozen
 *      node is owned by exactly one electrode.
 *
 * If a frozen node were owned by no electrode, `solveBasis` would hold it at
 * zero in every basis solution and no combination could ever reproduce its
 * intended potential. The assertion below enforces condition 2, because the
 * relaxation silently freezes the domain rim whether or not anything was
 * painted there - so an unpainted rim would otherwise impose an invisible
 * grounded box that no voltage setting could override.
 *
 * @returns {{basis: Float64Array[], reports: object[]}}
 */
export function solveBasis(grid, opts = {}) {
  const n = grid.nz * grid.nr;
  const basis = [];
  const reports = [];

  assertFrozenNodesAreOwned(grid);

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
