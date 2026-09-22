/**
 * Space charge - the beam's repulsion of itself.
 *
 * Why this is not pairwise Coulomb
 * --------------------------------
 * IonTrace is axisymmetric. A trajectory drawn in the meridional plane is not
 * a single ion: it is the cross-section of a RING of charge at radius r,
 * carrying its share of the beam current all the way around the azimuth.
 * Computing q1 q2 / 4 pi eps0 d^2 between two such rays would be the force
 * between two point charges, which is not the force between two rings, and
 * would also break the rotational symmetry the field solve depends on.
 *
 * The correct treatment for a long axisymmetric beam is Gauss's law. Take a
 * cylinder of radius r and length L about the axis. By symmetry E is purely
 * radial on its curved surface and the flat ends contribute nothing, so
 *
 *     E_r(r) . 2 pi r L = Q_enclosed / eps0 = lambda_enc(r) . L / eps0
 *
 *     E_r(r) = lambda_enc(r) / (2 pi eps0 r)
 *
 * where lambda_enc(r) is the line charge density enclosed within radius r.
 * Only the charge inside r matters; a uniform shell outside contributes
 * exactly nothing, which is the cylindrical shell theorem.
 *
 * For a beam of current I moving at axial speed v_z, the line charge density
 * is
 *
 *     lambda = I / v_z
 *
 * so a slow beam is a dense one. This is why space charge bites hardest where
 * the optics decelerate the beam, and it is the reason the model is driven by
 * a beam CURRENT - a real instrument parameter - rather than by the number of
 * rays drawn, which is a display setting and has no physics in it.
 *
 * Ray weighting
 * -------------
 * Each ray carries a fixed share w_i of the beam current, fixed at launch from
 * the assumed initial current density and conserved thereafter. The enclosed
 * fraction for ray i is then the sum of the shares of every ray currently
 * inside it, plus half its own - a ring does not exert a net force on itself,
 * and taking half places the ray at the middle of the annulus it represents.
 *
 * Because rays are re-sorted by radius each step rather than assumed laminar,
 * this stays valid after the beam crosses over at a focus, where the ordering
 * genuinely changes.
 *
 * What is neglected
 * -----------------
 *  - The beam's own magnetic field. Moving charges attract magnetically, and
 *    the self-force is reduced by a factor (1 - beta^2). For keV ions
 *    beta ~ 1e-4, so this is a part in 1e8 and is genuinely ignorable.
 *  - Image charges in the surrounding electrodes, which partially shield the
 *    space charge. This matters for a beam that fills the bore, and makes
 *    IonTrace's space-charge defocusing an OVERESTIMATE in that regime.
 *  - Axial (longitudinal) space charge. The long-beam approximation assumes
 *    the beam is much longer than it is wide, which is true for a continuous
 *    beam and false for a short bunch. Bunches are not modelled.
 *  - Poisson is not re-solved with the charge density; the beam field is
 *    added analytically to the solved electrode field. That superposition is
 *    exact for the free-space part and is what drops the image charges above.
 */

import { VACUUM_PERMITTIVITY } from './constants.js';

/** Coulomb constant 1 / (4 pi eps0), in V m / C. */
const COULOMB_K = 1 / (4 * Math.PI * VACUUM_PERMITTIVITY);

/**
 * Pairwise Coulomb field at each particle from every other particle.
 *
 * This is the DISCRETE alternative to the ring model below: each simulated
 * particle is a point charge and feels every other one directly,
 *
 *     E_i = (1 / 4 pi eps0) SUM_{j != i}  w q_j (r_i - r_j) / |r_i - r_j|^3
 *
 * Use it for a bunch, a cloud, or any case where the ions are genuinely
 * countable. Use the ring model for a continuous beam, where a trajectory
 * stands for a whole ring of charge rather than one ion.
 *
 * Macro-weighting
 * ---------------
 * `weight` is how many real ions each simulated particle stands for. With
 * weight 1 the calculation is literally N ions and the repulsion is, for any
 * N you can draw on screen, far too small to see - nine elementary charges
 * spread over millimetres produce a field of order 1e-3 V/m against electrode
 * fields of 1e4. That is the correct answer, not a defect, and it is why
 * particle codes weight.
 *
 * A macroparticle of weight w has charge wq and mass wm, so q/m is unchanged
 * and the electrode force is unaffected; only the mutual force scales, and it
 * scales linearly in w. The w real ions inside one macroparticle do not repel
 * each other in this model, which is the standard particle-in-cell
 * approximation.
 *
 * Softening
 * ---------
 * The 1/r^2 force diverges as two particles approach, and with a finite time
 * step a close pass would fling them apart with energy that came from nowhere.
 * `softening` replaces |d|^3 with (|d|^2 + eps^2)^{3/2} - Plummer softening -
 * which bounds the force at short range. It is a real approximation: close
 * encounters, and therefore collisional relaxation of the bunch, are
 * suppressed. It does not affect the long-range behaviour that drives beam
 * expansion.
 *
 * Cost is O(N^2). For the particle counts a browser will draw that is nothing,
 * but it is the reason this does not scale to a real PIC simulation.
 *
 * @param {number[]} xs        Transverse positions, metres.
 * @param {number[]} zs        Axial positions, metres.
 * @param {number[]} charges   Charge of each particle, coulombs.
 * @param {number} weight      Real ions represented per particle.
 * @param {number} softening   Plummer softening length, metres.
 * @returns {{Ex: number[], Ez: number[]}} Field at each particle, V/m.
 */
export function coulombField(xs, zs, charges, weight = 1, softening = 0) {
  const n = xs.length;
  const Ex = new Array(n).fill(0);
  const Ez = new Array(n).fill(0);
  if (n < 2 || weight === 0) return { Ex, Ez };

  const eps2 = softening * softening;

  for (let i = 0; i < n; i++) {
    let ex = 0;
    let ez = 0;
    for (let j = 0; j < n; j++) {
      if (j === i) continue;
      const dx = xs[i] - xs[j];
      const dz = zs[i] - zs[j];
      const d2 = dx * dx + dz * dz + eps2;
      if (d2 <= 0) continue;
      // 1 / d^3, via d^2 and its square root.
      const inv = 1 / (d2 * Math.sqrt(d2));
      const c = COULOMB_K * weight * charges[j] * inv;
      ex += c * dx;
      ez += c * dz;
    }
    Ex[i] = ex;
    Ez[i] = ez;
  }

  return { Ex, Ez };
}

/** 1 / (2 pi eps0), the recurring factor in the cylindrical Gauss result. */
const INV_TWO_PI_EPS0 = 1 / (2 * Math.PI * VACUUM_PERMITTIVITY);

/**
 * Fraction of the beam current carried by each ray, from its launch radius.
 *
 * A ray at radius r represents the annulus halfway to its neighbours on each
 * side, so for a beam of uniform current density its share is proportional to
 * that annulus's area. The shares sum to 1.
 *
 * The uniform-density assumption is the only profile IonTrace offers. A real
 * source is rarely uniform, and a peaked profile concentrates more current at
 * small radius and defocuses more strongly.
 *
 * @param {number[]} radii Launch radii in metres, any order.
 * @returns {number[]} Shares in the same order, summing to 1.
 */
export function currentShares(radii) {
  const n = radii.length;
  if (n === 0) return [];
  if (n === 1) return [1];

  // Rays at the SAME radius are the same ring, not separate ones.
  //
  // parallelBeam launches rays at signed offsets, so a beam spanning -R..+R
  // produces each non-zero radius twice - the two halves of one ring seen in
  // the meridional plane. Treating them as two distinct rings would give them
  // different annuli and count the ring's current twice. Radii are grouped
  // first, each group gets one annulus, and its share is split evenly among
  // the rays that make it up, so the ring's total current is counted once.
  const groups = groupByRadius(radii);
  const outerEdge = groups[groups.length - 1].radius;

  const shares = new Array(n).fill(0);
  let total = 0;

  for (let k = 0; k < groups.length; k++) {
    const r = groups[k].radius;
    const inner = k === 0 ? 0 : (groups[k - 1].radius + r) / 2;
    const outer = k === groups.length - 1 ? outerEdge : (r + groups[k + 1].radius) / 2;
    // Annulus area, up to the common factor of pi.
    const area = Math.max(0, outer * outer - inner * inner);
    const per = area / groups[k].members.length;
    for (const i of groups[k].members) shares[i] = per;
    total += area;
  }

  // A beam launched entirely on the axis has no area anywhere; share it out
  // evenly rather than dividing by zero.
  if (total <= 0) return new Array(n).fill(1 / n);
  return shares.map((s) => s / total);
}

/**
 * Indices grouped by exactly equal radius, in increasing radius order.
 *
 * Exact equality is the right test rather than a tolerance. The case that
 * matters is a ray and its mirror image, and the arithmetic makes those
 * exactly antisymmetric: the radial field is evaluated at |x|, so the
 * acceleration at +x is the exact negation of that at -x, and negation is
 * exact in floating point. Two rays that are merely close in radius really
 * are different rings and should be treated as such.
 */
function groupByRadius(radii) {
  const order = radii.map((_, i) => i).sort((a, b) => radii[a] - radii[b]);
  const groups = [];
  for (const i of order) {
    const last = groups[groups.length - 1];
    if (last && radii[i] === last.radius) last.members.push(i);
    else groups.push({ radius: radii[i], members: [i] });
  }
  return groups;
}

/**
 * Radial space-charge field at each ray, in V/m.
 *
 * Positive means outward. For a positive beam this is outward (defocusing);
 * for a negative beam `lineChargeDensity` is negative and the field points
 * inward, but the force qE is still outward. Sign handling therefore needs no
 * special case.
 *
 * @param {number[]} radii              Current radius of each ray, metres.
 * @param {number[]} shares             Current share of each ray (sums to 1).
 * @param {number} lineChargeDensity    lambda = I / v_z, in C/m (signed).
 * @param {number} [softening]          Radius below which the field is taken
 *                                      as zero. On the axis E_r is exactly
 *                                      zero by symmetry; the floor also keeps
 *                                      the 1/r from exploding for a ray that
 *                                      is merely very close to it.
 * @returns {number[]} E_r per ray, in the order given.
 */
export function spaceChargeField(radii, shares, lineChargeDensity, softening = 0) {
  const n = radii.length;
  const out = new Array(n).fill(0);
  if (n === 0 || lineChargeDensity === 0) return out;

  // Rays at equal radius are one ring and must feel an identical field, so
  // the enclosed charge is accumulated per GROUP rather than per ray. Walking
  // them individually would hand the second of a mirrored pair the first's
  // charge as though it were enclosed, which both breaks mirror symmetry and
  // invents a force between the two halves of a single ring.
  const groups = groupByRadius(radii);

  let enclosed = 0;
  for (const group of groups) {
    let groupShare = 0;
    for (const i of group.members) groupShare += shares[i];

    // Half of its own share: a ring exerts no net force on itself, and the
    // ray sits in the middle of the annulus it stands for.
    const fraction = enclosed + groupShare / 2;
    enclosed += groupShare;

    const r = group.radius;
    if (r > softening && r > 0) {
      const E = (lineChargeDensity * fraction * INV_TWO_PI_EPS0) / r;
      for (const i of group.members) out[i] = E;
    }
  }

  return out;
}

/**
 * Line charge density of a beam of current `currentA` travelling at `speedMs`.
 *
 * lambda = I / v. A beam that is brought to rest has, in this description,
 * infinite charge density - physically the charge piles up, and in a real
 * device that is precisely where a virtual cathode forms and the model stops
 * being valid. `minSpeed` caps it rather than returning Infinity, and callers
 * should treat a capped result as "outside the model's range" rather than as
 * an answer.
 */
export function lineChargeDensity(currentA, speedMs, minSpeed = 1) {
  const v = Math.max(Math.abs(speedMs), minSpeed);
  return currentA / v;
}
