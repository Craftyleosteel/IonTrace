/**
 * Reduced-Hessian refinement, and the null space of a beamline.
 *
 * What a Hessian method can and cannot do here
 * --------------------------------------------
 * Transmission is a COUNT. It is flat almost everywhere and jumps by one when
 * an ion clears an aperture it previously hit. It has no derivative, let alone
 * a second one, so no Newton method of any kind applies to it - which is why
 * the search in optimize.js is a bracketed scan and why that is not a
 * shortcoming to be fixed.
 *
 * But once the beam is through, a second question starts: of the many settings
 * that transmit everything, which delivers the best beam? THAT is smooth. The
 * spot size at the target moves continuously with the voltages as long as no
 * ion is near an edge, so it has a gradient and a Hessian, and second-order
 * information is exactly what coordinate descent lacks.
 *
 * So the two stages answer different questions, and neither replaces the
 * other: the scan acquires transmission, this refines quality.
 *
 * Why the null space matters
 * --------------------------
 * A beamline's knobs are not independent. A lens and the deflector behind it
 * trade off: strengthen one, weaken the other, and the beam at the target
 * barely notices. Along that combination the merit surface is a flat-bottomed
 * valley, and the Hessian has a near-zero eigenvalue pointing along it.
 *
 * That causes two problems and reveals one fact.
 *
 *   - Coordinate descent crawls along a valley that lies at an angle to its
 *     axes, taking many tiny steps where one diagonal step would do.
 *   - A plain Newton step DIVIDES by the curvature, so a near-zero eigenvalue
 *     produces an enormous, meaningless jump along the direction that matters
 *     least. Finite-difference noise in that eigenvalue is amplified without
 *     limit.
 *
 * The fix is the same in both cases: eigendecompose the Hessian, take the
 * Newton step only in the subspace where the curvature is genuinely positive,
 * and leave the rest alone. The directions left alone are the numerical null
 * space, and they are worth reporting rather than discarding - each one is a
 * combination of knobs the beam does not care about, which is a real and
 * useful statement about the instrument.
 *
 * Scaling
 * -------
 * Everything here is done in scaled coordinates: each knob is divided by a
 * natural size for it, so a deflector at 40 V and a lens at 300 V are
 * comparable. Without that the Hessian's eigenvalues mix units, its
 * eigenvectors are meaningless, and "near zero" has no scale to be near zero
 * against.
 */

import { toLocal, forwardOf } from './frames.js';
import { flyBeam } from './integrator.js';
import { applyKnob, readKnobs, writeKnobs } from './optimize.js';

/**
 * Eigenvalues and eigenvectors of a small symmetric matrix, by cyclic Jacobi.
 *
 * Jacobi rather than anything cleverer because the matrices here are a handful
 * of knobs square, and Jacobi is accurate on exactly that: it is backward
 * stable, needs no shifts or deflation, and gives orthogonal eigenvectors even
 * for repeated eigenvalues - which is the case that matters, since a pair of
 * equally sloppy directions is precisely what a degenerate beamline produces.
 *
 * @param {number[][]} input symmetric, n x n
 * @returns {{values: number[], vectors: number[][]}} vectors[k] is the k-th
 *          eigenvector; values[k] its eigenvalue.
 */
export function symmetricEigen(input) {
  const n = input.length;
  const a = input.map((row) => row.slice());
  const v = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))
  );

  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0;
    for (let p = 0; p < n; p++) {
      for (let q = p + 1; q < n; q++) off += a[p][q] * a[p][q];
    }
    if (off <= 1e-300) break;

    for (let p = 0; p < n; p++) {
      for (let q = p + 1; q < n; q++) {
        if (Math.abs(a[p][q]) < 1e-300) continue;
        const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const sign = theta >= 0 ? 1 : -1;
        const t = sign / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1);
        const s = t * c;

        for (let k = 0; k < n; k++) {
          const akp = a[k][p];
          const akq = a[k][q];
          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
        }
        for (let k = 0; k < n; k++) {
          const apk = a[p][k];
          const aqk = a[q][k];
          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }
        for (let k = 0; k < n; k++) {
          const vkp = v[k][p];
          const vkq = v[k][q];
          v[k][p] = c * vkp - s * vkq;
          v[k][q] = s * vkp + c * vkq;
        }
      }
    }
  }

  const values = a.map((row, i) => row[i]);
  const vectors = values.map((_, j) => v.map((row) => row[j]));

  // Largest curvature first, so "the null space" is the tail of the list.
  const order = values.map((_, i) => i).sort((i, j) => values[j] - values[i]);
  return {
    values: order.map((i) => values[i]),
    vectors: order.map((i) => vectors[i]),
  };
}

/**
 * How good the beam is where it is wanted.
 *
 * Mean squared transverse offset at the target, which is the natural quadratic
 * measure and is smooth in the voltages. Returns `null` if the beam is not
 * fully delivered: a merit that counts a lost ion as "a tighter beam" would
 * reward throwing the edge of the beam away, and every step here is rejected
 * unless transmission holds.
 */
/**
 * Where a trajectory crosses a plane, by interpolation.
 *
 * This is the difference between a merit function that can be differentiated
 * and one that cannot. Taking the last recorded point instead measures the
 * beam wherever the integrator's adaptive step happened to leave off, which
 * moves by a fraction of a step as the voltages change - so a finite
 * difference over a small voltage change reads that jitter rather than the
 * physics. On a first attempt it produced a Hessian whose largest eigenvalue
 * was negative and equal in size to its smallest: a pure numerical saddle,
 * entirely an artefact of where the steps landed.
 */
function crossingOf(points, frame) {
  const o = frame.o;
  const f = forwardOf(frame);
  const along = (p) =>
    (p.x - o[0]) * f[0] + ((p.y ?? 0) - o[1]) * f[1] + (p.z - o[2]) * f[2];

  for (let i = points.length - 1; i > 0; i--) {
    const b = along(points[i]);
    const a = along(points[i - 1]);
    if (b >= 0 && a < 0) {
      const g = a === b ? 0 : -a / (b - a);
      const p = points[i - 1];
      const q = points[i];
      return [
        p.x + g * (q.x - p.x),
        (p.y ?? 0) + g * ((q.y ?? 0) - (p.y ?? 0)),
        p.z + g * (q.z - p.z),
      ];
    }
  }
  return null;
}

export function beamQuality(beamline, makeIons, opts = {}) {
  const { target = beamline.mainEnd, want = null, ...flight } = opts;
  const ions = makeIons();
  const count = ions.length;
  const { tracks } = flyBeam(beamline, ions, { cfl: 0.05, maxSteps: 200000, ...flight });
  const frame = beamline.endFrame(target);

  let arrived = 0;
  let sum = 0;
  for (const t of tracks) {
    if (t.stop !== 'exited') continue;
    const last = t.points[t.points.length - 1];
    if (target) {
      const end = beamline.endNearest(last.x, last.y ?? 0, last.z);
      if (!end || end.element !== target.element || end.port !== target.port) continue;
    }
    const at = crossingOf(t.points, frame);
    if (!at) continue;
    arrived++;
    const [x, y] = toLocal(frame, at);
    sum += x * x + y * y;
  }

  const needed = want ?? count;
  return {
    arrived,
    count,
    merit: arrived >= needed && arrived > 0 ? sum / arrived : null,
    rms: arrived > 0 ? Math.sqrt(sum / arrived) : null,
  };
}

/** A natural size for a knob, used to make the Hessian dimensionless. */
function scaleOf(knob) {
  if (knob.seed) return Math.abs(knob.seed);
  const span = (knob.hi ?? knob.max) - (knob.lo ?? knob.min);
  return Math.max(Math.abs(span) / 6, knob.step);
}

/**
 * Refine a set of voltages by a Newton step in the well-conditioned subspace.
 *
 * @param {object} beamline
 * @param {() => object[]} makeIons
 * @param {object[]} knobs from `tunableKnobs`
 * @param {object} [options]
 * @param {number} [options.iterations=3]
 * @param {number} [options.probe=0.08]  finite-difference step, in knob scales
 * @param {number} [options.floor=0.02]  curvature below this fraction of the
 *        largest counts as null
 * @returns {Promise<object>} what moved, what did not, and why
 */
export async function refineNullSpace(beamline, makeIons, knobs, options = {}) {
  const {
    iterations = 3,
    probe = 0.08,
    floor = 0.02,
    flight = {},
    onProgress,
    shouldStop,
  } = options;

  const n = knobs.length;
  const start = readKnobs(beamline, knobs);
  const scales = knobs.map(scaleOf);
  const base = { merit: null, evaluations: 0, nullSpace: [], steps: 0 };
  if (n < 1) return { ...base, values: start, start, improved: false };

  /** Set the knobs from scaled offsets u, clamped to what each allows. */
  const put = (u) => {
    knobs.forEach((k, i) => {
      const raw = start[i] + u[i] * scales[i];
      const lo = Math.min(k.min, k.max);
      const hi = Math.max(k.min, k.max);
      applyKnob(beamline, k, Math.min(hi, Math.max(lo, raw)));
    });
  };

  let evaluations = 0;
  const first = beamQuality(beamline, makeIons, flight);
  evaluations++;
  const want = first.arrived;
  if (!want) {
    writeKnobs(beamline, knobs, start);
    return { ...base, values: start, start, improved: false, evaluations };
  }

  /** Merit at scaled offsets u; Infinity if the beam is no longer delivered. */
  const at = (u) => {
    put(u);
    const q = beamQuality(beamline, makeIons, { ...flight, want });
    evaluations++;
    return q.merit === null ? Infinity : q.merit;
  };

  let u = new Array(n).fill(0);
  let merit = at(u);
  const meritAtStart = merit;
  let nullSpace = [];
  let steps = 0;

  for (let iter = 0; iter < iterations; iter++) {
    if (shouldStop?.()) break;
    if (!Number.isFinite(merit)) break;

    /*
      Gradient and Hessian by central differences.

      Each knob gets its OWN step, shrunk until the merit is finite on both
      sides of it. That is not fussiness: a deflector transmits over a window
      about a tenth wide, so a probe of eight per cent either way can step
      clean out of it, the merit becomes "the beam is lost", and the whole
      second stage gives up before it has started - which is exactly what it
      did before this. The probe has to fit inside the region where the merit
      is smooth, and how wide that is differs from knob to knob.

      A knob whose probe never fits is left out of this iteration: no
      derivative of it can be measured, so it is reported as a direction
      nothing is known about rather than assigned a made-up curvature.
    */
    const hs = new Array(n).fill(probe);
    const plus = new Array(n).fill(Infinity);
    const minus = new Array(n).fill(Infinity);
    const usable = [];
    for (let i = 0; i < n; i++) {
      if (shouldStop?.()) break;
      for (let tries = 0; tries < 5; tries++) {
        const a = u.slice();
        a[i] += hs[i];
        const b = u.slice();
        b[i] -= hs[i];
        const fa = at(a);
        const fb = at(b);
        if (Number.isFinite(fa) && Number.isFinite(fb)) {
          plus[i] = fa;
          minus[i] = fb;
          usable.push(i);
          break;
        }
        hs[i] /= 2;
      }
    }
    if (usable.length === 0) break;

    const m = usable.length;
    const g = usable.map((i) => (plus[i] - minus[i]) / (2 * hs[i]));
    const H = Array.from({ length: m }, () => new Array(m).fill(0));
    for (let k = 0; k < m; k++) {
      const i = usable[k];
      H[k][k] = (plus[i] - 2 * merit + minus[i]) / (hs[i] * hs[i]);
    }
    let broke = false;
    for (let k = 0; k < m && !broke; k++) {
      for (let l = k + 1; l < m; l++) {
        const i = usable[k];
        const j = usable[l];
        const mk = (si, sj) => {
          const a = u.slice();
          a[i] += si * hs[i];
          a[j] += sj * hs[j];
          return at(a);
        };
        const pp = mk(1, 1);
        const pm = mk(1, -1);
        const mp = mk(-1, 1);
        const mm = mk(-1, -1);
        if (![pp, pm, mp, mm].every(Number.isFinite)) {
          // A corner of the stencil is outside the smooth region. Treat the
          // pair as uncoupled rather than reading a cross-derivative off a
          // discontinuity.
          H[k][l] = H[l][k] = 0;
          continue;
        }
        H[k][l] = H[l][k] = (pp - pm - mp + mm) / (4 * hs[i] * hs[j]);
      }
    }
    if (broke) break;

    const { values, vectors } = symmetricEigen(H);
    // Measured against the STIFFEST curvature present, so "flat" is relative to
    // how sharply this merit surface curves at all rather than to an absolute
    // number in units nobody chose.
    const biggest = Math.max(...values.map(Math.abs), 1e-300);
    const cut = floor * biggest;

    /*
      The step, assembled direction by direction.

      Dividing by a near-zero eigenvalue is what makes an unguarded Newton step
      fly off along the direction that matters least, so those directions are
      left out - they are the null space, and there is nothing to gain along a
      flat one anyway.

      Directions of genuine NEGATIVE curvature are a different case and must
      not be treated the same way. A beamline that is not yet near its optimum
      is usually at a saddle - one combination of knobs improves the beam,
      another makes it worse - and skipping the downhill direction leaves the
      method stalled there. Measured before this: a lens-and-deflector column
      3 mm from its best took a 0.2 V step and stopped. Dividing by |lambda|
      keeps the step a descent direction while still scaling it by how sharply
      the merit is curving.
    */
    const step = new Array(n).fill(0);
    nullSpace = [];
    for (let k = 0; k < m; k++) {
      const lam = values[k];
      const vec = vectors[k];
      const gk = vec.reduce((s, c, i) => s + c * g[i], 0);
      if (Math.abs(lam) > cut) {
        for (let i = 0; i < m; i++) step[usable[i]] -= (gk / Math.abs(lam)) * vec[i];
      } else {
        nullSpace.push({
          curvature: lam,
          relative: lam / biggest,
          slope: gk,
          // Written back in volts per knob, so it reads as "this much of this
          // one against that much of that one".
          direction: usable.map((i, l) => ({
            label: knobs[i].label,
            weight: vec[l] * scales[i],
          })),
        });
      }
    }

    const size = Math.hypot(...step);
    if (!Number.isFinite(size) || size === 0) break;

    // Backtracking line search. A Newton step is only as good as the quadratic
    // it came from, and near an aperture edge that quadratic is a guess.
    let took = false;
    for (let shrink = 1; shrink >= 1 / 32; shrink /= 2) {
      if (shouldStop?.()) break;
      const trial = u.map((x, i) => x + shrink * step[i]);
      const m = at(trial);
      if (Number.isFinite(m) && m < merit * (1 - 1e-6)) {
        u = trial;
        merit = m;
        took = true;
        steps++;
        break;
      }
    }

    await onProgress?.({ iteration: iter + 1, merit, evaluations, nullSpace });
    if (!took) break;
  }

  put(u);
  const final = beamQuality(beamline, makeIons, flight);
  evaluations++;
  // Only keep the refinement if it really is better and still delivers.
  const improved = final.merit !== null && final.merit < meritAtStart * (1 - 1e-6);
  if (!improved) writeKnobs(beamline, knobs, start);

  return {
    values: readKnobs(beamline, knobs),
    start,
    improved,
    steps,
    evaluations,
    meritBefore: meritAtStart,
    meritAfter: improved ? final.merit : meritAtStart,
    rmsBefore: Math.sqrt(meritAtStart),
    rmsAfter: Math.sqrt(improved ? final.merit : meritAtStart),
    nullSpace,
  };
}
