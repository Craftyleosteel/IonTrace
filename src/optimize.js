/**
 * Voltage tuning - search the solved fields for the settings that transmit
 * the most beam.
 *
 * Why this is cheap enough to be a button
 * ---------------------------------------
 * Every electrode voltage in this simulator is a scale factor on a stored unit
 * solution: the Laplace problem is solved once per electrode when the element
 * is built, and `setVoltage` only changes a multiplier (see docs/PHYSICS.md
 * section 3.3). So an optimiser may try a thousand voltage combinations
 * without ever touching the solver. The entire cost is flying the beam, and
 * that is the same cost as pressing Fly.
 *
 * That is not true of the geometry parameters. Changing an aperture radius or
 * an electrode thickness needs a fresh solve, so this optimiser does not touch
 * them - it tunes voltages only, which is also what an operator at a real
 * instrument can do.
 *
 * What is being maximised
 * -----------------------
 * Transmission is an integer count, so as an objective on its own it is flat
 * almost everywhere: most voltages transmit nothing, and the search has no
 * idea which way to move. The score below is that count plus two much smaller
 * terms that break the flatness:
 *
 *   1. partial credit for how far a lost ion got down the column, so a setting
 *      that destroys the beam in the last element beats one that destroys it in
 *      the first;
 *   2. a tie-break on the size of the surviving beam at the exit, so among
 *      settings that all transmit everything the tightest one wins.
 *
 * Both are scaled so they can never outweigh a single transmitted ion. The
 * reported figure of merit is always the honest count.
 *
 * The search
 * ----------
 * Coordinate descent with a coarse-to-fine bracketed scan on each knob. Not
 * gradient descent, and deliberately not: the objective is discontinuous
 * (an ion either clears an aperture or does not), so a derivative is
 * meaningless and a scan is what finds the operating window. Each knob gets a
 * coarse sweep of its whole range, then successively narrower sweeps around
 * the best point. Knobs are then revisited, because they interact - matching a
 * lens into a deflector changes which deflector voltage is best.
 */

import { flyBeam } from './integrator.js';
import { ELEMENT_TYPES } from './elements/index.js';
import { toLocal } from './frames.js';
import { refineNullSpace } from './reduced.js';

/**
 * Which parameters this optimiser is allowed to move.
 *
 * Keyed by element type, listing the parameter keys that scale a stored unit
 * solution rather than changing the geometry. Every one of these is declared
 * `rebuild: false` in the element registry, and `tunableKnobs` asserts that,
 * so the two cannot drift apart.
 */
export const TUNABLE = {
  aperture: ['voltage'],
  einzel: ['voltage'],
  bender: ['voltage'],
  // The RF quadrupole is deliberately NOT here, and this is the one entry
  // worth arguing about. Its voltages do change transmission, so they look
  // tunable - but a mass filter's job is selectivity, and the setting that
  // transmits the most beam is the one that filters nothing. Asked to
  // maximise transmission with the RF amplitude in its hands, the search
  // turns the RF down until the rods stop selecting: measured on a
  // lens-filter-deflector column, it dropped 250 V to 40 V and called it an
  // improvement. That is the optimiser working correctly on the wrong
  // objective. Tuning a filter means tuning (a, q) for a target mass, which
  // is a different search with a different goal, and it does not belong
  // behind a button labelled "best transmission".
};

/**
 * How far either side of an element's own suggested voltage to search.
 *
 * Symmetric about zero rather than about the suggestion, because the
 * suggestion is a magnitude: a negative ion wants the opposite polarity, and
 * the search should be able to find that without being told.
 */
const SEED_SPAN = 3;

/** The step `scoreBeamline` flies at when the caller names none. */
const DEFAULT_SCAN_CFL = 0.05;

/**
 * The knobs available on a beamline, with the range each may be searched over.
 *
 * Hard limits come from the element registry, so the optimiser can never
 * propose a voltage the user could not have typed in themselves. The range it
 * actually sweeps can be narrower, and for a deflector it should be: the
 * registry has to allow tens of kilovolts, because real deflectors run there,
 * while a 50 eV beam wants a few tens of volts. Sweeping the full range at any
 * affordable number of samples would step straight over the answer.
 *
 * So where an element can say what voltage it expects - the deflector can,
 * from a closed form - the sweep is centred on that instead. Give this
 * function the ion the column is being tuned for and it will ask.
 *
 * @param {object} beamline
 * @param {{energy?: number, charge?: number}} [ion] the ion, for elements that
 *        can suggest a voltage for it. Omit and every knob gets its full range.
 */
export function tunableKnobs(beamline, ion = null) {
  const knobs = [];
  beamline.elements.forEach((element, index) => {
    const type = element.typeKey;
    const spec = ELEMENT_TYPES[type];
    if (!spec) return;
    for (const key of TUNABLE[type] ?? []) {
      const field = spec.fields.find((f) => f.key === key);
      if (!field) continue;
      if (field.rebuild) {
        throw new Error(
          `${type}.${key} is listed as tunable but needs a re-solve; ` +
            'the optimiser only moves parameters that scale a stored solution'
        );
      }

      /*
        The element's own estimate of what this knob wants, from the same
        closed form the control's range and starting value come from. Taken
        from the registry's `scale` rather than from a method on the element,
        so every field that has one gets it - a lens wants -6 T/q and a filter
        a Mathieu q of 0.38, just as a deflector wants its matched voltage.

        Missing it is not cosmetic. Without a seed the knob's natural size
        falls back to a sixth of its full range, which for a lens allowed
        +/- 40 kV is 13 kV - so the sweep steps in kilovolts around an answer
        of a few hundred volts, and the finite-difference probe that the
        second stage uses is meaningless.
      */
      let seed = null;
      if (ion && typeof field.scale === 'function') {
        const v = field.scale(element.params, ion);
        if (Number.isFinite(v) && v !== 0) seed = Math.abs(v);
      }

      const span = seed === null ? null : seed * SEED_SPAN;
      knobs.push({
        index,
        type,
        key,
        label: `${element.label ?? spec.label} · ${field.label}`,
        unit: field.unit ?? '',
        min: field.min,
        max: field.max,
        // The range actually swept. Never wider than the registry allows.
        lo: span === null ? field.min : Math.max(field.min, -span),
        hi: span === null ? field.max : Math.min(field.max, span),
        seed,
        step: field.step,
        value: element.params[key] ?? 0,
      });
    }
  });
  return knobs;
}

/** Apply a knob's value to its element. */
export function applyKnob(beamline, knob, value) {
  const element = beamline.elements[knob.index];
  element.params[knob.key] = value;
  // Elements expose a setter for their principal voltage; the rest of the
  // tunable parameters are read straight out of `params` on every field
  // evaluation, so writing the parameter is enough.
  if (knob.key === 'voltage' && typeof element.setVoltage === 'function') {
    element.setVoltage(value);
  }
  // With fringe fields on, the field the ion actually flies through belongs to
  // a column solve spanning several elements, which keeps its own copy of the
  // voltages. Still fast adjust - the search does not re-solve anything.
  beamline.syncRuns?.();
}

/** Current value of every knob, for saving and restoring a trial. */
export function readKnobs(beamline, knobs) {
  return knobs.map((k) => beamline.elements[k.index].params[k.key] ?? 0);
}

/** Write a whole vector of knob values back. */
export function writeKnobs(beamline, knobs, values) {
  knobs.forEach((k, i) => applyKnob(beamline, k, values[i]));
}

/**
 * Fly a beam and score the result.
 *
 * `makeIons` must return a FRESH array each call: ions carry their own state
 * and are mutated by the integrator, so reusing them would score the second
 * trial against ions that have already flown.
 */
export function scoreBeamline(beamline, makeIons, opts = {}) {
  const ions = makeIons();
  const n = ions.length;
  if (n === 0) return { score: 0, transmitted: 0, count: 0, exitRadius: null };

  /*
    Transmitted WHERE.

    On a straight column that question does not arise. On a branching one it
    is the whole question: a deflector with hardware on both exits transmits
    everything at zero volts, straight out the back, and a search told merely
    to maximise transmission would discover that and switch the deflector off.
    That is not tuning a beamline, it is unplugging it.

    So an ion counts only if it leaves by the end being aimed at - the main
    line unless the caller names another.
  */
  const { target = beamline.mainEnd, ...flightOpts } = opts;
  const { tracks } = flyBeam(beamline, ions, {
    cfl: 0.05,
    maxSteps: 200000,
    ...flightOpts,
  });

  const last = beamline.elements[beamline.elements.length - 1];
  const path = last?.zEnd || 1;
  const exit = beamline.exitFrame;
  let transmitted = 0;
  let progress = 0;
  let sumR2 = 0;

  for (const t of tracks) {
    const end = t.points[t.points.length - 1];
    const arrived =
      t.stop === 'exited' &&
      (!target ||
        (() => {
          const e = beamline.endNearest(end.x, end.y ?? 0, end.z);
          return e && e.element === target.element && e.port === target.port;
        })());

    if (arrived) {
      transmitted++;
      // Transverse to the EXIT axis, not to the global z. After a bend those
      // are different directions, and the distance from the origin is mostly
      // the bend offset - a number that says nothing about the beam.
      const [ex, ey] = toLocal(exit, [end.x, end.y ?? 0, end.z]);
      sumR2 += ex * ex + ey * ey;
    } else if (t.stop === 'exited') {
      // Out of the instrument, but not where it was wanted. No credit, and no
      // partial credit either - it did not get "most of the way" anywhere.
    } else {
      // How far down the column it got, as a fraction of the path length.
      // `locate` returns null for a point already outside every element, which
      // for a lost ion means it left sideways - no credit for that.
      const hit = beamline.locate([end.x, end.y ?? 0, end.z]);
      progress += hit ? Math.min(1, hit.element.zEnd / path) : 0;
    }
  }

  const exitRadius = transmitted ? Math.sqrt(sumR2 / transmitted) : null;
  const scale = last?.bore || 0.005;

  // Partial credit and the tie-break are both capped below 1/n, so neither can
  // ever outweigh one more transmitted ion.
  const partial = (0.4 / n) * (progress / n);
  const tightness = exitRadius === null ? 0 : (0.4 / n) * Math.exp(-exitRadius / scale);

  return {
    score: transmitted / n + partial + tightness,
    transmitted,
    count: n,
    exitRadius,
  };
}

/** Round a value onto the knob's step grid and clamp it to the knob's range. */
function snap(knob, value) {
  const stepped = Math.round(value / knob.step) * knob.step;
  return Math.min(knob.max, Math.max(knob.min, stepped));
}

/**
 * Optimise a set of knobs for transmission.
 *
 * Runs as an async generator-driven loop so a browser can stay responsive and
 * the caller can cancel: `onProgress` is awaited, so yielding to the event
 * loop there is enough, and `shouldStop` is consulted between evaluations.
 *
 * @param {object} beamline
 * @param {() => object[]} makeIons   fresh ions for each trial
 * @param {object[]} knobs            from `tunableKnobs`
 * @param {object} [options]
 * @param {number} [options.passes=2]      sweeps over the whole knob list
 * @param {number} [options.coarse=13]     samples in the first sweep of a knob
 * @param {number} [options.refine=7]      samples in each narrowing sweep
 * @param {boolean} [options.polish=true]  run the reduced-Hessian second stage
 * @param {number} [options.levels=3]      narrowing sweeps after the coarse one
 * @param {number} [options.shrink=0.25]   bracket width kept at each level
 * @param {object} [options.flight]        options passed to the integrator
 * @param {(p: object) => any} [options.onProgress]
 * @param {() => boolean} [options.shouldStop]
 */
export async function optimizeVoltages(beamline, makeIons, knobs, options = {}) {
  const {
    passes = 2,
    coarse = 13,
    refine: refineSamples = 7,
    levels = 3,
    polish = true,
    scanSpeed = 4,
    flight = {},
    onProgress,
    shouldStop,
  } = options;

  /*
    The scan flies at a coarser time step than the answer is finally judged at.

    Its job is to RANK settings, not to measure one, and ranking survives a
    much coarser integration than measuring does. Measured on a lens-and-
    deflector column over a sweep of lens voltages: at two, four and eight
    times the default step the order of the settings is unchanged and the
    scores agree to four decimal places, while a flight drops from 70 ms to
    36, 17 and 8.

    The setting the scan chooses is then re-scored at full fidelity before
    anything is reported or refined, so the coarse step never reaches the
    answer - only the search path to it.

    Fewer IONS would also be cheaper and is not done: the same sweep reorders
    when the beam is thinned to five, because the tie-break is a mean over
    whichever ions are present and a different sample means a different mean.
  */
  const scanFlight = { ...flight, cfl: (flight.cfl ?? DEFAULT_SCAN_CFL) * scanSpeed };

  const start = readKnobs(beamline, knobs);
  let best = readKnobs(beamline, knobs);
  let bestResult = scoreBeamline(beamline, makeIons, scanFlight);
  let evaluations = 1;
  let cancelled = false;

  /*
    Try every knob at its own suggested value AT ONCE, before moving anything
    one at a time.

    Coordinate descent cannot find a setting that only works when two knobs are
    right together, because it never has both right at the same time: sweeping
    the lens while the deflector is off transmits nothing whatever the lens
    does, so the sweep learns nothing and moves on. Measured on a
    lens-and-deflector column starting from zero volts, the one-at-a-time
    search found nothing at all - and every element already knows roughly what
    it wants, so the combination of those guesses costs one flight to try.
  */
  if (knobs.some((k) => k.seed != null)) {
    const seeded = knobs.map((k, i) => (k.seed != null ? snap(k, k.seed) : start[i]));
    writeKnobs(beamline, knobs, seeded);
    const result = scoreBeamline(beamline, makeIons, scanFlight);
    evaluations++;
    if (result.score > bestResult.score) {
      bestResult = result;
      best = seeded;
    }
    writeKnobs(beamline, knobs, best);
  }

  // An upper bound, not an estimate: the seeds add at most three samples to
  // each knob's first level, and repeats of a value already tried are skipped
  // rather than re-flown. So the count can finish early but never overrun,
  // which is the direction a progress bar should err in.
  const perKnob = coarse + 3 + refineSamples * levels;
  const total = 1 + passes * knobs.length * perKnob;

  const report = async (knob) => {
    if (!onProgress) return;
    await onProgress({
      evaluations,
      total,
      fraction: Math.min(1, evaluations / total),
      transmitted: bestResult.transmitted,
      count: bestResult.count,
      knob: knob?.label ?? null,
      values: best.slice(),
    });
  };

  /*
    Settings already flown.

    Coordinate descent revisits the same point often: every narrowing level
    re-samples its own bracket centre, and a second pass re-walks ranges the
    first pass already covered. The flight is deterministic - the same ions
    through the same fields - so the answer is too, and re-flying it is pure
    waste. Keyed on the whole knob vector, because that is what the score
    depends on.
  */
  const flown = new Map();

  /** Try one value on one knob, keeping it only if it scores better. */
  const trial = (knob, ki, value) => {
    applyKnob(beamline, knob, value);
    const key = readKnobs(beamline, knobs).join(',');
    let result = flown.get(key);
    if (!result) {
      result = scoreBeamline(beamline, makeIons, scanFlight);
      evaluations++;
      flown.set(key, result);
    }
    if (result.score > bestResult.score) {
      bestResult = result;
      best = readKnobs(beamline, knobs);
      best[ki] = value;
      return true;
    }
    return false;
  };

  outer: for (let pass = 0; pass < passes && !cancelled; pass++) {
    for (let ki = 0; ki < knobs.length; ki++) {
      const knob = knobs[ki];

      // Every trial of this knob starts from the best setting of the others,
      // which is what makes this coordinate descent rather than independent
      // one-at-a-time scans.
      writeKnobs(beamline, knobs, best);

      let lo = knob.lo ?? knob.min;
      let hi = knob.hi ?? knob.max;
      let centre = best[ki];

      for (let level = 0; level <= levels; level++) {
        const samples = level === 0 ? coarse : refineSamples;
        const spacing = (hi - lo) / (samples - 1 || 1);

        const values = [];
        for (let s = 0; s < samples; s++) values.push(lo + spacing * s);
        if (level === 0) {
          // Always try the element's own estimate, both polarities, and the
          // setting already in place - explicitly, rather than trusting the
          // sweep to land on them. Whether an evenly spaced scan happens to
          // include the seed depends on the sample count, which is not
          // something the physics should hang on.
          if (knob.seed != null) values.push(knob.seed, -knob.seed);
          values.push(knob.value);
        }

        const tried = new Set();
        for (const raw of values) {
          if (shouldStop?.()) {
            cancelled = true;
            break outer;
          }
          const value = snap(knob, raw);
          if (tried.has(value)) continue;
          tried.add(value);
          if (trial(knob, ki, value)) centre = best[ki];
          writeKnobs(beamline, knobs, best);
        }
        await report(knob);

        // Refine within one sample spacing of the best point. That is the
        // right width by construction: a scan at spacing d can only have
        // missed the optimum by less than d, so bracketing +/- d is what makes
        // the next level able to find what this one stepped over. A fixed
        // fraction of the range would not - on a knob whose range is three
        // orders of magnitude wider than its answer, it never closes in.
        lo = Math.max(knob.min, centre - spacing);
        hi = Math.min(knob.max, centre + spacing);
        if (hi - lo < knob.step) break;
      }
    }
  }

  writeKnobs(beamline, knobs, best);

  // Back to full fidelity. Everything above ranked settings against each other
  // at a coarse step; what gets reported, and what the refinement works from,
  // is measured at the step the caller asked for.
  if (scanSpeed !== 1) {
    bestResult = scoreBeamline(beamline, makeIons, flight);
    evaluations++;
  }

  /*
    Second stage: once the beam is through, tighten it.

    The scan answers "does it get there", which is a discontinuous question and
    needs a scan. It cannot answer "of the settings that work, which is best",
    because the difference between them is a smooth quantity and a scan on a
    grid of one step per knob steps straight over it. That is a Newton problem,
    and it is where second-order information earns its keep - a lens and the
    deflector behind it trade off, so the optimum lies in a valley at an angle
    to both axes, which is the case coordinate descent is worst at.

    Skipped unless the beam is fully delivered: with ions still being lost the
    merit is not smooth, and a Hessian of it would be fiction.
  */
  let refinement = null;
  if (
    polish &&
    !cancelled &&
    knobs.length >= 1 &&
    bestResult.count > 0 &&
    bestResult.transmitted === bestResult.count
  ) {
    refinement = await refineNullSpace(beamline, makeIons, knobs, {
      flight,
      shouldStop,
      onProgress: async (p) => {
        evaluations += 0;
        await onProgress?.({
          evaluations,
          total,
          fraction: 1,
          transmitted: bestResult.transmitted,
          count: bestResult.count,
          knob: 'tightening the beam',
          values: readKnobs(beamline, knobs),
          refining: p,
        });
      },
    });
    evaluations += refinement.evaluations;
    if (refinement.improved) {
      best = readKnobs(beamline, knobs);
      bestResult = scoreBeamline(beamline, makeIons, flight);
      evaluations++;
    }
  }

  writeKnobs(beamline, knobs, best);
  const improved = bestResult.score > 0 && best.some((v, i) => v !== start[i]);

  return {
    values: best,
    start,
    improved,
    cancelled,
    evaluations,
    transmitted: bestResult.transmitted,
    count: bestResult.count,
    exitRadius: bestResult.exitRadius,
    refinement,
  };
}
