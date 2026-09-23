/**
 * Element registry.
 *
 * One entry per kind of optic, each describing how to build it and which of
 * its parameters are editable. The UI is generated from this rather than
 * hard-coded, so adding an element means adding it here and nowhere else -
 * which is the point of the whole architecture.
 *
 * `rebuild` says whether changing a parameter needs a fresh Laplace solve.
 * Voltages never do: they are a weighted sum of solutions already computed,
 * which is why they are instant and geometry is not.
 */

import { createDrift, DRIFT_DEFAULTS } from './drift.js';
import { createAperture, APERTURE_DEFAULTS } from './aperture.js';
import { createEinzel, EINZEL_ELEMENT_DEFAULTS } from './einzel.js';
import { createQuadrupole, QUADRUPOLE_DEFAULTS, amplitudeForQ } from './quadrupole.js';
import { createBender, BENDER_DEFAULTS, matchedVoltage } from './bender.js';

/** @typedef {{key: string, label: string, unit?: string, min: number, max: number, step: number, rebuild: boolean, help?: string, scale?: Function}} FieldSpec */

/**
 * Toolbar glyphs, drawn on a 32x18 canvas with the beam running left to right.
 *
 * Each one shows the element's actual cross-section rather than an abstract
 * symbol, so the toolbar reads as a picture of what will be placed: two rails
 * for a drift, a plate with a hole, three cylinders, four rods, a curve.
 */
const ICONS = {
  drift: '<path d="M2 4h28M2 14h28"/>',
  aperture: '<path d="M16 1v6M16 11v6"/><path d="M2 9h28" opacity=".35"/>',
  einzel: '<path d="M3 4h7M3 14h7M13 4h6M13 14h6M22 4h7M22 14h7"/>',
  quadrupole:
    '<circle cx="9" cy="5" r="2.6"/><circle cx="23" cy="5" r="2.6"/>' +
    '<circle cx="9" cy="13" r="2.6"/><circle cx="23" cy="13" r="2.6"/>',
  // Four curved electrodes in a grounded box, with the beam entering on one
  // axis and leaving on the perpendicular one.
  bender:
    '<rect x="9.5" y="1.5" width="15" height="15" rx="1" opacity=".3"/>' +
    '<path d="M21.7 7.29A5 5 0 0 0 18.71 4.3"/>' +
    '<path d="M15.29 4.3A5 5 0 0 0 12.3 7.29"/>' +
    '<path d="M12.3 10.71A5 5 0 0 0 15.29 13.7"/>' +
    '<path d="M18.71 13.7A5 5 0 0 0 21.7 10.71"/>' +
    '<path d="M2 9h15V1" opacity=".55"/>',
};

export const ELEMENT_TYPES = {
  drift: {
    icon: ICONS.drift,
    label: 'Drift',
    blurb: 'Field-free tube. Separates active elements so their fringe fields do not overlap.',
    create: createDrift,
    defaults: DRIFT_DEFAULTS,
    fields: [
      { key: 'length', label: 'Length', unit: 'mm', min: 2, max: 120, step: 1, rebuild: true },
      { key: 'bore', label: 'Bore', unit: 'mm', min: 2, max: 20, step: 0.5, rebuild: true },
    ],
  },

  aperture: {
    icon: ICONS.aperture,
    label: 'Aperture plate',
    blurb:
      'A charged plate with a hole. Not a passive opening — the equipotentials bulge through it, so it acts as a lens.',
    create: createAperture,
    defaults: APERTURE_DEFAULTS,
    fields: [
      { key: 'voltage', label: 'Potential', unit: 'V', min: -3000, max: 3000, step: 25, rebuild: false },
      { key: 'bore', label: 'Hole radius', unit: 'mm', min: 1, max: 12, step: 0.5, rebuild: true },
      { key: 'thickness', label: 'Thickness', unit: 'mm', min: 0.5, max: 10, step: 0.5, rebuild: true },
      { key: 'margin', label: 'Drift margin', unit: 'mm', min: 4, max: 40, step: 1, rebuild: true },
    ],
  },

  einzel: {
    icon: ICONS.einzel,
    label: 'Einzel lens',
    blurb:
      'Three coaxial cylinders, outer two grounded. Does no net work on a transmitted ion, and focuses for either polarity.',
    create: createEinzel,
    defaults: EINZEL_ELEMENT_DEFAULTS,
    fields: [
      // Electrostatic optics depends only on E/q, so what a lens does is set
      // by the centre potential as a multiple of T/q and by nothing else —
      // measured identical at 50 eV and 500 eV. Six times T/q is comfortably
      // inside the range that transmits everything (about 1 to 10 times) and
      // well clear of the over-focusing that starts near 20.
      //
      // Divided by the SIGNED charge, so the polarity follows the ion. A
      // negative ion meeting a negative centre electrode is decelerated, not
      // accelerated, and at six times its own energy it cannot climb the
      // barrier at all — it is reflected, and the lens transmits nothing.
      { key: 'voltage', label: 'Centre potential', unit: 'V', min: -40000, max: 40000, step: 50, rebuild: false, scale: (params, ion) => (-6 * ion.energy) / (ion.charge || 1) },
      { key: 'boreRadius', label: 'Bore', unit: 'mm', min: 3, max: 10, step: 0.5, rebuild: true },
      { key: 'centreLength', label: 'Centre length', unit: 'mm', min: 5, max: 40, step: 1, rebuild: true },
      { key: 'gap', label: 'Gap', unit: 'mm', min: 1, max: 12, step: 0.5, rebuild: true },
    ],
  },

  bender: {
    icon: ICONS.bender,
    label: 'Quadrupole deflector',
    blurb:
      'Four curved electrodes in a grounded box, at +V and −V on the diagonals, turning the beam ninety degrees. Its exit faces a different way from its entrance, so everything after it turns too — this is what makes the column a path rather than a line.',
    create: createBender,
    defaults: BENDER_DEFAULTS,
    fields: [
      // The hard limits have to reach tens of kilovolts, because real
      // deflectors run there — a 27 keV device with thin electrodes wants
      // 45 kV. But a 50 eV beam is matched at forty volts, and no single
      // fixed slider can serve both: across ±60 kV one pixel is over 500 V,
      // so the answer for a low-energy beam is unreachable by dragging. Hence
      // `scale`: the slider spans a few times the matched voltage for the ion
      // actually in the source, and follows it when the beam changes.
      // Signed by the charge: the matched voltage is a magnitude, but which
      // diagonal carries it decides which way the beam turns, and a negative
      // ion at a positive ion's polarity is steered into the wall.
      { key: 'voltage', label: 'Electrode voltage', unit: 'V', min: -60000, max: 60000, step: 1, rebuild: false, scale: (params, ion) => matchedVoltage(params, ion.energy, Math.abs(ion.charge) || 1) * Math.sign(ion.charge || 1), help: 'Applied as +V and −V on opposite diagonals. A new deflector arrives already matched to the beam.' },
      // Any angle, not just the four right angles it used to offer. The roll
      // is a rotation about the beam applied when the field is evaluated, so
      // an arbitrary angle costs nothing and needs no re-solve — there was
      // never a reason for the quarter-turn steps beyond the slider that used
      // to set it.
      { key: 'bendPlane', label: 'Rotation about the beam', unit: '°', min: -360, max: 360, step: 5, rebuild: false, help: '0° turns the beam left, 90° turns it down, 180° right, 270° up — and anything in between bends into a plane at that angle. The element is unchanged; it is simply rolled.' },
      { key: 'apertureRadius', label: 'Aperture radius r₀', unit: 'mm', min: 4, max: 40, step: 0.5, rebuild: true, help: 'Centre to the concave electrode faces. The matched voltage goes as (r₀/a)², so this and the two below set the operating voltage between them.' },
      { key: 'electrodeThickness', label: 'Electrode thickness', unit: 'mm', min: 0.5, max: 20, step: 0.5, rebuild: true },
      { key: 'boxClearance', label: 'Box clearance', unit: 'mm', min: 0.5, max: 20, step: 0.5, rebuild: true, help: 'Electrode backs to the grounded box. r₀ plus these two is the half-width a.' },
      { key: 'gapAngle', label: 'Aperture gap', unit: '°', min: 5, max: 40, step: 1, rebuild: true, help: 'How far short of each axis the electrodes stop. This is the beam’s way in and out.' },
      { key: 'height', label: 'Vertical aperture', unit: 'mm', min: 4, max: 60, step: 1, rebuild: true },
    ],
  },

  quadrupole: {
    icon: ICONS.quadrupole,
    label: 'Quadrupole',
    blurb:
      'Four rods, RF and DC. Converges in one transverse plane and diverges in the other at every instant; the RF is what makes it net-focusing in both.',
    create: createQuadrupole,
    defaults: QUADRUPOLE_DEFAULTS,
    fields: [
      // Stability depends on mass, so a fixed amplitude is stable only for
      // the ion it was written for. This one is set from the Mathieu q the
      // current ion would sit at — see amplitudeForQ.
      { key: 'rfAmplitude', label: 'RF amplitude', unit: 'V', min: 0, max: 20000, step: 10, rebuild: false, scale: (params, ion) => amplitudeForQ(params, ion.mass, ion.charge), help: 'V, zero-to-peak. Set to 0 for a DC quadrupole.' },
      { key: 'dcVoltage', label: 'DC offset', unit: 'V', min: -500, max: 500, step: 5, rebuild: false, help: 'U. The ratio U/V is what makes a filter selective.' },
      { key: 'frequency', label: 'Frequency', unit: 'MHz', min: 0.1, max: 5, step: 0.05, rebuild: false },
      { key: 'phase', label: 'Entry phase', unit: '°', min: 0, max: 360, step: 5, rebuild: false, help: 'RF phase when the simulation starts. Real transmission depends on it.' },
      { key: 'length', label: 'Rod length', unit: 'mm', min: 10, max: 200, step: 5, rebuild: false },
      { key: 'fieldRadius', label: 'Field radius r₀', unit: 'mm', min: 1.5, max: 10, step: 0.25, rebuild: true },
      { key: 'rodRadius', label: 'Rod radius', unit: 'mm', min: 1.5, max: 12, step: 0.1, rebuild: true, help: 'The ratio 1.1487 × r₀ cancels the 12-pole for round rods.' },
    ],
  },
};

/** Build an element of `type` with the given parameters. */
export function createElement(type, params = {}, solverOpts = {}) {
  const spec = ELEMENT_TYPES[type];
  if (!spec) throw new Error(`Unknown element type "${type}"`);
  const element = spec.create(params, solverOpts);
  element.typeKey = type;
  return element;
}

/** True if changing `key` on `type` requires a fresh Laplace solve. */
export function needsRebuild(type, key) {
  const spec = ELEMENT_TYPES[type];
  const f = spec?.fields.find((x) => x.key === key);
  // Unknown keys are treated as structural, which is the safe direction to
  // err in: a needless re-solve is slow, a skipped one is wrong.
  return f ? f.rebuild : true;
}

/** How far either side of its natural value a scaled slider should reach. */
const SLIDER_SPAN = 3;

/** Roughly this many positions on a slider, before rounding to a round step. */
const SLIDER_STEPS = 400;

/** The largest round number (1, 2, 2.5 or 5 x a power of ten) at most `x`. */
function niceStep(x) {
  if (!(x > 0)) return 1;
  const decade = 10 ** Math.floor(Math.log10(x));
  const f = x / decade;
  return (f >= 5 ? 5 : f >= 2.5 ? 2.5 : f >= 2 ? 2 : 1) * decade;
}

/**
 * The range a slider should offer for a field, which is not always the range
 * the parameter is allowed to take.
 *
 * Most fields are the same either way. A deflector's voltage is not: its hard
 * limits span ±60 kV because real deflectors run there, while the beam in
 * front of it may want forty volts. Offering the full span would put the
 * answer less than a pixel wide and make the control useless in exactly the
 * case a beginner meets first.
 *
 * So a field may declare `scale(params, ion)` - the value it naturally sits
 * near - and gets a slider a few times that instead. Two rules keep it safe:
 * the range is always clipped to the hard limits, so no slider can propose a
 * setting the parameter does not allow; and it is always widened to contain
 * the value already set, so re-rendering after a beam change can never
 * silently clamp a voltage the user chose deliberately.
 *
 * @param {FieldSpec} field
 * @param {object} params  the element's current parameters
 * @param {{energy: number, charge: number}} [ion]
 */
export function fieldRange(field, params, ion = null) {
  const base = { min: field.min, max: field.max, step: field.step };
  if (!field.scale || !ion) return base;

  const scale = Math.abs(field.scale(params, ion));
  if (!Number.isFinite(scale) || scale === 0) return base;
  return scaledSlider(field, scale, params[field.key]);
}

/**
 * A slider a few times as wide as `scale`, clipped to the field's hard limits
 * and widened to contain `set` if there is one.
 *
 * Shared with `startingParams` so the two agree by construction: a starting
 * value that did not land on its own slider's step would be rounded by the
 * browser the moment the panel drew it, and the element would quietly hold a
 * different voltage from the one on screen.
 */
function scaledSlider(field, scale, set) {
  // Symmetric about zero: the natural value is a magnitude, and a negative
  // ion wants the opposite polarity.
  let hi = Math.min(field.max, scale * SLIDER_SPAN);
  let lo = Math.max(field.min, -hi);

  if (Number.isFinite(set)) {
    hi = Math.min(field.max, Math.max(hi, set));
    lo = Math.max(field.min, Math.min(lo, set));
  }

  return { min: lo, max: hi, step: niceStep((hi - lo) / SLIDER_STEPS) };
}

/**
 * Parameters a newly placed element should start with, given the beam.
 *
 * A deflector at zero volts is not a deflector - it is a box the beam flies
 * into. Dropping one in from the toolbar and watching nothing happen is a
 * worse first impression than the element deserves, and the matched voltage
 * is known in closed form, so a new one arrives already set for the ion in
 * the source. Everything else keeps its own defaults.
 *
 * Only the initial value. Nothing here ever moves a voltage the user has set:
 * change the beam energy afterwards and the deflector stays where it is, and
 * the readout says how far off it now is - which is the energy selectivity of
 * the device, and worth seeing rather than hiding.
 */
export function startingParams(type, ion = null) {
  const spec = ELEMENT_TYPES[type];
  if (!spec || !ion) return {};
  const out = {};
  for (const field of spec.fields) {
    if (!field.scale) continue;
    const v = field.scale(spec.defaults, ion);
    if (!Number.isFinite(v) || v === 0) continue;
    // Snapped to the step the field will use for this value, so the box shows
    // exactly what the element holds. `toPrecision` clears the float dust that
    // multiplying a step back out leaves behind - without it a 4 u ion gets a
    // quadrupole reading 9.950000000000001 V.
    const { step } = scaledSlider(field, Math.abs(v), null);
    out[field.key] = Number((Math.round(v / step) * step).toPrecision(12));
  }
  return out;
}
