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
import { createQuadrupole, QUADRUPOLE_DEFAULTS } from './quadrupole.js';
import { createBender, BENDER_DEFAULTS } from './bender.js';

/** @typedef {{key: string, label: string, unit?: string, min: number, max: number, step: number, rebuild: boolean, help?: string}} FieldSpec */

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
      { key: 'voltage', label: 'Centre potential', unit: 'V', min: -4000, max: 4000, step: 50, rebuild: false },
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
      // Step 1 V, not something rounder. The range has to reach tens of
      // kilovolts because real deflectors run there, but a 50 eV beam is
      // matched at a few tens of volts and its whole transmitting window is
      // about ten per cent wide — a 25 V step would step clean over it.
      { key: 'voltage', label: 'Electrode voltage', unit: 'V', min: -20000, max: 20000, step: 1, rebuild: false, help: 'Applied as +V and −V on opposite diagonals. The matched value for the current ion is shown below; Tune searches the solved field for the value that actually transmits best.' },
      { key: 'bendPlane', label: 'Bend plane', unit: '°', min: 0, max: 270, step: 90, rebuild: false, help: '0° turns the beam horizontally, 90° vertically. Same element, same solve — it is simply rolled about the beam.' },
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
      { key: 'rfAmplitude', label: 'RF amplitude', unit: 'V', min: 0, max: 2000, step: 10, rebuild: false, help: 'V, zero-to-peak. Set to 0 for a DC quadrupole.' },
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
