/**
 * IonTrace front end.
 *
 * Wiring only: every physical quantity on screen is computed by the modules in
 * this directory, and nothing here adjusts a number to make a picture look
 * better.
 *
 * Two cost paths are kept visible to the user because they are the point of
 * the architecture:
 *
 *   changing a voltage   -> fast adjust, a weighted sum of stored solutions
 *   changing geometry    -> rebuild that element and re-solve Laplace
 *
 * and only the element that changed is ever re-solved, never the whole line.
 */

import { Beamline } from './beamline.js';
import {
  ELEMENT_TYPES,
  createElement,
  needsRebuild,
  fieldRange,
  startingParams,
} from './elements/index.js';
import { MATHIEU_Q_LIMIT } from './elements/quadrupole.js';
import { discBeam, focalCrossing } from './ion.js';
import { createFlight, kineticEnergy } from './integrator.js';
import { tunableKnobs, optimizeVoltages, TUNABLE } from './optimize.js';
import {
  joulesToEV,
  mToMm,
  mmToM,
  ELEMENTARY_CHARGE,
  ATOMIC_MASS_UNIT,
} from './constants.js';
import { NO_ELECTRODE } from './grid.js';
import { toGlobal, toLocal, forwardOf } from './frames.js';

const el = (id) => document.getElementById(id);

const canvas = el('scene');
const ctx = canvas.getContext('2d');
const cross = el('cross');
const crossCtx = cross.getContext('2d');
const statusEl = el('status');
const readoutEl = el('readout');
const trackEl = el('track');
const inspectorEl = el('inspector');
const addPanel = el('addPanel');
const toolsEl = el('tools');
const beamPanel = el('beamPanel');
const autoAlignBtn = el('autoAlign');
const optimizeBtn = el('optimize');
const tuneNote = el('tuneNote');
const scaleNote = el('scaleNote');
const flyButton = el('fly');

const inputs = {
  mass: el('mass'),
  charge: el('charge'),
  energy: el('energy'),
  rays: el('rays'),
  beamRadius: el('beamRadius'),
  divergence: el('divergence'),
  repulsion: el('repulsion'),
  beamCurrent: el('beamCurrent'),
  ionsPerParticle: el('ionsPerParticle'),
  method: el('method'),
  cfl: el('cfl'),
  showField: el('showField'),
  showContours: el('showContours'),
};

/** Controls that change only how the scene is drawn, never the physics. */
const DISPLAY_INPUTS = ['showField', 'showContours'];

/* ------------------------------------------------------------------ */
/* colour                                                              */
/* ------------------------------------------------------------------ */

function cssVar(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

function cssRGB(name) {
  const hex = cssVar(name).replace('#', '');
  const full = hex.length === 3 ? hex.split('').map((c) => c + c).join('') : hex;
  return [
    parseInt(full.slice(0, 2), 16),
    parseInt(full.slice(2, 4), 16),
    parseInt(full.slice(4, 6), 16),
  ];
}

const mix = (a, b, t) => [
  Math.round(a[0] + (b[0] - a[0]) * t),
  Math.round(a[1] + (b[1] - a[1]) * t),
  Math.round(a[2] + (b[2] - a[2]) * t),
];

/**
 * Diverging colour for a normalised potential in [-1, 1].
 *
 * Two hues with a neutral midpoint, so sign is carried by hue and magnitude by
 * saturation, and zero potential reads as "nothing" rather than as a third
 * category. The square root stretches the low-magnitude end, where the fringe
 * fields that actually do the focusing live.
 */
function divergingColour(t, neg, zero, pos) {
  const m = Math.min(1, Math.sqrt(Math.abs(t)));
  return t < 0 ? mix(zero, neg, m) : mix(zero, pos, m);
}

/* ------------------------------------------------------------------ */
/* state                                                               */
/* ------------------------------------------------------------------ */

let beamline = null;

/**
 * What the inspector is showing.
 *
 * `null` means nothing is selected, and the panel collapses to the add menu
 * and the Fly button. The ion SOURCE is a selectable object in its own right,
 * so the beam's settings live where every other object's do rather than in a
 * panel that is always on screen.
 */
let selection = null; // null | {kind:'source'} | {kind:'element', index}

const selectedIndex = () => (selection?.kind === 'element' ? selection.index : -1);

let view = null; // the top pane's transform, for hit-testing
let views = []; // every pane currently drawn
let trajectories = [];
let stats = {};
let flight = null;
let flightSpec = null;
let flightOpts = null;
let flightStarted = 0;
let animation = { running: false, frame: 0, perFrame: 40 };
let stale = false;

const readNumber = (input, fallback) => {
  const v = parseFloat(input.value);
  return Number.isFinite(v) ? v : fallback;
};

/**
 * A starting column that shows what the thing does without any setup.
 *
 * A lens and a corner. The corner is the point: a column that turns is what
 * distinguishes this from a ray diagram, and leaving it out of the first
 * thing anyone sees hid the most interesting element behind a toolbar button.
 *
 * Chosen for how forgiving it is rather than how impressive. Measured with
 * the shipped disc beam, it transmits all nine ions with the lens anywhere
 * from -600 V to 0 and the deflector anywhere from 32 to 42 V, so the first
 * thing a newcomer changes does not collapse it. A filter here instead would
 * look better and behave far worse: with a quadrupole in the line the same
 * beam swings between 0 and 6 of 9 over 50 V of lens, because the filter's
 * acceptance is sharp and not monotonic. That is real, and worth meeting
 * deliberately rather than on load.
 */
function defaultBeamline() {
  return new Beamline([
    createElement('drift', { length: 12, bore: 5 }),
    createElement('einzel', {
      gridStep: 0.5,
      voltage: -300, // six times T/q for the 50 eV beam in the source
      boreRadius: 5,
      housingRadius: 14,
      entryDrift: 15,
      exitDrift: 15,
    }),
    createElement('drift', { length: 14, bore: 5 }),
    // Matched for a 100 u, 50 eV singly-charged ion: V0 = 1.8556 (T/q)(r0/a)^2
    // is 39.8 V, and the solved field turns that through 89 degrees.
    createElement('bender', { voltage: 40 }),
    createElement('drift', { length: 30, bore: 5 }),
  ]);
}

/* ------------------------------------------------------------------ */
/* element editing                                                     */
/* ------------------------------------------------------------------ */

/** Rebuild one element from its current parameters, re-solving its field. */
function rebuildElement(index) {
  const old = beamline.elements[index];
  const rebuilt = createElement(old.typeKey, old.params);
  beamline.replace(index, rebuilt);
  return rebuilt;
}

/**
 * Apply a parameter change to an element.
 *
 * Whether this costs a Laplace solve is decided by the element registry, not
 * guessed here. Voltages and RF settings only rescale solutions that already
 * exist; anything that moves metal changes the boundary and must be re-solved.
 */
function setParam(index, key, value) {
  const element = beamline.elements[index];
  element.params[key] = value;

  if (needsRebuild(element.typeKey, key)) {
    statusEl.classList.add('busy');
    rebuildElement(index);
    statusEl.classList.remove('busy');
  } else if (key === 'voltage' && element.setVoltage) {
    element.setVoltage(value);
  } else {
    // Quadrupole drive parameters are read straight off `params` at
    // evaluation time, so there is nothing to recompute at all.
    beamline.layout();
  }

  markStale();
  renderTrack();
  renderInspector();
  render();
  drawReadout();
}

function addElement(type) {
  const after = selectedIndex();
  const index = after >= 0 ? after + 1 : beamline.elements.length;
  statusEl.classList.add('busy');
  // Matched to the beam that is in the source right now, so a deflector
  // dropped in from the toolbar actually deflects.
  beamline.add(createElement(type, startingParams(type, beamSpec())), index);
  statusEl.classList.remove('busy');
  selection = { kind: 'element', index };
  afterStructureChange();
}

function removeElement(index) {
  if (beamline.elements.length <= 1) return;
  beamline.remove(index);
  selection = null;
  afterStructureChange();
}

function moveElement(index, delta) {
  if (!beamline.move(index, delta)) return;
  selection = { kind: 'element', index: index + delta };
  afterStructureChange();
}

function afterStructureChange() {
  markStale();
  renderTrack();
  renderInspector();
  render();
  drawReadout();
}

/* ------------------------------------------------------------------ */
/* direct manipulation on the canvas                                   */
/* ------------------------------------------------------------------ */

/**
 * Pointer interaction with the beamline drawing.
 *
 * Elements are dragged along the axis to reorder them, and the beam source is
 * dragged to set its radius. The column is a contiguous sequence rather than
 * a free layout - elements sit end to end with no gaps - so dragging along z
 * means "put this one somewhere else in the order", and the drop indicator
 * shows where it will land.
 */
let drag = null;

/** Canvas pixel position of a pointer event, in CSS pixels. */
function pointerPos(e) {
  const rect = canvas.getBoundingClientRect();
  return { px: e.clientX - rect.left, py: e.clientY - rect.top };
}

/**
 * The world point under the cursor, on the y = 0 plane.
 *
 * The view is a top view of the global x-z plane, so a screen position maps
 * back to a point in that plane exactly. Hit testing then happens in world
 * coordinates, which is what lets it keep working after the column bends.
 */
function paneAt(py) {
  return views.find((v) => py >= v.top && py < v.top + v.height) ?? views[0] ?? view;
}

function worldAt(px, py) {
  const T = paneAt(py);
  if (!T) return [0, 0, 0];
  return T.unproject(px, py);
}

/** The element under a world point, or -1. */
function elementIndexAtWorld(g) {
  for (let i = 0; i < beamline.elements.length; i++) {
    const e = beamline.elements[i];
    const l = toLocal(e.frame, g);
    if (!e.contains(l[0], l[1], l[2])) continue;
    if (Math.hypot(l[0], l[1]) <= e.outerRadius * 1.25) return i;
  }
  return -1;
}

/** Whether the cursor is on one of the source handles. */
function onSourceHandle(px, py) {
  const first = beamline.elements[0];
  const T = paneAt(py);
  if (!first || !T) return false;
  const radius = mmToM(readNumber(inputs.beamRadius, 1));
  for (const sign of [-1, 1]) {
    const [hx, hy] = T.project(toGlobal(first.frame, across(T, sign * radius)));
    if (Math.hypot(px - hx, py - hy) < 10) return true;
  }
  return false;
}

/**
 * Distance along the reference path nearest a world point.
 *
 * Used to decide where a dragged element would land. Once the column can bend
 * there is no coordinate to compare against, so the drop position comes from
 * the closest point on the orbit itself.
 */
function nearestPathDistance(g) {
  let best = Infinity;
  let bestS = 0;
  for (const e of beamline.elements) {
    const n = e.curved ? 16 : 2;
    for (let k = 0; k <= n; k++) {
      const f = k / n;
      const p = toGlobal(
        e.frame,
        e.pathPoint ? e.pathPoint(f) : [0, 0, f * e.length]
      );
      const d = Math.hypot(p[0] - g[0], p[2] - g[2]);
      if (d < best) {
        best = d;
        bestS = e.zStart + f * e.length;
      }
    }
  }
  return bestS;
}

/** Where an element dropped at this path distance would be inserted. */
function dropIndexAtS(s) {
  for (let i = 0; i < beamline.elements.length; i++) {
    const e = beamline.elements[i];
    if (s < e.zStart + e.length / 2) return i;
  }
  return beamline.elements.length;
}

function select(next) {
  selection = next;
  renderTrack();
  renderInspector();
  render();
}

canvas.addEventListener('pointerdown', (e) => {
  if (!beamline) return;
  const { px, py } = pointerPos(e);

  if (onSourceHandle(px, py)) {
    drag = { kind: 'beam' };
    canvas.setPointerCapture(e.pointerId);
    select({ kind: 'source' });
    return;
  }

  const index = elementIndexAtWorld(worldAt(px, py));
  if (index < 0) {
    // Clicking empty space deselects, which is what collapses the panel back
    // to the add menu and the Fly button.
    select(null);
    return;
  }

  select({ kind: 'element', index });
  drag = { kind: 'element', from: index, startPx: px, startPy: py, moved: false };
  canvas.setPointerCapture(e.pointerId);
  canvas.style.cursor = 'grabbing';
});

canvas.addEventListener('pointermove', (e) => {
  if (!beamline || !view) return;
  const { px, py } = pointerPos(e);

  if (!drag) {
    canvas.style.cursor = onSourceHandle(px, py)
      ? 'grab'
      : elementIndexAtWorld(worldAt(px, py)) >= 0
        ? 'grab'
        : 'default';
    return;
  }

  if (drag.kind === 'beam') {
    // Distance from the axis, measured in the source's own frame so the
    // handles keep working if the first element is nudged, and along the
    // transverse axis of whichever pane is being dragged in.
    const first = beamline.elements[0];
    const T = paneAt(py);
    const l = toLocal(first.frame, worldAt(px, py));
    const mm = Math.max(0.1, Math.min(6, Math.abs(mToMm(l[T.axis]))));
    inputs.beamRadius.value = mm.toFixed(1);
    syncOutputs();
    renderTrack();
    markStale();
    render();
    return;
  }

  if (Math.hypot(px - drag.startPx, py - drag.startPy) > 5) drag.moved = true;
  drag.dropIndex = dropIndexAtS(nearestPathDistance(worldAt(px, py)));
  render();
});

function endDrag(e) {
  if (!drag) return;
  const finished = drag;
  drag = null;
  canvas.style.cursor = 'default';
  try {
    canvas.releasePointerCapture(e.pointerId);
  } catch {
    /* already released */
  }

  if (finished.kind === 'element' && finished.moved) {
    let target = finished.dropIndex ?? finished.from;
    // Removing the element first shifts everything after it down one.
    if (target > finished.from) target -= 1;
    if (target !== finished.from) {
      const [moved] = beamline.elements.splice(finished.from, 1);
      beamline.elements.splice(target, 0, moved);
      beamline.layout();
      selection = { kind: 'element', index: target };
      afterStructureChange();
      return;
    }
  }
  render();
  drawReadout();
}

canvas.addEventListener('pointerup', endDrag);
canvas.addEventListener('pointercancel', endDrag);

/* ------------------------------------------------------------------ */
/* track and inspector                                                 */
/* ------------------------------------------------------------------ */

const escapeHtml = (s) =>
  String(s).replace(/[&<>"]/g, (c) => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;',
  })[c]);

/** One-line summary of what an element is currently set to. */
function summarise(e) {
  const p = e.params;
  switch (e.typeKey) {
    case 'drift':
      return `${p.length} mm · ⌀${p.bore * 2} mm`;
    case 'aperture':
      return `${p.voltage} V · ⌀${p.bore * 2} mm`;
    case 'einzel':
      return `${p.voltage} V · ⌀${p.boreRadius * 2} mm`;
    case 'quadrupole':
      return p.rfAmplitude === 0
        ? `DC ${p.dcVoltage} V · ${p.length} mm`
        : `${p.rfAmplitude} V @ ${p.frequency} MHz · ${p.length} mm`;
    case 'bender':
      return `90° ${p.bendPlane === 0 ? 'horizontal' : 'vertical'} · ±${p.voltage.toFixed(
        0
      )} V · r₀ ${p.apertureRadius} mm`;
    default:
      return '';
  }
}

function renderTrack() {
  const idx = selectedIndex();
  const chips = beamline.elements
    .map((e, i) => {
      const last = beamline.elements.length - 1;
      const off = Math.hypot(e.align.dx, e.align.dy);
      return `
        <li class="chip ${i === idx ? 'sel' : ''} chip-${e.typeKey}">
          <button class="chip-body" data-act="select" data-index="${i}">
            <span class="chip-name">${escapeHtml(e.label)}${
              off > 0 ? ' <span class="nudged" title="Misaligned">off</span>' : ''
            }</span>
            <span class="chip-meta">${escapeHtml(summarise(e))}</span>
            <span class="chip-len">${mToMm(e.length).toFixed(0)} mm</span>
          </button>
          <span class="chip-tools">
            <button data-act="left" data-index="${i}" ${i === 0 ? 'disabled' : ''}
              title="Move upstream" aria-label="Move ${escapeHtml(e.label)} upstream">◀</button>
            <button data-act="right" data-index="${i}" ${i === last ? 'disabled' : ''}
              title="Move downstream" aria-label="Move ${escapeHtml(e.label)} downstream">▶</button>
            <button data-act="remove" data-index="${i}"
              ${beamline.elements.length <= 1 ? 'disabled' : ''}
              title="Remove" aria-label="Remove ${escapeHtml(e.label)}">×</button>
          </span>
        </li>`;
    })
    .join('');

  // The source is part of the column as far as selection goes.
  trackEl.innerHTML =
    `<li class="chip chip-source ${selection?.kind === 'source' ? 'sel' : ''}">
       <button class="chip-body" data-act="source">
         <span class="chip-name">Ion source</span>
         <span class="chip-meta">${escapeHtml(summariseBeam())}</span>
         <span class="chip-len">entrance</span>
       </button>
     </li>` + chips;

  autoAlignBtn.hidden = !beamline.misaligned;
}

function summariseBeam() {
  return (
    `${readNumber(inputs.mass, 100)} u · ${readNumber(inputs.charge, 1)}+ · ` +
    `${readNumber(inputs.energy, 50)} eV`
  );
}

/**
 * The placement toolbar.
 *
 * Generated from the element registry, so a new kind of optic appears here
 * the moment it is registered - there is no list of element types anywhere in
 * the UI to forget to update. Each button both clicks to append and drags to
 * place at a chosen point on the beamline.
 */
function renderTools() {
  toolsEl.innerHTML = Object.entries(ELEMENT_TYPES)
    .map(
      ([type, spec]) => `
        <button class="tool" data-add="${type}" draggable="true"
                title="${escapeHtml(spec.blurb)} — click to append, or drag onto the beamline">
          <svg viewBox="0 0 32 18" aria-hidden="true">${spec.icon ?? ''}</svg>
          <span>${escapeHtml(spec.label)}</span>
        </button>`
    )
    .join('');
}

/**
 * Show exactly the panel that belongs to the current selection.
 *
 * Nothing selected: the add menu and the Fly button, and nothing else. An
 * element or the source: that object's settings. Panels are hidden rather
 * than destroyed, so the inputs keep their values and the simulation never
 * depends on what is on screen.
 */
function renderInspector() {
  const isSource = selection?.kind === 'source';
  const e = beamline.elements[selectedIndex()];

  addPanel.hidden = selection !== null;
  beamPanel.hidden = !isSource;
  inspectorEl.hidden = selection === null;

  if (isSource) {
    inspectorEl.innerHTML = `
      <h2>Ion source ${backButton()}</h2>
      <p class="hint">
        Ions fill a <strong>disc</strong>, not a line. A beam launched along one
        axis would sit on a symmetry plane of every element and stay there — the
        motion would look flat because the source was, not because the physics
        is.
      </p>`;
    return;
  }

  if (!e) {
    inspectorEl.hidden = true;
    return;
  }

  const spec = ELEMENT_TYPES[e.typeKey];
  const ion = beamSpec();
  const rows = spec.fields
    .map((f) => {
      const value = e.params[f.key];
      // Typed, not dragged. A deflector's voltage spans four orders of
      // magnitude across the beams this simulator handles, and no slider
      // serves both a forty-volt setting and a forty-kilovolt one - on a
      // range wide enough for the second, the first is inside a pixel. The
      // arrows still step by a sensible amount for the beam in front of it.
      const r = fieldRange(f, e.params, ion);
      const instant = f.rebuild
        ? ''
        : '<span class="instant" title="No re-solve needed">fast</span>';
      return `
        <label class="field field-typed">
          <span class="field-label">
            ${escapeHtml(f.label)}${instant}
            <span class="unit">${escapeHtml(f.unit ?? '')}</span>
          </span>
          <input type="number" data-param="${f.key}" data-rebuild="${f.rebuild ? 1 : 0}"
                 min="${f.min}" max="${f.max}" step="${r.step}" value="${value}" />
          ${f.help ? `<span class="field-help">${escapeHtml(f.help)}</span>` : ''}
        </label>`;
    })
    .join('');

  inspectorEl.innerHTML = `
    <h2>${escapeHtml(e.label)} ${backButton()}</h2>
    <p class="hint">${escapeHtml(spec.blurb)}</p>
    ${rows}
    ${benderReadout(e)}
    ${quadrupoleReadout(e)}
    ${tuneRow(e)}
    ${alignmentRows(e)}
    <div class="row-actions">
      <button data-act="left" data-index="${selectedIndex()}">◀ Upstream</button>
      <button data-act="right" data-index="${selectedIndex()}">Downstream ▶</button>
      <button data-act="remove" data-index="${selectedIndex()}">Remove</button>
    </div>`;
}

/**
 * Tune this element alone.
 *
 * Offered only where there is something to tune. Tuning one element is not the
 * same as tuning the column - a deflector set in isolation is set for the beam
 * the elements before it happen to deliver - so the toolbar button that moves
 * every voltage together is named right next to it.
 */
function tuneRow(e) {
  const keys = TUNABLE[e.typeKey];
  if (!keys?.length) return '';
  const what = keys.length === 1 ? 'this voltage' : 'these voltages';
  return `
    <div class="row-actions tune-row">
      <button data-act="tune" data-index="${selectedIndex()}">Tune ${escapeHtml(
        e.label
      )}</button>
      <span class="field-help">
        Searches ${what} for the best transmission, leaving everything else
        alone. Optimise voltages, in the toolbar, moves the whole column at once.
      </span>
    </div>`;
}

const backButton = () =>
  '<button class="back" data-act="deselect" title="Deselect">✕</button>';

/**
 * Misalignment controls.
 *
 * A real beamline is never perfectly aligned, and seeing what a fraction of a
 * millimetre does to transmission is worth more than any amount of prose
 * about it. Offsets move this element only - its neighbours stay on their own
 * mounts, because that is how they are actually bolted down.
 */
function alignmentRows(e) {
  const mm = (v) => (v * 1e3).toFixed(2);
  return `
    <details class="sub" ${
      e.align.dx || e.align.dy ? 'open' : ''
    }>
      <summary>Alignment</summary>
      <label class="field field-typed">
        <span class="field-label">Offset x<span class="unit">mm</span></span>
        <input type="number" data-align="dx" min="-3" max="3" step="0.05"
               value="${mm(e.align.dx)}" />
      </label>
      <label class="field field-typed">
        <span class="field-label">Offset y<span class="unit">mm</span></span>
        <input type="number" data-align="dy" min="-3" max="3" step="0.05"
               value="${mm(e.align.dy)}" />
        <span class="field-help">
          Moves this element only. The ones after it stay on their own mounts,
          because each is bolted down independently.
        </span>
      </label>
    </details>`;
}

/**
 * Matched electrode voltage for the current ion, beside the control that sets
 * it.
 *
 * The tolerance here is measured, not assumed. A small beam through a
 * deflector with r0/a = 0.95 is fully transmitted from 0.9 to 1.1 times the
 * ideal matched value and dies outside that, so the window is about a tenth
 * either way - wide enough that a figure "off the matched value" is not by
 * itself bad news, and narrow enough to be worth showing.
 */
const BENDER_WINDOW = 0.1;

function benderReadout(e) {
  if (e.typeKey !== 'bender') return '';
  const V = e.matchedVoltage(
    readNumber(inputs.energy, 50),
    Math.abs(readNumber(inputs.charge, 1)) || 1
  );
  const set = e.params.voltage;
  const off = V === 0 ? 0 : Math.abs((set - V) / V);
  return `
    <div class="mathieu ${off <= BENDER_WINDOW ? 'ok' : 'bad'}">
      <span>matched = ${V.toFixed(0)} V</span>
      <span>set = ${set.toFixed(0)} V</span>
      <span class="verdict">${
        off <= BENDER_WINDOW
          ? `${(off * 100).toFixed(0)} % off the ideal value — inside the transmitting window`
          : `${(off * 100).toFixed(0)} % off — outside the window, the beam lands on an electrode`
      }</span>
      <button class="mini" data-act="match">Use matched voltage</button>
    </div>`;
}

/**
 * Mathieu parameters for the selected quadrupole and the current ion.
 *
 * These, not the voltages, are what decide whether an ion is transmitted, so
 * they belong next to the controls that set them.
 */
/**
 * Below about this many RF cycles in the rods, stability stops meaning much.
 *
 * Mathieu stability is an asymptotic property of the equation: it says where
 * the motion stays bounded for ever, not what happens to an ion that crosses
 * in a handful of periods. Such an ion can be thrown out with a perfectly
 * respectable (a, q). This project has been caught by it once already — a
 * 70 mm filter at 1.2 MHz gave 8.6 cycles and lost eight ions in nine, with
 * nothing in the stability numbers to show for it — so the cycle count is
 * shown beside them rather than left to be inferred.
 */
const RF_CYCLES_MIN = 15;

function quadrupoleReadout(e) {
  if (e.typeKey !== 'quadrupole') return '';
  const mass = readNumber(inputs.mass, 100);
  const charge = readNumber(inputs.charge, 1);
  const { a, q } = e.mathieu(mass, Math.abs(charge) || 1);
  const stable = Math.abs(q) < MATHIEU_Q_LIMIT && Math.abs(a) < 0.237;

  // Cycles the ion sees crossing the rods, at its entrance speed.
  const speed = Math.sqrt(
    (2 * readNumber(inputs.energy, 50) * ELEMENTARY_CHARGE * (Math.abs(charge) || 1)) /
      (mass * ATOMIC_MASS_UNIT)
  );
  const cycles = speed > 0 ? (e.length / speed) * e.params.frequency * 1e6 : 0;
  const brief = cycles < RF_CYCLES_MIN;

  return `
    <div class="mathieu ${stable ? 'ok' : 'bad'}">
      <span>a = ${a.toFixed(4)}</span>
      <span>q = ${q.toFixed(4)}</span>
      <span class="verdict">${
        stable
          ? 'inside the first stability region'
          : `outside it — q limit is ${MATHIEU_Q_LIMIT}`
      }</span>
    </div>
    <div class="mathieu ${brief ? 'bad' : 'ok'}">
      <span>${cycles.toFixed(1)} RF cycles</span>
      <span class="verdict">${
        brief
          ? 'too few to rely on — stability is asymptotic, and a short crossing ' +
            'loses ions whatever (a, q) says. Lengthen the rods or raise the frequency.'
          : 'long enough for the stability numbers above to mean something'
      }</span>
    </div>`;
}

/* ------------------------------------------------------------------ */
/* flight                                                              */
/* ------------------------------------------------------------------ */

function markStale() {
  if (trajectories.length === 0) return;
  stale = true;
  flyButton.classList.add('stale');
}

/** The ion the source is set to produce, without its spatial distribution. */
function beamSpec() {
  return {
    mass: readNumber(inputs.mass, 100),
    charge: readNumber(inputs.charge, 1),
    energy: readNumber(inputs.energy, 50),
    z: 0.2,
  };
}

/**
 * A fresh beam matching the source settings.
 *
 * Fresh every call, deliberately: ions carry their own state and the
 * integrator mutates them, so handing the same array to two flights would fly
 * the second one from wherever the first ended up. The optimiser calls this
 * hundreds of times.
 */
function makeBeam(count = Math.round(readNumber(inputs.rays, 9))) {
  return discBeam({
    ...beamSpec(),
    count,
    radius: readNumber(inputs.beamRadius, 1.5),
    divergence: readNumber(inputs.divergence, 0),
  });
}

function startFlight() {
  const spec = beamSpec();
  const count = Math.round(readNumber(inputs.rays, 9));

  let ions;
  try {
    // A DISC, not a line. Placing every ion on the x axis would put the whole
    // beam on a symmetry plane of every element here, where it would stay for
    // ever - the motion would look two-dimensional because the source was,
    // not because the physics is.
    ions = makeBeam(count);
  } catch (err) {
    trajectories = [];
    flight = null;
    stats.error = err.message;
    return;
  }
  stats.error = null;
  stats.flown = false;

  const repulsion = inputs.repulsion.value;
  const beamCurrent = readNumber(inputs.beamCurrent, 0) * 1e-6;
  const ionsPerParticle = 10 ** readNumber(inputs.ionsPerParticle, 6);

  stats.repulsion = repulsion;
  stats.particles = ions.length;

  const opts = {
    method: inputs.method.value,
    cfl: readNumber(inputs.cfl, 0.05),
    repulsion,
    beamCurrent,
    ionsPerParticle,
    recordEvery: 4,
    maxSteps: 400000,
  };

  flightSpec = spec;
  flightOpts = opts;
  flightStarted = performance.now();
  flight = createFlight(beamline, ions, opts);
  trajectories = flight.tracks;
  for (let i = 0; i < trajectories.length; i++) {
    trajectories[i].start = ions[i];
    trajectories[i].focus = null;
  }
}

function finishFlight() {
  stats.flyMs = performance.now() - flightStarted;
  for (const t of trajectories) t.focus = focalCrossing(t.points);

  const through = trajectories.filter((t) => t.stop === 'exited');
  stats.transmitted = through.length;
  stats.reflected = trajectories.filter((t) => t.stop === 'reflected').length;
  stats.struck = trajectories.filter((t) => t.stop === 'electrode').length;
  stats.total = trajectories.length;
  stats.drift = Math.max(0, ...trajectories.map((t) => t.energyDrift));
  stats.steps = flight?.steps ?? 0;

  // Where the surviving beam ends up, which is the number a user of a column
  // actually wants. A focal length is only meaningful for a single lens.
  if (through.length) {
    let worst = 0;
    let atIn = 0;
    let radius = 0;
    for (const t of through) {
      const first = t.points[0];
      const last = t.points[t.points.length - 1];
      radius = Math.max(radius, Math.hypot(last.x, last.y ?? 0));
      const kIn = joulesToEV(kineticEnergy(first));
      const kOut = joulesToEV(kineticEnergy(last));
      if (Math.abs(kOut - kIn) > Math.abs(worst)) {
        worst = kOut - kIn;
        atIn = kIn;
      }
    }
    stats.exitRadius = radius;
    stats.keIn = atIn;
    stats.worstWork = worst;
  } else {
    stats.exitRadius = null;
    stats.keIn = null;
    stats.worstWork = null;
  }

  stats.flown = true;
}

function fly() {
  cancelAnimationFrame(animation.frame);
  animation.running = false;

  try {
    startFlight();
  } catch (err) {
    readoutEl.innerHTML = stat('Error', 'failed', '', err.message, true);
    console.error(err);
    render();
    return;
  }

  stale = false;
  flyButton.classList.remove('stale');

  if (!flight) {
    render();
    drawReadout();
    return;
  }

  // Pace the display, not the physics. Getting this estimate wrong changes
  // only how long the animation takes, never where the ions go.
  const estimate = beamline.length / (readNumber(inputs.cfl, 0.05) * beamline.lengthScale);
  animation.perFrame = Math.max(1, Math.round(estimate / 120));
  animation.running = true;

  setFlyLabel('Flying…', 'ions in flight');
  drawReadout();
  tick();
}

/**
 * Advance the flight by a chunk and draw it, once per frame.
 *
 * The ions are genuinely being integrated here, not replayed. The chunk size
 * is a display choice with no effect on the trajectory, because the time step
 * is chosen from each ion's own state rather than from wall-clock time.
 */
function tick() {
  if (!flight) return;

  try {
    flight.advance(animation.perFrame);
  } catch (err) {
    animation.running = false;
    readoutEl.innerHTML = stat('Error', 'failed', '', err.message, true);
    console.error(err);
    return;
  }

  render();

  if (flight.done) {
    animation.running = false;
    finishFlight();
    drawReadout();
    render();
    setFlyLabel('Fly ions', 'launch the beam');
    return;
  }

  animation.frame = requestAnimationFrame(tick);
}

function setFlyLabel(label, hint) {
  flyButton.innerHTML =
    `<span class="fly-label">${label}</span><span class="fly-hint">${hint}</span>`;
}

/* ------------------------------------------------------------------ */
/* voltage tuning                                                      */
/* ------------------------------------------------------------------ */

/**
 * Search electrode voltages for the best transmission.
 *
 * Only voltages: every one of them is a multiplier on a field that was solved
 * when the element was built, so a few hundred trials cost a few hundred beam
 * flights and no solver time at all. The geometry is left alone, which is also
 * the constraint an operator at a real instrument works under.
 *
 * The search runs in the page rather than a worker, so it must hand the event
 * loop back or the button never repaints and the tab appears to hang. It does
 * that in `onProgress`, which `optimizeVoltages` awaits.
 */
let tuning = null;

function setTuneNote(text, bad = false) {
  tuneNote.textContent = text ?? '';
  tuneNote.hidden = !text;
  tuneNote.classList.toggle('warn', Boolean(bad));
}

async function runTuner(knobs, button, what) {
  if (tuning) {
    // A press during a search stops it rather than starting a competing one.
    // Pressing a *different* tune button does the same, so say so - otherwise
    // it looks as though the button did nothing.
    tuning.stop = true;
    if (tuning.button !== button) {
      setTuneNote(`Stopping the search already running; press again to tune ${what}.`);
    }
    return;
  }
  if (knobs.length === 0) {
    setTuneNote('Nothing to tune — this column has no adjustable voltages.', true);
    return;
  }

  cancelAnimationFrame(animation.frame);
  animation.running = false;

  const original = button.textContent;
  tuning = { stop: false, button };
  button.classList.add('busy');
  optimizeBtn.disabled = button !== optimizeBtn;
  flyButton.disabled = true;
  setTuneNote(
    `Tuning ${what} — ${knobs.length} voltage${knobs.length === 1 ? '' : 's'}, ` +
      'every trial flies the whole beam. Press the button again to stop and keep the best so far.'
  );

  let result;
  try {
    result = await optimizeVoltages(beamline, () => makeBeam(), knobs, {
      flight: { cfl: readNumber(inputs.cfl, 0.05), method: inputs.method.value },
      shouldStop: () => tuning.stop,
      onProgress: async (p) => {
        button.textContent = `${Math.round(p.fraction * 100)}% · ${p.transmitted}/${p.count} — cancel`;
        // Yields to the browser: without this the search would run to
        // completion before a single frame was painted.
        await new Promise((r) => setTimeout(r, 0));
      },
    });
  } catch (err) {
    console.error(err);
    setTuneNote(`Tuning failed: ${err.message}`, true);
    return;
  } finally {
    tuning = null;
    button.textContent = original;
    button.classList.remove('busy');
    optimizeBtn.disabled = false;
    flyButton.disabled = false;
  }

  const changed = knobs
    .map((k, i) => ({ k, from: result.start[i], to: result.values[i] }))
    .filter((c) => c.from !== c.to);

  const settings = changed
    .map((c) => `${c.k.label} ${c.from} → ${c.to} ${c.k.unit}`)
    .join('; ');
  const headline =
    `${result.transmitted}/${result.count} through` +
    (result.exitRadius === null
      ? ''
      : `, beam ${mToMm(result.exitRadius).toFixed(2)} mm at the exit`);

  setTuneNote(
    result.cancelled
      ? `Stopped after ${result.evaluations} trials — keeping the best: ${headline}.` +
          (settings ? ` ${settings}` : '')
      : changed.length === 0
        ? `${result.evaluations} trials: ${headline}. Nothing beat the settings already in place.`
        : `${result.evaluations} trials: ${headline}. ${settings}`,
    result.transmitted === 0
  );

  renderInspector();
  renderTrack();
  render();
  // The trajectories on screen were flown at the old voltages, so say so
  // rather than leaving a picture that no longer matches the settings.
  markStale();
  drawReadout();
}

/* ------------------------------------------------------------------ */
/* rendering                                                           */
/* ------------------------------------------------------------------ */

/**
 * Screen transform: a top view of the global x-z plane, fitted to whatever
 * the column actually occupies.
 *
 * Once a bender is in the line there is no axis to lay along the screen, so
 * the view frames the beamline's bounding box instead. Both axes carry the
 * SAME scale - a bent path drawn with unequal scales would show bend angles
 * that are not the bend angles - and the fit is recomputed whenever the
 * column changes shape.
 */
function makeTransform(width, paneHeight, paneTop, plane, scale, bounds, padding = 16) {
  const b = bounds;
  const margin = beamline.radiusLimit * 1.4;
  const minZ = b.minZ - margin;
  // The across-screen coordinate is global x in the top view and global y in
  // the side view; everything else about the two panes is identical.
  const lo = (plane === 'top' ? b.minX : b.minY) - margin;
  const hi = (plane === 'top' ? b.maxX : b.maxY) + margin;
  const axis = plane === 'top' ? 0 : 1;

  const spanZ = Math.max(1e-6, b.maxZ + margin - minZ);
  const spanT = Math.max(1e-6, hi - lo);

  const offX = padding + (width - 2 * padding - spanZ * scale) / 2;
  const offY = paneTop + (paneHeight - spanT * scale) / 2;

  return {
    plane,
    axis,
    top: paneTop,
    height: paneHeight,
    sx: (z) => offX + (z - minZ) * scale,
    sy: (t) => offY + (hi - t) * scale,
    project: (p) => [offX + (p[2] - minZ) * scale, offY + (hi - p[axis]) * scale],
    /** Screen position back to a world point in this pane's plane. */
    unproject: (px, py) => {
      const t = hi - (py - offY) / scale;
      const z = minZ + (px - offX) / scale;
      return axis === 0 ? [t, 0, z] : [0, t, z];
    },
    scale,
  };
}

/**
 * The scale both panes share.
 *
 * They must match, or the same element would be drawn at two sizes and the
 * eye would read a bend angle that is not there. Both axes within a pane also
 * carry that scale, so angles on screen are true angles.
 */
function fitScale(width, paneHeight, bounds, panes, padding = 16) {
  const margin = beamline.radiusLimit * 1.4;
  const spanZ = Math.max(1e-6, bounds.maxZ - bounds.minZ + 2 * margin);
  const spanX = Math.max(1e-6, bounds.maxX - bounds.minX + 2 * margin);
  const spanY = Math.max(1e-6, bounds.maxY - bounds.minY + 2 * margin);
  const spanT = panes === 1 ? spanX : Math.max(spanX, spanY);
  return Math.min((width - 2 * padding) / spanZ, (paneHeight - 2 * padding) / spanT);
}

function drawElementField(e, T) {
  // Only elements with an axisymmetric (z, r) map have something meaningful
  // to paint in this view. A quadrupole's solve lives in the transverse
  // plane and a bender's in its own bend frame; their slices belong in the
  // cross-section inset, not here.
  if (!e.grid || !e.field || e.typeKey === 'quadrupole' || e.typeKey === 'bender') return;

  const { grid, field } = e;
  const { nz, nr } = grid;

  let maxAbs = 0;
  for (let k = 0; k < field.phi.length; k++) {
    const a = Math.abs(field.phi[k]);
    if (a > maxAbs) maxAbs = a;
  }
  if (maxAbs === 0) return;

  const neg = cssRGB('--pot-neg');
  const zero = cssRGB('--pot-zero');
  const pos = cssRGB('--pot-pos');

  const rows = 2 * nr - 1;
  const off = document.createElement('canvas');
  off.width = nz;
  off.height = rows;
  const offCtx = off.getContext('2d');
  const img = offCtx.createImageData(nz, rows);

  for (let row = 0; row < rows; row++) {
    const j = Math.abs(nr - 1 - row);
    for (let i = 0; i < nz; i++) {
      const [r, g, b] = divergingColour(field.phi[j * nz + i] / maxAbs, neg, zero, pos);
      const p = (row * nz + i) * 4;
      img.data[p] = r;
      img.data[p + 1] = g;
      img.data[p + 2] = b;
      img.data[p + 3] = 255;
    }
  }
  offCtx.putImageData(img, 0, 0);

  // Painted in the element's own frame rather than screen-aligned, so a
  // misaligned or post-bend element carries its field map with it.
  ctx.save();
  withElementTransform(e, T, () => {
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(
      off,
      grid.z0 * T.scale,
      -grid.rMax * T.scale,
      grid.zLength * T.scale,
      2 * grid.rMax * T.scale
    );
  });
  ctx.restore();
}

/**
 * Run `body` with the canvas placed in an element's own frame.
 *
 * Inside, local axial runs along screen +x and local transverse along screen
 * -y, at the view's scale - so an element can draw itself as though it were
 * at the origin facing along the axis, wherever it has actually ended up.
 */
function withElementTransform(e, T, body) {
  const [ox, oy] = T.project(e.frame.o);
  const f = forwardOf(e.frame);
  // Screen x is global z; screen y is minus whichever transverse axis this
  // pane shows. The forward direction therefore appears as (fz, -f[axis]).
  const angle = Math.atan2(-f[T.axis], f[2]);
  ctx.translate(ox, oy);
  ctx.rotate(angle);
  body();
}

function drawContours(e, T) {
  if (!e.grid || !e.field || e.typeKey === 'quadrupole' || e.typeKey === 'bender') return;
  const { grid, field } = e;
  const { nz, nr } = grid;
  const phi = field.phi;

  let maxAbs = 0;
  for (let k = 0; k < phi.length; k++) maxAbs = Math.max(maxAbs, Math.abs(phi[k]));
  if (maxAbs === 0) return;

  const levels = [];
  for (let n = 1; n <= 9; n++) {
    levels.push((maxAbs * n) / 10, (-maxAbs * n) / 10);
  }

  ctx.save();
  withElementTransform(e, T, () => {});
  ctx.strokeStyle = cssVar('--gridline');
  ctx.globalAlpha = 0.85;
  ctx.lineWidth = 1;
  ctx.beginPath();

  const gx = (i) => e.grid.zAt(i) * T.scale;
  const gy = (r) => -r * T.scale;

  for (const level of levels) {
    for (let j = 0; j < nr - 1; j++) {
      for (let i = 0; i < nz - 1; i++) {
        const v00 = phi[j * nz + i];
        const v10 = phi[j * nz + i + 1];
        const v01 = phi[(j + 1) * nz + i];
        const v11 = phi[(j + 1) * nz + i + 1];
        const pts = [];
        const cross2 = (a, b, ia, ja, ib, jb) => {
          if ((a - level) * (b - level) >= 0) return;
          const s = (level - a) / (b - a);
          pts.push([ia + (ib - ia) * s, ja + (jb - ja) * s]);
        };
        cross2(v00, v10, i, j, i + 1, j);
        cross2(v10, v11, i + 1, j, i + 1, j + 1);
        cross2(v01, v11, i, j + 1, i + 1, j + 1);
        cross2(v00, v01, i, j, i, j + 1);
        if (pts.length < 2) continue;
        for (let p = 0; p + 1 < pts.length; p += 2) {
          const [a, b] = [pts[p], pts[p + 1]];
          for (const sign of [1, -1]) {
            ctx.moveTo(gx(a[0]), gy(sign * grid.rAt(a[1])));
            ctx.lineTo(gx(b[0]), gy(sign * grid.rAt(b[1])));
          }
        }
      }
    }
  }
  ctx.stroke();
  ctx.restore();
}

/**
 * Electrodes, as polygons in global space.
 *
 * Filled as paths rather than axis-aligned rectangles, because once the
 * column can bend an element's metal is no longer parallel to anything.
 */
function drawElectrodes(T) {
  ctx.save();
  for (const shape of beamline.outline()) {
    ctx.fillStyle = shape.ghost
      ? cssVar('--electrode-ghost')
      : shape.wall
        ? cssVar('--axis')
        : cssVar('--electrode');
    ctx.globalAlpha = shape.ghost ? 0.35 : shape.wall ? 0.5 : 1;
    // Axisymmetric metal looks the same in any plane containing the axis, so
    // the side view draws the quarter-turned copy rather than an edge-on
    // sliver of the x-z one.
    const corners =
      T.plane === 'side' && shape.cornersRolled ? shape.cornersRolled : shape.corners;
    ctx.beginPath();
    corners.forEach((c, i) => {
      const [px, py] = T.project(c);
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    });
    ctx.closePath();
    ctx.fill();
  }
  ctx.restore();
}

/** The reference orbit, which is a curve as soon as a bender is in the line. */
function drawReferencePath(T) {
  const pts = beamline.centreLine();
  if (pts.length < 2) return;
  ctx.save();
  ctx.strokeStyle = cssVar('--axis');
  ctx.lineWidth = 1;
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  pts.forEach((p, i) => {
    const [px, py] = T.project(p);
    if (i === 0) ctx.moveTo(px, py);
    else ctx.lineTo(px, py);
  });
  ctx.stroke();
  ctx.restore();
}

/**
 * Where one element ends and the next begins.
 *
 * Worth showing, because each was solved as a separate problem and the joins
 * are where that approximation lives. Drawn as a short bar across the
 * element's entrance face rather than a full-height line, since the faces are
 * no longer parallel once the column bends.
 */
function drawBoundaries(T) {
  ctx.save();
  ctx.strokeStyle = cssVar('--gridline');
  ctx.lineWidth = 1;
  ctx.setLineDash([2, 3]);
  for (const e of beamline.elements) {
    const r = e.outerRadius * 1.1;
    const a = T.project(toGlobal(e.frame, across(T, r)));
    const b = T.project(toGlobal(e.frame, across(T, -r)));
    ctx.beginPath();
    ctx.moveTo(a[0], a[1]);
    ctx.lineTo(b[0], b[1]);
    ctx.stroke();
  }
  ctx.restore();
}

function drawTrajectories(T) {
  const trajColour = cssVar('--traj');
  const halo = cssVar('--surface-1');
  const dim = stale ? 0.28 : 1;

  for (const pass of ['halo', 'line']) {
    ctx.save();
    ctx.strokeStyle = pass === 'halo' ? halo : trajColour;
    ctx.lineWidth = pass === 'halo' ? 4 : 2;
    ctx.globalAlpha = (pass === 'halo' ? 0.55 : 1) * dim;
    ctx.lineJoin = 'round';
    ctx.lineCap = 'round';
    for (const traj of trajectories) {
      if (traj.points.length < 2) continue;
      ctx.beginPath();
      for (let n = 0; n < traj.points.length; n++) {
        const p = traj.points[n];
        const [px, py] = T.project([p.x, p.y ?? 0, p.z]);
        if (n === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      ctx.stroke();
    }
    ctx.restore();
  }

  if (animation.running) {
    ctx.save();
    ctx.fillStyle = trajColour;
    ctx.strokeStyle = halo;
    ctx.lineWidth = 1.5;
    for (const traj of trajectories) {
      if (!traj.active) continue;
      const p = traj.state;
      const [px, py] = T.project([p.x, p.y ?? 0, p.z]);
      ctx.beginPath();
      ctx.arc(px, py, 3.5, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
    }
    ctx.restore();
  }

  ctx.save();
  ctx.globalAlpha = dim;
  ctx.fillStyle = cssVar('--electrode-edge');
  for (const traj of trajectories) {
    if (traj.stop !== 'electrode' || traj.active) continue;
    const p = traj.points[traj.points.length - 1];
    const [px, py] = T.project([p.x, p.y ?? 0, p.z]);
    ctx.beginPath();
    ctx.arc(px, py, 3, 0, Math.PI * 2);
    ctx.fill();
  }
  ctx.restore();
}

/** Outline around the selected element, so the inspector's subject is obvious. */
function drawSelection(T) {
  const e = beamline.elements[selectedIndex()];
  if (!e) return;
  const r = e.outerRadius * 1.12;
  const end = e.curved ? 0 : e.length;
  const corners = [
    across(T, r, 0),
    across(T, -r, 0),
    across(T, -r, end),
    across(T, r, end),
  ];

  ctx.save();
  ctx.strokeStyle = cssVar('--accent');
  ctx.fillStyle = cssVar('--accent');
  ctx.lineWidth = 2;

  if (e.curved) {
    // A bent element has no rectangle to outline, so trace its own orbit.
    ctx.globalAlpha = 0.65;
    ctx.beginPath();
    for (let k = 0; k <= 24; k++) {
      const [px, py] = T.project(toGlobal(e.frame, e.pathPoint(k / 24)));
      if (k === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    }
    ctx.lineWidth = Math.max(3, e.bore * 2 * T.scale);
    ctx.stroke();
  } else {
    ctx.globalAlpha = drag?.kind === 'element' && drag.from === selectedIndex() ? 0.22 : 0.1;
    ctx.beginPath();
    corners.forEach((c, i) => {
      const [px, py] = T.project(toGlobal(e.frame, c));
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    });
    ctx.closePath();
    ctx.fill();
    ctx.globalAlpha = 0.55;
    ctx.stroke();
  }
  ctx.restore();
}

/**
 * The ion source: draggable handles at the entrance showing the beam radius.
 *
 * Placed in the first element's frame, so it follows the column's entrance
 * wherever that is.
 */
function drawSource(T) {
  const first = beamline.elements[0];
  if (!first) return;
  const radius = mmToM(readNumber(inputs.beamRadius, 1));

  ctx.save();
  ctx.strokeStyle = cssVar('--traj');
  ctx.fillStyle = cssVar('--traj');
  ctx.lineWidth = 2;

  const at = (tx, tz) => T.project(toGlobal(first.frame, across(T, tx, tz)));
  const lead = Math.max(mmToM(2), radius * 0.8);

  ctx.beginPath();
  let p = at(radius, lead);
  ctx.moveTo(p[0], p[1]);
  p = at(radius, 0);
  ctx.lineTo(p[0], p[1]);
  p = at(-radius, 0);
  ctx.lineTo(p[0], p[1]);
  p = at(-radius, lead);
  ctx.lineTo(p[0], p[1]);
  ctx.stroke();

  for (const sign of [-1, 1]) {
    const [hx, hy] = at(sign * radius, 0);
    ctx.beginPath();
    ctx.arc(hx, hy, 4.5, 0, Math.PI * 2);
    ctx.fill();
  }
  ctx.restore();
}

/** Where a dragged or dropped element would land. */
function drawDropIndicator(T) {
  if (!drag || drag.dropIndex == null) return;
  if (drag.kind === 'element' && !drag.moved) return;

  const i = drag.dropIndex;
  const at =
    i >= beamline.elements.length
      ? { frame: beamline.exitFrame, r: beamline.radiusLimit }
      : { frame: beamline.elements[i].frame, r: beamline.elements[i].outerRadius };

  const a = T.project(toGlobal(at.frame, across(T, at.r * 1.4)));
  const b = T.project(toGlobal(at.frame, across(T, -at.r * 1.4)));

  ctx.save();
  ctx.strokeStyle = cssVar('--accent');
  ctx.lineWidth = 3;
  ctx.setLineDash([6, 4]);
  ctx.beginPath();
  ctx.moveTo(a[0], a[1]);
  ctx.lineTo(b[0], b[1]);
  ctx.stroke();
  ctx.restore();
}

/** A scale bar, since the view now fits itself rather than using a fixed span. */
function drawScale(T, width, height) {
  // Choose a round length that comes out a sensible size on screen.
  const targetPx = width * 0.18;
  const targetMm = mToMm(targetPx / T.scale);
  const nice = [1, 2, 5, 10, 20, 50, 100, 200, 500];
  const barMm = nice.reduce((a, b) => (Math.abs(b - targetMm) < Math.abs(a - targetMm) ? b : a));
  const barPx = mmToM(barMm) * T.scale;

  ctx.save();
  ctx.fillStyle = cssVar('--text-muted');
  ctx.strokeStyle = cssVar('--text-muted');
  ctx.lineWidth = 1;
  ctx.font = '11px ui-monospace, monospace';
  ctx.textBaseline = 'bottom';
  ctx.textAlign = 'center';

  const y = height - 14;
  const x0 = 14;
  ctx.beginPath();
  ctx.moveTo(x0, y);
  ctx.lineTo(x0 + barPx, y);
  ctx.moveTo(x0, y - 4);
  ctx.lineTo(x0, y + 4);
  ctx.moveTo(x0 + barPx, y - 4);
  ctx.lineTo(x0 + barPx, y + 4);
  ctx.stroke();
  ctx.fillText(`${barMm} mm`, x0 + barPx / 2, y - 5);
  ctx.restore();
}

/**
 * Where an ion sits in the transverse plane of the element it is passing
 * through, which is what a profile means once the column can turn.
 *
 * Plotting the global x and y instead only works while the beam travels along
 * z. After a right-hand bend the beam runs along -x, so global x is now the
 * direction of travel and global z is transverse - a profile drawn in (x, y)
 * shows the beam smeared across the screen in proportion to how far it has
 * flown, which is not a beam profile at all.
 *
 * An ion past the end of the column is measured against the exit frame, the
 * natural continuation of the last element's axis.
 */
function transverseAt(p) {
  const g = [p.x, p.y ?? 0, p.z];
  const hit = beamline.locate(g);
  const l = toLocal(hit ? hit.element.frame : beamline.exitFrame, g);
  return { x: l[0], y: l[1], element: hit?.element ?? null };
}

/**
 * Beam cross-section, looking down the axis.
 *
 * The main view is the x-z plane, which cannot show that the ions now move in
 * three dimensions. In a quadrupole they emphatically do: the field converges
 * in one transverse plane while diverging in the other, so a beam that looks
 * well behaved from the side can be being pulled into a line seen end-on.
 *
 * Every ion is drawn against the axis of the element it is in, so the profile
 * keeps meaning what it says after the beam turns a corner.
 */
function drawCrossSection() {
  const size = 132;
  const dpr = window.devicePixelRatio || 1;
  cross.width = Math.round(size * dpr);
  cross.height = Math.round(size * dpr);
  cross.style.width = `${size}px`;
  cross.style.height = `${size}px`;
  crossCtx.setTransform(dpr, 0, 0, dpr, 0, 0);

  crossCtx.clearRect(0, 0, size, size);
  crossCtx.fillStyle = cssVar('--surface-1');
  crossCtx.fillRect(0, 0, size, size);

  // Each ion against the axis of the element it is in. Computed once: locating
  // a point is a walk over the column, and this runs every frame.
  const spots = [];
  for (const traj of trajectories) {
    const p = traj.active ? traj.state : traj.points[traj.points.length - 1];
    if (!p) continue;
    spots.push({ ...transverseAt(p), live: Boolean(traj.active) });
  }

  // The aperture the leading ion is inside. Found by locating it in space
  // rather than by its z, which stops being a path coordinate once the column
  // bends.
  const leader = spots.find((s) => s.live) ?? spots[0];
  const elementHere = leader?.element ?? null;

  /*
    Framed on the aperture the beam is in, not on the widest in the line. A
    deflector's bore is nineteen millimetres against a drift's five, so a
    single fixed scale would draw the beam as a dot for most of its flight.
    The frame never crops the beam, so an ion outside its aperture is still
    visible - that is exactly the moment worth seeing.
  */
  const reach = Math.max(0, ...spots.map((s) => Math.hypot(s.x, s.y)));
  const bore = elementHere?.bore ?? Math.max(...beamline.elements.map((e) => e.bore));
  const limit = Math.max(bore, reach) * 1.15;
  const c = size / 2;
  const k = (size / 2 - 8) / limit;

  crossCtx.save();
  crossCtx.strokeStyle = cssVar('--gridline');
  crossCtx.lineWidth = 1;
  crossCtx.beginPath();
  crossCtx.moveTo(8, c);
  crossCtx.lineTo(size - 8, c);
  crossCtx.moveTo(c, 8);
  crossCtx.lineTo(c, size - 8);
  crossCtx.stroke();

  if (elementHere) {
    crossCtx.strokeStyle = cssVar('--electrode');
    crossCtx.globalAlpha = 0.6;
    crossCtx.beginPath();
    crossCtx.arc(c, c, elementHere.bore * k, 0, Math.PI * 2);
    crossCtx.stroke();
    crossCtx.globalAlpha = 1;
  }
  crossCtx.restore();

  crossCtx.save();
  crossCtx.fillStyle = cssVar('--traj');
  crossCtx.strokeStyle = cssVar('--surface-1');
  crossCtx.lineWidth = 1.2;
  for (const s of spots) {
    crossCtx.globalAlpha = s.live ? 1 : 0.35;
    crossCtx.beginPath();
    crossCtx.arc(c + s.x * k, c - s.y * k, 3, 0, Math.PI * 2);
    crossCtx.fill();
    crossCtx.stroke();
  }
  crossCtx.restore();

  // The scale is stated, because it is no longer fixed.
  crossCtx.fillStyle = cssVar('--text-muted');
  crossCtx.font = '10px ui-monospace, monospace';
  crossCtx.textAlign = 'center';
  crossCtx.fillText(
    elementHere
      ? `${elementHere.label} · ±${mToMm(limit).toFixed(0)} mm`
      : 'beam cross-section',
    c,
    size - 4
  );
}

function render() {
  if (!beamline || beamline.elements.length === 0) return;

  const cssWidth = canvas.parentElement.clientWidth;
  const bounds = beamline.bounds();
  const margin = beamline.radiusLimit * 1.4;

  // A side elevation appears only when the column actually leaves the
  // horizontal plane. For a straight or horizontally-bent line the top view
  // says everything, and a second empty pane would be wasted space.
  const panes = beamline.usesVerticalPlane ? 2 : 1;
  const spanZ = bounds.maxZ - bounds.minZ + 2 * margin;
  const spanT =
    panes === 1
      ? bounds.maxX - bounds.minX + 2 * margin
      : Math.max(
          bounds.maxX - bounds.minX + 2 * margin,
          bounds.maxY - bounds.minY + 2 * margin
        );
  const paneHeight = Math.max(
    150,
    Math.min(320, Math.round((cssWidth * spanT) / Math.max(1e-9, spanZ)))
  );
  const cssHeight = paneHeight * panes;

  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.round(cssWidth * dpr);
  canvas.height = Math.round(cssHeight * dpr);
  canvas.style.height = `${cssHeight}px`;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

  ctx.clearRect(0, 0, cssWidth, cssHeight);
  ctx.fillStyle = cssVar('--surface-1');
  ctx.fillRect(0, 0, cssWidth, cssHeight);

  const scale = fitScale(cssWidth, paneHeight, bounds, panes);
  const planes = panes === 1 ? ['top'] : ['top', 'side'];
  views = planes.map((plane, i) =>
    makeTransform(cssWidth, paneHeight, i * paneHeight, plane, scale, bounds)
  );
  view = views[0];

  for (const T of views) drawPane(T, cssWidth);

  if (panes === 2) {
    ctx.save();
    ctx.strokeStyle = cssVar('--border');
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, paneHeight + 0.5);
    ctx.lineTo(cssWidth, paneHeight + 0.5);
    ctx.stroke();
    ctx.restore();
  }

  drawScale(views[0], cssWidth, cssHeight);
  drawCrossSection();

  scaleNote.hidden = false;
  scaleNote.textContent = beamline.misaligned
    ? 'True scale · elements are misaligned'
    : 'True scale — both panes share one scale';
  scaleNote.classList.toggle('warn', beamline.misaligned);
}

/**
 * A local point displaced across the beam, in whichever transverse direction
 * this pane shows.
 *
 * Handles, selection outlines and element boundaries are all drawn "across
 * the beam", and which axis that means depends on the pane. Going through one
 * helper keeps every such marker honest in both views.
 */
function across(T, t, along = 0) {
  return T.axis === 0 ? [t, 0, along] : [0, t, along];
}

/** Everything that belongs in one projection. */
function drawPane(T, width) {
  ctx.save();
  ctx.beginPath();
  ctx.rect(0, T.top, width, T.height);
  ctx.clip();

  if (inputs.showField.checked) {
    for (const e of beamline.elements) drawElementField(e, T);
  }
  if (inputs.showContours.checked) {
    for (const e of beamline.elements) drawContours(e, T);
  }
  drawBoundaries(T);
  drawReferencePath(T);
  drawElectrodes(T);
  drawSelection(T);
  drawTrajectories(T);
  drawSource(T);
  drawDropIndicator(T);

  // Say which plane this is, since the two look alike for a straight column.
  ctx.fillStyle = cssVar('--text-muted');
  ctx.font = '10px ui-monospace, monospace';
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillText(T.plane === 'top' ? 'top  (x–z)' : 'side (y–z)', 8, T.top + 6);
  ctx.restore();
}

/* ------------------------------------------------------------------ */
/* readout                                                             */
/* ------------------------------------------------------------------ */

function stat(label, value, suffix = '', note = '', warn = false) {
  return `
    <div class="stat">
      <span class="stat-label">${label}</span>
      <span class="stat-value${warn ? ' warn' : ''}">${value}${
        suffix ? `<span class="suffix">${suffix}</span>` : ''
      }</span>
      ${note ? `<span class="stat-note">${note}</span>` : ''}
    </div>`;
}

function drawReadout() {
  if (stats.error) {
    readoutEl.innerHTML = stat('Beam', 'invalid', '', escapeHtml(stats.error), true);
    return;
  }

  const warnings = beamline.warnings;

  if (!stats.flown) {
    readoutEl.innerHTML =
      stat(
        'Beam',
        animation.running ? 'flying…' : 'not flown',
        '',
        animation.running ? 'integrating the ions now' : 'press Fly to launch the beam'
      ) + warningStats(warnings);
    return;
  }

  const loaded = stats.repulsion && stats.repulsion !== 'none';
  const rf = beamline.shortestPeriod !== null;
  const fate = [
    `${stats.transmitted} through`,
    stats.reflected ? `${stats.reflected} reflected` : null,
    stats.struck ? `${stats.struck} on metal` : null,
  ]
    .filter(Boolean)
    .join(' · ');

  const driftPct = stats.drift * 100;

  readoutEl.innerHTML =
    [
      stat(
        'Transmitted',
        `${stats.transmitted}/${stats.total}`,
        '',
        fate,
        stats.transmitted === 0
      ),
      stats.exitRadius === null
        ? stat('Exit radius', '—', '', 'nothing transmitted')
        : stat(
            'Exit radius',
            mToMm(stats.exitRadius).toFixed(2),
            'mm',
            'largest surviving ion at the end of the column'
          ),
      stats.keIn === null
        ? stat('Net work', '—', '', 'nothing transmitted')
        : stat(
            'Net work',
            stats.worstWork.toFixed(3),
            'eV',
            loaded || rf
              ? `on ${stats.keIn.toFixed(0)} eV in · expected: the field does work on the ion`
              : `worst of ${stats.transmitted} transmitted, on ${stats.keIn.toFixed(0)} eV in`,
            !loaded && !rf && Math.abs(stats.worstWork) > 0.02 * stats.keIn
          ),
      // Energy is a constant of the motion only in a static field with no
      // ion-ion interaction. An RF element does work on the ion by design -
      // that is how a quadrupole confines it at all - so a large figure there
      // is the physics, not the integrator, and flagging it would be
      // flagging correct behaviour. It is only a numerical diagnostic when
      // the column is static and the beam non-interacting.
      stat(
        'Energy drift',
        driftPct < 0.01 ? '<0.01' : driftPct.toFixed(2),
        '%',
        rf
          ? 'expected · a time-dependent field does work, so energy is not conserved'
          : loaded
            ? 'expected · ion–ion potential energy is not counted in ½mv² + qφ'
            : 'worst ½mv² + qφ deviation · grid quality, not step size',
        !rf && !loaded && stats.drift > 0.05
      ),
      stat(
        'Column',
        mToMm(beamline.length).toFixed(0),
        'mm',
        `${beamline.elements.length} elements · ${stats.steps.toLocaleString()} steps · ` +
          `${stats.flyMs.toFixed(0)} ms`
      ),
    ].join('') + warningStats(warnings);
}

function warningStats(warnings) {
  return warnings
    .slice(0, 3)
    .map((w) => stat('Model warning', '!', '', escapeHtml(w), true))
    .join('');
}

/* ------------------------------------------------------------------ */
/* wiring                                                              */
/* ------------------------------------------------------------------ */

/**
 * Readouts beside an input.
 *
 * Only one survives. Every other control now holds its own number in a box
 * you can type into, so repeating it alongside would just be the same figure
 * twice. This one is different: what is typed is an exponent, and the count it
 * stands for is the thing worth reading.
 */
const OUTPUTS = {
  ionsPerParticle: (v) => {
    const n = 10 ** parseFloat(v);
    return n < 10 ? n.toFixed(1) : n.toExponential(1).replace('e+', 'e');
  },
};

function syncOutputs() {
  for (const [id, fmt] of Object.entries(OUTPUTS)) {
    const out = el(`${id}Out`);
    if (out) out.textContent = fmt(inputs[id].value);
  }
  const model = inputs.repulsion.value;
  for (const node of document.querySelectorAll('[data-model]')) {
    node.hidden = node.dataset.model !== model;
  }
}

for (const [id, input] of Object.entries(inputs)) {
  const displayOnly = DISPLAY_INPUTS.includes(id);
  const event = input.type === 'range' ? 'input' : 'change';
  input.addEventListener(event, () => {
    syncOutputs();
    if (!displayOnly) markStale();
    // Mathieu numbers and the bender's matched voltage depend on the ion, so
    // the inspector follows the beam.
    if (id === 'mass' || id === 'charge' || id === 'energy') {
      renderTrack();
      renderInspector();
    }
    render();
    if (!displayOnly) drawReadout();
  });
}

// Escape deselects, which collapses the panel back to the add menu.
window.addEventListener('keydown', (e) => {
  if (e.key === 'Escape' && selection !== null) select(null);
});

/** Actions that can appear in either the track or the inspector. */
function handleAction(act, index) {
  switch (act) {
    case 'select':
      select({ kind: 'element', index });
      return true;
    case 'source':
      select({ kind: 'source' });
      return true;
    case 'deselect':
      select(null);
      return true;
    case 'left':
      moveElement(index, -1);
      return true;
    case 'right':
      moveElement(index, 1);
      return true;
    case 'remove':
      removeElement(index);
      return true;
    case 'match': {
      // Put the bender on the matched voltage for the current ion, snapped to
      // the slider's own step so the control shows exactly what was set.
      const e = beamline.elements[selectedIndex()];
      if (e?.typeKey !== 'bender') return true;
      const V = e.matchedVoltage(
        readNumber(inputs.energy, 50),
        Math.abs(readNumber(inputs.charge, 1)) || 1
      );
      const spec = ELEMENT_TYPES.bender.fields.find((f) => f.key === 'voltage');
      const { step } = fieldRange(spec, e.params, beamSpec());
      setParam(selectedIndex(), 'voltage', Math.round(V / step) * step);
      return true;
    }
    case 'tune': {
      const e = beamline.elements[index];
      if (!e) return true;
      const knobs = tunableKnobs(beamline, beamSpec()).filter((k) => k.index === index);
      // The button itself is inside the inspector, which this rerenders on
      // completion, so look it up now rather than holding a stale node.
      const button = inspectorEl.querySelector('button[data-act="tune"]');
      if (button) runTuner(knobs, button, `the ${e.label}`);
      return true;
    }
    default:
      return false;
  }
}

trackEl.addEventListener('click', (e) => {
  const button = e.target.closest('button[data-act]');
  if (button) handleAction(button.dataset.act, Number(button.dataset.index));
});

inspectorEl.addEventListener('click', (e) => {
  const button = e.target.closest('button[data-act]');
  if (button) handleAction(button.dataset.act, Number(button.dataset.index));
});

toolsEl.addEventListener('click', (e) => {
  const button = e.target.closest('button[data-add]');
  if (button) addElement(button.dataset.add);
});

toolsEl.addEventListener('dragstart', (e) => {
  const button = e.target.closest('button[data-add]');
  if (!button) return;
  e.dataTransfer.setData('text/iontrace-element', button.dataset.add);
  e.dataTransfer.effectAllowed = 'copy';
});

/* Dropping a tool onto the beamline places it where it was dropped, rather
   than at the end. The insertion point comes from the closest point on the
   reference orbit, so it works after the column bends. */
canvas.addEventListener('dragover', (e) => {
  e.preventDefault();
  e.dataTransfer.dropEffect = 'copy';
  const { px, py } = pointerPos(e);
  drag = { kind: 'insert', dropIndex: dropIndexAtS(nearestPathDistance(worldAt(px, py))) };
  render();
});

canvas.addEventListener('dragleave', () => {
  if (drag?.kind === 'insert') {
    drag = null;
    render();
  }
});

canvas.addEventListener('drop', (e) => {
  e.preventDefault();
  const type = e.dataTransfer.getData('text/iontrace-element');
  const index = drag?.dropIndex ?? beamline.elements.length;
  drag = null;
  if (!type || !ELEMENT_TYPES[type]) {
    render();
    return;
  }
  statusEl.classList.add('busy');
  beamline.add(createElement(type, startingParams(type, beamSpec())), index);
  statusEl.classList.remove('busy');
  selection = { kind: 'element', index };
  afterStructureChange();
});

autoAlignBtn.addEventListener('click', () => {
  beamline.autoAlign();
  afterStructureChange();
});

optimizeBtn.addEventListener('click', () => {
  runTuner(tunableKnobs(beamline, beamSpec()), optimizeBtn, 'the whole column');
});

/**
 * Read a typed field, or null if it does not yet hold a usable number.
 *
 * A field being typed into passes through states like "", "-" and "1e" on the
 * way to a value. None of those is a voltage, and none should be applied -
 * committing them would clamp the box to a limit mid-keystroke and fight the
 * person typing. Out-of-range values are clamped rather than rejected, so a
 * typo cannot put an element somewhere the sliders could not have.
 */
function readTyped(input) {
  if (input.value.trim() === '') return null;
  const v = parseFloat(input.value);
  if (!Number.isFinite(v)) return null;
  const lo = parseFloat(input.min);
  const hi = parseFloat(input.max);
  return Math.min(Number.isFinite(hi) ? hi : v, Math.max(Number.isFinite(lo) ? lo : v, v));
}

/**
 * Apply what a typed inspector field says.
 *
 * `live` distinguishes a keystroke from a commit. A voltage only rescales
 * stored solutions, so it can follow every keystroke and the picture updates
 * as the number is typed. Geometry needs a fresh Laplace solve, and re-solving
 * once per character while someone types "15" would solve for 1 first - slow,
 * and briefly wrong. Those wait for the field to be committed: blur, Enter, or
 * the stepper arrows.
 */
function applyTypedField(input, live) {
  const key = input.dataset.param;
  if (key) {
    if (live && input.dataset.rebuild === '1') return;
    const value = readTyped(input);
    if (value === null) return;
    setParam(selectedIndex(), key, value);
    return;
  }

  const alignKey = input.dataset.align;
  if (!alignKey) return;
  const mm = readTyped(input);
  if (mm === null) return;
  const element = beamline.elements[selectedIndex()];
  if (!element) return;
  // Misalignment costs no re-solve: the element is unchanged, only where it
  // sits. That is the whole reason it can be dragged smoothly.
  element.align = { ...element.align, [alignKey]: mmToM(mm) };
  beamline.layout();
  markStale();
  renderTrack();
  render();
  drawReadout();
}

// Inspector fields edit the selected element.
inspectorEl.addEventListener('input', (e) => {
  const input = e.target.closest('input[data-param], input[data-align]');
  if (input) applyTypedField(input, true);
});

inspectorEl.addEventListener('change', (e) => {
  const input = e.target.closest('input[data-param], input[data-align]');
  if (!input) return;
  applyTypedField(input, false);
  // Show what was actually taken. A value outside the element's limits is
  // clamped rather than refused, and the box should say so rather than
  // displaying a number the simulation is not using.
  const taken = readTyped(input);
  if (taken !== null && String(taken) !== input.value) input.value = taken;
});

flyButton.addEventListener('click', fly);

window.addEventListener('keydown', (e) => {
  if (e.key !== ' ' && e.key !== 'Enter') return;
  const tag = document.activeElement?.tagName;
  if (tag === 'INPUT' || tag === 'SELECT' || tag === 'BUTTON') return;
  e.preventDefault();
  fly();
});

let resizeTimer;
window.addEventListener('resize', () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(render, 100);
});

window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', render);

/* ------------------------------------------------------------------ */
/* start                                                               */
/* ------------------------------------------------------------------ */

statusEl.classList.add('busy');
beamline = defaultBeamline();
statusEl.classList.remove('busy');
selection = null;

syncOutputs();
renderTools();
renderTrack();
renderInspector();
render();
drawReadout();
fly();
