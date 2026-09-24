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

import { Beamline, exitsOf } from './beamline.js';
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
import { serialise, restore } from './scene.js';
import { topic } from './help.js';
import { LESSONS, allSteps } from './tutorial.js';
import {
  joulesToEV,
  mToMm,
  mmToM,
  ELEMENTARY_CHARGE,
  ATOMIC_MASS_UNIT,
} from './constants.js';
import { NO_ELECTRODE } from './grid.js';
import { toGlobal, toLocal, forwardOf, compose } from './frames.js';

const el = (id) => document.getElementById(id);

const canvas = el('scene');
const ctx = canvas.getContext('2d');
const cross = el('cross');
const crossCtx = cross.getContext('2d');
const statusEl = el('status');
const readoutEl = el('readout');
const inspectorEl = el('inspector');
const addPanel = el('addPanel');
const toolsEl = el('tools');
const beamPanel = el('beamPanel');
const autoAlignBtn = el('autoAlign');
const optimizeBtn = el('optimize');
const tuneNote = el('tuneNote');
const scaleNote = el('scaleNote');
const flyButton = el('fly');
const canvasFrame = el('canvasFrame');
const viewFly = el('viewFly');
const fullscreenBtn = el('fullscreen');
const fringeToggle = el('fringe');
const fringeNote = el('fringeNote');
const repulsionNote = el('repulsionNote');
const flowEl = el('flow');
const readoutHint = el('readoutHint');
const fieldPanel = el('fieldPanel');
const trackHint = el('trackHint');
const saveBtn = el('save');
const loadBtn = el('load');
const loadFile = el('loadFile');

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
  showLines: el('showLines'),
};

/** Controls that change only how the scene is drawn, never the physics. */
const DISPLAY_INPUTS = ['showField', 'showContours', 'showLines'];

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
/** A geometry rebuild waiting on an idle timer. See GEOMETRY_DELAY. */
let geometryTimer = null;

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

  // A column solve holds its own copy of the voltages. Re-applying them is
  // fast adjust over that grid, not a re-solve.
  beamline.syncRuns();

  markStale();
  renderTrack();
  refreshInspector();
  render();
  drawReadout();
}

/**
 * Place an element.
 *
 * `attach` names the exit it hangs from. Without one it goes on the first free
 * exit of whatever is selected, or the end of the line if nothing is - which
 * is what "add this" means when there is only one line to add it to.
 */
function addElement(type, attach = null) {
  const after = selectedIndex();
  const index = after >= 0 ? after + 1 : beamline.elements.length;
  statusEl.classList.add('busy');
  // Matched to the beam that is in the source right now, so a deflector
  // dropped in from the toolbar actually deflects.
  const placed = beamline.add(
    createElement(type, startingParams(type, beamSpec())),
    index,
    attach
  );
  statusEl.classList.remove('busy');
  selection = { kind: 'element', index: beamline.elements.indexOf(placed) };
  afterStructureChange();
}

/** The exit a new element should hang from, given what is selected. */
function defaultAttach() {
  const e = beamline.elements[selectedIndex()];
  if (!e) return null;
  const free = exitsOf(e).find((x) => !beamline.childAt(e, x.port));
  return free ? { parent: e, port: free.port } : null;
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
  describeFringe();
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

/**
 * Is the cursor on a drawn trajectory?
 *
 * Tested in screen pixels rather than in world coordinates, so the target is
 * the same size however far the view is zoomed out - which is what makes a
 * thin line clickable at all. Points are compared directly rather than by
 * distance to each segment: they are recorded densely enough that the gap
 * between them is well under the tolerance.
 */
function onTrajectory(px, py, tol = 7) {
  const T = paneAt(py);
  if (!T) return false;
  for (const t of trajectories) {
    const pts = t.points;
    if (!pts) continue;
    // Every few points is plenty for a hit test and keeps this cheap on a
    // beam of twenty thousand recorded steps.
    const stride = Math.max(1, Math.floor(pts.length / 400));
    for (let i = 0; i < pts.length; i += stride) {
      const p = pts[i];
      const [sx, sy] = T.project([p.x, p.y ?? 0, p.z]);
      if (Math.abs(sx - px) <= tol && Math.abs(sy - py) <= tol) return true;
    }
  }
  return false;
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


function select(next) {
  // A geometry rebuild waiting on a timer belongs to the element that was
  // selected when it was scheduled, so it is dropped rather than applied to
  // whatever is selected now.
  cancelGeometry();
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

  // A click on the beam itself selects the beam. The ions are the most
  // clickable thing on screen and used not to be clickable at all, so their
  // settings could only be reached through the one box at the far left of the
  // beamline panel.
  if (onTrajectory(px, py)) {
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

  // Selection only. Re-hanging an element is done in the flow chart, which is
  // a picture of the structure being edited rather than a picture of hardware
  // drawn to scale - so nothing there overlaps, and what you aim at is what
  // you get.
  select({ kind: 'element', index });
});

canvas.addEventListener('pointermove', (e) => {
  if (!beamline || !view) return;
  const { px, py } = pointerPos(e);

  if (!drag) {
    canvas.style.cursor = onSourceHandle(px, py)
      ? 'grab'
      : elementIndexAtWorld(worldAt(px, py)) >= 0 || onTrajectory(px, py)
        ? 'pointer'
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

});

function endDrag(e) {
  if (!drag) return;
  drag = null;
  canvas.style.cursor = 'default';
  try {
    canvas.releasePointerCapture(e.pointerId);
  } catch {
    /* already released */
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

/**
 * What an element is currently set to, as SHORT lines.
 *
 * Two of them at most, because an SVG `<text>` does not wrap - it runs on out
 * of the box and over whatever is beside it. A deflector's settings came to
 * 33 characters and a funnel's to 35, in a box that holds about 26, so both
 * spilled across the diagram.
 *
 * So the settings are split rather than truncated wherever possible: a clipped
 * line hides a number, and the numbers are the reason the line is there. The
 * clip in `fit` is a backstop for a value nobody anticipated being long, not
 * the normal case.
 */
function summarise(e) {
  const p = e.params;
  switch (e.typeKey) {
    case 'drift':
      return [`⌀${p.bore * 2} mm bore`];
    case 'aperture':
      return [`${p.voltage} V · ⌀${p.bore * 2} mm`];
    case 'einzel':
      return [`${p.voltage} V · ⌀${p.boreRadius * 2} mm`];
    case 'quadrupole':
      return p.rfAmplitude === 0
        ? [`DC ${p.dcVoltage} V`]
        : [`${p.rfAmplitude} V @ ${p.frequency} MHz`, `r₀ ${p.fieldRadius} mm`];
    case 'multipole':
      return [
        `${e.poles} rods · r₀ ${p.fieldRadius} mm`,
        `${p.rfAmplitude} V @ ${p.frequency} MHz`,
      ];
    case 'funnel':
      return [
        `${p.rings} rings · ⌀${p.entryRadius * 2}→${p.exitRadius * 2} mm`,
        `${p.rfAmplitude} V @ ${p.frequency} MHz`,
      ];
    case 'bender':
      return [
        `±${p.voltage.toFixed(0)} V · r₀ ${p.apertureRadius} mm`,
        `90° ${p.bendPlane === 0 ? 'horizontal' : 'vertical'}`,
      ];
    default:
      return [];
  }
}

/** Clip a string to `n` characters, marking that something was cut. */
const fit = (s, n) => (s.length <= n ? s : `${s.slice(0, n - 1)}…`);

/**
 * The beamline panel.
 *
 * One view: the flow chart. A column used to be a list and rendered as a list,
 * but once a deflector can have hardware on three of its exits the shape on
 * screen has to be the shape of the thing - a nested list puts two branches
 * one above the other, which reads as a sequence and is a lie about where the
 * beam can go.
 */
function renderTrack() {
  autoAlignBtn.hidden = !beamline.misaligned;
  renderFlow();
}

/* ------------------------------------------------------------------ */
/* flow chart                                                          */
/* ------------------------------------------------------------------ */

/*
  Flow-chart geometry, in CSS pixels.

  Sized to be read rather than to fit. The panel scrolls in both directions, so
  a long column runs off to the right and a branching one runs off the bottom,
  and neither is shrunk to squeeze into the panel - which would make the chart
  least legible exactly when there is most in it.
*/
const FLOW = { w: 228, h: 66, gapX: 62, gapY: 26, pad: 18 };

/*
  How much text a row holds, in characters.

  Derived from the box, not guessed. The name row runs from the left padding to
  the remove control, 182 px, at a 14 px proportional face whose average glyph
  is near 0.55 em. The settings rows have the full inner width, 200 px, at
  11.5 px monospace where every glyph is 0.6 em.

  Conservative on purpose. An SVG text element neither wraps nor clips; going
  over does not look crowded, it looks like the box has sprung a leak across
  whatever is beside it. Checked against every element at the extremes of every
  field it displays: the longest settings line fills 69 % of its row.
*/
const NAME_CHARS = 23;
const SUB_CHARS = 28;

/**
 * Lay the column out as a tidy tree.
 *
 * Column is depth from the source, so the beam runs left to right. Row is
 * assigned by walking the leaves in order and giving each the next lane; a
 * node with children sits at the mean of theirs, which keeps a junction
 * centred between the two lines leaving it instead of stuck on one of them.
 */
function flowLayout() {
  const pos = new Map();
  const sockets = new Map();
  let lane = 0;

  /*
    Every exit gets a lane of its own, whether or not anything is bolted to it.

    Allocating lanes only to exits that HAVE children was wrong in a way the
    picture made obvious: a deflector has three ways out, so with one branch
    built the other two sockets were both placed at the same half-row offset
    from their parent and drawn on top of each other, three labels in one
    illegible pile. An empty exit takes up room on a real beamline and it takes
    up room here.
  */
  const place = (e, col) => {
    const rows = [];
    for (const exit of exitsOf(e)) {
      const child = beamline.childAt(e, exit.port);
      if (child) {
        rows.push(place(child, col + 1));
      } else {
        const row = lane++;
        sockets.set(`${beamline.elements.indexOf(e)}:${exit.port}`, { col: col + 1, row, exit });
        rows.push(row);
      }
    }
    // Centred on what leaves it, so a junction sits between its branches
    // rather than level with one of them.
    const row = rows.length ? rows.reduce((a, b) => a + b, 0) / rows.length : lane++;
    pos.set(e, { col, row });
    return row;
  };

  for (const r of beamline.roots()) place(r, 1);
  // The source occupies column zero, level with whatever it feeds.
  const first = beamline.roots()[0];
  return {
    pos,
    sockets,
    sourceRow: first ? pos.get(first).row : 0,
    lanes: Math.max(lane, 1),
  };
}

const flowX = (col) => FLOW.pad + col * (FLOW.w + FLOW.gapX);
const flowY = (row) => FLOW.pad + row * (FLOW.h + FLOW.gapY);

/** An elbow from the right edge of one box to the left edge of another. */
function flowLink(x1, y1, x2, y2) {
  const mid = x1 + FLOW.gapX / 2;
  return `M${x1} ${y1} H${mid} V${y2} H${x2}`;
}

function renderFlow() {
  // Keep the scroll position across a redraw. Selecting a box, nudging a
  // voltage or flying the beam all redraw this, and having the panel jump back
  // to the top left every time would make a long column unusable.
  const { scrollLeft, scrollTop } = flowEl;
  const { pos, sockets, sourceRow, lanes } = flowLayout();
  // Wide enough for the sockets too, which sit one column past their parent.
  const cols = Math.max(
    2,
    ...[...pos.values()].map((p) => p.col + 1),
    ...[...sockets.values()].map((s) => s.col + 1)
  );
  const width = flowX(cols) + FLOW.pad;
  const height = flowY(lanes) + FLOW.pad;
  const idx = selectedIndex();

  const parts = [];

  // The source, and its link into the first element.
  const sy = flowY(sourceRow) + FLOW.h / 2;
  parts.push(
    `<g class="node source ${selection?.kind === 'source' ? 'sel' : ''}" data-act="source">
       <rect x="${flowX(0)}" y="${flowY(sourceRow)}" width="${FLOW.w}" height="${FLOW.h}" rx="8"/>
       <text x="${flowX(0) + 12}" y="${flowY(sourceRow) + 19}">Ion source</text>
       <text class="sub" x="${flowX(0) + 12}" y="${flowY(sourceRow) + 34}">${escapeHtml(
         summariseBeam()
       )}</text>
     </g>`
  );
  const root = beamline.roots()[0];
  if (root) {
    const p = pos.get(root);
    parts.push(
      `<path class="link" d="${flowLink(
        flowX(0) + FLOW.w,
        sy,
        flowX(p.col),
        flowY(p.row) + FLOW.h / 2
      )}"/>`
    );
  }

  for (const e of beamline.elements) {
    const p = pos.get(e);
    if (!p) continue;
    const x = flowX(p.col);
    const y = flowY(p.row);
    const i = beamline.elements.indexOf(e);

    for (const exit of exitsOf(e)) {
      const child = beamline.childAt(e, exit.port);
      const junction = exitsOf(e).length > 1;
      if (child) {
        const c = pos.get(child);
        parts.push(
          `<path class="link" d="${flowLink(
            x + FLOW.w,
            y + FLOW.h / 2,
            flowX(c.col),
            flowY(c.row) + FLOW.h / 2
          )}"/>`
        );
        if (junction) {
          parts.push(
            `<text class="port" x="${x + FLOW.w + 6}" y="${
              flowY(c.row) + FLOW.h / 2 - 5
            }">${escapeHtml(exit.label)}</text>`
          );
        }
      } else {
        // An unused exit, on its own lane, drawn as a socket to start a line
        // from. Placed in the column a child would occupy so the chart reads
        // the same whether an exit is filled or not.
        const s = sockets.get(`${i}:${exit.port}`);
        if (!s) continue;
        const ex = flowX(s.col);
        const ey = flowY(s.row) + FLOW.h / 2;
        parts.push(
          `<path class="link open" d="${flowLink(x + FLOW.w, y + FLOW.h / 2, ex, ey)}"/>`,
          `<g class="socket ${
            pendingPort?.parent === e && pendingPort?.port === exit.port ? 'armed' : ''
          }" data-act="port" data-index="${i}" data-port="${exit.port}">
             <rect x="${ex}" y="${ey - 15}" width="${FLOW.w * 0.8}" height="30" rx="15"/>
             <text x="${ex + 14}" y="${ey + 5}">+ ${escapeHtml(
               junction ? exit.label.toLowerCase() : 'add'
             )}</text>
           </g>`
        );
        if (junction) {
          parts.push(
            `<text class="port" x="${x + FLOW.w + 8}" y="${ey - 20}">${escapeHtml(
              exit.label
            )}</text>`
          );
        }
      }
    }

    /*
      Three rows, each one line and none of them wrapping.

      The length is right-aligned on the BOTTOM row rather than beside the
      name. Sharing the name's row cost it 40 px, which is what pushed
      "Quadrupole deflector" over the edge; the settings lines only fill about
      two thirds of their row, so down there it has room to spare and nothing
      to collide with.
    */
    const nudged = Math.hypot(e.align.dx, e.align.dy) > 0;
    const lines = summarise(e).slice(0, 2);
    const rows = [y + 41, y + 57];
    parts.push(
      `<g class="node ${i === idx ? 'sel' : ''} node-${e.typeKey}"
          data-act="select" data-index="${i}">
         <rect class="hit" x="${x}" y="${y}" width="${FLOW.w}" height="${FLOW.h}" rx="9"/>
         <text class="name" x="${x + 14}" y="${y + 23}"
           >${escapeHtml(fit(e.label + (nudged ? ' (off axis)' : ''), NAME_CHARS))}</text>
         ${lines
           .map(
             (t, r) =>
               `<text class="sub" x="${x + 14}" y="${rows[r]}">${escapeHtml(
                 fit(t, SUB_CHARS)
               )}</text>`
           )
           .join('')}
         <text class="len" text-anchor="end" x="${x + FLOW.w - 14}" y="${y + 57}"
           >${mToMm(e.length).toFixed(0)} mm</text>
         <g class="kill" data-act="remove" data-index="${i}">
           <circle cx="${x + FLOW.w - 16}" cy="${y + 16}" r="9"/>
           <text x="${x + FLOW.w - 16}" y="${y + 20}">×</text>
         </g>
       </g>`
    );
  }

  flowEl.innerHTML =
    `<svg viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" role="img">` +
    parts.join('') +
    '</svg>';

  flowEl.scrollLeft = scrollLeft;
  flowEl.scrollTop = scrollTop;
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
  // Both of these are about the beam rather than the column, so they appear
  // with it: the physics switches above the ion settings, the flight readout
  // below them.
  fieldPanel.hidden = !isSource;
  readoutEl.hidden = !isSource;
  readoutHint.hidden = isSource;

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
            ${f.help ? helpButton(f.help, f.label) : ''}
          </span>
          <input type="number" data-param="${f.key}" data-rebuild="${f.rebuild ? 1 : 0}"
                 min="${f.min}" max="${f.max}" step="${r.step}" value="${value}" />
        </label>`;
    })
    .join('');

  inspectorEl.innerHTML = `
    <h2>
      ${escapeHtml(e.label)}
      ${topic(e.typeKey) ? helpTopicButton(e.typeKey, e.label) : ''}
      ${backButton()}
    </h2>
    <p class="hint">${escapeHtml(spec.blurb)}</p>
    ${rows}
    <div id="derived">${bendDirectionRow(e)}${benderReadout(e)}${quadrupoleReadout(e)}${multipoleReadout(
      e
    )}${funnelReadout(e)}</div>
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

/**
 * Update the inspector without rebuilding it.
 *
 * `renderInspector` replaces the panel's HTML, which destroys the very input
 * the user is typing into - the field loses focus and the caret on every
 * keystroke, so nothing longer than one character can be typed. Anything that
 * fires while a field is being edited must come through here instead: it
 * writes to the existing nodes rather than replacing them.
 *
 * The value of a field is deliberately not written back while it has focus.
 * Doing so fights the person typing: "-3" would be clamped to the minimum and
 * rewritten under the caret before the rest of the number arrived.
 */
function refreshInspector() {
  const e = beamline.elements[selectedIndex()];
  if (!e) return;
  const spec = ELEMENT_TYPES[e.typeKey];
  const ion = beamSpec();

  for (const f of spec.fields) {
    const input = inspectorEl.querySelector(`input[data-param="${f.key}"]`);
    if (!input) continue;
    const r = fieldRange(f, e.params, ion);
    input.min = r.min;
    input.max = r.max;
    input.step = r.step;
    if (document.activeElement !== input) input.value = e.params[f.key];
  }

  const derived = inspectorEl.querySelector('#derived');
  if (derived) {
    derived.innerHTML = bendDirectionRow(e) + benderReadout(e) + quadrupoleReadout(e) + multipoleReadout(e) + funnelReadout(e);
  }
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

/**
 * Where a deflector sends the beam.
 *
 * Two things decide that, and asking someone to set them separately is asking
 * them to know the convention: the roll says which plane the bend happens in,
 * and the VOLTAGE says whether it happens at all. Off, the beam goes straight
 * through the box and out the far side. So "send it left" means roll to zero
 * AND put the electrodes on their matched voltage, and "straight on" means
 * zero volts whatever the roll.
 *
 * One click does both. The angle box stays for a bend that is not one of the
 * four square directions, and the voltage box for tuning it by hand.
 *
 * The straight option is listed last and set apart, because it is a different
 * kind of answer: it is the deflector doing nothing rather than doing
 * something in a direction.
 */
/*
  The four directions, as a roll and a polarity.

  Two independent things aim a deflector, and the pairing is what makes the
  buttons work. The ROLL says which plane it bends in - horizontal or vertical
  - and the SIGN of the voltage says which of the two ways within that plane.
  So Left and Right are the same hardware orientation at opposite polarity,
  and they leave through DIFFERENT PORTS: a column can have a line bolted to
  both and switch between them without anything moving. Up and Down are the
  same pair, rolled a quarter turn.
*/
const BEND_DIRECTIONS = [
  { label: 'Left', glyph: '←', deg: 0, sign: 1 },
  { label: 'Down', glyph: '↓', deg: 90, sign: 1 },
  { label: 'Right', glyph: '→', deg: 0, sign: -1 },
  { label: 'Up', glyph: '↑', deg: 90, sign: -1 },
];

/** The compass bearing a roll and a polarity send the beam on. */
const bearingOf = (deg, sign) => ((((deg + (sign < 0 ? 180 : 0)) % 360) + 360) % 360);

/** Where the beam is actually going, given the voltage and the roll. */
function bendState(e) {
  const matched = e.matchedVoltage(
    readNumber(inputs.energy, 50),
    Math.abs(readNumber(inputs.charge, 1)) || 1
  );
  const set = e.params.voltage;
  const bending = matched !== 0 && Math.abs(set / matched) > 0.5;
  return {
    matched,
    bending,
    // Compared as a bearing rather than as a roll, because a roll of 180 at
    // one polarity and a roll of 0 at the other put the beam in exactly the
    // same place - and a button that fails to light up for a setting it would
    // have produced is worse than no button.
    bearing: bearingOf(e.params.bendPlane, set < 0 ? -1 : 1),
  };
}

function bendDirectionRow(e) {
  if (e.typeKey !== 'bender') return '';
  const { bearing, bending } = bendState(e);
  const buttons = BEND_DIRECTIONS.map(
    (d, i) =>
      `<button class="dir ${
        bending && bearingOf(d.deg, d.sign) === bearing ? 'sel' : ''
      }" data-act="send" data-index="${i}"
               title="Send the beam ${d.label.toLowerCase()} — sets the bend plane and the matched voltage of the right polarity">
         <span class="dir-glyph">${d.glyph}</span>${escapeHtml(d.label)}
       </button>`
  ).join('');
  return `
    <div class="dirs">${buttons}</div>
    <button class="dir dir-wide ${bending ? '' : 'sel'}" data-act="send" data-index="-1"
            title="Switch the deflector off — the beam passes straight through">
      <span class="dir-glyph">⇢</span>Straight through
    </button>`;
}

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

/** RF periods an ion of the current beam sees crossing an element. */
function rfCycles(e) {
  const mass = readNumber(inputs.mass, 100);
  const charge = Math.abs(readNumber(inputs.charge, 1)) || 1;
  const speed = Math.sqrt(
    (2 * readNumber(inputs.energy, 50) * ELEMENTARY_CHARGE * charge) /
      (mass * ATOMIC_MASS_UNIT)
  );
  return speed > 0 ? (e.length / speed) * e.params.frequency * 1e6 : 0;
}

/**
 * What governs a multipole guide: the depth of its effective potential well,
 * and whether the approximation that produced it applies.
 */
function multipoleReadout(e) {
  if (e.typeKey !== 'multipole') return '';
  const mass = readNumber(inputs.mass, 100);
  const charge = Math.abs(readNumber(inputs.charge, 1)) || 1;
  const { depth, q, valid } = e.trapping(mass, charge);
  const cycles = rfCycles(e);
  const brief = cycles < RF_CYCLES_MIN;
  return `
    <div class="mathieu ${valid ? 'ok' : 'bad'}">
      <span>well ${depth.toFixed(2)} eV</span>
      <span>q = ${q.toFixed(3)}</span>
      <span class="verdict">${
        valid
          ? 'the drive is fast enough for an effective potential to mean something'
          : `q above 0.3 — the effective potential stops describing the motion, so ` +
            'this depth is not what the ion actually feels'
      }</span>
    </div>
    <div class="mathieu ${brief ? 'bad' : 'ok'}">
      <span>${cycles.toFixed(1)} RF cycles</span>
      <span class="verdict">${
        brief
          ? 'too few to guide on — lengthen the rods or raise the frequency'
          : 'long enough to guide on'
      }</span>
    </div>`;
}

/**
 * What governs an ion funnel here, which is not what governs a real one.
 */
function funnelReadout(e) {
  if (e.typeKey !== 'funnel') return '';
  const mass = readNumber(inputs.mass, 100);
  const charge = Math.abs(readNumber(inputs.charge, 1)) || 1;
  const well = e.wellAt(e.rings - 2, mass, charge);
  const cycles = rfCycles(e);
  const brief = cycles < RF_CYCLES_MIN;
  const push = e.params.dcEntry - e.params.dcExit;
  return `
    <div class="mathieu ${brief ? 'bad' : 'ok'}">
      <span>${cycles.toFixed(1)} RF cycles</span>
      <span>wall ${well.toFixed(1)} eV</span>
      <span>DC push ${push.toFixed(0)} V</span>
      <span class="verdict">${
        brief
          ? 'too few cycles for the RF wall to mean anything — raise the frequency'
          : 'the drive is fast enough for the RF wall to hold'
      }</span>
    </div>
    <div class="mathieu bad">
      <span class="verdict">
        No buffer gas. A real funnel works at 1–30 mbar and relies on collisions
        to damp the ions into the well; without that a <em>deeper</em> wall
        transmits <em>worse</em>, because the heating it causes has nowhere to
        go. Measured here: 1.1 eV of wall passes everything, 17 eV passes
        nothing. So this reproduces a funnel’s field, not its ability to
        collect a warm cloud.
      </span>
    </div>`;
}

/* ------------------------------------------------------------------ */
/* flight                                                              */
/* ------------------------------------------------------------------ */

function markStale() {
  if (trajectories.length === 0) return;
  stale = true;
  flyButton.classList.add('stale');
  viewFly.classList.add('stale');
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
  // Landing on a detector's collecting surface is a strike, but it is the one
  // strike that means the experiment worked, so it is counted separately.
  const struck = trajectories.filter((t) => t.stop === 'electrode');
  stats.detected = struck.filter((t) => {
    const p = t.points[t.points.length - 1];
    return beamline.detected(p.x, p.y ?? 0, p.z);
  }).length;
  stats.struck = struck.length - stats.detected;
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
  viewFly.classList.remove('stale');

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
  viewFly.textContent = label;
}

/* ------------------------------------------------------------------ */
/* full screen                                                         */
/* ------------------------------------------------------------------ */

const isFullscreen = () => document.fullscreenElement === canvasFrame;

function toggleFullscreen() {
  if (isFullscreen()) document.exitFullscreen();
  else canvasFrame.requestFullscreen?.().catch((err) => console.warn(err));
}

/*
  The canvas is sized from its container, and going full screen changes that
  container's size without resizing the window - so the redraw has to be hung
  on the fullscreen event rather than on `resize`.
*/
document.addEventListener('fullscreenchange', () => {
  viewFly.hidden = !isFullscreen();
  fullscreenBtn.textContent = isFullscreen() ? 'Exit full screen' : 'Full screen';
  render();
  drawReadout();
});

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

/**
 * The combinations of knobs the beam turns out not to care about.
 *
 * These come out of the Hessian's near-zero eigenvalues, and they are worth
 * saying out loud: each one is a statement about the instrument rather than
 * about the search. "Raising the lens by this while nudging the deflector by
 * that changes nothing" tells an operator which knob is redundant, and tells
 * anyone reading a tuned setting why it is not unique.
 */
function describeNullSpace(nullSpace) {
  if (!nullSpace?.length) return '';
  const say = (d) =>
    d.direction
      // Only the knobs that carry real weight in this direction.
      .filter((c) => Math.abs(c.weight) > 1e-9)
      .sort((a, b) => Math.abs(b.weight) - Math.abs(a.weight))
      .slice(0, 3)
      .map(
        (c) =>
          `${c.weight >= 0 ? '+' : '−'}${Math.abs(c.weight).toPrecision(3)} ` +
          `${c.label.split(' · ')[0]}`
      )
      .join(' with ');
  return (
    ` The beam is insensitive to ${nullSpace.length === 1 ? 'one combination' : `${nullSpace.length} combinations`}` +
    ` of these voltages — ${nullSpace.map(say).join('; ')} — so the setting above is not unique.`
  );
}

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

  // What the second stage did, if it ran. Said separately because it answers a
  // different question from the scan, and the flat directions it finds are a
  // statement about the instrument rather than about the search.
  const ref = result.refinement;
  const polish = !ref
    ? ''
    : (ref.improved
        ? ` Reduced-Hessian refinement tightened the beam from ${mToMm(ref.rmsBefore).toFixed(
            2
          )} to ${mToMm(ref.rmsAfter).toFixed(2)} mm in ${ref.steps} Newton ${
            ref.steps === 1 ? 'step' : 'steps'
          }.`
        : ' Reduced-Hessian refinement found nothing better: the scan had already ' +
          'reached the best setting nearby.') + describeNullSpace(ref.nullSpace);

  setTuneNote(
    (result.cancelled
      ? `Stopped after ${result.evaluations} trials — keeping the best: ${headline}.` +
        (settings ? ` ${settings}` : '')
      : changed.length === 0
        ? `${result.evaluations} trials: ${headline}. Nothing beat the settings already in place.`
        : `${result.evaluations} trials: ${headline}. ${settings}`) + polish,
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
 * Below what field strength a line has no direction worth following.
 *
 * Taken from the MEAN field over the element, not the peak. The peak sits on a
 * sharp electrode rim, where the field genuinely diverges as the grid is
 * refined (§11 of docs/PHYSICS.md) - so a threshold set as a fraction of it is
 * set by a singularity, and scales with the grid step rather than with the
 * physics. Measured: it came out high enough to stop every line within a
 * millimetre of the metal, leaving the bore and the drift blank.
 *
 * A mean is still pulled up by that rim, but only in proportion to how many
 * nodes are near it, which is few.
 */
function fieldFloor(field) {
  let sum = 0;
  let n = 0;
  for (let k = 0; k < field.Ez.length; k++) {
    const m = Math.hypot(field.Ez[k], field.Er[k]);
    if (m > 0) {
      sum += m;
      n++;
    }
  }
  return n ? (sum / n) * 0.02 : 0;
}

/**
 * Field lines: curves everywhere tangent to E.
 *
 * Traced rather than contoured. Starting from a seed, each step moves one
 * fraction of a grid cell along the unit field vector, which is the definition
 * of a field line and needs no marching-squares machinery. Both directions are
 * followed from each seed, so a line runs from the positive metal it leaves to
 * the negative metal it lands on rather than stopping at an arbitrary point.
 *
 * Seeded on a coarse lattice rather than on the electrodes. Electrode surfaces
 * are where the field is strongest, so seeding there crowds every line into
 * the gaps and leaves the interesting part - what escapes into the drift -
 * empty. A lattice covers both, and duplicate lines through the same cell are
 * dropped so the strong regions do not end up solid black.
 *
 * Axisymmetric elements only. A quadrupole's field lies in the transverse
 * plane, not this one, and a deflector's bend plane is drawn edge-on, so in
 * both cases tracing in the r-z slice would draw something that is not there.
 */
function drawFieldLines(e, T) {
  if (!e.grid || !e.field) return;
  // A quadrupole's or multipole's field lies in the plane TRANSVERSE to the
  // beam, so it does not appear in either of these views at all - drawing
  // anything for it here would be drawing a field that is not in this plane.
  // Those are handled in the cross-section instead.
  if (e.typeKey === 'quadrupole' || e.typeKey === 'multipole') return;
  if (e.typeKey === 'bender') return drawBenderFieldLines(e, T);
  const { grid, field } = e;
  const h = grid.step * 0.5;
  const maxSteps = Math.round((grid.zLength / h) * 1.5);

  const floor = fieldFloor(field);
  if (!(floor > 0)) return;

  // One line per cell of this lattice, so density is even rather than
  // following field strength.
  const cell = grid.step * 7;
  const claimed = new Set();
  const key = (z, r) => `${Math.round(z / cell)},${Math.round(r / cell)}`;

  const lines = [];
  for (let z = grid.z0 + cell; z < grid.z0 + grid.zLength; z += cell) {
    for (let r = grid.rMin + cell * 0.5; r < grid.rMax; r += cell) {
      if (claimed.has(key(z, r))) continue;
      const pts = [];
      for (const sense of [1, -1]) {
        let pz = z;
        let pr = r;
        const side = [];
        for (let n = 0; n < maxSteps; n++) {
          const { Ez, Er } = field.fieldAt(pz, pr);
          const mag = Math.hypot(Ez, Er);
          // A field line has no direction where there is no field.
          if (!(mag > floor)) break;
          side.push([pz, pr]);
          claimed.add(key(pz, pr));
          pz += (sense * h * Ez) / mag;
          pr += (sense * h * Er) / mag;
          if (pz < grid.z0 || pz > grid.z0 + grid.zLength) break;
          if (pr < grid.rMin || pr > grid.rMax) break;
          if (field.strikes(pr, 0, pz)) break;
        }
        // The backward half is walked outwards from the seed, so it has to be
        // reversed before being joined to the forward half.
        if (sense === 1) pts.push(...side.reverse());
        else pts.push(...side.slice(1));
      }
      if (pts.length > 3) lines.push(pts);
    }
  }
  if (lines.length === 0) return;

  ctx.save();
  withElementTransform(e, T, () => {});
  ctx.strokeStyle = cssVar('--accent');
  ctx.globalAlpha = 0.45;
  ctx.lineWidth = 1;
  ctx.beginPath();
  for (const pts of lines) {
    for (const sign of [1, -1]) {
      pts.forEach((p, i) => {
        const x = p[0] * T.scale;
        const y = -sign * p[1] * T.scale;
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      });
    }
  }
  ctx.stroke();
  ctx.restore();
}

/**
 * Field lines inside a quadrupole deflector.
 *
 * Worth a separate routine because a deflector is the one element whose field
 * and whose trajectory share a plane, and that plane is the one being drawn.
 * Everything else here is either axisymmetric, so its field lies in the r-z
 * slice, or transverse, so it does not appear in this view at all. A deflector
 * is neither: its field is a quadrupole lying flat in the bend plane, and the
 * lines running from the positive electrodes to the negative ones ARE the
 * thing that turns the beam.
 *
 * Drawn only in the pane that shows that plane. Rolled to bend vertically, the
 * field is in y-z and drawing it on the top view would be drawing a field
 * edge-on and pretending it was in the page.
 */
function drawBenderFieldLines(e, T) {
  const { grid, field } = e;
  const roll = ((e.params.bendPlane ?? 0) * Math.PI) / 180;
  const cosR = Math.cos(roll);
  const sinR = Math.sin(roll);

  // Does this pane show the plane the bend happens in? The pane's transverse
  // axis is x for the top view and y for the side.
  const inThisPane = Math.abs(T.axis === 0 ? cosR : sinR);
  if (inThisPane < 0.7) return;

  const floor = fieldFloor(field);
  if (!(floor > 0)) return;

  const extent = -grid.z0;
  const a = mmToM(
    e.params.apertureRadius + e.params.electrodeThickness + e.params.boxClearance
  );
  const h = grid.step * 0.5;
  const maxSteps = Math.round((4 * extent) / h);
  const cell = grid.step * 8;
  const claimed = new Set();
  const key = (z, x) => `${Math.round(z / cell)},${Math.round(x / cell)}`;

  /** Bend-plane (Z, X) to a point in the element's own frame. */
  const toLocal3 = (Z, X) => [X * cosR, X * sinR, Z + a];

  const lines = [];
  for (let Z = -extent + cell; Z < extent; Z += cell) {
    for (let X = -extent + cell; X < extent; X += cell) {
      if (claimed.has(key(Z, X))) continue;
      if (field.strikes(X, 0, Z)) continue;
      const pts = [];
      for (const sense of [1, -1]) {
        let pz = Z;
        let px = X;
        const side = [];
        for (let n = 0; n < maxSteps; n++) {
          const { Ez, Er } = field.fieldAt(pz, px);
          const mag = Math.hypot(Ez, Er);
          if (!(mag > floor)) break;
          side.push([pz, px]);
          claimed.add(key(pz, px));
          pz += (sense * h * Ez) / mag;
          px += (sense * h * Er) / mag;
          if (Math.abs(pz) > extent || Math.abs(px) > extent) break;
          if (field.strikes(px, 0, pz)) break;
        }
        if (sense === 1) pts.push(...side.reverse());
        else pts.push(...side.slice(1));
      }
      if (pts.length > 3) lines.push(pts);
    }
  }
  if (lines.length === 0) return;

  ctx.save();
  ctx.strokeStyle = cssVar('--accent');
  ctx.globalAlpha = 0.45;
  ctx.lineWidth = 1;
  ctx.beginPath();
  for (const pts of lines) {
    pts.forEach(([Z, X], i) => {
      const [sx, sy] = T.project(toGlobal(e.frame, toLocal3(Z, X)));
      if (i === 0) ctx.moveTo(sx, sy);
      else ctx.lineTo(sx, sy);
    });
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
  ctx.save();
  ctx.strokeStyle = cssVar('--axis');
  ctx.lineWidth = 1;
  ctx.setLineDash([4, 4]);
  // One stroke per route. Drawn as a single polyline, the jump from the end of
  // one branch back to the start of the next appears as a line that is not
  // there - and on a switched column it cuts straight across the diagram.
  for (const pts of beamline.centreLines()) {
    if (pts.length < 2) continue;
    ctx.beginPath();
    pts.forEach((p, i) => {
      const [px, py] = T.project(p);
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    });
    ctx.stroke();
  }
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
    ctx.globalAlpha = 0.1;
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
  /*
    How tall each pane may be.

    True scale is kept by `fitScale`, so this is a budget rather than a
    setting: give the panes more room and the whole column is drawn larger,
    up to the point where the aspect ratio of the beamline itself is the
    limit. The old ceiling of 320 px was well below that for any column
    shorter than about a metre, which is all of them.

    Full screen gets the display; otherwise a little over half the window,
    which leaves the beamline track and the readout visible beneath.
  */
  const budget = isFullscreen()
    ? window.innerHeight - 16
    : Math.min(760, Math.round(window.innerHeight * 0.58));
  const paneHeight = Math.max(
    150,
    Math.min(
      Math.floor(budget / panes),
      Math.round((cssWidth * spanT) / Math.max(1e-9, spanZ))
    )
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
  if (inputs.showLines.checked) {
    for (const e of beamline.elements) drawFieldLines(e, T);
  }
  drawBoundaries(T);
  drawReferencePath(T);
  drawElectrodes(T);
  drawSelection(T);
  drawTrajectories(T);
  drawSource(T);

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
    stats.detected ? `${stats.detected} detected` : null,
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
  describeRepulsion();
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
    // Mathieu numbers, the bender's matched voltage and the step a control
    // nudges by all depend on the ion, so the inspector follows the beam.
    // Refreshed in place rather than rebuilt: a rebuild would pull the focus
    // out of whatever is being typed into.
    if (id === 'mass' || id === 'charge' || id === 'energy') {
      renderTrack();
      refreshInspector();
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
    case 'send': {
      // `index` picks a direction here, not an element; -1 means "switch it
      // off and let the beam through".
      const i = selectedIndex();
      const e = beamline.elements[i];
      if (e?.typeKey !== 'bender') return true;
      const { matched } = bendState(e);
      if (index < 0) {
        setParam(i, 'voltage', 0);
      } else {
        const d = BEND_DIRECTIONS[index];
        // Plane first, then polarity: both are fast adjusts, but the readout
        // in between should never show a bend aimed at the old direction.
        setParam(i, 'bendPlane', d.deg);
        const spec = ELEMENT_TYPES.bender.fields.find((f) => f.key === 'voltage');
        const { step } = fieldRange(spec, e.params, beamSpec());
        setParam(i, 'voltage', d.sign * Math.round(matched / step) * step);
      }
      renderInspector();
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

/**
 * The exit a toolbar click will attach to, if one has been chosen.
 *
 * Clicking an empty socket in the track arms it; the next element placed goes
 * there and it disarms. Without this, starting a second branch would need
 * drag-and-drop onto a diagram where the two lines may overlap.
 */
let pendingPort = null;

const TRACK_HINT =
  'Click a box to select it, or the ions to see the beam. Drag a box onto ' +
  'another to re-hang it, and drag from the toolbar onto a socket to add.';

/* ------------------------------------------------------------------ */
/* dragging in the flow chart                                          */
/* ------------------------------------------------------------------ */

/**
 * Re-hanging an element by dragging its box.
 *
 * This used to be done on the diagram, where it was a poor fit: the diagram
 * draws a column to scale, so two branches can overlap, a long drift can be a
 * hundred times the size of the plate beside it, and what you are aiming at is
 * a picture of hardware rather than a picture of the structure you are
 * editing. The flow chart is the structure, drawn at a legible size with one
 * box per element - so that is where the editing gesture belongs, and the
 * diagram is left to do the one thing it is good at.
 */
let flowDrag = null;

/** The element or socket under a client point, if any. */
function flowTargetAt(clientX, clientY) {
  const hit = document.elementFromPoint(clientX, clientY)?.closest('[data-act]');
  if (!hit) return null;
  const index = Number(hit.dataset.index);
  if (hit.dataset.act === 'port') {
    return { kind: 'port', parent: beamline.elements[index] ?? null, port: hit.dataset.port };
  }
  if (hit.dataset.act === 'select' && beamline.elements[index]) {
    return { kind: 'element', element: beamline.elements[index] };
  }
  return null;
}

flowEl.addEventListener('pointerdown', (e) => {
  const node = e.target.closest('[data-act="select"]');
  if (!node || e.target.closest('[data-act="remove"]')) return;
  const element = beamline.elements[Number(node.dataset.index)];
  if (!element) return;
  flowDrag = { element, startX: e.clientX, startY: e.clientY, moved: false };
});

flowEl.addEventListener('pointermove', (e) => {
  if (!flowDrag) return;
  if (Math.hypot(e.clientX - flowDrag.startX, e.clientY - flowDrag.startY) > 5) {
    flowDrag.moved = true;
    flowEl.classList.add('dragging');
  }
});

flowEl.addEventListener('pointerup', (e) => {
  const drag = flowDrag;
  flowDrag = null;
  flowEl.classList.remove('dragging');
  if (!drag || !drag.moved) return;

  const target = flowTargetAt(e.clientX, e.clientY);
  if (!target) return;

  // Onto a socket: hang it on that exact exit. Onto a box: hang it on that
  // element's first free exit, or push whatever is there downstream.
  const parent = target.kind === 'port' ? target.parent : target.element;
  if (!parent || parent === drag.element) return;
  const port =
    target.kind === 'port'
      ? target.port
      : (exitsOf(parent).find((x) => !beamline.childAt(parent, x.port)) ?? exitsOf(parent)[0]).port;

  if (beamline.reparent(drag.element, parent, port)) {
    selection = { kind: 'element', index: beamline.elements.indexOf(drag.element) };
    afterStructureChange();
  }
});

flowEl.addEventListener('pointercancel', () => {
  flowDrag = null;
  flowEl.classList.remove('dragging');
});

/* Dropping a new element from the toolbar onto a socket or a box. */
flowEl.addEventListener('dragover', (e) => {
  if (!e.dataTransfer.types.includes('text/iontrace-element')) return;
  e.preventDefault();
  e.dataTransfer.dropEffect = 'copy';
});

flowEl.addEventListener('drop', (e) => {
  const type = e.dataTransfer.getData('text/iontrace-element');
  if (!type || !ELEMENT_TYPES[type]) return;
  e.preventDefault();

  const target = flowTargetAt(e.clientX, e.clientY);
  const parent = target ? (target.kind === 'port' ? target.parent : target.element) : null;
  const attach = parent
    ? {
        parent,
        port:
          target.kind === 'port'
            ? target.port
            : (exitsOf(parent).find((x) => !beamline.childAt(parent, x.port)) ??
                exitsOf(parent)[0]).port,
      }
    : pendingPort;
  armPort(null);
  addElement(type, attach);
});

/**
 * Arm an exit for the next element placed, and say so.
 *
 * Clicking a socket used to change a border and nothing else, which is not
 * enough to tell anyone that the click worked, let alone what to do next: the
 * button appeared to be broken. It now names the exit it armed and says where
 * the element comes from, and the toolbar lights up to match.
 */
function armPort(next) {
  pendingPort = next;
  for (const b of document.querySelectorAll('[data-act="port"]')) {
    const mine =
      next &&
      beamline.elements[Number(b.dataset.index)] === next.parent &&
      b.dataset.port === next.port;
    b.classList.toggle('armed', Boolean(mine));
  }
  toolsEl.classList.toggle('awaiting', Boolean(next));

  if (!next) {
    trackHint.textContent = TRACK_HINT;
    trackHint.classList.remove('armed');
    return;
  }
  const exit = exitsOf(next.parent).find((x) => x.port === next.port);
  const where =
    exitsOf(next.parent).length > 1
      ? `the ${next.parent.label}’s ${exit.label.toLowerCase()} exit`
      : `after the ${next.parent.label}`;
  trackHint.textContent = `Pick an element from the toolbar above to put it ${where}, or drag one onto the socket.`;
  trackHint.classList.add('armed');
}

/**
 * A click in either view of the beamline.
 *
 * The list and the flow chart show the same tree and answer to the same
 * gestures, so they share one handler rather than growing two that drift
 * apart. `[data-act]` is on a <button> in one and a <g> in the other, which is
 * the only difference that reaches here.
 */
function onStructureClick(e) {
  const port = e.target.closest('[data-act="port"]');
  if (port) {
    const parent = beamline.elements[Number(port.dataset.index)] ?? null;
    const already = pendingPort?.parent === parent && pendingPort?.port === port.dataset.port;
    armPort(already ? null : { parent, port: port.dataset.port });
    return;
  }
  const hit = e.target.closest('[data-act]');
  if (hit) handleAction(hit.dataset.act, Number(hit.dataset.index));
}

flowEl.addEventListener('click', onStructureClick);

inspectorEl.addEventListener('click', (e) => {
  const button = e.target.closest('button[data-act]');
  if (button) handleAction(button.dataset.act, Number(button.dataset.index));
});

toolsEl.addEventListener('click', (e) => {
  const button = e.target.closest('button[data-add]');
  if (!button) return;
  const attach = pendingPort ?? defaultAttach();
  armPort(null);
  addElement(button.dataset.add, attach);
});

toolsEl.addEventListener('dragstart', (e) => {
  const button = e.target.closest('button[data-add]');
  if (!button) return;
  e.dataTransfer.setData('text/iontrace-element', button.dataset.add);
  e.dataTransfer.effectAllowed = 'copy';
});

/*
  Nothing is dropped onto the diagram any more; it is dropped onto the flow
  chart. The diagram draws a column to scale, which is the right thing for
  seeing where ions go and the wrong thing for aiming at: two branches can
  overlap, and a long drift can be a hundred times the size of the plate beside
  it, so what you hit is not reliably what you meant. Editing happens on the
  picture of the structure; the diagram shows the result.
*/

/* ------------------------------------------------------------------ */
/* saving and loading                                                  */
/* ------------------------------------------------------------------ */

/** Everything outside the beamline object that belongs in a saved file. */
function sceneExtras() {
  return {
    source: {
      ...beamSpec(),
      rays: readNumber(inputs.rays, 9),
      beamRadius: readNumber(inputs.beamRadius, 1),
      divergence: readNumber(inputs.divergence, 0),
    },
    repulsion: inputs.repulsion.value,
    beamCurrent: readNumber(inputs.beamCurrent, 5),
    ionsPerParticle: readNumber(inputs.ionsPerParticle, 6),
  };
}

saveBtn.addEventListener('click', () => {
  const data = serialise(beamline, sceneExtras());
  const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  const stamp = new Date().toISOString().slice(0, 16).replace(/[:T]/g, '-');
  a.download = `iontrace-${stamp}.json`;
  a.click();
  // Revoked on the next turn of the loop, once the click has been handled.
  setTimeout(() => URL.revokeObjectURL(url), 0);
  setTuneNote(
    `Saved ${beamline.elements.length} elements as ${a.download}. The file holds the ` +
      'column and the beam, not the solved fields — loading it re-solves.'
  );
});

loadBtn.addEventListener('click', () => loadFile.click());

loadFile.addEventListener('change', async () => {
  const file = loadFile.files?.[0];
  loadFile.value = ''; // so the same file can be loaded twice running
  if (!file) return;

  let data;
  try {
    data = JSON.parse(await file.text());
  } catch (err) {
    setTuneNote(`Could not read ${file.name}: ${err.message}`, true);
    return;
  }

  statusEl.classList.add('busy');
  let restored;
  try {
    restored = restore(data);
  } catch (err) {
    statusEl.classList.remove('busy');
    setTuneNote(`Could not rebuild that column: ${err.message}`, true);
    return;
  }

  if (!restored.elements.length) {
    statusEl.classList.remove('busy');
    setTuneNote(`Nothing loaded from ${file.name}. ${restored.problems.join(' ')}`, true);
    return;
  }

  // The source first, so the elements are laid out against the right ion.
  const s = restored.source;
  if (s) {
    if (s.mass != null) inputs.mass.value = s.mass;
    if (s.charge != null) inputs.charge.value = s.charge;
    if (s.energy != null) inputs.energy.value = s.energy;
    if (s.rays != null) inputs.rays.value = s.rays;
    if (s.beamRadius != null) inputs.beamRadius.value = s.beamRadius;
    if (s.divergence != null) inputs.divergence.value = s.divergence;
  }
  const phys = restored.physics ?? {};
  if (phys.repulsion) inputs.repulsion.value = phys.repulsion;
  if (phys.beamCurrent != null) inputs.beamCurrent.value = phys.beamCurrent;
  if (phys.ionsPerParticle != null) inputs.ionsPerParticle.value = phys.ionsPerParticle;
  fringeToggle.checked = Boolean(phys.fringe);

  beamline = new Beamline();
  beamline.fringe = Boolean(phys.fringe);
  beamline.adopt(restored.elements);
  statusEl.classList.remove('busy');

  selection = null;
  armPort(null);
  syncOutputs();
  describeFringe();
  afterStructureChange();
  setTuneNote(
    `Loaded ${restored.elements.length} elements from ${file.name}.` +
      (restored.problems.length ? ` ${restored.problems.join(' ')}` : ''),
    restored.problems.length > 0
  );
});

autoAlignBtn.addEventListener('click', () => {
  beamline.autoAlign();
  afterStructureChange();
});

/**
 * Solve neighbouring elements together, or each alone.
 *
 * This is a physics setting, not a display one: it changes the field the ions
 * fly through, so it costs a solve and invalidates the drawn trajectories.
 */
function applyFringe() {
  statusEl.classList.add('busy');
  beamline.setFringe(fringeToggle.checked);
  statusEl.classList.remove('busy');
  describeFringe();
  markStale();
  render();
  drawReadout();
}

/**
 * What the repulsion setting is doing, said plainly.
 *
 * Off is a real choice and a reasonable default - most ion-optics work is
 * single-particle, and repulsion costs a lockstep integration - but it is not
 * a neutral one, and a beam that sails down a metre of drift without widening
 * is otherwise just puzzling.
 */
function describeRepulsion() {
  const model = inputs.repulsion.value;
  repulsionNote.textContent =
    model === 'none'
      ? 'Ions do not see each other, so a beam keeps its width down a field-free ' +
        'drift for ever. Real beams do not: switch this on to let them push apart.'
      : model === 'beam'
        ? 'Each ray is a ring of charge and the force comes from Gauss’s law on the ' +
          'current it encloses, so the spreading is set by the current rather than ' +
          'by how many rays are drawn.'
        : 'Every ion pushes on every other by Coulomb’s law, each standing in for ' +
          'the number of real ions below.';
}

function describeFringe() {
  if (!fringeToggle.checked) {
    fringeNote.textContent =
      'Each element is solved alone behind grounded end faces. Its field stops ' +
      'at its own boundary — which also means those faces shield it, whether ' +
      'or not anything real is there.';
    return;
  }
  const runs = beamline.runs;
  const spanned = runs.reduce((n, r) => n + (r.to - r.from + 1), 0);
  fringeNote.textContent = runs.length
    ? `${spanned} of ${beamline.elements.length} elements solved together in ` +
      `${runs.length} ${runs.length === 1 ? 'group' : 'groups'}. Fields now reach ` +
      'into their neighbours, and a grounded plate stops them. Quadrupoles and ' +
      'deflectors end a group: their own housings genuinely do terminate the field.'
    : 'No axisymmetric elements to solve together — a quadrupole and a deflector ' +
      'are each enclosed already.';
}

fringeToggle.addEventListener('change', applyFringe);

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
 * How long a geometry edit waits before the element is rebuilt.
 *
 * Geometry changes the boundary, so it needs a fresh Laplace solve, and that
 * is not free: measured, a drift rebuilds in 0.2 ms, an aperture plate in 5, a
 * lens in 24, a deflector in 62 and an RF quadrupole in 122. Re-solving on
 * every keystroke would stall on the slower ones and would also solve for the
 * wrong thing on the way - typing "15" means solving for 1 first.
 *
 * So geometry edits are applied after a short idle instead of being held back
 * until the field is left. Long enough that a burst of typing or a held-down
 * stepper arrow collapses into one solve; short enough that the hardware
 * visibly follows the number.
 */
const GEOMETRY_DELAY = 140;

/** Cancel a geometry rebuild that has not happened yet. */
function cancelGeometry() {
  clearTimeout(geometryTimer);
  geometryTimer = null;
}

/**
 * Apply what a typed inspector field says.
 *
 * `live` distinguishes a keystroke from a commit. A voltage only rescales
 * stored solutions, so it follows every keystroke directly. Geometry takes the
 * scheduled path above, so the drawing still follows what is being typed - the
 * plate really does get thicker as the number goes up - without a solve per
 * character.
 */
function applyTypedField(input, live) {
  const key = input.dataset.param;
  if (key) {
    const value = readTyped(input);
    if (value === null) return;
    const geometry = input.dataset.rebuild === '1';
    if (geometry && live) {
      // Pinned to the element selected NOW: the rebuild happens later, and by
      // then the selection may have moved on.
      const index = selectedIndex();
      cancelGeometry();
      geometryTimer = setTimeout(() => {
        geometryTimer = null;
        if (index === selectedIndex()) setParam(index, key, value);
      }, GEOMETRY_DELAY);
      return;
    }
    cancelGeometry();
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
viewFly.addEventListener('click', fly);
fullscreenBtn.addEventListener('click', toggleFullscreen);

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
/* help                                                                */
/* ------------------------------------------------------------------ */

const helpPop = el('helpPop');
const helpPopTitle = el('helpPopTitle');
const helpPopBody = el('helpPopBody');
const helpPopMore = el('helpPopMore');

/**
 * A question mark carrying its own text.
 *
 * Used for the per-parameter explanations, which live beside the parameter in
 * the element registry rather than in the topic list - they are about one
 * number on one element, and hoisting them somewhere central would separate
 * them from the thing they describe.
 */
function helpButton(text, about = '') {
  return `<button class="help-btn" type="button" data-help-text="${escapeHtml(text)}"
    data-help-title="${escapeHtml(about)}" aria-label="About ${escapeHtml(about || 'this')}">?</button>`;
}

/** A question mark pointing at a shared topic. */
function helpTopicButton(key, about = '') {
  return `<button class="help-btn" type="button" data-help="${escapeHtml(key)}"
    aria-label="About ${escapeHtml(about || key)}">?</button>`;
}

let helpAnchor = null;

function hideHelp() {
  helpPop.hidden = true;
  helpAnchor?.setAttribute('aria-expanded', 'false');
  helpAnchor = null;
}

/**
 * Put the popover beside its button, then pull it back on screen.
 *
 * Preferred position is just below and left-aligned, which reads naturally
 * next to a control. The clamp matters more than the preference: these buttons
 * sit at the right-hand edge of a narrow sidebar, so the natural position is
 * very often off screen, and a help popover that opens where it cannot be read
 * is worse than none.
 */
function placeHelp(button) {
  const pad = 8;
  const r = button.getBoundingClientRect();
  helpPop.style.left = '0px';
  helpPop.style.top = '0px';
  helpPop.hidden = false;
  const w = helpPop.offsetWidth;
  const h = helpPop.offsetHeight;

  let left = r.left;
  left = Math.min(left, window.innerWidth - w - pad);
  left = Math.max(pad, left);

  // Below unless that runs off the bottom and there is more room above.
  let top = r.bottom + 6;
  if (top + h > window.innerHeight - pad && r.top - h - 6 > pad) top = r.top - h - 6;
  top = Math.max(pad, Math.min(top, window.innerHeight - h - pad));

  helpPop.style.left = `${Math.round(left)}px`;
  helpPop.style.top = `${Math.round(top)}px`;
}

function showHelp(button) {
  const key = button.dataset.help;
  const t = key ? topic(key) : null;
  const title = t?.title ?? button.dataset.helpTitle ?? 'About this';
  const paras = t ? t.body : [button.dataset.helpText ?? ''];

  helpPopTitle.textContent = title;
  helpPopBody.innerHTML = paras
    .filter(Boolean)
    .map((p) => `<p>${escapeHtml(p)}</p>`)
    .join('');
  if (t?.more) {
    helpPopMore.textContent = `The full argument: ${t.more}`;
    helpPopMore.hidden = false;
  } else {
    helpPopMore.hidden = true;
  }

  helpAnchor?.setAttribute('aria-expanded', 'false');
  helpAnchor = button;
  button.setAttribute('aria-expanded', 'true');
  placeHelp(button);
}

/*
  One listener for every question mark in the document, present or future.

  The inspector rebuilds its HTML whenever the selection changes, so its help
  buttons are not the same elements from one moment to the next. Delegating
  from the document is what lets a button that did not exist when this ran
  still work - and `preventDefault` matters because several of these sit inside
  a <label>, where a click would otherwise be forwarded to the input.
*/
document.addEventListener('click', (ev) => {
  const btn = ev.target.closest?.('.help-btn');
  if (btn) {
    ev.preventDefault();
    ev.stopPropagation();
    if (helpAnchor === btn) hideHelp();
    else showHelp(btn);
    return;
  }
  if (!helpPop.hidden && !ev.target.closest?.('.help-pop')) hideHelp();
});

el('helpPopClose').addEventListener('click', hideHelp);
window.addEventListener('keydown', (ev) => {
  if (ev.key === 'Escape' && !helpPop.hidden) hideHelp();
});
window.addEventListener('resize', () => {
  if (helpAnchor) placeHelp(helpAnchor);
});

/* ------------------------------------------------------------------ */
/* tutorial                                                            */
/* ------------------------------------------------------------------ */

const tabBeamline = el('tabBeamline');
const tabTutorial = el('tabTutorial');
const panelControls = el('panelControls');
const panelTutorial = el('panelTutorial');
const layoutEl = document.querySelector('.layout');

/**
 * Swap the sidebar between the controls and the tutorial.
 *
 * Only the sidebar: the viewport keeps the diagram, the beamline chart and the
 * readout on screen throughout, because every lesson builds a real column and
 * flies it. A tutorial that covered the thing it was describing would be
 * asking the reader to take its word for the result.
 *
 * The layout widens in tutorial mode - 290 px is right for a stack of number
 * inputs and too narrow for prose.
 */
function showTab(which) {
  const tut = which === 'tutorial';
  panelControls.hidden = tut;
  panelTutorial.hidden = !tut;
  layoutEl.classList.toggle('tutoring', tut);
  tabBeamline.setAttribute('aria-selected', String(!tut));
  tabTutorial.setAttribute('aria-selected', String(tut));
  tabBeamline.classList.toggle('on', !tut);
  tabTutorial.classList.toggle('on', tut);
  if (tut) renderTutorial();
}

tabBeamline.addEventListener('click', () => showTab('beamline'));
tabTutorial.addEventListener('click', () => showTab('tutorial'));

const STEPS = allSteps();
let stepAt = 0;

/**
 * What a tutorial step is allowed to do.
 *
 * Everything here is something a person could do with the ordinary controls,
 * which is the point: a step cannot demonstrate behaviour the interface cannot
 * reproduce, so following the tutorial teaches the actual program rather than
 * a private back door into the model.
 */
const tutorialApi = {
  build(specs) {
    statusEl.classList.add('busy');
    const made = specs.map((s) => createElement(s.type, s.params ?? {}));
    for (let i = 1; i < made.length; i++) {
      made[i].from = { parent: made[i - 1], port: exitsOf(made[i - 1])[0].port };
    }
    if (made.length) made[0].from = { parent: null, port: 'out' };
    beamline = new Beamline();
    beamline.fringe = fringeToggle.checked;
    beamline.adopt(made);
    statusEl.classList.remove('busy');
    selection = null;
    armPort(null);
    afterStructureChange();
  },

  beam(spec) {
    for (const [k, v] of Object.entries(spec)) {
      if (inputs[k]) inputs[k].value = String(v);
    }
    syncOutputs();
    markStale();
    render();
  },

  physics(spec) {
    if (spec.repulsion !== undefined) inputs.repulsion.value = spec.repulsion;
    if (spec.beamCurrent !== undefined) inputs.beamCurrent.value = String(spec.beamCurrent);
    if (spec.fringe !== undefined) {
      fringeToggle.checked = Boolean(spec.fringe);
      statusEl.classList.add('busy');
      beamline.setFringe(fringeToggle.checked);
      statusEl.classList.remove('busy');
      describeFringe();
    }
    syncOutputs();
    markStale();
    render();
  },

  view(spec) {
    for (const [k, v] of Object.entries(spec)) {
      const key = `show${k[0].toUpperCase()}${k.slice(1)}`;
      if (inputs[key]) inputs[key].checked = Boolean(v);
    }
    render();
  },

  select(target) {
    if (target === null) select(null);
    else if (target === 'source') select({ kind: 'source' });
    else select({ kind: 'element', index: target });
  },

  fly() {
    fly();
  },

  async optimise() {
    await runTuner(tunableKnobs(beamline, beamSpec()), optimizeBtn, 'the whole column');
  },
};

function renderTutorial() {
  const here = STEPS[stepAt];
  if (!here) return;
  const { lesson, step, lessonIndex, stepIndex } = here;

  el('tutLessons').innerHTML = LESSONS.map((l, i) => {
    const first = STEPS.findIndex((s) => s.lesson === l);
    return `<li>
      <button type="button" class="tut-lesson${i === lessonIndex ? ' on' : ''}"
              data-step="${first}">
        <span class="tut-lesson-title">${escapeHtml(l.title)}</span>
        <span class="tut-lesson-sum">${escapeHtml(l.summary)}</span>
      </button>
    </li>`;
  }).join('');

  el('tutCrumb').textContent = `${lesson.title} · step ${stepIndex + 1} of ${lesson.steps.length}`;
  el('tutTitle').textContent = step.title;
  el('tutBody').innerHTML = step.body.map((p) => `<p>${escapeHtml(p)}</p>`).join('');

  const doBtn = el('tutDo');
  doBtn.hidden = !step.action;
  if (step.action) doBtn.textContent = step.action.label;

  const whyBtn = el('tutWhy');
  whyBtn.hidden = !step.topic;
  if (step.topic) whyBtn.dataset.help = step.topic;

  el('tutProgress').textContent = `${stepAt + 1} / ${STEPS.length}`;
  el('tutPrev').disabled = stepAt === 0;
  el('tutNext').disabled = stepAt === STEPS.length - 1;
}

el('tutLessons').addEventListener('click', (ev) => {
  const b = ev.target.closest('.tut-lesson');
  if (!b) return;
  stepAt = Number(b.dataset.step);
  renderTutorial();
});

el('tutPrev').addEventListener('click', () => {
  stepAt = Math.max(0, stepAt - 1);
  renderTutorial();
});
el('tutNext').addEventListener('click', () => {
  stepAt = Math.min(STEPS.length - 1, stepAt + 1);
  renderTutorial();
});

el('tutDo').addEventListener('click', async () => {
  const step = STEPS[stepAt]?.step;
  if (!step?.action) return;
  const btn = el('tutDo');
  btn.disabled = true;
  try {
    await step.action.run(tutorialApi);
  } catch (err) {
    setTuneNote(`That step could not run: ${err.message}`, true);
  } finally {
    btn.disabled = false;
    renderTutorial();
  }
});

/* ------------------------------------------------------------------ */
/* start                                                               */
/* ------------------------------------------------------------------ */

statusEl.classList.add('busy');
beamline = defaultBeamline();
statusEl.classList.remove('busy');
selection = null;

syncOutputs();
describeFringe();
renderTools();
armPort(null);
showTab('beamline');
renderTutorial();
renderTrack();
renderInspector();
render();
drawReadout();
fly();
