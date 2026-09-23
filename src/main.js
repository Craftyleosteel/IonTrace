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
import { ELEMENT_TYPES, createElement, needsRebuild } from './elements/index.js';
import { MATHIEU_Q_LIMIT } from './elements/quadrupole.js';
import { makeIon, parallelBeam, focalCrossing } from './ion.js';
import { createFlight, kineticEnergy } from './integrator.js';
import { joulesToEV, mToMm, mmToM } from './constants.js';
import { NO_ELECTRODE } from './grid.js';

const el = (id) => document.getElementById(id);

const canvas = el('scene');
const ctx = canvas.getContext('2d');
const cross = el('cross');
const crossCtx = cross.getContext('2d');
const statusEl = el('status');
const readoutEl = el('readout');
const trackEl = el('track');
const addersEl = el('adders');
const inspectorEl = el('inspector');
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
  zoom: el('zoom'),
  showField: el('showField'),
  showContours: el('showContours'),
};

/** Controls that change only how the scene is drawn, never the physics. */
const DISPLAY_INPUTS = ['showField', 'showContours', 'zoom'];

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
let selected = 0;
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
 * The einzel is set where it matches the beam into the quadrupole's
 * acceptance rather than merely where it focuses: at -400 V it over-focuses
 * and one ion in nine survives, at 0 V the beam is too wide entering the rods
 * and three survive, and at -150 V seven do. Matching one element to the next
 * is most of what building a column is.
 */
function defaultBeamline() {
  return new Beamline([
    createElement('drift', { length: 12, bore: 5 }),
    createElement('einzel', {
      gridStep: 0.5,
      voltage: -150,
      boreRadius: 5,
      housingRadius: 14,
      entryDrift: 8,
      exitDrift: 8,
    }),
    createElement('drift', { length: 14, bore: 5 }),
    // 150 mm at 2 MHz gives a 50 eV ion about thirty RF cycles in the rods.
    // That number matters more than it looks: stability is an asymptotic
    // property of the Mathieu equation, and an ion that crosses in a handful
    // of cycles can be thrown out whatever its (a, q) says. A shorter or
    // faster-crossing quadrupole loses ions for that reason alone.
    createElement('quadrupole', {
      gridStep: 0.25,
      length: 150,
      rfAmplitude: 250,
      frequency: 2,
    }),
    createElement('drift', { length: 25, bore: 5 }),
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
  const index = beamline.elements.length
    ? Math.min(selected + 1, beamline.elements.length)
    : 0;
  statusEl.classList.add('busy');
  const created = createElement(type);
  beamline.add(created, index);
  statusEl.classList.remove('busy');
  selected = index;
  afterStructureChange();
}

function removeElement(index) {
  if (beamline.elements.length <= 1) return;
  beamline.remove(index);
  selected = Math.max(0, Math.min(selected, beamline.elements.length - 1));
  afterStructureChange();
}

function moveElement(index, delta) {
  if (!beamline.move(index, delta)) return;
  selected = index + delta;
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
    default:
      return '';
  }
}

function renderTrack() {
  trackEl.innerHTML = beamline.elements
    .map((e, i) => {
      const last = beamline.elements.length - 1;
      return `
        <li class="chip ${i === selected ? 'sel' : ''} chip-${e.typeKey}"
            data-index="${i}">
          <button class="chip-body" data-act="select" data-index="${i}">
            <span class="chip-name">${escapeHtml(e.label)}</span>
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
}

function renderAdders() {
  addersEl.innerHTML = Object.entries(ELEMENT_TYPES)
    .map(
      ([type, spec]) =>
        `<button class="add" data-add="${type}" title="${escapeHtml(spec.blurb)}">
           + ${escapeHtml(spec.label)}
         </button>`
    )
    .join('');
}

function renderInspector() {
  const e = beamline.elements[selected];
  if (!e) {
    inspectorEl.innerHTML = '<h2>Element</h2><p class="hint">Nothing selected.</p>';
    return;
  }
  const spec = ELEMENT_TYPES[e.typeKey];

  const rows = spec.fields
    .map((f) => {
      const value = e.params[f.key];
      const instant = f.rebuild ? '' : '<span class="instant" title="No re-solve needed">fast</span>';
      return `
        <label class="field">
          <span class="field-label">
            ${escapeHtml(f.label)}${instant}
            <span class="unit">${escapeHtml(f.unit ?? '')}</span>
          </span>
          <input type="range" data-param="${f.key}"
                 min="${f.min}" max="${f.max}" step="${f.step}" value="${value}" />
          <output data-out="${f.key}">${value}</output>
          ${f.help ? `<span class="field-help">${escapeHtml(f.help)}</span>` : ''}
        </label>`;
    })
    .join('');

  inspectorEl.innerHTML = `
    <h2>${escapeHtml(e.label)}</h2>
    <p class="hint">${escapeHtml(spec.blurb)}</p>
    ${rows}
    ${quadrupoleReadout(e)}`;
}

/**
 * Mathieu parameters for the selected quadrupole and the current ion.
 *
 * These, not the voltages, are what decide whether an ion is transmitted, so
 * they belong next to the controls that set them.
 */
function quadrupoleReadout(e) {
  if (e.typeKey !== 'quadrupole') return '';
  const mass = readNumber(inputs.mass, 100);
  const charge = readNumber(inputs.charge, 1);
  const { a, q } = e.mathieu(mass, Math.abs(charge) || 1);
  const stable = Math.abs(q) < MATHIEU_Q_LIMIT && Math.abs(a) < 0.237;
  return `
    <div class="mathieu ${stable ? 'ok' : 'bad'}">
      <span>a = ${a.toFixed(4)}</span>
      <span>q = ${q.toFixed(4)}</span>
      <span class="verdict">${
        stable
          ? 'inside the first stability region'
          : `outside it — q limit is ${MATHIEU_Q_LIMIT}`
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

function startFlight() {
  const spec = {
    mass: readNumber(inputs.mass, 100),
    charge: readNumber(inputs.charge, 1),
    energy: readNumber(inputs.energy, 50),
    z: 0.2,
  };
  const count = Math.round(readNumber(inputs.rays, 9));
  const maxOffset = readNumber(inputs.beamRadius, 1.5);
  const divergence = readNumber(inputs.divergence, 0);

  let ions;
  try {
    ions = parallelBeam({ ...spec, count, maxOffset });
    // A real source is not perfectly collimated. Divergence is applied as a
    // linear fan so the outermost ion gets the full angle, which is the usual
    // way a beam's emittance is sketched.
    if (divergence !== 0 && count > 1) {
      ions = ions.map((ion, i) => {
        const frac = (2 * i) / (count - 1) - 1;
        return makeIon({
          ...spec,
          x: mToMm(ion.x),
          angle: frac * divergence,
        });
      });
    }
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
/* rendering                                                           */
/* ------------------------------------------------------------------ */

/**
 * Screen transform.
 *
 * The axial and transverse scales are equal by default, because stretching the
 * transverse axis makes trajectory angles unreadable. A long column is a very
 * thin strip at true scale, so the zoom control exists - and whenever it is
 * not 1, the view is labelled, because an unlabelled exaggerated plot is a
 * misleading one.
 */
function makeTransform(width, height, zoom) {
  const half = height / 2;
  const scale = width / beamline.length;
  return {
    sx: (z) => z * scale,
    sy: (x) => half - x * scale * zoom,
    scale,
    zoom,
  };
}

/** Half-height the view needs, in metres, at the given zoom. */
function viewHalfHeight(zoom) {
  return beamline.radiusLimit / zoom;
}

function drawElementField(e, T, width, height, zoom) {
  // Only elements with an axisymmetric (z, r) map have something meaningful
  // to paint in this view. A quadrupole's solve lives in the transverse
  // plane; its slice is drawn in the cross-section inset instead.
  if (!e.grid || !e.field || e.typeKey === 'quadrupole') return;

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

  const x0 = T.sx(e.zStart);
  const x1 = T.sx(e.zEnd);
  const yTop = T.sy(grid.rMax);
  const yBot = T.sy(-grid.rMax);
  ctx.imageSmoothingEnabled = true;
  ctx.imageSmoothingQuality = 'high';
  ctx.drawImage(off, x0, yTop, x1 - x0, yBot - yTop);
}

function drawContours(e, T) {
  if (!e.grid || !e.field || e.typeKey === 'quadrupole') return;
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
  ctx.strokeStyle = cssVar('--gridline');
  ctx.globalAlpha = 0.85;
  ctx.lineWidth = 1;
  ctx.beginPath();

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
            ctx.moveTo(T.sx(e.zStart + grid.zAt(a[0])), T.sy(sign * grid.rAt(a[1])));
            ctx.lineTo(T.sx(e.zStart + grid.zAt(b[0])), T.sy(sign * grid.rAt(b[1])));
          }
        }
      }
    }
  }
  ctx.stroke();
  ctx.restore();
}

function drawElectrodes(T) {
  ctx.save();
  for (const r of beamline.outline()) {
    ctx.fillStyle = r.ghost
      ? cssVar('--electrode-ghost')
      : r.wall
        ? cssVar('--axis')
        : cssVar('--electrode');
    ctx.globalAlpha = r.ghost ? 0.35 : r.wall ? 0.5 : 1;
    for (const sign of [1, -1]) {
      const yA = T.sy(sign * r.r0);
      const yB = T.sy(sign * r.r1);
      ctx.fillRect(
        T.sx(r.z0),
        Math.min(yA, yB),
        Math.max(1, T.sx(r.z1) - T.sx(r.z0)),
        Math.abs(yB - yA)
      );
    }
  }
  ctx.restore();
}

function drawBoundaries(T, height) {
  // Where one element ends and the next begins. Worth showing, because each
  // was solved as a separate problem and the joins are where that
  // approximation lives.
  ctx.save();
  ctx.strokeStyle = cssVar('--gridline');
  ctx.lineWidth = 1;
  ctx.setLineDash([2, 3]);
  ctx.beginPath();
  for (const e of beamline.elements) {
    const x = T.sx(e.zStart);
    ctx.moveTo(x, 0);
    ctx.lineTo(x, height);
  }
  ctx.stroke();
  ctx.restore();
}

function drawAxis(T, width) {
  ctx.save();
  ctx.strokeStyle = cssVar('--axis');
  ctx.lineWidth = 1;
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  ctx.moveTo(0, T.sy(0));
  ctx.lineTo(width, T.sy(0));
  ctx.stroke();
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
        if (n === 0) ctx.moveTo(T.sx(p.z), T.sy(p.x));
        else ctx.lineTo(T.sx(p.z), T.sy(p.x));
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
      ctx.beginPath();
      ctx.arc(T.sx(p.z), T.sy(p.x), 3.5, 0, Math.PI * 2);
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
    ctx.beginPath();
    ctx.arc(T.sx(p.z), T.sy(p.x), 3, 0, Math.PI * 2);
    ctx.fill();
  }
  ctx.restore();
}

function drawScale(T, width, height) {
  ctx.save();
  ctx.fillStyle = cssVar('--text-muted');
  ctx.font = '11px ui-monospace, monospace';
  ctx.textBaseline = 'bottom';
  const totalMm = mToMm(beamline.length);
  const stepMm = totalMm > 400 ? 100 : totalMm > 200 ? 50 : totalMm > 80 ? 20 : 10;
  for (let zmm = 0; zmm <= totalMm + 1e-9; zmm += stepMm) {
    const x = T.sx(mmToM(zmm));
    ctx.fillRect(x, height - 10, 1, 5);
    ctx.textAlign = zmm === 0 ? 'left' : 'center';
    ctx.fillText(`${zmm.toFixed(0)}`, x, height - 12);
  }
  ctx.textAlign = 'right';
  ctx.fillText('z / mm', width - 6, height - 12);
  ctx.restore();
}

/**
 * Beam cross-section, looking down the axis.
 *
 * The main view is the x-z plane, which cannot show that the ions now move in
 * three dimensions. In a quadrupole they emphatically do: the field converges
 * in one transverse plane while diverging in the other, so a beam that looks
 * well behaved from the side can be being pulled into a line seen end-on.
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

  // Framed on the widest bore in the line, so the scale does not jump about.
  const limit = Math.max(...beamline.elements.map((e) => e.bore)) * 1.15;
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

  // The aperture the ions are currently inside.
  const here = trajectories.find((t) => t.active) ?? trajectories[0];
  const zNow = here ? (here.active ? here.state.z : here.points[here.points.length - 1].z) : 0;
  const elementHere = beamline.elementAt(zNow);
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
  for (const traj of trajectories) {
    const p = traj.active ? traj.state : traj.points[traj.points.length - 1];
    if (!p) continue;
    crossCtx.globalAlpha = traj.active ? 1 : 0.35;
    crossCtx.beginPath();
    crossCtx.arc(c + p.x * k, c - (p.y ?? 0) * k, 3, 0, Math.PI * 2);
    crossCtx.fill();
    crossCtx.stroke();
  }
  crossCtx.restore();

  crossCtx.fillStyle = cssVar('--text-muted');
  crossCtx.font = '10px ui-monospace, monospace';
  crossCtx.textAlign = 'center';
  crossCtx.fillText(
    elementHere ? escapeHtml(elementHere.label) : 'beam cross-section',
    c,
    size - 4
  );
}

function render() {
  if (!beamline || beamline.elements.length === 0) return;

  const zoom = readNumber(inputs.zoom, 1);
  const cssWidth = canvas.parentElement.clientWidth;
  const halfH = viewHalfHeight(zoom);
  const aspect = beamline.length / (2 * beamline.radiusLimit) * zoom;
  const cssHeight = Math.max(170, Math.min(480, Math.round(cssWidth / aspect)));

  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.round(cssWidth * dpr);
  canvas.height = Math.round(cssHeight * dpr);
  canvas.style.height = `${cssHeight}px`;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

  ctx.clearRect(0, 0, cssWidth, cssHeight);
  ctx.fillStyle = cssVar('--surface-1');
  ctx.fillRect(0, 0, cssWidth, cssHeight);

  // The transverse scale that makes the visible half-height match the view.
  const T = {
    sx: (z) => (z / beamline.length) * cssWidth,
    sy: (x) => cssHeight / 2 - (x / halfH) * (cssHeight / 2),
  };

  if (inputs.showField.checked) {
    for (const e of beamline.elements) drawElementField(e, T, cssWidth, cssHeight, zoom);
  }
  if (inputs.showContours.checked) {
    for (const e of beamline.elements) drawContours(e, T);
  }
  drawBoundaries(T, cssHeight);
  drawAxis(T, cssWidth);
  drawElectrodes(T);
  drawTrajectories(T);
  drawScale(T, cssWidth, cssHeight);
  drawCrossSection();

  // An exaggerated transverse axis is legitimate but must never be silent.
  const exaggerated = Math.abs(zoom - 1) > 1e-9;
  scaleNote.hidden = !exaggerated;
  if (exaggerated) {
    scaleNote.textContent = `Transverse scale exaggerated ${zoom.toFixed(1)}× — angles are not true`;
  }
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

const OUTPUTS = {
  rays: (v) => v,
  beamRadius: (v) => parseFloat(v).toFixed(1),
  divergence: (v) => parseFloat(v).toFixed(1),
  beamCurrent: (v) => parseFloat(v).toFixed(1),
  ionsPerParticle: (v) => {
    const n = 10 ** parseFloat(v);
    return n < 10 ? n.toFixed(1) : n.toExponential(1).replace('e+', 'e');
  },
  cfl: (v) => parseFloat(v).toFixed(2),
  zoom: (v) => parseFloat(v).toFixed(1),
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
    // Mathieu numbers depend on the ion, so the inspector follows the beam.
    if (id === 'mass' || id === 'charge') renderInspector();
    render();
    if (!displayOnly) drawReadout();
  });
}

// Track: select, reorder, remove.
trackEl.addEventListener('click', (e) => {
  const button = e.target.closest('button[data-act]');
  if (!button) return;
  const index = Number(button.dataset.index);
  switch (button.dataset.act) {
    case 'select':
      selected = index;
      renderTrack();
      renderInspector();
      break;
    case 'left':
      moveElement(index, -1);
      break;
    case 'right':
      moveElement(index, 1);
      break;
    case 'remove':
      removeElement(index);
      break;
  }
});

addersEl.addEventListener('click', (e) => {
  const button = e.target.closest('button[data-add]');
  if (button) addElement(button.dataset.add);
});

// Inspector sliders edit the selected element.
inspectorEl.addEventListener('input', (e) => {
  const input = e.target.closest('input[data-param]');
  if (!input) return;
  const key = input.dataset.param;
  const value = parseFloat(input.value);
  const out = inspectorEl.querySelector(`[data-out="${key}"]`);
  if (out) out.textContent = value;
  setParam(selected, key, value);
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
selected = 1;

syncOutputs();
renderAdders();
renderTrack();
renderInspector();
render();
drawReadout();
fly();
