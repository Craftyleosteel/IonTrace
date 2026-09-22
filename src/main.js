/**
 * IonTrace front end.
 *
 * Wiring only: every physical quantity on screen is computed by the modules in
 * this directory, and nothing here adjusts a number to make a picture look
 * better. The two distinct cost paths are kept visible to the user because
 * they are the point of the architecture:
 *
 *   changing a voltage   -> fast adjust, a weighted sum of stored solutions
 *   changing geometry    -> rebuild the potential array and re-solve Laplace
 *
 * Drawing preserves the true aspect ratio of the lens. Ion optics figures
 * often stretch the radial axis to make deflections legible, which makes
 * trajectory angles unreadable; an honest picture is worth the thin strip.
 */

import { buildEinzelLens } from './geometries/einzel.js';
import { makeIon, parallelBeam, focalCrossing } from './ion.js';
import { flyBeam, flyIon, kineticEnergy } from './integrator.js';
import { joulesToEV, mToMm, mmToM } from './constants.js';
import { NO_ELECTRODE } from './grid.js';

/* ------------------------------------------------------------------ */
/* element lookup                                                      */
/* ------------------------------------------------------------------ */

const el = (id) => document.getElementById(id);

const canvas = el('scene');
const ctx = canvas.getContext('2d');
const statusEl = el('status');
const readoutEl = el('readout');

const inputs = {
  vCentre: el('vCentre'),
  vOuter: el('vOuter'),
  mass: el('mass'),
  charge: el('charge'),
  energy: el('energy'),
  rays: el('rays'),
  beamRadius: el('beamRadius'),
  boreRadius: el('boreRadius'),
  centreLength: el('centreLength'),
  gridStep: el('gridStep'),
  method: el('method'),
  cfl: el('cfl'),
  beamCurrent: el('beamCurrent'),
  showField: el('showField'),
  showContours: el('showContours'),
};

const flyButton = el('fly');

/** Controls whose change invalidates the solved field. */
const GEOMETRY_INPUTS = ['boreRadius', 'centreLength', 'gridStep'];

/**
 * Controls that change only how the scene is drawn.
 *
 * Everything else changes the physics, which makes any trajectory already on
 * screen stale: it was flown through a different field, or by a different
 * beam. Those are marked rather than silently redrawn, because a trajectory
 * that does not correspond to the settings beside it is simply wrong.
 */
const DISPLAY_INPUTS = ['showField', 'showContours'];

/* ------------------------------------------------------------------ */
/* colour                                                              */
/* ------------------------------------------------------------------ */

/** Read a CSS custom property and parse it as [r, g, b]. */
function cssRGB(name) {
  const raw = getComputedStyle(document.documentElement)
    .getPropertyValue(name)
    .trim();
  const hex = raw.replace('#', '');
  const full =
    hex.length === 3
      ? hex.split('').map((c) => c + c).join('')
      : hex;
  return [
    parseInt(full.slice(0, 2), 16),
    parseInt(full.slice(2, 4), 16),
    parseInt(full.slice(4, 6), 16),
  ];
}

function cssVar(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
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
/* simulation state                                                    */
/* ------------------------------------------------------------------ */

let model = null; // { grid, field, geometry }
let trajectories = [];
let stats = {};

/**
 * Animation state.
 *
 * `head` is a sample index, and because flyBeam advances every ion on a
 * shared time step, the same index is the same INSTANT for every ray. The
 * moving dots are therefore a real snapshot of where the bunch is, not a
 * drawing convenience.
 */
let animation = { head: Infinity, total: 0, running: false, frame: 0 };

/** True when the drawn trajectories no longer match the current settings. */
let stale = false;

function readNumber(input, fallback) {
  const v = parseFloat(input.value);
  return Number.isFinite(v) ? v : fallback;
}

/** Rebuild the potential array and solve every electrode's basis solution. */
function rebuild() {
  const t0 = performance.now();
  model = buildEinzelLens({
    boreRadius: readNumber(inputs.boreRadius, 6),
    centreLength: readNumber(inputs.centreLength, 20),
    gridStep: readNumber(inputs.gridStep, 0.5),
  });
  stats.solveMs = performance.now() - t0;
  stats.warnings = model.warnings;
}

/** Apply the current voltages. Cheap: no relaxation happens here. */
function applyVoltages() {
  const t0 = performance.now();
  const outer = readNumber(inputs.vOuter, 0);
  model.field.setVoltages({
    housing: 0,
    entrance: outer,
    centre: readNumber(inputs.vCentre, -2000),
    exit: outer,
  });
  stats.adjustMs = performance.now() - t0;

  // The domain end faces are part of the grounded housing, so they act as
  // solid 0 V plates across the aperture. That is harmless while the outer
  // cylinders are also at 0 V, but biasing them puts a real potential
  // difference across the entry and exit drifts and invents a field of
  // several kV/m where the actual instrument has none. The ion, meanwhile,
  // flies through those faces as if they were open. Flag the inconsistency
  // rather than pretending the result means something.
  stats.voltageWarning =
    outer !== 0
      ? `Entrance/exit biased to ${outer} V against grounded end faces: the ` +
        'drift regions carry a spurious accelerating field. See PHYSICS.md §2.1.'
      : null;
}

/** Fly the beam through the current field. */
function runFlight() {
  const { field, geometry } = model;
  const spec = {
    mass: readNumber(inputs.mass, 100),
    charge: readNumber(inputs.charge, 1),
    energy: readNumber(inputs.energy, 1000),
    z: 0.5,
  };
  const count = Math.round(readNumber(inputs.rays, 9));
  const maxOffset = readNumber(inputs.beamRadius, 4);

  let rays;
  try {
    rays = parallelBeam({ ...spec, count, maxOffset });
  } catch (err) {
    trajectories = [];
    stats.error = err.message;
    return;
  }
  stats.error = null;

  // Slider is in microamps; the physics is in amperes.
  const beamCurrent = readNumber(inputs.beamCurrent, 0) * 1e-6;
  stats.beamCurrent = beamCurrent;

  const opts = {
    method: inputs.method.value,
    cfl: readNumber(inputs.cfl, 0.05),
    beamCurrent,
    // Integrate at full resolution but keep every fourth point for drawing.
    // The energy diagnostic is evaluated on every step regardless, so this
    // costs nothing physically; it only avoids stroking tens of thousands of
    // line segments that land on the same pixels.
    recordEvery: 4,
  };

  const t0 = performance.now();
  // The whole beam flies together. With space charge on it has to: the force
  // on each ion depends on where the others are at that instant, so they
  // cannot be taken one at a time.
  const { tracks } = flyBeam(field, rays, opts);
  trajectories = tracks.map((t, i) => ({
    ...t,
    start: rays[i],
    focus: focalCrossing(t.points),
  }));
  stats.flyMs = performance.now() - t0;

  const lensCentre = (geometry.bounds.z3 + geometry.bounds.z4) / 2; // mm
  const fromCentre = (zMetres) => mToMm(zMetres) - lensCentre;

  // Only a forward exit is transmission. A reflected ion leaves through the
  // entrance face, which is a real and interesting result, but it is not the
  // beam getting through.
  const transmitted = trajectories.filter((t) => t.stop === 'exited');
  stats.transmitted = transmitted.length;
  stats.reflected = trajectories.filter((t) => t.stop === 'reflected').length;
  stats.struck = trajectories.filter((t) => t.stop === 'electrode').length;
  stats.total = trajectories.length;
  stats.drift = Math.max(0, ...trajectories.map((t) => t.energyDrift));

  // The paraxial focus comes from a dedicated probe ray close to the axis,
  // not from an average over the drawn rays. Averaging mixes the marginal and
  // paraxial foci, so the reported "focal length" would shift with the beam
  // radius and the ray count - neither of which is a property of the lens.
  //
  // The probe flies WITHOUT space charge, deliberately. This number is a
  // property of the optics, so it should not move when the beam current is
  // turned up; what the loaded beam actually does is reported separately
  // below. A probe ray flown inside the beam would also be unrepresentative,
  // since at 1 % of the bore it encloses almost no current and would feel
  // almost no self-field.
  const probe = flyIon(field, makeIon({ ...spec, x: 0.01 * geometry.boreRadius }), {
    ...opts,
    beamCurrent: 0,
  });
  const probeFocus = focalCrossing(probe.points);
  stats.paraxial =
    probe.stop === 'exited' && probeFocus
      ? { mm: fromCentre(probeFocus.z), extrapolated: probeFocus.extrapolated }
      : null;

  // What the loaded beam actually does. Space charge goes as 1/r, so a
  // sufficiently dense beam can never be brought to a point: it reaches a
  // minimum radius and expands again. Reporting a focal length in that regime
  // would be reporting something that does not exist.
  const outermost = transmitted
    .slice()
    .sort((a, b) => Math.abs(b.start.x) - Math.abs(a.start.x))[0];
  stats.waist = null;
  if (outermost) {
    const lensEnd = mmToM(geometry.bounds.z6);
    let waist = Infinity;
    let waistZ = 0;
    for (const p of outermost.points) {
      if (p.z < lensEnd) continue;
      if (Math.abs(p.x) < waist) {
        waist = Math.abs(p.x);
        waistZ = p.z;
      }
    }
    if (Number.isFinite(waist)) {
      stats.waist = {
        radiusMm: mToMm(waist),
        atMm: fromCentre(waistZ),
        crossed: outermost.focus !== null && !outermost.focus.extrapolated,
      };
    }
  }

  // Spherical aberration: how far short of the paraxial focus the outermost
  // transmitted ray crosses. Negative means under-corrected, which is what
  // every round electrostatic lens does.
  const marginal = transmitted
    .filter((t) => t.focus)
    .sort((a, b) => Math.abs(b.start.x) - Math.abs(a.start.x))[0];
  stats.aberration =
    marginal && stats.paraxial
      ? fromCentre(marginal.focus.z) - stats.paraxial.mm
      : null;

  // The einzel lens's defining property is that it does no net work, so the
  // figure that matters is the WORST departure across the beam, not one
  // arbitrary ray's.
  if (transmitted.length > 0) {
    let worst = 0;
    let atIn = 0;
    for (const t of transmitted) {
      const kIn = joulesToEV(kineticEnergy(t.points[0]));
      const kOut = joulesToEV(kineticEnergy(t.points[t.points.length - 1]));
      if (Math.abs(kOut - kIn) > Math.abs(worst)) {
        worst = kOut - kIn;
        atIn = kIn;
      }
    }
    stats.keIn = atIn;
    stats.worstWork = worst;
  } else {
    stats.keIn = null;
    stats.worstWork = null;
  }
}

/* ------------------------------------------------------------------ */
/* rendering                                                           */
/* ------------------------------------------------------------------ */

/** Screen transform for the current model and canvas size. */
function makeTransform(width, height) {
  const { grid } = model;
  const zLen = grid.zLength;
  const rMax = grid.rLength;
  return {
    sx: (z) => ((z - grid.z0) / zLen) * width,
    sy: (r) => height / 2 - (r / rMax) * (height / 2),
    zLen,
    rMax,
  };
}

/**
 * Paint the potential as a diverging field, mirrored about the axis.
 *
 * Built at grid resolution into an ImageData and then scaled up, which lets
 * the browser interpolate smoothly instead of showing one hard-edged rectangle
 * per node.
 */
function drawPotential(width, height) {
  const { grid, field } = model;
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

  const rows = 2 * nr - 1; // mirrored: -rMax .. +rMax
  const off = document.createElement('canvas');
  off.width = nz;
  off.height = rows;
  const offCtx = off.getContext('2d');
  const img = offCtx.createImageData(nz, rows);

  for (let row = 0; row < rows; row++) {
    // Row 0 is r = +rMax at the top, so j counts down to the axis and back up.
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
  ctx.imageSmoothingEnabled = true;
  ctx.imageSmoothingQuality = 'high';
  ctx.drawImage(off, 0, 0, width, height);
}

/**
 * Equipotential contours by marching squares on the solved grid.
 *
 * Contours are the most informative single overlay in ion optics: an ion is
 * deflected perpendicular to them, so where they bow across the bore is
 * exactly where focusing happens.
 */
function drawContours(width, height, T) {
  const { grid, field } = model;
  const { nz, nr } = grid;
  const phi = field.phi;

  let maxAbs = 0;
  for (let k = 0; k < phi.length; k++) maxAbs = Math.max(maxAbs, Math.abs(phi[k]));
  if (maxAbs === 0) return;

  const levels = [];
  const N = 11;
  for (let n = 1; n <= N; n++) {
    const frac = n / (N + 1);
    levels.push(maxAbs * frac, -maxAbs * frac);
  }

  ctx.save();
  ctx.strokeStyle = cssVar('--gridline');
  ctx.globalAlpha = 0.85;
  ctx.lineWidth = 1;
  ctx.beginPath();

  const gz = (i) => T.sx(grid.zAt(i));

  for (const level of levels) {
    for (let j = 0; j < nr - 1; j++) {
      for (let i = 0; i < nz - 1; i++) {
        const v00 = phi[j * nz + i];
        const v10 = phi[j * nz + i + 1];
        const v01 = phi[(j + 1) * nz + i];
        const v11 = phi[(j + 1) * nz + i + 1];

        // Edge crossings, in grid coordinates.
        const pts = [];
        const cross = (a, b, ia, ja, ib, jb) => {
          if ((a - level) * (b - level) >= 0) return;
          const t = (level - a) / (b - a);
          pts.push([ia + (ib - ia) * t, ja + (jb - ja) * t]);
        };
        cross(v00, v10, i, j, i + 1, j); // bottom
        cross(v10, v11, i + 1, j, i + 1, j + 1); // right
        cross(v01, v11, i, j + 1, i + 1, j + 1); // top
        cross(v00, v01, i, j, i, j + 1); // left

        if (pts.length < 2) continue;
        // With four crossings the cell is ambiguous; pairing them in order is
        // a standard resolution and is visually indistinguishable at this
        // contour density.
        for (let p = 0; p + 1 < pts.length; p += 2) {
          const [a, b] = [pts[p], pts[p + 1]];
          for (const sign of [1, -1]) {
            ctx.moveTo(gz(a[0]), T.sy(sign * grid.rAt(a[1])));
            ctx.lineTo(gz(b[0]), T.sy(sign * grid.rAt(b[1])));
          }
        }
      }
    }
  }

  ctx.stroke();
  ctx.restore();
}

/**
 * Draw electrodes from the grid itself rather than from the geometry
 * description, so any painted geometry renders without extra code.
 *
 * Nodes are merged into horizontal runs before filling. Drawing one rectangle
 * per node leaves hairline seams between them at fractional scales.
 */
function drawElectrodes(width, height, T) {
  const { grid } = model;
  const { nz, nr } = grid;
  const h = grid.step;

  ctx.save();
  ctx.fillStyle = cssVar('--electrode');

  for (let j = 0; j < nr; j++) {
    let runStart = -1;
    for (let i = 0; i <= nz; i++) {
      const isEl =
        i < nz &&
        grid.electrodeId[j * nz + i] !== NO_ELECTRODE &&
        // The domain end faces close the Laplace problem but are not hardware;
        // drawing them would put a wall across the beam's entrance and exit.
        !(i === 0 && grid.openFaces.zMin) &&
        !(i === nz - 1 && grid.openFaces.zMax);

      if (isEl && runStart === -1) runStart = i;
      if (!isEl && runStart !== -1) {
        const z0 = grid.zAt(runStart) - h / 2;
        const z1 = grid.zAt(i - 1) + h / 2;
        const r0 = grid.rAt(j) - h / 2;
        const r1 = grid.rAt(j) + h / 2;
        for (const sign of [1, -1]) {
          const yA = T.sy(sign * r0);
          const yB = T.sy(sign * r1);
          ctx.fillRect(
            T.sx(z0),
            Math.min(yA, yB),
            T.sx(z1) - T.sx(z0),
            Math.abs(yB - yA)
          );
        }
        runStart = -1;
      }
    }
  }

  ctx.restore();
}

/** The optic axis, drawn as a recessive hairline. */
function drawAxis(width, height, T) {
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

/**
 * Draw trajectories.
 *
 * Each ray is stroked twice: once in the surface colour as a slightly wider
 * halo, then in the trajectory hue. The halo separates rays that overlap and
 * keeps them legible over both poles of the diverging field, which is the
 * secondary encoding the palette's colour-vision margin relies on.
 */
function drawTrajectories(width, height, T) {
  const trajColour = cssVar('--traj');
  const halo = cssVar('--surface-1');

  // A stale path was flown through a different field or a different beam, so
  // it is drawn faintly rather than removed: the shape is still useful
  // context while a slider is moving, but it must not read as the answer.
  const dim = stale ? 0.28 : 1;
  const head = animation.head;

  for (const pass of ['halo', 'line']) {
    ctx.save();
    ctx.strokeStyle = pass === 'halo' ? halo : trajColour;
    ctx.lineWidth = pass === 'halo' ? 4 : 2;
    ctx.globalAlpha = (pass === 'halo' ? 0.55 : 1) * dim;
    ctx.lineJoin = 'round';
    ctx.lineCap = 'round';

    for (const traj of trajectories) {
      const last = Math.min(traj.points.length, head);
      if (last < 2) continue;
      ctx.beginPath();
      for (let n = 0; n < last; n++) {
        const p = traj.points[n];
        const x = T.sx(p.z);
        const y = T.sy(p.x);
        if (n === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }
    ctx.restore();
  }

  // The bunch itself, while it is still in flight. Every ray shares a time
  // base, so these dots are all the same instant - a real snapshot of where
  // the ions are, which is the thing a "fly" button is for.
  if (animation.running) {
    ctx.save();
    ctx.fillStyle = trajColour;
    ctx.strokeStyle = halo;
    ctx.lineWidth = 1.5;
    for (const traj of trajectories) {
      const n = Math.min(traj.points.length, head) - 1;
      if (n < 0 || n >= traj.points.length) continue;
      // Only ions still flying at this instant get a marker.
      if (head > traj.points.length) continue;
      const p = traj.points[n];
      ctx.beginPath();
      ctx.arc(T.sx(p.z), T.sy(p.x), 3.5, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
    }
    ctx.restore();
  }

  // Mark where a ray ended on metal. A strike is a real result, not a failure
  // to render, so it gets an explicit symbol rather than a line that stops.
  ctx.save();
  ctx.globalAlpha = dim;
  ctx.fillStyle = cssVar('--electrode-edge');
  for (const traj of trajectories) {
    if (traj.stop !== 'electrode') continue;
    if (head < traj.points.length) continue;
    const p = traj.points[traj.points.length - 1];
    ctx.beginPath();
    ctx.arc(T.sx(p.z), T.sy(p.x), 3, 0, Math.PI * 2);
    ctx.fill();
  }
  ctx.restore();
}

/** Axial scale, in millimetres. */
function drawScale(width, height, T) {
  const { grid } = model;
  ctx.save();
  ctx.fillStyle = cssVar('--text-muted');
  ctx.font = '11px ui-monospace, monospace';
  ctx.textBaseline = 'bottom';

  const totalMm = mToMm(grid.zLength);
  const stepMm = totalMm > 200 ? 50 : totalMm > 80 ? 20 : 10;
  for (let zmm = 0; zmm <= totalMm + 1e-9; zmm += stepMm) {
    const x = T.sx(zmm * 1e-3);
    ctx.fillRect(x, height - 10, 1, 5);
    ctx.textAlign = zmm === 0 ? 'left' : 'center';
    ctx.fillText(`${zmm.toFixed(0)}`, x, height - 12);
  }
  ctx.textAlign = 'right';
  ctx.fillText('z / mm', width - 6, height - 12);
  ctx.restore();
}

function render() {
  if (!model) return;

  const cssWidth = canvas.parentElement.clientWidth;
  const { grid } = model;
  const aspect = grid.zLength / (2 * grid.rLength);
  const cssHeight = Math.max(200, Math.round(cssWidth / aspect));

  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.round(cssWidth * dpr);
  canvas.height = Math.round(cssHeight * dpr);
  canvas.style.height = `${cssHeight}px`;

  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssWidth, cssHeight);
  ctx.fillStyle = cssVar('--surface-1');
  ctx.fillRect(0, 0, cssWidth, cssHeight);

  const T = makeTransform(cssWidth, cssHeight);

  if (inputs.showField.checked) drawPotential(cssWidth, cssHeight);
  if (inputs.showContours.checked) drawContours(cssWidth, cssHeight, T);
  drawAxis(cssWidth, cssHeight, T);
  drawElectrodes(cssWidth, cssHeight, T);
  drawTrajectories(cssWidth, cssHeight, T);
  drawScale(cssWidth, cssHeight, T);
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
    readoutEl.innerHTML = stat('Beam', 'invalid', '', stats.error, true);
    return;
  }

  // Measured from a probe ray near the axis, so it is a property of the lens
  // rather than of the drawn beam.
  const focal =
    stats.paraxial === null
      ? stat(
          'Paraxial focus',
          '—',
          '',
          stats.reflected > 0
            ? 'beam is reflected — this is an ion mirror, not a lens'
            : 'probe ray does not converge'
        )
      : stat(
          'Paraxial focus',
          stats.paraxial.mm.toFixed(1),
          'mm',
          'from lens centre' +
            (stats.paraxial.extrapolated ? ' · extrapolated beyond the grid' : '')
        );

  const aberration =
    stats.aberration === null
      ? stat('Spherical aberration', '—', '', 'needs a transmitted off-axis ray')
      : stat(
          'Spherical aberration',
          stats.aberration.toFixed(2),
          'mm',
          `outermost ray crosses ${Math.abs(stats.aberration).toFixed(2)} mm ` +
            (stats.aberration < 0 ? 'short (under-corrected)' : 'long')
        );

  // The einzel lens's defining property: entrance and exit are at the same
  // potential, so a transmitted ion must leave with the energy it arrived
  // with. Reported as the worst case across the beam.
  // With space charge on, a non-zero net work is CORRECT, not an error: the
  // beam's own field does real work on its ions as it expands, converting the
  // bunch's electrostatic energy into transverse kinetic energy. The einzel
  // lens's no-net-work property is a property of the LENS, and it only holds
  // for a beam whose ions do not interact. Flagging it as a discrepancy in
  // that case would be flagging correct physics.
  const loaded = stats.beamCurrent > 0;
  const energy =
    stats.keIn === null
      ? stat('Net work', '—', '', 'nothing transmitted')
      : stat(
          'Net work',
          stats.worstWork.toFixed(3),
          'eV',
          loaded
            ? `on ${stats.keIn.toFixed(0)} eV in · expected: the beam's own field does work`
            : `worst of ${stats.transmitted} transmitted, on ${stats.keIn.toFixed(0)} eV in · an einzel lens should do none`,
          !loaded && Math.abs(stats.worstWork) > 0.01 * stats.keIn
        );

  // What the loaded beam does, which is not the same as what the lens does
  // once space charge is on.
  const waist =
    stats.waist === null
      ? stat('Beam waist', '—', '', 'no transmitted ray to measure')
      : stats.waist.crossed
        ? stat(
            'Beam waist',
            'crossover',
            '',
            `outer ray reaches the axis at ${stats.waist.atMm.toFixed(1)} mm`
          )
        : stat(
            'Beam waist',
            stats.waist.radiusMm.toFixed(2),
            'mm',
            `at ${stats.waist.atMm.toFixed(1)} mm · space charge prevents a point focus`,
            true
          );

  const fate = [
    `${stats.transmitted} through`,
    stats.reflected ? `${stats.reflected} reflected` : null,
    stats.struck ? `${stats.struck} on metal` : null,
  ]
    .filter(Boolean)
    .join(' · ');

  const driftPct = stats.drift * 100;

  readoutEl.innerHTML = [
    focal,
    stats.beamCurrent > 0 ? waist : aberration,
    energy,
    stat('Transmitted', `${stats.transmitted}/${stats.total}`, '', fate),
    stat(
      'Energy drift',
      driftPct < 0.01 ? '<0.01' : driftPct.toFixed(2),
      '%',
      'worst ½mv² + qφ deviation · grid quality, not step size',
      stats.drift > 0.02
    ),
    stat(
      'Field solve',
      stats.solveMs.toFixed(0),
      'ms',
      `voltage change ${stats.adjustMs.toFixed(1)} ms · flight ${stats.flyMs.toFixed(0)} ms`
    ),
  ].join('');

  if (stats.voltageWarning) {
    readoutEl.innerHTML += stat('Model warning', '!', '', stats.voltageWarning, true);
  }
  if (stats.warnings?.length) {
    readoutEl.innerHTML += stat('Geometry warning', '!', '', stats.warnings[0], true);
  }
}

/* ------------------------------------------------------------------ */
/* update cycle                                                        */
/* ------------------------------------------------------------------ */

let framePending = false;
let resolvePending = false;
let debounceTimer;

/**
 * Schedule an update.
 *
 * The two paths are paced differently because they cost differently. A
 * voltage or beam change is a fast adjust plus a re-fly, cheap enough to run
 * on the next frame while a slider is still moving. A geometry change
 * re-solves Laplace's equation for every electrode, so it is debounced until
 * the slider settles - otherwise dragging the bore radius would queue one full
 * relaxation per pixel.
 *
 * @param {boolean} resolveField Whether the potential array must be rebuilt.
 */
function update(resolveField) {
  resolvePending = resolvePending || resolveField;

  if (resolveField) {
    statusEl.classList.add('busy');
    clearTimeout(debounceTimer);
    debounceTimer = setTimeout(runUpdate, 140);
    return;
  }

  if (framePending) return;
  framePending = true;
  requestAnimationFrame(runUpdate);
}

function runUpdate() {
  framePending = false;
  clearTimeout(debounceTimer);

  const doResolve = resolvePending || !model;
  resolvePending = false;

  try {
    if (doResolve) rebuild();
    applyVoltages();
    render();
    drawReadout();
  } catch (err) {
    readoutEl.innerHTML = stat('Error', 'failed', '', err.message, true);
    console.error(err);
  } finally {
    statusEl.classList.remove('busy');
  }
}

/* ------------------------------------------------------------------ */
/* flying                                                              */
/* ------------------------------------------------------------------ */

/** Mark the drawn trajectories as no longer matching the settings. */
function markStale() {
  if (trajectories.length === 0) return;
  stale = true;
  flyButton.classList.add('stale');
}

/**
 * Launch the beam.
 *
 * The flight is computed in full first and then revealed frame by frame. The
 * alternative - integrating one step per frame - would tie the physics time
 * step to the display refresh rate, which is a good way to get a different
 * trajectory on a different monitor.
 */
function fly() {
  cancelAnimationFrame(animation.frame);

  flyButton.disabled = true;
  statusEl.classList.add('busy');

  requestAnimationFrame(() => {
    try {
      if (!model) rebuild();
      applyVoltages();
      runFlight();

      stale = false;
      flyButton.classList.remove('stale');

      const total = Math.max(...trajectories.map((t) => t.points.length), 0);
      animation = { head: 0, total, running: total > 0, frame: 0 };

      drawReadout();
      animate();
    } catch (err) {
      readoutEl.innerHTML = stat('Error', 'failed', '', err.message, true);
      console.error(err);
      render();
    } finally {
      flyButton.disabled = false;
      statusEl.classList.remove('busy');
    }
  });
}

/** Reveal the computed flight over roughly a second and a half. */
function animate() {
  const FRAMES = 90;
  const perFrame = Math.max(1, Math.ceil(animation.total / FRAMES));

  const tick = () => {
    animation.head += perFrame;
    if (animation.head >= animation.total) {
      animation.head = Infinity;
      animation.running = false;
      render();
      return;
    }
    render();
    animation.frame = requestAnimationFrame(tick);
  };

  animation.frame = requestAnimationFrame(tick);
}

/* ------------------------------------------------------------------ */
/* wiring                                                              */
/* ------------------------------------------------------------------ */

const OUTPUTS = {
  vCentre: (v) => v,
  vOuter: (v) => v,
  rays: (v) => v,
  beamRadius: (v) => parseFloat(v).toFixed(1),
  beamCurrent: (v) => parseFloat(v).toFixed(1),
  boreRadius: (v) => parseFloat(v).toFixed(1),
  centreLength: (v) => v,
  cfl: (v) => parseFloat(v).toFixed(2),
};

function syncOutputs() {
  for (const [id, fmt] of Object.entries(OUTPUTS)) {
    const out = el(`${id}Out`);
    if (out) out.textContent = fmt(inputs[id].value);
  }
}

for (const [id, input] of Object.entries(inputs)) {
  const needsResolve = GEOMETRY_INPUTS.includes(id);
  const displayOnly = DISPLAY_INPUTS.includes(id);
  const event = input.type === 'range' ? 'input' : 'change';
  input.addEventListener(event, () => {
    syncOutputs();
    // Anything that is not purely cosmetic invalidates the flight already on
    // screen. It stays visible, dimmed, until the beam is flown again.
    if (!displayOnly) markStale();
    update(needsResolve);
  });
}

flyButton.addEventListener('click', fly);

// Space is the obvious key for "go", but only when the user is not part-way
// through typing a number into one of the fields.
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

// Repaint on a theme change so canvas colours follow the CSS tokens.
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', render);

syncOutputs();
update(true);
// Fly once on load so the page is not empty, and so the button's effect is
// obvious before it is pressed.
fly();
