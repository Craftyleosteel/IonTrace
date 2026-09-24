/**
 * Quadrupole deflector - four curved electrodes that turn the beam through a
 * right angle.
 *
 * This is the device a real beamline uses to fold its path: four electrodes
 * filling the corners of a grounded box, with concave faces forming a
 * near-circular aperture, held at +V, -V, +V, -V around the square. The beam
 * enters through the gap on one axis and leaves through the gap on the
 * perpendicular one.
 *
 * The shape of the electrodes
 * ---------------------------
 * Each one is a SQUARE BLOCK filling its quadrant, with a circular arc of
 * radius r0 bitten out of the inner corner and a straight-sided channel of
 * width `channelWidth` separating it from each neighbour. Four grounded posts
 * stand on the diagonals, in the corners, and each block is cut back to clear
 * them.
 *
 * Only one of those details reaches the beam, and it is not the obvious one
 * ------------------------------------------------------------------------
 * The aperture is bounded by the four electrode arcs and the four channel
 * mouths, and an electrode is a CONDUCTOR. Nothing behind its surface can
 * influence the field in front of it. That single fact settles most of the
 * geometry:
 *
 *   - Filling the corners changes nothing inside. The arc is still at +-V and
 *     the aperture cannot tell what is behind it. The blocks are the right
 *     shape because that is the instrument, not because the metal does
 *     anything.
 *   - The corner posts change nothing inside, for the same reason: they stand
 *     behind the blocks, which shield them completely. Measured at the shipped
 *     proportions, the potential 8 mm off axis on the diagonal is 37.28 V with
 *     posts and 37.28 V without - identical to every digit the test prints.
 *     They are structure, not field shaping, and the test below pins that.
 *   - The CHANNEL WIDTH is what reaches the beam, because it is the one choice
 *     that changes the arcs themselves: a wider channel cuts them shorter, and
 *     a shorter arc is a weaker quadrupole.
 *
 * The arc coverage, and what it costs
 * -----------------------------------
 * Worth making quantitative, because it is the whole voltage calibration.
 * Each arc spans a half-angle alpha = pi/4 - psi0 about its diagonal, with
 * psi0 = asin(w/r0) the angle at which the channel wall cuts the circle.
 * Projecting the boundary onto sin(2.theta) gives
 *
 *     phi = (4V/pi) sin(2.alpha) (r/r0)^2 sin(2.theta)
 *         = F . V . 2XZ/r0^2,        F = (4/pi) cos(2.psi0),
 *
 * so F is the factor by which the real electrode beats - or misses - the
 * ideal normalisation the matched voltage assumes, and the voltage needed
 * scales as 1/F. Fully covered boundary (w -> 0): F = 4/pi = 1.273. At the
 * shipped 10.5 mm channel: psi0 = 16.0 deg, F = 1.080. At a 16 mm channel:
 * psi0 = 24.9 deg, F = 0.822 - a THIRD more voltage for the same bend.
 *
 * That estimate treats the open channel mouths as sitting at zero, which they
 * do not, so it is a guide to the scaling rather than a number to quote. It
 * was enough to explain a suite-wide failure: widening the default channel
 * from an effective 5.24 mm half-width to 8 mm dropped F by 24 per cent, and
 * every test that flew at the ideal voltage under-bent into an electrode.
 *
 * Why it is not a sector bender
 * -----------------------------
 * A cylindrical sector deflects by pushing the beam radially inward with a
 * field that is everywhere perpendicular to the orbit. A quadrupole deflector
 * does something else entirely: it sets up a two-dimensional QUADRUPOLE field
 * IN the bend plane and lets it couple the two axes. Writing the potential as
 *
 *     phi = C x z
 *
 * with the beam entering along z and leaving along x, the equations of motion
 * are coupled rather than central:
 *
 *     m x'' = -qC z        m z'' = -qC x
 *
 * so the axial velocity is traded for transverse velocity as the ion crosses.
 * The trajectory is NOT a circular arc.
 *
 * Diagonalising with u = x + z and w = x - z separates them:
 *
 *     u'' = -b u     (oscillatory)
 *     w'' = +b w     (exponential),       b = qC/m
 *
 * The exponential branch is the reason this device is touchy: half of the
 * motion grows rather than oscillates, so an ion that is off-energy or
 * off-axis diverges from the design orbit instead of merely lagging it.
 *
 * The matched voltage
 * -------------------
 * Requiring the ion to arrive at the exit aperture with BOTH the right
 * position and the right direction fixes the design point. Writing s for half
 * the transit phase, the oscillatory branch demands v/(a.omega) = cot(s) and
 * the exponential branch demands v/(a.omega) = tanh(s), so
 *
 *     cot(s) = tanh(s)      =>      s = 0.937552,  theta = 1.875104
 *
 * whence (a.omega/v)^2 = coth^2(s) = 1.85565 and
 *
 *     V0 = 1.8556 (T/q) (r0/a)^2
 *
 * with T the kinetic energy, r0 the aperture radius and a the half-width of
 * the field region.
 *
 * There is a second root, and it is a trap
 * ----------------------------------------
 * cot(s) = tanh(s) has further roots, the next near s = 2.347 (k = 0.9641),
 * and direct integration of the coupled equations confirms it is also a clean
 * ninety degrees - the ion simply takes a longer way round, bending the other
 * way. It is not usable here. On that branch the oscillatory amplitude is
 * 1.427 a rather than 1.241 a, which carries the orbit out to rho = 1.01 a:
 * past the electrode faces at r0, which is necessarily less than a. The ion
 * lands on an electrode instead of reaching the exit. The short branch peaks
 * at rho = 0.88 a and clears the aperture with room to spare, so that is the
 * one this element is designed around.
 *
 * What the formula assumes, and therefore what to expect
 * ------------------------------------------------------
 * It assumes the ideal quadrupole potential everywhere inside the box and
 * nothing outside it. The solved field is neither: real electrodes subtend
 * finite arcs, the grounded box shapes the field near the apertures, grounded
 * posts interrupt the diagonals, and the field does not stop abruptly at the
 * entrance plane. Every one of those pushes the voltage that turns the beam
 * through exactly ninety degrees away from V0, by a few per cent, in a
 * direction that depends on the proportions.
 *
 * What the two changes above do to that number, in opposite directions:
 * filling the electrode corners strengthens the field per volt, so less
 * voltage is needed; the grounded posts weaken it, so more is. Which wins at
 * a given set of proportions is a question for a measurement, and the
 * calibration table that used to sit here was measured on the previous
 * electrode shape - thin annular arcs with wedge-shaped gaps - so it has been
 * removed rather than left to mislead. Re-measure it with `npm test`, or with
 * a sweep of V/V0 against the turn angle, before quoting numbers.
 *
 * So the matched voltage is a starting point of the right size, not the final
 * answer - which is exactly why a real deflector of this kind is followed by a
 * correction lens, and why this element ships with a tuner that searches the
 * solved field for the voltage that actually maximises transmission.
 */

import { PotentialArray, PLANAR } from '../grid.js';
import { solveBasis } from '../laplace.js';
import { Field } from '../field.js';
import { mmToM } from '../constants.js';
import { yawFrame, compose, translation, rollFrame, inverse } from '../frames.js';

export const BENDER_DEFAULTS = {
  apertureRadius: 19, // mm, r0 - centre to the concave electrode faces
  electrodeThickness: 9, // mm, metal from the arc to the flat back face
  boxClearance: 1, // mm, gap between electrode backs and the grounded box
  /*
    Full width of the four straight beam channels.

    This is the parameter that sets the operating voltage, so its default is
    not a round number chosen for looks. The channel edge cuts the electrode
    arc short, and the arc's angular extent is what sets the quadrupole
    strength (see the form factor in the header). 10.5 mm puts the edge at
    asin(5.25/19) = 16.0 degrees, the same coverage as the calibrated default
    this element shipped with before the electrodes were reshaped - so the
    matched voltage carries across the change instead of moving 30 per cent.
  */
  channelWidth: 10.5, // mm
  cornerSize: 5, // mm, the grounded corner posts (0 leaves them out)
  height: 30, // mm, aperture perpendicular to the bend plane
  voltage: 0, // V on each electrode (+V and -V on the diagonals)
  bendPlane: 0, // degrees of roll: 0 bends horizontally, 90 vertically
  gridStep: 0.5, // mm, in the bend plane
};

/** The first root of cot(s) = tanh(s), and the constant it gives. */
export const DEFLECTOR_DESIGN_PHASE = 0.937552;
export const DEFLECTOR_CONSTANT = Math.tanh(DEFLECTOR_DESIGN_PHASE) ** -2; // 1.85565

/**
 * Ideal matched voltage: V0 = k (T/q) (r0/a)^2.
 *
 * `a` is the half-width of the field region - the distance from the centre to
 * the entrance aperture - not the electrode radius. The two differ by the
 * electrode thickness and the box clearance, and the ratio enters squared, so
 * the distinction is worth a factor of two in practice.
 */
export function matchedVoltage(params, energyEV, chargeStates = 1) {
  const p = { ...BENDER_DEFAULTS, ...params };
  const r0 = p.apertureRadius;
  const a = r0 + p.electrodeThickness + p.boxClearance;
  const q = Math.abs(chargeStates);
  if (q === 0 || a === 0) return 0;
  return DEFLECTOR_CONSTANT * (energyEV / q) * (r0 / a) ** 2;
}

export function createBender(params = {}, solverOpts = {}) {
  const p = { ...BENDER_DEFAULTS, ...params };
  const warnings = [];

  const r0 = mmToM(p.apertureRadius);
  const outer = mmToM(p.apertureRadius + p.electrodeThickness);
  const a = mmToM(p.apertureRadius + p.electrodeThickness + p.boxClearance);
  const halfHeight = mmToM(p.height) / 2;
  const w = mmToM(p.channelWidth) / 2; // channel HALF-width, used throughout
  const clear = mmToM(p.boxClearance);
  const post = mmToM(Math.max(0, p.cornerSize));

  if (p.channelWidth <= 0) {
    throw new Error('The deflector needs beam channels to get ions in and out');
  }
  if (w >= outer) {
    throw new Error('Beam channels this wide leave no electrode between them');
  }
  if (w >= r0) {
    warnings.push(
      `Channels ${p.channelWidth} mm wide are wider than the ${2 * p.apertureRadius} mm ` +
        'aperture, so the concave faces are cut away entirely and the aperture is ' +
        'square. The field will be a poor quadrupole.'
    );
  }
  /*
    Square blocks put a flat back face parallel to the box wall along its whole
    width; the old curved backs only came that close at a single point. If the
    clearance is a step or less, the electrode node and the box node are
    NEIGHBOURS - both fixed, with nothing relaxing between them - so the gap
    carries no solved field at all and the electrode is, as far as the solve is
    concerned, sitting against the box. It is not an error and nothing detects
    it downstream, which is the reason to say so here.
  */
  if (p.boxClearance < 2 * p.gridStep) {
    warnings.push(
      `A box clearance of ${p.boxClearance} mm is under two grid steps, so the gap ` +
        'behind each electrode has no free nodes in it and is not resolved at all. ' +
        'Widen the clearance or refine the grid.'
    );
  }

  // Bend plane, as a roll about the beam direction. One solve serves every
  // plane: the field is the same, seen from a rotated frame.
  const roll = (p.bendPlane * Math.PI) / 180;
  const cosR = Math.cos(roll);
  const sinR = Math.sin(roll);
  const intoBend = (x, y) => [cosR * x + sinR * y, -sinR * x + cosR * y];
  const outOfBend = (u, v) => [cosR * u - sinR * v, sinR * u + cosR * v];

  /**
   * The solve lives in the BEND PLANE, and the beam travels inside it.
   *
   * This is the only element here where that is so. A lens is solved in the
   * r-z plane it is symmetric about; a mass filter in the plane transverse to
   * the beam. Here the two-dimensional field and the trajectory occupy the
   * SAME plane, which is precisely why the axes couple.
   *
   * Node count forced odd so the box is exactly symmetric about the centre -
   * without that the four electrodes do not see the same enclosure and the
   * deflector is lopsided.
   */
  const step = mmToM(p.gridStep);
  const half = Math.max(6, Math.round(a / step));
  const extent = half * step;
  const grid = new PotentialArray({
    nz: 2 * half + 1,
    nr: 2 * half + 1,
    step,
    symmetry: PLANAR,
    z0: -extent,
    r0: -extent,
  });

  // The grid's z axis is the transverse bend coordinate, not the beam axis,
  // so its end faces are box walls rather than apertures. Leaving them open -
  // the default, which suits every element solved in the r-z plane - would
  // make the entrance wall invisible to the collision test. The real holes in
  // the box are carved out by `inAperture` below.
  grid.openFaces = { zMin: false, zMax: false };

  const box = grid.addElectrode('box');
  const poleA = grid.addElectrode('poleA'); // the x.z > 0 diagonal
  const poleB = grid.addElectrode('poleB'); // the x.z < 0 diagonal
  grid.paintEnclosure(box);

  /*
    The corner posts.

    They sit in the corners of the box, on the diagonals, which is the one
    place a tie rod can run without crossing a beam channel. Being grounded
    they cannot touch the electrodes, so each block's diagonal corner is cut
    back by a clearance to make room.

    Neither the post nor that cut changes the field the beam sees. Both are
    behind the electrode surface, and the block shields the aperture from them
    completely - which is measured, not assumed: see the `corner posts` tests.
    They are here because the instrument has them and because an ion that
    wanders into the corner should hit grounded metal rather than an electrode.
    `cornerSize: 0` removes them and restores the full square block.
  */
  const postInner = extent - post;
  const blockCut = postInner - clear; // block corners stop here

  if (post > 0 && blockCut <= w) {
    throw new Error('Corner posts this large leave no electrode between the channels');
  }

  /**
   * One electrode: a square block filling its quadrant, with a concave
   * circular face of radius r0 bitten out of the inner corner, straight-sided
   * beam channels of half-width `w` along both axes, and its diagonal corner
   * cut back to clear the grounded post.
   *
   * The concave face is what makes the field quadrupolar near the centre; the
   * filled corners are what make it strong. An annular arc of the same inner
   * radius leaves the quadrant behind it empty, and that emptiness is a region
   * the grounded box reaches into.
   *
   * The grid's "z" axis carries the axial coordinate Z and its "r" axis the
   * transverse X, both measured from the deflector centre.
   */
  const tol = step * 1e-6;
  const quadrant = (sx, sz) => (Z, X) => {
    if (Math.sign(X) !== sx || Math.sign(Z) !== sz) return false;
    const ax = Math.abs(X);
    const az = Math.abs(Z);
    if (ax < w || az < w) return false; // the beam channels
    if (ax > outer + tol || az > outer + tol) return false; // the flat back faces
    if (Math.hypot(X, Z) < r0 - tol) return false; // the concave face
    if (post > 0 && ax > blockCut && az > blockCut) return false; // room for the post
    return true;
  };

  // Sign convention: a POSITIVE voltage bends a POSITIVE ion towards -x. That
  // choice is what ties the electrode polarity to the exit face, and the two
  // are not independent - with the diagonals swapped the same ion turns the
  // other way and leaves through the +x wall instead, which this element has
  // no aperture in. A device wired the other way round, as the published one
  // in the reference drawing is, is the mirror image of this and identical
  // physics.
  const paintedA =
    grid.paint(poleA, quadrant(1, 1)) + grid.paint(poleA, quadrant(-1, -1));
  const paintedB =
    grid.paint(poleB, quadrant(-1, 1)) + grid.paint(poleB, quadrant(1, -1));
  if (paintedA === 0 || paintedB === 0) {
    throw new Error('Deflector electrodes covered no grid nodes; check the geometry');
  }

  // The posts, grounded along with the box they stand in the corners of. They
  // get no voltage of their own, so they add no basis solution and cost
  // nothing to paint - but they are metal, and the ion collision test and the
  // Laplace solve both see them.
  if (post > 0) {
    grid.paint(box, (Z, X) => Math.abs(X) >= postInner && Math.abs(Z) >= postInner);
  }

  const { basis, reports } = solveBasis(grid, solverOpts);
  if (reports.some((r) => !r.converged)) {
    warnings.push('Laplace solve did not converge; this field is not trustworthy.');
  }

  // Unit solution, per volt on each electrode. Everything the element does is
  // this map scaled, so changing the voltage never re-solves.
  const unit = new Field(grid, basis);
  unit.setVoltages({ box: 0, poleA: -1, poleB: 1 });

  // Nominal path: a quarter turn of radius `a` from the entrance aperture to
  // the exit one. The REAL trajectory is not this curve - it is not a
  // circular arc at all - but the nominal orbit is what the element's
  // placement and the drawn reference line are built on.
  const length = (Math.PI / 2) * a;

  /** Local point -> deflector-centred bend coordinates. */
  function centred(x, y, z) {
    const [bx, by] = intoBend(x, y);
    return { X: bx, Z: z - a, h: by };
  }

  /**
   * The beam holes in the grounded box.
   *
   * The solve paints the box all the way round, because a grounded plane with
   * a hole in it is very nearly a grounded plane and solving it that way keeps
   * the Dirichlet boundary closed. But the holes are real, and an ion passing
   * through one must not be counted as hitting metal. Without this the failure
   * is total rather than subtle: every ion is destroyed the instant it reaches
   * the entrance plane, at every voltage.
   *
   * The hole is as wide as the channel it feeds, and no guesswork is needed to
   * say how wide that is: the channels are straight-sided, so the hole is
   * exactly `channelWidth` across. That is the point of a straight channel
   * over a wedge - the clear width is the same where the beam enters, where it
   * passes the electrodes and where it leaves, so one number describes it.
   *
   * There are FOUR of them, one on each face. That is not generosity, it is
   * the symmetry of the device: the electrodes stop short of both axes in both
   * directions, so the box has a clear channel along each. It is also what
   * makes the thing a switch rather than a corner. Turned off, the beam goes
   * in one face and straight out of the opposite one; turned on, it leaves
   * through a side. Cutting only the two holes the bend needs would wall off
   * the straight path and quietly destroy every ion that took it.
   */
  const holeHalf = w;
  const rim = step * 1.5;
  const inAperture = (X, Z) =>
    (Math.abs(Z) > extent - rim && Math.abs(X) < holeHalf) ||
    (Math.abs(X) > extent - rim && Math.abs(Z) < holeHalf);

  return {
    type: 'bender',
    label: 'Quadrupole deflector',
    params: p,
    length,
    bore: r0,
    // Deliberately no `clearBore`. The electrodes wrap around the beam rather
    // than lying outside a cylinder about the element's axis, so a point can
    // be well within r0 of that axis and still be buried in metal - a point at
    // 45 degrees, for instance, sits on an electrode face while being only
    // 14 mm from the axis of a 19 mm aperture. Distance from the axis is no
    // guide here and the real test has to run.
    outerRadius: Math.max(extent, halfHeight),
    lengthScale: step,
    shortestPeriod: null,
    warnings,
    grid,
    field: unit,
    curved: true,

    /** A quarter turn, conjugated by the roll so the beam is not twisted. */
    exitTransform: compose(
      compose(rollFrame(roll), compose(translation(-a, 0, a), yawFrame(Math.PI / 2))),
      inverse(rollFrame(roll))
    ),

    /**
     * THREE ways out, which is what makes this element a junction.
     *
     * The box has an aperture on each of its four faces, because the
     * electrodes stop short of both axes in both directions. One is the
     * entrance. The other three are all reachable, and which one the beam
     * takes is decided by the voltage alone:
     *
     *     +V   bent one way      out the -x face
     *      0   straight through  out the +z face
     *     -V   bent the other    out the +x face
     *
     * That is not a modelling convenience, it is what these devices are for: a
     * switch that sends a beam down one of three lines without moving any
     * hardware. Reversing the polarity mirrors the whole problem in x, so the
     * counter-bend is the same physics at the same matched magnitude - there
     * is no second voltage to find.
     *
     * A column can therefore have hardware bolted to all three, and the beam
     * goes where the field sends it. A port with nothing attached is simply an
     * opening the beam may leave through.
     *
     * The straight path is 2a long - in at one face, out at the opposite one -
     * while each bent path is a quarter arc of radius a. Different distances,
     * so each exit carries its own.
     */
    exits: [
      {
        port: 'bend',
        label: 'Bent',
        length,
        // The quarter arc, so the drawn path curves through the box rather
        // than cutting the corner.
        path: (f) => {
          const t = (f * Math.PI) / 2;
          const [x, y] = outOfBend(-a + a * Math.cos(t), 0);
          return [x, y, a * Math.sin(t)];
        },
        transform: compose(
          compose(rollFrame(roll), compose(translation(-a, 0, a), yawFrame(Math.PI / 2))),
          inverse(rollFrame(roll))
        ),
      },
      {
        port: 'straight',
        label: 'Straight through',
        length: 2 * a,
        transform: translation(0, 0, 2 * a),
      },
      {
        port: 'counter',
        label: 'Bent the other way',
        length,
        // The mirror image of the bend arc, which is exactly what reversing
        // the polarity produces.
        path: (f) => {
          const t = (f * Math.PI) / 2;
          const [x, y] = outOfBend(a - a * Math.cos(t), 0);
          return [x, y, a * Math.sin(t)];
        },
        transform: compose(
          compose(rollFrame(roll), compose(translation(a, 0, a), yawFrame(-Math.PI / 2))),
          inverse(rollFrame(roll))
        ),
      },
    ],

    pathPoint(f) {
      const t = (f * Math.PI) / 2;
      const [x, y] = outOfBend(-a + a * Math.cos(t), 0);
      return [x, y, a * Math.sin(t)];
    },

    setVoltage(v) {
      // Nothing to recompute: fieldAt scales the stored unit map.
      p.voltage = v;
    },

    contains(x, y, z) {
      const { X, Z, h } = centred(x, y, z);
      // A tolerance of one grid step: the entrance plane sits exactly on the
      // box edge, and without it an ion arriving there would belong to no
      // element at all and be reported as having wandered out of the column.
      const edge = extent + step;
      return Math.abs(X) <= edge && Math.abs(Z) <= edge && Math.abs(h) <= halfHeight * 1.5;
    },

    fieldAt(x, y, z) {
      const { X, Z } = centred(x, y, z);
      if (Math.abs(X) > extent || Math.abs(Z) > extent) return { Ex: 0, Ey: 0, Ez: 0 };
      // The grid's axes are (Z, X); there is no field along the height for an
      // ideal deflector, because the electrodes are uniform in that direction.
      const { Ez: eZ, Er: eX } = unit.fieldAt(Z, X);
      const [ex, ey] = outOfBend(p.voltage * eX, 0);
      return { Ex: ex, Ey: ey, Ez: p.voltage * eZ };
    },

    potentialAt(x, y, z) {
      const { X, Z } = centred(x, y, z);
      if (Math.abs(X) > extent || Math.abs(Z) > extent) return 0;
      return p.voltage * unit.potentialAt(Z, X);
    },

    strikes(x, y, z) {
      const { X, Z, h } = centred(x, y, z);
      if (Math.abs(h) > halfHeight) return true;
      if (inAperture(X, Z)) return false;
      if (Math.abs(X) > extent || Math.abs(Z) > extent) return true;
      // Field.strikes takes (transverse, height, axial). The grid's axial axis
      // carries Z and its transverse axis carries X, so they go in that order -
      // the reverse of fieldAt(Z, X), which takes (axial, transverse).
      return unit.strikes(X, 0, Z);
    },

    /** Ideal matched voltage for this ion, in volts on each electrode. */
    matchedVoltage(energyEV, chargeStates) {
      return matchedVoltage(p, energyEV, chargeStates);
    },

    /**
     * Four curved electrodes plus the box, as polygons in the element's own
     * frame. The concave faces are drawn as real arcs, because their shape is
     * what makes the field quadrupolar rather than a detail of the picture.
     */
    shapes() {
      // A point in the deflector's own centred coordinates, put back into the
      // element's frame - rolled into the bend plane, and shifted so the
      // centre sits one half-width along the beam.
      const corner = (X, Z) => {
        const [x, y] = outOfBend(X, 0);
        return [x, y, a + Z];
      };

      const out = [];
      const arcSteps = 12;

      /*
        The angle at which the concave face meets a channel wall. The wall is
        the straight line X = w, so the face runs from where that line cuts the
        circle round to its mirror image at Z = w. Wider channels eat the face
        from both ends; once the channel is wider than the aperture there is no
        face left and the inner corner is simply square.
      */
      const faceEnds = w < r0;
      const psi0 = faceEnds ? Math.asin(w / r0) : 0;
      const span = Math.PI / 2 - 2 * psi0;
      const notched = post > 0 && blockCut < outer;

      for (const [sx, sz] of [[1, 1], [-1, 1], [1, -1], [-1, -1]]) {
        const at = (X, Z) => corner(sx * X, sz * Z);
        const pts = [];

        if (faceEnds) {
          // The concave face, from the channel wall at X = w round to Z = w.
          for (let k = 0; k <= arcSteps; k++) {
            const psi = psi0 + (k / arcSteps) * span;
            pts.push(at(r0 * Math.sin(psi), r0 * Math.cos(psi)));
          }
        } else {
          pts.push(at(w, w));
        }

        // Out along one channel wall, round the flat back, home along the
        // other - with the diagonal corner cut square where a post stands.
        pts.push(at(outer, w));
        if (notched) {
          pts.push(at(outer, blockCut), at(blockCut, blockCut), at(blockCut, outer));
        } else {
          pts.push(at(outer, outer));
        }
        pts.push(at(w, outer));

        out.push({ points: pts });
      }

      // The posts themselves, drawn as the grounded structure they are.
      if (post > 0) {
        for (const [sx, sz] of [[1, 1], [-1, 1], [1, -1], [-1, -1]]) {
          const at = (X, Z) => corner(sx * X, sz * Z);
          out.push({
            points: [
              at(postInner, postInner),
              at(extent, postInner),
              at(extent, extent),
              at(postInner, extent),
            ],
            wall: true,
          });
        }
      }

      // The grounded box, as four thin walls.
      const t = Math.max(step * 2, a * 0.03);
      const wall = (x0, x1, z0, z1) => ({
        points: [corner(x0, z0), corner(x1, z0), corner(x1, z1), corner(x0, z1)],
        wall: true,
      });
      // All four faces in two pieces each, with the beam hole between them, so
      // the picture shows the same apertures the collision test uses.
      for (const near of [true, false]) {
        const a0 = near ? -extent : extent - t;
        const a1 = near ? -extent + t : extent;
        out.push(wall(-extent, -holeHalf, a0, a1)); // a Z face
        out.push(wall(holeHalf, extent, a0, a1));
        out.push(wall(a0, a1, -extent, -holeHalf)); // an X face
        out.push(wall(a0, a1, holeHalf, extent));
      }

      return out;
    },
  };
}
