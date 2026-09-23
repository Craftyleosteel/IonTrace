/**
 * Rigid placements in space.
 *
 * A beamline stops being a line the moment it contains a bender, so elements
 * cannot be located by a single axial coordinate any more. Each one carries a
 * FRAME instead: an origin and a set of local axes, saying where it sits and
 * which way it faces.
 *
 * Local coordinates are the same for every element, which is what keeps them
 * interchangeable:
 *
 *     local +z   the direction the beam travels on entry
 *     local +x   one transverse direction ("horizontal", the bend plane)
 *     local +y   the other ("vertical", the bend axis)
 *
 * so an element never needs to know where it is. It is handed a point in its
 * own frame, and returns a field in its own frame; the beamline does the
 * transforming. That is the whole reason a bender can be dropped into the
 * middle of a column without every other element having to care.
 *
 * A frame is stored as an origin `o` and a 3x3 matrix `m` in row-major order
 * whose COLUMNS are the local axes expressed in global coordinates. The matrix
 * is a rotation - orthonormal, determinant +1 - so its inverse is its
 * transpose, and transforming a point back into local coordinates costs the
 * same as transforming it out.
 */

/** The identity placement: at the origin, looking along +z. */
export function identityFrame() {
  return { o: [0, 0, 0], m: [1, 0, 0, 0, 1, 0, 0, 0, 1] };
}

/** A pure translation along the local axes. */
export function translation(x, y, z) {
  return { o: [x, y, z], m: [1, 0, 0, 0, 1, 0, 0, 0, 1] };
}

/**
 * Rotation about the local y axis - a bend in the horizontal plane.
 *
 * Positive `angle` turns the beam towards local -x, which is the convention
 * the bender element uses for a positive bend.
 */
export function yawFrame(angle) {
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  return { o: [0, 0, 0], m: [c, 0, -s, 0, 1, 0, s, 0, c] };
}

/** Rotation about the local x axis - a bend in the vertical plane. */
export function pitchFrame(angle) {
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  return { o: [0, 0, 0], m: [1, 0, 0, 0, c, s, 0, -s, c] };
}

/** A point given in `frame`'s local coordinates, expressed globally. */
export function toGlobal(frame, l) {
  const { o, m } = frame;
  return [
    o[0] + m[0] * l[0] + m[1] * l[1] + m[2] * l[2],
    o[1] + m[3] * l[0] + m[4] * l[1] + m[5] * l[2],
    o[2] + m[6] * l[0] + m[7] * l[1] + m[8] * l[2],
  ];
}

/**
 * A global point expressed in `frame`'s local coordinates.
 *
 * The inverse of a rotation is its transpose, so this is the same cost as the
 * forward transform and carries no accumulated error from inverting.
 */
export function toLocal(frame, g) {
  const { o, m } = frame;
  const dx = g[0] - o[0];
  const dy = g[1] - o[1];
  const dz = g[2] - o[2];
  return [
    m[0] * dx + m[3] * dy + m[6] * dz,
    m[1] * dx + m[4] * dy + m[7] * dz,
    m[2] * dx + m[5] * dy + m[8] * dz,
  ];
}

/**
 * A VECTOR - a field, a velocity - given locally, expressed globally.
 *
 * Distinct from `toGlobal` because a vector has no position: it rotates with
 * the frame but is not displaced by it. Applying the origin to a field would
 * be a category error that happens to compile.
 */
export function vectorToGlobal(frame, v) {
  const { m } = frame;
  return [
    m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
    m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
    m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
  ];
}

/** A global vector expressed in local coordinates. */
export function vectorToLocal(frame, v) {
  const { m } = frame;
  return [
    m[0] * v[0] + m[3] * v[1] + m[6] * v[2],
    m[1] * v[0] + m[4] * v[1] + m[7] * v[2],
    m[2] * v[0] + m[5] * v[1] + m[8] * v[2],
  ];
}

/**
 * `inner` placed relative to `outer`: the composition outer ∘ inner.
 *
 * Chaining these along a column is what makes elements snap together - each
 * one starts where the last one ended, facing wherever the last one left the
 * beam facing.
 */
export function compose(outer, inner) {
  const a = outer.m;
  const b = inner.m;
  const m = new Array(9);
  for (let r = 0; r < 3; r++) {
    for (let c = 0; c < 3; c++) {
      m[r * 3 + c] =
        a[r * 3] * b[c] + a[r * 3 + 1] * b[3 + c] + a[r * 3 + 2] * b[6 + c];
    }
  }
  return { o: toGlobal(outer, inner.o), m };
}

/**
 * A small misalignment, as a real beamline has.
 *
 * Offsets are transverse displacements of the element from where its mount
 * nominally puts it; tilts are angular errors. They are applied to the
 * element's OWN frame and are deliberately not propagated downstream - a
 * misaligned lens does not move the ones after it, because each is mounted
 * independently. Propagating them would model a bent bench rather than a
 * misaligned element.
 *
 * @param {object} align {dx, dy} in metres, {tiltX, tiltY} in radians.
 */
export function misalignment({ dx = 0, dy = 0, tiltX = 0, tiltY = 0 } = {}) {
  const shifted = translation(dx, dy, 0);
  if (tiltX === 0 && tiltY === 0) return shifted;
  return compose(shifted, compose(pitchFrame(tiltX), yawFrame(tiltY)));
}

/** True when a placement is the identity, to within round-off. */
export function isIdentity(frame, tol = 1e-12) {
  const id = identityFrame();
  for (let i = 0; i < 3; i++) if (Math.abs(frame.o[i] - id.o[i]) > tol) return false;
  for (let i = 0; i < 9; i++) if (Math.abs(frame.m[i] - id.m[i]) > tol) return false;
  return true;
}

/** The direction the beam is travelling in this frame, globally. */
export function forwardOf(frame) {
  const { m } = frame;
  return [m[2], m[5], m[8]];
}
