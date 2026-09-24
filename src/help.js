/**
 * The help topics behind every question mark in the interface.
 *
 * Why this file exists at all
 * ---------------------------
 * The explanations were already written - they are the docstrings at the top
 * of element.js, column.js, optimize.js and the rest, and they are the most
 * carefully checked prose in the project because they are where the physics is
 * justified. But they were readable only by someone with the source open, and
 * the person who most needs to know why a deflector wants 1.22 times its
 * closed-form voltage is the person staring at the deflector, not at the file.
 *
 * So this is not a second set of documentation to keep in step with the first.
 * It is the same explanations, in the shortest form that survives being read
 * in a popover, with the long derivation left where it belongs and pointed at.
 *
 * What a topic owes the reader
 * ----------------------------
 * Three things, in this order: what the control DOES, what it costs or breaks,
 * and the number that makes it concrete. The last is what separates this from
 * a tooltip. "Sets the deflector voltage" helps nobody; "the closed form runs
 * 20 to 40 per cent low because it knows only the aperture, and the tuner is
 * what finds the rest" is the thing worth knowing.
 *
 * Numbers quoted here are MEASURED against this solver, not quoted from a
 * textbook, and where a measurement has not been made the text says so rather
 * than rounding a guess into authority.
 */

/**
 * @typedef {object} Topic
 * @property {string} title  Short heading for the popover.
 * @property {string[]} body One paragraph per entry.
 * @property {string} [more] Where the full argument lives.
 */

/** @type {Record<string, Topic>} */
export const HELP = {
  /* ---------------------------------------------------------------- */
  /* how the simulation works                                          */
  /* ---------------------------------------------------------------- */

  laplace: {
    title: 'How the field is solved',
    body: [
      'Every electrode is painted onto a grid, and the potential between them is found by solving Laplace’s equation — the statement that in empty space the potential at a point is the average of its neighbours. The solver sweeps the grid until that is true everywhere.',
      'Nothing about the beam enters this. The field is a property of the metal and its voltages alone, which is why it can be computed once and reused for every ion.',
      'Axisymmetric elements are solved on a two-dimensional (z, r) grid and rotation supplies the third dimension exactly. A quadrupole is not axisymmetric, so it is solved on a transverse (x, y) grid instead, and a deflector in the plane its beam bends in.',
    ],
    more: 'docs/PHYSICS.md §2',
  },

  fastAdjust: {
    title: 'Why changing a voltage is instant',
    body: [
      'Each electrode is solved once at one volt, with every other electrode grounded. Because Laplace’s equation is linear, the field at any voltage is just those solutions scaled and added.',
      'So moving a voltage slider is a multiplication, not a re-solve. Changing the geometry is a re-solve, which is why those controls feel different.',
      'This is exact over voltages on a fixed set of conductors. It is not valid over geometry: adding a grounded plate next to a lens changes that lens’s field completely, and no amount of adding solutions together will show it.',
    ],
    more: 'docs/PHYSICS.md §3.2',
  },

  integrator: {
    title: 'How the ions are flown',
    body: [
      'Motion follows m·d²r/dt² = q·E(r) in three dimensions, integrated with Runge–Kutta 4 by default. Velocity Verlet is offered because it conserves energy better over long flights; RK4 is more accurate per step.',
      'The step is adaptive: it shrinks where the field varies quickly. Step fraction sets how cautious that is — smaller is more accurate and slower. 0.05 is a good default.',
      'Magnetic forces, buffer-gas collisions, image charges and relativistic correction are all absent. For ions at the energies modelled here the relativistic error is far below the grid error.',
    ],
    more: 'docs/PHYSICS.md §4',
  },

  /* ---------------------------------------------------------------- */
  /* the physics switches                                              */
  /* ---------------------------------------------------------------- */

  repulsion: {
    title: 'Ion–ion repulsion',
    body: [
      'Off, every ion flies as if alone — fast, and right for a weak beam. The two other settings make the ions see each other’s charge.',
      'Beam treats the ions as a continuous current, which is what a real steady beam is. Its spreading matches the closed-form envelope equation to about 0.1 %, so it is the one to trust for a DC beam.',
      'Coulomb treats each flown ion as a clump of real ions and sums the forces between them. It is the honest model for a short pulse, and much slower: the cost grows as the square of the ion count.',
    ],
    more: 'docs/PHYSICS.md §6',
  },

  fringe: {
    title: 'Fringe fields',
    body: [
      'Off, every element is solved alone behind grounded end faces. Those faces are a numerical device — Laplace’s equation needs a closed boundary — but they behave like a grounded plate a few millimetres from the hardware, so they do not merely omit the fringe field, they shield it.',
      'On, neighbouring axisymmetric elements share one grid and the field flows between them. Measured on an aperture plate: with its own caps in place the field outside it reads zero, and it is really 16 % of the plate potential at the boundary, falling to 0.5 % twenty millimetres out.',
      'This changes the field the ions fly through, not just the picture. It costs a Laplace solve per run of elements.',
    ],
    more: 'src/column.js',
  },

  shielding: {
    title: 'What stops a fringe field',
    body: [
      'Inside a grounded pipe of radius R, a disturbance dies along the axis as exp(−2.405·z/R). The decay length is R/2.405 and depends on nothing else — not the bore of whatever made the field, not its length, not its voltage.',
      'So a narrow grounded tube is an excellent shield: 5 mm of bore gives a decay length of 2.1 mm, and a fringe is down by a factor of 1000 within 15 mm. A wide housing is a poor one — 16 mm gives 6.7 mm, and the same 15 mm barely knocks it down by a factor of 9.',
      'Which is why the answer to “how do I keep this field off my detector” is usually an aperture at ground, not more distance.',
    ],
    more: 'docs/PHYSICS.md §10',
  },

  /* ---------------------------------------------------------------- */
  /* tools                                                             */
  /* ---------------------------------------------------------------- */

  optimise: {
    title: 'Optimise voltages',
    body: [
      'Searches every electrode voltage in the column for the combination that delivers the most beam, then refines the best one it found. It is cheap enough to be a button because a voltage is only a multiplier on an already-solved field.',
      'The search is in two stages. A coarse sweep finds the region that transmits at all — without it, sweeping a lens while the deflector is off transmits nothing and the search has no gradient to follow. A Newton step in the flat directions then tightens the beam without losing any of it.',
      'It never returns something worse than it started with: if nothing improves, the original voltages go back.',
    ],
    more: 'docs/PHYSICS.md §11',
  },

  tune: {
    title: 'Tuning one element',
    body: [
      'Searches this element’s own voltages only, leaving the rest of the column alone. Useful when you know which element is wrong.',
      'Optimise voltages, in the toolbar, moves everything at once — better when elements trade off against each other, which a lens and the deflector behind it certainly do.',
    ],
  },

  saveLoad: {
    title: 'Saving and loading',
    body: [
      'A column is saved as JSON: it is a tree, every browser parses it without a library, and you can open the file and read it. A spreadsheet row cannot hold a branching column without either flattening the branches away or sprawling into mostly empty columns.',
      'The file holds the description, not the solution — element types, their settings, how they hang off one another, the ion and the physics switches. No fields, no trajectories. Loading re-solves, so an old file gets today’s physics rather than yesterday’s cached numbers.',
      'A setting the file has that this version no longer does is reported rather than dropped silently, because a column that loads looking fine while describing different hardware is the worst available outcome.',
    ],
    more: 'src/scene.js',
  },

  flowChart: {
    title: 'The beamline chart',
    body: [
      'The column as a tree. A deflector has three ways out, so a beamline is not a list — it branches, and this is where the branching is visible.',
      'Drag a box to move that element and everything below it onto another exit. Click an empty socket to start a new line there.',
      'An exit with nothing attached is not an error. It is an opening the beam may leave through, which is exactly what a switch is for.',
    ],
  },

  align: {
    title: 'Misalignment',
    body: [
      'Moves this element only. The ones after it stay where they are, because each is bolted to its own mount — an offset here does not shift the rest of the line.',
      'This is how you find out what a real assembly tolerance costs you. A tenth of a millimetre on a strong lens is not a small effect.',
    ],
  },

  /* ---------------------------------------------------------------- */
  /* what the view shows                                               */
  /* ---------------------------------------------------------------- */

  potentialMap: {
    title: 'Potential map',
    body: [
      'Colour is the potential: negative one way, positive the other. It shows where the high and low ground is, which is what an ion rolls down.',
      'Where there is no potential there is no colour — the map fades out rather than stopping at an edge. That matters because the region solved is much larger than the hardware in it, so a map painted right to its boundary drew a rectangle around every element that had nothing physical about it.',
      'The dashed lines across the beam are the other half of the same idea: they mark where one solve ends and the next begins, so they appear only between elements solved separately. Turn on fringe fields and the ones inside a shared solve disappear, because there is no longer a seam there.',
    ],
  },

  equipotentials: {
    title: 'Equipotentials',
    body: [
      'Lines of constant potential. They bunch together where the field is strong, and an ion crossing them gains or loses exactly the energy of the interval it crossed, regardless of the path it took.',
    ],
  },

  fieldLines: {
    title: 'Field lines',
    body: [
      'Lines along E, running from positive metal to negative and crossing equipotentials at right angles. They are the clearest way to see a fringe field reaching out of an element — or being stopped by a grounded surface.',
      'They are not trajectories. An ion has inertia and generally does not follow one.',
      'Whether a line runs on into the next element tells you something real. With fringe fields off, each element is solved behind its own grounded end faces, so the field genuinely stops at the boundary and so does the line. Turn fringe fields on and neighbouring elements share one solve — the map runs unbroken and the lines carry through the joins, because there is no longer a boundary there to stop at.',
      'A deflector still breaks the chain either side of it: it is not a body of revolution, so it cannot share a grid with its neighbours.',
    ],
  },

  /* ---------------------------------------------------------------- */
  /* elements                                                          */
  /* ---------------------------------------------------------------- */

  drift: {
    title: 'Drift',
    body: [
      'A grounded tube. No field inside, so an ion crosses it in a straight line — unless the ions repel each other, in which case the beam spreads here and nowhere else is it so easy to see.',
      'Drifts are also what give neighbouring elements room. An element needs roughly three bore radii of margin before its own end caps stop clipping its fringe field.',
    ],
  },

  aperture: {
    title: 'Aperture plate',
    body: [
      'A plate with a hole, at a voltage. Charged, it is a lens — a weak one, and the simplest there is. At ground it is a shield, and a very effective one.',
      'Its field reaches roughly one bore radius either side, so it needs several radii of margin before its own boundary starts clipping it.',
    ],
  },

  tube: {
    title: 'Biased tube',
    body: [
      'A single cylinder at a DC voltage — an einzel lens with its two grounded guard cylinders taken away. Inside it the ion is at a different kinetic energy, which is the point of biasing it.',
      'Across the whole element the energy is unchanged, because the domain is closed by grounded faces at both ends: the ion starts and finishes at the same potential, so the exit fringe returns whatever the entrance fringe took. It focuses; it does not accelerate.',
      'That is a fact about the hardware, not the model. A tube only changes a beam’s energy if what comes after it sits at a different potential — so an accelerating stage needs this tube and a downstream element, solved together with fringe fields on.',
    ],
  },

  electrode: {
    title: 'Conductor',
    body: [
      'One piece of metal, any shape of revolution, at a DC voltage. Four numbers set the shape: thin in z is a plate, long in z is a tube, an inner radius of zero is a solid disc, and a narrow band of radius is a ring.',
      'It holds one conductor rather than a list of them because the beamline already composes things. Stack several with fringe fields on and they are solved together on one grid — the same answer a multi-electrode element would give, except each piece can be selected, tuned, dragged and saved on its own.',
      'Being a body of revolution is a real restriction. It can be a ring but never a pair of rods or a slit, because anything with corners in the transverse plane could not share an r-z solve with its neighbours — which is most of the point of being able to stack them.',
    ],
  },

  einzel: {
    title: 'Einzel lens',
    body: [
      'Three tubes: grounded, live, grounded. An ion leaves at the energy it arrived with, because the potential it ends in equals the potential it started in — the middle electrode only changes where it goes, never how fast.',
      'That is what makes it the workhorse of electrostatic optics: it focuses without touching the beam energy, so the tuning downstream does not shift when you change it.',
    ],
  },

  bender: {
    title: 'Quadrupole deflector',
    body: [
      'Four electrodes around a square aperture, at +V and −V on the diagonals. It does not push the beam round a corner like a sector — it sets up a quadrupole field in the plane the beam travels in, and lets the two axes trade velocity.',
      'Half of that motion grows exponentially rather than oscillating, which is why a deflector is touchy: an off-energy ion does not merely lag the design orbit, it leaves it. That is also what makes it an energy filter rather than merely a corner.',
      'The closed-form voltage runs 20–40 % low, because it knows only the aperture while a real ion is also kicked by the entrance and exit channels. It always transmits, but it does not turn a right angle — the tuner finds the voltage that does.',
    ],
    more: 'docs/PHYSICS.md §9',
  },

  benderPorts: {
    title: 'A deflector is a switch',
    body: [
      'There is an opening on all four faces, so three of them are ways out. Which one the beam takes is decided by the voltage alone: positive bends it one way, zero sends it straight through, negative bends it the other.',
      'That is what these devices are for — sending a beam down one of three lines without moving any hardware. Attach a line to each port and the voltage chooses between them.',
    ],
  },

  quadrupole: {
    title: 'Quadrupole mass filter',
    body: [
      'Four rods with a DC and an RF voltage on them. Only ions in a band of mass-to-charge survive the crossing; everything else grows without bound and hits a rod.',
      'Stability is set by the Mathieu parameters a and q. The useful edge is q = 0.90803 — past it nothing is stable at any mass. The default works at q = 0.38, comfortably inside.',
    ],
    more: 'docs/PHYSICS.md §8',
  },

  multipole: {
    title: 'Multipole guide',
    body: [
      'An RF-only guide. The rapidly oscillating field averages into an effective potential well that pushes ions toward the axis — deeper with more voltage, shallower with more mass.',
      'More poles make the well flatter in the middle and steeper at the wall, so a hexapole or octopole holds a wider, colder beam than a quadrupole does, and confines it less tightly.',
      'Without buffer gas to take energy out, a deeper well is not automatically better: nothing damps the radial motion the guide puts in.',
    ],
  },

  funnel: {
    title: 'Ion funnel',
    body: [
      'A stack of rings of decreasing bore, carrying an RF voltage that repels ions from the walls and a DC gradient that pushes them along. It squeezes a wide, diffuse cloud into a narrow beam.',
      'A real funnel works in a gas, where collisions damp the radial motion the RF puts in. There is no gas here, so a deeper RF well transmits worse rather than better — measured, 1.1 eV gives 9 of 9 and 17 eV gives none.',
    ],
  },

  detector: {
    title: 'Ion detector',
    body: [
      'A biased collector. A grounded entrance aperture with an active surface a few kilovolts below it, so the field through the aperture pulls ions in and accelerates them onto the surface rather than waiting to be hit.',
      'It distinguishes what it catches from what it loses: an ion landing on the active face is a count, one landing on the housing is a loss. Both are strikes, and reporting them as one number would hide the difference the element exists to draw.',
      'Solved alone, its field stops at its own boundary. With fringe fields on it reaches back up the column, which is what “sucks ions in” actually means.',
    ],
  },

  /* ---------------------------------------------------------------- */
  /* the beam                                                          */
  /* ---------------------------------------------------------------- */

  beamSpec: {
    title: 'The ion',
    body: [
      'Mass in unified atomic mass units, charge in elementary charges, kinetic energy in electron volts. Everything downstream scales off these three.',
      'A negative charge is allowed and is not a cosmetic change: every lens and deflector that worked for a positive ion now does the opposite, so a column tuned for one will not pass the other.',
    ],
  },

  divergence: {
    title: 'Divergence',
    body: [
      'Each ion is aimed outward at an angle proportional to its distance from the axis, so the beam expands as a cone from a waist at the source. Zero is perfectly collimated; the outermost ion gets the full angle.',
      'A real beam always has some. Tuning a column against a perfectly parallel beam flatters it.',
    ],
  },
};

/** A topic by key, or null. */
export function topic(key) {
  return Object.prototype.hasOwnProperty.call(HELP, key) ? HELP[key] : null;
}

/** Every topic key, for tests and for the tutorial's index. */
export function topicKeys() {
  return Object.keys(HELP);
}
