/**
 * The guided walkthrough.
 *
 * Why the steps DO things
 * -----------------------
 * A tutorial that only describes the interface is a manual with worse
 * navigation. Every step here can build the column it is talking about, set the
 * switch it is about to explain, and fly the beam - so the thing on screen is
 * always the thing being discussed, and the reader can stop following and start
 * poking at any point without losing their place.
 *
 * That also keeps the tutorial honest. A step whose `run` builds a column that
 * does not transmit is a step making a claim the simulator does not support,
 * and the test suite flies every one of them for exactly that reason. Prose can
 * drift away from the code; a column that has to fly cannot.
 *
 * What a lesson owes the reader
 * -----------------------------
 * Something they could not have guessed from the control's name. "Set the lens
 * to -300 V" is a manual; "the ion leaves at the energy it arrived with, so
 * focusing here does not detune anything downstream" is the reason einzel
 * lenses are everywhere. Where a number appears it is one measured against this
 * solver, not a plausible-looking round figure.
 *
 * The api handed to `run`
 * -----------------------
 * Deliberately small, and all of it things a person could do by hand with the
 * normal controls: build a column, set the beam, flip the physics switches,
 * choose what the view draws, select something, fly, tune. A step cannot reach
 * past the interface into the model, so nothing here can demonstrate behaviour
 * the controls cannot reproduce.
 */

/**
 * @typedef {object} TutorialApi
 * @property {(specs: {type: string, params?: object}[]) => void} build
 * @property {(spec: object) => void} beam
 * @property {(spec: object) => void} physics
 * @property {(spec: object) => void} view
 * @property {(target: number | 'source' | null) => void} select
 * @property {() => void} fly
 * @property {() => Promise<void>} optimise
 */

const DRIFT = (length = 25, bore = 6) => ({ type: 'drift', params: { length, bore } });

export const LESSONS = [
  /* ---------------------------------------------------------------- */
  {
    id: 'first',
    title: 'A first beamline',
    summary: 'Build a column, fly ions through it, and read what came out.',
    steps: [
      {
        title: 'What this simulator is doing',
        body: [
          'IonTrace solves the electric field from the shapes and voltages of metal, then integrates m·a = q·E to fly ions through it. Nothing is animated or approximated from a formula — the trajectories are the result of the field.',
          'That means it can be wrong in an interesting way rather than an obvious one: a wrong field produces a perfectly plausible-looking beam. It is why nearly every number in this program is checked against a closed form somewhere in the test suite.',
        ],
        topic: 'laplace',
      },
      {
        title: 'The simplest column',
        body: [
          'Every beamline starts with a drift — a grounded tube with no field inside it. An ion crosses it in a straight line.',
          'Press the button below to build one, then Fly ions. The trajectories should be perfectly straight, which is the most boring and most reassuring result available.',
        ],
        action: {
          label: 'Build a drift and fly',
          run: (api) => {
            api.build([DRIFT(60, 8)]);
            api.beam({ mass: 100, charge: 1, energy: 50, rays: 9, beamRadius: 1.5 });
            api.physics({ repulsion: 'none', fringe: false });
            api.fly();
          },
        },
        topic: 'drift',
      },
      {
        title: 'Reading the result',
        body: [
          'Click the beam in the diagram to see what it did — how many ions got through, where they ended up, and how wide the bundle is.',
          'The ion source is selectable like any element, so its settings live where everything else’s do. Mass, charge and energy are the three numbers everything downstream scales off.',
        ],
        action: {
          label: 'Select the beam',
          run: (api) => api.select('source'),
        },
        topic: 'beamSpec',
      },
    ],
  },

  /* ---------------------------------------------------------------- */
  {
    id: 'focus',
    title: 'Focusing the beam',
    summary: 'The einzel lens, and why it does not change the beam energy.',
    steps: [
      {
        title: 'A lens made of three tubes',
        body: [
          'An einzel lens is grounded, live, grounded. An ion is accelerated into the middle electrode and decelerated by exactly the same amount on the way out, so it leaves at the energy it arrived with.',
          'That is the whole reason it is the workhorse of electrostatic optics: it changes where the beam goes without changing how fast it travels, so tuning it does not detune everything downstream.',
        ],
        action: {
          label: 'Build a lens and fly',
          run: (api) => {
            api.build([
              DRIFT(20, 6),
              { type: 'einzel', params: { voltage: -300, boreRadius: 6, housingRadius: 15 } },
              DRIFT(60, 6),
            ]);
            api.beam({ mass: 100, charge: 1, energy: 50, rays: 9, beamRadius: 1.5, divergence: 0 });
            api.fly();
          },
        },
        topic: 'einzel',
      },
      {
        title: 'Moving the focus',
        body: [
          'Select the lens and change its voltage. The focus moves — more voltage, shorter focal length — and the beam energy at the exit does not move at all.',
          'Notice how quickly it responds. A voltage change is not a re-solve: each electrode was solved once at one volt, and any other voltage is that solution multiplied. Changing a dimension is a re-solve, and feels different.',
        ],
        action: {
          label: 'Select the lens',
          run: (api) => api.select(1),
        },
        topic: 'fastAdjust',
      },
      {
        title: 'Give it something harder',
        body: [
          'A perfectly parallel beam flatters a lens. Real beams diverge, and a lens that looked sharp will show its aberrations once the outer rays come in at an angle.',
          'This sets a 2° half-angle: each ion is aimed outward in proportion to its distance from the axis. Watch the edge rays cross at a different place from the middle ones — that is spherical aberration, and no amount of tuning removes it.',
        ],
        action: {
          label: 'Add divergence and fly',
          run: (api) => {
            api.beam({ divergence: 2 });
            api.fly();
          },
        },
        topic: 'divergence',
      },
    ],
  },

  /* ---------------------------------------------------------------- */
  {
    id: 'seeing',
    title: 'Seeing the field',
    summary: 'Three ways of drawing the same solved field, each answering a different question.',
    steps: [
      {
        title: 'Potential, equipotentials, field lines',
        body: [
          'The colour map shows where the high and low ground is. Equipotentials are the contour lines of that map — they bunch where the field is strong. Field lines run along E, at right angles to the equipotentials.',
          'Turn all three on together once, to see how they relate, then leave on whichever answers the question you have.',
        ],
        action: {
          label: 'Show everything',
          run: (api) => {
            api.view({ field: true, contours: true, lines: true });
            api.select(1);
          },
        },
        topic: 'fieldLines',
      },
      {
        title: 'Field lines are not trajectories',
        body: [
          'This is the mistake worth not making. A field line is the path a particle would take if it had no inertia and always moved along the local force.',
          'An ion has plenty of inertia. It crosses field lines constantly, and in a deflector it does something a field line never does — it trades forward velocity for sideways velocity and turns a corner.',
        ],
      },
    ],
  },

  /* ---------------------------------------------------------------- */
  {
    id: 'bend',
    title: 'Turning a corner',
    summary: 'The quadrupole deflector, the matched voltage, and why it needs tuning.',
    steps: [
      {
        title: 'Not a sector bender',
        body: [
          'A sector pushes the beam round a corner with a field perpendicular to the orbit, and the path is a circular arc. A quadrupole deflector does something else: it sets up a quadrupole field in the plane the beam travels in, and lets the two axes trade velocity. The path is not an arc.',
          'Diagonalise it and half the motion oscillates while the other half grows exponentially. That growing half is why a deflector is touchy — an off-energy ion does not lag the design orbit, it leaves it — and why the device doubles as an energy filter.',
        ],
        action: {
          label: 'Build a deflector and fly',
          run: (api) => {
            api.build([
              DRIFT(20, 6),
              { type: 'bender', params: {} },
              DRIFT(40, 6),
            ]);
            api.beam({ mass: 100, charge: 1, energy: 50, rays: 9, beamRadius: 1.0, divergence: 0 });
            api.fly();
          },
        },
        topic: 'bender',
      },
      {
        title: 'The closed form is a starting point',
        body: [
          'The deflector arrives at its matched voltage — the value derived from requiring the ion to reach the exit with both the right position and the right direction.',
          'Measured against the solved field, that number runs 20 to 40 % low at every set of proportions. It is not an error in the algebra: the derivation knows only the aperture, while a real ion is also kicked by the entrance and exit channels on the way through, and those kicks push the same way as the bend.',
          'It always transmits. It just does not turn a right angle.',
        ],
      },
      {
        title: 'Let the tuner find the rest',
        body: [
          'Tuning searches the solved field for the voltage that actually delivers the beam, which is the honest answer to a closed form that is known to be approximate.',
          'This is cheap because a deflector voltage is only a multiplier on an already-solved field — the search flies ions many times but solves Laplace once.',
        ],
        action: {
          label: 'Tune the column',
          run: async (api) => {
            await api.optimise();
          },
        },
        topic: 'optimise',
      },
    ],
  },

  /* ---------------------------------------------------------------- */
  {
    id: 'branch',
    title: 'Branching the beam',
    summary: 'A deflector has three ways out, so a column is a tree rather than a list.',
    steps: [
      {
        title: 'Three ports, one voltage',
        body: [
          'The deflector’s box has an opening on all four faces. One is the entrance; the other three are all reachable, and the voltage alone decides which the beam takes — positive bends one way, zero goes straight through, negative bends the other.',
          'That is what these devices are for: sending a beam down one of three lines without moving any hardware.',
        ],
        action: {
          label: 'Build a switch',
          run: (api) => {
            api.build([
              DRIFT(20, 6),
              { type: 'bender', params: { voltage: 0 } },
              DRIFT(40, 6),
            ]);
            api.beam({ mass: 100, charge: 1, energy: 50, rays: 9, beamRadius: 1.0 });
            api.fly();
          },
        },
        topic: 'benderPorts',
      },
      {
        title: 'The chart is where branches live',
        body: [
          'A branching column is not a list, so the Beamline chart draws it as a tree. Each exit gets its own lane, and an exit with nothing on it is an empty socket rather than an error.',
          'Click a socket to start a new line there. Drag a box to move that element, and everything below it, onto a different exit.',
        ],
        topic: 'flowChart',
      },
    ],
  },

  /* ---------------------------------------------------------------- */
  {
    id: 'charge',
    title: 'Ions that see each other',
    summary: 'Space charge, and where it shows up first.',
    steps: [
      {
        title: 'Why the beam spreads in a drift',
        body: [
          'With repulsion off, ions fly as if alone. Turn it on and a beam in a field-free drift does the one thing it obviously should: it spreads, because the ions push each other apart with nothing to hold them in.',
          'A drift is where this is easiest to see, precisely because nothing else is happening there.',
        ],
        action: {
          label: 'Build a long drift, repulsion on',
          run: (api) => {
            api.build([DRIFT(150, 12)]);
            api.beam({ mass: 100, charge: 1, energy: 10, rays: 9, beamRadius: 1.0, divergence: 0 });
            api.physics({ repulsion: 'beam', beamCurrent: 20 });
            api.fly();
          },
        },
        topic: 'repulsion',
      },
      {
        title: 'Which model to believe',
        body: [
          'Beam treats the ions as a continuous current, which is what a steady DC beam is. Its spreading matches the closed-form envelope equation to about 0.1 %.',
          'Coulomb treats each flown ion as a clump of real ones and sums the forces between them — the right model for a short pulse, and much slower, because the cost grows as the square of the ion count.',
          'Lower the energy or raise the current and the effect grows quickly: a slow beam spends longer pushing itself apart.',
        ],
      },
    ],
  },

  /* ---------------------------------------------------------------- */
  {
    id: 'fringe',
    title: 'Fringe fields and shielding',
    summary: 'What the grounded end faces are hiding, and how metal really stops a field.',
    steps: [
      {
        title: 'The boundary is not neutral',
        body: [
          'Laplace’s equation needs a closed boundary, so each element is solved inside grounded end faces. But a grounded plate a few millimetres from the hardware is exactly what stops a fringe field — so solving an element alone does not merely omit its fringe, it shields it.',
          'Turn fringe fields on and neighbouring elements share one grid, so the field flows between them instead of stopping at an artefact.',
        ],
        action: {
          label: 'Build an aperture, fringe on',
          run: (api) => {
            api.build([
              DRIFT(25, 6),
              { type: 'aperture', params: { voltage: -300 } },
              DRIFT(40, 6),
            ]);
            api.beam({ mass: 100, charge: 1, energy: 50, rays: 9, beamRadius: 1.0 });
            api.physics({ fringe: true });
            api.view({ field: true, lines: true });
            api.fly();
          },
        },
        topic: 'fringe',
      },
      {
        title: 'How far a field really reaches',
        body: [
          'Inside a grounded pipe of radius R a disturbance dies as exp(−2.405·z/R). The decay length is R/2.405 and depends on nothing else — not the bore that made the field, not its voltage.',
          'So a 5 mm bore gives a decay length of 2.1 mm and kills a fringe within 15 mm. A 16 mm housing gives 6.7 mm and barely dents it over the same distance. The answer to keeping a field off something is an aperture at ground, not more distance.',
        ],
        topic: 'shielding',
      },
    ],
  },

  /* ---------------------------------------------------------------- */
  {
    id: 'keep',
    title: 'Catching and keeping',
    summary: 'A detector that pulls ions in, and saving the column you built.',
    steps: [
      {
        title: 'A detector is an electrode',
        body: [
          'At a few kilovolts below ground, a detector does not wait to be hit — its field reaches through the entrance aperture and pulls ions in, which is what collection efficiency means.',
          'It also distinguishes what it catches from what it loses: landing on the active face is a count, landing on the housing is not. Both are strikes, and one number for both would hide the difference the element exists to draw.',
        ],
        action: {
          label: 'Build a column with a detector',
          run: (api) => {
            api.build([
              DRIFT(20, 6),
              { type: 'einzel', params: { voltage: -300, boreRadius: 6, housingRadius: 15 } },
              DRIFT(30, 6),
              { type: 'detector', params: { voltage: -3000 } },
            ]);
            api.beam({ mass: 100, charge: 1, energy: 50, rays: 9, beamRadius: 1.0 });
            api.fly();
          },
        },
        topic: 'detector',
      },
      {
        title: 'Keep what you built',
        body: [
          'Save writes the column to a JSON file — the description, not the solution. Element types, their settings, how they hang together, the ion, the physics switches. No fields and no trajectories.',
          'Loading re-solves, which takes a moment and means an old file gets today’s physics rather than yesterday’s cached numbers.',
        ],
        topic: 'saveLoad',
      },
    ],
  },
];

/** A lesson by id, or null. */
export function lesson(id) {
  return LESSONS.find((l) => l.id === id) ?? null;
}

/** Every step, flattened, with its lesson — for tests and for progress. */
export function allSteps() {
  return LESSONS.flatMap((l, li) =>
    l.steps.map((s, si) => ({ lesson: l, lessonIndex: li, step: s, stepIndex: si }))
  );
}
