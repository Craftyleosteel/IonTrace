# IonTrace

Open source ion trajectory simulation that runs in a browser. Solves Laplace's
equation on a potential array, then integrates Newton's second law through the
resulting field.

Built for physical chemists and ion opticians who want to see where their ions
go without installing anything.

**Status: v0.3.** Electrostatic and RF optics in vacuum, three-dimensional
trajectories, optional space charge, and beamlines assembled from drifts,
aperture plates, einzel lenses, quadrupoles and quadrupole deflectors — so a
column can turn a corner rather than only run straight. Voltages can be tuned
for transmission at the press of a button. The physics that is present is
validated; the physics that is absent is
[documented](docs/PHYSICS.md#5-what-is-deliberately-absent).

---

## Building a beamline

The main view is an editable column. Add elements, reorder them, remove them,
and fly ions through the result:

| Element | What it is |
|---|---|
| **Drift** | Field-free tube. Separates active elements so their fringe fields do not overlap. |
| **Aperture plate** | A charged plate with a hole. Not a passive opening — the equipotentials bulge through it, so it acts as a lens. |
| **Einzel lens** | Three coaxial cylinders, outer two grounded. Does no net work on a transmitted ion, and focuses for either polarity. |
| **Quadrupole** | Four rods, RF and DC. Converges in one transverse plane and diverges in the other at every instant; the RF is what makes it net-focusing in both. |
| **Quadrupole deflector** | Four curved electrodes in a grounded box, at $\pm V$ on the diagonals, turning the beam 90°. Not a sector: the field is a quadrupole *in* the bend plane, so the two axes couple and the path is not an arc. Rolls to any angle, so the column can turn left, down, or anywhere between. |

Elements carry their own local field solve, so **only the element you change
re-solves**, and voltages never re-solve at all.

That last point is what makes **Optimise voltages** practical. Since a voltage
is only a multiplier on a solved field, searching a few hundred combinations
costs a few hundred beam flights and no solver time — so the tuner can hunt for
the settings that transmit the most beam while you watch. Tune one element from
its own panel, or the whole column from the toolbar.

Elements arrive **already set for the ion in the source** — a lens at $-6\,T/q$,
a deflector at its matched voltage, a filter at a Mathieu $q$ of 0.38 — so
placing one and pressing Fly does something. All three follow the beam, and two
of them flip polarity for a negative ion; see
[PHYSICS.md §11](docs/PHYSICS.md). Every value is **typed rather than dragged**,
because a deflector's voltage spans four orders of magnitude across the beams
this handles and no slider serves both ends of that.

## What it computes

$$\nabla^2\phi = 0 \quad\text{between electrodes}, \qquad
m\ddot{\mathbf{r}} = q\mathbf{E}, \qquad \mathbf{E} = -\nabla\phi$$

**Ions move in three dimensions.** That is forced by the quadrupole: its field
depends on the azimuth, so an ion in it does not stay in a plane containing the
axis.

**Fields are still solved in two dimensions**, but which two depends on the
element. Einzel lenses, aperture plates and drifts are axisymmetric, so
$\partial\phi/\partial\theta = 0$ exactly and rotation supplies the third
dimension. A quadrupole's rods are uniform along $z$, so it is solved on a
*transverse* $(x, y)$ grid instead. A deflector's electrodes are uniform
perpendicular to the bend, so it is solved in the bend plane itself — the one
case where the field and the trajectory share a plane, which is exactly why its
axes couple. None of these is an approximation of a 3D solve; each is a
symmetry being used.

**Space charge is included.** The beam repels itself, driven by a beam
*current* rather than by how many rays you happen to draw. Note that a ray
here is not a point charge — under rotational symmetry it is a **ring** of
charge at radius $r$, so the force follows from Gauss's law on the enclosed
current, not from pairwise Coulomb:

$$E_r(r) = \frac{\lambda_{\text{enc}}(r)}{2\pi\varepsilon_0 r},
\qquad \lambda = \frac{I}{v_z}$$

Because the force on each ion depends on where the others are at that instant,
the whole beam is integrated in lockstep on a shared time step.

**RF fields are included.** A quadrupole's rod pairs are driven at
$W(t) = U + V\cos(\Omega t + \varphi)$, and because the pairs are always
driven antisymmetrically the field is exactly linear in $W$ - one solved map,
scaled by a coefficient that varies in time, with no re-solve. Mathieu $a$ and
$q$ are reported next to the controls that set them.

**Fringe fields are optional, and the toggle is a physics setting.** By default
each element is solved alone behind grounded end faces — which does not merely
omit its fringe field, it *shields* it, with a plate that is a numerical device
rather than hardware. Switch fringe fields on and stretches of axisymmetric
elements share one grid: fields reach into their neighbours, and a grounded
plate stops them.

That last part is why it needs a shared grid rather than a sum. A grounded
electrode contributes nothing to a superposition — it is at zero volts — yet it
changes a neighbouring lens's field completely, because it changes the
*boundary* of that lens's problem. Superposition is exact over voltages on
fixed conductors; it is not valid over geometry.

Inside a grounded pipe of radius $R$ a fringe dies as $e^{-2.405\,z/R}$, so a
narrow drift tube is already an excellent shield and a wide housing is a poor
one — measured against the solver to within 2–3 %. An einzel barely notices
(its outer cylinders are grounded, so it is nearly its own Faraday cage: 0.00 %
change); a charged aperture plate notices a lot (2.6 % on axis, and 3 of 9 ions
transmitted becomes 9 of 9). See [PHYSICS.md §10](docs/PHYSICS.md).

**Not included:** magnetic forces, buffer-gas collisions, image charges,
quadrupole fringe fields, relativistic correction. Each omission and its
consequences are listed in [docs/PHYSICS.md](docs/PHYSICS.md).

## Architecture

The design follows SIMION's, because that architecture is what makes electrode
geometry composable:

1. **Paint electrodes onto a grid.** `grid.paint(id, (z, r) => boolean)` marks
   nodes as fixed-potential boundary. Any geometry expressible as a predicate
   over $(z, r)$ is a valid element.
2. **Solve once per electrode.** Laplace's equation is solved with each
   electrode at 1 V and all others at 0 V, giving one unit solution each.
3. **Superpose to get any voltage set.** Because Laplace's equation is linear
   and all boundaries are Dirichlet,
   $\phi = \sum_i V_i \phi_i$. Changing a voltage costs one pass over the grid
   instead of a fresh relaxation — SIMION calls this *fast adjust*, and it is
   why the voltage sliders in the UI are instant while the geometry sliders
   pause to re-solve.

Adding a new element means describing its metal, not writing new field code.

## Running it

No build step, no dependencies. The site is plain static files.

ES modules will not load over `file://`, so it needs to be served over HTTP:

```sh
node serve.js          # then open http://localhost:8080
```

Or use any static server — the VS Code **Live Server** extension, `python3 -m
http.server`, or GitHub Pages.

### Tests

```sh
node tests/run-node.js
```

or open `tests/index.html` in a browser. Both run the same modules the site
loads, so there is no chance of the tested code and the shipped code diverging.

> If you have VS Code but no Node on `PATH`, its bundled Electron will run as
> Node:
> `ELECTRON_RUN_AS_NODE=1 "<path to>/Code.exe" tests/run-node.js`

## Validation

Every layer is checked against something outside itself, never against its own
output:

| Layer | Checked against |
|---|---|
| Constants | CODATA; textbook speed of a 1 eV electron |
| Laplace stencils | closed-form harmonic functions, including on the singular axis |
| Convergence | second order confirmed empirically against $1/\sqrt{r^2+z^2}$ |
| Fast adjust | superposition vs an independent direct solve |
| Conductor surfaces | full surface field recovered, not half |
| Absolute scale | mm→m pinned to literal metres; painted geometry matches the spec |
| Integrators | analytic simple harmonic motion; orders 4 and 2 confirmed |
| Space charge | closed-form uniform-beam field; cylindrical shell theorem; zero at zero current |
| Einzel lens | no net work, mirror symmetry, positive spherical aberration, reflection reported as reflection, focus independent of ion mass |
| Quadrupole | closed-form $(x^2-y^2)/r_0^2$ potential; four-fold symmetry; Mathieu $a$, $q$ and their scalings; stable ion transmits and unstable one is lost |
| Quadrupole deflector | solved potential is bilinear $\phi \propto XZ$; matched voltage derived from $\cot s = \tanh s$ and confirmed by integrating the coupled equations; measured transmission window |
| Voltage tuner | recovers a working deflector voltage from zero, and lands within 15 % of the independently derived matched value |
| Starting values | every element type transmits when placed, for ions from 4 u at 10 eV to 1000 u at 2 keV and for both polarities; `amplitudeForQ` inverts the Mathieu relation |
| Folded columns | a column bent through two right angles transmits; no branch of it claims another's ions; the beam is measured transversely to the axis it is actually on |
| Fringe fields | decay inside a grounded pipe matches $e^{-2.405z/R}$; a grounded plate measurably shields a charged one; a column solve agrees with the isolated solve where it should (an einzel, to 0.00 %) |
| Beamline | coordinate translation into placed elements; step size taken from the most demanding element; live chunking cannot change a trajectory |

The lens tests check properties the *real device* has, so they fail for
physical reasons rather than because an output number moved. The sharpest of
them is mass-independence: electrostatic optics depends only on $E/q$, so
1 u and 10 000 u ions at the same energy must focus at the same point (only the
flight time differs, as $\sqrt{m}$). A stray factor of $u$, of $e$, or of 1000
anywhere in the unit chain breaks it immediately.

This suite was written against an adversarial review. Four independent agents
audited the solver, the integrator, the units and ion-optics realism, and the
test suite itself — the last by mutation testing, deliberately breaking the
physics to find which tests failed to notice. `docs/PHYSICS.md` records what
they found, including the parts that are still wrong.

## Known limitations

Read [docs/PHYSICS.md §12](docs/PHYSICS.md) before trusting a number. The
headline ones:

- **Transmission is systematically pessimistic.** Electrode strikes resolve to
  the nearest grid node, shrinking every aperture by $h/2$ — a 6 mm bore models
  as 5.75 mm at the default resolution.
- **The lens converges at ≈ $O(h^{1.5})$, not $O(h^2)$**, limited by the field
  singularity at the sharp 90° electrode rims, where $|E| \sim \rho^{-1/3}$ and
  the peak field does not converge at all.
- **The energy-drift readout is a grid-quality number, not an integrator one.**
  It is flat in the time step and gauge-dependent.
- **Biasing the entrance or exit electrode** puts a spurious field across the
  drift regions, because the domain end faces act as grounded plates. The UI
  flags it.
- **Quadrupole fringe fields are absent.** The rods have a hard edge, so
  transmission through one is optimistic — entrance fringe loss is a real and
  well-known effect in mass filters.
- **Elements are solved in isolation**, so the field where two live elements
  meet is not a true solution for the pair. Leave a drift between them; the
  readout warns when you do not.
- **The tuner optimises what it is given.** It maximises transmission of the
  beam the source is currently set to produce, and a column tuned for nine ions
  at 1.5 mm is not necessarily tuned for a wider or more divergent one. It is
  also a local search: coordinate descent finds the best setting reachable from
  where it starts, which for a column with several interacting elements is not
  guaranteed to be the global optimum.
- **The tuner will not touch an RF quadrupole**, on purpose. Maximum
  transmission through a mass filter means switching the filter off, so its
  voltages are not the tuner's business — see
  [PHYSICS.md §9.5](docs/PHYSICS.md).

## Layout

```
index.html              the simulator
css/style.css
src/
  constants.js          SI constants and the only unit conversions
  grid.js               PotentialArray: the paintable domain
  laplace.js            SOR solver, planar + cylindrical stencils
  field.js              fast adjust, E = -grad phi, interpolation
  integrator.js         RK4 and velocity Verlet, adaptive step
  ion.js                practical units -> SI; beam factories
  spacecharge.js        ring/Gauss and discrete Coulomb models
  frames.js             rigid placements: where each element sits
  beamline.js           the column: layout, lookup, bounds, alignment
  optimize.js           voltage search for maximum transmission
  elements/             drift, aperture, einzel, quadrupole, deflector
    index.js            the registry the UI is generated from
  main.js               UI wiring and canvas rendering
tests/                  physics validation suite
docs/PHYSICS.md         equations, assumptions, and every omission
```

## Contributing

The priority is physics correctness over features. A pull request that adds an
element or a force term should come with tests that would fail if the physics
were wrong — a conservation law, a symmetry, or a closed-form limit — rather
than tests that pin current output.

If you find a physics error, please open an issue even if you are not sure. A
wrong trajectory that looks plausible is the failure mode this project most
needs help catching.

## Licence

IonTrace was created by Paco Navarro and is released under the MIT licence. See
[LICENSE](LICENSE).

You are free to use, copy, modify and redistribute it for any purpose,
including reshaping it to fit your own experimental setup and publishing the
result. The one condition is the MIT one: keep the copyright notice naming Paco
Navarro with any copy or substantial portion of the source you pass on.

If IonTrace contributes to something you publish - a paper, preprint, thesis,
poster or technical report - please cite it and link back to this repository:

> Navarro, P. (2026). IonTrace: open source ion trajectory simulation in the
> browser. https://github.com/Craftyleosteel/IonTrace

Citation metadata lives in [CITATION.cff](CITATION.cff), so GitHub's "Cite this
repository" button will generate BibTeX or APA for you. The citation is a
request rather than a licence condition; the attribution notice in the source
is the part the licence actually requires.
