# IonTrace

Open source ion trajectory simulation that runs in a browser. Solves Laplace's
equation on a potential array, then integrates Newton's second law through the
resulting field.

Built for physical chemists and ion opticians who want to see where their ions
go without installing anything.

**Status: v0.1.** Electrostatic optics in vacuum, one geometry (einzel lens).
The physics that is present is validated; the physics that is absent is
[documented](docs/PHYSICS.md#5-what-is-deliberately-absent).

---

## What it computes

$$\nabla^2\phi = 0 \quad\text{between electrodes}, \qquad
m\ddot{\mathbf{r}} = q\mathbf{E}, \qquad \mathbf{E} = -\nabla\phi$$

Geometry is axisymmetric, so the solve is two-dimensional in $(z, r)$ and
rotation supplies the third dimension exactly — not as an approximation.

**Space charge is included.** The beam repels itself, driven by a beam
*current* rather than by how many rays you happen to draw. Note that a ray
here is not a point charge — under rotational symmetry it is a **ring** of
charge at radius $r$, so the force follows from Gauss's law on the enclosed
current, not from pairwise Coulomb:

$$E_r(r) = \frac{\lambda_{\text{enc}}(r)}{2\pi\varepsilon_0 r},
\qquad \lambda = \frac{I}{v_z}$$

Because the force on each ion depends on where the others are at that instant,
the whole beam is integrated in lockstep on a shared time step.

**Not included:** magnetic forces, buffer-gas collisions, RF/time-dependent
fields, image charges, relativistic correction. Each omission and its
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

Read [docs/PHYSICS.md §7](docs/PHYSICS.md) before trusting a number. The
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
  geometries/einzel.js  the first element
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

MIT. See [LICENSE](LICENSE).
