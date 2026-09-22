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

**Not included:** magnetic forces, buffer-gas collisions, space charge,
RF/time-dependent fields, relativistic correction. Each omission and its
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
| Integrators | analytic simple harmonic motion; orders 4 and 2 confirmed |
| Einzel lens | no net work, mirror symmetry, positive spherical aberration, reflection above the barrier |

The lens tests check properties the *real device* has, so they fail for
physical reasons rather than because an output number moved.

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
