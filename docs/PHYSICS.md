# The physics in IonTrace

This document states exactly what IonTrace computes, what it assumes, and what
it leaves out. It is meant to be adversarial towards the code: if the code and
this document disagree, one of them is a bug.

Version 0.1 models **electrostatic ion optics in vacuum**. Nothing else.

---

## 1. Governing equations

### 1.1 The field

Between electrodes there is no free charge, so Gauss's law

$$\nabla \cdot \mathbf{E} = \rho/\varepsilon_0$$

reduces with $\rho = 0$ and $\mathbf{E} = -\nabla\phi$ to Laplace's equation:

$$\nabla^2 \phi = 0$$

Each electrode surface is a Dirichlet boundary at its applied potential. A
conductor in electrostatic equilibrium is an equipotential, which is why a
fixed value on every electrode node is the correct boundary condition and not
an approximation.

### 1.2 The motion

$$m\frac{d^2\mathbf{r}}{dt^2} = q\,\mathbf{E}(\mathbf{r})$$

That is the complete right-hand side. There is no magnetic term, no drag, no
ion–ion interaction and no relativistic correction. Section 5 lists what each
omission costs.

### 1.3 Conserved quantities

With a static field the total energy

$$E = \tfrac{1}{2}mv^2 + q\phi(\mathbf{r})$$

is a constant of the motion. IonTrace tracks it on every flight and reports the
worst deviation, because a single number that must not change is the most
informative diagnostic available: a unit slip, a mis-scaled field, a bad
interpolation and an unstable time step all show up in it.

For an axisymmetric field, angular momentum about the axis
$L_z = m r^2 \dot\theta$ is also conserved. IonTrace currently launches ions
with $\dot\theta = 0$, so $L_z = 0$ and the trajectory stays in a plane
containing the axis for all time. This is what makes a 2D integration exact
rather than approximate — see §3.3.

---

## 2. Geometry and symmetry

The solved domain is two-dimensional in $(z, r)$ with rotational symmetry about
$r = 0$.

**This is not a 2D approximation of a 3D problem.** For a rotationally
symmetric electrode set, $\partial\phi/\partial\theta = 0$ exactly, so the 3D
Laplace equation collapses to a 2D one with no information lost. The third
dimension is recovered by rotation. A genuinely 3D geometry — a quadrupole rod
set, a slotted electrode, anything with azimuthal structure — cannot be
represented this way and is out of scope for this version.

A `planar` mode also exists, treating the transverse coordinate as a Cartesian
$y$ with no symmetry assumption. It is used by the test suite for cases with
closed-form solutions and is not exposed in the UI.

### 2.1 The enclosure

An unbounded Dirichlet problem has no unique solution on a finite grid, so the
domain is closed by a conducting box. This is a real modelling decision with
real consequences: the box is at a fixed potential and therefore shields, so
placing it too close to the electrodes distorts the field. `buildEinzelLens`
warns when the enclosure sits within one bore radius of the cylinders.

The two **end faces** are a special case. They must carry a fixed potential for
the solve to be well posed, but physically they are where the beam enters and
leaves. They are marked as open faces (`grid.openFaces`): electrically they are
boundary nodes, but an ion reaching one is reported as having left the modelled
region rather than as having struck an electrode. The outer radial wall is
*not* open — in the einzel geometry that surface is the grounded housing, which
is real metal an ion can genuinely hit.

---

## 3. Numerical methods

### 3.1 Discretising Laplace's equation

Uniform grid of step $h$. Planar geometry gives the standard five-point stencil

$$\phi_{i,j} = \tfrac{1}{4}\left(\phi_{i+1,j} + \phi_{i-1,j} + \phi_{i,j+1} + \phi_{i,j-1}\right)$$

Cylindrical geometry must also discretise the $\frac{1}{r}\frac{\partial\phi}{\partial r}$
term. Off the axis, at $r = jh$ with $j > 0$:

$$\phi_{i,j} = \frac{1}{4}\left[\phi_{i,j+1}\left(1 + \frac{1}{2j}\right) + \phi_{i,j-1}\left(1 - \frac{1}{2j}\right) + \phi_{i+1,j} + \phi_{i-1,j}\right]$$

**On the axis the term is singular** and must be replaced by its limit.
Symmetry forces $\partial\phi/\partial r = 0$ at $r = 0$, so by L'Hôpital's rule
$\frac{1}{r}\frac{\partial\phi}{\partial r} \to \frac{\partial^2\phi}{\partial r^2}$
and Laplace's equation becomes

$$2\frac{\partial^2\phi}{\partial r^2} + \frac{\partial^2\phi}{\partial z^2} = 0$$

Mirroring the ghost node $\phi_{i,-1} = \phi_{i,1}$ gives

$$\phi_{i,0} = \frac{1}{6}\left(4\phi_{i,1} + \phi_{i+1,0} + \phi_{i-1,0}\right)$$

A mis-derived axis stencil is one of the easiest errors to make here and one of
the hardest to see, because the solution still looks plausible. The test suite
checks it directly against $\phi = z^2 - r^2/2$, which is exactly harmonic in
axisymmetric coordinates and which the stencil reproduces to solver tolerance.

Both stencils are second-order accurate. `converges at second order in the grid
step` verifies this empirically against $1/\sqrt{r^2+z^2}$, a smooth harmonic
function the stencil *cannot* reproduce exactly, so a genuine discretisation
error appears and can be measured.

### 3.2 Relaxation and fast adjust

The linear system is solved by successive over-relaxation with

$$\omega_{\text{opt}} = \frac{2}{1 + \sqrt{1 - \rho^2}}, \qquad
\rho = \frac{\cos(\pi/n_z) + \cos(\pi/n_r)}{2}$$

Because Laplace's equation is linear and all boundary conditions are Dirichlet,
the solution for arbitrary applied voltages is a superposition of per-electrode
unit solutions:

$$\phi(\mathbf{x}) = \sum_i V_i\, \phi_i(\mathbf{x})$$

where $\phi_i$ solves the problem with electrode $i$ at 1 V and all others at
0 V. This is SIMION's "fast adjust". The relaxation is paid once per electrode;
changing a voltage afterwards costs one multiply–add per node.

**The superposition is only valid because every basis solution satisfies the
same homogeneous outer boundary condition.** The grounded enclosure sits at 0 V
in all of them. If one basis solution had a different outer boundary, summing
them would be wrong. `superposes basis solutions to the same answer as a direct
solve` checks the claim numerically rather than trusting the argument.

### 3.3 Field interpolation

$\mathbf{E} = -\nabla\phi$ is evaluated by central differences **at the grid
nodes**, and those nodal field components are then bilinearly interpolated to
the ion's position.

The order matters. Interpolating $\phi$ and differentiating the interpolant is
cheaper but produces an $\mathbf{E}$ that jumps discontinuously across cell
boundaries, injecting spurious impulses that an integrator accumulates into
drifting energy. Differencing first and interpolating second gives a continuous
$\mathbf{E}$.

On the axis, $E_r$ is set to exactly zero rather than differenced. This is not
a numerical convenience: rotational symmetry leaves a radial field at $r = 0$
with no direction to point in, so $E_r(0, z) = 0$ is exact, and imposing it
prevents round-off from deflecting an on-axis ion.

**Known inconsistency.** The force uses interpolated nodal $\mathbf{E}$, while
the energy diagnostic uses interpolated $\phi$. These two are consistent only to
$O(h^2)$, so the reported energy drift contains a contribution that is *not*
integrator error and does not shrink when the time step shrinks. It does shrink
when the grid is refined, and the test `reduces energy drift as the grid is
refined` exists specifically to confirm that this is the explanation. An
interpolation scheme that is exactly conservative — bicubic $\phi$ with its
analytic gradient — would remove it, at the cost of more code. It has not been
done for this version.

### 3.4 Trajectory integration

Two integrators, kept deliberately:

| Method | Order | Symplectic | Notes |
|---|---|---|---|
| Runge–Kutta 4 | 4 | no | Default. Highest accuracy per step; the method SIMION uses. Energy error grows slowly and secularly. |
| Velocity Verlet | 2 | yes | Energy error oscillates about zero rather than accumulating. Kept for future periodic-field (trapping) work. |

Running a case through both and comparing is evidence about the trajectory
rather than about either method's internal consistency, so the suite does that.

**Measuring their order requires care.** Integrating an oscillator for a whole
number of periods and comparing position lands on a turning point where
$dx/dt = 0$; a phase error $\delta$ then shifts $x$ by only $A\delta^2/2$, and
the measurement reports the wrong thing entirely — RK4's amplitude decay (5th
order) and the square of Verlet's phase error (4th order). The suite samples at
a generic phase instead, where phase error enters linearly, and recovers 4 and
2 as expected.

### 3.5 Time step

Adaptive, with two limits and the smaller winning:

- **travel**: $|v|\,\Delta t \le C h$, so the ion cannot skip over grid cells
  and miss the field structure between them;
- **acceleration**: $\tfrac{1}{2}|a|\,\Delta t^2 \le C h$, which takes over near
  rest, where the travel limit alone would permit an unbounded step.

$C$ (the "step fraction") defaults to 0.05.

### 3.6 Flight termination

A flight ends when the ion leaves the domain (`exited`), strikes an electrode
that is not on an open face (`electrode`), or exhausts its step or time budget.
Electrode strikes are resolved to the nearest grid node, so the reported impact
point carries an uncertainty of order $h/2$.

---

## 4. Units

SI throughout the physics: metres, seconds, kilograms, coulombs, volts.
Conversion happens only at the boundary, in `constants.js` and `makeIon`.
Nothing downstream of `makeIon` sees a dalton, an electronvolt, a millimetre or
a degree.

| Constant | Value | Source |
|---|---|---|
| $e$ | 1.602176634 × 10⁻¹⁹ C | exact, SI 2019 |
| $u$ | 1.66053906892 × 10⁻²⁷ kg | CODATA 2022 |
| $m_e$ | 9.1093837139 × 10⁻³¹ kg | CODATA 2022 |
| $\varepsilon_0$ | 8.8541878188 × 10⁻¹² F/m | CODATA 2022 |
| $c$ | 299792458 m/s | exact |

### 4.1 A note on the relativistic check

`relativisticError` reports the fractional error of treating a particle as
Newtonian. The obvious expression,

$$\frac{(\gamma - 1)mc^2 - \tfrac{1}{2}mv^2}{(\gamma-1)mc^2}$$

must not be evaluated as written. For an ion $\beta \sim 10^{-4}$, so
$\gamma - 1 \sim 10^{-8}$: forming it as $\gamma$ minus one discards eight of
the sixteen available digits, and the subtraction in the numerator cancels
most of what survives. The result is wrong in its first significant figure.

Writing $s = \sqrt{1-\beta^2}$, the expression reduces exactly to

$$\frac{\beta^2 (2 + s)}{2(1 + s)}$$

with no cancellation anywhere. It tends to $\tfrac{3}{4}\beta^2$ as
$\beta \to 0$ and stays exact to $\beta \to 1$. The test asserts against the
leading-order form, which is how the original error was caught.

---

## 5. What is deliberately absent

Each of these is a real effect that IonTrace does not model. Results are only
trustworthy where the corresponding term is negligible.

| Omitted | Equation term | When it matters |
|---|---|---|
| Magnetic force | $q\,\mathbf{v}\times\mathbf{B}$ | Any magnetic sector, ICR cell, or fringe field from a nearby magnet. |
| Buffer-gas collisions | drag + stochastic kicks | Any trap or guide with He/N₂ at > ~10⁻⁴ mbar. Dominates ion motion in collisional cooling. |
| Space charge | ion–ion Coulomb | Dense clouds and high beam currents. Causes emittance growth and, in traps, frequency shifts. |
| RF / time-dependent fields | $\phi(\mathbf{r}, t)$ | Every Paul trap, quadrupole filter and ion funnel. This build solves a static field only. |
| Image charge | induced surface charge | Very close electrode approach; small for typical bore radii. |
| Relativistic correction | $\gamma$ | Electrons above ~10 keV. Reported by `relativisticError`, warned on above 0.1 %. |
| Surface effects | patch potentials, roughness | Real instruments at high precision. |
| Secondary emission | — | Ion–surface impact; flights simply end at metal. |

The absence of time-dependent fields is the largest gap for the ion-trapping
audience. The architecture anticipates it: the fast-adjust superposition
$\phi = \sum_i V_i \phi_i$ already makes an RF field a matter of making $V_i$ a
function of time, with no re-solve. The integrator would need its time step
tied to the RF period rather than to the grid step, and the symplectic Verlet
option exists for exactly that reason.

---

## 6. Validation

`tests/physics.test.js` — run in a browser at `tests/index.html`, or under Node.
Every test is written to fail for a physical reason, not to lock in current
output.

| Layer | Checked against |
|---|---|
| Constants | CODATA values; textbook speed of a 1 eV electron (5.93 × 10⁵ m/s) |
| Planar stencil | $\phi = z^2 - y^2$, exactly harmonic |
| Cylindrical stencil | $\phi = z^2 - r^2/2$, exactly harmonic, including the axis |
| Convergence | $1/\sqrt{r^2+z^2}$; second order confirmed empirically |
| Solver correctness | Laplace residual; maximum principle (no interior extrema) |
| Fast adjust | superposition vs an independent direct solve |
| Axis symmetry | $E_r(0,z) = 0$ to machine zero |
| Integrators | analytic simple harmonic motion; orders 4 and 2 confirmed; exactness in a uniform field |
| Symplecticity | Verlet energy bounded over 300 oscillation periods |
| Einzel lens | no net work; mirror symmetry; focusing; positive spherical aberration; monotonic focal length vs bias; reflection above the barrier; RK4 vs Verlet agreement |

The einzel-lens tests are the ones that matter most to a user, because they
check properties the *real device* has. An einzel lens does no net work on a
transmitted ion, its trajectories are mirror symmetric, and its outer rays
focus before its inner ones. Those hold regardless of how good the numerics
are, so a violation is a bug no matter how plausible the picture looks.

---

## 7. Known limitations of this version

1. Electrode surfaces are staircase approximations on the grid. SIMION
   mitigates this with surface-enhanced refinement; IonTrace does not. Expect
   field error near sharply curved electrodes to be first order in $h$ there,
   even though the interior solution is second order.
2. The energy diagnostic is inconsistent with the force at $O(h^2)$ — §3.3.
3. Electrode strikes are resolved only to the nearest node.
4. Space charge, collisions and time-dependent fields are absent — §5.
5. Only one geometry (`einzel`) ships. The potential-array architecture is
   geometry-agnostic; `paint()` accepts any predicate over $(z, r)$, so adding
   an element means describing its metal, not writing new field code.
