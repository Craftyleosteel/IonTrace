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
$L_z = m r^2 \dot\theta$ is conserved, because there is no $E_\theta$ to
produce a torque. In IonTrace $L_z$ is not merely conserved but **identically
zero**, and structurally so: `IonState` has no azimuthal slot and `makeIon`
sets only $v_x$ and $v_z$, so azimuthal velocity cannot be introduced. The
orbit therefore stays in a plane containing the axis for all time and there is
no centrifugal term. This is what makes the 2D integration exact rather than
approximate.

Note the consequence: a real beam that *does* carry angular momentum — from a
magnetic lens upstream, a skewed source, or off-axis extraction — cannot be
represented at all, and nothing warns. Absence of $E_\theta$ conserves $L_z$;
it does not make it zero.

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

The two **end faces** are a special case. They carry a fixed potential, and an
ion reaching one is reported as having left the modelled region rather than as
having struck an electrode (`grid.openFaces`). The outer radial wall is *not*
open — in the einzel geometry that surface is the grounded housing, which is
real metal an ion can genuinely hit.

Two honest caveats about that choice:

**Dirichlet end faces are a convenience, not a necessity.** A mixed
Dirichlet/Neumann problem is uniquely solvable so long as Dirichlet data exists
somewhere on the boundary, which it does (the housing and the three cylinders).
Homogeneous Neumann $\partial\phi/\partial z = 0$ on the end faces would be
equally well posed and would be the better model of a drift tube continuing
past the domain. IonTrace uses Dirichlet because it is simpler, not because it
is required.

**The end faces behave as solid grounded plates across the aperture.** This is
harmless while the outer cylinders are also at 0 V — 15 mm of drift is enough
for the entry field to vanish, verified by extending it to 120 mm and seeing no
change. It stops being harmless if the entrance or exit electrode is biased:
with `entrance = 100 V`, the 0 V plate at $z = 0$ creates a real potential
difference across the entry drift and an accelerating field of several kV/m
where the actual instrument has none. Meanwhile the trajectory code flies ions
straight through that same face as though it were open. The field and the
flight then disagree about what is there. The UI flags this when the outer
electrodes are biased away from zero; the proper fix is Neumann end faces.

There is also a solver invariant worth stating: the relaxation never updates
rim nodes, so any rim node not owned by an electrode would sit frozen at 0 V —
an invisible grounded box that no voltage setting could override.
`solveBasis` refuses to run in that case rather than returning a plausible
wrong field.

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

Both stencils are second-order accurate **for a smooth solution on a smooth
boundary**. `converges at second order in the grid step` verifies this
empirically against $1/\sqrt{r^2+z^2}$, a smooth harmonic function the stencil
*cannot* reproduce exactly, so a genuine discretisation error appears and can
be measured.

That qualifier matters, and §8 gives the measured rate for the actual lens,
which is closer to 1.5 than to 2. The limiter is not the stencil.

### 3.2 Relaxation and fast adjust

The linear system is solved by successive over-relaxation with

$$\omega_{\text{opt}} = \frac{2}{1 + \sqrt{1 - \rho^2}}, \qquad
\rho = \frac{\cos\!\big(\pi/(n_z-1)\big) + \cos\!\big(\pi/(n_r-1)\big)}{2}$$

The denominators are node *spacings*, not node counts. This is an estimate in
any case: it is derived for an empty rectangle and ignores interior
electrodes, and on the shipped geometry it costs roughly 30 % more sweeps than
an empirically optimal $\omega$. That is a speed matter only — the converged
solution is identical.

Because Laplace's equation is linear and all boundary conditions are Dirichlet,
the solution for arbitrary applied voltages is a superposition of per-electrode
unit solutions:

$$\phi(\mathbf{x}) = \sum_i V_i\, \phi_i(\mathbf{x})$$

where $\phi_i$ solves the problem with electrode $i$ at 1 V and all others at
0 V. This is SIMION's "fast adjust". The relaxation is paid once per electrode;
changing a voltage afterwards costs one multiply–add per node.

**What makes the superposition valid is not that anything is grounded.** In
IonTrace the enclosure is itself a registered electrode, so its own basis
solution holds the enclosure at 1 V and nothing in that solve is at zero — and
superposition still holds exactly. The real conditions are:

1. the set of Dirichlet-frozen nodes is **identical** in every basis solve, so
   the geometry must not change between them; and
2. those frozen nodes are **partitioned** by the electrodes — every frozen node
   is owned by exactly one electrode.

If a frozen node were owned by no electrode, `solveBasis` would hold it at zero
in every basis solution and no combination of voltages could ever reproduce its
intended potential. Condition 2 is asserted at build time for that reason.

What would genuinely break superposition: a floating (charge-constrained)
electrode, space charge, a field-dependent permittivity, or any geometry change
between solves. `superposes basis solutions to the same answer as a direct
solve` checks the result numerically rather than trusting the argument.

### 3.3 Field interpolation

$\mathbf{E} = -\nabla\phi$ is evaluated by central differences **at the grid
nodes**, and those nodal field components are then bilinearly interpolated to
the ion's position.

The order matters, though **not** for the reason it is tempting to give.
Interpolating $\phi$ and differentiating the interpolant produces an
$\mathbf{E}$ that jumps across cell boundaries, and it is easy to assert that
those jumps pump energy. They do not: that scheme gives
$\mathbf{E} = -\nabla\phi_{\text{bilinear}}$ *exactly*, so the work along any
path is exactly $-\Delta\phi_{\text{bilinear}}$ and the energy is exactly
conserved. A discontinuous but conservative force cannot pump energy. Its real
defects are accuracy — $\mathbf{E}$ is only first-order and piecewise constant
normal to each edge — and a jumpy force that upsets adaptive stepping.

The chosen scheme trades that for the opposite property. Interpolating $E_z$
and $E_r$ independently gives a **continuous** field that is **not curl-free**:
$\partial E_r/\partial z \neq \partial E_z/\partial r$ within a cell, so no
potential exists whose gradient it is, and there is no exactly conserved
energy. Measured circulation around a sub-cell loop in the einzel fringe is
$6.8\times10^{-2}$ V against a 68.6 V drop across that cell — a part in
$10^{3}$. It is the better trade, but it should be understood as a trade.

On the axis, $E_r$ is set to exactly zero rather than differenced. This is not
a numerical convenience: rotational symmetry leaves a radial field at $r = 0$
with no direction to point in, so $E_r(0, z) = 0$ is exact, and imposing it
prevents round-off from deflecting an on-axis ion.

**Conductor surfaces need one-sided differences.** A central difference taken
*at* a node lying on an electrode reaches one node into the metal, where the
potential is pinned at the electrode value. It returns
$(\phi_{\text{vac}} - V)/2h$ where the true surface derivative is
$(\phi_{\text{vac}} - V)/h$ — exactly half. Being a factor rather than a
truncation term, it does not shrink with refinement (measured ratio
0.446 → 0.485 → 0.496 as $h$ falls), and bilinear interpolation then spreads
the halved value a full cell into the vacuum, leaving a fixed $-25\%$ error
half a cell off the metal. That is exactly where aperture-grazing rays fly, and
the error is one-signed, so it does not cancel along a trajectory. Such nodes
therefore use the one-sided difference on the vacuum side.

**Known inconsistency.** The force uses interpolated nodal $\mathbf{E}$, while
the energy diagnostic uses interpolated $\phi$. The *pointwise* mismatch between
them is $O(h)$ in the smooth interior, and $O(1)$ — non-convergent — next to
electrode rims, where it tracks the diverging corner field. What converges at
second order is the quantity that actually matters: the mismatch is oscillatory
with near-zero mean over a cell traversal, so it largely cancels and the
*path-accumulated* energy drift falls as $O(h^2)$. Measured:
$1.10\times10^{-2} \to 2.54\times10^{-3} \to 5.37\times10^{-4} \to
1.40\times10^{-4}$ for $h = 1, 0.5, 0.25, 0.125$ mm.

Two consequences for the reported energy drift:

- It is a **grid-quality** number, not an integrator-quality one. It is flat to
  within 2 % across a 40× change in time step, so it cannot detect a degraded
  integrator. `keeps energy drift independent of the time step` asserts that
  flatness deliberately, and `reduces energy drift at second order` asserts the
  $h^2$ scaling. Together they pin the explanation rather than the value.
- It is normalised by $|E_0| = |KE + q\phi|$, which is **gauge-dependent**:
  shifting every electrode by a constant leaves the field and the trajectory
  identical but changes the reported percentage. Treat it as an order of
  magnitude, not a figure of merit.

An interpolation scheme that is exactly conservative — bicubic $\phi$ with its
analytic gradient — would remove the inconsistency entirely, at the cost of more
code. It has not been done for this version.

### 3.4 Trajectory integration

Two integrators, kept deliberately:

| Method | Order | Symplectic | Notes |
|---|---|---|---|
| Runge–Kutta 4 | 4 | no | Default. Highest accuracy per step; the method SIMION uses. Energy error grows slowly and secularly. |
| Velocity Verlet | 2 | in practice | Energy error oscillates about zero rather than accumulating. Kept for future periodic-field (trapping) work. |

Verlet's symplecticity carries an asterisk. Formally it requires
$\mathbf{F} = -\nabla U$, and the interpolated field is not curl-free, so no
such $U$ exists and the map is not strictly symplectic. Empirically the
non-conservative part is oscillatory and cancels: over 3000 oscillation periods
on a real grid field, Verlet's relative energy error stays bounded in
$[-2.2\times10^{-2}, -5.6\times10^{-3}]$ with no trend while RK4's ramps
monotonically from $-0.11$ to $-0.49$. The behaviour the table claims is real;
the guarantee behind it is not exact, and it would be the first thing to break
on a coarse grid or a very long trapping run.

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

$C$ (the "step fraction") defaults to 0.05, exported as `DEFAULT_CFL` so the
code, this document and the UI cannot drift apart.

Both limits bound **displacement**, not time, and that has consequences worth
knowing. There is no upper bound on $\Delta t$ itself: an ion nearly at rest in
a nearly field-free region is given an enormous step — 6.6 ms for an ion
balanced at the top of a +1000 V barrier, against a device transit time of
microseconds. The trajectory is not wrong, because the displacement limit still
holds, but the reported flight time is meaningless and `maxTime` becomes a
post-hoc detector rather than a bound. Relative velocity change is also
unbounded: whenever the acceleration limit binds, the step permits the speed to
at least triple. Neither limit is an accuracy control — they prevent the ion
skipping field structure, nothing more. Step-size accuracy is not currently
estimated at all, and the energy drift figure cannot supply it (§3.3).

### 3.6 Flight termination

A flight ends in one of these states:

| `stop` | Meaning |
|---|---|
| `exited` | Left forwards through the far face. **Transmitted.** |
| `reflected` | Came back out of the entrance face. **Not transmitted.** |
| `electrode` | Struck metal, including the outer housing wall. |
| `time-limit` / `step-limit` | Ran out of budget. |

The distinction between `exited` and `reflected` is not cosmetic. An einzel
lens biased above the beam energy is a working ion mirror, and reporting a
reflected ion as transmitted lets its backwards axis crossing masquerade as a
focal length — a negative one, for a device that is not a lens at all.

**Electrode strikes are biased, not merely uncertain.** Resolution is to the
nearest grid node, which shrinks every aperture by exactly $h/2$ in one
direction. Measured on a field-free lens, the largest transmitted radius is
5.499 / 5.749 / 5.874 mm at $h$ = 1 / 0.5 / 0.25 mm against a true 6 mm bore —
always `bore − h/2`. At the default 0.5 mm step the modelled bore is 5.75 mm,
a 4 % radius deficit and an 8 % acceptance-area deficit. **Transmission figures
from IonTrace are therefore systematically pessimistic**, and the bias shrinks
only linearly in $h$.

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
| Space charge — **now modelled**, see §7 | ion–ion Coulomb | — |
| RF / time-dependent fields | $\phi(\mathbf{r}, t)$ | Every Paul trap, quadrupole filter and ion funnel. This build solves a static field only. |
| Image charge | induced surface charge | Very close electrode approach; small for typical bore radii. |
| Relativistic correction | $\gamma$ | Reported by `relativisticError`; warned on above 0.1 %, which for an electron is **341 eV** ($\beta = 0.037$). A 10 keV electron is already 2.9 % off. Ions are safe: a 1 keV, 100 u ion is off by $1.6\times10^{-8}$. |
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
| Solver convergence | Laplace residual (iteration only — see below); maximum principle |
| Fast adjust | superposition vs an independent direct solve; unowned rim refused |
| Conductor surfaces | full surface field recovered, not half; zero field inside metal |
| Absolute scale | mm→m pinned to literal metres; a known field in V/m; painted geometry matches the specification |
| Integrators | analytic simple harmonic motion; orders 4 and 2 confirmed; exactness in a uniform field |
| Symplecticity | Verlet energy bounded over 300 oscillation periods |
| Focus detection | first forward crossing; exact axis landing; no focus for on-axis, diverging or reflected rays; extrapolation beyond the grid |
| Einzel lens | no net work; mirror symmetry; positive spherical aberration; monotonic focal length vs bias; reflection reported as reflection; mass-independence of the focus; anion at mirrored polarity |
| Energy drift | $O(h^2)$ in the grid **and** flat in the time step |

Two cautions about what this table does *not* claim.

**`maxResidual` does not validate the stencil.** It evaluates the same
expressions the relaxation iterates, so a mis-derived stencil reproduces its own
error and the residual still falls to round-off — demonstrated with wrong radial
weights of $1 \pm 1/j$, which gives 9.4×10⁻² error against the analytic solution
while the residual reads 3.4×10⁻¹⁴. Likewise the maximum principle holds *by
construction* for any convex-weighted stencil, right or wrong. The stencils are
validated instead by the closed-form harmonic tests, which do reject that error.

**NaN defeats comparison-based checks.** `NaN > x` is false, so an accumulator
written `if (err > worst) worst = err` silently discards non-finite values and
reports zero error on a solution that has blown up entirely. Both the solver's
convergence test and the suite's error accumulator use `Math.max`, which
propagates NaN, and `assertFinite` guards the arrays directly.

The einzel-lens tests are the ones that matter most to a user, because they
check properties the *real device* has. An einzel lens does no net work on a
transmitted ion, its trajectories are mirror symmetric, its outer rays focus
before its inner ones, and its focus is independent of ion mass at fixed energy
because electrostatic optics depends only on $E/q$. Those hold regardless of how
good the numerics are, so a violation is a bug no matter how plausible the
picture looks.

The mass-independence test is also the sharpest check on the unit chain in the
suite: a stray factor of $u$, of $e$, or of 1000 anywhere between `makeIon` and
the force would break it immediately, and it is external to the code in a way
that ratio and ordering tests are not.

---

## 7. Space charge

The beam repels itself. This is the one term from §5 that **is** modelled.

### 7.1 Why it is not pairwise Coulomb

A trajectory drawn in the meridional plane is not one ion. Under rotational
symmetry it is the cross-section of a **ring** of charge at radius $r$,
carrying its share of the beam current all the way round the azimuth.
Computing $q_1q_2/4\pi\varepsilon_0 d^2$ between two such rays would be the
force between two point charges, which is not the force between two rings, and
would break the symmetry the whole field solve rests on.

The correct treatment for a long axisymmetric beam is Gauss's law. On a
cylinder of radius $r$ and length $L$ about the axis, $\mathbf{E}$ is purely
radial on the curved surface and the flat ends contribute nothing, so

$$E_r(r)\,2\pi r L = \frac{\lambda_{\text{enc}}(r)\,L}{\varepsilon_0}
\qquad\Longrightarrow\qquad
E_r(r) = \frac{\lambda_{\text{enc}}(r)}{2\pi\varepsilon_0 r}$$

Only charge **inside** $r$ matters. A uniform shell outside contributes exactly
nothing — the cylindrical shell theorem, which the test suite checks directly.

For a uniform beam of radius $R$ this gives a field rising linearly from the
axis, $E_r = \lambda r / 2\pi\varepsilon_0 R^2$, and that closed form is what
the implementation is validated against.

### 7.2 Driven by current, not by ray count

$$\lambda = \frac{I}{v_z}$$

A slow beam is a dense one, which is why space charge bites hardest where an
optic decelerates the beam. The input is therefore a **beam current** — a real
instrument parameter — and not the number of rays, which is a display setting
with no physics in it. Doubling the drawn rays must not double the repulsion,
and does not.

Each ray carries a fixed share $w_i$ of the current, set at launch from the
assumed initial current density (uniform, the only profile offered) in
proportion to the annulus it represents. The enclosed fraction for ray $i$ is
the sum of the shares of every ray currently inside it plus **half its own** —
a ring exerts no net force on itself, and the half places the ray in the middle
of its annulus.

Rays are re-sorted by radius every step rather than assumed laminar, so the
model stays valid after the beam crosses over, where the ordering genuinely
changes. **Rays at equal radius are one ring**: a beam launched over signed
offsets produces each radius twice, and treating those as two rings would
double the beam's current.

### 7.3 Lockstep integration

With space charge the ions can no longer be flown one at a time. The force on
each depends on where all the others are *at that instant*, so `flyBeam`
advances the whole beam on a shared time step — the smallest any active ion
asks for. An ion that strikes metal or leaves the domain stops contributing,
which is correct: it is no longer part of the beam.

The self-field is **frozen across each step** rather than re-evaluated at each
Runge–Kutta stage, because the other ions have no defined position at an
intermediate stage. This is the standard particle-in-cell treatment and it
costs accuracy: the space-charge part of the motion is effectively second order
even though the electrode part remains fourth.

### 7.4 What the beam does

Two regimes, both real:

- **Weak space charge.** The crossover survives and moves downstream as
  current rises, because the repulsion opposes the lens.
- **Strong space charge.** There is **no point focus at all**. The self-field
  goes as $1/r$, so as the beam converges the repulsion diverges. The beam
  reaches a minimum radius — a waist — and expands again, and that waist grows
  with current. A tool that insisted on reporting a focal length here would be
  reporting something that does not exist, so the readout switches to the
  waist instead.

Note also that with space charge on, the einzel lens's **no-net-work property
no longer holds**, and that is correct rather than a numerical failure: the
beam's own field does real work on its ions as it expands, converting the
bunch's electrostatic energy into transverse kinetic energy.

### 7.5 What space charge still neglects

| Neglected | Consequence |
|---|---|
| The beam's own magnetic field | Moving charges attract magnetically; the self-force is reduced by $(1-\beta^2)$. For keV ions $\beta\sim10^{-4}$, so this is a part in $10^8$ — genuinely ignorable. |
| Image charges in the electrodes | Surrounding metal partially shields the space charge. IonTrace's defocusing is therefore an **overestimate** for a beam that fills the bore. |
| Longitudinal space charge | The long-beam approximation assumes the beam is much longer than it is wide — true for a continuous beam, false for a short bunch. **Bunches are not modelled.** |
| Non-uniform current density | Only a uniform profile is offered. A peaked profile concentrates current at small radius and defocuses more strongly. |
| Self-consistent Poisson | The beam field is added analytically to the solved electrode field rather than Poisson being re-solved with $\rho$. Exact for the free-space part; this is what drops the image charges above. |
| Virtual cathode formation | $\lambda = I/v_z$ diverges as the beam is brought to rest. The implementation caps it, and a capped result means "outside the model's range", not an answer. |

---

## 8. Known limitations of this version

**1. The lens does not converge at second order, and the reason is not
staircasing.** Measured Richardson orders for `buildEinzelLens` at
$h$ = 1 / 0.5 / 0.25 / 0.125 mm: on-axis $\phi$ near the gap ≈ 1.5, peak
on-axis $|E_z|$ ≈ 1.4. There is in fact **no staircase at all** — at the
default parameters every electrode edge lands exactly on a node
(`rAt(12) === mmToM(6)` is exactly true), so the boundary is represented
perfectly.

The limiter is the **sharp 90° conducting rim** at each cylinder end. The
vacuum wedge there subtends $3\pi/2$, so $\phi - \phi_{\text{edge}} \sim
\rho^{2/3}$ and $|E| \sim \rho^{-1/3}$: the field at the rim **diverges** as
the grid is refined rather than converging. Measured $|E|$ one node inside the
bore from the corner: 47.0 → 60.4 → 77.8 → 99.8 kV/m, a ratio of 1.285 per
halving against the predicted $2^{1/3} = 1.26$. That is the edge singularity,
confirmed.

This is a property of the *geometry*, not of the solver, and it would be
present with a perfectly body-fitted mesh. Real hardware has a finite edge
radius. Consequences: the global rate is ≈ $O(h^{1.5})$, the focal length
carries ≈ 0.2 % error at the default grid and ≈ 0.6 % at the 1 mm setting, and
any future field-strength or breakdown readout must not quote the peak $|E|$,
which is meaningless.

2. Nodal $\mathbf{E}$ is *not* generally second order even where $\phi$ is:
   central-differencing an $O(h^2)$-accurate $\phi$ amplifies the error by
   $1/h$. Measured order for on-axis $E_z$ against a smooth analytic harmonic
   is ≈ 1.6–1.75.
3. The energy diagnostic is inconsistent with the force — §3.3. It is a grid
   diagnostic, not an integrator one, and it is gauge-dependent.
4. Electrode strikes are resolved only to the nearest node, biasing every
   aperture inward by $h/2$ — §3.6.
5. The SOR stopping test is "largest nodal change per sweep", which
   underestimates the true iteration error by a factor that grows with grid
   size (≈ 8× at $n$ = 33, ≈ 74× at $n$ = 257). At the default tolerance of
   $10^{-9}$ on a 1 V basis this still leaves only ~$10^{-7}$ V of error, but
   it is not the guarantee the name suggests.
6. `optimalOmega` is derived for an empty rectangle and ignores interior
   electrodes; it costs roughly 30 % more sweeps than the empirical optimum on
   the shipped geometry. A speed matter only — the converged answer is the same.
7. Space charge, collisions and time-dependent fields are absent — §5.
8. Only one geometry (`einzel`) ships. The potential-array architecture is
   geometry-agnostic; `paint()` accepts any predicate over $(z, r)$, so adding
   an element means describing its metal, not writing new field code.
