# The physics in IonTrace

This document states exactly what IonTrace computes, what it assumes, and what
it leaves out. It is meant to be adversarial towards the code: if the code and
this document disagree, one of them is a bug.

Version 0.2 models **ion optics in vacuum**: electrostatic and radio-frequency
fields, three-dimensional trajectories, optional space charge, and beamlines
assembled from interchangeable elements.

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

$$m\frac{d^2\mathbf{r}}{dt^2}
= q\Big[\mathbf{E}_{\text{electrodes}}(\mathbf{r}, t)
      + \mathbf{E}_{\text{self}}(\mathbf{r})\Big]$$

in three Cartesian components. The electrode term may depend on time (§8); the
self term is the beam's own space charge and is optional (§7). There is still
no magnetic term, no drag and no relativistic correction — §5 lists what each
remaining omission costs.

### 1.3 Conserved quantities, and when they are not

In a **static** field with **no** ion–ion interaction, the total energy

$$E = \tfrac{1}{2}mv^2 + q\phi(\mathbf{r})$$

is a constant of the motion, and IonTrace tracks it as its single most
informative diagnostic: a unit slip, a mis-scaled field, a bad interpolation
and an unstable time step all show up in it.

Both qualifiers matter, and the readout states which case applies rather than
reporting a large number as though something were wrong:

| Condition | Is $\tfrac{1}{2}mv^2 + q\phi$ conserved? |
|---|---|
| Static field, no repulsion | **Yes** — a deviation is numerical error |
| Repulsion on | No — the tracked quantity omits the ions' mutual potential energy, which converts into kinetic energy as the beam expands (§7) |
| RF element present | No — a time-dependent field does work on the ion, which is how a quadrupole confines it at all (§8.7) |

Angular momentum about the axis, $L_z = m r^2\dot\theta$, is conserved in an
axisymmetric field because there is no $E_\theta$ to produce a torque. It is
**not** generally zero: an ion launched with azimuthal velocity keeps it. What
is true is that an ion launched with $L_z = 0$ in a purely axisymmetric column
stays in a plane containing the axis for ever, which is why a lens can be
understood in two dimensions. A quadrupole breaks the symmetry outright, so
$L_z$ is not conserved there at all.

---

## 2. Geometry and symmetry

**Ions move in three dimensions.** The state is $(x, y, z)$ with velocity
$(v_x, v_y, v_z)$, and that is forced by the quadrupole: its field depends on
the azimuth, so an ion in it does not stay in a plane containing the axis.

**Fields are still solved in two dimensions**, but which two depends on the
element:

| Element | Solved on | Why it is exact |
|---|---|---|
| Einzel lens, aperture plate, drift | $(z, r)$, axisymmetric | $\partial\phi/\partial\theta = 0$ exactly, so the 3D Laplace equation collapses to 2D with nothing lost. Rotation recovers the third dimension. |
| Quadrupole | $(x, y)$, transverse | The rods are uniform along $z$, so the potential does not depend on it. |

Neither is an approximation of a 3D solve; each is a symmetry being used.

A planar launch stays planar in an axisymmetric element — there is no
azimuthal field to take the ion out of the plane — which is why the
two-dimensional picture remains exact for a lens, and why moving to 3D changed
nothing for the cases that never needed it. The test suite checks this
explicitly.

### 2.0 Beamlines

Elements are composed into a column. Each is placed at its own $z$, and a field
query finds the element containing that point and delegates to it in local
coordinates. Adding a new kind of optic means describing its metal and its
field; the composition layer does not change.

**By default, elements are solved in isolation, not as one system.** Each is a
separate Dirichlet problem with grounded end faces, so its fringe field is
confined to its own footprint and stops abruptly at the boundary.

That default can be turned off. With fringe fields enabled, stretches of
axisymmetric elements are solved together on one grid and neighbouring
electrodes do see each other — including the case that matters most, a grounded
plate shielding a charged one. **§10** covers what changes, by how much, and
why superposition cannot substitute for it. The rest of this section describes
the isolated model, which remains the default and is what the accuracy figures
below refer to.

**How good is that?** Measured, by building the same column both ways and
comparing the on-axis potential:

| What changes | On-axis error |
|---|---|
| Splitting a lens into drift + lens + drift rather than one wider lens | 0.004 V out of 300 V |
| Element margin 0.8 bore radii | ≈ 19 % |
| Element margin 1.6 bore radii | ≈ 12 % |
| Element margin 2.4 bore radii | ≈ 6.6 % |
| Element margin 3.2 bore radii | ≈ 3.5 % |

The grounded faces *between* elements cost almost nothing. What costs accuracy
is each element's own margin clipping its **own** fringe field, and that error
decays roughly exponentially in margin / bore.

So the rule is not "leave a drift between elements" — it is **give each element
about three bore radii of its own margin**. Elements warn when they do not have
it. A drift between two live elements still helps, because a real grounded
drift tube genuinely shields, which is exactly what the isolation approximation
is pretending.

Solving the whole column at once in three dimensions is what a general-purpose
code does, and it would destroy the property that makes this interactive:
only the element you changed re-solves, and voltages never re-solve at all.
Solving an *axisymmetric run* at once is far cheaper — the symmetry is still
being used, so it is one more $r$–$z$ problem rather than a 3D one — and that
is what §10 does. Voltages stay fast adjust even there.

**A related trap, now fixed and worth recording.** Node coordinates are
accumulated as `z0 + i·h`, while geometry is written as `mm × 10⁻³`. Those two
routes to the same number differ in the last bit, and which way they differ
depends on where the element sits in the grid. Comparing them with a bare
`<=` therefore painted an electrode one node shorter or longer purely because
of its absolute position — so the *same* element solved to two different fields
in two different beamlines. On a 15 mm centre electrode at 0.4 mm resolution
that single node moved the on-axis potential by 7.8 V out of 300, because it
sits where the field changes at 20 V/mm. `grid.spans()` exists to make these
comparisons tolerant, and a test asserts that an element paints identically
wherever it is placed.

### 2.0.1 A column is a tree, and a deflector is a junction

A quadrupole deflector has **two ways out**. At its matched voltage it turns
the beam through a right angle; turned off, the beam goes straight through the
box and out the opposite face. That is not a modelling convenience — it is how
these devices are used, as a switch that sends a beam down one of two lines
without moving any hardware. The box has apertures on all four faces, which is
the four-fold symmetry of the electrodes rather than a special case.

So elements form a tree: each carries `from = {parent, port}`, and a deflector
can have hardware bolted to both of its exits. **Which branch the ions take is
decided by the field, not by the topology.** Measured on a switch with a line
on each exit: 9 of 9 out the straight branch at 0 V, 9 of 9 round the bend at
the matched 40 V, and nothing out of either at half or one and a half times it
— that last being the honest answer, since a partly-deflected beam lands on the
electrodes.

Three consequences worth stating:

**"Transmitted" needs a destination.** On a straight column the question does
not arise. On a branching one it is the whole question: a deflector with
hardware on both exits transmits everything at zero volts, straight out the
back. A search told merely to maximise transmission would discover that and
switch the deflector off — which is not tuning a beamline, it is unplugging it.
So the column names a **main end**, reached by always taking an element's first
exit, and the optimiser counts only ions that arrive there.

**Path distance no longer identifies a place.** Two elements on different
branches sit at the same distance from the source. `zStart` is still distance
along *that element's own branch*, which is what progress metrics need, but
anything that has to point at a location uses the frame, not the number.

**A run of shared field solves is a chain, not a slice.** Two branches leaving
one junction are adjacent in the element array and point in different
directions; painting both onto one $r$–$z$ grid would put one branch's hardware
on top of the other's. §10's runs follow the tree.

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
| RF / time-dependent fields — **now modelled**, see §8 | $\phi(\mathbf{r}, t)$ | — |
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

The beam repels itself. This is the one term from §5 that **is** modelled, and
there are **two models**, because they answer different questions.

| Model | A trajectory is… | Driven by | Use for |
|---|---|---|---|
| `coulomb` | one ion, feeling every other directly | ions per particle | a countable bunch or cloud |
| `beam` | a **ring** of charge | beam current | a continuous beam |
| `none` | non-interacting | — | lens characterisation |

Picking the wrong one is a physics error, not a preference. §7.1 explains why.

### 7.0 Discrete Coulomb

Each particle is a point charge:

$$\mathbf{E}_i = \frac{1}{4\pi\varepsilon_0}\sum_{j\neq i}
\frac{w\,q_j\,(\mathbf{r}_i-\mathbf{r}_j)}{|\mathbf{r}_i-\mathbf{r}_j|^3}$$

**Macro-weighting.** $w$ is how many real ions each simulated particle stands
for. At $w=1$ the calculation is literally $N$ ions, and for any $N$ you can
draw on screen the repulsion is far too small to see — nine elementary charges
spread over millimetres give a field of order $10^{-3}$ V/m against electrode
fields of $10^{4}$. That is the correct answer, not a defect, and it is why
particle codes weight. A macroparticle of weight $w$ has charge $wq$ and mass
$wm$, so $q/m$ and the electrode force are unchanged; only the mutual force
scales, linearly in $w$. The $w$ real ions inside one macroparticle do not
repel each other — the standard particle-in-cell approximation.

**Softening.** The $1/r^2$ force diverges as two particles approach, and with a
finite time step a close pass would fling them apart with energy from nowhere.
Plummer softening replaces $|d|^3$ with $(|d|^2+\epsilon^2)^{3/2}$, bounding the
force at short range while leaving the long range untouched. It is a real
approximation: close encounters, and therefore collisional relaxation of the
bunch, are suppressed.

Cost is $O(N^2)$ — trivial at the particle counts a browser will draw, and the
reason this does not scale to a real PIC simulation.

**The energy diagnostic does not apply.** IonTrace tracks
$\tfrac{1}{2}mv^2 + q\phi$ for the *electrode* field only. The ions' mutual
potential energy is not in it, and that energy is genuinely converted into
kinetic energy as the bunch expands, so a large reported "drift" with repulsion
on is the physics working rather than the integrator failing. A conserved total
for the interacting system — which would have to include
$\sum_{i<j} k\,w\,q_iq_j/d_{ij}$ and be evaluated for the system rather than
per ion — is not implemented. With repulsion on, the drift figure carries no
information about numerical quality.

### 7.1 Why the beam model is not pairwise Coulomb

For a **continuous** beam, a trajectory drawn in the meridional plane is not
one ion. Under rotational symmetry it is the cross-section of a **ring** of
charge at radius $r$, carrying its share of the beam current all the way round
the azimuth. Computing $q_1q_2/4\pi\varepsilon_0 d^2$ between two such rays
would be the force between two point charges, which is not the force between
two rings, and would break the symmetry the whole field solve rests on.

This is the distinction between the two models. Use `coulomb` when the
trajectories really are individual ions you could count; use `beam` when each
one stands for a continuous stream of them.

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

## 8. The quadrupole and RF fields

### 8.1 Why it is not axisymmetric

The ideal linear quadrupole potential is

$$\phi(x, y, t) = \big(U + V\cos\Omega t\big)\,\frac{x^2 - y^2}{r_0^2}$$

positive along $x$ and negative along $y$. The azimuth is the whole point, so
this cannot live on an axisymmetric grid at all. It is solved instead on a
**transverse** $(x, y)$ grid, with real round rods painted rather than ideal
hyperbolae — so the solution carries the higher multipoles a real rod set has.
The 12-pole in particular is why the rod-to-field radius ratio $1.1487$ is an
engineering choice and not a detail.

### 8.2 Every instant is defocusing in one plane

At any moment the field converges in one transverse plane and diverges in the
other. A **DC** quadrupole is therefore always unstable in one plane — it is a
singlet, and needs a partner of opposite sign to make a net-focusing doublet,
which is how accelerator FODO cells work. What makes an **RF** quadrupole
different is that the sign alternates faster than the ion can escape, so the
time-averaged force is inward in both planes for the right parameters.

### 8.3 Fast adjust with a time-dependent coefficient

The two rod pairs are always driven antisymmetrically, at $+W(t)$ and $-W(t)$,
so the field is exactly linear in $W$ and a single unit solution suffices:

$$\mathbf{E}(x, y, t) = W(t)\,\mathbf{E}_{\text{unit}}(x, y),
\qquad W(t) = U + V\cos(\Omega t + \varphi)$$

Nothing is approximated by that factorisation. It is the same superposition
used everywhere else, with a coefficient that happens to vary in time, and it
means the RF costs no re-solve at all.

### 8.4 Stability

$$a = \frac{8zeU}{m r_0^2 \Omega^2}, \qquad q = \frac{4zeV}{m r_0^2 \Omega^2}$$

These, not the voltages, decide whether an ion is transmitted. For an RF-only
filter ($a = 0$) the first stability region ends at $q = 0.90803$. Since
$q \propto 1/m$, each mass sits at a different $q$ at fixed RF, which is what
makes the device a mass filter.

**Stability is asymptotic.** It is a property of the Mathieu equation in the
limit of many oscillations, so an ion that crosses the rods in a handful of RF
cycles can be lost whatever its $(a, q)$ says. The shipped default gives about
thirty cycles; at eight, most of the beam is thrown out despite sitting
comfortably inside the stability region. This is real, not numerical, and is a
known constraint on short filters.

### 8.5 The step size must resolve the RF

A time-dependent field adds a third limit to §3.5: $\Delta t \le C\,T_{RF}$,
twenty steps per period at the default. Without it, a slow ion could take a
step spanning several field reversals and integrate a field it never
experienced — an unstable trajectory would look perfectly confined. The
Runge–Kutta stages are evaluated at $t$, $t+\Delta t/2$, $t+\Delta t/2$ and
$t+\Delta t$, which for a static field is redundant and for an RF one is the
difference between fourth order and second.

### 8.6 Hard edge, and what it costs

The field is uniform along the rods and zero outside them. Real fringe fields
at the rod ends are **not** modelled, and they matter: in a real mass filter
the entrance fringe is a well-known cause of transmission loss, because an ion
crosses a region where the RF is still ramping up. So transmission through this
element is **optimistic**, and an off-axis ion's potential energy jumps at the
rod ends because the potential is discontinuous there. The transverse dynamics
inside the rods, which is what sets stability, are unaffected.

Softening the edge with a longitudinal envelope $g(z)$ would be worse, not
better. $g(z)(x^2-y^2)$ does not satisfy Laplace's equation unless $g'' = 0$,
so it would trade a known approximation for a field that is not a field.

### 8.7 Stability is asymptotic, and a short filter is not stable

Mathieu stability says where the motion stays bounded **for ever**. It says
nothing about an ion that crosses the rods in a handful of RF periods, and such
an ion can be thrown out with a perfectly respectable $(a, q)$.

This is not a fine point. Measured: a 70 mm filter at 1.2 MHz gives a 100 u,
50 eV ion **8.6 cycles** and loses eight ions in nine, with nothing in the
stability numbers to show for it. The same filter at 150 mm and 2 MHz gives
30.5 cycles and transmits all nine. The cycle count,

$$N = \frac{L}{v}\,f = \frac{L f}{\sqrt{2T/m}},$$

is therefore shown in the interface beside $a$ and $q$, and flagged below
about 15. It depends on the ion's **energy**, which $a$ and $q$ do not — two
beams with identical Mathieu parameters can behave completely differently if
one crosses in five cycles and the other in fifty.

### 8.8 Energy is not conserved, and should not be

With a time-dependent field, $\tfrac{1}{2}mv^2 + q\phi$ is **not** a constant
of the motion: the RF does work on the ion, which is how a quadrupole confines
it at all. The energy-drift diagnostic of §1.3 therefore carries no information
about numerical quality once an RF element is in the line, and the readout says
so rather than reporting a large figure as though something were wrong.

---

## 8.9 Multipole guides and ion funnels

Both are RF devices, and neither is governed by Mathieu stability. What governs
them is the **effective potential** — Dehmelt's result for an ion in a field
oscillating faster than the ion can follow:

$$U^* = \frac{q^2 E_0^2}{4 m \Omega^2}$$

with $E_0$ the local field *amplitude*. IonTrace evaluates it from the solved
field rather than from an ideal form, so the depth quoted is for the device as
drawn.

### 8.9.1 Why more poles guide better

The ideal $2n$-pole potential is $\phi = V (r/r_0)^n \cos n\theta$, so
$|E| \propto r^{n-1}$ and therefore $U^* \propto r^{2n-2}$: **flat across the
middle and steep at the rods**. An ion spends most of its time where the RF is
weak, so it is heated less and the guide accepts a far wider mass range than a
quadrupole. The cost is a weak restoring force near the axis, which is why
multipoles make good traps and poor filters.

Recovered from the solved field of real round rods, fitting $|E| \propto r^k$
near the rim:

| Rods | $k$ measured | $n-1$ |
|---|---|---|
| 4 | 1.02 | 1 |
| 6 | 2.01 | 2 |
| 8 | 3.02 | 3 |
| 12 | 5.02 | 5 |

Fitted **near the rim** deliberately. For a 12-pole the multipole term at a
third of the field radius is a quarter of a per cent of the rim field, so
anything else present — and round rods always leave something — swamps it.
Measured at $r = 0.3r_0$ the same fit returns 3.97 rather than 5, which is a
statement about round rods, not about the solver.

Against the closed form for an ideal quadrupole, the solved well comes out
**9 % shallower**: round rods at the classic $1.1487 r_0$ cancel the 12-pole
term but do not reproduce a hyperbola, and the field at $r_0$ is correspondingly
weaker.

### 8.9.2 The funnel, and the gas that is not here

An ion funnel is a stack of rings with a shrinking aperture, adjacent rings in
opposite RF phase, and a DC gradient down the stack. It is solved in $r$–$z$
like a lens, with each ring its own electrode so that the DC ramp and the RF
alternation are two weighted readings of one basis:

$$\phi(r,z,t) = \phi_{\mathrm{DC}}(r,z) + \cos(\Omega t)\,\phi_{\mathrm{RF}}(r,z)$$

Both terms are fast adjusts; changing the gradient or the amplitude never
re-solves.

**A boundary condition that silently destroyed the beam.** With the domain's
end faces grounded, a DC ramp painted on the rings becomes a *hill*: zero at
the entrance, 13 V inside, zero again at the exit. A 5 eV beam was reflected in
its entirety — 0 of 9 through at every beam width — and nothing about that
looked like a boundary-condition error. Each end face now takes the potential
of the ring beside it, which is what a funnel embedded in a longer line
actually looks like, and the ramp falls monotonically from end to end.

**Without buffer gas, a deeper well transmits worse.** This is the honest core
of it. Measured on the shipped stack at 4 MHz:

| RF wall depth | Through |
|---|---|
| 0.3 eV | 9/9 |
| 1.1 eV | 9/9 |
| 4.3 eV | 3/9 |
| 17 eV | 0/9 |

The wall confines, but the RF heating that comes with it has nowhere to go, so
the ion climbs out. In a real funnel — 1 to 30 mbar — collisions carry that
energy away continuously, and that is not a detail of the device, it is the
mechanism. At the frequency a real funnel runs at, 0.5–1 MHz, a 100 u ion
crosses this stack in under six RF cycles, too few for an effective potential
to mean anything at all, and transmission is zero.

So this element reproduces a funnel's **field**, and the guiding a cold
near-axis beam gets from it. It does not reproduce a funnel's ability to
collect a warm, diffuse cloud, which is the job funnels exist to do. Its
defaults are chosen to work in vacuum and are therefore *not* the settings the
hardware is operated at. The interface says so rather than letting a
transmission figure stand unqualified.

---

## 9. The quadrupole deflector

### 9.1 A quadrupole *in* the bend plane

The bender in this simulator is not a cylindrical sector. It is the device a
real beamline uses to fold its path: four curved electrodes filling the corners
of a grounded box, at $+V$ and $-V$ on the diagonals, with the beam entering
through the gap on one axis and leaving through the gap on the perpendicular
one.

Each electrode is a **square block** filling its quadrant, with a circular arc
of radius $r_0$ cut out of the inner corner and a **straight-sided channel** of
width $2w$ between it and each neighbour. Four grounded posts stand on the
diagonals, in the corners of the box, and each block is cut back to clear them.

### 9.1.1 Only the arcs reach the beam

The aperture is bounded by the four electrode arcs and the four channel mouths,
and an electrode is a **conductor**: nothing behind its surface can influence
the field in front of it. That one fact settles most of the geometry, and it is
worth stating because it is counter-intuitive enough to have been got wrong
here once already.

- **Filling the corners changes nothing inside.** The arc is still at $\pm V$
  and the aperture cannot tell what is behind it. The blocks are that shape
  because that is the instrument.
- **The posts change nothing inside**, for the same reason — they stand behind
  the blocks. Measured at the shipped proportions, the potential 8 mm off axis
  on the diagonal is **37.28 V with posts and 37.28 V without**, identical to
  every digit the test prints. They are structure, not field shaping.
- **The channel width is what reaches the beam**, because it is the only choice
  here that changes the arcs themselves: a wider channel cuts them shorter, and
  a shorter arc is a weaker quadrupole.

### 9.1.2 Arc coverage sets the voltage

Each arc spans a half-angle $\alpha = \pi/4 - \psi_0$ about its diagonal, where
$\psi_0 = \arcsin(w/r_0)$ is the angle at which the channel wall cuts the
circle. Projecting the boundary onto $\sin 2\theta$,

$$\phi = \frac{4V}{\pi}\sin 2\alpha \left(\frac{r}{r_0}\right)^{2}\sin 2\theta
       = F\,V\,\frac{2XZ}{r_0^{2}}, \qquad F = \frac{4}{\pi}\cos 2\psi_0,$$

so $F$ would be the factor by which a real electrode beats — or misses — the
ideal normalisation $V_0$ assumes, with the voltage needed scaling as $1/F$.

**It over-predicts threefold, and the reason matters more than the model.**
Between a 10.5 mm channel and a 16 mm one, $F$ says the field should fall to
0.762. Measured on the solved field: **0.915**.

The model assumes the gaps sit at zero, as a gap in a driven boundary would. A
channel is not a gap. Its two walls belong to *neighbouring* blocks, which are
at opposite polarity, so the slot is lined with $+V$ facing $-V$ and goes on
driving field into the aperture long after the arc has ended. Widening the
channel moves those walls apart; it does not remove them. So the aperture field
is far less sensitive to channel width than arc coverage alone suggests, and
$F$ is useful only as a statement of direction.

### 9.1.3 The channel is a deflector in its own right

That $+V$ facing $-V$ across a slot is a **parallel-plate deflector**, and the
beam must cross one to get in. Treating the slot as uniform, of length
$L = $ `electrodeThickness` and half-gap $w$, the sideways throw at the matched
voltage is

$$d \;\approx\; \frac{k}{4}\left(\frac{r_0}{a}\right)^{2}\frac{L^{2}}{w},
\qquad k = 1.8556,$$

growing with the **square** of electrode depth and falling only as $1/w$. At the
shipped proportions with a 10.5 mm channel that is 3.1 mm of throw against
5.25 mm of clearance — and the measured transmission was 4 of 9.

This is real behaviour rather than a modelling artefact: a deflector with deep
channels really does steer its own beam into its entrance, which is why
instruments put a grounded aperture plate across the mouth. This model has no
such plate, so the default channel width is chosen against $d$ rather than
against $F$ — the field barely notices the difference, and the acceptance does.
A thin-electrode geometry has a short slot and never sees the effect at all,
which is why the calibration tests in §9.2 pin their own channel width.

The distinction is not cosmetic. A sector deflects with a field everywhere
perpendicular to the orbit, and the trajectory is a circular arc. A quadrupole
deflector sets up a two-dimensional quadrupole potential **in the plane the
beam travels in**,

$$\phi = C\,x\,z,$$

which couples the two axes rather than acting centrally:

$$m\ddot{x} = -qC z, \qquad m\ddot{z} = -qC x.$$

Axial velocity is traded for transverse velocity as the ion crosses, and the
path is **not** an arc. Diagonalising with $u = x+z$ and $w = x-z$ separates
them, at the cost of making plain what kind of device this is:

$$\ddot{u} = -b\,u \quad\text{(oscillatory)}, \qquad
  \ddot{w} = +b\,w \quad\text{(exponential)}, \qquad b = qC/m.$$

Half of the motion **grows** rather than oscillating. That is why a deflector
is touchy: an ion that is off-energy or off-axis does not merely lag the design
orbit, it departs from it exponentially. It is also what makes the device an
energy filter rather than merely a corner.

### 9.2 The matched voltage, and its evil twin

Requiring the ion to reach the exit aperture with both the right position and
the right direction closes the problem. With $s$ half the transit phase, the
oscillatory branch demands $v/(a\omega) = \cot s$ and the exponential branch
demands $v/(a\omega) = \tanh s$, so the design condition is

$$\cot s = \tanh s \;\Longrightarrow\; s = 0.937552,\quad
  \left(\tfrac{a\omega}{v}\right)^2 = \coth^2 s = 1.85565,$$

whence, with $T$ the kinetic energy, $r_0$ the aperture radius and $a$ the
half-width of the field region,

$$V_0 = 1.8556\,\frac{T}{q}\left(\frac{r_0}{a}\right)^2 .$$

That equation has further roots, the next near $s = 2.347$ giving $k = 0.9641$,
and direct integration of the coupled equations confirms it is *also* a clean
ninety degrees — the ion takes a longer way round, bending the other way. **It
is not usable.** On that branch the oscillatory amplitude is $1.427a$ rather
than $1.241a$, carrying the orbit out to $\rho = 1.01a$: past the electrode
faces at $r_0$, which is necessarily smaller than $a$. The ion lands on an
electrode instead of reaching the exit. The short branch peaks at $0.88a$ and
clears the aperture. Picking the wrong root is not a small error, and it does
not announce itself as one — the algebra is equally valid either way, and only
integrating the orbit tells them apart.

### 9.3 What the solved field does instead

The formula assumes the ideal quadrupole potential everywhere inside the box
and nothing outside it. The solved field is neither: electrodes subtend finite
arcs, the grounded box shapes the field near the apertures, grounded posts
interrupt the diagonals, and the field does not stop abruptly at the entrance
plane. Each of those moves the voltage that turns the beam through exactly 90°
away from $V_0$ by a few per cent.

Measured against the solved field, with a nine-ion beam of 1.5 mm radius at
100 u and 50 eV. An ion counts only if it leaves **by the bend port** — the
deflector has four, and with nothing attached to the straight one an ion that
sails through unbent also leaves the column, so counting exits alone would
score an *off* deflector as perfectly tuned:

| $r_0/a$ | turns 90° at | all nine bent over |
|---|---|---|
| 0.655 (shipped) | $V/V_0 = 1.22$ | 0.84 – 1.38 |
| 0.792 | 1.20 | 0.86 – 1.34 |
| 0.905 | 1.30 | 0.92 – 1.46 |
| 0.950 | 1.38 | 0.96 – 1.56 |

The closed form is systematically **low**, by 20–40 %, and never right. That is
not a defect in the algebra — it is §9.1.3 arriving: the ideal derivation knows
only the aperture, while a real ion is also kicked by the entrance and exit
channels on its way through. The window is wide enough (±25 % or better) that
$V_0$ always transmits, which is what makes it a usable starting point, but it
does not turn a right angle and the tuner is what finds the voltage that does.

The channel is the dominant control, at fixed $r_0/a = 0.655$:

| channel | turns 90° at | all nine bent over |
|---|---|---|
| 10.5 mm | $V/V_0 = 0.86$ | 0.70 – 0.96 |
| 13 mm | 1.02 | 0.76 – 1.14 |
| 15 mm (shipped) | 1.22 | 0.84 – 1.38 |
| 18 mm | 1.48 | 0.96 – 1.70 |

A *narrower* channel needs *less* voltage, and by a lot — far more than the
8.5 % the aperture field moves across this range. The reason is that the
channel kick points **the same way as the bend**: the $+V$ wall is the one on
the outside of the turn, so the slot starts steering the ion before it reaches
the aperture and the two effects add. Narrowing the channel strengthens that
kick as $1/w$ while barely touching the quadrupole, which is why $V_{90}$ moves
so much more than $F$ does.

What does not depend on that measurement: the closed form is the right *size*,
not the final answer, and the voltage tuner (§10) searches the solved field for
the value that actually maximises transmission. That is why the tuner exists.

### 9.4 The box has holes in it

The Laplace problem needs the domain closed, so the grounded box is painted all
the way round the solve. The beam apertures are real, though, and the collision
test carves them back out: a clear channel of half-width $w$ on each of the
four faces, the same $w$ that separates neighbouring electrodes. Because the
channels are straight-sided, that is one number rather than a width that
depends on where along the channel it is measured. Treating the closed rim as
metal destroys every ion on the entrance plane.

This is the same approximation as the open end faces of §2.1, and it has the
same justification: a grounded plane with a hole in it is very nearly a
grounded plane. Note that the grid's default open faces are the *wrong* ones
here — this element and the mass filter both solve in a plane of their own, so
for them the grid's $z$ ends are ordinary walls, and both close them
explicitly.

### 9.5 Rolling it into any plane

The bend plane is a roll about the beam, applied when the field is evaluated
rather than baked into the solve, so one solved deflector serves every
orientation and an arbitrary angle costs exactly what a quarter turn costs.
0° turns the beam left, 90° down, 180° right, 270° up, and anything between
bends into the plane at that azimuth — verified to nine figures against the
exit direction, with full transmission at 30°, 45° and 137°.

The roll is **conjugated**: the transverse axes are rotated into the bend
plane and back out again. Without that, a vertical deflector would also rotate
"up" into "sideways" for every element downstream, which is not what bending a
beamline does to the hardware bolted after it.

### 9.6 Two things a bend breaks that a straight column hides

**Element lookup.** Elements answer `contains` on their *axial* extent alone,
deliberately: an ion inside an element's length but outside its bore is still
that element's business, and should be reported as striking its wall rather
than as having wandered out of the column. In a straight line that is
unambiguous, because axial position identifies an element uniquely. It is not
once the path folds. A column bent through two right angles runs its last
drift back alongside its first, inside the first one's axial range but eighty
millimetres off its axis — and every ion entering the last drift was reported
as striking the wall of the first. The symptom is indistinguishable from
physics: the beam stops, at a plausible place, for a plausible reason.

So an element claims a point only if the point is in its **free space**.
Whichever element holds it in vacuum wins. Failing that it may still be
claimed as a strike, but only within 1.6 times its own outer radius — enough
to cover a step that overshoots a wall and the corners of a box whose diagonal
exceeds its stated radius, not enough to reach another branch of a folded
column. Beyond every element's envelope the point belongs to nothing, which is
what lets an ion at the exit plane be called an exit rather than a crash.

**The transverse plane.** A beam profile means the offsets perpendicular to
the direction of travel. Global $x$ and $y$ are those only while the beam runs
along $z$: after a right-hand bend the beam travels along $-x$, so global $x$
is now *along* the beam and global $z$ is transverse. A profile plotted in
global $(x, y)$ smears across the screen in proportion to distance flown. Each
ion is therefore measured in the frame of the element it is passing through,
and an ion past the end against the exit frame.

### 9.7 Tuning

Because every electrode voltage is a multiplier on a stored unit solution
(§3.2), a search over voltages costs beam flights and no solver time at all.
`src/optimize.js` exploits that: coordinate descent with a coarse-to-fine
bracketed scan on each knob, revisiting knobs because they interact.

**The objective.** Transmission is an integer count and therefore flat almost
everywhere: most voltages transmit nothing, and the search has no direction to
move in. Two much smaller terms break the flatness — partial credit for how far
a lost ion travelled down the column, and a tie-break on the size of the
surviving beam transverse to the **exit** axis. Both are capped below $1/N$, so
neither can outweigh a single transmitted ion; the reported figure is always
the honest count.

Gradient descent would be the wrong tool. The objective is genuinely
discontinuous — an ion either clears an aperture or does not — so a derivative
is meaningless and a scan is what finds the operating window.

**Where it looks.** A deflector's slider has to reach tens of kilovolts,
because real deflectors run there, while a 50 eV beam is matched at a few tens
of volts and its whole window is about a tenth wide. Any affordable sweep of
the full range steps straight over the answer. So where an element can say what
voltage it expects — the deflector can, from §9.2 — the sweep is centred on
that estimate instead, spanning $\pm 3 V_0$ so the opposite polarity a negative
ion needs is still reachable. The estimate itself, both polarities, and the
setting already in place are all tried explicitly, rather than left to whether
an evenly spaced scan happens to land on them.

Each refinement then brackets $\pm$ one sample spacing around the best point.
That width is the right one by construction: a scan at spacing $d$ can only
have missed the optimum by less than $d$. A fixed fraction of the range would
not close in at all on a knob whose range is three orders of magnitude wider
than its answer.

**A second stage, and what a Hessian can and cannot do here.** Transmission is
a count: flat almost everywhere, jumping by one when an ion clears an aperture
it previously hit. It has no derivative, let alone a second one, so no Newton
method of any kind applies to it — the scan is not a stopgap.

But once the beam is through, a different question starts: of the many settings
that transmit everything, which delivers the best beam? That *is* smooth. The
spot size at the target moves continuously with the voltages while no ion is
near an edge, so it has a gradient and a Hessian — and second-order information
is exactly what coordinate descent lacks. `src/reduced.js` finite-differences
both, eigendecomposes the Hessian by cyclic Jacobi, and takes a Newton step.

Three things make that work rather than merely sound plausible:

**Scaled coordinates.** Each knob is divided by a natural size taken from the
same closed form that sets its starting value — $-6T/q$ for a lens, the matched
voltage for a deflector, $q = 0.38$ for a filter. Without that the Hessian's
eigenvalues mix units and "near zero" has nothing to be near zero against. A
lens with no declared scale fell back to a sixth of its permitted range, 13 kV
for an answer of a few hundred volts, and the resulting Hessian was noise.

**Measuring at the target plane.** Each ion's crossing is interpolated onto the
exit plane rather than read off its last recorded point. The last point sits
wherever the adaptive step happened to leave off, which shifts as the voltages
change — so a finite difference over a small change reads that jitter rather
than the physics. Before this, the Hessian's largest eigenvalue came out
negative and equal in magnitude to its smallest: a pure numerical saddle.

**Treating three kinds of curvature differently.** A plain Newton step divides
by the curvature, so a near-zero eigenvalue produces an enormous jump along the
direction that matters least.

| Curvature | What it means | What is done |
|---|---|---|
| $\lambda > \tau$ | a real minimum in this direction | Newton step, $-g_k/\lambda$ |
| $\lambda < -\tau$ | downhill — a saddle | step with $-g_k/\lvert\lambda\rvert$, still a descent direction |
| $\lvert\lambda\rvert \le \tau$ | flat — the null space | left alone, and **reported** |

with $\tau$ a small fraction of the stiffest curvature present. Skipping
negative curvature instead of using $\lvert\lambda\rvert$ leaves the method
stalled at saddles, which is where an untuned beamline normally sits: measured,
a column 3 mm from its best took a 0.2 V step and stopped.

**The null space is the interesting output.** A beamline's knobs are not
independent — a lens and the deflector behind it trade off, so the merit surface
has a flat-bottomed valley and the Hessian a near-zero eigenvalue along it. That
direction is a combination of voltages the beam does not care about, which says
the tuned setting is not unique and which knob is redundant. It is reported
rather than discarded.

Measured on a lens-and-deflector column starting from zero volts: the scan alone
reaches an exit beam of 0.344 mm rms and the refinement takes it to 0.327 mm.
Modest, because the scan already gets close — and the larger gain from this work
came from something the investigation turned up rather than from the Newton step
itself. Coordinate descent could not find *any* working setting from zero,
because sweeping the lens while the deflector is off transmits nothing whatever
the lens does, so the sweep learns nothing and moves on. Every element already
knows roughly what it wants, so the combination of all those guesses is now
tried first, for the cost of one flight: 0 of 9 became 9 of 9.

**What makes it fast enough to watch.** Profiled rather than guessed at. One
flight of nine ions runs about 8,400 integrator steps, which produce 377,000
field evaluations and **603,000 element lookups** — four field queries per ion
per step plus a strike test, each of which was scanning the whole column and
calling both `contains` and `strikes` on every element it passed. Five changes,
measured on a lens-and-deflector search from zero volts (12.2 s originally):

| Change | Why it is free | After |
|---|---|---|
| Element lookup remembers the last element it found | Almost every query is about a point a fraction of a millimetre from the previous one, or another ion in the same beam. The hint is re-checked before it is trusted, and trajectories come out identical to the last bit. | 9.4 s |
| Settings already flown are not re-flown | Coordinate descent re-samples its own bracket centre at every level and re-walks ranges on a second pass. The flight is deterministic, so the answer is. Cut 141 evaluations to 96. | 6.7 s |
| The scan ranks at a coarser time step | Its job is to rank settings, not measure one, and ranking survives a much coarser integration. The winner is re-scored at full fidelity before anything is reported or refined. | 3.2 s |
| Points inside a clear aperture skip the strike test | Three quarters of strike tests are about a point comfortably inside a bore whose metal starts there, where one comparison settles it instead of a grid lookup. Elements declare a `clearBore`; a deflector cannot and does not. | 2.9 s |
| The scan's step raised from 4× to 8× | Ranking still exact on two columns including an RF one. | **2.4 s** |

The coarse step is a genuine trade-off, so it was measured before being taken.
Over a sweep of lens voltages on two columns — one ending in a deflector, one
containing an RF mass filter, which is the harder case — at 4×, 8× and 16× the
default step the **ordering of every setting is unchanged and the transmission
counts are identical**, while a flight drops from 70 ms to 17 and 8. What does
move, in the fifth decimal, is the tie-break between settings that all transmit
fully; those go on to the refinement stage, which works at full fidelity and is
a local search anyway.

### 10.7.1 Three things that looked faster and were not

Worth recording, because each is the obvious next idea:

**Velocity Verlet for the scan.** RK4 evaluates the field four times per step,
Verlet once, so it looks like a fourfold saving. It is not. Verlet's smaller
stability limit makes the adaptive stepper take more steps, so the real gain is
28 %; and at the coarse step it **reorders the sweep**, swapping two settings
whose scores differ in the fourth decimal. A cheaper integrator that changes
the answer is not a faster search.

**Fewer ions.** The same sweep reorders when the beam is thinned from nine to
five, because the tie-break is a mean over whichever ions are present and a
different sample is a different mean.

**Dropping the second coordinate pass.** The reasoning was that a second pass
exists to catch knob interactions and the Newton stage now handles those
properly through the off-diagonal Hessian terms. On a two-knob column that
holds — one pass and two give an identical 0.327 mm. On a three-knob column it
does not: two passes reach 0.847 mm against one pass's 0.936 mm. The pass stays.

**What it will not touch.** Only voltages: geometry needs a fresh solve, and
leaving it alone is the constraint an operator at a real instrument works under
anyway. And *not the RF quadrupole*, which is the one exclusion worth arguing
about. Its voltages do change transmission — but a mass filter's job is
selectivity, and the setting that transmits the most beam is the one that
filters nothing. Given the RF amplitude to play with, the search turns it down
until the rods stop selecting: measured on a lens–filter–deflector column, it
dropped 250 V to 40 V and counted it an improvement. That is the optimiser
working correctly on the wrong objective. Tuning a filter means tuning $(a, q)$
for a target mass, a different search with a different goal, and it does not
belong behind a button labelled "best transmission".

---

## 10. Fringe fields and column solves

### 10.1 An isolated element does not merely omit its fringe — it shields it

Every element is normally solved on its own grid, closed by grounded faces at
each end. Those faces are a numerical device: Laplace's equation needs a closed
boundary. But physically they are a **grounded plate a few millimetres from the
hardware**, and a grounded plate is precisely what stops a fringe field. The
isolated model therefore does not simply truncate the fringe; it terminates it,
with the wrong thing, at an arbitrary distance.

Turning fringe fields on solves stretches of neighbouring elements on one grid.
The internal faces vanish, and the field flows out of each element into
whatever is actually next to it.

### 10.2 Why superposition cannot do this

The cheap alternative — solve each element alone, then add the fields where
they overlap — cannot show shielding, and shielding is the point.

A grounded electrode contributes **nothing** to a superposition: it is at zero
volts, so its unit solution enters the sum multiplied by zero. Yet placing one
next to a charged lens changes that lens's field completely, because it changes
the **boundary** of the lens's own problem. Charge rearranges on both.

Superposition is exact over *electrode voltages on a fixed set of conductors* —
that is the fast adjust this whole simulator rests on (§3.2). It is not valid
over *geometry*. Adding a conductor is a different problem, not a larger sum.
So: one grid, all the metal on it, one solve.

### 10.3 What sets the decay: the pipe, not the source

Inside a grounded pipe of radius $R$, a disturbance at one end dies along the
axis as the lowest Bessel mode,

$$\phi(z) \sim e^{-j_{01} z / R}, \qquad j_{01} = 2.40483,$$

so the decay length is $R/2.405$ and depends on **nothing else** — not the bore
of the element that made the field, not its length, not its voltage. Measured
against this solver, fitted decay lengths past an einzel across housing radii
of 10–20 mm land within 3 % of $R/2.405$; measured downstream of a charged
aperture plate in drift tubes of 12 mm and 8 mm bore, within 2 %.

The practical consequences are worth stating plainly:

| Tube bore | Decay length | Field left after 15 mm |
|---|---|---|
| 3 mm | 1.2 mm | $4\times10^{-6}$ |
| 6 mm | 2.5 mm | $2\times10^{-3}$ |
| 12 mm | 5.0 mm | 0.05 |

A **narrow grounded drift is already an excellent shield**; a wide housing is a
poor one. This is why the answer to "how do I stop this field reaching my
detector" is usually an aperture plate at ground, not more distance.

### 10.4 How much it matters, element by element

Measured on axis at the middle of each element, isolated solve versus column
solve, with 12 mm drift tubes either side:

| Element | Difference | Why |
|---|---|---|
| Einzel lens | **0.00 %** | Its outer two cylinders are grounded, so it is very nearly its own Faraday cage. Only 0.285 V escapes past the exit cylinder of a lens sitting at 1946 V on axis — a part in seven thousand. |
| Aperture plate | **2.6 %** | A single charged plate with grounded space on both sides. Nothing shields it but distance. |

The trajectory consequence is larger than the potential difference suggests.
A 50 eV, nine-ion beam through a $-500$ V plate into a 12 mm tube: **3 of 9
transmitted with a 10.3 mm exit radius** when the plate is solved alone,
**9 of 9 with a 4.6 mm radius** when solved with its neighbours. Truncating the
field cuts the second half of the plate's lens action, so the ion leaves with a
radial kick that was never meant to stand on its own.

### 10.5 What a run is, and where it stops

Only axisymmetric elements can share an $r$–$z$ grid, so a run covers
consecutive drifts, aperture plates and einzel lenses and stops at a quadrupole
or a deflector. That is not much of a loss: both of those sit inside grounded
housings that genuinely terminate the field rather than merely appearing to.
A run also stops at a misaligned element, which is no longer a body of
revolution about the column axis.

Runs of nothing but drift are skipped: no electrode in them can ever hold a
voltage, so the solution is zero everywhere and the elements' own empty fields
say the same thing for free.

The run's own two end faces remain a numerical device — the problem still has
to be closed — but they now sit at the ends of the *run* rather than of every
element. How much potential is still alive there is measured after each solve
and warned about above 1 % of the largest electrode voltage.

Cost, for a browser: the shipped starting column re-solves in **34 ms**; a
seven-element lens stack, 606 × 36 nodes and 13 electrodes, in **226 ms**.
Voltages remain fast adjust over the shared grid, so the tuner is unaffected.

### 10.6 The quadrupole's hard edge, and a correction

§8.6 claims that softening the quadrupole's hard edge with a longitudinal
envelope "would trade a known approximation for a field that is not a field",
because $g(z)(x^2-y^2)$ does not satisfy Laplace's equation unless $g''=0$.
The premise is right and the conclusion is wrong. The naive product is not a
field, but it is the first term of one that is:

$$\phi = a(z)(x^2-y^2) - \frac{a''(z)}{12}(x^4-y^4) + O(a^{(4)}),$$

which satisfies $\nabla^2\phi = 0$ to the order shown — the quartic term exists
precisely to cancel the $a''(x^2-y^2)$ left over by the quadratic one. This is
the standard soft-edge multipole expansion of beam optics, and it is a
legitimate way to model a fringe field. It is not implemented here; the
quadrupole still has a hard edge, and its transmission is still optimistic.
The claim that it *could not* be done is withdrawn.

---

## 11. Starting values

An element placed from the toolbar arrives set for **the ion currently in the
source**, not for whichever ion the defaults happened to be written for. This
is a usability decision with a physics justification, and getting it wrong is
not a cosmetic failure: before it, three of the five element types did nothing
useful when placed. A deflector arrived at zero volts — so it did not deflect
at all, and the beam flew straight into the far wall of its box. An einzel
arrived forty times too strong for the default beam and transmitted three ions
in nine. A quadrupole arrived at $q = 1.83$, outside the first stability region
entirely, and transmitted none.

Each starting value comes from a closed form, so it follows the beam:

| Element | Set from | Why that quantity |
|---|---|---|
| Einzel lens | $V = -6\,T/q$ | Electrostatic optics depends only on $E/q$, so the lens is fixed by the centre potential as a multiple of $T/q$ and by nothing else — measured identical at 50 eV and 500 eV. Everything from $-1$ to $-10$ transmits fully; over-focusing sets in near $-20$. |
| Quadrupole deflector | $V_0 = 1.8556\,(T/q)(r_0/a)^2$ | §9.2. |
| RF quadrupole | $V$ giving $q = 0.38$ | Inverting $q = 4zeV/(mr_0^2\Omega^2)$. Not the apex of the stability region, which is where resolution is best and transmission most fragile — 0.38 is confining and well clear of the $q = 0.908$ cut-off. |

Two of these are **signed by the charge**, and that is the case a magnitude
gets wrong. A negative ion meeting the centre electrode of a lens set for a
positive one is decelerated rather than accelerated, cannot climb the barrier,
and is reflected; a deflector at the wrong polarity steers it into the wall
instead of round the corner. An RF filter has no polarity to get wrong, because
the drive reverses every half cycle regardless.

Nothing here ever moves a value the user has set. Change the beam energy
afterwards and every element stays exactly where it was, with the readout
saying how far off it now is — which for a deflector *is* its energy
selectivity, and worth seeing rather than hiding.

### 11.1 Why the controls are typed, not dragged

The same scaling argument decides the control. A deflector's voltage spans four
orders of magnitude across the beams this simulator handles: 8 V for a 10 eV
ion, 45 kV for a 27 keV one. A slider wide enough for the second puts the first
inside a single pixel — on a ±60 kV range, one pixel is over 500 V. The number
is typed instead, which is exact at every scale. The stepper arrows still nudge
by an amount chosen for the beam in front of the element, taken from the same
closed forms above.

---

## 12. Known limitations of this version

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
7. Collisions, magnetic forces and image charges are absent — §5. Quadrupole
   fringe fields are absent — §8.6.
8. Elements are solved in isolation rather than as one column, so the field
   where two live elements meet is not a true solution for the pair — §2.0.
9. Four element types ship: drift, aperture plate, einzel lens and quadrupole.
   The architecture is geometry-agnostic - `paint()` accepts any predicate over
   the solved plane - so adding an element means describing its metal and its
   field, not touching the composition layer.
