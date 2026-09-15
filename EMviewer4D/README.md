# Switch-on bench

The battery comes on over a few nanoseconds, and you watch the field find out.

Open `index.html` in a browser. No build step, no dependencies, no WebGL; works over
`file://`.

The two steady-state benches — [`../EMviewer`](../EMviewer/) in a plane and
[`../EMviewer3D`](../EMviewer3D/) in space — answer *what is the field once everything has
settled*. Everything in them is instantaneous: Coulomb and Biot–Savart with no retardation
anywhere, which is the right thing for a DC circuit and hides the question this one asks.
None of the three depends on the others.

## Why it exists

A circuit diagram implies the load learns about the battery the moment the switch closes.
It doesn't. Three separate delays show up here, on three different scales:

- **The field spreads at c from the battery.** Before the switch there is no field at all;
  then a front of it sweeps outward, and the colour map goes exactly black outside — not
  dim, black. The dashed sphere is that boundary, radius ct about the battery, and nothing
  crosses it.
- **The current lags the EMF.** The loop is an inductor, and the trace in the rail shows
  the two curves failing to coincide: the EMF arrives in a couple of nanoseconds, the
  current takes L/R to follow. The gap between the curves is what builds the field.
- **Energy arrives last.** In **E×B** mode nothing flows at a point until both fields
  exist there, so the load cannot be fed faster than light can carry the news.

Turn the wire resistance up and L/R shortens; open the loop out and L grows. The numbers
panel prints L/R against the time light takes to run along the wire, which is the ratio
that decides whether the model below is on firm ground.

## How it's computed

**Causality first.** A quasi-static solve on its own is *acausal*: it redistributes charge
along the whole wire the instant the battery moves, so every segment starts radiating at
once. Measured, that put field at a probe at up to 226% of its settled value 1.75 ns before
news from the battery could have got there, and charge on stretches of wire the light had
not reached. So each segment's response is delayed by its own light time from the battery,
`δᵢ = |cᵢ − battery| / c`. Nothing can then reach a probe sooner than
`(|battery→segment| + |segment→probe|)/c`, which is never less than `|battery→probe|/c`.

Two things must be repaired after delaying, and they are why the source is tabulated rather
than evaluated in closed form:

- **Neutrality.** The un-delayed charge sums to zero at a common instant; delayed, it does
  not — by up to 40% of the charge present, on the solenoid. The excess is removed *in
  proportion to |λ| itself*, which restores the sum to zero exactly while staying zero
  wherever the charge is still zero. A uniform subtraction would have put charge on
  segments the news had not reached, which is the very thing being fixed.
- **Continuity.** The current can no longer be a closed-form expression, because ∂λ/∂t now
  varies with each segment's own delay. It is re-derived by integrating ∂I/∂s = −∂λ/∂t
  along the wire. While part of the loop is still dark the constant of integration is fixed
  by requiring zero current there — there is no closed path yet, so nothing can circulate.
  Once the loop is lit the lumped I(t) takes over, blended across the moment it closes.

**The source.** Ohm's law along the wire, with the induced part of the field kept:

```
dV/ds  =  e(t) g_emf(s)  +  I(t) g_R(s)  +  İ(t) g_L(s)
```

with `g_L = -a(s)`, the tangential vector potential per unit current. Integrate once around
the loop and the three terms are the circuit equation,

```
e(t)  =  I R  +  İ L ,        L = ∮ a(s) ds
```

which is also exactly the condition for V to come back to itself after a lap. **L is not
assumed**: it is the same geometric kernel the charge solve uses, weighted by `t·t′`
instead of by 1. For a circular loop it comes out within a percent or so of
μ₀R[ln(8R/a) − 2], which is checked in the tests and appears nowhere in the code.

Because `dV/ds` is linear in (e, I, İ), so is the charge density. **Three solves of one
matrix** give basis densities, and from then on λ(s,t) is arithmetic — which is why moving
the clock costs nothing and the thing animates. The current then cannot be uniform along
the wire: continuity,

```
∂λ/∂t + ∂I/∂s = 0
```

fixes how it varies, and that variation is what radiates. Both are carried exactly rather
than left inconsistent.

**The fields.** Jefimenko's equations for a filament — Maxwell's equations solved for a
given source history, retardation included:

```
E = (1/4πε₀) Σ len [ λ/R² + λ̇/(cR) ] R̂  −  İ t̂ /(c²R)
B = (μ₀/4π)  Σ len [ I/R² + İ/(cR) ] t̂ × R̂
```

all evaluated at `t_r = t − R/c`. Causality is never imposed; it falls out of the lookup.

## What is approximated

**The response is a delayed quasi-static one, not a self-consistent solve.** The *shape* of
λ at each instant is still the instantaneous boundary-value solution, and I(t) still comes
from a lumped circuit equation; both are then delayed. That is sound when the switch-on
time is long compared with the time light takes to cross the loop, and it is being
stretched when it is not. The fields are then *exact for that source* — but this is not a
time-domain integral equation and does not pretend to be.

The numbers panel prints both times so you can see which regime you are in. With the
default rectangle, L/R is about 14 ns against 5.6 ns for light along the wire, which is a
comfortable margin. Drag the switch-on time down to a few hundred picoseconds and you are
looking at an illustration rather than a derivation.

**Everything else.** The LED is a plain resistor. The wire is a filament of radius `a`.
Each segment is a point element at its centre, which is also where its retarded time is
taken. There is no radiation reaction and no skin effect.

## What is on screen

The **colour map on the slice** is where the wavefront reads most clearly, so it is on by
default. Every colour ramp starts at exactly the background colour: on the steady benches
the bottom of the ramp was a dark tint, harmless because no point is ever unreached there,
but here a tinted floor would paint the unreached region as "a little bit of field" and
wipe out the boundary the whole bench is for.

The colour range is fixed to the **settled** field, not to the present instant. A range
recomputed each frame would rescale as the field grows and every moment would look equally
bright, which is the one thing that must not happen.

The dashed sphere has radius **ct** about the **battery** — the only thing that changes,
and so the only place news can start from. It is a boundary rather than a guide: the field
is strictly zero outside it, which is checked in the tests.

**E, B and E×B can all be shown at once**, and are independent toggles rather than a
one-of-three picker — **E** and **B** are on by default. Arrows and field lines layer
happily, each in its own ramp. The colour map cannot: two of them over one plane is mud, so
it takes the first field switched on.

Showing several fields is much cheaper than several times the cost, because one field
evaluation serves all of them — the colour ranges are gathered in a single sweep. Tracing
is the exception, since each field's lines follow their own path, so lines get shorter when
more than one is running.

**Arrows** sit on a lattice through the volume, and their length is a pixel count converted
back to a world length at that point's depth, so the projection foreshortens them — an
arrow pointing at the camera looks short, which is the depth cue. Full strength is 29 px at
1×, about a quarter of the on-screen gap between neighbouring lattice points.

Both sliders are logarithmic with 1× at the middle of their travel — arrow size 0.1× to
10×, line density 0.2× to 5×. Linear, they would spend most of the travel above 1× and
leave the small end unreachable.

**The movie draws exactly what a stopped frame draws.** Only dragging a handle and
scrubbing the clock get a coarser picture, because those have to answer the hand. Playback
does not: it is the thing worth looking at, so it pays full price and runs at whatever frame
rate that costs.

Which is not many frames a second. **Field line density costs more than it looks** — lines
pack in three dimensions, so halving the spacing is about seven times the lines, and the
seed lattice is what sets it, since once the separation rule is loose enough to stop binding
the count just tracks the lattice. Measured in the browser with **E** and **B** both up,
arrows and charge marks on:

| density | a frame | while dragging |
|---|---|---|
| 0.45× | 127 ms | 38 ms |
| 0.75× | 364 ms | 58 ms |
| **1.00×** (default) | 779 ms | 153 ms |

1× is about double the line count of the version before it. Turn it down, or show one
field at a time, if you would rather the movie were smooth than detailed.

Those absolute numbers are from the machine the table was taken on and the code as it stood
before the wire and sheet seeds; scale them by the fifth and the two fifths above. The
ratios are what travel between machines — a faster one here runs a 1× frame in about 130 ms.

**The clock keeps real time regardless.** It advances by the wall time each frame actually
took, so at the default 1 ns/s a 780 ms frame moves it 0.78 ns and the transient still plays
out over its true number of seconds, however few frames that is. The step is capped at two
seconds, which is only there to stop a backgrounded tab jumping the clock when it comes
back; the cap used to be 0.1 s, invisible while frames were 20 ms and a throttle to an
eighth speed once they were not.

**Field lines** are traced every frame, the movie included. Dragging and scrubbing get a
coarser seed lattice, shorter lines and a plain Euler step instead of the midpoint rule,
which halves the field evaluations that are the whole cost of tracing; the lattice is the
same shape either way, so lines do not jump about as the quality changes. Seeds in the dark
find no field and produce nothing, so the family grows outward on its own as the front
passes.

Seeds come from three places, and the lattice alone is not enough. It fills a **volume**,
which it has to: a line can be anywhere. But the region the loop encloses is a **sheet**
through that volume, with no volume of its own, so the lattice landed in it only by
accident — three seeds out of a hundred and twenty-five at the default density, and when
the lattice side came out even, which is what the moving picture uses, not one. The inside
of the loop drew nothing, and drew nothing most reliably while the movie was playing. That
is the worst place to lose: the **E** lines running from the + charge across to the − charge
are a good part of what the bench is for.

So the wire is seeded too — a ring of points a few radii off it, and the cone from the
loop's centroid over it, which is a surface the loop spans whatever it is bent into. No
notion of "the plane of the loop" appears anywhere, which is as well, since the helix and
saddle presets have not got one. The ring alone does not do it for **E×B**, which points
*along* the wire where it is close to it — energy running down the line — so a line started
there follows the wire out and never crosses the middle; that is what the cone is for.

This costs what it draws and no more: about a fifth more time for a stopped frame, which
buys about an eighth more lines, and about two fifths more on the dragging path, which had
no interior lines at all before. The per-line cost is unchanged.

The families are taken round-robin in proportion to their lengths rather than one after
another. The cap on lines cuts whatever is still in the list, so end-to-end let the first
family eat the whole budget: at 5× the wire ring took all 240 lines and the volume got
nothing at all, and the density slider had stopped adding coverage and was just packing
lines round the wire.

**How close two lines may come scales with the field.** A field line is a flux tube, and
there is no charge off the wire, so equal flux per line means the number crossing unit area
goes as |**F**| and the spacing as |**F**|^(−½) — which is the old convention that where the
lines crowd, the field is strong. A rule at a fixed fraction of the scene destroys it,
thinning hardest exactly where the crowding is the physics. Scaling the rule moves **E**
from |E|^0.37 to |E|^0.54 here, **B** from 0.33 to 0.48 and **E×B** from 0.02 to 0.43, for
the same line count and the same time. The sibling bench in
[`../EMviewer3D`](../EMviewer3D/) reaches 0.90 to 1.01 on the same measure; this one is
denser and its clamp binds harder, so it gets less of the way there.

**Wire the news has not reached is drawn dark.** That is the honest way to show where ct
has got to *along the wire*: it is continuous and exact, where the charge marks are
quantised and cannot be.

**Surface charge** marks sit at fixed positions, one per equal quantum of the *settled*
charge. What changes with time is how filled each one is: |λ| there now over |λ| there
once settled. They fade up in place, in the order the news reaches them, which is outward
from the battery.

Placing them by integrating |λ| at the current instant is the obvious thing and is what
this did at first. It is wrong, and visibly so: the point where the running total first
reaches one quantum moves inward as λ rises everywhere, so every mark slides toward
wherever the integration started. Anchored at the battery and run outward both ways, it
reads unmistakably as charge flowing *into* the battery, which is not happening — the
drift is an artefact of quantising a quantity that is growing. Fixing the positions removes
the motion entirely, and the count still carries the charge, because each mark is one
quantum drawn at its fill fraction.

The marks do not reach the front, and should not: λ is exactly zero there and takes the
whole switch-on time to come up, during which the front travels most of the way round the
loop. Measured, the freshly-lit wire behind the front carries a few thousandths of one
quantum. The dimmed wire is what marks the front.

## Controls

| | |
|---|---|
| Play / Restart | run the clock; the slider scrubs it |
| Switch-on time | how long the battery takes to come up |
| Slow motion | nanoseconds of circuit time per second of yours; 1 ns/s by default |
| drag the background | orbit |
| scroll / pinch | zoom |
| drag a white handle | reshape the wire, in the plane of the screen |
| hover | E, B, S and potential where the cursor's line of sight crosses the slice |

## Layout

```
index.html              markup
styles.css
src/solver4d.js         geometry, inductance, the circuit march, retarded fields (UMD)
src/view3d.js           orbit camera and projection (UMD)
src/app.js              canvas rendering, the clock, UI
test/physics4d.test.js  regression tests, no dependencies
```

## Tests

```
node test/physics4d.test.js
```

Checked against things known independently of this code:

- **Self-inductance** against μ₀R[ln(8R/a) − 2] for a thin ring, at three radii and wire
  gauges — within 2%.
- **The circuit march**: the residual of L İ + R I = e(t) along the whole history, taken
  with a five-point derivative so the test's own differencing does not dominate; and
  convergence in step size, which moves the answer by 4×10⁻¹¹ when the step is quartered.
- **Continuity**: ∂I/∂s + ∂λ/∂t around the whole loop, to 0.01% of peak λ̇; the analytic
  time derivatives against numerical ones; and neutrality at every instant.
- **Causality**: the field at a probe is *exactly* zero until the battery's light time and
  non-zero just after; a probe twice as far waits twice as long; and — the point of
  delaying the source — nothing arrives at the *nearest-wire* light time, which for the
  test probe is 1.8 ns earlier.
- **Radiation**: |B| falls as 1/R in the far field during the switch — measured in the
  plane of the loop, since the loop axis is the null of magnetic dipole radiation — and the
  settled axial field goes back to the 1/R³ of a dipole.
- **The steady state**: late in the history, λ and I settle to the DC solution, and the
  field matches what `../EMviewer3D` computes to machine precision when that sibling is
  present. It is not a dependency; the check skips if it is absent.
- **Field-line coverage**, which is about the drawing rather than the field. Random points
  are thrown into the scene and each is asked how far it is to the nearest drawn line of
  its own field. Points inside the loop are told apart from points outside by the solid
  angle the loop subtends at them — ±2π on the spanning surface, zero far off — so it needs
  no plane and works on the presets that have not got one, and both populations are held to
  the same band of distances from the wire so the comparison is of seeding and not of
  falloff. Inside must come out covered like outside, no point may be stranded, and where a
  line does pass close its tangent must agree with the field there to a degree or so — in
  sign as well as axis, since a whole line runs with +F from end to end. This is the check
  that the missing interior lines would have failed, and no test of the field could: nothing
  about the field was wrong.

## Ideas not built yet

- A genuinely self-consistent time-domain integral equation, which would remove the
  quasi-static assumption and let the switch-on time go arbitrarily short.
- The flux of **S** through a surface enclosing the LED, integrated over the transient
  against the energy the battery has delivered.
- Switching the battery *off*, which is where the inductive kick lives.
