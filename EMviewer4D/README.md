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

- **The field spreads at c.** Before the switch there is no field at all; then a front of
  it sweeps outward. At any moment the region that knows about the battery is bounded, and
  the colour map on the slice goes exactly black outside it — not dim, black. That boundary
  is not drawn or imposed: every history is zero before t = 0, so a point whose retarded
  time has not yet gone positive sums to nothing.
- **The current lags the EMF.** The loop is an inductor, and the trace in the rail shows
  the two curves failing to coincide: the EMF arrives in a couple of nanoseconds, the
  current takes L/R to follow. The gap between the curves is what builds the field.
- **Energy arrives last.** In **E×B** mode nothing flows at a point until both fields
  exist there, so the load cannot be fed faster than light can carry the news.

Turn the wire resistance up and L/R shortens; open the loop out and L grows. The numbers
panel prints L/R against the time light takes to run along the wire, which is the ratio
that decides whether the model below is on firm ground.

## How it's computed

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

**The source model is quasi-static.** I(t) comes from a lumped circuit equation using an
instantaneous inductance, and λ from an instantaneous boundary-value problem. That is sound
when the switch-on time is long compared with the time light takes to cross the loop, and
it is being stretched when it is not. The fields are then *exact for that source* — but
this is not a self-consistent time-domain integral equation and does not pretend to be.

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

The dashed sphere is a guide of radius **ct** about the origin. The circuit is not a point,
so news from its nearest part arrives a little sooner than the sphere suggests.

**Surface charge** spends one mark per equal quantum of charge, with the quantum fixed to
the settled state — so the marks thin out towards t = 0 rather than always filling the
wire, and the count is the charge at every moment.

There are no traced field lines here. Tracing is far too slow to redo every frame, and a
field line of a field that is still arriving is a confusing object anyway; the steady
benches are the place for those.

## Controls

| | |
|---|---|
| Play / Restart | run the clock; the slider scrubs it |
| Switch-on time | how long the battery takes to come up |
| Slow motion | nanoseconds of circuit time per second of yours |
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
- **Causality**: the field at a probe is *exactly* zero until R/c and non-zero just after,
  and a probe twice as far waits twice as long.
- **Radiation**: |B| falls as 1/R in the far field during the switch — measured in the
  plane of the loop, since the loop axis is the null of magnetic dipole radiation — and the
  settled axial field goes back to the 1/R³ of a dipole.
- **The steady state**: late in the history, λ and I settle to the DC solution, and the
  field matches what `../EMviewer3D` computes to machine precision when that sibling is
  present. It is not a dependency; the check skips if it is absent.

## Ideas not built yet

- A genuinely self-consistent time-domain integral equation, which would remove the
  quasi-static assumption and let the switch-on time go arbitrarily short.
- The flux of **S** through a surface enclosing the LED, integrated over the transient
  against the energy the battery has delivered.
- Switching the battery *off*, which is where the inductive kick lives.
