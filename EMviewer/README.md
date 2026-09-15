# Field bench

An interactive map of **E**, **B** and **S = E×B/μ₀** everywhere in the plane around a
battery–wire–LED circuit, where the wire can be grabbed and reshaped and the fields
re-solve as you drag.

Open `index.html` in a browser. No build step, no dependencies, works over `file://`.

## Why it exists

The usual circuit diagram says nothing about the space *around* the wire, which is where
almost everything interesting is happening. Two facts that are hard to believe until you
watch them:

- Along a resistanceless wire, **E** meets the metal at exactly a right angle. There is no
  tangential component, so no work is done on the carriers, and the interior field is zero
  even though current flows. Push the wire-resistance slider off zero and the field tilts
  over: that tilt *is* E = ρJ.
- The energy does not travel inside the wire. **S** leaves the battery sideways into empty
  space, glides along parallel to the wire without entering it, and turns inward at the LED.
  The wire steers the energy; the field carries it.

## How it's computed

The wire is a closed filament carrying an unknown line charge density λ(s). Ohm's law fixes
the potential along the filament up to a constant — flat along ideal wire, sloping across
the LED, stepping up across the battery — which turns the whole thing into a boundary-value
problem:

```
Σᵢ Mⱼᵢ λᵢ  −  C  =  V_target(sⱼ)      for every segment j
Σᵢ λᵢ Lᵢ           =  0                (charge neutrality)
```

`M` is the single-layer Coulomb operator, regularised at the wire radius `a`, with the
self-term `2·asinh(L/2a)` for a charged rod. Solve the bordered (N+1)×(N+1) system for λ and
the offset C, then **E** follows from Coulomb and **B** from Biot–Savart (exact
finite-segment formula near the wire, point-element approximation beyond 4 segment lengths).

Inside the conductor the filament model has nothing to say, but Ohm's law has everything:
the interior field is `E = ρJ` tangentially, so `fieldAt` blends to that within the wire
radius, and **B** is masked by `d²/a²` to recover the correct linear rise from the axis.

Because the circuit is planar, **E** lies in the plane and **B** is purely ±ẑ, so **S** is
in-plane too and the whole picture is honestly two-dimensional.

Nothing in the rendering is hand-drawn. The field-line shapes are the solution.

**The lines are evenly spaced, and that is deliberate.** The old convention is that line
density shows field strength, and where the field really is two-dimensional it is a theorem:
lines are flux tubes, so if each carries equal flux the number crossing unit length goes as
|**E**|. The bench in space earns it — there is no charge off the wire, so no divergence off
the wire, and [`../EMviewer3D`](../EMviewer3D/) holds its **E** lines to |E|^0.9 on purpose.

It is not available here, and not for want of trying. The kernel is the three-dimensional
1/r summed over a loop that happens to lie flat, so what is on screen is a *plane through* a
three-dimensional field. The in-plane divergence is −∂E_z/∂z, which does not vanish: flux
leaves the plane. Measured, 6.4% of the |**E**| around a 2 cm square, falling like the square
does — a finite divergence, not a rounding error, and it is checked in the tests. In-plane
lines are therefore not flux tubes: they converge and diverge for reasons that have nothing
to do with |**E**|.

So the spacing here is even and strength is carried by colour instead — on the arrows, and
on the background map. Forcing the spacing to track |**E**| would look like the theorem while
meaning nothing, which is worse than leaving it alone.

**Surface charge** turns λ — the thing actually solved for above — into the picture it
deserves. Rather than scale a glyph by λ, the overlay spends one mark per equal quantum of
charge, so the count of marks along any stretch of wire *is* the charge on that stretch and
marks per unit length is the density itself. Because the solve constrains the loop to be
neutral, the + and − marks always come out equal in number. This is the mechanism behind
both facts above: the field arrows leave the + beads and arrive at the − beads, and pushing
wire resistance off zero visibly drags the charge into a gradient along the whole run
instead of piling it at the terminals. It is on by default, here and on both sibling
benches: the charge is where the **E** lines begin and end, and without it the picture is
an unexplained tangle.

## Layout

```
index.html            markup
styles.css
src/solver.js         geometry, the charge solve, field evaluation (UMD: browser + node)
src/app.js            canvas rendering and UI
test/physics.test.js  regression tests, no dependencies
```

## Tests

```
node test/physics.test.js
```

Checks the solver against things known independently — Ohm's law, μ₀I/2πd, the direction of
energy flow — rather than against itself:

| check | result |
|---|---|
| boundary condition residual | 3.6e-15 V |
| charge neutrality | 3.1e-27 C |
| E inside ideal wire | 0 |
| E ∥ component just outside ideal wire | 0.5% of E⊥ |
| E inside LED vs IR/ℓ | 33.1 vs 33.1 V/m |
| B near wire vs μ₀I/2πd | 0.717 vs 0.667 µT |
| S just outside ideal wire | parallel, 1.2% inward leakage |
| S at the LED / at the battery | inward / outward |
| E inside resistive wire vs ρI | 0.95 vs 0.95 V/m |
| N=140 vs N=380 field at 2 cm | 0.5% |

The residual 0.5–1% anisotropies are filament discretisation, not physics; they shrink with N.

## Controls

| | |
|---|---|
| drag a white handle | reshape the wire |
| double-click canvas | insert a handle |
| alt-click / long-press a handle | delete it |
| hover anywhere | live E, B, S and potential at that point |

The view is fitted to the window, so the whole loop stays framed at any size or
aspect; 1 px = 1 mm is the ceiling, reached on a wide screen. The default loop is
about 1.7 m of wire.

## Modelling choices worth knowing

- The LED is a plain resistor. There is no forward-voltage knee and no diode behaviour, so
  reversing the battery just reverses everything.
- The wire is a filament of radius `a`; the interior is imposed analytically rather than
  resolved. Very close to the wire (under about a segment length) the field is smoothed.
- Everything is static. Displacement current, radiation and the switch-on transient are all
  absent — this is the DC steady state only.
- The wire radius slider goes up to 12 mm mostly so the interior is large enough to look at.

## Ideas not built yet

- Flux integral of **S** through a surface enclosing the LED, checked against I²R.
- Two loops, to show mutual inductance geometry.
- A switch-on transient, which needs a time-dependent solver and would answer a different
  question: where the momentum comes from before the steady state exists.
