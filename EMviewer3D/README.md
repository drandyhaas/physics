# Field bench 3D

An interactive map of **E**, **B** and **S = E×B/μ₀** around a battery–wire–LED circuit
**in space** — the loop no longer has to lie in a plane, and neither does the field.
Orbit the scene, reshape the wire, and everything re-solves as you drag.

Open `index.html` in a browser. No build step, no dependencies, no WebGL; works over
`file://`.

This is the three-dimensional counterpart of [`../EMviewer`](../EMviewer/), which does the
same physics for a circuit confined to a plane. Neither depends on the other.

## Why it exists

The planar bench had to answer "what does **B** look like?" with ⊗ and ⊙ marks, because a
planar circuit's magnetic field points purely in and out of the screen. That is a true
statement about a special case, and it hides the thing worth seeing: **B** is a vector
field whose lines close into rings threaded on the wire. Give the circuit a third
dimension and they show up as what they are.

- Switch to **B** with the flat rectangle and the slice at z = 0. Arrows point straight up
  through the middle of the loop and straight down outside it, and the field lines arch
  over the wire and back — one continuous set of closed rings, with the wire threaded
  through them.
- Load the **solenoid** and the same field gathers into a bundle down the axis, splaying
  out at the ends and looping back around the outside. Nothing in the code knows what a
  solenoid is; it is five turns of the same filament.
- Tilt a flat loop out of the plane and every number stays put. The `tilt` preset is the
  rectangle rotated 40° about x, and its loop length and peak charge density agree with
  the flat one to nine figures — a rotation is not physics.

## What is on screen

The field is drawn **throughout the volume**, not on a surface. Field lines are seeded on
a 4×4×4 lattice filling a box around the circuit, so the picture has a front and a back
rather than being a plane seen at an angle.

**Arrows ride on the field lines.** A field line's tangent is the field direction, so an
arrow anchored on a line already points the right way, and no extra field evaluation is
needed to orient it. They are spaced by screen distance, so their density stays put as you
zoom, and the spacing scales with the size slider so large arrows do not collide. Arrows
and lines switch independently: arrows alone give the flow without the curves.

Arrow length is a pixel count converted back to a world length at that point's depth, so
the projection foreshortens it — an arrow pointing at the camera looks short, which is the
depth cue.

Seeds are taken in lattice order, and one is rejected where an existing line already
passes — a minimum separation of 0.12 of the scene radius. Order is what keeps the spacing
even: seeding strongest-field-first puts the lines at intervals that look arbitrary.
Because the order is regular, the separation rule only ever thins, leaving a regular subset
of the lattice. It thins where thinning is needed and nowhere else: the solenoid drops from
51 lines to 35, the saddle from 60 to 44, the flat loop stays at 28.

Tracing stops on two conditions besides running out of steps or leaving the box. A line
that returns to its own seed has **closed**, and is finished — without that test a ring is
traced past its seed and then traced again backwards, drawing it two or three times over.
And a line whose path length exceeds four times its own bounding-box diagonal has stopped
going anywhere: near a curved wire a **B** line does not close but winds helically around
the wire, spending every step adding length inside a small volume. A full circular ring
scores 2.2 by that measure and is untouched; the winders reach 5 or 6.

Lines still crowd inside the solenoid, and should. Field lines bunch where the flux is
concentrated, which in a coil is the entire point; that density is the physics and not a
placement artefact.

Weak field is not culled. The bottom of the colour ramp is the background colour and
opacity falls with it, so arrows in weak field fade out on their own: what thins out is
the field thinning out, not a threshold. The colour range itself comes from a coarse
sampling of the volume — from the geometry rather than from whichever lines happen to be
traced, so one colour keeps meaning one field strength while the circuit is dragged.

The **slice plane** remains, with a narrower job — a flat cut for the colour map, and the
surface the readout probes, since a cursor in a 3D scene is a ray rather than a point.
With `z = 0` and a flat loop it is the planar field bench, seen from wherever you happen
to be standing.

**Surface charge** spends one mark per equal quantum of charge, so the count of marks
along a stretch of wire is the charge on it. The loop is held neutral by the solve, so
there are always as many + as −.

## How it's computed

The wire is a closed filament in space carrying an unknown line charge density λ(s).
Ohm's law fixes the potential along it up to a constant, which is a boundary-value
problem:

```
Σᵢ Mⱼᵢ λᵢ  −  C  =  V_target(sⱼ)      for every segment j
Σᵢ λᵢ Lᵢ           =  0                (charge neutrality)
```

The Coulomb kernel 1/r was always three-dimensional, so **that solve is the planar one,
unchanged** — the tests check this directly, and λ agrees with the 2D solver to 6×10⁻¹⁵
for a circuit that lies in a plane. What is new is **B**. Instead of one component it is
the full vector field of each finite straight filament,

```
B = (μ₀I/4π) (L × r₁) [ L·r₁/|r₁| − L·r₂/|r₂| ] / |L × r₁|²
```

with r₁, r₂ from the segment ends to the field point. Beyond a few segment lengths the
point element is indistinguishable and much cheaper, so it takes over there.

Unlike **E**, **B** is not softened at the wire radius. The interior is already handled by
taking the field down linearly from the axis inside the metal, and softening on top of
that biases the field low well outside it — 6% at four wire radii, which the brute-force
test caught.

## Rendering

There is no WebGL and no scene graph. `view3d.js` is an orbit camera (azimuth, elevation,
distance, field of view) and a perspective projection; `app.js` projects everything by
hand and paints back to front onto a 2D canvas. Two ordering rules do all the work:

- The **wire and its charge beads** go into one list sorted by depth, so a solenoid passes
  behind itself correctly. Each segment carries its own dark halo, which is what makes a
  near part of the wire read as being in front of a far part.
- **Field lines and slice arrows** are split against a single depth threshold — the depth
  of the loop's centroid — and drawn in two passes, one before the wire and one after.
  Per-segment sorting of thousands of thin strokes would cost more than it shows.

Focal length comes from the shorter side of the canvas, so the view contains the circuit
at any window shape rather than pushing it out through the sides of a tall one.

## Controls

| | |
|---|---|
| drag the background | orbit |
| scroll / pinch | zoom |
| drag a white handle | reshape the wire, in the plane of the screen at that handle's depth |
| double-click the wire | insert a handle |
| hover | E, B, S and potential where the cursor's line of sight crosses the slice plane |

Dragging a handle moves it in the plane of the screen, so which way it goes depends on
where you are standing. Orbit first, then drag.

## Layout

```
index.html            markup
styles.css
src/solver3d.js       geometry, the charge solve, field evaluation (UMD: browser + node)
src/view3d.js         orbit camera and projection (UMD: browser + node)
src/app.js            canvas rendering and UI
test/physics3d.test.js  regression tests, no dependencies
```

## Tests

```
node test/physics3d.test.js
```

The solver is checked against things known independently rather than against itself:

- **Brute-force Biot–Savart.** Every segment is subdivided into 60 point elements and
  summed, with no code shared with the formula being tested. Agreement is better than
  0.1% in magnitude and 10⁻⁶ in direction.
- **Analytic limits.** The field beside a long run against μ₀I/2πd, and a circular loop
  against μ₀I/2R at its centre and μ₀IR²/2(R²+z²)^{3/2} on its axis — all within 1%.
- **Reduction to the planar solver.** When `../EMviewer` is present as a sibling, a planar
  circuit is solved both ways and compared: λ to 6×10⁻¹⁵, E to 4×10⁻¹⁶, B_z to 4×10⁻¹³,
  and the in-plane components of B are checked to vanish. The sibling is not a
  dependency; the check skips if it is not there.
- **Ohm's law at the surface.** Just outside resistive wire the tangential field is ρI to
  0.3%, as it must be, since E_t is continuous across the surface.

## Modelling choices worth knowing

- The LED is a plain resistor, with no forward-voltage knee and no diode behaviour.
- The wire is a filament of radius `a`; the interior is imposed analytically rather than
  resolved.
- Everything is static — no displacement current, no radiation, no switch-on transient.
- Field lines stop at a bounding sphere rather than closing, so a line that would run a
  long way before returning is drawn as far as it gets.

## Ideas not built yet

- Flux of **S** through a surface enclosing the LED, checked against I²R.
- Two separate loops, for mutual inductance geometry.
- Depth-sorted field lines, if a scene ever needs them badly enough to pay for it.
