# The Worldsheet and the Copper Plane

A seven-page broadside on a single claim:

> The quantum probability amplitude of interacting **strings** — extended
> one-dimensional objects sweeping out a two-dimensional worldsheet — is,
> integrand for integrand, the same mathematical object as the statistics of
> **thermal fluctuations of electrical power** across a two-dimensional
> resistive **copper plane**, with each particle's momentum becoming a current
> injected at a point on the board's one-dimensional edge.

**The claim holds.** It is the Coulomb-gas representation of the Koba–Nielsen
amplitude, and it is not new: it is essentially the analogue model of Fairlie
and Nielsen (1970). What this document adds is that every step is *measured*
rather than cited — including the substitution joining the two sides, which is
the step that actually carries the claim — and that the places where the
one-line version misleads are named rather than left to be noticed.

📄 **[string_copper_correspondence.pdf](string_copper_correspondence.pdf)**

| page | what is on it |
| --- | --- |
| 1 | the statement, the identity, and the chain of results it rests on |
| 2 | the graphic: one harmonic field, read two ways |
| 3 | the mechanism in five steps, each with its measurement |
| 4 | three consequences, and the evidence table |
| 5 | six things the one-line version hides |
| 6 | where gravity enters — two routes, both of them geometry |
| 7 | what the analogy says gravity *is*, and the one thing it cannot say |

## Why it works

A two-dimensional sheet has a logarithmic Green's function, so the power
dissipated by a Kirchhoff-balanced set of contacts on its insulating rim is

    P = −(2R_s/π) Σ_(i<j) I_i I_j ln|z_i − z_j|  +  κ Σ_i I_i²

and Johnson–Nyquist noise weights a fluctuation by `exp(−P/P_N)` with
`P_N = 8k_BT·Δf` for a noise bandwidth Δf. Exponentiating turns that sum of
logs into a product of powers — which is the Koba–Nielsen factor — with

    α′ = R_s / πP_N

the Regge slope as a sheet resistance in units of noise power. Fixing three
contacts and sliding the rest gives Euler's Beta integral: the Veneziano
amplitude.

**The string is visible on page 2, and it is the level curve.** Mandelstam's map
`ρ = Σ α_i ln(z − ζ_i)` has `Re ρ = τ` and `Im ρ = σ`, and that is the same
logarithmic potential the sheet carries — so reading the currents as light-cone
momenta, `τ = −(π/R_s)V`: the voltage *is* the string's time, a level curve is
the string at one instant, and the flow lines are position along it. An
equipotential meets an insulating rim at 90° (measured: the radial/tangential
gradient ratio falls 1.71 → 0.200 →
0.0189 → 0.0019 as |w| → 1), which is
precisely an open string's free-endpoint condition. Nothing shrinks at a
puncture: the conformal map is compressing an infinite stretch of τ into a
point, and length is not conformally meaningful on a worldsheet.

## Where gravity enters (pages 6–7)

Open strings give gauge fields — their lightest massless state is spin 1.
Gravity needs **closed** strings, whose lightest massless state is spin 2, and a
massless spin-2 field must couple to everything's stress-energy. On a worldsheet
open versus closed is literally rim versus interior, so there are two routes and
both are geometry:

- **Move a contact off the rim.** Measured factor of two between the two
  propagators: 1.99689. But a current injection is the wrong *kind* of
  insertion — it makes the tachyon. A graviton operator carries the same
  monopole factor *and* a factor bilinear in ∂X.
- **Punch a hole.** Each hole is one loop, so a washer is the one-loop diagram.
  Cut and roll it and it is a cylinder of length ln(b/a). Sliced by circles each
  slice is a closed string crossing rim to rim; sliced radially each slice is an
  open string going once around. Same washer.

And **a washer's resistance, in units of its sheet resistance, is the annulus's
conformal modulus** — the one shape parameter conformal maps cannot remove, and
the parameter a one-loop amplitude integrates over:

    R / R_s  =  (1/2π) ln(b/a)  =  t

Measured at a=120, b=600: `R = 0.256748` against `t = 0.256150`
— 0.23%, converging at first order (ratios
1.918, 1.787, 1.595), since the staircased rim
is the limiting error. The sheet itself converges at *second* order
(ratios 3.986, 3.816).

### What it teaches about gravity

1. **Gravity is not optional.** One washer, two slicings — the open-string loop
   and the closed-string exchange are the same object. A theory of gauge fields
   on a boundary is not consistent without closed strings inside it.
2. **Gravity reads the flow's anisotropy, not its strength.** Across D planes
   the graviton's bilinear is J^μJ^ν. Its *trace* is the total Joule
   dissipation — that is the **dilaton**. The *traceless* remainder, how the
   dissipation is distributed among the planes, is the graviton. Hence
   universality: any current in any plane contributes. Hence also why a graviton
   needs ≥2 planes — a symmetric traceless tensor in one dimension is zero.
3. **Its long range is the high-resistance end** of the modulus integral.

### And what it cannot say

This is the **bosonic** string. Its lightest closed state is not the graviton
but the closed *tachyon*, and the massless level holds a dilaton and a B-field
beside the graviton — so the long-cylinder limit is not graviton exchange alone.
More basically the copper is the *worldsheet*, so gravity arrives as a string
**state**: the analogy gives the particle-physics face of gravity and nothing of
the geometric one. No metric, no curvature, no Einstein equation, and no sense
in which a board bends the room it sits in.

## What the one-line version hides (page 5)

1. **One plane carries one spacetime dimension.** The obstruction is **rank**,
   not sign: 2α′p_i·p_j has rank up to D−1 while I_iI_j has rank one.
   (Indefiniteness is *not* the reason and is not general — at s=−2, t=−1.5 the
   same matrix is positive semi-definite.)
2. **The integral is still divergent**, and SL(2,ℝ) invariance needs **both**
   Σp = 0 and α′p² = 1: the exponent at puncture i is 2α′p_i² − 2 − 2α′p_i·(Σp).
   Breaking the shell alone gives 0.0652; breaking conservation
   alone gives 0.0362.
3. **The self-energy is a regulator, not the mass shell** — κΣI² is a constant
   per plane, α′p² = 1 is a condition on the full D-dimensional momentum.
4. **Rim versus interior is open versus closed.**
5. **This is the worldsheet CFT, not holography** — and not a least-resistance
   path either, though not for the obvious reason: Rayleigh's minimum is exactly
   what sits in the exponent, and what the amplitude integrates over is the
   contact *positions*.
6. **The two readings use different charges** — light-cone momentum on page 2,
   the whole momentum vector in the Koba–Nielsen weight.

Two errors worth flagging because they recur in popular retellings: the power
delivered by current sources carries **no factor of ½** (that belongs to
assembling static charges, and inserting it halves α′), and `exp(−P)` alone is
**dimensionally incomplete** — a power is not an exponent, which is where P_N
comes from.

## Files

| file | what it does |
| --- | --- |
| `string_copper_correspondence.pdf` | the document |
| `verify_string_copper.py` | verifies every claim; writes `measured.json` |
| `measured.json` | the measured values, consumed by the generator |
| `string_copper_correspondence.py` | builds the PDF (refuses to run without `measured.json`) |
| `pdfcanvas.py` | a dependency-free vector-PDF canvas: shapes, a tiny TeX-ish maths typesetter, paragraph flow with inline maths |

## Reproducing

```bash
python3 verify_string_copper.py          # ~2 min; prints every check, writes measured.json
python3 string_copper_correspondence.py  # builds the PDF from those values
```

The generator reads `measured.json`, so **the document cannot quote a number
the verification run did not produce.** A `--quick` run is marked as such and
the generator refuses it. The layout carries build-time assertions, so a column
or table that would overflow fails the build instead of shipping.

### What the verification actually does

- The sheet is a real lattice of unit resistors (so `R_s = 1 Ω/sq` exactly),
  and its insulating-edge Laplacian is diagonalised by the DCT-II, so `L V = I`
  is solved **exactly** — no iterative tolerance to trust.
- The rim Green's function is measured at three sheet sizes to show
  **convergence** to `−R_s/π`.
- The power law is fitted over 300 random neutral current patterns
  with the self-energy term kept explicit, so it cannot absorb a wrong
  coefficient.
- **The bridge is measured.** The substitution joining the two sides used to be
  the one step nothing checked. Now: take the measured lattice dissipation,
  subtract the fitted self-energy, divide by a noise scale, and compare against
  the string-side product at α′ = R_s/πP_N. Worst implied error on P over six
  configurations: 0.114% — the lattice error, not
  a defect in the dictionary. (The ratio itself deviates ~1%, because
  `ratio = exp(−δP/P_N)` and the exponential amplifies.)
- Rayleigh's least-dissipation principle: a divergence-free loop of strength ε
  costs exactly `4ε²`.
- The rim's Neumann condition is measured, not asserted.
- The string side is checked against closed forms: the Veneziano integral
  against Euler's B, the three cyclic orderings against the three Mandelstam
  channels, `B(k+1,n−k+1) = 1/[(n+1)C(n,k)]` on 38 integer
  pairs, `aB(a,1.5) → 1` as `a → 0`, and SL(2,ℝ) invariance under four Möbius
  maps with each condition broken separately.
- The washer is a separate Dirichlet problem (conjugate gradient on a masked
  lattice) at four resolutions.

### Requirements

macOS, Python 3, and:

```bash
pip install numpy scipy pyobjc-framework-Quartz pyobjc-framework-Cocoa
```

`pdfcanvas.py` draws through Quartz and lays out text through AppKit, so it
needs no LaTeX, no matplotlib and no font bundling — but it is macOS-only. The
figure's flow field and string slices are the real harmonic solution for the
unit disk, `V(w) = −(R_s/π) Σ_i I_i ln|w − ζ_i|`, traced by RK4; nothing in the
figure is a sketch.

## References

- Euler to Goldbach, 8 January 1730 — the integral interpolating the factorial
- J. B. Johnson, *Phys. Rev.* **32** (1928) 97; H. Nyquist, *Phys. Rev.* **32** (1928) 110
- G. Veneziano, *Nuovo Cimento* **A57** (1968) 190
- Z. Koba and H. B. Nielsen, *Nucl. Phys.* **B10** (1969) 633
- D. B. Fairlie and H. B. Nielsen, *An analogue model for KSV theory*, *Nucl. Phys.* **B20** (1970) 637
- A. M. Polyakov, *Quantum geometry of bosonic strings*, *Phys. Lett.* **B103** (1981) 207
