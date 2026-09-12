# The Worldsheet and the Copper Plane

A four-page broadside on a single claim:

> The quantum probability amplitude of interacting **strings** — extended
> one-dimensional objects sweeping out a two-dimensional worldsheet — is,
> integrand for integrand, the same mathematical object as the statistics of
> **thermal fluctuations of electrical power** across a two-dimensional
> resistive **copper plane**, with each particle's momentum becoming a current
> injected at a point on the board's one-dimensional edge.

**The claim holds.** It is the Coulomb-gas representation of the Koba–Nielsen
amplitude, and it is not a new observation: it is essentially the analogue model
of Fairlie and Nielsen (1970). What this document adds is that every step is
*measured* rather than cited, and that the four places where the one-line
version misleads are stated explicitly.

📄 **[string_copper_correspondence.pdf](string_copper_correspondence.pdf)** —
page 1 the statement, page 2 the graphic, page 3 the mechanism, page 4 the
caveats and the evidence table, page 5 where gravity enters.

## Why it works

A two-dimensional sheet has a logarithmic Green's function. So the power
dissipated by a Kirchhoff-balanced set of contacts on its insulating rim is

    P = −(2R_s/π) Σ_{i<j} I_i I_j ln|z_i − z_j|  +  κ Σ_i I_i²

and Johnson–Nyquist noise weights a fluctuation by `exp(−P/P_N)` with
`P_N = 4k_BT/τ`. That exponential turns the sum of logs into a product of
powers — which is exactly the Koba–Nielsen factor `∏|z_i−z_j|^(2α′p_i·p_j)`,
with

    α′ = R_s τ / 4π k_B T

the Regge slope as a sheet resistance measured in units of noise power. Fixing
three contacts and sliding the rest then gives the Euler Beta function: the
Veneziano amplitude.

## What the short version hides

1. **One plane carries one spacetime dimension.** `p_i·p_j` is a Minkowski
   inner product, so the exponent matrix is indefinite (eigenvalues −3, 0, 5, 6
   at a four-tachyon point); a real current can only make `I_iI_j`, which is
   rank one and positive semi-definite. *D* dimensions need *D* superposed
   planes, the timelike one carrying imaginary current.
2. **The integral as written is divergent** — it still contains the infinite
   volume of SL(2,ℝ). Three punctures must be held fixed.
3. **The self-energy and the mass shell are the same object**: a pad's finite
   radius on one side, normal ordering and `α′p² = 1` on the other.
4. **This is the *worldsheet* CFT, not AdS/CFT.** The copper plane is the
   worldsheet, not spacetime. And the amplitude is not a least-resistance path:
   Rayleigh's principle gives only the single most probable configuration,
   whereas the amplitude sums over all of them.

## Where gravity enters (page 5)

Closed strings carry a massless spin-2 state — the graviton. Open strings do
not; their first massless state is a spin-1 gauge boson. So the rim/interior
distinction *is* the gauge-theory/gravity divide, and the measured factor of
two between the two propagators is that divide showing up in a lattice Green's
function.

Two things have to be said carefully, and the page says both:

- **A current injection is not a graviton.** `exp(ip·X)` is a monopole source
  and its state is the tachyon. A graviton vertex operator is bilinear in the
  worldsheet field, so it injects *no net current at all* — it couples to
  something quadratic in the local current, which in copper is the Joule
  dissipation density `R_s|J|²`.
- **Punching a hole is the cleaner route.** Copper with a hole has two rims:
  that is the annulus, which is the one-loop open-string diagram. Slice it by
  circles and each slice is a closed string propagating rim to rim (the
  graviton channel); slice it radially and each slice is an open string going
  once around the loop.

And the part that is exactly true and measurable: **a washer's resistance, in
units of its sheet resistance, is the annulus's conformal modulus** — the one
shape parameter conformal maps cannot remove, and the very parameter a
one-loop string amplitude integrates over.

    R / R_s  =  (1/2π) ln(b/a)  =  t

Measured at a=120, b=600: `R = 0.256748` against `t = 0.256150` — 0.23%, and
first order in the lattice spacing (the staircased rim is the only error, and
the error halves each time the rim is resolved twice as finely). Since the
cylinder's conformal length is `ln(b/a) = 2πR/R_s`, resistance is *separation*:
a high-resistance washer is the long-cylinder limit where only the lightest
closed states survive the trip, i.e. long-distance graviton exchange. As
b/a → 1 the open-string loop takes over instead.

What page 5 does **not** claim: that a circuit board gravitates beyond the pull
of its mass, that any measurement on copper could test string theory, or that
this is holography. The copper is the worldsheet, so this graviton is a string
state on an auxiliary two-dimensional surface — not curvature of the room the
board sits in.

Two errors worth flagging, because they are easy to make and both appear in
popular retellings: the power delivered by current sources carries **no factor
of ½** (that belongs to the energy of assembling static charges, and inserting
it halves `α′`), and `exp(−P)` alone is **dimensionally incomplete** — a power
is not an exponent, which is where `P_N` comes from.

## Files

| file | what it does |
| --- | --- |
| `string_copper_correspondence.pdf` | the document |
| `verify_string_copper.py` | verifies every claim; writes `measured.json` |
| `measured.json` | the measured values, consumed by the generator |
| `string_copper_correspondence.py` | builds the PDF (refuses to run without `measured.json`) |
| `pdfcanvas.py` | a small dependency-free vector-PDF canvas: shapes, a tiny TeX-ish maths typesetter, and paragraph flow with inline maths |

## Reproducing

```bash
python3 verify_string_copper.py          # ~75 s; prints every check, writes measured.json
python3 string_copper_correspondence.py  # builds the PDF from those values
```

The generator reads `measured.json`, so **the document cannot quote a number
the verification run did not produce.** A `--quick` verification run is marked
as such and the generator refuses it.

### What the verification actually does

Nothing is taken on citation:

- The resistive sheet is a real lattice of unit resistors (so `R_s = 1 Ω/sq`
  exactly). Its insulating-edge Laplacian is diagonalised by the DCT-II, so
  `L V = I` is solved **exactly** — there is no iterative tolerance to trust.
- The rim Green's function is measured at three sheet sizes to show
  **convergence** to `−R_s/π`, since the only error sources are the lattice at
  short range and wall images at long range.
- The power law is fitted over 300 random neutral current patterns with the
  self-energy term kept explicitly, so it cannot absorb a wrong coefficient.
- Rayleigh's least-dissipation principle is checked by adding a
  divergence-free loop and confirming the cost is exactly `4ε²`.
- The washer is solved as a separate Dirichlet problem (conjugate gradient on
  a masked lattice, inner rim at 1 V and outer at 0 V) at four resolutions, to
  show `R/R_s` converging on the annulus modulus.
- The string side is checked against closed forms: the Veneziano integral
  against Euler's B, the three cyclic orderings against the three Mandelstam
  channels, `B(k+1,n−k+1) = 1/[(n+1)C(n,k)]` on 38 integer pairs, and SL(2,ℝ)
  invariance of the full measure under four Möbius maps — on the mass shell
  (ratio 1.000000000000) and off it (0.0570), which shows the gauge symmetry
  *is* the mass-shell condition.

### Requirements

macOS, Python 3, and:

```bash
pip install numpy scipy pyobjc-framework-Quartz pyobjc-framework-Cocoa
```

`pdfcanvas.py` draws through Quartz and lays out text through AppKit, so it
needs no LaTeX, no matplotlib and no font bundling — but it is macOS-only.
The figure's flow field is the real harmonic solution for the unit disk,
`V(w) = −(R_s/π) Σ_i I_i ln|w − ζ_i|`, traced by RK4; it is not a sketch.

## References

- Euler to Goldbach, 8 January 1730 — the integral representation interpolating the factorial
- J. B. Johnson, *Phys. Rev.* **32** (1928) 97; H. Nyquist, *Phys. Rev.* **32** (1928) 110
- G. Veneziano, *Nuovo Cimento* **A57** (1968) 190
- Z. Koba and H. B. Nielsen, *Nucl. Phys.* **B10** (1969) 633
- D. B. Fairlie and H. B. Nielsen, *An analogue model for KSV theory*, *Nucl. Phys.* **B20** (1970) 637
- A. M. Polyakov, *Quantum geometry of bosonic strings*, *Phys. Lett.* **B103** (1981) 207
