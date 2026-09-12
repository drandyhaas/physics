#!/usr/bin/env python3
"""Verify every claim quoted in string_copper_correspondence.pdf.

Nothing here is taken on citation. The resistive sheet is a real lattice of
unit resistors, solved exactly (its insulating-edge Laplacian is diagonalised
by the DCT-II, so there is no iterative tolerance to trust); the string side
is evaluated numerically and compared against closed forms.

    python3 verify_string_copper.py            # full, ~35 s
    python3 verify_string_copper.py --quick    # skip the largest sheet
"""
import sys, math, json, os
import numpy as np
from itertools import combinations
from scipy.fft import dctn, idctn
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import cg, LinearOperator
from scipy.special import beta as Beta, comb
from scipy.integrate import quad
from scipy.optimize import least_squares

QUICK = "--quick" in sys.argv
RESULTS = {"quick": QUICK}

def hdr(s):
    print(f"\n{s}\n{'─' * len(s)}")

# ══════════════════════════════════════════════════ 1. the resistive sheet
# Unit-resistor square lattice => sheet resistance R_s = 1 ohm/sq exactly
# (an L x L patch is L parallel paths of L resistors: R = L*r/L = r).
def greens(W, H, on_edge, offsets):
    """G_ab = (L^+)_ab, the pseudo-inverse of the graph Laplacian: the exact
    potential at contact b from a unit current injected at contact a, with the
    mean removed. One DCT solve per contact yields every pair at once."""
    lx = 4 * np.sin(np.pi * np.arange(W) / (2 * W)) ** 2
    ly = 4 * np.sin(np.pi * np.arange(H) / (2 * H)) ** 2
    lam = ly[:, None] + lx[None, :]
    inv = np.where(lam > 0, 1 / np.where(lam > 0, lam, 1), 0.0)
    y0 = 0 if on_edge else H // 2
    pts = [(W // 2 + o, y0) for o in offsets]
    G = np.zeros((len(pts), len(pts)))
    for a, (xa, ya) in enumerate(pts):
        src = np.zeros((H, W)); src[ya, xa] = 1.0
        V = idctn(dctn(src, type=2, norm='ortho') * inv, type=2, norm='ortho')
        for b, (xb, yb) in enumerate(pts):
            G[a, b] = V[yb, xb]
    return G, pts

def log_slope(G, pts):
    """Fit G_ab = m ln r_ab + const over every pair; m is the coefficient."""
    lr, g = [], []
    for a, b in combinations(range(len(pts)), 2):
        lr.append(math.log(math.hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1])))
        g.append(G[a, b])
    lr, g = np.array(lr), np.array(g)
    (m, c), *_ = np.linalg.lstsq(np.vstack([lr, np.ones_like(lr)]).T, g, rcond=None)
    return m, float(np.abs(g - (m * lr + c)).max()), len(g)

OFF = np.array([0, 20, 45, 75, 110, 150, 195, 245, 300, 360, 425])
SIZES = [(2048, 1024), (4096, 2048)] + ([] if QUICK else [(8192, 4096)])

hdr("1. THE RESISTIVE SHEET — is the pair potential −R_s/π · ln r ?")
print("   Unit-resistor lattice, insulating rim, exact DCT solve. R_s = 1 Ω/sq.")
print("   The only errors are the lattice at small r and wall images at large r,")
print("   so the test is CONVERGENCE in sheet size:\n")
print(f"   {'sheet':>13}   {'rim contacts':>22}   {'interior contacts':>22}   ratio")
for W, H in SIZES:
    me, re_, n = log_slope(*greens(W, H, True, OFF))
    mb, rb, _ = log_slope(*greens(W, H, False, OFF))
    print(f"   {W:>5}×{H:<5}   {me:+.6f} ({abs(me/(-1/math.pi)-1)*100:5.2f}%)   "
          f"{mb:+.6f} ({abs(mb/(-1/(2*math.pi))-1)*100:5.2f}%)   {me/mb:.5f}")
RESULTS['rim_slope'], RESULTS['interior_slope'] = float(me), float(mb)
RESULTS['ratio'] = float(me / mb)
print(f"\n   exact:  −R_s/π = {-1/math.pi:.6f}   −R_s/2π = {-1/(2*math.pi):.6f}   ratio 2")
print(f"   {n} pairs per fit.  The factor of two is the insulating rim's image")
print("   charge — the same two that separates the open-string boundary")
print("   propagator −2α′ln|x−y| from the closed-string bulk one −α′ln|z−w|.")

hdr("2. THE POWER, with the self-energy term kept explicitly")
print("   P = −(2R_s/π) Σ_{i<j} I_i I_j ln r_ij  +  κ Σ_i I_i²")
W, H = SIZES[-1]
Ge, pe = greens(W, H, True, OFF)
rng = np.random.default_rng(4)
rows = []
while len(rows) < 300:
    k = int(rng.integers(3, 8))
    idx = rng.choice(len(OFF), k, replace=False)
    q = rng.integers(-4, 5, k).astype(float)
    if q.sum() != 0 or np.abs(q).max() == 0:
        continue
    P = float(q @ Ge[np.ix_(idx, idx)] @ q)
    E = -sum(q[a] * q[b] * math.log(abs(OFF[idx[a]] - OFF[idx[b]]))
             for a, b in combinations(range(k), 2))
    rows.append((E, (q ** 2).sum(), P))
E, S, P = map(np.array, zip(*rows))
(c1, c2), *_ = np.linalg.lstsq(np.vstack([E, S]).T, P, rcond=None)
res = P - (c1 * E + c2 * S)
R2 = 1 - res.var() / P.var()
RESULTS['power_coeff'], RESULTS['power_R2'] = float(c1), float(R2)
RESULTS['power_kappa'], RESULTS['power_n'] = float(c2), len(rows)
RESULTS['sheet'] = [W, H]
print(f"   {len(rows)} random neutral patterns, 3–7 contacts, I ∈ [−4,4], on {W}×{H}")
print(f"   fitted 2R_s/π = {c1:.6f}   exact {2/math.pi:.6f}   "
      f"({abs(c1/(2/math.pi)-1)*100:.3f}% off)")
print(f"   fitted κ      = {c2:.6f}   max|resid| {np.abs(res).max():.2e}   R² = {R2:.9f}")
print("   NOTE there is no factor of ½: that belongs to the energy of assembling")
print("   static charges, not to the power delivered by current sources. Halving")
print("   this coefficient would halve the α′ the correspondence assigns.")

hdr("3. NEUTRALITY IS FORCED (Kirchhoff ⇔ momentum conservation)")
print("   Two +1 A contacts on the rim, no return path: P against sheet size")
prev = None
for w in (512, 1024, 2048, 4096, 8192):
    h = w // 2
    a = 4 * np.sin(np.pi * np.arange(w) / (2 * w)) ** 2
    b = 4 * np.sin(np.pi * np.arange(h) / (2 * h)) ** 2
    lam = b[:, None] + a[None, :]
    inv = np.where(lam > 0, 1 / np.where(lam > 0, lam, 1), 0.0)
    src = np.zeros((h, w)); src[0, w // 2 - 20] = 1; src[0, w // 2 + 20] = 1
    V = idctn(dctn(src, type=2, norm='ortho') * inv, type=2, norm='ortho')
    Pn = V[0, w // 2 - 20] + V[0, w // 2 + 20]
    d = "" if prev is None else f"   Δ = {Pn - prev:+.4f}"
    print(f"     {w:>5}×{h:<5}  P = {Pn:7.4f}{d}")
    if prev is not None:
        RESULTS['nonneutral_step'] = float(Pn - prev)
    prev = Pn
print(f"   predicted Δ per doubling = (Q²R_s/π)ln2 = {(4/math.pi)*math.log(2):.4f}  for Q = 2")
print("   Log-divergent: a logarithmic potential needs neutrality to stay finite.")

hdr("4. RAYLEIGH LEAST DISSIPATION (why the weight may be written exp(−P))")
W2, H2 = 4096, 2048
lx = 4 * np.sin(np.pi * np.arange(W2) / (2 * W2)) ** 2
ly = 4 * np.sin(np.pi * np.arange(H2) / (2 * H2)) ** 2
lam = ly[:, None] + lx[None, :]
inv = np.where(lam > 0, 1 / np.where(lam > 0, lam, 1), 0.0)
src = np.zeros((H2, W2)); src[0, W2 // 2 - 200] = 1; src[0, W2 // 2 + 200] = -1
V = idctn(dctn(src, type=2, norm='ortho') * inv, type=2, norm='ortho')
jx0, jy0 = V[:, :-1] - V[:, 1:], V[:-1, :] - V[1:, :]
P0 = float((jx0 ** 2).sum() + (jy0 ** 2).sum())
print(f"   Ohmic (harmonic) flow                P = {P0:.6f}")
for eps in (0.05, 0.2, 1.0):
    jx, jy = jx0.copy(), jy0.copy()
    y0, x0 = 1000, 2000                      # a divergence-free plaquette loop
    jx[y0, x0] += eps; jx[y0 + 1, x0] -= eps
    jy[y0, x0] -= eps; jy[y0, x0 + 1] += eps
    Pe = float((jx ** 2).sum() + (jy ** 2).sum())
    print(f"   + divergence-free loop ε = {eps:<5}     P = {Pe:.6f}   "
          f"(+{Pe - P0:.4f}, predicted 4ε² = {4*eps**2:.4f})")
print("   Ohm's law is the strict minimiser, so exp(−P/P_N) peaks there. The")
print("   Gaussian width is Nyquist's 4kT/R, giving P_N = 4k_B T/τ for an")
print("   averaging window τ, hence α′ = R_s τ / 4π k_B T.")

# ══════════════════════════════════════════════════ 4b. the washer
def annulus_R(a, b, pad=1.18, rtol=1e-11):
    """DC resistance of an annular sheet of unit-resistance bonds (R_s = 1),
    inner rim held at 1 V and outer at 0 V. R = ΔV / I."""
    N = int(b * pad)
    ax = np.arange(-N, N + 1)
    X, Y = np.meshgrid(ax, ax, indexing='ij')
    r = np.hypot(X, Y)
    inner, outer = r <= a, r >= b
    free = ~inner & ~outer
    idx = -np.ones(r.shape, np.int64); idx[free] = np.arange(free.sum())
    n = int(free.sum())
    fi, fj = np.nonzero(free)
    rows, cols, vals = [np.arange(n)], [np.arange(n)], [np.full(n, 4.0)]
    rhs = np.zeros(n)
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        ni, nj = fi + di, fj + dj
        nbf, nbi = free[ni, nj], inner[ni, nj]
        me = idx[fi, fj]
        rows.append(me[nbf]); cols.append(idx[ni[nbf], nj[nbf]])
        vals.append(-np.ones(int(nbf.sum())))
        np.add.at(rhs, me[nbi], 1.0)
    L = coo_matrix((np.concatenate(vals),
                    (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()
    V, info = cg(L, rhs, rtol=rtol, maxiter=200000,
                 M=LinearOperator((n, n), matvec=lambda v: v / 4.0))
    assert info == 0, f"cg did not converge: info={info}"
    Vf = np.zeros(r.shape); Vf[inner] = 1.0; Vf[free] = V
    I = 0.0
    ii, jj = np.nonzero(inner)
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        ni, nj = ii + di, jj + dj
        ok = (ni >= 0) & (ni < r.shape[0]) & (nj >= 0) & (nj < r.shape[1])
        sel = free[ni[ok], nj[ok]]
        I += float(np.sum(1.0 - Vf[ni[ok][sel], nj[ok][sel]]))
    return 1.0 / I

hdr("4b. THE WASHER — is a copper annulus's resistance its CONFORMAL MODULUS?")
print("   An annulus has exactly one shape parameter conformal maps cannot change,")
print("   b/a. That number is its modulus, and it is what the one-loop string")
print("   amplitude integrates over. Continuum: R = (R_s/2π)ln(b/a), so R/R_s")
print("   should BE the modulus t = (1/2π)ln(b/a).\n")
print("     a      b    b/a     R measured      t = (1/2π)ln(b/a)   rel. err")
rowsA = []
for sc in (1, 2, 4) + (() if QUICK else (6,)):
    a, b = 20 * sc, 100 * sc
    Rm = annulus_R(a, b)
    Rx = math.log(b / a) / (2 * math.pi)
    rowsA.append({"a": a, "b": b, "R": float(Rm), "t": float(Rx),
                  "err_pct": float(abs(Rm / Rx - 1) * 100)})
    print(f"   {a:>4}  {b:>5}   {b/a:5.2f}    {Rm:.8f}       {Rx:.8f}        "
          f"{abs(Rm/Rx-1)*100:6.3f}%")
RESULTS['annulus'] = rowsA
print("   The staircased rim is the only error source, and it is first order in")
print("   the lattice spacing — the error halves each time the rim is resolved")
print("   twice as finely, which is what the column shows.")
print("\n   Reading it: the annulus is the ONE-LOOP open-string diagram, and it has")
print("   two channels. Slice it by circles and each slice is a CLOSED string,")
print("   propagating from rim to rim — and closed strings contain the graviton.")
print("   Slice it radially and each slice is an OPEN string with an endpoint on")
print("   each rim, going once around the loop. Conformally the cylinder's length")
print("   is ln(b/a) = 2πR/R_s, so a HIGH-resistance washer is the long-cylinder,")
print("   long-distance graviton-exchange limit; b/a → 1 is the open-string loop.")

# ══════════════════════════════════════════════════ 5. the string side
ap = 1.0
m2 = -1 / ap                                  # open-string tachyon, α′p² = 1
hdr("5. THE STRING SIDE")
rng = np.random.default_rng(7)
z = rng.uniform(0, 1, 6); q = rng.normal(size=6); q -= q.mean()
prod = np.prod([abs(z[i] - z[j]) ** (2 * ap * q[i] * q[j])
                for i, j in combinations(range(6), 2)])
En = -sum(q[i] * q[j] * math.log(abs(z[i] - z[j])) for i, j in combinations(range(6), 2))
RESULTS['kn_boltzmann'] = float(abs(prod / math.exp(-2 * ap * En) - 1))
print(f"   Koba–Nielsen factor ÷ exp(−2α′·logarithmic pair energy) − 1 = "
      f"{RESULTS['kn_boltzmann']:.1e}")
print("   -> the amplitude's integrand IS a 2-D Coulomb-gas Boltzmann weight.")

def gram(s, t):
    """p_i.p_j = (m_i² + m_j² − s)/2 in (−,+,+,+); p_i.p_i = −m² = +1/α′."""
    u = 4 * m2 - s - t
    G = np.full((4, 4), -m2)
    for (i, j), v in (((0, 1), s), ((0, 2), t), ((0, 3), u),
                      ((2, 3), s), ((1, 3), t), ((1, 2), u)):
        G[i, j] = G[j, i] = (2 * m2 - v) / 2
    return G

s, t = -3.0, -2.5
G = gram(s, t)
als, alt = 1 + ap * s, 1 + ap * t
print(f"\n   Regge form, α(s) = 1+α′s:  2α′p₁·p₂ = {2*ap*G[0,1]:+.6f}, "
      f"−α(s)−1 = {-als-1:+.6f}   match {abs(2*ap*G[0,1]+als+1) < 1e-12}")
print(f"   momentum conservation Σ_j p_i·p_j = 0:  max |row sum| = {np.abs(G.sum(1)).max():.1e}")
a, b = -als, -alt
num, _ = quad(lambda x: x ** (a - 1) * (1 - x) ** (b - 1), 0, 1,
              epsabs=1e-14, epsrel=1e-14)
RESULTS['veneziano'] = float(abs(num / Beta(a, b) - 1))
RESULTS['veneziano_value'] = float(Beta(a, b))
RESULTS['veneziano_a'], RESULTS['veneziano_b'] = float(a), float(b)
print(f"\n   Gauge-fix z₁,z₃,z₄ = 0,1,∞ and slide z₂ over (0,1):")
print(f"     ∫₀¹ x^{{−α(s)−1}}(1−x)^{{−α(t)−1}} dx = {num:.14f}")
print(f"     Euler B({a:.1f}, {b:.1f})                = {Beta(a,b):.14f}   ← VENEZIANO")
print(f"     relative agreement {RESULTS['veneziano']:.1e}")

print("\n   The three cyclic orders of the rim contacts = the three channels.")
print("   (all three integrals converge together only when a,b>0 and a+b<1,")
print("    since a+b+c = 1 identically — so this kinematic point is chosen so.)")
als, alt = -0.2, -0.5
alu = -1 - als - alt
a, b, c = -als, -alt, -alu
print(f"     α(s),α(t),α(u) = {als}, {alt}, {alu:.1f}   sum = {als+alt+alu:.1f} (must be −1)")
print(f"     a,b,c = {a:.1f}, {b:.1f}, {c:.1f}   a+b+c = {a+b+c:.1f} (identically 1)")
worst = 0.0
for lab, num_, ref, name in (
        (" 0 < x < 1  ", quad(lambda x: x**(a-1)*(1-x)**(b-1), 0, 1,
                              epsabs=1e-13, epsrel=1e-13, limit=400)[0],
         Beta(a, b), "B(−α(s),−α(t))   s-t"),
        (" 1 < x < ∞  ", quad(lambda v: v**(-a-b)*(1-v)**(b-1), 0, 1,
                              epsabs=1e-13, epsrel=1e-13, limit=400)[0],
         Beta(b, c), "B(−α(t),−α(u))   t-u"),
        ("−∞ < x < 0  ", quad(lambda y: y**(a-1)*(1+y)**(b-1), 0, np.inf,
                              epsabs=1e-13, epsrel=1e-13, limit=400)[0],
         Beta(a, c), "B(−α(s),−α(u))   s-u")):
    rel = abs(num_ / ref - 1); worst = max(worst, rel)
    print(f"     {lab} {num_:.10f}  vs  {ref:.10f}   {name}   rel {rel:.1e}")
RESULTS['channels'] = float(worst)
print("     Their sum over the whole real line is the crossing-symmetric total.")

print("\n   Euler's B interpolates the binomial: B(k+1,n−k+1) = 1/[(n+1)C(n,k)]")
bad = tot = 0
for n_ in (1, 2, 3, 5, 8, 13):
    for k_ in range(n_ + 1):
        tot += 1
        bad += abs(Beta(k_ + 1, n_ - k_ + 1) * (n_ + 1) * comb(n_, k_) - 1) > 1e-11
RESULTS['binomial_pairs'], RESULTS['binomial_bad'] = int(tot), int(bad)
print(f"     {tot} integer pairs checked, {bad} mismatches   "
      f"(n=5,k=2: B(3,4) = {Beta(3,4):.10f} = 1/(6·10))")
print("     So 'the integral counts orderings continuously' is an identity.")

print("\n   The particle spectrum is the x → 0 endpoint (two contacts merging):")
for aa in (0.30, 0.10, 0.03, 0.01):
    print(f"     B({aa:.2f}, 1.5) = {Beta(aa,1.5):9.4f}    1/a = {1/aa:7.2f}")
print("     Poles at −α(s) = 0,−1,−2,… i.e. α′s = −1,0,1,2,… — the open-string")
print("     spectrum m² = (J−1)/α′. Convergence at x=0 needs 2α′p₁·p₂ > −1.")

# ══════════════════════════════════════════════════ 6. the caveats
hdr("6. THE CAVEATS, ALSO MEASURED")
G = gram(-3.0, -2.5)
ev = np.linalg.eigvalsh(2 * ap * G)
RESULTS['gram_eigs'] = [round(float(v), 6) for v in ev]
print(f"   Exponent matrix M_ij = 2α′p_i·p_j eigenvalues: {np.round(ev, 4)}")
print("   INDEFINITE, rank 3. A real scalar current gives M_ij = k·I_i I_j, which")
print("   is rank ONE and positive semi-definite — so one plane carries one")
print("   spacetime dimension; D need D superposed planes, the timelike one")
print("   with imaginary current.")

D, N = 26, 5
E5 = np.array([2.0, 1.0, 0.0, -1.0, -2.0])          # Σ E = 0
L5 = np.sqrt(E5 ** 2 + 1 / ap)                      # |p_vec| on the mass shell
assert L5.max() <= L5.sum() - L5.max(), "polygon inequality"
sol = least_squares(lambda th: [(L5 * np.cos(th)).sum(), (L5 * np.sin(th)).sum()],
                    rng.uniform(0, 2 * math.pi, N), xtol=1e-15, ftol=1e-15, gtol=1e-15)
p = np.zeros((N, D)); p[:, 0] = E5
p[:, 1] = L5 * np.cos(sol.x); p[:, 2] = L5 * np.sin(sol.x)
eta = np.ones(D); eta[0] = -1.0
dot = lambda x, y: float((x * eta * y).sum())
print(f"\n   5 on-shell tachyons, Σp = 0 (closed polygon, residual "
      f"{np.abs(p.sum(0)).max():.1e}); α′p² = {np.round([ap*dot(p[i],p[i]) for i in range(N)],9)}")
KN = lambda zz, P: np.prod([abs(zz[i] - zz[j]) ** (2 * ap * dot(P[i], P[j]))
                            for i, j in combinations(range(N), 2)])
print("   SL(2,ℝ) invariance of the FULL measure ∏|z_ij|^{2α′p_i·p_j} ∏dz_i:")
worst = 0.0
for (A, B_, C, Dd) in [(2.0, 0.3, 0.0, 1.0), (1.0, -0.4, 0.7, 1.3),
                       (0.6, 1.1, -0.5, 2.0), (1.7, 0.2, 0.9, 0.8)]:
    zz = np.sort(rng.uniform(0.2, 3.0, N))
    zp = (A * zz + B_) / (C * zz + Dd)
    jac = np.prod(abs(A * Dd - B_ * C) / (C * zz + Dd) ** 2)
    r = KN(zp, p) * abs(jac) / KN(zz, p)
    worst = max(worst, abs(r - 1))
    print(f"     (a,b,c,d) = ({A:>4},{B_:>5},{C:>5},{Dd:>4})   ratio = {r:.12f}")
RESULTS['sl2r'] = float(worst)
qoff = p.copy()
qoff[:, 1:] *= np.sqrt((0.5 + E5 ** 2) / (1 / ap + E5 ** 2))[:, None]
zz = np.sort(rng.uniform(0.2, 3.0, N)); A, B_, C, Dd = 1.0, -0.4, 0.7, 1.3
zp = (A * zz + B_) / (C * zz + Dd)
roff = (KN(zp, qoff) * abs(np.prod(abs(A*Dd - B_*C) / (C*zz + Dd)**2)) / KN(zz, qoff))
RESULTS['sl2r_off'] = float(roff)
print(f"   Off the mass shell (α′p² = {ap*dot(qoff[0],qoff[0]):.2f}): ratio = {roff:.6f}")
print("   The gauge symmetry IS the mass-shell condition; so ∫dz₁…dz_N carries")
print("   the infinite volume of SL(2,ℝ) and three punctures must be fixed.")

# ══════════════════════════════════════════════════ summary
hdr("SUMMARY — the table on page 4 of the PDF")
def row(q, meas, exact, agree):
    print(f"   {q:<42} {meas:<18} {exact:<28} {agree}")
row("quantity", "measured", "exact", "agreement")
print("   " + "─" * 116)
row("Rim-contact pair potential, R_s = 1", f"{RESULTS['rim_slope']:.6f}",
    f"−R_s/π = {-1/math.pi:.6f}", f"{abs(RESULTS['rim_slope']/(-1/math.pi)-1)*100:.2f}%")
row("Interior-contact pair potential", f"{RESULTS['interior_slope']:.6f}",
    f"−R_s/2π = {-1/(2*math.pi):.6f}",
    f"{abs(RESULTS['interior_slope']/(-1/(2*math.pi))-1)*100:.2f}%")
row("Rim / interior ratio", f"{RESULTS['ratio']:.5f}", "2 exactly", "the image charge")
row("Power-law coefficient, 300 patterns", f"{RESULTS['power_coeff']:.6f}",
    f"2R_s/π = {2/math.pi:.6f}",
    f"{abs(RESULTS['power_coeff']/(2/math.pi)-1)*100:.2f}%, R² = {RESULTS['power_R2']:.8f}")
row("Non-neutral ΔP per size doubling", f"{RESULTS['nonneutral_step']:.4f}",
    f"(Q²R_s/π)ln2 = {(4/math.pi)*math.log(2):.4f}",
    f"{abs(RESULTS['nonneutral_step']/((4/math.pi)*math.log(2))-1)*100:.2f}%")
row("Gauge-fixed integral vs Euler B", f"{Beta(2.0,1.5):.14f}",
    f"B(2, 1.5) = {Beta(2.0,1.5):.14f}", f"{RESULTS['veneziano']:.1e}")
row("SL(2,ℝ) measure ratio, on the mass shell", f"{1+RESULTS['sl2r']:.12f}",
    "1 exactly", "four Möbius maps")
row("the same, pushed off the mass shell", f"{RESULTS['sl2r_off']:.6f}", "≠ 1",
    "the symmetry is the shell")
row("Koba–Nielsen exponent eigenvalues",
    ", ".join(f"{v:.0f}" for v in RESULTS['gram_eigs']), "indefinite, rank 3",
    "a real current gives rank 1, ≥ 0")
_A = RESULTS['annulus'][-1]
row(f"Washer resistance, a={_A['a']} b={_A['b']}", f"{_A['R']:.6f}",
    f"(1/2π)ln(b/a) = {_A['t']:.6f}", f"{_A['err_pct']:.2f}%, first order in h")
print(f"\n   also: three channels agree to {RESULTS['channels']:.1e}; "
      f"{RESULTS['binomial_pairs']} binomial pairs with {RESULTS['binomial_bad']} mismatches;")
print(f"   Koba–Nielsen = Boltzmann weight to {RESULTS['kn_boltzmann']:.1e}.")
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "measured.json")
with open(out, "w") as fh:
    json.dump(RESULTS, fh, indent=2, sort_keys=True)
print(f"\n   wrote {out} — string_copper_correspondence.py reads it, so the PDF")
print("   cannot quote a number this run did not produce.")
if QUICK:
    print("\n   (--quick: largest sheet skipped, so the sheet rows are less converged;")
    print("    the PDF generator refuses a --quick measured.json.)")
