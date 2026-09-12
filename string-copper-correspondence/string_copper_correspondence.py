#!/usr/bin/env python3
"""Build string_copper_correspondence.pdf.

A three-page broadside on the Coulomb-gas identity between the open-string
Koba-Nielsen amplitude and Johnson-Nyquist power fluctuations on a 2-D
resistive sheet. Every number quoted on page 3 comes from
verify_string_copper.py; the figure’s flow field is the real harmonic
solution for the unit disk, not a sketch.
"""
import json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pdfcanvas import Canvas, font_report

W, H = 792.0, 612.0                 # US Letter, landscape
M = 62.0                            # margin

PAPER = (0.991, 0.986, 0.977)
INK   = (0.098, 0.098, 0.110)
INK2  = (0.310, 0.310, 0.340)
MUTED = (0.520, 0.520, 0.550)
FAINT = (0.660, 0.655, 0.640)
RULE  = (0.830, 0.818, 0.786)
COOL  = (0.153, 0.337, 0.435)       # V < 0  (sinks)
WARM  = (0.733, 0.396, 0.141)       # V > 0  (sources)
MIDC  = (0.972, 0.960, 0.936)
STR   = (0.216, 0.247, 0.404)       # string-side accent
CU    = (0.659, 0.341, 0.118)       # copper-side accent

# ------------------------------------------------------------ measured values
# Read from measured.json so the document cannot quote a number that
# verify_string_copper.py did not actually produce. Regenerate with:
#     python3 verify_string_copper.py
_MJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "measured.json")
if not os.path.exists(_MJ):
    sys.exit("measured.json is missing — run: python3 verify_string_copper.py")
with open(_MJ) as _fh:
    MEAS = json.load(_fh)
if MEAS.get("quick"):
    sys.exit("measured.json came from a --quick run; rerun verify without it")

def num(x, d=6):
    """Fixed-point with a real minus sign (U+2212), not a hyphen."""
    return f"{x:.{d}f}".replace("-", "\u2212")

def sci(x, sig=1):
    """Returns MATHS MARKUP ($...$), so the exponent is typeset by STIXGeneral
    rather than by whatever face CoreText finds for a superscript digit."""
    e = int(math.floor(math.log10(abs(x))))
    man = f"{x / 10 ** e:.{sig}f}".rstrip("0").rstrip(".")
    return f"${man}\u00d710^{{{str(e).replace('-', '\u2212')}}}$"

def pct(x, ref, d=2):
    return f"{abs(x / ref - 1) * 100:.{d}f}%"

RIM, INT = MEAS["rim_slope"], MEAS["interior_slope"]
RIM_X, INT_X = -1 / math.pi, -1 / (2 * math.pi)
PWR, PWR_X = MEAS["power_coeff"], 2 / math.pi
NN, NN_X = MEAS["nonneutral_step"], (4 / math.pi) * math.log(2)
EIGS = ", ".join(f"{v:.0f}".replace("-", "\u2212") for v in MEAS["gram_eigs"])

# ------------------------------------------------------------------ field
# Unit disk, insulating rim, N point contacts ON the rim with sum I = 0.
# Then V(w) = -(R_s/pi) sum_i I_i ln|w - zeta_i| solves it EXACTLY: each term
# has constant normal derivative 1/2 on |w| = 1, so neutrality is precisely
# what makes the rest of the rim insulating.  R_s = 1 throughout.
ANG = np.deg2rad(np.array([158.0, 104.0, 30.0, -48.0, -126.0]))
CUR = np.array([+2.0, -1.0, +2.0, -2.0, -1.0])       # sum = 0
ZET = np.exp(1j * ANG)
assert abs(CUR.sum()) < 1e-12

def potential(w):
    w = np.asarray(w, dtype=complex)
    return -(1 / np.pi) * (CUR * np.log(np.abs(w[..., None] - ZET))).sum(-1)

def current(w):
    """J = -grad V / R_s, as a complex number (J_x + i J_y)."""
    w = np.asarray(w, dtype=complex)
    return np.conj((1 / np.pi) * (CUR / (w[..., None] - ZET)).sum(-1))

def trace(w0, rot=0.0, sign=1.0, h=0.005, nmax=5000, rmax=0.9965, stop=0.016,
          closed=True):
    """RK4 along the (rotated) unit current direction. rot=pi/2 -> equipotential."""
    rr = complex(math.cos(rot), math.sin(rot))
    def d(u):
        j = current(np.array([u]))[0] * rr
        m = abs(j)
        return 0j if m < 1e-14 else sign * j / m
    pts, w = [w0], w0
    for _ in range(nmax):
        k1 = d(w)
        if k1 == 0: break
        k2 = d(w + h * k1 / 2); k3 = d(w + h * k2 / 2); k4 = d(w + h * k3)
        step = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if step == 0: break
        w = w + h * step
        pts.append(w)
        if abs(w) > rmax: break
        if np.min(np.abs(w - ZET)) < stop: break
        if closed and len(pts) > 60 and abs(w - w0) < 1.2 * h: break
    return np.array(pts)

def streamlines(total=26):
    """Equal-flux seeding: a rim contact spreads its current over a half-disk,
    so equal flux means equal angle about the inward normal."""
    pos = CUR.clip(min=0).sum()
    out = []
    for i, q in enumerate(CUR):
        if q <= 0:
            continue
        n = max(2, int(round(total * q / pos)))
        base = math.atan2(-ZET[i].imag, -ZET[i].real)          # inward normal
        for k in range(n):
            a = base - math.pi / 2 + (k + 0.5) / n * math.pi
            w0 = ZET[i] + 0.055 * complex(math.cos(a), math.sin(a))
            if abs(w0) >= 0.997:
                continue
            p = trace(w0)
            if len(p) > 25:
                out.append(p)
    return out

def equipotentials(sl, count=11):
    """Seed orthogonal curves along the longest streamline."""
    src = max(sl, key=len)
    out = []
    for f in np.linspace(0.06, 0.94, count):
        w0 = src[int(f * (len(src) - 1))]
        if abs(w0) > 0.93 or np.min(np.abs(w0 - ZET)) < 0.10:
            continue
        a = trace(w0, rot=math.pi / 2, sign=+1.0)
        b = trace(w0, rot=math.pi / 2, sign=-1.0)
        if len(a) + len(b) > 40:
            out.append(np.concatenate([b[::-1], a]))
    return out

def rim_crossings(c, rr=0.988, nth=2400):
    """Angles on a near-rim circle where V = c. Every level curve of V meets the
    insulating rim (twice), so this finds both endpoints of every string slice."""
    th = np.linspace(0, 2 * np.pi, nth, endpoint=False)
    f = potential(rr * np.exp(1j * th)) - c
    out = []
    for k in range(nth):
        a, b = f[k], f[(k + 1) % nth]
        if not (a * b < 0):
            continue
        t0, t1 = th[k], th[k] + 2 * np.pi / nth
        f0 = a
        for _ in range(56):
            tm = 0.5 * (t0 + t1)
            fm = potential(np.array([rr * np.exp(1j * tm)]))[0] - c
            if f0 * fm < 0:
                t1 = tm
            else:
                t0, f0 = tm, fm
        out.append(rr * np.exp(0.5j * (t0 + t1)))
    return out

def string_slices(levels):
    """Each returned polyline is the open string at one instant of light-cone
    time: a level curve of V, with both ends free on the rim. Mandelstam’s map
    is rho = Σ α_i ln(z − ζ_i), so τ = Re rho = −(π/R_s)V and σ = Im rho."""
    out = []
    for c in levels:
        seeds = rim_crossings(c)
        used = [False] * len(seeds)
        for i, w0 in enumerate(seeds):
            if used[i]:
                continue
            best = None
            for sg in (+1.0, -1.0):                      # pick the inward branch
                probe = trace(w0, rot=math.pi / 2, sign=sg, h=0.004, nmax=4,
                              rmax=0.999, stop=0.004, closed=False)
                if len(probe) > 1 and abs(probe[-1]) < abs(probe[0]):
                    best = sg
            if best is None:
                continue
            path = trace(w0, rot=math.pi / 2, sign=best, h=0.004, nmax=9000,
                         rmax=0.988, stop=0.004, closed=False)
            if len(path) < 12:
                continue
            for j, wj in enumerate(seeds):              # its other end
                if j != i and abs(wj - path[-1]) < 0.05:
                    used[j] = True
            used[i] = True
            out.append((c, path))
    return out

def field_bitmap(n=460):
    """Premultiplied RGBA of V over the disk; row 0 = top. Clamped to a
    quantile so the log divergence at the contacts does not blow out the map."""
    ax = np.linspace(-1, 1, n)
    X, Y = np.meshgrid(ax, -ax)                 # -ax => row 0 is +y (top)
    Wg = X + 1j * Y
    rad = np.abs(Wg)
    V = np.zeros_like(rad)
    ok = rad < 0.9995
    V[ok] = potential(Wg[ok])
    core = rad < 0.93
    vmax = np.quantile(np.abs(V[core]), 0.93) or 1.0
    t = np.tanh(1.75 * np.clip(V / vmax, -1.6, 1.6)) / np.tanh(1.75)
    neg, mid, pos = np.array(COOL), np.array(MIDC), np.array(WARM)
    col = np.where(t[..., None] < 0,
                   mid + (neg - mid) * (-t[..., None]),
                   mid + (pos - mid) * (t[..., None]))
    a = (rad < 0.9975).astype(float)
    rgba = np.concatenate([col * a[..., None], a[..., None]], axis=-1)
    return (np.clip(rgba, 0, 1) * 255).astype(np.uint8).tobytes(), n

# ------------------------------------------------------------------ helpers
def headline(c, s, x, y, w, size=9.5, col=INK, role='bodyb'):
    """A run-in heading that may contain $maths$."""
    return c.flow(s, x, y, w, role, size, size * 1.32, col, just=False)

def eyebrow(c, s, x, y, col=MUTED, size=7.4):
    c.text(s.upper(), x, y, 'sansm', size, col, 'l', 1.55)

def rule(c, x1, x2, y, col=RULE, lw=0.6):
    c.line(x1, y, x2, y, col, lw)

def zigzag(c, p, q, n=6, amp=2.2, col=INK, lw=0.8):
    px, py = p; qx, qy = q
    L = math.hypot(qx - px, qy - py)
    if L < 1e-6: return
    ux, uy = (qx - px) / L, (qy - py) / L
    nx, ny = -uy, ux
    pts = [p]
    for k in range(1, 2 * n):
        t = (k / (2 * n)) * L
        s = amp * (1 if k % 2 else -1)
        pts.append((px + ux * t + nx * s, py + uy * t + ny * s))
    pts.append(q)
    c.polyline(pts, stroke=col, lw=lw)

def midarrow(c, pts, frac, col, size=3.6, lw=0.0):
    i = max(1, min(len(pts) - 1, int(frac * (len(pts) - 1))))
    (x0, y0), (x1, y1) = pts[i - 1], pts[i]
    c.arrowhead(x1, y1, math.atan2(y1 - y0, x1 - x0), size, col, slim=0.5)

def equation(c, cx, y, s=1.0):
    """The display identity, laid out in two passes so it can be centred."""
    BIG, SM, LIM = 27 * s, 14.2 * s, 6.9 * s
    pieces = []
    def add(kind, *a): pieces.append((kind, a))
    add('int', ); add('m', "\\rm{e}^{−P/P_N}\\,\\rm{d}z_1⋯\\rm{d}z_N"); add('gap', 13 * s)
    add('m', "="); add('gap', 13 * s)
    add('int', ); add('prod', ); add('m', "|z_i − z_j|^{2α′p_i·p_j}\\,\\rm{d}z_1⋯\\rm{d}z_N")
    tot = 0.0
    for k, a in pieces:
        if k == 'int':    tot += c.math_width("∫", BIG) + 2.5 * s
        elif k == 'prod': tot += max(c.math_width("∏", BIG * 0.86),
                                     c.math_width("1≤i<j≤N", LIM)) + 5 * s
        elif k == 'm':    tot += c.math_width(a[0], SM)
        else:             tot += a[0]
    x = cx - tot / 2
    for k, a in pieces:
        if k == 'int':
            x += c.math("∫", x, y - 5.5 * s, BIG, INK) + 2.5 * s
        elif k == 'prod':
            pw = c.math_width("∏", BIG * 0.86); lw_ = c.math_width("1≤i<j≤N", LIM)
            cw = max(pw, lw_) + 5 * s
            c.math("∏", x + (cw - pw) / 2, y - 1.5 * s, BIG * 0.86, INK)
            c.math("1≤i<j≤N", x + cw / 2, y - 17.0 * s, LIM, INK2, 'c')
            x += cw
        elif k == 'm':
            x += c.math(a[0], x, y, SM, INK)
        else:
            x += a[0]
    return tot

# ------------------------------------------------------------------ page 1
STATEMENT = (
    "Model an interacting particle not as a zero-dimensional point but as an extended "
    "one-dimensional string, sweeping out a continuous two-dimensional worldsheet in "
    "spacetime. Its quantum probability amplitude is then, integrand for integrand, the "
    "same mathematical object as the statistics of thermal fluctuations of electrical "
    "power across a two-dimensional resistive copper plane. In the correspondence the "
    "momenta of the incoming and outgoing particles become discrete electrical currents, "
    "injected and extracted at coordinate locations $z_i$ strictly along the "
    "one-dimensional outer boundary of the board. A two-dimensional sheet has a "
    "logarithmic Green’s function, so the power dissipated by any Kirchhoff-balanced set "
    "of edge currents is a sum of pairwise logarithms of the contact separations; and "
    "Johnson–Nyquist noise weights each microscopic thermal state by the exponential of "
    "its dissipated heat. Integrating those thermal states by sliding the injection "
    "points continuously along the boundary therefore reproduces the quantum path "
    "integral of the interacting strings, governed by the mathematical equivalence:")

PROVENANCE = [
    ("1730", "Euler", "the integral interpolating the factorial, to Goldbach, 8 January.", INK),
    ("1928", "Johnson · Nyquist", "measure, then explain, resistor noise. Phys. Rev. 32, 97 & 110.", CU),
    ("1968", "Veneziano", "the crossing-symmetric four-point amplitude. Nuovo Cimento A57, 190.", STR),
    ("1969", "Koba · Nielsen", "generalise it to N particles — the product above. Nucl. Phys. B10, 633.", STR),
    ("1970", "Fairlie · Nielsen", "the analogue model — the two sides joined. Nucl. Phys. B20, 637.", INK),
    ("1981", "Polyakov", "makes the worldsheet a 2-D statistical system. Phys. Lett. B103, 207.", STR),
]

def page1(c):
    c.page()
    c.rect(0, 0, W, H, fill=PAPER)
    eyebrow(c, "a verified analogy  ·  string amplitudes and Johnson noise", M, H - 56)
    c.text("The Worldsheet and the Copper Plane", M, H - 86, 'bodyb', 27.5, INK)
    c.text("A string scattering amplitude is the thermal-noise statistics of a resistive sheet",
           M, H - 106, 'bodyi', 12.8, INK2)
    rule(c, M, W - M, H - 120, INK, 0.9)

    colw = 322.0
    c.flow_columns(STATEMENT, [M, M + colw + 24], 470, colw, 'body', 10.0, 14.3, INK)

    # the identity
    c.rrect(M, 250, W - 2 * M, 86, 3, fill=(0.976, 0.968, 0.951), stroke=RULE, lw=0.6)
    equation(c, W / 2, 300, 1.0)
    for x, lab, subt, col in (
        (M + 168, "read as a copper plane",
         "the statistical weight of the dissipated power", CU),
        (W - M - 168, "read as a string worldsheet",
         "the Koba–Nielsen factor of the N-point amplitude", STR)):
        c.text(lab.upper(), x, 232, 'sansm', 7.2, col, 'c', 1.3)
        c.text(subt, x, 219, 'bodyi', 9.3, INK2, 'c')
    c.math("P \\rm{\\,dissipated power}\\;\\;\\; P_N \\rm{\\,the Nyquist power scale (p.\\,3)}"
           "\\;\\;\\; α′ \\rm{\\,the Regge slope (p.\\,3)}\\;\\;\\; N \\rm{\\,contacts}",
           W / 2, 204, 7.4, MUTED, 'c')
    c.text("an identity between INTEGRANDS; only N−3 of the N positions may actually "
           "be integrated (p. 4)", W / 2, 192, 'bodyi', 8.4, MUTED, 'c')

    # provenance — a named, published chain, not a new observation
    eyebrow(c, "not a new observation — the chain of results it rests on", M, 170)
    rule(c, M, W - M, 162, RULE, 0.6)
    for k, (yr, who, what, col) in enumerate(PROVENANCE):
        x = M + (k // 3) * (colw + 24)
        y = 144 - (k % 3) * 25
        c.text(yr, x, y, 'sansb', 9.2, INK, 'l', 0.4)
        wy = c.measure(yr, 'sansb', 9.2, 0.4)
        c.text(who, x + wy + 7, y, 'bodyb', 9.2, col)
        c.flow(what, x, y - 11.5, colw, 'body', 8.3, 11.0, INK2, just=False)
    rule(c, M, W - M, 62)
    c.text("Page 1 of 7  ·  the statement", M, 50, 'sansm', 7.2, MUTED, 'l', 1.2)
    c.text("Every number on pages 3, 4, 6 and 7 was measured, not cited",
           W - M, 50, 'sansm', 7.2, MUTED, 'r', 1.2)
    c.endpage()
# ------------------------------------------------------------------ page 2
def to_pt(w, cx, cy, r):
    return (cx + r * w.real, cy + r * w.imag)

def runs_inside(zs, rmax=0.985):
    """Split a complex polyline into the contiguous pieces inside |w| < rmax."""
    out, cur = [], []
    for z in zs:
        if abs(z) < rmax:
            cur.append(z)
        else:
            if len(cur) > 1: out.append(cur)
            cur = []
    if len(cur) > 1: out.append(cur)
    return out

def pad(c, cx, cy, r, ang, cur_, idx):
    """A rim pad, drawn straddling the board edge as a real edge pad does."""
    n = complex(math.cos(ang), math.sin(ang))
    P = (cx + r * n.real, cy + r * n.imag)
    c.save()
    c.ctx and None
    from Quartz import CGContextTranslateCTM, CGContextRotateCTM
    CGContextTranslateCTM(c.ctx, P[0], P[1]); CGContextRotateCTM(c.ctx, ang)
    c.rrect(-3.4, -6.6, 7.6, 13.2, 1.5, fill=(0.769, 0.478, 0.235), stroke=(0.40, 0.21, 0.09), lw=0.6)
    c.restore()
    src = cur_ > 0
    a0, a1 = (r + 26.0, r + 7.0) if src else (r + 7.0, r + 26.0)
    c.arrow(cx + a0 * n.real, cy + a0 * n.imag,
            cx + a1 * n.real, cy + a1 * n.imag, CU, 1.05, 4.8)
    lx = cx + (r + 34.0) * n.real; ly = cy + (r + 34.0) * n.imag
    c.math(f"I_{idx}", lx, ly - 3.2, 9.4, CU, 'c')
    ox = cx + (r + 47.0) * n.real; oy = cy + (r + 47.0) * n.imag
    c.text("in" if src else "out", ox, oy - 2.4, 'sansm', 6.6,
           CU if src else MUTED, 'c', 0.8)

def zlabel(c, cx, cy, r, ang, idx, col):
    n = complex(math.cos(ang), math.sin(ang))
    lx, ly = cx + (r - 14.0) * n.real, cy + (r - 14.0) * n.imag
    c.circle(lx, ly, 6.8, fill=(PAPER[0], PAPER[1], PAPER[2], 0.78))
    c.math(f"z_{idx}", lx, ly - 3.0, 8.6, col, 'c')

def string_legend(c, cx, ytop, w):
    """A ZOOM on the rim at one puncture. Locally the rim is straight and the
    level curves are semicircles about the puncture, so this is where the
    apparent collapse can be shown and explained."""
    c.text("DETAIL — THE RIM AT ONE PUNCTURE", cx, ytop, 'sansb', 6.4, MUTED, 'c', 1.25)
    y0 = ytop - 46                          # the rim, straight at this magnification
    px = cx - 46
    c.line(cx - 76, y0, cx + 76, y0, (STR[0], STR[1], STR[2], 0.95), 1.6)
    for rr in (6.0, 12.0, 19.0, 27.0):
        pts = [(px + rr * math.cos(t), y0 + rr * math.sin(t))
               for t in np.linspace(0, math.pi, 40)]
        c.polyline(pts, stroke=(STR[0], STR[1], STR[2], 0.92), lw=1.05)
        for e in (pts[0], pts[-1]):
            c.circle(e[0], e[1], 1.5, fill=STR)
    c.circle(px, y0, 2.6, fill=STR)
    c.math("z_i", px, y0 - 11.5, 7.6, STR, 'c')
    # tau increases AWAY from a source, i.e. radially across the arcs -- not
    # along the rim, which is what the arrow used to say
    a = math.radians(55.0)
    c.arrow(px + 7 * math.cos(a), y0 + 7 * math.sin(a),
            px + 33 * math.cos(a), y0 + 33 * math.sin(a),
            (0.10, 0.10, 0.12), 0.8, 4.0)
    c.math("τ", px + 40 * math.cos(a), y0 + 40 * math.sin(a) - 3, 8.0, INK, 'c')
    c.flow("A bold curve is the string at one instant, its two ends free on the rim. "
           "They crowd into the puncture — the marked point where a particle enters — "
           "because an infinite stretch of τ is squeezed into it. Nothing shrinks.",
           cx - w / 2, y0 - 20, w, 'bodyi', 7.9, 9.8, INK2, just=False)


def vignette_quantum(c, x, y, w, h):
    c.rrect(x, y, w, h, 2.5, fill=(0.977, 0.972, 0.963), stroke=RULE, lw=0.6)
    rng = np.random.default_rng(5)
    for amp, al in ((5.6, 0.20), (3.2, 0.34), (1.1, 1.0)):
        ph = rng.uniform(0, 2 * math.pi, 6)
        def dz(u, v):
            return (amp * (0.55 * math.sin(1.9 * u + ph[0]) + 0.45 * math.sin(2.7 * v + ph[1])),
                    amp * (0.50 * math.sin(2.3 * v + ph[2]) + 0.50 * math.sin(1.6 * u + ph[3])))
        col = (STR[0], STR[1], STR[2], al)
        for v in np.linspace(0.18, 0.82, 4):
            pts = []
            for u in np.linspace(0.08, 0.92, 40):
                d = dz(u * 6, v * 6)
                pts.append((x + u * w + d[0], y + v * h + d[1]))
            c.polyline(pts, stroke=col, lw=0.75 if al == 1.0 else 0.6)
        for u in np.linspace(0.10, 0.90, 7):
            pts = []
            for v in np.linspace(0.14, 0.86, 30):
                d = dz(u * 6, v * 6)
                pts.append((x + u * w + d[0], y + v * h + d[1]))
            c.polyline(pts, stroke=col, lw=0.55 if al == 1.0 else 0.45)

def vignette_thermal(c, x, y, w, h):
    c.rrect(x, y, w, h, 2.5, fill=(0.977, 0.972, 0.963), stroke=RULE, lw=0.6)
    rng = np.random.default_rng(11)
    nx, ny = 7, 4
    xs = np.linspace(x + 16, x + w - 16, nx)
    ys = np.linspace(y + 15, y + h - 15, ny)
    zig = {(1, 1), (3, 0), (4, 2), (2, 2), (5, 1)}
    for j, yy in enumerate(ys):
        for i, xx in enumerate(xs):
            if i < nx - 1:
                if (i, j) in zig:
                    zigzag(c, (xx + 3, yy), (xs[i + 1] - 3, yy), 5, 1.9,
                           (CU[0], CU[1], CU[2], 0.95), 0.75)
                else:
                    c.line(xx, yy, xs[i + 1], yy, (0.62, 0.58, 0.54), 0.55)
            if j < ny - 1:
                c.line(xx, yy, xx, ys[j + 1], (0.62, 0.58, 0.54), 0.55)
    for j, yy in enumerate(ys):
        for i, xx in enumerate(xs):
            c.circle(xx, yy, 1.0, fill=(0.40, 0.37, 0.34))
            if rng.random() < 0.62:
                a = rng.uniform(0, 2 * math.pi); L = rng.uniform(4.5, 9.0)
                c.arrow(xx, yy, xx + L * math.cos(a), yy + L * math.sin(a),
                        (COOL[0], COOL[1], COOL[2], 0.85), 0.7, 3.0)

def page2(c, bmp, n, sl, eq, slices):
    c.page()
    c.rect(0, 0, W, H, fill=PAPER)
    eyebrow(c, "the two complementary viewpoints", M, H - 50)
    c.text("One integrand, two readings", M, H - 76, 'bodyb', 20.5, INK)
    c.text("The same harmonic field on the same disk, with the same five boundary points",
           M, H - 93, 'bodyi', 10.6, INK2)
    rule(c, M, W - M, H - 105, INK, 0.8)

    R = 72.0
    CL, CR, CY = 202.0, 590.0, 344.0

    for cx, title, subt, col in ((CL, "The string worldsheet",
                                  "a 1-D string sweeping a 2-D surface", STR),
                                 (CR, "The resistive copper plane",
                                  "a 2-D sheet dissipating heat", CU)):
        c.text(title, cx, 494, 'bodyb', 12.4, col, 'c')
        c.text(subt, cx, 482, 'bodyi', 9.5, INK2, 'c')

    # ---- shared field, drawn identically on both sides
    for cx in (CL, CR):
        c.save(); c.clip_circle(cx, CY, R); c.alpha(0.93)
        c.image_rgba(bmp, n, n, cx - R, CY - R, 2 * R, 2 * R)
        c.restore()

    # ---- left: the SAME two curve families, emphasis swapped. Here the level
    # curves are bold — each one is the string at one instant of light-cone time
    # — and the flow lines, which mark position ALONG the string, are faint.
    c.save(); c.clip_circle(CL, CY, R * 0.999)
    for pth in sl:
        c.polyline([to_pt(z, CL, CY, R) for z in pth],
                   stroke=(0.30, 0.30, 0.36, 0.22), lw=0.4)
    for lev, pth in slices:
        pts = [to_pt(z, CL, CY, R) for z in pth]
        c.polyline(pts, stroke=(STR[0], STR[1], STR[2], 0.92), lw=1.05)
    for lev, pth in slices:                 # the string’s two FREE endpoints,
        for e in (to_pt(pth[0], CL, CY, R),  # drawn INSIDE the clip so none escapes
                  to_pt(pth[-1], CL, CY, R)):
            c.circle(e[0], e[1], 1.4, fill=STR)
    c.restore()
    c.circle(CL, CY, R, stroke=(STR[0], STR[1], STR[2], 0.95), lw=1.6)
    for i, a in enumerate(ANG):
        n_ = complex(math.cos(a), math.sin(a))
        src = CUR[i] > 0
        a0, a1 = ((R + 26.0, R + 7.0) if src else (R + 7.0, R + 26.0))
        c.arrow(CL + a0 * n_.real, CY + a0 * n_.imag,
                CL + a1 * n_.real, CY + a1 * n_.imag, STR, 1.05, 4.8)
        lx = CL + (R + 34.0) * n_.real; ly = CY + (R + 34.0) * n_.imag
        c.math(f"p_{i+1}", lx, ly - 3.2, 9.4, STR, 'c')
        ox = CL + (R + 47.0) * n_.real; oy = CY + (R + 47.0) * n_.imag
        c.text("in" if src else "out", ox, oy - 2.4, 'sansm', 6.6,
               STR if src else MUTED, 'c', 0.8)
        p_ = to_pt(np.exp(1j * a), CL, CY, R)
        c.circle(p_[0], p_[1], 2.6, fill=STR)
        zlabel(c, CL, CY, R, a, i + 1, (0.16, 0.18, 0.30))

    # ---- right: equipotentials, streamlines, board edge, pads
    c.save(); c.clip_circle(CR, CY, R * 0.999)
    for p in eq:
        c.polyline([to_pt(z, CR, CY, R) for z in p],
                   stroke=(0.24, 0.22, 0.18, 0.30), lw=0.45)
    for p in sl:
        pts = [to_pt(z, CR, CY, R) for z in p]
        c.polyline(pts, stroke=(0.13, 0.12, 0.13, 0.80), lw=0.72)
        for f in (0.34, 0.68):
            midarrow(c, pts, f, (0.13, 0.12, 0.13, 0.85), 3.5)
    c.restore()
    c.circle(CR, CY, R, stroke=(0.24, 0.13, 0.05), lw=1.7)
    c.circle(CR, CY, R - 3.2, stroke=(0.24, 0.13, 0.05, 0.35), lw=0.5)
    for i, a in enumerate(ANG):
        pad(c, CR, CY, R, a, CUR[i], i + 1)
        zlabel(c, CR, CY, R, a, i + 1, (0.30, 0.16, 0.06))

    # ---- the propagators, one under each disk
    c.math("⟨X^μ(z_i)X^ν(z_j)⟩ = −2α′η^{μν}\\,\\rm{ln}\\,|z_i − z_j|",
           CL, 220, 9.2, (0.16, 0.18, 0.30), 'c')
    c.math("V(w) = −(R_s/π) Σ_i I_i \\,\\rm{ln}\\,|w − ζ_i|",
           CR, 220, 9.2, (0.30, 0.16, 0.06), 'c')
    c.text("τ runs copper → teal, the way the current flows", CL, 208,
           'sansm', 6.8, MUTED, 'c', 1.1)
    c.text("the potential — exact, because the currents sum to zero", CR, 208,
           'sansm', 6.8, MUTED, 'c', 1.1)

    # ---- the centre column
    CX = 397.0
    string_legend(c, CX + 4, 478, 176)
    c.text("↔", CX, CY + 14, 'math', 30, INK, 'c')
    c.text("ONE INTEGRAND", CX, CY - 10, 'sansb', 7.4, INK, 'c', 1.5)
    c.text("TWO READINGS", CX, CY - 21, 'sansb', 7.4, INK, 'c', 1.5)
    rows = [("p_i", "I_i"), ("τ", "−(π/R_s)V"), ("α′", "R_s/πP_N"),
            ("\\rm{e}^{−S}", "\\rm{e}^{−P/P_N}"), ("Σp_i = 0", "ΣI_i = 0")]
    yy = CY - 44
    for a, b in rows:
        c.math(a, CX - 12, yy, 8.8, (0.16, 0.18, 0.30), 'r')
        c.text("↔", CX, yy, 'math', 8.6, MUTED, 'c')
        c.math(b, CX + 12, yy, 8.8, (0.30, 0.16, 0.06), 'l')
        yy -= 13.2

    # ---- the microscopic halves
    VW, VH, VY = 196.0, 72.0, 118.0
    vignette_quantum(c, CL - VW / 2, VY, VW, VH)
    vignette_thermal(c, CR - VW / 2, VY, VW, VH)
    for cx, head, cap, col in (
        (CL, "Microscopic fluctuations: quantum",
         "The embedding $X^μ$ fluctuates; each history weighted $\\rm{exp}(−S)$.", STR),
        (CR, "Microscopic fluctuations: thermal",
         "Every element rattles; each pattern weighted $\\rm{exp}(−P/P_N)$.", CU)):
        c.text(head, cx, VY - 12, 'sansb', 7.6, col, 'c', 1.0)
        cw_ = c.math_width(cap.replace('$', ''), 8.9)
        c.flow(cap, cx - cw_ / 2 - 6, VY - 24, cw_ + 14, 'bodyi', 8.9, 11, INK2, just=False)

    rule(c, M, W - M, 88)
    nr = ("Read each current as a light-cone momentum: Mandelstam’s map "
          "$ρ = Σ_i α_i \\rm{ln}(z − ζ_i)$ gives $τ = −(π/R_s)V$, so the voltage IS the "
          "string’s time, and $σ = \\rm{Im}\\,ρ$ is position along it.")
    nl = ("Colour is that same field on both sides — copper positive, teal negative. "
          "The flow lines are the current density $J = −∇V/R_s$; they are also the "
          "lines of constant $σ$.")
    for xx, t_ in ((M, nr), (M + 346, nl)):   # nr explains the LEFT panel
        n_ = len(c.wrap(t_, 'body', 8.8, 322))
        assert n_ <= 3, f"page-2 note column wraps to {n_} lines, into the folio"
        c.flow(t_, xx, 76, 322, 'body', 8.8, 11.8, INK2)
    c.text("Page 2 of 7  ·  the graphic", M, 40, 'sansm', 7.2, MUTED, 'l', 1.2)
    c.endpage()


# ------------------------------------------------------------------ page 3
def numbered(c, n, x, y, col):
    c.circle(x + 5.4, y + 3.0, 6.8, fill=(col[0], col[1], col[2], 0.13))
    c.text(str(n), x + 5.4, y, 'sansb', 8.2, col, 'c')

def steps_column(c, x, y, w, title, items, col, start=1):
    """Each item is a dict: head, body, optional disp (a centred display line),
    optional body2, optional chk."""
    eyebrow(c, title, x, y); y -= 19
    for k, it in enumerate(items, start):
        numbered(c, k, x, y, col)
        y = headline(c, it['head'], x + 18, y, w - 18) + 12.5
        yy = c.flow(it['body'], x + 18, y - 12.8, w - 18, 'body', 9.1, 12.4, INK2)
        if it.get('disp'):
            yy -= 4
            c.math(it['disp'], x + 18 + (w - 18) / 2, yy, it.get('dsize', 9.4), INK, 'c')
            yy -= 15
        if it.get('body2'):
            yy = c.flow(it['body2'], x + 18, yy - 1, w - 18, 'body', 9.1, 12.4, INK2)
        if it.get('chk'):
            c.text("✓", x + 18, yy + 1.0, 'sans', 8.0, (0.16, 0.45, 0.30))
            dx = c.measure("✓  ", 'sans', 8.0)
            yy = c.flow(it['chk'], x + 18 + dx, yy + 1.0, w - 18 - dx, 'bodyi', 8.5,
                        11.6, (0.16, 0.42, 0.29), just=False)
        y = yy - 10.0
    return y


STEPS_L = [
    dict(head="The string side is a product of powers.",
         body="An open string’s vertex operators — the insertions standing for each "
              "incoming or outgoing particle — sit on the one-dimensional boundary of the "
              "worldsheet, where the embedding has the propagator "
              "$⟨X^μ(x)X^ν(y)⟩ = −2α′η^{μν}\\,\\rm{ln}\\,|x − y|$. Contracting N of them "
              "turns that logarithm into a product, because the Gaussian average of a "
              "product of exponentials is the exponential of the sum of pairwise "
              "propagators:",
         disp="∏_{i<j}|z_i − z_j|^{2α′p_i·p_j}",
         body2="This is the bosonic string, and its N = 0 state is the tachyon "
               "($m^2 = −1/α′$, the bottom of the ladder in consequence 3).",
         chk="the product is exp of a logarithmic pair energy, to machine precision"),
    dict(head="A sheet’s rim Green’s function is $−R_s/π$ times a log.",
         body="A two-dimensional sheet has a logarithmic Green’s function, and a contact ON "
              "an insulating rim sees twice the interior value, because its image charge "
              "sits on top of it. That factor of two is not bookkeeping: it is the same two "
              "that separates the open-string boundary propagator from the closed-string "
              "bulk one.",
         chk=f"measured {num(RIM)} against $−R_s/π = {num(RIM_X)}$, {pct(RIM, RIM_X)}; "
             f"the error falls by {MEAS['rim_err_ratios'][0]:.1f}× per doubling, i.e. second order"),
    dict(head="Neutral rim currents dissipate a sum of logs.",
         body="Feeding that into $P = Σ_i I_i V(z_i)$ gives",
         disp="P = −(2R_s/π) Σ_{i<j} I_i I_j \\,\\rm{ln}\\,|z_i − z_j| + κ Σ_i I_i^2",
         body2="with $κ$ the self-energy coefficient, set by a pad’s finite radius. There is "
               "no factor of one half here — that belongs to the energy of assembling static "
               "charges, not to the power delivered by current sources, and inserting it "
               "would halve $α′$.",
         chk=f"fitted {num(PWR)} against $2R_s/π = {num(PWR_X)}$ over {MEAS['power_n']} "
             f"random neutral patterns, R² = {MEAS['power_R2']:.8f}"),
]
STEPS_R = [
    dict(head="Nyquist supplies the exponential — and its scale.",
         body="Thermal noise makes the contact currents a Gaussian whose quadratic form is "
              "exactly the dissipated power, so a fluctuation’s weight is "
              "$\\rm{exp}(−P/P_N)$ with $P_N = 8k_BTΔf$ for a noise bandwidth $Δf$ "
              "(equivalently $4k_BT/t$ for an averaging window $t$). Writing "
              "$\\rm{exp}(−P)$ alone is dimensionally incomplete: a power is not an "
              "exponent. Exponentiating step 3 now turns the sum of logs into a product of "
              "powers —",
         disp="\\rm{e}^{−P/P_N} = ∏_{i<j}|z_i − z_j|^{(2R_s/πP_N)I_iI_j}",
         dsize=9.0,
         body2="— and matching that exponent against $2α′p_i·p_j$ gives "
               "$α′ = R_s/πP_N$: the Regge slope, the constant fixing the string’s mass "
               "spectrum, as a sheet resistance in units of noise power. Its units are one "
               "over current squared, so $α′p^2$ is dimensionless.",
         chk=f"the substitution measured end to end on 6 configurations — worst implied "
             f"error on P, {MEAS['bridge_implied_P_err_pct']:.2f}%, which is the lattice error of step 2"),
    dict(head="Gauge-fix three contacts, then slide the rest.",
         body="Fixing $z_1,z_3,z_4 = 0,1,∞$ and integrating $z_2$ over $(0,1)$ collapses the "
              "N = 4 product to Euler’s Beta integral, with $α(s) = 1+α′s$ and $s$, $t$ the "
              "Mandelstam invariants — the squared collision energy and momentum transfer:",
         disp="∫_0^1 x^{−α(s)−1}(1−x)^{−α(t)−1}\\rm{d}x = B(−α(s),−α(t)) = "
              "Γ(−α(s))Γ(−α(t)) / Γ(−α(s)−α(t))",
         dsize=8.8,
         body2="where $B(a,b) = ∫_0^1 x^{a−1}(1−x)^{b−1}\\rm{d}x$, and $Γ$ is the integral "
               "Euler sent Goldbach in 1730 to interpolate the factorial. That is the "
               "Veneziano amplitude.",
         chk=f"the integral equals Euler’s B to {sci(MEAS['veneziano'])}"),
]
CONSEQ = [
    ("The three orderings are the three channels.",
     "Four contacts admit three distinct cyclic orders around the rim. Those are the three "
     "integration ranges of $x$, and they give $B(−α(s),−α(t))$, $B(−α(t),−α(u))$ and "
     "$B(−α(s),−α(u))$ — the s-t, t-u and s-u channels, whose sum is the crossing-symmetric total.",
     f"all three verified to {sci(MEAS['channels'])} or better"),
    ("The Beta function interpolates a binomial.",
     "$B(k+1, n−k+1) = 1/[(n+1)C(n,k)]$ exactly — so B is the continuous version of a "
     "combinatorial quantity, though it enters as the reciprocal. That is a separate fact "
     "from the three orderings, and not the reason there are three.",
     f"{MEAS['binomial_pairs']} integer pairs, {MEAS['binomial_bad']} mismatches"),
    ("The particle spectrum is two pads touching.",
     "The poles come from the endpoint $x → 0$, where two contacts merge. In the Γ form they "
     "are visible directly: $Γ(−α(s))$ blows up at $α(s) = 0,1,2,…$, i.e. $α′s = −1,0,1,…$, "
     "which is the open-string spectrum $m^2 = (J−1)/α′$.",
     f"$aB(a,1.5) → 1$ as $a → 0$: {', '.join(f'{v:.3f}' for v in MEAS['pole_residue'])}"),
]

def page3(c):
    c.page()
    c.rect(0, 0, W, H, fill=PAPER)
    eyebrow(c, "how the identity actually works", M, H - 50)
    c.text("Five steps, each one measured", M, H - 76, 'bodyb', 20.5, INK)
    c.text("Every ✓ below is a number produced by verify_string_copper.py, not a citation",
           M, H - 93, 'bodyi', 10.6, INK2)
    rule(c, M, W - M, H - 105, INK, 0.8)
    colw = 322.0; x2 = M + colw + 24
    yl = steps_column(c, M, H - 128, colw, "the chain, in order", STEPS_L, CU, 1)
    yr = steps_column(c, x2, H - 128, colw, "and its two closing moves", STEPS_R, STR, 4)
    assert min(yl, yr) > 72, f"page 3 steps end at {min(yl, yr):.1f}, into the folio"
    rule(c, M, W - M, 62)
    c.text("Page 3 of 7  ·  the mechanism", M, 50, 'sansm', 7.2, MUTED, 'l', 1.2)
    c.endpage()

# ------------------------------------------------------------------ page 4
CAVEAT_L = [
    ("One plane carries one spacetime dimension.",
     "The obstruction is RANK, not sign. $2α′p_i·p_j$ has rank up to D−1 — rank 3 at the "
     f"kinematic point $s = −3, t = −2.5, α′ = 1$, where its eigenvalues are {EIGS} — while "
     "a real scalar current can only produce $I_iI_j$, which is rank ONE. (Indefiniteness "
     "is not the reason, and is not even general: at $s = −2, t = −1.5$ the same matrix is "
     "positive semi-definite.) One copper plane therefore realises a single target "
     "dimension; D dimensions need D superposed planes, the timelike one carrying imaginary "
     "current. Reading a current as a momentum vector is the one step in the statement that "
     "cannot be taken literally."),
    ("The integral as written is still divergent.",
     "The measure is exactly $SL(2,ℝ)$-invariant — measured ratio "
     f"{num(1 + MEAS['sl2r'], 12)} under four Möbius maps — and it needs BOTH conditions: "
     "the exponent collected at puncture $i$ is $2α′p_i^2 − 2 − 2α′p_i·(Σp)$. Breaking the "
     f"mass shell alone gives {num(MEAS['sl2r_off'], 4)}; breaking momentum conservation "
     f"alone gives {num(MEAS['sl2r_noncons'], 4)}. Given $Σp = 0$ it is $α′p^2 = 1$ that "
     "does the work — that condition is exactly “the vertex operator has conformal weight "
     "one”. Either way $∫\\rm{d}z_1…\\rm{d}z_N$ carries the infinite volume of the gauge "
     "group, so three punctures must be held fixed and only N−3 integrated: an identity "
     "between integrands first, and between amplitudes only after."),
    ("The self-energy is a regulator, not the mass shell.",
     "The $i = j$ terms diverge on both sides: on the sheet the cut-off is a pad’s finite "
     "radius, on the worldsheet it is normal ordering. But they are not the same object. "
     "$κ Σ I_i^2$ is a multiplicative constant, per plane; $α′p^2 = 1$ is a CONDITION on the "
     "full D-dimensional momentum, and no single plane can enforce a constraint that sums "
     "over D of them. The two are related regularisations, not two names for one thing."),
]
CAVEAT_R = [
    ("Rim versus interior is open versus closed.",
     "The factor of two between a contact on the insulating rim and one in the interior — "
     f"measured as {MEAS['ratio']:.5f} — is the same factor of two between the open-string "
     "boundary propagator $−2α′\\,\\rm{ln}\\,|x − y|$ and the closed-string bulk one "
     "$−α′\\,\\rm{ln}\\,|z − w|$. Moving a contact off the edge into the copper is not a small "
     "perturbation of the model: it converts an open-string insertion into a closed-string one."),
    ("This is the worldsheet CFT, and not a least-resistance path.",
     "The two-dimensional theory here is the string’s own worldsheet, and the copper plane is "
     "that worldsheet — not spacetime. It is not the boundary CFT of AdS/CFT, where the dual "
     "theory lives on the boundary of the spacetime the gravity fills; conflating those is the "
     "most inviting misreading of the statement. Nor is the amplitude a least-resistance "
     "path — though not for the obvious reason. Rayleigh’s minimum is exactly what sits in "
     "the exponent: for a given set of contact currents the weight is $\\rm{exp}(−P/P_N)$ "
     "with $P$ the OHMIC dissipation, the fluctuations about it contributing only a "
     "source-independent factor that cancels. What the amplitude then does is integrate over "
     "the contact POSITIONS. Neither step is a minimisation over paths."),
    ("The two readings use different charges.",
     "The strings on page 2 come from Mandelstam’s map, which reads each current as a "
     "light-cone momentum $p^+$ — one component. The Koba–Nielsen weight reads the same "
     "current as the whole momentum vector. Two identifications, not one."),
]
EVIDENCE = [
    (f"Rim-contact pair potential, $R_s = 1$", num(RIM), f"$−R_s/π = {num(RIM_X)}$", pct(RIM, RIM_X)),
    ("Interior-contact pair potential", num(INT), f"$−R_s/2π = {num(INT_X)}$", pct(INT, INT_X)),
    ("Rim / interior ratio", f"{MEAS['ratio']:.5f}", "2 exactly", "the image charge"),
    (f"Power-law coefficient, {MEAS['power_n']} patterns", num(PWR),
     f"$2R_s/π = {num(PWR_X)}$", f"{pct(PWR, PWR_X)}, R² = {MEAS['power_R2']:.8f}"),
    ("Non-neutral ΔP per size doubling", num(NN, 4),
     f"$(Q^2R_s/π)\\rm{{ln}}2 = {num(NN_X, 4)}$", pct(NN, NN_X)),
    ("Gauge-fixed integral vs Euler B", f"{MEAS['veneziano_value']:.14f}",
     f"$B({MEAS['veneziano_a']:.0f}, {MEAS['veneziano_b']})"
     f" = {MEAS['veneziano_value']:.14f}$", sci(MEAS['veneziano'])),
    ("$SL(2,ℝ)$ measure ratio, on the mass shell", num(1 + MEAS['sl2r'], 12),
     "1 exactly", "four Möbius maps"),
    ("the same, pushed off the mass shell", num(MEAS['sl2r_off']), "≠ 1",
     "the symmetry is the shell"),
    ("Koba–Nielsen exponent eigenvalues", EIGS, "indefinite, rank 3",
     "a real current gives rank 1, ≥ 0"),
    (f"Washer resistance, $a = {MEAS['annulus'][-1]['a']}$, $b = {MEAS['annulus'][-1]['b']}$",
     f"{MEAS['annulus'][-1]['R']:.6f}",
     f"$(1/2π)\\rm{{ln}}(b/a) = {MEAS['annulus'][-1]['t']:.6f}$",
     f"{MEAS['annulus'][-1]['err_pct']:.2f}%, first order in h"),
]
REFS = ("Euler to Goldbach, 8 Jan 1730  ·  J. B. Johnson, Phys. Rev. 32 (1928) 97  ·  "
        "H. Nyquist, Phys. Rev. 32 (1928) 110  ·  G. Veneziano, Nuovo Cimento A57 (1968) 190  ·  "
        "Z. Koba & H. B. Nielsen, Nucl. Phys. B10 (1969) 633  ·  D. B. Fairlie & H. B. Nielsen, "
        "Nucl. Phys. B20 (1970) 637  ·  A. M. Polyakov, Phys. Lett. B103 (1981) 207")

def page4(c):
    c.page()
    c.rect(0, 0, W, H, fill=PAPER)
    eyebrow(c, "what follows from it, and what was measured", M, H - 50)
    c.text("Three consequences, and the evidence", M, H - 76, 'bodyb', 20.5, INK)
    c.text("Each consequence is precise rather than poetic; each number below is "
           "produced by verify_string_copper.py", M, H - 93, 'bodyi', 10.6, INK2)
    rule(c, M, W - M, H - 105, INK, 0.8)

    y = H - 128
    cw = (W - 2 * M - 40) / 3
    ends = []
    for k, (head, body, chk) in enumerate(CONSEQ):
        x = M + k * (cw + 20)
        yh = headline(c, head, x, y, cw, 9.5) + 12.4
        yy = c.flow(body, x, yh - 12.8, cw, 'body', 9.05, 12.3, INK2)
        c.text("✓", x, yy - 0.5, 'sans', 8.0, (0.16, 0.45, 0.30))
        dx = c.measure("✓  ", 'sans', 8.0)
        ends.append(c.flow(chk, x + dx, yy - 0.5, cw - dx, 'bodyi', 8.5, 11.5,
                           (0.16, 0.42, 0.29), just=False))

    # the table follows the band rather than being pinned to the foot, so the
    # slack falls in the bottom margin instead of as a hole mid-page
    ytab = min(ends) - 30
    rule(c, M, W - M, ytab + 18)
    eyebrow(c, "what was measured", M, ytab + 4)
    c.text("the sheet is a lattice of unit resistors solved exactly by DCT; the washer "
           "row is a conjugate-gradient solve", W - M, ytab + 4, 'bodyi', 8.6, MUTED, 'r')
    # measured content needs 150/66/106/111 pt, so give each column room for it
    cols = [M, M + 216, M + 322, M + 470]
    cwid = [200.0, 90.0, 130.0, W - M - (M + 470)]
    ytab -= 14
    for lab, cx in (("quantity", cols[0]), ("measured", cols[1]),
                    ("exact", cols[2]), ("agreement", cols[3])):
        c.text(lab.upper(), cx, ytab, 'sansb', 6.6, MUTED, 'l', 1.2)
    ytab -= 5
    rule(c, M, W - M, ytab, RULE, 0.5)
    ytab -= 11.2
    for k, row in enumerate(EVIDENCE):
        if k % 2 == 0:
            c.rect(M - 5, ytab - 3.2, W - 2 * M + 10, 12.8, fill=(0.972, 0.966, 0.953))
        for ci, (txt, role, col) in enumerate((
                (row[0], 'body', INK), (None, None, None),
                (row[2], 'body', INK2), (row[3], 'bodyi', MUTED))):
            if txt is None:
                continue
            n_ = len(c.wrap(txt, role, 8.5 if ci != 3 else 8.3, cwid[ci]))
            assert n_ == 1, f"page 4 cell {ci} wraps to {n_} lines: {txt!r}"
            c.flow(txt, cols[ci], ytab, cwid[ci], role, 8.5 if ci != 3 else 8.3,
                   11, col, just=False)
        c.math(row[1], cols[1], ytab, 8.5, (0.20, 0.30, 0.22))
        ytab -= 12.8
    last = ytab + 12.8
    rule(c, M, W - M, last - 20)
    c.flow(REFS, M, last - 32, W - 2 * M, 'body', 7.4, 9.8, MUTED, just=False)
    assert last - 54 > 34, f"page 4 runs to {last-54:.1f}"
    c.text("Page 4 of 7  ·  consequences and evidence", M, last - 54,
           'sansm', 7.2, MUTED, 'l', 1.2)
    c.endpage()


# ------------------------------------------------------------------ page 5
def page5(c):
    c.page()
    c.rect(0, 0, W, H, fill=PAPER)
    eyebrow(c, "where the analogy is exact, and where it stops", M, H - 50)
    c.text("Six things the one-line version hides", M, H - 76, 'bodyb', 20.5, INK)
    c.text("None of these breaks the identity; each one bounds what it may be read to mean",
           M, H - 93, 'bodyi', 10.6, INK2)
    rule(c, M, W - M, H - 105, INK, 0.8)
    colw = 322.0; x2 = M + colw + 24
    n = 0
    for x, items in ((M, CAVEAT_L), (x2, CAVEAT_R)):
        y = H - 128
        for head, body in items:
            n += 1
            numbered(c, n, x, y, CU if x == M else STR)
            y = headline(c, head, x + 18, y, colw - 18) + 12.5
            y = c.flow(body, x + 18, y - 12.8, colw - 18, 'body', 9.05, 12.3, INK2) - 12.0
        assert y > 70, f"page 5 caveat column ends at {y:.1f}, into the folio"
    rule(c, M, W - M, 62)
    c.text("Page 5 of 7  ·  the caveats", M, 50, 'sansm', 7.2, MUTED, 'l', 1.2)
    c.endpage()



# ------------------------------------------------------------------ page 6
def washer(c, cx, cy, ro, ri, mode):
    """An annular sheet. mode: 'resistor' | 'closed' | 'open'. All three are the
    SAME object, so the rims are drawn identically in all three."""
    c.circle(cx, cy, ro, fill=(0.878, 0.702, 0.549, 0.55))
    c.circle(cx, cy, ri, fill=PAPER)
    ang = lambda d, rr: (cx + rr * math.cos(math.radians(d)),
                         cy + rr * math.sin(math.radians(d)))
    if mode == 'resistor':
        for d in range(0, 360, 30):
            c.arrow(*ang(d, ri + 2.5), *ang(d, ro - 2.5), (CU[0], CU[1], CU[2], 0.9), 0.8, 4.0)
    elif mode == 'closed':
        for k, f in enumerate((0.14, 0.34, 0.54, 0.74, 0.94)):
            rr = ri + f * (ro - ri)
            bold = (k == 2)
            c.dash([] if bold else [1.6, 1.5])
            c.circle(cx, cy, rr, stroke=(STR[0], STR[1], STR[2], 0.95 if bold else 0.5),
                     lw=1.5 if bold else 0.6)
            c.dash([])
        c.arrow(*ang(62, ro + 4), *ang(62, ro + 17), (0.10, 0.10, 0.12), 0.9, 4.2)
        c.math("τ", *ang(62, ro + 23), 8.6, INK, 'c')
    else:
        for d in range(0, 360, 30):
            bold = (d == 30)
            p0, p1 = ang(d, ri), ang(d, ro)
            c.line(*p0, *p1, (STR[0], STR[1], STR[2], 0.95 if bold else 0.5),
                   1.5 if bold else 0.6)
            for e in (p0, p1):
                c.circle(e[0], e[1], 1.5 if bold else 1.0, fill=STR)
        c.arc(cx, cy, ro + 12, math.radians(44), math.radians(104),
              stroke=(0.10, 0.10, 0.12), lw=0.9)
        tip = ang(106, ro + 12)
        c.arrowhead(tip[0], tip[1], math.radians(196), 4.2, (0.10, 0.10, 0.12))
        c.math("τ", *ang(74, ro + 23), 8.6, INK, 'c')
    c.circle(cx, cy, ri, stroke=(0.24, 0.13, 0.05), lw=1.2)
    c.circle(cx, cy, ro, stroke=(0.24, 0.13, 0.05), lw=1.2)


OPEN_CLOSED = [
    ("the insertion sits", "on the 1-D rim", "in the 2-D interior"),
    ("so the string is", "open", "closed"),
    ("its propagator", "$\u22122α′\\,\\rm{ln}\\,|x \u2212 y|$",
                        "$\u2212α′\\,\\rm{ln}\\,|z \u2212 w|$"),
    ("lightest massless state", "spin 1 — a gauge boson", "spin 2 — the graviton"),
]

LEAD6 = (
    "Everything so far has been OPEN strings: contacts on the rim, insertions on the "
    "boundary. Their lightest massless state is spin 1 — a gauge boson, like a photon — so "
    "there is no gravity anywhere in it. Gravity is carried by CLOSED strings, whose lightest "
    "massless state is spin 2; and a massless spin-2 field has no choice but to couple to "
    "everything’s stress-energy, which is what gravity is. On a worldsheet, open versus "
    "closed is literally rim versus interior. So there are exactly two ways to put gravity "
    "into a sheet of copper, and both are changes to its GEOMETRY rather than new ingredients.")


def page6(c):
    c.page()
    c.rect(0, 0, W, H, fill=PAPER)
    eyebrow(c, "where gravity enters", M, H - 50)
    c.text("Punch a hole in the copper", M, H - 76, 'bodyb', 20.5, INK)
    c.text("Two routes, and both of them are geometry", M, H - 93, 'bodyi', 10.6, INK2)
    rule(c, M, W - M, H - 105, INK, 0.8)
    ylead = c.flow(LEAD6, M, H - 126, W - 2 * M, 'body', 9.6, 13.0, INK)
    rule(c, M, W - M, ylead + 3, RULE, 0.6)

    colw = 322.0; x2 = M + colw + 24
    ytop = ylead - 13

    eyebrow(c, "route one — move a contact off the rim", M, ytop, CU)
    y = ytop - 19
    cw = (92.0, 112.0, 116.0)
    xs = (M, M + cw[0], M + cw[0] + cw[1])
    c.text("A CONTACT ON THE RIM", xs[1], y - 11, 'sansb', 6.4, MUTED, 'l', 1.1)
    c.text("ONE IN THE INTERIOR", xs[2], y - 11, 'sansb', 6.4, MUTED, 'l', 1.1)
    y -= 16
    rule(c, M, M + sum(cw), y, RULE, 0.5)
    y -= 11
    for k, (lab, a_, b_) in enumerate(OPEN_CLOSED):
        if k % 2 == 0:
            c.rect(M - 4, y - 3.4, sum(cw) + 8, 12.6, fill=(0.972, 0.966, 0.953))
        c.text(lab, xs[0], y, 'body', 8.5, INK2)
        c.flow(a_, xs[1], y, cw[1] - 6, 'body', 8.5, 11, INK, just=False)
        c.flow(b_, xs[2], y, cw[2] - 6, 'bodyb', 8.5, 11, INK, just=False)
        y -= 12.6
    yL = c.flow(
        f"The factor of two measured between those two propagators — {MEAS['ratio']:.5f} — is "
        "that divide showing up in a lattice Green’s function. But a current injection is the "
        "WRONG KIND of insertion. $\\rm{e}^{ip·X}$ is a monopole, and the state it makes is "
        "the tachyon. A graviton operator $ε_{μν}∂X^μ∂X^ν\\rm{e}^{ip·X}$ carries that same "
        "monopole factor — it injects current like anything else — but multiplied by a factor "
        "bilinear in $∂X$, which is the current density. So route one says WHERE gravity "
        "lives, not how to source it.",
        M, y - 8, colw, 'body', 9.05, 12.3, INK2)

    eyebrow(c, "route two — punch a hole", x2, ytop, CU)
    y = c.flow(
        "Change the sheet’s TOPOLOGY instead. Copper with a hole has two rims, and in string "
        "perturbation theory each hole is one LOOP — so a washer is the one-loop diagram. Cut "
        "it open and roll it up and it is a cylinder of length $\\rm{ln}(b/a)$. It can be "
        "sliced two ways, exactly as page 2 sliced the disk, and neither is more correct:",
        x2, ytop - 19, colw, 'body', 9.05, 12.3, INK2)

    RO, RI = 34.0, 12.5
    cyc = y - 62
    cols3 = (x2 + 52, x2 + 161, x2 + 270)
    for cx_, mode in zip(cols3, ('resistor', 'closed', 'open')):
        washer(c, cx_, cyc, RO, RI, mode)
    caps = [("as a resistor", "radial flow, radius a out to b"),
            ("sliced by circles", "each slice a CLOSED string, rim to rim"),
            ("sliced radially", "each slice an OPEN string, once round")]
    ylab = cyc - RO - 22
    for cx_, (h_, t_) in zip(cols3, caps):
        c.text(h_, cx_, ylab, 'sansb', 6.9, INK, 'c', 1.0)
        c.flow(t_, cx_ - 53, ylab - 11, 106, 'bodyi', 7.9, 10.2, INK2, just=False)
    yR = c.flow(
        "Both describe the same washer — which is the whole point, and page 7 takes it up. "
        "Conformally that cylinder’s length is $\\rm{ln}(b/a) = 2πR/R_s$, so resistance is "
        "SEPARATION between the rims. Raise it and the closed string travels far, heavy "
        "states die off exponentially and only the lightest survive the trip; lower it "
        "toward $b/a → 1$ and the open-string loop takes over instead.",
        x2, ylab - 34, colw, 'body', 9.05, 12.3, INK2)
    assert min(yL, yR) > 74, f"page 6 columns end at {min(yL, yR):.1f}"
    rule(c, M, W - M, 62)
    c.text("Page 6 of 7  ·  where gravity enters", M, 50, 'sansm', 7.2, MUTED, 'l', 1.2)
    c.endpage()


# ------------------------------------------------------------------ page 7
LESSONS = [
    ("Gravity is not optional.",
     "One washer, two slicings: the open-string loop and the closed-string exchange are the "
     "SAME object, so you cannot have the first without the second. A theory of open strings "
     "— of gauge fields on a boundary — is not consistent without closed strings in its "
     "interior. On this reading gravity is not an extra force anyone put in; it is a "
     "consistency requirement of the gauge sector. The copper picture makes that visible: "
     "you cannot slice a washer only one way."),
    ("Gravity reads the flow’s anisotropy, not its strength.",
     "Across D superposed planes the graviton’s bilinear $∂X^μ∂X^ν$ is $J^μJ^ν$, a matrix of "
     "per-plane current densities, and it splits in two. The TRACE, $Σ_μ (J^μ)^2$, is the "
     "total Joule dissipation, and that part belongs to the DILATON. The TRACELESS remainder "
     "— how the dissipation is distributed among the planes — is the graviton. So gravity "
     "couples to the anisotropy of the flow rather than to how much heat is being made, and "
     "it couples wherever ANY current flows, whichever contact it came from. That is why "
     "gravity is universal. It also shows why a graviton needs at least two planes: a "
     "symmetric traceless tensor in one dimension is identically zero."),
    ("Its long range is the high-resistance end.",
     "The amplitude integrates over the washer’s shape, and the far end of that range — "
     "where the ratio and the resistance both run to infinity and the cylinder is long — is "
     "where only the lightest closed states survive the crossing. That limit is the "
     "long-distance force. Gravity’s reach is not a separate fact about gravity; it is which "
     "end of the shape integral you happen to stand at."),
]


def page7(c):
    c.page()
    c.rect(0, 0, W, H, fill=PAPER)
    eyebrow(c, "what the analogy teaches about gravity", M, H - 50)
    c.text("Three things it says gravity is", M, H - 76, 'bodyb', 20.5, INK)
    c.text("and one thing it cannot say", M, H - 93, 'bodyi', 10.6, INK2)
    rule(c, M, W - M, H - 105, INK, 0.8)

    cw = (W - 2 * M - 44) / 3
    ends = []
    for k, (head, body) in enumerate(LESSONS):
        x = M + k * (cw + 22)
        numbered(c, k + 1, x, H - 128, STR)
        y = headline(c, head, x + 18, H - 128, cw - 18) + 12.4
        ends.append(c.flow(body, x + 18, y - 12.8, cw - 18, 'body', 9.05, 12.3, INK2))

    y = min(ends) - 20
    colw = 322.0; x2 = M + colw + 24
    A = MEAS['annulus'][-1]
    yl = headline(c, "What a “modulus” is, and why it is here.", M, y, colw) + 12
    yl = c.flow(
        "A conformal map may stretch and bend a shape but never shear it. Any washer maps "
        "onto any other with the same ratio $b/a$ — and onto NONE with a different one — so "
        "a washer has exactly one shape parameter conformal maps cannot remove. That "
        "leftover number is its modulus, $t = (1/2π)\\rm{ln}(b/a)$, and an amplitude must "
        "integrate over it because it sums over every shape the surface can take: page 3's "
        f"gauge-fixing, one loop up. Measured, $a = {A['a']}$ and $b = {A['b']}$ give "
        f"$R = {A['R']:.6f}$ against $t = {A['t']:.6f}$, {A['err_pct']:.2f}% — so an "
        "ohmmeter on a washer reads a coordinate on the space of shapes.",
        M, yl - 12.8, colw, 'body', 9.05, 12.3, INK2)

    yr = headline(c, "What it cannot say.", x2, y, colw) + 12
    yr = c.flow(
        "This is the BOSONIC string, and that matters here. Its lightest closed state is not "
        "the graviton but the closed TACHYON, and the massless level holds a dilaton and a "
        "B-field beside the graviton — so the long-cylinder limit is not graviton exchange "
        "alone, nor even dominated by it. More basically: the copper is the WORLDSHEET, so "
        "gravity arrives as a string STATE, a particle that gets exchanged. The analogy "
        "gives the particle-physics face of gravity and nothing of the geometric one — no "
        "metric, no curvature, no Einstein equation.",
        x2, yr - 12.8, colw, 'body', 9.05, 12.3, INK2)

    yb = min(yl, yr) - 10
    rule(c, M, W - M, yb + 11)
    eyebrow(c, "reading all of this correctly", M, yb - 2)
    yb -= 20
    cw3 = (W - 2 * M - 40) / 3
    panes = [
        ("EXACT — AN IDENTITY", (0.16, 0.42, 0.29),
         "The integrand identity; the rim-versus-interior factor of two; $R/R_s$ as the "
         "annulus modulus. Theorems about one body of mathematics, measured here."),
        ("ANALOGY — A TRANSLATION", CU,
         "The page-2 dictionary translates; it does not identify. A current is one number, "
         "a momentum a spacetime vector — one plane, one dimension (p. 5)."),
        ("NOT CLAIMED", (0.46, 0.24, 0.28),
         "That a board gravitates beyond the pull of its mass; that measuring copper could "
         "test string theory; or that this is holography (p. 5)."),
    ]
    for k, (h_, col, t_) in enumerate(panes):
        x = M + k * (cw3 + 20)
        c.text(h_, x, yb, 'sansb', 7.0, col, 'l', 1.2)
        c.flow(t_, x, yb - 13, cw3, 'body', 8.6, 11.4, INK2)
    assert yb - 13 - 3 * 11.4 > 62, f"page 7 strip runs to {yb-13-3*11.4:.1f}"
    rule(c, M, W - M, 62)
    c.text("Page 7 of 7  ·  what it teaches", M, 50, 'sansm', 7.2, MUTED, 'l', 1.2)
    c.endpage()


# ------------------------------------------------------------------ main
def main():
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "string_copper_correspondence.pdf")
    bmp, n = field_bitmap()
    sl = streamlines(); eq = equipotentials(sl)
    slices = string_slices([1.25, 0.80, 0.50, 0.33, 0.20, 0.04, -0.16, -0.42, -0.85])
    print(f"  field     {n}×{n} px, {len(sl)} flow lines, {len(eq)} equipotentials, "
          f"{len(slices)} string slices")
    c = Canvas(out, W, H, {
        "kCGPDFContextTitle": "The Worldsheet and the Copper Plane",
        "kCGPDFContextAuthor": "Andy Haas",
        "kCGPDFContextCreator": "string_copper_correspondence.py",
        "kCGPDFContextSubject": "The Koba-Nielsen amplitude as Johnson-Nyquist power "
                                "fluctuations on a two-dimensional resistive sheet",
    })
    page1(c); page2(c, bmp, n, sl, eq, slices)
    page3(c); page4(c); page5(c); page6(c); page7(c)
    c.close()
    print(f"  wrote     {out}  ({os.path.getsize(out)/1024:.0f} kB)")

if __name__ == "__main__":
    main()
