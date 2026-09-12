#!/usr/bin/env python3
"""Build string_copper_correspondence.pdf.

A three-page broadside on the Coulomb-gas identity between the open-string
Koba-Nielsen amplitude and Johnson-Nyquist power fluctuations on a 2-D
resistive sheet. Every number quoted on page 3 comes from
verify_string_copper.py; the figure's flow field is the real harmonic
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

_SUP = str.maketrans("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")

def num(x, d=6):
    """Fixed-point with a real minus sign (U+2212), not a hyphen."""
    return f"{x:.{d}f}".replace("-", "\u2212")

def sci(x, sig=1):
    e = int(math.floor(math.log10(abs(x))))
    man = f"{x / 10 ** e:.{sig}f}".rstrip("0").rstrip(".")
    return f"{man}\u00d710{str(e).translate(_SUP)}"

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

def trace(w0, rot=0.0, sign=1.0, h=0.005, nmax=5000, rmax=0.9965, stop=0.016):
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
        if len(pts) > 60 and abs(w - w0) < 1.2 * h: break      # closed loop
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
    add('int', ); add('m', "\\rm{e}^{−P}\\,\\rm{d}z_1⋯\\rm{d}z_N"); add('gap', 13 * s)
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
    "logarithmic Green's function, so the power dissipated by any Kirchhoff-balanced set "
    "of edge currents is a sum of pairwise logarithms of the contact separations; and "
    "Johnson–Nyquist noise weights each microscopic thermal state by the exponential of "
    "its dissipated heat. Integrating those thermal states by sliding the injection "
    "points continuously along the boundary therefore reproduces the quantum path "
    "integral of the interacting strings, governed by the mathematical equivalence:")

PROVENANCE = [
    ("1730", "Euler", "the integral that interpolates the factorial, to Goldbach, 8 January."),
    ("1928", "Johnson · Nyquist", "measure, then explain, resistor noise. Phys. Rev. 32, 97 & 110."),
    ("1968", "Veneziano", "the crossing-symmetric four-point amplitude. Nuovo Cimento A57, 190."),
    ("1969", "Koba · Nielsen", "generalise it to N particles — the product above. Nucl. Phys. B10, 633."),
    ("1970", "Fairlie · Nielsen", "the analogue model: the amplitude as a flow's dissipation. B20, 637."),
    ("1981", "Polyakov", "makes the worldsheet a 2-D statistical system. Phys. Lett. B103, 207."),
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
    for x, lab, sub, col in (
        (M + 168, "read as a copper plane",
         "the Boltzmann weight of the dissipated power", CU),
        (W - M - 168, "read as a string worldsheet",
         "the Koba–Nielsen factor of the N-point amplitude", STR)):
        c.text(lab.upper(), x, 232, 'sansm', 7.2, col, 'c', 1.3)
        c.text(sub, x, 219, 'bodyi', 9.3, INK2, 'c')

    # provenance — a named, published chain, not a new observation
    eyebrow(c, "not a new observation — the chain of results it rests on", M, 194)
    rule(c, M, W - M, 186, RULE, 0.6)
    for k, (yr, who, what) in enumerate(PROVENANCE):
        x = M + (k // 3) * (colw + 24)
        y = 166 - (k % 3) * 25
        c.text(yr, x, y, 'sansb', 9.2, INK, 'l', 0.4)
        wy = c.measure(yr, 'sansb', 9.2, 0.4)
        c.text(who, x + wy + 7, y, 'bodyb', 9.2, CU if k < 3 else STR)
        c.flow(what, x, y - 11.5, colw, 'body', 8.3, 11.0, INK2, just=False)
    rule(c, M, W - M, 80)
    c.text("Page 1 of 5  ·  the statement", M, 66, 'sansm', 7.2, MUTED, 'l', 1.2)
    c.text("Every number on pages 3 to 5 was measured, not cited",
           W - M, 66, 'sansm', 7.2, MUTED, 'r', 1.2)
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

def conformal_mesh():
    """Images of the upper-half-plane coordinate grid under the Cayley map
    w = (z-i)/(z+i). The worldsheet's conformal structure is all it has."""
    cay = lambda z: (z - 1j) / (z + 1j)
    fam = []
    for a in (0.10, 0.24, 0.52, 1.1, 2.4, 6.0):
        t = np.linspace(-80, 80, 1400)
        fam.append(cay(t + 1j * a))
    for b in (-7.0, -2.8, -1.2, -0.45, 0.0, 0.45, 1.2, 2.8, 7.0):
        t = np.concatenate([np.geomspace(0.0015, 80, 1400)])
        fam.append(cay(b + 1j * t))
    return [seg for f in fam for seg in runs_inside(f)]

def ribbon(c, cx, cy, r, ang, L=40.0, wd=9.0, incoming=False, label="p", idx=1):
    """An external open string: a semi-infinite strip, conformally shrunk to a
    single boundary point. Transverse bows are the string at successive times."""
    n = complex(math.cos(ang), math.sin(ang))
    t = complex(-n.imag, n.real)
    P = (cx + r * n.real, cy + r * n.imag)
    def at(d, u):
        return (P[0] + d * n.real + u * t.real, P[1] + d * n.imag + u * t.imag)
    A1, B1 = at(L, wd), at(L, -wd)
    cA1, cA2 = at(0.30 * L, 0.22 * wd), at(0.70 * L, 0.88 * wd)
    cB1, cB2 = at(0.30 * L, -0.22 * wd), at(0.70 * L, -0.88 * wd)
    c.bezpath([('m', *P), ('c', *cA1, *cA2, *A1), ('l', *B1),
               ('c', *cB2, *cB1, *P), ('z',)],
              fill=(0.905, 0.898, 0.925, 0.34), stroke=(STR[0], STR[1], STR[2], 0.92), lw=0.9)
    for k, s in enumerate((0.26, 0.46, 0.68, 0.92)):
        hw = s * wd
        pts = []
        for u in np.linspace(-1, 1, 26):
            bow = (0.17 * hw * math.cos(math.pi * u / 2)
                   + 0.13 * hw * math.sin(2.4 * math.pi * u))
            pts.append(at(s * L + bow, u * hw))
        bold = (k == 3)
        c.polyline(pts, stroke=(STR[0], STR[1], STR[2], 0.95 if bold else 0.72),
                   lw=1.15 if bold else 0.7)
        for e in (pts[0], pts[-1]):          # the open string's two endpoints
            c.circle(e[0], e[1], 1.35 if bold else 1.0, fill=STR)
    # the momentum arrow, in for a source and out for a sink
    a0, a1 = 0.34 * L, L + 10.0
    if incoming:
        c.arrow(*at(a1, 0), *at(a0, 0), STR, 0.95, 4.6)
    else:
        c.arrow(*at(a0, 0), *at(a1, 0), STR, 0.95, 4.6)
    lx, ly = at(L + 18.0, 0)
    c.math(f"{label}_{idx}", lx, ly - 3.2, 9.4, STR, 'c')

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
    c.text("+" if src else "−", lx, ly + 8.6, 'sansb', 8.4, CU if src else COOL, 'c')

def zlabel(c, cx, cy, r, ang, idx, col):
    n = complex(math.cos(ang), math.sin(ang))
    lx, ly = cx + (r - 14.0) * n.real, cy + (r - 14.0) * n.imag
    c.circle(lx, ly, 7.4, fill=(PAPER[0], PAPER[1], PAPER[2], 0.86))
    c.math(f"z_{idx}", lx, ly - 3.0, 8.6, col, 'c')

def string_legend(c, cx, ytop, w):
    """A pinched strip like the ones on the disk, with the string called out.
    `ytop` is the top of the whole block; it grows downward."""
    L, wd = 80.0, 12.0
    x0 = cx - w / 2 + 8
    yc = ytop - 21                                    # the strip's axis
    at = lambda d, u: (x0 + d, yc + u)
    # the sweep direction, ABOVE the strip so it cannot sit on the caption
    c.arrow(*at(0.16 * L, wd + 8.5), *at(L, wd + 8.5),
            (STR[0], STR[1], STR[2], 0.8), 0.7, 3.6)
    c.text("τ", *at(L * 0.56, wd + 12.0), 'mathi', 7.8, STR, 'c')
    c.bezpath([('m', *at(0, 0)),
               ('c', *at(0.30 * L, 0.22 * wd), *at(0.70 * L, 0.88 * wd), *at(L, wd)),
               ('l', *at(L, -wd)),
               ('c', *at(0.70 * L, -0.88 * wd), *at(0.30 * L, -0.22 * wd), *at(0, 0)),
               ('z',)],
              fill=(0.905, 0.898, 0.925, 0.34),
              stroke=(STR[0], STR[1], STR[2], 0.92), lw=0.9)
    for k, s in enumerate((0.26, 0.46, 0.68, 0.92)):
        hw = s * wd
        pts = []
        for u in np.linspace(-1, 1, 26):
            bow = (0.17 * hw * math.cos(math.pi * u / 2)
                   + 0.13 * hw * math.sin(2.4 * math.pi * u))
            pts.append(at(s * L + bow, u * hw))
        bold = (k == 3)
        c.polyline(pts, stroke=(STR[0], STR[1], STR[2], 0.95 if bold else 0.70),
                   lw=1.3 if bold else 0.72)
        for e in (pts[0], pts[-1]):
            c.circle(e[0], e[1], 1.5 if bold else 1.05, fill=STR)
    c.circle(*at(0, 0), 2.0, fill=STR)                # the puncture it shrinks to
    c.math("z_i", x0 - 10, yc - 3, 7.6, STR, 'r')
    c.text("each rung is the 1-D string at one instant", cx, ytop - 42, 'bodyi', 7.9, INK2, 'c')
    c.text("the strip is the 2-D worldsheet it sweeps", cx, ytop - 51.5, 'bodyi', 7.9, INK2, 'c')


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

def page2(c, bmp, n, sl, eq, mesh):
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

    # ---- left: conformal mesh, rim, punctures, string strips
    c.save(); c.clip_circle(CL, CY, R * 0.999)
    for seg in mesh:
        c.polyline([to_pt(z, CL, CY, R) for z in seg],
                   stroke=(0.30, 0.30, 0.36, 0.26), lw=0.4)
    c.restore()
    c.circle(CL, CY, R, stroke=(STR[0], STR[1], STR[2], 0.95), lw=1.5)
    for i, a in enumerate(ANG):
        ribbon(c, CL, CY, R, a, incoming=(CUR[i] > 0), idx=i + 1)
    for i, a in enumerate(ANG):
        p = to_pt(np.exp(1j * a), CL, CY, R)
        c.circle(p[0], p[1], 2.5, fill=STR)
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
    c.text("the worldsheet field, in spacetime", CL, 208, 'sansm', 6.8, MUTED, 'c', 1.1)
    c.text("the potential — exact, because ΣI = 0", CR, 208, 'sansm', 6.8, MUTED, 'c', 1.1)

    # ---- the centre column
    CX = 397.0
    string_legend(c, CX + 6, 462, 138)
    c.text("↔", CX, CY + 14, 'math', 30, INK, 'c')
    c.text("ONE INTEGRAND", CX, CY - 10, 'sansb', 7.4, INK, 'c', 1.5)
    c.text("TWO READINGS", CX, CY - 21, 'sansb', 7.4, INK, 'c', 1.5)
    rows = [("p_i", "I_i"), ("X^μ(z)", "V(z)"), ("α′", "R_s/πP_N"),
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

    rule(c, M, W - M, 82)
    note = ("Colour is one scalar field, drawn identically on both sides — copper positive, teal negative: at left "
            "the embedding carrying the sheet into spacetime, at right the electrostatic potential. The flow lines "
            "are the true current density $J = −∇V/R_s$, traced by RK4 through the exact disk solution.")
    nl = len(c.wrap(note, 'body', 8.8, W - 2 * M))
    assert nl <= 2, f"page-2 note wraps to {nl} lines and would hit the folio"
    c.flow(note, M, 70, W - 2 * M, 'body', 8.8, 11.8, INK2)
    c.text("Page 2 of 5  ·  the graphic", M, 42, 'sansm', 7.2, MUTED, 'l', 1.2)
    c.endpage()


# ------------------------------------------------------------------ page 3
def numbered(c, n, x, y, col):
    c.circle(x + 5.4, y + 3.0, 6.8, fill=(col[0], col[1], col[2], 0.13))
    c.text(str(n), x + 5.4, y, 'sansb', 8.2, col, 'c')

def steps_column(c, x, y, w, title, items, col, start=1):
    eyebrow(c, title, x, y); y -= 19
    for k, (head, body, chk) in enumerate(items, start):
        numbered(c, k, x, y, col)
        y = headline(c, head, x + 18, y, w - 18) + 12.5
        yy = c.flow(body, x + 18, y - 12.8, w - 18, 'body', 9.1, 12.4, INK2)
        if chk:
            c.text("✓", x + 18, yy + 1.0, 'sans', 8.0, (0.16, 0.45, 0.30))
            dx = c.measure("✓  ", 'sans', 8.0)
            yy = c.flow(chk, x + 18 + dx, yy + 1.0, w - 18 - dx, 'bodyi', 8.5,
                        11.6, (0.16, 0.42, 0.29), just=False)
        y = yy - 10.0
    return y

STEPS_L = [
    ("The string side is a product of powers.",
     "An open string's vertex operators sit on the one-dimensional boundary of the "
     "worldsheet, where the embedding has the propagator "
     "$⟨X^μ(x)X^ν(y)⟩ = −2α′η^{μν}\\,\\rm{ln}\\,|x − y|$. Contracting N tachyon operators "
     "$\\rm{e}^{ip·X}$ turns that logarithm into the Koba–Nielsen factor "
     "$∏|z_i − z_j|^{2α′p_i·p_j}$.",
     f"the product is exp of a logarithmic pair energy, to {sci(MEAS['kn_boltzmann'])}"),
    ("A sheet's rim Green's function is $−R_s/π$ times a log.",
     "A two-dimensional sheet has a logarithmic Green's function, and a contact ON an "
     "insulating rim sees twice the interior value, because its image charge sits on top "
     "of it. That factor of two is not bookkeeping: it is the same two that separates the "
     "open-string boundary propagator from the closed-string bulk one.",
     f"measured {num(RIM)} against $−R_s/π = {num(RIM_X)}$, {pct(RIM, RIM_X)}, "
     f"converging as the sheet grows"),
    ("Neutral rim currents dissipate a sum of logs.",
     "Feeding that into $P = Σ_i I_i V(z_i)$ gives "
     "$P = −(2R_s/π) Σ_{i<j} I_i I_j \\,\\rm{ln}\\,|z_i − z_j| + κ Σ_i I_i^2$. There is no "
     "factor of one half here — that belongs to the energy of assembling static charges, "
     "not to the power delivered by current sources, and inserting it halves $α′$.",
     f"fitted {num(PWR)} against $2R_s/π = {num(PWR_X)}$ over {MEAS['power_n']} random "
     f"neutral patterns, R² = {MEAS['power_R2']:.8f}"),
]
STEPS_R = [
    ("Nyquist supplies the exponential — and its scale.",
     "Thermal noise makes the contact currents a Gaussian whose quadratic form is exactly "
     "the dissipated power, so a fluctuation's weight is $\\rm{exp}(−P/P_N)$ with "
     "$P_N = 4k_BT/τ$ for an averaging window $τ$. Writing $\\rm{exp}(−P)$ by itself is "
     "dimensionally incomplete: a power is not an exponent. Substituting step 3 then gives "
     "the Koba–Nielsen factor outright, with $α′ = R_sτ/4πk_BT$ — the Regge slope as a "
     "sheet resistance in units of noise power, with the right units of one over current squared.",
     "Ohmic flow strictly minimises P: a divergence-free loop of strength ε costs exactly 4ε²"),
    ("Gauge-fix three contacts, then slide the rest.",
     "Fixing $z_1,z_3,z_4 = 0,1,∞$ and integrating the last over $(0,1)$ collapses the "
     "N = 4 product to $∫x^{−α(s)−1}(1−x)^{−α(t)−1}\\rm{d}x = B(−α(s),−α(t))$, the "
     "Veneziano amplitude, with $α(s) = 1+α′s$. Euler's factorial-interpolation integral "
     "is the string amplitude.",
     f"the integral equals Euler's B to {sci(MEAS['veneziano'])}"),
]
CONSEQ = [
    ("The three orderings are the three channels.",
     "Four contacts admit three distinct cyclic orders around the rim. Those are the three "
     "integration ranges of $x$, and they give $B(−α(s),−α(t))$, $B(−α(t),−α(u))$ and "
     "$B(−α(s),−α(u))$ — the s-t, t-u and s-u channels, whose sum is the crossing-symmetric total.",
     f"all three verified to {sci(MEAS['channels'])} or better"),
    ("“Counting continuously” is an identity.",
     "$B(k+1, n−k+1) = 1/[(n+1)C(n,k)]$ exactly, so the Beta function really does "
     "interpolate the binomial coefficient. The integral counts arrangements the way a "
     "binomial does, with a continuous index rather than a discrete one.",
     f"{MEAS['binomial_pairs']} integer pairs, {MEAS['binomial_bad']} mismatches"),
    ("The particle spectrum is two pads touching.",
     "The s-channel poles come from the endpoint $x → 0$, where two contacts merge and the "
     "dissipated power diverges. Convergence needs $2α′p_1·p_2 > −1$; the poles sit at "
     "$α′s = −1, 0, 1, 2, …$, which is the open-string spectrum $m^2 = (J−1)/α′$.",
     "B(a,b) → 1/a confirmed as a → 0"),
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
    y = min(yl, yr) - 2
    rule(c, M, W - M, y + 10)
    eyebrow(c, "three consequences that are precise rather than poetic", M, y - 4)
    y -= 24
    cw = (W - 2 * M - 40) / 3
    for k, (head, body, chk) in enumerate(CONSEQ):
        x = M + k * (cw + 20)
        yh = headline(c, head, x, y, cw, 9.3) + 12.3
        yy = c.flow(body, x, yh - 12.6, cw, 'body', 8.9, 12.1, INK2)
        c.text("✓", x, yy - 0.5, 'sans', 8.0, (0.16, 0.45, 0.30))
        dx = c.measure("✓  ", 'sans', 8.0)
        c.flow(chk, x + dx, yy - 0.5, cw - dx, 'bodyi', 8.4, 11.4,
               (0.16, 0.42, 0.29), just=False)
    rule(c, M, W - M, 62)
    c.text("Page 3 of 5  ·  the mechanism", M, 50, 'sansm', 7.2, MUTED, 'l', 1.2)
    c.endpage()

# ------------------------------------------------------------------ page 4
CAVEAT_L = [
    ("One plane carries one spacetime dimension.",
     "$p_i·p_j$ is a Minkowski inner product, so the matrix of exponents is indefinite — "
     f"for four on-shell tachyons its eigenvalues are {EIGS}. A real scalar current can "
     "only produce $I_iI_j$, which is rank one and positive semi-definite. One copper plane "
     "therefore realises a single target dimension; D dimensions need D superposed planes, "
     "the timelike one carrying imaginary current. Reading a current as a momentum vector is "
     "the one step in the statement that cannot be taken literally."),
    ("The integral as written is still divergent.",
     "The measure is exactly $SL(2,ℝ)$-invariant — measured ratio "
     f"{num(1 + MEAS['sl2r'], 12)} under four Möbius maps — and it is invariant precisely "
     "because $α′p^2 = 1$; pushed off the mass shell the ratio becomes "
     f"{num(MEAS['sl2r_off'], 4)}. So $∫\\rm{{d}}z_1…\\rm{{d}}z_N$ contains the infinite "
     "volume of the gauge group. Three punctures must be held fixed and only N−3 integrated, "
     "which makes the equation an identity between integrands first and between amplitudes "
     "only after."),
    ("The self-energy and the mass shell are one thing.",
     "The $i = j$ terms diverge. On the sheet that is cut off by the finite radius of a real "
     "pad; on the worldsheet it is removed by normal ordering, and what normal ordering leaves "
     "behind is $α′p^2 = 1$. Regulator and mass shell are the same object under two names."),
]
CAVEAT_R = [
    ("Rim versus interior is open versus closed.",
     "The factor of two between a contact on the insulating rim and one in the interior — "
     f"measured as {MEAS['ratio']:.5f} — is the same factor of two between the open-string "
     "boundary propagator $−2α′\\,\\rm{ln}\\,|x−y|$ and the closed-string bulk one "
     "$−α′\\,\\rm{ln}\\,|z−w|$. Moving a contact off the edge into the copper is not a small "
     "perturbation of the model: it converts an open-string insertion into a closed-string one."),
    ("This is the worldsheet CFT, not holography.",
     "The two-dimensional theory here is the string's own worldsheet, and the copper plane is "
     "that worldsheet — not spacetime. It is not the boundary CFT of AdS/CFT, where the dual "
     "theory lives on the boundary of the spacetime the gravity fills; conflating the two is "
     "the most inviting misreading of the statement."),
    ("Nor is it a least-resistance path.",
     "Rayleigh's principle picks out only the single most probable configuration; the amplitude "
     "sums over all of them, which is the entire content of the exponential. And only the "
     "tachyonic integrand is this clean — excited states and superstrings carry further factors."),
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
    eyebrow(c, "where the analogy is exact, and where it stops", M, H - 50)
    c.text("Five things the one-line version hides", M, H - 76, 'bodyb', 20.5, INK)
    c.text("None of these breaks the identity; each one bounds what it may be read to mean",
           M, H - 93, 'bodyi', 10.6, INK2)
    rule(c, M, W - M, H - 105, INK, 0.8)
    colw = 322.0; x2 = M + colw + 24
    ends = []
    for x, items in ((M, CAVEAT_L), (x2, CAVEAT_R)):
        y = H - 126
        for head, body in items:
            y = headline(c, head, x, y, colw) + 12.5
            y = c.flow(body, x, y - 12.8, colw, 'body', 9.05, 12.25, INK2) - 10.0
        ends.append(y)

    ytab = min(min(ends) - 16, 212)
    rule(c, M, W - M, ytab + 18)
    eyebrow(c, "what was measured", M, ytab + 4)
    c.text("independent of any citation — the sheet is a lattice of unit resistors, solved exactly",
           W - M, ytab + 4, 'bodyi', 8.6, MUTED, 'r')
    cols = [M, M + 302, M + 452, M + 610]
    ytab -= 14
    for lab, cx in (("quantity", cols[0]), ("measured", cols[1]),
                    ("exact", cols[2]), ("agreement", cols[3])):
        c.text(lab.upper(), cx, ytab, 'sansb', 6.6, MUTED, 'l', 1.2)
    ytab -= 5
    rule(c, M, W - M, ytab, RULE, 0.5)
    ytab -= 11.2
    for k, row in enumerate(EVIDENCE):
        if k % 2 == 0:
            c.rect(M - 5, ytab - 3.6, W - 2 * M + 10, 12.6, fill=(0.972, 0.966, 0.953))
        c.flow(row[0], cols[0], ytab, 300, 'body', 8.5, 11, INK, just=False)
        c.math(row[1], cols[1], ytab, 8.5, (0.20, 0.30, 0.22))
        c.flow(row[2], cols[2], ytab, 156, 'body', 8.5, 11, INK2, just=False)
        c.text(row[3], cols[3], ytab, 'bodyi', 8.3, MUTED)
        ytab -= 12.3
    last = ytab + 12.3            # the loop steps once past the final row
    assert last > 62, f"page 4 evidence table's last row sits at {last:.1f}, "\
                      f"which runs into the references at 58"
    rule(c, M, W - M, 58)
    c.flow(REFS, M, 47, W - 2 * M - 128, 'body', 7.4, 9.8, MUTED, just=False)
    c.text("Page 4 of 5  ·  the caveats", W - M, 47, 'sansm', 7.2, MUTED, 'r', 1.2)
    c.endpage()

# ------------------------------------------------------------------ page 5
def washer(c, cx, cy, ro, ri, mode):
    """An annular sheet. mode: 'resistor' | 'closed' | 'open'."""
    c.circle(cx, cy, ro, fill=(0.878, 0.702, 0.549, 0.55))
    c.circle(cx, cy, ri, fill=PAPER)
    ang = lambda d, rr: (cx + rr * math.cos(math.radians(d)),
                         cy + rr * math.sin(math.radians(d)))
    if mode == 'resistor':
        for d in range(0, 360, 30):
            c.arrow(*ang(d, ri + 2.5), *ang(d, ro - 2.5), (CU[0], CU[1], CU[2], 0.9), 0.8, 4.0)
        c.circle(cx, cy, ri, stroke=(0.24, 0.13, 0.05), lw=1.5)
        c.circle(cx, cy, ro, stroke=(0.24, 0.13, 0.05), lw=1.5)
        c.math("a", cx, cy - 3.0, 8.2, INK, 'c')
        c.line(*ang(158, ro + 3), *ang(158, ro + 9), INK2, 0.5)
        c.math("b", *ang(158, ro + 15), 8.2, INK, 'c')
    elif mode == 'closed':
        for k, f in enumerate((0.14, 0.34, 0.54, 0.74, 0.94)):
            rr = ri + f * (ro - ri)
            bold = (k == 2)
            c.dash([] if bold else [1.6, 1.5])
            c.circle(cx, cy, rr, stroke=(STR[0], STR[1], STR[2], 0.95 if bold else 0.5),
                     lw=1.5 if bold else 0.6)
            c.dash([])
        c.circle(cx, cy, ri, stroke=(0.24, 0.13, 0.05, 0.5), lw=0.8)
        c.circle(cx, cy, ro, stroke=(0.24, 0.13, 0.05, 0.5), lw=0.8)
        c.arrow(*ang(62, ri + 1.5), *ang(62, ro + 9), (0.10, 0.10, 0.12), 0.9, 4.2)
        c.math("τ", *ang(62, ro + 15), 8.6, INK, 'c')
    else:
        for k, d in enumerate(range(0, 360, 30)):
            bold = (d == 30)
            p0, p1 = ang(d, ri), ang(d, ro)
            c.line(*p0, *p1, (STR[0], STR[1], STR[2], 0.95 if bold else 0.5),
                   1.5 if bold else 0.6)
            for e in (p0, p1):
                c.circle(e[0], e[1], 1.5 if bold else 1.0, fill=STR)
        c.circle(cx, cy, ri, stroke=(0.24, 0.13, 0.05, 0.5), lw=0.8)
        c.circle(cx, cy, ro, stroke=(0.24, 0.13, 0.05, 0.5), lw=0.8)
        c.arc(cx, cy, ro + 9, math.radians(44), math.radians(104),
              stroke=(0.10, 0.10, 0.12), lw=0.9)
        tip = ang(106, ro + 9)
        c.arrowhead(tip[0], tip[1], math.radians(106 + 90), 4.2, (0.10, 0.10, 0.12))
        c.math("τ", *ang(124, ro + 14), 8.6, INK, 'c')


STRIP_TOP = 158.0

OPEN_CLOSED = [
    ("the insertion sits", "on the 1-D rim", "in the 2-D interior"),
    ("so the string is", "open", "closed"),
    ("its propagator", "$−2α′\\,\\rm{ln}\\,|x−y|$", "$−α′\\,\\rm{ln}\\,|z−w|$"),
    ("lightest massless state", "spin 1 — a gauge boson", "spin 2 — the graviton"),
]

def page5(c):
    c.page()
    c.rect(0, 0, W, H, fill=PAPER)
    eyebrow(c, "where gravity enters — and how far that goes", M, H - 50)
    c.text("Punch a hole in the copper", M, H - 76, 'bodyb', 20.5, INK)
    c.text("Closed strings carry the graviton, and closed means the interior, not the rim",
           M, H - 93, 'bodyi', 10.6, INK2)
    rule(c, M, W - M, H - 105, INK, 0.8)

    colw = 322.0; x2 = M + colw + 24
    # ---- left: the open/closed table
    y = H - 126
    y = headline(c, "Gravity is a closed-string state.", M, y, colw) + 12
    cw = (92.0, 112.0, 116.0)
    xs = (M, M + cw[0], M + cw[0] + cw[1])
    for lab, xx in (("", xs[0]), ("A CONTACT ON THE RIM", xs[1]),
                    ("ONE IN THE INTERIOR", xs[2])):
        if lab:
            c.text(lab, xx, y - 12, 'sansb', 6.4, MUTED, 'l', 1.1)
    y -= 17
    rule(c, M, M + sum(cw), y, RULE, 0.5)
    y -= 11
    for k, (lab, a_, b_) in enumerate(OPEN_CLOSED):
        if k % 2 == 0:
            c.rect(M - 4, y - 3.6, sum(cw) + 8, 12.8, fill=(0.972, 0.966, 0.953))
        c.text(lab, xs[0], y, 'body', 8.5, INK2)
        c.flow(a_, xs[1], y, cw[1] - 6, 'body', 8.5, 11, INK, just=False)
        c.flow(b_, xs[2], y, cw[2] - 6, 'bodyb', 8.5, 11, INK, just=False)
        y -= 12.8
    y -= 6
    y = c.flow(
        f"The factor of two I measured between those two propagators — "
        f"{MEAS['ratio']:.5f} — is therefore not a curiosity of lattice Green's "
        "functions. It is the gauge-theory/gravity divide. Open strings give the "
        "gauge fields that live on a boundary; closed strings give gravity in the "
        "bulk. Sliding a current contact off the board edge and into the copper "
        "crosses exactly that line.",
        M, y, colw, 'body', 9.05, 12.3, INK2) - 14

    y = headline(c, "But injecting a current is not a graviton.", M, y, colw) + 12
    yendL = c.flow(
        "The insertion this document analyses is $\\rm{e}^{ip·X}$: a point charge, a "
        "monopole current source — and the state it creates is the tachyon, not the "
        "graviton. A graviton vertex operator is $ε_{μν}∂X^μ∂X^ν\\rm{e}^{ip·X}$ — one derivative holomorphic, "
        "one antiholomorphic — which is BILINEAR in the worldsheet field. Since $X$ "
        "plays the part of the "
        "potential, $∂X$ is the current density — so a graviton insertion injects no "
        "net current at all. It couples to something quadratic in the local current, "
        "which in copper is the Joule dissipation density $R_s|J|^2$. That is the same "
        "structure as the worldsheet stress tensor, and it is the honest version of "
        "“the graviton couples to stress-energy”.",
        M, y - 12.8, colw, 'body', 9.05, 12.3, INK2)
    assert yendL > STRIP_TOP + 4, f"page 5 left column ends at {yendL:.1f}"

    # ---- right: the washer
    y = H - 126
    y = headline(c, "A washer is the one-loop diagram.", x2, y, colw) + 12
    y = c.flow(
        "Copper with a hole has two rims, and that shape is the annulus — the one-loop "
        "open-string diagram. It can be sliced two ways, neither more correct:",
        x2, y - 12.8, colw, 'body', 9.05, 12.3, INK2) - 10

    RO, RI = 32.0, 12.0
    cyc = y - RO - 4
    cols3 = (x2 + 52, x2 + 161, x2 + 270)
    for cx_, mode in zip(cols3, ('resistor', 'closed', 'open')):
        washer(c, cx_, cyc, RO, RI, mode)
    caps = [("as a resistor", "radial flow, from a out to b"),
            ("sliced by circles", "each slice is a closed string, propagating rim to rim"),
            ("sliced radially", "each slice is an open string, once around the loop")]
    ylab = cyc - RO - 13
    for cx_, (h_, t_) in zip(cols3, caps):
        c.text(h_, cx_, ylab, 'sansb', 6.9, INK, 'c', 1.0)
        c.flow(t_, cx_ - 51, ylab - 11, 102, 'bodyi', 7.9, 10.2, INK2, just=False)
    y = ylab - 37

    y = headline(c, "What a “modulus” is.", x2, y, colw) + 12
    A = MEAS['annulus'][-1]
    y = c.flow(
        "A conformal map may stretch and bend a shape, but never shear it. Any washer "
        "maps onto any other with the same ratio $b/a$ — so a washer has exactly one "
        "shape parameter that conformal maps cannot remove. That leftover number is "
        "called its modulus, and an amplitude must "
        "integrate over it, because it sums over every shape the surface can take. Here "
        f"the modulus is $t = (1/2π)\\rm{{ln}}(b/a)$; on a lattice of unit resistors "
        f"$a = {A['a']}$, $b = {A['b']}$ give $R = {A['R']:.6f}$ against "
        f"$t = {A['t']:.6f}$ — {A['err_pct']:.2f}%.",
        x2, y - 12.8, colw, 'body', 9.05, 12.3, INK2) - 6
    yend = c.flow(
        "Conformally the cylinder's length is $\\rm{ln}(b/a) = 2πR/R_s$: resistance is "
        "separation. So a HIGH-resistance washer is the long-cylinder limit, where only "
        "the lightest closed states survive — graviton exchange between the rims; as "
        "$b/a → 1$ the open-string loop takes over instead.",
        x2, y, colw, 'body', 9.05, 12.3, INK2)
    assert yend > STRIP_TOP + 4, (
        f"page 5 right column ends at {yend:.1f}, below the strip at {STRIP_TOP}")

    # ---- the interpretation strip: the point of the page
    yb = STRIP_TOP - 16
    rule(c, M, W - M, yb + 16)
    eyebrow(c, "reading this page correctly", M, yb + 2)
    yb -= 20
    cw3 = (W - 2 * M - 40) / 3
    panes = [
        ("EXACT — an identity", (0.16, 0.42, 0.29),
         "The integrand identity itself; the rim-versus-interior factor of two; and "
         "$R/R_s$ as the annulus modulus. Theorems about one body of mathematics — and "
         "measured here, not cited."),
        ("ANALOGY — a dictionary", CU,
         "A current is one number, a momentum a spacetime vector: one plane, one "
         "dimension (p. 4). The copper is the WORLDSHEET — this graviton is a state on "
         "an auxiliary surface, not curvature of the room."),
        ("NOT CLAIMED", (0.46, 0.24, 0.28),
         "That a board gravitates beyond the pull of its mass; that measuring copper "
         "could test string theory; or that this is holography. A resistor network is "
         "not a quantum theory with a Regge slope."),
    ]
    for k, (h_, col, t_) in enumerate(panes):
        x = M + k * (cw3 + 20)
        c.text(h_, x, yb, 'sansb', 7.0, col, 'l', 1.2)
        c.flow(t_, x, yb - 13, cw3, 'body', 8.6, 11.7, INK2)
    rule(c, M, W - M, 62)
    c.text("Page 5 of 5  ·  gravity", M, 50, 'sansm', 7.2, MUTED, 'l', 1.2)
    c.text("A washer's ohmmeter reading is a measure on moduli space",
           W - M, 50, 'sansm', 7.2, MUTED, 'r', 1.2)
    c.endpage()


# ------------------------------------------------------------------ main
def main():
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "string_copper_correspondence.pdf")
    bmp, n = field_bitmap()
    sl = streamlines(); eq = equipotentials(sl); mesh = conformal_mesh()
    print(f"  field     {n}×{n} px, {len(sl)} streamlines, {len(eq)} equipotentials, "
          f"{len(mesh)} mesh arcs")
    c = Canvas(out, W, H, {
        "kCGPDFContextTitle": "The Worldsheet and the Copper Plane",
        "kCGPDFContextAuthor": "Andy Haas",
        "kCGPDFContextCreator": "string_copper_correspondence.py",
        "kCGPDFContextSubject": "The Koba-Nielsen amplitude as Johnson-Nyquist power "
                                "fluctuations on a two-dimensional resistive sheet",
    })
    page1(c); page2(c, bmp, n, sl, eq, mesh); page3(c); page4(c); page5(c)
    c.close()
    print(f"  wrote     {out}  ({os.path.getsize(out)/1024:.0f} kB)")

if __name__ == "__main__":
    main()
