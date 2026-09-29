"""Solve a toy of the whole route's geometry LP (three lanes on nine column rails,
pinned at both ends, lifted over a pad box) and write its primal positions and
dual forces as a TikZ picture, ../figs/lp_cartoon.tex.  Rows as in whole_geo:
pitch struts between same-layer neighbours, the pad box as a wall, length by
tangent cuts of sqrt(1+s^2)-1 (a string at unit tension), bends at a small weight.

    python3 gen/lp_cartoon.py
"""
import math
import os

import numpy as np
from scipy.optimize import linprog

K, N, P, B = 9, 3, 1.0, 1.45          # columns, lanes, pitch, box top (clearance included)
S = 1.35                              # the lanes' spacing at their pinned ends: room to spare there
BOX = range(3, 6)                     # the columns the box spans
W_BEND = 0.25
W_REF = 0.02                          # a weak pull of every bead toward its lane's straight reference
TANG = [0.25, 0.5, 1.0, 2.0]

nv = 0


def new():
    global nv
    nv += 1
    return nv - 1


o = [[new() for k in range(K)] for n in range(N)]
d = [[new() for k in range(K - 1)] for n in range(N)]      # length excess per segment
b = [[new() for k in range(K - 2)] for n in range(N)]      # |bend| per interior vertex
c = np.zeros(nv)
for n in range(N):
    for k in range(K - 1):
        c[d[n][k]] = 1.0
    for k in range(K - 2):
        c[b[n][k]] = W_BEND
    for k in range(1, K - 1):
        c[o[n][k]] = W_REF             # the reference is n*S, at or below every bead

A, rhs, tag = [], [], []


def row(coefs, r, t):
    a = np.zeros(nv)
    for i, v in coefs:
        a[i] += v
    A.append(a)
    rhs.append(r)
    tag.append(t)


for n in range(N):
    for k in range(K - 1):
        for t in TANG:
            f, fp = math.sqrt(1 + t * t) - 1, t / math.sqrt(1 + t * t)
            for sg in (1, -1):        # d >= f(t) + f'(t) (s - t),  s = sg (o[k+1] - o[k])
                row([(d[n][k], -1), (o[n][k + 1], sg * fp), (o[n][k], -sg * fp)], fp * t - f, None)
    for k in range(1, K - 1):
        for sg in (1, -1):
            row([(b[n][k - 1], -1), (o[n][k + 1], sg), (o[n][k], -2 * sg), (o[n][k - 1], sg)], 0, None)
for k in range(K):                    # pitch struts: o[n+1] - o[n] >= P
    for n in range(N - 1):
        row([(o[n][k], 1), (o[n + 1][k], -1)], -P, ('strut', n, k))
for k in BOX:                         # the box as a wall under lane 0
    row([(o[0][k], -1)], -B, ('wall', 0, k))

bounds = [(None, None)] * nv
for n in range(N):
    bounds[o[n][0]] = bounds[o[n][K - 1]] = (n * S, n * S)    # pinned ends
    for k in range(K - 1):
        bounds[d[n][k]] = (0, None)
    for k in range(K - 2):
        bounds[b[n][k]] = (0, None)

res = linprog(c, A_ub=np.array(A), b_ub=np.array(rhs), bounds=bounds, method='highs')
assert res.status == 0, res.message
y = -res.ineqlin.marginals
X = [[res.x[o[n][k]] for k in range(K)] for n in range(N)]
forces = [(t, y[i]) for i, t in enumerate(tag) if t and y[i] > 1e-7]
print('cost', round(res.fun, 4))
print('offsets', [[round(v, 3) for v in r] for r in X])
print('forces', [(t, round(v, 3)) for t, v in forces])

L = [r'\begin{tikzpicture}[x=0.78cm,y=0.62cm,font=\scriptsize]']
for k in range(K):
    L.append(r'\draw[gray!45,dashed,line width=0.3pt] (%d,-0.55)--(%d,3.8);' % (k, k))
L.append(r'\fill[gray!35] (%.2f,-0.5) rectangle (%.2f,%.2f);' % (BOX[0] - 0.35, BOX[-1] + 0.35, B - 0.12))
L.append(r'\node[font=\tiny] at (%.2f,%.2f) {pad box};' % ((BOX[0] + BOX[-1]) / 2, (B - 0.62) / 2))
for n in range(N):
    L.append(r'\draw[Fcu!85!black,line width=1.1pt] '
             + '--'.join('(%d,%.3f)' % (k, X[n][k]) for k in range(K)) + ';')
fmax = max(v for _, v in forces)
for (kind, n, k), v in forces:
    w = 0.5 + 2.2 * v / fmax
    if kind == 'strut':
        L.append(r'\draw[magenta!80!black,line width=%.2fpt] (%d,%.3f)--(%d,%.3f);'
                 % (w, k, X[n][k] + 0.07, k, X[n + 1][k] - 0.07))
        L.append(r'\node[magenta!70!black,right,font=\tiny,inner sep=1pt] at (%d,%.3f) {%.2f};'
                 % (k, (X[n][k] + X[n + 1][k]) / 2, v))
    else:
        # the wall's push, drawn from inside the box up to the bead it holds
        L.append(r'\draw[orange!90!black,line width=%.2fpt,-{Stealth[length=4pt]}] (%d,%.3f)--(%d,%.3f);'
                 % (w, k, B - 0.85, k, X[0][k] - 0.09))
        side = 'right' if k < (BOX[0] + BOX[-1]) / 2 else 'left'   # labels on the box's inner side
        L.append(r'\node[orange!80!black,%s,font=\tiny,inner sep=1.5pt] at (%d,%.3f) {%.2f};'
                 % (side, k, B - 0.6, v))
k0 = 1                                # a strut with room to spare carries nothing
L.append(r'\draw[gray,dotted,line width=0.6pt] (%d,%.3f)--(%d,%.3f);' % (k0, X[0][k0] + 0.07, k0, X[1][k0] - 0.07))
L.append(r'\node[gray,right,font=\tiny,inner sep=1pt] at (%d,%.3f) {0};' % (k0, (X[0][k0] + X[1][k0]) / 2))
for n in range(N):
    for k in range(K):
        fill = 'black' if k in (0, K - 1) else 'white'
        L.append(r'\filldraw[fill=%s,line width=0.5pt] (%d,%.3f) circle (2.2pt);' % (fill, k, X[n][k]))
L.append(r'\node[above,font=\tiny] at (6,%.3f) {$o_{n,k}$};' % (X[2][6] + 0.1))
L.append(r'\node[below,font=\tiny,gray] at (0,-0.55) {pinned};')
L.append(r'\node[below,font=\tiny,gray] at (%d,-0.55) {pinned};' % (K - 1))
L.append(r'\end{tikzpicture}')
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figs', 'lp_cartoon.tex')
open(out, 'w').write('\n'.join(L) + '\n')
print('total force on the box', round(sum(v for (t, _, _), v in forces if t == 'wall'), 3))
