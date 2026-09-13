/*
 * solver.js — electrostatics + magnetostatics of a deformable planar circuit.
 *
 * The wire is treated as a thin closed filament carrying an unknown line charge
 * density lambda(s). Ohm's law fixes the potential along the filament up to an
 * additive constant, which turns the problem into a boundary-value problem:
 *
 *     sum_i M_ji lambda_i  -  C  =  V_target(s_j)        for every segment j
 *     sum_i lambda_i L_i               =  0              (charge neutrality)
 *
 * Solve that for lambda and C, then E comes from Coulomb and B from Biot-Savart.
 * No field shape anywhere in this file is drawn by hand.
 *
 * Usable both as a browser global (window.FB) and as a CommonJS module.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.FB = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const EPS0 = 8.8541878128e-12;
  const KC = 1 / (4 * Math.PI * EPS0);
  const MU0 = 4e-7 * Math.PI;

  const clamp = (v, a, b) => (v < a ? a : v > b ? b : v);
  const lerp = (a, b, t) => a + (b - a) * t;

  /* ---------------------------------------------------------------
     geometry
     --------------------------------------------------------------- */

  // Closed Catmull-Rom through the control points. Returns the dense polyline
  // plus the dense index of each control point, so components can be anchored.
  function densePath(ctrl, per) {
    const n = ctrl.length, pts = [], anchor = [];
    for (let i = 0; i < n; i++) {
      const p0 = ctrl[(i - 1 + n) % n], p1 = ctrl[i],
            p2 = ctrl[(i + 1) % n], p3 = ctrl[(i + 2) % n];
      anchor.push(pts.length);
      for (let k = 0; k < per; k++) {
        const t = k / per, t2 = t * t, t3 = t2 * t;
        pts.push([
          0.5 * (2 * p1.x + (-p0.x + p2.x) * t + (2 * p0.x - 5 * p1.x + 4 * p2.x - p3.x) * t2 +
                 (-p0.x + 3 * p1.x - 3 * p2.x + p3.x) * t3),
          0.5 * (2 * p1.y + (-p0.y + p2.y) * t + (2 * p0.y - 5 * p1.y + 4 * p2.y - p3.y) * t2 +
                 (-p0.y + 3 * p1.y - 3 * p2.y + p3.y) * t3)
        ]);
      }
    }
    return { pts, anchor };
  }

  // Resample a closed polyline to N points of equal arc-length spacing.
  function resample(dense, N) {
    const M = dense.length, cum = new Float64Array(M + 1);
    for (let i = 1; i <= M; i++) {
      const p = dense[i - 1], q = dense[i % M];
      cum[i] = cum[i - 1] + Math.hypot(q[0] - p[0], q[1] - p[1]);
    }
    const L = cum[M], out = new Float64Array(N * 2);
    let j = 0;
    for (let i = 0; i < N; i++) {
      const s = L * i / N;
      while (j < M - 1 && cum[j + 1] < s) j++;
      const span = cum[j + 1] - cum[j];
      const f = span > 0 ? (s - cum[j]) / span : 0;
      const p = dense[j], q = dense[(j + 1) % M];
      out[i * 2] = p[0] + (q[0] - p[0]) * f;
      out[i * 2 + 1] = p[1] + (q[1] - p[1]) * f;
    }
    return { pts: out, L, cum };
  }

  /* ---------------------------------------------------------------
     dense linear solve (Gaussian elimination, partial pivoting)
     --------------------------------------------------------------- */
  function gauss(A, b, n) {
    for (let c = 0; c < n; c++) {
      let piv = c, best = Math.abs(A[c * n + c]);
      for (let r = c + 1; r < n; r++) {
        const v = Math.abs(A[r * n + c]);
        if (v > best) { best = v; piv = r; }
      }
      if (piv !== c) {
        for (let k = c; k < n; k++) {
          const t = A[c * n + k]; A[c * n + k] = A[piv * n + k]; A[piv * n + k] = t;
        }
        const t = b[c]; b[c] = b[piv]; b[piv] = t;
      }
      const d = A[c * n + c];
      if (Math.abs(d) < 1e-300) continue;
      for (let r = c + 1; r < n; r++) {
        const f = A[r * n + c] / d;
        if (f === 0) continue;
        for (let k = c; k < n; k++) A[r * n + k] -= f * A[c * n + k];
        b[r] -= f * b[c];
      }
    }
    const x = new Float64Array(n);
    for (let r = n - 1; r >= 0; r--) {
      let acc = b[r];
      for (let k = r + 1; k < n; k++) acc -= A[r * n + k] * x[k];
      x[r] = acc / (A[r * n + r] || 1e-300);
    }
    return x;
  }

  /* ---------------------------------------------------------------
     the circuit solve
     --------------------------------------------------------------- */
  const WIRE = 0, BATTERY = 1, LED = 2;

  /**
   * opts: { ctrl, N, per, radius, emf, rLED, rho }
   *   ctrl   array of {x,y} in metres
   *   N      number of filament segments
   *   per    spline samples per control-point interval
   *   radius wire radius in metres (also the regularisation length)
   *   emf    battery EMF in volts
   *   rLED   LED resistance in ohms
   *   rho    wire resistance per metre
   */
  function solveCircuit(opts) {
    const ctrl = opts.ctrl;
    const N = opts.N | 0, per = opts.per | 0, a = opts.radius;
    const a2 = a * a;

    const dense = densePath(ctrl, per);
    const rs = resample(dense.pts, N);
    const L = rs.L;
    if (!(L > 0.02)) return null;

    // components ride on control point 0 (battery) and the opposite one (LED)
    const midIdx = Math.floor(ctrl.length / 2);
    const sBat = rs.cum[Math.min(dense.anchor[0], rs.cum.length - 1)];
    const sLED = rs.cum[Math.min(dense.anchor[midIdx], rs.cum.length - 1)];
    const compLen = clamp(L * 0.09, 0.02, 0.09);

    // per-segment geometry, in flat typed arrays for fast field loops
    const cx = new Float64Array(N), cy = new Float64Array(N);
    const ax = new Float64Array(N), ay = new Float64Array(N);
    const bx = new Float64Array(N), by = new Float64Array(N);
    const tx = new Float64Array(N), ty = new Float64Array(N);
    const len = new Float64Array(N), sArr = new Float64Array(N);
    const type = new Uint8Array(N);

    let sAcc = 0;
    for (let i = 0; i < N; i++) {
      const px = rs.pts[i * 2], py = rs.pts[i * 2 + 1];
      const qx = rs.pts[((i + 1) % N) * 2], qy = rs.pts[((i + 1) % N) * 2 + 1];
      const dx = qx - px, dy = qy - py;
      const l = Math.hypot(dx, dy) || 1e-9;
      ax[i] = px; ay[i] = py; bx[i] = qx; by[i] = qy;
      cx[i] = (px + qx) / 2; cy[i] = (py + qy) / 2;
      tx[i] = dx / l; ty[i] = dy / l; len[i] = l;
      sArr[i] = sAcc + l / 2;
      sAcc += l;
    }

    const wrap = d => { d = ((d % L) + L) % L; return Math.min(d, L - d); };
    let batLen = 0, ledLen = 0, wireLen = 0;
    for (let i = 0; i < N; i++) {
      if (wrap(sArr[i] - sBat) < compLen / 2) { type[i] = BATTERY; batLen += len[i]; }
      else if (wrap(sArr[i] - sLED) < compLen / 2) { type[i] = LED; ledLen += len[i]; }
      else { type[i] = WIRE; wireLen += len[i]; }
    }
    if (batLen < 1e-6 || ledLen < 1e-6) return null;

    const Rwire = opts.rho * wireLen;
    const Rtot = opts.rLED + Rwire;
    const I = opts.emf / Rtot;              // conventional current, +s direction

    // prescribed potential gradient, and the resulting interior field E_t = -dV/ds
    const Et = new Float64Array(N);
    const dVds = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      if (type[i] === BATTERY) dVds[i] = opts.emf / batLen;
      else if (type[i] === LED) dVds[i] = -I * opts.rLED / ledLen;
      else dVds[i] = -I * opts.rho;
      Et[i] = -dVds[i];
    }
    const Vt = new Float64Array(N);
    let V = 0;
    for (let i = 0; i < N; i++) { Vt[i] = V + dVds[i] * len[i] / 2; V += dVds[i] * len[i]; }

    // assemble the (N+1) x (N+1) bordered system and solve
    const n = N + 1;
    const A = new Float64Array(n * n), rhs = new Float64Array(n);
    for (let j = 0; j < N; j++) {
      const jx = cx[j], jy = cy[j], base = j * n;
      for (let i = 0; i < N; i++) {
        if (i === j) {
          A[base + i] = KC * 2 * Math.asinh(len[i] / (2 * a));   // self-potential of a
        } else {                                                  // charged rod of radius a
          const dx = jx - cx[i], dy = jy - cy[i];
          A[base + i] = KC * len[i] / Math.sqrt(dx * dx + dy * dy + a2);
        }
      }
      A[base + N] = -1;
      rhs[j] = Vt[j];
    }
    for (let i = 0; i < N; i++) A[N * n + i] = len[i];
    rhs[N] = 0;

    const x = gauss(A, rhs, n);
    const lam = new Float64Array(N);
    let peak = 0;
    for (let i = 0; i < N; i++) { lam[i] = x[i]; if (Math.abs(lam[i]) > peak) peak = Math.abs(lam[i]); }

    return {
      N, L, a, I, Rwire, Rtot, batLen, ledLen, wireLen, compLen,
      V0: x[N], peak, closure: V,
      cx, cy, ax, ay, bx, by, tx, ty, len, s: sArr, type, lam, Et,
      outline: rs.pts,
      WIRE, BATTERY, LED
    };
  }

  /* ---------------------------------------------------------------
     field evaluation at an arbitrary point
     --------------------------------------------------------------- */

  // Inside the conductor the filament model says nothing useful, but Ohm's law
  // says everything: E = rho * J, tangential, which is exactly Et. Blend to it.
  // Inside the conductor B rises linearly from the axis, hence the d^2/a^2 mask.
  function fieldAt(sol, px, py) {
    const N = sol.N, a = sol.a, a2 = a * a;
    const { cx, cy, ax, ay, bx, by, tx, ty, len, lam } = sol;
    const kb = MU0 * sol.I / (4 * Math.PI);
    let ex = 0, ey = 0, bz = 0, v = 0, dmin = Infinity, near = 0;

    for (let i = 0; i < N; i++) {
      const dx = px - cx[i], dy = py - cy[i];
      const r2 = dx * dx + dy * dy;
      if (r2 < dmin) { dmin = r2; near = i; }

      const rr = r2 + a2, inv = 1 / Math.sqrt(rr);
      const kl = KC * lam[i] * len[i];
      const c = kl * inv / rr;
      ex += c * dx; ey += c * dy;
      v += kl * inv;

      if (r2 > 16 * len[i] * len[i]) {
        // far field: a point element is plenty
        const r3 = r2 * Math.sqrt(r2);
        bz += kb * len[i] * (tx[i] * dy - ty[i] * dx) / r3;
      } else {
        // near field: exact finite straight filament
        const Ax = ax[i] - px, Ay = ay[i] - py, Bx = bx[i] - px, By = by[i] - py;
        const la = Math.sqrt(Ax * Ax + Ay * Ay + 1e-16), lb = Math.sqrt(Bx * Bx + By * By + 1e-16);
        const cr = Ax * By - Ay * Bx, dt = Ax * Bx + Ay * By;
        const den = la * lb * (la * lb + dt);
        if (den > 1e-24) bz += kb * cr * (la + lb) / den;
      }
    }

    const d = Math.sqrt(dmin);
    if (d < a) bz *= (d * d) / a2;
    if (d < a * 1.35) {
      const w = clamp((a * 1.35 - d) / (a * 0.7), 0, 1);
      ex = lerp(ex, sol.Et[near] * tx[near], w);
      ey = lerp(ey, sol.Et[near] * ty[near], w);
    }
    return { ex, ey, bz, v, d, near };
  }

  /* ---------------------------------------------------------------
     grid sampling (the rendering hot path)
     --------------------------------------------------------------- */

  /**
   * opts: { w, h, cell, ppm, needE, needB }
   * w,h in pixels; the grid covers the rectangle centred on the world origin.
   */
  function buildGrid(sol, opts) {
    const w = opts.w, h = opts.h, ppm = opts.ppm;
    const needE = opts.needE !== false, needB = opts.needB !== false;
    const nx = Math.max(8, Math.floor(w / opts.cell));
    const ny = Math.max(8, Math.floor(h / opts.cell));
    const x0 = -w / 2 / ppm, y0 = -h / 2 / ppm;
    const dx = (w / ppm) / (nx - 1), dy = (h / ppm) / (ny - 1);

    const N = sol.N, a = sol.a, a2 = a * a;
    const { cx, cy, ax, ay, bx, by, tx, ty, len, lam, Et } = sol;
    const kb = MU0 * sol.I / (4 * Math.PI);

    const E = new Float32Array(nx * ny * 2);
    const B = new Float32Array(nx * ny);
    const D = new Float32Array(nx * ny);

    // precompute per-segment constants hoisted out of the inner loop
    const kl = new Float64Array(N), l2 = new Float64Array(N);
    for (let i = 0; i < N; i++) { kl[i] = KC * lam[i] * len[i]; l2[i] = 16 * len[i] * len[i]; }

    for (let j = 0; j < ny; j++) {
      const py = y0 + j * dy;
      for (let i = 0; i < nx; i++) {
        const px = x0 + i * dx;
        let ex = 0, ey = 0, bz = 0, dmin = Infinity, near = 0;

        for (let k = 0; k < N; k++) {
          const ddx = px - cx[k], ddy = py - cy[k];
          const r2 = ddx * ddx + ddy * ddy;
          if (r2 < dmin) { dmin = r2; near = k; }

          if (needE) {
            const rr = r2 + a2;
            const c = kl[k] / (rr * Math.sqrt(rr));
            ex += c * ddx; ey += c * ddy;
          }
          if (needB) {
            if (r2 > l2[k]) {
              bz += kb * len[k] * (tx[k] * ddy - ty[k] * ddx) / (r2 * Math.sqrt(r2));
            } else {
              const Ax = ax[k] - px, Ay = ay[k] - py, Bx = bx[k] - px, By = by[k] - py;
              const la = Math.sqrt(Ax * Ax + Ay * Ay + 1e-16), lb = Math.sqrt(Bx * Bx + By * By + 1e-16);
              const cr = Ax * By - Ay * Bx, dt = Ax * Bx + Ay * By;
              const den = la * lb * (la * lb + dt);
              if (den > 1e-24) bz += kb * cr * (la + lb) / den;
            }
          }
        }

        const d = Math.sqrt(dmin);
        if (needB && d < a) bz *= (d * d) / a2;
        if (needE && d < a * 1.35) {
          const wgt = clamp((a * 1.35 - d) / (a * 0.7), 0, 1);
          ex = lerp(ex, Et[near] * tx[near], wgt);
          ey = lerp(ey, Et[near] * ty[near], wgt);
        }

        const idx = j * nx + i;
        E[idx * 2] = ex; E[idx * 2 + 1] = ey; B[idx] = bz; D[idx] = d;
      }
    }
    return { nx, ny, x0, y0, dx, dy, E, B, D, w, h, ppm };
  }

  // Bilinear sample of the displayed vector field. mode 'E' or 'S'.
  function sampleVec(G, mode, x, y, out) {
    const fi = (x - G.x0) / G.dx, fj = (y - G.y0) / G.dy;
    if (fi < 0 || fj < 0 || fi > G.nx - 1.001 || fj > G.ny - 1.001) { out[0] = 0; out[1] = 0; return false; }
    const i = fi | 0, j = fj | 0, u = fi - i, v = fj - j;
    let vx = 0, vy = 0;
    for (let dj = 0; dj < 2; dj++) for (let di = 0; di < 2; di++) {
      const wgt = (di ? u : 1 - u) * (dj ? v : 1 - v), k = (j + dj) * G.nx + (i + di);
      if (mode === 'S') { const b = G.B[k]; vx += wgt * G.E[k * 2 + 1] * b / MU0; vy += wgt * -G.E[k * 2] * b / MU0; }
      else { vx += wgt * G.E[k * 2]; vy += wgt * G.E[k * 2 + 1]; }
    }
    out[0] = vx; out[1] = vy;
    return true;
  }

  function magAt(G, mode, k) {
    if (mode === 'B') return G.B[k];
    if (mode === 'S') { const b = G.B[k]; return Math.hypot(G.E[k * 2 + 1] * b, G.E[k * 2] * b) / MU0; }
    return Math.hypot(G.E[k * 2], G.E[k * 2 + 1]);
  }

  // S = E x B / mu0, in the plane, given a point sample.
  function poynting(f) {
    return [f.ey * f.bz / MU0, -f.ex * f.bz / MU0];
  }

  const PRESETS = {
    rect: [[-.26, -.15], [0, -.15], [.26, -.15], [.26, 0], [.26, .15], [0, .15], [-.26, .15], [-.26, 0]],
    circle: (() => { const a = []; for (let i = 0; i < 10; i++) { const t = i / 10 * 2 * Math.PI; a.push([.2 * Math.cos(t), .2 * Math.sin(t)]); } return a; })(),
    long: [[-.34, -.07], [-.1, -.07], [.12, -.07], [.34, -.07], [.34, .07], [.12, .07], [-.1, .07], [-.34, .07]],
    kink: [[-.3, -.16], [-.05, -.02], [.12, -.19], [.32, -.05], [.3, .14], [.05, .03], [-.12, .2], [-.32, .08]],
    folded: [[-.3, -.05], [-.05, -.16], [.28, -.1], [.3, -.02], [.05, 0], [-.05, .02], [.3, .04], [.28, .13]]
  };

  return {
    EPS0, KC, MU0, WIRE, BATTERY, LED, PRESETS,
    clamp, lerp, densePath, resample, gauss,
    solveCircuit, fieldAt, buildGrid, sampleVec, magAt, poynting
  };
});
