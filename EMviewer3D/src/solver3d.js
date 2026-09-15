/*
 * solver3d.js — electrostatics + magnetostatics of a deformable circuit in 3D.
 *
 * Same physics as the planar field bench, with the plane restriction lifted. The
 * wire is a thin closed filament in space carrying an unknown line charge density
 * lambda(s). Ohm's law fixes the potential along it up to an additive constant,
 * which is a boundary-value problem:
 *
 *     sum_i M_ji lambda_i  -  C  =  V_target(s_j)        for every segment j
 *     sum_i lambda_i L_i               =  0              (charge neutrality)
 *
 * The Coulomb kernel 1/r was always three-dimensional, so that solve carries over
 * unchanged; what changes is that B is now a full vector rather than a single
 * component, and the geometry is a space curve rather than a plane curve.
 *
 * Usable both as a browser global (window.FB3) and as a CommonJS module.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.FB3 = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const EPS0 = 8.8541878128e-12;
  const KC = 1 / (4 * Math.PI * EPS0);
  const MU0 = 4e-7 * Math.PI;

  const clamp = (v, a, b) => (v < a ? a : v > b ? b : v);
  const lerp = (a, b, t) => a + (b - a) * t;

  /* ---------------------------------------------------------------
     small vector helpers (plain objects; the hot loops use flat arrays)
     --------------------------------------------------------------- */
  const v3 = (x, y, z) => ({ x: x, y: y, z: z });
  const add = (a, b) => v3(a.x + b.x, a.y + b.y, a.z + b.z);
  const sub = (a, b) => v3(a.x - b.x, a.y - b.y, a.z - b.z);
  const scale = (a, s) => v3(a.x * s, a.y * s, a.z * s);
  const dot = (a, b) => a.x * b.x + a.y * b.y + a.z * b.z;
  const cross = (a, b) => v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
  const norm = a => Math.sqrt(dot(a, a));
  const unit = a => { const n = norm(a) || 1e-300; return scale(a, 1 / n); };

  /* ---------------------------------------------------------------
     geometry
     --------------------------------------------------------------- */

  // Closed Catmull-Rom through the control points, in space. Returns the dense
  // polyline plus the dense index of each control point, so components anchor.
  function densePath(ctrl, per) {
    const n = ctrl.length, pts = [], anchor = [];
    for (let i = 0; i < n; i++) {
      const p0 = ctrl[(i - 1 + n) % n], p1 = ctrl[i],
            p2 = ctrl[(i + 1) % n], p3 = ctrl[(i + 2) % n];
      anchor.push(pts.length);
      for (let k = 0; k < per; k++) {
        const t = k / per, t2 = t * t, t3 = t2 * t;
        const c = (a, b, c2, d) =>
          0.5 * (2 * b + (-a + c2) * t + (2 * a - 5 * b + 4 * c2 - d) * t2 +
                 (-a + 3 * b - 3 * c2 + d) * t3);
        pts.push([c(p0.x, p1.x, p2.x, p3.x), c(p0.y, p1.y, p2.y, p3.y), c(p0.z, p1.z, p2.z, p3.z)]);
      }
    }
    return { pts, anchor };
  }

  // Resample a closed space polyline to N points of equal arc-length spacing.
  function resample(dense, N) {
    const M = dense.length, cum = new Float64Array(M + 1);
    for (let i = 1; i <= M; i++) {
      const p = dense[i - 1], q = dense[i % M];
      cum[i] = cum[i - 1] + Math.hypot(q[0] - p[0], q[1] - p[1], q[2] - p[2]);
    }
    const L = cum[M], out = new Float64Array(N * 3);
    let j = 0;
    for (let i = 0; i < N; i++) {
      const s = L * i / N;
      while (j < M - 1 && cum[j + 1] < s) j++;
      const span = cum[j + 1] - cum[j];
      const f = span > 0 ? (s - cum[j]) / span : 0;
      const p = dense[j], q = dense[(j + 1) % M];
      out[i * 3]     = p[0] + (q[0] - p[0]) * f;
      out[i * 3 + 1] = p[1] + (q[1] - p[1]) * f;
      out[i * 3 + 2] = p[2] + (q[2] - p[2]) * f;
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
   *   ctrl   array of {x,y,z} in metres
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

    const midIdx = Math.floor(ctrl.length / 2);
    const sBat = rs.cum[Math.min(dense.anchor[0], rs.cum.length - 1)];
    const sLED = rs.cum[Math.min(dense.anchor[midIdx], rs.cum.length - 1)];
    const compLen = clamp(L * 0.09, 0.02, 0.09);

    const cx = new Float64Array(N), cy = new Float64Array(N), cz = new Float64Array(N);
    const ax = new Float64Array(N), ay = new Float64Array(N), az = new Float64Array(N);
    const bx = new Float64Array(N), by = new Float64Array(N), bz = new Float64Array(N);
    const tx = new Float64Array(N), ty = new Float64Array(N), tz = new Float64Array(N);
    const len = new Float64Array(N), sArr = new Float64Array(N);
    const type = new Uint8Array(N);

    let sAcc = 0;
    for (let i = 0; i < N; i++) {
      const px = rs.pts[i * 3], py = rs.pts[i * 3 + 1], pz = rs.pts[i * 3 + 2];
      const k = ((i + 1) % N) * 3;
      const qx = rs.pts[k], qy = rs.pts[k + 1], qz = rs.pts[k + 2];
      const dx = qx - px, dy = qy - py, dz = qz - pz;
      const l = Math.hypot(dx, dy, dz) || 1e-9;
      ax[i] = px; ay[i] = py; az[i] = pz;
      bx[i] = qx; by[i] = qy; bz[i] = qz;
      cx[i] = (px + qx) / 2; cy[i] = (py + qy) / 2; cz[i] = (pz + qz) / 2;
      tx[i] = dx / l; ty[i] = dy / l; tz[i] = dz / l; len[i] = l;
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
    const I = opts.emf / Rtot;

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

    const n = N + 1;
    const A = new Float64Array(n * n), rhs = new Float64Array(n);
    for (let j = 0; j < N; j++) {
      const jx = cx[j], jy = cy[j], jz = cz[j], base = j * n;
      for (let i = 0; i < N; i++) {
        if (i === j) {
          A[base + i] = KC * 2 * Math.asinh(len[i] / (2 * a));
        } else {
          const dx = jx - cx[i], dy = jy - cy[i], dz = jz - cz[i];
          A[base + i] = KC * len[i] / Math.sqrt(dx * dx + dy * dy + dz * dz + a2);
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

    // axis-aligned bounds, so the camera and the slice can frame themselves
    let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < N; i++) {
      if (cx[i] < lo[0]) lo[0] = cx[i]; if (cx[i] > hi[0]) hi[0] = cx[i];
      if (cy[i] < lo[1]) lo[1] = cy[i]; if (cy[i] > hi[1]) hi[1] = cy[i];
      if (cz[i] < lo[2]) lo[2] = cz[i]; if (cz[i] > hi[2]) hi[2] = cz[i];
    }

    return {
      N, L, a, I, Rwire, Rtot, batLen, ledLen, wireLen, compLen,
      V0: x[N], peak, closure: V, lo, hi,
      cx, cy, cz, ax, ay, az, bx, by, bz, tx, ty, tz,
      len, s: sArr, type, lam, Et,
      outline: rs.pts,
      WIRE, BATTERY, LED
    };
  }

  /* ---------------------------------------------------------------
     field evaluation at an arbitrary point in space
     --------------------------------------------------------------- */

  /*
   * E is Coulomb over the segments, softened at the wire radius. B is the exact
   * field of a finite straight filament,
   *
   *     B = (mu0 I / 4pi) (L x r1) [ L.r1/|r1| - L.r2/|r2| ] / |L x r1|^2
   *
   * with r1, r2 the vectors from the segment ends to the field point and L the
   * segment. Beyond a few segment lengths the point element is indistinguishable
   * and much cheaper, so it takes over there.
   *
   * Unlike E, B is not softened at the wire radius. The d < a mask below already
   * takes the interior down linearly from the axis, and softening on top of that
   * biases the field low well outside the metal -- 6% at four wire radii.
   */
  function fieldAt(sol, px, py, pz) {
    const N = sol.N, a = sol.a, a2 = a * a;
    const { cx, cy, cz, ax, ay, az, bx, by, bz, tx, ty, tz, len, lam } = sol;
    const kb = MU0 * sol.I / (4 * Math.PI);
    let ex = 0, ey = 0, ez = 0, Bx = 0, By = 0, Bz = 0, v = 0, dmin = Infinity, near = 0;

    for (let i = 0; i < N; i++) {
      const dx = px - cx[i], dy = py - cy[i], dz = pz - cz[i];
      const r2 = dx * dx + dy * dy + dz * dz;
      if (r2 < dmin) { dmin = r2; near = i; }

      const rr = r2 + a2, inv = 1 / Math.sqrt(rr);
      const kl = KC * lam[i] * len[i];
      const c = kl * inv / rr;
      ex += c * dx; ey += c * dy; ez += c * dz;
      v += kl * inv;

      if (r2 > 16 * len[i] * len[i]) {
        const r3 = r2 * Math.sqrt(r2);
        const k = kb * len[i] / r3;
        Bx += k * (ty[i] * dz - tz[i] * dy);
        By += k * (tz[i] * dx - tx[i] * dz);
        Bz += k * (tx[i] * dy - ty[i] * dx);
      } else {
        const Lx = bx[i] - ax[i], Ly = by[i] - ay[i], Lz = bz[i] - az[i];
        const r1x = px - ax[i], r1y = py - ay[i], r1z = pz - az[i];
        const r2x = px - bx[i], r2y = py - by[i], r2z = pz - bz[i];
        const gx = Ly * r1z - Lz * r1y, gy = Lz * r1x - Lx * r1z, gz = Lx * r1y - Ly * r1x;
        const den = gx * gx + gy * gy + gz * gz;
        if (den > 1e-40) {
          const n1 = Math.sqrt(r1x * r1x + r1y * r1y + r1z * r1z + 1e-16);
          const n2 = Math.sqrt(r2x * r2x + r2y * r2y + r2z * r2z + 1e-16);
          const f = ((Lx * r1x + Ly * r1y + Lz * r1z) / n1 -
                     (Lx * r2x + Ly * r2y + Lz * r2z) / n2) / den;
          Bx += kb * gx * f; By += kb * gy * f; Bz += kb * gz * f;
        }
      }
    }

    const d = Math.sqrt(dmin);
    if (d < a) { const m = (d * d) / a2; Bx *= m; By *= m; Bz *= m; }
    if (d < a * 1.35) {
      const w = clamp((a * 1.35 - d) / (a * 0.7), 0, 1);
      ex = lerp(ex, sol.Et[near] * tx[near], w);
      ey = lerp(ey, sol.Et[near] * ty[near], w);
      ez = lerp(ez, sol.Et[near] * tz[near], w);
    }
    return { ex, ey, ez, bx: Bx, by: By, bz: Bz, v, d, near };
  }

  // Poynting vector S = E x B / mu0 from a field sample.
  function poynting(f) {
    return {
      x: (f.ey * f.bz - f.ez * f.by) / MU0,
      y: (f.ez * f.bx - f.ex * f.bz) / MU0,
      z: (f.ex * f.by - f.ey * f.bx) / MU0
    };
  }

  // The vector a given mode draws, and its magnitude.
  function modeVec(f, mode) {
    if (mode === 'E') return { x: f.ex, y: f.ey, z: f.ez };
    if (mode === 'B') return { x: f.bx, y: f.by, z: f.bz };
    return poynting(f);
  }

  /* ---------------------------------------------------------------
     sampling a plane (the rendering hot path)
     --------------------------------------------------------------- */

  /**
   * opts: { o, u, v, half, n, mode }
   *   o     plane origin {x,y,z}
   *   u,v   orthonormal in-plane axes
   *   half  half-extent of the sampled square, metres
   *   n     samples per side
   * Returns the mode vector and its magnitude at every node, plus the distance
   * to the wire so the caller can cull samples inside the metal.
   */
  function buildSlice(sol, opts) {
    const n = opts.n | 0, half = opts.half, o = opts.o, u = opts.u, w = opts.v;
    const step = (2 * half) / (n - 1);
    const V = new Float32Array(n * n * 3), M = new Float32Array(n * n), D = new Float32Array(n * n);
    let mx = 0;
    for (let j = 0; j < n; j++) {
      const b = -half + j * step;
      for (let i = 0; i < n; i++) {
        const t = -half + i * step;
        const px = o.x + u.x * t + w.x * b;
        const py = o.y + u.y * t + w.y * b;
        const pz = o.z + u.z * t + w.z * b;
        const f = fieldAt(sol, px, py, pz);
        const q = modeVec(f, opts.mode);
        const k = (j * n + i);
        V[k * 3] = q.x; V[k * 3 + 1] = q.y; V[k * 3 + 2] = q.z;
        const m = Math.hypot(q.x, q.y, q.z);
        M[k] = m; D[k] = f.d;
        if (f.d > sol.a && m > mx) mx = m;
      }
    }
    return { n, half, o, u, v: w, step, V, M, D, max: mx, mode: opts.mode };
  }

  /*
   * Sampling a volume. The app uses this to set its colour range: taking the
   * range from geometry rather than from whichever field lines happen to be
   * traced keeps one colour meaning one field strength as the circuit is
   * dragged about.
   *
   * Points are nudged off the exact lattice by a deterministic jitter, which
   * keeps the samples from lining up with whatever symmetry the circuit has and
   * biasing the percentiles. It is a pure function of the lattice indices, so
   * the same circuit always gives the same range.
   */
  function jitter(i, j, k, s) {
    const v = Math.sin(i * 127.1 + j * 311.7 + k * 74.7 + s * 39.3) * 43758.5453;
    return v - Math.floor(v) - 0.5;
  }

  /**
   * opts: { half, n, mode }
   * n^3 samples spanning [-half, half] on each axis, centred on the origin.
   * Returns each sample's position, the mode vector there, its magnitude, and
   * the distance to the wire so the caller can drop samples inside the metal.
   */
  function buildVolume(sol, opts) {
    const n = opts.n | 0, half = opts.half;
    const step = (2 * half) / (n - 1), jit = step * 0.34;
    const T = n * n * n;
    const P = new Float32Array(T * 3), V = new Float32Array(T * 3);
    const M = new Float32Array(T), D = new Float32Array(T);
    let mx = 0;
    for (let k = 0; k < n; k++) {
      for (let j = 0; j < n; j++) {
        for (let i = 0; i < n; i++) {
          const id = (k * n + j) * n + i;
          const px = -half + i * step + jitter(i, j, k, 0) * jit;
          const py = -half + j * step + jitter(i, j, k, 1) * jit;
          const pz = -half + k * step + jitter(i, j, k, 2) * jit;
          const f = fieldAt(sol, px, py, pz);
          const q = modeVec(f, opts.mode);
          P[id*3] = px; P[id*3+1] = py; P[id*3+2] = pz;
          V[id*3] = q.x; V[id*3+1] = q.y; V[id*3+2] = q.z;
          const m = Math.hypot(q.x, q.y, q.z);
          M[id] = m; D[id] = f.d;
          if (f.d > sol.a && m > mx) mx = m;
        }
      }
    }
    return { n, half, step, P, V, M, D, max: mx, mode: opts.mode };
  }

  /* ---------------------------------------------------------------
     presets
     --------------------------------------------------------------- */
  const rect = [[-.26,-.15,0],[0,-.15,0],[.26,-.15,0],[.26,0,0],[.26,.15,0],[0,.15,0],[-.26,.15,0],[-.26,0,0]];
  const circle = (() => { const a = []; for (let i = 0; i < 10; i++) { const t = i / 10 * 2 * Math.PI; a.push([.2 * Math.cos(t), .2 * Math.sin(t), 0]); } return a; })();
  // A planar loop lifted out of the plane, to show that nothing was special
  // about z = 0: the same rectangle, rotated 40 degrees about the x axis.
  const tilt = (() => {
    const c = Math.cos(0.7), s = Math.sin(0.7);
    return rect.map(p => [p[0], p[1] * c, p[1] * s]);
  })();
  // Five turns of a solenoid, closed by a return leg up one side. The coil winds
  // about z rather than about a horizontal axis so that the default viewpoint
  // looks across it and it reads as a spring; wound about x it is seen close to
  // end-on, and the turns fold over each other into an unreadable tangle.
  const helix = (() => {
    const a = [], turns = 5, perTurn = 6, R = 0.1, n = turns * perTurn;
    for (let i = 0; i < n; i++) {
      const t = i / perTurn * 2 * Math.PI;
      a.push([R * Math.cos(t), R * Math.sin(t), -0.2 + 0.4 * i / (n - 1)]);
    }
    a.push([0.27, 0, 0.24], [0.31, 0, 0], [0.27, 0, -0.24]);
    return a;
  })();
  // A saddle: the rectangle with opposite corners pushed out of plane.
  const saddle = rect.map((p, i) => [p[0], p[1], 0.13 * Math.sin(i / rect.length * 4 * Math.PI)]);

  const PRESETS = { rect, circle, tilt, helix, saddle };

  /*
   * Seeded from the slice so the controls stay honest -- move the plane and you
   * reseed -- but integrated in three dimensions, so a line is free to leave the
   * plane immediately and generally does. B lines close into rings around the
   * wire; E lines run from positive surface charge to negative.
   */
  function trace(sol, mode, p0, sign, step, maxSteps, bound2) {
    const pts = [[p0.x, p0.y, p0.z, 0]];      // x, y, z, |F| at that point
    let x = p0.x, y = p0.y, z = p0.z, closed = false, run = 0;
    let lox = x, loy = y, loz = z, hix = x, hiy = y, hiz = z;
    const shut = (step * 0.9) * (step * 0.9);
    for (let i = 0; i < maxSteps; i++) {
      const f1 = fieldAt(sol, x, y, z);
      if (f1.d < sol.a * 1.1) break;              // ran into the metal
      const v1 = modeVec(f1, mode);
      const m1 = Math.hypot(v1.x, v1.y, v1.z);
      if (!(m1 > 0)) break;
      pts[pts.length - 1][3] = m1;             // we are standing on the last point
      const h = step * sign;
      const mx = x + v1.x / m1 * h * 0.5, my = y + v1.y / m1 * h * 0.5, mz = z + v1.z / m1 * h * 0.5;
      const f2 = fieldAt(sol, mx, my, mz);
      const v2 = modeVec(f2, mode);
      const m2 = Math.hypot(v2.x, v2.y, v2.z);
      if (!(m2 > 0)) break;
      x += v2.x / m2 * h; y += v2.y / m2 * h; z += v2.z / m2 * h;
      pts.push([x, y, z, m1]);
      // A ring that comes back to its seed is finished; tracing on would redraw it.
      if (i > 6) {
        const dx = x - p0.x, dy = y - p0.y, dz = z - p0.z;
        if (dx*dx + dy*dy + dz*dz < shut) { closed = true; break; }
      }
      // Near a curved wire a B line does not close: it winds helically around the
      // wire, drifting along it. Left to run it spends every step adding length
      // inside a small volume and paints a dense scribble. So stop a line once it
      // has covered more than a few times its own extent -- a full circular ring
      // scores 2.2 by this measure and is untouched, while the winders reach 5-6.
      run += step;
      lox = Math.min(lox, x); hix = Math.max(hix, x);
      loy = Math.min(loy, y); hiy = Math.max(hiy, y);
      loz = Math.min(loz, z); hiz = Math.max(hiz, z);
      if (i > 15) {
        const ex = hix - lox, ey = hiy - loy, ez = hiz - loz;
        if (run > 4 * Math.sqrt(ex*ex + ey*ey + ez*ez)) break;
      }
      if (x*x + y*y + z*z > bound2) break;
    }
    return { pts, closed };
  }

  /*
   * Seeds sit on a regular lattice through the VOLUME and are taken in that
   * order, so the lines fill the space instead of fanning out of one plane.
   *
   * Order is what keeps the spacing even. Seeding strongest-field-first is field
   * order rather than spatial order, and letting every point of a traced line
   * block further seeds lets one long B ring fence off an awkward region so the
   * next seed goes wherever the hole happens to be. The hash survives only as a
   * duplicate guard.
   */
  /*
   * Where a line is allowed to start.
   *
   * The lattice fills a VOLUME, and it has to: a line can be anywhere, and a
   * lattice through the cube is the only thing that finds one without knowing
   * where to look. But the region a loop encloses is a SHEET through that cube,
   * of no volume at all, so the lattice lands in it only by accident -- and
   * with GRID even, as it is here, there is no layer at z = 0 to land in at
   * all, so the flat presets drew no line inside the loop, ever. For E that is
   * the worst place to lose: the lines that run from the + surface charge
   * across to the - charge are the ones the bench is for.
   *
   * The wire is the one curve guaranteed to lie on the edge of that sheet, so
   * seed a ring of points a few radii off it and the sheet is covered by
   * construction -- with no notion of "the plane of the loop" needed anywhere,
   * which is as well, since the solenoid and saddle presets have not got one.
   * It is also where E lines begin, E being the one field here with ends.
   * Measured by the coverage audit in the tests, it helps every field and
   * costs under a tenth more lines, so it runs for all of them.
   */
  function seedPoints(sol, GRID, span) {
    const ring = [], cone = [], lattice = [];
    const off = sol.a * 5, stride = Math.max(1, Math.round(sol.N / 16));
    for (let i = 0; i < sol.N; i += stride) {
      const tx = sol.tx[i], ty = sol.ty[i], tz = sol.tz[i];
      // any unit vector across the wire; the branch only avoids the degenerate one
      let ux, uy, uz;
      if (Math.abs(tz) < 0.9) { ux = -ty; uy = tx; uz = 0; }
      else                    { ux = 0; uy = -tz; uz = ty; }
      const n = Math.hypot(ux, uy, uz) || 1; ux /= n; uy /= n; uz /= n;
      const vx = ty*uz - tz*uy, vy = tz*ux - tx*uz, vz = tx*uy - ty*ux;
      for (let k = 0; k < 4; k++) {
        const c = Math.cos(k * Math.PI / 2), w = Math.sin(k * Math.PI / 2);
        ring.push({ x: sol.cx[i] + off*(ux*c + vx*w),
                    y: sol.cy[i] + off*(uy*c + vy*w),
                    z: sol.cz[i] + off*(uz*c + vz*w) });
      }
    }
    // and on the sheet itself: the cone over the wire from the loop's own
    // centroid, which is a surface spanned by the loop whatever it is bent
    // into. The ring above is not enough on its own for E x B, which points
    // ALONG the wire where it is close to it -- energy running down the line --
    // so a line started there follows the wire out and never crosses the middle.
    let gx = 0, gy = 0, gz = 0;
    for (let i = 0; i < sol.N; i++) { gx += sol.cx[i]; gy += sol.cy[i]; gz += sol.cz[i]; }
    gx /= sol.N; gy /= sol.N; gz /= sol.N;
    for (let i = 0; i < sol.N; i += stride) for (let k = 1; k <= 3; k++) {
      const f = k / 4;
      cone.push({ x: sol.cx[i] + (gx - sol.cx[i])*f,
                  y: sol.cy[i] + (gy - sol.cy[i])*f,
                  z: sol.cz[i] + (gz - sol.cz[i])*f });
    }
    for (let k = 0; k < GRID; k++) for (let j = 0; j < GRID; j++) for (let i = 0; i < GRID; i++)
      lattice.push({ x: (2 * (i + 0.5) / GRID - 1) * span,
                     y: (2 * (j + 0.5) / GRID - 1) * span,
                     z: (2 * (k + 0.5) / GRID - 1) * span });
    return [ring, cone, lattice];
  }

  // Take the families round-robin, each in proportion to its own length, so
  // that whatever thins them thins all three alike rather than eating the first
  // one whole. Nothing caps the count on this bench, but the separation rule is
  // order-dependent -- first come, first served -- so the order still decides
  // which lines survive it.
  function roundRobin(fams) {
    const out = [], at = fams.map(() => 0);
    let total = 0;
    for (const f of fams) total += f.length;
    for (let n = 0; n < total; n++) {
      let pick = -1, behind = Infinity;
      for (let i = 0; i < fams.length; i++) {
        if (at[i] >= fams[i].length) continue;
        const frac = at[i] / fams[i].length;
        if (frac < behind) { behind = frac; pick = i; }
      }
      out.push(fams[pick][at[pick]++]);
    }
    return out;
  }

  function buildStreams(sol, opts) {
    const mode = opts.mode, R = opts.scene;
    const GRID = 4;            // seed lattice side; GRID^3 candidates
    const step = R * 0.022, maxSteps = 120, bound2 = (R * 2.3) * (R * 2.3);
    /* --- how close two lines may come ------------------------------------
       A field line is a flux tube. E has no divergence anywhere off the wire --
       all the charge is ON the wire -- which is checked in the tests, so if
       every line carries the same flux then the number crossing unit area goes
       as |F| and the spacing between them as |F|^(-1/2). That IS the statement
       "where the lines crowd, the field is strong", and it is the one piece of
       quantitative information the line PATTERN carries, as against the paths
       themselves.

       A separation rule at a fixed fraction of the scene destroys it, because
       it thins hardest exactly where the crowding is the physics. Measured on
       the settled rectangle, the fixed rule held the line density to |F|^0.76
       for E, |F|^0.49 for B and |F|^0.37 for E x B, where flux says |F|^1.

       The clamp is not a fudge. |F| runs away at the wire and to nothing far
       out, and an unclamped rule answers that with a hairball against the metal
       and an empty far field -- the whole line budget spent where the picture
       was already legible. A bit over a factor of two either side of the median
       covers about twenty in |F|, which is the range actually on screen. */
    const dup = R * 0.155;     // the scaled rule packs tighter than a fixed one
                               // at the same nominal spacing, so the base opens
                               // out to hold the line count
    const RMIN = 0.45, RMAX = 2.2;

    /* Candidates and the field at each, gathered first. That costs a field
       sweep per candidate -- well under a millisecond against the tens that
       tracing costs -- and buys the typical |F| over the scene, which is what
       the rule has to be measured against. */
    const cand = [];
    let ref = 0;
    {
      const mags = [];
      for (const fam of seedPoints(sol, GRID, R * 1.15)) {
        const keep = [];
        for (const p of fam) {
          const f = fieldAt(sol, p.x, p.y, p.z);
          if (f.d < sol.a * 2.5) continue;               // seed is in the metal
          const q = modeVec(f, mode);
          const m = Math.hypot(q.x, q.y, q.z);
          if (!(m > 0)) continue;                        // nothing to follow
          keep.push({ x: p.x, y: p.y, z: p.z, m });
          mags.push(m);
        }
        cand.push(keep);
      }
      if (!mags.length) return [];
      mags.sort((x, y) => x - y);
      ref = mags[mags.length >> 1];       // median: |F| spans decades at the wire
    }
    const sep = m => dup * clamp(Math.sqrt(ref/m), RMIN, RMAX);

    /* The radius varies from point to point now, so a set of occupied cell keys
       will not do it: the test has to be a real distance. Cells are the LARGEST
       radius the rule can ask for, which makes the 27 neighbours enough to cover
       any of them, and each cell keeps the line points that fell in it. */
    const cell = dup * RMAX;
    const grid = new Map();
    const ckey = (x, y, z) => Math.floor(x/cell) + ',' + Math.floor(y/cell) + ',' + Math.floor(z/cell);
    const drop = (x, y, z) => {
      const k = ckey(x, y, z);
      let arr = grid.get(k);
      if (!arr) grid.set(k, arr = []);
      arr.push(x, y, z);
    };
    const crowded = (x, y, z, r) => {
      const r2 = r*r;
      const i = Math.floor(x/cell), j = Math.floor(y/cell), k = Math.floor(z/cell);
      for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) for (let c = -1; c <= 1; c++) {
        const arr = grid.get((i+a) + ',' + (j+b) + ',' + (k+c));
        if (!arr) continue;
        for (let n = 0; n < arr.length; n += 3) {
          const dx = arr[n]-x, dy = arr[n+1]-y, dz = arr[n+2]-z;
          if (dx*dx + dy*dy + dz*dz < r2) return true;
        }
      }
      return false;
    };

    const out = [];
    for (const p of roundRobin(cand)) {
      if (crowded(p.x, p.y, p.z, sep(p.m))) continue;   // already on a line drawn
      const fwd = trace(sol, mode, p, +1, step, maxSteps, bound2);
      let line;
      if (fwd.closed) {
        line = fwd.pts;          // came back to the seed; tracing the other
      } else {                   // way would only redraw the same ring
        const back = trace(sol, mode, p, -1, step, maxSteps, bound2).pts;
        back.reverse(); back.pop();
        line = back.concat(fwd.pts);
      }
      if (line.length < 8) continue;
      for (const w of line) drop(w[0], w[1], w[2]);
      out.push(line);
    }
    return out;
  }

  return {
    EPS0, KC, MU0, WIRE, BATTERY, LED, PRESETS,
    clamp, lerp, v3, add, sub, scale, dot, cross, norm, unit,
    densePath, resample, gauss,
    solveCircuit, fieldAt, buildSlice, buildVolume, poynting, modeVec,
    trace, seedPoints, buildStreams
  };
});
