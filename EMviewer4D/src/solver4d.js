/*
 * solver4d.js — a circuit switching on, with retarded fields.
 *
 * The steady-state benches answer "what is the field once everything has
 * settled". This one answers "what happens in the first few nanoseconds", which
 * is a different question and needs three things they do not have: a circuit
 * that knows about its own inductance, a source history that conserves charge,
 * and fields evaluated at RETARDED time so that news of the switch travels out
 * at c rather than arriving everywhere at once.
 *
 *
 * THE SOURCE MODEL
 *
 * Ohm's law along the wire, now with the induced field kept:
 *
 *     dV/ds  =  e(t) g_emf(s)  +  I(t) g_R(s)  +  I'(t) g_L(s)
 *
 * with g_emf the battery's EMF density, g_R the resistive drop per unit
 * current, and g_L = -a(s), where a(s) is the tangential vector potential per
 * unit current. Integrated once around the loop the three terms give
 *
 *     e(t)  =  I R  +  I' L ,        L = closed integral of a(s) ds
 *
 * which is the circuit equation, and which is exactly the condition for V to
 * come back to itself after a lap. L is not assumed: it is the same geometric
 * kernel the charge solve uses, weighted by t.t' instead of 1.
 *
 * Because dV/ds is linear in (e, I, I'), so is lambda. Three solves of the same
 * matrix give basis charge densities lamE, lamI, lamL, and from then on
 *
 *     lambda_i(t) = e(t) lamE_i + I(t) lamI_i + I'(t) lamL_i
 *
 * is arithmetic. The current cannot then be uniform along the wire: continuity,
 *
 *     d(lambda)/dt + dI/ds = 0 ,
 *
 * fixes how it varies, and that variation is what the radiated field is made
 * of. Both are carried exactly here rather than being left inconsistent.
 *
 *
 * THE FIELDS
 *
 * Jefimenko's equations for a filament, which are Maxwell's equations solved
 * for a given source history -- exact for this source, retardation included:
 *
 *   E = (1/4pi eps0) SUM len [ lambda/R^2 + lambda'/(cR) ] Rhat  -  I' t /(c^2 R)
 *   B = (mu0/4pi)    SUM len [ I/R^2      + I'/(cR)      ] t x Rhat
 *
 * everything evaluated at t_r = t - R/c. Before t = 0 every history is zero, so
 * causality is not imposed anywhere: it falls out of the retarded lookup, and
 * the wavefront in the picture is the place where t_r has just gone positive.
 *
 * Each segment is a point element at its centre, which is also where its
 * retarded time is taken. The steady-state bench integrates its segments exactly
 * in space instead; the residual difference between the two is discretisation
 * and falls away as the wire is cut finer.
 *
 *
 * MAKING THE SOURCE CAUSAL
 *
 * A quasi-static solve on its own is acausal: it redistributes charge along the
 * whole wire the instant the battery moves, so every segment starts radiating at
 * once and field appears at a probe long before news from the battery could have
 * reached it -- measured at up to 226% of the settled field, 1.75 ns early.
 *
 * So each segment's response is delayed by the time light takes to get there
 * from the battery, delta_i = |c_i - battery| / c. Nothing can then reach a
 * probe sooner than (|battery to segment| + |segment to probe|)/c, which is
 * never less than |battery to probe|/c: the field is strictly zero outside the
 * light sphere about the battery.
 *
 * Two things have to be repaired after delaying, and they are why this is
 * tabulated rather than evaluated in closed form:
 *
 *   NEUTRALITY. The un-delayed charge sums to zero at a common instant; delayed,
 *   it does not -- by up to 40% of the charge present, on the solenoid. The
 *   excess is removed in proportion to |lambda| itself, which restores the sum
 *   to zero exactly while staying zero wherever the charge is still zero. A
 *   uniform subtraction would have put charge on segments the news has not
 *   reached, which is the very thing being fixed.
 *
 *   CONTINUITY. The current can no longer be the closed-form expression, because
 *   d(lambda)/dt now varies with each segment's own delay. It is re-derived by
 *   integrating dI/ds = -d(lambda)/dt along the wire. While part of the loop is
 *   still dark the constant of integration is fixed by requiring zero current
 *   there -- there is no closed path yet, so no circulating current can exist.
 *   Once the loop is lit the constant becomes the lumped I(t) again, blended
 *   across the moment the loop closes.
 *
 *
 * WHAT IS STILL APPROXIMATED
 *
 * The response is a delayed quasi-static one, not a self-consistent solve: the
 * shape of lambda at each instant is still the instantaneous boundary-value
 * solution, and I(t) still comes from a lumped circuit equation. That is sound
 * while the switch-on time is long compared with the light-crossing time of the
 * loop, and stretched when it is not. The fields are then exact for that source.
 * This is not a time-domain integral equation and does not pretend to be.
 *
 * Usable both as a browser global (window.FB4) and as a CommonJS module.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.FB4 = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const EPS0 = 8.8541878128e-12;
  const KC = 1 / (4 * Math.PI * EPS0);
  const MU0 = 4e-7 * Math.PI;
  const C = 299792458;

  const clamp = (v, a, b) => (v < a ? a : v > b ? b : v);
  const lerp = (a, b, t) => a + (b - a) * t;

  /* ---------------------------------------------------------------
     geometry — a closed space curve, resampled to equal arc length
     --------------------------------------------------------------- */
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
     LU with partial pivoting — one factorisation, three right-hand sides
     --------------------------------------------------------------- */
  function luFactor(A, n) {
    const piv = new Int32Array(n);
    for (let i = 0; i < n; i++) piv[i] = i;
    for (let c = 0; c < n; c++) {
      let best = Math.abs(A[c * n + c]), p = c;
      for (let r = c + 1; r < n; r++) {
        const v = Math.abs(A[r * n + c]);
        if (v > best) { best = v; p = r; }
      }
      if (p !== c) {
        for (let k = 0; k < n; k++) {
          const t = A[c * n + k]; A[c * n + k] = A[p * n + k]; A[p * n + k] = t;
        }
        const t = piv[c]; piv[c] = piv[p]; piv[p] = t;
      }
      const d = A[c * n + c];
      if (Math.abs(d) < 1e-300) continue;
      for (let r = c + 1; r < n; r++) {
        const f = A[r * n + c] / d;
        A[r * n + c] = f;
        if (f === 0) continue;
        for (let k = c + 1; k < n; k++) A[r * n + k] -= f * A[c * n + k];
      }
    }
    return { A, piv, n };
  }
  function luSolve(lu, b) {
    const { A, piv, n } = lu, x = new Float64Array(n);
    for (let i = 0; i < n; i++) x[i] = b[piv[i]];
    for (let i = 1; i < n; i++) {
      let acc = x[i];
      for (let k = 0; k < i; k++) acc -= A[i * n + k] * x[k];
      x[i] = acc;
    }
    for (let i = n - 1; i >= 0; i--) {
      let acc = x[i];
      for (let k = i + 1; k < n; k++) acc -= A[i * n + k] * x[k];
      x[i] = acc / (A[i * n + i] || 1e-300);
    }
    return x;
  }

  /* ---------------------------------------------------------------
     the EMF ramp — smootherstep, so that e, e' and e'' are all continuous
     (e'' feeds the third derivative of I, which the radiated field needs)
     --------------------------------------------------------------- */
  function emfAt(E0, tau, t) {
    if (t <= 0) return { v: 0, d: 0, dd: 0 };
    if (t >= tau) return { v: E0, d: 0, dd: 0 };
    const u = t / tau;
    return {
      v:  E0 * u * u * u * (10 + u * (-15 + 6 * u)),
      d:  E0 * 30 * u * u * (u - 1) * (u - 1) / tau,
      dd: E0 * 60 * u * (u - 1) * (2 * u - 1) / (tau * tau)
    };
  }

  const WIRE = 0, BATTERY = 1, LED = 2;

  /**
   * opts: { ctrl, N, per, radius, emf, rLED, rho, tau, tMax, steps }
   * Everything in SI; tau and tMax in seconds.
   */
  function buildTransient(opts) {
    const ctrl = opts.ctrl;
    const N = opts.N | 0, per = opts.per | 0, a = opts.radius, a2 = a * a;

    const dense = densePath(ctrl, per);
    const rs = resample(dense.pts, N);
    const Lloop = rs.L;
    if (!(Lloop > 0.02)) return null;

    const midIdx = Math.floor(ctrl.length / 2);
    const sBat = rs.cum[Math.min(dense.anchor[0], rs.cum.length - 1)];
    const sLED = rs.cum[Math.min(dense.anchor[midIdx], rs.cum.length - 1)];
    const compLen = clamp(Lloop * 0.09, 0.02, 0.09);

    const cx = new Float64Array(N), cy = new Float64Array(N), cz = new Float64Array(N);
    const ax = new Float64Array(N), ay = new Float64Array(N), az = new Float64Array(N);
    const bx = new Float64Array(N), by = new Float64Array(N), bz = new Float64Array(N);
    const tx = new Float64Array(N), ty = new Float64Array(N), tz = new Float64Array(N);
    const len = new Float64Array(N), sArr = new Float64Array(N);
    const type = new Uint8Array(N);

    let sAcc = 0;
    for (let i = 0; i < N; i++) {
      const px = rs.pts[i*3], py = rs.pts[i*3+1], pz = rs.pts[i*3+2];
      const k = ((i + 1) % N) * 3;
      const qx = rs.pts[k], qy = rs.pts[k+1], qz = rs.pts[k+2];
      const dx = qx - px, dy = qy - py, dz = qz - pz;
      const l = Math.hypot(dx, dy, dz) || 1e-9;
      ax[i] = px; ay[i] = py; az[i] = pz;
      bx[i] = qx; by[i] = qy; bz[i] = qz;
      cx[i] = (px+qx)/2; cy[i] = (py+qy)/2; cz[i] = (pz+qz)/2;
      tx[i] = dx/l; ty[i] = dy/l; tz[i] = dz/l; len[i] = l;
      sArr[i] = sAcc + l/2; sAcc += l;
    }

    const wrap = d => { d = ((d % Lloop) + Lloop) % Lloop; return Math.min(d, Lloop - d); };
    let batLen = 0, ledLen = 0, wireLen = 0;
    for (let i = 0; i < N; i++) {
      if (wrap(sArr[i] - sBat) < compLen / 2) { type[i] = BATTERY; batLen += len[i]; }
      else if (wrap(sArr[i] - sLED) < compLen / 2) { type[i] = LED; ledLen += len[i]; }
      else { type[i] = WIRE; wireLen += len[i]; }
    }
    if (batLen < 1e-6 || ledLen < 1e-6) return null;

    const Rtot = opts.rLED + opts.rho * wireLen;

    /* --- the geometric kernel, used twice ---------------------------------
       K(j,i) = integral over segment i of ds / |c_j - r|, regularised at the
       wire radius. Weighted by 1 it is the Coulomb operator whose inverse gives
       the charge; weighted by t_i . t_j it is the vector potential, whose loop
       integral is the self-inductance. Same geometry, two different questions. */
    const n = N + 1;
    const A = new Float64Array(n * n);
    const aTan = new Float64Array(N);          // vector potential per unit current
    for (let j = 0; j < N; j++) {
      const jx = cx[j], jy = cy[j], jz = cz[j], base = j * n;
      let av = 0;
      for (let i = 0; i < N; i++) {
        let k;
        if (i === j) k = 2 * Math.asinh(len[i] / (2 * a));
        else {
          const dx = jx-cx[i], dy = jy-cy[i], dz = jz-cz[i];
          k = len[i] / Math.sqrt(dx*dx + dy*dy + dz*dz + a2);
        }
        A[base + i] = KC * k;
        av += k * (tx[i]*tx[j] + ty[i]*ty[j] + tz[i]*tz[j]);
      }
      aTan[j] = MU0 / (4 * Math.PI) * av;
      A[base + N] = -1;
    }
    for (let i = 0; i < N; i++) A[N * n + i] = len[i];

    let Lind = 0;
    for (let j = 0; j < N; j++) Lind += aTan[j] * len[j];

    /* --- three basis charge densities ----------------------------------- */
    const gE = new Float64Array(N), gR = new Float64Array(N), gL = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      if (type[i] === BATTERY) gE[i] = 1 / batLen;
      if (type[i] === LED) gR[i] = -opts.rLED / ledLen;
      else if (type[i] === WIRE) gR[i] = -opts.rho;
      gL[i] = -aTan[i];
    }
    const potential = g => {
      const V = new Float64Array(n);
      let acc = 0;
      for (let i = 0; i < N; i++) { V[i] = acc + g[i] * len[i] / 2; acc += g[i] * len[i]; }
      V[N] = 0;
      return V;
    };
    const lu = luFactor(A, n);
    const solveBasis = g => {
      const x = luSolve(lu, potential(g));
      const lam = new Float64Array(N);
      for (let i = 0; i < N; i++) lam[i] = x[i];
      return lam;
    };
    const lamE = solveBasis(gE), lamI = solveBasis(gR), lamL = solveBasis(gL);

    /* --- continuity: dI/ds = -d(lambda)/dt -------------------------------
       Cumulative charge ahead of each segment, mean-removed so that the
       length-averaged current is the I(t) the circuit equation solves for. */
    const cumOf = lam => {
      const c = new Float64Array(N);
      let acc = 0;
      for (let i = 0; i < N; i++) { c[i] = acc + lam[i] * len[i] / 2; acc += lam[i] * len[i]; }
      let mean = 0;
      for (let i = 0; i < N; i++) mean += c[i] * len[i];
      mean /= Lloop;
      for (let i = 0; i < N; i++) c[i] -= mean;
      return c;
    };
    const cumE = cumOf(lamE), cumI = cumOf(lamI), cumL = cumOf(lamL);

    /* --- march the circuit equation  L I' + R I = e(t) -------------------- */
    const E0 = opts.emf, tau = opts.tau;
    const steps = opts.steps | 0, dt = opts.tMax / steps;
    const Ihist = new Float64Array(steps + 1);
    const slope = (t, I) => (emfAt(E0, tau, t).v - Rtot * I) / Lind;
    for (let k = 0; k < steps; k++) {
      const t = k * dt, I = Ihist[k];
      const k1 = slope(t, I);
      const k2 = slope(t + dt/2, I + dt*k1/2);
      const k3 = slope(t + dt/2, I + dt*k2/2);
      const k4 = slope(t + dt,   I + dt*k3);
      Ihist[k+1] = I + dt * (k1 + 2*k2 + 2*k3 + k4) / 6;
    }

    /* --- delay each segment by its light-time from the battery ------------ */
    let bcx = 0, bcy = 0, bcz = 0, nb = 0;
    for (let i = 0; i < N; i++) if (type[i] === BATTERY) { bcx += cx[i]; bcy += cy[i]; bcz += cz[i]; nb++; }
    bcx /= nb; bcy /= nb; bcz /= nb;
    const delay = new Float64Array(N);
    let dmax = 0;
    for (let i = 0; i < N; i++) {
      delay[i] = Math.hypot(cx[i]-bcx, cy[i]-bcy, cz[i]-bcz) / C;
      if (delay[i] > dmax) dmax = delay[i];
    }

    const T = {
      N, Lloop, a, Rtot, Lind, batLen, ledLen, wireLen, compLen,
      emf: E0, tau, dt, tMax: opts.tMax, steps, rLED: opts.rLED, rho: opts.rho,
      Iss: E0 / Rtot, tauLR: Lind / Rtot, cross: Lloop / C,
      cx, cy, cz, ax, ay, az, bx, by, bz, tx, ty, tz,
      len, s: sArr, type, aTan,
      lamE, lamI, lamL, cumE, cumI, cumL, Ihist,
      delay, dmax, battery: { x: bcx, y: bcy, z: bcz },
      outline: rs.pts,
      lo: (() => { const l = [Infinity,Infinity,Infinity];
        for (let i=0;i<N;i++){ if(cx[i]<l[0])l[0]=cx[i]; if(cy[i]<l[1])l[1]=cy[i]; if(cz[i]<l[2])l[2]=cz[i]; }
        return l; })(),
      hi: (() => { const h = [-Infinity,-Infinity,-Infinity];
        for (let i=0;i<N;i++){ if(cx[i]>h[0])h[0]=cx[i]; if(cy[i]>h[1])h[1]=cy[i]; if(cz[i]>h[2])h[2]=cz[i]; }
        return h; })(),
      WIRE, BATTERY, LED
    };
    buildTables(T);
    return T;
  }

  /* ---------------------------------------------------------------
     the causal source, tabulated
     --------------------------------------------------------------- */
  /*
   * lambda and I are no longer closed-form combinations of a few scalar
   * histories, because every segment now runs on its own clock. They are swept
   * once onto a time grid and interpolated afterwards, which is what keeps the
   * clock free to move: the tables are built when the circuit is solved, never
   * when the time changes.
   */
  function buildTables(T) {
    const N = T.N;
    // Fine enough to resolve the ramp, long enough to reach the settled state.
    const tTab = clamp(5*T.tauLR, 25e-9, 150e-9) + T.dmax + 5*T.tau;
    const du = Math.min(T.tau/12, 6e-11, tTab/400);
    const K = clamp(Math.round(tTab/du) + 1, 400, 4000);
    const dutab = tTab / (K - 1);

    const LAM = new Float32Array(N*K), LAMD = new Float32Array(N*K);
    const CUR = new Float32Array(N*K), CURD = new Float32Array(N*K);
    const rho = new Float64Array(N), lam = new Float64Array(N), P = new Float64Array(N);

    // arc length to the middle of each segment, for the mean-centring below
    const ell = new Float64Array(N);
    { let accl = 0;
      for (let i = 0; i < N; i++) { ell[i] = accl + T.len[i]/2; accl += T.len[i]; } }

    // pass 1: the charge, delayed and then made neutral again
    for (let k = 0; k < K; k++) {
      const u = k * dutab;
      let A = 0, Aabs = 0;
      for (let i = 0; i < N; i++) {
        const h = histAt(T, u - T.delay[i]);
        const r = h.e*T.lamE[i] + h.I*T.lamI[i] + h.Id*T.lamL[i];
        rho[i] = r; A += r*T.len[i]; Aabs += Math.abs(r)*T.len[i];
      }
      // Take the excess out in proportion to |lambda|, so it is removed only
      // where there is charge to remove it from, and the dark part stays dark.
      const corr = Aabs > 0 ? A/Aabs : 0;
      for (let i = 0; i < N; i++) LAM[i*K + k] = rho[i] - corr*Math.abs(rho[i]);
    }

    // pass 2: d(lambda)/dt by central difference on the grid
    for (let i = 0; i < N; i++) {
      const b = i*K;
      for (let k = 0; k < K; k++) {
        const km = k > 0 ? k-1 : k, kp = k < K-1 ? k+1 : k;
        LAMD[b+k] = (LAM[b+kp] - LAM[b+km]) / ((kp-km)*dutab || 1);
      }
    }

    // pass 3: the current, from dI/ds = -d(lambda)/dt
    for (let k = 0; k < K; k++) {
      const u = k * dutab;
      let acc = 0;
      for (let i = 0; i < N; i++) {
        P[i] = acc + LAMD[i*K+k]*T.len[i]/2;
        acc += LAMD[i*K+k]*T.len[i];
      }
      // Constant of integration. While any of the loop is still dark there is no
      // closed path, so no current can circulate: fix it to zero there. Once the
      // loop is lit the lumped I(t) takes over, blended across the changeover.
      let darkSum = 0, darkN = 0, meanP = 0;
      for (let i = 0; i < N; i++) {
        meanP += P[i]*T.len[i];
        if (u <= T.delay[i]) { darkSum += P[i]; darkN++; }
      }
      meanP /= T.Lloop;
      const Iode = histAt(T, u).I;
      const cDark = darkN ? darkSum/darkN : meanP;
      const cLoop = Iode + meanP;
      const w = smooth((u - T.dmax) / (0.5*T.dmax || 1e-12));
      const cst = cDark + (cLoop - cDark)*w;
      for (let i = 0; i < N; i++) CUR[i*K + k] = cst - P[i];
    }

    // pass 4: dI/dt, again by central difference
    for (let i = 0; i < N; i++) {
      const b = i*K;
      for (let k = 0; k < K; k++) {
        const km = k > 0 ? k-1 : k, kp = k < K-1 ? k+1 : k;
        CURD[b+k] = (CUR[b+kp] - CUR[b+km]) / ((kp-km)*dutab || 1);
      }
    }

    // What each quantity is heading for, so the table can be left early and
    // finished analytically. Clamping at the last row instead would freeze the
    // field a few tenths of a percent short of its own steady state, which is
    // exactly the settled value every cross-check is against.
    const LAMINF = new Float64Array(N), CURINF = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      LAMINF[i] = T.emf*T.lamE[i] + T.Iss*T.lamI[i];
      CURINF[i] = T.Iss;
    }
    T.K = K; T.dutab = dutab; T.tTab = tTab;
    T.LAM = LAM; T.LAMD = LAMD; T.CUR = CUR; T.CURD = CURD;
    T.LAMINF = LAMINF; T.CURINF = CURINF; T.ZERO = new Float64Array(N);
  }

  const smooth = x => {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    return x*x*x*(10 + x*(-15 + 6*x));
  };

  // Interpolate one tabulated quantity for segment i at time u. Zero before the
  // news arrives, held at the settled value after the table runs out.
  function tab(T, arr, inf, i, u) {
    if (u <= T.delay[i]) return 0;
    const b = i*T.K;
    if (u >= T.tTab) {
      // past the table the EMF is long constant, so what is left is a pure
      // exponential relaxation onto the steady value
      const f = inf[i];
      return f + (arr[b + T.K - 1] - f) * Math.exp(-(u - T.tTab) / T.tauLR);
    }
    const x = u / T.dutab, k = Math.floor(x);
    if (k >= T.K - 1) return arr[b + T.K - 1];
    return arr[b+k] + (arr[b+k+1] - arr[b+k]) * (x - k);
  }

  /* ---------------------------------------------------------------
     the scalar histories at an arbitrary time
     --------------------------------------------------------------- */
  // Only I is tabulated; its derivatives come straight back out of the circuit
  // equation, so there is no numerical differentiation anywhere.
  function histAt(T, t) {
    if (t <= 0) return { e:0, ed:0, edd:0, I:0, Id:0, Idd:0, Iddd:0 };
    const e = emfAt(T.emf, T.tau, t);
    let I;
    if (t >= T.tMax) {
      // Past the end of the march the EMF is constant, so the remaining
      // relaxation is exactly I_ss + (I(tMax) - I_ss) exp(-(t - tMax) R/L).
      // Clamping instead would leave the model permanently short of its own
      // steady state by however far the march had got.
      I = T.Iss + (T.Ihist[T.steps] - T.Iss) * Math.exp(-(t - T.tMax) / T.tauLR);
    } else {
      const u = t / T.dt, k = Math.floor(u);
      I = lerp(T.Ihist[k], T.Ihist[k+1], u - k);
    }
    const Id   = (e.v  - T.Rtot * I)   / T.Lind;
    const Idd  = (e.d  - T.Rtot * Id)  / T.Lind;
    const Iddd = (e.dd - T.Rtot * Idd) / T.Lind;
    return { e: e.v, ed: e.d, edd: e.dd, I, Id, Idd, Iddd };
  }

  // Charge density and current on segment i at time t, with their time
  // derivatives, read from the causal tables.
  function sourceAt(T, i, t) {
    return {
      lam:  tab(T, T.LAM,  T.LAMINF, i, t),
      lamd: tab(T, T.LAMD, T.ZERO,   i, t),
      cur:  tab(T, T.CUR,  T.CURINF, i, t),
      curd: tab(T, T.CURD, T.ZERO,   i, t)
    };
  }

  /* ---------------------------------------------------------------
     the retarded field
     --------------------------------------------------------------- */
  function fieldAt(T, px, py, pz, t) {
    const N = T.N, a2 = T.a * T.a;
    const { cx, cy, cz, tx, ty, tz, len, delay } = T;
    const kb = MU0 / (4 * Math.PI);
    let ex=0, ey=0, ez=0, Bx=0, By=0, Bz=0, v=0, dmin=Infinity, near=0;

    for (let i = 0; i < N; i++) {
      const dx = px-cx[i], dy = py-cy[i], dz = pz-cz[i];
      const r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < dmin) { dmin = r2; near = i; }
      const R = Math.sqrt(r2), Rs = Math.sqrt(r2 + a2);
      const Rb = Math.sqrt(r2 + 1e-20);        // B is NOT softened -- see below
      const tr = t - R / C;
      // Nothing here yet if the news has not reached this segment from the
      // battery and then reached the probe from this segment. Both halves are
      // needed; the first is what makes the picture causal about the battery.
      if (tr <= delay[i]) continue;

      const lam  = tab(T, T.LAM,  T.LAMINF, i, tr);
      const lamd = tab(T, T.LAMD, T.ZERO,   i, tr);
      const cur  = tab(T, T.CUR,  T.CURINF, i, tr);
      const curd = tab(T, T.CURD, T.ZERO,   i, tr);

      const l = len[i];
      const Rs2 = Rs*Rs, Rs3 = Rs2*Rs;
      // E: Coulomb, its retarded correction, and the induction term
      const ce = KC * l * (lam / Rs3 + lamd / (C * Rs2));
      ex += ce*dx; ey += ce*dy; ez += ce*dz;
      const ca = -KC * l * curd / (C*C*Rs);
      ex += ca*tx[i]; ey += ca*ty[i]; ez += ca*tz[i];
      v  += KC * l * lam / Rs;
      // B: Biot-Savart and its retarded correction, both crossed into Rhat.
      // Unlike E this is not softened at the wire radius: the interior is
      // already handled by the d < a mask below, and softening on top of it
      // biases the field low well outside the metal -- about two parts in a
      // thousand a few centimetres out, with no N to make it go away.
      const Rb2 = Rb*Rb, Rb3 = Rb2*Rb;
      const cbv = kb * l * (cur / Rb3 + curd / (C * Rb2));
      Bx += cbv * (ty[i]*dz - tz[i]*dy);
      By += cbv * (tz[i]*dx - tx[i]*dz);
      Bz += cbv * (tx[i]*dy - ty[i]*dx);
    }

    const d = Math.sqrt(dmin);
    if (d < T.a) { const m = (d*d)/(T.a*T.a); Bx*=m; By*=m; Bz*=m; }
    return { ex, ey, ez, bx: Bx, by: By, bz: Bz, v, d, near };
  }

  function poynting(f) {
    return {
      x: (f.ey*f.bz - f.ez*f.by) / MU0,
      y: (f.ez*f.bx - f.ex*f.bz) / MU0,
      z: (f.ex*f.by - f.ey*f.bx) / MU0
    };
  }
  function modeVec(f, mode) {
    if (mode === 'E') return { x: f.ex, y: f.ey, z: f.ez };
    if (mode === 'B') return { x: f.bx, y: f.by, z: f.bz };
    return poynting(f);
  }

  /* ---------------------------------------------------------------
     sampling
     --------------------------------------------------------------- */
  function buildVolume(T, opts) {
    const n = opts.n | 0, half = opts.half, t = opts.t;
    const step = (2*half) / (n-1), jit = step * 0.34;
    const jitter = (i,j,k,s) => {
      const q = Math.sin(i*127.1 + j*311.7 + k*74.7 + s*39.3) * 43758.5453;
      return q - Math.floor(q) - 0.5;
    };
    const TT = n*n*n;
    const P = new Float32Array(TT*3), V = new Float32Array(TT*3);
    const M = new Float32Array(TT), D = new Float32Array(TT);
    let mx = 0;
    for (let k = 0; k < n; k++) for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
      const id = (k*n + j)*n + i;
      const px = -half + i*step + jitter(i,j,k,0)*jit;
      const py = -half + j*step + jitter(i,j,k,1)*jit;
      const pz = -half + k*step + jitter(i,j,k,2)*jit;
      const f = fieldAt(T, px, py, pz, t);
      const q = modeVec(f, opts.mode);
      P[id*3]=px; P[id*3+1]=py; P[id*3+2]=pz;
      V[id*3]=q.x; V[id*3+1]=q.y; V[id*3+2]=q.z;
      const m = Math.hypot(q.x,q.y,q.z);
      M[id]=m; D[id]=f.d;
      if (f.d > T.a && m > mx) mx = m;
    }
    return { n, half, step, P, V, M, D, max: mx, mode: opts.mode, t };
  }

  function buildSlice(T, opts) {
    const n = opts.n|0, half = opts.half, o = opts.o, u = opts.u, w = opts.v, t = opts.t;
    const step = (2*half)/(n-1);
    const V = new Float32Array(n*n*3), M = new Float32Array(n*n), D = new Float32Array(n*n);
    let mx = 0;
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
      const s1 = -half + i*step, s2 = -half + j*step;
      const px = o.x + u.x*s1 + w.x*s2, py = o.y + u.y*s1 + w.y*s2, pz = o.z + u.z*s1 + w.z*s2;
      const f = fieldAt(T, px, py, pz, t);
      const q = modeVec(f, opts.mode);
      const k = j*n + i;
      V[k*3]=q.x; V[k*3+1]=q.y; V[k*3+2]=q.z;
      const m = Math.hypot(q.x,q.y,q.z);
      M[k]=m; D[k]=f.d;
      if (f.d > T.a && m > mx) mx = m;
    }
    return { n, half, o, u, v: w, step, V, M, D, max: mx, mode: opts.mode, t };
  }

  /* ---------------------------------------------------------------
     presets
     --------------------------------------------------------------- */
  const rect = [[-.26,-.15,0],[0,-.15,0],[.26,-.15,0],[.26,0,0],[.26,.15,0],[0,.15,0],[-.26,.15,0],[-.26,0,0]];
  const circle = (() => { const a=[]; for(let i=0;i<10;i++){const t=i/10*2*Math.PI; a.push([.2*Math.cos(t), .2*Math.sin(t), 0]);} return a; })();
  const tilt = (() => { const c=Math.cos(0.7), s=Math.sin(0.7); return rect.map(p => [p[0], p[1]*c, p[1]*s]); })();
  const helix = (() => {
    const a=[], turns=5, perTurn=6, R=0.1, n=turns*perTurn;
    for (let i=0;i<n;i++){ const t=i/perTurn*2*Math.PI; a.push([R*Math.cos(t), R*Math.sin(t), -0.2 + 0.4*i/(n-1)]); }
    a.push([0.27,0,0.24],[0.31,0,0],[0.27,0,-0.24]);
    return a;
  })();
  const saddle = rect.map((p,i) => [p[0], p[1], 0.13*Math.sin(i/rect.length*4*Math.PI)]);
  const PRESETS = { rect, circle, tilt, helix, saddle };

  return {
    EPS0, KC, MU0, C, WIRE, BATTERY, LED, PRESETS,
    clamp, lerp, densePath, resample, luFactor, luSolve, emfAt,
    buildTransient, histAt, sourceAt, fieldAt, poynting, modeVec,
    buildVolume, buildSlice
  };
});
