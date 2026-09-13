/*
 * app.js — canvas rendering and UI for the 3D field bench.
 * All physics lives in solver3d.js and all camera maths in view3d.js; this file
 * only decides how to draw them. There is no WebGL and no dependency: the scene
 * is projected by hand and painted back to front onto a 2D canvas, which is
 * enough for a wire, some arrows and a few dozen field lines.
 */
(function () {
  'use strict';

  const FB = window.FB3, V = window.V3;
  const { MU0, clamp, lerp } = FB;

  /* ---------------- number formatting ---------------- */
  const PREFIX = [[1e9,'G'],[1e6,'M'],[1e3,'k'],[1,''],[1e-3,'m'],[1e-6,'µ'],[1e-9,'n'],[1e-12,'p'],[1e-15,'f']];
  function fmt(v, unit) {
    if (!isFinite(v)) return '—';
    const a = Math.abs(v);
    if (a === 0) return '0 ' + unit;
    for (const [s, p] of PREFIX) {
      if (a >= s * 0.999) {
        const x = v / s;
        const d = Math.abs(x) >= 100 ? 0 : Math.abs(x) >= 10 ? 1 : 2;
        return x.toFixed(d) + ' ' + p + unit;
      }
    }
    return v.toExponential(2) + ' ' + unit;
  }

  /* ---------------- colour ramps ---------------- */
  const hex = h => [parseInt(h.slice(1,3),16), parseInt(h.slice(3,5),16), parseInt(h.slice(5,7),16)];
  function ramp(stops) {
    const lut = new Uint8Array(256 * 3);
    for (let i = 0; i < 256; i++) {
      const t = i / 255;
      let k = 0; while (k < stops.length - 2 && t > stops[k+1][0]) k++;
      const [t0, c0] = stops[k], [t1, c1] = stops[k+1];
      const f = clamp((t - t0) / (t1 - t0 || 1), 0, 1);
      lut[i*3]   = Math.round(lerp(c0[0], c1[0], f));
      lut[i*3+1] = Math.round(lerp(c0[1], c1[1], f));
      lut[i*3+2] = Math.round(lerp(c0[2], c1[2], f));
    }
    return lut;
  }
  const RAMP = {
    E: ramp([[0,hex('#05161d')],[.28,hex('#0d4f63')],[.55,hex('#2e9fb0')],[.8,hex('#8fe3e0')],[1,hex('#f4fffd')]]),
    S: ramp([[0,hex('#150f06')],[.3,hex('#7a3d10')],[.62,hex('#e8a33d')],[.85,hex('#ffd98a')],[1,hex('#fff8ea')]]),
    B: ramp([[0,hex('#1b1030')],[.3,hex('#4a2a7a')],[.58,hex('#9a5ad0')],[.82,hex('#d9a8f0')],[1,hex('#fdf2ff')]])
  };
  const LINE_TINT = { E: '160,214,220', B: '196,160,230', S: '232,190,120' };
  const rampCss = (mode, t) => {
    const i = clamp(Math.round(t * 255), 0, 255) * 3, L = RAMP[mode];
    return 'rgb(' + L[i] + ',' + L[i+1] + ',' + L[i+2] + ')';
  };

  /* ---------------- state ---------------- */
  const S = {
    ctrl: FB.PRESETS.rect.map(p => ({ x: p[0], y: p[1], z: p[2] })),
    mode: 'E',
    arrows: true, heat: false, lines: true, charge: false, frame: true,
    emf: 9, rLED: 100, rho: 0, radius: 0.003,
    cam: { az: 38, el: 22, dist: 1.15, fov: 45, target: { x: 0, y: 0, z: 0 } },
    slice: { axis: 'z', off: 0 },
    spin: false,
    drag: -1, hover: -1, dragging: false, orbiting: false, quality: 'high',
    sol: null, sl: null, vol: null, streams: null, range: null, scene: 0.4
  };

  const $ = id => document.getElementById(id);
  const cv = $('cv');
  const ctx = cv.getContext('2d', { alpha: false });
  let W = 0, H = 0, DPR = 1, F = null;      // F is the current camera frame

  /* ---------------- the slice plane ---------------- */
  function planeOf() {
    const n = S.slice.axis === 'x' ? { x:1, y:0, z:0 }
            : S.slice.axis === 'y' ? { x:0, y:1, z:0 }
                                   : { x:0, y:0, z:1 };
    const o = { x: n.x * S.slice.off, y: n.y * S.slice.off, z: n.z * S.slice.off };
    const [u, v] = V.planeAxes(n);
    return { n, o, u, v };
  }

  // How big the world is, from the circuit itself: everything else scales off it.
  function sceneRadius(sol) {
    let r = 0.05;
    for (let k = 0; k < 3; k++) r = Math.max(r, Math.abs(sol.lo[k]), Math.abs(sol.hi[k]));
    return r;
  }

  /* ---------------- colour range ---------------- */
  // Fields run over decades, so the ramp is logarithmic between two percentiles
  // of what is actually on the slice, ignoring samples inside the metal.
  function updateRange(vol, sol) {
    const vals = [];
    const T = vol.n * vol.n * vol.n;
    for (let k = 0; k < T; k++) if (vol.D[k] > sol.a * 1.2 && vol.M[k] > 0) vals.push(vol.M[k]);
    if (vals.length < 8) { S.range = null; return; }
    vals.sort((a, b) => a - b);
    const at = p => vals[clamp(Math.floor(p * (vals.length - 1)), 0, vals.length - 1)];
    let lo = at(0.30), hi = at(0.985);
    if (!(hi > lo * 1.0001)) hi = lo * 1.0001;
    S.range = { lo, hi, llo: Math.log(lo), lhi: Math.log(hi) };
  }
  const toT = m => {
    const r = S.range;
    if (!r || !(m > 0)) return 0;
    return clamp((Math.log(m) - r.llo) / (r.lhi - r.llo), 0, 1);
  };

  /* ---------------- field lines in space ---------------- */
  /*
   * Seeded from the slice so the controls stay honest -- move the plane and you
   * reseed -- but integrated in three dimensions, so a line is free to leave the
   * plane immediately and generally does. B lines close into rings around the
   * wire; E lines run from positive surface charge to negative.
   */
  function trace(sol, mode, p0, sign, step, maxSteps, bound2) {
    const pts = [[p0.x, p0.y, p0.z]];
    let x = p0.x, y = p0.y, z = p0.z, closed = false, run = 0;
    let lox = x, loy = y, loz = z, hix = x, hiy = y, hiz = z;
    const shut = (step * 0.9) * (step * 0.9);
    for (let i = 0; i < maxSteps; i++) {
      const f1 = FB.fieldAt(sol, x, y, z);
      if (f1.d < sol.a * 1.1) break;              // ran into the metal
      const v1 = FB.modeVec(f1, mode);
      const m1 = Math.hypot(v1.x, v1.y, v1.z);
      if (!(m1 > 0)) break;
      const h = step * sign;
      const mx = x + v1.x / m1 * h * 0.5, my = y + v1.y / m1 * h * 0.5, mz = z + v1.z / m1 * h * 0.5;
      const f2 = FB.fieldAt(sol, mx, my, mz);
      const v2 = FB.modeVec(f2, mode);
      const m2 = Math.hypot(v2.x, v2.y, v2.z);
      if (!(m2 > 0)) break;
      x += v2.x / m2 * h; y += v2.y / m2 * h; z += v2.z / m2 * h;
      pts.push([x, y, z]);
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
  const GRID = 4;              // seed lattice side; GRID^3 candidates
  function buildStreams(sol, mode) {
    const R = S.scene;
    const step = R * 0.022, maxSteps = 120, bound2 = (R * 2.3) * (R * 2.3);
    // Minimum separation between lines, as a fraction of the scene. Seeds are
    // rejected where an existing line already passes, so this both removes
    // duplicates and keeps neighbours apart. It is only ever a thinning rule:
    // because seeds are taken in lattice order, what survives is a regular
    // subset of the lattice rather than a greedy fill of whatever space is left.
    // Measured over the presets, raising it from 0.05 thins the crowded ones and
    // leaves the sparse ones alone: the 5-turn solenoid goes from 51 lines to 35
    // and its share of sampled points lying within 0.04R of another line from
    // 16.6% to 12.3%, the saddle from 60 to 44, while the flat loop stays at 28
    // because its lines were never close to begin with.
    //
    // The solenoid does not reach zero and should not. Field lines crowd where
    // the flux is concentrated, which inside a coil is the whole point; that
    // density is the physics, not a placement artefact.
    const dup = R * 0.12;
    const occupied = new Set();
    const key = (x, y, z) => Math.floor(x/dup) + ',' + Math.floor(y/dup) + ',' + Math.floor(z/dup);
    const near = (x, y, z) => {
      const i = Math.floor(x/dup), j = Math.floor(y/dup), k = Math.floor(z/dup);
      for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) for (let c = -1; c <= 1; c++)
        if (occupied.has((i+a) + ',' + (j+b) + ',' + (k+c))) return true;
      return false;
    };

    const out = [];
    const span = R * 1.15;
    for (let k = 0; k < GRID; k++) {
      for (let j = 0; j < GRID; j++) {
        for (let i = 0; i < GRID; i++) {
          const p = {
            x: (2 * (i + 0.5) / GRID - 1) * span,
            y: (2 * (j + 0.5) / GRID - 1) * span,
            z: (2 * (k + 0.5) / GRID - 1) * span
          };
          const f = FB.fieldAt(sol, p.x, p.y, p.z);
          if (f.d < sol.a * 2.5) continue;                 // seed is in the metal
          const q = FB.modeVec(f, mode);
          if (!(Math.hypot(q.x, q.y, q.z) > 0)) continue;  // nothing to follow
          if (near(p.x, p.y, p.z)) continue;               // already on a line drawn
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
          for (const w of line) occupied.add(key(w[0], w[1], w[2]));
          out.push(line);
        }
      }
    }
    return out;
  }

  /* ---------------- drawing ---------------- */

  // Depth fade: far things get fainter, so the eye can order them.
  const fade = d => clamp(1.25 - d / (S.cam.dist * 1.7), 0.12, 1);

  function drawHeat(sl, proj) {
    const n = sl.n;
    for (let j = 0; j < n - 1; j++) {
      for (let i = 0; i < n - 1; i++) {
        const a = proj[j*n+i], b = proj[j*n+i+1], c = proj[(j+1)*n+i+1], d = proj[(j+1)*n+i];
        if (!a || !b || !c || !d) continue;
        const m = (sl.M[j*n+i] + sl.M[j*n+i+1] + sl.M[(j+1)*n+i+1] + sl.M[(j+1)*n+i]) / 4;
        ctx.fillStyle = rampCss(S.mode, toT(m));
        ctx.beginPath();
        ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]);
        ctx.lineTo(c[0], c[1]); ctx.lineTo(d[0], d[1]);
        ctx.closePath(); ctx.fill();
      }
    }
  }

  function drawSliceFrame(pl, h) {
    const c = [];
    for (const [a, b] of [[-1,-1],[1,-1],[1,1],[-1,1]]) {
      const p = V.project(F,
        pl.o.x + pl.u.x*h*a + pl.v.x*h*b,
        pl.o.y + pl.u.y*h*a + pl.v.y*h*b,
        pl.o.z + pl.u.z*h*a + pl.v.z*h*b);
      if (!p) return;
      c.push(p);
    }
    ctx.beginPath();
    ctx.moveTo(c[0][0], c[0][1]);
    for (let i = 1; i < 4; i++) ctx.lineTo(c[i][0], c[i][1]);
    ctx.closePath();
    ctx.strokeStyle = 'rgba(120,170,185,.32)'; ctx.lineWidth = 1; ctx.setLineDash([5, 5]);
    ctx.stroke(); ctx.setLineDash([]);
  }

  /*
   * One arrow per volume sample, throughout the space rather than on a surface.
   * Arrow length is a fraction of the distance between neighbouring samples, so
   * arrows keep the same relationship to the lattice at any zoom instead of
   * growing into each other. Foreshortening is left alone: an arrow pointing at
   * the camera is meant to look short, that is the depth cue. Weak field needs
   * no culling -- the low end of the ramp is the background colour, so those
   * arrows fade out on their own.
   */
  function drawArrows(vol, sol, wantFar, split) {
    const T = vol.n * vol.n * vol.n;
    for (let k = 0; k < T; k++) {
      const m = vol.M[k];
      if (!(m > 0) || vol.D[k] < sol.a * 1.3) continue;
      const px = vol.P[k*3], py = vol.P[k*3+1], pz = vol.P[k*3+2];
      const p = V.project(F, px, py, pz);
      if (!p) continue;
      if ((p[2] > split) !== wantFar) continue;
      const tt = toT(m);
      const cell = vol.step * F.focal / p[2];          // sample spacing, in pixels
      const pix = clamp(cell * (0.3 + 0.45 * tt), 3.5, 44);
      const w = pix * p[2] / F.focal;                  // world length giving that many pixels
      const vx = vol.V[k*3] / m, vy = vol.V[k*3+1] / m, vz = vol.V[k*3+2] / m;
      const q = V.project(F, px + vx*w, py + vy*w, pz + vz*w);
      if (!q) continue;
      const dx = q[0] - p[0], dy = q[1] - p[1];
      const L = Math.hypot(dx, dy);
      ctx.strokeStyle = ctx.fillStyle = rampCss(S.mode, tt);
      ctx.globalAlpha = fade(p[2]) * (0.16 + 0.84 * tt);
      ctx.lineWidth = 1.4;
      ctx.beginPath(); ctx.moveTo(p[0], p[1]); ctx.lineTo(q[0], q[1]); ctx.stroke();
      if (L > 4) {                                     // head, only when there is room
        const ux = dx / L, uy = dy / L, hs = clamp(L * 0.38, 2.5, 5.5);
        ctx.beginPath();
        ctx.moveTo(q[0], q[1]);
        ctx.lineTo(q[0] - ux*hs - uy*hs*0.55, q[1] - uy*hs + ux*hs*0.55);
        ctx.lineTo(q[0] - ux*hs + uy*hs*0.55, q[1] - uy*hs - ux*hs*0.55);
        ctx.closePath(); ctx.fill();
      }
      ctx.globalAlpha = 1;
    }
  }

  // Field lines are split against a single depth threshold rather than sorted
  // per segment: the wire is one compact object, so "in front of the loop" or
  // "behind it" is the only ordering that actually reads.
  function drawStreams(streams, wantFar, split) {
    if (!streams) return;
    const tint = LINE_TINT[S.mode];
    ctx.lineWidth = 1.25; ctx.lineJoin = 'round'; ctx.lineCap = 'round';
    for (const line of streams) {
      let run = null, acc = 0, cnt = 0;
      const flush = () => {
        if (run && run.length > 1) {
          ctx.strokeStyle = 'rgba(' + tint + ',' + (fade(acc / cnt) * 0.5).toFixed(3) + ')';
          ctx.beginPath();
          ctx.moveTo(run[0][0], run[0][1]);
          for (let i = 1; i < run.length; i++) ctx.lineTo(run[i][0], run[i][1]);
          ctx.stroke();
        }
        run = null; acc = 0; cnt = 0;
      };
      for (const w of line) {
        const p = V.project(F, w[0], w[1], w[2]);
        if (!p || (p[2] > split) !== wantFar) { flush(); continue; }
        if (!run) run = [];
        run.push(p); acc += p[2]; cnt++;
      }
      flush();
    }
  }

  /*
   * The wire and its charge beads go into one list sorted back to front. This is
   * the one place true depth order matters: a solenoid has to pass behind itself.
   */
  function drawCircuit(sol) {
    const items = [];
    const N = sol.N;
    for (let i = 0; i < N; i++) {
      const a = V.project(F, sol.ax[i], sol.ay[i], sol.az[i]);
      const b = V.project(F, sol.bx[i], sol.by[i], sol.bz[i]);
      if (!a || !b) continue;
      items.push({ d: (a[2] + b[2]) / 2, kind: 0, a, b, t: sol.type[i] });
    }
    if (S.charge) for (const c of chargeMarks(sol)) {
      const p = V.project(F, c.x, c.y, c.z);
      if (p) items.push({ d: p[2], kind: 1, p, sign: c.sign });
    }
    items.sort((p, q) => q.d - p.d);

    const rWorld = Math.max(sol.a, 0.0016);
    for (const it of items) {
      if (it.kind === 0) {
        const wpx = clamp(2 * rWorld * F.focal / it.d, 1.6, 60);
        const t = it.t;
        // The halo is what makes a near part of the wire read as passing in front
        // of a far part, so it is drawn per segment rather than once for the run.
        // Its caps must be butt: round caps reach half a halo-width past each end,
        // which on segments this short punches a bite out of both neighbours and
        // turns the wire into a string of beads.
        ctx.lineCap = 'butt';
        ctx.strokeStyle = 'rgba(3,12,16,.92)'; ctx.lineWidth = wpx + 4.5;
        ctx.beginPath(); ctx.moveTo(it.a[0], it.a[1]); ctx.lineTo(it.b[0], it.b[1]); ctx.stroke();
        ctx.lineCap = 'round';
        ctx.strokeStyle = t === FB.BATTERY ? '#7fe6d2' : t === FB.LED ? '#ffb765' : '#c6d7da';
        ctx.lineWidth = t === FB.WIRE ? wpx : wpx + 3;
        ctx.beginPath(); ctx.moveTo(it.a[0], it.a[1]); ctx.lineTo(it.b[0], it.b[1]); ctx.stroke();
      } else {
        const R = clamp(0.006 * F.focal / it.d, 2.4, 9), h = R * 0.52;
        ctx.beginPath(); ctx.arc(it.p[0], it.p[1], R, 0, 6.2832);
        ctx.fillStyle = it.sign > 0 ? '#ff8a63' : '#6cc8ff'; ctx.fill();
        ctx.strokeStyle = 'rgba(4,14,19,.9)'; ctx.lineWidth = 1.2; ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(it.p[0] - h, it.p[1]); ctx.lineTo(it.p[0] + h, it.p[1]);
        if (it.sign > 0) { ctx.moveTo(it.p[0], it.p[1] - h); ctx.lineTo(it.p[0], it.p[1] + h); }
        ctx.strokeStyle = '#06171d'; ctx.lineWidth = 1.5; ctx.stroke();
      }
    }
  }

  /*
   * One mark per equal quantum of charge, so marks per unit length is the charge
   * density itself, and + and - come out equal in number because the solve holds
   * the loop neutral. Offset radially away from the loop centroid so the beads
   * sit on the outside of the metal rather than down its axis.
   */
  function chargeMarks(sol) {
    const N = sol.N, lam = sol.lam, len = sol.len;
    let Q = 0;
    for (let i = 0; i < N; i++) Q += Math.abs(lam[i]) * len[i];
    if (!(Q > 0)) return [];
    const M = clamp(Math.round(sol.L / (S.scene * 0.065)), 16, 240);
    const quantum = Q / M;
    let gx = 0, gy = 0, gz = 0;
    for (let i = 0; i < N; i++) { gx += sol.cx[i]; gy += sol.cy[i]; gz += sol.cz[i]; }
    gx /= N; gy /= N; gz /= N;
    const off = sol.a + S.scene * 0.022;

    const out = [];
    let acc = 0, next = quantum * 0.5;
    for (let i = 0; i < N; i++) {
      const q = Math.abs(lam[i]) * len[i];
      while (next < acc + q) {
        const f = (next - acc) / q;
        const x = sol.ax[i] + (sol.bx[i] - sol.ax[i]) * f;
        const y = sol.ay[i] + (sol.by[i] - sol.ay[i]) * f;
        const z = sol.az[i] + (sol.bz[i] - sol.az[i]) * f;
        let nx = x - gx, ny = y - gy, nz = z - gz;
        const d = nx*sol.tx[i] + ny*sol.ty[i] + nz*sol.tz[i];
        nx -= d*sol.tx[i]; ny -= d*sol.ty[i]; nz -= d*sol.tz[i];
        const nn = Math.hypot(nx, ny, nz) || 1;
        out.push({ x: x + nx/nn*off, y: y + ny/nn*off, z: z + nz/nn*off, sign: lam[i] >= 0 ? 1 : -1 });
        next += quantum;
      }
      acc += q;
    }
    return out;
  }

  function drawHandles() {
    const list = [];
    S.ctrl.forEach((p, i) => {
      const s = V.project(F, p.x, p.y, p.z);
      if (s) list.push({ s, i });
    });
    list.sort((a, b) => b.s[2] - a.s[2]);
    for (const { s, i } of list) {
      const act = i === S.hover || i === S.drag;
      const shrink = clamp(14 / S.ctrl.length, 0.5, 1);   // a 5-turn coil has 33
      const r = clamp(0.008 * F.focal / s[2], 3.5, 11) * shrink * (act ? 1.35 : 1);
      ctx.beginPath(); ctx.arc(s[0], s[1], r, 0, 7);
      ctx.fillStyle = act ? '#ffffff' : 'rgba(240,252,252,.8)'; ctx.fill();
      ctx.lineWidth = 2; ctx.strokeStyle = 'rgba(5,20,26,.85)'; ctx.stroke();
    }
  }

  // A corner triad, so which way is up is never a guess.
  function drawTriad() {
    const cx = W / 2, cy = H - 34, len = 26;
    const axes = [['x', 1,0,0, '#ff8a63'], ['y', 0,1,0, '#7fe6d2'], ['z', 0,0,1, '#a98af0']];
    ctx.font = '500 10.5px "IBM Plex Sans", sans-serif';
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    for (const [label, x, y, z, col] of axes) {
      const sx = (x*F.right.x + y*F.right.y + z*F.right.z) * len;
      const sy = -(x*F.up.x + y*F.up.y + z*F.up.z) * len;
      ctx.strokeStyle = col; ctx.lineWidth = 1.6; ctx.globalAlpha = .85;
      ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(cx + sx, cy + sy); ctx.stroke();
      ctx.fillStyle = col;
      ctx.fillText(label, cx + sx * 1.28, cy + sy * 1.28);
      ctx.globalAlpha = 1;
    }
  }

  /* ---------------- the frame ---------------- */
  function render() {
    ctx.fillStyle = '#05141a'; ctx.fillRect(0, 0, W, H);
    const sol = S.sol, sl = S.sl;
    if (!sol) return;
    F = V.frame(S.cam, W, H);

    // everything is ordered against the depth of the loop itself
    let gx = 0, gy = 0, gz = 0;
    for (let i = 0; i < sol.N; i++) { gx += sol.cx[i]; gy += sol.cy[i]; gz += sol.cz[i]; }
    const split = V.depthOf(F, gx / sol.N, gy / sol.N, gz / sol.N);

    if (sl && S.heat) {
      const n = sl.n, proj = new Array(n * n);
      for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
        const t = -sl.half + i * sl.step, b = -sl.half + j * sl.step;
        proj[j*n+i] = V.project(F,
          sl.o.x + sl.u.x*t + sl.v.x*b,
          sl.o.y + sl.u.y*t + sl.v.y*b,
          sl.o.z + sl.u.z*t + sl.v.z*b);
      }
      drawHeat(sl, proj);
    }
    if (S.frame) drawSliceFrame(planeOf(), S.scene * 1.55);

    if (S.lines) drawStreams(S.streams, true, split);
    if (S.vol && S.arrows) drawArrows(S.vol, sol, true, split);

    drawCircuit(sol);

    if (S.lines) drawStreams(S.streams, false, split);
    if (S.vol && S.arrows) drawArrows(S.vol, sol, false, split);

    drawHandles();
    drawTriad();
    drawLegend();
  }

  /* ---------------- recompute ---------------- */
  let pending = false;
  function recompute() {
    if (pending) return;
    pending = true;
    requestAnimationFrame(() => {
      pending = false;
      const hi = S.quality === 'high';
      const sol = FB.solveCircuit({
        ctrl: S.ctrl, N: hi ? 180 : 120, per: hi ? 44 : 26,
        radius: S.radius, emf: S.emf, rLED: S.rLED, rho: S.rho
      });
      if (!sol) { ctx.fillStyle = '#05141a'; ctx.fillRect(0, 0, W, H); return; }
      S.sol = sol;
      S.scene = sceneRadius(sol);
      S.vol = FB.buildVolume(sol, {
        half: S.scene * 1.25,
        n: hi ? 9 : 6,
        mode: S.mode
      });
      updateRange(S.vol, sol);
      // The slice is sampled only for the heatmap now; the outline needs just
      // the plane, and the readout evaluates the field where the ray lands.
      if (S.heat) {
        const pl = planeOf();
        S.sl = FB.buildSlice(sol, {
          o: pl.o, u: pl.u, v: pl.v,
          half: S.scene * 1.55,
          n: hi ? 40 : 20,
          mode: S.mode
        });
      } else S.sl = null;
      S.streams = (hi && S.lines) ? buildStreams(sol, S.mode) : null;
      render();
      updateFacts(sol);
    });
  }
  function settle() {
    S.quality = 'high';
    recompute();
  }

  /* ---------------- layout ---------------- */
  function resize() {
    const r = cv.getBoundingClientRect();
    if (r.width < 1 || r.height < 1) return;
    DPR = Math.min(window.devicePixelRatio || 1, 2);
    W = Math.round(r.width); H = Math.round(r.height);
    cv.width = Math.round(W*DPR); cv.height = Math.round(H*DPR);
    ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
    recompute();
  }

  /* ---------------- UI ---------------- */
  const MODE_NOTE = {
    E: 'The electric field of the surface charge the circuit builds up. Along ideal wire it meets the metal at right angles — no tangential push, no work done. At the LED it tips over and dives in.',
    B: 'The magnetic field of the current. In three dimensions it is a full vector: field lines close into rings threaded on the wire, and a solenoid gathers them into a bundle down its axis.',
    S: 'Energy flow, S = E×B/μ₀. It leaves the battery sideways into empty space, glides along outside the wire, and turns inward wherever power is being dissipated. The wire guides the energy; it does not carry it.'
  };
  const LEG_CAP = { E: 'Field strength |E|', B: 'Field strength |B|', S: 'Energy flux |S|' };
  const LEG_UNIT = { E: 'V/m', B: 'T', S: 'W/m²' };

  function drawLegend() {
    const r = S.range;
    $('leg-cap').textContent = LEG_CAP[S.mode];
    if (!r) { $('leg-lo').textContent = $('leg-hi').textContent = '—'; return; }
    $('leg-lo').textContent = fmt(r.lo, LEG_UNIT[S.mode]);
    $('leg-hi').textContent = fmt(r.hi, LEG_UNIT[S.mode]);
    const stops = [];
    for (let i = 0; i <= 6; i++) stops.push(rampCss(S.mode, i / 6) + ' ' + Math.round(i / 6 * 100) + '%');
    $('leg-bar').style.background = 'linear-gradient(90deg,' + stops.join(',') + ')';
  }

  function updateFacts(sol) {
    const P = sol.I * sol.I * S.rLED;
    $('f-I').textContent = fmt(sol.I, 'A');
    $('f-P').textContent = fmt(P, 'W');
    $('f-L').textContent = (sol.L * 100).toFixed(1) + ' cm';
    $('f-Ew').textContent = fmt(S.rho * sol.I, 'V/m');
    $('f-El').textContent = fmt(sol.I * S.rLED / sol.ledLen, 'V/m');
    $('f-q').textContent = fmt(sol.peak, 'C/m');
    $('f-n').textContent = S.ctrl.length;
  }

  function syncCam() {
    $('s-az').value = Math.round(S.cam.az);
    $('s-el').value = Math.round(S.cam.el);
    $('s-dist').value = Math.round(S.cam.dist * 100);
    $('s-fov').value = Math.round(S.cam.fov);
    $('v-az').textContent = Math.round(S.cam.az) + '°';
    $('v-el').textContent = Math.round(S.cam.el) + '°';
    $('v-dist').textContent = S.cam.dist.toFixed(2) + ' m';
    $('v-fov').textContent = Math.round(S.cam.fov) + '°';
  }

  // mode
  $('modeseg').addEventListener('click', e => {
    const b = e.target.closest('button'); if (!b) return;
    S.mode = b.dataset.mode;
    [...$('modeseg').children].forEach(x => x.setAttribute('aria-pressed', String(x === b)));
    $('modenote').textContent = MODE_NOTE[S.mode];
    settle();
  });

  // slice axis
  $('axisseg').addEventListener('click', e => {
    const b = e.target.closest('button'); if (!b) return;
    S.slice.axis = b.dataset.axis;
    [...$('axisseg').children].forEach(x => x.setAttribute('aria-pressed', String(x === b)));
    settle();
  });

  // checkboxes
  [['c-arrows','arrows'], ['c-heat','heat'], ['c-lines','lines'],
   ['c-charge','charge'], ['c-frame','frame']].forEach(([id, key]) => {
    $(id).addEventListener('change', e => {
      S[key] = e.target.checked;
      const needsData = (key === 'lines' && e.target.checked && !S.streams) ||
                        (key === 'heat' && e.target.checked && !S.sl);
      if (needsData) settle(); else render();
    });
  });
  $('c-spin').addEventListener('change', e => { S.spin = e.target.checked; if (S.spin) tick(); });

  // circuit sliders
  const SLIDERS = [
    ['s-emf','v-emf', v => { S.emf = v; return v.toFixed(1) + ' V'; }],
    ['s-led','v-led', v => { S.rLED = v; return v.toFixed(0) + ' Ω'; }],
    ['s-rho','v-rho', v => { S.rho = v; return v.toFixed(0) + ' Ω/m'; }],
    ['s-rad','v-rad', v => { S.radius = v/1000; return v.toFixed(1) + ' mm'; }]
  ];
  for (const [sid, vid, apply] of SLIDERS) {
    const el = $(sid);
    const run = () => { $(vid).textContent = apply(parseFloat(el.value)); };
    run();
    el.addEventListener('input', () => { run(); S.quality = 'low'; recompute(); bounce(); });
  }

  // slice offset, in centimetres on the slider
  $('s-off').addEventListener('input', e => {
    S.slice.off = parseFloat(e.target.value) / 1000;
    $('v-off').textContent = (S.slice.off * 100).toFixed(1) + ' cm';
    S.quality = 'low'; recompute(); bounce();
  });

  // camera sliders
  const CAM = [['s-az','az',1], ['s-el','el',1], ['s-dist','dist',0.01], ['s-fov','fov',1]];
  for (const [sid, key, scale] of CAM) {
    $(sid).addEventListener('input', e => {
      S.cam[key] = parseFloat(e.target.value) * scale;
      syncCam(); render();
    });
  }

  const VIEWS = {
    iso:   { az: 38,  el: 22 },
    front: { az: -90, el: 0 },
    top:   { az: -90, el: 89 },
    side:  { az: 0,   el: 0 },
    reset: { az: 38,  el: 22, dist: 1.15, fov: 45 }
  };
  $('viewchips').addEventListener('click', e => {
    const b = e.target.closest('button'); if (!b) return;
    Object.assign(S.cam, VIEWS[b.dataset.view]);
    syncCam(); render();
  });

  $('chips').addEventListener('click', e => {
    const b = e.target.closest('button'); if (!b) return;
    S.ctrl = FB.PRESETS[b.dataset.p].map(p => ({ x: p[0], y: p[1], z: p[2] }));
    settle();
  });

  $('railtoggle').addEventListener('click', () => {
    const open = $('rail').classList.toggle('open');
    $('railtoggle').setAttribute('aria-expanded', String(open));
    $('railtoggle').textContent = open ? 'Hide' : 'Controls';
  });

  let rt = null;
  const bounce = () => { clearTimeout(rt); rt = setTimeout(settle, 220); };

  /* ---------------- pointer ---------------- */
  const local = e => {
    const r = cv.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top];
  };
  function hitHandle(sx, sy) {
    let best = -1, bd = 18 * 18;
    for (let i = 0; i < S.ctrl.length; i++) {
      const p = V.project(F, S.ctrl[i].x, S.ctrl[i].y, S.ctrl[i].z);
      if (!p) continue;
      const d = (p[0]-sx)*(p[0]-sx) + (p[1]-sy)*(p[1]-sy);
      if (d < bd) { bd = d; best = i; }
    }
    return best;
  }

  const pointers = new Map();
  let last = null, dragDepth = 0, pinch0 = 0, dist0 = 0;

  cv.addEventListener('pointerdown', e => {
    cv.setPointerCapture(e.pointerId);
    pointers.set(e.pointerId, local(e));
    if (pointers.size === 2) {                       // pinch to zoom
      const [a, b] = [...pointers.values()];
      pinch0 = Math.hypot(a[0]-b[0], a[1]-b[1]); dist0 = S.cam.dist;
      S.dragging = false; S.orbiting = false;
      return;
    }
    const [sx, sy] = local(e);
    const i = hitHandle(sx, sy);
    if (i >= 0) {
      S.drag = i; S.dragging = true;
      dragDepth = V.depthOf(F, S.ctrl[i].x, S.ctrl[i].y, S.ctrl[i].z);
      S.quality = 'low'; S.streams = null;
    } else {
      S.orbiting = true; cv.classList.add('grabbing');
    }
    last = [sx, sy];
  });

  cv.addEventListener('pointermove', e => {
    const cur = local(e);
    if (pointers.has(e.pointerId)) pointers.set(e.pointerId, cur);

    if (pointers.size === 2 && pinch0 > 0) {
      const [a, b] = [...pointers.values()];
      const d = Math.hypot(a[0]-b[0], a[1]-b[1]);
      S.cam.dist = clamp(dist0 * pinch0 / Math.max(d, 1), 0.3, 3);
      syncCam(); render();
      return;
    }

    if (S.dragging && S.drag >= 0) {
      const w = V.unproject(F, cur[0], cur[1], dragDepth);
      const lim = S.scene * 2.2;
      S.ctrl[S.drag] = { x: clamp(w.x, -lim, lim), y: clamp(w.y, -lim, lim), z: clamp(w.z, -lim, lim) };
      recompute(); return;
    }
    if (S.orbiting && last) {
      S.cam.az = ((S.cam.az - (cur[0] - last[0]) * 0.42 + 180) % 360 + 360) % 360 - 180;
      S.cam.el = clamp(S.cam.el + (cur[1] - last[1]) * 0.35, -V.EL_MAX, V.EL_MAX);
      last = cur; syncCam(); render(); return;
    }

    // hover: light up a handle, and read the field where the ray meets the slice
    const h = hitHandle(cur[0], cur[1]);
    if (h !== S.hover) { S.hover = h; render(); }
    probe(cur[0], cur[1]);
  });

  function endPointer(e) {
    pointers.delete(e.pointerId);
    if (pointers.size < 2) pinch0 = 0;
    if (S.dragging) { S.dragging = false; S.drag = -1; bounce(); }
    S.orbiting = false;
    cv.classList.remove('grabbing');
    last = null;
  }
  cv.addEventListener('pointerup', endPointer);
  cv.addEventListener('pointercancel', endPointer);
  cv.addEventListener('pointerleave', () => {
    S.hover = -1;
    ['ro-e','ro-b','ro-s','ro-v'].forEach(id => { $(id).textContent = '—'; });
    render();
  });

  cv.addEventListener('wheel', e => {
    e.preventDefault();
    S.cam.dist = clamp(S.cam.dist * Math.exp(e.deltaY * 0.0012), 0.3, 3);
    syncCam(); render();
  }, { passive: false });

  // Double-click adds a handle, on the segment nearest the click.
  cv.addEventListener('dblclick', e => {
    const [sx, sy] = local(e);
    if (!S.sol) return;
    let best = -1, bd = Infinity;
    for (let i = 0; i < S.sol.N; i++) {
      const p = V.project(F, S.sol.cx[i], S.sol.cy[i], S.sol.cz[i]);
      if (!p) continue;
      const d = (p[0]-sx)*(p[0]-sx) + (p[1]-sy)*(p[1]-sy);
      if (d < bd) { bd = d; best = i; }
    }
    if (best < 0 || bd > 30*30) return;
    // insert between the two control points the nearest segment falls between
    const s = S.sol.s[best] / S.sol.L, k = Math.floor(s * S.ctrl.length);
    const a = S.ctrl[k], b = S.ctrl[(k+1) % S.ctrl.length];
    S.ctrl.splice(k + 1, 0, { x: (a.x+b.x)/2, y: (a.y+b.y)/2, z: (a.z+b.z)/2 });
    settle();
  });

  /* ---------------- probe ---------------- */
  // A screen point is a ray, not a place, so the readout reports where that ray
  // pierces the slice plane. Move the plane and you move the probe surface.
  function probe(sx, sy) {
    if (!S.sol || !F) return;
    const pl = planeOf();
    const far = V.unproject(F, sx, sy, F.focal);
    const dx = far.x - F.eye.x, dy = far.y - F.eye.y, dz = far.z - F.eye.z;
    const den = dx*pl.n.x + dy*pl.n.y + dz*pl.n.z;
    if (Math.abs(den) < 1e-9) return;
    const t = ((pl.o.x - F.eye.x)*pl.n.x + (pl.o.y - F.eye.y)*pl.n.y + (pl.o.z - F.eye.z)*pl.n.z) / den;
    if (t <= 0) return;
    const p = { x: F.eye.x + dx*t, y: F.eye.y + dy*t, z: F.eye.z + dz*t };
    const f = FB.fieldAt(S.sol, p.x, p.y, p.z);
    const s = FB.poynting(f);
    $('ro-e').textContent = fmt(Math.hypot(f.ex, f.ey, f.ez), 'V/m');
    $('ro-b').textContent = fmt(Math.hypot(f.bx, f.by, f.bz), 'T');
    $('ro-s').textContent = fmt(Math.hypot(s.x, s.y, s.z), 'W/m²');
    $('ro-v').textContent = fmt(f.v - S.sol.V0, 'V');
  }

  /* ---------------- spin ---------------- */
  let spinning = false;
  function tick() {
    if (spinning) return;
    spinning = true;
    let t0 = performance.now();
    const step = now => {
      if (!S.spin) { spinning = false; bounce(); return; }
      const dt = Math.min(now - t0, 60); t0 = now;
      S.cam.az = ((S.cam.az + dt * 0.012 + 180) % 360 + 360) % 360 - 180;
      syncCam(); render();
      requestAnimationFrame(step);
    };
    requestAnimationFrame(step);
  }

  /* ---------------- go ---------------- */
  window.addEventListener('resize', () => {
    S.quality = 'low'; resize(); bounce();
  });
  $('modenote').textContent = MODE_NOTE.E;
  $('v-off').textContent = '0.0 cm';
  syncCam();
  resize();
})();
