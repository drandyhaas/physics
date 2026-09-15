/*
 * app.js — rendering and UI for the switch-on bench.
 * All physics lives in solver4d.js and all camera maths in view3d.js. There is
 * no WebGL and no dependency: the scene is projected by hand and painted back to
 * front onto a 2D canvas.
 *
 * The one structural difference from the steady-state benches is that there are
 * now two costs, and they are kept apart. Solving the circuit is expensive and
 * depends only on the geometry and the component values, so it happens when
 * those change. Sampling the field is cheap-ish and depends on the clock, so it
 * happens every frame. Moving the clock never re-solves anything.
 */
(function () {
  'use strict';

  const FB = window.FB4, V = window.V3;
  const { clamp, lerp, C } = FB;
  const NS = 1e-9;

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
    const lut = new Uint8Array(256*3);
    for (let i = 0; i < 256; i++) {
      const t = i/255;
      let k = 0; while (k < stops.length-2 && t > stops[k+1][0]) k++;
      const [t0,c0] = stops[k], [t1,c1] = stops[k+1];
      const f = clamp((t-t0)/(t1-t0 || 1), 0, 1);
      lut[i*3]   = Math.round(lerp(c0[0], c1[0], f));
      lut[i*3+1] = Math.round(lerp(c0[1], c1[1], f));
      lut[i*3+2] = Math.round(lerp(c0[2], c1[2], f));
    }
    return lut;
  }
  // Every ramp starts at the background colour exactly. On the steady benches the
  // bottom of the ramp was a dark tint, which was harmless because there is no
  // such thing as an unreached point there. Here there is: outside the light
  // cone the field is exactly zero, and a tinted floor paints that region as
  // "a little bit of field" and wipes out the very boundary this bench is for.
  const BG = '#05141a';
  const RAMP = {
    E: ramp([[0,hex(BG)],[.28,hex('#0d4f63')],[.55,hex('#2e9fb0')],[.8,hex('#8fe3e0')],[1,hex('#f4fffd')]]),
    S: ramp([[0,hex(BG)],[.3,hex('#7a3d10')],[.62,hex('#e8a33d')],[.85,hex('#ffd98a')],[1,hex('#fff8ea')]]),
    B: ramp([[0,hex(BG)],[.3,hex('#4a2a7a')],[.58,hex('#9a5ad0')],[.82,hex('#d9a8f0')],[1,hex('#fdf2ff')]])
  };
  const LINE_TINT = { E: '160,214,220', B: '196,160,230', S: '232,190,120' };
  const rampCss = (mode, t) => {
    const i = clamp(Math.round(t*255), 0, 255)*3, L = RAMP[mode];
    return 'rgb(' + L[i] + ',' + L[i+1] + ',' + L[i+2] + ')';
  };

  /* ---------------- state ---------------- */
  const S = {
    ctrl: FB.PRESETS.rect.map(p => ({ x:p[0], y:p[1], z:p[2] })),
    modes: { E: true, B: true, S: false },
    arrows: true, heat: false, charge: true, front: true, frame: true, lines: true,
    arrowSize: 1.0,
    emf: 9, rLED: 100, rho: 0, radius: 0.003,
    tau: 1e-9, t: 0, tEnd: 40e-9, playing: false, speed: 50, lineDensity: 1,
    cam: { az: 38, el: 22, dist: 1.5, fov: 45, target: { x:0, y:0, z:0 } },
    slice: { axis: 'z', off: 0 },
    drag: -1, hover: -1, dragging: false, orbiting: false, scrubbing: false,
    T: null, vols: {}, sls: {}, streams: {}, marks: null, ranges: {}, scene: 0.4
  };

  const $ = id => document.getElementById(id);
  const cv = $('cv');
  const ctx = cv.getContext('2d', { alpha: false });
  const tcv = $('trace');
  const tctx = tcv.getContext('2d');
  let W = 0, H = 0, DPR = 1, F = null;

  /* ---------------- the slice plane ---------------- */
  function planeOf() {
    const n = S.slice.axis === 'x' ? {x:1,y:0,z:0} : S.slice.axis === 'y' ? {x:0,y:1,z:0} : {x:0,y:0,z:1};
    const o = { x:n.x*S.slice.off, y:n.y*S.slice.off, z:n.z*S.slice.off };
    const [u, v] = V.planeAxes(n);
    return { n, o, u, v };
  }
  function sceneRadius(T) {
    let r = 0.05;
    for (let k = 0; k < 3; k++) r = Math.max(r, Math.abs(T.lo[k]), Math.abs(T.hi[k]));
    return r;
  }

  /* ---------------- colour range ----------------
     Fixed to the SETTLED field, not to the present instant. A range recomputed
     each frame would rescale itself as the field grows, and the switch-on -- the
     one thing this bench exists to show -- would be invisible: everything would
     look the same brightness at every moment. */
  function updateRange(T) {
    const half = S.scene * 1.8, n = 9, step = (2*half)/(n-1);
    const tLate = 900e-9;
    const bag = { E: [], B: [], S: [] };
    for (let k = 0; k < n; k++) for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
      const x = -half + i*step, y = -half + j*step, z = -half + k*step;
      const f = FB.fieldAt(T, x, y, z, tLate + Math.hypot(x,y,z)/C);
      if (f.d < T.a*1.2) continue;
      // one field evaluation serves all three, which is why several modes at
      // once costs far less than several times as much
      for (const m of ['E','B','S']) {
        const q = FB.modeVec(f, m);
        const v = Math.hypot(q.x, q.y, q.z);
        if (v > 0) bag[m].push(v);
      }
    }
    S.ranges = {};
    for (const m of ['E','B','S']) {
      const vals = bag[m];
      if (vals.length < 8) continue;
      vals.sort((a,b) => a-b);
      const at = p => vals[clamp(Math.floor(p*(vals.length-1)), 0, vals.length-1)];
      let lo = at(0.18), hi = at(0.985);
      if (!(hi > lo*1.0001)) hi = lo*1.0001;
      S.ranges[m] = { lo, hi, llo: Math.log(lo), lhi: Math.log(hi) };
    }
  }
  const toT = (mode, m) => {
    const r = S.ranges[mode];
    if (!r || !(m > 0)) return 0;
    return clamp((Math.log(m) - r.llo) / (r.lhi - r.llo), 0, 1);
  };
  const active = () => ['E','B','S'].filter(m => S.modes[m]);

  /* ---------------- drawing ---------------- */
  const fade = d => clamp(1.25 - d/(S.cam.dist*1.7), 0.12, 1);

  function drawHeat(sl, proj, mode) {
    const n = sl.n;
    for (let j = 0; j < n-1; j++) for (let i = 0; i < n-1; i++) {
      const a = proj[j*n+i], b = proj[j*n+i+1], c = proj[(j+1)*n+i+1], d = proj[(j+1)*n+i];
      if (!a || !b || !c || !d) continue;
      const m = (sl.M[j*n+i] + sl.M[j*n+i+1] + sl.M[(j+1)*n+i+1] + sl.M[(j+1)*n+i]) / 4;
      ctx.fillStyle = rampCss(mode, toT(mode, m));
      ctx.beginPath();
      ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]); ctx.lineTo(c[0],c[1]); ctx.lineTo(d[0],d[1]);
      ctx.closePath(); ctx.fill();
    }
  }

  function drawSliceFrame(pl, h) {
    const c = [];
    for (const [a,b] of [[-1,-1],[1,-1],[1,1],[-1,1]]) {
      const p = V.project(F, pl.o.x+pl.u.x*h*a+pl.v.x*h*b, pl.o.y+pl.u.y*h*a+pl.v.y*h*b, pl.o.z+pl.u.z*h*a+pl.v.z*h*b);
      if (!p) return;
      c.push(p);
    }
    ctx.beginPath(); ctx.moveTo(c[0][0], c[0][1]);
    for (let i = 1; i < 4; i++) ctx.lineTo(c[i][0], c[i][1]);
    ctx.closePath();
    ctx.strokeStyle = 'rgba(120,170,185,.3)'; ctx.lineWidth = 1; ctx.setLineDash([5,5]);
    ctx.stroke(); ctx.setLineDash([]);
  }

  function drawArrows(vol, T, mode, wantFar, split) {
    const TT = vol.n*vol.n*vol.n;
    // 29 px at full strength and a 1x slider, a quarter of the on-screen gap
    // between neighbouring lattice points at the default zoom.
    const base = 29 * S.arrowSize;
    for (let k = 0; k < TT; k++) {
      const m = vol.M[k];
      if (!(m > 0) || vol.D[k] < T.a*1.3) continue;
      const px = vol.P[k*3], py = vol.P[k*3+1], pz = vol.P[k*3+2];
      const p = V.project(F, px, py, pz);
      if (!p) continue;
      if ((p[2] > split) !== wantFar) continue;
      const tt = toT(mode, m);
      const pix = base * (0.35 + 0.65*tt);
      const wl = pix * p[2] / F.focal;
      const vx = vol.V[k*3]/m, vy = vol.V[k*3+1]/m, vz = vol.V[k*3+2]/m;
      const q = V.project(F, px+vx*wl, py+vy*wl, pz+vz*wl);
      if (!q) continue;
      const dx = q[0]-p[0], dy = q[1]-p[1], L = Math.hypot(dx, dy);
      ctx.strokeStyle = ctx.fillStyle = rampCss(mode, tt);
      ctx.globalAlpha = fade(p[2]) * (0.16 + 0.84*tt);
      ctx.lineWidth = 1.4;
      ctx.beginPath(); ctx.moveTo(p[0],p[1]); ctx.lineTo(q[0],q[1]); ctx.stroke();
      if (L > 4) {
        const ux = dx/L, uy = dy/L, hs = clamp(L*0.38, 2.5, base*0.3);
        ctx.beginPath();
        ctx.moveTo(q[0],q[1]);
        ctx.lineTo(q[0]-ux*hs-uy*hs*0.55, q[1]-uy*hs+ux*hs*0.55);
        ctx.lineTo(q[0]-ux*hs+uy*hs*0.55, q[1]-uy*hs-ux*hs*0.55);
        ctx.closePath(); ctx.fill();
      }
      ctx.globalAlpha = 1;
    }
  }

  // Centred on the battery, because that is the only thing that changes: the
  // rest of the circuit cannot respond until the news gets to it. The field is
  // strictly zero outside this sphere -- it is a boundary, not a guide.
  function drawFront(T, t) {
    const r = C * t;
    if (!(r > 0)) return;
    const c = V.project(F, T.battery.x, T.battery.y, T.battery.z);
    if (!c) return;
    const rp = r * F.focal / c[2];
    if (!(rp > 2) || rp > 20000) return;
    ctx.beginPath(); ctx.arc(c[0], c[1], rp, 0, 6.2832);
    ctx.strokeStyle = 'rgba(232,163,61,.45)'; ctx.lineWidth = 1.4;
    ctx.setLineDash([7,6]); ctx.stroke(); ctx.setLineDash([]);
  }

  function drawCircuit(T, t) {
    const items = [];
    const N = T.N;
    for (let i = 0; i < N; i++) {
      const a = V.project(F, T.ax[i], T.ay[i], T.az[i]);
      const b = V.project(F, T.bx[i], T.by[i], T.bz[i]);
      if (!a || !b) continue;
      items.push({ d:(a[2]+b[2])/2, kind:0, a, b, t:T.type[i], lit: T.delay[i] < t });
    }
    if (S.charge && S.marks) for (const m of S.marks) {
      const st = markState(T, m, t);
      if (st.fill < 0.02) continue;                 // nothing has arrived here yet
      const p = V.project(F, m.x, m.y, m.z);
      if (p) items.push({ d:p[2], kind:1, p, sign:st.sign, fill:st.fill });
    }
    items.sort((p,q) => q.d - p.d);

    const rWorld = Math.max(T.a, 0.0016);
    for (const it of items) {
      if (it.kind === 0) {
        const wpx = clamp(2*rWorld*F.focal/it.d, 1.6, 60);
        ctx.lineCap = 'butt';
        ctx.strokeStyle = 'rgba(3,12,16,.92)'; ctx.lineWidth = wpx + 4.5;
        ctx.beginPath(); ctx.moveTo(it.a[0],it.a[1]); ctx.lineTo(it.b[0],it.b[1]); ctx.stroke();
        ctx.lineCap = 'round';
        // Wire the news has not reached yet is drawn dark. This is the honest
        // way to show where ct has got to along the wire: it is continuous and
        // exact, where the charge marks are quantised and cannot be. There is
        // almost no charge just behind the front anyway -- lambda is zero at the
        // front and takes the whole switch-on time to come up, by which point
        // the front has run most of the way round -- so a mark placed there
        // would claim a whole quantum where there is half a percent of one.
        ctx.strokeStyle = it.lit
          ? (it.t === FB.BATTERY ? '#7fe6d2' : it.t === FB.LED ? '#ffb765' : '#c6d7da')
          : (it.t === FB.BATTERY ? '#2c4f49' : it.t === FB.LED ? '#4a3a22' : '#33444a');
        ctx.lineWidth = it.t === FB.WIRE ? wpx : wpx + 3;
        ctx.beginPath(); ctx.moveTo(it.a[0],it.a[1]); ctx.lineTo(it.b[0],it.b[1]); ctx.stroke();
      } else {
        const R = clamp(0.006*F.focal/it.d, 2.4, 9), h = R*0.52;
        ctx.globalAlpha = it.fill;
        ctx.beginPath(); ctx.arc(it.p[0], it.p[1], R, 0, 6.2832);
        ctx.fillStyle = it.sign > 0 ? '#ff8a63' : '#6cc8ff'; ctx.fill();
        ctx.strokeStyle = 'rgba(4,14,19,.9)'; ctx.lineWidth = 1.2; ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(it.p[0]-h, it.p[1]); ctx.lineTo(it.p[0]+h, it.p[1]);
        if (it.sign > 0) { ctx.moveTo(it.p[0], it.p[1]-h); ctx.lineTo(it.p[0], it.p[1]+h); }
        ctx.strokeStyle = '#06171d'; ctx.lineWidth = 1.5; ctx.stroke();
        ctx.globalAlpha = 1;
      }
    }
  }

  /*
   * Charge marks.
   *
   * Positions are fixed, computed once from the SETTLED charge, one mark per
   * equal quantum of it. What varies with time is how filled each mark is:
   * |lambda| there now, over |lambda| there once settled.
   *
   * Placing them by integrating |lambda| at the current instant instead -- the
   * obvious thing, and what this did at first -- makes every mark slide toward
   * whatever point the integration starts from as the charge grows, because the
   * place where the running total first reaches one quantum moves inward as
   * lambda rises everywhere. Anchored at the battery and run outward both ways
   * it reads unmistakably as charge flowing INTO the battery, which is not
   * happening. The drift is an artefact of quantising a growing quantity.
   *
   * Fixing the positions removes the motion entirely: marks fade up in place, in
   * the order the news reaches them, which is outward from the battery. The
   * count still carries the charge -- each mark is one quantum of settled charge
   * and is drawn at its fill fraction, so the total ink is the charge present.
   */
  function buildMarks(T) {
    const N = T.N;
    const ref = new Float64Array(N);
    let Qss = 0;
    for (let i = 0; i < N; i++) {
      ref[i] = T.emf*T.lamE[i] + T.Iss*T.lamI[i];
      Qss += Math.abs(ref[i]) * T.len[i];
    }
    if (!(Qss > 0)) return [];
    const M = clamp(Math.round(T.Lloop/(S.scene*0.065)), 16, 240);
    const quantum = Qss / M;

    let gx=0, gy=0, gz=0;
    for (let i = 0; i < N; i++) { gx += T.cx[i]; gy += T.cy[i]; gz += T.cz[i]; }
    gx /= N; gy /= N; gz /= N;
    const off = T.a + S.scene*0.022;

    const out = [];
    let acc = 0, next = quantum*0.5;
    for (let i = 0; i < N; i++) {
      const q = Math.abs(ref[i]) * T.len[i];
      while (next < acc + q) {
        const f = (next - acc)/q;
        const x = T.ax[i] + (T.bx[i]-T.ax[i])*f;
        const y = T.ay[i] + (T.by[i]-T.ay[i])*f;
        const z = T.az[i] + (T.bz[i]-T.az[i])*f;
        let nx = x-gx, ny = y-gy, nz = z-gz;
        const d = nx*T.tx[i] + ny*T.ty[i] + nz*T.tz[i];
        nx -= d*T.tx[i]; ny -= d*T.ty[i]; nz -= d*T.tz[i];
        const nn = Math.hypot(nx,ny,nz) || 1;
        out.push({ i, x: x+nx/nn*off, y: y+ny/nn*off, z: z+nz/nn*off, ref: Math.abs(ref[i]) });
        next += quantum;
      }
      acc += q;
    }
    return out;
  }

  // How full each mark is at this instant, and which sign it is carrying.
  function markState(T, m, t) {
    const lam = FB.sourceAt(T, m.i, t).lam;
    return { fill: m.ref > 0 ? clamp(Math.abs(lam)/m.ref, 0, 1) : 0, sign: lam >= 0 ? 1 : -1 };
  }
  function drawStreams(streams, mode, wantFar, split) {
    if (!streams) return;
    const tint = LINE_TINT[mode];
    ctx.lineWidth = 1.25; ctx.lineJoin = 'round'; ctx.lineCap = 'round';
    for (const line of streams) {
      let run = null, acc = 0, cnt = 0;
      const flush = () => {
        if (run && run.length > 1) {
          ctx.strokeStyle = 'rgba(' + tint + ',' + (fade(acc/cnt)*0.55).toFixed(3) + ')';
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

  function drawHandles() {
    const list = [];
    S.ctrl.forEach((p,i) => { const s = V.project(F, p.x, p.y, p.z); if (s) list.push({s,i}); });
    list.sort((a,b) => b.s[2] - a.s[2]);
    const shrink = clamp(14/S.ctrl.length, 0.5, 1);
    for (const {s,i} of list) {
      const act = i === S.hover || i === S.drag;
      const r = clamp(0.008*F.focal/s[2], 3.5, 11) * shrink * (act ? 1.35 : 1);
      ctx.beginPath(); ctx.arc(s[0], s[1], r, 0, 7);
      ctx.fillStyle = act ? '#fff' : 'rgba(240,252,252,.8)'; ctx.fill();
      ctx.lineWidth = 2; ctx.strokeStyle = 'rgba(5,20,26,.85)'; ctx.stroke();
    }
  }

  function drawTriad() {
    const cx = W/2, cy = H-34, len = 26;
    ctx.font = '500 10.5px "IBM Plex Sans", sans-serif';
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    for (const [label,x,y,z,col] of [['x',1,0,0,'#ff8a63'],['y',0,1,0,'#7fe6d2'],['z',0,0,1,'#a98af0']]) {
      const sx = (x*F.right.x + y*F.right.y + z*F.right.z)*len;
      const sy = -(x*F.up.x + y*F.up.y + z*F.up.z)*len;
      ctx.strokeStyle = col; ctx.lineWidth = 1.6; ctx.globalAlpha = .85;
      ctx.beginPath(); ctx.moveTo(cx,cy); ctx.lineTo(cx+sx, cy+sy); ctx.stroke();
      ctx.fillStyle = col; ctx.fillText(label, cx+sx*1.28, cy+sy*1.28);
      ctx.globalAlpha = 1;
    }
  }

  /* ---------------- the EMF / current trace ---------------- */
  function drawTrace() {
    const T = S.T;
    const w = tcv.width, h = tcv.height;
    tctx.setTransform(1,0,0,1,0,0);
    tctx.fillStyle = '#12333d'; tctx.fillRect(0,0,w,h);
    if (!T) return;
    const pad = 10, x0 = pad, x1 = w-pad, y0 = pad, y1 = h-pad-14;
    const X = t => x0 + (t/S.tEnd)*(x1-x0);
    const Y = f => y1 - f*(y1-y0);
    tctx.strokeStyle = 'rgba(120,170,185,.25)'; tctx.lineWidth = 1;
    tctx.beginPath(); tctx.moveTo(x0,y1); tctx.lineTo(x1,y1); tctx.stroke();

    const curve = (get, col, width) => {
      tctx.beginPath();
      for (let i = 0; i <= 160; i++) {
        const t = S.tEnd*i/160, v = get(t);
        const px = X(t), py = Y(v);
        i ? tctx.lineTo(px,py) : tctx.moveTo(px,py);
      }
      tctx.strokeStyle = col; tctx.lineWidth = width; tctx.stroke();
    };
    curve(t => FB.emfAt(T.emf, T.tau, t).v / T.emf, 'rgba(127,230,210,.9)', 2);
    curve(t => FB.histAt(T, t).I / T.Iss, '#e8a33d', 2.2);

    const cx = X(S.t);
    tctx.strokeStyle = 'rgba(255,255,255,.55)'; tctx.lineWidth = 1;
    tctx.beginPath(); tctx.moveTo(cx,y0-2); tctx.lineTo(cx,y1+4); tctx.stroke();

    tctx.font = '500 17px "IBM Plex Sans", sans-serif'; tctx.textBaseline = 'top';
    tctx.textAlign = 'left';
    tctx.fillStyle = 'rgba(127,230,210,.9)'; tctx.fillText('EMF', x0, y1+8);
    tctx.fillStyle = '#e8a33d'; tctx.fillText('current', x0+52, y1+8);
    tctx.textAlign = 'right';
    tctx.fillStyle = 'rgba(139,176,184,.9)';
    tctx.fillText((S.tEnd/NS).toFixed(0) + ' ns', x1, y1+8);
  }

  /* ---------------- the frame ---------------- */
  function render() {
    ctx.fillStyle = '#05141a'; ctx.fillRect(0,0,W,H);
    const T = S.T;
    if (!T) return;
    F = V.frame(S.cam, W, H);

    let gx=0, gy=0, gz=0;
    for (let i = 0; i < T.N; i++) { gx += T.cx[i]; gy += T.cy[i]; gz += T.cz[i]; }
    const split = V.depthOf(F, gx/T.N, gy/T.N, gz/T.N);

    const on = active();
    // The colour map can only carry one field at a time -- two of them over the
    // same plane is mud -- so it takes the first one switched on. Arrows and
    // field lines layer perfectly well, each in its own ramp.
    const heatMode = on[0];
    if (S.heat && heatMode && S.sls[heatMode]) {
      const sl = S.sls[heatMode], n = sl.n, proj = new Array(n*n);
      for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
        const a = -sl.half + i*sl.step, b = -sl.half + j*sl.step;
        proj[j*n+i] = V.project(F, sl.o.x+sl.u.x*a+sl.v.x*b, sl.o.y+sl.u.y*a+sl.v.y*b, sl.o.z+sl.u.z*a+sl.v.z*b);
      }
      drawHeat(sl, proj, heatMode);
    }
    if (S.frame) drawSliceFrame(planeOf(), S.scene*2.6);
    for (const m of on) {
      if (S.lines) drawStreams(S.streams[m], m, true, split);
      if (S.arrows && S.vols[m]) drawArrows(S.vols[m], T, m, true, split);
    }
    drawCircuit(T, S.t);
    for (const m of on) {
      if (S.lines) drawStreams(S.streams[m], m, false, split);
      if (S.arrows && S.vols[m]) drawArrows(S.vols[m], T, m, false, split);
    }
    if (S.front) drawFront(T, S.t);
    drawHandles();
    drawTriad();
    drawLegend();
    drawTrace();
    updateClock();
  }

  /* ---------------- solve (slow) and sample (per frame) ---------------- */
  function solve() {
    const T = FB.buildTransient({
      ctrl: S.ctrl, N: 140, per: 40, radius: S.radius,
      emf: S.emf, rLED: S.rLED, rho: S.rho,
      tau: S.tau, tMax: 200*NS, steps: 8000
    });
    if (!T) return;
    S.T = T;
    S.scene = sceneRadius(T);
    S.tEnd = clamp(4*T.tauLR, 20*NS, 120*NS);
    $('s-t').max = String(Math.round(S.tEnd/NS*10));
    updateRange(T);
    S.marks = buildMarks(T);
    updateFacts();
    sample();
  }

  function sample() {
    const T = S.T;
    if (!T) return;
    // Only hand interactions get the coarse path. Playback draws exactly what a
    // stopped frame draws -- the movie is the thing worth looking at, and a
    // slower clock is a fairer price than a different picture.
    const fast = S.dragging || S.scrubbing;
    const on = active();
    S.vols = {}; S.sls = {}; S.streams = {};
    if (S.arrows) for (const m of on)
      S.vols[m] = FB.buildVolume(T, { half: S.scene*1.8, n: fast ? 7 : 9, mode: m, t: S.t });
    if (S.heat && on.length) {
      const pl = planeOf(), m = on[0];
      S.sls[m] = FB.buildSlice(T, {
        o: pl.o, u: pl.u, v: pl.v, half: S.scene*2.6,
        n: fast ? 34 : 60, mode: m, t: S.t
      });
    }
    // Tracing is the one cost that really does multiply with the number of
    // fields shown, because each follows its own path. Lines get shorter when
    // more than one is running so a frame stays inside its budget.
    if (S.lines) for (const m of on)
      S.streams[m] = FB.buildStreams(T, { mode: m, t: S.t, fast: fast, nModes: on.length,
                                          scene: S.scene, density: S.lineDensity });
    render();
  }

  /* ---------------- layout ---------------- */
  function resize() {
    const r = cv.getBoundingClientRect();
    if (r.width < 1 || r.height < 1) return;
    DPR = Math.min(window.devicePixelRatio || 1, 2);
    W = Math.round(r.width); H = Math.round(r.height);
    cv.width = Math.round(W*DPR); cv.height = Math.round(H*DPR);
    ctx.setTransform(DPR,0,0,DPR,0,0);
    render();
  }

  /* ---------------- UI ---------------- */
  const MODE_NOTE = {
    E: 'The electric field. Before the switch there is none; then a front of it sweeps outward at c, and behind the front the surface charge settles into the pattern the steady bench draws.',
    B: 'The magnetic field of the current. It cannot appear until the current does, and the current is held back by the loop’s own inductance — so B lags the EMF by L/R, and the news of it still travels at c.',
    S: 'Energy flow, S = E×B/μ₀. Zero until both fields exist at a point, which is why the load learns about the battery no faster than light can carry the news.'
  };
  const LEG_CAP = { E:'Field strength |E|', B:'Field strength |B|', S:'Energy flux |S|' };
  const LEG_LABEL = { E:'E', B:'B', S:'E\u00d7B' };
  const LEG_UNIT = { E:'V/m', B:'T', S:'W/m²' };

  function drawLegend() {
    const on = active();
    const box = $('legend-rows');
    if (!on.length) { box.innerHTML = '<div class="cap">nothing selected</div>'; return; }
    let html = '';
    for (const m of on) {
      const r = S.ranges[m];
      const stops = [];
      for (let i = 0; i <= 6; i++) stops.push(rampCss(m, i/6) + ' ' + Math.round(i/6*100) + '%');
      html += '<div class="legrow">' +
        '<div class="cap">' + LEG_CAP[m] + '</div>' +
        '<div class="bar" style="background:linear-gradient(90deg,' + stops.join(',') + ')"></div>' +
        '<div class="ends"><span>' + (r ? fmt(r.lo, LEG_UNIT[m]) : '—') + '</span>' +
        '<span>' + (r ? fmt(r.hi, LEG_UNIT[m]) : '—') + '</span></div></div>';
    }
    box.innerHTML = html;
  }

  function updateClock() {
    $('t-now').textContent = (S.t/NS).toFixed(2) + ' ns';
    $('t-sub').textContent = 'light has reached ' + (C*S.t).toFixed(2) + ' m';
    $('v-t').textContent = (S.t/NS).toFixed(2) + ' ns';
    $('s-t').value = String(Math.round(S.t/NS*10));
    const T = S.T;
    if (!T) return;
    const h = FB.histAt(T, S.t);
    $('f-e').textContent = fmt(h.e, 'V');
    $('f-I').textContent = fmt(h.I, 'A');
  }

  function updateFacts() {
    const T = S.T;
    if (!T) return;
    $('f-Iss').textContent = fmt(T.Iss, 'A');
    $('f-L').textContent = fmt(T.Lind, 'H');
    $('f-lr').textContent = (T.tauLR/NS).toFixed(2) + ' ns';
    $('f-cross').textContent = (T.cross/NS).toFixed(2) + ' ns';
    $('f-len').textContent = (T.Lloop*100).toFixed(1) + ' cm';
  }

  function syncCam() {
    $('s-az').value = Math.round(S.cam.az);
    $('s-el').value = Math.round(S.cam.el);
    $('s-dist').value = Math.round(S.cam.dist*100);
    $('s-fov').value = Math.round(S.cam.fov);
    $('v-az').textContent = Math.round(S.cam.az) + '°';
    $('v-el').textContent = Math.round(S.cam.el) + '°';
    $('v-dist').textContent = S.cam.dist.toFixed(2) + ' m';
    $('v-fov').textContent = Math.round(S.cam.fov) + '°';
  }

  const nsPerSec = () => 0.1 * Math.pow(10, S.speed/50);   // 0 -> 0.1, 50 -> 1, 100 -> 10 ns/s

  // Independent toggles, not a one-of-three picker: the fields layer.
  function syncModes() {
    [...$('modeseg').children].forEach(x =>
      x.setAttribute('aria-pressed', String(!!S.modes[x.dataset.mode])));
    const on = active();
    $('modenote').innerHTML = on.length
      ? on.map(m => '<b>' + LEG_LABEL[m] + '</b> ' + MODE_NOTE[m]).join('<br><br>')
      : 'Nothing selected.';
  }
  $('modeseg').addEventListener('click', e => {
    const b = e.target.closest('button'); if (!b) return;
    const m = b.dataset.mode;
    // never leave all three off; the last one on stays on
    if (S.modes[m] && active().length === 1) return;
    S.modes[m] = !S.modes[m];
    syncModes(); sample();
  });
  $('axisseg').addEventListener('click', e => {
    const b = e.target.closest('button'); if (!b) return;
    S.slice.axis = b.dataset.axis;
    [...$('axisseg').children].forEach(x => x.setAttribute('aria-pressed', String(x === b)));
    sample();
  });
  [['c-arrows','arrows'],['c-heat','heat'],['c-lines','lines'],['c-charge','charge'],['c-front','front'],['c-frame','frame']]
    .forEach(([id,key]) => $(id).addEventListener('change', e => { S[key] = e.target.checked; sample(); }));

  $('btn-play').addEventListener('click', () => setPlaying(!S.playing));
  $('btn-restart').addEventListener('click', () => { S.t = 0; setPlaying(true); });
  function setPlaying(v) {
    // No-op when nothing changes. Scrubbing calls this on every input event, and
    // without the guard each one resampled the whole scene at FULL quality --
    // before the scrub flag was even set -- and then again at the fast one. Two
    // samples per frame, the expensive one wasted.
    if (S.playing === v) return;
    S.playing = v;
    $('btn-play').textContent = v ? 'Pause' : 'Play';
    $('btn-play').setAttribute('aria-pressed', String(v));
    if (v) tick();
    else sample();
  }

  let scrubTimer = null;
  $('s-t').addEventListener('input', e => {
    S.t = parseFloat(e.target.value)/10*NS;
    S.scrubbing = true;          // before setPlaying, so any sample it does is cheap
    setPlaying(false);
    sample();
    clearTimeout(scrubTimer);
    scrubTimer = setTimeout(() => { S.scrubbing = false; sample(); }, 200);
  });
  $('s-tau').addEventListener('input', e => {
    S.tau = parseFloat(e.target.value)/10*NS;
    $('v-tau').textContent = (S.tau/NS).toFixed(1) + ' ns';
    bounce();
  });
  $('s-speed').addEventListener('input', e => {
    S.speed = parseFloat(e.target.value);
    $('v-speed').textContent = nsPerSec().toFixed(1) + ' ns/s';
  });
  // 0.2x to 5x, logarithmic and symmetric about 1x at the middle
  const ldensOf = v => Math.pow(5, (v - 50)/50);
  $('s-ldens').addEventListener('input', e => {
    S.lineDensity = ldensOf(parseFloat(e.target.value));
    $('v-ldens').textContent = S.lineDensity.toFixed(2) + '×';
    sample();
  });

  const asizeOf = v => Math.pow(10, (v - 50)/50);      // 0 -> 0.1x, 50 -> 1x, 100 -> 10x
  const asizeText = a => (a < 0.995 ? a.toFixed(2) : a < 9.95 ? a.toFixed(2) : a.toFixed(1)) + '×';
  $('s-asize').addEventListener('input', e => {
    S.arrowSize = asizeOf(parseFloat(e.target.value));
    $('v-asize').textContent = asizeText(S.arrowSize);
    render();
  });
  $('s-off').addEventListener('input', e => {
    S.slice.off = parseFloat(e.target.value)/1000;
    $('v-off').textContent = (S.slice.off*100).toFixed(1) + ' cm';
    sample();
  });

  const SLIDERS = [
    ['s-emf','v-emf', v => { S.emf = v; return v.toFixed(1) + ' V'; }],
    ['s-led','v-led', v => { S.rLED = v; return v.toFixed(0) + ' Ω'; }],
    ['s-rho','v-rho', v => { S.rho = v; return v.toFixed(0) + ' Ω/m'; }],
    ['s-rad','v-rad', v => { S.radius = v/1000; return v.toFixed(1) + ' mm'; }]
  ];
  for (const [sid,vid,apply] of SLIDERS) {
    const el = $(sid);
    const run = () => { $(vid).textContent = apply(parseFloat(el.value)); };
    run();
    el.addEventListener('input', () => { run(); bounce(); });
  }

  const CAM = [['s-az','az',1],['s-el','el',1],['s-dist','dist',0.01],['s-fov','fov',1]];
  for (const [sid,key,scale] of CAM)
    $(sid).addEventListener('input', e => { S.cam[key] = parseFloat(e.target.value)*scale; syncCam(); render(); });

  const VIEWS = {
    iso:{az:38,el:22}, front:{az:-90,el:0}, top:{az:-90,el:89}, side:{az:0,el:0},
    reset:{az:38,el:22,dist:1.5,fov:45}
  };
  $('viewchips').addEventListener('click', e => {
    const b = e.target.closest('button'); if (!b) return;
    Object.assign(S.cam, VIEWS[b.dataset.view]); syncCam(); render();
  });
  $('chips').addEventListener('click', e => {
    const b = e.target.closest('button'); if (!b) return;
    S.ctrl = FB.PRESETS[b.dataset.p].map(p => ({x:p[0],y:p[1],z:p[2]}));
    solve();
  });
  $('railtoggle').addEventListener('click', () => {
    const open = $('rail').classList.toggle('open');
    $('railtoggle').setAttribute('aria-expanded', String(open));
    $('railtoggle').textContent = open ? 'Hide' : 'Controls';
  });

  let rt = null;
  const bounce = () => { clearTimeout(rt); rt = setTimeout(solve, 200); };

  /* ---------------- pointer ---------------- */
  const local = e => { const r = cv.getBoundingClientRect(); return [e.clientX-r.left, e.clientY-r.top]; };
  function hitHandle(sx, sy) {
    let best = -1, bd = 18*18;
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
    if (pointers.size === 2) {
      const [a,b] = [...pointers.values()];
      pinch0 = Math.hypot(a[0]-b[0], a[1]-b[1]); dist0 = S.cam.dist;
      S.dragging = false; S.orbiting = false; return;
    }
    const [sx,sy] = local(e);
    const i = hitHandle(sx,sy);
    if (i >= 0) {
      S.drag = i; S.dragging = true;
      dragDepth = V.depthOf(F, S.ctrl[i].x, S.ctrl[i].y, S.ctrl[i].z);
    } else { S.orbiting = true; cv.classList.add('grabbing'); }
    last = [sx,sy];
  });
  cv.addEventListener('pointermove', e => {
    const cur = local(e);
    if (pointers.has(e.pointerId)) pointers.set(e.pointerId, cur);
    if (pointers.size === 2 && pinch0 > 0) {
      const [a,b] = [...pointers.values()];
      S.cam.dist = clamp(dist0*pinch0/Math.max(Math.hypot(a[0]-b[0], a[1]-b[1]), 1), 0.3, 4);
      syncCam(); render(); return;
    }
    if (S.dragging && S.drag >= 0) {
      const w = V.unproject(F, cur[0], cur[1], dragDepth);
      const lim = S.scene*2.2;
      S.ctrl[S.drag] = { x:clamp(w.x,-lim,lim), y:clamp(w.y,-lim,lim), z:clamp(w.z,-lim,lim) };
      render(); return;
    }
    if (S.orbiting && last) {
      S.cam.az = ((S.cam.az - (cur[0]-last[0])*0.42 + 180) % 360 + 360) % 360 - 180;
      S.cam.el = clamp(S.cam.el + (cur[1]-last[1])*0.35, -V.EL_MAX, V.EL_MAX);
      last = cur; syncCam(); render(); return;
    }
    const h = hitHandle(cur[0], cur[1]);
    if (h !== S.hover) { S.hover = h; render(); }
    probe(cur[0], cur[1]);
  });
  function endPointer(e) {
    pointers.delete(e.pointerId);
    if (pointers.size < 2) pinch0 = 0;
    if (S.dragging) { S.dragging = false; S.drag = -1; solve(); }
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
    S.cam.dist = clamp(S.cam.dist*Math.exp(e.deltaY*0.0012), 0.3, 4);
    syncCam(); render();
  }, { passive:false });

  function probe(sx, sy) {
    if (!S.T || !F) return;
    const pl = planeOf();
    const far = V.unproject(F, sx, sy, F.focal);
    const dx = far.x-F.eye.x, dy = far.y-F.eye.y, dz = far.z-F.eye.z;
    const den = dx*pl.n.x + dy*pl.n.y + dz*pl.n.z;
    if (Math.abs(den) < 1e-9) return;
    const t = ((pl.o.x-F.eye.x)*pl.n.x + (pl.o.y-F.eye.y)*pl.n.y + (pl.o.z-F.eye.z)*pl.n.z)/den;
    if (t <= 0) return;
    const p = { x:F.eye.x+dx*t, y:F.eye.y+dy*t, z:F.eye.z+dz*t };
    const f = FB.fieldAt(S.T, p.x, p.y, p.z, S.t);
    const s = FB.poynting(f);
    $('ro-e').textContent = fmt(Math.hypot(f.ex,f.ey,f.ez), 'V/m');
    $('ro-b').textContent = fmt(Math.hypot(f.bx,f.by,f.bz), 'T');
    $('ro-s').textContent = fmt(Math.hypot(s.x,s.y,s.z), 'W/m²');
    $('ro-v').textContent = fmt(f.v, 'V');
  }

  /* ---------------- the clock ---------------- */
  let running = false;
  function tick() {
    if (running) return;
    running = true;
    let t0 = performance.now();
    const step = now => {
      if (!S.playing) { running = false; sample(); return; }
      // Capped generously rather than tightly: the cap is only there to stop a
      // backgrounded tab jumping the clock on return, and a full-quality frame
      // legitimately takes most of a second.
      const dtWall = Math.min((now - t0)/1000, 2); t0 = now;
      S.t += dtWall * nsPerSec() * NS;
      if (S.t >= S.tEnd) { S.t = S.tEnd; S.playing = false; $('btn-play').textContent = 'Play'; }
      sample();
      requestAnimationFrame(step);
    };
    requestAnimationFrame(step);
  }

  /* ---------------- go ---------------- */
  window.addEventListener('resize', () => resize());
  syncModes();
  $('v-off').textContent = '0.0 cm';
  $('v-tau').textContent = (S.tau/NS).toFixed(1) + ' ns';
  $('v-speed').textContent = nsPerSec().toFixed(1) + ' ns/s';
  $('v-asize').textContent = asizeText(S.arrowSize);
  $('v-ldens').textContent = S.lineDensity.toFixed(2) + '×';
  syncCam();
  const r0 = cv.getBoundingClientRect();
  W = Math.round(r0.width) || 900; H = Math.round(r0.height) || 700;
  solve();
  resize();
})();
