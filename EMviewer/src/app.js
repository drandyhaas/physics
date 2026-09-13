/*
 * app.js — canvas rendering and UI for the field bench.
 * All physics lives in solver.js; this file only decides how to draw it.
 */
(function () {
  'use strict';

  const { MU0, clamp, lerp } = FB;
  // Pixels per metre. Fitted to the stage in resize() so the whole circuit stays
  // framed whatever the screen is; 1 px = 1 mm is only the large-monitor ceiling.
  let PPM = 1000;
  const VIEW_HW = 0.40;   // half-width of the world box kept on screen, metres
  const VIEW_HH = 0.24;   // half-height of the same box (covers every preset + margin)
  const PPM_MAX = 1000;

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
    B: ramp([[0,hex('#2ee8cc')],[.32,hex('#126b68')],[.5,hex('#08191f')],[.68,hex('#7a2a58')],[1,hex('#ff6fb5')]])
  };

  /* ---------------- state ---------------- */
  const S = {
    ctrl: FB.PRESETS.rect.map(p => ({ x: p[0], y: p[1] })),
    mode: 'E',
    heat: false, arrows: true, lines: true,
    emf: 9, rLED: 100, rho: 0, radius: 0.003,
    drag: -1, hover: -1, dragging: false, quality: 'high',
    sol: null, grid: null, range: { E: null, B: null, S: null }
  };

  const $ = id => document.getElementById(id);
  const cv = $('cv');
  const ctx = cv.getContext('2d', { alpha: false });
  let W = 0, H = 0, DPR = 1;

  const toScreen = (x, y) => [W/2 + x*PPM, H/2 - y*PPM];
  const toWorld  = (sx, sy) => [(sx - W/2)/PPM, (H/2 - sy)/PPM];

  /* ---------------- heatmap ---------------- */
  const heatCanvas = document.createElement('canvas');
  function drawHeat(G, mode) {
    const { nx, ny } = G;
    heatCanvas.width = nx; heatCanvas.height = ny;
    const hctx = heatCanvas.getContext('2d');
    const img = hctx.createImageData(nx, ny), d = img.data;
    const rng = S.range[mode], lo = rng.lo, hi = rng.hi, lut = RAMP[mode];

    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const k = j*nx + i;
        let t;
        if (mode === 'B') {
          const b = G.B[k];
          const u = clamp((Math.log10(Math.abs(b) + 1e-18) - lo) / (hi - lo || 1), 0, 1);
          t = 0.5 + (b < 0 ? -1 : 1) * u * 0.5;
        } else {
          t = clamp((Math.log10(FB.magAt(G, mode, k) + 1e-18) - lo) / (hi - lo || 1), 0, 1);
        }
        const c = clamp(Math.round(t * 255), 0, 255) * 3;
        const p = ((ny - 1 - j) * nx + i) * 4;      // grid j=0 is world bottom
        d[p] = lut[c]; d[p+1] = lut[c+1]; d[p+2] = lut[c+2]; d[p+3] = 255;
      }
    }
    hctx.putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(heatCanvas, 0, 0, W, H);
  }

  /* ---------------- arrows ---------------- */
  function drawArrows(G, mode) {
    if (mode === 'B') return drawBGlyphs(G);
    const step = W < 560 ? 30 : 26;
    const rng = S.range[mode], out = [0, 0];
    const BUCKETS = 6;
    const paths = Array.from({ length: BUCKETS }, () => new Path2D());

    for (let sy = step/2; sy < H; sy += step) {
      for (let sx = step/2; sx < W; sx += step) {
        const [wx, wy] = toWorld(sx, sy);
        if (!FB.sampleVec(G, mode, wx, wy, out)) continue;
        const m = Math.hypot(out[0], out[1]);
        if (m < 1e-15) continue;
        const t = clamp((Math.log10(m + 1e-18) - rng.lo) / ((rng.hi - rng.lo) || 1), 0, 1);
        if (t < 0.04) continue;
        const len = step * (0.3 + 0.62 * t);
        const ux = out[0]/m, uy = -out[1]/m;               // screen y is flipped
        const hx = sx + ux*len/2, hy = sy + uy*len/2;
        const P = paths[Math.min(BUCKETS-1, Math.floor(t * BUCKETS))];
        P.moveTo(sx - ux*len/2, sy - uy*len/2); P.lineTo(hx, hy);
        const hs = clamp(len * 0.3, 2.5, 7);
        P.moveTo(hx, hy); P.lineTo(hx - ux*hs + uy*hs*0.55, hy - uy*hs - ux*hs*0.55);
        P.moveTo(hx, hy); P.lineTo(hx - ux*hs - uy*hs*0.55, hy - uy*hs + ux*hs*0.55);
      }
    }
    const lut = RAMP[mode];
    ctx.lineCap = 'round'; ctx.lineJoin = 'round';
    for (let b = 0; b < BUCKETS; b++) {
      const c = Math.round((b + 0.7) / BUCKETS * 255) * 3;
      ctx.strokeStyle = `rgba(${lut[c]},${lut[c+1]},${lut[c+2]},${0.55 + 0.4*b/BUCKETS})`;
      ctx.lineWidth = 1 + b * 0.22;
      ctx.stroke(paths[b]);
    }
  }

  // B is perpendicular to the screen, so it gets ⊙ / ⊗ marks rather than arrows.
  function drawBGlyphs(G) {
    const step = W < 560 ? 46 : 40, rng = S.range.B;
    for (let sy = step/2; sy < H; sy += step) {
      for (let sx = step/2; sx < W; sx += step) {
        const [wx, wy] = toWorld(sx, sy);
        const fi = (wx - G.x0)/G.dx, fj = (wy - G.y0)/G.dy;
        if (fi < 0 || fj < 0 || fi > G.nx-1 || fj > G.ny-1) continue;
        const b = G.B[Math.round(fj)*G.nx + Math.round(fi)];
        const t = clamp((Math.log10(Math.abs(b) + 1e-18) - rng.lo) / ((rng.hi - rng.lo) || 1), 0, 1);
        if (t < 0.12) continue;
        const r = 2 + t*6, col = b > 0 ? '#ff8ec6' : '#5ff0d8';
        ctx.strokeStyle = col; ctx.lineWidth = 1 + t; ctx.globalAlpha = 0.35 + 0.55*t;
        ctx.beginPath(); ctx.arc(sx, sy, r, 0, 7); ctx.stroke();
        if (b > 0) { ctx.beginPath(); ctx.arc(sx, sy, Math.max(0.8, r*0.22), 0, 7); ctx.fillStyle = col; ctx.fill(); }
        else { const q = r*0.68; ctx.beginPath(); ctx.moveTo(sx-q,sy-q); ctx.lineTo(sx+q,sy+q); ctx.moveTo(sx+q,sy-q); ctx.lineTo(sx-q,sy+q); ctx.stroke(); }
        ctx.globalAlpha = 1;
      }
    }
  }

  /* ---------------- evenly spaced field lines ---------------- */
  function drawLines(G, mode) {
    if (mode === 'B') return;
    const sep = W < 560 ? 26 : 22;
    const gw = Math.ceil(W/sep), gh = Math.ceil(H/sep);
    const occ = new Uint8Array(gw*gh);
    const out = [0, 0], stepLen = 3.2/PPM, MAX = 320;
    const lines = [];

    const cellFree = (sx, sy) => {
      const i = (sx/sep)|0, j = (sy/sep)|0;
      return i >= 0 && j >= 0 && i < gw && j < gh && !occ[j*gw + i];
    };
    const mark = (sx, sy) => {
      const i = (sx/sep)|0, j = (sy/sep)|0;
      if (i >= 0 && j >= 0 && i < gw && j < gh) occ[j*gw + i] = 1;
    };

    const seeds = [];
    for (let j = 0; j < gh; j++) for (let i = 0; i < gw; i++)
      seeds.push([(i + 0.5)*sep, (j + 0.5)*sep]);
    // interleave so seeding spreads out rather than sweeping row by row
    seeds.sort((a, b) => ((a[0]*13.7 + a[1]*7.3) % 1) - ((b[0]*13.7 + b[1]*7.3) % 1));

    for (const seed of seeds) {
      if (!cellFree(seed[0], seed[1])) continue;
      const pts = [];
      for (const dir of [1, -1]) {
        let [wx, wy] = toWorld(seed[0], seed[1]);
        const run = [];
        for (let n = 0; n < MAX; n++) {
          if (!FB.sampleVec(G, mode, wx, wy, out)) break;
          let m = Math.hypot(out[0], out[1]);
          if (m < 1e-15) break;
          const k1x = out[0]/m*dir, k1y = out[1]/m*dir;        // RK2 midpoint
          if (!FB.sampleVec(G, mode, wx + k1x*stepLen*0.5, wy + k1y*stepLen*0.5, out)) break;
          m = Math.hypot(out[0], out[1]);
          if (m < 1e-15) break;
          wx += out[0]/m*dir*stepLen; wy += out[1]/m*dir*stepLen;
          const sc = toScreen(wx, wy);
          if (sc[0] < -20 || sc[1] < -20 || sc[0] > W+20 || sc[1] > H+20) break;
          const fi = (wx - G.x0)/G.dx, fj = (wy - G.y0)/G.dy;  // stop at the conductor
          if (fi >= 0 && fj >= 0 && fi < G.nx-1 && fj < G.ny-1 &&
              G.D[Math.round(fj)*G.nx + Math.round(fi)] < S.radius*1.05) break;
          if (n > 3 && !cellFree(sc[0], sc[1])) break;
          run.push(sc);
        }
        if (dir === 1) pts.push(...run); else pts.unshift(...run.reverse());
      }
      if (pts.length < 6) continue;
      for (const p of pts) mark(p[0], p[1]);
      lines.push(pts);
    }

    ctx.lineCap = 'round'; ctx.lineJoin = 'round';
    ctx.strokeStyle = mode === 'S' ? 'rgba(255,222,160,0.30)' : 'rgba(190,238,240,0.26)';
    ctx.lineWidth = 1.05;
    ctx.beginPath();
    for (const pts of lines) {
      ctx.moveTo(pts[0][0], pts[0][1]);
      for (let i = 1; i < pts.length; i++) ctx.lineTo(pts[i][0], pts[i][1]);
    }
    ctx.stroke();
  }

  /* ---------------- the circuit ---------------- */
  function drawCircuit(sol) {
    const pts = sol.outline, N = sol.N;
    const rpx = Math.max(2.2, S.radius * PPM);

    ctx.beginPath();
    let p = toScreen(pts[0], pts[1]);
    ctx.moveTo(p[0], p[1]);
    for (let i = 1; i < N; i++) { p = toScreen(pts[i*2], pts[i*2+1]); ctx.lineTo(p[0], p[1]); }
    ctx.closePath();
    ctx.lineJoin = 'round';
    ctx.strokeStyle = 'rgba(3,12,16,0.92)'; ctx.lineWidth = rpx*2 + 5; ctx.stroke();
    ctx.strokeStyle = '#c6d7da';            ctx.lineWidth = rpx*2;     ctx.stroke();

    const runOf = t => { const a = []; for (let i = 0; i < N; i++) if (sol.type[i] === t) a.push(i); return a; };
    const led = runOf(FB.LED), bat = runOf(FB.BATTERY);
    const P = sol.I * sol.I * S.rLED;

    const strokeRun = (run, w, col, cap) => {
      if (!run.length) return;
      ctx.beginPath();
      let a = toScreen(sol.ax[run[0]], sol.ay[run[0]]);
      ctx.moveTo(a[0], a[1]);
      for (const i of run) { const b = toScreen(sol.bx[i], sol.by[i]); ctx.lineTo(b[0], b[1]); }
      ctx.lineCap = cap || 'round'; ctx.lineWidth = w; ctx.strokeStyle = col; ctx.stroke();
    };

    if (led.length) {
      const mid = led[Math.floor(led.length/2)];
      const c = toScreen(sol.cx[mid], sol.cy[mid]);
      const bright = clamp(Math.sqrt(P/0.25), 0, 1);
      const R = 42 + 70*bright;
      const gr = ctx.createRadialGradient(c[0], c[1], 0, c[0], c[1], R);
      gr.addColorStop(0, `rgba(255,225,170,${0.36*bright + 0.05})`);
      gr.addColorStop(1, 'rgba(255,200,120,0)');
      ctx.fillStyle = gr; ctx.beginPath(); ctx.arc(c[0], c[1], R, 0, 7); ctx.fill();
      strokeRun(led, rpx*2 + 9, '#2a1608');
      strokeRun(led, rpx*2 + 3, `rgb(255,${Math.round(170+70*bright)},${Math.round(70+90*bright)})`);
    }

    if (bat.length) {
      strokeRun(bat, rpx*2 + 11, '#071b21', 'butt');
      strokeRun(bat, rpx*2 + 6, '#7fe6d2', 'butt');
      const last = bat[bat.length-1], first = bat[0], gm = bat[Math.floor(bat.length/2)];
      const nx = -sol.ty[gm], ny = sol.tx[gm];
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.font = '600 13px "IBM Plex Sans", sans-serif'; ctx.fillStyle = '#bfeee4';
      let q = toScreen(sol.bx[last] + nx*0.016, sol.by[last] + ny*0.016); ctx.fillText('+', q[0], q[1]);
      q = toScreen(sol.ax[first] + nx*0.016, sol.ay[first] + ny*0.016);   ctx.fillText('−', q[0], q[1]);
      q = toScreen(sol.cx[gm] - nx*0.022, sol.cy[gm] - ny*0.022);
      ctx.font = '500 11px "IBM Plex Sans", sans-serif'; ctx.fillStyle = 'rgba(150,220,210,.8)';
      ctx.fillText(S.emf.toFixed(1) + ' V', q[0], q[1]);
    }

    // current-direction ticks
    ctx.fillStyle = 'rgba(255,215,150,.85)';
    const TICKS = 26;
    for (let n = 0; n < TICKS; n++) {
      const i = Math.floor(n/TICKS*N);
      const c = toScreen(sol.cx[i], sol.cy[i]);
      const tx = sol.tx[i], ty = sol.ty[i], s = 3.4;
      ctx.beginPath();
      ctx.moveTo(c[0] + tx*s, c[1] - ty*s);
      ctx.lineTo(c[0] - tx*s*0.5 - ty*s*0.6, c[1] + ty*s*0.5 - tx*s*0.6);
      ctx.lineTo(c[0] - tx*s*0.5 + ty*s*0.6, c[1] + ty*s*0.5 + tx*s*0.6);
      ctx.closePath(); ctx.fill();
    }
  }

  function drawHandles() {
    S.ctrl.forEach((p, i) => {
      const c = toScreen(p.x, p.y);
      const act = i === S.hover || i === S.drag;
      ctx.beginPath(); ctx.arc(c[0], c[1], act ? 8 : 6, 0, 7);
      ctx.fillStyle = act ? '#ffffff' : 'rgba(240,252,252,.78)'; ctx.fill();
      ctx.lineWidth = 2; ctx.strokeStyle = 'rgba(5,20,26,.85)'; ctx.stroke();
    });
  }

  /* ---------------- render ---------------- */
  function updateRange(G, mode) {
    const n = G.nx * G.ny, vals = new Float64Array(n);
    for (let k = 0; k < n; k++) vals[k] = Math.log10(Math.abs(FB.magAt(G, mode, k)) + 1e-18);
    const sorted = Array.from(vals).sort((a, b) => a - b);
    const lo = sorted[Math.floor(n*0.22)], hi = sorted[Math.floor(n*0.985)];
    const prev = S.range[mode];
    const nl = prev ? lerp(prev.lo, lo, 0.5) : lo;
    const nh = prev ? lerp(prev.hi, hi, 0.5) : hi;
    S.range[mode] = { lo: nl, hi: Math.max(nh, nl + 0.4) };
  }

  function drawLegend(mode) {
    const rng = S.range[mode]; if (!rng) return;
    const lut = RAMP[mode], stops = [];
    for (let i = 0; i <= 8; i++) {
      const c = Math.round(i/8*255)*3;
      stops.push(`rgb(${lut[c]},${lut[c+1]},${lut[c+2]}) ${i/8*100}%`);
    }
    $('leg-bar').style.background = `linear-gradient(90deg, ${stops.join(',')})`;
    if (mode === 'B') {
      $('leg-cap').textContent = 'B out of the screen';
      $('leg-lo').textContent = '−' + fmt(Math.pow(10, rng.hi), 'T');
      $('leg-hi').textContent = '+' + fmt(Math.pow(10, rng.hi), 'T');
    } else {
      $('leg-cap').textContent = mode === 'S' ? 'Energy flow |S|' : 'Field strength |E|';
      const u = mode === 'S' ? 'W/m²' : 'V/m';
      $('leg-lo').textContent = fmt(Math.pow(10, rng.lo), u);
      $('leg-hi').textContent = fmt(Math.pow(10, rng.hi), u);
    }
  }

  function render() {
    ctx.fillStyle = '#05141a'; ctx.fillRect(0, 0, W, H);
    const sol = S.sol, G = S.grid;
    if (!sol || !G) return;
    if (S.heat) drawHeat(G, S.mode);
    if (S.lines && S.quality === 'high') drawLines(G, S.mode);
    if (S.arrows) drawArrows(G, S.mode);
    drawCircuit(sol);
    drawHandles();
    drawLegend(S.mode);
  }

  function updateFacts(sol) {
    $('f-I').textContent  = fmt(sol.I, 'A');
    $('f-P').textContent  = fmt(sol.I*sol.I*S.rLED, 'W');
    $('f-L').textContent  = (sol.L*100).toFixed(1) + ' cm';
    $('f-Ew').textContent = fmt(sol.I*S.rho, 'V/m');
    $('f-El').textContent = fmt(sol.I*S.rLED/sol.ledLen, 'V/m');
    $('f-Eb').textContent = fmt(-S.emf/sol.batLen, 'V/m');
    $('f-q').textContent  = fmt(sol.peak, 'C/m');
    $('nh').textContent   = S.ctrl.length;
  }

  /* ---------------- recompute pipeline ---------------- */
  let pending = false;
  function recompute() {
    if (pending) return;
    pending = true;
    requestAnimationFrame(() => {
      pending = false;
      const hi = S.quality === 'high';
      const sol = FB.solveCircuit({
        ctrl: S.ctrl, N: hi ? 260 : 140, per: hi ? 60 : 34,
        radius: S.radius, emf: S.emf, rLED: S.rLED, rho: S.rho
      });
      if (!sol) { ctx.fillStyle = '#05141a'; ctx.fillRect(0, 0, W, H); return; }
      S.sol = sol;
      S.grid = FB.buildGrid(sol, {
        w: W, h: H, ppm: PPM,
        cell: hi ? (W > 900 ? 5 : 6) : 12,
        needE: S.mode !== 'B',
        needB: S.mode !== 'E'
      });
      updateRange(S.grid, S.mode);
      render();
      updateFacts(sol);
    });
  }
  function settle() {
    if (S.quality === 'high') return;
    S.quality = 'high';
    recompute();
  }

  /* ---------------- layout ---------------- */
  function resize() {
    const r = cv.getBoundingClientRect();
    if (r.width < 1 || r.height < 1) return;   // stage not laid out yet
    DPR = Math.min(window.devicePixelRatio || 1, 2);
    W = Math.round(r.width);
    H = Math.round(r.height);
    // Contain-fit the world box: the narrower axis sets the scale, so a tall phone,
    // a squat landscape window and a wide monitor all show the whole circuit.
    PPM = Math.min(PPM_MAX, W / (2*VIEW_HW), H / (2*VIEW_HH));
    cv.width = Math.round(W*DPR); cv.height = Math.round(H*DPR);
    ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
    recompute();
  }

  /* ---------------- UI ---------------- */
  const MODE_NOTE = {
    E: 'The electric field of the surface charge the circuit builds up. Along ideal wire it meets the metal at right angles — no tangential push, no work done. At the LED it tips over and dives in.',
    B: 'The magnetic field of the current. In the plane of the circuit it points straight in and out of the screen: ⊗ away from you, ⊙ towards you.',
    S: 'Energy flow, S = E×B/μ₀. It leaves the battery sideways into empty space, glides along outside the wire, and turns inward wherever power is being dissipated. The wire guides the energy; it does not carry it.'
  };
  function setMode(m) {
    S.mode = m;
    for (const k of ['E', 'B', 'S']) $('m-' + k).setAttribute('aria-pressed', String(k === m));
    $('modenote').textContent = MODE_NOTE[m];
    recompute();                      // the grid may need a component it skipped
  }
  ['E', 'B', 'S'].forEach(k => $('m-' + k).addEventListener('click', () => setMode(k)));

  function bindSlider(id, valId, apply, show) {
    const el = $(id);
    const upd = () => {
      const v = parseFloat(el.value);
      apply(v);
      $(valId).textContent = show(v);
      S.quality = 'low'; recompute();
      clearTimeout(el._t); el._t = setTimeout(settle, 180);
    };
    el.addEventListener('input', upd);
    $(valId).textContent = show(parseFloat(el.value));
  }
  bindSlider('s-emf', 'v-emf', v => S.emf = v,           v => v.toFixed(2) + ' V');
  bindSlider('s-led', 'v-led', v => S.rLED = v,          v => v.toFixed(0) + ' Ω');
  bindSlider('s-rho', 'v-rho', v => S.rho = v,           v => v.toFixed(0) + ' Ω/m');
  bindSlider('s-rad', 'v-rad', v => S.radius = v/1000,   v => v.toFixed(1) + ' mm');

  [['c-heat','heat'], ['c-arrows','arrows'], ['c-lines','lines']].forEach(([id, key]) => {
    $(id).addEventListener('change', e => { S[key] = e.target.checked; render(); });
  });

  document.querySelectorAll('[data-preset]').forEach(b => {
    b.addEventListener('click', () => {
      S.ctrl = FB.PRESETS[b.dataset.preset].map(p => ({ x: p[0], y: p[1] }));
      recompute();
    });
  });

  $('railtoggle').addEventListener('click', () => {
    const open = $('rail').classList.toggle('open');
    $('railtoggle').setAttribute('aria-expanded', String(open));
    $('railtoggle').textContent = open ? 'Hide' : 'Controls';
  });

  /* ---------------- pointer ---------------- */
  function hitHandle(sx, sy) {
    for (let i = 0; i < S.ctrl.length; i++) {
      const c = toScreen(S.ctrl[i].x, S.ctrl[i].y);
      if (Math.hypot(c[0] - sx, c[1] - sy) < 16) return i;
    }
    return -1;
  }
  const local = e => {
    const r = cv.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top];
  };

  let lpTimer = null;
  cv.addEventListener('pointerdown', e => {
    const [sx, sy] = local(e);
    const i = hitHandle(sx, sy);
    if (i < 0) return;
    if (e.altKey && S.ctrl.length > 4) { S.ctrl.splice(i, 1); recompute(); return; }
    S.drag = i; S.dragging = true; S.quality = 'low';
    cv.setPointerCapture(e.pointerId);
    lpTimer = setTimeout(() => {
      if (S.drag === i && S.ctrl.length > 4) {
        S.ctrl.splice(i, 1); S.drag = -1; S.dragging = false; recompute();
      }
    }, 650);
  });

  cv.addEventListener('pointermove', e => {
    const [sx, sy] = local(e);
    if (S.dragging && S.drag >= 0) {
      clearTimeout(lpTimer);
      const [wx, wy] = toWorld(sx, sy);
      S.ctrl[S.drag].x = clamp(wx, -W/2/PPM + 0.01, W/2/PPM - 0.01);
      S.ctrl[S.drag].y = clamp(wy, -H/2/PPM + 0.01, H/2/PPM - 0.01);
      recompute();
      return;
    }
    const h = hitHandle(sx, sy);
    if (h !== S.hover) { S.hover = h; render(); }
    cv.style.cursor = h >= 0 ? 'grab' : 'crosshair';
    updateReadout(sx, sy);
  });

  function endDrag() {
    clearTimeout(lpTimer);
    if (S.dragging) { S.dragging = false; S.drag = -1; settle(); }
  }
  cv.addEventListener('pointerup', endDrag);
  cv.addEventListener('pointercancel', endDrag);
  cv.addEventListener('pointerleave', () => { $('rhead').textContent = 'Point at the field'; });

  cv.addEventListener('dblclick', e => {
    const [sx, sy] = local(e);
    const [wx, wy] = toWorld(sx, sy);
    let best = 0, bd = Infinity;
    for (let i = 0; i < S.ctrl.length; i++) {
      const a = S.ctrl[i], b = S.ctrl[(i + 1) % S.ctrl.length];
      const d = Math.hypot((a.x + b.x)/2 - wx, (a.y + b.y)/2 - wy);
      if (d < bd) { bd = d; best = i; }
    }
    S.ctrl.splice(best + 1, 0, { x: wx, y: wy });
    recompute();
  });

  function updateReadout(sx, sy) {
    if (!S.sol) return;
    const [wx, wy] = toWorld(sx, sy);
    const f = FB.fieldAt(S.sol, wx, wy);
    const Sv = FB.poynting(f);
    const t = S.sol.type[f.near];
    $('rhead').textContent = f.d < S.radius
      ? (t === FB.WIRE ? 'Inside the plain wire' : t === FB.LED ? 'Inside the LED' : 'Inside the battery')
      : (f.d < S.radius*3 ? 'Just outside the conductor' : 'In the space around the circuit');
    $('ro-e').textContent = fmt(Math.hypot(f.ex, f.ey), 'V/m');
    $('ro-b').textContent = fmt(f.bz, 'T');
    $('ro-s').textContent = fmt(Math.hypot(Sv[0], Sv[1]), 'W/m²');
    $('ro-v').textContent = fmt(f.v - S.sol.V0, 'V');
  }

  /* ---------------- go ---------------- */
  let rt = null;
  window.addEventListener('resize', () => {
    S.quality = 'low'; resize();
    clearTimeout(rt); rt = setTimeout(settle, 200);
  });
  $('modenote').textContent = MODE_NOTE.E;
  resize();
})();
