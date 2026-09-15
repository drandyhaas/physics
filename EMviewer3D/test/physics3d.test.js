/*
 * Physics regression tests for the 3D field bench. No dependencies.
 *   node test/physics3d.test.js
 *
 * These check the solver against things known independently: Ohm's law, the
 * infinite-wire field, the on-axis field of a current loop, a brute-force
 * Biot-Savart integral, and -- when the planar field bench is present as a
 * sibling -- the 2D solver it generalises.
 */
const FB = require('../src/solver3d.js');
const { MU0 } = FB;

let failures = 0;
function check(name, cond, extra) {
  if (!cond) failures++;
  console.log((cond ? '  ok   ' : '  FAIL ') + name + (extra ? '   [' + extra + ']' : ''));
}
function section(t) { console.log('\n' + t); }
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);

const base = {
  ctrl: FB.PRESETS.rect.map(p => ({ x: p[0], y: p[1], z: p[2] })),
  N: 260, per: 60, radius: 0.003, emf: 3, rLED: 100, rho: 0
};
const build = over => FB.solveCircuit(Object.assign({}, base, over));
const ctrlOf = pts => pts.map(p => ({ x: p[0], y: p[1], z: p[2] }));

/* ------------------------------------------------------------------ */
section('Bookkeeping');
let spread = 0;
let sol = build({});
check('solve succeeds', !!sol);
check('loop closes in potential', Math.abs(sol.closure) < 1e-12, 'dV=' + sol.closure.toExponential(2));
check('current is emf/R', Math.abs(sol.I - 3 / 100) < 1e-12, 'I=' + sol.I);
let Q = 0; for (let i = 0; i < sol.N; i++) Q += sol.lam[i] * sol.len[i];
check('circuit is charge neutral', Math.abs(Q) < 1e-24, 'Q=' + Q.toExponential(2));
check('segments are near enough equal length', (() => {
  let mn = Infinity, mx = 0;
  for (let i = 0; i < sol.N; i++) { mn = Math.min(mn, sol.len[i]); mx = Math.max(mx, sol.len[i]); }
  spread = (mx - mn) / mx;
  return spread < 0.02;   // chords fall a little short of the arc where it bends
})(), 'spread ' + (100 * spread).toFixed(3) + '%');
check('tangents are unit', (() => {
  for (let i = 0; i < sol.N; i++) {
    const n = Math.hypot(sol.tx[i], sol.ty[i], sol.tz[i]);
    if (Math.abs(n - 1) > 1e-12) return false;
  }
  return true;
})());

/* ------------------------------------------------------------------ */
section('Biot-Savart against brute-force integration');
// Subdivide every segment into K point elements and sum. This tests the exact
// finite-filament formula in fieldAt with no shared code between the two.
function bruteB(sol, p, K) {
  const k = MU0 * sol.I / (4 * Math.PI);
  let Bx = 0, By = 0, Bz = 0;
  for (let i = 0; i < sol.N; i++) {
    const Lx = sol.bx[i] - sol.ax[i], Ly = sol.by[i] - sol.ay[i], Lz = sol.bz[i] - sol.az[i];
    for (let j = 0; j < K; j++) {
      const t = (j + 0.5) / K;
      const qx = sol.ax[i] + Lx * t, qy = sol.ay[i] + Ly * t, qz = sol.az[i] + Lz * t;
      const dx = p[0] - qx, dy = p[1] - qy, dz = p[2] - qz;
      const r2 = dx * dx + dy * dy + dz * dz, r3 = r2 * Math.sqrt(r2);
      const ex = Lx / K, ey = Ly / K, ez = Lz / K;
      Bx += k * (ey * dz - ez * dy) / r3;
      By += k * (ez * dx - ex * dz) / r3;
      Bz += k * (ex * dy - ey * dx) / r3;
    }
  }
  return [Bx, By, Bz];
}
for (const [label, p] of [['near the wire', [0.26, 0.02, 0.012]],
                          ['off the plane', [0.1, 0.05, 0.09]],
                          ['far away',      [0.6, 0.4, 0.5]]]) {
  const f = FB.fieldAt(sol, p[0], p[1], p[2]);
  const b = bruteB(sol, p, 60);
  const mine = Math.hypot(f.bx, f.by, f.bz), ref = Math.hypot(b[0], b[1], b[2]);
  const cosang = (f.bx * b[0] + f.by * b[1] + f.bz * b[2]) / (mine * ref);
  check('|B| matches brute force ' + label, rel(mine, ref) < 2e-3,
        'mine=' + mine.toExponential(3) + ' ref=' + ref.toExponential(3));
  check('B direction matches ' + label, cosang > 0.9999, 'cos=' + cosang.toFixed(6));
}

/* ------------------------------------------------------------------ */
section('Analytic limits');
// Long thin loop: close to a straight run, so B -> mu0 I / 2 pi d beside it.
const longLoop = build({
  ctrl: ctrlOf([[-.6,-.02,0],[-.2,-.02,0],[.2,-.02,0],[.6,-.02,0],
                [.6,.02,0],[.2,.02,0],[-.2,.02,0],[-.6,.02,0]]),
  N: 400
});
for (const d of [0.05, 0.1]) {
  // beside the middle of the lower run, far from its ends
  const f = FB.fieldAt(longLoop, 0, -0.02 - d, 0);
  const got = Math.hypot(f.bx, f.by, f.bz);
  // the return run 0.04 m beyond carries current the other way, so it subtracts
  const want = MU0 * longLoop.I / (2 * Math.PI) * (1 / d - 1 / (d + 0.04));
  check('B beside a long run at d=' + d + ' m', rel(got, want) < 0.06,
        'got=' + got.toExponential(3) + ' want=' + want.toExponential(3) +
        ' (' + (100 * rel(got, want)).toFixed(1) + '%)');
}
// Circular loop: B at the centre is mu0 I / 2R on the axis.
const ring = build({ ctrl: ctrlOf(FB.PRESETS.circle), N: 400 });
{
  const f = FB.fieldAt(ring, 0, 0, 0);
  const got = Math.hypot(f.bx, f.by, f.bz);
  const want = MU0 * ring.I / (2 * 0.2);
  check('B at the centre of a ring is mu0 I / 2R', rel(got, want) < 0.02,
        'got=' + got.toExponential(3) + ' want=' + want.toExponential(3) +
        ' (' + (100 * rel(got, want)).toFixed(2) + '%)');
  check('that field is along the ring axis',
        Math.abs(f.bz) / got > 0.999, '|Bz|/|B|=' + (Math.abs(f.bz) / got).toFixed(6));
}
// On-axis: B(z) = mu0 I R^2 / 2 (R^2+z^2)^{3/2}
for (const z of [0.1, 0.3]) {
  const f = FB.fieldAt(ring, 0, 0, z);
  const got = Math.hypot(f.bx, f.by, f.bz);
  const R = 0.2, want = MU0 * ring.I * R * R / (2 * Math.pow(R * R + z * z, 1.5));
  check('B on the ring axis at z=' + z + ' m', rel(got, want) < 0.02,
        'got=' + got.toExponential(3) + ' want=' + want.toExponential(3) +
        ' (' + (100 * rel(got, want)).toFixed(2) + '%)');
}

/* ------------------------------------------------------------------ */
section('The field at an ideal wire');
// Just outside resistanceless wire, E must have no tangential component.
function probeOutside(sol, frac, off) {
  const i = Math.floor(sol.N * frac);
  // a normal: away from the loop centroid, with the tangent projected out
  let gx = 0, gy = 0, gz = 0;
  for (let k = 0; k < sol.N; k++) { gx += sol.cx[k]; gy += sol.cy[k]; gz += sol.cz[k]; }
  gx /= sol.N; gy /= sol.N; gz /= sol.N;
  let nx = sol.cx[i] - gx, ny = sol.cy[i] - gy, nz = sol.cz[i] - gz;
  const d = nx * sol.tx[i] + ny * sol.ty[i] + nz * sol.tz[i];
  nx -= d * sol.tx[i]; ny -= d * sol.ty[i]; nz -= d * sol.tz[i];
  const n = Math.hypot(nx, ny, nz); nx /= n; ny /= n; nz /= n;
  const f = FB.fieldAt(sol, sol.cx[i] + nx * off, sol.cy[i] + ny * off, sol.cz[i] + nz * off);
  const mag = Math.hypot(f.ex, f.ey, f.ez);
  return { i, f, mag, Et: f.ex * sol.tx[i] + f.ey * sol.ty[i] + f.ez * sol.tz[i] };
}
// Probe just clear of the interior blend zone (1.35a), so what is measured is the
// Coulomb solution rather than the imposed interior field.
const OFF = 0.0045;
{
  const p = probeOutside(sol, 0.62, OFF);   // a plain-wire stretch
  check('E meets ideal wire at a right angle',
        Math.abs(p.Et) / p.mag < 0.01, '|Et|/|E|=' + (Math.abs(p.Et) / p.mag).toFixed(4));
}
{
  const r = build({ rho: 60 });
  const p = probeOutside(r, 0.62, OFF);
  const ideal = probeOutside(sol, 0.62, OFF);
  // E_t is continuous across the surface, so just outside the metal it is still rho*I
  check('resistive wire tilts the field over', rel(Math.abs(p.Et), 60 * r.I) < 0.05,
        'Et=' + p.Et.toExponential(3) + ' rhoI=' + (60 * r.I).toExponential(3));
  check('an ideal wire does not',
        Math.abs(p.Et) / p.mag > 5 * Math.abs(ideal.Et) / ideal.mag,
        'tangential fraction ' + (Math.abs(ideal.Et) / ideal.mag).toFixed(4) +
        ' ideal vs ' + (Math.abs(p.Et) / p.mag).toFixed(4) + ' resistive');
  // inside the metal Ohm's law should hold: E_t = rho * I
  const i = Math.floor(r.N * 0.62);
  const f = FB.fieldAt(r, r.cx[i], r.cy[i], r.cz[i]);
  const Et = f.ex * r.tx[i] + f.ey * r.ty[i] + f.ez * r.tz[i];
  check('inside the metal E_t = rho I', rel(Math.abs(Et), 60 * r.I) < 0.02,
        'Et=' + Et.toExponential(3) + ' rhoI=' + (60 * r.I).toExponential(3));
}

/* ------------------------------------------------------------------ */
section('Energy flow');
{
  // S should point inward at the LED: the load is where energy arrives.
  const s9 = build({ emf: 9 });
  let i = -1; for (let k = 0; k < s9.N; k++) if (s9.type[k] === FB.LED) { i = k; break; }
  const mid = i + 3;
  let gx = 0, gy = 0; for (let k = 0; k < s9.N; k++) { gx += s9.cx[k]; gy += s9.cy[k]; }
  gx /= s9.N; gy /= s9.N;
  let nx = s9.cx[mid] - gx, ny = s9.cy[mid] - gy;
  const d = nx * s9.tx[mid] + ny * s9.ty[mid];
  nx -= d * s9.tx[mid]; ny -= d * s9.ty[mid];
  const n = Math.hypot(nx, ny); nx /= n; ny /= n;
  const off = 0.01;
  const f = FB.fieldAt(s9, s9.cx[mid] + nx * off, s9.cy[mid] + ny * off, s9.cz[mid]);
  const S = FB.poynting(f);
  check('S points into the LED', S.x * nx + S.y * ny < 0,
        'S.n=' + (S.x * nx + S.y * ny).toExponential(2) + ' W/m2');
}

/* ------------------------------------------------------------------ */
section('Non-planar circuits');
for (const name of ['tilt', 'helix', 'saddle', 'circle']) {
  const s = build({ ctrl: ctrlOf(FB.PRESETS[name]) });
  let q = 0; for (let i = 0; i < s.N; i++) q += s.lam[i] * s.len[i];
  check(name.padEnd(7) + ' solves and stays neutral',
        !!s && Math.abs(q) < 1e-23 && isFinite(s.peak),
        'L=' + (100 * s.L).toFixed(1) + 'cm  peak=' + s.peak.toExponential(2) + ' C/m');
}
{
  // The solenoid should put a strong axial field down its middle.
  const h = build({ ctrl: ctrlOf(FB.PRESETS.helix), emf: 9 });
  const inside = FB.fieldAt(h, 0, 0, 0);
  const outside = FB.fieldAt(h, 0, 0.34, 0);
  const bi = Math.hypot(inside.bx, inside.by, inside.bz);
  const bo = Math.hypot(outside.bx, outside.by, outside.bz);
  check('solenoid field is axial inside', Math.abs(inside.bz) / bi > 0.9,
        '|Bz|/|B|=' + (Math.abs(inside.bz) / bi).toFixed(3));
  check('solenoid field is stronger inside than out', bi > 8 * bo,
        'inside=' + bi.toExponential(2) + ' outside=' + bo.toExponential(2) +
        ' ratio=' + (bi / bo).toFixed(1));
}
{
  // A loop tilted out of plane is the same physics rotated: peak charge density
  // and loop length must be unchanged from the flat rectangle.
  const flat = build({});
  const tip = build({ ctrl: ctrlOf(FB.PRESETS.tilt) });
  check('tilting a loop changes no scalar', rel(tip.peak, flat.peak) < 1e-9 && rel(tip.L, flat.L) < 1e-9,
        'peak ' + flat.peak.toExponential(4) + ' vs ' + tip.peak.toExponential(4));
}

/* ------------------------------------------------------------------ */
section('Convergence in segment count');
{
  const a = build({ N: 140 }), b = build({ N: 380 });
  const fa = FB.fieldAt(a, 0.1, 0.3, 0.05), fb = FB.fieldAt(b, 0.1, 0.3, 0.05);
  const ma = Math.hypot(fa.ex, fa.ey, fa.ez), mb = Math.hypot(fb.ex, fb.ey, fb.ez);
  check('E away from the wire is resolution independent', rel(ma, mb) < 0.02,
        'N=140 ' + ma.toFixed(2) + ' vs N=380 ' + mb.toFixed(2) + ' V/m');
  const ba = Math.hypot(fa.bx, fa.by, fa.bz), bb = Math.hypot(fb.bx, fb.by, fb.bz);
  check('B away from the wire is resolution independent', rel(ba, bb) < 0.02,
        'N=140 ' + ba.toExponential(3) + ' vs N=380 ' + bb.toExponential(3) + ' T');
}

/* ------------------------------------------------------------------ */
// The planar bench is a sibling project, not a dependency. When it is there,
// check that this solver reduces to it for a circuit that lies in a plane.
section('Reduction to the planar solver');
let FB2 = null;
try { FB2 = require('../../EMviewer/src/solver.js'); } catch (e) { /* standalone */ }
if (!FB2) {
  console.log('  skip  planar field bench not present alongside');
} else {
  const o = { N: 260, per: 60, radius: 0.003, emf: 3, rLED: 100, rho: 0 };
  const s2 = FB2.solveCircuit(Object.assign({ ctrl: FB2.PRESETS.rect.map(p => ({ x: p[0], y: p[1] })) }, o));
  const s3 = FB.solveCircuit(Object.assign({ ctrl: ctrlOf(FB.PRESETS.rect) }, o));
  let worstLam = 0;
  for (let i = 0; i < s2.N; i++) worstLam = Math.max(worstLam, rel(s3.lam[i], s2.lam[i]));
  check('charge density agrees with the 2D solve', worstLam < 1e-9,
        'worst relative difference ' + worstLam.toExponential(2));
  let worstE = 0, worstB = 0;
  for (const p of [[0.1, 0.3], [0.35, 0.05], [0, 0], [0.27, 0.16]]) {
    const f2 = FB2.fieldAt(s2, p[0], p[1]);
    const f3 = FB.fieldAt(s3, p[0], p[1], 0);
    worstE = Math.max(worstE, rel(Math.hypot(f3.ex, f3.ey), Math.hypot(f2.ex, f2.ey)));
    worstB = Math.max(worstB, rel(Math.abs(f3.bz), Math.abs(f2.bz)));
    // in the plane of a planar circuit the in-plane B components must vanish
    const leak = Math.hypot(f3.bx, f3.by) / Math.abs(f3.bz);
    if (leak > 1e-9) { failures++; console.log('  FAIL in-plane B leaks at ' + p); }
  }
  check('E agrees with the 2D solve', worstE < 1e-9, 'worst ' + worstE.toExponential(2));
  check('B_z agrees with the 2D solve', worstB < 1e-3, 'worst ' + worstB.toExponential(2));
  check('B has no in-plane part for a planar circuit', true);
}

/* ------------------------------------------------------------------ */
section('Field-line coverage: a random point has a line of its own nearby');
/*
 * The lines ARE the picture, and their failure mode is not being wrong but
 * being absent: a region draws nothing and reads as no field rather than as no
 * line. The seeds fill a volume; the region a loop encloses is a sheet through
 * that volume with no volume of its own, so the lattice landed in it only by
 * accident -- and with GRID even it never landed in it at all, so the E lines
 * running from the + charge across to the - charge were simply not drawn.
 *
 * Nothing about the field was wrong, so no field test could catch it. So audit
 * the drawing from the other end: throw random points into the scene and ask,
 * of each, how far it is to the nearest drawn line of that same field. Judge it
 * against the region OUTSIDE the loop, which was never in doubt, with both
 * populations held to the same distance from the wire so the comparison is of
 * seeding and not of falloff.
 */
{
  // A deterministic stream, so a failure is the same failure next run.
  const rng = seed => { let s = seed >>> 0;
    return () => { s = (s*1664525 + 1013904223) >>> 0; return s / 4294967296; }; };

  // Van Oosterom-Strackee, over a fan of the wire polygon: the solid angle the
  // loop subtends at a point. It is +-2pi on the surface the loop spans and
  // falls to zero far away, so |Omega| > pi says "inside the loop" without any
  // reference to a plane -- which the tilted and non-planar presets lack.
  const solidAngle = (sol, p) => {
    const a0 = sol.cx[0]-p.x, a1 = sol.cy[0]-p.y, a2 = sol.cz[0]-p.z;
    const na = Math.hypot(a0, a1, a2);
    let om = 0;
    for (let i = 1; i < sol.N-1; i++) {
      const b0 = sol.cx[i]-p.x,   b1 = sol.cy[i]-p.y,   b2 = sol.cz[i]-p.z;
      const c0 = sol.cx[i+1]-p.x, c1 = sol.cy[i+1]-p.y, c2 = sol.cz[i+1]-p.z;
      const nb = Math.hypot(b0,b1,b2), nc = Math.hypot(c0,c1,c2);
      const det = a0*(b1*c2-b2*c1) - a1*(b0*c2-b2*c0) + a2*(b0*c1-b1*c0);
      const ab = a0*b0+a1*b1+a2*b2, ac = a0*c0+a1*c1+a2*c2, bc = b0*c0+b1*c1+b2*c2;
      om += 2*Math.atan2(det, na*nb*nc + ab*nc + ac*nb + bc*na);
    }
    return om;
  };
  const median = a => { const b = a.slice().sort((x,y) => x-y); return b[b.length >> 1]; };

  for (const preset of ['rect', 'circle', 'tilt']) {
    const sol = build({ ctrl: ctrlOf(FB.PRESETS[preset]) });
    let R = 0.05;
    for (let k = 0; k < 3; k++) R = Math.max(R, Math.abs(sol.lo[k]), Math.abs(sol.hi[k]));

    // Points in a ball around the loop, held to a band of distances from the
    // wire: near enough that a line ought to pass, far enough to be out of the
    // metal. Inside and outside are then being asked the same question.
    const pts = [];
    const rnd = rng(20240917);
    while (pts.length < 1500) {
      const x = (2*rnd()-1)*1.15*R, y = (2*rnd()-1)*1.15*R, z = (2*rnd()-1)*1.15*R;
      if (Math.hypot(x, y, z) > 1.15*R) continue;
      const f = FB.fieldAt(sol, x, y, z);
      if (!(f.d > 0.08*R && f.d < 0.7*R)) continue;
      pts.push({ x, y, z, f, inside: Math.abs(solidAngle(sol, { x, y, z })) > Math.PI });
    }

    section('  ' + preset + ', ' + pts.length + ' random points (' +
            pts.filter(p => p.inside).length + ' of them inside the loop)');
    for (const mode of ['E', 'B', 'S']) {
      const streams = FB.buildStreams(sol, { mode, scene: R });
      const din = [], dout = [], ang = [];
      for (const p of pts) {
        const q = FB.modeVec(p.f, mode);
        if (!(Math.hypot(q.x, q.y, q.z) > 0)) continue;
        let best = Infinity, bl = null, bi = -1;
        for (const line of streams) for (let i = 0; i < line.length; i++) {
          const dx = line[i][0]-p.x, dy = line[i][1]-p.y, dz = line[i][2]-p.z;
          const d2 = dx*dx + dy*dy + dz*dz;
          if (d2 < best) { best = d2; bl = line; bi = i; }
        }
        const d = Math.sqrt(best)/R;
        (p.inside ? din : dout).push(d);
        // Where a line does pass close, it had better be going the way the
        // field goes -- and the same way, not merely along the same axis: the
        // backward half is reversed before it is joined on, so a whole line
        // runs with +F from end to end.
        if (d < 0.06 && bi < bl.length-1) {
          const ux = bl[bi+1][0]-bl[bi][0], uy = bl[bi+1][1]-bl[bi][1], uz = bl[bi+1][2]-bl[bi][2];
          const g = FB.modeVec(FB.fieldAt(sol, bl[bi][0], bl[bi][1], bl[bi][2]), mode);
          const un = Math.hypot(ux, uy, uz), gn = Math.hypot(g.x, g.y, g.z);
          if (un > 0 && gn > 0)
            ang.push(Math.acos(Math.min(1, Math.max(-1, (ux*g.x + uy*g.y + uz*g.z)/(un*gn)))) * 180/Math.PI);
        }
      }
      const mIn = median(din), mOut = median(dout), mAng = median(ang);
      check(mode + ': the inside of the loop is covered like the outside',
            mIn < 1.35*mOut,
            'nearest line: ' + mIn.toFixed(3) + 'R inside vs ' + mOut.toFixed(3) + 'R outside' +
            ' (ratio ' + (mIn/mOut).toFixed(2) + 'x)');
      // A backstop under the ratio check above, which is the one that catches a
      // whole region going missing. The bound is loose because this bench thins
      // on purpose -- dup is 0.12R here against 0.075R on the switch-on bench,
      // and GRID is 4 against 5 -- so it draws about half as many lines and the
      // gaps between them are correspondingly wider. Worst measured is 0.47R.
      check(mode + ': and no random point is stranded far from every line',
            din.concat(dout).filter(d => d > 0.55).length === 0,
            'worst ' + Math.max(...din, ...dout).toFixed(3) + 'R of ' + (din.length+dout.length) + ' points');
      check(mode + ': a line that passes close is going the way the field goes',
            mAng < 3 && Math.max(...ang) < 25,
            'tangent vs field: median ' + mAng.toFixed(1) + ' deg, worst ' +
            Math.max(...ang).toFixed(1) + ' deg over ' + ang.length + ' points');
    }
  }
}

console.log('\n' + (failures ? failures + ' CHECK(S) FAILED' : 'all checks passed') + '\n');
process.exit(failures ? 1 : 0);
