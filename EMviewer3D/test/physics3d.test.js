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

console.log('\n' + (failures ? failures + ' CHECK(S) FAILED' : 'all checks passed') + '\n');
process.exit(failures ? 1 : 0);
