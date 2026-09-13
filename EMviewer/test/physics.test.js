/*
 * Physics regression tests. No dependencies.
 *   node test/physics.test.js
 *
 * These check the solver against things we know independently: Ohm's law,
 * the infinite-wire field, and the direction of energy flow.
 */
const FB = require('../src/solver.js');
const { MU0 } = FB;

let failures = 0;
function check(name, cond, extra) {
  if (!cond) failures++;
  console.log((cond ? '  ok   ' : '  FAIL ') + name + (extra ? '   [' + extra + ']' : ''));
}
function section(t) { console.log('\n' + t); }

const base = {
  ctrl: FB.PRESETS.rect.map(p => ({ x: p[0], y: p[1] })),
  N: 260, per: 60, radius: 0.003, emf: 3, rLED: 100, rho: 0
};
const build = over => FB.solveCircuit(Object.assign({}, base, over));

// probe just outside the conductor, offset along the outward normal
function probe(sol, frac, off) {
  const i = Math.floor(sol.N * frac);
  const nx = -sol.ty[i], ny = sol.tx[i];
  const f = FB.fieldAt(sol, sol.cx[i] + nx * off, sol.cy[i] + ny * off);
  return { i, nx, ny, f, Et: f.ex * sol.tx[i] + f.ey * sol.ty[i], En: f.ex * nx + f.ey * ny };
}
const firstOf = (sol, t) => { for (let i = 0; i < sol.N; i++) if (sol.type[i] === t) return i; return -1; };

/* ------------------------------------------------------------------ */
section('Bookkeeping');
let sol = build({});
check('solve succeeds', !!sol);
check('loop closes in potential', Math.abs(sol.closure) < 1e-12, 'dV=' + sol.closure.toExponential(2));
check('current is emf/R', Math.abs(sol.I - 3 / 100) < 1e-12, 'I=' + sol.I);
let Q = 0; for (let i = 0; i < sol.N; i++) Q += sol.lam[i] * sol.len[i];
check('circuit is charge neutral', Math.abs(Q) < 1e-24, 'Q=' + Q.toExponential(2));

// the boundary condition is what defines the answer, so verify it directly
let maxres = 0;
for (let j = 0; j < sol.N; j++) {
  let v = 0;
  for (let i = 0; i < sol.N; i++) {
    if (i === j) v += sol.lam[i] * FB.KC * 2 * Math.asinh(sol.len[i] / (2 * sol.a));
    else {
      const dx = sol.cx[j] - sol.cx[i], dy = sol.cy[j] - sol.cy[i];
      v += sol.lam[i] * FB.KC * sol.len[i] / Math.sqrt(dx * dx + dy * dy + sol.a * sol.a);
    }
  }
  let Vt = 0;
  for (let i = 0; i < j; i++) Vt += -sol.Et[i] * sol.len[i];
  Vt += -sol.Et[j] * sol.len[j] / 2;
  maxres = Math.max(maxres, Math.abs(v - sol.V0 - Vt));
}
check('boundary condition satisfied', maxres < 1e-10, 'max residual ' + maxres.toExponential(2) + ' V');

/* ------------------------------------------------------------------ */
section('Ideal wire: the field inside is zero, the field outside is perpendicular');
const iw = Math.floor(sol.N * 0.25);
const inside = FB.fieldAt(sol, sol.cx[iw], sol.cy[iw]);
check('E = 0 inside ideal wire', Math.hypot(inside.ex, inside.ey) < 1e-9,
  '|E|=' + Math.hypot(inside.ex, inside.ey).toExponential(2) + ' V/m');
let p = probe(sol, 0.25, sol.a * 1.8);
check('E meets ideal wire at right angles', Math.abs(p.Et / p.En) < 0.05,
  'Et/En=' + Math.abs(p.Et / p.En).toFixed(4));

/* ------------------------------------------------------------------ */
section('Ohm\'s law holds inside the components');
const li = firstOf(sol, FB.LED), bi = firstOf(sol, FB.BATTERY);
const inLED = FB.fieldAt(sol, sol.cx[li], sol.cy[li]);
const expLED = sol.I * base.rLED / sol.ledLen;
check('E inside LED = IR/length', Math.abs(Math.hypot(inLED.ex, inLED.ey) - expLED) / expLED < 0.02,
  Math.hypot(inLED.ex, inLED.ey).toFixed(1) + ' vs ' + expLED.toFixed(1) + ' V/m');
check('E inside LED drives the current', inLED.ex * sol.tx[li] + inLED.ey * sol.ty[li] > 0);
const inBat = FB.fieldAt(sol, sol.cx[bi], sol.cy[bi]);
check('E inside battery opposes the current', inBat.ex * sol.tx[bi] + inBat.ey * sol.ty[bi] < 0,
  'Et=' + (inBat.ex * sol.tx[bi] + inBat.ey * sol.ty[bi]).toFixed(1) + ' V/m');

/* ------------------------------------------------------------------ */
section('Magnetic field');
const dOff = sol.a * 3;
const pb = probe(sol, 0.25, dOff);
const Bexp = MU0 * sol.I / (2 * Math.PI * dOff);
check('B matches mu0 I / 2 pi d near the wire', Math.abs(Math.abs(pb.f.bz) - Bexp) / Bexp < 0.15,
  (Math.abs(pb.f.bz) * 1e6).toFixed(3) + ' vs ' + (Bexp * 1e6).toFixed(3) + ' uT');
check('B vanishes on the wire axis', Math.abs(inside.bz) < Bexp * 0.05);
// far-field point approximation must agree with the exact near-field formula
{
  const i = Math.floor(sol.N * 0.25);
  const far = FB.fieldAt(sol, sol.cx[i] - sol.ty[i] * 0.05, sol.cy[i] + sol.tx[i] * 0.05);
  check('B sign is consistent between near and far branches',
    Math.sign(far.bz) === Math.sign(pb.f.bz), 'near=' + pb.f.bz.toExponential(2) + ' far=' + far.bz.toExponential(2));
}

/* ------------------------------------------------------------------ */
section('Poynting flux goes where the power goes');
{
  const i = pb.i;
  const S = FB.poynting(pb.f);
  const along = S[0] * sol.tx[i] + S[1] * sol.ty[i];
  const into = -(S[0] * pb.nx + S[1] * pb.ny);
  check('S runs parallel to ideal wire, not into it', Math.abs(into / along) < 0.06,
    'along=' + along.toFixed(2) + ' into=' + into.toFixed(3) + ' W/m^2');
}
{
  const nx = -sol.ty[li], ny = sol.tx[li];
  const S = FB.poynting(FB.fieldAt(sol, sol.cx[li] + nx * sol.a * 1.8, sol.cy[li] + ny * sol.a * 1.8));
  check('S flows into the LED', -(S[0] * nx + S[1] * ny) > 0,
    'inward ' + (-(S[0] * nx + S[1] * ny)).toFixed(1) + ' W/m^2');
}
{
  const nx = -sol.ty[bi], ny = sol.tx[bi];
  const S = FB.poynting(FB.fieldAt(sol, sol.cx[bi] + nx * sol.a * 1.8, sol.cy[bi] + ny * sol.a * 1.8));
  check('S flows out of the battery', S[0] * nx + S[1] * ny > 0,
    'outward ' + (S[0] * nx + S[1] * ny).toFixed(1) + ' W/m^2');
}

/* ------------------------------------------------------------------ */
section('Resistive wire: the interior field switches on');
const rho = 60;
sol = build({ rho });
const ri = Math.floor(sol.N * 0.25);
{
  const f = FB.fieldAt(sol, sol.cx[ri], sol.cy[ri]);
  const exp = sol.I * rho;
  check('E inside resistive wire = rho_linear * I', Math.abs(Math.hypot(f.ex, f.ey) - exp) / exp < 0.03,
    Math.hypot(f.ex, f.ey).toFixed(2) + ' vs ' + exp.toFixed(2) + ' V/m');
  const q = probe(sol, 0.25, sol.a * 1.8);
  check('E outside now tilts along the wire', Math.abs(q.Et) > 0.4 * sol.I * rho,
    'Et=' + q.Et.toFixed(2) + ' V/m, interior ' + exp.toFixed(2));
  const nx = -sol.ty[ri], ny = sol.tx[ri];
  const S = FB.poynting(FB.fieldAt(sol, sol.cx[ri] + nx * sol.a * 1.8, sol.cy[ri] + ny * sol.a * 1.8));
  check('S now also flows into the plain wire', -(S[0] * nx + S[1] * ny) > 0,
    'inward ' + (-(S[0] * nx + S[1] * ny)).toFixed(3) + ' W/m^2');
}
{
  const zero = build({ rho: 0 });
  check('zero-resistance limit recovers zero interior field',
    Math.hypot(FB.fieldAt(zero, zero.cx[ri], zero.cy[ri]).ex, 0) < 1e-9);
}

/* ------------------------------------------------------------------ */
section('Convergence in segment count');
{
  const coarse = build({ N: 140, per: 34 });
  const fine = build({ N: 380, per: 80 });
  const pc = probe(coarse, 0.25, 0.02).En, pf = probe(fine, 0.25, 0.02).En;
  check('field away from the wire is resolution independent', Math.abs(pc - pf) / Math.abs(pf) < 0.03,
    'N=140 ' + pc.toFixed(2) + ' vs N=380 ' + pf.toFixed(2) + ' V/m');
}

/* ------------------------------------------------------------------ */
section('All shape presets solve');
for (const k of Object.keys(FB.PRESETS)) {
  const s = build({ ctrl: FB.PRESETS[k].map(p => ({ x: p[0], y: p[1] })) });
  check(k.padEnd(7) + ' L=' + (s.L * 100).toFixed(1) + 'cm  I=' + (s.I * 1000).toFixed(2) + 'mA',
    !!s && isFinite(s.peak) && s.peak > 0);
}

/* ------------------------------------------------------------------ */
section('Timing (render budget is 16 ms/frame while dragging)');
for (const q of [{ N: 140, per: 34, cell: 10, label: 'drag ' }, { N: 260, per: 60, cell: 5, label: 'settled' }]) {
  const s0 = build({ N: q.N, per: q.per });
  let t = Date.now();
  for (let i = 0; i < 5; i++) build({ N: q.N, per: q.per });
  const tSolve = (Date.now() - t) / 5;
  t = Date.now();
  FB.buildGrid(s0, { w: 1000, h: 640, cell: q.cell, ppm: 1000, needE: true, needB: false });
  const tE = Date.now() - t;
  t = Date.now();
  FB.buildGrid(s0, { w: 1000, h: 640, cell: q.cell, ppm: 1000, needE: true, needB: true });
  const tEB = Date.now() - t;
  console.log(`  ${q.label}  solve ${tSolve.toFixed(0)}ms   grid(E) ${tE}ms   grid(E+B) ${tEB}ms`);
}

console.log('\n' + (failures ? failures + ' FAILURE(S)' : 'all checks passed'));
process.exit(failures ? 1 : 0);
