/*
 * Physics regression tests for the switch-on bench. No dependencies.
 *   node test/physics4d.test.js
 *
 * The things worth doubting here are the inductance, the circuit march, whether
 * the source conserves charge, and whether the field really is retarded. Each
 * is checked against something known independently of this code.
 */
const F = require('../src/solver4d.js');
const { MU0, C } = F;

let failures = 0;
function check(name, cond, extra) {
  if (!cond) failures++;
  console.log((cond ? '  ok   ' : '  FAIL ') + name + (extra ? '   [' + extra + ']' : ''));
}
function section(t) { console.log('\n' + t); }
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);
const ctrlOf = p => p.map(q => ({ x: q[0], y: q[1], z: q[2] }));

const base = {
  N: 200, per: 44, radius: 0.003, emf: 9, rLED: 100, rho: 0,
  tau: 2e-9, tMax: 60e-9, steps: 600
};
const build = over => F.buildTransient(Object.assign({ ctrl: ctrlOf(F.PRESETS.rect) }, base, over));

/* ------------------------------------------------------------------ */
section('Self-inductance, against the thin-ring formula');
// A circular loop of radius R and wire radius a has L = mu0 R [ln(8R/a) - 2].
// Nothing about that formula appears in the solver: L there is the loop integral
// of the same geometric kernel the charge solve uses, weighted by t.t'.
for (const [R, a] of [[0.2, 0.003], [0.2, 0.001], [0.3, 0.004]]) {
  const ring = [];
  for (let i = 0; i < 14; i++) { const t = i/14*2*Math.PI; ring.push([R*Math.cos(t), R*Math.sin(t), 0]); }
  const T = build({ ctrl: ctrlOf(ring), radius: a, N: 400 });
  const want = MU0 * R * (Math.log(8*R/a) - 2);
  check('ring R=' + R + ' m, a=' + (a*1000) + ' mm', rel(T.Lind, want) < 0.06,
        'got ' + (T.Lind*1e6).toFixed(3) + ' uH, formula ' + (want*1e6).toFixed(3) +
        ' uH (' + (100*rel(T.Lind, want)).toFixed(1) + '%)');
}

/* ------------------------------------------------------------------ */
section('The circuit march');
{
  const T = build({});
  // Residual of L I' + R I = e(t), with I' taken from the marched history by
  // central difference -- so this tests the march, not the formula it came from.
  let worst = 0, at = 0;
  for (let k = 2; k < T.steps - 2; k++) {
    const t = k * T.dt;
    // five-point stencil: a central difference carries O(dt^2) error, which on
    // a 2 ns ramp is bigger than anything the march itself gets wrong
    const Id = (-T.Ihist[k+2] + 8*T.Ihist[k+1] - 8*T.Ihist[k-1] + T.Ihist[k-2]) / (12 * T.dt);
    const e = F.emfAt(T.emf, T.tau, t).v;
    const res = Math.abs(T.Lind * Id + T.Rtot * T.Ihist[k] - e) / T.emf;
    if (res > worst) { worst = res; at = t; }
  }
  check('L I\' + R I = e(t) along the whole march', worst < 2e-3,
        'worst residual ' + worst.toExponential(2) + ' of EMF, at ' + (at*1e9).toFixed(1) + ' ns');

  check('current starts from zero', T.Ihist[0] === 0);
  // the march ends at 60 ns, which is only about four L/R -- so compare against
  // where the analytic solution is at 60 ns, not against the asymptote
  const wantEnd = T.Iss * (1 - Math.exp(-(T.tMax - T.tau/2) / T.tauLR));
  check('current tracks the analytic rise', rel(T.Ihist[T.steps], wantEnd) < 0.01,
        (T.Ihist[T.steps]*1000).toFixed(2) + ' mA vs ' + (wantEnd*1000).toFixed(2) +
        ' mA at ' + (T.tMax*1e9) + ' ns (asymptote ' + (T.Iss*1000).toFixed(1) + ')');
  // and the march must be converged: halving the step must not move the answer
  const fine = build({ steps: base.steps * 4 });
  check('march is converged in step size', rel(T.Ihist[T.steps], fine.Ihist[fine.steps]) < 1e-4,
        'steps=' + base.steps + ' vs ' + (base.steps*4) + ': ' +
        rel(T.Ihist[T.steps], fine.Ihist[fine.steps]).toExponential(2));
  // after a step, I should rise as 1 - exp(-t/(L/R)); with a 2 ns ramp the
  // late-time approach is still exponential with that time constant
  const k1 = Math.round(20e-9 / T.dt), k2 = Math.round(40e-9 / T.dt);
  const decay = (T.Iss - T.Ihist[k2]) / (T.Iss - T.Ihist[k1]);
  const want = Math.exp(-(20e-9) / T.tauLR);
  check('approach is exponential with time constant L/R', rel(decay, want) < 0.05,
        'ratio over 20 ns ' + decay.toFixed(4) + ' vs exp(-20/' + (T.tauLR*1e9).toFixed(2) + ') = ' + want.toFixed(4));
}

/* ------------------------------------------------------------------ */
section('The source conserves charge');
{
  const T = build({});
  // continuity on a filament:  d(lambda)/dt + dI/ds = 0
  let worst = 0;
  for (const t of [0.5e-9, 1.5e-9, 3e-9, 8e-9]) {
    let scale = 0;
    for (let i = 0; i < T.N; i++) scale = Math.max(scale, Math.abs(F.sourceAt(T, i, t).lamd));
    for (let i = 0; i < T.N; i++) {
      const j = (i + 1) % T.N;
      const si = F.sourceAt(T, i, t), sj = F.sourceAt(T, j, t);
      const dIds = (sj.cur - si.cur) / (0.5 * (T.len[i] + T.len[j]));
      const lamd = 0.5 * (si.lamd + sj.lamd);
      worst = Math.max(worst, Math.abs(dIds + lamd) / scale);
    }
  }
  check('dI/ds = -d(lambda)/dt around the whole loop', worst < 0.02,
        'worst ' + (100*worst).toFixed(2) + '% of peak d(lambda)/dt');

  // and the analytic time derivative must match a numerical one
  const h = 1e-12; let wd = 0, sc = 0;
  for (let i = 0; i < T.N; i += 7) {
    const a = F.sourceAt(T, i, 3e-9 - h), b = F.sourceAt(T, i, 3e-9 + h);
    const num = (b.lam - a.lam) / (2*h);
    const ana = F.sourceAt(T, i, 3e-9).lamd;
    wd = Math.max(wd, Math.abs(num - ana)); sc = Math.max(sc, Math.abs(ana));
  }
  check('analytic d(lambda)/dt matches a numerical one', wd/sc < 0.02,
        'worst ' + (100*wd/sc).toFixed(2) + '%');

  // the loop stays neutral at every instant
  let qw = 0;
  for (const t of [0.5e-9, 2e-9, 10e-9, 50e-9]) {
    let q = 0;
    for (let i = 0; i < T.N; i++) q += F.sourceAt(T, i, t).lam * T.len[i];
    qw = Math.max(qw, Math.abs(q));
  }
  check('loop is neutral at every instant', qw < 1e-22, 'worst |Q| = ' + qw.toExponential(2) + ' C');
}

/* ------------------------------------------------------------------ */
section('The field is retarded');
{
  const T = build({});
  const p = [1.4, 0.9, 0.7];
  let dmin = Infinity;
  for (let i = 0; i < T.N; i++)
    dmin = Math.min(dmin, Math.hypot(p[0]-T.cx[i], p[1]-T.cy[i], p[2]-T.cz[i]));
  const arrive = dmin / C;
  const mag = t => { const f = F.fieldAt(T, p[0], p[1], p[2], t);
                     return Math.hypot(f.ex, f.ey, f.ez); };
  check('field is exactly zero before the news could arrive',
        mag(arrive * 0.98) === 0, 'probe at ' + dmin.toFixed(3) + ' m, arrival ' + (arrive*1e9).toFixed(3) + ' ns');
  check('field is non-zero shortly after', mag(arrive * 1.05) > 0);
  // the arrival time must track the distance, not just be early
  const q = [2.8, 1.8, 1.4];
  let dq = Infinity;
  for (let i = 0; i < T.N; i++)
    dq = Math.min(dq, Math.hypot(q[0]-T.cx[i], q[1]-T.cy[i], q[2]-T.cz[i]));
  const magq = t => { const f = F.fieldAt(T, q[0], q[1], q[2], t);
                      return Math.hypot(f.ex, f.ey, f.ez); };
  check('a probe twice as far waits twice as long',
        magq(dq/C * 0.98) === 0 && magq(dq/C * 1.05) > 0 && dq > 1.9 * dmin,
        'near ' + dmin.toFixed(2) + ' m, far ' + dq.toFixed(2) + ' m');
}

/* ------------------------------------------------------------------ */
section('Radiation: the far field falls off as 1/R');
{
  // Static fields go as 1/R^2, the radiated part as 1/R. Sampling each radius at
  // its own retarded moment compares the same instant of the source.
  const T = build({ tau: 1e-9 });
  const tPeak = 0.5e-9;
  // In the PLANE of the loop. Magnetic dipole radiation goes as sin(theta) from
  // the dipole axis, so the loop axis is its null -- measuring the falloff there
  // reads the near field and gives a slope between -1 and -2.
  const dir = [1, 0, 0];
  const at = R => {
    const f = F.fieldAt(T, dir[0]*R, dir[1]*R, dir[2]*R, tPeak + R/C);
    return Math.hypot(f.bx, f.by, f.bz);
  };
  const R1 = 200, R2 = 400;
  const slope = Math.log(at(R2)/at(R1)) / Math.log(R2/R1);
  check('|B| ~ 1/R in the far field during the switch', Math.abs(slope + 1) < 0.1,
        'log-log slope ' + slope.toFixed(3) + ' (radiation -1, static -2)');
  // and once settled there is no radiation left, so it goes back to 1/R^2... R^3
  // on the axis of a loop, which is the dipole result
  const late = R => {
    const f = F.fieldAt(T, 0, 0, R, 300e-9);
    return Math.hypot(f.bx, f.by, f.bz);
  };
  const ls = Math.log(late(8)/late(4)) / Math.log(2);
  check('settled, the axial field is a dipole falling as 1/R^3', Math.abs(ls + 3) < 0.06,
        'log-log slope ' + ls.toFixed(3));
}

/* ------------------------------------------------------------------ */
section('Late time reproduces the steady state');
{
  const T = build({});
  const tLate = 900e-9;                       // many L/R past the end of the march
  // the charge density must match the DC boundary-value solution
  let worst = 0, peak = 0;
  for (let i = 0; i < T.N; i++) {
    const lam = F.sourceAt(T, i, tLate).lam;
    const dc  = T.emf * T.lamE[i] + T.Iss * T.lamI[i];
    worst = Math.max(worst, Math.abs(lam - dc)); peak = Math.max(peak, Math.abs(dc));
  }
  check('charge density settles to the DC solution', worst/peak < 1e-5,
        'worst ' + (100*worst/peak).toExponential(2) + '%');
  let cw = 0, cp = 0;
  for (let i = 0; i < T.N; i++) {
    const cur = F.sourceAt(T, i, tLate).cur;
    cw = Math.max(cw, Math.abs(cur - T.Iss)); cp = T.Iss;
  }
  check('current settles to emf/R everywhere along the wire', cw/cp < 1e-5,
        'worst ' + (100*cw/cp).toExponential(2) + '%');
}

/* ------------------------------------------------------------------ */
// The 3D steady-state bench is a sibling project, not a dependency. When it is
// present, the settled field here must be the field it computes.
section('Late time reproduces the steady-state bench');
let FB3 = null;
try { FB3 = require('../../EMviewer3D/src/solver3d.js'); } catch (e) { /* standalone */ }
if (!FB3) console.log('  skip  3D field bench not present alongside');
else {
  const probes = [[0.1,0.3,0.05],[0.35,0.05,0],[0,0,0.2],[0.27,0.16,-0.1]];
  const compare = N => {
    const o = { N, per: 44, radius: 0.003, emf: 9, rLED: 100, rho: 0 };
    const T = build(o);
    const S = FB3.solveCircuit(Object.assign({ ctrl: ctrlOf(FB3.PRESETS.rect) }, o));
    let we = 0, wb = 0;
    for (const p of probes) {
      const R = Math.hypot(p[0],p[1],p[2]);
      const f4 = F.fieldAt(T, p[0], p[1], p[2], 900e-9 + R/F.C);
      const f3 = FB3.fieldAt(S, p[0], p[1], p[2]);
      we = Math.max(we, rel(Math.hypot(f4.ex,f4.ey,f4.ez), Math.hypot(f3.ex,f3.ey,f3.ez)));
      wb = Math.max(wb, rel(Math.hypot(f4.bx,f4.by,f4.bz), Math.hypot(f3.bx,f3.by,f3.bz)));
    }
    return { we, wb };
  };
  const c200 = compare(200);
  check('settled E matches the steady-state bench', c200.we < 1e-6,
        'worst ' + c200.we.toExponential(2) + ' relative');
  check('settled B matches the steady-state bench', c200.wb < 1e-6,
        'worst ' + c200.wb.toExponential(2) + ' relative');

  // Close to the wire the two part company by construction: the steady-state
  // bench integrates each segment exactly in space, while this one is a point
  // element at the segment centre, because that is where its retarded time is
  // taken too. Refining the space integral alone would be inconsistent, so the
  // difference stays -- but it has to be discretisation, which means this
  // solver's own answer must converge as the wire is cut finer. Measured
  // against itself, since the bench switches to point elements inside four
  // segment lengths and stops being an independent reference there.
  const nearB = N => {
    const T = build({ N, per: 44, radius: 0.003, emf: 9, rLED: 100, rho: 0 });
    const p = [0.272, 0, 0];                 // a fixed 12 mm outside the right-hand run
    const R = Math.hypot(p[0], p[1], p[2]);
    const f = F.fieldAt(T, p[0], p[1], p[2], 900e-9 + R/F.C);
    return Math.hypot(f.bx, f.by, f.bz);
  };
  const b1 = nearB(200), b2 = nearB(400), b3 = nearB(800);
  const d1 = Math.abs(b1 - b2), d2 = Math.abs(b2 - b3);
  check('B near the wire converges as the wire is cut finer', d2 < d1,
        '200->400 moves it ' + (100*d1/b3).toFixed(3) + '%, 400->800 ' + (100*d2/b3).toFixed(3) + '%');
  check('and it converges quadratically', d1 / d2 > 3.2,
        'successive differences shrink ' + (d1/d2).toFixed(1) + 'x (quadratic is 4x)');
}

console.log('\n' + (failures ? failures + ' CHECK(S) FAILED' : 'all checks passed') + '\n');
process.exit(failures ? 1 : 0);
