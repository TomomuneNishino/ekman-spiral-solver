// Browser front end.
//
// The physics lives entirely in solver.js / analytic.js / diagnostics.js /
// upwelling.js / checks.js. This file only reads the controls, calls those
// modules and draws the result.

import { Integrator, solveSteady, cellCentres, ekmanDepth, coriolis, OMEGA_EARTH } from './solver.js';
import * as A from './analytic.js';
import * as D from './diagnostics.js';
import { solveUpwelling } from './upwelling.js';
import { runChecks } from './checks.js';
import { Frame, span, depthColor, palette } from './plot.js';

const $ = (id) => document.getElementById(id);
const REDUCED = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

const DEFAULTS = { lat: 45, taum: 0.1, taud: 0, nu: -2, hd: 6, nz: 120, steps: 120, bot: 'free', sch: 'cn' };
// The spin-up starts paused at t = 0, so the button label and the animation
// always agree and the layer is seen growing from rest on the first press.
const state = { ...DEFAULTS, tab: 'spiral', frame: 0, playing: false, reveal: REDUCED ? 1 : 0 };

let P = null;        // physical parameters
let steady = null;   // steady solution and its diagnostics
let run = null;      // spin-up history
let pump = null;     // Ekman pumping solution

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------
function readControls() {
  for (const k of Object.keys(DEFAULTS)) {
    const el = $(k);
    state[k] = el.type === 'range' ? parseFloat(el.value) : el.value;
  }
  const f = coriolis(state.lat) || coriolis(1);
  const nu = 10 ** state.nu;
  const dir = (state.taud * Math.PI) / 180;
  const p = {
    f, nu, rho0: 1025,
    taux: state.taum * Math.cos(dir),
    tauy: state.taum * Math.sin(dir),
    bottom: state.bot,
  };
  p.DE = ekmanDepth(p);
  p.H = state.hd * p.DE;
  p.N = Math.round(state.nz);
  p.inertial = (2 * Math.PI) / Math.abs(p.f);
  p.dt = p.inertial / state.steps;
  return p;
}

const COMPASS = ['東', '東北東', '北東', '北北東', '北', '北北西', '北西', '西北西',
  '西', '西南西', '南西', '南南西', '南', '南南東', '南東', '東南東'];

function syncLabels() {
  const lat = state.lat;
  $('lat-v').textContent = lat === 0 ? '赤道 (f → 0)' : `${lat > 0 ? '北緯' : '南緯'} ${Math.abs(lat)}°`;
  $('taum-v').textContent = `${state.taum.toFixed(2)} N/m²`;
  $('taud-v').textContent = `${COMPASS[Math.round(state.taud / 22.5) % 16]}向き`;
  $('nu-v').textContent = `${(10 ** state.nu).toExponential(1)} m²/s`;
  // Depth in Ekman depths is the meaningful number, but metres are what a
  // reader can picture, so show both.
  $('hd-v').textContent = P ? `${state.hd} D_E ＝ ${(state.hd * P.DE).toFixed(0)} m`
    : `${state.hd} D_E`;
  $('nz-v').textContent = `${state.nz}`;
  $('steps-v').textContent = `${state.steps}`;
}

// ---------------------------------------------------------------------------
// Computation
// ---------------------------------------------------------------------------
function recompute() {
  P = readControls();
  syncLabels();

  const W = solveSteady(P);
  const z = Array.from(cellCentres(P));

  // The analytic solution is needed at two different resolutions. Errors must be
  // measured where the solver actually holds values, i.e. at cell centres. But a
  // reference curve sampled at cell centres turns into a visible polygon on a
  // coarse grid, which reads as though the closed-form spiral itself were
  // piecewise linear. So the drawn curve gets its own dense sampling, spanning
  // the full domain from the surface down to the floor.
  const nFine = 400;
  const zFine = [];
  for (let k = 0; k <= nFine; k++) zFine.push(-(k / nFine) * P.H);

  steady = {
    W, z, zFine,
    refFine: A.steadyProfile(P, zFine),
    ref: A.steadyProfile(P, z),
    d: D.diagnose(W, P),
    Mex: A.transportClassical(P),
    surf: A.surfaceVelocity(P),
    Pex: A.steadyPower(P),
  };

  // Spin-up over five inertial periods, storing scalars every step and velocity
  // profiles at a subsample suitable for animation.
  const nsteps = Math.round(5 * state.steps);
  const keep = Math.max(1, Math.round(nsteps / 220));
  const it = new Integrator(P, { scheme: state.sch, dt: P.dt });
  const h = {
    t: [0], Mx: [0], My: [0], Mxe: [0], Mye: [0], E: [0],
    Pw: [0], eps: [0], cumP: [0], cumE: [0], resM: [], resE: [], tRes: [],
    frames: [{ re: it.W.re.slice(), im: it.W.im.slice(), t: 0 }],
    blewUp: false,
  };
  let cumP = 0, cumE = 0;
  for (let n = 0; n < nsteps; n++) {
    const { Wold, Wmid, Wnew } = it.step();
    const Pw = D.windPowerPhysical(Wmid, P);
    const eps = D.dissipationPhysical(Wmid, P);
    cumP += Pw * it.dt; cumE += eps * it.dt;
    const M = D.transport(it.W, P);
    const Me = A.spinupTransport(P, it.t);
    h.t.push(it.t);
    h.Mx.push(M.re); h.My.push(M.im);
    h.Mxe.push(Me.re); h.Mye.push(Me.im);
    h.E.push(D.energy(it.W, P));
    h.Pw.push(Pw); h.eps.push(eps);
    h.cumP.push(cumP); h.cumE.push(cumE);
    h.tRes.push(it.t);
    h.resM.push(D.momentumResidual(Wold, Wnew, Wmid, P, it.dt).relative);
    h.resE.push(Math.abs(D.energyResidual(Wold, Wnew, Wmid, P, it.dt).relative));
    if ((n + 1) % keep === 0) h.frames.push({ re: it.W.re.slice(), im: it.W.im.slice(), t: it.t });
    if (!Number.isFinite(M.re) || Math.abs(M.re) > 1e12) { h.blewUp = true; break; }
  }
  // Same reasoning for the exact spin-up trajectory: at 10 steps per inertial
  // period the per-step samples would draw a decagon rather than a circle.
  const tEnd = h.t[h.t.length - 1] || 1;
  h.tFine = []; h.MxeFine = []; h.MyeFine = [];
  for (let k = 0; k <= 600; k++) {
    const tt = (k / 600) * tEnd;
    const Me = A.spinupTransport(P, tt);
    h.tFine.push(tt); h.MxeFine.push(Me.re); h.MyeFine.push(Me.im);
  }

  run = h;
  state.frame = Math.min(state.frame, run.frames.length - 1);
  $('tslider').max = String(run.frames.length - 1);

  pump = solveUpwelling(
    { ...P, taux: 0, tauy: 0 },
    { Ly: 4.0e6, Ny: 96, tau0: Math.hypot(P.taux, P.tauy) },
  );

  renderReadout();
  render();
}

// ---------------------------------------------------------------------------
// Readout strip
// ---------------------------------------------------------------------------
/** Readable number that stays readable when a run diverges. */
const num = (v, d = 1) =>
  (!Number.isFinite(v) ? '∞'
    : Math.abs(v) >= 1e5 || (v !== 0 && Math.abs(v) < 1e-2) ? v.toExponential(2)
    : v.toFixed(d));

const hhmm = (s) => (s < 7200 ? `${(s / 60).toFixed(0)} 分`
  : s < 2 * 86400 ? `${(s / 3600).toFixed(1)} 時間` : `${(s / 86400).toFixed(1)} 日`);

function renderReadout() {
  const { d, Mex } = steady;
  const speed = Math.hypot(d.surface.re, d.surface.im);
  const defl = A.deflectionDeg(d.surface, P);
  const Mmag = Math.hypot(d.M.re, d.M.im);
  const Mdir = (Math.atan2(d.M.im, d.M.re) * 180) / Math.PI;
  const rows = [
    ['エクマン深さ D<sub>E</sub>', `${P.DE.toFixed(1)}<small> m</small>`],
    ['海面流速', `${(speed * 100).toFixed(1)}<small> cm/s</small>`],
    ['風からの偏角', `${defl.toFixed(1)}<small>°（${defl < 0 ? '右' : '左'}へ）</small>`],
    ['鉛直積分輸送', `${Mmag.toFixed(0)}<small> kg/(m·s)</small>`],
    ['輸送の向き', `${Mdir.toFixed(0)}<small>°（東を 0°）</small>`],
    ['風の仕事率', `${(d.P * 1e3).toFixed(2)}<small> mW/m²</small>`],
    ['慣性周期', `${hhmm(P.inertial)}`],
    ['セル幅', `${(P.H / P.N / P.DE).toFixed(3)}<small> D_E</small>`],
  ];
  $('readout').innerHTML = rows
    .map(([k, v]) => `<div><dt>${k}</dt><dd>${v}</dd></div>`).join('');

  const relM = Math.hypot(d.M.re - Mex.re, d.M.im - Mex.im) / Math.hypot(Mex.re, Mex.im);
  const kb = document.getElementById('key-bot');
  if (kb) kb.textContent = `水深 ${P.H.toFixed(0)} m`;

  // The colour key is built from depthColor() rather than hard-coded in CSS, so
  // it can never drift out of step with what the figures actually draw. Ticks
  // are placed at whole Ekman depths, which is the scale the spiral lives on.
  const ramp = $('key-ramp');
  if (ramp) {
    const stops = [];
    for (let k = 0; k <= 24; k++) stops.push(`${depthColor(k / 24)} ${((k / 24) * 100).toFixed(1)}%`);
    ramp.style.background = `linear-gradient(90deg, ${stops.join(',')})`;
  }
  const ticks = $('key-ticks');
  if (ticks) {
    const step = Math.max(1, Math.ceil(state.hd / 6));
    let html = '';
    for (let d = step; d * P.DE < P.H * 0.94; d += step) {
      html += `<span style="left:${((d * P.DE) / P.H) * 100}%">${d} D_E</span>`;
    }
    ticks.innerHTML = html;
  }

  const dtMax = (P.H / P.N) ** 2 / (2 * P.nu);
  let note = `深さ ${P.H.toFixed(0)} m を ${P.N} 分割、Δt = ${P.dt.toFixed(0)} s。` +
    `輸送と τ/f の差は ${relM.toExponential(1)}。`;
  if (state.sch === 'euler') {
    note += ` 前進Euler法の安定条件は Δt &lt; Δz²/2ν<sub>z</sub> = ${dtMax.toFixed(0)} s。`;
  }
  if (run && run.blewUp) {
    note += ` <b style="color:var(--alert)">この設定では計算が発散したため途中で停止した。</b>`;
  }
  $('rail-note').innerHTML = note;
}

// ---------------------------------------------------------------------------
// Figures
// ---------------------------------------------------------------------------
function hodograph(canvas, W, { ghost = null, reveal = 1, showAnalytic = true, title = true } = {}) {
  const p = palette();
  const zs = steady.z;
  const us = [], vs = [], ts = [];
  const nShow = Math.max(2, Math.round(P.N * reveal));
  for (let k = 0; k < nShow; k++) { us.push(W.re[k]); vs.push(W.im[k]); ts.push(-zs[k] / P.H); }

  // The wind arrow is sized from the data, so the axes can then be fitted to
  // everything that will actually be drawn instead of to a symmetric box
  // around the origin, which would leave the spiral in one small corner.
  const gx = ghost ? Array.from(ghost.re) : [];
  const gy = ghost ? Array.from(ghost.im) : [];
  let R = 0;
  for (const arr of [us, gx]) for (const v of arr) R = Math.max(R, Math.abs(v));
  for (const arr of [vs, gy]) for (const v of arr) R = Math.max(R, Math.abs(v));
  R = R || 1e-6;
  const tm = Math.hypot(P.taux, P.tauy) || 1;
  const ax = (P.taux / tm) * R * 0.95, ay = (P.tauy / tm) * R * 0.95;

  const s0 = D.surfaceVelocityExtrapolated(W, P);
  const [xlo, xhi] = span([us, gx, [0, ax, title ? s0.re : 0]], 0.16);
  const [ylo, yhi] = span([vs, gy, [0, ay, title ? s0.im : 0]], 0.16);

  const f = new Frame(canvas, { pad: { l: 60, r: 30, t: 16, b: 42 } });
  f.begin({
    xlim: [xlo, xhi], ylim: [ylo, yhi], equal: true,
    xlabel: '東向き成分 u [m/s]', ylabel: '北向き成分 v [m/s]',
  });

  if (ghost) f.line(gx, gy, { color: p.rule, width: 5, alpha: 0.95, label: '定常解' });
  if (showAnalytic) {
    f.line(steady.refFine.map((c) => c.re), steady.refFine.map((c) => c.im),
      { color: p.muted, width: 2.2, dash: [6, 5], label: '解析解' });
  }

  // A hodograph is a plot of vector tips. Drawing the vectors themselves, from
  // the origin out to the curve at a sequence of depths, is what turns it back
  // into something you can read as "the current at this depth points there".
  // They go underneath the curve so the spiral stays the dominant shape.
  const surf = D.surfaceVelocityExtrapolated(W, P);
  if (title) {
    const dzc = P.H / P.N;
    const maxDep = Math.min(3, state.hd);
    const seen = new Set();
    for (let n = 6; n >= 1; n--) {
      const k = Math.min(nShow - 1, Math.max(0, Math.round((maxDep * (n / 6) * P.DE) / dzc - 0.5)));
      if (seen.has(k)) continue;
      seen.add(k);
      f.arrow(0, 0, W.re[k], W.im[k],
        { color: depthColor(-zs[k] / P.H), width: 2.6, head: 12, alpha: 0.9 });
    }
    // The surface vector is drawn to z = 0 itself, not to the topmost cell
    // centre half a cell below it, so the angle it makes with the wind is the
    // 45 degrees the theory predicts rather than 45 degrees minus half a cell.
    f.arrow(0, 0, surf.re, surf.im, { color: depthColor(0), width: 3.8, head: 16 });
  }

  f.rampLine(us, vs, ts, { width: 4, label: '数値解（色は深さ）' });
  f.points([title ? surf.re : us[0]], [title ? surf.im : vs[0]], { color: depthColor(0), r: 5.5 });
  f.arrow(0, 0, ax, ay, { color: p.accent, width: 3.6, head: 16, label: '風応力 τ' });

  if (title) {
    // The 45 degrees is the headline result, so put it on the figure rather
    // than only in the caption.
    const defl = A.deflectionDeg(surf, P);
    f.angleMark(0, 0, { x: P.taux / tm, y: P.tauy / tm }, { x: surf.re, y: surf.im }, {
      r: 66, color: p.accent, label: `${Math.abs(defl).toFixed(1)}°`,
    });
    f.annotate(surf.re, surf.im, '海面の流れ', { dx: 11, color: p.ink });
    for (const dep of [1, 2, 3]) {
      const k = Math.round((dep * P.DE) / (P.H / P.N) - 0.5);
      if (k > 0 && k < nShow) {
        f.points([W.re[k]], [W.im[k]], { color: depthColor(-zs[k] / P.H), r: 4.4 });
        // The marker carries the depth colour; the text does not, or it would
        // vanish into the background once the ramp goes dark.
        f.annotate(W.re[k], W.im[k], `${dep} D_E`, { dx: 9, color: p.muted });
      }
    }
  }
  // Sample everything that has been drawn, arrow included, so the legend lands
  // in whichever corner is actually free for this wind direction.
  const px = [...us, ...gx, surf.re], py = [...vs, ...gy, surf.im];
  for (let k = 0; k <= 10; k++) { px.push((ax * k) / 10); py.push((ay * k) / 10); }
  f.legend({ corner: f.pickCorner(px, py) });
}

/** Make a sequence of angles continuous by removing 360-degree jumps. */
function unwrapDeg(a) {
  const out = [a[0]];
  for (let k = 1; k < a.length; k++) {
    let v = a[k];
    while (v - out[k - 1] > 180) v -= 360;
    while (v - out[k - 1] < -180) v += 360;
    out.push(v);
  }
  return out;
}

/** Keep about `want` evenly spaced entries so reference curves stay visible under the data. */
function thin(arrays, want = 42) {
  const n = arrays[0].length;
  const step = Math.max(1, Math.ceil(n / want));
  return arrays.map((a) => a.filter((_, k) => k % step === 0 || k === n - 1));
}

function drawSpiral() {
  hodograph($('c-hodo'), steady.W, { reveal: state.reveal });

  const defl = A.deflectionDeg(steady.d.surface, P);
  let worst = 0, scale = 0;
  for (let k = 0; k < P.N; k++) {
    worst = Math.max(worst, Math.hypot(steady.W.re[k] - steady.ref[k].re,
      steady.W.im[k] - steady.ref[k].im));
    scale = Math.max(scale, Math.hypot(steady.ref[k].re, steady.ref[k].im));
  }
  $('cap-hodo').innerHTML =
    `色つきの線が数値解で、色は海面から水深 ${P.H.toFixed(0)} m までの深さを表す。` +
    `矢印は各深さでの流速ベクトルで、原点から曲線上の点へ伸びている。` +
    `海面の流れは風から<b>${Math.abs(defl).toFixed(1)}° ${defl < 0 ? '右' : '左'}</b>へ逸れ、` +
    `そこから 1 エクマン深さ潜るごとにさらに 1 ラジアンずつ回る。` +
    `螺旋全体を足し合わせた輸送は、風の<b>ちょうど 90° ${P.f > 0 ? '右' : '左'}</b>を向く。` +
    `数値解と解析解の差は、どの深さでも海面流速の` +
    `<b>${((worst / scale) * 100).toPrecision(2)}%</b> を超えない。`;

  // u and v against depth
  {
    const f = new Frame($('c-prof'));
    const zz = steady.z.map((z) => z / P.DE);
    const us = Array.from(steady.W.re), vs = Array.from(steady.W.im);
    const zzF = steady.zFine.map((z) => z / P.DE);
    const [lo, hi] = span([us, vs, steady.refFine.map((c) => c.re), steady.refFine.map((c) => c.im)]);
    f.begin({ xlim: [lo, hi], ylim: [-state.hd, 0], xlabel: '流速 [m/s]', ylabel: '深さ z / D_E' });
    f.line(steady.refFine.map((c) => c.re), zzF, { color: '#7EE7E8', width: 5, alpha: 0.4, label: '解析解 u, v' });
    f.line(steady.refFine.map((c) => c.im), zzF, { color: '#F0B429', width: 5, alpha: 0.4 });
    const [tu, tv, tz] = thin([us, vs, zz]);
    f.points(tu, tz, { color: '#7EE7E8', r: 2.6, label: '数値解 u（東西）' });
    f.points(tv, tz, { color: '#F0B429', r: 2.6, label: '数値解 v（南北）' });
    f.legend({ corner: 'br' });
  }

  // speed against depth, logarithmic
  {
    const f = new Frame($('c-prof2'), { pad: { l: 60, r: 18, t: 16, b: 40 } });
    const zz = steady.z.map((z) => z / P.DE);
    const zzF = steady.zFine.map((z) => z / P.DE);
    const tiny = 1e-12;
    const sp = [];
    for (let k = 0; k < P.N; k++) {
      sp.push(Math.max(tiny, Math.hypot(steady.W.re[k], steady.W.im[k])));
    }
    const spa = steady.refFine.map((c) => Math.max(tiny, Math.hypot(c.re, c.im)));
    const smax = Math.max(...sp, ...spa);
    f.begin({
      xlim: [smax * Math.exp(-state.hd) * 0.25, smax * 1.6], ylim: [-state.hd, 0], xlog: true,
      xlabel: '流速の大きさ [m/s]', ylabel: '深さ z / D_E',
    });
    // Pure exponential decay, anchored at the analytic surface speed.
    f.line(zzF.map((z) => spa[0] * Math.exp(z)), zzF,
      { color: palette().muted, width: 1.5, dash: [5, 4], label: '純粋な exp(z / D_E)' });
    f.line(spa, zzF, { color: '#F0B429', width: 5, alpha: 0.45, label: '解析解' });
    const [ts, tz] = thin([sp, zz]);
    f.points(ts, tz, { color: '#7EE7E8', r: 2.6, label: '数値解' });
    f.legend({ corner: 'br' });
  }

  // rotation angle against depth
  {
    const f = new Frame($('c-prof3'), { pad: { l: 60, r: 18, t: 16, b: 40 } });
    const zz = steady.z.map((z) => z / P.DE);
    const zzF = steady.zFine.map((z) => z / P.DE);
    const num = unwrapDeg(steady.z.map((_, k) =>
      A.deflectionDeg({ re: steady.W.re[k], im: steady.W.im[k] }, P)));
    const ana = unwrapDeg(steady.refFine.map((c) => A.deflectionDeg(c, P)));
    // Ideal spiral: -45 deg at the surface, one radian per Ekman depth.
    const sgn = Math.sign(P.f);
    const ideal = zzF.map((z) => sgn * (-45 + z * (180 / Math.PI)));
    const [lo, hi] = span([num, ana, ideal]);
    f.begin({ xlim: [lo, hi], ylim: [-state.hd, 0], xlabel: '風からの角度 [度]', ylabel: '深さ z / D_E' });
    f.line(ideal, zzF, { color: palette().muted, width: 1.5, dash: [5, 4], label: '1 D_E あたり 1 ラジアン' });
    f.line(ana, zzF, { color: '#F0B429', width: 5, alpha: 0.45, label: '解析解' });
    const [tn, tz] = thin([num, zz]);
    f.points(tn, tz, { color: '#7EE7E8', r: 2.6, label: '数値解' });
    f.legend({ corner: 'tr' });
    $('cap-prof3').innerHTML =
      `風応力の向きを 0° とした角度。±180° で折り返さないように連続化してある。` +
      `海面の ${A.deflectionDeg(steady.d.surface, P).toFixed(1)}° から始まり、` +
      `深さに対して線形に、1 エクマン深さあたりちょうど 1 ラジアン（57.3°）ずつ回っていく。`;
  }
}

function drawSpinup() {
  const fr = run.frames[Math.min(state.frame, run.frames.length - 1)];
  hodograph($('c-anim'), fr, { ghost: steady.W, showAnalytic: false, title: false });

  const i = Math.min(
    run.t.length - 1,
    Math.round((fr.t / (run.t[run.t.length - 1] || 1)) * (run.t.length - 1)),
  );

  // Trajectory of the depth-integrated transport
  {
    const f = new Frame($('c-traj'), { pad: { l: 66, r: 20, t: 16, b: 42 } });
    const [xlo, xhi] = span([run.Mx, run.MxeFine, [0, steady.Mex.re]], 0.14);
    const [ylo, yhi] = span([run.My, run.MyeFine, [0, steady.Mex.im]], 0.14);
    f.begin({
      xlim: [xlo, xhi], ylim: [ylo, yhi], equal: true,
      xlabel: '東西成分 M_x [kg/(m·s)]', ylabel: '南北成分 M_y [kg/(m·s)]',
    });
    f.line(run.MxeFine, run.MyeFine, { color: '#F0B429', width: 5, alpha: 0.45, label: '厳密解' });
    f.line(run.Mx.slice(0, i + 1), run.My.slice(0, i + 1), { color: '#7EE7E8', width: 1.8, label: '数値解' });
    f.points([steady.Mex.re], [steady.Mex.im], { color: '#F0B429', r: 5, label: 'エクマン輸送 τ/f' });
    f.points([run.Mx[i]], [run.My[i]], { color: '#E9F2F4', r: 4 });
    f.legend({ corner: f.pickCorner([...run.MxeFine, steady.Mex.re], [...run.MyeFine, steady.Mex.im]) });
  }

  const relerr = Math.hypot(run.Mx[i] - run.Mxe[i], run.My[i] - run.Mye[i]) /
    (Math.hypot(steady.Mex.re, steady.Mex.im) || 1);
  $('cap-traj').innerHTML =
    `静止した海に t = 0 から風を吹かせ続けると、輸送は定常のエクマン輸送のまわりを回る慣性振動を残す。` +
    `海底を応力なしにすると深さ平均流には摩擦が働かないので、この振動は<b>いつまでも減衰しない</b>。` +
    `実際の海洋でも、強風のあとの漂流ブイの軌跡に、慣性周期でループを描く運動が観測されている。` +
    `現時点での数値解と厳密解の差は相対で` +
    `<b>${Number.isFinite(relerr) ? relerr.toExponential(2) : '発散'}</b>。`;

  {
    const f = new Frame($('c-mts'), { pad: { l: 68, r: 16, t: 14, b: 38 } });
    const th = run.t.map((t) => t / P.inertial);
    const [lo, hi] = span([run.Mx, run.My]);
    f.begin({ xlim: [0, th[th.length - 1]], ylim: [lo, hi], xlabel: '時間 / 慣性周期', ylabel: '輸送 M [kg/(m·s)]' });
    const thF = run.tFine.map((t) => t / P.inertial);
    f.line(thF, run.MxeFine, { color: '#E9F2F4', width: 5, alpha: 0.28, label: '厳密解' });
    f.line(thF, run.MyeFine, { color: '#E9F2F4', width: 5, alpha: 0.28 });
    f.line(th, run.Mx, { color: '#7EE7E8', width: 1.8, label: '数値解 M_x' });
    f.line(th, run.My, { color: '#F0B429', width: 1.8, label: '数値解 M_y' });
    f.legend({ corner: 'tr' });
  }

  $('clock').textContent = `t = ${hhmm(fr.t)}（慣性周期の ${(fr.t / P.inertial).toFixed(2)} 倍）`;
  $('tslider').value = String(state.frame);
}

function drawBudget() {
  const th = run.t.map((t) => t / P.inertial);
  const balance = run.cumP.map((c, i) => c - run.cumE[i]);

  {
    const f = new Frame($('c-energy'), { pad: { l: 70, r: 16, t: 16, b: 40 } });
    const [lo, hi] = span([balance, run.E]);
    f.begin({
      xlim: [0, th[th.length - 1]], ylim: [Math.min(0, lo), hi],
      xlabel: '時間 / 慣性周期', ylabel: 'エネルギー [J/m²]',
    });
    f.line(th, balance, { color: '#F0B429', width: 5, alpha: 0.65, label: '∫風の仕事 − ∫散逸' });
    f.line(th, run.E, { color: '#7EE7E8', width: 1.8, label: '層内の運動エネルギー' });
    f.legend({ corner: 'br' });
    const n = th.length - 1;
    $('cap-energy').innerHTML = run.blewUp
      ? `この設定では計算が発散しているので、以下の数値は物理としては意味を持たない。` +
        `それでも 2 本の曲線は重なったままである。収支の恒等式は、スキームが正しかろうと` +
        `間違っていようと、それが計算した値について成り立つからである。` +
        `<b>収支が閉じることは、コードが自分の方程式を忠実に解いている証明であって、` +
        `その方程式が正しい答えを与える証明ではない。</b>後者を確かめるために解析解との比較がある。`
      : `静止した海に t = 0 から風を吹かせ続けたときの、層内のエネルギーの時間発展である。` +
        `この間に風は <b>${num(run.cumP[n], 0)} J/m²</b> の仕事をし、粘性は ` +
        `<b>${num(run.cumE[n], 0)} J/m²</b> の力学的エネルギーを熱に変えて取り去った。差の ` +
        `<b>${num(run.cumP[n] - run.cumE[n])} J/m²</b> は、層内に実際に存在する` +
        `運動エネルギーと一致する。2 本の曲線が重なっているのがそれである。`;
  }

  {
    const f = new Frame($('c-resid'), { pad: { l: 70, r: 16, t: 16, b: 40 } });
    const eps = 1e-18;
    const rm = run.resM.map((v) => Math.max(v, eps));
    const re = run.resE.map((v) => Math.max(v, eps));
    const top = Math.max(...rm, ...re, 1e-16);
    const th2 = run.tRes.map((t) => t / P.inertial);
    f.begin({
      xlim: [0, th2[th2.length - 1] || 1], ylim: [1e-18, Math.max(top * 30, 1e-12)], ylog: true,
      xlabel: '時間 / 慣性周期', ylabel: '相対残差',
    });
    f.line(th2, rm, { color: '#7EE7E8', width: 1.8, label: '運動量収支' });
    f.line(th2, re, { color: '#F0B429', width: 1.8, label: 'エネルギー収支' });
    f.line([0, th2[th2.length - 1] || 1], [2.2e-16, 2.2e-16],
      { color: palette().muted, width: 1.4, dash: [4, 4], label: '倍精度の機械イプシロン' });
    f.legend({ corner: 'tr' });
  }

  const worst = Math.max(...run.resM, ...run.resE);
  $('cap-resid').innerHTML =
    `運動量保存とエネルギー保存は、差分方程式の解析解が厳密に満たす等式である。` +
    `風応力をセル界面のフラックスの値として与えているため、セルを足し合わせると内部の応力が` +
    `厳密に打ち消え、残差は打切り誤差ではなく丸め誤差の水準に留まる。` +
    `この計算での最悪値は <b>${worst.toExponential(2)}</b>` +
    (state.sch === 'euler'
      ? '。前進Euler法でも収支は同じように閉じることに注意。収支が閉じるのは、'
        + 'コードが自分の方程式を忠実に解いていることの証明にすぎない。'
      : '。');

  {
    const f = new Frame($('c-power'), { pad: { l: 74, r: 16, t: 14, b: 38 } });
    const [lo, hi] = span([run.Pw, run.eps, [steady.Pex]]);
    f.begin({ xlim: [0, th[th.length - 1]], ylim: [Math.min(0, lo), hi], xlabel: '時間 / 慣性周期', ylabel: '仕事率 [W/m²]' });
    f.line([0, th[th.length - 1]], [steady.Pex, steady.Pex],
      { color: palette().muted, width: 1.4, dash: [5, 4], label: '解析値 |τ|²/(ρ₀ f D_E)' });
    f.line(th, run.Pw, { color: '#F0B429', width: 2, label: '風の仕事率' });
    f.line(th, run.eps, { color: '#E8705A', width: 2, label: '粘性散逸' });
    f.legend({ corner: 'br' });
  }
}

function drawPump() {
  const ykm = Array.from(pump.y, (v) => (v - 2.0e6) / 1000);
  const tx = pump.stresses.map((s) => s.taux);

  {
    // Two quantities in different units, so each gets its own axis: transport on
    // the left, stress on the right. Both ranges are made symmetric about zero
    // so that the two zero levels coincide and the curves can be read together.
    const f = new Frame($('c-pump1'), { pad: { l: 64, r: 16, t: 16, b: 40 } });
    const V = Array.from(pump.V);
    const va = Math.max(...V.map(Math.abs), 1e-12) * 1.08;
    const ta = Math.max(...tx.map(Math.abs), 1e-12) * 1.08;
    f.begin({
      xlim: [ykm[0], ykm[ykm.length - 1]], ylim: [-va, va],
      xlabel: '南北方向の位置 [km]', ylabel: '北向き輸送 V_E [m²/s]',
      y2lim: [-ta, ta], y2label: '東向き風応力 τ_x [N/m²]',
      xends: ['南 \u2190', '\u2192 北'],
    });
    f.line(ykm, tx, { color: '#F0B429', width: 2.2, label: '東向き風応力 τ_x（右軸）', axis: 'y2' });
    f.line(ykm, V, { color: '#7EE7E8', width: 2.2, label: '北向き輸送 V_E（左軸）' });
    f.legend({ corner: 'tr' });
  }

  {
    const f = new Frame($('c-pump2'), { pad: { l: 70, r: 16, t: 16, b: 40 } });
    const wy = Array.from(pump.we, (v) => v * 86400 * 365);
    const we = Array.from(pump.weExact, (v) => v * 86400 * 365);
    const [lo, hi] = span([wy, we]);
    f.begin({
      xlim: [ykm[0], ykm[ykm.length - 1]], ylim: [lo, hi],
      xlabel: '南北方向の位置 [km]', ylabel: '鉛直流速 w_e [m/年]',
      xends: ['南 \u2190', '\u2192 北'],
    });
    f.line(ykm, we, { color: '#F0B429', width: 5.5, alpha: 0.55, label: '解析解 curl(τ/f)/ρ₀' });
    f.line(ykm, wy, { color: '#63C9A8', width: 1.8, label: '数値解' });
    f.legend({ corner: 'tr' });
  }

  const peak = Math.max(...Array.from(pump.weExact, Math.abs)) * 86400 * 365;
  $('cap-pump').innerHTML =
    `正が湧昇（下から汲み上がる向き）。最大で <b>${peak.toFixed(0)} m/年</b>で、` +
    `実際の海で観測されるエクマン湧昇と同じオーダーである。南北に周期境界を課しているので、` +
    `1 周期にわたって足すと湧昇と沈降は厳密に相殺しなければならない。正味の流量は` +
    `<b>${pump.net === 0 ? 'ちょうどゼロ' : Math.abs(pump.net).toExponential(1) + ' m²/s'}</b>になっている。`;
}

// ---------------------------------------------------------------------------
// Verification tab
// ---------------------------------------------------------------------------
/** The parameter set a verification run used, spelled out for the reader. */
function verifyParams(cfg) {
  const s = Math.max(-1, Math.min(1, cfg.f / (2 * OMEGA_EARTH)));
  const lat = (Math.asin(s) * 180) / Math.PI;
  const mag = Math.hypot(cfg.taux, cfg.tauy);
  const dir = ((Math.atan2(cfg.tauy, cfg.taux) * 180) / Math.PI + 360) % 360;
  return [
    `${lat >= 0 ? '北緯' : '南緯'} ${Math.abs(lat).toFixed(0)}°（f = ${cfg.f.toExponential(3)} s⁻¹）`,
    `渦粘性 ν_z = ${cfg.nu.toExponential(1)} m²/s`,
    `密度 ρ₀ = ${cfg.rho0} kg/m³`,
    `風応力 |τ| = ${mag.toFixed(3)} N/m²（${COMPASS[Math.round(dir / 22.5) % 16]}向き）`,
    `海底 ${cfg.bottom === 'noslip' ? '粘着' : '応力なし'}`,
    `エクマン深さ D_E = ${cfg.DE.toFixed(2)} m`,
    `慣性周期 ${cfg.inertialHours.toFixed(1)} 時間`,
  ].join(' ・ ');
}

function renderVerification(r, live = false, seconds = null) {
  const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;');
  const parts = [];
  parts.push(
    `<div class="summary ${r.passed === r.total ? 'pass' : 'fail'}">` +
    `${r.total} 項目中 ${r.passed} 項目が合格` +
    `<div style="font:400 13px/1.8 var(--jp);color:var(--haze);margin-top:6px">` +
    `<b style="color:var(--foam);font-weight:500">` +
    `${live ? '左のスライダーの値で実行' : '基準ケースで実行'}</b>` +
    (seconds == null ? '' : `（計算 ${seconds.toFixed(1)} 秒）`) + '<br>' +
    esc(verifyParams(r.config)) + '<br>' +
    '格子数 N・計算領域の深さ H・時間刻み Δt は項目ごとに変えている（収束次数を測る項目は' +
    'これらを系統的に振る）。' +
    (live
      ? '海底の条件は、解析解との比較が定義できる「応力なし」に固定して走る。'
      : 'この設定は固定なので、結果は誰が実行しても同じになる。') +
    `</div></div>`,
  );
  for (const g of r.groups) {
    let html = `<div class="vgroup"><h3>${esc(g.title)}</h3>`;
    if (g.columns) {
      html += '<div class="tablewrap"><table><thead><tr>' +
        g.columns.map((c) => `<th>${esc(c)}</th>`).join('') +
        '</tr></thead><tbody>' +
        g.rows.map((row) => '<tr>' + row.map((v) => `<td>${esc(v)}</td>`).join('') + '</tr>').join('') +
        '</tbody></table></div>';
    }
    for (const n of g.notes) html += `<p class="vnote">${esc(n)}</p>`;
    for (const c of g.checks) {
      html += `<div class="vcheck"><span class="badge ${c.pass ? 'pass' : 'fail'}">` +
        `${c.pass ? '合格' : '不合格'}</span><span>${esc(c.name)}<br>` +
        `<code>結果 ${esc(c.value)} · 期待 ${esc(c.expected)}</code></span></div>`;
    }
    $('vout').insertAdjacentHTML('beforeend', html + '</div>');
    html = '';
  }
  $('vout').insertAdjacentHTML('afterbegin', parts.join(''));
}

// ---------------------------------------------------------------------------
// Rendering and events
// ---------------------------------------------------------------------------
function render() {
  if (!steady) return;
  if (state.tab === 'spiral') drawSpiral();
  else if (state.tab === 'spinup') drawSpinup();
  else if (state.tab === 'budget') drawBudget();
  else if (state.tab === 'pump') drawPump();
}

function selectTab(tab) {
  state.tab = tab;
  for (const b of $('tabs').children) b.setAttribute('aria-selected', String(b.dataset.tab === tab));
  for (const id of ['spiral', 'spinup', 'budget', 'pump', 'verify']) {
    $(`panel-${id}`).hidden = id !== tab;
  }
  render();
}

for (const k of Object.keys(DEFAULTS)) {
  $(k).addEventListener('input', () => { state.reveal = 1; recompute(); });
}
$('reset').addEventListener('click', () => {
  for (const [k, v] of Object.entries(DEFAULTS)) $(k).value = String(v);
  state.frame = 0; state.reveal = 1;
  recompute();
});
$('tabs').addEventListener('click', (e) => {
  const b = e.target.closest('button[data-tab]');
  if (b) selectTab(b.dataset.tab);
});
const syncPlayLabel = () => {
  $('play').textContent = state.playing ? '一時停止' : '再生';
};
$('play').addEventListener('click', () => {
  state.playing = !state.playing;
  syncPlayLabel();
});
$('tslider').addEventListener('input', (e) => {
  state.playing = false;
  syncPlayLabel();
  state.frame = parseInt(e.target.value, 10);
  drawSpinup();
});
/**
 * Run the suite and render it. `base` is undefined for the fixed reference case
 * and an explicit parameter set for the "use the sliders" mode. The bottom
 * condition is pinned to free-slip there: several of the reference solutions the
 * suite compares against (the τ/f transport, the spin-up transport) assume no
 * bottom stress, so a no-slip run would report failures that are physics rather
 * than bugs.
 */
function runVerification(btn, base, live) {
  const buttons = [$('runv'), $('runv-live')];
  for (const b of buttons) b.disabled = true;
  const label = btn.textContent;
  btn.textContent = '計算中…';
  $('vout').innerHTML = '';
  // Yield once so the button repaints before the synchronous run starts.
  setTimeout(() => {
    const t0 = performance.now();
    const r = runChecks(base);
    renderVerification(r, live, (performance.now() - t0) / 1000);
    for (const b of buttons) b.disabled = false;
    btn.textContent = label;
  }, 30);
}
$('runv').addEventListener('click', () => runVerification($('runv'), undefined, false));
$('runv-live').addEventListener('click', () => runVerification($('runv-live'), {
  f: P.f, nu: P.nu, rho0: P.rho0, taux: P.taux, tauy: P.tauy, bottom: 'free',
}, true));
window.addEventListener('resize', () => render());

// ---------------------------------------------------------------------------
// Animation loop: one reveal of the steady spiral on load, then the spin-up
// playback whenever that tab is open.
// ---------------------------------------------------------------------------
let last = performance.now();
function tick(now) {
  const dt = Math.min(0.05, (now - last) / 1000);
  last = now;
  if (state.reveal < 1) {
    state.reveal = Math.min(1, state.reveal + dt / 0.9);
    if (state.tab === 'spiral') drawSpiral();
  } else if (state.tab === 'spinup' && state.playing && run) {
    state.frame = (state.frame + 1) % run.frames.length;
    drawSpinup();
  }
  requestAnimationFrame(tick);
}

recompute();
selectTab('spiral');
syncPlayLabel();
requestAnimationFrame(tick);

// Canvas text is rasterised at draw time, so the very first frame would use a
// fallback face if the web fonts have not finished loading. Redraw once they have.
if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => render());
