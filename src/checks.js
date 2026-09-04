// The verification suite, as data.
//
// Both `node tests/verify.mjs` and the "Verification" tab of the web page call
// runChecks() and only differ in how they render the result. That is deliberate:
// the table published in the README is then a statement about the code that is
// actually deployed, not about a separate throwaway script.

import { Integrator, solveSteady, cellCentres, ekmanDepth, coriolis } from './solver.js';
import * as A from './analytic.js';
import * as D from './diagnostics.js';
import { solveUpwelling } from './upwelling.js';
import { abs, sub } from './complex.js';

export const DEFAULT_BASE = {
  f: coriolis(45),
  nu: 1.0e-2,
  rho0: 1025,
  taux: 0.1, // eastward: the westerlies blow from the west
  tauy: 0.0,
  bottom: 'free',
};

const e3 = (x) => x.toExponential(3);
const l2 = (num, ref) => {
  let s = 0;
  for (let k = 0; k < ref.length; k++) s += (num.re[k] - ref[k].re) ** 2 + (num.im[k] - ref[k].im) ** 2;
  return Math.sqrt(s / ref.length);
};
const linf = (num, ref) => {
  let m = 0;
  for (let k = 0; k < ref.length; k++) m = Math.max(m, Math.hypot(num.re[k] - ref[k].re, num.im[k] - ref[k].im));
  return m;
};
const order = (prev, cur) => (prev == null ? null : Math.log2(prev / cur));
const ord = (o) => (o == null ? '--' : o.toFixed(3));

export function runChecks(BASE = DEFAULT_BASE) {
  const DE = ekmanDepth(BASE);
  const groups = [];
  const checks = [];
  let id = 0;
  const add = (group, name, value, expected, pass) => {
    id += 1;
    const c = { id, name, value: String(value), expected, pass: !!pass };
    checks.push(c);
    group.checks.push(c);
    return c;
  };
  const G = (title, columns) => {
    const g = { title, columns: columns || null, rows: [], notes: [], checks: [] };
    groups.push(g);
    return g;
  };

  // -- 1. Grid convergence of the steady solution ---------------------------
  {
    const g = G('定常解の格子収束（格子を細かくすると解析解に近づくか）',
      ['格子数 N', 'Δz / D_E', 'L2 誤差 [m/s]', '次数', 'L∞ 誤差 [m/s]', '次数']);
    g.notes.push(
      '比較対象は、無限深の古典解ではなく、シミュレーターが解いているのと同じ有限深さ境界値問題の' +
      '厳密解である。こうしないと、領域を有限で打ち切ったことによる誤差が離散化誤差に化けて、' +
      '収束次数が汚れてしまう。');
    const H = 8 * DE;
    let p2 = null, pi = null, last = null;
    for (const N of [24, 48, 96, 192, 384, 768]) {
      const p = { ...BASE, H, N };
      const W = solveSteady(p);
      const ref = A.steadyProfile(p, Array.from(cellCentres(p)));
      const e2 = l2(W, ref), ei = linf(W, ref);
      const o2 = order(p2, e2), oi = order(pi, ei);
      g.rows.push([N, (H / N / DE).toFixed(4), e3(e2), ord(o2), e3(ei), ord(oi)]);
      p2 = e2; pi = ei; last = o2;
    }
    add(g, '空間精度の次数（格子幅を半分にすると誤差が 1/4 になるか）', last.toFixed(4), '2.00 ± 0.05', Math.abs(last - 2) < 0.05);
  }

  // -- 2. Surface deflection -------------------------------------------------
  {
    // The surface current sits 45 degrees to the right of the wind in the
    // northern hemisphere and 45 degrees to the left in the southern one.
    const defl0 = -45 * Math.sign(BASE.f);
    const g = G(`海面の流れは風の 45° ${BASE.f >= 0 ? '右' : '左'}向きか`,
      ['格子数 N', '偏角 [度]', '誤差の絶対値 [度]', '次数']);
    const H = 12 * DE;
    let prev = null, lastOrder = null, d1 = null, d0 = null;
    for (const N of [48, 96, 192, 384, 768]) {
      const p = { ...BASE, H, N };
      const s = D.surfaceVelocityExtrapolated(solveSteady(p), p);
      const defl = A.deflectionDeg(s, p);
      const err = Math.abs(defl - defl0);
      const o = order(prev, err);
      g.rows.push([N, defl.toFixed(6), e3(err), ord(o)]);
      prev = err; lastOrder = o; d0 = d1; d1 = defl;
    }
    const rich = (4 * d1 - d0) / 3; // second-order Richardson extrapolation
    g.notes.push(
      `最も細かい 2 つの格子から 2 次の Richardson 外挿を行うと ${rich.toFixed(8)} 度になる。` +
      '有限の格子で得た値から、格子幅ゼロの極限を推定する操作である。');
    g.notes.push(
      '「偏角の収束次数」は、格子数 N を倍にするたびに −45° からの誤差が何分の 1 になるかを' +
      '表す指標である。次数 2 は誤差が 1/4 になること、すなわち海面付近の離散化が' +
      '設計どおり 2 次精度で効いていることを意味する。');
    add(g, '海面の偏角（Richardson 外挿値）', rich.toFixed(6) + ' 度',
      defl0.toFixed(6), Math.abs(rich - defl0) < 1e-5);
    add(g, '偏角の収束次数', lastOrder.toFixed(4), '2.00 ± 0.10',
      Math.abs(lastOrder - 2) < 0.1);
  }

  // -- 3. Depth-integrated transport ----------------------------------------
  {
    const g = G('鉛直積分したエクマン輸送（層全体で運ばれる水の量）', ['', '東西成分 M_x [kg/(m·s)]', '南北成分 M_y [kg/(m·s)]']);
    const p = { ...BASE, H: 8 * DE, N: 256 };
    const M = D.transport(solveSteady(p), p);
    const Mex = A.transportClassical(p);
    const rel = abs(sub(M, Mex)) / abs(Mex);
    g.rows.push(['数値解', M.re.toExponential(6), M.im.toExponential(6)]);
    g.rows.push(['解析解 (τ_y − iτ_x)/f', Mex.re.toExponential(6), Mex.im.toExponential(6)]);
    g.notes.push(
      `輸送の向きは東から ${((Math.atan2(M.im, M.re) * 180) / Math.PI).toFixed(4)} 度、` +
      'つまり東向きの風のちょうど 90° 右である。海底を応力なしにすると、境界を横切る応力は' +
      '上端の風応力と下端のゼロだけになる。そのため、この等式は格子の粗さにも領域の深さにも' +
      'よらず、離散のレベルで厳密に成り立つ。');
    add(g, '輸送が (τ_y − iτ_x)/f と一致する', e3(rel), '< 1e-12（厳密に成立）', rel < 1e-12);

    let worst = 0;
    for (const N of [8, 16, 32, 64, 128, 512]) {
      for (const hf of [3, 5, 8, 12]) {
        const q = { ...BASE, H: hf * DE, N };
        worst = Math.max(worst, abs(sub(D.transport(solveSteady(q), q), Mex)) / abs(Mex));
      }
    }
    add(g, '輸送が格子数 N と深さ H によらない', e3(worst), '24 通りすべてで < 1e-11', worst < 1e-11);
  }

  // -- 4. Domain truncation with a no-slip bottom ----------------------------
  {
    const g = G('海底を粘着条件にしたときの領域打ち切り誤差（有限の深さで切ってよいか）',
      ['H / D_E', '輸送の相対誤差', '予測値 2exp(−H/D_E)', '有限深さ厳密解との差']);
    g.notes.push(
      '海底に粘着条件を課すと底面摩擦が運動量を抜き取るため、輸送は古典値より sech(λH) だけ' +
      '足りなくなる。この誤差が深さとともに指数関数的に減衰することが、本来は無限に深い海を' +
      '有限の領域で打ち切ってよい根拠になる。');
    let worstRatio = 0, worstFinite = 0;
    for (const hf of [1, 2, 3, 4, 5, 6, 8]) {
      const p = { ...BASE, H: hf * DE, N: Math.max(64, 32 * hf), bottom: 'noslip' };
      const M = D.transport(solveSteady(p), p);
      const Mex = A.transportClassical(p);
      const rel = abs(sub(M, Mex)) / abs(Mex);
      const pred = 2 * Math.exp(-hf);
      const fin = abs(sub(M, A.transportExact(p))) / abs(Mex);
      if (hf >= 3) {
        worstRatio = Math.max(worstRatio, Math.abs(rel / pred - 1));
        worstFinite = Math.max(worstFinite, fin);
      }
      g.rows.push([hf, e3(rel), e3(pred), e3(fin)]);
    }
    add(g, '打ち切り誤差が 2exp(−H/D_E) に従う', e3(worstRatio), 'H ≥ 3 D_E で相対 1% 以内', worstRatio < 0.01);
    add(g, '有限深さの厳密式と一致する', e3(worstFinite), '< 1e-4（空間離散化誤差 O(Δz²)）', worstFinite < 1e-4);
  }

  // -- 5. Discrete budget residuals -----------------------------------------
  {
    const g = G('立ち上がり計算中の離散収支の残差（何かが勝手に増減していないか）', ['診断量', '計算全体での最悪値']);
    const p = { ...BASE, H: 8 * DE, N: 200 };
    const it = new Integrator(p, { scheme: 'cn' });
    const Pref = A.steadyPower(p);
    let maxMom = 0, maxEne = 0, maxCor = 0, minEps = Infinity;
    const tEnd = 6 * 86400;
    while (it.t < tEnd) {
      const { Wold, Wmid, Wnew } = it.step();
      maxMom = Math.max(maxMom, D.momentumResidual(Wold, Wnew, Wmid, p, it.dt).relative);
      const er = D.energyResidual(Wold, Wnew, Wmid, p, it.dt);
      maxEne = Math.max(maxEne, Math.abs(er.relative));
      maxCor = Math.max(maxCor, Math.abs(D.coriolisWork(Wmid, p)) / Pref);
      minEps = Math.min(minEps, er.eps);
    }
    g.rows.push(['運動量収支の残差（相対）', e3(maxMom)]);
    g.rows.push(['エネルギー収支の残差（相対）', e3(maxEne)]);
    g.rows.push(['コリオリ力の仕事 / 風の入力', e3(maxCor)]);
    g.rows.push(['粘性散逸の最小値 [W/m²]', e3(minEps)]);
    g.notes.push(
      `${(tEnd / 86400).toFixed(0)} 日間を ${it.dt.toFixed(1)} 秒刻みで ${it.nstep} ステップ計算した。` +
      'これらは離散化した方程式が自分自身の収支を守っているかを測る量なので、' +
      '丸め誤差の水準に留まり続けるべきものである。時間とともに増えるようなら、' +
      'それは近似の限界ではなくバグである。');
    add(g, '運動量収支が閉じる', e3(maxMom), '< 1e-12', maxMom < 1e-12);
    add(g, 'エネルギー収支が閉じる', e3(maxEne), '< 1e-10', maxEne < 1e-10);
    add(g, 'コリオリ力が仕事をしない', e3(maxCor), '< 1e-14', maxCor < 1e-14);
    add(g, '粘性散逸が負にならない', e3(minEps), 'すべてのステップで 0 以上', minEps >= 0);
  }

  // -- 6. Time convergence ---------------------------------------------------
  {
    const g = G('慣性振動の厳密解に対する時間収束（時間刻みを細かくすると厳密解に近づくか）',
      ['ステップ数', 'Δt [s]', '輸送 M の相対誤差', '次数']);
    const p = { ...BASE, H: 8 * DE, N: 128 };
    const tEnd = 2.5 * ((2 * Math.PI) / Math.abs(p.f));
    const Mex = A.spinupTransport(p, tEnd);
    g.notes.push(
      '運動方程式を層全体で鉛直積分すると粘性項が完全に消え、dM/dt = T − ifM だけが残る。' +
      'その解 −i(T/f)(1 − exp(−ift)) には ν も H も格子も含まれないので、' +
      '時間積分の誤差だけを切り出して測ることができる。');
    let prev = null, lastOrder = null;
    for (const nsteps of [50, 100, 200, 400, 800, 1600]) {
      const dt = tEnd / nsteps;
      const it = new Integrator(p, { scheme: 'cn', dt });
      for (let n = 0; n < nsteps; n++) it.step();
      const rel = abs(sub(D.transport(it.W, p), Mex)) / abs(Mex);
      const o = order(prev, rel);
      g.rows.push([nsteps, dt.toFixed(2), e3(rel), ord(o)]);
      prev = rel; lastOrder = o;
    }
    add(g, '時間精度の次数（時間刻みを半分にすると誤差が 1/4 になるか）', lastOrder.toFixed(4), '2.00 ± 0.05', Math.abs(lastOrder - 2) < 0.05);

    const ref = [];
    for (const N of [16, 64, 256]) {
      const q = { ...BASE, H: 8 * DE, N };
      const itq = new Integrator(q, { scheme: 'cn', dt: tEnd / 200 });
      for (let n = 0; n < 200; n++) itq.step();
      ref.push(abs(sub(D.transport(itq.W, q), A.spinupTransport(q, tEnd))) / abs(Mex));
    }
    const spread = (Math.max(...ref) - Math.min(...ref)) / Math.max(...ref);
    g.notes.push(`実際、格子数 N = 16, 64, 256 のいずれでも誤差は ${ref.map((x) => x.toExponential(4)).join('、')} と一致する。`);
    add(g, '誤差が鉛直格子の細かさによらない', e3(spread), '< 1e-9', spread < 1e-9);
  }

  // -- 7. Steady energy balance ---------------------------------------------
  {
    const g = G('定常状態でのエネルギー収支（風が入れる仕事と粘性が奪う仕事の釣り合い）',
      ['格子数 N', '風の入力 P [W/m²]', '粘性散逸 ε [W/m²]', '|P − ε| / P', 'P と解析値の差']);
    const H = 12 * DE;
    let worstBalance = 0, lastErr = null, prevErr = null, lastOrder = null;
    for (const N of [64, 128, 256, 512, 1024]) {
      const p = { ...BASE, H, N };
      const W = solveSteady(p);
      const Pd = D.windPowerDiscrete(W, p);
      const ed = D.dissipation(W, p);
      const err = Math.abs(D.windPowerPhysical(W, p) - A.steadyPower(p)) / A.steadyPower(p);
      const bal = Math.abs(Pd - ed) / Pd;
      worstBalance = Math.max(worstBalance, bal);
      prevErr = lastErr; lastErr = err;
      if (prevErr != null) lastOrder = order(prevErr, err);
      g.rows.push([N, Pd.toExponential(6), ed.toExponential(6), e3(bal), e3(err)]);
    }
    const p = { ...BASE, H, N: 1024 };
    const W = solveSteady(p);
    const Pp = D.windPowerPhysical(W, p), ep = D.dissipationPhysical(W, p);
    const Pex = A.steadyPower(p);
    g.notes.push(
      'ここでの「定常」は、同じ風が吹き続けて流れが時間変化しなくなった状態（∂W/∂t = 0）を指す。' +
      `定常での解析的な入力は |τ|²/(ρ₀ f D_E) = ${Pex.toExponential(6)} W/m² である。` +
      'この式に隠れている cos(45°) は、上で測った海面の偏角そのものなので、' +
      '2 つの結果は同じ事実を別の角度から見ているにすぎない。');
    g.notes.push(
      `海面直下の半セル分を両辺に加えると、N = 1024 で P = ${Pp.toExponential(6)}、` +
      `ε = ${ep.toExponential(6)} W/m² となり完全に一致する。`);
    g.notes.push(
      '「P が |τ|²/(ρ₀ f D_E) に収束する」の次数も同じ読み方で、格子数 N を倍にすると' +
      '風の仕事率が解析値へ 1/4 ずつ近づくことを表す。');
    add(g, '定常状態で風の入力と粘性散逸が等しい（離散）', e3(worstBalance), '< 1e-12', worstBalance < 1e-12);
    add(g, '海面の半セル分を補正した入力と散逸が一致する', e3(Math.abs(Pp - ep) / Pex), '< 1e-12', Math.abs(Pp - ep) / Pex < 1e-12);
    add(g, 'P が |τ|²/(ρ₀ f D_E) に収束する', lastOrder.toFixed(4), '2.00 ± 0.10', Math.abs(lastOrder - 2) < 0.1);
  }

  // -- 8. Rotation must not create energy ------------------------------------
  {
    const g = G('回転がエネルギーを作ってはいけない（陰解法と陽解法の比較）', ['スキーム', 'エネルギーの相対変化', '理論予測']);
    const p = { ...BASE, nu: 0, taux: 0, tauy: 0, H: 100, N: 40 };
    const nsteps = 500;
    const dt = (5 * ((2 * Math.PI) / Math.abs(p.f))) / nsteps;
    const out = {};
    for (const scheme of ['cn', 'euler']) {
      const it = new Integrator(p, { scheme, dt });
      for (let k = 0; k < p.N; k++) { it.W.re[k] = 0.1; it.W.im[k] = 0.0; }
      const E0 = D.energy(it.W, p);
      for (let n = 0; n < nsteps; n++) it.step();
      out[scheme] = D.energy(it.W, p) / E0 - 1;
    }
    const predicted = (1 + (p.f * dt) ** 2) ** nsteps - 1;
    g.rows.push(['Crank–Nicolson法', e3(out.cn), '厳密に 0']);
    g.rows.push(['前進Euler法', e3(out.euler), e3(predicted)]);
    g.notes.push(
      '粘性も風もない状況では、コリオリ力は速度ベクトルの向きを変えるだけなので、' +
      'エネルギーは一定でなければならない。Crank–Nicolson法が施す Cayley 変換は絶対値が' +
      '厳密に 1 である。一方で前進Euler法は (1 − ifΔt) を掛けるが、その絶対値は 1 を超えるため、' +
      'エネルギーを勝手に作り出す。しかもその増え方は理論予測と 14 桁一致する。' +
      'このシミュレーターを陰解法にした理由がこれである。');
    add(g, 'Crank–Nicolson法はコリオリ回転だけを解いたときエネルギーを保存する', e3(Math.abs(out.cn)), '< 1e-12', Math.abs(out.cn) < 1e-12);
    add(g, 'Euler法による |W| の増加が |1 − ifΔt|^(2n) − 1 と一致する',
      e3(Math.abs(out.euler - predicted) / predicted), '< 1e-3',
      Math.abs(out.euler - predicted) / predicted < 1e-3);
  }

  // -- 9. Ekman pumping ------------------------------------------------------
  {
    const g = G('緯度依存の風が駆動するエクマン湧昇（層の底に生じる鉛直流）',
      ['南北格子数 Ny', 'w_e の相対誤差', '次数']);
    const base = { ...BASE, H: 8 * DE, N: 64 };
    const Ly = 4.0e6, tau0 = 0.1;
    let prev = null, lastOrder = null, ampl = 0, net = 0;
    for (const Ny of [32, 64, 128, 256]) {
      const r = solveUpwelling(base, { Ly, Ny, tau0 });
      let err = 0; ampl = 0;
      for (let j = 0; j < Ny; j++) {
        err = Math.max(err, Math.abs(r.we[j] - r.weExact[j]));
        ampl = Math.max(ampl, Math.abs(r.weExact[j]));
      }
      const rel = err / ampl;
      const o = order(prev, rel);
      g.rows.push([Ny, e3(rel), ord(o)]);
      prev = rel; lastOrder = o; net = r.net;
    }
    g.notes.push(
      `鉛直流速の最大値は ${(ampl * 86400 * 365).toFixed(1)} m/年で、実際に観測される` +
      'エクマン湧昇と同じオーダーである。南北に周期境界を課しているので湧昇と沈降は' +
      '厳密に相殺しなければならず、正味の流量は ' + net.toExponential(3) + ' m²/s になっている。');
    g.notes.push(
      '「w_e が curl(τ/f)/ρ₀ に収束する」の次数は、南北の格子数 Ny を倍にすると湧昇速度の' +
      '誤差が 1/4 になることを表す。ここで測っているのは鉛直方向ではなく、輸送の南北微分に' +
      '使っている中心差分の精度である。');
    add(g, 'w_e が curl(τ/f)/ρ₀ に収束する', lastOrder.toFixed(4), '2.00 ± 0.05', Math.abs(lastOrder - 2) < 0.05);
    add(g, '湧昇と沈降が相殺する', e3(Math.abs(net / (ampl * Ly))), '< 1e-14',
      Math.abs(net / (ampl * Ly)) < 1e-14);
  }

  return {
    config: {
      f: BASE.f, nu: BASE.nu, rho0: BASE.rho0, taux: BASE.taux, tauy: BASE.tauy,
      bottom: BASE.bottom,
      DE, inertialHours: (2 * Math.PI) / Math.abs(BASE.f) / 3600,
    },
    groups,
    checks,
    passed: checks.filter((c) => c.pass).length,
    total: checks.length,
  };
}
