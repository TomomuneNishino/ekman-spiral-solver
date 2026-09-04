#!/usr/bin/env node
// Verification suite for the Ekman layer solver.
//
//   node tests/verify.mjs           formatted report
//   node tests/verify.mjs --json    machine-readable results
//
// All of the work happens in src/checks.js, which the browser front end also
// imports, so the numbers below describe the deployed solver rather than a
// separate throwaway script.

import { runChecks, DEFAULT_BASE } from '../src/checks.js';

const WIDTH = 84;

// Terminal columns, not code points: CJK characters occupy two cells, so
// padding by string length would misalign every table containing Japanese.
const isWide = (c) => {
  const p = c.codePointAt(0);
  return (p >= 0x1100 && p <= 0x115f) || (p >= 0x2e80 && p <= 0xa4cf)
    || (p >= 0xac00 && p <= 0xd7a3) || (p >= 0xf900 && p <= 0xfaff)
    || (p >= 0xfe30 && p <= 0xfe6f) || (p >= 0xff00 && p <= 0xff60)
    || (p >= 0xffe0 && p <= 0xffe6) || (p >= 0x20000 && p <= 0x3fffd);
};
const cols = (s) => [...String(s)].reduce((n, c) => n + (isWide(c) ? 2 : 1), 0);
const padStart = (s, n) => ' '.repeat(Math.max(0, n - cols(s))) + s;

// Japanese has no inter-word spaces, so wrapping has to count columns and
// avoid leaving closing punctuation stranded at the start of a line.
const NO_LINE_START = '。、）」』】〕》,.):;!?';
function wrap(text, width) {
  const out = [];
  let line = '';
  for (const ch of [...String(text)]) {
    if (cols(line) + cols(ch) > width && !NO_LINE_START.includes(ch)) {
      out.push(line);
      line = ch === ' ' ? '' : ch;
    } else {
      line += ch;
    }
  }
  if (line.trim()) out.push(line);
  return out;
}

const r = runChecks(DEFAULT_BASE);
const rule = (ch = '=') => console.log(ch.repeat(WIDTH));

rule();
console.log('海洋エクマン層シミュレーター — 検証スイート');
rule();
console.log(`  緯度 北緯 45°         f      = ${r.config.f.toExponential(4)} 1/s`);
console.log(`  鉛直渦粘性係数        nu_z   = ${r.config.nu.toExponential(4)} m^2/s`);
console.log(`  基準密度              rho0   = ${r.config.rho0} kg/m^3`);
console.log(`  風応力                tau    = (${r.config.taux}, ${r.config.tauy}) N/m^2  東向き`);
console.log(`  エクマン深さ          D_E    = ${r.config.DE.toFixed(3)} m`);
console.log(`  慣性周期              2 pi/f = ${r.config.inertialHours.toFixed(3)} h`);
console.log('');

for (const [i, g] of r.groups.entries()) {
  const head = `-- ${i + 1}. ${g.title} `;
  console.log(head + '-'.repeat(Math.max(0, WIDTH - cols(head))));
  if (g.columns) {
    const all = [g.columns, ...g.rows.map((row) => row.map(String))];
    const w = g.columns.map((_, c) => Math.max(...all.map((row) => cols(row[c]))));
    const fmt = (row) => '   ' + row.map((v, c) => padStart(String(v), w[c] + 2)).join('');
    console.log(fmt(g.columns));
    console.log('   ' + '-'.repeat(w.reduce((a, b) => a + b + 2, 0)));
    for (const row of g.rows) console.log(fmt(row.map(String)));
  }
  for (const n of g.notes) {
    console.log('');
    for (const line of wrap(n, WIDTH - 6)) console.log('   ' + line);
  }
  console.log('');
  for (const c of g.checks) {
    console.log(`   [${c.pass ? '合格' : '不合格'}] ${c.id}. ${c.name}`);
    console.log(`          結果 ${c.value}   期待 ${c.expected}`);
  }
  console.log('');
}

rule();
console.log(`  ${r.total} 項目中 ${r.passed} 項目が合格`);
rule();

if (process.argv.includes('--json')) console.log(JSON.stringify(r, null, 2));
process.exit(r.passed === r.total ? 0 : 1);
