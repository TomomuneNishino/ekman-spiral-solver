// Minimal complex-number helpers.
//
// Two representations are used in this project:
//   * `{re, im}` objects  -- readable, used for closed-form analytic expressions
//   * paired Float64Array -- fast, used inside the tridiagonal solver
//
// Keeping the analytic formulas in object form makes them easy to check against
// the equations in the write-up; the hot loop never allocates.

export const C = (re, im = 0) => ({ re, im });

export const add = (a, b) => ({ re: a.re + b.re, im: a.im + b.im });
export const sub = (a, b) => ({ re: a.re - b.re, im: a.im - b.im });
export const mul = (a, b) => ({
  re: a.re * b.re - a.im * b.im,
  im: a.re * b.im + a.im * b.re,
});

export function div(a, b) {
  // Smith's algorithm: avoids overflow when |b| is very large or very small.
  if (Math.abs(b.re) >= Math.abs(b.im)) {
    const r = b.im / b.re;
    const d = b.re + b.im * r;
    return { re: (a.re + a.im * r) / d, im: (a.im - a.re * r) / d };
  }
  const r = b.re / b.im;
  const d = b.re * r + b.im;
  return { re: (a.re * r + a.im) / d, im: (a.im * r - a.re) / d };
}

export const scale = (a, s) => ({ re: a.re * s, im: a.im * s });
export const conj = (a) => ({ re: a.re, im: -a.im });
export const abs = (a) => Math.hypot(a.re, a.im);
export const arg = (a) => Math.atan2(a.im, a.re);

export function exp(a) {
  const m = Math.exp(a.re);
  return { re: m * Math.cos(a.im), im: m * Math.sin(a.im) };
}

// Real part of conj(a) * b -- the discrete inner product used for work and energy.
export const dot = (a, b) => a.re * b.re + a.im * b.im;

// i * a
export const imul = (a) => ({ re: -a.im, im: a.re });
