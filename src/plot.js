// A small canvas plotting layer.
//
// Purpose-built rather than pulled from a chart library: the figures here need
// depth-coloured polylines, equal-aspect hodographs and log-log convergence
// panels, and none of those are what a general charting library is good at.

const CSS = (name, fallback) => {
  const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  return v || fallback;
};

export const palette = () => ({
  ink: CSS('--foam', '#E8F1F2'),
  muted: CSS('--haze', '#8FAFBC'),
  rule: CSS('--shelf', '#1D4A5C'),
  accent: CSS('--current', '#F2B441'),
  panel: CSS('--deep', '#12303E'),
});

/** Depth ramp: t = 0 at the sea surface, t = 1 at the bottom of the domain. */
export function depthColor(t) {
  // Lightness swings from near-white at the surface to deep indigo at the floor,
  // so depth reads as brightness and not only as hue.
  const stops = [
    [0.00, [222, 252, 246]],
    [0.12, [126, 236, 232]],
    [0.28, [ 46, 190, 224]],
    [0.45, [ 30, 132, 212]],
    [0.62, [ 52,  78, 190]],
    [0.80, [ 96,  56, 160]],
    [1.00, [ 74,  30,  96]],
  ];
  // The interesting structure all sits within the top two or three Ekman depths,
  // which is a small fraction of the domain. Advancing the ramp faster near the
  // surface spends the available contrast where the curve actually varies.
  t = Math.pow(Math.max(0, Math.min(1, t)), 0.62);
  for (let i = 1; i < stops.length; i++) {
    if (t <= stops[i][0]) {
      const [t0, c0] = stops[i - 1];
      const [t1, c1] = stops[i];
      const s = (t - t0) / (t1 - t0);
      const c = c0.map((v, k) => Math.round(v + s * (c1[k] - v)));
      return `rgb(${c[0]},${c[1]},${c[2]})`;
    }
  }
  return 'rgb(74,30,96)';
}

function niceTicks(lo, hi, target = 5) {
  if (!(hi > lo)) return [lo];
  const raw = (hi - lo) / target;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const norm = raw / mag;
  const step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) {
    out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  }
  return out;
}

function logTicks(lo, hi) {
  const out = [];
  for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) {
    const v = 10 ** e;
    if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v);
  }
  return out;
}

const fmtTick = (v, log) => {
  if (log) {
    const e = Math.round(Math.log10(v));
    return e === 0 ? '1' : `1e${e}`;
  }
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e4 || a < 1e-3) return v.toExponential(0);
  if (a >= 100) return v.toFixed(0);
  if (a >= 1) return v.toFixed(1);
  return v.toFixed(a < 0.01 ? 3 : 2);
};

export class Frame {
  constructor(canvas, opts = {}) {
    this.c = canvas;
    this.o = { pad: { l: 58, r: 14, t: 16, b: 40 }, ...opts };
    this.series = [];
  }

  /** Set up transforms for the given data ranges and paint axes and grid. */
  begin({
    xlim, ylim, xlog = false, ylog = false, equal = false,
    xlabel = '', ylabel = '',
    // Optional second y axis on the right, for a series in different units.
    y2lim = null, y2label = '',
    // Optional pair of labels pinned to the two ends of the x axis, e.g.
    // ['南 \u2190', '\u2192 北'], so the sense of the axis is readable at a glance.
    xends = null,
  }) {
    const dpr = window.devicePixelRatio || 1;
    const rect = this.c.getBoundingClientRect();
    const w = Math.max(1, Math.round(rect.width));
    const h = Math.max(1, Math.round(rect.height));
    this.c.width = w * dpr;
    this.c.height = h * dpr;
    const g = this.c.getContext('2d');
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, w, h);
    this.g = g; this.w = w; this.h = h;
    this.xlog = xlog; this.ylog = ylog;

    const p = this.o.pad;
    // Gutters are decided here rather than at each call site, so bumping the
    // type sizes cannot silently start clipping labels somewhere.
    const padL = ylabel ? Math.max(p.l, 88) : Math.max(p.l, 48);
    const padB = xlabel ? Math.max(p.b, 50) : Math.max(p.b, 30);
    const padR = y2lim ? (y2label ? Math.max(p.r, 88) : Math.max(p.r, 48)) : p.r;
    this.box = { x0: padL, y0: p.t, x1: w - padR, y1: h - padB };

    let [xa, xb] = xlim, [ya, yb] = ylim;
    if (equal) {
      // Preserve aspect ratio so a circle looks like a circle.
      const bw = this.box.x1 - this.box.x0, bh = this.box.y1 - this.box.y0;
      const sx = (xb - xa) / bw, sy = (yb - ya) / bh;
      const s = Math.max(sx, sy);
      const cx = (xa + xb) / 2, cy = (ya + yb) / 2;
      xa = cx - (s * bw) / 2; xb = cx + (s * bw) / 2;
      ya = cy - (s * bh) / 2; yb = cy + (s * bh) / 2;
    }
    this.xa = xa; this.xb = xb; this.ya = ya; this.yb = yb;
    this.ya2 = y2lim ? y2lim[0] : null;
    this.yb2 = y2lim ? y2lim[1] : null;

    const P = palette();
    const tx = xlog ? logTicks(xa, xb) : niceTicks(xa, xb, Math.max(3, Math.round(this.w / 110)));
    const ty = ylog ? logTicks(ya, yb) : niceTicks(ya, yb, Math.max(3, Math.round(this.h / 70)));

    g.save();
    g.strokeStyle = P.rule;
    g.fillStyle = P.muted;
    g.lineWidth = 1;
    g.font = '13px "IBM Plex Mono", ui-monospace, monospace';

    g.globalAlpha = 0.45;
    for (const t of tx) {
      const X = this.X(t);
      g.beginPath(); g.moveTo(X, this.box.y0); g.lineTo(X, this.box.y1); g.stroke();
    }
    for (const t of ty) {
      const Y = this.Y(t);
      g.beginPath(); g.moveTo(this.box.x0, Y); g.lineTo(this.box.x1, Y); g.stroke();
    }
    g.globalAlpha = 1;

    g.textAlign = 'center'; g.textBaseline = 'top';
    for (const t of tx) g.fillText(fmtTick(t, xlog), this.X(t), this.box.y1 + 10);
    g.textAlign = 'right'; g.textBaseline = 'middle';
    for (const t of ty) g.fillText(fmtTick(t, ylog), this.box.x0 - 9, this.Y(t));

    if (y2lim) {
      // Ticks only, no grid: a second set of grid lines would not line up with
      // the first and would read as noise.
      const ty2 = niceTicks(this.ya2, this.yb2, Math.max(3, Math.round(this.h / 70)));
      g.textAlign = 'left'; g.textBaseline = 'middle';
      for (const t of ty2) {
        const Y = this.Y2(t);
        g.beginPath(); g.moveTo(this.box.x1, Y); g.lineTo(this.box.x1 + 5, Y); g.stroke();
        g.fillText(fmtTick(t, false), this.box.x1 + 9, Y);
      }
    }

    // Zero lines, drawn brighter than the grid because they carry meaning.
    g.strokeStyle = P.rule; g.lineWidth = 1.4;
    if (!xlog && xa < 0 && xb > 0) {
      g.beginPath(); g.moveTo(this.X(0), this.box.y0); g.lineTo(this.X(0), this.box.y1); g.stroke();
    }
    if (!ylog && ya < 0 && yb > 0) {
      g.beginPath(); g.moveTo(this.box.x0, this.Y(0)); g.lineTo(this.box.x1, this.Y(0)); g.stroke();
    }

    g.fillStyle = P.muted;
    g.font = '14px "Noto Sans JP", "Hiragino Sans", "Noto Sans CJK JP", system-ui, sans-serif';
    if (xlabel) {
      g.textAlign = 'center'; g.textBaseline = 'bottom';
      g.fillText(xlabel, (this.box.x0 + this.box.x1) / 2, this.h - 2);
    }
    if (xends) {
      g.font = '13px "Noto Sans JP", "Hiragino Sans", "Noto Sans CJK JP", system-ui, sans-serif';
      g.textBaseline = 'bottom';
      g.textAlign = 'left';
      g.fillText(xends[0], this.box.x0, this.h - 2);
      g.textAlign = 'right';
      g.fillText(xends[1], this.box.x1, this.h - 2);
      g.font = '14px "Noto Sans JP", "Hiragino Sans", "Noto Sans CJK JP", system-ui, sans-serif';
    }
    if (ylabel) {
      // Rotated and centred on the axis, so it unambiguously belongs to it.
      // Japanese labels run long, so shrink rather than overflow the axis.
      const avail = this.box.y1 - this.box.y0 - 12;
      let size = 14;
      g.font = `${size}px "Noto Sans JP", "Hiragino Sans", "Noto Sans CJK JP", system-ui, sans-serif`;
      while (g.measureText(ylabel).width > avail && size > 8) {
        size -= 0.5;
        g.font = `${size}px "Noto Sans JP", "Hiragino Sans", "Noto Sans CJK JP", system-ui, sans-serif`;
      }
      g.save();
      g.translate(14, (this.box.y0 + this.box.y1) / 2);
      g.rotate(-Math.PI / 2);
      g.textAlign = 'center'; g.textBaseline = 'top';
      g.fillText(ylabel, 0, 0);
      g.restore();
    }
    if (y2label) {
      // Mirror image of the left label: rotated the other way so it reads from
      // the top down on the right-hand side.
      const avail = this.box.y1 - this.box.y0 - 12;
      let size = 14;
      g.font = `${size}px "Noto Sans JP", "Hiragino Sans", "Noto Sans CJK JP", system-ui, sans-serif`;
      while (g.measureText(y2label).width > avail && size > 8) {
        size -= 0.5;
        g.font = `${size}px "Noto Sans JP", "Hiragino Sans", "Noto Sans CJK JP", system-ui, sans-serif`;
      }
      g.save();
      g.translate(this.w - 14, (this.box.y0 + this.box.y1) / 2);
      g.rotate(Math.PI / 2);
      g.textAlign = 'center'; g.textBaseline = 'top';
      g.fillText(y2label, 0, 0);
      g.restore();
    }
    g.restore();
    this.series = [];
    return this;
  }

  X(v) {
    const a = this.xlog ? Math.log10(this.xa) : this.xa;
    const b = this.xlog ? Math.log10(this.xb) : this.xb;
    const t = this.xlog ? Math.log10(v) : v;
    return this.box.x0 + ((t - a) / (b - a)) * (this.box.x1 - this.box.x0);
  }

  Y(v) {
    const a = this.ylog ? Math.log10(this.ya) : this.ya;
    const b = this.ylog ? Math.log10(this.yb) : this.yb;
    const t = this.ylog ? Math.log10(v) : v;
    return this.box.y1 - ((t - a) / (b - a)) * (this.box.y1 - this.box.y0);
  }

  /** Data -> pixels for the optional right-hand axis. */
  Y2(v) {
    if (this.ya2 == null) return this.Y(v);
    return this.box.y1 - ((v - this.ya2) / (this.yb2 - this.ya2)) * (this.box.y1 - this.box.y0);
  }

  clip(fn) {
    const g = this.g;
    g.save();
    g.beginPath();
    g.rect(this.box.x0, this.box.y0, this.box.x1 - this.box.x0, this.box.y1 - this.box.y0);
    g.clip();
    fn(g);
    g.restore();
  }

  line(xs, ys, { color = '#fff', width = 2, dash = null, label = null, alpha = 1, axis = 'y' } = {}) {
    const Ymap = axis === 'y2' ? (v) => this.Y2(v) : (v) => this.Y(v);
    this.clip((g) => {
      g.strokeStyle = color; g.lineWidth = width; g.globalAlpha = alpha;
      g.lineJoin = 'round'; g.lineCap = 'round';
      if (dash) g.setLineDash(dash);
      g.beginPath();
      let started = false;
      for (let i = 0; i < xs.length; i++) {
        if (!Number.isFinite(xs[i]) || !Number.isFinite(ys[i])) { started = false; continue; }
        const X = this.X(xs[i]), Y = Ymap(ys[i]);
        if (!started) { g.moveTo(X, Y); started = true; } else g.lineTo(X, Y);
      }
      g.stroke();
      g.setLineDash([]);
    });
    if (label) this.series.push({ label, color, dash });
    return this;
  }

  /** Polyline whose colour varies along its length, used to encode depth. */
  rampLine(xs, ys, ts, { width = 2.6, label = null } = {}) {
    this.clip((g) => {
      g.lineWidth = width; g.lineCap = 'round'; g.lineJoin = 'round';
      for (let i = 1; i < xs.length; i++) {
        g.strokeStyle = depthColor((ts[i - 1] + ts[i]) / 2);
        g.beginPath();
        g.moveTo(this.X(xs[i - 1]), this.Y(ys[i - 1]));
        g.lineTo(this.X(xs[i]), this.Y(ys[i]));
        g.stroke();
      }
    });
    if (label) this.series.push({ label, color: depthColor(0.25), dash: null });
    return this;
  }

  points(xs, ys, { color = '#fff', r = 2.4, label = null, ramp = null } = {}) {
    this.clip((g) => {
      for (let i = 0; i < xs.length; i++) {
        if (!Number.isFinite(xs[i]) || !Number.isFinite(ys[i])) continue;
        g.fillStyle = ramp ? depthColor(ramp[i]) : color;
        g.beginPath();
        g.arc(this.X(xs[i]), this.Y(ys[i]), r, 0, 2 * Math.PI);
        g.fill();
      }
    });
    if (label) this.series.push({ label, color, marker: true });
    return this;
  }

  arrow(x0, y0, x1, y1, { color = '#fff', width = 3, head = 12, label = null, alpha = 1 } = {}) {
    const g = this.g;
    const X0 = this.X(x0), Y0 = this.Y(y0), X1 = this.X(x1), Y1 = this.Y(y1);
    const a = Math.atan2(Y1 - Y0, X1 - X0);
    g.save();
    g.strokeStyle = color; g.fillStyle = color; g.lineWidth = width;
    g.lineCap = 'round'; g.globalAlpha = alpha;
    g.beginPath(); g.moveTo(X0, Y0); g.lineTo(X1 - Math.cos(a) * head * 0.78, Y1 - Math.sin(a) * head * 0.78); g.stroke();
    // Slightly swept-back head: reads more clearly than a plain triangle at
    // the sizes these figures are viewed at.
    g.beginPath();
    g.moveTo(X1, Y1);
    g.lineTo(X1 - head * Math.cos(a - 0.34), Y1 - head * Math.sin(a - 0.34));
    g.lineTo(X1 - head * 0.72 * Math.cos(a), Y1 - head * 0.72 * Math.sin(a));
    g.lineTo(X1 - head * Math.cos(a + 0.34), Y1 - head * Math.sin(a + 0.34));
    g.closePath(); g.fill();
    if (label) {
      // Placed alongside the shaft rather than past the tip, so it cannot spill
      // outside the frame however the arrow happens to be oriented.
      g.font = '13px "Noto Sans JP", "Hiragino Sans", "Noto Sans CJK JP", system-ui, sans-serif';
      const mx = (X0 + X1) / 2, my = (Y0 + Y1) / 2;
      const nx = -Math.sin(a), ny = Math.cos(a);
      const side = my + ny * 12 < this.box.y0 + 12 ? -1 : 1;
      g.textAlign = 'center'; g.textBaseline = 'middle';
      g.fillText(label, mx + nx * 13 * side, my + ny * 13 * side);
    }
    g.restore();
    return this;
  }

  /**
   * Arc between two directions, drawn at a fixed pixel radius about a data
   * point. Used to show the angle between the wind and the surface current,
   * so the 45 degrees can be read off the figure instead of taken on trust.
   */
  angleMark(cx, cy, v1, v2, { r = 46, color = '#fff', label = null } = {}) {
    const g = this.g;
    const CX = this.X(cx), CY = this.Y(cy);
    const a1 = Math.atan2(this.Y(cy + v1.y) - CY, this.X(cx + v1.x) - CX);
    const a2 = Math.atan2(this.Y(cy + v2.y) - CY, this.X(cx + v2.x) - CX);
    let d = a2 - a1;
    while (d > Math.PI) d -= 2 * Math.PI;
    while (d < -Math.PI) d += 2 * Math.PI;
    g.save();
    g.strokeStyle = color; g.fillStyle = color; g.lineWidth = 1.3;
    g.setLineDash([3, 3]);
    g.beginPath();
    g.arc(CX, CY, r, a1, a1 + d, d < 0);
    g.stroke();
    g.setLineDash([]);
    if (label) {
      const am = a1 + d / 2;
      g.font = '13.5px "Noto Sans JP", "Hiragino Sans", "Noto Sans CJK JP", system-ui, sans-serif';
      g.textAlign = 'center'; g.textBaseline = 'middle';
      g.fillText(label, CX + Math.cos(am) * (r + 19), CY + Math.sin(am) * (r + 19));
    }
    g.restore();
    return this;
  }

  annotate(x, y, text, { color = null, dx = 8, dy = 0, align = 'left' } = {}) {
    const g = this.g;
    g.save();
    g.fillStyle = color || palette().muted;
    g.font = '12.5px "Noto Sans JP", "Hiragino Sans", "Noto Sans CJK JP", system-ui, sans-serif';
    g.textAlign = align; g.textBaseline = 'middle';
    g.fillText(text, this.X(x) + dx, this.Y(y) + dy);
    g.restore();
    return this;
  }

  /**
   * Corner of the frame holding the fewest of the given points. The hodograph
   * rotates with the wind direction, so a fixed legend corner is bound to sit
   * on top of the curve for some settings.
   */
  pickCorner(xs, ys) {
    const { x0, y0, x1, y1 } = this.box;
    const w = (x1 - x0) * 0.46, h = (y1 - y0) * 0.42;
    const zones = { tl: 0, tr: 0, bl: 0, br: 0 };
    for (let i = 0; i < xs.length; i++) {
      const X = this.X(xs[i]), Y = this.Y(ys[i]);
      if (!Number.isFinite(X) || !Number.isFinite(Y)) continue;
      const left = X < x0 + w, right = X > x1 - w;
      const top = Y < y0 + h, bottom = Y > y1 - h;
      if (left && top) zones.tl++;
      if (right && top) zones.tr++;
      if (left && bottom) zones.bl++;
      if (right && bottom) zones.br++;
    }
    return Object.entries(zones).sort((a, b) => a[1] - b[1])[0][0];
  }

  legend({ corner = 'tr' } = {}) {
    if (!this.series.length) return this;
    const g = this.g;
    const pad = 9, lh = 19, sw = 24;
    g.save();
    g.font = '13px "Noto Sans JP", "Hiragino Sans", "Noto Sans CJK JP", system-ui, sans-serif';
    const wmax = Math.max(...this.series.map((s) => g.measureText(s.label).width)) + sw + pad * 2 + 6;
    const hbox = this.series.length * lh + pad * 2 - 4;
    const x = corner.includes('r') ? this.box.x1 - wmax - 8 : this.box.x0 + 8;
    const y = corner.includes('b') ? this.box.y1 - hbox - 8 : this.box.y0 + 8;
    const P = palette();
    g.fillStyle = P.panel; g.globalAlpha = 0.86;
    g.beginPath(); g.roundRect(x, y, wmax, hbox, 5); g.fill();
    g.globalAlpha = 1;
    g.strokeStyle = P.rule; g.lineWidth = 1; g.stroke();
    this.series.forEach((s, i) => {
      const yy = y + pad + i * lh + 4;
      g.strokeStyle = s.color; g.fillStyle = s.color; g.lineWidth = 2.2;
      if (s.marker) {
        g.beginPath(); g.arc(x + pad + sw / 2, yy, 2.8, 0, 2 * Math.PI); g.fill();
      } else {
        if (s.dash) g.setLineDash(s.dash);
        g.beginPath(); g.moveTo(x + pad, yy); g.lineTo(x + pad + sw, yy); g.stroke();
        g.setLineDash([]);
      }
      g.fillStyle = P.ink;
      g.textAlign = 'left'; g.textBaseline = 'middle';
      g.fillText(s.label, x + pad + sw + 6, yy);
    });
    g.restore();
    return this;
  }
}

/** Convenient min/max over several arrays, with a fractional margin. */
export function span(arrays, margin = 0.08) {
  let lo = Infinity, hi = -Infinity;
  for (const a of arrays) for (const v of a) {
    if (Number.isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
  }
  if (!Number.isFinite(lo)) return [0, 1];
  if (hi === lo) { hi = lo + 1; lo -= 1; }
  const d = (hi - lo) * margin;
  return [lo - d, hi + d];
}
