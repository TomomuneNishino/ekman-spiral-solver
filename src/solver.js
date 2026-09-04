// Numerical solver for the wind-driven oceanic Ekman layer.
//
// Governing equation, written in flux (conservative) form for the horizontal
// momentum per unit volume rho0*u, rho0*v:
//
//     d(rho0 u)/dt = d(tau_xz)/dz + rho0 f v
//     d(rho0 v)/dt = d(tau_yz)/dz - rho0 f u
//     tau_xz = rho0 nu dz(u),   tau_yz = rho0 nu dz(v)
//
// With the complex velocity W = u + i v and the complex stress T = tau_x + i tau_y
// this collapses to a single scalar equation:
//
//     dW/dt = nu * d2W/dz2 - i f W
//
// Boundary conditions
//   surface  z = 0 : rho0 nu dz(W) = T          (wind stress, imposed as a face flux)
//   bottom  z = -H : dz(W) = 0  (free)  or  W = 0  (noslip)
//
// GRID (staggered finite volume)
//   Cell centres  z_k     = -(k + 1/2) dz,  k = 0 .. N-1   (k = 0 is the top cell)
//   Cell faces    z_{j}   = -j dz,          j = 0 .. N      (j = 0 is the sea surface)
//   Velocities live at cell centres, stresses live at faces.
//
// Imposing the wind stress directly as the value of the j = 0 face flux (rather
// than through a ghost cell) is what makes the discrete momentum and energy
// budgets close to machine precision: summing the cell equations telescopes and
// leaves only the two boundary fluxes.

/**
 * @typedef {Object} Params
 * @property {number} f      Coriolis parameter [1/s]
 * @property {number} nu     vertical eddy viscosity [m^2/s]
 * @property {number} rho0   reference density [kg/m^3]
 * @property {number} taux   zonal wind stress [N/m^2]
 * @property {number} tauy   meridional wind stress [N/m^2]
 * @property {number} H      domain depth [m]
 * @property {number} N      number of cells
 * @property {'free'|'noslip'} bottom  bottom boundary condition
 */

export const OMEGA_EARTH = 7.2921e-5; // Earth's rotation rate [1/s]

/** Coriolis parameter at a given latitude in degrees. */
export const coriolis = (latitudeDeg) =>
  2 * OMEGA_EARTH * Math.sin((latitudeDeg * Math.PI) / 180);

/** Ekman depth D_E = sqrt(2 nu / f). */
export const ekmanDepth = (p) => Math.sqrt((2 * p.nu) / Math.abs(p.f));

/** Cell-centre depths, top cell first. */
export function cellCentres(p) {
  const dz = p.H / p.N;
  const z = new Float64Array(p.N);
  for (let k = 0; k < p.N; k++) z[k] = -(k + 0.5) * dz;
  return z;
}

/**
 * Tridiagonal complex operator L such that dW/dt = L W + S.
 * Row k reads  a_k W_{k-1} + b_k W_k + c_k W_{k+1}.
 * S is nonzero only in the top cell, where it carries the wind stress.
 */
export function buildOperator(p) {
  const { N, nu, f, H, rho0, taux, tauy, bottom } = p;
  const dz = H / N;
  const r = nu / (dz * dz);

  const aR = new Float64Array(N), aI = new Float64Array(N);
  const bR = new Float64Array(N), bI = new Float64Array(N);
  const cR = new Float64Array(N), cI = new Float64Array(N);

  for (let k = 0; k < N; k++) {
    aR[k] = k === 0 ? 0 : r;
    cR[k] = k === N - 1 ? 0 : r;
    // Interior cells see two active faces; the top cell's upper face carries the
    // prescribed wind stress and so contributes nothing to the operator.
    let diag = -2 * r;
    if (k === 0) diag = -r;
    if (k === N - 1) diag = bottom === 'noslip' ? -3 * r : -r;
    bR[k] = diag;
    bI[k] = -f; // Coriolis: pure rotation, no contribution to the real part
  }

  // Source from the surface stress: (tau_top - tau_bottom)/(rho0 dz) with tau_top = T.
  const sR = new Float64Array(N), sI = new Float64Array(N);
  sR[0] = taux / (rho0 * dz);
  sI[0] = tauy / (rho0 * dz);

  return { aR, aI, bR, bI, cR, cI, sR, sI, dz, N };
}

/**
 * Complex Thomas algorithm. Solves the tridiagonal system in place-free fashion
 * and returns {re, im}. Arrays a, b, c, d all have length n.
 */
export function thomas(aR, aI, bR, bI, cR, cI, dR, dI, n) {
  const cpR = new Float64Array(n), cpI = new Float64Array(n);
  const dpR = new Float64Array(n), dpI = new Float64Array(n);

  const cdiv = (xr, xi, yr, yi) => {
    if (Math.abs(yr) >= Math.abs(yi)) {
      const t = yi / yr, den = yr + yi * t;
      return [(xr + xi * t) / den, (xi - xr * t) / den];
    }
    const t = yr / yi, den = yr * t + yi;
    return [(xr * t + xi) / den, (xi * t - xr) / den];
  };

  let [q0r, q0i] = cdiv(cR[0], cI[0], bR[0], bI[0]);
  cpR[0] = q0r; cpI[0] = q0i;
  let [p0r, p0i] = cdiv(dR[0], dI[0], bR[0], bI[0]);
  dpR[0] = p0r; dpI[0] = p0i;

  for (let k = 1; k < n; k++) {
    // m = b_k - a_k * cp_{k-1}
    const mr = bR[k] - (aR[k] * cpR[k - 1] - aI[k] * cpI[k - 1]);
    const mi = bI[k] - (aR[k] * cpI[k - 1] + aI[k] * cpR[k - 1]);
    const [qr, qi] = cdiv(cR[k], cI[k], mr, mi);
    cpR[k] = qr; cpI[k] = qi;
    const nr = dR[k] - (aR[k] * dpR[k - 1] - aI[k] * dpI[k - 1]);
    const ni = dI[k] - (aR[k] * dpI[k - 1] + aI[k] * dpR[k - 1]);
    const [pr, pi] = cdiv(nr, ni, mr, mi);
    dpR[k] = pr; dpI[k] = pi;
  }

  const xR = new Float64Array(n), xI = new Float64Array(n);
  xR[n - 1] = dpR[n - 1]; xI[n - 1] = dpI[n - 1];
  for (let k = n - 2; k >= 0; k--) {
    xR[k] = dpR[k] - (cpR[k] * xR[k + 1] - cpI[k] * xI[k + 1]);
    xI[k] = dpI[k] - (cpR[k] * xI[k + 1] + cpI[k] * xR[k + 1]);
  }
  return { re: xR, im: xI };
}

/**
 * Steady state of the discrete system, obtained by solving L W = -S directly.
 * This isolates the spatial truncation error from any time-stepping error.
 */
export function solveSteady(p) {
  const op = buildOperator(p);
  const { N } = op;
  const dR = new Float64Array(N), dI = new Float64Array(N);
  for (let k = 0; k < N; k++) { dR[k] = -op.sR[k]; dI[k] = -op.sI[k]; }
  return thomas(op.aR, op.aI, op.bR, op.bI, op.cR, op.cI, dR, dI, N);
}

/**
 * Time integrator.
 *
 * scheme = 'cn'    Crank-Nicolson. Unconditionally stable, second order, and --
 *                  crucially for this problem -- the Cayley transform it applies
 *                  to the Coriolis term has amplification factor exactly 1, so
 *                  rotation neither creates nor destroys energy.
 * scheme = 'euler' Forward Euler / FTCS. Included only for comparison: its
 *                  Coriolis amplification factor is |1 - i f dt| > 1, so it
 *                  spuriously injects energy at O(f^2 dt^2).
 */
export class Integrator {
  constructor(p, { scheme = 'cn', dt = null } = {}) {
    this.p = { ...p };
    this.scheme = scheme;
    this.op = buildOperator(p);
    this.dt = dt ?? this.suggestedDt();
    this.t = 0;
    this.nstep = 0;
    this.W = { re: new Float64Array(p.N), im: new Float64Array(p.N) }; // start from rest
    if (scheme === 'cn') this.#prepareCN();
  }

  /** dt suggestion: resolve the inertial period; respect diffusive stability for Euler. */
  suggestedDt() {
    const inertial = (2 * Math.PI) / Math.abs(this.p.f);
    const dzs = (this.p.H / this.p.N) ** 2 / (2 * this.p.nu);
    return this.scheme === 'cn'
      ? inertial / 200
      : 0.4 * Math.min(dzs, inertial / 200);
  }

  #prepareCN() {
    const { N } = this.op;
    const h = this.dt / 2;
    // Left-hand side I - (dt/2) L, stored once and reused every step.
    this.lhs = {
      aR: new Float64Array(N), aI: new Float64Array(N),
      bR: new Float64Array(N), bI: new Float64Array(N),
      cR: new Float64Array(N), cI: new Float64Array(N),
    };
    for (let k = 0; k < N; k++) {
      this.lhs.aR[k] = -h * this.op.aR[k];
      this.lhs.aI[k] = -h * this.op.aI[k];
      this.lhs.bR[k] = 1 - h * this.op.bR[k];
      this.lhs.bI[k] = -h * this.op.bI[k];
      this.lhs.cR[k] = -h * this.op.cR[k];
      this.lhs.cI[k] = -h * this.op.cI[k];
    }
  }

  /** Apply L to a state, returning L W (without the source term). */
  applyL(W) {
    const { N, aR, aI, bR, bI, cR, cI } = this.op;
    const outR = new Float64Array(N), outI = new Float64Array(N);
    for (let k = 0; k < N; k++) {
      let re = bR[k] * W.re[k] - bI[k] * W.im[k];
      let im = bR[k] * W.im[k] + bI[k] * W.re[k];
      if (k > 0) {
        re += aR[k] * W.re[k - 1] - aI[k] * W.im[k - 1];
        im += aR[k] * W.im[k - 1] + aI[k] * W.re[k - 1];
      }
      if (k < N - 1) {
        re += cR[k] * W.re[k + 1] - cI[k] * W.im[k + 1];
        im += cR[k] * W.im[k + 1] + cI[k] * W.re[k + 1];
      }
      outR[k] = re; outI[k] = im;
    }
    return { re: outR, im: outI };
  }

  /**
   * Advance one step. Returns the time-centred state used by the step, which is
   * the state at which the discrete budgets close exactly:
   *   Crank-Nicolson -> (W^n + W^{n+1}) / 2
   *   forward Euler  -> W^n
   */
  step() {
    const { N, sR, sI } = this.op;
    const Wold = { re: this.W.re.slice(), im: this.W.im.slice() };

    if (this.scheme === 'cn') {
      const h = this.dt / 2;
      const LW = this.applyL(this.W);
      const rR = new Float64Array(N), rI = new Float64Array(N);
      for (let k = 0; k < N; k++) {
        rR[k] = this.W.re[k] + h * LW.re[k] + this.dt * sR[k];
        rI[k] = this.W.im[k] + h * LW.im[k] + this.dt * sI[k];
      }
      this.W = thomas(
        this.lhs.aR, this.lhs.aI, this.lhs.bR, this.lhs.bI,
        this.lhs.cR, this.lhs.cI, rR, rI, N,
      );
    } else {
      const LW = this.applyL(this.W);
      const nR = new Float64Array(N), nI = new Float64Array(N);
      for (let k = 0; k < N; k++) {
        nR[k] = this.W.re[k] + this.dt * (LW.re[k] + sR[k]);
        nI[k] = this.W.im[k] + this.dt * (LW.im[k] + sI[k]);
      }
      this.W = { re: nR, im: nI };
    }

    this.t += this.dt;
    this.nstep += 1;

    if (this.scheme === 'cn') {
      const mR = new Float64Array(N), mI = new Float64Array(N);
      for (let k = 0; k < N; k++) {
        mR[k] = 0.5 * (Wold.re[k] + this.W.re[k]);
        mI[k] = 0.5 * (Wold.im[k] + this.W.im[k]);
      }
      return { Wold, Wmid: { re: mR, im: mI }, Wnew: this.W };
    }
    return { Wold, Wmid: Wold, Wnew: this.W };
  }

  /** Advance until at least `tEnd`, returning the number of steps taken. */
  runTo(tEnd) {
    let n = 0;
    while (this.t < tEnd - 1e-12) { this.step(); n++; }
    return n;
  }
}
