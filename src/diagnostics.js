// Discrete diagnostics and budget residuals.
//
// The point of this module is that every quantity below is derived from the
// *discrete* equations, not from the continuous ones. The residuals it returns
// should therefore vanish to round-off, not merely to truncation error. Any
// residual that grows with resolution or with time is a bug, not a modelling
// approximation, which makes these functions the primary correctness check.
//
// Summing the cell momentum equations telescopes the interior face stresses and
// leaves
//     dM/dt = T - tau_bottom - i f M,          M = rho0 dz sum_k W_k
//
// Taking Re(conj(W_k) * (cell equation)) and summing gives
//     dE/dt = Re(conj(T) W_0) - epsilon,       E = (rho0/2) dz sum_k |W_k|^2
//     epsilon = rho0 nu / dz * sum_faces |W_{j-1} - W_j|^2  >= 0
// with the Coriolis term contributing exactly zero because Re(-i f |W|^2) = 0.

/** Depth-integrated momentum M = rho0 dz sum W_k  [kg/(m s)]. */
export function transport(W, p) {
  const dz = p.H / p.N;
  let re = 0, im = 0;
  for (let k = 0; k < p.N; k++) { re += W.re[k]; im += W.im[k]; }
  return { re: p.rho0 * dz * re, im: p.rho0 * dz * im };
}

/** Kinetic energy per unit area E = (rho0/2) dz sum |W_k|^2  [J/m^2]. */
export function energy(W, p) {
  const dz = p.H / p.N;
  let s = 0;
  for (let k = 0; k < p.N; k++) s += W.re[k] * W.re[k] + W.im[k] * W.im[k];
  return 0.5 * p.rho0 * dz * s;
}

/**
 * Stress on the bottom face.
 * free:   zero by construction.
 * noslip: rho0 nu * 2 W_{N-1} / dz, using the half-cell distance to the wall.
 */
export function bottomStress(W, p) {
  if (p.bottom !== 'noslip') return { re: 0, im: 0 };
  const dz = p.H / p.N;
  const c = (2 * p.rho0 * p.nu) / dz;
  return { re: c * W.re[p.N - 1], im: c * W.im[p.N - 1] };
}

/**
 * Rate of work done by the wind, in the form that makes the discrete energy
 * budget close exactly: Re(conj(T) * W_0), with W_0 the top *cell centre*.
 */
export function windPowerDiscrete(W, p) {
  return p.taux * W.re[0] + p.tauy * W.im[0];
}

/**
 * Second-order estimate of the true surface velocity, W(0) ~ W_0 + (dz/2) W'(0),
 * using the exactly known surface gradient W'(0) = T/(rho0 nu). Used when
 * comparing against the analytic surface current and the 45-degree result.
 */
export function surfaceVelocityExtrapolated(W, p) {
  const dz = p.H / p.N;
  const c = dz / (2 * p.rho0 * p.nu);
  return { re: W.re[0] + c * p.taux, im: W.im[0] + c * p.tauy };
}

/** Rate of work using the extrapolated surface velocity. */
export function windPowerPhysical(W, p) {
  const s = surfaceVelocityExtrapolated(W, p);
  return p.taux * s.re + p.tauy * s.im;
}

/**
 * Viscous dissipation, written as a sum of squares over the interior faces so
 * that it is manifestly non-negative. The no-slip wall adds one more term.
 */
export function dissipation(W, p) {
  const dz = p.H / p.N;
  const c = (p.rho0 * p.nu) / dz;
  let s = 0;
  for (let j = 1; j <= p.N - 1; j++) {
    const dr = W.re[j - 1] - W.re[j];
    const di = W.im[j - 1] - W.im[j];
    s += dr * dr + di * di;
  }
  let eps = c * s;
  if (p.bottom === 'noslip') {
    eps += 2 * c * (W.re[p.N - 1] ** 2 + W.im[p.N - 1] ** 2);
  }
  return eps;
}

/**
 * Work done by the Coriolis force, computed as written rather than assumed.
 * Should be exactly zero in floating point, term by term.
 */
export function coriolisWork(W, p) {
  const dz = p.H / p.N;
  let s = 0;
  for (let k = 0; k < p.N; k++) {
    // u * (f v) + v * (-f u)
    s += W.re[k] * (p.f * W.im[k]) + W.im[k] * (-p.f * W.re[k]);
  }
  return p.rho0 * dz * s;
}

/**
 * Momentum budget residual over one step.
 *   R = (M^{n+1} - M^n)/dt - (T - tau_bottom - i f M)
 * evaluated at the time-centred state Wmid returned by Integrator.step().
 */
export function momentumResidual(Wold, Wnew, Wmid, p, dt) {
  const Mo = transport(Wold, p);
  const Mn = transport(Wnew, p);
  const Mm = transport(Wmid, p);
  const tb = bottomStress(Wmid, p);
  const lhsR = (Mn.re - Mo.re) / dt;
  const lhsI = (Mn.im - Mo.im) / dt;
  // -i f M  ->  real: f*M.im, imag: -f*M.re
  const rhsR = p.taux - tb.re + p.f * Mm.im;
  const rhsI = p.tauy - tb.im - p.f * Mm.re;
  const scaleM = Math.hypot(p.taux, p.tauy) || 1;
  return {
    re: lhsR - rhsR,
    im: lhsI - rhsI,
    relative: Math.hypot(lhsR - rhsR, lhsI - rhsI) / scaleM,
  };
}

/**
 * Energy budget residual over one step.
 *   R = (E^{n+1} - E^n)/dt - (P - epsilon)
 * again evaluated at the time-centred state.
 */
export function energyResidual(Wold, Wnew, Wmid, p, dt) {
  const dE = (energy(Wnew, p) - energy(Wold, p)) / dt;
  const P = windPowerDiscrete(Wmid, p);
  const eps = dissipation(Wmid, p);
  const scaleE = Math.abs(P) + Math.abs(eps) || 1;
  return { absolute: dE - (P - eps), relative: (dE - (P - eps)) / scaleE, P, eps, dE };
}

/**
 * Dissipation with the surface half-cell included.
 *
 * The face sum above runs over interior faces only, so it covers the strip
 * z in [-H + dz/2, -dz/2] and misses the half cell just below the sea surface.
 * That strip has a known gradient, W'(0) = T/(rho0 nu), and contributes
 *
 *     rho0 nu (dz/2) |W'(0)|^2 = (dz/2) |T|^2 / (rho0 nu)
 *
 * which is exactly the amount by which the discrete wind power Re(conj(T) W_0)
 * -- evaluated at the top cell centre rather than at the surface -- undercounts
 * the physical wind power. The two corrections are the same number, so adding
 * it to both sides preserves the exact discrete balance while making each side
 * a second-order approximation to its continuous counterpart.
 */
export function dissipationPhysical(W, p) {
  const dz = p.H / p.N;
  const T2 = p.taux * p.taux + p.tauy * p.tauy;
  return dissipation(W, p) + (dz / 2) * (T2 / (p.rho0 * p.nu));
}

/** Convenience bundle for the UI. */
export function diagnose(W, p) {
  return {
    M: transport(W, p),
    E: energy(W, p),
    P: windPowerPhysical(W, p),
    eps: dissipation(W, p),
    surface: surfaceVelocityExtrapolated(W, p),
    coriolisWork: coriolisWork(W, p),
  };
}
