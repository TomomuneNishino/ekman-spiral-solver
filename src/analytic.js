// Closed-form solutions used as the reference for verification.
//
// Steady state:  nu W'' = i f W,  so W ~ exp(lambda z) with
//     lambda = sqrt(i f / nu) = (1 + i) / D_E,   D_E = sqrt(2 nu / f)
//
// Three reference solutions are provided, and it matters which one is used for
// which test:
//   * semiInfiniteProfile -- the textbook Ekman spiral on z in (-inf, 0]
//   * steadyProfile       -- the exact solution of the *same* boundary value
//                            problem the code solves, on the finite domain
//                            [-H, 0]. Grid-convergence tests must use this one,
//                            otherwise the domain-truncation error masquerades
//                            as a discretization error.
//   * spinupTransport     -- the exact time-dependent depth-integrated transport

import { C, add, sub, mul, div, scale, exp, abs, arg, conj } from './complex.js';
import { ekmanDepth } from './solver.js';

/**
 * lambda = sqrt(i f / nu) = (1 + i sign(f)) / D_E.
 *
 * The sign of f carries the hemisphere: in the north (f > 0) the spiral turns
 * clockwise with depth and the surface current sits 45 degrees to the right of
 * the wind; in the south every rotation is mirrored.
 */
export function lambda(p) {
  const dE = ekmanDepth(p);
  return C(1 / dE, Math.sign(p.f) / dE);
}

const stress = (p) => C(p.taux, p.tauy);

/**
 * Textbook Ekman spiral on a semi-infinite ocean:
 *   W(z) = T e^{-i pi/4} / (rho0 sqrt(f nu)) * e^{lambda z}
 */
export function semiInfiniteProfile(p, z) {
  const lam = lambda(p);
  const amp = surfaceVelocity(p);
  return z.map((zz) => mul(amp, exp(scale(lam, zz))));
}

/**
 * Exact steady solution on the finite domain [-H, 0].
 *
 *   free BC:    W = T cosh(lambda (z+H)) / (rho0 nu lambda sinh(lambda H))
 *   noslip BC:  W = T sinh(lambda (z+H)) / (rho0 nu lambda cosh(lambda H))
 *
 * Both are rewritten in terms of exp(lambda z) and exp(-2 lambda H) so that
 * nothing overflows when H >> D_E.
 */
export function steadyProfile(p, z) {
  const lam = lambda(p);
  const T = stress(p);
  const pref = div(T, scale(mul(C(p.nu, 0), lam), p.rho0)); // T / (rho0 nu lambda)
  const e2H = exp(scale(lam, -2 * p.H));
  const sign = p.bottom === 'noslip' ? -1 : 1;
  const den = p.bottom === 'noslip' ? add(C(1, 0), e2H) : sub(C(1, 0), e2H);

  return z.map((zz) => {
    const ez = exp(scale(lam, zz));                    // e^{lambda z}
    const eh = exp(scale(lam, -(zz + 2 * p.H)));       // e^{-lambda (z + 2H)}
    const num = add(ez, scale(eh, sign));
    return mul(pref, div(num, den));
  });
}

/**
 * Exact depth-integrated Ekman transport rho0 * int W dz on the finite domain.
 *
 * For the free (zero-stress) bottom this is exactly (tau_y - i tau_x)/f for
 * *any* H -- the truncation of the domain does not contaminate it at all,
 * because the only stresses crossing the boundaries are T at the top and zero
 * at the bottom. For the no-slip bottom there is a residual bottom stress and
 * the transport is short by a factor sech(lambda H) ~ 2 exp(-H/D_E).
 */
export function transportExact(p) {
  const T = stress(p);
  const ideal = scale(mul(C(0, -1), T), 1 / p.f); // -i T / f
  if (p.bottom !== 'noslip') return ideal;
  const lam = lambda(p);
  const e = exp(scale(lam, -p.H));
  const sech = div(scale(e, 2), add(C(1, 0), exp(scale(lam, -2 * p.H))));
  return mul(ideal, sub(C(1, 0), sech));
}

/** The classical result: (tau_y - i tau_x) / f, i.e. 90 degrees right of the wind. */
export function transportClassical(p) {
  return C(p.tauy / p.f, -p.taux / p.f);
}

/**
 * Exact depth-integrated transport during spin-up from rest.
 *
 * Integrating the momentum equation over the layer removes the viscous term
 * entirely (zero stress at the bottom), leaving dM/dt = T - i f M, hence
 *
 *     M(t) = -i (T/f) (1 - e^{-i f t})
 *
 * This is independent of nu, of H and of the grid, so it tests the time
 * integration in isolation. Physically it is an inertial oscillation spiralling
 * onto the steady Ekman transport.
 */
export function spinupTransport(p, t) {
  const T = stress(p);
  const decay = sub(C(1, 0), exp(C(0, -p.f * t)));
  return scale(mul(mul(C(0, -1), T), decay), 1 / p.f);
}

/**
 * Steady rate of work done by the wind on the Ekman layer, per unit area.
 *
 *   P = Re(conj(T) W(0)) = |T|^2 cos(pi/4) / (rho0 sqrt(f nu)) = |T|^2 / (rho0 f D_E)
 *
 * The cos(pi/4) is exactly the 45-degree surface deflection, so this number and
 * the 45-degree result are two faces of the same fact. In steady state the wind
 * input is balanced by viscous dissipation, P = epsilon.
 */
export function steadyPower(p) {
  const T2 = p.taux * p.taux + p.tauy * p.tauy;
  return T2 / (p.rho0 * Math.abs(p.f) * ekmanDepth(p));
}

/**
 * Surface velocity of the semi-infinite spiral,
 *   W(0) = T exp(-i sign(f) pi/4) / (rho0 sqrt(|f| nu)),
 * i.e. 45 degrees to the right of the wind in the northern hemisphere and 45
 * degrees to the left in the southern.
 */
export function surfaceVelocity(p) {
  const s = Math.sign(p.f);
  return scale(
    mul(stress(p), C(Math.SQRT1_2, -s * Math.SQRT1_2)),
    1 / (p.rho0 * Math.sqrt(Math.abs(p.f) * p.nu)),
  );
}

/** Deflection of the surface current from the wind stress, in degrees. */
export function deflectionDeg(W0, p) {
  let d = ((arg(W0) - arg(stress(p))) * 180) / Math.PI;
  while (d > 180) d -= 360;
  while (d < -180) d += 360;
  return d;
}

export { abs, arg, conj, mul, div, add, sub, scale, C };
