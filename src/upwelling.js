// Ekman pumping and suction.
//
// The one-dimensional spiral above assumes a horizontally uniform wind. Let the
// zonal stress instead vary with latitude, as the real zonal-wind belts do:
//
//     tau_x(y) = tau0 cos(2 pi y / Ly),   tau_y = 0
//
// Nothing in the momentum equations couples neighbouring columns -- there is no
// horizontal advection and no horizontal viscosity -- so every column is still
// an independent one-dimensional problem and can be solved with exactly the
// same routine. What changes is the diagnostic: the depth-integrated transport
// now has a horizontal divergence, and continuity converts that divergence into
// a vertical velocity at the base of the Ekman layer,
//
//     w_e = d(U_E)/dx + d(V_E)/dy
//
// with U_E + i V_E = (1/(rho0 f)) (tau_y - i tau_x) so that, for f constant,
//
//     w_e = -(1/(rho0 f)) d(tau_x)/dy = (tau0 / (rho0 f)) (2 pi / Ly) sin(2 pi y / Ly)
//
// Positive w_e is upwelling. The domain is periodic in y, which makes the
// closed-loop constraint that upwelling and downwelling balance exactly,
// integral of w_e dy = 0, a genuine test rather than an approximation.

import { solveSteady } from './solver.js';
import { transport } from './diagnostics.js';

/** Uniformly spaced periodic y grid, cell centred. */
export function yGrid(Ly, Ny) {
  const y = new Float64Array(Ny);
  const dy = Ly / Ny;
  for (let j = 0; j < Ny; j++) y[j] = (j + 0.5) * dy;
  return y;
}

/** tau_x(y) = tau0 cos(2 pi y / Ly), tau_y = 0. */
export function stressProfile(y, Ly, tau0) {
  return Array.from(y, (yy) => ({
    taux: tau0 * Math.cos((2 * Math.PI * yy) / Ly),
    tauy: 0,
  }));
}

/**
 * Solve one steady column per latitude and diagnose the Ekman pumping.
 * Returns transports (per unit density, i.e. m^2/s) and w_e in m/s.
 */
export function solveUpwelling(base, { Ly, Ny, tau0 }) {
  const y = yGrid(Ly, Ny);
  const dy = Ly / Ny;
  const stresses = stressProfile(y, Ly, tau0);

  const U = new Float64Array(Ny);
  const V = new Float64Array(Ny);
  const columns = [];

  for (let j = 0; j < Ny; j++) {
    const p = { ...base, taux: stresses[j].taux, tauy: stresses[j].tauy };
    const W = solveSteady(p);
    const M = transport(W, p);
    U[j] = M.re / p.rho0; // volume transport per unit width [m^2/s]
    V[j] = M.im / p.rho0;
    columns.push(W);
  }

  // Centred difference on the periodic grid; d/dx vanishes because nothing
  // depends on x.
  const we = new Float64Array(Ny);
  for (let j = 0; j < Ny; j++) {
    const jp = (j + 1) % Ny;
    const jm = (j - 1 + Ny) % Ny;
    we[j] = (V[jp] - V[jm]) / (2 * dy);
  }

  // Analytic pumping, w_e = (1/rho0) { d/dx (tau_y/f) - d/dy (tau_x/f) }.
  const weExact = new Float64Array(Ny);
  for (let j = 0; j < Ny; j++) {
    weExact[j] =
      ((tau0 * (2 * Math.PI)) / (base.rho0 * base.f * Ly)) *
      Math.sin((2 * Math.PI * y[j]) / Ly);
  }

  // Net vertical volume flux through the base of the layer, per unit length in x.
  let net = 0;
  for (let j = 0; j < Ny; j++) net += we[j] * dy;

  return { y, dy, stresses, U, V, we, weExact, net, columns };
}
