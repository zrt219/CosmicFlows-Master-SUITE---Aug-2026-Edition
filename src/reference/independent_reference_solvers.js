/**
 * @file independent_reference_solvers.js
 * @description Independent analytical and simple numerical reference implementations
 * for cross-verifying optimized production engines (Phases 101, 102).
 * 
 * Includes:
 * - Reference direct trilinear interpolator
 * - Reference classical fixed-step RK4 integrator
 * - Reference central finite difference gradient
 * - Reference direct determinant and Newton root evaluator
 * 
 * @module reference/independent_reference_solvers
 */

export class ReferenceTrilinearInterpolator {
  static eval(grid, buffer, x, y, z) {
    const [gx, gy, gz] = grid.coordToGridIndex(x, y, z);
    const ix = Math.max(0, Math.min(grid.nx - 2, Math.floor(gx)));
    const iy = Math.max(0, Math.min(grid.ny - 2, Math.floor(gy)));
    const iz = Math.max(0, Math.min(grid.nz - 2, Math.floor(gz)));

    const tx = gx - ix;
    const ty = gy - iy;
    const tz = gz - iz;

    const c000 = buffer[grid.index(ix, iy, iz)];
    const c100 = buffer[grid.index(ix + 1, iy, iz)];
    const c010 = buffer[grid.index(ix, iy + 1, iz)];
    const c110 = buffer[grid.index(ix + 1, iy + 1, iz)];
    const c001 = buffer[grid.index(ix, iy, iz + 1)];
    const c101 = buffer[grid.index(ix + 1, iy, iz + 1)];
    const c011 = buffer[grid.index(ix, iy + 1, iz + 1)];
    const c111 = buffer[grid.index(ix + 1, iy + 1, iz + 1)];

    const c00 = c000 * (1 - tx) + c100 * tx;
    const c10 = c010 * (1 - tx) + c110 * tx;
    const c01 = c001 * (1 - tx) + c101 * tx;
    const c11 = c011 * (1 - tx) + c111 * tx;

    const c0 = c00 * (1 - ty) + c10 * ty;
    const c1 = c01 * (1 - ty) + c11 * ty;

    return c0 * (1 - tz) + c1 * tz;
  }
}

export class ReferenceRK4Integrator {
  static step(velocityFn, x, y, z, dt) {
    const k1 = velocityFn(x, y, z);
    const k2 = velocityFn(x + 0.5 * dt * k1[0], y + 0.5 * dt * k1[1], z + 0.5 * dt * k1[2]);
    const k3 = velocityFn(x + 0.5 * dt * k2[0], y + 0.5 * dt * k2[1], z + 0.5 * dt * k2[2]);
    const k4 = velocityFn(x + dt * k3[0], y + dt * k3[1], z + dt * k3[2]);

    const xNext = x + (dt / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    const yNext = y + (dt / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    const zNext = z + (dt / 6.0) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);

    return [xNext, yNext, zNext];
  }
}
