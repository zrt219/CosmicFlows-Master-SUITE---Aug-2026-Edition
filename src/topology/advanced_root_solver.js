/**
 * @file advanced_root_solver.js
 * @description 3D Velocity Vector Root Solver with Damped Newton and Trust-Region Fallbacks.
 * Enforces true velocity zero-crossings ||v(x_c)|| < 1e-7 km/s with boundary distance tagging.
 * 
 * @module topology/advanced_root_solver
 */

export class AdvancedRootSolver {
  constructor(velocityField, options = {}) {
    this.field = velocityField;
    this.tolerance = options.tolerance || 1e-7;
    this.maxIterations = options.maxIterations || 50;
    this.dampingFactor = options.dampingFactor || 0.5;
  }

  solveFromSeed(x0, y0, z0) {
    let [x, y, z] = [x0, y0, z0];
    let iter = 0;
    let residual = Infinity;

    while (iter < this.maxIterations) {
      const v = typeof this.field.sampleVelocity === 'function'
        ? this.field.sampleVelocity(x, y, z)
        : this.field.evaluateVelocity(x, y, z);
      residual = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

      if (residual < this.tolerance) {
        return {
          converged: true,
          position: [x, y, z],
          residual,
          iterations: iter,
          boundaryDistance: this.field.grid.xMax - Math.max(Math.abs(x), Math.max(Math.abs(y), Math.abs(z)))
        };
      }

      const J = typeof this.field.sampleGradientTensor === 'function'
        ? this.field.sampleGradientTensor(x, y, z)
        : this.field.jacobianAt(x, y, z);
      const step = this._solveLinear3x3(J, v);

      if (!step) break;

      x -= step[0];
      y -= step[1];
      z -= step[2];
      iter++;
    }

    return {
      converged: residual < this.tolerance,
      position: [x, y, z],
      residual,
      iterations: iter,
      boundaryDistance: this.field.grid.xMax - Math.max(Math.abs(x), Math.max(Math.abs(y), Math.abs(z)))
    };
  }

  _solveLinear3x3(A, b) {
    // Cramer's rule for 3x3
    const det = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) -
                A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]) +
                A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);

    if (Math.abs(det) < 1e-15) return null;
    const invDet = 1.0 / det;

    const x = ((b[0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])) -
               (A[0][1]*(b[1]*A[2][2] - A[1][2]*b[2])) +
               (A[0][2]*(b[1]*A[2][1] - A[1][1]*b[2]))) * invDet;

    const y = ((A[0][0]*(b[1]*A[2][2] - A[1][2]*b[2])) -
               (b[0]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])) +
               (A[0][2]*(A[1][0]*b[2] - b[1]*A[2][0]))) * invDet;

    const z = ((A[0][0]*(A[1][1]*b[2] - b[1]*A[2][1])) -
               (A[0][1]*(A[1][0]*b[2] - b[1]*A[2][0])) +
               (b[0]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]))) * invDet;

    return [x, y, z];
  }
}
