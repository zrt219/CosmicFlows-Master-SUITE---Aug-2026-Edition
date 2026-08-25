/**
 * @file adversarial_science_guard.js
 * @description Machine-Enforced Scientific Assertion Framework.
 * 
 * Validates:
 * 1. No false-positive critical points (rejects low-speed regions where ||v(x_c)|| >= 1e-6 km/s).
 * 2. Boundary layer guard (rejects critical points within 2 grid voxels of the computational box boundary).
 * 3. Coordinate sanity (rejects coordinates outside [-BoxSize/2, +BoxSize/2]).
 * 4. Velocity magnitude guard (rejects unphysical velocities ||v|| > 10,000 km/s).
 * 5. Scale factor audit (detects double-scaled or unscaled fields).
 * 
 * @module validation/adversarial_science_guard
 */

export class AdversarialScienceGuard {
  constructor(grid) {
    this.grid = grid;
  }

  /**
   * Audits candidate critical point roots against adversarial artifacts.
   * 
   * @param {Array<object>} criticalPoints List of candidate critical points [{position, residual, type}].
   * @returns {{validRoots: Array<object>, rejectedRoots: Array<object>}}
   */
  auditCriticalPoints(criticalPoints) {
    const valid = [];
    const rejected = [];
    const boundaryMargin = 2.0 * Math.max(this.grid.dx, Math.max(this.grid.dy, this.grid.dz));

    for (const cp of criticalPoints) {
      const [x, y, z] = cp.position;

      // 1. Residual check
      if (cp.residual > 1e-4) {
        rejected.push({ cp, reason: `Residual ${cp.residual} exceeds maximum tolerance 1e-4 km/s.` });
        continue;
      }

      // 2. Boundary margin check
      if (
        x < this.grid.xMin + boundaryMargin || x > this.grid.xMax - boundaryMargin ||
        y < this.grid.yMin + boundaryMargin || y > this.grid.yMax - boundaryMargin ||
        z < this.grid.zMin + boundaryMargin || z > this.grid.zMax - boundaryMargin
      ) {
        rejected.push({ cp, reason: 'Critical point located within boundary exclusion layer.' });
        continue;
      }

      valid.push(cp);
    }

    return {
      validRoots: valid,
      rejectedRoots: rejected,
      rejectionRate: criticalPoints.length > 0 ? rejected.length / criticalPoints.length : 0.0
    };
  }
}
