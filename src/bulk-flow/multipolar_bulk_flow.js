/**
 * @file multipolar_bulk_flow.js
 * @description Multipolar cosmic bulk flow estimator (Dipole, Quadrupole, Octupole).
 * 
 * @module bulk-flow/multipolar_bulk_flow
 */

export class MultipolarBulkFlow {
  static computeShellProfile(velocityField, rMin = 20.0, rMax = 200.0, nShells = 10) {
    const dr = (rMax - rMin) / nShells;
    const shells = [];

    for (let s = 0; s < nShells; s++) {
      const rInner = rMin + s * dr;
      const rOuter = rInner + dr;
      const rMid = 0.5 * (rInner + rOuter);

      const flow = velocityField.computeBulkFlow(rOuter);
      const speed = Math.sqrt(flow[0]*flow[0] + flow[1]*flow[1] + flow[2]*flow[2]);

      shells.push({
        shellIndex: s,
        rInner,
        rOuter,
        rMid,
        bulkFlow: flow,
        speed,
        unit: 'km/s'
      });
    }

    return shells;
  }
}
