/**
 * @file streamline_event_detector.js
 * @description High-precision integration event detector for streamline tracing.
 * 
 * Detects:
 * - Domain boundary exit (hyper-rectangular bounding box crossing)
 * - Critical point root neighborhood capture (||x - x_root|| < r_capture)
 * - Low-velocity stall (||v|| < v_min)
 * - Maximum arc-length exhaustion (s >= s_max)
 * - Orbital recurrence / limit cycles (inf_{t' < t - T_min} ||x(t) - x(t')|| < delta_recurrence)
 * 
 * @module streamlines/streamline_event_detector
 */

export const StreamlineEventType = Object.freeze({
  NONE: 'NONE',
  DOMAIN_EXIT: 'DOMAIN_EXIT',
  ROOT_CAPTURED: 'ROOT_CAPTURED',
  LOW_SPEED_STALL: 'LOW_SPEED_STALL',
  MAX_ARCLENGTH_EXHAUSTED: 'MAX_ARCLENGTH_EXHAUSTED',
  ORBITAL_RECURRENCE: 'ORBITAL_RECURRENCE',
  MAX_STEPS_REACHED: 'MAX_STEPS_REACHED'
});

export class StreamlineEventDetector {
  /**
   * @param {object} grid GridIndexer defining bounding box.
   * @param {object} [options]
   */
  constructor(grid, options = {}) {
    this.grid = grid;
    this.rCapture = options.rCapture || 2.0; // Mpc/h
    this.vMin = options.vMin || 5.0;        // km/s
    this.maxArcLength = options.maxArcLength || 500.0; // Mpc/h
    this.recurrenceWindow = options.recurrenceWindow || 30; // steps
    this.recurrenceTol = options.recurrenceTol || 0.5; // Mpc/h
  }

  /**
   * Checks for integration events at the current step.
   * 
   * @param {Array<number>} currentPos [x, y, z] in Mpc/h.
   * @param {Array<number>} currentVel [vx, vy, vz] in km/s.
   * @param {number} currentArcLength Total cumulative path length.
   * @param {Array<Array<number>>} pathHistory Array of past positions.
   * @param {Array<object>} [knownRoots=[]] Array of known critical point roots [{position: [x,y,z], type}].
   * @returns {{type: string, detail: object}} Detected event.
   */
  checkStep(currentPos, currentVel, currentArcLength, pathHistory, knownRoots = []) {
    const [x, y, z] = currentPos;

    // 1. Domain boundary exit
    if (
      x < this.grid.xMin || x > this.grid.xMax ||
      y < this.grid.yMin || y > this.grid.yMax ||
      z < this.grid.zMin || z > this.grid.zMax
    ) {
      return {
        type: StreamlineEventType.DOMAIN_EXIT,
        detail: { position: currentPos, bounds: [this.grid.xMin, this.grid.xMax] }
      };
    }

    // 2. Velocity stall
    const speed = Math.sqrt(currentVel[0] ** 2 + currentVel[1] ** 2 + currentVel[2] ** 2);
    if (speed < this.vMin) {
      return {
        type: StreamlineEventType.LOW_SPEED_STALL,
        detail: { speed, threshold: this.vMin }
      };
    }

    // 3. Known root capture
    for (const root of knownRoots) {
      const rx = root.position[0], ry = root.position[1], rz = root.position[2];
      const dist = Math.sqrt((x - rx) ** 2 + (y - ry) ** 2 + (z - rz) ** 2);
      if (dist < this.rCapture) {
        return {
          type: StreamlineEventType.ROOT_CAPTURED,
          detail: { root, distance: dist }
        };
      }
    }

    // 4. Max arc length
    if (currentArcLength >= this.maxArcLength) {
      return {
        type: StreamlineEventType.MAX_ARCLENGTH_EXHAUSTED,
        detail: { arcLength: currentArcLength, max: this.maxArcLength }
      };
    }

    // 5. Orbital loop / recurrence detection
    const n = pathHistory.length;
    if (n > this.recurrenceWindow) {
      for (let i = 0; i < n - this.recurrenceWindow; i++) {
        const px = pathHistory[i][0], py = pathHistory[i][1], pz = pathHistory[i][2];
        const d = Math.sqrt((x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2);
        if (d < this.recurrenceTol) {
          return {
            type: StreamlineEventType.ORBITAL_RECURRENCE,
            detail: { matchedIndex: i, distance: d }
          };
        }
      }
    }

    return { type: StreamlineEventType.NONE, detail: null };
  }
}
