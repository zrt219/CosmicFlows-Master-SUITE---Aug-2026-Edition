/**
 * @file supergalactic_frame.js
 * @description Astrometric modeling of the de Vaucouleurs Supergalactic coordinate frame,
 * landmark supercluster coordinates, Local Sheet structure, and cosmic web landmarks
 * for the CosmicFlows-4 Research Workbench.
 * 
 * Reference:
 * - de Vaucouleurs et al. (1976) Second Reference Catalogue of Bright Galaxies (RC2)
 * - Tully et al. (2014) Nature 513, 71-73 (Laniakea supercluster)
 * - Dupuy & Courtois (2023) A&A 678, A176 (CosmicFlows-4 watershed boundaries)
 * 
 * @module coordinates/supergalactic_frame
 */

import { ASTROMETRIC_CONSTANTS } from './canonical_frame.js';
import { SupergalacticPosition } from './scientific_types.js';

/**
 * Authoritative Supergalactic Landmarks and Major Cosmological Structures.
 * Positions in Supergalactic Cartesian coordinates (SGX, SGY, SGZ) in units of Mpc/h.
 * @readonly
 */
export const SUPERGALACTIC_LANDMARKS = Object.freeze({
  VIRGO_CLUSTER: Object.freeze({
    name: 'Virgo Cluster Core (M87 / NGC 4486)',
    sgx: -3.4,
    sgy: 15.9,
    sgz: -0.5,
    distanceMpcOverH: 16.3,
    cz_kms: 1284.0,
    massMsunOverH: 1.2e14,
    description: 'Central gravitational core of the Local Supercluster'
  }),
  GREAT_ATTRACTOR: Object.freeze({
    name: 'Norma Cluster / Great Attractor Core (ACO 3627)',
    sgx: -35.0,
    sgy: 15.0,
    sgz: -10.0,
    distanceMpcOverH: 48.0,
    cz_kms: 4871.0,
    massMsunOverH: 5.0e15,
    description: 'Central gravitational attractor of the Laniakea basin of attraction'
  }),
  CENTAURUS_CLUSTER: Object.freeze({
    name: 'Centaurus Cluster (ACO 3526)',
    sgx: -25.0,
    sgy: 22.0,
    sgz: -5.0,
    distanceMpcOverH: 33.7,
    cz_kms: 3041.0,
    massMsunOverH: 2.0e14,
    description: 'Major cluster in the Centaurus-Hydra wall'
  }),
  HYDRA_CLUSTER: Object.freeze({
    name: 'Hydra Cluster (ACO 1060)',
    sgx: -30.0,
    sgy: 5.0,
    sgz: -15.0,
    distanceMpcOverH: 33.9,
    cz_kms: 3777.0,
    massMsunOverH: 2.5e14,
    description: 'Southern component of Hydra-Centaurus supercluster'
  }),
  PERSEUS_PISCES: Object.freeze({
    name: 'Perseus Cluster (ACO 426) / Perseus-Pisces Spine',
    sgx: 52.0,
    sgy: -20.0,
    sgz: -12.0,
    distanceMpcOverH: 56.5,
    cz_kms: 5366.0,
    massMsunOverH: 1.0e15,
    description: 'Prominent filamentary supercluster spine opposing Laniakea'
  }),
  COMA_CLUSTER: Object.freeze({
    name: 'Coma Cluster (ACO 1656)',
    sgx: 0.0,
    sgy: 74.0,
    sgz: 9.0,
    distanceMpcOverH: 74.5,
    cz_kms: 6925.0,
    massMsunOverH: 1.5e15,
    description: 'Great Wall anchor and rich galaxy cluster'
  }),
  SHAPLEY_CONCENTRATION: Object.freeze({
    name: 'Shapley Supercluster Core (ACO 3558)',
    sgx: -120.0,
    sgy: 65.0,
    sgz: -25.0,
    distanceMpcOverH: 139.0,
    cz_kms: 14500.0,
    massMsunOverH: 1.0e16,
    description: 'Most massive gravitationally bound concentration in the local volume'
  }),
  HERCULES_SUPERCLUSTER: Object.freeze({
    name: 'Hercules Cluster Complex (ACO 2151/2152)',
    sgx: 18.0,
    sgy: 46.0,
    sgz: 32.0,
    distanceMpcOverH: 58.7,
    cz_kms: 10800.0,
    massMsunOverH: 8.0e14,
    description: 'Hercules basin attractor'
  }),
  FORNAX_CLUSTER: Object.freeze({
    name: 'Fornax Cluster (NGC 1399)',
    sgx: -14.0,
    sgy: -12.0,
    sgz: 0.5,
    distanceMpcOverH: 18.5,
    cz_kms: 1425.0,
    massMsunOverH: 7.0e13,
    description: 'Nearby low-mass galaxy cluster in the southern sky'
  }),
  PAVO_INDUS: Object.freeze({
    name: 'Pavo-Indus Supercluster Core',
    sgx: -45.0,
    sgy: -30.0,
    sgz: -8.0,
    distanceMpcOverH: 54.7,
    cz_kms: 4200.0,
    massMsunOverH: 5.0e14,
    description: 'Southern filamentary wall connecting to Centaurus'
  })
});

/**
 * SupergalacticFrame provides methods for analyzing positions relative to the Supergalactic Plane,
 * Local Sheet, and known cosmological structures.
 */
export class SupergalacticFrame {
  /**
   * Evaluates the perpendicular distance of a point to the Supergalactic Plane (Z = 0).
   * @param {SupergalacticPosition|Array<number>} pos
   * @returns {number} Distance in Mpc/h (|SGZ|)
   */
  static distanceToSupergalacticPlane(pos) {
    const z = pos.z ?? pos[2];
    return Math.abs(z);
  }

  /**
   * Determines whether a position lies inside the Local Sheet (Tully et al. 2008).
   * The Local Sheet is a planar structure of thickness ~1.5 Mpc/h containing the Milky Way, Andromeda, and Virgo filament.
   * @param {SupergalacticPosition|Array<number>} pos
   * @param {number} [sheetHalfThickness=1.5]
   * @returns {boolean}
   */
  static isInLocalSheet(pos, sheetHalfThickness = 1.5) {
    const x = pos.x ?? pos[0];
    const y = pos.y ?? pos[1];
    const z = pos.z ?? pos[2];
    const inPlaneDistance = Math.sqrt(x * x + y * y);
    return Math.abs(z) <= sheetHalfThickness && inPlaneDistance <= 10.0;
  }

  /**
   * Finds the nearest major astronomical landmark to a given position.
   * @param {SupergalacticPosition|Array<number>} pos
   * @returns {{ landmark: Object, key: string, distanceMpcOverH: number }}
   */
  static findNearestLandmark(pos) {
    const x = pos.x ?? pos[0];
    const y = pos.y ?? pos[1];
    const z = pos.z ?? pos[2];

    let closestKey = null;
    let closestDist = Infinity;
    let closestLandmark = null;

    for (const [key, landmark] of Object.entries(SUPERGALACTIC_LANDMARKS)) {
      const dx = x - landmark.sgx;
      const dy = y - landmark.sgy;
      const dz = z - landmark.sgz;
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (dist < closestDist) {
        closestDist = dist;
        closestKey = key;
        closestLandmark = landmark;
      }
    }

    return {
      key: closestKey,
      landmark: closestLandmark,
      distanceMpcOverH: closestDist
    };
  }

  /**
   * Creates a SupergalacticPosition instance for a named landmark.
   * @param {string} landmarkKey
   * @returns {SupergalacticPosition}
   */
  static getLandmarkPosition(landmarkKey) {
    const landmark = SUPERGALACTIC_LANDMARKS[landmarkKey];
    if (!landmark) {
      throw new Error(`SupergalacticFrame: Unknown landmark key '${landmarkKey}'.`);
    }
    return new SupergalacticPosition(landmark.sgx, landmark.sgy, landmark.sgz);
  }
}
