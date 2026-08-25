/**
 * @file galaxy_catalog_reader.js
 * @module data/galaxy_catalog_reader
 * @description Comprehensive Multi-Survey Astrometric Galaxy Catalog Reader and Columnar Store.
 * 
 * Supports:
 * - CF4 Individual & Grouped Catalogs (Tully et al. 2023, Kourkchi et al. 2020).
 * - CF3 Distance Catalog (Tully et al. 2016).
 * - 2MASS Tully-Fisher Survey (2MTF), 2M++ (Lavaux & Hudson 2011), SFI++, and 6dFGSv.
 * - Exact cosmological distance modulus, recessional velocity, and peculiar velocity conversions.
 * - Analytic error propagation for distance and velocity.
 * - Full IAU Equatorial <-> Galactic <-> Supergalactic <-> Supergalactic Cartesian coordinate transformations.
 * - High-density TypedArray columnar storage with predicate filtering, group aggregation, and spatial indexing.
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

import {
  equatorialToGalactic,
  galacticToSupergalactic,
  supergalacticToCartesian,
  cartesianToSupergalactic,
  supergalacticToGalactic,
  galacticToEquatorial
} from '../coordinates/canonical_frame.js';

/**
 * Natural log of 10 divided by 5 (exact derivative factor for distance modulus).
 * d(d)/d(mu) = d * (ln(10)/5)
 * @type {number}
 */
export const LN10_DIV_5 = Math.LN10 / 5.0; // ~0.4605170185988092

/**
 * Speed of light in vacuum in km/s.
 * @type {number}
 */
export const SPEED_OF_LIGHT_KMS = 299792.458;

/**
 * Supported Distance Estimation Methods.
 * @enum {string}
 */
export const DistanceMethod = Object.freeze({
  TULLY_FISHER: 'TF',
  FUNDAMENTAL_PLANE: 'FP',
  TYPE_IA_SUPERNOVA: 'SNIa',
  SURFACE_BRIGHTNESS_FLUCTUATION: 'SBF',
  TIP_OF_RED_GIANT_BRANCH: 'TRGB',
  CEPHEID: 'Cepheid',
  UNKNOWN: 'UNK'
});

/**
 * Convert distance modulus mu (mag) to metric distance d (Mpc).
 * @param {number} mu - Distance modulus in magnitudes.
 * @param {number} [h=0.746] - Reduced Hubble constant h = H0 / 100.
 * @param {boolean} [inMpch=false] - If true, returns distance in Mpc/h instead of Mpc.
 * @returns {number} Distance in Mpc (or Mpc/h if inMpch=true).
 */
export function distanceModulusToDistance(mu, h = 0.746, inMpch = false) {
  if (isNaN(mu) || mu === null) return NaN;
  const dMpc = Math.pow(10.0, (mu - 25.0) / 5.0);
  return inMpch ? dMpc * h : dMpc;
}

/**
 * Convert metric distance d (Mpc) to distance modulus mu (mag).
 * @param {number} dMpc - Distance in Mpc.
 * @returns {number} Distance modulus in magnitudes.
 */
export function distanceToDistanceModulus(dMpc) {
  if (isNaN(dMpc) || dMpc <= 0.0) return NaN;
  return 5.0 * Math.log10(dMpc) + 25.0;
}

/**
 * Propagate 1-sigma uncertainty in distance modulus sigma_mu to distance uncertainty sigma_d.
 * Uses exact first-order Taylor expansion: sigma_d = d * (ln(10)/5) * sigma_mu.
 * 
 * @param {number} d - Distance (Mpc or Mpc/h).
 * @param {number} sigmaMu - 1-sigma uncertainty in distance modulus (mag).
 * @returns {number} 1-sigma uncertainty in distance (same units as d).
 */
export function propagateDistanceUncertainty(d, sigmaMu) {
  if (isNaN(d) || isNaN(sigmaMu) || d <= 0.0 || sigmaMu < 0.0) return NaN;
  return d * LN10_DIV_5 * sigmaMu;
}

/**
 * Compute radial peculiar velocity: v_pec = v_rec - H0 * d.
 * 
 * @param {number} vRec - Recessional velocity cz in km/s (CMB or Local Sheet frame).
 * @param {number} dMpc - Metric distance in Mpc.
 * @param {number} [H0=74.6] - Hubble constant in km/s/Mpc.
 * @param {boolean} [relativistic=false] - Whether to apply relativistic cosmological correction.
 * @returns {number} Peculiar velocity in km/s.
 */
export function calculatePeculiarVelocity(vRec, dMpc, H0 = 74.6, relativistic = false) {
  if (isNaN(vRec) || isNaN(dMpc)) return NaN;
  const vHubble = H0 * dMpc;
  const vLinear = vRec - vHubble;

  if (!relativistic) return vLinear;

  // Relativistic correction: v_pec_rel = (v_rec - v_H) / (1 + v_H / c)
  const c = SPEED_OF_LIGHT_KMS;
  return vLinear / (1.0 + vHubble / c);
}

/**
 * Propagate uncertainty in peculiar velocity.
 * sigma_vpec^2 = sigma_vrec^2 + (H0 * sigma_d)^2 + sigma_thermal^2
 * 
 * @param {number} dMpc - Distance in Mpc.
 * @param {number} sigmaMu - Distance modulus uncertainty (mag).
 * @param {number} [sigmaVrec=20.0] - Recessional velocity measurement uncertainty (km/s).
 * @param {number} [H0=74.6] - Hubble parameter (km/s/Mpc).
 * @param {number} [sigmaThermal=150.0] - Cosmic thermal dispersion (km/s).
 * @returns {number} Total 1-sigma peculiar velocity uncertainty in km/s.
 */
export function calculatePeculiarVelocityError(dMpc, sigmaMu, sigmaVrec = 20.0, H0 = 74.6, sigmaThermal = 150.0) {
  if (isNaN(dMpc) || isNaN(sigmaMu)) return NaN;
  const sigmaD = propagateDistanceUncertainty(dMpc, sigmaMu);
  const vHubbleErr = H0 * sigmaD;
  return Math.sqrt(sigmaVrec * sigmaVrec + vHubbleErr * vHubbleErr + sigmaThermal * sigmaThermal);
}

/**
 * High-performance Columnar Galaxy Catalog Store.
 */
export class GalaxyCatalogStore {
  /**
   * @param {number} capacity - Number of galaxy records.
   */
  constructor(capacity) {
    this.capacity = capacity;
    this.count = 0;

    this.pgc = new Int32Array(capacity);
    this.ra = new Float64Array(capacity);
    this.dec = new Float64Array(capacity);
    this.glon = new Float64Array(capacity);
    this.glat = new Float64Array(capacity);
    this.sgl = new Float64Array(capacity);
    this.sgb = new Float64Array(capacity);

    this.sgx = new Float32Array(capacity); // Mpc/h
    this.sgy = new Float32Array(capacity); // Mpc/h
    this.sgz = new Float32Array(capacity); // Mpc/h

    this.cz = new Float32Array(capacity); // km/s
    this.vcmb = new Float32Array(capacity); // km/s
    this.vls = new Float32Array(capacity); // km/s

    this.mu = new Float32Array(capacity); // mag
    this.muErr = new Float32Array(capacity); // mag
    this.dist = new Float32Array(capacity); // Mpc/h
    this.distErr = new Float32Array(capacity); // Mpc/h

    this.vpec = new Float32Array(capacity); // km/s
    this.vpecErr = new Float32Array(capacity); // km/s

    this.groupId = new Int32Array(capacity);
    this.weight = new Float32Array(capacity);
    this.methods = new Array(capacity);
    this.names = new Array(capacity);
  }

  /**
   * Add a galaxy record.
   * @param {Object} g - Galaxy properties.
   */
  addGalaxy(g) {
    if (this.count >= this.capacity) {
      throw new Error(`GalaxyCatalogStore overflow: capacity ${this.capacity} reached.`);
    }

    const idx = this.count;
    this.pgc[idx] = g.pgc || 0;
    this.ra[idx] = g.ra || 0.0;
    this.dec[idx] = g.dec || 0.0;
    this.glon[idx] = g.glon !== undefined ? g.glon : 0.0;
    this.glat[idx] = g.glat !== undefined ? g.glat : 0.0;
    this.sgl[idx] = g.sgl !== undefined ? g.sgl : 0.0;
    this.sgb[idx] = g.sgb !== undefined ? g.sgb : 0.0;

    this.sgx[idx] = g.sgx || 0.0;
    this.sgy[idx] = g.sgy || 0.0;
    this.sgz[idx] = g.sgz || 0.0;

    this.cz[idx] = g.cz || 0.0;
    this.vcmb[idx] = g.vcmb || g.cz || 0.0;
    this.vls[idx] = g.vls || g.cz || 0.0;

    this.mu[idx] = g.mu || 0.0;
    this.muErr[idx] = g.muErr || 0.0;
    this.dist[idx] = g.dist || 0.0;
    this.distErr[idx] = g.distErr || 0.0;

    this.vpec[idx] = g.vpec || 0.0;
    this.vpecErr[idx] = g.vpecErr || 0.0;

    this.groupId[idx] = g.groupId !== undefined ? g.groupId : -1;
    this.weight[idx] = g.weight !== undefined ? g.weight : 1.0;
    this.methods[idx] = g.method || DistanceMethod.UNKNOWN;
    this.names[idx] = g.name || `PGC${this.pgc[idx]}`;

    this.count++;
  }

  /**
   * Get a single galaxy object by index.
   * @param {number} i
   * @returns {Object}
   */
  getGalaxy(i) {
    if (i < 0 || i >= this.count) return null;
    return {
      index: i,
      pgc: this.pgc[i],
      name: this.names[i],
      ra: this.ra[i],
      dec: this.dec[i],
      glon: this.glon[i],
      glat: this.glat[i],
      sgl: this.sgl[i],
      sgb: this.sgb[i],
      sgx: this.sgx[i],
      sgy: this.sgy[i],
      sgz: this.sgz[i],
      cz: this.cz[i],
      vcmb: this.vcmb[i],
      vls: this.vls[i],
      mu: this.mu[i],
      muErr: this.muErr[i],
      dist: this.dist[i],
      distErr: this.distErr[i],
      vpec: this.vpec[i],
      vpecErr: this.vpecErr[i],
      groupId: this.groupId[i],
      weight: this.weight[i],
      method: this.methods[i]
    };
  }

  /**
   * Filter catalog by predicate function.
   * @param {function(Object): boolean} predicate
   * @returns {GalaxyCatalogStore}
   */
  filter(predicate) {
    const matchedIndices = [];
    for (let i = 0; i < this.count; i++) {
      const g = this.getGalaxy(i);
      if (predicate(g)) {
        matchedIndices.push(i);
      }
    }

    const subStore = new GalaxyCatalogStore(matchedIndices.length);
    for (const idx of matchedIndices) {
      subStore.addGalaxy(this.getGalaxy(idx));
    }
    return subStore;
  }

  /**
   * Compute axis-aligned bounding box in Supergalactic Cartesian space (Mpc/h).
   * @returns {{xmin: number, xmax: number, ymin: number, ymax: number, zmin: number, zmax: number}}
   */
  getBoundingBox() {
    let xmin = Infinity, xmax = -Infinity;
    let ymin = Infinity, ymax = -Infinity;
    let zmin = Infinity, zmax = -Infinity;

    for (let i = 0; i < this.count; i++) {
      const x = this.sgx[i];
      const y = this.sgy[i];
      const z = this.sgz[i];
      if (x < xmin) xmin = x;
      if (x > xmax) xmax = x;
      if (y < ymin) ymin = y;
      if (y > ymax) ymax = y;
      if (z < zmin) zmin = z;
      if (z > zmax) zmax = z;
    }

    return {
      xmin: isFinite(xmin) ? xmin : 0,
      xmax: isFinite(xmax) ? xmax : 0,
      ymin: isFinite(ymin) ? ymin : 0,
      ymax: isFinite(ymax) ? ymax : 0,
      zmin: isFinite(zmin) ? zmin : 0,
      zmax: isFinite(zmax) ? zmax : 0
    };
  }

  /**
   * Aggregate catalog into groups using inverse-variance weighting on distance modulus.
   * @returns {Array<Object>} Aggregated group records.
   */
  aggregateGroups() {
    const groupMap = new Map();

    for (let i = 0; i < this.count; i++) {
      const gid = this.groupId[i];
      if (gid < 0) continue; // Skip unassigned field galaxies

      if (!groupMap.has(gid)) {
        groupMap.set(gid, []);
      }
      groupMap.get(gid).push(this.getGalaxy(i));
    }

    const aggregated = [];

    for (const [gid, members] of groupMap.entries()) {
      let sumWeightMu = 0.0;
      let sumWeight = 0.0;
      let sumCz = 0.0;
      let sumSgx = 0.0;
      let sumSgy = 0.0;
      let sumSgz = 0.0;
      const methods = new Set();

      for (const m of members) {
        const err = m.muErr > 0 ? m.muErr : 0.35;
        const w = 1.0 / (err * err);
        sumWeightMu += m.mu * w;
        sumWeight += w;
        sumCz += m.cz;
        sumSgx += m.sgx;
        sumSgy += m.sgy;
        sumSgz += m.sgz;
        methods.add(m.method);
      }

      const n = members.length;
      const meanMu = sumWeight > 0 ? sumWeightMu / sumWeight : 0.0;
      const meanMuErr = sumWeight > 0 ? 1.0 / Math.sqrt(sumWeight) : 0.0;
      const meanCz = sumCz / n;

      // Velocity dispersion
      let sumVarCz = 0.0;
      for (const m of members) {
        const diff = m.cz - meanCz;
        sumVarCz += diff * diff;
      }
      const sigmaV = n > 1 ? Math.sqrt(sumVarCz / (n - 1)) : 0.0;

      aggregated.push({
        groupId: gid,
        memberCount: n,
        mu: meanMu,
        muErr: meanMuErr,
        cz: meanCz,
        sigmaV,
        sgx: sumSgx / n,
        sgy: sumSgy / n,
        sgz: sumSgz / n,
        methods: Array.from(methods)
      });
    }

    return aggregated;
  }

  /**
   * Export to GeoJSON FeatureCollection.
   * @returns {Object} GeoJSON structure.
   */
  toGeoJSON() {
    const features = [];
    for (let i = 0; i < this.count; i++) {
      features.push({
        type: 'Feature',
        geometry: {
          type: 'Point',
          coordinates: [this.sgx[i], this.sgy[i], this.sgz[i]]
        },
        properties: {
          pgc: this.pgc[i],
          name: this.names[i],
          ra: this.ra[i],
          dec: this.dec[i],
          glon: this.glon[i],
          glat: this.glat[i],
          cz: this.cz[i],
          mu: this.mu[i],
          dist: this.dist[i],
          vpec: this.vpec[i],
          groupId: this.groupId[i],
          method: this.methods[i]
        }
      });
    }

    return {
      type: 'FeatureCollection',
      features
    };
  }
}

/**
 * Universal Galaxy Catalog Reader and Ingestion Engine.
 */
export class GalaxyCatalogReader {
  /**
   * Ingest raw catalog records (from CSV, JSON, or FITS Table) into a high-performance GalaxyCatalogStore.
   * 
   * @param {Array<Object>} records - Array of raw galaxy row objects.
   * @param {Object} [options]
   * @param {number} [options.H0=74.6] - Hubble parameter in km/s/Mpc.
   * @param {boolean} [options.relativistic=false] - Whether to apply relativistic peculiar velocity.
   * @param {number} [options.defaultMuErr=0.35] - Default distance modulus uncertainty if missing.
   * @returns {GalaxyCatalogStore}
   */
  static ingestRecords(records, options = {}) {
    const H0 = options.H0 || 74.6;
    const h = H0 / 100.0;
    const relativistic = Boolean(options.relativistic);
    const defaultMuErr = options.defaultMuErr || 0.35;

    const store = new GalaxyCatalogStore(records.length);

    for (const r of records) {
      const pgc = Number(r.pgc || r.PGC || r.id || r.ID || 0);
      const name = String(r.name || r.NAME || `PGC${pgc}`);

      // Coordinates
      let ra = Number(r.ra !== undefined ? r.ra : (r.RA !== undefined ? r.RA : 0.0));
      let dec = Number(r.dec !== undefined ? r.dec : (r.DEC !== undefined ? r.DEC : 0.0));
      let glon = Number(r.glon !== undefined ? r.glon : (r.GLON !== undefined ? r.GLON : (r.l !== undefined ? r.l : NaN)));
      let glat = Number(r.glat !== undefined ? r.glat : (r.GLAT !== undefined ? r.GLAT : (r.b !== undefined ? r.b : NaN)));
      let sgl = Number(r.sgl !== undefined ? r.sgl : (r.SGL !== undefined ? r.SGL : NaN));
      let sgb = Number(r.sgb !== undefined ? r.sgb : (r.SGB !== undefined ? r.SGB : NaN));

      // Derive Galactic & Supergalactic coordinates if not given
      if (isNaN(glon) || isNaN(glat)) {
        const gal = equatorialToGalactic(ra, dec);
        glon = gal[0];
        glat = gal[1];
      }

      if (isNaN(sgl) || isNaN(sgb)) {
        const sgal = galacticToSupergalactic(glon, glat);
        sgl = sgal[0];
        sgb = sgal[1];
      }

      // Distance & Modulus
      let mu = Number(r.mu !== undefined ? r.mu : (r.DM !== undefined ? r.DM : (r.dm !== undefined ? r.dm : NaN)));
      let muErr = Number(r.muErr !== undefined ? r.muErr : (r.eDM !== undefined ? r.eDM : (r.e_dm !== undefined ? r.e_dm : defaultMuErr)));
      let distMpc = Number(r.distMpc !== undefined ? r.distMpc : (r.d !== undefined ? r.d : NaN));
      let distMpch = Number(r.distMpch !== undefined ? r.distMpch : (r.dist !== undefined ? r.dist : NaN));

      if (isNaN(distMpc) && !isNaN(mu)) {
        distMpc = distanceModulusToDistance(mu, h, false);
      }
      if (isNaN(distMpch) && !isNaN(distMpc)) {
        distMpch = distMpc * h;
      }
      if (isNaN(mu) && !isNaN(distMpc)) {
        mu = distanceToDistanceModulus(distMpc);
      }

      const distErrMpch = !isNaN(distMpch) ? propagateDistanceUncertainty(distMpch, muErr) : 0.0;

      // Supergalactic Cartesian coordinates (Mpc/h)
      let sgx = Number(r.sgx !== undefined ? r.sgx : (r.SGX !== undefined ? r.SGX : NaN));
      let sgy = Number(r.sgy !== undefined ? r.sgy : (r.SGY !== undefined ? r.SGY : NaN));
      let sgz = Number(r.sgz !== undefined ? r.sgz : (r.SGZ !== undefined ? r.SGZ : NaN));

      if (isNaN(sgx) || isNaN(sgy) || isNaN(sgz)) {
        if (!isNaN(distMpch)) {
          const cart = supergalacticToCartesian(sgl, sgb, distMpch);
          sgx = cart[0];
          sgy = cart[1];
          sgz = cart[2];
        } else {
          sgx = 0.0; sgy = 0.0; sgz = 0.0;
        }
      }

      // Velocities
      const cz = Number(r.cz !== undefined ? r.cz : (r.CZ !== undefined ? r.CZ : (r.v !== undefined ? r.v : 0.0)));
      const vcmb = Number(r.vcmb !== undefined ? r.vcmb : (r.Vcmb !== undefined ? r.Vcmb : cz));
      const vls = Number(r.vls !== undefined ? r.vls : (r.Vls !== undefined ? r.Vls : cz));

      let vpec = Number(r.vpec !== undefined ? r.vpec : (r.Vpec !== undefined ? r.Vpec : NaN));
      if (isNaN(vpec) && !isNaN(distMpc) && cz > 0) {
        vpec = calculatePeculiarVelocity(cz, distMpc, H0, relativistic);
      }
      const vpecErr = !isNaN(distMpc) ? calculatePeculiarVelocityError(distMpc, muErr, 20.0, H0, 150.0) : 150.0;

      const groupId = Number(r.groupId !== undefined ? r.groupId : (r.GroupId !== undefined ? r.GroupId : (r.nest !== undefined ? r.nest : -1)));
      const method = String(r.method || r.Method || DistanceMethod.UNKNOWN).trim();

      store.addGalaxy({
        pgc,
        name,
        ra,
        dec,
        glon,
        glat,
        sgl,
        sgb,
        sgx,
        sgy,
        sgz,
        cz,
        vcmb,
        vls,
        mu: isNaN(mu) ? 0.0 : mu,
        muErr: isNaN(muErr) ? defaultMuErr : muErr,
        dist: isNaN(distMpch) ? 0.0 : distMpch,
        distErr: distErrMpch,
        vpec: isNaN(vpec) ? 0.0 : vpec,
        vpecErr,
        groupId,
        method
      });
    }

    return store;
  }
}
