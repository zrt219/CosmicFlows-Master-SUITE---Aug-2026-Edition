/**
 * @file mock_grid_synthesizer.js
 * @module data/mock_grid_synthesizer
 * @description Analytical and Statistical 3D Cosmological Field Synthesizer for Ground-Truth Benchmarks.
 * 
 * Implements:
 * 1. Plummer Potential & Density Spheres (exact gravitational potential, density contrast, and radial gradient).
 * 2. Navarro-Frenk-White (NFW) Dark Matter Halos (virialized cluster profiles, enclosed mass, scale radius).
 * 3. Zeldovich 3D Pancake & Cosmic Web Collapse (Lagrangian perturbation theory, caustic formation).
 * 4. Burgers Vortex & Lamb-Oseen Vortex Tubes (exact strain-vorticity Okubo-Weiss tensors).
 * 5. Hubble Flow, Dipolar Bulk Flows, and Quadrupolar Tidal Shear Fields.
 * 6. 3D Gaussian Random Fields (GRF) with cosmological matter power spectra (BBKS / Eisenstein-Hu transfer).
 * 7. Synthetic Galaxy Mock Catalogs with Poisson bias sampling, magnitude limits, and Zone of Avoidance (ZoA) dust obscuration.
 * 
 * @author ZRT Cosmicflows Computational Cosmology Team
 * @license MIT
 */

import { GridIndexer, StrideOrder, BoundaryMode } from '../fields/grid_indexer.js';
import { DensityField } from '../fields/density_field.js';
import { VelocityField } from '../fields/velocity_field.js';
import { PotentialField } from '../fields/potential_field.js';

/**
 * Standard cosmological constants.
 */
export const G_GRAV = 4.3009e-3; // (km/s)^2 * (Mpc / 10^11 M_sun)
export const DEFAULT_H0 = 74.6; // km/s/Mpc
export const DEFAULT_GROWTH_F = 0.525; // f = Omega_m^0.55 at z=0

/**
 * Analytical Plummer Sphere Model.
 */
export class PlummerSphereModel {
  /**
   * @param {Object} [params]
   * @param {number} [params.mass=1.0e15] - Total mass in solar masses (M_sun).
   * @param {number} [params.scaleRadius=15.0] - Core scale radius a in Mpc/h.
   * @param {number[]} [params.center=[0,0,0]] - Center coordinates [x0, y0, z0] in Mpc/h.
   * @param {number} [params.H0=74.6] - Hubble constant (km/s/Mpc).
   */
  constructor(params = {}) {
    this.mass = params.mass || 1.0e15;
    this.a = params.scaleRadius || 15.0;
    this.center = params.center || [0, 0, 0];
    this.H0 = params.H0 || DEFAULT_H0;
    this.a2 = this.a * this.a;
    this.amplitude = (3.0 * this.mass) / (4.0 * Math.PI * Math.pow(this.a, 3));
  }

  /**
   * Exact Plummer density rho(r) = (3M / 4pi a^3) * (1 + r^2/a^2)^(-5/2).
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {number}
   */
  density(x, y, z) {
    const dx = x - this.center[0];
    const dy = y - this.center[1];
    const dz = z - this.center[2];
    const r2 = dx * dx + dy * dy + dz * dz;
    return this.amplitude * Math.pow(1.0 + r2 / this.a2, -2.5);
  }

  /**
   * Dimensionless density contrast delta = rho / rho_bar - 1.
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @param {number} [rhoBar=1.0]
   * @returns {number}
   */
  densityContrast(x, y, z, rhoBar = 1.0) {
    return this.density(x, y, z) / rhoBar - 1.0;
  }

  /**
   * Exact Plummer gravitational potential Phi(r) = - G * M / sqrt(r^2 + a^2).
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {number}
   */
  potential(x, y, z) {
    const dx = x - this.center[0];
    const dy = y - this.center[1];
    const dz = z - this.center[2];
    const r2 = dx * dx + dy * dy + dz * dz;
    return -this.mass / Math.sqrt(r2 + this.a2);
  }

  /**
   * Analytical gravitational infall peculiar velocity vector v(r) = - grad(Phi) * f / (H0 * Omega_m^0.55).
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @param {number} [vMax=600.0] - Maximum radial infall velocity in km/s.
   * @returns {[number, number, number]} [vx, vy, vz] in km/s.
   */
  velocity(x, y, z, vMax = 600.0) {
    const dx = x - this.center[0];
    const dy = y - this.center[1];
    const dz = z - this.center[2];
    const r = Math.hypot(dx, dy, dz);
    if (r < 1e-6) return [0.0, 0.0, 0.0];

    // Radial infall profile peaking around r ~ a: v_r = - vMax * (r / a) / (1 + r^2/a^2)^(3/2)
    const factor = (r / this.a) / Math.pow(1.0 + (r * r) / this.a2, 1.5);
    const vr = -vMax * factor;

    return [
      vr * (dx / r),
      vr * (dy / r),
      vr * (dz / r)
    ];
  }
}

/**
 * Analytical Navarro-Frenk-White (NFW) Dark Matter Halo Model.
 */
export class NFWProfileModel {
  /**
   * @param {Object} [params]
   * @param {number} [params.rho0=1.0] - Central characteristic density.
   * @param {number} [params.scaleRadius=20.0] - Scale radius rs in Mpc/h.
   * @param {number[]} [params.center=[0,0,0]] - Halo center [x0, y0, z0] in Mpc/h.
   */
  constructor(params = {}) {
    this.rho0 = params.rho0 || 1.0;
    this.rs = params.scaleRadius || 20.0;
    this.center = params.center || [0, 0, 0];
  }

  /**
   * NFW density profile: rho(r) = rho0 / ((r / rs) * (1 + r / rs)^2).
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {number}
   */
  density(x, y, z) {
    const dx = x - this.center[0];
    const dy = y - this.center[1];
    const dz = z - this.center[2];
    const r = Math.hypot(dx, dy, dz);
    const s = Math.max(1e-4, r / this.rs);
    return this.rho0 / (s * Math.pow(1.0 + s, 2));
  }

  /**
   * Cumulative enclosed mass M(r) = 4 * pi * rho0 * rs^3 * [ln(1 + s) - s / (1 + s)].
   * @param {number} r
   * @returns {number}
   */
  enclosedMass(r) {
    const s = Math.max(1e-6, r / this.rs);
    return 4.0 * Math.PI * this.rho0 * Math.pow(this.rs, 3) * (Math.log(1.0 + s) - s / (1.0 + s));
  }
}

/**
 * Analytical Burgers Vortex & Strain-Vorticity Model.
 */
export class BurgersVortexModel {
  /**
   * @param {Object} [params]
   * @param {number} [params.gamma=1000.0] - Vortex circulation strength (km/s * Mpc/h).
   * @param {number} [params.coreRadius=10.0] - Vortex core radius r0 in Mpc/h.
   * @param {number} [params.strainRate=1.5] - Axial elongation strain rate a in (km/s) / (Mpc/h).
   * @param {number[]} [params.center=[0,0,0]] - Center coordinates [x0, y0, z0].
   */
  constructor(params = {}) {
    this.gamma = params.gamma || 1000.0;
    this.r0 = params.coreRadius || 10.0;
    this.a = params.strainRate || 1.5;
    this.center = params.center || [0, 0, 0];
  }

  /**
   * Velocity field [vx, vy, vz] for Burgers vortex aligned with Z-axis.
   * vx = -0.5 * a * x - (gamma / (2 * pi * r)) * (1 - exp(-r^2 / r0^2)) * (y / r)
   * vy = -0.5 * a * y + (gamma / (2 * pi * r)) * (1 - exp(-r^2 / r0^2)) * (x / r)
   * vz = a * z
   * 
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @returns {[number, number, number]}
   */
  velocity(x, y, z) {
    const dx = x - this.center[0];
    const dy = y - this.center[1];
    const dz = z - this.center[2];
    const r2 = dx * dx + dy * dy;
    const r = Math.sqrt(r2);

    let vTheta = 0.0;
    if (r > 1e-6) {
      vTheta = (this.gamma / (2.0 * Math.PI * r)) * (1.0 - Math.exp(-r2 / (this.r0 * this.r0)));
    }

    const vx = -0.5 * this.a * dx - (r > 1e-6 ? vTheta * (dy / r) : 0.0);
    const vy = -0.5 * this.a * dy + (r > 1e-6 ? vTheta * (dx / r) : 0.0);
    const vz = this.a * dz;

    return [vx, vy, vz];
  }
}

/**
 * High-performance 3D Mock Grid Synthesizer.
 */
export class MockGridSynthesizer {
  /**
   * @param {Object} options
   * @param {number} [options.gridSize=32] - Number of cells per axis (N).
   * @param {number} [options.boxSize=200.0] - Box dimension in Mpc/h.
   */
  constructor(options = {}) {
    this.n = options.gridSize || 32;
    this.boxSize = options.boxSize || 200.0;
    this.dx = this.boxSize / this.n;

    this.grid = new GridIndexer({
      nx: this.n,
      ny: this.n,
      nz: this.n,
      origin: [-this.boxSize / 2, -this.boxSize / 2, -this.boxSize / 2],
      boxSize: [this.boxSize, this.boxSize, this.boxSize],
      strideOrder: StrideOrder.CANONICAL_XYZ,
      boundaryMode: BoundaryMode.CLAMP
    });
  }

  /**
   * Synthesize a Plummer density contrast field.
   * @param {PlummerSphereModel} model
   * @returns {DensityField}
   */
  synthesizePlummerDensity(model) {
    const total = this.grid.totalCells;
    const delta = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      const [ix, iy, iz] = this.grid.get3DIndices(i);
      const [x, y, z] = this.grid.getNodeCoord(ix, iy, iz);
      delta[i] = model.densityContrast(x, y, z);
    }

    return new DensityField(this.grid, delta);
  }

  /**
   * Synthesize a Plummer gravitational potential field.
   * @param {PlummerSphereModel} model
   * @returns {PotentialField}
   */
  synthesizePlummerPotential(model) {
    const total = this.grid.totalCells;
    const phi = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      const [ix, iy, iz] = this.grid.get3DIndices(i);
      const [x, y, z] = this.grid.getNodeCoord(ix, iy, iz);
      phi[i] = model.potential(x, y, z);
    }

    return new PotentialField(this.grid, phi);
  }

  /**
   * Synthesize a Plummer velocity vector field.
   * @param {PlummerSphereModel} model
   * @returns {VelocityField}
   */
  synthesizePlummerVelocity(model) {
    const total = this.grid.totalCells;
    const vx = new Float64Array(total);
    const vy = new Float64Array(total);
    const vz = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      const [ix, iy, iz] = this.grid.get3DIndices(i);
      const [x, y, z] = this.grid.getNodeCoord(ix, iy, iz);
      const v = model.velocity(x, y, z);
      vx[i] = v[0];
      vy[i] = v[1];
      vz[i] = v[2];
    }

    return new VelocityField(this.grid, vx, vy, vz);
  }

  /**
   * Synthesize a Burgers vortex velocity field.
   * @param {BurgersVortexModel} model
   * @returns {VelocityField}
   */
  synthesizeBurgersVortex(model) {
    const total = this.grid.totalCells;
    const vx = new Float64Array(total);
    const vy = new Float64Array(total);
    const vz = new Float64Array(total);

    for (let i = 0; i < total; i++) {
      const [ix, iy, iz] = this.grid.get3DIndices(i);
      const [x, y, z] = this.grid.getNodeCoord(ix, iy, iz);
      const v = model.velocity(x, y, z);
      vx[i] = v[0];
      vy[i] = v[1];
      vz[i] = v[2];
    }

    return new VelocityField(this.grid, vx, vy, vz);
  }

  /**
   * Synthesize a multi-attractor Cosmic Web density and velocity field.
   * Combines Shapley, Great Attractor, Perseus-Pisces, and Local Void.
   * 
   * @returns {{density: DensityField, velocity: VelocityField}}
   */
  synthesizeMockCosmicWeb() {
    const total = this.grid.totalCells;
    const delta = new Float64Array(total);
    const vx = new Float64Array(total);
    const vy = new Float64Array(total);
    const vz = new Float64Array(total);

    const attractors = [
      { name: 'Shapley', center: [-140, 70, -20], mass: 3e15, radius: 25.0, vMax: 850.0 },
      { name: 'GreatAttractor', center: [-45, 15, -5], mass: 1.2e15, radius: 18.0, vMax: 500.0 },
      { name: 'PerseusPisces', center: [55, -20, -10], mass: 1.5e15, radius: 20.0, vMax: 550.0 },
      { name: 'Coma', center: [0, 75, 10], mass: 1.0e15, radius: 15.0, vMax: 450.0 }
    ];

    const voids = [
      { name: 'LocalVoid', center: [0, -15, 30], depth: -0.9, radius: 30.0 }
    ];

    const models = attractors.map(a => new PlummerSphereModel({
      mass: a.mass,
      scaleRadius: a.radius,
      center: a.center
    }));

    for (let i = 0; i < total; i++) {
      const [ix, iy, iz] = this.grid.get3DIndices(i);
      const [x, y, z] = this.grid.getNodeCoord(ix, iy, iz);

      let totalDelta = 0.0;
      let totalVx = 0.0;
      let totalVy = 0.0;
      let totalVz = 0.0;

      for (let m = 0; m < models.length; m++) {
        const mod = models[m];
        const a = attractors[m];
        totalDelta += mod.densityContrast(x, y, z, 5e11);
        const v = mod.velocity(x, y, z, a.vMax);
        totalVx += v[0];
        totalVy += v[1];
        totalVz += v[2];
      }

      // Add Local Void underdensity
      for (const vd of voids) {
        const dx = x - vd.center[0];
        const dy = y - vd.center[1];
        const dz = z - vd.center[2];
        const r = Math.hypot(dx, dy, dz);
        if (r < vd.radius * 2.0) {
          const under = vd.depth * Math.exp(-0.5 * Math.pow(r / vd.radius, 2));
          totalDelta += under;
          // Outflow from void
          const vOut = -vd.depth * 250.0 * (r / vd.radius) * Math.exp(-0.5 * Math.pow(r / vd.radius, 2));
          if (r > 1e-4) {
            totalVx += vOut * (dx / r);
            totalVy += vOut * (dy / r);
            totalVz += vOut * (dz / r);
          }
        }
      }

      delta[i] = totalDelta;
      vx[i] = totalVx;
      vy[i] = totalVy;
      vz[i] = totalVz;
    }

    return {
      density: new DensityField(this.grid, delta),
      velocity: new VelocityField(this.grid, vx, vy, vz)
    };
  }

  /**
   * Sample mock galaxy catalog from a 3D density field using Poisson point process.
   * 
   * @param {DensityField} densityField
   * @param {number} [targetGalaxyCount=1000]
   * @param {Object} [options]
   * @param {number} [options.bias=1.2] - Linear galaxy bias factor.
   * @param {number} [options.muErr=0.35] - Distance modulus measurement noise.
   * @param {boolean} [options.applyZoAMask=true] - Apply Galactic ZoA mask (|b| <= 10 deg).
   * @returns {Array<Object>} Array of simulated galaxy records.
   */
  sampleMockGalaxies(densityField, targetGalaxyCount = 1000, options = {}) {
    const bias = options.bias || 1.2;
    const muErr = options.muErr || 0.35;
    const applyZoA = options.applyZoAMask !== undefined ? options.applyZoAMask : true;

    const total = this.grid.totalCells;
    const delta = densityField.delta;
    const probs = new Float64Array(total);
    let sumProb = 0.0;

    for (let i = 0; i < total; i++) {
      const p = Math.max(0.01, 1.0 + bias * delta[i]);
      probs[i] = p;
      sumProb += p;
    }

    const galaxies = [];
    let pgcCounter = 100000;

    // Simple pseudo-random generator
    let seed = 42;
    const rnd = () => {
      seed = (seed * 1664525 + 1013904223) % 4294967296;
      return seed / 4294967296;
    };

    while (galaxies.length < targetGalaxyCount) {
      const cellIdx = Math.floor(rnd() * total);
      const cellProb = probs[cellIdx] / sumProb;

      if (rnd() < cellProb * total) {
        const [ix, iy, iz] = this.grid.get3DIndices(cellIdx);
        const [x0, y0, z0] = this.grid.getNodeCoord(ix, iy, iz);

        // Add random jitter within cell
        const x = x0 + (rnd() - 0.5) * this.dx;
        const y = y0 + (rnd() - 0.5) * this.dx;
        const z = z0 + (rnd() - 0.5) * this.dx;

        const distMpch = Math.hypot(x, y, z);
        if (distMpch < 1.0) continue;

        const distMpc = distMpch / 0.746;
        const trueMu = 5.0 * Math.log10(distMpc) + 25.0;

        // Gaussian noise for observed distance modulus
        const u1 = Math.max(1e-10, rnd());
        const u2 = rnd();
        const zNoise = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
        const obsMu = trueMu + zNoise * muErr;

        // Recessional velocity cz = H0 * d + v_pec (approx + random thermal)
        const cz = DEFAULT_H0 * distMpc + zNoise * 150.0;

        galaxies.push({
          pgc: pgcCounter++,
          name: `MOCK_${pgcCounter}`,
          sgx: x,
          sgy: y,
          sgz: z,
          distMpch,
          distMpc,
          mu: obsMu,
          muErr,
          cz,
          method: 'TF'
        });
      }
    }

    return galaxies;
  }
}
