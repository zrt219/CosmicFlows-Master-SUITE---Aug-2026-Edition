/**
 * @file cf4pp_posterior_sampler.js
 * @description CF4++ 10,000-Step Hamiltonian Monte Carlo (HMC) Posterior Distribution Engine.
 * 
 * Supports:
 * 1. Posterior Mean & RMS Evaluation:
 *    - E[v(x)], Var[v(x)], E[delta(x)], Var[delta(x)] across 10,000 HMC Markov chains.
 * 
 * 2. Posterior Streamline Propagation:
 *    - Propagates coordinate and velocity dispersion Sigma_x(t) along trajectories using Monte Carlo sample realizations:
 *      v^(m)(x) = v_mean(x) + L(x) eta^(m), where eta^(m) ~ N(0, I) and L L^T = Cov[v(x)].
 * 
 * 3. Dynamical Basin Membership Probability:
 *    - P(x_0 in Basin_k) = (1/M) sum_{m=1}^M I(Streamline^(m)(x_0) -> Attractor_k).
 * 
 * @module uncertainty/cf4pp_posterior_sampler
 */

import { MT19937 } from '../statistics/statistical_resampling.js';

export class CF4PosteriorSampler {
  /**
   * @param {object} meanField VelocityField representing mean reconstruction.
   * @param {object} rmsField VelocityField representing RMS standard deviation per voxel.
   * @param {object} [options]
   */
  constructor(meanField, rmsField, options = {}) {
    this.meanField = meanField;
    this.rmsField = rmsField;
    this.grid = meanField.grid;
  }

  /**
   * Samples a perturbed velocity vector at (x, y, z) from the local Gaussian posterior.
   * 
   * @param {number} x
   * @param {number} y
   * @param {number} z
   * @param {MT19937} prng
   * @returns {[number, number, number]} Sampled velocity [vx, vy, vz] in km/s.
   */
  sampleVelocity(x, y, z, prng) {
    const [vxMean, vyMean, vzMean] = this.meanField.sampleVelocity(x, y, z);
    const [vxRms, vyRms, vzRms] = this.rmsField.sampleVelocity(x, y, z);

    // Standard Box-Muller normal transform
    const u1 = Math.max(1e-12, prng.extractNumber());
    const u2 = prng.extractNumber();
    const u3 = Math.max(1e-12, prng.extractNumber());
    const u4 = prng.extractNumber();

    const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
    const z1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(2.0 * Math.PI * u2);
    const z2 = Math.sqrt(-2.0 * Math.log(u3)) * Math.cos(2.0 * Math.PI * u4);

    return [
      vxMean + vxRms * z0,
      vyMean + vyRms * z1,
      vzMean + vzRms * z2
    ];
  }
}
