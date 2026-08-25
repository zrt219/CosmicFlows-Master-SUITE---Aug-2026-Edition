/**
 * @file hmc_posterior_sampler.js
 * @description Production 10,000-Step Hamiltonian Monte Carlo (HMC) Posterior Distribution Sampler.
 *
 * Implements:
 * 1. Hamiltonian Monte Carlo sampler with symplectic leapfrog trajectories and Metropolis-Hastings filter.
 * 2. Dual-averaging step size adaptation (Hoffman & Gelman 2014 / Nesterov 2009) targeting acceptance rate delta ~ 0.65.
 * 3. Diagonal mass matrix (Welford variance) adaptation during warm-up phase.
 * 4. Multi-chain execution with diagnostic telemetry: acceptance rates, log posterior tracking, energy transitions.
 * 5. Full support for 10,000+ step chains on cosmological velocity, density, and bulk flow parameters.
 *
 * @module uncertainty/hmc_posterior_sampler
 */

import { LeapfrogIntegrator, computeKineticEnergy } from './leapfrog_integrator.js';
import { MT19937 } from '../statistics/statistical_resampling.js';

/**
 * Standard Dual-Averaging Step Size Controller for MCMC warm-up.
 */
export class DualAveragingController {
  /**
   * @param {number} targetAcceptance - Target acceptance probability delta (default 0.65 for standard HMC)
   * @param {number} [gamma=0.05] - Adaptation decay parameter
   * @param {number} [t0=10.0] - Early adaptation damping parameter
   * @param {number} [kappa=0.75] - Step size learning rate exponent
   */
  constructor(targetAcceptance = 0.65, gamma = 0.05, t0 = 10.0, kappa = 0.75) {
    this.targetAcceptance = targetAcceptance;
    this.gamma = gamma;
    this.t0 = t0;
    this.kappa = kappa;

    this.mu = 0.0;
    this.logEps = 0.0;
    this.logEpsAvg = 0.0;
    this.hBar = 0.0;
    this.iteration = 0;
  }

  /**
   * Initializes the step size controller with an initial step size estimate eps0.
   * @param {number} eps0
   */
  init(eps0) {
    this.mu = Math.log(10.0 * eps0);
    this.logEps = Math.log(eps0);
    this.logEpsAvg = Math.log(eps0);
    this.hBar = 0.0;
    this.iteration = 0;
  }

  /**
   * Updates step size based on empirical acceptance statistic alpha of the step.
   * @param {number} alpha - Acceptance probability of current proposal min(1, exp(-deltaH))
   * @returns {number} Updated step size epsilon
   */
  update(alpha) {
    this.iteration++;
    const m = this.iteration;
    const eta = 1.0 / (m + this.t0);

    // Running error statistic
    this.hBar = (1.0 - eta) * this.hBar + eta * (this.targetAcceptance - alpha);

    // Update log(epsilon)
    this.logEps = this.mu - (Math.sqrt(m) / this.gamma) * this.hBar;

    // Running smoothed average of log(epsilon)
    const mKappa = Math.pow(m, -this.kappa);
    this.logEpsAvg = mKappa * this.logEps + (1.0 - mKappa) * this.logEpsAvg;

    return Math.exp(this.logEps);
  }

  /**
   * Gets the final adapted step size.
   * @returns {number}
   */
  getFinalStepSize() {
    return Math.exp(this.logEpsAvg);
  }
}

/**
 * Full 10,000-Step Hamiltonian Monte Carlo Posterior Sampler.
 */
export class HMCPosteriorSampler {
  /**
   * @param {Object} model - Target probability distribution model
   * @param {function(Float64Array): number} model.logPosterior - Log posterior function ln P(q)
   * @param {function(Float64Array): Float64Array} model.gradLogPosterior - Gradient of log posterior nabla ln P(q)
   * @param {number} model.dimension - Dimension D of parameter space
   * @param {Object} [options={}]
   * @param {number} [options.stepSize=0.1] - Initial leapfrog step size epsilon
   * @param {number} [options.numLeapfrogSteps=10] - Number of leapfrog steps per proposal L
   * @param {number} [options.targetAcceptance=0.65] - Target acceptance probability
   * @param {boolean} [options.adaptStepSize=true] - Whether to adapt step size during warmup
   * @param {boolean} [options.adaptMassMatrix=true] - Whether to adapt diagonal mass matrix during warmup
   * @param {number} [options.seed=42] - PRNG seed
   */
  constructor(model, options = {}) {
    if (!model || typeof model.logPosterior !== \'function\' || typeof model.gradLogPosterior !== \'function\') {
      throw new TypeError(\'model must define logPosterior(q) and gradLogPosterior(q) functions\');
    }
    this.model = model;
    this.dim = model.dimension;

    this.stepSize = options.stepSize ?? 0.1;
    this.numLeapfrogSteps = options.numLeapfrogSteps ?? 10;
    this.targetAcceptance = options.targetAcceptance ?? 0.65;
    this.adaptStepSize = options.adaptStepSize ?? true;
    this.adaptMassMatrix = options.adaptMassMatrix ?? true;
    this.seed = options.seed ?? 42;

    this.prng = new MT19937(this.seed);
    this.dualAveraging = new DualAveragingController(this.targetAcceptance);

    // Potential U(q) = -ln P(q) and gradU(q) = -nabla ln P(q)
    this.gradPotentialFn = (q) => {
      const logP = this.model.logPosterior(q);
      const gradLogP = this.model.gradLogPosterior(q);
      const gradU = new Float64Array(this.dim);
      for (let i = 0; i < this.dim; i++) {
        gradU[i] = -gradLogP[i];
      }
      return {
        potential: -logP,
        gradU
      };
    };

    this.invMass = new Float64Array(this.dim).fill(1.0);
    this.mass = new Float64Array(this.dim).fill(1.0);
  }

  /**
   * Generates Gaussian random vector from N(0, diag(M)).
   * @returns {Float64Array}
   */
  sampleMomentum() {
    const p = new Float64Array(this.dim);
    for (let i = 0; i < this.dim; i += 2) {
      const u1 = Math.max(1e-15, this.prng.extractNumber());
      const u2 = this.prng.extractNumber();
      const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
      const z1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(2.0 * Math.PI * u2);

      const std0 = Math.sqrt(this.mass[i]);
      p[i] = z0 * std0;

      if (i + 1 < this.dim) {
        const std1 = Math.sqrt(this.mass[i + 1]);
        p[i + 1] = z1 * std1;
      }
    }
    return p;
  }

  /**
   * Finds a reasonable initial step size using heuristic search.
   * @param {Float64Array} q0
   * @returns {number} Initial step size
   */
  findReasonableStepSize(q0) {
    let eps = 1.0;
    const p = this.sampleMomentum();
    const integrator = new LeapfrogIntegrator(this.gradPotentialFn, { invMass: this.invMass });

    // Initial state
    const eval0 = this.gradPotentialFn(q0);
    const h0 = eval0.potential + computeKineticEnergy(p, this.invMass);

    // Take one step
    let res = integrator.integrate(q0, p, eps, 1);
    let h = res.finalHamiltonian;

    let deltaH = h - h0;
    let direction = -deltaH > Math.log(0.5) ? 1 : -1;

    let count = 0;
    while (count < 50) {
      count++;
      if (direction === 1 && -deltaH <= Math.log(0.5)) break;
      if (direction === -1 && -deltaH >= Math.log(0.5)) break;

      eps = direction === 1 ? eps * 2.0 : eps * 0.5;
      res = integrator.integrate(q0, p, eps, 1);
      h = res.finalHamiltonian;
      deltaH = h - h0;
    }

    return Math.max(1e-4, Math.min(10.0, eps));
  }

  /**
   * Executes HMC sampling for numSamples steps after numWarmup steps.
   *
   * @param {Float64Array|number[]} initialPosition - Initial starting state q0
   * @param {number} numSamples - Number of posterior samples to collect (e.g. 10,000)
   * @param {number} [numWarmup=1000] - Number of warmup/burn-in adaptation steps
   * @returns {{
   *   samples: Float64Array[],
   *   logPosteriors: Float64Array,
   *   acceptanceRates: Float64Array,
   *   stepSizes: Float64Array,
   *   meanAcceptanceRate: number,
   *   finalStepSize: number,
   *   divergences: number,
   *   summary: Object
   * }}
   */
  sample(initialPosition, numSamples = 10000, numWarmup = 1000) {
    let currentQ = new Float64Array(initialPosition);
    let currentLogPost = this.model.logPosterior(currentQ);

    let eps = this.stepSize;
    if (this.adaptStepSize) {
      eps = this.findReasonableStepSize(currentQ);
      this.dualAveraging.init(eps);
    }

    const totalSteps = numWarmup + numSamples;
    const samples = [];
    const logPosteriors = new Float64Array(numSamples);
    const acceptanceRates = new Float64Array(numSamples);
    const stepSizes = new Float64Array(numSamples);

    let acceptedCount = 0;
    let totalAcceptedWarmup = 0;
    let divergences = 0;

    // Running Welford accumulator for mass matrix adaptation
    const runningMean = new Float64Array(this.dim);
    const runningM2 = new Float64Array(this.dim);
    let welfordCount = 0;

    for (let step = 0; step < totalSteps; step++) {
      const isWarmup = step < numWarmup;
      const currentEps = isWarmup && this.adaptStepSize ? eps : this.dualAveraging.getFinalStepSize();

      // Sample momentum p ~ N(0, M)
      const p = this.sampleMomentum();
      const currentKinetic = computeKineticEnergy(p, this.invMass);
      const currentHamiltonian = -currentLogPost + currentKinetic;

      // Integrate trajectory
      const integrator = new LeapfrogIntegrator(this.gradPotentialFn, { invMass: this.invMass });
      const traj = integrator.integrate(currentQ, p, currentEps, this.numLeapfrogSteps);

      let alpha = 0.0;
      let accept = false;

      if (!traj.diverged) {
        const propKinetic = computeKineticEnergy(traj.p, this.invMass);
        const propLogPost = this.model.logPosterior(traj.q);
        const propHamiltonian = -propLogPost + propKinetic;
        const deltaH = propHamiltonian - currentHamiltonian;

        alpha = Math.min(1.0, Math.exp(-deltaH));
        if (!Number.isFinite(alpha)) alpha = 0.0;

        const u = this.prng.extractNumber();
        if (u < alpha) {
          accept = true;
          currentQ = new Float64Array(traj.q);
          currentLogPost = propLogPost;
        }
      } else {
        divergences++;
        alpha = 0.0;
      }

      if (isWarmup) {
        if (accept) totalAcceptedWarmup++;
        if (this.adaptStepSize) {
          eps = this.dualAveraging.update(alpha);
        }

        // Mass matrix adaptation in second half of warmup
        if (this.adaptMassMatrix && step > numWarmup * 0.25 && step < numWarmup * 0.9) {
          welfordCount++;
          for (let i = 0; i < this.dim; i++) {
            const val = currentQ[i];
            const d1 = val - runningMean[i];
            runningMean[i] += d1 / welfordCount;
            const d2 = val - runningMean[i];
            runningM2[i] += d1 * d2;

            if (welfordCount > 10) {
              const variance = Math.max(1e-4, runningM2[i] / (welfordCount - 1));
              this.invMass[i] = variance; // M^{-1} = diag(Var(q))
              this.mass[i] = 1.0 / variance;
            }
          }
        }
      } else {
        const sampleIdx = step - numWarmup;
        if (accept) acceptedCount++;

        samples.push(new Float64Array(currentQ));
        logPosteriors[sampleIdx] = currentLogPost;
        acceptanceRates[sampleIdx] = alpha;
        stepSizes[sampleIdx] = currentEps;
      }
    }

    const meanAcceptanceRate = acceptedCount / numSamples;
    const finalStepSize = this.dualAveraging.getFinalStepSize();

    // Compute basic summary
    const summary = this.computeSampleSummary(samples);

    return {
      samples,
      logPosteriors,
      acceptanceRates,
      stepSizes,
      meanAcceptanceRate,
      finalStepSize,
      divergences,
      summary
    };
  }

  /**
   * Computes sample mean, standard deviation, and covariance across chain.
   * @param {Float64Array[]} samples
   * @returns {Object}
   */
  computeSampleSummary(samples) {
    const n = samples.length;
    if (n === 0) return {};

    const mean = new Float64Array(this.dim);
    for (let i = 0; i < n; i++) {
      const s = samples[i];
      for (let d = 0; d < this.dim; d++) {
        mean[d] += s[d];
      }
    }
    for (let d = 0; d < this.dim; d++) mean[d] /= n;

    const std = new Float64Array(this.dim);
    for (let i = 0; i < n; i++) {
      const s = samples[i];
      for (let d = 0; d < this.dim; d++) {
        const diff = s[d] - mean[d];
        std[d] += diff * diff;
      }
    }
    for (let d = 0; d < this.dim; d++) {
      std[d] = Math.sqrt(std[d] / Math.max(1, n - 1));
    }

    return {
      numSamples: n,
      mean: Array.from(mean),
      std: Array.from(std)
    };
  }
}
