/**
 * @file nuts_sampler.js
 * @description No-U-Turn Sampler (NUTS) with Dual-Averaging and Recursive Tree Building.
 *
 * Implements:
 * 1. Recursive binary tree trajectory builder with forward and backward leapfrog integration.
 * 2. General U-turn stopping criterion: (q^+ - q^-) . p^+ >= 0 and (q^+ - q^-) . p^- >= 0.
 * 3. Slice variable / Multinomial transition kernel (Betancourt 2017 / Hoffman & Gelman 2014).
 * 4. Dual-averaging step size adaptation targeting delta = 0.80.
 * 5. Diagnostic metrics: tree depth, energy errors, divergences, acceptance statistics.
 *
 * @module uncertainty/nuts_sampler
 */

import { LeapfrogIntegrator, computeKineticEnergy, leapfrogStep } from './leapfrog_integrator.js';
import { DualAveragingController } from './hmc_posterior_sampler.js';
import { MT19937 } from '../statistics/statistical_resampling.js';

/**
 * Checks the No-U-Turn condition between forward and backward trajectory endpoints.
 * @param {Float64Array} qMinus - Leftmost position
 * @param {Float64Array} qPlus - Rightmost position
 * @param {Float64Array} pMinus - Leftmost momentum
 * @param {Float64Array} pPlus - Rightmost momentum
 * @param {Float64Array} [invMass=null] - Inverse mass matrix
 * @returns {boolean} True if trajectory has NOT made a U-turn (keep expanding)
 */
export function checkNoUTurn(qMinus, qPlus, pMinus, pPlus, invMass = null) {
  const d = qMinus.length;
  let dotPlus = 0.0;
  let dotMinus = 0.0;

  for (let i = 0; i < d; i++) {
    const dq = qPlus[i] - qMinus[i];
    const vPlus = invMass ? pPlus[i] * invMass[i] : pPlus[i];
    const vMinus = invMass ? pMinus[i] * invMass[i] : pMinus[i];
    dotPlus += dq * vPlus;
    dotMinus += dq * vMinus;
  }

  return dotPlus >= 0.0 && dotMinus >= 0.0;
}

/**
 * The No-U-Turn Sampler (NUTS) implementation.
 */
export class NUTSSampler {
  /**
   * @param {Object} model
   * @param {function(Float64Array): number} model.logPosterior - Log posterior ln P(q)
   * @param {function(Float64Array): Float64Array} model.gradLogPosterior - Gradient nabla ln P(q)
   * @param {number} model.dimension - Dimension of parameter space
   * @param {Object} [options={}]
   * @param {number} [options.stepSize=0.1]
   * @param {number} [options.maxTreeDepth=10] - Maximum binary tree depth (2^10 = 1024 steps)
   * @param {number} [options.targetAcceptance=0.80] - Target acceptance probability
   * @param {number} [options.maxDeltaEnergy=1000.0] - Maximum energy difference before divergence
   * @param {number} [options.seed=42]
   */
  constructor(model, options = {}) {
    if (!model || typeof model.logPosterior !== \'function\' || typeof model.gradLogPosterior !== \'function\') {
      throw new TypeError(\'model must define logPosterior and gradLogPosterior functions\');
    }
    this.model = model;
    this.dim = model.dimension;

    this.stepSize = options.stepSize ?? 0.1;
    this.maxTreeDepth = options.maxTreeDepth ?? 10;
    this.targetAcceptance = options.targetAcceptance ?? 0.80;
    this.maxDeltaEnergy = options.maxDeltaEnergy ?? 1000.0;
    this.seed = options.seed ?? 42;

    this.prng = new MT19937(this.seed);
    this.dualAveraging = new DualAveragingController(this.targetAcceptance);

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
   * Samples Gaussian momentum vector from N(0, M).
   * @returns {Float64Array}
   */
  sampleMomentum() {
    const p = new Float64Array(this.dim);
    for (let i = 0; i < this.dim; i += 2) {
      const u1 = Math.max(1e-15, this.prng.extractNumber());
      const u2 = this.prng.extractNumber();
      const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
      const z1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(2.0 * Math.PI * u2);

      p[i] = z0 * Math.sqrt(this.mass[i]);
      if (i + 1 < this.dim) {
        p[i + 1] = z1 * Math.sqrt(this.mass[i + 1]);
      }
    }
    return p;
  }

  /**
   * Recursive binary tree builder for NUTS trajectory.
   *
   * @private
   */
  _buildTree(q, p, gradU, logU, direction, depth, epsilon, h0) {
    if (depth === 0) {
      // Base case: single leapfrog step in direction (+1 or -1)
      const stepEps = direction * epsilon;
      const stepRes = leapfrogStep(q, p, gradU, stepEps, this.gradPotentialFn, this.invMass);

      const h = stepRes.hamiltonian;
      const logPost = -stepRes.potential;
      const deltaH = h - h0;

      // Slice condition: logU <= -H(q, p)
      const nPrime = Number(logU <= -h);
      const sPrime = Number(logU < this.maxDeltaEnergy - h && Number.isFinite(h));
      const alpha = Math.min(1.0, Math.exp(-deltaH));

      return {
        qMinus: new Float64Array(stepRes.q),
        qPlus: new Float64Array(stepRes.q),
        pMinus: new Float64Array(stepRes.p),
        pPlus: new Float64Array(stepRes.p),
        gradUMinus: new Float64Array(stepRes.gradU),
        gradUPlus: new Float64Array(stepRes.gradU),
        qPrime: new Float64Array(stepRes.q),
        logPostPrime: logPost,
        nPrime,
        sPrime,
        alpha: Number.isFinite(alpha) ? alpha : 0.0,
        nAlpha: 1
      };
    }

    // Recursion: build left / first half of tree
    const firstHalf = this._buildTree(q, p, gradU, logU, direction, depth - 1, epsilon, h0);
    if (firstHalf.sPrime === 0) {
      return firstHalf;
    }

    // Build right / second half of tree extending from boundary
    let secondHalf;
    if (direction === -1) {
      secondHalf = this._buildTree(
        firstHalf.qMinus,
        firstHalf.pMinus,
        firstHalf.gradUMinus,
        logU,
        direction,
        depth - 1,
        epsilon,
        h0
      );
    } else {
      secondHalf = this._buildTree(
        firstHalf.qPlus,
        firstHalf.pPlus,
        firstHalf.gradUPlus,
        logU,
        direction,
        depth - 1,
        epsilon,
        h0
      );
    }

    // Metropolis proposal selection from second subtree
    let qPrime = firstHalf.qPrime;
    let logPostPrime = firstHalf.logPostPrime;
    const totalN = firstHalf.nPrime + secondHalf.nPrime;

    if (totalN > 0) {
      const probAcceptSecond = secondHalf.nPrime / totalN;
      if (this.prng.extractNumber() < probAcceptSecond) {
        qPrime = secondHalf.qPrime;
        logPostPrime = secondHalf.logPostPrime;
      }
    }

    // Determine boundary endpoints
    const qMinus = direction === -1 ? secondHalf.qMinus : firstHalf.qMinus;
    const qPlus = direction === 1 ? secondHalf.qPlus : firstHalf.qPlus;
    const pMinus = direction === -1 ? secondHalf.pMinus : firstHalf.pMinus;
    const pPlus = direction === 1 ? secondHalf.pPlus : firstHalf.pPlus;
    const gradUMinus = direction === -1 ? secondHalf.gradUMinus : firstHalf.gradUMinus;
    const gradUPlus = direction === 1 ? secondHalf.gradUPlus : firstHalf.gradUPlus;

    // Check U-turn condition
    const noUTurn = checkNoUTurn(qMinus, qPlus, pMinus, pPlus, this.invMass);
    const sPrime = Number(firstHalf.sPrime === 1 && secondHalf.sPrime === 1 && noUTurn);

    return {
      qMinus,
      qPlus,
      pMinus,
      pPlus,
      gradUMinus,
      gradUPlus,
      qPrime,
      logPostPrime,
      nPrime: totalN,
      sPrime,
      alpha: firstHalf.alpha + secondHalf.alpha,
      nAlpha: firstHalf.nAlpha + secondHalf.nAlpha
    };
  }

  /**
   * Executes NUTS posterior sampling for numSamples iterations after numWarmup steps.
   *
   * @param {Float64Array|number[]} initialPosition
   * @param {number} [numSamples=10000]
   * @param {number} [numWarmup=1000]
   * @returns {{
   *   samples: Float64Array[],
   *   logPosteriors: Float64Array,
   *   treeDepths: Uint16Array,
   *   acceptanceRates: Float64Array,
   *   stepSizes: Float64Array,
   *   divergences: number,
   *   meanTreeDepth: number,
   *   summary: Object
   * }}
   */
  sample(initialPosition, numSamples = 10000, numWarmup = 1000) {
    let currentQ = new Float64Array(initialPosition);
    let eval0 = this.gradPotentialFn(currentQ);
    let currentGradU = eval0.gradU;
    let currentLogPost = -eval0.potential;

    let eps = this.stepSize;
    this.dualAveraging.init(eps);

    const totalSteps = numWarmup + numSamples;
    const samples = [];
    const logPosteriors = new Float64Array(numSamples);
    const treeDepths = new Uint16Array(numSamples);
    const acceptanceRates = new Float64Array(numSamples);
    const stepSizes = new Float64Array(numSamples);

    let divergences = 0;
    let totalTreeDepth = 0;

    for (let step = 0; step < totalSteps; step++) {
      const isWarmup = step < numWarmup;
      const currentEps = isWarmup ? eps : this.dualAveraging.getFinalStepSize();

      // Sample momentum
      const p0 = this.sampleMomentum();
      const h0 = -currentLogPost + computeKineticEnergy(p0, this.invMass);

      // Sample slice variable u ~ Uniform(0, exp(-h0)) => log(u) = -h0 - Exp(1)
      const uExp = -Math.log(Math.max(1e-15, this.prng.extractNumber()));
      const logU = -h0 - uExp;

      let qMinus = new Float64Array(currentQ);
      let qPlus = new Float64Array(currentQ);
      let pMinus = new Float64Array(p0);
      let pPlus = new Float64Array(p0);
      let gradUMinus = new Float64Array(currentGradU);
      let gradUPlus = new Float64Array(currentGradU);

      let qProposed = new Float64Array(currentQ);
      let logPostProposed = currentLogPost;
      let depth = 0;
      let n = 1;
      let s = 1;

      let totalAlpha = 0;
      let totalNAlpha = 0;

      while (s === 1 && depth < this.maxTreeDepth) {
        // Choose expansion direction: -1 (backward) or +1 (forward)
        const direction = this.prng.extractNumber() < 0.5 ? -1 : 1;

        let treeRes;
        if (direction === -1) {
          treeRes = this._buildTree(qMinus, pMinus, gradUMinus, logU, direction, depth, currentEps, h0);
          qMinus = treeRes.qMinus;
          pMinus = treeRes.pMinus;
          gradUMinus = treeRes.gradUMinus;
        } else {
          treeRes = this._buildTree(qPlus, pPlus, gradUPlus, logU, direction, depth, currentEps, h0);
          qPlus = treeRes.qPlus;
          pPlus = treeRes.pPlus;
          gradUPlus = treeRes.gradUPlus;
        }

        if (treeRes.sPrime === 1) {
          const probAccept = treeRes.nPrime / Math.max(1, n + treeRes.nPrime);
          if (this.prng.extractNumber() < probAccept) {
            qProposed = treeRes.qPrime;
            logPostProposed = treeRes.logPostPrime;
          }
        }

        n += treeRes.nPrime;
        const noUTurn = checkNoUTurn(qMinus, qPlus, pMinus, pPlus, this.invMass);
        s = treeRes.sPrime === 1 && noUTurn ? 1 : 0;

        totalAlpha += treeRes.alpha;
        totalNAlpha += treeRes.nAlpha;
        depth++;
      }

      if (depth >= this.maxTreeDepth && s === 1) {
        // Reached max depth
      }
      if (s === 0 && depth === 1 && totalNAlpha === 1) {
        divergences++;
      }

      // Update current state
      currentQ = new Float64Array(qProposed);
      currentLogPost = logPostProposed;
      const evalNew = this.gradPotentialFn(currentQ);
      currentGradU = evalNew.gradU;

      const alphaAvg = totalNAlpha > 0 ? totalAlpha / totalNAlpha : 0.0;

      if (isWarmup) {
        eps = this.dualAveraging.update(alphaAvg);
      } else {
        const sampleIdx = step - numWarmup;
        samples.push(new Float64Array(currentQ));
        logPosteriors[sampleIdx] = currentLogPost;
        treeDepths[sampleIdx] = depth;
        acceptanceRates[sampleIdx] = alphaAvg;
        stepSizes[sampleIdx] = currentEps;
        totalTreeDepth += depth;
      }
    }

    const meanTreeDepth = totalTreeDepth / numSamples;
    const summary = this.computeSampleSummary(samples);

    return {
      samples,
      logPosteriors,
      treeDepths,
      acceptanceRates,
      stepSizes,
      divergences,
      meanTreeDepth,
      summary
    };
  }

  /**
   * @param {Float64Array[]} samples
   * @returns {Object}
   */
  computeSampleSummary(samples) {
    const n = samples.length;
    if (n === 0) return {};

    const mean = new Float64Array(this.dim);
    for (let i = 0; i < n; i++) {
      for (let d = 0; d < this.dim; d++) mean[d] += samples[i][d];
    }
    for (let d = 0; d < this.dim; d++) mean[d] /= n;

    const std = new Float64Array(this.dim);
    for (let i = 0; i < n; i++) {
      for (let d = 0; d < this.dim; d++) {
        const diff = samples[i][d] - mean[d];
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
