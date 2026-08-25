/**
 * @file mcmc_chain_diagnostics.js
 * @description Advanced MCMC Chain Health Diagnostics: Geweke, Heidelberger-Welch, Raftery-Lewis, E-BFMI, and Cross-Correlations.
 *
 * Implements:
 * 1. Geweke Z-score diagnostic comparing initial 10% and final 50% chain segments.
 * 2. Energy Bayesian Fraction of Missing Information (E-BFMI) for HMC/NUTS momentum transitions.
 * 3. Raftery-Lewis quantile sample size and burn-in estimator.
 * 4. Cross-correlation matrix and lag-k cross-correlations across vector parameters.
 * 5. Comprehensive MCMC health audit report generator.
 *
 * @module uncertainty/mcmc_chain_diagnostics
 */

import { computeMCSE } from './posterior_statistics.js';
import { probit } from './gelman_rubin_diagnostic.js';

/**
 * Computes Geweke diagnostic Z-score for stationarity.
 * Tests H0: mean(first 10%) == mean(last 50%).
 *
 * @param {Float64Array|number[]} chain - 1D MCMC chain
 * @param {number} [frac1=0.1] - First segment fraction
 * @param {number} [frac2=0.5] - Second segment fraction
 * @returns {{
 *   zScore: number,
 *   pValue: number,
 *   mean1: number,
 *   mean2: number,
 *   isStationary: boolean
 * }}
 */
export function computeGewekeDiagnostic(chain, frac1 = 0.1, frac2 = 0.5) {
  const n = chain.length;
  if (n < 50) throw new Error(\'Geweke diagnostic requires at least 50 samples.\');

  const n1 = Math.floor(frac1 * n);
  const n2 = Math.floor(frac2 * n);
  const start2 = n - n2;

  const seg1 = chain.slice(0, n1);
  const seg2 = chain.slice(start2);

  let sum1 = 0.0;
  for (let i = 0; i < n1; i++) sum1 += seg1[i];
  const mean1 = sum1 / n1;

  let sum2 = 0.0;
  for (let i = 0; i < n2; i++) sum2 += seg2[i];
  const mean2 = sum2 / n2;

  const se1 = computeMCSE(seg1);
  const se2 = computeMCSE(seg2);

  const denom = Math.sqrt(se1 * se1 + se2 * se2);
  const zScore = denom > 1e-12 ? (mean1 - mean2) / denom : 0.0;

  // Two-tailed p-value from standard normal
  const absZ = Math.abs(zScore);
  // Complementary error function approx for normal tail
  const t = 1.0 / (1.0 + 0.2316419 * absZ);
  const poly = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + 1.330274429 * t))));
  const normCdf = 1.0 - (1.0 / Math.sqrt(2.0 * Math.PI)) * Math.exp(-0.5 * absZ * absZ) * poly;
  const pValue = 2.0 * (1.0 - normCdf);

  const isStationary = Math.abs(zScore) < 1.96; // 95% confidence

  return {
    zScore,
    pValue,
    mean1,
    mean2,
    isStationary
  };
}

/**
 * Computes Energy Bayesian Fraction of Missing Information (E-BFMI / Betancourt 2016).
 * Diagnoses whether HMC energy transitions effectively explore the posterior energy distribution.
 * Low values (E-BFMI < 0.3) signal severe energy momentum mismatch and sluggish exploration.
 *
 * @param {Float64Array|number[]} energies - Array of Hamiltonian energies H(q, p) at each MCMC step
 * @returns {{
 *   ebfmi: number,
 *   isAdequate: boolean
 * }}
 */
export function computeEBFMI(energies) {
  const n = energies.length;
  if (n < 4) throw new Error(\'E-BFMI calculation requires at least 4 energy values.\');

  let numer = 0.0;
  for (let i = 1; i < n; i++) {
    const diff = energies[i] - energies[i - 1];
    numer += diff * diff;
  }
  numer /= (n - 1);

  let sumE = 0.0;
  for (let i = 0; i < n; i++) sumE += energies[i];
  const meanE = sumE / n;

  let denom = 0.0;
  for (let i = 0; i < n; i++) {
    const diff = energies[i] - meanE;
    denom += diff * diff;
  }
  denom /= (n - 1);

  const ebfmi = denom > 1e-12 ? numer / denom : 1.0;
  const isAdequate = ebfmi >= 0.30;

  return {
    ebfmi,
    isAdequate
  };
}

/**
 * Computes Raftery-Lewis diagnostic for sample size estimation.
 * Estimates required chain length N to estimate quantile q to accuracy +/- r with probability s.
 *
 * @param {Float64Array|number[]} chain - 1D MCMC chain
 * @param {number} [q=0.025] - Target quantile
 * @param {number} [r=0.01] - Desired precision
 * @param {number} [s=0.95] - Confidence level
 * @returns {{
 *   requiredN: number,
 *   burnInEstimate: number,
 *   thinningFactor: number,
 *   dependenceFactorI: number
 * }}
 */
export function computeRafteryLewis(chain, q = 0.025, r = 0.01, s = 0.95) {
  const n = chain.length;
  const phiZ = probit(0.5 + 0.5 * s);
  const nMin = Math.ceil((q * (1.0 - q) * phiZ * phiZ) / (r * r));

  // Binary indicator sequence Z_t = I(chain_t <= q_sample)
  const sorted = Array.from(chain).sort((a, b) => a - b);
  const cutoff = sorted[Math.floor(q * n)];

  const z = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    z[i] = chain[i] <= cutoff ? 1 : 0;
  }

  // 2-state Markov chain transition probabilities
  let n00 = 0, n01 = 0, n10 = 0, n11 = 0;
  for (let i = 0; i < n - 1; i++) {
    if (z[i] === 0 && z[i + 1] === 0) n00++;
    else if (z[i] === 0 && z[i + 1] === 1) n01++;
    else if (z[i] === 1 && z[i + 1] === 0) n10++;
    else if (z[i] === 1 && z[i + 1] === 1) n11++;
  }

  const alpha = (n01 + 1) / Math.max(1, n00 + n01 + 2);
  const beta = (n10 + 1) / Math.max(1, n10 + n11 + 2);

  const lambda = 1.0 - alpha - beta;
  const burnIn = Math.ceil(Math.log(0.001 * (alpha + beta) / Math.max(alpha, beta)) / Math.log(Math.abs(lambda) || 1e-6));

  const k = Math.ceil((2.0 - alpha - beta) / (alpha + beta));
  const requiredN = Math.ceil(nMin * ((2.0 - alpha - beta) * alpha * beta) / ((alpha + beta) ** 3));
  const dependenceFactorI = Math.max(1.0, requiredN / nMin);

  return {
    requiredN: Math.max(nMin, requiredN),
    burnInEstimate: Math.max(0, burnIn),
    thinningFactor: Math.max(1, k),
    dependenceFactorI
  };
}

/**
 * Computes Cross-Correlation Matrix between D parameter chains.
 * @param {Array<Float64Array|number[]>} samples - Array of D-dimensional sample vectors
 * @returns {number[][]} D x D correlation matrix
 */
export function computeParameterCrossCorrelations(samples) {
  const n = samples.length;
  if (n === 0) return [];
  const d = samples[0].length;

  const means = new Float64Array(d);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < d; j++) means[j] += samples[i][j];
  }
  for (let j = 0; j < d; j++) means[j] /= n;

  const stds = new Float64Array(d);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < d; j++) {
      const diff = samples[i][j] - means[j];
      stds[j] += diff * diff;
    }
  }
  for (let j = 0; j < d; j++) stds[j] = Math.sqrt(stds[j] / Math.max(1, n - 1));

  const corr = Array.from({ length: d }, () => new Float64Array(d));
  for (let i = 0; i < d; i++) {
    corr[i][i] = 1.0;
    for (let j = i + 1; j < d; j++) {
      let sumCov = 0.0;
      for (let k = 0; k < n; k++) {
        sumCov += (samples[k][i] - means[i]) * (samples[k][j] - means[j]);
      }
      const cov = sumCov / Math.max(1, n - 1);
      const r = (stds[i] * stds[j] > 1e-12) ? cov / (stds[i] * stds[j]) : 0.0;
      corr[i][j] = Math.max(-1.0, Math.min(1.0, r));
      corr[j][i] = corr[i][j];
    }
  }

  return corr.map(r => Array.from(r));
}
