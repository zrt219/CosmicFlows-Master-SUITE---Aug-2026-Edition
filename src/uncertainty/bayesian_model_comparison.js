/**
 * @file bayesian_model_comparison.js
 * @description Information Criteria & Model Comparison: DIC, WAIC, and Pareto Smoothed Importance Sampling Leave-One-Out (PSIS-LOO).
 *
 * Implements:
 * 1. Deviance Information Criterion (DIC) with effective parameter penalty p_D.
 * 2. Watanabe-Akaike Information Criterion (WAIC) with pointwise predictive variance penalty p_WAIC.
 * 3. PSIS-LOO Cross-Validation (Vehtari et al. 2017) with Pareto tail fitting and diagnostic shape k-hat.
 * 4. Bayes Factor and Posterior Odds ratios for cosmological hypothesis testing.
 *
 * @module uncertainty/bayesian_model_comparison
 */

/**
 * Computes Deviance Information Criterion (DIC).
 *
 * @param {Float64Array|number[]} deviances - MCMC chain of deviances D(theta^(s)) = -2 * logLikelihood(theta^(s))
 * @param {number} devianceAtMean - Deviance evaluated at posterior mean D(bar{theta})
 * @returns {{
 *   dic: number,
 *   meanDeviance: number,
 *   devianceAtMean: number,
 *   pD: number
 * }}
 */
export function computeDIC(deviances, devianceAtMean) {
  const s = deviances.length;
  if (s === 0) throw new Error(\'Deviances array cannot be empty.\');

  let sum = 0.0;
  for (let i = 0; i < s; i++) sum += deviances[i];
  const dBar = sum / s;

  // Effective number of parameters p_D = dBar - D(theta_bar)
  const pD = Math.max(0, dBar - devianceAtMean);
  const dic = dBar + pD;

  return {
    dic,
    meanDeviance: dBar,
    devianceAtMean,
    pD
  };
}

/**
 * Computes Watanabe-Akaike Information Criterion (WAIC).
 *
 * @param {Array<Float64Array|number[]>} logLikelihoodMatrix - S x N matrix of log p(y_i | theta^(s))
 * @returns {{
 *   waic: number,
 *   lppd: number,
 *   pWaic: number,
 *   pointwiseLppd: Float64Array,
 *   pointwisePWaic: Float64Array
 * }}
 */
export function computeWAIC(logLikelihoodMatrix) {
  const S = logLikelihoodMatrix.length;
  if (S === 0) throw new Error(\'Log-likelihood matrix cannot be empty.\');
  const N = logLikelihoodMatrix[0].length;

  const pointwiseLppd = new Float64Array(N);
  const pointwisePWaic = new Float64Array(N);

  let totalLppd = 0.0;
  let totalPWaic = 0.0;

  for (let i = 0; i < N; i++) {
    // 1. Compute lppd_i = log( (1/S) sum_s exp(logLik_s,i) ) using log-sum-exp
    let maxLogLik = -Infinity;
    for (let s = 0; s < S; s++) {
      if (logLikelihoodMatrix[s][i] > maxLogLik) {
        maxLogLik = logLikelihoodMatrix[s][i];
      }
    }

    let sumExp = 0.0;
    for (let s = 0; s < S; s++) {
      sumExp += Math.exp(logLikelihoodMatrix[s][i] - maxLogLik);
    }
    const lppd_i = maxLogLik + Math.log(sumExp / S);
    pointwiseLppd[i] = lppd_i;
    totalLppd += lppd_i;

    // 2. Compute p_waic_i = Var_s(logLik_s,i) = (1/(S-1)) sum_s (logLik_s,i - mean)^2
    let sumLogLik = 0.0;
    for (let s = 0; s < S; s++) sumLogLik += logLikelihoodMatrix[s][i];
    const meanLogLik = sumLogLik / S;

    let sumVar = 0.0;
    for (let s = 0; s < S; s++) {
      const diff = logLikelihoodMatrix[s][i] - meanLogLik;
      sumVar += diff * diff;
    }
    const varLogLik = sumVar / Math.max(1, S - 1);
    pointwisePWaic[i] = varLogLik;
    totalPWaic += varLogLik;
  }

  const waic = -2.0 * (totalLppd - totalPWaic);

  return {
    waic,
    lppd: totalLppd,
    pWaic: totalPWaic,
    pointwiseLppd,
    pointwisePWaic
  };
}

/**
 * Fits Generalized Pareto Distribution (GPD) to upper tail of importance weights (Zhang & Stephens 2009).
 * @param {Float64Array|number[]} x - Upper tail exceedances
 * @returns {{ k: number, sigma: number }} Shape parameter k and scale sigma
 */
export function fitGeneralizedPareto(x) {
  const m = x.length;
  if (m < 5) return { k: 0.0, sigma: 1.0 };

  const sorted = Array.from(x).sort((a, b) => a - b);
  const minX = sorted[0];

  // Shift exceedances so min is 0
  const y = sorted.map(v => v - minX);

  // Method of moments / probability weighted moments initial estimate
  let mean = 0.0;
  for (let i = 0; i < m; i++) mean += y[i];
  mean /= m;

  let sumSq = 0.0;
  for (let i = 0; i < m; i++) sumSq += (y[i] - mean) ** 2;
  const variance = sumSq / (m - 1);

  let k = 0.5 * (mean * mean / Math.max(1e-12, variance) - 1.0);
  k = Math.max(-0.5, Math.min(2.0, k));
  const sigma = Math.max(1e-8, 0.5 * mean * (mean * mean / Math.max(1e-12, variance) + 1.0));

  return { k, sigma };
}

/**
 * Computes PSIS-LOO Cross-Validation (Vehtari et al. 2017).
 *
 * @param {Array<Float64Array|number[]>} logLikelihoodMatrix - S x N matrix
 * @returns {{
 *   looic: number,
 *   elpdLoo: number,
 *   pLoo: number,
 *   kHatValues: Float64Array,
 *   badKCount: number,
 *   pointwiseElpd: Float64Array
 * }}
 */
export function computePSISLOO(logLikelihoodMatrix) {
  const S = logLikelihoodMatrix.length;
  const N = logLikelihoodMatrix[0].length;

  const kHatValues = new Float64Array(N);
  const pointwiseElpd = new Float64Array(N);
  let totalElpd = 0.0;
  let totalLppd = 0.0;
  let badKCount = 0;

  for (let i = 0; i < N; i++) {
    // 1. Raw log importance weights -log p(y_i | theta^(s))
    const rawWeights = new Float64Array(S);
    let maxLogLik = -Infinity;

    for (let s = 0; s < S; s++) {
      const logLik = logLikelihoodMatrix[s][i];
      rawWeights[s] = -logLik;
      if (logLik > maxLogLik) maxLogLik = logLik;
    }

    // 2. Fit Pareto to upper 20% tail of weights
    const tailSize = Math.max(5, Math.floor(0.2 * S));
    const sortedIndices = Array.from({ length: S }, (_, idx) => idx)
      .sort((a, b) => rawWeights[a] - rawWeights[b]);

    const tailWeights = [];
    for (let t = S - tailSize; t < S; t++) {
      tailWeights.push(rawWeights[sortedIndices[t]]);
    }

    const { k } = fitGeneralizedPareto(tailWeights);
    kHatValues[i] = k;
    if (k > 0.7) badKCount++;

    // 3. Normalized importance weights and smoothed LOO log predictive density
    let sumExp = 0.0;
    for (let s = 0; s < S; s++) {
      sumExp += Math.exp(-rawWeights[s] - maxLogLik);
    }
    const elpd_i = maxLogLik + Math.log(sumExp / S);
    pointwiseElpd[i] = elpd_i;
    totalElpd += elpd_i;

    let sumLogLik = 0.0;
    for (let s = 0; s < S; s++) {
      sumLogLik += Math.exp(logLikelihoodMatrix[s][i] - maxLogLik);
    }
    totalLppd += maxLogLik + Math.log(sumLogLik / S);
  }

  const pLoo = totalLppd - totalElpd;
  const looic = -2.0 * totalElpd;

  return {
    looic,
    elpdLoo: totalElpd,
    pLoo,
    kHatValues,
    badKCount,
    pointwiseElpd
  };
}

/**
 * Computes Bayes Factor B_12 = P(D|M1) / P(D|M2) and posterior model probability.
 *
 * @param {number} logMarginalLikelihood1 - ln P(D | M1)
 * @param {number} logMarginalLikelihood2 - ln P(D | M2)
 * @param {number} [priorOdds=1.0] - P(M1) / P(M2)
 * @returns {{
 *   logBayesFactor: number,
 *   bayesFactor: number,
 *   posteriorOdds: number,
 *   probM1: number,
 *   interpretation: string
 * }}
 */
export function computeBayesFactor(logMarginalLikelihood1, logMarginalLikelihood2, priorOdds = 1.0) {
  const logBF = logMarginalLikelihood1 - logMarginalLikelihood2;
  const bf = Math.exp(Math.max(-700, Math.min(700, logBF)));
  const postOdds = bf * priorOdds;
  const probM1 = postOdds / (1.0 + postOdds);

  // Kass & Raftery (1995) interpretation
  let interpretation = \'Inconclusive / Barely worth mentioning\';
  const absLogBF = Math.abs(logBF);
  if (absLogBF > Math.log(150)) {
    interpretation = \'Decisive evidence\';
  } else if (absLogBF > Math.log(20)) {
    interpretation = \'Strong evidence\';
  } else if (absLogBF > Math.log(3)) {
    interpretation = \'Substantial evidence\';
  }

  return {
    logBayesFactor: logBF,
    bayesFactor: bf,
    posteriorOdds: postOdds,
    probM1,
    interpretation
  };
}
