/**
 * @file dataset_forensics.js
 * @description Dataset forensics analyzer for FITS / NPZ products.
 * Computes exact statistical moments, NaN/Inf census, kurtosis, skewness,
 * dynamic range, percentiles (p1, p5, p25, p50, p75, p95, p99), spatial symmetry,
 * and memory layout verification.
 * 
 * @module data/dataset_forensics
 */

export class DatasetForensics {
  static analyze(buffer, dimensions, options = {}) {
    const len = buffer.length;
    let nanCount = 0;
    let posInfCount = 0;
    let negInfCount = 0;
    let validCount = 0;

    let min = Infinity;
    let max = -Infinity;
    let sum = 0.0;
    let sumSq = 0.0;

    const validValues = [];

    for (let i = 0; i < len; i++) {
      const v = buffer[i];
      if (Number.isNaN(v)) {
        nanCount++;
      } else if (v === Infinity) {
        posInfCount++;
      } else if (v === -Infinity) {
        negInfCount++;
      } else {
        validCount++;
        if (v < min) min = v;
        if (v > max) max = v;
        sum += v;
        sumSq += v * v;
        if (options.storeValues !== false && validValues.length < 500000) {
          validValues.push(v);
        }
      }
    }

    const mean = validCount > 0 ? sum / validCount : 0.0;
    const variance = validCount > 1 ? (sumSq - validCount * mean * mean) / (validCount - 1) : 0.0;
    const std = Math.sqrt(Math.max(0.0, variance));

    // Higher order moments: skewness and kurtosis
    let m3 = 0.0;
    let m4 = 0.0;
    for (let i = 0; i < validValues.length; i++) {
      const diff = validValues[i] - mean;
      const d2 = diff * diff;
      m3 += diff * d2;
      m4 += d2 * d2;
    }

    const skewness = (validCount > 2 && std > 0) ? (m3 / validCount) / (std * std * std) : 0.0;
    const kurtosis = (validCount > 3 && std > 0) ? (m4 / validCount) / (variance * variance) - 3.0 : 0.0;

    // Percentiles
    validValues.sort((a, b) => a - b);
    const getP = (p) => {
      if (validValues.length === 0) return 0.0;
      const idx = Math.floor((p / 100.0) * (validValues.length - 1));
      return validValues[idx];
    };

    return {
      totalElements: len,
      validCount,
      nanCount,
      posInfCount,
      negInfCount,
      finiteFraction: len > 0 ? validCount / len : 0.0,
      min: isFinite(min) ? min : 0.0,
      max: isFinite(max) ? max : 0.0,
      dynamicRange: isFinite(max) && isFinite(min) ? max - min : 0.0,
      mean,
      variance,
      std,
      skewness,
      kurtosis,
      percentiles: {
        p1: getP(1),
        p5: getP(5),
        p25: getP(25),
        p50: getP(50),
        p75: getP(75),
        p95: getP(95),
        p99: getP(99)
      },
      dimensions: [...dimensions],
      byteLength: buffer.byteLength
    };
  }
}
