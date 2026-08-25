/**
 * @file watershed_comparator.js
 * @description Quantitative Segmentation Comparator for Cosmicflows Watershed Basins.
 * Computes Jaccard index, Dice similarity coefficient, Precision, Recall, Centroid displacement,
 * Volume ratio, Transition matrix, and Cohen's Kappa between published reference and derived basins.
 *
 * References:
 * - Dupuy & Courtois (2023) MNRAS 527: Stability and Transition Matrix of Cosmic Basins.
 */

import { VoxelGrid, getBasinById, computeBasinStatistics } from './watershed_classifier.js';

/**
 * @typedef {Object} BasinPairMetrics
 * @property {number} basinId - Numeric basin ID
 * @property {string} name - Basin name
 * @property {number} countRef - Reference voxel count |A|
 * @property {number} countComp - Compared voxel count |B|
 * @property {number} countOverlap - Overlap voxel count |A ∩ B|
 * @property {number} countUnion - Union voxel count |A ∪ B|
 * @property {number} jaccardIndex - Jaccard similarity coefficient |A ∩ B| / |A ∪ B|
 * @property {number} diceCoefficient - Dice similarity 2|A ∩ B| / (|A| + |B|)
 * @property {number} precision - Precision |A ∩ B| / |B|
 * @property {number} recall - Recall |A ∩ B| / |A|
 * @property {number} f1Score - Harmonic mean of precision and recall
 * @property {number} volumeRefMpc3 - Volume in Mpc/h^3
 * @property {number} volumeCompMpc3 - Volume in Mpc/h^3
 * @property {number} volumeRatio - Volume ratio V_comp / V_ref
 * @property {number} centroidDisplacementMpc - Spatial centroid displacement in Mpc/h
 * @property {number} centroidDisplacementKms - Spatial centroid displacement in km/s
 * @property {[number, number, number]} centroidRefMpc - Reference centroid [x, y, z]
 * @property {[number, number, number]} centroidCompMpc - Compared centroid [x, y, z]
 */

/**
 * @typedef {Object} ComparisonResult
 * @property {Record<number, BasinPairMetrics>} perBasinMetrics - Metrics for each basin ID
 * @property {number[][]} confusionMatrix - Full N_basins x N_basins intersection counts C[i][j]
 * @property {number[][]} transitionMatrixPct - Row-normalized transition percentages P(B_j | A_i)
 * @property {number[]} basinIds - Ordered list of evaluated basin IDs
 * @property {number} macroJaccard - Unweighted mean of Jaccard indices across active basins (excluding background 0)
 * @property {number} macroDice - Unweighted mean of Dice coefficients across active basins
 * @property {number} microJaccard - Global volume-weighted Jaccard index
 * @property {number} microDice - Global volume-weighted Dice coefficient
 * @property {number} overallAccuracy - Proportion of identically labeled voxels
 * @property {number} cohensKappa - Cohen's Kappa inter-rater agreement statistic
 */

/**
 * Computes comprehensive quantitative segmentation metrics comparing two VoxelGrid segmentations.
 *
 * @param {VoxelGrid} gridRef - Reference / Ground Truth segmentation (e.g. Individual CF4)
 * @param {VoxelGrid} gridComp - Compared / Test segmentation (e.g. Grouped CF4 or Derived)
 * @returns {ComparisonResult}
 */
export function compareSegmentations(gridRef, gridComp) {
  if (gridRef.totalVoxels !== gridComp.totalVoxels) {
    throw new Error(`Dimension mismatch: gridRef has ${gridRef.totalVoxels} voxels, gridComp has ${gridComp.totalVoxels}`);
  }

  const total = gridRef.totalVoxels;
  const dxMpc = gridRef.metadata.dxMpc;
  const voxelVolMpc3 = dxMpc * dxMpc * dxMpc;
  const h0 = gridRef.metadata.hubbleH0;

  // Step 1: Collect unique labels across both grids
  const uniqueLabelsSet = new Set();
  for (let i = 0; i < total; i++) {
    uniqueLabelsSet.add(gridRef.data[i]);
    uniqueLabelsSet.add(gridComp.data[i]);
  }
  const basinIds = Array.from(uniqueLabelsSet).sort((a, b) => a - b);
  const maxLabel = Math.max(...basinIds);

  // Initialize confusion matrix
  const confusion = Array.from({ length: maxLabel + 1 }, () => new Float64Array(maxLabel + 1));
  const countsRef = new Float64Array(maxLabel + 1);
  const countsComp = new Float64Array(maxLabel + 1);

  // Compute spatial moments for centroid displacement
  const momentsRef = {};
  const momentsComp = {};
  for (const id of basinIds) {
    momentsRef[id] = { sumX: 0, sumY: 0, sumZ: 0, count: 0 };
    momentsComp[id] = { sumX: 0, sumY: 0, sumZ: 0, count: 0 };
  }

  let identicalVoxels = 0;

  for (let idx = 0; idx < total; idx++) {
    const lRef = gridRef.data[idx];
    const lComp = gridComp.data[idx];

    confusion[lRef][lComp]++;
    countsRef[lRef]++;
    countsComp[lComp]++;

    if (lRef === lComp) {
      identicalVoxels++;
    }

    const [ix, iy, iz] = gridRef.getGridCoords(idx);
    const { x, y, z } = gridRef.gridToWorld(ix, iy, iz, 'Mpc');

    if (momentsRef[lRef]) {
      momentsRef[lRef].sumX += x;
      momentsRef[lRef].sumY += y;
      momentsRef[lRef].sumZ += z;
      momentsRef[lRef].count++;
    }

    if (momentsComp[lComp]) {
      momentsComp[lComp].sumX += x;
      momentsComp[lComp].sumY += y;
      momentsComp[lComp].sumZ += z;
      momentsComp[lComp].count++;
    }
  }

  // Step 2: Per-basin metrics
  const perBasinMetrics = {};
  let sumJaccard = 0;
  let sumDice = 0;
  let activeBasinCount = 0;
  let globalIntersection = 0;
  let globalUnion = 0;

  for (const id of basinIds) {
    const a = countsRef[id];
    const b = countsComp[id];
    const inter = confusion[id][id];
    const un = a + b - inter;

    const jaccard = un > 0 ? inter / un : (a === 0 && b === 0 ? 1.0 : 0.0);
    const dice = (a + b) > 0 ? (2.0 * inter) / (a + b) : (a === 0 && b === 0 ? 1.0 : 0.0);
    const precision = b > 0 ? inter / b : (inter === 0 && a === 0 ? 1.0 : 0.0);
    const recall = a > 0 ? inter / a : (inter === 0 && b === 0 ? 1.0 : 0.0);
    const f1Score = (precision + recall) > 0 ? (2 * precision * recall) / (precision + recall) : 0.0;

    const volRef = a * voxelVolMpc3;
    const volComp = b * voxelVolMpc3;
    const volRatio = volRef > 0 ? volComp / volRef : (volComp === 0 ? 1.0 : Infinity);

    // Centroids
    const mR = momentsRef[id];
    const mC = momentsComp[id];

    const cRef = mR && mR.count > 0 ? [mR.sumX / mR.count, mR.sumY / mR.count, mR.sumZ / mR.count] : [0, 0, 0];
    const cComp = mC && mC.count > 0 ? [mC.sumX / mC.count, mC.sumY / mC.count, mC.sumZ / mC.count] : [0, 0, 0];

    const dx = cRef[0] - cComp[0];
    const dy = cRef[1] - cComp[1];
    const dz = cRef[2] - cComp[2];
    const dispMpc = (mR && mR.count > 0 && mC && mC.count > 0) ? Math.hypot(dx, dy, dz) : 0;
    const dispKms = dispMpc * h0;

    const basinDef = getBasinById(id);

    perBasinMetrics[id] = {
      basinId: id,
      name: basinDef.name,
      countRef: a,
      countComp: b,
      countOverlap: inter,
      countUnion: un,
      jaccardIndex: jaccard,
      diceCoefficient: dice,
      precision,
      recall,
      f1Score,
      volumeRefMpc3: volRef,
      volumeCompMpc3: volComp,
      volumeRatio: volRatio,
      centroidDisplacementMpc: dispMpc,
      centroidDisplacementKms: dispKms,
      centroidRefMpc: cRef,
      centroidCompMpc: cComp
    };

    if (id !== 0 && (a > 0 || b > 0)) {
      sumJaccard += jaccard;
      sumDice += dice;
      activeBasinCount++;
      globalIntersection += inter;
      globalUnion += un;
    }
  }

  // Step 3: Transition Matrix (% row probabilities)
  const transitionMatrixPct = [];
  for (let i = 0; i <= maxLabel; i++) {
    const row = [];
    const rowTotal = countsRef[i];
    for (let j = 0; j <= maxLabel; j++) {
      row.push(rowTotal > 0 ? (confusion[i][j] / rowTotal) * 100.0 : 0.0);
    }
    transitionMatrixPct.push(row);
  }

  // Step 4: Cohen's Kappa & Aggregate Accuracy
  const overallAccuracy = identicalVoxels / total;

  let expectedAgreement = 0;
  for (const id of basinIds) {
    expectedAgreement += (countsRef[id] / total) * (countsComp[id] / total);
  }
  const cohensKappa = (1.0 - expectedAgreement) > 1e-12
    ? (overallAccuracy - expectedAgreement) / (1.0 - expectedAgreement)
    : 1.0;

  const macroJaccard = activeBasinCount > 0 ? sumJaccard / activeBasinCount : 1.0;
  const macroDice = activeBasinCount > 0 ? sumDice / activeBasinCount : 1.0;
  const microJaccard = globalUnion > 0 ? globalIntersection / globalUnion : 1.0;
  const microDice = (globalIntersection * 2) / (globalUnion + globalIntersection || 1);

  // Convert confusion to regular arrays
  const confusionArray = confusion.map(row => Array.from(row));

  return {
    perBasinMetrics,
    confusionMatrix: confusionArray,
    transitionMatrixPct,
    basinIds,
    macroJaccard,
    macroDice,
    microJaccard,
    microDice,
    overallAccuracy,
    cohensKappa
  };
}
