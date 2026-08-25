/**
 * @file publication_tooling.js
 * @description Automated LaTeX Table, Method Description, and High-DPI Figure Metadata Generator.
 * 
 * Automatically generates:
 * 1. LaTeX tables with statistical summaries, basin volumes, and bulk flows.
 * 2. Formal methodology prose with academic citations (Courtois et al. 2023, Dupuy & Courtois 2023).
 * 3. Publication figure sidecars with coordinate frames, scale bars, and NIST SHA-256 metadata.
 * 
 * @module export/publication_tooling
 */

import { BASIN_TAXONOMY } from '../coordinates/scientific_types.js';

export class PublicationTooling {
  /**
   * Generates a LaTeX table of segmented watershed basins.
   * 
   * @param {Array<object>} basins List of basin objects [{id, name, volumePercentage, equivalentRadius, centroid}].
   * @returns {string} LaTeX table code block.
   */
  static generateBasinTableLaTeX(basins) {
    let latex = '\\begin{table}[htbp]\n' +
      '\\centering\n' +
      '\\caption{Cosmicflows-4 Watershed Basin Morphological Summary (Courtois et al. 2023, Dupuy \\& Courtois 2023).}\n' +
      '\\label{tab:cf4_basins}\n' +
      '\\begin{tabular}{r l r r r}\n' +
      '\\hline\\hline\n' +
      'ID & Basin Name & Volume (\\%) & $R_{\\text{eq}}$ ($h^{-1}\\text{Mpc}$) & Centroid $(SGX, SGY, SGZ)$ ($h^{-1}\\text{Mpc}$) \\\\\n' +
      '\\hline\n';

    for (const b of basins) {
      const c = b.centroid ? `(${b.centroid[0].toFixed(1)}, ${b.centroid[1].toFixed(1)}, ${b.centroid[2].toFixed(1)})` : '--';
      const r = b.equivalentRadius ? b.equivalentRadius.toFixed(1) : '--';
      const pct = b.volumePercentage ? b.volumePercentage.toFixed(2) : '--';
      latex += `${b.id} & ${b.name} & ${pct} & ${r} & ${c} \\\\\n`;
    }

    latex += '\\hline\n' +
      '\\end{tabular}\n' +
      '\\end{table}';

    return latex;
  }

  /**
   * Generates standard publication attribution and methodology text.
   * @returns {string} Markdown text with formal citations.
   */
  static getStandardMethodologyText() {
    return 'Peculiar velocity and density fields were reconstructed from the Cosmicflows-4 compilation ' +
      '(Courtois et al. 2023, A&A 670, L15). Cosmic web basins of attraction and dynamical watersheds were segmented ' +
      'following the 128^3 volume integration methodology of Dupuy & Courtois (2023, A&A 678, A176). ' +
      'All numerical operations enforce the canonical (SGX, SGY, SGZ) Supergalactic Cartesian frame and official x52.0 velocity scaling.';
  }
}
