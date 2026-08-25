/**
 * @file watershed_segmentation.js
 * @description High-level Watershed Segmentation and Table A.1 Basin Manager.
 * 
 * @module watershed/watershed_segmentation
 */

import { BASIN_TAXONOMY, EXTENDED_BASIN_TAXONOMY } from './watershed_classifier.js';

export class WatershedSegmentation {
  constructor(options = {}) {
    this.reconstruction = options.reconstruction || 'cf4-ungrouped';
  }

  /**
   * Retrieves Table A.1 metadata for a basin ID.
   * @param {number} basinId - Integer basin identifier.
   * @param {object} [opts={}] - Optional reconstruction settings.
   * @returns {object} Basin metadata record.
   */
  getBasinMetadata(basinId, opts = {}) {
    const basin = (EXTENDED_BASIN_TAXONOMY && EXTENDED_BASIN_TAXONOMY[basinId]) || BASIN_TAXONOMY[basinId];
    if (!basin) {
      return {
        id: basinId,
        name: 'Basin-' + basinId,
        abbreviation: 'B' + basinId,
        description: 'Uncatalogued basin'
      };
    }
    return {
      id: basin.id,
      name: basin.name,
      abbreviation: basin.abbreviation,
      description: basin.description,
      attractorSGMpc: basin.attractorSGMpc,
      volume: opts.reconstruction === 'cf4-grouped' ? basin.volumeCF4Grp1e6 : basin.volumeCF4Ind1e6,
      hexColor: basin.hexColor
    };
  }
}
