/**
 * @file dopri5_integrator.js
 * @description Canonical export wrapper and alias module for Dormand-Prince 5(4) (DOPRI5) Integrator.
 * Provides uniform .step(), .trace(), and .integrate() interfaces alongside Hairer-Wanner
 * 4th-order dense continuous interpolation and FSAL stage recycling.
 *
 * @module integration/dopri5_integrator
 */

import {
  DOPRI5_COEFFS,
  dormandPrinceStepRaw,
  dormandPrinceStep,
  traceDormandPrinceStreamline,
  denseInterpolateDOPRI5,
  DormandPrinceIntegrator
} from './dormand_prince.js';

export {
  DOPRI5_COEFFS,
  dormandPrinceStepRaw,
  dormandPrinceStep,
  traceDormandPrinceStreamline,
  denseInterpolateDOPRI5,
  DormandPrinceIntegrator
};

export class DOPRI5Integrator extends DormandPrinceIntegrator {
  constructor(config = {}) {
    super(config);
  }

  integrate(fieldFn, seedPos, overrideOptions = {}) {
    return this.trace(fieldFn, seedPos, overrideOptions);
  }
}

export default DOPRI5Integrator;
