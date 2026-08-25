/**
 * @file adaptive_step_controller.js
 * @description Advanced adaptive step size controller for high-order ODE integrators.
 * Supports:
 *   - Standard I-controller (Integral)
 *   - Gustafsson & Hairer-Söderlind PI-controller (Proportional-Integral)
 *   - Hairer-Söderlind H211b PID-controller
 *   - Error norm evaluation with mixed absolute/relative component tolerances
 *   - Step doubling / halving safeguards and stiffness monitoring.
 *
 * @module integration/adaptive_step_controller
 */

export const StepControllerType = Object.freeze({
  INTEGRAL: 'I',
  PI: 'PI',
  PID: 'PID'
});

export function computeScaledErrorNorm(errVec, y0, y1, options = {}) {
  const atol = options.atol ?? 1e-6;
  const rtol = options.rtol ?? 1e-6;
  const d = errVec.length;

  let sumSq = 0.0;
  for (let i = 0; i < d; i++) {
    const atol_i = Array.isArray(atol) ? atol[i] : atol;
    const rtol_i = Array.isArray(rtol) ? rtol[i] : rtol;
    const yMax = Math.max(Math.abs(y0[i]), Math.abs(y1[i]));
    const sc_i = atol_i + yMax * rtol_i;
    const ratio = errVec[i] / sc_i;
    sumSq += ratio * ratio;
  }

  return Math.sqrt(sumSq / d);
}

export class AdaptiveStepController {
  constructor(config = {}) {
    this.order = config.order ?? 5;
    this.type = config.type ?? StepControllerType.PI;
    this.safety = config.safety ?? 0.9;
    this.facMin = config.facMin ?? 0.2;
    this.facMax = config.facMax ?? 5.0;
    this.minStep = config.minStep ?? 1e-8;
    this.maxStep = config.maxStep ?? 100.0;
    this.atol = config.atol ?? 1e-6;
    this.rtol = config.rtol ?? 1e-6;
    this.maxConsecutiveRejections = config.maxConsecutiveRejections ?? 10;

    const p = this.order;
    if (this.type === StepControllerType.PI) {
      this.kI = config.kI ?? (0.6 / p);
      this.kP = config.kP ?? (0.2 / p);
      this.kD = 0.0;
    } else if (this.type === StepControllerType.PID) {
      this.kI = config.kI ?? (0.7 / p);
      this.kP = config.kP ?? (0.4 / p);
      this.kD = config.kD ?? (0.1 / p);
    } else {
      this.kI = config.kI ?? (1.0 / p);
      this.kP = 0.0;
      this.kD = 0.0;
    }

    this.prevErrNorm = 1.0;
    this.prevPrevErrNorm = 1.0;
    this.prevStep = 0.5;
    this.consecutiveRejections = 0;
    this.acceptedSteps = 0;
    this.rejectedSteps = 0;
  }

  reset() {
    this.prevErrNorm = 1.0;
    this.prevPrevErrNorm = 1.0;
    this.prevStep = 0.5;
    this.consecutiveRejections = 0;
    this.acceptedSteps = 0;
    this.rejectedSteps = 0;
  }

  evaluateStep(errVec, y0, y1, currentStep) {
    const errNorm = computeScaledErrorNorm(errVec, y0, y1, { atol: this.atol, rtol: this.rtol });
    const accepted = errNorm <= 1.0;

    let factor;
    const safeErr = Math.max(errNorm, 1e-10);

    if (this.type === StepControllerType.PI) {
      factor = this.safety * Math.pow(1.0 / safeErr, this.kI) * Math.pow(this.prevErrNorm, this.kP);
    } else if (this.type === StepControllerType.PID) {
      factor = this.safety * Math.pow(1.0 / safeErr, this.kI) * Math.pow(this.prevErrNorm, this.kP) * Math.pow(1.0 / Math.max(this.prevPrevErrNorm, 1e-10), this.kD);
    } else {
      factor = this.safety * Math.pow(1.0 / safeErr, 1.0 / this.order);
    }

    factor = Math.max(this.facMin, Math.min(this.facMax, factor));

    if (accepted) {
      this.acceptedSteps++;
      this.consecutiveRejections = 0;
      this.prevPrevErrNorm = this.prevErrNorm;
      this.prevErrNorm = safeErr;
      this.prevStep = currentStep;

      let nextStep = currentStep * factor;
      nextStep = Math.max(this.minStep, Math.min(this.maxStep, nextStep));

      return {
        accepted: true,
        errNorm,
        nextStep,
        factor,
        retry: false
      };
    } else {
      this.rejectedSteps++;
      this.consecutiveRejections++;

      if (this.consecutiveRejections >= this.maxConsecutiveRejections) {
        return {
          accepted: false,
          errNorm,
          nextStep: Math.max(this.minStep, currentStep * 0.5),
          factor: 0.5,
          retry: false,
          aborted: true
        };
      }

      let retryStep = currentStep * Math.min(factor, 0.5);
      retryStep = Math.max(this.minStep, Math.min(this.maxStep, retryStep));

      return {
        accepted: false,
        errNorm,
        nextStep: retryStep,
        factor,
        retry: true
      };
    }
  }
}

export default AdaptiveStepController;
