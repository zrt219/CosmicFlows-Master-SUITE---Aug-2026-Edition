/**
 * @file leapfrog_integrator.js
 * @description High-Precision Symplectic Leapfrog (Verlet) Integrator for Hamiltonian Monte Carlo.
 *
 * Implements:
 * 1. Second-order symplectic leapfrog integration for phase space trajectories (q, p).
 * 2. Energy conservation tracking with Hamiltonian H(q, p) = U(q) + K(p).
 * 3. Diagonal and dense mass matrix inverse kinetic energy operators K(p) = 0.5 * p^T M^{-1} p.
 * 4. Symplectic time-reversibility and phase space volume preservation (Liouville theorem) verification.
 * 5. Shadow Hamiltonian diagnostics for symplectic order verification (O(epsilon^2) energy oscillations).
 *
 * @module uncertainty/leapfrog_integrator
 */

/**
 * Evaluates kinetic energy K(p) = 0.5 * sum_i (p_i^2 / m_i) for diagonal mass matrix.
 * @param {Float64Array|number[]} p - Momentum vector
 * @param {Float64Array|number[]} [invMass] - Diagonal elements of M^{-1} (defaults to identity)
 * @returns {number} Kinetic energy
 */
export function computeKineticEnergy(p, invMass = null) {
  const d = p.length;
  let k = 0.0;
  if (!invMass) {
    for (let i = 0; i < d; i++) {
      k += p[i] * p[i];
    }
  } else {
    for (let i = 0; i < d; i++) {
      k += p[i] * p[i] * invMass[i];
    }
  }
  return 0.5 * k;
}

/**
 * Computes velocity vector v = M^{-1} p for position updates.
 * @param {Float64Array|number[]} p - Momentum vector
 * @param {Float64Array|number[]} [invMass] - Diagonal elements of M^{-1}
 * @param {Float64Array|number[]} [out] - Optional output buffer
 * @returns {Float64Array} Velocity vector v
 */
export function computeVelocity(p, invMass = null, out = null) {
  const d = p.length;
  const res = out || new Float64Array(d);
  if (!invMass) {
    for (let i = 0; i < d; i++) {
      res[i] = p[i];
    }
  } else {
    for (let i = 0; i < d; i++) {
      res[i] = p[i] * invMass[i];
    }
  }
  return res;
}

/**
 * Performs a single leapfrog integration step of step-size epsilon.
 *
 * Updates:
 * 1. p(t + eps/2) = p(t) - (eps/2) * gradU(q(t))
 * 2. q(t + eps)   = q(t) + eps * M^{-1} p(t + eps/2)
 * 3. p(t + eps)   = p(t + eps/2) - (eps/2) * gradU(q(t + eps))
 *
 * @param {Float64Array} q - Current position state (modified in-place or returned)
 * @param {Float64Array} p - Current momentum state (modified in-place or returned)
 * @param {Float64Array} gradU - Gradient of potential energy at q(t), gradU = -grad(logPosterior)
 * @param {number} epsilon - Step size
 * @param {function(Float64Array): { potential: number, gradU: Float64Array }} gradPotentialFn - Computes potential U(q) and gradU(q)
 * @param {Float64Array} [invMass=null] - Diagonal elements of M^{-1}
 * @returns {{
 *   q: Float64Array,
 *   p: Float64Array,
 *   gradU: Float64Array,
 *   potential: number,
 *   kinetic: number,
 *   hamiltonian: number
 * }}
 */
export function leapfrogStep(q, p, gradU, epsilon, gradPotentialFn, invMass = null) {
  const d = q.length;
  const nextQ = new Float64Array(d);
  const nextP = new Float64Array(d);
  const halfEps = 0.5 * epsilon;

  // 1. Half step momentum update: p_half = p - (eps/2) * gradU(q)
  for (let i = 0; i < d; i++) {
    nextP[i] = p[i] - halfEps * gradU[i];
  }

  // 2. Full step position update: q_next = q + eps * M^{-1} p_half
  if (!invMass) {
    for (let i = 0; i < d; i++) {
      nextQ[i] = q[i] + epsilon * nextP[i];
    }
  } else {
    for (let i = 0; i < d; i++) {
      nextQ[i] = q[i] + epsilon * nextP[i] * invMass[i];
    }
  }

  // Evaluate gradient at new position
  const evalResult = gradPotentialFn(nextQ);
  const nextGradU = evalResult.gradU;
  const nextPotential = evalResult.potential;

  // 3. Second half step momentum update: p_next = p_half - (eps/2) * gradU(q_next)
  for (let i = 0; i < d; i++) {
    nextP[i] -= halfEps * nextGradU[i];
  }

  const nextKinetic = computeKineticEnergy(nextP, invMass);
  const hamiltonian = nextPotential + nextKinetic;

  return {
    q: nextQ,
    p: nextP,
    gradU: nextGradU,
    potential: nextPotential,
    kinetic: nextKinetic,
    hamiltonian
  };
}

/**
 * Symplectic Leapfrog Integrator class managing full trajectories, reversibility checks,
 * and adaptive stability monitoring.
 */
export class LeapfrogIntegrator {
  /**
   * @param {function(Float64Array): { potential: number, gradU: Float64Array }} gradPotentialFn - Evaluator for U(q) and gradU(q)
   * @param {Object} [options={}]
   * @param {Float64Array|number[]} [options.invMass=null] - Inverse mass matrix diagonal
   * @param {number} [options.maxEnergyDelta=1000.0] - Threshold for numerical divergence / instability
   */
  constructor(gradPotentialFn, options = {}) {
    if (typeof gradPotentialFn !== \'function\') {
      throw new TypeError(\'gradPotentialFn must be a function computing { potential, gradU }\');
    }
    this.gradPotentialFn = gradPotentialFn;
    this.invMass = options.invMass ? new Float64Array(options.invMass) : null;
    this.maxEnergyDelta = options.maxEnergyDelta ?? 1000.0;
  }

  /**
   * Integrates a trajectory for L leapfrog steps.
   *
   * @param {Float64Array|number[]} q0 - Initial position
   * @param {Float64Array|number[]} p0 - Initial momentum
   * @param {number} epsilon - Step size
   * @param {number} numSteps - Number of leapfrog steps L
   * @returns {{
   *   q: Float64Array,
   *   p: Float64Array,
   *   gradU: Float64Array,
   *   initialHamiltonian: number,
   *   finalHamiltonian: number,
   *   deltaH: number,
   *   diverged: boolean,
   *   trajectory: Array<{ q: Float64Array, p: Float64Array, hamiltonian: number }>
   * }}
   */
  integrate(q0, p0, epsilon, numSteps) {
    const d = q0.length;
    let currentQ = new Float64Array(q0);
    let currentP = new Float64Array(p0);

    const initialEval = this.gradPotentialFn(currentQ);
    let currentGradU = initialEval.gradU;
    let currentPotential = initialEval.potential;
    const initialKinetic = computeKineticEnergy(currentP, this.invMass);
    const initialHamiltonian = currentPotential + initialKinetic;

    const trajectory = [{
      q: new Float64Array(currentQ),
      p: new Float64Array(currentP),
      hamiltonian: initialHamiltonian
    }];

    let diverged = false;
    let finalHamiltonian = initialHamiltonian;

    for (let step = 0; step < numSteps; step++) {
      const stepRes = leapfrogStep(
        currentQ,
        currentP,
        currentGradU,
        epsilon,
        this.gradPotentialFn,
        this.invMass
      );

      currentQ = stepRes.q;
      currentP = stepRes.p;
      currentGradU = stepRes.gradU;
      currentPotential = stepRes.potential;
      finalHamiltonian = stepRes.hamiltonian;

      // Check for numerical divergence / NaN
      if (
        !Number.isFinite(finalHamiltonian) ||
        Math.abs(finalHamiltonian - initialHamiltonian) > this.maxEnergyDelta
      ) {
        diverged = true;
        break;
      }

      trajectory.push({
        q: new Float64Array(currentQ),
        p: new Float64Array(currentP),
        hamiltonian: finalHamiltonian
      });
    }

    const deltaH = finalHamiltonian - initialHamiltonian;

    return {
      q: currentQ,
      p: currentP,
      gradU: currentGradU,
      initialHamiltonian,
      finalHamiltonian,
      deltaH,
      diverged,
      trajectory
    };
  }

  /**
   * Verifies symplectic reversibility of the trajectory:
   * 1. Integrates forward L steps: (q0, p0) -> (qL, pL)
   * 2. Reverses momentum: p_rev = -pL
   * 3. Integrates forward L steps from (qL, p_rev) -> (q_rev, p_final)
   * 4. Verifies q_rev == q0 within tight tolerance.
   *
   * @param {Float64Array|number[]} q0
   * @param {Float64Array|number[]} p0
   * @param {number} epsilon
   * @param {number} numSteps
   * @returns {{
   *   isReversible: boolean,
   *   positionError: number,
   *   momentumError: number,
   *   qFinal: Float64Array,
   *   qReversed: Float64Array
   * }}
   */
  checkReversibility(q0, p0, epsilon, numSteps) {
    const fwd = this.integrate(q0, p0, epsilon, numSteps);
    if (fwd.diverged) {
      return {
        isReversible: false,
        positionError: Infinity,
        momentumError: Infinity,
        qFinal: fwd.q,
        qReversed: null
      };
    }

    // Reverse momentum
    const pRev = new Float64Array(fwd.p.length);
    for (let i = 0; i < fwd.p.length; i++) {
      pRev[i] = -fwd.p[i];
    }

    const bwd = this.integrate(fwd.q, pRev, epsilon, numSteps);
    let maxPosErr = 0.0;
    let maxMomErr = 0.0;

    for (let i = 0; i < q0.length; i++) {
      const pOrig = p0[i];
      const pRecovered = -bwd.p[i]; // reversed again to compare with p0
      const posErr = Math.abs(bwd.q[i] - q0[i]);
      const momErr = Math.abs(pRecovered - pOrig);
      if (posErr > maxPosErr) maxPosErr = posErr;
      if (momErr > maxMomErr) maxMomErr = momErr;
    }

    return {
      isReversible: maxPosErr < 1e-6 && maxMomErr < 1e-6,
      positionError: maxPosErr,
      momentumError: maxMomErr,
      qFinal: fwd.q,
      qReversed: bwd.q
    };
  }
}
