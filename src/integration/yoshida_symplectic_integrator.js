/**
 * @file yoshida_symplectic_integrator.js
 * @description Symplectic Hamiltonian Integrators (Verlet, Yoshida 4th-Order, and Yoshida 6th-Order).
 * Preserves the symplectic 2-form dp ^ dq, phase-space volume (Liouville theorem), and guarantees
 * bounded energy oscillations without secular energy drift over astronomical timescales.
 *
 * Applicable to separable Hamiltonian dynamical systems:
 *   H(q, p) = T(p) + V(q) = (1/2) p^T M^{-1} p + V(q)
 * Equations of motion:
 *   dq/dt = dH/dp = M^{-1} p
 *   dp/dt = -dH/dq = -grad V(q) = F(q)
 *
 * References:
 * - Yoshida, H. (1990). "Construction of higher order symplectic integrators".
 *   Physics Letters A, 150(5-7), 262-268.
 * - Hairer, E., Lubich, C., & Wanner, G. (2006). "Geometric Numerical Integration:
 *   Structure-Preserving Algorithms for Ordinary Differential Equations". Springer.
 * - Suzuki, M. (1992). "General theory of fractal path integrals with applications
 *   to many-body theories and statistical physics". J. Math. Phys. 32, 400.
 *
 * @module integration/yoshida_symplectic_integrator
 */

import { evaluateField, isWithinBounds, TerminationReason, IntegrationDirection } from './rk4_classical.js';

export const YOSHIDA_COEFFS = Object.freeze({
  order4: {
    w1: 1.0 / (2.0 - Math.cbrt(2.0)),
    w0: -Math.cbrt(2.0) / (2.0 - Math.cbrt(2.0)),
    c: [
      0.5 / (2.0 - Math.cbrt(2.0)),
      0.5 * (1.0 - Math.cbrt(2.0)) / (2.0 - Math.cbrt(2.0)),
      0.5 * (1.0 - Math.cbrt(2.0)) / (2.0 - Math.cbrt(2.0)),
      0.5 / (2.0 - Math.cbrt(2.0))
    ],
    d: [
      1.0 / (2.0 - Math.cbrt(2.0)),
      -Math.cbrt(2.0) / (2.0 - Math.cbrt(2.0)),
      1.0 / (2.0 - Math.cbrt(2.0)),
      0.0
    ]
  },
  order6: {
    w1: -0.117767998417887e1,
    w2: 0.235573213359357e0,
    w3: 0.784513610477560e0,
    w0: 1.315682772563908e0,
    weights: [
      0.784513610477560e0,
      0.235573213359357e0,
      -0.117767998417887e1,
      1.315682772563908e0,
      -0.117767998417887e1,
      0.235573213359357e0,
      0.784513610477560e0
    ]
  },
  order8: {
    weights: [
      0.13020248308889008088e1,
      0.56116958178800328403e0,
      -0.38947796313452054140e0,
      0.15884190725516494080e0,
      -0.39590389413865757700e0,
      0.18453964097831570700e1,
      -0.25814578228248483320e0,
      -0.25048809405445690020e1
    ]
  }
});

export function verletStep(forceFn, q, p, t, dt, mass = 1.0) {
  const invM = 1.0 / mass;
  const dtHalf = 0.5 * dt;

  const qHalf = [
    q[0] + dtHalf * p[0] * invM,
    q[1] + dtHalf * p[1] * invM,
    q[2] + dtHalf * p[2] * invM
  ];

  const f = evaluateField(forceFn, t + dtHalf, qHalf);
  const pNext = [
    p[0] + dt * f[0],
    p[1] + dt * f[1],
    p[2] + dt * f[2]
  ];

  const qNext = [
    qHalf[0] + dtHalf * pNext[0] * invM,
    qHalf[1] + dtHalf * pNext[1] * invM,
    qHalf[2] + dtHalf * pNext[2] * invM
  ];

  return {
    qNext,
    pNext,
    tNext: t + dt,
    force: f
  };
}

export function yoshida4Step(forceFn, q, p, t, dt, mass = 1.0) {
  const { w1, w0 } = YOSHIDA_COEFFS.order4;

  let curQ = [q[0], q[1], q[2]];
  let curP = [p[0], p[1], p[2]];
  let curT = t;

  const weights = [w1, w0, w1];

  for (let i = 0; i < 3; i++) {
    const stepDt = weights[i] * dt;
    const res = verletStep(forceFn, curQ, curP, curT, stepDt, mass);
    curQ = res.qNext;
    curP = res.pNext;
    curT += stepDt;
  }

  return {
    qNext: curQ,
    pNext: curP,
    tNext: t + dt
  };
}

export function yoshida6Step(forceFn, q, p, t, dt, mass = 1.0) {
  const weights = YOSHIDA_COEFFS.order6.weights;
  let curQ = [q[0], q[1], q[2]];
  let curP = [p[0], p[1], p[2]];
  let curT = t;

  for (let i = 0; i < weights.length; i++) {
    const stepDt = weights[i] * dt;
    const res = yoshida4Step(forceFn, curQ, curP, curT, stepDt, mass);
    curQ = res.qNext;
    curP = res.pNext;
    curT += stepDt;
  }

  return {
    qNext: curQ,
    pNext: curP,
    tNext: t + dt
  };
}

export function computeHamiltonianEnergy(q, p, potentialFn, mass = 1.0) {
  const pSq = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
  const kinetic = 0.5 * pSq / mass;
  const potential = typeof potentialFn === 'function' ? potentialFn(q) : 0.0;
  return kinetic + potential;
}

export function computeAngularMomentum(q, p) {
  return [
    q[1] * p[2] - q[2] * p[1],
    q[2] * p[0] - q[0] * p[2],
    q[0] * p[1] - q[1] * p[0]
  ];
}

export function computeRungeLenzVector(q, p, GM, mass = 1.0) {
  const L = computeAngularMomentum(q, p);
  const pxL = [
    p[1] * L[2] - p[2] * L[1],
    p[2] * L[0] - p[0] * L[2],
    p[0] * L[1] - p[1] * L[0]
  ];

  const r = Math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
  if (r < 1e-12) return pxL;

  const factor = (GM * mass * mass) / r;
  return [
    pxL[0] - factor * q[0],
    pxL[1] - factor * q[1],
    pxL[2] - factor * q[2]
  ];
}

export function computeShadowHamiltonian(q, p, potentialFn, forceFn, dt, mass = 1.0) {
  const H0 = computeHamiltonianEnergy(q, p, potentialFn, mass);
  const f = evaluateField(forceFn, 0, q);
  const fSq = f[0] * f[0] + f[1] * f[1] + f[2] * f[2];
  const H_corr = (dt * dt / (24.0 * mass)) * fSq;
  return H0 + H_corr;
}

export function verifySymplecticJacobian(stepFn, forceFn, q0, p0, dt, mass = 1.0, delta = 1e-6) {
  const z0 = [...q0, ...p0];
  const J = [];

  for (let col = 0; col < 6; col++) {
    const zPlus = [...z0];
    const zMinus = [...z0];
    zPlus[col] += delta;
    zMinus[col] -= delta;

    const qPlus = zPlus.slice(0, 3);
    const pPlus = zPlus.slice(3, 6);
    const qMinus = zMinus.slice(0, 3);
    const pMinus = zMinus.slice(3, 6);

    const outPlus = stepFn(forceFn, qPlus, pPlus, 0.0, dt, mass);
    const outMinus = stepFn(forceFn, qMinus, pMinus, 0.0, dt, mass);

    const zOutPlus = [...outPlus.qNext, ...outPlus.pNext];
    const zOutMinus = [...outMinus.qNext, ...outMinus.pNext];

    const column = [];
    for (let row = 0; row < 6; row++) {
      column.push((zOutPlus[row] - zOutMinus[row]) / (2.0 * delta));
    }
    J.push(column);
  }

  const jacobian = [];
  for (let r = 0; r < 6; r++) {
    jacobian[r] = [];
    for (let c = 0; c < 6; c++) {
      jacobian[r][c] = J[c][r];
    }
  }

  let det = 1.0;
  const A = jacobian.map(row => [...row]);
  const n = 6;

  for (let i = 0; i < n; i++) {
    let pivot = i;
    for (let j = i + 1; j < n; j++) {
      if (Math.abs(A[j][i]) > Math.abs(A[pivot][i])) {
        pivot = j;
      }
    }

    if (Math.abs(A[pivot][i]) < 1e-12) {
      det = 0.0;
      break;
    }

    if (pivot !== i) {
      const tmp = A[i];
      A[i] = A[pivot];
      A[pivot] = tmp;
      det = -det;
    }

    det *= A[i][i];
    const pivotVal = A[i][i];

    for (let j = i + 1; j < n; j++) {
      const factor = A[j][i] / pivotVal;
      for (let k = i; k < n; k++) {
        A[j][k] -= factor * A[i][k];
      }
    }
  }

  return {
    jacobian,
    determinant: det,
    symplecticResidual: Math.abs(det - 1.0)
  };
}

export function traceSymplecticTrajectory(forceFn, seedQ, seedP, options = {}) {
  const {
    order = 4,
    dt = 0.1,
    maxSteps = 2000,
    mass = 1.0,
    potentialFn = null,
    GM = 100.0,
    bounds = null,
    direction = IntegrationDirection.FORWARD
  } = options;

  const dtSigned = direction === IntegrationDirection.BACKWARD ? -Math.abs(dt) : Math.abs(dt);
  const stepFn = order === 6 ? yoshida6Step : (order === 2 ? verletStep : yoshida4Step);

  let curQ = [Number(seedQ[0]), Number(seedQ[1]), Number(seedQ[2])];
  let curP = [Number(seedP[0]), Number(seedP[1]), Number(seedP[2])];
  let curT = 0.0;

  const positions = [[...curQ]];
  const momenta = [[...curP]];
  const times = [curT];
  const energies = [];
  const angularMomenta = [];
  const arcLengths = [0.0];

  let cumulativeArc = 0.0;
  let terminationReason = TerminationReason.MAX_STEPS;

  if (potentialFn) {
    energies.push(computeHamiltonianEnergy(curQ, curP, potentialFn, mass));
  }
  angularMomenta.push(computeAngularMomentum(curQ, curP));

  for (let step = 0; step < maxSteps; step++) {
    if (!isWithinBounds(curQ, bounds)) {
      terminationReason = TerminationReason.DOMAIN_EXIT;
      break;
    }

    try {
      const res = stepFn(forceFn, curQ, curP, curT, dtSigned, mass);
      const nextQ = res.qNext;
      const nextP = res.pNext;

      const dx = nextQ[0] - curQ[0];
      const dy = nextQ[1] - curQ[1];
      const dz = nextQ[2] - curQ[2];
      const stepDist = Math.sqrt(dx * dx + dy * dy + dz * dz);
      cumulativeArc += stepDist;

      curQ = nextQ;
      curP = nextP;
      curT = res.tNext;

      positions.push([...curQ]);
      momenta.push([...curP]);
      times.push(curT);
      arcLengths.push(cumulativeArc);

      if (potentialFn) {
        energies.push(computeHamiltonianEnergy(curQ, curP, potentialFn, mass));
      }
      angularMomenta.push(computeAngularMomentum(curQ, curP));
    } catch (err) {
      terminationReason = TerminationReason.NUMERICAL_FAILURE;
      break;
    }
  }

  let energyDriftMax = 0.0;
  let energyDriftRelMax = 0.0;
  if (energies.length > 1) {
    const e0 = energies[0];
    for (let i = 1; i < energies.length; i++) {
      const dE = Math.abs(energies[i] - e0);
      if (dE > energyDriftMax) energyDriftMax = dE;
      const rel = Math.abs(e0) > 1e-12 ? dE / Math.abs(e0) : dE;
      if (rel > energyDriftRelMax) energyDriftRelMax = rel;
    }
  }

  return {
    order,
    positions,
    momenta,
    times,
    arcLengths,
    energies,
    angularMomenta,
    totalSteps: positions.length - 1,
    totalArcLength: cumulativeArc,
    energyDriftMax,
    energyDriftRelMax,
    terminationReason
  };
}

export class YoshidaSymplecticIntegrator {
  constructor(config = {}) {
    this.order = config.order ?? 4;
    this.dt = config.dt ?? 0.1;
    this.mass = config.mass ?? 1.0;
    this.maxSteps = config.maxSteps ?? 2000;
    this.potentialFn = config.potentialFn ?? null;
    this.bounds = config.bounds ?? null;
    this.direction = config.direction ?? IntegrationDirection.FORWARD;
  }

  step(forceFn, q, p, t = 0.0, dt = this.dt, mass = this.mass) {
    if (this.order === 6) {
      return yoshida6Step(forceFn, q, p, t, dt, mass);
    } else if (this.order === 2) {
      return verletStep(forceFn, q, p, t, dt, mass);
    }
    return yoshida4Step(forceFn, q, p, t, dt, mass);
  }

  trace(forceFn, seedQ, seedP, overrideOptions = {}) {
    const opts = {
      order: this.order,
      dt: this.dt,
      mass: this.mass,
      maxSteps: this.maxSteps,
      potentialFn: this.potentialFn,
      bounds: this.bounds,
      direction: this.direction,
      ...overrideOptions
    };
    return traceSymplecticTrajectory(forceFn, seedQ, seedP, opts);
  }

  integrate(forceFn, seedQ, seedP, overrideOptions = {}) {
    return this.trace(forceFn, seedQ, seedP, overrideOptions);
  }
}

export default YoshidaSymplecticIntegrator;
