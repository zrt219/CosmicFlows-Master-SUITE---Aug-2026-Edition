"""
Analytic Vector Fields and Differential Invariant Verification Test Suite
=======================================================================
Mathematical test fixtures and numerical validation for canonical 3D vector fields:
  1. Constant laminar flow (zero divergence, zero vorticity, no critical points)
  2. Linear 3D sink field v(x) = -A x (attractor at origin, all negative eigenvalues)
  3. Linear 3D source field v(x) = +A x (repeller at origin, all positive eigenvalues)
  4. Saddle-Filament field (two negative eigenvalues, one positive)
  5. Saddle-Wall field (one negative eigenvalue, two positive)
  6. Pure rotational vortex field (non-zero curl, testing vorticity guards)
  7. 3D Plummer gravitational cluster model (smooth inflow with softened core)

Provides complete, exact mathematical implementations, high-order finite-difference stencils,
Cash-Karp RK45 streamline integrators, Jacobi / QR eigenvalue decomposition, and Newton-Raphson
critical point locators.
"""

import math
import numpy as np
import pytest
from typing import Callable, Tuple, List, Dict, Any, Optional


# =============================================================================
# 1. MATHEMATICAL VECTOR FIELD FIXTURES (ANALYTICAL FORMULATIONS)
# =============================================================================

class LaminarFlowField:
    """
    Constant Laminar Flow Field:
      v(x, y, z) = (u_0, v_0, w_0) = constant vector.
    Invariants:
      div(v) = 0
      curl(v) = (0, 0, 0)
      J(x) = 0 (3x3 null matrix)
      No critical points in R^3 when ||v_0|| > 0.
      Exact streamline: x(t) = x_0 + v_0 * t.
    """
    def __init__(self, u0: float = 250.0, v0: float = -150.0, w0: float = 75.0):
        self.v0 = np.array([u0, v0, w0], dtype=np.float64)
        if np.linalg.norm(self.v0) < 1e-12:
            raise ValueError("Laminar velocity vector must be non-zero for laminar flow fixture.")

    def velocity(self, x: np.ndarray) -> np.ndarray:
        return np.copy(self.v0)

    def analytical_divergence(self, x: np.ndarray) -> float:
        return 0.0

    def analytical_vorticity(self, x: np.ndarray) -> np.ndarray:
        return np.zeros(3, dtype=np.float64)

    def analytical_jacobian(self, x: np.ndarray) -> np.ndarray:
        return np.zeros((3, 3), dtype=np.float64)

    def analytical_streamline(self, x0: np.ndarray, t: float) -> np.ndarray:
        return x0 + self.v0 * t


class LinearSinkField:
    """
    Linear 3D Sink / Attractor Field:
      v(x) = -A x, where A = diag(lambda_1, lambda_2, lambda_3), lambda_i > 0.
    Invariants:
      Critical point at x* = (0, 0, 0).
      Jacobian J(x) = -A = diag(-lambda_1, -lambda_2, -lambda_3).
      All eigenvalues mu_i = -lambda_i < 0 (asymptotically stable sink / node).
      div(v) = -(lambda_1 + lambda_2 + lambda_3) < 0 everywhere.
      curl(v) = (0, 0, 0) everywhere.
      Exact streamline: x_i(t) = x_{i,0} * exp(-lambda_i * t).
    """
    def __init__(self, l1: float = 1.5, l2: float = 2.0, l3: float = 0.8):
        if l1 <= 0 or l2 <= 0 or l3 <= 0:
            raise ValueError("Sink field eigenvalues lambda_i must be strictly positive.")
        self.lambdas = np.array([l1, l2, l3], dtype=np.float64)
        self.A = np.diag(self.lambdas)

    def velocity(self, x: np.ndarray) -> np.ndarray:
        return -self.A @ np.asarray(x, dtype=np.float64)

    def analytical_divergence(self, x: np.ndarray) -> float:
        return float(-np.sum(self.lambdas))

    def analytical_vorticity(self, x: np.ndarray) -> np.ndarray:
        return np.zeros(3, dtype=np.float64)

    def analytical_jacobian(self, x: np.ndarray) -> np.ndarray:
        return -np.copy(self.A)

    def analytical_eigenvalues(self) -> np.ndarray:
        return -np.copy(self.lambdas)

    def analytical_streamline(self, x0: np.ndarray, t: float) -> np.ndarray:
        return np.asarray(x0, dtype=np.float64) * np.exp(-self.lambdas * t)


class LinearSourceField:
    """
    Linear 3D Source / Repeller Field:
      v(x) = +A x, where A = diag(lambda_1, lambda_2, lambda_3), lambda_i > 0.
    Invariants:
      Critical point at x* = (0, 0, 0).
      Jacobian J(x) = +A = diag(+lambda_1, +lambda_2, +lambda_3).
      All eigenvalues mu_i = +lambda_i > 0 (unstable source / repeller).
      div(v) = +(lambda_1 + lambda_2 + lambda_3) > 0 everywhere.
      curl(v) = (0, 0, 0) everywhere.
      Exact streamline: x_i(t) = x_{i,0} * exp(+lambda_i * t).
    """
    def __init__(self, l1: float = 1.1, l2: float = 1.4, l3: float = 1.9):
        if l1 <= 0 or l2 <= 0 or l3 <= 0:
            raise ValueError("Source field eigenvalues lambda_i must be strictly positive.")
        self.lambdas = np.array([l1, l2, l3], dtype=np.float64)
        self.A = np.diag(self.lambdas)

    def velocity(self, x: np.ndarray) -> np.ndarray:
        return self.A @ np.asarray(x, dtype=np.float64)

    def analytical_divergence(self, x: np.ndarray) -> float:
        return float(np.sum(self.lambdas))

    def analytical_vorticity(self, x: np.ndarray) -> np.ndarray:
        return np.zeros(3, dtype=np.float64)

    def analytical_jacobian(self, x: np.ndarray) -> np.ndarray:
        return np.copy(self.A)

    def analytical_eigenvalues(self) -> np.ndarray:
        return np.copy(self.lambdas)

    def analytical_streamline(self, x0: np.ndarray, t: float) -> np.ndarray:
        return np.asarray(x0, dtype=np.float64) * np.exp(self.lambdas * t)


class SaddleFilamentField:
    """
    Saddle-Filament Field (Cosmic Filament Spine):
      v(x) = (-lambda_1 x, -lambda_2 y, +lambda_3 z)
      with lambda_1, lambda_2, lambda_3 > 0.
    Invariants:
      Two negative eigenvalues (xy-plane collapse) and one positive eigenvalue (z-axis outflow along filament spine).
      Critical point at (0, 0, 0).
      Jacobian J = diag(-lambda_1, -lambda_2, +lambda_3).
      Eigenvalues: mu_1 = -lambda_1 < 0, mu_2 = -lambda_2 < 0, mu_3 = +lambda_3 > 0.
      div(v) = -lambda_1 - lambda_2 + lambda_3.
      curl(v) = (0, 0, 0).
      Exact streamline: x(t) = x0*e^(-l1 t), y(t) = y0*e^(-l2 t), z(t) = z0*e^(+l3 t).
    """
    def __init__(self, l1: float = 1.0, l2: float = 1.2, l3: float = 0.8):
        if l1 <= 0 or l2 <= 0 or l3 <= 0:
            raise ValueError("Eigenvalue parameters must be positive.")
        self.l1 = float(l1)
        self.l2 = float(l2)
        self.l3 = float(l3)
        self.A = np.diag([-self.l1, -self.l2, +self.l3])

    def velocity(self, x: np.ndarray) -> np.ndarray:
        return self.A @ np.asarray(x, dtype=np.float64)

    def analytical_divergence(self, x: np.ndarray) -> float:
        return float(-self.l1 - self.l2 + self.l3)

    def analytical_vorticity(self, x: np.ndarray) -> np.ndarray:
        return np.zeros(3, dtype=np.float64)

    def analytical_jacobian(self, x: np.ndarray) -> np.ndarray:
        return np.copy(self.A)

    def analytical_eigenvalues(self) -> np.ndarray:
        return np.array([-self.l1, -self.l2, +self.l3], dtype=np.float64)

    def analytical_streamline(self, x0: np.ndarray, t: float) -> np.ndarray:
        x0 = np.asarray(x0, dtype=np.float64)
        return np.array([
            x0[0] * math.exp(-self.l1 * t),
            x0[1] * math.exp(-self.l2 * t),
            x0[2] * math.exp(+self.l3 * t)
        ], dtype=np.float64)


class SaddleWallField:
    """
    Saddle-Wall Field (Cosmic Sheet / Zeldovich Pancake):
      v(x) = (-lambda_1 x, +lambda_2 y, +lambda_3 z)
      with lambda_1, lambda_2, lambda_3 > 0.
    Invariants:
      One negative eigenvalue (compression onto wall at x=0) and two positive eigenvalues (expansion within yz-plane).
      Critical point at (0, 0, 0).
      Jacobian J = diag(-lambda_1, +lambda_2, +lambda_3).
      Eigenvalues: mu_1 = -lambda_1 < 0, mu_2 = +lambda_2 > 0, mu_3 = +lambda_3 > 0.
      div(v) = -lambda_1 + lambda_2 + lambda_3.
      curl(v) = (0, 0, 0).
    """
    def __init__(self, l1: float = 1.5, l2: float = 0.7, l3: float = 0.9):
        if l1 <= 0 or l2 <= 0 or l3 <= 0:
            raise ValueError("Eigenvalue parameters must be positive.")
        self.l1 = float(l1)
        self.l2 = float(l2)
        self.l3 = float(l3)
        self.A = np.diag([-self.l1, +self.l2, +self.l3])

    def velocity(self, x: np.ndarray) -> np.ndarray:
        return self.A @ np.asarray(x, dtype=np.float64)

    def analytical_divergence(self, x: np.ndarray) -> float:
        return float(-self.l1 + self.l2 + self.l3)

    def analytical_vorticity(self, x: np.ndarray) -> np.ndarray:
        return np.zeros(3, dtype=np.float64)

    def analytical_jacobian(self, x: np.ndarray) -> np.ndarray:
        return np.copy(self.A)

    def analytical_eigenvalues(self) -> np.ndarray:
        return np.array([-self.l1, +self.l2, +self.l3], dtype=np.float64)

    def analytical_streamline(self, x0: np.ndarray, t: float) -> np.ndarray:
        x0 = np.asarray(x0, dtype=np.float64)
        return np.array([
            x0[0] * math.exp(-self.l1 * t),
            x0[1] * math.exp(+self.l2 * t),
            x0[2] * math.exp(+self.l3 * t)
        ], dtype=np.float64)


class RotationalVortexField:
    """
    Pure Rotational Vortex Field (Rigid / Beltrami rotation):
      v(x, y, z) = omega x r = (-omega_z * y + omega_y * z,
                                +omega_z * x - omega_x * z,
                                -omega_y * x + omega_x * y)
      For omega = (0, 0, omega_0): v(x, y, z) = (-omega_0 * y, +omega_0 * x, 0).
    Invariants:
      Divergence div(v) = 0 identically (strictly incompressible).
      Vorticity curl(v) = 2 * omega = (2*omega_x, 2*omega_y, 2*omega_z) everywhere.
      Jacobian J = [[0, -omega_z, omega_y], [omega_z, 0, -omega_x], [-omega_y, omega_x, 0]] (antisymmetric / skew-symmetric).
      Eigenvalues: 0 and +/- i * ||omega|| (pure imaginary pairs, zero trace).
      Closed circular streamlines with period T = 2*pi / ||omega||.
    """
    def __init__(self, omega_vector: Tuple[float, float, float] = (0.0, 0.0, 2.5)):
        self.omega = np.array(omega_vector, dtype=np.float64)
        self.omega_norm = float(np.linalg.norm(self.omega))
        if self.omega_norm < 1e-12:
            raise ValueError("Vorticity vector must be non-zero for rotational vortex fixture.")

    def velocity(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=np.float64)
        return np.cross(self.omega, x)

    def analytical_divergence(self, x: np.ndarray) -> float:
        return 0.0

    def analytical_vorticity(self, x: np.ndarray) -> np.ndarray:
        return 2.0 * self.omega

    def analytical_jacobian(self, x: np.ndarray) -> np.ndarray:
        wx, wy, wz = self.omega
        return np.array([
            [0.0, -wz, wy],
            [wz, 0.0, -wx],
            [-wy, wx, 0.0]
        ], dtype=np.float64)

    def analytical_streamline(self, x0: np.ndarray, t: float) -> np.ndarray:
        """
        Rodrigues' rotation formula for exact rotational streamline:
          x(t) = x0*cos(omega*t) + (k x x0)*sin(omega*t) + k*(k . x0)*(1 - cos(omega*t))
        where k = omega / ||omega||.
        """
        x0 = np.asarray(x0, dtype=np.float64)
        k = self.omega / self.omega_norm
        theta = self.omega_norm * t
        return (x0 * math.cos(theta) +
                np.cross(k, x0) * math.sin(theta) +
                k * np.dot(k, x0) * (1.0 - math.cos(theta)))


class PlummerGravitationalCluster:
    """
    3D Plummer Gravitational Cluster Model (Smooth Inflow with Softened Core):
      Potential: Phi(r) = - G * M / sqrt(r^2 + eps^2)
      Velocity Field (inflow): v(x) = -nabla Phi(x) = - G * M * (x - x_c) / (|x - x_c|^2 + eps^2)^(3/2)
    Invariants:
      Smooth C^inf everywhere with core softening length eps > 0.
      Center convergence: as x -> x_c, v(x) -> 0 smoothly (no singularity, no jump discontinuity).
      Peak inflow velocity occurs at r_peak = eps / sqrt(2) approx 0.7071 * eps.
      Analytical Divergence: div(v) = - 3 * G * M * eps^2 / (r^2 + eps^2)^(5/2) < 0 everywhere.
      Analytical Vorticity: curl(v) = (0, 0, 0) identically (irrotational conservative potential flow).
      Jacobian: J_ij(x) = - (G*M / (r^2 + eps^2)^(3/2)) * delta_ij + (3 * G * M * (x_i - x_{c,i})(x_j - x_{c,j}) / (r^2 + eps^2)^(5/2)).
      At center x = x_c: J = - (G*M / eps^3) * I_{3x3}, triple degenerate negative eigenvalues mu = -G*M / eps^3.
    """
    def __init__(self, center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                 GM: float = 1.0e6,
                 epsilon: float = 500.0):
        if epsilon <= 0:
            raise ValueError("Plummer core softening epsilon must be strictly positive.")
        if GM <= 0:
            raise ValueError("Cluster gravitational strength GM must be strictly positive.")
        self.center = np.array(center, dtype=np.float64)
        self.GM = float(GM)
        self.eps = float(epsilon)
        self.eps_sq = self.eps * self.eps

    def displacement(self, x: np.ndarray) -> np.ndarray:
        return np.asarray(x, dtype=np.float64) - self.center

    def velocity(self, x: np.ndarray) -> np.ndarray:
        r_vec = self.displacement(x)
        r_sq = float(np.dot(r_vec, r_vec))
        denom = (r_sq + self.eps_sq) ** 1.5
        return -(self.GM / denom) * r_vec

    def analytical_divergence(self, x: np.ndarray) -> float:
        r_vec = self.displacement(x)
        r_sq = float(np.dot(r_vec, r_vec))
        denom = (r_sq + self.eps_sq) ** 2.5
        return float(-3.0 * self.GM * self.eps_sq / denom)

    def analytical_vorticity(self, x: np.ndarray) -> np.ndarray:
        return np.zeros(3, dtype=np.float64)

    def analytical_jacobian(self, x: np.ndarray) -> np.ndarray:
        r_vec = self.displacement(x)
        r_sq = float(np.dot(r_vec, r_vec))
        r2_eps2 = r_sq + self.eps_sq
        denom_iso = r2_eps2 ** 1.5
        denom_ani = r2_eps2 ** 2.5
        iso_term = -(self.GM / denom_iso) * np.eye(3, dtype=np.float64)
        ani_term = (3.0 * self.GM / denom_ani) * np.outer(r_vec, r_vec)
        return iso_term + ani_term

    def peak_radius(self) -> float:
        return self.eps / math.sqrt(2.0)

    def peak_velocity_magnitude(self) -> float:
        r = self.peak_radius()
        r_sq = r * r
        return (self.GM * r) / ((r_sq + self.eps_sq) ** 1.5)


# =============================================================================
# 2. HIGH-PRECISION NUMERICAL DIFFERENTIAL OPERATORS & ALGORITHMS
# =============================================================================

def compute_numerical_jacobian(field_fn: Callable[[np.ndarray], np.ndarray],
                               x: np.ndarray,
                               h: float = 1e-5,
                               stencil_order: int = 4) -> np.ndarray:
    """
    Computes the 3x3 Jacobian matrix J_ij = dv_i / dx_j using 2nd or 4th-order central finite differences.
    """
    x = np.asarray(x, dtype=np.float64)
    J = np.zeros((3, 3), dtype=np.float64)
    
    if stencil_order == 4:
        # 4th order central difference coefficients: (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
        for j in range(3):
            e_j = np.zeros(3, dtype=np.float64)
            e_j[j] = 1.0
            f_p2 = field_fn(x + 2.0 * h * e_j)
            f_p1 = field_fn(x + 1.0 * h * e_j)
            f_m1 = field_fn(x - 1.0 * h * e_j)
            f_m2 = field_fn(x - 2.0 * h * e_j)
            J[:, j] = (-f_p2 + 8.0 * f_p1 - 8.0 * f_m1 + f_m2) / (12.0 * h)
    else:
        # Standard 2nd order central difference: (f(x+h) - f(x-h)) / (2h)
        for j in range(3):
            e_j = np.zeros(3, dtype=np.float64)
            e_j[j] = 1.0
            f_p = field_fn(x + h * e_j)
            f_m = field_fn(x - h * e_j)
            J[:, j] = (f_p - f_m) / (2.0 * h)
            
    return J


def compute_numerical_divergence(field_fn: Callable[[np.ndarray], np.ndarray],
                                 x: np.ndarray,
                                 h: float = 1e-5) -> float:
    """
    Computes divergence div(v) = Tr(J) = dv_x/dx + dv_y/dy + dv_z/dz via 4th-order finite difference.
    """
    J = compute_numerical_jacobian(field_fn, x, h=h, stencil_order=4)
    return float(np.trace(J))


def compute_numerical_vorticity(field_fn: Callable[[np.ndarray], np.ndarray],
                                x: np.ndarray,
                                h: float = 1e-5) -> np.ndarray:
    """
    Computes curl / vorticity curl(v) = (dv_z/dy - dv_y/dz, dv_x/dz - dv_z/dx, dv_y/dx - dv_x/dy).
    """
    J = compute_numerical_jacobian(field_fn, x, h=h, stencil_order=4)
    curl_x = J[2, 1] - J[1, 2]
    curl_y = J[0, 2] - J[2, 0]
    curl_z = J[1, 0] - J[0, 1]
    return np.array([curl_x, curl_y, curl_z], dtype=np.float64)


def jacobi_eigenvalue_diagonalization(A: np.ndarray,
                                      max_iter: int = 100,
                                      tol: float = 1e-15) -> Tuple[np.ndarray, np.ndarray]:
    """
    Exact Classical Jacobi Eigenvalue Algorithm for real symmetric 3x3 matrices.
    Returns sorted eigenvalues and corresponding orthogonal eigenvector columns.
    """
    n = A.shape[0]
    if n != 3 or A.shape[1] != 3:
        raise ValueError("Jacobi diagonalization expects 3x3 matrix.")
    
    # Symmetrize to guard against floating-point asymmetries
    A_mat = 0.5 * (A + A.T)
    V = np.eye(3, dtype=np.float64)
    
    for _ in range(max_iter):
        # Find maximum off-diagonal element
        p, q = 0, 1
        max_val = abs(A_mat[0, 1])
        if abs(A_mat[0, 2]) > max_val:
            p, q = 0, 2
            max_val = abs(A_mat[0, 2])
        if abs(A_mat[1, 2]) > max_val:
            p, q = 1, 2
            max_val = abs(A_mat[1, 2])
            
        if max_val < tol:
            break
            
        # Compute Jacobi rotation angle
        app = A_mat[p, p]
        aqq = A_mat[q, q]
        apq = A_mat[p, q]
        
        theta = 0.5 * math.atan2(2.0 * apq, (aqq - app))
        c = math.cos(theta)
        s = math.sin(theta)
        
        # Apply Givens / Jacobi similarity rotation
        J_rot = np.eye(3, dtype=np.float64)
        J_rot[p, p] = c
        J_rot[q, q] = c
        J_rot[p, q] = s
        J_rot[q, p] = -s
        
        A_mat = J_rot.T @ A_mat @ J_rot
        V = V @ J_rot
        
    eigenvalues = np.diag(A_mat)
    # Sort descending
    sort_idx = np.argsort(eigenvalues)[::-1]
    return eigenvalues[sort_idx], V[:, sort_idx]


def newton_raphson_critical_point_3d(field_fn: Callable[[np.ndarray], np.ndarray],
                                     x_init: np.ndarray,
                                     tol: float = 1e-12,
                                     max_iter: int = 50) -> Tuple[np.ndarray, bool, int]:
    """
    3D Newton-Raphson vector root locator for critical points v(x*) = 0.
      x_{k+1} = x_k - J(x_k)^(-1) * v(x_k)
    Includes backtracking line search for global convergence.
    """
    x = np.asarray(x_init, dtype=np.float64).copy()
    
    for iteration in range(max_iter):
        v = field_fn(x)
        norm_v = np.linalg.norm(v)
        if norm_v < tol:
            return x, True, iteration
            
        J = compute_numerical_jacobian(field_fn, x, h=1e-6, stencil_order=4)
        cond = np.linalg.cond(J)
        if cond > 1e14 or np.isnan(cond):
            # Regularize singular or near-singular Jacobian
            J_reg = J + 1e-8 * np.eye(3)
            delta = np.linalg.solve(J_reg, v)
        else:
            delta = np.linalg.solve(J, v)
            
        # Backtracking line search
        step_alpha = 1.0
        x_cand = x - step_alpha * delta
        while np.linalg.norm(field_fn(x_cand)) >= norm_v and step_alpha > 1e-4:
            step_alpha *= 0.5
            x_cand = x - step_alpha * delta
            
        x = x_cand
        
    v_final = field_fn(x)
    converged = np.linalg.norm(v_final) < tol
    return x, converged, max_iter


def rk45_cash_karp_integrate(field_fn: Callable[[np.ndarray], np.ndarray],
                             x0: np.ndarray,
                             t_end: float,
                             h_init: float = 0.05,
                             tol: float = 1e-8,
                             max_steps: int = 5000) -> Tuple[np.ndarray, np.ndarray]:
    """
    Cash-Karp Embedded Runge-Kutta 4(5) Streamline Integrator.
    Butcher tableau with adaptive stepsize control and local truncation error estimation.
    """
    # Butcher Tableau coefficients for Cash-Karp RK45
    c = np.array([0.0, 1.0/5.0, 3.0/10.0, 3.0/5.0, 1.0, 7.0/8.0])
    a = [
        [],
        [1.0/5.0],
        [3.0/40.0, 9.0/40.0],
        [3.0/10.0, -9.0/10.0, 6.0/5.0],
        [-11.0/54.0, 5.0/2.0, -70.0/27.0, 35.0/27.0],
        [1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0]
    ]
    b5 = np.array([37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0])
    b4 = np.array([2825.0/27648.0, 0.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 1.0/4.0])
    e_err = b5 - b4
    
    t = 0.0
    x = np.asarray(x0, dtype=np.float64).copy()
    h = h_init
    
    t_history = [t]
    x_history = [x.copy()]
    
    steps = 0
    while t < t_end and steps < max_steps:
        steps += 1
        if t + h > t_end:
            h = t_end - t
            
        k = np.zeros((6, 3), dtype=np.float64)
        k[0] = field_fn(x)
        for i in range(1, 6):
            x_stage = x.copy()
            for j in range(i):
                x_stage += h * a[i][j] * k[j]
            k[i] = field_fn(x_stage)
            
        # 5th-order candidate step and truncation error
        x5 = x + h * np.sum(b5[:, None] * k, axis=0)
        err = h * np.sum(e_err[:, None] * k, axis=0)
        err_norm = np.linalg.norm(err) / (tol * (1.0 + np.linalg.norm(x)))
        
        if err_norm <= 1.0 or h < 1e-12:
            # Step accepted
            t += h
            x = x5
            t_history.append(t)
            x_history.append(x.copy())
            
        # Adaptive step size adjustment
        if err_norm > 0:
            scale = 0.9 * (err_norm ** -0.2)
            scale = max(0.1, min(5.0, scale))
            h = h * scale
        else:
            h = h * 2.0
            
    return np.array(t_history), np.array(x_history)


# =============================================================================
# 3. UNIT TESTS FOR CONSTANT LAMINAR FLOW
# =============================================================================

class TestConstantLaminarFlow:
    """Test suite for constant laminar flow fields."""

    @pytest.fixture
    def laminar(self):
        return LaminarFlowField(u0=300.0, v0=-200.0, w0=100.0)

    def test_laminar_zero_divergence_across_domain(self, laminar):
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1200.0, -3400.0, 5600.0]),
            np.array([-5000.0, 8000.0, -12000.0]),
            np.array([9999.0, 9999.0, 9999.0])
        ]
        for pt in points:
            ana_div = laminar.analytical_divergence(pt)
            num_div = compute_numerical_divergence(laminar.velocity, pt, h=1e-4)
            assert ana_div == 0.0
            assert abs(num_div) < 1e-11, f"Expected 0 divergence, got {num_div} at {pt}"

    def test_laminar_zero_vorticity_across_domain(self, laminar):
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([500.0, -500.0, 500.0]),
            np.array([-1500.0, 2500.0, -3500.0])
        ]
        for pt in points:
            ana_curl = laminar.analytical_vorticity(pt)
            num_curl = compute_numerical_vorticity(laminar.velocity, pt, h=1e-4)
            np.testing.assert_array_equal(ana_curl, np.zeros(3))
            np.testing.assert_allclose(num_curl, np.zeros(3), atol=1e-11)

    def test_laminar_null_jacobian(self, laminar):
        pt = np.array([123.4, 567.8, -910.1])
        J_ana = laminar.analytical_jacobian(pt)
        J_num = compute_numerical_jacobian(laminar.velocity, pt, h=1e-4)
        np.testing.assert_array_equal(J_ana, np.zeros((3, 3)))
        np.testing.assert_allclose(J_num, np.zeros((3, 3)), atol=1e-11)

    def test_laminar_no_critical_points(self, laminar):
        for x_init in [np.zeros(3), np.array([100.0, -100.0, 200.0])]:
            _, converged, _ = newton_raphson_critical_point_3d(laminar.velocity, x_init, max_iter=20)
            assert not converged, "Constant laminar flow must not possess any critical points."

    def test_laminar_rk45_streamline_matches_exact_trajectory(self, laminar):
        x0 = np.array([10.0, -20.0, 30.0])
        t_end = 5.0
        t_hist, x_hist = rk45_cash_karp_integrate(laminar.velocity, x0, t_end=t_end, tol=1e-10)
        x_exact = laminar.analytical_streamline(x0, t_end)
        np.testing.assert_allclose(x_hist[-1], x_exact, rtol=1e-8, atol=1e-8)


# =============================================================================
# 4. UNIT TESTS FOR LINEAR 3D SINK FIELD (ATTRACTOR)
# =============================================================================

class TestLinear3DSinkField:
    """Test suite for linear 3D sink field v(x) = -A x."""

    @pytest.fixture
    def sink(self):
        return LinearSinkField(l1=2.5, l2=1.8, l3=0.9)

    def test_sink_critical_point_at_origin(self, sink):
        v_origin = sink.velocity(np.zeros(3))
        np.testing.assert_array_equal(v_origin, np.zeros(3))

    def test_sink_all_negative_eigenvalues(self, sink):
        evals_ana = sink.analytical_eigenvalues()
        assert np.all(evals_ana < 0.0), f"All sink eigenvalues must be strictly negative, got {evals_ana}"
        
        J_ana = sink.analytical_jacobian(np.zeros(3))
        evals_jacobi, _ = jacobi_eigenvalue_diagonalization(J_ana)
        np.testing.assert_allclose(evals_jacobi, np.sort(evals_ana)[::-1], atol=1e-12)

    def test_sink_strictly_negative_divergence(self, sink):
        expected_div = -(2.5 + 1.8 + 0.9)
        for pt in [np.zeros(3), np.array([10.0, -20.0, 30.0]), np.array([-50.0, 50.0, -50.0])]:
            ana_div = sink.analytical_divergence(pt)
            num_div = compute_numerical_divergence(sink.velocity, pt, h=1e-4)
            assert ana_div == pytest.approx(expected_div, 1e-12)
            assert num_div == pytest.approx(expected_div, 1e-7)
            assert num_div < 0.0

    def test_sink_zero_vorticity(self, sink):
        pt = np.array([12.0, -34.0, 56.0])
        num_curl = compute_numerical_vorticity(sink.velocity, pt, h=1e-4)
        np.testing.assert_allclose(num_curl, np.zeros(3), atol=1e-10)

    def test_sink_newton_raphson_convergence_to_origin(self, sink):
        x_init = np.array([845.0, -1290.0, 342.0])
        x_root, converged, iters = newton_raphson_critical_point_3d(sink.velocity, x_init)
        assert converged is True
        assert iters <= 5
        np.testing.assert_allclose(x_root, np.zeros(3), atol=1e-11)

    def test_sink_exponential_decay_streamlines(self, sink):
        x0 = np.array([100.0, 100.0, 100.0])
        t_end = 2.0
        t_hist, x_hist = rk45_cash_karp_integrate(sink.velocity, x0, t_end=t_end, tol=1e-9)
        x_exact = sink.analytical_streamline(x0, t_end)
        np.testing.assert_allclose(x_hist[-1], x_exact, rtol=1e-6, atol=1e-6)
        assert np.linalg.norm(x_hist[-1]) < np.linalg.norm(x0)


# =============================================================================
# 5. UNIT TESTS FOR LINEAR 3D SOURCE FIELD (REPELLER)
# =============================================================================

class TestLinear3DSourceField:
    """Test suite for linear 3D source field v(x) = +A x."""

    @pytest.fixture
    def source(self):
        return LinearSourceField(l1=1.2, l2=1.5, l3=2.1)

    def test_source_critical_point_at_origin(self, source):
        v_origin = source.velocity(np.zeros(3))
        np.testing.assert_array_equal(v_origin, np.zeros(3))

    def test_source_all_positive_eigenvalues(self, source):
        evals_ana = source.analytical_eigenvalues()
        assert np.all(evals_ana > 0.0), f"All source eigenvalues must be strictly positive, got {evals_ana}"
        
        J_ana = source.analytical_jacobian(np.zeros(3))
        evals_jacobi, _ = jacobi_eigenvalue_diagonalization(J_ana)
        np.testing.assert_allclose(evals_jacobi, np.sort(evals_ana)[::-1], atol=1e-12)

    def test_source_strictly_positive_divergence(self, source):
        expected_div = 1.2 + 1.5 + 2.1
        for pt in [np.zeros(3), np.array([5.0, 5.0, 5.0])]:
            ana_div = source.analytical_divergence(pt)
            num_div = compute_numerical_divergence(source.velocity, pt, h=1e-4)
            assert ana_div == pytest.approx(expected_div, 1e-12)
            assert num_div == pytest.approx(expected_div, 1e-7)
            assert num_div > 0.0

    def test_source_zero_vorticity(self, source):
        pt = np.array([-7.0, 14.0, -21.0])
        num_curl = compute_numerical_vorticity(source.velocity, pt, h=1e-4)
        np.testing.assert_allclose(num_curl, np.zeros(3), atol=1e-10)

    def test_source_exponential_growth_streamlines(self, source):
        x0 = np.array([1.0, 1.0, 1.0])
        t_end = 1.5
        t_hist, x_hist = rk45_cash_karp_integrate(source.velocity, x0, t_end=t_end, tol=1e-9)
        x_exact = source.analytical_streamline(x0, t_end)
        np.testing.assert_allclose(x_hist[-1], x_exact, rtol=1e-6, atol=1e-6)
        assert np.linalg.norm(x_hist[-1]) > np.linalg.norm(x0)


# =============================================================================
# 6. UNIT TESTS FOR SADDLE-FILAMENT FIELD
# =============================================================================

class TestSaddleFilamentField:
    """Test suite for Saddle-Filament field (two negative, one positive eigenvalue)."""

    @pytest.fixture
    def filament(self):
        return SaddleFilamentField(l1=1.0, l2=1.4, l3=0.8)

    def test_filament_eigenvalue_signature(self, filament):
        evals = filament.analytical_eigenvalues()
        neg_evals = evals[evals < 0]
        pos_evals = evals[evals > 0]
        assert len(neg_evals) == 2, "Saddle filament must have exactly 2 negative eigenvalues."
        assert len(pos_evals) == 1, "Saddle filament must have exactly 1 positive eigenvalue."

    def test_filament_jacobian_and_diagonalization(self, filament):
        J_ana = filament.analytical_jacobian(np.zeros(3))
        evals_jacobi, _ = jacobi_eigenvalue_diagonalization(J_ana)
        expected = np.array([0.8, -1.0, -1.4], dtype=np.float64)
        np.testing.assert_allclose(evals_jacobi, expected, atol=1e-12)

    def test_filament_xy_contraction_and_z_expansion(self, filament):
        x0 = np.array([10.0, 10.0, 1.0])
        t_end = 2.0
        t_hist, x_hist = rk45_cash_karp_integrate(filament.velocity, x0, t_end=t_end, tol=1e-9)
        x_final = x_hist[-1]
        
        # x and y should contract exponentially
        assert abs(x_final[0]) < abs(x0[0])
        assert abs(x_final[1]) < abs(x0[1])
        # z should expand exponentially along the filament spine
        assert abs(x_final[2]) > abs(x0[2])

    def test_filament_zero_curl(self, filament):
        pt = np.array([25.0, -15.0, 30.0])
        num_curl = compute_numerical_vorticity(filament.velocity, pt, h=1e-4)
        np.testing.assert_allclose(num_curl, np.zeros(3), atol=1e-10)


# =============================================================================
# 7. UNIT TESTS FOR SADDLE-WALL FIELD
# =============================================================================

class TestSaddleWallField:
    """Test suite for Saddle-Wall field (one negative, two positive eigenvalues)."""

    @pytest.fixture
    def wall(self):
        return SaddleWallField(l1=2.0, l2=0.8, l3=1.1)

    def test_wall_eigenvalue_signature(self, wall):
        evals = wall.analytical_eigenvalues()
        neg_evals = evals[evals < 0]
        pos_evals = evals[evals > 0]
        assert len(neg_evals) == 1, "Saddle wall must have exactly 1 negative eigenvalue."
        assert len(pos_evals) == 2, "Saddle wall must have exactly 2 positive eigenvalues."

    def test_wall_jacobian_and_diagonalization(self, wall):
        J_ana = wall.analytical_jacobian(np.zeros(3))
        evals_jacobi, _ = jacobi_eigenvalue_diagonalization(J_ana)
        expected = np.array([1.1, 0.8, -2.0], dtype=np.float64)
        np.testing.assert_allclose(evals_jacobi, expected, atol=1e-12)

    def test_wall_x_compression_and_yz_sheet_expansion(self, wall):
        x0 = np.array([10.0, 1.0, 1.0])
        t_end = 2.0
        t_hist, x_hist = rk45_cash_karp_integrate(wall.velocity, x0, t_end=t_end, tol=1e-9)
        x_final = x_hist[-1]
        
        # x contracts onto the wall plane (x=0)
        assert abs(x_final[0]) < abs(x0[0])
        # y and z expand along the 2D wall sheet
        assert abs(x_final[1]) > abs(x0[1])
        assert abs(x_final[2]) > abs(x0[2])

    def test_wall_zero_curl(self, wall):
        pt = np.array([-40.0, 50.0, -10.0])
        num_curl = compute_numerical_vorticity(wall.velocity, pt, h=1e-4)
        np.testing.assert_allclose(num_curl, np.zeros(3), atol=1e-10)


# =============================================================================
# 8. UNIT TESTS FOR PURE ROTATIONAL VORTEX FIELD
# =============================================================================

class TestRotationalVortexField:
    """Test suite for pure rotational vortex fields."""

    @pytest.fixture
    def vortex_z(self):
        return RotationalVortexField(omega_vector=(0.0, 0.0, 3.0))

    @pytest.fixture
    def vortex_arbitrary(self):
        return RotationalVortexField(omega_vector=(1.0, -2.0, 2.0))

    def test_vortex_zero_divergence_everywhere(self, vortex_z, vortex_arbitrary):
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([10.0, -20.0, 30.0]),
            np.array([-50.0, 50.0, -50.0])
        ]
        for vtx in [vortex_z, vortex_arbitrary]:
            for pt in points:
                num_div = compute_numerical_divergence(vtx.velocity, pt, h=1e-4)
                assert abs(num_div) < 1e-10, f"Vortex field must be strictly divergence-free, got {num_div}"

    def test_vortex_analytical_curl_matches_numerical_vorticity(self, vortex_z, vortex_arbitrary):
        points = [
            np.array([5.0, 5.0, 5.0]),
            np.array([-15.0, 25.0, 35.0])
        ]
        for vtx in [vortex_z, vortex_arbitrary]:
            for pt in points:
                ana_curl = vtx.analytical_vorticity(pt)
                num_curl = compute_numerical_vorticity(vtx.velocity, pt, h=1e-4)
                np.testing.assert_allclose(num_curl, ana_curl, rtol=1e-6, atol=1e-6)

    def test_vortex_circulation_stokes_theorem(self, vortex_z):
        """
        Verify Stokes' theorem: contour line integral oint v . dl = iint (curl v) . dA = 2 * omega_0 * pi * R^2.
        """
        radius = 5.0
        n_samples = 1000
        d_theta = 2.0 * math.pi / n_samples
        circulation = 0.0
        
        for i in range(n_samples):
            theta = i * d_theta
            # Path position: x = R*cos(th), y = R*sin(th), z = 0
            pos = np.array([radius * math.cos(theta), radius * math.sin(theta), 0.0])
            # Differential displacement: dl = (-R*sin(th), R*cos(th), 0) * d_theta
            dl = np.array([-radius * math.sin(theta), radius * math.cos(theta), 0.0]) * d_theta
            v = vortex_z.velocity(pos)
            circulation += float(np.dot(v, dl))
            
        expected_circulation = 2.0 * 3.0 * math.pi * (radius ** 2)  # curl_z = 2*3.0, area = pi*R^2
        assert circulation == pytest.approx(expected_circulation, rel=1e-4)

    def test_vortex_closed_circular_streamlines(self, vortex_z):
        x0 = np.array([10.0, 0.0, 0.0])
        period = 2.0 * math.pi / 3.0
        t_hist, x_hist = rk45_cash_karp_integrate(vortex_z.velocity, x0, t_end=period, tol=1e-10)
        # After 1 full period, streamline returns to initial seed x0
        np.testing.assert_allclose(x_hist[-1], x0, rtol=1e-5, atol=1e-5)


# =============================================================================
# 9. UNIT TESTS FOR 3D PLUMMER GRAVITATIONAL CLUSTER MODEL
# =============================================================================

class Test3DPlummerGravitationalCluster:
    """Test suite for 3D Plummer gravitational cluster model."""

    @pytest.fixture
    def cluster(self):
        return PlummerGravitationalCluster(center=(1000.0, -2000.0, 500.0), GM=5.0e7, epsilon=800.0)

    def test_plummer_center_convergence_no_singularity(self, cluster):
        # Velocity at exact center must be strictly zero (no step discontinuity or 1/r singularity)
        v_center = cluster.velocity(cluster.center)
        np.testing.assert_array_equal(v_center, np.zeros(3))
        
        # Test smooth continuous convergence as delta -> 0
        deltas = [100.0, 10.0, 1.0, 0.1, 0.001]
        v_norms = []
        for d in deltas:
            pt = cluster.center + np.array([d, 0.0, 0.0])
            v_norms.append(np.linalg.norm(cluster.velocity(pt)))
            
        # Verify strict monotonic decrease to 0
        for i in range(len(v_norms) - 1):
            assert v_norms[i] > v_norms[i+1], "Velocity magnitude must decrease smoothly as r -> 0."

    def test_plummer_peak_inflow_velocity_radius(self, cluster):
        r_peak = cluster.peak_radius()
        assert r_peak == pytest.approx(800.0 / math.sqrt(2.0), rel=1e-9)
        
        v_peak_theory = cluster.peak_velocity_magnitude()
        pt_peak = cluster.center + np.array([r_peak, 0.0, 0.0])
        v_peak_actual = np.linalg.norm(cluster.velocity(pt_peak))
        assert v_peak_actual == pytest.approx(v_peak_theory, rel=1e-9)
        
        # Test radius slightly inside and slightly outside peak
        pt_inside = cluster.center + np.array([r_peak * 0.8, 0.0, 0.0])
        pt_outside = cluster.center + np.array([r_peak * 1.2, 0.0, 0.0])
        assert np.linalg.norm(cluster.velocity(pt_inside)) < v_peak_actual
        assert np.linalg.norm(cluster.velocity(pt_outside)) < v_peak_actual

    def test_plummer_analytical_divergence_strictly_negative(self, cluster):
        points = [
            cluster.center,
            cluster.center + np.array([400.0, 400.0, 400.0]),
            cluster.center + np.array([2000.0, -1000.0, 3000.0])
        ]
        for pt in points:
            ana_div = cluster.analytical_divergence(pt)
            num_div = compute_numerical_divergence(cluster.velocity, pt, h=1e-2)
            assert ana_div < 0.0, "Plummer divergence must be strictly negative everywhere (inflow)."
            assert num_div == pytest.approx(ana_div, rel=1e-4)

    def test_plummer_analytical_vorticity_is_zero(self, cluster):
        pt = cluster.center + np.array([500.0, -300.0, 700.0])
        ana_curl = cluster.analytical_vorticity(pt)
        num_curl = compute_numerical_vorticity(cluster.velocity, pt, h=1e-2)
        np.testing.assert_array_equal(ana_curl, np.zeros(3))
        np.testing.assert_allclose(num_curl, np.zeros(3), atol=1e-8)

    def test_plummer_jacobian_matrix_and_center_degeneracy(self, cluster):
        # At cluster center, J = -(GM / eps^3) * I_3
        J_center = cluster.analytical_jacobian(cluster.center)
        expected_diag = -cluster.GM / (cluster.eps ** 3)
        np.testing.assert_allclose(J_center, expected_diag * np.eye(3), atol=1e-12)
        
        # Test general point against numerical finite difference Jacobian
        pt = cluster.center + np.array([600.0, -400.0, 300.0])
        J_ana = cluster.analytical_jacobian(pt)
        J_num = compute_numerical_jacobian(cluster.velocity, pt, h=1e-3, stencil_order=4)
        np.testing.assert_allclose(J_num, J_ana, rtol=1e-4, atol=1e-4)

    def test_plummer_streamline_inflow_into_cluster_core(self, cluster):
        x0 = cluster.center + np.array([3000.0, 2000.0, -1500.0])
        t_end = 50.0
        t_hist, x_hist = rk45_cash_karp_integrate(cluster.velocity, x0, t_end=t_end, tol=1e-8)
        
        # Distance to center must monotonically decrease
        dists = [np.linalg.norm(x - cluster.center) for x in x_hist]
        assert dists[-1] < dists[0], "Streamline must converge towards cluster core."
        for i in range(len(dists) - 1):
            assert dists[i+1] <= dists[i] + 1e-12, "Streamline inflow distance must decrease monotonically."


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
