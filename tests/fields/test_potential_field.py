"""
tests/fields/test_potential_field.py
Rigorous automated Pytest suite for Cosmological Potential & Gravitational Clustering Engine.

Tests:
1. Spectral 3D FFT Poisson Solver: Harmonic standing wave exact potential and discrete Laplacian residual.
2. Analytical Plummer Sphere: Potential Phi(r) and acceleration g(r) accuracy, isolated Green's FFT convolution solver.
3. Analytical Hernquist Profile: Potential and acceleration matching.
4. Cosmological Tidal Tensor & Jacobi Diagonalizer: Eigensolver orthonormality, eigenvalue ordering, trace invariance.
5. Gravitational Potential Energy Density: Local energy density and total integrated gravitational energy W.
6. Multigrid / SOR Poisson Solver: Convergence and agreement with spectral solver.
7. Acceleration vector buffers and finite-difference gradient stencils.
"""

import json
import math
import os
import subprocess
import numpy as np
import pytest

NODE_EXEC = "node"
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

def run_node_snippet(code: str) -> dict:
    """Executes an inline ES6 Node.js snippet and returns parsed JSON output."""
    wrapped_code = f"""
    {code}
    """
    proc = subprocess.run(
        [NODE_EXEC, "--input-type=module", "-e", wrapped_code],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Node execution failed with code {proc.returncode}:\nSTDERR:\n{proc.stderr}\nSTDOUT:\n{proc.stdout}")
    return json.loads(proc.stdout)


# ============================================================================
# 1. SPECTRAL 3D FFT POISSON SOLVER TESTS
# ============================================================================

class TestSpectralPoissonSolver:
    def test_harmonic_standing_wave_exact_solution(self):
        """Verify 3D FFT Poisson solver reproduces exact analytical solution for harmonic standing wave."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { PotentialField } from './src/fields/potential_field.js';

        const N = 32;
        const L = 100.0;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [L, L, L]
        });

        const kx = (2.0 * Math.PI) / L;
        const ky = (2.0 * Math.PI) / L;
        const kz = (2.0 * Math.PI) / L;
        const amplitude = 2.5;
        const fourPiG = 1.0;

        const harmonic = PotentialField.fromAnalyticHarmonic(grid, { kx, ky, kz, amplitude, fourPiG });
        const solved = PotentialField.solvePoissonFFT(harmonic.densityField, { fourPiG, kernel: 'continuous' });

        let maxAbsError = 0.0;
        let maxRelError = 0.0;
        let l2Error = 0.0;
        let l2Exact = 0.0;

        for (let i = 0; i < grid.totalCells; i++) {
          const exactVal = harmonic.potentialField.phi[i];
          const solvedVal = solved.phi[i];
          const diff = Math.abs(solvedVal - exactVal);

          if (diff > maxAbsError) maxAbsError = diff;
          if (Math.abs(exactVal) > 1e-6) {
            const rel = diff / Math.abs(exactVal);
            if (rel > maxRelError) maxRelError = rel;
          }

          l2Error += diff * diff;
          l2Exact += exactVal * exactVal;
        }

        const relativeL2 = Math.sqrt(l2Error / l2Exact);

        console.log(JSON.stringify({
          maxAbsError,
          relativeL2
        }));
        """
        res = run_node_snippet(code)
        assert res["relativeL2"] < 0.15, f"Harmonic wave FFT Poisson solver relative L2 error too high: {res['relativeL2']}"

    def test_poisson_discrete_laplacian_residual(self):
        """Verify discrete Laplacian nabla^2 Phi matches 4 * pi * G * delta(x) at interior voxels."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { PotentialField } from './src/fields/potential_field.js';

        const N = 32;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100]
        });

        const kx = (4.0 * Math.PI) / 100.0;
        const ky = (4.0 * Math.PI) / 100.0;
        const kz = (4.0 * Math.PI) / 100.0;
        const harmonic = PotentialField.fromAnalyticHarmonic(grid, { kx, ky, kz, amplitude: 1.0, fourPiG: 1.5 });
        const solved = PotentialField.solvePoissonFFT(harmonic.densityField, { fourPiG: 1.5, kernel: 'discrete' });

        let maxResidual = 0.0;
        for (let iz = 2; iz < N - 2; iz++) {
          for (let iy = 2; iy < N - 2; iy++) {
            for (let ix = 2; ix < N - 2; ix++) {
              const lap = solved.laplacianAtVoxel(ix, iy, iz);
              const idx = grid.index(ix, iy, iz);
              const source = 1.5 * harmonic.densityField.delta[idx];
              const res = Math.abs(lap - source);
              if (res > maxResidual) maxResidual = res;
            }
          }
        }

        console.log(JSON.stringify({ maxResidual }));
        """
        res = run_node_snippet(code)
        assert res["maxResidual"] < 1e-4, f"Discrete Laplacian Poisson residual too high: {res['maxResidual']}"


# ============================================================================
# 2. ANALYTICAL PLUMMER & HERNQUIST SPHERE TESTS
# ============================================================================

class TestAnalyticCosmicSpheres:
    def test_plummer_potential_and_acceleration_accuracy(self):
        """Verify Plummer profile potential and acceleration evaluation against exact analytic formulas."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { PotentialField } from './src/fields/potential_field.js';

        const N = 64;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-40, -40, -40],
          boxSize: [80, 80, 80]
        });

        const M = 1.0e14;
        const b = 8.0;
        const G = 4.30091e-9;

        const plummer = PotentialField.fromAnalyticPlummer(grid, { mass: M, scaleRadius: b, G, center: [0, 0, 0] });

        // Sample test points across resolved radial range (r >= dx)
        const testRadii = [4.0, 8.0, 12.0, 20.0, 30.0];
        const errors = [];

        for (const r of testRadii) {
          const testX = r / Math.sqrt(3.0);
          const testY = r / Math.sqrt(3.0);
          const testZ = r / Math.sqrt(3.0);

          const phiInterp = plummer.potentialField.evaluatePotential(testX, testY, testZ);
          const phiExact = plummer.exactPotential(testX, testY, testZ);
          const relPhiErr = Math.abs(phiInterp - phiExact) / Math.abs(phiExact);

          const gInterp = plummer.potentialField.evaluateAcceleration(testX, testY, testZ, '14point');
          const gExact = plummer.exactAcceleration(testX, testY, testZ);
          const gMagInterp = Math.sqrt(gInterp[0]**2 + gInterp[1]**2 + gInterp[2]**2);
          const gMagExact = Math.sqrt(gExact[0]**2 + gExact[1]**2 + gExact[2]**2);
          const relGErr = Math.abs(gMagInterp - gMagExact) / Math.abs(gMagExact);

          errors.push({ r, relPhiErr, relGErr });
        }

        console.log(JSON.stringify({ errors }));
        """
        res = run_node_snippet(code)
        for err in res["errors"]:
            assert err["relPhiErr"] < 0.05, f"Plummer potential interpolation error at r={err['r']} exceeds 5%: {err['relPhiErr']}"
            assert err["relGErr"] < 0.25, f"Plummer acceleration interpolation error at r={err['r']} exceeds 25%: {err['relGErr']}"

    def test_hernquist_potential_and_acceleration(self):
        """Verify Hernquist profile potential and acceleration fields."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { PotentialField } from './src/fields/potential_field.js';

        const N = 32;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100]
        });

        const M = 2.0e14;
        const a = 10.0;
        const G = 4.30091e-9;

        const hernquist = PotentialField.fromAnalyticHernquist(grid, { mass: M, scaleRadius: a, G, center: [0, 0, 0] });

        const x = 12.0, y = 15.0, z = -10.0;
        const phiInterp = hernquist.potentialField.evaluatePotential(x, y, z);
        const phiExact = hernquist.exactPotential(x, y, z);
        const relPhiErr = Math.abs(phiInterp - phiExact) / Math.abs(phiExact);

        const gInterp = hernquist.potentialField.evaluateAcceleration(x, y, z);
        const gExact = hernquist.exactAcceleration(x, y, z);

        const gErr = Math.sqrt((gInterp[0] - gExact[0])**2 + (gInterp[1] - gExact[1])**2 + (gInterp[2] - gExact[2])**2) /
                     Math.sqrt(gExact[0]**2 + gExact[1]**2 + gExact[2]**2);

        console.log(JSON.stringify({ relPhiErr, gErr }));
        """
        res = run_node_snippet(code)
        assert res["relPhiErr"] < 0.01, f"Hernquist potential error: {res['relPhiErr']}"
        assert res["gErr"] < 0.02, f"Hernquist acceleration error: {res['gErr']}"

    def test_isolated_greens_function_plummer_solver(self):
        """Verify Isolated Green's function FFT convolution solver matches analytic Plummer potential."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { PotentialField } from './src/fields/potential_field.js';

        const N = 32;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100]
        });

        const M = 1.0e14;
        const b = 10.0;
        const G = 4.30091e-9;

        const plummer = PotentialField.fromAnalyticPlummer(grid, { mass: M, scaleRadius: b, G, center: [0, 0, 0] });
        const isolated = PotentialField.solvePoissonIsolated(plummer.densityField, { G, rhoBar: 1.0, softening: b * 0.1 });

        // Compare potential in interior core region (r < 30 Mpc)
        let l2Diff = 0.0;
        let l2Exact = 0.0;
        let sampleCount = 0;

        for (let iz = 4; iz < N - 4; iz++) {
          for (let iy = 4; iy < N - 4; iy++) {
            for (let ix = 4; ix < N - 4; ix++) {
              const [x, y, z] = grid.gridIndexToCoord(ix, iy, iz);
              const r = Math.sqrt(x*x + y*y + z*z);
              if (r > 30.0) continue;

              const exactVal = plummer.exactPotential(x, y, z);
              const solvedVal = isolated.phi[grid.index(ix, iy, iz)];
              const diff = solvedVal - exactVal;

              l2Diff += diff * diff;
              l2Exact += exactVal * exactVal;
              sampleCount++;
            }
          }
        }

        const relativeL2 = Math.sqrt(l2Diff / l2Exact);
        console.log(JSON.stringify({ relativeL2, sampleCount }));
        """
        res = run_node_snippet(code)
        assert res["relativeL2"] < 0.05, f"Isolated Green's solver relative L2 error: {res['relativeL2']}"


# ============================================================================
# 3. TIDAL TENSOR & JACOBI DIAGONALIZER TESTS
# ============================================================================

class TestTidalTensorAndEigenvalues:
    def test_tidal_tensor_symmetry_and_jacobi_eigensolver(self):
        """Verify Tidal tensor is symmetric, Jacobi solver eigenvalues are sorted and eigenvectors orthonormal."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { PotentialField } from './src/fields/potential_field.js';

        const N = 32;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100]
        });

        const plummer = PotentialField.fromAnalyticPlummer(grid, { mass: 1e14, scaleRadius: 10.0, center: [0, 0, 0] });

        // Sample voxel away from center
        const ix = 12, iy = 14, iz = 18;
        const H = plummer.potentialField.computeHessianAtVoxel(ix, iy, iz);
        const symmErr01 = Math.abs(H[0][1] - H[1][0]);
        const symmErr02 = Math.abs(H[0][2] - H[2][0]);
        const symmErr12 = Math.abs(H[1][2] - H[2][1]);

        const eigenFull = plummer.potentialField.computeTidalEigenvaluesAtVoxel(ix, iy, iz, false);
        const [l1, l2, l3] = eigenFull.eigenvalues;
        const [e1, e2, e3] = eigenFull.eigenvectors;

        // Check eigenvalue order
        const isSorted = (l1 <= l2) && (l2 <= l3);

        // Check eigenvector orthonormality: e_i . e_j = delta_ij
        const dot11 = e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2];
        const dot22 = e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2];
        const dot33 = e3[0]*e3[0] + e3[1]*e3[1] + e3[2]*e3[2];
        const dot12 = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
        const dot13 = e1[0]*e3[0] + e1[1]*e3[1] + e1[2]*e3[2];
        const dot23 = e2[0]*e3[0] + e2[1]*e3[1] + e2[2]*e3[2];

        // Traceless tidal tensor test
        const eigenTraceless = plummer.potentialField.computeTidalEigenvaluesAtVoxel(ix, iy, iz, true);
        const traceTraceless = eigenTraceless.eigenvalues[0] + eigenTraceless.eigenvalues[1] + eigenTraceless.eigenvalues[2];

        console.log(JSON.stringify({
          symmErr: Math.max(symmErr01, symmErr02, symmErr12),
          isSorted,
          normErr: Math.max(Math.abs(dot11 - 1), Math.abs(dot22 - 1), Math.abs(dot33 - 1)),
          orthoErr: Math.max(Math.abs(dot12), Math.abs(dot13), Math.abs(dot23)),
          traceTraceless
        }));
        """
        res = run_node_snippet(code)
        assert res["symmErr"] < 1e-12, f"Hessian not symmetric: {res['symmErr']}"
        assert res["isSorted"] is True, "Eigenvalues not sorted in ascending order"
        assert res["normErr"] < 1e-12, f"Eigenvector not normalized: {res['normErr']}"
        assert res["orthoErr"] < 1e-12, f"Eigenvectors not orthogonal: {res['orthoErr']}"
        assert abs(res["traceTraceless"]) < 1e-12, f"Traceless tidal tensor trace non-zero: {res['traceTraceless']}"


# ============================================================================
# 4. GRAVITATIONAL POTENTIAL ENERGY DENSITY TESTS
# ============================================================================

class TestGravitationalEnergy:
    def test_gravitational_energy_integration(self):
        """Verify total gravitational potential energy W = int 0.5 * rho * Phi dV is negative and finite."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { PotentialField } from './src/fields/potential_field.js';

        const N = 32;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-50, -50, -50],
          boxSize: [100, 100, 100]
        });

        const M = 1.0e14;
        const b = 10.0;
        const G = 4.30091e-9;
        const plummer = PotentialField.fromAnalyticPlummer(grid, { mass: M, scaleRadius: b, G, center: [0, 0, 0] });

        const totalW = plummer.potentialField.computeTotalGravitationalPotentialEnergy(plummer.densityField);
        const energyBuffer = plummer.potentialField.computeEnergyDensityBuffer(plummer.densityField);

        // All local energy densities for attractive gravity should be <= 0
        let maxLocalW = -Infinity;
        for (let i = 0; i < energyBuffer.length; i++) {
          if (energyBuffer[i] > maxLocalW) maxLocalW = energyBuffer[i];
        }

        // Theoretical Plummer binding energy W_exact = -3 * pi / 32 * G * M^2 / b
        const wExact = -(3.0 * Math.PI / 32.0) * (G * M * M) / b;
        const relWErr = Math.abs(totalW - wExact) / Math.abs(wExact);

        console.log(JSON.stringify({
          totalW,
          wExact,
          relWErr,
          maxLocalW
        }));
        """
        res = run_node_snippet(code)
        assert res["totalW"] < 0, f"Total gravitational energy must be negative: {res['totalW']}"
        assert res["maxLocalW"] <= 0, f"Local energy density must be non-positive: {res['maxLocalW']}"
        # Box truncation at 5 * b will capture > 90% of Plummer energy
        assert res["relWErr"] < 0.15, f"Integrated Plummer energy relative error: {res['relWErr']}"


# ============================================================================
# 5. MULTIGRID / SOR POISSON SOLVER TESTS
# ============================================================================

class TestMultigridPoissonSolver:
    def test_multigrid_sor_convergence_and_fft_agreement(self):
        """Verify Multigrid/SOR Poisson solver converges and matches Spectral FFT solver."""
        code = """
        import { GridIndexer, BoundaryMode } from './src/fields/grid_indexer.js';
        import { PotentialField } from './src/fields/potential_field.js';

        const N = 16;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [50, 50, 50],
          boundaryMode: BoundaryMode.PERIODIC
        });

        const harmonic = PotentialField.fromAnalyticHarmonic(grid, { amplitude: 1.0, fourPiG: 1.0 });
        const solvedFFT = PotentialField.solvePoissonFFT(harmonic.densityField, { fourPiG: 1.0, kernel: 'discrete' });
        const solvedSOR = PotentialField.solvePoissonMultigrid(harmonic.densityField, { fourPiG: 1.0, maxIter: 800, tol: 1e-7, omega: 1.5 });

        // Normalize mean to compare shape
        let meanFFT = 0, meanSOR = 0;
        for (let i = 0; i < grid.totalCells; i++) {
          meanFFT += solvedFFT.phi[i];
          meanSOR += solvedSOR.phi[i];
        }
        meanFFT /= grid.totalCells;
        meanSOR /= grid.totalCells;

        let maxDiff = 0.0;
        for (let i = 0; i < grid.totalCells; i++) {
          const diff = Math.abs((solvedSOR.phi[i] - meanSOR) - (solvedFFT.phi[i] - meanFFT));
          if (diff > maxDiff) maxDiff = diff;
        }

        console.log(JSON.stringify({ maxDiff }));
        """
        res = run_node_snippet(code)
        assert res["maxDiff"] < 0.05, f"SOR vs FFT Poisson solution deviation too large: {res['maxDiff']}"


# ============================================================================
# 6. ACCELERATION BUFFERS & FINITE DIFFERENCE STENCILS
# ============================================================================

class TestAccelerationBuffers:
    def test_precomputed_acceleration_buffers(self):
        """Verify getAccelerationBuffers precomputes consistent gx, gy, gz fields."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { PotentialField } from './src/fields/potential_field.js';

        const N = 16;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-25, -25, -25],
          boxSize: [50, 50, 50]
        });

        const plummer = PotentialField.fromAnalyticPlummer(grid, { mass: 1e14, scaleRadius: 5.0 });
        const { gx, gy, gz } = plummer.potentialField.getAccelerationBuffers();

        // Verify interior voxel acceleration matching
        const ix = 8, iy = 8, iz = 10;
        const idx = grid.index(ix, iy, iz);
        const [gvx, gvy, gvz] = plummer.potentialField.computeAccelerationAtVoxel(ix, iy, iz);

        const diffX = Math.abs(gx[idx] - gvx);
        const diffY = Math.abs(gy[idx] - gvy);
        const diffZ = Math.abs(gz[idx] - gvz);

        console.log(JSON.stringify({
          bufferLength: gx.length,
          totalCells: grid.totalCells,
          diffX, diffY, diffZ
        }));
        """
        res = run_node_snippet(code)
        assert res["bufferLength"] == res["totalCells"]
        assert res["diffX"] < 30.0 and res["diffY"] < 30.0 and res["diffZ"] < 30.0
