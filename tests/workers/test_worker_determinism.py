"""
Multi-Threaded Worker Determinism Verification Test Suite
=========================================================
Rigorous mathematical verification that parallel/multi-threaded cosmological computation
workers produce bitwise and machine-precision identical numerical results compared to serial
reference executions.

Covers:
  1. Multi-threaded Streamline Integration (Batch RK45 Cash-Karp integration).
  2. Multi-threaded 3D Differential Tensor Grid Evaluation (Divergence, Vorticity, Jacobian).
  3. Multi-threaded Jacobi Eigenvalue Diagonalization across spatial tensor fields.
  4. Multi-threaded Cosmic Watershed Basin Segmentation & Voronoi partition.
  5. Multi-threaded Bulk Flow Minimum-Variance Multipole Estimation.
  6. Chunk-size invariance, thread safety, and concurrency repeatability.
"""

import math
import hashlib
import numpy as np
import pytest
from typing import List, Tuple, Dict, Any, Callable
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor


# =============================================================================
# 1. CORE NUMERICAL COMPUTATION KERNELS FOR WORKERS
# =============================================================================

def multi_cluster_velocity(x: np.ndarray) -> np.ndarray:
    """
    Analytic multi-cluster cosmological velocity field (Shapley, Great Attractor, Dipole Repeller).
    v(x) = sum_i K_i * (x - x_i) / (|x - x_i|^2 + eps_i^2)^(3/2)
    """
    clusters = [
        # (x, y, z, GM_strength, eps)
        (np.array([7200.0, -8600.0, -2400.0]), -2.3e8, 1000.0),  # Shapley Core (attractor)
        (np.array([-3500.0, 1500.0, -800.0]),  -1.2e8, 800.0),   # Great Attractor (attractor)
        (np.array([-10000.0, 10000.0, 12000.0]), +1.5e8, 1200.0) # Dipole Repeller (repeller)
    ]
    v = np.zeros(3, dtype=np.float64)
    x_arr = np.asarray(x, dtype=np.float64)
    for c_pos, K, eps in clusters:
        r_vec = x_arr - c_pos
        r_sq = np.dot(r_vec, r_vec)
        denom = (r_sq + eps * eps) ** 1.5
        v += (K / denom) * r_vec
    return v


def integrate_single_streamline(seed: np.ndarray,
                                t_end: float = 30.0,
                                h_init: float = 0.5,
                                tol: float = 1e-7,
                                max_steps: int = 1000) -> np.ndarray:
    """
    Cash-Karp RK45 integration of a single streamline. Returns trajectory array.
    """
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
    x = np.asarray(seed, dtype=np.float64).copy()
    h = h_init
    traj = [x.copy()]
    
    steps = 0
    while t < t_end and steps < max_steps:
        steps += 1
        if t + h > t_end:
            h = t_end - t
            
        k = np.zeros((6, 3), dtype=np.float64)
        k[0] = multi_cluster_velocity(x)
        for i in range(1, 6):
            x_stage = x.copy()
            for j in range(i):
                x_stage += h * a[i][j] * k[j]
            k[i] = multi_cluster_velocity(x_stage)
            
        x5 = x + h * np.sum(b5[:, None] * k, axis=0)
        err = h * np.sum(e_err[:, None] * k, axis=0)
        err_norm = np.linalg.norm(err) / (tol * (1.0 + np.linalg.norm(x)))
        
        if err_norm <= 1.0 or h < 1e-10:
            t += h
            x = x5
            traj.append(x.copy())
            
        if err_norm > 0:
            scale = 0.9 * (err_norm ** -0.2)
            scale = max(0.1, min(5.0, scale))
            h = h * scale
        else:
            h = h * 2.0
            
    return np.array(traj)


def compute_grid_voxel_tensor(point: np.ndarray, h: float = 1e-3) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Computes (divergence, curl, Jacobi-diagonalized eigenvalues) at a spatial coordinate.
    """
    x = np.asarray(point, dtype=np.float64)
    J = np.zeros((3, 3), dtype=np.float64)
    
    # 4th-order finite difference Jacobian
    for j in range(3):
        e_j = np.zeros(3, dtype=np.float64)
        e_j[j] = 1.0
        f_p2 = multi_cluster_velocity(x + 2.0 * h * e_j)
        f_p1 = multi_cluster_velocity(x + 1.0 * h * e_j)
        f_m1 = multi_cluster_velocity(x - 1.0 * h * e_j)
        f_m2 = multi_cluster_velocity(x - 2.0 * h * e_j)
        J[:, j] = (-f_p2 + 8.0 * f_p1 - 8.0 * f_m1 + f_m2) / (12.0 * h)
        
    div_v = float(np.trace(J))
    curl_v = np.array([J[2, 1] - J[1, 2], J[0, 2] - J[2, 0], J[1, 0] - J[0, 1]], dtype=np.float64)
    
    # Rate of strain / symmetric tidal tensor: S = 0.5 * (J + J.T)
    S = 0.5 * (J + J.T)
    
    # Exact Jacobi eigenvalue solver
    A = S.copy()
    for _ in range(50):
        # find max off diagonal
        p, q = 0, 1
        mv = abs(A[0, 1])
        if abs(A[0, 2]) > mv:
            p, q = 0, 2
            mv = abs(A[0, 2])
        if abs(A[1, 2]) > mv:
            p, q = 1, 2
            mv = abs(A[1, 2])
        if mv < 1e-15:
            break
        theta = 0.5 * math.atan2(2.0 * A[p, q], A[q, q] - A[p, p])
        c = math.cos(theta)
        s = math.sin(theta)
        R = np.eye(3)
        R[p, p] = c
        R[q, q] = c
        R[p, q] = s
        R[q, p] = -s
        A = R.T @ A @ R
        
    evals = np.sort(np.diag(A))[::-1]
    return div_v, curl_v, evals


def watershed_classify_voxel(point: np.ndarray,
                             attractor_centers: List[np.ndarray],
                             step_size: float = 50.0,
                             max_steps: int = 200) -> int:
    """
    Watershed cosmic basin segmentation: Traces streamline along flow field
    to identify destination attractor basin index.
    """
    x = np.asarray(point, dtype=np.float64).copy()
    for _ in range(max_steps):
        # Check proximity to attractor centers
        for idx, center in enumerate(attractor_centers):
            if np.linalg.norm(x - center) < 800.0:
                return idx
        v = multi_cluster_velocity(x)
        norm_v = np.linalg.norm(v)
        if norm_v < 1e-6:
            break
        x += (v / norm_v) * step_size
        
    # Nearest neighbor fallback
    dists = [np.linalg.norm(x - c) for c in attractor_centers]
    return int(np.argmin(dists))


def compute_bulk_flow_estimator(catalog_chunk: np.ndarray) -> np.ndarray:
    """
    Minimum-variance bulk flow dipole estimator for a chunk of galaxy positions and velocities.
    catalog_chunk: Nx6 array of [x, y, z, vx, vy, vz].
    Returns (3,) bulk flow dipole vector (V_X, V_Y, V_Z).
    """
    if len(catalog_chunk) == 0:
        return np.zeros(3, dtype=np.float64)
    # Sum weighted velocities
    positions = catalog_chunk[:, :3]
    velocities = catalog_chunk[:, 3:]
    r = np.linalg.norm(positions, axis=1, keepdims=True) + 1e-6
    unit_r = positions / r
    # Line of sight peculiar velocity: v_los = v . n
    v_los = np.sum(velocities * unit_r, axis=1, keepdims=True)
    # Weighted projection tensor: A_ij = sum n_i n_j
    A = np.zeros((3, 3), dtype=np.float64)
    b = np.zeros(3, dtype=np.float64)
    for i in range(len(catalog_chunk)):
        n = unit_r[i]
        A += np.outer(n, n)
        b += n * v_los[i, 0]
    # Invert to find bulk flow
    V_bulk = np.linalg.solve(A + 1e-8 * np.eye(3), b)
    return V_bulk


# =============================================================================
# 2. WORKER DISPATCHERS: SERIAL VS MULTI-THREADED
# =============================================================================

def run_streamlines_serial(seeds: np.ndarray) -> List[np.ndarray]:
    return [integrate_single_streamline(s) for s in seeds]


def run_streamlines_threaded(seeds: np.ndarray, num_workers: int = 4, chunk_size: int = 10) -> List[np.ndarray]:
    results = [None] * len(seeds)
    def worker(idx_start, idx_end):
        for i in range(idx_start, idx_end):
            results[i] = integrate_single_streamline(seeds[i])

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for i in range(0, len(seeds), chunk_size):
            end = min(len(seeds), i + chunk_size)
            futures.append(executor.submit(worker, i, end))
        for f in futures:
            f.result()
    return results


def run_grid_tensors_serial(coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(coords)
    divs = np.zeros(n, dtype=np.float64)
    curls = np.zeros((n, 3), dtype=np.float64)
    evals = np.zeros((n, 3), dtype=np.float64)
    for i in range(n):
        d, c, e = compute_grid_voxel_tensor(coords[i])
        divs[i] = d
        curls[i] = c
        evals[i] = e
    return divs, curls, evals


def run_grid_tensors_threaded(coords: np.ndarray, num_workers: int = 4, chunk_size: int = 16) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(coords)
    divs = np.zeros(n, dtype=np.float64)
    curls = np.zeros((n, 3), dtype=np.float64)
    evals = np.zeros((n, 3), dtype=np.float64)
    
    def process_chunk(idx_start, idx_end):
        for i in range(idx_start, idx_end):
            d, c, e = compute_grid_voxel_tensor(coords[i])
            divs[i] = d
            curls[i] = c
            evals[i] = e

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for i in range(0, n, chunk_size):
            end = min(n, i + chunk_size)
            futures.append(executor.submit(process_chunk, i, end))
        for f in futures:
            f.result()
    return divs, curls, evals


# =============================================================================
# 3. UNIT TESTS FOR WORKER DETERMINISM
# =============================================================================

class TestWorkerDeterminism:
    """Multi-threaded worker determinism verification test suite."""

    @pytest.fixture
    def seed_streamlines(self):
        np.random.seed(42)
        # Generate 100 fixed seed positions
        return np.random.uniform(-5000.0, 5000.0, size=(100, 3))

    @pytest.fixture
    def grid_coordinates(self):
        # 8x8x8 coordinate grid
        lin = np.linspace(-3000.0, 3000.0, 8)
        xx, yy, zz = np.meshgrid(lin, lin, lin, indexing='ij')
        return np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

    @pytest.fixture
    def galaxy_mock_catalog(self):
        np.random.seed(1337)
        n_galaxies = 1200
        pos = np.random.normal(0.0, 4000.0, size=(n_galaxies, 3))
        vel = np.array([multi_cluster_velocity(p) for p in pos])
        # Add small thermal dispersion
        vel += np.random.normal(0.0, 50.0, size=(n_galaxies, 3))
        return np.hstack([pos, vel])

    def test_streamline_serial_vs_parallel_exact_identity(self, seed_streamlines):
        """
        Verify that multi-threaded streamline batch integration reproduces serial
        integration with 0 error across all trajectory coordinates.
        """
        serial_results = run_streamlines_serial(seed_streamlines)
        
        for n_workers in [2, 4, 8]:
            parallel_results = run_streamlines_threaded(seed_streamlines, num_workers=n_workers, chunk_size=15)
            assert len(parallel_results) == len(serial_results)
            
            for i in range(len(seed_streamlines)):
                traj_serial = serial_results[i]
                traj_parallel = parallel_results[i]
                
                assert len(traj_serial) == len(traj_parallel), f"Trajectory step count mismatch on streamline {i}"
                # Strict machine-precision identity
                np.testing.assert_array_equal(
                    traj_parallel, traj_serial,
                    err_msg=f"Parallel streamline {i} differed from serial reference (workers={n_workers})"
                )

    def test_streamline_chunk_size_invariance(self, seed_streamlines):
        """
        Verify that varying chunk sizes (1, 5, 25, 50) does not alter results.
        """
        ref_results = run_streamlines_threaded(seed_streamlines, num_workers=4, chunk_size=1)
        for chunk_size in [5, 20, 50]:
            chunked_results = run_streamlines_threaded(seed_streamlines, num_workers=4, chunk_size=chunk_size)
            for i in range(len(seed_streamlines)):
                np.testing.assert_array_equal(chunked_results[i], ref_results[i])

    def test_grid_differential_tensors_parallel_identity(self, grid_coordinates):
        """
        Verify divergence, curl, and Jacobi eigenvalues calculated across worker threads
        match serial computation bit-for-bit.
        """
        div_s, curl_s, eval_s = run_grid_tensors_serial(grid_coordinates)
        
        for n_workers in [2, 4]:
            div_p, curl_p, eval_p = run_grid_tensors_threaded(grid_coordinates, num_workers=n_workers, chunk_size=32)
            
            np.testing.assert_array_equal(div_p, div_s, err_msg="Divergence grid parallel mismatch")
            np.testing.assert_array_equal(curl_p, curl_s, err_msg="Curl grid parallel mismatch")
            np.testing.assert_array_equal(eval_p, eval_s, err_msg="Eigenvalue grid parallel mismatch")

    def test_watershed_segmentation_parallel_determinism(self, grid_coordinates):
        """
        Verify watershed basin classification returns identical basin assignments across threads.
        """
        attractor_centers = [
            np.array([7200.0, -8600.0, -2400.0]), # Shapley
            np.array([-3500.0, 1500.0, -800.0])   # GA
        ]
        
        # Serial basin assignment
        serial_labels = np.array([watershed_classify_voxel(pt, attractor_centers) for pt in grid_coordinates])
        
        # Multi-threaded basin assignment
        def parallel_watershed(coords, num_workers=4):
            labels = np.zeros(len(coords), dtype=np.int32)
            def worker(i_start, i_end):
                for i in range(i_start, i_end):
                    labels[i] = watershed_classify_voxel(coords[i], attractor_centers)
            with ThreadPoolExecutor(max_workers=num_workers) as ex:
                futures = []
                for i in range(0, len(coords), 32):
                    futures.append(ex.submit(worker, i, min(len(coords), i + 32)))
                for f in futures:
                    f.result()
            return labels

        parallel_labels = parallel_watershed(grid_coordinates, num_workers=4)
        np.testing.assert_array_equal(parallel_labels, serial_labels)
        
        # Dice and Jaccard similarity between serial and parallel segmentations must be exactly 1.0
        for basin_id in [0, 1]:
            mask_s = (serial_labels == basin_id)
            mask_p = (parallel_labels == basin_id)
            intersection = np.sum(mask_s & mask_p)
            dice = 2.0 * intersection / (np.sum(mask_s) + np.sum(mask_p))
            jaccard = intersection / np.sum(mask_s | mask_p)
            assert dice == 1.0
            assert jaccard == 1.0

    def test_bulk_flow_estimator_deterministic_reduction(self, galaxy_mock_catalog):
        """
        Verify minimum-variance bulk flow multipole estimator determinism across partitions.
        """
        # Global serial calculation
        v_bulk_serial = compute_bulk_flow_estimator(galaxy_mock_catalog)
        assert np.linalg.norm(v_bulk_serial) > 0.0
        
        # Partitioned worker reduction
        def parallel_bulk_flow(catalog, num_partitions=4):
            # Each worker computes partial normal equations (A_k, b_k)
            n_total = len(catalog)
            chunk_len = math.ceil(n_total / num_partitions)
            partials = []
            
            def worker(chunk):
                positions = chunk[:, :3]
                velocities = chunk[:, 3:]
                r = np.linalg.norm(positions, axis=1, keepdims=True) + 1e-6
                unit_r = positions / r
                v_los = np.sum(velocities * unit_r, axis=1, keepdims=True)
                A = np.zeros((3, 3), dtype=np.float64)
                b = np.zeros(3, dtype=np.float64)
                for i in range(len(chunk)):
                    n = unit_r[i]
                    A += np.outer(n, n)
                    b += n * v_los[i, 0]
                return A, b
                
            with ThreadPoolExecutor(max_workers=num_partitions) as ex:
                futures = []
                for i in range(0, n_total, chunk_len):
                    futures.append(ex.submit(worker, catalog[i:i+chunk_len]))
                for f in futures:
                    partials.append(f.result())
                    
            A_total = sum(p[0] for p in partials)
            b_total = sum(p[1] for p in partials)
            return np.linalg.solve(A_total + 1e-8 * np.eye(3), b_total)

        v_bulk_parallel = parallel_bulk_flow(galaxy_mock_catalog, num_partitions=4)
        np.testing.assert_allclose(v_bulk_parallel, v_bulk_serial, rtol=1e-13, atol=1e-13)

    def test_worker_thread_safety_and_memory_isolation(self, seed_streamlines):
        """
        Execute 10 concurrent runs and verify bitwise SHA-256 hash identity across all runs.
        """
        hashes = []
        for _ in range(10):
            res = run_streamlines_threaded(seed_streamlines, num_workers=4, chunk_size=10)
            # Concatenate all coordinates and hash bytes
            all_bytes = np.concatenate(res).tobytes()
            h = hashlib.sha256(all_bytes).hexdigest()
            hashes.append(h)
            
        assert len(set(hashes)) == 1, "Thread safety violation: Non-deterministic results detected across runs!"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
