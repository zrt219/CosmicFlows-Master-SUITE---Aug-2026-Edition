/**
 * Shared utility module for trilinear interpolation of a velocity grid.
 */

self.trilinear = (function() {
    /**
     * Compute the flat index into the Float32Array for the (ix, iy, iz) voxel's vx component.
     * @param {number} ix - X index (0 to N-1)
     * @param {number} iy - Y index (0 to N-1)
     * @param {number} iz - Z index (0 to N-1)
     * @param {number} N - Grid dimension (e.g., 128)
     * @returns {number} The flat index.
     */
    function voxelIndex(ix, iy, iz, N) {
        return ((ix * N * N) + (iy * N) + iz) * 3;
    }

    /**
     * Trilinear interpolation of velocity at a given point.
     * @param {Float32Array} grid - Velocity grid [N][N][N][3]
     * @param {number} N - Grid dimension (e.g., 128)
     * @param {number} halfExtent - Grid half-extent in km/s (e.g., 16000)
     * @param {number} sgx - Supergalactic X coordinate
     * @param {number} sgy - Supergalactic Y coordinate
     * @param {number} sgz - Supergalactic Z coordinate
     * @returns {{vx: number, vy: number, vz: number}} Interpolated velocity
     */
    function trilinearInterpolate(grid, N, halfExtent, sgx, sgy, sgz) {
        // Map coordinates from [-halfExtent, halfExtent] to [0, N-1]
        let fx = (sgx + halfExtent) / (2 * halfExtent) * (N - 1);
        let fy = (sgy + halfExtent) / (2 * halfExtent) * (N - 1);
        let fz = (sgz + halfExtent) / (2 * halfExtent) * (N - 1);

        // Clamp to grid boundaries
        if (fx < 0) fx = 0; if (fx > N - 1) fx = N - 1;
        if (fy < 0) fy = 0; if (fy > N - 1) fy = N - 1;
        if (fz < 0) fz = 0; if (fz > N - 1) fz = N - 1;

        let ix0 = Math.floor(fx);
        let iy0 = Math.floor(fy);
        let iz0 = Math.floor(fz);

        let ix1 = Math.min(ix0 + 1, N - 1);
        let iy1 = Math.min(iy0 + 1, N - 1);
        let iz1 = Math.min(iz0 + 1, N - 1);

        let tx = fx - ix0;
        let ty = fy - iy0;
        let tz = fz - iz0;

        let inv_tx = 1.0 - tx;
        let inv_ty = 1.0 - ty;
        let inv_tz = 1.0 - tz;

        let w000 = inv_tx * inv_ty * inv_tz;
        let w100 = tx * inv_ty * inv_tz;
        let w010 = inv_tx * ty * inv_tz;
        let w110 = tx * ty * inv_tz;
        let w001 = inv_tx * inv_ty * tz;
        let w101 = tx * inv_ty * tz;
        let w011 = inv_tx * ty * tz;
        let w111 = tx * ty * tz;

        let i000 = voxelIndex(ix0, iy0, iz0, N);
        let i100 = voxelIndex(ix1, iy0, iz0, N);
        let i010 = voxelIndex(ix0, iy1, iz0, N);
        let i110 = voxelIndex(ix1, iy1, iz0, N);
        let i001 = voxelIndex(ix0, iy0, iz1, N);
        let i101 = voxelIndex(ix1, iy0, iz1, N);
        let i011 = voxelIndex(ix0, iy1, iz1, N);
        let i111 = voxelIndex(ix1, iy1, iz1, N);

        let vx = w000 * grid[i000]     + w100 * grid[i100]     + w010 * grid[i010]     + w110 * grid[i110]     + w001 * grid[i001]     + w101 * grid[i101]     + w011 * grid[i011]     + w111 * grid[i111];
        let vy = w000 * grid[i000 + 1] + w100 * grid[i100 + 1] + w010 * grid[i010 + 1] + w110 * grid[i110 + 1] + w001 * grid[i001 + 1] + w101 * grid[i101 + 1] + w011 * grid[i011 + 1] + w111 * grid[i111 + 1];
        let vz = w000 * grid[i000 + 2] + w100 * grid[i100 + 2] + w010 * grid[i010 + 2] + w110 * grid[i110 + 2] + w001 * grid[i001 + 2] + w101 * grid[i101 + 2] + w011 * grid[i011 + 2] + w111 * grid[i111 + 2];

        return { vx: vx, vy: vy, vz: vz };
    }

    /**
     * Compute the velocity gradient (Jacobian matrix) at a given point using 6-point central finite differencing.
     * @param {Float32Array} grid - Velocity grid [N][N][N][3]
     * @param {number} N - Grid dimension (e.g., 128)
     * @param {number} halfExtent - Grid half-extent in km/s (e.g., 16000)
     * @param {number} sgx - Supergalactic X coordinate
     * @param {number} sgy - Supergalactic Y coordinate
     * @param {number} sgz - Supergalactic Z coordinate
     * @param {number} [delta=250] - Step size for finite difference in km/s
     * @returns {number[][]} 3x3 matrix [[dvx/dx, dvx/dy, dvx/dz], [dvy/dx, dvy/dy, dvy/dz], [dvz/dx, dvz/dy, dvz/dz]]
     */
    function trilinearGradient(grid, N, halfExtent, sgx, sgy, sgz, delta) {
        if (delta === undefined) {
            delta = 250;
        }

        let v_px = trilinearInterpolate(grid, N, halfExtent, sgx + delta, sgy, sgz);
        let v_nx = trilinearInterpolate(grid, N, halfExtent, sgx - delta, sgy, sgz);
        let v_py = trilinearInterpolate(grid, N, halfExtent, sgx, sgy + delta, sgz);
        let v_ny = trilinearInterpolate(grid, N, halfExtent, sgx, sgy - delta, sgz);
        let v_pz = trilinearInterpolate(grid, N, halfExtent, sgx, sgy, sgz + delta);
        let v_nz = trilinearInterpolate(grid, N, halfExtent, sgx, sgy, sgz - delta);

        let den = 2.0 * delta;

        let dvx_dx = (v_px.vx - v_nx.vx) / den;
        let dvy_dx = (v_px.vy - v_nx.vy) / den;
        let dvz_dx = (v_px.vz - v_nx.vz) / den;

        let dvx_dy = (v_py.vx - v_ny.vx) / den;
        let dvy_dy = (v_py.vy - v_ny.vy) / den;
        let dvz_dy = (v_py.vz - v_ny.vz) / den;

        let dvx_dz = (v_pz.vx - v_nz.vx) / den;
        let dvy_dz = (v_pz.vy - v_nz.vy) / den;
        let dvz_dz = (v_pz.vz - v_nz.vz) / den;

        return [
            [dvx_dx, dvx_dy, dvx_dz],
            [dvy_dx, dvy_dy, dvy_dz],
            [dvz_dx, dvz_dy, dvz_dz]
        ];
    }

    return {
        voxelIndex: voxelIndex,
        trilinearInterpolate: trilinearInterpolate,
        trilinearGradient: trilinearGradient
    };
})();
