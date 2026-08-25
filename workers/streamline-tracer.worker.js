importScripts('./trilinear.js');

// Cash-Karp Butcher Tableau Constants
const a21 = 1/5;
const a31 = 3/40, a32 = 9/40;
const a41 = 3/10, a42 = -9/10, a43 = 6/5;
const a51 = -11/54, a52 = 5/2, a53 = -70/27, a54 = 35/27;
const a61 = 1631/55296, a62 = 175/512, a63 = 575/13824, a64 = 44275/110592, a65 = 253/4096;

// 5th order weights
const b1 = 37/378, b2 = 0, b3 = 250/621, b4 = 125/594, b5 = 0, b6 = 512/1771;
// 4th order weights
const b1s = 2825/27648, b2s = 0, b3s = 18575/48384, b4s = 13525/55296, b5s = 277/14336, b6s = 1/4;

self.onmessage = function(e) {
    const data = e.data;
    if (data.type === 'trace') {
        try {
            processTrace(data);
        } catch (err) {
            self.postMessage({ type: 'error', message: err.message });
        }
    }
};

function processTrace(params) {
    const { gridBuffer, N, halfExtent, seeds, maxSteps, stepRange, tolerance, attractors } = params;
    const grid = new Float32Array(gridBuffer);
    
    const results = [];
    let completed = 0;
    
    // Helper to get velocity from trilinear.js
    // Assumes a function `trilinearInterpolate(x, y, z, grid, N, halfExtent)` is available globally
    const getVel = (x, y, z) => {
        if (typeof trilinearInterpolate === 'function') {
            const v = trilinearInterpolate(x, y, z, grid, N, halfExtent);
            return Array.isArray(v) ? v : [v.x || 0, v.y || 0, v.z || 0];
        }
        return [0, 0, 0];
    };

    const sqr = (v) => v * v;
    const mag = (vx, vy, vz) => Math.sqrt(sqr(vx) + sqr(vy) + sqr(vz));
    
    for (const seed of seeds) {
        let x = seed.sgx;
        let y = seed.sgy;
        let z = seed.sgz;
        let h = stepRange[0];
        
        const waypoints = [];
        let steps = 0;
        let maxError = 0;
        let domainExit = false;
        
        const startX = x, startY = y, startZ = z;
        
        while (steps < maxSteps) {
            if (steps % 5 === 0) {
                waypoints.push({ x, y, z });
            }
            
            // Check domain limits
            if (Math.abs(x) > halfExtent || Math.abs(y) > halfExtent || Math.abs(z) > halfExtent) {
                domainExit = true;
                break;
            }
            
            const v1 = getVel(x, y, z);
            if (mag(v1[0], v1[1], v1[2]) < 5.0) {
                break; // Convergence
            }
            
            // Compute 6 stages (k1..k6)
            const k1x = v1[0] * h, k1y = v1[1] * h, k1z = v1[2] * h;
            
            const v2 = getVel(x + a21*k1x, y + a21*k1y, z + a21*k1z);
            const k2x = v2[0] * h, k2y = v2[1] * h, k2z = v2[2] * h;
            
            const v3 = getVel(x + a31*k1x + a32*k2x, y + a31*k1y + a32*k2y, z + a31*k1z + a32*k2z);
            const k3x = v3[0] * h, k3y = v3[1] * h, k3z = v3[2] * h;
            
            const v4 = getVel(x + a41*k1x + a42*k2x + a43*k3x, y + a41*k1y + a42*k2y + a43*k3y, z + a41*k1z + a42*k2z + a43*k3z);
            const k4x = v4[0] * h, k4y = v4[1] * h, k4z = v4[2] * h;
            
            const v5 = getVel(x + a51*k1x + a52*k2x + a53*k3x + a54*k4x, y + a51*k1y + a52*k2y + a53*k3y + a54*k4y, z + a51*k1z + a52*k2z + a53*k3z + a54*k4z);
            const k5x = v5[0] * h, k5y = v5[1] * h, k5z = v5[2] * h;
            
            const v6 = getVel(x + a61*k1x + a62*k2x + a63*k3x + a64*k4x + a65*k5x, y + a61*k1y + a62*k2y + a63*k3y + a64*k4y + a65*k5y, z + a61*k1z + a62*k2z + a63*k3z + a64*k4z + a65*k5z);
            const k6x = v6[0] * h, k6y = v6[1] * h, k6z = v6[2] * h;
            
            // 5th order solution
            const dx5 = b1*k1x + b2*k2x + b3*k3x + b4*k4x + b5*k5x + b6*k6x;
            const dy5 = b1*k1y + b2*k2y + b3*k3y + b4*k4y + b5*k5y + b6*k6y;
            const dz5 = b1*k1z + b2*k2z + b3*k3z + b4*k4z + b5*k5z + b6*k6z;
            
            // 4th order solution
            const dx4 = b1s*k1x + b2s*k2x + b3s*k3x + b4s*k4x + b5s*k5x + b6s*k6x;
            const dy4 = b1s*k1y + b2s*k2y + b3s*k3y + b4s*k4y + b5s*k5y + b6s*k6y;
            const dz4 = b1s*k1z + b2s*k2z + b3s*k3z + b4s*k4z + b5s*k5z + b6s*k6z;
            
            // Error estimate: |y5 - y4|
            const err = mag(dx5 - dx4, dy5 - dy4, dz5 - dz4);
            if (err > maxError) maxError = err;
            
            // Step size adaptation & rejection
            if (err > tolerance && h > stepRange[0]) {
                h = h * Math.max(0.2, 0.84 * Math.pow(tolerance / err, 0.25));
                if (h < stepRange[0]) h = stepRange[0];
                continue; // retry step with smaller h
            }
            
            // Accept step
            x += dx5;
            y += dy5;
            z += dz5;
            
            // Calculate next step size
            let h_new = h;
            if (err > 0) {
                h_new = h * Math.min(5, Math.max(0.2, 0.84 * Math.pow(tolerance / err, 0.25)));
            } else {
                h_new = h * 5;
            }
            
            h = Math.min(stepRange[1], Math.max(stepRange[0], h_new));
            steps++;
        }
        
        // Ensure terminal point is recorded if not exactly on a multiple of 5
        if (steps % 5 !== 0) {
             waypoints.push({ x, y, z });
        }
        
        // Classify terminal position
        let nearestBasin = 'Unknown';
        let minDistSq = Infinity;
        if (attractors && attractors.length > 0) {
            for (const att of attractors) {
                const distSq = sqr(x - att.sgx) + sqr(y - att.sgy) + sqr(z - att.sgz);
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    nearestBasin = att.name;
                }
            }
        }
        
        const displacement = mag(x - startX, y - startY, z - startZ);
        
        results.push({
            id: seed.id,
            waypoints,
            terminal: { sgx: x, sgy: y, sgz: z },
            basin: nearestBasin,
            displacement,
            steps,
            maxError,
            domainExit
        });
        
        completed++;
        if (completed % 50 === 0) {
            self.postMessage({ type: 'progress', completed, total: seeds.length });
        }
    }
    
    self.postMessage({ type: 'complete', results });
}
