importScripts('./trilinear.js');

self.onmessage = function(e) {
  try {
    const { type, catalog, N, halfExtent, cosmology } = e.data;
    
    if (type !== 'build') return;
    
    // Stage 1: Cloud-in-Cell (CIC) density field
    postProgress('cic', 0);
    const N3 = N * N * N;
    const density = new Float32Array(N3);
    
    // Grid geometry
    const L = 2 * halfExtent;
    const dx = L / (N - 1);
    
    for (let i = 0; i < catalog.length; i++) {
      const gal = catalog[i];
      // Compute grid coordinates
      const gx = (gal.sgx + halfExtent) / L * (N - 1);
      const gy = (gal.sgy + halfExtent) / L * (N - 1);
      const gz = (gal.sgz + halfExtent) / L * (N - 1);
      
      const ix = Math.floor(gx);
      const iy = Math.floor(gy);
      const iz = Math.floor(gz);
      
      const dx1 = gx - ix;
      const dy1 = gy - iy;
      const dz1 = gz - iz;
      const dx0 = 1.0 - dx1;
      const dy0 = 1.0 - dy1;
      const dz0 = 1.0 - dz1;
      
      // Distribute mass to 8 surrounding voxels (trilinear weights)
      if (ix >= 0 && ix < N - 1 && iy >= 0 && iy < N - 1 && iz >= 0 && iz < N - 1) {
        density[ix * N * N + iy * N + iz] += dx0 * dy0 * dz0;
        density[(ix + 1) * N * N + iy * N + iz] += dx1 * dy0 * dz0;
        density[ix * N * N + (iy + 1) * N + iz] += dx0 * dy1 * dz0;
        density[ix * N * N + iy * N + (iz + 1)] += dx0 * dy0 * dz1;
        density[(ix + 1) * N * N + (iy + 1) * N + iz] += dx1 * dy1 * dz0;
        density[(ix + 1) * N * N + iy * N + (iz + 1)] += dx1 * dy0 * dz1;
        density[ix * N * N + (iy + 1) * N + (iz + 1)] += dx0 * dy1 * dz1;
        density[(ix + 1) * N * N + (iy + 1) * N + (iz + 1)] += dx1 * dy1 * dz1;
      }
      
      if (i % 1000 === 0) {
        postProgress('cic', (i / catalog.length) * 20);
      }
    }
    
    let totalMass = 0;
    for (let i = 0; i < N3; i++) {
      totalMass += density[i];
    }
    const meanDensity = totalMass / N3;
    
    const delta = new Float32Array(N3);
    for (let i = 0; i < N3; i++) {
      delta[i] = density[i] / meanDensity - 1;
    }
    postProgress('cic', 30);
    
    // Stage 2: Wiener-filter velocity reconstruction
    postProgress('wiener', 30);
    const velocityGrid = new Float32Array(N3 * 3);
    const kernelRadius = 20; // cells
    const kernelRadius2 = kernelRadius * kernelRadius;
    const dV = dx * dx * dx;
    const prefactor = (cosmology.H0 * cosmology.f) / (4 * Math.PI);
    
    const sigma = 4.0;
    const sigma2 = sigma * sigma;
    
    // Real-space convolution
    for (let x = 0; x < N; x++) {
      for (let y = 0; y < N; y++) {
        for (let z = 0; z < N; z++) {
          const idx = (x * N * N + y * N + z) * 3;
          let vx = 0, vy = 0, vz = 0;
          
          const startX = Math.max(0, x - kernelRadius);
          const endX = Math.min(N - 1, x + kernelRadius);
          const startY = Math.max(0, y - kernelRadius);
          const endY = Math.min(N - 1, y + kernelRadius);
          const startZ = Math.max(0, z - kernelRadius);
          const endZ = Math.min(N - 1, z + kernelRadius);
          
          for (let nx = startX; nx <= endX; nx++) {
            for (let ny = startY; ny <= endY; ny++) {
              for (let nz = startZ; nz <= endZ; nz++) {
                if (nx === x && ny === y && nz === z) continue;
                
                const r2_grid = (nx - x) * (nx - x) + (ny - y) * (ny - y) + (nz - z) * (nz - z);
                
                if (r2_grid <= kernelRadius2) {
                  const dxx = (nx - x) * dx;
                  const dyy = (ny - y) * dx;
                  const dzz = (nz - z) * dx;
                  
                  const r2 = r2_grid * dx * dx;
                  const r = Math.sqrt(r2);
                  const r3 = r2 * r;
                  
                  // Gaussian smoothing term (acts as Wiener filter here)
                  const smooth = Math.exp(-r2_grid / (2 * sigma2));
                  const nidx = nx * N * N + ny * N + nz;
                  
                  const term = (delta[nidx] * dV * smooth) / r3;
                  vx += prefactor * term * dxx;
                  vy += prefactor * term * dyy;
                  vz += prefactor * term * dzz;
                }
              }
            }
          }
          velocityGrid[idx] = vx;
          velocityGrid[idx + 1] = vy;
          velocityGrid[idx + 2] = vz;
        }
      }
      
      // Update progress every few planes
      if (x % 4 === 0) {
        postProgress('wiener', 30 + (x / N) * 50);
      }
    }
    
    postProgress('wiener', 80);
    
    // Stage 3: Normalization
    postProgress('normalize', 80);
    
    // Find cell closest to observer at (0,0,0)
    const obsX = Math.floor((0 + halfExtent) / L * (N - 1));
    const obsY = Math.floor((0 + halfExtent) / L * (N - 1));
    const obsZ = Math.floor((0 + halfExtent) / L * (N - 1));
    
    let obsVx = 0, obsVy = 0, obsVz = 0;
    if (obsX >= 0 && obsX < N && obsY >= 0 && obsY < N && obsZ >= 0 && obsZ < N) {
       const oIdx = (obsX * N * N + obsY * N + obsZ) * 3;
       obsVx = velocityGrid[oIdx];
       obsVy = velocityGrid[oIdx + 1];
       obsVz = velocityGrid[oIdx + 2];
    }
    const obsSpeed = Math.sqrt(obsVx * obsVx + obsVy * obsVy + obsVz * obsVz) || 1;
    
    // Target CMB dipole scale
    const targetSpeed = 627.0;
    const scale = targetSpeed / obsSpeed;
    
    let sumSpeed = 0;
    let maxSpeed = 0;
    let sumSqSpeed = 0;
    
    for (let i = 0; i < N3; i++) {
      const idx = i * 3;
      velocityGrid[idx] *= scale;
      velocityGrid[idx + 1] *= scale;
      velocityGrid[idx + 2] *= scale;
      
      const v2 = velocityGrid[idx] * velocityGrid[idx] + 
                 velocityGrid[idx + 1] * velocityGrid[idx + 1] + 
                 velocityGrid[idx + 2] * velocityGrid[idx + 2];
      const speed = Math.sqrt(v2);
      sumSpeed += speed;
      sumSqSpeed += v2;
      if (speed > maxSpeed) maxSpeed = speed;
    }
    
    const meanSpeed = sumSpeed / N3;
    const rmsSpeed = Math.sqrt(sumSqSpeed / N3);
    const bulkFlow = targetSpeed;
    
    postProgress('normalize', 100);
    
    self.postMessage({
      type: 'complete',
      grid: velocityGrid.buffer,
      stats: { meanSpeed, maxSpeed, rmsSpeed, bulkFlow }
    }, [velocityGrid.buffer]);
    
  } catch (err) {
    self.postMessage({ type: 'error', message: err.message || String(err) });
  }
};

function postProgress(stage, percent) {
  self.postMessage({ type: 'progress', stage, percent: Math.round(percent) });
}
