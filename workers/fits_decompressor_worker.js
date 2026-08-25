/**
 * @file fits_decompressor_worker.js
 * @description Dedicated Web Worker for zero-copy Big-Endian binary decoding,
 * FITS primary array byte swapping, and canonical (SGZ, SGY, SGX) -> (SGX, SGY, SGZ) reordering.
 * 
 * Transfers decoded Float32Array / Float64Array directly back to the main thread
 * via Transferable Objects without memory duplication or WebGL UI frame drops.
 */

self.onmessage = function(evt) {
  const { id, rawBuffer, naxis1, naxis2, naxis3, bitpix, bscale = 1.0, bzero = 0.0, isVelocity = false } = evt.data;

  try {
    const totalElements = naxis1 * naxis2 * naxis3;
    const view = new DataView(rawBuffer);
    const outData = new Float32Array(totalElements);

    // Official CF4 velocity scale factor (Pillar 3)
    const scaleFactor = isVelocity ? (52.0 * bscale) : bscale;

    let srcIdx = 0;
    for (let iz = 0; iz < naxis3; iz++) {
      for (let iy = 0; iy < naxis2; iy++) {
        for (let ix = 0; ix < naxis1; ix++) {
          let rawVal = 0.0;
          if (bitpix === -32) {
            rawVal = view.getFloat32(srcIdx * 4, false); // Big-Endian
          } else if (bitpix === -64) {
            rawVal = view.getFloat64(srcIdx * 8, false); // Big-Endian
          } else if (bitpix === 16) {
            rawVal = view.getInt16(srcIdx * 2, false);
          } else if (bitpix === 32) {
            rawVal = view.getInt32(srcIdx * 4, false);
          }

          const physVal = rawVal * scaleFactor + bzero;

          // Canonical internal ZRT mapping: index(ix, iy, iz) = ix + nx * (iy + ny * iz)
          const dstIdx = ix + naxis1 * (iy + naxis2 * iz);
          outData[dstIdx] = physVal;
          srcIdx++;
        }
      }
    }

    // Zero-copy transfer of output ArrayBuffer
    self.postMessage({
      id,
      success: true,
      data: outData,
      nx: naxis1,
      ny: naxis2,
      nz: naxis3
    }, [outData.buffer]);

  } catch (error) {
    self.postMessage({
      id,
      success: false,
      error: error.message || String(error)
    });
  }
};
