// spherex-engine/backend/src/services/s3RangeFetcher.ts

import https from 'https';
import {
  calibrateAndMaskCutout,
  decodeBigEndianFlags32,
  decodeBigEndianFloat32,
} from '../core/fitsParser';
import {
  computeContiguousRowSpanByteRange,
  skyToDetectorPixel,
} from '../core/wcsProjection';
import {
  HydratedNativeCutout,
  SidecarExposureEntry,
  ViewportAnalysisRequest,
} from '../types/astro.types';
import { synthesize24RowBigEndianHduSlices } from './syntheticMefGenerator';

/**
 * Persistent HTTP/1.1 Keep-Alive connection pool configured for high-concurrency
 * S3 byte-range reads (320 concurrent sockets as specified in Section 1.2).
 */
export const s3HttpKeepAliveAgent = new https.Agent({
  keepAlive: true,
  keepAliveMsecs: 30000,
  maxSockets: 320,
  maxFreeSockets: 64,
  timeout: 2500,
});

/**
 * Executes a single HTTP/1.1 GET request with a Range: bytes=start-end header.
 */
async function fetchRemoteS3ByteRange(
  url: string,
  rangeHeaderValue: string
): Promise<Buffer> {
  return new Promise((resolve, reject) => {
    const req = https.get(
      url,
      {
        agent: s3HttpKeepAliveAgent,
        headers: {
          Range: rangeHeaderValue,
          Accept: 'application/octet-stream',
        },
        timeout: 2400,
      },
      (res) => {
        if (res.statusCode !== 206 && res.statusCode !== 200) {
          reject(new Error(`S3 Range GET failed with HTTP ${res.statusCode}`));
          res.resume();
          return;
        }
        const chunks: Buffer[] = [];
        res.on('data', (chunk: Buffer) => chunks.push(chunk));
        res.on('end', () => resolve(Buffer.concat(chunks)));
      }
    );
    req.on('error', reject);
    req.on('timeout', () => {
      req.destroy();
      reject(new Error('S3 Range GET timed out (exceeded 2400 ms SLA)'));
    });
  });
}

/**
 * Slices a 24-row x 2040-column detector span down to the local 24 x 24 column window
 * centered on iMin0Based.
 */
function sliceColumnsFloat32(
  fullRowSpan: Float32Array,
  rowCount: number,
  detectorWidth: number,
  iMin0Based: number,
  colCount: number
): Float32Array {
  const out = new Float32Array(rowCount * colCount);
  for (let r = 0; r < rowCount; r++) {
    const srcOffset = r * detectorWidth + iMin0Based;
    const dstOffset = r * colCount;
    for (let c = 0; c < colCount; c++) {
      out[dstOffset + c] = fullRowSpan[srcOffset + c];
    }
  }
  return out;
}

function sliceColumnsUint32(
  fullRowSpan: Uint32Array,
  rowCount: number,
  detectorWidth: number,
  iMin0Based: number,
  colCount: number
): Uint32Array {
  const out = new Uint32Array(rowCount * colCount);
  for (let r = 0; r < rowCount; r++) {
    const srcOffset = r * detectorWidth + iMin0Based;
    const dstOffset = r * colCount;
    for (let c = 0; c < colCount; c++) {
      out[dstOffset + c] = fullRowSpan[srcOffset + c];
    }
  }
  return out;
}

/**
 * Fetches and hydrates a single exposure's 24-row span across HDU 1 (SCI),
 * HDU 2 (FLAGS), HDU 3 (VARIANCE), and HDU 4 (ZODI) (783,360 bytes total).
 */
export async function fetchAndHydrateExposureCutout(params: {
  exposure: SidecarExposureEntry;
  raDeg: number;
  decDeg: number;
  rejectBitmask: number;
  preset: NonNullable<ViewportAnalysisRequest['simulation_preset']>;
}): Promise<{
  cutout: HydratedNativeCutout;
  bytesFetched: number;
}> {
  const { exposure, raDeg, decDeg, rejectBitmask, preset } = params;
  const DETECTOR_WIDTH = 2040;
  const DETECTOR_HEIGHT = 2040;

  // Project target (ra, dec) to 0-based detector coordinates
  const proj = skyToDetectorPixel(raDeg, decDeg, exposure.wcs);

  // Compute 24-row contiguous span for each of the 4 image HDUs
  const sciRange = computeContiguousRowSpanByteRange({
    j0BasedCenter: proj.j0BasedInt,
    rollAngleRad: exposure.wcs.rollAngleRad,
    hduDataStartByte: exposure.hduDataOffsets.sciStartByte,
  });

  const flagsRange = computeContiguousRowSpanByteRange({
    j0BasedCenter: proj.j0BasedInt,
    rollAngleRad: exposure.wcs.rollAngleRad,
    hduDataStartByte: exposure.hduDataOffsets.flagsStartByte,
  });

  const varRange = computeContiguousRowSpanByteRange({
    j0BasedCenter: proj.j0BasedInt,
    rollAngleRad: exposure.wcs.rollAngleRad,
    hduDataStartByte: exposure.hduDataOffsets.varianceStartByte,
  });

  const zodiRange = computeContiguousRowSpanByteRange({
    j0BasedCenter: proj.j0BasedInt,
    rollAngleRad: exposure.wcs.rollAngleRad,
    hduDataStartByte: exposure.hduDataOffsets.zodiStartByte,
  });

  let sciBuf: Buffer;
  let flagsBuf: Buffer;
  let varBuf: Buffer;
  let zodiBuf: Buffer;

  if (process.env.SPHEREX_S3_MODE === 'LIVE') {
    [sciBuf, flagsBuf, varBuf, zodiBuf] = await Promise.all([
      fetchRemoteS3ByteRange(exposure.s3Uri, sciRange.rangeHeaderValue),
      fetchRemoteS3ByteRange(exposure.s3Uri, flagsRange.rangeHeaderValue),
      fetchRemoteS3ByteRange(exposure.s3Uri, varRange.rangeHeaderValue),
      fetchRemoteS3ByteRange(exposure.s3Uri, zodiRange.rangeHeaderValue),
    ]);
  } else {
    // Synthesize authentic Big-Endian IEEE-754 24-row buffers in memory
    const synth = synthesize24RowBigEndianHduSlices({
      exposure,
      centerRaDeg: raDeg,
      centerDecDeg: decDeg,
      jMin0Based: sciRange.jMin0Based,
      rowCount: sciRange.rowCount,
      preset,
    });
    sciBuf = synth.sciBigEndian;
    flagsBuf = synth.flagsBigEndian;
    varBuf = synth.varianceBigEndian;
    zodiBuf = synth.zodiBigEndian;
  }

  const bytesFetched =
    sciBuf.byteLength + flagsBuf.byteLength + varBuf.byteLength + zodiBuf.byteLength; // 783,360 bytes

  // Decode Big-Endian IEEE-754 buffers to Little-Endian TypedArrays via native Buffer.swap32()
  const fullSci = decodeBigEndianFloat32(sciBuf, 1.0, 0.0);
  const fullFlags = decodeBigEndianFlags32(flagsBuf, 0);
  const fullVar = decodeBigEndianFloat32(varBuf, 1.0, 0.0);
  const fullZodi = decodeBigEndianFloat32(zodiBuf, 1.0, 0.0);

  // Trim horizontally to the local 24-column window around i0BasedInt
  const colCount = sciRange.rowCount; // 24x24 local native patch
  let iMin0Based = Math.max(0, proj.i0BasedInt - sciRange.rDetRowsHalf);
  if (iMin0Based + colCount > DETECTOR_WIDTH) {
    iMin0Based = DETECTOR_WIDTH - colCount;
  }

  const localSci = sliceColumnsFloat32(
    fullSci,
    sciRange.rowCount,
    DETECTOR_WIDTH,
    iMin0Based,
    colCount
  );
  const localFlags = sliceColumnsUint32(
    fullFlags,
    sciRange.rowCount,
    DETECTOR_WIDTH,
    iMin0Based,
    colCount
  );
  const localVar = sliceColumnsFloat32(
    fullVar,
    sciRange.rowCount,
    DETECTOR_WIDTH,
    iMin0Based,
    colCount
  );
  const localZodi = sliceColumnsFloat32(
    fullZodi,
    sciRange.rowCount,
    DETECTOR_WIDTH,
    iMin0Based,
    colCount
  );

  // Apply Zodiacal subtraction, Bit-2 (0x0004) variance x2.25 inflation, and bitmask
  const calibrated = calibrateAndMaskCutout({
    rawSci: localSci,
    rawFlags: localFlags,
    rawVarOrErr: localVar,
    rawZodi: localZodi,
    varianceIsSigma: exposure.hduDataOffsets.varianceIsSigma,
    rejectBitmask,
  });

  return {
    bytesFetched,
    cutout: {
      exposure,
      rowStart0Based: sciRange.jMin0Based,
      rowCount: sciRange.rowCount,
      colStart0Based: iMin0Based,
      colCount,
      targetDetectorX0Based: proj.x0BasedContinuous,
      targetDetectorY0Based: proj.y0BasedContinuous,
      fluxMjySr: calibrated.fluxMjySr,
      varianceMjySr2: calibrated.varianceMjySr2,
      flags: localFlags,
      validMask: calibrated.validMask,
      saturatedFraction: calibrated.saturatedFraction,
      hasPersistenceFlag: calibrated.hasPersistenceFlag,
      hasGhostFlag: calibrated.hasGhostFlag,
    },
  };
}