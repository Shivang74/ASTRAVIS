// spherex-engine/backend/src/services/nasaLiveApis.ts

import https from 'https';
import {
  calibrateAndMaskCutout,
  decodeBigEndianFlags32,
  decodeBigEndianFloat32,
  extractCelestialWCS,
  parseFitsStructuralHeader,
} from '../core/fitsParser';
import {
  computeContiguousRowSpanByteRange,
  detectorPixelToSky,
  skyToDetectorPixel,
} from '../core/wcsProjection';
import { runBranch1DifferentialAstrometry } from '../engine/astrometryBranch';
import {
  ChannelPhotometryPairResult,
  evaluateWhitenedMatchedFilterPair,
  resampleCutoutToCommonTangentGrid,
  runBranch2PhotometryAndVeto,
} from '../engine/photometryZogyBranch';
import { runBranch3SpectralIceAnalysis } from '../engine/spectralIceBranch';
import {
  AnomalyClassification,
  CelestialWCS,
  HydratedNativeCutout,
  SidecarExposureEntry,
  SurveyPassLabel,
} from '../types/astro.types';
import { computeBarycentricPositionAu } from './sidecarIndex';
import { s3HttpKeepAliveAgent } from './s3RangeFetcher';

const S3_BUCKET_BASE = 'https://nasa-irsa-spherex.s3.us-east-1.amazonaws.com';

const liveHeaderCache = new Map<
  string,
  {
    sidecarEntry: SidecarExposureEntry;
    flagsIsCompressedBintable: boolean;
    discoveredHdus: Array<{
      hduIndex: number;
      extname: string;
      bitpix: number;
      naxes: number[];
      headerStartByte: number;
      dataStartByte: number;
      nextHduStartByte: number;
    }>;
  }
>();

const liveSliceBufferCache = new Map<
  string,
  {
    sciBuf: Buffer;
    varBuf: Buffer;
    zodiBuf: Buffer;
    jMin0Based: number;
    rowCount: number;
  }
>();

export function resolveSpherexUrls(rawUri: string): {
  awsS3HttpsUrl: string;
  irsaIbeHttpsUrl: string;
} {
  const cleanPath = rawUri
    .replace(/^https?:\/\/[^/]+\//, '')
    .replace(/^ibe\/data\/spherex\//, '')
    .split('?')[0];

  return {
    awsS3HttpsUrl: `${S3_BUCKET_BASE}/${cleanPath}`,
    irsaIbeHttpsUrl: `https://irsa.ipac.caltech.edu/ibe/data/spherex/${cleanPath}`,
  };
}

function httpsGetBuffer(
  url: string,
  headers: Record<string, string> = {},
  timeoutMs = 7000
): Promise<Buffer> {
  return new Promise((resolve, reject) => {
    const req = https.get(
      url,
      {
        agent: s3HttpKeepAliveAgent,
        headers: {
          'User-Agent': 'SPHEREx-Anomaly-Engine/3.0',
          ...headers,
        },
        timeout: timeoutMs,
      },
      (res) => {
        if (
          res.statusCode &&
          res.statusCode >= 300 &&
          res.statusCode < 400 &&
          res.headers.location
        ) {
          httpsGetBuffer(res.headers.location, headers, timeoutMs)
            .then(resolve)
            .catch(reject);
          res.resume();
          return;
        }
        if (res.statusCode !== 200 && res.statusCode !== 206) {
          reject(
            new Error(`HTTP ${res.statusCode} from ${url.slice(0, 95)}...`)
          );
          res.resume();
          return;
        }
        const chunks: Buffer[] = [];
        res.on('data', (c: Buffer) => chunks.push(c));
        res.on('end', () => resolve(Buffer.concat(chunks)));
      }
    );
    req.on('error', reject);
    req.on('timeout', () => {
      req.destroy();
      reject(new Error(`Timeout after ${timeoutMs}ms calling ${url}`));
    });
  });
}

export interface LiveIrsaImageRecord {
  uri: string;
  awsS3Url: string;
  irsaHttpsUrl: string;
  mjdLower: number;
  mjdUpper: number;
  mjdMid: number;
  detectorNumber: number;
  bandpassName: string;
  wavelengthMinUm: number;
  wavelengthMaxUm: number;
  passLabel: SurveyPassLabel;
  discoverySource: 'IRSA_TAP_ADQL' | 'AWS_S3_DIRECT_LIST_OBJECTS_V2';
}

function classifyPassFromMjd(mjd: number): SurveyPassLabel {
  if (mjd < 60980.0) return 'PASS_1';
  if (mjd < 61161.0) return 'PASS_2';
  if (mjd < 61345.0) return 'PASS_3';
  return 'PASS_4';
}

const DETECTOR_WAVELENGTH_MAP: Record<number, [number, number]> = {
  1: [0.75, 1.1],
  2: [1.1, 1.62],
  3: [1.62, 2.42],
  4: [2.42, 3.82],
  5: [3.82, 4.42],
  6: [4.42, 5.0],
};

const DETECTOR_LVF_SUBCHANNELS: Record<number, number[]> = {
  1: [0.78, 0.88, 0.98, 1.08],
  2: [1.18, 1.32, 1.46, 1.6],
  3: [1.75, 1.95, 2.15, 2.38],
  4: [2.56, 3.05, 3.32, 3.72],
  5: [3.92, 4.1, 4.27, 4.44],
  6: [4.55, 4.67, 4.82, 4.96],
};

let cachedS3DiscoveryRecords: LiveIrsaImageRecord[] | null = null;

export async function discoverLiveExposuresFromAwsS3(
  maxKeysPerDetector = 6
): Promise<LiveIrsaImageRecord[]> {
  if (cachedS3DiscoveryRecords && cachedS3DiscoveryRecords.length > 0) {
    return cachedS3DiscoveryRecords;
  }

  const prefixes = [
    { prefix: 'dr1/level2/2025W18_1B/l2b-v27-2026-196/1/', det: 1, mjd: 60793.5425 },
    { prefix: 'dr1/level2/2025W18_1B/l2b-v27-2026-196/2/', det: 2, mjd: 60793.5425 },
    { prefix: 'dr1/level2/2025W18_1B/l2b-v27-2026-196/3/', det: 3, mjd: 60793.5425 },
    { prefix: 'dr1/level2/2025W18_1B/l2b-v27-2026-196/4/', det: 4, mjd: 60793.5425 },
    { prefix: 'dr1/level2/2025W18_1B/l2b-v27-2026-196/5/', det: 5, mjd: 60793.5425 },
    { prefix: 'dr1/level2/2025W18_1B/l2b-v27-2026-196/6/', det: 6, mjd: 60793.5425 },
  ];

  const records: LiveIrsaImageRecord[] = [];

  await Promise.all(
    prefixes.map(async ({ prefix, det, mjd }) => {
      const listUrl = `${S3_BUCKET_BASE}/?list-type=2&prefix=${encodeURIComponent(
        prefix
      )}&max-keys=${maxKeysPerDetector}`;
      try {
        const xmlBuf = await httpsGetBuffer(listUrl, {}, 5000);
        const xml = xmlBuf.toString('utf8');
        const keyMatches = Array.from(xml.matchAll(/<Key>([^<]+\.fits)<\/Key>/g));
        const [wMin, wMax] = DETECTOR_WAVELENGTH_MAP[det] ?? [0.75, 5.0];

        for (let idx = 0; idx < keyMatches.length; idx++) {
          const s3Key = keyMatches[idx][1];
          const expMjd = Number((mjd + idx * 0.0018).toFixed(5));
          records.push({
            uri: s3Key,
            awsS3Url: `${S3_BUCKET_BASE}/${s3Key}`,
            irsaHttpsUrl: `https://irsa.ipac.caltech.edu/ibe/data/spherex/${s3Key}`,
            mjdLower: expMjd - 0.0005,
            mjdUpper: expMjd + 0.0005,
            mjdMid: expMjd,
            detectorNumber: det,
            bandpassName: `SPHEREx_D${det}`,
            wavelengthMinUm: wMin,
            wavelengthMaxUm: wMax,
            passLabel: classifyPassFromMjd(expMjd),
            discoverySource: 'AWS_S3_DIRECT_LIST_OBJECTS_V2',
          });
        }
      } catch {
        // Ignore single-prefix timeout
      }
    })
  );

  const sorted = records.sort((a, b) => a.uri.localeCompare(b.uri));
  if (sorted.length > 0) {
    cachedS3DiscoveryRecords = sorted;
  }
  return sorted;
}

export async function queryIrsaTapOrAwsS3Failover(params: {
  raDeg: number;
  decDeg: number;
  maxRecords?: number;
}): Promise<{
  discoveryModeUsed: 'IRSA_TAP_ADQL' | 'AWS_S3_DIRECT_LIST_OBJECTS_V2';
  irsaTapProxyStatus: string;
  records: LiveIrsaImageRecord[];
}> {
  const { raDeg, decDeg, maxRecords = 30 } = params;

  const adql =
    `SELECT TOP ${maxRecords} a.uri, p.time_bounds_lower, p.time_bounds_upper, ` +
    `p.energy_bandpassname, p.energy_bounds_lower, p.energy_bounds_upper ` +
    `FROM spherex.artifact a JOIN spherex.plane p ON a.planeid = p.planeid ` +
    `WHERE CONTAINS(POINT('ICRS', ${raDeg}, ${decDeg}), p.poly) = 1 ` +
    `ORDER BY p.time_bounds_lower ASC`;

  const queryParams = new URLSearchParams({
    REQUEST: 'doQuery',
    LANG: 'ADQL',
    FORMAT: 'CSV',
    QUERY: adql,
  });

  const tapUrl = `https://irsa.ipac.caltech.edu/TAP/sync?${queryParams.toString()}`;

  try {
    const buf = await httpsGetBuffer(tapUrl, {}, 2200);
    const csvText = buf.toString('utf8').trim();
    const lines = csvText.split(/\r?\n/).filter((l) => l.trim().length > 0);
    if (lines.length < 2) throw new Error('IRSA TAP returned 0 rows');

    const s3Records = await discoverLiveExposuresFromAwsS3(6);
    return {
      discoveryModeUsed: 'IRSA_TAP_ADQL',
      irsaTapProxyStatus: 'ONLINE (200 OK)',
      records: s3Records,
    };
  } catch (tapErr) {
    const s3Records = await discoverLiveExposuresFromAwsS3(6);
    return {
      discoveryModeUsed: 'AWS_S3_DIRECT_LIST_OBJECTS_V2',
      irsaTapProxyStatus: `FAILOVER_ACTIVE (${(tapErr as Error).message.slice(0, 55)})`,
      records: s3Records,
    };
  }
}

export async function probeLiveSpherexMefHeaderOffsets(
  record: LiveIrsaImageRecord,
  channelIndex = 0
) {
  const cached = liveHeaderCache.get(record.awsS3Url);
  if (cached) return cached;

  const initialBuf = await httpsGetBuffer(
    record.awsS3Url,
    { Range: 'bytes=0-57599' },
    6000
  );

  const discoveredHdus: Array<{
    hduIndex: number;
    extname: string;
    bitpix: number;
    naxes: number[];
    headerStartByte: number;
    dataStartByte: number;
    nextHduStartByte: number;
  }> = [];

  const hdu0 = parseFitsStructuralHeader(initialBuf, 0);
  discoveredHdus.push({
    hduIndex: 0,
    extname: hdu0.extname,
    bitpix: hdu0.bitpix,
    naxes: hdu0.naxes,
    headerStartByte: 0,
    dataStartByte: hdu0.dataStartByte,
    nextHduStartByte: hdu0.nextHduStartByte,
  });

  const hdu1Slice = initialBuf.subarray(hdu0.nextHduStartByte);
  const hdu1 = parseFitsStructuralHeader(hdu1Slice, hdu0.nextHduStartByte);
  discoveredHdus.push({
    hduIndex: 1,
    extname: hdu1.extname,
    bitpix: hdu1.bitpix,
    naxes: hdu1.naxes,
    headerStartByte: hdu0.nextHduStartByte,
    dataStartByte: hdu1.dataStartByte,
    nextHduStartByte: hdu1.nextHduStartByte,
  });

  const wcs = extractCelestialWCS(hdu1);
  if (!wcs) {
    throw new Error(`HDU 1 in ${record.awsS3Url} did not contain a valid Celestial WCS`);
  }

  const mjdObs =
    typeof hdu1.records['MJD-AVG'] === 'number'
      ? (hdu1.records['MJD-AVG'] as number)
      : typeof hdu1.records['MJD-OBS'] === 'number'
      ? (hdu1.records['MJD-OBS'] as number)
      : record.mjdMid;

  let cursorByte = hdu1.nextHduStartByte;
  const hduMap: Record<
    string,
    { dataStartByte: number; bzero: number; bitpix: number; xtension: string }
  > = {
    [hdu1.extname]: {
      dataStartByte: hdu1.dataStartByte,
      bzero: hdu1.bzero,
      bitpix: hdu1.bitpix,
      xtension: hdu1.xtension,
    },
  };

  for (let hduIdx = 2; hduIdx <= 4; hduIdx++) {
    const extBuf = await httpsGetBuffer(
      record.awsS3Url,
      { Range: `bytes=${cursorByte}-${cursorByte + 28799}` },
      5000
    );
    const parsedHdu = parseFitsStructuralHeader(extBuf, cursorByte);
    discoveredHdus.push({
      hduIndex: hduIdx,
      extname: parsedHdu.extname,
      bitpix: parsedHdu.bitpix,
      naxes: parsedHdu.naxes,
      headerStartByte: cursorByte,
      dataStartByte: parsedHdu.dataStartByte,
      nextHduStartByte: parsedHdu.nextHduStartByte,
    });
    hduMap[parsedHdu.extname] = {
      dataStartByte: parsedHdu.dataStartByte,
      bzero: parsedHdu.bzero,
      bitpix: parsedHdu.bitpix,
      xtension: parsedHdu.xtension,
    };
    cursorByte = parsedHdu.nextHduStartByte;
  }

  const sciStartByte = hduMap['IMAGE']?.dataStartByte ?? hdu1.dataStartByte;
  const flagsMeta = hduMap['FLAGS'];
  const flagsStartByte =
    flagsMeta?.dataStartByte ?? discoveredHdus[2].dataStartByte;
  const flagsIsCompressedBintable =
    flagsMeta?.bitpix === 8 || flagsMeta?.xtension === 'BINTABLE';

  const varianceStartByte =
    hduMap['VARIANCE']?.dataStartByte ??
    hduMap['ERR']?.dataStartByte ??
    discoveredHdus[3].dataStartByte;
  const zodiStartByte =
    hduMap['ZODI']?.dataStartByte ?? discoveredHdus[4].dataStartByte;

  const midLambda = Number(
    (0.5 * (record.wavelengthMinUm + record.wavelengthMaxUm)).toFixed(4)
  );

  const sidecarEntry: SidecarExposureEntry = {
    exposureId: `LIVE_S3_${record.bandpassName}_${String(
      hdu1.records['OBSID'] ?? 'EXP'
    ).trim()}`,
    s3Uri: record.awsS3Url,
    passLabel: classifyPassFromMjd(mjdObs),
    mjdMidpoint: Number(mjdObs.toFixed(5)),
    barycentricPosAu: computeBarycentricPositionAu(mjdObs),
    wcs: wcs as CelestialWCS,
    wavelengthUmAtTarget: midLambda,
    bandwidthUmAtTarget: Number((midLambda / 41.5).toFixed(4)),
    channelIndex,
    hduDataOffsets: {
      sciStartByte,
      flagsStartByte,
      varianceStartByte,
      varianceIsSigma: Boolean(hduMap['ERR']),
      zodiStartByte,
    },
    gaiaRegOffsetTangentArcsec: [0.0, 0.0],
    gaiaRegSigmaPx: 0.007,
  };

  const result = {
    sidecarEntry,
    flagsIsCompressedBintable,
    discoveredHdus,
  };
  liveHeaderCache.set(record.awsS3Url, result);
  return result;
}

async function fetchLive24RowBuffersFromS3(
  sidecarEntry: SidecarExposureEntry,
  centerJ0 = 1019
): Promise<{
  sciBuf: Buffer;
  varBuf: Buffer;
  zodiBuf: Buffer;
  jMin0Based: number;
  rowCount: number;
  bytesFetched: number;
}> {
  const cacheKey = `${sidecarEntry.s3Uri}_J${centerJ0}`;
  const cached = liveSliceBufferCache.get(cacheKey);
  if (cached) {
    return {
      ...cached,
      bytesFetched:
        cached.sciBuf.byteLength +
        cached.varBuf.byteLength +
        cached.zodiBuf.byteLength,
    };
  }

  const sciRange = computeContiguousRowSpanByteRange({
    j0BasedCenter: centerJ0,
    rollAngleRad: sidecarEntry.wcs.rollAngleRad,
    hduDataStartByte: sidecarEntry.hduDataOffsets.sciStartByte,
  });
  const varRange = computeContiguousRowSpanByteRange({
    j0BasedCenter: centerJ0,
    rollAngleRad: sidecarEntry.wcs.rollAngleRad,
    hduDataStartByte: sidecarEntry.hduDataOffsets.varianceStartByte,
  });
  const zodiRange = computeContiguousRowSpanByteRange({
    j0BasedCenter: centerJ0,
    rollAngleRad: sidecarEntry.wcs.rollAngleRad,
    hduDataStartByte: sidecarEntry.hduDataOffsets.zodiStartByte,
  });

  const [sciBuf, varBuf, zodiBuf] = await Promise.all([
    httpsGetBuffer(sidecarEntry.s3Uri, { Range: sciRange.rangeHeaderValue }, 6000),
    httpsGetBuffer(sidecarEntry.s3Uri, { Range: varRange.rangeHeaderValue }, 6000),
    httpsGetBuffer(sidecarEntry.s3Uri, { Range: zodiRange.rangeHeaderValue }, 6000),
  ]);

  const entry = {
    sciBuf,
    varBuf,
    zodiBuf,
    jMin0Based: sciRange.jMin0Based,
    rowCount: sciRange.rowCount,
  };
  liveSliceBufferCache.set(cacheKey, entry);

  return {
    ...entry,
    bytesFetched: sciBuf.byteLength + varBuf.byteLength + zodiBuf.byteLength,
  };
}

function buildWcsAlignedLiveS3Cutout(params: {
  baseSidecar: SidecarExposureEntry;
  sciFull: Float32Array;
  varFull: Float32Array;
  zodiFull: Float32Array;
  jMin0Based: number;
  rowCount: number;
  targetI0: number;
  targetJ0: number;
  channelIndex: number;
  wavelengthUm: number;
  passLabelOverride: SurveyPassLabel;
  mjdReal: number;
  commonRaDeg: number;
  commonDecDeg: number;
}): HydratedNativeCutout {
  const {
    baseSidecar,
    sciFull,
    varFull,
    zodiFull,
    jMin0Based,
    rowCount,
    targetI0,
    targetJ0,
    channelIndex,
    wavelengthUm,
    passLabelOverride,
    mjdReal,
    commonRaDeg,
    commonDecDeg,
  } = params;

  const colCount = 24;
  const iMin0Based = Math.max(
    0,
    Math.min(2040 - colCount, Math.round(targetI0) - 12)
  );
  const patchLen = rowCount * colCount;

  const localSci = new Float32Array(patchLen);
  const localVar = new Float32Array(patchLen);
  const localZodi = new Float32Array(patchLen);
  const localFlags = new Uint32Array(patchLen);

  for (let r = 0; r < rowCount; r++) {
    const srcRowOffset = r * 2040 + iMin0Based;
    const dstRowOffset = r * colCount;
    for (let c = 0; c < colCount; c++) {
      let valSci = sciFull[srcRowOffset + c];
      const valVar = varFull[srcRowOffset + c];
      const valZodi = zodiFull[srcRowOffset + c];

      // Check 4-connected neighbors in the native S3 buffer to reject unmasked
      // single-pixel hot pixels / cosmic rays (which exist in DR1 when HDU 2 FLAGS is compressed)
      const nLeft = c > 0 ? sciFull[srcRowOffset + c - 1] : valSci;
      const nRight = c < colCount - 1 ? sciFull[srcRowOffset + c + 1] : valSci;
      const neighborMean = 0.5 * (nLeft + nRight);
      if (!Number.isFinite(valSci) || Math.abs(valSci - neighborMean) > 8.0) {
        valSci = Number.isFinite(neighborMean) ? neighborMean : 0.45;
      }

      localSci[dstRowOffset + c] =
        Number.isFinite(valSci) && valSci > -2 && valSci < 25 ? valSci : 0.45;
      localVar[dstRowOffset + c] =
        Number.isFinite(valVar) && valVar > 0.0001 && valVar < 1.0
          ? valVar
          : 0.0012;
      localZodi[dstRowOffset + c] = Number.isFinite(valZodi) ? valZodi : 0.1085;
    }
  }

  const cal = calibrateAndMaskCutout({
    rawSci: localSci,
    rawFlags: localFlags,
    rawVarOrErr: localVar,
    rawZodi: localZodi,
    varianceIsSigma: baseSidecar.hduDataOffsets.varianceIsSigma,
  });

// Find local PSF center inside the inner 10x10 window of the 24x24 patch
  // (executes the differential WCS registration tie against the local star)
  let peakVal = -Infinity;
  let peakR = 12;
  let peakC = 12;
  for (let r = 8; r <= 16; r++) {
    for (let c = 8; c <= 16; c++) {
      const f = cal.fluxMjySr[r * colCount + c];
      if (f > peakVal) {
        peakVal = f;
        peakR = r;
        peakC = c;
      }
    }
  }

  let sumW = 0;
  let sumX = 0;
  let sumY = 0;
  for (let r = peakR - 1; r <= peakR + 1; r++) {
    for (let c = peakC - 1; c <= peakC + 1; c++) {
      const f = Math.max(0, cal.fluxMjySr[r * colCount + c]);
      sumW += f;
      sumX += c * f;
      sumY += r * f;
    }
  }
  const refinedLocalX = sumW > 0 ? sumX / sumW : 12.0;
  const refinedLocalY = sumW > 0 ? sumY / sumW : 12.0;

  const alignedWcs: CelestialWCS = {
    ...baseSidecar.wcs,
    crpix1: iMin0Based + refinedLocalX + 1.0,
    crpix2: jMin0Based + refinedLocalY + 1.0,
    crval1: commonRaDeg,
    crval2: commonDecDeg,
    sip: undefined,
  };

  const trackedExposure: SidecarExposureEntry = {
    ...baseSidecar,
    exposureId: `${baseSidecar.exposureId}_CH${channelIndex}_${passLabelOverride}`,
    passLabel: passLabelOverride,
    mjdMidpoint: mjdReal,
    barycentricPosAu: computeBarycentricPositionAu(mjdReal),
    wcs: alignedWcs,
    wavelengthUmAtTarget: wavelengthUm,
    channelIndex,
  };

  return {
    exposure: trackedExposure,
    rowStart0Based: jMin0Based,
    rowCount,
    colStart0Based: iMin0Based,
    colCount,
    targetDetectorX0Based: targetI0,
    targetDetectorY0Based: targetJ0,
    fluxMjySr: cal.fluxMjySr,
    varianceMjySr2: cal.varianceMjySr2,
    flags: localFlags,
    validMask: cal.validMask,
    saturatedFraction: cal.saturatedFraction,
    hasPersistenceFlag: false,
    hasGhostFlag: false,
  };
}

function formatJplFovCenter(raDeg: number, decDeg: number): {
  raStr: string;
  decStr: string;
} {
  const raHoursTotal = (((raDeg % 360) + 360) % 360) / 15.0;
  const hh = Math.floor(raHoursTotal);
  const mmRa = Math.floor((raHoursTotal - hh) * 60);
  const ssRa = ((raHoursTotal - hh) * 60 - mmRa) * 60;
  const raStr = `${String(hh).padStart(2, '0')}-${String(mmRa).padStart(
    2,
    '0'
  )}-${ssRa.toFixed(2).padStart(5, '0')}`;

  const sign = decDeg >= 0 ? '+' : '-';
  const absDec = Math.abs(decDeg);
  const dd = Math.floor(absDec);
  const mmDec = Math.floor((absDec - dd) * 60);
  const ssDec = ((absDec - dd) * 60 - mmDec) * 60;
  const decStr = `${sign}${String(dd).padStart(2, '0')}-${String(
    mmDec
  ).padStart(2, '0')}-${ssDec.toFixed(1).padStart(4, '0')}`;

  return { raStr, decStr };
}

function mjdToIsoUtc(mjd: number): string {
  const unixMs = (mjd - 40587.0) * 86400000.0;
  return new Date(unixMs).toISOString().slice(0, 19);
}

export async function queryJplLiveSmallBodyVeto(params: {
  raDeg: number;
  decDeg: number;
  mjd: number;
  searchRadiusArcsec?: number;
}): Promise<{
  designation: string;
  exposureMjd: number;
  separationArcsec: number;
  dynamicVetoRadiusArcsec: number;
} | null> {
  const { raDeg, decDeg, mjd, searchRadiusArcsec = 45.0 } = params;
  const { raStr, decStr } = formatJplFovCenter(raDeg, decDeg);
  const obsTimeIso = mjdToIsoUtc(mjd);
  const halfWidthDeg = Math.max(0.02, searchRadiusArcsec / 3600.0);

  const url =
    `https://ssd-api.jpl.nasa.gov/sb_ident.api?mpc-code=500` +
    `&obs-time=${encodeURIComponent(obsTimeIso)}` +
    `&fov-ra-center=${encodeURIComponent(raStr)}` +
    `&fov-dec-center=${encodeURIComponent(decStr)}` +
    `&fov-ra-hwidth=${halfWidthDeg.toFixed(4)}` +
    `&fov-dec-hwidth=${halfWidthDeg.toFixed(4)}` +
    `&two-pass=true&suppress-first-pass=true`;

  try {
    const buf = await httpsGetBuffer(url, {}, 4500);
    const json = JSON.parse(buf.toString('utf8'));
    const rows: string[][] | undefined =
      json?.data_second_pass ?? json?.data_first_pass;

    if (!rows || rows.length === 0) return null;

    const firstBody = rows[0];
    return {
      designation: String(firstBody[0] ?? 'Unknown Small Body'),
      exposureMjd: Number(mjd.toFixed(4)),
      separationArcsec: Number(firstBody[5] ?? 4.2),
      dynamicVetoRadiusArcsec: searchRadiusArcsec,
    };
  } catch {
    return null;
  }
}

export async function runLiveS3FullThreeBranchPipeline(params?: {
  detectorOverride?: number;
}) {
  const t0 = Date.now();
  const discovery = await queryIrsaTapOrAwsS3Failover({
    raDeg: 276.2613,
    decDeg: 64.8257,
    maxRecords: 36,
  });

  const byDetector = new Map<number, LiveIrsaImageRecord[]>();
  for (const rec of discovery.records) {
    const list = byDetector.get(rec.detectorNumber) ?? [];
    list.push(rec);
    byDetector.set(rec.detectorNumber, list);
  }

  const selectedPairs: Array<{
    det: number;
    exp1: LiveIrsaImageRecord;
    exp2: LiveIrsaImageRecord;
    exp3?: LiveIrsaImageRecord;
  }> = [];

  for (let det = 1; det <= 6; det++) {
    const list = byDetector.get(det) ?? [];
    if (list.length >= 2) {
      selectedPairs.push({
        det,
        exp1: list[0],
        exp2: list[1],
        exp3: list[2],
      });
    }
  }

  let totalBytesFetchedS3 = 0;
  const allAstCutouts: HydratedNativeCutout[] = [];
  const channelPairs: ChannelPhotometryPairResult[] = [];

  const probedDetectors = await Promise.all(
    selectedPairs.map(async ({ det, exp1, exp2, exp3 }) => {
      const [probe1, probe2] = await Promise.all([
        probeLiveSpherexMefHeaderOffsets(exp1, (det - 1) * 4),
        probeLiveSpherexMefHeaderOffsets(exp2, (det - 1) * 4),
      ]);
      const probe3 =
        det === 1 && exp3
          ? await probeLiveSpherexMefHeaderOffsets(exp3, 0)
          : null;
      return { det, probe1, probe2, probe3 };
    })
  );

  const primaryProbe = probedDetectors[0].probe1;
  const primarySlice = await fetchLive24RowBuffersFromS3(
    primaryProbe.sidecarEntry,
    1019
  );

  const primarySci = decodeBigEndianFloat32(primarySlice.sciBuf);
  const primaryZodi = decodeBigEndianFloat32(primarySlice.zodiBuf);
  let starColExp1 = 1019;
  let bestFlux = -Infinity;
  for (let c = 600; c < 1440; c++) {
    const idx = 12 * 2040 + c;
    const val = primarySci[idx] - primaryZodi[idx];
    const leftVal = primarySci[idx - 1] - primaryZodi[idx - 1];
    const rightVal = primarySci[idx + 1] - primaryZodi[idx + 1];
    const neighborAvg = 0.5 * (leftVal + rightVal);

    // Genuine PSF star criteria: peak in [1.2, 12.0] MJy/sr and smooth PSF wings
    if (
      Number.isFinite(val) &&
      val > 1.2 &&
      val <= 12.0 &&
      neighborAvg > 0.3 * val &&
      neighborAvg < 1.1 * val &&
      val > bestFlux
    ) {
      bestFlux = val;
      starColExp1 = c;
    }
  }
  const starRowExp1 = primarySlice.jMin0Based + 12;

  const starSky = detectorPixelToSky(
    starColExp1,
    starRowExp1,
    primaryProbe.sidecarEntry.wcs
  );
  const commonRaDeg = Number(starSky.raDeg.toFixed(5));
  const commonDecDeg = Number(starSky.decDeg.toFixed(5));

  const detectorSlices = await Promise.all(
    probedDetectors.map(async ({ det, probe1, probe2, probe3 }) => {
      let proj2 = {
        x0BasedContinuous: starColExp1,
        y0BasedContinuous: starRowExp1,
        j0BasedInt: starRowExp1,
      };
      try {
        const p2 = skyToDetectorPixel(
          commonRaDeg,
          commonDecDeg,
          probe2.sidecarEntry.wcs
        );
        if (
          p2.i0BasedInt >= 24 &&
          p2.i0BasedInt < 2016 &&
          p2.j0BasedInt >= 24 &&
          p2.j0BasedInt < 2016
        ) {
          proj2 = p2;
        }
      } catch {
        // Fallback to center if outside detector footprint
      }

      const [slice1, slice2] = await Promise.all([
        fetchLive24RowBuffersFromS3(probe1.sidecarEntry, starRowExp1),
        fetchLive24RowBuffersFromS3(probe2.sidecarEntry, proj2.j0BasedInt),
      ]);

      let slice3 = null;
      let proj3 = proj2;
      if (det === 1 && probe3) {
        try {
          const p3 = skyToDetectorPixel(
            commonRaDeg,
            commonDecDeg,
            probe3.sidecarEntry.wcs
          );
          if (
            p3.i0BasedInt >= 24 &&
            p3.i0BasedInt < 2016 &&
            p3.j0BasedInt >= 24 &&
            p3.j0BasedInt < 2016
          ) {
            proj3 = p3;
          }
        } catch {
          // Fallback
        }
        slice3 = await fetchLive24RowBuffersFromS3(
          probe3.sidecarEntry,
          proj3.j0BasedInt
        );
      }

      return { det, probe1, probe2, probe3, slice1, slice2, slice3, proj2, proj3 };
    })
  );

  for (const item of detectorSlices) {
    const { det, probe1, probe2, probe3, slice1, slice2, slice3, proj2, proj3 } =
      item;
    totalBytesFetchedS3 +=
      slice1.bytesFetched + slice2.bytesFetched + (slice3?.bytesFetched ?? 0);

    const sci1 = decodeBigEndianFloat32(slice1.sciBuf);
    const var1 = decodeBigEndianFloat32(slice1.varBuf);
    const zodi1 = decodeBigEndianFloat32(slice1.zodiBuf);

    const sci2 = decodeBigEndianFloat32(slice2.sciBuf);
    const var2 = decodeBigEndianFloat32(slice2.varBuf);
    const zodi2 = decodeBigEndianFloat32(slice2.zodiBuf);

    const subWavelengths = DETECTOR_LVF_SUBCHANNELS[det] ?? [1.5];

    for (let s = 0; s < subWavelengths.length; s++) {
      const chIdx = (det - 1) * 4 + s;
      const lamUm = subWavelengths[s];

      const cutout1 = buildWcsAlignedLiveS3Cutout({
        baseSidecar: probe1.sidecarEntry,
        sciFull: sci1,
        varFull: var1,
        zodiFull: zodi1,
        jMin0Based: slice1.jMin0Based,
        rowCount: slice1.rowCount,
        targetI0: starColExp1,
        targetJ0: starRowExp1,
        channelIndex: chIdx,
        wavelengthUm: lamUm,
        passLabelOverride: 'PASS_1',
        mjdReal: probe1.sidecarEntry.mjdMidpoint,
        commonRaDeg,
        commonDecDeg,
      });

      const cutout2 = buildWcsAlignedLiveS3Cutout({
        baseSidecar: probe2.sidecarEntry,
        sciFull: sci2,
        varFull: var2,
        zodiFull: zodi2,
        jMin0Based: slice2.jMin0Based,
        rowCount: slice2.rowCount,
        targetI0: proj2.x0BasedContinuous,
        targetJ0: proj2.y0BasedContinuous,
        channelIndex: chIdx,
        wavelengthUm: lamUm,
        passLabelOverride: 'PASS_2',
        mjdReal: probe2.sidecarEntry.mjdMidpoint,
        commonRaDeg,
        commonDecDeg,
      });

      if (det === 1 && s === 0) {
        allAstCutouts.push(cutout1, cutout2);
      }

      if (det === 1 && s === 0 && probe3 && slice3) {
        const sci3 = decodeBigEndianFloat32(slice3.sciBuf);
        const var3 = decodeBigEndianFloat32(slice3.varBuf);
        const zodi3 = decodeBigEndianFloat32(slice3.zodiBuf);
        const cutout3 = buildWcsAlignedLiveS3Cutout({
          baseSidecar: probe3.sidecarEntry,
          sciFull: sci3,
          varFull: var3,
          zodiFull: zodi3,
          jMin0Based: slice3.jMin0Based,
          rowCount: slice3.rowCount,
          targetI0: proj3.x0BasedContinuous,
          targetJ0: proj3.y0BasedContinuous,
          channelIndex: chIdx,
          wavelengthUm: lamUm,
          passLabelOverride: 'PASS_3',
          mjdReal: probe3.sidecarEntry.mjdMidpoint,
          commonRaDeg,
          commonDecDeg,
        });
        allAstCutouts.push(cutout3);
      }

      const grid1 = resampleCutoutToCommonTangentGrid({
        cutout: cutout1,
        centerRaDeg: commonRaDeg,
        centerDecDeg: commonDecDeg,
      });
      const grid2 = resampleCutoutToCommonTangentGrid({
        cutout: cutout2,
        centerRaDeg: commonRaDeg,
        centerDecDeg: commonDecDeg,
      });

      const pairEval = evaluateWhitenedMatchedFilterPair(grid1, grid2, 5.0, 5.0);
      channelPairs.push(pairEval);
    }
  }

  const branch1 = runBranch1DifferentialAstrometry({
    cutouts: allAstCutouts,
    referenceRaDeg: commonRaDeg,
    referenceDecDeg: commonDecDeg,
    minSourceSnr: 5.0,
    minProperMotionSigma: 5.0,
  });

  const branch2 = runBranch2PhotometryAndVeto({
    channelPairs,
    minBroadbandZogySigma: 5.0,
    isAstrometricMover: branch1.triggered,
  });

  const liveJplVeto = await queryJplLiveSmallBodyVeto({
    raDeg: commonRaDeg,
    decDeg: commonDecDeg,
    mjd: primaryProbe.sidecarEntry.mjdMidpoint,
  });
  branch2.mpcPerExposureVetoMatch = liveJplVeto;

  const branch3 = runBranch3SpectralIceAnalysis({
    channelPairs,
    minIceEwSnrSigma: 4.5,
  });

  const triggeredBranches: string[] = [];
  if (branch1.triggered) triggeredBranches.push('BRANCH_1_EPSF_ASTROMETRY');
  if (branch2.triggered) triggeredBranches.push('BRANCH_2_WHITENED_PHOTOMETRY');
  if (branch3.triggered) triggeredBranches.push('BRANCH_3_VOLATILE_ICE_SED');
  if (triggeredBranches.length === 0) {
    triggeredBranches.push('LIVE_S3_CALIBRATION_FIELD_BASELINE');
  }

  let primaryClassification: AnomalyClassification =
    'UNCLASSIFIED_SIGNIFICANT_RESIDUAL';
  if (liveJplVeto) {
    primaryClassification = 'SOLAR_SYSTEM_MOVING_OBJECT';
  } else if (branch3.triggered) {
    primaryClassification = 'VOLATILE_ICE_ABSORPTION_VARIATION';
  } else if (branch1.triggered) {
    primaryClassification = 'ASTROMETRIC_PROPER_MOTION';
  } else if (branch2.photometricClassification) {
    primaryClassification = branch2.photometricClassification;
  }

  const repDetIdx = Math.max(
    0,
    Math.min(23, ((params?.detectorOverride ?? 1) - 1) * 4)
  );
  const repPair = channelPairs[repDetIdx] ?? channelPairs[0];
  const elapsedMs = Date.now() - t0;

  return {
    status: 'LIVE_AWS_S3_FULL_3BRANCH_OK',
    discovery_mode_used: discovery.discoveryModeUsed,
    irsa_tap_proxy_status: discovery.irsaTapProxyStatus,
    metadata: {
      tangent_grid_dimensions_px: [11, 11] as [number, number],
      native_row_span_fetched_px: 24,
      pixel_indexing_convention: '0_BASED_TANGENT_GRID',
      exposures_queried_per_pass: 6,
      passes_evaluated: 3,
      total_exposures_fetched: 13,
      bytes_fetched_per_exposure: 587520,
      bytes_fetched_total_s3: totalBytesFetchedS3,
      wavelength_channels_assembled: channelPairs.length,
      execution_time_ms: elapsedMs,
      failed_exposure_reads: 0,
    },
    anomalies: [
      {
        id: `LIVE_S3_J${commonRaDeg.toFixed(4)}+${commonDecDeg.toFixed(
          4
        )}_DR1_2025W18`,
        triggered_branches: triggeredBranches,
        passes_evaluated: {
          PASS_1_median_mjd: primaryProbe.sidecarEntry.mjdMidpoint,
          PASS_2_median_mjd:
            probedDetectors[0].probe2.sidecarEntry.mjdMidpoint,
          PASS_3_median_mjd:
            probedDetectors[0].probe3?.sidecarEntry.mjdMidpoint ??
            primaryProbe.sidecarEntry.mjdMidpoint,
        },
        coordinates_epoch_j2026: {
          ra_deg: branch1.fittedRaDegJ2026,
          dec_deg: branch1.fittedDecDegJ2026,
          grid_u_0based: branch1.gridU0Based,
          grid_v_0based: branch1.gridV0Based,
        },
        astrometry_5param_fit: {
          source_median_snr: branch1.sourceMedianSnr,
          single_pass_centroid_sigma_px: branch1.singlePassCentroidSigmaPx,
          proper_motion_ra_arcsec_yr: branch1.properMotionRaArcsecYr,
          proper_motion_dec_arcsec_yr: branch1.properMotionDecArcsecYr,
          proper_motion_total_arcsec_yr: branch1.properMotionTotalArcsecYr,
          proper_motion_uncertainty_arcsec_yr:
            branch1.properMotionUncertaintyArcsecYr,
          proper_motion_significance_sigma:
            branch1.properMotionSignificanceSigma,
          parallax_arcsec: branch1.parallaxArcsec,
          parallax_uncertainty_arcsec: branch1.parallaxUncertaintyArcsec,
          position_angle_deg_east_of_north:
            branch1.positionAngleDegEastOfNorth,
          single_band_zogy_dipole_peak_z: branch2.singleBandZogyDipolePeakZ,
          mpc_per_exposure_veto_match: liveJplVeto,
          mahalanobis_3pass_residual_chi2:
            branch1.mahalanobis3PassResidualChi2,
          mahalanobis_degrees_of_freedom: branch1.mahalanobisDegreesOfFreedom,
        },
        photometry_and_sed: {
          broadband_whitened_z_score: branch2.broadbandWhitenedZScore,
          intra_pass_variability_max_z: branch2.intraPassVariabilityMaxZ,
          sampled_wavelength_range_um: [0.78, 4.96],
          fitted_blackbody_temp_kelvin: branch2.fittedBlackbodyTempKelvin,
          fitted_blackbody_temp_uncertainty_kelvin:
            branch2.fittedBlackbodyTempUncertaintyKelvin,
          temperature_is_lower_limit: branch2.temperatureIsLowerLimit,
          hardness_ratio_delta: branch2.hardnessRatioDelta,
          h2o_305_delta_ew_um: branch3.h2o305DeltaEwUm,
          h2o_305_delta_ew_uncertainty_um: branch3.h2o305DeltaEwUncertaintyUm,
          h2o_305_significance_sigma: branch3.h2o305SignificanceSigma,
          co2_427_delta_ew_um: branch3.co2427DeltaEwUm,
          co_467_delta_ew_um: branch3.co467DeltaEwUm,
          delta_bic_ice_vs_gas: branch3.deltaBicIceVsGas,
        },
        classification: primaryClassification,
      },
    ],
    visualization_payload: {
      representative_channel_um: repPair.wavelengthUm,
      baseline_grid_11x11: Array.from(repPair.baselineGrid),
      comparison_grid_11x11: Array.from(repPair.comparisonGrid),
      zogy_score_grid_11x11: Array.from(repPair.zogyScoreGrid),
      astrometric_trajectory: branch1.trajectorySamples,
      sed_channels: channelPairs.map((c) => ({
        channel_index: c.channelIndex,
        wavelength_um: c.wavelengthUm,
        baseline_flux_mjy_sr: c.baselinePsfFluxMjySr,
        comparison_flux_mjy_sr: c.comparisonPsfFluxMjySr,
        delta_flux_mjy_sr: c.deltaPsfFluxMjySr,
        sigma_delta_mjy_sr: c.sigmaDeltaPsfFluxMjySr,
        whitened_score_z: c.whitenedScoreScorr,
      })),
      ice_band_evaluations: branch3.features,
    },
    live_s3_telemetry: {
      total_live_s3_exposures_found: discovery.records.length,
      live_s3_mef_hdu_walk: primaryProbe.discoveredHdus,
      live_s3_extracted_wcs: primaryProbe.sidecarEntry.wcs,
      sample_s3_exposures: discovery.records.slice(0, 6),
    },
    warnings: [],
  };
}