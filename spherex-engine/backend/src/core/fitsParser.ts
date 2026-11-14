// spherex-engine/backend/src/core/fitsParser.ts

import {
  CelestialWCS,
  FitsStructuralHeader,
  QUALITY_FLAG_BITS,
  SipPolynomialCoeffs,
} from '../types/astro.types';

export class EndCardNotFoundError extends Error {
  constructor(public readonly bytesScanned: number) {
    super(
      `FITS END card not found within ${bytesScanned} bytes; retry with expanded byte range.`
    );
    this.name = 'EndCardNotFoundError';
  }
}

/**
 * Parses an 80-char FITS value field, preserving slashes inside single-quoted strings,
 * unescaping doubled single quotes (''), and normalizing Fortran 'D' exponents.
 */
export function parseFitsValueField(rawField: string): string | number | boolean {
  const trimmed = rawField.trim();
  if (trimmed.startsWith("'")) {
    let strVal = '';
    for (let i = 1; i < trimmed.length; i++) {
      if (trimmed[i] === "'") {
        if (trimmed[i + 1] === "'") {
          strVal += "'";
          i++;
          continue;
        }
        return strVal.trimEnd();
      }
      strVal += trimmed[i];
    }
    return strVal.trimEnd();
  }

  const valuePart = trimmed.split('/')[0].trim();
  if (valuePart === 'T') return true;
  if (valuePart === 'F') return false;

  const normalizedNum = valuePart.replace(/([0-9])([dD])([+-]?[0-9])/g, '$1e$3');
  const num = Number(normalizedNum);
  if (!Number.isNaN(num) && normalizedNum !== '') return num;
  return valuePart;
}

/**
 * Stage 1: Parses structural FITS metadata and computes exact HDU byte boundaries.
 * Safe for all HDU types including HDU 0 (NAXIS=0) and HDU 5 (BINTABLE).
 */
export function parseFitsStructuralHeader(
  buffer: Buffer,
  hduOffset = 0
): FitsStructuralHeader {
  const CARD_SIZE = 80;
  const BLOCK_SIZE = 2880;
  const records: Record<string, string | number | boolean> = {};
  let lastStringKey = '';
  let cardCount = 0;
  let endFound = false;

  for (let offset = 0; offset + CARD_SIZE <= buffer.length; offset += CARD_SIZE) {
    cardCount++;
    const card = buffer.toString('ascii', offset, offset + CARD_SIZE);
    const keyword = card.slice(0, 8).trim();

    if (keyword === 'END') {
      endFound = true;
      break;
    }

    if (
      keyword === 'CONTINUE' &&
      lastStringKey &&
      typeof records[lastStringKey] === 'string'
    ) {
      const contVal = String(parseFitsValueField(card.slice(10)));
      records[lastStringKey] =
        String(records[lastStringKey]).replace(/&$/, '') + contVal;
      continue;
    }

    if (card.slice(8, 10) === '= ') {
      const parsed = parseFitsValueField(card.slice(10));
      records[keyword] = parsed;
      lastStringKey = typeof parsed === 'string' ? keyword : '';
    }
  }

  if (!endFound) {
    throw new EndCardNotFoundError(buffer.length);
  }

  const requireNum = (key: string): number => {
    if (typeof records[key] !== 'number') {
      throw new Error(`Missing mandatory FITS structural keyword: ${key}`);
    }
    return records[key] as number;
  };

  const bitpix = requireNum('BITPIX');
  const naxis = requireNum('NAXIS');
  const naxes: number[] = [];
  let elementCount = naxis > 0 ? 1 : 0;

  for (let k = 1; k <= naxis; k++) {
    const axisLen = requireNum(`NAXIS${k}`);
    naxes.push(axisLen);
    elementCount *= axisLen;
  }

  const pcount =
    typeof records['PCOUNT'] === 'number' ? (records['PCOUNT'] as number) : 0;
  const gcount =
    typeof records['GCOUNT'] === 'number' ? (records['GCOUNT'] as number) : 1;

  const headerBytes = Math.ceil((cardCount * CARD_SIZE) / BLOCK_SIZE) * BLOCK_SIZE;
  const dataStartByte = hduOffset + headerBytes;
  const rawDataBytes =
    (Math.abs(bitpix) / 8) * gcount * (pcount + elementCount);
  const paddedDataBytes = Math.ceil(rawDataBytes / BLOCK_SIZE) * BLOCK_SIZE;

  return {
    xtension: String(
      records['XTENSION'] ?? (records['SIMPLE'] === true ? 'PRIMARY' : 'UNKNOWN')
    ),
    extname: String(records['EXTNAME'] ?? 'PRIMARY'),
    bitpix,
    naxis,
    naxes,
    pcount,
    gcount,
    bscale:
      typeof records['BSCALE'] === 'number' ? (records['BSCALE'] as number) : 1.0,
    bzero:
      typeof records['BZERO'] === 'number' ? (records['BZERO'] as number) : 0.0,
    dataStartByte,
    nextHduStartByte: dataStartByte + paddedDataBytes,
    records,
  };
}

/**
 * Stage 2: Conditionally extracts 2D Celestial WCS parameters (including SIP polynomials)
 * only for 2D image HDUs. Returns null for NAXIS=0 (PRIMARY) or BINTABLE extensions.
 */
export function extractCelestialWCS(
  header: FitsStructuralHeader
): CelestialWCS | null {
  if (header.naxis < 2 || header.xtension === 'BINTABLE') return null;
  const r = header.records;
  const ctype1 = String(r['CTYPE1'] ?? '');
  if (!ctype1.startsWith('RA--')) return null;

  const req = (key: string): number => {
    if (typeof r[key] !== 'number') {
      throw new Error(`Missing mandatory WCS keyword ${key} in ${header.extname}`);
    }
    return r[key] as number;
  };

  let cd1_1 = r['CD1_1'] as number | undefined;
  let cd1_2 = (r['CD1_2'] as number | undefined) ?? 0;
  let cd2_1 = (r['CD2_1'] as number | undefined) ?? 0;
  let cd2_2 = r['CD2_2'] as number | undefined;

  if (cd1_1 === undefined || cd2_2 === undefined) {
    const cdelt1 = req('CDELT1');
    const cdelt2 = req('CDELT2');
    const pc1_1 = (r['PC1_1'] as number | undefined) ?? 1;
    const pc1_2 = (r['PC1_2'] as number | undefined) ?? 0;
    const pc2_1 = (r['PC2_1'] as number | undefined) ?? 0;
    const pc2_2 = (r['PC2_2'] as number | undefined) ?? 1;
    cd1_1 = cdelt1 * pc1_1;
    cd1_2 = cdelt1 * pc1_2;
    cd2_1 = cdelt2 * pc2_1;
    cd2_2 = cdelt2 * pc2_2;
  }

  let sip: SipPolynomialCoeffs | undefined;
  if (typeof r['A_ORDER'] === 'number' && typeof r['B_ORDER'] === 'number') {
    const aOrder = r['A_ORDER'] as number;
    const bOrder = r['B_ORDER'] as number;
    const apOrder =
      typeof r['AP_ORDER'] === 'number' ? (r['AP_ORDER'] as number) : undefined;
    const bpOrder =
      typeof r['BP_ORDER'] === 'number' ? (r['BP_ORDER'] as number) : undefined;

    const a: Record<string, number> = {};
    const b: Record<string, number> = {};
    const ap: Record<string, number> = {};
    const bp: Record<string, number> = {};

    for (let p = 0; p <= Math.max(aOrder, bOrder, apOrder ?? 0, bpOrder ?? 0); p++) {
      for (let q = 0; q <= Math.max(aOrder, bOrder, apOrder ?? 0, bpOrder ?? 0); q++) {
        const key = `${p}_${q}`;
        if (typeof r[`A_${key}`] === 'number') a[key] = r[`A_${key}`] as number;
        if (typeof r[`B_${key}`] === 'number') b[key] = r[`B_${key}`] as number;
        if (typeof r[`AP_${key}`] === 'number') ap[key] = r[`AP_${key}`] as number;
        if (typeof r[`BP_${key}`] === 'number') bp[key] = r[`BP_${key}`] as number;
      }
    }
    sip = { aOrder, bOrder, apOrder, bpOrder, a, b, ap, bp };
  }

  return {
    crpix1: req('CRPIX1'),
    crpix2: req('CRPIX2'),
    crval1: req('CRVAL1'),
    crval2: req('CRVAL2'),
    cd1_1,
    cd1_2,
    cd2_1,
    cd2_2,
    rollAngleRad: Math.atan2(-cd1_2, cd2_2),
    sip,
  };
}

/**
 * Fast Big-Endian IEEE-754 Float32 hydration using native C++ Buffer.swap32().
 */
export function decodeBigEndianFloat32(
  buffer: Buffer,
  bscale = 1.0,
  bzero = 0.0
): Float32Array {
  const aligned = Buffer.from(buffer);
  aligned.swap32();
  const rawFloats = new Float32Array(
    aligned.buffer,
    aligned.byteOffset,
    aligned.byteLength / 4
  );
  if (bscale === 1.0 && bzero === 0.0) return rawFloats;

  const scaled = new Float32Array(rawFloats.length);
  for (let i = 0; i < rawFloats.length; i++) {
    scaled[i] = Number.isNaN(rawFloats[i]) ? NaN : bzero + bscale * rawFloats[i];
  }
  return scaled;
}

/**
 * Decodes Big-Endian 32-bit integer flag arrays (HDU 2: FLAGS), supporting unsigned BZERO=2^31.
 */
export function decodeBigEndianFlags32(buffer: Buffer, bzero = 0): Uint32Array {
  const aligned = Buffer.from(buffer);
  aligned.swap32();
  const intView = new Int32Array(
    aligned.buffer,
    aligned.byteOffset,
    aligned.byteLength / 4
  );
  const out = new Uint32Array(intView.length);
  for (let i = 0; i < intView.length; i++) {
    out[i] = (intView[i] + bzero) >>> 0;
  }
  return out;
}

/**
 * Applies Equation 4.2: Zodiacal subtraction, Bit-2 (0x0004) 2.25x variance inflation,
 * ERR-to-VARIANCE squaring, and strict bitmask filtering (default 0xFFFB).
 */
export function calibrateAndMaskCutout(params: {
  rawSci: Float32Array;
  rawFlags: Uint32Array;
  rawVarOrErr: Float32Array;
  rawZodi: Float32Array;
  varianceIsSigma: boolean;
  rejectBitmask?: number;
}): {
  fluxMjySr: Float32Array;
  varianceMjySr2: Float32Array;
  validMask: Uint8Array;
  saturatedFraction: number;
  hasPersistenceFlag: boolean;
  hasGhostFlag: boolean;
} {
  const {
    rawSci,
    rawFlags,
    rawVarOrErr,
    rawZodi,
    varianceIsSigma,
    rejectBitmask = QUALITY_FLAG_BITS.DEFAULT_REJECT_MASK,
  } = params;

  const n = rawSci.length;
  const fluxMjySr = new Float32Array(n);
  const varianceMjySr2 = new Float32Array(n);
  const validMask = new Uint8Array(n);

  let saturatedCount = 0;
  let hasPersistenceFlag = false;
  let hasGhostFlag = false;

  for (let i = 0; i < n; i++) {
    const flag = rawFlags[i];
    if ((flag & QUALITY_FLAG_BITS.OVERFLOW_SATURATED) !== 0) saturatedCount++;
    if ((flag & QUALITY_FLAG_BITS.PERSISTENCE_ARTIFACT) !== 0) hasPersistenceFlag = true;
    if ((flag & QUALITY_FLAG_BITS.OPTICAL_GHOST_OR_STREAK) !== 0) hasGhostFlag = true;

    const fRaw = rawSci[i];
    const zodi = rawZodi[i];
    let vRaw = rawVarOrErr[i];
    if (varianceIsSigma) vRaw = vRaw * vRaw;

    // Bit 2 (0x0004, NONLINEAR_RESPONSE) inflates variance by 1.5^2 = 2.25x
    if ((flag & QUALITY_FLAG_BITS.NONLINEAR_RESPONSE) !== 0) {
      vRaw *= 2.25;
    }

    const isCleanFlag = (flag & rejectBitmask) === 0;
    const isFiniteVal =
      Number.isFinite(fRaw) && Number.isFinite(zodi) && Number.isFinite(vRaw) && vRaw > 0;

    if (isCleanFlag && isFiniteVal) {
      fluxMjySr[i] = fRaw - zodi;
      varianceMjySr2[i] = vRaw;
      validMask[i] = 1;
    } else {
      fluxMjySr[i] = 0;
      varianceMjySr2[i] = Infinity;
      validMask[i] = 0;
    }
  }

  return {
    fluxMjySr,
    varianceMjySr2,
    validMask,
    saturatedFraction: n > 0 ? saturatedCount / n : 0,
    hasPersistenceFlag,
    hasGhostFlag,
  };
}