// spherex-engine/backend/src/core/wcsProjection.ts

import { CelestialWCS, SipPolynomialCoeffs } from '../types/astro.types';

const DEG_TO_RAD = Math.PI / 180.0;
const RAD_TO_DEG = 180.0 / Math.PI;

function evalSipPoly(
  coeffs: Record<string, number>,
  order: number,
  u: number,
  v: number
): { val: number; dDu: number; dDv: number } {
  let val = 0;
  let dDu = 0;
  let dDv = 0;

  for (let p = 0; p <= order; p++) {
    for (let q = 0; q <= order - p; q++) {
      const c = coeffs[`${p}_${q}`];
      if (!c) continue;
      const uPow = p === 0 ? 1 : Math.pow(u, p);
      const vPow = q === 0 ? 1 : Math.pow(v, q);
      val += c * uPow * vPow;
      if (p > 0) dDu += c * p * Math.pow(u, p - 1) * vPow;
      if (q > 0) dDv += c * uPow * q * Math.pow(v, q - 1);
    }
  }
  return { val, dDu, dDv };
}

/**
 * Equation 4.1: Inverse Gnomonic (TAN) projection from sky (raDeg, decDeg)
 * to both 1-based FITS coordinates (xFits, yFits) and 0-based detector indices (x0, y0).
 */
export function skyToDetectorPixel(
  raDeg: number,
  decDeg: number,
  wcs: CelestialWCS
): {
  xiDeg: number;
  etaDeg: number;
  xFits1Based: number;
  yFits1Based: number;
  x0BasedContinuous: number;
  y0BasedContinuous: number;
  i0BasedInt: number;
  j0BasedInt: number;
} {
  const alpha = raDeg * DEG_TO_RAD;
  const delta = decDeg * DEG_TO_RAD;
  const alpha0 = wcs.crval1 * DEG_TO_RAD;
  const delta0 = wcs.crval2 * DEG_TO_RAD;

  const dAlpha = alpha - alpha0;
  const sinD = Math.sin(delta);
  const cosD = Math.cos(delta);
  const sinD0 = Math.sin(delta0);
  const cosD0 = Math.cos(delta0);

  const denom = sinD * sinD0 + cosD * cosD0 * Math.cos(dAlpha);
  if (denom <= 0) {
    throw new Error('Target coordinates lie >90 deg from tangent point.');
  }

  const xiDeg = RAD_TO_DEG * ((cosD * Math.sin(dAlpha)) / denom);
  const etaDeg =
    RAD_TO_DEG * ((sinD * cosD0 - cosD * sinD0 * Math.cos(dAlpha)) / denom);

  // Invert CD matrix: [u0, v0]^T = CD^-1 * [xiDeg, etaDeg]^T
  const det = wcs.cd1_1 * wcs.cd2_2 - wcs.cd1_2 * wcs.cd2_1;
  if (Math.abs(det) < 1e-18) {
    throw new Error('Singular WCS CD matrix.');
  }

  const u0 = (wcs.cd2_2 * xiDeg - wcs.cd1_2 * etaDeg) / det;
  const v0 = (-wcs.cd2_1 * xiDeg + wcs.cd1_1 * etaDeg) / det;

  let u = u0;
  let v = v0;

  // Apply SIP distortion inversion if present
  if (wcs.sip) {
    const sip: SipPolynomialCoeffs = wcs.sip;
    if (sip.apOrder !== undefined && sip.bpOrder !== undefined) {
      // Direct single-pass evaluation using inverse AP/BP polynomials
      const apEval = evalSipPoly(sip.ap, sip.apOrder, u0, v0);
      const bpEval = evalSipPoly(sip.bp, sip.bpOrder, u0, v0);
      u = u0 + apEval.val;
      v = v0 + bpEval.val;
    } else {
      // Exact 2D Newton-Raphson inversion using analytic Jacobian J_SIP
      for (let iter = 0; iter < 20; iter++) {
        const fA = evalSipPoly(sip.a, sip.aOrder, u, v);
        const gB = evalSipPoly(sip.b, sip.bOrder, u, v);

        const resU = u + fA.val - u0;
        const resV = v + gB.val - v0;

        if (Math.max(Math.abs(resU), Math.abs(resV)) < 1e-6) break;

        const j11 = 1 + fA.dDu;
        const j12 = fA.dDv;
        const j21 = gB.dDu;
        const j22 = 1 + gB.dDv;
        const jDet = j11 * j22 - j12 * j21;

        const du = (j22 * resU - j12 * resV) / jDet;
        const dv = (-j21 * resU + j11 * resV) / jDet;
        u -= du;
        v -= dv;
      }
    }
  }

  const xFits1Based = wcs.crpix1 + u;
  const yFits1Based = wcs.crpix2 + v;

  // Explicit 1-based FITS to 0-based detector coordinate conversion
  const x0BasedContinuous = xFits1Based - 1.0;
  const y0BasedContinuous = yFits1Based - 1.0;

  return {
    xiDeg,
    etaDeg,
    xFits1Based,
    yFits1Based,
    x0BasedContinuous,
    y0BasedContinuous,
    i0BasedInt: Math.round(xFits1Based) - 1,
    j0BasedInt: Math.round(yFits1Based) - 1,
  };
}

/**
 * Forward WCS transformation from 0-based detector coordinates (x0, y0)
 * back to sky coordinates (raDeg, decDeg) and local tangent offsets.
 */
export function detectorPixelToSky(
  x0Based: number,
  y0Based: number,
  wcs: CelestialWCS
): { raDeg: number; decDeg: number } {
  const xFits = x0Based + 1.0;
  const yFits = y0Based + 1.0;
  let u = xFits - wcs.crpix1;
  let v = yFits - wcs.crpix2;

  if (wcs.sip) {
    const fA = evalSipPoly(wcs.sip.a, wcs.sip.aOrder, u, v);
    const gB = evalSipPoly(wcs.sip.b, wcs.sip.bOrder, u, v);
    u += fA.val;
    v += gB.val;
  }

  const xiRad = (wcs.cd1_1 * u + wcs.cd1_2 * v) * DEG_TO_RAD;
  const etaRad = (wcs.cd2_1 * u + wcs.cd2_2 * v) * DEG_TO_RAD;

  const alpha0 = wcs.crval1 * DEG_TO_RAD;
  const delta0 = wcs.crval2 * DEG_TO_RAD;

  const rho = Math.hypot(xiRad, etaRad);
  if (rho < 1e-15) {
    return { raDeg: wcs.crval1, decDeg: wcs.crval2 };
  }

  const c = Math.atan(rho);
  const sinC = Math.sin(c);
  const cosC = Math.cos(c);
  const sinD0 = Math.sin(delta0);
  const cosD0 = Math.cos(delta0);

  const delta = Math.asin(cosC * sinD0 + (etaRad * sinC * cosD0) / rho);
  const alpha =
    alpha0 +
    Math.atan2(xiRad * sinC, rho * cosD0 * cosC - etaRad * sinD0 * sinC);

  return {
    raDeg: ((alpha * RAD_TO_DEG) % 360 + 360) % 360,
    decDeg: delta * RAD_TO_DEG,
  };
}

/**
 * Projects a North-Up, East-Left local tangent-plane coordinate (xiArcsec, etaArcsec)
 * centered on (centerRaDeg, centerDecDeg) directly onto 0-based detector coordinates (x0, y0).
 */
export function localTangentOffsetToDetectorPixel(
  centerRaDeg: number,
  centerDecDeg: number,
  dXiArcsecEast: number,
  dEtaArcsecNorth: number,
  wcs: CelestialWCS
): { x0Based: number; y0Based: number } {
  const decRad = centerDecDeg * DEG_TO_RAD;
  const raDeg =
    centerRaDeg + dXiArcsecEast / 3600.0 / Math.max(0.05, Math.cos(decRad));
  const decDeg = centerDecDeg + dEtaArcsecNorth / 3600.0;
  const proj = skyToDetectorPixel(raDeg, decDeg, wcs);
  return {
    x0Based: proj.x0BasedContinuous,
    y0Based: proj.y0BasedContinuous,
  };
}

/**
 * Equation 4.1 (Step 4): Computes the rotated-grid & kernel-padded contiguous row span
 * and exact HTTP byte-range headers for a 2040x2040 Float32/Int32 HDU.
 */
export function computeContiguousRowSpanByteRange(params: {
  j0BasedCenter: number;
  rollAngleRad: number;
  hduDataStartByte: number;
  detectorWidth?: number;
  detectorHeight?: number;
  gridSidePx?: number;
  lanczosRadiusPx?: number;
  sipMarginPx?: number;
}): {
  rDetRowsHalf: number;
  jMin0Based: number;
  jMax0Based: number;
  rowCount: number;
  byteStart: number;
  byteEnd: number;
  byteLength: number;
  rangeHeaderValue: string;
} {
  const {
    j0BasedCenter,
    rollAngleRad,
    hduDataStartByte,
    detectorWidth = 2040,
    detectorHeight = 2040,
    gridSidePx = 11,
    lanczosRadiusPx = 3,
    sipMarginPx = 1,
  } = params;

  const rotFactor =
    Math.abs(Math.cos(rollAngleRad)) + Math.abs(Math.sin(rollAngleRad));
  const rDetRowsHalf = Math.ceil(
    (gridSidePx / 2.0) * rotFactor + lanczosRadiusPx + sipMarginPx
  ); // 12 rows at 45 deg -> 24 rows total

  const totalRows = 2 * rDetRowsHalf;
  let jMin0Based = Math.max(0, j0BasedCenter - rDetRowsHalf);
  let jMax0Based = jMin0Based + totalRows - 1;

  if (jMax0Based >= detectorHeight) {
    jMax0Based = detectorHeight - 1;
    jMin0Based = Math.max(0, jMax0Based - totalRows + 1);
  }

  const rowCount = jMax0Based - jMin0Based + 1;
  const byteStart = hduDataStartByte + jMin0Based * detectorWidth * 4;
  const byteEnd = hduDataStartByte + (jMax0Based + 1) * detectorWidth * 4 - 1;
  const byteLength = byteEnd - byteStart + 1;

  return {
    rDetRowsHalf,
    jMin0Based,
    jMax0Based,
    rowCount,
    byteStart,
    byteEnd,
    byteLength,
    rangeHeaderValue: `bytes=${byteStart}-${byteEnd}`,
  };
}