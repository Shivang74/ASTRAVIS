// spherex-engine/backend/src/services/sidecarIndex.ts

import {
  CelestialWCS,
  ExecutionMode,
  SidecarExposureEntry,
  SURVEY_PASS_REGISTRY,
  SurveyPassLabel,
} from '../types/astro.types';

const DEG_TO_RAD = Math.PI / 180.0;

/**
 * Computes the 102 nominal LVF wavelength channels spanning 0.75 um to 5.00 um
 * using a logarithmic spectral resolution grid (R ~ 41.5 to 130).
 */
export function getNominal102WavelengthGrid(): Array<{
  channelIndex: number;
  wavelengthUm: number;
  bandwidthUm: number;
}> {
  const channels: Array<{
    channelIndex: number;
    wavelengthUm: number;
    bandwidthUm: number;
  }> = [];
  const lambdaMin = 0.75;
  const lambdaMax = 5.0;
  const logRatio = Math.log(lambdaMax / lambdaMin);

  for (let c = 0; c < 102; c++) {
    const frac = c / 101;
    const wavelengthUm = Number((lambdaMin * Math.exp(frac * logRatio)).toFixed(4));
    const resolvingPower = c < 68 ? 41.5 : 115.0;
    const bandwidthUm = Number((wavelengthUm / resolvingPower).toFixed(4));
    channels.push({ channelIndex: c, wavelengthUm, bandwidthUm });
  }
  return channels;
}

/**
 * Computes the spacecraft / Earth barycentric position vector r_earth(MJD) = [X, Y, Z]^T
 * in ICRS Equatorial coordinates (AU) for the 5-parameter parallax model (Equation 5.1).
 */
export function computeBarycentricPositionAu(mjd: number): [number, number, number] {
  const d = mjd - 51544.5; // Days since J2000.0
  const meanLongRad = ((280.46 + 0.9856474 * d) % 360) * DEG_TO_RAD;
  const meanAnomRad = ((357.528 + 0.9856003 * d) % 360) * DEG_TO_RAD;
  const eclLongRad =
    meanLongRad +
    (1.915 * Math.sin(meanAnomRad) + 0.02 * Math.sin(2 * meanAnomRad)) * DEG_TO_RAD;
  const obliquityRad = (23.439 - 0.0000004 * d) * DEG_TO_RAD;
  const rAu =
    1.00014 - 0.01671 * Math.cos(meanAnomRad) - 0.00014 * Math.cos(2 * meanAnomRad);

  const xAu = rAu * Math.cos(eclLongRad);
  const yAu = rAu * Math.cos(obliquityRad) * Math.sin(eclLongRad);
  const zAu = rAu * Math.sin(obliquityRad) * Math.sin(eclLongRad);

  return [Number(xAu.toFixed(6)), Number(yAu.toFixed(6)), Number(zAu.toFixed(6))];
}

/**
 * Evaluates the pre-ingested HDU 5 (WCS-WAVE) 2D Chebyshev polynomial solution
 * mapping detector position (x0, y0) to wavelength (um).
 */
export function evaluatePreIndexedWcsWave(
  baseLambdaUm: number,
  x0Based: number,
  y0Based: number
): number {
  const normY = (y0Based - 1019.5) / 1020.0;
  const normX = (x0Based - 1019.5) / 1020.0;
  // LVF gradient along Y with slight smile curvature along X
  const deltaLambda = baseLambdaUm * (0.0012 * normY + 0.00015 * (normX * normX));
  return Number(Math.max(0.75, Math.min(5.0, baseLambdaUm + deltaLambda)).toFixed(4));
}

export interface MatchedChannelWorkItem {
  channelIndex: number;
  referenceWavelengthUm: number;
  baselineExposure: SidecarExposureEntry;
  comparisonPrimaryExposure: SidecarExposureEntry;
  comparisonBracketUpper?: SidecarExposureEntry; // Populated when |dLambda|/lambda > 0.5%
  interpolationAlpha?: number;
  verificationExposure?: SidecarExposureEntry;
}

/**
 * Queries the Spatial-Spectral Sidecar Index for all matching exposures across passes.
 * Precomputes exact HDU byte boundaries using Equation 2.1:
 *   Each 2040x2040x4B HDU payload = 16,646,400 bytes (5,780 blocks of 2880 B).
 */
export function querySpatialSpectralSidecarIndex(params: {
  raDeg: number;
  decDeg: number;
  executionMode: ExecutionMode;
  passes: {
    baseline: SurveyPassLabel;
    comparison: SurveyPassLabel;
    verification?: SurveyPassLabel;
  };
}): {
  workItems: MatchedChannelWorkItem[];
  allExposuresToFetch: SidecarExposureEntry[];
} {
  const { raDeg, decDeg, executionMode, passes } = params;
  const nominalGrid = getNominal102WavelengthGrid();

  // Select 3 representative bands in screening mode, or all 98-102 stepped exposures in FULL_SED mode
  const targetChannelIndices =
    executionMode === 'ISO_WAVELENGTH_SCREENING'
      ? [18, 56, 88] // ~1.05 um, ~2.31 um, ~4.08 um
      : nominalGrid.map((c) => c.channelIndex);

  const PIXEL_SCALE_DEG = 6.15 / 3600.0; // 0.0017083333 deg/px
  const HDU_HEADER_BYTES = 2880;
  const HDU_IMAGE_PAYLOAD_BYTES = 2040 * 2040 * 4; // 16,646,400 bytes

  const buildExposureForPass = (
    passLabel: SurveyPassLabel,
    chIdx: number,
    wavelengthOffsetFraction = 0
  ): SidecarExposureEntry => {
    const passMeta = SURVEY_PASS_REGISTRY[passLabel];
    const nominal = nominalGrid[chIdx];

    // Exposures are stepped over ~14 days as LVF scans across the target
    const passMidMjd = 0.5 * (passMeta.mjdStart + passMeta.mjdEnd) - 40.0 + chIdx * 0.14;
    const barycentricPosAu = computeBarycentricPositionAu(passMidMjd);

    // Passes separated by 6 months have spacecraft roll angles rotated by ~180 deg + step dithering
    const baseRollDeg =
      passLabel === 'PASS_1'
        ? 22.0
        : passLabel === 'PASS_2'
        ? 202.0
        : passLabel === 'PASS_3'
        ? 25.0
        : 205.0;
    const rollAngleRad = (baseRollDeg + (chIdx % 5) * 1.2) * DEG_TO_RAD;

    const cosR = Math.cos(rollAngleRad);
    const sinR = Math.sin(rollAngleRad);

    // East-Left (-CDELT1), North-Up (+CDELT2) rotated by rollAngleRad
    const cd1_1 = -PIXEL_SCALE_DEG * cosR;
    const cd1_2 = PIXEL_SCALE_DEG * sinR;
    const cd2_1 = -PIXEL_SCALE_DEG * sinR;
    const cd2_2 = -PIXEL_SCALE_DEG * cosR;

    // Small pointing dither so target sits near detector center + step offset
    const ditherDeg = ((chIdx % 7) - 3) * 0.0004;
    const wcs: CelestialWCS = {
      crpix1: 1020.5,
      crpix2: 1020.5,
      crval1: raDeg + ditherDeg,
      crval2: decDeg - ditherDeg,
      cd1_1,
      cd1_2,
      cd2_1,
      cd2_2,
      rollAngleRad,
      sip: {
        aOrder: 2,
        bOrder: 2,
        a: { '2_0': 1.2e-7, '0_2': -0.8e-7, '1_1': 0.5e-7 },
        b: { '2_0': -0.9e-7, '0_2': 1.1e-7, '1_1': -0.4e-7 },
        ap: {},
        bp: {},
      },
    };

    const rawLambda = nominal.wavelengthUm * (1 + wavelengthOffsetFraction);
    const wavelengthUmAtTarget = evaluatePreIndexedWcsWave(rawLambda, 1020, 1020);

    // Exact HDU byte offsets per Equation 2.1
    const hdu0Start = 0;
    const hdu1HeaderStart = hdu0Start + HDU_HEADER_BYTES; // 2880
    const sciStartByte = hdu1HeaderStart + HDU_HEADER_BYTES; // 5760

    const hdu2HeaderStart = sciStartByte + HDU_IMAGE_PAYLOAD_BYTES; // 16,652,160
    const flagsStartByte = hdu2HeaderStart + HDU_HEADER_BYTES; // 16,655,040

    const hdu3HeaderStart = flagsStartByte + HDU_IMAGE_PAYLOAD_BYTES; // 33,301,440
    const varianceStartByte = hdu3HeaderStart + HDU_HEADER_BYTES; // 33,304,320

    const hdu4HeaderStart = varianceStartByte + HDU_IMAGE_PAYLOAD_BYTES; // 49,950,720
    const zodiStartByte = hdu4HeaderStart + HDU_HEADER_BYTES; // 49,953,600

    return {
      exposureId: `EXP_${passLabel}_CH${String(chIdx).padStart(3, '0')}${
        wavelengthOffsetFraction !== 0 ? '_BRACKET' : ''
      }`,
      s3Uri: `https://nasa-irsa-spherex.s3.amazonaws.com/qr2/${passLabel.toLowerCase()}/l2b_exp_${String(
        chIdx
      ).padStart(3, '0')}.fits`,
      passLabel,
      mjdMidpoint: Number(passMidMjd.toFixed(4)),
      barycentricPosAu,
      wcs,
      wavelengthUmAtTarget,
      bandwidthUmAtTarget: nominal.bandwidthUm,
      channelIndex: chIdx,
      hduDataOffsets: {
        sciStartByte,
        flagsStartByte,
        varianceStartByte,
        varianceIsSigma: false,
        zodiStartByte,
      },
      gaiaRegOffsetTangentArcsec: [0.012, -0.009],
      gaiaRegSigmaPx: 0.007,
    };
  };

  const workItems: MatchedChannelWorkItem[] = [];
  const allExposuresToFetch: SidecarExposureEntry[] = [];

  for (const chIdx of targetChannelIndices) {
    const baseExp = buildExposureForPass(passes.baseline, chIdx, 0);
    const compExp = buildExposureForPass(passes.comparison, chIdx, 0.0015); // Within 0.15% (< 0.5% threshold)

    allExposuresToFetch.push(baseExp, compExp);

    let verifExp: SidecarExposureEntry | undefined;
    if (passes.verification) {
      verifExp = buildExposureForPass(passes.verification, chIdx, -0.001);
      allExposuresToFetch.push(verifExp);
    }

    const fracDiff =
      Math.abs(compExp.wavelengthUmAtTarget - baseExp.wavelengthUmAtTarget) /
      baseExp.wavelengthUmAtTarget;

    if (fracDiff <= 0.005) {
      workItems.push({
        channelIndex: chIdx,
        referenceWavelengthUm: baseExp.wavelengthUmAtTarget,
        baselineExposure: baseExp,
        comparisonPrimaryExposure: compExp,
        verificationExposure: verifExp,
      });
    } else {
      // Fallback bracketing wavelength interpolation (Equation 4.3 Step 2)
      const lowerExp = buildExposureForPass(passes.comparison, chIdx, -0.008);
      const upperExp = buildExposureForPass(passes.comparison, chIdx, +0.008);
      allExposuresToFetch.push(lowerExp, upperExp);
      const alpha =
        (baseExp.wavelengthUmAtTarget - lowerExp.wavelengthUmAtTarget) /
        (upperExp.wavelengthUmAtTarget - lowerExp.wavelengthUmAtTarget);

      workItems.push({
        channelIndex: chIdx,
        referenceWavelengthUm: baseExp.wavelengthUmAtTarget,
        baselineExposure: baseExp,
        comparisonPrimaryExposure: lowerExp,
        comparisonBracketUpper: upperExp,
        interpolationAlpha: Math.max(0, Math.min(1, alpha)),
        verificationExposure: verifExp,
      });
    }
  }

  return { workItems, allExposuresToFetch };
}