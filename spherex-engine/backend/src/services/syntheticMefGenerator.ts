// spherex-engine/backend/src/services/syntheticMefGenerator.ts

import { localTangentOffsetToDetectorPixel } from '../core/wcsProjection';
import {
  QUALITY_FLAG_BITS,
  SidecarExposureEntry,
  ViewportAnalysisRequest,
} from '../types/astro.types';

/**
 * Evaluates Planck blackbody specific intensity B_nu(lambda, T) in MJy/sr (normalized scale).
 * Uses Equation 5.2: B_nu(lambda, T) proportional to (1 / lambda^3) / (exp(14387.77 / (lambda * T)) - 1).
 */
export function evaluateNormalizedPlanckBnu(
  wavelengthUm: number,
  tempKelvin: number,
  amplitudeScale = 1.0
): number {
  const x = 14387.7688 / (wavelengthUm * tempKelvin);
  if (x > 80) return 0;
  const bNu = (1.0 / Math.pow(wavelengthUm, 3)) * (1.0 / (Math.exp(x) - 1.0));
  // Normalize relative to 1.6 um at 2500 K
  const refX = 14387.7688 / (1.6 * 2500.0);
  const refBnu = (1.0 / Math.pow(1.6, 3)) * (1.0 / (Math.exp(refX) - 1.0));
  return amplitudeScale * (bNu / refBnu);
}

/**
 * Deterministic pseudo-random normal generator (Box-Muller) seeded by exposure & pixel index.
 */
function seededGaussian(seed: number): number {
  const s1 = Math.sin(seed * 12.9898 + 78.233) * 43758.5453;
  const u1 = Math.max(1e-7, s1 - Math.floor(s1));
  const s2 = Math.sin((seed + 1) * 93.9898 + 12.456) * 23421.6312;
  const u2 = Math.max(1e-7, s2 - Math.floor(s2));
  return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
}

/**
 * Synthesizes 24-row Big-Endian IEEE-754 binary buffers (195,840 bytes each) for
 * HDU 1 (SCI), HDU 2 (FLAGS), HDU 3 (VARIANCE), and HDU 4 (ZODI).
 */
export function synthesize24RowBigEndianHduSlices(params: {
  exposure: SidecarExposureEntry;
  centerRaDeg: number;
  centerDecDeg: number;
  jMin0Based: number;
  rowCount: number;
  preset: NonNullable<ViewportAnalysisRequest['simulation_preset']>;
}): {
  sciBigEndian: Buffer;
  flagsBigEndian: Buffer;
  varianceBigEndian: Buffer;
  zodiBigEndian: Buffer;
} {
  const { exposure, centerRaDeg, centerDecDeg, jMin0Based, rowCount, preset } = params;
  const width = 2040;
  const totalPixels = rowCount * width; // 24 * 2040 = 48,960 pixels (195,840 bytes)

  const sciBigEndian = Buffer.allocUnsafe(totalPixels * 4);
  const flagsBigEndian = Buffer.allocUnsafe(totalPixels * 4);
  const varianceBigEndian = Buffer.allocUnsafe(totalPixels * 4);
  const zodiBigEndian = Buffer.allocUnsafe(totalPixels * 4);

  const lambda = exposure.wavelengthUmAtTarget;
  const chIdx = exposure.channelIndex;

  // Compute source tangent-plane offset (dXiEastArcsec, dEtaNorthArcsec) and flux
  const tYears = (exposure.mjdMidpoint - 60842.318) / 365.25;
  const [X, Y, Z] = exposure.barycentricPosAu;
  const alpha0Rad = (centerRaDeg * Math.PI) / 180.0;
  const delta0Rad = (centerDecDeg * Math.PI) / 180.0;

  // Parallax factors P_alpha, P_delta (Equation 5.1 Step 2)
  const pAlpha = X * Math.sin(alpha0Rad) - Y * Math.cos(alpha0Rad);
  const pDelta =
    X * Math.cos(alpha0Rad) * Math.sin(delta0Rad) +
    Y * Math.sin(alpha0Rad) * Math.sin(delta0Rad) -
    Z * Math.cos(delta0Rad);

  let dXiArcsec = 1.1; // Base offset from cutout center (+0.18 px)
  let dEtaArcsec = 2.5; // Base offset (+0.41 px)
  let sourceFluxPeakMjySr = 12.0;

  if (preset === 'BROWN_DWARF_PROPER_MOTION') {
    // T = 1650 K L/T transition dwarf with mu_alpha* = -0.44"/yr, mu_delta = +0.73"/yr, parallax = 0.062"
    const muAlphaArcsecYr = -0.44;
    const muDeltaArcsecYr = 0.73;
    const parallaxArcsec = 0.062;
    dXiArcsec += muAlphaArcsecYr * tYears + parallaxArcsec * pAlpha;
    dEtaArcsec += muDeltaArcsecYr * tYears + parallaxArcsec * pDelta;
    sourceFluxPeakMjySr = evaluateNormalizedPlanckBnu(lambda, 1650, 14.5);
  } else if (preset === 'HOT_STELLAR_FLARE') {
    // Stationary M-dwarf (T=3100 K) that erupts with a T=14,500 K Rayleigh-Jeans flare in PASS_2
    const quiescentFlux = evaluateNormalizedPlanckBnu(lambda, 3100, 10.0);
    if (exposure.passLabel === 'PASS_2') {
      const flareFlux = evaluateNormalizedPlanckBnu(lambda, 14500, 18.0);
      sourceFluxPeakMjySr = quiescentFlux + flareFlux;
    } else {
      sourceFluxPeakMjySr = quiescentFlux;
    }
  } else if (preset === 'WATER_ICE_OUTGASSING') {
    // Stationary disk/comet with deep 3.05 um solid H2O ice trough appearing in PASS_2
    const baseContinuum = evaluateNormalizedPlanckBnu(lambda, 2200, 15.0);
    if (exposure.passLabel === 'PASS_2') {
      const h2oProfile = Math.exp(-Math.pow((lambda - 3.05) / 0.14, 2));
      const co2Profile = Math.exp(-Math.pow((lambda - 4.27) / 0.05, 2));
      const opticalDepth = 0.55 * h2oProfile + 0.35 * co2Profile;
      sourceFluxPeakMjySr = baseContinuum * Math.exp(-opticalDepth);
    } else {
      sourceFluxPeakMjySr = baseContinuum;
    }
  } else if (preset === 'MAIN_BELT_ASTEROID_VETO') {
    // Fast-moving asteroid present ONLY in PASS_2 channels 55-57 (intra-pass tracklet + MPC match)
    const quiescentBgStar = evaluateNormalizedPlanckBnu(lambda, 4200, 4.0);
    if (exposure.passLabel === 'PASS_2' && chIdx >= 55 && chIdx <= 57) {
      sourceFluxPeakMjySr = quiescentBgStar + 28.0;
    } else {
      sourceFluxPeakMjySr = quiescentBgStar;
    }
  }

  // Project true source sky position to 0-based continuous detector coordinates
  const detPos = localTangentOffsetToDetectorPixel(
    centerRaDeg,
    centerDecDeg,
    dXiArcsec,
    dEtaArcsec,
    exposure.wcs
  );

  const psfSigmaPx = 1.35 / 2.35482; // FWHM = 1.35 px -> sigma_PSF = 0.5733 px
  const zodiBase = 1.85 + 0.35 * Math.sin(lambda);
  const pixelVarianceBase = 0.04; // sigma = 0.20 MJy/sr

  const centerCol = Math.round(detPos.x0Based);

  for (let r = 0; r < rowCount; r++) {
    const detY = jMin0Based + r;
    const dy = detY - detPos.y0Based;
    const dy2 = dy * dy;

    for (let c = 0; c < width; c++) {
      const idx = r * width + c;
      const byteOffset = idx * 4;

      // Only compute detailed PSF + noise within +/- 20 columns of target to keep synthesis < 1 ms
      if (Math.abs(c - centerCol) <= 20) {
        const dx = c - detPos.x0Based;
        const r2 = dx * dx + dy2;
        const psfVal =
          r2 < 25.0 ? Math.exp(-0.5 * (r2 / (psfSigmaPx * psfSigmaPx))) : 0.0;

        const seed =
          exposure.mjdMidpoint * 1000 + detY * 2040 + c + exposure.channelIndex * 17;
        const noise = seededGaussian(seed) * Math.sqrt(pixelVarianceBase);

        // Set Bit 2 (NONLINEAR_RESPONSE = 0x0004) on one outer wing pixel to test variance inflation
        const flagVal =
          r === 3 && c === centerCol + 4 ? QUALITY_FLAG_BITS.NONLINEAR_RESPONSE : 0;

        const sciVal = zodiBase + sourceFluxPeakMjySr * psfVal + noise;

        sciBigEndian.writeFloatBE(sciVal, byteOffset);
        flagsBigEndian.writeUInt32BE(flagVal, byteOffset);
        varianceBigEndian.writeFloatBE(pixelVarianceBase, byteOffset);
        zodiBigEndian.writeFloatBE(zodiBase, byteOffset);
      } else {
        sciBigEndian.writeFloatBE(zodiBase, byteOffset);
        flagsBigEndian.writeUInt32BE(0, byteOffset);
        varianceBigEndian.writeFloatBE(pixelVarianceBase, byteOffset);
        zodiBigEndian.writeFloatBE(zodiBase, byteOffset);
      }
    }
  }

  return {
    sciBigEndian,
    flagsBigEndian,
    varianceBigEndian,
    zodiBigEndian,
  };
}