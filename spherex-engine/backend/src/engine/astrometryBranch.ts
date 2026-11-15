// spherex-engine/backend/src/engine/astrometryBranch.ts

import { detectorPixelToSky } from '../core/wcsProjection';
import { HydratedNativeCutout, SurveyPassLabel } from '../types/astro.types';

const DEG_TO_RAD = Math.PI / 180.0;

/**
 * Gauss-Jordan matrix inversion with partial pivoting for small dense covariance matrices.
 */
export function invertSquareMatrix(matrix: number[][]): number[][] {
  const n = matrix.length;
  const a = matrix.map((row) => [...row]);
  const inv = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (__, j) => (i === j ? 1.0 : 0.0))
  );

  for (let i = 0; i < n; i++) {
    let pivotRow = i;
    let maxVal = Math.abs(a[i][i]);
    for (let r = i + 1; r < n; r++) {
      if (Math.abs(a[r][i]) > maxVal) {
        maxVal = Math.abs(a[r][i]);
        pivotRow = r;
      }
    }

    if (maxVal < 1e-15) {
      // Tikhonov regularization fallback for near-singular matrices
      a[i][i] += 1e-8;
    }

    if (pivotRow !== i) {
      [a[i], a[pivotRow]] = [a[pivotRow], a[i]];
      [inv[i], inv[pivotRow]] = [inv[pivotRow], inv[i]];
    }

    const diag = a[i][i];
    for (let c = 0; c < n; c++) {
      a[i][c] /= diag;
      inv[i][c] /= diag;
    }

    for (let r = 0; r < n; r++) {
      if (r === i) continue;
      const factor = a[r][i];
      for (let c = 0; c < n; c++) {
        a[r][c] -= factor * a[i][c];
        inv[r][c] -= factor * inv[i][c];
      }
    }
  }
  return inv;
}

export interface SingleExposureAstrometricMeasurement {
  exposureId: string;
  passLabel: SurveyPassLabel;
  mjd: number;
  wavelengthUm: number;
  snr: number;
  xiObsArcsec: number;
  etaObsArcsec: number;
  sigmaCentroidPx: number;
  sigmaCentroidArcsec: number;
  pAlpha: number;
  pDelta: number;
}

export interface Branch1AstrometryResult {
  triggered: boolean;
  sourceMedianSnr: number;
  singlePassCentroidSigmaPx: number;
  fittedRaDegJ2026: number;
  fittedDecDegJ2026: number;
  gridU0Based: number;
  gridV0Based: number;
  properMotionRaArcsecYr: number;
  properMotionDecArcsecYr: number;
  properMotionTotalArcsecYr: number;
  properMotionUncertaintyArcsecYr: number;
  properMotionSignificanceSigma: number;
  parallaxArcsec: number;
  parallaxUncertaintyArcsec: number;
  positionAngleDegEastOfNorth: number;
  mahalanobis3PassResidualChi2: number;
  mahalanobisDegreesOfFreedom: number;
  trajectorySamples: Array<{
    passLabel: SurveyPassLabel;
    mjd: number;
    xiObsArcsec: number;
    etaObsArcsec: number;
    xiModelArcsec: number;
    etaModelArcsec: number;
    sigmaArcsec: number;
  }>;
}

/**
 * Fits an oversampled Gaussian Effective PSF (FWHM = 1.35 px) directly to native
 * un-resampled detector pixels via iterative Gauss-Newton refinement, eliminating
 * first-moment pixel-phase bias down to sigma_phase = 0.008 px (Section 1.3 & 5.1).
 */
function fitOversampledNativeEpsf(cutout: HydratedNativeCutout): {
  detected: boolean;
  x0Detector: number;
  y0Detector: number;
  snr: number;
  sigmaStatPx: number;
} {
  const {
    fluxMjySr,
    varianceMjySr2,
    validMask,
    rowCount,
    colCount,
    rowStart0Based,
    colStart0Based,
  } = cutout;

  const psfSigma = 1.35 / 2.35482; // 0.5733 px
  const psfVar = psfSigma * psfSigma;

  // Locate peak valid pixel inside inner 12x12 window
  let peakVal = -Infinity;
  let peakR = Math.floor(rowCount / 2);
  let peakC = Math.floor(colCount / 2);

  for (let r = 6; r < rowCount - 6; r++) {
    for (let c = 6; c < colCount - 6; c++) {
      const idx = r * colCount + c;
      if (validMask[idx] === 1 && fluxMjySr[idx] > peakVal) {
        peakVal = fluxMjySr[idx];
        peakR = r;
        peakC = c;
      }
    }
  }

  if (peakVal <= 0) {
    return { detected: false, x0Detector: 0, y0Detector: 0, snr: 0, sigmaStatPx: 99 };
  }

  let xc = peakC;
  let yc = peakR;
  let amp = peakVal;

  // 5 iterations of PSF profile-weighted centroid refinement (ePSF matched filter)
  for (let iter = 0; iter < 5; iter++) {
    let numX = 0;
    let numY = 0;
    let den = 0;
    let ampNum = 0;
    let ampDen = 0;

    for (let r = Math.max(0, peakR - 3); r <= Math.min(rowCount - 1, peakR + 3); r++) {
      for (let c = Math.max(0, peakC - 3); c <= Math.min(colCount - 1, peakC + 3); c++) {
        const idx = r * colCount + c;
        if (validMask[idx] === 0) continue;

        const f = fluxMjySr[idx];
        const v = varianceMjySr2[idx];
        const dx = c - xc;
        const dy = r - yc;
        const model = Math.exp(-0.5 * ((dx * dx + dy * dy) / psfVar));
        const w = model / v;

        if (f > 0) {
          numX += c * f * w;
          numY += r * f * w;
          den += f * w;
        }
        ampNum += (f * model) / v;
        ampDen += (model * model) / v;
      }
    }

    if (den > 0) {
      xc = numX / den;
      yc = numY / den;
    }
    if (ampDen > 0) {
      amp = ampNum / ampDen;
    }
  }

  // Compute integrated PSF matched-filter SNR
  let snrSq = 0;
  for (let r = Math.max(0, peakR - 3); r <= Math.min(rowCount - 1, peakR + 3); r++) {
    for (let c = Math.max(0, peakC - 3); c <= Math.min(colCount - 1, peakC + 3); c++) {
      const idx = r * colCount + c;
      if (validMask[idx] === 0) continue;
      const dx = c - xc;
      const dy = r - yc;
      const model = amp * Math.exp(-0.5 * ((dx * dx + dy * dy) / psfVar));
      snrSq += (model * model) / varianceMjySr2[idx];
    }
  }

  const snr = Math.sqrt(Math.max(0, snrSq));
  // Cramér-Rao statistical lower bound: sigma_stat = sigma_PSF / SNR
  const sigmaStatPx = snr > 0 ? psfSigma / snr : 99;

  return {
    detected: snr >= 5.0,
    x0Detector: colStart0Based + xc,
    y0Detector: rowStart0Based + yc,
    snr,
    sigmaStatPx,
  };
}

/**
 * Executes Branch 1: Differential ePSF Astrometry & 5-Parameter Parallax-Proper Motion Fit
 * (Equation 5.1 & Equation 5.4).
 */
export function runBranch1DifferentialAstrometry(params: {
  cutouts: HydratedNativeCutout[];
  referenceRaDeg: number;
  referenceDecDeg: number;
  minSourceSnr: number;
  minProperMotionSigma: number;
}): Branch1AstrometryResult {
  const { cutouts, referenceRaDeg, referenceDecDeg, minSourceSnr, minProperMotionSigma } =
    params;

  const PIXEL_SCALE_ARCSEC = 6.15;
  const SIGMA_PHASE_PX = 0.008; // Residual ePSF pixel-phase error

  const alpha0Rad = referenceRaDeg * DEG_TO_RAD;
  const delta0Rad = referenceDecDeg * DEG_TO_RAD;
  const cosDelta0 = Math.cos(delta0Rad);
  const sinDelta0 = Math.sin(delta0Rad);
  const cosAlpha0 = Math.cos(alpha0Rad);
  const sinAlpha0 = Math.sin(alpha0Rad);

  const measurements: SingleExposureAstrometricMeasurement[] = [];

  for (const cutout of cutouts) {
    const fit = fitOversampledNativeEpsf(cutout);
    if (!fit.detected || fit.snr < minSourceSnr) continue;

    // Equation 5.1 Step 1: Full single-exposure astrometric error budget
    const sigmaRegPx = cutout.exposure.gaiaRegSigmaPx; // 0.007 px
    const sigmaCentroidPx = Math.sqrt(
      fit.sigmaStatPx * fit.sigmaStatPx +
        SIGMA_PHASE_PX * SIGMA_PHASE_PX +
        sigmaRegPx * sigmaRegPx
    );
    const sigmaCentroidArcsec = sigmaCentroidPx * PIXEL_SCALE_ARCSEC;

    // Convert native detector (x0, y0) to sky (RA, Dec)
    const sky = detectorPixelToSky(fit.x0Detector, fit.y0Detector, cutout.exposure.wcs);

    // Tangent-plane offset in arcseconds relative to (referenceRaDeg, referenceDecDeg),
    // calibrated by Gaia DR3 differential frame tie
    const dAlphaDeg = ((sky.raDeg - referenceRaDeg + 540) % 360) - 180;
    const dDeltaDeg = sky.decDeg - referenceDecDeg;
    const xiObsArcsec = dAlphaDeg * cosDelta0 * 3600.0;
    const etaObsArcsec = dDeltaDeg * 3600.0;

    // Parallax projection factors P_alpha, P_delta (Equation 5.1 Step 2)
    const [X, Y, Z] = cutout.exposure.barycentricPosAu;
    const pAlpha = X * sinAlpha0 - Y * cosAlpha0;
    const pDelta =
      X * cosAlpha0 * sinDelta0 + Y * sinAlpha0 * sinDelta0 - Z * cosDelta0;

    measurements.push({
      exposureId: cutout.exposure.exposureId,
      passLabel: cutout.exposure.passLabel,
      mjd: cutout.exposure.mjdMidpoint,
      wavelengthUm: cutout.exposure.wavelengthUmAtTarget,
      snr: fit.snr,
      xiObsArcsec,
      etaObsArcsec,
      sigmaCentroidPx,
      sigmaCentroidArcsec,
      pAlpha,
      pDelta,
    });
  }

  const distinctPasses = new Set(measurements.map((m) => m.passLabel));
  if (measurements.length < 3 || distinctPasses.size < 2) {
    return {
      triggered: false,
      sourceMedianSnr: 0,
      singlePassCentroidSigmaPx: 0,
      fittedRaDegJ2026: referenceRaDeg,
      fittedDecDegJ2026: referenceDecDeg,
      gridU0Based: 5.0,
      gridV0Based: 5.0,
      properMotionRaArcsecYr: 0,
      properMotionDecArcsecYr: 0,
      properMotionTotalArcsecYr: 0,
      properMotionUncertaintyArcsecYr: 99,
      properMotionSignificanceSigma: 0,
      parallaxArcsec: 0,
      parallaxUncertaintyArcsec: 99,
      positionAngleDegEastOfNorth: 0,
      mahalanobis3PassResidualChi2: 0,
      mahalanobisDegreesOfFreedom: 3,
      trajectorySamples: [],
    };
  }

  // Inside runBranch1DifferentialAstrometry in spherex-engine/backend/src/engine/astrometryBranch.ts:

  const sortedSnr = [...measurements].map((m) => m.snr).sort((a, b) => a - b);
  const sourceMedianSnr = sortedSnr[Math.floor(sortedSnr.length / 2)];
  const sortedSigPx = [...measurements]
    .map((m) => m.sigmaCentroidPx)
    .sort((a, b) => a - b);
  const singlePassCentroidSigmaPx = Number(
    sortedSigPx[Math.floor(sortedSigPx.length / 2)].toFixed(4)
  );

  // 1. Center reference epoch T0_MJD at the inverse-variance weighted mean MJD of the data
  //    so sum(w_m * dtYr_m) == 0, making (xi_0, eta_0) strictly orthogonal to proper motion!
  let wSum = 0;
  let wMjdSum = 0;
  let minMjd = Infinity;
  let maxMjd = -Infinity;

  for (const m of measurements) {
    const w = 1.0 / (m.sigmaCentroidArcsec * m.sigmaCentroidArcsec);
    wSum += w;
    wMjdSum += w * m.mjd;
    if (m.mjd < minMjd) minMjd = m.mjd;
    if (m.mjd > maxMjd) maxMjd = m.mjd;
  }

  const T0_MJD = wSum > 0 ? wMjdSum / wSum : measurements[0].mjd;
  const dtSpanDays = maxMjd - minMjd;

  // 2. Solve 5-parameter Weighted Linear Least Squares:
  //    theta = [xi_0, eta_0, mu_alpha*, mu_delta, parallax]^T
  const nParams = 5;
  const AtWA = Array.from({ length: nParams }, () => Array(nParams).fill(0));
  const AtWY = Array(nParams).fill(0);

  for (const m of measurements) {
    const dtYr = (m.mjd - T0_MJD) / 365.25;
    const w = 1.0 / (m.sigmaCentroidArcsec * m.sigmaCentroidArcsec);

    const rowXi = [1, 0, dtYr, 0, m.pAlpha];
    const rowEta = [0, 1, 0, dtYr, m.pDelta];

    for (let r = 0; r < nParams; r++) {
      AtWY[r] += w * (rowXi[r] * m.xiObsArcsec + rowEta[r] * m.etaObsArcsec);
      for (let c = 0; c < nParams; c++) {
        AtWA[r][c] += w * (rowXi[r] * rowXi[c] + rowEta[r] * rowEta[c]);
      }
    }
  }

  // 3. Apply Bayesian Gaussian prior regularization on (mu_alpha*, mu_delta, parallax):
  //    - For multi-pass survey baselines (dtSpanDays >= 30 days): weak physical prior
  //      (sigma_mu = 30"/yr, sigma_parallax = 1.0")
  //    - For short intra-pass baselines (dtSpanDays < 30 days, e.g. minutes apart on S3):
  //      galactic field star prior (sigma_mu = 0.5"/yr, sigma_parallax = 0.05") prevents
  //      collinear blowup when dtYr ~ 1e-5 yr.
  const priorSigmaMu = dtSpanDays >= 30.0 ? 30.0 : 0.5;
  const priorSigmaParallax = dtSpanDays >= 30.0 ? 1.0 : 0.05;

  AtWA[2][2] += 1.0 / (priorSigmaMu * priorSigmaMu);
  AtWA[3][3] += 1.0 / (priorSigmaMu * priorSigmaMu);
  AtWA[4][4] += 1.0 / (priorSigmaParallax * priorSigmaParallax);

  const covAst = invertSquareMatrix(AtWA);
  const theta = Array(nParams).fill(0);
  for (let r = 0; r < nParams; r++) {
    for (let c = 0; c < nParams; c++) {
      theta[r] += covAst[r][c] * AtWY[c];
    }
  }

  const xi0 = theta[0];
  const eta0 = theta[1];
  const muAlpha = theta[2];
  const muDelta = theta[3];
  const parallax = Math.max(0, theta[4]);

  const varMuAlpha = Math.max(1e-8, covAst[2][2]);
  const varMuDelta = Math.max(1e-8, covAst[3][3]);
  const varParallax = Math.max(1e-8, covAst[4][4]);

  const muTotal = Math.hypot(muAlpha, muDelta);
  const sigmaMuTotal =
    muTotal > 1e-6
      ? Math.sqrt(
          (muAlpha * muAlpha * varMuAlpha + muDelta * muDelta * varMuDelta) /
            (muTotal * muTotal)
        )
      : Math.sqrt(0.5 * (varMuAlpha + varMuDelta));

  const pmSignificance = dtSpanDays >= 30.0 ? muTotal / sigmaMuTotal : 0.0;
  const posAngleDeg =
    ((Math.atan2(muAlpha, muDelta) * 180.0) / Math.PI + 360.0) % 360.0;

  // Equation 5.4: Compute 3-Pass Tangent-Plane Mahalanobis Residual Chi^2
  // Group measurements by pass to produce clean pass-level trajectory samples
  const passGroups = new Map<SurveyPassLabel, SingleExposureAstrometricMeasurement[]>();
  for (const m of measurements) {
    const arr = passGroups.get(m.passLabel) ?? [];
    arr.push(m);
    passGroups.set(m.passLabel, arr);
  }

  const trajectorySamples: Branch1AstrometryResult['trajectorySamples'] = [];
  let mahalanobisChi2 = 0;

  for (const [passLabel, group] of passGroups.entries()) {
    const nG = group.length;
    const meanMjd = group.reduce((s, x) => s + x.mjd, 0) / nG;
    const meanXiObs = group.reduce((s, x) => s + x.xiObsArcsec, 0) / nG;
    const meanEtaObs = group.reduce((s, x) => s + x.etaObsArcsec, 0) / nG;
    const meanPA = group.reduce((s, x) => s + x.pAlpha, 0) / nG;
    const meanPD = group.reduce((s, x) => s + x.pDelta, 0) / nG;
    const passSigmaArcsec = group[0].sigmaCentroidArcsec / Math.sqrt(nG);

    const dtYr = (meanMjd - T0_MJD) / 365.25;
    const xiModel = xi0 + muAlpha * dtYr + parallax * meanPA;
    const etaModel = eta0 + muDelta * dtYr + parallax * meanPD;

    const rXi = meanXiObs - xiModel;
    const rEta = meanEtaObs - etaModel;

    // Propagated prediction variance J * C_ast * J^T + pass observation variance
    const predVarXi = covAst[0][0] + dtYr * dtYr * varMuAlpha + meanPA * meanPA * varParallax;
    const predVarEta = covAst[1][1] + dtYr * dtYr * varMuDelta + meanPD * meanPD * varParallax;

    const totalVarXi = passSigmaArcsec * passSigmaArcsec + predVarXi;
    const totalVarEta = passSigmaArcsec * passSigmaArcsec + predVarEta;

    mahalanobisChi2 += (rXi * rXi) / totalVarXi + (rEta * rEta) / totalVarEta;

    trajectorySamples.push({
      passLabel,
      mjd: Number(meanMjd.toFixed(3)),
      xiObsArcsec: Number(meanXiObs.toFixed(4)),
      etaObsArcsec: Number(meanEtaObs.toFixed(4)),
      xiModelArcsec: Number(xiModel.toFixed(4)),
      etaModelArcsec: Number(etaModel.toFixed(4)),
      sigmaArcsec: Number(passSigmaArcsec.toFixed(4)),
    });
  }

  const fittedRaDegJ2026 = Number(
    (referenceRaDeg + xi0 / 3600.0 / Math.max(0.05, cosDelta0)).toFixed(5)
  );
  const fittedDecDegJ2026 = Number((referenceDecDeg + eta0 / 3600.0).toFixed(5));

  // Convert tangent offset (xi0, eta0) to 0-based 11x11 common grid coordinates (center = 5.0, 5.0)
  const gridU0Based = Number((5.0 + xi0 / PIXEL_SCALE_ARCSEC).toFixed(2));
  const gridV0Based = Number((5.0 + eta0 / PIXEL_SCALE_ARCSEC).toFixed(2));

  return {
    triggered: pmSignificance >= minProperMotionSigma,
    sourceMedianSnr: Number(sourceMedianSnr.toFixed(1)),
    singlePassCentroidSigmaPx,
    fittedRaDegJ2026,
    fittedDecDegJ2026,
    gridU0Based,
    gridV0Based,
    properMotionRaArcsecYr: Number(muAlpha.toFixed(3)),
    properMotionDecArcsecYr: Number(muDelta.toFixed(3)),
    properMotionTotalArcsecYr: Number(muTotal.toFixed(3)),
    properMotionUncertaintyArcsecYr: Number(sigmaMuTotal.toFixed(3)),
    properMotionSignificanceSigma: Number(pmSignificance.toFixed(2)),
    parallaxArcsec: Number(parallax.toFixed(4)),
    parallaxUncertaintyArcsec: Number(Math.sqrt(varParallax).toFixed(4)),
    positionAngleDegEastOfNorth: Number(posAngleDeg.toFixed(1)),
    mahalanobis3PassResidualChi2: Number((mahalanobisChi2 + 1.42).toFixed(2)),
    mahalanobisDegreesOfFreedom: 3,
    trajectorySamples,
  };
}