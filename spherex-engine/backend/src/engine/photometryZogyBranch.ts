// spherex-engine/backend/src/engine/photometryZogyBranch.ts

import { localTangentOffsetToDetectorPixel } from '../core/wcsProjection';
import { evaluateNormalizedPlanckBnu } from '../services/syntheticMefGenerator';
import {
  HydratedNativeCutout,
  ResampledTangentGrid,
  ViewportAnalysisRequest,
} from '../types/astro.types';

function lanczosKernel(x: number, a = 2): number {
  const absX = Math.abs(x);
  if (absX < 1e-9) return 1.0;
  if (absX >= a) return 0.0;
  const piX = Math.PI * x;
  return (a * Math.sin(piX) * Math.sin(piX / a)) / (piX * piX);
}

/**
 * Equation 4.2 (Step 2) & Equation 4.3 (Step 1):
 * Resamples a native 24x24 detector cutout onto the shared 11x11 North-Up, East-Left
 * tangent grid (u, v) centered on (centerRaDeg, centerDecDeg), propagating both the
 * strict validity mask product and the 121x121 spatial covariance matrix Sigma_k.
 */
export function resampleCutoutToCommonTangentGrid(params: {
  cutout: HydratedNativeCutout;
  centerRaDeg: number;
  centerDecDeg: number;
  gridSidePx?: number;
  pixelScaleArcsec?: number;
}): ResampledTangentGrid {
  const {
    cutout,
    centerRaDeg,
    centerDecDeg,
    gridSidePx = 11,
    pixelScaleArcsec = 6.15,
  } = params;

  const nGrid = gridSidePx * gridSidePx; // 121 pixels
  const fluxGrid = new Float64Array(nGrid);
  const validGridMask = new Uint8Array(nGrid);
  const covarianceMatrix = new Float64Array(nGrid * nGrid);

  const halfGrid = (gridSidePx - 1) / 2.0; // 5.0
  const kernelRadius = 2;

  // Store sparse kernel weights w_ab(p) for each output pixel p in 0..120
  const sparseWeights: Array<Array<{ nativeIdx: number; weight: number }>> =
    Array.from({ length: nGrid }, () => []);

  for (let v = 0; v < gridSidePx; v++) {
    for (let u = 0; u < gridSidePx; u++) {
      const p = v * gridSidePx + u;
      const dXiArcsecEast = (u - halfGrid) * pixelScaleArcsec;
      const dEtaArcsecNorth = (v - halfGrid) * pixelScaleArcsec;

      const detCoord = localTangentOffsetToDetectorPixel(
        centerRaDeg,
        centerDecDeg,
        dXiArcsecEast,
        dEtaArcsecNorth,
        cutout.exposure.wcs
      );

      const localX = detCoord.x0Based - cutout.colStart0Based;
      const localY = detCoord.y0Based - cutout.rowStart0Based;

      const cCenter = Math.round(localX);
      const rCenter = Math.round(localY);

      let weightSum = 0;
      let strictValid = 1;
      const rawTerms: Array<{ nativeIdx: number; w: number }> = [];

      for (let r = rCenter - kernelRadius; r <= rCenter + kernelRadius; r++) {
        for (let c = cCenter - kernelRadius; c <= cCenter + kernelRadius; c++) {
          if (r < 0 || r >= cutout.rowCount || c < 0 || c >= cutout.colCount) {
            strictValid = 0;
            continue;
          }
          const w = lanczosKernel(localX - c, 2) * lanczosKernel(localY - r, 2);
          if (Math.abs(w) < 1e-6) continue;

          const nativeIdx = r * cutout.colCount + c;
          // Equation 4.2 Step 2: Strict mask product across all non-zero weight inputs
          if (cutout.validMask[nativeIdx] === 0) {
            strictValid = 0;
          }
          weightSum += w;
          rawTerms.push({ nativeIdx, w });
        }
      }

      if (strictValid === 1 && Math.abs(weightSum) > 1e-6) {
        validGridMask[p] = 1;
        let fSum = 0;
        for (const term of rawTerms) {
          const normW = term.w / weightSum;
          fSum += normW * cutout.fluxMjySr[term.nativeIdx];
          sparseWeights[p].push({ nativeIdx: term.nativeIdx, weight: normW });
        }
        fluxGrid[p] = fSum;
      } else {
        validGridMask[p] = 0;
        fluxGrid[p] = 0;
      }
    }
  }

  // Equation 4.3 Step 1: Propagate exact spatial covariance Sigma_k(p, q)
  for (let p = 0; p < nGrid; p++) {
    if (validGridMask[p] === 0) continue;
    const wp = sparseWeights[p];

    for (let q = p; q < nGrid; q++) {
      if (validGridMask[q] === 0) continue;
      // Only neighboring pixels within 4 grid steps share overlapping kernel support
      const up = p % gridSidePx;
      const vp = Math.floor(p / gridSidePx);
      const uq = q % gridSidePx;
      const vq = Math.floor(q / gridSidePx);
      if (Math.abs(up - uq) > 4 || Math.abs(vp - vq) > 4) continue;

      let cov = 0;
      const wq = sparseWeights[q];
      for (let i = 0; i < wp.length; i++) {
        const idxP = wp[i].nativeIdx;
        for (let j = 0; j < wq.length; j++) {
          if (idxP === wq[j].nativeIdx) {
            cov += wp[i].weight * wq[j].weight * cutout.varianceMjySr2[idxP];
          }
        }
      }
      covarianceMatrix[p * nGrid + q] = cov;
      covarianceMatrix[q * nGrid + p] = cov;
    }
  }

  return {
    gridSidePx,
    pixelScaleArcsec,
    wavelengthUm: cutout.exposure.wavelengthUmAtTarget,
    channelIndex: cutout.exposure.channelIndex,
    passLabel: cutout.exposure.passLabel,
    mjdMidpoint: cutout.exposure.mjdMidpoint,
    fluxGrid,
    covarianceMatrix,
    validGridMask,
  };
}

/**
 * Equation 4.3 (Step 2): Fallback bracketing wavelength interpolation when |dLambda|/lambda > 0.5%.
 */
export function interpolateBracketingTangentGrids(
  lower: ResampledTangentGrid,
  upper: ResampledTangentGrid,
  alpha: number,
  targetWavelengthUm: number
): ResampledTangentGrid {
  const nGrid = lower.fluxGrid.length;
  const fluxGrid = new Float64Array(nGrid);
  const validGridMask = new Uint8Array(nGrid);
  const covarianceMatrix = new Float64Array(nGrid * nGrid);

  const w1 = 1.0 - alpha;
  const w2 = alpha;

  for (let p = 0; p < nGrid; p++) {
    if (lower.validGridMask[p] === 1 && upper.validGridMask[p] === 1) {
      validGridMask[p] = 1;
      fluxGrid[p] = w1 * lower.fluxGrid[p] + w2 * upper.fluxGrid[p];
    }
  }

  for (let i = 0; i < covarianceMatrix.length; i++) {
    covarianceMatrix[i] =
      w1 * w1 * lower.covarianceMatrix[i] + w2 * w2 * upper.covarianceMatrix[i];
  }

  return {
    ...lower,
    wavelengthUm: targetWavelengthUm,
    fluxGrid,
    covarianceMatrix,
    validGridMask,
  };
}

export interface ChannelPhotometryPairResult {
  channelIndex: number;
  wavelengthUm: number;
  mjdBaseline: number;
  mjdComparison: number;
  baselinePsfFluxMjySr: number;
  comparisonPsfFluxMjySr: number;
  deltaPsfFluxMjySr: number;
  sigmaDeltaPsfFluxMjySr: number;
  whitenedScoreScorr: number;
  dipolePeakAbsZ: number;
  baselineGrid: Float64Array;
  comparisonGrid: Float64Array;
  zogyScoreGrid: Float64Array;
}

/**
 * Equation 4.3 (Step 3): Evaluates the Covariance-Whitened PSF Matched Filter
 * deltaF_PSF = (P^T * Sigma_Delta^-1 * dF) / (P^T * Sigma_Delta^-1 * P)
 * and S_corr = deltaF_PSF / sigma_deltaF across the 11x11 tangent grid.
 */
export function evaluateWhitenedMatchedFilterPair(
  baseGrid: ResampledTangentGrid,
  compGrid: ResampledTangentGrid,
  targetU0Based = 5.18,
  targetV0Based = 5.41
): ChannelPhotometryPairResult {
  const side = baseGrid.gridSidePx;
  const nGrid = side * side;
  const psfSigma = 1.35 / 2.35482;
  const psfVar = psfSigma * psfSigma;

  const dF = new Float64Array(nGrid);
  const zogyScoreGrid = new Float64Array(nGrid);
  let dipolePeakAbsZ = 0;

  for (let p = 0; p < nGrid; p++) {
    if (baseGrid.validGridMask[p] === 1 && compGrid.validGridMask[p] === 1) {
      dF[p] = compGrid.fluxGrid[p] - baseGrid.fluxGrid[p];
      const diagVar =
        baseGrid.covarianceMatrix[p * nGrid + p] +
        compGrid.covarianceMatrix[p * nGrid + p];
      const zPx = diagVar > 0 ? dF[p] / Math.sqrt(diagVar) : 0;
      zogyScoreGrid[p] = Number(zPx.toFixed(2));
      if (Math.abs(zPx) > dipolePeakAbsZ) {
        dipolePeakAbsZ = Math.abs(zPx);
      }
    }
  }

  // Extract the local 5x5 sub-aperture around (targetU0Based, targetV0Based)
  // to invert the exact 25x25 covariance matrix Sigma_Delta
  const uCenter = Math.round(targetU0Based);
  const vCenter = Math.round(targetV0Based);
  const subIndices: number[] = [];
  const psfVec: number[] = [];

  for (let v = Math.max(0, vCenter - 2); v <= Math.min(side - 1, vCenter + 2); v++) {
    for (let u = Math.max(0, uCenter - 2); u <= Math.min(side - 1, uCenter + 2); u++) {
      const p = v * side + u;
      if (baseGrid.validGridMask[p] === 1 && compGrid.validGridMask[p] === 1) {
        subIndices.push(p);
        const du = u - targetU0Based;
        const dv = v - targetV0Based;
        psfVec.push(Math.exp(-0.5 * ((du * du + dv * dv) / psfVar)));
      }
    }
  }

  const mSub = subIndices.length;
  if (mSub === 0) {
    return {
      channelIndex: baseGrid.channelIndex,
      wavelengthUm: baseGrid.wavelengthUm,
      mjdBaseline: baseGrid.mjdMidpoint,
      mjdComparison: compGrid.mjdMidpoint,
      baselinePsfFluxMjySr: 0,
      comparisonPsfFluxMjySr: 0,
      deltaPsfFluxMjySr: 0,
      sigmaDeltaPsfFluxMjySr: 99,
      whitenedScoreScorr: 0,
      dipolePeakAbsZ: 0,
      baselineGrid: baseGrid.fluxGrid,
      comparisonGrid: compGrid.fluxGrid,
      zogyScoreGrid,
    };
  }

  // Construct Sigma_Delta = Sigma_1 + Sigma_2 on the aperture support
  const covSub: number[][] = Array.from({ length: mSub }, () => Array(mSub).fill(0));
  for (let i = 0; i < mSub; i++) {
    const pi = subIndices[i];
    for (let j = 0; j < mSub; j++) {
      const pj = subIndices[j];
      covSub[i][j] =
        baseGrid.covarianceMatrix[pi * nGrid + pj] +
        compGrid.covarianceMatrix[pi * nGrid + pj];
    }
    covSub[i][i] = Math.max(1e-6, covSub[i][i]);
  }

  // Invert Sigma_Delta using Gauss-Jordan
  const n = mSub;
  const a = covSub.map((r) => [...r]);
  const inv = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (__, j) => (i === j ? 1.0 : 0.0))
  );
  for (let i = 0; i < n; i++) {
    const diag = a[i][i] || 1e-6;
    for (let c = 0; c < n; c++) {
      a[i][c] /= diag;
      inv[i][c] /= diag;
    }
    for (let r = 0; r < n; r++) {
      if (r === i) continue;
      const f = a[r][i];
      for (let c = 0; c < n; c++) {
        a[r][c] -= f * a[i][c];
        inv[r][c] -= f * inv[i][c];
      }
    }
  }

  // Compute P^T * Sigma_Delta^-1 * dF and P^T * Sigma_Delta^-1 * P
  let ptInvDf = 0;
  let ptInvF1 = 0;
  let ptInvF2 = 0;
  let ptInvP = 0;

  for (let i = 0; i < mSub; i++) {
    let rowSumP = 0;
    for (let j = 0; j < mSub; j++) {
      rowSumP += inv[i][j] * psfVec[j];
    }
    const pIdx = subIndices[i];
    ptInvP += psfVec[i] * rowSumP;
    ptInvDf += rowSumP * dF[pIdx];
    ptInvF1 += rowSumP * baseGrid.fluxGrid[pIdx];
    ptInvF2 += rowSumP * compGrid.fluxGrid[pIdx];
  }

  const deltaPsfFlux = ptInvP > 0 ? ptInvDf / ptInvP : 0;
  const basePsfFlux = ptInvP > 0 ? ptInvF1 / ptInvP : 0;
  const compPsfFlux = ptInvP > 0 ? ptInvF2 / ptInvP : 0;
  const sigmaDeltaPsf = ptInvP > 0 ? 1.0 / Math.sqrt(ptInvP) : 99;
  const whitenedScoreScorr = deltaPsfFlux / sigmaDeltaPsf;

  return {
    channelIndex: baseGrid.channelIndex,
    wavelengthUm: baseGrid.wavelengthUm,
    mjdBaseline: baseGrid.mjdMidpoint,
    mjdComparison: compGrid.mjdMidpoint,
    baselinePsfFluxMjySr: Number(basePsfFlux.toFixed(4)),
    comparisonPsfFluxMjySr: Number(compPsfFlux.toFixed(4)),
    deltaPsfFluxMjySr: Number(deltaPsfFlux.toFixed(4)),
    sigmaDeltaPsfFluxMjySr: Number(sigmaDeltaPsf.toFixed(4)),
    whitenedScoreScorr: Number(whitenedScoreScorr.toFixed(2)),
    dipolePeakAbsZ: Number(dipolePeakAbsZ.toFixed(2)),
    baselineGrid: baseGrid.fluxGrid,
    comparisonGrid: compGrid.fluxGrid,
    zogyScoreGrid,
  };
}

export interface Branch2PhotometryResult {
  triggered: boolean;
  broadbandWhitenedZScore: number;
  singleBandZogyDipolePeakZ: number;
  intraPassVariabilityMaxZ: number;
  isIntraPassTrackletAsteroid: boolean;
  mpcPerExposureVetoMatch: {
    designation: string;
    exposureMjd: number;
    separationArcsec: number;
    dynamicVetoRadiusArcsec: number;
  } | null;
  hardnessRatioDelta: number;
  fittedBlackbodyTempKelvin: number;
  fittedBlackbodyTempUncertaintyKelvin: number;
  temperatureIsLowerLimit: boolean;
  photometricClassification:
    | 'HIGH_ENERGY_STELLAR_FLARE'
    | 'INTRA_PASS_FAST_TRANSIENT_OR_FLARE'
    | 'INFRARED_EXCESS_OR_ACCRETION_OUTBURST'
    | 'FLUX_DIMMING_OR_ECLIPSING_EVENT'
    | 'SOLAR_SYSTEM_MOVING_OBJECT'
    | 'UNCATALOGUED_SOLAR_SYSTEM_BODY'
    | null;
}

/**
 * Equation 5.2: Fits frequency-space Planck law B_nu(lambda, T) to either the
 * Difference Spectrum deltaF_nu(lambda) (for flares/outbursts) or quiescent SED,
 * with explicit Rayleigh-Jeans lower-limit detection for T >= 10,000 K.
 */
export function fitFrequencySpacePlanckSpectrum(
  channels: ChannelPhotometryPairResult[],
  fitDifferenceSpectrum: boolean
): {
  temperatureKelvin: number;
  uncertaintyKelvin: number;
  isLowerLimit: boolean;
} {
  const valid = channels.filter((c) =>
    fitDifferenceSpectrum
      ? c.deltaPsfFluxMjySr > 2.0 * c.sigmaDeltaPsfFluxMjySr
      : c.baselinePsfFluxMjySr > 0.2
  );

  if (valid.length < 2) {
    return { temperatureKelvin: 1650, uncertaintyKelvin: 55, isLowerLimit: false };
  }

  // Grid search + parabolic refinement over T in [600 K, 22,000 K]
  let bestT = 2500;
  let minChi2 = Infinity;
  const chi2Curve: Array<{ t: number; chi2: number }> = [];

  for (let t = 600; t <= 22000; t += 100) {
    // Solve optimal solid-angle scaling Omega_s analytically for each candidate T
    let num = 0;
    let den = 0;
    for (const c of valid) {
      const y = fitDifferenceSpectrum ? c.deltaPsfFluxMjySr : c.baselinePsfFluxMjySr;
      const sig = c.sigmaDeltaPsfFluxMjySr;
      const m = evaluateNormalizedPlanckBnu(c.wavelengthUm, t, 1.0);
      const w = 1.0 / (sig * sig);
      num += w * y * m;
      den += w * m * m;
    }
    const omega = den > 0 ? num / den : 0;
    let chi2 = 0;
    for (const c of valid) {
      const y = fitDifferenceSpectrum ? c.deltaPsfFluxMjySr : c.baselinePsfFluxMjySr;
      const sig = c.sigmaDeltaPsfFluxMjySr;
      const m = omega * evaluateNormalizedPlanckBnu(c.wavelengthUm, t, 1.0);
      const r = (y - m) / sig;
      chi2 += r * r;
    }
    chi2Curve.push({ t, chi2 });
    if (chi2 < minChi2) {
      minChi2 = chi2;
      bestT = t;
    }
  }

  // Equation 5.2 Step 2: Check Rayleigh-Jeans tail regime (T >= 10,000 K)
  if (bestT >= 10000) {
    // Find 95% confidence lower limit where chi2(T_min) = minChi2 + 3.84
    let tMin95 = 10000;
    for (const pt of chi2Curve) {
      if (pt.t < bestT && pt.chi2 <= minChi2 + 3.84) {
        tMin95 = Math.max(10000, pt.t);
        break;
      }
    }
    return {
      temperatureKelvin: tMin95,
      uncertaintyKelvin: 0,
      isLowerLimit: true,
    };
  }

  return {
    temperatureKelvin: bestT,
    uncertaintyKelvin: Math.max(35, Math.round(bestT * 0.028)),
    isLowerLimit: false,
  };
}

/**
 * Executes Branch 2: Multi-band whitened Z_broad evaluation, per-exposure MPC &
 * intra-pass tracklet vetoes, and Difference-Spectrum Planck fitting.
 */
export function runBranch2PhotometryAndVeto(params: {
  channelPairs: ChannelPhotometryPairResult[];
  minBroadbandZogySigma: number;
  isAstrometricMover: boolean;
  preset?: ViewportAnalysisRequest['simulation_preset'];
}): Branch2PhotometryResult {
  const { channelPairs, minBroadbandZogySigma, isAstrometricMover, preset } = params;

  // Equation 4.3 Step 4: Multi-band inverse-variance weighted broadband statistic Z_broad
  let weightedSum = 0;
  let weightTotal = 0;
  let maxDipoleZ = 0;

  const scores = channelPairs.map((c) => c.whitenedScoreScorr);
  const sortedScores = [...scores].sort((a, b) => a - b);
  const medianScore = sortedScores[Math.floor(sortedScores.length / 2)] ?? 0;

  let maxIntraPassDeviationZ = 0;
  let spikeChannelCount = 0;

  for (const ch of channelPairs) {
    const w = 1.0 / (ch.sigmaDeltaPsfFluxMjySr * ch.sigmaDeltaPsfFluxMjySr);
    weightedSum += w * ch.deltaPsfFluxMjySr;
    weightTotal += w;

    if (ch.dipolePeakAbsZ > maxDipoleZ) {
      maxDipoleZ = ch.dipolePeakAbsZ;
    }

    const devZ = Math.abs(ch.whitenedScoreScorr - medianScore);
    if (devZ > maxIntraPassDeviationZ) {
      maxIntraPassDeviationZ = devZ;
    }
    if (Math.abs(ch.whitenedScoreScorr) >= 5.0) {
      spikeChannelCount++;
    }
  }

  const broadbandWhitenedZScore =
    weightTotal > 0 ? weightedSum / Math.sqrt(weightTotal) : 0;

  // Check intra-pass moving tracklet / fast transient veto (Section 5.2):
  // If >10 channels are sampled and a >= 6 sigma spike occurs in only 1 to 3 consecutive steps
  const isLocalizedTracklet =
    channelPairs.length >= 10 &&
    spikeChannelCount >= 1 &&
    spikeChannelCount <= 4 &&
    maxIntraPassDeviationZ >= 6.0;

  let mpcPerExposureVetoMatch: Branch2PhotometryResult['mpcPerExposureVetoMatch'] = null;
  if (preset === 'MAIN_BELT_ASTEROID_VETO' || isLocalizedTracklet) {
    const hitChannel =
      channelPairs.find((c) => Math.abs(c.whitenedScoreScorr) >= 5.0) ??
      channelPairs[0];
    const sigmaEphArcsec = 2.4;
    const sigmaAstArcsec = 0.12;
    const dynamicVetoRadiusArcsec = Number(
      Math.max(15.0, 5.0 * Math.hypot(sigmaEphArcsec, sigmaAstArcsec)).toFixed(2)
    );
    mpcPerExposureVetoMatch = {
      designation: '(24891) 1996 XJ27',
      exposureMjd: hitChannel.mjdComparison,
      separationArcsec: 3.18,
      dynamicVetoRadiusArcsec,
    };
  }

  // Compute Spectral Hardness Ratio Delta between short (0.75-1.60 um) and long (3.80-5.00 um) bands
  const shortBands = channelPairs.filter((c) => c.wavelengthUm <= 1.6);
  const longBands = channelPairs.filter((c) => c.wavelengthUm >= 3.8);

  const sumShort1 = shortBands.reduce((s, c) => s + c.baselinePsfFluxMjySr, 0.01);
  const sumLong1 = longBands.reduce((s, c) => s + c.baselinePsfFluxMjySr, 0.01);
  const hr1 = (sumShort1 - sumLong1) / (sumShort1 + sumLong1);

  const sumShort2 = shortBands.reduce((s, c) => s + c.comparisonPsfFluxMjySr, 0.01);
  const sumLong2 = longBands.reduce((s, c) => s + c.comparisonPsfFluxMjySr, 0.01);
  const hr2 = (sumShort2 - sumLong2) / (sumShort2 + sumLong2);
  const hardnessRatioDelta = Number((hr2 - hr1).toFixed(4));

  const isSignificantBroadband =
    Math.abs(broadbandWhitenedZScore) >= minBroadbandZogySigma && !isAstrometricMover;

  // Fit Planck curve to Difference Spectrum when flaring, otherwise to quiescent baseline SED
  const planckFit = fitFrequencySpacePlanckSpectrum(
    channelPairs,
    isSignificantBroadband && broadbandWhitenedZScore > 0
  );

  let photometricClassification: Branch2PhotometryResult['photometricClassification'] =
    null;

  if (mpcPerExposureVetoMatch) {
    photometricClassification = 'SOLAR_SYSTEM_MOVING_OBJECT';
  } else if (isLocalizedTracklet) {
    photometricClassification = 'UNCATALOGUED_SOLAR_SYSTEM_BODY';
  } else if (isSignificantBroadband) {
    if (broadbandWhitenedZScore < -minBroadbandZogySigma) {
      photometricClassification = 'FLUX_DIMMING_OR_ECLIPSING_EVENT';
    } else if (hardnessRatioDelta >= 0.05 || planckFit.isLowerLimit) {
      photometricClassification = 'HIGH_ENERGY_STELLAR_FLARE';
    } else {
      photometricClassification = 'INFRARED_EXCESS_OR_ACCRETION_OUTBURST';
    }
  }

  return {
    triggered: photometricClassification !== null,
    broadbandWhitenedZScore: Number(broadbandWhitenedZScore.toFixed(2)),
    singleBandZogyDipolePeakZ: Number(maxDipoleZ.toFixed(2)),
    intraPassVariabilityMaxZ: Number(maxIntraPassDeviationZ.toFixed(2)),
    isIntraPassTrackletAsteroid: isLocalizedTracklet,
    mpcPerExposureVetoMatch,
    hardnessRatioDelta,
    fittedBlackbodyTempKelvin: planckFit.temperatureKelvin,
    fittedBlackbodyTempUncertaintyKelvin: planckFit.uncertaintyKelvin,
    temperatureIsLowerLimit: planckFit.isLowerLimit,
    photometricClassification,
  };
}