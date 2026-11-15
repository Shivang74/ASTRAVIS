// spherex-engine/backend/src/server.ts

import express, { Request, Response } from 'express';
import cors from 'cors';
import { z } from 'zod';
import {
  AnomalyClassification,
  HydratedNativeCutout,
  QUALITY_FLAG_BITS,
  SURVEY_PASS_REGISTRY,
  ViewportAnalysisRequest,
} from './types/astro.types';
import { querySpatialSpectralSidecarIndex } from './services/sidecarIndex';
import { fetchAndHydrateExposureCutout } from './services/s3RangeFetcher';
import { runBranch1DifferentialAstrometry } from './engine/astrometryBranch';
import {
  ChannelPhotometryPairResult,
  evaluateWhitenedMatchedFilterPair,
  interpolateBracketingTangentGrids,
  resampleCutoutToCommonTangentGrid,
  runBranch2PhotometryAndVeto,
} from './engine/photometryZogyBranch';
import { runBranch3SpectralIceAnalysis } from './engine/spectralIceBranch';
import { runLiveS3FullThreeBranchPipeline } from './services/nasaLiveApis';

const app = express();
app.use(cors());
app.use(express.json());

const RequestSchema = z.object({
  viewport: z.object({
    ra_deg: z.number().min(0).max(360),
    dec_deg: z.number().min(-90).max(90),
    cutout_side_arcsec: z.number().min(15).max(180).default(67.65),
    geometry: z.literal('SQUARE_TANGENT_GRID').default('SQUARE_TANGENT_GRID'),
  }),
  execution_mode: z
    .enum(['ISO_WAVELENGTH_SCREENING', 'FULL_SED_102_BAND'])
    .default('FULL_SED_102_BAND'),
  passes: z.object({
    baseline: z.enum(['PASS_1', 'PASS_2', 'PASS_3', 'PASS_4']).default('PASS_1'),
    comparison: z.enum(['PASS_1', 'PASS_2', 'PASS_3', 'PASS_4']).default('PASS_2'),
    verification: z
      .enum(['PASS_1', 'PASS_2', 'PASS_3', 'PASS_4'])
      .optional()
      .default('PASS_3'),
  }),
  thresholds: z
    .object({
      min_astrometric_source_snr: z.number().default(10.0),
      min_proper_motion_snr_sigma: z.number().default(5.0),
      min_broadband_zogy_sigma: z.number().default(5.0),
      min_ice_ew_snr_sigma: z.number().default(4.5),
      bad_pixel_bitmask: z
        .number()
        .default(QUALITY_FLAG_BITS.DEFAULT_REJECT_MASK),
    })
    .default({}),
  simulation_preset: z
    .enum([
      'BROWN_DWARF_PROPER_MOTION',
      'HOT_STELLAR_FLARE',
      'WATER_ICE_OUTGASSING',
      'MAIN_BELT_ASTEROID_VETO',
    ])
    .optional()
    .default('BROWN_DWARF_PROPER_MOTION'),
});

app.get('/api/v1/health', (_req: Request, res: Response) => {
  res.json({
    service: 'SPEC-ASTRO-MEF-03 Anomaly Detection Engine',
    revision: '3.0.0',
    status: 'HEALTHY',
    passes_registered: SURVEY_PASS_REGISTRY,
  });
});

app.get('/api/v1/live-irsa-search', async (req: Request, res: Response) => {
  const detectorOverride = Number(req.query.detector ?? 1);
  try {
    const result = await runLiveS3FullThreeBranchPipeline({
      detectorOverride,
    });
    res.status(200).json(result);
  } catch (err) {
    res.status(502).json({
      status: 'LIVE_S3_PIPELINE_ERROR',
      message: (err as Error).message,
    });
  }
});

app.post('/api/v1/spectral-delta', async (req: Request, res: Response) => {
  const startTimeMs = Date.now();

  const parsed = RequestSchema.safeParse(req.body);
  if (!parsed.success) {
    res.status(400).json({
      status: 'INVALID_REQUEST_PAYLOAD',
      errors: parsed.error.flatten(),
    });
    return;
  }

  const payload = parsed.data as ViewportAnalysisRequest;
  const { ra_deg, dec_deg } = payload.viewport;

  // Section 7.4: Check bright-star saturation exclusion zone (e.g., Vega at RA 279.2347, Dec +38.7837)
  const distToVegaArcsec =
    Math.hypot(
      (ra_deg - 279.2347) * Math.cos((dec_deg * Math.PI) / 180.0),
      dec_deg - 38.7837
    ) * 3600.0;

  if (distToVegaArcsec < 120.0) {
    res.status(422).json({
      status: 'BRIGHT_STAR_SATURATION_VETO',
      message:
        'Requested coordinates fall inside the bright-star saturation exclusion zone of Alpha Lyrae (Vega), where >50% of detector pixels fail Bit 1 (0x0002 OVERFLOW_SATURATED).',
      saturated_fraction: 0.78,
    });
    return;
  }

  try {
    // Step 1: Query Spatial-Spectral Sidecar Index (WCS & WCS-WAVE pre-evaluated)
    const { workItems, allExposuresToFetch } = querySpatialSpectralSidecarIndex({
      raDeg: ra_deg,
      decDeg: dec_deg,
      executionMode: payload.execution_mode,
      passes: payload.passes,
    });

    // Step 2: Concurrent 24-row byte-range fetches across all exposures
    const cutoutMap = new Map<string, HydratedNativeCutout>();
    let totalBytesFetchedS3 = 0;
    let failedExposureReads = 0;
    const warnings: string[] = [];

    const fetchPromises = allExposuresToFetch.map(async (exp) => {
      try {
        const { cutout, bytesFetched } = await fetchAndHydrateExposureCutout({
          exposure: exp,
          raDeg: ra_deg,
          decDeg: dec_deg,
          rejectBitmask: payload.thresholds.bad_pixel_bitmask,
          preset: payload.simulation_preset ?? 'BROWN_DWARF_PROPER_MOTION',
        });
        cutoutMap.set(exp.exposureId, cutout);
        totalBytesFetchedS3 += bytesFetched;
      } catch (err) {
        failedExposureReads++;
        warnings.push(
          `Failed reading exposure ${exp.exposureId} (ch=${exp.channelIndex}): ${
            (err as Error).message
          }`
        );
      }
    });

    await Promise.all(fetchPromises);

    const totalExposuresAttempted = allExposuresToFetch.length;
    const failureRatio = failedExposureReads / Math.max(1, totalExposuresAttempted);

    // Section 7.4: >15% exposure failure triggers HTTP 502
    if (failureRatio > 0.15) {
      res.status(502).json({
        status: 'INSUFFICIENT_SPECTRAL_COVERAGE',
        failed_exposure_reads: failedExposureReads,
        total_exposures_attempted: totalExposuresAttempted,
        warnings,
      });
      return;
    }

    const allHydratedCutouts = Array.from(cutoutMap.values());

    // Check if >50% of pixels are saturated across hydrated cutouts
    const meanSatFraction =
      allHydratedCutouts.reduce((s, c) => s + c.saturatedFraction, 0) /
      Math.max(1, allHydratedCutouts.length);
    if (meanSatFraction > 0.5) {
      res.status(422).json({
        status: 'BRIGHT_STAR_SATURATION_VETO',
        saturated_fraction: Number(meanSatFraction.toFixed(3)),
      });
      return;
    }

    // Step 3: Execute Branch 1 — Differential ePSF Astrometry & 5-Parameter Fit
    const branch1 = runBranch1DifferentialAstrometry({
      cutouts: allHydratedCutouts,
      referenceRaDeg: ra_deg,
      referenceDecDeg: dec_deg,
      minSourceSnr: payload.thresholds.min_astrometric_source_snr,
      minProperMotionSigma: payload.thresholds.min_proper_motion_snr_sigma,
    });

    // Step 4: Resample onto Common 11x11 Tangent Grid & Run Branch 2 (Whitened ZOGY Photometry)
    const channelPairs: ChannelPhotometryPairResult[] = [];
    for (const item of workItems) {
      const baseCutout = cutoutMap.get(item.baselineExposure.exposureId);
      const compPrimaryCutout = cutoutMap.get(
        item.comparisonPrimaryExposure.exposureId
      );
      if (!baseCutout || !compPrimaryCutout) continue;

      const baseGrid = resampleCutoutToCommonTangentGrid({
        cutout: baseCutout,
        centerRaDeg: ra_deg,
        centerDecDeg: dec_deg,
      });

      let compGrid = resampleCutoutToCommonTangentGrid({
        cutout: compPrimaryCutout,
        centerRaDeg: ra_deg,
        centerDecDeg: dec_deg,
      });

      if (
        item.comparisonBracketUpper &&
        item.interpolationAlpha !== undefined
      ) {
        const upperCutout = cutoutMap.get(
          item.comparisonBracketUpper.exposureId
        );
        if (upperCutout) {
          const upperGrid = resampleCutoutToCommonTangentGrid({
            cutout: upperCutout,
            centerRaDeg: ra_deg,
            centerDecDeg: dec_deg,
          });
          compGrid = interpolateBracketingTangentGrids(
            compGrid,
            upperGrid,
            item.interpolationAlpha,
            item.referenceWavelengthUm
          );
        }
      }

      const pairEval = evaluateWhitenedMatchedFilterPair(
        baseGrid,
        compGrid,
        branch1.gridU0Based,
        branch1.gridV0Based
      );
      channelPairs.push(pairEval);
    }

    const branch2 = runBranch2PhotometryAndVeto({
      channelPairs,
      minBroadbandZogySigma: payload.thresholds.min_broadband_zogy_sigma,
      isAstrometricMover: branch1.triggered,
      preset: payload.simulation_preset,
    });

    // Step 5: Execute Branch 3 — Solid-State Volatile Ice vs. Cool-Star Gas SED Analysis
    const branch3 = runBranch3SpectralIceAnalysis({
      channelPairs,
      minIceEwSnrSigma: payload.thresholds.min_ice_ew_snr_sigma,
    });

    // Determine final triggered branches and primary classification
    const triggeredBranches: string[] = [];
    if (branch1.triggered) triggeredBranches.push('BRANCH_1_EPSF_ASTROMETRY');
    if (branch2.triggered) triggeredBranches.push('BRANCH_2_WHITENED_PHOTOMETRY');
    if (branch3.triggered) triggeredBranches.push('BRANCH_3_VOLATILE_ICE_SED');

    let primaryClassification: AnomalyClassification =
      'UNCLASSIFIED_SIGNIFICANT_RESIDUAL';

    if (
      branch2.photometricClassification === 'SOLAR_SYSTEM_MOVING_OBJECT' ||
      branch2.photometricClassification === 'UNCATALOGUED_SOLAR_SYSTEM_BODY'
    ) {
      primaryClassification = branch2.photometricClassification;
    } else if (branch3.triggered) {
      primaryClassification = 'VOLATILE_ICE_ABSORPTION_VARIATION';
    } else if (branch1.triggered) {
      primaryClassification = 'ASTROMETRIC_PROPER_MOTION';
    } else if (branch2.photometricClassification) {
      primaryClassification = branch2.photometricClassification;
    }

    const executionTimeMs = Date.now() - startTimeMs;
    if (executionTimeMs > 2500) {
      res.status(504).json({
        status: 'GATEWAY_SLA_TIMEOUT',
        execution_time_ms: executionTimeMs,
      });
      return;
    }

    // Representative channel cutout for visual heatmap inspection (near 1.62 um or middle channel)
    const repPair =
      channelPairs[Math.floor(channelPairs.length / 2)] ?? channelPairs[0];

    const passesEvaluatedCount = payload.passes.verification ? 3 : 2;

    res.status(200).json({
      status: failedExposureReads > 0 ? 'PARTIAL' : 'COMPLETE',
      metadata: {
        tangent_grid_dimensions_px: [11, 11],
        native_row_span_fetched_px: 24,
        pixel_indexing_convention: '0_BASED_TANGENT_GRID',
        exposures_queried_per_pass: workItems.length,
        passes_evaluated: passesEvaluatedCount,
        total_exposures_fetched: allHydratedCutouts.length,
        bytes_fetched_per_exposure: 783360,
        bytes_fetched_total_s3: totalBytesFetchedS3,
        wavelength_channels_assembled: channelPairs.length,
        execution_time_ms: executionTimeMs,
        failed_exposure_reads: failedExposureReads,
      },
      anomalies: [
        {
          id: `ANOM_J${branch1.fittedRaDegJ2026.toFixed(4)}+${branch1.fittedDecDegJ2026.toFixed(
            4
          )}_${payload.passes.baseline}${payload.passes.comparison}_001`,
          triggered_branches: triggeredBranches,
          passes_evaluated: {
            PASS_1_median_mjd: 60842.318,
            PASS_2_median_mjd: 61025.812,
            PASS_3_median_mjd: 61208.114,
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
            mpc_per_exposure_veto_match: branch2.mpcPerExposureVetoMatch,
            mahalanobis_3pass_residual_chi2:
              branch1.mahalanobis3PassResidualChi2,
            mahalanobis_degrees_of_freedom:
              branch1.mahalanobisDegreesOfFreedom,
          },
          photometry_and_sed: {
            broadband_whitened_z_score: branch2.broadbandWhitenedZScore,
            intra_pass_variability_max_z: branch2.intraPassVariabilityMaxZ,
            sampled_wavelength_range_um: [0.75, 5.0],
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
        representative_channel_um: repPair?.wavelengthUm ?? 1.62,
        baseline_grid_11x11: repPair ? Array.from(repPair.baselineGrid) : [],
        comparison_grid_11x11: repPair ? Array.from(repPair.comparisonGrid) : [],
        zogy_score_grid_11x11: repPair ? Array.from(repPair.zogyScoreGrid) : [],
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
      warnings,
    });
  } catch (err) {
    res.status(500).json({
      status: 'INTERNAL_ENGINE_ERROR',
      message: (err as Error).message,
    });
  }
});

const PORT = Number(process.env.PORT ?? 4000);
app.listen(PORT, () => {
  console.log(
    `[SPEC-ASTRO-MEF-03] Anomaly Detection Engine listening on port ${PORT}`
  );
});