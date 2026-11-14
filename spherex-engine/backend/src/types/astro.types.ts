// spherex-engine/backend/src/types/astro.types.ts

export type SurveyPassLabel = 'PASS_1' | 'PASS_2' | 'PASS_3' | 'PASS_4';
export type ExecutionMode = 'ISO_WAVELENGTH_SCREENING' | 'FULL_SED_102_BAND';

export type AnomalyClassification =
  | 'ASTROMETRIC_PROPER_MOTION'
  | 'HIGH_ENERGY_STELLAR_FLARE'
  | 'INTRA_PASS_FAST_TRANSIENT_OR_FLARE'
  | 'INFRARED_EXCESS_OR_ACCRETION_OUTBURST'
  | 'FLUX_DIMMING_OR_ECLIPSING_EVENT'
  | 'VOLATILE_ICE_ABSORPTION_VARIATION'
  | 'SOLAR_SYSTEM_MOVING_OBJECT'
  | 'UNCATALOGUED_SOLAR_SYSTEM_BODY'
  | 'ARTIFACT_PERSISTENCE'
  | 'ARTIFACT_GHOST_STREAK'
  | 'UNCLASSIFIED_SIGNIFICANT_RESIDUAL';

export interface SurveyPassMetadata {
  epochId: 'E1' | 'E2' | 'E3' | 'E4';
  label: SurveyPassLabel;
  mjdStart: number; // Inclusive
  mjdEnd: number;   // Exclusive [mjdStart, mjdEnd)
  calendarWindow: string;
  status: 'COMPLETE' | 'IN_PROGRESS';
  coverageFraction: number;
}

export const SURVEY_PASS_REGISTRY: Record<SurveyPassLabel, SurveyPassMetadata> = {
  PASS_1: {
    epochId: 'E1',
    label: 'PASS_1',
    mjdStart: 60796.0,
    mjdEnd: 60980.0,
    calendarWindow: '2025-05-01 to 2025-11-01',
    status: 'COMPLETE',
    coverageFraction: 1.0,
  },
  PASS_2: {
    epochId: 'E2',
    label: 'PASS_2',
    mjdStart: 60980.0,
    mjdEnd: 61161.0,
    calendarWindow: '2025-11-01 to 2026-05-01',
    status: 'COMPLETE',
    coverageFraction: 1.0,
  },
  PASS_3: {
    epochId: 'E3',
    label: 'PASS_3',
    mjdStart: 61161.0,
    mjdEnd: 61345.0,
    calendarWindow: '2026-05-01 to 2026-11-01',
    status: 'COMPLETE',
    coverageFraction: 1.0,
  },
  PASS_4: {
    epochId: 'E4',
    label: 'PASS_4',
    mjdStart: 61345.0,
    mjdEnd: 61526.0,
    calendarWindow: '2026-11-01 to 2027-05-01',
    status: 'IN_PROGRESS',
    coverageFraction: 0.072,
  },
};

/**
 * Bitwise Quality Flags (HDU 2: FLAGS)
 * Default rejection mask is 0xFFFB (65531), which keeps Bit 2 (NONLINEAR_RESPONSE)
 * unmasked and inflates its variance by 2.25x (1.5x sigma).
 */
export const QUALITY_FLAG_BITS = {
  TRANSIENT_OR_COSMIC_RAY: 0x0001,
  OVERFLOW_SATURATED: 0x0002,
  NONLINEAR_RESPONSE: 0x0004,
  DEAD_OR_HOT_PIXEL: 0x0008,
  PERSISTENCE_ARTIFACT: 0x0010,
  OPTICAL_GHOST_OR_STREAK: 0x0020,
  TELEMETRY_OR_FLAT_ERROR: 0xffc0,
  DEFAULT_REJECT_MASK: 0xfffb, // 65531
} as const;

export interface FitsStructuralHeader {
  xtension: string;
  extname: string;
  bitpix: number;
  naxis: number;
  naxes: number[];
  pcount: number;
  gcount: number;
  bscale: number;
  bzero: number;
  dataStartByte: number;
  nextHduStartByte: number;
  records: Record<string, string | number | boolean>;
}

export interface SipPolynomialCoeffs {
  aOrder: number;
  bOrder: number;
  apOrder?: number;
  bpOrder?: number;
  a: Record<string, number>;  // e.g. "2_0": 1.2e-7
  b: Record<string, number>;
  ap: Record<string, number>;
  bp: Record<string, number>;
}

export interface CelestialWCS {
  crpix1: number; // 1-based FITS pixel center (nominally 1020.5)
  crpix2: number;
  crval1: number; // RA (deg)
  crval2: number; // Dec (deg)
  cd1_1: number;
  cd1_2: number;
  cd2_1: number;
  cd2_2: number;
  rollAngleRad: number;
  sip?: SipPolynomialCoeffs;
}

export interface SidecarExposureEntry {
  exposureId: string;
  s3Uri: string;
  passLabel: SurveyPassLabel;
  mjdMidpoint: number;
  barycentricPosAu: [number, number, number]; // [X, Y, Z] in ICRS AU
  wcs: CelestialWCS;
  wavelengthUmAtTarget: number;
  bandwidthUmAtTarget: number;
  channelIndex: number; // 0..101
  hduDataOffsets: {
    sciStartByte: number;
    flagsStartByte: number;
    varianceStartByte: number;
    varianceIsSigma: boolean; // True if EXTNAME='ERR'
    zodiStartByte: number;
  };
  gaiaRegOffsetTangentArcsec: [number, number]; // [dXi, dEta] from Gaia DR3 tie
  gaiaRegSigmaPx: number;
}

export interface HydratedNativeCutout {
  exposure: SidecarExposureEntry;
  rowStart0Based: number; // j_min
  rowCount: number;       // 24 rows
  colStart0Based: number; // i_min (local trim)
  colCount: number;       // 24 cols
  targetDetectorX0Based: number;
  targetDetectorY0Based: number;
  fluxMjySr: Float32Array;     // ZODI-subtracted calibrated flux
  varianceMjySr2: Float32Array; // Bit-2 inflated variance
  flags: Uint32Array;
  validMask: Uint8Array;       // 1 = valid, 0 = rejected
  saturatedFraction: number;
  hasPersistenceFlag: boolean;
  hasGhostFlag: boolean;
}

export interface ResampledTangentGrid {
  gridSidePx: number;          // 11
  pixelScaleArcsec: number;    // 6.15
  wavelengthUm: number;
  channelIndex: number;
  passLabel: SurveyPassLabel;
  mjdMidpoint: number;
  fluxGrid: Float64Array;      // Length 121 (11x11)
  covarianceMatrix: Float64Array; // 121 x 121 symmetric spatial covariance
  validGridMask: Uint8Array;   // Strict product of contributing native pixels
}

export interface ViewportAnalysisRequest {
  viewport: {
    ra_deg: number;
    dec_deg: number;
    cutout_side_arcsec: number;
    geometry: 'SQUARE_TANGENT_GRID';
  };
  execution_mode: ExecutionMode;
  passes: {
    baseline: SurveyPassLabel;
    comparison: SurveyPassLabel;
    verification?: SurveyPassLabel;
  };
  thresholds: {
    min_astrometric_source_snr: number;
    min_proper_motion_snr_sigma: number;
    min_broadband_zogy_sigma: number;
    min_ice_ew_snr_sigma: number;
    bad_pixel_bitmask: number;
  };
  simulation_preset?: 'BROWN_DWARF_PROPER_MOTION' | 'HOT_STELLAR_FLARE' | 'WATER_ICE_OUTGASSING' | 'MAIN_BELT_ASTEROID_VETO';
}