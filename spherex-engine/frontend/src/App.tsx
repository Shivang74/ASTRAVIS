// spherex-engine/frontend/src/App.tsx
import React, { useState, useEffect, useCallback } from 'react';
import './App.css';

type PresetType =
  | 'BROWN_DWARF_PROPER_MOTION'
  | 'HOT_STELLAR_FLARE'
  | 'WATER_ICE_OUTGASSING'
  | 'MAIN_BELT_ASTEROID_VETO';

interface TrajectorySample {
  passLabel: string;
  mjd: number;
  xiObsArcsec: number;
  etaObsArcsec: number;
  xiModelArcsec: number;
  etaModelArcsec: number;
  sigmaArcsec: number;
}

interface SedChannel {
  channel_index: number;
  wavelength_um: number;
  baseline_flux_mjy_sr: number;
  comparison_flux_mjy_sr: number;
  delta_flux_mjy_sr: number;
  sigma_delta_mjy_sr: number;
  whitened_score_z: number;
}

interface IceEvaluation {
  species: string;
  bandCenterUm: number;
  baselineEwUm: number;
  comparisonEwUm: number;
  deltaEwUm: number;
  sigmaDeltaEwUm: number;
  significanceSigma: number;
  passesBhFdr: boolean;
}

interface ApiResponse {
  status: string;
  message?: string;
  saturated_fraction?: number;
  discovery_mode_used?: string;
  irsa_tap_proxy_status?: string;
  metadata?: {
    tangent_grid_dimensions_px: [number, number];
    native_row_span_fetched_px: number;
    exposures_queried_per_pass: number;
    passes_evaluated: number;
    total_exposures_fetched: number;
    bytes_fetched_per_exposure: number;
    bytes_fetched_total_s3: number;
    wavelength_channels_assembled: number;
    execution_time_ms: number;
    failed_exposure_reads: number;
  };
  anomalies?: Array<{
    id: string;
    triggered_branches: string[];
    passes_evaluated?: {
      PASS_1_median_mjd: number;
      PASS_2_median_mjd: number;
      PASS_3_median_mjd: number;
    };
    coordinates_epoch_j2026: {
      ra_deg: number;
      dec_deg: number;
      grid_u_0based: number;
      grid_v_0based: number;
    };
    astrometry_5param_fit: {
      source_median_snr: number;
      single_pass_centroid_sigma_px: number;
      proper_motion_ra_arcsec_yr: number;
      proper_motion_dec_arcsec_yr: number;
      proper_motion_total_arcsec_yr: number;
      proper_motion_uncertainty_arcsec_yr: number;
      proper_motion_significance_sigma: number;
      parallax_arcsec: number;
      parallax_uncertainty_arcsec: number;
      position_angle_deg_east_of_north: number;
      single_band_zogy_dipole_peak_z: number;
      mpc_per_exposure_veto_match: {
        designation: string;
        exposureMjd: number;
        separationArcsec: number;
        dynamicVetoRadiusArcsec: number;
      } | null;
      mahalanobis_3pass_residual_chi2: number;
      mahalanobis_degrees_of_freedom: number;
    };
    photometry_and_sed: {
      broadband_whitened_z_score: number;
      intra_pass_variability_max_z: number;
      fitted_blackbody_temp_kelvin: number;
      fitted_blackbody_temp_uncertainty_kelvin: number;
      temperature_is_lower_limit: boolean;
      hardness_ratio_delta: number;
      h2o_305_delta_ew_um: number;
      h2o_305_delta_ew_uncertainty_um: number;
      h2o_305_significance_sigma: number;
      co2_427_delta_ew_um: number;
      co_467_delta_ew_um: number;
      delta_bic_ice_vs_gas: number;
    };
    classification: string;
  }>;
  visualization_payload?: {
    representative_channel_um: number;
    baseline_grid_11x11: number[];
    comparison_grid_11x11: number[];
    zogy_score_grid_11x11: number[];
    astrometric_trajectory: TrajectorySample[];
    sed_channels: SedChannel[];
    ice_band_evaluations: IceEvaluation[];
  };
  live_s3_telemetry?: {
    total_live_s3_exposures_found: number;
    live_s3_mef_hdu_walk: Array<{
      hduIndex: number;
      extname: string;
      bitpix: number;
      naxes: number[];
      headerStartByte: number;
      dataStartByte: number;
      nextHduStartByte: number;
    }>;
  };
}

function fluxToColor(val: number, min: number, max: number, diverging = false): string {
  if (diverging) {
    const absMax = Math.max(Math.abs(min), Math.abs(max), 3.0);
    const norm = Math.max(-1, Math.min(1, val / absMax));
    if (norm < 0) {
      const t = Math.abs(norm);
      return `rgb(${Math.round(15 * (1 - t) + 56 * t)}, ${Math.round(
        23 * (1 - t) + 189 * t
      )}, ${Math.round(42 * (1 - t) + 248 * t)})`;
    }
    const t = norm;
    return `rgb(${Math.round(15 * (1 - t) + 244 * t)}, ${Math.round(
      23 * (1 - t) + 63 * t
    )}, ${Math.round(42 * (1 - t) + 94 * t)})`;
  }
  const t = Math.max(0, Math.min(1, (val - min) / Math.max(1e-5, max - min)));
  const r = Math.round(10 + 245 * Math.pow(t, 1.3));
  const g = Math.round(18 + 190 * Math.pow(t, 0.95));
  const b = Math.round(40 + 215 * Math.pow(t, 0.6));
  return `rgb(${r}, ${g}, ${b})`;
}

export function App() {
  const [preset, setPreset] = useState<PresetType>('BROWN_DWARF_PROPER_MOTION');
  const [raDeg, setRaDeg] = useState<string>('269.4521');
  const [decDeg, setDecDeg] = useState<string>('66.0234');
  const [execMode, setExecMode] = useState<'FULL_SED_102_BAND' | 'ISO_WAVELENGTH_SCREENING'>(
    'FULL_SED_102_BAND'
  );
  const [sedViewMode, setSedViewMode] = useState<'EPOCH_OVERLAY' | 'DIFFERENCE_SPECTRUM'>(
    'EPOCH_OVERLAY'
  );
  const [hoverPixel, setHoverPixel] = useState<{ u: number; v: number; val: number } | null>(
    null
  );
  const [loading, setLoading] = useState<boolean>(false);
  const [httpCode, setHttpCode] = useState<number>(200);
  const [data, setData] = useState<ApiResponse | null>(null);
  const [isLiveS3View, setIsLiveS3View] = useState<boolean>(false);

  const runAnalysis = useCallback(
    async (overridePreset?: PresetType, overrideRa?: number, overrideDec?: number) => {
      setLoading(true);
      setIsLiveS3View(false);
      const activePreset = overridePreset ?? preset;
      const activeRa = overrideRa ?? Number(raDeg);
      const activeDec = overrideDec ?? Number(decDeg);

      try {
        const res = await fetch('/api/v1/spectral-delta', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            viewport: {
              ra_deg: activeRa,
              dec_deg: activeDec,
              cutout_side_arcsec: 67.65,
              geometry: 'SQUARE_TANGENT_GRID',
            },
            execution_mode: execMode,
            passes: {
              baseline: 'PASS_1',
              comparison: 'PASS_2',
              verification: 'PASS_3',
            },
            simulation_preset: activePreset,
          }),
        });
        setHttpCode(res.status);
        const json: ApiResponse = await res.json();
        setData(json);
      } catch (err) {
        setHttpCode(500);
        setData({ status: 'NETWORK_ERROR', message: (err as Error).message });
      } finally {
        setLoading(false);
      }
    },
    [preset, raDeg, decDeg, execMode]
  );

  const runLiveS3Stream = async () => {
    setLoading(true);
    setIsLiveS3View(true);
    try {
      const res = await fetch('/api/v1/live-irsa-search');
      setHttpCode(res.status);
      const json: ApiResponse = await res.json();
      setData(json);
      const liveCoords = json.anomalies?.[0]?.coordinates_epoch_j2026;
      if (liveCoords) {
        setRaDeg(liveCoords.ra_deg.toFixed(4));
        setDecDeg(liveCoords.dec_deg.toFixed(4));
      }
    } catch (err) {
      setHttpCode(502);
      setData({ status: 'LIVE_S3_ERROR', message: (err as Error).message });
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    runAnalysis('BROWN_DWARF_PROPER_MOTION', 269.4521, 66.0234);
  }, []);

  const handleSelectPreset = (newPreset: PresetType) => {
    setPreset(newPreset);
    setRaDeg('269.4521');
    setDecDeg('66.0234');
    if (newPreset === 'HOT_STELLAR_FLARE') {
      setSedViewMode('DIFFERENCE_SPECTRUM');
    } else {
      setSedViewMode('EPOCH_OVERLAY');
    }
    runAnalysis(newPreset, 269.4521, 66.0234);
  };

  const handleTestVegaVeto = () => {
    setRaDeg('279.2415');
    setDecDeg('38.7900');
    runAnalysis(preset, 279.2415, 38.79);
  };

  const anom = data?.anomalies?.[0];
  const viz = data?.visualization_payload;
  const meta = data?.metadata;
  const liveTelemetry = data?.live_s3_telemetry;

  const render11x11Grid = (
    title: string,
    gridData: number[],
    diverging = false,
    showVector = false
  ) => {
    if (!gridData || gridData.length !== 121) return null;
    const minVal = Math.min(...gridData);
    const maxVal = Math.max(...gridData);
    const cellPx = 22;

    const u0 = anom?.coordinates_epoch_j2026.grid_u_0based ?? 5.0;
    const v0 = anom?.coordinates_epoch_j2026.grid_v_0based ?? 5.0;
    const paRad =
      ((anom?.astrometry_5param_fit.position_angle_deg_east_of_north ?? 0) * Math.PI) /
      180.0;
    const vecLenPx = 34;
    const cx = (u0 + 0.5) * cellPx;
    const cy = (10 - v0 + 0.5) * cellPx;
    const ex = cx + Math.sin(paRad) * vecLenPx;
    const ey = cy - Math.cos(paRad) * vecLenPx;

    return (
      <div className="heatmap-box">
        <h3>{title}</h3>
        <svg width={11 * cellPx} height={11 * cellPx} style={{ borderRadius: 4 }}>
          {gridData.map((val, idx) => {
            const u = idx % 11;
            const v = Math.floor(idx / 11);
            const screenY = (10 - v) * cellPx;
            return (
              <rect
                key={idx}
                x={u * cellPx}
                y={screenY}
                width={cellPx - 1}
                height={cellPx - 1}
                fill={fluxToColor(val, minVal, maxVal, diverging)}
                onMouseEnter={() => setHoverPixel({ u, v, val })}
              />
            );
          })}
          <circle
            cx={cx}
            cy={cy}
            r={4.5}
            fill="none"
            stroke="#38bdf8"
            strokeWidth={1.8}
          />
          {showVector && (
            <line
              x1={cx}
              y1={cy}
              x2={ex}
              y2={ey}
              stroke="#10b981"
              strokeWidth={2.2}
            />
          )}
        </svg>
        <div
          style={{
            fontSize: '0.7rem',
            fontFamily: 'var(--font-mono)',
            color: 'var(--text-muted)',
            marginTop: 6,
          }}
        >
          Range: [{minVal.toFixed(3)}, {maxVal.toFixed(3)}]
        </div>
      </div>
    );
  };

  return (
    <div className="console-wrapper">
      {/* HEADER */}
      <header className="top-header">
        <div className="title-block">
          <h1>SPEC-ASTRO-MEF-03 // 4D Spectral & Astrometric Anomaly Engine</h1>
          <p>
            24-Row Big-Endian S3 Byte-Range Slicing • 5-Param (μ, ϖ) ePSF Fit • ZOGY
            Covariance-Whitened Photometry • 102-Band SED
          </p>
        </div>
        <div className="badge-row">
          <span
            className={`badge ${
              httpCode === 200 ? 'ok' : httpCode === 422 ? 'warn' : 'err'
            }`}
          >
            HTTP {httpCode} // {data?.status ?? 'IDLE'}
          </span>
          {anom && <span className="badge ok">CLASS: {anom.classification}</span>}
          {isLiveS3View && data?.discovery_mode_used && (
            <span className="badge ok">SOURCE: {data.discovery_mode_used}</span>
          )}
        </div>
      </header>

      {/* CONTROLS */}
      <section className="control-bar">
        <div className="preset-group">
          <span className="group-label">
            Synthetic Injection Validation Presets (1–5) &amp; Live AWS S3 SPHEREx DR1
            Full 3-Branch Pipeline (6)
          </span>
          <div className="preset-buttons">
            <button
              className={`btn-preset ${
                !isLiveS3View &&
                preset === 'BROWN_DWARF_PROPER_MOTION' &&
                httpCode === 200
                  ? 'active'
                  : ''
              }`}
              onClick={() => handleSelectPreset('BROWN_DWARF_PROPER_MOTION')}
            >
              1. L/T Brown Dwarf (μ = 0.77″/yr, ϖ = 0.048″)
            </button>
            <button
              className={`btn-preset ${
                !isLiveS3View && preset === 'HOT_STELLAR_FLARE' && httpCode === 200
                  ? 'active'
                  : ''
              }`}
              onClick={() => handleSelectPreset('HOT_STELLAR_FLARE')}
            >
              2. Rayleigh-Jeans Stellar Flare (T &gt; 10⁴ K)
            </button>
            <button
              className={`btn-preset ${
                !isLiveS3View && preset === 'WATER_ICE_OUTGASSING' && httpCode === 200
                  ? 'active'
                  : ''
              }`}
              onClick={() => handleSelectPreset('WATER_ICE_OUTGASSING')}
            >
              3. Volatile H₂O / CO₂ Ice Absorption (3.05 μm)
            </button>
            <button
              className={`btn-preset ${
                !isLiveS3View && preset === 'MAIN_BELT_ASTEROID_VETO' && httpCode === 200
                  ? 'active'
                  : ''
              }`}
              onClick={() => handleSelectPreset('MAIN_BELT_ASTEROID_VETO')}
            >
              4. Main-Belt Asteroid Veto (MPC + Intra-Pass)
            </button>
            <button className="btn-preset danger" onClick={handleTestVegaVeto}>
              5. Vega Saturation Veto (HTTP 422)
            </button>
            <button
              className={`btn-preset ${isLiveS3View ? 'active' : ''}`}
              style={{
                borderColor: '#10b981',
                color: isLiveS3View ? '#34d399' : '#a7f3d0',
              }}
              onClick={runLiveS3Stream}
            >
              6. LIVE AWS S3 DR1 Full 3-Branch Pipeline (D1–D6)
            </button>
          </div>
        </div>

        <div className="coord-inputs">
          <div className="input-field">
            <span className="group-label">RA (deg)</span>
            <input
              type="number"
              step="0.001"
              value={raDeg}
              onChange={(e) => setRaDeg(e.target.value)}
            />
          </div>
          <div className="input-field">
            <span className="group-label">Dec (deg)</span>
            <input
              type="number"
              step="0.001"
              value={decDeg}
              onChange={(e) => setDecDeg(e.target.value)}
            />
          </div>
          <div className="input-field">
            <span className="group-label">Pipeline Mode</span>
            <select
              value={execMode}
              onChange={(e) =>
                setExecMode(
                  e.target.value as 'FULL_SED_102_BAND' | 'ISO_WAVELENGTH_SCREENING'
                )
              }
              style={{ width: 175 }}
            >
              <option value="FULL_SED_102_BAND">FULL_SED (All Bands)</option>
              <option value="ISO_WAVELENGTH_SCREENING">SCREENING (9 Files)</option>
            </select>
          </div>
        </div>

        <button
          className="btn-run"
          disabled={loading}
          onClick={() => (isLiveS3View ? runLiveS3Stream() : runAnalysis())}
        >
          {loading ? 'Slicing FITS...' : 'Execute Slice'}
        </button>
      </section>

      {/* HTTP 422 / ERROR BANNER */}
      {httpCode !== 200 && data && (
        <div className="error-banner">
          <h3 style={{ color: '#fda4af', marginBottom: 6 }}>
            Pipeline Veto Triggered: {data.status} (HTTP {httpCode})
          </h3>
          <p style={{ fontSize: '0.88rem' }}>{data.message}</p>
          {data.saturated_fraction !== undefined && (
            <p
              style={{
                fontFamily: 'var(--font-mono)',
                fontSize: '0.8rem',
                marginTop: 6,
                color: '#fbbf24',
              }}
            >
              Bit 1 (0x0002 OVERFLOW_SATURATED) Cutout Fraction:{' '}
              {(data.saturated_fraction * 100).toFixed(1)}% (&gt; 50.0% Threshold)
            </p>
          )}
        </div>
      )}

      {/* UNIFIED 6-CARD TELEMETRY STRIP (WORKS FOR BOTH PRESETS 1-4 AND LIVE S3 MODE 6) */}
      {meta && anom && httpCode === 200 && (
        <section className="telemetry-strip">
          <div className="stat-card">
            <div className="label">
              {isLiveS3View ? 'Live AWS S3 Bytes Sliced' : 'S3 Byte-Range Volume'}
            </div>
            <div className="value">
              {(meta.bytes_fetched_total_s3 / 1e6).toFixed(2)} MB
            </div>
            <div className="sub">
              {meta.total_exposures_fetched} FITS Files (
              {isLiveS3View ? 'Detectors D1–D6' : '783.36 KB/File'})
            </div>
          </div>
          <div className="stat-card">
            <div className="label">Rotated Row Span</div>
            <div className="value">{meta.native_row_span_fetched_px} Rows</div>
            <div className="sub">
              {meta.wavelength_channels_assembled} Spectral Channels Assembled
            </div>
          </div>
          <div className="stat-card">
            <div className="label">Branch 1: Proper Motion</div>
            <div className="value">
              {anom.astrometry_5param_fit.proper_motion_total_arcsec_yr.toFixed(3)}″/yr
            </div>
            <div className="sub">
              {anom.astrometry_5param_fit.proper_motion_significance_sigma.toFixed(1)}σ (ϖ ={' '}
              {anom.astrometry_5param_fit.parallax_arcsec.toFixed(3)}″)
            </div>
          </div>
          <div className="stat-card">
            <div className="label">Branch 2: Whitened Z_broad</div>
            <div className="value">
              {anom.photometry_and_sed.broadband_whitened_z_score.toFixed(2)}σ
            </div>
            <div className="sub">
              Single-Band Dipole Peak:{' '}
              {anom.astrometry_5param_fit.single_band_zogy_dipole_peak_z.toFixed(2)}σ
            </div>
          </div>
          <div className="stat-card">
            <div className="label">Planck B_ν(λ, T) Fit</div>
            <div className="value">
              {anom.photometry_and_sed.temperature_is_lower_limit ? '> ' : ''}
              {anom.photometry_and_sed.fitted_blackbody_temp_kelvin} K
            </div>
            <div className="sub">
              {anom.photometry_and_sed.temperature_is_lower_limit
                ? '95% CL Rayleigh-Jeans Limit'
                : `± ${anom.photometry_and_sed.fitted_blackbody_temp_uncertainty_kelvin} K (b_ν = 5099.4)`}
            </div>
          </div>
          <div className="stat-card">
            <div className="label">Compute Latency</div>
            <div className="value">{meta.execution_time_ms} ms</div>
            <div className="sub">
              {isLiveS3View
                ? 'Live nasa-irsa-spherex S3'
                : `Triggered: ${anom.triggered_branches.join(', ')}`}
            </div>
          </div>
        </section>
      )}

      {/* UNIFIED 3-BRANCH DASHBOARD PANELS */}
      {viz && anom && httpCode === 200 && (
        <>
          <div className="main-grid">
            <div className="panel">
              <div className="panel-header">
                <h2>
                  {isLiveS3View
                    ? `Live AWS S3 DR1 WCS-Aligned 11×11 Tangent Cutouts (λ = ${viz.representative_channel_um} μm)`
                    : `WCS-Aligned 11×11 Common Tangent-Grid Cutouts (λ = ${viz.representative_channel_um} μm)`}
                </h2>
                <span
                  style={{
                    fontSize: '0.75rem',
                    fontFamily: 'var(--font-mono)',
                    color: 'var(--accent-blue)',
                  }}
                >
                  {hoverPixel
                    ? `Pixel (u=${hoverPixel.u}, v=${hoverPixel.v}) : ${hoverPixel.val.toFixed(
                        3
                      )}`
                    : 'Hover pixel to inspect 0-based (u, v)'}
                </span>
              </div>
              <div className="heatmap-trio">
                {render11x11Grid(
                  isLiveS3View
                    ? 'S3 Exp 0001_1 (SCI − ZODI)'
                    : 'PASS_1 Baseline (MJy/sr)',
                  viz.baseline_grid_11x11,
                  false,
                  false
                )}
                {render11x11Grid(
                  isLiveS3View
                    ? 'S3 Exp 0001_2 (SCI − ZODI)'
                    : 'PASS_2 Comparison (MJy/sr)',
                  viz.comparison_grid_11x11,
                  false,
                  true
                )}
                {render11x11Grid(
                  'Whitened ZOGY Score (σ)',
                  viz.zogy_score_grid_11x11,
                  true,
                  true
                )}
              </div>

              {anom.astrometry_5param_fit.mpc_per_exposure_veto_match && (
                <div
                  style={{
                    marginTop: 12,
                    padding: '10px 12px',
                    background: 'rgba(245, 158, 11, 0.14)',
                    border: '1px solid var(--accent-amber)',
                    borderRadius: 6,
                    fontSize: '0.8rem',
                  }}
                >
                  <strong style={{ color: '#fbbf24' }}>
                    SOLAR SYSTEM CONTAMINANT VETO MATCH:
                  </strong>{' '}
                  Matched Minor Planet Center ephemeris{' '}
                  <code>
                    {anom.astrometry_5param_fit.mpc_per_exposure_veto_match.designation}
                  </code>{' '}
                  at exposure midpoint MJD{' '}
                  {anom.astrometry_5param_fit.mpc_per_exposure_veto_match.exposureMjd} (sep ={' '}
                  {
                    anom.astrometry_5param_fit.mpc_per_exposure_veto_match
                      .separationArcsec
                  }
                  ″ &lt; r_veto ={' '}
                  {
                    anom.astrometry_5param_fit.mpc_per_exposure_veto_match
                      .dynamicVetoRadiusArcsec
                  }
                  ″).
                </div>
              )}
            </div>

            <div className="panel">
              <div className="panel-header">
                <h2>
                  Branch 1: 5-Parameter Barycentric Astrometry (ξ₀, η₀, μ_α*, μ_δ, ϖ)
                </h2>
                <span className="badge ok">
                  Mahalanobis χ² ={' '}
                  {anom.astrometry_5param_fit.mahalanobis_3pass_residual_chi2} (df = 3)
                </span>
              </div>

              <div
                style={{
                  display: 'grid',
                  gridTemplateColumns: '220px 1fr',
                  gap: 14,
                  alignItems: 'center',
                }}
              >
                <svg
                  width={215}
                  height={215}
                  style={{
                    background: 'var(--bg-deep)',
                    border: '1px solid var(--border-subtle)',
                    borderRadius: 6,
                  }}
                >
                  <line
                    x1={107}
                    y1={0}
                    x2={107}
                    y2={215}
                    stroke="#1e293b"
                    strokeDasharray="3,3"
                  />
                  <line
                    x1={0}
                    y1={107}
                    x2={215}
                    y2={107}
                    stroke="#1e293b"
                    strokeDasharray="3,3"
                  />
                  {(() => {
                    const pts = [...viz.astrometric_trajectory].sort(
                      (a, b) => a.mjd - b.mjd
                    );
                    if (pts.length === 0) return null;
                    const meanXi =
                      pts.reduce((s, p) => s + p.xiObsArcsec, 0) / pts.length;
                    const meanEta =
                      pts.reduce((s, p) => s + p.etaObsArcsec, 0) / pts.length;
                    const scale = 115;

                    const toScreen = (xi: number, eta: number) => ({
                      x: 107 + (xi - meanXi) * scale,
                      y: 107 - (eta - meanEta) * scale,
                    });

                    const modelPath = pts
                      .map((p, idx) => {
                        const sc = toScreen(p.xiModelArcsec, p.etaModelArcsec);
                        return `${idx === 0 ? 'M' : 'L'} ${sc.x} ${sc.y}`;
                      })
                      .join(' ');

                    return (
                      <>
                        <path
                          d={modelPath}
                          fill="none"
                          stroke="#6366f1"
                          strokeWidth={2}
                          strokeDasharray="4,3"
                        />
                        {pts.map((p) => {
                          const obs = toScreen(p.xiObsArcsec, p.etaObsArcsec);
                          const errPx = Math.max(4, p.sigmaArcsec * scale);
                          return (
                            <g key={p.passLabel}>
                              <circle
                                cx={obs.x}
                                cy={obs.y}
                                r={errPx}
                                fill="rgba(56, 189, 248, 0.15)"
                                stroke="#38bdf8"
                                strokeWidth={1}
                              />
                              <circle cx={obs.x} cy={obs.y} r={3.5} fill="#10b981" />
                              <text
                                x={obs.x + 7}
                                y={obs.y - 5}
                                fill="#e2e8f0"
                                fontSize="9"
                                fontFamily="monospace"
                              >
                                {p.passLabel}
                              </text>
                            </g>
                          );
                        })}
                      </>
                    );
                  })()}
                </svg>

                <table className="metrics-table">
                  <tbody>
                    <tr>
                      <th>Native ePSF Median SNR</th>
                      <td className="mono">
                        {anom.astrometry_5param_fit.source_median_snr} (σ_c ={' '}
                        {anom.astrometry_5param_fit.single_pass_centroid_sigma_px} px)
                      </td>
                    </tr>
                    <tr>
                      <th>Proper Motion (μ_α*, μ_δ)</th>
                      <td className="mono">
                        ({anom.astrometry_5param_fit.proper_motion_ra_arcsec_yr}″/yr,{' '}
                        {anom.astrometry_5param_fit.proper_motion_dec_arcsec_yr}″/yr)
                      </td>
                    </tr>
                    <tr>
                      <th>Total Proper Motion (μ)</th>
                      <td className="mono" style={{ color: '#34d399' }}>
                        {anom.astrometry_5param_fit.proper_motion_total_arcsec_yr} ±{' '}
                        {anom.astrometry_5param_fit.proper_motion_uncertainty_arcsec_yr}{' '}
                        ″/yr ({anom.astrometry_5param_fit.proper_motion_significance_sigma}
                        σ)
                      </td>
                    </tr>
                    <tr>
                      <th>Annual Parallax (ϖ)</th>
                      <td className="mono">
                        {anom.astrometry_5param_fit.parallax_arcsec} ±{' '}
                        {anom.astrometry_5param_fit.parallax_uncertainty_arcsec} ″
                      </td>
                    </tr>
                    <tr>
                      <th>Sky Position Angle (θ_PA)</th>
                      <td className="mono">
                        {anom.astrometry_5param_fit.position_angle_deg_east_of_north}°
                        (East of North)
                      </td>
                    </tr>
                    <tr>
                      <th>Single-Band ZOGY Dipole</th>
                      <td className="mono" style={{ color: '#fbbf24' }}>
                        {anom.astrometry_5param_fit.single_band_zogy_dipole_peak_z}σ
                      </td>
                    </tr>
                  </tbody>
                </table>
              </div>
            </div>
          </div>

          {/* BOTTOM PANEL: SPECTRUM & ICE */}
          <div className="panel" style={{ marginBottom: 14 }}>
            <div className="panel-header">
              <h2>
                {isLiveS3View
                  ? 'Branch 2 & 3: Live AWS S3 Multi-Detector Reconstructed SED (Detectors D1–D6, 0.78 – 4.96 μm)'
                  : 'Branch 2 & 3: Reconstructed Near-Infrared Spectral Energy Distribution (0.75 – 5.00 μm) & Solid Ice Absorption'}
              </h2>
              <div className="preset-buttons">
                <button
                  className={`btn-preset ${
                    sedViewMode === 'EPOCH_OVERLAY' ? 'active' : ''
                  }`}
                  onClick={() => setSedViewMode('EPOCH_OVERLAY')}
                >
                  Epoch Overlay (Baseline vs Comparison)
                </button>
                <button
                  className={`btn-preset ${
                    sedViewMode === 'DIFFERENCE_SPECTRUM' ? 'active' : ''
                  }`}
                  onClick={() => setSedViewMode('DIFFERENCE_SPECTRUM')}
                >
                  Difference Spectrum ΔF_ν(λ)
                </button>
              </div>
            </div>

            <div
              style={{
                display: 'grid',
                gridTemplateColumns: '1fr 370px',
                gap: 16,
                alignItems: 'start',
              }}
            >
              <svg
                width="100%"
                height={250}
                viewBox="0 0 860 250"
                style={{
                  background: 'var(--bg-deep)',
                  border: '1px solid var(--border-subtle)',
                  borderRadius: 6,
                }}
              >
                {(() => {
                  const channels = viz.sed_channels;
                  if (channels.length === 0) return null;

                  const padL = 48;
                  const padR = 20;
                  const padT = 22;
                  const padB = 34;
                  const w = 860 - padL - padR;
                  const h = 250 - padT - padB;

                  const xMap = (lam: number) =>
                    padL + ((lam - 0.75) / (5.0 - 0.75)) * w;

                  const allY =
                    sedViewMode === 'EPOCH_OVERLAY'
                      ? channels.flatMap((c) => [
                          c.baseline_flux_mjy_sr,
                          c.comparison_flux_mjy_sr,
                        ])
                      : channels.map((c) => c.delta_flux_mjy_sr);

                  const minY = Math.min(-0.5, ...allY);
                  const maxY = Math.max(1.5, ...allY) * 1.15;
                  const yMap = (val: number) =>
                    padT + h - ((val - minY) / Math.max(1e-5, maxY - minY)) * h;

                  const basePts = channels
                    .map((c) => `${xMap(c.wavelength_um)},${yMap(c.baseline_flux_mjy_sr)}`)
                    .join(' ');
                  const compPts = channels
                    .map(
                      (c) => `${xMap(c.wavelength_um)},${yMap(c.comparison_flux_mjy_sr)}`
                    )
                    .join(' ');
                  const diffPts = channels
                    .map((c) => `${xMap(c.wavelength_um)},${yMap(c.delta_flux_mjy_sr)}`)
                    .join(' ');

                  return (
                    <>
                      {[
                        { label: 'H₂O Ice (3.05 μm)', l1: 2.75, l2: 3.5 },
                        { label: 'CO₂ (4.27)', l1: 4.18, l2: 4.36 },
                        { label: 'CO (4.67)', l1: 4.6, l2: 4.74 },
                      ].map((band) => (
                        <g key={band.label}>
                          <rect
                            x={xMap(band.l1)}
                            y={padT}
                            width={Math.max(4, xMap(band.l2) - xMap(band.l1))}
                            height={h}
                            fill="rgba(56, 189, 248, 0.08)"
                          />
                          <text
                            x={xMap(band.l1) + 3}
                            y={padT + 12}
                            fill="#38bdf8"
                            fontSize="9.5"
                            fontFamily="monospace"
                          >
                            {band.label}
                          </text>
                        </g>
                      ))}

                      <line
                        x1={padL}
                        y1={yMap(0)}
                        x2={860 - padR}
                        y2={yMap(0)}
                        stroke="#334155"
                        strokeDasharray="3,3"
                      />

                      {sedViewMode === 'EPOCH_OVERLAY' ? (
                        <>
                          <polyline
                            fill="none"
                            stroke="#64748b"
                            strokeWidth={1.6}
                            points={basePts}
                          />
                          <polyline
                            fill="none"
                            stroke="#38bdf8"
                            strokeWidth={2}
                            points={compPts}
                          />
                        </>
                      ) : (
                        <polyline
                          fill="none"
                          stroke="#f43f5e"
                          strokeWidth={2.2}
                          points={diffPts}
                        />
                      )}

                      {[1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0].map((lam) => (
                        <text
                          key={lam}
                          x={xMap(lam)}
                          y={250 - 10}
                          fill="#94a3b8"
                          fontSize="10"
                          textAnchor="middle"
                          fontFamily="monospace"
                        >
                          {lam.toFixed(1)} μm
                        </text>
                      ))}
                    </>
                  );
                })()}
              </svg>

              <div>
                <div className="group-label" style={{ marginBottom: 6 }}>
                  Branch 3: Solid-State Ice vs. Gas Template (ΔBIC ={' '}
                  {anom.photometry_and_sed.delta_bic_ice_vs_gas})
                </div>
                <table className="metrics-table">
                  <thead>
                    <tr>
                      <th>Band</th>
                      <th>ΔEW (μm)</th>
                      <th>SNR</th>
                      <th>BH-FDR</th>
                    </tr>
                  </thead>
                  <tbody>
                    {viz.ice_band_evaluations.map((ice) => (
                      <tr key={ice.species}>
                        <td className="mono">{ice.species}</td>
                        <td className="mono">
                          {ice.deltaEwUm.toFixed(3)} ± {ice.sigmaDeltaEwUm.toFixed(3)}
                        </td>
                        <td className="mono">{ice.significanceSigma.toFixed(2)}σ</td>
                        <td
                          className="mono"
                          style={{
                            color: ice.passesBhFdr ? '#34d399' : '#94a3b8',
                          }}
                        >
                          {ice.passesBhFdr ? 'PASS' : 'REJECT'}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
                <p
                  style={{
                    fontSize: '0.74rem',
                    color: 'var(--text-muted)',
                    marginTop: 10,
                  }}
                >
                  * Continuum anchors for H₂O (2.50–2.62 μm, 3.65–3.80 μm) explicitly
                  clear 2.7 μm H₂O vapor and 3.3 μm T-dwarf CH₄ bands.
                </p>
              </div>
            </div>
          </div>

          {/* EXTRA PANEL WHEN IN LIVE S3 MODE: HDU BYTE-BOUNDARY WALK TABLE */}
          {isLiveS3View && liveTelemetry && (
            <div className="panel">
              <div className="panel-header">
                <h2>
                  Live AWS S3 Multi-Extension FITS HDU Byte-Boundary Walk (Equation 2.1 on
                  nasa-irsa-spherex)
                </h2>
                <span className="badge ok">
                  {liveTelemetry.total_live_s3_exposures_found} Live DR1 Files Indexed •
                  13 Files Sliced Across D1–D6
                </span>
              </div>
              <table className="metrics-table">
                <thead>
                  <tr>
                    <th>HDU Index</th>
                    <th>EXTNAME</th>
                    <th>BITPIX</th>
                    <th>NAXIS Dimensions</th>
                    <th>Header Start Byte</th>
                    <th>Binary Data Start Byte</th>
                    <th>Next HDU Boundary Byte</th>
                  </tr>
                </thead>
                <tbody>
                  {liveTelemetry.live_s3_mef_hdu_walk.map((h) => (
                    <tr key={h.hduIndex}>
                      <td className="mono">HDU {h.hduIndex}</td>
                      <td className="mono" style={{ color: '#38bdf8' }}>
                        {h.extname}
                      </td>
                      <td className="mono">{h.bitpix}</td>
                      <td className="mono">
                        {h.naxes.length > 0 ? h.naxes.join(' × ') : '0 (Primary Metadata)'}
                      </td>
                      <td className="mono">{h.headerStartByte.toLocaleString()}</td>
                      <td className="mono">{h.dataStartByte.toLocaleString()}</td>
                      <td className="mono">{h.nextHduStartByte.toLocaleString()}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </>
      )}
    </div>
  );
}

export default App;