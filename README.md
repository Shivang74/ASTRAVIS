# ASTRAVIS: Surgical 4D Spectral & Astrometric Anomaly Engine for NASA SPHEREx

**Specification Revision:** `SPEC-ASTRO-MEF-03` (v3.0 Production Grade)

**Data Sources:** Live AWS Open Data S3 (`nasa-irsa-spherex`), Caltech/IPAC IRSA TAP/SIA, NASA/JPL Solar System Dynamics (`sb_ident.api`), ESA Gaia DR3

**Tech Stack:** TypeScript, Node.js, Express, React 18, Vite, Custom SVG 4D Telemetry Engine

---

## 1. Overview & The "Data Gravity" Challenge

NASA's **SPHEREx** (Spectro-Photometer for the History of the Universe, Epoch of Reionization, and Ices Explorer) space telescope maps the entire sky every 6 months across **102 near-infrared spectral channels (0.75 to 5.00 um)**. Because SPHEREx captures light through Linear Variable Filters (LVFs), a single star's complete 102-band Spectral Energy Distribution (SED) is spread across roughly 100 separate exposures per survey pass.

Each Level-2 Multi-Extension FITS (MEF) exposure is **63.34 MB** and contains 6 Header Data Units (`PRIMARY`, `IMAGE`, `FLAGS`, `VARIANCE`, `ZODI`, and `WCS-WAVE`). Under traditional workflows, inspecting a single target across three 6-month survey passes requires downloading nearly 300 full FITS files (over **18.5 GB** of raw data) just to analyze an `11 x 11` pixel patch of sky.

**ASTRAVIS** solves this bottleneck by bringing surgical byte-range streaming and a decoupled 3-branch astrophysics engine directly to the browser:

- **99% Bandwidth Reduction via HTTP 206 Byte-Range Slicing:** Instead of downloading full 63.34 MB files, I probe FITS headers using 28 KB range reads (Equation 3.1), step cleanly over tile-compressed binary tables (`HDU 2: FLAGS`), and fetch only the exact **24 contiguous detector rows** (`195.84 KB` per HDU) needed to cover the rotated target cutout.
- **Resilient Dual-Source Cloud Discovery:** I query Caltech/IPAC IRSA's Table Access Protocol (`/TAP/sync`) and automatically fail over to direct **AWS S3 `ListObjectsV2`** discovery against `nasa-irsa-spherex.s3.us-east-1.amazonaws.com` whenever upstream university gateways experience downtime.
- **Three Parallel Physics Branches:** I simultaneously evaluate native sub-pixel astrometry (separating proper motion from Earth's annual parallax), covariance-whitened ZOGY difference photometry (with live NASA JPL asteroid vetoes), and 102-band solid-state volatile ice spectroscopy (`H2O`, `CO2`, `CO`).

---

## 2. System Architecture (`SPEC-ASTRO-MEF-03`)

```text
+------------------------------------------------------------------------------+
|                          CLIENT VISUALIZATION LAYER                          |
|  1. Emits Request: { ra_deg, dec_deg, cutout_side_arcsec, mode, passes }     |
+--------------------------------------+---------------------------------------+
                                       | HTTP POST (Application/JSON)
                                       v
+------------------------------------------------------------------------------+
|                       STATELESS COMPUTE ENGINE (API NODE)                    |
|  2. Queries Spatial-Spectral Sidecar Index / Live AWS S3 Bucket:             |
|     - Evaluates WCS (RA---TAN-SIP) & WCS-WAVE polynomials at (RA, Dec)       |
|     - Computes detector roll angle, wavelength, and 24-row byte offsets      |
|     - Mode A (Screening): 3 iso-wavelength pairs (or bracketing pairs)       |
|     - Mode B (Full SED): Multi-exposure sweep across 0.75 to 5.00 um         |
|     - Mode C (Live S3 DR1): Streams 13 live DR1 files across Detectors D1-D6 |
+--------------------------------------+---------------------------------------+
                                       | Persistent HTTP/1.1 Pool (320 Sockets)
                                       v
+------------------------------------------------------------------------------+
|                REMOTE CLOUD OBJECT STORAGE (nasa-irsa-spherex)               |
|  3. Streams 24-Row Spans (195.84 KB/HDU) for SCI, FLAGS, VARIANCE, ZODI      |
+--------------------------------------+---------------------------------------+
                                       | Big-Endian IEEE-754 Binary Streams
                                       v
+------------------------------------------------------------------------------+
|                       STATELESS COMPUTE ENGINE (API NODE)                    |
|  4. Swap32 Hydration -> BSCALE/BZERO -> Bit 2 Variance x2.25 -> Minus ZODI   |
|  5. Executes Three Independent Parallel Detection & Classification Branches: |
|     * Branch 1 (Differential ePSF Astrometry):                               |
|       Fits oversampled ePSF on native pixels + Gaia DR3 tie -> 5-Param Fit   |
|     * Branch 2 (Common-Grid ZOGY & Intra-Pass Photometry):                   |
|       Resamples to 11x11 (u,v) grid -> Whitened Score S_corr -> JPL Veto     |
|       -> Fits Planck B_nu(lambda, T) to Difference Spectrum (with R-J limit) |
|     * Branch 3 (Narrowband Solid Ice vs. Gas Template SED Analysis):         |
|       Evaluates Delta EW at 3.05, 4.27, 4.67 um with BH-FDR & Delta BIC      |
+--------------------------------------+---------------------------------------+
                                       | HTTP 200 OK (COMPLETE or PARTIAL JSON)
                                       v
+------------------------------------------------------------------------------+
|                          CLIENT VISUALIZATION LAYER                          |
|  6. Renders 11x11 WCS Cutouts, 5-Param Astrometric Radar & 102-Band SEDs     |
+------------------------------------------------------------------------------+
```

---

## 3. Core Mathematical & Physical Innovations

### 3.1 Multi-Extension FITS HDU Boundary Stepping (Equation 3.1)

Live SPHEREx Data Release 1 (`l2b-v27`) files store HDU 2 (`FLAGS`) as a tile-compressed binary table (`BITPIX = 8`, `NAXIS = [8, 2040]`) with a 13.3 MB `PCOUNT` heap. My two-stage parser (`parseFitsStructuralHeader` and `extractCelestialWCS`) explicitly accounts for `PCOUNT` and `GCOUNT`, so it never crashes on `NAXIS = 0` or `BINTABLE` extensions:

```text
O_data(n)  = O_HDU(n) + 2880 * ceil( (N_cards(n) * 80) / 2880 )
N_bytes(n) = ( |BITPIX(n)| / 8 ) * GCOUNT(n) * [ PCOUNT(n) + Product_{k=1..NAXIS} NAXIS_k(n) ]
O_HDU(n+1) = O_data(n) + 2880 * ceil( N_bytes(n) / 2880 )
```

### 3.2 Roll-Angle Padded 24-Row Span & Native `Buffer.swap32()` Hydration

An `11 x 11` North-Up tangent-plane cutout projected onto a detector rotated by roll angle theta spans `11 * (|cos(theta)| + |sin(theta)|)` rows (up to 15.56 rows at 45 deg). Adding a Lanczos-3 interpolation kernel radius (3 px per side) and a 1 px SIP distortion margin yields a required span of 24 contiguous detector rows (`24 * 2040 * 4 = 195,840` bytes per HDU).

I decode incoming Big-Endian IEEE-754 streams in-place via Node's C++ `Buffer.swap32()`, then subtract the foreground Zodiacal dust (HDU 4: `ZODI`). Pixels flagged with Bit 2 (`0x0004 NONLINEAR_RESPONSE`) under the default mask `0xFFFB` are retained, with their variance multiplied by `1.5^2 = 2.25`.

### 3.3 Branch 1: Differential ePSF Astrometry & 5-Parameter Parallax Solution

Standard difference-image pipelines miss slow-moving brown dwarfs and outer solar system bodies because a sub-pixel shift at moderate SNR produces difference lobes below 5 sigma. Branch 1 runs independently on every detected source (SNR >= 10):

- It fits an oversampled Effective PSF (FWHM = 1.35 px) directly to un-resampled detector pixels, reducing pixel-phase bias to `sigma_phase = 0.008 px`.
- It solves the 5-parameter barycentric trajectory in local tangent-plane coordinates (xi, eta), using Earth's barycentric position vector `r_earth(t) = [X, Y, Z]^T` in AU:

```text
xi(t_m)  = xi_0  + mu_alpha* * (t_m - t_0) + parallax * P_alpha(t_m)
eta(t_m) = eta_0 + mu_delta  * (t_m - t_0) + parallax * P_delta(t_m)
```

- I center the reference epoch `t_0` at the inverse-variance weighted mean MJD of the observations (`sum(w_m * dt_m) = 0`) and apply Bayesian Gaussian prior regularization on (mu, parallax) for short intra-pass baselines. This guarantees a well-conditioned normal matrix `A^T * W * A` across any cadence.

### 3.4 Branch 2: Covariance-Whitened ZOGY Photometry, Live JPL Veto & Difference-Spectrum Planck Fit

Resampling native detector pixels onto a shared `11 x 11` tangent grid correlates neighboring grid pixels. I propagate the exact `121 x 121` spatial covariance matrix `Sigma_Delta = Sigma_1 + Sigma_2` and evaluate point-source flux changes using the covariance-whitened PSF matched filter:

```text
Delta_F_PSF = (P^T * Sigma_Delta^-1 * Delta_F) / (P^T * Sigma_Delta^-1 * P)

S_corr = Delta_F_PSF / sqrt( (P^T * Sigma_Delta^-1 * P)^-1 )
```

- **Solar System Veto:** I query NASA JPL's live `sb_ident.api` at each exposure's midpoint `MJD_m` and flag step-to-step intra-pass tracklets, so main-belt asteroids are never misclassified as stellar flares.
- **Difference-Spectrum Planck Fitting:** I fit frequency-space blackbody intensity `B_nu(lambda, T)` (`b_nu = 5099.4 um*K`) directly to the Difference Spectrum `Delta_F_nu(lambda)`. When a hot flare (T >= 10,000 K) enters the Rayleigh-Jeans tail in the 0.75 to 5.00 um window, the engine reports a 95% confidence lower limit (`T_min,95% > 10,000 K`).

### 3.5 Branch 3: Solid-State Volatile Ice vs. Cool-Star Gas Templates

I evaluate continuum-subtracted Equivalent Widths (Delta EW) across three key astrobiological ice bands:

- **Water Ice (H2O at 3.05 um):** The blue continuum (2.50 to 2.62 um) and red continuum (3.65 to 3.80 um) are anchored outside the 2.7 um water vapor and 3.3 um T-dwarf methane (CH4) bands.
- **Carbon Dioxide Ice (CO2 at 4.27 um)** and **Carbon Monoxide Ice (CO at 4.67 um).**

I enforce Benjamini-Hochberg False Discovery Rate control (`alpha_FDR = 10^-4`, `|Delta EW| >= 4.5 sigma`) and require a Bayesian Information Criterion advantage `Delta BIC_(ice-gas) > 10` over substellar gas-phase templates.

---

## 4. Interactive Dashboard Modes (Buttons 1 to 6)

- **Preset 1: L/T Brown Dwarf** (mu = 0.77"/yr, parallax = 0.048"). Demonstrates how a sub-pixel mover that produces only a 4.25 sigma single-band ZOGY dipole is captured at 21.86 sigma by Branch 1's multi-exposure 5-parameter ePSF fit, while recovering its 1,700 K temperature.
- **Preset 2: Rayleigh-Jeans Stellar Flare** (T > 10^4 K). Simulates a rapid UV/optical flare on an M-dwarf, isolating the flaring exposures via the Difference Spectrum and reporting a 95% CL lower limit (> 10,000 K).
- **Preset 3: Volatile H2O / CO2 Ice Absorption** (3.05 um). Simulates a cometary outgassing or circumstellar disk event, triggering a PASS on the `H2O_305` and `CO2_427` Benjamini-Hochberg FDR table.
- **Preset 4: Main-Belt Asteroid Veto** (MPC + Intra-Pass). Simulates a fast-moving asteroid crossing a single LVF step in Pass 2, triggering the yellow Solar System Contaminant Veto banner.
- **Preset 5: Vega Saturation Veto** (HTTP 422). Points the viewport at Alpha Lyrae (Vega), where more than 50% of cutout pixels fail Bit 1 (`0x0002 OVERFLOW_SATURATED`), halting the pipeline with HTTP 422 Unprocessable Entity.
- **Mode 6: LIVE AWS S3 DR1 Full 3-Branch Pipeline (D1 to D6).** Connects live to `nasa-irsa-spherex.s3.us-east-1.amazonaws.com`, indexes 36 real Data Release 1 files, slices 7.64 MB across 13 files (Detectors D1 to D6), walks all 5 HDUs live, and fits a real 3,600 K North Ecliptic Pole field star at SNR = 157.8.

---

## 5. Quick Start (Local & GitHub Codespaces)

### Prerequisites

- Node.js v18+ or v20+
- npm v9+

### 1. Start the Backend Engine (Port 4000)

```bash
cd spherex-engine/backend
npm install
npm run dev
```

### 2. Start the Frontend Visualization Console (Port 5173)

Open a second terminal:

```bash
cd spherex-engine/frontend
npm install
npm run dev
```

Open `http://localhost:5173` (or click the forwarded port 5173 link in GitHub Codespaces). Vite automatically proxies `/api/*` requests to `http://localhost:4000`.

### 3. CLI Smoke Tests

```bash
# Test the 3-Branch Synthetic Injection Pipeline
curl -s -X POST http://localhost:4000/api/v1/spectral-delta \
  -H "Content-Type: application/json" \
  -d '{
    "viewport": { "ra_deg": 269.4521, "dec_deg": 66.0234, "cutout_side_arcsec": 67.65, "geometry": "SQUARE_TANGENT_GRID" },
    "execution_mode": "FULL_SED_102_BAND",
    "passes": { "baseline": "PASS_1", "comparison": "PASS_2", "verification": "PASS_3" },
    "simulation_preset": "BROWN_DWARF_PROPER_MOTION"
  }' | head -n 45

# Test the Live AWS S3 SPHEREx DR1 + NASA JPL Small-Body Pipeline
curl -s "http://localhost:4000/api/v1/live-irsa-search" | head -n 45
```

---

## 6. Repository Structure

```text
spherex-engine/
├── backend/
│   ├── package.json
│   ├── tsconfig.json
│   └── src/
│       ├── server.ts                      # Express API (HTTP 200/422/502/504 contracts)
│       ├── types/
│       │   └── astro.types.ts             # Pass registry, bitmasks, WCS & API schemas
│       ├── core/
│       │   ├── fitsParser.ts              # 2-stage FITS parser, Swap32 & Bit-2 calibration
│       │   └── wcsProjection.ts           # TAN-SIP 2D Newton-Raphson & 24-row span math
│       ├── services/
│       │   ├── sidecarIndex.ts            # Pre-indexed WCS-WAVE & barycentric ephemeris
│       │   ├── s3RangeFetcher.ts          # 320-socket HTTP/1.1 Keep-Alive Range fetcher
│       │   ├── syntheticMefGenerator.ts   # Big-Endian IEEE-754 24-row buffer synthesizer
│       │   └── nasaLiveApis.ts            # Live IRSA TAP, AWS S3 failover & JPL sb_ident
│       └── engine/
│           ├── astrometryBranch.ts        # Branch 1: Native ePSF & 5-param (mu, parallax)
│           ├── photometryZogyBranch.ts    # Branch 2: Common-grid covariance & Planck fit
│           └── spectralIceBranch.ts       # Branch 3: 102-band ice EW, BH-FDR & Delta BIC
└── frontend/
    ├── package.json
    ├── vite.config.ts                     # Port 5173 server with /api proxy to port 4000
    └── src/
        ├── App.tsx                        # 4D interactive heatmaps, radar & SED console
        └── App.css                        # Dark-mode astrophysics telemetry styling
```

---

## 7. NASA Data & Partner Resources Used

- **NASA IRSA SPHEREx AWS Open Data S3 Archive** (`nasa-irsa-spherex`): https://registry.opendata.aws/nasa-irsa-spherex/
- **NASA/IPAC Infrared Science Archive (IRSA) SPHEREx Mission Portal:** https://irsa.ipac.caltech.edu/Missions/spherex.html
- **NASA/IPAC IRSA Table Access Protocol (TAP) API:** https://irsa.ipac.caltech.edu/TAP
- **NASA/JPL Solar System Dynamics Small-Body Identification API (`sb_ident`):** https://ssd-api.jpl.nasa.gov/doc/sb_ident.html
- **ESA Gaia Data Release 3 (DR3) Astrometric Catalog:** https://www.cosmos.esa.int/web/gaia/dr3
- **IAU Minor Planet Center (MPC):** https://www.minorplanetcenter.net/

---

## 8. Use of Super Intelligence (SI) Disclosure

In accordance with NASA Space Apps Challenge rules, I transparently disclose that Large Language Model (SI) assistance (Google Gemini) was used during this project as an interactive mathematical sounding board, specification reviewer, and pair-programming assistant for generating TypeScript/React boilerplate, debugging normal-matrix conditioning, and formatting technical documentation. I directed, tested, and verified all system architecture decisions, live NASA AWS S3 byte-range integrations, JPL API verifications, and interactive validation tests. No NASA logos, flags, or mission insignias were generated or modified using SI tools.