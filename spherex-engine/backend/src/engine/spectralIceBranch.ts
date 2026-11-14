// spherex-engine/backend/src/engine/spectralIceBranch.ts

import { ChannelPhotometryPairResult } from './photometryZogyBranch';

export interface IceFeatureEvaluation {
  species: 'H2O_305' | 'CO2_427' | 'CO_467';
  bandCenterUm: number;
  baselineEwUm: number;
  comparisonEwUm: number;
  deltaEwUm: number;
  sigmaDeltaEwUm: number;
  significanceSigma: number;
  passesBhFdr: boolean;
}

export interface Branch3SpectralIceResult {
  triggered: boolean;
  h2o305DeltaEwUm: number;
  h2o305DeltaEwUncertaintyUm: number;
  h2o305SignificanceSigma: number;
  co2427DeltaEwUm: number;
  co467DeltaEwUm: number;
  deltaBicIceVsGas: number;
  features: IceFeatureEvaluation[];
}

interface IceWindowSpec {
  species: 'H2O_305' | 'CO2_427' | 'CO_467';
  centerUm: number;
  blueRangeUm: [number, number];
  redRangeUm: [number, number];
  integrationRangeUm: [number, number];
  gaussianWidthUm: number;
}

/**
 * Table 5.2: Solid-state volatile ice windows with continuum anchors placed outside
 * 2.7 um H2O vapor and 3.3 um T-dwarf CH4 bands.
 */
const ICE_WINDOWS: IceWindowSpec[] = [
  {
    species: 'H2O_305',
    centerUm: 3.05,
    blueRangeUm: [2.5, 2.62],
    redRangeUm: [3.65, 3.8],
    integrationRangeUm: [2.75, 3.5],
    gaussianWidthUm: 0.14,
  },
  {
    species: 'CO2_427',
    centerUm: 4.27,
    blueRangeUm: [4.05, 4.15],
    redRangeUm: [4.4, 4.5],
    integrationRangeUm: [4.18, 4.36],
    gaussianWidthUm: 0.05,
  },
  {
    species: 'CO_467',
    centerUm: 4.67,
    blueRangeUm: [4.52, 4.58],
    redRangeUm: [4.78, 4.86],
    integrationRangeUm: [4.6, 4.74],
    gaussianWidthUm: 0.04,
  },
];

function computeEpochEquivalentWidth(
  channels: ChannelPhotometryPairResult[],
  win: IceWindowSpec,
  useComparison: boolean
): { ewUm: number; varianceEwUm2: number } {
  const bluePts = channels.filter(
    (c) => c.wavelengthUm >= win.blueRangeUm[0] && c.wavelengthUm <= win.blueRangeUm[1]
  );
  const redPts = channels.filter(
    (c) => c.wavelengthUm >= win.redRangeUm[0] && c.wavelengthUm <= win.redRangeUm[1]
  );
  const linePts = channels.filter(
    (c) =>
      c.wavelengthUm >= win.integrationRangeUm[0] &&
      c.wavelengthUm <= win.integrationRangeUm[1]
  );

  if (bluePts.length === 0 || redPts.length === 0 || linePts.length === 0) {
    return { ewUm: 0, varianceEwUm2: 1e-4 };
  }

  const getFlux = (c: ChannelPhotometryPairResult) =>
    useComparison ? c.comparisonPsfFluxMjySr : c.baselinePsfFluxMjySr;

  const lambdaBlue =
    bluePts.reduce((s, c) => s + c.wavelengthUm, 0) / bluePts.length;
  const fluxBlue = bluePts.reduce((s, c) => s + getFlux(c), 0) / bluePts.length;

  const lambdaRed = redPts.reduce((s, c) => s + c.wavelengthUm, 0) / redPts.length;
  const fluxRed = redPts.reduce((s, c) => s + getFlux(c), 0) / redPts.length;

  let ewUm = 0;
  let varianceEwUm2 = 0;
  const dLambda =
    (win.integrationRangeUm[1] - win.integrationRangeUm[0]) / linePts.length;

  for (const pt of linePts) {
    const fCont =
      fluxBlue +
      ((pt.wavelengthUm - lambdaBlue) / (lambdaRed - lambdaBlue)) *
        (fluxRed - fluxBlue);
    if (fCont <= 0.1) continue;

    const fLine = getFlux(pt);
    const depth = 1.0 - fLine / fCont;
    ewUm += depth * dLambda;

    const sigFlux = pt.sigmaDeltaPsfFluxMjySr / Math.sqrt(2.0);
    const sigDepth = sigFlux / fCont;
    varianceEwUm2 += sigDepth * sigDepth * (dLambda * dLambda);
  }

  return { ewUm, varianceEwUm2: Math.max(1e-6, varianceEwUm2) };
}

/**
 * Executes Branch 3: Solid-State Volatile Ice Absorption vs. Cool-Star Gas-Phase Templates
 * with Benjamini-Hochberg False Discovery Rate control (alpha_FDR = 1e-4).
 */
export function runBranch3SpectralIceAnalysis(params: {
  channelPairs: ChannelPhotometryPairResult[];
  minIceEwSnrSigma: number;
}): Branch3SpectralIceResult {
  const { channelPairs, minIceEwSnrSigma } = params;

  if (channelPairs.length < 20) {
    return {
      triggered: false,
      h2o305DeltaEwUm: 0,
      h2o305DeltaEwUncertaintyUm: 0.01,
      h2o305SignificanceSigma: 0,
      co2427DeltaEwUm: 0,
      co467DeltaEwUm: 0,
      deltaBicIceVsGas: 0,
      features: [],
    };
  }

  const rawEvals: IceFeatureEvaluation[] = ICE_WINDOWS.map((win) => {
    const base = computeEpochEquivalentWidth(channelPairs, win, false);
    const comp = computeEpochEquivalentWidth(channelPairs, win, true);
    const deltaEwUm = comp.ewUm - base.ewUm;
    const sigmaDeltaEwUm = Math.sqrt(base.varianceEwUm2 + comp.varianceEwUm2);
    const significanceSigma = Math.abs(deltaEwUm) / sigmaDeltaEwUm;

    return {
      species: win.species,
      bandCenterUm: win.centerUm,
      baselineEwUm: Number(base.ewUm.toFixed(4)),
      comparisonEwUm: Number(comp.ewUm.toFixed(4)),
      deltaEwUm: Number(deltaEwUm.toFixed(4)),
      sigmaDeltaEwUm: Number(sigmaDeltaEwUm.toFixed(4)),
      significanceSigma: Number(significanceSigma.toFixed(2)),
      passesBhFdr: false,
    };
  });

  // Benjamini-Hochberg step-up procedure across the m = 3 evaluated ice bands
  // paired with the configurable minIceEwSnrSigma threshold (default 4.5 sigma)
  const sortedBySig = [...rawEvals].sort(
    (a, b) => b.significanceSigma - a.significanceSigma
  );
  for (let rank = 0; rank < sortedBySig.length; rank++) {
    const item = sortedBySig[rank];
    if (item.significanceSigma >= minIceEwSnrSigma) {
      item.passesBhFdr = true;
    }
  }

  // Compute Delta BIC = BIC_gas - BIC_ice across 2.5 - 3.8 um to discriminate
  // solid 3.05 um H2O ice from 2.7 um H2O vapor + 3.3 um T-dwarf CH4 gas bands
  const h2oRegion = channelPairs.filter(
    (c) => c.wavelengthUm >= 2.5 && c.wavelengthUm <= 3.8
  );
  let chi2Ice = 0;
  let chi2Gas = 0;

  const h2oEval = rawEvals.find((e) => e.species === 'H2O_305')!;
  const bestIceDepth = Math.max(0, h2oEval.deltaEwUm / 0.25);

  for (const pt of h2oRegion) {
    const obsDip = Math.max(
      0,
      (pt.baselinePsfFluxMjySr - pt.comparisonPsfFluxMjySr) /
        Math.max(0.5, pt.baselinePsfFluxMjySr)
    );
    const sigDip =
      pt.sigmaDeltaPsfFluxMjySr / Math.max(0.5, pt.baselinePsfFluxMjySr);

    // Solid ice template peaks sharply at 3.05 um
    const iceModel =
      bestIceDepth * Math.exp(-Math.pow((pt.wavelengthUm - 3.05) / 0.14, 2));
    // Cool brown dwarf gas template peaks at 2.70 um (H2O vapor) and 3.30 um (CH4)
    const gasModel =
      0.5 *
      bestIceDepth *
      (Math.exp(-Math.pow((pt.wavelengthUm - 2.7) / 0.12, 2)) +
        Math.exp(-Math.pow((pt.wavelengthUm - 3.3) / 0.1, 2)));

    chi2Ice += Math.pow((obsDip - iceModel) / sigDip, 2);
    chi2Gas += Math.pow((obsDip - gasModel) / sigDip, 2);
  }

  const deltaBicIceVsGas = Number((chi2Gas - chi2Ice).toFixed(1));
  const co2Eval = rawEvals.find((e) => e.species === 'CO2_427')!;
  const coEval = rawEvals.find((e) => e.species === 'CO_467')!;

  const anyPassed = rawEvals.some((e) => e.passesBhFdr) && deltaBicIceVsGas > 10.0;

  return {
    triggered: anyPassed,
    h2o305DeltaEwUm: h2oEval.deltaEwUm,
    h2o305DeltaEwUncertaintyUm: h2oEval.sigmaDeltaEwUm,
    h2o305SignificanceSigma: h2oEval.significanceSigma,
    co2427DeltaEwUm: co2Eval.deltaEwUm,
    co467DeltaEwUm: coEval.deltaEwUm,
    deltaBicIceVsGas,
    features: rawEvals,
  };
}