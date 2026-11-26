(* ::Package:: *)
(*
  Fisher-Information Inequality for Boolean Pairs
  Exact Mathematica / Wolfram Language symbolic verification file.

  Source tags: green clovers \CASmark{01},...,\CASmark{26}.
  Mathematica 13.x / 14.x recommended. All rational constants are EXACT.

  Instructions:
     Get["/path/to/fisher_boolean_pairs_checks.wl"]
  or wolframscript -file fisher_boolean_pairs_checks.wl.

  Each tag reports PASS, FAIL, or UNRESOLVED. A PASS certifies the specific
  rational identity, derivative, or semialgebraic claim documented here,
  NOT the non-algebraic mathematical reasoning surrounding that claim.

  The independent combinatorial/Fourier arguments, Taylor remainder
  comparisons, global range assumptions, scalar ODE comparison, and
  coverage of parameter regimes are NOT reduced to numerical checks.
*)

ClearAll["Global`*"];
$Assumptions = True;
$HistoryLength = 0;
failedLeaves = {};
SetAttributes[casLeaf, HoldRest];
casLeaf[id_String, assertion_] := Module[{ans},
  ans = Quiet @ Check[
    TimeConstrained[TrueQ[FullSimplify[assertion]], 180, $Aborted],
    $Failed
  ];
  Print["\\CASmark{", id, "} [green clover] : ", Which[TrueQ[ans], "PASS", ans === $Aborted, "TIMEOUT", True, "FAIL/UNRESOLVED"]];
  If[! TrueQ[ans], AppendTo[failedLeaves, id]];
  ans
];
(* This supplies exact real-quantifier elimination, not grid sampling. *)
ClearAll[allReal];
allReal[vars_List, premise_, conclusion_] :=
  Resolve[ForAll[vars, Implies[premise, conclusion]], Reals];

ClearAll[qfun, phi, quad, parity, dstar];
qfun[lam_, d_] := lam (1 - d)/((1 + d) (1 - lam d));
phi[lam_, d_] := qfun[lam, d]/(lam - d);
quad[lam_, d_] := -4 (1 - lam) d^2 + (lam (1 - lam) - 4) d + lam^2 + 3 lam;
parity[k_, r_] := k^2 r^(2 k - 2)/(1 - r^(2 k));
dstar[lam_] := 3 lam/4 + lam^2 (3 lam - 1)/12;

(* Benchmark: direct four-cell computation, included for completeness.
   This identity is not a separate blue leaf because the benchmark itself
   is a primary statement in the Introduction. *)
ClearAll[checkParityBenchmark];
checkParityBenchmark[] := FullSimplify[
  2 ((k r^(k - 1)/4)^2/((1 + r^k)/4)) +
  2 ((k r^(k - 1)/4)^2/((1 - r^k)/4)) == parity[k, r],
  0 < r < 1 && Element[k, Integers] && k >= 1
];
Print["Parity Fisher benchmark: ", checkParityBenchmark[]];

(* ===== Curvature and zero-covariance normalization ===== *)
(* \CASmark{01}: curvature at independence. *)
casLeaf["01",
  Together[p (1 - p) qq (1 - qq) * (
    1/(p qq) + 1/(p (1 - qq)) +
    1/((1 - p) qq) + 1/((1 - p) (1 - qq))
  ) - 1] == 0
];

(* \CASmark{02}: beta in terms of the two endpoint parameters. *)
ClearAll[bet];
bet[lam_, eta_] := lam eta/((lam + eta) (1 + lam eta));
casLeaf["02",
  Together[bet[lam, eta] (lam + 1/lam + eta + 1/eta) - 1] == 0
];

(* \CASmark{03}: sum of the four reciprocal cell probabilities. *)
ClearAll[wcurv];
wcurv[lam_, eta_, d_] := bet[lam, eta] (
  1/(lam - d) + 1/(1/lam - d) +
  1/(eta + d) + 1/(1/eta + d)
);
casLeaf["03",
  Together[wcurv[lam, eta, d] -
    lam eta (1 - d^2)/((lam - d) (1 - lam d) (eta + d) (1 + eta d))] == 0
];

(* \CASmark{04}: curvature factorization. *)
casLeaf["04",
  Together[wcurv[lam, eta, d] -
    phi[lam, d] eta (1 + d)^2/((eta + d) (1 + eta d))] == 0
];

(* \CASmark{05}: nonnegative difference that eliminates eta. *)
casLeaf["05",
  Expand[(eta + d) (1 + eta d) - eta (1 + d)^2 - d (1 - eta)^2] == 0 &&
  allReal[{eta, d}, 0 < eta <= 1 && 0 < d < 1,
    (eta + d) (1 + eta d) - eta (1 + d)^2 >= 0]
];

(* ===== Normalized curvature and spectral cutoffs ===== *)
(* \CASmark{06}: monotonicity of Q on 0<d<lambda<1. *)
casLeaf["06",
  Together[D[qfun[lam, d], lam] -
    (1 - d)/((1 + d) (1 - lam d)^2)] == 0 &&
  Together[D[qfun[lam, d], d] -
    lam (-2 + lam (1 + 2 d - d^2))/((1 + d)^2 (1 - lam d)^2)] == 0 &&
  allReal[{lam, d}, 0 < d < lam < 1,
    0 < qfun[lam, d] <= lam &&
    D[qfun[lam, d], lam] > 0 && D[qfun[lam, d], d] < 0]
];

(* \CASmark{07}: equality between the spectral test and the quadratic. *)
casLeaf["07",
  Together[
    1 - (((1 + lam)^2/4) - d^2) phi[lam, d] -
    (1 - lam) quad[lam, d]/(4 (1 + d) (lam - d) (1 - lam d))
  ] == 0
];

(* \CASmark{08}: three exact cutoff identities and strict positivity
   of the small-ceiling quadratic. *)
cc0 = 741/1000;
casLeaf["08",
  Together[quad[lam, cc0 lam]/lam -
    (363831 lam^2 - 113831 lam + 9000)/250000] == 0 &&
  Together[quad[lam, 2 lam - 1] - 14 lam (1 - lam)^2] == 0 &&
  Together[quad[lam, 3 lam/4] - lam^2 (3 lam - 1)/2] == 0 &&
  Discriminant[363831 lam^2 - 113831 lam + 9000, lam] < 0
];

(* \CASmark{09}: extended spectral cutoff and bounds on the derivative. *)
casLeaf["09",
  allReal[{lam, d}, 1/3 <= lam < 1 && 0 <= d <= lam,
    3 lam/4 <= dstar[lam] < lam &&
    -D[quad[lam, d], d] < 6 &&
    D[quad[lam, d], d] < 0] &&
  allReal[{lam}, 1/3 <= lam < 1, quad[lam, dstar[lam]] >= 0]
];

(* ===== Negative-slope algebra ===== *)
(* \CASmark{10}: derivative in lambda and the wide-gap quadratic. *)
ll0 = (1 + 7 d)/8;
casLeaf["10",
  Together[D[phi[lam, d], lam] -
    (-d (1 - d) (1 - lam^2)/((1 + d) (lam - d)^2 (1 - lam d)^2))] == 0 &&
  Together[(1 + d) (1 - ll0 d) - 4 ll0 (1 - d)^2 -
    ((1 - d)/8) (35 d^2 - 9 d + 4)] == 0 &&
  Discriminant[35 d^2 - 9 d + 4, d] < 0
];

(* \CASmark{11}: narrow-gap discriminant certificate. *)
llw = d + w (1 - d)/8;
ngap = (1 + d) (1 - llw d) -
  8 llw (1 - d) ((1 - d) - 4 w (1 - d)/8);
polyNeg = (72 - 41 w + 4 w^2) d^2 +
  (-48 + 47 w - 8 w^2) d + (8 - 8 w + 4 w^2);
casLeaf["11",
  Together[8 ngap/(1 - d) - polyNeg] == 0 &&
  Factor[Discriminant[polyNeg, d]] == w (32 w^2 + 385 w - 896) &&
  allReal[{w}, 0 <= w <= 1,
    72 - 41 w + 4 w^2 > 0 &&
    w (32 w^2 + 385 w - 896) <= 0]
];

(* ===== Positive-slope, small correlation ===== *)
(* \CASmark{12}: rational expression for derivative of log gap. *)
casLeaf["12",
  Together[D[x - Log[1 + x] - x^2/(2 + 4 x/3), x] -
    x^3/((1 + x) (3 + 2 x)^2)] == 0
];

(* \CASmark{13}: all exact coefficients and positivity of the small cubic. *)
smallCub[lam_] := (1 + cc0 lam) (1 - cc0 lam^2) -
  lam (1 - cc0 lam) (6 - 20 cc0 lam/3);
casLeaf["13",
  Together[smallCub[lam] - (
    1 - 5259 lam/1000 + 1729 lam^2/200 - 4209621 lam^3/1000000
  )] == 0 &&
  allReal[{lam}, 0 <= lam <= 1/2, smallCub[lam] > 0]
];

(* ===== Short scale ===== *)
(* \CASmark{14}: Q at d=2lambda-1 and endpoint comparison. *)
casLeaf["14",
  Together[(1 + lam) qfun[lam, 2 lam - 1] -
    (1 + lam)/(1 + 2 lam)] == 0 &&
  allReal[{lam}, 1/2 <= lam < 1,
    (1 + lam)/(1 + 2 lam) <= 3/4]
];

(* ===== Normalized parity benchmark ===== *)
(* \CASmark{15}: exact rational part of log(3) bound.
   The analytic bound log(3)>263/240 follows from its
   positive-term power series at 1/2, as written in the paper. *)
casLeaf["15",
  (263/240)^2 - 6/5 == 49/57600 &&
  FullSimplify[Log[3] > 263/240]
];

(* ===== Outer-ceiling affine estimates ===== *)
ClearAll[planeInc, planeDec];
planeInc[d_, z_] := (13/12) (1 - d) + z/4;
planeDec[d_, z_, aa_, bb_] :=
  (13/12) (1 - aa d) + (4/3 - 13 bb/12) z;
(* \CASmark{16}: affine-plane identities and rational norm bounds. *)
casLeaf["16",
  (3/5)^2 + 2^2/(20/3) == 24/25 &&
  (3/4)^2 + (5/3)^2/(55/8) == 1531/1584 &&
  Together[planeInc[d, 2 d/5] - planeDec[d, 2 d/5, 3/5, 2]] == 0 &&
  Together[planeInc[d, 2 d/5] - (65 - 59 d)/60] == 0 &&
  Together[planeInc[d, 3 d/8] - planeDec[d, 3 d/8, 3/4, 5/3]] == 0 &&
  Together[planeInc[d, 3 d/8] - (104 - 95 d)/96] == 0
];

(* \CASmark{17}: exact polynomial obtained in outer small ceiling. *)
outerSmall[lam_] := 81 - 531 lam + 690 lam^2 - 272 lam^3;
casLeaf["17",
  Together[
    60 (1 + 2 lam/3) (1 - 2 lam^2/3) -
    (20/3) lam (1 - 2 lam/3) (65 - 59 (2 lam/3)) -
    (20/27) outerSmall[lam]
  ] == 0 &&
  outerSmall[1/5] == 28/125 &&
  allReal[{lam}, 0 <= lam <= 1/5,
    D[outerSmall[lam], lam] < 0 && outerSmall[lam] > 0]
];

(* \CASmark{18}: exact polynomial obtained in outer large ceiling. *)
ClearAll[outerLargeExpr, outerLargeP];
outerLargeExpr[lam_] := 96 (1 + 5 lam/6) (1 - 5 lam^2/6) -
   (55/8) lam (1 - 5 lam/6) (104 - 95 (5 lam/6));
outerLargeP[y_] := 160 + 17952 y - 10854 y^2 - 5993 y^3;
casLeaf["18",
  Together[1440 outerLargeExpr[(y + 4)/5] - outerLargeP[y]] == 0 &&
  allReal[{y}, 0 <= y <= 1, outerLargeP[y] >= 160 + 1105 y > 0]
];

(* ===== Central-ceiling plane and concavity certificates ===== *)
(* Shared rational constants for the central-ceiling checks. *)
mm = 55/8;
kappa = 13/12;

(* \CASmark{19}: square-root concavity numerator, eta substitution
   after taking both derivatives with lambda fixed. *)
ClearAll[vaux, etaSub, jpoly];
vaux[lam_, d_] := (lam - d) qfun[lam, d];
etaSub[lam_, d_] := (lam - d)/(1 - lam d);
jpoly[d_, eta_] := -1 + 4 eta (1 - d) +
  (2 - 8 d) eta^2 - 4 (1 - d) eta^3 + 3 eta^4;
casLeaf["19",
  Together[
    2 vaux[lam, d] D[vaux[lam, d], {d, 2}] -
    D[vaux[lam, d], d]^2 -
    (d + etaSub[lam, d])^2/
      ((1 + d)^4 (1 + d etaSub[lam, d])^2) *
      jpoly[d, etaSub[lam, d]]
  ] == 0
];

(* \CASmark{20}: rational domain restrictions and negativity in both
   d>=1/4 and d<=1/4 intervals. This is exact QE (possibly slower).
   It does not prove the functional analytic differentiation step. *)
condCentral = 1/5 <= lam <= 4/5 && 37 lam/50 <= d < lam;
casLeaf["20",
  allReal[{lam, d}, condCentral,
    etaSub[lam, d] <= 1 - d &&
    Implies[d >= 1/4,
      jpoly[d, etaSub[lam, d]] <= -21/125] &&
    Implies[d <= 1/4,
      etaSub[lam, d] < 1/8 &&
      jpoly[d, etaSub[lam, d]] < -15/32]
  ]
];

(* \CASmark{21}: exact partial fractions for rational convexity. *)
ClearAll[ratConv];
ratConv[lam_, d_, aa_] := qfun[lam, d] (kappa - aa d);
casLeaf["21",
  Together[ratConv[lam, d, aa] - (
    -aa + 2 lam (aa + kappa)/((1 + lam) (1 + d)) +
    (aa - kappa lam) (1 - lam)/((1 + lam) (1 - lam d))
  )] == 0
];

(* \CASmark{22}: exact univariate and bivariate rational estimates
   supporting convexity for a=2/3, 8/13<=lambda<=4/5. *)
casLeaf["22",
  allReal[{lam}, 8/13 <= lam <= 4/5,
    (kappa lam - 2/3) (1 - lam) lam <= 4/125] &&
  allReal[{lam, d}, 8/13 <= lam <= 4/5 && 0 <= d <= 3/4,
    (1 + d)/(1 - lam d) <= 35/8 &&
    D[ratConv[lam, d, 2/3], {d, 2}] >= 0] &&
  (4/125) (35/8)^3 == 343/128 && 343/128 < 7/2
];

(* \CASmark{23}: central G at zero relative gap. *)
gZero[lam_] := 1 - mm qfun[lam, lam] (kappa - lam);
casLeaf["23",
  Together[gZero[lam] -
    (756 lam^2 - 523 lam + 96)/(96 (1 + lam)^2)] == 0 &&
  Discriminant[756 lam^2 - 523 lam + 96, lam] == -16775
];

(* \CASmark{24}: lower bounds for three square-root endpoint factors;
   we square only after verifying both sides are nonnegative. *)
ClearAll[radDiff];
radDiff[ss_, kk_, lam_] :=
  (1 - ss^2 lam^2) (1 - ss lam^2) - (1 - kk lam^2)^2;
casLeaf["24",
  mm/16 > (13/20)^2 && mm (13/50) > (4/3)^2 &&
  allReal[{lam}, 1/5 <= lam <= 4/5,
    1 - (11/12) lam^2 > 0 &&
    1 - (13/20) lam^2 > 0 &&
    radDiff[15/16, 11/12, lam] >= 0 &&
    radDiff[37/50, 13/20, lam] >= 0]
];

(* \CASmark{25}: coefficients and exact positivity of the three cubics.
   These three Reduce/Resolve checks replace the page-long expansions. *)
ClearAll[kfun, endpointCubic, p1, p2, p3];
kfun[ss_, lam_] := (1 + ss lam) (1 - ss lam^2);
endpointCubic[aa_, bb_, yy_, v0_, kk_, lam_] := Module[{ss = 1 - yy},
  kfun[ss, lam] - mm lam (1 - ss lam) (kappa - aa ss lam) +
    bb v0 lam (1 - kk lam^2)
];
p1[lam_] := 1 - (970048/153600) lam +
  (1918500/153600) lam^2 - (1090581/153600) lam^3;
p2[lam_] := 1 - (17504/3072) lam +
  (31770/3072) lam^2 - (17363/3072) lam^3;
p3[lam_] := 1 - (604950/120000) lam +
  (979575/120000) lam^2 - (496892/120000) lam^3;
casLeaf["25",
  Together[endpointCubic[1, 3/10, 1/16, 13/20, 11/12, lam] - p1[lam]] == 0 &&
  Together[endpointCubic[2/3, 5/4, 1/16, 13/20, 11/12, lam] - p2[lam]] == 0 &&
  Together[endpointCubic[2/3, 5/4, 13/50, 4/3, 13/20, lam] - p3[lam]] == 0 &&
  allReal[{lam}, 1/5 <= lam <= 4/5, p1[lam] > 0] &&
  allReal[{lam}, 1/5 <= lam <= 4/5, p2[lam] > 0] &&
  allReal[{lam}, 1/5 <= lam <= 4/5, p3[lam] > 0]
];

(* ===== Logarithmic path / Riccati denominator ===== *)
(* \CASmark{26}: derivative and positivity of R denominator. *)
ClearAll[dDen, rr];
dDen[uu_, ee_, oo_] := 13 ee/12 + 4 uu/3 - uu^2/(oo - 4 uu/5);
rr[uu_, ee_, oo_] := uu^2/dDen[uu, ee, oo];
casLeaf["26",
  Together[2 dDen[uu, ee, oo] - uu D[dDen[uu, ee, oo], uu] -
    (13 ee/6 + 4 uu/3 + (4/5) uu^3/(oo - 4 uu/5)^2)] == 0 &&
  Together[D[rr[uu, ee, oo], uu] -
    uu (2 dDen[uu, ee, oo] - uu D[dDen[uu, ee, oo], uu])/
     dDen[uu, ee, oo]^2] == 0 &&
  allReal[{ss, ee, oo}, 0 <= ss <= 2/5 && ee > 0 && oo > 0,
    13 ee/(12 oo) + ss (4/3 - ss/(1 - 4 ss/5)) > 0]
];

Print["----------------------------------------"];
Print["Total checks: 26"];
Print["Failed/unresolved clover marks: ", failedLeaves];
Print[If[failedLeaves === {}, "ALL SYMBOLIC CHECKS PASSED", "INSPECT RESULTS BEFORE RELYING ON CERTIFICATES"]];
