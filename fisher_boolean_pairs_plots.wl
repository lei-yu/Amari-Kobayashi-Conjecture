(* ::Package:: *)
(*
  Supplementary figure source for the Fisher-information Boolean-pairs paper.

  All six plots illustrate exact statements in the main companion file
  fisher_boolean_pairs_checks.wl. Plotting does not certify continuum claims:
  the proof certificates are the exact rational and quantifier-elimination
  statements tagged by the green clovers \CASmark{07}, \CASmark{08}, ... .

  To generate the figures in the directory ./figures, run:
     Get["/path/to/fisher_boolean_pairs_plots.wl"];
  or:
     wolframscript -file fisher_boolean_pairs_plots.wl

  The supplied figure PDF files are independently generated vector previews
  from the same formulas. Mathematica can regenerate/replace them here.
*)

ClearAll[casPlotBaseDir, casPlotOutDir, casExport,
  casQuad, casSmallN, casOuterSmall, casOuterLarge,
  casJ, casCurvature, casP1, casP2, casP3,
  casColors];

casPlotBaseDir = If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]], Directory[]];
casPlotOutDir = FileNameJoin[{casPlotBaseDir, "figures"}];
If[! DirectoryQ[casPlotOutDir], CreateDirectory[casPlotOutDir]];
casExport[file_String, p_] := Module[{target},
  target = FileNameJoin[{casPlotOutDir, file <> ".pdf"}];
  Export[target, p, "PDF"];
  Print["Exported figure: ", target];
  target
];
casColors = {
  RGBColor[39/255, 90/255, 141/255],
  RGBColor[34/255, 105/255, 68/255],
  RGBColor[138/255, 75/255, 145/255]
};


(* Supplementary Fig. S1: \CASmark{07}, \CASmark{08}, \CASmark{09}.
   x = d/lambda, lambda = 1/5, 1/2, 4/5. *)
casQuad[lam_, d_] := -4 (1 - lam) d^2 +
  (lam (1 - lam) - 4) d + lam^2 + 3 lam;
casExport["cas_spectral_cutoffs",
  Plot[Evaluate[{
      casQuad[1/5, x/5]/(1/5),
      casQuad[1/2, x/2]/(1/2),
      casQuad[4/5, (4/5) x]/(4/5)
    }], {x, 0, 1},
    Frame -> True, ImageSize -> 610, GridLines -> Automatic,
    FrameStyle -> Directive[Black, 11], LabelStyle -> Directive[Black, 12],
    PlotRangePadding -> Scaled[0.03],
    PlotStyle -> (Directive[#, AbsoluteThickness[2.3]] & /@ casColors),
    PlotRange -> {{0, 1}, {-1.6, 4.2}},
    FrameLabel -> {"d/lambda", "Q_lambda(d) / lambda"},
    PlotLegends -> {"lambda = 1/5", "lambda = 1/2", "lambda = 4/5"},
    Epilog -> {Directive[RGBColor[.7, .35, .17], Dashed, AbsoluteThickness[1.5]],
      Line[{{741/1000, -1.6}, {741/1000, 4.2}}]}
  ]
];

(* Supplementary Fig. S2: \CASmark{13}, positivity on 0<=lambda<=1/2. *)
casSmallN[lam_] := 1 - (5259/1000) lam + (1729/200) lam^2 -
  (4209621/1000000) lam^3;
casExport["cas_small_cubic",
  Plot[casSmallN[lam], {lam, 0, 1/2},
    Frame -> True, ImageSize -> 610, GridLines -> Automatic,
    FrameStyle -> Directive[Black, 11], LabelStyle -> Directive[Black, 12],
    PlotRangePadding -> Scaled[0.03],
    PlotStyle -> Directive[casColors[[1]], AbsoluteThickness[2.5]],
    PlotPoints -> 90, MaxRecursion -> 5,
    ScalingFunctions -> {"Linear", "Log"},
    PlotRange -> {{0, 1/2}, {10^-4, 3/2}},
    FrameLabel -> {"lambda", "N(lambda) [log scale]"}
  ]
];

(* Supplementary Fig. S3: \CASmark{17}, \CASmark{18}.
   The two curves are independently normalized by their initial values. *)
casOuterSmall[lam_] := 81 - 531 lam + 690 lam^2 - 272 lam^3;
casOuterLarge[t_] := 160 + 17952 t - 10854 t^2 - 5993 t^3;
casExport["cas_outer_cubics",
  Plot[{casOuterSmall[t/5]/81, casOuterLarge[t]/160}, {t, 0, 1},
    Frame -> True, ImageSize -> 610, GridLines -> Automatic,
    FrameStyle -> Directive[Black, 11], LabelStyle -> Directive[Black, 12],
    PlotRangePadding -> Scaled[0.03],
    PlotStyle -> {Directive[casColors[[2]], AbsoluteThickness[2.5]],
                  Directive[casColors[[1]], AbsoluteThickness[2.5]]},
    ScalingFunctions -> {"Linear", "Log"},
    PlotRange -> {{0, 1}, {10^-3, 100}},
    FrameLabel -> {"t", "Normalized values [log scale]"},
    PlotLegends -> {"P_small(t/5) / 81", "P_large(t) / 160"}
  ]
];

(* Supplementary Fig. S4: \CASmark{19}, \CASmark{20}.
   eta=(lambda-d)/(1-lambda d), d=lambda*s, s in [37/50,1]. *)
casJ[d_, eta_] := -1 + 4 eta (1 - d) + (2 - 8 d) eta^2 -
  4 (1 - d) eta^3 + 3 eta^4;
casExport["cas_central_J_negative",
  DensityPlot[
    casJ[lam s, (lam - lam s)/(1 - lam^2 s)],
    {lam, 1/5, 4/5}, {s, 37/50, 1},
    Frame -> True, ImageSize -> 610,
    FrameLabel -> {"lambda", "d/lambda"},
    LabelStyle -> Directive[Black, 12],
    PlotPoints -> 55, MaxRecursion -> 2,
    PlotRange -> {-1, -3/5}, ColorFunction -> (Blend[{RGBColor[.17, .18, .45], RGBColor[.11, .62, .61], RGBColor[.90, .85, .24]}, #] &),
    PlotLegends -> Automatic
  ]
];

(* Supplementary Fig. S5: \CASmark{22}.
   Exact second derivative from the partial fractions in the paper. *)
casCurvature[lam_, d_] := 4 lam (2/3 + 13/12)/
    ((1 + lam) (1 + d)^3) +
  2 (2/3 - 13 lam/12) (1 - lam) lam^2/
    ((1 + lam) (1 - lam d)^3);
casExport["cas_central_rational_convexity",
  DensityPlot[casCurvature[lam, d],
    {lam, 8/13, 4/5}, {d, 0, 3/4},
    Frame -> True, ImageSize -> 610,
    FrameLabel -> {"lambda", "d"},
    LabelStyle -> Directive[Black, 12],
    PlotPoints -> 55, MaxRecursion -> 2,
    PlotRange -> All, ColorFunction -> (Blend[{RGBColor[.10, .24, .43], RGBColor[.45, .48, .56], RGBColor[.95, .81, .34]}, #] &),
    PlotLegends -> Automatic
  ]
];

(* Supplementary Fig. S6: \CASmark{25}, each cubic positive for
   lambda in [1/5,4/5] by exact rational quantifier elimination. *)
casP1[lam_] := 1 - (970048/153600) lam +
  (1918500/153600) lam^2 - (1090581/153600) lam^3;
casP2[lam_] := 1 - (17504/3072) lam +
  (31770/3072) lam^2 - (17363/3072) lam^3;
casP3[lam_] := 1 - (604950/120000) lam +
  (979575/120000) lam^2 - (496892/120000) lam^3;
casExport["cas_central_cubics",
  Plot[Evaluate[{casP1[lam], casP2[lam], casP3[lam]}],
    {lam, 1/5, 4/5},
    Frame -> True, ImageSize -> 610, GridLines -> Automatic,
    FrameStyle -> Directive[Black, 11], LabelStyle -> Directive[Black, 12],
    PlotRangePadding -> Scaled[0.03],
    PlotStyle -> (Directive[#, AbsoluteThickness[2.5]] & /@ casColors),
    ScalingFunctions -> {"Linear", "Log"},
    PlotPoints -> 90, MaxRecursion -> 4,
    PlotRange -> {{1/5, 4/5}, {10^-3, 1}},
    FrameLabel -> {"lambda", "P_j(lambda) [log scale]"},
    PlotLegends -> {"P1", "P2", "P3"}
  ]
];
Print["All six supplementary figure exports attempted."];
