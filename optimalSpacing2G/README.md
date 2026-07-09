# optimalSpacing2G

Two-group extension of the `secrdesign` optimal spacing idea, searching over **absolute trap spacing (m)** to optimise a user-chosen criterion based on expected captures/recaptures for two groups with their own `D`, `lambda0`, `sigma`, and detection function.

## What it does
- Scales a baseline trap layout (`traps0`) to candidate spacings (metres), rebuilds a mask, and for each group computes expected captures (`En`) and recaptures (`Er`) with `Enrm`.
- Objective options:
  - `criterion = "sum_min"`: maximise `min(En1 + En2, Er1 + Er2)`.
  - `criterion = "all_min"`: maximise `min(En1, En2, Er1, Er2)`.
  - `criterion = "min_mean_CV"`: maximise `-mean(CV1, CV2)`.
  - `criterion = "min_max_CV"`: maximise `-max(CV1, CV2)`, the negative worst group CV.
  - `criterion = "max_mean_En2"`: maximise `En2plus1 + En2plus2`, where `En2plus` is from `secrdesign::En2()`.
  - `criterion = "max_min_En2"`: maximise `min(En2plus1, En2plus2)`.
- Group CV is computed as `CV_g = 1 / sqrt(min(En_g, Er_g))`.
- `En2plus_g` is the expected number of group `g` individuals detected at two or more detectors.
- Poisson only; no simulation block. Uses `secrdesign:::dfcast`/`replacedefaults` internally to harmonise detection functions/pars.
- Returns a data frame of per-spacing values and the spacing that optimises the chosen criterion.

## Usage
```r
library(secrdesign)
source("optimalSpacing2G/optimalSpacing2G.R")

grid <- make.grid(7, 7, spacing = 50)  # baseline traps (multi-catch)

# One-group reference (secrdesign native)
os1 <- optimalSpacing(
  D = 5,
  traps = grid,
  detectpar = list(lambda0 = 0.2, sigma = 40),
  noccasions = 5,
  plt = FALSE
)
os1$rotRSE$optimum.spacing  # e.g., ~57 m

# Two groups; identical params should recover the same spacing
os2 <- optimalSpacing2G(
  D1 = 5, D2 = 5,
  traps0 = grid,
  detectpar1 = list(lambda0 = 0.2, sigma = 40),
  detectpar2 = list(lambda0 = 0.2, sigma = 40),
  noccasions = 5,
  spacing_m = seq(40, 120, 1),
  criterion = "sum_min"
)
os2$optimum.spacing         # ~57 m, matches one-group result
head(os2$values)            # includes spacing, En1, En2, Er1, Er2, CV1, CV2, En2plus1, En2plus2, crit
```

## Notes
- `spacing_m` sets the candidate spacings (metres); choose the range/resolution you need. The optimum is selected from the evaluated candidate spacings with the maximum `crit`.
- All criteria maximise `crit`; `min_mean_CV` and `min_max_CV` use negative CV values so that maximising `crit` is equivalent to minimising CV.
- `mask` buffer defaults to `xsigma * max(sigma1, sigma2) + spacing`; adjust inside `getCrit2G_abs` if you prefer a different rule.
- Trap scaling is done by multiplying `x`/`y` and spacing attributes by the requested factor; detector type is preserved.
