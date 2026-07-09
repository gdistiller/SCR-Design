##############################################################################
## Two-group optimal spacing (numerical search over absolute spacings)
## Possible objectives: for two groups A and B
## "sum_min" = maximize min(En_A + En_B, Er_A + Er_B)
## "all_min" = maximize min(En_A, EnB, Er_A, Er_B)
## "min_mean_CV" = minimize mean(CV_A, CV_B)
## "min_max_CV" = minimize max(CV_A, CV_B)
## "max_mean_En2" = maximise mean(En2_A, En2_B) where En2 = n animals on >1 trap
## "max_min_En2" = maximise min(En2_A, En2_B)
##############################################################################

optimalSpacing2G <- function(
    D1, D2,
    traps0,
    detectpar1, detectpar2,
    noccasions,
    nrepeats = 1,
    detectfn1 = "HHN",
    detectfn2 = "HHN",
    xsigma = 4,
    spacing_m = seq(200, 4000, 200),
    criterion = c(
      "sum_min", "all_min", "min_mean_CV", "min_max_CV",
      "max_mean_En2", "max_min_En2"
    ),
    CF = 1.0,
    ...
) {

  valid_criteria <- c(
    "sum_min", "all_min", "min_mean_CV", "min_max_CV",
    "max_mean_En2", "max_min_En2"
  )
  if (missing(criterion)) {
    criterion <- valid_criteria[1]
  } else if (length(criterion) != 1 || !criterion %in% valid_criteria) {
    stop(
      "'criterion' must be one of: ",
      paste(sQuote(valid_criteria), collapse = ", "),
      call. = FALSE
    )
  }

  detectfn1 <- secr:::secr_valid.detectfn(detectfn1, valid = c(0, 1, 2, 14:19))
  detectfn2 <- secr:::secr_valid.detectfn(detectfn2, valid = c(0, 1, 2, 14:19))
  dfc1 <- secrdesign:::dfcast(detectfn1, detectpar1)
  dfc2 <- secrdesign:::dfcast(detectfn2, detectpar2)
  detectfn1 <- dfc1$detectfn
  detectfn2 <- dfc2$detectfn
  detectpar1 <- dfc1$detectpar
  detectpar2 <- dfc2$detectpar

  base_spacing <- spacing(traps0)

  args <- list(...)
  defaultmaskargs <- list(nx = 64, type = "trapbuffer")
  dotsmask <- args[names(args) %in% c("nx", "type", "poly", "poly.habitat")]
  maskargs_base <- secrdesign:::replacedefaults(defaultmaskargs, dotsmask)

  if (any(detector(traps0) == "single")) {
    warning("treating single-catch traps as multi-catch", call. = FALSE)
    detector(traps0) <- "multi"
  }

  values <- lapply(
    spacing_m,
    getCrit2G_abs,
    traps0,
    base_spacing,
    xsigma,
    detectpar1,
    detectpar2,
    D1,
    D2,
    noccasions,
    nrepeats,
    detectfn1,
    detectfn2,
    maskargs_base,
    CF,
    criterion
  )
  values <- do.call(rbind, values)
  values <- as.data.frame(values)
  names(values) <- c(
    "spacing", "En1", "En2", "Er1", "Er2",
    "CV1", "CV2", "En2plus1", "En2plus2", "crit"
  )

  opt <- if (!is.na(CF) && nrow(values) > 0) {
    ok <- is.finite(values$spacing) & is.finite(values$crit)
    if (any(ok)) {
      finite_values <- values[ok, , drop = FALSE]
      optrow <- finite_values[which.max(finite_values$crit), ]
      list(spacing = optrow$spacing, objective = optrow$crit)
    } else {
      list(spacing = NA, objective = NA)
    }
  } else {
    list(spacing = NA, objective = NA)
  }

  out <- list(
    values = values,
    optimum.spacing = opt$spacing,
    optimum.crit = opt$objective,
    maximum.crit = opt$objective,
    criterion = criterion,
    traps.base = traps0,
    detectpar1 = detectpar1,
    detectpar2 = detectpar2
  )
  class(out) <- c("optimalSpacing2G", "list")
  out
}

##############################################################################

getCrit2G_abs <- function(
    S,
    traps0,
    base_spacing,
    xsigma,
    detectpar1,
    detectpar2,
    D1,
    D2,
    noccasions,
    nrepeats,
    detectfn1,
    detectfn2,
    maskargs_base,
    CF,
    criterion
) {

  scalefac <- S / base_spacing
  trapS <- scale_traps(traps0, scalefac)

  maskargs <- maskargs_base
  maskargs$traps <- trapS
  maskargs$buffer <- xsigma * max(detectpar1$sigma, detectpar2$sigma) + S
  mask <- do.call(make.mask, maskargs)

  En1 <- Enrm(D1, trapS, mask, detectpar1, noccasions, detectfn1) * nrepeats
  En2 <- Enrm(D2, trapS, mask, detectpar2, noccasions, detectfn2) * nrepeats
  Er1 <- En1[2]
  Er2 <- En2[2]
  En1 <- En1[1]
  En2 <- En2[1]
  CV1 <- 1 / sqrt(min(En1, Er1))
  CV2 <- 1 / sqrt(min(En2, Er2))
  En2plus1 <- secrdesign::En2(D1, trapS, mask, detectpar1, noccasions, detectfn1)[2] * nrepeats
  En2plus2 <- secrdesign::En2(D2, trapS, mask, detectpar2, noccasions, detectfn2)[2] * nrepeats

  critval <- switch(
    criterion,
    sum_min = min(En1 + En2, Er1 + Er2),
    all_min = min(En1, En2, Er1, Er2),
    min_mean_CV = -mean(c(CV1, CV2)),
    min_max_CV = -max(CV1, CV2),
    max_mean_En2 = En2plus1 + En2plus2,
    max_min_En2 = min(En2plus1, En2plus2)
  )

  c(S, En1, En2, Er1, Er2, CV1, CV2, En2plus1, En2plus2, critval * CF)
}

scale_traps <- function(traps, scalefac) {
  tr <- traps
  tr$x <- tr$x * scalefac
  tr$y <- tr$y * scalefac
  for (nm in c("spacing", "spacex", "spacey")) {
    if (!is.null(attr(tr, nm))) {
      attr(tr, nm) <- attr(tr, nm) * scalefac
    }
  }
  tr
}

##############################################################################
