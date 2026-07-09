# run: environment(make.lacework) <- asNamespace("secr”) after modifying fn

make.lacework <- function (region, spacing = c(100, 20), times = NULL, 
                           origin = NULL, rotate = 0, 
                           radius = NULL, nintersections = NULL,
                           ndetectors = NULL, ndetectors_tol = 0,
                           detector = "multi", keep.design = TRUE) {
  spacing <- as.numeric(spacing)
  if (length(spacing) == 1) {
    if (is.null(times)) stop("make.lacework requires times if spacing length 1")
    b <- spacing[1]
    a <- b * times
  }
  else {
    if (length(spacing) != 2) {
      stop("invalid spacing 2-vector in make.lacework")
    }
    a <- spacing[1]
    b <- spacing[2]
  }
  whole_number <- function(x) {
    length(x) == 1 && is.finite(x) && x == floor(x)
  }
  if (!is.null(nintersections) &&
      (!whole_number(nintersections) || nintersections < 0)) {
    stop("nintersections must be a non-negative whole number")
  }
  if (!is.null(ndetectors) &&
      (!whole_number(ndetectors) || ndetectors < 1)) {
    stop("ndetectors must be a positive whole number")
  }
  if (!whole_number(ndetectors_tol) || ndetectors_tol < 0) {
    stop("ndetectors_tol must be a non-negative whole number")
  }
  if (!is.null(radius) && !is.null(ndetectors)) {
    stop("radius and ndetectors cannot both be set")
  }
  temporigin <- origin
  fraction <- 1.0  ## suspended code
  region <- boundarytoSF(region)
  if (rotate != 0) {
    centrexy <- sf::st_centroid(sf::st_as_sfc(sf::st_bbox(region)))
    centrexy <- sf::st_coordinates(centrexy)
    region <- sfrotate(region, degrees = -rotate, centrexy = centrexy, usecentroid = FALSE)
  }    
  
  bbox <- matrix(sf::st_bbox(region), ncol = 2)                   ## after rotation
  if (is.null(origin)) {
    origin <- bbox[,1]
    origin <- origin - runif(2)*a
  }

  make_design <- function(region) {
    bbox <- matrix(sf::st_bbox(region), ncol = 2)
    rx <- diff(bbox[1,])
    ry <- diff(bbox[2,])
    n1 <- ceiling((rx + a) / a) + 2
    n2 <- ceiling((ry + a) / b) + 2
    gridx <- make.grid(nx=n1, ny=n2, spacex = a, spacey = b,
                       originxy = origin, detector = detector, ID = 'numxb')
    metax <- data.frame(
      family = "x",
      line = round((gridx$x - origin[1]) / a),
      x = gridx$x,
      y = gridx$y
    )
    n1 <- ceiling((ry + a) / a) + 2
    n2 <- ceiling((rx + a) / b) + 2
    gridy <- make.grid(ny=n1, nx=n2, spacey = a, spacex = b,
                       originxy = origin, detector = detector, ID = 'numyb')
    metay <- data.frame(
      family = "y",
      line = round((gridy$y - origin[2]) / a),
      x = gridy$x,
      y = gridy$y
    )
    if (fraction < 1) {
      OKx <- ((gridx$y-origin[2]) %% a) < fraction * a
      OKy <- ((gridy$x-origin[1]) %% a) < fraction * a
      gridx <- subset(gridx,OKx)
      gridy <- subset(gridy,OKy)
      metax <- metax[OKx, , drop = FALSE]
      metay <- metay[OKy, , drop = FALSE]
    }
    grid <- rbind(gridx, gridy, renumber = FALSE, checkdetector = FALSE)
    meta <- rbind(metax, metay)
    dupl <- duplicated(round(grid,5))
    crossings <- subset(grid, dupl)
    crossingmeta <- data.frame(
      xline = round((crossings$x - origin[1]) / a),
      yline = round((crossings$y - origin[2]) / a),
      x = crossings$x,
      y = crossings$y
    )
    grid <- subset(grid, !dupl)
    meta <- meta[!dupl, , drop = FALSE]
    OKgrid <- pointsInPolygon(grid, region)
    OKcrossings <- pointsInPolygon(crossings, region)
    grid <- subset(grid, OKgrid)
    meta <- meta[OKgrid, , drop = FALSE]
    crossings <- subset(crossings, OKcrossings)
    crossingmeta <- crossingmeta[OKcrossings, , drop = FALSE]
    list(grid = grid, crossings = crossings, meta = meta,
         crossingmeta = crossingmeta)
  }

  choose_crossings <- function(crossingmeta, n) {
    if (n == 0) {
      return(list(xline = numeric(0), yline = numeric(0),
                  dimensions = c(nx = 0, ny = 0)))
    }
    xlines <- sort(unique(crossingmeta$xline))
    ylines <- sort(unique(crossingmeta$yline))
    centre <- c(mean(range(crossingmeta$xline)), mean(range(crossingmeta$yline)))
    present <- matrix(FALSE, nrow = length(xlines), ncol = length(ylines))
    present[cbind(match(crossingmeta$xline, xlines),
                  match(crossingmeta$yline, ylines))] <- TRUE
    ps <- rbind(0, cbind(0, apply(present, 2, cumsum)))
    ps <- t(apply(ps, 1, cumsum))
    rectsum <- function(ix, iy, nx, ny) {
      x2 <- ix + nx - 1
      y2 <- iy + ny - 1
      ps[x2 + 1, y2 + 1] - ps[ix, y2 + 1] - ps[x2 + 1, iy] + ps[ix, iy]
    }
    best <- NULL
    bestscore <- NULL
    better_score <- function(score, bestscore) {
      if (is.null(bestscore)) return(TRUE)
      firstdiff <- which(score != bestscore)[1]
      !is.na(firstdiff) && score[firstdiff] < bestscore[firstdiff]
    }

    for (nx in seq_along(xlines)) {
      for (ny in seq_along(ylines)) {
        product <- nx * ny
        if (product < n) next
        aspect <- max(nx, ny) / min(nx, ny)
        for (ix in seq_len(length(xlines) - nx + 1)) {
          xs <- xlines[ix:(ix + nx - 1)]
          for (iy in seq_len(length(ylines) - ny + 1)) {
            ys <- ylines[iy:(iy + ny - 1)]
            if (rectsum(ix, iy, nx, ny) != product) next
            rectanglecentre <- c(mean(range(xs)), mean(range(ys)))
            centred <- sum((rectanglecentre - centre)^2)
            aspect_ok <- aspect <= 2
            score <- c(!aspect_ok, product, abs(log(nx / ny)), centred)
            if (better_score(score, bestscore)) {
              bestscore <- score
              best <- list(xline = xs, yline = ys,
                           dimensions = c(nx = nx, ny = ny))
            }
          }
        }
      }
    }
    if (is.null(best)) {
      stop("Could not find a rectangular block of intersections")
    }
    best
  }

  subset_design <- function(design, keep_grid, keep_crossings) {
    list(
      grid = subset(design$grid, keep_grid),
      crossings = subset(design$crossings, keep_crossings),
      meta = design$meta[keep_grid, , drop = FALSE],
      crossingmeta = design$crossingmeta[keep_crossings, , drop = FALSE]
    )
  }

  design0 <- make_design(region)
  selected_dimensions <- NULL
  selected_xline <- NULL
  selected_yline <- NULL
  cluster_halfwidth <- NULL
  pre_drop_radius <- NULL
  pre_drop_ndetectors <- NULL
  dropped_ndetectors <- 0
  actual_nintersections <- nrow(design0$crossings)

  if (!is.null(nintersections)) {
    if (nrow(design0$crossings) < nintersections) {
      stop("Requested nintersections ", nintersections,
           " exceeds the maximum available ", nrow(design0$crossings))
    }
    selected <- choose_crossings(design0$crossingmeta, nintersections)
    keep_crossings <- design0$crossingmeta$xline %in% selected$xline &
      design0$crossingmeta$yline %in% selected$yline
    if (sum(keep_crossings) > nintersections) {
      candidate_crossings <- design0$crossingmeta[keep_crossings, , drop = FALSE]
      cropped_centroid <- colMeans(candidate_crossings[, c("x", "y")])
      distance_to_cropped_centroid <-
        (candidate_crossings$x - cropped_centroid[1])^2 +
        (candidate_crossings$y - cropped_centroid[2])^2
      drop_local <- order(distance_to_cropped_centroid, decreasing = TRUE)[
        seq_len(sum(keep_crossings) - nintersections)
      ]
      keep_index <- which(keep_crossings)
      keep_crossings[keep_index[drop_local]] <- FALSE
    }
    keep_grid <- rep(FALSE, nrow(design0$meta))
    if (sum(keep_crossings) > 0) {
      retained_crossings <- design0$crossingmeta[keep_crossings, , drop = FALSE]
      cluster_halfwidth <- a / 2
      for (xline in sort(unique(retained_crossings$xline))) {
        y_on_line <- retained_crossings$y[retained_crossings$xline == xline]
        on_line <- design0$meta$family == "x" & design0$meta$line == xline
        yrange <- range(y_on_line) + c(-cluster_halfwidth, cluster_halfwidth)
        keep_grid[on_line] <- design0$meta$y[on_line] >= yrange[1] &
          design0$meta$y[on_line] <= yrange[2]
      }
      for (yline in sort(unique(retained_crossings$yline))) {
        x_on_line <- retained_crossings$x[retained_crossings$yline == yline]
        on_line <- design0$meta$family == "y" & design0$meta$line == yline
        xrange <- range(x_on_line) + c(-cluster_halfwidth, cluster_halfwidth)
        keep_grid[on_line] <- design0$meta$x[on_line] >= xrange[1] &
          design0$meta$x[on_line] <= xrange[2]
      }
      demoted_crossings <- !keep_crossings & (
        design0$crossingmeta$xline %in% unique(retained_crossings$xline) |
          design0$crossingmeta$yline %in% unique(retained_crossings$yline)
      )
      if (any(demoted_crossings)) {
        for (i in which(demoted_crossings)) {
          cm <- design0$crossingmeta[i, , drop = FALSE]
          xline <- retained_crossings[retained_crossings$xline == cm$xline, , drop = FALSE]
          yline <- retained_crossings[retained_crossings$yline == cm$yline, , drop = FALSE]
          on_x_span <- nrow(xline) > 0 &&
            cm$y >= min(xline$y) - cluster_halfwidth &&
            cm$y <= max(xline$y) + cluster_halfwidth
          on_y_span <- nrow(yline) > 0 &&
            cm$x >= min(yline$x) - cluster_halfwidth &&
            cm$x <= max(yline$x) + cluster_halfwidth
          if (on_x_span || on_y_span) {
            gridrow <- design0$crossings[i, , drop = FALSE]
            rownames(gridrow) <- paste0(rownames(gridrow), ".det")
            metarow <- data.frame(
              family = if (on_x_span) "x" else "y",
              line = if (on_x_span) cm$xline else cm$yline,
              x = cm$x,
              y = cm$y
            )
            design0$grid <- rbind(design0$grid, gridrow,
                                  renumber = FALSE, checkdetector = FALSE)
            design0$meta <- rbind(design0$meta, metarow)
            keep_grid <- c(keep_grid, TRUE)
          }
        }
      }
    }
    design0 <- subset_design(design0, keep_grid, keep_crossings)
    selected_dimensions <- selected$dimensions
    selected_xline <- selected$xline
    selected_yline <- selected$yline
    actual_nintersections <- nrow(design0$crossings)
  }

  if (!is.null(ndetectors)) {
    if (nrow(design0$crossings) == 0) {
      stop("Cannot search radius because there are no intersections")
    }
    if (ndetectors < nrow(design0$crossings)) {
      stop("ndetectors must be at least the number of intersections")
    }
    d <- distancetotrap(design0$grid, design0$crossings)
    candidates <- sort(unique(c(0, d)))
    counts <- vapply(candidates, function(r) sum(d <= r), integer(1))
    eligible <- which(counts >= ndetectors)
    if (length(eligible) < 1) {
      stop("Could not achieve ndetectors ", ndetectors,
           "; maximum was ", max(counts), " detectors")
    }
    best <- eligible[which.min(counts[eligible] - ndetectors)]
    radius <- candidates[best]
    pre_drop_radius <- radius
    pre_drop_ndetectors <- counts[best]
  }

  grid <- design0$grid
  crossings <- design0$crossings
  if (!is.null(radius)) {
    grid <- subset(grid, distancetotrap(grid, crossings)<=radius)
  }
  if (!is.null(ndetectors) && nrow(grid) > ndetectors) {
    pre_drop_xy <- as.matrix(grid[, c("x", "y")])
    centroid <- colMeans(pre_drop_xy)
    cross_xy <- as.matrix(crossings[, c("x", "y")])
    distance_to_crossing <- as.matrix(dist(rbind(pre_drop_xy, cross_xy)))
    distance_to_crossing <- distance_to_crossing[
      seq_len(nrow(pre_drop_xy)),
      nrow(pre_drop_xy) + seq_len(nrow(cross_xy)),
      drop = FALSE
    ]
    cluster <- max.col(-distance_to_crossing, ties.method = "first")
    distance_to_centroid <- sqrt((pre_drop_xy[,1] - centroid[1])^2 +
                                   (pre_drop_xy[,2] - centroid[2])^2)
    droppable <- apply(round(distance_to_crossing, 5), 1, min) > 0
    keep <- rep(TRUE, nrow(grid))
    todrop <- nrow(grid) - ndetectors
    while (todrop > 0) {
      candidate <- vapply(seq_len(nrow(crossings)), function(k) {
        eligible <- which(keep & droppable & cluster == k)
        if (length(eligible) < 1) return(NA_integer_)
        eligible[which.max(distance_to_centroid[eligible])]
      }, integer(1))
      candidate <- candidate[!is.na(candidate)]
      if (length(candidate) < 1) {
        stop("Could not drop enough detectors without removing intersections")
      }
      candidate <- candidate[order(distance_to_centroid[candidate], decreasing = TRUE)]
      thisdrop <- head(candidate, todrop)
      keep[thisdrop] <- FALSE
      todrop <- todrop - length(thisdrop)
    }
    dropped_ndetectors <- sum(!keep)
    grid <- subset(grid, keep)
  }
  rownames(grid) <- sapply(lapply(strsplit(rownames(grid), ".", TRUE), rev), paste, collapse='.')
  if (rotate != 0) {
    grid <- secr::rotate(grid, degrees = rotate, centrexy = centrexy)
    crossings <- secr::rotate(crossings, degrees = rotate, centrexy = centrexy)
  }
  attr(grid, 'crossings') <- as.matrix(crossings)
  if (keep.design) {
    design <- list (
      'function' = 'make.lacework',
      region = region, 
      spacing = spacing,
      origin = temporigin, 
      rotate = rotate, 
      radius = radius,
      nintersections = nintersections,
      actual_nintersections = actual_nintersections,
      selected_dimensions = selected_dimensions,
      selected_xline = selected_xline,
      selected_yline = selected_yline,
      cluster_halfwidth = cluster_halfwidth,
      cluster_cell_rule = if (!is.null(cluster_halfwidth)) "nearest_selected_intersection" else NULL,
      ndetectors = ndetectors,
      ndetectors_tol = ndetectors_tol,
      pre_drop_radius = pre_drop_radius,
      pre_drop_ndetectors = pre_drop_ndetectors,
      dropped_ndetectors = dropped_ndetectors,
      ndetectors_target_met = is.null(ndetectors) || nrow(grid) == ndetectors,
      detector = detector)
    attr(grid, 'design') <- design
  }
  
  grid
}
