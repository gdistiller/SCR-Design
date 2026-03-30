# function to crop lacework design to N points 
#
# start with all points
# repeatedly remove one point at a time
# at each step, define the current centroid of the retained points
# compute each point’s distance from that centroid
# remove the farthest point
# repeat until exactly N points remain

crop_to_n <- function(x, y, N, xy_ratio = 1,
                                        tie_break = c("boundary", "random")) {
  stopifnot(length(x) == length(y))
  stopifnot(N >= 1, N <= length(x))
  stopifnot(xy_ratio > 0)
  
  tie_break <- match.arg(tie_break)
  
  pts <- data.frame(
    id = seq_along(x),
    x = x,
    y = y
  )
  
  keep <- rep(TRUE, nrow(pts))
  
  # Anisotropic squared distance to centroid.
  # If xy_ratio = width / height:
  #   larger xy_ratio means x-extent is allowed to be larger,
  #   so x-deviations should be penalized a bit less.
  anisotropic_dist2 <- function(x, y, cx, cy, xy_ratio) {
    ((x - cx) / xy_ratio)^2 + (y - cy)^2
  }
  
  while (sum(keep) > N) {
    cur <- pts[keep, , drop = FALSE]
    
    cx <- mean(cur$x)
    cy <- mean(cur$y)
    
    d2 <- anisotropic_dist2(cur$x, cur$y, cx, cy, xy_ratio)
    
    # candidate(s): farthest from current centroid
    max_d2 <- max(d2)
    cand_local <- which(d2 == max_d2)
    
    # tie-break if multiple points are equally far
    if (length(cand_local) > 1) {
      cand <- cur[cand_local, , drop = FALSE]
      
      if (tie_break == "boundary") {
        # Prefer removing the point closest to the current boundary
        # (more "outside-looking" in a box sense)
        xmin <- min(cur$x); xmax <- max(cur$x)
        ymin <- min(cur$y); ymax <- max(cur$y)
        
        margin <- pmin(
          cand$x - xmin,
          xmax - cand$x,
          cand$y - ymin,
          ymax - cand$y
        )
        
        remove_local <- cand_local[which.min(margin)]
      } else {
        remove_local <- sample(cand_local, 1)
      }
    } else {
      remove_local <- cand_local
    }
    
    remove_id <- cur$id[remove_local]
    keep[pts$id == remove_id] <- FALSE
  }
  
  kept <- pts[keep, , drop = FALSE]
  
  list(
    points = kept[, c("x","y")],
    n_kept = nrow(kept),
    centroid = c(x = mean(kept$x), y = mean(kept$y)),
    bbox = c(
      xmin = min(kept$x),
      xmax = max(kept$x),
      ymin = min(kept$y),
      ymax = max(kept$y)
    ),
    removed_ids = pts$id[!keep]
  )
}