# =============================================================================
# Shiny App: Management Zone Delimitation
# Fuzzy C-Means Clustering with FPI/NCE optimal k selection
# =============================================================================
# Requirements:
#   install.packages(c("shiny","shinythemes","shinyjs","e1071","vegan","gstat",
#                      "sf","terra","sp","ggplot2","ggrepel","patchwork","plotly",
#                      "DT","readr","stringr","shinycssloaders"),
#                    repos = "https://cran.r-project.org")
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(shinyjs)
  library(e1071)
  library(vegan)
  library(gstat)
  library(sf)
  library(terra)
  library(sp)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(plotly)
  library(DT)
  library(tidyr)
  library(shinycssloaders)
})

# ── Increase file upload limit (default is 5 MB — too small for raster files)
#    1 GB = 1e9 bytes | Adjust as needed for your largest raster
options(shiny.maxRequestSize = 1e9)   # 1 GB limit — covers most GeoTIFFs

# ── Diagnostic logger ────────────────────────────────────────────────────────
# Writes a timestamped trace of the pipeline to mz_debug.log next to this file.
# This is what lets us see WHERE Part 4 stalls on a given machine: run the app,
# reproduce the issue, then read mz_debug.log.  Set MZ_DEBUG=FALSE to silence.
MZ_VERSION <- "2026-06-03-resize-fix"
MZ_DEBUG   <- TRUE
MZ_LOGFILE <- file.path(getwd(), "mz_debug.log")
mzlog <- function(...) {
  if (!isTRUE(MZ_DEBUG)) return(invisible())
  line <- sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%OS2"),
                  paste0(..., collapse = ""))
  try(cat(line, file = MZ_LOGFILE, append = TRUE), silent = TRUE)
}

# ─────────────────────────────────────────────────────────────────────────────
# VALIDATION INDEX FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

# Fuzziness Performance Index (FPI) — McBratney & Moore (1985)
# FPI = 1 - (k*PC - 1) / (k - 1), where PC = sum(U^2) / n
# Optimal k = argmin(FPI) — minimum fuzzy overlap between zones
fpi_index <- function(U) {
  k  <- ncol(U)
  pc <- sum(U^2) / nrow(U)
  1 - ((k * pc) - 1) / (k - 1)
}

# Normalized Classification Entropy (NCE)
# NCE = -sum(u_ik * log(u_ik)) / (n * log(k))
# Optimal k = argmax(NCE) — maximum classification entropy
nce_index <- function(U) {
  k <- ncol(U)
  n <- nrow(U)
  U_safe <- pmax(U, .Machine$double.xmin)
  -sum(U_safe * log(U_safe)) / (n * log(k))
}

# Xie-Beni Index (XB) — lower is better
xie_beni <- function(X, U, centroids, m = 2) {
  n   <- nrow(X)
  k   <- ncol(U)
  d   <- ncol(X)
  Um  <- U^m
  dist_sq <- vapply(seq_len(k), function(j) {
    rowSums((X - matrix(centroids[j, ], nrow = n, ncol = d, byrow = TRUE))^2)
  }, numeric(n))
  Jm <- sum(Um * dist_sq)
  cdist_sq <- as.matrix(dist(centroids))^2
  diag(cdist_sq) <- NA
  min_dist_sq <- min(cdist_sq, na.rm = TRUE)
  Jm / (n * min_dist_sq)
}

# Partition Entropy (PE) — lower is better
pe_index <- function(U, m = 2) {
  n <- nrow(U)
  -sum(pmax(U, .Machine$double.xmin) * log(pmax(U, .Machine$double.xmin))) / n
}

# Fukuyama-Sugeno (FS) — lower is better
fs_index <- function(X, U, centroids, m = 2) {
  n   <- nrow(X)
  k   <- ncol(U)
  d   <- ncol(X)
  Um  <- U^m
  dist_sq <- vapply(seq_len(k), function(j) {
    rowSums((X - matrix(centroids[j, ], nrow = n, ncol = d, byrow = TRUE))^2)
  }, numeric(n))
  term1 <- sum(Um * dist_sq)
  w  <- colSums(Um)
  cbar <- colSums(t(centroids) * w) / sum(w)
  dist_cbar_sq <- rowSums((centroids - matrix(cbar, nrow = k, ncol = d, byrow = TRUE))^2)
  term1 - sum(w * dist_cbar_sq)
}

# Min-max normalization
normalize <- function(x) {
  rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (rng == 0) return(rep(0.5, length(x)))
  (x - min(x, na.rm = TRUE)) / rng
}

# ─────────────────────────────────────────────────────────────────────────────
# CORE ANALYSIS FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

runValidation <- function(pca_scores, k_range = 2:6, m = 2) {
  results <- lapply(k_range, function(k) {
    set.seed(42)
    fcm_out <- cmeans(pca_scores, centers = k, m = m,
                      method = "cmeans", iter.max = 300, verbose = FALSE)
    U <- fcm_out$membership
    C <- fcm_out$centers

    data.frame(
      k   = k,
      XB  = xie_beni(pca_scores, U, C, m = m),
      PE  = pe_index(U, m = m),
      FS  = fs_index(pca_scores, U, C, m = m),
      FPI = fpi_index(U),
      NCE = nce_index(U)
    )
  })
  validation <- do.call(rbind, results)
  # Add the *_rank columns and consensus_rank, matching what the standalone
  # script's evaluate_k() writes. Without these, the Validation CSV that the
  # Export button writes is missing the columns the Quarto report expects
  # (consensus_rank, used in section 5's validity table).
  validation$XB_rank <- rank(validation$XB, ties.method = "min")
  validation$FPI_rank <- rank(validation$FPI, ties.method = "min")
  validation$NCE_rank <- rank(validation$NCE, ties.method = "min")
  validation$consensus_rank <- validation$XB_rank +
    validation$FPI_rank + validation$NCE_rank
  validation
}

selectOptimalK <- function(df, method = "consensus") {
  # All indices are argmin.  "consensus" (and "auto" for backward
  # compatibility) sums the ranks of XB + FPI + NCE and picks the
  # k with the lowest consensus_rank — matches the reference.
  consensus <- {
    df$XB_rank  <- rank(df$XB,  ties.method = "min")
    df$FPI_rank <- rank(df$FPI, ties.method = "min")
    df$NCE_rank <- rank(df$NCE, ties.method = "min")
    df$consensus_rank <- df$XB_rank + df$FPI_rank + df$NCE_rank
    df$k[which.min(df$consensus_rank)]
  }
  switch(method,
    XB        = df$k[which.min(df$XB)],
    FPI       = df$k[which.min(df$FPI)],
    NCE       = df$k[which.min(df$NCE)],
    PE        = df$k[which.min(df$PE)],
    FS        = df$k[which.min(df$FS)],
    auto = , consensus = consensus
  )
}

runFCM <- function(pca_scores, k, m = 2) {
  set.seed(42)
  fcm_out <- cmeans(pca_scores, centers = k, m = m,
                    method = "cmeans", iter.max = 300, verbose = FALSE)
  list(
    membership  = fcm_out$membership,
    centroids   = fcm_out$centers,
    cluster_id  = max.col(fcm_out$membership)
  )
}

loadBoundary <- function(path) {
  ext <- tools::file_ext(path)
  if (ext == "gpkg") {
    st_read(path, quiet = TRUE)
  } else if (ext %in% c("shp", "json", "geojson")) {
    st_read(path, quiet = TRUE)
  } else {
    stop("Unsupported boundary format. Use .gpkg, .shp, or .geojson")
  }
}

loadRaster <- function(path) {
  rast(path)
}

loadPointsCSV <- function(path) {
  read.csv(path)
}

buildFeatureMatrix <- function(obs_raw, soil_masked, soil_vars, pca_thresh = 0.80,
                               point_crs = 4326) {
  # obs_raw  : data frame with X, Y columns (already renamed from col_lon/col_lat)
  # soil_masked: masked SpatRaster from terra
  # soil_vars  : character vector of selected variable names from CSV
  # pca_thresh : cumulative variance threshold (default 0.80 = 80%)
  # point_crs  : CRS of the observation points (EPSG code or CRS string;
  #               default 4326 for lat/lon; use raster CRS for UTM/projected data)

  # Convert observation points to sf
  obs_sf <- st_as_sf(obs_raw, coords = c("X", "Y"), crs = point_crs)
  obs_sf <- st_transform(obs_sf, crs(soil_masked))

  # Extract raster values at point locations.
  # NOTE: must be terra::extract — tidyr (loaded later) masks `extract` with an
  # S3 generic that has no SpatRaster method ("no applicable method" error).
  extracted <- terra::extract(soil_masked, obs_sf, ID = FALSE)

  # Use CSV soil variables (NOT raster layer names — those may differ)
  df_vars <- obs_raw[, soil_vars, drop = FALSE]

  # Overwrite with raster-extracted values where layer names match.
  # Direct name-based lookup avoids cbind column-name collisions that occur
  # when CSV columns and raster layers share names.
  for (v in names(soil_masked)) {
    if (v %in% names(df_vars) && v %in% names(extracted)) {
      df_vars[[v]] <- extracted[[v]]
    }
  }

  # Median imputation
  for (j in seq_along(soil_vars)) {
    if (any(is.na(df_vars[[j]]))) {
      med <- median(df_vars[[j]], na.rm = TRUE)
      if (is.na(med)) med <- 0.5
      df_vars[[j]][is.na(df_vars[[j]])] <- med
    }
  }

  df_norm <- as.data.frame(lapply(df_vars, normalize))
  stopifnot("NAs remain after normalisation" = !any(is.na(df_norm)))

  # PCA
  pca_raw  <- prcomp(df_norm, center = FALSE, scale. = FALSE)
  eig      <- pca_raw$sdev^2
  prop_var <- eig / sum(eig)
  cumvar   <- cumsum(prop_var)
  n_pc     <- which(cumvar >= pca_thresh)[1]
  if (is.na(n_pc)) n_pc <- length(eig)

  list(
    df_vars    = df_vars,
    df_norm    = df_norm,
    pca_scores = pca_raw$x[, seq_len(n_pc), drop = FALSE],
    pca_raw    = pca_raw,
    prop_var   = prop_var,
    cumvar     = cumvar,
    n_pc       = n_pc,
    obs_sf     = obs_sf   # return sf for later use in zone mapping
  )
}

buildZoneMap <- function(obs_sf, membership, boundary_terra, soil_masked, k,
                         method = "kriging", max_grid_cells = NULL) {
  # Mirrors delimit_management_zones.R::interpolate_membership.
  #
  # `max_grid_cells` is a performance switch:
  #   - NULL (default): kriging runs on ALL valid cells of the masked
  #     raster.  Bit-identical to the reference, but slow on large
  #     extents (~38s on the 1.5M-cell bundled raster).
  #   - <integer>: aggregate the raster to ~max_grid_cells, run kriging
  #     on the smaller grid, then bilinear-resample back to the full
  #     resolution.  ~10x faster, results equivalent within rounding.
  # The Shiny app passes 40000 so the UI doesn't freeze; the test
  # harness leaves it NULL to assert bit-identical equivalence.
  ref_rast    <- soil_masked[[1]]
  valid_cells <- which(!is.na(values(ref_rast, mat = FALSE)))

  if (length(valid_cells) == 0) {
    stop("No valid (non-NA) cells in the masked raster. Check boundary extent.")
  }

  if (is.null(max_grid_cells) || length(valid_cells) <= max_grid_cells) {
    target_rast  <- ref_rast
    target_valid <- valid_cells
    do_resample  <- FALSE
  } else {
    fact         <- max(1L, ceiling(sqrt(length(valid_cells) / max_grid_cells)))
    target_rast  <- aggregate(ref_rast, fact = fact, fun = mean, na.rm = TRUE)
    target_rast  <- mask(target_rast, boundary_terra)
    target_valid <- which(!is.na(values(target_rast, mat = FALSE)))
    do_resample  <- TRUE
  }

  grid_xy <- as.data.frame(xyFromCell(target_rast, target_valid))
  names(grid_xy) <- c("X", "Y")
  coordinates(grid_xy) <- ~ X + Y
  # PROJ4 string, not WKT — safe with modern sp/terra
  proj4string(grid_xy) <- CRS(crs(target_rast, proj = TRUE))

  layers <- vector("list", k)

  for (j in seq_len(k)) {
    # Attach the membership as a named column on the Spatial object so
    # gstat can resolve the formula LHS — required for both krige()
    # and idw() in this gstat build.
    pts <- obs_sf
    pts$membership_value <- membership[, j]
    pts_sp <- as(pts["membership_value"], "Spatial")

    predicted <- tryCatch({
      if (method == "kriging") {
        vg <- variogram(membership_value ~ 1, pts_sp)
        initial_model <- vgm(
          psill = max(stats::var(pts$membership_value), 1e-6),
          model = "Sph",
          range = max(diff(range(grid_xy@coords[, 1])),
                      diff(range(grid_xy@coords[, 2]))) / 3,
          nugget = 0
        )
        vg_fit <- fit.variogram(vg, initial_model)
        krige(membership_value ~ 1, pts_sp, grid_xy,
              model = vg_fit, debug.level = 0)$var1.pred
      } else {
        idw(membership_value ~ 1, pts_sp, grid_xy,
            idp = 2, debug.level = 0)$var1.pred
      }
    }, error = function(e) {
      message("    Falling back to IDW for z", j, ": ", conditionMessage(e))
      idw(membership_value ~ 1, pts_sp, grid_xy,
          idp = 2, debug.level = 0)$var1.pred
    })

    vals <- rep(NA_real_, ncell(target_rast))
    vals[target_valid] <- pmin(pmax(predicted, 0), 1)

    r <- setValues(target_rast, vals)
    r <- mask(r, boundary_terra)
    if (do_resample) {
      # Upsample to the original raster's resolution.  mask() again
      # because bilinear can leak values just outside the boundary.
      r <- resample(r, ref_rast, method = "bilinear")
      r <- mask(r, boundary_terra)
    }
    names(r) <- paste0("Memb_z", j)
    layers[[j]] <- r
  }

  rast(layers)
}

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────

ui <- fluidPage(

  useShinyjs(),

  tags$head(
    tags$link(rel = "stylesheet", type = "text/css",
              href = "https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap"),
    tags$style(HTML("
      * { font-family: 'Inter', sans-serif; }
      body { background: #f0f2f5; margin: 0; }

      /* ── Sidebar ── */
      .sidebar {
        position: fixed; top: 0; left: 0; width: 260px; height: 100vh;
        background: linear-gradient(180deg, #1a252f 0%, #2c3e50 100%);
        color: white; z-index: 1000; overflow-y: auto;
        box-shadow: 2px 0 10px rgba(0,0,0,0.15);
      }
      .sidebar-logo {
        padding: 24px 20px; border-bottom: 1px solid rgba(255,255,255,0.1);
        font-size: 18px; font-weight: 700; letter-spacing: 0.5px;
      }
      .sidebar-logo span { color: #3498db; }
      .sidebar-nav { padding: 10px 0; }
      .nav-item {
        padding: 12px 20px; cursor: pointer; transition: all 0.2s;
        border-left: 3px solid transparent; font-size: 14px; font-weight: 400;
        display: flex; align-items: center; gap: 10px;
      }
      .nav-item:hover { background: rgba(255,255,255,0.08); }
      .nav-item.active {
        background: rgba(52,152,219,0.2); border-left-color: #3498db;
        font-weight: 600;
      }
      .nav-item i { font-size: 16px; width: 20px; text-align: center; }

      /* ── Main content ── */
      .main-content { margin-left: 260px; min-height: 100vh; }

      /* ── Header bar ── */
      .top-bar {
        background: white; padding: 16px 28px; display: flex;
        justify-content: space-between; align-items: center;
        box-shadow: 0 2px 8px rgba(0,0,0,0.06); position: sticky; top: 0; z-index: 100;
      }
      .top-bar h1 { margin: 0; font-size: 22px; font-weight: 700; color: #1a252f; }
      .top-bar-meta { font-size: 13px; color: #95a5a6; }

      /* ── Content area ── */
      .content-area { padding: 24px 28px; }

      /* ── Cards ── */
      .card {
        background: white; border-radius: 12px; box-shadow: 0 2px 12px rgba(0,0,0,0.06);
        margin-bottom: 20px; overflow: hidden;
      }
      .card-header {
        padding: 16px 20px; border-bottom: 1px solid #eee;
        display: flex; justify-content: space-between; align-items: center;
        background: #fafbfc;
      }
      .card-title { font-size: 15px; font-weight: 600; color: #1a252f; margin: 0; }
      .card-body { padding: 20px; }

      /* ── Section divider ── */
      .section-title {
        font-size: 13px; font-weight: 700; text-transform: uppercase;
        letter-spacing: 1px; color: #3498db; margin-bottom: 12px;
      }

      /* ── File input boxes ── */
      .file-box {
        background: #f8f9fa; border: 2px dashed #dde1e6; border-radius: 10px;
        padding: 20px; text-align: center; transition: all 0.2s; cursor: pointer;
        margin-bottom: 12px;
      }
      .file-box:hover { border-color: #3498db; background: #eaf4fb; }
      .file-box.loaded { border-color: #27ae60; background: #eafaf1; }
      .file-box-icon { font-size: 28px; color: #bdc3c7; margin-bottom: 8px; }
      .file-box.loaded .file-box-icon { color: #27ae60; }
      .file-box-label { font-size: 14px; color: #7f8c8d; font-weight: 600; }
      .file-box-status { font-size: 12px; color: #95a5a6; margin-top: 4px; }

      /* ── Settings inputs ── */
      .setting-group { margin-bottom: 16px; }
      .setting-label { font-size: 12px; font-weight: 600; color: #5d6d7e; margin-bottom: 6px; text-transform: uppercase; letter-spacing: 0.5px; }
      .form-control-sm { border-radius: 8px; font-size: 13px; }
      select.form-control { border-radius: 8px; font-size: 13px; }

      /* ── Buttons ── */
      .btn-primary-custom {
        background: linear-gradient(135deg, #3498db, #2980b9); color: white;
        border: none; border-radius: 8px; padding: 10px 20px; font-size: 14px;
        font-weight: 600; width: 100%; cursor: pointer; transition: all 0.2s;
        box-shadow: 0 4px 12px rgba(52,152,219,0.3);
      }
      .btn-primary-custom:hover { transform: translateY(-1px); box-shadow: 0 6px 16px rgba(52,152,219,0.4); }
      .btn-primary-custom:disabled { opacity: 0.6; cursor: not-allowed; transform: none; }

      .btn-secondary-custom {
        background: white; color: #3498db; border: 2px solid #3498db;
        border-radius: 8px; padding: 10px 20px; font-size: 14px; font-weight: 600;
        width: 100%; cursor: pointer; transition: all 0.2s;
      }
      .btn-secondary-custom:hover { background: #eaf4fb; }

      /* ── Metric cards ── */
      .metric-row { display: flex; gap: 12px; margin-bottom: 16px; }
      .metric-card {
        flex: 1; background: linear-gradient(135deg, #667eea, #764ba2);
        color: white; padding: 16px; border-radius: 12px; text-align: center;
      }
      .metric-card.green  { background: linear-gradient(135deg, #11998e, #38ef7d); }
      .metric-card.orange { background: linear-gradient(135deg, #f093fb, #f5576c); }
      .metric-card.amber  { background: linear-gradient(135deg, #f5af19, #f12711); }
      .metric-value { font-size: 28px; font-weight: 700; }
      .metric-label { font-size: 11px; opacity: 0.85; margin-top: 4px; text-transform: uppercase; letter-spacing: 0.5px; }

      /* ── Alert / info boxes ── */
      .info-msg {
        background: #eaf4fb; border-left: 4px solid #3498db;
        padding: 10px 14px; border-radius: 6px; font-size: 13px;
        color: #2c3e50; margin-bottom: 12px;
      }
      .success-msg {
        background: #eafaf1; border-left: 4px solid #27ae60;
        padding: 10px 14px; border-radius: 6px; font-size: 13px;
        color: #1e8449; margin-bottom: 12px;
      }
      .warning-msg {
        background: #fef9e7; border-left: 4px solid #f39c12;
        padding: 10px 14px; border-radius: 6px; font-size: 13px;
        color: #9a7d0a; margin-bottom: 12px;
      }

      /* ── Data table styles ── */
      .dataTables_wrapper { font-size: 13px !important; }
      table.dataTable thead th { background: #f8f9fa; font-weight: 700; color: #1a252f; }
      .DT-table { border: none !important; }
      table.dataTable tbody tr:hover { background: #f0f7ff !important; }

      /* ── Plotly ── */
      .plotly .modebar { display: none !important; }
      .modebar { display: none !important; }

      /* ── Progress bar ── */
      .loading-overlay { text-align: center; padding: 40px; color: #95a5a6; }

      /* ── Membership grid ── */
      .memb-grid { display: grid; gap: 16px; }

      /* ── Help tab ── */
      .help-hero { background: linear-gradient(135deg, #f4f7fb 0%, #eaf4fb 100%);
                   border-left: 4px solid #3498db; }
      .help-hero .card-body { padding: 20px 24px; }
      .help-file-box { background: #f8f9fa; border: 2px dashed #dde1e6;
                       border-radius: 10px; padding: 18px 16px;
                       height: 100%; text-align: center;
                       transition: all 0.2s; }
      .help-file-box:hover { border-color: #3498db; background: #eaf4fb; }
      .help-file-icon { font-size: 32px; color: #3498db; }
      .help-file-ext { color: #5d6d7e; font-size: 12px; font-weight: 600;
                       letter-spacing: 0.3px; }
      .help-figure { width: 100%; height: auto; border-radius: 8px;
                     border: 1px solid #eee; display: block; }
      .help-step-grid { display: grid; grid-template-columns: repeat(3, 1fr);
                        gap: 12px; margin-top: 18px; }
      .help-step { background: #f8f9fa; border-radius: 8px; padding: 10px 12px;
                   display: flex; align-items: center; gap: 10px; }
      .help-step-num { background: #3498db; color: white; width: 28px;
                       height: 28px; border-radius: 50%; display: flex;
                       align-items: center; justify-content: center;
                       font-weight: 700; font-size: 13px; flex-shrink: 0; }
      .help-caption { padding: 8px 4px 0 4px; font-size: 13px;
                      color: #2c3e50; }
      .help-figure-tall { max-width: 520px; width: 100%; height: auto;
                          margin: 0 auto; display: block; }
    "))
  ),

  # ── Sidebar ──────────────────────────────────────────────────────────────
  div(class = "sidebar",
    div(class = "sidebar-logo",
      tags$img(src = "mz_logo.svg", height = "36px", width = "36px",
               style = "vertical-align: middle; margin-right: 10px;",
               alt = "MZ Delimitation"),
      tags$span("MZ Delimitation")
    ),
    div(class = "sidebar-nav",
      div(class = "nav-item active", id = "nav-help",     `data-tab` = "tab-help",
        icon("book"), "User Guide"),
      div(class = "nav-item",        id = "nav-data",     `data-tab` = "tab-data",
        icon("database"), "1. Data Input"),
      div(class = "nav-item",        id = "nav-validate", `data-tab` = "tab-validate",
        icon("chart-line"), "2. Validation"),
      div(class = "nav-item",        id = "nav-cluster",  `data-tab` = "tab-cluster",
        icon("layer-group"), "3. Clustering"),
      div(class = "nav-item",        id = "nav-maps",     `data-tab` = "tab-maps",
        icon("map"), "4. Zone Maps"),
      div(class = "nav-item",        id = "nav-stats",    `data-tab` = "tab-stats",
        icon("table"), "5. Statistics"),
      div(class = "nav-item",        id = "nav-export",   `data-tab` = "tab-export",
        icon("download"), "6. Export")
    )
  ),

  # ── Main content ─────────────────────────────────────────────────────────
  div(class = "main-content",

    # Tab 1: Data Input ────────────────────────────────────────────────────────
    div(id = "tab-data", style = "display:none;",
      div(class = "top-bar",
        div(h1("Data Input")),
        div(class = "top-bar-meta", "Step 1 of 6")
      ),
      div(class = "content-area",
        fluidRow(
          # Left: File inputs
          column(8,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Source Data Files")
              ),
              div(class = "card-body",
                div(class = "section-title", "Study Area Boundary"),
                fileInput("input_boundary", NULL,
                  accept = c(".gpkg", ".shp", ".geojson", ".json"),
                  buttonLabel = icon("folder-open"),
                  placeholder = "Select boundary file (.gpkg / .shp / .geojson)..."),
                uiOutput("boundary_status"),

                div(style = "margin-top: 24px;", class = "section-title", "Soil Properties Raster"),
                fileInput("input_raster", NULL,
                  accept = c(".tif", ".tiff", ".asc", ".bil"),
                  buttonLabel = icon("image"),
                  placeholder = "Select soil raster (.tif / .asc / .bil)..."),
                uiOutput("raster_status"),

                div(style = "margin-top: 24px;", class = "section-title", "Observation Points"),
                fileInput("input_points", NULL,
                  accept = c(".csv", ".xlsx"),
                  buttonLabel = icon("map-marker-alt"),
                  placeholder = "Select observation CSV..."),
                uiOutput("points_status")
              )
            ),

            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Variable Selection"),
                div(style = "font-size:12px; color:#7f8c8d;",
                    "Detected from your CSV after upload")
              ),
              div(class = "card-body",
                uiOutput("soil_vars_ui"),
                uiOutput("var_col_preview"),
                div(class = "info-msg",
                  "Variables will appear after you upload the observation CSV.
                   Select at least 3 numeric columns for clustering.")
              )
            )
          ),

          # Right: Settings
          column(4,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Analysis Settings")
              ),
              div(class = "card-body",
                div(class = "setting-group",
                  tags$label(class = "setting-label", "Coordinate Columns"),
                  textInput("col_lon", NULL, value = "longitude", placeholder = "X / longitude column"),
                  textInput("col_lat", NULL, value = "latitude",  placeholder = "Y / latitude column")
                ),
                div(class = "setting-group",
                  tags$label(class = "setting-label", "Point Coordinate System"),
                  selectInput("point_crs_mode", NULL,
                    choices = c(
                      "Auto-detect" = "auto",
                      "EPSG:4326 — WGS84 (lat/lon)" = "4326",
                      "Same as raster CRS" = "raster"
                    ),
                    selected = "auto"
                  ),
                  helpText("Auto-detect: lat/lon for small values, raster CRS for UTM/large values",
                           style = "font-size:11px; color:#95a5a6;")
                ),
                div(class = "setting-group",
                  tags$label(class = "setting-label", "Fuzziness exponent (m)"),
                  numericInput("fcm_m", NULL, value = 2, min = 1.1, max = 3, step = 0.1),
                  helpText("m = 2 is standard; higher m = fuzzier partitions", style = "font-size:11px; color:#95a5a6;")
                ),
                div(class = "setting-group",
                  tags$label(class = "setting-label", "PCA Variance Threshold (%)"),
                  numericInput("pca_thresh", NULL, value = 80, min = 50, max = 99, step = 1),
                  helpText("PCs retained to reach this % cumulative variance", style = "font-size:11px; color:#95a5a6;")
                )
              )
            ),

            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Data Summary")
              ),
              div(class = "card-body",
                uiOutput("data_summary")
              )
            )
          )
        ),

        fluidRow(
          column(12,
            div(class = "btn-primary-custom", id = "btn_load_data",
              icon("arrow-right"), " Load Data & Continue to Validation ",
              style = "max-width: 320px; display: inline-block; margin-top: 8px;")
          )
        )
      )
    ),

    # Tab 2: Validation ───────────────────────────────────────────────────────
    div(id = "tab-validate", style = "display:none;",
      div(class = "top-bar",
        div(h1("Validation — Optimal k Selection")),
        div(class = "top-bar-meta", "Step 2 of 6")
      ),
      div(class = "content-area",
        fluidRow(
          column(4,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Validation Settings")
              ),
              div(class = "card-body",
                div(class = "setting-group",
                  tags$label(class = "setting-label", "Min clusters (k)"),
                  numericInput("k_min", NULL, value = 2, min = 2, max = 10)
                ),
                div(class = "setting-group",
                  tags$label(class = "setting-label", "Max clusters (k)"),
                  numericInput("k_max", NULL, value = 6, min = 3, max = 15)
                ),
                hr(),
                div(class = "setting-group",
                  tags$label(class = "setting-label", "Optimal k Method"),
                  selectInput("opt_method", NULL,
                    choices = c(
                      "FPI — Fuzziness Performance Index" = "FPI",
                      "NCE — Normalized Classification Entropy" = "NCE",
                      "XB — Xie-Beni Index" = "XB",
                      "PE — Partition Entropy" = "PE",
                      "FS — Fukuyama-Sugeno Index" = "FS",
                      "Auto (FPI + NCE consensus)" = "auto"
                    ),
                    selected = "FPI"
                  )
                ),
                div(class = "info-msg",
                  HTML("<strong>FPI / NCE:</strong> argmax (peak)<br>
                        <strong>XB / PE / FS:</strong> argmin (valley)"))
              )
            ),

            div(class = "card",
              div(class = "card-body",
                div(class = "metric-row",
                  div(class = "metric-card",
                    div(class = "metric-value", textOutput("opt_k_display")),
                    div(class = "metric-label", "Optimal k")
                  ),
                  div(class = "metric-card green",
                    div(class = "metric-value", textOutput("opt_fpi_display")),
                    div(class = "metric-label", "FPI")
                  ),
                  div(class = "metric-card orange",
                    div(class = "metric-value", textOutput("opt_nce_display")),
                    div(class = "metric-label", "NCE")
                  )
                ),
                div(class = "metric-row",
                  div(class = "metric-card amber",
                    div(class = "metric-value", textOutput("opt_xb_display")),
                    div(class = "metric-label", "Xie-Beni")
                  )
                )
              )
            ),

            div(class = "btn-primary-custom", id = "btn_run_validation",
              icon("play"), " Run Validation ", style = "margin-top: 12px;")
          ),

          column(8,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Validation Index Plot"),
                div(style = "font-size:12px; color:#7f8c8d;",
                  textOutput("validation_method_label"))
              ),
              div(class = "card-body",
                withSpinner(plotlyOutput("plot_validation", height = "320px"), type = 4)
              )
            ),
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "All Indices Comparison")
              ),
              div(class = "card-body",
                withSpinner(plotlyOutput("plot_all_indices", height = "280px"), type = 4)
              )
            ),
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Validation Table"),
                div(style = "font-size:12px; color:#7f8c8d;", "★ = best per index")
              ),
              div(class = "card-body",
                DT::dataTableOutput("table_validation")
              )
            )
          )
        ),

        fluidRow(
          column(12,
            div(style = "margin-top: 16px;",
              actionLink("link_back_data", icon("arrow-left"), " Back to Data"),
              div(style = "float:right;",
                actionLink("link_next_cluster", "Continue to Clustering ", icon("arrow-right")))
            )
          )
        )
      )
    ),

    # Tab 3: Clustering ────────────────────────────────────────────────────────
    div(id = "tab-cluster", style = "display:none;",
      div(class = "top-bar",
        div(h1("Fuzzy C-Means Clustering")),
        div(class = "top-bar-meta", "Step 3 of 6")
      ),
      div(class = "content-area",
        fluidRow(
          column(3,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Cluster Settings")
              ),
              div(class = "card-body",
                div(class = "setting-group",
                  tags$label(class = "setting-label", "Number of Zones (k)"),
                  numericInput("manual_k", NULL, value = NA, min = 2, max = 15),
                  helpText("Leave blank to use optimal k from validation", style = "font-size:11px; color:#95a5a6;")
                ),
                div(class = "success-msg",
                  textOutput("cluster_k_info"))
              )
            ),
            div(class = "btn-primary-custom", id = "btn_run_clustering",
              icon("cogs"), " Run FCM Clustering ", style = "margin-top: 12px;")
          ),

          column(9,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Cluster Distribution")
              ),
              div(class = "card-body",
                withSpinner(plotlyOutput("plot_cluster_dist", height = "260px"), type = 4)
              )
            )
          ),

          column(12,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Zone Assignments")
              ),
              div(class = "card-body",
                DT::dataTableOutput("table_assignments")
              )
            )
          )
        ),

        fluidRow(
          column(12,
            div(style = "margin-top: 16px;",
              actionLink("link_back_validate", icon("arrow-left"), " Back to Validation"),
              div(style = "float:right;",
                actionLink("link_next_maps", "Continue to Zone Maps ", icon("arrow-right")))
            )
          )
        )
      )
    ),

    # Tab 4: Zone Maps ────────────────────────────────────────────────────────
    div(id = "tab-maps", style = "display:none;",
      div(class = "top-bar",
        div(h1("Zone Maps & Membership Surfaces")),
        div(class = "top-bar-meta", "Step 4 of 6")
      ),
      div(class = "content-area",
        fluidRow(
          column(8,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Management Zone Map"),
                actionButton("btn_regen_maps", "Regenerate",
                             icon = icon("rotate"),
                             class = "btn btn-sm",
                             style = "float:right; padding:2px 10px; font-size:12px;")
              ),
              div(class = "card-body",
                uiOutput("zone_status_ui"),
                plotOutput("plot_zone_map", height = "420px"),
                div(class = "info-msg", "Hard zone assignment: each pixel assigned to zone with highest membership value.")
              )
            )
          ),
          column(4,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Zone Summary")
              ),
              div(class = "card-body",
                uiOutput("zone_summary_ui")
              )
            )
          ),
          column(12,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Fuzzy Membership Surfaces"),
                div(style = "font-size:12px; color:#7f8c8d;",
                    "Bright = strong membership | Dark = weak membership")
              ),
              div(class = "card-body",
                plotOutput("plot_membership", height = "460px")
              )
            )
          )
        ),

        fluidRow(
          column(12,
            div(style = "margin-top: 16px;",
              actionLink("link_back_cluster", icon("arrow-left"), " Back to Clustering"),
              div(style = "float:right;",
                actionLink("link_next_stats", "Continue to Statistics ", icon("arrow-right")))
            )
          )
        )
      )
    ),

    # Tab 5: Statistics ────────────────────────────────────────────────────────
    div(id = "tab-stats", style = "display:none;",
      div(class = "top-bar",
        div(h1("Zone Statistics & ANOVA")),
        div(class = "top-bar-meta", "Step 5 of 6")
      ),
      div(class = "content-area",
        fluidRow(
          column(6,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Mean Soil Properties by Zone")
              ),
              div(class = "card-body",
                DT::dataTableOutput("table_zone_means")
              )
            )
          ),
          column(6,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "ANOVA — Zone Differences"),
                div(style = "font-size:12px; color:#7f8c8d;", "*** p<0.001 | ** p<0.01 | * p<0.05")
              ),
              div(class = "card-body",
                DT::dataTableOutput("table_anova")
              )
            )
          ),
          column(12,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Zone Means Comparison")
              ),
              div(class = "card-body",
                withSpinner(plotlyOutput("plot_zone_means", height = "320px"), type = 4)
              )
            )
          )
        ),

        fluidRow(
          column(12,
            div(style = "margin-top: 16px;",
              actionLink("link_back_maps", icon("arrow-left"), " Back to Zone Maps"),
              div(style = "float:right;",
                actionLink("link_next_export", "Continue to Export ", icon("arrow-right")))
            )
          )
        )
      )
    ),

    # Tab 6: Export ────────────────────────────────────────────────────────────
    div(id = "tab-export", style = "display:none;",
      div(class = "top-bar",
        div(h1("Export Results")),
        div(class = "top-bar-meta", "Step 6 of 6")
      ),
      div(class = "content-area",
        fluidRow(
          column(6,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Download Data")
              ),
              div(class = "card-body",
                div(style = "display: flex; flex-direction: column; gap: 12px;",
                  actionButton("btn_export_validation", "Download Validation CSV",
                               icon("download"), style = "background:#2c3e50; color:white; width:100%;"),
                  actionButton("btn_export_zone_stats", "Download Zone Statistics CSV",
                               icon("download"), style = "background:#2c3e50; color:white; width:100%;"),
                  actionButton("btn_export_assignments", "Download Point Assignments CSV",
                               icon("download"), style = "background:#2c3e50; color:white; width:100%;"),
                  actionButton("btn_export_anova", "Download ANOVA Results CSV",
                               icon("download"), style = "background:#2c3e50; color:white; width:100%;")
                )
              )
            ),
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Download Maps")
              ),
              div(class = "card-body",
                div(style = "display: flex; flex-direction: column; gap: 12px;",
                  actionButton("btn_export_zone_tif", "Download Zones GeoTIFF (.tif)",
                               icon("layer-group"), style = "background:#16a085; color:white; width:100%;"),
                  actionButton("btn_export_zone_png", "Download Zone Map PNG (300 DPI)",
                               icon("image"), style = "background:#16a085; color:white; width:100%;"),
                  div(style = "font-size:12px; color:#7f8c8d;",
                      "GeoTIFF keeps the CRS for GIS; PNG is a 300-DPI publication figure.")
                )
              )
            )
          ),
          column(6,
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Session Summary")
              ),
              div(class = "card-body",
                uiOutput("session_summary")
              )
            ),
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Generate HTML Report")
              ),
              div(class = "card-body",
                div(style = "font-size:12px; color:#7f8c8d; margin-bottom:10px;",
                    "Renders report/mz_report.qmd with the current session's results. ",
                    "Writes the missing outputs (including the transition mask) to ",
                    "outputs/ and runs the Quarto CLI; the resulting HTML is saved to ",
                    "report/mz_report.html."),
                div(class = "setting-group",
                    tags$label(class = "setting-label", "Author"),
                    textInput("report_author", NULL, value = "MZ Analysis",
                              placeholder = "e.g. C. Carbajal")
                ),
                div(class = "setting-group",
                    tags$label(class = "setting-label", "Study area name"),
                    textInput("report_study_area", NULL, value = "Study Area",
                              placeholder = "e.g. INIA Test Field")
                ),
                actionButton("btn_export_report", "Generate Report",
                             icon("file-alt"),
                             style = "background:#2980b9; color:white; width:100%; margin-top:6px;"),
                div(style = "font-size:11px; color:#95a5a6; margin-top:6px;",
                    "Renders report/mz_report.qmd with the current session's results. ",
                    "Writes the missing outputs (CSVs, zone map, validation index plot, ",
                    "transition mask) to outputs/ and runs the Quarto CLI; the resulting ",
                    "HTML lands at report/mz_report.html. Requires Quarto CLI on PATH ",
                    "(~3 s for the bundled data).")
              )
            ),
            div(class = "card",
              div(class = "card-header",
                span(class = "card-title", "Restart")
              ),
              div(class = "card-body",
                actionButton("btn_restart", "Start New Analysis",
                             icon("refresh"), style = "background:#e74c3c; color:white; width:100%;")
              )
            )
          )
        ),

        fluidRow(
          column(12,
            div(style = "margin-top: 16px;",
              actionLink("link_back_stats", icon("arrow-left"), " Back to Statistics")
            )
          )
        )
      )
    ),

    # Tab 0: User Guide ───────────────────────────────────────────────────────
    # Hand-built, in-app quick-start guide.  Uses the existing SVG diagrams
    # in www/assets/figures/ and a couple of sample-output PNGs — no
    # external markdown file, no runtime markdown package needed.  Keep it
    # basics-only: data contract, the 6 steps, a sample-output panel
    # (zone map + interpolated membership surfaces), and a few tips.
    # The full reference lives in the standalone
    # delimit_management_zones_USER_GUIDE.md for users who want depth.
    # This is the FIRST tab in the sidebar and the default-active one, so
    # the user lands on the guide the moment the app starts.
    div(id = "tab-help",
      div(class = "top-bar",
        div(h1("User Guide")),
        div(class = "top-bar-meta", "Quick start & reference")
      ),
      div(class = "content-area",

        # ── Hero ───────────────────────────────────────────────────────────
        div(class = "card help-hero",
          div(class = "card-body",
            tags$img(src = "mz_logo.svg",
                     style = "height:64px; width:64px; float:left;
                              margin-right:20px; margin-top:4px;"),
            h2(style = "margin:0 0 6px 0; color:#1a252f;",
               "MZ Delimitation"),
            div(style = "color:#3498db; font-weight:600; font-size:13px;
                        text-transform:uppercase; letter-spacing:1px;
                        margin-bottom:10px;",
               "Fuzzy C-Means management-zone delineation"),
            p(style = "color:#5d6d7e; font-size:14px; line-height:1.55;
                      margin-bottom:0;",
              "Upload a boundary, a soil-property raster stack, and a
              point CSV with measured values.  Pick three or more soil
              variables.  The app runs PCA → Fuzzy C-Means → kriged
              membership surfaces, then validates zones with FPI / NCE /
              XB and exports CSV / GeoTIFF / PNG.  Six steps, ~3 minutes
              end-to-end on the bundled SoilGrids demo data.")
          )
        ),

        # ── Data you need ──────────────────────────────────────────────────
        div(class = "card",
          div(class = "card-header",
            span(class = "card-title",
              icon("database"), " Data you need (3 files)"),
            div(style = "font-size:12px; color:#7f8c8d;",
                "All three are required before Step 1 can finish")
          ),
          div(class = "card-body",
            fluidRow(
              column(4, div(class = "help-file-box",
                div(class = "help-file-icon", icon("draw-polygon")),
                h4("Boundary", style = "margin:8px 0 4px 0;"),
                div(class = "help-file-ext", ".gpkg  ·  .shp  ·  .geojson"),
                p(style = "color:#5d6d7e; font-size:13px; margin:6px 0 0 0;",
                  "Polygon defining the field.  The raster is cropped and
                  masked to this extent.")
              )),
              column(4, div(class = "help-file-box",
                div(class = "help-file-icon", icon("image")),
                h4("Soil raster stack", style = "margin:8px 0 4px 0;"),
                div(class = "help-file-ext", ".tif  ·  .asc  ·  .bil"),
                p(style = "color:#5d6d7e; font-size:13px; margin:6px 0 0 0;",
                  "One band per soil property (BD, CEC, pH, SOC, …).
                  Layer names should match the CSV variable names when
                  possible so values can be cross-checked.")
              )),
              column(4, div(class = "help-file-box",
                div(class = "help-file-icon", icon("map-marker-alt")),
                h4("Observation points", style = "margin:8px 0 4px 0;"),
                div(class = "help-file-ext", ".csv"),
                p(style = "color:#5d6d7e; font-size:13px; margin:6px 0 0 0;",
                  "Longitude / latitude (or projected X / Y) plus the
                  numeric soil columns.  Include an ", code("fid"),
                  " column — it's used to join point assignments back to
                  your source data.")
              ))
            )
          )
        ),

        # ── The 6 steps ────────────────────────────────────────────────────
        div(class = "card",
          div(class = "card-header",
            span(class = "card-title",
              icon("layer-group"), " The 6-step workflow")),
          div(class = "card-body",
            tags$img(src = "assets/figures/fig1_workflow.svg",
                     class = "help-figure",
                     alt = "MZ Delimitation six-step workflow diagram"),
            div(class = "help-step-grid",
              div(class = "help-step",
                div(class = "help-step-num", "1"),
                div(strong("Data Input"), br(),
                    span(style = "color:#7f8c8d; font-size:12px;",
                         "Upload files, select variables, load & PCA"))),
              div(class = "help-step",
                div(class = "help-step-num", "2"),
                div(strong("Validation"), br(),
                    span(style = "color:#7f8c8d; font-size:12px;",
                         "FPI / NCE / XB / PE / FS over k = 2..K"))),
              div(class = "help-step",
                div(class = "help-step-num", "3"),
                div(strong("Clustering"), br(),
                    span(style = "color:#7f8c8d; font-size:12px;",
                         "Final FCM, hard zones, point assignments"))),
              div(class = "help-step",
                div(class = "help-step-num", "4"),
                div(strong("Zone Maps"), br(),
                    span(style = "color:#7f8c8d; font-size:12px;",
                         "Kriged membership surfaces + hard zone map"))),
              div(class = "help-step",
                div(class = "help-step-num", "5"),
                div(strong("Statistics"), br(),
                    span(style = "color:#7f8c8d; font-size:12px;",
                         "Per-zone means, ANOVA, comparison plot"))),
              div(class = "help-step",
                div(class = "help-step-num", "6"),
                div(strong("Export"), br(),
                    span(style = "color:#7f8c8d; font-size:12px;",
                         "CSV tables, GeoTIFF, 300-DPI PNG figure")))
            )
          )
        ),

        # ── Per-step procedure details ──────────────────────────────────────
        # One card per step that has a process diagram in
        # www/assets/figures/.  These complement the 6-step overview above
        # by showing the actual algorithm / data flow for each step.  The
        # final "Sample output" card below shows the end product (zone
        # map + membership surfaces) so the reader sees the full picture
        # from "how it works" to "what you get".

        # Step 1 — Data Input & Loading
        div(class = "card",
          div(class = "card-header",
            span(class = "card-title",
              icon("database"), " Step 1 — Data Input & Loading"),
            div(style = "font-size:12px; color:#7f8c8d;",
                "CRS alignment, raster extraction, median imputation,
                 min–max normalisation, PCA")
          ),
          div(class = "card-body",
            fluidRow(
              column(7,
                tags$img(src = "assets/figures/fig2_data_pipeline.svg",
                         class = "help-figure",
                         alt = "MZ data pipeline diagram"),
                div(class = "help-caption",
                  strong("Data pipeline"),
                  br(), span(style = "color:#7f8c8d;",
                    "Boundary + raster + CSV flow through CRS alignment,
                    point extraction, and feature matrix construction."))
              ),
              column(5,
                tags$img(src = "assets/figures/fig3_data_loading.svg",
                         class = "help-figure-tall",
                         alt = "Data loading detail"),
                div(class = "help-caption",
                  strong("Data loading detail"),
                  br(), span(style = "color:#7f8c8d;",
                    "Per-layer raster values overwrite matching CSV
                    columns; NAs are median-imputed, then min–max
                    normalised to [0, 1] before PCA."))
              )
            )
          )
        ),

        # Step 2 — Validation
        div(class = "card",
          div(class = "card-header",
            span(class = "card-title",
              icon("chart-line"), " Step 2 — Validation"),
            div(style = "font-size:12px; color:#7f8c8d;",
                "FPI / NCE / XB / PE / FS across k, consensus ranking")
          ),
          div(class = "card-body",
            tags$img(src = "assets/figures/fig4_validation.svg",
                     class = "help-figure-tall",
                     alt = "Validation indices diagram"),
            div(class = "help-caption",
              strong("Validation indices"),
              br(), span(style = "color:#7f8c8d;",
                "Five indices (FPI, NCE, XB, PE, FS) are computed for each
                k in 2..K.  Each picks an argmin (or argmax) and the
                consensus rank (XB + FPI + NCE) is the default tie-breaker."))
          )
        ),

        # Step 3 — Clustering
        div(class = "card",
          div(class = "card-header",
            span(class = "card-title",
              icon("layer-group"), " Step 3 — Clustering"),
            div(style = "font-size:12px; color:#7f8c8d;",
                "Fuzzy C-Means in PCA space, hard zone = argmax membership")
          ),
          div(class = "card-body",
            tags$img(src = "assets/figures/fig5_clustering.svg",
                     class = "help-figure-tall",
                     alt = "Fuzzy C-Means clustering diagram"),
            div(class = "help-caption",
              strong("Fuzzy C-Means"),
              br(), span(style = "color:#7f8c8d;",
                "Each point has a membership vector of length k (sums to 1).
                Hard cluster id is argmax of that vector.  Fuzziness
                exponent m = 2 is standard; higher m = softer boundaries."))
          )
        ),

        # Step 4 — Zone Maps & Interpolation
        div(class = "card",
          div(class = "card-header",
            span(class = "card-title",
              icon("map"), " Step 4 — Zone Maps & Interpolation"),
            div(style = "font-size:12px; color:#7f8c8d;",
                "Ordinary kriging of membership surfaces, hard zone map,
                 transition mask")
          ),
          div(class = "card-body",
            tags$img(src = "assets/figures/fig6_interpolation.svg",
                     class = "help-figure-tall",
                     alt = "Interpolation / kriging diagram"),
            div(class = "help-caption",
              strong("Kriging the membership surfaces"),
              br(), span(style = "color:#7f8c8d;",
                "For each zone, an ordinary-kriging spherical variogram
                is fit to the observation points and predicted onto the
                raster grid.  Hard zone = argmax across k surfaces.  IDW
                is used automatically if the variogram fails to converge."))
          )
        ),

        # ── Sample output ──────────────────────────────────────────────────
        # Only the zone map and the kriged membership surfaces — the two
        # outputs that come out of Step 4 (Zone Maps).  Bars/PCA/validation
        # plots are intentionally excluded; they live in Steps 2 and 5
        # of the running app and would just duplicate what the user sees
        # when they actually run the analysis.
        div(class = "card",
          div(class = "card-header",
            span(class = "card-title",
              icon("image"), " Sample output — Zone maps"),
            div(style = "font-size:12px; color:#7f8c8d;",
                "From the bundled SoilGrids demo (k = 3, m = 2)")
          ),
          div(class = "card-body",
            fluidRow(
              column(6,
                tags$img(src = "assets/figures/fig6b_example_zone_map.png",
                         class = "help-figure",
                         alt = "Example management zone map"),
                div(class = "help-caption",
                  strong("Management zone map"),
                  br(), span(style = "color:#7f8c8d;",
                    "Hard zone (k = 3).  Each pixel assigned to the
                    cluster with the highest fuzzy membership."))
              ),
              column(6,
                tags$img(src = "assets/figures/fig6c_example_membership_surfaces.png",
                         class = "help-figure",
                         alt = "Example membership surfaces"),
                div(class = "help-caption",
                  strong("Interpolated membership surfaces"),
                  br(), span(style = "color:#7f8c8d;",
                    "One panel per zone — the kriged membership
                    surfaces.  Bright = high membership, dark = low.
                    Overlap between panels = fuzziness."))
              )
            )
          )
        ),

        # ── Tips ───────────────────────────────────────────────────────────
        div(class = "card",
          div(class = "card-header",
            span(class = "card-title",
              icon("lightbulb"), " Quick tips"))
          ,
          div(class = "card-body",
            tags$ul(style = "color:#2c3e50; font-size:14px; line-height:1.7;",
              tags$li(strong("Variables: "),
                "Pick 3-8 soil variables that are conceptually distinct
                (e.g. BD, pH, SOC).  Highly correlated pairs (Sand +
                Silt) don't add information once PCA runs."),
              tags$li(strong("PCA threshold: "),
                "80% is a good default.  Lower it to 60-70% if you want
                a sharper separation, raise to 90% if interpretation
                matters more than speed."),
              tags$li(strong("Optimal k: "),
                "The default ", code("FPI"), " selection is fine for
                most cases.  If FPI and NCE disagree, switch the method
                in the validation tab and re-run."),
              tags$li(strong("CRS: "),
                "Auto-detect handles UTM and lat/lon, but if the raster
                is in a national projection and your points are
                lat/lon, set ", code("Point Coordinate System"),
                " to ", code("EPSG:4326"), " explicitly."),
              tags$li(strong("First run is slow: "),
                "Step 4 (kriging) is the slow one — expect 5-30 s on a
                typical field-sized raster.  Larger rasters are
                downsampled internally so the UI never freezes.")
            )
          )
        ),

        fluidRow(
          column(12,
            div(style = "margin-top: 16px;",
              actionLink("link_back_export", icon("arrow-left"), " Back to Export")
            )
          )
        )
      )
    )
  ),

  # ── Navigation JS ───────────────────────────────────────────────────────────
  tags$script(HTML("
    $(document).ready(function() {
      $('.nav-item').on('click', function() {
        var tabId = $(this).data('tab');
        $('.nav-item').removeClass('active');
        $(this).addClass('active');
        $('[id^=\"tab-\"]').hide();
        $('#' + tabId).show();
        window.scrollTo(0, 0);
        // CRITICAL: plots rendered while their tab was display:none come back
        // at size 0 and never repaint on their own.  Firing a resize forces
        // Shiny to re-measure the now-visible output and redraw it at full
        // size — this is what stops the Zone Maps spinner from spinning
        // forever even though the server finished rendering.
        //
        // For DataTables the resize is not enough on its own: the DT widget
        // also caches its column widths at init, so we explicitly call
        // columns.adjust() on every visible DataTable when the new tab
        // becomes visible.  This is what makes the Statistics tab tables
        // (table_zone_means / table_anova) actually appear after the tab
        // is shown — without it the eager-rendered DT widgets stay at
        // width 0 and look blank.
        setTimeout(function() {
          $(window).trigger('resize');
          if ($.fn.DataTable) {
            $('table.dataTable').each(function() {
              var dt = $(this).DataTable();
              if (dt && dt.columns && dt.columns.adjust) {
                dt.columns.adjust();
              }
            });
          }
        }, 250);
      });

      // Wire the custom <div> action buttons to Shiny.
      // Plain divs do not create Shiny input bindings, so without this the
      // observeEvent(input$btn_*) handlers on the server never fire.
      $(document).on('click', '#btn_load_data, #btn_run_validation, #btn_run_clustering', function() {
        Shiny.setInputValue($(this).attr('id'), Math.random(), {priority: 'event'});
      });
    });
  "))
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

  mzlog("===== NEW SESSION ===== app version: ", MZ_VERSION,
        " | terra ", as.character(packageVersion("terra")),
        " | log: ", MZ_LOGFILE)

  rv <- reactiveValues(
    boundary       = NULL,
    soil_masked    = NULL,
    obs_sf         = NULL,
    soil_vars      = NULL,
    pca_scores     = NULL,
    pca_raw        = NULL,
    prop_var       = NULL,
    n_pc           = NULL,
    df_vars        = NULL,
    df_norm        = NULL,
    validation_df  = NULL,
    selected_k     = NULL,
    fcm_result     = NULL,
    zone_stack     = NULL,
    zone_hard      = NULL,
    k_final        = NULL,
    # ── CSV column detection ──
    csv_columns    = NULL,   # all column names from uploaded CSV
    numeric_cols   = NULL,   # detected numeric columns from CSV
    coord_cols     = NULL    # detected coordinate columns
  )

  # ── Navigation helpers ─────────────────────────────────────────────────────
  observeEvent(TRUE, {
    observe({shinyjs::onclick("link_back_data",      function(e) { shinyjs::runjs("$('.nav-item').removeClass('active'); $('#nav-data').addClass('active'); $('[id^=\"tab-\"]').hide(); $('#tab-data').show(); window.scrollTo(0,0);")})})
    observe({shinyjs::onclick("link_next_cluster",   function(e) { shinyjs::runjs("$('.nav-item').removeClass('active'); $('#nav-cluster').addClass('active'); $('[id^=\"tab-\"]').hide(); $('#tab-cluster').show(); window.scrollTo(0,0);")})})
    observe({shinyjs::onclick("link_back_validate",  function(e) { shinyjs::runjs("$('.nav-item').removeClass('active'); $('#nav-validate').addClass('active'); $('[id^=\"tab-\"]').hide(); $('#tab-validate').show(); window.scrollTo(0,0);")})})
    observe({shinyjs::onclick("link_next_maps",      function(e) { shinyjs::runjs("$('.nav-item').removeClass('active'); $('#nav-maps').addClass('active'); $('[id^=\"tab-\"]').hide(); $('#tab-maps').show(); window.scrollTo(0,0); setTimeout(function(){$(window).trigger('resize');},250);")})})
    observe({shinyjs::onclick("link_back_cluster",   function(e) { shinyjs::runjs("$('.nav-item').removeClass('active'); $('#nav-cluster').addClass('active'); $('[id^=\"tab-\"]').hide(); $('#tab-cluster').show(); window.scrollTo(0,0);")})})
    observe({shinyjs::onclick("link_next_stats",     function(e) { shinyjs::runjs("$('.nav-item').removeClass('active'); $('#nav-stats').addClass('active'); $('[id^=\"tab-\"]').hide(); $('#tab-stats').show(); window.scrollTo(0,0); setTimeout(function(){$(window).trigger('resize'); if(window.jQuery&&jQuery.fn.DataTable){jQuery('table.dataTable').each(function(){var dt=jQuery(this).DataTable(); if(dt&&dt.columns&&dt.columns.adjust) dt.columns.adjust();});} },250);")})})
    observe({shinyjs::onclick("link_back_maps",      function(e) { shinyjs::runjs("$('.nav-item').removeClass('active'); $('#nav-maps').addClass('active'); $('[id^=\"tab-\"]').hide(); $('#tab-maps').show(); window.scrollTo(0,0); setTimeout(function(){$(window).trigger('resize');},250);")})})
    observe({shinyjs::onclick("link_next_export",    function(e) { shinyjs::runjs("$('.nav-item').removeClass('active'); $('#nav-export').addClass('active'); $('[id^=\"tab-\"]').hide(); $('#tab-export').show(); window.scrollTo(0,0);")})})
    observe({shinyjs::onclick("link_back_stats",     function(e) { shinyjs::runjs("$('.nav-item').removeClass('active'); $('#nav-stats').addClass('active'); $('[id^=\"tab-\"]').hide(); $('#tab-stats').show(); window.scrollTo(0,0); setTimeout(function(){$(window).trigger('resize'); if(window.jQuery&&jQuery.fn.DataTable){jQuery('table.dataTable').each(function(){var dt=jQuery(this).DataTable(); if(dt&&dt.columns&&dt.columns.adjust) dt.columns.adjust();});} },250);")})})
    observe({shinyjs::onclick("link_back_export",    function(e) { shinyjs::runjs("$('.nav-item').removeClass('active'); $('#nav-export').addClass('active'); $('[id^=\"tab-\"]').hide(); $('#tab-export').show(); window.scrollTo(0,0);")})})
  }, once = TRUE)

  # Default state before CSV upload
  output$soil_vars_ui <- renderUI({
    div(style = "color: #bdc3c7; font-style: italic; font-size: 13px;",
      "Upload the observation CSV first to see available variables.")
  })
  output$var_col_preview <- renderUI(NULL)

  # ── CSV column detection on upload ─────────────────────────────────────────
  observeEvent(input$input_points, {
    req(input$input_points$datapath)

    tryCatch({
      df_raw <- read.csv(input$input_points$datapath, nrows = 10, header = TRUE)

      # Detect numeric columns (exclude coordinates)
      all_cols <- names(df_raw)

      # Find potential coordinate columns
      coord_patterns <- c("lon", "lat", "x$", "y$", "easting", "northing",
                         "long", "coord", "geometry", "fid", "id")
      coord_mask <- sapply(all_cols, function(col) {
        any(sapply(coord_patterns, function(p) grepl(p, col, ignore.case = TRUE)))
      })

      # Numeric detection from sample
      num_mask <- sapply(df_raw, function(col) is.numeric(col) || is.integer(col))

      # Candidate numeric vars: numeric columns that are NOT coordinate columns
      coord_cols_detected <- all_cols[coord_mask]
      numeric_candidates   <- all_cols[num_mask & !coord_mask]

      rv$csv_columns  <- all_cols
      rv$coord_cols  <- coord_cols_detected
      rv$numeric_cols <- numeric_candidates

      # Auto-suggest coordinate columns
      lon_candidates <- all_cols[grepl("lon|long|east|x$", all_cols, ignore.case = TRUE)]
      lat_candidates <- all_cols[grepl("lat|y$|north", all_cols, ignore.case = TRUE)]

      if (length(lon_candidates) > 0) updateTextInput(session, "col_lon", value = lon_candidates[1])
      if (length(lat_candidates) > 0) updateTextInput(session, "col_lat", value = lat_candidates[1])

      # Build checkbox UI for numeric variables
      output$soil_vars_ui <- renderUI({
        req(length(numeric_candidates) > 0)

        div(
          checkboxGroupInput("soil_vars", NULL,
            choices = setNames(numeric_candidates, numeric_candidates),
            selected = numeric_candidates[1:min(8, length(numeric_candidates))],
            inline = FALSE
          ),
          div(style = "margin-top: 8px;",
            actionLink("select_all_vars", "[Select All]", style = "font-size:12px;"),
            actionLink("deselect_all_vars", "[Deselect All]", style = "font-size:12px; margin-left: 10px;")
          )
        )
      })

      # Preview of numeric column data
      output$var_col_preview <- renderUI({
        req(numeric_candidates)
        preview_cols <- numeric_candidates[1:min(6, length(numeric_candidates))]
        preview_df   <- df_raw[, preview_cols, drop = FALSE]
        col_classes  <- sapply(preview_df, function(x) class(x)[1])

        tagList(
          div(style = "margin-top: 12px; font-size: 12px; color: #7f8c8d;",
            strong("Detected numeric columns (first 6 shown):"),
            HTML(paste0("<br>",
              paste(sprintf("<span style='font-family:monospace;'>%s</span> <span style='color:#bdc3c7;'>(%s)</span>",
                preview_cols, col_classes), collapse = " &nbsp; ")))
          )
        )
      })

      output$points_status <- renderUI({
        n_num   <- length(numeric_candidates)
        n_coord <- length(coord_cols_detected)
        tagList(
          div(class = "success-msg",
            icon("check-circle"), strong(basename(input$input_points$name)),
            br(), sprintf("%d rows, %d columns detected", nrow(df_raw), length(all_cols)),
            br(), sprintf("%d numeric variables, %d coordinate columns", n_num, n_coord)
          )
        )
      })

    }, error = function(e) {
      output$points_status <- renderUI(
        div(class = "warning-msg", icon("exclamation-triangle"),
          "Could not read CSV columns: ", e$message)
      )
    })
  })

  # ── Select / Deselect All ────────────────────────────────────────────────────
  observeEvent(input$select_all_vars, {
    updateCheckboxGroupInput(session, "soil_vars", selected = rv$numeric_cols)
  })
  observeEvent(input$deselect_all_vars, {
    updateCheckboxGroupInput(session, "soil_vars", selected = character(0))
  })

  # ── Data loading ──────────────────────────────────────────────────────────

  # Top-level Data Summary card.  Reads rv$obs_sf / rv$soil_vars / rv$n_pc /
  # rv$point_crs_label so it auto-refreshes as soon as the user clicks
  # "Load Data & Continue".  Before data is loaded, it shows a friendly
  # pointer instead of an empty card body (the previous version was
  # defined as a side-effect inside the btn_load_data observe, so the card
  # was blank for any user who uploaded files but hadn't clicked Load yet).
  output$data_summary <- renderUI({
    if (is.null(rv$obs_sf) || is.null(rv$soil_vars)) {
      return(div(class = "info-msg",
        icon("info-circle"),
        "Data not loaded yet. Upload the three files and click",
        strong(" Load Data & Continue "),
        "to populate this card."))
    }
    tagList(
      div(icon("check"), strong(nrow(rv$obs_sf)), "observation points"),
      div(icon("check"), strong(length(rv$soil_vars)), "soil variables selected"),
      div(icon("check"), strong(rv$n_pc), "PC(s) retained (≥",
          input$pca_thresh, "% var)"),
      div(icon("check"), "Point CRS: ",
          code(rv$point_crs_label %||% "—")),
      div(style = "font-size:11px; color:#7f8c8d; margin-top:4px;",
        "Variables: ", paste(rv$soil_vars, collapse = ", "))
    )
  })

  observeEvent(input$btn_load_data, {
    req(input$input_boundary$datapath, input$input_raster$datapath, input$input_points$datapath)

    # Validate: at least 3 variables must be selected
    selected_vars <- input$soil_vars
    if (length(selected_vars) < 3 || is.null(selected_vars)) {
      showNotification("Please select at least 3 soil variables from your CSV before loading data.",
                       type = "error", duration = 6)
      return()
    }

    shinyjs::disable("btn_load_data")
    shinyjs::html("btn_load_data", '<span style="color:white;"><i class="fa fa-spinner fa-spin"></i> Loading boundary...</span>', add = FALSE)

    tryCatch({

      # ── Stage 1: Boundary ────────────────────────────────────────────────
      output$boundary_status <- renderUI(
        div(class = "info-msg", HTML('<i class="fa fa-spinner fa-spin"></i>'), " Loading boundary...")
      )
      rv$boundary <- loadBoundary(input$input_boundary$datapath)
      output$boundary_status <- renderUI({
        div(class = "success-msg", icon("check"), "Boundary: ",
            basename(input$input_boundary$name), " — ",
            nrow(rv$boundary), "polygon(s)")
      })

      # ── Stage 2: Raster ─────────────────────────────────────────────────
      shinyjs::html("btn_load_data", '<span style="color:white;"><i class="fa fa-spinner fa-spin"></i> Loading raster...</span>', add = FALSE)
      output$raster_status <- renderUI(
        div(class = "info-msg", HTML('<i class="fa fa-spinner fa-spin"></i>'), " Loading raster (large file — please wait)...")
      )
      soil_stack   <- loadRaster(input$input_raster$datapath)

      # CRS alignment
      raster_crs   <- crs(soil_stack)
      boundary_crs <- st_crs(rv$boundary)
      if (!identical(as.character(raster_crs), as.character(boundary_crs))) {
        rv$boundary <- st_transform(rv$boundary, raster_crs)
      }

      boundary_terra <- vect(rv$boundary)
      rv$soil_masked <- crop(soil_stack, boundary_terra)
      rv$soil_masked <- mask(rv$soil_masked, boundary_terra)

      output$raster_status <- renderUI({
        div(class = "success-msg", icon("check"), "Raster: ",
            basename(input$input_raster$name), " — ",
            paste(names(rv$soil_masked), collapse = ", "))
      })

      # ── Stage 3: Points ──────────────────────────────────────────────────
      shinyjs::html("btn_load_data", '<span style="color:white;"><i class="fa fa-spinner fa-spin"></i> Loading points...</span>', add = FALSE)
      output$points_status <- renderUI(
        div(class = "info-msg", HTML('<i class="fa fa-spinner fa-spin"></i>'), " Loading observation points...")
      )
      obs_raw <- loadPointsCSV(input$input_points$datapath)
      names(obs_raw)[names(obs_raw) == input$col_lon] <- "X"
      names(obs_raw)[names(obs_raw) == input$col_lat] <- "Y"

      # ── Resolve point CRS ────────────────────────────────────────────────
      crs_mode <- input$point_crs_mode %||% "auto"
      if (crs_mode == "auto") {
        # Heuristic: if coordinates exceed geographic bounds (±180 / ±90),
        # they are projected (e.g., UTM). Otherwise assume WGS84 lat/lon.
        x_vals <- obs_raw[["X"]]
        y_vals <- obs_raw[["Y"]]
        is_proj <- max(abs(x_vals), na.rm = TRUE) > 180 ||
                   max(abs(y_vals), na.rm = TRUE) > 90
        if (is_proj) {
          point_crs <- crs(rv$soil_masked)       # use raster CRS (e.g., UTM zone)
          crs_label <- paste0("Projected (", as.character(crs(rv$soil_masked)), ")")
        } else {
          point_crs <- 4326
          crs_label <- "EPSG:4326 (WGS84)"
        }
      } else if (crs_mode == "raster") {
        point_crs <- crs(rv$soil_masked)
        crs_label <- paste0("Raster CRS (", as.character(crs(rv$soil_masked)), ")")
      } else {
        point_crs <- as.numeric(crs_mode)
        crs_label <- paste0("EPSG:", crs_mode)
      }
      rv$point_crs_label <- crs_label    # for display in data summary

      output$points_status <- renderUI({
        div(class = "success-msg", icon("check"), "Points: ",
            basename(input$input_points$name), " — ",
            nrow(obs_raw), "observations loaded")
      })

      # ── Stage 4: Feature matrix + PCA ────────────────────────────────────
      shinyjs::html("btn_load_data", '<span style="color:white;"><i class="fa fa-spinner fa-spin"></i> Computing PCA...</span>', add = FALSE)
      rv$soil_vars <- input$soil_vars

      fm <- buildFeatureMatrix(obs_raw, rv$soil_masked, rv$soil_vars,
                             pca_thresh = input$pca_thresh / 100,
                             point_crs  = point_crs)
      rv$pca_scores <- fm$pca_scores
      rv$pca_raw    <- fm$pca_raw
      rv$prop_var   <- fm$prop_var
      rv$n_pc       <- fm$n_pc
      rv$df_vars    <- fm$df_vars
      rv$df_norm    <- fm$df_norm
      rv$obs_sf     <- fm$obs_sf

      # ── Done ─────────────────────────────────────────────────────────────
      shinyjs::html("btn_load_data", '<i class="fa fa-arrow-right"></i> Load Data &amp; Continue to Validation ', add = FALSE)
      shinyjs::enable("btn_load_data")

      # NOTE: output$data_summary is defined at the TOP LEVEL of the server
      # (see below).  It is a regular reactive that reads rv$obs_sf,
      # rv$soil_vars, etc., so it auto-refreshes on its own — no need to
      # reassign it from inside this observe.  The previous side-effect
      # pattern (`output$data_summary <- renderUI({...})` here) is exactly
      # the anti-pattern the rest of the app was refactored to avoid, and
      # it left the Data Summary card empty before the first click of
      # "Load Data & Continue".

      showNotification(paste0("Data loaded! ", nrow(rv$obs_sf), " points, ",
                               length(rv$soil_vars), " variables, ",
                               rv$n_pc, " PCs. Proceed to Step 2."),
                       type = "message", duration = 5)

    }, error = function(e) {
      msg <- paste0("Load failed: ", e$message)
      output$boundary_status <- renderUI(div(class = "warning-msg", icon("exclamation-triangle"), msg))
      output$raster_status   <- renderUI(div(class = "warning-msg", icon("exclamation-triangle"), msg))
      output$points_status   <- renderUI(div(class = "warning-msg", icon("exclamation-triangle"), msg))
      shinyjs::html("btn_load_data", '<i class="fa fa-arrow-right"></i> Load Data &amp; Continue to Validation ', add = FALSE)
      shinyjs::enable("btn_load_data")
      showNotification(msg, type = "error", duration = 8)
    })
  })

  # ── Validation ────────────────────────────────────────────────────────────
  observeEvent(input$btn_run_validation, {
    req(rv$pca_scores)

    k_range <- seq(input$k_min, input$k_max)
    rv$validation_df <- runValidation(rv$pca_scores, k_range = k_range, m = input$fcm_m)
    rv$selected_k    <- selectOptimalK(rv$validation_df, input$opt_method)
    updateNumericInput(session, "manual_k", value = rv$selected_k)

    method <- input$opt_method

    # Metrics
    sel <- rv$validation_df[rv$validation_df$k == rv$selected_k, ]
    output$opt_k_display  <- renderText(paste0("k = ", rv$selected_k))
    output$opt_fpi_display <- renderText(sprintf("%.4f", sel$FPI))
    output$opt_nce_display <- renderText(sprintf("%.4f", sel$NCE))
    output$opt_xb_display  <- renderText(sprintf("%.4f", sel$XB))
    output$validation_method_label <- renderText(paste0("Method: ", method,
      if (method %in% c("FPI", "NCE", "auto")) " (argmax)" else " (argmin)"))

    # Main validation plot
    output$plot_validation <- renderPlotly({
      df <- rv$validation_df
      p <- ggplot(df, aes(x = .data[["k"]], y = .data[[method]])) +
        geom_line(color = "#3498db", linewidth = 1.2) +
        geom_point(color = "#2c3e50", size = 4) +
        geom_point(data = df[df$k == rv$selected_k, , drop = FALSE],
                   aes(x = .data[["k"]], y = .data[[method]]),
                   color = "#e74c3c", size = 7, shape = 21, stroke = 2) +
        labs(title = paste0(method, " Index by Number of Zones"),
             x = "Number of Zones (k)", y = method) +
        theme_minimal(base_size = 13) +
        theme(plot.title = element_text(face = "bold", color = "#1a252f"),
              plot.margin = margin(10,10,10,10))
      # ggplotly drops ggplot subtitles, so carry "Optimal k" in the plotly title
      ggplotly(p) %>% layout(
        hovermode = "x unified",
        title = list(text = paste0(method,
          " Index by Number of Zones<br><sup>Optimal k = ", rv$selected_k, "</sup>")))
    })

    # All indices comparison
    output$plot_all_indices <- renderPlotly({
      req(rv$validation_df)
      df_long <- tidyr::gather(rv$validation_df, key = "Index", value = "Value", -k)
      p <- ggplot(df_long, aes(x = k, y = Value, color = Index)) +
        geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
        facet_wrap(~Index, ncol = 3, scales = "free_y") +
        labs(title = "All Validation Indices", x = "k", y = "Value") +
        theme_minimal(base_size = 11) +
        theme(legend.position = "none",
              strip.background = element_rect(fill = "#ecf0f1"),
              plot.title = element_text(face = "bold", color = "#1a252f"))
      ggplotly(p) %>% layout(hovermode = "x unified")
    })

    # Validation table
    output$table_validation <- DT::renderDataTable({
      df <- rv$validation_df
      df$XB_rank  <- rank(df$XB,  ties.method = "min")
      df$FPI_rank <- rank(df$FPI, ties.method = "min")
      df$NCE_rank <- rank(df$NCE, ties.method = "min")
      df$consensus_rank <- df$XB_rank + df$FPI_rank + df$NCE_rank

      df$Best_XB  <- ifelse(df$k == df$k[which.min(df$XB)],  "★", "")
      df$Best_PE  <- ifelse(df$k == df$k[which.min(df$PE)],  "★", "")
      df$Best_FS  <- ifelse(df$k == df$k[which.min(df$FS)],  "★", "")
      df$Best_FPI <- ifelse(df$k == df$k[which.min(df$FPI)], "★", "")
      df$Best_NCE <- ifelse(df$k == df$k[which.min(df$NCE)], "★", "")
      df$Best_Consensus <- ifelse(df$k == df$k[which.min(df$consensus_rank)], "★", "")
      df$Sel      <- ifelse(df$k == rv$selected_k, "✓", "")

      DT::datatable(df, options = list(pageLength = 8, dom = "t", ordering = FALSE),
                    rownames = FALSE, caption = paste0("Selected k = ", rv$selected_k)) %>%
        DT::formatSignif(c("XB", "PE", "FS", "FPI", "NCE"), 4)
    })

    showNotification(paste0("Validation done. Optimal k = ", rv$selected_k, " (", input$opt_method, ")"), type = "message")
  })

  # ── Clustering ─────────────────────────────────────────────────────────────
  observeEvent(input$btn_run_clustering, {
    req(rv$pca_scores, rv$selected_k)

    # manual_k is a free-text numericInput: coerce to integer and bound to the data
    k <- if (!is.na(input$manual_k) && input$manual_k >= 2) {
      as.integer(round(input$manual_k))
    } else {
      rv$selected_k
    }
    k_max <- nrow(rv$pca_scores) - 1
    validate(
      need(k >= 2, "Number of zones (k) must be at least 2."),
      need(k <= k_max, sprintf(
        "Number of zones (k = %d) cannot exceed the number of observations - 1 (max %d).",
        k, k_max))
    )
    rv$k_final <- k

    # Guard the clustering: if FCM fails (e.g. fewer distinct data points than
    # k, usually because median-imputation collapsed points that fell outside
    # the raster), surface a clear error and STOP here.  Leaving rv$fcm_result
    # NULL is exactly what made Part 4 wait forever with a spinning loader.
    mzlog("CLUSTER: start runFCM k=", k, " m=", input$fcm_m,
          " pca_rows=", nrow(rv$pca_scores))
    fcm <- tryCatch(
      runFCM(rv$pca_scores, k = k, m = input$fcm_m),
      error = function(e) e)
    if (inherits(fcm, "error")) {
      mzlog("CLUSTER: FAILED — ", conditionMessage(fcm))
      showNotification(
        paste0("Clustering failed: ", conditionMessage(fcm),
               " — try a smaller k, or revisit your variable selection / boundary."),
        type = "error", duration = 12)
      return()
    }
    mzlog("CLUSTER: ok, fcm_result set")
    rv$fcm_result <- fcm
    cluster_id <- rv$fcm_result$cluster_id
    membership <- rv$fcm_result$membership

    output$cluster_k_info <- renderText(paste0("Running with k = ", k))

    # Cluster distribution — tabulate over 1..k so an empty hard cluster shows as 0
    dist_df <- data.frame(Zone = seq_len(k),
                          Count = tabulate(cluster_id, nbins = k))

    output$plot_cluster_dist <- renderPlotly({
      p <- ggplot(dist_df, aes(x = as.factor(Zone), y = Count, fill = as.factor(Zone))) +
        geom_col(width = 0.6) +
        geom_text(aes(label = Count), vjust = -0.3, size = 4) +
        scale_fill_brewer(palette = "Set2") +
        labs(title = paste0("Points per Zone (k = ", k, ")"),
             x = "Zone", y = "Number of Points") +
        theme_minimal(base_size = 13) +
        theme(legend.position = "none", plot.title = element_text(face = "bold"))
      ggplotly(p)
    })

    # Assignments table
    rv$obs_sf$Zone <- cluster_id
    if (!"fid" %in% names(rv$obs_sf)) rv$obs_sf$fid <- seq_len(nrow(rv$obs_sf))
    assign_df <- as.data.frame(st_drop_geometry(rv$obs_sf[, c("fid", "Zone")]))
    output$table_assignments <- DT::renderDataTable({
      DT::datatable(assign_df, options = list(pageLength = 10), rownames = FALSE)
    })

    showNotification(paste0("FCM clustering done — ", k, " zones"), type = "message")
  })

  # ── Zone Maps ──────────────────────────────────────────────────────────────
  # ARCHITECTURE NOTE — why the Zone Maps spinner used to spin forever:
  # The previous version created `output$plot_zone_map` (and the membership
  # outputs) as a SIDE EFFECT *inside* a plain observe().  An output created
  # that way only exists once the observe runs all the way to that line, so
  # if the heavy buildZoneMap()/kriging step was slow, hit an error, or a
  # req() upstream halted, `output$plot_zone_map` was never defined and the
  # withSpinner() over it kept spinning with no plot and no visible error —
  # exactly the "bucle"/loop the user reported.
  #
  # The fix: do the heavy work ONCE in a cached reactive (zoneMap), and
  # define every Zone Maps output at the TOP LEVEL of the server (below),
  # reading that reactive.  Now the outputs always exist, the spinner always
  # resolves, and any failure surfaces as a readable on-screen message
  # instead of an endless spinner.
  zone_palette <- c("#d7191c", "#fdae61", "#abdda4", "#2c7bb6", "#5e4fa2",
                    "#1a9850", "#fc8d59", "#91bfdb", "#7fc97f", "#beaed4")
  # Color ramp that gracefully extends the base palette if k > its length.
  zoneColors <- function(k) {
    if (k <= length(zone_palette)) zone_palette[seq_len(k)]
    else grDevices::colorRampPalette(zone_palette)(k)
  }

  # Tag a hard-zone raster with the FULL 1..k category table AND an explicit
  # value->colour table, then plot it WITHOUT a `col=` argument so terra uses
  # that table.  This is what makes each zone draw in its own colour: a plain
  # integer raster is rendered as CONTINUOUS, so terra stretches a `col` vector
  # across the data's actual [min,max].  When a zone id is absent — e.g. dropped
  # by modal aggregation for the on-screen map, or simply never the argmax — the
  # colours shift and stop matching the zone-% legend and the exported PNG (the
  # on-screen and full-res renders could even disagree).  Empirically (terra
  # 1.9.27) only a coltab pins value z -> zoneColors(k)[z] regardless of which
  # ids occur; setting factor levels alone does NOT (a col= vector is still
  # stretched).  The levels are kept purely for the GeoTIFF's "Zone N" labels.
  asCategoricalZones <- function(r, k) {
    levels(r)        <- data.frame(value = seq_len(k),
                                   Zone  = paste0("Zone ", seq_len(k)))
    terra::coltab(r) <- data.frame(value = seq_len(k), col = zoneColors(k))
    r
  }

  # Manual "Regenerate" trigger — lets the user force a fresh recompute from
  # the Zone Maps tab if anything ever looks stuck.
  maps_trigger <- reactiveVal(0L)
  observeEvent(input$btn_regen_maps, maps_trigger(isolate(maps_trigger()) + 1L))

  # Heavy zone-map computation, cached so it runs ONCE per clustering result
  # and is shared by every output below.  req() (before the tryCatch) keeps
  # the outputs in their "waiting" state until clustering is done; a genuine
  # failure is returned as ok = FALSE so the renderers can show a message
  # rather than dying and leaving the spinner stuck.
  zoneMap <- reactive({
    maps_trigger()  # depend on the manual trigger so the button forces a redo

    # Self-diagnosing prerequisites.  We deliberately DO NOT use req() here:
    # req() halts silently, which left the output value-less and the spinner
    # spinning forever with no explanation.  Instead we report exactly what is
    # missing so the renderers can show it on screen.
    miss <- c()
    if (is.null(rv$boundary))    miss <- c(miss, "boundary (Step 1)")
    if (is.null(rv$soil_masked)) miss <- c(miss, "raster (Step 1)")
    if (is.null(rv$obs_sf))      miss <- c(miss, "observation points (Step 1)")
    if (is.null(rv$fcm_result))  miss <- c(miss, "clustering result (Step 3)")
    if (length(miss) > 0) {
      mzlog("ZONEMAP: waiting — missing: ", paste(miss, collapse = ", "))
      return(list(ok = FALSE, waiting = TRUE,
                  message = paste0("Waiting for: ", paste(miss, collapse = ", "),
                                   ". Complete that step, then return here.")))
    }

    k          <- rv$k_final
    membership <- rv$fcm_result$membership
    cluster_id <- rv$fcm_result$cluster_id
    mzlog("ZONEMAP: prereqs ok | k=", k,
          " | soil_masked cells=", terra::ncell(rv$soil_masked),
          " | obs=", nrow(rv$obs_sf),
          " | raster_crs=", terra::crs(rv$soil_masked, describe = TRUE)$code)
    boundary_terra <- vect(rv$boundary)

    tryCatch({
      mzlog("ZONEMAP: buildZoneMap start")
      zone_stack <- buildZoneMap(rv$obs_sf, membership, boundary_terra,
                                 rv$soil_masked, k,
                                 method = "kriging", max_grid_cells = 40000)
      mzlog("ZONEMAP: buildZoneMap done, nlyr=", terra::nlyr(zone_stack))
      rv$zone_stack <- zone_stack  # full-res, kept for downstream export

      # Hard zone = layer index of the max membership per cell.
      # terra::which.max() is a vectorised C++ op (~0.1s on the bundled
      # stack); the old app()-with-R-closure took 25s+ and was a second
      # reason Part 4 appeared to hang.
      zone_hard <- mask(which.max(zone_stack), boundary_terra)
      names(zone_hard) <- "Zone"
      mzlog("ZONEMAP: which.max+mask done")

      # Downsample for display so terra::plot / ggplotly stay responsive.
      # Modal aggregation preserves the dominant zone label per block.
      DISPLAY_MAX_CELLS <- 500000
      zone_hard_disp <- if (ncell(zone_hard) > DISPLAY_MAX_CELLS) {
        aggregate(zone_hard,
                  fact = ceiling(sqrt(ncell(zone_hard) / DISPLAY_MAX_CELLS)),
                  fun = "modal", na.rm = TRUE)
      } else zone_hard

      # Pin both the full-res and the display rasters to the FULL 1..k category
      # table AFTER aggregation, so each zone id keeps its own colour even when
      # modal aggregation drops a minority zone — keeping the on-screen map, the
      # exported PNG and the zone-% legend in agreement (see asCategoricalZones).
      zone_hard      <- asCategoricalZones(zone_hard, k)
      zone_hard_disp <- asCategoricalZones(zone_hard_disp, k)
      rv$zone_hard   <- zone_hard  # full-res categorical zones, for GeoTIFF/PNG

      # Downsample the WHOLE membership stack for a single static panel plot.
      # Rendering all k surfaces as one server-side image (terra::plot) avoids
      # shipping k interactive plotly widgets to the browser — those k
      # client-side geom_raster widgets were the last element on this tab heavy
      # enough to make the page feel frozen ("loop") on some machines.
      MEMB_MAX_CELLS <- 120000
      memb_stack_disp <- if (ncell(zone_stack) > MEMB_MAX_CELLS) {
        aggregate(zone_stack,
                  fact = ceiling(sqrt(ncell(zone_stack) / MEMB_MAX_CELLS)),
                  fun = mean, na.rm = TRUE)
      } else zone_stack
      names(memb_stack_disp) <- paste0("Zone ", seq_len(k))

      mzlog("ZONEMAP: SUCCESS (ok=TRUE), display cells=",
            terra::ncell(zone_hard_disp))
      list(ok = TRUE, k = k, cluster_id = cluster_id,
           zone_hard_disp = zone_hard_disp, memb_stack_disp = memb_stack_disp)
    }, error = function(e) {
      mzlog("ZONEMAP: ERROR — ", conditionMessage(e))
      list(ok = FALSE, waiting = FALSE, message = conditionMessage(e))
    })
  })

  # Shared placeholder painter so the map / membership cards ALWAYS show a
  # readable message (never an endless spinner) when the result isn't ready.
  drawPlaceholder <- function(zm) {
    op <- graphics::par(mar = c(0, 0, 0, 0)); on.exit(graphics::par(op))
    plot.new()
    waiting <- isTRUE(zm$waiting)
    text(0.5, 0.58,
         if (waiting) "Zone maps not ready yet" else "Zone map could not be generated",
         cex = 1.25, font = 2, col = if (waiting) "#2c3e50" else "#c0392b")
    text(0.5, 0.44, zm$message %||% "", cex = 0.95, col = "#7f8c8d")
    if (!waiting)
      text(0.5, 0.34, "Click “Regenerate” to try again.",
           cex = 0.9, col = "#7f8c8d")
    invisible(NULL)
  }

  # Kick the computation off (and surface progress) as soon as clustering is
  # ready, instead of lazily on tab-open.  zoneMap() is cached, so the
  # outputs below reuse this result rather than recomputing.
  observeEvent(
    list(rv$fcm_result, rv$obs_sf, rv$soil_masked, rv$boundary, rv$k_final),
    {
      req(rv$fcm_result, rv$obs_sf, rv$soil_masked, rv$boundary)
      showNotification("Building zone maps (kriging)…",
                       type = "message", duration = NULL, id = "zone_build")
      zm <- zoneMap()
      if (isTRUE(zm$ok)) {
        showNotification(paste0("Zone maps ready — ", zm$k, " zones"),
                         type = "message", duration = 4, id = "zone_build")
      } else {
        showNotification(paste0("Zone map error: ", zm$message),
                         type = "error", duration = 10, id = "zone_build")
      }
    },
    ignoreInit = TRUE)

  # Main hard-zone map.  terra::plot() is instant regardless of raster size
  # (ggplotly + a discrete geom_raster fill chokes on 1.5M cells and was the
  # original cause of the never-resolving spinner).
  output$plot_zone_map <- renderPlot({
    mzlog("RENDER plot_zone_map: enter")
    zm <- zoneMap()
    if (!isTRUE(zm$ok)) {
      mzlog("RENDER plot_zone_map: placeholder")
      return(drawPlaceholder(zm))
    }
    mzlog("RENDER plot_zone_map: drawing map")
    # No col= here: zone_hard_disp carries a value->colour table (coltab) set by
    # asCategoricalZones(), which pins each zone id to its colour.  Passing a
    # col= vector instead would let terra stretch it across the data range and
    # mis-colour zones when an id is missing.  The colour key is zone_summary_ui.
    terra::plot(zm$zone_hard_disp,
                main = paste0("Management Zone Map  (k = ", zm$k,
                              ", Fuzzy C-Means, m = ", isolate(input$fcm_m), ")"),
                axes = FALSE, box = FALSE, legend = FALSE)
    mzlog("RENDER plot_zone_map: done")
  })

  # Zone summary — % of points per zone (tabulate over 1..k so empty zones
  # still show 0%).  When the zone map isn't ready, show a friendly pointer
  # instead of an empty card body — req() halts silently and used to leave
  # the card blank the whole time the user was on Step 1–3.
  output$zone_summary_ui <- renderUI({
    zm <- zoneMap()
    mzlog("RENDER zone_summary_ui: ok=", isTRUE(zm$ok),
          " waiting=", isTRUE(zm$waiting),
          " k=", if (isTRUE(zm$ok)) zm$k else "—")
    if (!isTRUE(zm$ok)) {
      waiting <- isTRUE(zm$waiting)
      return(div(class = if (waiting) "info-msg" else "warning-msg",
        icon(if (waiting) "hourglass-half" else "exclamation-triangle"),
        " ", zm$message %||%
          "Zone map not ready yet — complete Steps 1–3."))
    }
    k   <- zm$k
    cnt <- tabulate(zm$cluster_id, nbins = k)
    pct <- round(100 * cnt / sum(cnt), 1)
    cols <- zoneColors(k)
    tagList(lapply(seq_len(k), function(z) {
      div(style = paste0("padding:8px 12px; border-radius:8px; margin-bottom:6px; background:",
            cols[z], "; color:white; font-weight:600; font-size:13px;"),
          paste0("Zone ", z, ": ", pct[z], "% of pixels"))
    }))
  })

  # Eager-render the Zone Summary just like the Zone Map (above).  Without
  # this, the uiOutput placeholder sits empty while the tab is hidden, and
  # the content it gets on un-suspend can land in a not-yet-laid-out
  # container — leaving the card body blank even though zoneMap() says ok
  # (which is exactly the state the user's screenshot shows: Zone Map
  # rendered, Zone Summary empty).  The eager render matches what we already
  # do for plot_zone_map / plot_membership, so all three Zone Maps outputs
  # paint together the moment the tab becomes visible.
  outputOptions(output, "zone_summary_ui", suspendWhenHidden = FALSE)

  # Membership surfaces — ALL k panels rendered as ONE static server-side
  # image.  No per-zone plotly widgets, so nothing heavy is shipped to the
  # browser and the panel always paints in one shot.
  output$plot_membership <- renderPlot({
    zm <- zoneMap()
    if (!isTRUE(zm$ok)) return(drawPlaceholder(zm))
    terra::plot(zm$memb_stack_disp,
                col  = grDevices::hcl.colors(50, "Plasma", rev = TRUE),
                nc   = ceiling(sqrt(zm$k)),
                axes = FALSE, box = FALSE,
                mar  = c(1.5, 1.5, 2.5, 3.5))
  })

  # Plain-text status so the user always sees the state in words, never just a
  # spinning loader with no explanation.
  output$zone_status_ui <- renderUI({
    zm <- zoneMap()
    if (isTRUE(zm$ok)) {
      div(class = "success-msg", icon("check"),
          sprintf(" Zone maps ready — %d zones.", zm$k))
    } else if (isTRUE(zm$waiting)) {
      div(class = "info-msg", icon("hourglass-half"), " ", zm$message)
    } else {
      div(class = "warning-msg", icon("exclamation-triangle"),
          paste0(" Zone map could not be generated: ", zm$message))
    }
  })

  # Render the Part 4 plots even while their tab is hidden, so the image is
  # already present (and correctly sized after the resize trigger) the instant
  # the user opens Zone Maps — no half-painted, never-clearing output.
  outputOptions(output, "plot_zone_map",   suspendWhenHidden = FALSE)
  outputOptions(output, "plot_membership", suspendWhenHidden = FALSE)

  # ── Statistics ─────────────────────────────────────────────────────────────
  # ARCHITECTURE NOTE — same fix as Part 4 (Zone Maps), applied here:
  # The previous version created output$table_zone_means / output$table_anova /
  # output$plot_zone_means as SIDE EFFECTS *inside* a plain observe() guarded by
  # req().  When req() halted (clustering/df not ready yet) the outputs were
  # never defined, so the DT tables and the plotly spinner had nothing to
  # resolve to and spun forever — the "loop sensation".  The export buttons
  # worked because each recomputes the stats itself in its own observeEvent.
  #
  # The fix mirrors zoneMap(): do the work ONCE in a cached reactive
  # (zoneStats) that NEVER halts — it returns a self-diagnosing list
  # {ok, waiting, message, ...} — and define every Statistics output at the TOP
  # LEVEL reading that reactive.  Now the outputs always exist, always resolve,
  # and any "not ready" / failure surfaces as a readable message.
  zoneStats <- reactive({
    miss <- c()
    if (is.null(rv$fcm_result)) miss <- c(miss, "clustering result (Step 3)")
    if (is.null(rv$df_vars))    miss <- c(miss, "soil variables (Step 1)")
    if (is.null(rv$soil_vars))  miss <- c(miss, "variable names (Step 1)")
    if (length(miss) > 0) {
      return(list(ok = FALSE, waiting = TRUE,
                  message = paste0("Waiting for: ", paste(miss, collapse = ", "),
                                   ". Complete that step, then return here.")))
    }

    tryCatch({
      soil_vars <- rv$soil_vars
      df_stats  <- cbind(rv$df_vars, Zone = rv$fcm_result$cluster_id)

      # Zone means
      zone_means <- aggregate(. ~ Zone, data = df_stats, FUN = mean)
      zone_means$Zone <- paste0("Zone ", zone_means$Zone)

      # ANOVA — wrap aov() in tryCatch so a degenerate model (e.g. zero
      # residual df, single-value group) returns NA instead of failing the
      # whole reactive.  sapply() with USE.NAMES = FALSE guarantees a numeric
      # vector of length(soil_vars) regardless of intermediate failures.
      anova_res <- sapply(soil_vars, function(v) {
        fit <- aov(as.formula(paste(v, "~ Zone")), data = df_stats)
        tryCatch(summary(fit)[[1]]$"Pr(>F)"[1], error = function(e) NA_real_)
      }, USE.NAMES = FALSE)

      anova_df <- data.frame(Variable = soil_vars, `p-value` = anova_res,
                             check.names = FALSE)
      anova_df$Significant <- ifelse(is.na(anova_df$`p-value`), "ns",
                            ifelse(anova_df$`p-value` < 0.001, "***",
                            ifelse(anova_df$`p-value` < 0.01, "**",
                            ifelse(anova_df$`p-value` < 0.05, "*", "ns"))))

      zone_long <- tidyr::gather(zone_means, key = "Variable", value = "Value", -Zone)

      list(ok = TRUE, zone_means = zone_means, anova_df = anova_df,
           zone_long = zone_long)
    }, error = function(e) {
      list(ok = FALSE, waiting = FALSE, message = conditionMessage(e))
    })
  })

  # Small placeholder data frame so a DT table renders a readable "not ready"
  # message instead of perpetually showing its built-in loading state.
  statsPlaceholderDT <- function(zs) {
    DT::datatable(
      data.frame(Status = zs$message %||%
                   "Statistics not ready yet — complete the earlier steps."),
      options = list(dom = "t", ordering = FALSE), rownames = FALSE,
      colnames = "")
  }

  output$table_zone_means <- DT::renderDataTable({
    zs <- zoneStats()
    if (!isTRUE(zs$ok)) return(statsPlaceholderDT(zs))
    # Round only the numeric columns — Zone is the character label
    # ("Zone 1", …), and round() on a data frame with a non-numeric column
    # errors with "non-numeric-alike variable(s) in data frame: Zone".
    zm <- zs$zone_means
    num_cols <- vapply(zm, is.numeric, logical(1))
    zm[num_cols] <- round(zm[num_cols], 3)
    DT::datatable(zm, options = list(pageLength = 10),
                  rownames = FALSE, caption = "Mean soil properties per zone")
  })

  output$table_anova <- DT::renderDataTable({
    zs <- zoneStats()
    if (!isTRUE(zs$ok)) return(statsPlaceholderDT(zs))
    DT::datatable(zs$anova_df, options = list(pageLength = 10),
                  rownames = FALSE, caption = "ANOVA p-values — significance of zone differences") %>%
      DT::formatSignif("p-value", 4)
  })

  output$plot_zone_means <- renderPlotly({
    zs <- zoneStats()
    if (!isTRUE(zs$ok)) {
      return(plotly_empty(type = "scatter", mode = "markers") %>%
               layout(title = list(text = zs$message %||%
                 "Statistics not ready yet — complete the earlier steps.")))
    }
    p <- ggplot(zs$zone_long, aes(x = Variable, y = Value, fill = Zone, group = Zone)) +
      geom_col(position = "dodge", width = 0.7) +
      labs(title = "Soil Properties by Zone",
           x = "Soil Variable", y = "Mean Value (normalized)") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", color = "#1a252f"),
            legend.position = "right")
    ggplotly(p) %>% layout(title = list(text = paste0(
      "Soil Properties by Zone<br><sup>Mean values per zone (normalized scale)</sup>")))
  })

  # Pre-render every output while hidden so it is resolved the instant the
  # tab opens.  Both DT tables and the plot need this — a DataTable (or
  # plotly widget) created while its container was display:none bakes in
  # width 0 and never reflows on its own.  The JS tab-show handler below
  # calls DataTable.columns.adjust() (in addition to the window resize) so
  # the eager-rendered tables get their column widths fixed when the tab
  # becomes visible.  The heavy work is the shared zoneStats() reactive,
  # cached across all three outputs — so making it eager costs one compute
  # per (re)clustering, not three.  (table_validation / table_assignments
  # dodge the same bug for a different reason — they are created inside
  # their button observeEvent while their own tab is already visible.)
  outputOptions(output, "plot_zone_means",    suspendWhenHidden = FALSE)
  outputOptions(output, "table_zone_means",   suspendWhenHidden = FALSE)
  outputOptions(output, "table_anova",        suspendWhenHidden = FALSE)

  # ── Export ─────────────────────────────────────────────────────────────────
  # All exports write to outputs/ (matching the batch script's layout). The
  # directory is created on first use so the app does not pollute the
  # project until the user actually exports something.

  ensure_outputs_dir <- function() {
    if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)
  }

  observeEvent(input$btn_export_validation, {
    req(rv$validation_df)
    ensure_outputs_dir()
    write.csv(rv$validation_df, "outputs/mz_validation.csv", row.names = FALSE)
    showNotification("Validation CSV saved → outputs/mz_validation.csv",
                     type = "message")
  })

  observeEvent(input$btn_export_zone_stats, {
    req(rv$fcm_result, rv$df_vars)
    ensure_outputs_dir()
    df_stats <- cbind(rv$df_vars, Zone = rv$fcm_result$cluster_id)
    zone_means <- aggregate(. ~ Zone, data = df_stats, FUN = mean)
    write.csv(zone_means, "outputs/mz_zone_stats.csv", row.names = FALSE)
    showNotification("Zone statistics CSV saved → outputs/mz_zone_stats.csv",
                     type = "message")
  })

  observeEvent(input$btn_export_assignments, {
    req(rv$obs_sf)
    ensure_outputs_dir()
    pa <- rv$obs_sf
    if (!"fid" %in% names(pa)) pa$fid <- seq_len(nrow(pa))
    sel <- intersect(c("fid", "Zone"), names(pa))
    assign_df <- as.data.frame(st_drop_geometry(pa[, sel]))
    write.csv(assign_df, "outputs/mz_point_assignments.csv", row.names = FALSE)
    showNotification("Point assignments CSV saved → outputs/mz_point_assignments.csv",
                     type = "message")
  })

  observeEvent(input$btn_export_anova, {
    req(rv$fcm_result, rv$df_vars, rv$soil_vars)
    ensure_outputs_dir()
    df_stats <- cbind(rv$df_vars, Zone = rv$fcm_result$cluster_id)
    anova_res <- lapply(rv$soil_vars, function(v) {
      summary(aov(as.formula(paste(v, "~ Zone")), data = df_stats))[[1]]$"Pr(>F)"[1]
    })
    anova_df <- data.frame(Variable = rv$soil_vars, `p-value` = unlist(anova_res))
    write.csv(anova_df, "outputs/mz_anova.csv", row.names = FALSE)
    showNotification("ANOVA CSV saved → outputs/mz_anova.csv",
                     type = "message")
  })

  # ── Map exports ──────────────────────────────────────────────────────────────
  # Both reuse rv$zone_hard, the full-resolution hard-zone raster cached by
  # zoneMap() (Part 4). req() keeps the handler inert until the zone map exists,
  # and any failure surfaces as an on-screen error notification rather than a
  # silent crash — same defensive style as the rest of the app.

  # GeoTIFF — full-resolution integer zone raster, CRS preserved for GIS.
  observeEvent(input$btn_export_zone_tif, {
    if (is.null(rv$zone_hard)) {
      showNotification("Generate the zone map (Step 4) before exporting the GeoTIFF.",
                       type = "warning", duration = 6)
      return(invisible(NULL))
    }
    tryCatch({
      ensure_outputs_dir()
      out <- "outputs/mz_zones.tif"
      terra::writeRaster(rv$zone_hard, out, overwrite = TRUE,
                         datatype = "INT1U", NAflag = 255)
      mzlog("EXPORT: wrote GeoTIFF ", out, " | cells=", terra::ncell(rv$zone_hard))
      showNotification(paste0("Zones GeoTIFF saved → ", out), type = "message",
                       duration = 5)
    }, error = function(e) {
      mzlog("EXPORT: GeoTIFF error — ", conditionMessage(e))
      showNotification(paste0("GeoTIFF export failed: ", conditionMessage(e)),
                       type = "error", duration = 10)
    })
  })

  # PNG — 300-DPI publication figure of the hard-zone map, with a zone legend.
  observeEvent(input$btn_export_zone_png, {
    if (is.null(rv$zone_hard)) {
      showNotification("Generate the zone map (Step 4) before exporting the PNG.",
                       type = "warning", duration = 6)
      return(invisible(NULL))
    }
    tryCatch({
      ensure_outputs_dir()
      out  <- "outputs/mz_zone_map.png"
      k    <- rv$k_final
      cols <- zoneColors(k)
      grDevices::png(out, width = 8, height = 6, units = "in", res = 300)
      on.exit(grDevices::dev.off(), add = TRUE)
      # No col= : rv$zone_hard carries the value->colour table (coltab) set by
      # asCategoricalZones(), so each zone id draws in its fixed colour.  The
      # manual legend uses zoneColors(k) by index, which now matches exactly.
      terra::plot(rv$zone_hard,
                  main = paste0("Management Zone Map  (k = ", k,
                                ", Fuzzy C-Means, m = ", isolate(input$fcm_m), ")"),
                  axes = FALSE, box = FALSE, legend = FALSE)
      graphics::legend("topright", legend = paste0("Zone ", seq_len(k)),
                       fill = cols, border = NA, bty = "n", cex = 0.9)
      mzlog("EXPORT: wrote PNG ", out, " @300dpi | k=", k)
      showNotification(paste0("Zone map PNG (300 DPI) saved → ", out),
                       type = "message", duration = 5)
    }, error = function(e) {
      mzlog("EXPORT: PNG error — ", conditionMessage(e))
      showNotification(paste0("PNG export failed: ", conditionMessage(e)),
                       type = "error", duration = 10)
    })
  })

  # ── HTML report ─────────────────────────────────────────────────────────────
  # Self-contained: writes everything the Quarto report expects to outputs/
  # (including the transition mask the rest of the app does not persist),
  # then calls R/render_mz_report.R to render report/mz_report.qmd to HTML.
  # Mirrors the standalone batch pipeline's output layout (data in outputs/,
  # qmd in report/), so the same report works from either entry point.
  observeEvent(input$btn_export_report, {
    need <- list(
      validation_df = rv$validation_df,
      fcm_result    = rv$fcm_result,
      df_vars       = rv$df_vars,
      obs_sf        = rv$obs_sf,
      zone_hard     = rv$zone_hard,
      soil_vars     = rv$soil_vars
    )
    missing_keys <- names(need)[vapply(need, is.null, logical(1))]
    if (length(missing_keys) > 0) {
      showNotification(
        paste0("Cannot generate report: run steps 1–4 first (missing: ",
               paste(missing_keys, collapse = ", "), ")."),
        type = "warning", duration = 8
      )
      return(invisible(NULL))
    }

    withProgress(message = "Rendering HTML report…", value = 0, {
      tryCatch({
        incProgress(0.1, detail = "Writing outputs/")
        ensure_outputs_dir()

        # 1. Validation table
        write.csv(rv$validation_df, "outputs/mz_validation.csv",
                  row.names = FALSE)

        # 2. Per-zone means
        df_stats   <- cbind(rv$df_vars, Zone = rv$fcm_result$cluster_id)
        zone_means <- aggregate(. ~ Zone, data = df_stats, FUN = mean)
        write.csv(zone_means, "outputs/mz_zone_stats.csv", row.names = FALSE)

        # 3. ANOVA (p-values only — the .qmd recomputes η² from raw data)
        anova_res <- lapply(rv$soil_vars, function(v) {
          summary(aov(as.formula(paste(v, "~ Zone")), data = df_stats))[[1]]$"Pr(>F)"[1]
        })
        write.csv(data.frame(variable = rv$soil_vars,
                             p_value  = unlist(anova_res)),
                  "outputs/mz_anova.csv", row.names = FALSE)

        # 4. Point assignments (rich: fid, lon, lat, Zone, memberships)
        pa <- rv$obs_sf
        if (!"fid" %in% names(pa)) pa$fid <- seq_len(nrow(pa))
        coords <- st_coordinates(pa)
        pa_tbl <- data.frame(
          fid       = pa$fid,
          longitude = coords[, 1],
          latitude  = coords[, 2],
          Zone      = rv$fcm_result$cluster_id,
          round(as.data.frame(rv$fcm_result$membership), 6)
        )
        names(pa_tbl)[6:ncol(pa_tbl)] <- paste0("membership_z", seq_len(ncol(rv$fcm_result$membership)))
        write.csv(pa_tbl, "outputs/mz_point_assignments.csv", row.names = FALSE)

        # 5. Hard zone map
        terra::writeRaster(rv$zone_hard, "outputs/mz_zone_map.tif",
                           overwrite = TRUE, datatype = "INT1U", NAflag = 255)

        # 5a. Per-zone membership rasters. The Quarto report's uncertainty
        #     section (Shannon entropy) reads mz_membership_z*.tif. The
        #     standalone batch script writes these from the kriging output;
        #     the Shiny app previously did not. Persist them now so the
        #     report has everything it needs without the user having to
        #     re-run the batch pipeline.
        if (!is.null(rv$zone_stack)) {
          k_layers <- terra::nlyr(rv$zone_stack)
          for (j in seq_len(k_layers)) {
            terra::writeRaster(rv$zone_stack[[j]],
                               sprintf("outputs/mz_membership_z%d.tif", j),
                               overwrite = TRUE, datatype = "FLT4S",
                               NAflag = -9999)
          }
        }

        # 5b. PNG figure of the zone map (300 DPI, with legend)
        png_path <- "outputs/mz_zone_map.png"
        k        <- rv$k_final
        cols     <- zoneColors(k)
        grDevices::png(png_path, width = 8, height = 6, units = "in", res = 300)
        on.exit(grDevices::dev.off(), add = TRUE)
        terra::plot(rv$zone_hard,
                    main = paste0("Management Zone Map  (k = ", k,
                                  ", Fuzzy C-Means, m = ",
                                  isolate(input$fcm_m), ")"),
                    axes = FALSE, box = FALSE, legend = FALSE)
        graphics::legend("topright",
                         legend = paste0("Zone ", seq_len(k)),
                         fill = cols, border = NA, bty = "n", cex = 0.9)

        # 5c. Validation index plot
        vi_path <- "outputs/mz_validation_all_indices.png"
        scale01 <- function(x) {
          rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
          if (!is.finite(rng) || rng == 0) return(rep(0.5, length(x)))
          (x - min(x, na.rm = TRUE)) / rng
        }
        vplot_df <- rbind(
          data.frame(k = rv$validation_df$k, index = "XB",
                     value = scale01(rv$validation_df$XB)),
          data.frame(k = rv$validation_df$k, index = "FPI",
                     value = scale01(rv$validation_df$FPI)),
          data.frame(k = rv$validation_df$k, index = "NCE",
                     value = scale01(rv$validation_df$NCE)),
          data.frame(k = rv$validation_df$k, index = "PE",
                     value = scale01(rv$validation_df$PE)),
          data.frame(k = rv$validation_df$k, index = "FS",
                     value = scale01(rv$validation_df$FS))
        )
        vplot <- ggplot2::ggplot(vplot_df, ggplot2::aes(k, value, color = index)) +
          ggplot2::geom_line(linewidth = 0.8) +
          ggplot2::geom_point(size = 2) +
          ggplot2::scale_x_continuous(breaks = rv$validation_df$k) +
          ggplot2::labs(x = "Number of zones (k)",
                        y = "Scaled index value (lower = better)",
                        color = "Index",
                        title = "Cluster-validity indices") +
          ggplot2::theme_minimal(base_size = 12)
        ggplot2::ggsave(vi_path, vplot, width = 7, height = 4.5, dpi = 300)

        # 6. Transition mask — the standalone script writes this, the app
        #    does not. Compute from rv$zone_stack (full-res membership stack)
        #    at the same default 0.60 threshold.
        if (!is.null(rv$zone_stack)) {
          trans_thr <- 0.60
          transition_mask <- terra::app(rv$zone_stack, function(x) {
            if (all(is.na(x))) return(NA_real_)
            as.numeric(max(x, na.rm = TRUE) < trans_thr)
          })
          names(transition_mask) <- "Transition"
          terra::writeRaster(transition_mask, "outputs/mz_transition.tif",
                             overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        }

        incProgress(0.5, detail = "Staging inputs & running Quarto CLI")

        # Stage the user's uploaded input files to a stable, extension-
        # preserving location so the Quarto report (a separate process) can
        # read them. Shiny uploads live in session temp files WITHOUT their
        # original extension, which can confuse format auto-detection in
        # sf/terra; copying them out with the original name removes that risk
        # and lets the report's Inputs/PCA/ANOVA/entropy sections reflect the
        # actual data instead of the bundled demo.
        stage_dir <- file.path(tempdir(), "mz_report_inputs")
        dir.create(stage_dir, showWarnings = FALSE, recursive = TRUE)
        stage_upload <- function(upload, fallback_name) {
          if (is.null(upload) || is.null(upload$datapath) ||
              !file.exists(upload$datapath)) return(NULL)
          nm <- if (!is.null(upload$name) && nzchar(upload$name)) {
            upload$name
          } else {
            fallback_name
          }
          dest <- file.path(stage_dir, nm)
          ok <- file.copy(upload$datapath, dest, overwrite = TRUE)
          if (isTRUE(ok)) normalizePath(dest, mustWork = TRUE) else NULL
        }
        staged_boundary <- isolate(stage_upload(input$input_boundary, "boundary.gpkg"))
        staged_raster   <- isolate(stage_upload(input$input_raster,   "soil_predictions.tif"))
        staged_points   <- isolate(stage_upload(input$input_points,   "soilgrids_data.csv"))

        # Resolve absolute paths up front so the wrapper does not depend on
        # CWD. The wrapper is sourced to the GLOBAL env so the function is
        # reliably found on the next line (sourcing with local = TRUE inside
        # withProgress can drop the binding in some Shiny versions).
        project_root <- normalizePath(getwd(), mustWork = TRUE)
        qmd_path  <- file.path(project_root, "report", "mz_report.qmd")
        wrap_path <- file.path(project_root, "R", "render_mz_report.R")
        out_html_path <- file.path(project_root, "report", "mz_report.html")

        if (!file.exists(qmd_path)) {
          stop("mz_report.qmd not found at ", qmd_path)
        }
        if (!file.exists(wrap_path)) {
          stop("R/render_mz_report.R not found at ", wrap_path)
        }

        # Capture pre-render mtime so we can verify the file was actually
        # updated afterwards.
        html_mtime_before <- if (file.exists(out_html_path)) {
          file.info(out_html_path)$mtime
        } else {
          NA
        }

        source(wrap_path, local = FALSE)
        mzlog("REPORT: wrapper sourced from ", wrap_path)

        out_html <- render_mz_report(
          outputs_dir          = "outputs",
          data_dir             = "data",
          boundary_path        = staged_boundary,
          raster_path          = staged_raster,
          points_path          = staged_points,
          soil_vars            = paste(rv$soil_vars, collapse = ","),
          author               = isolate(input$report_author %||% "MZ Analysis"),
          study_area_name      = isolate(input$report_study_area %||% "Study Area"),
          k_range              = seq(isolate(input$k_min), isolate(input$k_max)),
          fuzziness            = isolate(input$fcm_m),
          pca_threshold        = isolate(input$pca_thresh) / 100,
          transition_threshold = 0.60,
          interpolation_method = "kriging",
          qmd_path             = qmd_path,
          verbose              = FALSE
        )

        incProgress(0.9, detail = "Verifying output")
        if (!file.exists(out_html)) {
          stop("Render reported success but output file is missing: ", out_html)
        }
        html_mtime_after <- file.info(out_html)$mtime
        if (!is.na(html_mtime_before) &&
            !is.na(html_mtime_after) &&
            html_mtime_after <= html_mtime_before) {
          stop("Render finished but the HTML mtime did not change ",
               "(before=", format(html_mtime_before), ", after=",
               format(html_mtime_after), "). Check that Quarto wrote to ",
               out_html, " and that no read-only file is in the way.")
        }
        size_kb <- round(file.info(out_html)$size / 1024, 1)
        mzlog("REPORT: rendered ", out_html,
              " | size=", size_kb, " KB",
              " | mtime=", format(html_mtime_after))

        incProgress(1.0, detail = "Done")
        showNotification(
          paste0("HTML report rendered → ", out_html,
                 "  (", size_kb, " KB)"),
          type = "message", duration = 6
        )
      }, error = function(e) {
        mzlog("REPORT: error — ", conditionMessage(e))
        showNotification(
          paste0("Report rendering failed: ", conditionMessage(e)),
          type = "error", duration = 15
        )
      })
    })
  })

  output$session_summary <- renderUI({
    k <- rv$k_final
    n_pts <- if (!is.null(rv$obs_sf)) nrow(rv$obs_sf) else "—"
    n_vars <- if (!is.null(rv$soil_vars)) length(rv$soil_vars) else "—"
    method <- input$opt_method

    tagList(
      div(strong("Optimal k:"), " ", rv$selected_k),
      div(strong("Method:"), " ", method),
      div(strong("Points:"), " ", n_pts),
      div(strong("Variables:"), " ", n_vars),
      div(strong("FPI:"), " ", sprintf("%.4f", rv$validation_df$FPI[rv$validation_df$k == rv$selected_k][1])),
      div(strong("NCE:"), " ", sprintf("%.4f", rv$validation_df$NCE[rv$validation_df$k == rv$selected_k][1])),
      div(strong("Xie-Beni:"), " ", sprintf("%.4f", rv$validation_df$XB[rv$validation_df$k == rv$selected_k][1]))
    )
  })

  observeEvent(input$btn_restart, {
    session$reload()
  })
}

# =============================================================================
# RUN
# =============================================================================

shinyApp(ui = ui, server = server)