# =============================================================================
# Management Zone Delimitation - Fuzzy C-Means Clustering
# =============================================================================
# Purpose:
#   Complete standalone workflow for delimiting agricultural management zones
#   from soil-property raster predictions and point observations.
#
# Inputs (data/):
#   agro_geo.gpkg              - Study-area boundary polygon
#   soil_predictions.tif        - Multi-layer raster of soil properties
#   soilgrids_data.csv        - Observation points with soil variables
#
# Outputs (outputs/):
#   Raster:  mz_zone_map.tif, mz_transition.tif, mz_membership_z*.tif
#   Tables:  mz_validation.csv, mz_zone_stats.csv, mz_anova.csv,
#            mz_point_assignments.csv
#   Figures: mz_validation_all_indices.png, mz_zone_map.png,
#            mz_pca_scores_by_zone.png
# =============================================================================

suppressPackageStartupMessages({
  library(e1071)
  library(gstat)
  library(sf)
  library(terra)
  library(sp)
  library(ggplot2)
  library(patchwork)
})

HAS_GGREPEL <- requireNamespace("ggrepel", quietly = TRUE)

# ---- Configuration -----------------------------------------------------------

DEFAULT_BASE_DIR <- "Delineation_MZ"
BASE_DIR <- if (dir.exists(DEFAULT_BASE_DIR)) DEFAULT_BASE_DIR else getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
OUT_DIR  <- file.path(BASE_DIR, "outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BOUNDARY_PATH <- file.path(DATA_DIR, "agro_geo.gpkg")
RASTER_PATH   <- file.path(DATA_DIR, "soil_predictions.tif")
POINTS_PATH   <- file.path(DATA_DIR, "soilgrids_data.csv")

SOIL_VARS <- c("BD", "CEC", "Fragm", "Sand", "Silt", "Clay",
               "N", "OCD", "pH", "SOC")
K_RANGE <- 2:6
PCA_THRESHOLD <- 0.80
FUZZINESS <- 2
TRANSITION_THRESHOLD <- 0.60

INTERPOLATION_METHOD <- tolower(Sys.getenv("MZ_INTERPOLATION_METHOD",
                                           unset = "kriging"))
if (!INTERPOLATION_METHOD %in% c("kriging", "idw")) {
  stop("MZ_INTERPOLATION_METHOD must be 'kriging' or 'idw'.")
}

set.seed(42)

# ---- Helper functions --------------------------------------------------------

normalize <- function(x) {
  rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (!is.finite(rng) || rng == 0) return(rep(0.5, length(x)))
  (x - min(x, na.rm = TRUE)) / rng
}

scale01 <- function(x) {
  rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (!is.finite(rng) || rng == 0) return(rep(0.5, length(x)))
  (x - min(x, na.rm = TRUE)) / rng
}

xie_beni <- function(X, U, centroids, m = 2) {
  n <- nrow(X)
  k <- ncol(U)
  d <- ncol(X)
  Um <- U^m

  dist_sq <- vapply(seq_len(k), function(j) {
    rowSums((X - matrix(centroids[j, ], nrow = n, ncol = d,
                        byrow = TRUE))^2)
  }, numeric(n))

  jm <- sum(Um * dist_sq)
  cdist_sq <- as.matrix(dist(centroids))^2
  diag(cdist_sq) <- NA
  min_dist_sq <- min(cdist_sq, na.rm = TRUE)
  jm / (n * min_dist_sq)
}

partition_entropy <- function(U) {
  -sum(pmax(U, .Machine$double.xmin) *
         log(pmax(U, .Machine$double.xmin))) / nrow(U)
}

fpi_index <- function(U) {
  k <- ncol(U)
  pc <- sum(U^2) / nrow(U)
  1 - ((k * pc) - 1) / (k - 1)
}

nce_index <- function(U) {
  partition_entropy(U) / log(ncol(U))
}

fs_index <- function(X, U, centroids, m = 2) {
  n <- nrow(X)
  k <- ncol(U)
  d <- ncol(X)
  Um <- U^m

  dist_sq <- vapply(seq_len(k), function(j) {
    rowSums((X - matrix(centroids[j, ], nrow = n, ncol = d,
                        byrow = TRUE))^2)
  }, numeric(n))

  term1 <- sum(Um * dist_sq)
  w <- colSums(Um)
  cbar <- colSums(t(centroids) * w) / sum(w)
  dist_cbar_sq <- rowSums((centroids - matrix(cbar, nrow = k, ncol = d,
                                              byrow = TRUE))^2)
  term1 - sum(w * dist_cbar_sq)
}

median_impute <- function(df) {
  for (nm in names(df)) {
    if (anyNA(df[[nm]])) {
      med <- median(df[[nm]], na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      df[[nm]][is.na(df[[nm]])] <- med
      cat("  Imputed ", nm, " with median ", round(med, 4), "\n", sep = "")
    }
  }
  df
}

fit_cmeans <- function(scores, k, m = 2) {
  set.seed(42)
  cmeans(scores, centers = k, m = m, method = "cmeans",
         iter.max = 300, verbose = FALSE)
}

evaluate_k <- function(scores, k_range, m = 2) {
  out <- lapply(k_range, function(k) {
    fcm <- fit_cmeans(scores, k, m = m)
    U <- fcm$membership
    C <- fcm$centers
    data.frame(
      k = k,
      XB = xie_beni(scores, U, C, m = m),
      FPI = fpi_index(U),
      NCE = nce_index(U),
      PE = partition_entropy(U),
      FS = fs_index(scores, U, C, m = m)
    )
  })

  validation <- do.call(rbind, out)
  validation$XB_rank <- rank(validation$XB, ties.method = "min")
  validation$FPI_rank <- rank(validation$FPI, ties.method = "min")
  validation$NCE_rank <- rank(validation$NCE, ties.method = "min")
  validation$consensus_rank <- validation$XB_rank +
    validation$FPI_rank + validation$NCE_rank
  validation
}

select_optimal_k <- function(validation) {
  validation$k[which.min(validation$consensus_rank)]
}

interpolate_membership <- function(obs_sf, membership, ref_rast, boundary_vect,
                                   method = "kriging") {
  valid_cells <- which(!is.na(values(ref_rast, mat = FALSE)))
  grid_xy <- as.data.frame(xyFromCell(ref_rast, valid_cells))
  names(grid_xy) <- c("X", "Y")
  coordinates(grid_xy) <- ~ X + Y
  proj4string(grid_xy) <- CRS(crs(ref_rast, proj = TRUE))

  layers <- vector("list", ncol(membership))

  for (j in seq_len(ncol(membership))) {
    cat("  Interpolating membership surface z", j, "\n", sep = "")

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

    vals <- rep(NA_real_, ncell(ref_rast))
    vals[valid_cells] <- pmin(pmax(predicted, 0), 1)

    r <- setValues(ref_rast, vals)
    r <- mask(r, boundary_vect)
    names(r) <- paste0("Memb_z", j)
    layers[[j]] <- r
  }

  rast(layers)
}

save_validation_plot <- function(validation, path) {
  plot_df <- rbind(
    data.frame(k = validation$k, index = "XB", value = scale01(validation$XB)),
    data.frame(k = validation$k, index = "FPI", value = scale01(validation$FPI)),
    data.frame(k = validation$k, index = "NCE", value = scale01(validation$NCE)),
    data.frame(k = validation$k, index = "PE", value = scale01(validation$PE)),
    data.frame(k = validation$k, index = "FS", value = scale01(validation$FS))
  )

  p <- ggplot(plot_df, aes(k, value, color = index)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = validation$k) +
    labs(x = "Number of zones (k)",
         y = "Scaled index value",
         color = "Index",
         title = "Cluster-validity indices") +
    theme_minimal(base_size = 12)

  ggsave(path, p, width = 7, height = 4.5, dpi = 300)
}

save_zone_plot <- function(zone_raster, path) {
  plot_raster <- zone_raster
  if (ncell(plot_raster) > 1000000) {
    fact <- ceiling(sqrt(ncell(plot_raster) / 1000000))
    plot_raster <- aggregate(plot_raster, fact = fact, fun = "modal",
                             na.rm = TRUE)
  }

  zone_df <- as.data.frame(plot_raster, xy = TRUE, na.rm = TRUE)
  names(zone_df) <- c("x", "y", "Zone")

  p <- ggplot(zone_df, aes(x, y, fill = factor(Zone))) +
    geom_raster() +
    coord_equal() +
    labs(fill = "Zone", title = "Management zone map") +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_blank())

  ggsave(path, p, width = 7, height = 7, dpi = 300)
}

# ---- 1. Load and align spatial data -----------------------------------------

cat("=== Step 1: Loading inputs ===\n")
boundary <- st_read(BOUNDARY_PATH, quiet = TRUE)
boundary <- st_make_valid(boundary)
soil_stack <- rast(RASTER_PATH)
points_raw <- read.csv(POINTS_PATH, stringsAsFactors = FALSE)

missing_cols <- setdiff(c("fid", "longitude", "latitude", SOIL_VARS),
                        names(points_raw))
if (length(missing_cols) > 0) {
  stop("Point CSV is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

cat("  Boundary features: ", nrow(boundary), "\n", sep = "")
cat("  Raster layers: ", paste(names(soil_stack), collapse = ", "), "\n", sep = "")
cat("  Point rows: ", nrow(points_raw), "\n", sep = "")

raster_crs <- crs(soil_stack)
if (!identical(as.character(st_crs(boundary)), as.character(raster_crs))) {
  boundary <- st_transform(boundary, raster_crs)
}
boundary_vect <- vect(boundary)
soil_masked <- crop(soil_stack, boundary_vect)
soil_masked <- mask(soil_masked, boundary_vect)

points_sf <- st_as_sf(points_raw, coords = c("longitude", "latitude"),
                      crs = 4326)
points_sf <- st_transform(points_sf, raster_crs)

# ---- 2. Extract raster values and build feature matrix -----------------------

cat("\n=== Step 2: Building feature matrix ===\n")
extracted <- terra::extract(soil_masked, vect(points_sf), ID = FALSE)
df_vars <- points_raw[, SOIL_VARS, drop = FALSE]

for (var in intersect(SOIL_VARS, names(extracted))) {
  use_raster_value <- is.finite(extracted[[var]])
  df_vars[[var]][use_raster_value] <- extracted[[var]][use_raster_value]
}

df_vars <- as.data.frame(lapply(df_vars, as.numeric))
df_vars <- median_impute(df_vars)
df_norm <- as.data.frame(lapply(df_vars, normalize))
stopifnot(!anyNA(df_norm))

cat("  Features: ", nrow(df_norm), " rows x ", ncol(df_norm), " columns\n",
    sep = "")

# ---- 3. PCA ------------------------------------------------------------------

cat("\n=== Step 3: PCA ===\n")
pca_raw <- prcomp(df_norm, center = FALSE, scale. = FALSE)
eig <- pca_raw$sdev^2
prop_var <- eig / sum(eig)
cumvar <- cumsum(prop_var)
n_pc <- which(cumvar >= PCA_THRESHOLD)[1]
if (is.na(n_pc)) n_pc <- length(eig)
pca_scores <- pca_raw$x[, seq_len(n_pc), drop = FALSE]

cat("  Retained PCs: ", n_pc, "\n", sep = "")
cat("  Cumulative variance: ", round(cumvar[n_pc] * 100, 1), "%\n", sep = "")

# ---- 4. Fuzzy C-means validation and final model -----------------------------

cat("\n=== Step 4: Fuzzy C-means validation ===\n")
validation <- evaluate_k(pca_scores, K_RANGE, m = FUZZINESS)
write.csv(validation, file.path(OUT_DIR, "mz_validation.csv"),
          row.names = FALSE)
print(validation[, c("k", "XB", "FPI", "NCE", "PE", "FS",
                     "consensus_rank")])

optimal_k <- select_optimal_k(validation)
cat("  Selected k: ", optimal_k, "\n", sep = "")

final_fcm <- fit_cmeans(pca_scores, optimal_k, m = FUZZINESS)
membership <- final_fcm$membership
zone_id <- max.col(membership)

point_assignments <- cbind(
  points_raw[, c("fid", "longitude", "latitude")],
  Zone = zone_id,
  round(as.data.frame(membership), 6)
)
names(point_assignments)[5:ncol(point_assignments)] <-
  paste0("membership_z", seq_len(optimal_k))
write.csv(point_assignments, file.path(OUT_DIR, "mz_point_assignments.csv"),
          row.names = FALSE)

# ---- 5. Membership surfaces and zone rasters --------------------------------

cat("\n=== Step 5: Interpolating membership surfaces ===\n")
ref_rast <- soil_masked[[1]]
zone_stack <- interpolate_membership(
  obs_sf = points_sf,
  membership = membership,
  ref_rast = ref_rast,
  boundary_vect = boundary_vect,
  method = INTERPOLATION_METHOD
)

zone_hard <- app(zone_stack, function(x) {
  if (all(is.na(x))) return(NA_real_)
  which.max(x)
})
names(zone_hard) <- "Zone"

transition_mask <- app(zone_stack, function(x) {
  if (all(is.na(x))) return(NA_real_)
  as.numeric(max(x, na.rm = TRUE) < TRANSITION_THRESHOLD)
})
names(transition_mask) <- "Transition"

cat("\n=== Step 6: Writing raster outputs ===\n")
writeRaster(zone_hard, file.path(OUT_DIR, "mz_zone_map.tif"),
            overwrite = TRUE, datatype = "INT1U", NAflag = 255)
writeRaster(transition_mask, file.path(OUT_DIR, "mz_transition.tif"),
            overwrite = TRUE, datatype = "INT1U", NAflag = 255)

for (j in seq_len(nlyr(zone_stack))) {
  writeRaster(zone_stack[[j]],
              file.path(OUT_DIR, paste0("mz_membership_z", j, ".tif")),
              overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
}

# ---- 6. Statistics and figures ----------------------------------------------

cat("\n=== Step 7: Zone statistics and figures ===\n")
zone_stats <- aggregate(df_vars, by = list(Zone = zone_id), FUN = mean)
write.csv(zone_stats, file.path(OUT_DIR, "mz_zone_stats.csv"),
          row.names = FALSE)

anova_results <- do.call(rbind, lapply(SOIL_VARS, function(var) {
  fit <- aov(df_vars[[var]] ~ factor(zone_id))
  p_val <- summary(fit)[[1]][["Pr(>F)"]][1]
  data.frame(variable = var, p_value = p_val)
}))
write.csv(anova_results, file.path(OUT_DIR, "mz_anova.csv"),
          row.names = FALSE)

save_validation_plot(validation,
                     file.path(OUT_DIR, "mz_validation_all_indices.png"))
save_zone_plot(zone_hard, file.path(OUT_DIR, "mz_zone_map.png"))

scores_df <- as.data.frame(pca_raw$x[, 1:min(2, ncol(pca_raw$x)), drop = FALSE])
scores_df$Zone <- factor(zone_id)
if (ncol(scores_df) >= 3) {
  p_scores <- ggplot(scores_df, aes(PC1, PC2, color = Zone)) +
    geom_point(size = 2.5) +
    labs(title = "PCA scores by management zone") +
    theme_minimal(base_size = 12)
  ggsave(file.path(OUT_DIR, "mz_pca_scores_by_zone.png"),
         p_scores, width = 6, height = 5, dpi = 300)
}

cat("\nDone.\n")
cat("  Selected zones: ", optimal_k, "\n", sep = "")
cat("  Outputs written to: ", OUT_DIR, "\n", sep = "")
