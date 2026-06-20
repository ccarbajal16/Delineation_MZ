# =============================================================================
# render_mz_report.R
# -----------------------------------------------------------------------------
# Wrapper around the Quarto CLI to render `report/mz_report.qmd` from R.
#
# Usage:
#   source("R/render_mz_report.R")
#   render_mz_report()                                         # default HTML
#   render_mz_report(format = "pdf", author = "My Name")
#   render_mz_report(outputs_dir = "outputs/run_2024_06_05")   # alt outputs
#
# Requires: Quarto CLI >= 1.3 (https://quarto.org/docs/get-started/).
# The `quarto` R package is NOT required; we shell out to the CLI directly.
# =============================================================================

#' Render the MZ delineation report
#'
#' @param outputs_dir Path to the directory containing the pipeline outputs
#'   (`mz_validation.csv`, `mz_zone_map.tif`, etc.). Default: `"outputs"`.
#' @param data_dir Path to the input data directory. Default: `"data"`. Used as a
#'   fallback for the input files when the explicit `*_path` arguments are not
#'   supplied (the batch-script entry point relies on this).
#' @param boundary_path Optional explicit path to the study-area boundary file.
#'   When supplied (e.g. by the Shiny app with the user's uploaded file), the
#'   report uses it instead of `<data_dir>/agro_geo.gpkg`. `NULL` keeps the
#'   batch-script behaviour.
#' @param raster_path Optional explicit path to the soil-property raster.
#'   Overrides `<data_dir>/soil_predictions.tif`. `NULL` keeps the default.
#' @param points_path Optional explicit path to the observation-points CSV.
#'   Overrides `<data_dir>/soilgrids_data.csv`. `NULL` keeps the default.
#' @param soil_vars Optional soil-variable names to analyze, as a character
#'   vector (or a single comma-separated string). When supplied, overrides the
#'   report's built-in default variable set so the report adapts to whatever
#'   columns the caller actually clustered on. `NULL` keeps the default.
#' @param format Output format. One of `"html"`, `"pdf"`. Default: `"html"`.
#' @param author Author name to print on the report cover.
#' @param study_area_name Short label for the study area (printed on cover).
#' @param k_range Integer vector of k values evaluated by the pipeline.
#' @param fuzziness FCM fuzziness exponent (m).
#' @param pca_threshold Cumulative-variance threshold for PCA retention.
#' @param transition_threshold Max-membership threshold for the transition mask.
#' @param interpolation_method One of `"kriging"`, `"idw"`.
#' @param seed Random seed used by the pipeline.
#' @param qmd_path Path to the Quarto source file. If `NULL` (default), the
#'   function locates `report/mz_report.qmd` relative to the project root.
#' @param output_file Output filename (without extension). If `NULL`, defaults
#'   to `mz_report.<format>` written next to the `.qmd`.
#' @param output_dir Directory to write the rendered file to. If `NULL`, the
#'   file is written next to the `.qmd` (default Quarto behavior).
#' @param quarto_bin Path to the Quarto executable. Default: `"quarto"`.
#' @param verbose If `TRUE`, stream Quarto's stdout/stderr to the R console.
#'
#' @return Invisibly returns the absolute path of the rendered file.
#' @export
render_mz_report <- function(
  outputs_dir = "outputs",
  data_dir = "data",
  boundary_path = NULL,
  raster_path = NULL,
  points_path = NULL,
  soil_vars = NULL,
  format = c("html", "pdf"),
  author = "MZ Analysis",
  study_area_name = "Study Area",
  k_range = 2:6,
  fuzziness = 2,
  pca_threshold = 0.80,
  transition_threshold = 0.60,
  interpolation_method = c("kriging", "idw"),
  seed = 42,
  qmd_path = NULL,
  output_file = NULL,
  output_dir = NULL,
  quarto_bin = "quarto",
  verbose = FALSE
) {
  format <- match.arg(format)
  interpolation_method <- match.arg(interpolation_method)

  # ---- Locate the .qmd if not provided -------------------------------------
  if (is.null(qmd_path)) {
    # Caller of this function is the project root (heuristic: walk up to find
    # `delineation_management_zones.R`).
    root <- try_find_project_root()
    if (is.null(root)) {
      stop("Could not locate the project root automatically. ",
           "Pass `qmd_path` explicitly.")
    }
    qmd_path <- file.path(root, "report", "mz_report.qmd")
  }
  if (!file.exists(qmd_path)) {
    stop("Quarto source not found: ", qmd_path)
  }

  # ---- Resolve absolute paths ---------------------------------------------
  qmd_dir <- normalizePath(dirname(qmd_path), mustWork = TRUE)
  qmd_path <- normalizePath(qmd_path, mustWork = TRUE)
  outputs_dir <- normalizePath(outputs_dir, mustWork = TRUE)
  data_dir <- normalizePath(data_dir, mustWork = TRUE)

  if (!dir.exists(outputs_dir)) {
    stop("outputs_dir does not exist: ", outputs_dir)
  }
  if (!dir.exists(data_dir)) {
    stop("data_dir does not exist: ", data_dir)
  }

  # ---- Check Quarto CLI is available --------------------------------------
  # Resolve the absolute path via Sys.which() so the call works even when
  # PATH inside the calling R session is sparser than the user's shell
  # (e.g. Shiny apps launched from a launcher that scrubs PATH).
  if (!file.exists(quarto_bin)) {
    resolved <- Sys.which("quarto")
    if (nzchar(resolved) && file.exists(resolved)) {
      quarto_bin <- resolved
    }
  }
  quarto_version <- tryCatch(
    system2(quarto_bin, "--version", stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  if (is.na(quarto_version) || !nzchar(quarto_version)) {
    stop("Quarto CLI not found on PATH (binary: '", quarto_bin, "'). ",
         "Install from https://quarto.org/docs/get-started/.")
  }

  # ---- Build -P parameter flags -------------------------------------------
  # We use repeated `-P key:value` flags rather than `--params-file`; the
  # latter is forwarded to pandoc which does not recognize it.
  p_flag <- function(key, val) {
    if (is.character(val)) {
      sprintf("-P %s:%s", key, shQuote(val))
    } else if (is.numeric(val) && length(val) > 1L) {
      arr <- paste(format(val, trim = TRUE, scientific = FALSE), collapse = ",")
      # Wrap list value in double quotes; Quarto's CLI splits on whitespace
      # so spaces inside the brackets break the YAML parser.
      sprintf('-P %s:"[%s]"', key, arr)
    } else if (is.numeric(val)) {
      sprintf("-P %s:%s", key, format(val, trim = TRUE, scientific = FALSE))
    } else {
      sprintf("-P %s:%s", key, as.character(val))
    }
  }

  # Optional input-path / variable params. Only emitted when the caller
  # supplies them, so the batch-script entry point (which relies on the demo
  # files living in data/ and the built-in variable set) keeps working
  # unchanged. soil_vars is collapsed to a comma-separated string to dodge
  # YAML/JSON array-quoting headaches on the command line; the .qmd splits it.
  opt_flag <- function(key, val) {
    if (is.null(val)) return(NULL)
    if (is.character(val) && length(val) == 1L && !nzchar(val)) return(NULL)
    if (identical(key, "soil_vars") && is.character(val) && length(val) > 1L) {
      val <- paste(val, collapse = ",")
    }
    p_flag(key, val)
  }

  param_flags <- c(
    p_flag("outputs_dir",          outputs_dir),
    p_flag("data_dir",             data_dir),
    opt_flag("boundary_path",      boundary_path),
    opt_flag("raster_path",        raster_path),
    opt_flag("points_path",        points_path),
    opt_flag("soil_vars",          soil_vars),
    p_flag("author",               author),
    p_flag("study_area_name",      study_area_name),
    p_flag("k_range",              as.integer(k_range)),
    p_flag("fuzziness",            fuzziness),
    p_flag("pca_threshold",        pca_threshold),
    p_flag("transition_threshold", transition_threshold),
    p_flag("interpolation_method", interpolation_method),
    p_flag("seed",                 seed)
  )

  # ---- Build the render command -------------------------------------------
  # By default, Quarto writes `<qmd-basename>.<format>` next to the .qmd.
  # We only pass `--output-dir` when the caller asked for a different one,
  # and we always resolve it to an absolute path. Quarto treats relative
  # `--output-dir` as relative to the .qmd's directory, which gives
  # surprising nesting (e.g. `report/outputs/...`).

  args <- c(
    "render", shQuote(qmd_path),
    "--to", format
  )
  if (!is.null(output_dir)) {
    out_dir <- normalizePath(output_dir, mustWork = FALSE)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    args <- c(args, "--output-dir", shQuote(out_dir))
  }

  # `output_file` is supported for compatibility but disabled by default.
  # Quarto's `--output NAME` (no extension) writes a file with that exact
  # name and a sibling `<NAME>_files/` resource dir, which is fragile.
  # We intentionally do NOT pass `--output`; use `output_dir` instead.
  if (!is.null(output_file)) {
    warning("`output_file` is ignored; use `output_dir` to relocate output.",
            call. = FALSE)
  }

  args <- c(args, param_flags)
  # NOTE: we deliberately do NOT pass --quiet. With --quiet, Quarto also
  # suppresses knitr's error traceback on stderr, which would make every
  # render failure invisible (the wrapper would only see a non-zero exit
  # status with no clue why). Output is captured to temp files below, so the
  # console stays clean regardless of verbosity.

  message("Rendering MZ report (format = ", format, ") ...")
  if (verbose) {
    message("Command: ", quarto_bin, " ",
            paste(args, collapse = " "))
  }

  # Capture BOTH channels to temp files. Quarto writes its progress to stderr
  # and knitr's error traceback (the bit we need on failure) to stderr as well.
  # The files are removed on exit; on failure their contents are folded into
  # the error message so the caller sees the real cause.
  stdout_log <- tempfile(fileext = "_quarto_stdout.log")
  stderr_log <- tempfile(fileext = "_quarto_stderr.log")
  on.exit({
    if (file.exists(stdout_log)) file.remove(stdout_log)
    if (file.exists(stderr_log)) file.remove(stderr_log)
  }, add = TRUE)

  result <- system2(quarto_bin, args, stdout = stdout_log, stderr = stderr_log)
  # system2 returns the exit status as the value directly when stdout/stderr
  # are file paths (no "status" attribute); it is carried as an attribute
  # only when stdout = TRUE. Handle both so a Quarto failure is always caught
  # instead of falling through to the (misleading) mtime check.
  status <- attr(result, "status")
  if (is.null(status)) status <- suppressWarnings(as.integer(result[1L]))
  read_log <- function(f) {
    if (file.exists(f) && file.size(f) > 0)
      paste(readLines(f, warn = FALSE), collapse = "\n")
    else ""
  }
  if (!is.na(status) && status != 0L) {
    combined <- trimws(paste(c(read_log(stdout_log), read_log(stderr_log)),
                             collapse = "\n"))
    if (!nzchar(combined)) {
      combined <- "(no output captured — re-run with verbose = TRUE for the command)"
    }
    stop("Quarto render failed with status ", status, ".\n",
         "--- quarto output ---\n", combined, "\n--- end output ---")
  }
  if (verbose) {
    out_txt <- read_log(stdout_log)
    if (nzchar(out_txt)) message(out_txt)
  }

  rendered <- if (is.null(output_dir)) {
    file.path(qmd_dir, paste0(tools::file_path_sans_ext(basename(qmd_path)),
                              ".", format))
  } else {
    out_dir <- normalizePath(output_dir, mustWork = FALSE)
    file.path(out_dir,
              paste0(tools::file_path_sans_ext(basename(qmd_path)), ".",
                     format))
  }
  if (!file.exists(rendered)) {
    stop("Render reported success but the output file is missing: ", rendered)
  }
  message("Done: ", rendered)
  invisible(normalizePath(rendered, mustWork = TRUE))
}

# Suppress "unused argument" notes for the no-longer-used helper.
# (kept as a comment marker — write_params_yaml was removed)

# ---- Internal helpers -------------------------------------------------------

#' Find the project root by looking upward for `delineation_management_zones.R`.
#' @return Absolute path or `NULL` if not found.
#' @keywords internal
try_find_project_root <- function() {
  d <- getwd()
  for (i in 0:8) {
    if (file.exists(file.path(d, "delineation_management_zones.R"))) {
      return(normalizePath(d, mustWork = TRUE))
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  NULL
}
