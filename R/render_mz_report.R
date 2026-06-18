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
#' @param data_dir Path to the input data directory. Default: `"data"`.
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

  param_flags <- c(
    p_flag("outputs_dir",          outputs_dir),
    p_flag("data_dir",             data_dir),
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
  if (!verbose) args <- c(args, "--quiet")

  message("Rendering MZ report (format = ", format, ") ...")
  if (verbose) {
    message("Command: ", quarto_bin, " ",
            paste(args, collapse = " "))
  }

  # Always capture stderr to a temp file so failures are diagnosable even
  # when the caller asked for quiet output. The file is removed on success.
  stderr_log <- tempfile(fileext = "_quarto_stderr.log")
  on.exit(if (file.exists(stderr_log)) file.remove(stderr_log), add = TRUE)

  result <- system2(quarto_bin, args, stdout = "", stderr = stderr_log)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) {
    err_msg <- if (file.exists(stderr_log) && file.size(stderr_log) > 0) {
      paste(readLines(stderr_log, warn = FALSE), collapse = "\n")
    } else {
      "(no stderr captured — re-run with verbose = TRUE for the command)"
    }
    stop("Quarto render failed with status ", status, ".\n",
         "--- stderr ---\n", err_msg, "\n--- end stderr ---")
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
