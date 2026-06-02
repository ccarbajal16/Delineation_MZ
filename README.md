# Management Zone Delimitation — Fuzzy C-Means Clustering

A standalone R workflow that delineates agricultural management zones from
multi-layer soil-property rasters and point observations using Fuzzy C-Means
(FCM) clustering, PCA dimensionality reduction, and spatial interpolation.

![Management zone map](outputs/mz_zone_map.png)

---

## Overview

![Workflow diagram](figures/standalone_fig1_batch_workflow.png)

The script reads three local data sources, builds a normalized soil-property
matrix, reduces dimensionality with PCA, evaluates k = 2–6 candidate zone
counts with five cluster-validity indices, fits the final FCM model, and
interpolates per-zone membership surfaces across the study area. All results
are written as reproducible GeoTIFF, CSV, and PNG artifacts.

---

## Repository layout

```
Delineation_MZ/
├── delineation_management_zones.R   # Main script
├── data/
│   ├── agro_geo.gpkg                # Study-area boundary polygon
│   ├── soil_predictions.tif         # Multi-layer soil-property raster
│   └── soilgrids_data.csv           # Point observations with soil variables
├── figures/                         # Workflow diagrams (this README)
└── outputs/                         # All script outputs (generated)
```

---

## Requirements

Install dependencies from an R console before running the script:

```r
install.packages(c("e1071", "gstat", "sf", "terra", "sp",
                   "ggplot2", "patchwork"))
# optional but recommended for labelled plots:
install.packages("ggrepel")
```

| Package | Role |
|---------|------|
| `terra` / `sf` | Raster and vector I/O, reprojection, masking |
| `e1071` | Fuzzy C-Means via `cmeans()` |
| `gstat` | Ordinary kriging and IDW interpolation |
| `ggplot2` / `patchwork` | Validation and zone map figures |

---

## Running the script

```bash
Rscript delineation_management_zones.R
```

The interpolation method defaults to **ordinary kriging**. Override with an
environment variable:

```bash
MZ_INTERPOLATION_METHOD=idw Rscript delineation_management_zones.R
```

---

## Workflow steps

### Step 1 — Load and align spatial data

The boundary polygon is read from `agro_geo.gpkg` and reprojected to match
the raster CRS. The raster stack (`soil_predictions.tif`) is cropped and
masked to the boundary. Point observations are loaded from
`soilgrids_data.csv` and transformed to the same CRS.

Soil variables used:

```
BD  CEC  Fragm  Sand  Silt  Clay  N  OCD  pH  SOC
```

---

### Step 2 — Build the feature matrix

Raster values are extracted at each observation point with `terra::extract`.
Where a valid raster value exists it overwrites the CSV value, so raster
predictions take precedence. Missing values are imputed with the per-variable
median. All variables are then normalized to [0, 1].

---

### Step 3 — PCA

`prcomp()` is applied to the normalized matrix. Principal components are
retained until cumulative explained variance reaches **80 %** (configurable
via `PCA_THRESHOLD`). The resulting score matrix is the input to FCM.

---

### Step 4 — Fuzzy C-Means validation and zone selection

![Validation logic](figures/standalone_fig2_validation_logic.png)

FCM (`e1071::cmeans`, fuzziness m = 2, 300 iterations) is fitted for each
candidate k in {2, 3, 4, 5, 6}. Five validity indices are computed:

| Index | Formula / description | Direction |
|-------|-----------------------|-----------|
| **XB** (Xie–Beni) | Intra-cluster compactness / inter-cluster separation | lower is better |
| **FPI** (Fuzzy Partition Index) | `1 − ((k·PC − 1)/(k − 1))` | lower is better |
| **NCE** (Normalized Class Entropy) | `PE / log(k)` | lower is better |
| **PE** (Partition Entropy) | Mean fuzzy entropy | reported |
| **FS** (Fukuyama–Sugeno) | Compactness minus scatter | reported |

The **consensus rank** (XB_rank + FPI_rank + NCE_rank) selects the optimal k.
Results are saved to `outputs/mz_validation.csv`:

| k | XB | FPI | NCE | PE | FS | consensus_rank |
|---|-----|-----|-----|----|----|----------------|
| 2 | 0.104 | 0.375 | 0.438 | 0.304 | −46.3 | 14 |
| 3 | 0.092 | 0.295 | 0.329 | 0.361 | −82.7 | 9 |
| 4 | 0.096 | 0.300 | 0.309 | 0.428 | −97.8 | 10 |
| **5** | **0.062** | **0.270** | **0.272** | **0.437** | **−107.4** | **5** |
| 6 | 0.110 | 0.262 | 0.252 | 0.452 | −125.1 | 7 |

**Selected k = 5** (lowest consensus rank).

The validation index plot is saved as
`outputs/mz_validation_all_indices.png`.

---

### Step 5 — Final FCM model

FCM is re-fitted with the selected k = 5. This yields:

- **U** — n × 5 membership matrix (each row sums to 1).
- **C** — 5 × p centroid matrix.
- **zone_id** — hard zone label per point (`max.col(U)`).

Point-level results are written to `outputs/mz_point_assignments.csv`
(columns: `fid`, `longitude`, `latitude`, `Zone`,
`membership_z1` … `membership_z5`).

---

### Step 6 — Interpolate membership surfaces

![Interpolation outputs](figures/standalone_fig3_interpolation_outputs.png)

For each of the k zones, the membership column is spatially interpolated
across all valid raster cells:

1. **Ordinary kriging** — variogram fitted with a spherical model; range
   initialized to one-third of the study-area extent.
2. **IDW fallback** (idp = 2) — activated automatically if kriging fails, or
   forced via the environment variable.

Interpolated values are clipped to [0, 1] and masked to the boundary.
Individual membership rasters are saved as
`outputs/mz_membership_z1.tif` … `outputs/mz_membership_z5.tif`.

**Hard zone map** — each cell is assigned to the zone with the highest
membership (`which.max`), saved as `outputs/mz_zone_map.tif` and
`outputs/mz_zone_map.png`.

**Transition mask** — cells where the maximum membership is below 0.60 are
flagged as transition areas and saved as `outputs/mz_transition.tif`.

---

### Step 7 — Zone statistics and figures

![Statistical outputs](figures/standalone_fig4_statistics.png)

**Zone means** (`outputs/mz_zone_stats.csv`) — per-variable mean of the
original (un-normalized) values grouped by hard zone assignment:

| Zone | BD | CEC | Fragm | Sand | Silt | Clay | N | OCD | pH | SOC |
|------|----|-----|-------|------|------|------|---|-----|----|-----|
| 1 | 1.37 | 13.4 | 15.4 | 78.6 | 12.9 | 8.5 | 3.19 | 14.3 | 7.48 | 13.9 |
| 2 | 1.37 | 13.8 | 16.9 | 78.8 | 12.7 | 8.4 | 3.09 | 14.2 | 7.68 | 14.0 |
| 3 | 1.38 | 14.3 | 15.7 | 77.4 | 13.6 | 9.0 | 2.42 | 13.3 | 7.66 | 12.5 |
| 4 | 1.38 | 13.8 | 17.4 | 78.6 | 12.9 | 8.5 | 3.43 | 14.3 | 7.71 | 14.6 |
| 5 | 1.38 | 14.0 | 17.8 | 79.3 | 12.6 | 8.1 | 4.39 | 15.8 | 7.69 | 13.6 |

**ANOVA** (`outputs/mz_anova.csv`) — one-way ANOVA per soil variable testing
whether zone membership explains significant variation (p < 0.05):

| Variable | p-value | Significant? |
|----------|---------|--------------|
| BD | 0.311 | |
| CEC | 0.049 | ✓ |
| Fragm | < 0.001 | ✓ |
| Sand | 0.007 | ✓ |
| Silt | 0.074 | |
| Clay | 0.072 | |
| N | < 0.001 | ✓ |
| OCD | 0.003 | ✓ |
| pH | < 0.001 | ✓ |
| SOC | 0.221 | |

Six of the ten soil variables show statistically significant differences among
zones, confirming that the five-zone solution captures meaningful spatial
variation in soil properties.

**PCA score plot** — PC1 vs PC2 coloured by zone, saved as
`outputs/mz_pca_scores_by_zone.png`.

---

## Outputs summary

| File | Type | Description |
|------|------|-------------|
| `mz_zone_map.tif` | GeoTIFF INT1U | Hard zone labels (1–k) |
| `mz_transition.tif` | GeoTIFF INT1U | 1 = transition cell (max membership < 0.60) |
| `mz_membership_z*.tif` | GeoTIFF FLT4S | Per-zone fuzzy membership [0, 1] |
| `mz_validation.csv` | CSV | Validity indices for k = 2–6 |
| `mz_zone_stats.csv` | CSV | Mean soil properties per zone |
| `mz_anova.csv` | CSV | ANOVA p-values per soil variable |
| `mz_point_assignments.csv` | CSV | Zone labels and memberships per observation |
| `mz_validation_all_indices.png` | PNG | Scaled validity indices vs k |
| `mz_zone_map.png` | PNG | Rendered hard zone map |
| `mz_pca_scores_by_zone.png` | PNG | PCA biplot coloured by zone |

---

## Configuration

All tunable parameters are declared at the top of the script:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `SOIL_VARS` | `BD CEC … SOC` | Variables extracted and clustered |
| `K_RANGE` | `2:6` | Candidate zone counts evaluated |
| `PCA_THRESHOLD` | `0.80` | Minimum cumulative variance for PC retention |
| `FUZZINESS` | `2` | FCM fuzziness exponent m |
| `TRANSITION_THRESHOLD` | `0.60` | Maximum membership below which a cell is a transition zone |
| `MZ_INTERPOLATION_METHOD` | `kriging` | `"kriging"` or `"idw"` (env var) |

Reproducibility is anchored by `set.seed(42)` applied globally and inside
`fit_cmeans()`.
