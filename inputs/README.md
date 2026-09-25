# Input data

This directory contains the small, versioned inputs required to reproduce the analysis. Large binaries, downloaded observations, and provider-controlled datasets remain ignored by Git.

## Required inputs

### Cropland Data Layer raster

Supply a single-band categorical CDL raster covering Texas. The preprocessing script identifies class `2` as cotton:

```text
path/to/cdl_texas.tif
```

Record the following before running the analysis:

- CDL year
- download date
- CropScape export or download URL
- geographic extent and projection
- any clipping, reprojection, or category filtering applied before export

The source raster was exported on February 28, 2025. It most likely represents the 2024 CDL—the latest annual CDL expected to have been available at that time—but the source-year field is absent from the retained metadata. Do not silently present the year as verified. Use the included `texas_cotton_area.geojson` as the authoritative fixed boundary for reproducing the reported model run.

### ENSO classification

The main pipeline expects:

```text
inputs/enso_years_1975_2025.csv
```

Required columns:

| Column | Description |
|---|---|
| `year` | Four-digit year |
| `phase` | `El Nino`, `La Nina`, or `Neutral` |

The exact 51-row table used by the workflow is included. The manuscript uses annual mean ONI thresholds of ≥0.5 for El Niño, ≤−0.5 for La Niña, and Neutral otherwise; this is a study-specific annual aggregation rather than NOAA’s official sustained three-month event definition.

## Downloaded during execution

The scripts retrieve NOAA GHCN-Daily observations, nClimGrid-Daily fields, gridMET fields, and Census cartographic boundaries. Keep their caches under `outputs/` or another ignored local directory.
