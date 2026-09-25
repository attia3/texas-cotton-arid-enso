# Texas Cotton ARID–ENSO Analysis

R workflow supporting the manuscript **“A Phenology-Based Agricultural Reference Index for Mapping Texas Cotton Drought Stress and Relative Yield.”** The analysis combines a Cropland Data Layer (CDL) cotton footprint, NOAA station and gridded weather data, a phenology-aware Agricultural Reference Index for Drought (ARID) model, a relative-yield response equation, and ENSO classifications for 1975–2025.

The repository includes the derived cotton-area boundary and the exact annual ENSO table used by the workflow. Large weather downloads, caches, and derived rasters are regenerated locally and are not stored in Git.

## Workflow

1. Convert CDL class 2 (cotton) to an approximately 5-km cotton-area footprint.
2. Discover and download a spatially thinned NOAA GHCN-Daily station network.
3. Download nClimGrid-Daily and gridMET inputs, bias-correct precipitation and temperature, calculate phase-specific ARID, and estimate annual relative yield and yield loss.
4. Produce manuscript figures and summary tables, including ENSO composites.

## Repository structure

```text
scripts/
  prepare_texas_cotton_area.R
  fetch_texas_cotton_noaa_network.R
  compute_arid_batch.R
  build_arid_texas_cotton_phenology_yield.R
  make_manuscript_figures.R
inputs/
  README.md
  enso_years_1975_2025.csv
  texas_cotton_area.geojson
examples/
  locations_template.csv
outputs/
  README.md
```

Large downloaded or derived rasters, NetCDF files, station observations, and model outputs are intentionally excluded from Git. See [`inputs/README.md`](inputs/README.md) for the required inputs and provenance fields.

## Requirements

- R 4.1 or newer
- R packages: `terra`, `sf`, `jsonlite`, `ggplot2`, `dplyr`, and `viridis`
- System geospatial libraries required by `sf` and `terra` (commonly GDAL, GEOS, and PROJ)
- A free [NOAA Climate Data Online token](https://www.ncdc.noaa.gov/cdo-web/token)

Install the R packages with:

```r
install.packages(c("terra", "sf", "jsonlite", "ggplot2", "dplyr", "viridis"))
```

## Run the analysis

Run commands from the repository root.

### 1. Prepare or reuse the cotton footprint

The derived boundary used by the workflow is included as `inputs/texas_cotton_area.geojson`. To recreate a boundary from a source CDL raster, run:

```bash
Rscript scripts/prepare_texas_cotton_area.R \
  path/to/cdl_texas.tif \
  inputs/texas_cotton_area.geojson \
  167
```

The aggregation factor of 167 converts 30-m CDL cells to an approximately 5-km occupancy grid. The supplied source raster was exported on February 28, 2025 and most likely represents the 2024 CDL, but the source-year field was not retained in its metadata. The included GeoJSON is therefore the authoritative analysis mask for reproducing the reported run.

### 2. Download the NOAA station network

Store the NOAA token in an environment variable, never in a script or committed file.

macOS/Linux:

```bash
export NOAA_CDO_TOKEN="your-token"
```

Windows PowerShell:

```powershell
$env:NOAA_CDO_TOKEN = "your-token"
```

Then run:

```bash
Rscript scripts/fetch_texas_cotton_noaa_network.R \
  inputs/texas_cotton_area.geojson \
  1975 2025 \
  outputs/texas_cotton_noaa \
  100
```

### 3. Run the gridded ARID and relative-yield workflow

The main script sources the included ARID functions from `scripts/compute_arid_batch.R` and uses the included annual ENSO table.

```bash
Rscript scripts/build_arid_texas_cotton_phenology_yield.R \
  outputs/texas_cotton_noaa/texas_cotton_station_catalog.csv \
  Texas \
  Texas_Cotton \
  1975 \
  2025 \
  outputs/texas_cotton_grid \
  inputs/texas_cotton_area.geojson \
  inputs/enso_years_1975_2025.csv
```

This stage downloads substantial gridded data and can require considerable disk space and runtime.

The standalone location-level ARID calculator can also be run with the cleaned template in `examples/locations_template.csv`:

```bash
Rscript scripts/compute_arid_batch.R \
  examples/locations_template.csv \
  outputs/arid_locations
```

### 4. Create manuscript figures and tables

```bash
Rscript scripts/make_manuscript_figures.R \
  outputs/texas_cotton_grid \
  inputs/texas_cotton_area.geojson \
  outputs/texas_cotton_noaa/texas_cotton_station_catalog.csv
```

## Primary data sources

- [USDA NASS CropScape/Cropland Data Layer](https://nassgeodata.gmu.edu/CropScape/)
- [NOAA Climate Data Online](https://www.ncei.noaa.gov/cdo-web/)
- [NOAA nClimGrid-Daily](https://www.ncei.noaa.gov/products/land-based-station/nclimgrid-daily)
- [gridMET](https://www.climatologylab.org/gridmet.html)
- [U.S. Census Bureau cartographic boundary files](https://www.census.gov/geographies/mapping-files/time-series/geo/cartographic-boundary.html)
- [NOAA Climate Prediction Center Oceanic Niño Index](https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php)

Users are responsible for following the source providers’ terms, documenting the precise data versions used, and citing each data product appropriately.

## Scientific scope

The annual outputs are modeled **relative** yield responses, not direct predictions of lint yield. Important assumptions include spatially uniform soil/crop parameters, transferred yield-response coefficients, station-network availability, fallback meteorological estimates for early years, and an annualized ENSO classification. These limitations should accompany any interpretation of the maps or summaries.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). 
## License

The code is released under the [MIT License](LICENSE). Input datasets remain subject to their respective providers’ terms and citation requirements.
