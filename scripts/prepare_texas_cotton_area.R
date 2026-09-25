#!/usr/bin/env Rscript

# Convert a categorical CDL raster into a coarse cotton-area polygon.
# CDL class 2 is cotton. The coarse mask prevents vectorizing millions of
# individual 30 m pixels before the statewide grid workflow begins.
#
# Usage:
#   Rscript scripts/prepare_texas_cotton_area.R \
#     "C:/.../clipped.TIF" inputs/texas_cotton_area.geojson 167

required_packages <- c("terra", "sf")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Please install required R package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

ensure_parent <- function(path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
}

prepare_cotton_area <- function(cdl_file, output_file, aggregate_factor = 167L) {
  if (!file.exists(cdl_file)) {
    stop("CDL raster was not found: ", cdl_file, call. = FALSE)
  }

  aggregate_factor <- as.integer(aggregate_factor)
  if (is.na(aggregate_factor) || aggregate_factor < 1L) {
    stop("aggregate_factor must be a positive integer.", call. = FALSE)
  }

  cdl <- terra::rast(cdl_file)
  if (terra::nlyr(cdl) != 1L) stop("Expected a single-band categorical CDL raster.", call. = FALSE)

  cdl_levels <- terra::levels(cdl)[[1]]
  if (!is.null(cdl_levels) && "Value" %in% names(cdl_levels)) {
    if (!any(as.integer(cdl_levels$Value) == 2L)) {
      stop("The raster categories do not contain CDL class 2 (Cotton).", call. = FALSE)
    }
  }

  message("Reading CDL raster: ", cdl_file)
  message("Raster resolution: ", paste(terra::res(cdl), collapse = " x "), " map units")
  message("Aggregating by factor ", aggregate_factor, " (approximately 5 km at 30 m source resolution).")

  cotton <- terra::ifel(cdl == 2, 1, 0)
  coarse <- terra::aggregate(cotton, fact = aggregate_factor, fun = max, na.rm = TRUE)
  coarse[is.infinite(coarse)] <- NA
  cotton_mask <- terra::ifel(coarse > 0, 1, NA)

  message("Converting coarse cotton mask to polygons.")
  polygons <- terra::as.polygons(cotton_mask, dissolve = TRUE, na.rm = TRUE)
  polygons <- sf::st_as_sf(polygons)
  names(polygons)[names(polygons) == names(polygons)[1]] <- "cotton"

  if (is.na(sf::st_crs(polygons))) {
    sf::st_crs(polygons) <- terra::crs(cdl, proj = TRUE)
  }
  polygons <- sf::st_make_valid(sf::st_transform(polygons, 4326))
  polygons <- polygons[!sf::st_is_empty(polygons), , drop = FALSE]
  polygons <- polygons[polygons$cotton > 0, , drop = FALSE]

  if (nrow(polygons) == 0) stop("No cotton pixels (class 2) were found in the CDL raster.", call. = FALSE)

  ensure_parent(output_file)
  if (file.exists(output_file)) unlink(output_file)
  sf::st_write(polygons, output_file, quiet = TRUE, append = FALSE)

  message("Wrote cotton-area boundary: ", normalizePath(output_file, winslash = "/", mustWork = FALSE))
  invisible(polygons)
}

args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0) {
  if (length(args) < 2L || length(args) > 3L) {
    stop("Usage: Rscript scripts/prepare_texas_cotton_area.R CDL_TIF OUTPUT_GEOJSON [AGGREGATE_FACTOR]", call. = FALSE)
  }
  prepare_cotton_area(
    cdl_file = args[[1]],
    output_file = args[[2]],
    aggregate_factor = if (length(args) == 3L) args[[3]] else 167L
  )
}
