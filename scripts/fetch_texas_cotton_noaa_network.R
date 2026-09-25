#!/usr/bin/env Rscript

# Discover and download a thinned NOAA GHCN-Daily station network for the
# Texas cotton area. Daily observations are downloaded one station-year at a
# time because NOAA CDO limits large daily-data date ranges.
#
# Usage:
#   Rscript scripts/fetch_texas_cotton_noaa_network.R \
#     inputs/texas_cotton_area.geojson \
#     1975 2025 \
#     outputs/texas_cotton_noaa \
#     100
#
# Required environment variable:
#   NOAA_CDO_TOKEN

base_url <- "https://www.ncei.noaa.gov/cdo-web/api/v2"
required_packages <- c("jsonlite", "sf")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Please install required R package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

token <- Sys.getenv("NOAA_CDO_TOKEN")
if (token == "") stop("Set NOAA_CDO_TOKEN in your environment before running this script.", call. = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x[[1]]) || identical(x[[1]], "")) y else x
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

as_query_string <- function(query) {
  query <- Filter(function(x) !is.null(x) && length(x) > 0 && !is.na(x[[1]]), query)
  if (length(query) == 0) return("")
  paste(
    names(query),
    vapply(query, function(x) utils::URLencode(as.character(x[[1]]), reserved = TRUE), character(1)),
    sep = "=", collapse = "&"
  )
}

cdo_get <- function(endpoint, query = list(), attempts = 4L) {
  query_string <- as_query_string(query)
  request_url <- paste0(base_url, endpoint, if (query_string != "") paste0("?", query_string) else "")
  last_error <- NULL

  for (attempt in seq_len(attempts)) {
    result <- tryCatch({
      con <- url(request_url, open = "rb", headers = c(token = token))
      on.exit(close(con), add = TRUE)
      raw <- readBin(con, what = "raw", n = 100000000)
      if (length(raw) == 0) stop("NOAA returned an empty response.")
      jsonlite::fromJSON(rawToChar(raw), simplifyDataFrame = TRUE)
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      NULL
    })

    if (!is.null(result)) {
      Sys.sleep(0.25)
      return(result)
    }
    Sys.sleep(min(30, 2 ^ attempt))
  }

  stop("NOAA CDO request failed after ", attempts, " attempts: ", last_error, call. = FALSE)
}

bind_results <- function(chunks) {
  if (length(chunks) == 0) return(data.frame(stringsAsFactors = FALSE))
  columns <- unique(unlist(lapply(chunks, names)))
  chunks <- lapply(chunks, function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    missing <- setdiff(columns, names(x))
    for (column in missing) x[[column]] <- NA
    x[, columns, drop = FALSE]
  })
  do.call(rbind, chunks)
}

cdo_get_all <- function(endpoint, query = list(), page_limit = 1000L) {
  offset <- 1L
  chunks <- list()

  repeat {
    payload <- cdo_get(endpoint, c(query, list(limit = page_limit, offset = offset)))
    batch <- payload$results
    if (is.null(batch) || length(batch) == 0) break
    batch <- as.data.frame(batch, stringsAsFactors = FALSE)
    chunks[[length(chunks) + 1L]] <- batch

    count <- suppressWarnings(as.integer(payload$metadata$resultset$count %||% NA_integer_))
    if (nrow(batch) < page_limit || is.na(count) || offset + nrow(batch) - 1L >= count) break
    offset <- offset + nrow(batch)
  }

  bind_results(chunks)
}

read_area <- function(area_file) {
  if (!file.exists(area_file)) stop("Cotton-area boundary was not found: ", area_file, call. = FALSE)
  area <- sf::st_read(area_file, quiet = TRUE)
  if (is.na(sf::st_crs(area))) sf::st_crs(area) <- 4326
  area <- sf::st_make_valid(sf::st_transform(area, 4326))
  area <- area[!sf::st_is_empty(area), , drop = FALSE]
  if (nrow(area) == 0) stop("Cotton-area boundary contains no geometries.", call. = FALSE)
  area
}

discover_stations <- function(area, start_date, end_date, max_stations = 100L) {
  bbox <- sf::st_bbox(area)
  extent <- paste(bbox[["ymin"]], bbox[["xmin"]], bbox[["ymax"]], bbox[["xmax"]], sep = ",")

  stations <- cdo_get_all(
    "/stations",
    query = list(
      datasetid = "GHCND",
      startdate = start_date,
      enddate = end_date,
      extent = extent,
      datatypeid = "PRCP,TMAX,TMIN",
      sortfield = "datacoverage",
      sortorder = "desc"
    )
  )
  if (nrow(stations) == 0) stop("NOAA returned no GHCND stations in the cotton-area bounding box.", call. = FALSE)

  required <- c("id", "name", "latitude", "longitude")
  missing <- setdiff(required, names(stations))
  if (length(missing) > 0) stop("NOAA station response is missing: ", paste(missing, collapse = ", "), call. = FALSE)

  stations$latitude <- as.numeric(stations$latitude)
  stations$longitude <- as.numeric(stations$longitude)
  stations <- stations[is.finite(stations$latitude) & is.finite(stations$longitude), , drop = FALSE]
  stations_sf <- sf::st_as_sf(stations, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  inside <- lengths(sf::st_intersects(stations_sf, area)) > 0L
  stations <- stations[inside, , drop = FALSE]
  if (nrow(stations) == 0) stop("No NOAA GHCND stations fall inside the cotton-area boundary.", call. = FALSE)

  # Retain at most one strong station per approximately 0.5-degree cell before
  # applying the final cap, which avoids over-representing dense urban clusters.
  stations$network_cell <- paste(floor(stations$latitude / 0.5), floor(stations$longitude / 0.5), sep = "_")
  coverage <- if ("datacoverage" %in% names(stations)) as.numeric(stations$datacoverage) else rep(0, nrow(stations))
  coverage[!is.finite(coverage)] <- -Inf
  stations$coverage_rank <- coverage
  stations <- stations[order(stations$network_cell, -stations$coverage_rank), , drop = FALSE]
  stations <- stations[!duplicated(stations$network_cell), , drop = FALSE]
  stations <- stations[order(-stations$coverage_rank, stations$id), , drop = FALSE]
  stations <- head(stations, max_stations)
  stations$location_id <- make.unique(gsub("[^A-Za-z0-9]+", "_", stations$id))
  stations
}

fetch_station_year <- function(station_id, year) {
  daily <- cdo_get_all(
    "/data",
    query = list(
      datasetid = "GHCND",
      stationid = station_id,
      startdate = sprintf("%04d-01-01", year),
      enddate = sprintf("%04d-12-31", year),
      units = "metric",
      includemetadata = "false",
      datatypeid = "PRCP,TMAX,TMIN,AWND,ADPT"
    )
  )
  if (nrow(daily) == 0) return(NULL)
  required <- c("date", "datatype", "value")
  if (length(setdiff(required, names(daily))) > 0) return(NULL)

  daily$date <- as.Date(daily$date)
  daily$value <- as.numeric(daily$value)
  wide <- reshape(daily[, c("date", "datatype", "value")], idvar = "date", timevar = "datatype", direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))
  names(wide)[names(wide) == "date"] <- "DATE"
  for (field in c("PRCP", "TMAX", "TMIN", "ADPT", "AWND")) {
    if (!field %in% names(wide)) wide[[field]] <- NA_real_
  }
  wide <- wide[order(wide$DATE), c("DATE", "PRCP", "TMAX", "TMIN", "ADPT", "AWND"), drop = FALSE]
  wide
}

fetch_station <- function(station, start_year, end_year, weather_dir) {
  chunks <- list()
  for (year in seq.int(start_year, end_year)) {
    message("  ", station$id, " ", year)
    chunk <- tryCatch(
      fetch_station_year(station$id, year),
      error = function(e) {
        warning("Skipping ", station$id, " ", year, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(chunk)) chunks[[length(chunks) + 1L]] <- chunk
  }
  if (length(chunks) == 0) return(NULL)

  data <- do.call(rbind, chunks)
  data <- data[!duplicated(data$DATE), , drop = FALSE]
  data$STATION <- station$id
  data$NAME <- station$name %||% NA_character_
  data$LATITUDE <- as.numeric(station$latitude)
  data$LONGITUDE <- as.numeric(station$longitude)
  data$ELEVATION <- as.numeric(station$elevation %||% NA_real_)
  data <- data[, c("STATION", "NAME", "LATITUDE", "LONGITUDE", "ELEVATION", "DATE", "ADPT", "AWND", "PRCP", "TMAX", "TMIN"), drop = FALSE]

  output_file <- file.path(weather_dir, paste0(station$location_id, "_NOAA_daily.csv"))
  write.csv(data, output_file, row.names = FALSE, na = "")
  output_file
}

run_download <- function(area_file, start_year, end_year, output_dir, max_stations = 100L) {
  start_year <- as.integer(start_year)
  end_year <- as.integer(end_year)
  if (!is.finite(start_year) || !is.finite(end_year) || start_year > end_year) stop("Invalid year range.", call. = FALSE)
  output_dir <- ensure_dir(output_dir)
  weather_dir <- ensure_dir(file.path(output_dir, "weather"))
  area <- read_area(area_file)
  stations <- discover_stations(area, sprintf("%04d-01-01", start_year), sprintf("%04d-12-31", end_year), as.integer(max_stations))

  station_file <- file.path(output_dir, "texas_cotton_noaa_stations.csv")
  write.csv(stations, station_file, row.names = FALSE, na = "")

  catalog <- list()
  for (i in seq_len(nrow(stations))) {
    station <- stations[i, , drop = FALSE]
    message("Downloading station ", i, "/", nrow(stations), ": ", station$id)
    weather_file <- fetch_station(station, start_year, end_year, weather_dir)
    if (is.null(weather_file)) {
      warning("No data downloaded for station ", station$id)
      next
    }
    catalog[[length(catalog) + 1L]] <- data.frame(
      station_id = station$id,
      weather_file = normalizePath(weather_file, winslash = "/", mustWork = FALSE),
      latitude = as.numeric(station$latitude),
      longitude = as.numeric(station$longitude),
      elevation = as.numeric(station$elevation %||% NA_real_),
      AWC = 0.09,
      DDC = 0.60,
      RCN = 58,
      RZD = 500,
      WUC = 0.096,
      stringsAsFactors = FALSE
    )
  }

  if (length(catalog) == 0) stop("No station weather files were successfully downloaded.", call. = FALSE)
  catalog <- do.call(rbind, catalog)
  catalog_file <- file.path(output_dir, "texas_cotton_station_catalog.csv")
  write.csv(catalog, catalog_file, row.names = FALSE, na = "")
  message("Wrote station catalog: ", catalog_file)
  invisible(catalog)
}

args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0) {
  if (length(args) < 4L || length(args) > 5L) {
    stop("Usage: Rscript scripts/fetch_texas_cotton_noaa_network.R AREA_GEOJSON START_YEAR END_YEAR OUTPUT_DIR [MAX_STATIONS]", call. = FALSE)
  }
  run_download(args[[1]], args[[2]], args[[3]], args[[4]], if (length(args) == 5L) args[[5]] else 100L)
}
