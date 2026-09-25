#!/usr/bin/env Rscript

# Build a Texas cotton-area phenology-based ARID grid from NOAA nClimGrid-Daily
# plus local station bias correction, then compute relative yield and relative
# yield loss using the cotton model in Equation 8.
#
# Usage:
#   Rscript scripts/build_arid_texas_cotton_phenology_yield.R \
#     inputs/texas_cotton_station_catalog.csv \
#     Texas \
#     Texas_Cotton \
#     1975 \
#     2025 \
#     outputs/texas_cotton_grid \
#     inputs/texas_cotton_area.geojson \
#     inputs/enso_years_1975_2025.csv

required_packages <- c("terra", "sf")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Please install required R package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

options(tigris_use_cache = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

script_file <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath("scripts", winslash = "/", mustWork = FALSE)
}

source(file.path(script_dir, "compute_arid_batch.R"))

neutral_ratio <- 1
temp_bias_default <- 0
precip_ratio_bounds <- c(0.25, 4)
cotton_phase_base_temp_c <- 15.56
cotton_phase_ttt <- c(60, 230, 225, 200, 370, 200)
cotton_phase_names <- c(
  "planting_emergence",
  "emergence_pinhead_square",
  "pinhead_square_first_bloom",
  "first_bloom_peak_bloom",
  "peak_bloom_first_open_boll",
  "first_open_boll_harvest"
)
cotton_phase_thresholds <- cumsum(cotton_phase_ttt)
yield_model_intercept <- 0.64
yield_phase_exponents <- c(0.01, -0.11, 0.16, 0.09, 0.06, 0.08)

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

normalize_file_value <- function(x) {
  y <- trimws(as.character(x))
  y <- gsub('^"+|"+$', "", y)
  y <- gsub("^'+|'+$", "", y)
  y
}

read_station_weather <- function(path) {
  raw <- read.csv(path, stringsAsFactors = FALSE)
  required <- c("DATE", "PRCP", "TMAX", "TMIN")
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0) {
    stop("Station weather file is missing required columns: ", paste(missing, collapse = ", "), " -> ", path)
  }

  dates <- as.Date(raw$DATE)
  weather <- data.frame(
    date = dates,
    year = as.integer(strftime(dates, "%Y")),
    month = as.integer(strftime(dates, "%m")),
    doy = as.integer(strftime(dates, "%j")),
    tmax = as.numeric(raw$TMAX),
    tmin = as.numeric(raw$TMIN),
    precipitation = as.numeric(raw$PRCP),
    stringsAsFactors = FALSE
  )

  weather <- weather[stats::complete.cases(weather[, c("date", "year", "month", "doy", "tmax", "tmin", "precipitation")]), ]
  weather
}

read_station_catalog <- function(path) {
  stations <- read.csv(path, stringsAsFactors = FALSE)
  required <- c("station_id", "weather_file", "latitude", "longitude", "elevation")
  missing <- setdiff(required, names(stations))
  if (length(missing) > 0) {
    stop("Station catalog is missing required columns: ", paste(missing, collapse = ", "))
  }

  stations$weather_file <- vapply(stations$weather_file, normalize_file_value, character(1))
  stations$latitude <- as.numeric(stations$latitude)
  stations$longitude <- as.numeric(stations$longitude)
  stations$elevation <- as.numeric(stations$elevation)
  stations
}

monthly_dates <- function(start_year, end_year) {
  seq(as.Date(sprintf("%d-01-01", start_year)), as.Date(sprintf("%d-12-01", end_year)), by = "1 month")
}

download_if_missing <- function(url, destination) {
  if (file.exists(destination)) {
    ok_existing <- tryCatch({
      test_raster <- terra::rast(destination)
      terra::nlyr(test_raster) > 0
    }, error = function(e) {
      FALSE
    })
    if (ok_existing) {
      return(destination)
    }
    unlink(destination)
  }

  options(timeout = max(600, getOption("timeout")))
  attempts <- 3
  ok <- FALSE

  for (attempt in seq_len(attempts)) {
    if (file.exists(destination)) {
      unlink(destination)
    }

    result <- tryCatch({
      withCallingHandlers(
        {
          download.file(
            url,
            destination,
            mode = "wb",
            quiet = TRUE,
            method = "libcurl"
          )
        },
        warning = function(w) {
          stop(conditionMessage(w), call. = FALSE)
        }
      )
      TRUE
    }, error = function(e) {
      FALSE
    })

    if (result && file.exists(destination)) {
      valid <- tryCatch({
        test_raster <- terra::rast(destination)
        terra::nlyr(test_raster) > 0
      }, error = function(e) {
        FALSE
      })

      if (valid) {
        ok <- TRUE
        break
      }
    }

    if (file.exists(destination)) {
      unlink(destination)
    }

    Sys.sleep(min(30, 5 * attempt))
  }

  if (!ok) {
    return(NULL)
  }

  destination
}

download_file_simple <- function(url, destination) {
  if (file.exists(destination)) {
    return(destination)
  }

  options(timeout = max(600, getOption("timeout")))
  attempts <- 3

  for (attempt in seq_len(attempts)) {
    if (file.exists(destination)) {
      unlink(destination)
    }

    ok <- tryCatch({
      download.file(url, destination, mode = "wb", quiet = TRUE, method = "libcurl")
      file.exists(destination) && file.info(destination)$size > 0
    }, error = function(e) {
      FALSE
    })

    if (isTRUE(ok)) {
      return(destination)
    }

    Sys.sleep(min(30, 5 * attempt))
  }

  NULL
}

download_nclimgrid_month <- function(year, month, cache_dir) {
  month_key <- sprintf("%04d%02d", year, month)
  month_dir <- file.path(cache_dir, month_key)
  ensure_dir(month_dir)

  nc_files <- list.files(month_dir, pattern = "\\.nc$", recursive = TRUE, full.names = TRUE)
  if (length(nc_files) > 0) {
    return(list(month_dir = month_dir, nc_files = nc_files))
  }

  base_url <- sprintf("https://www.ncei.noaa.gov/thredds/fileServer/nclimgrid-daily/%d", year)
  combined_name <- sprintf("ncdd-%s-grd-scaled.nc", month_key)
  combined_file <- download_if_missing(
    sprintf("%s/%s", base_url, combined_name),
    file.path(month_dir, combined_name)
  )

  if (is.null(combined_file)) {
    variable_names <- c("prcp", "tmin", "tmax")
    for (variable_name in variable_names) {
      variable_file <- sprintf("%s-%s-grd-scaled.nc", variable_name, month_key)
      result <- download_if_missing(
        sprintf("%s/%s", base_url, variable_file),
        file.path(month_dir, variable_file)
      )
      if (is.null(result)) {
        stop(
          "Could not download nClimGrid file for ",
          variable_name,
          " in ",
          month_key,
          "."
        )
      }
    }
  }

  nc_files <- list.files(month_dir, pattern = "\\.nc$", recursive = TRUE, full.names = TRUE)
  if (length(nc_files) == 0) {
    stop("No NetCDF files were downloaded for ", month_key)
  }

  list(month_dir = month_dir, nc_files = nc_files)
}

match_nc_files <- function(nc_files) {
  combined_hits <- nc_files[grepl("^ncdd-.*\\.nc$", basename(nc_files), ignore.case = TRUE)]
  if (length(combined_hits) > 0) {
    return(list(prcp = combined_hits[[1]], tmin = combined_hits[[1]], tmax = combined_hits[[1]]))
  }

  variables <- c(prcp = "prcp", tmin = "tmin", tmax = "tmax")
  result <- lapply(variables, function(name) {
    hits <- nc_files[grepl(name, basename(nc_files), ignore.case = TRUE)]
    if (length(hits) == 0) {
      stop("Could not locate a NetCDF file for variable: ", name)
    }
    hits[[1]]
  })
  result
}

read_nclimgrid_variable <- function(path, variable_name) {
  if (grepl("^ncdd-.*\\.nc$", basename(path), ignore.case = TRUE)) {
    terra::rast(path, subds = variable_name)
  } else {
    terra::rast(path)
  }
}

download_gridmet_year <- function(year, variable_name, cache_dir) {
  if (year < 1979) {
    return(NULL)
  }

  year_dir <- file.path(cache_dir, as.character(year))
  ensure_dir(year_dir)
  destination <- file.path(year_dir, sprintf("%s_%d.nc", variable_name, year))

  direct_urls <- c(
    sprintf("https://www.northwestknowledge.net/metdata/data/%s_%d.nc", variable_name, year),
    sprintf("https://www.northwestknowledge.net/metdata/data/permanent/%s_%d.nc", variable_name, year)
  )

  if (file.exists(destination)) {
    ok_existing <- tryCatch({
      test_raster <- terra::rast(destination)
      terra::nlyr(test_raster) > 0
    }, error = function(e) {
      FALSE
    })
    if (ok_existing) {
      return(destination)
    }
    unlink(destination)
  }

  for (url in direct_urls) {
    result <- download_if_missing(url, destination)
    if (!is.null(result)) {
      return(result)
    }
  }

  stop("Could not download gridMET file for ", variable_name, " in ", year, ".")
}

gridmet_variable_stack <- function(path) {
  terra::rast(path)
}

align_gridmet_to_template <- function(path, template_raster, layer_indices = NULL) {
  source <- terra::rast(path)
  if (!is.null(layer_indices)) {
    source <- source[[layer_indices]]
  }
  source <- terra::crop(source, terra::ext(template_raster), snap = "out")
  terra::resample(source, template_raster, method = "bilinear")
}

is_growing_season <- function(dates) {
  md <- format(dates, "%m-%d")
  md >= "05-01" & md <= "10-15"
}

specific_humidity_to_ea <- function(q, elevation_m) {
  pressure <- 101.3 * ((293 - 0.0065 * elevation_m) / 293)^5.26
  (q * pressure) / (0.622 + 0.378 * q)
}

ea_to_dewpoint <- function(ea) {
  log_term <- log(pmax(ea, 1e-6) / 0.6108)
  (237.3 * log_term) / (17.27 - log_term)
}

state_to_fips <- function(state_name) {
  state_map <- setNames(
    c(
      "01", "02", "04", "05", "06", "08", "09", "10", "11", "12", "13", "15",
      "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27",
      "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39",
      "40", "41", "42", "44", "45", "46", "47", "48", "49", "50", "51", "53",
      "54", "55", "56", "72"
    ),
    c(
      "alabama", "alaska", "arizona", "arkansas", "california", "colorado",
      "connecticut", "delaware", "district of columbia", "florida", "georgia",
      "hawaii", "idaho", "illinois", "indiana", "iowa", "kansas", "kentucky",
      "louisiana", "maine", "maryland", "massachusetts", "michigan", "minnesota",
      "mississippi", "missouri", "montana", "nebraska", "nevada", "new hampshire",
      "new jersey", "new mexico", "new york", "north carolina", "north dakota",
      "ohio", "oklahoma", "oregon", "pennsylvania", "rhode island",
      "south carolina", "south dakota", "tennessee", "texas", "utah", "vermont",
      "virginia", "washington", "west virginia", "wisconsin", "wyoming",
      "puerto rico"
    )
  )
  abbrev_map <- setNames(
    c(
      "01", "02", "04", "05", "06", "08", "09", "10", "11", "12", "13", "15",
      "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27",
      "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39",
      "40", "41", "42", "44", "45", "46", "47", "48", "49", "50", "51", "53",
      "54", "55", "56", "72"
    ),
    c(
      "al", "ak", "az", "ar", "ca", "co", "ct", "de", "dc", "fl", "ga", "hi",
      "id", "il", "in", "ia", "ks", "ky", "la", "me", "md", "ma", "mi", "mn",
      "ms", "mo", "mt", "ne", "nv", "nh", "nj", "nm", "ny", "nc", "nd", "oh",
      "ok", "or", "pa", "ri", "sc", "sd", "tn", "tx", "ut", "vt", "va", "wa",
      "wv", "wi", "wy", "pr"
    )
  )
  key <- tolower(trimws(state_name))
  state_map[[key]] %||% abbrev_map[[key]] %||% NA_character_
}

download_census_counties <- function() {
  zip_url <- "https://www2.census.gov/geo/tiger/GENZ2024/shp/cb_2024_us_county_500k.zip"
  zip_path <- file.path(tempdir(), "cb_2024_us_county_500k.zip")
  unzip_dir <- file.path(tempdir(), "cb_2024_us_county_500k")
  shp_path <- file.path(unzip_dir, "cb_2024_us_county_500k.shp")

  if (!file.exists(shp_path)) {
    ensure_dir(unzip_dir)
    if (!file.exists(zip_path)) {
      result <- download_file_simple(zip_url, zip_path)
      if (is.null(result)) {
        stop("Could not download Census county boundary file from ", zip_url)
      }
    }
    utils::unzip(zip_path, exdir = unzip_dir)
  }

  sf::st_read(shp_path, quiet = TRUE)
}

get_analysis_boundary <- function(state_name, region_name, boundary_file) {
  if (is.null(boundary_file) || boundary_file == "") {
    stop(
      "A cotton-growing-area boundary file is required for the Texas run. ",
      "Pass a GeoJSON, shapefile, or GeoPackage as the seventh argument."
    )
  }

  boundary <- sf::st_read(boundary_file, quiet = TRUE)
  if (!inherits(boundary, "sf")) {
    boundary <- tryCatch(
      sf::st_as_sf(boundary),
      error = function(e) {
        stop("Analysis-area boundary could not be converted to an sf object: ", boundary_file)
      }
    )
  }
  if (is.na(sf::st_crs(boundary))) {
    boundary <- sf::st_set_crs(boundary, 4326)
  }
  boundary <- sf::st_make_valid(sf::st_transform(boundary, 4326))
  if (nrow(boundary) == 0) {
    stop("The analysis-area boundary file contains no features: ", boundary_file)
  }
  boundary
}

extract_station_series <- function(raster, stations_sf, dates, variable_name) {
  extracted <- terra::extract(raster, terra::vect(stations_sf))
  values <- extracted[, -1, drop = FALSE]
  value_matrix <- as.matrix(values)
  data.frame(
    station_id = rep(stations_sf$station_id, times = length(dates)),
    date = rep(dates, each = nrow(stations_sf)),
    variable = variable_name,
    grid_value = as.numeric(value_matrix),
    stringsAsFactors = FALSE
  )
}

build_bias_tables <- function(monthly_data, station_weather, station_info) {
  obs_lookup <- do.call(rbind, lapply(names(station_weather), function(station_id) {
    x <- station_weather[[station_id]]
    data.frame(
      station_id = station_id,
      date = x$date,
      month = x$month,
      year = x$year,
      prcp_obs = x$precipitation,
      tmin_obs = x$tmin,
      tmax_obs = x$tmax,
      stringsAsFactors = FALSE
    )
  }))

  grid_lookup <- do.call(rbind, monthly_data)
  merged <- merge(grid_lookup, obs_lookup, by = c("station_id", "date"), all.x = TRUE, all.y = FALSE)
  merged$month <- as.integer(strftime(merged$date, "%m"))
  merged$year <- as.integer(strftime(merged$date, "%Y"))

  bias_rows <- list()
  active_rows <- list()

  for (station_id in unique(merged$station_id)) {
    station_data <- merged[merged$station_id == station_id, ]
    active_months <- unique(station_data[, c("year", "month")])
    active_months <- active_months[stats::complete.cases(active_months), ]
    if (nrow(active_months) > 0) {
      active_months$station_id <- station_id
      active_rows[[station_id]] <- active_months[, c("station_id", "year", "month")]
    }

    for (month in 1:12) {
      subset_month <- station_data[station_data$month == month, ]

      temp_tmax <- subset_month[stats::complete.cases(subset_month[, c("tmax_obs", "tmax_grid")]), ]
      temp_tmin <- subset_month[stats::complete.cases(subset_month[, c("tmin_obs", "tmin_grid")]), ]
      prcp <- subset_month[stats::complete.cases(subset_month[, c("prcp_obs", "prcp_grid")]), ]

      tmax_bias <- if (nrow(temp_tmax) > 0) mean(temp_tmax$tmax_obs - temp_tmax$tmax_grid, na.rm = TRUE) else temp_bias_default
      tmin_bias <- if (nrow(temp_tmin) > 0) mean(temp_tmin$tmin_obs - temp_tmin$tmin_grid, na.rm = TRUE) else temp_bias_default

      prcp_ratio <- neutral_ratio
      if (nrow(prcp) > 0) {
        denominator <- sum(prcp$prcp_grid, na.rm = TRUE)
        numerator <- sum(prcp$prcp_obs, na.rm = TRUE)
        if (is.finite(denominator) && denominator > 0) {
          prcp_ratio <- numerator / denominator
        }
      }
      prcp_ratio <- max(precip_ratio_bounds[[1]], min(precip_ratio_bounds[[2]], prcp_ratio))

      bias_rows[[paste(station_id, month, sep = "_")]] <- data.frame(
        station_id = station_id,
        month = month,
        tmax_bias = tmax_bias,
        tmin_bias = tmin_bias,
        prcp_ratio = prcp_ratio,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    monthly_bias = do.call(rbind, bias_rows),
    active_months = if (length(active_rows) > 0) do.call(rbind, active_rows) else data.frame()
  )
}

idw_surface <- function(cell_xy, station_xy, station_values, power = 2, default_value = 0) {
  valid <- is.finite(station_values)
  if (!any(valid)) {
    return(rep(default_value, nrow(cell_xy)))
  }

  pts <- station_xy[valid, , drop = FALSE]
  vals <- station_values[valid]

  if (nrow(pts) == 1) {
    return(rep(vals[[1]], nrow(cell_xy)))
  }

  distances <- outer(seq_len(nrow(cell_xy)), seq_len(nrow(pts)), Vectorize(function(i, j) {
    sqrt((cell_xy[i, 1] - pts[j, 1])^2 + (cell_xy[i, 2] - pts[j, 2])^2)
  }))

  exact_match <- distances == 0
  weights <- 1 / pmax(distances, 1e-9)^power
  weights <- weights / rowSums(weights)
  result <- as.numeric(weights %*% vals)

  if (any(exact_match)) {
    match_rows <- which(apply(exact_match, 1, any))
    for (row in match_rows) {
      result[row] <- vals[which(exact_match[row, ])[1]]
    }
  }

  result
}

vectorized_arid_step <- function(tmax, tmin, prcp, doy, year, lat_rad, elevation, state, params, solar_radiation = NULL, ws2 = NULL, ea = NULL) {
  days <- ifelse(year %% 4 == 0, 366, 365)
  tmean <- (tmax + tmin) / 2

  ws2_used <- if (is.null(ws2)) rep(2, length(tmax)) else ws2
  es <- (
    0.6108 * exp(17.27 * tmax / (tmax + 237.3)) +
      0.6108 * exp(17.27 * tmin / (tmin + 237.3))
  ) / 2
  slope <- (0.6108 * exp(17.27 * tmean / (tmean + 237.3)) * 4098) / (tmean + 237.3)^2
  psc <- 0.665 * 10^-3 * 101.3 * ((293 - 0.0065 * elevation) / 293)^5.26

  irdes <- 1 + 0.033 * cos(2 * pi * doy / days)
  sd <- 0.409 * sin(2 * pi * doy / days - 1.39)
  ssa <- acos(-tan(lat_rad) * tan(sd))
  extra <- 24 * 60 * 0.082 / pi * irdes *
    (ssa * sin(lat_rad) * sin(sd) + cos(lat_rad) * cos(sd) * sin(ssa))

  solar_used <- if (is.null(solar_radiation)) {
    estimate_solar_radiation(tmax = tmax, tmin = tmin, extra = extra)
  } else {
    solar_radiation
  }
  swr <- (1 - 0.23) * solar_used
  csr <- (0.75 + 2 * 10^-5 * elevation) * extra
  rrad <- solar_used / csr

  ea_used <- if (is.null(ea)) {
    0.6108 * exp(17.27 * tmin / (tmin + 237.3))
  } else {
    ea
  }
  lwr <- 4.903 * 10^-9 *
    ((tmax + 273.16)^4 + (tmin + 273.16)^4) / 2 *
    (0.34 - 0.14 * sqrt(ea_used)) *
    (1.35 * rrad - 0.35)
  nrad <- swr - lwr
  eto <- (
    0.408 * slope * nrad +
      psc * (900 / (tmean + 273)) * ws2_used * (es - ea_used)
  ) / (slope + psc * (1 + 0.34 * ws2_used))

  s <- 25400 / params$RCN - 254
  runoff <- ifelse(prcp > 0.2 * s, (prcp - 0.2 * s)^2 / (prcp + 0.8 * s), 0)
  cwbd <- prcp - runoff

  wbd <- cwbd + state$wat
  drainage <- ifelse(wbd / params$RZD > params$AWC, params$RZD * params$DDC * (wbd / params$RZD - params$AWC), 0)
  wad <- wbd - drainage
  transpiration <- pmin(params$WUC * wad, eto)
  wat <- wad - transpiration
  arid <- ifelse(eto > 0, 1 - transpiration / eto, NA_real_)

  list(
    eto = eto,
    arid = arid,
    wat = wat,
    precip = prcp,
    tmin = tmin,
    tmax = tmax,
    solar_radiation = solar_used,
    ws2 = ws2_used,
    ea = ea_used,
    tdew = ea_to_dewpoint(ea_used)
  )
}

fetch_oni_table <- function(local_enso_file = NULL) {
  if (!is.null(local_enso_file) && local_enso_file != "") {
    enso <- read.csv(local_enso_file, stringsAsFactors = FALSE)
    required <- c("year", "phase")
    missing <- setdiff(required, names(enso))
    if (length(missing) > 0) {
      stop("Local ENSO file is missing required columns: ", paste(missing, collapse = ", "))
    }
    enso$year <- as.integer(enso$year)
    return(enso[, c("year", "phase")])
  }

  urls <- c(
    "https://www.climate.gov/feeds/dashboard/dashboard-data-oni-graph",
    "https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php"
  )

  monthly <- NULL

  for (url in urls) {
    lines <- tryCatch(readLines(url, warn = FALSE), error = function(e) character(0))
    if (length(lines) == 0) {
      next
    }

    value_lines <- grep("^[0-9]{6},[-]?[0-9]+\\.?[0-9]*$", lines, value = TRUE)
    if (length(value_lines) > 0) {
      parts <- do.call(rbind, strsplit(value_lines, ",", fixed = TRUE))
      monthly <- data.frame(
        year = as.integer(substr(parts[, 1], 1, 4)),
        month = as.integer(substr(parts[, 1], 5, 6)),
        oni = as.numeric(parts[, 2]),
        stringsAsFactors = FALSE
      )
      monthly <- monthly[stats::complete.cases(monthly), ]
      if (nrow(monthly) > 0) {
        break
      }
    }
  }

  if (is.null(monthly) || nrow(monthly) == 0) {
    stop(
      "Could not retrieve ONI data from the configured NOAA sources. ",
      "Please retry with internet access or provide a local ENSO classification table."
    )
  }

  annual <- stats::aggregate(oni ~ year, data = monthly, FUN = function(x) mean(x, na.rm = TRUE))
  annual$phase <- ifelse(
    annual$oni >= 0.5,
    "El Nino",
    ifelse(annual$oni <= -0.5, "La Nina", "Neutral")
  )
  annual
}

annual_raster <- function(template_raster, values, valid_cells) {
  x <- template_raster[[1]]
  out <- x
  filled <- rep(NA_real_, terra::ncell(out))
  filled[valid_cells] <- values
  terra::values(out) <- filled
  out
}

phase_thermal_time <- function(tmax, tmin, base_temp_c = cotton_phase_base_temp_c) {
  tavg <- (tmax + tmin) / 2
  pmax(tavg - base_temp_c, 0)
}

compute_relative_yield <- function(one_minus_arid_by_phase) {
  safe_term <- pmax(one_minus_arid_by_phase, 1e-6)
  yield_model_intercept *
    safe_term[, 1]^yield_phase_exponents[[1]] *
    safe_term[, 2]^yield_phase_exponents[[2]] *
    safe_term[, 3]^yield_phase_exponents[[3]] *
    safe_term[, 4]^yield_phase_exponents[[4]] *
    safe_term[, 5]^yield_phase_exponents[[5]] *
    safe_term[, 6]^yield_phase_exponents[[6]]
}

write_phase_summary <- function(phenology_summary, phase_dir) {
  phase_stats <- stats::aggregate(
    phenology_summary[, c("phase_id", "county_mean_arid", "county_mean_one_minus_arid", "county_total_precip_mm", "county_mean_tmin_c", "county_mean_tmax_c", "county_mean_tdew_c", "county_mean_ws_ms", "county_mean_srad_mj_m2_day", "county_mean_phase_days")],
    by = list(enso_phase = phenology_summary$enso_phase, phenology_phase = phenology_summary$phenology_phase),
    FUN = mean,
    na.rm = TRUE
  )
  phase_stats <- phase_stats[order(phase_stats$enso_phase, phase_stats$phase_id), ]
  write.csv(phase_stats, file.path(phase_dir, "phenology_phase_summary.csv"), row.names = FALSE)
}

run_pipeline <- function(station_catalog_file, state_name, region_name, start_year, end_year, output_dir, boundary_file, local_enso_file) {
  output_dir <- ensure_dir(output_dir)
  cache_dir <- ensure_dir(file.path(output_dir, "cache_nclimgrid"))
  supplemental_cache_dir <- ensure_dir(file.path(output_dir, "cache_gridmet"))
  weather_dir <- ensure_dir(file.path(output_dir, "phenology_weather"))
  arid_dir <- ensure_dir(file.path(output_dir, "phenology_arid"))
  phase_dir <- ensure_dir(file.path(output_dir, "phenology_phase_means"))
  yield_dir <- ensure_dir(file.path(output_dir, "relative_yield_loss"))
  diagnostics_dir <- ensure_dir(file.path(output_dir, "diagnostics"))

  stations <- read_station_catalog(station_catalog_file)
  station_weather <- stats::setNames(lapply(stations$weather_file, read_station_weather), stations$station_id)

  params <- list(
    AWC = param_value(stations[1, ], "AWC"),
    DDC = param_value(stations[1, ], "DDC"),
    RCN = param_value(stations[1, ], "RCN"),
    RZD = param_value(stations[1, ], "RZD"),
    WUC = param_value(stations[1, ], "WUC")
  )

  county_boundary <- get_analysis_boundary(state_name, region_name, boundary_file)
  stations_sf <- sf::st_as_sf(stations, coords = c("longitude", "latitude"), crs = 4326)
  stations_sf$station_id <- stations$station_id

  month_index <- monthly_dates(start_year, end_year)

  message("Pass 1/2: building station-vs-grid bias tables.")
  bias_pass <- list()

  for (month_start in month_index) {
    month_start_date <- as.Date(month_start, origin = "1970-01-01")
    year <- as.integer(format(month_start_date, "%Y"))
    month <- as.integer(format(month_start_date, "%m"))
    if (!(month %in% 5:10)) {
      next
    }
    archive <- download_nclimgrid_month(year, month, cache_dir)
    matched <- match_nc_files(archive$nc_files)

    prcp_raster <- terra::mask(terra::crop(read_nclimgrid_variable(matched$prcp, "prcp"), terra::vect(county_boundary)), terra::vect(county_boundary))
    tmin_raster <- terra::mask(terra::crop(read_nclimgrid_variable(matched$tmin, "tmin"), terra::vect(county_boundary)), terra::vect(county_boundary))
    tmax_raster <- terra::mask(terra::crop(read_nclimgrid_variable(matched$tmax, "tmax"), terra::vect(county_boundary)), terra::vect(county_boundary))
    dates <- seq(month_start_date, by = "1 day", length.out = terra::nlyr(prcp_raster))
    season_keep <- is_growing_season(dates)
    if (!any(season_keep)) {
      next
    }
    prcp_raster <- prcp_raster[[which(season_keep)]]
    tmin_raster <- tmin_raster[[which(season_keep)]]
    tmax_raster <- tmax_raster[[which(season_keep)]]
    dates <- dates[season_keep]

    prcp_series <- extract_station_series(prcp_raster, stations_sf, dates, "prcp")
    tmin_series <- extract_station_series(tmin_raster, stations_sf, dates, "tmin")
    tmax_series <- extract_station_series(tmax_raster, stations_sf, dates, "tmax")

    month_bias <- Reduce(function(x, y) merge(x, y, by = c("station_id", "date"), all = TRUE), list(
      stats::setNames(prcp_series[, c("station_id", "date", "grid_value")], c("station_id", "date", "prcp_grid")),
      stats::setNames(tmin_series[, c("station_id", "date", "grid_value")], c("station_id", "date", "tmin_grid")),
      stats::setNames(tmax_series[, c("station_id", "date", "grid_value")], c("station_id", "date", "tmax_grid"))
    ))
    bias_pass[[format(month_start_date, "%Y-%m")]] <- month_bias
  }

  bias_tables <- build_bias_tables(bias_pass, station_weather, stations)
  write.csv(bias_tables$monthly_bias, file.path(diagnostics_dir, "monthly_station_bias.csv"), row.names = FALSE)
  write.csv(bias_tables$active_months, file.path(diagnostics_dir, "station_active_months.csv"), row.names = FALSE)

  oni_by_year <- fetch_oni_table(local_enso_file)
  write.csv(oni_by_year, file.path(output_dir, "enso_year_classification.csv"), row.names = FALSE)

  message("Pass 2/2: applying bias correction and computing phenology-based ARID grids.")
  template <- NULL
  template_values <- NULL
  county_xy <- NULL
  county_lat_rad <- NULL
  county_elevation <- NULL
  county_cell_count <- NULL
  state <- list()
  phenology_summary <- list()
  phenology_arid_files <- list()
  yield_summary <- list()
  diagnostics_rows <- list()

  current_year <- NULL
  phase_arid_sum <- NULL
  phase_arid_count <- NULL
  phase_prcp_sum <- NULL
  phase_tmin_sum <- NULL
  phase_tmax_sum <- NULL
  phase_tdew_sum <- NULL
  phase_ws_sum <- NULL
  phase_srad_sum <- NULL
  phase_day_count <- NULL

  initialize_phase_arrays <- function(cell_count) {
    list(
      arid_sum = matrix(0, nrow = cell_count, ncol = length(cotton_phase_names)),
      arid_count = matrix(0, nrow = cell_count, ncol = length(cotton_phase_names)),
      prcp_sum = matrix(0, nrow = cell_count, ncol = length(cotton_phase_names)),
      tmin_sum = matrix(0, nrow = cell_count, ncol = length(cotton_phase_names)),
      tmax_sum = matrix(0, nrow = cell_count, ncol = length(cotton_phase_names)),
      tdew_sum = matrix(0, nrow = cell_count, ncol = length(cotton_phase_names)),
      ws_sum = matrix(0, nrow = cell_count, ncol = length(cotton_phase_names)),
      srad_sum = matrix(0, nrow = cell_count, ncol = length(cotton_phase_names)),
      day_count = matrix(0, nrow = cell_count, ncol = length(cotton_phase_names))
    )
  }

  normalize_cell_vector <- function(values, label) {
    values <- as.numeric(values)
    valid_mask <- as.logical(template_values)
    if (length(values) == county_cell_count) {
      return(values)
    }
    if (!is.null(template) && length(values) == terra::ncell(template)) {
      return(values[valid_mask])
    }
    stop(
      "Unexpected ", label, " length: got ", length(values),
      "; expected ", county_cell_count, " valid cells or ",
      if (!is.null(template)) terra::ncell(template) else NA_integer_,
      " template cells. Check that the weather rasters share the template grid."
    )
  }

  finalize_year <- function(year_value) {
    if (is.null(year_value) || is.null(phase_day_count)) {
      return(NULL)
    }

    enso_row <- oni_by_year[oni_by_year$year == year_value, , drop = FALSE]
    enso_value <- if (nrow(enso_row) > 0) enso_row$phase[[1]] else "Unknown"
    phase_arid_means_year <- matrix(NA_real_, nrow = county_cell_count, ncol = length(cotton_phase_names))

    for (phase_id in seq_along(cotton_phase_names)) {
      phase_days <- phase_day_count[, phase_id]
      if (!any(phase_days > 0, na.rm = TRUE)) {
        next
      }

      phase_arid_mean <- phase_arid_sum[, phase_id] / pmax(phase_arid_count[, phase_id], 1)
      phase_tmin_mean <- phase_tmin_sum[, phase_id] / pmax(phase_days, 1)
      phase_tmax_mean <- phase_tmax_sum[, phase_id] / pmax(phase_days, 1)
      phase_tdew_mean <- phase_tdew_sum[, phase_id] / pmax(phase_days, 1)
      phase_ws_mean <- phase_ws_sum[, phase_id] / pmax(phase_days, 1)
      phase_srad_mean <- phase_srad_sum[, phase_id] / pmax(phase_days, 1)
      phase_arid_means_year[, phase_id] <- phase_arid_mean

      phase_name <- cotton_phase_names[[phase_id]]
      file_stub <- sprintf("%d_P%d_%s", year_value, phase_id, phase_name)

      arid_raster <- annual_raster(template, phase_arid_mean, template_values)
      precip_raster <- annual_raster(template, phase_prcp_sum[, phase_id], template_values)
      tmin_raster <- annual_raster(template, phase_tmin_mean, template_values)
      tmax_raster <- annual_raster(template, phase_tmax_mean, template_values)
      tdew_raster <- annual_raster(template, phase_tdew_mean, template_values)
      ws_raster <- annual_raster(template, phase_ws_mean, template_values)
      srad_raster <- annual_raster(template, phase_srad_mean, template_values)

      arid_file <- file.path(arid_dir, sprintf("ARID_%s.tif", file_stub))
      terra::writeRaster(arid_raster, arid_file, overwrite = TRUE)
      terra::writeRaster(precip_raster, file.path(weather_dir, sprintf("precip_mm_%s.tif", file_stub)), overwrite = TRUE)
      terra::writeRaster(tmin_raster, file.path(weather_dir, sprintf("tmin_c_%s.tif", file_stub)), overwrite = TRUE)
      terra::writeRaster(tmax_raster, file.path(weather_dir, sprintf("tmax_c_%s.tif", file_stub)), overwrite = TRUE)
      terra::writeRaster(tdew_raster, file.path(weather_dir, sprintf("tdew_c_%s.tif", file_stub)), overwrite = TRUE)
      terra::writeRaster(ws_raster, file.path(weather_dir, sprintf("wind_ms_%s.tif", file_stub)), overwrite = TRUE)
      terra::writeRaster(srad_raster, file.path(weather_dir, sprintf("srad_MJm2d_%s.tif", file_stub)), overwrite = TRUE)

      phenology_summary[[file_stub]] <<- data.frame(
        year = year_value,
        phase_id = phase_id,
        phenology_phase = phase_name,
        enso_phase = enso_value,
        county_mean_arid = mean(phase_arid_mean, na.rm = TRUE),
        county_mean_one_minus_arid = mean(1 - phase_arid_mean, na.rm = TRUE),
        county_total_precip_mm = mean(phase_prcp_sum[, phase_id], na.rm = TRUE),
        county_mean_tmin_c = mean(phase_tmin_mean, na.rm = TRUE),
        county_mean_tmax_c = mean(phase_tmax_mean, na.rm = TRUE),
        county_mean_tdew_c = mean(phase_tdew_mean, na.rm = TRUE),
        county_mean_ws_ms = mean(phase_ws_mean, na.rm = TRUE),
        county_mean_srad_mj_m2_day = mean(phase_srad_mean, na.rm = TRUE),
        county_mean_phase_days = mean(phase_days, na.rm = TRUE),
        stringsAsFactors = FALSE
      )

      phenology_arid_files[[file_stub]] <<- arid_file
    }

    valid_yield_rows <- rowSums(is.finite(phase_arid_means_year)) == length(cotton_phase_names)
    if (any(valid_yield_rows)) {
      one_minus_arid <- 1 - phase_arid_means_year
      relative_yield <- rep(NA_real_, county_cell_count)
      relative_yield[valid_yield_rows] <- compute_relative_yield(one_minus_arid[valid_yield_rows, , drop = FALSE])
      relative_yield_loss <- 1 - relative_yield

      yield_raster <- annual_raster(template, relative_yield, template_values)
      loss_raster <- annual_raster(template, relative_yield_loss, template_values)

      terra::writeRaster(
        yield_raster,
        file.path(yield_dir, sprintf("relative_yield_%d.tif", year_value)),
        overwrite = TRUE
      )
      terra::writeRaster(
        loss_raster,
        file.path(yield_dir, sprintf("relative_yield_loss_%d.tif", year_value)),
        overwrite = TRUE
      )

      yield_summary[[as.character(year_value)]] <<- data.frame(
        year = year_value,
        enso_phase = enso_value,
        county_mean_relative_yield = mean(relative_yield, na.rm = TRUE),
        county_mean_relative_yield_loss = mean(relative_yield_loss, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  for (month_start in month_index) {
    month_start_date <- as.Date(month_start, origin = "1970-01-01")
    year <- as.integer(format(month_start_date, "%Y"))
    month <- as.integer(format(month_start_date, "%m"))
    if (!(month %in% 5:10)) {
      next
    }

    if (!is.null(current_year) && year != current_year) {
      finalize_year(current_year)
    }

    if (is.null(current_year) || year != current_year) {
      current_year <- year
      if (is.null(template)) {
        archive <- download_nclimgrid_month(year, month, cache_dir)
        matched <- match_nc_files(archive$nc_files)
        template <- terra::mask(terra::crop(read_nclimgrid_variable(matched$prcp, "prcp"), terra::vect(county_boundary)), terra::vect(county_boundary))
        template <- template[[1]]
        valid_cells <- !is.na(terra::values(template))
        template_values <- valid_cells
        xy_all <- terra::xyFromCell(template, seq_len(terra::ncell(template)))
        county_xy <- xy_all[valid_cells, , drop = FALSE]
        county_lat_rad <- county_xy[, 2] * pi / 180
        county_cell_count <- nrow(county_xy)
        county_elevation <- idw_surface(
          cell_xy = county_xy,
          station_xy = sf::st_coordinates(stations_sf),
          station_values = stations$elevation,
          power = 2,
          default_value = stats::median(stations$elevation, na.rm = TRUE)
        )
      }

      phase_arrays <- initialize_phase_arrays(county_cell_count)
      phase_arid_sum <- phase_arrays$arid_sum
      phase_arid_count <- phase_arrays$arid_count
      phase_prcp_sum <- phase_arrays$prcp_sum
      phase_tmin_sum <- phase_arrays$tmin_sum
      phase_tmax_sum <- phase_arrays$tmax_sum
      phase_tdew_sum <- phase_arrays$tdew_sum
      phase_ws_sum <- phase_arrays$ws_sum
      phase_srad_sum <- phase_arrays$srad_sum
      phase_day_count <- phase_arrays$day_count
      state$wat <- rep(params$RZD * params$AWC, county_cell_count)
      state$cum_tt <- rep(0, county_cell_count)
    }

    archive <- download_nclimgrid_month(year, month, cache_dir)
    matched <- match_nc_files(archive$nc_files)

    prcp_raster <- terra::mask(terra::crop(read_nclimgrid_variable(matched$prcp, "prcp"), terra::vect(county_boundary)), terra::vect(county_boundary))
    tmin_raster <- terra::mask(terra::crop(read_nclimgrid_variable(matched$tmin, "tmin"), terra::vect(county_boundary)), terra::vect(county_boundary))
    tmax_raster <- terra::mask(terra::crop(read_nclimgrid_variable(matched$tmax, "tmax"), terra::vect(county_boundary)), terra::vect(county_boundary))
    dates <- seq(month_start_date, by = "1 day", length.out = terra::nlyr(prcp_raster))
    season_keep <- is_growing_season(dates)
    if (!any(season_keep)) {
      next
    }
    prcp_raster <- prcp_raster[[which(season_keep)]]
    tmin_raster <- tmin_raster[[which(season_keep)]]
    tmax_raster <- tmax_raster[[which(season_keep)]]
    dates <- dates[season_keep]

    gridmet_srad_path <- download_gridmet_year(year, "srad", supplemental_cache_dir)
    gridmet_vs_path <- download_gridmet_year(year, "vs", supplemental_cache_dir)
    gridmet_sph_path <- download_gridmet_year(year, "sph", supplemental_cache_dir)

    if (!is.null(gridmet_srad_path) && !is.null(gridmet_vs_path) && !is.null(gridmet_sph_path)) {
      # gridMET and nClimGrid use similar but non-identical 1/24-degree
      # grids. Resample all supplemental fields to the nClimGrid template
      # before applying the valid-cell mask. Only the current month's days
      # are read and resampled to avoid holding three full annual stacks.
      doy_match <- as.integer(format(dates, "%j"))
      gridmet_days <- terra::nlyr(terra::rast(gridmet_srad_path))
      if (gridmet_days == 365L && ((year %% 4L) == 0L)) {
        doy_match <- ifelse(doy_match > 59L, doy_match - 1L, doy_match)
      }
      srad_raster <- align_gridmet_to_template(gridmet_srad_path, template, doy_match)
      vs_raster <- align_gridmet_to_template(gridmet_vs_path, template, doy_match)
      sph_raster <- align_gridmet_to_template(gridmet_sph_path, template, doy_match)
    } else {
      srad_raster <- NULL
      vs_raster <- NULL
      sph_raster <- NULL
    }

    active_month <- bias_tables$active_months[
      bias_tables$active_months$year == year & bias_tables$active_months$month == month,
      ,
      drop = FALSE
    ]

    month_bias <- merge(
      stations[, c("station_id", "longitude", "latitude", "elevation")],
      bias_tables$monthly_bias[bias_tables$monthly_bias$month == month, ],
      by = "station_id",
      all.x = FALSE,
      all.y = TRUE
    )
    if (nrow(active_month) > 0) {
      month_bias <- month_bias[month_bias$station_id %in% active_month$station_id, , drop = FALSE]
    } else {
      month_bias <- month_bias[0, , drop = FALSE]
    }

    station_xy <- if (nrow(month_bias) > 0) as.matrix(month_bias[, c("longitude", "latitude")]) else matrix(numeric(0), ncol = 2)
    tmax_bias_surface <- idw_surface(county_xy, station_xy, month_bias$tmax_bias %||% numeric(0), default_value = temp_bias_default)
    tmin_bias_surface <- idw_surface(county_xy, station_xy, month_bias$tmin_bias %||% numeric(0), default_value = temp_bias_default)
    prcp_ratio_surface <- idw_surface(county_xy, station_xy, month_bias$prcp_ratio %||% numeric(0), default_value = neutral_ratio)

    station_points <- sf::st_coordinates(stations_sf)
    station_tmax_bias <- idw_surface(station_points, station_xy, month_bias$tmax_bias %||% numeric(0), default_value = temp_bias_default)
    station_tmin_bias <- idw_surface(station_points, station_xy, month_bias$tmin_bias %||% numeric(0), default_value = temp_bias_default)
    station_prcp_ratio <- idw_surface(station_points, station_xy, month_bias$prcp_ratio %||% numeric(0), default_value = neutral_ratio)

    raw_prcp_station <- extract_station_series(prcp_raster, stations_sf, dates, "prcp")
    raw_tmin_station <- extract_station_series(tmin_raster, stations_sf, dates, "tmin")
    raw_tmax_station <- extract_station_series(tmax_raster, stations_sf, dates, "tmax")

    corrected_station <- Reduce(function(x, y) merge(x, y, by = c("station_id", "date"), all = TRUE), list(
      stats::setNames(raw_prcp_station[, c("station_id", "date", "grid_value")], c("station_id", "date", "raw_prcp")),
      stats::setNames(raw_tmin_station[, c("station_id", "date", "grid_value")], c("station_id", "date", "raw_tmin")),
      stats::setNames(raw_tmax_station[, c("station_id", "date", "grid_value")], c("station_id", "date", "raw_tmax"))
    ))
    corrected_station <- merge(
      corrected_station,
      do.call(rbind, lapply(names(station_weather), function(station_id) {
        x <- station_weather[[station_id]]
        data.frame(
          station_id = station_id,
          date = x$date,
          obs_prcp = x$precipitation,
          obs_tmin = x$tmin,
          obs_tmax = x$tmax,
          stringsAsFactors = FALSE
        )
      })),
      by = c("station_id", "date"),
      all.x = TRUE
    )
    corrected_station$corrected_prcp <- corrected_station$raw_prcp * station_prcp_ratio[match(corrected_station$station_id, stations$station_id)]
    corrected_station$corrected_tmin <- corrected_station$raw_tmin + station_tmin_bias[match(corrected_station$station_id, stations$station_id)]
    corrected_station$corrected_tmax <- corrected_station$raw_tmax + station_tmax_bias[match(corrected_station$station_id, stations$station_id)]
    diagnostics_rows[[format(month_start_date, "%Y-%m")]] <- corrected_station

    prcp_values <- terra::values(prcp_raster)
    tmin_values <- terra::values(tmin_raster)
    tmax_values <- terra::values(tmax_raster)
    srad_values <- if (!is.null(srad_raster)) terra::values(srad_raster) else NULL
    vs_values <- if (!is.null(vs_raster)) terra::values(vs_raster) else NULL
    sph_values <- if (!is.null(sph_raster)) terra::values(sph_raster) else NULL

    for (layer_index in seq_len(ncol(prcp_values))) {
      prcp_day <- prcp_values[template_values, layer_index] * prcp_ratio_surface
      tmin_day <- tmin_values[template_values, layer_index] + tmin_bias_surface
      tmax_day <- tmax_values[template_values, layer_index] + tmax_bias_surface
      tmax_day <- pmax(tmax_day, tmin_day)

      if (!is.null(srad_values)) {
        solar_mj_day <- pmax(srad_values[template_values, layer_index], 0) * 0.0864
      } else {
        solar_mj_day <- NULL
      }

      if (!is.null(vs_values)) {
        ws_day <- pmax(vs_values[template_values, layer_index], 0)
      } else {
        ws_day <- NULL
      }

      if (!is.null(sph_values)) {
        ea_day <- specific_humidity_to_ea(pmax(sph_values[template_values, layer_index], 0), county_elevation)
      } else {
        ea_day <- NULL
      }

      step <- vectorized_arid_step(
        tmax = tmax_day,
        tmin = tmin_day,
        prcp = pmax(prcp_day, 0),
        doy = as.integer(format(dates[layer_index], "%j")),
        year = year,
        lat_rad = county_lat_rad,
        elevation = county_elevation,
        state = state,
        params = params,
        solar_radiation = solar_mj_day,
        ws2 = ws_day,
        ea = ea_day
      )

      state$wat <- step$wat
      step$arid <- normalize_cell_vector(step$arid, "ARID")
      step$precip <- normalize_cell_vector(step$precip, "precipitation")
      step$tmin <- normalize_cell_vector(step$tmin, "minimum temperature")
      step$tmax <- normalize_cell_vector(step$tmax, "maximum temperature")
      step$tdew <- normalize_cell_vector(step$tdew, "dewpoint")
      step$ws2 <- normalize_cell_vector(step$ws2, "wind speed")
      step$solar_radiation <- normalize_cell_vector(step$solar_radiation, "solar radiation")
      tt_day <- phase_thermal_time(step$tmax, step$tmin)
      active_cells <- state$cum_tt < cotton_phase_thresholds[[length(cotton_phase_thresholds)]]
      phase_ids <- pmin(findInterval(state$cum_tt, vec = cotton_phase_thresholds) + 1L, length(cotton_phase_names))
      valid <- !is.na(step$arid) & active_cells

      for (phase_id in seq_along(cotton_phase_names)) {
        phase_cells <- which(valid & phase_ids == phase_id)
        if (length(phase_cells) == 0L) {
          next
        }
        phase_arid_sum[phase_cells, phase_id] <- phase_arid_sum[phase_cells, phase_id] + step$arid[phase_cells]
        phase_arid_count[phase_cells, phase_id] <- phase_arid_count[phase_cells, phase_id] + 1
        phase_prcp_sum[phase_cells, phase_id] <- phase_prcp_sum[phase_cells, phase_id] + step$precip[phase_cells]
        phase_tmin_sum[phase_cells, phase_id] <- phase_tmin_sum[phase_cells, phase_id] + step$tmin[phase_cells]
        phase_tmax_sum[phase_cells, phase_id] <- phase_tmax_sum[phase_cells, phase_id] + step$tmax[phase_cells]
        phase_tdew_sum[phase_cells, phase_id] <- phase_tdew_sum[phase_cells, phase_id] + step$tdew[phase_cells]
        phase_ws_sum[phase_cells, phase_id] <- phase_ws_sum[phase_cells, phase_id] + step$ws2[phase_cells]
        phase_srad_sum[phase_cells, phase_id] <- phase_srad_sum[phase_cells, phase_id] + step$solar_radiation[phase_cells]
        phase_day_count[phase_cells, phase_id] <- phase_day_count[phase_cells, phase_id] + 1
      }

      state$cum_tt <- state$cum_tt + tt_day
    }
  }

  finalize_year(current_year)

  diagnostics <- do.call(rbind, diagnostics_rows)
  diagnostics_summary <- do.call(rbind, lapply(split(diagnostics, diagnostics$station_id), function(x) {
    data.frame(
      station_id = unique(x$station_id),
      prcp_mae_raw = mean(abs(x$obs_prcp - x$raw_prcp), na.rm = TRUE),
      prcp_mae_corrected = mean(abs(x$obs_prcp - x$corrected_prcp), na.rm = TRUE),
      tmin_mae_raw = mean(abs(x$obs_tmin - x$raw_tmin), na.rm = TRUE),
      tmin_mae_corrected = mean(abs(x$obs_tmin - x$corrected_tmin), na.rm = TRUE),
      tmax_mae_raw = mean(abs(x$obs_tmax - x$raw_tmax), na.rm = TRUE),
      tmax_mae_corrected = mean(abs(x$obs_tmax - x$corrected_tmax), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(diagnostics_summary, file.path(diagnostics_dir, "station_overlap_diagnostics.csv"), row.names = FALSE)

  phenology_summary_df <- do.call(rbind, phenology_summary)
  rownames(phenology_summary_df) <- NULL
  phenology_summary_df <- phenology_summary_df[order(phenology_summary_df$year, phenology_summary_df$phase_id), ]
  write.csv(phenology_summary_df, file.path(output_dir, "phenology_texas_cotton_summary.csv"), row.names = FALSE)
  write_phase_summary(phenology_summary_df, phase_dir)

  if (length(yield_summary) > 0) {
    yield_summary_df <- do.call(rbind, yield_summary)
    rownames(yield_summary_df) <- NULL
    yield_summary_df <- yield_summary_df[order(yield_summary_df$year), ]
    write.csv(yield_summary_df, file.path(output_dir, "relative_yield_summary.csv"), row.names = FALSE)
  }

  for (enso_name in c("El Nino", "La Nina", "Neutral")) {
    for (phase_id in seq_along(cotton_phase_names)) {
      phase_name <- cotton_phase_names[[phase_id]]
      phase_keys <- sprintf(
        "%d_P%d_%s",
        phenology_summary_df$year[phenology_summary_df$enso_phase == enso_name & phenology_summary_df$phase_id == phase_id],
        phase_id,
        phase_name
      )
      phase_keys <- phase_keys[phase_keys %in% names(phenology_arid_files)]
      if (length(phase_keys) == 0) {
        next
      }
      rasters <- terra::rast(unname(unlist(phenology_arid_files[phase_keys])))
      phase_mean <- terra::app(rasters, fun = mean, na.rm = TRUE)
      phase_file <- file.path(
        phase_dir,
        sprintf("%s_P%d_%s_mean_ARID.tif", gsub(" ", "_", enso_name), phase_id, phase_name)
      )
      terra::writeRaster(phase_mean, phase_file, overwrite = TRUE)
    }
  }

  invisible(list(
    phenology_summary = phenology_summary_df,
    diagnostics = diagnostics_summary
  ))
}

args <- commandArgs(trailingOnly = TRUE)
  if (sys.nframe() == 0) {
  if (length(args) != 8) {
    stop(
      "Usage: Rscript scripts/build_arid_texas_cotton_phenology_yield.R ",
      "inputs/texas_cotton_station_catalog.csv Texas Texas_Cotton 1975 2025 outputs/texas_cotton_grid cotton_area_boundary.geojson enso.csv"
    )
  }

  run_pipeline(
    station_catalog_file = args[[1]],
    state_name = args[[2]],
    region_name = args[[3]],
    start_year = as.integer(args[[4]]),
    end_year = as.integer(args[[5]]),
    output_dir = args[[6]],
    boundary_file = args[[7]],
    local_enso_file = args[[8]]
  )
}
