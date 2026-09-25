#!/usr/bin/env Rscript

# Batch ARID calculator for multiple Texas cotton locations.
#
# Usage:
#   Rscript scripts/compute_arid_batch.R inputs/locations.csv outputs/arid
#
# locations.csv columns:
#   location_id,weather_file,latitude,elevation
#
# Optional format column:
#   weather_format ("csv"/"noaa_daily" for NOAA daily files)
#
# Optional location-specific parameter columns:
#   AWC,DDC,RCN,RZD,WUC
#
# Weather file columns, with or without a header:
#   solar_radiation,tmax,tmin,tdew,precipitation,wind_speed_10m,doy,year

default_params <- list(
  AWC = 0.09,
  DDC = 0.60,
  RCN = 58,
  RZD = 500,
  WUC = 0.096
)

weather_names <- c(
  "solar_radiation",
  "tmax",
  "tmin",
  "tdew",
  "precipitation",
  "wind_speed_10m",
  "doy",
  "year"
)

to_doy <- function(date_value) {
  as.integer(strftime(as.Date(date_value), format = "%j"))
}

to_year <- function(date_value) {
  as.integer(strftime(as.Date(date_value), format = "%Y"))
}

estimate_solar_radiation <- function(tmax, tmin, extra, coastal = FALSE) {
  # FAO-56 temperature-range estimate for Rs when observations are unavailable.
  krs <- if (coastal) 0.19 else 0.16
  krs * sqrt(pmax(tmax - tmin, 0)) * extra
}

read_weather <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE)
  has_header <- grepl("[A-Za-z]", first_line)
  weather <- read.table(
    path,
    header = has_header,
    sep = "",
    stringsAsFactors = FALSE
  )

  if (ncol(weather) < length(weather_names)) {
    stop("Weather file must contain at least 8 columns: ", path)
  }

  weather <- weather[, seq_along(weather_names)]
  names(weather) <- weather_names
  weather
}

read_noaa_daily <- function(path) {
  raw <- read.csv(path, stringsAsFactors = FALSE)
  required <- c("DATE", "PRCP", "TMAX", "TMIN")
  missing <- setdiff(required, names(raw))

  if (length(missing) > 0) {
    stop("NOAA daily file is missing required columns: ", paste(missing, collapse = ", "))
  }

  data.frame(
    solar_radiation = NA_real_,
    tmax = as.numeric(raw$TMAX),
    tmin = as.numeric(raw$TMIN),
    tdew = if ("ADPT" %in% names(raw)) as.numeric(raw$ADPT) else NA_real_,
    precipitation = as.numeric(raw$PRCP),
    wind_speed_10m = if ("AWND" %in% names(raw)) as.numeric(raw$AWND) else NA_real_,
    doy = to_doy(raw$DATE),
    year = to_year(raw$DATE),
    stringsAsFactors = FALSE
  )
}

load_weather <- function(path, format = NULL) {
  selected_format <- if (is.null(format) || is.na(format) || format == "") {
    tolower(tools::file_ext(path))
  } else {
    tolower(format)
  }

  if (selected_format %in% c("csv", "noaa_daily")) {
    read_noaa_daily(path)
  } else {
    read_weather(path)
  }
}

param_value <- function(location, name) {
  value <- location[[name]]
  if (is.null(value) || length(value) == 0 || is.na(value) || value == "") {
    return(default_params[[name]])
  }
  as.numeric(value)
}

compute_arid <- function(weather, latitude, elevation, params = default_params) {
  n <- nrow(weather)
  eto <- numeric(n)
  wbd <- numeric(n)
  wat <- numeric(n)
  arid <- numeric(n)

  lat <- latitude * pi / 180
  psc <- 0.665 * 10^-3 * 101.3 * ((293 - 0.0065 * elevation) / 293)^5.26

  for (i in seq_len(n)) {
    days <- if (weather$year[i] %% 4 == 0) 366 else 365
    tmean <- (weather$tmax[i] + weather$tmin[i]) / 2

    ws2 <- weather$wind_speed_10m[i] * 4.87 / log(67.8 * 10 - 5.42)
    es <- (
      0.6108 * exp(17.27 * weather$tmax[i] / (weather$tmax[i] + 237.3)) +
        0.6108 * exp(17.27 * weather$tmin[i] / (weather$tmin[i] + 237.3))
    ) / 2
    slope <- (
      0.6108 * exp(17.27 * tmean / (tmean + 237.3)) * 4098
    ) / (tmean + 237.3)^2

    irdes <- 1 + 0.033 * cos(2 * pi * weather$doy[i] / days)
    sd <- 0.409 * sin(2 * pi * weather$doy[i] / days - 1.39)
    ssa <- acos(-tan(lat) * tan(sd))
    extra <- 24 * 60 * 0.082 / pi * irdes *
      (ssa * sin(lat) * sin(sd) + cos(lat) * cos(sd) * sin(ssa))
    if (is.na(weather$solar_radiation[i])) {
      weather$solar_radiation[i] <- estimate_solar_radiation(
        tmax = weather$tmax[i],
        tmin = weather$tmin[i],
        extra = extra
      )
    }
    swr <- (1 - 0.23) * weather$solar_radiation[i]
    csr <- (0.75 + 2 * 10^-5 * elevation) * extra
    rrad <- weather$solar_radiation[i] / csr

    ea <- if (!is.na(weather$tdew[i])) {
      0.6108 * exp(17.27 * weather$tdew[i] / (weather$tdew[i] + 237.3))
    } else {
      0.6108 * exp(17.27 * weather$tmin[i] / (weather$tmin[i] + 237.3))
    }
    lwr <- 4.903 * 10^-9 *
      ((weather$tmax[i] + 273.16)^4 + (weather$tmin[i] + 273.16)^4) / 2 *
      (0.34 - 0.14 * sqrt(ea)) *
      (1.35 * rrad - 0.35)
    nrad <- swr - lwr

    ws2_value <- if (!is.na(ws2)) ws2 else 2
    eto[i] <- (
      0.408 * slope * nrad +
        psc * (900 / (tmean + 273)) * ws2_value * (es - ea)
    ) / (slope + psc * (1 + 0.34 * ws2_value))

    s <- 25400 / params$RCN - 254
    runoff <- if (weather$precipitation[i] > 0.2 * s) {
      (weather$precipitation[i] - 0.2 * s)^2 / (weather$precipitation[i] + 0.8 * s)
    } else {
      0
    }

    cwbd <- weather$precipitation[i] - runoff
    water_after_transpiration <- if (i == 1) {
      params$RZD * params$AWC
    } else {
      wat[i - 1]
    }

    wbd[i] <- cwbd + water_after_transpiration
    drainage <- if (wbd[i] / params$RZD > params$AWC) {
      params$RZD * params$DDC * (wbd[i] / params$RZD - params$AWC)
    } else {
      0
    }

    water_after_drainage <- wbd[i] - drainage
    transpiration <- min(params$WUC * water_after_drainage, eto[i])
    wat[i] <- water_after_drainage - transpiration
    arid[i] <- if (eto[i] > 0) 1 - transpiration / eto[i] else NA_real_
  }

  data.frame(
    year = weather$year,
    doy = weather$doy,
    precipitation = weather$precipitation,
    solar_radiation = weather$solar_radiation,
    ETo = eto,
    WBD = wbd,
    WAT = wat,
    ARID = arid
  )
}

run_batch <- function(locations_file, output_dir) {
  locations <- read.csv(locations_file, stringsAsFactors = FALSE)
  required <- c("location_id", "weather_file", "latitude", "elevation")
  missing <- setdiff(required, names(locations))

  if (length(missing) > 0) {
    stop("locations.csv is missing required columns: ", paste(missing, collapse = ", "))
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  all_results <- list()

  for (i in seq_len(nrow(locations))) {
    location <- locations[i, ]
    params <- list(
      AWC = param_value(location, "AWC"),
      DDC = param_value(location, "DDC"),
      RCN = param_value(location, "RCN"),
      RZD = param_value(location, "RZD"),
      WUC = param_value(location, "WUC")
    )

    weather_format <- if ("weather_format" %in% names(locations)) location$weather_format else NULL
    weather <- load_weather(location$weather_file, weather_format)
    result <- compute_arid(
      weather = weather,
      latitude = as.numeric(location$latitude),
      elevation = as.numeric(location$elevation),
      params = params
    )
    result$location_id <- location$location_id
    result <- result[, c("location_id", setdiff(names(result), "location_id"))]

    output_file <- file.path(output_dir, paste0(location$location_id, "_ARID.csv"))
    write.csv(result, output_file, row.names = FALSE)
    all_results[[location$location_id]] <- result
  }

  combined <- do.call(rbind, all_results)
  write.csv(combined, file.path(output_dir, "TX_cotton_ARID_all_locations.csv"), row.names = FALSE)
  invisible(combined)
}

args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0) {
  if (length(args) != 2) {
    stop("Usage: Rscript scripts/compute_arid_batch.R inputs/locations.csv outputs/arid")
  }

  run_batch(args[[1]], args[[2]])
}
