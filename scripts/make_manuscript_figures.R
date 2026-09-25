#!/usr/bin/env Rscript

# Create manuscript figures and summary tables from the completed statewide
# Texas cotton ARID/phenology/yield run.

required_packages <- c("terra", "sf", "ggplot2", "dplyr", "viridis")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Please install required R package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript scripts/make_manuscript_figures.R OUTPUT_DIR COTTON_AREA_GEOJSON STATION_CATALOG_CSV", call. = FALSE)
}

output_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
area_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
station_catalog_file <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
figure_dir <- file.path(output_dir, "manuscript_figures")
table_dir <- file.path(output_dir, "manuscript_tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

theme_manuscript <- function() {
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.25),
      axis.title = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      strip.text = ggplot2::element_text(face = "bold", color = "black"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "grey25")
    )
}

save_plot <- function(plot, filename, width, height) {
  ggplot2::ggsave(
    filename = file.path(figure_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

area <- sf::st_read(area_file, quiet = TRUE)
if (is.na(sf::st_crs(area))) sf::st_crs(area) <- 4326
area <- sf::st_make_valid(sf::st_transform(area, 4326))
stations <- read.csv(station_catalog_file, stringsAsFactors = FALSE)
stations_sf <- sf::st_as_sf(stations, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

phase_summary <- read.csv(file.path(output_dir, "phenology_texas_cotton_summary.csv"), stringsAsFactors = FALSE)
yield_summary <- read.csv(file.path(output_dir, "relative_yield_summary.csv"), stringsAsFactors = FALSE)
enso <- read.csv(file.path(output_dir, "enso_year_classification.csv"), stringsAsFactors = FALSE)

phase_order <- c(
  "planting_emergence", "emergence_pinhead_square", "pinhead_square_first_bloom",
  "first_bloom_peak_bloom", "peak_bloom_first_open_boll", "first_open_boll_harvest"
)
phase_labels <- c(
  planting_emergence = "P1: Planting-emergence",
  emergence_pinhead_square = "P2: Emergence-pinhead square",
  pinhead_square_first_bloom = "P3: Pinhead square-first bloom",
  first_bloom_peak_bloom = "P4: First bloom-peak bloom",
  peak_bloom_first_open_boll = "P5: Peak bloom-first open boll",
  first_open_boll_harvest = "P6: First open boll-harvest"
)
phase_summary$phase_label <- factor(phase_labels[phase_summary$phenology_phase], levels = unname(phase_labels[phase_order]))
phase_summary$enso_phase <- factor(phase_summary$enso_phase, levels = c("El Nino", "La Nina", "Neutral"))
yield_summary$enso_phase <- factor(yield_summary$enso_phase, levels = c("El Nino", "La Nina", "Neutral"))

# Figure 1: cotton footprint and NOAA station network.
fig1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = area, fill = "#dcefd4", color = "#4f7d4f", linewidth = 0.25) +
  ggplot2::geom_sf(data = stations_sf, color = "#b21f35", size = 0.8, alpha = 0.75) +
  ggplot2::coord_sf(expand = FALSE) +
  ggplot2::labs(
    title = "Texas cotton analysis footprint and NOAA station network",
    subtitle = "CDL cotton pixels aggregated to the approximately 5-km analysis footprint",
    x = NULL, y = NULL
  ) +
  theme_manuscript() +
  ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
save_plot(fig1, "figure_1_cotton_area_stations.png", 7.2, 5.0)

# Figure 2: annual cotton-area mean relative yield loss.
enso_colors <- c("El Nino" = "#c43c39", "La Nina" = "#2b6ca3", "Neutral" = "#666666")
phase_means <- yield_summary |>
  dplyr::group_by(enso_phase) |>
  dplyr::summarise(mean_loss = mean(county_mean_relative_yield_loss, na.rm = TRUE), .groups = "drop")
fig2 <- ggplot2::ggplot(yield_summary, ggplot2::aes(year, county_mean_relative_yield_loss)) +
  ggplot2::geom_line(color = "#222222", linewidth = 0.45) +
  ggplot2::geom_point(ggplot2::aes(color = enso_phase), size = 1.5) +
  ggplot2::geom_hline(data = phase_means, ggplot2::aes(yintercept = mean_loss, color = enso_phase), linetype = "dashed", linewidth = 0.45, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = enso_colors) +
  ggplot2::scale_x_continuous(breaks = seq(1975, 2025, by = 5)) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  ggplot2::labs(
    title = "Annual modeled cotton relative yield loss",
    subtitle = "Points are colored by annual ENSO classification; dashed lines show ENSO-group means",
    x = "Year", y = "Relative yield loss", color = "ENSO phase"
  ) +
  theme_manuscript()
save_plot(fig2, "figure_2_annual_yield_loss.png", 7.2, 4.2)

# Build annual ENSO composite relative-yield-loss rasters.
loss_dir <- file.path(output_dir, "relative_yield_loss")
loss_files <- list.files(loss_dir, pattern = "^relative_yield_loss_[0-9]{4}\\.tif$", full.names = TRUE)
loss_years <- as.integer(sub("^relative_yield_loss_([0-9]{4})\\.tif$", "\\1", basename(loss_files)))
loss_lookup <- data.frame(year = loss_years, file = loss_files, stringsAsFactors = FALSE) |>
  dplyr::left_join(enso[, c("year", "phase")], by = "year")
composite_dir <- file.path(output_dir, "enso_composites")
dir.create(composite_dir, recursive = TRUE, showWarnings = FALSE)

composite_files <- list()
for (phase in c("El Nino", "La Nina", "Neutral")) {
  selected <- loss_lookup$file[loss_lookup$phase == phase]
  if (length(selected) == 0) next
  composite <- terra::app(terra::rast(selected), fun = mean, na.rm = TRUE)
  composite_file <- file.path(composite_dir, paste0("mean_relative_yield_loss_", gsub(" ", "_", phase), ".tif"))
  terra::writeRaster(composite, composite_file, overwrite = TRUE)
  composite_files[[phase]] <- composite_file
}

map_data <- do.call(rbind, lapply(names(composite_files), function(phase) {
  r <- terra::rast(composite_files[[phase]])
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "value"
  d$enso_phase <- phase
  d
}))
map_data$enso_phase <- factor(map_data$enso_phase, levels = c("El Nino", "La Nina", "Neutral"))
map_range <- range(map_data$value, na.rm = TRUE)
fig3 <- ggplot2::ggplot(map_data, ggplot2::aes(x, y, fill = value)) +
  ggplot2::geom_raster() +
  ggplot2::geom_sf(data = area, inherit.aes = FALSE, fill = NA, color = "grey15", linewidth = 0.25) +
  ggplot2::facet_wrap(~enso_phase, nrow = 1) +
  ggplot2::coord_sf(expand = FALSE) +
  viridis::scale_fill_viridis(option = "magma", direction = -1, limits = map_range, name = "Yield loss") +
  ggplot2::labs(title = "ENSO composites of modeled cotton relative yield loss", subtitle = "Mean annual 5-km rasters within each ENSO class, 1975-2025", x = NULL, y = NULL) +
  theme_manuscript() +
  ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
save_plot(fig3, "figure_3_enso_yield_loss_maps.png", 9.0, 4.2)

# Figure 4: ENSO departures from Neutral years.
neutral_raster <- terra::rast(composite_files[["Neutral"]])
difference_files <- list()
for (phase in c("El Nino", "La Nina")) {
  difference <- terra::rast(composite_files[[phase]]) - neutral_raster
  difference_file <- file.path(composite_dir, paste0("difference_", gsub(" ", "_", phase), "_minus_Neutral_yield_loss.tif"))
  terra::writeRaster(difference, difference_file, overwrite = TRUE)
  difference_files[[phase]] <- difference_file
}
diff_data <- do.call(rbind, lapply(names(difference_files), function(phase) {
  r <- terra::rast(difference_files[[phase]])
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "value"
  d$comparison <- paste(phase, "minus Neutral")
  d
}))
diff_limit <- max(abs(diff_data$value), na.rm = TRUE)
fig4 <- ggplot2::ggplot(diff_data, ggplot2::aes(x, y, fill = value)) +
  ggplot2::geom_raster() +
  ggplot2::geom_sf(data = area, inherit.aes = FALSE, fill = NA, color = "grey15", linewidth = 0.25) +
  ggplot2::facet_wrap(~comparison, nrow = 1) +
  ggplot2::coord_sf(expand = FALSE) +
  ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0, limits = c(-diff_limit, diff_limit), name = "Loss difference") +
  ggplot2::labs(title = "ENSO departures from Neutral relative yield loss", subtitle = "Positive values indicate greater modeled loss than Neutral years", x = NULL, y = NULL) +
  theme_manuscript() +
  ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
save_plot(fig4, "figure_4_enso_difference_maps.png", 7.2, 4.2)

# Figure 5: phase-specific ARID by ENSO class.
phase_heat <- phase_summary |>
  dplyr::group_by(enso_phase, phase_label) |>
  dplyr::summarise(mean_arid = mean(county_mean_arid, na.rm = TRUE), .groups = "drop")
fig5 <- ggplot2::ggplot(phase_heat, ggplot2::aes(phase_label, enso_phase, fill = mean_arid)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.35) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", mean_arid)), size = 3.1) +
  viridis::scale_fill_viridis(option = "cividis", limits = range(phase_heat$mean_arid, na.rm = TRUE), name = "Mean ARID") +
  ggplot2::labs(title = "Phase-specific cotton drought stress by ENSO class", subtitle = "Cotton-area mean phase ARID across 1975-2025", x = NULL, y = NULL) +
  theme_manuscript() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1))
save_plot(fig5, "figure_5_phase_enso_heatmap.png", 8.4, 3.8)

# Manuscript-ready tables.
data_sources <- data.frame(
  source = c("NOAA nClimGrid-Daily", "gridMET", "NOAA GHCN-Daily", "USDA Cropland Data Layer"),
  variables_or_role = c("Precipitation, minimum temperature, maximum temperature", "Solar radiation, wind speed, specific humidity", "Station precipitation and temperature for bias correction", "Cotton class 2 analysis footprint"),
  temporal_coverage = c("1975-2025", "1979-2025; fallback estimates for 1975-1978", "Station-dependent records within 1975-2025", "Layer supplied by the user"),
  spatial_resolution = c("Approximately 5 km", "Approximately 4 km", "Point observations", "30 m source, aggregated to approximately 5 km"),
  processing_role = c("Continuous gridded baseline", "Supplemental gridded weather", "Monthly bias correction of PRCP, TMIN, and TMAX", "Defines valid cotton analysis cells"),
  stringsAsFactors = FALSE
)
write.csv(data_sources, file.path(table_dir, "table_1_data_sources.csv"), row.names = FALSE)

phase_parameters <- data.frame(
  phase = paste0("P", 1:6),
  phenological_interval = unname(phase_labels[phase_order]),
  thermal_time_Cd = c(60, 230, 225, 200, 370, 200),
  equation_8_exponent = c(0.01, -0.11, 0.16, 0.09, 0.06, 0.08),
  stringsAsFactors = FALSE
)
write.csv(phase_parameters, file.path(table_dir, "table_2_phenological_parameters.csv"), row.names = FALSE)

enso_summary <- yield_summary |>
  dplyr::group_by(enso_phase) |>
  dplyr::summarise(
    years = dplyr::n(),
    mean_relative_yield = mean(county_mean_relative_yield, na.rm = TRUE),
    sd_relative_yield = sd(county_mean_relative_yield, na.rm = TRUE),
    mean_relative_yield_loss = mean(county_mean_relative_yield_loss, na.rm = TRUE),
    sd_relative_yield_loss = sd(county_mean_relative_yield_loss, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(enso_summary, file.path(table_dir, "table_3_enso_yield_summary.csv"), row.names = FALSE)

diagnostics <- read.csv(file.path(output_dir, "diagnostics", "station_overlap_diagnostics.csv"), stringsAsFactors = FALSE)
diagnostic_summary <- data.frame(
  variable = c("Precipitation", "Minimum temperature", "Maximum temperature"),
  mean_mae_raw = c(mean(diagnostics$prcp_mae_raw, na.rm = TRUE), mean(diagnostics$tmin_mae_raw, na.rm = TRUE), mean(diagnostics$tmax_mae_raw, na.rm = TRUE)),
  mean_mae_corrected = c(mean(diagnostics$prcp_mae_corrected, na.rm = TRUE), mean(diagnostics$tmin_mae_corrected, na.rm = TRUE), mean(diagnostics$tmax_mae_corrected, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
write.csv(diagnostic_summary, file.path(table_dir, "table_4_bias_correction_summary.csv"), row.names = FALSE)

message("Wrote figures to: ", figure_dir)
message("Wrote tables to: ", table_dir)
