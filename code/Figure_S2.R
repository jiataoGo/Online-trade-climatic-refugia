# Figure S2: taxon richness and mean suitable area.
# Rebuild data: Rscript code/Figure_S2.R --recalculate

# 1. Paths. Set the working directory to the repository root before running.
OUTPUT_DIR <- "results"
PROJECT_ROOT <- "data"

# 2. Set paths and styles. nat = naturalized non-invasive.
CONFIG <- list(
  root = PROJECT_ROOT,
  output_dir = OUTPUT_DIR,
  source_csv = file.path(PROJECT_ROOT, "Step4_Refugia_Full_Filtered_NotVeg.csv"),
  raster_dir = file.path(PROJECT_ROOT, "SDM_combined"),
  map_dir = file.path(PROJECT_ROOT, "China_map"),
  png_file = file.path(OUTPUT_DIR, "Figure_S2.png"),

  recalculate = "--recalculate" %in% commandArgs(trailingOnly = TRUE),
  expected_species = 748L,
  copy_selected = FALSE,
  write_geotiffs = TRUE,
  write_tables = TRUE,

  main_xlim = c(72, 143),
  main_ylim = c(16, 55),
  inset_xlim = c(105, 125.5),
  inset_ylim = c(3, 24),
  inset_placement = c(left = 0.79, bottom = 0.025, right = 0.995, top = 0.40),
  area_xlim = c(0, 2.06),
  figure_width_mm = 240,
  figure_height_mm = 205,
  dpi = 600,
  font_family = "Arial",
  base_font_size = 12,
  zero_grey = "#C9C9C9",
  nodata_grey = "#E9E9E9",
  richness_palette = c("#E0F1DA", "#B7DFAB", "#78C679", "#31A354", "#08783E", "#004D2A"),
  class_colors = c(inv = "#547EA0", nat = "#547EA0", nonnat = "#547EA0"),
  class_labels = c(inv = "Invasive", nat = "Naturalized non-invasive", nonnat = "Non-naturalized")
)

required_packages <- c("terra", "readr", "sf", "ggplot2", "patchwork", "ragg")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install the missing R package before running this script: ", package)
  }
}

# Load validated data or recalculate from rasters.
# PART 1. DATA PREPARATION ---------------------------------------------


config <- CONFIG
if (config$recalculate) {
  # Recalculate from source rasters.

  # 1. Check input paths.
  if (!file.exists(config$source_csv)) stop("The source CSV does not exist.")
  if (!dir.exists(config$raster_dir)) stop("The source raster directory does not exist.")
  output_dir <- config$output_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  classes <- c("inv", "nat", "nonnat")
  class_english <- c(inv = "Invasive", nat = "Naturalized non-invasive",
                     nonnat = "Non-naturalized")
  # 2. Read records and identify unique taxa.
  records <- as.data.frame(readr::read_csv(
    config$source_csv, col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE), stringsAsFactors = FALSE)
  species_column <- if ("species_dir" %in% names(records)) "species_dir" else "species"
  required_columns <- c(species_column, "invasion_label", "presence_cells", "absence_cells")
  missing_columns <- setdiff(required_columns, names(records))
  if (length(missing_columns) > 0L) {
    stop("Source CSV is missing columns: ", paste(missing_columns, collapse = ", "))
  }
  records$species_dir <- trimws(records[[species_column]])
  records$invasion_label <- trimws(records$invasion_label)
  if (!all(!is.na(records$species_dir) & nzchar(records$species_dir))) {
    stop("Species identifiers contain missing or empty values.")
  }
  if (!all(!is.na(records$invasion_label) & nzchar(records$invasion_label))) {
    stop("Invasion labels contain missing or empty values.")
  }
  if (!setequal(unique(records$invasion_label), classes)) {
    stop("Expected exactly three invasion labels: inv, nat and nonnat.")
  }
  species <- unique(records[c("species_dir", "invasion_label")])
  if (anyDuplicated(species$species_dir)) {
    stop("A species has conflicting invasion labels in the source CSV.")
  }
  species <- species[order(species$species_dir), , drop = FALSE]
  rownames(species) <- NULL
  if (nrow(species) != config$expected_species) {
    stop(sprintf("Expected %d unique species; found %d.",
                 config$expected_species, nrow(species)))
  }
  record_counts <- table(records$species_dir)
  species$csv_record_count <- as.integer(record_counts[species$species_dir])

  # Validate CSV cell counts against raster values.
  expected_counts <- unique(records[c("species_dir", "presence_cells", "absence_cells")])
  if (anyDuplicated(expected_counts$species_dir)) {
    stop("A species has conflicting presence/absence counts in the source CSV.")
  }
  for (column in c("presence_cells", "absence_cells")) {
    count_values <- suppressWarnings(as.numeric(expected_counts[[column]]))
    if (!all(is.finite(count_values) & count_values >= 0 & count_values == floor(count_values))) {
      stop("Invalid nonnegative integer values in ", column, ".")
    }
    expected_counts[[column]] <- count_values
  }
  count_index <- match(species$species_dir, expected_counts$species_dir)

  # 3. Match exact raster filenames to the CSV.
  files <- list.files(config$raster_dir, pattern = "_combined_binary[.]tif$",
                      recursive = TRUE, full.names = TRUE)
  files <- files[basename(dirname(files)) == "EnsemblePlot"]
  identifiers <- sub("_combined_binary[.]tif$", "", basename(files))
  if (anyDuplicated(identifiers)) {
    stop("Duplicate combined TIFF filenames were found in the source directory.")
  }
  species$source_path <- files[match(species$species_dir, identifiers)]
  species$file_status <- ifelse(is.na(species$source_path), "missing", "matched")
  extra <- !identifiers %in% species$species_dir
  if (any(species$file_status == "missing")) {
    if (config$write_tables) {
      readr::write_excel_csv(species[species$file_status == "missing", , drop = FALSE],
                             file.path(output_dir, "missing_species.csv"))
    }
    stop("Some requested combined TIFFs are missing; inspect missing_species.csv.")
  }

  # 4. Calculate geographic cell areas in km².
  template <- terra::rast(species$source_path[1])
  if (!isTRUE(terra::nlyr(template) == 1L && terra::is.lonlat(template))) {
    stop("The template must be a single-layer longitude/latitude raster.")
  }
  domain <- !is.na(terra::values(template, mat = FALSE))
  if (!any(domain)) stop("The template has no valid model cells.")
  cell_area <- terra::cellSize(template, mask = FALSE, unit = "km")
  areas <- terra::values(cell_area, mat = FALSE)
  if (!all(is.finite(areas) & areas > 0)) stop("Cell areas must be finite and positive.")
  domain_area <- sum(areas[domain])
  richness <- matrix(0L, nrow = terra::ncell(template), ncol = length(classes),
                     dimnames = list(NULL, classes))
  species$presence_cells <- species$absence_cells <- species$nodata_cells <- 0L
  species$suitable_area_km2 <- species$valid_area_km2 <- NA_real_
  species$mask_difference_cells <- 0L
  species$source_size_bytes <- file.info(species$source_path)$size
  species$source_md5 <- unname(tools::md5sum(species$source_path))
  if (anyNA(species$source_md5)) stop("Could not calculate every source-file MD5.")
  species$selected_copy <- NA_character_
  species$copy_md5_match <- NA
  species$csv_presence_cells <- expected_counts$presence_cells[count_index]
  species$csv_absence_cells <- expected_counts$absence_cells[count_index]
  if (config$copy_selected) {
    for (label in classes) {
      dir.create(file.path(output_dir, "selected_binary", label),
                 recursive = TRUE, showWarnings = FALSE)
    }
  }

  # 5. Aggregate rasters: 1 = suitable, 0 = unsuitable, NA = unmodelled.
  for (i in seq_len(nrow(species))) {
    raster <- terra::rast(species$source_path[i])
    species_name <- species$species_dir[i]
    if (!isTRUE(terra::nlyr(raster) == 1L &&
                terra::compareGeom(template, raster, stopOnError = FALSE))) {
      stop("Raster grid mismatch: ", species_name)
    }
    raster_values <- terra::values(raster, mat = FALSE)
    valid <- !is.na(raster_values)
    if (!all(raster_values[valid] %in% c(0, 1))) {
      stop("Non-binary raster values: ", species_name)
    }
    present <- valid & raster_values == 1
    class_index <- match(species$invasion_label[i], classes)
    richness[, class_index] <- richness[, class_index] + as.integer(present)
    species$presence_cells[i] <- sum(present)
    species$absence_cells[i] <- sum(valid & raster_values == 0)
    species$nodata_cells[i] <- sum(!valid)
    species$suitable_area_km2[i] <- sum(areas[present])
    species$valid_area_km2[i] <- sum(areas[valid])
    species$mask_difference_cells[i] <- sum(valid != domain)

    if (config$copy_selected) {
      destination <- file.path(output_dir, "selected_binary", species$invasion_label[i],
                               basename(species$source_path[i]))
      if (!file.copy(species$source_path[i], destination, overwrite = TRUE)) {
        stop("Could not copy the selected raster: ", species_name)
      }
      species$selected_copy[i] <- destination
      species$copy_md5_match[i] <- unname(tools::md5sum(destination)) == species$source_md5[i]
    }
    if (i %% 50L == 0L || i == nrow(species)) {
      message(sprintf("Validated and aggregated %d / %d species.", i, nrow(species)))
    }
  }
  species$csv_counts_match <- species$presence_cells == species$csv_presence_cells &
    species$absence_cells == species$csv_absence_cells

  # Require identical NoData masks.
  if (!all(species$mask_difference_cells == 0L)) {
    stop("Species have different NoData masks; aggregation stopped to protect missing-data semantics.")
  }
  if (config$copy_selected) {
    if (!isTRUE(all(species$copy_md5_match))) stop("A selected copy failed its MD5 check.")
    copied <- list.files(file.path(output_dir, "selected_binary"),
                         pattern = "_combined_binary[.]tif$", recursive = TRUE, full.names = TRUE)
    expected_keys <- paste(species$invasion_label, basename(species$source_path), sep = "/")
    copied_keys <- paste(basename(dirname(copied)), basename(copied), sep = "/")
    if (!isTRUE(length(copied_keys) == nrow(species) && setequal(copied_keys, expected_keys))) {
      stop("selected_binary contains unexpected or stale TIFFs; inspect it before rerunning.")
    }
  }
  if (!all(species$csv_counts_match)) {
    warning("Some CSV pixel counts differ from the source TIFFs; the analysis uses TIFF values.",
            call. = FALSE)
  }
  richness[!domain, ] <- NA_integer_
  total_richness <- rowSums(richness)
  class_presence <- richness > 0
  # Encode group overlap: inv = 1, nat = 2, nonnat = 4.
  overlap_code <- as.integer(class_presence[, 1]) +
    2L * as.integer(class_presence[, 2]) + 4L * as.integer(class_presence[, 3])

  # 6. Summarize union and per-taxon areas, including zeros.
  # SD measures variation among taxa.
  summary_rows <- vector("list", length(classes))
  for (j in seq_along(classes)) {
    label <- classes[j]
    present <- !is.na(overlap_code) & class_presence[, j]
    exclusive <- !is.na(overlap_code) & overlap_code == 2^(j - 1L)
    species_areas <- species$suitable_area_km2[species$invasion_label == label]
    summary_rows[[j]] <- data.frame(
      invasion_label = label,
      class_en = unname(class_english[label]), species_count = length(species_areas),
      occupied_cells = sum(present), union_area_km2 = sum(areas[present]),
      coverage_pct_of_valid_grid = 100 * sum(areas[present]) / domain_area,
      exclusive_area_km2 = sum(areas[exclusive]),
      max_species_per_grid = max(richness[, j], na.rm = TRUE),
      area_weighted_mean_species_per_grid = sum(richness[domain, j] * areas[domain]) / domain_area,
      summed_species_area_km2 = sum(species_areas),
      mean_species_area_km2 = mean(species_areas),
      median_species_area_km2 = stats::median(species_areas),
      sd_species_area_km2 = stats::sd(species_areas),
      min_species_area_km2 = min(species_areas),
      q25_species_area_km2 = unname(stats::quantile(species_areas, 0.25)),
      q75_species_area_km2 = unname(stats::quantile(species_areas, 0.75)),
      max_species_area_km2 = max(species_areas), zero_area_species = sum(species_areas == 0),
      stringsAsFactors = FALSE)
  }
  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL
  mean_columns <- c("invasion_label", "class_en", "species_count",
                    "mean_species_area_km2", "median_species_area_km2", "sd_species_area_km2",
                    "min_species_area_km2", "q25_species_area_km2", "q75_species_area_km2",
                    "max_species_area_km2", "zero_area_species")
  mean_area <- summary[mean_columns]
  mean_area$mean_pct_of_valid_grid <- 100 * mean_area$mean_species_area_km2 / domain_area

  overlap <- data.frame(
    bitmask = 0:7, inv_present = as.integer(bitwAnd(0:7, 1L) != 0L),
    nat_present = as.integer(bitwAnd(0:7, 2L) != 0L),
    nonnat_present = as.integer(bitwAnd(0:7, 4L) != 0L),
    combination = c("none", "inv_only", "nat_only", "inv_nat", "nonnat_only",
                     "inv_nonnat", "nat_nonnat", "inv_nat_nonnat"),
    cells = 0L, area_km2 = 0)
  for (i in seq_len(nrow(overlap))) {
    selected <- !is.na(overlap_code) & overlap_code == overlap$bitmask[i]
    overlap$cells[i] <- sum(selected)
    overlap$area_km2[i] <- sum(areas[selected])
  }
  overlap$pct_of_valid_grid <- 100 * overlap$area_km2 / domain_area

  # 7. Cross-check area totals.
  stopifnot(sum(overlap$cells) == sum(domain))
  stopifnot(abs(sum(overlap$area_km2) - domain_area) < 1e-6)
  for (j in seq_along(classes)) {
    selected <- overlap[[paste0(classes[j], "_present")]] == 1L
    stopifnot(abs(sum(overlap$area_km2[selected]) - summary$union_area_km2[j]) < 1e-6)
    stopifnot(abs(sum(richness[domain, j] * areas[domain]) -
                  summary$summed_species_area_km2[j]) < 1e-4)
  }
  stopifnot(all(abs(summary$mean_species_area_km2 * summary$species_count -
                     summary$summed_species_area_km2) < 1e-4))

  grid <- as.data.frame(terra::xyFromCell(template, which(domain)))
  grid$grid_id <- which(domain)
  grid$cell_area_km2 <- areas[domain]
  for (label in classes) grid[[paste0(label, "_species_count")]] <- richness[domain, label]
  grid$all_species_count <- total_richness[domain]
  grid$overlap_bitmask <- overlap_code[domain]
  overall <- data.frame(
    csv_rows = nrow(records), unique_species = nrow(species),
    selected_combined_files = nrow(species), available_combined_files = length(files),
    extra_files_excluded = sum(extra), excluded_species = paste(sort(identifiers[extra]), collapse = "; "),
    valid_grid_cells = sum(domain), valid_grid_area_km2 = domain_area,
    all_union_cells = sum(total_richness > 0, na.rm = TRUE),
    all_union_area_km2 = sum(areas[!is.na(total_richness) & total_richness > 0]),
    all_max_species_per_grid = max(total_richness, na.rm = TRUE),
    csv_pixel_count_mismatches = sum(!species$csv_counts_match),
    differing_masks = sum(species$mask_difference_cells > 0),
    all_copies_md5_verified = if (config$copy_selected) all(species$copy_md5_match) else NA,
    nrow = terra::nrow(template), ncol = terra::ncol(template),
    x_resolution_degrees = terra::res(template)[1], y_resolution_degrees = terra::res(template)[2],
    crs = terra::crs(template, proj = TRUE), csv_md5 = unname(tools::md5sum(config$source_csv)),
    r_version = R.version.string, terra_version = as.character(utils::packageVersion("terra")),
    readr_version = as.character(utils::packageVersion("readr")), stringsAsFactors = FALSE)

  # 8. Export richness rasters, preserving NoData.
  if (config$write_geotiffs) {
    raster_output <- file.path(output_dir, "rasters")
    dir.create(raster_output, recursive = TRUE, showWarnings = FALSE)
    if (nrow(species) >= 65535L) stop("Species richness exceeds the output integer range.")
    for (label in classes) {
      output_name <- paste0(label, "_species_richness")
      output_raster <- terra::setValues(terra::rast(template), richness[, label])
      names(output_raster) <- output_name
      terra::writeRaster(output_raster, file.path(raster_output, paste0(output_name, ".tif")),
                         overwrite = TRUE,
                         wopt = list(datatype = "INT2U", NAflag = 65535, gdal = "COMPRESS=DEFLATE"))
    }
  }

  # 9. Save summary tables and the plotting cache.
  results <- list(grid = grid, summary = summary, mean_area = mean_area, overall = overall,
                  extent = as.vector(terra::ext(template)), resolution = terra::res(template),
                  overlap = overlap)
  if (config$write_tables) {
    table_outputs <- list(
      species_areas_and_checks.csv = species,
      class_mean_species_area.csv = mean_area,
      grid_species_counts.csv = grid,
      overall_summary_and_qa.csv = overall)
    for (filename in names(table_outputs)) {
      readr::write_excel_csv(table_outputs[[filename]], file.path(output_dir, filename))
    }
    saveRDS(results, file.path(PROJECT_ROOT, "plot_data.rds"))
  }
  message(sprintf("Analysis complete: %d species, %d valid cells, %.2f km^2 model domain.",
                  nrow(species), sum(domain), domain_area))
} else {
  # Reuse cached results for styling.
  cache_file <- file.path(PROJECT_ROOT, "plot_data.rds")
  if (!file.exists(cache_file)) stop("Run once with recalculate = TRUE first.")
  results <- readRDS(cache_file)
  stopifnot(results$overall$unique_species == config$expected_species)
  if (unname(tools::md5sum(config$source_csv)) != results$overall$csv_md5) {
    stop("The source CSV has changed. Set recalculate = TRUE.")
  }
  message("Styling-only run: cached raster results are being used.")
}


# PART 2. FIGURE -------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(sf)
})
grid_data <- results$grid
area_summary <- results$summary

font_scale <- config$base_font_size / 10

# 3. Read display boundaries; area estimates are unchanged.
china_raw <- st_read(file.path(config$map_dir, "bou1_4p.shp"),
                     options = "ENCODING=GBK", quiet = TRUE)
south <- st_read(file.path(config$map_dir, "south_sea.shp"), quiet = TRUE)
if (is.na(st_crs(china_raw))) st_crs(china_raw) <- 4326
china_raw <- st_transform(china_raw, 4326)
south <- st_transform(south, 4326)
dash <- south[south$GBCODE == 61010 & south$BOU1_4M_ != 691, ]
stopifnot(nrow(dash) == 9L)
islets <- south[south$GBCODE != 61010, ]
land_geometry <- st_union(st_make_valid(st_set_crs(st_geometry(china_raw), NA)))
frame <- st_as_sfc(st_bbox(c(xmin = 70, ymin = 1, xmax = 145, ymax = 57)))
outside_geometry <- st_difference(frame, land_geometry)
land <- st_sf(geometry = st_set_crs(land_geometry, 4326))
outside <- st_sf(geometry = st_set_crs(outside_geometry, 4326))

# Convert boundaries to paths, preserving polygon holes.
geometry_layers <- list(land = land, outside = outside, islets = islets, dash = dash)
map_paths <- list()
for (layer_name in names(geometry_layers)) {
  coordinates <- as.data.frame(st_coordinates(geometry_layers[[layer_name]]))
  id_columns <- grep("^L[0-9]+$", names(coordinates), value = TRUE)
  coordinates$group <- 1L
  coordinates$subgroup <- 1L
  for (column_name in id_columns) {
    coordinates$subgroup <- paste(coordinates$subgroup, coordinates[[column_name]], sep = "_")
  }
  group_columns <- id_columns
  if (layer_name %in% c("land", "outside")) group_columns <- id_columns[-1]
  for (column_name in group_columns) {
    coordinates$group <- paste(coordinates$group, coordinates[[column_name]], sep = "_")
  }
  map_paths[[layer_name]] <- coordinates
}
land_path <- map_paths$land
outside_path <- map_paths$outside
islet_path <- map_paths$islets
dash_path <- map_paths$dash

# 4. Set colours and text styles.
classes <- c("nonnat", "nat", "inv")
labels <- config$class_labels
bar_colors <- config$class_colors
continuous <- config$richness_palette
zero_grey <- config$zero_grey
nodata_grey <- config$nodata_grey
publication_theme <- theme_classic(base_size = config$base_font_size, base_family = config$font_family) +
  theme(axis.title = element_text(size = 10 * font_scale),
    axis.text = element_text(size = 9 * font_scale, colour = "black"),
    axis.line = element_line(linewidth = .3),
    axis.ticks = element_line(linewidth = .3),
    axis.ticks.length = grid::unit(1.2, "mm"),
    plot.title = element_text(size = 10.5 * font_scale, face = "plain", hjust = .5),
    plot.subtitle = element_text(size = 9.5 * font_scale, hjust = .5),
    plot.tag = element_text(size = 13 * font_scale, face = "plain"),
    plot.tag.position = "topleft",
    legend.title = element_text(size = 10 * font_scale),
    legend.text = element_text(size = 9 * font_scale),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(2, 2, 2, 2, "mm"))
# 5. Draw maps at the original raster resolution.
number_of_rows <- results$overall$nrow
number_of_columns <- results$overall$ncol
full_grid <- expand.grid(
  x = results$extent[1] + (seq_len(number_of_columns) - .5) * results$resolution[1],
  y = results$extent[4] - (seq_len(number_of_rows) - .5) * results$resolution[2])
common_max <- max(grid_data[, paste0(classes, "_species_count")])

maps <- list()
for (class_index in seq_along(classes)) {
  class_name <- classes[class_index]
  raster_data <- full_grid
  raster_data$richness <- NA_real_
  raster_data$richness[grid_data$grid_id] <- grid_data[[paste0(class_name, "_species_count")]]
  raster_data$is_zero <- as.numeric(!is.na(raster_data$richness) & raster_data$richness == 0)
  raster_data$positive_richness <- ifelse(!is.na(raster_data$richness) & raster_data$richness > 0, raster_data$richness, NA_real_)
  view_plots <- list()
  for (view_name in c("main", "inset")) {
    inset <- view_name == "inset"
    map_plot <- ggplot() +
      geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
        rule = "evenodd", fill = nodata_grey, colour = NA) +
      geom_raster(data = raster_data, aes(x, y, alpha = is_zero), fill = zero_grey,
                  interpolate = FALSE, show.legend = FALSE) +
      scale_alpha_identity(guide = "none") +
      geom_raster(data = raster_data, aes(x, y, fill = positive_richness),
                  interpolate = FALSE, show.legend = !inset) +
      scale_fill_gradientn(colours = continuous, limits = c(1, common_max),
        breaks = c(1, 50, 100, 150, 200, 250, common_max),
        name = "Predicted taxon richness", na.value = "transparent",
        guide = if (inset) "none" else guide_colourbar(title.position = "top",
          title.hjust = .5, barwidth = grid::unit(85, "mm"),
          barheight = grid::unit(3.5, "mm"), ticks.colour = "white")) +
      geom_polygon(data = outside_path, aes(X, Y, group = group, subgroup = subgroup),
        rule = "evenodd", fill = "white", colour = NA) +
      geom_polygon(data = land_path, aes(X, Y, group = group, subgroup = subgroup),
        rule = "evenodd", fill = NA, colour = "#4D5658", linewidth = if (inset) .23 else .28) +
      geom_path(data = islet_path, aes(X, Y, group = group),
        colour = "#586260", linewidth = if (inset) .23 else .15) +
      geom_path(data = dash_path, aes(X, Y, group = group),
        colour = "#313A3B", linewidth = if (inset) .55 else .35) +
      labs(x = NULL, y = NULL)
    if (inset) {
      map_plot <- map_plot + coord_fixed(ratio = 1 / cos(13.5 * pi / 180),
            xlim = config$inset_xlim, ylim = config$inset_ylim, expand = FALSE) +
        theme_void(base_family = config$font_family) +
        theme(panel.border = element_rect(fill = NA, colour = "#7D8484", linewidth = .3),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(0, 0, 0, 0), legend.position = "none")
    } else {
      map_plot <- map_plot + coord_fixed(ratio = 1 / cos(35.5 * pi / 180),
            xlim = config$main_xlim, ylim = config$main_ylim, expand = FALSE) +
        scale_x_continuous(breaks = c(80, 100, 120), labels = paste0(c(80, 100, 120), "°E")) +
        scale_y_continuous(breaks = c(20, 35, 50), labels = paste0(c(20, 35, 50), "°N")) +
        publication_theme +
        theme(legend.position = "bottom", axis.line = element_blank(),
          panel.border = element_rect(fill = NA, colour = "#AEB5B5", linewidth = .3))
    }
    view_plots[[view_name]] <- map_plot
  }
  # Use a shared extent and scale, with South China Sea insets.
  main <- view_plots$main + labs(tag = LETTERS[class_index], title = labels[class_name],
    subtitle = paste0("n = ", area_summary$species_count[area_summary$invasion_label == class_name]))
  position <- config$inset_placement
  maps[[class_index]] <- main + inset_element(view_plots$inset,
    left = position[["left"]], bottom = position[["bottom"]],
    right = position[["right"]], top = position[["top"]],
    align_to = "panel", on_top = TRUE)
}
# 6. Plot mean suitable area.
area_summary$invasion_label <- factor(area_summary$invasion_label, levels = classes)
area_labels <- c(inv = "Invasive", nat = "Naturalized\nnon-invasive",
                 nonnat = "Non-naturalized")
mean_area <- ggplot(area_summary, aes(x = invasion_label,
                           y = mean_species_area_km2 / 1e6, fill = invasion_label)) +
  geom_col(width = .58, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", mean_species_area_km2 / 1e6)),
            vjust = -.6, size = 3.4 * font_scale, family = config$font_family) +
  scale_fill_manual(values = bar_colors, guide = "none") +
  scale_x_discrete(labels = area_labels) +
  scale_y_continuous(limits = config$area_xlim, breaks = seq(0, 2, .5), expand = c(0, 0)) +
  labs(tag = "D", x = NULL,
       y = "Mean suitable area per taxon (10⁶ km²)") +
  publication_theme +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 9 * font_scale, lineheight = 1.05),
        plot.margin = margin(2, 4, 2, 4, "mm"))
# 7. Arrange panels with the shared legend above.
figure <- wrap_plots(maps[[1]], maps[[2]], maps[[3]],
                    wrap_elements(full = mean_area), guide_area(),
                    design = "EE\nAB\nCD", heights = c(.20, 1, 1), guides = "collect") &
  theme(legend.position = "top")
# 8. Export PNG (600 dpi) and PDF.
dir.create(dirname(config$png_file), recursive = TRUE, showWarnings = FALSE)
ragg::agg_png(config$png_file,
  width = config$figure_width_mm, height = config$figure_height_mm,
  units = "mm", res = config$dpi, background = "white")
print(figure)
dev.off()
grDevices::cairo_pdf(file.path(config$output_dir, "Figure_S2.pdf"),
  width = config$figure_width_mm / 25.4, height = config$figure_height_mm / 25.4,
  family = config$font_family, bg = "white")
print(figure)
invisible(dev.off())
message("Figure S2 saved as PDF and PNG: ", config$png_file)
message("Nine South China Sea segments: ", paste(dash$BOU1_4M_, collapse = ", "))

message("Workflow finished. All current outputs are in: ", config$output_dir)
