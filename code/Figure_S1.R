# Figure S1: sale-form proportions by year and status.

# 1. Packages and paths -----------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
})
# Set the working directory to the repository root before running.
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)
input_dir <- "data"
csv_path <- file.path(input_dir, "Step4_Refugia_Full_Filtered_NotVeg.csv")

# 2. Group labels -----------------------------------
status_labels <- c(`1` = "Non-naturalized", `2` = "Naturalized non-invasive", `3` = "Invasive")
form_levels <- c("Seed only", "Live plant only", "Both", "Unclassified")

# 3. Read and validate records -----------------------------------
trade <- read_csv(csv_path, show_col_types = FALSE, col_types = cols(.default = col_character())) %>%
  transmute(year = factor(year_label, levels = c("2019", "2022")), plant_name_id = str_trim(accepted_plant_name_id),
    status_code = str_trim(invasion.status), seed = as.integer(propagule_type ==
      "seed"), adult = as.integer(propagule_type == "adult"))
stopifnot(nrow(trade) == 21351L, n_distinct(trade$plant_name_id) == 748L, !anyNA(trade$plant_name_id),
  all(trade$plant_name_id != ""), all(trade$status_code %in% names(status_labels)),
  all(trade$seed %in% 0:1), all(trade$adult %in% 0:1))
status_check <- trade %>%
  distinct(year, plant_name_id, status_code) %>%
  count(year, plant_name_id)
stopifnot(all(status_check$n == 1L))

# 4. Summarize sale forms -----------------------------------
species_data <- trade %>%
  mutate(invasion_status = factor(status_labels[status_code], levels = status_labels)) %>%
  group_by(year, invasion_status, plant_name_id) %>%
  summarise(seed = any(seed == 1), live_plant = any(adult == 1), .groups = "drop") %>%
  mutate(sale_category = factor(case_when(seed & live_plant ~ "Both", seed ~
    "Seed only", live_plant ~ "Live plant only", TRUE ~ "Unclassified"), levels = form_levels))
species_totals <- species_data %>%
  count(year, invasion_status, name = "n_species")
form_summary <- species_data %>%
  count(year, invasion_status, sale_category, name = "n_species") %>%
  group_by(year, invasion_status) %>%
  mutate(proportion = n_species/sum(n_species), label = scales::number(proportion *
    100, accuracy = 0.1)) %>%
  ungroup()
form_check <- form_summary %>%
  summarise(total = sum(proportion), .by = c(year, invasion_status))
stopifnot(all(abs(form_check$total - 1) < 1e-10))
# Compare all 18 category counts and proportions before replacing the source
# table.
previous_file <- file.path(output_dir, "Figure_S1_sale_form_source_data.csv")
if (file.exists(previous_file)) {
  previous <- read_csv(previous_file, show_col_types = FALSE) %>%
    mutate(across(c(year, invasion_status, sale_category), as.character))
  current <- form_summary %>%
    mutate(across(c(year, invasion_status, sale_category), as.character))
  comparison <- full_join(previous, current, by = c("year", "invasion_status",
    "sale_category"), suffix = c("_previous", "_current"))
  stopifnot(nrow(comparison) == 18L, !anyNA(comparison$n_species_previous), !anyNA(comparison$n_species_current),
    all(comparison$n_species_previous == comparison$n_species_current), all(abs(comparison$proportion_previous -
      comparison$proportion_current) < 1e-12))
}

# 5. Colours and theme -----------------------------------
npg <- (ggsci::pal_npg("nrc"))(10)
form_colours <- c(`Seed only` = npg[2], `Live plant only` = npg[3], Both = npg[5],
  Unclassified = "#BDBDBD")
my_theme <- theme_classic(base_size = 16, base_family = "serif") + theme(text = element_text(face = "plain"),
  strip.background = element_blank(), strip.text = element_text(face = "plain",
    size = 16), axis.text = element_text(colour = "black", size = 14), legend.text = element_text(size = 14),
  legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(10,
    10, 10, 10))

# 6. Draw the figure -----------------------------------
p_forms <- ggplot(form_summary, aes(year, proportion, fill = sale_category)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = label), size = 4.5, family = "serif", lineheight = 1,
    position = position_stack(vjust = 0.5, reverse = TRUE)) + geom_text(data = species_totals,
  aes(x = year, y = 1.055, label = paste0("N = ", n_species)), inherit.aes = FALSE,
  size = 4.5, family = "serif") + scale_fill_manual(values = form_colours) +
  scale_y_continuous(labels = scales::label_number(scale = 100), breaks = seq(0,
    1, 0.25), limits = c(0, 1.11), expand = c(0, 0)) + facet_wrap(~invasion_status,
  nrow = 1) + labs(x = NULL, y = "Taxa (%)") + my_theme

# 7. Export -----------------------------------
figure_base <- file.path(output_dir, "Figure_S1")
ggsave(paste0(figure_base, ".png"), p_forms, device = ragg::agg_png, width = 300,
  height = 115, units = "mm", dpi = 600, bg = "white")
ggsave(paste0(figure_base, ".pdf"), p_forms, device = grDevices::cairo_pdf, width = 300,
  height = 115, units = "mm", family = "serif", bg = "white")
write_csv(form_summary, file.path(output_dir, "Figure_S1_sale_form_source_data.csv"))
write_csv(species_totals, file.path(output_dir, "Figure_S1_group_totals.csv"))
message("Figure S1 completed: sale forms only; all 18 category values preserved.")
print(species_totals)
