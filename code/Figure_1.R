#=====================================================================
# Faceted Bubble Map Script - CNS Publication Quality (with South Sea Inset)
#=====================================================================
library(tidyverse)
library(sf)
library(ggsci)
library(ggpubr)
library(patchwork)

# Set the working directory to the repository root before running.

dir.create("results", showWarnings = FALSE)

# 1. Load Data
df <- read_csv("data/Step4_Refugia_Full_Filtered_NotVeg.csv")
glimpse(df)

# Prepare data
df$sale_form <- factor(df$propagule_type, levels = c("seed", "adult"))
df$year_label <- factor(df$year_label)
df$is_refugia_label <- ifelse(df$is_refugia == 1, "Refugia", "In-niche")
df$prop_label <- ifelse(df$sale_form == "seed", "Seed", "Live")

# 3. Load Maps
world_shp <- st_read("./data/worldmap/World_Land_Clean.shp", quiet = TRUE)
china_shp <- st_read("./data/China_map/bou2_4p.shp", quiet = TRUE)
south_sea_shp <- st_read("./data/China_map/south_sea.shp", quiet = TRUE)

if (is.na(st_crs(world_shp))) st_crs(world_shp) <- 4326
if (is.na(st_crs(china_shp))) st_crs(china_shp) <- 4326
if (is.na(st_crs(south_sea_shp))) st_crs(south_sea_shp) <- 4326

# 4. Aggregate
df_agg <- df %>%
  group_by(wcg_lng, wcg_lat, shop_id, is_refugia_label, year_label, prop_label) %>%
  summarise(n_species = n_distinct(taxon_name), .groups = "drop") %>%
  st_as_sf(coords = c("wcg_lng", "wcg_lat"), crs = 4326)

write_csv(st_drop_geometry(df_agg), "results/df_agg.csv")


# 5. CNS Palette
cns_colors <- c("Refugia" = "#5069A8", "In-niche" = "#A92E2E")

# 6. Annotation Data
annot_stats <- df %>%
  group_by(year_label, prop_label) %>%
  summarise(
    total_stores = n_distinct(shop_id),
    refugia_stores = n_distinct(shop_id[which(is_refugia == 1)]),
    refugia_pct = round(refugia_stores / total_stores * 100, 1),
    n_species = n_distinct(taxon_name),
    .groups = "drop"
  )

df_agg$prop_label <- factor(df_agg$prop_label, levels = c("Seed", "Live"))
annot_stats$prop_label <- factor(annot_stats$prop_label, levels = c("Seed", "Live"))
df_agg <- df_agg %>% arrange(factor(is_refugia_label, levels = c("Refugia", "In-niche")))

# 7. Create Main Map
p_main <- ggplot() +
  geom_sf(data = world_shp, fill = "#F8F9F9", color = "#D0D3D4", linewidth = 0.1) +
  geom_sf(data = china_shp, fill = "white", color = "#566573", linewidth = 0.3) +
  geom_sf(data = south_sea_shp, color = "#566573", linewidth = 0.3) +
  geom_sf(data = df_agg, aes(size = n_species, fill = is_refugia_label),
          shape = 21, alpha = 0.5, color = "black", stroke = 0.1) +
  scale_size_area(max_size = 12, name = "No. taxa", breaks = c(5, 10, 20, 40)) +
  scale_fill_manual(values = cns_colors, guide = "none") +
  coord_sf(xlim = c(73, 135), ylim = c(18, 54), expand = FALSE) +
  scale_x_continuous(breaks = seq(80, 130, 10)) +
  scale_y_continuous(breaks = seq(20, 50, 10)) +
  labs(y = "Latitude", x = "Longitude") +
  facet_grid(year_label ~ prop_label) +
  theme_bw() +
  theme(
    text = element_text(family = "serif", size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
    legend.position = c(0.56, 0.15),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


# Add stats to main map
p_main <- p_main +
  geom_text(data = annot_stats,
            aes(x = 75, y = 50,
                label = paste0("Nurseries: ", total_stores, "\n",
                             "Taxa: ", n_species, "\n",
                             "Refugia ratio: ", refugia_pct, "%")),
            hjust = 0, size = 4, family = "serif")



# 8. Create South Sea Inset Map
cat("Creating South Sea inset...\n")
p_inset <- ggplot() +
  geom_sf(data = world_shp, fill = "#F8F9F9", color = "#D0D3D4", linewidth = 0.1) +
  geom_sf(data = china_shp, fill = "white", color = "#566573", linewidth = 0.3) +
  geom_sf(data = south_sea_shp, color = "#566573", linewidth = 0.3) +
  coord_sf(xlim = c(106, 122), ylim = c(3, 25)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.background = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# 9. Combine with Patchwork
final_plot <- p_main +
  inset_element(p_inset, left = 0.88, bottom = 0.1, right = 0.95, top = 0.3, align_to = 'full')

# 10. Export.
ggsave("results/Figure_1.png", final_plot, device = ragg::agg_png,
       width = 10, height = 1013 / 126.7, units = "in", dpi = 600, bg = "white")
ggsave("results/Figure_1.pdf", final_plot, device = grDevices::cairo_pdf,
       width = 10, height = 1013 / 126.7, units = "in",
       family = "serif", bg = "white")
