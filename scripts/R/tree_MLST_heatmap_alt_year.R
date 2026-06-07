# ============================================================================
# Arbre MLST + heatmap Altitude / Année
# Entrées : results/phylogeny/MLST_iqtree.treefile
#           results/metadata_genaspe_species_final.tsv
#           data/metadata_genaspe.csv
# Sorties : tree_MLST_heatmap_altitude_annee.png/.pdf
# ============================================================================

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(ape, phangorn, tidyverse, ggtree, ggtext, ggnewscale, readr, tibble, grid)

in_tree <- "results/phylogeny/MLST_iqtree.treefile"
species_file <- "results/metadata_genaspe_species_final.tsv"
tick_metadata_file <- "data/metadata_genaspe.csv"
out_png <- "tree_altitude_year_heatmaps.png"
out_pdf <- "tree_altitude_year_heatmaps.pdf"

samples_to_exclude <- c("ET78", "EM2069", "ET43", "ET76", "ET79", "ET110", "ET44")
target_colors <- c(
  "B. afzelii" = "#E41A1C", "B. burgdorferi sensu stricto" = "#f5a4e1",
  "B. garinii" = "#228B22", "B. miyamotoi" = "#FF8C00",
  "B. valaisiana" = "#8B4513", "B. lusitaniae" = "#b103fc",
  "B. spielmanii" = "#00A6A6", "B. bavariensis" = "#7A4EAB", "NA" = "#999999"
)
heatmap_colors <- c(
  "400 m" = "#FEE8C8", "800 m" = "#FDBB84", "1200 m" = "#E34A33", "1500 m" = "#7F0000",
  "2023" = "#80B1D3", "2024" = "#FB8072", "2025" = "#B3DE69", "Inconnue" = "#BDBDBD", "Inconnu" = "#D9D9D9"
)

clean_species <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  str_squish(str_replace_all(x, "\\s*\\(candidat\\)", ""))
}

extract_species_vec <- function(x) {
  x <- clean_species(x)
  if (x %in% c("", "NA", "inconnu")) return(character(0))
  x <- str_remove(x, regex("^co-infection\\s*-\\s*", ignore_case = TRUE))
  sort(unique(str_trim(unlist(str_split(x, ",")))))
}

short_species <- function(x) {
  case_when(
    is.na(x) | x == "" ~ "NA",
    TRUE ~ str_replace(x, "^Borrelia\\s+", "B. ")
  )
}

italic_md <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\b(B\\.\\s+\\p{L}+)", "<i>\\1</i>")
  x <- str_replace_all(x, "\\b(Borrelia\\s+\\p{L}+)", "<i>\\1</i>")
  x
}

first_non_empty <- function(x, default = "Inconnu") {
  x <- x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) default else x[1]
}

parse_mlst_tips <- function(labels) {
  tibble(raw_label = labels) %>%
    mutate(label = sub("_.*$", "", raw_label), base_label = sub("-\\d+$", "", label))
}

load_mlst_species <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(base_label = as.character(`n°`), species_vec = map(species_final_retained, extract_species_vec)) %>%
    distinct(base_label, .keep_all = TRUE) %>%
    pmap_dfr(function(base_label, species_vec) {
      if (length(species_vec) <= 1) {
        tibble(label = base_label, base_label = base_label, species_display = ifelse(length(species_vec) == 0, "NA", species_vec[1]))
      } else {
        tibble(label = paste0(base_label, "-", seq_along(species_vec)), base_label = base_label, species_display = species_vec)
      }
    }) %>%
    mutate(species_short = short_species(species_display))
}

load_tick_meta <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      base_label = as.character(`n°`),
      site_code = str_extract(`Type échantillon`, "S\\d+"),
      altitude = recode(site_code, S12 = "400 m", S13 = "800 m", S14 = "1200 m", S15 = "1500 m", .default = "Inconnu"),
      annee = coalesce(str_extract(as.character(`date collecte`), "(19|20)\\d{2}"), "Inconnue")
    ) %>%
    group_by(base_label) %>%
    summarise(altitude = first_non_empty(altitude), annee = first_non_empty(annee, "Inconnue"), .groups = "drop")
}

tr <- read.tree(in_tree)
tips0 <- parse_mlst_tips(tr$tip.label)
tips_to_drop <- tips0 %>% filter(label %in% samples_to_exclude | base_label %in% samples_to_exclude) %>% pull(raw_label)
if (length(tips_to_drop) > 0) tr <- drop.tip(tr, tips_to_drop)

tips <- parse_mlst_tips(tr$tip.label)
tr$tip.label <- tips$label
tr <- tryCatch(phangorn::midpoint(ladderize(tr, right = TRUE)), error = function(e) ladderize(tr, right = TRUE))

meta <- tips %>%
  select(label, base_label) %>%
  left_join(load_mlst_species(species_file), by = c("label", "base_label")) %>%
  left_join(load_tick_meta(tick_metadata_file), by = "base_label") %>%
  mutate(
    species_short = coalesce(species_short, "NA"),
    altitude = coalesce(altitude, "Inconnu"),
    annee = coalesce(annee, "Inconnue"),
    tip_label = label
  )

alt_matrix <- meta %>%
  select(label, Altitude = altitude) %>%
  column_to_rownames("label") %>%
  as.data.frame()

year_matrix <- meta %>%
  select(label, Année = annee) %>%
  column_to_rownames("label") %>%
  as.data.frame()

color_values <- target_colors[names(target_colors) %in% unique(meta$species_short)]
missing_species <- setdiff(unique(meta$species_short), names(color_values))
if (length(missing_species) > 0) color_values <- c(color_values, setNames(rep("#666666", length(missing_species)), missing_species))

p <- ggtree(tr, linewidth = 0.6) %<+% meta + theme_tree2() +
  geom_tippoint(aes(color = species_short), size = 3, alpha = 0.9, na.rm = TRUE) +
  geom_tiplab(aes(label = tip_label), align = TRUE, linetype = "dotted", linewidth = 0.25, size = 2.4, offset = 0.002) +
  scale_color_manual(values = color_values, labels = italic_md(names(color_values)), name = "Espèce") +
  guides(color = guide_legend(order = 1, ncol = 2, title.position = "top", override.aes = list(size = 4)))

p <- gheatmap(
  p,
  alt_matrix,
  offset = 0.08,
  width = 0.035,
  colnames = TRUE,
  colnames_position = "top",
  colnames_angle = 45,
  colnames_offset_y = 0,
  font.size = 2.8,
  color = "white",
  hjust = 0
) +
  scale_fill_manual(
    values = c(
      "400 m" = "#FEE8C8",
      "800 m" = "#FDBB84",
      "1200 m" = "#E34A33",
      "1500 m" = "#7F0000",
      "Inconnu" = "#D9D9D9"
    ),
    breaks = c("400 m", "800 m", "1200 m", "1500 m", "Inconnu"),
    na.value = "transparent",
    name = "Altitude"
  ) +
  guides(
    fill = guide_legend(
      order = 2,
      nrow = 1,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0,
      label.position = "right"
    )
  )

p <- p + ggnewscale::new_scale_fill()

p <- gheatmap(
  p,
  year_matrix,
  offset = 0.135,
  width = 0.035,
  colnames = TRUE,
  colnames_position = "top",
  colnames_angle = 45,
  colnames_offset_y = 0,
  font.size = 2.8,
  color = "white",
  hjust = 0
) +
  scale_fill_manual(
    values = c(
      "2023" = "#80B1D3",
      "2024" = "#FB8072",
      "2025" = "#B3DE69",
      "Inconnue" = "#BDBDBD"
    ),
    na.value = "transparent",
    name = "Année"
  ) +
  guides(
    fill = guide_legend(
      order = 3,
      nrow = 1,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0,
      label.position = "right"
    )
  ) +
  theme(
    plot.margin = margin(15, 30, 30, 30),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.box = "vertical",
    legend.box.just = "left",
    legend.direction = "horizontal",
    legend.title.position = "top",
    legend.title.align = 0,
    legend.key = element_blank(),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(fill = "white", color = NA),
    legend.margin = margin(4, 6, 4, 6),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key.width = unit(5, "mm"),
    legend.spacing.y = unit(2.5, "mm"),
    legend.spacing.x = unit(2.5, "mm"),
    legend.text = ggtext::element_markdown(size = 8.5),
    legend.title = element_text(size = 9.5, hjust = 0)
  )

print(p)
ggsave(out_png, p, width = 18, height = 10, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(out_pdf, p, width = 18, height = 10, bg = "white", limitsize = FALSE)
cat("Fichiers sauvegardés :\n", out_png, "\n", out_pdf, "\n")
