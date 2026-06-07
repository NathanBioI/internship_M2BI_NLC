# ==============================================================================
# SCRIPT COMPLET : ARBRE ospC + HEATMAP ALTITUDE + ANNÉE
# VERSION A4 PORTRAIT - HEATMAP COLLÉE À L'ARBRE
# - arbre orientation normale
# - suppression des labels alignés qui créaient le grand espace
# - heatmap compacte à 2 colonnes : Altitude + Année
# - heatmap collée au bout des branches
# - légende au-dessus, hors de l'arbre
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Chargement des bibliothèques
# ------------------------------------------------------------------------------
if (!require("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  ape, phangorn, tidyverse, ggtree, ggtext,
  readr, lubridate, grid, tibble
)

# ------------------------------------------------------------------------------
# 2. Paramètres
# ------------------------------------------------------------------------------

in_tree <- "results/ospc/ospc_selection/selected_ospC_iqtree.treefile"

ospc_decisions_file <- "results/ospc/ospc_selection/ospc_sample_decisions.tsv"
species_file <- "results/metadata_genaspe_species_final.tsv"
tick_metadata_file <- "data/metadata_genaspe.csv"

out_svg <- "ospC_tree_A4_heatmap_collee.svg"
out_pdf <- "ospC_tree_A4_heatmap_collee.pdf"

DISPLAY_FULL_SAMPLE_LABEL <- FALSE
samples_to_exclude <- c("ET77", "ET78")

# ------------------------------------------------------------------------------
# 3. Fonctions utilitaires
# ------------------------------------------------------------------------------

italicize_species_only <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "NA"] <- "NA"
  
  idx <- !is.na(x) & x != "NA"
  
  x[idx] <- stringr::str_replace_all(
    x[idx],
    "\\b([A-Z]\\.\\s+\\p{L}+)\\b",
    "<i>\\1</i>"
  )
  
  x[idx] <- stringr::str_replace_all(
    x[idx],
    "\\b(Borrelia\\s+\\p{L}+)\\b",
    "<i>\\1</i>"
  )
  
  x
}

clean_species <- function(x) {
  if (is.na(x)) return(NA_character_)
  
  x %>%
    as.character() %>%
    str_replace_all("\\s*\\(candidat\\)", "") %>%
    str_squish()
}

extract_species_vec <- function(x) {
  if (is.na(x) || x == "" || x == "NA" || x == "inconnu") {
    return(character(0))
  }
  
  x <- str_remove(x, regex("^co-infection\\s*-\\s*", ignore_case = TRUE))
  x <- clean_species(x)
  
  parts <- str_split(x, ",")[[1]] %>%
    str_trim()
  
  parts <- parts[parts != "" & !is.na(parts)]
  
  sort(unique(parts))
}

first_non_na <- function(x, default = NA_character_) {
  x <- x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) return(default)
  x[1]
}

parse_ospc_tip_labels <- function(tip_labels) {
  tibble(raw_label = tip_labels) %>%
    mutate(
      sample_full = str_extract(raw_label, "^[^|]+"),
      sample_full = if_else(is.na(sample_full) | sample_full == "", raw_label, sample_full),
      base_label = str_extract(sample_full, "^[0-9]+"),
      base_label = if_else(is.na(base_label) | base_label == "", sample_full, base_label),
      decision_from_label = str_match(raw_label, "^[^|]+\\|([^|]+)\\|")[, 2],
      candidate_from_label = str_match(raw_label, "^[^|]+\\|[^|]+\\|([^|]+)")[, 2]
    )
}

# ------------------------------------------------------------------------------
# 4. Chargement et nettoyage de l'arbre ospC
# ------------------------------------------------------------------------------

if (!file.exists(in_tree)) {
  stop("Fichier arbre introuvable : ", in_tree)
}

tr <- read.tree(in_tree)

tree_tips_initial <- parse_ospc_tip_labels(tr$tip.label)

tips_to_drop <- tree_tips_initial %>%
  filter(
    raw_label %in% samples_to_exclude |
      sample_full %in% samples_to_exclude |
      base_label %in% samples_to_exclude
  ) %>%
  pull(raw_label)

if (length(tips_to_drop) > 0) {
  tr <- drop.tip(tr, tips_to_drop)
}

tree_tips <- parse_ospc_tip_labels(tr$tip.label)

if (DISPLAY_FULL_SAMPLE_LABEL) {
  tr$tip.label <- tree_tips$sample_full
  tree_tips <- tree_tips %>%
    mutate(label = sample_full)
} else {
  tr$tip.label <- tree_tips$base_label
  tree_tips <- tree_tips %>%
    mutate(label = base_label)
}

tr <- ladderize(tr, right = TRUE)

tr <- tryCatch(
  phangorn::midpoint(tr),
  error = function(e) {
    message("Midpoint rooting impossible, arbre conservé tel quel : ", e$message)
    tr
  }
)

# ------------------------------------------------------------------------------
# 5. Chargement des décisions ospC
# ------------------------------------------------------------------------------

if (file.exists(ospc_decisions_file)) {
  ospc_decisions <- read_tsv(ospc_decisions_file, show_col_types = FALSE) %>%
    mutate(
      sample = as.character(sample),
      base_label = str_extract(sample, "^[0-9]+"),
      base_label = if_else(is.na(base_label) | base_label == "", sample, base_label),
      selected_source = if_else(is.na(selected_source) | selected_source == "", "NA", selected_source),
      decision = if_else(is.na(decision) | decision == "", "NA", decision),
      selected_len = as.character(selected_len),
      selected_depth = as.character(selected_depth),
      selected_blast_species = if_else(
        is.na(selected_blast_species) | selected_blast_species == "",
        "NA",
        selected_blast_species
      ),
      mlst_species = if_else(
        is.na(mlst_species) | mlst_species == "",
        "NA",
        mlst_species
      )
    ) %>%
    select(
      base_label, sample, mlst_species, mlst_status, decision,
      selected_candidate_id, selected_template, selected_source,
      selected_len, selected_depth, selected_orf_aa,
      selected_blast_hit, selected_blast_species, notes
    ) %>%
    distinct(base_label, .keep_all = TRUE)
} else {
  warning("Fichier de décisions ospC introuvable : ", ospc_decisions_file)
  
  ospc_decisions <- tibble(
    base_label = unique(tree_tips$base_label),
    sample = NA_character_,
    mlst_species = NA_character_,
    mlst_status = NA_character_,
    decision = NA_character_,
    selected_candidate_id = NA_character_,
    selected_template = NA_character_,
    selected_source = NA_character_,
    selected_len = NA_character_,
    selected_depth = NA_character_,
    selected_orf_aa = NA_character_,
    selected_blast_hit = NA_character_,
    selected_blast_species = NA_character_,
    notes = NA_character_
  )
}

# ------------------------------------------------------------------------------
# 6. Métadonnées espèces MLST
# ------------------------------------------------------------------------------

if (file.exists(species_file)) {
  meta_raw <- read_tsv(species_file, show_col_types = FALSE)
  
  meta_base <- meta_raw %>%
    mutate(
      base_label = as.character(`n°`),
      species_display_raw = case_when(
        !is.na(species_final_retained) & species_final_retained != "NA" ~ as.character(species_final_retained),
        TRUE ~ "NA"
      ),
      species_display_raw = map_chr(species_display_raw, clean_species),
      species_vec = map(species_display_raw, extract_species_vec)
    ) %>%
    select(base_label, species_vec) %>%
    distinct(base_label, .keep_all = TRUE)
  
  meta_expanded <- map2_dfr(meta_base$base_label, meta_base$species_vec, \(base_label, species_vec) {
    if (length(species_vec) <= 1) {
      tibble(
        base_label = base_label,
        species_display = if (length(species_vec) == 0) "NA" else species_vec[1]
      )
    } else {
      tibble(
        base_label = base_label,
        species_display = paste(species_vec, collapse = ", ")
      )
    }
  }) %>%
    mutate(
      species_short = case_when(
        species_display == "NA" ~ "NA",
        str_detect(species_display, ",") ~ "Co-infection",
        TRUE ~ str_replace(species_display, "^([A-Za-z])[a-z]+\\s+", "\\1. ")
      ),
      species_short_md = map_chr(species_short, italicize_species_only)
    )
} else {
  warning("Fichier espèces introuvable : ", species_file)
  
  meta_expanded <- tibble(
    base_label = unique(tree_tips$base_label),
    species_display = "NA",
    species_short = "NA",
    species_short_md = "NA"
  )
}

# ------------------------------------------------------------------------------
# 7. Métadonnées tiques : altitude + année
# ------------------------------------------------------------------------------

if (file.exists(tick_metadata_file)) {
  tick_raw <- read_tsv(tick_metadata_file, show_col_types = FALSE)
  
  tick_meta <- tick_raw %>%
    mutate(
      base_label = as.character(`n°`),
      site_code = str_extract(`Type échantillon`, "S\\d+"),
      altitude = case_when(
        site_code == "S12" ~ "400 m",
        site_code == "S13" ~ "800 m",
        site_code == "S14" ~ "1200 m",
        site_code == "S15" ~ "1500 m",
        TRUE ~ "Inconnu"
      ),
      date_full = dmy(`date collecte`),
      annee = case_when(
        !is.na(date_full) ~ as.character(year(date_full)),
        TRUE ~ "Inconnue"
      )
    ) %>%
    filter(!is.na(base_label), base_label != "") %>%
    group_by(base_label) %>%
    summarise(
      altitude = first_non_na(altitude, "Inconnu"),
      annee = first_non_na(annee, "Inconnue"),
      .groups = "drop"
    )
} else {
  warning("Fichier métadonnées tiques introuvable : ", tick_metadata_file)
  
  tick_meta <- tibble(
    base_label = unique(tree_tips$base_label),
    altitude = "Inconnu",
    annee = "Inconnue"
  )
}

# ------------------------------------------------------------------------------
# 8. Fusion finale des métadonnées
# ------------------------------------------------------------------------------

meta <- tree_tips %>%
  transmute(
    label = label,
    raw_label = raw_label,
    sample_full = sample_full,
    base_label = base_label,
    decision_from_label = decision_from_label,
    candidate_from_label = candidate_from_label
  ) %>%
  left_join(ospc_decisions, by = "base_label") %>%
  left_join(meta_expanded, by = "base_label") %>%
  left_join(tick_meta, by = "base_label") %>%
  mutate(
    species_display = case_when(
      !is.na(species_display) & species_display != "NA" ~ species_display,
      !is.na(mlst_species) & mlst_species != "NA" ~ mlst_species,
      TRUE ~ "NA"
    ),
    species_short = case_when(
      is.na(species_short) | species_short == "NA" ~ str_replace(species_display, "^([A-Za-z])[a-z]+\\s+", "\\1. "),
      TRUE ~ species_short
    ),
    species_short = if_else(is.na(species_short) | species_short == "", "NA", species_short),
    species_short_md = map_chr(species_short, italicize_species_only),
    
    altitude = if_else(is.na(altitude), "Inconnu", altitude),
    annee = if_else(is.na(annee), "Inconnue", annee),
    
    decision = if_else(is.na(decision), decision_from_label, decision),
    decision = if_else(is.na(decision), "NA", decision),
    
    selected_source = if_else(is.na(selected_source), "NA", selected_source),
    selected_len = if_else(is.na(selected_len), "NA", selected_len),
    selected_depth = if_else(is.na(selected_depth), "NA", selected_depth),
    
    tip_label = if (isTRUE(DISPLAY_FULL_SAMPLE_LABEL)) sample_full else base_label
  )

# ------------------------------------------------------------------------------
# 9. Heatmap compacte : Altitude + Année
# ------------------------------------------------------------------------------

metadata_matrix <- meta %>%
  select(label, altitude, annee) %>%
  rename(
    Altitude = altitude,
    Année = annee
  ) %>%
  column_to_rownames("label") %>%
  as.data.frame()

# ------------------------------------------------------------------------------
# 10. Diagnostics
# ------------------------------------------------------------------------------

cat("\n--- Fichier arbre ospC ---\n")
print(in_tree)

cat("\n--- Nombre de tips dans l'arbre ospC ---\n")
print(length(tr$tip.label))

cat("\n--- Tips exclus de l'arbre ---\n")
print(tips_to_drop)

cat("\n--- Résumé espèces MLST sur l'arbre ospC ---\n")
print(table(meta$species_short, useNA = "ifany"))

cat("\n--- Résumé altitude ---\n")
print(table(meta$altitude, useNA = "ifany"))

cat("\n--- Résumé année ---\n")
print(table(meta$annee, useNA = "ifany"))

cat("\n--- Résumé décisions ospC ---\n")
print(table(meta$decision, useNA = "ifany"))

cat("\n--- Résumé source des séquences sélectionnées ---\n")
print(table(meta$selected_source, useNA = "ifany"))

# ------------------------------------------------------------------------------
# 11. Arbre de base
# ------------------------------------------------------------------------------

p <- ggtree(tr, linewidth = 0.55) %<+% meta + theme_tree2()

p <- p +
  geom_tippoint(
    aes(color = species_short),
    shape = 16,
    size = 2.3,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  geom_tiplab(
    aes(label = tip_label),
    align = FALSE,
    size = 1.9,
    hjust = 0,
    offset = 0,
    na.rm = TRUE
  )

# ------------------------------------------------------------------------------
# 12. Couleurs espèces
# ------------------------------------------------------------------------------

distinct_colors <- c(
  "B. afzelii" = "#E41A1C",
  "B. burgdorferi sensu stricto" = "#f5a4e1",
  "B. burgdorferi sensu lato" = "#984ea3",
  "B. garinii" = "#228B22",
  "B. miyamotoi" = "#FF8C00",
  "B. valaisiana" = "#8B4513",
  "B. lusitaniae" = "#b103fc",
  "B. spielmanii" = "#00A6A6",
  "B. bavariensis" = "#7A4EAB",
  "Co-infection" = "#000000",
  "NA" = "#999999"
)

present_species <- unique(meta$species_short)
present_species <- present_species[!is.na(present_species)]

color_values <- distinct_colors[names(distinct_colors) %in% present_species]
missing_species <- setdiff(present_species, names(color_values))

if (length(missing_species) > 0) {
  extra_cols <- setNames(rep("#666666", length(missing_species)), missing_species)
  color_values <- c(color_values, extra_cols)
}

species_labels_md <- setNames(
  italicize_species_only(names(color_values)),
  names(color_values)
)

p <- p +
  scale_color_manual(
    values = color_values,
    labels = species_labels_md,
    name = "Espèce MLST"
  )

# ------------------------------------------------------------------------------
# 13. Heatmap collée à l'arbre
# ------------------------------------------------------------------------------

metadata_colors <- c(
  "400 m" = "#FEE8C8",
  "800 m" = "#FDBB84",
  "1200 m" = "#E34A33",
  "1500 m" = "#7F0000",
  "Inconnu" = "#D9D9D9",
  "2022" = "#BEBADA",
  "2023" = "#80B1D3",
  "2024" = "#FB8072",
  "2025" = "#B3DE69",
  "2026" = "#FDB462",
  "Inconnue" = "#BDBDBD"
)

p <- gheatmap(
  p,
  metadata_matrix,
  offset = 0.006,
  width = 0.050,
  colnames = TRUE,
  colnames_position = "top",
  colnames_angle = 45,
  colnames_offset_y = 0.06,
  font.size = 2.0,
  color = "white",
  hjust = 0
) +
  scale_fill_manual(
    values = metadata_colors,
    breaks = c(
      "400 m", "800 m", "1200 m", "1500 m",
      "2023", "2024", "2025", "2026"
    ),
    na.value = "transparent",
    name = "Altitude / Année"
  )
p <- p +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  coord_cartesian(clip = "off")

# ------------------------------------------------------------------------------
# 14. Thème
# ------------------------------------------------------------------------------

p <- p +
  guides(
    color = guide_legend(
      order = 1,
      nrow = 2,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0,
      label.position = "right",
      override.aes = list(size = 3.5)
    ),
    fill = guide_legend(
      order = 2,
      nrow = 1,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0,
      label.position = "right"
    )
  ) +
  theme(
    legend.position = c(0.08, 0.985),
    legend.justification = c(0, 1),
    legend.box = "vertical",
    legend.box.just = "left",
    legend.direction = "horizontal",
    
    legend.title.position = "top",
    legend.title.align = 0,
    
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 1, 0),
    legend.key.width = unit(3.8, "mm"),
    legend.key.height = unit(3.8, "mm"),
    legend.spacing.y = unit(1.2, "mm"),
    legend.spacing.x = unit(1.8, "mm"),
    legend.text = ggtext::element_markdown(size = 6.8),
    legend.title = element_text(size = 7.8, hjust = 0),
    
    plot.margin = margin(5, 5, 5, 5)
  )

print(p)

# ------------------------------------------------------------------------------
# 15. Sauvegarde A4 portrait
# ------------------------------------------------------------------------------

#ggsave(out_svg, plot = p, width = 8.27, height = 11.69, limitsize = FALSE)
ggsave(out_pdf, plot = p, width = 8.27, height = 11.69, limitsize = FALSE)
#ggsave(out_png, plot = p, width = 8.27, height = 11.69, units = "in", dpi = 1200, bg = "white", limitsize = FALSE)
#ggsave(out_jpg, plot = p, width = 8.27, height = 11.69, units = "in", dpi = 1200, bg = "white", limitsize = FALSE)
cat("\n--- Fichiers sauvegardés ---\n")
cat(out_svg, "\n")
cat(out_pdf, "\n")
