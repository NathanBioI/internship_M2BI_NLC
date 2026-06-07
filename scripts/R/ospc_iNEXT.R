# ==============================================================================
# ospC iNEXT - analyse complète simplifiée + figures propres
# ------------------------------------------------------------------------------
# Objectif :
#   1. Repartir des sorties finales du pipeline ospC.
#   2. Regrouper les séquences ospC en groupes à ~95 % d'identité.
#   3. Lancer iNEXT sur la richesse en groupes ospC.
#   4. Produire uniquement les figures principales et une table de probabilités.
#
# Entrées attendues, avec détection automatique ancien/nouveau chemin :
#   results/ospc/ospc_selection/selected_ospC_consensus_oriented.mafft.fasta
#   results/ospc/ospc_selection/ospc_sample_decisions.tsv
#   results/metadata_genaspe_species_final.tsv
#
# Sorties principales :
#   results/ospc_iNEXT_final/figures/*.pdf
#   results/ospc_iNEXT_final/figures/*.png
#   results/ospc_iNEXT_final/ospc_iNEXT_probabilites_couverture.tsv
# ===============================================================================

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

required_packages <- c("tidyverse", "readr", "ape", "iNEXT", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(
    "Packages manquants : ", paste(missing_packages, collapse = ", "),
    "\nInstalle-les avec : install.packages(c(",
    paste0("'", missing_packages, "'", collapse = ", "), "))"
  )
}

library(tidyverse)
library(readr)
library(ape)
library(iNEXT)
library(ggplot2)

# ------------------------------------------------------------------------------
# 2. Paramètres
# ------------------------------------------------------------------------------

selection_dirs <- c(
  "results/ospc/ospc_selection",
  "results/ospc_selection"
)

selection_dir <- selection_dirs[file.exists(selection_dirs)][1]
if (is.na(selection_dir)) {
  stop("Aucun dossier de sélection ospC trouvé dans : ", paste(selection_dirs, collapse = ", "))
}

alignment_file <- file.path(selection_dir, "selected_ospC_consensus_oriented.mafft.fasta")
decisions_file <- file.path(selection_dir, "ospc_sample_decisions.tsv")
metadata_file <- "results/metadata_genaspe_species_final.tsv"

outdir <- "results/ospc_iNEXT_final"
figdir <- file.path(outdir, "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# Distance 0.05 avec ape::dist.dna(model = "raw") ≈ groupes à 95 % d'identité.
cluster_max_distance <- 0.05
cluster_threshold_label <- "95_id"

# Figures totales : espèces gardées si nombre d'échantillons suffisant.
min_samples_per_species_total <- 5

# Figures site / année : on se concentre sur les deux espèces les plus représentées.
species_focus <- c("Borrelia garinii", "Borrelia valaisiana")
sites_to_keep <- c("S12", "S13")
site_label_map <- c("S12" = "400 m (S12)", "S13" = "800 m (S13)")

min_samples_per_subcurve <- 5
min_clusters_per_subcurve <- 2
require_repeated_group_for_subcurve <- TRUE

q_values_inext <- c(0)
nboot_inext <- 200
set.seed(123)

species_order_total <- c(
  "Borrelia garinii",
  "Borrelia valaisiana",
  "Borrelia afzelii",
  "Borrelia burgdorferi sensu stricto"
)

species_order_focus <- c("Borrelia garinii", "Borrelia valaisiana")

species_expr_map <- c(
  "Borrelia garinii" = "italic(Borrelia~garinii)",
  "Borrelia valaisiana" = "italic(Borrelia~valaisiana)",
  "Borrelia afzelii" = "italic(Borrelia~afzelii)",
  "Borrelia burgdorferi sensu stricto" = "italic(Borrelia~burgdorferi)~'sensu stricto'"
)

site_curve_order <- c("Total", "400 m (S12)", "800 m (S13)")
site_colors <- c("Total" = "black", "400 m (S12)" = "#D55E00", "800 m (S13)" = "#0072B2")
site_fills <- c("Total" = "grey55", "400 m (S12)" = "#D55E00", "800 m (S13)" = "#0072B2")

# ------------------------------------------------------------------------------
# 3. Fonctions utilitaires
# ------------------------------------------------------------------------------

sanitize_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

ospc_expr <- function(prefix = "", suffix = "") {
  parse(text = paste0("\"", prefix, "\"*italic(osp)*plain(C)*\"", suffix, "\""))[[1]]
}

clean_species <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "_", " ")
  x <- str_squish(x)
  x[x %in% c("", "NA", "NaN", "nan", "NULL", "null")] <- NA_character_
  x
}

is_excluded_species <- function(x) {
  x <- clean_species(x)
  x_low <- str_to_lower(x)
  is.na(x) |
    x == "" |
    str_detect(x_low, "co-infection|coinfection|co infection") |
    str_detect(x_low, "inconnu|unknown|ambiguous|ambigu") |
    str_detect(x_low, "miyamotoi")
}

find_first_existing_col <- function(df, candidates, context = "table") {
  found <- candidates[candidates %in% colnames(df)]
  if (length(found) == 0) {
    stop(
      "Aucune colonne parmi ", paste(candidates, collapse = ", "),
      " trouvée dans ", context, ". Colonnes disponibles : ",
      paste(colnames(df), collapse = ", ")
    )
  }
  found[1]
}

extract_collection_year <- function(date_collecte, id_collecte = NA_character_) {
  date_collecte <- as.character(date_collecte)
  id_collecte <- as.character(id_collecte)
  year_from_date <- str_extract(date_collecte, "(19|20)\\d{2}")
  year_from_id <- str_match(id_collecte, "_((19|20)\\d{2})\\d{4}$")[, 2]
  year <- ifelse(!is.na(year_from_date), year_from_date, year_from_id)
  ifelse(is.na(year) | year == "", NA_character_, year)
}

make_incidence_freq_vector <- function(data_curve) {
  n_units <- n_distinct(data_curve$sample)

  incidence <- data_curve %>%
    group_by(ospC_cluster) %>%
    summarise(freq = n_distinct(sample), .groups = "drop") %>%
    arrange(desc(freq), ospC_cluster)

  vec <- c(n_units, incidence$freq)
  names(vec) <- c("n_sampling_units", incidence$ospC_cluster)
  vec
}

map_fasta_headers_to_samples <- function(fasta_headers, decisions) {
  sample_list <- unique(decisions$sample)

  out <- tibble(
    fasta_header = fasta_headers,
    sample = NA_character_,
    mapping_method = NA_character_
  )

  if ("selected_candidate_id" %in% colnames(decisions)) {
    m <- decisions %>%
      filter(!is.na(selected_candidate_id), selected_candidate_id != "") %>%
      select(sample, selected_candidate_id) %>%
      distinct()

    idx <- match(out$fasta_header, m$selected_candidate_id)
    hit <- !is.na(idx) & is.na(out$sample)
    out$sample[hit] <- m$sample[idx[hit]]
    out$mapping_method[hit] <- "exact_selected_candidate_id"
  }

  if ("selected_template" %in% colnames(decisions)) {
    m <- decisions %>%
      filter(!is.na(selected_template), selected_template != "") %>%
      select(sample, selected_template) %>%
      distinct()

    idx <- match(out$fasta_header, m$selected_template)
    hit <- !is.na(idx) & is.na(out$sample)
    out$sample[hit] <- m$sample[idx[hit]]
    out$mapping_method[hit] <- "exact_selected_template"
  }

  for (i in which(is.na(out$sample))) {
    parts <- str_split(out$fasta_header[i], fixed("|"))[[1]] %>% str_squish()
    hit <- intersect(parts, sample_list)
    if (length(hit) >= 1) {
      out$sample[i] <- hit[1]
      out$mapping_method[i] <- "pipe_field_match"
    }
  }

  for (i in which(is.na(out$sample))) {
    candidates <- c(
      sub("__.*$", "", out$fasta_header[i]),
      sub("\\|.*$", "", out$fasta_header[i])
    )
    hit <- intersect(candidates, sample_list)
    if (length(hit) >= 1) {
      out$sample[i] <- hit[1]
      out$mapping_method[i] <- "prefix_match"
    }
  }

  for (i in which(is.na(out$sample))) {
    possible <- sample_list[vapply(sample_list, function(s) startsWith(out$fasta_header[i], s), logical(1))]
    if (length(possible) >= 1) {
      possible <- possible[order(nchar(possible), decreasing = TRUE)]
      out$sample[i] <- possible[1]
      out$mapping_method[i] <- "longest_sample_prefix"
    }
  }

  out
}

parse_curve_id <- function(x) {
  x <- as.character(x)
  mat <- str_split_fixed(x, "\\|\\|", 2)
  tibble(
    curve_id = x,
    species = ifelse(mat[, 2] == "", x, mat[, 1]),
    curve_group = ifelse(mat[, 2] == "", "Total", mat[, 2])
  )
}

recode_inext_method <- function(x) {
  x_low <- str_to_lower(as.character(x))
  case_when(
    str_detect(x_low, "observ") ~ "Observé",
    str_detect(x_low, "extrap") ~ "Extrapolation",
    str_detect(x_low, "interp|raref") ~ "Raréfaction",
    TRUE ~ as.character(x)
  )
}

extract_inext_plot_table <- function(inext_object, mode = c("size_based", "coverage_based")) {
  mode <- match.arg(mode)
  if (mode == "size_based") {
    tbl <- as_tibble(inext_object$iNextEst$size_based)
    x_col <- "t"
  } else {
    tbl <- as_tibble(inext_object$iNextEst$coverage_based)
    x_col <- "SC"
  }

  parsed <- parse_curve_id(tbl$Assemblage)

  tbl %>%
    bind_cols(parsed %>% select(species, curve_group, curve_id)) %>%
    mutate(
      q = .data[["Order.q"]],
      x_value = .data[[x_col]],
      Method_label = recode_inext_method(Method),
      qD.LCL = pmax(qD.LCL, 0),
      species = clean_species(species),
      curve_group = as.character(curve_group)
    )
}

extract_inext_data_info_columns <- function(data_info) {
  list(
    assemblage_col = find_first_existing_col(data_info, c("Assemblage", "site", "Site", "assemblage"), "iNEXT DataInfo"),
    n_units_col = find_first_existing_col(data_info, c("n", "T", "t", "m", "U", "sample_size", "SampleSize"), "iNEXT DataInfo"),
    s_obs_col = find_first_existing_col(data_info, c("S.obs", "S.obs.", "S_obs", "Sobs", "Observed"), "iNEXT DataInfo"),
    coverage_col = find_first_existing_col(data_info, c("SC", "SC(n)", "SampleCoverage", "sample_coverage"), "iNEXT DataInfo")
  )
}

make_coverage_summary <- function(inext_object, analysis) {
  data_info <- as_tibble(inext_object$DataInfo)
  cols <- extract_inext_data_info_columns(data_info)
  parsed <- parse_curve_id(data_info[[cols$assemblage_col]])

  data_info %>%
    bind_cols(parsed %>% select(species, curve_group, curve_id)) %>%
    transmute(
      analysis = analysis,
      curve_id = curve_id,
      species = clean_species(species),
      curve_group = curve_group,
      n_sampling_units = .data[[cols$n_units_col]],
      n_groups_observed = .data[[cols$s_obs_col]],
      sample_coverage = .data[[cols$coverage_col]],
      probability_new_sample_known_group = sample_coverage,
      probability_new_sample_new_group = 1 - sample_coverage
    ) %>%
    arrange(analysis, species, curve_group)
}

make_inext_input_from_table <- function(df, name_col = "curve_id") {
  df %>%
    group_split(.data[[name_col]]) %>%
    set_names(map_chr(., ~ unique(.x[[name_col]]))) %>%
    map(make_incidence_freq_vector)
}

run_inext <- function(df, analysis_name) {
  message("Lancement iNEXT : ", analysis_name)
  iNEXT(
    make_inext_input_from_table(df, "curve_id"),
    q = q_values_inext,
    datatype = "incidence_freq",
    se = TRUE,
    conf = 0.95,
    nboot = nboot_inext
  )
}

# Résumé des courbes site/année avant iNEXT, pour éviter les courbes trop faibles.
filter_subcurves <- function(curve_raw, min_samples, min_clusters, require_repeated_group) {
  summary <- curve_raw %>%
    distinct(curve_id, species_final_retained, curve_group, sample, ospC_cluster) %>%
    count(curve_id, species_final_retained, curve_group, ospC_cluster, name = "incidence_frequency") %>%
    group_by(curve_id, species_final_retained, curve_group) %>%
    summarise(
      n_samples = n_distinct(curve_raw$sample[curve_raw$curve_id == first(curve_id)]),
      n_ospC_clusters = n_distinct(ospC_cluster),
      n_repeated_clusters = sum(incidence_frequency >= 2),
      .groups = "drop"
    ) %>%
    mutate(
      min_required = if_else(curve_group == "Total", min_samples_per_species_total, min_samples),
      retained = n_samples >= min_required &
        n_ospC_clusters >= min_clusters &
        (curve_group == "Total" | !require_repeated_group | n_repeated_clusters >= 1)
    )

  curve_raw %>% filter(curve_id %in% summary$curve_id[summary$retained])
}

add_species_lab <- function(df, species_order) {
  df %>%
    mutate(
      species = factor(species, levels = species_order),
      species_lab = species_expr_map[as.character(species)],
      species_lab = factor(species_lab, levels = species_expr_map[species_order])
    ) %>%
    filter(!is.na(species), !is.na(species_lab))
}

# ------------------------------------------------------------------------------
# 4. Fonctions graphiques, avec labels propres de script_2.R
# ------------------------------------------------------------------------------

save_plot <- function(plot, filename, width = 11, height = 7) {
  ggsave(file.path(figdir, paste0(filename, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(figdir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 300)
}

plot_total_faceted <- function(plot_df, x_label, title, subtitle, output_prefix) {
  plot_df <- plot_df %>% filter(q == 0) %>% add_species_lab(species_order_total)

  p <- ggplot(plot_df, aes(x = x_value, y = qD)) +
    geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL, group = curve_id), fill = "grey75", alpha = 0.45, colour = NA) +
    geom_line(aes(linetype = Method_label, group = interaction(curve_id, Method_label)), linewidth = 0.9, colour = "black") +
    geom_point(data = plot_df %>% filter(Method_label == "Observé"), size = 2.4, colour = "black", show.legend = FALSE) +
    facet_wrap(~ species_lab, scales = "free_x", ncol = 2, labeller = label_parsed) +
    scale_linetype_manual(
      values = c("Raréfaction" = "solid", "Observé" = "solid", "Extrapolation" = "dashed"),
      breaks = c("Raréfaction", "Extrapolation"),
      labels = c("Raréfaction", "Extrapolation"),
      drop = TRUE
    ) +
    labs(title = title, subtitle = subtitle, x = x_label, y = ospc_expr("Richesse estimée en groupes ", ""), linetype = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  save_plot(p, output_prefix, 11, 7)
  invisible(p)
}

plot_focus_faceted <- function(plot_df, x_label, title, subtitle, output_prefix, curve_order, manual_colors, manual_fills, width = 11, height = 6) {
  plot_df <- plot_df %>%
    filter(q == 0, species %in% species_order_focus) %>%
    add_species_lab(species_order_focus) %>%
    mutate(curve_group = factor(curve_group, levels = curve_order)) %>%
    filter(!is.na(curve_group))

  p <- ggplot(plot_df, aes(x = x_value, y = qD)) +
    geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL, group = curve_id, fill = curve_group), alpha = 0.16, colour = NA) +
    geom_line(aes(colour = curve_group, linetype = Method_label, group = interaction(curve_id, Method_label)), linewidth = 0.95) +
    geom_point(data = plot_df %>% filter(Method_label == "Observé"), aes(colour = curve_group), size = 2.6, show.legend = FALSE) +
    facet_wrap(~ species_lab, scales = "free_x", ncol = 2, labeller = label_parsed) +
    scale_colour_manual(values = manual_colors, drop = TRUE) +
    scale_fill_manual(values = manual_fills, drop = TRUE) +
    scale_linetype_manual(
      values = c("Raréfaction" = "solid", "Observé" = "solid", "Extrapolation" = "dashed"),
      breaks = c("Raréfaction", "Extrapolation"),
      labels = c("Raréfaction", "Extrapolation"),
      drop = TRUE
    ) +
    labs(title = title, subtitle = subtitle, x = x_label, y = ospc_expr("Richesse estimée en groupes ", ""), colour = "Courbe", fill = "Courbe", linetype = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  save_plot(p, output_prefix, width, height)
  invisible(p)
}

# ------------------------------------------------------------------------------
# 5. Lecture des données et reconstruction sample / espèce / site / année
# ------------------------------------------------------------------------------

message("Lecture de l'alignement : ", alignment_file)
aln <- ape::read.dna(alignment_file, format = "fasta")
fasta_headers <- rownames(aln)

message("Lecture des décisions ospC : ", decisions_file)
decisions <- read_tsv(decisions_file, show_col_types = FALSE) %>%
  mutate(sample = as.character(sample))

message("Lecture des métadonnées : ", metadata_file)
metadata_raw <- read_tsv(metadata_file, show_col_types = FALSE)

metadata_id_col <- find_first_existing_col(
  metadata_raw,
  c("sample", "Sample", "sample_id", "id", "ID", "n°", "n", "n_ADN", "n°ADN EPIA", "ID_ADN", "ID_COLLECTE"),
  "metadata_file"
)

date_col <- find_first_existing_col(
  metadata_raw,
  c("date collecte", "Date collecte", "date_collecte", "collection_date", "date"),
  "metadata_file"
)

id_collecte_col <- if ("ID_COLLECTE" %in% colnames(metadata_raw)) "ID_COLLECTE" else NA_character_
site_col <- if ("Type échantillon" %in% colnames(metadata_raw)) "Type échantillon" else id_collecte_col

if (!"species_final_retained" %in% colnames(metadata_raw)) {
  stop("La colonne species_final_retained est absente de ", metadata_file)
}

metadata <- metadata_raw %>%
  mutate(
    metadata_id_raw = as.character(.data[[metadata_id_col]]),
    sample_key = str_extract(metadata_id_raw, "\\d+"),
    site_raw = if (!is.na(site_col)) as.character(.data[[site_col]]) else NA_character_,
    id_collecte_raw = if (!is.na(id_collecte_col)) as.character(.data[[id_collecte_col]]) else NA_character_,
    site_code = str_extract(site_raw, "S12|S13|S14"),
    site_code = if_else(is.na(site_code), str_extract(id_collecte_raw, "S12|S13|S14"), site_code),
    site_label = case_when(
      site_code == "S12" ~ site_label_map[["S12"]],
      site_code == "S13" ~ site_label_map[["S13"]],
      TRUE ~ site_code
    ),
    collection_year = extract_collection_year(.data[[date_col]], id_collecte_raw),
    species_final_retained = clean_species(species_final_retained)
  ) %>%
  filter(!is.na(sample_key), sample_key != "") %>%
  distinct(sample_key, .keep_all = TRUE)

header_map <- map_fasta_headers_to_samples(fasta_headers, decisions)

unmapped <- header_map %>% filter(is.na(sample) | sample == "")
if (nrow(unmapped) > 0) {
  stop(nrow(unmapped), " headers FASTA n'ont pas pu être reliés à un sample.")
}

input_clean <- header_map %>%
  mutate(sample_key = str_extract(sample, "^\\d+")) %>%
  left_join(metadata, by = "sample_key") %>%
  filter(!is_excluded_species(species_final_retained)) %>%
  distinct(fasta_header, sample, species_final_retained, site_code, site_label, collection_year, .keep_all = TRUE)

if (nrow(input_clean) == 0) {
  stop("Aucune séquence ospC exploitable après filtrage des espèces.")
}

# ------------------------------------------------------------------------------
# 6. Distances ospC et clustering à 95 %
# ------------------------------------------------------------------------------

headers_keep <- input_clean$fasta_header
aln_keep <- aln[headers_keep, , drop = FALSE]

message("Nombre de séquences gardées : ", nrow(aln_keep))

dist_mat <- ape::dist.dna(aln_keep, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
if (any(is.na(dist_mat))) {
  stop("La matrice de distance ospC contient des NA. Vérifie l'alignement MAFFT.")
}

clusters_vec <- cutree(hclust(as.dist(dist_mat), method = "average"), h = cluster_max_distance)

cluster_table <- tibble(
  fasta_header = names(clusters_vec),
  ospC_cluster_number = as.integer(clusters_vec),
  ospC_cluster = paste0(cluster_threshold_label, "_cluster_", ospC_cluster_number)
) %>%
  left_join(input_clean, by = "fasta_header") %>%
  arrange(species_final_retained, site_code, collection_year, ospC_cluster_number, sample)

cluster_summary <- cluster_table %>%
  group_by(species_final_retained) %>%
  summarise(
    n_samples = n_distinct(sample),
    n_sequences = n_distinct(fasta_header),
    n_ospC_clusters = n_distinct(ospC_cluster),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples))

message("Résumé groupes ospC :")
print(cluster_summary)

# ------------------------------------------------------------------------------
# 7. iNEXT : total par espèce, sites S12/S13, années
# ------------------------------------------------------------------------------

coverage_outputs <- list()

species_keep_total <- cluster_summary %>%
  filter(n_samples >= min_samples_per_species_total) %>%
  pull(species_final_retained)

total_curves <- cluster_table %>%
  filter(species_final_retained %in% species_keep_total) %>%
  mutate(curve_group = "Total", curve_id = species_final_retained)

inext_total <- run_inext(total_curves, "total_species")
coverage_outputs[["total_species"]] <- make_coverage_summary(inext_total, "total_species")

total_size_plot <- extract_inext_plot_table(inext_total, "size_based")
total_coverage_plot <- extract_inext_plot_table(inext_total, "coverage_based")

site_total <- cluster_table %>%
  filter(species_final_retained %in% species_focus) %>%
  mutate(curve_group = "Total", curve_id = paste(species_final_retained, curve_group, sep = "||"))

site_sub <- cluster_table %>%
  filter(species_final_retained %in% species_focus, site_code %in% sites_to_keep) %>%
  mutate(curve_group = site_label, curve_id = paste(species_final_retained, curve_group, sep = "||"))

site_curves <- bind_rows(site_total, site_sub) %>%
  filter_subcurves(min_samples_per_subcurve, min_clusters_per_subcurve, require_repeated_group_for_subcurve)

inext_sites <- run_inext(site_curves, "total_S12_S13")
coverage_outputs[["total_S12_S13"]] <- make_coverage_summary(inext_sites, "total_S12_S13")

total_sites_size_plot <- extract_inext_plot_table(inext_sites, "size_based")
total_sites_coverage_plot <- extract_inext_plot_table(inext_sites, "coverage_based")

year_total <- cluster_table %>%
  filter(species_final_retained %in% species_focus) %>%
  mutate(curve_group = "Total", curve_id = paste(species_final_retained, curve_group, sep = "||"))

year_sub <- cluster_table %>%
  filter(species_final_retained %in% species_focus, !is.na(collection_year), collection_year != "") %>%
  mutate(curve_group = paste0("Année ", collection_year), curve_id = paste(species_final_retained, curve_group, sep = "||"))

year_curves <- bind_rows(year_total, year_sub) %>%
  filter_subcurves(min_samples_per_subcurve, min_clusters_per_subcurve, require_repeated_group_for_subcurve)

inext_years <- run_inext(year_curves, "total_years")
coverage_outputs[["total_years"]] <- make_coverage_summary(inext_years, "total_years")

total_years_size_plot <- extract_inext_plot_table(inext_years, "size_based")
total_years_coverage_plot <- extract_inext_plot_table(inext_years, "coverage_based")

probability_table <- bind_rows(coverage_outputs) %>%
  mutate(
    sample_coverage_pct = sample_coverage * 100,
    probability_new_sample_known_group_pct = probability_new_sample_known_group * 100,
    probability_new_sample_new_group_pct = probability_new_sample_new_group * 100
  )

write_tsv(probability_table, file.path(outdir, "ospc_iNEXT_probabilites_couverture.tsv"))

# ------------------------------------------------------------------------------
# 8. Figures principales, avec mise en forme propre
# ------------------------------------------------------------------------------

plot_total_faceted(
  plot_df = total_size_plot,
  x_label = "Nombre d'échantillons",
  title = ospc_expr("Raréfaction/extrapolation de la richesse en groupes ", ""),
  subtitle = ospc_expr("Analyse iNEXT ; groupes ", " à 95 % d'identité approximative ; ruban gris = IC 95 %"),
  output_prefix = "Figure11A_iNEXT_q0_total_sample_size_by_species"
)

plot_total_faceted(
  plot_df = total_coverage_plot,
  x_label = "Couverture d'échantillonnage",
  title = ospc_expr("Richesse en groupes ", " standardisée par couverture"),
  subtitle = ospc_expr("Analyse iNEXT ; groupes ", " à 95 % d'identité approximative ; ruban gris = IC 95 %"),
  output_prefix = "Figure11A_iNEXT_q0_total_coverage_based_by_species"
)

plot_focus_faceted(
  plot_df = total_sites_size_plot,
  x_label = "Nombre d'échantillons",
  title = ospc_expr("Raréfaction/extrapolation des groupes ", " par site"),
  subtitle = ospc_expr("B. garinii et B. valaisiana ; Total, S12/400 m et S13/800 m ; groupes ", " à 95 % ; rubans = IC 95 % iNEXT"),
  output_prefix = "Figure11B_iNEXT_q0_total_S12_S13_sample_size_garinii_valaisiana",
  curve_order = site_curve_order,
  manual_colors = site_colors,
  manual_fills = site_fills,
  width = 11,
  height = 6
)

plot_focus_faceted(
  plot_df = total_sites_coverage_plot,
  x_label = "Couverture d'échantillonnage",
  title = ospc_expr("Richesse en groupes ", " par site standardisée par couverture"),
  subtitle = ospc_expr("B. garinii et B. valaisiana ; Total, S12/400 m et S13/800 m ; groupes ", " à 95 % ; rubans = IC 95 % iNEXT"),
  output_prefix = "Figure11B_iNEXT_q0_total_S12_S13_coverage_based_garinii_valaisiana",
  curve_order = site_curve_order,
  manual_colors = site_colors,
  manual_fills = site_fills,
  width = 11,
  height = 6
)

year_groups_present <- total_years_size_plot %>% filter(q == 0) %>% distinct(curve_group) %>% pull(curve_group)
year_curve_order <- c("Total", sort(setdiff(year_groups_present, "Total")))
year_colors <- setNames(grDevices::hcl.colors(length(year_curve_order), palette = "Dark 3"), year_curve_order)
year_colors["Total"] <- "black"
year_fills <- year_colors
year_fills["Total"] <- "grey55"

plot_focus_faceted(
  plot_df = total_years_size_plot,
  x_label = "Nombre d'échantillons",
  title = ospc_expr("Raréfaction/extrapolation des groupes ", " par année de collecte"),
  subtitle = ospc_expr("B. garinii et B. valaisiana ; Total et années de collecte ; groupes ", " à 95 % ; rubans = IC 95 % iNEXT"),
  output_prefix = "Figure11C_iNEXT_q0_total_years_sample_size_garinii_valaisiana",
  curve_order = year_curve_order,
  manual_colors = year_colors,
  manual_fills = year_fills,
  width = 11,
  height = 6
)

plot_focus_faceted(
  plot_df = total_years_coverage_plot,
  x_label = "Couverture d'échantillonnage",
  title = ospc_expr("Richesse en groupes ", " par année standardisée par couverture"),
  subtitle = ospc_expr("B. garinii et B. valaisiana ; Total et années de collecte ; groupes ", " à 95 % ; rubans = IC 95 % iNEXT"),
  output_prefix = "Figure11C_iNEXT_q0_total_years_coverage_based_garinii_valaisiana",
  curve_order = year_curve_order,
  manual_colors = year_colors,
  manual_fills = year_fills,
  width = 11,
  height = 6
)

cat("\nAnalyse iNEXT ospC terminée.\n")
cat("Figures : ", figdir, "\n")
cat("Table de probabilités : ", file.path(outdir, "ospc_iNEXT_probabilites_couverture.tsv"), "\n")
