library(ape)
library(phangorn)
library(dplyr)
library(ggplot2)
library(ggtree)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(tibble)

# --- 1. Chargement de l'arbre ---
in_tree <- "results/phylogeny/MLST_iqtree.treefile"
tr <- read.tree(in_tree)

# Echantillons / tips à exclure, repris depuis l'autre script
samples_to_exclude <- c("ET78", "EM2069", "ET43", "ET76", "ET79", "ET110", "ET44")

tree_tips_initial <- tibble(raw_label = tr$tip.label) %>%
  mutate(
    label_no_suffix = sub("_.*$", "", raw_label),
    sample_id = sub("-\\d+$", "", label_no_suffix),
    split_idx = if_else(
      str_detect(label_no_suffix, "-\\d+$"),
      as.integer(str_extract(label_no_suffix, "(?<=-)\\d+$")),
      NA_integer_
    )
  )

# Suppression des tips correspondant aux échantillons exclus,
# que le nom soit sous forme ETxx ou ETxx-1 / ETxx-2.
tips_to_drop <- tree_tips_initial %>%
  filter(label_no_suffix %in% samples_to_exclude | sample_id %in% samples_to_exclude) %>%
  pull(raw_label)

if (length(tips_to_drop) > 0) {
  tr <- drop.tip(tr, tips_to_drop)
}

# Mise à jour des labels après exclusion
tree_tips <- tibble(raw_label = tr$tip.label) %>%
  mutate(
    label_no_suffix = sub("_.*$", "", raw_label),
    sample_id = sub("-\\d+$", "", label_no_suffix),
    split_idx = if_else(
      str_detect(label_no_suffix, "-\\d+$"),
      as.integer(str_extract(label_no_suffix, "(?<=-)\\d+$")),
      NA_integer_
    )
  )

tr$tip.label <- tree_tips$label_no_suffix
tr <- ladderize(tr, right = TRUE)
tr <- phangorn::midpoint(tr)

# --- 2. Métadonnées espèce ---
meta_raw <- readr::read_tsv(
  "results/metadata_genaspe_species_final.tsv",
  show_col_types = FALSE,
  name_repair = "unique"
)

clean_species <- function(x) {
  x %>%
    replace_na("") %>%
    str_replace_all("\\s*\\(candidat\\)", "") %>%
    str_squish()
}

extract_species_vec <- function(x) {
  x <- x %>% replace_na("") %>% str_squish()
  if (x == "" || x == "NA" || x == "inconnu") return(character(0))
  
  x <- str_remove(x, regex("^co-infection\\s*-\\s*", ignore_case = TRUE))
  x <- clean_species(x)
  
  parts <- str_split(x, ",", simplify = FALSE)[[1]] %>%
    str_trim() %>%
    discard(~ .x == "") %>%
    unique() %>%
    sort()
  
  parts
}

meta_base <- meta_raw %>%
  mutate(
    base_label = as.character(`n°`),
    species_display_raw = case_when(
      !is.na(species_final_retained) &
        species_final_retained != "" &
        species_final_retained != "NA" ~ species_final_retained,
      TRUE ~ "NA"
    ),
    species_display_raw = clean_species(species_display_raw),
    species_vec = map(species_display_raw, extract_species_vec)
  ) %>%
  select(base_label, species_display_raw, species_vec) %>%
  distinct(base_label, .keep_all = TRUE)

meta_expanded <- purrr::map2_dfr(
  meta_base$base_label,
  meta_base$species_vec,
  \(base_label, species_vec) {
    if (length(species_vec) == 0) {
      tibble(
        label = base_label,
        base_label = base_label,
        species_display = "NA"
      )
    } else if (length(species_vec) == 1) {
      tibble(
        label = base_label,
        base_label = base_label,
        species_display = species_vec[1]
      )
    } else {
      tibble(
        label = paste0(base_label, "-", seq_along(species_vec)),
        base_label = base_label,
        species_display = species_vec
      )
    }
  }
) %>%
  mutate(
    species_short = str_replace(species_display, "^([A-Za-z])[a-z]+\\s+", "\\1. ")
  )

# On retire aussi les échantillons exclus des métadonnées expansées,
# afin qu'ils ne ressortent pas dans les diagnostics comme "absents de l'arbre".
meta_expanded <- meta_expanded %>%
  filter(!(label %in% samples_to_exclude | base_label %in% samples_to_exclude))

# --- 3. Fusion des métadonnées ---
meta <- tree_tips %>%
  transmute(
    label = label_no_suffix,
    base_label = sample_id
  ) %>%
  left_join(meta_expanded, by = c("label", "base_label")) %>%
  mutate(
    species_display = case_when(
      !is.na(species_display) & species_display != "" ~ species_display,
      TRUE ~ "NA"
    ),
    species_short = case_when(
      !is.na(species_short) & species_short != "" ~ species_short,
      TRUE ~ "NA"
    ),
    tip_label = label
  ) %>%
  distinct(label, .keep_all = TRUE)

# --- 4. Présence / absence des marqueurs depuis mafft_input ---
mafft_dir <- "results/pubMLST_typing/mafft_input"

fasta_files <- list.files(
  mafft_dir,
  pattern = "\\.fasta$",
  full.names = TRUE
)

fasta_files <- fasta_files[!grepl("\\.aln\\.fasta$", fasta_files)]

target_loci <- c("clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA")
fasta_files <- fasta_files[basename(fasta_files) %in% paste0(target_loci, ".fasta")]

read_fasta_headers <- function(path) {
  x <- readLines(path, warn = FALSE)
  x[grepl("^>", x)] %>%
    sub("^>", "", .) %>%
    str_trim()
}

marker_long <- purrr::map_dfr(fasta_files, function(f) {
  locus <- sub("\\.fasta$", "", basename(f))
  headers <- read_fasta_headers(f)
  
  tibble(
    locus = locus,
    header = headers
  ) %>%
    mutate(
      label = sub("\\|.*$", "", header),
      base_label = sub("-\\d+$", "", label)
    ) %>%
    filter(!(label %in% samples_to_exclude | base_label %in% samples_to_exclude)) %>%
    count(label, locus, name = "value")
})

# Diagnostic conservé : indique uniquement s'il y avait plusieurs séquences
# pour un même locus, mais la heatmap reste en présence / absence.
marker_duplicates <- marker_long %>%
  filter(value > 1)

cat("\n--- Duplications de marqueurs avant binarisation (value > 1) ---\n")
print(marker_duplicates)

marker_matrix <- tibble(label = tr$tip.label) %>%
  left_join(
    marker_long %>%
      tidyr::pivot_wider(
        names_from = locus,
        values_from = value,
        values_fill = 0
      ),
    by = "label"
  )

for (loc in target_loci) {
  if (!loc %in% colnames(marker_matrix)) {
    marker_matrix[[loc]] <- 0L
  }
}

marker_matrix <- marker_matrix %>%
  select(label, all_of(target_loci)) %>%
  distinct(label, .keep_all = TRUE)

marker_df <- marker_matrix %>%
  column_to_rownames("label") %>%
  as.data.frame()

# Heatmap binaire : présent = case grise ; absent = case vide.
marker_df_presence <- marker_df
marker_df_presence[] <- lapply(
  marker_df_presence,
  function(x) ifelse(x > 0, "Présent", NA_character_)
)

# --- 5. Diagnostics ---
cat("\n--- Tips exclus de l'arbre ---\n")
print(tips_to_drop)

cat("\n--- Entrées métadonnées absentes de l'arbre ---\n")
print(setdiff(unique(meta_expanded$label), tr$tip.label))

cat("\n--- Tips de l'arbre sans métadonnées espèce associées ---\n")
print(setdiff(tr$tip.label, unique(meta$label)))

cat("\n--- Tips de l'arbre sans info marqueur ---\n")
print(setdiff(tr$tip.label, rownames(marker_df_presence)))

# --- 6. Graphique arbre ---
p <- ggtree(tr, linewidth = 0.6) %<+% meta + theme_tree2()

p <- p +
  geom_tippoint(
    aes(color = species_short),
    shape = 16,       # point rond partout
    size = 3,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  geom_tiplab(
    aes(label = tip_label),
    align = TRUE,
    linetype = "dotted",
    linewidth = 0.25,
    size = 2.8,
    offset = 0.002,
    na.rm = TRUE
  ) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(20, 10, 10, 10)
  ) +
  guides(
    color = guide_legend(
      title = "Espèce",
      order = 1,
      label.theme = element_text(face = "italic")
    )
  )

distinct_colors <- c(
  "B. afzelii" = "#E41A1C",
  "B. burgdorferi sensu stricto" = "#f5a4e1",
  "B. garinii" = "#228B22",
  "B. miyamotoi" = "#FF8C00",
  "B. valaisiana" = "#8B4513",
  "B. lusitaniae" = "#b103fc",
  "B. spielmanii" = "#00A6A6",
  "B. bavariensis" = "#7A4EAB",
  "NA" = "#999999"
)

present_species <- unique(meta$species_short)
color_values <- distinct_colors[names(distinct_colors) %in% present_species]

missing_species <- setdiff(present_species, names(color_values))
if (length(missing_species) > 0) {
  extra_cols <- setNames(rep("#666666", length(missing_species)), missing_species)
  color_values <- c(color_values, extra_cols)
}

p <- p + scale_color_manual(values = color_values)

# --- 7. Ajout des 8 colonnes de marqueurs à droite ---
p <- gheatmap(
  p,
  marker_df_presence,
  offset = 0.08,
  width = 0.35,
  colnames = TRUE,
  colnames_position = "top",
  colnames_angle = 90,
  colnames_offset_y = 0,
  font.size = 2.8,
  hjust = 0
) +
  scale_fill_manual(
    values = c("Présent" = "#BDBDBD"),
    na.value = "white",
    guide = "none"
  ) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.box = "horizontal",
    legend.key = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

print(p)

# Sauvegardes possibles
# ggsave("tree_MLST_heatmap_presence_absence.svg", plot = p, width = 12, height = 18, limitsize = FALSE)
ggsave("tree_MLST_heatmap_presence_absence.png", plot = p, width = 12, height = 18, dpi = 1200, limitsize = FALSE)
