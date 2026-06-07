# ============================================================================
# Arbres MLST et ospC face à face + liens entre échantillons communs
# Entrées : results/phylogeny/MLST_iqtree.treefile
#           results/ospc/ospc_selection/selected_ospC_iqtree.treefile
#           results/metadata_genaspe_species_final.tsv
#           data/metadata_genaspe.csv
# Sorties : face_to_face_MLST_ospC.png/.jpg/.svg/.tsv
# ============================================================================

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(ape, phangorn, tidyverse, ggtree, ggtext, readr, tibble)

mlst_tree_file <- "results/phylogeny/MLST_iqtree.treefile"
ospc_tree_file <- "results/ospc/ospc_selection/selected_ospC_iqtree.treefile"
species_file <- "results/metadata_genaspe_species_final.tsv"
tick_metadata_file <- "data/metadata_genaspe.csv"

out_png <- "face_to_face_MLST_ospC.png"
out_jpg <- "face_to_face_MLST_ospC.jpg"
out_svg <- "face_to_face_MLST_ospC.svg"
out_links <- "face_to_face_MLST_ospC.tsv"

mlst_samples_to_exclude <- c("ET78", "EM2069", "ET43", "ET76", "ET79", "ET110", "ET44")
ospc_samples_to_exclude <- c("ET77", "ET78")
mlst_species_to_exclude <- c("B. miyamotoi")
ospc_species_to_exclude <- c("Co-infection")

TREE_WIDTH <- 3.2
INNER_GAP <- 2.2
LABEL_OFFSET <- 0.035
PLOT_WIDTH <- 11.69
PLOT_HEIGHT <- 8.27
PLOT_DPI <- 300

species_colors <- c(
  "B. afzelii" = "#E41A1C", "B. burgdorferi sensu stricto" = "#f5a4e1",
  "B. garinii" = "#228B22", "B. miyamotoi" = "#FF8C00",
  "B. valaisiana" = "#8B4513", "B. lusitaniae" = "#b103fc",
  "Co-infection" = "#000000", "NA" = "#999999"
)
site_shapes <- c("400 m" = 16, "800 m" = 17, "1200 m" = 15, "1500 m" = 18, "Inconnu" = 1)

clean_species <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  str_squish(str_replace_all(x, "\\s*\\(candidat\\)", ""))
}

extract_species_vec <- function(x) {
  x <- clean_species(x)
  if (x %in% c("", "NA", "inconnu")) return(character(0))
  x <- str_remove(x, regex("^co-infection\\s*-\\s*", ignore_case = TRUE))
  sort(unique(str_trim(unlist(str_split(x, ",")))))
}

short_species <- function(x) case_when(is.na(x) | x == "" ~ "NA", TRUE ~ str_replace(x, "^Borrelia\\s+", "B. "))

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

rescale_num <- function(x, to, from = range(x, na.rm = TRUE)) {
  if (!is.finite(diff(from)) || diff(from) == 0) return(rep(mean(to), length(x)))
  (x - from[1]) / diff(from) * diff(to) + to[1]
}

parse_mlst_tips <- function(labels) {
  tibble(raw_label = labels) %>%
    mutate(label = sub("_.*$", "", raw_label), base_label = sub("-\\d+$", "", label))
}

parse_ospc_tips <- function(labels) {
  tibble(raw_label = labels) %>%
    mutate(
      sample_full = str_extract(raw_label, "^[^|]+"),
      sample_full = if_else(is.na(sample_full) | sample_full == "", raw_label, sample_full),
      base_label = str_extract(sample_full, "^[0-9]+"),
      base_label = if_else(is.na(base_label) | base_label == "", sample_full, base_label),
      label = base_label
    )
}

prepare_tree <- function(tree_file, parser, exclude) {
  tr <- read.tree(tree_file)
  tips0 <- parser(tr$tip.label)
  drop <- tips0 %>% filter(raw_label %in% exclude | label %in% exclude | base_label %in% exclude) %>% pull(raw_label)
  if (length(drop) > 0) tr <- drop.tip(tr, drop)
  tips <- parser(tr$tip.label)
  tr$tip.label <- tips$label
  tr <- tryCatch(phangorn::midpoint(ladderize(tr, right = TRUE)), error = function(e) ladderize(tr, right = TRUE))
  list(tree = tr, tips = parser(tr$tip.label))
}

load_species_meta <- function(path) {
  base <- read_tsv(path, show_col_types = FALSE) %>%
    transmute(base_label = as.character(`n°`), species_vec = map(species_final_retained, extract_species_vec)) %>%
    distinct(base_label, .keep_all = TRUE)

  mlst <- pmap_dfr(base, function(base_label, species_vec) {
    if (length(species_vec) <= 1) {
      tibble(label = base_label, base_label = base_label, species = ifelse(length(species_vec) == 0, "NA", species_vec[1]))
    } else {
      tibble(label = paste0(base_label, "-", seq_along(species_vec)), base_label = base_label, species = species_vec)
    }
  }) %>% mutate(species_short = short_species(species))

  ospc <- base %>%
    mutate(species = map_chr(species_vec, ~ if (length(.x) == 0) "NA" else if (length(.x) == 1) .x[1] else "Co-infection"),
           species_short = short_species(species)) %>%
    select(base_label, species_short)

  list(mlst = mlst %>% select(label, base_label, species_short), ospc = ospc)
}

load_tick_meta <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      base_label = as.character(`n°`),
      site_code = str_extract(`Type échantillon`, "S\\d+"),
      site = recode(site_code, S12 = "400 m", S13 = "800 m", S14 = "1200 m", S15 = "1500 m", .default = "Inconnu")
    ) %>%
    group_by(base_label) %>%
    summarise(site = first_non_empty(site), .groups = "drop")
}

annotate_tips <- function(tips, species_meta, tick_meta, by_label = TRUE) {
  joined <- if (by_label) {
    tips %>% left_join(species_meta, by = c("label", "base_label"))
  } else {
    tips %>% left_join(species_meta, by = "base_label")
  }
  joined %>%
    left_join(tick_meta, by = "base_label") %>%
    mutate(species_short = coalesce(species_short, "NA"), site = coalesce(site, "Inconnu"))
}

filter_tree_species <- function(tr, tips, species_exclude) {
  to_drop <- tips %>% filter(species_short %in% species_exclude) %>% pull(label) %>% intersect(tr$tip.label)
  if (length(to_drop) > 0) {
    tr <- drop.tip(tr, to_drop)
    tr <- ladderize(tr, right = TRUE)
    tips <- tips %>% filter(label %in% tr$tip.label)
  }
  list(tree = tr, tips = tips)
}

layout_tree <- function(tr, tips, side = c("left", "right"), y_to = c(1, length(tr$tip.label))) {
  side <- match.arg(side)
  d <- ggtree(tr)$data %>% as_tibble()
  x_max <- max(d$x, na.rm = TRUE); if (!is.finite(x_max) || x_max == 0) x_max <- 1
  tip_col <- ifelse(side == "left", -INNER_GAP / 2, INNER_GAP / 2)
  root_col <- ifelse(side == "left", tip_col - TREE_WIDTH, tip_col + TREE_WIDTH)

  d <- d %>%
    mutate(
      x0 = x,
      y = rescale_num(y, y_to),
      x = if (side == "left") root_col + x0 / x_max * TREE_WIDTH else root_col - x0 / x_max * TREE_WIDTH
    ) %>%
    left_join(tips %>% distinct(label, .keep_all = TRUE), by = "label")

  edges <- d %>%
    filter(!is.na(parent), node != parent) %>%
    select(node, parent, x_child = x, y_child = y) %>%
    left_join(d %>% select(parent = node, x_parent = x, y_parent = y), by = "parent")

  branches <- bind_rows(
    edges %>% transmute(x = x_parent, xend = x_child, y = y_child, yend = y_child),
    edges %>% group_by(parent, x_parent) %>% summarise(y = min(y_child), yend = max(y_child), .groups = "drop") %>%
      filter(y != yend) %>% transmute(x = x_parent, xend = x_parent, y = y, yend = yend)
  )

  tips_xy <- d %>%
    filter(isTip) %>%
    transmute(label, base_label, species_short, site, x_branch = x, x_tip = tip_col, y)

  terminals <- tips_xy %>% transmute(x = x_branch, xend = x_tip, y = y, yend = y)
  list(branches = branches, terminals = terminals, tips = tips_xy, x_max = x_max, root_col = root_col, tip_col = tip_col)
}

make_axis <- function(layout, side, y = 0.55) {
  brks <- pretty(c(0, layout$x_max), n = 4)
  brks <- brks[brks >= 0 & brks <= layout$x_max]
  map_x <- function(v) if (side == "left") layout$root_col + v / layout$x_max * TREE_WIDTH else layout$root_col - v / layout$x_max * TREE_WIDTH
  list(
    segments = bind_rows(
      tibble(x = map_x(0), xend = map_x(layout$x_max), y = y, yend = y),
      tibble(x = map_x(brks), xend = map_x(brks), y = y, yend = y - 0.16)
    ),
    labels = tibble(x = map_x(brks), y = y - 0.36, label = format(signif(brks, 3), scientific = FALSE, trim = TRUE))
  )
}

mlst <- prepare_tree(mlst_tree_file, parse_mlst_tips, mlst_samples_to_exclude)
ospc <- prepare_tree(ospc_tree_file, parse_ospc_tips, ospc_samples_to_exclude)
species_meta <- load_species_meta(species_file)
tick_meta <- load_tick_meta(tick_metadata_file)

mlst_ann <- annotate_tips(mlst$tips, species_meta$mlst, tick_meta, by_label = TRUE)
ospc_ann <- annotate_tips(ospc$tips, species_meta$ospc, tick_meta, by_label = FALSE)

mlst_f <- filter_tree_species(mlst$tree, mlst_ann, mlst_species_to_exclude)
ospc_f <- filter_tree_species(ospc$tree, ospc_ann, ospc_species_to_exclude)

n_y <- max(length(mlst_f$tree$tip.label), length(ospc_f$tree$tip.label))
mlst_lay <- layout_tree(mlst_f$tree, mlst_f$tips, "left", c(1, n_y))
ospc_lay <- layout_tree(ospc_f$tree, ospc_f$tips, "right", c(1, n_y))

branch_segments <- bind_rows(mlst_lay$branches, ospc_lay$branches)
terminal_segments <- bind_rows(mlst_lay$terminals, ospc_lay$terminals)
tip_points <- bind_rows(mlst_lay$tips %>% mutate(tree = "MLST"), ospc_lay$tips %>% mutate(tree = "ospC"))

links <- inner_join(
  tip_points %>% filter(tree == "MLST") %>% select(base_label, mlst_label = label, mlst_species = species_short, x_mlst = x_tip, y_mlst = y),
  tip_points %>% filter(tree == "ospC") %>% select(base_label, ospc_label = label, ospc_species = species_short, x_ospc = x_tip, y_ospc = y),
  by = "base_label"
) %>% mutate(link_species = mlst_species)

write_tsv(links, out_links)

labels_df <- tip_points %>%
  filter(base_label %in% unique(links$base_label)) %>%
  mutate(label_x = if_else(tree == "MLST", x_tip + LABEL_OFFSET, x_tip - LABEL_OFFSET), hjust = if_else(tree == "MLST", 0, 1))

axis_left <- make_axis(mlst_lay, "left")
axis_right <- make_axis(ospc_lay, "right")
axis_segments <- bind_rows(axis_left$segments, axis_right$segments)
axis_labels <- bind_rows(axis_left$labels, axis_right$labels)

present_species <- unique(c(tip_points$species_short, links$link_species))
color_values <- species_colors[names(species_colors) %in% present_species]
missing_species <- setdiff(present_species, names(color_values))
if (length(missing_species) > 0) color_values <- c(color_values, setNames(rep("#666666", length(missing_species)), missing_species))
shape_values <- site_shapes[names(site_shapes) %in% unique(tip_points$site)]

p <- ggplot() +
  geom_segment(data = branch_segments, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.34) +
  geom_segment(data = terminal_segments, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.24, linetype = "dotted", color = "grey82") +
  geom_segment(data = axis_segments, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.28, color = "grey20") +
  geom_text(data = axis_labels, aes(x = x, y = y, label = label), size = 2.3, color = "grey20") +
  geom_segment(data = links, aes(x = x_mlst, xend = x_ospc, y = y_mlst, yend = y_ospc, color = link_species), linewidth = 0.28, alpha = 0.55) +
  geom_point(data = tip_points, aes(x = x_tip, y = y, color = species_short, shape = site), size = 1.7, alpha = 0.96) +
  geom_text(data = labels_df, aes(x = label_x, y = y, label = label, hjust = hjust), size = 1.75) +
  annotate("text", x = -INNER_GAP / 2 - TREE_WIDTH / 2, y = n_y + 4, label = "MLST", fontface = "bold", size = 4) +
  annotate("text", x = INNER_GAP / 2 + TREE_WIDTH / 2, y = n_y + 4, label = "ospC", fontface = "bold.italic", size = 4) +
  scale_color_manual(values = color_values, labels = italic_md(names(color_values)), name = "Espèce") +
  scale_shape_manual(values = shape_values, name = "Altitude") +
  coord_cartesian(xlim = c(-INNER_GAP / 2 - TREE_WIDTH - 0.2, INNER_GAP / 2 + TREE_WIDTH + 0.2), ylim = c(-0.35, n_y + 6), clip = "off") +
  theme_void() +
  theme(
    legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
    legend.text = ggtext::element_markdown(size = 7), legend.title = element_text(size = 8),
    legend.key = element_blank(), plot.margin = margin(10, 10, 10, 10)
  )

print(p)
ggsave(out_png, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, units = "in", dpi = PLOT_DPI, bg = "white", limitsize = FALSE)
ggsave(out_jpg, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, units = "in", dpi = PLOT_DPI, bg = "white", limitsize = FALSE)
ggsave(out_svg, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, units = "in", limitsize = FALSE)
cat("Fichiers sauvegardés :\n", out_png, "\n", out_jpg, "\n", out_svg, "\n", out_links, "\n")
