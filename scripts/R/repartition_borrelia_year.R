# Répartition proportionnelle des espèces de Borrelia par année
# Entrée : results/table_infection.tsv
# Sorties : TSV + PNG + PDF

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, ggtext, scales)

# Paramètres

input_file <- "results/table_infection.tsv"

output_table <- "results/repartition_2023_2024.tsv"
output_plot_png <- "results/repartition_2023_2024.png"
output_plot_pdf <- "results/repartition_2023_2024.pdf"

species_cols <- c(
  "Borrelia afzelii",
  "Borrelia burgdorferi sensu stricto",
  "Borrelia garinii",
  "Borrelia lusitaniae",
  "Borrelia valaisiana",
  "Borrelia miyamotoi"
)

annees <- c("2023", "2024")
colors_year <- c("2023" = "#66c2a5", "2024" = "#fc8d62")

plot_width <- 22
plot_height <- 9.5
bar_height <- 0.32
dodge_offset <- 0.18
species_spacing <- 1.35
ic_tick_height <- 0.11

# Fonctions

format_species <- function(x) {
  x <- as.character(x)
  if_else(str_detect(x, "^Borrelia\\s+"), paste0("*", x, "*"), x)
}

binom_ci_pct <- function(x, n) {
  if (is.na(x) || is.na(n) || n == 0 || x == 0) return(c(NA_real_, NA_real_))
  binom.test(x, n)$conf.int * 100
}

format_label_ic <- function(prop, ci_inf, ci_sup, n) {
  if (is.na(ci_inf) || is.na(ci_sup)) {
    sprintf("**%.1f%%** (n = %d)", prop, n)
  } else {
    sprintf("**%.1f%%** [IC95%% : %.1f–%.1f] (n = %d)", prop, ci_inf, ci_sup, n)
  }
}

safe_presence <- function(x) {
  as.integer(replace_na(suppressWarnings(as.numeric(x)), 0))
}

# Lecture + comptage

dir.create("results", showWarnings = FALSE, recursive = TRUE)

counts <- read_tsv(input_file, show_col_types = FALSE) %>%
  mutate(annee = str_extract(as.character(`date collecte`), "20\\d{2}")) %>%
  filter(annee %in% annees) %>%
  mutate(across(all_of(species_cols), safe_presence)) %>%
  select(annee, all_of(species_cols)) %>%
  pivot_longer(all_of(species_cols), names_to = "espece", values_to = "presence") %>%
  group_by(espece, annee) %>%
  summarise(effectif = sum(presence), .groups = "drop") %>%
  complete(espece = species_cols, annee = annees, fill = list(effectif = 0)) %>%
  group_by(espece) %>%
  mutate(
    total_espece = sum(effectif),
    proportion_pct = if_else(total_espece > 0, effectif / total_espece * 100, 0)
  ) %>%
  ungroup() %>%
  filter(total_espece > 0) %>%
  rowwise() %>%
  mutate(
    ci = list(binom_ci_pct(effectif, total_espece)),
    ci_inf_pct = ci[[1]],
    ci_sup_pct = ci[[2]]
  ) %>%
  ungroup() %>%
  select(-ci)

write_tsv(counts, output_table)

effectif_total <- counts %>%
  distinct(espece, total_espece) %>%
  summarise(n = sum(total_espece), .groups = "drop") %>%
  pull(n)

# Positions graphiques

order_species <- counts %>%
  distinct(espece, total_espece) %>%
  arrange(desc(total_espece))

miyamotoi <- order_species %>% filter(str_detect(espece, regex("miyamoto", ignore_case = TRUE))) %>% pull(espece)
others <- order_species %>% filter(!str_detect(espece, regex("miyamoto", ignore_case = TRUE))) %>% pull(espece)
visual_order <- c(others, miyamotoi)

df_y <- tibble(
  espece = rev(visual_order),
  y_base = seq_along(visual_order) * species_spacing
)

plot_data <- counts %>%
  filter(effectif > 0) %>%
  left_join(df_y, by = "espece") %>%
  group_by(espece) %>%
  mutate(
    n_groups = n_distinct(annee),
    y_pos = case_when(
      n_groups == 1 ~ y_base,
      annee == "2023" ~ y_base + dodge_offset,
      annee == "2024" ~ y_base - dodge_offset,
      TRUE ~ y_base
    )
  ) %>%
  ungroup() %>%
  mutate(
    label_bar = mapply(format_label_ic, proportion_pct, ci_inf_pct, ci_sup_pct, effectif),
    xmin = 0,
    xmax = proportion_pct,
    ymin = y_pos - bar_height / 2,
    ymax = y_pos + bar_height / 2,
    ic_ymin = y_pos - ic_tick_height,
    ic_ymax = y_pos + ic_tick_height,
    label_x = pmin(if_else(is.na(ci_sup_pct), proportion_pct, ci_sup_pct) + 4, 112)
  )

plot_ic <- plot_data %>% filter(!is.na(ci_inf_pct), !is.na(ci_sup_pct))
labels_species <- df_y %>% mutate(espece_md = format_species(espece))
y_min <- min(df_y$y_base) - 0.75
y_max <- max(df_y$y_base) + 0.75

# Graphique

p <- ggplot() +
  geom_rect(
    data = plot_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = annee),
    color = "black", alpha = 0.9
  ) +
  geom_segment(
    data = plot_ic,
    aes(x = ci_inf_pct, xend = ci_sup_pct, y = y_pos, yend = y_pos),
    linewidth = 0.75, color = "black"
  ) +
  geom_segment(
    data = plot_ic,
    aes(x = ci_inf_pct, xend = ci_inf_pct, y = ic_ymin, yend = ic_ymax),
    linewidth = 0.75, color = "black"
  ) +
  geom_segment(
    data = plot_ic,
    aes(x = ci_sup_pct, xend = ci_sup_pct, y = ic_ymin, yend = ic_ymax),
    linewidth = 0.75, color = "black"
  ) +
  geom_richtext(
    data = plot_data,
    aes(x = label_x, y = y_pos, label = label_bar),
    hjust = 0, size = 4.4, color = "black", fill = NA, label.color = NA
  ) +
  geom_richtext(
    data = labels_species,
    aes(x = 0, y = y_base, label = espece_md),
    inherit.aes = FALSE, hjust = 1, nudge_x = -1.8,
    fill = NA, label.color = NA, size = 5.8, color = "black"
  ) +
  scale_fill_manual(values = colors_year, drop = FALSE) +
  coord_cartesian(xlim = c(0, 160), ylim = c(y_min, y_max), clip = "off") +
  scale_x_continuous(labels = label_percent(scale = 1), breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = NULL) +
  theme_minimal(base_size = 17) +
  labs(
    title = "Répartition proportionnelle de chaque espèce entre 2023 et 2024",
    subtitle = paste0("Pourcentages calculés au sein de chaque espèce ; IC95% binomiaux ; effectif total = ", effectif_total),
    x = "Proportion au sein de chaque espèce",
    y = NULL,
    fill = "Année"
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.title.x = element_text(size = 17, face = "bold", color = "black", margin = margin(t = 12)),
    plot.title = element_text(face = "bold", size = 21, color = "black"),
    plot.subtitle = element_text(size = 14, color = "black", margin = margin(b = 12)),
    legend.position = "top",
    legend.title = element_text(size = 15, face = "bold", color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 25, r = 520, b = 25, l = 360)
  )

print(p)
print(counts)

ggsave(output_plot_png, p, width = plot_width, height = plot_height, dpi = 300, bg = "white")
ggsave(output_plot_pdf, p, width = plot_width, height = plot_height, bg = "white")
