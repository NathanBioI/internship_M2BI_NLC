# Analyse statistique genASPE / Borrelia
# Entrée : genASPE_table_prevalence.tsv

pkgs <- c("readr", "dplyr", "tidyr", "stringr", "lubridate", "purrr", "tibble")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# Paramètres
tsv_path <- "results/table_infection.tsv"
out_dir <- "stats_genASPE_Borrelia"
dir.create(out_dir, showWarnings = FALSE)

exclude_site <- "S15"
alt <- c(S12 = 400, S13 = 800, S14 = 1200, S15 = 1600)

FISHER_WORKSPACE <- 2e8
FISHER_CONTROL_MULT <- 100

species_cols <- c(
  "Borrelia afzelii",
  "Borrelia burgdorferi sensu stricto",
  "Borrelia garinii",
  "Borrelia lusitaniae",
  "Borrelia valaisiana",
  "Borrelia miyamotoi"
)

# Fonctions 
num <- function(x) readr::parse_number(stringr::str_replace_all(as.character(x), ",", "."))

read_date <- function(x) {
  d <- suppressWarnings(lubridate::dmy(as.character(x)))
  i <- is.na(d)
  d[i] <- suppressWarnings(lubridate::ymd(as.character(x)[i]))
  as.Date(d)
}

# Saison
saison_mois <- function(d) dplyr::case_when(
  lubridate::month(d) %in% 3:5 ~ "Printemps",
  lubridate::month(d) %in% 6:8 ~ "Ete",
  lubridate::month(d) %in% 9:11 ~ "Automne",
  lubridate::month(d) %in% c(12, 1, 2) ~ "Hiver"
)

w <- function(x, file) readr::write_tsv(x, file.path(out_dir, file), na = "")

# Test exact de Fisher
# Si le tableau est 2x2, R renvoie aussi un odds ratio et son intervalle de confiance.
# Pour les tableaux r x c, il n'y a pas d'odds ratio unique : les colonnes sont mises à NA.
test_tab <- function(tab, analyse) {
  tab <- as.matrix(tab)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]

  if (all(dim(tab) == c(2, 2))) {
    ft <- stats::fisher.test(tab, simulate.p.value = FALSE)
    test <- "Fisher exact 2x2"
    p <- ft$p.value
    odds_ratio <- unname(ft$estimate)
    conf_low <- unname(ft$conf.int[1])
    conf_high <- unname(ft$conf.int[2])
  } else {
    ft <- stats::fisher.test(
      tab,
      simulate.p.value = FALSE,
      workspace = FISHER_WORKSPACE,
      control = list(mult = FISHER_CONTROL_MULT)
    )
    test <- "Fisher exact rxc"
    p <- ft$p.value
    odds_ratio <- NA_real_
    conf_low <- NA_real_
    conf_high <- NA_real_
  }

  tibble(
    analyse,
    test,
    p_value = p,
    odds_ratio = odds_ratio,
    conf_low = conf_low,
    conf_high = conf_high,
    n_total = sum(tab)
  )
}

assoc_test <- function(df, row_var, col_var, label) {
  d <- df %>% filter(!is.na(.data[[row_var]]), !is.na(.data[[col_var]]))
  test_tab(table(d[[row_var]], d[[col_var]]), label)
}

# Résidus
stdres <- function(df, factor_var, label) {
  d <- df %>% filter(!is.na(.data[[factor_var]]), !is.na(espece))
  tab <- table(d$espece, d[[factor_var]])
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]

  as.data.frame(as.table(suppressWarnings(stats::chisq.test(tab, correct = FALSE)$stdres))) %>%
    as_tibble() %>%
    setNames(c("espece", "niveau", "std_residual")) %>%
    mutate(analyse = paste(label, "espece x", factor_var), facteur = factor_var, .before = 1)
}

# Comparaisons deux à deux
pairwise_comp <- function(df, factor_var, label) {
  d <- df %>% filter(!is.na(.data[[factor_var]]), !is.na(espece))
  lev <- sort(unique(as.character(d[[factor_var]])))

  purrr::map_dfr(combn(lev, 2, simplify = FALSE), function(p) {
    dd <- d %>% filter(.data[[factor_var]] %in% p)
    test_tab(
      table(dd$espece, dd[[factor_var]]),
      paste(label, factor_var, paste(p, collapse = " vs "), sep = " - ")
    ) %>%
      mutate(facteur = factor_var, comparaison = paste(p, collapse = " vs "), .after = analyse)
  })
}

# Lecture du TSV
raw <- readr::read_tsv(
  tsv_path,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  name_repair = "unique"
) %>%
  select(where(~ !all(is.na(.x) | .x == "")))

# Mise en forme
ticks <- raw %>%
  filter(stringr::str_detect(as.character(`n°ADN EPIA`), "^ADN\\s*EPPAT")) %>%
  transmute(
    sample_id = as.character(`n°ADN EPIA`),
    site = stringr::str_extract(as.character(Site), "S[0-9]+"),
    date_collecte = read_date(`date collecte`),
    annee = factor(lubridate::year(date_collecte)),
    saison = factor(saison_mois(date_collecte), levels = c("Printemps", "Ete", "Automne", "Hiver")),
    across(all_of(species_cols), num)
  ) %>%
  mutate(
    site = factor(site, levels = names(alt)),
    altitude_m = unname(alt[as.character(site)]),
    site_test = factor(
      if_else(as.character(site) %in% exclude_site, NA_character_, as.character(site)),
      levels = setdiff(names(alt), exclude_site)
    )
  )

borrelia <- ticks %>%
  select(sample_id, site, site_test, altitude_m, date_collecte, annee, saison, all_of(species_cols)) %>%
  pivot_longer(all_of(species_cols), names_to = "espece", values_to = "n") %>%
  mutate(
    n = replace_na(as.integer(round(n)), 0L),
    complexe = if_else(espece == "Borrelia miyamotoi", "Borrelia miyamotoi", "Bbsl")
  ) %>%
  filter(n > 0) %>%
  uncount(n, .id = "copie") %>%
  mutate(
    espece = factor(espece, levels = species_cols),
    complexe = factor(complexe, levels = c("Bbsl", "Borrelia miyamotoi"))
  )

bbsl <- borrelia %>% filter(complexe == "Bbsl")
borrelia_site <- borrelia %>% filter(!is.na(site_test))
bbsl_site <- bbsl %>% filter(!is.na(site_test))

# Sorties
resume <- tibble(
  indicateur = c("Bacteries Borrelia identifiees", "Bacteries Bbsl identifiees", "Bacteries B. miyamotoi"),
  valeur = c(nrow(borrelia), nrow(bbsl), sum(borrelia$espece == "Borrelia miyamotoi"))
)

effectifs <- borrelia %>%
  count(espece, name = "n_bacteries") %>%
  mutate(prop = n_bacteries / sum(n_bacteries)) %>%
  arrange(desc(n_bacteries))

# Tests globaux
tests <- bind_rows(
  assoc_test(borrelia_site, "espece", "site_test", "Borrelia  : espece x site"),
  assoc_test(borrelia, "espece", "annee", "Borrelia  : espece x annee"),
  assoc_test(borrelia, "espece", "saison", "Borrelia  : espece x saison"),
  assoc_test(borrelia_site, "complexe", "site_test", "Bbsl/miyamotoi x site"),
  assoc_test(borrelia, "complexe", "annee", "Bbsl/miyamotoi x annee"),
  assoc_test(borrelia, "complexe", "saison", "Bbsl/miyamotoi x saison"),
  assoc_test(bbsl_site, "espece", "site_test", "Bbsl  : espece x site"),
  assoc_test(bbsl, "espece", "annee", "Bbsl  : espece x annee"),
  assoc_test(bbsl, "espece", "saison", "Bbsl  : espece x saison")
)

# Résidus
residus <- bind_rows(
  stdres(borrelia_site, "site_test", "Borrelia "),
  stdres(borrelia, "annee", "Borrelia "),
  stdres(borrelia, "saison", "Borrelia "),
  stdres(bbsl_site, "site_test", "Bbsl "),
  stdres(bbsl, "annee", "Bbsl "),
  stdres(bbsl, "saison", "Bbsl ")
)

# Comparaisons deux à deux entre niveaux de site, année et saison
pairs <- bind_rows(
  pairwise_comp(borrelia_site, "site_test", "Borrelia "),
  pairwise_comp(borrelia, "annee", "Borrelia "),
  pairwise_comp(borrelia, "saison", "Borrelia "),
  pairwise_comp(bbsl_site, "site_test", "Bbsl "),
  pairwise_comp(bbsl, "annee", "Bbsl "),
  pairwise_comp(bbsl, "saison", "Bbsl ")
)

# Écriture des fichiers
w(resume, "00_resume_bacteries.tsv")
w(borrelia, "01_bacteries_Borrelia_long.tsv")
w(bbsl, "02_bacteries_Bbsl_long.tsv")
w(effectifs, "03_effectifs_especes.tsv")
w(tests, "04_tests_globaux.tsv")
w(residus, "05_residus_standardises.tsv")
w(pairs, "06_pairwise_composition.tsv")

cat("\nTerminé","\n\n")
print(resume)
print(effectifs)
print(tests)
