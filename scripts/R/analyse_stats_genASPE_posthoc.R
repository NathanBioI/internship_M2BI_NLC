set.seed(123)

pkgs <- c("readr", "dplyr", "tidyr", "stringr", "purrr", "tibble")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

input_dir <- if (dir.exists("stats_genASPE_Borrelia_bacteries_only")) "stats_genASPE_Borrelia_bacteries_only" else "."
out_dir <- "stats_genASPE_posthoc_bacteries_only"
dir.create(out_dir, showWarnings = FALSE)

FISHER_WORKSPACE <- 2e8    # augmenter si fisher.test renvoie une erreur de workspace
FISHER_CONTROL_MULT <- 100 # augmente la place reservee aux chemins du calcul exact

r <- function(file) {
  path <- file.path(input_dir, file)
  if (!file.exists(path)) stop("Fichier introuvable : ", path)
  readr::read_tsv(path, show_col_types = FALSE)
}

w <- function(x, file) readr::write_tsv(x, file.path(out_dir, file), na = "")


test_tab <- function(tab, analyse) {
  tab <- as.matrix(tab)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]

  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(tibble(analyse, test = NA_character_, p_value = NA_real_,
                  n_total = sum(tab), min_expected = NA_real_))
  }

  # min_expected est conserve uniquement comme indicateur descriptif
  # des faibles effectifs attendus. Il ne sert plus a choisir entre chi2 et Fisher,
  # car tous les tableaux sont testes avec Fisher exact.
  chi0 <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  min_exp <- min(chi0$expected)

  if (all(dim(tab) == c(2, 2))) {
    p <- stats::fisher.test(tab, simulate.p.value = FALSE)$p.value
    test <- "Fisher exact 2x2"
  } else {
    p <- stats::fisher.test(
      tab,
      simulate.p.value = FALSE,
      workspace = FISHER_WORKSPACE,
      control = list(mult = FISHER_CONTROL_MULT)
    )$p.value
    test <- paste0("Fisher exact rxc workspace=", FISHER_WORKSPACE)
  }

  tibble(analyse, test, p_value = p, n_total = sum(tab),
         min_expected = min_exp)
}

global_test <- function(df, factor_var, label) {
  d <- df %>% filter(!is.na(.data[[factor_var]]), !is.na(espece))
  test_tab(table(d$espece, d[[factor_var]]), paste(label, "espece x", factor_var))
}

stdres <- function(df, factor_var, label) {
  d <- df %>% filter(!is.na(.data[[factor_var]]), !is.na(espece))
  tab <- table(d$espece, d[[factor_var]])
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]

  if (nrow(tab) < 2 || ncol(tab) < 2) return(tibble())

  as.data.frame(as.table(suppressWarnings(stats::chisq.test(tab, correct = FALSE)$stdres))) %>%
    as_tibble() %>%
    setNames(c("espece", "niveau", "std_residual")) %>%
    mutate(analyse = paste(label, "espece x", factor_var), facteur = factor_var, .before = 1)
}

pairwise_comp <- function(df, factor_var, label) {
  d <- df %>% filter(!is.na(.data[[factor_var]]), !is.na(espece))
  lev <- sort(unique(as.character(d[[factor_var]])))
  if (length(lev) < 2) return(tibble())

  purrr::map_dfr(combn(lev, 2, simplify = FALSE), function(p) {
    dd <- d %>% filter(.data[[factor_var]] %in% p)
    test_tab(table(dd$espece, dd[[factor_var]]),
             paste(label, factor_var, paste(p, collapse = " vs "), sep = " - ")) %>%
      mutate(facteur = factor_var, comparaison = paste(p, collapse = " vs "), .after = analyse)
  })
}

presence_bacterie <- function(df, factor_var, label) {
  d <- df %>% filter(!is.na(.data[[factor_var]]), !is.na(espece))

  grid <- tidyr::expand_grid(
    espece = sort(unique(as.character(d$espece))),
    niveau = sort(unique(as.character(d[[factor_var]])))
  )

  purrr::pmap_dfr(grid, function(espece, niveau) {
    tab <- table(
      espece_cible = d$espece == espece,
      niveau_cible = as.character(d[[factor_var]]) == niveau
    )

    test_tab(tab, paste(label, espece, "x", factor_var, niveau)) %>%
      mutate(
        jeu = label,
        espece = espece,
        facteur = factor_var,
        niveau = niveau,
        n_espece_niveau = sum(d$espece == espece & as.character(d[[factor_var]]) == niveau),
        n_espece_total = sum(d$espece == espece),
        n_niveau = sum(as.character(d[[factor_var]]) == niveau),
        prop_dans_niveau = n_espece_niveau / n_niveau,
        .after = analyse
      )
  })
}

borrelia <- r("stats_genASPE_Borrelia/01_bacteries_Borrelia_long.tsv")
bbsl <- r("stats_genASPE_Borrelia/02_bacteries_Bbsl_long.tsv")

borrelia_site <- borrelia %>% filter(!is.na(site_test))
bbsl_site <- bbsl %>% filter(!is.na(site_test))

tests_globaux <- bind_rows(
  global_test(borrelia_site, "site_test", "Borrelia 174"),
  global_test(borrelia, "annee", "Borrelia 174"),
  global_test(borrelia, "saison", "Borrelia 174"),
  global_test(bbsl_site, "site_test", "Bbsl 152"),
  global_test(bbsl, "annee", "Bbsl 152"),
  global_test(bbsl, "saison", "Bbsl 152")
)

residus <- bind_rows(
  stdres(borrelia_site, "site_test", "Borrelia 174"),
  stdres(borrelia, "annee", "Borrelia 174"),
  stdres(borrelia, "saison", "Borrelia 174"),
  stdres(bbsl_site, "site_test", "Bbsl 152"),
  stdres(bbsl, "annee", "Bbsl 152"),
  stdres(bbsl, "saison", "Bbsl 152")
)

presence <- bind_rows(
  presence_bacterie(borrelia_site, "site_test", "Borrelia 174"),
  presence_bacterie(borrelia, "annee", "Borrelia 174"),
  presence_bacterie(borrelia, "saison", "Borrelia 174"),
  presence_bacterie(bbsl_site, "site_test", "Bbsl 152"),
  presence_bacterie(bbsl, "annee", "Bbsl 152"),
  presence_bacterie(bbsl, "saison", "Bbsl 152")
)

pairwise <- bind_rows(
  pairwise_comp(borrelia_site, "site_test", "Borrelia 174"),
  pairwise_comp(borrelia, "annee", "Borrelia 174"),
  pairwise_comp(borrelia, "saison", "Borrelia 174"),
  pairwise_comp(bbsl_site, "site_test", "Bbsl 152"),
  pairwise_comp(bbsl, "annee", "Bbsl 152"),
  pairwise_comp(bbsl, "saison", "Bbsl 152")
)

w(tests_globaux, "01_tests_globaux.tsv")
w(residus, "02_residus_standardises.tsv")
w(presence, "03_tests_presence_especes_bacteries.tsv")
w(pairwise, "04_pairwise_composition.tsv")

cat("\nTermine. Dossier :", out_dir, "\n")
cat("Tests de presence : 03_tests_presence_especes_bacteries.tsv\n\n")
print(tests_globaux)
