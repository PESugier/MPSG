# Journal des décisions

| N° | Date | Décision | Appliquée à l'étape |
|----|------|----------|---------------------|
| 1 | 2026-11-26 | Tag `v1.1` = état du dépôt au moment de la publication | 1 (fait) |
| 2 | 2026-11-26 | Auteurs, ORCID et mainteneur : bloc D2 ci-dessous | 3 (création de `DESCRIPTION`) |
| 3 | 2026-11-26 | Documentation du jeu de données d'exemple : texte D3 ci-dessous | Étape de création du jeu de données |
| 4 | 2026-11-26 | Les données complètes comptent 3 707 SNPs (le papier indique 3 766, doublons en DL complet inclus). À signaler sobrement dans `paper/README.md` et dans la documentation du jeu de données | 2, puis création du jeu de données |
| 5 | 2026-11-26 | Fonctions héritées : ne conserver que les 3 fonctions appelées par les scripts du papier (`f_cpp_inv_var`, `f_selec_mod_metapsg_gv`, `f_selec_mod_metapsg_gv_parallel`), avec les mêmes arguments et les mêmes sorties, et un message de dépréciation | Étape de la nouvelle interface |
| 6 | 2026-11-26 | Supprimer le code C++ (RcppArmadillo) : inutile une fois l'inversion remplacée par une décomposition spectrale. Les résultats de référence viennent des fichiers du papier, pas d'une ré-exécution du code C++ | Étape d'optimisation |
| 7 | 2026-11-26 | Supprimer l'argument `same_var` (conservé seulement dans la signature de `f_cpp_inv_var`) et la fonction `f_grid_log10` (jamais appelée) | Étape des corrections |
| 8 | 2026-11-26 | Publier un erratum du papier sur GitHub (`paper/ERRATA.md`), après relecture par les co-auteurs. Points identifiés : facteur 1/n_s manquant dans le MIC (section « Hyper-parameters tuning » et annexe A) ; λ1 et λ2 inversés dans l'algorithme 1 ; signe du terme (ρ/2)‖U‖² du lagrangien augmenté ; nombre de SNPs | Étape dédiée, contenu à valider |
| 9 | 2026-11-26 | Étape 3 = empaquetage sans aucun changement d'algorithme. Le C++ est conservé jusqu'à l'étape d'optimisation (Rtools est présent). Dépendances inutilisées ou redondantes retirées (doSNOW, tictoc, abind, Matrix). Version de développement 1.1.9000. Seules les 3 fonctions héritées sont exportées | 3 |

## D2 — Auteurs (champ `Authors@R` de `DESCRIPTION`)

```r
Authors@R: c(
    person("Pierre-Emmanuel", "Sugier", , "pierre-emmanuel.sugier@u-pariscite.fr",
           role = c("aut", "cre"), comment = c(ORCID = "0000-0002-5846-1104")),
    person("Yazdan", "Asgari", role = "aut", comment = c(ORCID = "0000-0001-6993-6956")),
    person("Mohammed", "Sedki", role = "aut", comment = c(ORCID = "0000-0001-5489-9611")),
    person("Thérèse", "Truong", role = "aut", comment = c(ORCID = "0000-0002-2943-6786")),
    person("Benoit", "Liquet", role = "aut", comment = c(ORCID = "0000-0002-8136-2294"))
  )
```

## D3 — Documentation du jeu de données d'exemple (15 gènes)

```
Gene-level summary statistics for 255 SNPs in 15 genes, a subset of the data
analysed in Sugier et al. (2025). Individual-level data come from the CECILE
study, a French population-based case-control study on breast cancer (1,125
cases and 1,172 controls), and from the French case-control studies of the
EPITHYR consortium on thyroid cancer (CATHY, Young-thyr and E3N; 1,129 female
cases and 1,174 female controls). Only women of European ancestry were kept.
SNPs come from genes of 10 non-overlapping candidate pathways; within each gene,
one SNP of each pair in very high linkage disequilibrium (r2 > 0.98) was removed.
For each study and each gene, SNP effects (log odds ratios) and their covariance
matrix were estimated jointly with a Bayesian hierarchical logistic regression
with Gaussian priors (BhGLM package). Only summary statistics are distributed.

The full data contain 3,707 SNPs in 331 genes (the 3,766 SNPs reported in the
paper include duplicated SNPs in complete linkage disequilibrium, removed before
the analysis). They are available in the GitHub repository (folder paper/application).

We thank Pascal Guénel (PI of CECILE and CATHY) and Florent de Vathaire (PI of
Young-thyr) for providing the data.
```

## Bugs connus (à traiter à l'étape des corrections, pas avant)
- Boucle infinie quand un β̂ vaut exactement 0 (bornes de λ infinies dans `find_lambda_max`).
- `f_metapsg_gv_parallel_grid` échoue avec une seule valeur de α (`for (i in 2:n_alpha)`).
- `same_var` ignoré par les fonctions de grille (décision 7 : suppression).
- `f_grid_log10` erronée et jamais appelée (décision 7 : suppression).
- `paper/application/script_metapsg_BC_TC.r` : ligne 99 non commentée (erreur de syntaxe). Archive : décider s'il faut la corriger.
- `paper/application/README.md` : variable `pleio_genes` inexistante. Archive : idem.

## Décisions en attente (à trancher au moment de l'étape concernée)
- Critère d'arrêt de l'ADMM, tolérances par défaut, garde-fou sur rho
- Filtre ad hoc chi2 du papier : retrouver le script utilisé avant toute implémentation
- Grille de lambda par défaut (proposition : alpha = 0.1 à 0.9 par 0.1, eps = 1e-3, 30 valeurs)
