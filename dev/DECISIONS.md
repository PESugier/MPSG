# Journal des décisions

| N° | Date | Décision | Appliquée à l'étape |
|----|------|----------|---------------------|
| 1 | 2026-11-26 | Tag `v1.1` = état du dépôt au moment de la publication | 1 (fait) |
| 2 | 2026-11-26 | Auteurs, ORCID et mainteneur : bloc D2 ci-dessous | 3 (création de `DESCRIPTION`) |
| 3 | 2026-11-26 | Documentation du jeu de données d'exemple : texte D3 ci-dessous | Étape de création du jeu de données |
| 4 | 2026-11-26 | Les données complètes comptent 3 707 SNPs (le papier indique 3 766, doublons en DL complet inclus). À signaler sobrement dans `paper/README.md` et dans la documentation du jeu de données | 2, puis création du jeu de données |

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

## Décisions en attente
- Fonctions héritées à conserver (toutes / celles utilisées par les scripts du papier / aucune)
- Suppression du code C++
- Suppression de `same_var` et de `f_grid_log10`
- Critère d'arrêt, tolérances par défaut, garde-fou sur rho
- Intégration du filtre ad hoc chi2 du papier
- Grille de lambda par défaut
- Écarts papier/code (formule du MIC, algorithme 1) : signaler aux co-auteurs ?
