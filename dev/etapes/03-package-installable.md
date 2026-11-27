# Étape 3 — Rendre le dépôt installable comme package R (sans changer l'algorithme)

## Objectif
Transformer `R/` et `src/` en un package R installable et vérifiable, en gardant
**exactement** les calculs du code publié (tag `v1.1`), C++ compris. Deux
contrôles de non-régression :
- un test du package sur 15 gènes, comparé à une référence produite par le code
  original lui-même ;
- la reproduction complète du tableau 6 du papier (331 gènes), comparée aux
  résultats sauvegardés dans `paper/application/results_mpsg/`.

Fichiers fournis avec cette étape (non suivis au départ) :
`dev/verification/make_reference_15genes.R`, `dev/verification/ref_15genes_linux.rds`,
`dev/verification/reproduce_table6.R`, `tests/testthat.R`,
`tests/testthat/test-reference-15genes.R`. `dev/DECISIONS.md` est mis à jour.
**Ne pas modifier ces fichiers** : s'ils posent problème, s'arrêter et le signaler.

## Conventions
- Pour appeler R, utiliser `Rscript` et `Rcmd` (`Rcmd build`, `Rcmd check`,
  `Rcmd INSTALL`). Sous PowerShell, `R` est un alias d'`Invoke-History` : ne
  jamais taper `R CMD`.
- Préserver les fins de ligne existantes des fichiers modifiés (les fichiers de
  `R/` sont en CRLF).

## 0. Pré-vérifications (s'arrêter si l'une échoue)
1. Branche courante `dev`, contenant le commit de l'étape 2
   (« Move paper material… ») ; `origin/dev` existe et est à jour.
2. `git status` : seuls les fichiers fournis ci-dessus sont nouveaux, et
   `dev/DECISIONS.md` est modifié.
3. `Rscript --version` répond (R 4.4.2) ; `Rscript -e "cat(pkgbuild::has_build_tools())"` renvoie TRUE.
4. Reporter `Rscript -e ".libPaths()"`. Les installations se font dans la
   première bibliothèque accessible en écriture (normalement la bibliothèque
   utilisateur) ; si aucune ne l'est, s'arrêter et le signaler.

## 1. Packages R nécessaires
Installer depuis CRAN (`repos = "https://cloud.r-project.org"`) **uniquement ceux
qui manquent** parmi : roxygen2, Rcpp, RcppArmadillo, foreach, doParallel, abind.
Ne rien installer d'autre. Reporter ce qui a été installé et les versions.

## 2. Référence du code original produite sur cette machine
1. `Rscript dev/verification/make_reference_15genes.R dev/verification/ref_15genes_windows.rds`
   (le script extrait le code du tag `v1.1` et le compile avec sourceCpp).
2. Comparer `ref_15genes_windows.rds` à `ref_15genes_linux.rds` (produit par
   Claude sur Linux) pour chacune des 4 configurations : écart maximal sur
   `lambda`, `MIC` et sur les coefficients `beta_select` ; égalité exacte de
   `ng_hit`, `nv_hit` et `iter`. Reporter le tout. Si `ng_hit` ou `nv_hit`
   diffèrent, s'arrêter et le signaler.
3. Créer `tests/testthat/fixtures/` et y copier
   `dev/verification/ref_15genes_windows.rds` sous le nom `ref_15genes.rds`.

## 3. Modifications pour l'empaquetage
Aucune autre modification de code que celles-ci.

### 3.1 `DESCRIPTION` (nouveau, contenu exact)
```
Package: MPSG
Type: Package
Title: Meta-Analysis Model for Pleiotropy Selection with Group Structure
Version: 1.1.9000
Authors@R: c(
    person("Pierre-Emmanuel", "Sugier", , "pierre-emmanuel.sugier@u-pariscite.fr",
           role = c("aut", "cre"), comment = c(ORCID = "0000-0002-5846-1104")),
    person("Yazdan", "Asgari", role = "aut", comment = c(ORCID = "0000-0001-6993-6956")),
    person("Mohammed", "Sedki", role = "aut", comment = c(ORCID = "0000-0001-5489-9611")),
    person("Thérèse", "Truong", role = "aut", comment = c(ORCID = "0000-0002-2943-6786")),
    person("Benoit", "Liquet", role = "aut", comment = c(ORCID = "0000-0002-8136-2294"))
  )
Description: Detects cross-phenotype (pleiotropic) genetic associations from
    summary statistics of several studies. A penalised multivariate
    meta-analysis model with group structure selects, jointly across all
    studies, groups of variants (such as genes or pathways) and variants
    within groups. The model is fitted with an alternating direction method
    of multipliers (ADMM) algorithm and tuned with a modified information
    criterion. The method is described in Sugier et al. (2025)
    <doi:10.1093/biostatistics/kxaf037>.
License: GPL (>= 3)
URL: https://github.com/PESugier/MPSG
BugReports: https://github.com/PESugier/MPSG/issues
Depends:
    R (>= 4.1.0)
Imports:
    doParallel,
    foreach,
    parallel,
    Rcpp
LinkingTo:
    Rcpp,
    RcppArmadillo
Suggests:
    testthat (>= 3.0.0)
Config/testthat/edition: 3
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
```

### 3.2 `R/MPSG-package.R` (nouveau, contenu exact)
```r
#' MPSG: Meta-Analysis Model for Pleiotropy Selection with Group Structure
#'
#' Implementation of the method described in Sugier et al. (2025),
#' Biostatistics, 26(1), kxaf037, doi:10.1093/biostatistics/kxaf037.
#'
#' @keywords internal
#' @useDynLib MPSG, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
"_PACKAGE"
```

### 3.3 Lignes à retirer des fichiers R (et rien d'autre)
- `R/metapsg_gv_solver.r` : `library(Matrix)`, `library(Rcpp)`,
  `library(RcppArmadillo)`, `sourceCpp("arma_inv.cpp")`, `source("metapsg_sub.r")`.
- `R/metapsg_sub.r` : `sourceCpp("arma_inv.cpp")`, `library(abind)`.
- `R/metapsg_parallel.r` : `library(tictoc)`, `library(doParallel)`,
  `library(doSNOW)`, `source("metapsg_gv_solver.r")` (en tête de fichier), les
  deux lignes `registerDoSNOW(...)` (le backend doParallel enregistré juste
  avant reste en place), et les deux `source("metapsg_gv_solver.r")` situés dans
  les blocs `foreach`.
- Dans `R/metapsg_parallel.r`, ajouter l'argument `.packages = "MPSG"` aux deux
  appels `foreach(...)`, pour que les processus parallèles chargent le package
  au lieu de sourcer les fichiers.
- Les `print()` et les bugs connus (voir `dev/DECISIONS.md`) restent tels quels.

### 3.4 Fichiers à supprimer (`git rm`)
`src/RcppExports.o`, `src/arma_inv.o`, `src/METAPSG.dll` (binaires compilés),
`R/arma_inv.cpp` (doublon de `src/arma_inv.cpp`), `R/README.md` (un dossier
`R/` de package ne doit contenir que du code R).

### 3.5 `src/Makevars` et `src/Makevars.win`
Supprimer le bloc de commentaires et la ligne `CXX_STD = CXX11` (le CRAN demande
de la retirer). Garder les lignes `PKG_CXXFLAGS` et `PKG_LIBS` telles quelles.

### 3.6 Interface Rcpp
`Rscript -e "Rcpp::compileAttributes()"` pour régénérer `R/RcppExports.R` et
`src/RcppExports.cpp` avec le préfixe `_MPSG_` (au lieu de `_METAPSG_`).
`arma_inv` ne doit pas être exportée.

### 3.7 Documentation et exports
Ajouter un bloc roxygen avec `@export` uniquement pour les trois fonctions
utilisées par les scripts du papier : `f_cpp_inv_var`, `f_selec_mod_metapsg_gv`,
`f_selec_mod_metapsg_gv_parallel`. Chaque bloc contient : un titre, une phrase
indiquant qu'il s'agit de l'interface utilisée pour Sugier et al. (2025), la
description de chaque argument (d'après le code, sans inventer de
comportement), la description de la valeur retournée. Aucune modification du
corps des fonctions. Puis `Rscript -e "roxygen2::roxygenise()"` (génère
`NAMESPACE` et `man/`).

### 3.8 `.Rbuildignore` et `.gitignore` (nouveaux, contenu exact)
`.Rbuildignore` :
```
^CLAUDE\.md$
^dev$
^paper$
^.*\.Rproj$
^\.Rproj\.user$
```
`.gitignore` :
```
.Rproj.user
.Rhistory
.RData
.Ruserdata
src/*.o
src/*.so
src/*.dll
*.tar.gz
*.Rcheck/
```

## 4. Construction et vérification
1. `Rcmd build .` puis, avec la variable d'environnement `NOT_CRAN=true` (pour
   que le test de non-régression s'exécute), `Rcmd check --no-manual MPSG_1.1.9000.tar.gz`.
2. Objectif : 0 ERROR, 0 WARNING. Reporter **mot pour mot** chaque WARNING et
   chaque NOTE.
3. Seule correction autorisée sans demander : la NOTE « no visible binding for
   global variable 'i' » (déclarer `i <- NULL` avant le `foreach` concerné).
   Pour tout le reste, s'arrêter et signaler.

## 5. Reproduction complète du tableau 6
1. Installer le package : `Rcmd INSTALL .`
2. Lancer en arrière-plan, avec journal :
   `Rscript dev/verification/reproduce_table6.R > dev/verification/reproduce_table6.log 2>&1`
   (les deux analyses, pondérée et non pondérée, tournent l'une après l'autre ;
   cela peut prendre longtemps). Surveiller le journal jusqu'à la fin.
3. Reporter le contenu de `dev/verification/reproduce_table6_results.txt`,
   le nombre de cœurs utilisés et la durée. Si les gènes ou les nombres de SNPs
   diffèrent de la référence, ne rien corriger : le signaler.

## 6. Indexation
`git add` des fichiers créés, modifiés et supprimés, y compris
`dev/verification/` (scripts, deux fichiers .rds, résultats, journal) et
`tests/`. Ne pas indexer `*.tar.gz` ni `MPSG.Rcheck/` (ignorés par `.gitignore`).
Pas de commit.

## 7. Rapport
`dev/rapports/03-package-installable.md`, au format défini dans `CLAUDE.md`.
Message de commit proposé :
`Make the repository an installable R package (no change to the algorithm)`
