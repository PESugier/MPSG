# Étape 4 — Corrections des bugs connus

## Objectif
Corriger les bugs listés dans `dev/DECISIONS.md` (section « Bugs connus »),
décisions 7 et 10. Avec les réglages utilisés dans le papier, **les résultats ne
doivent pas changer** : le test de référence sur 15 gènes et la reproduction
complète du tableau 6 doivent redonner exactement les mêmes valeurs qu'à l'étape 3.

Fichiers fournis avec cette étape : `CLAUDE.md` (mis à jour : section
Environnement), `dev/DECISIONS.md` (mis à jour), `dev/etapes/04-corrections.md`,
`dev/etapes/04-reference.diff`, `tests/testthat/test-bugfixes.R`.
**Ne pas modifier le fichier de test fourni** ; en cas de problème, s'arrêter et
le signaler.

`dev/etapes/04-reference.diff` montre les modifications de `R/` attendues. Il a
été produit et testé par Claude sur une reconstruction du dépôt après l'étape 3.
Le contexte peut différer légèrement du vôtre (blocs roxygen, espaces) : ne pas
l'appliquer avec `git apply`, mais s'en servir comme référence pour faire les
modifications décrites ci-dessous, à la main.

## 0. Pré-vérifications (s'arrêter si l'une échoue)
1. Branche `dev` ; HEAD = commit de l'étape 3 (« Make the repository an
   installable R package… ») ; à jour avec `origin/dev`.
2. `git status` : `CLAUDE.md` et `dev/DECISIONS.md` modifiés ; nouveaux :
   `dev/etapes/04-corrections.md`, `dev/etapes/04-reference.diff`,
   `tests/testthat/test-bugfixes.R`. Rien d'autre.

## 1. Corrections dans `R/`
Règles : uniquement les modifications ci-dessous ; dans chaque appel de fonction
modifié, passer les arguments **par nom** ; conserver les fins de ligne CRLF.

### 1.1 Boucle infinie quand un β̂ vaut exactement 0 (`R/metapsg_sub.r`)
- Ajouter, juste avant `find_lambda_max_in_g_for_vsnorm`, la fonction :
```r
# Largest finite value of x (NA if there is none). Bounds are infinite or
# undefined when an observed effect is exactly 0; such values are ignored.
f_max_finite <- function(x){
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}
```
- Dans `find_lambda_max_in_g_for_vsnorm` et `find_lambda_max_in_g_for_gnorm` :
  `v_lambda_max[g] <- max(v_this_g, na.rm=TRUE)` devient
  `v_lambda_max[g] <- f_max_finite(v_this_g)`, et
  `result <- max(v_lambda_max)/coeff` devient
  `result <- f_max_finite(v_lambda_max)/coeff`.
  Sans effet quand aucun effet n'est nul : `max(..., na.rm = TRUE)` écartait déjà
  les NaN, seuls les Inf sont nouvellement écartés.

### 1.2 `find_lambda_max` (`R/metapsg_sub.r`) : remplacer la fonction entière par
```r
find_lambda_max <- function(sumstat, n_s, alpha=0.1, w_g, w_p, f_test=1.5, f_selec=1.2,
                            rho=1, t_inc=2, t_dec=2, mu=10, max_iter=1000, tol_abs=1e-03, tol_rel=1e-02,
                            max_halving=100){
  result <- 0
  solve_at <- function(lambda) metapsg_gv(sumstat, n_s, alpha=alpha, lambda=lambda, w_p=w_p, w_g=w_g,
                                          rho=rho, t_inc=t_inc, t_dec=t_dec, mu=mu, max_iter=max_iter,
                                          tol_abs=tol_abs, tol_rel=tol_rel)
  l_max <- suppressWarnings(min(c(find_lambda_max_in_g_for_gnorm(sumstat, w_g), find_lambda_max_in_g_for_vsnorm(sumstat, w_p)), na.rm=TRUE))
  if (!is.finite(l_max)) stop("Unable to compute a starting value for lambda: all observed effects are zero or undefined.")
  test <- solve_at(l_max)
  n_halving <- 0
  while (test$info[dim(test$info)[1],]$z_norm==0) {
    n_halving <- n_halving + 1
    if (n_halving > max_halving) stop("No value of lambda giving a non-empty model was found after ", max_halving, " halvings.")
    l_max <- l_max / 2
    test <- solve_at(l_max)
  }
  test <- solve_at(f_test*l_max)
  if(test$info[dim(test$info)[1],]$z_norm!=0){result <- f_test * l_max}else{result <- f_selec * l_max}
  return(result)
}
```
Les valeurs par défaut des réglages ADMM sont celles de `metapsg_gv` : le
comportement est inchangé quand l'utilisateur garde les réglages par défaut.

### 1.3 Transmettre les réglages ADMM à `find_lambda_max`
Dans `f_selec_mod_metapsg_gv` (`R/metapsg_gv_solver.r`) et dans le `foreach` de
`f_metapsg_gv_parallel_grid` (`R/metapsg_parallel.r`), l'appel à
`find_lambda_max(...)` reçoit en plus
`rho=rho, t_inc=t_inc, t_dec=t_dec, mu=mu, max_iter=max_iter, tol_abs=tol_abs, tol_rel=tol_rel`.

### 1.4 Une seule valeur de α (`R/metapsg_parallel.r`)
Dans `f_metapsg_gv_parallel_grid` : `for (i in 2:n_alpha){` devient
`for (i in seq_len(n_alpha)[-1]){`.

### 1.5 `same_var`
- `metapsg_gv` : retirer l'argument `same_var = FALSE` ; dans ses trois appels à
  `f_cpp_mrho_inv`, retirer `same_var=same_var`.
- `f_cpp_mrho_inv` : retirer l'argument `same_var` et la branche `if (same_var==TRUE)`,
  en gardant seulement la boucle `for (k in 1:K){result[,,k] <- arma_inv(...)}`
  (et les commentaires associés à la branche supprimée).
- `f_metapsg_gv_parallel_grid` : retirer l'argument `same_var = FALSE`, et
  `same_var=same_var` dans son appel depuis `f_selec_mod_metapsg_gv_parallel`.
- Les deux appels à `metapsg_gv` (dans `f_selec_mod_metapsg_gv` et dans le
  `foreach` de `f_selec_mod_metapsg_gv_parallel`) : retirer `same_var = FALSE`
  et passer tous les autres arguments par nom (voir le diff de référence).
- `f_selec_mod_metapsg_gv` et `f_selec_mod_metapsg_gv_parallel` (exportées) :
  **garder** l'argument `same_var = FALSE` dans la signature, et ajouter comme
  première ligne du corps :
  `  if (isTRUE(same_var)) warning("'same_var' is not supported and is ignored.")`
- Ne pas toucher à `f_cpp_inv_var` ni aux fonctions inutilisées `f_inv_var` et
  `f_mrho_inv` (nettoyage prévu plus tard).

### 1.6 Supprimer `f_grid_log10` (`R/metapsg_sub.r`)

## 2. Documentation
Dans les blocs roxygen de `f_selec_mod_metapsg_gv` et
`f_selec_mod_metapsg_gv_parallel` :
- `same_var` : « Not supported and ignored (a warning is issued if `TRUE`). Kept
  for compatibility with existing scripts. »
- retirer la mention selon laquelle `alpha` doit contenir au moins deux valeurs ;
- si la documentation dit que la recherche de lambda_max utilise les réglages
  par défaut, indiquer qu'elle utilise désormais les réglages ADMM fournis.

Puis `Rscript -e "roxygen2::roxygenise()"`. `NAMESPACE` ne doit pas changer.

## 3. Archive du papier (`paper/`)
1. `paper/application/script_metapsg_BC_TC.r`, ligne 99 : `RESULT FILE MANAGEMENT`
   devient `# RESULT FILE MANAGEMENT`.
2. `paper/application/README.md`, ligne 210 : `append(pleio_genes,` devient
   `append(pleio_genes_res,`.
3. Ajouter à la fin de `paper/README.md` (contenu exact) :
```markdown

## Changes to the archived material

The files of this folder are kept as they were at publication (tag `v1.1`),
except for two one-line fixes so that the archived code runs:

- `application/script_metapsg_BC_TC.r`: the line `RESULT FILE MANAGEMENT` was not
  commented out, which caused a syntax error.
- `application/README.md`: `append(pleio_genes, ...)` replaced by
  `append(pleio_genes_res, ...)`; the variable `pleio_genes` was not defined.
```
Vérifier avec `git diff` que seules ces lignes changent dans ces fichiers.

## 4. `.gitignore`
Ajouter la ligne `dev/verification/*.log`.

## 5. Construction et vérification
1. `Rcmd build .`, puis avec `NOT_CRAN=true` : `Rcmd check --no-manual MPSG_1.1.9000.tar.gz`.
2. Attendu : Status OK (0 ERROR, 0 WARNING, 0 NOTE). Tests : les 48 assertions du
   test de référence et les 5 de `test-bugfixes.R` passent (53 au total). Reporter
   mot pour mot toute WARNING ou NOTE ; ne rien corriger d'autre sans demander.

## 6. Reproduction complète du tableau 6 (contrôle de non-régression)
1. `Rcmd INSTALL .`
2. En arrière-plan :
   `Rscript dev/verification/reproduce_table6.R > dev/verification/reproduce_table6.log 2>&1`
3. À la fin : `git diff dev/verification/reproduce_table6_results.txt`. Seules les
   lignes de durée (`time:`, `total time:`) peuvent différer. Si une autre ligne
   change, ne rien corriger : s'arrêter et le signaler.

## 7. Indexation
`git add` des fichiers modifiés, créés et supprimés (y compris
`reproduce_table6_results.txt` s'il a changé). Le journal `.log` est ignoré.
Pas de commit.

## 8. Rapport
`dev/rapports/04-corrections.md`, au format défini dans `CLAUDE.md`.
Message de commit proposé : écrire le message ci-dessous dans le fichier
`dev/rapports/04-commit-message.txt` (il sera utilisé avec `git commit -F`) :
```
Fix known bugs (no change to results with the settings of the paper)

- find_lambda_max: ignore infinite bounds when an observed effect is exactly 0,
  and cap the number of halvings (fixes an infinite loop)
- find_lambda_max: use the ADMM settings given by the user
- parallel grid search: allow a single value of alpha
- same_var: removed from the computation; ignored with a warning in the
  legacy functions
- f_grid_log10: removed (wrong and never used)
- paper/: two one-line fixes so that the archived scripts run
```
