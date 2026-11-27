# Rapport — Étape 4 : corrections des bugs connus

Date : 2026-11-27 — branche `dev` — aucun commit effectué.

## 1. Résumé

Les bugs connus de `dev/DECISIONS.md` sont corrigés (décisions 7 et 10) :
- boucle infinie de `find_lambda_max` quand un β̂ vaut 0 ;
- réglages ADMM de l'utilisateur ignorés par `find_lambda_max` ;
- échec de la grille parallèle avec un seul α ;
- `same_var`, retiré du calcul et ignoré avec un avertissement dans les fonctions héritées ;
- `f_grid_log10`, supprimée ;
- deux erreurs d'une ligne dans les scripts archivés de `paper/`.

Les modifications de `R/` sont identiques, ligne à ligne, à celles du diff de référence `dev/etapes/04-reference.diff`.
`Rcmd check` se termine en **Status: OK**, avec 53 assertions réussies sur 53.
La reproduction du tableau 6 redonne exactement les mêmes résultats qu'à l'étape 3 : seules les durées changent.

## 2. Fichiers créés, modifiés ou déplacés

`git status --short` (après indexation) :

```
M  .gitignore
M  CLAUDE.md
M  R/metapsg_gv_solver.r
M  R/metapsg_parallel.r
M  R/metapsg_sub.r
M  dev/DECISIONS.md
A  dev/etapes/04-corrections.md
A  dev/etapes/04-reference.diff
A  dev/rapports/04-commit-message.txt
A  dev/rapports/04-corrections.md
M  dev/verification/reproduce_table6_results.txt
M  man/f_selec_mod_metapsg_gv.Rd
M  man/f_selec_mod_metapsg_gv_parallel.Rd
M  paper/README.md
M  paper/application/README.md
M  paper/application/script_metapsg_BC_TC.r
A  tests/testthat/test-bugfixes.R

```

Ignorés, non indexés : `MPSG.Rcheck/`, `MPSG_1.1.9000.tar.gz`, `dev/verification/reproduce_table6.log`, `src/MPSG.dll`, `src/*.o`.

`git diff --cached --stat -M` :

```
 .gitignore                                    |   1 +
 CLAUDE.md                                     |   9 +-
 R/metapsg_gv_solver.r                         |  17 ++--
 R/metapsg_parallel.r                          |  14 +--
 R/metapsg_sub.r                               |  52 ++++++-----
 dev/DECISIONS.md                              |  18 ++--
 dev/etapes/04-corrections.md                  | 176 ++++++++++++++++++++++++++++++++++++
 dev/etapes/04-reference.diff                  | 221 ++++++++++++++++++++++++++++++++++++++++++++++
 dev/rapports/04-commit-message.txt            |  10 +++
 dev/rapports/04-corrections.md                | 134 ++++++++++++++++++++++++++++
 dev/verification/reproduce_table6_results.txt |   6 +-
 man/f_selec_mod_metapsg_gv.Rd                 |   4 +-
 man/f_selec_mod_metapsg_gv_parallel.Rd        |   7 +-
 paper/README.md                               |  10 +++
 paper/application/README.md                   |   2 +-
 paper/application/script_metapsg_BC_TC.r      |   2 +-
 tests/testthat/test-bugfixes.R                |  34 +++++++
 17 files changed, 659 insertions(+), 58 deletions(-)

```

## 3. Vérifications effectuées

### Section 0 — Pré-vérifications : OK
- Branche `dev`, HEAD = `c6b36df` (« Make the repository an installable R package… ») = `origin/dev`.
- `git status` : `CLAUDE.md` et `dev/DECISIONS.md` modifiés. Nouveaux fichiers : `dev/etapes/04-corrections.md`, `dev/etapes/04-reference.diff` et `tests/testthat/test-bugfixes.R`. Rien d'autre.

### Section 1 — Corrections dans `R/`
- Méthode : remplacements littéraux par un script Perl. Chaque remplacement contrôle le nombre d'occurrences attendu (par exemple 3 appels à `f_cpp_mrho_inv`, 2 `max(v_this_g, na.rm=TRUE)`). Ce choix préserve le CRLF et les espaces de fin de ligne.
- Contrôle : j'ai extrait les lignes de code ajoutées et supprimées de `git diff R/` (hors lignes roxygen et lignes vides) et je les ai comparées à celles de `dev/etapes/04-reference.diff`. Les deux sont **identiques** (75 lignes).
- Les arguments sont passés par nom dans tous les appels modifiés.
- Les trois fichiers restent en CRLF pur (vérifié avec `file`).

### Section 2 — Documentation
- `same_var` est maintenant documenté ainsi : « Not supported and ignored (a warning is issued if `TRUE`). Kept for compatibility with existing scripts. » La fonction parallèle hérite de ce texte via `@inheritParams`.
- La mention « it must contain at least two values » est retirée de `alpha` (version parallèle).
- `roxygenise()` n'a réécrit que `man/f_selec_mod_metapsg_gv.Rd` et `man/f_selec_mod_metapsg_gv_parallel.Rd`. `NAMESPACE` et `DESCRIPTION` sont inchangés.

### Section 3 — Archive `paper/`
- `git diff` ne montre que les 3 modifications demandées : la ligne 99 du script, la ligne 210 de `application/README.md`, et le bloc ajouté à la fin de `paper/README.md`. Les fins de ligne d'origine sont conservées : CRLF dans le dossier `application/`, LF pour `paper/README.md`.
- Contrôle supplémentaire : `parse("paper/application/script_metapsg_BC_TC.r")` réussit désormais.

### Section 4 — `.gitignore`
- La ligne `dev/verification/*.log` est ajoutée. `git check-ignore` confirme que `reproduce_table6.log` est ignoré.

### Section 5 — Build et check
- `Rcmd build .` : OK.
- `NOT_CRAN=true Rcmd check --no-manual MPSG_1.1.9000.tar.gz` : **Status: OK**, aucune WARNING ni NOTE (0 occurrence dans `00check.log`).
- Tests : `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 53 ]` (48 assertions de référence + 5 de `test-bugfixes.R`).

### Section 6 — Reproduction du tableau 6
- `Rcmd INSTALL .` : OK. Le script a tourné en arrière-plan sur 11 cœurs, en 20 min 24 s, avec le code de sortie 0.
- `git diff dev/verification/reproduce_table6_results.txt` : **seules les lignes de durée diffèrent** :
  - poids adaptatifs : `time: 16.1 min` devient `14.7 min` ;
  - poids standard : `time: 5.4 min` devient `5.7 min` ;
  - total : `total time: 21.5 min` devient `20.4 min`.
- Toutes les autres lignes sont identiques : meilleurs modèles (α, λ, MIC), gènes, nombres de SNPs, et écarts maximaux aux β de référence (1.18e-14 et 7.91e-16).

## 4. Écarts, anomalies, questions ouvertes

1. **`f_cpp_mrho_inv` : commentaires supprimés.** J'ai supprimé les deux lignes commentées `#result <- apply(...)` de la branche retirée, comme dans le diff de référence. La ligne `#result <- do.call(abind::abind, ...)`, qui suit le bloc, est conservée.
2. **Mentions de `same_var` restantes** :
   - dans les signatures de `f_cpp_inv_var`, `f_inv_var` et `f_mrho_inv`, comme demandé (ces fonctions ne sont pas touchées) ;
   - dans le bloc commenté `# TEST` à la fin de `R/metapsg_parallel.r` (2 lignes), hors périmètre. Il faudra le nettoyer lors d'une prochaine étape.
3. **Documentation de lambda_max.** La documentation ne disait pas que la recherche de `lambda_max` utilise les réglages par défaut. Rien n'a donc été ajouté, conformément à la section 2. Faut-il préciser maintenant que les réglages ADMM fournis sont utilisés ?
4. **Lignes vides.** La suppression de `f_grid_log10` laisse deux lignes vides consécutives dans `R/metapsg_sub.r`, comme dans le diff de référence.
5. **Fins de ligne des `.Rd`.** Les `.Rd` régénérés sont en LF dans la copie de travail. Git avertit « LF will be replaced by CRLF » à cause de `core.autocrlf=true`. Sans effet sur le contenu indexé ; à traiter avec le `.gitattributes` prévu.
6. **Artefacts de build.** `src/*.o` et `src/MPSG.dll` sont présents après `Rcmd INSTALL .`. Ils sont ignorés par git. J'avais supprimé ceux produits par `roxygenise()` avant le build.
7. Détail technique : un `Rscript -e` contenant des guillemets imbriqués a échoué depuis Git Bash. Je suis passé par un script temporaire hors du dépôt, comme le recommande `CLAUDE.md`.

## 5. Message de commit proposé

Écrit dans `dev/rapports/04-commit-message.txt` (pour `git commit -F`) :

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
