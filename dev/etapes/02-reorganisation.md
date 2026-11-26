# Étape 2 — Réorganisation du dépôt (aucune modification de code)

## Objectif
Regrouper le matériel du papier dans `paper/` pour libérer la racine du dépôt,
qui accueillera le package à l'étape 3. Uniquement des déplacements (avec
`git mv`, pour conserver l'historique), un nouveau `paper/README.md` et la mise
à jour des liens du `README.md` racine.

## Pré-vérifications (s'arrêter et le signaler si l'une échoue)
1. `git status` : rien de modifié, seuls `CLAUDE.md` et `dev/` sont non suivis.
2. Le tag `v1.1` existe en local et sur origin :
   `git tag -l v1.1` puis `git ls-remote --tags origin v1.1`.
3. La branche courante est `main` et elle est à jour : `git fetch origin` puis
   `git status` doit indiquer « up to date with 'origin/main' ».
4. La racine contient `Application/`, `R/`, `README.md`, `Simulation_study/`,
   `src/` (en plus de `CLAUDE.md` et `dev/`).

## Actions
1. Créer la branche `dev` depuis `main` : `git switch -c dev`.
2. Créer le dossier `paper/`, puis :
   - `git mv Application paper/application`
   - `git mv Simulation_study paper/simulation_study`
3. Créer `paper/README.md` avec exactement le contenu de la section suivante.
4. Dans le `README.md` racine, mettre à jour uniquement les liens vers les
   anciens emplacements (dont l'URL absolue
   `https://github.com/PESugier/MPSG/blob/main/Application/README.md`) pour
   pointer vers les nouveaux, en liens relatifs (`paper/application/README.md`).
   Ne rien changer d'autre dans ce fichier.
5. Rechercher dans tout le dépôt les autres mentions de `Application/` ou de
   `Simulation_study/` et les lister dans le rapport, sans les modifier (les
   scripts du papier contiennent des chemins absolus du cluster : on n'y touche pas).
6. `git add CLAUDE.md dev/ paper/README.md README.md`

## Contenu de `paper/README.md`

````markdown
# Material of Sugier et al. (2025)

This folder contains the scripts, data and results of the paper

> Sugier P-E, Asgari Y, Sedki M, Truong T, Liquet B (2025). Meta-analysis models
> with group structure for pleiotropy detection at gene and variant level using
> summary statistics from multiple datasets. *Biostatistics*, 26(1), kxaf037.
> <https://doi.org/10.1093/biostatistics/kxaf037>

The state of the repository at the time of publication is available with the
git tag `v1.1`.

| Folder | Content |
|---|---|
| [`application/`](application/) | Summary statistics, annotation, analysis scripts and results of the breast and thyroid cancer application |
| [`simulation_study/`](simulation_study/) | Scripts of the simulation study and resulting tables |

Note on the application data: the summary statistics cover 3,707 SNPs in 331
genes. The 3,766 SNPs reported in the paper include duplicated SNPs (in complete
linkage disequilibrium) that were removed before the analysis.
````

## À ne pas faire
- Ne modifier aucun fichier de `R/`, `src/`, ni le contenu des fichiers déplacés.
- Ne pas supprimer les binaires de `src/` (étape 3).
- Ne pas corriger les bugs connus des scripts (étape 4).
- Pas de commit, pas de push.

## Vérifications finales
1. `git status` : les déplacements apparaissent comme des renommages
   (`renamed:`), pas comme suppression + ajout.
2. `git diff --cached --stat -M` et `git diff --cached -M --name-status` :
   tous les fichiers déplacés doivent être en `R100` (renommage sans modification).
3. `git diff --cached README.md` : seules les lignes de liens ont changé, et les
   fins de ligne sont inchangées.

## Rapport
`dev/rapports/02-reorganisation.md`, au format défini dans `CLAUDE.md`.
Message de commit proposé : `Move paper material (application, simulation study) to paper/`
