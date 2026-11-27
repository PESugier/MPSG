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

## Changes to the archived material

The files of this folder are kept as they were at publication (tag `v1.1`),
except for two one-line fixes so that the archived code runs:

- `application/script_metapsg_BC_TC.r`: the line `RESULT FILE MANAGEMENT` was not
  commented out, which caused a syntax error.
- `application/README.md`: `append(pleio_genes, ...)` replaced by
  `append(pleio_genes_res, ...)`; the variable `pleio_genes` was not defined.
