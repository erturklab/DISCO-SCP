# 3D-SCP: 3D single-cell proteomics of cleared mouse brain

Analysis code for a study of 3D single-cell proteomics of DISCO-cleared mouse brain (Ertürk lab, Helmholtz Zentrum München).

The code reproduces all quantitative panels of Figures 1–4 and Extended Data Figures 1–7: linking 3D morphology to single-cell proteomes, method performance, cell-type identity (Thy1-EGFP neurons vs Cx3cr1-EGFP microglia) and regional differences. Schematics, 3D registration and tracing images are not part of this repository.

## Contents

| Script | Content | Figures |
|---|---|---|
| `Scripts/00_build_master_dataset.R` | DIA-NN output → master protein table, sample metadata, contaminant flags | – |
| `Scripts/01_morphology_link.R` | Morphology of traced neurons vs depth-adjusted protein levels | Fig. 1, ED Figs. 1–2 |
| `Scripts/02_method_performance.R` | LC columns, input titration, sensitivity, compartments, batch concordance | Fig. 2, ED Fig. 3 |
| `Scripts/03_cell_type_identity.R` | Neurons vs microglia: detection, differential abundance, literature overlap, controls | Fig. 3, ED Figs. 4–5 |
| `Scripts/04_regional_analysis.R` | Regional differences per cell type: differential expression and GSEA (GO:BP) | Fig. 4, ED Figs. 6–7 |
| `Scripts/05_unpack_supplementary_data2.R` | Optional: writes the sheets of Supplementary Data 2 back into the files the scripts expect (replaces step 00) | – |
| `Scripts/shared/setup.R` | Paths, thresholds, colours, data loading (sourced by scripts 01–04) | – |
| `Scripts/shared/excluded_samples.R` | Samples removed from the whole study (duplicate neurons, invalid well) | – |
| `Scripts/shared/functions.R` | Shared functions (depth adjustment, safe Spearman, themes, source-data export) | – |
| `Scripts/compare_summary_numbers.R` | Optional: compares the summary-number files of two runs | – |
| `lmd_export/` | Python notebook that converts QuPath annotations into the Leica LMD cutting file (8-well strip); adapted from CosciaLab/Qupath_to_LMD | Methods (sample collection) |
| `install_packages.R` | Installs the required R packages | – |

Each figure script runs on its own once the master dataset exists (step 2 below).

## Requirements

* R 4.4.2 (Windows)
* R packages: run `source("install_packages.R")` once. The package versions of every run are saved in `Figure*_software_versions.txt`.
* Python (only for `lmd_export/`): see `lmd_export/README.md`.
* Font: figures are set in Arial. On Linux, install Arial or a metric-compatible font (e.g. Liberation Sans).
* Internet access for `02_method_performance.R` (UniProt, see Reproducibility notes).
* Tested on Windows with 32 GB RAM and an Intel Core i7-12800HX (12th Gen) processor. Scripts 00–04 together take about 6 minutes.

## Data

The data are not stored in this repository.

| Data | Where |
|---|---|
| Raw mass spectrometry files and DIA-NN output | PRIDE / ProteomeXchange, accession will be added here (upload in progress) |
| Processed master dataset, control-well screen and morphology table | Supplementary Data 2 of the paper (one Excel file, see "Two ways to reproduce") |

Published tables used for comparison (download from the publishers; not redistributed here). Put them in `Data/` under these names:

| File name | Source |
|---|---|
| `Sharma2015_SuppTable2.xlsx` | Sharma, K. et al. Cell type– and brain region–resolved mouse brain proteome. *Nat. Neurosci.* **18**, 1819–1831 (2015). Supplementary Table 2. https://doi.org/10.1038/nn.4160 |
| `Rayaprolu2022_SuppData9.xlsx` | Rayaprolu, S. et al. Cell type-specific biotin labeling in vivo resolves regional neuronal and astrocyte proteomic differences in mouse brain. *Nat. Commun.* **13**, 2927 (2022). Supplementary Data 9. https://doi.org/10.1038/s41467-022-30623-x |
| `Johnson2026_SI_Tables.xlsx` | Johnson, C. C. et al. Aspiration Patch Proteomics Enables CE–ESI–MS Proteotyping of Identified Single Neurons in Intact Brain Tissue. *Anal. Chem.* **98**, 28088–28096 (2026). Table S3. https://doi.org/10.1021/acs.analchem.6c03049 |
| `Lloyd2024_SuppTable_mmc2.xlsx` | Lloyd, A. F. et al. Deep proteomic analysis of microglia reveals fundamental biological differences between model systems. *Cell Rep.* **43**, 114908 (2024). Table S2 (sheet "Mouse"). https://doi.org/10.1016/j.celrep.2024.114908 |
| `1-s2.0-S0896627326003855-mmc4.xlsx` | Zhang, H. et al. Single-cell proteome atlas of aging mouse microglia reveals subpopulation-specific phagoproteome. *Neuron* https://doi.org/10.1016/j.neuron.2026.05.014 (2026). Table S3. |

Further inputs expected in `Data/`: for step 00 the DIA-NN protein matrix (`report_CombinedBrainEmpLib_05Jun2026_360_MBR.pg_matrix.tsv`) and the plate design table (`Plate_Design_Master_Table.xlsx`); for the figure scripts the HeLa test-injection statistics (`*200fL_200pL_QCs_combiLib_17Jun2026*stats*.tsv`), the DIA-NN report (`report_CombinedBrainEmpLib_05Jun2026_360_MBR.parquet`) and precursor matrix (`..._pr_matrix.tsv`). File names are set in `Scripts/shared/setup.R` and at the top of each script.

## How to run

1. Put the input files in `Data/`.
2. Build the master dataset (once): `source("Scripts/00_build_master_dataset.R")`.
3. Set the R working directory to the project folder (or `options(scp.project_root = "/path/to/project")`).
4. Run a figure script, e.g. `source("Scripts/04_regional_analysis.R")`.

### Two ways to reproduce

* **Option A, from the raw data:** download the PRIDE files into `Data/` and follow steps 1–4 above (script 00 builds the master tables).
* **Option B, from Supplementary Data 2:** put `Supplementary_Data_2.xlsx` into `Data/` and run `source("Scripts/05_unpack_supplementary_data2.R")`. It writes the four tables into `Results/master_dataset/` and `Morphological_analysis.xlsx` into `Data/`. Then skip step 2 and run the figure scripts. The precursor matrix, the DIA-NN report and the HeLa test-injection statistics are still needed from PRIDE.

Outputs are written to `Results/figures/<Figure>/`: figure panels (SVG), the tables behind every panel (CSV), `Source_Data_Figure*.xlsx`, `Figure*_summary_numbers.txt` (all numbers quoted in legends and text) and `Figure*_software_versions.txt`.

## Analysis summary

* **Cells**: contaminant-filtered protein groups per cell; ≥ 200 = usable cell, ≥ 600 = analysis set for tests and regional comparisons.
* **Depth adjustment**: per protein, log2 intensity ~ log2(protein groups per cell) [+ batch]; the residual is used for display and correlations. No normalisation, no imputation.
* **Statistics**: all tests are reported as nominal P and BH-adjusted q; morphology associations and pathway results are exploratory.

## Reproducibility notes

* Random seeds are set where randomness is used (`set.seed(42)` for fgsea, fixed jitter seeds).
* GO gene sets come from `msigdbr`; the version is written to `Figure4_summary_numbers.txt`.
* Subcellular annotation (Fig. 2d) is queried from UniProt and cached in `Data/uniprot_subcellular_location_cache.csv`. Keep the cache to get identical numbers; refreshing it uses the current UniProt release.
* Package versions of every run are saved in `Figure*_software_versions.txt`.

## License

R code: Apache License 2.0 (see `LICENSE`). `lmd_export/`: GPL-3.0 (see `lmd_export/LICENSE`), because it is adapted from GPL-3.0 code. Data deposits carry their own licenses.

## Contact

Eva Küsters, Ertürk lab, Helmholtz Zentrum München — questions and bug reports: please open an issue on GitHub.
