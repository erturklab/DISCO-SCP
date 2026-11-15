###############################################################################
# setup.R -- shared setup for 01_morphology_link.R .. 04_regional_analysis.R
#
# What this file does (in order):
#   1. Config: paths, thresholds, region names, samples excluded from the study;
#      sources functions.R (all function definitions)
#   2. Load the master protein table, rename its columns, drop excluded samples,
#      remove contaminants, build the log2 protein matrix
#   3. Optional precursor-level matrix (for the input-amount panels)
#   4. Colours and figure-size constants
#   5. Cell cohorts (cohort_meta) and per-cell protein counts
#   6. Output folders
# Functions (depth_residuals(), themes, input checks, ...) live in functions.R.
#
# Sourced at the top of each figure script; not run on its own.
#
# Run all scripts with the R working directory set to the project folder:
#   Data/      raw instrument files, literature reference tables
#   Scripts/   0N_*.R, 00_build_master_dataset.R;  Scripts/shared/  setup.R, functions.R
#   Results/   master_dataset/ (from 00_build_master_dataset.R) and figures/
# (Alternatively set options(scp.project_root = "/path/to/project") before sourcing.)
# Required R packages: see install_packages.R (run it once) and README.md.
###############################################################################

# ---- 1. Config ---------------------------------------------------------------

PROJECT_ROOT <- getOption("scp.project_root", default = normalizePath(getwd(), winslash = "/"))
if (!dir.exists(file.path(PROJECT_ROOT, "Data"))) {
  stop("No 'Data' folder found in ", PROJECT_ROOT, ". Set the R working directory to the project folder ",
       "(or options(scp.project_root = ...)) and try again.")
}


# All function definitions (input checks, depth_residuals(), themes, ...)
options(width = 250)   
source(file.path(PROJECT_ROOT, "Scripts", "shared", "functions.R"))

cran_pkgs <- c("tidyverse", "readxl", "ggrepel", "patchwork", "colorspace", "ggdist", "scales")
check_packages(cran_pkgs, "setup.R")
invisible(lapply(setdiff(cran_pkgs, "scales"), library, character.only = TRUE))   

# Figures are written as SVG (SVG_DEVICE in functions.R); tables as CSV.

# Literature reference tables (supplementary tables of published studies)
SHARMA_FILE  <- file.path(PROJECT_ROOT, "Data", "Sharma2015_SuppTable2.xlsx")
SHARMA_SHEET <- "For Suppl table 2 proteinGroups"
SHARMA_SKIP  <- 2

LLOYD_FILE         <- file.path(PROJECT_ROOT, "Data", "Lloyd2024_SuppTable_mmc2.xlsx")
JOHNSON_FILE       <- file.path(PROJECT_ROOT, "Data", "Johnson2026_SI_Tables.xlsx")     
RAYAPROLU2022_FILE <- file.path(PROJECT_ROOT, "Data", "Rayaprolu2022_SuppData9.xlsx")    
RAYAPROLU2022_SHEET <- "SuppData9"
RAYAPROLU2022_SKIP  <- 3                                                                   # title rows; header = row 4
NEURON2026_FILE    <- file.path(PROJECT_ROOT, "Data", "1-s2.0-S0896627326003855-mmc4.xlsx")  
PRECURSOR_PARQUET_FILE <- file.path(PROJECT_ROOT, "Data", "report_CombinedBrainEmpLib_05Jun2026_360_MBR.parquet")  # DIA-NN long-format report

OUT_DIR <- file.path(PROJECT_ROOT, "Results", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Two different cell-inclusion bars -- do not mix them up:
#   MIN_PROTEINS_VIABLE (200): "is this a usable cell at all". Applied once, centrally
#     (passes_qc_tier below), on the contaminant-filtered protein count.
#   MIN_PROTEINS (600): stricter bar for analyses that need depth per cell
#     (differential expression, regional comparisons, UpSet plots).
MIN_PROTEINS_VIABLE <- 200
MIN_PROTEINS        <- 600
MISS_THRESHOLD      <- 0.20   # minimum fraction of cells detecting a protein for the cell-type DE (03_cell_type_identity.R)
MIN_N_DEPTH_FIT     <- 15     # minimum number of cells detecting a protein to fit its depth regression

# Region names (single source of truth: rename a region here and every script follows).
REGIONS_THY1 <- c("Cingulate Cortex", "Somatosensory Cortex", "Hippocampal Region")
REGIONS_MG   <- c("Cerebellum Gray Matter", "Cerebellum White Matter", "Somatosensory Cortex")

# Short labels for axes and legends. The region is still called "Cingulate Cortex" in the
# data; in the manuscript it is named correctly - medial prefrontal cortex (mPFC). Only the display
# label changes, the data value stays.
REGION_ABBREV <- c(
  "Cerebellum Gray Matter"  = "CereGM",
  "Cerebellum White Matter" = "CereWM",
  "Somatosensory Cortex"    = "SoC",
  "Cingulate Cortex"        = "mPFC",
  "Hippocampal Region"      = "Hippo"
)
# Full written-out display names (legends, methods)
REGION_FULL_DISPLAY <- c(
  "Cerebellum Gray Matter"  = "Cerebellum gray matter",
  "Cerebellum White Matter" = "Cerebellum white matter",
  "Somatosensory Cortex"    = "Somatosensory cortex",
  "Cingulate Cortex"        = "Medial prefrontal cortex",
  "Hippocampal Region"      = "Hippocampal formation"
)

# Samples excluded from the whole study. They are removed from every table right after loading,
# so no analysis can pick them up by accident.
source(file.path(PROJECT_ROOT, "Scripts", "shared", "excluded_samples.R"))   # DUPLICATE_NEURON_SAMPLES, INVALID_WELL_SAMPLES, EXCLUDED_SAMPLES

# ---- 2. Load the master dataset ----------------------------------------------

message("\nLoading master dataset...")
MASTER_DATASET_DIR <- file.path(PROJECT_ROOT, "Results", "master_dataset")
req_files <- file.path(MASTER_DATASET_DIR,
                       c("3D_SCP_all_Data_Master_file.csv", "master_sample_metadata.csv", "contamination_analysis.csv"))
missing_files <- req_files[!file.exists(req_files)]
if (length(missing_files) > 0) {
  stop("Master dataset files not found:\n  ", paste(missing_files, collapse = "\n  "),
       "\nRun 00_build_master_dataset.R first.", call. = FALSE)
}

master_data              <- read_csv(req_files[1], show_col_types = FALSE)   # proteins x samples (+ annotation columns)
sample_metadata          <- read_csv(req_files[2], show_col_types = FALSE)   # one row per sample
contamination_annotation <- read_csv(req_files[3], show_col_types = FALSE)   # one row per gene, contaminant flags

# Required columns and unique identifiers (nothing is dropped or repaired silently)
require_columns(master_data, "gene_label", "Master protein table")
require_columns(sample_metadata, c("sample_id", "clean_sample_name", "cell_type", "tier", "region", "batch",
                                   "cell_num", "raw_filename", "n_proteins_detected"), "master_sample_metadata.csv")
require_columns(contamination_annotation, c("gene_label", "is_crap_or_contaminant"), "contamination_analysis.csv")
if (anyDuplicated(sample_metadata$sample_id) || anyDuplicated(sample_metadata$clean_sample_name)) {
  stop("master_sample_metadata.csv: sample_id and clean_sample_name must each be unique.", call. = FALSE)
}

# The master table names its sample columns with the original raw sample_id. Rename them to
# clean_sample_name, which is what the rest of the project uses.
annotation_cols <- c("gene_label", "Protein.Group", "Protein.Names",
                     "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences")
sample_cols <- setdiff(colnames(master_data), annotation_cols)                       # all columns that are samples
clean_name  <- setNames(sample_metadata$clean_sample_name, sample_metadata$sample_id) # raw name -> clean name
if (anyNA(clean_name[sample_cols])) {
  stop("Sample column(s) without an entry in sample_metadata: ",
       paste(sample_cols[is.na(clean_name[sample_cols])], collapse = ", "))
}
colnames(master_data)[match(sample_cols, colnames(master_data))] <- unname(clean_name[sample_cols])
not_measured <- setdiff(sample_metadata$clean_sample_name, colnames(master_data))
if (length(not_measured) > 0) {
  warning(length(not_measured), " sample(s) in the metadata have no column in the master table (no proteome): ",
          paste(head(not_measured, 5), collapse = ", "), if (length(not_measured) > 5) ", ...")
}
unknown_excluded <- setdiff(EXCLUDED_SAMPLES, sample_metadata$clean_sample_name)
if (length(unknown_excluded) > 0) {
  warning("EXCLUDED_SAMPLES not found in the metadata (check the names): ", paste(unknown_excluded, collapse = ", "))
}

# Remove excluded samples everywhere, once.
master_data     <- master_data %>% dplyr::select(-dplyr::any_of(EXCLUDED_SAMPLES))
sample_metadata <- sample_metadata %>% dplyr::filter(!clean_sample_name %in% EXCLUDED_SAMPLES)
message("  Excluded ", length(EXCLUDED_SAMPLES), " samples (duplicate neurons / invalid well); ",
        nrow(sample_metadata), " samples remain")

# Remove contaminants (flags come from 00_build_master_dataset.R: cRAP, skin, immunoglobulin,
# and the tissue-free PEN control screen).
clean_genes <- contamination_annotation %>% dplyr::filter(!is_crap_or_contaminant) %>% dplyr::pull(gene_label)
protein_df_filtered <- master_data %>%
  dplyr::select(-dplyr::any_of(setdiff(annotation_cols, "gene_label"))) %>%
  dplyr::filter(gene_label %in% clean_genes)
if (anyDuplicated(protein_df_filtered$gene_label)) {
  warning("Duplicated gene_label values in the master table -- name-based lookups use the first match.")
}

# Export for manual inspection (not read back by any script).
write_csv(protein_df_filtered, file.path(PROJECT_ROOT, "Results", "protein_data_contamination_filtered.csv"))

# Protein x sample matrices. log2 of the MaxLFQ intensity; NA = not detected. No imputation.
protein_mat_linear <- as.matrix(protein_df_filtered %>% dplyr::select(-gene_label))
rownames(protein_mat_linear) <- protein_df_filtered$gene_label
n_zero <- sum(protein_mat_linear == 0, na.rm = TRUE)
if (n_zero > 0) {   # a zero would become -Inf after log2; it means "not detected"
  message("  Note: ", n_zero, " zero intensities set to NA")
  protein_mat_linear[protein_mat_linear == 0] <- NA
}
if (!is.numeric(protein_mat_linear)) stop("The master protein table contains non-numeric intensity columns.", call. = FALSE)
if (any(protein_mat_linear < 0, na.rm = TRUE) || any(is.infinite(protein_mat_linear)) || any(is.nan(protein_mat_linear))) {
  stop("The master protein table contains negative, infinite or NaN intensities.", call. = FALSE)
}
protein_mat_log2 <- log2(protein_mat_linear)

# Number of proteins detected per sample (after contaminant removal). This is THE depth measure:
# used for cell inclusion and as depth covariate. Use this, not the raw count in sample_metadata.
n_proteins_lookup <- setNames(colSums(!is.na(protein_mat_log2)), colnames(protein_mat_log2))

# ---- 3. Precursor-level matrix (not MaxLFQ) -----------------------------------
# Used only for the input-amount panels (1x vs 10x); see load_precursor_matrix() in functions.R.
precursor_mat_log2 <- load_precursor_matrix(PRECURSOR_PARQUET_FILE)

# ---- 4. Colours and figure-size constants ------------------------------------
PAL_THY1      <- "#4A6C8C"
PAL_MICROGLIA <- "#C15A42"
PAL_BULK      <- "#D9932E"
PAL_SHARMA    <- "#7A8C4A"
PAL_ACCENT    <- "#8C6C7B"
PAL_DANGER    <- "#8C3D2E"
PAL_CAUTION   <- "#A69783"
PAL_BG        <- "#D4CDC0"
PAL_BG_LINE   <- "grey70"
PAL_CELL10    <- "#C15A42"
PAL_PRODA     <- "#5A7A5A"
PAL_MSQROB2   <- "#6B7A9E"
PAL_SUBCELL <- c(Mitochondria = "#B5563B", Nuclear = "#6B7A8C", Cytoskeleton = "#8C6C7B",
                 Membrane = "#D9932E", Cytosolic = "#7A8C4A", `Un-assigned` = "#A69783")

# Region colours, keyed by the data value ("Cingulate Cortex" = mPFC, "Hippocampal Region" = Hippo)
PAL_REGION_THY1 <- c(`Cingulate Cortex` = "#7BA7C7", `Somatosensory Cortex` = "#2F4A6B", `Hippocampal Region` = "#9B6B9E")
PAL_REGION_MG   <- c(`Somatosensory Cortex` = "#C14248", `Cerebellum Gray Matter` = "#C15A42", `Cerebellum White Matter` = "#C18042")

# Batch display labels and colours (for plots and legends)
BATCH_DISPLAY_LABELS <- c(`01Apr2026` = "Thy1-B1", `24Apr2026` = "Thy1-B2",
                          `20260529` = "Thy1-B3", `30Apr2026` = "MG-B1")
BATCH_PALETTE <- c(`01Apr2026` = colorspace::lighten(PAL_THY1, 0.35),
                   `24Apr2026` = PAL_THY1,
                   `20260529`  = colorspace::darken(PAL_THY1, 0.35),
                   `30Apr2026` = PAL_MICROGLIA)


# Journal specification: 180 mm wide figures, Arial, 7 pt text, nothing below 6 pt,
# bold axis titles. ggsave() takes inches, so sizes are written in mm and converted.
FIG_W_MM <- 180
# One rounding rule for every number printed in the summary files / text: .5 always rounds UP (R's round() and
# sprintf() round exact .5 to the even number, so 1042.5 could print as 1042 in one place and 1043 in another).
round_half_up <- function(x, digits = 0) floor(x * 10^digits + 0.5 + 1e-9) / 10^digits
PT_BASE  <- 7                                   # all text
PT_MIN   <- 6                                   # smallest allowed text (gene names, ticks)

# Cell-type labels used on plots (figures 3 and 4)
LAB_NEURON    <- "Thy1-EGFP neurons"
LAB_MICROGLIA <- "CX3CR1-EGFP microglia"

# ---- 5. Cell cohorts ----------------------------------------------------------

# cohort_meta: one row per sample (cell type, region, batch, tier, ...), plus
#   n_proteins_filtered = proteins detected after contaminant removal
#   passes_qc_tier      = n_proteins_filtered >= MIN_PROTEINS_VIABLE
# n_proteins_detected (from sample_metadata) is the RAW count before contaminant removal;
# use it only when raw yield itself is shown, never for cell inclusion or depth correction.
cohort_meta <- sample_metadata %>%
  dplyr::select(-sample_id) %>%
  dplyr::rename(sample_col = clean_sample_name, cell_type_label = cell_type) %>%
  dplyr::mutate(
    n_proteins_filtered = as.numeric(n_proteins_lookup[sample_col]),
    passes_qc_tier = !is.na(n_proteins_filtered) & n_proteins_filtered >= MIN_PROTEINS_VIABLE
  )
n_proteins_lookup_RAW_YIELD_ONLY <- setNames(cohort_meta$n_proteins_detected, cohort_meta$sample_col)

mg_singlecell_cols <- cohort_meta %>%
  dplyr::filter(tier == "tier_1x_single_cell", cell_type_label == "Microglia", passes_qc_tier) %>% dplyr::pull(sample_col)
thy1_bulk10_cols <- cohort_meta %>%
  dplyr::filter(tier == "tier_10x_area", passes_qc_tier) %>% dplyr::pull(sample_col)
thy1_singlecell_cols <- cohort_meta %>%
  dplyr::filter(tier == "tier_1x_single_cell", cell_type_label == "Thy1", passes_qc_tier) %>% dplyr::pull(sample_col)

if (length(thy1_singlecell_cols) == 0 || length(mg_singlecell_cols) == 0 || length(thy1_bulk10_cols) == 0) {
  stop("An analysis cohort is empty (Thy1 single cells: ", length(thy1_singlecell_cols), ", microglia single cells: ",
       length(mg_singlecell_cols), ", Thy1 10x areas: ", length(thy1_bulk10_cols),
       "). Check tier / cell_type values in master_sample_metadata.csv.", call. = FALSE)
}
message("  Microglia single cells: ", length(mg_singlecell_cols),
        "  |  Thy1 10x areas: ", length(thy1_bulk10_cols),
        "  |  Thy1 single cells: ", length(thy1_singlecell_cols))

# Tissue-free PEN_SC control wells (03_cell_type_identity.R, ED Fig. 4b-d). The same well is left out here as in the
# contaminant screen of 00_build_master_dataset.R (statistically confirmed outlier well; see its
# run_pen_contamination_screen()), so that the wells shown in the figure are the wells behind the
# contaminant flags. The two thresholds below must equal PEN_DETECTION_FRAC_MIN / PEN_PERCENTILE_MIN
# in 00_build_master_dataset.R; 03_cell_type_identity.R checks this against pen_sc_systematic_screen_ALL_PROTEINS.csv.
PEN_EXCLUDED_WELL_PATTERN <- "G4_474"          # matched against the raw file name
PEN_DETECTION_FRAC_MIN    <- 5 / 7             # a protein is flagged if seen in >= 5 of 7 control wells ...
PEN_PERCENTILE_MIN        <- 69                # ... at >= the 69th intensity percentile (dataset-wide)
pen_sc_cols <- cohort_meta %>%
  dplyr::filter(is_pen_sc, !stringr::str_detect(raw_filename, PEN_EXCLUDED_WELL_PATTERN)) %>%
  dplyr::pull(sample_col)
message("  PEN_SC control wells used: ", length(pen_sc_cols), " (outlier well '", PEN_EXCLUDED_WELL_PATTERN, "' left out)")

# Link between the precursor-matrix column names (parquet `Run`) and the clean sample names.
# The column keeps the name `sample_id`, as the figure scripts expect.
if (!is.null(precursor_mat_log2)) {
  sample_id_lookup <- sample_metadata %>%
    dplyr::mutate(sample_id = stringr::str_remove(raw_filename, "\\.d$")) %>%
    dplyr::select(sample_id, clean_sample_name) %>%
    dplyr::left_join(cohort_meta, by = c("clean_sample_name" = "sample_col"))
  message("  Precursor matrix: ", sum(sample_id_lookup$sample_id %in% colnames(precursor_mat_log2)),
          " of ", nrow(sample_id_lookup), " samples matched")
}

# ---- 6. Output folders -------------------------------------------------------------

# Figure numbering follows the manuscript: Fig 1 morphology link (01_morphology_link.R), Fig 2 method
# development (02_method_performance.R), Fig 3 cell identity (03_cell_type_identity.R), Fig 4 regional differences (04_regional_analysis.R).
DIR_F1          <- file.path(OUT_DIR, "Figure1_Morphology_link")
DIR_F1_SUPP     <- file.path(OUT_DIR, "Figure1_SUPPLEMENTARY")
DIR_F2          <- file.path(OUT_DIR, "Figure2_Method_performance")
DIR_F2_SUPP     <- file.path(OUT_DIR, "Figure2_SUPPLEMENTARY")
DIR_F3          <- file.path(OUT_DIR, "Figure3_Cell_identity")
DIR_F3_SUPP     <- file.path(OUT_DIR, "Figure3_SUPPLEMENTARY")
DIR_F4          <- file.path(OUT_DIR, "Figure4_Regional_analysis")
DIR_F4_SUPP     <- file.path(OUT_DIR, "Figure4_SUPPLEMENTARY")

invisible(lapply(list(DIR_F1, DIR_F1_SUPP, DIR_F2, DIR_F2_SUPP, DIR_F3, DIR_F3_SUPP, DIR_F4, DIR_F4_SUPP),
                 dir.create, showWarnings = FALSE, recursive = TRUE))