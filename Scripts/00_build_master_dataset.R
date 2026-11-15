###############################################################################
# 00_build_master_dataset.R
#
# Builds a shareable master dataset from the raw DIA-NN pg_matrix and the
# plate design table: clean, well-documented files for handing to any
# collaborator.
#
# Produces, in Results/master_dataset/:
#   1. 3D_SCP_all_Data_Master_file.csv  — primary file: gene identity +
#      one column per sample with RAW intensities. No filtering, no
#      normalisation, no imputation.
#   2. master_sample_metadata.csv       — one row per sample: region,
#      batch, cell type, tier, well/cell position, cell size.
#   3. contamination_analysis.csv       — this project's contaminant
#      calls, kept separate and clearly labelled as interpretation.
#   4. pen_sc_systematic_screen_ALL_PROTEINS.csv — full PEN_SC screen
#      audit trail (every protein, not just flagged ones).
#   5. README.md                        — describes all of the above.
#
# Expects the project laid out as:
#   Data/      — raw instrument files, plate design table, literature refs
#   Scripts/   — this file, 01-05 scripts, and shared/ (setup.R, functions.R)
#   Results/   — created if missing; master_dataset/ and figures/ written here
#
# Run with the R working directory set to the project folder (or set
# options(scp.project_root = "/path/to/project") first). Every path below is
# built from PROJECT_ROOT with file.path().
###############################################################################

PROJECT_ROOT <- getOption("scp.project_root", default = normalizePath(getwd(), winslash = "/"))

library(tidyverse)
library(readxl)

# ---- Config -----------------------------------------------------------

DATA_FILE   <- file.path(PROJECT_ROOT, "Data", "report_CombinedBrainEmpLib_05Jun2026_360_MBR.pg_matrix.tsv")
MASTER_FILE <- file.path(PROJECT_ROOT, "Data", "Plate_Design_Master_Table.xlsx")
OUT_DIR     <- file.path(PROJECT_ROOT, "Results", "master_dataset")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_PROTEINS_QC     <- 200    # informational flag only; not applied as a filter here
PEN_DETECTION_FRAC_MIN <- 5/7 # PEN_SC screen: min fraction of clean controls a gene must appear in
PEN_PERCENTILE_MIN     <- 69  # PEN_SC screen: min intensity percentile (dataset-wide) when detected

REGION_FULL <- c(CiC = "Cingulate Cortex", SoC = "Somatosensory Cortex",
                 Hippo = "Hippocampal Region", Cg = "Cerebellum Gray Matter",
                 Cw = "Cerebellum White Matter")

# Classic cRAP list (keratins/trypsin/hemoglobin/fibrinogen) and skin/
# cornified-envelope genes (Sharma non-overlap deep dive). PEN_SC-derived
# contaminants are NOT hardcoded here — they're recomputed fresh from the
# data every run by run_pen_contamination_screen() below.
CRAP_GENES_BASE <- c(
  "KRT1","KRT2","KRT3","KRT4","KRT5","KRT6A","KRT6B","KRT7","KRT8","KRT9",
  "KRT10","KRT13","KRT14","KRT15","KRT16","KRT17","KRT18","KRT19","KRT20",
  "Krt1","Krt2","Krt4","Krt5","Krt6a","Krt6b","Krt7","Krt8","Krt10",
  "Krt12","Krt13","Krt14","Krt15","Krt16","Krt17","Krt18","Krt19","Krt23",
  "Krt24","Krt25","Krt26","Krt27","Krt28","Krt31","Krt32","Krt33a","Krt33b",
  "Krt34","Krt35","Krt36","Krt42","Krt71","Krt72","Krt73","Krt74","Krt75",
  "Krt76","Krt77","Krt78","Krt79","Krt80","Krt82","Krt84","Krt85","Krt86",
  "Krt87","Krt90","Krt222","Krtap8-1","Krtcap2",
  "PRSS1","PRSS2","TRY1","TRY2",
  "Hbb-bs","Hbb-b1","Hbb-b2","Hba","Hba-a1","Hba-a2","Fga","Fgb","Fgg"
)
CRAP_GENES_SKIN_NEW <- c("Cdsn","Dsg1b","Pkp1","Sdr9c7","Serpinb5","Alox12b","Calml3")

# Non-filtering audit columns (recorded but never used to exclude a
# protein) — same patterns used in the Figure 3 fgsea audit, defined once
# here so they can't drift out of sync.
RIBO_PROTEASOME_PATTERN <- "^(Rpl|Rps|Psma|Psmb|Psmc|Psmd)"
BLOOD_IMMUNE_GENES <- c("Hmox1","Hbb.*","Hba.*","Alb","Apoa1","Apoc1","Apoc3","Ctsg","Mpo","Ltf",
                        "Serpinc1","Serpina1.*","Anxa3","Gc","Hpx","C3","C4b","Lbp","Fabp4","Cd36",
                        "G6pdx","Vwf","Serpinb1a","Serpinb1c",
                        "Cp","Itih.*","Hp","Ambp","Fga","Fgb","Fgg","Apoa2","Apob","Ttr","A2m","Orm1","Orm2","Trf")
BLOOD_IMMUNE_PATTERN <- paste0("^(", paste(BLOOD_IMMUNE_GENES, collapse = "|"), ")$")
MYELIN_OLIGO_PANEL <- c("Mbp", "Plp1", "Mog", "Cnp")

# ---- Step 1: parse sample columns & assign tiers -----------------------

#' Extract batch / well / region / cell-ID / control-status flags from a
#' raw DIA-NN sample column name. The trailing number is the global run index.
parse_sample_column <- function(col) {
  base <- basename(gsub("\\\\", "/", col))
  batch <- case_when(
    str_detect(col, "01Apr2026") ~ "01Apr2026",
    str_detect(col, "25Apr2026") ~ "24Apr2026",
    str_detect(col, "20260529")  ~ "20260529",
    str_detect(col, "30Apr2026") ~ "30Apr2026",
    TRUE ~ NA_character_
  )
  m_well <- str_match(base, "_([A-H])(\\d{1,2})_")
  m_cell <- str_match(base, "(CiC|SoC|Hippo|Cg|Cw)_Cell_(\\d+)")
  m_injection <- str_match(base, "_1_(\\d+)\\.d$")

  tibble(
    sample_col = col, base = base, batch = batch,
    well_row = m_well[, 2], well_col = suppressWarnings(as.integer(m_well[, 3])),
    region_code = m_cell[, 2], cell_num = suppressWarnings(as.integer(m_cell[, 3])),
    is_bulk     = str_detect(base, "_Bulk_") & !str_detect(str_to_upper(base), "PEN"),
    is_pen_sc   = str_detect(base, regex("PEN_SC", ignore_case = TRUE)),
    is_pen_bulk = str_detect(base, regex("PEN_Bulk", ignore_case = TRUE)),
    is_blank    = str_detect(base, regex("Blank", ignore_case = TRUE)),
    is_auditory = str_detect(base, regex("Auditory|Visual", ignore_case = TRUE)),
    injection_order = suppressWarnings(as.integer(m_injection[, 2]))
  )
}

#' Join parsed sample columns to the plate design table via a three-tier
#' fallback (type-matched well -> region+cell number -> plain well
#' position).
join_sample_metadata <- function(parsed, master) {
  by_well <- parsed %>%
    mutate(expected_type_cat = case_when(
      is_blank ~ "blank", is_pen_sc ~ "pen_sc", is_pen_bulk ~ "pen_bulk",
      is_bulk ~ "bulk", TRUE ~ "single_cell"
    )) %>%
    left_join(
      master %>%
        mutate(master_type_cat = case_when(
          str_detect(str_to_lower(sample_type), "blank") ~ "blank",
          str_detect(str_to_lower(sample_type), "pen_sc|pen sc") ~ "pen_sc",
          str_detect(str_to_lower(sample_type), "pen_bulk|pen bulk") ~ "pen_bulk",
          str_detect(str_to_lower(sample_type), "^bulk") ~ "bulk",
          TRUE ~ "single_cell"
        )) %>%
        select(batch, well_row, well_col, master_type_cat, sample_type, region,
               area_um2_final, area_based_cell_equivalent, include_in_study),
      by = c("batch", "well_row", "well_col", "expected_type_cat" = "master_type_cat")
    ) %>%
    select(-expected_type_cat)

  master_cell <- master %>% filter(!is.na(cell_num))
  by_cell <- parsed %>%
    left_join(master_cell %>% select(batch, region_code, cell_num, sample_type, region,
                                     area_um2_final, area_based_cell_equivalent, include_in_study),
             by = c("batch", "region_code", "cell_num"))

  merged <- by_well
  fallback_cols <- c("sample_type", "region", "area_um2_final",
                     "area_based_cell_equivalent", "include_in_study")

  # Fallback 1: region + cell number, for single cells the type-matched join missed.
  need_fallback <- is.na(merged$sample_type) & !is.na(parsed$cell_num)
  matched_via_fallback <- rep(FALSE, nrow(merged))
  if (any(need_fallback)) {
    fb <- by_cell %>% filter(sample_col %in% merged$sample_col[need_fallback])
    idx <- match(merged$sample_col[need_fallback], fb$sample_col)
    for (col in fallback_cols) merged[[col]][need_fallback] <- fb[[col]][idx]
    matched_via_fallback[need_fallback] <- TRUE
  }

  # Fallback 2 (last resort): plain well position, no sample-type disambiguation.
  # Rescues samples only — never overrides a successful type-matched join.
  still_missing <- is.na(merged$sample_type)
  matched_via_welltype_fallback <- rep(FALSE, nrow(merged))
  if (any(still_missing)) {
    well_only <- master %>%
      group_by(batch, well_row, well_col) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(batch, well_row, well_col, all_of(fallback_cols))
    wo <- parsed %>% filter(sample_col %in% merged$sample_col[still_missing]) %>%
      left_join(well_only, by = c("batch", "well_row", "well_col"))
    idx <- match(merged$sample_col[still_missing], wo$sample_col)
    for (col in fallback_cols) merged[[col]][still_missing] <- wo[[col]][idx]
    matched_via_welltype_fallback[still_missing] <- TRUE
    message("  NOTE: ", sum(still_missing), " sample(s) needed the well-position-only fallback; ",
            "recovered ", sum(still_missing & !is.na(merged$sample_type)),".")
  }

  # Safety check: exactly one output row per input sample column.
  if (nrow(merged) != nrow(parsed)) {
    dup_cols <- merged$sample_col[duplicated(merged$sample_col) | duplicated(merged$sample_col, fromLast = TRUE)]
    warning("Join produced ", nrow(merged), " rows from ", nrow(parsed), " raw samples -- ",
            length(unique(dup_cols)), " sample_col(s) matched more than once even after fallback. ",
            "Likely a duplicate batch+well+type row in Plate_Design_Master_Table.xlsx. Affected: ",
            paste(unique(dup_cols), collapse = "; "))
  }
  n_unmatched <- sum(is.na(merged$sample_type))
  if (n_unmatched > 0) {
    warning(n_unmatched, " sample(s) still unmatched after all fallbacks -- will fall through to ",
            "tier='unclassified' and be silently excluded downstream. Inspect: ",
            paste(merged$sample_col[is.na(merged$sample_type)], collapse = "; "))
  }

  merged %>% mutate(matched_via_region_cellnum_fallback = matched_via_fallback,
                    matched_via_welltype_fallback = matched_via_welltype_fallback)
}

#' Assign cell type and analysis tier. Control/blank status is checked
#' before include_in_study.
assign_tiers <- function(merged) {
  merged %>%
    mutate(
      cell_type = case_when(
        batch %in% c("01Apr2026", "24Apr2026", "20260529") ~ "Thy1",
        batch == "30Apr2026" ~ "Microglia",
        TRUE ~ NA_character_
      ),
      tier = case_when(
        is_pen_sc | is_pen_bulk ~ "PEN_control",
        is_blank ~ "blank",
        !is.na(include_in_study) & !include_in_study ~ "excluded_offtarget_region",
        is_auditory ~ "excluded_offtarget_region",
        sample_type == "Bulk" & batch == "01Apr2026" ~ "tier_10x_area",
        sample_type == "Bulk" & batch == "20260529"  ~ "tier_90x_area",
        sample_type == "Single_Cell" ~ "tier_1x_single_cell",
        TRUE ~ "unclassified"
      ),
      known_caveat = case_when(
        tier == "tier_90x_area" ~ paste("Within-run instrument degradation documented for this tier",
                                        "-- excluded from concordance/reproducibility figures elsewhere",
                                        "in this project, but included here unfiltered"),
        matched_via_welltype_fallback ~ paste("Metadata recovered via plain well-position fallback",
                                              "(no sample-type match found) -- verify sample_type/region"),
        TRUE ~ NA_character_
      )
    )
}

# ---- Step 2: PEN_SC contamination screen --------------------------------

#' Flag proteins as PEN_SC background contamination if detected in most/
#' all clean PEN_SC controls at high intensity (one outlier control well is
#' left out of the screen).
run_pen_contamination_screen <- function(mat_raw, sample_cols, out_dir) {
  pen_sc_cols <- sample_cols[str_detect(sample_cols, regex("PEN_SC", ignore_case = TRUE))]
  pen_sc_cols <- pen_sc_cols[!str_detect(pen_sc_cols, "G4_474")]  # outlier control well, left out
  message("  PEN_SC controls used (1 outlier excluded, G4): ", length(pen_sc_cols))

  if (length(pen_sc_cols) == 0) {
    message("  NOTE: no PEN_SC columns found -- screen skipped, no PEN-derived contaminants flagged.")
    results <- tibble(gene_label = character(), n_pen_detected = integer(),
                      pen_detection_frac = double(), pen_percentile = double(),
                      meets_pen_contamination_bar = logical())
  } else {
    pen_sub <- mat_raw[, pen_sc_cols, drop = FALSE]
    n_pen_detected <- rowSums(!is.na(pen_sub))
    pen_median_val <- apply(pen_sub, 1, median, na.rm = TRUE)

    all_signal_sorted <- sort(as.vector(mat_raw)[!is.na(as.vector(mat_raw))])
    pen_percentile <- ifelse(is.na(pen_median_val), NA_real_,
                             100 * findInterval(pen_median_val, all_signal_sorted) / length(all_signal_sorted))

    results <- tibble(
      gene_label = rownames(mat_raw),
      n_pen_detected = n_pen_detected,
      n_pen_total = length(pen_sc_cols),
      pen_detection_frac = n_pen_detected / length(pen_sc_cols),
      pen_percentile = pen_percentile
    ) %>%
      mutate(meets_pen_contamination_bar = !is.na(pen_percentile) &
               pen_detection_frac >= PEN_DETECTION_FRAC_MIN & pen_percentile >= PEN_PERCENTILE_MIN)

    message("  PEN_SC screen: ", sum(results$meets_pen_contamination_bar), " of ", nrow(results),
            " proteins flagged (>=", round(PEN_DETECTION_FRAC_MIN * length(pen_sc_cols)), "/",
            length(pen_sc_cols), " controls, >=", PEN_PERCENTILE_MIN, "th percentile).")
  }

  write_csv(results, file.path(out_dir, "pen_sc_systematic_screen_ALL_PROTEINS.csv"))
  results
}

# ---- Step 3: contamination analysis (this project's own judgment) ------

build_crap_lookup <- function(pen_screen_results) {
  pen_hits <- pen_screen_results %>% filter(meets_pen_contamination_bar) %>% arrange(desc(pen_percentile))
  tibble(
    gene_label = c(CRAP_GENES_BASE, CRAP_GENES_SKIN_NEW, pen_hits$gene_label),
    crap_reason = c(
      rep("keratin/trypsin/hemoglobin/fibrinogen (classic cRAP list)", length(CRAP_GENES_BASE)),
      rep("skin/cornified envelope (Sharma non-overlap deep dive)", length(CRAP_GENES_SKIN_NEW)),
      paste0("PEN_SC background (systematic screen: ", pen_hits$n_pen_detected, "/", pen_hits$n_pen_total,
             " controls, ", round(pen_hits$pen_percentile, 1), "th percentile)")
    )
  ) %>%
    distinct(gene_label, .keep_all = TRUE)  # in case a gene meets more than one category
}

#' Builds contamination_analysis.csv: the filtering flag
#' (is_crap_or_contaminant) plus non-filtering audit columns for
#' ribosomal/proteasome and blood/immune gene patterns.
build_contamination_analysis <- function(raw, crap_lookup) {
  raw %>%
    distinct(gene_label, .keep_all = TRUE) %>%
    select(gene_label) %>%
    left_join(crap_lookup, by = "gene_label") %>%
    mutate(
      is_igkv_family = str_detect(gene_label, "^Igkv"),
      is_crap_or_contaminant = !is.na(crap_reason) | is_igkv_family,
      crap_reason = if_else(is_igkv_family,
                            "Immunoglobulin kappa variable gene family (blood/immune contamination, matched by prefix)",
                            crap_reason),
      matches_ribo_proteasome_pattern = str_detect(gene_label, RIBO_PROTEASOME_PATTERN),
      matches_blood_immune_pattern = str_detect(gene_label, BLOOD_IMMUNE_PATTERN),
      is_myelin_oligo_panel_gene = gene_label %in% MYELIN_OLIGO_PANEL
    )
}

# ---- Step 4: sample metadata & primary intensity file -------------------

#' Build a human-readable sample name from plate-design facts (used as the
#' intensity matrix's column headers) and finalize sample_metadata.
finalize_sample_metadata <- function(merged, mat_raw) {
  n_proteins_raw <- colSums(!is.na(mat_raw))

  region_abbrev_lookup <- setNames(names(REGION_FULL), REGION_FULL)
  tier_abbrev_lookup <- c(
    tier_1x_single_cell = "1x", tier_10x_area = "10x", tier_90x_area = "90x",
    PEN_control = "PENctrl", blank = "Blank", excluded_offtarget_region = "OffTarget"
  )

  sample_metadata <- merged %>%
    mutate(
      n_proteins_detected = n_proteins_raw[sample_col],
      passes_qc_200protein_threshold = !is.na(n_proteins_detected) & n_proteins_detected >= MIN_PROTEINS_QC
    ) %>%
    rename(sample_id = sample_col, raw_filename = base) %>%
    mutate(
      region_abbrev = coalesce(region_abbrev_lookup[region], "NA"),
      tier_abbrev = coalesce(tier_abbrev_lookup[tier], "Unclass"),
      well_id = if_else(!is.na(well_row) & !is.na(well_col),
                        paste0(well_row, sprintf("%02d", well_col)), NA_character_),
      cell_id = if_else(!is.na(cell_num), paste0("Cell", sprintf("%03d", cell_num)), NA_character_),
      # Combine well + cell ID when both exist, rather than picking one.
      id_part = case_when(
        !is.na(well_id) & !is.na(cell_id) ~ paste0(well_id, "_", cell_id),
        !is.na(well_id) ~ well_id,
        !is.na(cell_id) ~ cell_id,
        TRUE ~ paste0("Row", row_number())
      ),
      clean_sample_name = paste(coalesce(cell_type, "NA"), region_abbrev, tier_abbrev, batch, id_part, sep = "_")
    ) %>%
    select(-region_abbrev, -tier_abbrev, -well_id, -cell_id, -id_part)

  # Resolve any name collisions with a numeric suffix rather than silently
  # overwriting a matrix column.
  dupe_mask <- duplicated(sample_metadata$clean_sample_name) |
    duplicated(sample_metadata$clean_sample_name, fromLast = TRUE)
  if (any(dupe_mask)) {
    sample_metadata <- sample_metadata %>%
      group_by(clean_sample_name) %>%
      mutate(.dup_group_n = n(), .dup_row_num = row_number()) %>%
      ungroup() %>%
      mutate(clean_sample_name = ifelse(.dup_group_n > 1,
                                        paste0(clean_sample_name, "_dup", .dup_row_num),
                                        clean_sample_name)) %>%
      select(-.dup_group_n, -.dup_row_num)
    message("  NOTE: ", sum(dupe_mask), " sample name collision(s) resolved with '_dupN' suffixes -- ",
            "check Plate_Design_Master_Table.xlsx for a duplicate batch/well row.")
  }

  sample_metadata %>%
    select(sample_id, clean_sample_name, raw_filename, batch, well_row, well_col, region, cell_type,
           sample_type, tier, cell_num, area_um2_final, area_based_cell_equivalent, include_in_study,
           is_bulk, is_pen_sc, is_pen_bulk, is_blank, is_auditory, injection_order,
           matched_via_region_cellnum_fallback, matched_via_welltype_fallback, n_proteins_detected,
           passes_qc_200protein_threshold, known_caveat) %>%
    arrange(batch, well_row, well_col)
}

#' Assemble the primary, judgment-free intensity file: protein identity +
#' one column per sample, raw intensities only. Columns keep their
#' ORIGINAL raw sample_id names (NOT renamed to clean_sample_name) -- this
#' file is meant to stay exactly traceable to the instrument output, with
#' no parsing/naming judgment applied. Use master_sample_metadata.csv's
#' sample_id <-> clean_sample_name columns to translate when needed;
#' setup.R does that translation downstream, for analysis use.
build_master_intensity_file <- function(raw, mat_raw, sample_metadata) {
  protein_identity <- raw %>%
    select(gene_label, Protein.Group, Protein.Names, First.Protein.Description,
           N.Sequences, N.Proteotypic.Sequences) %>%
    distinct(gene_label, .keep_all = TRUE)

  protein_identity %>%
    left_join(as.data.frame(mat_raw) %>% rownames_to_column("gene_label"), by = "gene_label")
}

# ---- Step 5: README ------------------------------------------------------

write_release_readme <- function(out_dir) {
  readme_text <- glue::glue('
# 3D-SCP Brain Proteomics — Master Dataset

Generated: {format(Sys.Date())}

## What this is

Single-cell and pooled-area (Bulk) proteomics of DISCO-cleared mouse brain
tissue, profiling Thy1-GFP-M neurons (Cingulate Cortex, Somatosensory
Cortex, Hippocampal Region) and CX3CR1-eGFP microglia (Somatosensory
Cortex, Cerebellum Gray Matter, Cerebellum White Matter), collected by
laser capture microdissection (LCM) and analysed by timsTOF DIA-NN mass
spectrometry.

Samples span four acquisition batches: 01Apr2026, 24Apr2026 (raw files
dated 25Apr2026), 30Apr2026 (CX3CR1/microglia), and 20260529.

## Files in this release

### 1. `3D_SCP_all_Data_Master_file.csv` — the primary, unbiased data file
One row per protein. **No analytical judgment calls in this file** — pure
gene identity and raw measurements only.

**Leading columns (protein identity — descriptive, not judgment):**
- `gene_label` — gene symbol (first gene if a protein group maps to
  multiple; falls back to the UniProt Protein.Group accession if no gene
  symbol is available)
- `Protein.Group`, `Protein.Names`, `First.Protein.Description` — original
  DIA-NN/UniProt annotation
- `N.Sequences`, `N.Proteotypic.Sequences` — peptide evidence counts

**Sample columns (one per sample, raw intensities):**
Named with the ORIGINAL raw DIA-NN sample identifier (`sample_id` in file
2) -- deliberately NOT renamed to a readable form here, so this file stays
exactly traceable to the instrument output with no parsing/naming
judgment applied to it. To work with human-readable sample names, join
against `master_sample_metadata.csv`\'s `sample_id` <-> `clean_sample_name`
columns.

**Values are raw, linear-scale MS intensities as reported by DIA-NN.**
DIA-NN encodes non-detection as `0`; this has been converted to a proper
empty/NA value. No contaminant filtering, no normalisation, no log
transform, and no imputation have been applied.

### 2. `master_sample_metadata.csv` — objective, measured facts per sample
One row per sample — **all samples included**, including blanks, PEN
controls, and off-target-region samples; nothing silently dropped.
- `sample_id` — original DIA-NN identifier; join key to the intensity
  matrix column headers in file 1
- `clean_sample_name` — readable identifier built from objective
  plate-design facts, for human use:
  ```
  {{CellType}}_{{RegionAbbrev}}_{{Tier}}_{{Batch}}_{{WellID_and/or_CellID}}
  ```
  - `CellType`: `Thy1` or `Microglia`
  - `RegionAbbrev`: `CiC` (Cingulate Cortex), `SoC` (Somatosensory Cortex),
    `Hippo` (Hippocampal Region), `Cg` (Cerebellum Gray Matter), `Cw`
    (Cerebellum White Matter)
  - `Tier`: `1x` (single cell), `10x`/`90x` (Bulk area tiers — see caveat #2
    below), `PENctrl`, `Blank`, `OffTarget`, or `Unclass`
  - `Batch`: acquisition batch, e.g. `01Apr2026`
  - `WellID_and/or_CellID`: plate well (e.g. `A03`), LCM cell number (e.g.
    `Cell014`), or both concatenated (`A03_Cell014`) when a sample has both

  Examples: `Thy1_CiC_1x_01Apr2026_A03`, `Microglia_Cg_1x_30Apr2026_Cell014`,
  `Thy1_CiC_1x_01Apr2026_A03_Cell014`
- `raw_filename` — base filename, for traceability back to the source data
- `batch`, `well_row`, `well_col`, `region`, `cell_type`, `sample_type`,
  `tier`, `cell_num` — plate design / experimental design facts
- **`area_um2_final`, `area_based_cell_equivalent`** — LCM-measured cell/
  collection size. Cell size is a known driver of single-cell proteome
  variation (Lanz et al. 2023, J Proteome Res), and correlates strongly
  with per-cell protein depth in this dataset — check this column early
  for any PCA or variance analysis.
- `is_bulk`, `is_pen_sc`, `is_pen_bulk`, `is_blank`, `is_auditory` — raw
  filename-parsed flags used to derive `tier`
- `matched_via_region_cellnum_fallback`, `matched_via_welltype_fallback`
  — TRUE if a sample\'s metadata came from a join fallback rather than
  the primary type-matched well join
- `n_proteins_detected` — protein count per sample, from the raw
  (unfiltered) intensity matrix
- `known_caveat` — free-text notes on documented technical issues (e.g.
  the 20260529 Bulk tier\'s within-run instrument degradation)

### 3. `contamination_analysis.csv` — this project\'s own interpretive work
Kept **deliberately separate** from the neutral data above. One row per
protein:
- `gene_label` — join key back to file 1
- `is_crap_or_contaminant` — TRUE if flagged on this project\'s
  cRAP/contaminant watch-list. **This reflects analytical judgment, not
  verified fact** — use it, ignore it, or redo the analysis yourself from
  the raw data in file 1.
- `crap_reason` — which category flagged it (classic cRAP list;
  skin/cornified-envelope; PEN_SC background, with the exact detection
  count and percentile recorded per gene; or Igkv immunoglobulin family).
  See `pen_sc_systematic_screen_ALL_PROTEINS.csv` for the full screen
  results across every protein, not just flagged ones.
- `matches_ribo_proteasome_pattern`, `matches_blood_immune_pattern`,
  `is_myelin_oligo_panel_gene` — non-filtering audit columns; these do
  NOT affect `is_crap_or_contaminant`.

### 4. This README

## Caveats worth knowing before you analyse this data

1. **`cell_type` is an experimental design label, not a purity claim.**
   The 01Apr2026 Bulk wells (`tier_10x_area`) are unsorted ~10-cell-
   equivalent LCM tissue punches from Thy1-targeted regions — labelled
   "Thy1" because that\'s the traced lineage, but the tissue physically
   contains whatever cell types sit in that area.
2. **The 20260529 Bulk tier (`tier_90x_area`) shows within-run instrument
   degradation** (protein counts decline across well position within a
   replicate block). Excluded from concordance/reproducibility figures
   elsewhere in this project but included, unfiltered, in this release.
3. **The regional design is not fully crossed**: microglia were only
   collected from Cerebellum Gray Matter, Cerebellum White Matter, and
   Somatosensory Cortex; Thy1 neurons only from Cingulate Cortex,
   Somatosensory Cortex, and Hippocampal Region. Somatosensory Cortex is
   the only region shared between both cell types.
4. **cRAP/contaminant flagging is this project\'s own working list** — a
   defensible starting point, not a definitive one. Inspect `crap_reason`
   and adjust for your own analysis.
5. **Off-target wells are excluded via `include_in_study`**, not filename
   text matching.

')
  writeLines(readme_text, file.path(out_dir, "README.md"))
}

# ---- Step 6: console summary ---------------------------------------------

print_sample_count_summary <- function(sample_metadata) {
  message("\n", strrep("=", 70), "\nSAMPLE COUNT SUMMARY -- check this against what you expect\n", strrep("=", 70))

  message("\nBy tier (all batches combined):")
  print(sample_metadata %>% count(tier, name = "n") %>% arrange(desc(n)), n = Inf)

  message("\nBy tier x batch:")
  print(sample_metadata %>% count(tier, batch) %>% pivot_wider(names_from = batch, values_from = n, values_fill = 0), n = Inf)

  message("\nSingle cells by cell type x region (QC not yet applied):")
  print(sample_metadata %>% filter(tier == "tier_1x_single_cell") %>%
          count(cell_type, region) %>% pivot_wider(names_from = region, values_from = n, values_fill = 0), n = Inf)

  message("\nSingle cells passing QC (>=", MIN_PROTEINS_QC, " proteins) by cell type x region:")
  print(sample_metadata %>% filter(tier == "tier_1x_single_cell", passes_qc_200protein_threshold) %>%
          count(cell_type, region) %>% pivot_wider(names_from = region, values_from = n, values_fill = 0), n = Inf)

  message("\nControls/excluded: Blanks=", sum(sample_metadata$tier == "blank"),
          " | PEN controls=", sum(sample_metadata$tier == "PEN_control"),
          " | Off-target excluded=", sum(sample_metadata$tier == "excluded_offtarget_region"),
          " | Unclassified (should be 0)=", sum(sample_metadata$tier == "unclassified"))
  message(strrep("=", 70), "\n")
}

# ---- Run pipeline ---------------------------------------------------------

message("\n[00] Building shareable master dataset...")

message("  Parsing sample columns & joining plate design...")
master <- read_excel(MASTER_FILE, sheet = "Master Well Table") %>%
  mutate(region_code = names(REGION_FULL)[match(region, REGION_FULL)])

header_line <- readLines(DATA_FILE, n = 1)
all_cols <- str_split(header_line, "\t")[[1]]
anno_cols <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description",
              "N.Sequences", "N.Proteotypic.Sequences")
sample_cols <- setdiff(all_cols, anno_cols)

parsed <- map_dfr(sample_cols, parse_sample_column)
merged <- join_sample_metadata(parsed, master) %>% assign_tiers()

message("  Loading raw intensity matrix...")
raw <- read_tsv(DATA_FILE, show_col_types = FALSE) %>%
  mutate(gene_label = coalesce(str_extract(Genes, "^[^;]+"), Protein.Group))
mat_raw <- raw %>% select(all_of(sample_cols)) %>% as.matrix()
rownames(mat_raw) <- raw$gene_label
mat_raw[mat_raw == 0] <- NA   # DIA-NN's non-detection placeholder -> NA

pen_screen_results <- run_pen_contamination_screen(mat_raw, sample_cols, OUT_DIR)
crap_lookup <- build_crap_lookup(pen_screen_results)

sample_metadata <- finalize_sample_metadata(merged, mat_raw)
write_csv(sample_metadata, file.path(OUT_DIR, "master_sample_metadata.csv"))
message("  Wrote master_sample_metadata.csv (", nrow(sample_metadata), " samples)")

offtarget_by_batch <- sample_metadata %>% filter(tier == "excluded_offtarget_region") %>% count(batch, name = "n")
message("  Off-target exclusion: ", sum(sample_metadata$tier == "excluded_offtarget_region"), " samples (by batch: ",
        paste0(offtarget_by_batch$batch, "=", offtarget_by_batch$n, collapse = ", "),
        ") -- sanity-check against expectation (16 in 01Apr2026); anything elsewhere is worth a closer look.")

message("  Building contamination analysis (kept separate -- this is judgment, not fact)...")
contamination_analysis <- build_contamination_analysis(raw, crap_lookup)
write_csv(contamination_analysis, file.path(OUT_DIR, "contamination_analysis.csv"))
message("  Wrote contamination_analysis.csv (", nrow(contamination_analysis), " proteins, ",
        sum(contamination_analysis$is_crap_or_contaminant), " flagged)")

message("  Writing primary data file (neutral: gene identity + raw intensities only)...")
master_combined <- build_master_intensity_file(raw, mat_raw, sample_metadata)
write_csv(master_combined, file.path(OUT_DIR, "3D_SCP_all_Data_Master_file.csv"))
message("  Wrote 3D_SCP_all_Data_Master_file.csv (", nrow(master_combined), " proteins x ",
        ncol(master_combined) - 6, " samples)")

message("  Writing README.md...")
write_release_readme(OUT_DIR)

print_sample_count_summary(sample_metadata)
message("  Master dataset release complete: ", OUT_DIR)
