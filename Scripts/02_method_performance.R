###############################################################################
# 02_method_performance.R -- 3D-SCP method development and sensitivity
#
# Generates the quantitative panels of Figure 2 and Extended Data (ED) Fig. 3:
#   a    LC column comparison (HeLa, 1 cell equivalent): protein groups and peak capacity
#   b    input-amount titration (HeLa, 0.1-100 cell equivalents)
#   c    sensitivity: Thy1 single cells vs 10-cell areas
#   d    subcellular compartment coverage of the detected proteins (UniProt)
#   e    protein detection overlap, single cells vs 10-cell areas (Venn)
#   ED Fig. 3a-e  protein yield by batch; batch-to-batch concordance; collection area and yield by
#                 cell type; precursor-level intensity scaling (10-cell area vs single cells);
#                 single-cell average vs 10-cell area by region
# Panel letters are added afterwards in a vector editor.
#
# Sections:
#   1. Setup
#   2. Cell counts and cohort flow
#   3. Fig. 2a, 2b: HeLa test injections
#   4. Fig. 2c: sensitivity
#   5. Fig. 2d: subcellular compartments
#   6. Fig. 2e: protein overlap (Venn)
#   7. Main Figure 2: assembly
#   8. Extended Data Fig. 3
#   9. Summary numbers
#
# Cell inclusion (contaminant-filtered protein groups per cell; MIN_PROTEINS_VIABLE and MIN_PROTEINS
# in setup.R): >= 200 for the descriptive yield panels (Fig. 2c-e, ED 3a, 3c, 3d, and 10-cell areas);
# >= 600 (analysis set) wherever correlations are computed (ED 3b, 3e) and in Figs. 3 and 4.
# Never the raw count before contaminant removal. 
#
# ED 3d and 3e compare intensities between 1 and 10 cells and therefore use the precursor-level
# matrix (raw Precursor.Quantity), not MaxLFQ: MaxLFQ assumes that most proteins do not change
# between the compared samples, which is violated by design here. Fig. 2a-e and ED 3a, 3c use
# protein-group counts or presence/absence only; ED 3b compares protein-level (MaxLFQ) averages
# between batches.
#
# Inputs : the master dataset and cohorts from setup.R; Data/*stats.tsv (HeLa injections);
#          the DIA-NN precursor report (needed for ED 3d and 3e); the UniProt web service
#          (results cached in Data/uniprot_subcellular_location_cache.csv).
# Outputs: Figure2_MAIN_assembled.svg, EDFig3_assembled.svg, the CSV tables of every panel,
#          Figure2_sample_membership.csv and Figure2_summary_numbers.txt.
###############################################################################
# Run from the project folder, or set options(scp.project_root = "/path/to/project") first.
source(file.path(getOption("scp.project_root", getwd()), "Scripts", "shared", "setup.R"))

###############################################################################
# 1. SETUP
# Main Figure 2 is 6.8 in wide (room at the side for extra artwork): 2 rows x 3 panels, each row
# FIG2_ROW_H_IN high. Significance stars are 11 pt. Text size and theme come from setup.R.
###############################################################################
FIG2_W_IN     <- 6.8
FIG2_ROW_H_IN <- 2.7
STAR_SIZE  <- pt_to_geom(11)                 # significance stars 10-12 pt

###############################################################################
# 2. CELL COUNTS AND COHORT FLOW
###############################################################################
cells_collected <- cohort_meta %>% dplyr::filter(tier == "tier_1x_single_cell")
n_cells_collected <- nrow(cells_collected)
collected_cols <- cells_collected$sample_col

# Protein counts per cell are the contaminant-filtered counts already stored in cohort_meta
# (n_proteins_filtered, from setup.R).
cols_in_mat <- intersect(collected_cols, colnames(protein_mat_log2))

single_cell_max_check <- cells_collected %>%
  dplyr::filter(passes_qc_tier, is.finite(n_proteins_filtered)) %>%
  dplyr::group_by(cell_type_label) %>%
  dplyr::slice_max(n_proteins_filtered, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(cell_type_label, sample_col, n_proteins_filtered)

n_proteins_ever_measured <- sum(rowSums(!is.na(protein_mat_log2[, cols_in_mat, drop = FALSE])) > 0)
cells_used <- cells_collected %>% dplyr::filter(is.finite(n_proteins_filtered), n_proteins_filtered >= MIN_PROTEINS)
n_cells_used <- nrow(cells_used)

# Empty/fully-failed wells: collected cells with zero detected proteins
# after contaminant filtering (not just below the 600-protein analysis
# threshold -- genuinely nothing detected at all). Also reports the
# broader low-count distribution for context, since "empty" vs. "just
# below threshold" is a real distinction worth seeing, not one fixed
# cutoff.
n_cells_empty <- sum(cells_collected$n_proteins_filtered == 0, na.rm = TRUE)
empty_well_breakdown <- cells_collected %>%
  dplyr::filter(is.finite(n_proteins_filtered)) %>%
  dplyr::mutate(count_band = dplyr::case_when(
    n_proteins_filtered == 0 ~ "0 (empty/failed)",
    n_proteins_filtered < 10 ~ "1-9",
    n_proteins_filtered < 50 ~ "10-49",
    n_proteins_filtered < MIN_PROTEINS ~ paste0("50-", MIN_PROTEINS - 1),
    TRUE ~ paste0(">=", MIN_PROTEINS)
  )) %>%
  dplyr::count(count_band, name = "n_cells")


###############################################################################
# 3. FIG. 2a, 2b: HeLa TEST INJECTIONS (no tissue cells)
###############################################################################
message("\n[Figure 2] Method development + data reliability...")

fresh_candidates <- list.files(file.path(PROJECT_ROOT, "Data"), pattern = "200fL_200pL_QCs_combiLib_17Jun2026.*stats.*\\.tsv$",
                               ignore.case = TRUE, full.names = TRUE)
if (length(fresh_candidates) == 0) {
  stop("HeLa titration stats file not found in ", file.path(PROJECT_ROOT, "Data"),
       " (pattern 200fL_200pL_QCs_combiLib_17Jun2026*stats*.tsv); it is needed for Fig 2a and 2b.", call. = FALSE)
} else {
  if (length(fresh_candidates) > 1) warning("Several titration stats files found; using the first: ", fresh_candidates[1])
  FRESH_TITRATION_FILE <- fresh_candidates[1]
  message("  Using fresh-cell titration file: ", FRESH_TITRATION_FILE)
  
  classify_fresh_sample <- function(path) {
    base <- basename(gsub("\\\\", "/", path))
    if (stringr::str_detect(base, "_13256")) return(c(NA, NA, NA))
    is_ht <- stringr::str_detect(path, "HT")
    column <- if (is_ht) "HT" else "PLUS"
    spd <- stringr::str_match(path, "(\\d+)SPD")[,2]
    load_group <- dplyr::case_when(
      stringr::str_starts(base, "25pg_QC_KT") ~ "0.1 cell",
      stringr::str_starts(base, "250pg_QC_23Apr2026") | stringr::str_starts(base, "250pg_QC_30Apr2026") ~ "1 cell",
      stringr::str_starts(base, "2_5ng_QC_23Apr") ~ "10 cell",
      stringr::str_starts(base, "12_5ng_QC_DPM_15May2026") ~ "50 cell",
      stringr::str_starts(base, "25ng_QC") ~ "100 cell",
      stringr::str_starts(base, "25pg_QC_17Apr2026") & is_ht ~ "0.1 cell",
      stringr::str_starts(base, "250pg_QC_17Apr2026") & is_ht ~ "1 cell",
      TRUE ~ NA_character_
    )
    c(load_group, column, spd)
  }
  
  fresh_stats <- read_tsv(FRESH_TITRATION_FILE, show_col_types = FALSE)
  classified <- purrr::map(fresh_stats$File.Name, classify_fresh_sample)
  fresh_stats$load_group <- purrr::map_chr(classified, 1)
  fresh_stats$column     <- purrr::map_chr(classified, 2)
  fresh_stats$spd        <- purrr::map_chr(classified, 3)
  fresh_clean <- fresh_stats %>% dplyr::filter(!is.na(load_group)) %>%
    dplyr::mutate(k_proteins = Proteins.Identified / 1000,
                  load_group = factor(load_group, levels = c("0.1 cell","1 cell","10 cell","50 cell","100 cell")))
  
  one_cell <- fresh_clean %>% dplyr::filter(load_group == "1 cell")
  col_summary <- one_cell %>% dplyr::group_by(column) %>%
    dplyr::summarise(mean_k = mean(k_proteins), sd_k = sd(k_proteins), n = dplyr::n(),
                     spd = paste(unique(spd), collapse=","), .groups = "drop")
  ht_mean <- col_summary$mean_k[col_summary$column == "HT"]
  plus_mean <- col_summary$mean_k[col_summary$column == "PLUS"]
  fold_increase <- plus_mean / ht_mean
  t_res <- t.test(k_proteins ~ column, data = one_cell)
  sig_stars <- dplyr::case_when(t_res$p.value < 0.001 ~ "***", t_res$p.value < 0.01 ~ "**",
                                t_res$p.value < 0.05 ~ "*", TRUE ~ "n.s.")
  col_summary$column <- factor(col_summary$column, levels = c("HT","PLUS"))
  one_cell$column <- factor(one_cell$column, levels = c("HT","PLUS"))
  bracket_y <- max(one_cell$k_proteins) * 1.12
  
  # Shared x-axis labels for the three 1-cell subplots: column + SPD + n
  # together, since SPD is fully confounded with column in this dataset
  # (HT: 40SPD only; PLUS: 20SPD only).
  spd_n_labels <- setNames(
    paste0(c("\u00b5PAC HT", "\u00b5PAC Neo Plus"), "\n",
           c("10 min", "45 min"), " gradient\nn = ",
           col_summary$n[match(c("HT","PLUS"), col_summary$column)]),
    c("HT","PLUS"))
  
  fig2a_columns <- ggplot(col_summary, aes(column, mean_k)) +
    geom_col(fill = "white", colour = "black", width = 0.55, linewidth = 0.6) +
    geom_errorbar(aes(ymin = mean_k - sd_k, ymax = mean_k + sd_k), width = 0.15) +
    geom_jitter(data = one_cell, aes(column, k_proteins), width = 0.08, size = 1.8, alpha = 0.7) +
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y) +
    annotate("text", x = 1.5, y = bracket_y * 1.03, label = sig_stars, size = STAR_SIZE) +
    annotate("label", x = 1.5, y = bracket_y * 1.22, label = paste0(round(fold_increase,2), "-fold"),
             fill = "#EDF0DE", colour = PAL_SHARMA, linewidth = 0, fontface = "bold", size = pt_to_geom(PT_BASE)) +
    scale_x_discrete(labels = spd_n_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +   # bars start exactly at 0; top space for the bracket/fold-label
    labs(x = NULL, y = "Protein groups detected (thousands)",
         title = "HeLa digest, 1 cell eq.") +
    THEME_JOURNAL()
  
  fwhm_summary <- one_cell %>% dplyr::group_by(column) %>%
    dplyr::summarise(mean_fwhm = mean(FWHM.RT), sd_fwhm = sd(FWHM.RT), n = dplyr::n(), .groups = "drop")
  fwhm_summary$column <- factor(fwhm_summary$column, levels = c("HT","PLUS"))
  fwhm_t <- t.test(FWHM.RT ~ column, data = one_cell)
  write_csv(fwhm_summary, file.path(DIR_F2, "2a_peak_width_summary.csv"))
  message("  Peak width: HT mean FWHM.RT=", round(fwhm_summary$mean_fwhm[fwhm_summary$column=="HT"], 4),
          "  PLUS mean FWHM.RT=", round(fwhm_summary$mean_fwhm[fwhm_summary$column=="PLUS"], 4),
          "  p=", signif(fwhm_t$p.value, 3))
  
  ms1_summary <- one_cell %>% dplyr::group_by(column) %>%
    dplyr::summarise(mean_ms1 = mean(MS1.Signal), sd_ms1 = sd(MS1.Signal), n = dplyr::n(), .groups = "drop")
  ms1_summary$column <- factor(ms1_summary$column, levels = c("HT","PLUS"))
  ms1_t <- t.test(MS1.Signal ~ column, data = one_cell)
  write_csv(ms1_summary, file.path(DIR_F2, "2a_ms1_signal_summary.csv"))
  message("  MS1 signal: HT mean=", signif(ms1_summary$mean_ms1[ms1_summary$column=="HT"], 3),
          "  PLUS mean=", signif(ms1_summary$mean_ms1[ms1_summary$column=="PLUS"], 3),
          "  p=", signif(ms1_t$p.value, 3))
  
  ###########################################################################
  # Peak capacity (Nc = 1 + gradient_time / FWHM.RT).
  # Gradient times below are nominal values (~10 min for HT, ~45 min for PLUS) -- not logged per-run in
  # stats.tsv, and not exact to the second.
  ###########################################################################
  GRADIENT_TIME_MIN <- c(HT = 10, PLUS = 45)   # nominal gradient length, minutes (per LC operator)
  one_cell <- one_cell %>%
    dplyr::mutate(gradient_time_min = GRADIENT_TIME_MIN[as.character(column)],
                  peak_capacity = 1 + gradient_time_min / FWHM.RT)
  peakcap_summary <- one_cell %>% dplyr::group_by(column) %>%
    dplyr::summarise(mean_peakcap = mean(peak_capacity), sd_peakcap = sd(peak_capacity), n = dplyr::n(), .groups = "drop")
  peakcap_summary$column <- factor(peakcap_summary$column, levels = c("HT", "PLUS"))
  peakcap_t <- t.test(peak_capacity ~ column, data = one_cell)
  peakcap_sig <- dplyr::case_when(peakcap_t$p.value < 0.001 ~ "***", peakcap_t$p.value < 0.01 ~ "**",
                                  peakcap_t$p.value < 0.05 ~ "*", TRUE ~ "n.s.")
  peakcap_bracket_y <- max(one_cell$peak_capacity) * 1.12
  
  fig2a_peakcap <- ggplot(peakcap_summary, aes(column, mean_peakcap)) +
    geom_col(fill = "white", colour = "black", width = 0.55, linewidth = 0.6) +
    geom_errorbar(aes(ymin = mean_peakcap - sd_peakcap, ymax = mean_peakcap + sd_peakcap), width = 0.15) +
    geom_jitter(data = one_cell, aes(column, peak_capacity), width = 0.08, size = 1.8, alpha = 0.7) +
    annotate("segment", x = 1, xend = 2, y = peakcap_bracket_y, yend = peakcap_bracket_y) +
    annotate("text", x = 1.5, y = peakcap_bracket_y * 1.03, label = peakcap_sig, size = STAR_SIZE) +
    scale_x_discrete(labels = spd_n_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +   # bars start exactly at 0
    labs(x = NULL, y = "Peak capacity (1 + gradient / FWHM)",
         title = "HeLa digest, 1 cell eq.") +
    THEME_JOURNAL()
  write_csv(peakcap_summary, file.path(DIR_F2, "2a_peak_capacity_summary.csv"))
  message("  Peak capacity (nominal gradient times, see code comment): HT mean=", round(peakcap_summary$mean_peakcap[peakcap_summary$column=="HT"], 1),
          "  PLUS mean=", round(peakcap_summary$mean_peakcap[peakcap_summary$column=="PLUS"], 1),
          "  p=", signif(peakcap_t$p.value, 3))
  
  write_csv(col_summary, file.path(DIR_F2, "2a_column_comparison_summary.csv"))
  
  message("  Fig 2b: input-amount titration...")
  plus_all <- fresh_clean %>% dplyr::filter(column == "PLUS")
  load_to_ce <- c("0.1 cell" = 0.1, "1 cell" = 1, "10 cell" = 10, "50 cell" = 50, "100 cell" = 100)
  load_summary_log <- plus_all %>% dplyr::group_by(load_group) %>%
    dplyr::summarise(mean_k = mean(k_proteins), sd_k = sd(k_proteins), n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(cell_equiv = load_to_ce[as.character(load_group)])
  
  set.seed(1)
  plus_all_log <- plus_all %>% dplyr::mutate(
    cell_equiv = load_to_ce[as.character(load_group)],
    # jitter range widened from +/-0.15 to +/-0.22 -- at the original width, replicate
    # points with near-identical values (e.g. within the 50/100-cell groups) sat almost
    # exactly on top of each other and were indistinguishable even with a point outline.
    cell_equiv_jit = cell_equiv * 2^stats::runif(dplyr::n(), -0.22, 0.22)
  )
  
  load_summary_log <- load_summary_log %>% dplyr::mutate(
    cap_xmin = cell_equiv / 1.15,
    cap_xmax = cell_equiv * 1.15
  )
  
  
  fig2b_titration <- ggplot(load_summary_log, aes(cell_equiv, mean_k)) +
    geom_segment(aes(x = cell_equiv, xend = cell_equiv, y = mean_k - sd_k, yend = mean_k + sd_k),
                 colour = "black", linewidth = 0.7) +
    geom_segment(aes(x = cap_xmin, xend = cap_xmax, y = mean_k + sd_k, yend = mean_k + sd_k),
                 colour = "black", linewidth = 0.7) +
    geom_segment(aes(x = cap_xmin, xend = cap_xmax, y = mean_k - sd_k, yend = mean_k - sd_k),
                 colour = "black", linewidth = 0.7) +
    geom_point(size = 2.3, colour = "black", shape = 21, fill = "white", stroke = 1, alpha = 0.9) +
    geom_point(data = plus_all_log, aes(cell_equiv_jit, k_proteins), size = 2.0, alpha = 0.8,
               shape = 21, fill = PAL_BULK, colour = "white", stroke = 0.9, inherit.aes = FALSE) +
    scale_x_log10(breaks = load_to_ce, labels = as.character(load_to_ce),
                  expand = expansion(mult = 0.15)) +
    labs(x = "HeLa digest (cell equivalents)", y = "Protein groups detected (thousands)") +
    THEME_JOURNAL() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # Diluted HeLa reference material (µPAC Neo Plus only). Open circle =
  # mean; colored points = raw replicates, jittered proportionally.
  
  write_csv(load_summary_log, file.path(DIR_F2, "2b_load_titration_summary.csv"))
}
###############################################################################
# 4. FIG. 2c: SENSITIVITY
# Thy1 single cells vs 10-cell areas (all cells with >= 200 protein groups; the dashed line marks the
# 600-protein analysis threshold used in Figs. 3 and 4).
###############################################################################
message("  Fig 2c: sensitivity (real tissue)...")
sens_rows <- list()
for (tier_name in c("tier_1x_single_cell", "tier_10x_area")) {
  sub <- cohort_meta %>% dplyr::filter(tier==tier_name, passes_qc_tier,
                                       tier_name!="tier_1x_single_cell" | cell_type_label=="Thy1") %>%
    dplyr::select(sample_col, batch)
  sens_rows[[tier_name]] <- tibble::tibble(tier=tier_name, sample_col=sub$sample_col, batch=sub$batch,
                                           n_proteins=n_proteins_lookup[sub$sample_col])
}
sens_df <- dplyr::bind_rows(sens_rows) %>%
  dplyr::mutate(tier_label = factor(dplyr::recode(tier, tier_1x_single_cell="Single cells", tier_10x_area="10-cell area"),
                                    levels=c("Single cells","10-cell area")),
                batch_label = dplyr::recode(batch, !!!BATCH_DISPLAY_LABELS))


sens_summary <- sens_df %>% dplyr::group_by(tier_label) %>%
  dplyr::summarise(median_n = median(n_proteins), n = dplyr::n(), .groups = "drop")
fold_sens <- sens_summary$median_n[sens_summary$tier_label == "10-cell area"] /
  sens_summary$median_n[sens_summary$tier_label == "Single cells"]
# How many of the >=200 cells also clear the 600-protein analysis threshold (Figs 3, 4)
sens_600 <- sens_df %>% dplyr::group_by(tier_label) %>%
  dplyr::summarise(n_600 = sum(n_proteins >= MIN_PROTEINS),
                   median_600 = median(n_proteins[n_proteins >= MIN_PROTEINS]), .groups = "drop")
fold_sens_600 <- sens_600$median_600[sens_600$tier_label == "10-cell area"] /
  sens_600$median_600[sens_600$tier_label == "Single cells"]
tier_n_labels <- setNames(
  paste0(levels(sens_df$tier_label), "\n(n = ",
         sens_summary$n[match(levels(sens_df$tier_label), sens_summary$tier_label)], "; ",
         sens_600$n_600[match(levels(sens_df$tier_label), sens_600$tier_label)], " \u2265600)"),
  levels(sens_df$tier_label))
bracket_y_G <- max(sens_df$n_proteins) * 1.12

fig2c_sensitivity <- ggplot(sens_df, aes(tier_label, n_proteins, fill=tier)) +
  geom_boxplot(outlier.shape=NA, alpha=0.6) +
  geom_jitter(aes(colour=batch_label), width=0.18, size=1.2, alpha=0.65, stroke = 0) +
  # dashed line = 600-protein threshold used for the analysis set in Figs 3 and 4;
  # everything shown here has >=200 proteins (descriptive QC bar)
  geom_hline(yintercept = MIN_PROTEINS, linetype = "dashed", colour = "grey30", linewidth = 0.45) +
  annotate("text", x = 2.55, y = MIN_PROTEINS, label = "600", hjust = 1, vjust = -0.5, size = pt_to_geom(PT_BASE), colour = "grey30") +
  annotate("segment", x = 1, xend = 2, y = bracket_y_G, yend = bracket_y_G) +
  annotate("label", x = 1.5, y = bracket_y_G * 1.10, label = paste0(round(fold_sens,2), "-fold (\u2265200)"),
           fill = "#EDF0DE", colour = PAL_SHARMA, linewidth = 0, fontface = "bold", size = pt_to_geom(PT_BASE)) +
  scale_fill_manual(values=c(tier_1x_single_cell=PAL_THY1, tier_10x_area=PAL_BULK), guide="none") +
  scale_colour_manual(values = setNames(BATCH_PALETTE, BATCH_DISPLAY_LABELS[names(BATCH_PALETTE)]), name = "Batch") +
  scale_x_discrete(labels = tier_n_labels) +
  labs(x=NULL, y="Protein groups detected") +
  THEME_JOURNAL() +
  theme(legend.position = c(0.02, 0.86), legend.justification = c(0, 1),
        legend.background = element_blank(), legend.key.size = grid::unit(2.5, "mm"),
        legend.spacing.y = grid::unit(0, "mm"))

###############################################################################
# 5. FIG. 2d: SUBCELLULAR COMPARTMENTS
# UniProt subcellular-location annotation, multi-label (a protein belongs to every compartment it lists; proteins
# without annotation are left out of the percentages). The panel is the cumulative view over all samples of a group
# (it describes the combined dataset).
###############################################################################
message("  Fig 2d: subcellular coverage (UniProt, cached)...")

UNIPROT_CACHE_FILE <- file.path(PROJECT_ROOT, "Data", "uniprot_subcellular_location_cache.csv")
UNIPROT_LOG_FILE   <- file.path(PROJECT_ROOT, "Data", "uniprot_query_log.csv")   # one row per request batch: date, UniProt release, counts
# UNIPROT_REFRESH: FALSE = use the cache and only ask UniProt for accessions not in it (same annotation every run);
# TRUE = ignore the cache and ask UniProt for ALL accessions again (the old cache is kept as *_before_refresh.csv).
# Set to TRUE once to rebuild the annotation from a single UniProt release, then set back to FALSE.
UNIPROT_REFRESH <- FALSE

sc_det <- rowSums(!is.na(protein_mat_log2[, thy1_singlecell_cols, drop=FALSE])) > 0
detected_singlecell <- rownames(protein_mat_log2)[sc_det]
bulk_det <- rowSums(!is.na(protein_mat_log2[, thy1_bulk10_cols, drop=FALSE])) > 0
detected_10x <- rownames(protein_mat_log2)[bulk_det]

# Look-up table gene symbol -> UniProt accession, taken from master_data. Protein.Group can list several
# accessions separated by ";", so only the first one is kept. One row per gene symbol.
gene_to_accession <- master_data %>%
  dplyr::select(gene_label, Protein.Group) %>%
  dplyr::mutate(accession = stringr::str_extract(Protein.Group, "^[^;]+")) %>%
  dplyr::select(gene_label, accession) %>%
  dplyr::distinct(gene_label, .keep_all = TRUE)

all_genes_needed <- union(detected_singlecell, detected_10x)
accessions_needed <- gene_to_accession %>%
  dplyr::filter(gene_label %in% all_genes_needed, !is.na(accession)) %>%
  dplyr::distinct(accession) %>% dplyr::pull(accession)

if (UNIPROT_REFRESH && file.exists(UNIPROT_CACHE_FILE)) {
  file.copy(UNIPROT_CACHE_FILE, sub("\\.csv$", "_before_refresh.csv", UNIPROT_CACHE_FILE), overwrite = TRUE)
  if (file.exists(UNIPROT_LOG_FILE)) file.copy(UNIPROT_LOG_FILE, sub("\\.csv$", "_before_refresh.csv", UNIPROT_LOG_FILE), overwrite = TRUE)
  file.remove(UNIPROT_CACHE_FILE); if (file.exists(UNIPROT_LOG_FILE)) file.remove(UNIPROT_LOG_FILE)
  message("  UNIPROT_REFRESH = TRUE: old cache kept as *_before_refresh.csv, all accessions are asked again")
}
if (file.exists(UNIPROT_CACHE_FILE)) {
  uniprot_cache_2d <- read_csv(UNIPROT_CACHE_FILE, show_col_types = FALSE,
                               col_types = cols(accession = col_character(), subcell_raw = col_character()))
  accessions_needed <- setdiff(accessions_needed, uniprot_cache_2d$accession)
} else {
  uniprot_cache_2d <- tibble::tibble(accession = character(), subcell_raw = character())
}


# query_uniprot_batch_2d(accs): asks UniProt for the subcellular-location comment of up to 100 accessions.
# Returns a tibble (accession, subcell_raw), or NULL if the request failed (network error, timeout,
# HTTP error or unexpected format); a failed batch is never written to the cache.
query_uniprot_batch_2d <- function(accs) {
  query <- paste0("accession:(", paste(accs, collapse = " OR "), ")")
  resp <- tryCatch(
    httr::GET("https://rest.uniprot.org/uniprotkb/search", httr::timeout(60),
              query = list(query = query, fields = "accession,cc_subcellular_location",
                           format = "tsv", size = length(accs))),
    error = function(e) e)
  if (inherits(resp, "error")) { warning("UniProt request failed: ", conditionMessage(resp), call. = FALSE); return(NULL) }
  if (httr::status_code(resp) != 200) { warning("UniProt returned HTTP ", httr::status_code(resp), call. = FALSE); return(NULL) }
  hdr <- function(h) { v <- httr::headers(resp)[[h]]; if (is.null(v)) NA_character_ else as.character(v) }
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  if (nchar(trimws(txt)) == 0) {
    df <- tibble::tibble(accession = character(), subcell_raw = character())
  } else {
    df <- readr::read_tsv(txt, show_col_types = FALSE)
    if (ncol(df) != 2) { warning("Unexpected UniProt response format (", ncol(df), " columns)", call. = FALSE); return(NULL) }
    colnames(df) <- c("accession", "subcell_raw")
    df$subcell_raw <- as.character(df$subcell_raw)
  }
  attr(df, "release") <- hdr("x-uniprot-release"); attr(df, "release_date") <- hdr("x-uniprot-release-date")
  df
}

failed_accessions_2d <- character()
if (length(accessions_needed) > 0) {
  check_packages("httr", "querying UniProt")
  batch_size <- 100
  batches <- split(accessions_needed, ceiling(seq_along(accessions_needed) / batch_size))
  message("  Querying UniProt in ", length(batches), " batch(es)...")
  query_log_2d <- list()
  new_results <- purrr::map_dfr(seq_along(batches), function(i) {
    res <- query_uniprot_batch_2d(batches[[i]])
    Sys.sleep(0.5)
    if (is.null(res)) { failed_accessions_2d <<- c(failed_accessions_2d, batches[[i]]); return(NULL) }
    rel <- attr(res, "release"); rel_date <- attr(res, "release_date")
    # A batch answer can be incomplete. Every accession missing from the answer is asked once more on its own;
    # only an accession that UniProt still does not return (obsolete or merged entry) is cached without an
    # annotation, so a partial answer can no longer silently turn proteins into "unannotated".
    unanswered <- setdiff(batches[[i]], res$accession)
    n_first <- nrow(res); recovered <- 0L; still_missing <- character()
    for (a in unanswered) {
      r <- query_uniprot_batch_2d(a); Sys.sleep(0.2)
      if (is.null(r)) failed_accessions_2d <<- c(failed_accessions_2d, a)       # request failed: not cached, retried next run
      else if (nrow(r) == 0) still_missing <- c(still_missing, a)              # UniProt really has no entry for it
      else { res <- dplyr::bind_rows(res, r); recovered <- recovered + 1L }
    }
    query_log_2d[[length(query_log_2d) + 1]] <<- tibble::tibble(
      date = as.character(Sys.Date()), uniprot_release = rel, uniprot_release_date = rel_date,
      n_requested = length(batches[[i]]), n_returned_first_try = n_first, n_recovered_by_single_query = recovered,
      n_not_returned = length(still_missing), n_with_location_text = sum(!is.na(res$subcell_raw) & nzchar(trimws(res$subcell_raw))))
    dplyr::bind_rows(res, tibble::tibble(accession = still_missing, subcell_raw = NA_character_))
  })
  uniprot_cache_2d <- dplyr::bind_rows(uniprot_cache_2d, new_results) %>% dplyr::distinct(accession, .keep_all = TRUE)
  write_csv(uniprot_cache_2d, UNIPROT_CACHE_FILE)
  new_log <- dplyr::bind_rows(query_log_2d)
  if (file.exists(UNIPROT_LOG_FILE)) new_log <- dplyr::bind_rows(read_csv(UNIPROT_LOG_FILE, show_col_types = FALSE, col_types = cols(.default = col_character())) %>%
                                                                   dplyr::mutate(dplyr::across(dplyr::starts_with("n_"), as.numeric)), new_log)
  write_csv(new_log, UNIPROT_LOG_FILE)
  writeLines(paste0("UniProt subcellular-location cache last extended on ", Sys.Date(), " (", nrow(new_results),
                    " accessions added; ", nrow(uniprot_cache_2d), " in total). UniProt release(s) of this extension: ",
                    paste(unique(stats::na.omit(dplyr::bind_rows(query_log_2d)$uniprot_release)), collapse = ", "),
                    ". Per-batch details: uniprot_query_log.csv"),
             file.path(DIR_F2, "uniprot_subcellular_location_cache_info.txt"))
}
if (length(failed_accessions_2d) > 0) {
  warning(length(failed_accessions_2d), " accession(s) could not be retrieved from UniProt and are counted as ",
          "unannotated (left out of the percentages) in this run (they are not cached; run the script again to retry).", call. = FALSE)
}
n_no_accession_2d <- sum(is.na(gene_to_accession$accession[gene_to_accession$gene_label %in% all_genes_needed]))
if (n_no_accession_2d > 0) message("  ", n_no_accession_2d, " detected protein group(s) have no accession and are unannotated (left out of the percentages)")

# Compartment annotation, MULTI-LABEL. Every location that UniProt lists for a protein is classified into one
# compartment by keyword (within one listed location, in this priority order: mitochondri -> Mitochondria;
# nucleus/nuclear -> Nuclear; cytoskeleton -> Cytoskeleton; membrane -> Membrane; cytoplasm/cytosol -> Cytosolic;
# anything else -> "Other (annotated)"), and a protein belongs to EVERY compartment it lists. Text after ";" in a
# location (e.g. "Multi-pass membrane protein") and "Note=" text are not locations and are ignored.
# Proteins WITHOUT any UniProt subcellular-location annotation (no accession, accession not returned, no comment)
# are not assigned and are left out of the denominator: percentages are "% of annotated protein groups". They can add
# up to more than 100% (a protein can be in several compartments) and describe annotation, not measured localisation.
COMPARTMENTS <- c("Mitochondria", "Nuclear", "Cytoskeleton", "Membrane", "Cytosolic", "Other (annotated)")
parse_locations <- function(text) {
  t <- stringr::str_remove(text, "^SUBCELLULAR LOCATION:\\s*")
  t <- stringr::str_remove(t, "\\s*Note=.*$")
  t <- stringr::str_remove_all(t, "\\{[^}]*\\}")                 # evidence tags
  t <- stringr::str_remove_all(t, "\\[[^]]*\\]:\\s*")            # isoform prefixes
  entries <- stringr::str_trim(unlist(stringr::str_split(t, "\\.\\s*")))
  entries <- stringr::str_trim(stringr::str_remove(entries, ";.*$"))
  entries[entries != ""]
}
classify_location <- function(entry) {
  t <- stringr::str_to_lower(entry)
  dplyr::case_when(
    stringr::str_detect(t, "mitochondri") ~ "Mitochondria",
    stringr::str_detect(t, "nucleus|nuclear") ~ "Nuclear",
    stringr::str_detect(t, "cytoskeleton") ~ "Cytoskeleton",
    stringr::str_detect(t, "membrane") ~ "Membrane",
    stringr::str_detect(t, "cytoplasm|cytosol") ~ "Cytosolic",
    TRUE ~ "Other (annotated)")
}
# compartments_from_text(): NULL = no annotation; otherwise the set of compartments the protein lists
compartments_from_text <- function(text) {
  if (is.na(text) || !nzchar(trimws(text))) return(NULL)
  entries <- parse_locations(text)
  if (length(entries) == 0) return("Other (annotated)")
  unique(classify_location(entries))
}
uniprot_comp_list <- setNames(purrr::map(uniprot_cache_2d$subcell_raw, compartments_from_text), uniprot_cache_2d$accession)
gene_acc       <- gene_to_accession$accession[match(rownames(protein_mat_log2), gene_to_accession$gene_label)]
gene_comps     <- uniprot_comp_list[gene_acc]
gene_annotated <- !vapply(gene_comps, is.null, logical(1)); names(gene_annotated) <- rownames(protein_mat_log2)
gene_comp_mat  <- t(vapply(gene_comps, function(cl) COMPARTMENTS %in% cl, logical(length(COMPARTMENTS))))   # proteins x compartments (TRUE/FALSE)
dimnames(gene_comp_mat) <- list(rownames(protein_mat_log2), COMPARTMENTS)

# Cumulative view: all protein groups detected in >= 1 sample of the group (describes the combined dataset, not a typical cell)
summarise_compartments <- function(genes, label) {
  annotated <- genes[gene_annotated[genes]]
  tibble::tibble(bucket = COMPARTMENTS, n_proteins = as.integer(colSums(gene_comp_mat[annotated, , drop = FALSE])),
                 n_annotated = length(annotated), n_detected = length(genes), group = label) %>%
    dplyr::mutate(pct = 100 * n_proteins / n_annotated)
}
sc_summary_2d   <- summarise_compartments(detected_singlecell, "Single cells (cumulative)")
bulk_summary_2d <- summarise_compartments(detected_10x, "10-cell area (cumulative)")
combined_2d <- dplyr::bind_rows(sc_summary_2d, bulk_summary_2d)
write_csv(combined_2d, file.path(DIR_F2, "2d_subcellular_classification_summary.csv"))

bucket_order <- sc_summary_2d %>% dplyr::arrange(dplyr::desc(pct)) %>% dplyr::pull(bucket)
bucket_order <- c(setdiff(bucket_order, "Other (annotated)"), "Other (annotated)")          # "Other" always last
combined_2d <- combined_2d %>% dplyr::mutate(bucket = factor(bucket, levels = rev(bucket_order)))
message("  Subcellular coverage (multi-label): 10x area ", length(detected_10x), " protein groups (", bulk_summary_2d$n_annotated[1],
        " annotated) across ", length(thy1_bulk10_cols), " wells  |  single cells ", length(detected_singlecell), " protein groups (",
        sc_summary_2d$n_annotated[1], " annotated) across ", length(thy1_singlecell_cols), " cells")

fig2d_compartments <- ggplot(combined_2d, aes(x = bucket, y = pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Single cells (cumulative)" = PAL_THY1, "10-cell area (cumulative)" = PAL_BULK),
                    labels = c("Single cells (cumulative)" = "Single cells", "10-cell area (cumulative)" = "10-cell areas"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +   # bars start exactly at 0 (this is the axis coord_flip puts at the bottom)
  coord_flip() +
  labs(x=NULL, y="Annotated protein groups (%)") +
  THEME_JOURNAL() +
  theme(legend.position = c(0.98, 0.25), legend.justification = c(1, 0.5),
        legend.background = element_blank(), legend.key.size = grid::unit(2.5, "mm")) +
  guides(fill = guide_legend(reverse = TRUE))

###############################################################################
# 6. FIG. 2e: PROTEIN OVERLAP (VENN)
###############################################################################
message("  Fig 2e: tier overlap (Venn diagram)...")
detected_set <- function(cols) rownames(protein_mat_log2)[rowSums(!is.na(protein_mat_log2[, cols, drop=FALSE])) > 0]

set_1x  <- detected_set(thy1_singlecell_cols)
set_10x <- detected_set(thy1_bulk10_cols)


circle_overlap_area <- function(d, r1, r2) {
  if (d >= r1 + r2) return(0)
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1)) +
    r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2)) -
    0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
}

n_10x_only_2e   <- length(setdiff(set_10x, set_1x))
n_1x_only_2e    <- length(setdiff(set_1x, set_10x))
n_overlap_2e    <- length(intersect(set_10x, set_1x))
n_total_10x_2e  <- n_10x_only_2e + n_overlap_2e
n_total_1x_2e   <- n_1x_only_2e + n_overlap_2e
n_union_2e      <- n_10x_only_2e + n_overlap_2e + n_1x_only_2e


if (n_1x_only_2e + n_overlap_2e == 0 || n_10x_only_2e + n_overlap_2e == 0) {
  stop("Venn diagram: no protein group detected in the single cells or in the 10-cell areas.", call. = FALSE)
}
stopifnot(n_union_2e == length(union(set_1x, set_10x)),
          n_total_1x_2e == length(set_1x), n_total_10x_2e == length(set_10x))   # counts add up to the sets

r_10x_2e <- sqrt(n_total_10x_2e / pi)
r_1x_2e  <- sqrt(n_total_1x_2e / pi)


# Distance between the circle centres so that the overlap area equals the shared count; the two
# degenerate cases (no overlap, or one set completely inside the other) have no root to solve for.
d_solved_2e <- if (n_overlap_2e == 0) {
  r_10x_2e + r_1x_2e
} else if (n_overlap_2e == min(n_total_10x_2e, n_total_1x_2e)) {
  abs(r_10x_2e - r_1x_2e)
} else {
  uniroot(function(d) circle_overlap_area(d, r_10x_2e, r_1x_2e) - n_overlap_2e,
          interval = c(abs(r_10x_2e - r_1x_2e) + 1e-6, r_10x_2e + r_1x_2e - 1e-6))$root
}
if (abs(circle_overlap_area(d_solved_2e, r_10x_2e, r_1x_2e) - n_overlap_2e) > 0.5) {
  stop("Venn diagram: the drawn overlap does not match the shared count.", call. = FALSE)
}

cy_10x_2e <- d_solved_2e / 2
cy_1x_2e  <- -d_solved_2e / 2

make_circle_df <- function(cx, cy, r, n = 200) {
  theta <- seq(0, 2 * pi, length.out = n)
  tibble::tibble(x = cx + r * cos(theta), y = cy + r * sin(theta))
}
circle_10x_df_2e <- make_circle_df(0, cy_10x_2e, r_10x_2e)
circle_1x_df_2e  <- make_circle_df(0, cy_1x_2e, r_1x_2e)

pct_10x_only_2e <- sprintf("%.1f", 100 * n_10x_only_2e / n_union_2e)   # one decimal so the three shares add to 100.0%
pct_1x_only_2e  <- sprintf("%.1f", 100 * n_1x_only_2e / n_union_2e)
pct_overlap_2e  <- sprintf("%.1f", 100 * n_overlap_2e / n_union_2e)


top_edge_10x_2e    <- cy_10x_2e + r_10x_2e
bottom_edge_10x_2e <- cy_10x_2e - r_10x_2e
top_edge_1x_2e     <- cy_1x_2e + r_1x_2e
bottom_edge_1x_2e  <- cy_1x_2e - r_1x_2e

y_10x_only_label_2e <- (top_edge_10x_2e + top_edge_1x_2e) / 2
y_1x_only_label_2e  <- (bottom_edge_1x_2e + bottom_edge_10x_2e) / 2
y_overlap_label_2e  <- (max(bottom_edge_10x_2e, bottom_edge_1x_2e) + min(top_edge_10x_2e, top_edge_1x_2e)) / 2

# Only the "shared" number sits inside (the two crescents are too thin to hold
# text without touching the circle outline). The two exclusive counts go OUTSIDE,
# directly under/over the set name, so nothing overlaps a line. All text 7 pt.
fmt_n <- function(n) format(n, big.mark = ",")
shared_label_2e <- tibble::tibble(
  x = 0, y = y_overlap_label_2e,
  label = paste0("Shared\n", fmt_n(n_overlap_2e), " (", pct_overlap_2e, "%)"))
setname_df_2e <- tibble::tibble(
  x = 0,
  y = c(top_edge_10x_2e + 0.05 * r_10x_2e, bottom_edge_1x_2e - 0.05 * r_1x_2e),
  vj = c(0, 1),
  label = c(paste0("10-cell areas only\n", fmt_n(n_10x_only_2e), " (", pct_10x_only_2e, "%)"),
            paste0("Single cells only\n", fmt_n(n_1x_only_2e), " (", pct_1x_only_2e, "%)")))

fig2e_venn <- ggplot2::ggplot() +
  ggplot2::geom_polygon(data = circle_10x_df_2e, ggplot2::aes(x, y),
                        fill = PAL_BULK, alpha = 0.35, colour = "grey30", linewidth = 0.5) +
  ggplot2::geom_polygon(data = circle_1x_df_2e, ggplot2::aes(x, y),
                        fill = PAL_THY1, alpha = 0.35, colour = "grey30", linewidth = 0.5) +
  ggplot2::geom_text(data = shared_label_2e, ggplot2::aes(x, y, label = label),
                     size = PT_BASE, size.unit = "pt", lineheight = 0.95) +
  ggplot2::geom_text(data = setname_df_2e, ggplot2::aes(x, y, label = label, vjust = vj),
                     size = PT_BASE, size.unit = "pt", lineheight = 0.95) +
  ggplot2::expand_limits(y = c(bottom_edge_1x_2e - 0.60 * r_1x_2e, top_edge_10x_2e + 0.60 * r_10x_2e)) +
  ggplot2::coord_equal() +
  ggplot2::theme_void()


set_size_summary_2e <- tibble::tibble(tier = c("Single cells", "10-cell area"),
                                      n_proteins_detected = c(length(set_1x), length(set_10x)))
write_csv(set_size_summary_2e, file.path(DIR_F2, "2e_tier_overlap_set_sizes.csv"))
message("  Tier overlap: 1x n=", length(set_1x), "  10x area n=", length(set_10x))

message("  Figure 2 main panels (a-e) saved to ", DIR_F2)

###############################################################################
# 7. MAIN FIGURE 2: ASSEMBLY
# Row 1 = protein count | peak capacity | input-amount titration
# Row 2 = sensitivity | subcellular compartments | Venn diagram
###############################################################################
# Bold category labels on all category axes.
CAT_X_BOLD <- ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold"))
CAT_Y_BOLD <- ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold"))
fig2a_columns              <- fig2a_columns + CAT_X_BOLD
fig2a_peakcap <- fig2a_peakcap + CAT_X_BOLD
fig2c_sensitivity <- fig2c_sensitivity + CAT_X_BOLD
fig2d_compartments <- fig2d_compartments + CAT_Y_BOLD   # compartment names (axis is flipped)
fig2_row1 <- (fig2a_columns | fig2a_peakcap | fig2b_titration) +
  patchwork::plot_layout(widths = rep(1, 3))
fig2_row2 <- (fig2c_sensitivity | fig2d_compartments | fig2e_venn) +
  patchwork::plot_layout(widths = c(1, 1, 1))

fig2_main_assembled <- (fig2_row1 / fig2_row2) +
  patchwork::plot_layout(heights = c(1.2, 1.2))


FIG2_EXPORT_WIDTH_IN  <- FIG2_W_IN
FIG2_EXPORT_HEIGHT_IN <- 2 * FIG2_ROW_H_IN

source_panel("Fig2a_protein_groups", fig2a_columns, "Protein groups (thousands) per HeLa injection at 1 cell equivalent, HT and Neo Plus configuration; mean and SD per configuration",
             columns = c("File.Name", "column", "spd", "k_proteins", "mean_k", "sd_k", "n"), table_names = c("mean_sd", "per_inj"))
source_panel("Fig2a_peak_capacity", fig2a_peakcap, "Peak capacity per HeLa injection, HT and Neo Plus configuration; mean and SD per configuration",
             columns = c("File.Name", "column", "spd", "gradient_time_min", "FWHM.RT", "peak_capacity", "mean_peakcap", "sd_peakcap", "n"), table_names = c("mean_sd", "per_inj"))
source_panel("Fig2b_titration", fig2b_titration, "Protein groups (thousands) per HeLa injection across 0.1-100 cell equivalents (Neo Plus); mean and SD per input",
             columns = c("File.Name", "load_group", "cell_equiv", "k_proteins", "mean_k", "sd_k", "n"), table_names = c("mean_sd", "per_inj"))
source_panel("Fig2c_sensitivity", fig2c_sensitivity, "Protein groups per Thy1 single cell and per 10-cell area (samples with >= 200 protein groups)",
             columns = c("sample_col", "tier_label", "batch", "batch_label", "n_proteins"))
source_panel("Fig2d_compartments", fig2d_compartments, "Subcellular compartment coverage (percent of annotated protein groups), cumulative over single cells and 10-cell areas")
source_panel("Fig2e_venn_counts", data.frame(set = c("shared", "10-cell areas only", "single cells only", "10-cell areas total", "single cells total"),
                                             protein_groups = c(n_overlap_2e, n_total_10x_2e - n_overlap_2e, n_total_1x_2e - n_overlap_2e, n_total_10x_2e, n_total_1x_2e)),
             "Counts behind the Venn diagram")
source_panel("Fig2e_venn_protein_groups", data.frame(gene = union(detected_singlecell, detected_10x),
                                                     single_cells = union(detected_singlecell, detected_10x) %in% detected_singlecell,
                                                     ten_cell_areas = union(detected_singlecell, detected_10x) %in% detected_10x),
             "Protein groups detected in at least one single cell and/or at least one 10-cell area")
ggsave(file.path(DIR_F2, "Figure2_MAIN_assembled.svg"), fig2_main_assembled,
       width = FIG2_EXPORT_WIDTH_IN, height = FIG2_EXPORT_HEIGHT_IN,
       units = "in", device = SVG_DEVICE)
message("  Figure 2 assembled layout saved as SVG (", FIG2_EXPORT_WIDTH_IN, " x ",
        FIG2_EXPORT_HEIGHT_IN, " in) -- insert at 100%, do not rescale.")

###############################################################################
# 8. EXTENDED DATA FIGURE 3: method reliability and instrument QC
###############################################################################
message("\n[ED Fig. 3] Method reliability & instrument QC...")

message("  ED Fig 3 panel B: batch-to-batch single-cell-average concordance (Thy1, per region)...")

MIN_FRAC_DETECTED_BATCH_PB <- 0.5
# Cells need >= MIN_PROTEINS (600) proteins here: this panel reports Pearson/Spearman correlations,
# which need the stricter bar used for the analysis set, not the 200-protein bar of the descriptive panels.
build_batch_single_cell_average <- function(cols) {
  if (length(cols) == 0) return(NULL)
  sub <- protein_mat_log2[, cols, drop = FALSE]
  frac_det <- rowMeans(!is.na(sub))
  keep <- frac_det >= MIN_FRAC_DETECTED_BATCH_PB
  rowMeans(sub[keep, , drop = FALSE], na.rm = TRUE)
}

thy1_batches_present <- cohort_meta %>% dplyr::filter(tier=="tier_1x_single_cell", cell_type_label=="Thy1") %>%
  dplyr::distinct(batch) %>% dplyr::pull(batch)
batch_pairs <- if (length(thy1_batches_present) >= 2) utils::combn(sort(thy1_batches_present), 2, simplify = FALSE) else list()

batch_corr_summary <- list()
for (region_val in REGIONS_THY1) {
  for (pair in batch_pairs) {
    b1 <- pair[1]; b2 <- pair[2]
    cols1 <- cohort_meta %>% dplyr::filter(tier=="tier_1x_single_cell", cell_type_label=="Thy1",
                                           region==region_val, batch==b1, n_proteins_lookup[sample_col] >= MIN_PROTEINS) %>% dplyr::pull(sample_col)
    cols2 <- cohort_meta %>% dplyr::filter(tier=="tier_1x_single_cell", cell_type_label=="Thy1",
                                           region==region_val, batch==b2, n_proteins_lookup[sample_col] >= MIN_PROTEINS) %>% dplyr::pull(sample_col)
    pb1 <- build_batch_single_cell_average(cols1); pb2 <- build_batch_single_cell_average(cols2)
    if (is.null(pb1) || is.null(pb2)) {
      message("  ED 3b: ", REGION_ABBREV[[region_val]], " ", b1, " vs ", b2, " skipped (no cells >= ", MIN_PROTEINS, " in one batch)")
      next
    }
    common_g <- intersect(names(pb1), names(pb2))
    if (length(common_g) < 3) {
      message("  ED 3b: ", REGION_ABBREV[[region_val]], " ", b1, " vs ", b2, " skipped (fewer than 3 shared proteins)")
      next
    }
    x <- pb1[common_g]; y <- pb2[common_g]
    if (sd(x) == 0 || sd(y) == 0) {
      message("  ED 3b: ", REGION_ABBREV[[region_val]], " ", b1, " vs ", b2, " skipped (constant values; correlation undefined)")
      next
    }
    r_p <- cor(x, y); r_s <- cor(x, y, method = "spearman")
    key <- paste0(region_val, "__", b1, "_vs_", b2)
    b1_label <- BATCH_DISPLAY_LABELS[[b1]]; b2_label <- BATCH_DISPLAY_LABELS[[b2]]
    # Fisher z transform for a closed-form 95% CI on Pearson r -- cheap,
    # appropriate for this sample size range, no bootstrap needed.
    z <- atanh(r_p); se_z <- 1 / sqrt(length(common_g) - 3)
    ci_lo <- tanh(z - 1.96 * se_z); ci_hi <- tanh(z + 1.96 * se_z)
    batch_corr_summary[[key]] <- tibble::tibble(region = region_val, batch_1 = b1, batch_2 = b2,
                                                comparison_label = paste0(REGION_ABBREV[[region_val]], ": ", b1_label, " vs ", b2_label),
                                                pearson_r = r_p, pearson_ci_lo = ci_lo, pearson_ci_hi = ci_hi,
                                                spearman_r = r_s, n = length(common_g),
                                                n_cells_1 = length(cols1), n_cells_2 = length(cols2))
  }
}
if (length(batch_corr_summary) > 0) {
  region_order <- dplyr::bind_rows(batch_corr_summary) %>%
    dplyr::group_by(region) %>% dplyr::summarise(mean_r = mean(pearson_r), .groups = "drop") %>%
    dplyr::arrange(mean_r) %>% dplyr::pull(region)
  summary_df <- dplyr::bind_rows(batch_corr_summary) %>%
    dplyr::mutate(region = factor(region, levels = region_order)) %>%
    dplyr::arrange(region, pearson_r) %>%
    dplyr::mutate(comparison_label = factor(comparison_label, levels = comparison_label))
  
  edfig3b_concordance <- ggplot(summary_df, aes(pearson_r, comparison_label, colour = region)) +
    geom_pointrange(aes(xmin = pearson_ci_lo, xmax = pearson_ci_hi), size = 0.35, fatten = 3) +
    geom_text(aes(x = pearson_ci_hi + 0.018, label = paste0("n = ", format(n, big.mark = ","))),
              hjust = 0, size = pt_to_geom(PT_MIN), colour = "grey35") +
    scale_colour_manual(values = PAL_REGION_THY1, guide = "none") +
    coord_cartesian(xlim = c(min(summary_df$pearson_ci_lo) - 0.02, 1.06)) +
    labs(x = "Pearson r (95% CI); n = protein groups", y = NULL) +
    THEME_JOURNAL() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  
  write_csv(summary_df, file.path(DIR_F2_SUPP, "EDFig3_b_batch_concordance_summary.csv"))
  message("  ED Fig 3 panel B saved (", length(batch_corr_summary), " region x batch-pair comparisons)")
} else {
  message("  NOTE: ED Fig 3 panel B skipped -- fewer than 2 Thy1 batches with QC-passed cells found.")
  edfig3b_concordance <- NULL
}

# ---- ED Fig. 3b control (table only, the figure is unchanged) ------------------------------------------
# Question: is the batch-to-batch correlation higher for the SAME region than for a DIFFERENT region?
# Every pair of (region, batch) groups from different batches is compared. If matched-region pairs are not
# clearly higher than different-region pairs, the correlation only shows that neurons look alike.
# pearson_r_centred repeats the comparison after subtracting each protein's average over all groups, so
# only the differences between groups remain (proteins detected in every group only).
groups_pb <- cohort_meta %>%
  dplyr::filter(tier == "tier_1x_single_cell", cell_type_label == "Thy1", region %in% REGIONS_THY1,
                n_proteins_lookup[sample_col] >= MIN_PROTEINS) %>%
  dplyr::group_by(region, batch) %>%
  dplyr::summarise(cols = list(sample_col), n_cells = dplyr::n(), .groups = "drop")
groups_pb$avg <- lapply(groups_pb$cols, build_batch_single_cell_average)
common_all <- Reduce(intersect, lapply(groups_pb$avg, names))
centred <- NULL
if (length(common_all) >= 10) {
  A <- sapply(groups_pb$avg, function(v) v[common_all]); centred <- A - rowMeans(A)
}
rows_cmp <- list()
for (i in seq_len(nrow(groups_pb))) for (j in seq_len(nrow(groups_pb))) {
  if (i >= j || groups_pb$batch[i] == groups_pb$batch[j]) next
  a <- groups_pb$avg[[i]]; b <- groups_pb$avg[[j]]
  g <- intersect(names(a), names(b)); if (length(g) < 3) next
  x <- a[g]; y <- b[g]; if (sd(x) == 0 || sd(y) == 0) next
  rows_cmp[[length(rows_cmp) + 1]] <- tibble::tibble(
    pair_type = if (groups_pb$region[i] == groups_pb$region[j]) "same region" else "different region",
    region_1 = REGION_ABBREV[[groups_pb$region[i]]], batch_1 = BATCH_DISPLAY_LABELS[[groups_pb$batch[i]]], n_cells_1 = groups_pb$n_cells[i],
    region_2 = REGION_ABBREV[[groups_pb$region[j]]], batch_2 = BATCH_DISPLAY_LABELS[[groups_pb$batch[j]]], n_cells_2 = groups_pb$n_cells[j],
    n_protein_groups = length(g), pearson_r = cor(x, y), spearman_r = cor(x, y, method = "spearman"),
    pearson_r_centred = if (is.null(centred)) NA_real_ else cor(centred[, i], centred[, j]),
    n_protein_groups_centred = if (is.null(centred)) NA_integer_ else length(common_all))
}
if (length(rows_cmp) > 0) {
  region_pair_cmp <- dplyr::bind_rows(rows_cmp) %>% dplyr::arrange(pair_type, dplyr::desc(pearson_r))
  write_csv(region_pair_cmp, file.path(DIR_F2_SUPP, "EDFig3_b_same_vs_different_region_correlations.csv"))
  cmp_sum <- region_pair_cmp %>% dplyr::group_by(pair_type) %>%
    dplyr::summarise(pairs = dplyr::n(), min_r = min(pearson_r), median_r = median(pearson_r), max_r = max(pearson_r),
                     median_r_centred = median(pearson_r_centred, na.rm = TRUE), .groups = "drop")
  message("  ED 3b control: same-region vs different-region correlations across batches")
  print(as.data.frame(cmp_sum), digits = 3)
} else {
  message("  NOTE: ED 3b control skipped -- no usable pairs of groups from different batches.")
}

batch_level_df <- cohort_meta %>%
  dplyr::filter(tier == "tier_1x_single_cell", cell_type_label == "Thy1", passes_qc_tier) %>%
  dplyr::mutate(n_proteins_filtered = as.numeric(n_proteins_lookup[sample_col])) %>%
  dplyr::filter(is.finite(n_proteins_filtered)) %>%
  dplyr::mutate(batch = factor(batch, levels = intersect(names(BATCH_DISPLAY_LABELS), unique(batch))))
message("  Cross-batch yield level (Thy1 single cells only): n=", nrow(batch_level_df), " across ",
        dplyr::n_distinct(batch_level_df$batch), " batches")
batch_medians <- batch_level_df %>% dplyr::group_by(batch) %>%
  dplyr::summarise(median_proteins = median(n_proteins_filtered), n = dplyr::n(), .groups = "drop")
print(batch_medians)


batch_fold_diff <- max(batch_medians$median_proteins) / min(batch_medians$median_proteins)
message("  Batch medians ranged from ", min(batch_medians$median_proteins), " to ",
        max(batch_medians$median_proteins), " (", round(batch_fold_diff, 2), "-fold difference)")
# No statistical test on this panel: descriptive only (medians and fold difference above).

edfig3a_yield_batch <- ggplot(batch_level_df, aes(batch, n_proteins_filtered, fill = batch)) +
  ggdist::stat_halfeye(
    adjust = 0.5, width = 0.45, .width = 0,
    justification = -0.56, point_colour = NA, alpha = 0.6
  ) +
  geom_boxplot(width = 0.42, outlier.shape = NA, alpha = 0.85) +
  geom_point(
    size = 1, alpha = 0.25,
    position = position_jitter(seed = 1, width = 0.16)
  ) +
  scale_fill_manual(values = setNames(rep(PAL_THY1, nlevels(batch_level_df$batch)), levels(batch_level_df$batch)), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +   # raincloud starts exactly at 0
  coord_cartesian(xlim = c(0.45, nlevels(batch_level_df$batch) + 0.45), clip = "off") +
  scale_x_discrete(labels = setNames(paste0(BATCH_DISPLAY_LABELS[as.character(batch_medians$batch)], "\n(n = ", batch_medians$n, ")"),
                                     as.character(batch_medians$batch))) +
  labs(x = NULL, y = "Protein groups detected") +
  THEME_JOURNAL()


if (!"area_um2_final" %in% colnames(cohort_meta)) {
  message("  SKIPPED: cell-area/yield-by-type panel -- 'area_um2_final' not found in cohort_meta/sample_metadata. ",
          "Confirm master_sample_metadata.csv includes this column (see 00_build_master_dataset.R README).")
  edfig3c_area <- edfig3c_yield <- NULL
} else {
  area_yield_df <- cohort_meta %>%
    dplyr::filter(tier == "tier_1x_single_cell", cell_type_label %in% c("Thy1", "Microglia"),
                  passes_qc_tier, is.finite(area_um2_final)) %>%
    dplyr::mutate(n_proteins_filtered = as.numeric(n_proteins_lookup[sample_col])) %>%
    dplyr::filter(is.finite(n_proteins_filtered)) %>%
    dplyr::mutate(cell_type_display = factor(dplyr::if_else(cell_type_label == "Thy1", "Thy1-EGFP neurons", "Cx3cr1-EGFP microglia"),
                                             levels = c("Thy1-EGFP neurons", "Cx3cr1-EGFP microglia")))
  n_by_type <- area_yield_df %>% dplyr::count(cell_type_display)
  type_n_labels <- setNames(paste0(c("Thy1-EGFP\nneurons", "Cx3cr1-EGFP\nmicroglia")[match(n_by_type$cell_type_display, c("Thy1-EGFP neurons", "Cx3cr1-EGFP microglia"))], "\n(n = ", n_by_type$n, ")"), n_by_type$cell_type_display)
  message("  Area/yield by cell type: n=", nrow(area_yield_df), " QC-passing single cells (",
          paste(sprintf("%s=%d", n_by_type$cell_type_display, n_by_type$n), collapse=", "), ")")
  
  # Area is reported descriptively only 
  area_yield_medians <- area_yield_df %>% dplyr::group_by(cell_type_display) %>%
    dplyr::summarise(median_area = median(area_um2_final), median_yield = median(n_proteins_filtered), .groups = "drop")
  print(area_yield_medians)
  message("  Area (descriptive only, no test): Thy1+ median=",
          round(area_yield_medians$median_area[area_yield_medians$cell_type_display=="Thy1-EGFP neurons"],1),
          "  Cx3cr1+ median=", round(area_yield_medians$median_area[area_yield_medians$cell_type_display=="Cx3cr1-EGFP microglia"],1))
  message("  Yield (descriptive only, no test): Thy1+ median=", area_yield_medians$median_yield[area_yield_medians$cell_type_display=="Thy1-EGFP neurons"],
          ", Cx3cr1+ median=", area_yield_medians$median_yield[area_yield_medians$cell_type_display=="Cx3cr1-EGFP microglia"])
  
  edfig3c_area <- ggplot(area_yield_df, aes(cell_type_display, area_um2_final, fill = cell_type_display)) +
    ggdist::stat_halfeye(
      adjust = 0.5, width = 0.45, .width = 0,
      justification = -0.56, point_colour = NA, alpha = 0.6
    ) +
    geom_boxplot(width = 0.42, outlier.shape = NA, alpha = 0.85) +
    geom_point(
      size = 1, alpha = 0.25,
      position = position_jitter(seed = 1, width = 0.16)
    ) +
    scale_fill_manual(values = c(`Thy1-EGFP neurons` = PAL_THY1, `Cx3cr1-EGFP microglia` = PAL_MICROGLIA), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +   # raincloud starts exactly at 0
    coord_cartesian(xlim = c(0.5, 2.5), clip = "off") +
    scale_x_discrete(labels = type_n_labels) +
    labs(x = NULL, y = "LMD collection area (\u00b5m\u00b2)") +
    THEME_JOURNAL()
  
  edfig3c_yield <- ggplot(area_yield_df, aes(cell_type_display, n_proteins_filtered, fill = cell_type_display)) +
    ggdist::stat_halfeye(
      adjust = 0.5, width = 0.45, .width = 0,
      justification = -0.56, point_colour = NA, alpha = 0.6
    ) +
    geom_boxplot(width = 0.42, outlier.shape = NA, alpha = 0.85) +
    geom_point(
      size = 1, alpha = 0.25,
      position = position_jitter(seed = 1, width = 0.16)
    ) +
    scale_fill_manual(values = c(`Thy1-EGFP neurons` = PAL_THY1, `Cx3cr1-EGFP microglia` = PAL_MICROGLIA), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +   # raincloud starts exactly at 0
    coord_cartesian(xlim = c(0.5, 2.5), clip = "off") +
    scale_x_discrete(labels = type_n_labels) +
    labs(x = NULL, y = "Protein groups detected") +
    THEME_JOURNAL()
  
}


message("  ED Fig 3 panel D: log-log scaling exponent (input amount vs. abundance) -- precursor-level, NOT MaxLFQ...")

if (is.null(precursor_mat_log2)) {
  message("  SKIPPED: precursor matrix not loaded -- see warning above.")
  edfig3d_scaling <- NULL
} else {
  scaling_1x_ids  <- sample_id_lookup %>% dplyr::filter(tier=="tier_1x_single_cell", cell_type_label=="Thy1", passes_qc_tier) %>% dplyr::pull(sample_id)   # same >=200 cells as Fig 2c (ED 3d is the intensity companion of 2c)
  scaling_10x_ids <- sample_id_lookup %>% dplyr::filter(tier=="tier_10x_area", passes_qc_tier) %>% dplyr::pull(sample_id)
  scaling_1x_ids  <- intersect(scaling_1x_ids, colnames(precursor_mat_log2))
  scaling_10x_ids <- intersect(scaling_10x_ids, colnames(precursor_mat_log2))
  message("  Precursor-level scaling: ", length(scaling_1x_ids), " 1x samples, ", length(scaling_10x_ids), " 10x samples matched")
  
  # Same minimum for BOTH sample types: >= MIN_PROTEINS_VIABLE (200) contaminant-filtered protein groups
  # (passes_qc_tier), i.e. the same samples as Fig. 2c. Checked, not assumed.
  samples_used <- sample_id_lookup %>% dplyr::filter(sample_id %in% c(scaling_1x_ids, scaling_10x_ids))
  stopifnot(nrow(samples_used) == length(scaling_1x_ids) + length(scaling_10x_ids),
            all(n_proteins_lookup[samples_used$clean_sample_name] >= MIN_PROTEINS_VIABLE))
  message("  CHECK passed: all ", nrow(samples_used), " samples of the scaling analysis (single cells and 10-cell areas) have >= ",
          MIN_PROTEINS_VIABLE, " protein groups")
  
  MIN_FRAC_DETECTED_SCALING <- 0.8
  # scaling_table(pr_mat): for every precursor detected in >= 80% of the samples of BOTH groups: the mean log2 intensity
  # in each group (over the samples where it was detected; nothing imputed), fold change = 2^(mean log2 in 10-cell
  # areas - mean log2 in single cells), exponent beta = (that difference) / log2(10). beta = 1 would be proportional to the
  # nominal 10-fold difference in cell number. fold_change is computed directly (not from beta) so that its mean is a real mean.
  scaling_table <- function(pr_mat) {
    m1  <- pr_mat[, scaling_1x_ids, drop = FALSE]
    m10 <- pr_mat[, scaling_10x_ids, drop = FALSE]
    keep <- rowMeans(!is.na(m1)) >= MIN_FRAC_DETECTED_SCALING & rowMeans(!is.na(m10)) >= MIN_FRAC_DETECTED_SCALING
    if (sum(keep) < 10) return(NULL)
    a <- rowMeans(m1[keep, , drop = FALSE], na.rm = TRUE)
    b <- rowMeans(m10[keep, , drop = FALSE], na.rm = TRUE)
    tibble::tibble(precursor = names(a), mean_log2_1x = a, mean_log2_10x = b,
                   fold_change = 2^(b - a), beta = (b - a) / log2(10)) %>% dplyr::filter(is.finite(beta))
  }
  scaling_df <- scaling_table(precursor_mat_log2)
  if (is.null(scaling_df)) {
    message("  NOTE: fewer than 10 precursors detected in >=", MIN_FRAC_DETECTED_SCALING*100,
            "% of both 1x and 10x samples — ED Fig 3d skipped.")
    edfig3d_scaling <- NULL
  } else {
    mean_beta <- mean(scaling_df$beta)
    median_beta <- median(scaling_df$beta)
    mean_fold_change <- mean(scaling_df$fold_change)
    median_fold_change <- median(scaling_df$fold_change)
    write_csv(scaling_df, file.path(DIR_F2, "EDFig3_d_scaling_exponent_per_precursor.csv"))
    message("  Scaling exponent (precursor-level): n=", nrow(scaling_df), " precursors, mean beta=", round(mean_beta,3),
            ", median beta=", round(median_beta,3))
    message("  Fold-change (precursor-level, 10-cell area vs single cells): mean=", round(mean_fold_change,3),
            "-fold, median=", round(median_fold_change,3), "-fold")
    # Effect of the duplicate-precursor rule (only exists if duplicate records carried different quantities)
    scaling_rule_check <- NULL
    first_row_mat <- attr(precursor_mat_log2, "first_row_matrix")
    if (!is.null(first_row_mat)) {
      old_df <- scaling_table(first_row_mat)
      if (!is.null(old_df)) {
        scaling_rule_check <- tibble::tibble(
          rule = c("first row in file order (old)", "lowest Q.Value (used now)"),
          n_precursors = c(nrow(old_df), nrow(scaling_df)),
          median_fold_change = c(median(old_df$fold_change), median_fold_change),
          median_beta = c(median(old_df$beta), median_beta))
        write_csv(scaling_rule_check, file.path(DIR_F2, "EDFig3_d_duplicate_rule_check.csv"))
        message("  Effect of the duplicate-precursor rule on ED Fig 3d:"); print(as.data.frame(scaling_rule_check))
      }
    }
    
    # Intensity violin + box of the matched precursors. Only the median fold change is printed on the panel
    # (no test: a test over >1,000 precursors is always "significant" and says nothing about effect size);
    # exponent, n and the calculation go into the legend.
    n_matched <- nrow(scaling_df)
    violin_df <- dplyr::bind_rows(
      tibble::tibble(intensity = scaling_df$mean_log2_1x, group = "Single cells"),
      tibble::tibble(intensity = scaling_df$mean_log2_10x, group = "10-cell area")
    ) %>% dplyr::mutate(group = factor(group, levels = c("Single cells", "10-cell area")))   # same order as Fig 2c: single cells left, 10-cell area right
    group_n_labels <- c(`Single cells` = paste0("Single cells\n(n = ", length(scaling_1x_ids), ")"),
                        `10-cell area` = paste0("10-cell areas\n(n = ", length(scaling_10x_ids), ")"))
    y_top <- max(violin_df$intensity, na.rm = TRUE); y_rng <- diff(range(violin_df$intensity, na.rm = TRUE))
    
    edfig3d_scaling <- ggplot(violin_df, aes(group, intensity, fill = group)) +
      geom_violin(alpha = 0.35, colour = NA) +
      geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.85) +
      annotate("segment", x = 1, xend = 2, y = y_top + 0.04 * y_rng, yend = y_top + 0.04 * y_rng) +
      annotate("label", x = 1.5, y = y_top + 0.11 * y_rng, label = paste0(sprintf("%.1f", median_fold_change), "-fold"),
               fill = "#EDF0DE", colour = PAL_SHARMA, linewidth = 0, fontface = "bold", size = pt_to_geom(PT_BASE)) +
      scale_fill_manual(values = c(`Single cells` = PAL_THY1, `10-cell area` = PAL_BULK), guide = "none") +
      scale_x_discrete(labels = group_n_labels) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
      labs(x = NULL, y = "Mean log2 precursor intensity") +
      THEME_JOURNAL()
  }
}
message("  ED Fig 3 panel E: single-cell average vs. 10x area concordance by region -- precursor-level, NOT MaxLFQ...")
# Supplementary only: a 10-cell area is not a purified population.
concordance_plots_kde <- list(); concordance_summary <- list(); concordance_src <- list()
if (is.null(precursor_mat_log2)) {
  message("  SKIPPED: concordance panel -- precursor matrix not loaded, see warning above.")
  edfig3e_concordance <- NULL
} else {
  for (region_val in REGIONS_THY1) {
    sc_ids <- sample_id_lookup %>% dplyr::filter(tier=="tier_1x_single_cell", cell_type_label=="Thy1",
                                                 region==region_val, passes_qc_tier,
                                                 n_proteins_lookup[clean_sample_name] >= MIN_PROTEINS) %>% dplyr::pull(sample_id)
    bk_ids <- sample_id_lookup %>% dplyr::filter(tier=="tier_10x_area", region==region_val, passes_qc_tier) %>% dplyr::pull(sample_id)
    sc_ids <- intersect(sc_ids, colnames(precursor_mat_log2))
    bk_ids <- intersect(bk_ids, colnames(precursor_mat_log2))
    if (length(sc_ids)==0 || length(bk_ids)==0) next
    sc_avg <- rowMeans(precursor_mat_log2[, sc_ids, drop=FALSE], na.rm=TRUE)
    bulk_mean  <- rowMeans(precursor_mat_log2[, bk_ids, drop=FALSE], na.rm=TRUE)
    shared <- is.finite(sc_avg) & is.finite(bulk_mean)
    x <- sc_avg[shared]; y <- bulk_mean[shared]
    if (length(x) < 3) {
      message("  ED 3e: ", REGION_ABBREV[[region_val]], " skipped (fewer than 3 shared precursors)")
      next
    }
    if (sd(x) == 0 || sd(y) == 0) {
      message("  ED 3e: ", REGION_ABBREV[[region_val]], " skipped (constant values; correlation undefined)")
      next
    }
    r_p <- cor(x, y); r_s <- cor(x, y, method="spearman")
    df <- tibble::tibble(x=x, y=y)
    concordance_src[[region_val]] <- tibble::tibble(region = REGION_ABBREV[[region_val]], precursor = names(x), single_cell_average_log2 = as.numeric(x), area_mean_log2 = as.numeric(y))
    # Title = region and cell/area numbers; r and rho are in EDFig3_e_concordance_by_region_summary.csv.
    plot_title <- paste0(REGION_ABBREV[[region_val]], " (", length(sc_ids), " cells, ", length(bk_ids), " areas)")
    plot_x <- "Single-cell average, log2 intensity"
    plot_y <- "10-cell area, log2 intensity"
    
    # Density-coloured scatter; density is normalised per panel so all 3 regions share one colour bar.
    dens_grid <- MASS::kde2d(df$x, df$y, n = 200, lims = c(range(df$x), range(df$y)))
    interp_density <- function(x0, y0, grid) {
      ix <- pmin(pmax(findInterval(x0, grid$x), 1), length(grid$x) - 1)
      iy <- pmin(pmax(findInterval(y0, grid$y), 1), length(grid$y) - 1)
      grid$z[cbind(ix, iy)]
    }
    df$point_density <- interp_density(df$x, df$y, dens_grid)
    df$density_norm <- df$point_density / max(df$point_density)
    
    concordance_plots_kde[[region_val]] <- ggplot(df, aes(x, y, colour = density_norm)) +
      geom_point_dense(size = 0.5, alpha = 0.7, stroke = 0) +
      scale_colour_viridis_c(option = "viridis", name = "Relative density", limits = c(0, 1), breaks = c(0, 0.5, 1),
                             guide = guide_colourbar(title.position = "left", title.vjust = 0.8,
                                                     barwidth = grid::unit(28, "mm"), barheight = grid::unit(2.5, "mm"))) +
      geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey40") +
      coord_cartesian(xlim = range(df$x), ylim = range(df$y)) +
      labs(title=plot_title, x=plot_x, y=plot_y) +
      THEME_JOURNAL()
    
    concordance_summary[[region_val]] <- tibble::tibble(region=region_val, pearson_r=r_p, spearman_r=r_s, n=length(x),
                                                        n_cells=length(sc_ids), n_areas=length(bk_ids))
  }
  if (length(concordance_plots_kde) == 0) {
    message("  SKIPPED: concordance panel -- no region has both single cells and 10-cell areas in the precursor matrix.")
    edfig3e_concordance <- NULL
  } else {
    write_csv(dplyr::bind_rows(concordance_summary), file.path(DIR_F2_SUPP, "EDFig3_e_concordance_by_region_summary.csv"))
    edfig3e_concordance <- wrap_plots(concordance_plots_kde, ncol = 3, nrow = 1) +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
  }
  # ED 3e legend: per-region n, Pearson r and Spearman rho are in EDFig3_e_concordance_by_region_summary.csv.
  # Precursor-level, not MaxLFQ (see header). Point colour = local point density (KDE, normalised per panel).
}

# ---- ED Fig. 3 assembly: 180 mm (7.09 in) wide; panels that were skipped are left out ----
#   row 1: a yield by batch | b batch concordance
#   row 2: c area + yield by cell type | d precursor scaling
#   row 3: e single-cell average vs 10-cell area, one panel per region
# A missing optional input (precursor matrix, collection area) removes only the panels that need it.
ED3_W_IN <- 7.09
ED3_ROW_H_IN <- c(2.5, 2.6, 2.9)   # row 3 is taller: it also holds the shared colour bar
ed3_row <- function(panels, widths) {   # panels side by side; NULL (skipped) panels dropped; NULL if none left
  keep <- !vapply(panels, is.null, logical(1))
  if (!any(keep)) return(NULL)
  Reduce(`|`, panels[keep]) + patchwork::plot_layout(widths = widths[keep])
}
ed3_bold_x <- function(p) if (is.null(p)) NULL else p + CAT_X_BOLD
row1_ed3 <- ed3_row(list(ed3_bold_x(edfig3a_yield_batch), edfig3b_concordance), c(1.2, 1.3))
row2_ed3 <- ed3_row(list(ed3_bold_x(edfig3c_area), ed3_bold_x(edfig3c_yield), ed3_bold_x(edfig3d_scaling)), c(1.25, 1.25, 0.9))
ed3_rows <- list(row1_ed3, row2_ed3, edfig3e_concordance)
ed3_keep <- !vapply(ed3_rows, is.null, logical(1))
ed3_assembled <- Reduce(`/`, ed3_rows[ed3_keep]) + patchwork::plot_layout(heights = ED3_ROW_H_IN[ed3_keep])
source_panel("EDFig3a_yield_by_batch", edfig3a_yield_batch, "Protein groups per Thy1 single cell by preparation batch (>= 200 protein groups)",
             columns = c("sample_col", "batch", "region", "n_proteins_filtered"))
source_panel("EDFig3b_batch_concordance", edfig3b_concordance, "Batch-to-batch concordance of single-cell averages (>= 600 protein groups) per region and batch pair")
source_panel("EDFig3c_area_and_yield", if (!is.null(edfig3c_area)) edfig3c_area$data, "LMD collection area (left panel) and protein groups (right panel) per cell, Thy1 neurons and Cx3cr1 microglia (>= 200 protein groups)",
             columns = c("sample_col", "cell_type_display", "batch", "area_um2_final", "n_proteins_filtered"))
source_panel("EDFig3d_scaling", if (!is.null(edfig3d_scaling)) scaling_df, "Per precursor (detected in >= 80% of samples of both groups): mean log2 intensity in single cells and 10-cell areas, fold change and scaling exponent (beta)")
source_panel("EDFig3e_concordance", if (length(concordance_src) > 0) dplyr::bind_rows(concordance_src), "Per precursor and region: single-cell average and 10-cell area mean (log2 intensity)")
ggsave(file.path(DIR_F2_SUPP, "EDFig3_assembled.svg"), ed3_assembled,
       width = ED3_W_IN, height = sum(ED3_ROW_H_IN[ed3_keep]), units = "in", device = SVG_DEVICE)
message("  ED Fig. 3 assembled (", sum(ed3_keep), " of 3 rows): ", ED3_W_IN, " x ", sum(ED3_ROW_H_IN[ed3_keep]), " in (insert at 100%)")
write_source_data(file.path(DIR_F2, "Source_Data_Figure2.xlsx"), "Source data for Figure 2 and Extended Data Figure 3")

message("  ED Fig. 3 (a-e) saved to ", DIR_F2_SUPP)

###############################################################################
# 9. SUMMARY NUMBERS -- every number needed for the Figure 2 / ED Figure 3 legends and
# results text, collected in one place. Printed to the console AND saved to
# Figure2_summary_numbers.txt in DIR_F2.
# Each block is guarded with exists(), since a panel may have been skipped.
# All protein counts = contaminant-filtered counts (see header).
###############################################################################
summary_lines <- character()
say <- function(...) { l <- paste0(...); message(l); summary_lines <<- c(summary_lines, l) }
med_iqr <- function(x) sprintf("median %.0f (IQR %.0f-%.0f)", round_half_up(median(x, na.rm = TRUE)),
                               round_half_up(quantile(x, 0.25, na.rm = TRUE)), round_half_up(quantile(x, 0.75, na.rm = TRUE)))
fmt_p <- function(p) formatC(p, format = "e", digits = 1)

say("\n==================== FIGURE 2 & EXTENDED DATA FIGURE 3 -- SUMMARY STATISTICS ====================")

# ---- WHICH SAMPLES ARE USED WHERE (audit; stops the script on any mismatch) ----
say("\n--- WHICH SAMPLES ARE USED WHERE (counts = contaminant-filtered protein groups) ---")
membership <- cohort_meta %>%
  dplyr::filter(tier %in% c("tier_1x_single_cell", "tier_10x_area")) %>%
  dplyr::transmute(
    sample_col, tier, cell_type_label, batch,
    n_protein_groups_filtered = as.numeric(n_proteins_lookup[sample_col]),
    collected   = TRUE,
    ge200       = passes_qc_tier,
    ge600       = !is.na(n_protein_groups_filtered) & n_protein_groups_filtered >= MIN_PROTEINS,
    in_fig2cde  = sample_col %in% c(thy1_singlecell_cols, thy1_bulk10_cols))
write_csv(membership, file.path(DIR_F2, "Figure2_sample_membership.csv"))

# guards: the cells drawn in panel c are exactly the cells used for d and e
sc_c   <- sens_df$sample_col[sens_df$tier == "tier_1x_single_cell"]
area_c <- sens_df$sample_col[sens_df$tier == "tier_10x_area"]
stopifnot(setequal(sc_c,   thy1_singlecell_cols),
          setequal(area_c, thy1_bulk10_cols),
          all(n_proteins_lookup[c(sc_c, area_c)] >= MIN_PROTEINS_VIABLE))
say("CHECK passed: panel c, d and e use the SAME samples (every one >=200 filtered protein groups)")

if (exists("one_cell") && exists("plus_all")) {
  say("Panels a, b, peak capacity: HeLa digest injections only (no tissue cells). ",
      "1 cell eq.: HT n=", sum(one_cell$column == "HT"), ", Neo Plus n=", sum(one_cell$column == "PLUS"),
      "; titration (Neo Plus) n=", nrow(plus_all), " runs")
}
for (ct in c("Thy1", "Microglia")) {
  m <- membership %>% dplyr::filter(tier == "tier_1x_single_cell", cell_type_label == ct)
  say(ct, " single cells: collected ", nrow(m), " -> >=200: ", sum(m$ge200), " -> >=600: ", sum(m$ge600))
}
m10 <- membership %>% dplyr::filter(tier == "tier_10x_area")
say("10-cell areas: collected ", nrow(m10), " -> >=200: ", sum(m10$ge200), " -> >=600: ", sum(m10$ge600))
say("Panel c (sensitivity): Thy1 single cells >=200 n=", length(sc_c), " (of which >=600: ",
    sum(n_proteins_lookup[sc_c] >= MIN_PROTEINS), ") + 10-cell areas >=200 n=", length(area_c),
    " (of which >=600: ", sum(n_proteins_lookup[area_c] >= MIN_PROTEINS), "); microglia NOT shown")
say("Panel c fold labels: ", round(fold_sens, 2), "-fold = median of ALL >=200 samples (n=", length(sc_c), " cells vs ",
    length(area_c), " areas; shown on panel). ", round(fold_sens_600, 2), "-fold = medians of the >=600 samples only (",
    sum(n_proteins_lookup[sc_c] >= MIN_PROTEINS), " cells vs ", sum(n_proteins_lookup[area_c] >= MIN_PROTEINS), " areas; LEGEND TEXT ONLY)")
say("Panel d (compartments) and e (Venn): same ", length(thy1_singlecell_cols), " Thy1 single cells (>=200) and ",
    length(thy1_bulk10_cols), " 10-cell areas (>=200); a protein group counts as detected if present in >=1 of these samples")
say("Microglia are not shown in main Fig 2; they enter only the cell-flow counts, ED Fig 3c (>=200) and Figs 3-4 (>=600)")
say("ED Fig 3a: Thy1 single cells >=200 | 3b: single-cell averages from >=600 cells | 3c: Thy1 + microglia >=200 | 3d, 3e: 10-cell areas vs single cells (see their sections)")
say("Sample-by-sample list written to Figure2_sample_membership.csv")

# ---- Cell flow and protein counts per cell (Fig 2c, text, Methods) ----
say("\n--- Single cells: cell flow (all protein counts contaminant-filtered) ---")
cc <- cells_collected %>% dplyr::filter(cell_type_label %in% c("Thy1", "Microglia"))
for (ct in c("Thy1", "Microglia")) {
  d <- cc %>% dplyr::filter(cell_type_label == ct)
  d200 <- d %>% dplyr::filter(is.finite(n_proteins_filtered), n_proteins_filtered >= MIN_PROTEINS_VIABLE)
  d600 <- d %>% dplyr::filter(is.finite(n_proteins_filtered), n_proteins_filtered >= MIN_PROTEINS)
  say(ct, ": collected ", nrow(d), " | >=", MIN_PROTEINS_VIABLE, " protein groups: ", nrow(d200),
      " (", round(100 * nrow(d200) / nrow(d), 1), "%) | >=", MIN_PROTEINS, " protein groups: ", nrow(d600),
      " (", round(100 * nrow(d600) / nrow(d), 1), "%)")
  say("   protein groups per cell, >=200 set: ", med_iqr(d200$n_proteins_filtered), "; max ", max(d200$n_proteins_filtered))
  say("   protein groups per cell, >=600 set: ", med_iqr(d600$n_proteins_filtered), "; max ", max(d600$n_proteins_filtered))
}
say("All single cells collected (Thy1 + Microglia + any other): ", n_cells_collected,
    " | used for analysis (>=", MIN_PROTEINS, "): ", n_cells_used,
    " (", round(100 * n_cells_used / n_cells_collected, 1), "%)")
say("Protein groups ever measured in any single cell: ", n_proteins_ever_measured)
say("Empty/failed wells (0 protein groups after filtering): ", n_cells_empty,
    " (", round(100 * n_cells_empty / n_cells_collected, 1), "% of collected)")
say("Count bands (all collected single cells):")
for (i in seq_len(nrow(empty_well_breakdown))) say("   ", empty_well_breakdown$count_band[i], ": ", empty_well_breakdown$n_cells[i])
say("Highest single-cell count (>=200 set), per cell type:")
for (i in seq_len(nrow(single_cell_max_check))) say("   ", single_cell_max_check$cell_type_label[i], ": ",
                                                    single_cell_max_check$n_proteins_filtered[i], " (", single_cell_max_check$sample_col[i], ")")

# ---- Fig 2a: column comparison ----
if (exists("fold_increase")) {
  say("\n--- Fig 2a: column comparison (HeLa reference material, 1 cell equivalent) ---")
  for (cl in c("HT", "PLUS")) {
    r <- col_summary[col_summary$column == cl, ]
    say(cl, ": ", round(r$mean_k * 1000), " +/- ", round(r$sd_k * 1000), " protein groups (mean +/- SD), n=", r$n,
        ", ", r$spd, " SPD, nominal gradient ", GRADIENT_TIME_MIN[[cl]], " min")
  }
  say("Fold increase HT -> PLUS: ", round(fold_increase, 2), " (Welch t-test P=", fmt_p(t_res$p.value), ")")
  for (cl in c("HT", "PLUS")) {
    r <- peakcap_summary[peakcap_summary$column == cl, ]
    f <- fwhm_summary[fwhm_summary$column == cl, ]
    say("Peak capacity ", cl, ": ", round(r$mean_peakcap, 1), " +/- ", round(r$sd_peakcap, 1), "; FWHM ", round(f$mean_fwhm, 4), " +/- ", round(f$sd_fwhm, 4), " min")
  }
  say("Peak capacity PLUS/HT: ", round(peakcap_summary$mean_peakcap[peakcap_summary$column == "PLUS"] /
                                         peakcap_summary$mean_peakcap[peakcap_summary$column == "HT"], 2),
      "-fold (P=", fmt_p(peakcap_t$p.value), "); peak width (FWHM) P=", fmt_p(fwhm_t$p.value))
}

# ---- Fig 2b: titration ----
if (exists("load_summary_log")) {
  say("\n--- Fig 2b: input-amount titration (HeLa, uPAC Neo Plus) ---")
  for (i in seq_len(nrow(load_summary_log))) {
    say(as.character(load_summary_log$load_group[i]), " (", load_summary_log$cell_equiv[i], " cell eq.): ",
        round(load_summary_log$mean_k[i] * 1000), " +/- ", round(load_summary_log$sd_k[i] * 1000), " protein groups, n=", load_summary_log$n[i])
  }
}

# ---- Fig 2c: sensitivity ----
if (exists("sens_summary")) {
  say("\n--- Fig 2c: sensitivity, Thy1 single cells vs 10-cell areas ---")
  say("Cells shown = all with >=", MIN_PROTEINS_VIABLE, " protein groups; dashed line = ", MIN_PROTEINS, " (analysis threshold)")
  for (tl in levels(sens_df$tier_label)) {
    v <- sens_df$n_proteins[sens_df$tier_label == tl]
    say(tl, " (>=200): n=", length(v), ", ", med_iqr(v), "; range ", min(v), "-", max(v),
        " | of these >=600: n=", sum(v >= MIN_PROTEINS), " (", round(100 * mean(v >= MIN_PROTEINS), 1), "%), ", med_iqr(v[v >= MIN_PROTEINS]))
  }
  say("Median fold difference, 10-cell area vs single cells: ", round(fold_sens, 2), " (>=200 set, as shown in panel); ",
      round(fold_sens_600, 2), " (>=600 set only)")
  say("Single-cell batches in panel: ", paste(names(table(sens_df$batch_label[sens_df$tier_label == "Single cells"])),
                                              table(sens_df$batch_label[sens_df$tier_label == "Single cells"]), collapse = "; "))
  if ("area_um2_final" %in% colnames(cohort_meta)) {
    a1 <- cohort_meta$area_um2_final[cohort_meta$tier == "tier_1x_single_cell" & cohort_meta$cell_type_label == "Thy1" & cohort_meta$passes_qc_tier]
    a10 <- cohort_meta$area_um2_final[cohort_meta$tier == "tier_10x_area" & cohort_meta$passes_qc_tier]
    say("Collected area (um^2): Thy1 single cells ", med_iqr(a1), "; 10-cell areas ", med_iqr(a10),
        "; ratio of medians ", round(median(a10, na.rm = TRUE) / median(a1, na.rm = TRUE), 1))
  }
}

# ---- Fig 2d: compartments ----
if (exists("combined_2d")) {
  say("\n--- Fig 2d: subcellular compartment coverage, multi-label (% of ANNOTATED protein groups; can add to > 100%) ---")
  ce <- combined_2d %>% dplyr::mutate(bucket = as.character(bucket)) %>% dplyr::arrange(group, dplyr::desc(pct))
  for (g in unique(ce$group)) say(g, ": ", ce$n_annotated[ce$group == g][1], " of ", ce$n_detected[ce$group == g][1], " detected protein groups have a UniProt location annotation (",
                                  round(100 * ce$n_annotated[ce$group == g][1] / ce$n_detected[ce$group == g][1], 1), "%); the rest is left out of the percentages")
  for (i in seq_len(nrow(ce))) say(ce$group[i], " | ", ce$bucket[i], ": ", ce$n_proteins[i], " protein groups (", round(ce$pct[i], 1), "% of annotated)")
  # UniProt record (what was asked, when, which release, what came back)
  acc_used <- unique(stats::na.omit(gene_to_accession$accession[gene_to_accession$gene_label %in% all_genes_needed]))
  cache_used <- uniprot_cache_2d[uniprot_cache_2d$accession %in% acc_used, ]
  n_text_used <- sum(!is.na(cache_used$subcell_raw) & nzchar(trimws(cache_used$subcell_raw)))
  say("UniProt record: annotation uses the FIRST (lead) accession of each protein group; ", length(acc_used), " accessions of the detected protein groups; ",
      n_text_used, " (", round(100 * n_text_used / length(acc_used), 1), "%) have a subcellular-location comment, ",
      length(acc_used) - n_text_used, " have none (no comment, not returned, or request failed); ", length(setdiff(acc_used, uniprot_cache_2d$accession)), " not in the cache")
  if (file.exists(UNIPROT_LOG_FILE)) {
    ulog <- read_csv(UNIPROT_LOG_FILE, show_col_types = FALSE, col_types = cols(.default = col_character()))
    say("UniProt record: query dates ", paste(unique(ulog$date), collapse = ", "), " | UniProt release(s) ", paste(unique(stats::na.omit(ulog$uniprot_release)), collapse = ", "),
        " (release date ", paste(unique(stats::na.omit(ulog$uniprot_release_date)), collapse = ", "), ") | ", sum(as.numeric(ulog$n_recovered_by_single_query), na.rm = TRUE),
        " accessions missing from a batch answer were recovered by single queries; ", sum(as.numeric(ulog$n_not_returned), na.rm = TRUE), " not returned by UniProt at all")
  } else {
    say("UniProt record: NO query log (the cache was built before logging existed, so the UniProt release is unknown). To record it: set UNIPROT_REFRESH <- TRUE once, run the script, set it back to FALSE.")
  }
  say("Total protein groups: single cells (cumulative) ", length(detected_singlecell), " in ", length(thy1_singlecell_cols),
      " cells; 10-cell areas (cumulative) ", length(detected_10x), " in ", length(thy1_bulk10_cols), " wells")
}

# ---- Fig 2e: Venn ----
if (exists("n_overlap_2e")) {
  say("\n--- Fig 2e: protein overlap, single cells vs 10-cell areas (Venn) ---")
  say("Shared: ", n_overlap_2e, " (", pct_overlap_2e, "% of union) | 10-cell only: ", n_10x_only_2e, " (", pct_10x_only_2e,
      "%) | single-cell only: ", n_1x_only_2e, " (", pct_1x_only_2e, "%) | union: ", n_union_2e)
  say("Share of single-cell protein groups also found in 10-cell areas: ", round(100 * n_overlap_2e / n_total_1x_2e, 1), "%")
  say("Share of 10-cell protein groups also found in single cells: ", round(100 * n_overlap_2e / n_total_10x_2e, 1), "%")
  say("Set sizes: single cells ", length(set_1x), "; 10-cell areas ", length(set_10x))
}

# ---- ED Fig 3 cohort check: thresholds, batches and counts vs Fig. 2 ----
say("\n--- ED Fig 3 COHORT CHECK (thresholds, batches, counts vs Fig. 2) ---")
say("Thresholds: ED3a >=200 | ED3b >=600 cells (per batch x region) | ED3c >=200 (both cell types) | ED3d >=200 (same samples as Fig 2c) | ED3e single cells >=600 vs 10-cell areas >=200")
coh <- membership %>% dplyr::filter(tier == "tier_1x_single_cell")
for (ct in c("Thy1", "Microglia")) {
  for (b in sort(unique(coh$batch[coh$cell_type_label == ct]))) {
    m <- coh %>% dplyr::filter(cell_type_label == ct, batch == b)
    say("   ", ct, " | batch ", b, " (", as.character(BATCH_DISPLAY_LABELS[[as.character(b)]]), "): collected ", nrow(m),
        " | >=200: ", sum(m$ge200), " | >=600: ", sum(m$ge600))
  }
}
if (exists("batch_level_df")) {
  say("ED3a plotted n = ", nrow(batch_level_df), " Thy1 cells (", paste(batch_medians$n, collapse = " + "), ") vs Fig 2c: ", length(thy1_singlecell_cols),
      ifelse(nrow(batch_level_df) == length(thy1_singlecell_cols), "  -> MATCH", "  -> MISMATCH, CHECK"))
}
if (exists("n_by_type")) {
  say("ED3c plotted n: ", paste(n_by_type$cell_type_display, n_by_type$n, collapse = "; "),
      " | expected (>=200 and finite collection area): Thy1 ", length(thy1_singlecell_cols), ", microglia ", length(mg_singlecell_cols),
      ifelse(all(n_by_type$n == c(length(thy1_singlecell_cols), length(mg_singlecell_cols))), "  -> MATCH", "  -> MISMATCH (cells without a recorded area?)"))
}
if (exists("scaling_1x_ids")) {
  say("ED3d samples: ", length(scaling_1x_ids), " Thy1 single cells + ", length(scaling_10x_ids), " 10-cell areas vs Fig 2c: ",
      length(thy1_singlecell_cols), " + ", length(thy1_bulk10_cols),
      ifelse(length(scaling_1x_ids) == length(thy1_singlecell_cols) && length(scaling_10x_ids) == length(thy1_bulk10_cols), "  -> MATCH", "  -> MISMATCH (cells missing from the precursor matrix?)"))
}
if (exists("batch_corr_summary") && length(batch_corr_summary) > 0) {
  say("ED3b: batch pairs are compared only within a region; smallest group = ",
      min(c(summary_df$n_cells_1, summary_df$n_cells_2)), " cells (>=600) -- see per-comparison lines below")
}

# ---- ED Fig 3a ----
if (exists("batch_medians")) {
  say("\n--- ED Fig 3a: protein yield by batch (Thy1 single cells, >=200, descriptive) ---")
  for (i in seq_len(nrow(batch_medians))) say(BATCH_DISPLAY_LABELS[[as.character(batch_medians$batch[i])]], ": median ",
                                              batch_medians$median_proteins[i], " protein groups, n=", batch_medians$n[i])
  say("Fold difference between batch medians: ", round(batch_fold_diff, 2))
}

# ---- ED Fig 3b ----
if (exists("summary_df")) {
  say("\n--- ED Fig 3b: batch-to-batch concordance (single-cell averages, >=600 cells, per region) ---")
  say(nrow(summary_df), " region x batch-pair comparisons; Pearson r range ", round(min(summary_df$pearson_r), 2), "-",
      round(max(summary_df$pearson_r), 2), "; Spearman rho range ", round(min(summary_df$spearman_r), 2), "-", round(max(summary_df$spearman_r), 2))
  for (i in seq_len(nrow(summary_df))) say("   ", as.character(summary_df$comparison_label[i]), ": r=", round(summary_df$pearson_r[i], 2),
                                           " (95% CI ", round(summary_df$pearson_ci_lo[i], 2), "-", round(summary_df$pearson_ci_hi[i], 2), "), n=", summary_df$n[i], " protein groups; cells >=600 per batch: ",
                                           summary_df$n_cells_1[i], " (", as.character(BATCH_DISPLAY_LABELS[[as.character(summary_df$batch_1[i])]]), ") vs ",
                                           summary_df$n_cells_2[i], " (", as.character(BATCH_DISPLAY_LABELS[[as.character(summary_df$batch_2[i])]]), ")")
}

# ---- ED Fig 3c ----
if (exists("area_yield_medians")) {
  say("\n--- ED Fig 3c: collection area and yield by cell type (>=200 cells) ---")
  for (i in seq_len(nrow(area_yield_medians))) say(as.character(area_yield_medians$cell_type_display[i]), ": median area ",
                                                   round(area_yield_medians$median_area[i]), " um^2; median yield ", area_yield_medians$median_yield[i],
                                                   " protein groups, n=", n_by_type$n[n_by_type$cell_type_display == area_yield_medians$cell_type_display[i]])
  say("Descriptive only, no test: cell type is confounded with mouse and batch (microglia come from a single batch)")
}

# ---- ED Fig 3d ----
if (exists("mean_beta")) {
  say("\n--- ED Fig 3d: precursor-level scaling, 10-cell area vs Thy1 single cells, both >=200 (same cells as Fig 2c), matched precursors ---")
  say("n=", nrow(scaling_df), " precursors (>=", MIN_FRAC_DETECTED_SCALING * 100, "% detected in both groups); median fold-change ",
      round(median_fold_change, 2), ", mean ", round(mean_fold_change, 2), "; median beta ", round(median_beta, 2),
      ", mean beta ", round(mean_beta, 2), " (beta=1 would be proportional); no test (descriptive)")
  say("Area ratio actually collected (10-cell areas / single cells, ratio of medians) is printed under Fig 2c above; the exponent is relative to the NOMINAL 10-fold difference in cell number")
  if (exists("scaling_rule_check") && !is.null(scaling_rule_check)) {
    say("Duplicate-precursor rule: median fold change ", paste(sprintf("%s = %.3f", scaling_rule_check$rule, scaling_rule_check$median_fold_change), collapse = " | "))
  }
  dupc <- attr(precursor_mat_log2, "duplicate_check")
  if (!is.null(dupc)) say("Duplicate Run + Precursor.Id records: ", dupc$n_keys_duplicated, " keys duplicated, ", dupc$n_keys_differing,
                          " with different quantities (max relative difference ", signif(dupc$max_relative_difference, 3), "); rule: ", dupc$rule)
}

# ---- ED Fig 3e ----
if (exists("concordance_summary") && length(concordance_summary) > 0) {
  say("\n--- ED Fig 3e: single-cell average (cells >=600) vs 10-cell area (>=200), precursor level ---")
  for (cs in concordance_summary) say(REGION_ABBREV[[cs$region]], ": Pearson r=", round(cs$pearson_r, 2),
                                      ", Spearman rho=", round(cs$spearman_r, 2), ", n=", cs$n, " precursors; ", cs$n_cells, " single cells (>=600) vs ", cs$n_areas, " 10-cell areas (>=200)")
}

say("\n====================================== END OF SUMMARY ======================================\n")
SW_PKGS <- c("tidyverse", "dplyr", "stringr", "readr", "ggplot2", "patchwork", "ggdist", "svglite", "ggrastr", "MASS", "httr", "arrow")
say("\nSoftware versions (all packages used by this script): ", software_versions(SW_PKGS))
write_software_versions(file.path(DIR_F2, "Figure2_software_versions.txt"), SW_PKGS)
writeLines(summary_lines, file.path(DIR_F2, "Figure2_summary_numbers.txt"))
message("Summary saved to ", file.path(DIR_F2, "Figure2_summary_numbers.txt"))