###############################################################################
# 03_cell_type_identity.R -- Cell identity (Figure 3) and Extended Data (ED) Fig. 4
#
# Thy1-EGFP neurons vs Cx3cr1-EGFP microglia, single cells with >= 600 protein groups.
# Main Figure 3:
#   a  marker detection in both cell types, adjusted for proteome depth (logistic regression)
#   b  differential expression, neurons vs microglia (Wilcoxon rank-sum on log2 MaxLFQ;
#      no imputation, no normalisation)
#   c  share of literature candidates that are higher in the expected cell type
#      (detection and abundance)
#   d  neuron proteins vs published neuron datasets (UpSet)
#   e  microglia proteins vs published microglia datasets (UpSet)
# ED Fig. 4:
#   a  missingness by abundance quintile (why nothing is imputed)
#   b  tissue-free control (PEN_SC) wells: protein groups per well
#   c  marker detection in the control wells
#   d  background-protein screen of the control wells, full proteome
#   e  microglia rank-abundance vs independent datasets
#   f  volcano plots of the three DE methods (Wilcoxon, proDA, msqrob2)
#   g  agreement of the log2 fold changes between the three methods
# ED Fig. 5: heatmap of the literature candidates, cell by cell.
# Panel letters are added afterwards in a vector editor.
#
# Sections:
#   1. Setup
#   2. Shared data preparation
#   3. Cell-type tests (used by panels a-c): Wilcoxon DE, detection test
#   4. Literature candidates
#   5. Fig. 3a: marker detection
#   6. Fig. 3b: volcano
#   7. Fig. 3c: candidate summary
#   8. Fig. 3d, 3e: UpSet plots
#   9. Main Figure 3: assembly
#  10. ED Fig. 4a: missingness
#  11. ED Fig. 4b-d: tissue-free control wells
#  12. ED Fig. 4e: rank-abundance
#  13. ED Fig. 4f, 4g: DE method comparison
#  14. ED Fig. 4: assembly
#  15. ED Fig. 5: heatmap of the literature candidates
#  16. Summary numbers
#  17. Source data (Excel)
#
# Cell inclusion (contaminant-filtered protein groups per cell; MIN_PROTEINS in setup.R):
# >= 600 (analysis set) for every cell-level analysis here; the rank-abundance panel (ED 4e) uses
# the >= 200 cohorts defined in setup.R. Never the raw count before contaminant removal.
#
# Depth: the number of proteins detected differs between cells and between the two cell types
# (microglia have deeper profiles). The detection tests therefore include log2(filtered protein
# count) as a covariate (panel a, panel c, the detection test). Batch is NOT in these models:
# microglia come from a single batch, so batch and cell type cannot be separated. The abundance test
# of panel b (Wilcoxon) does not use depth; ED 4 checks it against lm(log2 intensity ~ cell type +
# log2 depth) as a sensitivity analysis (section 13).
#
# Inputs : the master dataset and cohorts from setup.R; the literature tables listed in setup.R
#          (Sharma 2015, Rayaprolu 2022, Johnson 2026, Lloyd 2024, Neuron 2026 Table S3);
#          Results/master_dataset/pen_sc_systematic_screen_ALL_PROTEINS.csv (ED 4d);
#          the DIA-NN precursor matrix .tsv (msqrob2, ED 4f-g; file name set in section 1).
# Outputs: Figure3_MAIN_assembled.svg, EDFig4_assembled.svg, SuppFig_heatmap_literature_candidates.svg,
#          the CSV tables behind every panel, Figure3_cell_membership.csv, Figure3_summary_numbers.txt.
###############################################################################
# Run from the project folder, or set options(scp.project_root = "/path/to/project") first.
source(file.path(getOption("scp.project_root", getwd()), "Scripts", "shared", "setup.R"))

###############################################################################
# 1. SETUP
###############################################################################
check_packages(c("ComplexUpset", "data.table", "proDA", "QFeatures", "msqrob2", "MsCoreUtils",
                 "SummarizedExperiment", "S4Vectors"), "03_cell_type_identity.R")

# Literature tables (paths in setup.R) and the precursor matrix used by msqrob2
for (lit_file in c(SHARMA_FILE, LLOYD_FILE, JOHNSON_FILE, RAYAPROLU2022_FILE, NEURON2026_FILE)) {
  require_file(lit_file, "Literature table", "Put the file in the project's Data/ folder (names are set in setup.R).")
}
PR_MATRIX_FILE <- file.path(PROJECT_ROOT, "Data", "report_CombinedBrainEmpLib_05Jun2026_360_MBR_pr_matrix.tsv")
require_file(PR_MATRIX_FILE, "DIA-NN precursor matrix (pr_matrix.tsv, needed for msqrob2 in ED Fig. 4f-g)")

# Journal geometry. Main Figure 3: row 1 = a detection | b volcano | c candidate summary (each plotting
# area exactly FIG3_PANEL_H_IN high); row 2 = d neuron UpSet | e microglia UpSet. ED Fig. 4: four rows.
FIG3_W_IN       <- mm_to_in(FIG_W_MM)
FIG3_PANEL_H_IN <- 2.0
FIG3_ROW_H_IN   <- c(3.0, 3.2)
ED4_ROW_H_IN    <- rep(2.4, 4)              # 9.6 in in total (journal maximum is about 9.7 in)

# Proteins shown in panels a and b (and marked in the heatmap). Chosen from the literature candidate
# table written in section 4 (Figure3_literature_candidates_tests.csv).
FIG3_PROTEINS_NEURON    <- c("Cplx2", "Amph", "Camk2a", "Slc12a5", "Syp", "Necap1", "Cend1", "Ppfia3", "Uchl1")
FIG3_PROTEINS_MICROGLIA <- c("Sh3bp1", "Arhgdib", "Pygl", "Gbp5", "Ctsb")
# Rank-abundance panel (ED 4e): Aif1 + the proteins above + two further immune proteins + the
# canonical microglia markers; markers not detected in a dataset are drawn as open circles.
MICROGLIA_CANONICAL  <- c("P2ry12", "Tmem119", "Csf1r", "Cx3cr1", "Trem2", "Sall1", "Siglech", "Fcrls", "Itgam")
MICROGLIA_RANK_PANEL <- unique(c("Aif1", FIG3_PROTEINS_MICROGLIA, "Itgb2", "Cybb", MICROGLIA_CANONICAL))

# Literature rules (section 4) and UpSet rule (section 8)
LIT_MIN_CELLS   <- 3      # a candidate must be detected in >= 3 of our cells (both cell types pooled) to be testable
LIT_FDR         <- 0.05
LIT_LOG2FOLD    <- 1      # Sharma 2015: at least 2-fold higher in the target cell type
MIN_CELLS_UPSET <- 3      # a protein counts as detected in this study if seen in >= 3 cells

###############################################################################
# 2. SHARED DATA PREPARATION
# The analysis set: Thy1 and microglia single cells with >= MIN_PROTEINS protein groups.
###############################################################################
message("\n[Figure 3] Preparing the analysis set...")

de_meta <- cohort_meta %>%
  dplyr::filter(tier == "tier_1x_single_cell", cell_type_label %in% c("Thy1", "Microglia")) %>%
  dplyr::select(-dplyr::any_of("sample_id")) %>%
  dplyr::rename(cell_type = cell_type_label, sample_id = sample_col)
mat_de <- protein_mat_linear[, de_meta$sample_id]
pass_cells <- colSums(!is.na(mat_de)) >= MIN_PROTEINS
mat_de <- mat_de[, pass_cells]
de_cells <- de_meta[de_meta$sample_id %in% colnames(mat_de), ] %>%
  dplyr::mutate(batch = as.factor(batch))
mat_de <- mat_de[rowMeans(!is.na(mat_de)) >= MISS_THRESHOLD, ]
# No normalisation: plain log2 MaxLFQ
mat_de_log <- log2(mat_de)
thy1_ids <- de_cells %>% dplyr::filter(cell_type == "Thy1") %>% dplyr::pull(sample_id)
mg_ids   <- de_cells %>% dplyr::filter(cell_type == "Microglia") %>% dplyr::pull(sample_id)
message("  Analysis set: ", length(thy1_ids), " Thy1 neurons + ", length(mg_ids), " microglia (>= ", MIN_PROTEINS,
        " protein groups); log2 matrix (not normalised): ", nrow(mat_de_log), " proteins detected in >= ",
        100 * MISS_THRESHOLD, "% of these cells")

# Per-cell table for the detection tests: log2 depth, centred on the pooled median, so that predictions
# at depth_c = 0 are "at the pooled median depth". Microglia = reference level.
cell_df <- tibble::tibble(
  sample_id  = c(mg_ids, thy1_ids),
  cell_type  = factor(c(rep("Microglia", length(mg_ids)), rep("Thy1", length(thy1_ids))), levels = c("Microglia", "Thy1")),
  depth_log2 = log2(as.numeric(n_proteins_lookup[c(mg_ids, thy1_ids)]))
)
if (anyNA(cell_df$depth_log2)) stop("Missing protein counts (n_proteins_lookup) for cells of the analysis set.", call. = FALSE)
ref_depth_log2 <- median(cell_df$depth_log2)
cell_df$depth_c <- cell_df$depth_log2 - ref_depth_log2
message("  Reference depth for the detection tests (pooled median): ", round(2^ref_depth_log2), " protein groups")

# ---- Depth overlap check for Fig. 3a ----
# Fig. 3a predicts detection at the pooled median depth. This is only reliable if cells of both types
# have depths around that value: where does the reference depth lie within each cell type, and how much do the
# two cell types overlap in depth?
ref_depth_proteins <- 2^ref_depth_log2
depth_overlap <- tibble::tibble(cell_type = as.character(cell_df$cell_type), proteins = 2^cell_df$depth_log2) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(cells = dplyr::n(), depth_min = min(proteins), depth_q25 = quantile(proteins, 0.25), depth_median = median(proteins),
                   depth_q75 = quantile(proteins, 0.75), depth_max = max(proteins),
                   ref_depth_percentile = round(100 * mean(proteins <= ref_depth_proteins)),
                   ref_inside_range = ref_depth_proteins >= min(proteins) & ref_depth_proteins <= max(proteins),
                   ref_inside_IQR = ref_depth_proteins >= quantile(proteins, 0.25) & ref_depth_proteins <= quantile(proteins, 0.75), .groups = "drop")
overlap_lo <- max(depth_overlap$depth_min); overlap_hi <- min(depth_overlap$depth_max)
depth_overlap$share_cells_in_overlap_range <- vapply(depth_overlap$cell_type, function(ct) {
  p <- 2^cell_df$depth_log2[as.character(cell_df$cell_type) == ct]; round(100 * mean(p >= overlap_lo & p <= overlap_hi)) }, numeric(1))
write_csv(depth_overlap, file.path(DIR_F3_SUPP, "Figure3a_depth_overlap_check.csv"))
message("  Depth overlap check: reference depth ", round(ref_depth_proteins), " protein groups; ",
        paste(sprintf("%s: percentile %d, inside range %s, inside IQR %s", depth_overlap$cell_type, depth_overlap$ref_depth_percentile,
                      depth_overlap$ref_inside_range, depth_overlap$ref_inside_IQR), collapse = " | "),
        "; depth range shared by both cell types: ", round(overlap_lo), "-", round(overlap_hi))
p_depth_overlap <- ggplot(cell_df, aes(x = 2^depth_log2, fill = cell_type)) +
  geom_density(alpha = 0.5, colour = NA, adjust = 0.8) +
  geom_vline(xintercept = ref_depth_proteins, linetype = "dashed", colour = "grey20", linewidth = 0.4) +
  annotate("text", x = ref_depth_proteins, y = Inf, label = "reference depth", hjust = -0.08, vjust = 1.6, size = pt_to_geom(PT_BASE)) +
  scale_fill_manual(values = c(Thy1 = PAL_THY1, Microglia = PAL_MICROGLIA), labels = c(Thy1 = LAB_NEURON, Microglia = LAB_MICROGLIA), name = NULL) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Protein groups per cell (log scale)", y = "Density") +
  THEME_JOURNAL() + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
ggsave(file.path(DIR_F3_SUPP, "Figure3a_depth_overlap_check.svg"), p_depth_overlap, width = 3.4, height = 2.6, device = SVG_DEVICE)

# Detection matrix (contaminant-filtered proteins x analysis-set cells) and detection rate per cell type
det_mat     <- !is.na(protein_mat_log2[, cell_df$sample_id, drop = FALSE])
det_thy1    <- rowMeans(det_mat[, cell_df$cell_type == "Thy1",      drop = FALSE])
det_mg      <- rowMeans(det_mat[, cell_df$cell_type == "Microglia", drop = FALSE])
label_by_upper <- rownames(protein_mat_log2); names(label_by_upper) <- toupper(label_by_upper)
label_by_upper <- label_by_upper[!duplicated(names(label_by_upper))]      # UPPER-CASE name -> label used here

# Logistic regression of detection (yes/no) on cell type + log2 depth for one protein. Returns the fitted
# model, or NULL if the protein is detected in all or in none of the cells, or the fit fails.
fit_detection <- function(g) {
  if (!g %in% rownames(det_mat)) return(NULL)
  d <- cell_df
  d$detected <- as.integer(det_mat[g, ])
  if (length(unique(d$detected)) < 2) return(NULL)
  tryCatch(suppressWarnings(glm(detected ~ cell_type + depth_c, data = d, family = binomial)), error = function(e) NULL)
}
inv_logit <- function(x) 1 / (1 + exp(-x))

# Proteins chosen for the panels: check that they exist in the contaminant-filtered matrix
resolve_proteins <- function(x, what) {
  lab <- unname(label_by_upper[toupper(x)])
  for (i in which(is.na(lab))) message("  skipped ", x[i], " (", what, "): not in the contaminant-filtered protein matrix")
  unique(lab[!is.na(lab)])
}
neuron_markers   <- resolve_proteins(FIG3_PROTEINS_NEURON, "neuron")
microglia_markers <- resolve_proteins(FIG3_PROTEINS_MICROGLIA, "microglia")
all_markers <- c(neuron_markers, microglia_markers)
if (length(neuron_markers) == 0 || length(microglia_markers) == 0) stop("No neuron or no microglia proteins left to show; check FIG3_PROTEINS_* in section 1.", call. = FALSE)
message("  Neuron proteins   : ", paste(neuron_markers, collapse = ", "))
message("  Microglia proteins: ", paste(microglia_markers, collapse = ", "))

###############################################################################
# 3. CELL-TYPE TESTS (used by panels a-c)
###############################################################################

# ---- 3.1 Differential expression: Wilcoxon rank-sum on the non-imputed log2 data ----------------------
# Primary DE method of the manuscript. Per protein: cells where it was detected, at least 3 per cell type;
# log2FC = difference of medians (neurons - microglia); BH over all tested proteins.
message("\n[Figure 3b] Cell-type DE (Wilcoxon, non-imputed log2 data)...")
de_ct <- purrr::map_dfr(rownames(mat_de_log), function(g) {
  x <- mat_de_log[g, thy1_ids]; x <- x[!is.na(x)]
  y <- mat_de_log[g, mg_ids];   y <- y[!is.na(y)]
  if (length(x) < 3 || length(y) < 3) return(NULL)
  wt <- wilcox.test(x, y)
  tibble::tibble(gene = g, logFC = median(x) - median(y), P.Value = wt$p.value, n_a = length(x), n_b = length(y))
})
de_ct$adj.P.Val <- p.adjust(de_ct$P.Value, method = "BH")
de_ct <- de_ct %>% dplyr::arrange(P.Value)

# Contamination flags for every tested protein: proteins matching ribosomal/proteasome or blood/immune name
# patterns are flagged, not removed (the contaminant filter itself is applied in setup.R). The flags come from
# contamination_annotation (setup.R).
de_ct <- de_ct %>%
  dplyr::left_join(contamination_annotation %>%
                     dplyr::select(gene_label, matches_ribo_proteasome_pattern, matches_blood_immune_pattern),
                   by = c("gene" = "gene_label")) %>%
  dplyr::mutate(contamination_flagged = matches_ribo_proteasome_pattern | matches_blood_immune_pattern)
n_flagged_ribo  <- sum(de_ct$matches_ribo_proteasome_pattern, na.rm = TRUE)
n_flagged_blood <- sum(de_ct$matches_blood_immune_pattern, na.rm = TRUE)
n_ct_flagged_sig <- sum(de_ct$adj.P.Val < 0.05 & de_ct$contamination_flagged, na.rm = TRUE)
message("  ", sum(de_ct$adj.P.Val < 0.05, na.rm = TRUE), " / ", nrow(de_ct), " proteins at FDR q<0.05; ",
        n_ct_flagged_sig, " of those match ribosomal/proteasome or blood/immune name patterns (flagged, not removed)")
write_csv(de_ct, file.path(DIR_F3, "Figure3b_DE_wilcoxon_all_proteins.csv"))

# ---- 3.2 Detection test: detected yes/no ~ cell type + log2 depth, for every protein -----------------
# BH over all proteins detected in >= 10% of the cells of either cell type. z_thy1 > 0: detected more often in neurons.
message("  Detection test (logistic, detected ~ cell type + log2 depth)...")
detect_genes <- rownames(protein_mat_log2)[det_thy1 >= 0.10 | det_mg >= 0.10]
de_detect <- purrr::map_dfr(detect_genes, function(g) {
  fit <- fit_detection(g)
  if (is.null(fit) || !isTRUE(fit$converged)) return(NULL)
  cf <- suppressWarnings(summary(fit)$coefficients)
  if (!"cell_typeThy1" %in% rownames(cf)) return(NULL)
  est <- cf["cell_typeThy1", "Estimate"]; se <- cf["cell_typeThy1", "Std. Error"]
  if (!is.finite(se) || se >= 10) return(NULL)                                # unstable fit (very rarely detected protein)
  tibble::tibble(gene = g, est_thy1_vs_mg = est, se = se, z_thy1 = est / se,
                 P.Value = cf["cell_typeThy1", "Pr(>|z|)"],
                 detected_thy1_pct = 100 * det_thy1[[g]], detected_mg_pct = 100 * det_mg[[g]])
})
de_detect$adj.P.Val <- p.adjust(de_detect$P.Value, method = "BH")
write_csv(de_detect, file.path(DIR_F3, "Figure3_detection_test_logistic.csv"))
message("  Detection test: ", nrow(de_detect), " proteins, ", sum(de_detect$adj.P.Val < LIT_FDR), " different at FDR<", LIT_FDR,
        "; intensity test (Wilcoxon): ", nrow(de_ct), " proteins, ", sum(de_ct$adj.P.Val < LIT_FDR, na.rm = TRUE), " different")

###############################################################################
# 4. LITERATURE CANDIDATES
# Objective, rule-based list of proteins that published datasets assign to neurons or to microglia.
# Nothing here selects on differences in our data. Rules:
#   - Sharma 2015: quantified in BOTH isolated neurons and isolated microglia, at least 2-fold higher in
#     the target cell type
#   - independent datasets: neurons   -> Rayaprolu 2022 (Camk2a-TurboID, >= 2-fold, FDR < 0.05) OR
#                                        Johnson 2026 (single neurons), on top of Sharma
#                           microglia -> Lloyd 2024 (ex vivo, detected in >= 3 of 4 samples) AND
#                                        Neuron 2026 single-cell Table S3
#   - this study: only testability. Detected in at least LIT_MIN_CELLS cells of the analysis set (both cell types
#     pooled, so the filter does not look at the expected cell type). No selection on the size or direction of any difference.
# Every candidate is then tested in our data (detection and abundance, section 3) without filtering.
# Written: Figure3_literature_lists.csv (candidates + detection rates),
#          Figure3_literature_candidates_tests.csv (candidates + test results; the table to choose the
#          proteins of panels a and b from).
###############################################################################
message("\n[Figure 3] Literature candidates...")

# Sharma 2015: the intensity column of each isolated cell type is the first column block ("...4")
sharma_raw <- suppressMessages(readxl::read_excel(SHARMA_FILE, sheet = SHARMA_SHEET, skip = SHARMA_SKIP))   # duplicate-header renaming messages are expected
sharma_column <- function(tbl, keyword) {
  exact <- paste0("Isolated ", keyword, "...4")
  if (exact %in% colnames(tbl)) return(exact)
  cols <- colnames(tbl)[stringr::str_detect(colnames(tbl), paste0("^Isolated ", keyword))]
  if (length(cols) == 0) stop("Sharma table: no column starting with 'Isolated ", keyword, "'", call. = FALSE)
  suffix <- suppressWarnings(as.integer(stringr::str_extract(cols, "\\d+$")))
  if (all(is.na(suffix))) cols[1] else cols[which.min(suffix)]
}
first_gene <- function(x) toupper(stringr::str_extract(as.character(x), "^[^;]+"))
sharma_mg_col <- sharma_column(sharma_raw, "Microglia")
sharma_ne_col <- sharma_column(sharma_raw, "Neurons")
message("  Sharma columns used: '", sharma_ne_col, "' (neurons) and '", sharma_mg_col, "' (microglia)")
sharma_fold <- tibble::tibble(
  G   = first_gene(sharma_raw[["Gene names"]]),
  mgv = suppressWarnings(as.numeric(sharma_raw[[sharma_mg_col]])),
  nev = suppressWarnings(as.numeric(sharma_raw[[sharma_ne_col]]))
) %>%
  dplyr::filter(!is.na(G), G != "", !is.na(mgv), !is.na(nev)) %>%           # quantified in BOTH cell types
  dplyr::arrange(dplyr::desc(mgv)) %>% dplyr::distinct(G, .keep_all = TRUE) %>%
  dplyr::mutate(neuron = nev - mgv, microglia = mgv - nev)                  # difference toward the target cell type (log2 scale)

# Reference tables, read once and used for the candidate rules here and for the UpSet plots (section 8)
lloyd_raw  <- readxl::read_excel(LLOYD_FILE, sheet = "Mouse")
lloyd_n3   <- rowSums(sapply(lloyd_raw[grep("Copy number exMM", colnames(lloyd_raw), value = TRUE)],
                             function(x) suppressWarnings(as.numeric(x))) > 0, na.rm = TRUE)
set_lloyd  <- unique(stats::na.omit(first_gene(lloyd_raw[["Gene names"]][lloyd_n3 >= 3])))      # detected in >= 3 of 4 ex vivo samples
s3_raw     <- readxl::read_excel(NEURON2026_FILE, skip = 1)                                      # row 1 is the title
set_s3     <- unique(stats::na.omit(toupper(as.character(s3_raw[["Gene_name"]]))))
ray_raw    <- readxl::read_excel(RAYAPROLU2022_FILE, sheet = RAYAPROLU2022_SHEET, skip = RAYAPROLU2022_SKIP)
stopifnot(all(c("Gene Symbols", "Majority Protein IDs", "Log2 (difference)", "BH-FDR") %in% names(ray_raw)))
# Rayaprolu 2022 lists everything quantified incl. background: keep neuron-labelled = enriched >= 2-fold, BH-FDR < 0.05
ray_enriched <- !is.na(ray_raw[["Gene Symbols"]]) & as.numeric(ray_raw[["Log2 (difference)"]]) >= 1 & as.numeric(ray_raw[["BH-FDR"]]) < 0.05
set_ray    <- unique(stats::na.omit(toupper(trimws(as.character(ray_raw[["Gene Symbols"]][ray_enriched])))))
# Johnson 2026 (Anal Chem), Table S3: DA | PVN | 5-HT blocks side by side (UniProt acc. in cols 1/5/9, gene symbol in 3/7/11)
johnson_raw <- readxl::read_excel(JOHNSON_FILE, sheet = "Table S3", skip = 3, col_names = FALSE)
set_johnson <- unique(stats::na.omit(toupper(sub("[; ].*$", "", trimws(unlist(lapply(c(3, 7, 11), function(j) as.character(johnson_raw[[j]][-1]))))))))
set_johnson <- set_johnson[set_johnson != "" & set_johnson != "NA"]
message("  Rayaprolu 2022: ", sum(!is.na(ray_raw[["Gene Symbols"]])), " quantified rows, ", sum(ray_enriched),
        " enriched (log2FC>=1, FDR<0.05) | Johnson 2026 neurons: ", length(set_johnson), " proteins")

lit_lists <- purrr::map_dfr(c("neuron", "microglia"), function(cell) {
  # cells (both cell types pooled) in which each protein was detected; label-blind, so candidates are not
  # selected for being seen in their expected cell type
  det_target <- det_thy1 * sum(cell_df$cell_type == "Thy1") + det_mg * sum(cell_df$cell_type == "Microglia")
  names(det_target) <- toupper(names(det_target))
  s0 <- sharma_fold[sharma_fold[[cell]] >= LIT_LOG2FOLD, ]
  in_refs <- if (cell == "neuron") ((s0$G %in% set_ray) | (s0$G %in% set_johnson)) else (s0$G %in% set_lloyd & s0$G %in% set_s3)
  s1 <- s0[in_refs, ]
  d  <- det_target[s1$G]
  s2 <- s1[!is.na(d) & round(d) >= LIT_MIN_CELLS, ]
  message(sprintf("  %-9s: Sharma >=2-fold %5d -> also in the independent dataset(s) %4d -> detected in >= %d of our cells (pooled) %4d",
                  cell, nrow(s0), nrow(s1), LIT_MIN_CELLS, nrow(s2)))
  s2 <- s2[order(s2[[cell]], decreasing = TRUE), ]
  tibble::tibble(
    cell_type = cell, gene = unname(label_by_upper[s2$G]), rank_by_literature_fold = seq_len(nrow(s2)),
    literature_fold_Sharma = round(2^s2[[cell]], 1),
    detected_Thy1_pct = round(100 * det_thy1[unname(label_by_upper[s2$G])]),
    detected_microglia_pct = round(100 * det_mg[unname(label_by_upper[s2$G])])
  )
})
write_csv(lit_lists, file.path(DIR_F3, "Figure3_literature_lists.csv"))

# Every candidate with its result in our data (no selection on our data). clear_in_detection_test marks the
# candidates that are also clearly higher (FDR < 0.05) in the expected cell type (used by the heatmap, section 15).
lit_tests <- purrr::map_dfr(c("neuron", "microglia"), function(cell) {
  s <- if (cell == "neuron") 1 else -1
  lit_lists[lit_lists$cell_type == cell, ] %>%
    dplyr::left_join(de_detect %>% dplyr::select(gene, detection_test_q = adj.P.Val, z_thy1), by = "gene") %>%
    dplyr::left_join(de_ct %>% dplyr::select(gene, intensity_diff_Thy1_minus_microglia = logFC, intensity_test_q = adj.P.Val), by = "gene") %>%
    dplyr::transmute(cell_type = cell, gene, literature_fold_Sharma, detected_Thy1_pct, detected_microglia_pct,
                     detection_test_q, intensity_diff_Thy1_minus_microglia, intensity_test_q,
                     clear_in_detection_test = !is.na(detection_test_q) & detection_test_q < LIT_FDR & s * z_thy1 > 0,
                     also_clear_in_intensity_test = !is.na(intensity_test_q) & intensity_test_q < LIT_FDR & s * intensity_diff_Thy1_minus_microglia > 0)
}) %>% dplyr::arrange(cell_type, detection_test_q)
write_csv(lit_tests, file.path(DIR_F3, "Figure3_literature_candidates_tests.csv"))
message("  Literature candidates (no selection on our data): neurons ", sum(lit_tests$cell_type == "neuron"),
        ", microglia ", sum(lit_tests$cell_type == "microglia"), " -> ", file.path(DIR_F3, "Figure3_literature_candidates_tests.csv"))

# Are the proteins chosen for panels a and b literature candidates? (information only)
for (i in seq_along(all_markers)) {
  ct <- if (all_markers[i] %in% neuron_markers) "neuron" else "microglia"
  if (!all_markers[i] %in% lit_tests$gene[lit_tests$cell_type == ct])
    message("  NOTE: ", all_markers[i], " (", ct, ") is not in the literature candidate table")
}

###############################################################################
# 5. FIG. 3a: MARKER DETECTION, BOTH CELL TYPES (dumbbell)
# Raw detection percentages mix real cell-type differences with each cell type's typical proteome depth
# (microglia profiles are deeper). Per protein: logistic regression (detected ~ cell type + log2 depth),
# then the predicted detection probability of each cell type AT THE POOLED MEDIAN DEPTH. The 95% CI is
# computed on the link scale and back-transformed, so it stays inside 0-100%. Read the interval, not
# just the point: several proteins with different point estimates have overlapping intervals.
###############################################################################
message("\n[Figure 3a] Marker detection, depth-adjusted...")
message("  Cells: ", LAB_NEURON, " n=", length(thy1_ids), " | ", LAB_MICROGLIA, " n=", length(mg_ids))

newdata_ref <- tibble::tibble(cell_type = factor(c("Thy1", "Microglia"), levels = c("Microglia", "Thy1")), depth_c = 0)
detect_wide <- purrr::map_dfr(all_markers, function(g) {
  fit <- fit_detection(g)
  if (is.null(fit)) return(NULL)
  pl <- predict(fit, newdata = newdata_ref, type = "link", se.fit = TRUE)
  tibble::tibble(
    gene = g,
    marker_set = dplyr::if_else(g %in% neuron_markers, "Neuronal markers", "Microglial markers"),
    pct_thy1 = 100 * inv_logit(pl$fit[1]),
    pct_thy1_lo = 100 * inv_logit(pl$fit[1] - 1.96 * pl$se.fit[1]),
    pct_thy1_hi = 100 * inv_logit(pl$fit[1] + 1.96 * pl$se.fit[1]),
    pct_microglia = 100 * inv_logit(pl$fit[2]),
    pct_microglia_lo = 100 * inv_logit(pl$fit[2] - 1.96 * pl$se.fit[2]),
    pct_microglia_hi = 100 * inv_logit(pl$fit[2] + 1.96 * pl$se.fit[2])
  )
})
n_unmodeled <- length(all_markers) - nrow(detect_wide)
if (n_unmodeled > 0) {
  message("  NOTE: ", n_unmodeled, " protein(s) could not be modelled (detected in all or none of the cells, or the fit failed): ",
          paste(setdiff(all_markers, detect_wide$gene), collapse = ", "))
}
write_csv(detect_wide, file.path(DIR_F3, "Figure3a_detection_depth_adjusted.csv"))

detect_both <- dplyr::bind_rows(
  detect_wide %>% dplyr::transmute(gene, marker_set, cell_type = LAB_NEURON,
                                   pct_detected = pct_thy1, pct_lo = pct_thy1_lo, pct_hi = pct_thy1_hi),
  detect_wide %>% dplyr::transmute(gene, marker_set, cell_type = LAB_MICROGLIA,
                                   pct_detected = pct_microglia, pct_lo = pct_microglia_lo, pct_hi = pct_microglia_hi)
)
marker_order <- detect_both %>% dplyr::filter(cell_type == LAB_NEURON) %>%
  dplyr::arrange(marker_set, dplyr::desc(pct_detected)) %>% dplyr::pull(gene)
detect_wide <- detect_wide %>% dplyr::mutate(gene = factor(gene, levels = rev(marker_order)))
detect_both <- detect_both %>% dplyr::mutate(gene = factor(gene, levels = rev(marker_order)))

p3a_detection <- ggplot() +
  geom_segment(data = detect_wide, aes(x = pct_thy1, xend = pct_microglia, y = gene, yend = gene),
               colour = "grey70", linewidth = 0.6) +
  geom_errorbarh(data = detect_both, aes(xmin = pct_lo, xmax = pct_hi, y = gene, colour = cell_type),
                 height = 0.25, alpha = 0.5, show.legend = FALSE) +
  geom_point(data = detect_both, aes(pct_detected, gene, colour = cell_type), size = 1.8, alpha = 0.9) +
  scale_colour_manual(values = setNames(c(PAL_THY1, PAL_MICROGLIA), c(LAB_NEURON, LAB_MICROGLIA)), name = NULL) +
  facet_grid(marker_set ~ ., scales = "free_y", space = "free_y",
             labeller = ggplot2::as_labeller(c("Neuronal markers" = "Neuronal\nmarkers", "Microglial markers" = "Microglial\nmarkers"))) +   # two lines so the label fits
  labs(x = "Predicted detection (%), depth-adjusted, 95% CI", y = NULL) +
  THEME_JOURNAL() + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), panel.spacing.y = grid::unit(1.5, "mm"))

###############################################################################
# 6. FIG. 3b: VOLCANO, NEURONS vs MICROGLIA (Wilcoxon; the test itself is in section 3.1)
# Highlighted: the proteins chosen in section 1. Grey: all other tested proteins.
###############################################################################
message("\n[Figure 3b] Volcano...")
build_volcano <- function(de, neuron_set, microglia_set) {
  de_marked <- de %>%
    dplyr::mutate(
      marker_set = dplyr::case_when(gene %in% neuron_set ~ "Neuronal markers",
                                    gene %in% microglia_set ~ "Microglial markers",
                                    TRUE ~ NA_character_),
      label = dplyr::if_else(!is.na(marker_set), gene, NA_character_))
  ggplot(de_marked, aes(logFC, -log10(adj.P.Val))) +
    geom_point_dense(data = dplyr::filter(de_marked, is.na(marker_set)), colour = "grey80", size = 0.6, alpha = 0.4) +
    geom_point(data = dplyr::filter(de_marked, !is.na(marker_set)), aes(colour = marker_set), size = 1.5) +
    ggrepel::geom_text_repel(data = dplyr::filter(de_marked, !is.na(marker_set)), aes(label = label, colour = marker_set),
                             size = pt_to_geom(PT_MIN), max.overlaps = 30, fontface = "bold", show.legend = FALSE,
                             segment.size = 0.25, min.segment.length = 0.1) +
    scale_colour_manual(values = c("Neuronal markers" = PAL_THY1, "Microglial markers" = PAL_MICROGLIA), name = NULL) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +   # points/axis start exactly at 0
    scale_x_continuous(limits = c(-1, 1) * max(abs(de_marked$logFC), na.rm = TRUE)) +   # symmetric limits: log2FC = 0 in the middle
    labs(x = "log2FC (neurons \u2212 microglia)", y = "\u2212log10(adj. P)") +
    THEME_JOURNAL() + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
}
p3b_volcano <- build_volcano(de_ct, neuron_markers, microglia_markers)
not_testable <- setdiff(all_markers, de_ct$gene)
if (length(not_testable) > 0) message("  Not testable in panel b (fewer than 3 detections in one cell type): ", paste(not_testable, collapse = ", "))
message("  ", sum(de_ct$adj.P.Val < 0.05, na.rm = TRUE), " / ", nrow(de_ct), " proteins significant at FDR q<0.05")

###############################################################################
# 7. FIG. 3c: DO THE LITERATURE CANDIDATES BEHAVE AS EXPECTED IN OUR DATA?
# For every candidate (section 4), in OUR data:
#   detection: is it detected more often in the expected cell type, at the same proteome depth?
#              (logistic regression as in panel a, predicted at the pooled median depth)
#   abundance: is it more abundant in the expected cell type? (Wilcoxon, section 3.1; log2FC toward the
#              expected cell type)
# Proteins detected in >= 50% of the cells of BOTH types are reported separately: detection frequency
# cannot separate cell types for them (ceiling effect).
# Panel c shows, per cell type, the share of candidates that go the expected way, for the candidates
# that could be tested in BOTH readouts (same n for both bars).
###############################################################################
message("\n[Figure 3c] Check of all literature candidates...")
candidate_detection <- function(g) {
  fit <- fit_detection(g)
  if (is.null(fit)) return(NULL)
  pr <- inv_logit(predict(fit, newdata = newdata_ref, type = "link"))
  cf <- suppressWarnings(summary(fit)$coefficients)
  tibble::tibble(gene = g, pct_thy1 = 100 * pr[1], pct_microglia = 100 * pr[2],
                 est_thy1_vs_mg = unname(coef(fit)["cell_typeThy1"]),
                 det_P = if ("cell_typeThy1" %in% rownames(cf)) cf["cell_typeThy1", "Pr(>|z|)"] else NA_real_)
}
cand_det <- purrr::map_dfr(unique(lit_lists$gene), candidate_detection)

cand_all <- lit_lists %>%
  dplyr::select(cell_type, gene, literature_fold_Sharma) %>%
  dplyr::left_join(cand_det, by = "gene") %>%
  dplyr::left_join(de_ct %>% dplyr::select(gene, logFC, abu_P = P.Value), by = "gene") %>%
  dplyr::mutate(
    det_thy1_pct = 100 * det_thy1[gene], det_mg_pct = 100 * det_mg[gene],
    pct_expected = dplyr::if_else(cell_type == "neuron", pct_thy1, pct_microglia),
    pct_other    = dplyr::if_else(cell_type == "neuron", pct_microglia, pct_thy1),
    detection_gap = pct_expected - pct_other,
    detection_expected_way = dplyr::if_else(cell_type == "neuron", est_thy1_vs_mg > 0, est_thy1_vs_mg < 0),
    log2fc_toward_expected = dplyr::if_else(cell_type == "neuron", logFC, -logFC),
    abundance_expected_way = log2fc_toward_expected > 0,
    common_protein = det_thy1_pct >= 50 & det_mg_pct >= 50,
    shown_in_Fig3a = (cell_type == "neuron" & gene %in% neuron_markers) |
      (cell_type == "microglia" & gene %in% microglia_markers)
  )
# q values: the SAME proteome-wide tests as the volcano (3b) and the heatmap set (ED Fig. 5): BH over all tested proteins
# (section 3), not over the candidates only, so a candidate is "clear" by one rule everywhere. NA = not in that test
# (detection test: detected in >= 10% of the cells of either type; abundance test: see section 3.1).
cand_all <- cand_all %>%
  dplyr::left_join(de_detect %>% dplyr::select(gene, detection_q = adj.P.Val), by = "gene") %>%
  dplyr::left_join(de_ct %>% dplyr::select(gene, abundance_q = adj.P.Val), by = "gene")
write_csv(cand_all, file.path(DIR_F3, "Figure3c_candidates_detection_and_abundance.csv"))

summarise_readout <- function(d, readout) {
  way <- if (readout == "detection") d$detection_expected_way else d$abundance_expected_way
  q   <- if (readout == "detection") d$detection_q else d$abundance_q
  eff <- if (readout == "detection") d$detection_gap else d$log2fc_toward_expected
  ok <- !is.na(way)
  tibble::tibble(readout = readout, n_tested = sum(ok), n_expected_way = sum(way[ok]), pct_expected_way = round(100 * mean(way[ok]), 1),
                 median_effect = round(median(eff[ok]), 2),     # detection: percentage points; abundance: log2 fold change
                 n_clear_q05 = sum(way[ok] & q[ok] < 0.05, na.rm = TRUE))
}
cand_summary <- purrr::map_dfr(c("neuron", "microglia"), function(ct) {
  d <- cand_all[cand_all$cell_type == ct, ]
  dplyr::bind_rows(
    summarise_readout(d, "detection") %>% dplyr::mutate(group = "all candidates"),
    summarise_readout(d, "abundance") %>% dplyr::mutate(group = "all candidates"),
    summarise_readout(d[d$common_protein %in% TRUE, ], "detection") %>% dplyr::mutate(group = "detected in >=50% of cells of both types"),
    summarise_readout(d[d$common_protein %in% TRUE, ], "abundance") %>% dplyr::mutate(group = "detected in >=50% of cells of both types"),
    summarise_readout(d[!(d$common_protein %in% TRUE), ], "detection") %>% dplyr::mutate(group = "others"),
    summarise_readout(d[!(d$common_protein %in% TRUE), ], "abundance") %>% dplyr::mutate(group = "others")
  ) %>% dplyr::mutate(cell_type = ct, .before = 1)
})
write_csv(cand_summary, file.path(DIR_F3, "Figure3c_candidates_expected_way_summary.csv"))
print(as.data.frame(cand_summary), row.names = FALSE)

# Panel c: candidates tested in BOTH readouts
cand_both <- cand_all %>% dplyr::filter(!is.na(detection_expected_way), !is.na(abundance_expected_way))
cand_both_summary <- purrr::map_dfr(c("neuron", "microglia"), function(ct) {
  d <- cand_both[cand_both$cell_type == ct, ]
  dplyr::bind_rows(summarise_readout(d, "detection"), summarise_readout(d, "abundance")) %>%
    dplyr::mutate(cell_type = ct, group = "candidates tested in both readouts", .before = 1)
})
write_csv(cand_both_summary, file.path(DIR_F3, "Figure3c_data.csv"))
cand_bars <- cand_both_summary %>%
  dplyr::mutate(cell_type = factor(cell_type, c("neuron", "microglia"), c("Neurons", "Microglia")),
                readout = factor(tools::toTitleCase(readout), c("Detection", "Abundance")),
                label = paste0(n_expected_way, "/", n_tested))
cand_bars$key <- paste(cand_bars$cell_type, cand_bars$readout)
cand_fill <- c(`Neurons Detection` = PAL_THY1, `Neurons Abundance` = colorspace::lighten(PAL_THY1, 0.5),
               `Microglia Detection` = PAL_MICROGLIA, `Microglia Abundance` = colorspace::lighten(PAL_MICROGLIA, 0.5))
# Dark bar = detection, light bar = abundance (see the figure legend).
p3c_candidates <- ggplot(cand_bars, aes(cell_type, pct_expected_way, fill = key, group = readout)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 50, linetype = 2, linewidth = 0.3) +
  geom_text(aes(label = label), position = position_dodge(0.8), angle = 90, hjust = -0.12, vjust = 0.5, size = pt_to_geom(PT_MIN)) +
  scale_fill_manual(values = cand_fill, guide = "none") +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25), expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "Candidates higher in\nexpected cell type (%)") +
  THEME_JOURNAL() + theme(axis.text.x = element_text(face = "bold"))

###############################################################################
# 8. FIG. 3d, 3e: OVERLAP WITH PUBLISHED DATASETS (UpSet)
# d: neurons (this study vs Sharma 2015, Rayaprolu 2022, Johnson 2026); e: microglia (this study vs Sharma
# 2015, Zhang 2026 = Neuron 2026 Table S3, Lloyd 2024).
# How proteins are matched (same rules for both plots):
#  * by UniProt accession first (any accession shared), gene symbol second. Tables that only have gene
#    symbols (Zhang 2026) get accessions through a symbol -> accession dictionary built from all other
#    tables, so old/new symbols (IARS vs IARS1, SEPT9 vs SEPTIN9) no longer count as different proteins.
#  * reference lists are taken as the publication reports them; only Rayaprolu 2022 is filtered
#    (enriched >= 2-fold, FDR < 0.05), because its table lists everything quantified.
#  * this study: a protein counts as detected if seen in >= MIN_CELLS_UPSET cells of the analysis set.
###############################################################################
message("\n[Figure 3d, 3e] Overlap with published datasets (UpSet, matched by UniProt accession)...")

study_entries <- function(cols) {
  det_n <- rowSums(!is.na(protein_mat_linear[, cols, drop = FALSE]))
  keep  <- names(det_n)[det_n >= MIN_CELLS_UPSET]
  ann   <- master_data[match(keep, master_data$gene_label), ]
  list(entries = data.frame(acc = ann$Protein.Group, gene = keep, stringsAsFactors = FALSE),
       n_ge1 = sum(det_n >= 1), n_kept = length(keep), det_n = det_n)
}

make_overlap_upset <- function(df, study_label, fill_col, set_fill) {
  others <- setdiff(setdiff(names(df), "protein"), study_label)
  ints <- unlist(lapply(0:length(others), function(k) combn(others, k, simplify = FALSE)), recursive = FALSE)
  ints <- lapply(ints, function(cmb) c(study_label, cmb))               # only intersections that include this study
  ComplexUpset::upset(
    df, intersect = c(study_label, others), intersections = ints,
    name = "Dataset", width_ratio = 0.2, height_ratio = 0.35,
    base_annotations = list("Intersection size" = ComplexUpset::intersection_size(
      fill = fill_col, bar_number_threshold = 1.01, text = list(size = pt_to_geom(PT_MIN))) +
        ggplot2::scale_y_continuous(limits = c(0, 3000), breaks = seq(0, 3000, 1000), expand = ggplot2::expansion(mult = c(0, 0)))),   # same y axis (0-3000) in both UpSets
    themes = ComplexUpset::upset_modify_themes(list(
      "intersections_matrix" = ggplot2::theme(text = ggplot2::element_text(size = PT_BASE),
                                              axis.text.y = ggplot2::element_text(size = PT_BASE),
                                              axis.title = ggplot2::element_text(size = PT_BASE, face = "bold")),
      "overall_sizes"        = ggplot2::theme(text = ggplot2::element_text(size = PT_BASE),
                                              axis.text = ggplot2::element_text(size = PT_BASE),
                                              axis.title = ggplot2::element_text(size = PT_BASE, face = "bold")),
      "Intersection size"    = ggplot2::theme(text = ggplot2::element_text(size = PT_BASE),
                                              axis.text = ggplot2::element_text(size = PT_BASE),
                                              axis.title = ggplot2::element_text(size = PT_BASE, face = "bold")))),
    set_sizes = ComplexUpset::upset_set_size(geom = ggplot2::geom_bar(fill = set_fill, width = 0.6)) +
      ggplot2::ylab("Total protein\ncount") +
      ggplot2::scale_y_continuous(breaks = 2000, expand = ggplot2::expansion(mult = c(0, 0.06))))
}

# ---- Fig. 3e: microglia -------------------------------------------------------------------------------
mg_upset_cols <- cohort_meta %>%
  dplyr::filter(tier == "tier_1x_single_cell", cell_type_label == "Microglia",
                !is.na(n_proteins_filtered), n_proteins_filtered >= MIN_PROTEINS) %>% dplyr::pull(sample_col)
mg_study <- study_entries(mg_upset_cols)

.v <- suppressWarnings(as.numeric(sharma_raw[[sharma_mg_col]])); .k <- is.finite(.v)
sharma_mg_entries <- data.frame(acc = sharma_raw[["Majority protein IDs"]][.k], gene = sharma_raw[["Gene names"]][.k])

# Lloyd 2024: ex vivo microglia, detected in >= 1 of the samples (the candidate rule in section 4 requires >= 3 of 4)
exmm_cols_lloyd <- colnames(lloyd_raw)[stringr::str_detect(colnames(lloyd_raw), "exMM")]
.n <- rowSums(sapply(lloyd_raw[exmm_cols_lloyd], function(x) suppressWarnings(as.numeric(x)) > 0), na.rm = TRUE)
lloyd_entries <- data.frame(acc = lloyd_raw[["Majority protein IDs"]][.n >= 1], gene = lloyd_raw[["Gene names"]][.n >= 1])

zhang_entries <- data.frame(acc = NA_character_, gene = unique(stats::na.omit(as.character(s3_raw[["Gene_name"]]))))

mg_sets <- list("3D-SCP" = mg_study$entries, "Sharma et al. 2015" = sharma_mg_entries,
                "Zhang et al. 2026" = zhang_entries, "Lloyd et al. 2024" = lloyd_entries)
upset_df_mg <- build_overlap_table(mg_sets)
message("  Microglia UpSet: this study ", mg_study$n_kept, " proteins in >=", MIN_CELLS_UPSET, " cells (", mg_study$n_ge1, " in >=1 cell) | set sizes ",
        paste(names(mg_sets), colSums(upset_df_mg[names(mg_sets)]), sep = "=", collapse = ", "))
write_csv(upset_df_mg, file.path(DIR_F3, "Figure3e_microglia_literature_overlap_data.csv"))
p3e_upset_microglia <- make_overlap_upset(upset_df_mg, "3D-SCP", PAL_MICROGLIA, PAL_ACCENT)
message("  Fig. 3e saved (n=", nrow(upset_df_mg), " protein entries across all 4 sets)")

mg_study_rows <- upset_df_mg[upset_df_mg[["3D-SCP"]], ]
mg_other_cols <- setdiff(names(mg_sets), "3D-SCP")
n_mg_study_total   <- nrow(mg_study_rows)
n_mg_study_only    <- sum(rowSums(mg_study_rows[mg_other_cols]) == 0)
n_mg_overlap_any   <- n_mg_study_total - n_mg_study_only
n_mg_all_four      <- sum(rowSums(mg_study_rows[mg_other_cols]) == length(mg_other_cols))
mg_merge <- merged_entries_report(upset_df_mg)
message("  Microglia UpSet entries: ", mg_merge$n_entries, " -> ", mg_merge$n_proteins, " proteins (", mg_merge$n_dropped, " dropped without accession/gene; ", nrow(mg_merge$merged), " entries merged into shared proteins)")
if (nrow(mg_merge$merged) > 0) write_csv(mg_merge$merged, file.path(DIR_F3, "Figure3e_microglia_merged_entries.csv"))
message("  Fig 3e headline: ", n_mg_overlap_any, "/", n_mg_study_total, " (", round(100 * n_mg_overlap_any / n_mg_study_total, 1),
        "%) of 3D-SCP microglia proteins also in >=1 other dataset; ", n_mg_all_four, "/", n_mg_study_total, " (",
        round(100 * n_mg_all_four / n_mg_study_total, 1), "%) in all 4 datasets")

# ---- Fig. 3d: neurons ---------------------------------------------------------------------------------
thy1_upset_cols <- cohort_meta %>%
  dplyr::filter(tier=="tier_1x_single_cell", cell_type_label=="Thy1",
                !is.na(n_proteins_filtered), n_proteins_filtered >= MIN_PROTEINS) %>% dplyr::pull(sample_col)
neu_study <- study_entries(thy1_upset_cols)

.v <- suppressWarnings(as.numeric(sharma_raw[[sharma_ne_col]])); .k <- is.finite(.v)
sharma_ne_entries <- data.frame(acc = sharma_raw[["Majority protein IDs"]][.k], gene = sharma_raw[["Gene names"]][.k])

# Johnson et al. 2026 (Anal Chem), Table S3: DA | PVN | 5-HT blocks side by side (UniProt acc. in cols 1/5/9, gene symbol in 3/7/11)
johnson_entries <- dplyr::bind_rows(lapply(list(c(1, 3), c(5, 7), c(9, 11)), function(j)
  data.frame(acc = as.character(johnson_raw[[j[1]]][-1]), gene = as.character(johnson_raw[[j[2]]][-1]), stringsAsFactors = FALSE)))
johnson_entries <- johnson_entries[!(is.na(johnson_entries$acc) & is.na(johnson_entries$gene)), ]

# Rayaprolu et al. 2022 (Nat Commun), Suppl. Data 9: Camk2a-TurboID neuron-labelled = enriched >=2-fold, BH-FDR < 0.05
rayaprolu_entries <- data.frame(acc = ray_raw[["Majority Protein IDs"]][ray_enriched], gene = ray_raw[["Gene Symbols"]][ray_enriched])

neu_sets <- list("3D-SCP" = neu_study$entries, "Sharma et al. 2015" = sharma_ne_entries,
                 "Rayaprolu et al. 2022" = rayaprolu_entries, "Johnson et al. 2026" = johnson_entries)
upset_df_neu <- build_overlap_table(neu_sets)
message("  Neuron UpSet: this study ", neu_study$n_kept, " proteins in >=", MIN_CELLS_UPSET, " cells (", neu_study$n_ge1, " in >=1 cell) | set sizes ",
        paste(names(neu_sets), colSums(upset_df_neu[names(neu_sets)]), sep = "=", collapse = ", "))
write_csv(upset_df_neu, file.path(DIR_F3, "Figure3d_neuron_literature_overlap_data.csv"))
p3d_upset_neuron <- make_overlap_upset(upset_df_neu, "3D-SCP", PAL_THY1, PAL_SHARMA)
message("  Fig. 3d saved (n=", nrow(upset_df_neu), " protein entries across ", length(neu_sets), " sets)")

neu_study_rows <- upset_df_neu[upset_df_neu[["3D-SCP"]], ]
neu_other_cols <- setdiff(names(neu_sets), "3D-SCP")
n_neu_study_total <- nrow(neu_study_rows)
n_neu_study_only  <- sum(rowSums(neu_study_rows[neu_other_cols]) == 0)
n_neu_overlap_any <- n_neu_study_total - n_neu_study_only
n_neu_all         <- sum(rowSums(neu_study_rows[neu_other_cols]) == length(neu_other_cols))
neu_merge <- merged_entries_report(upset_df_neu)
message("  Neuron UpSet entries: ", neu_merge$n_entries, " -> ", neu_merge$n_proteins, " proteins (", neu_merge$n_dropped, " dropped without accession/gene; ", nrow(neu_merge$merged), " entries merged into shared proteins)")
if (nrow(neu_merge$merged) > 0) write_csv(neu_merge$merged, file.path(DIR_F3, "Figure3d_neuron_merged_entries.csv"))
message("  Fig 3d headline: ", n_neu_overlap_any, "/", n_neu_study_total, " (", round(100 * n_neu_overlap_any / n_neu_study_total, 1),
        "%) of 3D-SCP neuron proteins also in >=1 other dataset; ", n_neu_all, "/", n_neu_study_total, " (",
        round(100 * n_neu_all / n_neu_study_total, 1), "%) in all ", length(neu_sets), " datasets")

###############################################################################
# 9. MAIN FIGURE 3: ASSEMBLY
#   row 1:  a detection (10 of 24 columns) | b volcano (9) | c candidate summary (5)
#   row 2:  d neuron UpSet (12)            | e microglia UpSet (12)
# Row 1 plotting areas are exactly FIG3_PANEL_H_IN high (fix_panel_h). The UpSet panels are ComplexUpset
# objects (small patchworks themselves) and are used as single cells; their text size is set inside
# make_overlap_upset(). Check in the SVG that all text is 7 pt.
###############################################################################
message("\n[Figure 3] Assembling the main figure...")
fig3_design <- "
AAAAAAAAAABBBBBBBBBCCCCC
DDDDDDDDDDDDEEEEEEEEEEEE
"
text_7pt <- ggplot2::theme(text = ggplot2::element_text(size = PT_BASE), axis.text = ggplot2::element_text(size = PT_BASE),
                           axis.title = ggplot2::element_text(size = PT_BASE, face = "bold"))
fig3_row1_panels <- list(fix_panel_h(p3a_detection & text_7pt, FIG3_PANEL_H_IN),
                         fix_panel_h(p3b_volcano & text_7pt, FIG3_PANEL_H_IN),
                         fix_panel_h(p3c_candidates & text_7pt, FIG3_PANEL_H_IN))
fig3_main_assembled <- (fig3_row1_panels[[1]] + fig3_row1_panels[[2]] + fig3_row1_panels[[3]] +
                          p3d_upset_neuron + p3e_upset_microglia) +
  patchwork::plot_layout(design = fig3_design, heights = FIG3_ROW_H_IN) & text_7pt
ggsave(file.path(DIR_F3, "Figure3_MAIN_assembled.svg"), fig3_main_assembled,
       width = FIG3_W_IN, height = sum(FIG3_ROW_H_IN), units = "in", device = SVG_DEVICE)
message("  Figure 3 saved (", round(FIG3_W_IN, 2), " x ", sum(FIG3_ROW_H_IN), " in) -- insert at 100%, do not rescale.")

###############################################################################
# 10. ED FIG. 4a: MISSINGNESS BY ABUNDANCE QUINTILE (purely descriptive; no test or model)
# The mean fraction of missing values falls with protein abundance (mean per abundance quintile), the
# pattern expected for values missing because the signal is below detection: the reason nothing is
# imputed. Scoped to exactly the population of the DE analysis: single cells with >= MIN_PROTEINS
# contaminant-filtered protein groups, Thy1 + microglia only, all regions pooled.
###############################################################################
message("\n[ED Fig. 4a] Missingness by abundance quintile (descriptive)...")

# The same cells as the DE analysis: single-cell tier, >= MIN_PROTEINS, Thy1 + microglia only (no control wells).
missingness_diag_cols <- cohort_meta %>%
  dplyr::filter(tier == "tier_1x_single_cell", cell_type_label %in% c("Thy1", "Microglia"),
                !is.na(n_proteins_filtered), n_proteins_filtered >= MIN_PROTEINS) %>%   # analysis set: contaminant-filtered count >= 600 (same cells as the DE)
  dplyr::pull(sample_col)
missingness_diag_cols <- intersect(missingness_diag_cols, colnames(protein_mat_log2))
mat_missing_diag <- protein_mat_log2[, missingness_diag_cols, drop = FALSE]
message("  Using ", ncol(mat_missing_diag), " analysis-set single cells (>=600 contaminant-filtered proteins; single-cell tier only; ",
        "Thy1 + Microglia only; all regions pooled) x ", nrow(mat_missing_diag), " contamination-filtered proteins")

frac_missing <- rowMeans(is.na(mat_missing_diag))
mean_abundance_observed <- apply(mat_missing_diag, 1, function(x) mean(x, na.rm = TRUE))

# drop proteins with 0% or 100% missingness -- undefined/trivial for this diagnostic
keep_diag <- frac_missing > 0 & frac_missing < 1 & is.finite(mean_abundance_observed)
diag_df <- tibble::tibble(
  gene = rownames(mat_missing_diag)[keep_diag],
  mean_log2_abundance_observed = mean_abundance_observed[keep_diag],
  frac_missing = frac_missing[keep_diag]
)
message("  n=", nrow(diag_df), " proteins with intermediate missingness (excludes always-detected/never-detected)")
n_always_det <- sum(frac_missing == 0, na.rm = TRUE)    # detected in every analysis-set cell
n_never_det  <- sum(frac_missing == 1, na.rm = TRUE)    # detected in none of the analysis-set cells
message("  always detected: ", n_always_det, " | never detected: ", n_never_det, " | intermediate: ", nrow(diag_df))

diag_df <- diag_df %>% dplyr::mutate(abundance_quintile = dplyr::ntile(mean_log2_abundance_observed, 5))
quintile_summary <- diag_df %>%
  dplyr::group_by(abundance_quintile) %>%
  dplyr::summarise(mean_frac_missing = mean(frac_missing), n_proteins = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(abundance_quintile_label = factor(abundance_quintile,
                                                  labels = c("Q1 (low)","Q2","Q3","Q4","Q5 (high)")))
write_csv(quintile_summary, file.path(DIR_F3_SUPP, "EDFig4a_missingness_by_abundance_quintile.csv"))
print(quintile_summary)
message("  ED Fig. 4a quintile range: Q1 (lowest abundance) = ",
        round(100 * quintile_summary$mean_frac_missing[quintile_summary$abundance_quintile == 1], 1),
        "% missing (n=", quintile_summary$n_proteins[quintile_summary$abundance_quintile == 1],
        ")  |  Q5 (highest abundance) = ",
        round(100 * quintile_summary$mean_frac_missing[quintile_summary$abundance_quintile == 5], 1),
        "% missing (n=", quintile_summary$n_proteins[quintile_summary$abundance_quintile == 5], ")")

ed4a_missingness <- ggplot(quintile_summary, aes(abundance_quintile_label, mean_frac_missing)) +
  geom_col(fill = PAL_ACCENT, width = 0.6) +   # both cell types pooled, so the neutral colour
  geom_text(aes(label = paste0(round(mean_frac_missing * 100, 1), "%\n(n=", n_proteins, ")")),
            vjust = -0.3, size = pt_to_geom(PT_MIN), lineheight = 0.9) +
  labs(x = "Abundance quintile (observed values)", y = "Mean fraction missing") +
  scale_y_continuous(limits = c(0, max(quintile_summary$mean_frac_missing) * 1.25),
                     expand = expansion(mult = c(0, 0))) +   # bars start at 0; headroom for the labels is in the 1.25x limit
  THEME_JOURNAL() + theme(axis.text.x = element_text(face = "bold", angle = 30, hjust = 1))

write_csv(diag_df, file.path(DIR_F3_SUPP, "EDFig4a_missingness_per_protein.csv"))

###############################################################################
# 11. ED FIG. 4b-d: TISSUE-FREE CONTROL (PEN_SC) WELLS
# b  protein groups detected per control well (all control wells, incl. the outlier well; well names are
#    replaced by numbers 1-n, ordered by protein count)
# c  detection of the markers of panels a/b in the control wells used for the screen
# d  the background screen of the whole proteome: a protein is flagged as background if it is seen in
#    >= 5 of 7 control wells AND its median intensity in the controls is >= the 69th percentile of all
#    signals (thresholds PEN_* in setup.R; the screen itself is run in 00_build_master_dataset.R and
#    its flags feed the contaminant filter). The wells used (pen_sc_cols) are defined in setup.R.
###############################################################################
message("\n[ED Fig. 4b-d] Tissue-free control wells...")

# ---- ED 4b: one bar per well --------------------------------------------------------------------------
pen_wells_df <- cohort_meta %>%
  dplyr::filter(is_pen_sc, is.finite(n_proteins_detected)) %>%
  dplyr::mutate(well_label = paste0(dplyr::coalesce(well_row, ""), dplyr::coalesce(as.character(well_col), "")))
if (nrow(pen_wells_df) == 0) stop("No PEN_SC control wells found in cohort_meta (column is_pen_sc).", call. = FALSE)
pen_wells_df <- pen_wells_df %>%
  dplyr::arrange(n_proteins_detected) %>%
  dplyr::mutate(well_number = factor(dplyr::row_number(), levels = dplyr::row_number()))
message("  PEN_SC per well: n=", nrow(pen_wells_df), " control wells, range ", min(pen_wells_df$n_proteins_detected), "-",
        max(pen_wells_df$n_proteins_detected), " proteins (well identities anonymized to 1-", nrow(pen_wells_df),
        " on the plot; true labels: ", paste(pen_wells_df$well_label, collapse = ", "), ")")
ed4b_pen_wells <- ggplot(pen_wells_df, aes(well_number, n_proteins_detected)) +
  geom_col(width = 0.6, fill = PAL_ACCENT) +   # neutral colour: control wells
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +   # bars start at 0
  coord_flip() +
  labs(x = "Tissue/cell-free control well", y = "Protein groups detected") +
  THEME_JOURNAL()

# ---- ED 4c: marker detection in the control wells -------------------------------------------------------
# Only markers with non-zero control detection are drawn (the others are 0% and invisible); the full list,
# including the zeros, is in the CSV.
pen_marker_detection <- purrr::map_dfr(all_markers, function(g) {
  if (!g %in% rownames(protein_mat_log2)) return(NULL)
  tibble::tibble(gene = g,
                 marker_set = dplyr::if_else(g %in% neuron_markers, "Neuronal markers", "Microglial markers"),
                 pct_pen = 100 * mean(!is.na(protein_mat_log2[g, pen_sc_cols])))
})
write_csv(pen_marker_detection, file.path(DIR_F3_SUPP, "EDFig4c_marker_detection_PEN_controls.csv"))
pen_marker_shown <- pen_marker_detection %>% dplyr::filter(pct_pen > 0)
if (nrow(pen_marker_shown) == 0) stop("None of the markers is detected in the control wells: ED Fig. 4c would be empty.", call. = FALSE)
ed4c_pen_markers <- ggplot(pen_marker_shown %>% dplyr::mutate(gene = factor(gene, levels = marker_order)),
                           aes(pct_pen, gene, fill = marker_set)) +
  geom_col(width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c("Neuronal markers" = PAL_THY1, "Microglial markers" = PAL_MICROGLIA), guide = "none") +
  facet_grid(marker_set ~ ., scales = "free_y", space = "free_y",
             labeller = ggplot2::as_labeller(c("Neuronal markers" = "Neuronal\nmarkers", "Microglial markers" = "Microglial\nmarkers"))) +
  # axis labels as k/n controls (pct_pen is always a multiple of 100/n)
  scale_x_continuous(
    breaks = seq_along(pen_sc_cols) / length(pen_sc_cols) * 100,
    labels = paste0(seq_along(pen_sc_cols), "/", length(pen_sc_cols)),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(x = "Controls with protein detected", y = NULL) +
  THEME_JOURNAL() + theme(panel.spacing.y = grid::unit(1.5, "mm"))

# ---- ED 4d: background screen of the whole proteome -------------------------------------------------------
PEN_SCREEN_FILE <- file.path(MASTER_DATASET_DIR, "pen_sc_systematic_screen_ALL_PROTEINS.csv")
require_file(PEN_SCREEN_FILE, "PEN_SC contaminant screen table", "It is written by 00_build_master_dataset.R (step 2); rerun that script.")
pen_screen_full <- read_csv(PEN_SCREEN_FILE, show_col_types = FALSE)
require_columns(pen_screen_full, c("gene_label", "n_pen_detected", "n_pen_total", "pen_detection_frac", "pen_percentile",
                                   "meets_pen_contamination_bar"), "pen_sc_systematic_screen_ALL_PROTEINS.csv")
# The wells and thresholds of this script must be the ones behind the flags in the screen file
if (max(pen_screen_full$n_pen_total, na.rm = TRUE) != length(pen_sc_cols)) {
  stop("The screen in ", basename(PEN_SCREEN_FILE), " used ", max(pen_screen_full$n_pen_total, na.rm = TRUE),
       " control wells, but setup.R selects ", length(pen_sc_cols), " (PEN_EXCLUDED_WELL_PATTERN).", call. = FALSE)
}
flag_recomputed <- !is.na(pen_screen_full$pen_percentile) &
  pen_screen_full$pen_detection_frac >= PEN_DETECTION_FRAC_MIN & pen_screen_full$pen_percentile >= PEN_PERCENTILE_MIN
if (any(flag_recomputed != pen_screen_full$meets_pen_contamination_bar, na.rm = TRUE)) {
  stop("PEN_DETECTION_FRAC_MIN / PEN_PERCENTILE_MIN in setup.R do not reproduce the flags in ", basename(PEN_SCREEN_FILE),
       ". Make them equal to the values in 00_build_master_dataset.R.", call. = FALSE)
}
n_pen_total_proteins <- nrow(pen_screen_full)
n_pen_ever_detected  <- sum(pen_screen_full$n_pen_detected > 0, na.rm = TRUE)
n_pen_flagged        <- sum(pen_screen_full$meets_pen_contamination_bar, na.rm = TRUE)
message("  PEN_SC screen: ", n_pen_total_proteins, " proteins screened, ", n_pen_ever_detected, " seen in at least one control (",
        round(100 * n_pen_ever_detected / n_pen_total_proteins, 1), "%), ", n_pen_flagged, " flagged (",
        round(100 * n_pen_flagged / n_pen_total_proteins, 1), "% of all, ", round(100 * n_pen_flagged / n_pen_ever_detected, 1),
        "% of those ever detected)")

# Median intensity in the controls (y axis): taken from master_data BEFORE the contaminant filter, because
# protein_mat_log2 no longer contains the flagged proteins. The flags themselves come from 00_build_master_dataset.R.
mat_unfiltered <- master_data %>%
  dplyr::select(-dplyr::any_of(c("Protein.Group", "Protein.Names", "First.Protein.Description",
                                 "N.Sequences", "N.Proteotypic.Sequences")))
mat_unfiltered <- as.matrix(mat_unfiltered %>% dplyr::select(-gene_label)) %>% `rownames<-`(master_data$gene_label)
mat_unfiltered[mat_unfiltered == 0] <- NA
mat_unfiltered <- log2(mat_unfiltered)
pen_median_log2 <- apply(mat_unfiltered[, intersect(pen_sc_cols, colnames(mat_unfiltered)), drop = FALSE], 1,
                         function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else median(x) })
pen_screen_full <- pen_screen_full %>%
  dplyr::left_join(tibble::tibble(gene_label = names(pen_median_log2), pen_median_log2_intensity = unname(pen_median_log2)), by = "gene_label")
message("  ", sum(pen_screen_full$meets_pen_contamination_bar & is.finite(pen_screen_full$pen_median_log2_intensity), na.rm = TRUE), " / ",
        n_pen_flagged, " flagged proteins have a plottable intensity value")

# Intensity threshold: percentile of all signals before contaminant filtering, as in the screen, so the dashed line
# matches the flagged / not-flagged colours.
all_signal_vals <- as.vector(mat_unfiltered); all_signal_vals <- all_signal_vals[is.finite(all_signal_vals)]
pen_y_threshold <- stats::quantile(all_signal_vals, PEN_PERCENTILE_MIN / 100, na.rm = TRUE)

pen_plot_df <- pen_screen_full %>% dplyr::filter(n_pen_detected > 0, is.finite(pen_median_log2_intensity))
n_pen_controls_total <- max(pen_plot_df$n_pen_total, na.rm = TRUE)
# x axis: all fractions k/n (1/7 ... 7/7)
pen_frac_breaks <- seq_len(n_pen_controls_total) / n_pen_controls_total
pen_frac_labels <- paste0(seq_len(n_pen_controls_total), "/", n_pen_controls_total)
# The vertical threshold line sits halfway between 4/7 and 5/7, so jittered points never appear on the wrong side.
pen_frac_breaks_below <- pen_frac_breaks[pen_frac_breaks < PEN_DETECTION_FRAC_MIN]
pen_threshold_line_x <- if (length(pen_frac_breaks_below) > 0) (max(pen_frac_breaks_below) + PEN_DETECTION_FRAC_MIN) / 2 else PEN_DETECTION_FRAC_MIN - 0.05

ed4d_pen_screen <- ggplot(pen_plot_df, aes(pen_detection_frac, pen_median_log2_intensity, colour = meets_pen_contamination_bar)) +
  annotate("rect", xmin = pen_threshold_line_x, xmax = Inf, ymin = pen_y_threshold, ymax = Inf, fill = PAL_DANGER, alpha = 0.08) +
  geom_jitter(width = 0.025, height = 0, alpha = 0.35, size = 0.6, stroke = 0) +
  geom_hline(yintercept = pen_y_threshold, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = pen_threshold_line_x, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  scale_x_continuous(breaks = pen_frac_breaks, labels = pen_frac_labels) +
  scale_colour_manual(values = c(`TRUE` = PAL_DANGER, `FALSE` = "grey55"),
                      labels = c(`TRUE` = "Meets background-signal criteria", `FALSE` = "Below background-signal threshold"), name = NULL) +
  labs(x = "Controls with protein detected (fraction)", y = "Median log2 intensity in controls") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2, alpha = 1))) +   # legend dots larger than the small, transparent plot points
  THEME_JOURNAL() + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1), legend.background = element_blank(),
                          panel.grid.minor.x = element_blank(), legend.margin = margin(0, 0, 0, 0))

###############################################################################
# 12. ED FIG. 4e: RANK-ABUNDANCE OF MICROGLIAL MARKERS vs INDEPENDENT DATASETS
# Mean log2 abundance per protein (>= 200-protein cohorts of setup.R: microglia single cells, 10-cell
# areas) next to the Sharma 2015 microglia column, each dataset ranked and scaled 0-1 on its own:
# only the rank position within a dataset is meaningful, not the absolute values across datasets.
# Markers not detected in a dataset are drawn as open circles in a "not detected" column to the right
# of the 0-100% axis (their height carries no information).
###############################################################################
message("\n[ED Fig. 4e] Rank-abundance, microglia...")

sharma_microglia <- sharma_raw %>%
  dplyr::select(gene_raw = `Gene names`, intensity = dplyr::all_of(sharma_mg_col)) %>%
  dplyr::filter(!is.na(intensity)) %>%
  dplyr::mutate(gene = stringr::str_extract(gene_raw, "^[^;]+")) %>%
  dplyr::filter(!is.na(gene), gene != "") %>%
  dplyr::distinct(gene, .keep_all = TRUE)
# The Sharma column is already log2 LFQ intensity (range ~21-39), like our log2 values: used as it is, no further transform.
message("  Sharma microglia column '", sharma_mg_col, "': range ", signif(min(sharma_microglia$intensity), 3), " to ",
        signif(max(sharma_microglia$intensity), 3))
if (max(sharma_microglia$intensity) > 100) stop("The Sharma column is not on a log2 scale (maximum > 100); check SHARMA_* in setup.R.", call. = FALSE)

build_rank_curve <- function(values_named, dataset_label) {
  values_named <- values_named[is.finite(values_named)]
  ord <- order(values_named, decreasing = TRUE)
  values_named <- values_named[ord]
  n <- length(values_named)
  rank_pct <- 100 * (seq_len(n) - 1) / (n - 1)
  val_norm <- (values_named - min(values_named)) / (max(values_named) - min(values_named))
  tibble::tibble(gene = names(values_named), dataset = dataset_label,
                 rank = seq_len(n), rank_pct = rank_pct, value_norm = val_norm)
}

# Markers not detected in a dataset are drawn as hollow circles right of the dotted line, one column per dataset
# (single cells, 10-cell areas, Sharma 2015) and one row per marker. The row height carries no information.
SENTINEL_X <- c(104, 108.5, 113)      # one x position per dataset (same order as the dataset levels)
SENTINEL_LABEL_X <- 116
add_not_detected_sentinel <- function(rank_curves, full_marker_panel, dataset_levels) {
  detected_combos <- rank_curves %>% dplyr::filter(gene %in% full_marker_panel) %>%
    dplyr::distinct(gene, dataset) %>% dplyr::mutate(detected = TRUE)
  missing <- tidyr::expand_grid(gene = full_marker_panel, dataset = dataset_levels) %>%
    dplyr::mutate(dataset = factor(dataset, levels = dataset_levels)) %>%
    dplyr::left_join(detected_combos, by = c("gene", "dataset")) %>%
    dplyr::filter(is.na(detected))
  genes_missing <- full_marker_panel[full_marker_panel %in% missing$gene]          # panel order
  row_y <- setNames(seq(0.92, 0.05, length.out = max(2, length(genes_missing)))[seq_along(genes_missing)], genes_missing)
  missing %>% dplyr::arrange(dataset, gene) %>%
    dplyr::mutate(rank_pct = SENTINEL_X[as.integer(dataset)], value_norm = unname(row_y[gene]))
}

vals_1x  <- rowMeans(protein_mat_log2[, mg_singlecell_cols, drop=FALSE], na.rm=TRUE)
vals_10x <- rowMeans(protein_mat_log2[, thy1_bulk10_cols, drop=FALSE], na.rm=TRUE)
vals_sharma <- setNames(sharma_microglia$intensity, sharma_microglia$gene)

curve_1x     <- build_rank_curve(vals_1x[is.finite(vals_1x)], "Single cells (Cx3cr1-EGFP microglia)")
curve_10x    <- build_rank_curve(vals_10x[is.finite(vals_10x)], "10-cell areas (unsorted)")
curve_sharma <- build_rank_curve(vals_sharma, "Sharma et al. 2015 (MACS microglia)")
rank_curves <- dplyr::bind_rows(curve_1x, curve_10x, curve_sharma) %>%
  dplyr::mutate(dataset = factor(dataset, levels = c("Single cells (Cx3cr1-EGFP microglia)",
                                                     "10-cell areas (unsorted)",
                                                     "Sharma et al. 2015 (MACS microglia)")))

rank_marker_positions <- rank_curves %>% dplyr::filter(gene %in% MICROGLIA_RANK_PANEL)
not_in_any <- setdiff(MICROGLIA_RANK_PANEL, unique(rank_marker_positions$gene))
write_csv(rank_marker_positions, file.path(DIR_F3_SUPP, "EDFig4e_marker_positions_across_datasets.csv"))
message("  Markers not detected in ANY of the 3 datasets: ",
        if (length(not_in_any) > 0) paste(not_in_any, collapse=", ") else "none")

not_detected_mg <- add_not_detected_sentinel(rank_curves, MICROGLIA_RANK_PANEL, levels(rank_curves$dataset))
write_csv(not_detected_mg %>% dplyr::select(gene, dataset), file.path(DIR_F3_SUPP, "EDFig4e_markers_not_detected_per_dataset.csv"))
not_detected_labels <- not_detected_mg %>% dplyr::distinct(gene, value_norm)        # one label per marker
message("  ED 4e: ", nrow(not_detected_labels), " markers are missing in at least one dataset (one label each; ", nrow(not_detected_mg), " hollow circles)")

dataset_colours <- c("Single cells (Cx3cr1-EGFP microglia)"=PAL_MICROGLIA,
                     "10-cell areas (unsorted)"=PAL_BULK,
                     "Sharma et al. 2015 (MACS microglia)"=PAL_SHARMA)

# Short legend entries; the full description is in the figure legend.
dataset_labels <- c("Single cells (Cx3cr1-EGFP microglia)" = "Single microglia", "10-cell areas (unsorted)" = "10-cell areas",
                    "Sharma et al. 2015 (MACS microglia)" = "Sharma 2015 (MACS)")
ed4e_rank_abundance <- ggplot(rank_curves, aes(rank_pct, value_norm, colour = dataset)) +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  geom_point(data = rank_marker_positions, size = 1.5) +
  geom_point(data = not_detected_mg, shape = 1, size = 1.5, stroke = 0.6) +
  ggrepel::geom_text_repel(data = rank_marker_positions, aes(label = gene), size = pt_to_geom(PT_MIN),
                           max.overlaps = Inf, show.legend = FALSE,
                           force = 3, max.iter = 20000, box.padding = 0.25,
                           min.segment.length = 0, seed = 1) +
  # one label per marker that is not detected in at least one dataset
  geom_text(data = not_detected_labels, aes(x = SENTINEL_LABEL_X, y = value_norm, label = gene), hjust = 0,
            size = pt_to_geom(PT_MIN), colour = "grey30", inherit.aes = FALSE) +
  geom_vline(xintercept = 100, linetype = "dotted", colour = "grey60") +
  annotate("text", x = 102.5, y = 1.08, label = "not detected", hjust = 0, size = pt_to_geom(PT_MIN), colour = "grey40") +
  scale_colour_manual(values = dataset_colours, labels = dataset_labels, name = NULL) +
  coord_cartesian(xlim = c(0, 128), clip = "off") +
  labs(x = "Abundance rank percentile within dataset (0% = most abundant)",
       y = "Normalised abundance") +
  guides(colour = guide_legend(nrow = 1)) +
  THEME_JOURNAL() +
  theme(
    legend.position = "bottom", legend.justification = "left", legend.background = element_blank(),
    legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-2, 0, 0, 0),
    plot.margin = ggplot2::margin(t = 12, r = 28, b = 4, l = 4)
  )

###############################################################################
# 13. ED FIG. 4f, 4g: DE METHOD COMPARISON FOR THE CELL-TYPE CONTRAST
# Wilcoxon (primary, section 3.1) vs proDA (dropout-aware, protein level) vs msqrob2 (precursor level, robust
# summarisation to protein). All three are non-imputed, use the same cells and the same contaminant filter,
# and no normalisation. limma + imputation is deliberately not shown (ED Fig. 4a shows why imputation is
# inappropriate here). The question is whether three independent non-imputed approaches agree.
#   f  volcano plots of the three methods (same marker highlighting as Fig. 3b, shared axes)
#   g  pairwise agreement of the log2 fold changes
# Sensitivity (table and summary numbers only): lm(log2 intensity ~ cell type + log2 depth) on the same cells
# and proteins as the Wilcoxon test. Batch is not in the model: it is completely confounded with cell type.
###############################################################################
message("\n[ED Fig. 4f-g] DE method comparison: Wilcoxon vs proDA vs msqrob2...")

# ---- proDA on the same matrix and cells as the Wilcoxon test --------------------------------------------
message("  Running proDA...")
run_proda_celltype <- function(mat_log, thy1_ids, mg_ids) {
  mat_sub <- mat_log[, c(thy1_ids, mg_ids)]
  group <- factor(ifelse(colnames(mat_sub) %in% thy1_ids, "Thy1", "Microglia"), levels = c("Microglia", "Thy1"))
  fit <- proDA::proDA(mat_sub, design = ~ group)
  res <- proDA::test_diff(fit, contrast = "groupThy1")
  tibble::tibble(gene = res$name, proda_logFC = res$diff, proda_pval = res$pval, proda_padj = res$adj_pval)
}
res_proda_ct <- run_proda_celltype(mat_de_log, thy1_ids, mg_ids)
message("    ", sum(res_proda_ct$proda_padj < 0.05, na.rm = TRUE), " / ", nrow(res_proda_ct), " proteins at FDR q<0.05")

# ---- msqrob2: precursor level (DIA-NN precursor matrix), robust summarisation to protein ----------------
# Same contaminant screen as the protein-level methods (clean_genes), no normalisation, precursors seen in >= 3 cells.
run_msqrob2_celltype <- function(pr_matrix_file, sample_metadata, thy1_ids, mg_ids, clean_genes) {
  message("    Reading precursor matrix (large file, may take a minute)...")
  pr_raw <- data.table::fread(pr_matrix_file, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
  meta_cols <- c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description",
                 "Proteotypic", "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge", "Precursor.Id")
  raw_file_cols <- setdiff(colnames(pr_raw), meta_cols)
  name_lookup <- setNames(raw_file_cols, basename(gsub("\\\\", "/", raw_file_cols)))
  
  target_ids <- c(thy1_ids, mg_ids)
  meta_sub <- sample_metadata %>% dplyr::filter(clean_sample_name %in% target_ids)
  meta_sub$raw_basename <- basename(gsub("\\\\", "/", meta_sub$raw_filename))
  matched <- meta_sub %>% dplyr::filter(raw_basename %in% names(name_lookup))
  message("    Matched ", nrow(matched), " / ", length(target_ids), " cells to pr_matrix columns.")
  if (nrow(matched) < length(target_ids)) {
    stop("msqrob2: ", length(target_ids) - nrow(matched), " cell(s) of the analysis set have no column in the precursor matrix. Fix the raw-file names in the metadata so that all cells match.", call. = FALSE)
  }
  
  pr_sub <- pr_raw[, c("Precursor.Id", "Genes", name_lookup[matched$raw_basename]), drop = FALSE]
  colnames(pr_sub) <- c("Precursor.Id", "Genes", matched$clean_sample_name)
  pr_sub <- pr_sub[!is.na(pr_sub$Genes) & pr_sub$Genes != "" & !grepl(";", pr_sub$Genes), , drop = FALSE]
  n_before_contam <- nrow(pr_sub)
  pr_sub <- pr_sub[pr_sub$Genes %in% clean_genes, , drop = FALSE]      # SAME contamination screen as the protein-level methods
  message("    Contamination filter: ", n_before_contam, " -> ", nrow(pr_sub), " precursor rows")
  pr_sub <- pr_sub[!duplicated(pr_sub$Precursor.Id), , drop = FALSE]
  
  ecols <- matched$clean_sample_name
  coldata <- S4Vectors::DataFrame(
    row.names = ecols,
    group = factor(ifelse(ecols %in% thy1_ids, "Thy1", "Microglia"), levels = c("Microglia", "Thy1"))
  )
  pe <- QFeatures::readQFeatures(pr_sub, quantCols = ecols, fnames = "Precursor.Id", name = "peptideRaw")
  SummarizedExperiment::rowData(pe[["peptideRaw"]])$Genes <- pr_sub$Genes[match(rownames(pe[["peptideRaw"]]), pr_sub$Precursor.Id)]
  SummarizedExperiment::colData(pe) <- coldata[colnames(pe[["peptideRaw"]]), , drop = FALSE]
  pe <- QFeatures::zeroIsNA(pe, "peptideRaw")
  pe <- QFeatures::logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
  # no normalisation (consistent with the protein-level analysis)
  pep_mat <- SummarizedExperiment::assay(pe[["peptideLog"]])
  keep_rows <- rowSums(!is.na(pep_mat)) >= 3
  pe[["peptideLog"]] <- pe[["peptideLog"]][keep_rows, ]
  message("    Precursors after >=3-observation filter: ", sum(keep_rows))
  pe <- QFeatures::aggregateFeatures(pe, i = "peptideLog", fcol = "Genes", name = "geneLevel", fun = MsCoreUtils::robustSummary)
  pe <- msqrob2::msqrob(object = pe, i = "geneLevel", formula = ~ group)
  L <- msqrob2::makeContrast(contrasts = "groupThy1 = 0", parameterNames = c("(Intercept)", "groupThy1"))
  pe_test <- msqrob2::hypothesisTest(object = pe, i = "geneLevel", contrast = L)
  rd <- SummarizedExperiment::rowData(pe_test[["geneLevel"]])
  ht <- rd[["groupThy1"]]
  tibble::tibble(gene = rownames(rd), msqrob2_logFC = ht$logFC, msqrob2_pval = ht$pval, msqrob2_padj = ht$adjPval)
}
message("  Running msqrob2 (precursor level)...")
res_msqrob2_ct <- run_msqrob2_celltype(PR_MATRIX_FILE, sample_metadata, thy1_ids, mg_ids, clean_genes)
message("    ", sum(res_msqrob2_ct$msqrob2_padj < 0.05, na.rm = TRUE), " / ", nrow(res_msqrob2_ct), " proteins at FDR q<0.05")

# ---- Sensitivity: lm(log2 intensity ~ cell type + log2 depth), same cells and proteins as the Wilcoxon test ----
message("  Running the depth-adjusted linear model (sensitivity)...")
lm_ids      <- c(thy1_ids, mg_ids)
lm_is_thy1  <- c(rep(TRUE, length(thy1_ids)), rep(FALSE, length(mg_ids)))
lm_depth    <- log2(as.numeric(n_proteins_lookup[lm_ids]))
res_lm_ct <- purrr::map_dfr(rownames(mat_de_log), function(g) {
  y <- unname(mat_de_log[g, lm_ids]); ok <- !is.na(y)
  if (sum(ok & lm_is_thy1) < 3 || sum(ok & !lm_is_thy1) < 3) return(NULL)          # same rule as the Wilcoxon test
  d <- data.frame(y = y[ok], thy1 = as.numeric(lm_is_thy1[ok]), depth = lm_depth[ok])
  cf <- summary(lm(y ~ thy1 + depth, data = d))$coefficients
  tibble::tibble(gene = g, lm_logFC = cf["thy1", "Estimate"], lm_pval = cf["thy1", "Pr(>|t|)"], lm_depth_slope = cf["depth", "Estimate"])
})
res_lm_ct$lm_padj <- p.adjust(res_lm_ct$lm_pval, method = "BH")
message("    ", sum(res_lm_ct$lm_padj < 0.05, na.rm = TRUE), " / ", nrow(res_lm_ct), " proteins at FDR q<0.05")

message("  Merging Wilcoxon / proDA / msqrob2 results...")
# Same input for all methods: the contaminant-filtered proteins of de_ct (Fig. 3b).
wilcoxon_for_merge <- de_ct %>% dplyr::select(gene, wilcoxon_logFC = logFC, wilcoxon_pval = P.Value, wilcoxon_padj = adj.P.Val)
comparison_ct <- wilcoxon_for_merge %>%
  dplyr::full_join(res_proda_ct, by = "gene") %>%
  dplyr::full_join(res_msqrob2_ct, by = "gene") %>%
  dplyr::left_join(res_lm_ct, by = "gene")
write_csv(comparison_ct, file.path(DIR_F3_SUPP, "EDFig4_method_comparison_all_proteins.csv"))

sign_agreement_ct <- function(a, b) {
  both <- !is.na(a) & !is.na(b)
  tibble::tibble(n_both_tested = sum(both), pct_sign_agree = round(100 * mean(sign(a[both]) == sign(b[both])), 1))
}
message("  Sign agreement, wilcoxon vs proda: "); print(sign_agreement_ct(comparison_ct$wilcoxon_logFC, comparison_ct$proda_logFC))
message("  Sign agreement, wilcoxon vs msqrob2: "); print(sign_agreement_ct(comparison_ct$wilcoxon_logFC, comparison_ct$msqrob2_logFC))
message("  Sign agreement, proda vs msqrob2: "); print(sign_agreement_ct(comparison_ct$proda_logFC, comparison_ct$msqrob2_logFC))

# Global agreement: effect sizes (log2FC) of ALL proteins with finite values in both methods being compared.
# Pearson r (continuous effect size) and Spearman rho (robust to a few extreme points) are both reported.
concordance_pairs <- list(
  c("wilcoxon_logFC", "proda_logFC", "Wilcoxon", "proDA"),
  c("wilcoxon_logFC", "msqrob2_logFC", "Wilcoxon", "msqrob2"),
  c("proda_logFC", "msqrob2_logFC", "proDA", "msqrob2")
)

# Axis range per method (the same in every panel in which the method appears): symmetric around 0 and wide enough for all
# proteins. If a single extreme value would squash the rest (> 5x the 99.5th percentile), the axis stops at the 99.5th
# percentile and the few points beyond it are drawn as triangles at the edge.
method_range <- function(col) {
  v <- comparison_ct[[col]]; v <- v[is.finite(v)]
  mx <- max(abs(v)); q <- as.numeric(quantile(abs(v), 0.995))
  lim <- if (mx > 5 * q) q else mx
  if (mx > 5 * q) message("  ", col, ": axis limited to +/-", round(q * 1.05, 2), " (", sum(abs(v) > q * 1.05), " extreme value(s) shown as triangles)")
  c(-lim, lim) * 1.05
}
per_method_ranges <- list(
  wilcoxon_logFC = method_range("wilcoxon_logFC"),
  proda_logFC = method_range("proda_logFC"),
  msqrob2_logFC = method_range("msqrob2_logFC")
)

ed4g_panels <- purrr::map(concordance_pairs, function(pair) {
  col1 <- pair[1]; col2 <- pair[2]; lab1 <- pair[3]; lab2 <- pair[4]
  d <- comparison_ct %>% dplyr::filter(is.finite(.data[[col1]]), is.finite(.data[[col2]]))
  r_pearson <- cor(d[[col1]], d[[col2]], method = "pearson")
  r_spearman <- cor(d[[col1]], d[[col2]], method = "spearman")
  message("  Concordance ", lab1, " vs ", lab2, ": n=", nrow(d),
          ", Pearson r=", round(r_pearson, 3), ", Spearman rho=", round(r_spearman, 3))
  xlim_trim <- per_method_ranges[[col1]]
  ylim_trim <- per_method_ranges[[col2]]
  d <- d %>% dplyr::mutate(
    off_scale = .data[[col1]] < xlim_trim[1] | .data[[col1]] > xlim_trim[2] |
      .data[[col2]] < ylim_trim[1] | .data[[col2]] > ylim_trim[2],
    x_clip = pmin(pmax(.data[[col1]], xlim_trim[1]), xlim_trim[2]),
    y_clip = pmin(pmax(.data[[col2]], ylim_trim[1]), ylim_trim[2])
  )
  n_off <- sum(d$off_scale)
  if (n_off > 0) message("    ", n_off, " point(s) off the trimmed view (shown clipped at edge, triangle)")
  ggplot(d, aes(x_clip, y_clip)) +
    geom_point_dense(data = dplyr::filter(d, !off_scale), size = 0.4, alpha = 0.3, colour = "grey40", stroke = 0) +
    geom_point(data = dplyr::filter(d, off_scale), size = 1, alpha = 0.5, colour = "grey40", shape = 17) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    coord_cartesian(xlim = xlim_trim, ylim = ylim_trim) +
    labs(x = paste0(lab1, " log2FC"), y = paste0(lab2, " log2FC")) +
    ggtitle(paste0(lab1, " vs ", lab2)) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.5, size = pt_to_geom(PT_MIN),
             label = paste0("Spearman \u03c1 = ", round(r_spearman, 2), "\nn = ", format(nrow(d), big.mark = ","))) +
    THEME_JOURNAL()
})

ed4g_concordance <- ed4g_panels[[1]] | ed4g_panels[[2]] | ed4g_panels[[3]]

# Headline-marker robustness: do the methods agree for the proteins shown in Fig. 3a/b? (lm_padj: depth-adjusted sensitivity)
marker_robustness <- comparison_ct %>%
  dplyr::filter(gene %in% all_markers) %>%
  dplyr::mutate(
    wilcoxon_sig = !is.na(wilcoxon_padj) & wilcoxon_padj < 0.05,
    proda_sig    = !is.na(proda_padj) & proda_padj < 0.05,
    msqrob2_sig  = !is.na(msqrob2_padj) & msqrob2_padj < 0.05,
    n_methods_sig = wilcoxon_sig + proda_sig + msqrob2_sig,
    direction_consistent = dplyr::case_when(
      is.na(wilcoxon_logFC) | is.na(proda_logFC) | is.na(msqrob2_logFC) ~ NA,
      TRUE ~ (sign(wilcoxon_logFC) == sign(proda_logFC)) & (sign(wilcoxon_logFC) == sign(msqrob2_logFC))
    )
  ) %>%
  dplyr::select(gene, wilcoxon_padj, proda_padj, msqrob2_padj, lm_padj, n_methods_sig, direction_consistent)
write_csv(marker_robustness, file.path(DIR_F3_SUPP, "EDFig4_headline_marker_robustness.csv"))
print(marker_robustness)

# Sensitivity of the primary test to proteome depth: Wilcoxon (no depth) vs lm with log2 depth
sig_w <- !is.na(comparison_ct$wilcoxon_padj) & comparison_ct$wilcoxon_padj < 0.05
sig_l <- !is.na(comparison_ct$lm_padj) & comparison_ct$lm_padj < 0.05
lm_pair <- comparison_ct %>% dplyr::filter(is.finite(wilcoxon_logFC), is.finite(lm_logFC))
lm_check <- list(n_tested = nrow(res_lm_ct), sig_wilcoxon = sum(sig_w), sig_lm = sum(sig_l), sig_both = sum(sig_w & sig_l),
                 dir_agree_sig_w = mean(sign(comparison_ct$wilcoxon_logFC[sig_w]) == sign(comparison_ct$lm_logFC[sig_w]), na.rm = TRUE),
                 r_pearson = cor(lm_pair$wilcoxon_logFC, lm_pair$lm_logFC),
                 r_spearman = cor(lm_pair$wilcoxon_logFC, lm_pair$lm_logFC, method = "spearman"),
                 median_depth_slope = median(res_lm_ct$lm_depth_slope))
message("  Depth sensitivity: significant (q<0.05) Wilcoxon ", lm_check$sig_wilcoxon, " | lm+depth ", lm_check$sig_lm, " | both ", lm_check$sig_both,
        " | direction agrees among Wilcoxon-significant: ", round(100 * lm_check$dir_agree_sig_w, 1), "% | log2FC Pearson r = ",
        round(lm_check$r_pearson, 3), ", Spearman rho = ", round(lm_check$r_spearman, 3))


# ---- ED 4f: one volcano per method, same marker highlighting as Fig. 3b ----------------------------------------
# All three volcanoes share the axis limits of the Wilcoxon result (primary method). Points outside the limits
# (msqrob2 has a few) are not dropped: they are drawn as triangles at the edge, with exact values in
# EDFig4f_<method>_off_scale_points.csv.
build_volcano_method <- function(data, logfc_col, padj_col, method_label, shared_xlim, shared_ylim) {
  d <- data %>%
    dplyr::mutate(
      marker_set = dplyr::case_when(
        gene %in% neuron_markers ~ "Neuronal markers",
        gene %in% microglia_markers ~ "Microglial markers",
        TRUE ~ NA_character_
      ),
      label = dplyr::if_else(!is.na(marker_set), gene, NA_character_),
      logfc_val = .data[[logfc_col]],
      padj_val = .data[[padj_col]],
      neglog10p = -log10(padj_val)
    ) %>%
    dplyr::filter(is.finite(logfc_val), is.finite(padj_val))
  n_sig <- sum(d$padj_val < 0.05, na.rm = TRUE)
  
  # Points outside the shared axis limits are clipped to the edge (triangles), not dropped.
  d <- d %>%
    dplyr::mutate(
      off_scale_x = logfc_val < shared_xlim[1] | logfc_val > shared_xlim[2],
      off_scale_y = neglog10p > shared_ylim[2],
      off_scale = off_scale_x | off_scale_y,
      logfc_clipped = pmin(pmax(logfc_val, shared_xlim[1] * 0.97), shared_xlim[2] * 0.97),
      neglog10p_clipped = pmin(neglog10p, shared_ylim[2] * 0.97),
      scale_status = dplyr::if_else(off_scale, "Off-scale (clipped)", "In range")
    )
  n_off_scale <- sum(d$off_scale)
  if (n_off_scale > 0) message("  ", method_label, ": ", n_off_scale, " point(s) off-scale (outside shared axis range), shown clipped at the edge")
  
  # Shape is mapped in aes(), so the off-scale legend entry appears only in panels that have such points.
  p <- ggplot(d, aes(logfc_clipped, neglog10p_clipped)) +
    geom_point_dense(data = dplyr::filter(d, is.na(marker_set)), aes(shape = scale_status),
                     colour = "grey80", size = 0.5, alpha = 0.4, show.legend = FALSE) +
    geom_point(data = dplyr::filter(d, !is.na(marker_set)), aes(colour = marker_set, shape = scale_status), size = 1.4) +
    ggrepel::geom_text_repel(data = dplyr::filter(d, !is.na(marker_set)), aes(label = label, colour = marker_set),
                             size = pt_to_geom(PT_MIN), max.overlaps = 30, show.legend = FALSE,
                             min.segment.length = 0, segment.size = 0.25, segment.alpha = 0.5,
                             box.padding = 0.2, point.padding = 0.1, force = 0.6) +
    scale_colour_manual(values = c("Neuronal markers" = PAL_THY1, "Microglial markers" = PAL_MICROGLIA), name = NULL) +
    scale_shape_manual(values = c("In range" = 16, "Off-scale (clipped)" = 17), name = NULL) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +   # y axis starts at 0
    coord_cartesian(xlim = shared_xlim, ylim = shared_ylim) +
    labs(x = "log2FC (neurons \u2212 microglia)", y = "\u2212log10(adj. P)") +
    THEME_JOURNAL() +
    theme(legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(0, 0, 0, 0),
          plot.title = element_text(size = PT_BASE, face = "bold", hjust = 0.5))
  p <- p + ggtitle(method_label)           # panel title = method
  if (n_off_scale == 0) p <- p + guides(shape = "none")     # no shape legend when nothing is off-scale
  # legend symbols: circle = in range, triangle = off-scale (grey, so they do not suggest a marker set)
  if (n_off_scale > 0) p <- p + guides(shape = guide_legend(override.aes = list(shape = c(16, 17), colour = "grey30", size = 1.5)))
  if (n_off_scale > 0) {
    off_pts <- d %>% dplyr::filter(off_scale) %>% dplyr::arrange(dplyr::desc(abs(logfc_val)))
    write_csv(off_pts %>% dplyr::select(gene, logfc_val, padj_val), file.path(DIR_F3_SUPP, paste0("EDFig4f_", method_label, "_off_scale_points.csv")))
  }
  p
}

wilcoxon_logfc_finite <- comparison_ct$wilcoxon_logFC[is.finite(comparison_ct$wilcoxon_logFC)]
wilcoxon_neglog10p_finite <- -log10(comparison_ct$wilcoxon_padj[is.finite(comparison_ct$wilcoxon_padj)])
wilcoxon_neglog10p_finite <- wilcoxon_neglog10p_finite[is.finite(wilcoxon_neglog10p_finite)]
x_range <- range(wilcoxon_logfc_finite); x_margin <- diff(x_range) * 0.02
shared_xlim <- c(-1, 1) * (max(abs(x_range)) + x_margin)   # symmetric around 0 so log2FC = 0 is centred in all three volcanoes
shared_ylim <- c(0, max(wilcoxon_neglog10p_finite) * 1.02)
message("  ED 4f shared axis limits: x=[", round(shared_xlim[1],2), ", ", round(shared_xlim[2],2),
        "]  y=[0, ", round(shared_ylim[2],2), "]")

ed4f_wilcoxon <- build_volcano_method(comparison_ct, "wilcoxon_logFC", "wilcoxon_padj", "Wilcoxon", shared_xlim, shared_ylim)
ed4f_proda    <- build_volcano_method(comparison_ct, "proda_logFC", "proda_padj", "proDA", shared_xlim, shared_ylim)
ed4f_msqrob2  <- build_volcano_method(comparison_ct, "msqrob2_logFC", "msqrob2_padj", "msqrob2", shared_xlim, shared_ylim)

# one shared legend at the bottom
ed4f_volcano_trio <- (ed4f_wilcoxon | ed4f_proda | ed4f_msqrob2) +
  patchwork::plot_layout(guides = "collect") & ggplot2::theme(legend.position = "bottom")

###############################################################################
# 14. ED FIG. 4: ASSEMBLY
#   row 1: a missingness | b control wells | c marker detection in controls
#   row 2: d background screen | e rank-abundance (microglia)
#   row 3: f volcano plots of the three DE methods (one letter, three plots)
#   row 4: g log2FC agreement between the methods (one letter, three plots)
# Four equal rows = 9.6 in (journal maximum is about 9.7 in / 247 mm). Panel letters are added in a vector editor.
###############################################################################
message("\n[ED Fig. 4] Assembling the extended data figure...")
ed4_row1 <- (ed4a_missingness | ed4b_pen_wells | ed4c_pen_markers) + patchwork::plot_layout(widths = c(1, 1, 1.15))
ed4_row2 <- (ed4d_pen_screen | ed4e_rank_abundance) + patchwork::plot_layout(widths = c(1, 1.45))
ed4_assembled <- (ed4_row1 / ed4_row2 / ed4f_volcano_trio / ed4g_concordance) + patchwork::plot_layout(heights = ED4_ROW_H_IN)
ggsave(file.path(DIR_F3_SUPP, "EDFig4_assembled.svg"), ed4_assembled,
       width = FIG3_W_IN, height = sum(ED4_ROW_H_IN), units = "in", device = SVG_DEVICE)
message("  ED Fig. 4 saved (", round(FIG3_W_IN, 2), " x ", sum(ED4_ROW_H_IN), " in) -- insert at 100%, do not rescale.")

###############################################################################
# 15. ED FIG. 5: HEATMAP OF THE LITERATURE CANDIDATES
# Rows    : the literature candidates (section 4) that are ALSO clearly higher in the expected cell type by the
#           detection test; neuron proteins on top, microglia proteins below, strongest detection difference
#           first. Proteins shown in panels a/b are marked with *.
# Columns : all analysis-set single cells, neurons | microglia, ordered by number of proteins detected within
#           each cell type -- fewest on the left, most on the right (no depth bar: the values are depth-adjusted).
# Colour  : depth-adjusted residual (log2). Per protein, a straight line of log2 intensity against log2(number of
#           detected protein groups) is fitted across ALL analysis-set cells (neurons + microglia together) in which
#           the protein was detected; the residual is shown (0 = level expected for a cell with that many proteins).
#           Batch is NOT in this model: all microglia come from one batch, so a batch term would absorb the
#           neuron-vs-microglia difference. Grey = not detected (no imputation). Proteins detected in fewer than
#           MIN_N_DEPTH_FIT cells get no residual and are left out (reported below).
#           Numbers next to each name: % of neurons / % of microglia that detected the protein.
###############################################################################
message("\n[Supplementary heatmap] Literature candidates, cell by cell...")
hm_cells <- de_cells %>%
  dplyr::mutate(depth = as.numeric(n_proteins_lookup[sample_id]),
                cell_type_lab = factor(dplyr::if_else(cell_type == "Thy1", LAB_NEURON, LAB_MICROGLIA),
                                       levels = c(LAB_NEURON, LAB_MICROGLIA))) %>%
  dplyr::arrange(cell_type_lab, depth) %>%
  dplyr::group_by(cell_type_lab) %>% dplyr::mutate(cell_rank = dplyr::row_number()) %>% dplyr::ungroup()

hm_raw <- protein_mat_log2[, hm_cells$sample_id, drop = FALSE]                    # plain log2, no per-cell centring

hm_prot <- lit_tests %>% dplyr::filter(clear_in_detection_test) %>% dplyr::select(gene, cell_type, detection_test_q)
hm_prot <- hm_prot %>%
  dplyr::filter(gene %in% rownames(hm_raw)) %>%
  dplyr::mutate(group = factor(dplyr::if_else(cell_type == "neuron", "Neuronal markers", "Microglial markers"),
                               levels = c("Neuronal markers", "Microglial markers"))) %>%
  dplyr::arrange(group, detection_test_q)

# depth-only residuals (use_batch = FALSE, see the header); rows without a residual are reported, not hidden
hm_rows_before <- nrow(hm_prot)
hm_resid <- depth_residuals(hm_cells$sample_id, protein_mat = protein_mat_log2[hm_prot$gene, , drop = FALSE], use_batch = FALSE)
hm_has_resid <- rowSums(!is.na(hm_resid)) > 0
hm_dropped <- hm_prot$gene[!hm_prot$gene %in% rownames(hm_resid)[hm_has_resid]]
message("  Heatmap rows: ", hm_rows_before, " candidates -> ", hm_rows_before - length(hm_dropped), " with a residual; left out (detected in < ",
        MIN_N_DEPTH_FIT, " cells): ", if (length(hm_dropped)) paste(hm_dropped, collapse = ", ") else "none")
hm_prot <- hm_prot[!hm_prot$gene %in% hm_dropped, ]

hm_is_thy1 <- hm_cells$cell_type == "Thy1"
hm_pct_n <- rowMeans(!is.na(hm_raw[hm_prot$gene, hm_is_thy1,  drop = FALSE]))
hm_pct_m <- rowMeans(!is.na(hm_raw[hm_prot$gene, !hm_is_thy1, drop = FALSE]))
hm_shown <- all_markers
hm_prot$label <- paste0(hm_prot$gene, ifelse(hm_prot$gene %in% hm_shown, " *", ""),
                        "  (", round(100 * hm_pct_n), "/", round(100 * hm_pct_m), "%)")

hm_res <- hm_resid[hm_prot$gene, , drop = FALSE]                                       # depth-adjusted residuals, NA kept
hm_long <- as.data.frame(hm_res) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "sample_id", values_to = "residual") %>%
  dplyr::left_join(hm_cells %>% dplyr::select(sample_id, cell_type_lab, cell_rank), by = "sample_id") %>%
  dplyr::left_join(hm_prot %>% dplyr::select(gene, group, label), by = "gene") %>%
  dplyr::mutate(label = factor(label, levels = rev(hm_prot$label)))                   # first protein at the top
write_csv(as.data.frame(hm_res) %>% tibble::rownames_to_column("gene"),
          file.path(DIR_F3_SUPP, "SuppFig_heatmap_literature_candidates_residuals.csv"))

# ---- Supplementary Data: one table of ALL literature candidates (the list + detection + abundance + flags) -------
# One row per candidate and cell type. Every number quoted in the text can be recounted from this table:
#   Fig. 3c            = rows with tested_in_both_readouts = TRUE (expected-way columns)
#   ED Fig. 5          = rows with in_ED_Fig5 = TRUE (clear in the detection test and with a depth-adjusted value)
#   Fig. 3a,b          = rows with shown_in_Fig3ab = TRUE
# q values come from the proteome-wide tests (section 3), the same rule everywhere.
n_thy1_cells <- sum(cell_df$cell_type == "Thy1"); n_mg_cells <- sum(cell_df$cell_type == "Microglia")
ed5_key <- paste(hm_prot$cell_type, hm_prot$gene)
supp_candidates <- lit_lists %>%
  dplyr::select(cell_type, gene, rank_by_literature_fold, literature_fold_Sharma,
                detected_Thy1_pct, detected_microglia_pct) %>%
  dplyr::left_join(cand_all %>% dplyr::select(cell_type, gene, depth_adjusted_pct_Thy1 = pct_thy1, depth_adjusted_pct_microglia = pct_microglia,
                                              detection_gap_pct_points = detection_gap, detection_expected_way, detection_q,
                                              log2FC_toward_expected_cell_type = log2fc_toward_expected, abundance_expected_way, abundance_q,
                                              detected_in_50pct_of_both = common_protein, shown_in_Fig3ab = shown_in_Fig3a),
                   by = c("cell_type", "gene")) %>%
  dplyr::left_join(lit_tests %>% dplyr::select(cell_type, gene, clear_in_detection_test, also_clear_in_intensity_test), by = c("cell_type", "gene")) %>%
  dplyr::mutate(
    n_cells_detected_Thy1 = round(det_thy1[gene] * n_thy1_cells), n_cells_detected_microglia = round(det_mg[gene] * n_mg_cells),
    n_cells_detected_total = n_cells_detected_Thy1 + n_cells_detected_microglia,
    gene_upper = toupper(gene),
    also_in_Rayaprolu2022 = dplyr::if_else(cell_type == "neuron", gene_upper %in% set_ray, NA),
    also_in_Johnson2026 = dplyr::if_else(cell_type == "neuron", gene_upper %in% set_johnson, NA),
    also_in_Lloyd2024 = dplyr::if_else(cell_type == "microglia", gene_upper %in% set_lloyd, NA),
    also_in_Neuron2026_TableS3 = dplyr::if_else(cell_type == "microglia", gene_upper %in% set_s3, NA),
    tested_in_both_readouts = !is.na(detection_expected_way) & !is.na(abundance_expected_way),
    in_ED_Fig5 = paste(cell_type, gene) %in% ed5_key) %>%
  dplyr::select(-gene_upper) %>%
  dplyr::arrange(cell_type, rank_by_literature_fold)
write_csv(supp_candidates, file.path(DIR_F3_SUPP, "SupplementaryData_literature_candidates.csv"))
message("  Supplementary table of all literature candidates: ", nrow(supp_candidates), " rows (neurons ", sum(supp_candidates$cell_type == "neuron"),
        ", microglia ", sum(supp_candidates$cell_type == "microglia"), "); in both readouts ", sum(supp_candidates$tested_in_both_readouts),
        "; in ED Fig. 5 ", sum(supp_candidates$in_ED_Fig5))
# colour range = the full range of the finite residuals (nothing is clipped), symmetric around 0
HM_LIM <- ceiling(max(abs(hm_long$residual), na.rm = TRUE) * 2) / 2
hm_finite <- hm_long$residual[is.finite(hm_long$residual)]
message("  Heatmap residuals: ", length(hm_finite), " finite values, range ", round(min(hm_finite), 2), " to ", round(max(hm_finite), 2),
        "; ", sum(abs(hm_finite) > 2), " beyond +/-2; colour scale +/-", HM_LIM)

hm_main_panel <- ggplot(hm_long, aes(cell_rank, label, fill = residual)) +
  geom_tile(width = 1.06, height = 1.06) +   # slight overlap avoids hairline seams between tiles
  facet_grid(group ~ cell_type_lab, scales = "free", space = "free") +
  scale_fill_gradientn(colours = c("#053061","#2166AC","#4393C3","#92C5DE","white","#F4A582","#D6604D","#B2182B","#67001F"),
                       values = scales::rescale(c(-1, -0.625, -0.375, -0.175, 0, 0.175, 0.375, 0.625, 1), from = c(-1, 1)),   # same steps as 01_morphology_link.R
                       limits = c(-HM_LIM, HM_LIM), oob = scales::squish, na.value = "grey85",
                       name = "Depth-adjusted\nresidual (log2)") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "single cells, ordered by number of proteins detected (fewest to most)", y = NULL) +
  THEME_OPEN_AXIS(PT_BASE) +
  theme(axis.text.y = element_text(size = PT_MIN), axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), strip.text.x = element_text(face = "bold"), strip.text.y = element_text(face = "bold"),
        panel.spacing = grid::unit(0.3, "lines"), legend.key.height = grid::unit(0.5, "cm"))

hm_full <- hm_main_panel      # no depth bar above the heatmap: the values are already depth-adjusted
# Size: journal width (FIG_W_MM); 9 in high, up to 9.45 in if needed so that the row labels (PT_MIN) do not touch
# (top strip, x label and margins take about 0.6 in).
HM_W_IN <- mm_to_in(FIG_W_MM); HM_H_MAX_IN <- 9.45
HM_H_IN <- min(HM_H_MAX_IN, max(9, 0.6 + nrow(hm_prot) * PT_MIN / 72))
hm_row_pitch_pt <- (HM_H_IN - 0.6) * 72 / nrow(hm_prot)
message("  Heatmap size: ", round(HM_W_IN, 2), " x ", round(HM_H_IN, 2), " in for ", nrow(hm_prot), " rows (", round(hm_row_pitch_pt, 1), " pt per row)")
if (hm_row_pitch_pt < PT_MIN - 1e-9) {
  message("  NOTE: ", nrow(hm_prot), " rows do not fit at ", HM_H_MAX_IN, " in with ", PT_MIN, " pt text (about ", floor((HM_H_MAX_IN - 0.6) * 72 / PT_MIN),
          " rows fit): labels will touch. Show fewer rows, or split into neuron and microglia heatmaps.")
}
ggsave(file.path(DIR_F3_SUPP, "SuppFig_heatmap_literature_candidates.svg"), hm_full, width = HM_W_IN, height = HM_H_IN, device = SVG_DEVICE)
message("  Heatmap saved: ", nrow(hm_prot), " proteins (", sum(hm_prot$group == "Neuronal markers"), " neuron, ",
        sum(hm_prot$group == "Microglial markers"), " microglia) x ", nrow(hm_cells), " cells")

###############################################################################
# 16. SUMMARY NUMBERS
# Every number needed for the Figure 3 / ED Fig. 4 legends and the results text, printed and saved to
# Figure3_summary_numbers.txt (no new analyses; the objects come from the sections above).
###############################################################################
summary_lines <- character()
say <- function(...) { l <- paste0(...); message(l); summary_lines <<- c(summary_lines, l) }
fmt_p <- function(p) formatC(p, format = "e", digits = 1)
med_iqr <- function(x) sprintf("median %.0f (IQR %.0f-%.0f)", round_half_up(median(x, na.rm = TRUE)),
                               round_half_up(quantile(x, 0.25, na.rm = TRUE)), round_half_up(quantile(x, 0.75, na.rm = TRUE)))

say("==================== FIGURE 3 / ED FIGURE 4 SUMMARY ====================")
say("Rule: contaminants removed first, then protein-count cutoffs on the filtered count. ",
    ">=", MIN_PROTEINS, " = analysis set (DE, detection, missingness, UpSet); >=200 = descriptive panels (rank-abundance).")

# ---- which cells are used where -------------------------------------------------------------------------
say("\n--- WHICH CELLS ARE USED WHERE (sample IDs compared directly; counts = contaminant-filtered protein groups) ---")
sc <- cohort_meta %>% dplyr::filter(tier == "tier_1x_single_cell", cell_type_label %in% c("Thy1", "Microglia"))
membership3 <- sc %>% dplyr::transmute(
  sample_col, cell_type = cell_type_label, batch,
  n_protein_groups_filtered = as.numeric(n_proteins_lookup[sample_col]),
  ge200 = !is.na(n_protein_groups_filtered) & n_protein_groups_filtered >= MIN_PROTEINS_VIABLE,
  ge600 = !is.na(n_protein_groups_filtered) & n_protein_groups_filtered >= MIN_PROTEINS,
  in_a_b_c_DE_detection_candidates = sample_col %in% c(thy1_ids, mg_ids),
  in_UpSet_d_e                     = sample_col %in% c(thy1_upset_cols, mg_upset_cols),
  in_missingness_ED4a              = sample_col %in% missingness_diag_cols,
  in_rank_abundance_ED4e           = sample_col %in% c(thy1_singlecell_cols, mg_singlecell_cols))
write_csv(membership3, file.path(DIR_F3, "Figure3_cell_membership.csv"))
same_cells <- function(label, a, b) say("  ", label, ": ", length(a), " vs ", length(b), ifelse(setequal(a, b), "  -> SAME CELLS", "  -> DIFFERENT CELLS, CHECK"))
say("Analysis threshold: >=", MIN_PROTEINS, " protein groups. Neurons = ", LAB_NEURON, "; microglia = ", LAB_MICROGLIA, ".")
say("a, b, c (detection, volcano, candidate summary; one analysis set): neurons ", length(thy1_ids), ", microglia ", length(mg_ids))
say("d / e (UpSet, protein counted if seen in >=", MIN_CELLS_UPSET, " cells): neurons ", length(thy1_upset_cols), ", microglia ", length(mg_upset_cols))
say("Same-cells checks:")
same_cells("neurons, a/b/c vs d (UpSet)", thy1_ids, thy1_upset_cols)
same_cells("microglia, a/b/c vs e (UpSet)", mg_ids, mg_upset_cols)
same_cells("a/b/c vs ED 4a (missingness)", c(thy1_ids, mg_ids), missingness_diag_cols)
say("Cells with >=", MIN_PROTEINS, " filtered protein groups: neurons ", sum(membership3$ge600 & membership3$cell_type == "Thy1"),
    ", microglia ", sum(membership3$ge600 & membership3$cell_type == "Microglia"))
say("Rank-abundance (ED 4e) uses the >=", MIN_PROTEINS_VIABLE, " cohorts: neurons ", sum(membership3$ge200 & membership3$cell_type == "Thy1"),
    ", microglia ", sum(membership3$ge200 & membership3$cell_type == "Microglia"))
say("Tissue-free PEN control wells: ", length(pen_sc_cols), " used for the marker and screen panels (outlier well '", PEN_EXCLUDED_WELL_PATTERN, "' left out)")
say("Sample-by-sample list written to Figure3_cell_membership.csv")

# ---- cells ------------------------------------------------------------------------------------------------
say("\n--- Cells (single-cell tier) ---")
for (ct in c("Thy1", "Microglia")) {
  d <- cohort_meta %>% dplyr::filter(tier == "tier_1x_single_cell", cell_type_label == ct)
  say(ct, ": collected ", nrow(d), " | >=", MIN_PROTEINS_VIABLE, " filtered: ", sum(d$n_proteins_filtered >= MIN_PROTEINS_VIABLE, na.rm = TRUE),
      " | >=", MIN_PROTEINS, " filtered: ", sum(d$n_proteins_filtered >= MIN_PROTEINS, na.rm = TRUE))
}
say("Analysis set used for the DE: Thy1 n=", length(thy1_ids), ", microglia n=", length(mg_ids), " (total ", length(thy1_ids) + length(mg_ids), ")")
say("Filtered protein count in the analysis set: Thy1 ", med_iqr(n_proteins_lookup[thy1_ids]), " | microglia ", med_iqr(n_proteins_lookup[mg_ids]))
say("Depth difference Thy1 vs microglia (Wilcoxon): P=", fmt_p(wilcox.test(n_proteins_lookup[thy1_ids], n_proteins_lookup[mg_ids])$p.value))

# ---- ED 4a ------------------------------------------------------------------------------------------------
say("\n--- Missingness diagnostic (ED Fig. 4a) ---")
say("Cells: ", ncol(mat_missing_diag), " (analysis set) x ", nrow(mat_missing_diag), " contaminant-filtered proteins; ",
    nrow(diag_df), " with intermediate missingness")
say("Missing fraction: Q1 = ", round(100 * quintile_summary$mean_frac_missing[1], 1), "% | Q5 = ",
    round(100 * quintile_summary$mean_frac_missing[5], 1), "%")

# ---- Fig. 3b ------------------------------------------------------------------------------------------------
say("\n--- Cell-type DE (Fig. 3b, Wilcoxon) ---")
say("log2 matrix (not normalised): ", nrow(mat_de_log), " proteins (detected in >=", 100 * MISS_THRESHOLD, "% of cells) x ", ncol(mat_de_log), " cells")
say("Tested (>=3 cells per group): ", nrow(de_ct), " | q<0.05: ", sum(de_ct$adj.P.Val < 0.05, na.rm = TRUE),
    " (", round(100 * mean(de_ct$adj.P.Val < 0.05, na.rm = TRUE), 1), "%)")
sig <- de_ct %>% dplyr::filter(adj.P.Val < 0.05)
say("  higher in Thy1: ", sum(sig$logFC > 0), " | higher in microglia: ", sum(sig$logFC < 0),
    " | with |log2FC|>1: ", sum(abs(sig$logFC) > 1), " (Thy1 ", sum(sig$logFC > 1), ", microglia ", sum(sig$logFC < -1), ")")
say("  median |log2FC| of significant proteins: ", round(median(abs(sig$logFC)), 2))
say("Contamination audit (all tested proteins): ribosomal/proteasome pattern ", n_flagged_ribo, " | blood/immune pattern ", n_flagged_blood,
    " | significant AND ribo/blood-flagged ", n_ct_flagged_sig, " (flagged, not removed)")

say("\n--- Detection test (logistic, detected ~ cell type + log2 depth) ---")
say("Tested: ", nrow(de_detect), " proteins | q<0.05: ", sum(de_detect$adj.P.Val < 0.05),
    " (Thy1-higher ", sum(de_detect$adj.P.Val < 0.05 & de_detect$z_thy1 > 0), ", microglia-higher ", sum(de_detect$adj.P.Val < 0.05 & de_detect$z_thy1 < 0), ")")

# ---- Fig. 3a ------------------------------------------------------------------------------------------------
say("\n--- Depth-adjusted marker detection (Fig. 3a) ---")
say(LAB_NEURON, " n=", length(thy1_ids), " | ", LAB_MICROGLIA, " n=", length(mg_ids), " | reference depth (pooled median) = ", round(2^ref_depth_log2), " protein groups")
say("Markers modelled: ", nrow(detect_wide), " of ", length(all_markers))
for (ms in unique(detect_wide$marker_set)) {
  d <- detect_wide[detect_wide$marker_set == ms, ]
  say("  ", ms, " (n=", nrow(d), "): predicted detection Thy1 ", med_iqr(d$pct_thy1), "% | microglia ", med_iqr(d$pct_microglia), "%")
}

say("\n--- Depth overlap check for Fig. 3a (detection is predicted at the pooled median depth) ---")
say("Reference depth (pooled median of all analysis-set cells): ", round(ref_depth_proteins), " protein groups")
for (i in seq_len(nrow(depth_overlap))) {
  d <- depth_overlap[i, ]
  say("  ", d$cell_type, " (n=", d$cells, "): protein groups min ", round_half_up(d$depth_min), " | IQR ", round_half_up(d$depth_q25), "-", round_half_up(d$depth_q75), " | median ", round_half_up(d$depth_median), " | max ", round_half_up(d$depth_max),
      "; the reference depth is at the ", d$ref_depth_percentile, "th percentile of this cell type (inside the range: ", d$ref_inside_range, "; inside the IQR: ", d$ref_inside_IQR, ")")
}
say("  Depth range covered by BOTH cell types: ", round(overlap_lo), "-", round(overlap_hi), " protein groups; share of cells inside it: ",
    paste(sprintf("%s %d%%", depth_overlap$cell_type, depth_overlap$share_cells_in_overlap_range), collapse = ", "))
# ---- Fig. 3c ------------------------------------------------------------------------------------------------
say("\n--- Literature candidates (no selection on our data; Fig. 3c) ---")
say("Candidates passing the literature criteria and detected in >= ", LIT_MIN_CELLS, " of our cells (pooled): neurons ", sum(lit_lists$cell_type == "neuron"),
    " | microglia ", sum(lit_lists$cell_type == "microglia"))
for (i in seq_len(nrow(cand_summary))) {
  r <- cand_summary[i, ]
  say(sprintf("%-9s %-10s %-42s: %3d/%3d expected way (%.0f%%); median %s %s; %d clear at q<0.05 (proteome-wide test, as ED Fig. 5)",
              r$cell_type, r$readout, r$group, r$n_expected_way, r$n_tested, r$pct_expected_way,
              if (r$readout == "detection") "difference" else "log2FC", r$median_effect, r$n_clear_q05))
}
for (i in seq_len(nrow(cand_both_summary))) {
  r <- cand_both_summary[i, ]
  say(sprintf("%-9s %-10s %-42s: %3d/%3d expected way (%.0f%%)   <- Fig. 3c (same candidates for both bars; descriptive, no test)",
              r$cell_type, r$readout, r$group, r$n_expected_way, r$n_tested, r$pct_expected_way))
}
say("(detection: median difference in percentage points, depth-adjusted; abundance: median log2 fold change toward the expected cell type)")
say("Candidates also clearly higher in our data by the detection test (heatmap set): neurons ", sum(lit_tests$cell_type == "neuron" & lit_tests$clear_in_detection_test),
    " | microglia ", sum(lit_tests$cell_type == "microglia" & lit_tests$clear_in_detection_test))
say("Supplementary heatmap: depth-adjusted residuals (straight line on log2 protein count, batch NOT in the model); ", nrow(hm_prot), " of ", hm_rows_before,
    " candidate rows shown; ", length(hm_dropped), " left out because detected in < ", MIN_N_DEPTH_FIT, " analysis-set cells (no residual)",
    if (length(hm_dropped)) paste0(": ", paste(hm_dropped, collapse = ", ")) else "")
say("  colour scale +/-", HM_LIM, " = full range of the ", length(hm_finite), " finite residuals; ", sum(abs(hm_finite) > 2), " beyond +/-2; nothing clipped")
say("Markers shown in panels a/b: ", length(neuron_markers), " neuronal (", paste(neuron_markers, collapse = ", "),
    ") | ", length(microglia_markers), " microglial (", paste(microglia_markers, collapse = ", "), ")")
say("Shown markers vs literature list:")
marker_vs_literature <- function(g, ct) {
  o <- lit_tests[lit_tests$cell_type == ct & lit_tests$clear_in_detection_test, ]
  l <- lit_lists[lit_lists$cell_type == ct, ]
  if (g %in% o$gene) "candidate + clearly higher (detection test)" else if (g %in% l$gene) "candidate, NOT clearly higher (detection test)" else "NOT a literature candidate"
}
for (g in neuron_markers)    say("  neuron    ", sprintf("%-8s", g), ": ", marker_vs_literature(g, "neuron"))
for (g in microglia_markers) say("  microglia ", sprintf("%-8s", g), ": ", marker_vs_literature(g, "microglia"))

say("\n--- CANDIDATE DENOMINATORS (Fig. 3c) ---")
say("Detection = logistic model (detected yes/no ~ cell type + log2 proteome depth), predicted at the pooled median depth; 'expected way' = the model coefficient points to the expected cell type.")
say("Abundance = Wilcoxon test of Fig. 3b on log2 intensities; 'expected way' = log2 fold change points to the expected cell type. The two readouts are separate tests.")
for (ct in c("neuron", "microglia")) {
  d <- cand_all[cand_all$cell_type == ct, ]
  has_det <- !is.na(d$detection_expected_way); has_abu <- !is.na(d$abundance_expected_way)
  say(ct, ": candidates ", nrow(d), " | detection testable ", sum(has_det), " | abundance testable ", sum(has_abu),
      " | BOTH (Fig. 3c, same candidates for both bars) ", sum(has_det & has_abu),
      " | detection only ", sum(has_det & !has_abu), " | abundance only ", sum(!has_det & has_abu))
  if (any(!has_det & has_abu)) say("   abundance-only (no detection result; model not estimable): ", paste(d$gene[!has_det & has_abu], collapse = ", "))
}

# ---- Fig. 3d, 3e, ED 4e -----------------------------------------------------------------------------------
say("\n--- UpSet panels (>=", MIN_PROTEINS, " filtered cells) and rank-abundance (>=", MIN_PROTEINS_VIABLE, " filtered cells) ---")
say("UpSet cells: microglia n=", length(mg_upset_cols), " | Thy1 (all regions) n=", length(thy1_upset_cols))
say("UpSet: protein counted as detected in this study if seen in >=", MIN_CELLS_UPSET, " cells: microglia ", mg_study$n_kept, " (", mg_study$n_ge1,
    " in >=1 cell) | neurons ", neu_study$n_kept, " (", neu_study$n_ge1, " in >=1 cell)")
say("UpSet entries -> proteins (entries sharing an accession or gene symbol are one protein; list in Figure3e_microglia_merged_entries.csv / Figure3d_neuron_merged_entries.csv): microglia ", mg_merge$n_entries, " -> ", mg_merge$n_proteins, " (", mg_merge$n_dropped, " dropped without accession/gene) | neurons ", neu_merge$n_entries, " -> ", neu_merge$n_proteins, " (", neu_merge$n_dropped, " dropped)")
say("Microglia UpSet (Fig. 3e): ", n_mg_overlap_any, "/", n_mg_study_total, " (", round(100 * n_mg_overlap_any / n_mg_study_total, 1),
    "%) of 3D-SCP proteins also in >=1 other dataset; ", n_mg_all_four, " (", round(100 * n_mg_all_four / n_mg_study_total, 1), "%) in all 4")
say("Neuron UpSet (Fig. 3d, all Thy1): ", n_neu_overlap_any, "/", n_neu_study_total, " (", round(100 * n_neu_overlap_any / n_neu_study_total, 1),
    "%) of 3D-SCP proteins also in >=1 other dataset; ", n_neu_all, " (", round(100 * n_neu_all / n_neu_study_total, 1), "%) in all ", length(neu_sets), " datasets")
say("Rank-abundance cells (>=", MIN_PROTEINS_VIABLE, "): microglia single cells n=", length(mg_singlecell_cols), " | 10-cell areas n=", length(thy1_bulk10_cols))
say("Proteins in the rank curves: microglia single ", sum(is.finite(vals_1x)), " | 10-cell areas ", sum(is.finite(vals_10x)), " | Sharma microglia ", length(vals_sharma))
say("Detection of each microglia marker (rank-abundance panel; not a DE result):")
for (g in MICROGLIA_RANK_PANEL) {
  if (!g %in% rownames(protein_mat_log2)) { say("  ", sprintf("%-8s", g), ": not in the contaminant-filtered matrix"); next }
  say(sprintf("  %-8s microglia >=%d: %d/%d | microglia >=%d: %d/%d | Thy1 >=%d: %d/%d",
              g, MIN_PROTEINS_VIABLE, sum(!is.na(protein_mat_log2[g, mg_singlecell_cols])), length(mg_singlecell_cols),
              MIN_PROTEINS, sum(!is.na(protein_mat_log2[g, mg_ids])), length(mg_ids),
              MIN_PROTEINS, sum(!is.na(protein_mat_log2[g, thy1_ids])), length(thy1_ids)))
}

# ---- ED 4b-d ----------------------------------------------------------------------------------------------
say("\n--- Tissue-free control (PEN_SC) wells (ED Fig. 4b-d) ---")
say("Control wells in the per-well panel: ", nrow(pen_wells_df), " | proteins detected per well: ", min(pen_wells_df$n_proteins_detected),
    "-", max(pen_wells_df$n_proteins_detected), " (median ", median(pen_wells_df$n_proteins_detected), ")")
say("Wells used for the marker and screen panels (outlier well left out): ", length(pen_sc_cols))
say("Markers with non-zero detection in controls: ", sum(pen_marker_detection$pct_pen > 0), " of ", nrow(pen_marker_detection),
    " (", paste(pen_marker_detection$gene[pen_marker_detection$pct_pen > 0], collapse = ", "), ")")
say("Systematic screen: ", n_pen_total_proteins, " proteins screened | ever detected in a control ", n_pen_ever_detected,
    " (", round(100 * n_pen_ever_detected / n_pen_total_proteins, 1), "%) | meet the contamination bar ", n_pen_flagged,
    " (", round(100 * n_pen_flagged / n_pen_total_proteins, 1), "% of screened; ", round(100 * n_pen_flagged / n_pen_ever_detected, 1), "% of ever-detected)")
say("Denominators: protein groups in the master matrix / screened against the controls ", n_pen_total_proteins,
    " | analysed for missingness (contaminant-filtered) ", nrow(mat_missing_diag), " | difference ", n_pen_total_proteins - nrow(mat_missing_diag),
    " (of which flagged in the controls ", n_pen_flagged, "; other exclusions, e.g. contaminant lists, ", n_pen_total_proteins - nrow(mat_missing_diag) - n_pen_flagged, ")")
say("Missingness panel, ", ncol(mat_missing_diag), " cells x ", nrow(mat_missing_diag), " protein groups: always detected ", n_always_det,
    " | never detected ", n_never_det, " | intermediate (plotted) ", nrow(diag_df))
say("Control wells: ", sum(cohort_meta$is_pen_sc), " collected | ", length(pen_sc_cols), " used for the marker and screen panels")

# ---- ED 4f, 4g ----------------------------------------------------------------------------------------------
say("\n--- DE method comparison (ED Fig. 4f, 4g) ---")
say("Significant (q<0.05): Wilcoxon=", sum(comparison_ct$wilcoxon_padj < 0.05, na.rm = TRUE),
    " | proDA=", sum(comparison_ct$proda_padj < 0.05, na.rm = TRUE), " | msqrob2=", sum(comparison_ct$msqrob2_padj < 0.05, na.rm = TRUE),
    " | of ", nrow(comparison_ct), " proteins compared")
for (pair in concordance_pairs) {
  d <- comparison_ct %>% dplyr::filter(is.finite(.data[[pair[1]]]), is.finite(.data[[pair[2]]]))
  say("Agreement ", pair[3], " vs ", pair[4], ": n=", nrow(d), ", Pearson r=", round(cor(d[[pair[1]]], d[[pair[2]]]), 3),
      ", Spearman rho=", round(cor(d[[pair[1]]], d[[pair[2]]], method = "spearman"), 3))
}
say("Headline markers (n=", nrow(marker_robustness), "): ", sum(marker_robustness$n_methods_sig == 3, na.rm = TRUE),
    " significant in all 3 methods; ", sum(marker_robustness$direction_consistent, na.rm = TRUE), "/",
    sum(!is.na(marker_robustness$direction_consistent)), " direction-consistent")
say("Depth sensitivity of the Wilcoxon test: lm(log2 intensity ~ cell type + log2 depth) on the same ", lm_check$n_tested, " proteins")
say("   significant (q<0.05): Wilcoxon ", lm_check$sig_wilcoxon, " | lm+depth ", lm_check$sig_lm, " | both ", lm_check$sig_both,
    " | direction agrees among Wilcoxon-significant proteins: ", round(100 * lm_check$dir_agree_sig_w, 1), "%")
say("   log2FC Wilcoxon vs lm+depth: Pearson r=", round(lm_check$r_pearson, 3), ", Spearman rho=", round(lm_check$r_spearman, 3),
    " | median depth slope (log2 intensity per doubling of depth): ", round(lm_check$median_depth_slope, 3))

# ---- microglial identity markers ----------------------------------------------------------------------------
say("\n--- MICROGLIAL IDENTITY MARKERS in the analysis-set profiles (undetected markers are never shown as detected) ---")
marker_status <- function(g, cols) {
  if (!g %in% rownames(protein_mat_log2)) return("not in the contaminant-filtered matrix")
  k <- sum(!is.na(protein_mat_log2[g, cols]))
  if (k == 0) paste0("not detected (0/", length(cols), ")") else paste0(k, "/", length(cols))
}
say("Collective pattern (Ctsb, Arhgdib, Pygl): detected in microglia (n=", length(mg_ids), ") vs neurons (n=", length(thy1_ids), ")")
for (g in c("Ctsb", "Arhgdib", "Pygl")) say("   ", g, ": microglia ", marker_status(g, mg_ids), " | neurons ", marker_status(g, thy1_ids))
dd <- detect_wide[as.character(detect_wide$gene) %in% c("Ctsb", "Arhgdib", "Pygl"), ]
for (i in seq_len(nrow(dd))) say("   depth-adjusted predicted detection ", as.character(dd$gene[i]), ": microglia ", round(dd$pct_microglia[i]), "% vs neurons ", round(dd$pct_thy1[i]), "%")
say("Canonical microglial markers (status in the microglia profiles):")
for (g in c("Aif1", "P2ry12", "Tmem119", "Csf1r", "Cx3cr1", "Trem2", "Siglech", "Fcrls", "Sall1", "Itgam"))
  say("   ", g, if (g == "Aif1") " (Iba1)" else "", ": microglia ", marker_status(g, mg_ids), " | neurons ", marker_status(g, thy1_ids))

say("\n====================================== END OF SUMMARY ======================================")
SW_PKGS <- c("tidyverse", "dplyr", "stringr", "readxl", "ggplot2", "ggrepel", "patchwork", "ComplexUpset", "svglite", "ggrastr", "scales", "colorspace", "data.table", "proDA", "msqrob2", "QFeatures", "SummarizedExperiment", "MsCoreUtils", "S4Vectors")
say("\nSoftware versions (all packages used by this script): ", software_versions(SW_PKGS))
write_software_versions(file.path(DIR_F3, "Figure3_software_versions.txt"), SW_PKGS)
writeLines(summary_lines, file.path(DIR_F3, "Figure3_summary_numbers.txt"))
message("\nSummary saved to ", file.path(DIR_F3, "Figure3_summary_numbers.txt"))

###############################################################################
# 17. SOURCE DATA: Figure 3, Extended Data Figures 4 and 5 (one Excel file, one sheet per panel)
###############################################################################
source_panel("Fig3a_detection", detect_wide %>% dplyr::mutate(gene = as.character(gene)),
             "Marker proteins: predicted detection (%) at the pooled median depth in neurons (thy1) and microglia, with 95% CI (logistic regression, depth-adjusted)")
source_panel("Fig3b_volcano", de_ct %>% dplyr::transmute(gene, marker_set = dplyr::case_when(gene %in% neuron_markers ~ "Neuronal markers", gene %in% microglia_markers ~ "Microglial markers", TRUE ~ NA_character_),
                                                         log2FC_neurons_minus_microglia = logFC, P_value = P.Value, adj_P_value_BH = adj.P.Val, n_neurons = n_a, n_microglia = n_b, contamination_flagged),
             "Neurons vs microglia, all tested proteins (Wilcoxon test; log2FC = difference of medians; BH-adjusted P)")
source_panel("Fig3c_candidates", cand_both_summary, "Literature candidates tested in both readouts: number and percentage higher in the expected cell type, detection and abundance")
source_panel("Fig3d_upset_neurons", upset_df_neu, "Protein groups in this study and published neuronal datasets (TRUE = present); basis of the UpSet plot")
source_panel("Fig3e_upset_microglia", upset_df_mg, "Protein groups in this study and published microglia datasets (TRUE = present); basis of the UpSet plot")
source_panel("EDFig4a_missingness", quintile_summary, "Mean fraction of missing values by abundance quintile", columns = c("abundance_quintile_label", "mean_frac_missing", "n_proteins"))
source_panel("EDFig4b_control_wells", pen_wells_df, "Protein groups per tissue-free control well (wells anonymised to numbers, ordered by protein groups)", columns = c("well_number", "n_proteins_detected"))
source_panel("EDFig4c_control_markers", pen_marker_shown, "Marker proteins detected in tissue-free control wells")
source_panel("EDFig4d_control_screen", pen_plot_df, "Per protein: detection fraction and median intensity in tissue-free control wells")
source_panel("EDFig4e_rank_abundance", rank_curves, "Rank-abundance curves of single-cell microglia, 10-cell areas and published microglia data")
source_panel("EDFig4f_g_method_comparison", comparison_ct, "Neurons vs microglia, all proteins, three differential-abundance methods (Wilcoxon, proDA, msqrob2): log2FC and adjusted P")
source_panel("EDFig5_heatmap", hm_long %>% dplyr::transmute(protein = gene, marker_set = group, label = as.character(label), sample_id, cell_type = cell_type_lab, cell_rank, depth_adjusted_residual_log2 = residual),
             "Depth-adjusted residuals (log2) of literature candidates in single cells, ordered by proteome depth within cell type; empty = not detected")
write_source_data(file.path(DIR_F3, "Source_Data_Figure3.xlsx"), "Source data for Figure 3 and Extended Data Figures 4 and 5")