###############################################################################
# 01_morphology_link.R -- 3D-SCP links whole-organ structure to single-cell proteomes
#
# Generates the quantitative panels of Figure 1 and its supplement:
#   d    example-cell heatmap: 4 traced Thy1 neurons, depth-adjusted protein levels
#   e, f Eml4 vs soma volume, Rtcb vs neurite count (exploratory candidates)
#   ED Fig. 1a  heatmap of all traced neurons
#   ED Fig. 2 gene-wise correlation screens, by-region plots, tables
# The workflow schematic and the 3D registration/tracing images are made separately.
#
# Sections:
#   1. Setup (sourced from shared/setup.R)
#   2. Data preparation: match traced neurons to their proteomes
#   3. Data preparation: depth adjustment (residuals from setup.R)
#   4. Data preparation: gene panel and example cells
#   5. Heatmaps (Fig. 1d, ED Fig. 1a)
#   6. Gene-wise correlation screens (ED Fig. 2c, Supp. Tables)
#   7. Eml4 and Rtcb scatter plots (Fig. 1e, 1f)
#   8. Supplementary plots and tables (ED Fig. 2)
#   9. Summary numbers
#
# Depth adjustment = per protein, log2 intensity ~ log2(proteins detected) + batch (see
# depth_residuals() in setup.R). A residual of +1 means "2-fold higher than expected for a
# cell of this depth". The fit is done ONCE across all Thy1 single cells; every panel selects
# cells from that one residual matrix.
###############################################################################
# Run from the project folder, or set options(scp.project_root = "/path/to/project") first.
source(file.path(getOption("scp.project_root", getwd()), "Scripts", "shared", "setup.R"))

# Nominal (uncorrected) P value, formatted for plot annotations
fmt_p_nominal <- function(p) if (p < 0.001) "nominal P < 0.001" else sprintf("nominal P = %.3f", p)

# Minimum number of cells per protein
MIN_N_SOMA_SCREEN    <- 12   # soma-volume screen (all traced cells)
MIN_N_NEURITE_SCREEN <- 6    # neurite-count screen (SoC traced cells)
N_LABELLED_PER_DIRECTION <- 5   # volcano labels: top 5 proteins per direction (by P value)
SENS_MIN_N <- 5                 # minimum complete pairs for the extra numbers in the summary (per-region, >= 600 proteins)
# All correlations in this script go through safe_spearman() (functions.R): complete finite pairs,
# minimum n, >= 2 distinct values in both variables, warnings/errors caught, invalid results rejected.

###############################################################################
# 2. DATA PREPARATION: MATCH TRACED NEURONS TO THEIR PROTEOMES
#
# Morphology file: one row per traced neuron. Proteomics_cell_id has the form
# "Somatosensory_Cortex_Cell_7"; publication_id is the gap-free number used in the figures.
# Two rows were assigned the same ID ambiguously and are marked publication_id == "excluded";
# both are dropped. Only the columns used below are kept.
###############################################################################
MORPHOLOGY_FILE <- file.path(PROJECT_ROOT, "Data", "Morphological_analysis.xlsx")
MORPHOLOGY_SHEET <- "morphological_features_Gen3_Thy"
require_file(MORPHOLOGY_FILE, "Morphology file")
if (!MORPHOLOGY_SHEET %in% readxl::excel_sheets(MORPHOLOGY_FILE)) {
  stop("Sheet '", MORPHOLOGY_SHEET, "' not found in ", MORPHOLOGY_FILE, ". Sheets: ",
       paste(readxl::excel_sheets(MORPHOLOGY_FILE), collapse = ", "), call. = FALSE)
}
# Soma volume is analysed in calibrated volume (um^3). SOMA_VOLUME_COLUMN = the name of that column in the morphology
# sheet; inside the script it is called soma_volume.
# The column is found automatically: the one column starting with "soma_volume_" that is not the pixel column.
soma_candidates <- setdiff(grep("^soma_volume_", names(read_excel(MORPHOLOGY_FILE, sheet = MORPHOLOGY_SHEET, n_max = 0)), value = TRUE), "soma_volume_px")
if (length(soma_candidates) != 1) stop("Expected exactly one calibrated soma volume column (soma_volume_*, not _px); found: ",
                                       paste(soma_candidates, collapse = ", "), call. = FALSE)
SOMA_VOLUME_COLUMN <- soma_candidates
message("Soma volume column used: ", SOMA_VOLUME_COLUMN)
SOMA_UNIT          <- "\u00b5m\u00b3"
morph_raw <- read_excel(MORPHOLOGY_FILE, sheet = MORPHOLOGY_SHEET)
require_columns(morph_raw, c("Proteomics_cell_id", "publication_id", SOMA_VOLUME_COLUMN, "neurite_count"), "Morphology sheet")
if (!is.numeric(morph_raw[[SOMA_VOLUME_COLUMN]]) || !is.numeric(morph_raw$neurite_count)) {
  stop(SOMA_VOLUME_COLUMN, " and neurite_count must be numeric in the morphology sheet.", call. = FALSE)
}
morph_raw$soma_volume <- morph_raw[[SOMA_VOLUME_COLUMN]]
morph_clean <- morph_raw %>%
  dplyr::filter(publication_id != "excluded") %>%
  dplyr::select(Proteomics_cell_id, publication_id, soma_volume, neurite_count)
if (anyDuplicated(morph_clean$Proteomics_cell_id)) {
  stop("Duplicated Proteomics_cell_id in the morphology file: ",
       paste(unique(morph_clean$Proteomics_cell_id[duplicated(morph_clean$Proteomics_cell_id)]), collapse = ", "), call. = FALSE)
}
message("Morphology rows: ", nrow(morph_raw), " total, ", nrow(morph_clean), " after dropping the ambiguous pair")

# Build the same ID ("Region_Cell_N") for the proteomics samples of the traced batch (24Apr2026,
# Thy1) and join. cohort_meta (setup.R) already excludes duplicate neurons and invalid wells.
sub_meta <- cohort_meta %>%
  dplyr::filter(cell_type_label == "Thy1", batch == "24Apr2026", !is.na(cell_num)) %>%
  dplyr::mutate(constructed_id = paste0(gsub(" ", "_", region), "_Cell_", as.integer(cell_num)))
if (anyDuplicated(sub_meta$constructed_id)) {
  stop("Duplicated constructed IDs in the proteomics metadata (the join would duplicate rows): ",
       paste(unique(sub_meta$constructed_id[duplicated(sub_meta$constructed_id)]), collapse = ", "), call. = FALSE)
}

no_proteome <- setdiff(morph_clean$Proteomics_cell_id, sub_meta$constructed_id)
if (length(no_proteome) > 0) message("Traced neurons without a matching proteome: ", paste(no_proteome, collapse = ", "))

# merged_morph: one row per traced neuron with a proteome (>= 200 protein groups), morphology + cell metadata
merged_morph <- morph_clean %>%
  dplyr::inner_join(sub_meta, by = c("Proteomics_cell_id" = "constructed_id")) %>%
  dplyr::filter(passes_qc_tier)
if (nrow(merged_morph) == 0) stop("No traced neuron could be matched to a QC-passing proteome.", call. = FALSE)
n_not_single_cell <- sum(!merged_morph$sample_col %in% thy1_singlecell_cols)
merged_morph <- merged_morph %>% dplyr::filter(sample_col %in% thy1_singlecell_cols)
message("Traced neurons matched to a QC-passing Thy1 single-cell proteome: ", nrow(merged_morph),
        if (n_not_single_cell > 0) paste0(" (", n_not_single_cell, " dropped: not in the single-cell tier)") else "")

# Supplementary Data: every traced neuron (also those without a usable proteome) with its morphology, so that the
# soma-volume and neurite-count screens and Extended Data Fig. 1 and 2a can be recounted from one table.
morph_table <- morph_clean %>%
  dplyr::left_join(sub_meta %>% dplyr::select(constructed_id, sample_col, region, n_proteins_filtered, passes_qc_tier),
                   by = c("Proteomics_cell_id" = "constructed_id")) %>%
  dplyr::transmute(Proteomics_cell_id, publication_id, region,
                   soma_volume_um3 = soma_volume, neurite_count,
                   n_protein_groups_filtered = n_proteins_filtered,
                   proteome_at_least_200_protein_groups = !is.na(passes_qc_tier) & passes_qc_tier,
                   in_correlation_screens = !is.na(sample_col) & sample_col %in% merged_morph$sample_col)
write_csv(morph_table, file.path(DIR_F1_SUPP, "SuppTable_morphology_per_traced_neuron.csv"))
message("Morphology table: ", nrow(morph_table), " traced neurons, ", sum(morph_table$in_correlation_screens), " in the correlation screens")

###############################################################################
# 3. DATA PREPARATION: DEPTH ADJUSTMENT
# One fit across all Thy1 single cells (all regions, all batches), see setup.R.
###############################################################################
resid_thy1 <- depth_residuals(thy1_singlecell_cols)
n_usable_genes <- sum(rowSums(!is.na(resid_thy1)) > 0)
message("Depth fit: ", ncol(resid_thy1), " Thy1 single cells; ", n_usable_genes, " / ", nrow(resid_thy1),
        " proteins have a residual (detected in >= ", MIN_N_DEPTH_FIT, " cells)")
stopifnot(all(merged_morph$sample_col %in% colnames(resid_thy1)),
          all(merged_morph$n_proteins_filtered == n_proteins_lookup[merged_morph$sample_col]))

###############################################################################
# 4. DATA PREPARATION: GENE PANEL AND EXAMPLE CELLS
###############################################################################
# Neuron-function protein panel, grouped by function (one continuous heatmap; the groups are shown
# by gene order, a label strip and separator lines). Multi-subunit complexes are represented by
# one member (e.g. Atp6v0a1 for the vacuolar ATPase, Tuba1b for tubulins). None of the proteins is
# flagged as a contaminant or as non-brain tissue protein.
gene_groups <- list(
  "Synaptic vesicle/release" = c("Syn1","Syp","Syt1","Snap25","Snap91","Vamp2","Cplx2","Dnm1","Dnm3","Nsf","Rab8b"),
  "Postsynaptic scaffold"    = c("Homer1","Shank1","Shank2","Shank3","Dlg2","Gria2","Grin1"),
  "Cell adhesion"            = c("Ncam1","L1cam","Nlgn3","Nrxn3","Gpm6a"),
  "Cytoskeleton"             = c("Tuba1b","Map1b","Map2","Stmn1","Sptbn1","Septin7","Septin11","Dpysl2","Ctnnd2","Nefm","Pfn2","Rac1"),
  "Calcium/ion transport"    = c("Camk2a","Camk2b","Atp2b2","Atp2a2","Atp1b1","Atp6v0a1","Ppp3ca","Slc25a12","Hpca"),
  "Vesicle trafficking"      = c("Ap2b1","Cltc","Gdi1","Gnai1","Gna13","Gnb2"),
  "General/identity"         = c("Thy1","Uchl1","Eno2","Nptn","Gap43","Ank3","Mapk1")
)
stopifnot(!anyDuplicated(unlist(gene_groups)))
gene_panel_f1 <- unlist(gene_groups, use.names = FALSE)
gene_to_group <- setNames(rep(names(gene_groups), lengths(gene_groups)), gene_panel_f1)

# Example cells: the best-covered members of two visually chosen groups (Cells 1 and 10: different
# looking; Cells 7 and 14: similar looking). Chosen from the morphology groupings, not by searching for
# the largest differences. Labels use publication_id (Cell 14 is "Cell 13" in the figures).
final_cell_ids <- c("Somatosensory_Cortex_Cell_1", "Somatosensory_Cortex_Cell_7",
                    "Somatosensory_Cortex_Cell_10", "Somatosensory_Cortex_Cell_14")
final_labels   <- c("Cell 1", "Cell 7", "Cell 10", "Cell 13")
final_sample_cols <- merged_morph$sample_col[match(final_cell_ids, merged_morph$Proteomics_cell_id)]
if (anyNA(final_sample_cols)) stop("Example cell(s) not found in merged_morph: ", paste(final_cell_ids[is.na(final_sample_cols)], collapse = ", "))
names(final_sample_cols) <- final_labels

# Panel proteins with a residual in at least one example cell; ordered by group, then alphabetically
usable_genes_f1 <- gene_panel_f1[gene_panel_f1 %in% rownames(resid_thy1)]
usable_genes_f1 <- usable_genes_f1[rowSums(!is.na(resid_thy1[usable_genes_f1, final_sample_cols, drop = FALSE])) > 0]
gene_order_f1 <- unlist(lapply(gene_groups, function(g) sort(g[g %in% usable_genes_f1])), use.names = FALSE)
if (length(gene_order_f1) == 0) stop("None of the panel proteins has a residual in the example cells.", call. = FALSE)
not_in_panel <- setdiff(gene_panel_f1, gene_order_f1)
if (length(not_in_panel) > 0) message("Panel proteins without a residual in any example cell (not shown): ", paste(not_in_panel, collapse = ", "))
message("Panel: ", length(gene_order_f1), " of ", length(gene_panel_f1), " proteins have a residual in at least one example cell")

# Position of the groups along the x axis (label strip and separator lines)
# Positions are taken from the gene order that is actually drawn (gene_order_f1), so a group that lost
# all its proteins in the filtering simply does not exist here: no empty group, no repeated or stray
# separator line. Column k of the heatmap sits at x = k.
shown_runs    <- rle(unname(gene_to_group[gene_order_f1]))
stopifnot(!anyDuplicated(shown_runs$values))             # every group must be one continuous block of columns
group_sizes   <- shown_runs$lengths                      # columns per shown group
group_bounds  <- cumsum(group_sizes)                     # last column of each shown group
group_starts  <- group_bounds - group_sizes              # number of columns before each shown group
n_genes_f1    <- length(gene_order_f1)
group_label_df <- data.frame(group = shown_runs$values, n_genes = group_sizes,
                             centre = group_starts + (group_sizes + 1) / 2)   # middle of columns start+1 ... bound

# wrap_group_label(): breaks a group name into lines so that no line is wider than its own group of
# columns (a character of bold Arial is about 0.62 em wide). Breaks are allowed at spaces and after "/".
# panel_width_in = approximate width of the heatmap area in inches (a little too small is safe).
wrap_group_label <- function(label, n_cols, panel_width_in, pt = PT_MIN) {
  col_in    <- panel_width_in / n_genes_f1
  max_chars <- max(4, floor(0.95 * n_cols * col_in * 72 / (0.62 * pt)))
  lines <- gsub("/ ", "/", strwrap(gsub("/", "/ ", label), width = max_chars))
  if (max(nchar(lines)) > max_chars) {
    message("  NOTE: heatmap group label '", label, "' may be wider than its ", n_cols, " columns; shorten the name or widen the figure")
  }
  paste(lines, collapse = "\n")
}

###############################################################################
# 5. HEATMAPS (Fig. 1d, ED Fig. 1a)
###############################################################################
# Matrix (proteins x cells) -> long table
resid_to_long <- function(mat) {
  as.data.frame(mat) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "sample_col", values_to = "residual")
}

# One heatmap function for all heatmaps of this figure. long_df needs the columns gene (factor, in
# panel order), row (factor; the FIRST level is drawn at the BOTTOM) and residual; with
# facet_regions = TRUE also region_abbrev. Uses group_label_df / group_bounds from above.
#   lim = the colour scale runs from -lim to +lim (RESID_LIM, set below from the data); the function stops if a
#         displayed value would fall outside it, so nothing is clipped silently.
#   panel_width_in = approximate width of the heatmap area, only used to break the group names into lines.
#   strip_height_in = height of the group-name strip in inches.
plot_resid_heatmap <- function(long_df, lim, bg = "white", missing_colour = "black", border_colour = "black",
                               border_lw = 0.3, facet_regions = FALSE, y_text_size = PT_BASE,
                               y_text_bold = FALSE, strip_height_in = 0.3, panel_width_in = 4.9) {
  dark <- identical(bg, "black")
  fg   <- if (dark) "white" else "black"
  mid  <- if (dark) "grey85" else "white"
  stopifnot(all(abs(long_df$residual[is.finite(long_df$residual)]) <= lim + 1e-9))   # nothing is clipped
  
  label_df <- group_label_df
  label_df$label <- mapply(wrap_group_label, label_df$group, label_df$n_genes, MoreArgs = list(panel_width_in = panel_width_in))
  # The heatmap x axis is discrete (columns at 1..n, default padding 0.6 on both sides); the strip uses the
  # same range, so every group name sits exactly above the middle of its own columns.
  p_strip <- ggplot(label_df, aes(x = centre, y = 1, label = label)) +
    geom_text(colour = fg, fontface = "bold", size = pt_to_geom(PT_MIN), family = "Arial", lineheight = 0.85) +
    scale_x_continuous(limits = c(0.4, n_genes_f1 + 0.6), expand = c(0, 0)) +
    theme_void()
  
  # colour steps are fixed fractions of the range, so they move with `lim`
  p_heat <- ggplot(long_df, aes(x = gene, y = row, fill = residual)) +
    geom_tile(colour = border_colour, linewidth = border_lw) +
    geom_vline(xintercept = group_bounds[-length(group_bounds)] + 0.5, colour = fg, linewidth = 0.4) +
    scale_fill_gradientn(colours = c("#053061", "#2166AC", "#4393C3", "#92C5DE", mid, "#F4A582", "#D6604D", "#B2182B", "#67001F"),
                         values = scales::rescale(c(-1, -0.625, -0.375, -0.175, 0, 0.175, 0.375, 0.625, 1), from = c(-1, 1)),
                         limits = c(-lim, lim), oob = scales::squish,
                         na.value = missing_colour,
                         name = "Depth-adjusted\nresidual (log2)") +
    labs(x = NULL, y = NULL) +
    THEME_JOURNAL() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = PT_MIN, colour = fg),
          axis.text.y = element_text(size = y_text_size, face = if (y_text_bold) "bold" else "plain", colour = fg),
          legend.key.height = unit(4, "mm"), legend.key.width = unit(2.5, "mm"),
          panel.grid = element_blank())
  
  if (facet_regions) {
    p_heat <- p_heat +
      facet_grid(region_abbrev ~ ., scales = "free_y", space = "free_y", switch = "y") +
      theme(strip.placement = "outside", strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 0, face = "bold", size = PT_BASE),
            panel.spacing.y = unit(1.5, "mm"))
  }
  if (dark) {
    p_heat <- p_heat + theme(plot.background = element_rect(fill = "black", colour = NA),
                             panel.background = element_rect(fill = "black", colour = NA),
                             legend.background = element_rect(fill = "black", colour = NA),
                             legend.text = element_text(colour = "white"),
                             legend.title = element_text(colour = "white"))
    p_strip <- p_strip + theme(plot.background = element_rect(fill = "black", colour = NA))
  }
  bg_theme <- theme(plot.background = element_rect(fill = if (dark) "black" else "white", colour = NA),
                    plot.margin = if (facet_regions) margin(5.5, 5.5, 5.5, 5.5) else margin(0, 0, 0, 0))
  p_strip / p_heat +
    patchwork::plot_layout(heights = grid::unit(c(strip_height_in, 1), c("in", "null"))) +
    patchwork::plot_annotation(theme = bg_theme)
}

# ---- Colour range: the full range of the residuals shown (panel d cells are a subset of the ED Fig. 1a cells) ----
# One shared scale for both heatmaps, so the same colour means the same residual in Fig. 1d and ED Fig. 1a.
# The extremes are listed first and nothing is removed or clipped; missing values (protein not detected)
# keep their own colour and are NOT part of the colour bar (explain them in the figure legend).
resid_shown <- resid_to_long(resid_thy1[gene_order_f1, merged_morph$sample_col, drop = FALSE]) %>%
  dplyr::filter(is.finite(residual))
RESID_LIM <- ceiling(max(abs(resid_shown$residual)) * 2) / 2     # next 0.5 above the largest |residual|
message("Heatmap residuals shown: ", nrow(resid_shown), " finite values (range ", round(min(resid_shown$residual), 2), " to ",
        round(max(resid_shown$residual), 2), "); ", sum(abs(resid_shown$residual) > 4), " beyond +/-4 -> colour scale +/-", RESID_LIM)
resid_extreme <- resid_shown %>% dplyr::filter(abs(residual) > 4) %>% dplyr::arrange(dplyr::desc(abs(residual)))
if (nrow(resid_extreme) > 0) {
  print(as.data.frame(head(resid_extreme, 15)))
  write_csv(resid_extreme, file.path(DIR_F1_SUPP, "SuppTable_heatmap_residuals_beyond_4.csv"))
}

# ---- Panel d: 4 example cells, white and black background ----
plot_df <- resid_to_long(resid_thy1[gene_order_f1, final_sample_cols, drop = FALSE]) %>%
  dplyr::mutate(row  = factor(names(final_sample_cols)[match(sample_col, final_sample_cols)], levels = rev(final_labels)),
                gene = factor(gene, levels = gene_order_f1))

# Missing values (protein not detected) are black, blending into the background.
p1_white <- plot_resid_heatmap(plot_df, lim = RESID_LIM, bg = "white", y_text_bold = TRUE, panel_width_in = 4.9)
p1_black <- plot_resid_heatmap(plot_df, lim = RESID_LIM, bg = "black", y_text_bold = TRUE, panel_width_in = 4.9)
# Panel d: 6.4 x 2.0 inches
ggsave(file.path(DIR_F1, "Figure1_example_cells_heatmap_white.svg"), p1_white,
       width = 6.4, height = 2.0, device = SVG_DEVICE)
ggsave(file.path(DIR_F1, "Figure1_example_cells_heatmap_black.svg"), p1_black,
       width = 6.4, height = 2.0, device = SVG_DEVICE, bg = "black")
message("Panel d saved (white and black background) to ", DIR_F1)

# ---- ED Fig. 1a: all traced neurons, grouped by region ----
# EVERY traced neuron (all rows of the morphology table after dropping the ambiguous pair), also those without a
# proteome or with < 200 protein groups. Those neurons have no depth-adjusted values (they are not part of the
# depth fit or any analysis), so their row is completely grey; their label ends with an asterisk (*).
# Missing values are grey. Nothing is imputed.
tr_all <- morph_clean %>%
  dplyr::left_join(dplyr::select(sub_meta, constructed_id, sample_col, region), by = c("Proteomics_cell_id" = "constructed_id")) %>%
  dplyr::mutate(region = ifelse(is.na(region), gsub("_", " ", sub("_Cell_[0-9]+$", "", Proteomics_cell_id)), region),
                pub = as.integer(publication_id),
                region_abbrev = factor(unname(REGION_ABBREV[region]), levels = unname(REGION_ABBREV[REGIONS_THY1])),
                in_analysis = !is.na(sample_col) & sample_col %in% colnames(resid_thy1)) %>%
  dplyr::arrange(region_abbrev, pub)
if (anyNA(tr_all$region_abbrev)) stop("ED Fig. 1a: a traced neuron has a region that is not in REGION_ABBREV: ",
                                      paste(unique(tr_all$region[is.na(tr_all$region_abbrev)]), collapse = ", "), call. = FALSE)
if (sum(tr_all$in_analysis) != nrow(merged_morph)) stop("ED Fig. 1a: neurons with values (", sum(tr_all$in_analysis),
                                                        ") differ from merged_morph (", nrow(merged_morph), ").", call. = FALSE)
# publication_id is a gap-free 1..n sequence; if it is not unique across regions, prefix the region
tr_all$label <- if (!any(duplicated(tr_all$pub))) paste0("Cell ", tr_all$pub) else paste0(tr_all$region_abbrev, " ", tr_all$pub)
tr_all$label <- paste0(tr_all$label, ifelse(tr_all$in_analysis, "", " *"))

mat_ed <- matrix(NA_real_, nrow = length(gene_order_f1), ncol = nrow(tr_all), dimnames = list(gene_order_f1, tr_all$label))
mat_ed[, tr_all$label[tr_all$in_analysis]] <- resid_thy1[gene_order_f1, tr_all$sample_col[tr_all$in_analysis], drop = FALSE]
plot_df_ed <- resid_to_long(mat_ed) %>%
  dplyr::mutate(idx = match(sample_col, tr_all$label),
                row = factor(sample_col, levels = rev(tr_all$label)),     # first cell at the top
                region_abbrev = tr_all$region_abbrev[idx],
                gene = factor(gene, levels = gene_order_f1))

p_ed1a <- plot_resid_heatmap(plot_df_ed, lim = RESID_LIM, bg = "white", missing_colour = "grey80", border_colour = "grey92",
                             border_lw = 0.1, facet_regions = TRUE, y_text_size = PT_MIN, strip_height_in = 0.3, panel_width_in = 5.4)
ed1a_height_mm <- max(80, 2.4 * nrow(tr_all) + 40)   # ~2.4 mm per cell row at 6 pt
ggsave(file.path(DIR_F1_SUPP, "EDFig1a_all_traced_neurons_heatmap.svg"), p_ed1a,
       width = mm_to_in(FIG_W_MM), height = mm_to_in(ed1a_height_mm), device = SVG_DEVICE)
message("ED Fig. 1a saved: ", nrow(tr_all), " traced neurons x ", n_genes_f1, " proteins (", sum(tr_all$in_analysis),
        " with values, ", sum(!tr_all$in_analysis), " without: no proteome or < ", MIN_PROTEINS_VIABLE, " protein groups, marked *); per region: ",
        paste(names(table(tr_all$region_abbrev)), table(tr_all$region_abbrev), collapse = "; "))

###############################################################################
# 6. GENE-WISE CORRELATION SCREENS (ED Fig. 2c, Supp. Tables)
# For every protein: Spearman correlation between its depth-adjusted residual and a morphology
# measure across traced neurons. Per protein only the cells in which it was detected (a residual exists)
# are used; nothing is imputed. Every correlation goes through safe_spearman() (functions.R).
# Proteins that cannot be tested are NOT dropped silently: they are written, with the reason, to a
# "skipped proteins" table. P values are nominal; q = Benjamini-Hochberg over the TESTED proteins of the
# screen. Output tables: gene, spearman_r, p_nominal, q_BH, n (complete pairs), ties_present, warning.
###############################################################################
run_full_gene_correlation <- function(cell_meta, morph_col, min_n, resid = resid_thy1) {
  cols  <- intersect(cell_meta$sample_col, colnames(resid))
  morph <- setNames(cell_meta[[morph_col]], cell_meta$sample_col)[cols]
  if (length(cols) < min_n) stop(morph_col, ": only ", length(cols), " cells; fewer than min_n = ", min_n, call. = FALSE)
  if (length(unique(na.omit(morph))) < 2) stop(morph_col, " is constant in the cells given; a correlation is not defined.", call. = FALSE)
  per_protein <- purrr::map_dfr(rownames(resid), function(g) {
    dplyr::bind_cols(tibble::tibble(gene = g),
                     safe_spearman(resid[g, cols], morph, min_n, what = c("protein residual", morph_col)))
  })
  skipped <- per_protein %>% dplyr::filter(status == "skipped") %>%
    dplyr::select(gene, n_complete_pairs = n, reason)
  results <- per_protein %>% dplyr::filter(status == "tested") %>%
    dplyr::select(gene, spearman_r, p_nominal, n, ties_present, warning)
  if (nrow(results) == 0) stop(morph_col, ": no protein could be tested (minimum ", min_n, " complete pairs).", call. = FALSE)
  results$q_BH <- p.adjust(results$p_nominal, method = "BH")
  results <- results[, c("gene", "spearman_r", "p_nominal", "q_BH", "n", "ties_present", "warning")]
  message("  ", morph_col, ": ", nrow(results), " proteins tested, ", nrow(skipped), " skipped (", length(cols), " cells, minimum ", min_n, " pairs). Reasons: ",
          paste(names(table(skipped$reason)), table(skipped$reason), sep = " = ", collapse = "; "),
          "; P approximate because of ties in ", sum(results$ties_present), " tested proteins")
  list(results = dplyr::arrange(results, p_nominal), skipped = skipped)
}

traced_soc <- merged_morph %>% dplyr::filter(region == "Somatosensory Cortex")   # neurite count only varies enough in SoC (see section 7)
soma_screen    <- run_full_gene_correlation(merged_morph, "soma_volume", MIN_N_SOMA_SCREEN)
neurite_screen <- run_full_gene_correlation(traced_soc, "neurite_count", MIN_N_NEURITE_SCREEN)
soma_corr_table    <- soma_screen$results
neurite_corr_table <- neurite_screen$results
write_csv(soma_corr_table, file.path(DIR_F1_SUPP, "SuppTable_soma_volume_correlation_all_genes.csv"))
write_csv(neurite_corr_table, file.path(DIR_F1_SUPP, "SuppTable_neurite_count_correlation_all_genes.csv"))
write_csv(soma_screen$skipped, file.path(DIR_F1_SUPP, "SuppTable_soma_volume_correlation_skipped_proteins.csv"))
write_csv(neurite_screen$skipped, file.path(DIR_F1_SUPP, "SuppTable_neurite_count_correlation_skipped_proteins.csv"))
message("Screens: soma volume ", nrow(soma_corr_table), " proteins tested (", sum(soma_corr_table$p_nominal < 0.05), " nominal P<0.05, ",
        sum(soma_corr_table$q_BH < 0.1), " BH q<0.1); neurite count ", nrow(neurite_corr_table), " proteins tested (",
        sum(neurite_corr_table$p_nominal < 0.05), " nominal P<0.05, ", sum(neurite_corr_table$q_BH < 0.1), " BH q<0.1)")

# Volcano-style plot: the top N proteins per direction (by P value) are labelled
plot_corr_volcano <- function(corr_table, top_n, xlab) {
  d <- corr_table
  hl <- dplyr::bind_rows(d %>% dplyr::filter(spearman_r > 0) %>% dplyr::arrange(p_nominal) %>% head(top_n),
                         d %>% dplyr::filter(spearman_r < 0) %>% dplyr::arrange(p_nominal) %>% head(top_n))
  corr_table$is_hl <- corr_table$gene %in% hl$gene
  ggplot(corr_table, aes(spearman_r, -log10(p_nominal))) +
    geom_point(data = subset(corr_table, !is_hl), colour = "grey75", size = 0.6, alpha = 0.5) +
    geom_point(data = subset(corr_table, is_hl), colour = "#B2182B", size = 1.6) +
    ggrepel::geom_text_repel(data = subset(corr_table, is_hl), aes(label = gene), colour = "#B2182B",
                             fontface = "italic", size = pt_to_geom(PT_MIN), max.overlaps = 20, min.segment.length = 0,
                             segment.size = 0.3, box.padding = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    labs(x = xlab, y = "-log10(nominal P)") +
    THEME_JOURNAL()
}
p_soma_volcano    <- plot_corr_volcano(soma_corr_table, N_LABELLED_PER_DIRECTION, "Spearman r (vs. soma volume)")
p_neurite_volcano <- plot_corr_volcano(neurite_corr_table, N_LABELLED_PER_DIRECTION, "Spearman r (vs. neurite count)")

###############################################################################
# 7. EML4 AND RTCB SCATTER PLOTS (Fig. 1e, 1f)
# Two exploratory candidates. The statistics (rho, nominal P, BH q, n) are TAKEN FROM THE SCREEN TABLES
# of section 6, not recomputed, so the panels and the all-protein screens cannot disagree. They are
# reported in the legend, not drawn on the panels. Whether a candidate passes BH is checked and printed
# below (it is not assumed). The grey line is an ordinary least-squares line; the statistics are Spearman.
#   Eml4 vs soma volume: all traced neurons (3 regions). Rtcb vs neurite count: SoC only, because
#   neurite counts in mPFC and Hippo are close to zero (ceiling effect) and SoC is the one region
#   with usable variance.
###############################################################################
if (!all(c("Eml4", "Rtcb") %in% rownames(resid_thy1))) stop("Eml4 and/or Rtcb not in the residual matrix.", call. = FALSE)
eml4_df <- merged_morph %>%
  dplyr::mutate(eml4_residual = resid_thy1["Eml4", sample_col]) %>%
  dplyr::filter(!is.na(eml4_residual), !is.na(soma_volume))
rtcb_df <- traced_soc %>%
  dplyr::mutate(rtcb_residual = resid_thy1["Rtcb", sample_col]) %>%
  dplyr::filter(!is.na(rtcb_residual), !is.na(neurite_count))

# Statistics of one candidate = its row in the screen table; stops if the protein was not tested there
# or if the plotted cells differ from the cells the screen used.
screen_row_for <- function(tbl, g, n_plotted, what) {
  r <- tbl[tbl$gene == g, ]
  if (nrow(r) != 1) stop(g, " is not in the ", what, " screen table (not testable; see the skipped-proteins table).", call. = FALSE)
  if (r$n != n_plotted) stop(g, ": the scatterplot shows ", n_plotted, " cells but the ", what, " screen used ", r$n, ".", call. = FALSE)
  r
}
scr_eml4 <- screen_row_for(soma_corr_table, "Eml4", nrow(eml4_df), "soma-volume")
scr_rtcb <- screen_row_for(neurite_corr_table, "Rtcb", nrow(rtcb_df), "neurite-count")
cand_label <- function(r) sprintf("Spearman rho = %.2f, %s, BH q = %.2f, n = %d%s", r$spearman_r, fmt_p_nominal(r$p_nominal),
                                  r$q_BH, r$n, if (isTRUE(r$ties_present)) " (ties: P approximate)" else "")
eml4_label <- cand_label(scr_eml4)
rtcb_label <- cand_label(scr_rtcb)
message("Panel e (legend text): ", eml4_label)
message("Panel f (legend text): ", rtcb_label)
message("BH check: Eml4 q = ", signif(scr_eml4$q_BH, 3), ", Rtcb q = ", signif(scr_rtcb$q_BH, 3),
        if (scr_eml4$q_BH < 0.1 || scr_rtcb$q_BH < 0.1) "  -> at least one candidate PASSES q < 0.1: update the 'exploratory' wording"
        else "  -> neither passes q < 0.1: both stay exploratory candidates")

DARK_SCATTER_THEME <- theme(
  plot.background = element_rect(fill = "black", colour = NA),
  panel.background = element_rect(fill = "black", colour = NA),
  axis.text = element_text(colour = "white"),
  axis.title = element_text(colour = "white"),
  axis.line = element_line(colour = "white"),
  axis.ticks = element_line(colour = "white"),
  legend.background = element_rect(fill = "black", colour = NA),
  legend.text = element_text(colour = "white"),
  legend.title = element_text(colour = "white")
)

p_eml4 <- ggplot(eml4_df, aes(x = soma_volume, y = eml4_residual)) +
  geom_point(aes(colour = region), size = 1.2, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey70", linewidth = 0.5) +
  scale_colour_manual(values = PAL_REGION_THY1, labels = REGION_ABBREV, name = NULL) +
  labs(x = paste0("Soma volume (", SOMA_UNIT, ")"), y = "Eml4 depth-adjusted residual") +
  THEME_JOURNAL() + DARK_SCATTER_THEME +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(2.5, "mm"), legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-4, 0, 0, 0))
ggsave(file.path(DIR_F1, "Figure1_Eml4_vs_soma_volume.svg"), p_eml4,   # panel e: 2.7 x 1.9 in
       width = 2.7, height = 1.9, device = SVG_DEVICE, bg = "black")

p_rtcb <- ggplot(rtcb_df, aes(x = neurite_count, y = rtcb_residual)) +
  geom_point(colour = unname(PAL_REGION_THY1["Somatosensory Cortex"]), size = 1.2, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey70", linewidth = 0.5) +
  labs(x = "Neurite count", y = "Rtcb depth-adjusted residual") +
  THEME_JOURNAL() + DARK_SCATTER_THEME
ggsave(file.path(DIR_F1, "Figure1_Rtcb_vs_neurite_count.svg"), p_rtcb,   # panel f: 2.7 x 1.9 in
       width = 2.7, height = 1.9, device = SVG_DEVICE, bg = "black")
message("Panels e and f saved to ", DIR_F1)

###############################################################################
# 8. SUPPLEMENTARY PLOTS AND TABLES (ED Fig. 2)
###############################################################################
# Eml4 within each region: the trend is not a single-region effect
p_eml4_by_region <- ggplot(eml4_df, aes(x = soma_volume, y = eml4_residual, colour = region)) +
  geom_point(size = 1.2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  scale_colour_manual(values = PAL_REGION_THY1, labels = REGION_ABBREV, name = NULL) +
  labs(x = paste0("Soma volume (", SOMA_UNIT, ")"), y = "Eml4 depth-adjusted residual") +
  THEME_JOURNAL() +
  theme(legend.position = "right", legend.key.size = unit(2.5, "mm"))

# One row: soma-volume screen | neurite-count screen | Eml4 within region
p_ed_row <- (p_soma_volcano | p_neurite_volcano | p_eml4_by_region) +
  patchwork::plot_layout(widths = c(1, 1, 1.25))
ggsave(file.path(DIR_F1_SUPP, "EDFig2_screens_and_Eml4_by_region_ROW.svg"), p_ed_row,
       width = mm_to_in(FIG_W_MM), height = 2.5, device = SVG_DEVICE)

# Soma volume and neurite count by region: why neurite count is restricted to SoC (mPFC and Hippo have a
# ceiling effect) and why Hippo is treated as genuinely smaller in the soma-volume analysis
morph_by_region_plot <- function(y_col, y_lab) {
  ggplot(merged_morph, aes(x = region, y = .data[[y_col]], fill = region)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.5) +
    geom_jitter(width = 0.12, size = 0.9, alpha = 0.6) +
    scale_fill_manual(values = PAL_REGION_THY1, guide = "none") +
    scale_x_discrete(labels = REGION_ABBREV) +
    labs(x = NULL, y = y_lab) +
    THEME_JOURNAL()
}
p_confounds <- morph_by_region_plot("soma_volume", paste0("Soma volume (", SOMA_UNIT, ")")) | morph_by_region_plot("neurite_count", "Neurite count")
ggsave(file.path(DIR_F1_SUPP, "EDFig2_morphology_by_region_confounds.svg"), p_confounds,
       width = 5.0, height = 2.5, device = SVG_DEVICE)

# Reference tables: panel coverage per example cell; the protein panel with its functional groups
coverage_table <- purrr::map_dfr(names(final_sample_cols), function(lbl) {
  n_det <- sum(!is.na(resid_thy1[gene_order_f1, final_sample_cols[[lbl]]]))
  tibble::tibble(cell = lbl, genes_detected = n_det, panel_size = length(gene_order_f1),
                 pct_detected = round(100 * n_det / length(gene_order_f1)))
})
write_csv(coverage_table, file.path(DIR_F1_SUPP, "SuppTable_example_cell_panel_coverage.csv"))
write_csv(tibble::tibble(gene = gene_order_f1, functional_group = unname(gene_to_group[gene_order_f1])),
          file.path(DIR_F1_SUPP, "SuppTable_protein_panel_functional_groups.csv"))
message("Supplementary items saved to ", DIR_F1_SUPP)

###############################################################################
# 9. SUMMARY NUMBERS
# Every number needed for the Figure 1 legend, results text and Extended Data Fig. 2.
# Printed and saved to Figure1_summary_numbers.txt.
###############################################################################
summary_lines <- character()
say <- function(...) { l <- paste0(...); message(l); summary_lines <<- c(summary_lines, l) }
fmt_p <- function(p) formatC(p, format = "e", digits = 1)
# keep_half = TRUE: a median of exactly x.5 (halfway between two whole numbers) is printed as "x-(x+1)", e.g. 4-5 (use for whole-number counts such as neurites)
med_iqr <- function(x, keep_half = FALSE) {
  m <- median(x, na.rm = TRUE)
  m_txt <- if (keep_half && abs(m * 2 - round(m * 2)) < 1e-9 && abs(m - round(m)) > 1e-9) sprintf("%.0f-%.0f", floor(m), ceiling(m)) else sprintf("%.0f", round_half_up(m))
  sprintf("median %s (IQR %.0f-%.0f)", m_txt, round_half_up(quantile(x, 0.25, na.rm = TRUE)), round_half_up(quantile(x, 0.75, na.rm = TRUE)))
}
# spear(): same rules as the screens (safe_spearman); a case that cannot be tested says why instead of failing
spear <- function(x, y, min_n = SENS_MIN_N) {
  r <- safe_spearman(x, y, min_n)
  if (r$status != "tested") return(paste0("not testable (", r$reason, ")"))
  sprintf("rho=%.2f, nominal P=%s, n=%d%s", r$spearman_r, fmt_p(r$p_nominal), r$n, if (isTRUE(r$ties_present)) " (ties: P approximate)" else "")
}
screen_row <- function(tbl, g) { r <- tbl[tbl$gene == g, ]
if (nrow(r) == 0) "not in screen (not testable)" else sprintf("screen: rho=%.2f, nominal P=%s, BH q=%.2f, n=%d", r$spearman_r, fmt_p(r$p_nominal), r$q_BH, r$n) }

say("==================== FIGURE 1 -- SUMMARY STATISTICS ====================")
say("\n--- Data and method ---")
say("Protein matrix: protein_mat_log2 (contaminant-filtered, no normalisation, no imputation)")
say("Depth adjustment: depth_residuals() in setup.R = log2 intensity ~ log2(contaminant-filtered protein count) + batch; ",
    "fitted once across ", ncol(resid_thy1), " Thy1 single cells (all regions, all batches); ",
    n_usable_genes, " of ", nrow(resid_thy1), " proteins have a residual (detected in >= ", MIN_N_DEPTH_FIT, " cells)")
say("Cell inclusion: >= ", MIN_PROTEINS_VIABLE, " contaminant-filtered protein groups (passes_qc_tier in setup.R)")

say("\n--- Cells (morphology-traced Thy1 batch 24Apr2026) ---")
say("Morphology rows in file: ", nrow(morph_raw), "; after dropping the ambiguous pair: ", nrow(morph_clean))
say("Matched to a proteome with >=", MIN_PROTEINS_VIABLE, " protein groups (merged_morph): ", nrow(merged_morph))
say("  per region: ", paste(names(table(merged_morph$region)), table(merged_morph$region), collapse = "; "))
hi <- merged_morph$n_proteins_filtered >= MIN_PROTEINS
say("  with >=", MIN_PROTEINS, " proteins: ", sum(hi), " (per region: ",
    paste(names(table(merged_morph$region[hi])), table(merged_morph$region[hi]), collapse = "; "), ")")
say("Proteins per cell (contaminant-filtered), traced cells: ", med_iqr(merged_morph$n_proteins_filtered),
    "; range ", min(merged_morph$n_proteins_filtered), "-", max(merged_morph$n_proteins_filtered))

say("\n--- Morphology by region (ED Fig 2 confound panels) ---")
for (r in unique(merged_morph$region)) {
  d <- merged_morph[merged_morph$region == r, ]
  say(REGION_ABBREV[[r]], " (n=", nrow(d), "): soma volume ", med_iqr(d$soma_volume), " ", SOMA_UNIT, "; neurite count ", med_iqr(d$neurite_count, keep_half = TRUE))
}
say("Soma volume differs by region: Kruskal-Wallis P=", fmt_p(kruskal.test(soma_volume ~ region, data = merged_morph)$p.value))
say("Neurite count differs by region: Kruskal-Wallis P=", fmt_p(kruskal.test(neurite_count ~ region, data = merged_morph)$p.value))
hip <- merged_morph$region == "Hippocampal Region"
say("Soma volume, Hippo vs cortex (Wilcoxon): P=",
    fmt_p(wilcox.test(merged_morph$soma_volume[hip], merged_morph$soma_volume[!hip])$p.value))

say("\n--- Eml4 vs soma volume (Fig 1e, ED Fig 2d); exploratory candidate, nominal statistics ---")
say("All three regions: ", spear(eml4_df$soma_volume, eml4_df$eml4_residual))
say("Same relationship in the all-gene screen: ", screen_row(soma_corr_table, "Eml4"))
e600 <- eml4_df[eml4_df$n_proteins_filtered >= MIN_PROTEINS, ]
say("Only cells >=", MIN_PROTEINS, " proteins: ", spear(e600$soma_volume, e600$eml4_residual))
for (r in unique(eml4_df$region)) {
  d <- eml4_df[eml4_df$region == r, ]
  sl <- if (nrow(d) >= 2 && dplyr::n_distinct(d$soma_volume) >= 2) unname(coef(lm(eml4_residual ~ soma_volume, data = d))[2]) else NA_real_
  say("  within ", REGION_ABBREV[[r]], ": ", spear(d$soma_volume, d$eml4_residual),
      "; slope ", signif(sl, 2), " (", ifelse(is.na(sl), "not defined", ifelse(sl < 0, "negative", "positive")), ")")
}
say("Eml4 detected in ", nrow(eml4_df), " of ", nrow(merged_morph), " traced cells")

say("\n--- Rtcb vs neurite count, SoC only (Fig 1f); exploratory candidate, nominal statistics ---")
say("All SoC cells: ", spear(rtcb_df$neurite_count, rtcb_df$rtcb_residual))
say("Same relationship in the all-gene screen: ", screen_row(neurite_corr_table, "Rtcb"))
r600 <- rtcb_df[rtcb_df$n_proteins_filtered >= MIN_PROTEINS, ]
say("Only cells >=", MIN_PROTEINS, " proteins: ", spear(r600$neurite_count, r600$rtcb_residual))
say("Neurite count in SoC cells with Rtcb: ", med_iqr(rtcb_df$neurite_count, keep_half = TRUE),
    "; range ", min(rtcb_df$neurite_count), "-", max(rtcb_df$neurite_count))

say("\n--- Gene-wise screens (ED Fig 2c; Supp Tables) ---")
say("Soma volume: ", nrow(soma_corr_table), " proteins tested (>= ", MIN_N_SOMA_SCREEN, " complete pairs), ", nrow(soma_screen$skipped),
    " not tested (reasons: ", paste(names(table(soma_screen$skipped$reason)), table(soma_screen$skipped$reason), sep = " = ", collapse = "; "), "); ",
    sum(soma_corr_table$p_nominal < 0.05),
    " with nominal P<0.05 (uncorrected); ", sum(soma_corr_table$q_BH < 0.1), " with BH q<0.1; ",
    sum(soma_corr_table$p_nominal < 0.05 & soma_corr_table$spearman_r < 0), " negative / ",
    sum(soma_corr_table$p_nominal < 0.05 & soma_corr_table$spearman_r > 0), " positive among P<0.05")
say("Neurite count (SoC): ", nrow(neurite_corr_table), " proteins tested (>= ", MIN_N_NEURITE_SCREEN, " complete pairs), ", nrow(neurite_screen$skipped),
    " not tested (reasons: ", paste(names(table(neurite_screen$skipped$reason)), table(neurite_screen$skipped$reason), sep = " = ", collapse = "; "), "); ",
    sum(neurite_corr_table$p_nominal < 0.05),
    " with nominal P<0.05 (uncorrected); ", sum(neurite_corr_table$q_BH < 0.1), " with BH q<0.1; P approximate (ties) for ",
    sum(neurite_corr_table$ties_present), " of ", nrow(neurite_corr_table), " tested proteins")
say("Top 5 soma hits (by P): ", paste(sprintf("%s (r=%.2f)", head(soma_corr_table$gene, 5), head(soma_corr_table$spearman_r, 5)), collapse = ", "))
say("Top 5 neurite hits (by P): ", paste(sprintf("%s (r=%.2f)", head(neurite_corr_table$gene, 5), head(neurite_corr_table$spearman_r, 5)), collapse = ", "))

say("\n--- Example cells (Fig 1d) and protein panel ---")
say("Panel: ", length(gene_order_f1), " proteins shown of ", length(gene_panel_f1), "; ", length(gene_groups), " functional groups")
print(coverage_table)
say("Example-cell coverage (%): ", paste(coverage_table$cell, coverage_table$pct_detected, sep = "=", collapse = "; "))
ex <- merged_morph %>%
  dplyr::filter(Proteomics_cell_id %in% final_cell_ids) %>%
  dplyr::mutate(label = paste0("Cell ", publication_id),
                proteins_raw = n_proteins_detected,                       # before contaminant removal
                proteins_filtered = as.integer(n_proteins_filtered)) %>%  # after contaminant removal
  dplyr::select(label, Proteomics_cell_id, sample_col, soma_volume, neurite_count, proteins_raw, proteins_filtered) %>%
  dplyr::arrange(match(label, final_labels))
print(as.data.frame(ex))
# Table for Fig. 1a (cells 1, 7, 10, 13)
write_csv(ex %>% dplyr::transmute(Cell = label, `Protein groups (contaminant-filtered)` = proteins_filtered,
                                  !!paste0("Soma volume (", SOMA_UNIT, ")") := soma_volume, `Neurite count` = neurite_count,
                                  Record = Proteomics_cell_id),
          file.path(DIR_F1, "Figure1a_example_cell_metrics_table.csv"))
for (i in seq_len(nrow(ex))) {
  say(ex$label[i], " (", ex$Proteomics_cell_id[i], " | ", ex$sample_col[i], "): soma volume ", ex$soma_volume[i],
      " ", SOMA_UNIT, ", ", ex$neurite_count[i], " neurites, ", ex$proteins_raw[i], " proteins (raw) / ",
      ex$proteins_filtered[i], " proteins (contaminant-filtered)")
}

SW_PKGS <- c("tidyverse", "dplyr", "tidyr", "purrr", "readxl", "ggplot2", "ggrepel", "patchwork", "colorspace", "ggdist", "svglite", "ggrastr", "scales", "arrow")
say("\nSoftware versions (all packages used by this script): ", software_versions(SW_PKGS))
# ---- Source data (one Excel file: Figure 1, Extended Data Figure 1, Extended Data Fig. 2) ----
source_panel("Fig1a_example_cells", ex %>% dplyr::transmute(cell = label, record = Proteomics_cell_id, protein_groups_filtered = proteins_filtered,
                                                            soma_volume_um3 = soma_volume, neurite_count = neurite_count),
             "Example cells: protein groups (contaminant-filtered), soma volume and neurite count")
source_panel("Fig1d_heatmap", plot_df %>% dplyr::transmute(cell = row, protein = gene, functional_group = unname(gene_to_group[as.character(gene)]), depth_adjusted_residual_log2 = residual),
             "Depth-adjusted residuals (log2) of the protein panel in the 4 example cells; proteins not detected are missing rows")
source_panel("Fig1e_Eml4_soma_volume", eml4_df, "Eml4 depth-adjusted residual vs soma volume, all traced neurons with a proteome",
             columns = c("publication_id", "Proteomics_cell_id", "sample_col", "region", "soma_volume", "eml4_residual"))
source_panel("Fig1f_Rtcb_neurites", rtcb_df, "Rtcb depth-adjusted residual vs neurite count, Somatosensory Cortex neurons",
             columns = c("publication_id", "Proteomics_cell_id", "sample_col", "region", "neurite_count", "rtcb_residual"))
source_panel("EDFig1a_heatmap", plot_df_ed %>% dplyr::transmute(cell = row, region = region_abbrev, protein = gene, functional_group = unname(gene_to_group[as.character(gene)]), depth_adjusted_residual_log2 = residual),
             "Depth-adjusted residuals (log2) of the protein panel in all traced neurons; cells marked * have no usable proteome (empty rows); proteins not detected are missing rows")
source_panel("EDFig2_soma_volume_screen", soma_corr_table, "Spearman correlation of every protein's residual with soma volume (volcano plot)")
source_panel("EDFig2_neurite_screen", neurite_corr_table, "Spearman correlation of every protein's residual with neurite count, Somatosensory Cortex (volcano plot)")
source_panel("EDFig2_morphology_by_region", merged_morph, "Soma volume and neurite count per traced neuron and region (Eml4 within region uses the table of Fig. 1e)",
             columns = c("publication_id", "Proteomics_cell_id", "region", "soma_volume", "neurite_count"))
write_source_data(file.path(DIR_F1, "Source_Data_Figure1.xlsx"), "Source data for Figure 1 and Extended Data Figures 1 and 2")

write_software_versions(file.path(DIR_F1, "Figure1_software_versions.txt"), SW_PKGS)
writeLines(summary_lines, file.path(DIR_F1, "Figure1_summary_numbers.txt"))
message("Summary saved to ", file.path(DIR_F1, "Figure1_summary_numbers.txt"))