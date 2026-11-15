###############################################################################
# 04_regional_analysis.R -- Regional differences: Figure 4 and Extended Data (ED) Figs. 6, 7
#
# Thy1-EGFP neurons (3 regions) and Cx3cr1-EGFP microglia (3 regions), single cells with >= MIN_PROTEINS
# protein groups, one animal per cell type (cells are the units). All region pairs per cell type are tested.
#
# Figure 4   b, c  independent GO:BP pathways of neurons / microglia: bars (|NES|) and single-cell pathway heatmap
# ED Fig. 6  a, b  protein detection overlap across regions (UpSet), neurons / microglia
#            c     protein-level heatmap of the pathways, neurons
# ED Fig. 7        protein-level heatmap of the pathways, microglia
#
# Regional DE:  per protein lm(log2 intensity ~ region [+ slide] + log2 depth), variance moderated with limma,
#               BH per comparison. Slide is added when cells come from more than one slide.
# GSEA:         fgsea on GO:BP (msigdbr), ranked by the moderated t. Likely artifacts: >= 30% of the leading edge
#               is ribosomal/proteasomal, blood/immune or neuron/glia marker protein.
# Independent pathways: padj < GSEA_PADJ, not a likely artifact, not redundant (fgsea::collapsePathways).
# Pathway score of a cell: mean over the leading edge of the per-protein z-scored depth-adjusted residuals
#               (lm(y ~ region [+ slide] + depth), depth and slide removed, region effect kept), z-scored per pathway
#               across cells. Display only; the statistics are the GSEA.
#
# Inputs : setup.R
# Outputs: Figure4_Regional_analysis/Figure4_pathways_{neurons,microglia}.svg, Figure4_summary_numbers.txt,
#          Figure4_SUPPLEMENTARY/ (EDFig6_ab_upsets, EDFig6c_protein_heatmap_neurons, EDFig7_protein_heatmap_microglia,
#          Fig4_DE_*.csv, Fig4_GSEA_*.csv, Fig4_independent_pathways.csv, further tables)
###############################################################################
source(file.path(getOption("scp.project_root", getwd()), "Scripts", "shared", "setup.R"))

###############################################################################
# 1. SETUP
###############################################################################
check_packages(c("fgsea", "msigdbr", "limma", "ComplexUpset", "BiocParallel", "data.table", "patchwork"), "04_regional_analysis.R")
require_columns(contamination_annotation, c("gene_label", "matches_ribo_proteasome_pattern", "matches_blood_immune_pattern"),
                "contamination_analysis.csv")
library(ggplot2); library(patchwork)

DE_MIN_DETECTED          <- 0.5
DE_Q_REPORT              <- 0.1
GSEA_MIN_SIZE            <- 10
GSEA_MAX_SIZE            <- 500
GSEA_PADJ                <- 0.05
ARTIFACT_LE_FRACTION     <- 0.30
SCORE_MIN_GENES_FRACTION <- 0.5
UPSET_MIN_DETECTED       <- 0.2
UPSET_Y_MAX              <- 1750
UPSET_H_IN               <- 2.8
N_GENES_SHOWN            <- 3
ROW_H_MM                 <- 3.9
GENE_ROW_MM              <- 2.9
LABEL_MAX_CHARS          <- 46
HEAT_GAP_MM              <- 2.5

COMPARISONS <- tibble::tribble(
  ~cell_type,  ~region_a,                 ~region_b,
  "Thy1",      "Cingulate Cortex",        "Hippocampal Region",
  "Thy1",      "Somatosensory Cortex",    "Hippocampal Region",
  "Thy1",      "Cingulate Cortex",        "Somatosensory Cortex",
  "Microglia", "Cerebellum Gray Matter",  "Somatosensory Cortex",
  "Microglia", "Cerebellum White Matter", "Somatosensory Cortex",
  "Microglia", "Cerebellum Gray Matter",  "Cerebellum White Matter"
) %>%
  dplyr::mutate(label = paste(REGION_ABBREV[region_a], "vs", REGION_ABBREV[region_b]),
                id = paste(cell_type, label), file_id = gsub(" ", "_", id))

REGIONS_BY_TYPE    <- list(Thy1 = REGIONS_THY1, Microglia = REGIONS_MG)
REGION_COL_BY_TYPE <- list(Thy1      = setNames(PAL_REGION_THY1, REGION_ABBREV[names(PAL_REGION_THY1)]),
                           Microglia = setNames(PAL_REGION_MG,   REGION_ABBREV[names(PAL_REGION_MG)]))
HEAT_COLOURS <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "white", "#F4A582", "#D6604D", "#B2182B", "#67001F")
TITLE_BY_TYPE <- c(Thy1 = "Thy1-EGFP neurons", Microglia = "Cx3cr1-EGFP microglia")

###############################################################################
# 2. DATA PREPARATION
###############################################################################
message("\n[Figure 4] Preparing the analysis set...")

sc_meta <- cohort_meta %>%
  dplyr::filter(tier == "tier_1x_single_cell", cell_type_label %in% c("Thy1", "Microglia"),
                n_proteins_filtered >= MIN_PROTEINS,
                (cell_type_label == "Thy1" & region %in% REGIONS_THY1) | (cell_type_label == "Microglia" & region %in% REGIONS_MG)) %>%
  dplyr::filter(sample_col %in% colnames(protein_mat_log2))
if (anyNA(sc_meta$injection_order)) stop("Missing injection_order for cells of the analysis set.", call. = FALSE)
meta_by_type <- list(Thy1 = dplyr::filter(sc_meta, cell_type_label == "Thy1"),
                     Microglia = dplyr::filter(sc_meta, cell_type_label == "Microglia"))
message("  ", nrow(meta_by_type$Thy1), " Thy1 neurons, ", nrow(meta_by_type$Microglia), " microglia")
print(dplyr::count(sc_meta, cell_type_label, region, batch))

gene_flags <- contamination_annotation %>%
  dplyr::transmute(gene = gene_label,
                   flagged_ribo_proteasome = dplyr::coalesce(matches_ribo_proteasome_pattern, FALSE),
                   flagged_blood_immune    = dplyr::coalesce(matches_blood_immune_pattern, FALSE),
                   contamination_flagged   = flagged_ribo_proteasome | flagged_blood_immune)

###############################################################################
# 3. REGIONAL DIFFERENTIAL EXPRESSION
###############################################################################
message("\n[Figure 4] Regional differential expression...")

regional_de <- function(meta, region_a, region_b, mat = protein_mat_log2) {
  cells_a <- meta$sample_col[meta$region == region_a]
  cells_b <- meta$sample_col[meta$region == region_b]
  all_cells <- c(cells_a, cells_b)
  sub <- mat[, all_cells, drop = FALSE]
  detected_a <- rowMeans(!is.na(sub[, cells_a, drop = FALSE]))
  detected_b <- rowMeans(!is.na(sub[, cells_b, drop = FALSE]))
  genes <- rownames(sub)[detected_a >= DE_MIN_DETECTED & detected_b >= DE_MIN_DETECTED]
  
  design <- tibble::tibble(region = as.integer(all_cells %in% cells_a),
                           batch  = meta$batch[match(all_cells, meta$sample_col)],
                           depth  = log2(as.numeric(n_proteins_lookup[all_cells])))
  per_gene <- purrr::map_dfr(genes, function(g) {
    d <- design
    d$y <- as.numeric(sub[g, ])
    d <- d[!is.na(d$y), , drop = FALSE]
    fml <- if (dplyr::n_distinct(d$batch) > 1) y ~ region + factor(batch) + depth else y ~ region + depth
    fit <- tryCatch(lm(fml, data = d), error = function(e) NULL)
    if (is.null(fit) || fit$df.residual <= 0) return(NULL)
    coefs <- summary(fit)$coefficients
    if (!"region" %in% rownames(coefs)) return(NULL)
    s2 <- summary(fit)$sigma^2
    tibble::tibble(gene = g, log2FC = coefs["region", "Estimate"], s2 = s2, df = fit$df.residual,
                   unscaled_var = coefs["region", "Std. Error"]^2 / s2, n_a = sum(d$region == 1), n_b = sum(d$region == 0))
  })
  if (nrow(per_gene) < 10) stop("Regional DE ", region_a, " vs ", region_b, ": fewer than 10 proteins could be fitted.", call. = FALSE)
  if (any(!is.finite(per_gene$s2) | per_gene$s2 <= 0)) stop("Regional DE ", region_a, " vs ", region_b, ": a protein has zero residual variance.", call. = FALSE)
  
  sq <- limma::squeezeVar(per_gene$s2, per_gene$df)
  t_mod <- per_gene$log2FC / sqrt(sq$var.post * per_gene$unscaled_var)
  per_gene$pval <- 2 * stats::pt(-abs(t_mod), per_gene$df + sq$df.prior)
  per_gene$qval <- p.adjust(per_gene$pval, method = "BH")
  per_gene$t_mod <- t_mod
  per_gene %>% dplyr::select(gene, log2FC, t_mod, pval, qval, n_a, n_b) %>% dplyr::arrange(pval)
}

de_results <- list()
for (i in seq_len(nrow(COMPARISONS))) {
  cmp <- COMPARISONS[i, ]
  d <- regional_de(meta_by_type[[cmp$cell_type]], cmp$region_a, cmp$region_b) %>% dplyr::left_join(gene_flags, by = "gene")
  de_results[[cmp$id]] <- d
  message(sprintf("  %-27s: %4d proteins tested | q<%.1f: %d", cmp$id, nrow(d), DE_Q_REPORT, sum(d$qval < DE_Q_REPORT)))
  write_csv(d, file.path(DIR_F4_SUPP, paste0("Fig4_DE_", cmp$file_id, ".csv")))
}

###############################################################################
# 4. GSEA ON GO:BP
###############################################################################
message("\n[Figure 4] GSEA (GO:BP)...")

msig_go <- function(subcollection) {
  if (utils::packageVersion("msigdbr") >= "10.0.0") {
    msigdbr::msigdbr(species = "Mus musculus", collection = "C5", subcollection = subcollection)
  } else {
    msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = subcollection)
  }
}
go_bp <- msig_go("GO:BP")
go_pathways <- split(go_bp$gene_symbol, go_bp$gs_name)
message("  GO:BP gene sets: ", length(go_pathways), " (msigdbr ", as.character(utils::packageVersion("msigdbr")), ")")

NEURON_GLIA_MARKERS <- c("Nefl", "Nefm", "Nefh", "Ina", "Plp1", "Mbp", "Mog", "Mag", "Cnp", "Gfap", "Vim", "S100b",
                         "Tppp", "Sirt2", "Nfasc", "Bcan")
flag_ribo  <- setNames(gene_flags$flagged_ribo_proteasome, gene_flags$gene)
flag_blood <- setNames(gene_flags$flagged_blood_immune, gene_flags$gene)
leading_edge_share <- function(genes, flagged) round(100 * mean(flagged[genes] %in% TRUE), 1)
audit_leading_edge <- function(res) {
  res$pct_ribo_proteasome <- vapply(res$leadingEdge, leading_edge_share, numeric(1), flagged = flag_ribo)
  res$pct_blood_immune    <- vapply(res$leadingEdge, leading_edge_share, numeric(1), flagged = flag_blood)
  res$pct_neuron_glia     <- vapply(res$leadingEdge, function(g) round(100 * mean(g %in% NEURON_GLIA_MARKERS), 1), numeric(1))
  res$likely_artifact <- res$pct_ribo_proteasome >= 100 * ARTIFACT_LE_FRACTION | res$pct_blood_immune >= 100 * ARTIFACT_LE_FRACTION |
    res$pct_neuron_glia >= 100 * ARTIFACT_LE_FRACTION
  res$likely_artifact[is.na(res$likely_artifact)] <- FALSE
  res
}

run_gsea <- function(de_df, label) {
  ranks <- de_df %>% dplyr::filter(is.finite(t_mod)) %>% dplyr::arrange(dplyr::desc(t_mod))
  rank_vec <- setNames(ranks$t_mod, ranks$gene)
  set.seed(42)
  res <- fgsea::fgsea(pathways = go_pathways, stats = rank_vec, minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE, eps = 0,
                      BPPARAM = BiocParallel::SerialParam()) %>%
    tibble::as_tibble() %>% dplyr::arrange(pval) %>% audit_leading_edge()
  near <- res$pathway[!is.na(res$padj) & res$padj < 0.15]
  res$is_redundant <- FALSE
  if (length(near) > 1) {
    collapsed <- fgsea::collapsePathways(data.table::as.data.table(res[res$pathway %in% near, ]), pathways = go_pathways, stats = rank_vec)
    res$is_redundant <- res$pathway %in% near & !res$pathway %in% collapsed$mainPathways
  }
  message(sprintf("  %-27s: %d GO:BP terms tested, %d with padj<%.2f", label, nrow(res), sum(res$padj < GSEA_PADJ, na.rm = TRUE), GSEA_PADJ))
  res
}

gsea_results <- list()
for (i in seq_len(nrow(COMPARISONS))) {
  id <- COMPARISONS$id[i]
  gsea_results[[id]] <- run_gsea(de_results[[id]], id)
  write_csv(dplyr::mutate(gsea_results[[id]], leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")),
            file.path(DIR_F4_SUPP, paste0("Fig4_GSEA_", COMPARISONS$file_id[i], ".csv")))
}

###############################################################################
# 5. INDEPENDENT PATHWAYS
###############################################################################
pathway_label <- function(pathway) {
  x <- tolower(gsub("_", " ", sub("^GOBP_", "", pathway)))
  x <- gsub("\\bmrna\\b", "mRNA", gsub("\\btrna\\b", "tRNA", gsub("\\brna\\b", "RNA", gsub("\\bdna\\b", "DNA", gsub("\\batp\\b", "ATP", x)))))
  x <- ifelse(grepl("^(mRNA|tRNA)", x), x, paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))))
  ifelse(nchar(x) > LABEL_MAX_CHARS, paste0(substr(x, 1, LABEL_MAX_CHARS - 1), "…"), x)
}

all_gsea <- dplyr::bind_rows(lapply(seq_len(nrow(COMPARISONS)), function(i) {
  g <- gsea_results[[COMPARISONS$id[i]]]
  g$leadingEdge <- vapply(g$leadingEdge, paste, character(1), collapse = ";")
  g$cell_type <- COMPARISONS$cell_type[i]; g$comparison <- COMPARISONS$label[i]; g$comparison_id <- COMPARISONS$id[i]
  g$region_a <- COMPARISONS$region_a[i]; g$region_b <- COMPARISONS$region_b[i]
  g
}))
all_gsea$is_redundant[is.na(all_gsea$is_redundant)] <- FALSE

selection_counts <- all_gsea %>%
  dplyr::group_by(cell_type, comparison) %>%
  dplyr::summarise(tested = dplyr::n(),
                   significant = sum(padj < GSEA_PADJ, na.rm = TRUE),
                   artifact = sum(padj < GSEA_PADJ & likely_artifact, na.rm = TRUE),
                   redundant = sum(padj < GSEA_PADJ & !likely_artifact & is_redundant, na.rm = TRUE),
                   independent = sum(padj < GSEA_PADJ & !likely_artifact & !is_redundant, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(factor(cell_type, c("Thy1", "Microglia")), match(comparison, COMPARISONS$label))
write_csv(selection_counts, file.path(DIR_F4_SUPP, "Fig4_pathway_selection_counts.csv"))

indep <- all_gsea %>%
  dplyr::filter(!is.na(padj), padj < GSEA_PADJ, !likely_artifact, !is_redundant) %>%
  dplyr::mutate(higher_region = ifelse(NES > 0, REGION_ABBREV[region_a], REGION_ABBREV[region_b]),
                label = pathway_label(pathway),
                genes_text = vapply(strsplit(leadingEdge, ";", fixed = TRUE), function(v) paste(head(v, N_GENES_SHOWN), collapse = ", "), ""),
                comparison = factor(comparison, levels = COMPARISONS$label)) %>%
  dplyr::arrange(factor(cell_type, c("Thy1", "Microglia")), comparison, padj)
message("  Independent pathways: ", nrow(indep), " (", paste(names(table(indep$cell_type)), table(indep$cell_type), collapse = ", "), ")")
write_csv(dplyr::select(indep, cell_type, comparison, pathway, label, padj, NES, higher_region, size, leadingEdge),
          file.path(DIR_F4_SUPP, "Fig4_independent_pathways.csv"))

###############################################################################
# 6. DEPTH-ADJUSTED PROTEIN VALUES AND PATHWAY SCORES (display only)
###############################################################################
region_kept_matrix <- function(mat, cells, m) {
  depth <- log2(as.numeric(n_proteins_lookup[cells]))
  out <- matrix(NA_real_, nrow(mat), length(cells), dimnames = list(rownames(mat), cells))
  for (g in rownames(mat)) {
    y <- as.numeric(mat[g, cells]); ok <- !is.na(y)
    if (sum(ok) < MIN_N_DEPTH_FIT) next
    d <- data.frame(y = y[ok], depth = depth[ok], region = factor(m$region[ok]), batch = factor(m$batch[ok]))
    d$region <- droplevels(d$region); d$batch <- droplevels(d$batch)
    if (nlevels(d$region) < 2) next
    f <- if (nlevels(d$batch) > 1) lm(y ~ region + batch + depth, data = d) else lm(y ~ region + depth, data = d)
    rc <- grep("^region", names(coef(f)))
    out[g, ok] <- residuals(f) + as.numeric(model.matrix(f)[, rc, drop = FALSE] %*% coef(f)[rc])
  }
  out
}

score_set <- function(z, genes, n_min = SCORE_MIN_GENES_FRACTION * length(genes)) {
  zz <- z[intersect(genes, rownames(z)), , drop = FALSE]
  v <- colMeans(zz, na.rm = TRUE)
  v[colSums(!is.na(zz)) < n_min] <- NA
  v
}

score_pathways <- function(cell_type) {
  P <- dplyr::filter(indep, cell_type == !!cell_type)
  if (nrow(P) == 0) return(NULL)
  P$comparison <- droplevels(P$comparison)
  P$row <- paste(P$comparison, P$pathway)
  regions <- REGIONS_BY_TYPE[[cell_type]]
  m <- meta_by_type[[cell_type]]; m <- m[m$region %in% regions, ]
  le <- strsplit(P$leadingEdge, ";", fixed = TRUE)
  genes <- intersect(unique(unlist(le)), rownames(protein_mat_log2))
  resid <- region_kept_matrix(protein_mat_log2[genes, , drop = FALSE], m$sample_col, m)
  z <- t(scale(t(resid)))
  sc <- vapply(seq_len(nrow(P)), function(i) score_set(z, le[[i]]), numeric(nrow(m)))
  sc <- matrix(sc, nrow = nrow(m), dimnames = list(m$sample_col, P$row))
  list(P = P, meta = m, z = z, resid = resid, score_z = apply(sc, 2, function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)))
}

make_panel <- function(cell_type) {
  S <- score_pathways(cell_type)
  if (is.null(S)) return(NULL)
  P <- S$P; regions <- REGIONS_BY_TYPE[[cell_type]]
  m <- S$meta
  m <- m[order(match(m$region, regions), n_proteins_lookup[m$sample_col]), ]
  m$x <- seq_len(nrow(m))
  region_n <- table(REGION_ABBREV[m$region])
  region_lab_n <- setNames(paste0(names(region_n), "\n(n=", as.integer(region_n), ")"), names(region_n))
  m$region_f <- factor(region_lab_n[REGION_ABBREV[m$region]], levels = region_lab_n[REGION_ABBREV[regions]])
  Z <- S$score_z[m$sample_col, , drop = FALSE]
  
  heat <- as.data.frame(Z) %>% tibble::rownames_to_column("sample_col") %>%
    tidyr::pivot_longer(-sample_col, names_to = "row", values_to = "z") %>%
    dplyr::left_join(dplyr::select(m, sample_col, x, region_f), by = "sample_col") %>%
    dplyr::left_join(dplyr::select(P, row, comparison), by = "row")
  heat$row <- factor(heat$row, levels = rev(P$row))
  row_lab <- setNames(P$label, P$row)
  bars <- P; bars$row <- factor(bars$row, levels = rev(P$row)); bars$absNES <- abs(bars$NES); bars$blk <- "|NES|"
  
  build <- function(z_lim) {
    p_bar <- ggplot(bars, aes(absNES, row, fill = higher_region)) +
      geom_col(width = 0.72) +
      geom_text(aes(label = genes_text), hjust = 0, nudge_x = 0.04, size = pt_to_geom(PT_MIN), family = "Arial", fontface = "italic", colour = "grey35") +
      facet_grid(comparison ~ blk, scales = "free_y", space = "free_y", switch = "y",
                 labeller = labeller(comparison = function(x) gsub(" vs ", "\nvs\n", x))) +
      scale_fill_manual(values = REGION_COL_BY_TYPE[[cell_type]], name = "Higher in") +
      scale_y_discrete(labels = row_lab) +
      scale_x_continuous(limits = c(0, max(bars$absNES) * 2.3), breaks = 0:2, expand = c(0, 0)) +
      labs(x = "|NES|", y = NULL) +
      THEME_JOURNAL() +
      theme(strip.placement = "outside",
            strip.text.y.left = element_text(angle = 0, hjust = 0.5, face = "bold", size = PT_BASE, lineheight = 0.9),
            strip.background.y = element_blank(),
            strip.background.x = element_rect(fill = "white", colour = NA), strip.text.x = element_text(colour = "white"),
            panel.spacing.y = unit(1.5, "mm"))
    p_heat <- ggplot(heat, aes(x, row, fill = z)) +
      geom_tile(height = 0.86) +
      facet_grid(comparison ~ region_f, scales = "free", space = "free_y", labeller = labeller(comparison = function(x) rep("", length(x)))) +
      scale_fill_gradientn(colours = HEAT_COLOURS, limits = c(-z_lim, z_lim), breaks = seq(-z_lim, z_lim, by = 2), na.value = "grey85",
                           name = "Depth-adjusted pathway score (z)",
                           guide = guide_colourbar(barwidth = unit(24, "mm"), barheight = unit(2.2, "mm"), title.position = "top", title.hjust = 0.5)) +
      labs(x = NULL, y = NULL) +
      THEME_JOURNAL() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(),
            panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.5),
            strip.background.x = element_rect(fill = "grey88", colour = NA), strip.text.x = element_text(face = "bold", size = PT_BASE, lineheight = 0.9),
            strip.text.y = element_blank(), strip.background.y = element_blank(),
            panel.spacing.x = unit(HEAT_GAP_MM, "mm"), panel.spacing.y = unit(1.5, "mm"))
    (p_bar | p_heat) + plot_layout(widths = c(1.1, 1)) +
      plot_annotation(title = TITLE_BY_TYPE[[cell_type]], theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = PT_BASE))) &
      theme(legend.position = "bottom", legend.text = element_text(size = PT_BASE), legend.title = element_text(size = PT_BASE, face = "bold"),
            legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"), legend.key.size = unit(3, "mm"))
  }
  
  dep <- log2(as.numeric(n_proteins_lookup[m$sample_col]))
  depth_check <- do.call(rbind, lapply(regions, function(rg) {
    ii <- which(m$region == rg)
    data.frame(cell_type = cell_type, region = REGION_ABBREV[[rg]], n_cells = length(ii), pathway_row = P$row,
               rho_score_vs_depth = vapply(seq_len(ncol(Z)), function(j)
                 suppressWarnings(cor(Z[ii, j], dep[ii], method = "spearman", use = "pairwise.complete.obs")), numeric(1)))
  }))
  list(S = S, P = P, cells = m, build = build, z_need = max(abs(heat$z), na.rm = TRUE), n_rows = nrow(P), n_comp = nlevels(P$comparison),
       depth_check = depth_check, scores = cbind(m[, c("sample_col", "region")], as.data.frame(Z)))
}

save_figure <- function(fig, name, dir, w_in, h_in) {
  ggsave(file.path(dir, paste0(name, ".svg")), fig, width = w_in, height = h_in, units = "in", device = SVG_DEVICE)
  message("  ", name, ": ", round(w_in * 25.4), " x ", round(h_in * 25.4), " mm")
}

###############################################################################
# 7. PROTEIN-LEVEL HEATMAP (genes named beside the bars)
###############################################################################
make_gene_panel <- function(res, cell_type) {
  P <- res$P; m <- res$cells; R <- res$S$resid
  rows <- list()
  for (i in seq_len(nrow(P))) {
    g <- head(strsplit(P$leadingEdge[i], ";", fixed = TRUE)[[1]], N_GENES_SHOWN)
    for (k in seq_along(g)) {
      v <- if (g[k] %in% rownames(R)) as.numeric(R[g[k], m$sample_col]) else rep(NA_real_, nrow(m))
      rows[[length(rows) + 1]] <- data.frame(comparison = as.character(P$comparison[i]), pathway = i,
                                             pathway_lab = paste(strwrap(P$label[i], 34), collapse = "\n"), gene = g[k], key = paste0(i, "|", g[k]), gene_order = k,
                                             sample_col = m$sample_col, x = m$x, region_f = m$region_f, resid = v, stringsAsFactors = FALSE)
    }
  }
  D <- dplyr::bind_rows(rows)
  comps <- levels(P$comparison)
  build <- function(lim) {
    plots <- lapply(seq_along(comps), function(ci) {
      d <- D[D$comparison == comps[ci], ]
      d$pathway_lab <- factor(d$pathway_lab, levels = unique(d$pathway_lab[order(d$pathway)]))
      d$key <- factor(d$key, levels = rev(unique(d$key[order(d$pathway, d$gene_order)])))
      is_top <- ci == 1
      ggplot(d, aes(x, key, fill = pmax(pmin(resid, lim), -lim))) +
        geom_tile(height = 0.92) +
        facet_grid(pathway_lab ~ region_f, scales = "free", space = "free_y", switch = "y") +
        scale_y_discrete(labels = function(z) sub("^[0-9]+\\|", "", z)) +
        scale_fill_gradientn(colours = HEAT_COLOURS, limits = c(-lim, lim), breaks = seq(-lim, lim, by = if (lim >= 4) 2 else 1),
                             na.value = "grey85", name = "Depth-adjusted residual (log2)",
                             guide = guide_colourbar(barwidth = unit(24, "mm"), barheight = unit(2.2, "mm"), title.position = "top", title.hjust = 0.5)) +
        labs(x = NULL, y = NULL, title = comps[ci]) +
        THEME_JOURNAL() +
        theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(),
              axis.text.y = element_text(size = PT_MIN, face = "italic"),
              panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.4),
              strip.placement = "outside", strip.background.y = element_blank(),
              strip.text.y.left = element_text(angle = 0, hjust = 1, size = PT_MIN, lineheight = 0.9),
              strip.text.x = if (is_top) element_text(face = "bold", size = PT_BASE, lineheight = 0.9) else element_blank(),
              strip.background.x = if (is_top) element_rect(fill = "grey88", colour = NA) else element_blank(),
              plot.title = element_text(face = "bold", size = PT_BASE, hjust = 0), plot.title.position = "plot",
              panel.spacing.x = unit(HEAT_GAP_MM, "mm"), panel.spacing.y = unit(0.8, "mm"))
    })
    n_rows <- vapply(comps, function(cn) sum(D$comparison == cn & D$sample_col == D$sample_col[1]), 1)
    wrap_plots(plots, ncol = 1, heights = n_rows + c(2.5, rep(0.6, length(comps) - 1))) +
      plot_layout(guides = "collect") +
      plot_annotation(title = TITLE_BY_TYPE[[cell_type]], theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = PT_BASE))) &
      theme(legend.position = "bottom", legend.text = element_text(size = PT_BASE), legend.title = element_text(size = PT_BASE, face = "bold"),
            legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"))
  }
  list(build = build, need = as.numeric(max(abs(D$resid), na.rm = TRUE)),
       n_rows = length(unique(D$key)), n_path = nrow(P), n_comp = length(comps), data = D)
}

###############################################################################
# 8. DEPTH RANGE PER REGION
###############################################################################
depth_range_table <- function(res, cell_type) {
  m <- res$S$meta; cells <- m$sample_col
  do.call(rbind, lapply(REGIONS_BY_TYPE[[cell_type]], function(rg) { x <- n_proteins_lookup[cells[m$region == rg]]
  data.frame(cell_type = cell_type, region = REGION_ABBREV[[rg]], n_cells = length(x), protein_groups_min = min(x), median = median(x), max = max(x)) }))
}

###############################################################################
# 9. FIGURE 4 and ED FIGS. 6c, 7
###############################################################################
message("\n[Figure 4] Pathway heatmaps...")
res <- list(Thy1 = make_panel("Thy1"), Microglia = make_panel("Microglia"))
res <- res[!vapply(res, is.null, logical(1))]
Z_LIM <- 2 * ceiling(max(vapply(res, function(r) r$z_need, 1)) / 2)
for (ct in names(res)) {
  r <- res[[ct]]
  save_figure(r$build(Z_LIM), paste0("Figure4_pathways_", c(Thy1 = "neurons", Microglia = "microglia")[[ct]]), DIR_F4,
              mm_to_in(FIG_W_MM), mm_to_in(r$n_rows * ROW_H_MM + r$n_comp * 2 + 29))
}

message("\n[Figure 4] Protein-level heatmaps...")
gene_panels <- lapply(names(res), function(ct) make_gene_panel(res[[ct]], ct)); names(gene_panels) <- names(res)
G_LIM <- ceiling(max(vapply(gene_panels, function(g) g$need, 1)))
ED_NAME <- c(Thy1 = "EDFig6c_protein_heatmap_neurons", Microglia = "EDFig7_protein_heatmap_microglia")
for (ct in names(gene_panels)) {
  g <- gene_panels[[ct]]
  save_figure(g$build(G_LIM), ED_NAME[[ct]], DIR_F4_SUPP, mm_to_in(FIG_W_MM), mm_to_in(g$n_rows * GENE_ROW_MM + g$n_path * 0.8 + g$n_comp * 6 + 24))
}
write_csv(dplyr::bind_rows(lapply(gene_panels, function(g) g$data[, c("comparison", "pathway", "gene", "sample_col", "resid")])),
          file.path(DIR_F4_SUPP, "Fig4_protein_residuals_per_cell.csv"))
write_csv(dplyr::bind_rows(lapply(res, function(r) tidyr::pivot_longer(r$scores, -c(sample_col, region), names_to = "comparison_pathway", values_to = "pathway_score_z"))),
          file.path(DIR_F4_SUPP, "Fig4_independent_pathway_scores_per_cell.csv"))
depth_all <- dplyr::bind_rows(lapply(res, function(r) r$depth_check))
write_csv(depth_all, file.path(DIR_F4_SUPP, "Fig4_independent_pathway_score_vs_depth.csv"))
write_csv(dplyr::bind_rows(lapply(names(res), function(ct) depth_range_table(res[[ct]], ct))), file.path(DIR_F4_SUPP, "Fig4_depth_range_per_region.csv"))

###############################################################################
# 10. ED FIG. 6a, 6b: PROTEIN DETECTION OVERLAP (UpSet)
###############################################################################
message("\n[Figure 4] ED Fig. 6a, 6b...")

region_detection_sets <- function(meta) {
  sets <- lapply(unique(meta$region), function(r) {
    cells <- meta$sample_col[meta$region == r]
    names(which(rowMeans(!is.na(protein_mat_log2[, cells, drop = FALSE])) >= UPSET_MIN_DETECTED))
  })
  setNames(sets, unique(meta$region))
}
build_region_upset_df <- function(sets) {
  df <- tibble::tibble(protein = Reduce(union, sets))
  for (nm in names(sets)) df[[unname(REGION_ABBREV[nm])]] <- df$protein %in% sets[[nm]]
  df
}
make_region_upset <- function(df, fill_col) {
  sizes <- table(do.call(paste, c(df[setdiff(names(df), "protein")], sep = "_")))
  if (max(sizes) >= UPSET_Y_MAX) stop("An UpSet intersection is larger than the y axis (UPSET_Y_MAX).", call. = FALSE)
  axis_theme <- ggplot2::theme(text = ggplot2::element_text(size = PT_BASE), axis.text = ggplot2::element_text(size = PT_BASE),
                               axis.title = ggplot2::element_text(size = PT_BASE, face = "bold"))
  ComplexUpset::upset(
    df, intersect = setdiff(names(df), "protein"), name = "Region", width_ratio = 0.2, height_ratio = 0.5,
    base_annotations = list("Intersection size" = ComplexUpset::intersection_size(fill = fill_col, bar_number_threshold = 1.01, text = list(size = pt_to_geom(PT_MIN))) +
                              ggplot2::scale_y_continuous(limits = c(0, UPSET_Y_MAX), breaks = seq(0, 1600, 400), expand = ggplot2::expansion(mult = c(0, 0)))),
    themes = ComplexUpset::upset_modify_themes(list("intersections_matrix" = axis_theme, "overall_sizes" = axis_theme, "Intersection size" = axis_theme)),
    set_sizes = ComplexUpset::upset_set_size(geom = ggplot2::geom_bar(fill = PAL_ACCENT, width = 0.6)) +
      ggplot2::ylab("Total protein\ngroups") +
      ggplot2::scale_y_continuous(breaks = 2000, expand = ggplot2::expansion(mult = c(0, 0.06))))
}
region_sets <- list(Thy1 = region_detection_sets(meta_by_type$Thy1), Microglia = region_detection_sets(meta_by_type$Microglia))
upset_df <- lapply(region_sets, build_region_upset_df)
write_csv(upset_df$Thy1, file.path(DIR_F4_SUPP, "ED6a_region_overlap_neurons.csv"))
write_csv(upset_df$Microglia, file.path(DIR_F4_SUPP, "ED6b_region_overlap_microglia.csv"))
save_figure(make_region_upset(upset_df$Thy1, PAL_THY1) | make_region_upset(upset_df$Microglia, PAL_MICROGLIA),
            "EDFig6_ab_upsets", DIR_F4_SUPP, mm_to_in(FIG_W_MM), UPSET_H_IN)

###############################################################################
# 11. SUMMARY NUMBERS FOR LEGENDS AND TEXT
###############################################################################
summary_lines <- character()
say <- function(...) { l <- paste0(...); message(l); summary_lines <<- c(summary_lines, l) }
say_table <- function(df) for (l in capture.output(print(as.data.frame(df), row.names = FALSE))) say(l)
fmt_p <- function(p) formatC(p, format = "e", digits = 1)
lab <- function(x) unname(REGION_ABBREV[x])
cells_txt <- function(ct) { m <- meta_by_type[[ct]]; tb <- table(factor(lab(m$region), levels = lab(REGIONS_BY_TYPE[[ct]])))
paste0(nrow(m), " cells (", paste(names(tb), tb, sep = " ", collapse = ", "), ")") }
n_ind <- table(factor(indep$cell_type, c("Thy1", "Microglia")))
n_sig_all <- sum(selection_counts$significant); n_nonart <- sum(selection_counts$significant - selection_counts$artifact)

say("==================== FIGURE 4 / ED FIGS. 6, 7: NUMBERS FOR LEGENDS AND TEXT ====================")
say("Software: R ", getRversion(), " | msigdbr ", as.character(utils::packageVersion("msigdbr")), " (", length(go_pathways), " GO:BP sets) | fgsea ",
    as.character(utils::packageVersion("fgsea")), " | limma ", as.character(utils::packageVersion("limma")))

say("\n--- 1. Cells (analysis set: >= ", MIN_PROTEINS, " contaminant-filtered protein groups) ---")
for (ct in c("Thy1", "Microglia")) {
  m <- meta_by_type[[ct]]
  say(ct, ": ", cells_txt(ct))
  say_table(m %>% dplyr::group_by(region) %>%
              dplyr::summarise(cells = dplyr::n(), slides = dplyr::n_distinct(batch), depth_median = median(n_proteins_filtered),
                               depth_q1 = quantile(n_proteins_filtered, .25), depth_q3 = quantile(n_proteins_filtered, .75), .groups = "drop") %>%
              dplyr::mutate(region = lab(region)))
  tt <- as.data.frame.matrix(table(lab(m$region), m$batch)); say("  cells per region x slide:"); say_table(cbind(region = rownames(tt), tt))
  say("  depth differs between regions (Kruskal-Wallis): P = ", fmt_p(kruskal.test(n_proteins_filtered ~ region, data = m)$p.value))
}

say("\n--- 2. Regional differential expression (flagged proteins kept) ---")
de_overview <- purrr::map_dfr(COMPARISONS$id, function(id) {
  d <- de_results[[id]]
  tibble::tibble(comparison = id, tested = nrow(d), q05 = sum(d$qval < 0.05), q10 = sum(d$qval < DE_Q_REPORT),
                 q10_flagged = sum(d$qval < DE_Q_REPORT & d$contamination_flagged, na.rm = TRUE), q10_and_fc1 = sum(d$qval < DE_Q_REPORT & abs(d$log2FC) >= 1))
})
say_table(de_overview)

say("\n--- 3. Pathway selection (GSEA, GO:BP, padj < ", GSEA_PADJ, ") ---")
say_table(selection_counts)
say("All comparisons: ", n_sig_all, " significant -> ", n_nonart, " after removing likely artifacts -> ", nrow(indep), " independent (neurons ", n_ind[["Thy1"]],
    ", microglia ", n_ind[["Microglia"]], ")")

say("\n--- 4. Independent pathways shown in Fig. 4 (order as in the figure) ---")
say_table(indep %>% dplyr::transmute(cell_type, comparison, pathway = label, higher_in = higher_region, NES = round(NES, 2), padj = signif(padj, 2),
                                     set_size = size, leading_edge_shown = genes_text))
shown_missing <- unlist(lapply(gene_panels, function(g) {
  d <- g$data; all_na <- tapply(is.na(d$resid), d$key, all)
  unique(d$gene[d$key %in% names(all_na)[all_na]]) }))
say("Protein rows shown in ED Figs. 6c, 7: ", paste(vapply(gene_panels, function(g) g$n_rows, 1), collapse = " + "), "; rows with no value in any cell: ",
    if (length(shown_missing)) paste(shown_missing, collapse = ", ") else "none")

say("\n--- 5. Detection overlap (ED Fig. 6a, 6b): protein detected in >= ", 100 * UPSET_MIN_DETECTED, "% of a region's cells ---")
for (ct in c("Thy1", "Microglia")) {
  sets <- region_sets[[ct]]; core <- Reduce(intersect, sets); all_prot <- Reduce(union, sets)
  say(ct, ": ", length(all_prot), " proteins | detected in all ", length(sets), " regions: ", length(core), " (", round(100 * length(core) / length(all_prot), 1), "%)")
  for (r in names(sets)) {
    only <- setdiff(sets[[r]], unlist(sets[setdiff(names(sets), r)]))
    say("   ", lab(r), ": ", length(sets[[r]]), " detected | ", length(only), " only in this region")
  }
}

say("\n--- 6. Colour scales ---")
say("Pathway score (Fig. 4): z-score, shared scale +/-", Z_LIM, "; observed range ", paste(sprintf("%s %.2f", names(res), vapply(res, function(r) r$z_need, 1)), collapse = ", "), " (max |z|)")
say("Protein residual (ED Figs. 6c, 7): log2 units, shared scale +/-", G_LIM, " (covers the full range of the displayed values; nothing clipped)")
say("Gray = protein not detected in the cell (pathway score: fewer than ", 100 * SCORE_MIN_GENES_FRACTION, "% of the leading-edge proteins detected); nothing imputed.")

say("\n--- 7. Depth control (Spearman rho of pathway score with proteome depth, within region) ---")
say_table(depth_all %>% dplyr::group_by(cell_type, region) %>%
            dplyr::summarise(pathways = dplyr::n(), median_rho = round(median(rho_score_vs_depth, na.rm = TRUE), 2),
                             median_abs_rho = round(median(abs(rho_score_vs_depth), na.rm = TRUE), 2),
                             pct_abs_rho_above_0.3 = round(100 * mean(abs(rho_score_vs_depth) > 0.3, na.rm = TRUE)), .groups = "drop"))

SW_PKGS <- c("tidyverse", "dplyr", "ggplot2", "patchwork", "ComplexUpset", "svglite", "ggrastr", "scales", "data.table", "limma", "fgsea", "msigdbr", "BiocParallel")
say("\nSoftware versions: ", software_versions(SW_PKGS))
# ---- Source data (Figure 4, Extended Data Figures 6 and 7) ----
CT_NAME <- c(Thy1 = "neurons", Microglia = "microglia")
for (ct in names(res)) {
  r <- res[[ct]]
  source_panel(paste0("Fig4_pathways_", CT_NAME[[ct]]),
               dplyr::transmute(r$P, comparison = as.character(comparison), pathway, label, NES, padj, set_size = size, higher_region, leading_edge_shown = genes_text),
               paste0("Independent GO:BP pathways, ", CT_NAME[[ct]], ": normalised enrichment score (bar length is |NES|), adjusted P and the region with the higher level"))
  source_panel(paste0("Fig4_pathway_scores_", CT_NAME[[ct]]),
               r$scores %>% tidyr::pivot_longer(-c(sample_col, region), names_to = "row", values_to = "pathway_score_z") %>%
                 dplyr::left_join(dplyr::select(r$P, row, comparison, pathway, label), by = "row") %>%
                 dplyr::left_join(dplyr::select(r$cells, sample_col, cell_order = x), by = "sample_col") %>%
                 dplyr::transmute(sample_col, region = unname(REGION_ABBREV[region]), cell_order, protein_groups = as.numeric(n_proteins_lookup[sample_col]),
                                  comparison = as.character(comparison), pathway, label, pathway_score_z),
               paste0("Depth-adjusted pathway score (z) per ", CT_NAME[[ct]], " cell and pathway; cell_order = column position in the heatmap (by region, then proteome depth)"))
}
source_panel("EDFig6a_overlap_neurons", upset_df$Thy1, paste0("Protein groups detected in >= ", 100 * UPSET_MIN_DETECTED, "% of the neurons of a region (TRUE = present); basis of the UpSet plot"))
source_panel("EDFig6b_overlap_microglia", upset_df$Microglia, paste0("Protein groups detected in >= ", 100 * UPSET_MIN_DETECTED, "% of the microglia of a region (TRUE = present); basis of the UpSet plot"))
for (ct in names(gene_panels)) {
  source_panel(if (ct == "Thy1") "EDFig6c_protein_residuals" else "EDFig7_protein_residuals",
               gene_panels[[ct]]$data %>% dplyr::transmute(comparison, pathway_number = pathway, pathway = gsub("\n", " ", pathway_lab), gene, gene_order, sample_col,
                                                           region = sub("\n.*$", "", as.character(region_f)), cell_order = x, depth_adjusted_residual_log2 = resid),
               paste0("Depth-adjusted, region-kept protein residuals (log2) of the leading-edge proteins shown, ", CT_NAME[[ct]], "; empty = not detected"))
}
write_source_data(file.path(DIR_F4, "Source_Data_Figure4.xlsx"), "Source data for Figure 4 and Extended Data Figures 6 and 7")

write_software_versions(file.path(DIR_F4, "Figure4_software_versions.txt"), SW_PKGS)
writeLines(summary_lines, file.path(DIR_F4, "Figure4_summary_numbers.txt"))
message("\n[Figure 4] Done. Summary: ", file.path(DIR_F4, "Figure4_summary_numbers.txt"))