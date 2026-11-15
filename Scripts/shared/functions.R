###############################################################################
# functions.R -- function definitions shared by all scripts (no data is loaded here)
#
# Sourced automatically by setup.R. Contents:
#   1. Input checks:     check_packages(), require_file(), require_columns()
#   2. Data loading:     load_precursor_matrix()
#   3. Statistics:       depth_residuals() (the single depth-adjustment method), build_overlap_table()
#   4. Plotting:         THEME_OPEN_AXIS(), THEME_JOURNAL(), RAINCLOUD_LAYERS(), mm_to_in(), pt_to_geom(), fix_panel_h()
# Constants these functions use (PT_BASE, MIN_N_DEPTH_FIT, ...) are defined in setup.R.
###############################################################################

# ---- 1. Input checks ---------------------------------------------------------
# check_packages(pkgs, needed_for)
#   Stops with a clear message if any package in `pkgs` is not installed. Never installs anything.
check_packages <- function(pkgs, needed_for = "the analysis scripts") {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R package(s) needed for ", needed_for, ": ", paste(missing, collapse = ", "),
         "\nInstall them once with source(\"install_packages.R\") from the project folder (see README.md).",
         call. = FALSE)
  }
}
# Input checks: stop early with a message that says what is wrong and where.
# require_file(path, what, hint)  -- stops if the file does not exist; returns `path` invisibly.
# require_columns(df, cols, what) -- stops if `df` lacks any of `cols`; returns `df` invisibly.
require_file <- function(path, what, hint = NULL) {
  if (!file.exists(path)) stop(what, " not found: ", path, if (!is.null(hint)) paste0("\n", hint), call. = FALSE)
  invisible(path)
}
require_columns <- function(df, cols, what) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop(what, " is missing required column(s): ", paste(missing, collapse = ", "),
         "\n  columns found: ", paste(colnames(df), collapse = ", "), call. = FALSE)
  }
  invisible(df)
}

# ---- 2. Data loading ---------------------------------------------------------
# precursor_long_to_log2_matrix(pr_long): long table (Run, Precursor.Id, Precursor.Quantity; one row per
# run and precursor) -> precursors x runs matrix of log2 quantity; 0 = not detected -> NA.
precursor_long_to_log2_matrix <- function(pr_long) {
  pr_wide <- pr_long %>%
    dplyr::select(Run, Precursor.Id, Precursor.Quantity) %>%
    tidyr::pivot_wider(names_from = Run, values_from = Precursor.Quantity)
  pr_mat <- as.matrix(pr_wide %>% dplyr::select(-Precursor.Id))
  rownames(pr_mat) <- pr_wide$Precursor.Id
  pr_mat[pr_mat == 0] <- NA
  log2(pr_mat)
}

# resolve_precursor_duplicates(pr_long): a precursor should appear once per run. If a Run + Precursor.Id
# appears several times, this checks whether the rows carry the same Precursor.Quantity.
#   * same quantity in all duplicates -> keeping one row loses nothing;
#   * different quantities -> explicit rule instead of an arbitrary row: keep the row with the lowest
#     Q.Value (best identification confidence), then the lowest PG.Q.Value, then the highest quantity.
# Needs the columns Run, Precursor.Id, Precursor.Quantity, Q.Value, PG.Q.Value.
# Returns list(resolved = one row per Run + Precursor.Id, first_row = the table the OLD rule gave (first row in
# file order), check = one-row summary: n_rows, n_keys_duplicated, n_keys_differing, max_relative_difference, rule).
resolve_precursor_duplicates <- function(pr_long, tol = 1e-9) {
  key <- pr_long[, c("Run", "Precursor.Id")]
  is_dup <- duplicated(key) | duplicated(key, fromLast = TRUE)
  first_row <- dplyr::distinct(pr_long, Run, Precursor.Id, .keep_all = TRUE)
  if (!any(is_dup)) {
    check <- tibble::tibble(n_rows = nrow(pr_long), n_keys_duplicated = 0L, n_keys_differing = 0L,
                            max_relative_difference = 0, rule = "none needed (no duplicates)")
    return(list(resolved = first_row, first_row = first_row, check = check))
  }
  d <- pr_long[is_dup, ] %>%
    dplyr::group_by(Run, Precursor.Id) %>%
    dplyr::summarise(q_min = min(Precursor.Quantity), q_max = max(Precursor.Quantity), .groups = "drop") %>%
    dplyr::mutate(rel_diff = dplyr::if_else(q_max > 0, (q_max - q_min) / q_max, 0))
  n_diff <- sum(d$rel_diff > tol)
  resolved <- if (n_diff == 0) first_row else {
    pr_long %>% dplyr::arrange(Run, Precursor.Id, Q.Value, PG.Q.Value, dplyr::desc(Precursor.Quantity)) %>%
      dplyr::distinct(Run, Precursor.Id, .keep_all = TRUE)
  }
  check <- tibble::tibble(n_rows = nrow(pr_long), n_keys_duplicated = nrow(d), n_keys_differing = n_diff,
                          max_relative_difference = max(d$rel_diff),
                          rule = if (n_diff == 0) "identical quantities: keeping one row is safe"
                          else "different quantities: keep the row with the lowest Q.Value (then PG.Q.Value, then highest quantity)")
  list(resolved = resolved, first_row = first_row, check = check)
}

# load_precursor_matrix(parquet_file)
#   Input : path to the DIA-NN long-format parquet report (raw Precursor.Quantity).
#   Output: precursors x runs matrix of log2 Precursor.Quantity (NA = not detected), column names =
#           `Run`; or NULL if the file does not exist (the input-amount panels are then skipped).
#   Only used for the input-amount panels (1x vs 10x). 
load_precursor_matrix <- function(parquet_file) {
  message("Loading precursor-level matrix (raw Precursor.Quantity)...")
  if (!file.exists(parquet_file)) {
    message("  WARNING: precursor parquet file not found at ", parquet_file,
            " -- the input-amount panels will be skipped.")
    precursor_mat_log2 <- NULL
  } else {
    # The parquet file exists, so arrow is required (the panels are not skipped silently).
    check_packages("arrow", "reading the DIA-NN parquet report")
    
    pr_long <- arrow::read_parquet(
      parquet_file,
      col_select = c("Run", "Precursor.Id", "Precursor.Quantity",
                     "Q.Value", "PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "PG.PEP"))
    n_before_fdr <- nrow(pr_long)
    
    # Keep confident identifications only (standard DIA-NN report filters)
    pr_long <- pr_long %>%
      dplyr::filter(Q.Value <= 0.01, PG.Q.Value <= 0.01, Lib.Q.Value <= 0.01,
                    Lib.PG.Q.Value <= 0.01, PG.PEP <= 0.15)
    message("  Precursor rows: ", n_before_fdr, " -> ", nrow(pr_long), " after FDR filtering")
    
    # Duplicate Run + Precursor.Id records: checked and resolved by an explicit rule (see
    # resolve_precursor_duplicates() below), then long -> wide -> matrix; 0 = not detected -> NA; log2.
    # `Run` equals raw_filename in sample_metadata without the ".d" suffix.
    dup <- resolve_precursor_duplicates(pr_long)
    message("  Duplicate Run + Precursor.Id records: ", dup$check$n_keys_duplicated, " keys duplicated, ",
            dup$check$n_keys_differing, " with different quantities (max relative difference ",
            signif(dup$check$max_relative_difference, 3), "). Rule: ", dup$check$rule)
    precursor_mat_log2 <- precursor_long_to_log2_matrix(dup$resolved)
    attr(precursor_mat_log2, "duplicate_check") <- dup$check
    # only when the quantities differ: the matrix made with the OLD rule (first row in file order), so that the
    # effect of the rule on the precursor-level results can be shown
    if (dup$check$n_keys_differing > 0) attr(precursor_mat_log2, "first_row_matrix") <- precursor_long_to_log2_matrix(dup$first_row)
    message("  Precursor matrix: ", nrow(precursor_mat_log2), " precursors x ", ncol(precursor_mat_log2), " runs")
  }
  precursor_mat_log2
}

# ---- 3. Statistics -----------------------------------------------------------
# depth_residuals(cells, min_n, ...): depth adjustment used in every figure.
# For every protein, across the cells given:
#     log2 intensity ~ log2(proteins detected)   [+ batch, if the cells span more than one batch]
# The residual says how much higher or lower a cell is than expected for its depth (and batch).
# Returns a proteins x cells matrix of residuals (NA = not detected or no fit).
# Arguments: cells = sample names (columns of protein_mat); min_n = minimum number of cells with a
#   detection. protein_mat, meta_table and n_proteins default to the objects built in setup.R.
depth_residuals <- function(cells, min_n = MIN_N_DEPTH_FIT, protein_mat = protein_mat_log2,
                            meta_table = cohort_meta, n_proteins = n_proteins_lookup, use_batch = TRUE) {
  if (length(cells) == 0) stop("depth_residuals(): no cells given.", call. = FALSE)
  if (!all(cells %in% colnames(protein_mat))) {
    stop("depth_residuals(): cells not in the protein matrix: ", paste(setdiff(cells, colnames(protein_mat)), collapse = ", "), call. = FALSE)
  }
  meta <- meta_table[match(cells, meta_table$sample_col), ]
  stopifnot(!anyNA(meta$batch))
  d_all <- data.frame(depth = as.numeric(log2(n_proteins[cells])), batch = factor(meta$batch))
  res <- matrix(NA_real_, nrow(protein_mat), length(cells),
                dimnames = list(rownames(protein_mat), cells))
  for (g in rownames(protein_mat)) {
    y  <- protein_mat[g, cells]
    ok <- !is.na(y)
    if (sum(ok) < min_n) next
    d <- d_all[ok, , drop = FALSE]
    d$y <- y[ok]
    d$batch <- droplevels(d$batch)
    fit <- if (use_batch && nlevels(d$batch) > 1) lm(y ~ depth + batch, data = d) else lm(y ~ depth, data = d)
    res[g, ok] <- residuals(fit)
  }
  res
}

# safe_spearman(x, y, min_n, what): one robust Spearman correlation, used by every protein-wide screen,
# candidate-protein plot and summary/sensitivity number so that they all follow the same rules.
#   1. keeps only complete, finite pairs of x and y
#   2. needs at least `min_n` such pairs
#   3. needs at least 2 distinct values in x AND in y among those pairs
#   4. runs cor.test(method = "spearman") and catches errors; warnings are recorded, not hidden
#   5. accepts the result only if rho and P are finite and valid
# Returns a one-row tibble: n (complete pairs), spearman_r, p_nominal (NOT multiple-testing corrected),
# ties_present (TRUE: R cannot compute an exact P and uses an approximation), warning (text or NA),
# status ("tested" / "skipped") and reason (why it was skipped, else NA). Never stops the calling loop.
# `what` = names of x and y, used only in the skip reasons.
safe_spearman <- function(x, y, min_n, what = c("x", "y")) {
  res <- function(status, reason = NA_character_, n = NA_integer_, rho = NA_real_, p = NA_real_,
                  ties = NA, warn = NA_character_) {
    tibble::tibble(n = as.integer(n), spearman_r = rho, p_nominal = p, ties_present = ties,
                   warning = warn, status = status, reason = reason)
  }
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < min_n) return(res("skipped", paste0("fewer than ", min_n, " complete pairs"), n))
  if (length(unique(x)) < 2) return(res("skipped", paste0(what[1], " constant among the complete pairs"), n))
  if (length(unique(y)) < 2) return(res("skipped", paste0(what[2], " constant among the complete pairs"), n))
  ties <- anyDuplicated(x) > 0 || anyDuplicated(y) > 0
  warns <- character()
  ct <- tryCatch(
    withCallingHandlers(stats::cor.test(x, y, method = "spearman"),
                        warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)
  if (inherits(ct, "error")) return(res("skipped", paste0("cor.test failed: ", conditionMessage(ct)), n))
  rho <- unname(ct$estimate); p <- ct$p.value
  # For small n without ties R uses an Edgeworth approximation that can return P = 0 for strong correlations.
  # A P value of exactly 0 is impossible, so such a result is replaced by the t approximation (as used when ties are present).
  if (is.finite(p) && p <= 0) {
    ct2 <- tryCatch(suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE)), error = function(e) NULL)
    if (!is.null(ct2) && is.finite(ct2$p.value) && ct2$p.value > 0) p <- ct2$p.value
  }
  if (length(rho) != 1 || length(p) != 1 || !is.finite(rho) || !is.finite(p) || abs(rho) > 1 || p < 0 || p > 1) {
    return(res("skipped", "correlation estimate or P value not finite or not valid", n))
  }
  res("tested", NA_character_, n, rho, p, ties, if (length(warns)) paste(unique(warns), collapse = "; ") else NA_character_)
}

# fix_panel_h(p, h_in): makes the plotting area (inside the axes) of plot `p` exactly h_in inches high,
# so that y axes have the same length whatever the axis text, legend or title around the plot. The
# plot sits at the top of its patchwork cell. Facets share the height. Returns a patchwork element.
fix_panel_h <- function(p, h_in) {
  g <- ggplot2::ggplotGrob(p)
  rows <- unique(g$layout$t[grepl("^panel", g$layout$name)])
  w <- suppressWarnings(as.numeric(g$heights[rows])); if (anyNA(w) || sum(w) <= 0) w <- rep(1, length(rows))
  g$heights[rows] <- grid::unit(h_in * w / sum(w), "in")
  patchwork::wrap_elements(full = grid::grobTree(g, vp = grid::viewport(y = grid::unit(1, "npc"), just = "top",
                                                                        height = sum(g$heights))))
}

# build_overlap_table(sets): which proteins are found in which datasets (used for the UpSet plot).
#
# The problem: the same protein is often written differently in different tables.
# A simple match on one column would count such a protein twice.
#
# The idea: every row ("entry") gets one or more "tokens": its accession(s) and its gene symbol (written "G:SYMBOL").
# Two entries that share ANY token are treated as the same protein, and this chains: if A shares a token with B
# and B shares one with C, then A, B and C are one protein.
#
# Input  sets: named list, one data.frame per dataset, with columns
#              `acc`  (UniProt accession(s), separated by ';' or ',', may be NA) and
#              `gene` (gene symbol, may be NA).
# Output tibble: one row per protein, a label (`protein`), and one TRUE/FALSE column per dataset
#              (TRUE = this protein is in that dataset). The input entries are attached as attr(, "entries").
build_overlap_table <- function(sets) {
  # Step 1: stack all datasets into one long table (one row per entry, with the name of its dataset).
  ent <- dplyr::bind_rows(lapply(names(sets), function(nm) {
    s <- sets[[nm]]
    data.frame(dataset = nm, acc = as.character(s$acc), gene = as.character(s$gene), stringsAsFactors = FALSE)
  }))
  # Step 2: tidy the names. Gene: keep only the first symbol, upper case, no spaces; empty or "NA" becomes a real NA.
  ent$gene <- toupper(trimws(sub("[; ].*$", "", ent$gene)))
  ent$gene[is.na(ent$gene) | ent$gene %in% c("", "NA")] <- NA_character_
  # Accessions: split "P1;P2" into separate accessions and remove the isoform suffix ("P12345-2" -> "P12345").
  acc_list <- lapply(strsplit(ifelse(is.na(ent$acc), "", ent$acc), "[;,]"), function(a) {
    a <- trimws(a); a <- sub("-\\d+$", "", a); unique(a[a != "" & a != "NA"]) })
  # Step 3: fill in missing accessions. A dictionary "gene symbol -> accessions" is built from entries that have both.
  # Entries that have only a gene symbol take their accessions from it, so they can match entries from other tables.
  has_both <- which(!is.na(ent$gene) & lengths(acc_list) > 0)
  bridge <- lapply(split(unlist(acc_list[has_both]), rep(ent$gene[has_both], lengths(acc_list[has_both]))), unique)
  for (i in which(!is.na(ent$gene) & lengths(acc_list) == 0)) {
    b <- bridge[[ent$gene[i]]]; if (!is.null(b)) acc_list[[i]] <- b
  }
  # Step 4: the tokens of each entry = its accessions + its gene symbol (with the prefix "G:").
  # Entries without any token (no accession and no gene) cannot be matched and are dropped later (has_tok = FALSE).
  tok <- lapply(seq_len(nrow(ent)), function(i) c(acc_list[[i]], if (!is.na(ent$gene[i])) paste0("G:", ent$gene[i])))
  has_tok <- lengths(tok) > 0
  # Step 5: group entries that share a token. `edge` is a list of (entry, token) pairs. Every entry starts in its own
  # group (its row number). Each round, an entry takes the smallest group number among all entries that share a token
  # with it. This repeats until nothing changes; entries that end with the same number are the same protein.
  edge <- data.frame(entry = rep(seq_len(nrow(ent)), lengths(tok)), token = unlist(tok), stringsAsFactors = FALSE)
  comp <- seq_len(nrow(ent))
  repeat {                                         # (technical name: connected components by label propagation)
    tmin <- tapply(comp[edge$entry], edge$token, min)
    new  <- tapply(tmin[edge$token], edge$entry, min)
    nc <- comp; nc[as.integer(names(new))] <- as.integer(new)
    if (identical(nc, comp)) break
    comp <- nc
  }
  # Step 6: build the output table. `comp` = protein id; one TRUE/FALSE column per dataset.
  ent$comp <- comp
  ent$has_tok <- has_tok
  ent_ok <- ent[has_tok, ]
  out <- data.frame(comp = sort(unique(ent_ok$comp)))
  for (nm in names(sets)) out[[nm]] <- out$comp %in% ent_ok$comp[ent_ok$dataset == nm]
  # Step 7: label each protein with its first gene symbol, or "ACC:<accession>" if it has no symbol.
  first_gene <- tapply(ent_ok$gene, ent_ok$comp, function(g) { g <- g[!is.na(g)]; if (length(g)) g[1] else NA_character_ })
  first_acc  <- tapply(ent_ok$acc,  ent_ok$comp, function(a) { a <- a[!is.na(a)]; if (length(a)) sub("[;,].*$", "", a[1]) else NA_character_ })
  lab <- first_gene[as.character(out$comp)]
  lab[is.na(lab)] <- paste0("ACC:", first_acc[as.character(out$comp)][is.na(lab)])
  # Sanity message: a group much larger than the number of datasets suggests that unrelated proteins got merged.
  sizes <- table(ent_ok$comp)
  message("  build_overlap_table: ", nrow(ent_ok), " entries -> ", nrow(out), " proteins; largest merged group = ", max(sizes),
          " entries (", sum(sizes > 2 * length(sets)), " unusually large group(s))")
  res <- tibble::as_tibble(cbind(protein = unname(lab), out[names(sets)]))
  attr(res, "entries") <- ent                         # every input entry with its merged-protein id (comp) and has_tok
  res
}

# merged_entries_report(tab, study): explains why a dataset has more entries than proteins.
#   Entries are merged when they share an accession or gene symbol (two rows for the same protein), and
#   entries with neither an accession nor a gene symbol are dropped.
#   Input: `tab` = the result of build_overlap_table(); `study` = the dataset to report on.
#   Returns list(n_entries = input rows, n_proteins = proteins after merging, n_dropped = rows without any
#   accession or symbol, merged = table of the entries that were merged into the same protein).
merged_entries_report <- function(tab, study = "3D-SCP") {
  e <- attr(tab, "entries"); e <- e[e$dataset == study, , drop = FALSE]
  ok <- e[e$has_tok, , drop = FALSE]
  dup <- unique(ok$comp[duplicated(ok$comp)])
  m <- ok[ok$comp %in% dup, c("comp", "acc", "gene"), drop = FALSE]
  m <- m[order(m$comp), , drop = FALSE]
  list(n_entries = nrow(e), n_proteins = length(unique(ok$comp)), n_dropped = sum(!e$has_tok), merged = m)
}
# ---- 3b. Software versions ---------------------------------------------------
# software_versions(pkgs)       one line "R x.y.z | pkg 1.2.3 | ..." for the Methods; not-installed packages are marked.
# write_software_versions(file, pkgs)  writes R version, system, the package table and the full sessionInfo() to `file`.
software_versions <- function(pkgs) {
  ver <- vapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "NOT INSTALLED", character(1))
  paste0("R ", getRversion(), " | ", paste(pkgs, ver, sep = " ", collapse = " | "))
}
write_software_versions <- function(file, pkgs) {
  ver <- vapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "NOT INSTALLED", character(1))
  writeLines(c(paste0("Written: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
               R.version.string, paste0("Platform: ", R.version$platform), "",
               "Packages used by this script:", paste0("  ", format(pkgs), "  ", ver), "",
               "Full sessionInfo():", utils::capture.output(print(utils::sessionInfo()))), file)
  invisible(ver)
}

# ---- 4. Plotting -------------------------------------------------------------
# THEME_OPEN_AXIS(base_size, base_family): ggplot theme with open (L-shaped) axes.
THEME_OPEN_AXIS <- function(base_size = 10, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line  = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4)
    )
}

# mm_to_in(x): millimetres -> inches (ggsave() takes inches; figure sizes are written in mm).
mm_to_in <- function(x) x / 25.4
# pt_to_geom(x): font size in pt -> the `size` argument of geom_text() / geom_text_repel() (mm).
pt_to_geom <- function(x) x / ggplot2::.pt

# THEME_JOURNAL(): journal specification (Arial, PT_BASE pt text, bold axis titles) on top of THEME_OPEN_AXIS.
THEME_JOURNAL <- function() {
  THEME_OPEN_AXIS(PT_BASE) +
    ggplot2::theme(
      axis.title   = ggplot2::element_text(face = "bold", size = PT_BASE),
      axis.text    = ggplot2::element_text(size = PT_BASE),
      legend.text  = ggplot2::element_text(size = PT_BASE),
      legend.title = ggplot2::element_text(size = PT_BASE, face = "bold"),
      strip.text   = ggplot2::element_text(size = PT_BASE, face = "bold"),
      plot.title   = ggplot2::element_text(size = PT_BASE, face = "bold"),
      legend.key.size = grid::unit(3, "mm")
    )
}

# Raincloud plot layers: half-violin + boxplot + jittered points.
# (`dots_binwidth` is not used; it is kept so that older calls do not fail.)
RAINCLOUD_LAYERS <- function(dots_binwidth, point_size = 1.5,
                             point_alpha = 0.2, violin_alpha = 0.6,
                             violin_width = 0.6) {
  list(
    ggdist::stat_halfeye(adjust = 0.5, width = violin_width, justification = -0.3,
                         .width = 0, point_colour = NA, alpha = violin_alpha),
    ggplot2::geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.85),
    ggplot2::geom_point(size = point_size, alpha = point_alpha,
                        position = ggplot2::position_jitter(seed = 1, width = 0.1))
  )
}
# SVG export device for ggsave(device = ...); names Arial explicitly. Uses svglite when installed (text stays text),
# otherwise the built-in cairo SVG (text as outlines).
SVG_DEVICE <- function(file, width, height, ...) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = file, width = width, height = height, system_fonts = list(sans = "Arial"), fix_text_size = FALSE)
  } else {
    grDevices::svg(filename = file, width = width, height = height, family = "Arial")
  }
}

# geom_point_dense(): same arguments as geom_point(), for layers with thousands of points. With ggrastr installed the points are drawn
# as one 600 dpi bitmap (axes, text and boxes stay vector), so the SVG stays small; without it,
# ordinary points are drawn (and a note is printed once).
geom_point_dense <- function(..., raster.dpi = 600) {
  if (requireNamespace("ggrastr", quietly = TRUE)) return(ggrastr::geom_point_rast(..., raster.dpi = raster.dpi))
  if (is.null(getOption("scp.dense_note"))) { message("NOTE: package ggrastr not installed; dense point layers are drawn as vector points (install.packages('ggrastr'))"); options(scp.dense_note = TRUE) }
  ggplot2::geom_point(...)
}

# ---- Source data -----------------------------------------------------------------------------------------------------------
# source_panel("Fig2c", plot, "what the panel shows") stores the data behind a panel; write_source_data(file, title) writes all
# stored panels into ONE Excel file (one sheet per panel, plus a README sheet) and empties the store. A ggplot panel contributes
# the data frame it was built from (plot data and the data of its layers); a patchwork panel contributes each of its sub-plots;
# a data frame is stored as it is. NULL (a skipped panel) is ignored.
.SOURCE_STORE <- new.env()
.SOURCE_STORE$panels <- list()

.plot_tables <- function(p) {
  if (inherits(p, "patchwork")) {
    subs <- c(p$patches$plots, list(p))
    out <- list()
    for (s in subs) out <- c(out, .plot_tables(if (identical(s, p)) { q <- p; q$patches <- NULL; class(q) <- setdiff(class(q), "patchwork"); q } else s))
    return(out)
  }
  tabs <- list()
  if (is.data.frame(p$data) && nrow(p$data) > 0) tabs[[1]] <- p$data
  for (l in p$layers) if (is.data.frame(l$data) && nrow(l$data) > 0) tabs[[length(tabs) + 1]] <- l$data
  tabs <- tabs[!duplicated(lapply(tabs, function(d) paste(dim(d), paste(names(d), collapse = "|"), paste(utils::head(unlist(d[1, ]), 3), collapse = "|"))))]
  helper_cols <- c("x", "y", "xend", "yend", "xmin", "xmax", "ymin", "ymax", "xintercept", "yintercept", "intercept", "slope", "label")
  tabs[!vapply(tabs, function(d) nrow(d) <= 3 && all(names(d) %in% helper_cols), logical(1))]    # drop annotation helpers (brackets, reference lines, label positions)
}

.clean_table <- function(d) {
  if ("File.Name" %in% names(d)) d$File.Name <- sub("^.*[\\\\/]", "", d$File.Name)   # raw-file name without the local folder
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  d <- d[, vapply(d, function(v) is.atomic(v) || is.factor(v), logical(1)), drop = FALSE]
  d[] <- lapply(d, function(v) if (is.factor(v)) as.character(v) else v)
  d
}

# columns: optional character vector; only these columns are kept (a table without any of them is dropped).
source_panel <- function(name, x, description = "", columns = NULL, table_names = NULL) {
  if (is.null(x)) return(invisible(NULL))
  tabs <- if (inherits(x, "ggplot")) .plot_tables(x) else if (is.data.frame(x)) list(x) else stop("source_panel(): ", name, " is neither a plot nor a data frame")
  tabs <- lapply(tabs, .clean_table)
  if (!is.null(columns)) { tabs <- lapply(tabs, function(d) d[, intersect(columns, names(d)), drop = FALSE]); tabs <- tabs[vapply(tabs, ncol, 1L) > 0] }
  if (length(tabs) == 0) { message("  Source data: ", name, " has no data frame; not stored"); return(invisible(NULL)) }
  if (!is.null(table_names) && length(table_names) == length(tabs)) names(tabs) <- table_names
  .SOURCE_STORE$panels[[name]] <- list(description = description, tables = tabs)
  invisible(NULL)
}

write_source_data <- function(file, title = "Source data") {
  P <- .SOURCE_STORE$panels
  if (length(P) == 0) { message("  Source data: nothing stored"); return(invisible(NULL)) }
  sheets <- list(); readme <- data.frame(panel = character(), sheet = character(), rows = integer(), columns = character(), description = character())
  for (nm in names(P)) {
    tabs <- P[[nm]]$tables
    for (k in seq_along(tabs)) {
      tn <- if (!is.null(names(tabs))) names(tabs)[k] else as.character(k)
      sh <- substr(if (length(tabs) == 1) nm else paste0(nm, "_", tn), 1, 31)
      sheets[[sh]] <- tabs[[k]]
      readme[nrow(readme) + 1, ] <- list(nm, sh, nrow(tabs[[k]]), paste(names(tabs[[k]]), collapse = ", "), if (k == 1) P[[nm]]$description else paste0("same panel, table: ", tn))
    }
  }
  intro <- data.frame(info = c(title,
                               "One sheet per panel: the data used to draw the panel (column names are those of the figure script). Panels with several tables have one sheet per table.",
                               "Not detected = empty cell (no imputation). Intensities are log2 unless stated otherwise."))
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "README"); openxlsx::writeData(wb, "README", intro, colNames = FALSE)
    openxlsx::writeData(wb, "README", readme, startRow = 5, headerStyle = openxlsx::createStyle(textDecoration = "bold"))
    openxlsx::setColWidths(wb, "README", cols = 1:5, widths = c(14, 26, 8, 60, 70))
    for (sh in names(sheets)) { openxlsx::addWorksheet(wb, sh); openxlsx::writeData(wb, sh, sheets[[sh]]); openxlsx::freezePane(wb, sh, firstRow = TRUE) }
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    message("  Source data written: ", file, " (", length(sheets), " sheets)")
  } else {
    dir <- sub("\\.xlsx$", "_csv", file); dir.create(dir, showWarnings = FALSE)
    utils::write.csv(readme, file.path(dir, "README.csv"), row.names = FALSE)
    for (sh in names(sheets)) utils::write.csv(sheets[[sh]], file.path(dir, paste0(sh, ".csv")), row.names = FALSE)
    message("  Source data: package openxlsx not installed; ", length(sheets), " CSV files written to ", dir)
  }
  .SOURCE_STORE$panels <- list()
  invisible(file)
}