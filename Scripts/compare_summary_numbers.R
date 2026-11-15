# compare_summary_numbers.R -- checks that a re-run reproduces earlier results.
# Usage (from the project folder):
#   source("Scripts/compare_summary_numbers.R")
#   compare_summary_numbers("Results_baseline", "Results")
# Matching *_summary_numbers.txt files of the two folders are compared. Lines are compared as sets (order and
# line breaks do not matter); lines that only contain script names, software versions, query dates or legend/Methods
# text are ignored. Lines found only in the old or only in the new file are printed.

IGNORE_PATTERN <- paste(c("Software versions", "query dates", "helpers\\.R|setup\\.R", "^\\[(Fig|ED Fig|Example|Methods|Supplementary)",
                          "DRAFT TEXT", "^=+$", "^=+ .* =+$", "^ED Fig 3a:", "-+$"), collapse = "|")

clean_lines <- function(x) { x <- trimws(gsub("\\s+", " ", x)); x[nzchar(x) & !grepl(IGNORE_PATTERN, x)] }

compare_summary_numbers <- function(old_dir, new_dir) {
  old_files <- list.files(old_dir, pattern = "_summary_numbers\\.txt$", recursive = TRUE, full.names = TRUE)
  new_files <- list.files(new_dir, pattern = "_summary_numbers\\.txt$", recursive = TRUE, full.names = TRUE)
  if (!length(old_files)) stop("No *_summary_numbers.txt files found in ", old_dir)
  n_diff <- 0
  for (f in old_files) {
    g <- new_files[basename(new_files) == basename(f)]
    if (length(g) != 1) { message("Not found (or not unique) in new run: ", f); n_diff <- n_diff + 1; next }
    a <- clean_lines(readLines(f, warn = FALSE)); b <- clean_lines(readLines(g, warn = FALSE))
    only_old <- setdiff(a, b); only_new <- setdiff(b, a)
    if (!length(only_old) && !length(only_new)) { message("same numbers: ", basename(f)); next }
    message("DIFFERENT: ", basename(f))
    for (x in head(only_old, 30)) message("  only in old: ", substr(x, 1, 200))
    for (x in head(only_new, 30)) message("  only in new: ", substr(x, 1, 200))
    n_diff <- n_diff + 1
  }
  message(if (n_diff == 0) "All summary files contain the same numbers." else paste(n_diff, "file(s) differ."))
  invisible(n_diff)
}
