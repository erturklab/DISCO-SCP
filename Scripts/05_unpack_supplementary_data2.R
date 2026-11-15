###############################################################################
# 05_unpack_supplementary_data2.R
# For readers who start from "Supplementary Data 2" (OPTION B in its README sheet) instead of running
# 00_build_master_dataset.R. Writes the sheets back into the files the scripts expect:
#   Results/master_dataset/3D_SCP_all_Data_Master_file.csv, master_sample_metadata.csv,
#   contamination_analysis.csv, pen_sc_systematic_screen_ALL_PROTEINS.csv   and   Data/Morphological_analysis.xlsx
#
# How to use (RStudio): set the R working directory to the project folder, put Supplementary_Data_2.xlsx into
# Data/ (or change IN_FILE), then source() this file. Afterwards run the figure scripts 01-04.
# Packages: readr, openxlsx, readxl (installed by install_packages.R).
###############################################################################

suppressPackageStartupMessages({ library(readr); library(openxlsx); library(readxl) })

PROJECT_ROOT <- getOption("scp.project_root", default = normalizePath(getwd(), winslash = "/"))
IN_FILE     <- file.path(PROJECT_ROOT, "Data", "Supplementary_Data_2.xlsx")
MASTER_DIR  <- file.path(PROJECT_ROOT, "Results", "master_dataset")
MORPH_SHEET <- "morphological_features_Gen3_Thy"
MORPH_OUT   <- file.path(PROJECT_ROOT, "Data", "Morphological_analysis.xlsx")

CSV_SHEETS <- c(                                   # same sheet names as in Supplementary Data 2
  master_protein_table   = "3D_SCP_all_Data_Master_file.csv",
  sample_metadata        = "master_sample_metadata.csv",
  contaminant_annotation = "contamination_analysis.csv",
  pen_control_screen     = "pen_sc_systematic_screen_ALL_PROTEINS.csv"
)

if (!file.exists(IN_FILE)) stop("File not found: ", IN_FILE, "\nPut Supplementary_Data_2.xlsx into Data/.", call. = FALSE)
sheets <- excel_sheets(IN_FILE)
absent <- setdiff(c(names(CSV_SHEETS), MORPH_SHEET), sheets)
if (length(absent) > 0) stop("Sheets missing in the workbook: ", paste(absent, collapse = ", "), call. = FALSE)

dir.create(MASTER_DIR, showWarnings = FALSE, recursive = TRUE)
for (sh in names(CSV_SHEETS)) {
  x <- read_excel(IN_FILE, sheet = sh, guess_max = 100000, .name_repair = "minimal")
  write_csv(x, file.path(MASTER_DIR, CSV_SHEETS[[sh]]))
  message("  written: Results/master_dataset/", CSV_SHEETS[[sh]], "  (", nrow(x), " rows)")
}

morph <- as.data.frame(read_excel(IN_FILE, sheet = MORPH_SHEET, .name_repair = "minimal"))
wb <- createWorkbook(); addWorksheet(wb, MORPH_SHEET); writeData(wb, MORPH_SHEET, morph)
dir.create(dirname(MORPH_OUT), showWarnings = FALSE, recursive = TRUE)
saveWorkbook(wb, MORPH_OUT, overwrite = TRUE)
message("  written: Data/Morphological_analysis.xlsx  (", nrow(morph), " rows)")
message("Done. Next: run the figure scripts 01-04.")
