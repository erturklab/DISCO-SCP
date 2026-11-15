###############################################################################
# install_packages.R -- run ONCE, from the project folder:  source("install_packages.R")
# Installs the R packages needed by the analysis scripts that are not installed yet.
# The analysis scripts themselves never install anything.
#
# The package versions of every run are saved in Results/figures/<Figure>/Figure*_software_versions.txt.
###############################################################################
cran <- c("tidyverse", "readxl", "ggrepel", "patchwork", "colorspace", "ggdist", "scales", "arrow", "ComplexUpset", "msigdbr", "httr",
          "data.table", "openxlsx", "glue", "svglite", "ggrastr", "MASS")
bioc <- c("fgsea", "SummarizedExperiment", "QFeatures", "msqrob2", "proDA", "limma", "BiocParallel", "MsCoreUtils", "S4Vectors")

installed <- rownames(installed.packages())
cran_todo <- setdiff(cran, installed)
if (length(cran_todo)) install.packages(cran_todo)

bioc_todo <- setdiff(bioc, installed)
if (length(bioc_todo)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(bioc_todo, update = FALSE, ask = FALSE)
}
still_missing <- setdiff(c(cran, bioc), rownames(installed.packages()))
if (length(still_missing)) stop("Could not install: ", paste(still_missing, collapse = ", ")) else message("All packages installed.")
