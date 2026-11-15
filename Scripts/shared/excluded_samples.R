###############################################################################
# excluded_samples.R -- samples removed from the whole study
# Sourced by shared/setup.R so that all scripts use the same lists.
###############################################################################

# Lower-depth duplicate of a neuron that was collected twice
DUPLICATE_NEURON_SAMPLES <- c("Thy1_SoC_1x_24Apr2026_B08_Cell026", "Thy1_SoC_1x_24Apr2026_F08_Cell030",
                              "Thy1_CiC_1x_24Apr2026_G02_Cell015", "Thy1_Hippo_1x_24Apr2026_B10_Cell010")
# Invalid well: two cells were collected together into this one well by mistake (LCM).
INVALID_WELL_SAMPLES <- c("Thy1_CiC_1x_24Apr2026_B04_Cell026")
EXCLUDED_SAMPLES <- c(DUPLICATE_NEURON_SAMPLES, INVALID_WELL_SAMPLES)
