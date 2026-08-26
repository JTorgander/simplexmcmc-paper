# Initialize the repository by installing its R dependencies and downloading
# the fitted experiment models.
#
# Run from the repository root with:
#   Rscript init.R

source("Helpers/install_packages.R")
source("Helpers/download_experiments.R")

message("Repository initialization complete.")
