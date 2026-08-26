# Install the R packages required by the experiment scripts.
#
# Run from the repository root with:
#   Rscript install_packages.R

cran_repo <- "https://cloud.r-project.org"
stan_repo <- "https://mc-stan.org/r-packages/"

required_packages <- c(
  "bridgestan",
  "cmdstanr",
  "gridExtra",
  "jsonlite",
  "LaplacesDemon",
  "loo",
  "MASS",
  "nimble",
  "posterior",
  "scales",
  "showtext",
  "spikeslab",
  "tidyverse",
  "VGAM"
)

installed_packages <- rownames(installed.packages())
missing_packages <- setdiff(required_packages, installed_packages)

if (length(missing_packages) == 0L) {
  message("All required R packages are already installed.")
} else {
  message("Installing: ", paste(missing_packages, collapse = ", "))
  install.packages(
    missing_packages,
    repos = c(STAN = stan_repo, CRAN = cran_repo),
    dependencies = TRUE
  )
}

still_missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(still_missing) > 0L) {
  stop(
    "The following packages could not be installed: ",
    paste(still_missing, collapse = ", ")
  )
}

message("Package setup complete.")
message(
  "To compile and run Stan models, CmdStan must also be installed. ",
  "Check with cmdstanr::check_cmdstan_toolchain() and install with ",
  "cmdstanr::install_cmdstan()."
)
