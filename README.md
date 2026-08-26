# SimplexMCMC

This repository contains supplementary code and material for the SimplexMCMC
paper, including R and Stan code for running its three experiments.

- `experiment 1.R`–`experiment 3.R`: experiment entry points
- `ex1_utils.R`–`ex3_utils.R`: experiment-specific helpers
- `models.R` and `Models/`: R and Stan model definitions
- `download_experiments.R`: downloads the fitted experiment models

## Initialize

From the repository root, run:

```sh
Rscript init.R
```

This installs the required R packages and downloads the fitted models. Compiling
Stan models also requires a working CmdStan toolchain; follow the guidance printed
by the initialization script if it is not already installed.
