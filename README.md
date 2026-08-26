# SimplexMCMC

This repository contains supplementary code and material for the SimplexMCMC
paper, including R and Stan code for running its three experiments.

- `experiment 1.R`–`experiment 3.R`: experiment entry points
- `Helpers/`: experiment utilities and R model setup
- `Models/`: Stan model definitions
- `rdhmc.R`: reflective dynamic HMC implementation
- `init.R`: installs dependencies and downloads fitted experiment results

## Initialize

From the repository root, run:

```sh
Rscript init.R
```

This installs the required R packages and downloads the fitted results into
`Experiments/`. Compiling Stan models also requires a working CmdStan toolchain;
follow the guidance printed by the initialization script if needed.
