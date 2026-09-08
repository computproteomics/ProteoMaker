# ProteoMaker: Simulating proteoforms in bottom-up proteomics

[![bio.tools](https://img.shields.io/badge/bio.tools-proteomaker-blue.svg?labelColor=gray&logo=data:image/svg%2bxml;base64,PHN2ZyBpZD0idXVpZC02ZTIxYTIzOC04NWFmLTRlNDctYjI1OC05ZTEyZDg2MzJmYmUiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgdmlld0JveD0iMCAwIDQzMi44NCA0MzIuODQiPjxwYXRoIGQ9Ik01NS4zLDQyNy4yN2wtNDkuNzMtNDkuNzNjLTcuNDMtNy40My03LjQzLTE5LjQ3LDAtMjYuODlsMTMxLjYzLTEzMS42M2M3LjQzLTcuNDMsMTkuNDctNy40MywyNi44OSwwbDQ5LjczLDQ5LjczYzcuNDMsNy40Myw3LjQzLDE5LjQ3LDAsMjYuODlsLTEzMS42MywxMzEuNjNjLTcuNDMsNy40Mi0xOS40Niw3LjQyLTI2Ljg5LDBaIiBmaWxsPSIjZmZmIi8+PHBhdGggZD0iTTIyNy43MSwyNTUuNzFsLTcuMTgsNy4zNy02LjM0LDYuNC01MS4xMy01MC4xOSw2LjM0LTYuNCw3LjE4LTcuMzdjNi42NC02Ljc2LDE3LjQ1LTYuODYsMjQuMjEtLjIybDI2LjY5LDI2LjJjNi43Nyw2LjY0LDYuODcsMTcuNDUuMjMsMjQuMjFoMFoiIGZpbGw9IiNmZmYiLz48cGF0aCBkPSJNNDMwLjQsMTguNTNsLTE2LjExLTE2LjEiIGZpbGw9IiNmZmYiLz48cGF0aCBkPSJNMzQ4LjQ2LDYzLjM3bC0xMTkuNzQsMTE5LjczLTI0Ljk0LDI0Ljk0LDIxLjAxLDIxLjAxLDI0Ljk0LTI0Ljk0LDExOS43NC0xMTkuNzMiIGZpbGw9IiNmZmYiLz48cGF0aCBkPSJNMzY5LjQ0LDg0LjM4bDI3LjE2LDEuOSwzNi4yNC02NS4zTDQxMS44NSwwbC02NS4zLDM2LjIzLDEuOSwyNy4xNyIgZmlsbD0iI2ZmZiIvPjxwYXRoIGQ9Ik0xMTIuNjEsMTU1LjZjLTI5LjE0LDExLjgzLTYzLjgsNS45NC04Ny40NC0xNy43QzIuOTYsMTE1LjY5LTMuNTgsODMuNzcsNS41NCw1NS44MyIgZmlsbD0iI2ZmZiIvPjxwYXRoIGQ9Ik01Ny4xMiw0LjI2YzI3LjkxLTkuMTUsNTkuODctMi41OCw4Mi4wNywxOS42MywyMy42MSwyMy42MSwyOS41Myw1OC4yNSwxNy43LDg3LjM5IiBmaWxsPSIjZmZmIi8+PHBhdGggZD0iTTI3Ny4yMywzMjAuMjJjLTExLjgzLDI5LjE0LTUuOTQsNjMuOCwxNy43LDg3LjQ0LDIyLjIxLDIyLjIxLDU0LjEzLDI4Ljc1LDgyLjA3LDE5LjYzIiBmaWxsPSIjZmZmIi8+PHBhdGggZD0iTTMyMS41NiwyNzUuOTRjMjkuMTQtMTEuODMsNjMuNzgtNS45Miw4Ny4zOSwxNy43LDIyLjIxLDIyLjIxLDI4Ljc4LDU0LjE2LDE5LjYzLDgyLjA3IiBmaWxsPSIjZmZmIi8+PHBhdGggZD0iTTE2My45MSwyMDcuNTRMNDIuNyw5MS4yNmw1MC43MS00OC4xMiwxMzAuMjMsMTM1LjU3LTkuMjIsOS4xOC0xMy4wOSwxMi4yOC01Ljk4LTQuMTloMGMtLjI4LS4xOC00LjMtMi4yMy00LjYtMi4zNi0uMzctLjE2LTEuOTUtMS4xNi01LjctLjg0LTMuMTEuMjctNS43NCwxLjYyLTguNTcsMy45OGwtMTIuNTgsMTAuNzhoLjAxWiIgZmlsbD0iI2ZmZiIvPjxwYXRoIGQ9Ik0yMjQuMTgsMjY4LjI0bDExNS4zMywxMjAuNDQsNDguMzMtNTAuODktMTM0LjUyLTEyOS4zNC05LjIyLDkuMjYtMTIuMzQsMTMuMTMsNC4xNSw1Ljk1aDBjLjE3LjI4LDIuMiw0LjI4LDIuMzMsNC41OC4xNi4zNywxLjE1LDEuOTQuOCw1LjY4LS4yOCwzLjEtMS42NSw1Ljc0LTQuMDMsOC41OGwtMTAuODQsMTIuNjJoLjAxWiIgZmlsbD0iI2ZmZiIvPjwvc3ZnPg==)](https://bio.tools/proteomaker)
![Tool Type](https://img.shields.io/badge/tool%20type-Library%20%7C%20Web%20application-blue.svg?labelColor=gray)
[![Bridge](https://img.shields.io/badge/bridge-bio.tools%20%E2%86%92%20github-blue.svg?labelColor=orange&logo=data:image/svg%2bxml;base64,PHN2ZyBpZD0idXVpZC0zYTVmNzA5MS0xMzM2LTQxNGUtYjBmNS1jMDdkMmQwYWM2MDEiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgdmlld0JveD0iMCAwIDEwLjAxIDYuOCI+PHBhdGggZD0iTTUuNjYsNS42NmwtLjc0Ljc0Yy0uMjYuMjYtLjU5LjQtLjk1LjRzLS43LS4xNC0uOTYtLjRMLjUyLDMuOTFjLS4yMy0uMjMtLjM3LS41My0uMzktLjg1TDAsMS40NWMtLjAzLS4zNy4xLS43NC4zNi0xLjAxQy42NC4xNCwxLjA0LS4wMiwxLjQ1LjAxbDEuNjEuMTJjLjMyLjAzLjYyLjE3Ljg0LjM5bC42Mi42MWMuMTYuMTcuMTYuNDMsMCwuNTktLjE3LjE2LS40My4xNi0uNTksMGwtLjYxLS42MWMtLjA5LS4wOS0uMi0uMTQtLjMzLS4xNWwtMS42LS4xM2MtLjE2LS4wMS0uMzIuMDUtLjQyLjE3LS4xLjEtLjE1LjI1LS4xNC4zOWwuMTMsMS42MWMuMDEuMTIuMDYuMjQuMTUuMzJsMi40OSwyLjVjLjEuMDkuMjMuMTUuMzcuMTUuMTMsMCwuMjYtLjA2LjM2LS4xNWwuNzQtLjc0Yy4wMi4xNC4wOS4yOC4yLjM4LjEuMTEuMjQuMTguMzkuMloiIGZpbGw9IiNmZmYiLz48cGF0aCBkPSJNMTAsMS40NWwtLjEzLDEuNjFjLS4wMi4zMi0uMTYuNjItLjM5Ljg1bC0yLjQ5LDIuNDljLS4yNi4yNi0uNi40LS45Ni40LS4xMiwwLS4yNS0uMDItLjM3LS4wNS0uMjItLjA3LS4zNC0uMy0uMjgtLjUyLjA2LS4yMi4yOS0uMzQuNTEtLjI4LjA1LjAxLjEuMDIuMTQuMDIuMTQsMCwuMjctLjA2LjM3LS4xNWwyLjQ5LTIuNWMuMDktLjA4LjE1LS4yLjE1LS4zMmwuMTMtMS42MWMuMDEtLjE0LS4wNC0uMjktLjE0LS4zOS0uMS0uMTItLjI2LS4xOC0uNDItLjE3bC0xLjYuMTNjLS4xMy4wMS0uMjQuMDYtLjMzLjE1bC0xLjc3LDEuNzdjLS4wMy0uMTQtLjEtLjI3LS4yMS0uMzgtLjEtLjExLS4yNC0uMTgtLjM4LS4ybDEuNzgtMS43OGMuMjItLjIyLjUyLS4zNi44NC0uMzlsMS42MS0uMTJzLjA3LS4wMS4xLS4wMWMuMzgsMCwuNzQuMTYuOTkuNDQuMjYuMjcuMzkuNjQuMzYsMS4wMVoiIGZpbGw9IiNmZmYiLz48Y2lyY2xlIGN4PSI4LjA1IiBjeT0iMS45NiIgcj0iLjcxIiBmaWxsPSIjZmZmIi8+PGNpcmNsZSBjeD0iMS45NiIgY3k9IjEuOTYiIHI9Ii43MSIgZmlsbD0iI2ZmZiIvPjxwYXRoIGQ9Ik00LjUyLDUuMjVjLS4xMi4xMi0uMzEuMTUtLjQ2LjA5LS4wNS0uMDItLjA5LS4wNS0uMTMtLjA5bC0uMzMtLjMzYy0uMjUtLjI1LS4zOS0uNTktLjM5LS45NXMuMTQtLjcuMzktLjk1bC4zMS0uMzFjLjE2LS4xNi40Mi0uMTYuNTgsMCwuMTUuMTUuMTcuMzkuMDMuNTZsLS4zMy4zM2MtLjEuMS0uMTYuMjMtLjE2LjM3cy4wNi4yNy4xNi4zN2wuMzMuMzNjLjE2LjE3LjE2LjQyLDAsLjU4WiIgZmlsbD0iI2ZmZiIvPjxwYXRoIGQ9Ik02Ljc5LDMuOTdjMCwuMzYtLjE0LjctLjM5Ljk1bC0uMzMuMzNjLS4xNy4xNi0uNDMuMTYtLjU5LDAtLjE1LS4xNS0uMTYtLjM5LS4wMy0uNTVsLjM2LS4zNmMuMS0uMS4xNi0uMjMuMTYtLjM3cy0uMDYtLjI3LS4xNi0uMzdsLS4zLS4zYy0uMTYtLjE2LS4xNi0uNDIsMC0uNTkuMDgtLjA4LjE5LS4xMi4yOS0uMTIuMDcsMCwuMTQuMDIuMi4wNXMuMTEuMDguMTYuMTNjLjA4LjA4LjE1LjE1LjI0LjI0LjI1LjI1LjM5LjU5LjM5Ljk1WiIgZmlsbD0iI2ZmZiIvPjwvc3ZnPg==)](https://bio-tools.github.io/biohackathon2025/)

[![DOI](https://zenodo.org/badge/233879107.svg)](https://doi.org/10.5281/zenodo.22641240)

ProteoMaker is a platform for the  generation of an in-silico bottom-up proteomics data set with a ground truth on the level of proteoforms. 

All the parameters that are used to generate the data are described in [vignettes/Parameters.qmd](vignettes/Parameters.qmd) (render to HTML if you prefer). The script that runs the entire pipeline is [`RunSims.R`](inst/cmd/RunSims.R). Alternatively, you can use the [Vignette](vignettes/Vignette.html). The simulations with multiple parameters can be set up and run. ProteoMaker also provides comparison of the results with the ground truth using [benchmarking metrics](vignettes/Benchmarks.qmd) and visual comparison between simulated data sets.

You can also use the Shiny app: https://computproteomics.bmb.sdu.dk/app_direct/ProteoMaker/

The pipeline can be described as follows:

0) General functions to run the simulations: [`00_BatchRunFuncs.R`](R/00_BatchRunFuncs.R)
1) Generation of ground truth data at the proteoform level [`01_GenerateGroundTruth.R`](R/01_GenerateGroundTruth.R).
2) Digestion of the proteoforms from the ground truth: [`02_Digestion.R`](R/02_Digestion.R).
3) In silico MS run: [`03_MSRun.R`](R/03_MSRun.R).
4) Functions for data analysis from the peptide to proteins: [`04_DataAnalysis.R`](R/04_DataAnalysis.R).
5) Statistical testing: [`05_Statistics.R`](R/05_Statistics.R).
6) Benchmarking: [`06_Benchmarks.R`](R/06_Benchmarks.R).

## Installation

Install the package from GitHub with `pak`:
```
install.packages("pak")
pak::pak("computproteomics/ProteoMaker")
```

## Quick Start

1. Load the package in R:
   ```r
   library(ProteoMaker)
   ```
2. Run the default simulation (adjust the output folder as needed):
   ```r
   Param  <- def_param()
   Config <- set_proteomaker(resultFilePath = "results")
   Result <- run_sims(Param, Config)
   ```
   Intermediate `outputDataAnalysis_<hash>.RData` files and benchmark tables are written to `results/`.
3. Explore the outputs, for example:
   ```r
   #Explore output
   Benchmatrix <- matrix_benchmarks(Result, Config)
   visualize_benchmarks(Benchmatrix)
   ```
   or open the vignette for a full walkthrough.

## Repository Overview

| Path | Description |
| --- | --- |
| R/ | Core simulation, analysis, and benchmarking functions |
| inst/config/parameters.yaml | Default parameter definitions consumed by `def_param()` |
| inst/cmd/RunSims.R | Convenience script to run the full simulation pipeline |
| vignettes/ | Walkthroughs and usage examples |
| inst/img/ | Diagrams and figures (e.g. pipeline layout) |
| tests/ | Automated tests for key functionality |
| inst/shiny/ | Shiny interface for interactive configuration (if used) |

### Running full batches and benchmarking

Running the [Vignette](vignettes/Vignette.html) allows running full batches without having to re-run the data sets which have been built with the same set of parameters. In addition, the pipeline is run hierarchically to avoid repetitive execution of identical down-stream analysis. This is done via creating hashes of the parameter configurations and writing intermediate and final results into respective tables.

Re-running the full batch with different assessment of the benchmarking metrics will avoid re-running the data set generation and analysis, and thus should be superfast.

Important remarks:

- You need to always define _all_ parameters in the beginning of this script
- Be aware that changing downstream parameters (ground truth, digestion) can immensely increase the number of possible parameter settings
- Keep always the result files in the respective folder (`resultFilePath`) if you didn't change anything in the pipeline such as any of the methods in the sourced files. This will allow you to run the full batch without re-running the data set generation and analysis.
- Benchmarking results are skipped for data sets with fewer than 100 quantified proteins.
