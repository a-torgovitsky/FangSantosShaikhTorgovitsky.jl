# FangSantosShaikhTorgovitsky

Repo: https://github.com/a-torgovitsky/FangSantosShaikhTorgovitsky.jl

This package contains the code for the Monte Carlo simulations in "Inference for Large-Scale Linear Systems with Known Coefficients" by Fang, Santos, Shaikh, and Torgovitsky (2021, working paper).
A shortened version of the working paper (not including these simulations) is forthcoming in _Econometrica_.

The code is written in Julia, with plots and tables written in R.
If you want to apply the methodology, the [`lpinfer` R package](https://github.com/conroylau/lpinfer) is a more user-friendly way to do so.
This repository is just for replication purposes.

## Running the code

To run the code in serial:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate() # download the correct package versions
using FangSantosShaikhTorgovitsky
menu()
```

You should see this:
```julia
================================================================================
                        ███████ ███████ ███████ ████████
                        ██      ██      ██         ██
                        █████   ███████ ███████    ██
                        ██           ██      ██    ██
                        ██      ███████ ███████    ██
================================================================================
Inference for Large-Scale Systems of Linear Inequalities
Fang, Santos, Shaikh, and Torgovitsky (2021)
================================================================================
Which simulation would you like to run?
         1. Identified sets (Figure 3)
         2. Comparison with Wald in nearly point identified cases (Figure 4)
         3. Size of FSST in many designs with λᵇ (Table 1)
         4. Size of FSST in many designs with λʳ (Table 2)
         5. Power (Figure 5)
Enter choice:
```
Select from the options given.

Option 1 (identified sets) takes under a minute to run.

The others options are much more time-consuming and should be run in parallel.
To do this, start Julia with:

```sh
julia -p 40 # you likely want a different number
```

Then run the following commands:

```julia
using Distributed
using Pkg
Pkg.activate(".")
Pkg.instantiate() # download correct versions only on one worker
@everywhere using Pkg
@everywhere Pkg.activate(".")
# Avoid race condition by loading on a single worker first
using FangSantosShaikhTorgovitsky
@everywhere using FangSantosShaikhTorgovitsky
menu()
```

Then select from the options given as before.

## Plots and tables

Plots and tables were coded in R.
They are called from Julia automatically when using `menu()`.

## Packages, versions and requirements

- Julia 1.1.0
    - Instantiating as above should install the correct Julia package versions.
- `R` packages used for figures and tables:
  - `tidyverse` 1.3.1
  - `ggplot2` 3.3.5
  - `tikzDevice` 0.12.3.1
  - `ggpubr` 0.4.0
  - `numform` 0.6.4
    - This package is not on CRAN but can be installed via `devtools::install_github('trinker/numform')`
- Gurobi 9.1.2
