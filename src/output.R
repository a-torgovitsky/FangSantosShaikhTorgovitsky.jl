library("tidyverse")
library("ggplot2")
library("tikzDevice")
library("ggpubr")

# The script also uses the numform library to remove leading zeros, see
# https://stackoverflow.com/a/42578991/509489
# Not going to load it as a library since it has name conflicts with tidyverse,
# but here's how you install it:
#
# devtools::install_github('trinker/numform')

source("theme.R")
source("latex.R")
source("load_and_clean.R")
source("plot_nearid.R")
source("plot_power.R")
source("plot_idset.R")
source("table_size.R")
source("plot_errors.R")
