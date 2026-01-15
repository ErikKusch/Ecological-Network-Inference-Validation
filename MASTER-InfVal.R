#' ####################################################################### #
#' PROJECT: [InfVal; Master Script for Parametrising Entire Analyses]
#' CONTENTS:
#'  - Master File for Framework Execution
#'  DEPENDENCIES:
#'  -
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
rm(list = ls())
set.seed(42)

# PREAMBLE ================================================================
## Directories ------------------------------------------------------------
message("Registering Directories")
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "Data")
Dir.Models <- file.path(Dir.Base, "Models")
Dir.Exports <- file.path(Dir.Base, "Exports")
Dir.Concept <- file.path(Dir.Base, "Concept")
Dir.Scripts <- file.path(Dir.Base, "RScripts")
Dirs <- c(Dir.Data, Dir.Models, Dir.Exports, Dir.Concept)
CreateDir <- sapply(Dirs, function(x) if (!dir.exists(x)) dir.create(x))
if (!dir.exists(Dir.Data)) {
  dir.create(Dir.Data)
}

## Packages ---------------------------------------------------------------
message("Loading Libraries")
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}
package_vec <- c(
  ## Data Simulation Packages
  "parallel",
  "doParallel",
  "pbapply",
  # Inference Packages
  "netassoc", # for netassoc models
  "Hmsc", # for HMSC models
  "igraph", # for graph representation
  "devtools", # to install cooccur
  ## Result Packages
  "ggplot2", # for plotting
  # "tidybayes", # for plotting
  # "brms",
  # "rethinking",
  # "reshape2",
  "cowplot"
  # "scales",
  # "ggnewscale",
  # "ggpubr"
)
sapply(package_vec, install.load.package)

## NetSimVal Setup --------------------------------------------------------
if ("NetSimVal" %in% rownames(installed.packages()) == FALSE) {
  devtools::install_github("https://github.com/ErikKusch/NetSimVal")
}
library(NetSimVal)
package_vec <- c(package_vec, "NetSimVal")

## COOCCUR Setup ----------------------------------------------------------
if ("cooccur" %in% rownames(installed.packages()) == FALSE) {
  devtools::install_github("https://github.com/griffithdan/cooccur")
}
library(cooccur)
package_vec <- c(package_vec, "cooccur")

## Functionality -------------------------------------------------------
`%nin%` <- Negate(`%in%`)
source(file.path(Dir.Scripts, "Paper Plotting Functions.r"))

# RUN SET-UP ===========================================================
## Simulation Parameters -----------------------------------------------
# package loading
if ("NetSimVal" %in% rownames(installed.packages()) == FALSE) {
  devtools::install_github("https://github.com/ErikKusch/NetSimVal")
}
library(NetSimVal)

## global simulation parameters -----
n_runs <- 1e3 # number of networks to simulate and infer for

## Network Creation
n_spec <- 20 # number of species per network
NetworkType <- "Association" # type of links
Sparcity <- 0.1 # how many % of associations should be exactly 0
MaxStrength <- 1 # absolute maximum of interspecific links

## Initial Individual Creation
n_individuals <- 1e3 # number of individuals for initialisation
n_mode <- "each" # how to interpret the above number
Env_range <- c(0, 10) # environmental landscape range
Trait_sd <- 0.5 # standard deviation of traits per species

## Space Creation
Env_col <- Env_row <- 1e3
x_gradient <- function(x) x # function relating x-location to environmental value
y_gradient <- function(y) 0 # function relating y-location to environmental value

## Carrying Capacity Creation
k_range <- c(300, 300) # lower and upper limit of carrying capacities

## Simulation Parameters
d0 <- 0.4 # base death rate
b0 <- 0.6 # base birth rate
t_max <- 50 # simulation time
t_inter <- 1 # when to record data
Env_sd <- 1 # environmental maladaption SD, higher = more permissive environment
migration <- 0.5 # sd of 0 centred normal for relocation of offspring
migration_top <- 0.05 # point at which migration peaks
Effect_Dis <- 0.25 # distance at which link effect manifests
verbose <- TRUE # whether to produce simulation progress tracker in console

## run names
RunNames <- list("NoSpaceBinaryInterac", "BinaryInterac", "ContinuousInterac", "RareSpecies")

## Inference Settings --------------------------------------------------
### HMSC MCMC Settings
nSamples <- 7e3
thin <- 1
nWarmup <- round(nSamples * 0.3 * thin, 0)
nChains <- 4
### Gridding of data along each axis:
n_Grid <- 20

## Cluster for parallel computation ------------------------------------
# message("Registering Clusters")
# DesiredCores <- ifelse(n_runs > parallel::detectCores(), 2, 50) # parallel::detectCores()
# ncores <- ifelse(parallel::detectCores() > DesiredCores, DesiredCores, parallel::detectCores())
# cl <- parallel::makeCluster(
#   ncores
#   # , outfile = "Log.txt"
# )
# parallel::clusterExport(cl, varlist = ls(), envir = environment())
# clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))
# doSNOW::registerDoSNOW(cl)

# LOOPING OVER RUNS ======================================================
message("Actual Simulations")
# lapply(RunNames, function(RunName) {
for (RunName in RunNames) {
  # RunName <- RunNames[[1]]
  print(RunName)
  # parallel::clusterExport(cl, varlist = "RunName", envir = environment())
  ## METADATA FOR RUN WRITING ============================================
  MetaF <- file.path(Dir.Data, paste0("META-", RunName, ".RData"))
  if (file.exists(MetaF)) {
    stop("A run with this name has already been executed.")
  }

  ## DATA SIMULATION =====================================================
  source(file.path(Dir.Scripts, "DataSimulations.R"))

  ## ASSOCIATION INFERENCE ===============================================
  source(file.path(Dir.Scripts, "Inference.R"))

  ## POST-INFERENCE ANALYSES =============================================
  # source(file.path(Dir.Scripts, "Results.R"))
  # save.image(file = MetaF)
}
# }
# )

# NON-RESULT FIGURE CREATION ==============================================
source(file.path(Dir.Scripts, "Figure 1 - Demo of NetSimVal Capabilities.r"))
source(file.path(Dir.Scripts, "Figure 2 - Demo of SimulationRuns.r"))
