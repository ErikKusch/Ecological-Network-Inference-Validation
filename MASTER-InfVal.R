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
Dir.Scripts <- file.path(Dir.Base, "R Scripts")
Dirs <- c(Dir.Data, Dir.Models, Dir.Exports, Dir.Concept)
CreateDir <- sapply(Dirs, function(x) if(!dir.exists(x)) dir.create(x))
if(!dir.exists(Dir.Data)){dir.create(Dir.Data)}

## Packages ---------------------------------------------------------------
message("Loading Libraries")
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  ## Data Simulation Packages
  "parallel",
  "doParallel",
  "foreach",
  "doSNOW",
  "pbapply",
  "randcorr",
  "lubridate",
  # Inference Packages
  "Hmsc", # for HMSC models
  # "cooccur", # for COOCCUR models
  "igraph", # for graph representation
  # "devtools", # to install betalink
  # "betalink", # for computation of network dissimilarities
  ## Result Packages
  "ggplot2", # for plotting
  "tidybayes", # for plotting
  "brms", 
  # "rethinking", 
  "reshape2", 
  "cowplot", 
  "scales", 
  "ggnewscale"
)
sapply(package_vec, install.load.package)

## Functionality -------------------------------------------------------
`%nin%` <- Negate(`%in%`)

### Data Simulation Functions 
source(file.path(Dir.Scripts, "SimulationFrameworkFunctions.R"))

### Whole-Network Inference Accuracy Function
FUN_Matcomparison <- function(mat1, mat2){
  eq <- mat1==mat2 # avoid to later compute this twice
  100-round(sum(!eq, na.rm = TRUE)/sum(!is.na(eq))*100, 2) # get the percentage of non equal values
}

# RUN SET-UP ===========================================================
## Simulation Parameters -----------------------------------------------
n_runs <- 1e3 # number of networks to simulate and infer for
## Network Creation
n_spec = 30 # number of species per network
NetworkType = "Association" # type of links
Sparcity = 0 # how many % of associations should be exactly 0
MaxStrength = 20 # absolute maximum of interspecific links
## Initial Individual Creation
n_individuals = 5e2 # number of individuals for initialisation
n_mode = "each" # how to interpret the above number
Env_range = c(0, 10) # environmental landscape range
Trait_sd = 1 # standard deviation of traits per species
## Carrying Capacity Creation
k_range = c(300,300)
## Simulation Parameters
d0 = 0.4 # base death rate
b0 = 0.6 # base birth rate
env.xy = function(x = NULL, y = NULL){x} #environmental maladpation function
t_max = 20 # simulation time
t_inter = 1 # when to record data
Env_sd = 5 # environmental maladaption SD, higher = more permissive environment
migration = 0.5 # sd of 0 centred normal for relocation of offspring
Effect_Dis = 1 # distance at which link effect manifests
verbose = FALSE # whether to produce simulation progress tracker in console
## File Naming
RunName <- "Revision"

## Inference Settings --------------------------------------------------
### HMSC MCMC Settings
nSamples <- 7e3
thin <- 1
nWarmup <- round(nSamples*0.3*thin, 0)
nChains <- 4
### Gridding of data along each axis:
n_Grid <- 5

## Cluster for parallel computation ------------------------------------
message("Registering Clusters")
ncores <- 100 # ifelse(parallel::detectCores() > n_runs, n_runs, parallel::detectCores())
cl <- parallel::makeCluster(ncores
                            # , outfile = "Log.txt"
                            )
parallel::clusterExport(cl, varlist = ls(), envir = environment())
doSNOW::registerDoSNOW(cl)

# METADATA FOR RUN WRITING =============================================
MetaF <- file.path(Dir.Data, paste0("META-", RunName, ".RData"))
if(file.exists(MetaF)){stop("A run with this name has already been executed.")}

# DATA SIMULATION ======================================================
source(file.path(Dir.Scripts, "DataSimulations.R"))
save.image(file = MetaF)

# ASSOCIATION INFERENCE ================================================
source(file.path(Dir.Scripts, "Inference.R"))

# POST-INFERENCE ANALYSES ==============================================
source(file.path(Dir.Scripts, "Results.R"))
