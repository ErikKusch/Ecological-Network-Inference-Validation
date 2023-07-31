#' ####################################################################### #
#' PROJECT: [PhD; X - DATA SIMULATIONS] 
#' CONTENTS: 
#'  - Simulation Calls for CH4
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
rm(list = ls())

# PREAMBLE ================================================================
## Directories ------------------------------------------------------------
message("Registering Directories")
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "Data")
if(!dir.exists(Dir.Data)){dir.create(Dir.Data)}

## Packages ---------------------------------------------------------------
message("Loading Libraries")
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "parallel",
  "doParallel",
  "foreach",
  "doSNOW"
)
sapply(package_vec, install.load.package)

# SIMULATION FUNCTIONS =================================================
message("Registering Simulation Framework")
source("SimulationFrameworkFunctions.R")

# SIMULATION RUNS ======================================================
## number of simulations to run
n_runs <- 5e2

## cluster
message("Registering Clusters")
ncores <- ifelse(parallel::detectCores() > 50, 50, parallel::detectCores())
cl <- parallel::makeCluster(ncores)
doSNOW::registerDoSNOW(cl)
# doParallel::registerDoParallel(cl)

## running simulations
message("Running Simulation Framework")
for(Env_sd in c(1, 0.75, 2.5, 5, 10)){
  print(paste("Data simulation env_sd =", Env_sd))
  RunName <- paste("Association", Env_sd, sep = "_")
  pb <- txtProgressBar(max = n_runs, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  parallel::clusterExport(cl, varlist = c("Env_sd", "RunName"), envir = environment())
  foreach(ITER = 1:n_runs, .options.snow = opts) %dopar% {
    if(file.exists(file.path(Dir.Data, paste0(RunName, "_SIM_", ITER, ".RData")))){
      1+1
    }else{
      Simulation_Output <- FUN.SimulationFramework(
        seed = ITER,
        ## Network Creation
        n_spec = 20,
        NetworkType = "Association", # or "Association"
        Sparcity = 0.5,
        MaxStrength = 1,
        ## Initial Individual Creation
        n_individuals = 5e2,
        n_mode = "each", # or "total"
        Env_range = c(0, 10),
        Trait_sd = 1,
        ## Carrying Capacity Creation
        k_range = c(200,200),
        ## Simulation Parameters
        d0 = 0.4,
        b0 = 0.6,
        env.xy = function(x = NULL, y = NULL){x},
        t_max = 20,
        t_inter = 0.1,
        sd = 1,
        migration = 0.2,
        Effect_Dis = 0.5,
        verbose = FALSE
      )
      save(Simulation_Output, file = file.path(Dir.Data, paste0(RunName, "_SIM_", ITER, ".RData")))
    }
  }
  close(pb)
}

## closing cluster
stopCluster(cl) 
