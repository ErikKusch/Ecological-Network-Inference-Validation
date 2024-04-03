#' ####################################################################### #
#' PROJECT: [InfVal; Demographic Simulation of Interacting Species] 
#' CONTENTS: 
#'  - Simulation Calls across parallel cluster
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
message("Running Simulation Framework")

pb <- txtProgressBar(max = n_runs, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

foreach(ITER = 1:n_runs, .options.snow = opts) %dopar% {
  if(file.exists(file.path(Dir.Data, paste0(RunName, "_SIM_", ITER, ".RData")))){
    1+1
  }else{
    Simulation_Output <- FUN.SimulationFramework(
      seed = ITER,
      ## Network Creation
      n_spec = n_spec,
      NetworkType = NetworkType,
      Sparcity = Sparcity,
      MaxStrength = MaxStrength,
      ## Initial Individual Creation
      n_individuals = n_individuals,
      n_mode = n_mode,
      Env_range = Env_range,
      Trait_sd = Trait_sd,
      ## Carrying Capacity Creation
      k_range = k_range,
      ## Simulation Parameters
      d0 = d0,
      b0 = b0,
      env.xy = env.xy,
      t_max = t_max,
      t_inter = t_inter,
      sd = Env_sd,
      migration = migration,
      Effect_Dis = Effect_Dis,
      verbose = verbose
    )
    save(Simulation_Output, file = file.path(Dir.Data, paste0(RunName, "_SIM_", ITER, ".RData")))
  }
}
close(pb)

# SimExecCatch <- pbsapply(1:n_runs, cl = cl, FUN = function(ITER){
#   if(file.exists(file.path(Dir.Data, paste0(RunName, "_SIM_", ITER, ".RData")))){
#     1+1
#   }else{
#     Simulation_Output <- FUN.SimulationFramework(
#       seed = ITER,
#       ## Network Creation
#       n_spec = n_spec,
#       NetworkType = NetworkType,
#       Sparcity = Sparcity,
#       MaxStrength = MaxStrength,
#       ## Initial Individual Creation
#       n_individuals = n_individuals,
#       n_mode = n_mode,
#       Env_range = Env_range,
#       Trait_sd = Trait_sd,
#       ## Carrying Capacity Creation
#       k_range = k_range,
#       ## Simulation Parameters
#       d0 = d0,
#       b0 = b0,
#       env.xy = env.xy,
#       t_max = t_max,
#       t_inter = t_inter,
#       sd = Env_sd,
#       migration = migration,
#       Effect_Dis = Effect_Dis,
#       verbose = verbose
#     )
#     save(Simulation_Output, file = file.path(Dir.Data, paste0(RunName, "_SIM_", ITER, ".RData")))
#   }
# })

# ### lightening data load by retaining no intermediate simulation steps
# fs <- list.files(pattern = ".RData", Dir.Data, full.names = TRUE)
# fsize <- sapply(fs, FUN = function(x){file.size(x)})
# reduce_ls <- pbsapply(fs[fsize > 2e7], cl = cl, FUN = function(x){
#   load(x)
#   Simulation_Output$Simulation <- Simulation_Output$Simulation[c(1, 
#                                                                  (length(Simulation_Output$Simulation)-1), 
#                                                                  length(Simulation_Output$Simulation)
#                                                                  )
#                                                                ]
#   save(Simulation_Output, file = x)
# })