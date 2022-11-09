#' ####################################################################### #
#' PROJECT: [PhD; X - DATA SIMULATIONA] 
#' CONTENTS: 
#'  - Simulation Calls for CH4
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
rm(list = ls())

# PREAMBLE ================================================================
## Directories ------------------------------------------------------------
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "Data")
if(!dir.exists(Dir.Data)){dir.create(Dir.Data)}

## Packages ---------------------------------------------------------------
### CRAN ------------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "parallel",
  "doParallel",
  "foreach"
)
sapply(package_vec, install.load.package)

# SIMULATION FUNCTION ==================================================
FUN.SimulationFramework <- function(
    seed = 42,
    ## Network Creation
    n_spec = 10,
    NetworkType = "Interaction", # or "Association"
    Sparcity = 0.2,
    MaxStrength = 0.5,
    ## Initial Individual Creation
    n_individuals = 1e3,
    n_mode = "each", # or "total"
    Env_range = c(0, 10),
    Trait_sd = 1,
    ## Carrying Capacity Creation
    k_range = c(200,200),
    ## Simulation Parameters
    d0 = 0.2,
    b0 = 0.8,
    t_max = 10,
    t_inter = 0.1,
    sd = 1,
    migration = 0.1,
    Effect_Dis = 0.5,
    verbose = TRUE
){
  set.seed(seed)
  call_info <- match.call()
  
  message("#### Package Loading")
  install.load.package <- function(x) {
    if (!require(x, character.only = TRUE))
      install.packages(x, repos='http://cran.us.r-project.org')
    require(x, character.only = TRUE)
  }
  package_vec <- c(
    "randcorr", # for random correlation matrices to be used as adjacency matrices
    "igraph" # for representing adjacency matrices as graphs
  )
  sapply(package_vec, install.load.package)
  
  message("#### Network Creation")
  if(!(NetworkType %in% c("Interaction", "Association"))){stop('NetworkType needs to be either "Interaction" or "Association"')}
  Rand_corr <- randcorr(n_spec) # establish random correlation matrix
  if(NetworkType == "Interaction"){
    Rand_corr[lower.tri(Rand_corr)] <- randcorr(n_spec)[lower.tri(randcorr(n_spec))]
    Rand_corr[sample(which(!is.na(Rand_corr)), as.integer(sum(!is.na(Rand_corr))*Sparcity))] <- 0
    Rand_corr <- graph_from_adjacency_matrix(adjmatrix = Rand_corr,
                                             mode = "directed",
                                             weighted = TRUE,
                                             diag = FALSE)
  }
  if(NetworkType == "Association"){
    Rand_corr[lower.tri(Rand_corr)] <- NA # make into undirected adjacency matrix representation
    Rand_corr[sample(which(!is.na(Rand_corr)), as.integer(sum(!is.na(Rand_corr))*Sparcity))] <- 0
    Rand_corr <- graph_from_adjacency_matrix(adjmatrix = Rand_corr,
                                             mode = "undirected",
                                             weighted = TRUE,
                                             diag = FALSE)
  }
  E(Rand_corr)$weight <-  E(Rand_corr)$weight*MaxStrength
  Network_igraph <- Rand_corr
  
  message("#### Initial Individual Creation")
  if(!(n_mode %in% c("each", "total"))){stop('n_mode needs to be either "each" or "total"')}
  ## generating IDs and species memeberships
  if(n_mode == "each"){
    ID <-  1:(n_individuals*n_spec)
    Species <- rep(paste("Sp", 1:n_spec, sep = "_"), each = n_individuals)
  }else{
    ID <- 1:n_individuals
    Species <- sample(paste("Sp", 1:n_spec, sep = "_"), size = n_individuals, replace = TRUE)
  }
  ## creating data frame to hold individuals
  ID_df <- data.frame(ID = ID, Trait = NA, X = NA, Y = NA, Species = Species)
  ## trait mean generating
  Trait_means <- runif(n = n_spec, min = Env_range[1], max = Env_range[2])
  names(Trait_means) <- paste("Sp", 1:n_spec, sep = "_")
  ## sampling trait values and assigning to data frame
  for(x in names(Trait_means)){
    ID_df$Trait[ID_df$Species == x] <- rnorm(n = length(x %in% ID_df$Species), 
                                             mean = Trait_means[x], sd = Trait_sd)
  }
  ## individual locations
  ID_df$X <- runif(n = nrow(ID_df), min = Env_range[1], max = Env_range[2])
  ID_df$Y <- runif(n = nrow(ID_df), min = Env_range[1], max = Env_range[2])
  
  message("#### Carrying Capacity Creation")
  k_vec <- as.integer(runif(n = n_spec, min = k_range[1], max = k_range[2]))
  names(k_vec) <- paste("Sp", 1:n_spec, sep = "_")
  
  message("#### Simulation Run")
  ## helper functions
  dt <- function(b0, d0, k, N, Tr, Env, sd){d0 + N * exp((Tr-Env)^2)/sd * (b0 - d0)/k}
  env.xy <- function(x = NULL){x}
  ## list object to store individuals at ech time step
  ID_ls <- list(ID_df)
  ## Network as adjacency matrix
  Effect_Mat <- as.matrix(as_adjacency_matrix(Network_igraph, attr = "weight")) # columns affect rows
  rownames(Effect_Mat) <- colnames(Effect_Mat) <- names(k_vec)
  ## setting star times
  t <- 0 # start time at 1
  names(ID_ls)[length(ID_ls)] <- t
  
  ## progress bar
  pb <- txtProgressBar(min = 0, max = t_max, style = 3)
  
  ## simulation loop over time steps
  while(t < t_max){
    if(verbose){print(t)}
    ## vectors for storing birth and death probabilities for each individual
    birth_prob <- rep(NA, nrow(ID_df))
    death_prob <- rep(NA, nrow(ID_df))
    for(i in 1:nrow(ID_df)){
      birth_prob[i] <- b0
      ### identify abundances of each species within distance box around focal individual
      Abundances_i <- rep(0, ncol(Effect_Mat))
      names(Abundances_i) <- colnames(Effect_Mat)
      Abundances_obs <- table(ID_df[-i,][ID_df[-i, "X"] < ID_df[i,"X"]+Effect_Dis &
                                           ID_df[-i, "X"] > ID_df[i,"X"]-Effect_Dis &
                                           ID_df[-i, "Y"] < ID_df[i,"Y"]+Effect_Dis &
                                           ID_df[-i, "Y"] > ID_df[i,"Y"]-Effect_Dis,
                                         "Species"])
      Abundances_i[match(names(Abundances_obs), names(Abundances_i))] <- Abundances_obs
      ### Extract all effects that the focal species is subject to in the interaction matrix
      Effects_i <- Effect_Mat[which(rownames(Effect_Mat) == ID_df[i, "Species"]),]
      ### Weigh effects by abundances
      WeightedEffects_i <- Abundances_i * Effects_i
      ### Calculate final effect
      if(sum(Abundances_i) != 0){
        FinalEffect_i <- sum(WeightedEffects_i)/sum(Abundances_i) * -1 # times -1 to make negative effects register as higher death rate, divided by sum of abundances to create weighted mean of effects
      }else{
        FinalEffect_i <- 0
      }
      # print(paste(sum(WeightedEffects_i), 
      #             sum(Abundances_i), 
      #             FinalEffect_i, sep = " - "))
      ### Assign death probability, final effect alters entire death rate
      death_prob[i] <- dt(b0 = b0, d0 = d0,
                          k = k_vec[which(names(k_vec) == ID_df$Species[i])],
                          N = sum(ID_df$Species == ID_df$Species[i]), 
                          Tr = ID_df$Trait[i], 
                          Env = env.xy(x = ID_df$X[i]), sd = sd) + 
        FinalEffect_i
      death_prob[i] <- ifelse(death_prob[i]<0, 0, death_prob[i]) # make sure probabilities are never < 0
    }
    names(birth_prob) <- names(death_prob) <- ID_df$ID
    ## event identification
    EventSample_vec <- paste(rep(c("Birth", "Death"), each = nrow(ID_df)), 
                             names(birth_prob), sep="_")
    event <- sample(
      EventSample_vec,
      size = 1,
      prob = c(birth_prob, death_prob)
    )
    ## event evaluation
    event_eval <- strsplit(event, split = "_")
    event_ID <- event_eval[[1]][2]
    event_EV <- event_eval[[1]][1]
    if(event_EV == "Birth"){
      append_df <- ID_df[ID_df$ID == event_ID, ]
      append_df$ID <- max(ID_df$ID)+1
      movement.x <- rnorm(1, 0, migration)
      movement.y <- rnorm(1, 0, migration)
      newloc.x <- append_df$X+movement.x
      newloc.y <- append_df$Y+movement.y
      ## ensuring species don't disperse beyond the environmental limit
      newloc.x <- ifelse(newloc.x<Env_range[1], Env_range[1], newloc.x)
      newloc.x <- ifelse(newloc.x>Env_range[2], Env_range[2], newloc.x)
      newloc.y <- ifelse(newloc.y<Env_range[1], Env_range[1], newloc.y)
      newloc.y <- ifelse(newloc.y>Env_range[2], Env_range[2], newloc.y)
      append_df$X <- newloc.x
      append_df$Y <- newloc.y
      ID_df <- rbind(ID_df, append_df)
    }
    if(event_EV == "Death"){
      ID_df <- ID_df[ID_df$ID != event_ID, ]
    }
    ## Gillespie time
    ### identify by how much time advances
    tadvance <- rexp(1, rate = sum(c(birth_prob, death_prob))) 
    t <- t + tadvance
    ### record data only if interval is met
    if(t - as.numeric(names(ID_ls)[length(ID_ls)]) >= t_inter){
      if(verbose){message(t)}
      ID_ls <- c(ID_ls, list(ID_df))
      names(ID_ls)[length(ID_ls)] <- t
      saveobj <- list(Call = call_info, Network = Network_igraph, 
                      Traits = Trait_means, K = k_vec, Simulation = ID_ls)
      save(saveobj, 
           file = paste0("TEMP_SIM_", seed, ".RData"))
    }
    ## update progress
    setTxtProgressBar(pb, t)
    if(nrow(ID_df) == 0){warning("All species went extinct"); break}
  }
  unlink(paste0("TEMP_SIM_", seed, ".RData"))
  return(list(Call = call_info, Network = Network_igraph, 
              Traits = Trait_means, K = k_vec, Simulation = ID_ls))
}


# SIMULATION RUNS ======================================================
## cluster
ncores <- parallel::detectCores()
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

## number of simulations to run
n_runs <- 1e2

## Associations --------------------------------------------------------
RunName <- "NDD-RIM"
print(Sys.time())
foreach(ITER = 1:n_runs) %dopar% {
  Simulation_Output <- FUN.SimulationFramework(
    seed = ITER,
    ## Network Creation
    n_spec = 20,
    NetworkType = "Interaction", # or "Association"
    Sparcity = 0.5,
    MaxStrength = 1,
    ## Initial Individual Creation
    n_individuals = 4e2,
    n_mode = "each", # or "total"
    Env_range = c(0, 10),
    Trait_sd = 1,
    ## Carrying Capacity Creation
    k_range = c(200,200),
    ## Simulation Parameters
    d0 = 0.4,
    b0 = 0.6,
    t_max = 10,
    t_inter = 0.1,
    sd = 1,
    migration = 0.2,
    Effect_Dis = 0.5,
    verbose = TRUE
  )
  save(Simulation_Output, file = file.path(Dir.Data, paste0(RunName, "_SIM_", ITER, ".RData")))
}
print(Sys.time())

