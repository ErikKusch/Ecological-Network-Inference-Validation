#' ####################################################################### #
#' PROJECT: [PhD; 4 - DATA SIMULATIONS] 
#' CONTENTS: 
#'  - Data Simulation Framework Functions
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
#' These are the functions which make up our simulation framework.

# RANDOM NETWORK ===============================================================
#' @param n_spec Numeric. Number of species / size of node set.
#' @param NetworkType Character. Whether to create an interaction or association matrix. In an interaction matrix, columns act on rows. In an association matrix, entries are mirrored around the diagonal. Possible values: "Association" or "Interaction".
#' @param Sparcity Numeric. Proportion of unique associations or interactions which are exactly 0.
#' @param MaxStrength Numeric. Maximum value of associations or interactions.
#' @param seed Numeric. Seed for random processes.
#' @return An igraph object with association/interaction strength stored as "weight" attribute of edges.
FUN.RandNetwork <- function(n_spec = 20,
                            NetworkType = "Association", # or "Association"
                            Sparcity = 0.5,
                            MaxStrength = 1,
                            seed = 42){
  library(randcorr)
  library(igraph)
  
  set.seed(seed)
  
  if(!(NetworkType %in% c("Interaction", "Association"))){stop('NetworkType needs to be either "Interaction" or "Association"')}
  
  Rand_corr <- randcorr::randcorr(n_spec) # establish random correlation matrix
  if(NetworkType == "Interaction"){
    Rand_corr[lower.tri(Rand_corr)] <- randcorr::randcorr(n_spec)[lower.tri(randcorr::randcorr(n_spec))]
    Rand_corr[sample(which(!is.na(Rand_corr)), as.integer(sum(!is.na(Rand_corr))*Sparcity))] <- 0
    Rand_corr <- igraph::graph_from_adjacency_matrix(adjmatrix = Rand_corr,
                                                     mode = "directed",
                                                     weighted = TRUE,
                                                     diag = FALSE)
  }
  if(NetworkType == "Association"){
    Rand_corr[lower.tri(Rand_corr)] <- NA # make into undirected adjacency matrix representation
    diag(Rand_corr) <- NA
    Rand_corr[sample(which(!is.na(Rand_corr)), as.integer(sum(!is.na(Rand_corr))*Sparcity))] <- 0
    Rand_corr <- igraph::graph_from_adjacency_matrix(adjmatrix = Rand_corr,
                                                     mode = "undirected",
                                                     weighted = TRUE,
                                                     diag = FALSE)
  }
  E(Rand_corr)$weight <-  E(Rand_corr)$weight*MaxStrength
  return(Rand_corr)
}

# RANDOM INITIALISING INDIVIDUALS ==============================================
#' @param n_spec Numeric. Number of species / size of node set.
#' @param n_individuals Numeric. How many initialising individuals to randomly create.
#' @param n_mode Character. How to use n_individuals. If "each" the the number of initialising individuals will be equal to n_individual multiplied by n_spec. If "total", n_individuals are created with random species memberships within n_spec.
#' @param Env_range Numeric vector of length 2. Minimum and maximum values of environmental gradient.
#' @param Trait_sd Numeric. Standard deviation of trait distributions from which trait values of each individual are drawn.
#' @param seed Numeric. Seed for random processes.
#' @return A list of length 2 containing (1) ID_df - the data frame of initialising individuals with columns ID, Trait, X, Y, and Species, and (2) Trait_means - a named vector of random mean trait values for each species.
FUN.RandIndividuals <- function(n_spec = 10, 
                                n_individuals = 4e2,
                                n_mode = "each", # or "total"
                                Env_range = c(0, 10),
                                Trait_sd = 1,
                                seed = 42){
  set.seed(seed)
  
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
  
  ## return object
  return(list(ID_df = ID_df, Trait_means = Trait_means))
}

# RANDOM CARRYING CAPACITIES ===================================================
#' @param n_spec Numeric. Number of species / size of node set.
#' @param k_range Numeric vector of length 2. Minimum and maximum values of carrying capacity range.
#' @param seed Numeric. Seed for random processes.
#' @return A named vector of random carrying capacity values for each species.
FUN.RandK <- function(n_spec = 20, k_range = c(200,200), seed = 42){
  set.seed(seed)
  k_vec <- as.integer(runif(n = n_spec, min = k_range[1], max = k_range[2]))
  names(k_vec) <- paste("Sp", 1:n_spec, sep = "_")
  return(k_vec)
}

# SIMULATION FUNCTION ==========================================================
#' @param d0 Numeric. background death rate (gets altered by the environment and interactions).
#' @param b0 Numeric. background birth rate (remains constant).
#' @param env.xy Function with arguments X and Y. Translates X and Y coordinates into optimal local phenotype.
#' @param t_max Numeric. Maximum simulation time.
#' @param t_inter Numeric. Interval length at which to record simulation outputs (measured in simulation time).
#' @param sd Numeric. Habitat suitability in death rate function. Higher values allow individuals to persist in areas of greater environmental maladaptation.
#' @param migration Numeric. Standard deviation of 0-centred normal dsitribution from which natal dispersal is drawn.
#' @param Effect_Dis Numeric. Distance within which neighbouring individuals interact with a focal individual.
#' @param Network_igraph An igraph object with association/interaction strength stored as "weight" attribute of edges.
#' @param k_vec Named vector containing carrying capacity for each species.
#' @param ID_df Data frame of initialising individuals with columns ID, Trait, X, Y, and Species
#' @param Env_range Numeric vector of length 2. Minimum and maximum values of environmental gradient.
#' @param seed Numeric. Seed for random processes.
#' @param verbose Logical. Whether to print simulation time and sampling interval.
#' @return A list containing data frame object with the same columns as ID_df at each sampling interval defined via n_inter until t_max is reached.
SIM.Comp <- function(d0 = 0.4,
                     b0 = 0.6,
                     env.xy = function(x = NULL, y = NULL){x},
                     t_max = 10,
                     t_inter = 0.1,
                     sd = 2.5,
                     migration = 0.2,
                     Effect_Dis = 0.5,
                     Network_igraph,
                     k_vec, 
                     ID_df,
                     Env_range,
                     seed = 42,
                     verbose = TRUE, # whether to print progress in time as current time
                     RunName = ""
){
  call_info <- match.call()
  set.seed(seed)
  
  ## Network as adjacency matrix
  Effect_Mat <- igraph::as_adjacency_matrix(Network_igraph, attr = "weight") # columns affect rows
  rownames(Effect_Mat) <- colnames(Effect_Mat) <- names(k_vec)
  
  ## dynamic death rate initialisation
  ID_df <- Sim.d0Update(ID_df = ID_df, which = "Initial", 
                        env.xy = env.xy, d0 = d0, b0 = b0, sd = sd, 
                        Effect_Mat = Effect_Mat, k_vec = k_vec, 
                        Effect_Dis = Effect_Dis, seed = seed)
  
  ## list object to store individuals at each time step
  ID_ls <- list(ID_df)
  
  ## setting start times
  t <- 0 # start time at 1
  names(ID_ls)[length(ID_ls)] <- t
  ## progress bar
  if(!verbose){pb <- txtProgressBar(min = 0, max = t_max, style = 3)}
  TimeStart <- Sys.time()
  ## simulation loop over time steps
  while(t < t_max){
    # if(verbose){print(t)}
    ## vectors for storing birth and death probabilities for each individual
    birth_prob <- rep(b0, nrow(ID_df))
    death_prob <- ID_df$dt
    names(birth_prob) <- names(death_prob) <- ID_df$ID
    
    ## event identification
    EventSample_vec <- paste(rep(c("Birth", "Death"), each = nrow(ID_df)), 
                             names(birth_prob), sep="_")
    ProbSample_vec <- c(birth_prob, death_prob)
    if(any(ProbSample_vec == Inf)){
      event <- EventSample_vec[which(ProbSample_vec == Inf)[1]]
    }else{
      event <- sample(
        EventSample_vec,
        size = 1,
        prob = ProbSample_vec
      )
    }
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
      affected_row <- ID_df[nrow(ID_df),]
    }
    if(event_EV == "Death"){
      affected_row <- ID_df[ID_df$ID == as.numeric(event_ID),]
      ID_df <- ID_df[ID_df$ID != event_ID, ]
    }
    ## dynamic death rate recalculation
    ID_df <- Sim.d0Update(ID_df = ID_df, which = affected_row, event = event_EV,
                          env.xy = env.xy, d0 = d0, b0 = b0, sd = sd, 
                          Effect_Mat = Effect_Mat, k_vec = k_vec, 
                          Effect_Dis = Effect_Dis, seed = round(t, 3)*1e4+seed)
    ## Gillespie time
    ### identify by how much time advances
    tadvance <- rexp(1, rate = sum(c(birth_prob, death_prob))) 
    t <- t + tadvance
    ### record data only if interval is met
    if(t - as.numeric(names(ID_ls)[length(ID_ls)]) >= t_inter){
      # if(verbose){message(t)}
      ID_ls <- c(ID_ls, list(ID_df))
      names(ID_ls)[length(ID_ls)] <- t
      saveobj <- list(Call = call_info, Network = Network_igraph, K = k_vec, Simulation = ID_ls)
      save(saveobj, 
           file = paste0("TEMP_SIM_", RunName, "-", seed, ".RData"))
    }
    ## update progress
    if(!verbose){setTxtProgressBar(pb, t)}
    # estimator of finish
    if(verbose){
      TimeEnd <- Sys.time()
      PercRem <- ((t_max-t)/t_max)*100
      PercDone <- 100-PercRem
      if(PercDone > 100){PercDone <- 100}
      TimeElapsed <- difftime(TimeEnd, TimeStart, units = "secs")
      TimeRemaining <- (TimeElapsed/PercDone)*(100-PercDone)
      # EstimatedFinish <- TimeStart+TimeRemaining
      cat('\r', sprintf("%05s", format(round(PercDone, 2), nsmall = 2)), "%  - ", paste("Estimated time remaining =", round(lubridate::seconds_to_period(TimeRemaining), 0), "         "))
      flush.console()
    }
    if(nrow(ID_df) == 0){warning("All species went extinct"); break}
  }
  ID_ls <- c(ID_ls, list(ID_df))
  names(ID_ls)[length(ID_ls)] <- t
  unlink(paste0("TEMP_SIM_", RunName, "-", seed, ".RData"))
  return(ID_ls)
}

#' This is a helper function for updating of dynamic death rates in the simulation run
Sim.d0Update <- function(ID_df = ID_df, which = "Initial", event = NULL,
                         env.xy = env.xy, d0 = d0, b0 = b0, sd = sd, 
                         Effect_Mat = Effect_Mat, k_vec = k_vec, 
                         Effect_Dis = Effect_Dis, seed = seed){
  # set.seed(seed)
  
  ## dynamic death rate functions
  d0P <- function(N, b0, d0, k){N * (b0 - d0)/k}
  d0E <- function(Tr, Env, sd){exp((abs(Tr-Env)/sd)^2)}
  d0Omega <- function(Effect_Mat, Effect_Dis, ID_df, i){
    Abundances_i <- rep(0, ncol(Effect_Mat))
    names(Abundances_i) <- colnames(Effect_Mat)
    Abundances_obs <- table(ID_df[-i,][
      abs(ID_df[i,"X"]-ID_df[-i, "X"]) <= Effect_Dis &
        abs(ID_df[i,"Y"]-ID_df[-i, "Y"]) <= Effect_Dis,
      "Species"])
    Abundances_i[match(names(Abundances_obs), names(Abundances_i))] <- Abundances_obs
    ### Extract all effects that the focal species is subject to in the interaction matrix
    Effects_i <- Effect_Mat[which(rownames(Effect_Mat) == ID_df[i, "Species"]),]
    ### Weigh effects by abundances
    WeightedEffects_i <- Abundances_i * Effects_i
    ### Calculate final effect
    if(sum(Abundances_i) != 0){
      FinalEffect_i <- sum(WeightedEffects_i)/sum(Abundances_i)
    }else{
      FinalEffect_i <- 0
    }
    return(FinalEffect_i)
  }
  dt <- function(d0, d0P, d0E, d0Omega){d0 + d0P*d0E - d0Omega}
  ## queried d0 update
  if(which[1] == "Initial"){
    ## dynamic death rate components
    ### population dynamics
    N_vec <- table(ID_df$Species)
    d0P_vec <- d0P(N = N_vec, b0 = b0, d0 = d0, 
                   k = k_vec[match(names(k_vec), names(N_vec))])
    ID_df$d0P <- as.numeric(d0P_vec[match(ID_df$Species, names(d0P_vec))])
    ### environment
    ID_df$d0E <- d0E(Tr = ID_df$Trait, Env = env.xy(x = ID_df$X, y = ID_df$Y), 
                     sd = sd)
    ### Interactions
    ID_df$d0Omega <- sapply(1:nrow(ID_df), FUN = function(i){
      d0Omega(Effect_Mat = Effect_Mat, Effect_Dis = Effect_Dis,
              ID_df = ID_df, i = i)
    })
  }else{
    ## updating population dynamics
    if(sum(ID_df$Species == which$Species) > 0){
      d0P <- d0P(N = sum(ID_df$Species == which$Species), 
                 b0 = b0, d0 = d0, 
                 k = k_vec[names(k_vec) == which$Species])
      ID_df[ID_df$Species == which$Species, "d0P"] <- d0P
    }
    ## environment
    if(event == "Birth"){
      ID_df$d0E[ID_df$ID == which$ID] <- d0E(Tr = which$Trait, 
                                             Env = env.xy(x = which$X, 
                                                          y = which$Y), 
                                             sd = sd)
    }
    ## interactions
    ### individuals affected by addition or removal of other indiivudal
    newinterac <- ID_df[abs(ID_df$X-which$X) <= Effect_Dis &
                          abs(ID_df$Y-which$Y) <= Effect_Dis, "ID"]
    if(length(newinterac)>0){
      ID_df$d0Omega[ID_df$ID %in% newinterac] <- sapply(newinterac, FUN = function(i){
        d0Omega(Effect_Mat = Effect_Mat, Effect_Dis = Effect_Dis,
                ID_df = ID_df, i = i)
      })
    }
  }
  ID_df$dt <- dt(d0 = d0, 
                 d0P = ID_df$d0P, 
                 d0E = ID_df$d0E, 
                 d0Omega = ID_df$d0Omega)
  ID_df$dt[ID_df$dt < 0] <- 0 # make sure probabilities are never < 0
  return(ID_df)
}


# SIMULATION FRAMEWORK =========================================================
#' @param seed Numeric. Seed for random processes.
#' @param n_spec Numeric. Number of species / size of node set.
#' @param NetworkType Character. Whether to create an interaction or association matrix. In an interaction matrix, columns act on rows. In an association matrix, entries are mirrored around the diagonal. Possible values: "Association" or "Interaction".
#' @param Sparcity Numeric. Proportion of unique associations or interactions which are exactly 0.
#' @param MaxStrength Numeric. Maximum value of associations or interactions.
#' @param Env_range Numeric vector of length 2. Minimum and maximum values of environmental gradient.
#' @param Trait_sd Numeric. Standard deviation of trait disttributions from which trait values of each individual are drawn.
#' @param n_spec Numeric. Number of species / size of node set.
#' @param k_range Numeric vector of length 2. Minimum and maximum values of carrying capacity range.
#' @param d0 Numeric. background death rate (gets altered by the environment and interactions).
#' @param b0 Numeric. background birth rate (remains constant).
#' @param env.xy Function with arguments X and Y. Translates X and Y coordinates into optimal local phenotype.
#' @param t_max Numeric. Maximum simulation time.
#' @param t_inter Numeric. Interval length at which to record simulation outputs (measured in simulation time).
#' @param sd Numeric. Habitat suitability in death rate function. Higher values allow individuals to persist in areas of greater environmental maladaptation.
#' @param migration Numeric. Standard deviation of 0-centred normal dsitribution from which natal dispersal is drawn.
#' @param Effect_Dis Numeric. Distance within which neighbouring individuals interact with a focal individual.
#' @param verbose Logical. Whether to print simulation time and sampling interval.

#' @return A list containing data frame object with the same columns as ID_df at each sampling interval defined via n_inter until t_max is reached.
FUN.SimulationFramework <- function(
    seed = 42,
    ## Network Creation
    n_spec = 20,
    NetworkType = "Association", # or "Interaction"
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
    env.xy = function(x = NULL, y = NULL){x},
    t_max = 10,
    t_inter = 0.1,
    sd = 2.5,
    migration = 0.2,
    Effect_Dis = 0.5,
    verbose = TRUE,
    RunName = ""
){
  call_info <- match.call()
  
  # message("#### Package Loading")
  # install.load.package <- function(x) {
  #   if (!require(x, character.only = TRUE))
  #     install.packages(x, repos='http://cran.us.r-project.org')
  #   require(x, character.only = TRUE)
  # }
  # package_vec <- c(
  #   "randcorr", # for random correlation matrices to be used as adjacency matrices
  #   "igraph" # for representing adjacency matrices as graphs
  # )
  # sapply(package_vec, install.load.package)
  
  if(verbose){message("#### Network Creation")}
  Network_igraph <- FUN.RandNetwork(n_spec = n_spec, 
                                    NetworkType = NetworkType,
                                    Sparcity = Sparcity,
                                    MaxStrength = MaxStrength, 
                                    seed = seed)
  
  if(verbose){message("#### Initial Individual Creation")}
  RandInd_ls <- FUN.RandIndividuals(n_spec = length(V(Network_igraph)), 
                                    n_individuals = n_individuals,
                                    n_mode = n_mode,
                                    Env_range = Env_range,
                                    Trait_sd = Trait_sd,
                                    seed = seed)
  ID_df <- RandInd_ls[[1]]
  Trait_means <- RandInd_ls[[2]]
  
  if(verbose){message("#### Carrying Capacity Creation")}
  k_vec <- FUN.RandK(n_spec = length(V(Network_igraph)), 
                     k_range = k_range, 
                     seed = seed)
  
  if(verbose){message("#### Simulation Run")}
  Sim.Output <- SIM.Comp(d0 = d0, b0 = b0,
                         env.xy = env.xy,
                         k_vec = k_vec, 
                         t_max = t_max,
                         t_inter = t_inter,
                         sd = sd, 
                         ID_df = ID_df,
                         Env_range = Env_range,
                         migration = migration,
                         Effect_Dis = Effect_Dis,
                         Network_igraph = Network_igraph,
                         verbose = verbose,
                         seed = seed,
                         RunName = RunName
  )
  
  message("#### Returning Output")
  return(list(Call = call_info, Network = Network_igraph, 
              Traits = Trait_means, K = k_vec, Simulation = Sim.Output))
}
