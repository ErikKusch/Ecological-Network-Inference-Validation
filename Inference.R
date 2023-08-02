#' ####################################################################### #
#' PROJECT: [PhD; X - DATA SIMULATIONS] 
#' CONTENTS: 
#'  - Infer ecological networks from simulation outputs
#'  - Compute network dissimilarities
#'  DEPENDENCIES:
#'  - DataSimulations_HMSC.R must have been run
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
rm(list=ls())
set.seed(42)

## DIRECTORIES -------------------------------------------------------------
Dir.Base <- getwd() # read out the project directory
Dir.Data <- file.path(Dir.Base, "Data")
Dir.Models <- file.path(Dir.Base, "Models")
Dir.Exports <- file.path(Dir.Base, "Exports")
Dirs <- c(Dir.Data, Dir.Models, Dir.Exports)
CreateDir <- sapply(Dirs, function(x) if(!dir.exists(x)) dir.create(x))

## Packages ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "Hmsc", # for HMSC models
  "cooccur", # for COOCCUR models
  "igraph", # for graph representation
  "devtools", # to install betalink
  "betalink", # for computation of network dissimilarities
  "pbapply", # for parallel processing
  "doParallel", # for clustercalls
  "ggplot2", # for plotting
  "tidybayes" # for plotting
)
sapply(package_vec, install.load.package)

if("betalink" %in% rownames(installed.packages()) == FALSE){
  install_version("betalink", version = "2.2.1", repos = "http://cran.us.r-project.org")
}
library(betalink)

`%nin%` <- Negate(`%in%`)

## Bayes Settings ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
nSamples <- 7e3
# nSamples <- 2e3
thin <- 1
nWarmup <- round(nSamples*0.3*thin, 0)
nChains <- 4
n_Grid <- 10
message(paste0("thin = ",as.character(thin),"; samples = ",as.character(nSamples)))

# DATA =====================================================================
Data_fs <- list.files(Dir.Data, pattern = ".RData")
# Data_fs <- Data_fs[unlist(lapply(strsplit(Data_fs, split = "_"), "[[", 2)) == 1]

# INFERENCE FUNCTIONS ======================================================
FUN.Inference <- function(Simulation_Output = NULL, 
                          Dir.Exports = NULL,
                          Dir.Models = NULL,
                          ModelSave = FALSE,
                          Treatment_Iter = NULL,
                          Cores = 1){
  ### DATA READOUT ####
  ## Simulated individuals
  ID_df <- Simulation_Output$Simulation[[length(Simulation_Output$Simulation)]]
  
  ## Networks
  V(Simulation_Output$Network)$names <- paste0("Sp_", V(Simulation_Output$Network))
  Network_Weighted <- Simulation_Output$Network
  mode <- ifelse(is.directed(Network_Weighted), "directed", "undirected")
  Network_Realised <- as.matrix(as_adjacency_matrix(Network_Weighted, attr = "weight"))
  colnames(Network_Realised) <- rownames(Network_Realised) <- V(Network_Weighted)$names
  Network_Realised <- Network_Realised[colnames(Network_Realised) %in% unique(ID_df$Species),
                   colnames(Network_Realised) %in% unique(ID_df$Species)]
  
  ## Network of only survived species till end of simulation
  Network_SurvNonReal <- graph_from_adjacency_matrix(Network_Realised, mode = mode, weighted = TRUE)
  Network_Weighted <- Network_SurvNonReal
  E(Network_SurvNonReal)$weight <- ifelse(E(Network_SurvNonReal)$weight > 0, 1, -1)
  
  ## Network of only survived species till end of simulation and within reach of each other to actually interact
  # Figuring out trait differences between potentially interacting species
  SPTrait_df <- aggregate(ID_df, Trait ~ Species, FUN = mean)
  SPTrait_df$SD <- aggregate(ID_df, Trait ~ Species, FUN = sd)$Trait
  SPTrait_df <- SPTrait_df[match(colnames(Network_Realised), SPTrait_df$Species), ]
  SPTrait_mat <- abs(outer(SPTrait_df$Trait, SPTrait_df$Trait, '-'))
  colnames(SPTrait_mat) <- rownames(SPTrait_mat) <- colnames(Network_Realised)
  SPTraitSD_mat <- abs(outer(SPTrait_df$SD, SPTrait_df$SD, '+'))
  colnames(SPTraitSD_mat) <- rownames(SPTraitSD_mat) <- colnames(Network_Realised)
  TraitDiff_mat <- SPTrait_mat-SPTraitSD_mat
  
  # limitting to realised interactions
  Network_Realised[(TraitDiff_mat) > Simulation_Output$Call$sd+Simulation_Output$Call$Effect_Dis] <- 0 # anything greater apart in enviro pref than the interaction window (0.5) + environmental sd cannot be realised
  Network_SurvReal <- igraph::graph_from_adjacency_matrix(adjmatrix = Network_Realised,
                                                      mode = mode,
                                                      weighted = TRUE,
                                                      diag = FALSE)
  Network_WeightedReal <- Network_SurvReal
  E(Network_SurvReal)$weight <- ifelse(E(Network_SurvReal)$weight > 0, 1, -1)
  
  ### DATA PREPRATION ####
  # Y: Site X Species matrix
  ## make locational data into site X species matrix
  GridCoords <- seq(from = eval(Simulation_Output$Call[["Env_range"]])[1], 
                    to = eval(Simulation_Output$Call[["Env_range"]])[2], 
                    length = n_Grid+1)[-(n_Grid+1)]
  grids_df <- expand.grid(GridCoords, GridCoords)
  colnames(grids_df) <- c("X", "Y")
  # grids_df <- data.frame(X = GridCoords,
  #                        Y = 0)
  grids_df$GridID <- 1:nrow(grids_df)
  GridsID_vec <- lapply(1:nrow(ID_df),
                        FUN = function(x){
                          Xs <- which(ID_df[x, "X"] >= grids_df$X)
                          Ys <- which(ID_df[x, "Y"] >= grids_df$Y)
                          Xs[tail(which(Xs %in% Ys), 1)]
                        }
  )
  grids_df$X <- grids_df$X+diff(GridCoords)[1]/2
  grids_df$Y <- grids_df$Y+diff(GridCoords)[1]/2
  GridsID <- unlist(GridsID_vec)
  Pop_dfBASE <- data.frame(
    matrix(0, nrow = nrow(grids_df),
           ncol = length(unique(ID_df$Species))+1
    )
  )
  colnames(Pop_dfBASE) <- c("GridsID", sort(unique(ID_df$Species)))
  Pop_dfBASE$GridsID <- grids_df$GridID
  #### observed frequencies
  Poptab <- as.data.frame.matrix(table(GridsID, ID_df$Species))
  Poptab$GridsID <- as.numeric(rownames(Poptab))
  #### matching observed with base frame
  PoptabStore <- Pop_dfBASE
  PoptabStore[match(Poptab$GridsID, PoptabStore$GridsID), -1] <- Poptab[,-ncol(Poptab)]
  ### storing site X species matrix
  Y <- PoptabStore[,-1] # rownames are grid IDs
  
  # X: covariates to be used as predictors, If you don't have covariate data, indicate this by X=NULL
  X <- data.frame(
    # X = grids_df$X,
    # Y = grids_df$Y,
    Environment = grids_df$X
  )
  rownames(X) <- grids_df$GridID
  XNaive <- data.frame(
    # X = grids_df$X,
    # Y = grids_df$Y,
    Environment = rep(1, nrow(grids_df))
  )
  rownames(XNaive) <- grids_df$GridID
  
  # S: study design, including units of study and their possible coordinates, If you don't have variables that define the study design, indicate this by S=NULL
  S <- data.frame(GridID = grids_df$GridID)
  
  # Tr: species traits (note that T is a reserved word in R and that's why we use Tr); If you don't have trait data, indicate this by Tr=NULL.
  Traits <- aggregate(Trait ~ Species, data = ID_df, FUN = mean)
  Traits_vec <- Traits$Trait
  names(Traits_vec) <- Traits$Species
  Tr <- data.frame(Trait = Traits_vec[match(colnames(Y), names(Traits_vec))])
  
  # P: phylogenetic information given by taxonomical levels, e.g. order, family, genus, species; If TP does not have phylogenetic data (because you don't have such data at all, or because, it is given in tree-format, like is the case in this example), indicate this with P=NULL  
  P <- NULL 
  
  ### DATA CHECKS ####
  if(!is.numeric(as.matrix(Y)) || !is.logical(as.matrix(Y)) && !is.finite(sum(Y, na.rm=TRUE))){
    stop("Species data should be numeric and have finite values")}
  if(any(is.na(S))){stop("study design has NA values - not allowed for")}
  if(any(is.na(X))){stop("Covariate data has NA values - not allowed for")}
  if(any(is.na(P))){stop("P has NA values - not allowed for")}
  
  ### Model Specification ----
  ## Model Formulae
  XFormula <- as.formula("~ Environment")
  # XFormula2 <- as.formula("~ X+Y")
  TrFormula <- as.formula("~ Trait")
  ## StudyDesign
  studyDesign <- data.frame(GridID = as.factor(S$GridID))
  St <- studyDesign$GridID
  rL.site <- HmscRandomLevel(units = levels(St))
  xy <- data.frame(X = grids_df$X,
                   Y = grids_df$Y)
  rownames(xy) <- studyDesign[,1]
  rL.coords <- HmscRandomLevel(sData=xy)
  
  ## Models
  InformedMod <- Hmsc(Y = Y, XData = X,  XFormula = XFormula,
                      TrData = Tr, TrFormula = TrFormula,
                      distr = "poisson",
                      studyDesign = studyDesign,
                      ranLevels = {list("GridID" = rL.site)}
  )
  NaiveMod <- Hmsc(Y = Y, XData = XNaive,  XFormula = ~1,
                   TrData = NULL, TrFormula = NULL,
                   distr = "poisson",
                   studyDesign = studyDesign,
                   ranLevels = {list("GridID" = rL.site)}
  )
  models_ls <- list(InformedMod
                    # ,
                    # NaiveMod
                    )
  names(models_ls) <- c("Informed"
                        # ,"Naive"
                        )
  
  ### Modelling and Dissimilarity Computation ----
  # message("Modelling")
  for(Model_Iter in 1:length(models_ls)){
    Start_t <- Sys.time()
    print(Model_Iter)
    hmsc_model <- models_ls[[Model_Iter]]
    
    ModelPath <- file.path(Dir.Models, 
                           paste0(tools::file_path_sans_ext(Treatment_Iter), "_", 
                                  names(models_ls)[Model_Iter], ".Rdata"))
    
    if(file.exists(ModelPath)){
      load(ModelPath)
    }else{
      hmsc_model <- sampleMcmc(hmsc_model, samples = nSamples, thin = thin,
                               transient = nWarmup,
                               nChains = nChains,
                               nParallel = Cores
                               # nChains
      )
      if(ModelSave){save(hmsc_model, file = ModelPath)}
    }
    OmegaCor <- computeAssociations(hmsc_model)
    supportLevel <- 0.95
    me <- as.data.frame(OmegaCor[[1]]$mean)
    me <- cbind(hmsc_model$spNames,me)
    colnames(me)[1] <- ""
    po <- as.data.frame(OmegaCor[[1]]$support)
    po <- cbind(hmsc_model$spNames,po)
    colnames(po)[1] <- ""
    ne <- as.data.frame(1-OmegaCor[[1]]$support)
    ne <- cbind(hmsc_model$spNames,ne)
    colnames(ne)[1] <- ""
    vals <- list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
    
    ### Interaction/Association Matrix ----
    Interaction_mean <- vals$`Posterior mean`[,-1]
    Interaction_ProbPos <- vals$`Pr(x>0)`[,-1]
    Interaction_ProbNeg <- vals$`Pr(x<0)`[,-1]
    Partner2 <- c()
    for(i in 1:(length(colnames(Interaction_mean))-1)){
      Partner2 <- c(Partner2, colnames(Interaction_mean)[-c(1:i)])
    }
    Interactions_igraph <- data.frame(Partner1 = rep(rownames(Interaction_mean), 
                                                     times = (length(colnames(Interaction_mean))-1):0),
                                      Partner2 = Partner2,
                                      Inter_mean = t(Interaction_mean)[lower.tri(t(Interaction_mean), diag = FALSE)],
                                      Inter_ProbPos = t(Interaction_ProbPos)[lower.tri(t(Interaction_ProbPos), diag = FALSE)],
                                      Inter_ProbNeg = t(Interaction_ProbNeg)[lower.tri(t(Interaction_ProbNeg), diag = FALSE)]
    )
    Interactions_HMSC <- Interactions_igraph[order(abs(Interactions_igraph$Inter_mean), decreasing = TRUE), ]
    
    ### Network Dissimilarity ----
    Interactions_HMSC$Sig <- FALSE
    Interactions_HMSC$Sig[
      Interactions_HMSC$Inter_ProbPos >= 0.95 | 
        Interactions_HMSC$Inter_ProbNeg >= 0.95] <- TRUE
    
    HMSC_ig <- graph_from_data_frame(Interactions_HMSC[Interactions_HMSC$Sig == TRUE,], 
                                     directed = FALSE)
    Spec_vec <- unique(c(Interactions_HMSC$Partner1, Interactions_HMSC$Partner2))
    if(length(E(HMSC_ig)) != 0){
      E(HMSC_ig)$weight <- E(HMSC_ig)$Inter_mean
      E(HMSC_ig)$weight[E(HMSC_ig)$weight > 0] <- 1
      E(HMSC_ig)$weight[E(HMSC_ig)$weight < 0] <- -1
      origvert <- names(V(HMSC_ig))
      # HMSC_ig <- set.vertex.attribute(HMSC_ig, "names", value = names(V(HMSC_ig)))
      
      addvert <- Spec_vec[Spec_vec %nin% names(V(HMSC_ig))
      ]
      HMSC_ig <- add_vertices(HMSC_ig,
                              nv = length(addvert))
      HMSC_ig <- set.vertex.attribute(HMSC_ig, "name", value = c(origvert, addvert))
    }else{
      HMSC_ig <- make_empty_graph(n = length(Spec_vec))
      HMSC_ig <- set.vertex.attribute(HMSC_ig, "names", value = Spec_vec)
    }
    
    # V(Simulation_Output$Network)$name <- paste0("Sp_", V(Simulation_Output$Network))
    # E(Simulation_Output$Network)$weight <- ifelse(E(Simulation_Output$Network)$weight > 0, 1, -1)
    
    Graphs_ls <- list(HMSC = HMSC_ig, 
                      SurvNonReal = Network_SurvNonReal,
                      SurvReal = Network_SurvReal)
    betadiv <- network_betadiversity(N = Graphs_ls)[-3,]
    
    ### Export ----
    End_t <- Sys.time()
    models_ls[[Model_Iter]] <- list(
      # HMSC_model = hmsc_model,
      HMSC_associations = Interactions_HMSC,
      Graphs = Graphs_ls,
      Dissimilarity = betadiv,
      Duration = End_t - Start_t,
      grids = n_Grid
    )
    print(models_ls[[Model_Iter]]$Duration)
  } # model loop
  
  ## COCCUR
  mat_Iter <- Y
  mat_Iter[mat_Iter > 0] <- 1
  model_coccurr <- cooccur(mat = t(mat_Iter), type = "spp_site", 
                           thresh = FALSE, spp_names = TRUE)
  Interac_df <- effect.sizes(model_coccurr, standardized = TRUE)
  Interac_df$pLT <- prob.table(model_coccurr)$p_lt
  Interac_df$pGT <- prob.table(model_coccurr)$p_gt
  Interac_df$Sig <- Interac_df[,4] < 0.05 | Interac_df[,5] < 0.05
  colnames(Interac_df)[1:2] <- c("Partner1", "Partner2")
  
  COOCCUR_ig <- graph_from_data_frame(Interac_df[Interac_df$Sig == TRUE,], 
                                      directed = mode)
  Spec_vec <- unique(c(Interac_df$Partner1, Interac_df$Partner2))
  if(length(E(COOCCUR_ig)) != 0){
    E(COOCCUR_ig)$weight <- E(COOCCUR_ig)$effects
    E(COOCCUR_ig)$weight[E(COOCCUR_ig)$weight > 0] <- 1
    E(COOCCUR_ig)$weight[E(COOCCUR_ig)$weight < 0] <- -1
    origvert <- names(V(COOCCUR_ig))
    addvert <- Spec_vec[Spec_vec %nin% names(V(COOCCUR_ig))
    ]
    COOCCUR_ig <- add_vertices(COOCCUR_ig, 
                               nv = length(addvert))
    COOCCUR_ig <- set.vertex.attribute(COOCCUR_ig, "name", value = c(origvert, addvert))
  }else{
    COOCCUR_ig <- make_empty_graph(n = length(Spec_vec))
    COOCCUR_ig <- set.vertex.attribute(COOCCUR_ig, "name", value = Spec_vec)
  }
  
  # V(Simulation_Output$Network)$name <- paste0("Sp_", V(Simulation_Output$Network))
  # E(Simulation_Output$Network)$weight <- ifelse(E(Simulation_Output$Network)$weight > 0, 1, -1)
  Graphs_ls <- list(COOCCUR = COOCCUR_ig, 
                    SurvNonReal = Network_SurvNonReal,
                    SurvReal = Network_SurvReal)
  betadiv <- network_betadiversity(N = Graphs_ls)[-3, ]
  
  models_ls$COOCCUR <- list(COOCCUR_associations = Interac_df,
                            Graphs = Graphs_ls,
                            Dissimilarity = betadiv
  )
  
  ## FINAL OUTPUT
  # models_ls$TrueNetwork <- Network_Weighted
  models_ls$ID_df <- ID_df
  models_ls$Weighted <- list(NonReal = Network_Weighted,
                             Real = Network_WeightedReal)
  save(models_ls, file = file.path(Dir.Exports, Treatment_Iter)) 
  return(models_ls)
}

# ANALYSIS =================================================================
message("############ INFERENCE & NETWORK DISSIMILARITY COMPUTATION")

print("Registering Cluster")
nCores <- ifelse(parallel::detectCores()>25, 25, parallel::detectCores())
cl <- parallel::makeCluster(nCores) # for parallel pbapply functions
parallel::clusterExport(cl,
                        varlist = c("Data_fs", "nSamples", "thin", "nWarmup", "nChains",
                                    "%nin%", "Dir.Data", "Dir.Models", "Dir.Exports",
                                    "install.load.package", "package_vec", "n_Grid", "FUN.Inference"),
                        envir = environment()
)
clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))

print("Inference of Networks")
Inference_ls <- pblapply(1:length(Data_fs), 
                         cl = cl,
                         FUN = function(Treatment_Enum){
                           Treatment_Iter <- Data_fs[Treatment_Enum]
                           
                           if(file.exists(file.path(Dir.Exports, Treatment_Iter))){
                             load(file.path(Dir.Exports, Treatment_Iter))
                           }else{
                             load(file.path(Dir.Data, Treatment_Iter)) # loads list object "SimulationOutput"
                             if(nrow(Simulation_Output$Simulation[[length(Simulation_Output$Simulation)]]) == 0){
                               models_ls <- NA
                             }else{
                               models_ls <- FUN.Inference(Simulation_Output =  Simulation_Output,
                                                          Dir.Exports = Dir.Exports,
                                                          Dir.Models = Dir.Models,
                                                          Treatment_Iter = Treatment_Iter) 
                             }
                           }
                           
                           # reporting back to top
                           models_ls
                           
                         })
names(Inference_ls) <- Data_fs
# stop("Remove NA entries from Inference_ls - these are instances were the simulation resulted in extinction of all species")