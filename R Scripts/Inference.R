#' ####################################################################### #
#' PROJECT: [InfVal; HMSC-Inference of Pairwise Associations] 
#' CONTENTS: 
#'  - Infer ecological networks from simulation outputs
#'  - Compute network dissimilarities
#'  DEPENDENCIES:
#'  -
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
message("Running Association Inference")

# DATA =====================================================================
Data_fs <- list.files(Dir.Data, pattern = ".RData")
Data_fs <- Data_fs[grepl(RunName, Data_fs)]
Data_fs <- Data_fs[!grepl("META", Data_fs)]

# INFERENCE FUNCTION =======================================================
FUN.Inference <- function(Simulation_Output = NULL, 
                          Dir.Exports = NULL,
                          Dir.Models = NULL,
                          ModelSave = TRUE,
                          Treatment_Iter = NULL,
                          n_Grid = 10,
                          Cores = 1){
  ### DATA READOUT ####
  ## Simulated individuals
  ID_df <- Simulation_Output$Simulation[[length(Simulation_Output$Simulation)]]
  ID2_df <- Simulation_Output$Simulation[[(length(Simulation_Output$Simulation)-1)]]
  
  ## Networks
  V(Simulation_Output$Network)$names <- paste0("Sp_", V(Simulation_Output$Network))
  Network_Weighted <- Simulation_Output$Network
  mode <- ifelse(is.directed(Network_Weighted), "directed", "undirected")
  Network_Realised <- as.matrix(as_adjacency_matrix(Network_Weighted, attr = "weight"))
  colnames(Network_Realised) <- rownames(Network_Realised) <- V(Network_Weighted)$names
  Network_Realised <- Network_Realised[colnames(Network_Realised) %in% unique(ID_df$Species),
                   colnames(Network_Realised) %in% unique(ID_df$Species)]
  NonReal_mat <- Network_Realised
  diag(NonReal_mat) <- NA
  if(!is.directed(Simulation_Output$Network)){
    NonReal_mat[lower.tri(NonReal_mat) ] <- NA
  }
  
  ## Network of only survived species till end of simulation
  Network_SurvNonReal <- graph_from_adjacency_matrix(Network_Realised, mode = mode, weighted = TRUE)
  Network_Weighted <- Network_SurvNonReal
  E(Network_SurvNonReal)$weight <- ifelse(E(Network_SurvNonReal)$weight > 0, 1, -1)
  
  ## Network of only survived species till end of simulation and within reach of each other to actually interact
  # Figuring out trait differences between potentially interacting species
  SPTrait_df <- data.frame(Simulation_Output$Traits)
  SPTrait_df$Species <- rownames(SPTrait_df)
  colnames(SPTrait_df) <- c("Trait", "Species")
  # SPTrait_df <- aggregate(ID_df, Trait ~ Species, FUN = mean) # might want to consider whether to delimit this by initialising or final trait values
  # SPTrait_df$SD <- aggregate(ID_df, Trait ~ Species, FUN = sd)$Trait
  SPTrait_df$SD <- 1
  SPTrait_df <- SPTrait_df[match(colnames(Network_Realised), SPTrait_df$Species), ]
  SPTrait_mat <- abs(outer(SPTrait_df$Trait, SPTrait_df$Trait, '-'))
  colnames(SPTrait_mat) <- rownames(SPTrait_mat) <- colnames(Network_Realised)
  SPTraitSD_mat <- abs(outer(SPTrait_df$SD, SPTrait_df$SD, '+'))
  colnames(SPTraitSD_mat) <- rownames(SPTraitSD_mat) <- colnames(Network_Realised)
  TraitDiff_mat <- SPTrait_mat #-SPTraitSD_mat
  
  # limitting to realised interactions
  Network_Realised[(TraitDiff_mat) > eval(Simulation_Output$Call$sd)+eval(Simulation_Output$Call$Effect_Dis)] <- 0 # anything greater apart in enviro pref than the interaction window (0.5) + environmental sd cannot be realised
  Real_mat <- Network_Realised
  diag(Real_mat) <- NA
  if(!is.directed(Simulation_Output$Network)){
    Real_mat[lower.tri(Real_mat) ] <- NA
  }
  Network_SurvReal <- igraph::graph_from_adjacency_matrix(adjmatrix = Network_Realised,
                                                      mode = mode,
                                                      weighted = TRUE,
                                                      diag = FALSE)
  Network_WeightedReal <- Network_SurvReal
  E(Network_SurvReal)$weight <- ifelse(E(Network_SurvReal)$weight > 0, 1, -1)
  
  # True Matrices
  mat_ls <- list(NonRealised = NonReal_mat,
                 Realised = Real_mat)
  
  ### DATA PREPRATION ####
  # Y: Site X Species matrix
  ## make locational data into site X species matrix
  Fun.Gridding <- function(from = eval(Simulation_Output$Call[["Env_range"]])[1], 
                           to = eval(Simulation_Output$Call[["Env_range"]])[2],
                           n_Grid2 = n_Grid,
                           ID_df = NULL){
    GridCoords <- seq(from = from, 
                      to = to, 
                      length = n_Grid2+1)[-(n_Grid2+1)]
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
    PoptabStore <- PoptabStore[,-1] # rownames are grid IDs
    list(PoptabStore, grids_df)
  }
  grids_df <- Fun.Gridding(ID_df = ID_df)[[2]]
  Abundances <- Fun.Gridding(ID_df = ID_df)[[1]]
  Abundances2 <- Fun.Gridding(ID_df = ID2_df)[[1]]
  Performances <- Abundances-Abundances2
  Occurrences <- sign(Abundances)
  
  ### storing site X species matrix
  Y <- list(Occurrence = Occurrences,
            Abundance = Abundances,
            Performance = Performances)
  # print(Y)
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
  X <- list(Naive = XNaive,
            Informed = X)
  
  # S: study design, including units of study and their possible coordinates, If you don't have variables that define the study design, indicate this by S=NULL
  S <- data.frame(GridID = grids_df$GridID)
  
  # Tr: species traits (note that T is a reserved word in R and that's why we use Tr); If you don't have trait data, indicate this by Tr=NULL.
  Traits <- aggregate(Trait ~ Species, data = ID_df, FUN = mean)
  Traits_vec <- Traits$Trait
  names(Traits_vec) <- Traits$Species
  Tr <- data.frame(Trait = Traits_vec[match(colnames(Y$Occurrence), names(Traits_vec))])
  
  # P: phylogenetic information given by taxonomical levels, e.g. order, family, genus, species; If TP does not have phylogenetic data (because you don't have such data at all, or because, it is given in tree-format, like is the case in this example), indicate this with P=NULL  
  P <- NULL 
  
  ### DATA CHECKS ####
  test_ls <- lapply(Y, FUN = function(Y){if(!is.numeric(as.matrix(Y)) || !is.logical(as.matrix(Y)) && !is.finite(sum(Y, na.rm=TRUE))){
    stop("Species data should be numeric and have finite values")}})
  test_ls <- lapply(X, FUN = function(X){if(any(is.na(X))){stop("Covariate data has NA values - not allowed for")}})
  if(any(is.na(S))){stop("study design has NA values - not allowed for")}
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
  models_ls <- lapply(names(Y), FUN = function(z){
    if(z == "Occurrence"){
      InformedMod <- Hmsc(Y = Y[[z]], XData = X$Informed,  XFormula = XFormula,
                          TrData = Tr, TrFormula = TrFormula,
                          distr = "probit",
                          studyDesign = studyDesign,
                          ranLevels = {list("GridID" = rL.site)}
      )
      NaiveMod <- Hmsc(Y = Y[[z]], XData = X$Naive,  XFormula = ~1,
                       TrData = NULL, TrFormula = NULL,
                       distr = "probit",
                       studyDesign = studyDesign,
                       ranLevels = {list("GridID" = rL.site)}
      )
    }else{
      InformedMod <- Hmsc(Y = Y[[z]], XData = X$Informed,  XFormula = XFormula,
                          TrData = Tr, TrFormula = TrFormula,
                          distr = "poisson",
                          studyDesign = studyDesign,
                          ranLevels = {list("GridID" = rL.site)}
      )
      NaiveMod <- Hmsc(Y = Y[[z]], XData = X$Naive,  XFormula = ~1,
                       TrData = NULL, TrFormula = NULL,
                       distr = "poisson",
                       studyDesign = studyDesign,
                       ranLevels = {list("GridID" = rL.site)}
      )
      }
    models_ls <- list(InformedMod,NaiveMod)
    names(models_ls) <- c("Informed","Naive")
    models_ls
  })
  names(models_ls) <- names(Y)
  models_ls <- unlist(models_ls, recursive = FALSE)
  
  ### Modelling and Dissimilarity Computation ----
  # message("Modelling")
  HMSCmat_ls <- as.list(rep(NA, length(models_ls)))
  names(HMSCmat_ls) <- names(models_ls)
  for(Model_Iter in 1:length(models_ls)){
    Start_t <- Sys.time()
    print(names(models_ls)[Model_Iter])
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
      OmegaCor <- computeAssociations(hmsc_model)
      if(ModelSave){save(OmegaCor, file = ModelPath)}
    }
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
    
    if(length(E(HMSC_ig)) != 0){
      HMSC_ig <- permute(as.undirected(HMSC_ig), match(V(HMSC_ig)$name, colnames(Real_mat)))
      net_mat <- as_adjacency_matrix(HMSC_ig, attr = "weight",
                                     type = "upper", sparse = FALSE)
    }else{
      net_mat <- matrix(rep(0, length = length(V(HMSC_ig))^2), ncol = length(V(HMSC_ig)))
      colnames(net_mat) <- rownames(net_mat) <- sort(V(HMSC_ig)$name)
    }
    net_mat[lower.tri(net_mat)] <- NA
    diag(net_mat) <- NA
    colnames(net_mat) <- rownames(net_mat) <- V(HMSC_ig)$name
    HMSC_mat <- net_mat
    HMSCmat_ls[[Model_Iter]] <- HMSC_mat
    
    betadiv <- data.frame(i = rep(names(models_ls)[Model_Iter], 2),
                          j = c("SurvNonReal", "SurvReal"),
                          Accuracy = c(FUN_Matcomparison(mat1 = sign(NonReal_mat), mat2 = HMSC_mat),
                                       FUN_Matcomparison(mat1 = sign(Real_mat), mat2 = HMSC_mat))
                          )
  
    Graph <- HMSC_ig
    
    ### Export ----
    End_t <- Sys.time()
    models_ls[[Model_Iter]] <- list(
      # HMSC_model = hmsc_model,
      HMSC_associations = Interactions_HMSC,
      Graph = Graph,
      Dissimilarity = betadiv,
      Duration = End_t - Start_t,
      grids = n_Grid
    )
    print(models_ls[[Model_Iter]]$Duration)
  } # model loop
  
  # ## COOCCUR
  # mat_Iter <- Occurrences
  # mat_Iter[mat_Iter > 0] <- 1
  # model_coccurr <- cooccur(mat = t(mat_Iter), type = "spp_site", 
  #                          thresh = FALSE, spp_names = TRUE)
  # Interac_df <- effect.sizes(model_coccurr, standardized = TRUE)
  # Interac_df$pLT <- prob.table(model_coccurr)$p_lt
  # Interac_df$pGT <- prob.table(model_coccurr)$p_gt
  # Interac_df$Sig <- Interac_df[,4] < 0.05 | Interac_df[,5] < 0.05
  # colnames(Interac_df)[1:2] <- c("Partner1", "Partner2")
  # 
  # COOCCUR_ig <- graph_from_data_frame(Interac_df[Interac_df$Sig == TRUE,], 
  #                                     directed = mode)
  # Spec_vec <- unique(c(Interac_df$Partner1, Interac_df$Partner2))
  # if(length(E(COOCCUR_ig)) != 0){
  #   E(COOCCUR_ig)$weight <- E(COOCCUR_ig)$effects
  #   E(COOCCUR_ig)$weight[E(COOCCUR_ig)$weight > 0] <- 1
  #   E(COOCCUR_ig)$weight[E(COOCCUR_ig)$weight < 0] <- -1
  #   origvert <- names(V(COOCCUR_ig))
  #   addvert <- Spec_vec[Spec_vec %nin% names(V(COOCCUR_ig))
  #   ]
  #   COOCCUR_ig <- add_vertices(COOCCUR_ig, 
  #                              nv = length(addvert))
  #   COOCCUR_ig <- set.vertex.attribute(COOCCUR_ig, "name", value = c(origvert, addvert))
  # }else{
  #   COOCCUR_ig <- make_empty_graph(n = length(Spec_vec))
  #   COOCCUR_ig <- set.vertex.attribute(COOCCUR_ig, "name", value = Spec_vec)
  # }
  # 
  # # V(Simulation_Output$Network)$name <- paste0("Sp_", V(Simulation_Output$Network))
  # # E(Simulation_Output$Network)$weight <- ifelse(E(Simulation_Output$Network)$weight > 0, 1, -1)
  # COOCCUR_ig <- permute(as.undirected(COOCCUR_ig), match(V(COOCCUR_ig)$name, colnames(Real_mat)))
  # 
  # net_mat <- as_adjacency_matrix(as.undirected(COOCCUR_ig), attr = "weight",
  #                                type = "upper", sparse = FALSE)
  # net_mat[lower.tri(net_mat)] <- NA
  # diag(net_mat) <- NA
  # colnames(net_mat) <- rownames(net_mat) <- V(COOCCUR_ig)$name
  # COOCCUR_mat <- net_mat
  # 
  # Graph <- COOCCUR_ig
  # 
  # # betadiv <- network_betadiversity(N = Graphs_ls)[-3, ]
  # 
  # betadiv <- data.frame(i = c("COOCCUR", "COOCCUR"),
  #                       j = c("SurvNonReal", "SurvReal"),
  #                       Accuracy = c(FUN_Matcomparison(mat1 = sign(NonReal_mat), mat2 = COOCCUR_mat),
  #                                    FUN_Matcomparison(mat1 = sign(Real_mat), mat2 = COOCCUR_mat))
  # )
  # 
  # models_ls$COOCCUR <- list(COOCCUR_associations = Interac_df,
  #                           Graph = Graph,
  #                           Dissimilarity = betadiv
  # )
  # 
  ## FINAL OUTPUT
  models_ls$mats <- list(True = mat_ls,
                         HMSC = HMSCmat_ls
                         # ,
                         # COOCCUR = COOCCUR_mat
                         )
  models_ls$Graphs <- list(SurvNonReal = Network_SurvNonReal,
                           SurvReal = Network_SurvReal)
  models_ls$ID_df <- ID_df
  models_ls$Weighted <- list(NonReal = Network_Weighted,
                             Real = Network_WeightedReal)
  save(models_ls, file = file.path(Dir.Exports, Treatment_Iter)) 
  return(models_ls)
}

# ANALYSIS =================================================================
parallel::clusterExport(cl,
                        varlist = c("FUN.Inference"),
                        envir = environment()
)
clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))

Inference_ls <- pblapply(Data_fs, 
                         cl = cl,
                         FUN = function(Treatment_Iter){
                           print(Treatment_Iter)
                           if(file.exists(file.path(Dir.Exports, Treatment_Iter))){
                             load(file.path(Dir.Exports, Treatment_Iter))
                           }else{
                             load(file.path(Dir.Data, Treatment_Iter)) # loads list object "SimulationOutput"
                             if(nrow(Simulation_Output$Simulation[[length(Simulation_Output$Simulation)]]) == 0){
                               models_ls <- NA
                             }else{
                               models_ls <- FUN.Inference(Simulation_Output = Simulation_Output,
                                                          Dir.Exports = Dir.Exports,
                                                          Dir.Models = Dir.Models,
                                                          Treatment_Iter = Treatment_Iter,
                                                          n_Grid = n_Grid,
                                                          ModelSave = TRUE)
                             }
                           }
                           # reporting back to top
                           models_ls
                         })
names(Inference_ls) <- Data_fs