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
Data_fs <- Data_fs[!grepl("Environment", Data_fs)]
Data_fs <- Data_fs[startsWith(Data_fs, RunName)]

# INFERENCE FUNCTIONALITY ==================================================
## Gridding of simulation outputs and summary to abundances of species per cell
Fun.Gridding <- function(
    from = Env_range[1],
    to = Env_range[2],
    n_Grid2 = n_Grid,
    ID_df = NULL) {
  GridCoords <- seq(
    from = from,
    to = to,
    length = n_Grid2 + 1
  )[-(n_Grid2 + 1)]
  grids_df <- expand.grid(GridCoords, GridCoords)
  colnames(grids_df) <- c("X", "Y")
  # grids_df <- data.frame(X = GridCoords,
  #                        Y = 0)
  grids_df$GridID <- 1:nrow(grids_df)
  GridsID_vec <- lapply(1:nrow(ID_df),
    FUN = function(x) {
      Xs <- which(ID_df[x, "X"] >= grids_df$X)
      Ys <- which(ID_df[x, "Y"] >= grids_df$Y)
      Xs[tail(which(Xs %in% Ys), 1)]
    }
  )
  grids_df$X <- grids_df$X + diff(GridCoords)[1] / 2
  grids_df$Y <- grids_df$Y + diff(GridCoords)[1] / 2
  GridsID <- unlist(GridsID_vec)
  Pop_dfBASE <- data.frame(
    matrix(0,
      nrow = nrow(grids_df),
      ncol = length(unique(ID_df$Species)) + 1
    )
  )
  colnames(Pop_dfBASE) <- c("GridsID", sort(unique(ID_df$Species)))
  Pop_dfBASE$GridsID <- grids_df$GridID
  #### observed frequencies
  Poptab <- as.data.frame.matrix(table(GridsID, ID_df$Species))
  Poptab$GridsID <- as.numeric(rownames(Poptab))
  #### matching observed with base frame
  PoptabStore <- Pop_dfBASE
  PoptabStore[match(Poptab$GridsID, PoptabStore$GridsID), -1] <- Poptab[, -ncol(Poptab)]
  PoptabStore <- PoptabStore[, -1] # rownames are grid IDs
  list(PoptabStore, grids_df)
}

## Reducing network to only survived species component and making it a matrix
Fun.SurvNetwork <- function(Igraph, ID_df) {
  NetMat_Realised <- as.matrix(as_adjacency_matrix(Igraph, attr = "weight"))
  colnames(NetMat_Realised) <- rownames(NetMat_Realised) <- V(Network_Realised)$names
  NetMat_Realised <- NetMat_Realised[
    colnames(NetMat_Realised) %in% unique(ID_df$Species),
    colnames(NetMat_Realised) %in% unique(ID_df$Species)
  ]
  return(NetMat_Realised)
}

# ACTUAL INFERENCE =========================================================
pblapply(Data_fs, cl = 20, FUN = function(Treatment_Iter) {
  # Treatment_Iter <- Data_fs[1]
  FNAME <- file.path(Dir.Models, basename(Treatment_Iter))

  if (file.exists(FNAME)) {
    1 + 1
  } else {
    load(file.path(Dir.Data, paste0(tools::file_path_sans_ext(Treatment_Iter), "_Environment.RData"))) # loads "Env_mat"
    load(file.path(Dir.Data, Treatment_Iter)) # loads list objects "SimResult", "Network_igraph", "CarryingK_vec", "Niches_vec"

    ## Simulated individuals --------------------------------------------------
    ID_df <- SimResult[[length(SimResult)]]
    ID_Grid <- Fun.Gridding(ID_df = ID_df)

    ## Networks ---------------------------------------------------------------
    Network_True <- Network_igraph
    V(Network_True)$names <- paste0("Sp_", V(Network_True))

    ### realisation
    if (RunName %in% "NoSpaceBinaryInterac") {
      Network_Realised <- Network_True
    } else {
      Network_Realised <- NetSimVal::Val.Realise(Network_True, Trait_means = Niches_vec, Effect_Dis, Env_sd)
      V(Network_Realised)$names <- V(Network_True)$names
    }

    ### turning into matrices
    NetMat_ls <- lapply(list(True = Network_True, Real = Network_Realised), FUN = function(Network) {
      NetMat <- Fun.SurvNetwork(Network, ID_df)
      diag(NetMat) <- NA
      if (!is.directed(Network_True)) {
        NetMat[lower.tri(NetMat)] <- NA
      }
      NetMat
    })


    ## Inferences -------------------------------------------------------------
    ### COOCCCUR -------------------
    mat_Iter <- ID_Grid[[1]]
    mat_Iter <- mat_Iter[rowSums(mat_Iter) != 0, ]
    mat_Iter <- sign(mat_Iter)

    model_coccurr <- cooccur(
      mat = t(mat_Iter),
      type = "spp_site",
      thresh = FALSE,
      spp_names = TRUE
    )
    model_coccur <- data.frame(
      Effects = effect.sizes(model_coccurr, standardized = TRUE),
      Sig = prob.table(model_coccurr)$p_lt < 0.05 | prob.table(model_coccurr)$p_gt < 0.05
    )
    Cooccur_results <- list(Effects = NetMat_ls$True, Sig = NetMat_ls$True)
    Cooccur_results <- lapply(names(Cooccur_results), FUN = function(x) {
      colu <- tail(which(startsWith(colnames(model_coccur), prefix = x)), 1)
      for (i in seq_len(nrow(model_coccur))) {
        sp1 <- model_coccur$Effects.sp1[i]
        sp2 <- model_coccur$Effects.sp2[i]
        Cooccur_results[[x]][sp1, sp2] <- model_coccur[i, colu]
        Cooccur_results[[x]][sp2, sp1] <- model_coccur[i, colu]
      }
      Cooccur_results[[x]]
    })
    names(Cooccur_results) <- c("Effects", "Sig")
    Cooccur_results$Sig <- Cooccur_results$Sig == 1


    ### NETASSOC -------------------
    mat_Iter <- ID_Grid[[1]]
    mat_Iter <- mat_Iter[rowSums(mat_Iter) != 0, ]

    model_netassoc <- make_netassoc_network(
      obs = t(mat_Iter),
      plot = FALSE, verbose = FALSE
    )
    Netassoc_results <- list(
      Effects = model_netassoc$matrix_spsp_ses_all,
      Sig = model_netassoc$matrix_spsp_pvalue < 0.05
    )


    ### HMSC -------------------
    # Y: site X species abundance
    Y <- ID_Grid[[1]]

    # X: covariates to be used as predictors, If you don't have covariate data, indicate this by X=NULL
    ## defione coarse boundaries
    id_grid <- ID_Grid[[2]]
    dx <- diff(sort(unique(id_grid$X)))[1]
    dy <- diff(sort(unique(id_grid$Y)))[1]
    x_centers <- sort(unique(id_grid$X))
    y_centers <- sort(unique(id_grid$Y))
    x_breaks <- c(x_centers - dx / 2, max(x_centers) + dx / 2)
    y_breaks <- c(y_centers - dy / 2, max(y_centers) + dy / 2)
    ## make env_df with cells and corresponding coarse cells
    env_df <- expand.grid(
      X = as.numeric(rownames(Env_mat)),
      Y = as.numeric(colnames(Env_mat))
    )
    env_df$value <- as.vector(Env_mat)
    env_df$Xc <- x_centers[
      findInterval(env_df$X, x_breaks, rightmost.closed = TRUE)
    ]
    env_df$Yc <- y_centers[
      findInterval(env_df$Y, y_breaks, rightmost.closed = TRUE)
    ]
    env_df <- merge(
      env_df,
      id_grid,
      by.x = c("Xc", "Yc"),
      by.y = c("X", "Y"),
      all.x = FALSE
    )
    grid_means <- aggregate(value ~ GridID, env_df, mean)

    X <- data.frame(
      X = ID_Grid[[2]]$X,
      Y = ID_Grid[[2]]$Y,
      Environment = grid_means$value
    )

    XFormula <- as.formula("~ X+Y")
    if (length(unique(X$Environment)) == 1) {
      XFormula <- as.formula("~ X+Y")
    } else {
      XFormula <- as.formula("~ X+Y+Environment")
    }

    # S: study design, including units of study and their possible coordinates, If you don't have variables that define the study design, indicate this by S=NULL
    S <- data.frame(GridID = ID_Grid[[2]]$GridID)

    studyDesign <- data.frame(GridID = as.factor(S$GridID))
    St <- studyDesign$GridID
    rL.site <- HmscRandomLevel(units = levels(St))
    xy <- X[, c("X", "Y")]
    rownames(xy) <- studyDesign[, 1]
    rL.coords <- HmscRandomLevel(sData = xy)

    # Tr: species traits (note that T is a reserved word in R and that's why we use Tr); If you don't have trait data, indicate this by Tr=NULL.
    Traits <- aggregate(Trait ~ Species, data = ID_df, FUN = mean)
    Traits_vec <- Traits$Trait
    names(Traits_vec) <- Traits$Species
    Tr <- data.frame(Trait = Traits_vec[match(colnames(Y), names(Traits_vec))])

    if (length(unique(X$Environment)) == 1) {
      Tr <- TrFormula <- NULL
    } else {
      TrFormula <- as.formula("~ Trait")
    }


    # P: phylogenetic information given by taxonomical levels, e.g. order, family, genus, species; If TP does not have phylogenetic data (because you don't have such data at all, or because, it is given in tree-format, like is the case in this example), indicate this with P=NULL
    P <- NULL

    # MODEL
    HMSC_mod <- Hmsc(
      Y = Y, XData = X, XFormula = XFormula,
      TrData = Tr, TrFormula = TrFormula,
      distr = "poisson",
      studyDesign = studyDesign,
      ranLevels = {
        list("GridID" = rL.site)
      }
    )
    hmsc_model <- sampleMcmc(HMSC_mod,
      samples = nSamples, thin = thin,
      transient = nWarmup,
      nChains = nChains,
      nParallel = 1
      # nChains
    )
    OmegaCor <- computeAssociations(hmsc_model)

    # ESTIMATES
    # mean effect
    me <- as.data.frame(OmegaCor[[1]]$mean)
    diag(me) <- NA
    ## support for positive association
    po <- as.data.frame(OmegaCor[[1]]$support)
    diag(po) <- NA

    HMSC_results <- list(
      Effects = me,
      Sig = po > 0.95 | po < 0.05
    )

    ## Save Object -----------------------------------------------------------
    Inference_ls <- list(
      Cooccur = Cooccur_results,
      Netassoc = Netassoc_results,
      HMSC = HMSC_results
    )
    save(
      Inference_ls,
      file = FNAME
    )
  }
})
