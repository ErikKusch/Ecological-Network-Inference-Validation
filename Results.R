#' ####################################################################### #
#' PROJECT: [PhD; X - DATA SIMULATIONS] 
#' CONTENTS: 
#'  - Compare Inferred and Known Networks
#'  - Tally Error Rates for Inferred Networks
#'  DEPENDENCIES:
#'  - Inference.R
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
source("Inference.R")
package_vec <- c(package_vec, "brms", "rethinking", "reshape2", "cowplot", "scales")
sapply(package_vec, install.load.package)
Dir.Base <- getwd() # read out the project directory
Dir.Concept <- file.path(Dir.Base, "Concept")
Dirs <- c(Dir.Concept)
CreateDir <- sapply(Dirs, function(x) if(!dir.exists(x)) dir.create(x))

# BETA DIVERSITY ===========================================================
message("Comparison of Dissimilarities")
Dissimilarities <- do.call(rbind, 
                           pblapply(Inference_ls, FUN = function(x){
                             data.frame(
                               # Naive = x$Naive$Dissimilarity$OS, 
                               HMSC = x$Informed$Dissimilarity$OS,
                               COOCCUR = x$COOCCUR$Dissimilarity$OS
                             )
                           })
)
Dissimilarities_df <- data.frame(OS = c(Dissimilarities$COOCCUR,
                                        Dissimilarities$HMSC),
                                 Mode = rep(c("COOCCUR", "HMSC"), 
                                            each = nrow(Dissimilarities))
)

OS_gg <- ggplot(Dissimilarities_df, 
                aes(x = OS, y = factor(Mode, levels = c("COOCCUR", "HMSC")))) + 
  stat_halfeye() + 
  theme_bw() + labs(y = "")
ggsave(OS_gg, filename = file.path(Dir.Exports, "OS.png"), 
       width = 16, height = 9, units = "cm")

# DETECTION REGRESSION =====================================================
message("Detection of Positive, Negative, and Absent Associations")
Detection_ls <- pblapply(names(Inference_ls), FUN = function(SimName){
  x <- Inference_ls[[SimName]]
  ## matrices of networks
  ### true/known
  True_mat <- as_adjacency_matrix(x$Informed$Graphs$True, 
                                  attr = "weight", type = "upper", sparse = FALSE) 
  Weighted_mat <- as_adjacency_matrix(x$TrueNetwork, 
                                      attr = "weight", type = "upper", sparse = FALSE)
  colnames(Weighted_mat) <- rownames(Weighted_mat) <- colnames(True_mat)
  
  #### inferred
  HMSC_mat <- as_adjacency_matrix(x$Informed$Graphs$HMSC, 
                                      attr = "weight", type = "upper", sparse = FALSE)
  COOCCUR_mat <- as_adjacency_matrix(x$COOCCUR$Graphs$COOCCUR, 
                                   attr = "weight", type = "upper", sparse = FALSE)
  
  ### sorting
  True_mat <- True_mat[order(rownames(True_mat)),order(colnames(True_mat))]
  Weighted_mat <- Weighted_mat[order(rownames(Weighted_mat)),order(colnames(Weighted_mat))]
  COOCCUR_mat <- COOCCUR_mat[order(rownames(COOCCUR_mat)),order(colnames(COOCCUR_mat))]
  HMSC_mat <- HMSC_mat[order(rownames(HMSC_mat)),order(colnames(HMSC_mat))]
  
  InferenceMats_ls <- list(HMSC = HMSC_mat,
                           COOCCUR = COOCCUR_mat)
  
  ### limitting to species who did not go extinct in simulation
  True_mat <- True_mat[rownames(True_mat) %in% rownames(COOCCUR_mat),
                       colnames(True_mat) %in% colnames(COOCCUR_mat)]
  Weighted_mat <- Weighted_mat[rownames(Weighted_mat) %in% rownames(HMSC_mat),
                               colnames(Weighted_mat) %in% colnames(HMSC_mat)]
  TruePos_mat <- TrueNeg_mat <- TrueAbs_mat <- True_mat
  TruePos_mat[TruePos_mat != 1] <- NA
  TrueNeg_mat[TrueNeg_mat != -1] <- NA
  TrueAbs_mat[TrueAbs_mat != 0] <- NA
  
  ## environmental preferences of species
  Traits <- aggregate(Trait ~ Species, data = x$ID_df, FUN = mean)
  EnvDiff_mat <- True_mat
  for(i in rownames(EnvDiff_mat)){
    focaltrait <- Traits$Trait[Traits$Species == i]
    alltraitsorder <- Traits$Trait[match(Traits$Species, colnames(EnvDiff_mat))]
    EnvDiff_mat[i,] <- abs(focaltrait - alltraitsorder)
  }
  
  ## inference classification
  InferenceCalss_ls <- lapply(names(InferenceMats_ls), FUN = function(k){
    data.frame(
      CorrectPos = as.vector(as.numeric((TruePos_mat + InferenceMats_ls[[k]]) == 2)), # true positive
      CorrectNeg = as.vector(as.numeric((TrueNeg_mat + InferenceMats_ls[[k]]) == -2)), # true negative
      CorrectAbsent = as.vector(as.numeric(((TrueAbs_mat == 0) + (InferenceMats_ls[[k]] == 0) == 2))), # true absent
      Magnitude = abs(as.vector(Weighted_mat)),
      SignTrue = as.vector(sign(Weighted_mat)),
      SIgnInferred = as.vector(InferenceMats_ls[[k]]),
      EnvDiff = as.vector(EnvDiff_mat),
      Mode = k,
      Sim = SimName
    )
  })
  Inference_df <- do.call(rbind, InferenceCalss_ls)
  Inference_df
})
Detect_df <- do.call(rbind, Detection_ls)

## Correct Identification of Positive Associations -----------------------------
print("Models of Positive Associations")
if(file.exists(file.path(Dir.Exports, "Bayes_Model_Positive.RData"))){
  load(file.path(Dir.Exports, "Bayes_Model_Positive.RData"))
}else{
  Bayes_Model_Positive <- list(HMSC = NA, COOCCUR = NA)
  for(i in names(Bayes_Model_Positive)){
    model_df <- Detect_df[Detect_df$Mode == i, ]
    Bayes_Model_Positive[[i]] <- brm(formula = CorrectPos ~Magnitude * EnvDiff,
                                data = model_df,
                                family = bernoulli(link = "logit"),
                                warmup = nWarmup,
                                iter = nSamples,
                                chains = nChains,
                                cores = nChains,
                                seed = 42)
  }
  save(Bayes_Model_Positive, 
       file = file.path(Dir.Exports, "Bayes_Model_Positive.RData"))
}

## Correct Identification of Negative Associations -----------------------------
print("Models of Negative Associations")
if(file.exists(file.path(Dir.Exports, "Bayes_Model_Negative.RData"))){
  load(file.path(Dir.Exports, "Bayes_Model_Negative.RData"))
}else{
  Bayes_Model_Negative <- list(HMSC = NA, COOCCUR = NA)
  for(i in names(Bayes_Model_Negative)){
    model_df <- Detect_df[Detect_df$Mode == i, ]
    Bayes_Model_Negative[[i]] <- brm(formula = CorrectNeg ~Magnitude * EnvDiff,
                                     data = model_df,
                                     family = bernoulli(link = "logit"),
                                     warmup = nWarmup,
                                     iter = nSamples,
                                     chains = nChains,
                                     cores = nChains,
                                     seed = 42)
  }
  save(Bayes_Model_Negative, 
       file = file.path(Dir.Exports, "Bayes_Model_Negative.RData"))
}

## Correct Identification of Absent Associations -------------------------------
print("Models of Absent Associations")
if(file.exists(file.path(Dir.Exports, "Bayes_Model_Absent.RData"))){
  load(file.path(Dir.Exports, "Bayes_Model_Absent.RData"))
}else{
  Bayes_Model_Absent <- list(HMSC = NA, COOCCUR = NA)
  for(i in names(Bayes_Model_Absent)){
    model_df <- Detect_df[Detect_df$Mode == i, ]
    Bayes_Model_Absent[[i]] <- brm(formula = CorrectAbsent ~Magnitude * EnvDiff,
                                     data = model_df,
                                     family = bernoulli(link = "logit"),
                                     warmup = nWarmup,
                                     iter = nSamples,
                                     chains = nChains,
                                     cores = nChains,
                                     seed = 42)
  }
  save(Bayes_Model_Absent, 
       file = file.path(Dir.Exports, "Bayes_Model_Absent.RData"))
}

## Plotting --------------------------------------------------------------------
print("Plotting")
Plot_ls <- list(Positive = Bayes_Model_Positive,
                Negative = Bayes_Model_Negative,
                Absent = Bayes_Model_Absent
                )
ProbMat <- expand.grid(seq(from = 0, to = 1, by = 0.05),
                       seq(from = 0, to = 10, by = 0.5))
colnames(ProbMat) <- c("Magnitude", "EnvDiff")

for(k in names(Plot_ls)){
  Bayes_ls <- Plot_ls[[k]]
  
  HMSC_df <- posterior_samples(Bayes_ls[["HMSC"]])
  HMSC_df <- reshape2::melt(HMSC_df[,1:4])
  HMSC_df$Method <- "HMSC"
  COOCCUR_df <- posterior_samples(Bayes_ls[["COOCCUR"]])
  COOCCUR_df <- reshape2::melt(COOCCUR_df[,1:4])
  COOCCUR_df$Method <- "COOCCUR"
  
  post_plot <- rbind(HMSC_df, COOCCUR_df)
  post_plot$variable <- c("Intercept", "Magnitude", "Environmental \n Difference", "Interaction")[match(post_plot$variable, unique(post_plot$variable))]
  
  
  Coeff_gg <- ggplot(post_plot, aes(x = value, 
                                    y = factor(variable, levels = rev(c("Intercept", "Magnitude", "Environmental \n Difference", "Interaction"))))) + 
    stat_halfeye() + 
    labs(y = "Model Coefficient", x = "Logit Value") + 
    facet_wrap(~Method, ncol = 1) + 
    theme_bw()
  
  for(i in names(Bayes_ls)){
    post_df <- posterior_samples(Bayes_ls[[i]])
    prob_iter <- ProbMat
    prob_iter$Prob <- apply(prob_iter, MARGIN = 1, FUN = function(x){
      mean(inv_logit(post_df$b_Intercept + 
                       post_df$b_EnvDiff * as.numeric(x[2]) + 
                       (post_df$b_Magnitude + post_df$`b_Magnitude:EnvDiff` * as.numeric(x[2])) * as.numeric(x[1])
      )
      )
    })
    prob_iter$Method <- i
    if(exists("prob_plot")){
      prob_plot <- rbind(prob_plot, prob_iter)
    }else{
      prob_plot <- prob_iter
    }
  }
  
  
  Prob_gg <- ggplot(prob_plot, aes(x = EnvDiff, y = Magnitude, fill = Prob)) + 
    geom_tile() + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Probability of \n Indetification")) + 
    theme_bw() + 
    labs(y = "Association Magnitude", x = "Environmental Difference") +
    scale_fill_viridis_c(option = "E", limits = c(0,1)) + 
    facet_wrap(~ Method, ncol = 1)
  
  Plot_ls[[k]] <- cowplot::plot_grid(Coeff_gg, Prob_gg)
  
  ggsave(Plot_ls[[k]], 
         filename = file.path(Dir.Exports, paste0("Detection_", k, ".png")), 
         width = 30, height = 20, units = "cm")
  
  rm(prob_plot)
}

# ERROR RATES ==============================================================
message("Inference Error Rates")
ErrorRates_df <- do.call(rbind, 
                         pblapply(Inference_ls, FUN = function(x){
                           # matrix extraction
                           True_mat <- as_adjacency_matrix(x$Informed$Graphs$True, 
                                                           attr = "weight", type = "upper", sparse = FALSE) 
                           HMSC_mat <- as_adjacency_matrix(x$Informed$Graphs$HMSC, 
                                                               attr = "weight", type = "upper", sparse = FALSE)
                           COOCCUR_mat <- as_adjacency_matrix(x$COOCCUR$Graphs$COOCCUR, 
                                                            attr = "weight", type = "upper", sparse = FALSE)
                           # sorting
                           True_mat <- True_mat[order(rownames(True_mat)),order(colnames(True_mat))]
                           COOCCUR_mat <- COOCCUR_mat[order(rownames(COOCCUR_mat)),order(colnames(COOCCUR_mat))]
                           HMSC_mat <- HMSC_mat[order(rownames(HMSC_mat)),order(colnames(HMSC_mat))]
                           ### limitting to species who did not go extinct in simulation
                           True_mat <- True_mat[rownames(True_mat) %in% rownames(COOCCUR_mat),
                                                colnames(True_mat) %in% colnames(COOCCUR_mat)]
                           
                           # metrics
                           mat_ls <- list(COOCCUR = COOCCUR_mat,
                                          HMSC = HMSC_mat)
                           
                           do.call(rbind, lapply(names(mat_ls),
                                                 function(y){
                                                   data.frame(
                                                     Values = c(
                                                       # true positive associations; how many of the inferred positive associations are correctly identified as such?
                                                       TP = sum((True_mat + mat_ls[[y]]) == 2)/sum(mat_ls[[y]] == 1), 
                                                       # true negative associations; how many of the inferred negative associations are correctly identified as such?
                                                       TN = sum((True_mat + mat_ls[[y]]) == -2)/sum(mat_ls[[y]] == -1), 
                                                       # falsely inferred positive; how many of the inferred positive associations are incorrectly identified as such?
                                                       FP = sum((True_mat == 0) + (mat_ls[[y]] == 1) == 2)/sum(mat_ls[[y]] == 1), 
                                                       # falsely inferred negative; how many of the inferred negative associations are incorrectly identified as such?
                                                       FN = sum((True_mat == 0) + (mat_ls[[y]] == -1) == 2)/sum(mat_ls[[y]] == -1),
                                                       # true positive links which were not inferred; how many of the true positive associations were not inferred?
                                                       MP = sum((True_mat == 1) + (mat_ls[[y]] == 0) == 2)/sum(True_mat == 1), 
                                                       # true negative links which were not inferred; how many of the true negative associations were not inferred?
                                                       MN = sum((True_mat == -1) + (mat_ls[[y]] == 0) == 2)/sum(True_mat == -1),
                                                       # true absent; how many of the inferred absent interactions are correctly identified as such?
                                                       TA = sum((True_mat == 0) + (mat_ls[[y]] == 0) == 2)/sum(mat_ls[[y]] == 0),
                                                       # false absent; how many of the inferred absent interactions are correctly identified as such?
                                                       FA = sum((True_mat != 0) + (mat_ls[[y]] == 0) == 2)/sum(mat_ls[[y]] == 0)
                                                     ),
                                                     Metric = c("TP", "TN", "FP", "FN", "MP", "MN", "TA", "FA"),
                                                     # Identifier
                                                     Mode = y
                                                   )
                                                 }))
                         })
)

ErrorRates_gg <- ggplot(ErrorRates_df, 
                        aes(x = Values, 
                            y = factor(Mode, levels = c("COOCCUR", "HMSC")))) + 
  stat_halfeye() + 
  facet_wrap(~factor(Metric, levels = c("TP", "TN", "FP", "FN", "MP", "MN", "TA", "FA")),
             ncol = 2) + 
  theme_bw() + labs(x = "Detection Rate [%]", y = "")
ggsave(ErrorRates_gg, filename = file.path(Dir.Exports, "ErrorRates.png"), 
       width = 20, height = 24, units = "cm")

# CONCEPT VISUALISATION ========================================================
message("Conceptual Visualisation")
source("SimulationFrameworkFunctions.R")

## Data Generation -------------------------------------------------------------
print("Data Generation")
if(file.exists(file.path(Dir.Concept, "ConceptSim.RData"))){
  load(file.path(Dir.Concept, "ConceptSim.RData"))
}else{
  start <- Sys.time()
  Simulation_Output <- FUN.SimulationFramework(
    seed = 21,
    ## Network Creation
    n_spec = 5,
    NetworkType = "Association", # or "Association"
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
  end <- Sys.time()
  print(end-start)
  save(Simulation_Output, file = file.path(Dir.Concept, "ConceptSim.RData")) 
}

## Input Data Visualisation ----------------------------------------------------
print("Input Data Visualisation")

### starting constellation ----
first_df <- Simulation_Output$Simulation[[1]]
first_df$Species <- factor(as.numeric(gsub(first_df$Species, pattern = "Sp_", replacement = "")))
Traits_df <- data.frame(Traits = Simulation_Output$Traits,
                        Species = names(Simulation_Output$Traits))
Traits_df$Species <- factor(gsub(Traits_df$Species, 
                                 pattern = "Sp_", replacement = ""))
Initial_gg <- ggplot(first_df, 
                     aes(x = X, y = Y, col = Species, shape = Species)) + 
  geom_point() + 
  scale_shape_manual(values=1:nlevels(first_df$Species)) + 
  scale_color_viridis_d() + 
  geom_vline(data = Traits_df, 
             aes(xintercept = Traits, 
                 col = Species)) + 
  xlim(eval(Simulation_Output$Call[["Env_range"]])[1],
       eval(Simulation_Output$Call[["Env_range"]])[2]) + 
  theme_bw()

### ending constellation ----
ID_df <- Simulation_Output$Simulation[[length(Simulation_Output$Simulation)]]
last_df <- Simulation_Output$Simulation[[length(Simulation_Output$Simulation)]]
last_df$Species <- factor(as.numeric(gsub(last_df$Species, pattern = "Sp_", replacement = "")))
Traits_df <- aggregate(Trait ~ Species, data = ID_df, FUN = mean)
Traits_df$Species <- factor(as.numeric(gsub(Traits_df$Species, 
                                            pattern = "Sp_", replacement = "")))
Final_gg <- ggplot(last_df, 
                   aes(x = X, y = Y, col = Species, shape = Species)) + 
  geom_point() + 
  scale_shape_manual(values=1:nlevels(last_df$Species)) + 
  scale_color_viridis_d() + 
  geom_vline(data = Traits_df, 
             aes(xintercept = Trait, 
                 col = Species)) + 
  xlim(eval(Simulation_Output$Call[["Env_range"]])[1],
       eval(Simulation_Output$Call[["Env_range"]])[2]) +
  theme_bw()
leg <- get_legend(Final_gg + theme(legend.position = "bottom"))

### constellation plotting ----
Input_gg <- plot_grid(plot_grid(Initial_gg + theme(legend.position = "none"), 
                    Final_gg + theme(legend.position = "none"), 
                    nrow = 1),
          leg, rel_heights = c(1,0.1), ncol = 1)
ggsave(Input_gg, filename = file.path(Dir.Exports, "SpatialInputs.png"), 
       width = 30, height = 16, units = "cm")

### species-site matrices ----
## make locational data into site X species matrix
GridCoords <- seq(from = eval(Simulation_Output$Call[["Env_range"]])[1], 
                  to = eval(Simulation_Output$Call[["Env_range"]])[2], 
                  length = n_Grid+1)[-(n_Grid+1)]
grids_df <- expand.grid(GridCoords, GridCoords)
colnames(grids_df) <- c("X", "Y")
grids_df$GridID <- 1:nrow(grids_df)
GridsID_vec <- lapply(1:nrow(ID_df),
                      FUN = function(x){
                        Xs <- which(ID_df[x, "X"] >= grids_df$X)
                        Ys <- which(ID_df[x, "Y"] >= grids_df$Y)
                        Xs[tail(which(Xs %in% Ys), 1)]
                      }
)
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

### Plotting
grid_vis <- ggplot(last_df, 
                   aes(x = X, y = Y, col = Species, shape = Species)) + 
  geom_point() + 
  scale_shape_manual(values=1:nlevels(last_df$Species)) + 
  scale_color_viridis_d() + 
  geom_vline(xintercept = GridCoords+diff(GridCoords)[1]/2) + 
  geom_hline(yintercept = GridCoords+diff(GridCoords)[1]/2) + 
  xlim(eval(Simulation_Output$Call[["Env_range"]])[1],
       eval(Simulation_Output$Call[["Env_range"]])[2]) +
  theme_bw()

SPSite_df <- cbind(Y, (grids_df+diff(GridCoords)[1]/2)[,1:2])
loop_vec <- colnames(SPSite_df)[-((ncol(SPSite_df)-1):ncol(SPSite_df))]
plot_col <- viridis_pal()(length(loop_vec))
names(plot_col) <- loop_vec
plot_ls <- as.list(rep(NA, length = (length(loop_vec)+1)))
plot_ls[[1]] <- grid_vis
names(plot_ls) <- c("gridvis", loop_vec)
for(i in loop_vec){
  abund_df <- SPSite_df[, c(i, "X", "Y")]
  colnames(abund_df)[1] <- "Abund"
  plot_ls[[i]] <- ggplot(abund_df, aes(x = X, y = Y, fill = Abund)) + 
    geom_tile(color = "black", lwd = 0.5, linetype = 1) + 
    coord_fixed() +
    scale_fill_gradient(low = "grey",
                        high = plot_col[names(plot_col) == i]) + 
    guides(fill = guide_colourbar(barwidth = 2,
                                  barheight = 15,
                                  title = "Associatiuon")) + 
    theme_bw() + theme(legend.position = "none") + theme(plot.margin = unit(c(0,0,0,0), "cm")) + labs(y = "", x = "")
}

InputMatrices_gg <- plot_grid(plot_ls[[1]],
          plot_grid(plotlist = plot_ls[-1], nrow = 1),
          ncol = 1, rel_heights = c(3, 1.7), labels = c("A", "B")
          )
ggsave(InputMatrices_gg, filename = file.path(Dir.Exports, "InputMatrices.png"), 
       width = 30, height = 20, units = "cm")


## Abundance Through Time ------------------------------------------------------
print("Abundance Visualisation through Time")
Abund_time <- pblapply(names(Simulation_Output$Simulation), 
                       FUN = function(t){
                         ID_iter <- Simulation_Output$Simulation[[t]]
                         cbind(data.frame(table(ID_iter$Species)), t)
                       })
Abund_time <- do.call(rbind, Abund_time)
Abund_time$t <- as.numeric(Abund_time$t)
colnames(Abund_time)[1:2] <- c("Species", "Abundance")
Abund_time$Species <- factor(as.numeric(gsub(Abund_time$Species, 
                       pattern = "Sp_", replacement = "")))

AbundTime_gg <- ggplot(Abund_time, aes(x = t, y = Abundance, col = Species)) +
  geom_line(size = 1.5) + 
  scale_color_viridis_d() + 
  theme_bw()
ggsave(AbundTime_gg, filename = file.path(Dir.Exports, "AbundanceTime.png"), 
       width = 30, height = 20, units = "cm")  

## Spatial Gradient ------------------------------------------------------------
print("Spatial Gradient Visualisation")
env.xy <- function(x = NULL, y = NULL){x}
gridseq <- seq(from = eval(Simulation_Output$Call[["Env_range"]])[1],
               to = eval(Simulation_Output$Call[["Env_range"]])[2],
               length = 5e1)
gridmat <- expand.grid(gridseq, gridseq)
colnames(gridmat) <- c("x", "y")
gridmat$Phenotype <- apply(gridmat, MARGIN = 1, FUN = function(k){
  env.xy(x = k[1], y = k[2])
})

Gradient_gg <- ggplot(gridmat, aes(x = x, y = y, fill = Phenotype)) + 
  geom_tile() + 
  coord_fixed() + 
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Optimal \n Phenotype")) + 
  theme_bw() + theme(plot.margin = unit(c(0,0,0,0), "cm"))
ggsave(Gradient_gg, filename = file.path(Dir.Exports, "SpatialGradient.png"), 
       width = 20, height = 22, units = "cm")  

## Network Inference -----------------------------------------------------------
print("Network Inference")
source2 <- function(file, start, end, ...) {
  file.lines <- scan(file, what=character(), skip=start-1, nlines=end-start+1, sep='\n')
  file.lines.collapsed <- paste(file.lines, collapse='\n')
  source(textConnection(file.lines.collapsed), ...)
}
source2("Inference.R", 62, 313)

models_ls <- FUN.Inference(Simulation_Output, 
                           Dir.Models = Dir.Concept,
                           Dir.Exports = Dir.Concept,
                           Treatment_Iter = "Concept",
                           Cores = 1)

print("Input network matrix")
net_mat <- as_adjacency_matrix(Simulation_Output$Network, attr = "weight",
                               type = "upper", sparse = FALSE)
net_mat[lower.tri(net_mat)] <- NA
diag(net_mat) <- NA
colnames(net_mat) <- rownames(net_mat) <- V(Simulation_Output$Network)
edg_df1 <- melt(net_mat)
colnames(edg_df1) <- c("Partner 1", "Partner 2", "Strength")
TrueMat_gg <- ggplot(edg_df1, aes(x = `Partner 1`, y = `Partner 2`, fill = Strength)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) + 
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Associatiuon")) + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
  scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")

ggsave(TrueMat_gg, filename = file.path(Dir.Exports, "Matrix_True.png"), 
       width = 20, height = 22, units = "cm")  

print("HMSC network matrix")
if(length(E(models_ls$Informed$Graphs$HMSC)) != 0){
  net_mat <- as_adjacency_matrix(models_ls$Informed$Graphs$HMSC, attr = "weight",
                                 type = "upper", sparse = FALSE)
}else{
  net_mat <- matrix(rep(0, length = length(V(models_ls$Informed$Graphs$HMSC))^2), ncol = length(V(models_ls$Informed$Graphs$HMSC)))
  colnames(net_mat) <- rownames(net_mat) <- sort(V(models_ls$Informed$Graphs$HMSC)$name)
}

net_mat[lower.tri(net_mat)] <- NA
diag(net_mat) <- NA
colnames(net_mat) <- rownames(net_mat) <- V(models_ls$Informed$Graphs$HMSC)
edg_df <- melt(net_mat)
colnames(edg_df) <- c("Partner 1", "Partner 2", "Strength")
edg_df$Method <- "HMSC"
edg_df2 <- edg_df

print("COOCCUR network matrix")
net_mat <- as_adjacency_matrix(models_ls$COOCCUR$Graphs$COOCCUR, attr = "weight",
                               type = "upper", sparse = FALSE)
net_mat[lower.tri(net_mat)] <- NA
diag(net_mat) <- NA
colnames(net_mat) <- rownames(net_mat) <- V(models_ls$COOCCUR$Graphs$COOCCUR)
edg_df <- melt(net_mat)
colnames(edg_df) <- c("Partner 1", "Partner 2", "Strength")
edg_df$Method <- "COOCCUR"
edg_df2 <- rbind(edg_df2, edg_df)
# edg_df2$Strength <- as.numeric(edg_df2$Strength)
edg_df2$Correct <- NA
edg_df2$Correct[which((edg_df2$Strength == 1) + (sign(edg_df1$Strength) == 1) == 2)] <- 1
edg_df2$Correct[which((edg_df2$Strength == -1) + (sign(edg_df1$Strength) == -1) == 2)] <- 1
edg_df2$Correct[which((edg_df2$Strength == 0) + (sign(edg_df1$Strength) == 0) == 2)] <- 1
edg_df2$Correct[which(is.na(edg_df2$Correct) & !is.na(edg_df2$Strength))] <- 0
edg_df2$Correct <- factor(edg_df2$Correct)

InfMat_gg <- ggplot(edg_df2, 
                    aes(x = `Partner 1`, y = `Partner 2`, 
                        fill = Strength, shape = Correct)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) + 
  coord_fixed() +
  geom_point(size = 3) + 
  scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "",
                     guide = "none") +  
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Associatiuon")) + 
  theme_bw() + 
  facet_wrap(~Method) + 
  theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
  scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")

ggsave(InfMat_gg, filename = file.path(Dir.Exports, "Matrix_Inferred.png"), 
       width = 40, height = 22, units = "cm")  

