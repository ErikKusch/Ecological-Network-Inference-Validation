#' ####################################################################### #
#' PROJECT: [PhD; X - DATA SIMULATIONS] 
#' CONTENTS: 
#'  - Compare Inferred and Known Networks
#'  - Tally Error Rates for Inferred Networks
#'  DEPENDENCIES:
#'  - Inference.R
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =====================================================================
source("Inference.R")
# Inference_ls <- pblapply(
#   list.files(Dir.Exports, pattern = "Association_.*.RData", full.names = TRUE),
#   FUN = function(x){
#     load(x)
#     models_ls
#   })
# names(Inference_ls) <- list.files(Dir.Exports, pattern = "Association_.*.RData")
# OneSD <- unlist(lapply(strsplit(names(Inference_ls), split = "_"), "[[", 2)) == 1
# Inference_ls <- Inference_ls[OneSD]
package_vec <- c(package_vec, "brms", "rethinking", "reshape2", 
                 "cowplot", "scales", "ggnewscale")
sapply(package_vec, install.load.package)
Dir.Base <- getwd() # read out the project directory
Dir.Concept <- file.path(Dir.Base, "Concept")
Dirs <- c(Dir.Concept)
CreateDir <- sapply(Dirs, function(x) if(!dir.exists(x)) dir.create(x))

# INFERRED VS. INFERRED ========================================================
InfComp <- pblapply(Inference_ls, FUN = function(Sim){
  mats <- Sim$mats$HMSC
  # mats$COOCCUR <- Sim$mats$COOCCUR
  comp_df <- expand.grid(names(mats), names(mats))
  comp_df$Difference <- apply(comp_df, MARGIN = 1, FUN = function(z){
    FUN_Matcomparison(mats[[z[1]]], mats[[z[2]]])
  })
  comp_df
})
InfComp_df <- as.data.frame(do.call(rbind, InfComp))

cases <- as.character(unique(InfComp_df$Var1))
cases <- c(
  # "COOCCUR", 
           rev(cases[startsWith(cases, "O")]),
           rev(cases[startsWith(cases, "A")]), rev(cases[startsWith(cases, "P")]))

InfComp_df$Var1 <- factor(InfComp_df$Var1, levels = cases)
InfComp_df$Var2 <- factor(InfComp_df$Var2, levels = cases)

mean_df <- aggregate(Difference ~ Var1 + Var2, InfComp_df, FUN = mean)
nums <- nrow(mean_df)
names <- levels(mean_df$Var1)
mean_df <- reshape(mean_df, idvar = "Var1", timevar = "Var2", direction = "wide")[,-1]
colnames(mean_df) <- rownames(mean_df) <- names
# mean_df <- mean_df[match(rownames(mean_df), cases), match(colnames(mean_df), rev(cases))]
diag(mean_df) <- NA
mean_df[lower.tri(mean_df)] <- NA

sd_df <- aggregate(Difference ~ Var1 + Var2, InfComp_df, FUN = sd)
names <- levels(sd_df$Var1)
sd_df <- reshape(sd_df, idvar = "Var1", timevar = "Var2", direction = "wide")[,-1]
colnames(sd_df) <- rownames(sd_df) <- names
diag(sd_df) <- NA
sd_df[upper.tri(sd_df)] <- NA

plot_df <- rbind(as.data.frame(as.table(as.matrix(mean_df))),
                 as.data.frame(as.table(as.matrix(sd_df)))
)
plot_df$metric <- c(rep("Mean", nums), rep("SD", nums))

plot_df$Var1 <- gsub(as.character(plot_df$Var1), pattern = "Occurrence.", replacement = "[O]")
plot_df$Var1 <- gsub(as.character(plot_df$Var1), pattern = "Abundance.", replacement = "[A]")
plot_df$Var1 <- gsub(as.character(plot_df$Var1), pattern = "Performance.", replacement = "[P]")
plot_df$Var1 <- gsub(as.character(plot_df$Var1), pattern = "Informed", replacement = "Clim")
plot_df$Var1 <- gsub(as.character(plot_df$Var1), pattern = "Naive", replacement = "")

plot_df$Var2 <- gsub(as.character(plot_df$Var2), pattern = "Occurrence.", replacement = "[O]")
plot_df$Var2 <- gsub(as.character(plot_df$Var2), pattern = "Abundance.", replacement = "[A]")
plot_df$Var2 <- gsub(as.character(plot_df$Var2), pattern = "Performance.", replacement = "[P]")
plot_df$Var2 <- gsub(as.character(plot_df$Var2), pattern = "Informed", replacement = "Clim")
plot_df$Var2 <- gsub(as.character(plot_df$Var2), pattern = "Naive", replacement = "")

cases <- gsub(as.character(cases), pattern = "Occurrence.", replacement = "[O]")
cases <- gsub(as.character(cases), pattern = "Abundance.", replacement = "[A]")
cases <- gsub(as.character(cases), pattern = "Performance.", replacement = "[P]")
cases <- gsub(as.character(cases), pattern = "Informed", replacement = "Clim")
cases <- gsub(as.character(cases), pattern = "Naive", replacement = "")

ggplot(mapping = aes(x = factor(Var1, levels = cases), y = factor(Var2, levels = cases))) +
  geom_tile(data = plot_df[plot_df$metric == "Mean",], aes(fill = Freq)) +
  scale_fill_viridis_c(option = "E", na.value="transparent", name = "Mean", position = "left") +
  # scale_fill_gradient(low = "darkred", high = "forestgreen", na.value="transparent", 
  #                     name = "Mean", position = "left") +
  new_scale_fill() +
  geom_tile(data = plot_df[plot_df$metric == "SD",], aes(fill = Freq)) +
  # scale_fill_gradient(low = "beige", high = "darkblue", na.value="transparent", 
  #                     name = "SD", position = "right") + 
  scale_fill_viridis_c(option = "B", na.value="transparent", name = "SD", position = "right") +
  theme_bw() + 
  labs(x = "", y = "", title = "Inference Similarity")

ggsave(filename = file.path(Dir.Exports, "Fig_InfVInf.png"), 
       width = 16, height = 9, units = "cm")

# BETA DIVERSITY ===============================================================
message("Comparison of Dissimilarities")
Dissimilarities_df <- do.call(rbind, 
                              pblapply(Inference_ls, FUN = function(x){
                                binding_df <- do.call(rbind, lapply(x, "[[", "Dissimilarity"))
                                rownames(binding_df) <- c()
                                binding_df
                              })
)
Dissimilarities_df <- Dissimilarities_df[Dissimilarities_df$i != "COOCCUR", ]
Dissimilarities_df$j <- gsub(Dissimilarities_df$j, pattern = "SurvNonReal", replacement = "Full True Network")
Dissimilarities_df$j <- gsub(Dissimilarities_df$j, pattern = "SurvReal", replacement = "Realisable True Network")

# cases <- unique(Dissimilarities_df$i)
# cases <- rev(c("COOCCUR", rev(cases[startsWith(cases, "O")]), 
#   rev(cases[startsWith(cases, "A")]), rev(cases[startsWith(cases, "P")])))

Dissimilarities_df$i <- gsub(as.character(Dissimilarities_df$i), pattern = "Occurrence.", replacement = "[O]")
Dissimilarities_df$i <- gsub(as.character(Dissimilarities_df$i), pattern = "Abundance.", replacement = "[A]")
Dissimilarities_df$i <- gsub(as.character(Dissimilarities_df$i), pattern = "Performance.", replacement = "[P]")
Dissimilarities_df$i <- gsub(as.character(Dissimilarities_df$i), pattern = "Informed", replacement = "Clim")
Dissimilarities_df$i <- gsub(as.character(Dissimilarities_df$i), pattern = "Naive", replacement = "")

OS_gg <- ggplot(Dissimilarities_df[Dissimilarities_df$j == "Realisable True Network",], 
                aes(x = Accuracy, y = factor(i, levels = cases))
) + 
  stat_halfeye() +
  # geom_violin() + 
  facet_wrap(~j, scales = "free_x") + 
  theme_bw() + labs(y = "", x = "Inference Accuracy")
OS_gg
ggsave(OS_gg, filename = file.path(Dir.Exports, "Fig_OS.png"), 
       width = 16, height = 12, units = "cm")

# DETECTION REGRESSION =========================================================
message("Detection of Positive, Negative, and Absent Associations")
Detection_ls <- pblapply(names(Inference_ls), FUN = function(SimName){
  # SimName <- names(Inference_ls)[1]
  x <- Inference_ls[[SimName]]
  
  ## matrices
  True_mat <- abs(x$mats$True$Realised)
  Mats_ls <- x$mats$True
  InferenceMats_ls <- c(x$mats$HMSC, list(COOCCUR = x$mats$COOCCUR))
  
  ## trait diff
  Traits <- aggregate(Trait ~ Species, data = x$ID_df, FUN = mean)
  # Traits$Species <- gsub(Traits$Species, pattern = "Sp_", replacement = "")
  if(length(rownames(True_mat)[rownames(True_mat) %nin% Traits$Species]) > 0){
    Traits <- rbind(Traits, 
                    data.frame(
                      Species = rownames(True_mat)[rownames(True_mat) %nin% Traits$Species],
                      Trait = NA
                    ))
  }
  EnvDiff_mat <- True_mat
  for(i in rownames(EnvDiff_mat)){
    focaltrait <- Traits$Trait[Traits$Species == i]
    alltraitsorder <- Traits$Trait[match(colnames(EnvDiff_mat), Traits$Species)]
    EnvDiff_mat[i,] <- abs(focaltrait - alltraitsorder)
  }
  EnvDiff_mat[lower.tri(EnvDiff_mat)] <- NA
  diag(EnvDiff_mat) <- NA
  
  True_mat <- sign(Mats_ls[[1]])
  Weighted_mat <- Mats_ls[[1]]
  
  mode_ls <- lapply(names(InferenceMats_ls), FUN = function(k){
    data.frame(
      CorrectPos = as.numeric((True_mat + InferenceMats_ls[[k]]) == 2), # true positive
      CorrectNeg = as.numeric((True_mat + InferenceMats_ls[[k]]) == -2), # true negative
      CorrectAbsent = as.numeric(((True_mat == 0) + (InferenceMats_ls[[k]] == 0) == 2)), # true absent
      Magnitude = abs(as.vector(Weighted_mat)),
      SignTrue = as.vector(sign(Weighted_mat)),
      SignInferred = as.vector(InferenceMats_ls[[k]]),
      EnvDiff = as.vector(EnvDiff_mat),
      Mode = k,
      Sim = SimName
    )
  })
  do.call(rbind, mode_ls)
})
names(Detection_ls) <- names(Inference_ls)
Detection_df <- do.call(rbind, Detection_ls)
Detection_df <- Detection_df[Detection_df$Mode != "COOCCUR", ]

## detection classification
message("Detection of different correct associations")
DetectionModels_df <- lapply(X = unique(Detection_df$Mode), 
                             FUN = function(Model){
                               print(paste("########", Model))
                               model_df <- Detection_df[Detection_df$Mode == Model, ]
                               
                               ## Correct Identification of Positive Associations -----------------------------
                               print("Models of Positive Associations")
                               if(file.exists(file.path(Dir.Exports, paste0("Bayes_Model_Positive_", Model,".RData")))){
                                 load(file.path(Dir.Exports, paste0("Bayes_Model_Positive_", Model,".RData")))
                               }else{
                                 Bayes_Model_Positive <- brm(formula = CorrectPos ~Magnitude * EnvDiff,
                                                             data = model_df,
                                                             family = bernoulli(link = "logit"),
                                                             warmup = nWarmup,
                                                             iter = nSamples,
                                                             chains = nChains,
                                                             cores = nChains,
                                                             seed = 42)
                                 save(Bayes_Model_Positive, 
                                      file = file.path(Dir.Exports, paste0("Bayes_Model_Positive_", Model,".RData")))
                               }
                               
                               ## Correct Identification of Negative Associations -----------------------------
                               print("Models of Negative Associations")
                               if(file.exists(file.path(Dir.Exports, paste0("Bayes_Model_Negative_", Model,".RData")))){
                                 load(file.path(Dir.Exports, paste0("Bayes_Model_Negative_", Model,".RData")))
                               }else{
                                 Bayes_Model_Negative <- brm(formula = CorrectNeg ~Magnitude * EnvDiff,
                                                             data = model_df,
                                                             family = bernoulli(link = "logit"),
                                                             warmup = nWarmup,
                                                             iter = nSamples,
                                                             chains = nChains,
                                                             cores = nChains,
                                                             seed = 42)
                                 save(Bayes_Model_Negative, 
                                      file = file.path(Dir.Exports, paste0("Bayes_Model_Negative_", Model,".RData")))
                               }
                               
                               ## Correct Identification of Absent Associations -------------------------------
                               print("Models of Absent Associations")
                               if(file.exists(file.path(Dir.Exports, paste0("Bayes_Model_Absent_", Model,".RData")))){
                                 load(file.path(Dir.Exports, paste0("Bayes_Model_Absent_", Model,".RData")))
                               }else{
                                 Bayes_Model_Absent <- brm(formula = CorrectAbsent ~ Magnitude * EnvDiff,
                                                                data = model_df,
                                                                family = bernoulli(link = "logit"),
                                                                warmup = nWarmup,
                                                                iter = nSamples,
                                                                chains = nChains,
                                                                cores = nChains,
                                                                seed = 42)
                                 save(Bayes_Model_Absent,
                                      file = file.path(Dir.Exports, paste0("Bayes_Model_Absent_", Model,".RData")))
                               }
                               list(Pos = Bayes_Model_Positive,
                                    Neg = Bayes_Model_Negative,
                                    Abs = Bayes_Model_Absent)
                             })
names(DetectionModels_df) <- unique(Detection_df$Mode)

## inference classification
message("Inference of different associations (correct or not)")
InferenceModels_df <- lapply(X = unique(Detection_df$Mode), 
                             FUN = function(Model){
                               print(paste("########", Model))
                               model_df <- Detection_df[Detection_df$Mode == Model, ]
                               model_df <- model_df[!is.na(model_df$SignInferred), ]
                               model_df$Magnitude <- model_df$Magnitude*model_df$SignTrue
                               
                               ## Correct Identification of Positive Associations -----------------------------
                               print("Models of Positive Associations")
                               if(file.exists(file.path(Dir.Exports, paste0("BayesInf_Model_Positive_", Model,".RData")))){
                                 load(file.path(Dir.Exports, paste0("BayesInf_Model_Positive_", Model,".RData")))
                               }else{
                                 run_df <- model_df
                                 run_df$SignInferred[run_df$SignInferred != 1] <- 0
                                 Bayes_Model_Positive <- brm(formula = SignInferred ~ Magnitude * EnvDiff,
                                                             data = run_df,
                                                             family = bernoulli(link = "logit"),
                                                             warmup = nWarmup,
                                                             iter = nSamples,
                                                             chains = nChains,
                                                             cores = nChains,
                                                             seed = 42)
                                 save(Bayes_Model_Positive, 
                                      file = file.path(Dir.Exports, paste0("BayesInf_Model_Positive_", Model,".RData")))
                               }
                               
                               ## Correct Identification of Negative Associations -----------------------------
                               print("Models of Negative Associations")
                               if(file.exists(file.path(Dir.Exports, paste0("BayesInf_Model_Negative_", Model,".RData")))){
                                 load(file.path(Dir.Exports, paste0("BayesInf_Model_Negative_", Model,".RData")))
                               }else{
                                 run_df <- model_df
                                 run_df$SignInferred[run_df$SignInferred != -1] <- 0
                                 run_df$SignInferred <- abs(run_df$SignInferred)
                                 Bayes_Model_Negative <- brm(formula = SignInferred ~Magnitude * EnvDiff,
                                                             data = run_df,
                                                             family = bernoulli(link = "logit"),
                                                             warmup = nWarmup,
                                                             iter = nSamples,
                                                             chains = nChains,
                                                             cores = nChains,
                                                             seed = 42)
                                 save(Bayes_Model_Negative, 
                                      file = file.path(Dir.Exports, paste0("BayesInf_Model_Negative_", Model,".RData")))
                               }
                               
                               ## Correct Identification of Absent Associations -------------------------------
                               print("Models of Absent Associations")
                               if(file.exists(file.path(Dir.Exports, paste0("BayesInf_Model_Absent_", Model,".RData")))){
                                 load(file.path(Dir.Exports, paste0("BayesInf_Model_Absent_", Model,".RData")))
                               }else{
                                 run_df <- model_df
                                 run_df$SignInferred[run_df$SignInferred != 0] <- 99
                                 run_df$SignInferred[run_df$SignInferred == 0] <- 1
                                 run_df$SignInferred[run_df$SignInferred == 99] <- 0
                                 
                                 Bayes_Model_Absent <- brm(formula = SignInferred ~ Magnitude * EnvDiff,
                                                           data = run_df,
                                                           family = bernoulli(link = "logit"),
                                                           warmup = nWarmup,
                                                           iter = nSamples,
                                                           chains = nChains,
                                                           cores = nChains,
                                                           seed = 42)
                                 save(Bayes_Model_Absent,
                                      file = file.path(Dir.Exports, paste0("BayesInf_Model_Absent_", Model,".RData")))
                               }
                               list(Pos = Bayes_Model_Positive,
                                    Neg = Bayes_Model_Negative,
                                    Abs = Bayes_Model_Absent)
                             })
names(InferenceModels_df) <- unique(Detection_df$Mode)

## Plotting --------------------------------------------------------------------
message("Plotting")
FUN.BayesPlot <- function(Model_ls = DetectionModels_df,
                          which = "Occurrence.Informed",
                          colo = "Inference"){
  plot1_ls <- Model_ls[which]
  
  return_ls <- pblapply(names(plot1_ls), FUN = function(NAME_iter){
    plot_ls <- plot1_ls[[NAME_iter]]
    post_ls <- lapply(plot_ls, FUN = function(x){
      Y <- posterior_samples(x)
      Y[,-1] <- apply(Y[,-1], MARGIN = 2, FUN = function(x){
        exp(Y$b_Intercept + x) - 
          exp(Y$b_Intercept)
      })
      Y$b_Intercept <- exp(Y$b_Intercept)
      Y <- reshape2::melt(Y[,1:4])
      Y
    })
    
    ProbMat <- expand.grid(seq(from = 0, to = 10, length.out = 1e2),
                           seq(from = 0, to = 3, length.out = 1e2))
    if(colo == "Inference"){
      ProbMat <- expand.grid(seq(from = -10, to = 10, length.out = 1e2),
                             seq(from = 0, to = 3, length.out = 1e2))
    }
    
    
    colnames(ProbMat) <- c("Magnitude", "EnvDiff")
    
    prob_ls <- lapply(plot_ls, FUN = function(x){
      post_df <- posterior_samples(x)
      prob_iter <- ProbMat
      prob_iter$Prob <- apply(prob_iter, MARGIN = 1, FUN = function(x){
        mean(inv_logit(post_df$b_Intercept + 
                         post_df$b_EnvDiff * as.numeric(x[2]) + 
                         (post_df$b_Magnitude + post_df$`b_Magnitude:EnvDiff` * as.numeric(x[2])) * as.numeric(x[1])
        )
        )
      })
      prob_iter
    })
    
    plots_list <- lapply(names(post_ls), FUN = function(name_iter){
      post_plot <- post_ls[[name_iter]]
      post_plot$variable <- c("Intercept", "Magnitude", "Environmental \n Difference", "Interaction")[match(post_plot$variable, unique(post_plot$variable))]
      
      Coeff_gg <- ggplot(post_plot,
                         aes(x = value)) + 
        # geom_histogram() +
        stat_halfeyeh() +
        labs(y = "", x = "Model Coefficient") + 
        facet_wrap(~ factor(variable, levels = c("Intercept", "Magnitude", "Environmental \n Difference", "Interaction")), 
                   ncol = 2, scales = "free") + 
        geom_vline(xintercept = 0) + 
        theme_bw() 
      
      prob_plot <- prob_ls[[name_iter]]
      
      Prob_gg <- ggplot(prob_plot, aes(x = EnvDiff, y = Magnitude, fill = Prob)) + 
        geom_tile() + 
        guides(fill = guide_colourbar(barwidth = 2,
                                      barheight = 15,
                                      title = "Probability of \n Indetification")) + 
        theme_bw() + 
        labs(y = "Association Magnitude", x = "Environmental Difference")
      
      if(colo == "Inference"){
        Prob_gg <- Prob_gg + scale_fill_viridis_c(option = "A", limits = c(0,1))
      }else{
        Prob_gg <- Prob_gg + scale_fill_viridis_c(option = "E", limits = c(0,1))
      }
      cowplot::plot_grid(Coeff_gg, Prob_gg)
    })
    
    title <- ggdraw() + 
      draw_label(
        NAME_iter,
        fontface = 'bold',
        x = 0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
      )
    
    cowplot::plot_grid(title, 
                       cowplot::plot_grid(plotlist = plots_list, ncol = 1
                                          # , labels = names(post_ls)
                                          ), ncol = 1, rel_heights = c(0.1, 1)
                       )
    }
    )
  names(return_ls) <- names(plot1_ls)
  return_ls
}
oldw <- getOption("warn")
options(warn=-1)
suppressMessages({Detection_plots <- FUN.BayesPlot(Model_ls = DetectionModels_df, which = names(DetectionModels_df), colo = "Detection")})
suppressMessages({Inference_plots <- FUN.BayesPlot(Model_ls = InferenceModels_df, which = names(DetectionModels_df), colo = "Inference")})
options(warn = oldw)

pdf(file.path(Dir.Exports, "FIG_Inference.pdf"), width = 16, height = 22, onefile = TRUE)
Inference_plots
dev.off()

pdf(file.path(Dir.Exports, "FIG_Detection.pdf"), width = 16, height = 22, onefile = TRUE)
Detection_plots
dev.off()


# ERROR RATES ==================================================================
message("Inference Error Rates")
ErrorRates_ls <- pblapply(names(Inference_ls), FUN = function(SimName){
  # SimName <- names(Inference_ls)[[98]]
  x <- Inference_ls[[SimName]]
  
  ## matrices
  True_mat <- sign(x$mats$True$Realised)
  InferenceMats_ls <- c(x$mats$HMSC, list(COOCCUR = x$mats$COOCCUR))
  
  do.call(rbind, lapply(names(InferenceMats_ls),
                        function(y){
                          data.frame(
                            Values = c(
                              # true positive associations; how many of the inferred positive associations are correctly identified as such? - NaN means no positives in the inferred matrix
                              TP = sum((True_mat + InferenceMats_ls[[y]]) == 2, na.rm = TRUE)/
                                sum(InferenceMats_ls[[y]] == 1, na.rm = TRUE), 
                              # true negative associations; how many of the inferred negative associations are correctly identified as such? - NaN means no negatives in the inferred matrix
                              TN = sum((True_mat + InferenceMats_ls[[y]]) == -2, na.rm = TRUE)/
                                sum(InferenceMats_ls[[y]] == -1, na.rm = TRUE),
                              # true absent; how many of the inferred absent interactions are correctly identified as such? - NaN means no absents in the inferred matrix
                              TA = sum((True_mat == 0) + (InferenceMats_ls[[y]] == 0) == 2, na.rm = TRUE)/
                                sum(InferenceMats_ls[[y]] == 0, na.rm = TRUE),
                              
                              # falsely inferred positive; how many of the inferred positive associations are incorrectly identified as such?
                              FP = sum((True_mat != 1) + (InferenceMats_ls[[y]] == 1) == 2, na.rm = TRUE)/
                                sum(InferenceMats_ls[[y]] == 1, na.rm = TRUE), 
                              # falsely inferred negative; how many of the inferred negative associations are incorrectly identified as such?
                              FN = sum((True_mat != -1) + (InferenceMats_ls[[y]] == -1) == 2, na.rm = TRUE)/
                                sum(InferenceMats_ls[[y]] == -1, na.rm = TRUE),
                              # false absent; how many of the inferred absent interactions are correctly identified as such?
                              FA = sum((True_mat != 0) + (InferenceMats_ls[[y]] == 0) == 2, na.rm = TRUE)/
                                sum(InferenceMats_ls[[y]] == 0, na.rm = TRUE),
                              
                              # true positive links which were not inferred; how many of the true positive associations were not inferred?
                              MP = (sum((True_mat == 1) + (InferenceMats_ls[[y]] != 1) == 2, na.rm = TRUE)/
                                      sum(True_mat == 1, na.rm = TRUE)), 
                              # true negative links which were not inferred; how many of the true negative associations were not inferred?
                              MN = (sum((True_mat == -1) + (InferenceMats_ls[[y]] != -1) == 2, na.rm = TRUE)/
                                      sum(True_mat == -1, na.rm = TRUE)),
                              # true absent links which were not inferred; how many of the true absent associations were not inferred?
                              MA = (sum((True_mat == 0) + (InferenceMats_ls[[y]] != 0) == 2, na.rm = TRUE)/
                                      sum(True_mat == 0, na.rm = TRUE))
                              
                            ),
                            Metric = c("TP", "TN", "TA", "FP", "FN", "FA", "MP", "MN", "MA"),
                            # Identifier
                            Approach = y,
                            SimName = SimName
                          )
                        }))
})
ErrorRates_df <- do.call(rbind, ErrorRates_ls)
ErrorRates_df$Approach <- gsub(as.character(ErrorRates_df$Approach), pattern = "Occurrence.", replacement = "[O]")
ErrorRates_df$Approach <- gsub(as.character(ErrorRates_df$Approach), pattern = "Abundance.", replacement = "[A]")
ErrorRates_df$Approach <- gsub(as.character(ErrorRates_df$Approach), pattern = "Performance.", replacement = "[P]")
ErrorRates_df$Approach <- gsub(as.character(ErrorRates_df$Approach), pattern = "Informed", replacement = "Clim")
ErrorRates_df$Approach <- gsub(as.character(ErrorRates_df$Approach), pattern = "Naive", replacement = "")
ErrorRates_df <- ErrorRates_df[ErrorRates_df$Approach != "COOCCUR", ]
NonNA_count <- aggregate(Values ~ Approach + Metric, data = ErrorRates_df,
                         function(x) { sum(!is.na(x)) },
                         na.action = NULL)
colnames(NonNA_count)[3] <- "n"

ggplot(ErrorRates_df, aes(y = Values, x = Approach)) +
  geom_text(data = NonNA_count, aes(x = Approach, y = 1.1, label = n)) + 
  geom_boxplot() +
  facet_wrap(~ factor(Metric, levels = c("TP", "FP", "MP", "TN", "FN", "MN", "TA", "FA", "MA")),
             ncol = 3, scales = "free") +
  theme_bw() + labs(y = "Detection Rate [%]", x = "") + 
  lims(y = c(0, 1.1))
ggsave(filename = file.path(Dir.Exports, paste0("Fig_ErrorRates.png")),
       width = 32, height = 20, units = "cm")

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
    Sparcity = 0.1,
    MaxStrength = 10,
    ## Initial Individual Creation
    n_individuals = 4e2,
    n_mode = "each", # or "total"
    Env_range = c(0, 10),
    Trait_sd = 1,
    ## Carrying Capacity Creation
    k_range = c(300,300),
    ## Simulation Parameters
    d0 = 0.4,
    b0 = 0.6,
    t_max = 10,
    t_inter = 0.1,
    sd = 2.5,
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
ggsave(Input_gg, filename = file.path(Dir.Concept, "Fig_SpatialInputs.png"), 
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
ggsave(InputMatrices_gg, filename = file.path(Dir.Concept, "Fig_InputMatrices.png"), 
       width = 36, height = 20, units = "cm")

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
ggsave(AbundTime_gg, filename = file.path(Dir.Concept, "Fig_AbundanceTime.png"), 
       width = 30, height = 20, units = "cm")  

## Spatial Gradient ------------------------------------------------------------
print("Spatial Gradient Visualisation")
env.xy <- function(x = NULL, y = NULL){x}
gridseq <- seq(from = eval(Simulation_Output$Call[["Env_range"]])[1],
               to = eval(Simulation_Output$Call[["Env_range"]])[2],
               length = 1e1)
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
ggsave(Gradient_gg, filename = file.path(Dir.Concept, "Fig_SpatialGradient.png"), 
       width = 20, height = 22, units = "cm")  

## Network Inference -----------------------------------------------------------
print("Network Inference")
# source2 <- function(file, start, end, ...) {
#   file.lines <- scan(file, what=character(), skip=start-1, nlines=end-start+1, sep='\n')
#   file.lines.collapsed <- paste(file.lines, collapse='\n')
#   source(textConnection(file.lines.collapsed), ...)
# }
# source2("Inference.R", 62, 360)

models_ls <- FUN.Inference(Simulation_Output, 
                           Dir.Models = Dir.Concept,
                           ModelSave = TRUE,
                           Dir.Exports = Dir.Concept,
                           Treatment_Iter = "Concept",
                           Cores = 4)


## Network Realisation ---------------------------------------------------------
### True NonRealised ----
net_mat <- models_ls$mats$True$NonRealised
edg_df1 <- melt(net_mat)
colnames(edg_df1) <- c("Partner 1", "Partner 2", "Strength")
TrueMat_gg <- ggplot(edg_df1, aes(x = `Partner 1`, y = `Partner 2`, fill = Strength)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) + 
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Association")) + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
  scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")

### Trait Difference ----
SPTrait_df <- data.frame(Simulation_Output$Traits)
SPTrait_df$Species <- rownames(SPTrait_df)
colnames(SPTrait_df) <- c("Trait", "Species")
# SPTrait_df <- aggregate(ID_df, Trait ~ Species, FUN = mean) # might want to consider whether to delimit this by initialising or final trait values
# SPTrait_df$SD <- aggregate(ID_df, Trait ~ Species, FUN = sd)$Trait
SPTrait_df$SD <- 1
SPTrait_mat <- abs(outer(SPTrait_df$Trait, SPTrait_df$Trait, '-'))
colnames(SPTrait_mat) <- rownames(SPTrait_mat) <- SPTrait_df$Species
SPTraitSD_mat <- abs(outer(SPTrait_df$SD, SPTrait_df$SD, '+'))
colnames(SPTraitSD_mat) <- rownames(SPTraitSD_mat) <- SPTrait_df$Species

TraitDiff_mat <- SPTrait_mat#-SPTraitSD_mat
diag(TraitDiff_mat) <- NA
TraitDiff_mat[lower.tri(TraitDiff_mat)] <- NA

d2.df <- reshape2::melt(TraitDiff_mat, c("x", "y"), value.name = "z")
d2.df$Realised <- FALSE
d2.df$Realised[which(d2.df$z < Simulation_Output$Call$sd+Simulation_Output$Call$Effect_Dis)] <- TRUE

TraitDiff_gg <- ggplot(data = d2.df, aes(x = x, y = y, fill = z))+
  geom_tile(color = "black", lwd = 0.5, linetype = 1) + 
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Trait \n Difference")) + 
  scale_fill_viridis_c() +
  geom_point(aes(shape = Realised)) + 
  scale_shape_manual(values = c(26, 15)) + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
  labs(x = "Partner 1", y = "Partner 2") + guides(shape = "none") 

### True Realised ----
net_mat <- models_ls$mats$True$Realised
diag(net_mat) <- NA
colnames(net_mat) <- rownames(net_mat) <- V(Simulation_Output$Network)
True_mat <- net_mat
edg_df1 <- melt(net_mat)
colnames(edg_df1) <- c("Partner 1", "Partner 2", "Strength")
TrueMatReal_gg <- ggplot(edg_df1, aes(x = `Partner 1`, y = `Partner 2`, fill = Strength)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) + 
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Association")) + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
  scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")

### Plotting ----
ggsave(cowplot::plot_grid(TrueMat_gg, TraitDiff_gg, TrueMatReal_gg, ncol = 3, labels = "AUTO"),
       filename = file.path(Dir.Concept, "Fig_Realisation.png"), 
       width = 42, height = 11, units = "cm"
)

## Inference Visualisation -----------------------------------------------------
PrepMat <- function(matrix = models_ls$mats$COOCCUR, edg_df2 = edg_df1, name = "COOCCUR"){
  colnames(matrix) <- rownames(matrix) <- gsub(colnames(matrix), pattern = "Sp_", replacement = "")
  model_df <- melt(matrix)
  colnames(model_df) <- c("Partner 1", "Partner 2", "Strength")
  model_df$Approach <- name
  model_df$Correct <- model_df$Strength == sign(edg_df2$Strength)
  model_df
}

# COOCCUR_plot <- PrepMat()
# COOCCUR_plot <- ggplot(COOCCUR_plot, aes(x = `Partner 1`, y = `Partner 2`, fill = Strength)) +
#   geom_tile(color = "black", lwd = 0.5, linetype = 1) + 
#   coord_fixed() +
#   guides(fill = guide_colourbar(barwidth = 2,
#                                 barheight = 15,
#                                 title = "Inferred \n Association")) + 
#   geom_point(aes(shape = Correct)) + 
#   scale_shape_manual(values = c(26, 15)) + 
#   theme_bw() + 
#   theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
#   scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365") + 
#   guides(shape = "none") 

model_df <- do.call(rbind, lapply(names(models_ls$mats$HMSC), FUN = function(x){
  PrepMat(matrix = models_ls$mats$HMSC[[x]], name = x)
}))

cases <- as.character(unique(model_df$Approach))
cases <- c(
  # "COOCCUR", 
  rev(cases[startsWith(cases, "O")]),
  rev(cases[startsWith(cases, "A")]), rev(cases[startsWith(cases, "P")]))

HMSC_plot <- ggplot(model_df, aes(x = `Partner 1`, y = `Partner 2`, fill = Strength)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) + 
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Inferred \n Association")) + 
  geom_point(aes(shape = Correct)) + 
  scale_shape_manual(values = c(26, 15)) + 
  theme_bw() + 
  facet_wrap(~factor(Approach, levels = cases), ncol = 3, dir="v") +
  theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
  scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365") + 
  guides(shape = "none") 

cowplot::plot_grid(
  cowplot::plot_grid(
    TrueMatReal_gg
    # , 
    # COOCCUR_plot
    ),
  HMSC_plot + theme(plot.margin = unit(c(0,0,0,0), "cm")), 
  ncol = 1, rel_heights = c(2.5, 3), labels = "AUTO"
)

ggsave(filename = file.path(Dir.Concept, "Fig_Inference.png"), 
       width = 24, height = 29, units = "cm")
