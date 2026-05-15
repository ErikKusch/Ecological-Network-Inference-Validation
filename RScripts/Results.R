#' ####################################################################### #
#' PROJECT: [InfVal; Post-Inference Visualisation & Analyses]
#' CONTENTS:
#'  - Compare Inferred and Known Networks
#'  - Tally Error Rates for Inferred Networks
#'  DEPENDENCIES:
#'  - Inference.R
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
message("Compiling Results")

# FUN_Matcomparison <- function(mat1, mat2) {
#   eq <- sign(mat1) == sign(mat2) # avoid to later compute this twice
#   round(sum(eq, na.rm = TRUE) / sum(!is.na(eq)) * 100, 2) # get the percentage of equal values
# }

# INFERRED VS. INFERRED ========================================================
print("INFERRED VS. INFERRED")
InfComp <- pblapply(
  Inference_ls[-1], # we ignore the first because it is the toy example and is used further down for plotting
  FUN = function(Sim) {
    # Sim <- Inference_ls[[2]]
    methods_vec <- names(Sim$Inferrences)
    mats <- lapply(methods_vec, function(x) {
      Sim$Inferrences[[x]]$Effects * Sim$Inferrences[[x]]$Sig
    })
    names(mats) <- methods_vec
    comp_df <- expand.grid(names(mats), names(mats))
    comp_df$Similarity <- apply(comp_df, MARGIN = 1, FUN = function(z) {
      Val.Accuracy(mats[[z[1]]], mats[[z[2]]])
    })
    comp_df
  }
)
InfComp_df <- as.data.frame(do.call(rbind, InfComp))

cases <- as.character(unique(InfComp_df$Var1))
# cases <- c(
#   # "COOCCUR",
#   rev(cases[startsWith(cases, "O")]),
#   rev(cases[startsWith(cases, "A")]), rev(cases[startsWith(cases, "P")]))

InfComp_df$Var1 <- factor(InfComp_df$Var1, levels = c("Cooccur", "Netassoc", "HMSC"))
InfComp_df$Var2 <- factor(InfComp_df$Var2, levels = c("Cooccur", "Netassoc", "HMSC"))

mean_df <- aggregate(Similarity ~ Var1 + Var2, InfComp_df, FUN = mean)
nums <- nrow(mean_df)
names <- levels(mean_df$Var1)
mean_df <- reshape(mean_df, idvar = "Var1", timevar = "Var2", direction = "wide")[, -1]
colnames(mean_df) <- rownames(mean_df) <- names
# mean_df <- mean_df[match(rownames(mean_df), cases), match(colnames(mean_df), rev(cases))]
diag(mean_df) <- NA
mean_df[lower.tri(mean_df)] <- NA

sd_df <- aggregate(Similarity ~ Var1 + Var2, InfComp_df, FUN = sd)
names <- levels(sd_df$Var1)
sd_df <- reshape(sd_df, idvar = "Var1", timevar = "Var2", direction = "wide")[, -1]
colnames(sd_df) <- rownames(sd_df) <- names
diag(sd_df) <- NA
sd_df[upper.tri(sd_df)] <- NA

plot_df <- rbind(
  as.data.frame(as.table(as.matrix(mean_df))),
  as.data.frame(as.table(as.matrix(sd_df)))
)
plot_df$metric <- c(rep("Mean", nums), rep("SD", nums))

plot_df$Freq <- round(plot_df$Freq, 2)

# plot_df$Var1 <- gsub(as.character(plot_df$Var1), pattern = "Occurrence.", replacement = "[O]")
# plot_df$Var1 <- gsub(as.character(plot_df$Var1), pattern = "Abundance.", replacement = "[A]")
# plot_df$Var1 <- gsub(as.character(plot_df$Var1), pattern = "Performance.", replacement = "[P]")
# plot_df$Var1 <- gsub(as.character(plot_df$Var1), pattern = "Informed", replacement = "Clim")
# plot_df$Var1 <- gsub(as.character(plot_df$Var1), pattern = "Naive", replacement = "")

# plot_df$Var2 <- gsub(as.character(plot_df$Var2), pattern = "Occurrence.", replacement = "[O]")
# plot_df$Var2 <- gsub(as.character(plot_df$Var2), pattern = "Abundance.", replacement = "[A]")
# plot_df$Var2 <- gsub(as.character(plot_df$Var2), pattern = "Performance.", replacement = "[P]")
# plot_df$Var2 <- gsub(as.character(plot_df$Var2), pattern = "Informed", replacement = "Clim")
# plot_df$Var2 <- gsub(as.character(plot_df$Var2), pattern = "Naive", replacement = "")

# cases <- gsub(as.character(cases), pattern = "Occurrence.", replacement = "[O]")
# cases <- gsub(as.character(cases), pattern = "Abundance.", replacement = "[A]")
# cases <- gsub(as.character(cases), pattern = "Performance.", replacement = "[P]")
# cases <- gsub(as.character(cases), pattern = "Informed", replacement = "Clim")
# cases <- gsub(as.character(cases), pattern = "Naive", replacement = "")

InfInf_gg <- ggplot(mapping = aes(x = factor(Var1, levels = cases), y = factor(Var2, levels = cases))) +
  geom_tile(data = plot_df[plot_df$metric == "Mean", ], aes(fill = Freq)) +
  geom_label(data = plot_df[plot_df$metric == "Mean", ], aes(label = round(Freq, 2))) +
  scale_fill_viridis_c(option = "E", na.value = "transparent", name = "Mean", position = "left") +
  new_scale_fill() +
  geom_tile(data = plot_df[plot_df$metric == "SD", ], aes(fill = Freq)) +
  geom_label(data = plot_df[plot_df$metric == "SD", ], aes(label = round(Freq, 2))) +
  scale_fill_viridis_c(option = "B", na.value = "transparent", name = "SD", position = "right") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.key.width = unit(dev.size()[1] / 25, "inches"),
    legend.key.height = unit(dev.size()[1] / 20, "inches"),
    legend.box = "horizontal"
  ) +
  labs(
    x = "", y = ""
    # , title = "Inference Similarity"
  )

ggsave(InfInf_gg,
  filename = file.path(Dir.Exports, paste0(RunName, "-Fig_InfVInf.png")),
  width = 16, height = 9, units = "cm"
)

# WHOLE-NETWORK ACCURACY =======================================================
print("WHOLE-NETWORK ACCURACY")
Accuracy_df <- do.call(rbind, pblapply(Inference_ls[-1], FUN = function(Sim) {
  # Sim <- Inference_ls[[2]]

  ## matrices
  True_mat <- Sim$True$Realised
  methods_vec <- names(Sim$Inferrences)
  mats <- lapply(methods_vec, function(x) {
    Sim$Inferrences[[x]]$Effects * Sim$Inferrences[[x]]$Sig
  })

  data.frame(
    Accuracy = sapply(mats, FUN = function(mat) {
      Val.Accuracy(True_mat, mat)
    }),
    Approach = methods_vec
  )
}))


OS_gg <- ggplot(
  Accuracy_df,
  aes(x = round(Accuracy, 2), y = factor(Approach, levels = cases))
) +
  stat_halfeye() +
  # geom_violin() +
  # facet_wrap(~j, scales = "free_x") +
  lims(x = c(0, 100)) +
  theme_bw() +
  labs(y = "", x = "Inference Accuracy [%]")
OS_gg
ggsave(OS_gg, filename = file.path(Dir.Exports, paste0(RunName, "Fig_OS.png")), width = 16, height = 12, units = "cm")


# WithinCompare <- list(c("[O]", "[O]Clim"),
#                       c("[A]", "[A]Clim"),
#                       c("[P]", "[P]Clim"))
# OSBP_All_gg <- ggplot(Dissimilarities_df[Dissimilarities_df$j == "Realisable True Network",],
#                       aes(y = Accuracy, x = factor(i, levels = cases))
# ) +
#   geom_violin() +
#   geom_boxplot(width = 0.3) +
#   stat_compare_means(comparisons = WithinCompare,
#                      label.y = max(Dissimilarities_df[Dissimilarities_df$j == "Realisable True Network","Accuracy"])+4,
#                      label = "p.signif", hide.ns = TRUE) +
#   facet_wrap(~j, scales = "free_x") +
#   theme_bw() + labs(x = "", y = "Inference Accuracy [%]")

# AcrossCompare <- list(c("[O]", "[A]"),
#                       c("[O]", "[P]"),
#                       c("[A]", "[P]"))
# OSBP_NonClim_gg <- ggplot(Dissimilarities_df[(Dissimilarities_df$j == "Realisable True Network" &
#                                                 !grepl(pattern = "Clim", Dissimilarities_df$i)),],
#                           aes(y = Accuracy, x = factor(i, levels = cases))
# ) +
#   geom_violin() +
#   geom_boxplot(width = 0.3) +
#   stat_compare_means(comparisons = AcrossCompare,
#                      label = "p.signif", hide.ns = TRUE) +
#   facet_wrap(~j, scales = "free_x") +
#   theme_bw() + labs(x = "", y = "Inference Accuracy [%]")

# AcrossCCompare <- list(c("[O]Clim", "[A]Clim"),
#                        c("[O]Clim", "[P]Clim"),
#                        c("[A]Clim", "[P]Clim"))
# OSBP_Clim_gg <- ggplot(Dissimilarities_df[(Dissimilarities_df$j == "Realisable True Network" &
#                                              grepl(pattern = "Clim", Dissimilarities_df$i)),],
#                        aes(y = Accuracy, x = factor(i, levels = cases))
# ) +
#   geom_violin() +
#   geom_boxplot(width = 0.3) +
#   stat_compare_means(comparisons = AcrossCCompare,
#                      label = "p.signif", hide.ns = TRUE) +
#   facet_wrap(~j, scales = "free_x") +
#   theme_bw() + labs(x = "", y = "Inference Accuracy [%]")

# OSBP_gg <- plot_grid(OSBP_All_gg,
#                      plot_grid(OSBP_NonClim_gg, OSBP_Clim_gg, ncol = 2, labels = c("B", "C")),
#                      ncol = 1, labels = c("A", ""))
# OSBP_gg
# ggsave(OSBP_gg, filename = file.path(Dir.Exports, paste0(RunName, "Fig_OSBP.png")),
#        width = 16*2, height = 12*2, units = "cm")

# ERROR RATES ==================================================================
print("INFERENCE ERROR RATES")
ErrorRates_ls <- pblapply(Inference_ls[-1], FUN = function(Sim) {
  # Sim <- Inference_ls[[98]]

  methods_vec <- names(Sim$Inferrences)
  mats <- lapply(methods_vec, function(x) {
    Sim$Inferrences[[x]]$Effects * Sim$Inferrences[[x]]$Sig
  })
  names(mats) <- methods_vec

  df <- do.call(rbind, lapply(methods_vec, FUN = function(method) {
    df <- Val.ErrorRates(Sim$True$Realised, mats[[method]])
    df$Approach <- method
    df
  }))
  df
})
ErrorRates_df <- do.call(rbind, ErrorRates_ls)
# ErrorRates_df <- do.call(rbind, ErrorRates_ls)
# ErrorRates_df$Approach <- gsub(as.character(ErrorRates_df$Approach), pattern = "Occurrence.", replacement = "[O]")
# ErrorRates_df$Approach <- gsub(as.character(ErrorRates_df$Approach), pattern = "Abundance.", replacement = "[A]")
# ErrorRates_df$Approach <- gsub(as.character(ErrorRates_df$Approach), pattern = "Performance.", replacement = "[P]")
# ErrorRates_df$Approach <- gsub(as.character(ErrorRates_df$Approach), pattern = "Informed", replacement = "Clim")
# ErrorRates_df$Approach <- gsub(as.character(ErrorRates_df$Approach), pattern = "Naive", replacement = "")
# ErrorRates_df <- ErrorRates_df[ErrorRates_df$Approach != "COOCCUR", ]
# ErrorRates_df <- ErrorRates_df[ErrorRates_df$Metric %nin% c("FP", "FN", "FA"), ]

ErrorRates_df$Metric <- c(
  "TP" = "True Positives",
  "FP" = "False Positives",
  "MP" = "Missed Positives",
  "TN" = "True Negatives",
  "FN" = "False Negatives",
  "MN" = "Missed Negatives",
  "TA" = "True Absent",
  "FA" = "False Absent",
  "MA" = "Missed Absent"
)[ErrorRates_df$Metric]

NonNA_count <- aggregate(Values ~ Approach + Metric,
  data = ErrorRates_df,
  function(x) {
    sum(!is.na(x))
  },
  na.action = NULL
)
colnames(NonNA_count)[3] <- "n"

ER_gg <- ggplot(ErrorRates_df, aes(y = Values, x = factor(Approach, levels = cases))) +
  geom_text(data = NonNA_count, aes(x = factor(Approach, levels = cases), y = 110, label = n)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  facet_wrap(
    ~ factor(Metric, levels = c(
      "True Positives",
      "False Positives",
      "Missed Positives",
      "True Negatives",
      "False Negatives",
      "Missed Negatives",
      "True Absent",
      "False Absent",
      "Missed Absent"
    )),
    ncol = 3
  ) +
  theme_bw() +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red") +
  labs(y = "Detection Rate [%]", x = "") +
  lims(y = c(0, 140)) +
  scale_y_continuous(breaks = seq(0, 100, 10))
ER_gg
ggsave(ER_gg,
  filename = file.path(Dir.Exports, paste0(RunName, "Fig_ErrorRates.png")),
  width = 32, height = 20, units = "cm"
)

# BAEYSIAN MODELS ==============================================================
print("BAYESIAN MODEL PREPAPARTION")
# message("Detection of Positive, Negative, and Absent Associations")
Detection_ls <- pblapply(Inference_ls[-1], FUN = function(Sim) {
  # Sim <- Inference_ls[[2]]

  ## matrices
  True_mat <- Sim$True$Realised

  methods_vec <- names(Sim$Inferrences)
  mats <- lapply(methods_vec, function(x) {
    Sim$Inferrences[[x]]$Effects * Sim$Inferrences[[x]]$Sig
  })
  names(mats) <- methods_vec

  ## trait diff; need ID_df here!!!!
  Traits <- aggregate(Trait ~ Species, data = Sim$ID_df, FUN = mean)
  # Traits$Species <- gsub(Traits$Species, pattern = "Sp_", replacement = "")
  if (length(rownames(True_mat)[rownames(True_mat) %nin% Traits$Species]) > 0) {
    Traits <- rbind(
      Traits,
      data.frame(
        Species = rownames(True_mat)[rownames(True_mat) %nin% Traits$Species],
        Trait = NA
      )
    )
  }
  EnvDiff_mat <- True_mat
  for (i in rownames(EnvDiff_mat)) {
    focaltrait <- Traits$Trait[Traits$Species == i]
    alltraitsorder <- Traits$Trait[match(colnames(EnvDiff_mat), Traits$Species)]
    EnvDiff_mat[i, ] <- abs(focaltrait - alltraitsorder)
  }
  EnvDiff_mat[lower.tri(EnvDiff_mat)] <- NA
  diag(EnvDiff_mat) <- NA

  Weighted_mat <- True_mat
  True_mat <- sign(Weighted_mat)

  mode_ls <- lapply(names(mats), FUN = function(k) {
    # print(k)
    data.frame(
      CorrectPos = as.numeric((True_mat + sign(mats[[k]])) == 2), # true positive
      CorrectNeg = as.numeric((True_mat + sign(mats[[k]])) == -2), # true negative
      CorrectAbsent = as.numeric(((True_mat == 0) + (sign(mats[[k]]) == 0) == 2)), # true absent
      MagnitudeTrue = abs(as.vector(Weighted_mat)),
      MagnitudeInferred = abs(as.vector(as.matrix(mats[[k]]))),
      SignTrue = as.vector(sign(Weighted_mat)),
      SignInferred = sign(as.vector(as.matrix(mats[[k]]))),
      EnvDiff = as.vector(EnvDiff_mat),
      Mode = k
    )
  })
  do.call(rbind, mode_ls)
})
# names(Detection_ls) <- names(Inference_ls)
Detection_df <- do.call(rbind, Detection_ls)
Detection_df <- na.omit(Detection_df)
# Detection_df <- Detection_df[Detection_df$Mode != "COOCCUR", ]
# Detection_df <- Detection_df[Detection_df$Mode %in% c(
#   "Abundance.Informed",
#   "Abundance.Naive",
#   "Occurrence.Naive"
# ), ]

## detection classification
print("DETECTION REGRESSION - likelihood of correct inference")
# message("Detection of different correct associations")
DetectionModels_df <- pblapply(
  X = unique(Detection_df$Mode),
  # cl = length(unique(Detection_df$Mode)),
  FUN = function(Model) {
    # Model <- unique(Detection_df$Mode)[1]
    print(paste("########", Model))
    model_df <- Detection_df[Detection_df$Mode == Model, ]
    model_df$UnderlyingTrue <- model_df$SignTrue * model_df$MagnitudeTrue

    ## Correct Identification of Positive Associations -----------------------------
    # print("Models of Positive Associations")
    if (file.exists(file.path(Dir.Exports, paste0(RunName, "Bayes_Model_Positive_", Model, ".RData")))) {
      load(file.path(Dir.Exports, paste0(RunName, "Bayes_Model_Positive_", Model, ".RData")))
    } else {
      Bayes_Model_Positive <- brm(
        formula = CorrectPos ~ UnderlyingTrue * EnvDiff,
        data = model_df,
        family = bernoulli(link = "logit"),
        warmup = nWarmup,
        iter = nSamples,
        chains = nChains,
        cores = nChains,
        seed = 42
      )
      save(Bayes_Model_Positive,
        file = file.path(Dir.Exports, paste0(RunName, "Bayes_Model_Positive_", Model, ".RData"))
      )
    }

    ## Correct Identification of Negative Associations -----------------------------
    # print("Models of Negative Associations")
    if (file.exists(file.path(Dir.Exports, paste0(RunName, "Bayes_Model_Negative_", Model, ".RData")))) {
      load(file.path(Dir.Exports, paste0(RunName, "Bayes_Model_Negative_", Model, ".RData")))
    } else {
      Bayes_Model_Negative <- brm(
        formula = CorrectNeg ~ UnderlyingTrue * EnvDiff,
        data = model_df,
        family = bernoulli(link = "logit"),
        warmup = nWarmup,
        iter = nSamples,
        chains = nChains,
        cores = nChains,
        seed = 42
      )
      save(Bayes_Model_Negative,
        file = file.path(Dir.Exports, paste0(RunName, "Bayes_Model_Negative_", Model, ".RData"))
      )
    }

    # ## Correct Identification of Absent Associations -------------------------------
    # # print("Models of Absent Associations")
    # if(file.exists(file.path(Dir.Exports, paste0(RunName, "Bayes_Model_Absent_", Model,".RData")))){
    #   load(file.path(Dir.Exports, paste0(RunName, "Bayes_Model_Absent_", Model,".RData")))
    # }else{
    #   Bayes_Model_Absent <- brm(formula = CorrectAbsent ~ Magnitude * EnvDiff,
    #                             data = model_df,
    #                             family = bernoulli(link = "logit"),
    #                             warmup = nWarmup,
    #                             iter = nSamples,
    #                             chains = nChains,
    #                             cores = 1,
    #                             seed = 42)
    #   save(Bayes_Model_Absent,
    #        file = file.path(Dir.Exports, paste0(RunName, "Bayes_Model_Absent_", Model,".RData")))
    # }
    list(
      Pos = Bayes_Model_Positive,
      Neg = Bayes_Model_Negative
      # ,
      # Abs = Bayes_Model_Absent
    )
  }
)
names(DetectionModels_df) <- unique(Detection_df$Mode)

## inference classification
print("INFERENCE REGRESSION - likelihood of inferring a specific association")
InferenceModels_df <- pblapply(
  X = unique(Detection_df$Mode),
  # cl = length(unique(Detection_df$Mode)),
  FUN = function(Model) {
    print(paste("########", Model))
    model_df <- Detection_df[Detection_df$Mode == Model, ]
    model_df <- model_df[!is.na(model_df$SignInferred), ]
    model_df$UnderlyingTrue <- model_df$SignTrue * model_df$MagnitudeTrue

    ## Correct Identification of Positive Associations -----------------------------
    # print("Models of Positive Associations")
    if (file.exists(file.path(Dir.Exports, paste0(RunName, "BayesInf_Model_Positive_", Model, ".RData")))) {
      load(file.path(Dir.Exports, paste0(RunName, "BayesInf_Model_Positive_", Model, ".RData")))
    } else {
      run_df <- model_df
      run_df$SignInferred[run_df$SignInferred != 1] <- 0
      Bayes_Model_Positive <- brm(
        formula = SignInferred ~ UnderlyingTrue * EnvDiff,
        data = run_df,
        family = bernoulli(link = "logit"),
        warmup = nWarmup,
        iter = nSamples,
        chains = nChains,
        cores = nChains,
        seed = 42
      )
      save(Bayes_Model_Positive,
        file = file.path(Dir.Exports, paste0(RunName, "BayesInf_Model_Positive_", Model, ".RData"))
      )
    }

    ## Correct Identification of Negative Associations -----------------------------
    # print("Models of Negative Associations")
    if (file.exists(file.path(Dir.Exports, paste0(RunName, "BayesInf_Model_Negative_", Model, ".RData")))) {
      load(file.path(Dir.Exports, paste0(RunName, "BayesInf_Model_Negative_", Model, ".RData")))
    } else {
      run_df <- model_df
      run_df$SignInferred[run_df$SignInferred != -1] <- 0
      run_df$SignInferred <- abs(run_df$SignInferred)
      Bayes_Model_Negative <- brm(
        formula = SignInferred ~ UnderlyingTrue * EnvDiff,
        data = run_df,
        family = bernoulli(link = "logit"),
        warmup = nWarmup,
        iter = nSamples,
        chains = nChains,
        cores = nChains,
        seed = 42
      )
      save(Bayes_Model_Negative,
        file = file.path(Dir.Exports, paste0(RunName, "BayesInf_Model_Negative_", Model, ".RData"))
      )
    }

    # ## Correct Identification of Absent Associations -------------------------------
    # # print("Models of Absent Associations")
    # if(file.exists(file.path(Dir.Exports, paste0(RunName, "BayesInf_Model_Absent_", Model,".RData")))){
    #   load(file.path(Dir.Exports, paste0(RunName, "BayesInf_Model_Absent_", Model,".RData")))
    # }else{
    #   run_df <- model_df
    #   run_df$SignInferred[run_df$SignInferred != 0] <- 99
    #   run_df$SignInferred[run_df$SignInferred == 0] <- 1
    #   run_df$SignInferred[run_df$SignInferred == 99] <- 0
    #
    #   Bayes_Model_Absent <- brm(formula = SignInferred ~ Magnitude * EnvDiff,
    #                             data = run_df,
    #                             family = bernoulli(link = "logit"),
    #                             warmup = nWarmup,
    #                             iter = nSamples,
    #                             chains = nChains,
    #                             cores = 1,
    #                             seed = 42)
    #   save(Bayes_Model_Absent,
    #        file = file.path(Dir.Exports, paste0(RunName, "BayesInf_Model_Absent_", Model,".RData")))
    # }
    list(
      Pos = Bayes_Model_Positive,
      Neg = Bayes_Model_Negative
      # ,
      # Abs = Bayes_Model_Absent
    )
  }
)
names(InferenceModels_df) <- unique(Detection_df$Mode)

## Plotting --------------------------------------------------------------------
print("BAYESIAN MODEL PLOTTING")
FUN.BayesPlot <- function(Model_ls = DetectionModels_df,
                          which = "Occurrence.Informed",
                          colo = "Inference") {
  plot1_ls <- Model_ls[which]

  return_ls <- pblapply(names(plot1_ls), FUN = function(NAME_iter) {
    # NAME_iter <- names(plot1_ls)[1]
    plot_ls <- plot1_ls[[NAME_iter]]
    post_ls <- lapply(plot_ls, FUN = function(x) {
      Y <- posterior_samples(x)
      Y[, -1] <- apply(Y[, -1], MARGIN = 2, FUN = function(x) {
        exp(Y$b_Intercept + x) -
          exp(Y$b_Intercept)
      })
      Y$b_Intercept <- exp(Y$b_Intercept)
      Y <- reshape2::melt(Y[, 1:4])
      Y
    })

    ProbMat <- expand.grid(
      seq(from = -1, to = 1, length.out = 1e2),
      seq(from = 0, to = Env_sd + Effect_Dis, length.out = 1e2)
    )
    if (colo == "Inference") {
      ProbMat <- expand.grid(
        seq(from = -1, to = 1, length.out = 1e2),
        seq(from = 0, to = 20, length.out = 1e2)
      )
    }
    colnames(ProbMat) <- c("Magnitude", "EnvDiff")

    prob_ls <- lapply(plot_ls, FUN = function(x) {
      # prob_iter <- add_predicted_draws(newdata = ProbMat, object = x)
      post_df <- posterior_samples(x)
      colnames(post_df) <- c("b_Intercept", "b_Magnitude", "b_EnvDiff", "b_Magnitude:EnvDiff", "Intercept", "lprior", "lp__")
      prob_iter <- ProbMat
      prob_iter$Prob <- apply(prob_iter, MARGIN = 1, FUN = function(x) {
        mean(inv_logit(post_df$b_Intercept +
          post_df$b_EnvDiff * as.numeric(x[2]) +
          (post_df$b_Magnitude + post_df$`b_Magnitude:EnvDiff` * as.numeric(x[2])) * as.numeric(x[1])))
      })
      prob_iter
    })

    plots_list <- lapply(names(post_ls), FUN = function(name_iter) {
      post_plot <- post_ls[[name_iter]]
      post_plot$variable <- c("Intercept", "Magnitude", "Environmental \n Difference", "Interaction")[match(post_plot$variable, unique(post_plot$variable))]

      Coeff_gg <- ggplot(
        post_plot,
        aes(x = value)
      ) +
        # geom_histogram() +
        stat_halfeyeh() +
        labs(y = "", x = "Model Coefficient") +
        facet_wrap(~ factor(variable, levels = c("Intercept", "Magnitude", "Environmental \n Difference", "Interaction")),
          ncol = 2, scales = "free"
        ) +
        geom_vline(xintercept = 0) +
        theme_bw()

      prob_plot <- prob_ls[[name_iter]]

      Prob_gg <- ggplot(prob_plot, aes(x = EnvDiff, y = Magnitude, fill = Prob)) +
        geom_tile() +
        guides(fill = guide_colourbar(
          barwidth = 2,
          barheight = 15,
          title = "Probability of \n Indetification"
        )) +
        theme_bw() +
        labs(y = "Association Magnitude", x = "Environmental Difference")

      if (colo == "Inference") {
        Prob_gg <- Prob_gg + scale_fill_viridis_c(option = "A", limits = c(0, 1))
      } else {
        Prob_gg <- Prob_gg + scale_fill_viridis_c(option = "E", limits = c(0, 1))
      }
      cowplot::plot_grid(Coeff_gg, Prob_gg)
    })

    title <- ggdraw() +
      draw_label(
        NAME_iter,
        fontface = "bold",
        x = 0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
      )

    cowplot::plot_grid(title,
      cowplot::plot_grid(
        plotlist = plots_list, ncol = 1
        # , labels = names(post_ls)
      ),
      ncol = 1, rel_heights = c(0.1, 1)
    )
  })
  names(return_ls) <- names(plot1_ls)
  return_ls
}
oldw <- getOption("warn")
options(warn = -1)
suppressMessages({
  Detection_plots <- FUN.BayesPlot(Model_ls = DetectionModels_df, which = names(DetectionModels_df), colo = "Detection")
})
suppressMessages({
  Inference_plots <- FUN.BayesPlot(Model_ls = InferenceModels_df, which = names(InferenceModels_df), colo = "Inference")
})
options(warn = oldw)

pdf(file.path(Dir.Exports, paste0(RunName, "FIG_Detection.pdf")), width = 16, height = 22 * 2 / 3, onefile = TRUE)
print(Detection_plots)
dev.off()

pdf(file.path(Dir.Exports, paste0(RunName, "FIG_Inference.pdf")), width = 16, height = 22 * 2 / 3, onefile = TRUE)
print(Inference_plots)
dev.off()


# CONCEPT VISUALISATION ========================================================
print("CONCEPTUAL VISUALISATION")
# source("SimulationFrameworkFunctions.R")

Treatment_Iter <- Data_fs[1]
load(file.path(Dir.Data, paste0(paste0(strsplit(tools::file_path_sans_ext(Treatment_Iter), split = "_")[[1]][[1]], "_DEMO"), "_Environment.RData"))) # loads "Env_mat"
load(file.path(Dir.Data, Treatment_Iter)) # loads list objects "SimResult", "Network_igraph", "CarryingK_vec", "Niches_vec"

## Input Data Visualisation ----------------------------------------------------
print("Input Data Visualisation")

### starting constellation ----f
first_df <- SimResult[[1]]
first_df$Species <- factor(gsub(first_df$Species, pattern = "Sp_", replacement = ""))
Traits_df <- data.frame(
  Traits = Niches_vec,
  Species = names(Niches_vec)
)
Traits_df$Species <- factor(gsub(Traits_df$Species,
  pattern = "Sp_", replacement = ""
))

Initial_gg <- ggplot(
  first_df,
  aes(x = X, y = Y, col = Species, shape = Species)
) +
  geom_point() +
  scale_shape_manual(values = 1:nlevels(first_df$Species)) +
  scale_color_viridis_d() +
  geom_vline(
    data = Traits_df,
    aes(
      xintercept = Traits,
      col = Species
    )
  ) +
  xlim(
    as.numeric(head(colnames(Env_mat), 1)),
    as.numeric(tail(colnames(Env_mat), 1))
  ) +
  theme_bw()
Initial_gg

### ending constellation ----
GridCoords <- seq(
  from = as.numeric(head(colnames(Env_mat), 1)),
  to = as.numeric(tail(colnames(Env_mat), 1)),
  length = n_Grid + 1
) # [-(n_Grid + 1)]
ID_df <- SimResult[[length(SimResult)]]
last_df <- SimResult[[length(SimResult)]]
last_df$Species <- factor(gsub(last_df$Species, pattern = "Sp_", replacement = ""))
ID2_df <- SimResult[[length(SimResult) - 1]]
Traits_df <- aggregate(Trait ~ Species, data = ID_df, FUN = mean)
Traits_df$Species <- factor(gsub(Traits_df$Species,
  pattern = "Sp_", replacement = ""
))
Final_gg <- ggplot(
  last_df,
  aes(x = X, y = Y, col = Species, shape = Species)
) +
  geom_point() +
  scale_shape_manual(values = 1:nlevels(last_df$Species)) +
  scale_color_viridis_d() +
  # geom_vline(data = Traits_df,
  #            aes(xintercept = Trait,
  #                col = Species)) +
  xlim(
    as.numeric(head(colnames(Env_mat), 1)),
    as.numeric(tail(colnames(Env_mat), 1))
  ) +
  geom_vline(xintercept = GridCoords) +
  geom_hline(yintercept = GridCoords) +
  theme_bw()
Final_gg
leg <- get_legend(Final_gg + theme(legend.position = "bottom"))

### constellation plotting ----
Input_gg <- plot_grid(
  plot_grid(Initial_gg + theme(legend.position = "none"),
    Final_gg + theme(legend.position = "none"),
    nrow = 1, labels = "AUTO"
  ),
  leg,
  rel_heights = c(1, 0.1), ncol = 1
)
Input_gg
# ggsave(Input_gg, filename = file.path(Dir.Concept, "Fig_SpatialInputs.png"),
#        width = 30, height = 16, units = "cm")

### species-site matrices ----
## make locational data into site X species matrix
grids_df <- expand.grid(GridCoords, GridCoords)
colnames(grids_df) <- c("X", "Y")
grids_df$GridID <- 1:nrow(grids_df)

PoptabStore_ls <- pblapply(list(ID_df, ID2_df), FUN = function(ID_df) {
  GridsID_vec <- lapply(1:nrow(ID_df),
    FUN = function(x) {
      Xs <- which(ID_df[x, "X"] >= grids_df$X)
      Ys <- which(ID_df[x, "Y"] >= grids_df$Y)
      Xs[tail(which(Xs %in% Ys), 1)]
    }
  )
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
  ### storing site X species matrix
  Y <- PoptabStore[, -1] # rownames are grid IDs
})
Y <- PoptabStore_ls[[1]] # last simulation step

### Plotting
grid_vis <- ggplot(
  last_df,
  aes(x = X, y = Y, col = Species, shape = Species)
) +
  geom_point() +
  scale_shape_manual(values = 1:nlevels(last_df$Species)) +
  scale_color_viridis_d() +
  geom_vline(xintercept = GridCoords) +
  geom_hline(yintercept = GridCoords) +
  xlim(
    as.numeric(head(colnames(Env_mat), 1)),
    as.numeric(tail(colnames(Env_mat), 1))
  ) +
  theme_bw()
grid_vis

InputGrids_ls <- pblapply(c("Occurrence", "Abundance", "Performance"), FUN = function(Iter) {
  # print(Iter)
  Y_Iter <- Y
  if (Iter == "Occurrence") {
    Y_Iter <- sign(Y)
  }
  if (Iter == "Performance") {
    Y_Iter <- Y <- PoptabStore_ls[[2]] - PoptabStore_ls[[1]] # second-to-last-last
  }
  # Y_Iter[Y_Iter == 0] <- NA
  SPSite_df <- cbind(Y_Iter, (grids_df + diff(GridCoords)[1] / 2)[, 1:2])
  loop_vec <- colnames(SPSite_df)[-((ncol(SPSite_df) - 1):ncol(SPSite_df))]
  plot_col <- viridis_pal()(length(loop_vec))
  names(plot_col) <- loop_vec
  plot_ls <- as.list(rep(NA, length = (length(loop_vec) + 1)))
  plot_ls[[1]] <- grid_vis
  names(plot_ls) <- c("gridvis", loop_vec)
  for (i in loop_vec) {
    abund_df <- SPSite_df[, c(i, "X", "Y")]
    colnames(abund_df)[1] <- "Abund"
    if (Iter == "Occurrence") {
      abund_df$Abund <- factor(abund_df$Abund)
    }
    plot_ls[[i]] <- ggplot(abund_df, aes(x = X, y = Y, fill = Abund)) +
      geom_tile(color = "black", lwd = 0.5, linetype = 1) +
      coord_fixed() +
      theme_bw() +
      theme(legend.position = "bottom") + # "none"
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      labs(y = "", x = "")
    if (Iter == "Occurrence") {
      plot_ls[[i]] <- plot_ls[[i]] + scale_fill_manual(
        values = c(
          "white",
          as.character(plot_col[names(plot_col) == i])
        ),
        name = ""
      )
    }
    if (Iter == "Abundance") {
      plot_ls[[i]] <- plot_ls[[i]] + scale_fill_gradient(
        low = "white",
        high = plot_col[names(plot_col) == i]
      ) +
        guides(fill = guide_colourbar(
          barwidth = 8,
          barheight = 0.75,
          title = ""
        ))
    }
    if (Iter == "Performance") {
      plot_ls[[i]] <- plot_ls[[i]] + scale_fill_gradient2(
        low = "darkred",
        high = plot_col[names(plot_col) == i]
      ) +
        guides(fill = guide_colourbar(
          barwidth = 8,
          barheight = 0.75,
          title = ""
        ))
    }
    plot_ls[[i]] <- plot_ls[[i]] +
      theme(legend.box.spacing = unit(0, "pt")) # The spacing between the plotting area and the legend box (unit)
  }
  plot_grid(plotlist = plot_ls[-1], nrow = 1)
})

InputMatrices_gg <- plot_grid(Input_gg, # grid_vis
  plot_grid(plotlist = InputGrids_ls[-3], nrow = 2, labels = c("C", "D")),
  ncol = 1, rel_heights = c(2.15, 3), labels = c("A", "")
)
InputMatrices_gg
ggsave(InputMatrices_gg,
  filename = file.path(Dir.Concept, paste0(RunName, "_Fig_InputMatrices.png")),
  width = 36, height = 32, units = "cm"
)

## Abundance Through Time ------------------------------------------------------
print("Abundance Visualisation through Time")
Abund_time <- pblapply(names(SimResult),
  FUN = function(t) {
    ID_iter <- SimResult[[t]]
    cbind(data.frame(table(ID_iter$Species)), t)
  }
)
Abund_time <- do.call(rbind, Abund_time)
Abund_time$t <- as.numeric(Abund_time$t)
colnames(Abund_time)[1:2] <- c("Species", "Abundance")
Abund_time$Species <- factor(as.numeric(gsub(Abund_time$Species,
  pattern = "Sp_", replacement = ""
)))

AbundTime_gg <- ggplot(Abund_time, aes(x = t, y = Abundance, col = Species)) +
  geom_line(size = 1.5) +
  scale_color_viridis_d() +
  theme_bw()
ggsave(AbundTime_gg,
  filename = file.path(Dir.Concept, paste0(RunName, "_Fig_AbundanceTime.png")),
  width = 30, height = 20, units = "cm"
)

## Spatial Gradient ------------------------------------------------------------
print("Spatial Gradient Visualisation")
env.xy <- function(x = NULL, y = NULL) {
  x
}
gridseq <- seq(
  from = as.numeric(head(colnames(Env_mat), 1)),
  to = as.numeric(tail(colnames(Env_mat), 1)),
  length = 1e2
)
gridmat <- expand.grid(gridseq, gridseq)
colnames(gridmat) <- c("x", "y")
gridmat$Phenotype <- apply(gridmat, MARGIN = 1, FUN = function(k) {
  env.xy(x = k[1], y = k[2])
})

Gradient_gg <- ggplot(gridmat, aes(x = x, y = y, fill = Phenotype)) +
  geom_tile() +
  coord_fixed() +
  guides(fill = guide_colourbar(
    barwidth = 2,
    barheight = 15,
    title = "Optimal \n Phenotype"
  )) +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
ggsave(Gradient_gg,
  filename = file.path(Dir.Concept, paste0(RunName, "_Fig_SpatialGradient.png")),
  width = 20, height = 22, units = "cm"
)

## Network Inference -----------------------------------------------------------
print("Network Inference")

## Network Realisation ---------------------------------------------------------
### True NonRealised ----
net_mat <- Inference_ls[[1]]$True
edg_df1 <- melt(net_mat)
colnames(edg_df1) <- c("Partner 1", "Partner 2", "Strength")
TrueMat_gg <- ggplot(edg_df1, aes(x = `Partner 1`, y = `Partner 2`, fill = Strength)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) +
  coord_fixed() +
  guides(fill = guide_colourbar(
    barwidth = 2,
    barheight = 15,
    title = "Association"
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -20, hjust = 0)) +
  scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")

### Trait Difference ----
SPTrait_df <- data.frame(Niches_vec)
SPTrait_df$Species <- rownames(SPTrait_df)
colnames(SPTrait_df) <- c("Trait", "Species")
# SPTrait_df <- aggregate(ID_df, Trait ~ Species, FUN = mean) # might want to consider whether to delimit this by initialising or final trait values
# SPTrait_df$SD <- aggregate(ID_df, Trait ~ Species, FUN = sd)$Trait
SPTrait_df$SD <- Trait_sd
SPTrait_mat <- abs(outer(SPTrait_df$Trait, SPTrait_df$Trait, "-"))
colnames(SPTrait_mat) <- rownames(SPTrait_mat) <- SPTrait_df$Species
SPTraitSD_mat <- abs(outer(SPTrait_df$SD, SPTrait_df$SD, "+"))
colnames(SPTraitSD_mat) <- rownames(SPTraitSD_mat) <- SPTrait_df$Species

TraitDiff_mat <- SPTrait_mat #-SPTraitSD_mat
diag(TraitDiff_mat) <- NA
TraitDiff_mat[lower.tri(TraitDiff_mat)] <- NA

d2.df <- reshape2::melt(TraitDiff_mat, c("x", "y"), value.name = "z")
d2.df$Realised <- FALSE
d2.df$Realised[which(d2.df$z <= Env_sd + Effect_Dis)] <- TRUE

TraitDiff_gg <- ggplot(data = d2.df, aes(x = x, y = y, fill = z)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  coord_fixed() +
  guides(fill = guide_colourbar(
    barwidth = 2,
    barheight = 15,
    title = "Trait \n Difference"
  )) +
  scale_fill_viridis_c() +
  geom_point(aes(shape = Realised), col = "white", size = 5) +
  scale_shape_manual(values = c(26, 15)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -20, hjust = 0)) +
  labs(x = "Partner 1", y = "Partner 2") +
  guides(shape = "none")

### True Realised ----
net_mat <- Inference_ls[[1]]$True$Realised
diag(net_mat) <- NA
# True_mat <- net_mat
edg_df1 <- melt(net_mat)
colnames(edg_df1) <- c("Partner 1", "Partner 2", "Strength")
TrueMatReal_gg <- ggplot(edg_df1, aes(x = `Partner 1`, y = `Partner 2`, fill = Strength)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) +
  coord_fixed() +
  guides(fill = guide_colourbar(
    barwidth = 2,
    barheight = 15,
    title = "Association"
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -20, hjust = 0)) +
  scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")

### Plotting ----
ggsave(cowplot::plot_grid(TrueMat_gg, TraitDiff_gg, TrueMatReal_gg, ncol = 3, labels = "AUTO"),
  filename = file.path(Dir.Concept, paste0(RunName, "_Fig_Realisation.png")),
  width = 42, height = 11, units = "cm"
)

## Inference Visualisation -----------------------------------------------------
PrepMat <- function(matrix = models_ls$mats$COOCCUR, edg_df2 = edg_df1, name = "COOCCUR") {
  matrix <- as.matrix(matrix$Effects * matrix$Sig)
  class(matrix) <- "matrix"
  # colnames(matrix) <- rownames(matrix) <- gsub(colnames(matrix), pattern = "Sp_", replacement = "")
  model_df <- melt(matrix)
  colnames(model_df) <- c("Partner 1", "Partner 2", "Strength")
  model_df$Strength <- sign(model_df$Strength)
  model_df$Approach <- name
  model_df$Correct <- sign(model_df$Strength) == sign(edg_df2$Strength)
  model_df
}

model_df <- do.call(rbind, lapply(names(Inference_ls[[1]]$Inferrences), FUN = function(x) {
  # print(x)
  PrepMat(matrix = Inference_ls[[1]]$Inferrences[[x]], edg_df2 = edg_df1, name = x)
}))

cases <- as.character(unique(model_df$Approach))
# cases <- c(
#   # "COOCCUR",
#   rev(cases[startsWith(cases, "O")]),
#   rev(cases[startsWith(cases, "A")]), rev(cases[startsWith(cases, "P")])
# )

model_df$Strength <- factor(model_df$Strength, levels = c(-1, 0, 1))
Inference_plot <- ggplot(model_df, aes(x = `Partner 1`, y = `Partner 2`, fill = Strength)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) +
  coord_fixed() +
  guides(fill = guide_colourbar(
    barwidth = 2,
    barheight = 15,
    title = "Inferred \n Association"
  )) +
  geom_point(aes(shape = Correct)) +
  scale_shape_manual(values = c(26, 15)) +
  theme_bw() +
  facet_wrap(~ factor(Approach, levels = cases), ncol = 3, dir = "v") +
  theme(axis.text.x = element_text(angle = -20, hjust = 0)) +
  scale_fill_manual(
    values = c("-1" = "#5ab4ac", "0" = "white", "1" = "#d8b365"),
    name = "Association"
  ) +
  # scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365") +
  guides(shape = "none")

p <- cowplot::plot_grid(
  cowplot::plot_grid(
    TrueMatReal_gg
    # ,
    # COOCCUR_plot
  ),
  Inference_plot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
  ncol = 1, rel_heights = c(2.5, 3), labels = "AUTO"
)

ggsave(
  p,
  filename = file.path(Dir.Concept, paste0(RunName, "_Fig_Inference.png")),
  width = 24, height = 29, units = "cm"
)
