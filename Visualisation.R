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
Dir.Concept <- file.path(Dir.Base, "Concept")
Dirs <- c(Dir.Concept)
CreateDir <- sapply(Dirs, function(x) if(!dir.exists(x)) dir.create(x))

## Packages ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "Hmsc", # for HMSC models
  "igraph", # for graph representation
  "devtools", # to install betalink
  "ggplot2", # for plotting
  "tidybayes", # for plotting
  "reshape2", # for data reformatting
  "scales" # for gradient in network matrices
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

# SIMULATION FUNCTION ======================================================
source("SimulationFrameworkFunctions.R")

if(file.exists(file.path(Dir.Concept, "ConceptSim.RData"))){
  load(file.path(Dir.Concept, "ConceptSim.RData"))
}else{
  Simulation_Output <- FUN.SimulationFramework(
    seed = 42,
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
  save(Simulation_Output, file = file.path(Dir.Concept, "ConceptSim.RData")) 
}

# SPATIAL & HMSC INPUT VISUALISATION =======================================
ID_df <- Simulation_Output$Simulation[[length(Simulation_Output$Simulation)]]

### DATA PREPRATION ####
# Y: Site X Species matrix
## make locational data into site X species matrix
GridCoords <- seq(from = eval(Simulation_Output$Call[["Env_range"]])[1], 
                  to = eval(Simulation_Output$Call[["Env_range"]])[2], 
                  length = n_Grid+1)[-(n_Grid+1)]
# grids_df <- expand.grid(GridCoords, GridCoords)
# colnames(grids_df) <- c("X", "Y")
grids_df <- data.frame(X = GridCoords,
                       Y = 0)
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


Traits_df <- data.frame(Traits = Simulation_Output$Traits,
                        Species = names(Simulation_Output$Traits))
Traits_df$Species <- factor(gsub(Traits_df$Species, pattern = "Sp_", replacement = ""))

# starting constellation
first_df <- Simulation_Output$Simulation[[1]]
first_df$Species <- factor(as.numeric(gsub(first_df$Species, pattern = "Sp_", replacement = "")))
Initial_gg <- ggplot(first_df, 
                     aes(x = X, y = Y, col = Species, shape = Species)) + 
  geom_point() + 
  scale_shape_manual(values=1:nlevels(first_df$Species)) + 
  geom_vline(data = Traits_df, 
             aes(xintercept = Traits, 
                 col = Species)) + 
  xlim(eval(Simulation_Output$Call[["Env_range"]])[1],
       eval(Simulation_Output$Call[["Env_range"]])[2]) + 
  theme_bw()

# ending constellation
last_df <- Simulation_Output$Simulation[[length(Simulation_Output$Simulation)]]
last_df$Species <- factor(as.numeric(gsub(last_df$Species, pattern = "Sp_", replacement = "")))
Traits_df <- aggregate(Trait ~ Species, data = ID_df, FUN = mean)
Traits_df$Species <- factor(as.numeric(gsub(Traits_df$Species, pattern = "Sp_", replacement = "")))
Final_gg <- ggplot(last_df, 
                   aes(x = X, y = Y, col = Species, shape = Species)) + 
  geom_point() + 
  scale_shape_manual(values=1:nlevels(last_df$Species)) + 
  geom_vline(data = Traits_df, 
             aes(xintercept = Trait, 
                 col = Species)) + 
  xlim(eval(Simulation_Output$Call[["Env_range"]])[1],
       eval(Simulation_Output$Call[["Env_range"]])[2]) +
  theme_bw()

# transect abundances
Abund_long <- reshape(PoptabStore, 
                      direction = "long",
                      varying = list(names(PoptabStore)[-1]),
                      v.names = "Abundances",
                      idvar = c("GridsID"),
                      timevar = "Species",
                      times = names(PoptabStore)[-1])
Abund_long$Coord <- rep(GridCoords, length(unique(Abund_long$Species)))+
  (GridCoords[2]-GridCoords[1])/2
Abund_long$Species <- factor(as.numeric(gsub(Abund_long$Species, pattern = "Sp_", replacement = "")))

Abund_gg <- ggplot(Abund_long, aes(x = Coord, y = Abundances, fill = Species)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_bw()

# abundances through time
Abund_time <- pblapply(names(Simulation_Output$Simulation), 
                       FUN = function(t){
                         ID_iter <- Simulation_Output$Simulation[[t]]
                         cbind(data.frame(table(ID_iter$Species)), t)
                       })
Abund_time <- do.call(rbind, Abund_time)
Abund_time$t <- as.numeric(Abund_time$t)
colnames(Abund_time)[1:2] <- c("Species", "Abundance")

AbundTime_gg <- ggplot(Abund_time, aes(x = t, y = Abundance, col = Species)) + 
  geom_line() + 
  # ylim(0,250) +
  theme_bw()

# # spatial gradient
# message("spatial gradient visualisation")

# HMSC INFERENCE & NETWORKS ================================================
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
FullMod <- Hmsc(Y = Y, XData = X,  XFormula = XFormula,
                TrData = Tr, TrFormula = TrFormula,
                distr = "poisson",
                studyDesign = studyDesign,
                ranLevels = {list("GridID" = rL.site)}
)

# message("Modelling")
hmsc_model <- FullMod
hmsc_model <- sampleMcmc(hmsc_model, samples = nSamples, thin = thin,
                         transient = nWarmup,
                         nChains = nChains,
                         nParallel = nChains
)

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
Interactions_HMSC$Sig <- FALSE
Interactions_HMSC$Sig[
  Interactions_HMSC$Inter_ProbPos >= 0.95 | 
    Interactions_HMSC$Inter_ProbNeg >= 0.95] <- TRUE
HMSC_ig <- graph_from_data_frame(Interactions_HMSC[Interactions_HMSC$Sig == TRUE,], 
                                 directed = FALSE)
Spec_vec <- unique(c(Interactions_HMSC$Partner1, Interactions_HMSC$Partner2))
if(length(E(HMSC_ig)) != 0){
  E(HMSC_ig)$weight <- E(HMSC_ig)$Inter_mean
  # E(HMSC_ig)$weight[E(HMSC_ig)$weight > 0] <- 1
  # E(HMSC_ig)$weight[E(HMSC_ig)$weight < 0] <- -1
  origvert <- names(V(HMSC_ig))
  addvert <- Spec_vec[Spec_vec %nin% names(V(HMSC_ig))
  ]
  HMSC_ig <- add_vertices(HMSC_ig, 
                          nv = length(addvert))
  HMSC_ig <- set.vertex.attribute(HMSC_ig, "name", value = c(origvert, addvert))
}else{
  HMSC_ig <- make_empty_graph(n = length(Spec_vec))
  HMSC_ig <- set.vertex.attribute(HMSC_ig, "name", value = Spec_vec)
}

message("Input network matrix")
net_mat <- as_adjacency_matrix(Simulation_Output$Network, attr = "weight",
                               type = "upper", sparse = FALSE)
net_mat[lower.tri(net_mat)] <- NA
diag(net_mat) <- NA
colnames(net_mat) <- rownames(net_mat) <- V(Simulation_Output$Network)
edg_df <- melt(net_mat)
colnames(edg_df) <- c("Partner 1", "Partner 2", "Strength")
ggplot(edg_df, aes(x = `Partner 1`, y = `Partner 2`, fill = Strength)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) + 
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Associatiuon")) + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
  scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")


message("Inferred network matrix")
net_mat <- as_adjacency_matrix(HMSC_ig, attr = "weight",
                               type = "upper", sparse = FALSE)
net_mat[lower.tri(net_mat)] <- NA
diag(net_mat) <- NA
colnames(net_mat) <- rownames(net_mat) <- V(HMSC_ig)
edg_df <- melt(net_mat)
colnames(edg_df) <- c("Partner 1", "Partner 2", "Strength")
ggplot(edg_df, aes(x = `Partner 1`, y = `Partner 2`, fill = Strength)) +
  geom_tile(color = "black", lwd = 0.5, linetype = 1) + 
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Associatiuon")) + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle = -20, hjust = 0)) + 
  scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365")











































