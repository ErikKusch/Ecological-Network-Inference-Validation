#' ####################################################################### #
#' PROJECT: [PhD - Chapter 4]  
#' CONTENTS: 
#'  - Preparation of PFTC data for IF-REM application
#'  DEPENDENCIES:
#'  - None
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
rm(list = ls())

# PREAMBLE ================================================================
## Directories ------------------------------------------------------------
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "Data")
if(!dir.exists(Dir.Data)){dir.create(Dir.Data)}
Dir.Data.PFTC <- file.path(Dir.Data, "PFTC")
if(!dir.exists(Dir.Data.PFTC)){dir.create(Dir.Data.PFTC)}
Dir.Plots <- file.path(Dir.Base, "Plots")
if(!dir.exists(Dir.Plots)){dir.create(Dir.Plots)}

## Packages ---------------------------------------------------------------
### CRAN ------------------------------------------------------------------
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

## FUNCTIONALITY -----------------------------------------------------------
`%nin%` <- Negate(`%in%`)

hush <- function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

Sort.DF <- function(Data = NULL, Column = NULL){
  Data[order(Data[ , Column] ), ]
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# DATA RETRIEVAL ==========================================================
# library(igraph)
# Rand_igraph <- erdos.renyi.game(n = 10, # number of nodes
#                                 p.or.m = 8,
#                                 type = c("gnm"), # "gnp" ...  probability for drawing an edge between two arbitrary vertices; "gnm" number of edges in the graph
#                                 directed = FALSE,
#                                 loops = FALSE
#                                 )
# library(randcorr)
Rand_corr <- randcorr(20) # establish random correlation matrix
Rand_corr[lower.tri(Rand_corr)] <- NA # make into undirected adjacency matrix representation
Rand_corr[sample(which(!is.na(Rand_corr)), sum(!is.na(Rand_corr))*0.5)] <- 0

Rand_corr <- graph_from_adjacency_matrix(adjmatrix = Rand_corr,
                                         mode = "undirected",
                                         weighted = TRUE,
                                         diag = FALSE)
plot(Rand_corr)

E(Rand_corr)$weight
save(Rand_corr, file = file.path(Dir.Data, "RandomNetworkDevelBIG.RData"))
