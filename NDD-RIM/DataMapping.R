#' ####################################################################### #
#' PROJECT: [PhD; X - DATA SIMULATIONA] 
#' CONTENTS: 
#'  - Maps simulation results and transforms into data useable by NDD-RIM 
#'  DEPENDENCIES:
#'  - DataSimulations_NDD-RIM.R
#' AUTHOR: [Malyon Bimler]
#' ####################################################################### #


# PREAMBLE ================================================================
## Directories ------------------------------------------------------------
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "DataMaps")
if(!dir.exists(Dir.Data)){dir.create(Dir.Data)}

## Packages ---------------------------------------------------------------
library(stringr)
library(plotrix)  # to draw circles on plots

# TOY DATA ================================================================

# get list of simulation names - set up now but important for later
Simul.Names <- list.files('Data/', full.names = F)
Simul.Names <- stringr::str_remove_all(Simul.Names, '.RData')

# set colours for each species 
species.col <- RColorBrewer::brewer.pal(12, 'Set3')

# this can be iterated later on when we have more simulations 
load('Data/NDD-RIM_ToyData_SIM_42.RData')
# select last simulation step
datmap <- Simulation_Output$Simulation[[length(Simulation_Output$Simulation)]]

## CREATING OBSERVATIONS =================================================

# Set number of focal individuals to sample (aka observe)
n = 10   
# don't sample from the outer edge of the map
not.edge <- datmap[datmap$X > 1 & datmap$X < 9 & datmap$Y > 1 & datmap$Y < 9, ]
# focal individual ID 
foc.i <- sample(not.edge$ID, n)
foc.i.row <- rownames(datmap[datmap$ID %in% foc.i, ])
# this can later be modified to be more specific e.g. select 5 indivs of species 1 
# TO FIX: DON'T ALLOW OVERLAPPING NEIGHBOURHOODS

# get coordinates of each focal individual 
focal.coord <- datmap[datmap$ID %in% foc.i, c('X', "Y")]

# Set size of the interaction neighbourhood
size.in <- 0.5

# store observations in dataframe below
observations <- data.frame(matrix(NA, nrow = n, ncol = 4 + length(levels(as.factor(datmap$Species)))))
colnames(observations) <- c('Sim', 'focalID', 'focalSp', 'performance', levels(as.factor(datmap$Species)))

# For each individual observations, get the number of neighbours of each species 
for (f in 1:n) {
  xrange <- c(focal.coord[foc.i.row[f], 'X'] - size.in, focal.coord[foc.i.row[f], 'X'] + size.in)
  yrange <- c(focal.coord[foc.i.row[f], 'Y'] - size.in, focal.coord[foc.i.row[f], 'Y'] + size.in)
  neighbours <- datmap[datmap$X > xrange[1] & datmap$X < xrange[2] & 
                         datmap$Y > yrange[1] & datmap$Y < yrange[2], ]
  # insert results into df
  observations[f, ] <- c('ToySim',                       # Simulation ID - modify when working on multiple sims
                         foc.i[f],                       # focalID 
                         datmap[foc.i.row[f], ]$Species, # focal species 
                         NA,                             # individual performance
                         table(neighbours$Species))      # neighbour abundances
}


## MAP ====================================================================
# Add extra space to right of plot area
par(mar = c(2, 2, 2, 6)) 
# plot all individuals
plot(datmap$X, datmap$Y, 
     pch = 16, 
     col = species.col[as.factor(datmap$Species)],
     xlab = 'x coordinates', ylab = 'y coordinates', main = Simul.Names)
legend(
  "topright", 
  inset=c(-0.25,0),
  legend = paste("Color", levels(as.factor(datmap$Species))), # for readability of legend
  col = species.col,
  pch = 16, 
  cex = .7 
)
# remove the edges from the map
polygon(x = c(1, 9, 9, 1), y = c(1, 1, 9, 9), lty = 2, border = 'black')
# show focal individuals which are observed
points(datmap[foc.i.row, ]$X, datmap[foc.i.row, ]$Y, 
       pch = 21, col = 'black', bg = species.col[as.factor(datmap[foc.i.row, ]$Species)])
# and their interaction neighbourhood
mapply(draw.circle, datmap[foc.i.row, ]$X, datmap[foc.i.row, ]$Y, 0.5 , lty = 2)



