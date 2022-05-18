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
  
)
sapply(package_vec, install.load.package)

#### KrigR ------
if("KrigR" %in% rownames(installed.packages()) == FALSE){ # KrigR check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("https://github.com/ErikKusch/KrigR", force = TRUE)
}
library(KrigR) 
try(source("X - PersonalSettings.R")) # I do this here to specify number of cores and API credentials and am thus not sharing this file
# CDS API (needed for ERA5-Land downloads)
if(!exists("API_Key") | !exists("API_User")){ # CS API check: if CDS API credentials have not been specified elsewhere
  API_User <- readline(prompt = "Please enter your Climate Data Store API user number and hit ENTER.")
  API_Key <- readline(prompt = "Please enter your Climate Data Store API key number and hit ENTER.")
} # end of CDS API check
# NUMBER OF CORES
if(!exists("numberOfCores")){ # Core check: if number of cores for parallel processing has not been set yet
  numberOfCores <- as.numeric(readline(prompt = paste("How many cores do you want to allocate to these processes? Your machine has", parallel::detectCores())))
} # end of Core check

#### PhyloMaker ------
if("V.PhyloMaker" %in% rownames(installed.packages()) == FALSE){ # PhyloMaker check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("jinyizju/V.PhyloMaker")
}
library(V.PhyloMaker) 

#### Rethinking ------
if("rethinking" %in% rownames(installed.packages()) == FALSE){ # PhyloMaker check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("stan-dev/cmdstanr")
  install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))
  devtools::install_github("rmcelreath/rethinking")
}
library(rethinking)

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

Sort.DF <- function(Data = NULL, Column = NULL){
  Data[order(Data[ , Column] ), ]
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

## Sourcing ---------------------------------------------------------------
source("X - Functions_Data.R")

# DATA RETRIEVAL ==========================================================
if(!file.exists(file.path(Dir.Data.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"))){
download.file(url = "https://osf.io/hjpwt/download",
              destfile = file.path(Dir.Data.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"), mode = "wb")
}
Raw_df <- read.csv(file.path(Dir.Data.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")) # load pftc data
Raw_df$site <- paste(Raw_df$site, Raw_df$treatment, Raw_df$plot_id, sep = "_")

## Fitness Data Frame -----------------------------------------------------
Fitness_var <- "dry_mass_g" # defineby what metric we assess plant performance
Perf_df <- Raw_df[Raw_df$trait == Fitness_var, c("value", "taxon", "site")]

message("You now have the performance of each individual belonging to each species at each site in the object Perf_df.")
head(Perf_df)

## Neighbour Data Frame ---------------------------------------------------
Neigh_df <- stats::aggregate(individual_nr ~ site + taxon, data = Raw_df, FUN = max)
Neigh_df <- stats::reshape(data = Neigh_df, 
               idvar = "site", 
               timevar = "taxon",
               direction = "wide")
colnames(Neigh_df) <- gsub("individual_nr.", "", colnames(Neigh_df))
Neigh_df[is.na(Neigh_df)] <- 0

message("You now have the abundance of each species at each site in the object Neigh_df.")
head(Neigh_df)

## Climate Data Frame -----------------------------------------------------
ECV_vec <- c("2m_temperature", "volumetric_soil_water_layer_1", "total_precipitation")

if(!file.exists(file.path(Dir.Data.PFTC, "Metadata_df.csv"))){ # bioclimatic data not established yet
  if(!file.exists(file.path(Dir.Data.PFTC, "PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx"))){
    download.file(url = "https://osf.io/uk85w/download", 
                  destfile = "PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx", mode = "wb")
  }
  Metadata_df <- as.data.frame(readxl::read_xlsx("PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx")) # read metadata file, available here https://osf.io/uk85w/
  Metadata_df$SiteID <- with(Metadata_df, paste(Site, Treatment, PlotID, sep = "_")) # create ID column which combines site, treatment, and plotid
  Metadata_df$Year <- min(unique(Raw_df$year)) # set observation years to 2018 for all sites (that was the field season most of the data was collected throughout)
  colnames(Metadata_df)[5] <- "Lat"
  colnames(Metadata_df)[6] <- "Lon"
  
  ### THIS IS WHERE NEW AND DIFFERENT CLIMATE DATA NEEDS TO BE DOWNLOADED
  for(Variable_Iter in ECV_vec){
    if(startsWith(x = Variable_Iter, prefix = "total_")){
      PrecipFix <- TRUE
    }else{
      PrecipFix <- FALSE
    }
    ## Download Anylsis Data at Highest Resolution
    Era5Data_ras <- download_ERA(Variable = Variable_Iter,
                                 DateStart = "1998-01-01",
                                 DateStop = paste0(unique(Metadata_df$Year),"-12-31"),
                                 TResolution = "month",
                                 TStep = 1,
                                 API_Key = API_Key,
                                 API_User = API_User,
                                 Extent = Metadata_df,
                                 ID = "SiteID",
                                 Buffer = 0.1,
                                 Dir = Dir.Data.PFTC,
                                 FileName = paste0(Variable_Iter, "_Data"),
                                 Cores = Cores,
                                 TryDown = 42,
                                 SingularDL = TRUE,
                                 verbose = TRUE,
                                 PrecipFix = PrecipFix
    )
    ## Download Uncertainty Data
    Era5Uncert_ras <- download_ERA(Variable = Variable_Iter,
                                   DataSet = "era5",
                                   Type = "ensemble_members",
                                   DateStart = "1998-01-01",
                                   DateStop = paste0(unique(Metadata_df$Year),"-12-31"),
                                   TResolution = "month",
                                   TStep = 1,
                                   API_Key = API_Key,
                                   API_User = API_User,
                                   Extent = Metadata_df,
                                   Buffer = 0.5,
                                   ID = "SiteID",
                                   Dir = Dir.Data.PFTC,
                                   FileName = paste0(Variable_Iter, "_Uncert.nc"),
                                   Cores = parallel::detectCores(),
                                   TryDown = 42,
                                   SingularDL = TRUE,
                                   verbose = TRUE,
                                   PrecipFix = PrecipFix
    )
    Era5Uncert_ras <- stackApply(Era5Uncert_ras, rep(1:(nlayers(Era5Uncert_ras)/10), each = 10), sd)
    
    ## Extract Data to Locations
    Extract_df <- cbind(Metadata_df$Lon, Metadata_df$Lat)
    colnames(Extract_df) <- c("x", "y")
    Extract_sp <- SpatialPoints(Extract_df)
    Data_mat <- raster::extract(Era5Data_ras, Extract_sp)
    Uncert_mat <- raster::extract(Era5Uncert_ras, Extract_sp)
    
    ## Save Individual Dataframes of Time-Series
    # figuring out intervals
    Down_start <- lubridate::date("1998-01-01")
    Down_end <- lubridate::date(paste0(unique(Metadata_df$Year), "-12-31"))
    T_seq <- seq(Down_start, Down_end, by = "month")
    ## data frame and saving
    
    counter <- 1
    for(Site_Iter in Metadata_df$SiteID){
      FName <- file.path(Dir.Data.PFTC, paste0(Site_Iter, ".csv"))
      Raw_df <- data.frame(
        Date = T_seq,
        Data = Data_mat[counter, ],
        UC = Uncert_mat[counter, ]
      )
      colnames(Raw_df)[2:3] <- c(Variable_Iter, paste0(Variable_Iter, "_UC"))
      if(file.exists(FName)){ # if any of these already exists
        FName1 <- read.csv(FName)[-1] # read the already existing data frame
        Raw_df <- cbind(FName1, Raw_df)[,-ncol(FName1)-1] # append current data as columns
      }
      write.csv(Raw_df, file = FName)
      counter <- counter + 1
    }
    
    ## Save Mean Values to Original Data Source
    Metadata_df$XYZ <- NA
    colnames(Metadata_df)[ncol(Metadata_df)] <- Variable_Iter
    Metadata_df[, ncol(Metadata_df)] <- rowMeans(Data_mat)
    Metadata_df$XYZ <- NA
    colnames(Metadata_df)[ncol(Metadata_df)] <- paste0(Variable_Iter, "_SD")
    Metadata_df[, ncol(Metadata_df)] <- apply(Data_mat, 1, sd)
    Metadata_df$XYZ <- NA
    colnames(Metadata_df)[ncol(Metadata_df)] <- paste0(Variable_Iter, "_UC")
    Metadata_df[, ncol(Metadata_df)] <- rowMeans(Uncert_mat)
    
    ## save metadata
    write.csv(Metadata_df, file = file.path(Dir.Data.PFTC, "Metadata_df.csv"))
  }
}
Metadata_df <- read.csv(file.path(Dir.Data.PFTC, "Metadata_df.csv"))[-1]
Metadata_df <- Metadata_df[, -c(1:4, 8, 10)]
colnames(Metadata_df) <- gsub("siteid", "site", tolower(colnames(Metadata_df)))

message("You now have the elevation, geolocational data, and climate mean, standard deviation, and data uncertainty at each site in the object Metadata_df")
head(Metadata_df)

## Trait Similarity -------------------------------------------------------
Trait_df <- stats::aggregate(value ~ trait + taxon, 
                             data = Raw_df[Raw_df$trait != Fitness_var, ], 
                             FUN = mean, na.action = na.omit)
Trait_df <- stats::reshape(data = Trait_df, 
                          idvar = "taxon", 
                          timevar = "trait",
                          direction = "wide")
colnames(Trait_df) <- gsub("value.", "", colnames(Trait_df))

message("You now have the mean trait expression of multiple traits for each species across the entire study in the object Trait_df")
head(Trait_df)
