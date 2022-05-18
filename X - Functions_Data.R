#' ####################################################################### #
#' PROJECT: [PhD; X - DATA FUNCTIONALITY] 
#' CONTENTS: 
#'  - Functionality for retrieval of environmental data for observations
#'  - Functionality for phylogenetic distance calculation for differing levels of taxonomic resolution
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# CLIMATE DATA =============================================================
FUN.ClimData <- function(Variable = "2m_temperature",
                         Data = Metadata_df,
                         ID = "ID",
                         Year = "Year",
                         Cores = 1,
                         Shape = NULL,
                         Dir = getwd()){
  
  # Download Identification ----
  Data$Download <- 2 # create download query column, 2 indicates insufficient metadata data
  Data$Download[rowSums(!is.na(Data[,c("Lat","Lon", Year)])) == 3] <- 1 # 1 indicates sufficient metadata availability
  Data$Download[which(as.numeric(as.character(Data[,Year])) < 1981 & Data$Download == 1)] <- 3 # 3 indicates study end date predating first available year of era5-land data availability
  if(!is.null(Shape)){
    DownloadQueries <- which(Data$Download == 1)
    CoordinateCheck_df <- Data[DownloadQueries,5:6]
    OverCheck_df <- raster::extract(x = Shape, y = data.frame(CoordinateCheck_df$Lon, CoordinateCheck_df$Lat), method = "bilinear")
    OutofBounds <- which(is.na(OverCheck_df)) ## these are the points which have download query 1, but fall outside of era5-land shape
    Data$Download[DownloadQueries][OutofBounds] <- 4 # 4 indicates locations outside of land mask
  }
  DownsNeeded <- Data$Download == 1
  
  # Data Download ----
  if(length(unique(Data$Year)) == 1){ # if data set was collected in the same year
    Era5Data_ras <- download_ERA(Variable = Variable,
                                 DateStart = "1981-01-01",
                                 DateStop = paste0(unique(Data$Year),"-12-31"),
                                 TResolution = "year",
                                 TStep = unique(Data$Year)-1981+1,
                                 API_Key = API_Key,
                                 API_User = API_User,
                                 Extent = Data[DownsNeeded, ],
                                 ID = "SiteID",
                                 Buffer = 0.1,
                                 Dir = Dir.PFTC,
                                 FileName = paste0(Variable, "_Data"),
                                 Cores = Cores,
                                 TryDown = 42
    )
    Era5Uncert_ras <- download_ERA(Variable = Variable,
                                   DataSet = "era5",
                                   Type = "ensemble_members",
                                   DateStart = "1981-01-01",
                                   DateStop = paste0(unique(Data$Year),"-12-31"),
                                   TResolution = "year",
                                   TStep = unique(Data$Year)-1981+1,
                                   API_Key = API_Key,
                                   API_User = API_User,
                                   Extent = Data,
                                   Buffer = 0.5,
                                   ID = "SiteID",
                                   Dir = Dir.PFTC,
                                   FileName = paste0(Variable, "_Uncert.nc"),
                                   Cores = Cores,
                                   TryDown = 42
    )
    Era5Uncert_ras <- stackApply(Era5Uncert_ras, rep(1,10), sd)
    
    ###Extract data to locations here
    Extract_df <- cbind(Data$Lon, Data$Lat)
    colnames(Extract_df) <- c("x", "y")
    Extract_sp <- SpatialPoints(Extract_df)
    Data_vec <- raster::extract(Era5Data_ras, Extract_sp)
    Uncert_vec <- raster::extract(Era5Uncert_ras, Extract_sp)
    
    
  }else{ # if individual data points/locations come with individual end years
    for(Down_Iter in DownsNeeded){
      
    }
  }
  
  
  Data <- Data[ , -ncol(Data)]
  Data$Data <- Data_vec
  Data$Uncert <- Uncert_vec
  colnames(Data)[(ncol(Data)-1):ncol(Data)] <- paste0(Variable, c("_Data", "_Uncert"))
  return(Data)
}

# PHYLOGENETIC DISTANCE ====================================================
FUN.PhyloDist <- function(SpeciesNames = NULL, verbatim = TRUE){
  SpeciesNames <- unique(gsub(pattern = " " , replacement = "_", x = SpeciesNames)) # make sure that binary nomenclature is seperated by underscores
  Status <- rep(3, length(SpeciesNames)) # this tracks at which level data was matched with phylogeny, 3 indexes failure
  # Fully recognised species ----
  Recognised_Spec <- SpeciesNames[SpeciesNames %in% V.PhyloMaker::tips.info$species]
  Status[SpeciesNames %in% V.PhyloMaker::tips.info$species] <- 1 # match at species-level
  if(isTRUE(verbatim)){message(paste(length(Recognised_Spec), "are recognised at species level by V.Phylomaker"))}
  
  # Genus-Level recognised species ----
  Unrecognised_Spec <- SpeciesNames[SpeciesNames %nin% Recognised_Spec] ## these are the species we were given, but aren't resolved to species-level
  Reported_Gen <- unlist(lapply(strsplit(x = Unrecognised_Spec, split = "_"), '[[', 1)) # genuses of non-fully recongised species
  Known_Gen <- unlist(lapply(strsplit(x = V.PhyloMaker::tips.info$species, split = "_"), '[[', 1)) # all genuses contained within phylogeny megatree
  Recognised_Gen <- Reported_Gen[Reported_Gen %in% Known_Gen]
  Status[SpeciesNames %in% Unrecognised_Spec[Reported_Gen %in% Known_Gen]] <- 2 # genus-level match
  if(isTRUE(verbatim)){message(paste(length(Recognised_Gen), "are recognised at genus level by V.Phylomaker"))}
  
  # Non-recognised species ----
  Failed_Spec <- Unrecognised_Spec[Reported_Gen %nin% Known_Gen]
  if(isTRUE(verbatim)){message(paste(length(Failed_Spec), "are recognised neither at species or genus level by V.Phylomaker"))}
  if(length(Failed_Spec > 0) & isTRUE(verbatim)){
    message(paste("These species are:", paste(Failed_Spec, collapse = ", ")))
  }
  
  # Random Sampling of Genus-level species ----
  Phylo_Spec <- c()
  if(length(Recognised_Gen) > 0){
    Sample_Spec <- V.PhyloMaker::tips.info$species ## all species we can sample from
    Sample_Gen <- Known_Gen ## all genuses for matching
    Sample_Gen <- Known_Gen[Sample_Spec %nin% Recognised_Spec] ## remove already fully resolved species from sample spectrum
    Sample_Spec <- Sample_Spec[Sample_Spec %nin% Recognised_Spec] ## remove already fully resolved species from sample spectrum
    ## loop over all sampling genuses
    Phylo_Spec <- NA
    for(Samp_Iter in 1:length(Recognised_Gen)){ 
      if(sum(Sample_Gen %in% Recognised_Gen[Samp_Iter]) == 0){## there are not enough species left over for us to randomly sample one for such a record
        Phylo_Spec[Samp_Iter] <- NA
      }else{
        Phylo_Spec[Samp_Iter] <- sample(Sample_Spec[Sample_Gen %in% Recognised_Gen[Samp_Iter]], 1)
      }
      Sample_Gen <- Sample_Gen[Sample_Spec %nin% Phylo_Spec[Samp_Iter]]
      Sample_Spec <- Sample_Spec[Sample_Spec %nin% Phylo_Spec[Samp_Iter]]
    }
    Phylo_Spec <- c(Recognised_Spec, Phylo_Spec)
    names(Phylo_Spec) <- c(SpeciesNames[Status == 1], SpeciesNames[Status == 2])
    Phylo_Spec <- na.omit(Phylo_Spec)
  }else{
    Phylo_Spec <- Recognised_Spec
  }
  # Phylogeny Building ----
  ## prepare phylogeny building
  phylo_df <- data.frame(species = Phylo_Spec,
                         genus =  unlist(lapply(strsplit(Phylo_Spec, split = "_"), `[[`, 1)),
                         family =  V.PhyloMaker::tips.info$family[na.omit(base:: match(Phylo_Spec, V.PhyloMaker::tips.info$species))]
  )
  ## build phylogeny
  phylo_phylo <- V.PhyloMaker::phylo.maker(sp.list = phylo_df, scenarios = "S3")
  tree <- phylo_phylo$scenario.3
  tree$edge.length <-  tree$edge.length + 0.001
  
  # Phylogenetic Distance ----
  phylo_dist <- ape::cophenetic.phylo(tree) ## this looks wrong in the output
  colnames(phylo_dist) <- names(Phylo_Spec[match(colnames(phylo_dist), Phylo_Spec)])
  rownames(phylo_dist) <- names(Phylo_Spec[match(colnames(phylo_dist), Phylo_Spec)])
  
  # Export of Results -----
  return(phylo_dist) 
}