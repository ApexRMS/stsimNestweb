# Set environment variable TZ when running on AWS EC2 instance
Sys.setenv(TZ='UTC')

# Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

# Load libraries
library(rsyncrosim)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(terra))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(glmmTMB))

# Setup ----
progressBar(type = "message", message = "Preparing inputs...")

## Function definitions ----
# Define function to facilitate recoding a vector using a look-up table
lookup <- function(x, old, new){
  dplyr::recode(x, !!!set_names(new, old))
}

## Directories ----
# scriptDir <- dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])))
scriptDir <- ssimEnvironment()$PackageDirectory
spatialInputsDir <- file.path("Model-Inputs", "Spatial")
tabularDataDir <- file.path("Data", "Tabular")

## Connect to SyncroSim ----
myScenario <- scenario()

# Load relevant datasheets
# stsim
RunControl <- datasheet(myScenario, "stsim_RunControl")
Stratum <- datasheet(myScenario, "stsim_Stratum", includeKey = TRUE)
SecondaryStratum <- datasheet(myScenario, "stsim_SecondaryStratum", includeKey = TRUE)
StateClass <- datasheet(myScenario, "stsim_StateClass", includeKey = TRUE)
InitialConditionsSpatial <- datasheet(myScenario, "stsim_InitialConditionsSpatial")
OutputSpatialState <- datasheet(myScenario, "stsim_OutputSpatialState")

# stsim stock-flow
OutputSpatialStockGroup <- datasheet(myScenario, "stsim_OutputSpatialStockGroup")
StockType <- datasheet(myScenario, "stsim_StockType")
StockGroup <- datasheet(myScenario, "stsim_StockGroup")
StockTypeGroupMembership <- datasheet(myScenario, "stsim_StockTypeGroupMembership")

# stsimNetweb
SiteType <- datasheet(myScenario, "stsimNestweb_SiteType")
SpeciesID <- datasheet(myScenario, "stsimNestweb_Species", includeKey = TRUE) %>%
  pull(SpeciesId, name = Name)
Site <- datasheet(myScenario, "stsimNestweb_SiteValue") %>% slice(1)
OutputOptions <- datasheet(myScenario, "stsimNestweb_OutputOptions")
HabitatModel <- datasheet(myScenario, "stsimNestweb_HabitatModel")
InvalidHabitat <- datasheet(myScenario, "stsimNestweb_InvalidHabitat")
OutputHabitatAmount <- data.frame(
  Iteration = integer(0),
  Timestep = integer(0),
  StratumId = character(0),
  SecondaryStratumId = character(0),
  Site = character(0),
  Species = character(0),
  Amount = numeric(0)
)

OutputSpatialHabitat <- data.frame(
  Iteration = integer(0),
  Timestep = integer(0),
  Species = character(0),
  FileName = character(0)
)

OutputSpatialHabitatChange <- data.frame(
  Iteration = integer(0),
  Timestep = integer(0),
  Species = character(0),
  FileName = character(0)
)

## Setup Parameters ----
# Iterations
iterations <- seq(RunControl$MinimumIteration, RunControl$MaximumIteration)

# Timesteps
timestepsTabular <- if(OutputOptions$SummaryOutputHA){
  seq(RunControl$MinimumTimestep, RunControl$MaximumTimestep, by = OutputOptions$SummaryOutputHATimesteps) %>%
    c(RunControl$MaximumTimestep) %>%
    unique()
} else c()

timestepsSpatial <- if(OutputOptions$RasterOutputHA){
  seq(RunControl$MinimumTimestep, RunControl$MaximumTimestep, by = OutputOptions$RasterOutputHATimesteps) %>%
    c(RunControl$MaximumTimestep) %>%
    unique()
} else c()

timesteps <- c(timestepsTabular, timestepsSpatial) %>%
  unique() %>%
  sort()

# Species
species <- as.character(HabitatModel$Name)

# Species codes
# speciesCodes <- read_csv(file.path("D:/nestweb", tabularDataDir, "species-codes.csv"), show_col_types = FALSE)

# Invalid Habitat
invalidHabitatLookup <- InvalidHabitat %>%
  left_join(tibble(
    StratumId = c(Stratum$Name, rep(NA, length(Stratum$Name))),
    StratumId2 = rep(Stratum$Name, 2)),
    multiple = "all") %>%
  left_join(tibble(
    StateClassId = c(StateClass$Name, rep(NA, length(StateClass$Name))),
    StateClassId2 = rep(StateClass$Name, 2)),
    multiple = "all") %>%
  left_join(tibble(
    Species = c(names(SpeciesID), rep(NA, length(names(SpeciesID)))),
    Species2 = rep(names(SpeciesID), 2)),
    multiple = "all") %>%
  dplyr::select("Species" = "Species2",
                "StateClassId" = "StateClassId2",
                "StratumId" = "StratumId2") %>%
  mutate(HabitatMask = 0) %>%
  bind_rows(anti_join(expand_grid(Species = names(SpeciesID), StateClassId = StateClass$Name, StratumId = Stratum$Name), .)) %>%
  unique() %>%
  mutate(StateClassId = StateClassId %>% lookup(StateClass$Name, StateClass$StateClassId),
         StratumId = StratumId %>% lookup(Stratum$Name, Stratum$StratumId),
         StateClassStratumId = (StateClassId * 10) + StratumId) %>%
  select(Species, StateClassStratumId, HabitatMask)

# Square meter to hectare conversion
scaleFactor <- 0.0001

## Load data ----
# Load Stratum raster
stratum <- rast(InitialConditionsSpatial$StratumFileName)

# Load mean decay raster
meanDecay <- rast(file.path(scriptDir, "mean-decay.tif"))

# Create a template raster
templateRaster <- stratum
templateRaster[!is.na(templateRaster)] <- 1
names(templateRaster) <- "habitatSuitability"

# Pixel resolution
cellResolution <- templateRaster %>% res()

# Area of a pixel (square meter)
cellArea <- cellResolution[1]^2

# Get Strata and site values
StrataData <- data.frame(
  StratumId = rast(InitialConditionsSpatial$StratumFileName)[] %>%
    as.vector() %>%
    lookup(Stratum$StratumId, Stratum$Name),
  SecondaryStratumId = rast(InitialConditionsSpatial$SecondaryStratumFileName)[] %>%
    as.vector() %>%
    lookup(SecondaryStratum$SecondaryStratumId, SecondaryStratum$Name),
  Site = rast(Site$FileName)[] %>% as.vector() %>% lookup(SiteType$ID, SiteType$Name))

# Get mean values for other habitat model variables
rawNestwebData <- read_csv(file.path(scriptDir, "Habitat selection - full dataset (30.03.2022).csv"), show_col_types = FALSE) %>%
  mutate(Num_2BI = Num_2BI_FD + Num_2BI_Sx + Num_2BI_Pl) %>% 
  dplyr::select(Num_Trees, Num_2BI) %>% 
  summarise(MeanNumTrees = mean(Num_Trees),
            MeanNum2BI = mean(Num_2BI))

## Setup files and folders ----

e <- ssimEnvironment()
transferDir <- e$TransferDirectory

# Predict Habitat ----
progressBar(type = "message", message = "Running main code...")
progressBar(type = "begin", totalSteps = length(iterations) * length(timesteps) * length(species))

# Build parameter sampling table
# zzz: apply up2date() to models - get claude's help
models <- HabitatModel$ModelFileName %>%
  map(~{
    e <- new.env(parent = emptyenv())
    load(.x, envir = e)
    get(ls(e)[1], envir = e)
  }) %>%
  map(up2date) %>%
  set_names(HabitatModel$Name)

# Parameterize the sampling distribution for each parameter and model
parameterTable <- imap_dfr(
  models,
  ~{
    # Get list of model parameters
    parameters <- fixef(.x)[1]$cond %>% names()
    means <- fixef(.x)[1]$cond
    stdErrors <- summary(.x)[6]$coefficients$cond[,2] 
    
    tibble(
      species = .y,
      parameter = parameters,
      mean = means,
      sd = stdErrors)})

for(iteration in iterations){
  # Sample parameters
  parameterTable <- parameterTable %>% 
    mutate(value = map2_dbl(mean, sd, rnorm, n = 1))
  
  # For each species model in models, update the aspen and diameter coefficients with value in parameterTable
  for(i in seq(1:length(models))){
    
    # Get model name to filter paramer tables
    modelName <- names(models)[i]
    parameterTableIt <- parameterTable %>% filter(species == modelName)
    
    # Get species model and overwrite parameter coefficients
    speciesModel <- models[[i]]
    
    speciesModel$fit$par[names(speciesModel$fit$par) == "beta"] <- parameterTableIt %>% pull(value)
    speciesModel$fit$parfull[names(speciesModel$fit$parfull) == "beta"] <- parameterTableIt %>% pull(value)
    
    models[[i]] <- speciesModel
  }
  
  for(timestep in timesteps){
  
    # Load aspen cover and diameter
    aspenCover <- datasheet(myScenario, "stsim_OutputSpatialStockGroup",
                            lookupsAsFactors = FALSE) %>%
      filter(Iteration == iteration, Timestep == timestep,
             StockGroupId == "Aspen Cover (%) [Type]") %>%
      pull(Filename) %>%
      rast()
    
    # Convert aspen raster to proportion (rather than percentage)
    aspenCover[] <- aspenCover[]/100
    names(aspenCover) <- "Perc_At"
    
    diameter <- datasheet(myScenario, "stsim_OutputSpatialStockGroup",
                          lookupsAsFactors = FALSE) %>%
      filter(Iteration == iteration, Timestep == timestep,
             StockGroupId == "Diameter (cm) [Type]") %>%
      pull(Filename) %>%
      rast()
    names(diameter) <- "Median_DBH"
    
    # Convert TST raster to binary Y/N cut data
    cutReclasMatrix <- matrix(data = c(0,60,1), ncol = 3, nrow = 1)
    
    cutRaster <- datasheet(myScenario, "stsim_OutputSpatialTST",
                           lookupsAsFactors = FALSE) %>%
      filter(Iteration == iteration, Timestep == timestep) %>%
      pull(Filename) %>%
      rast()
    
    crs(cutRaster) <- crs(templateRaster)
    
    cut <- cutRaster %>% 
      classify(rcl = cutReclasMatrix, 
                     include.lowest = TRUE,
                     others = NA) %>% 
      buffer(width = cellResolution[1]) %>% 
      mask(mask = templateRaster)
    
    distanceToCut <- distance(cut, target = 0, unit = 'm')
    names(distanceToCut) <- "dist_to_cut"
    
    # Calculate distance to forest edge
    forestRcl <- matrix(c(20,30,40,41,42,43,44,50,60,70,90,
                          NA,NA,1,1,1,1,1,NA,NA,NA,NA),
                        ncol = 2, nrow = 11)
    
    stateClass <- datasheet(myScenario, "stsim_OutputSpatialState",
                            lookupsAsFactors = FALSE) %>%
      filter(Iteration == iteration, Timestep == timestep) %>%
      pull(Filename) %>%
      rast()

    crs(stateClass) <- crs(templateRaster)

    forest <- stateClass %>% classify(rcl = forestRcl)

    distanceToForest <- distance(forest, unit = 'm') %>%
      mask(templateRaster)
    names(distanceToForest) <- "edge_near"

    # Create dataframe of habitat suitability model inputs
    habitatSuitabilityDf <- data.frame(Perc_At = aspenCover[],
                                       Median_DBH = diameter[],
                                       edge_near = distanceToForest[],
                                       Num_Trees = rawNestwebData$MeanNumTrees,
                                       Num_2BI = rawNestwebData$MeanNum2BI,
                                       Mean_decay = meanDecay[],
                                       dist_to_cut = distanceToCut[],
                                       Site = StrataData$Site)

    # Combine state class and stratum raster to create habitat masking raster
    stateClassStratum <- ((stateClass * 10) + stratum) %>% suppressWarnings()
    
    for(aSpecies in species){
      # Create habitat mask
      reclassMatrix <- invalidHabitatLookup %>%
        filter(Species == aSpecies) %>%
        select(-Species) %>%
        as.matrix()

      habitatMask <- stateClassStratum %>%
        classify(rcl = reclassMatrix, others = NA)

      # Predict habitat suitability
      model <- models[[aSpecies]]
      habitatSuitabilityDf$pred <- predict(model, newdata = habitatSuitabilityDf, type = "response", allow.new.levels = TRUE)
      habitatSuitabilityDf$pred[is.nan(habitatSuitabilityDf$pred)] <- NA
      habitatSuitabilityDf$invalidHabitat <- habitatMask[] %>% as.vector()

      # Mask out invalid habitat
      habitatSuitabilityDf <- habitatSuitabilityDf %>%
        mutate(finalHabitat = case_when(!is.na(invalidHabitat) ~ invalidHabitat,
                                        is.na(invalidHabitat) ~ pred))

      # Output habitat raster
      if(timestep %in% timestepsSpatial) {
        hsFilename <- file.path(transferDir, str_c("hs.sp", SpeciesID[aSpecies], ".it", iteration, ".ts", timestep, ".tif"))

        # Create and write habitat raster
        rast(templateRaster, vals = habitatSuitabilityDf$finalHabitat) %>%
          writeRaster(hsFilename,
                      datatype = "FLT4S",
                      overwrite = TRUE,
                      NAflag = -9999)

        OutputSpatialHabitat <- rbind(OutputSpatialHabitat, data.frame(
          Iteration = as.integer(iteration),
          Timestep = as.integer(timestep),
          Species = aSpecies,
          FileName = hsFilename
        ))

        # Output habitat change raster
        if(OutputOptions$RasterOutputHAC){
          if(timestep != min(timestepsSpatial)){
            hscFilename <- file.path(transferDir, str_c("hsc.sp", SpeciesID[aSpecies], ".it", iteration, ".ts", timestep, ".tif"))

            habitatData <- data.frame(
              tsMin = rast(file.path(transferDir, str_c("hs.sp", SpeciesID[aSpecies], ".it", iteration, ".ts", min(timestepsSpatial), ".tif")))[] %>% as.vector(),
              tsCurrent = rast(hsFilename)[] %>% as.vector()) %>%
              mutate(difference = tsCurrent - tsMin)

            rast(templateRaster, vals = habitatData$difference) %>%
              writeRaster(hscFilename,
                          overwrite = TRUE,
                          NAflag = -9999)

            OutputSpatialHabitatChange <- rbind(OutputSpatialHabitatChange, data.frame(
              Iteration = as.integer(iteration),
              Timestep = as.integer(timestep),
              Species = aSpecies,
              FileName = hscFilename
            ))
          }
        }
      }
      
      # Calculate tabular output and append to 
      if(timestep %in% timestepsTabular) {
        OutputHabitatAmount <- bind_rows(
          OutputHabitatAmount, 
          habitatSuitabilityDf %>% 
            dplyr::select(Amount = finalHabitat) %>% 
            bind_cols(StrataData) %>% 
            filter(!is.na(Amount)) %>% 
            group_by(StratumId, SecondaryStratumId, Site) %>%
            summarise(Amount = sum(Amount) * cellArea * scaleFactor, .groups = "drop") %>% 
            mutate(
              Timestep = as.integer(timestep),
              Iteration = as.integer(iteration),
              Species = aSpecies
            ))}
    
      # Increment progress bar
      progressBar()
    }
  }
}

# Save outputs
saveDatasheet(myScenario, OutputSpatialHabitat, "stsimNestweb_OutputSpatialHabitat")
saveDatasheet(myScenario, OutputSpatialHabitatChange, "stsimNestweb_OutputSpatialHabitatChange")
saveDatasheet(myScenario, OutputHabitatAmount, "stsimNestweb_OutputHabitatAmount")


# Clean up ----

# Wrap up SyncroSim progress bar
progressBar("end")