# Set environment variable TZ when running on AWS EC2 instance
Sys.setenv(TZ='UTC')

# Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

# Load libraries
library(rsyncrosim)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(terra))
suppressPackageStartupMessages(library(MuMIn))
suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(unmarked))

# Setup ----
progressBar(type = "message", message = "Preparing inputs...")

## Function definitions ----
# Define function to facilitate recoding a vector using a look-up table
lookup <- function(x, old, new){
  dplyr::recode(x, !!!set_names(new, old))
}

## Connect to SyncroSim ----
myScenario <- scenario()

# Load relevant datasheets
# stsim
RunControl <- datasheet(myScenario, "stsim_RunControl")
Stratum <- datasheet(myScenario, "stsim_Stratum")
SecondaryStratum <- datasheet(myScenario, "stsim_SecondaryStratum")
StateClass <- datasheet(myScenario, "stsim_StateClass")
InitialConditionsSpatial <- datasheet(myScenario, "stsim_InitialConditionsSpatial")
OutputSpatialState <- datasheet(myScenario, "stsim_OutputSpatialState")

#stsimsf
OutputSpatialStockGroup <- datasheet(myScenario, "stsimsf_OutputSpatialStockGroup")
StockType <- datasheet(myScenario, "stsimsf_StockType")
StockGroup <- datasheet(myScenario, "stsimsf_StockGroup")
StockTypeGroupMembership <- datasheet(myScenario, "stsimsf_StockTypeGroupMembership")

# stsimNetweb
SiteType <- datasheet(myScenario, "stsimNestweb_SiteType")
SpeciesID <- datasheet(myScenario, "stsimNestweb_Species", includeKey = TRUE) %>% 
  pull(SpeciesID, name = Name)
Site <- datasheet(myScenario, "stsimNestweb_SiteValue")
OutputOptions <- datasheet(myScenario, "stsimNestweb_OutputOptions")
HabitatModel <- datasheet(myScenario, "stsimNestweb_HabitatModel")
InvalidHabitat <- datasheet(myScenario, "stsimNestweb_InvalidHabitat")
OutputHabitatAmount <- data.frame()

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
species <- HabitatModel$Name

# Invalid Habitat
invalidHabitatLookup <- InvalidHabitat %>% 
  left_join(tibble(
    StratumID = c(Stratum$Name, rep(NA, length(Stratum$Name))),
    StratumID2 = rep(Stratum$Name, 2)), 
    multiple = "all") %>% 
  left_join(tibble(
    StateClassID = c(StateClass$Name, rep(NA, length(StateClass$Name))),
    StateClassID2 = rep(StateClass$Name, 2)), 
    multiple = "all") %>% 
  left_join(tibble(
    Species = c(names(SpeciesID), rep(NA, length(names(SpeciesID)))),
    Species2 = rep(names(SpeciesID), 2)), 
    multiple = "all") %>% 
  dplyr::select("Species" = "Species2",
                "StateClassID" = "StateClassID2",
                "StratumID" = "StratumID2") %>% 
  mutate(HabitatMask = 0) %>% 
  bind_rows(anti_join(expand_grid(Species = names(SpeciesID), StateClassID = StateClass$Name, StratumID = Stratum$Name), .)) %>% 
  unique() %>% 
  mutate(StateClassID = StateClassID %>% lookup(StateClass$Name, StateClass$ID),
         StratumID = StratumID %>% lookup(Stratum$Name, Stratum$ID),
         StateClassStratumID = (StateClassID * 10) + StratumID) %>% 
  select(Species, StateClassStratumID, HabitatMask)

# Square meter to hectare conversion
scaleFactor <- 0.0001

## Load spatial data ----
# Load Stratum raster
stratum <- rast(InitialConditionsSpatial$StratumFileName)

# Create a template raster
templateRaster <- stratum
templateRaster[!is.na(templateRaster)] <- 1
names(templateRaster) <- "habitatSuitability"

# Pixel resolution
cellResolution <- templateRaster %>% res()

# Pixel count
#cellCount <- freq(templateRaster, value = 1)$count

# Area of a pixel (square meter)
cellArea <- cellResolution[1]^2

# Total area
# totalArea <- cellCount * cellArea

# Get Strata and site values
StrataData <- data.frame(
  StratumID = rast(InitialConditionsSpatial$StratumFileName)[] %>% as.vector(),
  SecondaryStratumID = rast(InitialConditionsSpatial$SecondaryStratumFileName)[] %>% as.vector(),
  Site = rast(Site$FileName)[] %>% as.vector() %>% lookup(SiteType$ID, SiteType$Name))

# Get area of unique strata and site combinations
# Area <- StrataData %>% 
#   count(StratumID, SecondaryStratumID, Site) %>% 
#   rename(Area = n) %>% 
#   mutate(Area = Area * cellArea * scaleFactor)
# 
# # Add Area to StrataData
# StrataData<- StrataData %>% 
#   left_join(Area)

## Setup files and folders ----

# Create temp folder, ensure it is empty
tempDir <- ssimEnvironment()$TempDirectory

spatialOutputDir <- file.path(tempDir, "SpatialOutputs") %>% normalizePath(mustWork = FALSE)

unlink(spatialOutputDir, recursive = TRUE, force = TRUE)

dir.create(spatialOutputDir)

# Generate filenames for potential outputs

## Handle empty values ----

# Main Code Here ----
progressBar(type = "message", message = "Running main code...")
progressBar(type = "begin", totalSteps = length(iterations) * length(timesteps) * length(species))

# Build parameter sampling table
modelNames <- map_chr(HabitatModel$ModelFileName, load)
for(m in HabitatModel$ModelFileName) load(m)
models <- modelNames %>% 
  map(get) %>% 
  set_names(HabitatModel$Name)
rm(list = modelNames)

# Parameterize the sampling distribution for each parameter and model
parameterTable <- imap_dfr(
  models,
  ~{
    # Connect clean parameter names to variable names in model fit
    parameters <- c('Aspen Cover' = 'scale(Perc_At)',
                    'Diameter' = 'scale(Median_DBH)')
    means <- coef(.x)[parameters]
    stdErrors <- .x$coefArray[, 'Std. Error', parameters] %>% 
      apply(MARGIN = 2, FUN = mean, na.rm = TRUE) # MARGIN = 2 is used to average within columns 
    
    tibble(
      species = .y,
      parameter = names(parameters),
      mean = means,
      sd = stdErrors)}) # Double check this assumption

for(iteration in iterations){
  # Sample parameters
  parameterTable <- parameterTable %>% 
    mutate(value = map2_dbl(mean, sd, rnorm, n = 1))
  
  for(timestep in timesteps){
    
    # Load aspen cover and diameter
    aspenCover <- datasheetSpatRaster(
      ssimObject = myScenario, 
      datasheet = "stsimsf_OutputSpatialStockGroup", 
      iteration = iteration, 
      timestep = timestep,
      filterColumn = "StockGroupID",
      filterValue = "Deciduous Canopy Cover (%) [Type]")
    
    # Convert aspen raster to proportion (rather than percentage)
    aspenCover[] <- aspenCover[]/100
    names(aspenCover) <- 'Perc_At'
    
    # diameter <- datasheetSpatRaster(
    #   ssimObject = myScenario, 
    #   datasheet = "stsimsf_OutputSpatialStockGroup", 
    #   iteration = iteration, 
    #   timestep = timestep,
    #   filterColumn = "StockGroupID",
    #   filterValue = "Diameter (cm) [Type]")
    
    # Convert TST raster to binary Y/N cut data
    # cut <- datasheetSpatRaster(
    #   ssimObject = myScenario, 
    #   datasheet = "stsim_OutputSpatialTST", 
    #   iteration = iteration, 
    #   timestep = timestep)[] %>% 
    #   as.vector() %>% 
    #   {case_when(. <= 60 ~ "Y", TRUE ~ "N")}
    
    # Create dataframe of habitat suitability model inputs
    habitatSuitabilityDf <- data.frame(Perc_At = aspenCover[], 
                                     # Median_DBH = diameter[],
                                       Median_DBH = 0,
                                       edge_near = 0,
                                       Num_2BI = 0,
                                       Mean_decay = 0,
                                       dist_to_cut = 0,
                                     # cut_harvest0 = cut, 
                                       cut_harvest0 = 'N', 
                                       Site = StrataData$Site)
    
    # Load state class raster
    stateClass <- datasheetSpatRaster(
      ssimObject = myScenario, 
      datasheet = "stsim_OutputSpatialState", 
      iteration = iteration, 
      timestep = timestep) 
    
    # Combine state class and stratum raster to create habitat masking raster
    stateClassStratum <- ((stateClass * 10) + stratum) %>% suppressWarnings()
    
    for(aSpecies in species){
      # Predict habitat suitability
      model <- models[[aSpecies]]
      habitatSuitabilityDf$pred <- predict(model, newdata = habitatSuitabilityDf, type = "response", allow.new.levels = TRUE)
      habitatSuitabilityDf$pred[is.nan(habitatSuitabilityDf$pred)] <- NA
      
      # Output habitat raster
      if(timestep %in% timestepsSpatial) {
        outputFilename <- file.path(spatialOutputDir, str_c("hs.sp", SpeciesID[aSpecies], ".it", iteration, ".ts", timestep, ".tif")) %>% 
          normalizePath(mustWork = FALSE)
        
        # Create habitat mask
        reclassMatrix <- invalidHabitatLookup %>% 
          filter(Species == aSpecies) %>%  
          select(-Species) %>% 
          as.matrix()
        
        habitatMask <- stateClassStratum %>% 
          classify(rcl = reclassMatrix)
        
       maskedHabitat <- data.frame(
         invalidHabitat = habitatMask[] %>% as.vector(),
         habitat = habitatSuitabilityDf$pred) %>% 
         mutate(finalHabitat = case_when(!is.na(invalidHabitat) ~ invalidHabitat,
                                         is.na(invalidHabitat) ~ habitat))
        
        # Create and write habitat raster
        rast(templateRaster, vals = maskedHabitat$finalHabitat) %>% 
          writeRaster(outputFilename,
                      overwrite = TRUE,
                      NAflag = -9999)
        
        # Output habitat change raster
        if(OutputOptions$RasterOutputHAC){
          if(timestep != min(timestepsSpatial)){
            outputFilename = file.path(spatialOutputDir, str_c("hsc.sp", SpeciesID[aSpecies], ".it", iteration, ".ts", timestep, ".tif")) %>% 
              normalizePath(mustWork = FALSE)
            
            habitatData <- data.frame(
              ts2016 = rast(file.path(spatialOutputDir, str_c("hs.sp", SpeciesID[aSpecies], ".it", iteration, ".ts", min(timestepsSpatial), ".tif")))[] %>% as.vector(),
              tsCurrent = rast(file.path(spatialOutputDir, str_c("hs.sp", SpeciesID[aSpecies], ".it", iteration, ".ts", timestep, ".tif")))[] %>% as.vector()) %>% 
              mutate(difference = tsCurrent - ts2016)
            
            rast(templateRaster, vals = habitatData$difference) %>% 
              writeRaster(outputFilename, 
                          overwrite = TRUE,
                          NAflag = -9999)
          }
        }
      }
      
      # Calculate tabular output and append to 
      if(timestep %in% timestepsTabular) {
        OutputHabitatAmount <- bind_rows(
          OutputHabitatAmount, 
          habitatSuitabilityDf %>% 
            dplyr::select(Amount = pred) %>% 
            bind_cols(StrataData) %>% 
            filter(!is.na(Amount)) %>% 
            group_by(StratumID, SecondaryStratumID, Site) %>% 
            summarise(Amount = sum(Amount) * cellArea * scaleFactor, .groups = "drop") %>% 
            mutate(
              Timestep = timestep,
              Iteration = iteration,
              Species = aSpecies
            ))}
      
      # Increment progress bar
      progressBar()
    }
  }
}

# Save spatial outputs
OutputSpatialHabitat <- tibble(FileName = list.files(spatialOutputDir, "hs\\..+tif", full.names = TRUE) %>% normalizePath()) %>%
  mutate(
    temp = basename(FileName),
    Iteration = temp %>% str_extract("it\\d+") %>% str_replace("it", "") %>% as.numeric(),
    Timestep = temp %>% str_extract("ts\\d+") %>% str_replace("ts", "") %>% as.numeric(),
    Species = temp %>% str_extract("sp\\d+") %>% str_replace("sp", "") %>% as.numeric(),
    Species = lookup(Species, SpeciesID, names(SpeciesID))) %>% 
  dplyr::select(-temp) %>% 
  as.data.frame()

saveDatasheet(myScenario, OutputSpatialHabitat, "stsimNestweb_OutputSpatialHabitat")

OutputSpatialHabitatChange <- tibble(FileName = list.files(spatialOutputDir, "hsc\\..+tif", full.names = TRUE) %>% normalizePath()) %>%
  mutate(
    temp = basename(FileName),
    Iteration = temp %>% str_extract("it\\d+") %>% str_replace("it", "") %>% as.numeric(),
    Timestep = temp %>% str_extract("ts\\d+") %>% str_replace("ts", "") %>% as.numeric(),
    Species = temp %>% str_extract("sp\\d+") %>% str_replace("sp", "") %>% as.numeric(),
    Species = lookup(Species, SpeciesID, names(SpeciesID))) %>% 
  dplyr::select(-temp) %>% 
  as.data.frame()

saveDatasheet(myScenario, OutputSpatialHabitatChange, "stsimNestweb_OutputSpatialHabitatChange")

# Save tabular output
OutputHabitatAmount <- OutputHabitatAmount %>% 
  mutate(
    StratumID = StratumID %>% lookup(Stratum$ID, Stratum$Name),
    SecondaryStratumID = SecondaryStratumID %>% lookup(SecondaryStratum$ID, SecondaryStratum$Name))
saveDatasheet(myScenario, OutputHabitatAmount, "stsimNestweb_OutputHabitatAmount")


# Clean up ----

# Wrap up SyncroSim progress bar
progressBar("end")