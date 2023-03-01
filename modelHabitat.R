# zzz:
# Load spatial data to disambiguating habitat suitability (strata, site)
# Create a datasheet for Site
# Build habitat suitability maps

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


## Connect to SyncroSim ----
myScenario <- scenario()

# Load relevant datasheets
# stsim
RunControl <- datasheet(myScenario, "stsim_RunControl")
Stratum <- datasheet(myScenario, "stsim_Stratum")
SecondaryStrtaum <- datasheet(myScenario, "stsim_SecondaryStratum")
StateClass <- datasheet(myScenario, "stsim_StateClass")
InitialConditionsSpatial <- datasheet(myScenario, "stsim_InitialConditionsSpatial")
OutputSpatialState <- datasheet(myScenario, "stsim_OutputSpatialState")

#stsimsf
OutputSpatialStockGroup <- datasheet(myScenario, "stsimsf_OutputSpatialStockGroup")
StockType <- datasheet(myScenario, "stsimsf_StockType")
StockGroup <- datasheet(myScenario, "stsimsf_StockGroup")
StockTypeGroupMembership <- datasheet(myScenario, "stsimsf_StockTypeGroupMembership")

# stsimNetweb
OutputOptions <- datasheet(myScenario, "stsimNestweb_OutputOptions")
OutputOptionsSpatial <- datasheet(myScenario, "stsimNestweb_OutputOptionsSpatial")
HabitatModel <- datasheet(myScenario, "stsimNestweb_HabitatModel")
SpeciesID <- datasheet(myScenario, "stsimNestweb_Species", includeKey = TRUE) %>% 
  pull(SpeciesID, name = Name)

## Setup Parameters ----
# Iterations 
iterations <- seq(RunControl$MinimumIteration, RunControl$MaximumIteration)

# Timesteps
timestepsTabular <- if(OutputOptions$SummaryOutputHA){
  seq(RunControl$MinimumTimestep, RunControl$MaximumTimestep, by = OutputOptions$SummaryOutputHATimesteps) %>% 
    c(RunControl$MaximumTimestep) %>% 
    unique()
} else c()

timestepsSpatial <- if(OutputOptionsSpatial$RasterOutputHA){
  seq(RunControl$MinimumTimestep, RunControl$MaximumTimestep, by = OutputOptionsSpatial$RasterOutputHATimesteps) %>% 
    c(RunControl$MaximumTimestep) %>% 
    unique()
} else c()

timesteps <- c(timestepsTabular, timestepsSpatial) %>% 
  unique() %>% 
  sort()

# Species
species <- HabitatModel$Name

# Square meter to hectare conversion
scaleFactor <- 0.0001

# Load template raster
templateRaster <- rast(InitialConditionsSpatial$StratumFileName)
templateRaster[!is.na(templateRaster)] <- 1
names(templateRaster) <- "habitatSuitability"

# Pixel resolution
cellResolution <- templateRaster %>% res()

# Pixel count
cellCount <- freq(templateRaster, value = 1)$count

# Area of a pixel (square meter)
cellArea <- cellResolution[1]^2

# Total area
totalArea <- cellCount * cellArea



## Setup files and folders ----

# Create temp folder, ensure it is empty
tempDir <- ssimEnvironment()$TempDirectory

# Generate filenames for potential outputs

## Handle empty values ----

## Function definitions ----

# Main Code Here zzz ----
progressBar(type = "message", message = "Running main code...")

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
    aspenCover <- datasheetRaster(
      ssimObject = myScenario, 
      datasheet = "stsimsf_OutputSpatialStockGroup", 
      iteration = iteration, 
      timestep = timestep,
      filterColumn = "StockGroupID",
      filterValue = "Aspen Cover (%) [Type]")
    
    # Convert aspen raster to proportion (rather than percentage)
    aspenCover[] <- aspenCover[]/100
    
    diameter <- datasheetRaster(
      ssimObject = myScenario, 
      datasheet = "stsimsf_OutputSpatialStockGroup", 
      iteration = iteration, 
      timestep = timestep,
      filterColumn = "StockGroupID",
      filterValue = "Diameter (cm) [Type]")
    
    # Create dataframe of habitat suitability model inputs
    habitatSuitabilityDf <- data.frame(Perc_At = aspenCover[], 
                                       Median_DBH = diameter[],
                                       edge_near = 0,
                                       Num_2BI = 0,
                                       Mean_decay = 0,
                                       dist_to_cut = 0,
                                       # zzz: how to handle categorical variables??
                                       cut_harvest0 = "N",
                                       Site = "YY")
    
    for(aSpecies in species){
      # Predict habitat suitability
      model <- models[[aSpecies]]
      habitatSuitabilityDf$pred <- predict(model, newdata = habitatSuitabilityDf, type = "response", allow.new.levels = TRUE)
      
      # Output raster
      outputFilename <- file.path(tempDir, str_c("hs.sp_", SpeciesID[aSpecies], ".it", iteration, ".ts", timestep, ".tif")) %>% 
        normalizePath(mustWork = FALSE)
      
      rast(templateRaster, vals = habitatSuitabilityDf$pred) %>% 
        writeRaster(outputFilename, overwrite = TRUE)
      
    }
  }
}

OutputSpatialHabitat <- tibble(FileName = list.files(tempDir, ".tif", full.names = TRUE)) #%>% 

# Clean up ----

# Wrap up SyncroSim progress bar
progressBar("end")