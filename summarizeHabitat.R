# Set environment variable TZ when running on AWS EC2 instance
Sys.setenv(TZ='UTC')

# Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

# Load libraries
library(rsyncrosim)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(terra))

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
OutputOptions <- datasheet(myScenario, "stsimNestweb_OutputOptions")
OutputSpatialHabitat <- datasheet(myScenario, "stsimNestweb_OutputSpatialHabitat")
SpeciesID <- datasheet(myScenario, "stsimNestweb_Species", includeKey = TRUE) %>% 
  pull(SpeciesID, name = Name)

## Setup Parameters ----
# Timesteps 
timesteps <- OutputSpatialHabitat$Timestep %>% 
  unique() %>% 
  sort()

# Species
species <- OutputSpatialHabitat$Species %>% 
  unique() %>% 
  as.character() %>% 
  sort()

## Setup files and folders ----

# Create temp folder, ensure it is empty
tempDir <- ssimEnvironment()$TempDirectory

spatialOutputDir <- file.path(tempDir, "SpatialOutputs") %>% normalizePath(mustWork = FALSE)

unlink(spatialOutputDir, recursive = TRUE, force = TRUE)

dir.create(spatialOutputDir)

# Main Code Here ----
if(OutputOptions$RasterOutputHAAverage) {
  progressBar(type = "message", message = "Running main code...")
  progressBar(type = "begin", totalSteps = length(timesteps) * length(species))
  
  for(timestep in timesteps){
    # Get all habitat suitability maps for a given timestep
    habitatSuitability <- datasheetSpatRaster(
      ssimObject = myScenario, 
      datasheet = "stsimNestweb_OutputSpatialHabitat", 
      timestep = timestep) 
    
    # Repeat for habitat suitability change maps
    # NB: The first timestep is excluded because no change raster is calculated
    if(OutputOptions$RasterOutputHACAverage) {
      if(timestep != min(timesteps)){
        # Get all habitat suitability change maps for a given timestep
        habitatSuitabilityChange <- datasheetSpatRaster(
          ssimObject = myScenario, 
          datasheet = "stsimNestweb_OutputSpatialHabitatChange",
          timestep = timestep) 
      }
    }
    
    for(aSpecies in species){
      # Determine output filename based on species and tiemstep
      outputHabitatFilename <- file.path(spatialOutputDir, str_c("hsa.sp", SpeciesID[aSpecies], ".ts", timestep, ".tif")) %>% 
        normalizePath(mustWork = FALSE)
      
      # Subset layers by species
      habitatLayerNames <- names(habitatSuitability) %>% 
        str_subset(str_c("sp", SpeciesID[aSpecies], "\\."))
  
      # Calculate spatial averages
      habitatSuitability[[habitatLayerNames]] %>% 
        mean() %>% 
        writeRaster(outputHabitatFilename, 
                    overwrite = TRUE,
                    NAflag = -9999)
      
      # Repeat for habitat suitability change
      if(OutputOptions$RasterOutputHACAverage) {
        if(timestep != min(timesteps)){
          # Determine output filename based on species and tiemstep
          outputHabitatChangeFilename <- file.path(spatialOutputDir, str_c("hsca.sp", SpeciesID[aSpecies], ".ts", timestep, ".tif")) %>% 
            normalizePath(mustWork = FALSE)
          
          # Subset layers by species
          habitatChangeLayerNames <- names(habitatSuitabilityChange) %>% 
            str_subset(str_c("sp", SpeciesID[aSpecies], "\\."))
          
          # Calculate spatial averages
          habitatSuitabilityChange[[habitatChangeLayerNames]] %>% 
            mean() %>% 
            writeRaster(outputHabitatChangeFilename, 
                        overwrite = TRUE,
                        NAflag = -9999)
        }
      }
      
      # Increment 
      progressBar()
    }
  }
  OutputSpatialHabitatAverage <- data.frame(
    FileName = list.files(spatialOutputDir, pattern = "hsa\\..+tif", full.names = TRUE) %>% 
      normalizePath(),
    Iteration = 1) %>% 
    mutate(
      temp = basename(FileName),
      Timestep = temp %>% str_extract("ts\\d+") %>% str_replace("ts", "") %>% as.numeric(),
      Species = temp %>% str_extract("sp\\d+") %>% str_replace("sp", "") %>% as.numeric(),
      Species = lookup(Species, SpeciesID, names(SpeciesID))) %>% 
    dplyr::select(-temp) %>% 
    as.data.frame()
  
  saveDatasheet(myScenario, OutputSpatialHabitatAverage, "stsimNestweb_OutputSpatialHabitatAverage")
  
  OutputSpatialHabitatChangeAverage <- data.frame(
    FileName = list.files(spatialOutputDir, pattern = "hsca\\..+tif", full.names = TRUE) %>% 
      normalizePath(),
    Iteration = 1) %>% 
    mutate(
      temp = basename(FileName),
      Timestep = temp %>% str_extract("ts\\d+") %>% str_replace("ts", "") %>% as.numeric(),
      Species = temp %>% str_extract("sp\\d+") %>% str_replace("sp", "") %>% as.numeric(),
      Species = lookup(Species, SpeciesID, names(SpeciesID))) %>% 
    dplyr::select(-temp) %>% 
    as.data.frame()
  
  saveDatasheet(myScenario, OutputSpatialHabitatChangeAverage, "stsimNestweb_OutputSpatialHabitatChangeAverage")
}
