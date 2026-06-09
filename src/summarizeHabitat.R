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
  pull(SpeciesId, name = Name)

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

e <- ssimEnvironment()
transferDir <- e$TransferDirectory

# Main Code Here ----
if(OutputOptions$RasterOutputHAAverage) {
  progressBar(type = "message", message = "Running main code...")
  progressBar(type = "begin", totalSteps = length(timesteps) * length(species))

  OutputSpatialHabitatAverage <- data.frame(
    Iteration = integer(0),
    Timestep = integer(0),
    Species = character(0),
    FileName = character(0)
  )

  OutputSpatialHabitatChangeAverage <- data.frame(
    Iteration = integer(0),
    Timestep = integer(0),
    Species = character(0),
    FileName = character(0)
  )

  for(timestep in timesteps){
    # Get all habitat suitability maps for a given timestep
    hsFiles <- datasheet(myScenario, "stsimNestweb_OutputSpatialHabitat",
                         lookupsAsFactors = FALSE) %>%
      filter(Timestep == timestep) %>%
      pull(FileName)
    habitatSuitability <- rast(hsFiles)
    names(habitatSuitability) <- tools::file_path_sans_ext(basename(hsFiles))

    # Repeat for habitat suitability change maps
    # NB: The first timestep is excluded because no change raster is calculated
    if(OutputOptions$RasterOutputHACAverage) {
      if(timestep != min(timesteps)){
        # Get all habitat suitability change maps for a given timestep
        hscFiles <- datasheet(myScenario, "stsimNestweb_OutputSpatialHabitatChange",
                              lookupsAsFactors = FALSE) %>%
          filter(Timestep == timestep) %>%
          pull(FileName)
        habitatSuitabilityChange <- rast(hscFiles)
        names(habitatSuitabilityChange) <- tools::file_path_sans_ext(basename(hscFiles))
      }
    }

    for(aSpecies in species){
      habitatLayerNames <- names(habitatSuitability) %>%
        str_subset(str_c("sp", SpeciesID[aSpecies], "\\."))

      # Determine output filename based on species and timestep
      outputHabitatFilename <- file.path(transferDir, str_c("hsa.sp", SpeciesID[aSpecies], ".ts", timestep, ".tif"))

      # Calculate spatial averages
      habitatSuitability[[habitatLayerNames]] %>%
        mean() %>%
        writeRaster(outputHabitatFilename,
                    overwrite = TRUE,
                    NAflag = -9999)

      OutputSpatialHabitatAverage <- rbind(OutputSpatialHabitatAverage, data.frame(
        Iteration = 1L,
        Timestep = as.integer(timestep),
        Species = aSpecies,
        FileName = outputHabitatFilename
      ))

      # Repeat for habitat suitability change
      if(OutputOptions$RasterOutputHACAverage) {
        if(timestep != min(timesteps)){
          # Determine output filename based on species and timestep
          outputHabitatChangeFilename <- file.path(transferDir, str_c("hsca.sp", SpeciesID[aSpecies], ".ts", timestep, ".tif"))

          # Subset layers by species
          habitatChangeLayerNames <- names(habitatSuitabilityChange) %>%
            str_subset(str_c("sp", SpeciesID[aSpecies], "\\."))

          # Calculate spatial averages
          habitatSuitabilityChange[[habitatChangeLayerNames]] %>%
            mean() %>%
            writeRaster(outputHabitatChangeFilename,
                        overwrite = TRUE,
                        NAflag = -9999)

          OutputSpatialHabitatChangeAverage <- rbind(OutputSpatialHabitatChangeAverage, data.frame(
            Iteration = 1L,
            Timestep = as.integer(timestep),
            Species = aSpecies,
            FileName = outputHabitatChangeFilename
          ))
        }
      }

      # Increment
      progressBar()
    }
  }

  saveDatasheet(myScenario, OutputSpatialHabitatAverage, "stsimNestweb_OutputSpatialHabitatAverage")
  saveDatasheet(myScenario, OutputSpatialHabitatChangeAverage, "stsimNestweb_OutputSpatialHabitatChangeAverage")
}
