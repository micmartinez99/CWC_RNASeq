# The purpose of this script is to take the raw hydrolyzed urolithin data
# obtained from Dr. Anthony Provatas and Nuoxi Fan and extract the relevant information.
# This data contains the levels of 9 urolithin metabolites (IsoA, UroA, UroB, UroC, UroD, UroE,
# UroM5, UroM6, and UroM7) at three timepoints (V1: Visit 1, V2: Before walnuts, V3: After walnuts)
# We are interested in the delta between visit 3 and visit 2 to assess the ability to produce
# urolithins. 

# Load libraries
library(tidyverse)
library(dplyr)

# Initialize a function to create output directories
newDir <- function(x) {
  if (!dir.exists(x)) {
    dir.create(x)
  }
}

# Initialize an output directory
newDir("Outputs/001_Tidy_Urolithins")

################################################################################

# Read in the raw urolithin Data
uros <- read.csv("Data/Urolithin_Data/Hydrolyzed_Urolithin_Data_Batches1_2.csv")
colnames(uros) <- uros[1,]
colnames(uros) <- trimws(colnames(uros))
uros <- uros[-1,]

# Set a new column in the dataframe that will hold information about what visit this is
# as well as patient ID
uros <- uros %>%
  mutate(visit = case_when(
    grepl("V1$", `Sample ID`) ~ "V1",
    grepl("V2$", `Sample ID`) ~ "V2",
    grepl("V3$", `Sample ID`) ~ "V3",
    TRUE ~ NA_character_)) %>%
  mutate(patient_id = paste0("p", sub("V.*", "", `Sample ID`)))
  

# Remove any rows that are visit 1
uros_23 <- uros[uros$visit != "V1",]

# Pivot longer the dataframe for easier manipulation
uros_23 <- uros_23 %>%
  pivot_longer(-c(`Sample ID`, visit, patient_id), names_to = "Metabolite", values_to = "quant")

# Get a vector of metabolites to iterate over
metabolites <- unique(uros_23$Metabolite)

# Initialize a list to store the delta values
delta_list <- vector("list", length(metabolites))
names(delta_list) <- metabolites

# Calculate delta values (V3-V2) for each metabolite
for (i in metabolites) {
  
  # Subset the dataframe based on metabolite i
  subset <- uros_23[uros_23$Metabolite == i,]
  
  # Pivot wider the data
  subset <- subset %>%
    pivot_wider(names_from = visit, values_from = quant, names_prefix = "timepoint_")
  
  # Extract the visit 2 information
  visit2 <- subset[,c(2,3,4)]
  visit2 <- na.omit(visit2)
  
  # Extract the visit 3 information
  visit3 <- subset[,c(2,3,5)]
  visit3 <- na.omit(visit3)
  
  # Append visit 3 values to the visit 2 dataframe 
  visit2$timepoint_V3 <- visit3$timepoint_V3
  
  # Re name the dataframe
  delta <- visit2
  delta$timepoint_V2 <- as.numeric(delta$timepoint_V2)
  delta$timepoint_V3 <- as.numeric(delta$timepoint_V3)
  
  # Convert quant columns to numeric and calculate the delta
  deltaVal <- delta %>%
    mutate(delta = timepoint_V3 - timepoint_V2)
  
  delta_list[[i]] <- deltaVal$delta
  
  # Save delta values as a csv file
  metaboliteName <- gsub(" ", "_", i)
  fileName <- paste(metaboliteName, "Delta_Values.csv", sep = "_")
  write.csv(deltaVal, paste("Outputs/001_Tidy_Urolithins", fileName, sep = "/"))
}

# Convert the list to a dataframe, set column names to metabolites
delta_values_all <- do.call(cbind, lapply(delta_list, as.data.frame))
colnames(delta_values_all) <- names(delta_list)
rownames(delta_values_all) <- uro_delta$patient_id

# Save as a csv
write.csv(delta_values_all, file = "Outputs/001_Tidy_Urolithins/All_Delta_Values.csv")

