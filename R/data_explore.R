# This is Josh's R Script

# Load packages
library(tidyverse)

# install.packages("skimr")
# install.packages(c("tidyverse", "skimr"))
library(skimr)

# Reading in data
biomass_data <- read.csv("data/biomass_conversions.csv")

# Filtering Frozen samples
biomass_data_frozen <- biomass_data %>%
  filter(Preservation.Method == "Frozen")

# Check out structure of biomass_data_frozen
str(biomass_data_frozen)
skim(biomass_data_frozen)

