biomass_data <- read.csv("data/biomass_conversions.csv")
biomass_data_frozen <- biomass_data %>%
  filter(Preservation.Method == "Frozen")
