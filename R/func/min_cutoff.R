# Inputs
# data <- biomass_data_taxaclean
# weight_measure_1 <- expr(wm_g)
# weight_measure_1 <- expr(dm_g)

min_cutoff <- function(data, 
                       weight_measure_1, 
                       weight_measure_2, 
                       cut_off = 0.001){  #What goes into a function
  # Where the function starts
  
  # Filtering at cutoff
  output <- data |> 
    filter({{weight_measure_1}} >= cut_off,
           {{weight_measure_2}} >= cut_off)
  
  # Give me output
  return(output)
}

# Test function
min_cutoff(data = biomass_data_taxaclean, dm_g, wm_g, 0.00001) 

# Check IRL
biomass_data_taxaclean |>
  filter(wm_g >= 0.001, dm_g >=0.001) |> nrow()

##########
.data = biomass_data_taxaclean
ref_weight_measure <- expr(wm_g)
other_weight_measure <- expr(dm_g)

# https://dplyr.tidyverse.org/articles/programming.html#data-masking
# https://brad-cannell.github.io/r_notes/tidy-evaluation.html
check_mass_diff <- function(data, 
                            ref_weight_measure, 
                            other_weight_measure){
  
  #reference <- data[[ref_weight_measure]]
  #other <- .data[[rlang::as_label(other_weight_measure)]]
  reference <- data |> 
    pull({{ref_weight_measure}})
  
  
  # if(reference >= other)
  #   rlang::warn("Reference weight measure is heavier than other")
  
  reference
  
}

check_mass_diff(biomass_data_taxaclean, wm_g)



biomass_data_taxaclean |> 
  filter(!wm_g >= dm_g)

