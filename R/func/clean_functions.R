library(tidyverse)
library(rlang)
library(janitor)

#' Excludes values over a certain cut off
#'
#' @param data 
#' @param weight_measure_1 
#' @param weight_measure_2 
#' @param cut_off 
#'
#' @return a data where weight_measure_1 weight_measure_2 are above cut_off value
#'
#' @examples
#' Test function
#' min_cutoff(data = biomass_data_taxaclean, afdm_g, wm_g) 

min_cutoff <- function(data, 
                       weight_measure_1, 
                       weight_measure_2, 
                       cut_off = 0.001){  #What goes into a function
  
  # Filtering at cutoff
  output <- data |> 
    dplyr::filter({{weight_measure_1}} >= cut_off,
           {{weight_measure_2}} >= cut_off)
  
  # Give me output
  return(output)
}

# Inputs
# data <- biomass_data_taxaclean
# weight_measure_1 <- expr(wm_g)
# weight_measure_1 <- expr(dm_g)

# Check IRL
# biomass_data_taxaclean |>
#   filter(wm_g >= 0.001, dm_g >=0.001) |> nrow()

#' Excluding light-heavier ones
#'
#' @param data
#' @param heavier_weight_measure 
#' @param lighter_weight_measure 
#'
#' @return
#'
#' @examples
#' check_mass_diff(biomass_data_taxaclean, wm_g, dm_g)
#'
#'biomass_data_taxaclean |> 
#'  drop_na(dm_g, afdm_g) |> 
#'  check_mass_diff(dm_g, afdm_g)


check_mass_diff <- function(data,  # Must not contain NAs for the two weight measures
                            heavier_weight_measure, 
                            lighter_weight_measure){
  # Pull out the weights independently
  heavier <- data |> 
    dplyr::pull({{heavier_weight_measure}})
  
  lighter <- data |> 
    dplyr::pull({{lighter_weight_measure}})
  
  # Identify where some heavier_weight_measure obs are lighter than lighter_weight_measure (indicate possible measure error)
  if(any(heavier < lighter)){
    num_obs_lighter <- length(which(heavier < lighter))
    rlang::warn(paste0("Excluding ", num_obs_lighter, " observations where heavier_weight_measure is less than lighter_weight_measure"))
    
  # Exclude the lighter ones
    data |> 
      dplyr::filter({{heavier_weight_measure}} >= {{lighter_weight_measure}})
  } else{
    message("No observations where heavier_weight_measure is less than lighter_weight_measure")
    return(data)
  } 

}

# data = biomass_data_taxaclean
# heavier_weight_measure <- expr(wm_g)
# lighter_weight_measure <- expr(dm_g)


taxa_notenough_obs <- function(preserved_method_data,
                               taxa_names){
  
  # Identify which ones I don't want
  notenough_taxa <- preserved_method_data |> 
    dplyr::group_by({{taxa_names}}) |> 
    dplyr::summarise(n = n()) |> 
    dplyr::filter(n < 10) |> 
    dplyr::pull({{taxa_names}})
  
  # Tells user which ones are exclude
  ntaxa <- length(notenough_taxa) # how many exclude
  y <- paste(notenough_taxa, collapse = ", ") # join taxa strings
  message(paste0(ntaxa, " taxa excluded: ",  y))
  
  # Exclude these from dataset
  output <- preserved_method_data |> 
    filter(! taxa_clean %in% notenough_taxa)
  
  return(output)
}





# https://dplyr.tidyverse.org/articles/programming.html#data-masking
# https://brad-cannell.github.io/r_notes/tidy-evaluation.html