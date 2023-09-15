library(tidyverse)
library(rlang)
library(janitor)

#' Excludes values under a certain cutoff mass
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

#' Excluding observations where dried samples are heavier than wet, or burnt are heavier than dried
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


#' Removing taxa with less than 10 observations
#'
#' @param preserved_method_data 
#' @param taxa_names 
#'
#' @return data including only taxa with at least 10 observations
#'
#' @examples
#' taxa_notenough_obs(step_3, taxa_clean)

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


#' Identifying observations outside given percentiles
#'
#' @param mass_measure 
#' @param lower_percentile
#' @param upper_percentile
#'
#' @return A logical vector indicating whether each observation of the specified mass measure is between the specified percentiles
#'
#' @examples
#' df |> mutate(keep_data_mass1_percentiles = biomass_percentiles(wm_g, 0.025, 0.975))

biomass_percentiles <- function(mass_measure, lower_percentile, upper_percentile){
    ifelse(mass_measure >= quantile(mass_measure, lower_percentile, na.rm = TRUE) 
         & mass_measure <= quantile(mass_measure, upper_percentile, na.rm = TRUE), TRUE, FALSE)
}


#' comparing leverage values to cutoff value, with inbuilt filter for taxa
#'
#' @param preservation_method_data 
#' @param taxon 
#' @param mass_measure_1 
#' @param mass_measure_2 
#'
#' @return logical vector indicating TRUE where leverage value is under the cutoff value for outlier designation, otherwise FALSE
#'
#' @examples
#' leverage_with_taxa(biomass_frozen_percentiles, "taxa_clean == 'Amphipoda'", "wm_g", "dm_g")

leverage_with_taxa <- function(preservation_method_data, taxon, mass_measure_1, mass_measure_2){
  # filtering for a particular taxon
  taxon <- rlang::parse_expr(taxon)
  data <- dplyr::filter(preservation_method_data, !!taxon)
  
  #running an OLS regression and calculating associated leverage values
  lm_object <- lm(eval(as.name(mass_measure_2)) ~ eval(as.name(mass_measure_1)), data = data)
  leverage_values <- hatvalues(lm_object)
  
  #calculating cutoff leverage value, based on Aguinis et al. (2013)
  predictors <- length(coef(lm_object))-1
  n <- nrow(lm_object$model)
  leverage_cutoff <- 2*(predictors+1)/n
  
  #creating logical vector indicating whether leverage value is under the cutoff
  output <- ifelse(leverage_values < leverage_cutoff, TRUE, FALSE)
  return(output)
}


#' comparing leverage values to cutoff value for outlier candidate designation
#'
#' @param preservation_method_data 
#' @param mass_measure_1 
#' @param mass_measure_2 
#'
#' @return logical vector indicating TRUE where leverage value is under the cutoff value, otherwise FALSE
#'
#' @examples
#' keep_data_leverage(biomass_frozen_percentiles, "wm_g", "dm_g")

keep_data_leverage <- function(preservation_method_data, mass_measure_1, mass_measure_2){
  #running an OLS regression and calculating associated leverage values
  lm_object <- lm(eval(as.name(mass_measure_2)) ~ eval(as.name(mass_measure_1)), data = preservation_method_data)
  leverage <- hatvalues(lm_object)
  
  #calculating cutoff leverage value, based on Aguinis et al. (2013)
  predictors <- length(coef(lm_object))-1
  n <- nrow(lm_object$model)
  leverage_cutoff <- 2*(predictors+1)/n
  
  #creating logical vector indicating whether leverage value is under the cutoff
  output <- tibble(leverage_values = leverage, keep_data_leverage = ifelse(leverage < leverage_cutoff, TRUE, FALSE))
  return(output)
}

#' comparing studentized deleted residuals to cutoff value
#'
#' @param preservation_method_data 
#' @param mass_measure_1 
#' @param mass_measure_2 
#'
#' @return logical vector indicating TRUE where studentized deleted residuals are within the cutoff range for outlier candidate designation, otherwise FALSE
#'
#' @examples
#' keep_data_stud_resid(biomass_frozen_percentiles, "wm_g", "dm_g")

keep_data_stud_resid <- function(preservation_method_data, mass_measure_1, mass_measure_2){
  #running an OLS regression and calculating associated studentized deleted residuals
  lm_object <- lm(eval(as.name(mass_measure_2)) ~ eval(as.name(mass_measure_1)), data = preservation_method_data)
  stud_resid <- rstudent(lm_object)
  
  #calculating cutoff value, based on Aguinis et al. (2013)
  predictors <- length(coef(lm_object))-1
  n <- nrow(lm_object$model)
  alpha_tdist <- 0.05/n
  df_tdist <- n-predictors-1
  critical_t <- qt(alpha_tdist/2, df_tdist, lower.tail = FALSE, log.p = FALSE)
  
  #creating logical vector indicating whether studentized deleted residual is within the cutoff range
  output <- ifelse(stud_resid < critical_t & stud_resid > -critical_t, TRUE, FALSE)
  return(output)
}


#' Comparing Cook's distance to cutoff value for outlier candidate designation
#'
#' @param preservation_method_data 
#' @param mass_measure_1 
#' @param mass_measure_2 
#'
#' @return logical vector indicating TRUE where leverage value is under the cutoff value, otherwise FALSE
#'
#' @examples
#' keep_data_cook(biomass_frozen_percentiles, "wm_g", "dm_g")

keep_data_cook <- function(preservation_method_data, mass_measure_1, mass_measure_2){
  #running an OLS regression and calculating associated cook's distance
  lm_object <- lm(eval(as.name(mass_measure_2)) ~ eval(as.name(mass_measure_1)), data = preservation_method_data)
  cooks_d <- cooks.distance(lm_object)
  
  #calculating cutoff value, based on Aguinis et al. (2013)
  predictors <- length(coef(lm_object))-1
  n <- nrow(lm_object$model)
  cooks_d_cutoff <- qf(0.5, predictors+1, n-predictors-1, lower.tail = FALSE)
  
  #creating logical vector indicating whether cook's distance is under the cutoff
  output <- ifelse(cooks_d < cooks_d_cutoff, TRUE, FALSE)
  return(output)
}


#' Comparing dffits value to cutoff value for outlier candidate designation
#'
#' @param preservation_method_data 
#' @param mass_measure_1 
#' @param mass_measure_2 
#'
#' @return logical vector indicating TRUE where dffits values are within the cutoff range for outlier candidate designation, otherwise FALSE
#'
#' @examples
#' keep_data_dffits(biomass_frozen_percentiles, "wm_g", "dm_g")

keep_data_dffits <- function(preservation_method_data, mass_measure_1, mass_measure_2){
  #running an OLS regression and calculating associated dffits values
  lm_object <- lm(eval(as.name(mass_measure_2)) ~ eval(as.name(mass_measure_1)), data = preservation_method_data)
  dffits <- dffits(lm_object)
  
  #calculating cutoff value, based on Aguinis et al. (2013)
  predictors <- length(coef(lm_object))-1
  n <- nrow(lm_object$model)
  dffits_cutoff <- 2*sqrt((predictors+1)/n)
  
  #creating logical vector indicating whether dffits value is within the cutoff range
  output <- ifelse(dffits < dffits_cutoff & dffits > -dffits_cutoff, TRUE, FALSE)
  return(output)
}


#' Comparing dfbetas value to cutoff value for outlier candidate designation
#'
#' @param preservation_method_data 
#' @param mass_measure_1 
#' @param mass_measure_2 
#'
#' @return list of two logical vectors, one for the intercept and one for the slope, indicating TRUE where dfbetas values are within the cutoff range for outlier candidate designation, otherwise FALSE
#'
#' @examples
#' keep_data_dfbetas(biomass_frozen_percentiles, "wm_g", "dm_g")
# keep_data_dfbetas(biomass_frozen_percentiles, "wm_g", "dm_g")
keep_data_dfbetas <- function(preservation_method_data, mass_measure_1, mass_measure_2){
  #running an OLS regression and calculating associated dfbetas values
  lm_object <- lm(eval(as.name(mass_measure_2)) ~ eval(as.name(mass_measure_1)), data = preservation_method_data)
  dfbetas <- dfbetas(lm_object)
  
  #calculating cutoff value, based on Aguinis et al. (2013)
  n <- nrow(lm_object$model)
  dfbetas_cutoff <- 2/sqrt(n)
  
  #creating logical vector indicating whether dfbetas values are within the cutoff range
  output <- ifelse(dfbetas < dfbetas_cutoff & dfbetas > -dfbetas_cutoff, TRUE, FALSE)
  return(output)
}


#Trying to make a function for checking which observations of a particular taxa are identified by different measures of outlyingess
outlier_check <- function(preservation_method_data, taxon, outlier_measure, keep_data_outlier_measure){ 
  taxon <- rlang::parse_expr(taxon)
  data <- dplyr::filter(preservation_method_data, !!taxon)
  data |>
  select(outlier_measure, keep_data_outlier_measure, obs_num) |>
  arrange(desc(outlier_measure))|> 
  print(n = Inf) |>
  tabyl(keep_data_outlier_measure)
}
# outlier_check(taxa_unnested_outlier_candidates_identified, "taxa_clean == 'Amphipoda'", leverage_values, keep_data_leverage)


# https://dplyr.tidyverse.org/articles/programming.html#data-masking
# https://brad-cannell.github.io/r_notes/tidy-evaluation.html
# https://stackoverflow.com/questions/72356097/how-to-pass-dplyr-filter-a-character-vector-as-part-of-a-function
# https://stackoverflow.com/questions/52856711/use-function-arguments-in-lm-formula-within-function-environment?rq=3
