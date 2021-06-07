# Metadata_Prep_Functions

#----------------------------------------------------------------------------
##--------------            Metadata lists             --------------
#----------------------------------------------------------------------------

motortest_vars <- c(
  "pole_test_s",
  "sticker_test_s",
  "x_of_steps",
  "errors_step",
  "beam_time_to_cross",
  "wire_hang_s",
  "hindlimb_score")

grouping_vars <- c(
  "cage_id",
  "diet",
  "genotype",
  "treatment_group",
  "description"
)

mouse_vars <- c(
  "diversity_shannon",
  "host_weight",
  "collection_timestamp",
  "donor_id"
)


#----------------------------------------------------------------------------
##--------------           METADATA PREP FUNCTION           --------------
#----------------------------------------------------------------------------

process_meta <- function(dat, verbose = F){

  env <- microbiome::meta(dat) %>% 
    dplyr::mutate(mouse_id = as.factor(mouse_id)) %>%
    dplyr::mutate(collection_timestamp = as.character(collection_timestamp)) %>% 
    dplyr::mutate(collection_timestamp = 
                    substr(collection_timestamp, 1, nchar(collection_timestamp)-6 )) %>% 
    dplyr::mutate(collection_timestamp = paste0(
      substr(collection_timestamp, 1, nchar(collection_timestamp)-2), "20",
      substr(collection_timestamp, nchar(collection_timestamp)-1, nchar(collection_timestamp)) ))  %>% 
    dplyr::mutate(collection_timestamp = as.Date(collection_timestamp, format = "%m/%d/%Y")) 
  
  if (verbose == T){
    print(str(env))
  }
  cat ("Metadata Processing complete \n")
  # load env into global enviornment
  assign("env",env,envir = .GlobalEnv)
  return(env)
}

#----------------------------------------------------------------------------
##--------------           TRIM METADATA FUNC           --------------
#----------------------------------------------------------------------------

trim_meta <- function(env, min_ratio){
  if (typeof(env) != "list"){
    stop("Please input preprocessed metdata list \n RUN process_meta() function first")
  }
  ###  Remove (Yes/No) questions with less than a certain % (min ratio) of samples as minority
  rmvlst_index <- c()
  rmvlst_names <- c()
  cat("List of Metadata below filtering threshold: \n")
  for (i in 1:length(env)) {
    p <- na.omit(env[, i])
    x <- count(p)
    if (length(x$x) == 2) {
      a <- min(x$freq)
      b <- sum(x$freq)
      if ((a / b) < min_ratio) {
        print(colnames(env)[i])
        rmvlst_names <- c(rmvlst_names, colnames(env)[i])
        rmvlst_index <- c(rmvlst_index, i)
      }
      ###  Remove questions with only a single factor response
    } else if (length(x$x) == 1) {
      print(colnames(env)[i])
      rmvlst_names <- c(rmvlst_names, colnames(env)[i])
      rmvlst_index <- c(rmvlst_index, i)
    }
  }
  if (!is.null(rmvlst_index)) {
    env <- env[-rmvlst_index]
  }
  cat("\nMetadata Filtering Complete: \n\n\n")
  assign("env",env,envir = .GlobalEnv)
}


