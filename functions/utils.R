## ------------------------------------------------
## Script name: utils.R
## Purpose of script: data processing and cleaning 
##                    functions functions used 
##                    throughout the the analysis
## Author: Anqi Zhu
## Date Created: 2023-05-26
## ------------------------------------------------
## Notes:
##   ReadSingleCRISPRscan - read single cell txt 
##   NumberCells - count number of cells original
##   Aggregate_Sample - aggregate samples
##   Split_Column - split the column of specified
##   logLikNB_more - check log likelihood function 
##                   of NBinom
## ------------------------------------------------

## Load up the packages ---------------------------

suppressMessages(require(stringr))
suppressMessages(require(dplyr))

## ------------------------------------------------

# This function reads one txt file and output a data table of all rows in the txt
ReadSingleCRISPRscan <- function(filename, includeCellBarcode = TRUE, col_remove, col_names) {
  if (includeCellBarcode) {
    CellBarcode <- gsub(".txt", "", unlist(strsplit(filename,split = "/"))[8])
    dat <- read.table(filename, header = FALSE, sep = "\t", col.names = paste0("V",seq_len(14)), fill = TRUE)
    dat <- dat[, -col_remove]
    colnames(dat) <- col_names
    dat[, "cell_barcode"] <- CellBarcode
  } else {
    dat <- read.table(filename, header = FALSE, sep = "\t", col.names = paste0("V",seq_len(14)), fill = TRUE)
    dat <- dat[, -col_remove]
    colnames(dat) <- col_names
  }
  dat
}

# This function counts the number of files within a folder that is identified by sample id
NumberCells <- function(root_dir, sample_id, ...){
  dirs <- list.dirs(root_dir, recursive= FALSE)
  matching_dirs <- dirs[grep(sample_id, dirs)]
  num_files <- length(list.files(matching_dirs))
  return(num_files)
}

# aggregate single cell files into sample level information
# parallel computing! usually use mc.cores specified by bash script
Aggregate_Sample <- function(root_dir, sample_id, col_remove, col_names, mc.cores = NA, ...){
  dirs <- list.dirs(root_dir, recursive= FALSE)
  matching_dirs <- dirs[grep(sample_id, dirs)]
  file_names <- list.files(matching_dirs, full.names = TRUE)
  
  if (!is.na(mc.cores)){
    res_list <- mclapply(file_names, FUN = ReadSingleCRISPRscan, col_remove = col_remove, col_names = col_names, mc.cores=mc.cores)
  }
  
  res_list <- mclapply(file_names, FUN = ReadSingleCRISPRscan, col_remove = col_remove, col_names = col_names)
  
  res_df <- do.call(rbind, res_list)
  res_df
}

# a util function to split column
Split_Column <- function(data_frame, col_name, delim){
  # delim required to be a regex
  library(tidyr)
  
  df_split <- tidyr::separate_rows(data_frame, {{ col_name }}, sep = delim)
  df_split
}

# check log likelihood function of NBINOM
logLikNB_more <- function(y, x, coef, beta, beta_other, param, offset) {
  if (is.numeric(coef)) coef=as.integer(coef)
  if (!is.null(beta_other)){
    stopifnot(!is.null(coef))
    if (!is.integer(coef)){
      coef <- match(coef, colnames(x))
    }
    xbeta <- x[,coef]*beta + x[,-coef] %*% beta_other + offset
  } else {
    xbeta <- x %*% beta + offset    
  }
  mean.hat <- exp(xbeta)
  dnbinom(y, mu=mean.hat, size=1/param, log=TRUE)
}

