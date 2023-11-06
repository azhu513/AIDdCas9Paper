## ------------------------------------------------
## Script name: model_helper.R
## Purpose of script: To process the data in this 
##                    project and apply the DESeq2
##                    workflow and apeglm shrinkage
##                    estimates
##
## Author: Anqi Zhu
## Date Created: 2023-05-31
## ------------------------------------------------
## Notes: 
##   imputeTR - impute missing total reads
##   check_imputeTR - check imputation result
##   GetDDS - form data objects for DESeq2 model 
##             and apply DESeq2 and apeglm and 
##             output results 
##
## ------------------------------------------------

## Load up the packages ---------------------------

suppressMessages(library(DESeq2))
suppressMessages(library(apeglm))

## ------------------------------------------------

# impute missing total reads
imputeTR <- function(x, name_pattern){
  col_idx <- grep(name_pattern, colnames(x))
  x_new <- x[, col_idx]
  out <- apply(x_new, 2, FUN = function(y) ifelse(y==0, max(y), y))
  colnames(out) <- colnames(x)[col_idx]
  out
}

# check imputation result of 'imputeTR'
check_imputeTR <- function(m) all(apply(m, 2, function(x) length(unique(x)) == 1))

# form data objects for DESeq2 model and apply DESeq2 and apeglm
GetDDS <- function(data, 
                   trt_var, 
                   trt_ref_level, 
                   other_cond_vars, 
                   row_group_var,
                   select_cols, 
                   sample_id_vars, 
                   value_vars, 
                   ct_name, 
                   tr_name, 
                   read_min_count, 
                   read_min_sample, 
                   collapse = FALSE,
                   collapse_var = NA, ...){
  
  newdata <- data 
  
  # get the columns indexes of selected columns
  column_names <- colnames(newdata)
  select_cols_index <- match(select_cols, column_names)
  
  # change the data to wide format
  dat_wide <- newdata  |> 
    select(all_of(select_cols_index)) |> 
    pivot_wider(names_from = sample_id_vars,
                names_sep = "_",
                values_from = value_vars, 
                values_fill = 0) |>  
    arrange(.data[[row_group_var]])
  
  # count matrix, imputations are still 0 assuming no detection means no count
  cts_full <- dat_wide[, grep(pattern = ct_name, x = colnames(dat_wide))]
  
  # The total reads should not be imputed with 0, 
  # instead, should be imputed by total reads of that target in that sample

    dat_wide_split <- dat_wide %>% ungroup() %>% group_by_at({{row_group_var}})
    trs_full_split_key <- group_keys(dat_wide_split)
    trs_full_split_ls <- group_split(dat_wide_split) 
    trs_full_split_imputed <- lapply(trs_full_split_ls, FUN = imputeTR, name_pattern = tr_name)
    
    ## Sanity checking
    if(!(all(sapply(trs_full_split_imputed, check_imputeTR)))) {
      stop("Error: total count imputation is not correct, some group/sample has different numbers of total count!")
    }
  
    row_group_var_key <- unique(dat_wide[,row_group_var])
    
    if (!all(row_group_var_key==trs_full_split_key)){
      stop("Error: number of groups do not match ")
    }
    
    trs_full <- do.call(rbind, trs_full_split_imputed)
    
  # get the row data for DESeq object
  rowdata <- dat_wide[,!(grepl(pattern = ct_name, x = colnames(dat_wide)) | grepl(tr_name, x = colnames(dat_wide)))]
  
  # get the coldata (or experiment design)
  if (!is.na(other_cond_vars)){
    if (!is.na(collapse_var)){
      cond_col_index <- match(c(sample_id_vars, collapse_var, other_cond_vars, trt_var), column_names)
    } else {
      cond_col_index <- match(c(sample_id_vars, other_cond_vars, trt_var), column_names)
  }
  } else {
    if (!is.na(collapse_var)){
      cond_col_index <- match(c(sample_id_vars, collapse_var, trt_var), column_names)
    } else {
      cond_col_index <- match(c(sample_id_vars, trt_var), column_names)
    }
  }
  
  conditions_df <- newdata %>% 
    select(all_of(cond_col_index)) %>%
    distinct(.keep_all = TRUE) 
  
  conditions_df <- data.frame(conditions_df)
  
  pos_index <- length(unlist(strsplit(colnames(cts_full)[1], "_")))
  
  cts_cols_id <- sapply(strsplit(colnames(cts_full), "_"), "[", pos_index)
  trs_cols_id <- sapply(strsplit(colnames(trs_full), "_"), "[", pos_index)
  
  if (!all(cts_cols_id==trs_cols_id)){
    stop("Column orders of count matrix and total read matrix are different!")
  }
  
  cond_sample_id <- pull(.data=conditions_df, 1)
  conditions_df <- conditions_df[match(cts_cols_id, cond_sample_id),] 
  coldata <- as.data.frame(conditions_df, stringsAsFactors = TRUE)
  coldata[[trt_var]] <- relevel(factor(coldata[[trt_var]]), ref=trt_ref_level)
  
  if (!is.na(other_cond_vars)) {
  coldata[[other_cond_vars]] <- factor(coldata[[other_cond_vars]])
  }
  
  if (!is.na(other_cond_vars)){
  dds <- DESeqDataSetFromMatrix(countData=cts_full, colData=coldata, design =  as.formula(paste0("~", paste(c(other_cond_vars, trt_var), collapse="+"))))
  } else {
  dds <- DESeqDataSetFromMatrix(countData=cts_full, colData=coldata, design =  as.formula(paste0("~",  trt_var))) 
  }
  
  mcols(dds) <- rowdata
  
  if (collapse) {
    
    stopifnot(!is.na(collapse_var))
    
    dds <- collapseReplicates( dds,
                               groupby = dds[[collapse_var]],
                               run = dds[[sample_id_vars]])
    
    coldata[[collapse_var]] <- factor(coldata[[collapse_var]])
    groupby = coldata[[collapse_var]]
    sp <- split(seq(along=groupby), groupby)
    
    trs_collapse <- sapply(sp, function(i) MatrixGenerics::rowSums(trs_full[, i, drop=FALSE]))
    
    mode(trs_collapse) <- "integer"

    norm_factors <- (trs_collapse+1)/ exp(rowMeans(log(trs_collapse+1)))
    
  } else {
    norm_factors <- (trs_full+1)/ exp(rowMeans(log(trs_full+1)))
  }
  
  normalizationFactors(dds) <- norm_factors
  
  idx <- rowSums(counts(dds, normalized=TRUE) >= read_min_count ) >= read_min_sample
  
  if (nrow(dds) < 50) {
    dds <- DESeq(dds[idx, ], fitType="mean", betaPrior = FALSE)
  } else {
    dds <- DESeq(dds[idx, ], betaPrior = FALSE)
  }
  res <- results(dds)
  resLFC <- lfcShrink(dds, coef = tail(resultsNames(dds), 1), res=res, type="apeglm", svalue=TRUE, returnList=TRUE)
  
  
  return(list("WideData" = dat_wide, 
              "dds" = dds,
              "res" = res,
              "resLFC" = resLFC))
}