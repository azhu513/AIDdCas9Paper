## ----------------------------------------------------
## Script name: data_process.R
## Purpose of script: process the txt files of single cell
##                    amplicon-seq from crisprscan workflow
## Author: Anqi Zhu
## Date Created: 2023-04-25
## ----------------------------------------------------
## Input: cell_barcode.txt from running crisprscan
## Output:    data frame digested from all cells from a
##            sample with cell barcode attached,
##            valid snv alleles in a separated data frame,
##            summary information of total reads per cell
##            per amplicon
## ----------------------------------------------------

## Load up the packages -------------------------------

suppressMessages(require(dplyr))
suppressMessages(require(tidyr))
suppressMessages(require(parallel))
suppressMessages(require(data.table))

## ----------------------------------------------------

## load up our functions into memory ------------------

source("./functions/utils.R")

## ----------------------------------------------------

## Batch submission configurations --------------------

# specifies that you will be manually setting the cores
options(future.availablecores.methods = "mc.cores")
# SLURM_CPUS_PER_TASK is the amount of cores specified in the job environment
options(mc.cores = Sys.getenv("SLURM_CPUS_PER_TASK"))

# get the argument from command line
args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]

## thresholds for data cleaning-----------------------
thresholds <- list("allele_count_threshold" = 3,
                           "allele_freq_threshold" = 0.05,
                           "total_count_threshold" = 10,
                           "mutation_count_threshold" = 2,
                           "amino_count_threshold" = 3,
                           "bulk_read_min_count"= 10,
                           "bulk_read_min_sample" = 3,
                           "sc_min_count"= 5,
                           "sc_min_sample"=3)

allele_count_threshold <- thresholds[['allele_count_threshold']]
allele_freq_threshold <- thresholds[['allele_freq_threshold']]
total_count_threshold <- thresholds[['total_count_threshold']]
mutation_count_threshold <- thresholds[['mutation_count_threshold']]
## ----------------------------------------------------

# arguments for reading the txt file
remove_idx = c(1,2,3,9,10)
col_names = c("Target", "Position", "Allele", "Allele_Count", "Allele_Freq", "Total_Count", "Allele_Seq", "Mut_Pos", "Ref_Seq")

CNV_Dup_amplicons <- "^EGFR_|^hEGFR_522$"

# Function to process MissionBio single cell data to remove low frequency alleles, and remove amplicons with total read counts lower than the threshold specified
# Also, split the alleles into single mutations if there are multiple mutations on the same allele and sum the read counts up by mutation
#
# Note that for the simplicity of the code I baked in column names that only works for my column names for MissionBio data, change it when use it on another set of data

ProcessSingleCellFile <- function(file_name,
                                  col_remove, col_names,
                                  amplicon_include,
                                  amplicon_exclude,
                                  amplicon_name,
                                  allele_threshold,
                                  amplicon_threshold,
                                  split_col_name,
                                  split_delim,
                                  filter_col_name_allele,
                                  filter_col_name_amplicon,
                                  allele_sum_value_count,
                                  allele_sum_value_freq,
                                  mut_count_threshold,
                                  joinby,
                                  ...){

  # read the text file per cell, remove unnecessary columns, and assign column names
  # also, remove NTC and CNV amplicons
  dat_all <- ReadSingleCRISPRscan(file_name, includeCellBarcode = TRUE, col_remove, col_names)

  dat <- dat_all |>
    filter(!grepl(amplicon_exclude, {{ amplicon_name }})) |>
    mutate({{ split_col_name }} := sub(pattern = "\\|$", replacement="", x = {{ split_col_name }}))

  # count WT alleles per amplicon target region and append the count and freq to the two data frames
  WT_count_freq <- dat |>
    filter({{ split_col_name }}=="WT") |>
    select({{ amplicon_name }}, {{ allele_sum_value_count }}, {{ allele_sum_value_freq }}) |>
    rename(WT_Count = {{ allele_sum_value_count }}, WT_Freq = {{ allele_sum_value_freq }})

  # split alleles into single mutations
  singlemut_dat <- tidyr::separate_rows(dat, {{ split_col_name }}, sep = split_delim)

  # process singlemut_dat and sum up the read counts by amplicon and allele
  singlemut_dat <- singlemut_dat |>
    group_by({{ amplicon_name }}, {{ split_col_name }}) |>
    mutate(mutation_count = sum({{ allele_sum_value_count }})) |>
    filter( {{ filter_col_name_amplicon }} > amplicon_threshold) |>
    distinct({{ amplicon_name }}, {{ split_col_name }}, mutation_count, .keep_all = TRUE) |>
    ungroup()

  # the invalid single nt mutations are those with low read counts
  invalid_mutations_df <- singlemut_dat |>
    filter(mutation_count < mut_count_threshold)

  amps <- invalid_mutations_df |> distinct({{ amplicon_name }})

  # get the valid mutations and form the final OUTPUT of singlemut data
  singlemut_dat <- singlemut_dat |>
    left_join(WT_count_freq, by = joinby) |>
    filter(mutation_count >= mut_count_threshold)

  # process dat
  # remove those noise mutations fron invalid snv
  # consolidate alleles and counts/frequencies
  # and then remove low count alleles and amplicons with low total counts

  for (i in seq_len(nrow(amps))){
    amp <- as.character(amps[i,])
    muts <- invalid_mutations_df |> filter({{ amplicon_name }} == amp) |> distinct({{ split_col_name }})
    for (j in seq_len(nrow(muts))) {
      mut <- as.character(muts[j,])

      dat <- dat |>
        mutate({{ split_col_name }} := ifelse({{ amplicon_name }} == amp, gsub(pattern = mut, replacement="", x = {{ split_col_name }}), {{ split_col_name }})) |>
        mutate({{ split_col_name }} := sub(pattern = "\\|$", replacement="", x = {{ split_col_name }})) |>
        mutate({{ split_col_name }} := sub(pattern = "^\\|", replacement="", x = {{ split_col_name }})) |>
        filter({{ split_col_name }} != "")
    }
  }

  dat <- dat |>
    group_by({{ amplicon_name }}, {{ split_col_name }}) |>
    mutate(allele_count_new = sum({{ allele_sum_value_count }}),
           allele_freq_new = sum({{ allele_sum_value_freq }})) |>
    filter(allele_count_new >= allele_threshold ) |>
    filter({{ filter_col_name_amplicon }} > amplicon_threshold) |>
    ungroup() |>
    left_join(WT_count_freq, by = joinby) |>
    distinct({{ amplicon_name }}, {{ split_col_name }}, allele_count_new, allele_freq_new, .keep_all = TRUE)


  return(list(singlemut_dat,
              dat,
              dat_all))
}

defaultFun <- function(file_name){
  res <- file_name |> ProcessSingleCellFile(col_remove = remove_idx, col_names = col_names,
                                            amplicon_include = NA,
                                            amplicon_exclude = CNV_Dup_amplicons,
                                            amplicon_name = Target,
                                            allele_threshold = allele_count_threshold,
                                            amplicon_threshold = total_count_threshold,
                                            split_col_name = Allele,
                                            split_delim = "\\|",
                                            filter_col_name_allele = Allele_Count,
                                            filter_col_name_amplicon = Total_Count,
                                            allele_sum_value_count = Allele_Count,
                                            allele_sum_value_freq = Allele_Freq,
                                            mut_count_threshold = mutation_count_threshold,
                                            joinby = "Target")

  res
}

root_dir = "crisprbe/All_crispr-scan_full"
dirs <- list.dirs(root_dir, recursive= FALSE)
matching_dirs <- dirs[grep(pattern = sample_id, x = dirs)]
file_names <- list.files(matching_dirs, full.names = TRUE)

res_list <- mclapply(file_names, FUN = defaultFun)

data_full <- do.call(bind_rows, lapply(res_list, "[[", 2)) # which are not full full they are processed with alleles and amplicons removed

singlemut_full <- do.call(bind_rows, lapply(res_list, "[[", 1)) # which are not full full they are processed with alleles and amplicons removed
data_all_full <- do.call(bind_rows, lapply(res_list, "[[", 3)) # real full data, nothing filtered

target_cell_total_reads <- data.frame(data_all_full) |>
  distinct(Target, cell_barcode, Total_Count)

## create summary data frame group by cell and amplicon target region
# number of alleles, valid alleles by count threshold and allele frequency
#
cell_total <- data.frame(data_full) |>
  distinct(Target, cell_barcode, .keep_all = TRUE) |>
  group_by(cell_barcode) |>
  summarize(cell_total_read = sum(Total_Count))

Target_total_reads <- data.frame(data_full) |>
  distinct(Target, cell_barcode, .keep_all = TRUE) |>
  group_by(Target) |>
  summarize(target_total_read = sum(Total_Count))

summary_count <- data.frame(data_full) |>
  left_join(cell_total, by = c('cell_barcode')) |>
  filter(allele_count_new >= allele_count_threshold ) |>
  filter(Total_Count > total_count_threshold) |>
  group_by(cell_barcode) |>
  mutate(n_amplicon = n_distinct(Target),
         amplicon_readcount_proportion = Total_Count/cell_total_read) |>
  group_by(Target, .add=TRUE) |>
  mutate(n_allele_total = n_distinct(Allele),
         n_allele_valid_by_count = sum(allele_count_new> allele_count_threshold),
         n_allele_valid_by_freq = sum(allele_freq_new > allele_freq_threshold),
         n_allele_valid_by_both = ifelse(Total_Count <= 60, n_allele_valid_by_freq, n_allele_valid_by_count),
         max_count = max(allele_count_new),
         max_freq =  max(allele_freq_new),
         all_WT = (n_allele_valid_by_both == 1 & ((Total_Count <= 60 & WT_Freq == max_freq) | (Total_Count > 60 & WT_Count == max_count)))) |>
  ungroup()

count_cell <- singlemut_full |>
  filter(Total_Count > total_count_threshold) |>
  ungroup() |>
  group_by(Target) |>
  mutate(n_total_cells = n_distinct(cell_barcode)) |>
  group_by(Allele, .add = TRUE) |>
  mutate(n_allele_cells = n_distinct(cell_barcode))  |>
  ungroup() |>
  distinct(Target, Allele, .keep_all = TRUE) |>
  select(-c(cell_barcode, Allele_Count, Allele_Freq, Total_Count, mutation_count, WT_Count, WT_Freq)) |>
  arrange(Target, Allele)

write.table(data_full, file = paste0("scratch/data/", set_id, "/sample_", sample_id, ".txt"), row.names=FALSE, sep="\t", quote = FALSE)

write.table(singlemut_full, file = paste0("scratch/data/", set_id, "/single_mut_sample_", sample_id, ".txt"), row.names=FALSE, sep="\t", quote = FALSE)

write.table(target_cell_total_reads, file = paste0("scratch/data/", set_id, "/target_cell_total_count_", sample_id, ".txt"), row.names=FALSE, sep="\t", quote = FALSE)

write.table(data_all_full, file = paste0("scratch/data/", set_id, "/data_all_full_", sample_id, ".txt"), row.names=FALSE, sep="\t", quote = FALSE)

save(summary_count, count_cell, file = paste0("data/", set_id, "/MissionBio_data_summary_sample_", sample_id, ".RData"))

## After this step, run MissionBio_annotate_raw_data.sh
