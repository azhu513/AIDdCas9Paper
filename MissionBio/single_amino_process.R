## ----------------------------------------------------
## Script name: single_amino_process.R
## Purpose of script: process the processed data files
##                    of single cell amplicon-seq to
##                    have single amino level informat-
##                    ion
## Author: Anqi Zhu
## Date Created: 2023-04-25
## ----------------------------------------------------
## Input: sample_", sample_id, ".txt.annotated" from
##        running data_process.R and annotate_raw_data
## Output:    data frame of single amino acids change
##            at single cell level
## ----------------------------------------------------
## Load up the packages -------------------------------

suppressMessages(require(dplyr))
suppressMessages(require(tidyr))
suppressMessages(require(parallel))
suppressMessages(require(tidyverse))

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
sample_id <- args[1] # sample id

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

dat_annotate <- read.table(file = paste0("scratch/data/", set_id, "/sample_", sample_id, ".txt.annotated"), sep="\t", header = F)

data_all_full <- read.table(file = paste0("scratch/data/", set_id, "/data_all_full_", sample_id, ".txt"), sep="\t", header = T)

col_names = c("Target", "Position", "Allele", "Allele_Count", "Allele_Freq", "Total_Count", "Allele_Seq", "Mut_Pos", "Ref_Seq", "cell_barcode", "allele_count_new",	"allele_freq_new", "WT_Count",	"WT_Freq", "Location", "Gene", "Consequence", "Amino_acids", "IMPACT", "protein_position")

colnames(dat_annotate) <- col_names

dat_total_cells <- data_all_full |>
  ungroup() |>
  group_by(Target) |>
  summarize(n_total_cells = n_distinct(cell_barcode)) |>
  ungroup()

dat_frame <- dat_annotate |>
  filter(grepl(pattern = 'frameshift_variant|inframe_deletion|inframe_insertion|protein_altering_variant|protein_altering_variant', x = Consequence)) |>
  mutate(amino_change = paste(protein_position, Amino_acids),
         amino_count = allele_count_new) |>
  filter(amino_count >= amino_count_threshold) |>
  select(Target, Position, Allele, allele_count_new, Total_Count, cell_barcode, amino_change, Consequence, IMPACT, amino_count)

dat <- dat_annotate |>
  filter(Amino_acids!='-' & Amino_acids!='' & protein_position!='-' & Allele != 'WT')  |>
  filter(!grepl(pattern = 'frameshift_variant|inframe_deletion|inframe_insertion|protein_altering_variant|protein_altering_variant', x = Consequence) & !grepl(pattern = 'synonymous_variant', x = Consequence)) |>
  separate(Amino_acids, sep = "/", into = c("before", "after")) |>
  separate(protein_position, sep = "-", into = c("start_pos", "end_pos"))

dat_single <- dat |>
  filter(is.na(end_pos))

dat_mult <- dat |>
  filter(!is.na(end_pos))

dat_mult_split = list()

for (i in seq_len(nrow(dat_mult))){
  before_new = unlist(strsplit(as.character(dat_mult[i, "before"]), ""))
  after_new = unlist(strsplit(as.character(dat_mult[i, "after"]), ""))

  max_length = max(length(before_new), length(after_new))
  before_new_padded = c(before_new, rep(NA, max_length - length(before_new)))
  after_new_padded = c(after_new, rep(NA, max_length - length(after_new)))

  dat_mult_split[[i]] = data.frame(
    Target = dat_mult[i, "Target"],
    Position = dat_mult[i, "Position"],
    Allele = dat_mult[i, "Allele"],
    cell_barcode = dat_mult[i, "cell_barcode"],
    before_new_padded, after_new_padded)
}

dat_mult_df = do.call(rbind, dat_mult_split)

dat_mult_pos <- dat_mult |>
  mutate(positions = map2(start_pos, end_pos, seq)) |>
  unnest(positions)

dat_mult_df = dat_mult_df |>
  group_by(Target, Position, Allele, cell_barcode)
key1 = group_keys(dat_mult_df)
dat_mult_df = dat_mult_df |>
  group_split()

dat_mult_pos = dat_mult_pos |>
  group_by(Target, Position, Allele, cell_barcode)
key2 = group_keys(dat_mult_pos)
dat_mult_pos = dat_mult_pos |>
  group_split()

all(key1 == key2)

align_list = list()
for (i in seq_len(nrow(key1))) {
  if (nrow(dat_mult_pos[[i]])==nrow(dat_mult_df[[i]])){
    pos = dat_mult_pos[[i]] |>
      select(-c(before,after))
    selected = dat_mult_df[[i]] |>
      select(before = before_new_padded, after = after_new_padded)
    align_list[[i]] = cbind(pos, selected)
  }
}
dat_mult_aligned_df = do.call(rbind, align_list)

final_1 = dat_mult_aligned_df |>
  select(Target, Position, Allele, allele_count_new, Total_Count, cell_barcode, positions, before, after, Consequence, IMPACT) |>
  mutate(after = ifelse(is.na(after), "", after)) |>
  filter(before!=after)

final_2 = dat_single |>
  select(Target, Position, Allele, allele_count_new, Total_Count, cell_barcode, positions=start_pos, before, after, Consequence, IMPACT) |>
  mutate(after = ifelse(is.na(after), "", after))

dat_ <- rbind(final_1, final_2) |>
  mutate(amino_change = paste0(positions, " ", before, "/", after)) |>
  group_by(Target, Position, cell_barcode, amino_change) |>
  mutate(amino_count = sum(allele_count_new)) |>
  filter(amino_count >= amino_count_threshold) |>
  distinct(Target, Position, cell_barcode, amino_change, amino_count, .keep_all = TRUE) |>
  select(-c(before, after, positions)) |>
  ungroup()

count_cell <- rbind(dat_, dat_frame) |>
  group_by(Target, amino_change) |>
  mutate(n_amino_cells = n_distinct(cell_barcode))  |>
  ungroup() |>
  distinct(Target, amino_change, .keep_all = TRUE) |>
  left_join(dat_total_cells, by = "Target")

dat_full <- rbind(dat_, dat_frame)

write.table(count_cell, file = paste0("scratch/data/", set_id, "/single_amino_cell_count_sample_", sample_id, ".txt"), row.names=FALSE, sep="\t", quote = FALSE)

write.table(dat_full, file = paste0("scratch/data/", set_id, "/single_amino_sample_", sample_id, ".txt"), row.names=FALSE, sep="\t", quote = FALSE)

write.table(dat_total_cells, file = paste0("scratch/data/", set_id, "/single_amino_sample_total_cells_", sample_id, ".txt"), row.names=FALSE, sep="\t", quote = FALSE)
