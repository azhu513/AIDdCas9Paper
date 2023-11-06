## ----------------------------------------------------
## Script name: sc_cell_count_model.R
## Purpose of script: model the cell count instead of
##                    read count in Bayes GLM
## Author: Anqi Zhu
## Date Created: 2023-04-25
## ----------------------------------------------------
## Notes: calculate fold change in cell population with
##        single mutation level data
##
## ----------------------------------------------------

## Load up the packages -------------------------------

suppressMessages(require(dplyr))
suppressMessages(require(tidyr))

## ----------------------------------------------------

## load up our functions into memory ------------------

source("./functions/model_helper.R")

## ----------------------------------------------------

## Batch submission configurations --------------------

# specifies that you will be manually setting the cores
options(future.availablecores.methods = "mc.cores")
# SLURM_CPUS_PER_TASK is the amount of cores specified in the job environment
options(mc.cores = Sys.getenv("SLURM_CPUS_PER_TASK"))

## thresholds -----------------------------------------
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
amino_count_threshold <- thresholds[['amino_count_threshold']]
read_min_count <- thresholds[['sc_min_count']]
read_min_sample <- thresholds[['sc_min_sample']]

## ----------------------------------------------------

# sample information
samples <- c(paste0("0", seq(1,9,1)), seq(10, 16))
cond_df <- data.frame(samID = samples,
                      condition = c(rep('uninfected', 4), rep("DOX", 6), rep("Erl", 3), rep("Osi", 3)),
                      sample_pair = c(0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 1, 2, 3, 1, 2, 3))

# gather data
data_long <- list()
for (i in seq_along(samples)){
  sample_id <- samples[i]
  dat <- read.table(paste0("scratch/data/", set_id, "/single_amino_cell_count_sample_", sample_id, ".txt"), sep = "\t", header = TRUE)

  full <- dat |>
    ungroup() |>
    # select(Target, Position, Allele, amino_change, n_amino_cells, n_total_cells) |>
    distinct(Target, amino_change, .keep_all = TRUE) |>
    mutate(samID = sample_id)

  data_long[[i]] <- full
  rm(dat, full)
}

data_long_df <- do.call(rbind, data_long)
rm(data_long)

# add sample information
data_long_df <- data_long_df |> left_join(cond_df, by = "samID")

# remove extremely low coverage regions

data_long_df <- data_long_df |>
  filter(Target != 'hBARF_69_OT1' & Target != 'hBRAF_69_OT5' & Target != 'hEGFR_529_OT2' & Target != 'hEGFR_522_OT7'& Target != 'hEGFR_522_OT3' & Target != 'hBRAF_48_OT1' )

# identify technical replicate sequencing samples
data_long_df$sample <- paste(data_long_df$condition, data_long_df$sample_pair)

# in this analysis, we will use the single mutation data from split full

data_Erl <- data_long_df |> filter(samID %in% c('05', '07', '09', '06', '08', '10', '11', '12', '13'))
data_Osi <- data_long_df |> filter(samID %in% c('05', '07', '09', '06', '08', '10', '14', '15', '16'))

DDS_ls_Erl <- GetDDS(data = data_Erl, trt_var = "condition", trt_ref_level = "DOX", other_cond_vars = "sample_pair", row_group_var = "Target", select_cols = c("samID", "Target", "Position", "amino_change", "n_amino_cells", "n_total_cells"), sample_id_vars = "samID", value_vars = c("n_amino_cells", "n_total_cells"), ct_name = "n_amino_cells", tr_name = "n_total_cells", read_min_count = read_min_count, read_min_sample = read_min_sample, collapse=TRUE, collapse_var = "sample")

DDS_ls_Osi <- GetDDS(data = data_Osi, trt_var = "condition", trt_ref_level = "DOX", other_cond_vars = "sample_pair", row_group_var = "Target", select_cols = c("samID", "Target", "Position", "amino_change", "n_amino_cells", "n_total_cells"), sample_id_vars = "samID", value_vars = c("n_amino_cells", "n_total_cells"), ct_name = "n_amino_cells", tr_name = "n_total_cells", read_min_count = read_min_count, read_min_sample = read_min_sample, collapse=TRUE, collapse_var = "sample")

res_Erl_full <- data.frame(cbind(
  rowData(DDS_ls_Erl[['dds']])[,c("Target", "Position", "amino_change")],
  DDS_ls_Erl[['res']],
  DDS_ls_Erl[['resLFC']]$res[, c("log2FoldChange", "svalue")],
  DDS_ls_Erl[['resLFC']]$fit$fsr)) |>
  mutate(baseMean_corrected = baseMean/1.5,
         shrink_log2FoldChange = log2FoldChange.1,
         local_fsr = conditionErl) |>
  arrange(desc(shrink_log2FoldChange)) |>
  select(-c(log2FoldChange, log2FoldChange.1, conditionErl))

res_Osi_full <- data.frame(cbind(
  rowData(DDS_ls_Osi[['dds']])[,c("Target", "Position", "amino_change")],
  DDS_ls_Osi[['res']],
  DDS_ls_Osi[['resLFC']]$res[, c("log2FoldChange", "svalue")],
  DDS_ls_Osi[['resLFC']]$fit$fsr)) |>
  mutate(baseMean_corrected = baseMean/1.5,
         shrink_log2FoldChange = log2FoldChange.1,
         local_fsr = conditionOsi) |>
  arrange(desc(shrink_log2FoldChange)) |>
  select(-c(log2FoldChange, log2FoldChange.1, conditionOsi))

write.csv(res_Erl_full, file = paste0("output/", set_id, "_MissionBio_lollipop_res_amino_Erl", Sys.Date(),".csv"))
write.csv(res_Osi_full, file = paste0("output/", set_id, "_MissionBio_lollipop_res_amino_Osi", Sys.Date(), ".csv"))

