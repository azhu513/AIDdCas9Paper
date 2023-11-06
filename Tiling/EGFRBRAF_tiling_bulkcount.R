## ---------------------------------------------------
## Script name: EGFRBRAF_tiling_bulkcount.R
## Purpose of script: Model the bulk count of
##                    amlicon-seq in the EGFR/BRAF
##                    tiling experiment
## Author: Anqi Zhu
## Date Created: 2023-03-12
## Date Modified: 2023-05-22
## ---------------------------------------------------

## Load up the packages ------------------------------

suppressMessages(require(dplyr))
suppressMessages(require(tidyr))
suppressMessages(require(ggplot2))

## ------------------------------------------------

## load up our functions into memory --------------

source("../functions/utils.R")
source("../functions/model_helper")

## ------------------------------------------------

## load previously processed data
load("data/dataset1_clean.RData")

## further split time replicate to time and replicate
EGFRtarget_keep$time <- gsub(pattern = "[^0-9.-]", "", EGFRtarget_keep$time_replicate)
BRAFtarget_keep$time <- gsub(pattern = "[^0-9.-]", "", BRAFtarget_keep$time_replicate)

EGFRtarget_keep$rep <- ifelse(EGFRtarget_keep$time_replicate != "", substr(EGFRtarget_keep$time_replicate, 3,3), "")
BRAFtarget_keep$rep <- ifelse(BRAFtarget_keep$time_replicate != "", substr(BRAFtarget_keep$time_replicate, 3,3), "")

## Check predrug for every time point, and replicates

#### EGFR has predrug and both replicates for all time points
EGFRtarget_keep <- EGFRtarget_keep %>% group_by(time_replicate, treatment)
table(group_keys(EGFRtarget_keep))
EGFRtarget_keep <- EGFRtarget_keep %>% ungroup()

#### BRAF has predrug for all time points and only replicate a for D0 and only replicate b for D3
BRAFtarget_keep <- BRAFtarget_keep %>% group_by(time_replicate, treatment)
table(group_keys(BRAFtarget_keep))
BRAFtarget_keep <- BRAFtarget_keep %>% ungroup()

## a pairwise comparison between time points for pre-drug samples by

## EGFR
plt_conditions <- EGFRtarget_keep %>%
  select(SamId, SampleName, time_replicate, time, rep, treatment) %>%
  distinct(.keep_all = TRUE)

EGFRtarget_keep  %>%
  filter(treatment=="PreDrug") %>%
  mutate(log_count = log(AlleleCount),
         target_allele = paste0(Target, "_", Allele)) %>%
  filter(log_count < 6) %>%
  select(time_replicate, Target, target_allele, log_count) %>%
  tidyr::pivot_wider(names_from=time_replicate, values_from=log_count) %>%
  tibble::column_to_rownames('target_allele')  %>%
  GGally::ggpairs(aes(color=Target), columns = c("D6a", "D6b", "D9a", "D9b", "D3a", "D3b", "D0a", "D0b"), title = "EGFR log allele count")

## BRAF

plt_conditions <- BRAFtarget_keep %>%
  select(SamId, SampleName, time_replicate, time, rep, treatment) %>%
  distinct(.keep_all = TRUE)

BRAFtarget_keep  %>%
  filter(treatment=="PreDrug") %>%
  mutate(log_count = log(AlleleCount),
         target_allele = paste0(Target, "_", Allele)) %>%
  filter(log_count < 6) %>%
  select(time_replicate, Target, target_allele, log_count) %>%
  tidyr::pivot_wider(names_from=time_replicate, values_from=log_count) %>%
  tibble::column_to_rownames('target_allele')  %>%
  GGally::ggpairs(aes(color=Target), columns = c("D6a", "D6b", "D9a", "D9b", "D3a", "D3b", "D0a", "D0b"), title = "log allele count")

## Form the dataset for DESeq2 for each target region, due to the reason that each samples are paired with target region
# region_dat_list <- large_data %>% group_by(Target) %>% group_split()
#
# Compare the drug vs predrug condition for each drug (and dose)
large_data <- rbind(EGFRtarget_keep, BRAFtarget_keep) |> ungroup()

for (trt in c("Erlot100", "Erlot500", "Osi020", "Osi100")){
  dat <- large_data |>
    filter(treatment == "PreDrug"|treatment == trt) |>
    mutate(sample_id = paste0(treatment, "-", time_replicate))

  dds <- GetDDS(data = dat, trt_var = "treatment", trt_ref_level = "PreDrug", other_cond_vars = "time_replicate", row_group_var = "Target", select_cols = c("sample_id", "Target", "Position", "Allele", "AlleleCount", "TotalReads", "Type", "AminoChange", "Impact", "ProteinPosition"), sample_id_vars = "sample_id", value_vars = c("AlleleCount", "TotalReads"), ct_name = "AlleleCount", tr_name = "TotalReads", read_min_cout = 5, read_min_sample = 3)

  res_full <- data.frame(cbind(
    rowData(dds[['dds']])[,c("Target", "Position", "Allele", "Amino_acids", "protein_position")],
    dds[['res']],
    dds[['resLFC']]$res[, c("log2FoldChange", "svalue")],
    dds[['resLFC']]$fit$fsr)) |>
    mutate(shrink_log2FoldChange = log2FoldChange.1) |>
    arrange(desc(shrink_log2FoldChange)) |>
    select(-c(log2FoldChange, log2FoldChange.1))

  write.csv(res_full, file = paste0("./output/EGFRBRAF_tiling_bulk_", trt, ".csv"), quote = T, row.names = F)
}




