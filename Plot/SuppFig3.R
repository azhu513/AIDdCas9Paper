## ------------------------------------------------
## Script name: SuppFig3.R
## Purpose of script: Supplementary Figure 3
##
## Author: Kristel Dorighi
## ------------------------------------------------

## Load up the packages ---------------------------

library(tidyverse)
library(ggpubr)
library(reshape)
library(readr)
library(dplyr)

## ------------------------------------------------

#Time course of base editor induction with DOX

#load data
D <- read_tsv(file = "AIDdCas9_timecourse.tsv")
S <- read.csv(file = "AIDdCas9_timecourse_SampleInfo.csv", header = TRUE)

SD <- merge(D, S, by.x = "SamId", by.y = "SamID")
SD1 <- SD %>% filter(Target.x == Target.y)
SD2 <- SD1 %>% select(SamId, Allele:Percent, TotalReads, AminoChange, ProteinPosition, DOX, Replicate, Target.x)
SD2 <- SD2 %>% dplyr::rename(sgRNA = Target.x)



#count the number of mutations per allele

# split on |
mutations <- strsplit(
  SD2$Allele,
  "\\|"
)

SD2 <- SD2[rep(1:nrow(SD2), sapply(mutations, length)), ]
SD2$mutation <- unlist(mutations)

# add location info 
SD2$location <- as.numeric(gsub("\\D", "", SD2$mutation))
SD2$substitution <- gsub("[[:digit:]]", "", SD2$mutation)
SD2$ref_nt <- sapply(SD2$substitution, function(x) {
  strsplit(x, ">")[[1]][1]
})
SD2$mut_nt <- sapply(SD2$substitution, function(x) {
  strsplit(x, ">")[[1]][2]
})

# look at substitutions and WT alleles only:
SD2s <- SD2[grepl(">", SD2$mutation), ]
SD2w <- SD2[grepl("WT", SD2$Allele), ]

SD2s <- rbind(SD2s, SD2w)

#Get info on C>T vs 'other'

SD2s$MutClass <- "other"
SD2s$MutClass[SD2s$substitution == "C>T"] <- "C>T"
SD2s$MutClass[SD2s$substitution == "C>A"] <- "C>A"
SD2s$MutClass[SD2s$substitution == "C>G"] <- "C>G"
SD2s$MutClass[SD2s$substitution == "G>A"] <- "G>A"
SD2s$MutClass[SD2s$substitution == "G>C"] <- "G>C"
SD2s$MutClass[SD2s$substitution == "G>T"] <- "G>T"
SD2s$MutClass[SD2s$substitution == "WT"] <- "WT"




#calculate mutation frequency at each location for each type of mutation
SA <- SD2s %>% 
  select(sgRNA, DOX, Replicate, AlleleCount, TotalReads, mutation:MutClass) %>% 
  group_by(sgRNA, DOX, Replicate, location, MutClass, TotalReads) %>% 
  summarize(MutCount = sum(AlleleCount)) %>% 
  transmute(MutFreq = MutCount/TotalReads)


#add in location information
TF <- read.delim(file = "sgRNA_PAMs.txt", header = TRUE, sep = "\t")
TFs <- TF %>% select(Name, Position, Strand, Sequence, PAM)
TFSA <- merge(SA, TFs, by.x = "sgRNA", by.y = "Name")
TFSA$PAMrel <- TFSA$location - TFSA$Position


#change DOX names from D1DOX to D1 DOX and use only sgRNA_19
TFSA1 <- TFSA

TFSA1["DOX"][TFSA1["DOX"]== "noDOX"] <- "no DOX"
TFSA1["DOX"][TFSA1["DOX"]== "D1DOX"] <- "D1 DOX"
TFSA1["DOX"][TFSA1["DOX"]== "D3DOX"] <- "D3 DOX"
TFSA1["DOX"][TFSA1["DOX"]== "D6DOX"] <- "D6 DOX"
TFSA1["DOX"][TFSA1["DOX"]== "D9DOX"] <- "D9 DOX"


#Plots (update DOX value to get each individual plot)

#CD81-26
TFSA2 <- TFSA1 %>% 
  filter(PAMrel >= -30 & PAMrel <= -9) %>%
  filter(MutClass != "other") %>%
  filter(sgRNA == "sgCD81ex5_26", DOX == "no DOX")
fig <- ggbarplot(
  TFSA2, x = "PAMrel", y = "MutFreq", add = "mean_se", error.plot = "errorbar", fill = "MutClass"
) +
  scale_fill_manual("MutClass", values = c("C>A"="#94D2BD","C>G"="#0A9396","C>T"="#005F73", "G>A"="#BB3E03", "G>C"="#EE9B00", "G>T"="#E9D8A6"))+
  labs(x = "Distance from PAM", y = "Mutation Frequency") + ylim(0,0.7)
fig + theme_classic() + rotate_x_text(90)


#CD81-19
TFSA2 <- TFSA1 %>% 
  filter(PAMrel >= -30 & PAMrel <= -9) %>%
  filter(MutClass != "other") %>%
  filter(sgRNA == "sgCD81ex5_19", DOX == "no DOX")
fig <- ggbarplot(
  TFSA2, x = "PAMrel", y = "MutFreq", add = "mean_se", error.plot = "errorbar", fill = "MutClass"
) +
  scale_fill_manual("MutClass", values = c("C>A"="#94D2BD","C>G"="#0A9396","C>T"="#005F73", "G>A"="#BB3E03", "G>C"="#EE9B00", "G>T"="#E9D8A6"))+
  labs(x = "Distance from PAM", y = "Mutation Frequency") + ylim(0,0.7)
fig + theme_classic() + rotate_x_text(90)





#calculate the number of mutations per allele
SA <- SD2s %>% 
  select(Allele, sgRNA, DOX, Replicate, AlleleCount, TotalReads) %>%
  group_by(Allele, sgRNA, DOX, Replicate, AlleleCount, TotalReads) %>% 
  mutate(MutNumber = n()) %>%
  unique() %>%
  ungroup()


#assign WT alleles 0 and condense alleles with more than 6 mutations
SA$MutNumber[SA$Allele == "WT"] <- "0"
SA$MutNumber[SA$MutNumber >= "6"] <- "6+"

#determine a mutation frequency for each allele in each condition
SA1 <- SA %>% 
  group_by(Allele, MutNumber, Replicate, sgRNA, DOX, TotalReads) %>% 
  summarize(MutCount = sum(AlleleCount)) %>% 
  transmute(MutFreq = MutCount/TotalReads) %>%
  ungroup()

#Add up the mutation frequencies of alleles with different numbers of mutations
SA2 <- SA1 %>% group_by(MutNumber, sgRNA, DOX, Replicate) %>%
  summarize(TotalFreq = sum(MutFreq))

#change DOX names from D1DOX to D1 DOX
SA4 <- SA2

SA4["DOX"][SA4["DOX"]== "noDOX"] <- "no DOX"
SA4["DOX"][SA4["DOX"]== "D1DOX"] <- "D1 DOX"
SA4["DOX"][SA4["DOX"]== "D3DOX"] <- "D3 DOX"
SA4["DOX"][SA4["DOX"]== "D6DOX"] <- "D6 DOX"
SA4["DOX"][SA4["DOX"]== "D9DOX"] <- "D9 DOX"


#add error bars
SA7 <- SA4

SA7$DOX <- factor(SA7$DOX, levels=c("no DOX","D1 DOX", "D3 DOX","D6 DOX", "D9 DOX"))
SA7$MutNumber <- factor(SA7$MutNumber, levels=c("0","1", "2","3", "4", "5", "6+"))

#excluding WT
SA8 <- SA7 %>% filter(MutNumber != "0")

#black & white
fig <- ggbarplot(SA8, y = "TotalFreq", facet.by = c("sgRNA", "DOX"), fill = "MutNumber", add = "mean_se", error.plot = "errorbar", position = position_stack(), 
                 palette = c("#F8F9FA","#DEE2E6","#ADB5BD", "#6C757D", "#495057", "#212529"))
fig + theme_classic()



#number of alleles per dox time point
SDtop <- SD2 %>% filter(Percent > 0.001 & Allele != "WT") 
SDcount <- SDtop %>% group_by(sgRNA, DOX, Replicate) %>% summarise(N = n())
write.csv(SDcount, file = "DOXTimecourse_alleleCount.csv")


