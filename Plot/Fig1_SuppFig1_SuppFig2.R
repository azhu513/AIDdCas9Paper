## ------------------------------------------------
## Script name: Fig1_SuppFig1_SuppFig2
## Purpose of script: Figure 1 & Supplementary 
##                    Figure 1 & 2
##
## Author: Kristel Dorighi
## ------------------------------------------------

## Load up the packages ---------------------------

library(tidyverse)
library(ggpubr)
library(reshape)

## ------------------------------------------------

## load datsets and sample info and join

D <- read_tsv(file = "DUAL_ABE_CBE_tabdelim.tsv")
S <- read.csv(file = "DUAL_ABE_CBE_SampleInfo.csv", header = TRUE)
SD <- left_join(D, S, by = "SamId")

#filter for the base editors that we want to use in this analysis (ie. AIDdCas9, BE4, AIDCas9n2xUGI, AIDCas9n)
#filter so that we're only looking at alignments for the specific sgRNAs used in each sample or parental (ie. none)


SD1 <- SD %>% filter(
    BaseEditor == "AID-Cas9n-2xUGI" | 
    BaseEditor == "AID-Cas9n" |
    BaseEditor == "BE4" |
    BaseEditor == "AID-dCas9"
) %>%
  filter(
    sgRNA == Target | 
      sgRNA == "none" |
      sgRNA == "NTC_39"
  )


#select only the columns we need for the analysis
SD2 <- SD1 %>% select(Allele:Percent, TotalReads, AminoChange, ProteinPosition, SamId, BaseEditor, Target, sgRNA, Replicate, Sequence)

#make a copy to use for subsequent analysis of amino acid change and variant type
SD4 <- SD1 %>% select(Allele:Percent, TotalReads, AminoChange, ProteinPosition, SamId, BaseEditor, Target, sgRNA, Replicate, Type)


#Percent Editing efficiency at each nucleotide

#split alleles into nucleotide change and frequency

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

# look at substitutions only:
SD2s <- SD2[grepl(">", SD2$mutation), ]

#Get info on C>T vs 'other'
SD2s$MutClass <- "other"
SD2s$MutClass[SD2s$substitution == "C>T"] <- "C>T"
SD2s$MutClass[SD2s$substitution == "C>A"] <- "C>A"
SD2s$MutClass[SD2s$substitution == "C>G"] <- "C>G"
SD2s$MutClass[SD2s$substitution == "G>A"] <- "G>A"
SD2s$MutClass[SD2s$substitution == "G>C"] <- "G>C"
SD2s$MutClass[SD2s$substitution == "G>T"] <- "G>T"

#calculate mutation frequency at each location according to the MutClass for each base editor, sgRNA and replicate
SA <- SD2s %>% 
  select(Target, sgRNA, BaseEditor, Replicate, AlleleCount, TotalReads, mutation:MutClass) %>% 
  group_by(Target, sgRNA, BaseEditor, Replicate, location, MutClass, TotalReads) %>% 
  summarize(MutCount = sum(AlleleCount)) %>% 
  transmute(MutFreq = MutCount/TotalReads)

#compute position relative to PAM
TF <- read.delim(file = "sgRNA_PAMs.txt", header = TRUE, sep = "\t")
TFs <- TF %>% select(Name, Position)
TFSA <- merge(SA, TFs, by.x = "Target", by.y = "Name")
TFSA$PAMrel <- TFSA$location - TFSA$Position


#Plots with standard error for error bars and fixed Y-axis to 0.7 (update BaseEditor and sgRNA to make each plot)
SA1 <- TFSA %>% 
  filter(BaseEditor == "AID-dCas9", sgRNA == "sgCD81ex5_25")
SA2 <- SA1 %>%
  filter(PAMrel >= -30 & PAMrel <= 0)

fig <- ggbarplot(
  SA2, x = "PAMrel", y = "MutFreq", add = "mean_se", error.plot = "errorbar", fill = "MutClass"
) +
  scale_fill_manual("MutClass", values = c("C>A"="#94D2BD","C>G"="#0A9396","C>T"="#005F73", "G>A"="#BB3E03", "G>C"="#EE9B00", "G>T"="#E9D8A6", "other"="gray"))+
  labs(x = "Distance from PAM", y = "Mutation Frequency") + ylim(0,0.7)
fig + theme_classic() + rotate_x_text(90)


#both C->X and G->X
#scale - 500x500
SA1 <- TFSA %>% 
  filter(PAMrel >= -30 & PAMrel <= 0) %>%
  filter(BaseEditor == "BE4", sgRNA == "NTC_39")
fig <- ggbarplot(
  SA1, x = "PAMrel", y = "MutFreq", add = "mean_se", error.plot = "errorbar", fill = "MutClass"
) +
  scale_fill_manual("MutClass", values = c("C>A"="#94D2BD","C>G"="#0A9396","C>T"="#005F73", "G>A"="#BB3E03", "G>C"="#EE9B00", "G>T"="#E9D8A6", "other"="gray"))+
  labs(x = "Distance from PAM", y = "Mutation Frequency") + ylim(0,0.7)
fig + theme_classic() + rotate_x_text(90)


#for NTC
TFSA_NTC <- TFSA %>% 
  filter(sgRNA == "NTC_39") %>% group_by(BaseEditor, Replicate, MutClass, location) %>% 
  transmute(NTCPercent = mean(MutFreq)) %>%
  mutate(ampPos = location - 2395403)

range(TFSA_NTC$location)
2395403 2395504

SA1 <- TFSA_NTC %>% 
  #filter(PAMrel >= -30 & PAMrel <= 0) %>%
  filter(BaseEditor == "BE4")
fig <- ggbarplot(
  SA1, x = "ampPos", y = "NTCPercent", add = "mean_se", error.plot = "errorbar", fill = "MutClass"
) +
  scale_fill_manual("MutClass", values = c("C>A"="#94D2BD","C>G"="#0A9396","C>T"="#005F73", "G>A"="#BB3E03", "G>C"="#EE9B00", "G>T"="#E9D8A6", "other"="gray"))+
  labs(x = "Position in Amplicon", y = "Percent Edited") + ylim(0,0.7)
fig + theme_classic() + rotate_x_text(90)


#select Cs - 16
SA1 <- TFSA %>% 
  filter(MutClass == "C>T" | MutClass == "C>A" | MutClass == "C>G", sgRNA == "sgCD81ex5_16", BaseEditor == "AID-dCas9")
SA2 <- SA1 %>%
  filter(PAMrel >= -30 & PAMrel <= -5)
fig <- ggbarplot(
  SA2, x = "PAMrel", y = "MutFreq", add = "mean_se", error.plot = "errorbar", fill = "MutClass", 
) +
  scale_fill_manual("MutClass", values = c("C>A"="#94D2BD","C>G"="#0A9396","C>T"="#005F73"))+
  labs(x = "Distance from PAM", y = "Mutation Frequency") + ylim(0,0.7) 
fig + theme_classic() 

#select Cs - 19
SA1 <- TFSA %>% 
  filter(MutClass == "C>T" | MutClass == "C>A" | MutClass == "C>G", sgRNA == "sgCD81ex5_19", BaseEditor == "AID-dCas9")
SA2 <- SA1 %>%
  filter(PAMrel >= -25 & PAMrel <= -9)
fig <- ggbarplot(
  SA2, x = "PAMrel", y = "MutFreq", add = "mean_se", error.plot = "errorbar", fill = "MutClass", 
) +
  scale_fill_manual("MutClass", values = c("C>A"="#94D2BD","C>G"="#0A9396","C>T"="#005F73"))+
  labs(x = "Distance from PAM", y = "Mutation Frequency") + ylim(0,0.7) 
fig + theme_classic() 

#select Cs - 25
SA1 <- TFSA %>% 
  filter(MutClass == "C>T" | MutClass == "C>A" | MutClass == "C>G", sgRNA == "sgCD81ex5_25", BaseEditor == "BE4")
SA2 <- SA1 %>%
  filter(PAMrel >= -25 & PAMrel <= -7)
fig <- ggbarplot(
  SA2, x = "PAMrel", y = "MutFreq", add = "mean_se", error.plot = "errorbar", fill = "MutClass", 
) +
  scale_fill_manual("MutClass", values = c("C>A"="#94D2BD","C>G"="#0A9396","C>T"="#005F73"))+
  labs(x = "Distance from PAM", y = "Mutation Frequency") + ylim(0,0.7) 
fig + theme_classic() 

#select Cs - 26
SA1 <- TFSA %>% 
  filter(MutClass == "C>T" | MutClass == "C>A" | MutClass == "C>G", sgRNA == "sgCD81ex5_26", BaseEditor == "AID-dCas9")
SA2 <- SA1 %>%
  filter(PAMrel >= -25 & PAMrel <= -5)
fig <- ggbarplot(
  SA2, x = "PAMrel", y = "MutFreq", add = "mean_se", error.plot = "errorbar", fill = "MutClass", 
) +
  scale_fill_manual("MutClass", values = c("C>A"="#94D2BD","C>G"="#0A9396","C>T"="#005F73"))+
  labs(x = "Distance from PAM", y = "Mutation Frequency") + ylim(0,0.7) 
fig + theme_classic() 


#select Gs - 16
SA1 <- TFSA %>% 
  filter(MutClass == "G>T" | MutClass == "G>C" | MutClass == "G>A", sgRNA == "sgCD81ex5_16", BaseEditor == "BE4")
SA2 <- SA1 %>%
  filter(PAMrel >= -35 & PAMrel <= -10)
fig <- ggbarplot(
  SA2, x = "PAMrel", y = "MutFreq", add = "mean_se", error.plot = "errorbar", fill = "MutClass", 
) +
  scale_fill_manual("MutClass", values = c("G>A"="#BB3E03", "G>C"="#EE9B00", "G>T"="#E9D8A6"))+
  labs(x = "Distance from PAM", y = "Mutation Frequency") + ylim(0,0.25) 
fig + theme_classic() 

#select Gs - 19
SA1 <- TFSA %>% 
  filter(MutClass == "G>T" | MutClass == "G>C" | MutClass == "G>A", sgRNA == "sgCD81ex5_19", BaseEditor == "BE4")
SA2 <- SA1 %>%
  filter(PAMrel >= -35 & PAMrel <= -10)
fig <- ggbarplot(
  SA2, x = "PAMrel", y = "MutFreq", add = "mean_se", error.plot = "errorbar", fill = "MutClass", 
) +
  scale_fill_manual("MutClass", values = c("G>A"="#BB3E03", "G>C"="#EE9B00", "G>T"="#E9D8A6"))+
  labs(x = "Distance from PAM", y = "Mutation Frequency") + ylim(0,0.25) 
fig + theme_classic() 

#select Gs - 25
SA1 <- TFSA %>% 
  filter(MutClass == "G>T" | MutClass == "G>C" | MutClass == "G>A", sgRNA == "sgCD81ex5_25", BaseEditor == "BE4")
SA2 <- SA1 %>%
  filter(PAMrel >= -35 & PAMrel <= -10)
fig <- ggbarplot(
  SA2, x = "PAMrel", y = "MutFreq", add = "mean_se", error.plot = "errorbar", fill = "MutClass", 
) +
  scale_fill_manual("MutClass", values = c("G>A"="#BB3E03", "G>C"="#EE9B00", "G>T"="#E9D8A6"))+
  labs(x = "Distance from PAM", y = "Mutation Frequency") + ylim(0,0.25) 
fig + theme_classic() 

#select Gs - 26
SA1 <- TFSA %>% 
  filter(MutClass == "G>T" | MutClass == "G>C" | MutClass == "G>A", sgRNA == "sgCD81ex5_26", BaseEditor == "BE4")
SA2 <- SA1 %>%
  filter(PAMrel >= -35 & PAMrel <= -10)
fig <- ggbarplot(
  SA2, x = "PAMrel", y = "MutFreq", add = "mean_se", error.plot = "errorbar", fill = "MutClass", 
) +
  scale_fill_manual("MutClass", values = c("G>A"="#BB3E03", "G>C"="#EE9B00", "G>T"="#E9D8A6"))+
  labs(x = "Distance from PAM", y = "Mutation Frequency") + ylim(0,0.25) 
fig + theme_classic() 

                       
#Analysis of top alleles

SDseq <- SD2 %>% select(sgRNA, BaseEditor, Allele, Percent, AminoChange, ProteinPosition, Sequence) %>% 
  group_by(sgRNA, BaseEditor, Allele, Sequence) %>%
  transmute(meanPercent = mean(Percent)) %>%
  distinct()

#print out files with all alleles for each sgRNA
C25 <- SDseq %>% filter(sgRNA == "sgCD81ex5_25") %>%
  pivot_wider(
    names_from = BaseEditor,
    values_from = meanPercent
  )
write.csv(C25, file = "DualABECBE_sgCD81ex5_25_withseqs.csv")

C26 <- SDseq %>% filter(sgRNA == "sgCD81ex5_26") %>%
  pivot_wider(
    names_from = BaseEditor,
    values_from = meanPercent
  )
write.csv(C26, file = "DualABECBE_sgCD81ex5_26_withseqs.csv")

C19 <- SDseq %>% filter(sgRNA == "sgCD81ex5_19") %>%
  pivot_wider(
    names_from = BaseEditor,
    values_from = meanPercent
  )
write.csv(C19, file = "DualABECBE_sgCD81ex5_19_withseqs.csv")

C16 <- SDseq %>% filter(sgRNA == "sgCD81ex5_16") %>%
  pivot_wider(
    names_from = BaseEditor,
    values_from = meanPercent
  )
write.csv(C16, file = "DualABECBE_sgCD81ex5_16_withseqs.csv")




#Analysis of Allele Type

#Categorize edits into groups
#Catergories: C>T only, contains C>A/C>G edits but no G edits, G>N edits, 
#mixed alleles are categoriezed by the latest assignment

SD3AC <- SD2

SD3AC$AlleleClass <- "other"
SD3AC$AlleleClass[grepl("WT", SD3AC$Allele)] <- "WT"
SD3AC$AlleleClass[grepl("C>T", SD3AC$Allele)] <- "C>T only"
SD3AC$AlleleClass[grepl("C>G", SD3AC$Allele)] <- "C>G/A"
SD3AC$AlleleClass[grepl("C>A", SD3AC$Allele)] <- "C>G/A"
SD3AC$AlleleClass[grepl("G>C", SD3AC$Allele)] <- "G>N"
SD3AC$AlleleClass[grepl("G>T", SD3AC$Allele)] <- "G>N"
SD3AC$AlleleClass[grepl("G>A", SD3AC$Allele)] <- "G>N"
SD3AC$AlleleClass[grepl("D", SD3AC$Allele) | grepl("I", SD3AC$Allele)] <- "Indel"


table(SD3AC$BaseEditor, SD3AC$AlleleClass)

#determine the percent of total reads in each category
#sum the allele counts to get the total reads for each replicate
#sum within each category to get the frequency of each category

ACounts <- SD3AC %>% select(sgRNA, Target, BaseEditor, AlleleCount, Replicate, AlleleClass) %>% 
  group_by(sgRNA, Target, BaseEditor, Replicate, AlleleClass) %>%
  transmute(sumAlleleCount = sum(AlleleCount)) %>%
  distinct()

TCounts <- ACounts %>% 
  group_by(sgRNA, Target, BaseEditor, Replicate) %>%
  mutate(totalReadCount = sum(sumAlleleCount)) %>%
  distinct()

TCounts$Percent <- TCounts$sumAlleleCount / TCounts$totalReadCount

#Allele classes with CD81 sgRNAs
#Plots with error bars

#sg-19
AClass <- TCounts %>% 
  filter(sgRNA == "sgCD81ex5_19")
AClass$BaseEditor <- factor(AClass$BaseEditor, levels = c("BE4", "AID-Cas9n-2xUGI", "AID-Cas9n", "AID-dCas9"))
AClass$AlleleClass <- factor(AClass$AlleleClass, levels = c("WT", "G>N", "C>G/A", "C>T only", "Indel", "other"))
fig <- ggbarplot(
  AClass, x = "BaseEditor", y = "Percent", add = "mean_se", error.plot = "errorbar", fill = "AlleleClass"
) +
  scale_fill_manual("AlleleClass", values = c("WT" = "white","C>G/A"="#94D2BD","C>T only"="#005F73", "G>N"="#EE9B00", "Indel" = "grey", "other" = "#5F6369"))+
  labs(y = "Percent")
fig + theme_classic() + rotate_x_text(90)

#sg-26
AClass <- TCounts %>% 
  filter(sgRNA == "sgCD81ex5_26")
AClass$BaseEditor <- factor(AClass$BaseEditor, levels = c("BE4", "AID-Cas9n-2xUGI", "AID-Cas9n", "AID-dCas9"))
AClass$AlleleClass <- factor(AClass$AlleleClass, levels = c("WT", "G>N", "C>G/A", "C>T only", "Indel", "other"))
fig <- ggbarplot(
  AClass, x = "BaseEditor", y = "Percent", add = "mean_se", error.plot = "errorbar", fill = "AlleleClass"
) +
  scale_fill_manual("AlleleClass", values = c("WT" = "white","C>G/A"="#94D2BD","C>T only"="#005F73", "G>N"="#EE9B00", "Indel" = "grey", "other" = "#5F6369"))+
  labs(y = "Percent")
fig + theme_classic() + rotate_x_text(90)

#sg-16
AClass <- TCounts %>% 
  filter(sgRNA == "sgCD81ex5_16")
AClass$BaseEditor <- factor(AClass$BaseEditor, levels = c("BE4", "AID-Cas9n-2xUGI", "AID-Cas9n", "AID-dCas9"))
AClass$AlleleClass <- factor(AClass$AlleleClass, levels = c("WT", "G>N", "C>G/A", "C>T only", "Indel", "other"))
fig <- ggbarplot(
  AClass, x = "BaseEditor", y = "Percent", add = "mean_se", error.plot = "errorbar", fill = "AlleleClass"
) +
  scale_fill_manual("AlleleClass", values = c("WT" = "white","C>G/A"="#94D2BD","C>T only"="#005F73", "G>N"="#EE9B00", "Indel" = "grey", "other" = "#5F6369"))+
  labs(y = "Percent")
fig + theme_classic() + rotate_x_text(90)

#sg-25
AClass <- TCounts %>% 
  filter(sgRNA == "sgCD81ex5_25")
AClass$BaseEditor <- factor(AClass$BaseEditor, levels = c("BE4", "AID-Cas9n-2xUGI", "AID-Cas9n", "AID-dCas9"))
AClass$AlleleClass <- factor(AClass$AlleleClass, levels = c("WT", "G>N", "C>G/A", "C>T only", "Indel", "other"))
fig <- ggbarplot(
  AClass, x = "BaseEditor", y = "Percent", add = "mean_se", error.plot = "errorbar", fill = "AlleleClass"
) +
  scale_fill_manual("AlleleClass", values = c("WT" = "white","C>G/A"="#94D2BD","C>T only"="#005F73", "G>N"="#EE9B00", "Indel" = "grey", "other" = "#5F6369"))+
  labs(y = "Percent")
fig + theme_classic() + rotate_x_text(90)


#Allele classes with the non-targeting control (NTC-39)

#take the mean editing efficiency across all four targets per replicate
AClass <- TCounts %>% 
  filter(sgRNA == "NTC_39") %>% group_by(BaseEditor, Replicate, AlleleClass) %>% 
  summarize(NTCPercent = mean(Percent))
AClass$BaseEditor <- factor(AClass$BaseEditor, levels = c("BE4", "AID-Cas9n-2xUGI", "AID-Cas9n", "AID-dCas9"))
AClass$AlleleClass <- factor(AClass$AlleleClass, levels = c("WT", "G>N", "C>G/A", "C>T only", "Indel", "other"))
fig <- ggbarplot(
  AClass, x = "BaseEditor", y = "NTCPercent", add = "mean_se", error.plot = "errorbar", fill = "AlleleClass"
) +
  scale_fill_manual("AlleleClass", values = c("WT" = "white","C>G/A"="#94D2BD","C>T only"="#005F73", "G>N"="#EE9B00", "Indel" = "grey", "other" = "#5F6369"))+
  labs(y = "Percent")
fig + theme_classic() + rotate_x_text(90)




#Analysis of amino acid changes
#for all alleles greater than 0.1% frequency


SDtop <- SD4 %>% filter(Percent > 0.001 & Allele != "WT" & sgRNA != "none") 

SDaa <- SDtop %>% 
  group_by(BaseEditor, sgRNA, Target, Replicate, AminoChange, ProteinPosition, TotalReads, Type) %>% 
  transmute(AACount = sum(AlleleCount)) %>%
  mutate(AAFreq = AACount/TotalReads) %>%
  unique()

SDaa$TypeS <- "none"
SDaa$TypeS[SDaa$Type == "inframe_deletion" ] <- "indel (inframe)"
SDaa$TypeS[SDaa$Type == "LD" ] <- "indel (frameshift)"
SDaa$TypeS[SDaa$Type == "inframe_insertion"] <- "indel (inframe)"
SDaa$TypeS[SDaa$Type == "frameshift_variant"] <- "indel (frameshift)"
SDaa$TypeS[SDaa$Type == "protein_altering_variant"] <- "indel (inframe)"
SDaa$TypeS[SDaa$Type == "synonymous_variant"] <- "synonymous variant"
SDaa$TypeS[SDaa$Type == "missense_variant"] <- "missense variant"
SDaa$TypeS[SDaa$Type == "stop_gained"] <- "stop gained"
SDaa$TypeS[SDaa$Type == "coding_sequence_variant,intron_variant"] <- "missense variant"
SDaa$TypeS[SDaa$Type == "splice_acceptor_variant"] <- "missense variant"

AAChangeCount <- SDaa %>% filter(sgRNA != "NTC_39") %>% group_by(BaseEditor, sgRNA, Replicate, TypeS) %>% summarise(N = n())

#Plots

AAChangeCount$BaseEditor <- factor(AAChangeCount$BaseEditor, levels = c("BE4", "AID-Cas9n-2xUGI", "AID-Cas9n", "AID-dCas9"))
AAChangeCount$TypeS <- factor(AAChangeCount$TypeS, levels = c("missense variant", "stop gained", "synonymous variant","indel (inframe)", "indel (frameshift)"))

#16
AAChangeCount1 <- AAChangeCount %>% filter(sgRNA == "sgCD81ex5_16")

fig <- ggbarplot(
  AAChangeCount1, x = "BaseEditor", y = "N", add = "mean_se", error.plot = "errorbar", fill = "TypeS"
) +
  scale_fill_manual("TypeS", values = c("missense variant"="#0A9396","stop gained"="#E9D8A6", "synonymous variant"="white", "indel (inframe)" = "grey", "indel (frameshift)" = "#5F6369"))+
  labs(y = "total number")
fig + theme_classic() + rotate_x_text(90) + ylim(0,65)

#19
AAChangeCount1 <- AAChangeCount %>% filter(sgRNA == "sgCD81ex5_19")

fig <- ggbarplot(
  AAChangeCount1, x = "BaseEditor", y = "N", add = "mean_se", error.plot = "errorbar", fill = "TypeS"
) +
  scale_fill_manual("TypeS", values = c("missense variant"="#0A9396","stop gained"="#E9D8A6", "synonymous variant"="white", "indel (inframe)" = "grey", "indel (frameshift)" = "#5F6369"))+
  labs(y = "total number")
fig + theme_classic() + rotate_x_text(90) + ylim(0,65)

#25
AAChangeCount1 <- AAChangeCount %>% filter(sgRNA == "sgCD81ex5_25")

fig <- ggbarplot(
  AAChangeCount1, x = "BaseEditor", y = "N", add = "mean_se", error.plot = "errorbar", fill = "TypeS"
) +
  scale_fill_manual("TypeS", values = c("missense variant"="#0A9396","stop gained"="#E9D8A6", "synonymous variant"="white", "indel (inframe)" = "grey", "indel (frameshift)" = "#5F6369"))+
  labs(y = "total number")
fig + theme_classic() + rotate_x_text(90) + ylim(0,65)

#26
AAChangeCount1 <- AAChangeCount %>% filter(sgRNA == "sgCD81ex5_26")

fig <- ggbarplot(
  AAChangeCount1, x = "BaseEditor", y = "N", add = "mean_se", error.plot = "errorbar", fill = "TypeS"
) +
  scale_fill_manual("TypeS", values = c("missense variant"="#0A9396","stop gained"="#E9D8A6", "synonymous variant"="white", "indel (inframe)" = "grey", "indel (frameshift)" = "#5F6369"))+
  labs(y = "number of unique peptide sequences")
fig + theme_classic() + rotate_x_text(90)+ ylim(0,65)

