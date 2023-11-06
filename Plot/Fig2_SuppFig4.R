## ------------------------------------------------
## Script name: Fig2_SuppFig4
## Purpose of script: Figure 2 and Supplementary 
##                    Figure 4
##
## Author: Kristel Dorighi
## ------------------------------------------------

## Load up the packages ---------------------------
library(tidyverse)
library(ggpubr)
library(reshape)
## ------------------------------------------------

#CD81 and OR5M9 sgRNA tiling heatmap

#load data
T <- read.csv(file = "AIDdCas9_OR5M9_CD81.csv", header = TRUE) %>%
  select(SamId, Target, Allele:Percent, TotalReads)
K <- read.csv(file = "AIDdCas9_OR5M9_CD81_SampleInfo.csv", header = TRUE)
TK <- left_join(T, K, by = "SamId")

#filter so that we're only looking at alignments for the specific sgRNAs used in each sample or parental
TK1 <- TK %>%
  filter(
    sgRNA == Target | 
      sgRNA == "PC9 parental" |
      sgRNA == "NTC_39" |
      sgRNA == "NTC_40"
  )

SD2 <- TK1

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


#calculate mutation frequency at each location according to the MutClass for each base editor, sgRNA and take the mean of the replicates
SA <- SD2s %>% 
  select(Target, sgRNA, Replicate, AlleleCount, TotalReads, mutation:MutClass) %>% 
  group_by(Target, sgRNA, Replicate, location, MutClass, TotalReads) %>% 
  summarize(MutCount = sum(AlleleCount)) %>% 
  transmute(MutFreq = MutCount/TotalReads) %>%
  group_by(Target, sgRNA, MutClass, location) %>%
  summarise(meanMutFreq = mean(MutFreq)) %>%
  ungroup()


#select for C>X or G>X only and separate into CD81 and OR5M9
SA1 <- SA %>% filter(Target == sgRNA,
                     MutClass != "other",
                     grepl("CD81", sgRNA)) %>% select(sgRNA, location, MutClass, meanMutFreq)

SO1 <- SA %>% filter(Target == sgRNA,
                     MutClass != "other",
                     grepl("OR5", sgRNA)) %>% select(sgRNA, location, MutClass, meanMutFreq)


#for CD81
#update location to begin at 1
SA1$Position <- SA1$location - 2395401

#sum all the C>X and G>X
#add a new column C>X or G>X
SA1$MutSim <- SA1$MutClass
SA1$MutSim[grepl("C>", SA1$MutClass)] <- "C>X"
SA1$MutSim[grepl("G>", SA1$MutClass)] <- "G>X"

SA2 <- SA1 %>% 
  select(sgRNA, location, meanMutFreq, Position, MutSim) %>% 
  group_by(sgRNA, location, Position, MutSim) %>% 
  summarize(MutSum = sum(meanMutFreq)) %>%
  ungroup()


#convert G>X rows to negative also fill NA as 0
SAG <- SA2 %>% filter(MutSim == "G>X")
SAC <- SA2 %>% filter(MutSim == "C>X")

SAG1 <- SAG 
SAG1$MutSum <- SAG$MutSum * -1

SA3 <- rbind(SAC, SAG1)
SA4 <- pivot_wider(SA3, 
                   names_from = sgRNA,
                   values_from = MutSum,
                   values_fill = 0) %>% select(Position, contains("CD81"))

SA5 <- as.data.frame(t(SA4[,-1]))
colnames(SA5) <- SA4$Position
SA5$sgRNA <- rownames(SA5)
SA6 <- SA5 %>% pivot_longer(cols = !sgRNA, names_to = "Position", values_to = "MutSum" )
SA6$Position <- as.numeric(SA6$Position)

fig <- ggplot(SA6, aes(Position, sgRNA, fill = MutSum, )) +
  geom_tile() + scale_fill_gradient2(low = "#BB3E03", mid = "white", high = "#005F73", na.value = 'white') 

fig + theme_classic()


#for OR5M9

#update location to begin at 1
SO1$Position <- SO1$location - 56462581

#sum all the C>X and G>X
#add a new column C>X or G>X
SO1$MutSim <- SO1$MutClass
SO1$MutSim[grepl("C>", SO1$MutClass)] <- "C>X"
SO1$MutSim[grepl("G>", SO1$MutClass)] <- "G>X"

SO2 <- SO1 %>% 
  select(sgRNA, location, meanMutFreq, Position, MutSim) %>% 
  group_by(sgRNA, location, Position, MutSim) %>% 
  summarize(MutSum = sum(meanMutFreq)) %>%
  ungroup()

#convert G>X rows to negative also fill NA as 0
SOG <- SO2 %>% filter(MutSim == "G>X")
SOC <- SO2 %>% filter(MutSim == "C>X")

SOG1 <- SOG 
SOG1$MutSum <- SOG$MutSum * -1

SO3 <- rbind(SOC, SOG1)

SO4 <- pivot_wider(SO3, 
                   names_from = sgRNA,
                   values_from = MutSum,
                   values_fill = 0) %>% select(Position, contains("OR5"))

SO5 <- as.data.frame(t(SO4[,-1]))
colnames(SO5) <- SO4$Position

SO5$sgRNA <- rownames(SO5)
SO6 <- SO5 %>% pivot_longer(cols = !sgRNA, names_to = "Position", values_to = "MutSum" )
SO6$Position <- as.numeric(SO6$Position)

fig <- ggplot(SO6, aes(Position, sgRNA, fill = MutSum, )) +
  geom_tile() + scale_fill_gradient2(low = "#BB3E03", mid = "white", high = "#005F73", na.value = 'white') 

fig + theme_classic()




#Coverage tiling
#separate into CD81 and OR5M9, take only C and G mutations and restrict from PAM of first sgRNA to PAM of last sgRNA

CovCD81 <- SA %>% filter(Target == sgRNA, MutClass != "other", location > 2395411, location < 2395498,
                         grepl("CD81", sgRNA)) %>% select(sgRNA, location, MutClass, meanMutFreq)

CovOR5M9 <- SA %>% filter(Target == sgRNA, MutClass != "other", location > 56462587, location < 56462764,
                          grepl("OR5", sgRNA)) %>% select(sgRNA, location, MutClass, meanMutFreq)


#CD81
CovCD81$Position <- CovCD81$location - 2395411

#all edited nucleotides
length(unique(CovCD81$Position))
42


#OR5M9
CovOR5M9$Position <- CovOR5M9$location - 56462587

#all
length(unique(CovOR5M9$Position))
74
