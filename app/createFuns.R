library(shiny)
library(data.table)
library(tidyverse)
library(stringr)
library(stringdist)
library(wesanderson)
setwd('~/neocovid-app/app')

colors <- c('#7ed6df','#535c68','#badc58','#30336b',
            '#f6e58d','#6ab04c','#22a6b3','#f9ca24',
            '#4834d4','#95afc0','#e056fd','#ff7979',
            '#be2edd','#c7ecee','#f0932b','#eb4d4b',
            '#130f40','#ffbe76','#dff9fb','#686de0')
##### SARS COV 2 PROTEIN OUTLINE #####
proteinSequences <- fread('data/UP000464024.tsv')
proteinSequences2 <- proteinSequences %>% mutate(prot_start = 1) %>% arrange(GenomePos) %>%
  mutate(polar_end = cumsum(Length) + 100*GenomePos,
         polar_start = polar_end - Length + 1,
         Gene = factor(Gene, levels = Gene),
         Protein = factor(Protein, levels = Protein))
proteinSequences3 <- proteinSequences2 %>% mutate(Chain = strsplit(Chain, split = " CHAIN")) %>%
  unnest(Chain) %>% mutate(Chain = gsub('\\b[CHAIN]* ', '', Chain)) %>%
  mutate(pos_start = as.numeric(gsub("\\b(\\d+)\\.\\.(\\d+); .+", "\\1", Chain)),
         pos_end = as.numeric(gsub("\\b(\\d+)\\.\\.(\\d+); .+", "\\2", Chain))) %>%
  group_by(Gene, pos_start, pos_end) %>%
  mutate(Chain = gsub("\\b\\d+\\.\\.\\d+;.+note=(.+)", "\\1", Chain),
         Chain = strsplit(Chain, split = ";")[[1]][1],
         Chain = gsub("\"", "", Chain),
         Chain = factor(Chain, levels = Chain),
         chain_polar_start = pos_start + polar_start - 1,
         chain_polar_end = pos_end + polar_start - 1)
#
ggplot(data = NULL) +
  geom_hline(yintercept = -1, colour = 'grey80', size = 2) +
  # Proteins
  geom_rect(data = proteinSequences3, aes(xmin = polar_start, xmax = polar_end, ymin = -1.5, ymax = -0.5, fill = Gene)) +
  scale_fill_manual(values = sample(colors, 12)) +
  # Chains
  geom_rect(data = filter(proteinSequences3, polar_start != chain_polar_start | polar_end != chain_polar_end),
            aes(xmin = chain_polar_start, xmax = chain_polar_end, ymin = -1.75, ymax = -0.25), fill = 'grey80', alpha = 0.3) +
  #
  scale_x_continuous(limits = c(0, max(proteinSequences3$polar_end))) + scale_y_continuous(limits = c(-10,3)) +
  coord_polar() +
  theme_void()

predictions <- fread('~/projects/covid/results/tables/SupplementaryTable3_FilteredAntigenBindingPredictions.tsv')
# predictions <- fread('~/projects/covid/classii/results/Table_FilteredClassIIAntigenBindingPredictions.tsv')
epitopes <- unique(predictions$Epitope.Seq)
epiPositions <- data.table(Peptide = epitopes) %>% group_by(Peptide) %>%
  mutate(findPep = ifelse(length(grep(Peptide, proteinSequences2$Sequence))>0, grep(Peptide, proteinSequences2$Sequence), NA)) %>% 
  # filter(!is.na(findPep)) %>%
  mutate(Gene = proteinSequences2$Gene[findPep],
         Position = str_locate(proteinSequences2$Sequence[findPep], Peptide)[1],
         Position_end = str_locate(proteinSequences2$Sequence[findPep], Peptide)[2])
foundEpis <- epiPositions %>% filter(!is.na(findPep)) %>% mutate(findPep = NULL)
notFoundEpis <- epiPositions %>% filter(is.na(findPep))
# #
# m =3
# getMore <- epiPositions %>% filter(is.na(findPep)) %>%
#   mutate(findPep = ifelse(length(agrep(Peptide, proteinSequences2$Sequence, max = m))>0, 
#                           agrep(Peptide, proteinSequences2$Sequence, max =m), NA)) %>%
#   filter(!is.na(findPep))
getMore <- lapply(notFoundEpis$Peptide, function(x){
  all_dist <- lapply(foundEpis$Peptide, function(y){
    adist(x, y)
  }) %>% unlist
  if(length(which(all_dist<2))>0) {
    d <- data.frame(Peptide = x, closestPeptide = foundEpis$Peptide[which(all_dist<2)])
    e <- left_join(d, foundEpis, by = c('closestPeptide' = 'Peptide'))
    return(e)
  } else {
    d <- data.frame(Peptide = x, closestPeptide = as.character(NA))
    return(d)
  }
})
moreEpis <- Reduce(full_join, getMore) %>% filter(!is.na(closestPeptide)) %>% 
  full_join(foundEpis)
prepEpis <- moreEpis %>% left_join(dplyr::select(proteinSequences2, Gene, Protein, polar_start, polar_end)) %>% 
  group_by(Peptide) %>% 
  mutate(peptide.polar_start = polar_start + Position, peptide.polar_end = polar_start + Position_end)

plotEpis <- prepEpis %>% arrange(peptide.polar_start) %>%
  mutate(all_pos = list(seq(peptide.polar_start, peptide.polar_end, by = 1))) %>% 
  unnest(cols = c(all_pos)) %>%
  group_by(all_pos) %>% mutate(n = 1:n()) %>% 
  group_by(Peptide, closestPeptide, Gene, Protein, Position, Position_end, polar_start, polar_end, peptide.polar_start, peptide.polar_end) %>% 
  summarise(y = max(n)) %>% ungroup %>% arrange(peptide.polar_start) %>%
  mutate(y_new = ifelse(y == lag(y, default = 0), y-1, y))

# ggplot(data = NULL) +
#   geom_hline(yintercept = -1, colour = 'grey80', size = 2) +
#   # Proteins
#   geom_rect(data = proteinSequences3, aes(xmin = polar_start, xmax = polar_end, ymin = -1.5, ymax = -0.5, fill = Gene)) +
#   scale_fill_manual(values = sample(colors, 12)) +
#   # Chains
#   geom_rect(data = filter(proteinSequences3, polar_start != chain_polar_start | polar_end != chain_polar_end),
#             aes(xmin = chain_polar_start, xmax = chain_polar_end, ymin = -1.75, ymax = -0.25), fill = 'grey80', alpha = 0.3) +
#   # Peptides
#   geom_rect(data = plotEpis, aes(xmin = peptide.polar_start, xmax = peptide.polar_end, ymin = y_new, ymax = y_new + 0.5, fill = Gene)) +
#   #
#   scale_x_continuous(limits = c(0, max(proteinSequences3$polar_end))) + 
#   scale_y_continuous(limits = c(-10,max(plotEpis$y_new)+2)) +
#   coord_polar() +
#   theme_void()
summPredictions <- predictions %>% 
  group_by(Epitope.Seq, HLA.Allele) %>%
  summarise(`Median Affinity (nM)` = median(Median.Score, na.rm = T),
            `Best Score (nM)` = median(Best.Score, na.rm = T),
            `Best Score Method` = unique(Best.Score.Method, na.rm = T),
            `Predicted Stability (NetMHCstabpan)` = median(Predicted.Stability, na.rm = T),
            `Half Life` = median(Half.Life, na.rm = T),
            `Stability Rank` = median(Stability.Rank, na.rm = T))
summPredictions <- dplyr::select(plotEpis, Peptide, Gene, Protein, Position) %>%
  full_join(summPredictions, by = c("Peptide" = "Epitope.Seq"))
# write.table(summPredictions, file = "~/projects/covid/results_si/int/summariseClassIpredictions.tsv", sep = "\t", quote=F, row.names = F)

# CLASS II 

predictionsII <- fread('~/projects/covid/classii/results/Table_FilteredClassIIAntigenBindingPredictions.tsv')
epitopesII <- unique(predictionsII$`Epitope Seq`) # 2,547
epiPositionsII <- data.table(Peptide = epitopesII) %>% group_by(Peptide) %>%
  mutate(findPep = ifelse(length(grep(Peptide, proteinSequences2$Sequence))>0, 
                          grep(Peptide, proteinSequences2$Sequence), NA)) %>% 
  # filter(!is.na(findPep)) %>%
  mutate(Gene = proteinSequences2$Gene[findPep],
         Position = str_locate(proteinSequences2$Sequence[findPep], Peptide)[1],
         Position_end = str_locate(proteinSequences2$Sequence[findPep], Peptide)[2])
foundEpisII <- epiPositionsII %>% filter(!is.na(findPep)) %>% mutate(findPep = NULL) # 2,348
notFoundEpisII <- epiPositionsII %>% filter(is.na(findPep)) # 199
# #
# m =3
# getMore <- epiPositions %>% filter(is.na(findPep)) %>%
#   mutate(findPep = ifelse(length(agrep(Peptide, proteinSequences2$Sequence, max = m))>0, 
#                           agrep(Peptide, proteinSequences2$Sequence, max =m), NA)) %>%
#   filter(!is.na(findPep))
getMoreII <- lapply(notFoundEpisII$Peptide, function(x){
  all_dist <- lapply(foundEpisII$Peptide, function(y){
    adist(x, y)
  }) %>% unlist
  if(length(which(all_dist<2))>0) {
    d <- data.frame(Peptide = x, closestPeptide = foundEpis$Peptide[which(all_dist<2)])
    e <- left_join(d, foundEpis, by = c('closestPeptide' = 'Peptide'))
    return(e)
  } else {
    d <- data.frame(Peptide = x, closestPeptide = as.character(NA))
    return(d)
  }
})
moreEpisII <- Reduce(full_join, getMoreII) %>% filter(!is.na(closestPeptide)) %>% 
  full_join(foundEpisII) #2,440
prepEpisII <- moreEpisII %>% left_join(dplyr::select(proteinSequences2, Gene, Protein, polar_start, polar_end)) %>% 
  group_by(Peptide) %>% 
  mutate(peptide.polar_start = polar_start + Position, peptide.polar_end = polar_start + Position_end)

plotEpisII <- prepEpisII %>% arrange(peptide.polar_start) %>%
  mutate(all_pos = list(seq(peptide.polar_start, peptide.polar_end, by = 1))) %>% 
  unnest(cols = c(all_pos)) %>%
  group_by(all_pos) %>% mutate(n = 1:n()) %>% 
  group_by(Peptide, closestPeptide, Gene, Protein, Position, Position_end, polar_start, polar_end, peptide.polar_start, peptide.polar_end) %>% 
  summarise(y = max(n)) %>% ungroup %>% arrange(peptide.polar_start) %>%
  mutate(y_new = ifelse(y == lag(y, default = 0), y-1, y))

summPredictionsII <- predictionsII %>% 
  group_by(`Epitope Seq`, `HLA Allele`) %>%
  summarise(`Median Affinity (nM)` = median(`Median Score`, na.rm = T),
            `Best Score (nM)` = median(`Best Score`, na.rm = T),
            `Best Score Method` = unique(`Best Score Method`, na.rm = T)[1]
            # `Predicted Stability (NetMHCstabpan)` = median(`Predicted Stability`, na.rm = T),
            # `Half Life` = median(`Half Life`, na.rm = T),
            # `Stability Rank` = median(`Stability Rank`, na.rm = T)
            )
summPredictionsII <- dplyr::select(plotEpisII , Peptide, Gene, Protein, Position) %>%
  full_join(summPredictionsII, by = c("Peptide" = "Epitope Seq"))
# write.table(summPredictionsII, file = "~/projects/covid/results_si/int/summariseClassIIpredictions.tsv", sep = "\t", quote=F, row.names = F)

##### HLA View #####
summPredictions2 <- dplyr::select(plotEpis, Peptide, Gene, Protein, Position, peptide.polar_start, peptide.polar_end) %>%
  full_join(summPredictions)
summPredictions2II <- dplyr::select(plotEpisII, Peptide, Gene, Protein, Position, peptide.polar_start, peptide.polar_end) %>%
  full_join(summPredictionsII)
pickHLA = "HLA-A*02:02"
getHlaFamily <- gsub("(HLA-[A-C]\\*\\d+):\\d+", "\\1", pickHLA)
filtHlaFamily <- summPredictions2 %>% ungroup %>% mutate(HlaFamily = gsub("(HLA-[A-C]\\*\\d+):\\d+", "\\1", HLA.Allele)) %>%
  filter(HlaFamily == getHlaFamily & HLA.Allele != pickHLA)
summHlaFamily <- filtHlaFamily %>% group_by(Peptide, peptide.polar_start) %>% summarise(n = n())
filterPreds <- filter(summPredictions2, HLA.Allele == pickHLA)
filterPredsPlusHlaFamily <- filter(summPredictions2, HLA.Allele == pickHLA | HLA.Allele %in% filtHlaFamily$HLA.Allele)

ggplot(data = NULL) +
  geom_hline(data = NULL, yintercept = -25, colour = 'grey80', size = 2) +
  # Peptides from the same HLA family members
  geom_segment(data = summHlaFamily, aes(x = peptide.polar_start, xend = peptide.polar_start,
                                         y = 510, yend = Inf, colour = n)) +
  scale_colour_gradientn(colours = wes_palette('Zissou1', 50, type = 'continuous')) +
  # Chains
  geom_rect(data = proteinSequences3, aes(fill = Gene, xmin = chain_polar_start, xmax = chain_polar_end,
                                          ymin = -50, ymax = 0),
            # ymin = -1.75, ymax = -0.25),
            fill = 'grey80', alpha = 0.3) +
  # Proteins
  geom_rect(data = proteinSequences3, aes(fill = Gene, xmin = polar_start, xmax = polar_end,
                                          ymin = -40, ymax = -10)) +
                                          # ymin = -1.5, ymax = -0.5)) +
  scale_fill_manual(values = sample(colors, 12)) +
  # Peptides from specific HLA type
  geom_col(data = filterPreds, aes(fill = Gene, x = peptide.polar_start, 
                                    y = 500-`Median Affinity (nM)`), width = 5) +
  #
  scale_x_continuous(limits = c(min(proteinSequences3$polar_start), max(proteinSequences3$polar_end))) +
  scale_y_continuous(limits = c(-50,550), breaks = c(0, 100, 200, 300, 400, 500, 550),
                     labels = c(500,400,300,200,100,0, paste0('Other ',getHlaFamily," Alleles"))) +
  labs(y = "Predicted Median Affinity (nM)", x = "Position in SARS-CoV-2",
       colour = "N Other Alleles") +
  theme()

#####
# fulldf <- fread('data/tb_for_katie/fulldf.tsv')
freqcountry <- fread('data/tb_for_katie/freq_country_V2.tsv')
perccountry <- fread('data/tb_for_katie/percent_country.tsv')
fulldf <- fread('data/tb_for_katie/fulldf.txt')
#
# download.file("http://thematicmapping.org/downloads/TM_WORLD_BORDERS_SIMPL-0.3.zip" , destfile="data/world_shape_file.zip")
library(rgdal)
world <- readOGR(dsn = "data/",
                 layer = "TM_WORLD_BORDERS_SIMPL-0.3",
                 verbose = FALSE)
# setdiff(fulldf2$World_Country, world$NAME)
# grep('Samoa', world$NAME, value = TRUE)
fulldf <- fulldf %>% mutate(NAME = World_Country)
summ_pops <- fulldf %>% group_by(NAME) %>% 
  summarise(PICK = sample(1:length(colors), 1, replace = TRUE))
world@data <- world@data %>% left_join(summ_pops) %>%
  mutate(PICK = factor(PICK))

#####
plotEpisI = plotEpis %>% mutate(Class = "I")
summPredictionsI = summPredictions %>% mutate(Class = "I") %>% rename("HLA Allele" = "HLA.Allele")
summPredictions2I = summPredictions2 %>% mutate(Class = "I") %>% rename("HLA Allele" = "HLA.Allele")

plotEpis <- plotEpisII %>% mutate(Class = "II") %>% 
  full_join(plotEpisI) %>%
  mutate(Class = factor(Class, levels = c('I','II')))
summPredictions <- summPredictionsII %>% mutate(Class = "II") %>% 
  full_join(summPredictionsI) %>%
  mutate(Class = factor(Class, levels = c('I','II')),
         `HLA Allele` = gsub("HLA-", "", `HLA Allele`))
summPredictions2 <- summPredictions2II %>% mutate(Class = "II") %>% full_join(summPredictions2I) %>%
  mutate(Class = factor(Class, levels = c('I','II')),
         `HLA Allele` = gsub("HLA-", "", `HLA Allele`))

###
iedb <- fread('data/iedb20200602.epitope_full_v3.csv', skip = 1)
uniqueIedbAnt <- unique(iedb$Description)

# getMoreII <- lapply(unique(summPredictions$Peptide), function(x){
#   find <- lapply(uniqueIedbAnt, function(y){
#     if ()
#     grepl(x, y)
#   }) %>% unlist
#   find_filt <- Filter(Negate(is.null), find)
#   all_dist <- lapply(uniqueIedbAnt, function(y){
#     adist(x, y)
#   }) %>% unlist
#   if (length(find_filt)>0) {
#     return()
#   } else if(length(which(all_dist<2))>0) {
#     d <- data.frame(Peptide = x, closestPeptide = foundEpis$Peptide[which(all_dist<2)])
#     e <- left_join(d, foundEpis, by = c('closestPeptide' = 'Peptide'))
#     return(e)
#   } else {
#     d <- data.frame(Peptide = x, closestPeptide = as.character(NA))
#     return(d)
#   }
# })
###
tm <- fulldf %>% group_by(Country, World_Country) %>% summarise(n=n()) %>%
  dplyr::select(Country, World_Country)
countryPopFreq <- freqcountry %>% left_join(tm) %>% 
  dplyr::select(World_Country, Region, Locus, Superfamily, Allele, Country_Frequency, country_n)
colnames(countryPopFreq) <- c('Country','Region', 'Locus','HLA Superfamily', 'HLA Allele',
                              "Allele Frequency", "N Individuals (Sampled)")
#

save(world, fulldf, countryPopFreq,
     colors,
     proteinSequences3, 
     plotEpis, summPredictions, summPredictions2,
     # uniqueIedbAnt,
     file = 'data/forApp.Rda')




