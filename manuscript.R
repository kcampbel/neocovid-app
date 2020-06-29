library(data.table)
library(tidyverse)
library(broom)
library(patchwork)
setwd('~/neocovid-app/app')
load('data/forApp.Rda')

colors <- c('#7ed6df','#535c68','#badc58','#30336b',
            '#f6e58d','#6ab04c','#22a6b3','#f9ca24',
            '#4834d4','#95afc0','#e056fd','#ff7979',
            '#be2edd','#c7ecee','#f0932b','#eb4d4b',
            '#130f40','#ffbe76','#dff9fb','#686de0') #colors[c(8,4,12)]

##### ~~~~~ SUMMAR of pMHCs <500nM ~~~~~ #####
# Number of unique peptides
predictions[, list(alleles = unique(Peptide)), by = 'Class'][, list(n = .N), by = 'Class']
# Number of Alleles associated with predictions
predictions[, list(alleles = unique(hlaalleles)), by = 'Class'][, list(n = .N), by = 'Class']

# Overlap of Class I and II peptides
cIpeptides <- unique(predictions[Class == 'I']$Peptide)
cIIpeptides <- unique(predictions[Class == 'II']$Peptide)

# Class I that is nested in a Class II peptide
subString <- lapply(cIpeptides, function(cI){ any(grepl(cI, cIIpeptides)) })
table(unlist(subString))
# Class II that contain a Class I peptide
subString <- lapply(cIIpeptides, function(cII){ any(unlist(lapply(cIpeptides, function(cI) { grepl(cI, cII) } ))) })
table(unlist(subString))

# Peptides per gene
geneLengths <- unique(plotProteins[, c('Gene','Length')])
set.seed(771)
geneLengths$color <- sample(colors, 12, replace = FALSE)
#
peptideByGene <- merge(predictions[, list(nAlleles = .N), by = c('Class','Gene','Peptide', 'peptide.polar_start')], geneLengths, by = "Gene")
nByGenePerClass <- peptideByGene[, list(nPeptides = .N), by = c('Class','Gene','Length','color')][!is.na(Gene)]
nByGene <- nByGenePerClass[, list(nPeptides = sum(nPeptides)), by = c('Gene','Length','color')]
glance(lm(nByGene$nPeptides ~ nByGene$Length))

compareNPeptidesPerGene <- rbind(nByGene, nByGenePerClass, fill = TRUE)
tm <- geneLengths[, list(Class = c('I','II',NA)), by = c("Gene","Length")]
compareNPeptidesPerGene <- data.table(left_join(tm, compareNPeptidesPerGene))
compareNPeptidesPerGene[, nPeptides := ifelse(is.na(nPeptides), 0, nPeptides)]
ggplot(compareNPeptidesPerGene,
       aes(x = Length, y = nPeptides+1, colour = Class, group = Class)) +
  geom_smooth(method = 'lm', alpha = 0.4) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = colors[c(8,4)], na.value = colors[12]) +
  # scale_x_continuous(trans = 'log2') + scale_y_continuous(trans = 'log2') +
  # scale_x_log10() + scale_y_log10() +
  labs(x = "Length (AA)", y = "N Peptides") +
  theme(text = element_text(size = 12), title = element_text(size = 12))
ggsave('../figs/FigureS1.pdf', height = 4, width = 4)
#

peptideByGene$ClassX <- paste0("Class ", peptideByGene$Class)
dummy <- peptideByGene[, list(y = max(nAlleles)*1.05), by = c('Class', 'ClassX','Gene')]
a <- ggplot(data = peptideByGene, aes(colour = Gene)) +
  # Peptides
  geom_segment(aes(x = peptide.polar_start-100, xend = peptide.polar_start-100, 
               y = 0, yend = nAlleles)) +
  scale_color_manual(values = geneLengths$color, na.value = 'grey80') +
  #
  geom_blank(data = dummy, aes(y = y)) +
  scale_x_continuous(limits = c(0, max(plotProteins$polar_end))) +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "N Alleles", colour = "Viral Gene") +
  facet_grid(ClassX ~ ., scales = 'free_y') +
  guides(colour = guide_legend(nrow = 1)) +
  theme_bw() + 
  theme(text = element_text(size = 12), title = element_text(size =12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.position = 'none')
# 
b <- ggplot(data = plotProteins, aes(fill = Gene)) +
  geom_hline(yintercept = -25, colour = 'grey80', size = 2) +
  # Proteins
  geom_rect(aes(xmin = polar_start-100, xmax = polar_end-100, 
                ymin = -40, ymax = -10)) +
  scale_fill_manual(values = geneLengths$color, na.value = 'grey80') +
  # Chains
  # geom_rect(aes(xmin = chain_polar_start-100, xmax = chain_polar_end-100, 
                # ymin = -50, ymax = 0), fill = 'grey80', alpha = 0.3) +
  # geom_point(data = forPlot, aes(x = peptide.polar_start, y = 500-`Median Affinity (nM)`, colour = Gene, shape = hlagene), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = geneLengths$color, na.value = 'grey80', guide = FALSE) +
  #
  scale_x_continuous(limits = c(0, max(plotProteins$polar_end))) +
  scale_y_continuous(limits = c(-50, 0)) +
  labs(y = "Predicted Binding Affinity (nM)", x = "Position in Viral Proteome",
       colour = "Viral Gene") +
  facet_grid(' ' ~ ., scales = 'free_y', space = 'free_y') +
  guides(fill = guide_legend(nrow = 1)) +
  theme_classic() +
  theme(text = element_text(size = 10), title = element_text(size = 12),
        strip.background = element_rect(colour = NA, fill = 'white'),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = 'bottom', legend.justification = 'left')
#
a / b  + plot_layout(heights = c(5,1))
ggsave('../figs/Figure1B.pdf', height = 4, width = 9.56)


##### ~~~~~ IEDB COMPARISON ~~~~~ #####
setkey(comparePublished, Label)
iedb <- comparePublished['IEDB']
setkey(iedb, Comparison)
table(iedb['Exact Match', list(peplen = nchar(`Published Epitope`)), by = 'Published Epitope']$peplen)
table(iedb['Close Match', list(peplen = nchar(`Published Epitope`)), by = 'Published Epitope']$peplen)
table(iedb['Nested', list(peplen = nchar(`Published Epitope`)), by = 'Published Epitope']$peplen)
table(iedb[Comparison %in% c('Close Match','Nested'), list(peplen = nchar(`Published Epitope`)), by = 'Published Epitope']$peplen)
iedb[Comparison %in% c('Exact Match','Close Match','Nested'), list(nPeps = length(unique(`Published Epitope`))), by = `HLA Restriction`]
table(iedb[Comparison %in% c('Exact Match','Close Match','Nested'), list(nPeps = paste0(list(unique(`HLA Restriction`)))), by = `Published Epitope`]$nPeps)
#
setkey(predictions, Peptide)
otherHlas <- predictions[Peptide %in% iedb$Peptide, list(nHLAs = .N, uniqueHLAs = list(`HLA Allele`)), by = c('Peptide', 'Class')]
otherHlas[, list(min = min(nHLAs), median = median(nHLAs), max = max(nHLAs)), by = 'Class']

##### ~~~~~ Grifoni COMPARISON ~~~~~ #####
setkey(comparePublished, Source)
grifoni <- comparePublished['Grifoni Cell Host & Microbe 2020']
table(grifoni[, list(peplen = nchar(`Published Epitope`)), by = 'Published Epitope']$peplen)
table(grifoni[, list(peplen = nchar(Peptide)), by = 'Peptide']$peplen)
table(grifoni$Gene)
table(predictions[! Peptide %in% c(grifoni$Peptide, iedb$Peptide), list(n = .N), by = c('Peptide','Class')]$Class)

##### ~~~~~ HLA profile ~~~~~ #####
nByHLA <- predictions[, list(nPeptides = .N), by = c('Class','HLA Allele')]
nByHLA[, list(min = min(nPeptides), median = median(nPeptides), max = max(nPeptides)), by = c('Class')]
predictions[, list(nPeptides = .N), by = c('Class','hlagene')]




