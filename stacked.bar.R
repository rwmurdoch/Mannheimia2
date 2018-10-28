### this script is crap ###
### this script has been replaced with stacked.bar.2 and is no longer useful at all

library("magrittr")
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)

#making an additional metadata column that combines CvI and Site
variable1 = as.character(get_variable(physeq.f, "C_vs_I"))
variable2 = as.character(get_variable(physeq.f, "Site"))
sample_data(physeq.f)$NewPastedVar <- mapply(paste0, variable1, variable2, 
                                            collapse = "_")
merge_samples(physeq.f, "NewPastedVar")

#makes a new physeq object with only phylum information
physeq_phylum <- tax_glom(physeq.f, taxrank = "phylum")

#aggregates by infection and site
physeq_phylum_CS = merge_samples(physeq_phylum, "NewPastedVar")
sample_data(physeq_phylum_CS)$NewPastedVar <- factor(sample_names(physeq_phylum_CS))
physeq_phylum_CS = transform_sample_counts(physeq_phylum_CS, function(x) {x/sum(x)} )
physeq_phylum_CS = psmelt(physeq_phylum_CS)
physeq_phylum_CS = filter(physeq_phylum_CS, Abundance > 0.02) 
physeq_phylum_CS = arrange(physeq_phylum_CS, phylum) 



#make a pallete
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

#plotting crude CvsI
ggplot(physeq_phylum_CS, aes(x = C_vs_I, y = Abundance, fill = phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum composition of all samples by site and infection status") 



#save for later
ggplot(physeq_phylum_C, aes(x = Study_day, y = Abundance, fill = phylum)) + 
    geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n")  +
  ggtitle("Phylum Composition of Lake Erie \n Bacterial Communities by Sampling Site") 
