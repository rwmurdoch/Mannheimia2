### 180904 This script works finally, in so far as it will test diff abundance across days 2-6
### then filter out taxa that are already divergent at day 1
### what I learned from this is that after false discovery rate is factored in, very few differential taxa remain

# TO DO: format and save the tables
# consider funenling this to a taxatimeplot for genera of interest

## This is an attempt, after many rounds of analysis, to come up with a better screen for differences
## introduce a control for starting conditions by dividing by the day 1 abundance
## currently, it excludes day 1 for the comparisons as an initial screen
## the pipeline then should move to the stacked.bar.2 script which is then directed at taxa of interest by site

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
library("DESeq2")

## Deseq does not use direct comparisons, it calculates curves and deviations
## changing to relative abundance not necessary

#agglomerate to genus level, only run once, this takes a long time
#physeq.genus = tax_glom(physeq.f, taxrank = "genus")
#physeq.genus.26 = subset_samples(physeq.genus, Study_day != 1)
#physeq.genus.1 = subset_samples(physeq.genus, Study_day == 1)
#type RAW to see non-normalized result
#type FILTER to see result with wondy day 1 taxa removed

titleD = "TS days 2-6: Infected / Control"
target.set = subset_samples(physeq.genus.26, Site == "TS")
target.set.1 = subset_samples(physeq.genus.1, Site == "TS")
UP = 0.5 #sets upper and lower limits for log2 filtration
LOW = -0.5


diagdds = phyloseq_to_deseq2(target.set, ~ C_vs_I)

# method to deal with too many zeros, skip if not needed
geoMeans = apply(counts(diagdds, normalized = FALSE), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

# always use this
diagdds = DESeq(diagdds, test="Wald", fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.1
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(target.set)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$phylum = factor(as.character(sigtab$phylum), levels=names(x))

# family order
#x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
#x = sort(x, TRUE)
#sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
RAW <- ggplot(sigtab, aes(x=genus, y=log2FoldChange, color=family)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))  + ggtitle(titleD)

#write.csv(sigtab, "TS.days26.diffabund.csv")


####### this is aimed at normalizing to starting abundance, pretty difficult to deal with and NOT COMPLETE
# the general strategy is to first filter the whole genusglom taxtable...
# removing all genera where C and I are different to start with
# but how to define different? Attempt to calculate ratios
# instead of filtering out everything with alpha value less than 0.5...
# filter out everything outside of -0.5 to 0.5 log2foldchange

# this is simple definition of a function that will replace 0 with 1
rep.zero = function(x){
  ifelse(x == 0, 1, x)
}

# substitute the zeros with 1's, leaving other counts the same
target.set.2 <- transform_sample_counts(target.set.1, rep.zero)
diagdds.1 = phyloseq_to_deseq2(target.set.2, ~ C_vs_I)

# method to deal with too many zeros, skip if not needed
#geoMeans = apply(counts(diagdds, normalized = FALSE), 1, gm_mean)
#diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

diagdds.1 = DESeq(diagdds.1, fitType="local")


res.1 = results(diagdds.1, cooksCutoff = FALSE)

sigtab.1 = res.1[which(res.1$log2FoldChange < UP), ] 
sigtab.2 = sigtab.1[which(sigtab.1$log2FoldChange > LOW), ]
#sigtab.2 = sigtab.2[which(sigtab.2$)] #hold this thought
sigtab.2 = cbind(as(sigtab.2, "data.frame"), as(tax_table(target.set.2)[rownames(sigtab.2), ], "matrix"))
head(sigtab.2)

dim(sigtab.2)

include <- rownames(sigtab.2)


############################################################################################################
##### this next part filters the deseq2 results table from days2-6 using the include hash list##############
############################################################################################################

sigtab.f = subset(sigtab, rownames(sigtab) %in% rownames(sigtab.2))

head(sigtab.f)

dim(sigtab.f)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab.f$log2FoldChange, sigtab.f$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.f$phylum = factor(as.character(sigtab.f$phylum), levels=names(x))

# family order
#x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
#x = sort(x, TRUE)
#sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
FILTER <- ggplot(sigtab.f, aes(x=genus, y=log2FoldChange, color=family)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))  + ggtitle(titleD)
