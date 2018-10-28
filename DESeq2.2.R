## this is a second attempt to run DeSeq2 on the phyloseq data
## http://joey711.github.io/phyloseq-extensions/DESeq2.html
## results are largely the same, but I think this is more appropriate overall

##I am going to agglomerate to genus level, I don't think OTU level is relevant
##https://github.com/joey711/phyloseq/issues/683

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#physeq.genus = tax_glom(physeq.f, "genus")

library("DESeq2")
titleD = "PH day 4: Infected / Control"
physeq.NL.day4 = subset_samples(physeq.genus, Site == "NL")
physeq.NL.day4 = subset_samples(physeq.NL.day4, Study_day!=1)
diagdds = phyloseq_to_deseq2(physeq.NL.day4, ~ C_vs_I)
diagdds = DESeq(diagdds, test="Wald", fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq.NL.day4)[rownames(sigtab), ], "matrix"))
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
x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
ggplot(sigtab, aes(x=genus, y=log2FoldChange, color=phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))  + ggtitle(titleD)
