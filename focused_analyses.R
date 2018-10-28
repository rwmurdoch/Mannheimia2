#subset starting from site subsets to only those from infected animals on days 1 and 2
physeq.NL.I = subset_samples(physeq.NL, C_vs_I=="I")
physeq.NL.I.day12 = subset_samples(physeq.NL.I, Study_day%in%c(1,2))
physeq.NS.I = subset_samples(physeq.NS, C_vs_I=="I")
physeq.NS.I.day12 = subset_samples(physeq.NS.I, Study_day%in%c(1,2))
physeq.PH.I = subset_samples(physeq.PH, C_vs_I=="I")
physeq.PH.I.day12 = subset_samples(physeq.PH.I, Study_day%in%c(1,2))
physeq.TS.I = subset_samples(physeq.TS, C_vs_I=="I")
physeq.TS.I.day12 = subset_samples(physeq.TS.I, Study_day%in%c(1,2))

plot_richness(physeq.NL.I, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
p <- plot_richness(physeq.NL.I, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
(div.all.site <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                                  alpha = 0.1)) + ggtitle("Overall Alpha-diversity") + geom_smooth(aes(col=Tag_Number))

plot_richness(physeq.NS.I, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
p <- plot_richness(physeq.NS.I, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
(div.all.site <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                                  alpha = 0.1)) + ggtitle("Overall Alpha-diversity") + geom_smooth(aes(col=Tag_Number))

plot_richness(physeq.PH.I, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
p <- plot_richness(physeq.PH.I, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
(div.all.site <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                                  alpha = 0.1)) + ggtitle("Overall Alpha-diversity") + geom_smooth(aes(col=Tag_Number))

plot_richness(physeq.TS.I, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
p <- plot_richness(physeq.TS.I, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
(div.all.site <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                                  alpha = 0.1)) + ggtitle("Overall Alpha-diversity") + geom_smooth(aes(col=Tag_Number))

## Adonis test
library("vegan")

# permanova, test DAY 1 VS DAY 2
dist.wuni.physeq.NL.I.day12 <- phyloseq::distance(physeq.NL.I.day12, method = "wunifrac")
dist.wuni.physeq.NS.I.day12 <- phyloseq::distance(physeq.NS.I.day12, method = "wunifrac")
dist.wuni.physeq.PH.I.day12 <- phyloseq::distance(physeq.PH.I.day12, method = "wunifrac")
dist.wuni.physeq.TS.I.day12 <- phyloseq::distance(physeq.TS.I.day12, method = "wunifrac")

adonis.physeq.NL.I.day12 <- data.frame(sample_data(physeq.NL.I.day12))
adonis.physeq.NS.I.day12 <- data.frame(sample_data(physeq.NS.I.day12))
adonis.physeq.PH.I.day12 <- data.frame(sample_data(physeq.PH.I.day12))
adonis.physeq.TS.I.day12 <- data.frame(sample_data(physeq.TS.I.day12))

adonis(dist.wuni.physeq.NL.I.day12 ~ Study_day, data = adonis.physeq.NL.I.day12)
adonis(dist.wuni.physeq.NS.I.day12 ~ Study_day, data = adonis.physeq.NS.I.day12)
adonis(dist.wuni.physeq.PH.I.day12 ~ Study_day, data = adonis.physeq.PH.I.day12)
adonis(dist.wuni.physeq.TS.I.day12 ~ Study_day, data = adonis.physeq.TS.I.day12)

#ordinations
ord.d = ordinate(physeq.NL.I, method = "CCA", distance = "wunifrac")
g <- plot_ordination(physeq, ord.d, color = "Tag_Number") + 
  geom_point(mapping = aes()) + 
  ggtitle("NL infected, CCA: WUniFrac") +
  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) 

g + geom_smooth(aes(col=Tag_Number), method = "lm", se=F)

ord.d = ordinate(physeq.NS.I, method = "CCA", distance = "wunifrac")
g <- plot_ordination(physeq, ord.d, color = "Tag_Number") + 
  geom_point(mapping = aes()) + 
  ggtitle("NS infected, CCA: WUniFrac") +
  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) 

g + geom_smooth(aes(col=Tag_Number), method = "lm", se=F)

ord.d = ordinate(physeq.NS, method = "CCA", distance = "wunifrac")
g <- plot_ordination(physeq, ord.d, color = "Tag_Number") + 
  geom_point(mapping = aes()) + 
  ggtitle("NS infected, CCA: WUniFrac") +
  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) 

g + geom_smooth(aes(col=C_vs_I), method = "lm", se=F)

ord.d = ordinate(physeq.NS, method = "CCA", distance = "wunifrac")
g <- plot_ordination(physeq.NS, ord.d, color = "C_vs_I") + 
  geom_point(mapping = aes()) + 
  ggtitle("NS infected, CCA: WUniFrac") +
  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) 

g + geom_smooth(aes(col=Tag_Number), method="lm")

ord.d = ordinate(physeq.PH.I.day12, method = "CCA", distance = "wunifrac")
plot_ordination(physeq, ord.d, color = "Study_day") + 
  geom_point(mapping = aes()) + 
  ggtitle("PH days 1 and 2, CCA: WUniFrac") +
  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) 

ord.d = ordinate(physeq.TS.I.day12, method = "CCA", distance = "wunifrac")
plot_ordination(physeq, ord.d, color = "Study_day") + 
  geom_point(mapping = aes()) + 
  ggtitle("TS days 1 and 2, CCA: WUniFrac") +
  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) 

#subset starting from site subsets to only those from NONINFECTED animals on days 1 and 2
physeq.NL.C = subset_samples(physeq.NL, C_vs_I=="C")
physeq.NL.C.day12 = subset_samples(physeq.NL.C, Study_day%in%c(1,2))
physeq.NS.C = subset_samples(physeq.NS, C_vs_I=="C")
physeq.NS.C.day12 = subset_samples(physeq.NS.C, Study_day%in%c(1,2))
physeq.PH.C = subset_samples(physeq.PH, C_vs_I=="C")
physeq.PH.C.day12 = subset_samples(physeq.PH.C, Study_day%in%c(1,2))
physeq.TS.C = subset_samples(physeq.TS, C_vs_I=="C")
physeq.TS.C.day12 = subset_samples(physeq.TS.C, Study_day%in%c(1,2))

plot_richness(physeq.NL.C, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
p <- plot_richness(physeq.NL.C, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
(div.all.site <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                                  alpha = 0.1)) + ggtitle("Overall Alpha-diversity") + geom_smooth(aes(col=Tag_Number))

plot_richness(physeq.NS.C, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
p <- plot_richness(physeq.NS.C, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
(div.all.site <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                                  alpha = 0.1)) + ggtitle("Overall Alpha-diversity") + geom_smooth(aes(col=Tag_Number))

plot_richness(physeq.PH.C, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
p <- plot_richness(physeq.PH.C, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
(div.all.site <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                                  alpha = 0.1)) + ggtitle("Overall Alpha-diversity") + geom_smooth(aes(col=Tag_Number))

plot_richness(physeq.TS.C, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
p <- plot_richness(physeq.TS.C, x= "Study_day", "Tag_Number", measures=c("Chao1", "Shannon"))
(div.all.site <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                                  alpha = 0.1)) + ggtitle("Overall Alpha-diversity") + geom_smooth(aes(col=Tag_Number))
