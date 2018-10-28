ord.d = ordinate(physeq.NS, method = "CCA", distance = "wunifrac")
g <- plot_ordination(physeq, ord.d, color = "Tag_Number", axes = c(1,4)) + 
  geom_point(mapping = aes()) + 
  ggtitle("NS infected, CCA: WUniFrac") +
  scale_colour_brewer(type="qual", palette="Paired") +
  scale_fill_brewer(type="qual", palette="Dark2") +
  geom_point(size=5) 

g
g + geom_smooth(aes(col=Tag_Number), method = "auto", se=F)



study_day.NS = sample_data(physeq.NS)$Study_day

ord.d = ordinate(physeq.NS, method = "CCA", distance = "wunifrac")
g <- plot_ordination(physeq.NS, ord.d, color = "Tag_Number", axes = c(1,2), shape = "C_vs_I", label = "Study_day") + 
  ggtitle("NS infected, CCA: WUniFrac") +
  scale_colour_brewer(type="qual", palette="Spectral") +
  geom_point(size=3) + theme_dark() +
  geom_text(mapping = aes(label = Study_day), size = 5, vjust = 1.5) 

g

## attempting a constraqined vectorized ordination as shown by
## http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#constrained_ordinations

cap.ord = ordinate(physeq.NS, method = "CAP", distance = "bray", formula = ~ Study_day + Farm)
g <- plot_ordination(physeq.NS, cap.ord, color = "Tag_Number", axes = c(1,2), shape = "C_vs_I", label = "Study_day") + 
  ggtitle("NS infected, CCA: WUniFrac") +
  scale_colour_brewer(type="qual", palette="Spectral") +
  geom_point(size=3) + theme_dark() +
  geom_text(mapping = aes(label = Study_day), size = 5, vjust = 1.5) 

g

#now we add vectors

arrowmat = vegan::scores(cap.ord, display = "bp")
arrowdf = data.frame(labels = rownames(arrowmat), arrowmat)

arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# make a new graphic

g + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )

## checking other axes from pCCA above

g <- plot_ordination(physeq.NS, cap.ord, color = "Tag_Number", axes = c(5,6), shape = "C_vs_I", label = "Study_day") + 
  ggtitle("NS infected, CCA: WUniFrac") +
  scale_colour_brewer(type="qual", palette="Spectral") +
  geom_point(size=3) + theme_dark() +
  geom_text(mapping = aes(label = Study_day), size = 5, vjust = 1.5) 

g


### this starts to produce a useful figure:


ord.d = ordinate(physeq.NS, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.NS, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I", label = "Study_day") + 
  ggtitle("NS, NMDS, bray-curtis") +
  scale_fill_distiller(type="seq", palette="Spectral") +
  geom_point(size=3) + theme_dark() +
  geom_text(mapping = aes(label = Study_day), size = 5, vjust = 1.5) 

g

ord.d = ordinate(physeq.NL, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.NL, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I", label = "Study_day") + 
  ggtitle("NL, NMDS, bray-curtis") +
  scale_fill_distiller(type="seq", palette="Spectral") +
  geom_point(size=3) + theme_dark() +
  geom_text(mapping = aes(label = Study_day), size = 5, vjust = 1.5) 

g

ord.d = ordinate(physeq.PH, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.PH, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I", label = "Study_day") + 
  ggtitle("PH, NMDS, bray-curtis") +
  scale_fill_distiller(type="seq", palette="Spectral") +
  geom_point(size=3) + theme_dark() +
  geom_text(mapping = aes(label = Study_day), size = 5, vjust = 1.5) 

g

ord.d = ordinate(physeq.TS, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.TS, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I") + 
  ggtitle("TS, NMDS, bray-curtis") + 
  scale_fill_distiller(type="seq", palette="Spectral") +
  geom_point(size=3) + theme_dark() +
  geom_text(mapping = aes(label = Tag_Number), size = 3, vjust = 1.5) +
  geom_text(mapping = aes(label = Study_day), size = 3, vjust = 3)
g

## trying to do a more dramatic color scale for day
## this test works, has both size and a weird color scale

ord.d = ordinate(physeq.TS, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.TS, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I") + 
  ggtitle("TS, NMDS, bray-curtis") + 
  geom_point(aes(size = Study_day)) + scale_colour_gradientn(colours = terrain.colors(10)) + 
  theme_dark() +
  geom_text(mapping = aes(label = Tag_Number), size = 3, vjust = 1.5) +
  geom_text(mapping = aes(label = Study_day), size = 3, vjust = 3)
g



## improving on previoius, now with a bold color scale only



ord.d = ordinate(physeq.NL, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.NL, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I") + 
  ggtitle("NL, NMDS, bray-curtis") + 
  geom_point(aes(size = 3)) + scale_colour_gradient2(midpoint = 3.5, low = "red", mid = "white", high = "blue") + 
  theme_dark() +
  geom_text(mapping = aes(label = Tag_Number), size = 3, vjust = 1.5) +
  geom_text(mapping = aes(label = Study_day), size = 3, vjust = 3)
g

ord.d = ordinate(physeq.NS, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.NS, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I") + 
  ggtitle("NS, NMDS, bray-curtis") + 
  geom_point(aes(size = 3)) + scale_colour_gradient2(midpoint = 3.5, low = "red", mid = "white", high = "blue") + 
  theme_dark() +
  geom_text(mapping = aes(label = Tag_Number), size = 3, vjust = 1.5) +
  geom_text(mapping = aes(label = Study_day), size = 3, vjust = 3)
g

ord.d = ordinate(physeq.PH, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.PH, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I") + 
  ggtitle("PH, NMDS, bray-curtis") + 
  geom_point(aes(size = 3)) + scale_colour_gradient2(midpoint = 3.5, low = "red", mid = "white", high = "blue") + 
  theme_dark() +
  geom_text(mapping = aes(label = Tag_Number), size = 3, vjust = 1.5) +
  geom_text(mapping = aes(label = Study_day), size = 3, vjust = 3)
g

ord.d = ordinate(physeq.TS, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.TS, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I") + 
  ggtitle("TS, NMDS, bray-curtis") + 
  geom_point(aes(size = 3)) + scale_colour_gradient2(midpoint = 3.5, low = "red", mid = "white", high = "blue") + 
  theme_dark() +
  geom_text(mapping = aes(label = Tag_Number), size = 3, vjust = 1.5) +
  geom_text(mapping = aes(label = Study_day), size = 3, vjust = 3)
g



#### attempting to convert to a partially constrained CCA, removing effect of tag number

#first have to
ord.d = cca.phyloseq(physeq.TS, formula = Tag_Number, "CCA")
