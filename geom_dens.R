plot_ordination(physeq.NL, ordNL.CCA, color = "C_vs_I", axes=c(1,2)) + 
  geom_point(mapping = aes()) +
  ggtitle("NL: CCA: WUniFrac") +  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=4) 

plot_ordination(physeq.NL, ordNL.CCA, color = "C_vs_I", axes=c(1,2)) + 
  geom_point(mapping = aes()) +
  ggtitle("NL: CCA: WUniFrac") +  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=4) + geom_density2d()

pNL = plot_ordination(physeq.NL, ordNL.CCA, color = "C_vs_I", axes=c(1,2)) + 
  geom_point(mapping = aes()) +
  ggtitle("NL: CCA: WUniFrac") +  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=4) 

p0 = ggplot(pNL$data, pNL$mapping) + geom_density_2d()
p0
