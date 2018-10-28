## Facet wrapping C vs I by site with day shading ###

ord.d = ordinate(physeq.NL, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.NL, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I") + 
  ggtitle("NL, NMDS, bray-curtis") + facet_wrap("C_vs_I") +
  geom_point(aes(size = 3)) + scale_colour_gradient2(midpoint = 3.5, low = "red", mid = "white", high = "blue") + 
  theme_dark() +
  geom_text(mapping = aes(label = Tag_Number), size = 2.8, vjust = 1.9)
  
g

ord.d = ordinate(physeq.NS, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.NS, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I") + 
  ggtitle("NS, NMDS, bray-curtis") + facet_wrap("C_vs_I") +
  geom_point(aes(size = 3)) + scale_colour_gradient2(midpoint = 3.5, low = "red", mid = "white", high = "blue") + 
  theme_dark() +
  geom_text(mapping = aes(label = Tag_Number), size = 2.8, vjust = 1.9)
  
g

ord.d = ordinate(physeq.PH, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.PH, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I") + 
  ggtitle("PH, NMDS, bray-curtis") + facet_wrap("C_vs_I") +
  geom_point(aes(size = 3)) + scale_colour_gradient2(midpoint = 3.5, low = "red", mid = "white", high = "blue") + 
  theme_dark() +
  geom_text(mapping = aes(label = Tag_Number), size = 2.8, vjust = 1.9)
  
g

ord.d = ordinate(physeq.TS, method = "NMDS", distance = "bray")
g <- plot_ordination(physeq.TS, ord.d, color = "Study_day", axes = c(1,2), shape = "C_vs_I") + 
  ggtitle("TS, NMDS, bray-curtis") + facet_wrap("C_vs_I") +
  geom_point(aes(size = 3)) + scale_colour_gradient2(midpoint = 3.5, low = "red", mid = "white", high = "blue") + 
  theme_dark() +
  geom_text(mapping = aes(label = Tag_Number), size = 2.8, vjust = 1.9)
  
g
