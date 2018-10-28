######## THIS SCRIPT IS NOT USED #############

# this follows on taxa.time.plot.R but is directed instead at a single site
taxa.diff = c("Selenomonas_3","Citrobacter","Succiniclasticum","Haemophilus","Comamonas","Mannheimia","Selenomonas_3")
genus.glom = tax_glom(physeq.f, taxrank = "genus")
physeq.diff = subset_taxa(genus.glom, genus == "Selenomonas_3"|genus =="Citrobacter"|genus=="Succiniclasticum"
                          |genus=="Haemophilus"|genus=="Comamonas"|genus=="Mannheimia")


## this script allows for plotting taxa over time with display of every sample and..
## display of a moving average with variance shadow

# define your physeq or subset here
set <- physeq.NS


# define your target taxon here
TAXON <- "Mannheimia"
# define chart title here
TITLE = "Relative abundance of Mannheimia in NS"
# enter a name for the linear scale plot output file
PLOT = "Mannheimia.NL"

## You will have to alter manually the variables you want to facet by in the ggplot call at the end of the script

physeq.pasturellaceae = subset_taxa(set, genus == TAXON)

### now I need to subset by infection status and turn each into relative abundance as derived from sum of 
### counts in each sample

# first turn the taxon-specific OTU tables into dataframes

OTU.P = as(otu_table(physeq.pasturellaceae), "matrix")
OTU.P = as.data.frame(OTU.P)
P.sums = as.data.frame(colSums(OTU.P))

# then get the entire study as an OTU data frame and pull out sums for each sample
OTU.Site = as(otu_table(set), "matrix")
OTU.Site = as.data.frame(OTU.Site)
Site.sums = as.data.frame(colSums(OTU.Site))

#and divide the taxon counts by total sequence sums by sample

OTU.P.all.rel = P.sums / Site.sums
#replace zeros with a minimal value so as to allow for geom smoothing
OTU.P.all.rel[OTU.P.all.rel==0] <- 0.000001

# now I have full Taxon OTU tables in relative abundance format

### moving towards a by-sample plot of sum Pasteurellaceae by day


# just add the rel Pasteurellaceae to sub-setted metadata tables
metadata.Site = sample_data(set)
metadata.Site$Study_day <- as.numeric(sample_data(set)$Study_day) #make sure day is numeric, helps with the plotting here
OTU.P.all.rel.sum.meta = cbind(metadata.Site, OTU.P.all.rel)
colnames(OTU.P.all.rel.sum.meta) [6] = "Relative.Abundance"

# now I have added rel-abundance of Pasteurellaceae to the metadata file
# I have to figure a nice wqy to chart this out

library(ggplot2)
library(scales)

my.color = "#6494e0"

df = OTU.P.all.rel.sum.meta
options(scipen=999)
theme_set(theme_bw())
ggo <- ggplot(df, aes(x=Study_day, y=Relative.Abundance)) + 
  geom_point(color = my.color) + 
  geom_smooth(method = "auto", se=TRUE, level = 0.66) + 
  xlim(c(0, 15)) + 
  ylim(c(-.2, 1)) +
  scale_y_log10(limits=c(0.0000001,1)) +
  labs(subtitle="66% confidence interval", 
       y="log10 relative abundance", 
       x="Study Day", 
       title=TITLE) +
  facet_grid(rows=vars(C_vs_I), cols=vars(Site))
plot(ggo)


### this is an attemp to overlay instead of facet C_vs_I

gglog <- ggplot(df, aes(x=Study_day, y=Relative.Abundance, color = C_vs_I)) + 
  geom_smooth(method = "auto", se=TRUE, level = 0.66) +
  xlim(c(1, 6)) + 
  ylim(c(-.2, 1)) +
  scale_y_log10(limits=c(0.0000001,1)) +
  labs(subtitle="66% confidence interval", 
       y="log10 relative abundance", 
       title=TITLE) +
  facet_grid(cols=vars(Site))


### this is same as above but with linear y-axis
png(PLOT, x=600, y=600)

gglin <- ggplot(df, aes(x=Study_day, y=Relative.Abundance, color = C_vs_I)) + 
  geom_smooth(method = "auto", se=TRUE, level = 0.66) +
  xlim(c(1, 6)) + 
  ylim(c(-.2, 0.75)) +
  labs(subtitle="66% confidence interval", 
       y="relative abundance", 
       x="Study Day", 
       title=TITLE) +
  facet_grid(cols=vars(Site))

plot(gglin)
dev.off()


