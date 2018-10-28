## Initial installation, probably don't bother with over and over ##
## very large package, long time to install ##

## method 1 ##
## does not work ##

#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')

## dependency; tibble? ##
#install.packages("tidyverse")

## method 2 ##
#source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
#local = TRUE)




library("phyloseq")
packageVersion("phyloseq")

library("ggplot2")
packageVersion("ggplot2")

library("scales")
packageVersion("scales")

library("grid")
packageVersion("grid")

## this is maybe not the best idea ##
theme_set(theme_bw())

plot_richness(physeq, x = "Site") + geom_boxplot()
plot_richness(physeq, x = "Study_day") + geom_boxplot()
plot_richness(physeq, x = "C_vs_I") + geom_boxplot()
plot_richness(physeq, x = "Tag_Number") + geom_boxplot()

plot_richness(physeq, x = "Site") 
plot_richness(physeq, x = "Study_day") 
plot_richness(physeq, x = "C_vs_I") 
plot_richness(physeq, x = "Tag_Number") 

ord

