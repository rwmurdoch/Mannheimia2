#otu	<- read.table(file = "phyloseq/otu_table.txt",	header = TRUE)
#head(otu)

#tax	<- read.table(file="phyloseq/taxonomy.tsv", sep	=	'\t',	header = TRUE)
#head(tax)

#merged_file	<-	merge(otu, tax, by.x = c("OTUID"), by.y=c("OTUID"))
#head(merged_file)

#write.table(merged_file, file	=	"combined_otu_tax",	sep = '\t', col.names = 
#  TRUE,	row.names	=	FALSE)

#	It	seems	tedious	but	you	need	to	open	the	merged	.txt	file	in	excel	and	split	
#into	two	files:	one	for	taxonomy	(containing	only	the	columns	OTUID	and	
#taxonomic	info)	and	the	other	for	the	OTU	matrix	(containing	only	OTUID	and	
#abundances	in	each	sample).	Note:	for	the	taxonomy	file,	you	need	to	use	data	
#—>	text-to-columns	in	Excel	and	separate	on	semicolon	to	get	columns	for	
#kingdom,	phylum,	class,	etc…	once	you	make	these	two	separate	files	in	excel,	
#save	each	as	a	.csv	

library("ggplot2")
library("phyloseq")
library("ape")

otu_table = read.csv(file = "otu_table.csv", sep=",", row.names=1)
otu_table = as.matrix(otu_table)

taxonomy = read.csv("taxonomy_new.csv", sep=",", row.names=1)
taxonomy = as.matrix(taxonomy)

#metadata2 is the original
#metadata3.txt has new values for day
metadata = read.table("Sample_metadata3.txt", row.names=1, header = TRUE, stringsAsFactors = FALSE)
#this is new (180716), to avoid tagnumber and day being treated as integers
metadata$Tag_Number <- as.character(metadata$Tag_Number)
#181012 first try this new day as an integer
#metadata$Study_day <- as.character(metadata$Study_day)
metadata$C_vs_I <- as.factor(metadata$C_vs_I)

phy_tree = read_tree("phyloseq/tree.nwk")

OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata)
#	(tree	was	already	imported	as	a	phyloseq	object)

#	check	that	your	OTU	names	are	consistent	across	objects
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)

#	make	sure	files	have	the	same	sample	names	
sample_names(OTU)
sample_names(META)

#	merge	into	one	phyloseq	object
physeq	=	phyloseq(OTU,	TAX,	META,	phy_tree)
physeq





