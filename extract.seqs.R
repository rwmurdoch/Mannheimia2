### this script is aimed at extracting a subset of seqs in conjunction with the OTU table

## First, the taxon that you are interested in make a sub-table and turn into an OTU table with hashes as seq identifiers
library(phyloseq)

# write the desired output file name
output.file = "Pasteurellaceae.OTU.seq.table.csv"
output.fasta = "Pasteurellaceae.seqs.fasta"

# designate your target taxon
target.seq.set <- subset_taxa(physeq.f, family == "Pasteurellaceae")

# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(target.seq.set), "matrix")
# transpose if necessary
if(taxa_are_rows(physeq.f)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)
OTUdf = as.data.frame(t(OTUdf))

library(data.table)
OTUdt = setDT(OTUdf, keep.rownames = TRUE)[]
colnames(OTUdt)[1] <- "seq_name"

## Next, import the hash to seq table
## the target is a multifasta from the qiime2 seqs.qzv

library("Biostrings")

all.seqs <- readDNAStringSet("Sam_dada2_seqs.fasta")
seq_name = names(all.seqs)
sequence = paste(all.seqs)
all.seqs <- data.frame(seq_name, sequence)

## next we do a simple bind, bind the sequence onto the OTUdt object
OTU.seq.table <- dplyr::left_join(OTUdt, all.seqs, by = "seq_name")

## now we have to extract and bind the taxonomy
tax.lookup <- as.data.frame(tax_table(target.seq.set))
tax.lookup <- setDT(tax.lookup, keep.rownames = TRUE)[]
colnames(tax.lookup)[1] <- "seq_name"
OTU.seq.table <- dplyr::left_join(OTU.seq.table, tax.lookup, by = "seq_name")

## and finally write the table
write.csv(OTU.seq.table, output.file)

## next we also need a multifasta with informative sequence name
# pull over just species, seq_name, and sequence
seqs <- OTU.seq.table[,c("seq_name","species","sequence")]
seqs$name <- paste(seqs$species,seqs$seq_name,sep="_")
seqs <- seqs[,c(4,3)]
install.packages("seqRFLP")
library("seqRFLP")
df.fasta = dataframe2fas(seqs, file="df.fasta")
write(df.fasta,output.fasta)
