library(ggpubr)
library(tidyverse)

# importing other tables 
otus_ITS_uparse_R1 <-
  read.delim(
    "otu_table_ITS_UPARSE_R1.txt",
    header = TRUE,
    row.names = 1,)

metadata_ITS_uparse_R1 <-
  read.delim(
    "root_fun_map_soybean_2021_1.txt",
    row.names = 1,
    header = TRUE,)

otus_seq_ITS_uparse_R1 <-
  readDNAStringSet(
    "otus_R1.fasta",
    format = "fasta",
    seek.first.rec = TRUE,
    use.names = TRUE)

ITS_taxonomy<-read.delim("constax_taxonomy.txt",
                         header=TRUE, 
                         row.names=1)

physeq_ITS_uparse <- phyloseq(otu_table(otus_ITS_uparse_R1, taxa_are_rows = TRUE),
                              sample_data(metadata_ITS_uparse_R1),
                              tax_table(as.matrix(ITS_taxonomy)),
                              otus_seq_ITS_uparse_R1) 

physeq_ITS_uparse
str(physeq_ITS_uparse)
head(sample_data(physeq_ITS_uparse))
tax_table(physeq_ITS_uparse)[tax_table(physeq_ITS_uparse)==""]<- NA
head(tax_table(physeq_ITS_uparse))

# checking the phyloseq object
sort(unique(as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom)) # not everything is Fungi
# remove non fungi 
#remove chloropast,mitochondria,cyanobacteria (Moriah added code below - from Ried's script)
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Phylum!="Chloroplast")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Class!="Chloroplast")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Order!="Chloroplast")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Family!="Chloroplast")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Genus!="Chloroplast")
sort(unique(as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom))

physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Phylum!="Mitochondria")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Class!="Mitochondria")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Order!="Mitochondria")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Family!="Mitochondria")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Genus!="Mitochondria")
sort(unique(as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom))

physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Phylum!="Cynobacteria")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Class!="Cynobacteria")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Order!="Cynobacteria")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Family!="Cynobacteria")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Genus!="Cynobacteria")
sort(unique(as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom))

physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Kingdom!="Anthophyta")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Kingdom!="Rhizaria_1")
physeq_ITS_uparse <- subset_taxa(physeq_ITS_uparse, Kingdom!="Viridiplantae_1")
sort(unique(as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom))
# end of Moriah's additiong
nrow(as.data.frame(tax_table(physeq_ITS_uparse))[as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom!="Fungi",])

df_fungi_uparse <-
  as.data.frame(as.matrix(sample_data(physeq_ITS_uparse)))
df_fungi_uparse$LibrarySize <- sample_sums(physeq_ITS_uparse)
df_fungi_uparse <-
  df_fungi_uparse[order(df_fungi_uparse$LibrarySize), ]
df_fungi_uparse$Index <-
  seq(nrow(df_fungi_uparse)) # sample numbering
df_fungi_uparse

# Delete samples with a mean of less than 1000 (Moriah added this)
samplesover1000_all <- subset_samples(physeq_ITS_uparse, sample_sums(physeq_ITS_uparse) > 1000)
any(taxa_sums(samplesover1000_all) == 0)
set.seed(81)
rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))

# reorder the factor - (Moriah wants to know what are these levels???)
df_fungi_uparse$Site <- factor(
  df_fungi_uparse$Site,
  levels = c(
    "Citta_di_Castello",
    "Ripabianca",
    "San_Giovanni_dAsso",
    "control"))

ggarrange(
  ggplot(data=df_fungi_uparse, aes(x=Index, y=LibrarySize, color=is.neg)) +
    geom_point(alpha =0.7) +
    labs(title="Fungi", x="Sample number", y="Read number") + 
    theme_classic() +
    scale_colour_manual("Negative control", values = c("grey", "red")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)),
  ggplot(df_fungi_uparse, aes(x = LibrarySize)) + # Histogram of sample read counts
    geom_histogram(color = "indianred", fill = "indianred", binwidth = 1000) +
    theme_classic() +
    #facet_grid(~Treatment, scales = "free_x", space="free_x") +
    labs(title=NULL, x="Read number", y="Sample number") +
    theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)),
  labels = c("A", "B", "C", "D"),
  widths = c(1,1,1,1),
  align = "v", ncol = 2, nrow = 2,
  common.legend = TRUE,
  legend = "bottom") -> lib_size_all

glib_size_all
    