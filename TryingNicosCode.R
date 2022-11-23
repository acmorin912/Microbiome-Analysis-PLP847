# importing other tables 
otus_ITS_uparse_R1 <-
  read.delim(
    "otu_table_ITS_UPARSE_R1.txt",
    header = TRUE,
    row.names = 1,)

metadata_ITS_uparse_R1 <-
  read.delim(
    "root_fun_map_soybean_2021.txt",
    row.names = 1,
    header = TRUE,)

otus_seq_ITS_uparse_R1 <-
  readDNAStringSet(
    "otus_R1.fasta",
    format = "fasta",
    seek.first.rec = TRUE,
    use.names = TRUE)

physeq_ITS_uparse <- phyloseq(otu_table(otus_ITS_uparse_R1, taxa_are_rows = TRUE),
                              sample_data(metadata_ITS_uparse_R1),
                              tax_table(as.matrix(taxonomy_ITS08_filt)),
                              otus_seq_ITS_uparse_R1) 

physeq_ITS_uparse
str(physeq_ITS_uparse)
head(sample_data(physeq_ITS_uparse))
tax_table(physeq_ITS_uparse)[tax_table(physeq_ITS_uparse)==""]<- NA
head(tax_table(physeq_ITS_uparse))

# checking the phyloseq object
sort(unique(as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom)) # everything is Fungi
nrow(as.data.frame(tax_table(physeq_ITS_uparse))[as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom!="Fungi",])

df_fungi_uparse <-
  as.data.frame(as.matrix(sample_data(physeq_ITS_uparse)))
df_fungi_uparse$LibrarySize <- sample_sums(physeq_ITS_uparse)
df_fungi_uparse <-
  df_fungi_uparse[order(df_fungi_uparse$LibrarySize), ]
df_fungi_uparse$Index <-
  seq(nrow(df_fungi_uparse)) # sample numbering
df_fungi_uparse

# reorder the factor
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

lib_size_all
    