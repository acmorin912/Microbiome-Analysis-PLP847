library(phyloseq)
library(Biostrings)
library(tidyverse)
library(vegan)
library(decontam)
library(ggpubr)
#import
mapping <- read.delim("FIXED_FIXED_root_fun_map_soybean_2021.txt", row.names = 1, header= TRUE, sep = "\t")
head(mapping)

otutable<- read.delim("otu_table_ITS_UPARSE_R1.txt", row.names = 1, header= TRUE, sep = "\t")
head(otutable)

taxonomy<- read.delim("constax_taxonomy.txt", row.names = 1, header= TRUE, sep= "\t")
head(taxonomy)

sequences<- readDNAStringSet("otus_R1.fasta", format = "fasta", seek.first.rec = TRUE, use.names = TRUE)
sequences
taxonomy[, "Kingdom"] <- as.factor(gsub("_1", "", taxonomy[, "Kingdom"]))
taxonomy[, "Phylum"] <- as.factor(gsub("_1", "", taxonomy[, "Phylum"]))
taxonomy[, "Class"] <- as.factor(gsub("_1", "", taxonomy[, "Class"]))
taxonomy[, "Order"] <- as.factor(gsub("_1", "",taxonomy[, "Order"]))
taxonomy[, "Family"] <- as.factor(gsub("_1", "", taxonomy[, "Family"]))
taxonomy[, "Genus"] <- as.factor(gsub("_1", "", taxonomy[, "Genus"]))
taxonomy[, "Species"] <- as.factor(gsub("_1", "", taxonomy[, "Species"]))

table(taxonomy$High_level_taxonomy)
table(taxonomy$Kingdom)

untarget <- c("Ichthyosporia","Metazoa", "Protista", 
              "Rhizaria", "Viridiplantae", "Stramenopila")

apply(taxonomy, 2, function(x) which(x %in% untarget))

###ISSUES HERE ###
taxonomy_bad <- 
  subset(
    taxonomy,
    High_level_taxonomy%in%untarget |
      Kingdom%in%untarget  &
      HL_hit_query_cover<60 |
      HL_hit_percent_id<60  )

dim(taxonomy_bad)
head(taxonomy_bad)

taxonomy_filt <-
  taxonomy[!(rownames(taxonomy) %in% rownames(taxonomy_bad)), ]
dim(taxonomy_filt)
### Proceeding WITHOUT the above. come back if needed! ###

#Phyloseq object
physeq<- phyloseq(otu_table(otutable, taxa_are_rows = TRUE), sample_data(mapping), tax_table(as.matrix(taxonomy)), sequences)

physeq
head(sample_data(physeq))
tax_table(physeq)[tax_table(physeq)==""]<- NA
head(tax_table(physeq))

#decontamination, revist? Possible error?
count(as.data.frame(as.matrix(sample_data(
  physeq))), Description, BarcodePlate)

table(physeq@sam_data$BarcodePlate)
# adding column for controls
physeq@sam_data$is.neg <- physeq@sam_data$Description == "Control"
head(sample_data(physeq))

contaminants <-
  isContaminant(
    physeq,
    method = "prevalence",
    batch = "BarcodePlate",
    neg = "is.neg",
    threshold = 0.4)

table(contaminants$contaminant) 


# remove contaminants ---------------------------------------------------
remove_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
  
physeq_filt <-
    remove_taxa(physeq, rownames(subset(
      contaminants, contaminant%in%c("TRUE"))))
  
physeq_filt
  
  
  # plot depth --------------------------------------------------
physeq_filt@sam_data$LibrarySize <- as.numeric(as.character(sample_sums(physeq_filt)))
physeq_filt@sam_data$Index <- as.numeric(as.character(seq(nrow(physeq_filt@sam_data))))
head(physeq_filt@sam_data)
  
  
  # function to plot depth
PlotDepth <- function(physeq){
df <-
      as(sample_data(physeq), "matrix")
df <- 
      as.data.frame(df)
    # reconvert to numeric
df$LibrarySize <- as.numeric(as.character(df$LibrarySize))
df$Index <- as.numeric(as.character(df$Index))
    # order
df <- df[order(df$LibrarySize), ]
df$Index <- seq(nrow(df))
    # inspect
ggplot(data=df, aes(x=Index, y=LibrarySize, color=is.neg)) +
      geom_point(alpha =0.7) -> plot_dist
    return(plot_dist)  
  }

PlotDepth(physeq_filt) +
    labs(title="OTU", 
         subtitle = "Samples read depth", 
         x="Sample index", 
         y="Read number")
  # remove control samp[les  
physeq_clean <-
    subset_samples(physeq_filt, is.neg%in%c("FALSE"))
  
otu_table(physeq_clean) <-
    otu_table(physeq_clean)[which(rowSums(otu_table(physeq_clean)) > 0), ]
  
physeq_clean
  # rarefying dataset at even depth -----------------------------------------------------------------------
set.seed(2018)
sort(colSums(otu_table(physeq_clean)))
  
physeq_rare = rarefy_even_depth(physeq_clean, sample.size = 8000, 
                                  rngseed = FALSE, replace = TRUE, 
                                  trimOTUs = TRUE, verbose = TRUE)
  
otu_table(physeq_rare) <- otu_table(physeq_rare)[which(rowSums(otu_table(physeq_rare)) >= 1),]
physeq_rare
  
colSums(otu_table(physeq_rare))
any(taxa_sums(physeq_rare) == 0)
  
physeq_rare
  
  
  # Diversity analysis ------------------------------------------
  
  # Barplot -----------------------------------------------------
ExtractBar <- function(physeq, last_gen, Rank){
    print(head(tax_table(physeq)))
    top30 <- 
    names(sort(taxa_sums(physeq), TRUE)[1:last_gen]) # get 30 top genera
 physeq_gen <- 
    prune_taxa(top30, physeq)
    print(head(physeq_gen@sam_data))
    df_bar <- 
    physeq_gen %>%
    tax_glom(taxrank = Rank) %>%                     
    transform_sample_counts(function(x) {x/sum(x)} ) %>% 
    psmelt() %>%                                         
    filter(Abundance > 0.04) %>%                         
     arrange(get(Rank))                                     
    print(levels(df_bar[,Rank]))
    return(df_bar)
  }
  
PlotBar <- function(df, Rank){
  barplot <-
  ggplot(df, aes(x = Description,
                     y = Abundance, 
                     fill = get(Rank))) + 
  geom_bar(stat= "identity", position = "stack") +
  guides(fill = guide_legend(ncol = 1, title = "Taxa") )
  return(barplot)
}
  
head(physeq_rare@sam_data)
physeq_rare@sam_data$Description
PlotBar(ExtractBar(physeq_rare, 30, "Genus"), "Genus") 
  
#Getting alpha diversity
plot_richness(physeq_rare, x="Description")
