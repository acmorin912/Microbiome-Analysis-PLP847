#Ashlynn Morin, HPCC R, Fungal Roots Analysis
#New code for 2021 data
#to get here, you will need to set your working directory as one of mine for all the data.
#PATHWAY FOR WORKING DIRECTORY: /mnt/research/bonito_lab/Morin/fungal_roots_2021/Needed data now
#hopefully accessible? If not I will add to my github as well :)
#If you work on it, please save it with the date, time, and your initals (Also make a note where you started working versus ended working)
#I am following the following resource from Reid Longley: https://github.com/longleyr/Management-of-Soybean-Code-and-Files/blob/master/Updated%20Network%20Example.R"

#Get Phyloseq (since install ____ isnt working)
source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)

#Getting Biostrings (since install ______ isn't working)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

#Getting "SpiecEasi"
install.packages(devtools)
library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)

#make sure packages there
library("Biostrings")
library("phyloseq")
library("ggplot2")
library("SpiecEasi")
library("igraph")
library("vegan")
library("stringi")
library("rhdf5")
library("zlibbioc")
library("S4Vectors")
library("yaml")
library("colorspace")
library("indicspecies")

#Time to get to work
#Create phyloseq ITS_OTU table
ITS_otus<- read.delim("otu_table_ITS_UPARSE_R1.txt",
                      row.names=1) 
head(ITS_otus)

#ITS_OTU_PHY
ITS_otus_phy <-otu_table(ITS_otus,
                         taxa_are_rows = TRUE)
ITS_otus_phy

#ITS_Metadata (WILL NEED TO UPDATE THIS NAME!!!!)
ITS_metadata <-read.delim("root_fun_map_soybean_2021.txt",
                          row.names=1)
ITS_metadata
ITS_metadata_phy <-sample_data(ITS_metadata)

#ITS taxonomy
  ITS_taxonomy<- read.delim("otu_taxonomy.blast",
                            header= TRUE,
                            row.names =NULL,
                            sep= ",")
ITS_taxonomy
ITS_taxonomy2<-ITS_taxonomy
ITS_taxonomy2
ITS_taxonomy_phy <- tax_table(as.matrix(ITS_taxonomy))
ITS_taxonomy_phy

#Sequences
ITS_sequences<- readDNAStringSet("otus_R1.fasta", format="fasta", seek.first.rec = TRUE, use.names = TRUE)
ITS_sequences

phyloseq_object_Fungi <- phyloseq(ITS_otus_phy, ITS_metadata_phy, ITS_taxonomy_phy, ITS_sequences)
#^^ISSUES HERE
#11/17/22 11:12am, stopped working AM