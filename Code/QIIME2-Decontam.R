####################
# Authors: Ashleigh Fleming Marshall
# ashleigh.marshall@ioz.ac.uk | ashleigh.marshall.16@ucl.ac.uk
####################

### MICROBIOME, DISEASE, AND HATCHING FAILURE

##Clear workspace
rm(list=ls())

#### INSTALL AND LOAD PACKAGES #### DO NOT UPDATE ANY OF THESE PACKAGES!!!
library("dada2"); packageVersion("dada2") #‘1.24.0’
library("phyloseq"); packageVersion("phyloseq") #‘1.40.0’
library("vegan"); packageVersion("vegan") #‘2.6.4’
library("Biostrings"); packageVersion("Biostrings") #‘2.64.1’
library("ggplot2"); packageVersion("ggplot2") #‘3.4.0’
library("decontam"); packageVersion("decontam") #'1.16.0'
library("qiime2R"); packageVersion("qiime2R") #‘0.99.6’
library("tidyverse"); packageVersion("tidyverse") #'1.3.2'
library("biomformat"); packageVersion("biomformat") #'1.24.0'
library("DESeq2"); packageVersion("DESeq2") #‘1.36.0'
library("dendextend"); packageVersion("dendextend") #‘1.16.0’
library("ggview"); packageVersion("ggview") #‘0.1.0’
#library("viridis"); packageVersion("viridis") #‘0.6.2’

setwd("Data/DataAnalyses")

#Adapted from https://github.com/yanxianl/Li_AqFl1-Microbiota_2021/blob/master/code/03_filtering.Rmd

#Import metadata
mtd <- read_tsv(("sample-metadata-new-noNA-removemock_forR_wEnvironmentalData_2.tsv"), comment = "#q2")
colnames(mtd)

#Import feature table
table_250.256_silva_138.1.filtered <- read_qza("250.256_silva_138.1_table-withphyla-no-unassigned-no-mitochondria-no-chloroplast-no-archaea-no-eukaryota.qza")

#Convert to count tables
count_tab_250.256_silva_138.1.filtered <- table_250.256_silva_138.1.filtered$data %>% as.data.frame() 

#Import taxonomies
taxonomy_silva_138.1_250.256 <- read_qza("silva_138.1_taxonomy_250.256.qza")

#Convert format to work with decontam
tax_tab_silva_138.1_250.256 <- taxonomy_silva_138.1_250.256$data %>% 
  as.data.frame() %>%
  mutate(Taxon = gsub("D_0", "k", Taxon), Taxon = gsub("D_1", "p", Taxon),
         Taxon = gsub("D_2", "c", Taxon), Taxon = gsub("D_3", "o", Taxon),
         Taxon = gsub("D_4", "f", Taxon), Taxon = gsub("D_5", "g", Taxon),
         Taxon = gsub("D_6", "s", Taxon)) %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  select(-Confidence) #Warning message: Expected 7 pieces. Missing pieces filled with `NA` in 4525 rows [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, ...] = this matches the number of NAs in the species column

#Remove likely contaminants based on prevalence in control samples
#Applying it to the data filtered to remove non-target length sequences (filtered to 250-256 nts) and filtered to remove anything not classified to phyla level, anything unassigned at kingdom level, and any mitochondria, chloroplasts, archae, and eukaryota
colnames(count_tab_250.256_silva_138.1.filtered)
negative_controls_250.256_silva_138.1 <- count_tab_250.256_silva_138.1.filtered[,c(237:239)]
vector_for_decontam_250.256_silva_138.1 <- c(rep(FALSE,236),rep(TRUE,3),rep(FALSE,16))
contam_df_250.256_silva_138.1 <- isContaminant(t(count_tab_250.256_silva_138.1.filtered), neg=vector_for_decontam_250.256_silva_138.1, method="prevalence") 
table(contam_df_250.256_silva_138.1$contaminant) #Identified 6 contaminants
contam_df_250.256_silva_138.1_0.5 <- isContaminant(t(count_tab_250.256_silva_138.1.filtered), neg=vector_for_decontam_250.256_silva_138.1, method="prevalence",threshold=0.5) 
table(contam_df_250.256_silva_138.1_0.5$contaminant) #Identified 18 contaminants
contam_df_250.256_silva_138.1_0.2 <- isContaminant(t(count_tab_250.256_silva_138.1.filtered), neg=vector_for_decontam_250.256_silva_138.1, method="prevalence",threshold=0.2) 
table(contam_df_250.256_silva_138.1_0.2$contaminant) #Identified 9 contaminants
contam_df_250.256_silva_138.1_0.3 <- isContaminant(t(count_tab_250.256_silva_138.1.filtered), neg=vector_for_decontam_250.256_silva_138.1, method="prevalence",threshold=0.3) 
table(contam_df_250.256_silva_138.1_0.3$contaminant) #Identified 12 contaminants
contam_df_250.256_silva_138.1_0.4 <- isContaminant(t(count_tab_250.256_silva_138.1.filtered), neg=vector_for_decontam_250.256_silva_138.1, method="prevalence",threshold=0.4) 
table(contam_df_250.256_silva_138.1_0.4$contaminant) #Identified 15 contaminants
contam_asvs_250.256_silva_138.1 <- row.names(contam_df_250.256_silva_138.1[contam_df_250.256_silva_138.1$contaminant == TRUE, ])
contaminants_250.256_silva_138.1 <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_250.256_silva_138.1, ] #Same as identified in DADA2 pipeline
contam_asvs_250.256_silva_138.1_0.5 <- row.names(contam_df_250.256_silva_138.1_0.5[contam_df_250.256_silva_138.1_0.5$contaminant == TRUE, ])
contaminants_250.256_silva_138.1_0.5 <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_250.256_silva_138.1_0.5, ] #Same as identified in DADA2 pipeline
contam_asvs_250.256_silva_138.1_0.2 <- row.names(contam_df_250.256_silva_138.1_0.2[contam_df_250.256_silva_138.1_0.2$contaminant == TRUE, ])
contaminants_250.256_silva_138.1_0.2 <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_250.256_silva_138.1_0.2, ] #Same as identified in DADA2 pipeline
contam_asvs_250.256_silva_138.1_0.3 <- row.names(contam_df_250.256_silva_138.1_0.3[contam_df_250.256_silva_138.1_0.3$contaminant == TRUE, ])
contaminants_250.256_silva_138.1_0.3 <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_250.256_silva_138.1_0.3, ] #Same as identified in DADA2 pipeline with one extra
contam_asvs_250.256_silva_138.1_0.4 <- row.names(contam_df_250.256_silva_138.1_0.4[contam_df_250.256_silva_138.1_0.4$contaminant == TRUE, ])
contaminants_250.256_silva_138.1_0.4 <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_250.256_silva_138.1_0.4, ] #Same as identified in DADA2 pipeline
#write.csv(contaminants_250.256_silva_138.1_0.4,"DecontamOutputs/contaminants_250.256_silva_138.1_0.4.csv")

# Remove the contaminants - cut-off at 0.4 (see below)
# Making new count table
count_tab_250.256_silva_138.1.filtered_no_contam <- count_tab_250.256_silva_138.1.filtered[!row.names(count_tab_250.256_silva_138.1.filtered) %in% contam_asvs_250.256_silva_138.1_0.4, ] #250-256 nts, non-target-taxa filtered using silva 138.1
dim(count_tab_250.256_silva_138.1.filtered_no_contam) #6556 255
#saveRDS(count_tab_250.256_silva_138.1.filtered_no_contam, file = "count_tab_250.256_silva_138.1.filtered_no_contam.RDS") 

# Making new taxonomy table
tax_tab_silva_138.1_250.256_no_contam <- tax_tab_silva_138.1_250.256[!row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_250.256_silva_138.1_0.4, ]
dim(tax_tab_silva_138.1_250.256_no_contam) #6954 7

###

# Trying with phyloseq method instead (see https://benjjneb.github.io/decontam/vignettes/decontam_intro.html#putting-it-all-together)

# Make phyloseq object with data with non-target-length and non-target taxa removed (with Silva 138.1)
ps_138.1 <- phyloseq(sample_data(column_to_rownames(mtd, "SampleID")),
               otu_table(as.matrix(count_tab_250.256_silva_138.1.filtered), taxa_are_rows = TRUE),
               tax_table(as.matrix(tax_tab_silva_138.1_250.256))) #6571 taxa and 255 samples
ps_138.1 #6571 taxa and 255 samples
head(sample_data(ps_138.1))

df_138.1 <- as.data.frame(sample_data(ps_138.1)) # Put sample_data into a ggplot-friendly data.frame
df_138.1$LibrarySize <- sample_sums(ps_138.1)
df_138.1 <- df_138.1[order(df_138.1$LibrarySize),]
df_138.1$Index <- seq(nrow(df_138.1))
ggplot(data=df_138.1, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
class(df_138.1)

# Remove mock communities
rownames(sample_data(ps_138.1))
sample_or_neg <- c("Sample","Negative Control")
ps_138.1_nomock <- subset_samples(ps_138.1, Sample %in% sample_or_neg)
rownames(sample_data(ps_138.1_nomock)) #Successfully removed mock communities

sample_data(ps_138.1_nomock)$is.neg <- sample_data(ps_138.1_nomock)$Sample_or_Control == "Control Sample"
contamdf.prev_138.1_nomock <- isContaminant(ps_138.1_nomock, method="prevalence", neg="is.neg")
table(contamdf.prev_138.1_nomock$contaminant) #Identified 5 potential contaminants
which(contamdf.prev_138.1_nomock$contaminant) #95  655 2803 4253 4738
contamdf.prev02_138.1_nomock <- isContaminant(ps_138.1_nomock, method="prevalence", neg="is.neg", threshold=0.2)
table(contamdf.prev02_138.1_nomock$contaminant) #Identified 9 potential contaminants
contamdf.prev03_138.1_nomock <- isContaminant(ps_138.1_nomock, method="prevalence", neg="is.neg", threshold=0.3)
table(contamdf.prev03_138.1_nomock$contaminant) #Identified 11 potential contaminants
contamdf.prev04_138.1_nomock <- isContaminant(ps_138.1_nomock, method="prevalence", neg="is.neg", threshold=0.4)
table(contamdf.prev04_138.1_nomock$contaminant) #Identified 15 potential contaminants
contamdf.prev05_138.1_nomock <- isContaminant(ps_138.1_nomock, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05_138.1_nomock$contaminant) #Identified 18 potential contaminants
which(contamdf.prev05_138.1_nomock$contaminant) #95  305  322  655 2191 2302 2803 2875 3003 3149 3188 3199 4094 4253 4738 4920 5353 6334
contam_asvs_0.1_ps_138.1_nomock <- row.names(contamdf.prev_138.1_nomock[contamdf.prev_138.1_nomock$contaminant == TRUE, ])
contaminants_0.1_ps_138.1_nomock <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.1_ps_138.1_nomock, ] #Same as identified in other methods
contam_asvs_0.2_ps_138.1_nomock <- row.names(contamdf.prev02_138.1_nomock[contamdf.prev02_138.1_nomock$contaminant == TRUE, ])
contaminants_0.2_ps_138.1_nomock <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.2_ps_138.1_nomock, ] #Same as identified in other methods
contam_asvs_0.3_ps_138.1_nomock <- row.names(contamdf.prev03_138.1_nomock[contamdf.prev03_138.1_nomock$contaminant == TRUE, ])
contaminants_0.3_ps_138.1_nomock <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.3_ps_138.1_nomock, ] #Same as identified in other methods
contam_asvs_0.4_ps_138.1_nomock <- row.names(contamdf.prev04_138.1_nomock[contamdf.prev04_138.1_nomock$contaminant == TRUE, ])
contaminants_0.4_ps_138.1_nomock <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.4_ps_138.1_nomock, ] #Same as identified in other methods
contam_asvs_0.5_ps_138.1_nomock <- row.names(contamdf.prev05_138.1_nomock[contamdf.prev05_138.1_nomock$contaminant == TRUE, ])
contaminants_0.5_ps_138.1_nomock <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.5_ps_138.1_nomock, ] #Same as identified in other methods
#write.csv(contaminants_0.5_ps_138.1_nomock,"DecontamOutputs/contaminants_0.5_ps_138.1_nomock.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
ps_138.1.pa_nomock <- transform_sample_counts(ps_138.1_nomock, function(abund) 1*(abund>0))
ps_138.1.pa.neg_nomock <- prune_samples(sample_data(ps_138.1.pa_nomock)$Sample_or_Control == "Control Sample", ps_138.1.pa_nomock)
ps_138.1.pa.pos_nomock <- prune_samples(sample_data(ps_138.1.pa_nomock)$Sample_or_Control == "True Sample", ps_138.1.pa_nomock)
# Make data.frame of prevalence in positive and negative samples
df.pa_138.1_nomock <- data.frame(pa.pos=taxa_sums(ps_138.1.pa.pos_nomock), pa.neg=taxa_sums(ps_138.1.pa.neg_nomock),
                               contaminant=contamdf.prev05_138.1_nomock$contaminant)
ggplot(data=df.pa_138.1_nomock, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Check with only kaki and only hihi (mocks removed)
colnames(sample_data(ps_138.1_nomock))
sample_data(ps_138.1_nomock)$species
ps_138.1_nomock_kaki <- subset_samples(ps_138.1_nomock, !species == "Hihi")
rownames(sample_data(ps_138.1_nomock_kaki)) 
ps_138.1_nomock_hihi <- subset_samples(ps_138.1_nomock, !species == "Kaki")
rownames(sample_data(ps_138.1_nomock_hihi)) 

sample_data(ps_138.1_nomock_kaki)$is.neg <- sample_data(ps_138.1_nomock_kaki)$Sample_or_Control == "Control Sample"
contamdf.prev_138.1_kaki <- isContaminant(ps_138.1_nomock_kaki, method="prevalence", neg="is.neg")
table(contamdf.prev_138.1_kaki$contaminant) #Identified 6 potential contaminants
which(contamdf.prev_138.1_kaki$contaminant) #655 3003 3149 3188 4253 4738
contamdf.prev02_138.1_kaki <- isContaminant(ps_138.1_nomock_kaki, method="prevalence", neg="is.neg", threshold=0.2)
table(contamdf.prev02_138.1_kaki$contaminant) #Identified 13 potential contaminants
contamdf.prev03_138.1_kaki <- isContaminant(ps_138.1_nomock_kaki, method="prevalence", neg="is.neg", threshold=0.3)
table(contamdf.prev03_138.1_kaki$contaminant) #Identified 14 potential contaminants
contamdf.prev04_138.1_kaki <- isContaminant(ps_138.1_nomock_kaki, method="prevalence", neg="is.neg", threshold=0.4)
table(contamdf.prev04_138.1_kaki$contaminant) #Identified 15 potential contaminants
contamdf.prev05_138.1_kaki <- isContaminant(ps_138.1_nomock_kaki, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05_138.1_kaki$contaminant) #Identified 18 potential contaminants
which(contamdf.prev05_138.1_kaki$contaminant) #95  305  322  655 1241 2191 2302 2875 3003 3149 3188 3199 4094 4253 4738 4920 5353 6334
contam_asvs_0.1_ps_138.1_kaki <- row.names(contamdf.prev_138.1_kaki[contamdf.prev_138.1_kaki$contaminant == TRUE, ])
contaminants_0.1_ps_138.1_kaki <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.1_ps_138.1_kaki, ] 
contam_asvs_0.2_ps_138.1_kaki <- row.names(contamdf.prev02_138.1_kaki[contamdf.prev02_138.1_kaki$contaminant == TRUE, ])
contaminants_0.2_ps_138.1_kaki <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.2_ps_138.1_kaki, ] 
contam_asvs_0.3_ps_138.1_kaki <- row.names(contamdf.prev03_138.1_kaki[contamdf.prev03_138.1_kaki$contaminant == TRUE, ])
contaminants_0.3_ps_138.1_kaki <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.3_ps_138.1_kaki, ] 
contam_asvs_0.4_ps_138.1_kaki <- row.names(contamdf.prev04_138.1_kaki[contamdf.prev04_138.1_kaki$contaminant == TRUE, ])
contaminants_0.4_ps_138.1_kaki <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.4_ps_138.1_kaki, ] 
contam_asvs_0.5_ps_138.1_kaki <- row.names(contamdf.prev05_138.1_kaki[contamdf.prev05_138.1_kaki$contaminant == TRUE, ])
contaminants_0.5_ps_138.1_kaki <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.5_ps_138.1_kaki, ] 
#write.csv(contaminants_0.5_ps_138.1_kaki,"DecontamOutputs/contaminants_0.5_ps_138.1_kaki.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
ps_138.1.pa_kaki <- transform_sample_counts(ps_138.1_nomock_kaki, function(abund) 1*(abund>0))
ps_138.1.pa.neg_kaki <- prune_samples(sample_data(ps_138.1.pa_kaki)$Sample_or_Control == "Control Sample", ps_138.1.pa_kaki)
ps_138.1.pa.pos_kaki <- prune_samples(sample_data(ps_138.1.pa_kaki)$Sample_or_Control == "True Sample", ps_138.1.pa_kaki)
# Make data.frame of prevalence in positive and negative samples
df.pa_138.1_kaki <- data.frame(pa.pos=taxa_sums(ps_138.1.pa.pos_kaki), pa.neg=taxa_sums(ps_138.1.pa.neg_kaki),
                          contaminant=contamdf.prev_138.1_kaki$contaminant)
ggplot(data=df.pa_138.1_kaki, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

sample_data(ps_138.1_nomock_hihi)$is.neg <- sample_data(ps_138.1_nomock_hihi)$Sample_or_Control == "Control Sample"
contamdf.prev_138.1_hihi <- isContaminant(ps_138.1_nomock_hihi, method="prevalence", neg="is.neg")
table(contamdf.prev_138.1_hihi$contaminant) #Identified 4 potential contaminants
which(contamdf.prev_138.1_hihi$contaminant) #95 2803 4253 4738
contamdf.prev02_138.1_hihi <- isContaminant(ps_138.1_nomock_hihi, method="prevalence", neg="is.neg", threshold=0.2)
table(contamdf.prev02_138.1_hihi$contaminant) #Identified 8 potential contaminants
contamdf.prev03_138.1_hihi <- isContaminant(ps_138.1_nomock_hihi, method="prevalence", neg="is.neg", threshold=0.3)
table(contamdf.prev03_138.1_hihi$contaminant) #Identified 9 potential contaminants
contamdf.prev04_138.1_hihi <- isContaminant(ps_138.1_nomock_hihi, method="prevalence", neg="is.neg", threshold=0.4)
table(contamdf.prev04_138.1_hihi$contaminant) #Identified 12 potential contaminants
contamdf.prev05_138.1_hihi <- isContaminant(ps_138.1_nomock_hihi, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05_138.1_hihi$contaminant) #Identified 16 potential contaminants
which(contamdf.prev05_138.1_hihi$contaminant) #95  322  655 2191 2302 2803 3003 3149 3188 3199 4094 4253 4738 4920 5353 6334
contam_asvs_0.1_ps_138.1_hihi <- row.names(contamdf.prev_138.1_hihi[contamdf.prev_138.1_hihi$contaminant == TRUE, ])
contaminants_0.1_ps_138.1_hihi <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.1_ps_138.1_hihi, ] 
contam_asvs_0.2_ps_138.1_hihi <- row.names(contamdf.prev02_138.1_hihi[contamdf.prev02_138.1_hihi$contaminant == TRUE, ])
contaminants_0.2_ps_138.1_hihi <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.2_ps_138.1_hihi, ] 
contam_asvs_0.3_ps_138.1_hihi <- row.names(contamdf.prev03_138.1_hihi[contamdf.prev03_138.1_hihi$contaminant == TRUE, ])
contaminants_0.3_ps_138.1_hihi <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.3_ps_138.1_hihi, ] 
contam_asvs_0.4_ps_138.1_hihi <- row.names(contamdf.prev04_138.1_hihi[contamdf.prev04_138.1_hihi$contaminant == TRUE, ])
contaminants_0.4_ps_138.1_hihi <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.4_ps_138.1_hihi, ] 
contam_asvs_0.5_ps_138.1_hihi <- row.names(contamdf.prev05_138.1_hihi[contamdf.prev05_138.1_hihi$contaminant == TRUE, ])
contaminants_0.5_ps_138.1_hihi <- tax_tab_silva_138.1_250.256[row.names(tax_tab_silva_138.1_250.256) %in% contam_asvs_0.5_ps_138.1_hihi, ] 
#write.csv(contaminants_0.5_ps_138.1_hihi,"DecontamOutputs/contaminants_0.5_ps_138.1_hihi.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
ps_138.1.pa_hihi <- transform_sample_counts(ps_138.1_nomock_hihi, function(abund) 1*(abund>0))
ps_138.1.pa.neg_hihi <- prune_samples(sample_data(ps_138.1.pa_hihi)$Sample_or_Control == "Control Sample", ps_138.1.pa_hihi)
ps_138.1.pa.pos_hihi <- prune_samples(sample_data(ps_138.1.pa_hihi)$Sample_or_Control == "True Sample", ps_138.1.pa_hihi)
# Make data.frame of prevalence in positive and negative samples
df.pa_138.1_hihi <- data.frame(pa.pos=taxa_sums(ps_138.1.pa.pos_hihi), pa.neg=taxa_sums(ps_138.1.pa.neg_hihi),
                               contaminant=contamdf.prev_138.1_hihi$contaminant)
ggplot(data=df.pa_138.1_hihi, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

## Remove contaminants

# Remove 15 contaminants identified with probability threshold of 0.4 in all samples. Going with 0.4 cut-off as majority seem likely to be contaminants based on biological likelihood and previous identification as contaminants in other studies, plus doesn't exclude any ASVs not identified as contaminants in each species separately
# Same 15 found with and without mock communities, but want to keep mock communities for now use the data with mocks.
ps_138.1 #6571 taxa and 255 samples
ps_138.1_no_contam_0.4 <- prune_taxa(!contamdf.prev04_138.1_nomock$contaminant, ps_138.1)
ps_138.1_no_contam_0.4 #6556 taxa and 255 samples = successfully removed the 15 contaminants

# Export the feature table without the contaminants to be used again in QIIME2
otu_table(ps_138.1_no_contam_0.4) %>%
  as("matrix") %>%
  make_biom() %>%
  write_biom("filtered-silva-138.1-table-decontam_0.4.biom")
