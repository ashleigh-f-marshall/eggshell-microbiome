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
library("fantaxtic"); packageVersion("fantaxtic") #‘0.2.0’

setwd("Data/DataAnalyses")

#Adapted from https://github.com/yanxianl/Li_AqFl1-Microbiota_2021/blob/master/code/03_filtering.Rmd

#Import metadata
mtd <- read_tsv(("sample-metadata-new-noNA-removemock_forR_wEnvironmentalData_2.tsv"), comment = "#q2")

#Import feature table and convert to count table
table_250.256_silva_138.1.filtered_decontam_noMCnoNTC <- read_qza("filtered-silva-138.1-table-decontam_0.4_noMCnoNTC.qza")
count_tab_250.256_silva_138.1.filtered_decontam_noMCnoNTC <- table_250.256_silva_138.1.filtered_decontam_noMCnoNTC$data %>% as.data.frame()
dim(count_tab_250.256_silva_138.1.filtered_decontam_noMCnoNTC) #6528 249

#Import taxonomies
taxonomy_silva_138.1_250.256_filtered_decontam_noMCnoNTC <- read_qza("silva_138.1_taxonomy_250.256-filtered-decontam_noMCnoNTC.qza")

tax_tab_silva_138.1_250.256_filtered_decontam_noMCnoNTC <- taxonomy_silva_138.1_250.256_filtered_decontam_noMCnoNTC$data %>% 
  as.data.frame() %>%
  mutate(Taxon = gsub("D_0", "k", Taxon), Taxon = gsub("D_1", "p", Taxon),
         Taxon = gsub("D_2", "c", Taxon), Taxon = gsub("D_3", "o", Taxon),
         Taxon = gsub("D_4", "f", Taxon), Taxon = gsub("D_5", "g", Taxon),
         Taxon = gsub("D_6", "s", Taxon)) %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  select(-Confidence) #Warning message: Expected 7 pieces. Missing pieces filled with `NA` in 4525 rows [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, ...] = this matches the number of NAs in the species column
dim(tax_tab_silva_138.1_250.256_filtered_decontam_noMCnoNTC) #6528 7

test_physeq <- qza_to_phyloseq(features="filtered-silva-138.1-table-decontam_0.4_noMCnoNTC.qza",
                               tree="rooted-tree_allfilter_decontam_noMCnoNTC.qza",
                               taxonomy="silva_138.1_taxonomy_250.256-filtered-decontam_noMCnoNTC.qza")
test_physeq_sd <- phyloseq(sample_data(column_to_rownames(mtd, "SampleID")))
physeq_tree_noMCnoNTC <- merge_phyloseq(test_physeq,test_physeq_sd)
physeq_tree_noMCnoNTC #6528 taxa and 249 samples

#### DIFFERENTIAL ABUNDANCE TESTING ####

physeq_tree_noMCnoNTC_kaki <- subset_samples(physeq_tree_noMCnoNTC, !species == "Hihi")
rownames(sample_data(physeq_tree_noMCnoNTC_kaki)) 
physeq_tree_noMCnoNTC_hihi <- subset_samples(physeq_tree_noMCnoNTC, !species == "Kaki")
rownames(sample_data(physeq_tree_noMCnoNTC_hihi))

## DESeq2 (https://astrobiomike.github.io/amplicon/dada2_workflow_ex)
all_deseq <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC, ~species) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
all_deseq <- DESeq(all_deseq)
deseq_res_hihi_vs_kaki <- results(all_deseq, alpha=0.05, contrast=c("species", "Hihi", "Kaki"))
summary(deseq_res_hihi_vs_kaki) #Out of 3207 ASVs with nonzero total read count and adj-p < 0.5 there are 6 (0.19%) increased (i.e. greater proportion in hihi than kaki) when comparing hihi to kaki, and 29 (0.9%) decreased (i.e. at a lower count abundance in the hihi than in the kaki samples) = 35 significant in total
#With adj-p <0.01 4 increased and 27 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_deseq_res_hihi_vs_kaki <- deseq_res_hihi_vs_kaki[which(deseq_res_hihi_vs_kaki$padj < 0.05), ]
summary(sigtab_res_deseq_res_hihi_vs_kaki) #Table only contains ASVs considered significantly differentially abundant

count_tab_clean_phy <- otu_table(count_tab_250.256_silva_138.1.filtered_decontam_noMCnoNTC, taxa_are_rows=T)
taxonomy_silva_138.1_250.256_noMCnoNTC <- read_qza("silva_138.1_taxonomy_250.256-filtered-decontam_noMCnoNTC.qza")
tax_tab <- taxonomy_silva_138.1_250.256_noMCnoNTC$data %>% 
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  select(-Confidence) #Warning message: Expected 7 pieces. Missing pieces filled with `NA` in 4252 rows [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, ...].
dim(tax_tab) #6528 7
tax_tab <- as.matrix(tax_tab)
tax_tab_phy <- tax_table(tax_tab)

ASV_physeq <- phyloseq(count_tab_clean_phy, tax_tab_phy, mtd)
tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))] <- gsub(tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))], pattern = "d__", replacement = "")
tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))] <- gsub(tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))], pattern = " p__", replacement = "")
tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))] <- gsub(tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))], pattern = " c__", replacement = "")
tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))] <- gsub(tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))], pattern = " o__", replacement = "")
tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))] <- gsub(tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))], pattern = " f__", replacement = "")
tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))] <- gsub(tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))], pattern = " g__", replacement = "")
tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))] <- gsub(tax_table(ASV_physeq)[, colnames(tax_table(ASV_physeq))], pattern = " s__", replacement = "")

head(tax_table(ASV_physeq))
ASV_physeq_rename <- name_na_taxa(ASV_physeq, na_label = "Unidentified <tax> (<rank>)")
head(tax_table(ASV_physeq_rename))
ASV_physeq_rename_genera <- ASV_physeq_rename
ASV_physeq_rename_families <- ASV_physeq_rename
ASV_physeq_rename_phyla <- ASV_physeq_rename
tax_table(ASV_physeq_rename_genera) <- tax_table(ASV_physeq_rename)[,1:6]
tax_table(ASV_physeq_rename_families) <- tax_table(ASV_physeq_rename)[,1:5]
tax_table(ASV_physeq_rename_phyla) <- tax_table(ASV_physeq_rename)[,1:2]

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_deseq_res_hihi_vs_kaki_with_tax <- cbind(as(sigtab_res_deseq_res_hihi_vs_kaki, "data.frame"), as(tax_table(ASV_physeq_rename)[row.names(sigtab_res_deseq_res_hihi_vs_kaki), ], "matrix"))
dim(sigtab_res_deseq_res_hihi_vs_kaki_with_tax) #35 13

#Sort by the log2FoldChange column
sigtab_res_deseq_res_hihi_vs_kaki_with_tax[order(sigtab_res_deseq_res_hihi_vs_kaki_with_tax$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the 35 significant ASVs = negative log2FoldChange means lower count abundance in hihi vs. kaki, positive log2FoldChange means greater count abundance in hihi vs. kaki
#write.csv(sigtab_res_deseq_res_hihi_vs_kaki_with_tax,"DifferentialAbundance/DESeq2_Results/deseq_results_hihi.vs.kaki_qiime2.csv")

dim(table(tax_table(physeq_tree_noMCnoNTC)[,1])) #1 kingdom
dim(table(tax_table(physeq_tree_noMCnoNTC)[,2])) #35 phyla = matches collapsed table from QIIME2
dim(table(tax_table(physeq_tree_noMCnoNTC)[,3])) #86 classes
dim(table(tax_table(physeq_tree_noMCnoNTC)[,4])) #214 orders
dim(table(tax_table(physeq_tree_noMCnoNTC)[,5])) #379 families
dim(table(tax_table(physeq_tree_noMCnoNTC)[,6])) #807 genera = matches collapsed table from QIIME2 except no g__ and doesn't separate if have differences at higher classifications (i.e. all g__uncultured counted together)
dim(table(tax_table(physeq_tree_noMCnoNTC)[,7])) #592 species
genera_list_fromR <- table(tax_table(physeq_tree_noMCnoNTC)[,6])
#write.csv(genera_list_fromR,"genera_list_fromR.csv")
as.character(unique(tax_table(physeq_tree_noMCnoNTC)[,6])) #808 (includes NA)

physeq_tree_noMCnoNTC_phylum <- tax_glom(physeq_tree_noMCnoNTC,"Phylum")
physeq_tree_noMCnoNTC_phylum #35 taxa and 249 samples
physeq_tree_noMCnoNTC_phylum_NArmFalse <- tax_glom(physeq_tree_noMCnoNTC,"Phylum",NArm=FALSE)
physeq_tree_noMCnoNTC_phylum_NArmFalse #35 taxa and 249 samples
tax_table(physeq_tree_noMCnoNTC_phylum) #NAs used to replace Class, Order, Family, Genus, Species
physeq_tree_noMCnoNTC_family <- tax_glom(physeq_tree_noMCnoNTC,"Family")
physeq_tree_noMCnoNTC_family #398 taxa and 249 samples
physeq_tree_noMCnoNTC_family_NArmFalse <- tax_glom(physeq_tree_noMCnoNTC,"Family",NArm=FALSE)
physeq_tree_noMCnoNTC_family_NArmFalse #424 taxa and 249 samples
tax_table(physeq_tree_noMCnoNTC_family) #NAs used to replace Genus, Species
physeq_tree_noMCnoNTC_genus <- tax_glom(physeq_tree_noMCnoNTC,"Genus")
physeq_tree_noMCnoNTC_genus #910 taxa and 249 samples = leaves out any taxa not identified at genus level
physeq_tree_noMCnoNTC_genus_NArmFalse <- tax_glom(physeq_tree_noMCnoNTC,"Genus",NArm=FALSE)
physeq_tree_noMCnoNTC_genus_NArmFalse #1013 taxa and 249 samples
tax_table(physeq_tree_noMCnoNTC_genus) #NAs used to replace Species
genera_list_tax_glom <- tax_table(physeq_tree_noMCnoNTC_genus)
genera_list_tax_glom_NArmFalse <- tax_table(physeq_tree_noMCnoNTC_genus_NArmFalse) #Anything not identified at genus is still there but labelled NA
#write.csv(genera_list_tax_glom,"genera_list_tax_glom.csv")
#write.csv(genera_list_tax_glom_NArmFalse,"genera_list_tax_glom_NArmFalse.csv")

#From https://github.com/joey711/phyloseq/issues/850 --> renaming any NAs to be named for the last identified taxa. Only need to do for genera and family because already filtered to remove anything not identified at phylum level
tax.clean <- data.frame(tax_table(physeq_tree_noMCnoNTC_genus_NArmFalse))
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}
tax_table(physeq_tree_noMCnoNTC_genus_NArmFalse) <- as.matrix(tax.clean)

tax.clean_fam <- data.frame(tax_table(physeq_tree_noMCnoNTC_family_NArmFalse))
for (i in 1:7){ tax.clean_fam[,i] <- as.character(tax.clean_fam[,i])}
tax.clean_fam[is.na(tax.clean_fam)] <- ""
for (i in 1:nrow(tax.clean_fam)){
  if (tax.clean_fam[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean_fam[i,1], sep = "")
    tax.clean_fam[i, 2:7] <- kingdom
  } else if (tax.clean_fam[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean_fam[i,2], sep = "")
    tax.clean_fam[i, 3:7] <- phylum
  } else if (tax.clean_fam[i,4] == ""){
    class <- paste("Class_", tax.clean_fam[i,3], sep = "")
    tax.clean_fam[i, 4:7] <- class
  } else if (tax.clean_fam[i,5] == ""){
    order <- paste("Order_", tax.clean_fam[i,4], sep = "")
    tax.clean_fam[i, 5:7] <- order
  } else if (tax.clean_fam[i,6] == ""){
    family <- paste("Family_", tax.clean_fam[i,5], sep = "")
    tax.clean_fam[i, 6:7] <- family
  } else if (tax.clean_fam[i,7] == ""){
    tax.clean_fam$Species[i] <- paste("Genus",tax.clean_fam$Genus[i], sep = "_")
  }
}
tax_table(physeq_tree_noMCnoNTC_family_NArmFalse) <- as.matrix(tax.clean_fam)

physeq_tree_noMCnoNTC_genus_NArmFalse_kaki <- subset_samples(physeq_tree_noMCnoNTC_genus_NArmFalse, !species == "Hihi")
rownames(sample_data(physeq_tree_noMCnoNTC_genus_NArmFalse_kaki)) 
physeq_tree_noMCnoNTC_genus_NArmFalse_kaki #1013 taxa and 92 samples
physeq_tree_noMCnoNTC_genus_NArmFalse_hihi <- subset_samples(physeq_tree_noMCnoNTC_genus_NArmFalse, !species == "Kaki")
rownames(sample_data(physeq_tree_noMCnoNTC_genus_NArmFalse_hihi)) 
physeq_tree_noMCnoNTC_genus_NArmFalse_hihi #1013 taxa and 157 samples

physeq_tree_noMCnoNTC_family_NArmFalse_kaki <- subset_samples(physeq_tree_noMCnoNTC_family_NArmFalse, !species == "Hihi")
rownames(sample_data(physeq_tree_noMCnoNTC_family_NArmFalse_kaki)) 
physeq_tree_noMCnoNTC_family_NArmFalse_kaki #424 taxa and 92 samples
physeq_tree_noMCnoNTC_family_NArmFalse_hihi <- subset_samples(physeq_tree_noMCnoNTC_family_NArmFalse, !species == "Kaki")
rownames(sample_data(physeq_tree_noMCnoNTC_family_NArmFalse_hihi)) 
physeq_tree_noMCnoNTC_family_NArmFalse_hihi #424 taxa and 157 samples

physeq_tree_noMCnoNTC_phylum_kaki <- subset_samples(physeq_tree_noMCnoNTC_phylum, !species == "Hihi")
rownames(sample_data(physeq_tree_noMCnoNTC_phylum_kaki)) 
physeq_tree_noMCnoNTC_phylum_kaki #35 taxa and 92 samples
physeq_tree_noMCnoNTC_phylum_hihi <- subset_samples(physeq_tree_noMCnoNTC_phylum, !species == "Kaki")
rownames(sample_data(physeq_tree_noMCnoNTC_phylum_hihi)) 
physeq_tree_noMCnoNTC_phylum_hihi #35 taxa and 157 samples

## DESeq2 (https://astrobiomike.github.io/amplicon/dada2_workflow_ex)
#Grouped by genera level
all_deseq_genus <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_genus_NArmFalse, ~species) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
all_deseq_genus <- DESeq(all_deseq_genus)
deseq_res_hihi_vs_kaki_genus <- results(all_deseq_genus, alpha=0.05, contrast=c("species", "Hihi", "Kaki"))
summary(deseq_res_hihi_vs_kaki_genus) #Out of 778 genera with nonzero total read count and adj-p < 0.5 there are 2 (0.26%) increased (i.e. greater proportion in hihi than kaki) when comparing hihi to kaki, and 110 (14%) decreased (i.e. at a lower count abundance in the hihi than in the kaki samples) = 112 significant in total
#With adj-p <0.01 2 (0.26%) increased and 65 decreased (8.4%)

#Subset this table to only include ones that pass p < 0.05
sigtab_res_deseq_res_hihi_vs_kaki_genus <- deseq_res_hihi_vs_kaki_genus[which(deseq_res_hihi_vs_kaki_genus$padj < 0.05), ]
summary(sigtab_res_deseq_res_hihi_vs_kaki_genus) #Table only contains genera considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_deseq_res_hihi_vs_kaki_with_tax_genus <- cbind(as(sigtab_res_deseq_res_hihi_vs_kaki_genus, "data.frame"), as(tax_table(ASV_physeq_rename_genera)[row.names(sigtab_res_deseq_res_hihi_vs_kaki_genus), ], "matrix"))
dim(sigtab_res_deseq_res_hihi_vs_kaki_with_tax_genus) #112

#Sort by the log2FoldChange column
sigtab_res_deseq_res_hihi_vs_kaki_with_tax_genus[order(sigtab_res_deseq_res_hihi_vs_kaki_with_tax_genus$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the 112 significant genera = negative log2FoldChange means lower count abundance in hihi vs. kaki, positive log2FoldChange means greater count abundance in hihi vs. kaki
#write.csv(sigtab_res_deseq_res_hihi_vs_kaki_with_tax_genus,"DifferentialAbundance/DESeq2_Results/deseq_results_hihi.vs.kaki_qiime2_genus.csv")

#Grouped by family level
all_deseq_family <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_family, ~species) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
all_deseq_family <- DESeq(all_deseq_family)
deseq_res_hihi_vs_kaki_family <- results(all_deseq_family, alpha=0.05, contrast=c("species", "Hihi", "Kaki"))
summary(deseq_res_hihi_vs_kaki_family) #Out of 343 familes with nonzero total read count and adj-p < 0.5 there are 2 (0.58%) increased (i.e. greater proportion in hihi than kaki) when comparing hihi to kaki, and 90 (26%) decreased (i.e. at a lower count abundance in the hihi than in the kaki samples) = 92 significant in total
#With adj-p <0.01 2 increased and 59 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_deseq_res_hihi_vs_kaki_family <- deseq_res_hihi_vs_kaki_family[which(deseq_res_hihi_vs_kaki_family$padj < 0.05), ]
summary(sigtab_res_deseq_res_hihi_vs_kaki_family) #Table only contains families considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_deseq_res_hihi_vs_kaki_with_tax_family <- cbind(as(sigtab_res_deseq_res_hihi_vs_kaki_family, "data.frame"), as(tax_table(ASV_physeq_rename_families)[row.names(sigtab_res_deseq_res_hihi_vs_kaki_family), ], "matrix"))
dim(sigtab_res_deseq_res_hihi_vs_kaki_with_tax_family) #92

#Sort by the log2FoldChange column
sigtab_res_deseq_res_hihi_vs_kaki_with_tax_family[order(sigtab_res_deseq_res_hihi_vs_kaki_with_tax_family$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the 92 significant genera = negative log2FoldChange means lower count abundance in hihi vs. kaki, positive log2FoldChange means greater count abundance in hihi vs. kaki
#write.csv(sigtab_res_deseq_res_hihi_vs_kaki_with_tax_family,"DifferentialAbundance/DESeq2_Results/deseq_results_hihi.vs.kaki_qiime2_family.csv")

#Grouped by phyla level
all_deseq_phylum <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_phylum, ~species) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
all_deseq_phylum <- DESeq(all_deseq_phylum)
deseq_res_hihi_vs_kaki_phylum <- results(all_deseq_phylum, alpha=0.05, contrast=c("species", "Hihi", "Kaki"))
summary(deseq_res_hihi_vs_kaki_phylum) #Out of 32 phyla with nonzero total read count and adj-p < 0.5 there is 1 (3.1%) increased (i.e. greater proportion in hihi than kaki) when comparing hihi to kaki, and 0 (0%) decreased (i.e. at a lower count abundance in the hihi than in the kaki samples) = 1 significant in total
#With adj-p <0.01 1 increased and 0 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_deseq_res_hihi_vs_kaki_phylum <- deseq_res_hihi_vs_kaki_phylum[which(deseq_res_hihi_vs_kaki_phylum$padj < 0.05), ]
summary(sigtab_res_deseq_res_hihi_vs_kaki_phylum) #Table only contains phyla considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_deseq_res_hihi_vs_kaki_with_tax_phylum <- cbind(as(sigtab_res_deseq_res_hihi_vs_kaki_phylum, "data.frame"), as(tax_table(ASV_physeq_rename_phyla)[row.names(sigtab_res_deseq_res_hihi_vs_kaki_phylum), ], "matrix"))
dim(sigtab_res_deseq_res_hihi_vs_kaki_with_tax_phylum) #1

#Sort by the log2FoldChange column
sigtab_res_deseq_res_hihi_vs_kaki_with_tax_phylum[order(sigtab_res_deseq_res_hihi_vs_kaki_with_tax_phylum$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the significant phylum = negative log2FoldChange means lower count abundance in hihi vs. kaki, positive log2FoldChange means greater count abundance in hihi vs. kaki
#write.csv(sigtab_res_deseq_res_hihi_vs_kaki_with_tax_phylum,"DifferentialAbundance/DESeq2_Results/deseq_results_hihi.vs.kaki_qiime2_phylum.csv")

## DESeq2 (https://astrobiomike.github.io/amplicon/dada2_workflow_ex)

## Kaki only

# Wild vs. captive 

# ASV level
kaki_deseq_wild_vs_captive <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_kaki, ~wild.captive) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_wild_vs_captive <- DESeq(kaki_deseq_wild_vs_captive)
kaki_deseq_res_wild_vs_captive <- results(kaki_deseq_wild_vs_captive, alpha=0.05, contrast=c("wild.captive", "Wild", "Captive"))
summary(kaki_deseq_res_wild_vs_captive) #Out of 2662 ASVs with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in wild than captive) when comparing wild to captive, and 0 (0.0%) decreased (i.e. at a lower count abundance in the wild than in the captive samples) = 0 significant in total
#With adj-p <0.01 0 increased and 0 decreased

# Genus level
kaki_deseq_genus_wild_vs_captive <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_genus_NArmFalse_kaki, ~wild.captive) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_genus_wild_vs_captive <- DESeq(kaki_deseq_genus_wild_vs_captive)
kaki_deseq_genus_res_wild_vs_captive <- results(kaki_deseq_genus_wild_vs_captive, alpha=0.05, contrast=c("wild.captive", "Wild", "Captive"))
summary(kaki_deseq_genus_res_wild_vs_captive) #Out of 710 genera with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in wild than captive) when comparing wild to captive, and 0 (0.0%) decreased (i.e. at a lower count abundance in the wild than in the captive samples) = 0 significant in total
#With adj-p <0.01 0 increased and 0 decreased

# Phylum level
kaki_deseq_phylum_wild_vs_captive <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_phylum_kaki, ~wild.captive) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_phylum_wild_vs_captive <- DESeq(kaki_deseq_phylum_wild_vs_captive)
kaki_deseq_phylum_res_wild_vs_captive <- results(kaki_deseq_phylum_wild_vs_captive, alpha=0.05, contrast=c("wild.captive", "Wild", "Captive"))
summary(kaki_deseq_phylum_res_wild_vs_captive) #Out of 32 phyla with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in wild than captive) when comparing wild to captive, and 2 (6.2%) decreased (i.e. at a lower count abundance in the wild than in the captive samples) = 2 significant in total
#With adj-p <0.01 0 increased and 0 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_kaki_deseq_phylum_res_wild_vs_captive <- kaki_deseq_phylum_res_wild_vs_captive[which(kaki_deseq_phylum_res_wild_vs_captive$padj < 0.05), ]
summary(sigtab_res_kaki_deseq_phylum_res_wild_vs_captive) #Table only contains phyla considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_kaki_deseq_phylum_res_wild_vs_captive_with_tax <- cbind(as(sigtab_res_kaki_deseq_phylum_res_wild_vs_captive, "data.frame"), as(tax_table(ASV_physeq_rename_phyla)[row.names(sigtab_res_kaki_deseq_phylum_res_wild_vs_captive), ], "matrix"))
dim(sigtab_res_kaki_deseq_phylum_res_wild_vs_captive_with_tax) #2

#Sort by the log2FoldChange column
sigtab_res_kaki_deseq_phylum_res_wild_vs_captive_with_tax[order(sigtab_res_kaki_deseq_phylum_res_wild_vs_captive_with_tax$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the significant phylum = negative log2FoldChange means lower count abundance in Wild vs. Captive, positive log2FoldChange means greater count abundance in Wild vs. Captive
#write.csv(sigtab_res_kaki_deseq_phylum_res_wild_vs_captive_with_tax,"DifferentialAbundance/DESeq2_Results/deseq_results_kaki_wild_vs_captive_qiime2_phylum.csv")

# Male experience

# ASV level
kaki_deseq_male_experience <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_kaki, ~male_experienced.naive_social) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_male_experience <- DESeq(kaki_deseq_male_experience)
kaki_deseq_res_male_experience <- results(kaki_deseq_male_experience, alpha=0.05, contrast=c("male_experienced.naive_social", "Experienced", "Naive"))
summary(kaki_deseq_res_male_experience) #Out of 2535 ASVs with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in Experienced than Naive) when comparing Experienced to Naive, and 0 (0.0%) decreased (i.e. at a lower count abundance in the Experienced than in the Naive samples) = 0 significant in total
#With adj-p <0.01 0 increased and 0 decreased

# Genus level
kaki_deseq_genus_male_experience <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_genus_NArmFalse_kaki, ~male_experienced.naive_social) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_genus_male_experience <- DESeq(kaki_deseq_genus_male_experience)
kaki_deseq_genus_res_male_experience <- results(kaki_deseq_genus_male_experience, alpha=0.05, contrast=c("male_experienced.naive_social", "Experienced", "Naive"))
summary(kaki_deseq_genus_res_male_experience) #Out of 695 genera with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in Experienced than Naive) when comparing Experienced to Naive, and 0 (0.0%) decreased (i.e. at a lower count abundance in the Experienced than in the Naive samples) = 0 significant in total
#With adj-p <0.01 0 increased and 0 decreased

# Phylum level
kaki_deseq_phylum_male_experience <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_phylum_kaki, ~male_experienced.naive_social) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_phylum_male_experience <- DESeq(kaki_deseq_phylum_male_experience)
kaki_deseq_phylum_res_male_experience <- results(kaki_deseq_phylum_male_experience, alpha=0.05, contrast=c("male_experienced.naive_social", "Experienced", "Naive"))
summary(kaki_deseq_phylum_res_male_experience) #Out of 32 phyla with nonzero total read count and adj-p < 0.5 there are 2 (6.2%) increased (i.e. greater proportion in Experienced than Naive) when comparing Experienced to Naive, and 0 (0.0%) decreased (i.e. at a lower count abundance in the Experienced than in the Naive samples) = 2 significant in total
#With adj-p <0.01 2 increased and 0 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_kaki_deseq_phylum_res_male_experience <- kaki_deseq_phylum_res_male_experience[which(kaki_deseq_phylum_res_male_experience$padj < 0.05), ]
summary(sigtab_res_kaki_deseq_phylum_res_male_experience) #Table only contains phyla considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_kaki_deseq_phylum_res_male_experience_with_tax <- cbind(as(sigtab_res_kaki_deseq_phylum_res_male_experience, "data.frame"), as(tax_table(ASV_physeq_rename_phyla)[row.names(sigtab_res_kaki_deseq_phylum_res_male_experience), ], "matrix"))
dim(sigtab_res_kaki_deseq_phylum_res_male_experience_with_tax) #2

#Sort by the log2FoldChange column
sigtab_res_kaki_deseq_phylum_res_male_experience_with_tax[order(sigtab_res_kaki_deseq_phylum_res_male_experience_with_tax$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the significant phylum = negative log2FoldChange means lower count abundance in Experienced vs. Naive, positive log2FoldChange means greater count abundance in Experienced vs. Naive
#write.csv(sigtab_res_kaki_deseq_phylum_res_male_experience_with_tax,"DifferentialAbundance/DESeq2_Results/deseq_results_kaki_male_experience_qiime2_phylum.csv")

# Male wild- or captive-laid

# ASV level
kaki_deseq_male_wild.captive <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_kaki, ~male_lay_wild.captive) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_male_wild.captive <- DESeq(kaki_deseq_male_wild.captive)
kaki_deseq_res_male_wild.captive <- results(kaki_deseq_male_wild.captive, alpha=0.05, contrast=c("male_lay_wild.captive", "Wild", "Captive"))
summary(kaki_deseq_res_male_wild.captive) #Out of 2551 ASVs with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in Wild than Captive) when comparing Wild to Captive, and 2 (0.078%) decreased (i.e. at a lower count abundance in the Wild than in the Captive samples) = 2 significant in total
#With adj-p <0.01 0 increased and 2 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_kaki_deseq_res_male_wild.captive <- kaki_deseq_res_male_wild.captive[which(kaki_deseq_res_male_wild.captive$padj < 0.05), ]
summary(sigtab_res_kaki_deseq_res_male_wild.captive) #Table only contains ASVs considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_kaki_deseq_res_male_wild.captive_with_tax <- cbind(as(sigtab_res_kaki_deseq_res_male_wild.captive, "data.frame"), as(tax_table(ASV_physeq)[row.names(sigtab_res_kaki_deseq_res_male_wild.captive), ], "matrix"))
dim(sigtab_res_kaki_deseq_res_male_wild.captive_with_tax) #2

#Sort by the log2FoldChange column
sigtab_res_kaki_deseq_res_male_wild.captive_with_tax[order(sigtab_res_kaki_deseq_res_male_wild.captive_with_tax$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the significant phylum = negative log2FoldChange means lower count abundance in Wild vs. Captive, positive log2FoldChange means greater count abundance in Wild vs. Captive
#write.csv(sigtab_res_kaki_deseq_res_male_wild.captive_with_tax,"DifferentialAbundance/DESeq2_Results/deseq_results_kaki_male_wild.captive_qiime2.csv")

# Genus level
kaki_deseq_genus_male_wild.captive <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_genus_NArmFalse_kaki, ~male_lay_wild.captive) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_genus_male_wild.captive <- DESeq(kaki_deseq_genus_male_wild.captive)
kaki_deseq_genus_res_male_wild.captive <- results(kaki_deseq_genus_male_wild.captive, alpha=0.05, contrast=c("male_lay_wild.captive", "Wild", "Captive"))
summary(kaki_deseq_genus_res_male_wild.captive) #Out of 699 genera with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in Wild than Captive) when comparing Wild to Captive, and 1 (0.14%) decreased (i.e. at a lower count abundance in the Wild than in the Captive samples) = 1 significant in total
#With adj-p <0.01 0 increased and 1 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_kaki_deseq_genus_res_male_wild.captive <- kaki_deseq_genus_res_male_wild.captive[which(kaki_deseq_genus_res_male_wild.captive$padj < 0.05), ]
summary(sigtab_res_kaki_deseq_genus_res_male_wild.captive) #Table only contains ASVs considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_kaki_deseq_genus_res_male_wild.captive_with_tax <- cbind(as(sigtab_res_kaki_deseq_genus_res_male_wild.captive, "data.frame"), as(tax_table(ASV_physeq_rename_genera)[row.names(sigtab_res_kaki_deseq_genus_res_male_wild.captive), ], "matrix"))
dim(sigtab_res_kaki_deseq_genus_res_male_wild.captive_with_tax) #1

#Sort by the log2FoldChange column
sigtab_res_kaki_deseq_genus_res_male_wild.captive_with_tax[order(sigtab_res_kaki_deseq_genus_res_male_wild.captive_with_tax$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the significant phylum = negative log2FoldChange means lower count abundance in Wild vs. Captive, positive log2FoldChange means greater count abundance in Experienced vs. Naive
#write.csv(sigtab_res_kaki_deseq_genus_res_male_wild.captive_with_tax,"DifferentialAbundance/DESeq2_Results/deseq_results_kaki_male_wild.captive_qiime2_genus.csv")

# Phylum level
kaki_deseq_phylum_male_wild.captive <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_phylum_kaki, ~male_lay_wild.captive) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_phylum_male_wild.captive <- DESeq(kaki_deseq_phylum_male_wild.captive)
kaki_deseq_phylum_res_male_wild.captive <- results(kaki_deseq_phylum_male_wild.captive, alpha=0.05, contrast=c("male_lay_wild.captive", "Wild", "Captive"))
summary(kaki_deseq_phylum_res_male_wild.captive) #Out of 31 phyla with nonzero total read count and adj-p < 0.5 there is 1 (3.2%) increased (i.e. greater proportion in Wild than Captive) when comparing Wild to Captive, and 0 (0.0%) decreased (i.e. at a lower count abundance in the Wild than in the Captive samples) = 1 significant in total
#With adj-p <0.01 0 increased and 0 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_kaki_deseq_phylum_res_male_wild.captive <- kaki_deseq_phylum_res_male_wild.captive[which(kaki_deseq_phylum_res_male_wild.captive$padj < 0.05), ]
summary(sigtab_res_kaki_deseq_phylum_res_male_wild.captive) #Table only contains phyla considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_kaki_deseq_phylum_res_male_wild.captive_with_tax <- cbind(as(sigtab_res_kaki_deseq_phylum_res_male_wild.captive, "data.frame"), as(tax_table(ASV_physeq_rename_phyla)[row.names(sigtab_res_kaki_deseq_phylum_res_male_wild.captive), ], "matrix"))
dim(sigtab_res_kaki_deseq_phylum_res_male_wild.captive_with_tax) #1

#Sort by the log2FoldChange column
sigtab_res_kaki_deseq_phylum_res_male_wild.captive_with_tax[order(sigtab_res_kaki_deseq_phylum_res_male_wild.captive_with_tax$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the significant phylum = negative log2FoldChange means lower count abundance in Wild vs. Captive, positive log2FoldChange means greater count abundance in Wild vs. Captive
#write.csv(sigtab_res_kaki_deseq_phylum_res_male_wild.captive_with_tax,"DifferentialAbundance/DESeq2_Results/deseq_results_kaki_male_wild.captive_qiime2_phyla.csv")

# Hatching success

# ASV level
kaki_deseq_hatched.unhatched <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_kaki, ~hatched.unhatched) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_hatched.unhatched <- DESeq(kaki_deseq_hatched.unhatched)
kaki_deseq_res_hatched.unhatched <- results(kaki_deseq_hatched.unhatched, alpha=0.05, contrast=c("hatched.unhatched", "Hatched", "Unhatched"))
summary(kaki_deseq_res_hatched.unhatched) #Out of 2705 ASVs with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in hatched than unhatched) when comparing hatched to unhatched, and 0 (0.0%) decreased (i.e. at a lower count abundance in the hatched than in the unhatched samples) = 0 significant in total
#With adj-p <0.01 0 increased and 0 decreased

# Genus level
kaki_deseq_genus_hatched.unhatched <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_genus_NArmFalse_kaki, ~hatched.unhatched) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_genus_hatched.unhatched <- DESeq(kaki_deseq_genus_hatched.unhatched)
kaki_deseq_genus_res_hatched.unhatched <- results(kaki_deseq_genus_hatched.unhatched, alpha=0.05, contrast=c("hatched.unhatched", "Hatched", "Unhatched"))
summary(kaki_deseq_genus_res_hatched.unhatched) #Out of 720 genera with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in hatched than unhatched) when comparing hatched to unhatched, and 0 (0.0%) decreased (i.e. at a lower count abundance in the hatched than in the unhatched samples) = 0 significant in total
#With adj-p <0.01 0 increased and 0 decreased

# Family level
kaki_deseq_family_hatched.unhatched <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_family_NArmFalse_kaki, ~hatched.unhatched) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_family_hatched.unhatched <- DESeq(kaki_deseq_family_hatched.unhatched)
kaki_deseq_family_res_hatched.unhatched <- results(kaki_deseq_family_hatched.unhatched, alpha=0.05, contrast=c("hatched.unhatched", "Hatched", "Unhatched"))
summary(kaki_deseq_family_res_hatched.unhatched) #Out of 349 genera with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in hatched than unhatched) when comparing hatched to unhatched, and 0 (0.0%) decreased (i.e. at a lower count abundance in the hatched than in the unhatched samples) = 0 significant in total
#With adj-p <0.01 0 increased and 0 decreased

# Phylum level
kaki_deseq_phylum_hatched.unhatched <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_phylum_kaki, ~hatched.unhatched) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
kaki_deseq_phylum_hatched.unhatched <- DESeq(kaki_deseq_phylum_hatched.unhatched)
kaki_deseq_phylum_res_hatched.unhatched <- results(kaki_deseq_phylum_hatched.unhatched, alpha=0.05, contrast=c("hatched.unhatched", "Hatched", "Unhatched"))
summary(kaki_deseq_phylum_res_hatched.unhatched) #Out of 32 phyla with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in hatched than unhatched) when comparing hatched to unhatched, and 0 (0.0%) decreased (i.e. at a lower count abundance in the hatched than in the unhatched samples) = 0 significant in total
#With adj-p <0.01 0 increased and 0 decreased

#

## Hihi only

# Female experience

# ASV level
hihi_deseq_female_experience <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_hihi, ~female_experienced.naive_social) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
hihi_deseq_female_experience <- DESeq(hihi_deseq_female_experience)
hihi_deseq_res_female_experience <- results(hihi_deseq_female_experience, alpha=0.05, contrast=c("female_experienced.naive_social", "Experienced", "Naive"))
summary(hihi_deseq_res_female_experience) #Out of 452 ASVs with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in Experienced than Naive) when comparing Experienced to Naive, and 2 (0.44%) decreased (i.e. at a lower count abundance in the Experienced than in the Naive samples) = 2 significant in total
#With adj-p <0.01 0 increased and 1 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_hihi_deseq_res_female_experience <- hihi_deseq_res_female_experience[which(hihi_deseq_res_female_experience$padj < 0.05), ]
summary(sigtab_res_hihi_deseq_res_female_experience) #Table only contains ASVs considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_hihi_deseq_res_female_experience_with_tax <- cbind(as(sigtab_res_hihi_deseq_res_female_experience, "data.frame"), as(tax_table(ASV_physeq_rename)[row.names(sigtab_res_hihi_deseq_res_female_experience), ], "matrix"))
dim(sigtab_res_hihi_deseq_res_female_experience_with_tax) #2

#Sort by the log2FoldChange column
sigtab_res_hihi_deseq_res_female_experience_with_tax[order(sigtab_res_hihi_deseq_res_female_experience_with_tax$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the significant phylum = negative log2FoldChange means lower count abundance in Experienced vs. Naive, positive log2FoldChange means greater count abundance in Experienced vs. Naive
#write.csv(sigtab_res_hihi_deseq_res_female_experience_with_tax,"DifferentialAbundance/DESeq2_Results/deseq_results_hihi_male_female_experience_qiime2.csv")

# Genus level
hihi_deseq_genus_female_experience <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_genus_NArmFalse_hihi, ~female_experienced.naive_social) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
hihi_deseq_genus_female_experience <- DESeq(hihi_deseq_genus_female_experience)
hihi_deseq_genus_res_female_experience <- results(hihi_deseq_genus_female_experience, alpha=0.05, contrast=c("female_experienced.naive_social", "Experienced", "Naive"))
summary(hihi_deseq_genus_res_female_experience) #Out of 273 genera with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in Experienced than Naive) when comparing Experienced to Naive, and 1 (0.37%) decreased (i.e. at a lower count abundance in the Experienced than in the Naive samples) = 1 significant in total
#With adj-p <0.01 0 increased and 1 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_hihi_deseq_genus_res_female_experience <- hihi_deseq_genus_res_female_experience[which(hihi_deseq_genus_res_female_experience$padj < 0.05), ]
summary(sigtab_res_hihi_deseq_genus_res_female_experience) #Table only contains ASVs considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_hihi_deseq_genus_res_female_experience_with_tax <- cbind(as(sigtab_res_hihi_deseq_genus_res_female_experience, "data.frame"), as(tax_table(ASV_physeq_rename_genera)[row.names(sigtab_res_hihi_deseq_genus_res_female_experience), ], "matrix"))
dim(sigtab_res_hihi_deseq_genus_res_female_experience_with_tax) #1

#Sort by the log2FoldChange column
sigtab_res_hihi_deseq_genus_res_female_experience_with_tax[order(sigtab_res_hihi_deseq_genus_res_female_experience_with_tax$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the significant phylum = negative log2FoldChange means lower count abundance in Experienced vs. Naive, positive log2FoldChange means greater count abundance in Experienced vs. Naive
#write.csv(sigtab_res_hihi_deseq_genus_res_female_experience_with_tax,"DifferentialAbundance/DESeq2_Results/deseq_results_hihi_male_female_experience_qiime2_genus.csv")

# Phylum level
hihi_deseq_phylum_female_experience <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_phylum_hihi, ~female_experienced.naive_social) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
hihi_deseq_phylum_female_experience <- DESeq(hihi_deseq_phylum_female_experience)
hihi_deseq_phylum_res_female_experience <- results(hihi_deseq_phylum_female_experience, alpha=0.05, contrast=c("female_experienced.naive_social", "Experienced", "Naive"))
summary(hihi_deseq_phylum_res_female_experience) #Out of 18 phyla with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in Experienced than Naive) when comparing Experienced to Naive, and 0 (0.0%) decreased (i.e. at a lower count abundance in the Experienced than in the Naive samples) = 0 significant in total
#With adj-p <0.01 0 increased and 0 decreased

# Hatching success

#Need to remove NAs
physeq_tree_noMCnoNTC_hihi_hatched <- subset_samples(physeq_tree_noMCnoNTC_hihi, hatched.unhatched == "Hatched")
physeq_tree_noMCnoNTC_hihi_unhatched <- subset_samples(physeq_tree_noMCnoNTC_hihi, hatched.unhatched == "Unhatched")
physeq_tree_noMCnoNTC_hihi_HS <- merge_phyloseq(physeq_tree_noMCnoNTC_hihi_hatched,physeq_tree_noMCnoNTC_hihi_unhatched)
sample_data(physeq_tree_noMCnoNTC_hihi_HS)$hatched.unhatched

# ASV level
hihi_deseq_hatched.unhatched <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_hihi_HS, ~hatched.unhatched) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
hihi_deseq_hatched.unhatched <- DESeq(hihi_deseq_hatched.unhatched)
hihi_deseq_res_hatched.unhatched <- results(hihi_deseq_hatched.unhatched, alpha=0.05, contrast=c("hatched.unhatched", "Hatched", "Unhatched"))
summary(hihi_deseq_res_hatched.unhatched) #Out of 332 ASVs with nonzero total read count and adj-p < 0.5 there are 1 (0.3%) increased (i.e. greater proportion in hatched than unhatched) when comparing hatched to unhatched, and 0 (0.0%) decreased (i.e. at a lower count abundance in the hatched than in the unhatched samples) = 0 significant in total
#With adj-p <0.01 1 increased and 0 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_hihi_deseq_res_hatched.unhatched <- hihi_deseq_res_hatched.unhatched[which(hihi_deseq_res_hatched.unhatched$padj < 0.05), ]
summary(sigtab_res_hihi_deseq_res_hatched.unhatched) #Table only contains ASVs considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_hihi_deseq_res_hatched.unhatched_with_tax <- cbind(as(sigtab_res_hihi_deseq_res_hatched.unhatched, "data.frame"), as(tax_table(ASV_physeq_rename)[row.names(sigtab_res_hihi_deseq_res_hatched.unhatched), ], "matrix"))
dim(sigtab_res_hihi_deseq_res_hatched.unhatched_with_tax) #1

#Sort by the log2FoldChange column
sigtab_res_hihi_deseq_res_hatched.unhatched_with_tax[order(sigtab_res_hihi_deseq_res_hatched.unhatched_with_tax$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the significant phylum = negative log2FoldChange means lower count abundance in hatched vs. unhatched, positive log2FoldChange means greater count abundance in hatched vs. unhatched
#write.csv(sigtab_res_hihi_deseq_res_hatched.unhatched_with_tax,"DifferentialAbundance/DESeq2_Results/deseq_results_hihi_hatched.unhatched_qiime2.csv")

#Need to remove NAs
physeq_tree_noMCnoNTC_genus_NArmFalse_hihi_hatched <- subset_samples(physeq_tree_noMCnoNTC_genus_NArmFalse_hihi, hatched.unhatched == "Hatched")
physeq_tree_noMCnoNTC_genus_NArmFalse_hihi_unhatched <- subset_samples(physeq_tree_noMCnoNTC_genus_NArmFalse_hihi, hatched.unhatched == "Unhatched")
physeq_tree_noMCnoNTC_genus_NArmFalse_hihi_HS <- merge_phyloseq(physeq_tree_noMCnoNTC_genus_NArmFalse_hihi_hatched,physeq_tree_noMCnoNTC_genus_NArmFalse_hihi_unhatched)
sample_data(physeq_tree_noMCnoNTC_genus_NArmFalse_hihi_HS)$hatched.unhatched

# Genus level
hihi_deseq_genus_hatched.unhatched <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_genus_NArmFalse_hihi_HS, ~hatched.unhatched) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
hihi_deseq_genus_hatched.unhatched <- DESeq(hihi_deseq_genus_hatched.unhatched)
hihi_deseq_genus_res_hatched.unhatched <- results(hihi_deseq_genus_hatched.unhatched, alpha=0.05, contrast=c("hatched.unhatched", "Hatched", "Unhatched"))
summary(hihi_deseq_genus_res_hatched.unhatched) #Out of 219 genera with nonzero total read count and adj-p < 0.5 there are 4 (1.8%) increased (i.e. greater proportion in hatched than unhatched) when comparing hatched to unhatched, and 0 (0.0%) decreased (i.e. at a lower count abundance in the hatched than in the unhatched samples) = 0 significant in total
#With adj-p <0.01 4 increased and 0 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_hihi_deseq_genus_res_hatched.unhatched <- hihi_deseq_genus_res_hatched.unhatched[which(hihi_deseq_genus_res_hatched.unhatched$padj < 0.05), ]
summary(sigtab_res_hihi_deseq_genus_res_hatched.unhatched) #Table only contains ASVs considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_hihi_deseq_genus_res_hatched.unhatched_with_tax <- cbind(as(sigtab_res_hihi_deseq_genus_res_hatched.unhatched, "data.frame"), as(tax_table(ASV_physeq_rename_genera)[row.names(sigtab_res_hihi_deseq_genus_res_hatched.unhatched), ], "matrix"))
dim(sigtab_res_hihi_deseq_genus_res_hatched.unhatched_with_tax) #4

#Sort by the log2FoldChange column
sigtab_res_hihi_deseq_genus_res_hatched.unhatched_with_tax[order(sigtab_res_hihi_deseq_genus_res_hatched.unhatched_with_tax$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the significant phylum = negative log2FoldChange means lower count abundance in hatched vs. unhatched, positive log2FoldChange means greater count abundance in hatched vs. unhatched
#write.csv(sigtab_res_hihi_deseq_genus_res_hatched.unhatched_with_tax,"DifferentialAbundance/DESeq2_Results/deseq_results_hihi_hatched.unhatched_qiime2_genus.csv")

#Need to remove NAs
physeq_tree_noMCnoNTC_family_NArmFalse_hihi_hatched <- subset_samples(physeq_tree_noMCnoNTC_family_NArmFalse_hihi, hatched.unhatched == "Hatched")
physeq_tree_noMCnoNTC_family_NArmFalse_hihi_unhatched <- subset_samples(physeq_tree_noMCnoNTC_family_NArmFalse_hihi, hatched.unhatched == "Unhatched")
physeq_tree_noMCnoNTC_family_NArmFalse_hihi_HS <- merge_phyloseq(physeq_tree_noMCnoNTC_family_NArmFalse_hihi_hatched,physeq_tree_noMCnoNTC_family_NArmFalse_hihi_unhatched)
sample_data(physeq_tree_noMCnoNTC_family_NArmFalse_hihi_HS)$hatched.unhatched

# Family level
hihi_deseq_family_hatched.unhatched <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_family_NArmFalse_hihi_HS, ~hatched.unhatched) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
hihi_deseq_family_hatched.unhatched <- DESeq(hihi_deseq_family_hatched.unhatched)
hihi_deseq_family_res_hatched.unhatched <- results(hihi_deseq_family_hatched.unhatched, alpha=0.05, contrast=c("hatched.unhatched", "Hatched", "Unhatched"))
summary(hihi_deseq_family_res_hatched.unhatched) #Out of 124 families with nonzero total read count and adj-p < 0.5 there are 2 (1.6%) increased (i.e. greater proportion in hatched than unhatched) when comparing hatched to unhatched, and 0 (0.0%) decreased (i.e. at a lower count abundance in the hatched than in the unhatched samples) = 0 significant in total
#With adj-p <0.01 2 increased and 0 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_hihi_deseq_family_res_hatched.unhatched <- hihi_deseq_family_res_hatched.unhatched[which(hihi_deseq_family_res_hatched.unhatched$padj < 0.05), ]
summary(sigtab_res_hihi_deseq_family_res_hatched.unhatched) #Table only contains ASVs considered significantly differentially abundant

#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_hihi_deseq_family_res_hatched.unhatched_with_tax <- cbind(as(sigtab_res_hihi_deseq_family_res_hatched.unhatched, "data.frame"), as(tax_table(ASV_physeq_rename_families)[row.names(sigtab_res_hihi_deseq_family_res_hatched.unhatched), ], "matrix"))
dim(sigtab_res_hihi_deseq_family_res_hatched.unhatched_with_tax) #2

#Sort by the log2FoldChange column
sigtab_res_hihi_deseq_family_res_hatched.unhatched_with_tax[order(sigtab_res_hihi_deseq_family_res_hatched.unhatched_with_tax$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the significant phylum = negative log2FoldChange means lower count abundance in hatched vs. unhatched, positive log2FoldChange means greater count abundance in hatched vs. unhatched
#write.csv(sigtab_res_hihi_deseq_family_res_hatched.unhatched_with_tax,"DifferentialAbundance/DESeq2_Results/deseq_results_hihi_hatched.unhatched_qiime2_family.csv")

#Need to remove NAs
physeq_tree_noMCnoNTC_phylum_hihi_hatched <- subset_samples(physeq_tree_noMCnoNTC_phylum_hihi, hatched.unhatched == "Hatched")
physeq_tree_noMCnoNTC_phylum_hihi_unhatched <- subset_samples(physeq_tree_noMCnoNTC_phylum_hihi, hatched.unhatched == "Unhatched")
physeq_tree_noMCnoNTC_phylum_hihi_HS <- merge_phyloseq(physeq_tree_noMCnoNTC_phylum_hihi_hatched,physeq_tree_noMCnoNTC_phylum_hihi_unhatched)
sample_data(physeq_tree_noMCnoNTC_phylum_hihi_HS)$hatched.unhatched

# Phylum level
hihi_deseq_phylum_hatched.unhatched <- phyloseq_to_deseq2(physeq_tree_noMCnoNTC_phylum_hihi_HS, ~hatched.unhatched) #Warning message: In DESeqDataSet(se, design = design, ignoreRank) :some variables in design formula are characters, converting to factors
hihi_deseq_phylum_hatched.unhatched <- DESeq(hihi_deseq_phylum_hatched.unhatched)
hihi_deseq_phylum_res_hatched.unhatched <- results(hihi_deseq_phylum_hatched.unhatched, alpha=0.05, contrast=c("hatched.unhatched", "Hatched", "Unhatched"))
summary(hihi_deseq_phylum_res_hatched.unhatched) #Out of 18 phyla with nonzero total read count and adj-p < 0.5 there are 0 (0.0%) increased (i.e. greater proportion in hatched than unhatched) when comparing hatched to unhatched, and 0 (0.0%) decreased (i.e. at a lower count abundance in the hatched than in the unhatched samples) = 0 significant in total
#With adj-p <0.01 0 increased and 0 decreased
