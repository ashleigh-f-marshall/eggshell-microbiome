####################
# Authors: Ashleigh Fleming Marshall
# ashleigh.marshall@ioz.ac.uk | ashleigh.marshall.16@ucl.ac.uk
####################

### MICROBIOME, DISEASE, AND HATCHING FAILURE

##Clear workspace
rm(list=ls())

## Figures for microbiome chapters

# Load packages
library("ggplot2"); packageVersion("ggplot2") #‘3.4.0’
library("microbial"); packageVersion("microbial") #‘0.0.20’
library("qiime2R"); packageVersion("qiime2R") #‘0.99.6’
library("phyloseq"); packageVersion("phyloseq") #‘1.40.0’
library("ggtext"); packageVersion("ggtext") #‘0.1.2’
library("patchwork"); packageVersion("patchwork") #'1.1.2'
library("ggview"); packageVersion("ggview") #‘0.1.0’
library("tidyverse"); packageVersion("tidyverse") #‘1.3.2’
library("RColorBrewer"); packageVersion("RColorBrewer") #‘1.1.3’
library("scales"); packageVersion("scales") #‘1.2.1’
library("ggrepel"); packageVersion("ggrepel") #‘0.9.2’
library("glue"); packageVersion("glue") #'1.6.2'
library("ggalt"); packageVersion("ggalt") #'0.4.0’
library("colorBlindness"); packageVersion("colorBlindness") #'0.1.9’
library("ggsignif"); packageVersion("ggsignif") #'0.6.4’
library("vegan"); packageVersion("vegan") #‘2.6.4’
library("ggforce"); packageVersion("ggforce") #‘0.4.1’
library("grafify"); packageVersion("grafify") #‘3.0.1’
library("ggh4x"); packageVersion("ggh4x") #‘0.2.3’

#for(i in 1:14){
#  print(hue_pal()(i))
#}

#show_col(viridis_pal()(3))

#display.brewer.all(colorblindFriendly = TRUE)

setwd("Data/DataAnalyses")

# Import metadata
mtd <- read_tsv(("sample-metadata-new-noNA-removemock_forR_wEnvironmentalData_2.tsv"), comment = "#q2")

# Import feature table and convert to count table
table_250.256_silva_138.1.filtered_decontam_noMCnoNTC <- read_qza("filtered-silva-138.1-table-decontam_0.4_noMCnoNTC.qza")
count_tab_250.256_silva_138.1.filtered_decontam_noMCnoNTC <- table_250.256_silva_138.1.filtered_decontam_noMCnoNTC$data %>% as.data.frame() 

#Import taxonomies
taxonomy_silva_138.1_250.256_filtered_decontam_noMCnoNTC <- read_qza("silva_138.1_taxonomy_250.256-filtered-decontam_noMCnoNTC.qza")

tax_tab_silva_138.1_250.256_filtered_decontam <- taxonomy_silva_138.1_250.256_filtered_decontam_noMCnoNTC$data %>% 
  as.data.frame() %>%
  mutate(Taxon = gsub("D_0", "k", Taxon), Taxon = gsub("D_1", "p", Taxon),
         Taxon = gsub("D_2", "c", Taxon), Taxon = gsub("D_3", "o", Taxon),
         Taxon = gsub("D_4", "f", Taxon), Taxon = gsub("D_5", "g", Taxon),
         Taxon = gsub("D_6", "s", Taxon)) %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  select(-Confidence) #Warning message: Expected 7 pieces. Missing pieces filled with `NA` in 4525 rows [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, ...] = this matches the number of NAs in the species column
dim(tax_tab_silva_138.1_250.256_filtered_decontam) #6528 7

test_physeq <- qza_to_phyloseq(features="filtered-silva-138.1-table-decontam_0.4_noMCnoNTC.qza",
                               tree="rooted-tree_allfilter_decontam_noMCnoNTC.qza",
                               taxonomy="silva_138.1_taxonomy_250.256-filtered-decontam_noMCnoNTC.qza")
test_physeq_sd <- phyloseq(sample_data(column_to_rownames(mtd, "SampleID")))
physeq_tree_noMCnoNTC <- merge_phyloseq(test_physeq,test_physeq_sd)
physeq_tree_noMCnoNTC #6528 taxa and 249 samples

##### PCoA (ordination plots) from Riffomonas tutorials
# Load data
bray_mat <- read_tsv("all_exported_diversity_metrics/bray_curtis_distance_matrix/distance-matrix.tsv")
jaccard_mat <- read_tsv("all_exported_diversity_metrics/jaccard_distance_matrix/distance-matrix.tsv")
unweighted_mat <- read_tsv("all_exported_diversity_metrics/unweighted_unifrac_distance_matrix/distance-matrix.tsv")
weighted_mat <- read_tsv("all_exported_diversity_metrics/weighted_unifrac_distance_matrix/distance-matrix.tsv")
colnames(bray_mat)[1] <- "SampleID"
colnames(jaccard_mat)[1] <- "SampleID"
colnames(unweighted_mat)[1] <- "SampleID"
colnames(weighted_mat)[1] <- "SampleID"
bray_mat_tbl <- as_tibble(bray_mat)
jaccard_mat_tbl <- as_tibble(jaccard_mat)
unweighted_mat_tbl <- as_tibble(unweighted_mat)
weighted_mat_tbl <- as_tibble(weighted_mat)

mtd_tbl <- as_tibble(mtd)
colnames(mtd_tbl)[9] <- "study_species"
mtd_tbl_noMCnoNTC <- mtd_tbl %>%
  filter(!is.na(female_combo))
colnames(mtd_tbl_noMCnoNTC)[29] <- "female_age_cohort.sq"
colnames(mtd_tbl_noMCnoNTC)[45] <- "male_age_cohort.sq"
colnames(mtd_tbl_noMCnoNTC)[62] <- "first_egg_lay"
colnames(mtd_tbl_noMCnoNTC)[79] <- "re.est_first_egg_lay"
colnames(mtd_tbl_noMCnoNTC)[89] <- "clutch.hatch.percent"
mtd_tbl_noMCnoNTC_hihi <- mtd_tbl_noMCnoNTC %>%
  filter(study_species=="Hihi")
mtd_tbl_noMCnoNTC_kaki <- mtd_tbl_noMCnoNTC %>%
  filter(study_species=="Kaki")
mtd_tbl_short <- mtd_tbl %>%
  select(SampleID,clutch.nest_id_corrected,study_species,female_combo,wild.captive,hatched.unhatched)
mtd_tbl_short$SampleID
mtd_tbl_short_noMCnoNTC <- mtd_tbl_short %>%
  filter(!is.na(female_combo))
mtd_tbl_short_noMCnoNTC$SampleID
mtd_tbl_short_noMCnoNTC$study_species<-gsub("Kaki","Kakī",mtd_tbl_short_noMCnoNTC$study_species)
mtd_tbl_short_noMCnoNTC$female_combo<-gsub("\\.","/",mtd_tbl_short_noMCnoNTC$female_combo)
mtd_tbl_short_noMCnoNTC$clutch.nest_id_corrected<-gsub("\\.","/",mtd_tbl_short_noMCnoNTC$clutch.nest_id_corrected)

dist_matrix_species_weighted <- weighted_mat %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_species_unweighted <- unweighted_mat %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_species_jaccard <- jaccard_mat %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_species_bray <- bray_mat %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()

mtd_tbl_noMCnoNTC_kaki %>% arrange(SampleID)

#adonis2(dist_matrix_kaki_weighted~female_combo + nest_location + swab_day_parentalInc + clutch_number, data = mtd_tbl_noMCnoNTC_kaki)
#adonis2(dist_matrix_kaki_unweighted~female_combo + nest_location + swab_day_parentalInc + clutch_number, data = mtd_tbl_noMCnoNTC_kaki)
#adonis2(dist_matrix_kaki_jaccard~female_combo + nest_location + swab_day_parentalInc + clutch_number, data = mtd_tbl_noMCnoNTC_kaki)
#adonis2(dist_matrix_kaki_bray~female_combo + nest_location + swab_day_parentalInc + clutch_number,  method = "bray", data = mtd_tbl_noMCnoNTC_kaki)

#adonis2(dist_matrix_hihi_weighted~male_age_cohort.sq + mean_rain_laytoswab + male_age_cohort + plate_ID, data = mtd_tbl_noMCnoNTC_hihi)
#adonis2(dist_matrix_hihi_unweighted~male_age_cohort.sq + mean_rain_laytoswab + male_age_cohort + plate_ID, data = mtd_tbl_noMCnoNTC_hihi)
#adonis2(dist_matrix_hihi_jaccard~male_age_cohort.sq + mean_rain_laytoswab + male_age_cohort + plate_ID, data = mtd_tbl_noMCnoNTC_hihi)
#adonis2(dist_matrix_hihi_bray~male_age_cohort.sq + mean_rain_laytoswab + male_age_cohort + plate_ID, data = mtd_tbl_noMCnoNTC_hihi)

# PCoA
pcoa_species_weighted <- cmdscale(dist_matrix_species_weighted, eig = TRUE, add= TRUE)
pcoa_species_unweighted <- cmdscale(dist_matrix_species_unweighted, eig = TRUE, add= TRUE)
pcoa_species_jaccard <- cmdscale(dist_matrix_species_jaccard, eig = TRUE, add= TRUE)
pcoa_species_bray <- cmdscale(dist_matrix_species_bray, eig = TRUE, add= TRUE)

positions_species_weighted <- pcoa_species_weighted$points
positions_species_unweighted <- pcoa_species_unweighted$points
positions_species_jaccard <- pcoa_species_jaccard$points
positions_species_bray <- pcoa_species_bray$points

colnames(positions_species_weighted) <- c("pcoa1","pcoa2")
colnames(positions_species_unweighted) <- c("pcoa1","pcoa2")
colnames(positions_species_jaccard) <- c("pcoa1","pcoa2")
colnames(positions_species_bray) <- c("pcoa1","pcoa2")

percent_explained_species_weighted <- 100 * pcoa_species_weighted$eig / sum(pcoa_species_weighted$eig)
percent_explained_species_unweighted <- 100 * pcoa_species_unweighted$eig / sum(pcoa_species_unweighted$eig)
percent_explained_species_jaccard <- 100 * pcoa_species_jaccard$eig / sum(pcoa_species_jaccard$eig)
percent_explained_species_bray <- 100 * pcoa_species_bray$eig / sum(pcoa_species_bray$eig)

pretty_pe_species_weighted <- format(round(percent_explained_species_weighted[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_species_unweighted <- format(round(percent_explained_species_unweighted[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_species_jaccard <- format(round(percent_explained_species_jaccard[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_species_bray <- format(round(percent_explained_species_bray[1:2], digits = 1), nsmall = 1, trim = TRUE)

labs_species_weighted <- c(glue("PCo 1 ({pretty_pe_species_weighted[1]}%)"),
          glue("PCo 2 ({pretty_pe_species_weighted[2]}%)"))
labs_species_unweighted <- c(glue("PCo 1 ({pretty_pe_species_unweighted[1]}%)"),
                           glue("PCo 2 ({pretty_pe_species_unweighted[2]}%)"))
labs_species_jaccard <- c(glue("PCo 1 ({pretty_pe_species_jaccard[1]}%)"),
                           glue("PCo 2 ({pretty_pe_species_jaccard[2]}%)"))
labs_species_bray <- c(glue("PCo 1 ({pretty_pe_species_bray[1]}%)"),
                           glue("PCo 2 ({pretty_pe_species_bray[2]}%)"))

PCoA_species_weighted <- positions_species_weighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = study_species )) +
  geom_point() +
  scale_color_discrete(name = "Study Species") +
  geom_encircle(expand = 0, aes(fill = study_species), alpha = 0.3) + 
  scale_fill_discrete(name = "Study Species") +
  labs(x = labs_species_weighted[1], y = labs_species_weighted[2]) +
  theme_bw() +
  ggtitle("Weighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
PCoA_species_weighted
PCoA_species_unweighted <- positions_species_unweighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = study_species )) +
  geom_point() +
  scale_color_discrete(name = "Study Species") +
  geom_encircle(expand = 0, aes(fill = study_species), alpha = 0.3) + 
  scale_fill_discrete(name = "Study Species") +
  labs(x = labs_species_unweighted[1], y = labs_species_unweighted[2]) +
  theme_bw() +
  ggtitle("Unweighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
PCoA_species_unweighted
PCoA_species_jaccard <- positions_species_jaccard %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = study_species )) +
  geom_point() +
  scale_color_discrete(name = "Study Species") +
  geom_encircle(expand = 0, aes(fill = study_species), alpha = 0.3) + 
  scale_fill_discrete(name = "Study Species") +
  labs(x = labs_species_jaccard[1], y = labs_species_jaccard[2]) +
  theme_bw() +
  ggtitle("Jaccard Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
PCoA_species_jaccard
PCoA_species_bray <- positions_species_bray %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = study_species )) +
  geom_point() +
  scale_color_discrete(name = "Study Species") +
  geom_encircle(expand = 0, aes(fill = study_species), alpha = 0.3) + 
  scale_fill_discrete(name = "Study Species") +
  labs(x = labs_species_bray[1], y = labs_species_bray[2]) +
  theme_bw() +
  ggtitle("Bray-Curtis Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
PCoA_species_bray

# Hihi only
mtd_tbl_short_noMCnoNTC_hihi <- mtd_tbl_short_noMCnoNTC %>%
  filter(study_species=="Hihi")
mtd_tbl_short_noMCnoNTC_hihi_noHSNA <- mtd_tbl_short_noMCnoNTC_hihi %>%
  filter(!is.na(hatched.unhatched))

bray_mat_hihi <- read_tsv("hihi_exported_diversity_metrics/bray_curtis_distance_matrix/distance-matrix.tsv")
jaccard_mat_hihi <- read_tsv("hihi_exported_diversity_metrics/jaccard_distance_matrix/distance-matrix.tsv")
unweighted_mat_hihi <- read_tsv("hihi_exported_diversity_metrics/unweighted_unifrac_distance_matrix/distance-matrix.tsv")
weighted_mat_hihi <- read_tsv("hihi_exported_diversity_metrics/weighted_unifrac_distance_matrix/distance-matrix.tsv")
colnames(bray_mat_hihi)[1] <- "SampleID"
colnames(jaccard_mat_hihi)[1] <- "SampleID"
colnames(unweighted_mat_hihi)[1] <- "SampleID"
colnames(weighted_mat_hihi)[1] <- "SampleID"
bray_mat_hihi_tbl <- as_tibble(bray_mat_hihi)
jaccard_mat_hihi_tbl <- as_tibble(jaccard_mat_hihi)
unweighted_mat_hihi_tbl <- as_tibble(unweighted_mat_hihi)
weighted_mat_hihi_tbl <- as_tibble(weighted_mat_hihi)

#dist_matrix_hihi <- weighted_mat_hihi %>%
#  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
#  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by =c("b" = "SampleID")) %>%
#  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by = "SampleID") %>%
#  select(SampleID, b, distances) %>%
#  pivot_wider(names_from = "b" , values_from = "distances") %>%
#  select(-SampleID) %>%
#  as.dist()

dist_matrix_hihi_weighted <- weighted_mat_hihi %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_hihi_unweighted <- unweighted_mat_hihi %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_hihi_jaccard <- jaccard_mat_hihi %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_hihi_bray <- bray_mat_hihi %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()

dist_matrix_hihi_noHSNA_weighted <- weighted_mat_hihi %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi_noHSNA, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi_noHSNA, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_hihi_noHSNA_unweighted <- unweighted_mat_hihi %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi_noHSNA, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi_noHSNA, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_hihi_noHSNA_jaccard <- jaccard_mat_hihi %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi_noHSNA, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi_noHSNA, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_hihi_noHSNA_bray <- bray_mat_hihi %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi_noHSNA, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi_noHSNA, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()

# PCoA
#pcoa_hihi <- cmdscale(dist_matrix_hihi, eig = TRUE, add= TRUE)
#positions_hihi <- pcoa_hihi$points
#colnames(positions_hihi) <- c("pcoa1","pcoa2")

#percent_explained_hihi <- 100 * pcoa_hihi$eig / sum(pcoa_hihi$eig)
#pretty_pe_hihi <- format(round(percent_explained_hihi[1:2], digits = 1), nsmall = 1, trim = TRUE)

#labs_hihi <- c(glue("PCo 1 ({pretty_pe_hihi[1]}%)"),
#               glue("PCo 2 ({pretty_pe_hihi[2]}%)"))
#positions_hihi %>%
#  as_tibble(rownames = "SampleID") %>%
#  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by = "SampleID") %>%
#  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
#  geom_point() +
#  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.3) + 
  #  geom_polygon(aes(fill = female_combo), alpha = 0.3) +
#  labs(x = labs_hihi[1], y = labs_hihi[2]) +
#  theme_bw()

pcoa_hihi_weighted <- cmdscale(dist_matrix_hihi_weighted, eig = TRUE, add= TRUE)
pcoa_hihi_unweighted <- cmdscale(dist_matrix_hihi_unweighted, eig = TRUE, add= TRUE)
pcoa_hihi_jaccard <- cmdscale(dist_matrix_hihi_jaccard, eig = TRUE, add= TRUE)
pcoa_hihi_bray <- cmdscale(dist_matrix_hihi_bray, eig = TRUE, add= TRUE)
pcoa_hihi_noHSNA_weighted <- cmdscale(dist_matrix_hihi_noHSNA_weighted, eig = TRUE, add= TRUE)
pcoa_hihi_noHSNA_unweighted <- cmdscale(dist_matrix_hihi_noHSNA_unweighted, eig = TRUE, add= TRUE)
pcoa_hihi_noHSNA_jaccard <- cmdscale(dist_matrix_hihi_noHSNA_jaccard, eig = TRUE, add= TRUE)
pcoa_hihi_noHSNA_bray <- cmdscale(dist_matrix_hihi_noHSNA_bray, eig = TRUE, add= TRUE)

positions_hihi_weighted <- pcoa_hihi_weighted$points
positions_hihi_unweighted <- pcoa_hihi_unweighted$points
positions_hihi_jaccard <- pcoa_hihi_jaccard$points
positions_hihi_bray <- pcoa_hihi_bray$points
positions_hihi_noHSNA_weighted <- pcoa_hihi_noHSNA_weighted$points
positions_hihi_noHSNA_unweighted <- pcoa_hihi_noHSNA_unweighted$points
positions_hihi_noHSNA_jaccard <- pcoa_hihi_noHSNA_jaccard$points
positions_hihi_noHSNA_bray <- pcoa_hihi_noHSNA_bray$points

colnames(positions_hihi_weighted) <- c("pcoa1","pcoa2")
colnames(positions_hihi_unweighted) <- c("pcoa1","pcoa2")
colnames(positions_hihi_jaccard) <- c("pcoa1","pcoa2")
colnames(positions_hihi_bray) <- c("pcoa1","pcoa2")
colnames(positions_hihi_noHSNA_weighted) <- c("pcoa1","pcoa2")
colnames(positions_hihi_noHSNA_unweighted) <- c("pcoa1","pcoa2")
colnames(positions_hihi_noHSNA_jaccard) <- c("pcoa1","pcoa2")
colnames(positions_hihi_noHSNA_bray) <- c("pcoa1","pcoa2")

percent_explained_hihi_weighted <- 100 * pcoa_hihi_weighted$eig / sum(pcoa_hihi_weighted$eig)
percent_explained_hihi_unweighted <- 100 * pcoa_hihi_unweighted$eig / sum(pcoa_hihi_unweighted$eig)
percent_explained_hihi_jaccard <- 100 * pcoa_hihi_jaccard$eig / sum(pcoa_hihi_jaccard$eig)
percent_explained_hihi_bray <- 100 * pcoa_hihi_bray$eig / sum(pcoa_hihi_bray$eig)
percent_explained_hihi_noHSNA_weighted <- 100 * pcoa_hihi_noHSNA_weighted$eig / sum(pcoa_hihi_noHSNA_weighted$eig)
percent_explained_hihi_noHSNA_unweighted <- 100 * pcoa_hihi_noHSNA_unweighted$eig / sum(pcoa_hihi_noHSNA_unweighted$eig)
percent_explained_hihi_noHSNA_jaccard <- 100 * pcoa_hihi_noHSNA_jaccard$eig / sum(pcoa_hihi_noHSNA_jaccard$eig)
percent_explained_hihi_noHSNA_bray <- 100 * pcoa_hihi_noHSNA_bray$eig / sum(pcoa_hihi_noHSNA_bray$eig)

pretty_pe_hihi_weighted <- format(round(percent_explained_hihi_weighted[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_hihi_unweighted <- format(round(percent_explained_hihi_unweighted[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_hihi_jaccard <- format(round(percent_explained_hihi_jaccard[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_hihi_bray <- format(round(percent_explained_hihi_bray[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_hihi_noHSNA_weighted <- format(round(percent_explained_hihi_noHSNA_weighted[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_hihi_noHSNA_unweighted <- format(round(percent_explained_hihi_noHSNA_unweighted[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_hihi_noHSNA_jaccard <- format(round(percent_explained_hihi_noHSNA_jaccard[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_hihi_noHSNA_bray <- format(round(percent_explained_hihi_noHSNA_bray[1:2], digits = 1), nsmall = 1, trim = TRUE)

labs_hihi_weighted <- c(glue("PCo 1 ({pretty_pe_hihi_weighted[1]}%)"),
               glue("PCo 2 ({pretty_pe_hihi_weighted[2]}%)"))
labs_hihi_unweighted <- c(glue("PCo 1 ({pretty_pe_hihi_unweighted[1]}%)"),
                      glue("PCo 2 ({pretty_pe_hihi_unweighted[2]}%)"))
labs_hihi_jaccard <- c(glue("PCo 1 ({pretty_pe_hihi_jaccard[1]}%)"),
                      glue("PCo 2 ({pretty_pe_hihi_jaccard[2]}%)"))
labs_hihi_bray <- c(glue("PCo 1 ({pretty_pe_hihi_bray[1]}%)"),
                      glue("PCo 2 ({pretty_pe_hihi_bray[2]}%)"))
labs_hihi_noHSNA_weighted <- c(glue("PCo 1 ({pretty_pe_hihi_noHSNA_weighted[1]}%)"),
                               glue("PCo 2 ({pretty_pe_hihi_noHSNA_weighted[2]}%)"))
labs_hihi_noHSNA_unweighted <- c(glue("PCo 1 ({pretty_pe_hihi_noHSNA_unweighted[1]}%)"),
                                 glue("PCo 2 ({pretty_pe_hihi_noHSNA_unweighted[2]}%)"))
labs_hihi_noHSNA_jaccard <- c(glue("PCo 1 ({pretty_pe_hihi_noHSNA_jaccard[1]}%)"),
                              glue("PCo 2 ({pretty_pe_hihi_noHSNA_jaccard[2]}%)"))
labs_hihi_noHSNA_bray <- c(glue("PCo 1 ({pretty_pe_hihi_noHSNA_bray[1]}%)"),
                           glue("PCo 2 ({pretty_pe_hihi_noHSNA_bray[2]}%)"))

PCoA_hihi_hatch_weighted <- positions_hihi_noHSNA_weighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi_noHSNA, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = hatched.unhatched )) +
  geom_point() +
  scale_color_discrete(name = "Hatching Success") +
  geom_encircle(expand = 0, aes(fill = hatched.unhatched), alpha = 0.3) + 
  scale_fill_discrete(name = "Hatching Success") +
  labs(x = labs_hihi_noHSNA_weighted[1], y = labs_hihi_noHSNA_weighted[2]) +
  theme_bw() +
  ggtitle("Weighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
PCoA_hihi_hatch_unweighted <- positions_hihi_noHSNA_unweighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi_noHSNA, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = hatched.unhatched )) +
  geom_point() +
  scale_color_discrete(name = "Hatching Success") +
  geom_encircle(expand = 0, aes(fill = hatched.unhatched), alpha = 0.3) + 
  scale_fill_discrete(name = "Hatching Success") +
  labs(x = labs_hihi_noHSNA_unweighted[1], y = labs_hihi_noHSNA_unweighted[2]) +
  theme_bw() +
  ggtitle("Unweighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
PCoA_hihi_hatch_jaccard <- positions_hihi_noHSNA_jaccard %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi_noHSNA, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = hatched.unhatched )) +
  geom_point() +
  scale_color_discrete(name = "Hatching Success") +
  geom_encircle(expand = 0, aes(fill = hatched.unhatched), alpha = 0.3) + 
  scale_fill_discrete(name = "Hatching Success") +
  labs(x = labs_hihi_noHSNA_jaccard[1], y = labs_hihi_noHSNA_jaccard[2]) +
  theme_bw() +
  ggtitle("Jaccard Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
PCoA_hihi_hatch_bray <- positions_hihi_noHSNA_bray %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi_noHSNA, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = hatched.unhatched )) +
  geom_point() +
  scale_color_discrete(name = "Hatching Success") +
  geom_encircle(expand = 0, aes(fill = hatched.unhatched), alpha = 0.3) + 
  scale_fill_discrete(name = "Hatching Success") +
  labs(x = labs_hihi_noHSNA_bray[1], y = labs_hihi_noHSNA_bray[2]) +
  theme_bw() +
  ggtitle("Bray-Curtis Distance") +
  facet_grid(row = vars("Hihi"), scale = "free_y", space = "free") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text.y = element_text(size = 14, face="bold"))

PCoA_hihi_hatch_weighted
PCoA_hihi_hatch_unweighted
PCoA_hihi_hatch_jaccard
PCoA_hihi_hatch_bray

all_hihi_hatch_PCoA <- (PCoA_hihi_hatch_weighted | PCoA_hihi_hatch_unweighted | PCoA_hihi_hatch_jaccard | PCoA_hihi_hatch_bray) + plot_layout(guides = 'collect')
ggview(all_hihi_hatch_PCoA, units = "px", height = 2000, width = 5000)

PCoA_hihi_female_weighted <- positions_hihi_weighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
  geom_point() +
  scale_color_discrete(name = "Hihi Dam ID") +
  geom_encircle(expand = 0, aes(color = female_combo), alpha = 1) + 
  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.2) + 
  scale_fill_discrete(name = "Hihi Dam ID") +
  labs(x = labs_hihi_weighted[1], y = labs_hihi_weighted[2]) +
  theme_bw() +
  ggtitle("Weighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)) +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_hihi_female_unweighted <- positions_hihi_unweighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
  geom_point() +
  scale_color_discrete(name = "Hihi Dam ID") +
  geom_encircle(expand = 0, aes(color = female_combo), alpha = 1) + 
  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.2) + 
  scale_fill_discrete(name = "Hihi Dam ID") +
  labs(x = labs_hihi_unweighted[1], y = labs_hihi_unweighted[2]) +
  theme_bw() +
  ggtitle("Unweighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_hihi_female_jaccard <- positions_hihi_jaccard %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
  geom_point() +
  scale_color_discrete(name = "Hihi Dam ID") +
  geom_encircle(expand = 0, aes(color = female_combo), alpha = 1) + 
  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.2) + 
  scale_fill_discrete(name = "Hihi Dam ID") +
  labs(x = labs_hihi_jaccard[1], y = labs_hihi_jaccard[2]) +
  theme_bw() +
  ggtitle("Jaccard Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_hihi_female_bray <- positions_hihi_bray %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_hihi, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
  geom_point() +
  scale_color_discrete(name = "Hihi Dam ID") +
  geom_encircle(expand = 0, aes(color = female_combo), alpha = 1) + 
  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.2) + 
  scale_fill_discrete(name = "Hihi Dam ID") +
  labs(x = labs_hihi_bray[1], y = labs_hihi_bray[2]) +
  theme_bw() +
  ggtitle("Bray-Curtis Distance") +
#  facet_grid(row = vars("Hihi"), scale = "free_y", space = "free") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text.y = element_text(size = 14, face="bold")) +
  guides(color=guide_legend(ncol=2, byrow=FALSE))

PCoA_hihi_female_weighted
PCoA_hihi_female_unweighted
PCoA_hihi_female_jaccard
PCoA_hihi_female_bray

all_hihi_female_PCoA <- (PCoA_hihi_female_weighted | PCoA_hihi_female_unweighted | PCoA_hihi_female_jaccard | PCoA_hihi_female_bray) + plot_layout(guides = 'collect')
ggview(all_hihi_female_PCoA, units = "px", height = 2000, width = 5000)
all_hihi_female_PCoA

# Kakī only
mtd_tbl_short_noMCnoNTC_kaki <- mtd_tbl_short_noMCnoNTC %>%
  filter(study_species=="Kakī")
mtd_tbl_short_noMCnoNTC_kaki_repeats <- filter(mtd_tbl_short_noMCnoNTC_kaki, female_combo %in% c("GBKG/RBK","BKOBK/WY","BKBKW/BKBK","BKYO/GW","BKYO/GBK","BKBKY/WBK"))

bray_mat_kaki <- read_tsv("kaki_exported_diversity_metrics/bray_curtis_distance_matrix/distance-matrix.tsv")
jaccard_mat_kaki <- read_tsv("kaki_exported_diversity_metrics/jaccard_distance_matrix/distance-matrix.tsv")
unweighted_mat_kaki <- read_tsv("kaki_exported_diversity_metrics/unweighted_unifrac_distance_matrix/distance-matrix.tsv")
weighted_mat_kaki <- read_tsv("kaki_exported_diversity_metrics/weighted_unifrac_distance_matrix/distance-matrix.tsv")
colnames(bray_mat_kaki)[1] <- "SampleID"
colnames(jaccard_mat_kaki)[1] <- "SampleID"
colnames(unweighted_mat_kaki)[1] <- "SampleID"
colnames(weighted_mat_kaki)[1] <- "SampleID"
bray_mat_kaki_tbl <- as_tibble(bray_mat_kaki)
jaccard_mat_kaki_tbl <- as_tibble(jaccard_mat_kaki)
unweighted_mat_kaki_tbl <- as_tibble(unweighted_mat_kaki)
weighted_mat_kaki_tbl <- as_tibble(weighted_mat_kaki)

dist_matrix_kaki_weighted <- weighted_mat_kaki %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_kaki_unweighted <- unweighted_mat_kaki %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_kaki_jaccard <- jaccard_mat_kaki %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_kaki_bray <- bray_mat_kaki %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()

# PCoA
pcoa_kaki_weighted <- cmdscale(dist_matrix_kaki_weighted, eig = TRUE, add= TRUE)
pcoa_kaki_unweighted <- cmdscale(dist_matrix_kaki_unweighted, eig = TRUE, add= TRUE)
pcoa_kaki_jaccard <- cmdscale(dist_matrix_kaki_jaccard, eig = TRUE, add= TRUE)
pcoa_kaki_bray <- cmdscale(dist_matrix_kaki_bray, eig = TRUE, add= TRUE)

positions_kaki_weighted <- pcoa_kaki_weighted$points
positions_kaki_unweighted <- pcoa_kaki_unweighted$points
positions_kaki_jaccard <- pcoa_kaki_jaccard$points
positions_kaki_bray <- pcoa_kaki_bray$points

colnames(positions_kaki_weighted) <- c("pcoa1","pcoa2")
colnames(positions_kaki_unweighted) <- c("pcoa1","pcoa2")
colnames(positions_kaki_jaccard) <- c("pcoa1","pcoa2")
colnames(positions_kaki_bray) <- c("pcoa1","pcoa2")

percent_explained_kaki_weighted <- 100 * pcoa_kaki_weighted$eig / sum(pcoa_kaki_weighted$eig)
percent_explained_kaki_unweighted <- 100 * pcoa_kaki_unweighted$eig / sum(pcoa_kaki_unweighted$eig)
percent_explained_kaki_jaccard <- 100 * pcoa_kaki_jaccard$eig / sum(pcoa_kaki_jaccard$eig)
percent_explained_kaki_bray <- 100 * pcoa_kaki_bray$eig / sum(pcoa_kaki_bray$eig)

pretty_pe_kaki_weighted <- format(round(percent_explained_kaki_weighted[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_kaki_unweighted <- format(round(percent_explained_kaki_unweighted[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_kaki_jaccard <- format(round(percent_explained_kaki_jaccard[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_kaki_bray <- format(round(percent_explained_kaki_bray[1:2], digits = 1), nsmall = 1, trim = TRUE)

labs_kaki_weighted <- c(glue("PCo 1 ({pretty_pe_kaki_weighted[1]}%)"),
               glue("PCo 2 ({pretty_pe_kaki_weighted[2]}%)"))
labs_kaki_unweighted <- c(glue("PCo 1 ({pretty_pe_kaki_unweighted[1]}%)"),
               glue("PCo 2 ({pretty_pe_kaki_unweighted[2]}%)"))
labs_kaki_jaccard <- c(glue("PCo 1 ({pretty_pe_kaki_jaccard[1]}%)"),
               glue("PCo 2 ({pretty_pe_kaki_jaccard[2]}%)"))
labs_kaki_bray <- c(glue("PCo 1 ({pretty_pe_kaki_bray[1]}%)"),
               glue("PCo 2 ({pretty_pe_kaki_bray[2]}%)"))

PCoA_kaki_hatch_weighted <- positions_kaki_weighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = hatched.unhatched )) +
  geom_point() +
  scale_color_discrete(name = "Hatching Success") +
  geom_encircle(expand = 0, aes(fill = hatched.unhatched), alpha = 0.3) + 
  scale_fill_discrete(name = "Hatching Success") +  labs(x = labs_kaki_weighted[1], y = labs_kaki_weighted[2]) +
  theme_bw() +
  ggtitle("Weighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
PCoA_kaki_hatch_unweighted <- positions_kaki_unweighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = hatched.unhatched )) +
  geom_point() +
  scale_color_discrete(name = "Hatching Success") +
  geom_encircle(expand = 0, aes(fill = hatched.unhatched), alpha = 0.3) + 
  scale_fill_discrete(name = "Hatching Success") +  labs(x = labs_kaki_unweighted[1], y = labs_kaki_unweighted[2]) +
  theme_bw() +
  ggtitle("Unweighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
PCoA_kaki_hatch_jaccard <- positions_kaki_jaccard %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = hatched.unhatched )) +
  geom_point() +
  scale_color_discrete(name = "Hatching Success") +
  geom_encircle(expand = 0, aes(fill = hatched.unhatched), alpha = 0.3) + 
  scale_fill_discrete(name = "Hatching Success") +  labs(x = labs_kaki_jaccard[1], y = labs_kaki_jaccard[2]) +
  theme_bw() +
  ggtitle("Jaccard Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
PCoA_kaki_hatch_bray <- positions_kaki_bray %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = hatched.unhatched )) +
  geom_point() +
  scale_color_discrete(name = "Hatching Success") +
  geom_encircle(expand = 0, aes(fill = hatched.unhatched), alpha = 0.3) + 
  scale_fill_discrete(name = "Hatching Success") +  
  labs(x = labs_kaki_bray[1], y = labs_kaki_bray[2]) +
  theme_bw() +
  ggtitle("Bray-Curtis Distance") +
  facet_grid(row = vars("Kakī"), scale = "free_y", space = "free") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text.y = element_text(size = 14, face="bold"))

PCoA_kaki_hatch_weighted
PCoA_kaki_hatch_unweighted
PCoA_kaki_hatch_jaccard
PCoA_kaki_hatch_bray

all_kaki_hatch_PCoA <- (PCoA_kaki_hatch_weighted | PCoA_kaki_hatch_unweighted | PCoA_kaki_hatch_jaccard | PCoA_kaki_hatch_bray) + plot_layout(guides = 'collect')
ggview(all_kaki_hatch_PCoA, units = "px", height = 2000, width = 5000)

all_hatch_PCoA <- (PCoA_hihi_hatch_weighted | PCoA_hihi_hatch_unweighted | PCoA_hihi_hatch_jaccard | PCoA_hihi_hatch_bray) / 
  (PCoA_kaki_hatch_weighted | PCoA_kaki_hatch_unweighted | PCoA_kaki_hatch_jaccard | PCoA_kaki_hatch_bray) + plot_layout(guides = 'collect')
ggview(all_hatch_PCoA, units = "px", height = 4000, width = 5000)

PCoA_kaki_female_weighted <- positions_kaki_weighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Dam ID") +
  geom_encircle(expand = 0, aes(color = female_combo), alpha = 1) + 
  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.2) + 
  scale_fill_discrete(name = "Kakī Dam ID") +  labs(x = labs_kaki_weighted[1], y = labs_kaki_weighted[2]) +
  theme_bw() +
  ggtitle("Weighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.justification = "left") +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_kaki_female_unweighted <- positions_kaki_unweighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Dam ID") +
  geom_encircle(expand = 0, aes(color = female_combo), alpha = 1) + 
  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.2) + 
  scale_fill_discrete(name = "Kakī Dam ID") +  labs(x = labs_kaki_unweighted[1], y = labs_kaki_unweighted[2]) +
  theme_bw() +
  ggtitle("Unweighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.justification = "left") +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_kaki_female_jaccard <- positions_kaki_jaccard %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Dam ID") +
  geom_encircle(expand = 0, aes(color = female_combo), alpha = 1) + 
  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.2) + 
  scale_fill_discrete(name = "Kakī Dam ID") +  labs(x = labs_kaki_jaccard[1], y = labs_kaki_jaccard[2]) +
  theme_bw() +
  ggtitle("Jaccard Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.justification = "left") +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_kaki_female_bray <- positions_kaki_bray %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Dam ID") +
  geom_encircle(expand = 0, aes(color = female_combo), alpha = 1) + 
  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.2) + 
  scale_fill_discrete(name = "Kakī Dam ID") +  
  labs(x = labs_kaki_bray[1], y = labs_kaki_bray[2]) +
  theme_bw() +
  ggtitle("Bray-Curtis Distance") +
#  facet_grid(row = vars("Kakī"), scale = "free_y", space = "free") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.justification = "left",
        strip.text.y = element_text(size = 14, face="bold")) + 
  guides(color=guide_legend(ncol=2, byrow=FALSE))

PCoA_kaki_female_weighted
PCoA_kaki_female_unweighted
PCoA_kaki_female_jaccard
PCoA_kaki_female_bray

all_kaki_female_PCoA <- (PCoA_kaki_female_weighted | PCoA_kaki_female_unweighted | PCoA_kaki_female_jaccard | PCoA_kaki_female_bray) + plot_layout(guides = 'collect')
ggview(all_kaki_female_PCoA, units = "px", height = 2000, width = 5000)

all_female_PCoA <- (PCoA_hihi_female_weighted | PCoA_hihi_female_unweighted | PCoA_hihi_female_jaccard | PCoA_hihi_female_bray) / 
  (PCoA_kaki_female_weighted | PCoA_kaki_female_unweighted | PCoA_kaki_female_jaccard | PCoA_kaki_female_bray) + plot_layout(guides = 'collect')
ggview(all_female_PCoA, units = "px", height = 4000, width = 5000)
all_female_PCoA

PCoA_kaki_clutchID_weighted <- positions_kaki_weighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = clutch.nest_id_corrected)) +
  geom_point() +
  scale_color_discrete(name = "Kakī Clutch ID") +
  geom_encircle(expand = 0, aes(color = clutch.nest_id_corrected), alpha = 1) + 
  geom_encircle(expand = 0, aes(fill = clutch.nest_id_corrected), alpha = 0.2) + 
  scale_fill_discrete(name = "Kakī Clutch ID") +  labs(x = labs_kaki_weighted[1], y = labs_kaki_weighted[2]) +
  theme_bw() +
  ggtitle("Weighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.justification = "left") +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_kaki_clutchID_unweighted <- positions_kaki_unweighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = clutch.nest_id_corrected )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Clutch ID") +
  geom_encircle(expand = 0, aes(color = clutch.nest_id_corrected), alpha = 1) + 
  geom_encircle(expand = 0, aes(fill = clutch.nest_id_corrected), alpha = 0.2) + 
  scale_fill_discrete(name = "Kakī Clutch ID") +  labs(x = labs_kaki_unweighted[1], y = labs_kaki_unweighted[2]) +
  theme_bw() +
  ggtitle("Unweighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.justification = "left") +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_kaki_clutchID_jaccard <- positions_kaki_jaccard %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = clutch.nest_id_corrected )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Clutch ID") +
  geom_encircle(expand = 0, aes(color = clutch.nest_id_corrected), alpha = 1) + 
  geom_encircle(expand = 0, aes(fill = clutch.nest_id_corrected), alpha = 0.2) + 
  scale_fill_discrete(name = "Kakī Clutch ID") +  labs(x = labs_kaki_jaccard[1], y = labs_kaki_jaccard[2]) +
  theme_bw() +
  ggtitle("Jaccard Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.justification = "left") +
  guides(color=guide_legend(ncol=2, byrow=FALSE))

#labs_text <- positions_kaki_bray %>%
#  as_tibble(rownames = "SampleID") %>%
#  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>% 
#  group_by(clutch.nest_id_corrected) %>%
#  filter(pcoa2 == max(pcoa2))
PCoA_kaki_clutchID_bray <- positions_kaki_bray %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = clutch.nest_id_corrected )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Clutch ID") +
  geom_encircle(expand = 0, aes(color = clutch.nest_id_corrected), alpha = 1) + 
  geom_encircle(expand = 0, aes(fill = clutch.nest_id_corrected), alpha = 0.2) + 
  scale_fill_discrete(name = "Kakī Clutch ID") +  
  labs(x = labs_kaki_bray[1], y = labs_kaki_bray[2]) +
  theme_bw() +
  ggtitle("Bray-Curtis Distance") +
#  geom_text(data = labs_text, 
#            aes(pcoa1, pcoa2, label = clutch.nest_id_corrected)) +
#  facet_grid(row = vars("Kakī"), scale = "free_y", space = "free") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.justification = "left",
        strip.text.y = element_text(size = 14, face="bold")) + 
  guides(color=guide_legend(ncol=2, byrow=FALSE))

PCoA_kaki_clutchID_weighted
PCoA_kaki_clutchID_unweighted
PCoA_kaki_clutchID_jaccard
PCoA_kaki_clutchID_bray

all_kaki_clutchID_PCoA <- (PCoA_kaki_clutchID_weighted | PCoA_kaki_clutchID_unweighted | PCoA_kaki_clutchID_jaccard | PCoA_kaki_clutchID_bray) + plot_layout(guides = 'collect')
ggview(all_kaki_clutchID_PCoA, units = "px", height = 2000, width = 5000)

all_kaki_female.clutchID_PCoA <- ((PCoA_kaki_female_weighted | PCoA_kaki_female_unweighted | PCoA_kaki_female_jaccard | PCoA_kaki_female_bray) + plot_layout(guides = 'collect')) / 
  ((PCoA_kaki_clutchID_weighted | PCoA_kaki_clutchID_unweighted | PCoA_kaki_clutchID_jaccard | PCoA_kaki_clutchID_bray) + plot_layout(guides = 'collect')) & theme(legend.justification = "left")
ggview(all_kaki_female.clutchID_PCoA, units = "px", height = 4000, width = 5000)
all_kaki_female.clutchID_PCoA

PCoA_kaki_wildcap_weighted <- positions_kaki_weighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = wild.captive )) +
  geom_point() +
#  scale_color_discrete(name = "Egg Captive- or Wild-laid") +
  scale_color_manual(values=c("limegreen","purple"), name = "Egg Captive- or Wild-laid") +
  geom_encircle(expand = 0, aes(fill = wild.captive), alpha = 0.3) + 
#  scale_fill_discrete(name = "Egg Captive- or Wild-laid") +  
  scale_fill_manual(values=c("limegreen","purple"),labels=c("Captive","Wild"),name = "Egg Captive- or Wild-laid") +
  labs(x = labs_kaki_weighted[1], y = labs_kaki_weighted[2]) +
  theme_bw() +
  ggtitle("Weighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
PCoA_kaki_wildcap_unweighted <- positions_kaki_unweighted %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = wild.captive )) +
  geom_point() +
  #  scale_color_discrete(name = "Egg Captive- or Wild-laid") +
  scale_color_manual(values=c("limegreen","purple"), name = "Egg Captive- or Wild-laid") +
  geom_encircle(expand = 0, aes(fill = wild.captive), alpha = 0.3) + 
  #  scale_fill_discrete(name = "Egg Captive- or Wild-laid") +  
  scale_fill_manual(values=c("limegreen","purple"),labels=c("Captive","Wild"),name = "Egg Captive- or Wild-laid") +
  labs(x = labs_kaki_unweighted[1], y = labs_kaki_unweighted[2]) +
  theme_bw() +
  ggtitle("Unweighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
PCoA_kaki_wildcap_jaccard <- positions_kaki_jaccard %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = wild.captive )) +
  geom_point() +
  #  scale_color_discrete(name = "Egg Captive- or Wild-laid") +
  scale_color_manual(values=c("limegreen","purple"), name = "Egg Captive- or Wild-laid") +
  geom_encircle(expand = 0, aes(fill = wild.captive), alpha = 0.3) + 
  #  scale_fill_discrete(name = "Egg Captive- or Wild-laid") +  
  scale_fill_manual(values=c("limegreen","purple"),labels=c("Captive","Wild"),name = "Egg Captive- or Wild-laid") + 
  labs(x = labs_kaki_jaccard[1], y = labs_kaki_jaccard[2]) +
  theme_bw() +
  ggtitle("Jaccard Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
PCoA_kaki_wildcap_bray <- positions_kaki_bray %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = wild.captive )) +
  geom_point() +
  #  scale_color_discrete(name = "Egg Captive- or Wild-laid") +
  scale_color_manual(values=c("limegreen","purple"), name = "Egg Captive- or Wild-laid") +
  geom_encircle(expand = 0, aes(fill = wild.captive), alpha = 0.3) + 
  #  scale_fill_discrete(name = "Egg Captive- or Wild-laid") +  
  scale_fill_manual(values=c("limegreen","purple"),labels=c("Captive","Wild"),name = "Egg Captive- or Wild-laid") +
  labs(x = labs_kaki_bray[1], y = labs_kaki_bray[2]) +
  theme_bw() +
  ggtitle("Bray-Curtis Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

all_kaki_wildcap_PCoA <- (PCoA_kaki_wildcap_weighted | PCoA_kaki_wildcap_unweighted | PCoA_kaki_wildcap_jaccard | PCoA_kaki_wildcap_bray) + plot_layout(guides = 'collect')
ggview(all_kaki_wildcap_PCoA, units = "px", height = 2000, width = 5000)

dist_matrix_kaki_weighted_repeats <- weighted_mat_kaki %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_kaki_unweighted_repeats <- unweighted_mat_kaki %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_kaki_jaccard_repeats <- jaccard_mat_kaki %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()
dist_matrix_kaki_bray_repeats <- bray_mat_kaki %>%
  pivot_longer(cols = -SampleID, names_to = "b", values_to ="distances") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by =c("b" = "SampleID")) %>%
  select(SampleID, b, distances) %>%
  pivot_wider(names_from = "b" , values_from = "distances") %>%
  select(-SampleID) %>%
  as.dist()

# PCoA
pcoa_kaki_weighted_repeats <- cmdscale(dist_matrix_kaki_weighted_repeats, eig = TRUE, add= TRUE)
pcoa_kaki_unweighted_repeats <- cmdscale(dist_matrix_kaki_unweighted_repeats, eig = TRUE, add= TRUE)
pcoa_kaki_jaccard_repeats <- cmdscale(dist_matrix_kaki_jaccard_repeats, eig = TRUE, add= TRUE)
pcoa_kaki_bray_repeats <- cmdscale(dist_matrix_kaki_bray_repeats, eig = TRUE, add= TRUE)

positions_kaki_weighted_repeats <- pcoa_kaki_weighted_repeats$points
positions_kaki_unweighted_repeats <- pcoa_kaki_unweighted_repeats$points
positions_kaki_jaccard_repeats <- pcoa_kaki_jaccard_repeats$points
positions_kaki_bray_repeats <- pcoa_kaki_bray_repeats$points

colnames(positions_kaki_weighted_repeats) <- c("pcoa1","pcoa2")
colnames(positions_kaki_unweighted_repeats) <- c("pcoa1","pcoa2")
colnames(positions_kaki_jaccard_repeats) <- c("pcoa1","pcoa2")
colnames(positions_kaki_bray_repeats) <- c("pcoa1","pcoa2")

percent_explained_kaki_weighted_repeats <- 100 * pcoa_kaki_weighted_repeats$eig / sum(pcoa_kaki_weighted_repeats$eig)
percent_explained_kaki_unweighted_repeats <- 100 * pcoa_kaki_unweighted_repeats$eig / sum(pcoa_kaki_unweighted_repeats$eig)
percent_explained_kaki_jaccard_repeats <- 100 * pcoa_kaki_jaccard_repeats$eig / sum(pcoa_kaki_jaccard_repeats$eig)
percent_explained_kaki_bray_repeats <- 100 * pcoa_kaki_bray_repeats$eig / sum(pcoa_kaki_bray_repeats$eig)

pretty_pe_kaki_weighted_repeats <- format(round(percent_explained_kaki_weighted_repeats[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_kaki_unweighted_repeats <- format(round(percent_explained_kaki_unweighted_repeats[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_kaki_jaccard_repeats <- format(round(percent_explained_kaki_jaccard_repeats[1:2], digits = 1), nsmall = 1, trim = TRUE)
pretty_pe_kaki_bray_repeats <- format(round(percent_explained_kaki_bray_repeats[1:2], digits = 1), nsmall = 1, trim = TRUE)

labs_kaki_weighted_repeats <- c(glue("PCo 1 ({pretty_pe_kaki_weighted_repeats[1]}%)"),
                        glue("PCo 2 ({pretty_pe_kaki_weighted_repeats[2]}%)"))
labs_kaki_unweighted_repeats <- c(glue("PCo 1 ({pretty_pe_kaki_unweighted_repeats[1]}%)"),
                          glue("PCo 2 ({pretty_pe_kaki_unweighted_repeats[2]}%)"))
labs_kaki_jaccard_repeats <- c(glue("PCo 1 ({pretty_pe_kaki_jaccard_repeats[1]}%)"),
                       glue("PCo 2 ({pretty_pe_kaki_jaccard_repeats[2]}%)"))
labs_kaki_bray_repeats <- c(glue("PCo 1 ({pretty_pe_kaki_bray_repeats[1]}%)"),
                    glue("PCo 2 ({pretty_pe_kaki_bray_repeats[2]}%)"))

PCoA_kaki_female_weighted_repeats <- positions_kaki_weighted_repeats %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Dam ID") +
  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.3) + 
  scale_fill_discrete(name = "Kakī Dam ID") +  labs(x = labs_kaki_weighted_repeats[1], y = labs_kaki_weighted_repeats[2]) +
  theme_bw() +
  ggtitle("Weighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_kaki_female_unweighted_repeats <- positions_kaki_unweighted_repeats %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Dam ID") +
  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.3) + 
  scale_fill_discrete(name = "Kakī Dam ID") +  labs(x = labs_kaki_unweighted_repeats[1], y = labs_kaki_unweighted_repeats[2]) +
  theme_bw() +
  ggtitle("Unweighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_kaki_female_jaccard_repeats <- positions_kaki_jaccard_repeats %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Dam ID") +
  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.3) + 
  scale_fill_discrete(name = "Kakī Dam ID") +  labs(x = labs_kaki_jaccard_repeats[1], y = labs_kaki_jaccard_repeats[2]) +
  theme_bw() +
  ggtitle("Jaccard Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_kaki_female_bray_repeats <- positions_kaki_bray_repeats %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = female_combo )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Dam ID") +
  geom_encircle(expand = 0, aes(fill = female_combo), alpha = 0.3) + 
  scale_fill_discrete(name = "Kakī Dam ID") +  
  labs(x = labs_kaki_bray_repeats[1], y = labs_kaki_bray_repeats[2]) +
  theme_bw() +
  ggtitle("Bray-Curtis Distance") +
  facet_grid(row = vars("Kakī"), scale = "free_y", space = "free") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text.y = element_text(size = 14, face="bold")) + 
  guides(color=guide_legend(ncol=2, byrow=FALSE))

PCoA_kaki_clutchID_weighted_repeats <- positions_kaki_weighted_repeats %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = clutch.nest_id_corrected )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Dam ID") +
  geom_encircle(expand = 0, aes(fill = clutch.nest_id_corrected), alpha = 0.3) + 
  scale_fill_discrete(name = "Kakī Dam ID") +  labs(x = labs_kaki_weighted_repeats[1], y = labs_kaki_weighted_repeats[2]) +
  theme_bw() +
  ggtitle("Weighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_kaki_clutchID_unweighted_repeats <- positions_kaki_unweighted_repeats %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = clutch.nest_id_corrected )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Dam ID") +
  geom_encircle(expand = 0, aes(fill = clutch.nest_id_corrected), alpha = 0.3) + 
  scale_fill_discrete(name = "Kakī Dam ID") +  labs(x = labs_kaki_unweighted_repeats[1], y = labs_kaki_unweighted_repeats[2]) +
  theme_bw() +
  ggtitle("Unweighted UniFrac Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_kaki_clutchID_jaccard_repeats <- positions_kaki_jaccard_repeats %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = clutch.nest_id_corrected )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Dam ID") +
  geom_encircle(expand = 0, aes(fill = clutch.nest_id_corrected), alpha = 0.3) + 
  scale_fill_discrete(name = "Kakī Dam ID") +  labs(x = labs_kaki_jaccard_repeats[1], y = labs_kaki_jaccard_repeats[2]) +
  theme_bw() +
  ggtitle("Jaccard Distance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(color=guide_legend(ncol=2, byrow=FALSE))
PCoA_kaki_clutchID_bray_repeats <- positions_kaki_bray_repeats %>%
  as_tibble(rownames = "SampleID") %>%
  inner_join(., mtd_tbl_short_noMCnoNTC_kaki_repeats, by = "SampleID") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = clutch.nest_id_corrected )) +
  geom_point() +
  scale_color_discrete(name = "Kakī Dam ID") +
  geom_encircle(expand = 0, aes(fill = clutch.nest_id_corrected), alpha = 0.3) + 
  scale_fill_discrete(name = "Kakī Dam ID") +  
  labs(x = labs_kaki_bray_repeats[1], y = labs_kaki_bray_repeats[2]) +
  theme_bw() +
  ggtitle("Bray-Curtis Distance") +
  facet_grid(row = vars("Kakī"), scale = "free_y", space = "free") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text.y = element_text(size = 14, face="bold")) + 
  guides(color=guide_legend(ncol=2, byrow=FALSE))


##### Alpha Diversity Plots
# Hihi vs. kaki
shannon_all <- read_tsv("all_exported_diversity_metrics/shannon_vector/alpha-diversity.tsv")
evenness_all <- read_tsv("all_exported_diversity_metrics/evenness_vector/alpha-diversity.tsv")
faith_pd_all <- read_tsv("all_exported_diversity_metrics/faith_pd_vector/alpha-diversity.tsv")
observed_features_all <- read_tsv("all_exported_diversity_metrics/observed_features_vector/alpha-diversity.tsv")
colnames(shannon_all)[1] <- "SampleID"
colnames(evenness_all)[1] <- "SampleID"
colnames(faith_pd_all)[1] <- "SampleID"
colnames(observed_features_all)[1] <- "SampleID"

all_alpha_data <- mtd_tbl_short_noMCnoNTC %>%
  inner_join(., shannon_all, by = "SampleID") %>%
  inner_join(., evenness_all, by = "SampleID") %>%
  inner_join(., faith_pd_all, by = "SampleID") %>%
  inner_join(., observed_features_all, by = "SampleID")

all_alpha_shannon_plot <- all_alpha_data %>%
  ggplot(aes(x = study_species, y= shannon_entropy, fill = study_species, color = study_species, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab("Alpha Diversity Measure") +
  xlab(NULL) +
  scale_fill_discrete(name = "Study Species") +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Shannon's Diversity Index") +
  stat_signif(comparisons = list(c("Hihi","Kakī")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
all_alpha_shannon_plot
all_alpha_observed_plot <- all_alpha_data %>%
  ggplot(aes(x = study_species, y= observed_features, fill = study_species, color = study_species, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_discrete(name = "Study Species") +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Observed Features") +
  stat_signif(comparisons = list(c("Hihi","Kakī")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
all_alpha_observed_plot
all_alpha_faith_plot <- all_alpha_data %>%
  ggplot(aes(x = study_species, y= faith_pd, fill = study_species, color = study_species, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_discrete(name = "Study Species") +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Faith's Phylogenetic Diversity") +
  stat_signif(comparisons = list(c("Hihi","Kakī")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
all_alpha_faith_plot
all_alpha_evenness_plot <- all_alpha_data %>%
  ggplot(aes(x = study_species, y= pielou_evenness, fill = study_species, color = study_species, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_discrete(name = "Study Species") +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Pielou's Evenness") +
  stat_signif(comparisons = list(c("Hihi","Kakī")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
all_alpha_evenness_plot

all_alpha_all_plot <- (all_alpha_shannon_plot | all_alpha_observed_plot | all_alpha_faith_plot | all_alpha_evenness_plot) + plot_layout(guides = 'collect')
ggview(all_alpha_all_plot, units = "px", height = 2000, width = 5000)

all_alpha.beta_plots <- (all_alpha_shannon_plot | all_alpha_observed_plot | all_alpha_faith_plot | all_alpha_evenness_plot) / 
  (PCoA_species_weighted | PCoA_species_unweighted | PCoA_species_jaccard | PCoA_species_bray) + plot_layout(guides = 'collect')
ggview(all_alpha.beta_plots, units = "px", height = 4000, width = 5000)
all_alpha.beta_plots

# Hihi only
shannon_hihi <- read_tsv("hihi_exported_diversity_metrics/shannon_vector/alpha-diversity.tsv")
evenness_hihi <- read_tsv("hihi_exported_diversity_metrics/evenness_vector/alpha-diversity.tsv")
faith_pd_hihi <- read_tsv("hihi_exported_diversity_metrics/faith_pd_vector/alpha-diversity.tsv")
observed_features_hihi <- read_tsv("hihi_exported_diversity_metrics/observed_features_vector/alpha-diversity.tsv")
colnames(shannon_hihi)[1] <- "SampleID"
colnames(evenness_hihi)[1] <- "SampleID"
colnames(faith_pd_hihi)[1] <- "SampleID"
colnames(observed_features_hihi)[1] <- "SampleID"

hihi_alpha_data <- mtd_tbl_short_noMCnoNTC_hihi_noHSNA %>%
  inner_join(., shannon_hihi, by = "SampleID") %>%
  inner_join(., evenness_hihi, by = "SampleID") %>%
  inner_join(., faith_pd_hihi, by = "SampleID") %>%
  inner_join(., observed_features_hihi, by = "SampleID")

hihi_alpha_shannon_plot <- hihi_alpha_data %>%
  ggplot(aes(x = hatched.unhatched, y= shannon_entropy, fill = hatched.unhatched, color = hatched.unhatched, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab("Alpha Diversity Measure") +
  xlab(NULL) +
  scale_fill_discrete(name = "Hatching Success") +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Shannon's Diversity Index") +
  stat_signif(comparisons = list(c("Hatched","Unhatched")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

hihi_alpha_observed_plot <- hihi_alpha_data %>%
  ggplot(aes(x = hatched.unhatched, y= observed_features, fill = hatched.unhatched, color = hatched.unhatched, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_discrete(name = "Hatching Success") +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Observed Features") +
  stat_signif(comparisons = list(c("Hatched","Unhatched")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

hihi_alpha_faith_plot <- hihi_alpha_data %>%
  ggplot(aes(x = hatched.unhatched, y= faith_pd, fill = hatched.unhatched, color = hatched.unhatched, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_discrete(name = "Hatching Success") +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Faith's Phylogenetic Diversity") +
  stat_signif(comparisons = list(c("Hatched","Unhatched")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

hihi_alpha_evenness_plot <- hihi_alpha_data %>%
  ggplot(aes(x = hatched.unhatched, y= pielou_evenness, fill = hatched.unhatched, color = hatched.unhatched, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_discrete(name = "Hatching Success") +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Pielou's Evenness") +
  stat_signif(comparisons = list(c("Hatched","Unhatched")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  facet_grid(row = vars("Hihi"), scale = "free_y", space = "free") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.y = element_text(size = 14, face="bold"))

hihi_alpha_all_plot <- (hihi_alpha_shannon_plot | hihi_alpha_observed_plot | hihi_alpha_faith_plot | hihi_alpha_evenness_plot) + plot_layout(guides = 'collect')
ggview(hihi_alpha_all_plot, units = "px", height = 2000, width = 5000)

# Kakī only
shannon_kaki <- read_tsv("kaki_exported_diversity_metrics/shannon_vector/alpha-diversity.tsv")
evenness_kaki <- read_tsv("kaki_exported_diversity_metrics/evenness_vector/alpha-diversity.tsv")
faith_pd_kaki <- read_tsv("kaki_exported_diversity_metrics/faith_pd_vector/alpha-diversity.tsv")
observed_features_kaki <- read_tsv("kaki_exported_diversity_metrics/observed_features_vector/alpha-diversity.tsv")
colnames(shannon_kaki)[1] <- "SampleID"
colnames(evenness_kaki)[1] <- "SampleID"
colnames(faith_pd_kaki)[1] <- "SampleID"
colnames(observed_features_kaki)[1] <- "SampleID"

kaki_alpha_data <- mtd_tbl_short_noMCnoNTC_kaki %>%
  inner_join(., shannon_kaki, by = "SampleID") %>%
  inner_join(., evenness_kaki, by = "SampleID") %>%
  inner_join(., faith_pd_kaki, by = "SampleID") %>%
  inner_join(., observed_features_kaki, by = "SampleID")

kaki_alpha_shannon_plot <- kaki_alpha_data %>%
  ggplot(aes(x = hatched.unhatched, y= shannon_entropy, fill = hatched.unhatched, color = hatched.unhatched, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab("Alpha Diversity Measure") +
  xlab(NULL) +
  scale_fill_discrete(name = "Hatching Success") +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Shannon's Diversity Index") +
  stat_signif(comparisons = list(c("Hatched","Unhatched")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

kaki_alpha_observed_plot <- kaki_alpha_data %>%
  ggplot(aes(x = hatched.unhatched, y= observed_features, fill = hatched.unhatched, color = hatched.unhatched, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_discrete(name = "Hatching Success") +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Observed Features") +
  stat_signif(comparisons = list(c("Hatched","Unhatched")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

kaki_alpha_faith_plot <- kaki_alpha_data %>%
  ggplot(aes(x = hatched.unhatched, y= faith_pd, fill = hatched.unhatched, color = hatched.unhatched, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_discrete(name = "Hatching Success") +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Faith's Phylogenetic Diversity") +
  stat_signif(comparisons = list(c("Hatched","Unhatched")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

kaki_alpha_evenness_plot <- kaki_alpha_data %>%
  ggplot(aes(x = hatched.unhatched, y= pielou_evenness, fill = hatched.unhatched, color = hatched.unhatched, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_discrete(name = "Hatching Success") +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Pielou's Evenness")  +
  stat_signif(comparisons = list(c("Hatched","Unhatched")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  facet_grid(row = vars("Kakī"), scale = "free_y", space = "free") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.y = element_text(size = 14, face="bold"))

kaki_alpha_all_plot <- (kaki_alpha_shannon_plot | kaki_alpha_observed_plot | kaki_alpha_faith_plot | kaki_alpha_evenness_plot) + plot_layout(guides = 'collect')
ggview(kaki_alpha_all_plot, units = "px", height = 2000, width = 5000)

all_alpha_hatch_plots <- (hihi_alpha_shannon_plot | hihi_alpha_observed_plot | hihi_alpha_faith_plot | hihi_alpha_evenness_plot) / 
  (kaki_alpha_shannon_plot | kaki_alpha_observed_plot | kaki_alpha_faith_plot | kaki_alpha_evenness_plot) + plot_layout(guides = 'collect')
ggview(all_alpha_hatch_plots, units = "px", height = 4000, width = 5000)

#all_alpha_beta_plots <- ((((PCoA_hihi_hatch_weighted | PCoA_hihi_hatch_unweighted | PCoA_hihi_hatch_jaccard | PCoA_hihi_hatch_bray) / 
#  (PCoA_kaki_hatch_weighted | PCoA_kaki_hatch_unweighted | PCoA_kaki_hatch_jaccard | PCoA_kaki_hatch_bray))  + plot_layout(guides = 'collect')) /
#  (((hihi_alpha_shannon_plot | hihi_alpha_observed_plot | hihi_alpha_faith_plot | hihi_alpha_evenness_plot) / 
#  (kaki_alpha_shannon_plot | kaki_alpha_observed_plot | kaki_alpha_faith_plot | kaki_alpha_evenness_plot)) + plot_layout(guides = 'collect')))

all_alpha_beta_plots <- ((hihi_alpha_shannon_plot | hihi_alpha_observed_plot | hihi_alpha_faith_plot | hihi_alpha_evenness_plot) / 
                           (kaki_alpha_shannon_plot | kaki_alpha_observed_plot | kaki_alpha_faith_plot | kaki_alpha_evenness_plot) /
                           (PCoA_hihi_hatch_weighted | PCoA_hihi_hatch_unweighted | PCoA_hihi_hatch_jaccard | PCoA_hihi_hatch_bray) / 
                           (PCoA_kaki_hatch_weighted | PCoA_kaki_hatch_unweighted | PCoA_kaki_hatch_jaccard | PCoA_kaki_hatch_bray)) + plot_layout(guides = 'collect') &
   theme(legend.position='bottom')

#all_alpha_beta_plots <- ((PCoA_hihi_hatch_weighted | PCoA_hihi_hatch_unweighted | PCoA_hihi_hatch_jaccard | PCoA_hihi_hatch_bray) / 
#                           (PCoA_kaki_hatch_weighted | PCoA_kaki_hatch_unweighted | PCoA_kaki_hatch_jaccard | PCoA_kaki_hatch_bray) + plot_layout(guides = 'collect')) /
#                           ((hihi_alpha_shannon_plot | hihi_alpha_observed_plot | hihi_alpha_faith_plot | hihi_alpha_evenness_plot) / 
#                           (kaki_alpha_shannon_plot | kaki_alpha_observed_plot | kaki_alpha_faith_plot | kaki_alpha_evenness_plot) + plot_layout(guides = 'collect'))

#all_alpha_beta_plots <- (all_hatch_PCoA) / (all_alpha_plots)

ggview(all_alpha_beta_plots, units = "px", height = 8000, width = 5000)
all_alpha_beta_plots

kaki_wildcap_alpha_shannon_plot <- kaki_alpha_data %>%
  ggplot(aes(x = wild.captive, y= shannon_entropy, fill = wild.captive, color = wild.captive, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab("Alpha Diversity Measure") +
  xlab(NULL) +
  #  scale_fill_discrete(name = "Egg Captive- or Wild-laid") +
  scale_fill_manual(values=c("limegreen","purple"),labels=c("Captive","Wild"),name = "Egg Captive- or Wild-laid") +
  scale_color_manual(values=c("limegreen","purple")) +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Shannon's Diversity Index") +
  stat_signif(comparisons = list(c("Captive","Wild")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) + #should be kruskal.test but won't work and wilcox.test seems to give same results
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
kaki_wildcap_alpha_shannon_plot
kaki_wildcap_alpha_observed_plot <- kaki_alpha_data %>%
  ggplot(aes(x = wild.captive, y= observed_features, fill = wild.captive, color = wild.captive, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab(NULL) +
  xlab(NULL) +
  #  scale_fill_discrete(name = "Egg Captive- or Wild-laid") +
  scale_fill_manual(values=c("limegreen","purple"),labels=c("Captive","Wild"),name = "Egg Captive- or Wild-laid") +
  scale_color_manual(values=c("limegreen","purple")) +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Observed Features") +
  stat_signif(comparisons = list(c("Captive","Wild")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
kaki_wildcap_alpha_observed_plot
kaki_wildcap_alpha_faith_plot <- kaki_alpha_data %>%
  ggplot(aes(x = wild.captive, y= faith_pd, fill = wild.captive, color = wild.captive, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab(NULL) +
  xlab(NULL) +
  #  scale_fill_discrete(name = "Egg Captive- or Wild-laid") +
  scale_fill_manual(values=c("limegreen","purple"),labels=c("Captive","Wild"),name = "Egg Captive- or Wild-laid") +
  scale_color_manual(values=c("limegreen","purple")) +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Faith's Phylogenetic Diversity") +
  stat_signif(comparisons = list(c("Captive","Wild")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
kaki_wildcap_alpha_faith_plot
kaki_wildcap_alpha_evenness_plot <- kaki_alpha_data %>%
  ggplot(aes(x = wild.captive, y= pielou_evenness, fill = wild.captive, color = wild.captive, alpha = 0.3)) +
  geom_boxplot(outlier.alpha = 1) +
  geom_jitter(shape = 16, color = "black", alpha = 1) + 
  ylab(NULL) +
  xlab(NULL) +
#  scale_fill_discrete(name = "Egg Captive- or Wild-laid") +
  scale_fill_manual(values=c("limegreen","purple"),labels=c("Captive","Wild"),name = "Egg Captive- or Wild-laid") +
  scale_color_manual(values=c("limegreen","purple")) +
  theme_bw() +
  guides(color = "none", alpha = "none") +
  ggtitle("Pielou's Evenness")  +
  stat_signif(comparisons = list(c("Captive","Wild")), test = "wilcox.test", map_signif_level = TRUE, color = "black", alpha = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
kaki_wildcap_alpha_evenness_plot

all_kaki_wildcap_alpha <- (kaki_wildcap_alpha_shannon_plot | kaki_wildcap_alpha_observed_plot | kaki_wildcap_alpha_faith_plot | kaki_wildcap_alpha_evenness_plot) + plot_layout(guides = 'collect')
ggview(all_kaki_wildcap_alpha, units = "px", height = 2000, width = 5000)

all_kaki_wildcap_alpha.beta_plots <- (kaki_wildcap_alpha_shannon_plot | kaki_wildcap_alpha_observed_plot | kaki_wildcap_alpha_faith_plot | kaki_wildcap_alpha_evenness_plot) / 
  (PCoA_kaki_wildcap_weighted | PCoA_kaki_wildcap_unweighted | PCoA_kaki_wildcap_jaccard | PCoA_kaki_wildcap_bray) + plot_layout(guides = 'collect')
ggview(all_kaki_wildcap_alpha.beta_plots, units = "px", height = 4000, width = 5000)
all_kaki_wildcap_alpha.beta_plots

##### Relative Abundance Plots #####
# Load data
rel.abund_phylum <- read.csv("Phyla_RelAbundance_df.csv")
rel.abund_family <- read.csv("Family_RelAbundance_df.csv")
rel.abund_genus <- read.csv("Genus_RelAbundance_df.csv")
rel.abund_phylum_top5 <- read.csv("Phyla_RelAbundance_df_top5.csv")
rel.abund_family_top5 <- read.csv("Family_RelAbundance_df_top5.csv")
rel.abund_genus_top5 <- read.csv("Genus_RelAbundance_df_top5.csv")

# Top 12 taxa overall
phyla_plot <- rel.abund_phylum %>%
  mutate(Phylum = fct_relevel(Phylum,"Other","Gemmatimonadota","Deinococcota","Chloroflexi","Bdellovibrionota","Verrucomicrobiota","Planctomycetota","Acidobacteriota","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
  ggplot(aes(x = Group, y = Rel.abund.percent, fill = Phylum)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance (%)") +
  ggtitle("Relative abundance of phyla") +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
phyla_plot

family_plot <- rel.abund_family %>%
  mutate(Family = fct_relevel(Family,"Other","Microbacteriaceae","Sanguibacteraceae","Micrococcaceae","Sphingomonadaceae","Xanthomonadaceae","Moraxellaceae","Staphylococcaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
  ggplot(aes(x = Group, y = Rel.abund.percent, fill = Family)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance (%)") +
  ggtitle("Relative abundance of families") +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
family_plot

my.labels_genera <- c("Other","Comamonadaceae genus",expression(italic("Sphingomonas")),expression(italic("Streptococcus")),expression(italic("Sanguibacter-Flavimobilis")),expression(italic("Acinetobacter")),expression(italic("Stenotrophomonas")),expression(italic("Staphylococcus")),expression(italic("Janthinobacterium")),expression(italic("Pseudoalteromonas")),expression(italic("Halomonas")),expression(italic("Ralstonia")),expression(italic("Pseudomonas")))
genera_plot <- rel.abund_genus %>%
  mutate(Genus = fct_relevel(Genus,"Other","Comamonadaceae genus","Sphingomonas","Streptococcus","Sanguibacter-Flavimobilis","Acinetobacter","Stenotrophomonas","Staphylococcus","Janthinobacterium","Pseudoalteromonas","Halomonas","Ralstonia","Pseudomonas")) %>%
  ggplot(aes(x = Group, y = Rel.abund.percent, fill = Genus)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance (%)") +
  ggtitle("Relative abundance of genera") +
  theme_classic() +
  scale_fill_discrete(name="Genus",labels=my.labels_genera)  +
    theme(legend.text.align = 0,
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))
genera_plot

# Top 5 in each
my.labels_phylum_top5 <- c("Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")
phyla_plot_top5 <- rel.abund_phylum_top5 %>%
  mutate(Phylum = fct_relevel(Phylum,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
  ggplot(aes(x = Group, y = Rel.abund.percent, fill = Phylum)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance (%)") +
  #ggtitle("Relative abundance of phyla") +
  theme_classic() +
  scale_fill_discrete(name="Phylum",labels=my.labels_phylum_top5) +
  scale_fill_manual(values=c("grey","#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"),labels=my.labels_phylum_top5) +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.justification = "left")
phyla_plot_top5

my.labels_family_top5 <- c("Other","Xanthomonadaceae","Staphylococcaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")
family_plot_top5 <- rel.abund_family_top5 %>%
  mutate(Family = fct_relevel(Family,"Other","Xanthomonadaceae","Staphylococcaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
  ggplot(aes(x = Group, y = Rel.abund.percent, fill = Family)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance (%)") +
  #ggtitle("Relative abundance of families") +
  theme_classic() +
  scale_fill_discrete(name="Family",labels=my.labels_family_top5) +
  scale_fill_manual(values=c("grey","#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FB61D7"),labels=my.labels_family_top5) +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.justification = "left")
family_plot_top5

my.labels_genera_top5 <- c("Other",expression(italic("Stenotrophomonas")),expression(italic("Staphylococcus")),expression(italic("Janthinobacterium")),expression(italic("Pseudoalteromonas")),expression(italic("Halomonas")),expression(italic("Ralstonia")),expression(italic("Pseudomonas")))
genera_plot_top5 <- rel.abund_genus_top5 %>%
  mutate(Genus = fct_relevel(Genus,"Other","Stenotrophomonas","Staphylococcus","Janthinobacterium","Pseudoalteromonas","Halomonas","Ralstonia","Pseudomonas")) %>%
  ggplot(aes(x = Group, y = Rel.abund.percent, fill = Genus)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance (%)") +
  #ggtitle("Relative abundance of genera") +
  theme_classic() +
  scale_fill_discrete(name="Genus",labels=my.labels_genera_top5) +
  scale_fill_manual(values=c("grey","#F8766D","#C49A00","#53B400","#00C094","#00B6EB","#A58AFF","#FB61D7"),labels=my.labels_genera_top5) +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.text.align = 0,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.justification = "left")
genera_plot_top5

allplot <- genera_plot_top5/
  family_plot_top5 /
  phyla_plot_top5
allplot

#ggview(allplot, units = "px", height = 3600, width = 2200)

##### Differential Abundance Plots #####
### DESeq2 Plots

# Hihi vs. kakī
deseq_ASV <- read.csv("DifferentialAbundance/DESeq2_Results/deseq_results_hihi.vs.kaki_qiime2.csv")
deseq_genus <- read.csv("DifferentialAbundance/DESeq2_Results/deseq_results_hihi.vs.kaki_qiime2_genus.csv")
deseq_family <- read.csv("DifferentialAbundance/DESeq2_Results/deseq_results_hihi.vs.kaki_qiime2_family.csv")
deseq_phylum <- read.csv("DifferentialAbundance/DESeq2_Results/deseq_results_hihi.vs.kaki_qiime2_phylum.csv")
dim(deseq_ASV) #35 14
dim(deseq_genus) #112 13
dim(deseq_family) #92 12
dim(deseq_phylum) #1 9
head(deseq_ASV)

## ASVs
ASV_phyla_labels <- c("Actinobacteriota","Firmicutes","Proteobacteria")
ASV_deseq2_plot <- ggplot(deseq_ASV, aes(x=log2FoldChange, y=X, fill=Phylum)) + 
  geom_bar(stat="identity") +
  xlim(-30,30) +
  ylab("ASV") +
  scale_fill_discrete(name="Family",labels=my.labels_family_top5) +
  scale_fill_manual(values=c("#E38900","#00BFC4","#FB61D7"),labels=ASV_phyla_labels)
ASV_deseq2_plot

deseq_ASV %>%
  count(Species)

deseq_ASV[1,]$Species <- "Pseudochrobactrum sp. [1]"
deseq_ASV[2,]$Species <- "Pseudochrobactrum sp. [2]"
deseq_ASV[3,]$Species <- "Sphingomonas sp."
deseq_ASV[4,]$Species <- "Acetobacteraceae bacterium"
deseq_ASV[5,]$Species <- "Devosia sp."
deseq_ASV[6,]$Species <- "Arthrobacter sp."
deseq_ASV[7,]$Species <- "Paeniglutamicibacter sp."
deseq_ASV[8,]$Species <- "Unidentified Micrococcaceae (Family) [1]"
deseq_ASV[9,]$Species <- "Unidentified Micrococcaceae (Family) [2]"
deseq_ASV[11,]$Species <- "Microbacterium sp. [1]"
deseq_ASV[12,]$Species <- "Microbacterium sp. [2]"
deseq_ASV[13,]$Species <- "Leucobacter sp."
deseq_ASV[14,]$Species <- "Pseudoclavibacter sp."
deseq_ASV[15,]$Species <- "Sanguibacter-Flavimobilis sp."
deseq_ASV[16,]$Species <- "Rhodococcus sp."
deseq_ASV[17,]$Species <- "Bacillus sp. [1]"
deseq_ASV[18,]$Species <- "Bacillus sp. [2]"
deseq_ASV[19,]$Species <- "Exiguobacterium sp."
deseq_ASV[20,]$Species <- "Stenotrophomonas rhizophila"
deseq_ASV[21,]$Species <- "Herbaspirillum sp."
deseq_ASV[22,]$Species <- "Verticiella sp."
deseq_ASV[23,]$Species <- "Janthinobacterium sp."
deseq_ASV[24,]$Species <- "Massilia sp."
deseq_ASV[26,]$Species <- "Pseudoalteromonas sp."
deseq_ASV[27,]$Species <- "Pseudomonas sp. [1]"
deseq_ASV[28,]$Species <- "Pseudomonas sp. [2]"
deseq_ASV[29,]$Species <- "Pseudomonas sp. [3]"
deseq_ASV[30,]$Species <- "Halomonas sp. [1]"
deseq_ASV[31,]$Species <- "Halomonas sp. [2]"
deseq_ASV[32,]$Species <- "Halomonas sp. [3]"
deseq_ASV[33,]$Species <- "Halomonas sp. [4]"
deseq_ASV[34,]$Species <- "Halomonas sp. [5]"
deseq_ASV[35,]$Species <- "Pseudomonas sp. [4]"

deseq_ASV %>%
  count(Species)

ASV_deseq2_plot <- ggplot(deseq_ASV, aes(x=log2FoldChange, y=factor(Species,levels=rev(levels(factor(Species)))), fill=Phylum)) + 
  geom_bar(stat="identity") +
  xlim(-30,30) +
  ylab("ASV") +
  scale_fill_manual(values=c("#E38900","#00BFC4","#FB61D7"),labels=ASV_phyla_labels) +
  scale_y_discrete(breaks=c("Acetobacteraceae bacterium","Arthrobacter sp.","Bacillus sp. [1]","Bacillus sp. [2]","Devosia sp.","Exiguobacterium sp.","Halomonas sp. [1]","Halomonas sp. [2]","Halomonas sp. [3]","Halomonas sp. [4]","Halomonas sp. [5]","Herbaspirillum sp.","Janthinobacterium sp.","Leucobacter sp.","Massilia sp.","Microbacterium sp. [1]","Microbacterium sp. [2]","Paeniglutamicibacter sp.","Pseudoalteromonas sp.","Pseudochrobactrum sp. [1]","Pseudochrobactrum sp. [2]","Pseudoclavibacter sp.","Pseudomonas sp. [1]","Pseudomonas sp. [2]","Pseudomonas sp. [3]","Pseudomonas sp. [4]","Rhodococcus sp.","Sanguibacter-Flavimobilis sp.","Sphingomonas sp.","Stenotrophomonas rhizophila","Unidentified Comamonadaceae (Family)","Unidentified Microbacteriaceae (Family)","Unidentified Micrococcaceae (Family) [1]","Unidentified Micrococcaceae (Family) [2]","Verticiella sp."),
                   labels=c("Acetobacteraceae bacterium",expression(paste(italic("Arthrobacter"), " sp.")),expression(paste(italic("Bacillus"), " sp. [1]")),expression(paste(italic("Bacillus"), " sp. [2]")),expression(paste(italic("Devosia"), " sp.")),expression(paste(italic("Exiguobacterium"), " sp.")),expression(paste(italic("Halomonas"), " sp. [1]")),expression(paste(italic("Halomonas"), " sp. [2]")),expression(paste(italic("Halomonas"), " sp. [3]")),expression(paste(italic("Halomonas"), " sp. [4]")),expression(paste(italic("Halomonas"), " sp. [5]")),expression(paste(italic("Herbaspirillum"), " sp.")),expression(paste(italic("Janthinobacterium"), " sp.")),expression(paste(italic("Leucobacter"), " sp.")),expression(paste(italic("Massilia"), " sp.")),expression(paste(italic("Microbacterium"), " sp. [1]")),expression(paste(italic("Microbacterium"), " sp. [2]")),expression(paste(italic("Paeniglutamicibacter"), " sp.")),expression(paste(italic("Pseudoalteromonas"), " sp.")),expression(paste(italic("Pseudochrobactrum"), " sp. [1]")),expression(paste(italic("Pseudochrobactrum"), " sp. [2]")),expression(paste(italic("Pseudoclavibacter"), " sp.")),expression(paste(italic("Pseudomonas"), " sp. [1]")),expression(paste(italic("Pseudomonas"), " sp. [2]")),expression(paste(italic("Pseudomonas"), " sp. [3]")),expression(paste(italic("Pseudomonas"), " sp. [4]")),expression(paste(italic("Rhodococcus"), " sp.")),expression(paste(italic("Sanguibacter-Flavimobilis"), " sp.")),expression(paste(italic("Sphingomonas"), " sp.")),expression(italic("Stenotrophomonas rhizophila")),"Unidentified Comamonadaceae (Family)","Unidentified Microbacteriaceae (Family)","Unidentified Micrococcaceae (Family) [1]","Unidentified Micrococcaceae (Family) [2]",expression(paste(italic("Verticiella"), " sp."))))
ASV_deseq2_plot

## Genera
deseq_genus[deseq_genus$Genus=="uncultured",]

deseq_genus[1,]$Genus <- "Uncultured Pirellulaceae (Family)"
deseq_genus[2,]$Genus <- "Uncultured Gemmataceae (Family)"
deseq_genus[5,]$Genus <- "Uncultured Isosphaeraceae (Family)"
deseq_genus[6,]$Genus <- "Uncultured Vicinamibacterales (Order)"
deseq_genus[13,]$Genus <- "Uncultured Gaiellales (Order)"
deseq_genus[14,]$Genus <- "Uncultured Acidobacteriales (Order)"
deseq_genus[26,]$Genus <- "Uncultured Rhizobiales_Incertae_Sedis (Family)"
deseq_genus[46,]$Genus <- "Uncultured Chitinophagaceae (Family)"
deseq_genus[48,]$Genus <- "Uncultured Saprospiraceae (Family)"
deseq_genus[75,]$Genus <- "Uncultured Solirubrobacteraceae (Family)"
deseq_genus[77,]$Genus <- "Uncultured Microtrichales (Order)"
deseq_genus[91,]$Genus <- "Uncultured Gemmatimonadaceae (Family)"
deseq_genus[96,]$Genus <- "Uncultured Neisseriaceae (Family)"

deseq_genus %>%
  count(Genus)

genus_deseq2_plot <- ggplot(deseq_genus, aes(x=log2FoldChange, y=factor(Genus,levels=rev(levels(factor(Genus)))), fill=Phylum)) + 
  geom_bar(stat="identity") +
  xlim(-10,10) +
  ylab("Genus") +
  scale_y_discrete(breaks=c("Acidothermus","Acinetobacter","Actinopolyspora","alphaI_cluster","Anaerococcus","Aquisphaera","Arcicella","Azotobacter","Bacillus","Bdellovibrio","Blastococcus","Brevundimonas","Brochothrix","Bryobacter","Candidatus_Alysiosphaera","Candidatus_Udaeobacter","Chryseobacterium","Chthoniobacter","Clostridium_sensu_stricto_1","Conexibacter","Corynebacterium","Crinalium_SAG_22.89","Cutibacterium","Cyanobium_PCC-6307","Deinococcus","Devosia","Dyadobacter","Ellin6067","Exiguobacterium","Ferruginibacter","Fibrella","Finegoldia","Flavobacterium","Fonticella","Gaiella","Gemmatimonas","Geodermatophilus","Haliangium","Halomonas","Herbaspirillum","Hymenobacter","IMCC26256","Janthinobacterium","KD4-96","Klenkia","Leptotrichia","Leucobacter","Longivirga","Massilia","Methylocystis","Methylotenera","Microbacterium","Modestobacter","Mycobacterium","Neisseria","P3OB-42","Paeniglutamicibacter","PeM15","Peptoniphilus","Polymorphobacter","Prevotella_7","Pseudoalteromonas","Pseudochrobactrum","Pseudoclavibacter","Pseudomonas","Ralstonia","Rhizorhapis","Rhodobacter","Rhodoblastus","Rhodococcus","Rhodoferax","Rhodomicrobium","Roseiarcus","Roseomonas","Sanguibacter-Flavimobilis","SC-I-84","Sediminibacterium","Segetibacter","Sphingomonas","Sphingorhabdus","Sporosarcina","Staphylococcus","Stenotrophomonas","Streptococcus","Tychonema_CCAP_1459-11B","Uncultured Acidobacteriales (Order)","Uncultured Chitinophagaceae (Family)","Uncultured Gaiellales (Order)","Uncultured Gemmataceae (Family)","Uncultured Gemmatimonadaceae (Family)","Uncultured Isosphaeraceae (Family)","Uncultured Microtrichales (Order)","Uncultured Neisseriaceae (Family)","Uncultured Pirellulaceae (Family)","Uncultured Rhizobiales_Incertae_Sedis (Family)","Uncultured Saprospiraceae (Family)","Uncultured Solirubrobacteraceae (Family)","Uncultured Vicinamibacterales (Order)","Unidentified Bacillales (Order)","Unidentified Blastocatellaceae (Family)","Unidentified Comamonadaceae (Family)","Unidentified Gammaproteobacteria (Class)","Unidentified Intrasporangiaceae (Family)","Unidentified Isosphaeraceae (Family)","Unidentified Kineosporiaceae (Family)","Unidentified Microbacteriaceae (Family)","Unidentified Micrococcaceae (Family)","Unidentified Pasteurellaceae (Family)","Unidentified Planococcaceae (Family)","Unidentified Xanthobacteraceae (Family)","Ureaplasma","Verticiella"),
                   labels=c(expression(italic("Acidothermus")),expression(italic("Acinetobacter")),expression(italic("Actinopolyspora")),expression(italic("alphaI cluster")),expression(italic("Anaerococcus")),expression(italic("Aquisphaera")),expression(italic("Arcicella")),expression(italic("Azotobacter")),expression(italic("Bacillus")),expression(italic("Bdellovibrio")),expression(italic("Blastococcus")),expression(italic("Brevundimonas")),expression(italic("Brochothrix")),expression(italic("Bryobacter")),expression(italic("Candidatus Alysiosphaera")),expression(italic("Candidatus Udaeobacter")),expression(italic("Chryseobacterium")),expression(italic("Chthoniobacter")),expression(italic("Clostridium sensu stricto 1")),expression(italic("Conexibacter")),expression(italic("Corynebacterium")),expression(italic("Crinalium SAG 22.89")),expression(italic("Cutibacterium")),expression(italic("Cyanobium PCC-6307")),expression(italic("Deinococcus")),expression(italic("Devosia")),expression(italic("Dyadobacter")),expression(italic("Ellin6067")),expression(italic("Exiguobacterium")),expression(italic("Ferruginibacter")),expression(italic("Fibrella")),expression(italic("Finegoldia")),expression(italic("Flavobacterium")),expression(italic("Fonticella")),expression(italic("Gaiella")),expression(italic("Gemmatimonas")),expression(italic("Geodermatophilus")),expression(italic("Haliangium")),expression(italic("Halomonas")),expression(italic("Herbaspirillum")),expression(italic("Hymenobacter")),expression(italic("IMCC26256")),expression(italic("Janthinobacterium")),expression(italic("KD4-96")),expression(italic("Klenkia")),
expression(italic("Leptotrichia")),expression(italic("Leucobacter")),expression(italic("Longivirga")),expression(italic("Massilia")),expression(italic("Methylocystis")),expression(italic("Methylotenera")),expression(italic("Microbacterium")),expression(italic("Modestobacter")),expression(italic("Mycobacterium")),expression(italic("Neisseria")),expression(italic("P3OB-42")),expression(italic("Paeniglutamicibacter")),expression(italic("PeM15")),expression(italic("Peptoniphilus")),expression(italic("Polymorphobacter")),expression(italic("Prevotella 7")),expression(italic("Pseudoalteromonas")),expression(italic("Pseudochrobactrum")),expression(italic("Pseudoclavibacter")),expression(italic("Pseudomonas")),expression(italic("Ralstonia")),expression(italic("Rhizorhapis")),expression(italic("Rhodobacter")),expression(italic("Rhodoblastus")),expression(italic("Rhodococcus")),expression(italic("Rhodoferax")),expression(italic("Rhodomicrobium")),expression(italic("Roseiarcus")),expression(italic("Roseomonas")),expression(italic("Sanguibacter-Flavimobilis")),expression(italic("SC-I-84")),expression(italic("Sediminibacterium")),expression(italic("Segetibacter")),expression(italic("Sphingomonas")),expression(italic("Sphingorhabdus")),expression(italic("Sporosarcina")),expression(italic("Staphylococcus")),expression(italic("Stenotrophomonas")),expression(italic("Streptococcus")),expression(italic("Tychonema CCAP 1459-11B")),"Uncultured Acidobacteriales (Order)","Uncultured Chitinophagaceae (Family)","Uncultured Gaiellales (Order)","Uncultured Gemmataceae (Family)","Uncultured Gemmatimonadaceae (Family)","Uncultured Isosphaeraceae (Family)","Uncultured Microtrichales (Order)","Uncultured Neisseriaceae (Family)","Uncultured Pirellulaceae (Family)","Uncultured Rhizobiales_Incertae_Sedis (Family)","Uncultured Saprospiraceae (Family)","Uncultured Solirubrobacteraceae (Family)","Uncultured Vicinamibacterales (Order)","Unidentified Bacillales (Order)","Unidentified Blastocatellaceae (Family)","Unidentified Comamonadaceae (Family)","Unidentified Gammaproteobacteria (Class)","Unidentified Intrasporangiaceae (Family)","Unidentified Isosphaeraceae (Family)","Unidentified Kineosporiaceae (Family)","Unidentified Microbacteriaceae (Family)","Unidentified Micrococcaceae (Family)","Unidentified Pasteurellaceae (Family)","Unidentified Planococcaceae (Family)","Unidentified Xanthobacteraceae (Family)",expression(italic("Ureaplasma")),expression(italic("Verticiella"))))
genus_deseq2_plot + geom_text(aes(x = -5,y="Staphylococcus",label="More abundant in kakī")) + geom_text(aes(x = 5,y="Staphylococcus",label="More abundant in hihi"))

ggview(genus_deseq2_plot, units = "px", height = 3600, width = 600)

## Family
deseq_family %>%
  count(Family)
deseq_family[deseq_family$Family=="uncultured",]

deseq_family[3,]$Family <- "Uncultured Vicinamibacterales (Order)"
deseq_family[9,]$Family <- "Uncultured Gaiellales (Order)"
deseq_family[11,]$Family <- "Uncultured Acidobacteriales (Order)"
deseq_family[52,]$Family <- "Uncultured Frankiales (Order)"
deseq_family[59,]$Family <- "Uncultured Microtrichales (Order)"

family_deseq2_plot <- ggplot(deseq_family, aes(x=log2FoldChange, y=factor(Family,levels=rev(levels(factor(Family)))), fill=Phylum)) + 
  geom_bar(stat="identity") +
  xlim(-10,10) +
  ylab("Family")
family_deseq2_plot

## Phylum
phylum_deseq2_plot <- ggplot(deseq_phylum, aes(x=log2FoldChange, y=Phylum, fill=Phylum)) + 
  geom_bar(stat="identity") +
  xlim(-2,2) +
  ylab("Phylum") +
  scale_fill_manual(values=c("#00BFC4"),labels="Firmicutes")
phylum_deseq2_plot

all_deseq_plots <- (ASV_deseq2_plot/plot_spacer()/phylum_deseq2_plot)|family_deseq2_plot|genus_deseq2_plot
all_deseq_plots

layout <- "
AACCDD
AACCDD
##CCDD
BBCCDD
"
ASV_deseq2_plot + phylum_deseq2_plot + genus_deseq2_plot + family_deseq2_plot +
  plot_layout(design = layout)

ASV_deseq2_gen.fam_plot <- genus_deseq2_plot + family_deseq2_plot + plot_layout(guides = "collect")
ggview(ASV_deseq2_gen.fam_plot, units = "px", height = 3000, width = 3000)

# Hihi hatched vs. unhatched
deseq_hihi_hatch_ASV <- read.csv("DifferentialAbundance/DESeq2_Results/deseq_results_hihi_hatched.unhatched_qiime2.csv")
deseq_hihi_hatch_genus <- read.csv("DifferentialAbundance/DESeq2_Results/deseq_results_hihi_hatched.unhatched_qiime2_genus.csv")
deseq_hihi_hatch_family <- read.csv("DifferentialAbundance/DESeq2_Results/deseq_results_hihi_hatched.unhatched_qiime2_family.csv")
dim(deseq_hihi_hatch_ASV) #1 14
dim(deseq_hihi_hatch_genus) #4 13
dim(deseq_hihi_hatch_family) #2 12
head(deseq_hihi_hatch_ASV)

## ASVs

deseq_hihi_hatch_ASV[1,]$Species <- "Haemophilus sp."

deseq_hihi_hatch_ASV %>%
  count(Species)
max(deseq_hihi_hatch_ASV$log2FoldChange)

ASV_deseq2_hihi_hatch_plot <- ggplot(deseq_hihi_hatch_ASV, aes(x=log2FoldChange, y=factor(Species,levels=rev(levels(factor(Species)))), fill=Phylum)) + 
  geom_bar(stat="identity") +
  xlim(-25,25) +
  ylab("ASV") +
#  scale_fill_manual(values=c("#E38900","#00BFC4","#FB61D7"),labels=ASV_phyla_labels) +
  scale_y_discrete(breaks=c("Haemophilus sp."),
                   labels=c(expression(paste(italic("Haemophilus"), " sp."))))
ASV_deseq2_hihi_hatch_plot

## Genera
deseq_hihi_hatch_genus %>%
  count(Genus)

max(deseq_hihi_hatch_genus$log2FoldChange)

genus_deseq2_hihi_hatch_plot <- ggplot(deseq_hihi_hatch_genus, aes(x=log2FoldChange, y=factor(Genus,levels=rev(levels(factor(Genus)))), fill=Phylum)) + 
  geom_bar(stat="identity") +
  xlim(-25,25) +
  ylab("Genus") +
  scale_y_discrete(breaks=c("Actinomyces","Haemophilus","Rhodococcus","Sphingobium"),
                   labels=c(expression(italic("Actinomyces")),expression(italic("Haemophilus")),expression(italic("Rhodococcus")),expression(italic("Sphingobium"))))
genus_deseq2_hihi_hatch_plot

## Family
deseq_hihi_hatch_family %>%
  count(Family)

max(deseq_hihi_hatch_family$log2FoldChange)

family_deseq2_hihi_hatch_plot <- ggplot(deseq_hihi_hatch_family, aes(x=log2FoldChange, y=factor(Family,levels=rev(levels(factor(Family)))), fill=Phylum)) + 
  geom_bar(stat="identity") +
  xlim(-25,25) +
  ylab("Family")
family_deseq2_hihi_hatch_plot

deseq_hihi_hatch_ASV
deseq_hihi_hatch_genus$Species <- ""
deseq_hihi_hatch_family$Genus <- ""
deseq_hihi_hatch_family$Species <- ""

all_hihi_hatch <- rbind(deseq_hihi_hatch_ASV,deseq_hihi_hatch_genus,deseq_hihi_hatch_family)
all_hihi_hatch$taxa.name[1] <- all_hihi_hatch$Species[1]
all_hihi_hatch$taxa.name[2:5] <- all_hihi_hatch$Genus[2:5]
all_hihi_hatch$taxa.name[6:7] <- all_hihi_hatch$Family[6:7]
all_hihi_hatch$taxa.level[1] <- "ASV"
all_hihi_hatch$taxa.level[2:5] <- "Genus"
all_hihi_hatch$taxa.level[6:7] <- "Family"

ASV_phyla_labels_2 <- c("Actinobacteriota","Bacteroidota","Proteobacteria")
all_deseq2_hihi_hatch_plot <- ggplot(all_hihi_hatch, aes(x=log2FoldChange, y=factor(taxa.name,levels=rev(levels(factor(taxa.name)))), fill=Phylum)) + 
  geom_bar(stat="identity") +
  xlim(-25,25) +
  ylab(NULL) +
  scale_fill_manual(values=c("#882255","#117733","#6699CC"),labels=ASV_phyla_labels_2) +
  scale_y_discrete(breaks=c("Haemophilus sp.","Actinomyces","Haemophilus","Rhodococcus","Sphingobium","Chitinophagaceae","Nocardiaceae"),
                   labels=c(expression(paste(italic("Haemophilus"), " sp.")),expression(italic("Actinomyces")),expression(italic("Haemophilus")),expression(italic("Rhodococcus")),expression(italic("Sphingobium")),"Chitinophagaceae","Nocardiaceae")) +
  facet_grid(row = vars(factor(taxa.level,levels=c("ASV","Genus","Family"))), scale = "free_y", space = "free", switch="y") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        plot.margin = unit(c(1,0,0,0),"cm"))
all_deseq2_hihi_hatch_plot
cvdPlot(all_deseq2_hihi_hatch_plot)

### ANCOM Plots
ASV_ancom <- read.table(file = "DifferentialAbundance/ANCOM_volcano_ASV.tsv", sep = '\t', header = TRUE)
genus_ancom <- read.table(file = "DifferentialAbundance/ANCOM_volcano_genus.tsv", sep = '\t', header = TRUE)
family_ancom <- read.table(file = "DifferentialAbundance/ANCOM_volcano_family.tsv", sep = '\t', header = TRUE)
phylum_ancom <- read.table(file = "DifferentialAbundance/ANCOM_volcano_phylum.tsv", sep = '\t', header = TRUE)
head(ASV_ancom)

ASV_ancom$Species_parse <- gsub("Unidentified", "Unid.", ASV_ancom$Species_parse)
ASV_ancom$Species_parse <- gsub("Uncultured", "Uncult.", ASV_ancom$Species_parse)

max(ASV_ancom$W)

ASV_volcano <- ggplot(ASV_ancom,aes(x=clr,y=W,color=reject_null_hypothesis,shape=clr>0)) +
  geom_point(size = 4) +
#  ggtitle("Differentially abundant ASVs") +
  xlab("Centered log-ratio (CLR)") +
  ylab("W") +
  ylim(NA,85) +
  scale_color_manual(values = c("grey","#F8766D"),
                    name = "Differentially Abundant",
                    breaks = c("FALSE","TRUE"),
                    labels = c("False","True")) +
  scale_shape_manual(values = c(17, 16),
                     name = NULL,
                     labels = c("Hihi > Kakī","Kakī > Hihi")) +
  geom_text_repel(data=ASV_ancom %>% filter(reject_null_hypothesis=="TRUE"),
            aes(label=Species_parse),
            colour = "black",
            size = 4.5,
            max.overlaps = Inf,
            box.padding = unit(0.2, "lines"),
            point.padding = unit(0.2, "lines"),
            parse=TRUE,
            segment.colour = "grey") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        plot.title = element_text(size = 16, face = "bold"))
ASV_volcano

genus_ancom$genera_parse <- gsub("Unidentified", "Unid.", genus_ancom$genera_parse)
genus_ancom$genera_parse <- gsub("Uncultured", "Uncult.", genus_ancom$genera_parse)

max(genus_ancom$W)

genus_volcano <- ggplot(genus_ancom,aes(x=clr,y=W,color=reject_null_hypothesis,shape=clr>0)) +
  geom_point(size = 4) +
#  ggtitle("Differentially abundant genera") +
  xlab("Centered log-ratio (CLR)") +
  ylab("W") +
  ylim(NA,125) +
  scale_color_manual(values = c("grey","#F8766D"),
                     name = "Differentially Abundant",
                     breaks = c("FALSE","TRUE"),
                     labels = c("False","True")) +
  scale_shape_manual(values = c(17, 16),
                     name = NULL,
                     labels = c("Hihi > Kakī","Kakī > Hihi")) +
  geom_text_repel(data=genus_ancom %>% filter(reject_null_hypothesis=="TRUE"),
                  aes(label=genera_parse),
                  colour = "black",
                  size = 4.5,
                  max.overlaps = Inf,
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"),
                  parse=TRUE,
                  segment.colour = "grey") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        plot.title = element_text(size = 16, face = "bold"))
genus_volcano

family_ancom$family <- gsub("Unidentified", "Unid.", family_ancom$family)

max(family_ancom$W)

family_ancom
family_volcano <- ggplot(family_ancom,aes(x=clr,y=W,color=reject_null_hypothesis,shape=clr>0)) +
  geom_point(size = 4) +
#  ggtitle("Differentially abundant families") +
  xlab("Centered log-ratio (CLR)") +
  ylab("W") +
  ylim(NA,95) +
  scale_color_manual(values = c("grey","#F8766D"),
                     name = "Differentially Abundant",
                     breaks = c("FALSE","TRUE"),
                     labels = c("False","True")) +
  scale_shape_manual(values = c(17, 16),
                     name = NULL,
                     labels = c("Hihi > Kakī","Kakī > Hihi")) +
  geom_text_repel(data=family_ancom %>% filter(reject_null_hypothesis=="TRUE"),
                  aes(label=family),
                  colour = "black",
                  size = 4.5,
                  max.overlaps = Inf,
                  box.padding = unit(0.8, "lines"),
                  point.padding = unit(0.4, "lines"),
                  segment.colour = "grey") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        plot.title = element_text(size = 16, face = "bold"))
family_volcano

phylum_ancom
phylum_volcano <- ggplot(phylum_ancom,aes(x=clr,y=W,color=reject_null_hypothesis,shape=clr>0)) +
  geom_point(size = 4) +
#  ggtitle("Differentially abundant phyla") +
  xlab("Centered log-ratio (CLR)") +
  ylab("W") +
  scale_color_manual(values = c("grey","#F8766D"),
                     name = "Differentially Abundant",
                     breaks = c("FALSE","TRUE"),
                     labels = c("False","True")) +
  scale_shape_manual(values = c(17, 16),
                     name = NULL,
                     labels = c("Hihi > Kakī","Kakī > Hihi")) +
  geom_text_repel(data=phylum_ancom %>% filter(reject_null_hypothesis=="TRUE"),
                  aes(label=phyla),
                  colour = "black",
                  size = 4.5,
                  max.overlaps = Inf,
                  segment.colour = "grey") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        plot.title = element_text(size = 16, face = "bold"))
phylum_volcano

all_volcano <- (ASV_volcano | genus_volcano) /
  (family_volcano | phylum_volcano) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
all_volcano
ggview(all_volcano, units = "px", height = 2000, width = 2000)

##### Following riffomonas tutorials #####
# Import metadata
mtd <- read_tsv(("sample-metadata-new-noNA-removemock_forR_wEnvironmentalData_2.tsv"), comment = "#q2")
mtd_df <- as.data.frame(mtd)
head(mtd_df)
colnames(mtd_df)[1] <- "sample_id"
colnames(mtd_df)[9] <- "study_species"
mtd_df <- mtd_df %>%
  select(sample_id,study_species,wild.captive,hatched.unhatched)

mtd_df$hatched.unhatched
mtd_df_NArm <- mtd_df %>%
  drop_na(hatched.unhatched)
mtd_df_NArm$hatched.unhatched

# Import data for count table
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
otu_counts_data <- read_tsv("data_for_otu_counts.txt")
colnames(otu_counts_data)[1:5]

otu_counts <- otu_counts_data %>%
  pivot_longer(-sample_id, names_to="otu", values_to = "count")

#Import taxonomies
taxonomy_silva_138.1_250.256_filtered_decontam_noMCnoNTC <- read_qza("silva_138.1_taxonomy_250.256-filtered-decontam_noMCnoNTC.qza")

tax_tab_silva_138.1_250.256_filtered_decontam <- taxonomy_silva_138.1_250.256_filtered_decontam_noMCnoNTC$data %>% 
  as.data.frame() %>%
  mutate(Taxon = gsub("D_0", "k", Taxon), Taxon = gsub("D_1", "p", Taxon),
         Taxon = gsub("D_2", "c", Taxon), Taxon = gsub("D_3", "o", Taxon),
         Taxon = gsub("D_4", "f", Taxon), Taxon = gsub("D_5", "g", Taxon),
         Taxon = gsub("D_6", "s", Taxon)) %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  select(-Confidence) #Warning message: Expected 7 pieces. Missing pieces filled with `NA` in 4525 rows [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, ...] = this matches the number of NAs in the species column
dim(tax_tab_silva_138.1_250.256_filtered_decontam) #6528 7

taxonomy <- tax_tab_silva_138.1_250.256_filtered_decontam
taxonomy <- tibble::rownames_to_column(taxonomy,"otu")
taxonomy[,2] <- gsub(taxonomy[,2],pattern = "d__", replacement = "")
taxonomy[,3] <- gsub(taxonomy[,3],pattern = " p__", replacement = "")
taxonomy[,4] <- gsub(taxonomy[,4],pattern = " c__", replacement = "")
taxonomy[,5] <- gsub(taxonomy[,5],pattern = " o__", replacement = "")
taxonomy[,6] <- gsub(taxonomy[,6],pattern = " f__", replacement = "")
taxonomy[,7] <- gsub(taxonomy[,7],pattern = " g__", replacement = "")
taxonomy[,8] <- gsub(taxonomy[,8],pattern = " s__", replacement = "")

taxonomy_NArm <- taxonomy

for (i in 2:nrow(taxonomy_NArm)){
  if (is.na(taxonomy_NArm[i,3])){
    kingdom <- paste("Kingdom_", taxonomy_NArm[i,2], sep = "")
    taxonomy_NArm[i, 3:8] <- kingdom
  } else if (is.na(taxonomy_NArm[i,4])){
    phylum <- paste("Phylum_", taxonomy_NArm[i,3], sep = "")
    taxonomy_NArm[i, 4:8] <- phylum
  } else if (is.na(taxonomy_NArm[i,5])){
    class <- paste("Class_", taxonomy_NArm[i,4], sep = "")
    taxonomy_NArm[i, 5:8] <- class
  } else if (is.na(taxonomy_NArm[i,6])){
    order <- paste("Order_", taxonomy_NArm[i,5], sep = "")
    taxonomy_NArm[i, 6:8] <- order
  } else if (is.na(taxonomy_NArm[i,7])){
    family <- paste("Family_", taxonomy_NArm[i,6], sep = "")
    taxonomy_NArm[i, 7:8] <- family
  } else if (is.na(taxonomy_NArm[i,8])){
    taxonomy_NArm$Species[i] <- paste("Genus",taxonomy_NArm$Genus[i], sep = "_")
  }
}
taxonomy$Genus[1:20]
taxonomy_NArm$Genus[1:20]

otu_rel_abund <- inner_join(mtd_df, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy_NArm, by="otu") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "otu"),
               names_to="level",
               values_to="taxon")

otu_rel_abund_NArm <- inner_join(mtd_df_NArm, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy_NArm, by="otu") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "otu"),
               names_to="level",
               values_to="taxon")

colnames(otu_rel_abund)
colnames(otu_rel_abund_NArm)

#saveRDS(otu_rel_abund,file="otu_rel_abund.rds")
#saveRDS(otu_rel_abund_NArm,file="otu_rel_abund_NArm.rds")

## Grouped stacked barcharts for groups
# Hihi vs. Kakī
phylum_rel_abund <- otu_rel_abund %>%
  filter(level=="Phylum") %>%
  group_by(study_species, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(study_species, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

#phylum_rel_abund %>%
#  arrange(desc(mean_rel_abund))

phylum_pool <- phylum_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 0.2,
            mean = mean(mean_rel_abund),
            .groups="drop")

phylum_rel_abund$study_species<-gsub("Kaki","Kakī",phylum_rel_abund$study_species)

#riff_hihi.vs.kakī_rel.abund_phyla_plot <- inner_join(phylum_rel_abund, phylum_pool, by="taxon") %>%
#  mutate(taxon = if_else(pool, "Other", taxon)) %>%
#  group_by(study_species, taxon) %>%
#  summarize(mean_rel_abund = sum(mean_rel_abund),
#            mean = min(mean),
#            .groups="drop") %>%
#  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
#  ggplot(aes(x=study_species, y=mean_rel_abund, fill=taxon)) +
#  geom_col() +
#  scale_fill_manual(name="Phylum",
#                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
##                    values=c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D","grey")) +
#                    values=c(brewer.pal(5,"Dark2"),"grey")) +
#  scale_x_discrete(labels=c("Hihi","Kakī")) +
#  scale_y_continuous(expand=c(0, 0)) +
#  labs(x="Study Species",
#       y="Mean Relative Abundance (%)") +
#  theme_classic() +
#  theme(axis.text.x = element_markdown(),
#        legend.text = element_markdown(),
#        legend.key.size = unit(10, "pt"),
#        legend.justification = "left")
#riff_hihi.vs.kakī_rel.abund_phyla_plot

riff_hihi.vs.kakī_rel.abund_phyla_plot <- inner_join(phylum_rel_abund, phylum_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(study_species, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
  ggplot(aes(x=study_species, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Phylum",
                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
                    #                    values=c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D","grey")) +
                    values=c(brewer.pal(5,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~study_species, scale = "free_x", space = "free", switch="x") +
  force_panelsizes(cols = c(1, 92/157)) +
  labs(x="Study Species",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_hihi.vs.kakī_rel.abund_phyla_plot

cvdPlot(riff_hihi.vs.kakī_rel.abund_phyla_plot)

family_rel_abund <- otu_rel_abund %>%
  filter(level=="Family") %>%
  group_by(study_species, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(study_species, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

family_pool <- family_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 2,
            mean = mean(mean_rel_abund),
            .groups="drop")

family_rel_abund$study_species<-gsub("Kaki","Kakī",family_rel_abund$study_species)

#riff_hihi.vs.kakī_rel.abund_family_plot <- inner_join(family_rel_abund, family_pool, by="taxon") %>%
#  mutate(taxon = if_else(pool, "Other", taxon)) %>%
#  group_by(study_species, taxon) %>%
#  summarize(mean_rel_abund = sum(mean_rel_abund),
#            mean = min(mean),
#            .groups="drop") %>%
#  mutate(taxon = fct_relevel(taxon,"Other","Moraxellaceae","Xanthomonadaceae","Staphylococcaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
#  ggplot(aes(x=study_species, y=mean_rel_abund, fill=taxon)) +
#  geom_col() +
#  scale_fill_manual(name="Family",
#                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Pseudoalteromonadaceae","Oxalobacteraceae","Staphylococcaceae","Xanthomonadaceae","Moraxellaceae","Other"),
##                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
#                    values=c(brewer.pal(8,"Dark2"),"grey")) +
#  scale_x_discrete(labels=c("Hihi","Kakī")) +
#  scale_y_continuous(expand=c(0, 0)) +
#  labs(x="Study Species",
#       y="Mean Relative Abundance (%)") +
#  theme_classic() +
#  theme(axis.text.x = element_markdown(),
#        legend.text = element_markdown(),
#        legend.key.size = unit(10, "pt"),
#        legend.justification = "left")
#riff_hihi.vs.kakī_rel.abund_family_plot

riff_hihi.vs.kakī_rel.abund_family_plot <- inner_join(family_rel_abund, family_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(study_species, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Moraxellaceae","Xanthomonadaceae","Staphylococcaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
  ggplot(aes(x=study_species, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Family",
                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Pseudoalteromonadaceae","Oxalobacteraceae","Staphylococcaceae","Xanthomonadaceae","Moraxellaceae","Other"),
                    #                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~study_species, scale = "free_x", space = "free", switch="x") +
  force_panelsizes(cols = c(1, 92/157)) +
  labs(x="Study Species",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_hihi.vs.kakī_rel.abund_family_plot

genus_rel_abund <- otu_rel_abund %>%
  filter(level=="Genus") %>%
  group_by(study_species, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(study_species, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "^(\\S*)$", "*\\1*"))

genus_pool <- genus_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1.55,
            mean = mean(mean_rel_abund),
            .groups="drop")

genus_rel_abund$study_species<-gsub("Kaki","Kakī",genus_rel_abund$study_species)

#riff_hihi.vs.kakī_rel.abund_genus_plot <- inner_join(genus_rel_abund, genus_pool, by="taxon") %>%
#  mutate(taxon = if_else(pool, "Other", taxon)) %>%
#  group_by(study_species, taxon) %>%
#  summarize(mean_rel_abund = sum(mean_rel_abund),
#            mean = min(mean),
#            .groups="drop") %>%
#  mutate(taxon = fct_relevel(taxon,"Other","*Acinetobacter*","*Stenotrophomonas*","*Staphylococcus*","*Janthinobacterium*","*Pseudoalteromonas*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
#  ggplot(aes(x=study_species, y=mean_rel_abund, fill=taxon)) +
#  geom_col() +
#  scale_fill_manual(name="Genus",
#                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Pseudoalteromonas*","*Janthinobacterium*","*Staphylococcus*","*Stenotrophomonas*","*Acinetobacter*","Other"),
##                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
#                    values=c(brewer.pal(8,"Dark2"),"grey")) +
#  scale_x_discrete(labels=c("Hihi","Kakī")) +
#  scale_y_continuous(expand=c(0, 0)) +
#  labs(x="Study Species",
#       y="Mean Relative Abundance (%)") +
#  theme_classic() +
#  theme(axis.text.x = element_markdown(),
#        legend.text = element_markdown(),
#        legend.key.size = unit(10, "pt"),
#        legend.justification = "left")
#riff_hihi.vs.kakī_rel.abund_genus_plot

riff_hihi.vs.kakī_rel.abund_genus_plot <- inner_join(genus_rel_abund, genus_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(study_species, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","*Acinetobacter*","*Stenotrophomonas*","*Staphylococcus*","*Janthinobacterium*","*Pseudoalteromonas*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
  ggplot(aes(x=study_species, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Genus",
                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Pseudoalteromonas*","*Janthinobacterium*","*Staphylococcus*","*Stenotrophomonas*","*Acinetobacter*","Other"),
                    #                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~study_species, scale = "free_x", space = "free", switch="x") +
  force_panelsizes(cols = c(1, 92/157)) +
  labs(x="Study Species",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_hihi.vs.kakī_rel.abund_genus_plot

# Hihi - hatched vs. unhatched
phylum_hihi_rel_abund <- otu_rel_abund_NArm %>%
  filter(level=="Phylum") %>%
  filter(study_species=="Hihi") %>%
  group_by(hatched.unhatched, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

phylum_hihi_rel_abund_df <- as.data.frame(phylum_hihi_rel_abund)

phylum_hihi_rel_abund_df[phylum_hihi_rel_abund_df$hatched.unhatched=="Hatched",] %>%
  arrange(desc(mean_rel_abund))

phylum_hihi_rel_abund[phylum_hihi_rel_abund$hatched.unhatched=="Unhatched",] %>%
  arrange(desc(mean_rel_abund))

phylum_hihi_rel_abund[phylum_hihi_rel_abund$taxon=="Cyanobacteria",]

phylum_hihi_pool <- phylum_hihi_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 0.25,
            mean = mean(mean_rel_abund),
            .groups="drop")

#riff_hihi_hatched.vs.unhatched_rel.abund_phyla_plot <- inner_join(phylum_hihi_rel_abund, phylum_hihi_pool, by="taxon") %>%
#  mutate(taxon = if_else(pool, "Other", taxon)) %>%
#  group_by(hatched.unhatched, taxon) %>%
#  summarize(mean_rel_abund = sum(mean_rel_abund),
#            mean = min(mean),
#            .groups="drop") %>%
#  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
#  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
#  geom_col() +
#  scale_fill_manual(name=NULL,
#                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
#                    values=c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D","grey")) +
#  scale_y_continuous(expand=c(0, 0)) +
#  labs(x="Hatching Success",
#       y="Mean Relative Abundance (%)") +
#  theme_classic() +
#  theme(axis.text.x = element_markdown(),
#        legend.text = element_markdown(),
#        legend.key.size = unit(10, "pt"))
#riff_hihi_hatched.vs.unhatched_rel.abund_phyla_plot

riff_hihi_hatched.vs.unhatched_rel.abund_phyla_plot <- inner_join(phylum_hihi_rel_abund, phylum_hihi_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Phylum",
                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
#                    values=c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D","grey")) +
                    values=c(brewer.pal(5,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x") +
  labs(x="Hatching Success",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_hihi_hatched.vs.unhatched_rel.abund_phyla_plot

family_hihi_rel_abund <- otu_rel_abund_NArm %>%
  filter(level=="Family") %>%
  filter(study_species=="Hihi") %>%
  group_by(hatched.unhatched, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

family_hihi_rel_abund_df <- as.data.frame(family_hihi_rel_abund)

family_hihi_rel_abund_df[family_hihi_rel_abund_df$hatched.unhatched=="Unhatched",] %>%
  arrange(desc(mean_rel_abund))

family_hihi_rel_abund[family_hihi_rel_abund$hatched.unhatched=="Unhatched",] %>%
  arrange(desc(mean_rel_abund))

family_hihi_rel_abund[family_hihi_rel_abund$taxon=="Micrococcaceae",]

family_hihi_pool <- family_hihi_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1.8,
            mean = mean(mean_rel_abund),
            .groups="drop")

#riff_hihi_hatched.vs.unhatched_rel.abund_family_plot <- inner_join(family_hihi_rel_abund, family_hihi_pool, by="taxon") %>%
#  mutate(taxon = if_else(pool, "Other", taxon)) %>%
#  group_by(hatched.unhatched, taxon) %>%
#  summarize(mean_rel_abund = sum(mean_rel_abund),
#            mean = min(mean),
#            .groups="drop") %>%
#  mutate(taxon = fct_relevel(taxon,"Other","Oxalobacteraceae","Moraxellaceae","Staphylococcaceae","Orbaceae","Pseudoalteromonadaceae","Pseudomonadaceae","Halomonadaceae","Burkholderiaceae")) %>%
#  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
#  geom_col() +
#  scale_fill_manual(name=NULL,
#                    breaks=c("Burkholderiaceae","Halomonadaceae","Pseudomonadaceae","Pseudoalteromonadaceae","Orbaceae","Staphylococcaceae","Moraxellaceae","Oxalobacteraceae","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
#  scale_y_continuous(expand=c(0, 0)) +
#  labs(x="Hatching Success",
#       y="Mean Relative Abundance (%)") +
#  theme_classic() +
#  theme(axis.text.x = element_markdown(),
#        legend.text = element_markdown(),
#        legend.key.size = unit(10, "pt"))
#riff_hihi_hatched.vs.unhatched_rel.abund_family_plot

riff_hihi_hatched.vs.unhatched_rel.abund_family_plot <- inner_join(family_hihi_rel_abund, family_hihi_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Oxalobacteraceae","Moraxellaceae","Staphylococcaceae","Orbaceae","Pseudoalteromonadaceae","Pseudomonadaceae","Halomonadaceae","Burkholderiaceae")) %>%
  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Family",
                    breaks=c("Burkholderiaceae","Halomonadaceae","Pseudomonadaceae","Pseudoalteromonadaceae","Orbaceae","Staphylococcaceae","Moraxellaceae","Oxalobacteraceae","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x") +
  labs(x="Hatching Success",
       y="Mean Relative Abundance (%)") +
  theme_classic()  +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_hihi_hatched.vs.unhatched_rel.abund_family_plot

genus_hihi_rel_abund <- otu_rel_abund_NArm %>%
  filter(level=="Genus") %>%
  filter(study_species=="Hihi") %>%
  group_by(hatched.unhatched, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "^(\\S*)$", "*\\1*"))

genus_hihi_rel_abund_df <- as.data.frame(genus_hihi_rel_abund)

genus_hihi_rel_abund_df[genus_hihi_rel_abund_df$hatched.unhatched=="Unhatched",] %>%
  arrange(desc(mean_rel_abund))

genus_hihi_rel_abund[genus_hihi_rel_abund$hatched.unhatched=="Unhatched",] %>%
  arrange(desc(mean_rel_abund))

genus_hihi_rel_abund[genus_hihi_rel_abund$taxon=="*Pseudochrobactrum*",]

genus_hihi_pool <- genus_hihi_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1.6,
            mean = mean(mean_rel_abund),
            .groups="drop")

#riff_hihi_hatched.vs.unhatched_rel.abund_genus_plot <- inner_join(genus_hihi_rel_abund, genus_hihi_pool, by="taxon") %>%
#  mutate(taxon = if_else(pool, "Other", taxon)) %>%
#  group_by(hatched.unhatched, taxon) %>%
#  summarize(mean_rel_abund = sum(mean_rel_abund),
#            mean = min(mean),
#            .groups="drop") %>%
#  mutate(taxon = fct_relevel(taxon,"Other","*Acinetobacter*","*Janthinobacterium*","*Staphylococcus*","*Orbus*","*Pseudoalteromonas*","*Pseudomonas*","*Halomonas*","*Ralstonia*")) %>%
#  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
#  geom_col() +
#  scale_fill_manual(name=NULL,
#                    breaks=c("*Ralstonia*","*Halomonas*","*Pseudomonas*","*Pseudoalteromonas*","*Orbus*","*Staphylococcus*","*Janthinobacterium*","*Acinetobacter*","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
#  scale_y_continuous(expand=c(0, 0)) +
#  labs(x="Hatching Success",
#       y="Mean Relative Abundance (%)") +
#  theme_classic() +
#  theme(axis.text.x = element_markdown(),
#        legend.text = element_markdown(),
#        legend.key.size = unit(10, "pt"))
#riff_hihi_hatched.vs.unhatched_rel.abund_genus_plot

riff_hihi_hatched.vs.unhatched_rel.abund_genus_plot <- inner_join(genus_hihi_rel_abund, genus_hihi_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","*Acinetobacter*","*Janthinobacterium*","*Staphylococcus*","*Orbus*","*Pseudoalteromonas*","*Pseudomonas*","*Halomonas*","*Ralstonia*")) %>%
  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Genus",
                    breaks=c("*Ralstonia*","*Halomonas*","*Pseudomonas*","*Pseudoalteromonas*","*Orbus*","*Staphylococcus*","*Janthinobacterium*","*Acinetobacter*","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x") +
  labs(x="Hatching Success",
       y="Mean Relative Abundance (%)") +
  theme_classic()  +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_hihi_hatched.vs.unhatched_rel.abund_genus_plot

all_relabund_hihi <- riff_hihi_hatched.vs.unhatched_rel.abund_phyla_plot / riff_hihi_hatched.vs.unhatched_rel.abund_family_plot / riff_hihi_hatched.vs.unhatched_rel.abund_genus_plot
ggview(all_relabund_hihi, units = "px", height = 3600, width = 2200)

# Kakī - hatched vs. unhatched
phylum_kaki_rel_abund <- otu_rel_abund %>%
  filter(level=="Phylum") %>%
  filter(study_species=="Kaki") %>%
  group_by(hatched.unhatched, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

phylum_kaki_rel_abund_df <- as.data.frame(phylum_kaki_rel_abund)

phylum_kaki_rel_abund_df[phylum_kaki_rel_abund_df$hatched.unhatched=="Hatched",] %>%
  arrange(desc(mean_rel_abund))

phylum_kaki_rel_abund[phylum_kaki_rel_abund$hatched.unhatched=="Unhatched",] %>%
  arrange(desc(mean_rel_abund))

phylum_kaki_pool <- phylum_kaki_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 0.2,
            mean = mean(mean_rel_abund),
            .groups="drop")

#riff_kakī_hatched.vs.unhatched_rel.abund_phyla_plot <- inner_join(phylum_kaki_rel_abund, phylum_kaki_pool, by="taxon") %>%
#  mutate(taxon = if_else(pool, "Other", taxon)) %>%
#  group_by(hatched.unhatched, taxon) %>%
#  summarize(mean_rel_abund = sum(mean_rel_abund),
#            mean = min(mean),
#            .groups="drop") %>%
#  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
#  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
#  geom_col() +
#  scale_fill_manual(name=NULL,
#                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
#                    values=c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D","grey")) +
#  scale_y_continuous(expand=c(0, 0)) +
#  labs(x="Hatching Success",
#       y="Mean Relative Abundance (%)") +
#  theme_classic() +
#  theme(axis.text.x = element_markdown(),
#        legend.text = element_markdown(),
#        legend.key.size = unit(10, "pt"))
#riff_kakī_hatched.vs.unhatched_rel.abund_phyla_plot

riff_kakī_hatched.vs.unhatched_rel.abund_phyla_plot <- inner_join(phylum_kaki_rel_abund, phylum_kaki_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Phylum",
                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
#                    values=c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D","grey")) +
                    values=c(brewer.pal(5,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x") +
  labs(x="Hatching Success",
       y="Mean Relative Abundance (%)") +
  theme_classic()  +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_kakī_hatched.vs.unhatched_rel.abund_phyla_plot

family_kaki_rel_abund <- otu_rel_abund %>%
  filter(level=="Family") %>%
  filter(study_species=="Kaki") %>%
  group_by(hatched.unhatched, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

family_kaki_rel_abund_df <- as.data.frame(family_kaki_rel_abund)

family_kaki_rel_abund_df[family_kaki_rel_abund_df$hatched.unhatched=="Unhatched",] %>%
  arrange(desc(mean_rel_abund))

family_kaki_rel_abund[family_kaki_rel_abund$hatched.unhatched=="Unhatched",] %>%
  arrange(desc(mean_rel_abund))

family_kaki_rel_abund[family_kaki_rel_abund$taxon=="Moraxellaceae",]

family_kaki_pool <- family_kaki_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1.5,
            mean = mean(mean_rel_abund),
            .groups="drop")

#riff_kakī_hatched.vs.unhatched_rel.abund_family_plot <- inner_join(family_kaki_rel_abund, family_kaki_pool, by="taxon") %>%
#  mutate(taxon = if_else(pool, "Other", taxon)) %>%
#  group_by(hatched.unhatched, taxon) %>%
#  summarize(mean_rel_abund = sum(mean_rel_abund),
#            mean = min(mean),
#            .groups="drop") %>%
#  mutate(taxon = fct_relevel(taxon,"Other","Micrococcaceae","Sanguibacteraceae","Pseudoalteromonadaceae","Oxalobacteraceae","Xanthomonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
#  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
#  geom_col() +
#  scale_fill_manual(name=NULL,
#                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Xanthomonadaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Sanguibacteraceae","Micrococcaceae","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
#  scale_y_continuous(expand=c(0, 0)) +
#  labs(x="Hatching Success",
#       y="Mean Relative Abundance (%)") +
#  theme_classic() +
#  theme(axis.text.x = element_markdown(),
#        legend.text = element_markdown(),
#        legend.key.size = unit(10, "pt"))
#riff_kakī_hatched.vs.unhatched_rel.abund_family_plot

riff_kakī_hatched.vs.unhatched_rel.abund_family_plot <- inner_join(family_kaki_rel_abund, family_kaki_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Micrococcaceae","Sanguibacteraceae","Pseudoalteromonadaceae","Oxalobacteraceae","Xanthomonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Family",
                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Xanthomonadaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Sanguibacteraceae","Micrococcaceae","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x") +
  labs(x="Hatching Success",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_kakī_hatched.vs.unhatched_rel.abund_family_plot

genus_kaki_rel_abund <- otu_rel_abund %>%
  filter(level=="Genus") %>%
  filter(study_species=="Kaki") %>%
  group_by(hatched.unhatched, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")%>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "^(\\S*)$", "*\\1*"))

genus_kaki_rel_abund_df <- as.data.frame(genus_kaki_rel_abund)

genus_kaki_rel_abund_df[genus_kaki_rel_abund_df$hatched.unhatched=="Hatched",] %>%
  arrange(desc(mean_rel_abund))

genus_kaki_rel_abund[genus_kaki_rel_abund$hatched.unhatched=="Unhatched",] %>%
  arrange(desc(mean_rel_abund))

genus_kaki_rel_abund[genus_kaki_rel_abund$taxon=="*Pseudochrobactrum*",]

genus_kaki_pool <- genus_kaki_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1.1,
            mean = mean(mean_rel_abund),
            .groups="drop")

#riff_kakī_hatched.vs.unhatched_rel.abund_genus_plot <- inner_join(genus_kaki_rel_abund, genus_kaki_pool, by="taxon") %>%
#  mutate(taxon = if_else(pool, "Other", taxon)) %>%
#  group_by(hatched.unhatched, taxon) %>%
#  summarize(mean_rel_abund = sum(mean_rel_abund),
#            mean = min(mean),
#            .groups="drop") %>%
#  mutate(taxon = fct_relevel(taxon,"Other","*Pseudochrobactrum*","*Sanguibacter-Flavimobilis*","*Pseudoalteromonas*","*Stenotrophomonas*","*Janthinobacterium*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
#  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
#  geom_col() +
#  scale_fill_manual(name=NULL,
#                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Janthinobacterium*","*Stenotrophomonas*","*Pseudoalteromonas*","*Sanguibacter-Flavimobilis*","*Pseudochrobactrum*","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
#  scale_y_continuous(expand=c(0, 0)) +
#  labs(x="Hatching Success",
#       y="Mean Relative Abundance (%)") +
#  theme_classic() +
#  theme(axis.text.x = element_markdown(),
#        legend.text = element_markdown(),
#        legend.key.size = unit(10, "pt"))
#riff_kakī_hatched.vs.unhatched_rel.abund_genus_plot

riff_kakī_hatched.vs.unhatched_rel.abund_genus_plot <- inner_join(genus_kaki_rel_abund, genus_kaki_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","*Pseudochrobactrum*","*Sanguibacter-Flavimobilis*","*Pseudoalteromonas*","*Stenotrophomonas*","*Janthinobacterium*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Genus",
                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Janthinobacterium*","*Stenotrophomonas*","*Pseudoalteromonas*","*Sanguibacter-Flavimobilis*","*Pseudochrobactrum*","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x") +
  labs(x="Hatching Success",
       y="Mean Relative Abundance (%)") +
  theme_classic()  +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_kakī_hatched.vs.unhatched_rel.abund_genus_plot

all_relabund_kakī <- riff_kakī_hatched.vs.unhatched_rel.abund_phyla_plot / riff_kakī_hatched.vs.unhatched_rel.abund_family_plot / riff_kakī_hatched.vs.unhatched_rel.abund_genus_plot
ggview(all_relabund_kakī, units = "px", height = 3600, width = 2200)

# Kakī - wild vs. captive
phylum_kaki_wildcap_rel_abund <- otu_rel_abund %>%
  filter(level=="Phylum") %>%
  filter(study_species=="Kaki") %>%
  group_by(wild.captive, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(wild.captive, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

phylum_kaki_wildcap_pool <- phylum_kaki_wildcap_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 0.3,
            mean = mean(mean_rel_abund),
            .groups="drop")

phylum_kaki_wildcap_rel_abund_df <- as.data.frame(phylum_kaki_wildcap_rel_abund)

phylum_kaki_wildcap_rel_abund_df[phylum_kaki_wildcap_rel_abund_df$wild.captive=="Wild",] %>%
  arrange(desc(mean_rel_abund))

phylum_kaki_wildcap_rel_abund[phylum_kaki_wildcap_rel_abund$wild.captive=="Wild",] %>%
  arrange(desc(mean_rel_abund))

#riff_kakī_wild.vs.captive_rel.abund_phyla_plot <- inner_join(phylum_kaki_wildcap_rel_abund, phylum_kaki_wildcap_pool, by="taxon") %>%
#  mutate(taxon = if_else(pool, "Other", taxon)) %>%
#  group_by(wild.captive, taxon) %>%
#  summarize(mean_rel_abund = sum(mean_rel_abund),
#            mean = min(mean),
#            .groups="drop") %>%
#  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
#  ggplot(aes(x=wild.captive, y=mean_rel_abund, fill=taxon)) +
#  geom_col() +
#  scale_fill_manual(name=NULL,
#                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
##                   values=c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D","grey")) +
#                    values=c(brewer.pal(5,"Dark2"),"grey")) +
#  scale_y_continuous(expand=c(0, 0)) +
#  labs(x="Egg Captive- or Wild-laid",
#       y="Mean Relative Abundance (%)") +
#  theme_classic() +
#  theme(axis.text.x = element_markdown(),
#        legend.text = element_markdown(),
#        legend.key.size = unit(10, "pt"))
#riff_kakī_wild.vs.captive_rel.abund_phyla_plot

riff_kakī_wild.vs.captive_rel.abund_phyla_plot <- inner_join(phylum_kaki_wildcap_rel_abund, phylum_kaki_wildcap_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(wild.captive, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
  ggplot(aes(x=wild.captive, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name=NULL,
                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
                    #                   values=c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D","grey")) +
                    values=c(brewer.pal(5,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~wild.captive, scale = "free_x", space = "free", switch="x") +
  labs(x="Egg Captive- or Wild-laid",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_kakī_wild.vs.captive_rel.abund_phyla_plot

family_kaki_wildcap_rel_abund <- otu_rel_abund %>%
  filter(level=="Family") %>%
  filter(study_species=="Kaki") %>%
  group_by(wild.captive, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(wild.captive, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

family_kaki_wildcap_pool <- family_kaki_wildcap_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1.64,
            mean = mean(mean_rel_abund),
            .groups="drop")

family_kaki_wildcap_rel_abund_df <- as.data.frame(family_kaki_wildcap_rel_abund)

family_kaki_wildcap_rel_abund_df[family_kaki_wildcap_rel_abund_df$wild.captive=="Wild",] %>%
  arrange(desc(mean_rel_abund))

family_kaki_wildcap_rel_abund[family_kaki_wildcap_rel_abund$wild.captive=="Wild",] %>%
  arrange(desc(mean_rel_abund))

#riff_kakī_wild.vs.captive_rel.abund_family_plot <- inner_join(family_kaki_wildcap_rel_abund, family_kaki_wildcap_pool, by="taxon") %>%
#  mutate(taxon = if_else(pool, "Other", taxon)) %>%
#  group_by(wild.captive, taxon) %>%
#  summarize(mean_rel_abund = sum(mean_rel_abund),
#            mean = min(mean),
#            .groups="drop") %>%
#  mutate(taxon = fct_relevel(taxon,"Other","Sanguibacteraceae","Staphylococcaceae","Pseudoalteromonadaceae","Oxalobacteraceae","Xanthomonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
#  ggplot(aes(x=wild.captive, y=mean_rel_abund, fill=taxon)) +
#  geom_col() +
#  scale_fill_manual(name=NULL,
#                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Xanthomonadaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Staphylococcaceae","Sanguibacteraceae","Other"),
##                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
#                     values=c(brewer.pal(8,"Dark2"),"grey")) +
#  scale_y_continuous(expand=c(0, 0)) +
#  labs(x="Egg Captive- or Wild-laid",
#       y="Mean Relative Abundance (%)") +
#  theme_classic() +
#  theme(axis.text.x = element_markdown(),
#        legend.text = element_markdown(),
#        legend.key.size = unit(10, "pt"))
#riff_kakī_wild.vs.captive_rel.abund_family_plot

riff_kakī_wild.vs.captive_rel.abund_family_plot <- inner_join(family_kaki_wildcap_rel_abund, family_kaki_wildcap_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(wild.captive, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Sanguibacteraceae","Staphylococcaceae","Pseudoalteromonadaceae","Oxalobacteraceae","Xanthomonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
  ggplot(aes(x=wild.captive, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name=NULL,
                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Xanthomonadaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Staphylococcaceae","Sanguibacteraceae","Other"),
                    #                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~wild.captive, scale = "free_x", space = "free", switch="x") +
  labs(x="Egg Captive- or Wild-laid",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_kakī_wild.vs.captive_rel.abund_family_plot

genus_kaki_wildcap_rel_abund <- otu_rel_abund %>%
  filter(level=="Genus") %>%
  filter(study_species=="Kaki") %>%
  group_by(wild.captive, sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(wild.captive, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")%>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "^(\\S*)$", "*\\1*"))

genus_kaki_wildcap_pool <- genus_kaki_wildcap_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1.5,
            mean = mean(mean_rel_abund),
            .groups="drop")

genus_kaki_wildcap_rel_abund_df <- as.data.frame(genus_kaki_wildcap_rel_abund)

genus_kaki_wildcap_rel_abund_df[genus_kaki_wildcap_rel_abund_df$wild.captive=="Wild",] %>%
  arrange(desc(mean_rel_abund))

genus_kaki_wildcap_rel_abund[genus_kaki_wildcap_rel_abund$wild.captive=="Wild",] %>%
  arrange(desc(mean_rel_abund))

#riff_kakī_wild.vs.captive_rel.abund_genus_plot <- inner_join(genus_kaki_wildcap_rel_abund, genus_kaki_wildcap_pool, by="taxon") %>%
#  mutate(taxon = if_else(pool, "Other", taxon)) %>%
#  group_by(wild.captive, taxon) %>%
#  summarize(mean_rel_abund = sum(mean_rel_abund),
#            mean = min(mean),
#            .groups="drop") %>%
#  mutate(taxon = fct_relevel(taxon,"Other","*Sanguibacter-Flavimobilis*","*Staphylococcus*","*Pseudoalteromonas*","*Stenotrophomonas*","*Janthinobacterium*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
#  ggplot(aes(x=wild.captive, y=mean_rel_abund, fill=taxon)) +
#  geom_col() +
#  scale_fill_manual(name=NULL,
#                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Janthinobacterium*","*Stenotrophomonas*","*Pseudoalteromonas*","*Staphylococcus*","*Sanguibacter-Flavimobilis*","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
#                    values=c(brewer.pal(8,"Dark2"),"grey")) +
#  scale_y_continuous(expand=c(0, 0)) +
#  labs(x="Egg Captive- or Wild-laid",
#       y="Mean Relative Abundance (%)") +
#  theme_classic() +
#  theme(axis.text.x = element_markdown(),
#        legend.text = element_markdown(),
#        legend.key.size = unit(10, "pt"))
#riff_kakī_wild.vs.captive_rel.abund_genus_plot

riff_kakī_wild.vs.captive_rel.abund_genus_plot <- inner_join(genus_kaki_wildcap_rel_abund, genus_kaki_wildcap_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(wild.captive, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","*Sanguibacter-Flavimobilis*","*Staphylococcus*","*Pseudoalteromonas*","*Stenotrophomonas*","*Janthinobacterium*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
  ggplot(aes(x=wild.captive, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name=NULL,
                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Janthinobacterium*","*Stenotrophomonas*","*Pseudoalteromonas*","*Staphylococcus*","*Sanguibacter-Flavimobilis*","Other"),
                    #                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~wild.captive, scale = "free_x", space = "free", switch="x") +
  labs(x="Egg Captive- or Wild-laid",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_kakī_wild.vs.captive_rel.abund_genus_plot

## Grouped stacked barcharts for individuals
# Hihi vs. Kakī
phylum_rel_abund_ind <- otu_rel_abund %>%
  filter(level=="Phylum") %>%
  group_by(study_species, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop")

phylum_pool_ind <- phylum_rel_abund_ind %>%
  group_by(study_species, taxon) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 0.2,
            mean = mean(mean),
            .groups="drop")

phylum_sample_order <- phylum_rel_abund_ind %>%
  filter(taxon == "Proteobacteria") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample_id, order)

phylum_rel_abund_ind$study_species<-gsub("Kaki","Kakī",phylum_rel_abund_ind$study_species)

riff_hihi.vs.kakī_rel.abund.ind_phyla_plot <- inner_join(phylum_rel_abund_ind, phylum_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, study_species, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
  inner_join(., phylum_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Phylum",
                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
#                    values=c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D","grey")) +
                    values=c(brewer.pal(5,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~study_species, scale = "free_x", space = "free", switch="x") +
  labs(x="Study Species",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_hihi.vs.kakī_rel.abund.ind_phyla_plot

family_rel_abund_ind <- otu_rel_abund %>%
  filter(level=="Family") %>%
  group_by(study_species, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop")

#saveRDS(family_rel_abund_ind,"family_rel_abund_ind.rds")

family_pool_ind <- family_rel_abund_ind %>%
  group_by(study_species, taxon) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 2,
            mean = mean(mean),
            .groups="drop")

family_sample_order <- family_rel_abund_ind %>%
  filter(taxon == "Pseudomonadaceae") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample_id, order)

family_rel_abund_ind$study_species<-gsub("Kaki","Kakī",family_rel_abund_ind$study_species)

riff_hihi.vs.kakī_rel.abund.ind_family_plot <- inner_join(family_rel_abund_ind, family_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, study_species, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Moraxellaceae","Xanthomonadaceae","Staphylococcaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
  inner_join(., family_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Family",
                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Pseudoalteromonadaceae","Oxalobacteraceae","Staphylococcaceae","Xanthomonadaceae","Moraxellaceae","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~study_species, scale = "free_x", space = "free", switch="x") +
  labs(x="Study Species",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_hihi.vs.kakī_rel.abund.ind_family_plot

genus_rel_abund_ind <- otu_rel_abund %>%
  filter(level=="Genus") %>%
  group_by(study_species, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "^(\\S*)$", "*\\1*"))

#saveRDS(genus_rel_abund_ind,"genus_rel_abund_ind.rds")

genus_pool_ind <- genus_rel_abund_ind %>%
  group_by(study_species, taxon) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 1.55,
            mean = mean(mean),
            .groups="drop")

genus_sample_order <- genus_rel_abund_ind %>%
  filter(taxon == "*Pseudomonas*") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample_id, order)

genus_rel_abund_ind$study_species<-gsub("Kaki","Kakī",genus_rel_abund_ind$study_species)

riff_hihi.vs.kakī_rel.abund.ind_genus_plot <- inner_join(genus_rel_abund_ind, genus_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, study_species, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","*Acinetobacter*","*Stenotrophomonas*","*Staphylococcus*","*Janthinobacterium*","*Pseudoalteromonas*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
  inner_join(., genus_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Genus",
                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Pseudoalteromonas*","*Janthinobacterium*","*Staphylococcus*","*Stenotrophomonas*","*Acinetobacter*","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~study_species, scale = "free_x", space = "free", switch="x") +
  labs(x="Study Species",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_hihi.vs.kakī_rel.abund.ind_genus_plot

hihi.vs.kaki_relabund_all <- (riff_hihi.vs.kakī_rel.abund_phyla_plot | riff_hihi.vs.kakī_rel.abund.ind_phyla_plot) /
  (riff_hihi.vs.kakī_rel.abund_family_plot | riff_hihi.vs.kakī_rel.abund.ind_family_plot) /
  (riff_hihi.vs.kakī_rel.abund_genus_plot | riff_hihi.vs.kakī_rel.abund.ind_genus_plot)

ggview(hihi.vs.kaki_relabund_all, units = "px", height = 4000, width = 3000)
hihi.vs.kaki_relabund_all

# Hihi - hatched vs. unhatched
phylum_hihi_rel_abund_ind <- otu_rel_abund_NArm %>%
  filter(level=="Phylum") %>%
  filter(study_species=="Hihi") %>%
  group_by(hatched.unhatched, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop")

phylum_hihi_pool_ind <- phylum_hihi_rel_abund_ind %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 0.25,
            mean = mean(mean),
            .groups="drop")

phylum_hihi_sample_order <- phylum_hihi_rel_abund_ind %>%
  filter(taxon == "Proteobacteria") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample_id, order)

riff_hihi_hatched.vs.unhatched_rel.abund.ind_phyla_plot <- inner_join(phylum_hihi_rel_abund_ind, phylum_hihi_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, hatched.unhatched, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
  inner_join(., phylum_hihi_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Phylum",
                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
#                    values=c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D","grey")) +
                    values=c(brewer.pal(5,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x") +
  labs(x="Hatching Success",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_hihi_hatched.vs.unhatched_rel.abund.ind_phyla_plot

family_hihi_rel_abund_ind <- otu_rel_abund_NArm %>%
  filter(level=="Family") %>%
  filter(study_species=="Hihi") %>%
  group_by(hatched.unhatched, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop")

family_hihi_pool_ind <- family_hihi_rel_abund_ind %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 1.8,
            mean = mean(mean),
            .groups="drop")

family_hihi_sample_order <- family_hihi_rel_abund_ind %>%
  filter(taxon == "Burkholderiaceae") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample_id, order)

riff_hihi_hatched.vs.unhatched_rel.abund.ind_family_plot <- inner_join(family_hihi_rel_abund_ind, family_hihi_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, hatched.unhatched, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Oxalobacteraceae","Moraxellaceae","Staphylococcaceae","Orbaceae","Pseudoalteromonadaceae","Pseudomonadaceae","Halomonadaceae","Burkholderiaceae")) %>%
  inner_join(., family_hihi_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Family",
                    breaks=c("Burkholderiaceae","Halomonadaceae","Pseudomonadaceae","Pseudoalteromonadaceae","Orbaceae","Staphylococcaceae","Moraxellaceae","Oxalobacteraceae","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x") +
  labs(x="Hatching Success",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_hihi_hatched.vs.unhatched_rel.abund.ind_family_plot

genus_hihi_rel_abund_ind <- otu_rel_abund_NArm %>%
  filter(level=="Genus") %>%
  filter(study_species=="Hihi") %>%
  group_by(hatched.unhatched, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "^(\\S*)$", "*\\1*"))

genus_hihi_pool_ind <- genus_hihi_rel_abund_ind %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 1.6,
            mean = mean(mean),
            .groups="drop")

genus_hihi_sample_order <- genus_hihi_rel_abund_ind %>%
  filter(taxon == "*Ralstonia*") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample_id, order)

riff_hihi_hatched.vs.unhatched_rel.abund.ind_genus_plot <- inner_join(genus_hihi_rel_abund_ind, genus_hihi_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, hatched.unhatched, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","*Acinetobacter*","*Janthinobacterium*","*Staphylococcus*","*Orbus*","*Pseudoalteromonas*","*Pseudomonas*","*Halomonas*","*Ralstonia*")) %>%
  inner_join(., genus_hihi_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Genus",
                    breaks=c("*Ralstonia*","*Halomonas*","*Pseudomonas*","*Pseudoalteromonas*","*Orbus*","*Staphylococcus*","*Janthinobacterium*","*Acinetobacter*","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x") +
  labs(x="Hatching Success",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_hihi_hatched.vs.unhatched_rel.abund.ind_genus_plot

# Kakī - hatched vs. unhatched
phylum_kaki_rel_abund_ind <- otu_rel_abund %>%
  filter(level=="Phylum") %>%
  filter(study_species=="Kaki") %>%
  group_by(hatched.unhatched, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop")

phylum_kaki_pool_ind <- phylum_kaki_rel_abund_ind %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 0.2,
            mean = mean(mean),
            .groups="drop")

phylum_kaki_sample_order <- phylum_kaki_rel_abund_ind %>%
  filter(taxon == "Proteobacteria") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample_id, order)

riff_kakī_hatched.vs.unhatched_rel.abund.ind_phyla_plot <- inner_join(phylum_kaki_rel_abund_ind, phylum_kaki_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, hatched.unhatched, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
  inner_join(., phylum_kaki_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Phylum",
                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
#                    values=c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D","grey")) +
                    values=c(brewer.pal(5,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x") +
  labs(x="Hatching Success",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_kakī_hatched.vs.unhatched_rel.abund.ind_phyla_plot

family_kaki_rel_abund_ind <- otu_rel_abund %>%
  filter(level=="Family") %>%
  filter(study_species=="Kaki") %>%
  group_by(hatched.unhatched, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop")

family_kaki_pool_ind <- family_kaki_rel_abund_ind %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 1.5,
            mean = mean(mean),
            .groups="drop")

family_kaki_sample_order <- family_kaki_rel_abund_ind %>%
  filter(taxon == "Pseudomonadaceae") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample_id, order)

riff_kakī_hatched.vs.unhatched_rel.abund.ind_family_plot <- inner_join(family_kaki_rel_abund_ind, family_kaki_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, hatched.unhatched, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Micrococcaceae","Sanguibacteraceae","Pseudoalteromonadaceae","Oxalobacteraceae","Xanthomonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
  inner_join(., family_kaki_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Family",
                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Xanthomonadaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Sanguibacteraceae","Micrococcaceae","Other"),
 #                   values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x") +
  labs(x="Hatching Success",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_kakī_hatched.vs.unhatched_rel.abund.ind_family_plot

genus_kaki_rel_abund_ind <- otu_rel_abund %>%
  filter(level=="Genus") %>%
  filter(study_species=="Kaki") %>%
  group_by(hatched.unhatched, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "^(\\S*)$", "*\\1*"))

genus_kaki_pool_ind <- genus_kaki_rel_abund_ind %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 1.1,
            mean = mean(mean),
            .groups="drop")

genus_kaki_sample_order <- genus_kaki_rel_abund_ind %>%
  filter(taxon == "*Pseudomonas*") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample_id, order)

riff_kakī_hatched.vs.unhatched_rel.abund.ind_genus_plot <- inner_join(genus_kaki_rel_abund_ind, genus_kaki_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, hatched.unhatched, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","*Pseudochrobactrum*","*Sanguibacter-Flavimobilis*","*Pseudoalteromonas*","*Stenotrophomonas*","*Janthinobacterium*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
  inner_join(., genus_kaki_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Genus",
                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Janthinobacterium*","*Stenotrophomonas*","*Pseudoalteromonas*","*Sanguibacter-Flavimobilis*","*Pseudochrobactrum*","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x") +
  labs(x="Hatching Success",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_kakī_hatched.vs.unhatched_rel.abund.ind_genus_plot

all_relabund_hihi_ind <- riff_hihi_hatched.vs.unhatched_rel.abund.ind_phyla_plot / riff_hihi_hatched.vs.unhatched_rel.abund.ind_family_plot / riff_hihi_hatched.vs.unhatched_rel.abund.ind_genus_plot
ggview(all_relabund_hihi_ind, units = "px", height = 3600, width = 2200)
all_relabund_kakī_ind <- riff_kakī_hatched.vs.unhatched_rel.abund.ind_phyla_plot / riff_kakī_hatched.vs.unhatched_rel.abund.ind_family_plot / riff_kakī_hatched.vs.unhatched_rel.abund.ind_genus_plot
ggview(all_relabund_kakī_ind, units = "px", height = 3600, width = 2200)

hihi_relabund_all <- (riff_hihi_hatched.vs.unhatched_rel.abund_phyla_plot | riff_hihi_hatched.vs.unhatched_rel.abund.ind_phyla_plot) /
  (riff_hihi_hatched.vs.unhatched_rel.abund_family_plot | riff_hihi_hatched.vs.unhatched_rel.abund.ind_family_plot) /
  (riff_hihi_hatched.vs.unhatched_rel.abund_genus_plot | riff_hihi_hatched.vs.unhatched_rel.abund.ind_genus_plot)

kaki_relabund_all <- (riff_kakī_hatched.vs.unhatched_rel.abund_phyla_plot | riff_kakī_hatched.vs.unhatched_rel.abund.ind_phyla_plot) /
  (riff_kakī_hatched.vs.unhatched_rel.abund_family_plot | riff_kakī_hatched.vs.unhatched_rel.abund.ind_family_plot) /
  (riff_kakī_hatched.vs.unhatched_rel.abund_genus_plot | riff_kakī_hatched.vs.unhatched_rel.abund.ind_genus_plot)

ggview(hihi_relabund_all, units = "px", height = 4000, width = 3000)
ggview(kaki_relabund_all, units = "px", height = 4000, width = 3000)

# Kakī - wild vs. captive
phylum_kaki_wildcap_rel_abund_ind <- otu_rel_abund %>%
  filter(level=="Phylum") %>%
  filter(study_species=="Kaki") %>%
  group_by(wild.captive, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop")

phylum_kaki_wildcap_pool_ind <- phylum_kaki_wildcap_rel_abund_ind %>%
  group_by(wild.captive, taxon) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 0.3,
            mean = mean(mean),
            .groups="drop")

phylum_kaki_wildcap_sample_order <- phylum_kaki_wildcap_rel_abund_ind %>%
  filter(taxon == "Proteobacteria") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample_id, order)

riff_kakī_wild.vs.captive_rel.abund.ind_phyla_plot <- inner_join(phylum_kaki_wildcap_rel_abund_ind, phylum_kaki_wildcap_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, wild.captive, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
  inner_join(., phylum_kaki_wildcap_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name=NULL,
                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
 #                   values=c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D","grey")) +
                    values=c(brewer.pal(5,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~wild.captive, scale = "free_x", space = "free", switch="x") +
  labs(x="Egg Captive- or Wild-laid",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_kakī_wild.vs.captive_rel.abund.ind_phyla_plot

family_kaki_wildcap_rel_abund_ind <- otu_rel_abund %>%
  filter(level=="Family") %>%
  filter(study_species=="Kaki") %>%
  group_by(wild.captive, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop")

family_kaki_wildcap_pool_ind <- family_kaki_wildcap_rel_abund_ind %>%
  group_by(wild.captive, taxon) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 1.64,
            mean = mean(mean),
            .groups="drop")

family_kaki_wildcap_sample_order <- family_kaki_wildcap_rel_abund_ind %>%
  filter(taxon == "Pseudomonadaceae") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample_id, order)

riff_kakī_wild.vs.captive_rel.abund.ind_family_plot <- inner_join(family_kaki_wildcap_rel_abund_ind, family_kaki_wildcap_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, wild.captive, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Sanguibacteraceae","Staphylococcaceae","Pseudoalteromonadaceae","Oxalobacteraceae","Xanthomonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
  inner_join(., family_kaki_wildcap_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name=NULL,
                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Xanthomonadaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Staphylococcaceae","Sanguibacteraceae","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~wild.captive, scale = "free_x", space = "free", switch="x") +
  labs(x="Egg Captive- or Wild-laid",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_kakī_wild.vs.captive_rel.abund.ind_family_plot

genus_kaki_wildcap_rel_abund_ind <- otu_rel_abund %>%
  filter(level=="Genus") %>%
  filter(study_species=="Kaki") %>%
  group_by(wild.captive, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon,
                             "^(\\S*)$", "*\\1*"))

genus_kaki_wildcap_pool_ind <- genus_kaki_wildcap_rel_abund_ind %>%
  group_by(wild.captive, taxon) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 1.5,
            mean = mean(mean),
            .groups="drop")

genus_kaki_wildcap_sample_order <- genus_kaki_wildcap_rel_abund_ind %>%
  filter(taxon == "*Pseudomonas*") %>%
  arrange(desc(rel_abund)) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample_id, order)

riff_kakī_wild.vs.captive_rel.abund.ind_genus_plot <- inner_join(genus_kaki_wildcap_rel_abund_ind, genus_kaki_wildcap_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, wild.captive, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","*Sanguibacter-Flavimobilis*","*Staphylococcus*","*Pseudoalteromonas*","*Stenotrophomonas*","*Janthinobacterium*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
  inner_join(., genus_kaki_wildcap_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name=NULL,
                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Janthinobacterium*","*Stenotrophomonas*","*Pseudoalteromonas*","*Staphylococcus*","*Sanguibacter-Flavimobilis*","Other"),
#                    values=c("#FF61CC","#C77CFF","#00A9FF","#00BFC4","#00BE67","#7CAE00","#CD9600","#F8766D","grey")) +
                    values=c(brewer.pal(8,"Dark2"),"grey")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~wild.captive, scale = "free_x", space = "free", switch="x") +
  labs(x="Egg Captive- or Wild-laid",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_text(size = 12),
        axis.title=element_text(size=12))
riff_kakī_wild.vs.captive_rel.abund.ind_genus_plot

kakī_wild.vs.captive_relabund_all <- (riff_kakī_wild.vs.captive_rel.abund_phyla_plot | riff_kakī_wild.vs.captive_rel.abund.ind_phyla_plot) /
  (riff_kakī_wild.vs.captive_rel.abund_family_plot | riff_kakī_wild.vs.captive_rel.abund.ind_family_plot) /
  (riff_kakī_wild.vs.captive_rel.abund_genus_plot | riff_kakī_wild.vs.captive_rel.abund.ind_genus_plot)

ggview(kakī_wild.vs.captive_relabund_all, units = "px", height = 4000, width = 3000)
kakī_wild.vs.captive_relabund_all

# Hihi/Kakī hatched vs. unhatched family level only

plot_grafify_palette(palette="safe")
plot_grafify_palette(palette="kelly")
plot_grafify_palette(palette="muted")

"Burkholderiaceae" = safe_wine_#661100
"Halomonadaceae" = safe_skyblue_#6699CC
"Pseudomonadaceae" = safe_green_#117733
"Pseudoalteromonadaceae" = safe_reddish_#882255
"Orbaceae" = safe_bluegreen_#44AA99
"Staphylococcaceae" = safe_violet_#332288
"Moraxellaceae" = safe_bush_#999933
"Oxalobacteraceae" = safe_purple_#AA4499
"Xanthomonadaceae" = safe_yellow_#DDCC77
"Sanguibacteraceae" = safe_blue_#88CCEE
"Micrococcaceae" = safe_red_#CC6677
"Other" = grey

pretty_hihi_label <- c("Hatched" = "Hatched<br>(*n* = 84)",
                       "Unhatched" = "Unhatched<br>(*n* = 22)")
riff_hihi_hatched.vs.unhatched_rel.abund_family_plot_NEW <- inner_join(family_hihi_rel_abund, family_hihi_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Oxalobacteraceae","Moraxellaceae","Staphylococcaceae","Orbaceae","Pseudoalteromonadaceae","Pseudomonadaceae","Halomonadaceae","Burkholderiaceae")) %>%
  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Family",
                    breaks=c("Burkholderiaceae","Halomonadaceae","Pseudomonadaceae","Pseudoalteromonadaceae","Orbaceae","Staphylococcaceae","Moraxellaceae","Oxalobacteraceae","Other"),
                    values=c("#661100","#6699CC","#117733","#882255","#44AA99","#332288","#999933","#AA4499","#999999")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x",
             labeller = labeller(hatched.unhatched = pretty_hihi_label)) +
  force_panelsizes(cols = c(1, 22/84)) +
  labs(x="Hatching Success",
       y="Mean Relative Abundance (%)") +
  theme_classic()  +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_markdown(size = 12),
        axis.title=element_text(size=12))
riff_hihi_hatched.vs.unhatched_rel.abund_family_plot_NEW

riff_hihi_hatched.vs.unhatched_rel.abund.ind_family_plot_NEW <- inner_join(family_hihi_rel_abund_ind, family_hihi_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, hatched.unhatched, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Oxalobacteraceae","Moraxellaceae","Staphylococcaceae","Orbaceae","Pseudoalteromonadaceae","Pseudomonadaceae","Halomonadaceae","Burkholderiaceae")) %>%
  inner_join(., family_hihi_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Family",
                    breaks=c("Burkholderiaceae","Halomonadaceae","Pseudomonadaceae","Pseudoalteromonadaceae","Orbaceae","Staphylococcaceae","Moraxellaceae","Oxalobacteraceae","Other"),
                    values=c("#661100","#6699CC","#117733","#882255","#44AA99","#332288","#999933","#AA4499","#999999")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x",
             labeller = labeller(hatched.unhatched = pretty_hihi_label)) +
  labs(x="Hatching Success",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_markdown(size = 12),
        axis.title=element_text(size=12))
riff_hihi_hatched.vs.unhatched_rel.abund.ind_family_plot_NEW

pretty_kaki_label <- c("Hatched" = "Hatched<br>(*n* = 73)",
                       "Unhatched" = "Unhatched<br>(*n* = 19)")
riff_kakī_hatched.vs.unhatched_rel.abund_family_plot_NEW <- inner_join(family_kaki_rel_abund, family_kaki_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(hatched.unhatched, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Micrococcaceae","Sanguibacteraceae","Pseudoalteromonadaceae","Oxalobacteraceae","Xanthomonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
  ggplot(aes(x=hatched.unhatched, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Family",
                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Xanthomonadaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Sanguibacteraceae","Micrococcaceae","Other"),
                    values=c("#117733","#661100","#6699CC","#DDCC77","#AA4499","#882255","#88CCEE","#CC6677","#999999")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x",
             labeller = labeller(hatched.unhatched = pretty_kaki_label)) +
  force_panelsizes(cols = c(1, 19/73)) +
  labs(x="Hatching Success",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_markdown(size = 12),
        axis.title=element_text(size=12))
riff_kakī_hatched.vs.unhatched_rel.abund_family_plot_NEW

riff_kakī_hatched.vs.unhatched_rel.abund.ind_family_plot_NEW <- inner_join(family_kaki_rel_abund_ind, family_kaki_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, hatched.unhatched, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Micrococcaceae","Sanguibacteraceae","Pseudoalteromonadaceae","Oxalobacteraceae","Xanthomonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
  inner_join(., family_kaki_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Family",
                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Xanthomonadaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Sanguibacteraceae","Micrococcaceae","Other"),
                    values=c("#117733","#661100","#6699CC","#DDCC77","#AA4499","#882255","#88CCEE","#CC6677","#999999")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~hatched.unhatched, scale = "free_x", space = "free", switch="x",
             labeller = labeller(hatched.unhatched = pretty_kaki_label)) +
  labs(x="Hatching Success",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_markdown(size = 12),
        axis.title=element_text(size=12))
riff_kakī_hatched.vs.unhatched_rel.abund.ind_family_plot_NEW

hihi.kaki_hatched.vs.unhatched_family <- (riff_hihi_hatched.vs.unhatched_rel.abund_family_plot_NEW | riff_hihi_hatched.vs.unhatched_rel.abund.ind_family_plot_NEW) /
  (riff_kakī_hatched.vs.unhatched_rel.abund_family_plot_NEW | riff_kakī_hatched.vs.unhatched_rel.abund.ind_family_plot_NEW)
ggview(hihi.kaki_hatched.vs.unhatched_family, units = "px", height = 3000, width = 3000)
hihi.kaki_hatched.vs.unhatched_family

# Hihi vs. kaki - genus level only

pretty_hihi.vs.kaki_label <- c("Hihi" = "Hihi<br>(*n* = 157)",
                               "Kakī" = "Kakī<br>(*n* = 92)")
riff_hihi.vs.kakī_rel.abund_genus_plot_NEW<- inner_join(genus_rel_abund, genus_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(study_species, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","*Acinetobacter*","*Stenotrophomonas*","*Staphylococcus*","*Janthinobacterium*","*Pseudoalteromonas*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
  ggplot(aes(x=study_species, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Genus",
                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Pseudoalteromonas*","*Janthinobacterium*","*Staphylococcus*","*Stenotrophomonas*","*Acinetobacter*","Other"),
                    values=c("#882255","#999933","#44AA99","#AA4499","#332288","#117733","#DDCC77","#CC6677","#999999")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~study_species, scale = "free_x", space = "free", switch="x",
             labeller = labeller(study_species = pretty_hihi.vs.kaki_label)) +
  force_panelsizes(cols = c(1, 92/157)) +
  labs(x="Study Species",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_markdown(size = 12),
        axis.title=element_text(size=12))
riff_hihi.vs.kakī_rel.abund_genus_plot_NEW

riff_hihi.vs.kakī_rel.abund.ind_genus_plot_NEW <- inner_join(genus_rel_abund_ind, genus_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, study_species, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","*Acinetobacter*","*Stenotrophomonas*","*Staphylococcus*","*Janthinobacterium*","*Pseudoalteromonas*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
  inner_join(., genus_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Genus",
                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Pseudoalteromonas*","*Janthinobacterium*","*Staphylococcus*","*Stenotrophomonas*","*Acinetobacter*","Other"),
                    values=c("#882255","#999933","#44AA99","#AA4499","#332288","#117733","#DDCC77","#CC6677","#999999")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~study_species, scale = "free_x", space = "free", switch="x",
             labeller = labeller(study_species = pretty_hihi.vs.kaki_label)) +
  labs(x="Study Species",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_markdown(size = 12),
        axis.title=element_text(size=12))
riff_hihi.vs.kakī_rel.abund.ind_genus_plot_NEW

hihi.vs.kaki_genus <- (riff_hihi.vs.kakī_rel.abund_genus_plot_NEW | riff_hihi.vs.kakī_rel.abund.ind_genus_plot_NEW)
ggview(hihi.vs.kaki_genus, units = "px", height = 2000, width = 3000)
hihi.vs.kaki_genus

# Kakī - wild vs. captive

plot_grafify_palette(palette="safe")

pretty_kaki_wildcap_label <- c("Captive" = "Captive-laid<br>(*n* = 23)",
                               "Wild" = "Wild-laid<br>(*n* = 69)")

riff_kakī_wild.vs.captive_rel.abund_phyla_plot_NEW <- inner_join(phylum_kaki_wildcap_rel_abund, phylum_kaki_wildcap_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(wild.captive, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
  ggplot(aes(x=wild.captive, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_grafify(palette = "safe") +
  scale_fill_manual(name="Phylum",
                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
                    values=c("#882255","#999933","#44AA99","#AA4499","#332288","#999999")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~wild.captive, scale = "free_x", space = "free", switch="x",
             labeller = labeller(wild.captive = pretty_kaki_wildcap_label)) +
  force_panelsizes(cols = c(23/69, 1)) +
  labs(x="Egg Captive- or Wild-laid",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_markdown(size = 12),
        axis.title=element_text(size=12))
riff_kakī_wild.vs.captive_rel.abund_phyla_plot_NEW

#882255
#999933
#44AA99
#AA4499
#332288
#117733
#DDCC77
#CC6677

riff_kakī_wild.vs.captive_rel.abund_family_plot_NEW <- inner_join(family_kaki_wildcap_rel_abund, family_kaki_wildcap_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(wild.captive, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Sanguibacteraceae","Staphylococcaceae","Pseudoalteromonadaceae","Oxalobacteraceae","Xanthomonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
  ggplot(aes(x=wild.captive, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Family",
                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Xanthomonadaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Staphylococcaceae","Sanguibacteraceae","Other"),
                    values=c("#882255","#999933","#44AA99","#AA4499","#332288","#117733","#DDCC77","#CC6677","#999999")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~wild.captive, scale = "free_x", space = "free", switch="x",
             labeller = labeller(wild.captive = pretty_kaki_wildcap_label)) +
  force_panelsizes(cols = c(23/69, 1)) +
  labs(x="Egg Captive- or Wild-laid",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_markdown(size = 12),
        axis.title=element_text(size=12))
riff_kakī_wild.vs.captive_rel.abund_family_plot_NEW

riff_kakī_wild.vs.captive_rel.abund_genus_plot_NEW <- inner_join(genus_kaki_wildcap_rel_abund, genus_kaki_wildcap_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(wild.captive, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","*Sanguibacter-Flavimobilis*","*Staphylococcus*","*Pseudoalteromonas*","*Stenotrophomonas*","*Janthinobacterium*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
  ggplot(aes(x=wild.captive, y=mean_rel_abund, fill=taxon)) +
  geom_col(width = 1.5) +
  scale_fill_manual(name="Genus",
                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Janthinobacterium*","*Stenotrophomonas*","*Pseudoalteromonas*","*Staphylococcus*","*Sanguibacter-Flavimobilis*","Other"),
                    values=c("#882255","#999933","#44AA99","#AA4499","#332288","#117733","#DDCC77","#CC6677","#999999")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~wild.captive, scale = "free_x", space = "free", switch="x",
             labeller = labeller(wild.captive = pretty_kaki_wildcap_label)) +
  force_panelsizes(cols = c(23/69, 1)) +
  labs(x="Egg Captive- or Wild-laid",
       y="Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        legend.position = "none",
        strip.text.x = element_markdown(size = 12),
        axis.title=element_text(size=12))
riff_kakī_wild.vs.captive_rel.abund_genus_plot_NEW

riff_kakī_wild.vs.captive_rel.abund.ind_phyla_plot_NEW <- inner_join(phylum_kaki_wildcap_rel_abund_ind, phylum_kaki_wildcap_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, wild.captive, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Cyanobacteria","Bacteroidota","Firmicutes","Actinobacteriota","Proteobacteria")) %>%
  inner_join(., phylum_kaki_wildcap_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Phylum",
                    breaks=c("Proteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Cyanobacteria","Other"),
                    values=c("#882255","#999933","#44AA99","#AA4499","#332288","#999999")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~wild.captive, scale = "free_x", space = "free", switch="x",
             labeller = labeller(wild.captive = pretty_kaki_wildcap_label)) +
  labs(x="Egg Captive- or Wild-laid",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_markdown(size = 12),
        axis.title=element_text(size=12))
riff_kakī_wild.vs.captive_rel.abund.ind_phyla_plot_NEW

riff_kakī_wild.vs.captive_rel.abund.ind_family_plot_NEW <- inner_join(family_kaki_wildcap_rel_abund_ind, family_kaki_wildcap_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, wild.captive, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","Sanguibacteraceae","Staphylococcaceae","Pseudoalteromonadaceae","Oxalobacteraceae","Xanthomonadaceae","Halomonadaceae","Burkholderiaceae","Pseudomonadaceae")) %>%
  inner_join(., family_kaki_wildcap_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Family",
                    breaks=c("Pseudomonadaceae","Burkholderiaceae","Halomonadaceae","Xanthomonadaceae","Oxalobacteraceae","Pseudoalteromonadaceae","Staphylococcaceae","Sanguibacteraceae","Other"),
                    values=c("#882255","#999933","#44AA99","#AA4499","#332288","#117733","#DDCC77","#CC6677","#999999")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~wild.captive, scale = "free_x", space = "free", switch="x",
             labeller = labeller(wild.captive = pretty_kaki_wildcap_label)) +
  labs(x="Egg Captive- or Wild-laid",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_markdown(size = 12),
        axis.title=element_text(size=12))
riff_kakī_wild.vs.captive_rel.abund.ind_family_plot_NEW

riff_kakī_wild.vs.captive_rel.abund.ind_genus_plot_NEW <- inner_join(genus_kaki_wildcap_rel_abund_ind, genus_kaki_wildcap_pool_ind, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, wild.captive, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  mutate(taxon = fct_relevel(taxon,"Other","*Sanguibacter-Flavimobilis*","*Staphylococcus*","*Pseudoalteromonas*","*Stenotrophomonas*","*Janthinobacterium*","*Halomonas*","*Ralstonia*","*Pseudomonas*")) %>%
  inner_join(., genus_kaki_wildcap_sample_order, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id),
         sample_id = fct_reorder(sample_id, order)) %>%
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(width = 1) +
  scale_fill_manual(name="Genus",
                    breaks=c("*Pseudomonas*","*Ralstonia*","*Halomonas*","*Janthinobacterium*","*Stenotrophomonas*","*Pseudoalteromonas*","*Staphylococcus*","*Sanguibacter-Flavimobilis*","Other"),
                    values=c("#882255","#999933","#44AA99","#AA4499","#332288","#117733","#DDCC77","#CC6677","#999999")) +
  scale_y_continuous(expand=c(0, 0)) +
  facet_grid(~wild.captive, scale = "free_x", space = "free", switch="x",
             labeller = labeller(wild.captive = pretty_kaki_wildcap_label)) +
  labs(x="Egg Captive- or Wild-laid",
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(size=12),
        legend.title = element_text(size=14),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.justification = "left",
        strip.text.x = element_markdown(size = 12),
        axis.title=element_text(size=12))
riff_kakī_wild.vs.captive_rel.abund.ind_genus_plot_NEW

kakī_wild.vs.captive_relabund_all_NEW <- (riff_kakī_wild.vs.captive_rel.abund_phyla_plot_NEW | riff_kakī_wild.vs.captive_rel.abund.ind_phyla_plot_NEW) /
  (riff_kakī_wild.vs.captive_rel.abund_family_plot_NEW | riff_kakī_wild.vs.captive_rel.abund.ind_family_plot_NEW) /
  (riff_kakī_wild.vs.captive_rel.abund_genus_plot_NEW | riff_kakī_wild.vs.captive_rel.abund.ind_genus_plot_NEW)

ggview(kakī_wild.vs.captive_relabund_all_NEW, units = "px", height = 4000, width = 3000)
kakī_wild.vs.captive_relabund_all_NEW

#fake_df <- data.frame(letter = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T"),
#                              Counts = c(20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1))
#fake_df_4 <- fake_df[fake_df$Counts>=17,]
#fake_df_5 <- fake_df[fake_df$Counts>=16,]
#fake_df_6 <- fake_df[fake_df$Counts>=15,]
#fake_df_7 <- fake_df[fake_df$Counts>=14,]
#fake_df_8 <- fake_df[fake_df$Counts>=13,]
#ggplot(fake_df_8, aes(x = letter, y = Counts, fill = letter)) +
#  geom_bar(stat = "identity") +
#  scale_fill_grafify(palette = "safe",ColSeq = TRUE)

adonis_figure_data <- read.csv("adonis_figure_data.csv")
adonis_figure_data_hihi <- adonis_figure_data %>%
  filter(adonis_figure_data$species=="hihi")
adonis_figure_data_kaki <- adonis_figure_data %>%
  filter(adonis_figure_data$species=="kakī")

adonis_figure_data_hihi$beta_metric_order = adonis_figure_data_hihi$beta_metric
adonis_figure_data_hihi$beta_metric_order = factor(adonis_figure_data_hihi$beta_metric_order, levels=c("Weighted UniFrac","Unweighted UniFrac","Jaccard","Bray-Curtis"))
adonis_figure_hihi <- ggplot(adonis_figure_data_hihi, aes(x = fct_rev(factor(factor, level = c("Clutch/dam/sire ID","Dam experienced vs. naïve","Sire age","Sire age^2","Sire total social partners","Sire total genetic partners","Mean daily temp.","Mean daily rainfall","Overall min temp.","Overall max temp.","Sequencing plate ID","Freezer delay"))), y = R.2)) +
  geom_bar(stat = "identity", fill = "#37abc8") + coord_flip() +
  facet_grid(~as.factor(beta_metric_order)) +
  ylab(bquote(R^2)) +
  xlab(NULL) + 
  ylim(0,0.4) +
  geom_vline(xintercept = c(2.5,6.5),linetype="dotted") +
  scale_x_discrete(labels = c("Freezer delay","Sequencing plate ID","Overall max temp.","Overall min temp.","Mean daily rainfall","Mean daily temp.","Sire total genetic partners","Sire total social partners",expression("Sire age"^2),"Sire age","Dam experienced vs. naïve","Clutch/dam/sire ID")) +
  geom_text(aes(label = significance_star), hjust = -0.2, vjust = 0.65) +
  theme_classic() +
  theme(strip.text.x = element_text(face = "bold"))
adonis_figure_hihi

adonis_figure_data_kaki$beta_metric_order = adonis_figure_data_kaki$beta_metric
adonis_figure_data_kaki$beta_metric_order = factor(adonis_figure_data_kaki$beta_metric_order, levels=c("Weighted UniFrac","Unweighted UniFrac","Jaccard","Bray-Curtis"))
adonis_figure_kaki <- ggplot(adonis_figure_data_kaki, aes(x = fct_rev(factor(factor, level = c("Clutch ID","Dam/sire ID","Sire experienced vs. naïve","Sire age","Sire age^2","Sire total social partners","Clutch size","Clutch number","Est. days parental incubation at swab","Mean daily temp.","Mean daily rainfall","Overall min temp.","Overall max temp.","Nest location","Egg captive- vs. wild-laid","Sire captive- or wild-laid","Freezer delay"))), y = R.2)) +
  geom_bar(stat = "identity", fill = "#ff9955") + coord_flip() +
  facet_grid(~as.factor(beta_metric_order)) +
  ylab(bquote(R^2)) +
  xlab(NULL) + 
  ylim(0,0.8) +
  geom_vline(xintercept = c(1.5,3.5,8.5),linetype="dotted") +
  scale_x_discrete(labels = c("Freezer delay","Sire captive- or wild-laid","Egg captive- vs. wild-laid","Nest location","Overall max temp.","Overall min temp.","Mean daily rainfall","Mean daily temp.","Est. days parental incubation at swab","Clutch number","Clutch size","Sire total social partners",expression("Sire age"^2),"Sire age","Sire experienced vs. naïve","Dam/sire ID","Clutch ID")) +
  geom_text(aes(label = significance_star), hjust = -0.2, vjust = 0.65) +
  theme_classic() +
  theme(strip.text.x = element_text(face = "bold"))
adonis_figure_kaki

adonis_figure_hihi +
  adonis_figure_kaki + plot_layout(ncol = 1, heights = c(1, 1.3))

##### ALL PLOTS #####
phyla_plot
family_plot
genera_plot
phyla_plot_top5
family_plot_top5
genera_plot_top5
ASV_deseq2_plot
genus_deseq2_plot
family_deseq2_plot
phylum_deseq2_plot
ASV_volcano
genus_volcano
family_volcano
phylum_volcano
riff_hihi.vs.kakī_rel.abund_phyla_plot
riff_hihi.vs.kakī_rel.abund_family_plot
riff_hihi.vs.kakī_rel.abund_genus_plot
riff_hihi_hatched.vs.unhatched_rel.abund_phyla_plot
riff_hihi_hatched.vs.unhatched_rel.abund_family_plot
riff_hihi_hatched.vs.unhatched_rel.abund_genus_plot
riff_kakī_hatched.vs.unhatched_rel.abund_phyla_plot
riff_kakī_hatched.vs.unhatched_rel.abund_family_plot
riff_kakī_hatched.vs.unhatched_rel.abund_genus_plot
riff_kakī_wild.vs.captive_rel.abund_phyla_plot
riff_kakī_wild.vs.captive_rel.abund_family_plot
riff_kakī_wild.vs.captive_rel.abund_genus_plot
riff_hihi.vs.kakī_rel.abund.ind_phyla_plot
riff_hihi.vs.kakī_rel.abund.ind_family_plot
riff_hihi.vs.kakī_rel.abund.ind_genus_plot
riff_hihi_hatched.vs.unhatched_rel.abund.ind_phyla_plot
riff_hihi_hatched.vs.unhatched_rel.abund.ind_family_plot
riff_hihi_hatched.vs.unhatched_rel.abund.ind_genus_plot
riff_kakī_hatched.vs.unhatched_rel.abund.ind_phyla_plot
riff_kakī_hatched.vs.unhatched_rel.abund.ind_family_plot
riff_kakī_hatched.vs.unhatched_rel.abund.ind_genus_plot
riff_kakī_wild.vs.captive_rel.abund.ind_phyla_plot
riff_kakī_wild.vs.captive_rel.abund.ind_family_plot
riff_kakī_wild.vs.captive_rel.abund.ind_genus_plot
