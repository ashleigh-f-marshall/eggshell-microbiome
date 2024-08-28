####################
# Authors: Ashleigh Fleming Marshall
# ashleigh.marshall@ioz.ac.uk | ashleigh.marshall.16@ucl.ac.uk
####################

### MICROBIOME, DISEASE, AND HATCHING FAILURE

##Clear workspace
rm(list=ls())

#### Install and load packages ####

library("phyloseq"); packageVersion("phyloseq") #‘1.40.0’
library("vegan"); packageVersion("vegan") #‘2.6.4’’
library("tidyverse"); packageVersion("tidyverse") #'1.3.2'
library("ggplot2"); packageVersion("ggplot2") #'3.4.0'
library("ggtext"); packageVersion("ggtext") #'0.1.2'
library("ggview"); packageVersion("ggview") #‘0.1.0’
library("lme4"); packageVersion("lme4") #‘1.1.31’
library("MuMIn"); packageVersion("MuMIn") #‘1.47.1’
library("rcompanion"); packageVersion("rcompanion") #‘2.4.18’
library("greybox"); packageVersion("greybox") #‘1.0.7’
library("ltm"); packageVersion("ltm") #‘1.2.0’
library("performance"); packageVersion("performance") #'0.10.2'
library("AICcmodavg"); packageVersion("AICcmodavg") #‘2.3.1’

setwd("Data/DataAnalyses")

#### Read in and process data ####

mtd <- read_tsv(("sample-metadata-new-noNA-removemock_forR_wEnvironmentalData_2.tsv"), comment = "#q2")
colnames(mtd)[29] <- "female_age_cohort.sq"
colnames(mtd)[45] <- "male_age_cohort.sq"
colnames(mtd)[62] <- "first_egg_lay"
colnames(mtd)[79] <- "re.est_first_egg_lay"
colnames(mtd)[89] <- "clutch.hatch.percent"
mtd_short <- mtd %>%
  dplyr::select(SampleID,plate_ID,species,wild.captive,clutch.nest_id_corrected,nest_location,clutch_number,clutch_size,hatched_total,female_combo,female_age_cohort,female_age_cohort.sq,female_lay_wild.captive,female_experienced.naive_social,female_total.partners_social,female_total.partners_genetic,Male_combo,male_age_cohort,male_age_cohort.sq,male_lay_wild.captive,male_experienced.naive_social,male_total.partners_social,male_total.partners_genetic,first_egg_lay,hatched,re.est_first_egg_lay,hatched.unhatched,failure_stage,COD,development_day.13,clutch.hatch.percent,freezer_delay,swab_day_parentalInc,mean_temp_laytoswab,min_temp_laytoswab,max_temp_laytoswab,mean_rain_laytoswab)
				
mtd_hihi <- mtd %>%
  filter(species=="Hihi")
mtd_kaki <- mtd %>%
  filter(species=="Kaki")

mtd_kaki %>%
  count(Male_combo)

mtd_short_hihi <- mtd_short %>%
  filter(species=="Hihi")
mtd_short_kaki <- mtd_short %>%
  filter(species=="Kaki")

shannon_hihi <- read_tsv("hihi_exported_diversity_metrics/shannon_vector/alpha-diversity.tsv")
colnames(shannon_hihi)[1] <- "SampleID"
evenness_hihi <- read_tsv("hihi_exported_diversity_metrics/evenness_vector/alpha-diversity.tsv")
colnames(evenness_hihi)[1] <- "SampleID"
faith_pd_hihi <- read_tsv("hihi_exported_diversity_metrics/faith_pd_vector/alpha-diversity.tsv")
colnames(faith_pd_hihi)[1] <- "SampleID"
observed_features_hihi <- read_tsv("hihi_exported_diversity_metrics/observed_features_vector/alpha-diversity.tsv")
colnames(observed_features_hihi)[1] <- "SampleID"

hihi_alpha <- shannon_hihi %>%
  inner_join(., evenness_hihi, by = "SampleID") %>%
  inner_join(., faith_pd_hihi, by = "SampleID") %>%
  inner_join(., observed_features_hihi, by = "SampleID")
colnames(hihi_alpha)

shannon_kaki <- read_tsv("kaki_exported_diversity_metrics/shannon_vector/alpha-diversity.tsv")
colnames(shannon_kaki)[1] <- "SampleID"
evenness_kaki <- read_tsv("kaki_exported_diversity_metrics/evenness_vector/alpha-diversity.tsv")
colnames(evenness_kaki)[1] <- "SampleID"
faith_pd_kaki <- read_tsv("kaki_exported_diversity_metrics/faith_pd_vector/alpha-diversity.tsv")
colnames(faith_pd_kaki)[1] <- "SampleID"
observed_features_kaki <- read_tsv("kaki_exported_diversity_metrics/observed_features_vector/alpha-diversity.tsv")
colnames(observed_features_kaki)[1] <- "SampleID"

kaki_alpha <- shannon_kaki %>%
  inner_join(., evenness_kaki, by = "SampleID") %>%
  inner_join(., faith_pd_kaki, by = "SampleID") %>%
  inner_join(., observed_features_kaki, by = "SampleID")
colnames(kaki_alpha)

genus_rel_abund_ind <- readRDS("genus_rel_abund_ind.rds")
family_rel_abund_ind <- readRDS("family_rel_abund_ind.rds")
colnames(genus_rel_abund_ind)[1] <- "species"
colnames(genus_rel_abund_ind)[2] <- "SampleID"
colnames(family_rel_abund_ind)[1] <- "species"
colnames(family_rel_abund_ind)[2] <- "SampleID"
genus_rel_abund_ind_hihi <- genus_rel_abund_ind %>%
  filter(species == "Hihi")
genus_rel_abund_ind_kaki <- genus_rel_abund_ind %>%
  filter(species == "Kaki")
family_rel_abund_ind_hihi <- family_rel_abund_ind %>%
  filter(species == "Hihi")
family_rel_abund_ind_kaki <- family_rel_abund_ind %>%
  filter(species == "Kakī")

hihi_rel_abund_Staphylococcus <-  filter(genus_rel_abund_ind_hihi, taxon == "Staphylococcus")
colnames(hihi_rel_abund_Staphylococcus)[4] <- "rel_abund_Staphylococcus"
hihi_rel_abund_Staphylococcus <- dplyr::select(hihi_rel_abund_Staphylococcus, -c(taxon,species))
hihi_rel_abund_Streptococcus <-  filter(genus_rel_abund_ind_hihi, taxon == "Streptococcus")
colnames(hihi_rel_abund_Streptococcus)[4] <- "rel_abund_Streptococcus"
hihi_rel_abund_Streptococcus <- dplyr::select(hihi_rel_abund_Streptococcus, -c(taxon,species))
hihi_rel_abund_Enterococcus <-  filter(genus_rel_abund_ind_hihi, taxon == "Enterococcus")
colnames(hihi_rel_abund_Enterococcus)[4] <- "rel_abund_Enterococcus"
hihi_rel_abund_Enterococcus <- dplyr::select(hihi_rel_abund_Enterococcus, -c(taxon,species))
hihi_rel_abund_Campylobacter <-  filter(genus_rel_abund_ind_hihi, taxon == "Campylobacter")
colnames(hihi_rel_abund_Campylobacter)[4] <- "rel_abund_Campylobacter"
hihi_rel_abund_Campylobacter <- dplyr::select(hihi_rel_abund_Campylobacter, -c(taxon,species))
hihi_rel_abund_Mycoplasma <-  filter(genus_rel_abund_ind_hihi, taxon == "Mycoplasma")
colnames(hihi_rel_abund_Mycoplasma)[4] <- "rel_abund_Mycoplasma"
hihi_rel_abund_Mycoplasma <- dplyr::select(hihi_rel_abund_Mycoplasma, -c(taxon,species))
hihi_rel_abund_Escherichia <-  filter(genus_rel_abund_ind_hihi, taxon == "Escherichia-Shigella")
colnames(hihi_rel_abund_Escherichia)[4] <- "rel_abund_Escherichia"
hihi_rel_abund_Escherichia <- dplyr::select(hihi_rel_abund_Escherichia, -c(taxon,species))
hihi_rel_abund_Pseudomonas <-  filter(genus_rel_abund_ind_hihi, taxon == "Pseudomonas")
colnames(hihi_rel_abund_Pseudomonas)[4] <- "rel_abund_Pseudomonas"
hihi_rel_abund_Pseudomonas <- dplyr::select(hihi_rel_abund_Pseudomonas, -c(taxon,species))
hihi_rel_abund_Neisseria <-  filter(genus_rel_abund_ind_hihi, taxon == "Neisseria")
colnames(hihi_rel_abund_Neisseria)[4] <- "rel_abund_Neisseria"
hihi_rel_abund_Neisseria <- dplyr::select(hihi_rel_abund_Neisseria, -c(taxon,species))
hihi_rel_abund_Enterobacteriaceae <-  filter(family_rel_abund_ind_hihi, taxon == "Enterobacteriaceae")
colnames(hihi_rel_abund_Enterobacteriaceae)[4] <- "rel_abund_Enterobacteriaceae"
hihi_rel_abund_Enterobacteriaceae <- dplyr::select(hihi_rel_abund_Enterobacteriaceae, -c(taxon,species))

hihi_rel_abund_data <- hihi_rel_abund_Staphylococcus %>%
  inner_join(., hihi_rel_abund_Streptococcus, by = "SampleID") %>%
  inner_join(., hihi_rel_abund_Enterococcus, by = "SampleID") %>%
  inner_join(., hihi_rel_abund_Campylobacter, by = "SampleID") %>%
  inner_join(., hihi_rel_abund_Mycoplasma, by = "SampleID") %>%
  inner_join(., hihi_rel_abund_Escherichia, by = "SampleID") %>%
  inner_join(., hihi_rel_abund_Pseudomonas, by = "SampleID") %>%
  inner_join(., hihi_rel_abund_Neisseria, by = "SampleID") %>%
  inner_join(., hihi_rel_abund_Enterobacteriaceae, by = "SampleID")
colnames(hihi_rel_abund_data)

kaki_rel_abund_Staphylococcus <-  filter(genus_rel_abund_ind_kaki, taxon == "Staphylococcus")
colnames(kaki_rel_abund_Staphylococcus)[4] <- "rel_abund_Staphylococcus"
kaki_rel_abund_Staphylococcus <- dplyr::select(kaki_rel_abund_Staphylococcus, -c(taxon,species))
kaki_rel_abund_Streptococcus <-  filter(genus_rel_abund_ind_kaki, taxon == "Streptococcus")
colnames(kaki_rel_abund_Streptococcus)[4] <- "rel_abund_Streptococcus"
kaki_rel_abund_Streptococcus <- dplyr::select(kaki_rel_abund_Streptococcus, -c(taxon,species))
kaki_rel_abund_Enterococcus <-  filter(genus_rel_abund_ind_kaki, taxon == "Enterococcus")
colnames(kaki_rel_abund_Enterococcus)[4] <- "rel_abund_Enterococcus"
kaki_rel_abund_Enterococcus <- dplyr::select(kaki_rel_abund_Enterococcus, -c(taxon,species))
kaki_rel_abund_Campylobacter <-  filter(genus_rel_abund_ind_kaki, taxon == "Campylobacter")
colnames(kaki_rel_abund_Campylobacter)[4] <- "rel_abund_Campylobacter"
kaki_rel_abund_Campylobacter <- dplyr::select(kaki_rel_abund_Campylobacter, -c(taxon,species))
kaki_rel_abund_Mycoplasma <-  filter(genus_rel_abund_ind_kaki, taxon == "Mycoplasma")
colnames(kaki_rel_abund_Mycoplasma)[4] <- "rel_abund_Mycoplasma"
kaki_rel_abund_Mycoplasma <- dplyr::select(kaki_rel_abund_Mycoplasma, -c(taxon,species))
kaki_rel_abund_Escherichia <-  filter(genus_rel_abund_ind_kaki, taxon == "Escherichia-Shigella")
colnames(kaki_rel_abund_Escherichia)[4] <- "rel_abund_Escherichia"
kaki_rel_abund_Escherichia <- dplyr::select(kaki_rel_abund_Escherichia, -c(taxon,species))
kaki_rel_abund_Pseudomonas <-  filter(genus_rel_abund_ind_kaki, taxon == "Pseudomonas")
colnames(kaki_rel_abund_Pseudomonas)[4] <- "rel_abund_Pseudomonas"
kaki_rel_abund_Pseudomonas <- dplyr::select(kaki_rel_abund_Pseudomonas, -c(taxon,species))
kaki_rel_abund_Neisseria <-  filter(genus_rel_abund_ind_kaki, taxon == "Neisseria")
colnames(kaki_rel_abund_Neisseria)[4] <- "rel_abund_Neisseria"
kaki_rel_abund_Neisseria <- dplyr::select(kaki_rel_abund_Neisseria, -c(taxon,species))
kaki_rel_abund_Enterobacteriaceae <-  filter(family_rel_abund_ind_kaki, taxon == "Enterobacteriaceae")
colnames(kaki_rel_abund_Enterobacteriaceae)[4] <- "rel_abund_Enterobacteriaceae"
kaki_rel_abund_Enterobacteriaceae <- dplyr::select(kaki_rel_abund_Enterobacteriaceae, -c(taxon,species))

kaki_rel_abund_data <- kaki_rel_abund_Staphylococcus %>%
  inner_join(., kaki_rel_abund_Streptococcus, by = "SampleID") %>%
  inner_join(., kaki_rel_abund_Enterococcus, by = "SampleID") %>%
  inner_join(., kaki_rel_abund_Campylobacter, by = "SampleID") %>%
  inner_join(., kaki_rel_abund_Mycoplasma, by = "SampleID") %>%
  inner_join(., kaki_rel_abund_Escherichia, by = "SampleID") %>%
  inner_join(., kaki_rel_abund_Pseudomonas, by = "SampleID") %>%
  inner_join(., kaki_rel_abund_Neisseria, by = "SampleID") %>%
  inner_join(., kaki_rel_abund_Enterobacteriaceae, by = "SampleID")
colnames(kaki_rel_abund_data)

hihi_data <- mtd_short_hihi %>%
  inner_join(., hihi_alpha, by = "SampleID") %>%
  inner_join(., hihi_rel_abund_data, by = "SampleID")
colnames(hihi_data)

kaki_data <- mtd_short_kaki %>%
  inner_join(., kaki_alpha, by = "SampleID") %>%
  inner_join(., kaki_rel_abund_data, by = "SampleID")
colnames(kaki_data)

kaki_collection_weight <- read.csv("kaki_collection_weight.csv")

kaki_data <- inner_join(kaki_data, kaki_collection_weight, by = "SampleID")
colnames(kaki_data)

hihi_data$clutch.hatch.prop <- hihi_data$hatched_total/hihi_data$clutch_size
kaki_data$clutch.hatch.prop <- kaki_data$hatched_total/kaki_data$clutch_size

hihi_data %>%
  count(hatched.unhatched) #84 hatched, 22 unhatched
(84/(84+22))*100 #79.24528
kaki_data %>%
  count(hatched.unhatched) #73 hatched, 19 unhatched
(73/(73+19))*100 #79.34783
kaki_data[kaki_data$wild.captive=="Wild",] %>%
  count(hatched.unhatched) #52 hatched, 17 unhatched
(52/(52+17))*100 #75.36232
kaki_data[kaki_data$wild.captive=="Captive",] %>%
  count(hatched.unhatched) #21 hatched, 2 unhatched
(21/(21+2))*100 #91.30435
mean(hihi_data$clutch.hatch.prop)*100 #75.15924
mean(kaki_data$clutch.hatch.prop)*100 #77.17391
mean(kaki_data[kaki_data$wild.captive=="Wild",]$clutch.hatch.prop)*100 #72.46377
mean(kaki_data[kaki_data$wild.captive=="Captive",]$clutch.hatch.prop)*100 #91.30435

res_1 <- prop.test(x = c(52,21), n = c(69,23), alternative = "two.sided", correct = TRUE)
res_1 #p = 0.1808
sqrt(res_1$statistic) #1.338251 = z score
res_2 <- prop.test(x = c(52,21), n = c(69,23), alternative = "two.sided", correct = FALSE)
res_2 #p = 0.1019
sqrt(res_2$statistic) #1.63564

#### Simple regression of kaki incubation length and relative abundance of Pseudomonas ####

hist(kaki_data$rel_abund_Pseudomonas)
ggplot(kaki_data,aes(swab_day_parentalInc, rel_abund_Pseudomonas)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x = "Estimated days of parental incubation at swab", y = expression(paste("Relative abundance of ",italic("Pseudomonas")))) +
  theme_minimal()
cor(kaki_data$rel_abund_Pseudomonas,kaki_data$swab_day_parentalInc) #-0.4723612
cor.test(kaki_data$rel_abund_Pseudomonas,kaki_data$swab_day_parentalInc) #-0.4723612, p-value = 1.995e-06 = low-moderate
cor(kaki_data$rel_abund_Pseudomonas,kaki_data$swab_day_parentalInc, method = "spearman") #-0.2705131
kaki_incubation_pseudomonas <- lm(rel_abund_Pseudomonas ~ swab_day_parentalInc, data = kaki_data)
plot(kaki_incubation_pseudomonas)
summary(kaki_incubation_pseudomonas)

#### Check for correlations between variables ####

## Continuous vs. continuous = Pearson Correlation
# Dam age vs. sire age - Hihi
cor.test(hihi_data$female_age_cohort, hihi_data$male_age_cohort) #0.0395515, p-value = 0.6229 = negligible correlation (not sig?)
cor.test(hihi_data$female_age_cohort.sq, hihi_data$male_age_cohort.sq) #0.003155689, p-value = 0.9687 = negligible correlation (not sig?)
plot(hihi_data$female_age_cohort, hihi_data$male_age_cohort)
plot(hihi_data$female_age_cohort.sq, hihi_data$male_age_cohort.sq)

# Dam age vs. sire age - Kaki
cor.test(kaki_data$female_age_cohort, kaki_data$male_age_cohort) #0.6057583, p-value = 4.027e-10 = moderate correlation
cor.test(kaki_data$female_age_cohort.sq, kaki_data$male_age_cohort.sq) #0.4796891, p-value = 2.259e-06 = low-moderate correlation
plot(kaki_data$female_age_cohort, kaki_data$male_age_cohort)
plot(kaki_data$female_age_cohort.sq, kaki_data$male_age_cohort.sq)

# Dam age vs. dam total social partners - Hihi
cor.test(hihi_data$female_age_cohort, hihi_data$female_total.partners_social) #0.8045378, p-value < 2.2e-16 = high correlation
cor.test(hihi_data$female_age_cohort.sq, hihi_data$female_total.partners_social) #0.8084775, p-value < 2.2e-16 = high correlation

# Dam age vs. dam total social partners - Kaki
cor.test(kaki_data$female_age_cohort, kaki_data$female_total.partners_social) #0.77124, p-value < 2.2e-16 = high correlation
cor.test(kaki_data$female_age_cohort.sq, kaki_data$female_total.partners_social) #0.804769, p-value < 2.2e-16 = high correlation

# Dam age vs. dam total genetic partners - Hihi
cor.test(hihi_data$female_age_cohort, hihi_data$female_total.partners_genetic) #0.9340186, p-value < 2.2e-16 = very high correlation
cor.test(hihi_data$female_age_cohort.sq, hihi_data$female_total.partners_genetic) #0.9168213, p-value < 2.2e-16 = very high correlation

# Sire age vs. Sire total social partners - Hihi
cor.test(hihi_data$male_age_cohort, hihi_data$male_total.partners_social) #0.7396615, p-value < 2.2e-16 = high correlation
cor.test(hihi_data$male_age_cohort.sq, hihi_data$male_total.partners_social) #0.8037683, p-value < 2.2e-16 = high correlation

# Sire age vs. Sire total social partners - Kaki
cor.test(kaki_data$male_age_cohort, kaki_data$male_total.partners_social) #0.609879, p-value = 1.097e-10 = moderate correlation
cor.test(kaki_data$male_age_cohort.sq, kaki_data$male_total.partners_social) #0.6082335, p-value = 1.27e-10 = moderate correlation

# Sire age vs. Sire total genetic partners - Hihi
cor.test(hihi_data$male_age_cohort, hihi_data$male_total.partners_genetic) #0.8025488, p-value < 2.2e-16 = high correlation
cor.test(hihi_data$male_age_cohort.sq, hihi_data$male_total.partners_genetic) #0.7941576, p-value < 2.2e-16 = high correlation

# Egg weight vs. days of parental incubation at collection - Kaki
cor.test(kaki_data$collection_egg_weight, kaki_data$swab_day_parentalInc) #-0.2683281, p-value = 0.01101 = weak correlation

# Temperature variables - Hihi
cor.test(hihi_data$min_temp_laytoswab, hihi_data$max_temp_laytoswab) #0.6821413, p-value < 2.2e-16 = moderate correlation
cor.test(hihi_data$min_temp_laytoswab, hihi_data$mean_temp_laytoswab) #0.9698433, p-value < 2.2e-16 = very high correlation
cor.test(hihi_data$max_temp_laytoswab, hihi_data$mean_temp_laytoswab) #0.7627425, p-value < 2.2e-16 = high correlation

# Temperature variables - Kaki
cor.test(kaki_data$min_temp_laytoswab, kaki_data$max_temp_laytoswab) #0.1612957, p-value = 0.1245 = weak correlation
cor.test(kaki_data$min_temp_laytoswab, kaki_data$mean_temp_laytoswab) #0.7247286, p-value = 3.169e-16 = moderate correlation
cor.test(kaki_data$max_temp_laytoswab, kaki_data$mean_temp_laytoswab) #0.7636094, p-value < 2.2e-16 = high correlation

## Categorical vs. categorical = Chi-square Tests (significance and Cramer's V values)
# Dam experienced/naive vs. sire experienced.naive - Hihi
tbl_1 <-table(hihi_data$female_experienced.naive_social,hihi_data$male_experienced.naive_social)
chi2_1 <- chisq.test(tbl_1, correct=F)
chi2_1 # p < 0.001 = not independent
sqrt(chi2_1$statistic/sum(tbl_1)) #Cramer's V = 0.3532031
cramerV(hihi_data$female_experienced.naive_social,hihi_data$male_experienced.naive_social) #Cramer V = 0.3532 = moderate association
cramer(hihi_data$female_experienced.naive_social,hihi_data$male_experienced.naive_social) #Cramer's V = 0.3318 = moderate association

# Dam experienced/naive vs. sire experienced.naive - Kaki
tbl_2 <-table(kaki_data$female_experienced.naive_social,kaki_data$male_experienced.naive_social)
chi2_2 <- chisq.test(tbl_2, correct=F)
chi2_2 # p < 0.001 = not independent
sqrt(chi2_2$statistic/sum(tbl_2)) #Cramer's V = 0.719195
cramerV(kaki_data$female_experienced.naive_social,kaki_data$male_experienced.naive_social) #Cramer V = 0.7034 = strong association
cramer(kaki_data$female_experienced.naive_social,kaki_data$male_experienced.naive_social) #Cramer's V = 0.6909 = strong association

# Nest location vs. female ID - Kaki
tbl_3 <-table(kaki_data$nest_location,kaki_data$female_combo)
chi2_3 <- chisq.test(tbl_3, correct=F)
chi2_3 # p-value < 2.2e-16 = not independent
sqrt(chi2_3$statistic/sum(tbl_3)) #Cramer's V = 3.510916
cramerV(kaki_data$nest_location,kaki_data$female_combo) #Cramer V = 0.9738 = strong association
cramer(kaki_data$nest_location,kaki_data$female_combo) #Cramer's V = 0.9288 = strong association

# Egg captive/wild vs. dam captive/wild laid - Kaki
tbl_4 <-table(kaki_data$wild.captive,kaki_data$female_lay_wild.captive)
chi2_4 <- chisq.test(tbl_4, correct=F)
chi2_4 # p-value = 0.0001379 = not independent
sqrt(chi2_4$statistic/sum(tbl_4)) #Cramer's V = 0.4063597
cramerV(kaki_data$wild.captive,kaki_data$female_lay_wild.captive) #Cramer V = 0.3974 = moderate association
cramer(kaki_data$wild.captive,kaki_data$female_lay_wild.captive) #Cramer's V = 0.3652 = moderate association

# Egg captive/wild vs. sire captive/wild laid - Kaki
tbl_5 <-table(kaki_data$wild.captive,kaki_data$male_lay_wild.captive)
chi2_5 <- chisq.test(tbl_5, correct=F)
chi2_5 # p-value = 7.889e-05 = not independent
sqrt(chi2_5$statistic/sum(tbl_5)) #Cramer's V = 0.411581
cramerV(kaki_data$wild.captive,kaki_data$male_lay_wild.captive) #Cramer V = 0.4116 = moderate association
cramer(kaki_data$wild.captive,kaki_data$male_lay_wild.captive) #Cramer's V = 0.3725 = moderate association

# Dam captive/wild laid vs. sire captive/wild laid - Kaki
tbl_6 <-table(kaki_data$female_lay_wild.captive,kaki_data$male_lay_wild.captive)
chi2_6 <- chisq.test(tbl_6, correct=F)
chi2_6 # p-value = 1.942e-05 = not independent
sqrt(chi2_6$statistic/sum(tbl_6)) #Cramer's V = 0.4553419
cramerV(kaki_data$female_lay_wild.captive,kaki_data$male_lay_wild.captive) #Cramer V = 0.4453 = moderate association
cramer(kaki_data$female_lay_wild.captive,kaki_data$male_lay_wild.captive) #Cramer's V = 0.4182 = moderate association

## Continuous vs. categorical = Point biserial correlation
# Dam age vs. dam experience - Hihi
hihi_data$dam_exp_dummy <- ifelse(hihi_data$female_experienced.naive_social == "Experienced", 1, 0)
cor.test(hihi_data$dam_exp_dummy, hihi_data$female_age_cohort) #0.755615, p-value < 2.2e-16 = high correlation
cor.test(hihi_data$dam_exp_dummy, hihi_data$female_age_cohort.sq) #0.6078867, p-value < 2.2e-16 = moderate correlation

# Sire age vs. sire experience - Hihi
hihi_data$sire_exp_dummy <- ifelse(hihi_data$male_experienced.naive_social == "Experienced", 1, 0)
cor.test(hihi_data$sire_exp_dummy, hihi_data$male_age_cohort) #0.6630874, p-value < 2.2e-16 = moderate correlation
cor.test(hihi_data$sire_exp_dummy, hihi_data$male_age_cohort.sq) #0.5127647, p-value = 6.618e-12 = moderate correlation

# Dam age vs. dam experience - Kaki
kaki_data$dam_exp_dummy <- ifelse(kaki_data$female_experienced.naive_social == "Experienced", 1, 0)
cor.test(kaki_data$dam_exp_dummy, kaki_data$female_age_cohort) #0.506982, p-value < 4.644e-07 = moderate correlation
cor.test(kaki_data$dam_exp_dummy, kaki_data$female_age_cohort.sq) #0.3700054, p-value = 0.000388 = low correlation

# Sire age vs. sire experience - Kaki
kaki_data$sire_exp_dummy <- ifelse(kaki_data$male_experienced.naive_social == "Experienced", 1, 0)
cor.test(kaki_data$sire_exp_dummy, kaki_data$male_age_cohort) #0.7088779, p-value < 2.654e-15 = high correlation
cor.test(kaki_data$sire_exp_dummy, kaki_data$male_age_cohort.sq) #0.6018522, p-value = 2.219e-10 = moderate correlation

## Running GLMMs
colnames(hihi_data)

## Hihi

hihi_data_bin_NArm <- na.omit(hihi_data[,c("hatched","female_combo","shannon_entropy","pielou_evenness","faith_pd","observed_features","rel_abund_Staphylococcus","rel_abund_Enterococcus","rel_abund_Pseudomonas","rel_abund_Enterobacteriaceae","clutch_size","female_age_cohort","female_age_cohort.sq","male_age_cohort","male_age_cohort.sq","mean_temp_laytoswab","mean_rain_laytoswab")])
hihi_data_prop_NArm <- na.omit(hihi_data[,c("clutch.hatch.prop","female_combo","shannon_entropy","pielou_evenness","faith_pd","observed_features","rel_abund_Staphylococcus","rel_abund_Enterococcus","rel_abund_Pseudomonas","rel_abund_Enterobacteriaceae","clutch_size","female_age_cohort","female_age_cohort.sq","male_age_cohort","male_age_cohort.sq","clutch_size","mean_temp_laytoswab","mean_rain_laytoswab")])

# Hihi binary hatching success
hihi_bin_model_null <- glmer(hatched ~ 1 + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))

hihi_bin_model_shannon <- glmer(hatched ~ scale(shannon_entropy) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_shannon)
hihi_bin_model_pielou <- glmer(hatched ~ scale(pielou_evenness) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_pielou)
hihi_bin_model_faith <- glmer(hatched ~ scale(faith_pd) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_faith)
hihi_bin_model_observed <- glmer(hatched ~ scale(observed_features) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_observed) 
hihi_bin_model_Staphylococcus <- glmer(hatched ~ scale(rel_abund_Staphylococcus) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_Staphylococcus) 
hihi_bin_model_Enterococcus <- glmer(hatched ~ scale(rel_abund_Enterococcus) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_Enterococcus) 
hihi_bin_model_Pseudomonas <- glmer(hatched ~ scale(rel_abund_Pseudomonas) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_Pseudomonas) 
hihi_bin_model_Enterobacteriaceae <- glmer(hatched ~ scale(rel_abund_Enterobacteriaceae) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_Enterobacteriaceae) 
hihi_bin_model_clutch.size <- glmer(hatched ~ scale(clutch_size) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_clutch.size) 
hihi_bin_model_dam.age <- glmer(hatched ~ scale(female_age_cohort) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_dam.age)
hihi_bin_model_dam.age.sq <- glmer(hatched ~ scale(female_age_cohort.sq) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_dam.age.sq) 
hihi_bin_model_sire.age <- glmer(hatched ~ scale(male_age_cohort) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_sire.age) 
hihi_bin_model_sire.age.sq <- glmer(hatched ~ scale(male_age_cohort.sq) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_sire.age.sq) 
hihi_bin_model_mean.temp <- glmer(hatched ~ scale(mean_temp_laytoswab) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_mean.temp) 
hihi_bin_model_mean.rain <- glmer(hatched ~ scale(mean_rain_laytoswab) + (1|female_combo),data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_mean.rain) 

anova(hihi_bin_model_null,hihi_bin_model_shannon) # better = hihi_bin_model_shannon (significant) *
anova(hihi_bin_model_null,hihi_bin_model_pielou) # better = hihi_bin_model_pielou (not significant) --> near sig.
anova(hihi_bin_model_null,hihi_bin_model_faith) # better = hihi_bin_model_null (not significant)
anova(hihi_bin_model_null,hihi_bin_model_observed) # better = hihi_bin_model_null (not significant)
anova(hihi_bin_model_null,hihi_bin_model_Staphylococcus) # better = hihi_bin_model_Staphylococcus (not significant)
anova(hihi_bin_model_null,hihi_bin_model_Enterococcus) # better = hihi_bin_model_Enterococcus (not significant)
anova(hihi_bin_model_null,hihi_bin_model_Pseudomonas) # better = hihi_bin_model_null (not significant)
anova(hihi_bin_model_null,hihi_bin_model_Enterobacteriaceae) # better = hihi_bin_model_null (not significant)
anova(hihi_bin_model_null,hihi_bin_model_clutch.size) # better = hihi_bin_model_null (not significant)
anova(hihi_bin_model_null,hihi_bin_model_dam.age) # better = hihi_bin_model_dam.age (significant) *
anova(hihi_bin_model_null,hihi_bin_model_dam.age.sq) # better = hihi_bin_model_dam.age.sq (significant) *
anova(hihi_bin_model_null,hihi_bin_model_sire.age) # better = hihi_bin_model_null (not significiant)
anova(hihi_bin_model_null,hihi_bin_model_sire.age.sq) # better = hihi_bin_model_null (not significant)
anova(hihi_bin_model_null,hihi_bin_model_mean.temp) # better = hihi_bin_model_null (not significant)
anova(hihi_bin_model_null,hihi_bin_model_mean.rain) # better = hihi_bin_model_null (not significant)

anova(hihi_bin_model_dam.age, hihi_bin_model_dam.age.sq) # identical?

hihi_bin_model_full.sq <- glmer(hatched ~ scale(shannon_entropy) + scale(female_age_cohort.sq) + (1|female_combo), data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_full.sq)
check_overdispersion(hihi_bin_model_full.sq) # No overdispersion detected
hihi_bin_model_full <- glmer(hatched ~ scale(shannon_entropy) + scale(female_age_cohort) + (1|female_combo), data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_full)
check_overdispersion(hihi_bin_model_full) # No overdispersion detected

anova(hihi_bin_model_full,hihi_bin_model_full.sq) # better = hihi_bin_model_full.sq (not significant)
anova(hihi_bin_model_full.sq,hihi_bin_model_null) # better = hihi_bin_model_full.sq (significant)
anova(hihi_bin_model_full.sq,hihi_bin_model_shannon) # better = hihi_bin_model_full.sq (significant)
anova(hihi_bin_model_full.sq,hihi_bin_model_dam.age) # better = hihi_bin_model_full.sq (significant)
anova(hihi_bin_model_shannon,hihi_bin_model_dam.age) # better = hihi_bin_model_shannon (not significant)
anova(hihi_bin_model_full.sq,hihi_bin_model_dam.age.sq) # better = hihi_bin_model_full.sq (significant)
anova(hihi_bin_model_shannon,hihi_bin_model_dam.age.sq) # better = hihi_bin_model_shannon (not significant)

hihi_bin_model_int_1 <- glmer(hatched ~ scale(shannon_entropy) * scale(female_age_cohort.sq) + (1|female_combo), data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_int_1)
check_overdispersion(hihi_bin_model_int_1) # No overdispersion detected
hihi_bin_model_int_2 <- glmer(hatched ~ scale(shannon_entropy) * scale(female_age_cohort) + (1|female_combo), data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_int_2)
check_overdispersion(hihi_bin_model_int_2) # No overdispersion detected
hihi_bin_model_int_3 <- glmer(hatched ~ scale(shannon_entropy) * scale(mean_temp_laytoswab) + scale(female_age_cohort.sq) + (1|female_combo), data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_int_3)
check_overdispersion(hihi_bin_model_int_3) # No overdispersion detected
hihi_bin_model_int_4 <- glmer(hatched ~ scale(shannon_entropy) * scale(mean_temp_laytoswab) + scale(female_age_cohort) + (1|female_combo), data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_int_4)
check_overdispersion(hihi_bin_model_int_4) # No overdispersion detected
hihi_bin_model_int_5 <- glmer(hatched ~ scale(shannon_entropy) * scale(mean_rain_laytoswab) + scale(female_age_cohort.sq) + (1|female_combo), data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_int_5)
check_overdispersion(hihi_bin_model_int_5) # No overdispersion detected
hihi_bin_model_int_6 <- glmer(hatched ~ scale(shannon_entropy) * scale(mean_rain_laytoswab) + scale(female_age_cohort) + (1|female_combo), data=hihi_data_bin_NArm,family=binomial(link="logit"))
summary(hihi_bin_model_int_6)
check_overdispersion(hihi_bin_model_int_6) # No overdispersion detected

anova(hihi_bin_model_full,hihi_bin_model_int_1) #better = hihi_bin_model_full (not significant)
anova(hihi_bin_model_full,hihi_bin_model_int_2) #better = hihi_bin_model_full (not significant)
anova(hihi_bin_model_full,hihi_bin_model_int_3) #better = hihi_bin_model_full (not significant)
anova(hihi_bin_model_full,hihi_bin_model_int_4) #better = hihi_bin_model_full (not significant)
anova(hihi_bin_model_full,hihi_bin_model_int_5) #better = hihi_bin_model_full (not significant)
anova(hihi_bin_model_full,hihi_bin_model_int_6) #better = hihi_bin_model_full (not significant)

hihi_bin_models <- list(hihi_bin_model_null,hihi_bin_model_shannon,hihi_bin_model_pielou,hihi_bin_model_faith,hihi_bin_model_observed,hihi_bin_model_Staphylococcus,hihi_bin_model_Enterococcus,hihi_bin_model_Pseudomonas,hihi_bin_model_Enterobacteriaceae,hihi_bin_model_clutch.size,hihi_bin_model_dam.age,hihi_bin_model_dam.age.sq,hihi_bin_model_sire.age,hihi_bin_model_sire.age.sq,hihi_bin_model_mean.temp,hihi_bin_model_mean.rain,hihi_bin_model_full, hihi_bin_model_full.sq, hihi_bin_model_int_1, hihi_bin_model_int_2, hihi_bin_model_int_3,hihi_bin_model_int_4,hihi_bin_model_int_5,hihi_bin_model_int_6)
hihi_bin_model_names <- c("hihi_bin_model_null","hihi_bin_model_shannon","hihi_bin_model_pielou","hihi_bin_model_faith","hihi_bin_model_observed","hihi_bin_model_Staphylococcus","hihi_bin_model_Enterococcus","hihi_bin_model_Pseudomonas","hihi_bin_model_Enterobacteriaceae","hihi_bin_model_clutch.size","hihi_bin_model_dam.age","hihi_bin_model_dam.age.sq","hihi_bin_model_sire.age","hihi_bin_model_sire.age.sq",'hihi_bin_model_mean.temp',"hihi_bin_model_mean.rain","hihi_bin_model_full","hihi_bin_model_full.sq", "hihi_bin_model_int_1", "hihi_bin_model_int_2", "hihi_bin_model_int_3","hihi_bin_model_int_4","hihi_bin_model_int_5","hihi_bin_model_int_6")

hihi_bin_models_ranked <- aictab(cand.set = hihi_bin_models, modnames = hihi_bin_model_names, sort = TRUE, second.ord=T)
hihi_bin_models_ranked[1,]
summary(hihi_bin_model_full.sq)
confint(hihi_bin_model_full.sq)
hihi_bin_models_best_Cum.Wt <- hihi_bin_models_ranked[which(hihi_bin_models_ranked$Cum.Wt<=0.95),]
hihi_bin_models_best_AICc <- hihi_bin_models_ranked[which(hihi_bin_models_ranked$Delta_AICc<=2),]

write.csv(hihi_bin_models_ranked,"hihi_bin_models_ranked.csv")

hihi_bin_mod.sel <- model.sel(hihi_bin_models)
hihi_bin_mod.sel_delta2 <- subset(hihi_bin_mod.sel, delta <2)
hihi_bin_mod.avg <- model.avg(hihi_bin_mod.sel,subset=delta <2)
summary(hihi_bin_mod.avg)
sw(hihi_bin_mod.avg) # relative importance
confint(hihi_bin_mod.avg,full=TRUE) # full-model averaged 95% confidence intervals

# Hihi proportional clutch hatching success
hihi_prop_model_null <- glmer(clutch.hatch.prop ~ 1 + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))

hihi_prop_model_shannon <- glmer(clutch.hatch.prop ~ scale(shannon_entropy) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_shannon)
hihi_prop_model_pielou <- glmer(clutch.hatch.prop ~ scale(pielou_evenness) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_pielou)
hihi_prop_model_faith <- glmer(clutch.hatch.prop ~ scale(faith_pd) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_faith)
hihi_prop_model_observed <- glmer(clutch.hatch.prop ~ scale(observed_features) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_observed)
hihi_prop_model_Staphylococcus <- glmer(clutch.hatch.prop ~ scale(rel_abund_Staphylococcus) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_Staphylococcus) 
hihi_prop_model_Enterococcus <- glmer(clutch.hatch.prop ~ scale(rel_abund_Enterococcus) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_Enterococcus) 
hihi_prop_model_Pseudomonas <- glmer(clutch.hatch.prop ~ scale(rel_abund_Pseudomonas) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_Pseudomonas) 
hihi_prop_model_Enterobacteriaceae <- glmer(clutch.hatch.prop ~ scale(rel_abund_Enterobacteriaceae) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_Enterobacteriaceae) 
hihi_prop_model_clutch.size <- glmer(clutch.hatch.prop ~ scale(clutch_size) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_clutch.size)
hihi_prop_model_dam.age <- glmer(clutch.hatch.prop ~ scale(female_age_cohort) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_dam.age)
hihi_prop_model_dam.age.sq <- glmer(clutch.hatch.prop ~ scale(female_age_cohort.sq) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_dam.age.sq)
hihi_prop_model_sire.age <- glmer(clutch.hatch.prop ~ scale(male_age_cohort) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_sire.age)
hihi_prop_model_sire.age.sq <- glmer(clutch.hatch.prop ~ scale(male_age_cohort.sq) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_sire.age.sq)
hihi_prop_model_mean.temp <- glmer(clutch.hatch.prop ~ scale(mean_temp_laytoswab) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_mean.temp)
hihi_prop_model_mean.rain <- glmer(clutch.hatch.prop ~ scale(mean_rain_laytoswab) + (1|female_combo),data=hihi_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(hihi_prop_model_mean.rain)

anova(hihi_prop_model_null,hihi_prop_model_shannon) # better = hihi_prop_model_null (not significant)
anova(hihi_prop_model_null,hihi_prop_model_pielou) # better = hihi_prop_model_null (not significant)
anova(hihi_prop_model_null,hihi_prop_model_faith) # better = hihi_prop_model_null (not significant)
anova(hihi_prop_model_null,hihi_prop_model_observed) # better = hihi_prop_model_null (not significant)
anova(hihi_prop_model_null,hihi_prop_model_Staphylococcus) # better = hihi_prop_model_null (not significant)
anova(hihi_prop_model_null,hihi_prop_model_Enterococcus) # better = hihi_prop_model_null (not significant)
anova(hihi_prop_model_null,hihi_prop_model_Pseudomonas) # better = hihi_prop_model_null (not significant)
anova(hihi_prop_model_null,hihi_prop_model_Enterobacteriaceae) # better = hihi_prop_model_null (not significant)
anova(hihi_prop_model_null,hihi_prop_model_clutch.size) # better = hihi_prop_model_null (not significant) --> near sig.
anova(hihi_prop_model_null,hihi_prop_model_dam.age) # better = hihi_prop_model_dam (not significant) --> very close to sig.
anova(hihi_prop_model_null,hihi_prop_model_dam.age.sq) # better = hihi_prop_model_dam.age.sq (not significant) --> near sig.
anova(hihi_prop_model_null,hihi_prop_model_sire.age) # better = hihi_prop_model_null (not significant)
anova(hihi_prop_model_null,hihi_prop_model_sire.age.sq) # better = hihi_prop_model_null (not significant)
anova(hihi_prop_model_null,hihi_prop_model_mean.temp) # better = hihi_prop_model_null (not significant)
anova(hihi_prop_model_null,hihi_prop_model_mean.rain) # better = hihi_prop_model_null (not significant)

hihi_prop_models <- list(hihi_prop_model_null,hihi_prop_model_shannon,hihi_prop_model_pielou,hihi_prop_model_faith,hihi_prop_model_observed,hihi_prop_model_Staphylococcus,hihi_prop_model_Enterococcus,hihi_prop_model_Pseudomonas,hihi_prop_model_Enterobacteriaceae,hihi_prop_model_clutch.size,hihi_prop_model_dam.age,hihi_prop_model_dam.age.sq,hihi_prop_model_sire.age,hihi_prop_model_sire.age.sq,hihi_prop_model_mean.temp,hihi_prop_model_mean.rain)
hihi_prop_model_names <- c("hihi_prop_model_null","hihi_prop_model_shannon","hihi_prop_model_pielou","hihi_prop_model_faith","hihi_prop_model_observed","hihi_prop_model_Staphylococcus","hihi_prop_model_Enterococcus","hihi_prop_model_Pseudomonas","hihi_prop_model_Enterobacteriaceae","hihi_prop_model_clutch.size","hihi_prop_model_dam.age","hihi_prop_model_dam.age.sq","hihi_prop_model_sire.age","hihi_prop_model_sire.age.sq",'hihi_prop_model_mean.temp',"hihi_prop_model_mean.rain")

hihi_prop_models_ranked <- aictab(cand.set = hihi_prop_models, modnames = hihi_prop_model_names, sort = TRUE, second.ord=T)
hihi_prop_models_ranked[1,]
summary(hihi_prop_model_dam.age)
confint(hihi_prop_model_dam.age)
hihi_prop_models_best_Cum.Wt <- hihi_prop_models_ranked[which(hihi_prop_models_ranked$Cum.Wt<=0.95),]
hihi_prop_models_best_AICc <- hihi_prop_models_ranked[which(hihi_prop_models_ranked$Delta_AICc<=2),]
hihi_prop_models_best_AICc_list <- list(hihi_prop_models_best_AICc$Modnames)

write.csv(hihi_prop_models_ranked,"hihi_prop_models_ranked.csv")

hihi_prop_mod.sel <- model.sel(hihi_prop_models)
hihi_prop_mod.sel_delta2 <- subset(hihi_prop_mod.sel, delta <2)
hihi_prop_mod.avg <- model.avg(hihi_prop_mod.sel,subset=delta <2)
summary(hihi_prop_mod.avg)
sw(hihi_prop_mod.avg) # relative importance
confint(hihi_prop_mod.avg,full=TRUE) # full-model averaged 95% confidence intervals

## Kakī
kaki_data_bin_NArm <- na.omit(kaki_data[,c("clutch.nest_id_corrected","hatched","female_combo","shannon_entropy","pielou_evenness","faith_pd","observed_features","rel_abund_Staphylococcus","rel_abund_Enterococcus","rel_abund_Pseudomonas","rel_abund_Enterobacteriaceae","clutch_size","female_age_cohort","female_age_cohort.sq","clutch_number","collection_egg_weight","swab_day_parentalInc","wild.captive","female_lay_wild.captive","male_lay_wild.captive","mean_temp_laytoswab","mean_rain_laytoswab")])
kaki_data_prop_NArm <- na.omit(kaki_data[,c("clutch.nest_id_corrected","clutch.hatch.prop","female_combo","shannon_entropy","pielou_evenness","faith_pd","observed_features","rel_abund_Staphylococcus","rel_abund_Enterococcus","rel_abund_Pseudomonas","rel_abund_Enterobacteriaceae","clutch_size","female_age_cohort","female_age_cohort.sq","clutch_number","collection_egg_weight","swab_day_parentalInc","wild.captive","female_lay_wild.captive","male_lay_wild.captive","clutch_size","mean_temp_laytoswab","mean_rain_laytoswab")])

# Kakī binary hatching success
kaki_bin_model_null <- glmer(hatched ~ 1 + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))

kaki_bin_model_shannon <- glmer(hatched ~ scale(shannon_entropy) + (1|female_combo), data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_shannon) # Warning: Failed to converge
kaki_bin_model_pielou <- glmer(hatched ~ scale(pielou_evenness) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_pielou) # Warning: Failed to converge
kaki_bin_model_faith <- glmer(hatched ~ scale(faith_pd) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_faith)
kaki_bin_model_observed <- glmer(hatched ~ scale(observed_features) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_observed)
kaki_bin_model_Staphylococcus <- glmer(hatched ~ scale(rel_abund_Staphylococcus) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit")) # Warning: Failed to converge
summary(kaki_bin_model_Staphylococcus) 
kaki_bin_model_Enterococcus <- glmer(hatched ~ scale(rel_abund_Enterococcus) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_Enterococcus) # Warning: Model is nearly unidentifiable: large eigenvalue ratio
kaki_bin_model_Pseudomonas <- glmer(hatched ~ scale(rel_abund_Pseudomonas) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_Pseudomonas) # Warning: Failed to converge
kaki_bin_model_Enterobacteriaceae <- glmer(hatched ~ scale(rel_abund_Enterobacteriaceae) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_Enterobacteriaceae) 
kaki_bin_model_clutch.size <- glmer(hatched ~ scale(clutch_size) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_clutch.size)
kaki_bin_model_dam.age <- glmer(hatched ~ scale(female_age_cohort) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_dam.age)
kaki_bin_model_dam.age.sq <- glmer(hatched ~ scale(female_age_cohort.sq) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_dam.age.sq)
kaki_bin_model_clutch.no <- glmer(hatched ~ scale(clutch_number) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_clutch.no) # Warning: boundary (singular) fit
kaki_bin_model_egg.weight <- glmer(hatched ~ scale(collection_egg_weight) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_egg.weight) # Warning: Failed to converge
kaki_bin_model_incub <- glmer(hatched ~ scale(swab_day_parentalInc) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_incub) # Warning: boundary (singular) fit
kaki_bin_model_egg.wild.captive <- glmer(hatched ~ wild.captive + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_egg.wild.captive) # Warning: boundary (singular) fit
#kaki_bin_model_dam.wild.captive <- glmer(hatched ~ female_lay_wild.captive + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
#summary(kaki_bin_model_dam.wild.captive) # Warning: boundary (singular) fit
#kaki_bin_model_sire.wild.captive <- glmer(hatched ~ male_lay_wild.captive + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
#summary(kaki_bin_model_sire.wild.captive)
kaki_bin_model_mean.temp <- glmer(hatched ~ scale(mean_temp_laytoswab) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_mean.temp)
kaki_bin_model_mean.rain <- glmer(hatched ~ scale(mean_rain_laytoswab) + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit"))
summary(kaki_bin_model_mean.rain)

anova(kaki_bin_model_null,kaki_bin_model_shannon) # better = kaki_bin_model_null (not significant)
anova(kaki_bin_model_null,kaki_bin_model_pielou) # better = kaki_bin_model_null (not significant)
anova(kaki_bin_model_null,kaki_bin_model_faith) # better = kaki_bin_model_faith (not significant)
anova(kaki_bin_model_null,kaki_bin_model_observed) # better = kaki_bin_model_observed (not significant) --> very near sig.
anova(kaki_bin_model_null,kaki_bin_model_Staphylococcus) # better = kaki_bin_model_Staphylococcus (not significant)
anova(kaki_bin_model_null,kaki_bin_model_Enterococcus) # better = kaki_bin_model_Enterococcus (not significant) --> near sig.
anova(kaki_bin_model_null,kaki_bin_model_Pseudomonas) # better = kaki_bin_model_null (not significant)
anova(kaki_bin_model_null,kaki_bin_model_Enterobacteriaceae) # better = kaki_bin_model_null (not significant)
anova(kaki_bin_model_null,kaki_bin_model_clutch.size) # better = kaki_bin_model_null (not significant)
anova(kaki_bin_model_null,kaki_bin_model_dam.age) # better = kaki_bin_model_null (not significant)
anova(kaki_bin_model_null,kaki_bin_model_dam.age.sq) # better = kaki_bin_model_null (not significant)
anova(kaki_bin_model_null,kaki_bin_model_clutch.no) # better = kaki_bin_model_clutch.no (not significant)
anova(kaki_bin_model_null,kaki_bin_model_egg.weight) # better = kaki_bin_model_null (not significant)
anova(kaki_bin_model_null,kaki_bin_model_incub) # better = kaki_bin_model_incub (significant) *
anova(kaki_bin_model_null,kaki_bin_model_egg.wild.captive) # better = kaki_bin_model_egg.wild.captive (significant) *
#anova(kaki_bin_model_null,kaki_bin_model_dam.wild.captive) # better = kaki_bin_model_dam.wild.captive (significant) *
#anova(kaki_bin_model_null,kaki_bin_model_sire.wild.captive) # better = kaki_bin_model_null (not significant)
anova(kaki_bin_model_null,kaki_bin_model_mean.temp) # better = kaki_bin_model_null (not significant)
anova(kaki_bin_model_null,kaki_bin_model_mean.rain) # better = kaki_bin_model_null (not significant)

kaki_bin_model_full <- glmer(hatched ~ scale(swab_day_parentalInc) + wild.captive + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit")) #Warning: boundary (singular) fit
check_overdispersion(kaki_bin_model_full) # No overdispersion detected

anova(kaki_bin_model_full, kaki_bin_model_incub) # better = kaki_bin_model_incub (not significant)
anova(kaki_bin_model_full, kaki_bin_model_egg.wild.captive) # better = kaki_bin_model_full (significant)

kaki_bin_model_int_1 <- glmer(hatched ~ scale(swab_day_parentalInc) * wild.captive + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit")) #Warning: boundary (singular) fit
check_overdispersion(kaki_bin_model_int_1) # No overdispersion detected
kaki_bin_model_int_2 <- glmer(hatched ~ scale(swab_day_parentalInc) * rel_abund_Pseudomonas + (1|female_combo),data=kaki_data_bin_NArm,family=binomial(link="logit")) #Warning: boundary (singular) fit
check_overdispersion(kaki_bin_model_int_2) # No overdispersion detected

kaki_bin_models <- list(kaki_bin_model_null,kaki_bin_model_shannon,kaki_bin_model_pielou,kaki_bin_model_faith,kaki_bin_model_observed,kaki_bin_model_Staphylococcus,kaki_bin_model_Enterococcus,kaki_bin_model_Pseudomonas,kaki_bin_model_Enterobacteriaceae,kaki_bin_model_clutch.size,kaki_bin_model_dam.age,kaki_bin_model_dam.age.sq,kaki_bin_model_clutch.no,kaki_bin_model_egg.weight,kaki_bin_model_incub,kaki_bin_model_egg.wild.captive,kaki_bin_model_mean.temp,kaki_bin_model_mean.rain,kaki_bin_model_full, kaki_bin_model_int_1, kaki_bin_model_int_2)
kaki_bin_model_names <- c("kaki_bin_model_null","kaki_bin_model_shannon","kaki_bin_model_pielou","kaki_bin_model_faith","kaki_bin_model_observed","kaki_bin_model_Staphylococcus","kaki_bin_model_Enterococcus","kaki_bin_model_Pseudomonas","kaki_bin_model_Enterobacteriaceae","kaki_bin_model_clutch.size","kaki_bin_model_dam.age","kaki_bin_model_dam.age.sq","kaki_bin_model_clutch.no","kaki_bin_model_egg.weight","kaki_bin_model_incub","kaki_bin_model_egg.wild.captive","kaki_bin_model_mean.temp","kaki_bin_model_mean.rain","kaki_bin_model_full", "kaki_bin_model_int_1", "kaki_bin_model_int_2")

kaki_bin_models_ranked <- aictab(cand.set = kaki_bin_models, modnames = kaki_bin_model_names, sort = TRUE, second.ord=T)
kaki_bin_models_ranked
kaki_bin_models_ranked[1,]
summary(kaki_bin_model_incub)
confint(kaki_bin_model_incub)
kaki_bin_models_best_Cum.Wt <- kaki_bin_models_ranked[which(kaki_bin_models_ranked$Cum.Wt<=0.95),]
kaki_bin_models_best_AICc <- kaki_bin_models_ranked[which(kaki_bin_models_ranked$Delta_AICc<=2),]
kaki_bin_models_best_AICc

write.csv(kaki_bin_models_ranked,"kaki_bin_models_ranked.csv")

kaki_bin_mod.sel <- model.sel(kaki_bin_models)
kaki_bin_mod.sel_delta2 <- subset(kaki_bin_mod.sel, delta <2)
kaki_bin_mod.avg <- model.avg(kaki_bin_mod.sel,subset=delta <2)
summary(kaki_bin_mod.avg)
sw(kaki_bin_mod.avg) # relative importance
confint(kaki_bin_mod.avg,full=TRUE) # full-model averaged 95% confidence intervals

# Kakī proportional clutch hatching success
kaki_prop_model_null <- glmer(clutch.hatch.prop ~ 1 + (1|female_combo), data=kaki_data_prop_NArm, weights=clutch_size, family=binomial(link="logit")) # Warning: Failed to converge

kaki_prop_model_shannon <- glmer(clutch.hatch.prop ~ scale(shannon_entropy) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_shannon)
kaki_prop_model_pielou <- glmer(clutch.hatch.prop ~ scale(pielou_evenness) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_pielou)
kaki_prop_model_faith <- glmer(clutch.hatch.prop ~ scale(faith_pd) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_faith)
kaki_prop_model_observed <- glmer(clutch.hatch.prop ~ scale(observed_features) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_observed)
kaki_prop_model_Staphylococcus <- glmer(clutch.hatch.prop ~ scale(rel_abund_Staphylococcus) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_Staphylococcus) 
kaki_prop_model_Enterococcus <- glmer(clutch.hatch.prop ~ scale(rel_abund_Enterococcus) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_Enterococcus) 
kaki_prop_model_Pseudomonas <- glmer(clutch.hatch.prop ~ scale(rel_abund_Pseudomonas) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_Pseudomonas) 
kaki_prop_model_Enterobacteriaceae <- glmer(clutch.hatch.prop ~ scale(rel_abund_Enterobacteriaceae) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_Enterobacteriaceae) # Warning: Failed to converge
kaki_prop_model_clutch.size <- glmer(clutch.hatch.prop ~ scale(clutch_size) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_clutch.size) # Warning: Failed to converge
kaki_prop_model_dam.age <- glmer(clutch.hatch.prop ~ scale(female_age_cohort) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_dam.age)
kaki_prop_model_dam.age.sq <- glmer(clutch.hatch.prop ~ scale(female_age_cohort.sq) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_dam.age.sq)
kaki_prop_model_clutch.no <- glmer(clutch.hatch.prop ~ scale(clutch_number) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_clutch.no)
kaki_prop_model_egg.weight <- glmer(clutch.hatch.prop ~ scale(collection_egg_weight) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_egg.weight)
kaki_prop_model_incub <- glmer(clutch.hatch.prop ~ scale(swab_day_parentalInc) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_incub)
kaki_prop_model_egg.wild.captive <- glmer(clutch.hatch.prop ~ wild.captive + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_egg.wild.captive)
#kaki_prop_model_dam.wild.captive <- glmer(clutch.hatch.prop ~ female_lay_wild.captive + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
#summary(kaki_prop_model_dam.wild.captive)
#kaki_prop_model_sire.wild.captive <- glmer(clutch.hatch.prop ~ male_lay_wild.captive + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
#summary(kaki_prop_model_sire.wild.captive)
kaki_prop_model_mean.temp <- glmer(clutch.hatch.prop ~ scale(mean_temp_laytoswab) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_mean.temp)
kaki_prop_model_mean.rain <- glmer(clutch.hatch.prop ~ scale(mean_rain_laytoswab) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_mean.rain)

anova(kaki_prop_model_null,kaki_prop_model_shannon) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_pielou) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_faith) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_observed) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_Staphylococcus) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_Enterococcus) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_Pseudomonas) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_Enterobacteriaceae) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_clutch.size) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_dam.age) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_dam.age.sq) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_clutch.no) # better = kaki_prop_model_clutch.no (significant) *
anova(kaki_prop_model_null,kaki_prop_model_egg.weight) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_incub) # better = kaki_prop_model_incub (not significant)
anova(kaki_prop_model_null,kaki_prop_model_egg.wild.captive) # better = kaki_prop_model_egg.wild.captive (not significant) --> near sig.
#anova(kaki_prop_model_null,kaki_prop_model_dam.wild.captive) # better = kaki_prop_model_dam.wild.captive (significant) *
#anova(kaki_prop_model_null,kaki_prop_model_sire.wild.captive) # better = kaki_prop_model_null (not significant)
anova(kaki_prop_model_null,kaki_prop_model_mean.temp) # better = kaki_prop_model_mean.temp (significant) *
anova(kaki_prop_model_null,kaki_prop_model_mean.rain) # better = kaki_prop_model_null (not significant)

kaki_prop_model_full <- glmer(clutch.hatch.prop ~ scale(clutch_number) + scale(mean_temp_laytoswab) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
summary(kaki_prop_model_full)
check_overdispersion(kaki_prop_model_full) # No overdispersion detected

anova(kaki_prop_model_full, kaki_prop_model_clutch.no) # better = kaki_prop_model_full (not significant)
anova(kaki_prop_model_full, kaki_prop_model_mean.temp) # better = kaki_prop_model_full (not significant)

kaki_prop_model_int_1 <- glmer(clutch.hatch.prop ~ scale(clutch_number) * scale(mean_temp_laytoswab) + (1|female_combo),data=kaki_data_prop_NArm,weights=clutch_size,family=binomial(link="logit"))
check_overdispersion(kaki_prop_model_int_1) # No overdispersion detected

anova(kaki_prop_model_full,kaki_prop_model_int_1) # better = kaki_prop_model_full (not significant)

kaki_prop_models <- list(kaki_prop_model_null,kaki_prop_model_shannon,kaki_prop_model_pielou,kaki_prop_model_faith,kaki_prop_model_observed,kaki_prop_model_Staphylococcus,kaki_prop_model_Enterococcus,kaki_prop_model_Pseudomonas,kaki_prop_model_Enterobacteriaceae,kaki_prop_model_clutch.size,kaki_prop_model_dam.age,kaki_prop_model_dam.age.sq,kaki_prop_model_clutch.no,kaki_prop_model_egg.weight,kaki_prop_model_incub,kaki_prop_model_egg.wild.captive,kaki_prop_model_mean.temp,kaki_prop_model_mean.rain,kaki_prop_model_full, kaki_prop_model_int_1)
kaki_prop_model_names <- c("kaki_prop_model_null","kaki_prop_model_shannon","kaki_prop_model_pielou","kaki_prop_model_faith","kaki_prop_model_observed","kaki_prop_model_Staphylococcus","kaki_prop_model_Enterococcus","kaki_prop_model_Pseudomonas","kaki_prop_model_Enterobacteriaceae","kaki_prop_model_clutch.size","kaki_prop_model_dam.age","kaki_prop_model_dam.age.sq","kaki_prop_model_clutch.no","kaki_prop_model_egg.weight","kaki_prop_model_incub","kaki_prop_model_egg.wild.captive","kaki_prop_model_mean.temp","kaki_prop_model_mean.rain","kaki_prop_model_full", "kaki_prop_model_int_1")

kaki_prop_models_ranked <- aictab(cand.set = kaki_prop_models, modnames = kaki_prop_model_names, sort = TRUE, second.ord=T)
kaki_prop_models_ranked[1,]
summary(kaki_prop_model_full)
confint(kaki_prop_model_full)
kaki_prop_models_best_Cum.Wt <- kaki_prop_models_ranked[which(kaki_prop_models_ranked$Cum.Wt<=0.95),]
kaki_prop_models_best_AICc <- kaki_prop_models_ranked[which(kaki_prop_models_ranked$Delta_AICc<=2),]
kaki_prop_models_best_AICc_list <- list(kaki_prop_models_best_AICc$Modnames)

write.csv(kaki_prop_models_ranked,"kaki_prop_models_ranked.csv")

kaki_prop_mod.sel <- model.sel(kaki_prop_models)
kaki_prop_mod.sel_delta2 <- subset(kaki_prop_mod.sel, delta <2)
kaki_prop_mod.avg <- model.avg(kaki_prop_mod.sel,subset=delta <2)
summary(kaki_prop_mod.avg)
sw(kaki_prop_mod.avg) # relative importance
confint(kaki_prop_mod.avg,full=TRUE) # full-model averaged 95% confidence intervals







