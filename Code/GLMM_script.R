##Clear workspace
rm(list=ls())

#### INSTALL AND LOAD PACKAGES #### DO NOT UPDATE ANY OF THESE PACKAGES!!!
#library("dada2"); packageVersion("dada2") #‘1.24.0’
library("phyloseq"); packageVersion("phyloseq") #‘1.40.0’
#library("vegan"); packageVersion("vegan") #‘2.6.4’
#library("Biostrings"); packageVersion("Biostrings") #‘2.64.1’
#library("ggplot2"); packageVersion("ggplot2") #‘3.4.0’
#library("decontam"); packageVersion("decontam") #'1.16.0'
library("qiime2R"); packageVersion("qiime2R") #‘0.99.6’
library("tidyverse"); packageVersion("tidyverse") #'1.3.2'
#library("biomformat"); packageVersion("biomformat") #'1.24.0'
#library("DESeq2"); packageVersion("DESeq2") #‘1.36.0'
#library("dendextend"); packageVersion("dendextend") #‘1.16.0’
library("ggview"); packageVersion("ggview") #‘0.1.0’
library("lme4"); packageVersion("lme4") #‘1.1.31’
library("MuMIn"); packageVersion("MuMIn") #‘1.47.1’
library("stringr"); packageVersion("stringr") #‘1.4.1’
library("car"); packageVersion("car") #'3.1.1'
library("performance"); packageVersion("performance") #'0.10.1'

setwd("Data/DataAnalyses")
mtd <- read_tsv(("sample-metadata-new-noNA-removemock_forR_wEnvironmentalData_2.tsv"), comment = "#q2")

#Import feature table and convert to count table
table_250.256_silva_138.1.filtered_decontam_noMCnoNTC <- read_qza("filtered-silva-138.1-table-decontam_0.4_noMCnoNTC.qza")
count_tab_250.256_silva_138.1.filtered_decontam_noMCnoNTC <- table_250.256_silva_138.1.filtered_decontam_noMCnoNTC$data %>% as.data.frame() 

#Import taxonomies
taxonomy_silva_138.1_250.256_filtered_decontam <- read_qza("silva_138.1_taxonomy_250.256-filtered-decontam.qza")

#Convert format to work with decontam
tax_tab_silva_138.1_250.256_filtered_decontam <- taxonomy_silva_138.1_250.256_filtered_decontam$data %>% 
  as.data.frame() %>%
  mutate(Taxon = gsub("D_0", "k", Taxon), Taxon = gsub("D_1", "p", Taxon),
         Taxon = gsub("D_2", "c", Taxon), Taxon = gsub("D_3", "o", Taxon),
         Taxon = gsub("D_4", "f", Taxon), Taxon = gsub("D_5", "g", Taxon),
         Taxon = gsub("D_6", "s", Taxon)) %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  select(-Confidence) #Warning message: Expected 7 pieces. Missing pieces filled with `NA` in 4252 rows [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, ...] = this matches the number of NAs in the species column

physeq_feature.tree.tax <- qza_to_phyloseq(features="filtered-silva-138.1-table-decontam_0.4_noMCnoNTC.qza",
                               tree="rooted-tree_allfilter_decontam_noMCnoNTC.qza",
                               taxonomy="silva_138.1_taxonomy_250.256-filtered-decontam.qza")
physeq_meta <- phyloseq(sample_data(column_to_rownames(mtd, "SampleID")))
physeq_tree_noMCnoNTC <- merge_phyloseq(physeq_feature.tree.tax,physeq_meta)
physeq_tree_noMCnoNTC #6528 taxa and 249 samples
sample_names(physeq_tree_noMCnoNTC)

##Following van Veelen (2018)

#Subset samples into kaki and hihi
physeq_tree_noMCnoNTC_hihi <- subset_samples(physeq_tree_noMCnoNTC, !species == "Kaki")
physeq_tree_noMCnoNTC_hihi #6528 taxa and 157 samples
physeq_tree_noMCnoNTC_kaki <- subset_samples(physeq_tree_noMCnoNTC, !species == "Hihi")
physeq_tree_noMCnoNTC_kaki #6528 taxa and 92 samples

#Rarefy data based on QIIME2 alpha rarefaction curves to keep all samples
physeq_tree_noMCnoNTC_rarefied = rarefy_even_depth(physeq_tree_noMCnoNTC,rngseed=1,sample.size=4896,replace=F)
physeq_tree_noMCnoNTC_rarefied #5839 taxa and 249 samples
physeq_tree_noMCnoNTC_hihi_rarefied = rarefy_even_depth(physeq_tree_noMCnoNTC_hihi,rngseed=1,sample.size=4896,replace=F)
physeq_tree_noMCnoNTC_hihi_rarefied #3055 taxa and 157 samples
physeq_tree_noMCnoNTC_kaki_rarefied = rarefy_even_depth(physeq_tree_noMCnoNTC_kaki,rngseed=1,sample.size=14396 ,replace=F)
physeq_tree_noMCnoNTC_kaki_rarefied #3922 taxa and 92 samples
sample_data(physeq_tree_noMCnoNTC_kaki_rarefied)

#### Alpha diversity metrics ####

# Create an estimate_richness file for each subset
# Plot a suite of Alpha diversity richness metrics for full ASV table
sample_data(physeq_tree_noMCnoNTC_rarefied)
# Estimate richness eggshells of rarefied data
richness_all.eggs_rarefied <- estimate_richness(physeq_tree_noMCnoNTC_rarefied) #adds X to front of kaki rownames and changes all - to .
richness_hihi.eggs_rarefied <- estimate_richness(physeq_tree_noMCnoNTC_hihi_rarefied) #changes all - to .
richness_kaki.eggs_rarefied <- estimate_richness(physeq_tree_noMCnoNTC_kaki_rarefied) #adds X to front of kaki rownames and changes all - to .

# Fix rownames: cut the "X" and leave the second part of the split
row.names(richness_kaki.eggs_rarefied) <- sapply(strsplit(rownames(richness_kaki.eggs_rarefied),"X"),"[",2)

# Can't work out how to get rid of X in all eggshells but doesn't matter because need to deal with hihi and kakī seperately from here because used different rarefaction depths

## Prepare data

mtd.df <- as.data.frame(mtd) #convert metadata to dataframe
mtd.df <- mtd.df[order(mtd.df$SampleID),] #re-order the metadata to match the richness dataframe
mtd.df_hihi <- subset(mtd.df, mtd.df$species == "Hihi") #subset metadata dataframe to just hihi samples
mtd.df_kaki <- subset(mtd.df, mtd.df$species == "Kaki") #subset metadata dataframe to just kaki samples

richness_kaki.eggs_rarefied_df <- data.frame(names=row.names(richness_kaki.eggs_rarefied),richness_kaki.eggs_rarefied)
rownames(richness_kaki.eggs_rarefied_df) <- NULL #match format to metadata dataframe
dim(richness_kaki.eggs_rarefied_df) #92 10
richness_hihi.eggs_rarefied_df <- data.frame(names=row.names(richness_hihi.eggs_rarefied),richness_hihi.eggs_rarefied)
rownames(richness_hihi.eggs_rarefied_df) <- NULL #match format to metadata dataframe
dim(richness_hihi.eggs_rarefied_df) #157 10

# Select just the metadata columns needed for each dataset
mtd.df_hihi_data <- mtd.df_hihi[,c("SampleID","plate_ID","clutch.nest_id_corrected","clutch_size","hatched_total","female_cohort","female_combo","female_age_cohort","female_experienced.naive_social","female_first.season_social","female_total.partners_social","female_total.partners_genetic","male_cohort","Male_combo","male_age_cohort","male_experienced.naive_social","male_first.season_social","male_total.partners_social","male_total.partners_genetic","1st_egg_lay","1st_egg_lay_mintemp","1st_egg_lay_maxtemp","1st_egg_lay_rain","hatched","hatched.unhatched","clutch.hatch%","swab_date")]
dim(mtd.df_hihi_data) #157 34

mtd.df_kaki_data <- mtd.df_kaki[,c("SampleID","plate_ID","wild.captive","clutch.nest_id_corrected","nest_location","clutch_number","clutch_size","hatched_total","female_cohort","female_combo","female_age_cohort","female_lay_wild.captive","female_hatch_wild.captive","female_experienced.naive_social","female_first.season_social","female_total.partners_social","male_cohort","Male_combo","male_age_cohort","male_lay_wild.captive","male_hatch_wild.captive","male_experienced.naive_social","male_first.season_social","male_total.partners_social","1st_egg_lay","1st_egg_lay_mintemp","1st_egg_lay_maxtemp","1st_egg_lay_rain","hatched","hatched.unhatched","clutch.hatch%","swab_date")]
dim(mtd.df_kaki_data) #92 37

names(mtd.df_hihi_data)[names(mtd.df_hihi_data) == "clutch.hatch%"] <- "clutch.hatch.percent"
names(mtd.df_kaki_data)[names(mtd.df_kaki_data) == "clutch.hatch%"] <- "clutch.hatch.percent"

mtd.df_hihi_data$clutch.hatch.prop <- mtd.df_hihi_data$hatched_total/mtd.df_hihi_data$clutch_size
mtd.df_kaki_data$clutch.hatch.prop <- mtd.df_kaki_data$hatched_total/mtd.df_kaki_data$clutch_size

richness_hihi.eggs_rarefied_df$SampleID <- mtd.df_hihi_data$SampleID
richness_kaki.eggs_rarefied_df$SampleID <- mtd.df_kaki_data$SampleID

# Add calculated alpha diversity metrics as columns in the dataframes
final.data_hihi <- merge(mtd.df_hihi_data,richness_hihi.eggs_rarefied_df,by="SampleID")
final.data_kaki <- merge(mtd.df_kaki_data,richness_kaki.eggs_rarefied_df,by="SampleID")

# Add column for days into breeding season clutch was initiated = date 1st clutch was initiated = 0 (because Julian date won't work for kaki because it crosses 2019/20); kaki first clutch initiated = 11/10/2019, hihi first clutch initiated = 01/11/2019
final.data_hihi[,"FirstLay"] <- "2019/11/01"
final.data_hihi$`1st_egg_lay` <- format(as.Date(final.data_hihi$`1st_egg_lay`,format="%d/%m/%Y"),"%Y/%m/%d")
final.data_hihi[,"SeasonDay"] <- as.numeric(difftime(final.data_hihi$`1st_egg_lay`,final.data_hihi$FirstLay,units=c("days"),tz="UTC"))
final.data_hihi[,"SeasonDay"] 

final.data_kaki[,"FirstLay"] <- "2019/10/11"
final.data_kaki$`1st_egg_lay` <- format(as.Date(final.data_kaki$`1st_egg_lay`,format="%d/%m/%Y"),"%Y/%m/%d")
final.data_kaki[,"SeasonDay"] <- as.numeric(difftime(final.data_kaki$`1st_egg_lay`,final.data_kaki$FirstLay,units=c("days"),tz="UTC"))
final.data_kaki[,"SeasonDay"]

# Add column for days into season egg was swabbed = using date 1st clutch was initiated as start of season (because Julian date won't work for kaki because it crosses 2019/20)
final.data_hihi$swab_date <- format(as.Date(final.data_hihi$swab_date,format="%d/%m/%Y"),"%Y/%m/%d")
final.data_hihi[,"SwabDay"] <- as.numeric(difftime(final.data_hihi$swab_date,final.data_hihi$FirstLay,units=c("days"),tz="UTC"))
final.data_hihi[,"SwabDay"] 

final.data_kaki$swab_date <- format(as.Date(final.data_kaki$swab_date,format="%d/%m/%Y"),"%Y/%m/%d")
final.data_kaki[,"SwabDay"] <- as.numeric(difftime(final.data_kaki$swab_date,final.data_kaki$FirstLay,units=c("days"),tz="UTC"))
final.data_kaki[,"SwabDay"] 

# Save the final datasets
#write.csv(final.data_hihi,"final.data_hihi_forglmm.csv")
#write.csv(final.data_kaki,"final.data_kaki_forglmm.csv")

## Hihi
# Hatched vs. unhatched (binary)
hihi_bin_model_sat_df <- na.omit(final.data_hihi[,c("hatched","Shannon","female_experienced.naive_social","male_experienced.naive_social","clutch_size","female_age_cohort","female_total.partners_social","female_total.partners_genetic","male_age_cohort","male_total.partners_social","male_total.partners_genetic","1st_egg_lay_mintemp","1st_egg_lay_maxtemp","1st_egg_lay_rain","SeasonDay","clutch.nest_id_corrected","female_combo","plate_ID")])
hihi_bin_model_sat <- glmer(hatched~
                                    scale(Shannon) +
                                    female_experienced.naive_social +
                                    male_experienced.naive_social +
                                    scale(clutch_size) +
#                                    scale(female_age_cohort^2) +
                                   scale(female_total.partners_social) +
#                                   scale(female_total.partners_genetic) +
                                    scale(male_age_cohort) +
                                   scale(male_total.partners_social) +
#                                   scale(male_total.partners_genetic) +
                                    scale(`1st_egg_lay_mintemp`) +
                                    scale(`1st_egg_lay_maxtemp`) +
                                    scale(`1st_egg_lay_rain`) +
                                    scale(SeasonDay) +
                                    (1|female_combo) +
                                    (1|plate_ID),
                                  family=binomial(link="logit"),
                                  data=hihi_bin_model_sat_df,
                                  na.action="na.fail")
#Warnings: boundary (singular) fit: see help('isSingular')
summary(hihi_bin_model_sat)
car::vif(hihi_bin_model_sat) # Use VIF < 5 as threshold
check_overdispersion(hihi_bin_model_sat) # p > 0.05 = No overdispersion detected
hihi_bin_model_selection <- dredge(hihi_bin_model_sat,trace=2)
nrow(hihi_bin_model_selection) #2048
#saveRDS(hihi_bin_model_selection,"hihi_bin_model_selection.rds")
#readRDS("hihi_bin_model_selection.rds")
hihi_bin_top_model <- get.models(hihi_bin_model_selection,subset = 1)[[1]] #Warning: boundary (singular) fit: see help('isSingular')
hihi_bin_top_model
car::vif(hihi_bin_top_model)
summary(hihi_bin_top_model)
hihi_bin_averaged_models <-(model.avg(hihi_bin_model_selection,subset = delta <=2))
summary(hihi_bin_averaged_models)
sw(hihi_bin_averaged_models) # relative importance
confint(hihi_bin_averaged_models,full=TRUE) # full-model averaged 95% confidence intervals

# Clutch hatching success (proportion) - need to add clutch size as weights to use binomial family
hihi_prop_model_sat_df <- na.omit(final.data_hihi[,c("clutch.hatch.prop","Shannon","female_experienced.naive_social","male_experienced.naive_social","clutch_size","female_age_cohort","female_total.partners_social","female_total.partners_genetic","male_age_cohort","male_total.partners_social","male_total.partners_genetic","1st_egg_lay_mintemp","1st_egg_lay_maxtemp","1st_egg_lay_rain","SeasonDay","clutch.nest_id_corrected","female_combo","plate_ID")])
hihi_prop_model_sat <- glmer(clutch.hatch.prop~
                            scale(Shannon) +
                            female_experienced.naive_social +
                            male_experienced.naive_social +
                            scale(clutch_size) +
#                            scale(female_age_cohort^2) +
                            scale(female_total.partners_social) +
#                            scale(female_total.partners_genetic) +
                            scale(male_age_cohort) +
                            scale(male_total.partners_social) +
                            scale(male_total.partners_genetic) +
                            scale(`1st_egg_lay_mintemp`) +
                            scale(`1st_egg_lay_maxtemp`) +
                            scale(`1st_egg_lay_rain`) +
                            scale(SeasonDay) + 
                            (1|female_combo) +
                            (1|plate_ID),family=binomial(link="logit"),
                            weights=clutch_size,
                            data=hihi_prop_model_sat_df,
                            na.action="na.fail")
#Warnings: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :Model failed to converge with max|grad| = 0.0050977 (tol = 0.002, component 1)
summary(hihi_prop_model_sat)
car::vif(hihi_prop_model_sat) # Use VIF < 5 as threshold
check_overdispersion(hihi_prop_model_sat) # # p > 0.05 = No overdispersion detected
hihi_prop_model_selection <- dredge(hihi_prop_model_sat,trace=2)
nrow(hihi_prop_model_selection) #4096
#saveRDS(hihi_prop_model_selection,"hihi_prop_model_selection.rds")
#readRDS("hihi_prop_model_selection.rds")
hihi_prop_top_model <- get.models(hihi_prop_model_selection,subset = 1)[[1]] #Warning: boundary (singular) fit: see help('isSingular')
hihi_prop_top_model
car::vif(hihi_prop_top_model)
summary(hihi_prop_top_model)
hihi_prop_averaged_models <- (model.avg(hihi_prop_model_selection,subset = delta <=2))
summary(hihi_prop_averaged_models)
sw(hihi_prop_averaged_models) # relative importance
confint(hihi_prop_averaged_models,full=TRUE) # full-model averaged 95% confidence intervals

## Kaki
# Hatched vs. unhatched (binary)
kaki_bin_model_sat_df <- na.omit(final.data_kaki[,c("hatched","Shannon","wild.captive","nest_location","female_lay_wild.captive","male_lay_wild.captive","female_experienced.naive_social","male_experienced.naive_social","clutch_number","clutch_size","female_age_cohort","female_total.partners_social","male_age_cohort","male_total.partners_social","1st_egg_lay_mintemp","1st_egg_lay_maxtemp","1st_egg_lay_rain","SeasonDay","clutch.nest_id_corrected","female_combo","plate_ID")])
kaki_bin_model_sat <- glmer(hatched~
                              scale(Shannon) + 
                              wild.captive +
#                              nest_location +
                              female_lay_wild.captive +
                              male_lay_wild.captive +
#                              female_experienced.naive_social +
                              male_experienced.naive_social +
#                              scale(clutch_number) +
                              scale(clutch_size) +
#                              scale(female_age_cohort) +
                              scale(female_total.partners_social) +
#                              scale(male_age_cohort) +
                              scale(male_total.partners_social) +
                              scale(`1st_egg_lay_mintemp`) +
#                              scale(`1st_egg_lay_maxtemp`) +
                              scale(`1st_egg_lay_rain`) +
                              scale(SeasonDay) +
                              (1|female_combo) + 
                              (1|plate_ID),
                            family=binomial(link="logit"),
                            data=kaki_bin_model_sat_df,
                            na.action="na.fail")
#Warnings: boundary (singular) fit: see help('isSingular')
summary(kaki_bin_model_sat)
car::vif(kaki_bin_model_sat) # Use VIF < 5 as threshold 
check_overdispersion(kaki_bin_model_sat) # # p > 0.05 = No overdispersion detected
kaki_bin_model_selection <- dredge(kaki_bin_model_sat,trace=2)
nrow(kaki_bin_model_selection) #2048
#saveRDS(kaki_bin_model_selection,"kaki_bin_model_selection.rds")
#readRDS("kaki_bin_model_selection.rds")
kaki_bin_top_model <- get.models(kaki_bin_model_selection,subset = 1)[[1]] #Warning: boundary (singular) fit: see help('isSingular')
kaki_bin_top_model
car::vif(kaki_bin_top_model)
summary(kaki_bin_top_model)
kaki_bin_averaged_models <- model.avg(kaki_bin_model_selection,subset = delta <=2)
summary(kaki_bin_averaged_models)
sw(kaki_bin_averaged_models) # relative importance
confint(kaki_bin_averaged_models,full=TRUE) # full-model averaged 95% confidence intervals

# Clutch hatching success (proportion) - need to add clutch size as weights to use binomial family
kaki_prop_model_sat_df <- na.omit(final.data_kaki[,c("clutch.hatch.prop","Shannon","wild.captive","nest_location","female_lay_wild.captive","male_lay_wild.captive","female_experienced.naive_social","male_experienced.naive_social","clutch_number","clutch_size","female_age_cohort","female_total.partners_social","male_age_cohort","male_total.partners_social","1st_egg_lay_mintemp","1st_egg_lay_maxtemp","1st_egg_lay_rain","SeasonDay","clutch.nest_id_corrected","female_combo","plate_ID")])
kaki_prop_model_sat <- glmer(clutch.hatch.prop~
                               scale(Shannon) + 
                               wild.captive +
#                               nest_location +
                               female_lay_wild.captive +
                               male_lay_wild.captive +
#                               female_experienced.naive_social +
                               male_experienced.naive_social +
#                               scale(clutch_number) +
                               scale(clutch_size) +
#                               scale(female_age_cohort) +
                               scale(female_total.partners_social) +
#                               scale(male_age_cohort) +
                               scale(male_total.partners_social) +
                               scale(`1st_egg_lay_mintemp`) +
#                               scale(`1st_egg_lay_maxtemp`) +
                               scale(`1st_egg_lay_rain`) +
                               scale(SeasonDay) +
                               (1|female_combo) + 
                               (1|plate_ID),
                             family=binomial(link="logit"),
                             weights=clutch_size,
                             data=kaki_prop_model_sat_df,
                             na.action="na.fail")
#Warnings: boundary (singular) fit: see help('isSingular')
summary(kaki_prop_model_sat)
car::vif(kaki_prop_model_sat) # Use VIF < 5 as threshold 
check_overdispersion(kaki_prop_model_sat) # p > 0.05 = No overdispersion detected
kaki_prop_model_selection <- dredge(kaki_prop_model_sat,trace=2)
nrow(kaki_prop_model_selection) #2048
#saveRDS(kaki_prop_model_selection,"kaki_prop_model_selection.rds")
#readRDS("kaki_prop_model_selection.rds")
kaki_prop_top_model <- get.models(kaki_prop_model_selection,subset = 1)[[1]] #Warning: boundary (singular) fit: see help('isSingular')
kaki_prop_top_model
car::vif(kaki_prop_top_model)
summary(kaki_prop_top_model)
kaki_prop_averaged_models <- (model.avg(kaki_prop_model_selection,subset = delta <=2))
summary(kaki_prop_averaged_models)
sw(kaki_prop_averaged_models) # relative importance
confint(kaki_prop_averaged_models,full=TRUE) # full-model averaged 95% confidence intervals

options(na.action = "na.omit")

# Null model = just the random factor(s)
hihi_model_null_1 <- glmer(hatched ~ 1 + (1|female_combo),data=final.data_hihi,family=binomial)
hihi_model_null_2 <- glmer(hatched ~ 1 + (1|female_combo) + (1|plate_ID),data=final.data_hihi,family=binomial)
anova(hihi_model_null_1,hihi_model_null_2) # hihi_model_null_1 better (lower AIC)
hihi_model_full <- glmer(hatched ~ 
                           Shannon + 
                           (female_age_cohort^2) + 
                           clutch_size +
                           (1|female_combo),
                         data=final.data_hihi,family=binomial)
anova(hihi_model_null_1,hihi_model_full) # hihi_model_full model better (lower AIC)
summary(hihi_model_full) # clutch size not significant
hihi_model_1 <- glmer(hatched ~ 
                           Shannon + 
                           (female_age_cohort^2) + 
#                           clutch_size +
                           (1|female_combo),
                         data=final.data_hihi,family=binomial) 
anova(hihi_model_full,hihi_model_1) # hihi_model_1 better (lower AIC)
hihi_model_2 <- glmer(hatched ~ 
                        Shannon + 
                        female_age_cohort + 
                        #clutch_size +
                        (1|female_combo),
                      data=final.data_hihi,family=binomial) 
anova(hihi_model_1,hihi_model_2) # both models exactly the same??




#####

#### ALPHA AND BETA DIVERSITY ####
##https://micca.readthedocs.io/en/latest/phyloseq.html

sample_data(physeq_tree_noMCnoNTC)

tab <- (otu_table(physeq_tree_noMCnoNTC))
class(tab) <- "matrix"
rarecurve(t(tab),step=50,cex=0.5)

physeq_tree_noMCnoNTC_rarefied = rarefy_even_depth(physeq_tree_noMCnoNTC,rngseed=1,sample.size=4896,replace=F) 
#`set.seed(1)` was used to initialize repeatable random subsampling
#689OTUs were removed because they are no longer present in any sample after random subsampling

plot_bar(physeq_tree_noMCnoNTC_rarefied,fill="Phylum")
plot_bar(physeq_tree_noMCnoNTC_rarefied, fill="Phylum") + facet_wrap(~species, scales="free_x", nrow=1)
physeq_tree_noMCnoNTC_rarefied_phylum_NArmFalse <- tax_glom(physeq_tree_noMCnoNTC_rarefied,"Phylum",NArm=FALSE)
physeq_tree_noMCnoNTC_rarefied_phylum_NArmFalse
plot_bar(physeq_tree_noMCnoNTC_rarefied_phylum_NArmFalse, fill="Phylum") + facet_wrap(~species, scales="free_x", nrow=1)
ggview(units = "px", height = 1800, width = 10000)

plot_richness(physeq_tree_noMCnoNTC_rarefied, x="species", color="species", measures=c("Observed"))
plot_richness(physeq_tree_noMCnoNTC_rarefied, x="species", color="species", measures=c("Shannon"))
plot_richness(physeq_tree_noMCnoNTC_rarefied, x="species", color="species", measures=c("Simpson"))
plot_richness(physeq_tree_noMCnoNTC_rarefied, x="species", measures=c("Observed", "Shannon","Simpson")) + geom_boxplot()
rich = estimate_richness(physeq_tree_noMCnoNTC_rarefied)
rich
pairwise.wilcox.test(rich$Observed, sample_data(physeq_tree_noMCnoNTC_rarefied)$species)
pairwise.wilcox.test(rich$Shannon, sample_data(physeq_tree_noMCnoNTC_rarefied)$species)
pairwise.wilcox.test(rich$Simpson, sample_data(physeq_tree_noMCnoNTC_rarefied)$species)

# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(physeq_tree_noMCnoNTC_rarefied, method="unifrac", weighted=F)
ordination = ordinate(physeq_tree_noMCnoNTC_rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(physeq_tree_noMCnoNTC_rarefied, ordination, color="species") + theme(aspect.ratio=1)

adonis2(wunifrac_dist ~ sample_data(physeq_tree_noMCnoNTC_rarefied)$species)




