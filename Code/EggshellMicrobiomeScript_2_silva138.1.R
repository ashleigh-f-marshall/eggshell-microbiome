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
library("DECIPHER"); packageVersion("DECIPHER") #‘2.24.0’
library("tidyverse"); packageVersion("tidyverse") #‘1.3.2’
library("DESeq2"); packageVersion("DESeq2") #‘1.36.0'
library("dendextend"); packageVersion("dendextend") #‘1.16.0’
library("viridis"); packageVersion("viridis") #‘0.6.2’
library("ANCOMBC"); packageVersion("ANCOMBC") #'1.6.4'

#### LOAD AND PREPARE DATA ####
setwd("Data/210212_Data")
path <- "Data/210212_Data"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

#saveRDS(fnFs, file = "outputs/fnFs.RDS") 
#saveRDS(fnRs, file = "outputs/fnRs.RDS") 
fnFs <- readRDS("outputs/fnFs.RDS")
fnRs <- readRDS("outputs/fnRs.RDS")
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#### VISUALISE DATA ####

# Quality profiles of forward reads
plotQualityProfile(fnFs[1:2]) #just shows first 2

# Quality profiles of reverse reads
plotQualityProfile(fnRs[1:2]) #just shows first 2

#### FILTER AND TRIM
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#saveRDS(filtFs, file = "outputs/filtFs.RDS") 
#saveRDS(filtRs, file = "outputs/filtRs.RDS")
filtFs <- readRDS("outputs/filtFs.RDS")
filtRs <- readRDS("outputs/filtRs.RDS")

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])

#saveRDS(out, file = "outputs/out.RDS")
out <- readRDS("outputs/out.RDS")

# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE) #100488960 total bases in 418704 reads from 15 samples will be used for learning the error rates
errR <- learnErrors(filtRs, multithread=TRUE) #100055200 total bases in 625345 reads from 20 samples will be used for learning the error rates.
plotErrors(errF, nominalQ=TRUE) #Looks OK
plotErrors(errR, nominalQ=TRUE) #Looks OK

#saveRDS(errF, file = "outputs/errF.RDS")
#saveRDS(errR, file = "outputs/errR.RDS")
errF <- readRDS("outputs/errF.RDS")
errR <- readRDS("outputs/errR.RDS")

##From https://astrobiomike.github.io/amplicon/dada2_workflow_ex
#derep_forward <- derepFastq(filtFs, verbose=TRUE)
#names(derep_forward) <- sample.names # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
#derep_reverse <- derepFastq(filtRs, verbose=TRUE)
#names(derep_reverse) <- sample.names
#dada_forward <- dada(derep_forward, err=errF, pool="pseudo")
# dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder
#dada_reverse <- dada(derep_reverse, err=errR, pool="pseudo")
# dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE) #255 samples - default pool = FALSE (pooling can increase sensitivity to sequence variants that may be present at very low frequencies in multiple samples, but independent inference is still very accurate and is less prone to reporting types of false-positives like contaminants)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE) #255 samples - default pool = FALSE
dadaFs[[1]] ##58 sequence variants were inferred from 6575 input unique sequences
dadaRs[[1]] ##48 sequence variants were inferred from 4965 input unique sequences

#saveRDS(dadaFs, file = "outputs/dadaFs.RDS") 
#saveRDS(dadaRs, file = "outputs/dadaRs.RDS") 
dadaFs <- readRDS("outputs/dadaFs.RDS")
dadaRs <- readRDS("outputs/dadaRs.RDS")

#### MERGE PAIRED READS ####
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE) #by default merged sequences only output if forward and reverse overlap by at least 12 bases and are identical in overlapping region
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#saveRDS(mergers, file = "outputs/mergers.RDS") 
mergers <- readRDS("outputs/mergers.RDS")

#### CONSTRUCT SEQUENCE TABLE ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab) #255 samples, 8507 ASVs
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab))) #lengths between 240 and 388, mainly clustered around 253 but ~10% at 240 (trimmed length of forward reads?)
reads.by.len <- tapply(colSums(seqtab), nchar(colnames(seqtab)), sum)
plot(as.numeric(names(reads.by.len)), reads.by.len, xlab="Length", ylab="Total Reads") #Looks ok? Highest point around 253, only small bump at 240
# Remove non-target-length sequences
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]
dim(seqtab2) #255 samples, 7532 ASVs
table(nchar(getSequences(seqtab2)))

#saveRDS(seqtab, file = "outputs/seqtab.RDS") 
#saveRDS(seqtab2, file = "outputs/seqtab2.RDS") 
seqtab <- readRDS("outputs/seqtab.RDS")
seqtab2 <- readRDS("outputs/seqtab2.RDS")

#### REMOVE CHIMERAS ####
#seqtab.nochim_full <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) #Identified 718 bimeras out of 8507 input sequences = 8.4%
#dim(seqtab.nochim_full) #255 samples, 7789 ASVs
#table(nchar(getSequences(seqtab.nochim_full)))
#sum(seqtab.nochim_full)/sum(seqtab) #0.9932448 - chimeras only account for around 0.7% of the merged sequence reads

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE) #Identified 501 bimeras out of 7532 input sequences = 6.7%
dim(seqtab.nochim) #255 samples, 7031 ASVs
table(nchar(getSequences(seqtab.nochim)))
sum(seqtab.nochim)/sum(seqtab2) #0.9970991 - chimeras only account for around 0.3% of the merged sequence reads

#saveRDS(seqtab.nochim, file = "outputs/seqtab.nochim.RDS")
seqtab.nochim <- readRDS("outputs/seqtab.nochim.RDS")

#### TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track_full <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim_full),finalpercent=round(rowSums(seqtab.nochim_full)/out[,1]*100,1))
colnames(track_full) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","final%")
rownames(track_full) <- sample.names
head(track_full) #Looks good
track_full[c("MC-1","MC-2","MC-3"),]
mean(round(rowSums(seqtab.nochim_full)/out[,1]*100,1)) #88.03%
#write.csv(track_full,"Tracked_Reads_full.csv")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim),finalpercent=round(rowSums(seqtab.nochim)/out[,1]*100,1))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","final%")
rownames(track) <- sample.names
head(track) #Looks good
track[c("MC-1","MC-2","MC-3"),]
mean(round(rowSums(seqtab.nochim)/out[,1]*100,1)) #86.47%
#write.csv(track,"Tracked_Reads.csv")

#Evaluate accuracy of ASV recognition
unqs.mock_1 <- seqtab.nochim["MC-1",]
unqs.mock_1 <- sort(unqs.mock_1[unqs.mock_1>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock_1), "sample sequences present in the Mock community.\n") #DADA2 inferred 22 sample sequences present in the Mock community
unqs.mock_2 <- seqtab.nochim["MC-2",]
unqs.mock_2 <- sort(unqs.mock_2[unqs.mock_2>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock_2), "sample sequences present in the Mock community.\n") #DADA2 inferred 22 sample sequences present in the Mock community
unqs.mock_3 <- seqtab.nochim["MC-3",]
unqs.mock_3 <- sort(unqs.mock_3[unqs.mock_3>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock_3), "sample sequences present in the Mock community.\n") #DADA2 inferred 25 sample sequences present in the Mock community

unqs.mock_1_full <- seqtab.nochim_full["MC-1",]
unqs.mock_1_full <- sort(unqs.mock_1_full[unqs.mock_1_full>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock_1_full), "sample sequences present in the Mock community.\n") #DADA2 inferred 23 sample sequences present in the Mock community
unqs.mock_2_full <- seqtab.nochim_full["MC-2",]
unqs.mock_2_full <- sort(unqs.mock_2_full[unqs.mock_2_full>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock_2_full), "sample sequences present in the Mock community.\n") #DADA2 inferred 23 sample sequences present in the Mock community
unqs.mock_3_full <- seqtab.nochim_full["MC-3",]
unqs.mock_3_full <- sort(unqs.mock_3_full[unqs.mock_3_full>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock_3_full), "sample sequences present in the Mock community.\n") #DADA2 inferred 25 sample sequences present in the Mock community

##No reference sequences available for mock community but should be 25 ASVs post analysis

#### ASSIGN TAXONOMY #####
taxa_silva_sp <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)

taxa_silva_sp.print <- taxa_silva_sp # Removing sequence rownames for display only
rownames(taxa_silva_sp.print) <- NULL
head(taxa_silva_sp.print)

#saveRDS(taxa_silva_sp, file = "outputs/taxa_silva_sp.RDS")
taxa_silva_sp <- readRDS("outputs/taxa_silva_sp.RDS")
dim(taxa_silva_sp) #7031 7

#write.csv(taxa_silva_sp, "taxa_silva_sp.csv")

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("SILVA_SSU_r138_2019.RData")
ids_top <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ids_both <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE) # use all processors
#saveRDS(ids_top, file = "outputs/ids_top.RDS")
#saveRDS(ids_both, file = "outputs/ids_both.RDS")
ids_top <- readRDS("outputs/ids_top.RDS")
ids_both <- readRDS("outputs/ids_both.RDS")
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid_top <- t(sapply(ids_top, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid_top) <- ranks; rownames(taxid_top) <- getSequences(seqtab.nochim)

taxid_both <- t(sapply(ids_both, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid_both) <- ranks; rownames(taxid_both) <- getSequences(seqtab.nochim)

#write.csv(taxid_top, "taxa_decipher_top.csv")
#write.csv(taxid_both, "taxa_decipher_both.csv")

#saveRDS(taxid, file = "outputs/taxid.RDS")
taxid <- readRDS("outputs/taxid.RDS")

#Remove non-target sequences (chloroplasts, mitochondria, and archaea) N.B. eukaryota already removed when reference databases reformatted for DADA2 compatibility
dim(seqtab.nochim) #255 7031
dim(taxa_silva_sp) #7031 6

is.chloro_silva_sp <- taxa_silva_sp[,"Order"] %in% "Chloroplast"
table(is.chloro_silva_sp) #241 chloroplasts
seqtab.nochloro_silva_sp <- seqtab.nochim[,!is.chloro_silva_sp]
dim(seqtab.nochloro_silva_sp) #255 6790
taxa_silva_sp.nochloro <- taxa_silva_sp[!is.chloro_silva_sp,]
dim(taxa_silva_sp.nochloro) #6781 7

is.mito_silva_sp <- taxa_silva_sp.nochloro[,"Family"] %in% "Mitochondria"
table(is.mito_silva_sp) #129 mitochondria
seqtab.nochloro.nomito_silva_sp <- seqtab.nochloro_silva_sp[,!is.mito_silva_sp]
dim(seqtab.nochloro.nomito_silva_sp) #255 6661
taxa_silva_sp.nochloro.nomito <- taxa_silva_sp.nochloro[!is.mito_silva_sp,]
dim(taxa_silva_sp.nochloro.nomito) #6661 7

is.archaea_silva_sp <- taxa_silva_sp.nochloro.nomito[,"Kingdom"] %in% "Archaea"
table(is.archaea_silva_sp) #10 archaea
seqtab.nochloro.nomito.noarchaea_silva_sp <- seqtab.nochloro.nomito_silva_sp[,!is.archaea_silva_sp]
dim(seqtab.nochloro.nomito.noarchaea_silva_sp) #255 6651
taxa_silva_sp.nochloro.nomito.noarchaea <- taxa_silva_sp.nochloro.nomito[!is.archaea_silva_sp,]
dim(taxa_silva_sp.nochloro.nomito.noarchaea) #6651 7

is.phylaNA_silva_sp <- is.na(taxa_silva_sp.nochloro.nomito.noarchaea[,"Phylum"])
table(is.phylaNA_silva_sp) #26 NA
seqtab.nochloro.nomito.noarchaea.nophylaNA_silva_sp <- seqtab.nochloro.nomito.noarchaea_silva_sp[,!is.phylaNA_silva_sp]
dim(seqtab.nochloro.nomito.noarchaea.nophylaNA_silva_sp) #255 6625
taxa_silva_sp.nochloro.nomito.noarchaea.nophylaNA <- taxa_silva_sp.nochloro.nomito.noarchaea[!is.phylaNA_silva_sp,]
dim(taxa_silva_sp.nochloro.nomito.noarchaea.nophylaNA) #6625 7

seqtab_clean_silva_sp <- seqtab.nochloro.nomito.noarchaea.nophylaNA_silva_sp #255 6625
taxa_clean_silva_sp <- taxa_silva_sp.nochloro.nomito.noarchaea.nophylaNA #6625 7

#### EXTRACT PRODUCTS OF DADA2 ####

#Give headers easier names
asv_seqs_silva_sp <- colnames(seqtab_clean_silva_sp)
asv_headers_silva_sp <- vector(dim(seqtab_clean_silva_sp)[2], mode="character")

for (i in 1:dim(seqtab_clean_silva_sp)[2]) {
  asv_headers_silva_sp[i] <- paste(">ASV", i, sep="_")
}

#Make and write fasta of final ASV sequences
asv_fasta_silva_sp <- c(rbind(asv_headers_silva_sp, asv_seqs_silva_sp))
#write(asv_fasta_silva_sp, "ASVs_silva_sp.fa")

#Make and write count table
asv_tab_silva_sp <- t(seqtab_clean_silva_sp)
row.names(asv_tab_silva_sp) <- sub(">", "", asv_headers_silva_sp)
#write.table(asv_tab_silva_sp, "ASVs_counts_silva_sp.tsv", sep="\t", quote=F, col.names=NA)
#write.csv(asv_tab_silva_sp, "ASVs_counts_silva_sp.csv")

#Make and write taxonomy table
colnames(taxa_clean_silva_sp) #"Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
rownames(taxa_clean_silva_sp)
asv_tax_silva_sp <- taxa_clean_silva_sp
rownames(asv_tax_silva_sp) <- gsub(pattern=">", replacement="", x=asv_headers_silva_sp)
colnames(asv_tax_silva_sp)
#write.table(asv_tax_silva_sp, "ASVs_taxonomy_silva_sp.tsv", sep = "\t", quote=F, col.names=NA)
#write.csv(asv_tax_silva_sp, "ASVs_taxonomy_silva_sp.csv")

#Remove likely contaminants based on prevalence in control samples
colnames(asv_tab_silva_sp)
negative_controls_silva_sp <- asv_tab_silva_sp[,c(237:239)]
vector_for_decontam_silva_sp <- c(rep(FALSE,236),rep(TRUE,3),rep(FALSE,16))
contam_df_silva_sp <- isContaminant(t(asv_tab_silva_sp), neg=vector_for_decontam_silva_sp, method="prevalence")
table(contam_df_silva_sp$contaminant) #Identified 6 contaminants
contam_df_silva_sp_0.5 <- isContaminant(t(asv_tab_silva_sp), neg=vector_for_decontam_silva_sp, method="prevalence",threshold=0.5)
table(contam_df_silva_sp_0.5$contaminant) #Identified 18 contaminants
contam_df_silva_sp_0.2 <- isContaminant(t(asv_tab_silva_sp), neg=vector_for_decontam_silva_sp, method="prevalence",threshold=0.2)
table(contam_df_silva_sp_0.2$contaminant) #Identified 9 contaminants
contam_df_silva_sp_0.3 <- isContaminant(t(asv_tab_silva_sp), neg=vector_for_decontam_silva_sp, method="prevalence",threshold=0.3)
table(contam_df_silva_sp_0.3$contaminant) #Identified 11 contaminants
contam_df_silva_sp_0.4 <- isContaminant(t(asv_tab_silva_sp), neg=vector_for_decontam_silva_sp, method="prevalence",threshold=0.4)
table(contam_df_silva_sp_0.4$contaminant) #Identified 15 contaminants
contam_asvs_silva_sp <- row.names(contam_df_silva_sp[contam_df_silva_sp$contaminant == TRUE, ])
silva_sp_contaminants <-asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp, ]
contam_asvs_silva_sp_0.5 <- row.names(contam_df_silva_sp_0.5[contam_df_silva_sp_0.5$contaminant == TRUE, ])
silva_sp_contaminants_0.5 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_0.5, ]
contam_asvs_silva_sp_0.2 <- row.names(contam_df_silva_sp_0.2[contam_df_silva_sp_0.2$contaminant == TRUE, ])
silva_sp_contaminants_0.2 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_0.2, ]
contam_asvs_silva_sp_0.3 <- row.names(contam_df_silva_sp_0.3[contam_df_silva_sp_0.3$contaminant == TRUE, ])
silva_sp_contaminants_0.3 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_0.3, ]
contam_asvs_silva_sp_0.4 <- row.names(contam_df_silva_sp_0.4[contam_df_silva_sp_0.4$contaminant == TRUE, ])
silva_sp_contaminants_0.4 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_0.4, ]
#write.csv(silva_sp_contaminants,"silva_sp_contaminants.csv")
#write.csv(silva_sp_contaminants_0.5,"silva_sp_contaminants_0.5.csv")
#write.csv(silva_sp_contaminants_0.2,"silva_sp_contaminants_0.2.csv")
#write.csv(silva_sp_contaminants_0.3,"silva_sp_contaminants_0.3.csv")
#write.csv(silva_sp_contaminants_0.4,"silva_sp_contaminants_0.4.csv")

# Checking for differences if remove mock communities before running decontam
colnames(asv_tab_silva_sp)
mock_communities_silva_sp <- asv_tab_silva_sp[,c(234:236)]
asv_tab_silva_sp_nomock <- asv_tab_silva_sp[,c(1:233,237:255)]
colnames(asv_tab_silva_sp_nomock)
negative_controls_silva_sp_nomock <- asv_tab_silva_sp_nomock[,c(234:236)]
vector_for_decontam_silva_sp_nomock <- c(rep(FALSE,233),rep(TRUE,3),rep(FALSE,16))
contam_df_silva_sp_nomock <- isContaminant(t(asv_tab_silva_sp_nomock), neg=vector_for_decontam_silva_sp_nomock, method="prevalence")
table(contam_df_silva_sp_nomock$contaminant) #Identified 5 contaminants
contam_df_silva_sp_nomock_0.5 <- isContaminant(t(asv_tab_silva_sp_nomock), neg=vector_for_decontam_silva_sp_nomock, method="prevalence",threshold=0.5)
table(contam_df_silva_sp_nomock_0.5$contaminant) #Identified 18 contaminants
contam_df_silva_sp_nomock_0.2 <- isContaminant(t(asv_tab_silva_sp_nomock), neg=vector_for_decontam_silva_sp_nomock, method="prevalence",threshold=0.2)
table(contam_df_silva_sp_nomock__0.2$contaminant) #Identified 9 contaminants
contam_df_silva_sp_nomock_0.3 <- isContaminant(t(asv_tab_silva_sp_nomock), neg=vector_for_decontam_silva_sp_nomock, method="prevalence",threshold=0.3)
table(contam_df_silva_sp_nomock_0.3$contaminant) #Identified 11 contaminants
contam_df_silva_sp_nomock_0.4 <- isContaminant(t(asv_tab_silva_sp_nomock), neg=vector_for_decontam_silva_sp_nomock, method="prevalence",threshold=0.4)
table(contam_df_silva_sp_nomock_0.4$contaminant) #Identified 15 contaminants
contam_asvs_silva_sp_nomock <- row.names(contam_df_silva_sp_nomock[contam_df_silva_sp_nomock$contaminant == TRUE, ])
silva_sp_nomock_contaminants <-asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_nomock, ] #Identical to with mock but didn't pick up ASV_48
contam_asvs_silva_sp_nomock_0.5 <- row.names(contam_df_silva_sp_nomock_0.5[contam_df_silva_sp_nomock_0.5$contaminant == TRUE, ])
silva_sp_nomock_contaminants_0.5 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_nomock_0.5, ] #Identical to with mock
contam_asvs_silva_sp_nomock_0.2 <- row.names(contam_df_silva_sp_nomock_0.2[contam_df_silva_sp_nomock_0.2$contaminant == TRUE, ])
silva_sp_nomock_contaminants_0.2 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_nomock_0.2, ] #Identical to with mock
contam_asvs_silva_sp_nomock_0.3 <- row.names(contam_df_silva_sp_nomock_0.3[contam_df_silva_sp_nomock_0.3$contaminant == TRUE, ])
silva_sp_nomock_contaminants_0.3 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_nomock_0.3, ] #Identical to with mock
contam_asvs_silva_sp_nomock_0.4 <- row.names(contam_df_silva_sp_nomock_0.4[contam_df_silva_sp_nomock_0.4$contaminant == TRUE, ])
silva_sp_nomock_contaminants_0.4 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_nomock_0.4, ] #Identical to with mock

hist(contam_df_silva_sp_nomock$p, 100,)

metadata_all <- read.table("Data/EggshellMicrobiomeProject.csv",header=TRUE,sep=",",row.names=1)
ps <- phyloseq(otu_table(seqtab_clean_silva_sp,taxa_are_rows=FALSE),
               sample_data(metadata_all),
               tax_table(taxa_clean_silva_sp))
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contam_df_silva_sp_0.4$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

## Running decontam on kakī and hihi samples separately (as recommended in Davis et al. 2018)
colnames(asv_tab_silva_sp)
asv_tab_silva_sp_kaki_NTC <- asv_tab_silva_sp[,c(1:92,237:239)] #Select just kakī and NTC, leave out mock
colnames(asv_tab_silva_sp_kaki_NTC)
asv_tab_silva_sp_hihi_NTC <- asv_tab_silva_sp[,c(93:233,237:255)] #Select just hihi and NTC, leave out mock
colnames(asv_tab_silva_sp_hihi_NTC)

colnames(asv_tab_silva_sp_kaki_NTC)
negative_controls_silva_sp_kaki_NTC <- asv_tab_silva_sp_kaki_NTC[,c(93:95)]
vector_for_decontam_silva_sp_kaki_NTC <- c(rep(FALSE,92),rep(TRUE,3))
contam_df_silva_sp_kaki_NTC <- isContaminant(t(asv_tab_silva_sp_kaki_NTC), neg=vector_for_decontam_silva_sp_kaki_NTC, method="prevalence")
table(contam_df_silva_sp_kaki_NTC$contaminant) #Identified 6 contaminants
contam_df_silva_sp_kaki_NTC_0.5 <- isContaminant(t(asv_tab_silva_sp_kaki_NTC), neg=vector_for_decontam_silva_sp_kaki_NTC, method="prevalence",threshold=0.5)
table(contam_df_silva_sp_kaki_NTC_0.5$contaminant) #Identified 18 contaminants
contam_df_silva_sp_kaki_NTC_0.2 <- isContaminant(t(asv_tab_silva_sp_kaki_NTC), neg=vector_for_decontam_silva_sp_kaki_NTC, method="prevalence",threshold=0.2)
table(contam_df_silva_sp_kaki_NTC_0.2$contaminant) #Identified 13 contaminants
contam_df_silva_sp_kaki_NTC_0.3 <- isContaminant(t(asv_tab_silva_sp_kaki_NTC), neg=vector_for_decontam_silva_sp_kaki_NTC, method="prevalence",threshold=0.3)
table(contam_df_silva_sp_kaki_NTC_0.3$contaminant) #Identified 14 contaminants
contam_df_silva_sp_kaki_NTC_0.4 <- isContaminant(t(asv_tab_silva_sp_kaki_NTC), neg=vector_for_decontam_silva_sp_kaki_NTC, method="prevalence",threshold=0.4)
table(contam_df_silva_sp_kaki_NTC_0.4$contaminant) #Identified 15 contaminants
contam_asvs_silva_sp_kaki_NTC <- row.names(contam_df_silva_sp_kaki_NTC[contam_df_silva_sp_kaki_NTC$contaminant == TRUE, ])
silva_sp_kaki_NTC_contaminants <-asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_kaki_NTC, ]
contam_asvs_silva_sp_kaki_NTC_0.5 <- row.names(contam_df_silva_sp_kaki_NTC_0.5[contam_df_silva_sp_kaki_NTC_0.5$contaminant == TRUE, ])
silva_sp_kaki_NTC_contaminants_0.5 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_kaki_NTC_0.5, ]
contam_asvs_silva_sp_kaki_NTC_0.2 <- row.names(contam_df_silva_sp_kaki_NTC_0.2[contam_df_silva_sp_kaki_NTC_0.2$contaminant == TRUE, ])
silva_sp_kaki_NTC_contaminants_0.2 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_kaki_NTC_0.2, ]
contam_asvs_silva_sp_kaki_NTC_0.3 <- row.names(contam_df_silva_sp_kaki_NTC_0.3[contam_df_silva_sp_kaki_NTC_0.3$contaminant == TRUE, ])
silva_sp_kaki_NTC_contaminants_0.3 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_kaki_NTC_0.3, ]
contam_asvs_silva_sp_kaki_NTC_0.4 <- row.names(contam_df_silva_sp_kaki_NTC_0.4[contam_df_silva_sp_kaki_NTC_0.4$contaminant == TRUE, ])
silva_sp_kaki_NTC_contaminants_0.4 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_kaki_NTC_0.4, ]
#write.csv(silva_sp_kaki_NTC_contaminants,"silva_sp_kaki_NTC_contaminants.csv")
#write.csv(silva_sp_kaki_NTC_contaminants_0.5,"silva_sp_kaki_NTC_contaminants_0.5.csv")
#write.csv(silva_sp_kaki_NTC_contaminants_0.2,"silva_sp_kaki_NTC_contaminants_0.2.csv")
#write.csv(silva_sp_kaki_NTC_contaminants_0.3,"silva_sp_kaki_NTC_contaminants_0.3.csv")
#write.csv(silva_sp_kaki_NTC_contaminants_0.4,"silva_sp_kaki_NTC_contaminants_0.4.csv")

hist(contam_df_silva_sp_kaki_NTC$p, 100,)

colnames(asv_tab_silva_sp_hihi_NTC)
negative_controls_silva_sp_hihi_NTC <- asv_tab_silva_sp_hihi_NTC[,c(142:144)]
vector_for_decontam_silva_sp_hihi_NTC <- c(rep(FALSE,141),rep(TRUE,3),rep(FALSE,16))
contam_df_silva_sp_hihi_NTC <- isContaminant(t(asv_tab_silva_sp_hihi_NTC), neg=vector_for_decontam_silva_sp_hihi_NTC, method="prevalence")
table(contam_df_silva_sp_hihi_NTC$contaminant) #Identified 4 contaminants
contam_df_silva_sp_hihi_NTC_0.5 <- isContaminant(t(asv_tab_silva_sp_hihi_NTC), neg=vector_for_decontam_silva_sp_hihi_NTC, method="prevalence",threshold=0.5)
table(contam_df_silva_sp_hihi_NTC_0.5$contaminant) #Identified 16 contaminants
contam_df_silva_sp_hihi_NTC_0.2 <- isContaminant(t(asv_tab_silva_sp_hihi_NTC), neg=vector_for_decontam_silva_sp_hihi_NTC, method="prevalence",threshold=0.2)
table(contam_df_silva_sp_hihi_NTC_0.2$contaminant) #Identified 8 contaminants
contam_df_silva_sp_hihi_NTC_0.3 <- isContaminant(t(asv_tab_silva_sp_hihi_NTC), neg=vector_for_decontam_silva_sp_hihi_NTC, method="prevalence",threshold=0.3)
table(contam_df_silva_sp_hihi_NTC_0.3$contaminant) #Identified 9 contaminants
contam_df_silva_sp_hihi_NTC_0.4 <- isContaminant(t(asv_tab_silva_sp_hihi_NTC), neg=vector_for_decontam_silva_sp_hihi_NTC, method="prevalence",threshold=0.4)
table(contam_df_silva_sp_hihi_NTC_0.4$contaminant) #Identified 12 contaminants
contam_asvs_silva_sp_hihi_NTC <- row.names(contam_df_silva_sp_hihi_NTC[contam_df_silva_sp_hihi_NTC$contaminant == TRUE, ])
silva_sp_hihi_NTC_contaminants <-asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_hihi_NTC, ]
contam_asvs_silva_sp_hihi_NTC_0.5 <- row.names(contam_df_silva_sp_hihi_NTC_0.5[contam_df_silva_sp_hihi_NTC_0.5$contaminant == TRUE, ])
silva_sp_hihi_NTC_contaminants_0.5 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_hihi_NTC_0.5, ]
contam_asvs_silva_sp_hihi_NTC_0.2 <- row.names(contam_df_silva_sp_hihi_NTC_0.2[contam_df_silva_sp_hihi_NTC_0.2$contaminant == TRUE, ])
silva_sp_hihi_NTC_contaminants_0.2 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_hihi_NTC_0.2, ]
contam_asvs_silva_sp_hihi_NTC_0.3 <- row.names(contam_df_silva_sp_hihi_NTC_0.3[contam_df_silva_sp_hihi_NTC_0.3$contaminant == TRUE, ])
silva_sp_hihi_NTC_contaminants_0.3 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_hihi_NTC_0.3, ]
contam_asvs_silva_sp_hihi_NTC_0.4 <- row.names(contam_df_silva_sp_hihi_NTC_0.4[contam_df_silva_sp_hihi_NTC_0.4$contaminant == TRUE, ])
silva_sp_hihi_NTC_contaminants_0.4 <- asv_tax_silva_sp[row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_hihi_NTC_0.4, ]
#write.csv(silva_sp_hihi_NTC_contaminants,"silva_sp_hihi_NTC_contaminants.csv")
#write.csv(silva_sp_hihi_NTC_contaminants_0.5,"silva_sp_hihi_NTC_contaminants_0.5.csv")
#write.csv(silva_sp_hihi_NTC_contaminants_0.2,"silva_sp_hihi_NTC_contaminants_0.2.csv")
#write.csv(silva_sp_hihi_NTC_contaminants_0.3,"silva_sp_hihi_NTC_contaminants_0.3.csv")
#write.csv(silva_sp_hihi_NTC_contaminants_0.4,"silva_sp_hihi_NTC_contaminants_0.4.csv")

hist(contam_df_silva_sp_hihi_NTC$p, 100,)

## Going with 0.4 cut-off as majority seem likely to be contaminants based on biological likelihood and previous identification as contaminants in other studies, plus doesn't exclude any ASVs not identified as contaminants in each species separately

# Making new fasta file
contam_indices_silva_sp_0.4 <- which(asv_fasta_silva_sp %in% paste0(">", contam_asvs_silva_sp_0.4))
dont_want_silva_sp_0.4 <- sort(c(contam_indices_silva_sp_0.4, contam_indices_silva_sp_0.4 + 1))
asv_fasta_no_contam_silva_sp_0.4 <- asv_fasta_silva_sp[- dont_want_silva_sp_0.4]

# Making new count table
asv_tab_no_contam_silva_sp_0.4 <- asv_tab_silva_sp[!row.names(asv_tab_silva_sp) %in% contam_asvs_silva_sp_0.4, ] #6610 255

# Making new taxonomy table
asv_tax_no_contam_silva_sp_0.4 <- asv_tax_silva_sp[!row.names(asv_tax_silva_sp) %in% contam_asvs_silva_sp_0.4, ] #6610 7

## Writing them out to files
#write(asv_fasta_no_contam_silva_sp_0.4, "ASVs-no-contam_silva_sp_0.4.fa")
#write.table(asv_tab_no_contam_silva_sp_0.4, "ASVs_counts-no-contam_silva_sp_0.4.tsv", sep="\t", quote=F, col.names=NA)
#write.table(asv_tax_no_contam_silva_sp_0.4, "ASVs_taxonomy-no-contam_silva_sp_0.4.tsv", sep="\t", quote=F, col.names=NA)


#### STATISTICAL ANALYSES #####

rm(list=ls())

#### INSTALL AND LOAD PACKAGES #### DO NOT UPDATE ANY OF THESE PACKAGES!!!
library("dada2"); packageVersion("dada2") #‘1.24.0’
library("phyloseq"); packageVersion("phyloseq") #‘1.40.0’
library("vegan"); packageVersion("vegan") #‘2.6.4’
library("Biostrings"); packageVersion("Biostrings") #‘2.64.1’
library("ggplot2"); packageVersion("ggplot2") #‘3.4.0’
library("decontam"); packageVersion("decontam") #'1.16.0'
library("DECIPHER"); packageVersion("DECIPHER") #‘2.24.0’
library("tidyverse"); packageVersion("tidyverse") #‘1.3.2’
library("DESeq2"); packageVersion("DESeq2") #‘1.36.0'
library("dendextend"); packageVersion("dendextend") #‘1.16.0’
library("viridis"); packageVersion("viridis") #‘0.6.2’
library("ANCOMBC"); packageVersion("ANCOMBC") #'1.6.4'

#### LOAD AND PREPARE DATA ####
setwd("Data/210212_Data")

count_tab <- read.table("ASVs_counts-no-contam_silva_sp_0.4.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")

#Remove mock communities and negative controls
dim(count_tab) #6610 255
colnames(count_tab)
mock_and_negative <- count_tab[,c(234:239)]
count_tab_clean <- count_tab[,c(1:233,240:255)]
dim(count_tab_clean) #6610 249
colnames(count_tab_clean)
#write.csv(count_tab_clean,"count_tab_clean.csv")

#Loading in taxonomy
tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam_silva_sp_0.4.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
dim(tax_tab) #6610 7

#Load metadata
#metadata_all <- read.table("Data/EggshellMicrobiomeProject.csv",header=TRUE,sep=",",row.names=1)
metadata_all <- read.table("Data/EggshellMicrobiomeProject_forR.txt",header=TRUE,sep=",",row.names=1) #Using text file because prevents sample names changing to dates
head(metadata_all) #Samples as rows, variables as columns
dim(metadata_all) #255 83
#Remove mocks and negative controls just in case
rownames(metadata_all)
mock_and_negative_meta <- metadata_all[c(234:239),]
metadata_all_clean <- metadata_all[c(1:233,240:255),]
dim(metadata_all_clean) #249 83

#Make subsets for kaki and hihi with and without negatives and mocks
metadata_hihi <- metadata_all[c(93:255),]
metadata_kaki <- metadata_all[c(1:92,234:239),]
metadata_hihi_clean <- metadata_all[c(93:233,240:255),]
metadata_kaki_clean <- metadata_all[c(1:92),]
dim(metadata_hihi) #163 83
dim(metadata_kaki) #98 83
dim(metadata_hihi_clean) #157 83
dim(metadata_kaki_clean) #92 83

count_tab_hihi <- count_tab[,c(93:255)]
count_tab_kaki <- count_tab[,c(1:92,234:239)]
count_tab_hihi_clean <- count_tab[,c(93:233,240:255)]
count_tab_kaki_clean <- count_tab[,c(1:92)]
dim(count_tab_hihi) #6610 163
dim(count_tab_kaki) #6610 98
dim(count_tab_hihi_clean) #6610 157
dim(count_tab_kaki_clean) #6610 92

## BETA DIVERSITY (https://astrobiomike.github.io/amplicon/dada2_workflow_ex)

#Need to normalise our samples to account for differences in sampling depth. Subsampling to lowest samples depth or turning counts into proportions = generally shunned. Going to use variance stabilising transformation instead using DESeq2 (see McMurdie and Holmes 2014)

## All samples
#Make a DESeq2 object
deseq_counts <- DESeqDataSetFromMatrix(count_tab_clean, colData = metadata_all_clean, design = ~species) #N.B. need to include something for colData and design here to make it work, but don't matter for purposes of transforming the data
#Warning message: In DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
#Pull out transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)
#Calculating Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))

## Hierarchical clustering
#Use Euclidean distance matrix to make and plot a hierarchical clustering of samples
euc_clust <- hclust(euc_dist, method="ward.D2")
plot(euc_clust)

#Add a colour column to use in this plot
metadata_all_clean_colour <- metadata_all_clean
metadata_all_clean_colour$colour <- metadata_all_clean_colour$species
head(metadata_all_clean_colour)
metadata_all_clean_colour$colour<-gsub("Hihi","red",metadata_all_clean_colour$colour)
metadata_all_clean_colour$colour<-gsub("Kaki","blue",metadata_all_clean_colour$colour)

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(metadata_all_clean_colour$colour[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols
plot(euc_dend, ylab="VST Euc. dist.") #Very strong hihi/kaki clusters - looks like some samples of each in the wrong clusters? Looks like same samples that could see in boxplot in qiime2 - hihi seem ok, but 19-7, 19-8, 19-9, 19-10, 19-160, 19-161, 19-162, 19-163, 19-164, 19-165, 19-166 (4 different entire clutches) spread throughout and within hihi cluster

colnames(metadata_all_clean)
## Ordination
#Making a phyloseq object with the transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
metadata_all_clean_colour_phy <- sample_data(metadata_all_clean_colour)
vst_physeq <- phyloseq(vst_count_phy, metadata_all_clean_colour_phy)

#Generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues #Can scale the axes according to their magnitude of separating apart the samples
plot_ordination(vst_physeq, vst_pcoa, color="species") + 
  geom_point(size=1) + labs(col="species") + 
  geom_text(aes(label=rownames(metadata_all_clean_colour), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(metadata_all_clean_colour$colour[order(metadata_all_clean_colour$species)])) + 
  theme_bw() + theme(legend.position="none") #Same hihi/kaki clusters as seen in dendogram

#Making a subset with these non-clustering kaki samples removed just in case
metadata_all_clean_refined <- metadata_all_clean[!(row.names(metadata_all_clean) %in% c("19-7", "19-8", "19-9", "19-10", "19-160", "19-161", "19-162", "19-163", "19-164", "19-165", "19-166")),]
rownames(metadata_all_clean_refined)
dim(metadata_all_clean_refined) #238 83
colnames(count_tab_clean)
count_tab_clean_refined <- count_tab_clean[,!names(count_tab_clean) %in% c("19-7", "19-8", "19-9", "19-10", "19-160", "19-161", "19-162", "19-163", "19-164", "19-165", "19-166")]
colnames(count_tab_clean_refined)
dim(count_tab_clean_refined) #6610 238

## Hihi only
#Make a DESeq2 object
deseq_counts_hihi <- DESeqDataSetFromMatrix(count_tab_hihi_clean, colData = metadata_hihi_clean, design = ~clutch.nest_id_corrected) #N.B. need to include something for colData and design here to make it work, but don't matter for purposes of transforming the data
#Note: levels of factors in the design contain characters other than letters, numbers, '_' and '.'. It is recommended (but not required) to use only letters, numbers, and delimiters '_' or '.', as these are safe characters for column names in R. [This is a message, not a warning or an error]
#Warning message: In DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors
deseq_counts_vst_hihi <- varianceStabilizingTransformation(deseq_counts_hihi)
#Pull out transformed table
vst_trans_count_tab_hihi <- assay(deseq_counts_vst_hihi)
#Calculating Euclidean distance matrix
euc_dist_hihi <- dist(t(vst_trans_count_tab_hihi))

## Hierarchical clustering
#Use Euclidean distance matrix to make and plot a hierarchical clustering of samples
euc_clust_hihi <- hclust(euc_dist_hihi, method="ward.D2")
plot(euc_clust_hihi)

#Add a colour column to use in this plot
metadata_hihi_clean_colour <- metadata_hihi_clean
metadata_hihi_clean_colour$colour <- metadata_hihi_clean_colour$hatched.unhatched
head(metadata_hihi_clean_colour)
metadata_hihi_clean_colour$colour<-gsub("Unhatched","red",metadata_hihi_clean_colour$colour)
metadata_hihi_clean_colour$colour<-gsub("Hatched","blue",metadata_hihi_clean_colour$colour)

euc_dend_hihi <- as.dendrogram(euc_clust_hihi, hang=0.1)
dend_cols_hihi <- as.character(metadata_hihi_clean_colour$colour[order.dendrogram(euc_dend_hihi)])
labels_colors(euc_dend_hihi) <- dend_cols_hihi
plot(euc_dend_hihi, ylab="VST Euc. dist.") #No obvious clustering for hatched/unhatched

## Ordination
#Making a phyloseq object with the transformed table
vst_count_phy_hihi <- otu_table(vst_trans_count_tab_hihi, taxa_are_rows=T)
metadata_hihi_clean_colour_phy <- sample_data(metadata_hihi_clean_colour)
vst_physeq_hihi <- phyloseq(vst_count_phy_hihi, metadata_hihi_clean_colour_phy)

#Generating and visualizing the PCoA with phyloseq
vst_pcoa_hihi <- ordinate(vst_physeq_hihi, method="MDS", distance="euclidean")
eigen_vals_hihi <- vst_pcoa_hihi$values$Eigenvalues #Can scale the axes according to their magnitude of separating apart the samples
plot_ordination(vst_physeq_hihi, vst_pcoa_hihi, color="hatched.unhatched") + 
  geom_point(size=1) + labs(col="hatched.unhatched") + 
  geom_text(aes(label=rownames(metadata_hihi_clean_colour), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(metadata_hihi_clean_colour$colour[order(metadata_hihi_clean_colour$hatched.unhatched)])) + 
  theme_bw() + theme(legend.position="none") #No obvious clustering for hatched/unhatched, B22-1A-3 stands out

## Kakī only
#Taking out the samples clustering with hihi to use in comparisons just in case
metadata_kaki_clean_refined <- metadata_kaki_clean[!(row.names(metadata_kaki_clean) %in% c("19-7", "19-8", "19-9", "19-10", "19-160", "19-161", "19-162", "19-163", "19-164", "19-165", "19-166")),]
rownames(metadata_kaki_clean_refined)
dim(metadata_kaki_clean_refined) #81 83
count_tab_kaki_clean_refined <- count_tab_kaki_clean[,!names(count_tab_kaki_clean) %in% c("19-7", "19-8", "19-9", "19-10", "19-160", "19-161", "19-162", "19-163", "19-164", "19-165", "19-166")]
colnames(count_tab_kaki_clean_refined)
dim(count_tab_kaki_clean_refined) #6610 92

#metadata_kaki_clean <- metadata_kaki_clean_refined #Used to check what happens if non-clustered kaki removed
#count_tab_kaki_clean <- count_tab_kaki_clean_refined #Used to check what happens if non-clustered kaki removed

#Make a DESeq2 object
deseq_counts_kaki <- DESeqDataSetFromMatrix(count_tab_kaki_clean, colData = metadata_kaki_clean, design = ~clutch.nest_id_corrected) #N.B. need to include something for colData and design here to make it work, but don't matter for purposes of transforming the data
#Note: levels of factors in the design contain characters other than letters, numbers, '_' and '.'. It is recommended (but not required) to use only letters, numbers, and delimiters '_' or '.', as these are safe characters for column names in R. [This is a message, not a warning or an error]
#Warning message: In DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors
deseq_counts_vst_kaki <- varianceStabilizingTransformation(deseq_counts_kaki)
#Pull out transformed table
vst_trans_count_tab_kaki <- assay(deseq_counts_vst_kaki)
#Calculating Euclidean distance matrix
euc_dist_kaki <- dist(t(vst_trans_count_tab_kaki))

## Hierarchical clustering
#Use Euclidean distance matrix to make and plot a hierarchical clustering of samples
euc_clust_kaki <- hclust(euc_dist_kaki, method="ward.D2")
plot(euc_clust_kaki)
colnames (metadata_kaki)
#Add a colour column to use in this plot
metadata_kaki_clean_colour <- metadata_kaki_clean
metadata_kaki_clean_colour$colour <- metadata_kaki_clean_colour$hatched.unhatched
head(metadata_kaki_clean_colour)
metadata_kaki_clean_colour$colour<-gsub("Unhatched","red",metadata_kaki_clean_colour$colour)
metadata_kaki_clean_colour$colour<-gsub("Hatched","blue",metadata_kaki_clean_colour$colour)
metadata_kaki_clean_colour.wc <- metadata_kaki_clean
metadata_kaki_clean_colour.wc$colour <- metadata_kaki_clean_colour.wc$wild.captive
head(metadata_kaki_clean_colour.wc)
metadata_kaki_clean_colour.wc$colour<-gsub("Captive","red",metadata_kaki_clean_colour.wc$colour)
metadata_kaki_clean_colour.wc$colour<-gsub("Wild","blue",metadata_kaki_clean_colour.wc$colour)

euc_dend_kaki <- as.dendrogram(euc_clust_kaki, hang=0.1)
euc_dend_kaki.wc <- as.dendrogram(euc_clust_kaki, hang=0.1)
dend_cols_kaki <- as.character(metadata_kaki_clean_colour$colour[order.dendrogram(euc_dend_kaki)])
dend_cols_kaki.wc <- as.character(metadata_kaki_clean_colour.wc$colour[order.dendrogram(euc_dend_kaki.wc)])
labels_colors(euc_dend_kaki) <- dend_cols_kaki
labels_colors(euc_dend_kaki.wc) <- dend_cols_kaki.wc
plot(euc_dend_kaki, ylab="VST Euc. dist.") #No apparent clustering for hatched/unhatched
plot(euc_dend_kaki.wc, ylab="VST Euc. dist.") #Maybe some slight clustering for wild/captive

## Ordination
#Making a phyloseq object with the transformed table
vst_count_phy_kaki <- otu_table(vst_trans_count_tab_kaki, taxa_are_rows=T)
metadata_kaki_clean_colour_phy <- sample_data(metadata_kaki_clean_colour)
metadata_kaki_clean_colour.wc_phy <- sample_data(metadata_kaki_clean_colour.wc)
vst_physeq_kaki <- phyloseq(vst_count_phy_kaki, metadata_kaki_clean_colour_phy)
vst_physeq_kaki.wc <- phyloseq(vst_count_phy_kaki, metadata_kaki_clean_colour.wc_phy)

#Generating and visualizing the PCoA with phyloseq
vst_pcoa_kaki <- ordinate(vst_physeq_kaki, method="MDS", distance="euclidean")
vst_pcoa_kaki.wc <- ordinate(vst_physeq_kaki.wc, method="MDS", distance="euclidean")
eigen_vals_kaki <- vst_pcoa_kaki$values$Eigenvalues #Can scale the axes according to their magnitude of separating apart the samples
eigen_vals_kaki.wc <- vst_pcoa_kaki.wc$values$Eigenvalues #Can scale the axes according to their magnitude of separating apart the samples
plot_ordination(vst_physeq_kaki, vst_pcoa_kaki, color="hatched.unhatched") + 
  geom_point(size=1) + labs(col="hatched.unhatched") + 
  geom_text(aes(label=rownames(metadata_kaki_clean_colour), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(metadata_kaki_clean_colour$colour[order(metadata_kaki_clean_colour$hatched.unhatched)])) + 
  theme_bw() + theme(legend.position="none") #Maybe some slight clustering, 19-85 and 19-60 appear to be outliers
plot_ordination(vst_physeq_kaki.wc, vst_pcoa_kaki.wc, color="wild.captive") + 
  geom_point(size=1) + labs(col="wild.captive") + 
  geom_text(aes(label=rownames(metadata_kaki_clean_colour.wc), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(metadata_kaki_clean_colour.wc$colour[order(metadata_kaki_clean_colour.wc$wild.captive)])) + 
  theme_bw() + theme(legend.position="none") #Maybe some slight clustering for hatched/unhatched, 19-85 and 19-60 appear to be outliers

## ALPHA DIVERSITY (https://astrobiomike.github.io/amplicon/dada2_workflow_ex)

## Rarefaction curves - note = rarefaction curves should not be used to estimate total richness of a sample, but can be useful in exploring the data
rarecurve(t(count_tab_clean), step=100, col=metadata_all_clean_colour$colour, lwd=2, ylab="ASVs", label=F) #Use the t() to transpose the table because the vegan package expects rows to be samples and columns to be ASVs
#Kaki (blue) seem to generally have greater richness (unique number of sequences recovered) than hihi
abline(v=(min(rowSums(t(count_tab_clean))))) #Add a vertical line at the fewest seqs in any sample

## Richness and diversity estimates
#Plotting Chao1 richness estimates and Shannon diversity values.
#Create a phyloseq object using  un-transformed count table
count_tab_clean_phy <- otu_table(count_tab_clean, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
ASV_physeq <- phyloseq(count_tab_clean_phy, tax_tab_phy, metadata_all_clean_colour_phy)

#Call the plot_richness() function on phyloseq object
plot_richness(ASV_physeq, color="species", measures=c("Observed", "Shannon")) + 
  scale_color_manual(values=unique(metadata_all_clean_colour$colour[order(metadata_all_clean_colour$species)])) +
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) #Very clear clustering of hihi and kaki in both Observed features and Shannon plots

#Plot grouping by species and colouring with hatched.unhatched
plot_richness(ASV_physeq, x="species", color="hatched.unhatched", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(metadata_all_clean_colour$colour[order(metadata_all_clean_colour$species)])) +
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) #Nothing super obvious, maybe slight lower alpha diversity in unhatched eggs compared to hatched eggs. Hihi lower Chao1 and higher Shannon diversity than kaki --> hihi lower richness, but higher evenness


## TAXONOMIC SUMMARIES (https://astrobiomike.github.io/amplicon/dada2_workflow_ex)
# Note: advises to pull out ASVs that become important in the story and BLAST them, and maybe integrate into phylogenetic trees to get a more robust idea of who they are most closely related to

#Make a count table that has summed all ASVs that were in the same phylum
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Phylum")) 
#Make a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="Phylum"))[,"Phylum"]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

#Already filtered out sequences not assigned taxonomy at phylum level so don't need to account for these

#Remove the Proteobacteria, so we can next add them back in broken down by class
temp_major_taxa_counts_tab <- phyla_counts_tab[!row.names(phyla_counts_tab) %in% "Proteobacteria", ]

class_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Class")) #Make count table broken down by class (contains classes beyond the # Proteobacteria too at this point)
class_tax_phy_tab <- tax_table(tax_glom(ASV_physeq, taxrank="Class")) #Make a table that holds the phylum and class level info

phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("Phylum"=phy_tmp_vec, "Class"=class_tmp_vec, row.names = rows_tmp)
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$Phylum == "Proteobacteria", "Class"]) #Make a vector of just the Proteobacteria classes

rownames(class_counts_tab) <- as.vector(class_tax_tab$Class) #Changing the row names like above so that they correspond to the taxonomy, rather than an ASV identifier
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ] #Make a table of the counts of the Proteobacterial classes

#There may be some sequences that were resolved to the level # of Proteobacteria, but not any further, which would be missing from class table - we can find the sum of them by subtracting the proteo class count table # from just the Proteobacteria row from the original phylum-level count table
proteo_no_class_annotated_counts <- phyla_counts_tab[row.names(phyla_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)

#Combine the tables:
major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria"=proteo_no_class_annotated_counts)

#Check we didn't miss any other sequences, we can compare the column sums to see if they are the same, if "TRUE", we know nothing fell through the cracks
identical(colSums(major_taxa_counts_tab), colSums(count_tab_clean)) #TRUE!!!

#Generate a proportions table for summarizing:
major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)
#Check the dimensions of this table at this point
dim(major_taxa_proportions_tab) #37 249 --> currently 37 phyla, may be too busy in the figure
colnames(major_taxa_proportions_tab)
rownames(major_taxa_proportions_tab)

#Keeping taxa that make up greater than 2% in any individual sample
temp_filt_major_taxa_proportions_tab <- as.data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 2, ]) #Need to use as.data.frame to stop R changing the sample names (adding X in front of number and changing - to .)
dim(temp_filt_major_taxa_proportions_tab) #12 249 --> drops to 12, much easier to plot
colnames(temp_filt_major_taxa_proportions_tab)

#Add the rest of the taxa (below the 2% cut-off) to a row called "Other" (which will also keep  totals at 100%)
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)

#Make plots using the proportions table

#First make a copy of table that's safe for manipulating
filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab

#Make the data easier to use with ggplot
filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_tab_for_plot) #Add a column of the taxa names so that it is within the table
filt_major_taxa_proportions_tab_for_plot.g <- pivot_longer(filt_major_taxa_proportions_tab_for_plot, !Major_Taxa, names_to = "Sample", values_to = "Proportion") %>% data.frame() #Transform the table into narrow, or long, format
head(filt_major_taxa_proportions_tab_for_plot.g)
dim(filt_major_taxa_proportions_tab_for_plot.g) #3237 3
head(filt_major_taxa_proportions_tab_for_plot)

#Want a table with "color" and "species" of each sample to merge into our plotting table so we can use that more easily in our plotting function
sample_info_for_merge<-data.frame("Sample"=row.names(metadata_all_clean_colour), "Species"=metadata_all_clean_colour$species, "Colour"=metadata_all_clean_colour$colour, stringsAsFactors=F)
filt_major_taxa_proportions_tab_for_plot.g2 <- merge(filt_major_taxa_proportions_tab_for_plot.g, sample_info_for_merge) #Merge this with the plotting table we just made
head(filt_major_taxa_proportions_tab_for_plot.g)
head(sample_info_for_merge)
head(filt_major_taxa_proportions_tab_for_plot.g2) #Worked!
#Stacked bar charts for each taxon per sample:
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples") #Works but too many samples to be meaningful

#Boxplot with each box as a major taxon
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(Major_Taxa, Proportion)) +
  geom_jitter(aes(color=factor(Species), shape=factor(Species)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_tab_for_plot.g2$Colour[order(filt_major_taxa_proportions_tab_for_plot.g2$Species)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank()) +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="All samples")

#Boxplot with each box as a major taxon - species independently

#Subset our plotting table
hihi_sample_IDs <- row.names(metadata_all_clean_colour)[metadata_all_clean_colour$species == "Hihi"]
kaki_sample_IDs <- row.names(metadata_all_clean_colour)[metadata_all_clean_colour$species == "Kaki"]
filt_major_taxa_proportions_hihi_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% hihi_sample_IDs, ]
filt_major_taxa_proportions_kaki_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% kaki_sample_IDs, ]

#Hihi samples
ggplot(filt_major_taxa_proportions_hihi_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(Species), shape=factor(Species)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_hihi_only_tab_for_plot.g$Colour[order(filt_major_taxa_proportions_hihi_only_tab_for_plot.g$Species)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) + # moved legend to top 
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Hihi samples only")

#Kaki samples
ggplot(filt_major_taxa_proportions_kaki_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(Species), shape=factor(Species)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_kaki_only_tab_for_plot.g$Colour[order(filt_major_taxa_proportions_kaki_only_tab_for_plot.g$Species)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) + # moved legend to top 
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Kakī samples only")

#Pie charts
hihi_sample_major_taxa_proportion_tab <- filt_major_taxa_proportions_hihi_only_tab_for_plot.g[, c(1:3)] %>% pivot_wider(names_from = Major_Taxa, values_from = Proportion) %>% column_to_rownames("Sample") %>% t() %>% data.frame()
kaki_sample_major_taxa_proportion_tab <- filt_major_taxa_proportions_kaki_only_tab_for_plot.g[, c(1:3)] %>% pivot_wider(names_from = Major_Taxa, values_from = Proportion) %>% column_to_rownames("Sample") %>% t() %>% data.frame()
hihi_sample_summed_major_taxa_proportions_vec <- rowSums(hihi_sample_major_taxa_proportion_tab)
kaki_sample_summed_major_taxa_proportions_vec <- rowSums(kaki_sample_major_taxa_proportion_tab)
hihi_sample_major_taxa_summary_tab <- data.frame("Major_Taxa"=names(hihi_sample_summed_major_taxa_proportions_vec), "Proportion"=hihi_sample_summed_major_taxa_proportions_vec, row.names=NULL)
kaki_sample_major_taxa_summary_tab <- data.frame("Major_Taxa"=names(kaki_sample_summed_major_taxa_proportions_vec), "Proportion"=kaki_sample_summed_major_taxa_proportions_vec, row.names=NULL)
#Hihi only
ggplot(data.frame(hihi_sample_major_taxa_summary_tab), aes(x="Hihi samples", y=Proportion, fill=Major_Taxa)) + 
  geom_bar(width=1, stat="identity") +
  coord_polar("y") +
  scale_fill_viridis(discrete=TRUE) +
  ggtitle("Hihi samples only") +
  theme_void() +
  theme(plot.title = element_text(hjust=0.5), legend.title=element_blank())
#Kaki only
ggplot(data.frame(kaki_sample_major_taxa_summary_tab), aes(x="Kakī samples", y=Proportion, fill=Major_Taxa)) + 
  geom_bar(width=1, stat="identity") +
  coord_polar("y") +
  scale_fill_viridis(discrete=TRUE) +
  ggtitle("Kakī samples only") +
  theme_void() +
  theme(plot.title = element_text(hjust=0.5), legend.title=element_blank())

## BETADISPER AND PERMUTATIONAL ANOVA (https://astrobiomike.github.io/amplicon/dada2_workflow_ex)
#Permutational ANOVA test used to see if any available information is indicative of community structure
#Can test if there is a significant different between sample types
#Function 'adonis' can tell us if there is a statistical difference between groups, but has an assumption which must be checked with 'betadisper' (needs to be sufficient level of homogeneity of dispersion within groups, if not adonis can be unreliable)

#All samples
anova(betadisper(euc_dist, metadata_all_clean$species)) #p = 3.055e-15 --> SIGNIFICANT = means there is a difference between group dispersions, so can't trust the results of an adonis (permutational anova) because the assumption of homogeneous within-group dispersions is not met
#Kaki - wild vs. captive
anova(betadisper(euc_dist_kaki, metadata_kaki_clean$wild.captive)) #p = 0.5145 --> NOT SIGNIFICANT
#Kaki - hatched vs. unhatched
anova(betadisper(euc_dist_kaki, metadata_kaki_clean$hatched.unhatched)) #p = 0.05069 --> NOT SIGNIFICANT (only just!!!)
#Hihi - hatched vs. unhatched
anova(betadisper(euc_dist_hihi, metadata_hihi_clean$hatched.unhatched)) #p = 0.4688 --> NOT SIGNIFICANT
#Not significant means can test if the groups host statistically different communities based on adonis (meet assumption)
adonis2(euc_dist_kaki~metadata_kaki_clean$wild.captive) #p = 0.004 - significant difference in microbial communities in wild vs. captive kaki
adonis2(euc_dist_kaki~metadata_kaki_clean$hatched.unhatched) #p = 0.984 - no significant difference in microbial communities in hatched vs. unhatched kaki
metadata_hihi_clean_noHatchNA <- metadata_hihi_clean[-which(is.na(metadata_hihi_clean$hatched.unhatched)),]
HatchNA <- metadata_hihi_clean[is.na(metadata_hihi_clean$hatched.unhatched),]
row.names(HatchNA)
class(vst_trans_count_tab_hihi)
vst_trans_count_tab_hihi_noHatchNA <- vst_trans_count_tab_hihi[,!colnames(vst_trans_count_tab_hihi) %in% c("AB-3-1","AB-3-2","AB-3-3","AB-3-4","B1-19-1","B1-19-2","B1-19-3","B1-19-4","B1-28-1","B1-28-2","B21-13-1","B21-13-2","B21-13-3","B22-25-1","B22-25-3","B22-25-4","B22-28-1","B22-28-2","B22-28-3","B22-28-4","B4-3-1","B4-3-2","B4-3-3","B4-3-4","BH-2-1","BH-2-3","BH-2-4","CH-1-1","CH-1-2","CH-1-3","CH-1-4","KB-1-1","KB-1-2","KB-1-3","KB-1-4","LHV-2-1","LHV-2-2","LHV-2-3","LHV-2-4","SC-4-1","SC-4-2","SC-4-3","SC-4-4","SD-3-1","SD-3-2","SD-5-5","SD-5-6","SD-5-7","SV-4-1","SV-4-2","SV-4-3")]
colnames(vst_trans_count_tab_hihi_noHatchNA)
euc_dist_hihi_noHatchNA <- dist(t(vst_trans_count_tab_hihi_noHatchNA))
adonis2(euc_dist_hihi_noHatchNA~metadata_hihi_clean_noHatchNA$hatched.unhatched) #p = 0.898 - no significant difference in microbial communities in hatched vs. unhatched hihi

#Plot
#Generate and visualize the PCoA with phyloseq
kaki_eigen_vals <- vst_pcoa_kaki.wc$values$Eigenvalues #Scale the axes according to their magnitude of separating apart the samples
#Make new ordination of just kaki with adonis statistic
plot_ordination(vst_physeq_kaki.wc, vst_pcoa_kaki.wc, color="wild.captive") + 
  geom_point(size=1) + labs(col="wild.captive") + 
  geom_text(aes(label=rownames(metadata_kaki_clean_colour.wc), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=50, y=128, label="Wild vs captive") +
  annotate("text", x=50, y=122, label="Permutational ANOVA = 0.003") + 
  coord_fixed(sqrt(kaki_eigen_vals[2]/kaki_eigen_vals[1])) + ggtitle("PCoA - kakī only") + 
  scale_color_manual(values=unique(metadata_kaki_clean_colour.wc$colour[order(metadata_kaki_clean_colour.wc$wild.captive)])) + 
  theme_bw() + theme(legend.position="none")

## DIFFERENTIAL ABUNDANCE ANALYSIS WITH DESEQ2

# Can use DESeq2 to find out which ASVs (and possibly which taxa) are contributing to that difference
# QUESTION: Can I use DESeq2 to compare between kaki and hihi when it didn't meet the assumptions for adonis (not significant according to betadisper)

#Making a kaki-only phyloseq object of non-transformed values (as that is what DESeq2 operates on)
kaki_count_phy <- otu_table(count_tab_clean[, colnames(count_tab_clean) %in% kaki_sample_IDs], taxa_are_rows=T)
#write.csv(kaki_count_phy,"kaki_count_phy.csv")
class(kaki_count_phy)
kaki_count_phy_no0 <- kaki_count_phy[rowSums(kaki_count_phy)>0,] #Removing ASVs with zero counts in all samples to check if makes a difference - doesn't
dim(kaki_count_phy_no0) #4064   92
kaki_count_physeq <- phyloseq(kaki_count_phy, metadata_kaki_clean_colour.wc_phy)
dim(kaki_count_phy) #6610   92
kaki_count_physeq #6610   92
otu_table(kaki_count_physeq)
#Convert phyloseq object to a deseq object
kaki_deseq <- phyloseq_to_deseq2(kaki_count_physeq, ~wild.captive)
dim(kaki_deseq) #6610   92
#Run deseq standard analysis:
kaki_deseq <- DESeq(kaki_deseq,test="Wald",fitType="parametric")
kaki_deseq #6610   92

####FROM HERE - NOT WORKING?

#Pull out results table - specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_wild_vs_captive <- results(kaki_deseq, alpha=0.05, contrast=c("wild.captive", "Wild", "Captive"))
summary(deseq_res_wild_vs_captive) #Out of 2670 ASVs with adj-p <0.05 shows nothing increased or decreased?


# this tells us out of ~1,800 ASVs, with adj-p < 0.01, there are 7 increased when comparing altered basalts to glassy basalts, and about 6 decreased
# "decreased" in this case means at a lower count abundance in the altered basalts than in the glassy basalts, and "increased" means greater proportion in altered than in glassy
# remember, this is done with a drastically reduced dataset, which is hindering the capabilities here quite a bit i'm sure

# let's subset this table to only include these that pass our specified significance level
sigtab_res_deseq_altered_vs_glassy <- deseq_res_altered_vs_glassy[which(deseq_res_altered_vs_glassy$padj < 0.01), ]

# now we can see this table only contains those we consider significantly differentially abundant
summary(sigtab_res_deseq_altered_vs_glassy) 

# next let's stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_deseq_altered_vs_glassy_with_tax <- cbind(as(sigtab_res_deseq_altered_vs_glassy, "data.frame"), as(tax_table(ASV_physeq)[row.names(sigtab_res_deseq_altered_vs_glassy), ], "matrix"))

# and now let's sort that table by the baseMean column
sigtab_deseq_altered_vs_glassy_with_tax[order(sigtab_deseq_altered_vs_glassy_with_tax$baseMean, decreasing=T), ]

# this puts a sequence derived from a Rhizobiales at the second to highest (first is unclassified) that was detected in ~7 log2fold greater abundance in the glassy basalts than in the highly altered basalts


#Make a phyloseq object of non-transformed values (as that is what DESeq2 operates on)
all_count_physeq <- phyloseq(count_tab_clean_phy, metadata_all_clean_colour_phy)
#write.csv(count_tab_clean_phy,"all_count_physeq.csv")
#Convert phyloseq object to a deseq object
all_deseq <- phyloseq_to_deseq2(all_count_physeq, ~species)
#Run deseq standard analysis
all_deseq <- DESeq(all_deseq)

#Pull out results table - specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_hihi_vs_kaki <- results(all_deseq, alpha=0.05, contrast=c("species", "Hihi", "Kaki"))
summary(deseq_res_hihi_vs_kaki) #Out of 3230 ASVs, with adj-p < 0.5 there are 6 increased (i.e. greater proportion in hihi than kaki) when comparing hihi to kaki, and 29 decreased (i.e. at a lower count abundance in the hihi than in the kaki samples) = 35 significant in total
#With adj-p <0.01 4 increased and 27 decreased

#Subset this table to only include ones that pass p < 0.05
sigtab_res_deseq_res_hihi_vs_kaki <- deseq_res_hihi_vs_kaki[which(deseq_res_hihi_vs_kaki$padj < 0.05), ]
summary(sigtab_res_deseq_res_hihi_vs_kaki) #Table only contains ASVs considered significantly differentially abunadant
#Stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_res_deseq_res_hihi_vs_kaki_with_tax <- cbind(as(sigtab_res_deseq_res_hihi_vs_kaki, "data.frame"), as(tax_table(ASV_physeq)[row.names(sigtab_res_deseq_res_hihi_vs_kaki), ], "matrix"))
dim(sigtab_res_deseq_res_hihi_vs_kaki_with_tax) #35 13

#Sort by the log2FoldChange column
sigtab_res_deseq_res_hihi_vs_kaki_with_tax[order(sigtab_res_deseq_res_hihi_vs_kaki_with_tax$log2FoldChange, decreasing=T), ] #Shows the taxonomy of the 35 significant ASVs = negative log2FoldChange means lower count abundance in hihi vs. kaki, positive log2FoldChange means greater count abundance in hihi vs. kaki

#write.csv(sigtab_res_deseq_res_hihi_vs_kaki_with_tax,"deseq_results_dada2.csv")

## Same as in QIIME2 pipeline

#From http://joey711.github.io/phyloseq-extensions/DESeq2.html
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab_res_deseq_res_hihi_vs_kaki_with_tax$log2FoldChange, sigtab_res_deseq_res_hihi_vs_kaki_with_tax$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_res_deseq_res_hihi_vs_kaki_with_tax$Phylum = factor(as.character(sigtab_res_deseq_res_hihi_vs_kaki_with_tax$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab_res_deseq_res_hihi_vs_kaki_with_tax$log2FoldChange, sigtab_res_deseq_res_hihi_vs_kaki_with_tax$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_res_deseq_res_hihi_vs_kaki_with_tax$Genus = factor(as.character(sigtab_res_deseq_res_hihi_vs_kaki_with_tax$Genus), levels=names(x))
ggplot(sigtab_res_deseq_res_hihi_vs_kaki_with_tax, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

## ANCOM-BC (https://microbiome.github.io/course_2021_radboud/differential-abundance-analysis.html)

# Requires phylosec format
#ancombc_res = ancombc2(all_count_physeq,  p_adj_method = "fdr", 
#                       tax_level = "Species",
#                       fix_formula = "species",
#                       lib_cut = 0, 
#                       group = "species", 
#                       struc_zero = TRUE, 
#                       neg_lb = TRUE, 
#                       alpha = 0.05, 
#                       global = TRUE)
#ancombc_res$res$diff_speciesKaki




#### HAND-OFF TO PHYLOSEQ
theme_set(theme_bw())

metadata_all <- read.table("Data/EggshellMicrobiomeProject.csv",header=TRUE,sep=",",row.names=1)
metadata_all_noNA <- read.table("Data/EggshellMicrobiomeProject_noNA.csv",header=TRUE,sep=",",row.names=1)
metadata_hihi <- read.table("Data/EggshellMicrobiomeProject_hihi.csv",header=TRUE,sep=",",row.names=1)
metadata_hihi_noNA <- read.table("Data/EggshellMicrobiomeProject_hihi_noNA.csv",header=TRUE,sep=",",row.names=1)
metadata_kaki <- read.table("Data/EggshellMicrobiomeProject_kaki.csv",header=TRUE,sep=",",row.names=1)
metadata_kaki_noNA <- read.table("Data/EggshellMicrobiomeProject_kaki_noNA.csv",header=TRUE,sep=",",row.names=1)

ps <- phyloseq(otu_table(seqtab.nochim,taxa_are_rows=FALSE),
               sample_data(metadata_all),
               tax_table(taxa_gg))
sample_names(ps)
ps <- prune_samples(sample_names(ps) != "MC-1", ps) # Remove mock community samples
ps <- prune_samples(sample_names(ps) != "MC-2", ps) # Remove mock community samples
ps <- prune_samples(sample_names(ps) != "MC-3", ps) # Remove mock community samples
sample_names(ps)
#ps <- prune_samples(sample_names(ps) != "NTC-1", ps) # Remove mock community samples
#ps <- prune_samples(sample_names(ps) != "NTC-2", ps) # Remove mock community samples
#ps <- prune_samples(sample_names(ps) != "NTC-3", ps) # Remove mock community samples
#sample_names(ps)

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

plot_richness(ps, x="species", measures=c("Shannon", "Simpson"), color="wild.captive")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="species", title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="species", fill="Family") + facet_wrap(~species, scales="free_x") #not working
