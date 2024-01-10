setwd("/home/eteleeb/projects")
#library("iCluster")
library("iClusterPlus")
library("lattice")
library("gplots")
library(ggplot2)
library(GenomicRanges)
library(gridExtra)
library(parallel)
library(reshape2)
library(data.table)
library(gplots)
library(ggsignif)
library(dplyr)
library(VennDiagram)
library(RColorBrewer)
library(venn)
library(enrichR)
library(clusterProfiler)
library(org.Hs.eg.db)
library("AnnotationDbi")
library(gage)
library(pathview)
library(gageData)

## good colors c("#00AFBB", "#E7B800")

# set up directories 
dir.create('omics_integration/iCluster_output')
dir.create('omics_integration/iCluster_output/Figures')

############################### read all required datasets  ############################### 
## read clinical data
#mend.clin <- read.csv('/40/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/03.-Phenotype/2019_09_25_WashU_MendelianVsSporadics_clinical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F)
mend.clin <- read.csv('/home/data/WashU_Data/bulkRNASeq/201703_MendelianVsSporadics/03.-Phenotype/2019_09_25_WashU_MendelianVsSporadics_clinical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F) 

#sunshine.clin <- read.csv('/40/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/03.-Phenotype/2020_10_13_WashU_SUNSHINE_clinical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F) 
sunshine.clin <- read.csv('/home/data/WashU_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/03.-Phenotype/2020_10_13_WashU_SUNSHINE_clinical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F) 

sunshine.clin$Subj_ID <- gsub("MAP_", "", sunshine.clin$Subj_ID)
sunshine.clin$Subj_ID <- gsub("MAC_", "", sunshine.clin$Subj_ID)
sunshine.clin$Subj_ID <- gsub("DIAN_", "", sunshine.clin$Subj_ID)
clin.shared.cols <- intersect(colnames(sunshine.clin), colnames(mend.clin))
mend.clin <- mend.clin[, clin.shared.cols]
mend.clin$cohort <- 'M'
sunshine.clin <- sunshine.clin[, clin.shared.cols]
sunshine.clin$cohort <- 'S'
clin.data <- unique(rbind(mend.clin, sunshine.clin))
clin.data[!is.na(clin.data$Status) & clin.data$Status=="ADAD_carrier", "Status"] = "ADAD_Carrier"
clin.data <- unique(clin.data[, !colnames(clin.data) %in% c("study","OtherNPDx", "cohort")])
dim(clin.data)

## read count & expression 
mend.cts <- as.data.frame(fread('/home/data/WashU_Data/bulkRNASeq/201703_MendelianVsSporadics/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_count_matrix_2nd_read_strand.tsv', header =T, check.names = F, stringsAsFactors = F))
sunshine.cts <- as.data.frame(fread('/home/data/WashU_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_count_matrix_2nd_read_strand.tsv', header =T, check.names = F, stringsAsFactors = F))
mend.tpm <- as.data.frame(fread('/home/data/WashU_Data/bulkRNASeq/201703_MendelianVsSporadics/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies/gene_quant_matrix_TPM.tsv',header=T, check.names=F, stringsAsFactors = F))
sunshine.tpm <- as.data.frame(fread('/home/data/WashU_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies/gene_quant_matrix_TPM.tsv',header=T, check.names=F, stringsAsFactors = F))
mend.fpkm <- as.data.frame(fread('/home/data/WashU_Data/bulkRNASeq/201703_MendelianVsSporadics/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_fpkm_matrix_2nd_read_strand.tsv',header=T, check.names=F, stringsAsFactors = F))
sunshine.fpkm <- as.data.frame(fread('/home/data/WashU_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_fpkm_matrix_2nd_read_strand.tsv', header =T, check.names = F, stringsAsFactors = F))
meta.cols <- sunshine.cts[, c('GeneID', 'GeneName', 'GeneBiotype')]

## read technical data 
#mend.tech <- read.csv('/40/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/03.-Phenotype/2020_09_11_WashU_MendalianVsSporadics_technical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F)
mend.tech <- read.csv('/home/data/WashU_Data//bulkRNASeq/201703_MendelianVsSporadics/03.-Phenotype/2020_09_11_WashU_MendalianVsSporadics_technical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F)
#sunshine.tech <- read.csv('/40/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/03.-Phenotype/2020_09_16_WashU_SUNSHINE_technical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F)
sunshine.tech <- read.csv('/home/data/WashU_Data//bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/03.-Phenotype/2020_09_16_WashU_SUNSHINE_technical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F)

##################################### read exp data ##################################### 
rnaseq <- as.data.frame(fread('omics_integration/data/brain_exp_matrix_FPKM_2nd_read_strand.tsv', header =T, stringsAsFactors = F, check.names = F))
#rnaseq <- read.table('omics_integration/data/brain_exp_matrix_Salmon.tsv', header =T, stringsAsFactors = F, check.names = F)
mean.expr = rowMeans(rnaseq[, -c(1,2)], na.rm = T)
rnaseq = rnaseq[order(mean.expr, decreasing=T),]
rnaseq = rnaseq[!duplicated(rnaseq[["GeneName"]]),]
rownames(rnaseq) <- rnaseq$GeneName
rnaseq$GeneID <- NULL
rnaseq$GeneName <- NULL
## filter lowly expressed genes 
#keep <- rowSums(rnaseq <= 0.2) >= floor(length(colnames(rnaseq))*0.75)
#rnaseq <- rnaseq[!keep, ]

########### after batch correction expression 
# bc.rnaseq <- read.table('batch_correction/TPM/HK_genes/rnaseq_after_correction_nromalized_counts_k2.tsv', header =T, stringsAsFactors = F, check.names = F)
# bc.rnaseq$GeneID = rownames(bc.rnaseq)
# bc.rnaseq <- merge (bc.rnaseq, rnaseq[, c('GeneID', 'GeneName')])
# rownames(bc.rnaseq) <- bc.rnaseq$GeneName
# bc.rnaseq$GeneID <- NULL
# bc.rnaseq$GeneName <- NULL
# s.cols <- colnames(bc.rnaseq)[grepl("HW", colnames(bc.rnaseq))]
# m.cols <- colnames(bc.rnaseq)[!colnames(bc.rnaseq) %in% s.cols]
# mend.tech <- read.csv('/40/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/03.-Phenotype/2020_09_11_WashU_MendalianVsSporadics_technical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F)
# sunshine.tech <- read.csv('/40/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/03.-Phenotype/2020_09_16_WashU_SUNSHINE_technical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F)
# for (i in 1:ncol(bc.rnaseq)) {
#   col.name <- colnames(bc.rnaseq)[i]
#   if (col.name %in% m.cols) {
#     subj.id <- mend.tech[mend.tech$Sample_ID==col.name, 'Subj_ID']
#     colnames(bc.rnaseq)[colnames(bc.rnaseq)==col.name] <- subj.id
#   } else {
#     subj.id <- sunshine.tech[sunshine.tech$Sample_Name==col.name, 'Subj_ID']
#     colnames(bc.rnaseq)[colnames(bc.rnaseq)==col.name] <- subj.id
#   }
# }
# colnames(bc.rnaseq) <- gsub("MAP_", "ID_", colnames(bc.rnaseq))
# colnames(bc.rnaseq) <- gsub("MAC_", "ID_", colnames(bc.rnaseq))
# colnames(bc.rnaseq) <- gsub("DIAN_", "ID_", colnames(bc.rnaseq))
# r.data.cols <- colnames(bc.rnaseq) 
# bc.rnaseq = bc.rnaseq[apply(bc.rnaseq, 1, var) > 0,]
# rnaseq = bc.rnaseq
#############################################################################################

################################ read proteomics data ######################################
# proteomics clinical data
prt.phen.data <-  read.csv('omics_integration/data/Proteomics/v4_brain_afterQC_phenoFile-2020.csv', header =T, check.names = F, stringsAsFactors = F)
#prt.data <- read.table('omics_integration/data/brain_protein_matrix_v2.tsv', header =T, sep="\t", stringsAsFactors = F, check.names = F)
# ------------------------------------------------------------------------------------
## combat corrected file 
prt.data <- readRDS('/home/data/WashU_Data/Metabolomics_Data/02.-Processed_Data/201910_Metabolon_BrainTissue_CULEBRA/04-kizer_analysis/multi-omics/data/model_data_ad_combat.rds')
prt.col.data <- colnames(prt.data)[!colnames(prt.data) %in% c('ExtIdentifier', 'MAPID', 'Status', 'age_at_death', 'sex')]
prt.data <- as.data.frame(t(prt.data[, prt.col.data]))
#prt.data <- as.data.frame(t(readRDS('omics_integration/data/Proteomics/02-ComBat_433x1092.rds')))
prt.data$SOMAseqID <- rownames(prt.data)
org.prt.data <- read.csv('omics_integration/data/Proteomics/20200512_brain_afterQC_proteomic_exprs.csv', header =T, check.names = F, stringsAsFactors = F)[,c(1,5,8)]
prt.data <- merge(prt.data, org.prt.data)
prt.data$SOMAseqID <- NULL

# SampleId MAPID ExtIdentifier
exid <- colnames(prt.data)[grepl("EXID", colnames(prt.data))]
for (i in 1:length(exid)) {
  id <- exid[i]
  ## extract the sample id & mapid
  map_id <- prt.phen.data[prt.phen.data$ExtIdentifier == id, 'MAPID']
  if (is.na(map_id)) { stop() }
  colnames(prt.data)[colnames(prt.data) == id] <- paste0('ID_', map_id)
}
# ------------------------------------------------------------------------------------

p.data.cols <- colnames(prt.data)[grepl('ID_', colnames(prt.data))]
mean.expr = rowMeans(prt.data[, p.data.cols], na.rm = T)
prt.data = prt.data[order(mean.expr, decreasing=T),]
prt.data = prt.data[!duplicated(prt.data[["EntrezGeneSymbol"]]),]
rownames(prt.data) <- prt.data$EntrezGeneSymbol
prt.data$EntrezGeneID <- NULL
prt.data$EntrezGeneSymbol <- NULL
prt.data$Target <- NULL

## handle missing values by excluding the ones with NA in > 20%
prt.data$na_count <- apply(prt.data, 1, function(x) sum(is.na(x)))
prt.data <- prt.data[prt.data$na_count < length(colnames(prt.data)[colnames(prt.data) !='na_count']) * 0.20, ]
prt.data$na_count <- NULL
## replace all NA with the mean value
for(i in 1:nrow(prt.data)){
  indx = which(is.na(prt.data[i,]))
  row.min <- rowMeans(prt.data[i,], na.rm = T)
  prt.data[i, indx] <- row.min
}

## fix missing genes 
missing.genes <- read.table('omics_integration/data/missung_genes_prt.tsv', header =T, stringsAsFactors = F, check.names = F, sep="\t")
for (k in 1:nrow(missing.genes[1:22, ])) {
  gg <- missing.genes$org.name[k]
  aliase.g <- missing.genes$aliase.name[k]
  rownames(prt.data)[rownames(prt.data) ==gg] = aliase.g
}

#############################################################################################
## read metabolomics data 
metab <- read.table('omics_integration/data/QCd_metabolites_for_RNAseq_brains.csv', header =T, sep=',', stringsAsFactors = F, check.names = F)
metab.data.cols <- colnames(metab)[14:ncol(metab)]
metab <- metab[, c('COMP.ID', metab.data.cols)]
metab$COMP.ID <- paste0('COMPID_', metab$COMP.ID)
rownames(metab) <- metab$COMP.ID
metab$COMP.ID <- NULL
colnames(metab) <- paste0('ID_', colnames(metab))
## check the number of missing values - exclude the ones with >=20%
metab$na_count <- apply(metab, 1, function(x) sum(is.na(x)))
metab <- metab[metab$na_count < length(metab.data.cols) * 0.20, ]
metab$na_count <- NULL
## replace all NA with the lowest value
for(i in 1:nrow(metab)){
  indx = which(is.na(metab[i,]))
  r.min <- rowMeans(metab[i,], na.rm = T)
  metab[i, indx] <- r.min
}

##############################################################################
### read methylation data
# meth <- as.data.frame(fread('omics_integration/data//beta_values_methylation_abdallah.csv'))
# colnames(meth)[colnames(meth)=="V1"] = "CG_ID"
# meth.key <- read.table('omics_integration/data/phenotype_key.txt', header =T, sep="\t", stringsAsFactors = F, check.names = F)
# meth.cols = as.data.frame(colnames(meth)[-1])
# colnames(meth.cols) = "Sample_Name"
# meth.cols = merge(meth.cols, meth.key, sort =F)
# colnames(meth) = c('CG_ID', as.character(meth.cols$UniquePhenoID))
# rownames(meth) <- meth$CG_ID
# meth$CG_ID <- NULL
# colnames(meth) <- gsub("MAP", "ID", colnames(meth))
# colnames(meth) <- gsub("DIAN", "ID", colnames(meth))
##############################################################################

## extract shared samples 
shared.samples <- intersect(colnames(rnaseq), colnames(prt.data))
shared.samples <- intersect(shared.samples, colnames(metab))
#shared.samples <- intersect(shared.samples, colnames(meth))
## extract shared genes
shared.genes <- intersect(rownames(rnaseq), rownames(prt.data))

## extract shared samples 
rnaseq <- rnaseq[, shared.samples]
prt.data <- prt.data[, shared.samples]
metab <- metab[, shared.samples]
#meth <- meth[, shared.samples]

## extract shared genes data 
#rnaseq <- rnaseq[rownames(rnaseq) %in% shared.genes, ]
#prt.data <- prt.data[rownames(prt.data) %in% shared.genes, ]

## check all samples in the same order 
all(colnames(rnaseq)==colnames(prt.data))
all(colnames(prt.data)==colnames(metab))
#all(colnames(metab)==colnames(meth))

## make a list of omic datasets
omics.data <- list(rnaseq=rnaseq, proteomic= prt.data, metab=metab)
## check the dimensions of all datasets 
sapply(omics.data, dim)

## remove failed/flagged samples 
flagged.samples <- c(mend.tech[mend.tech$QC_Stats=='Fail', 'Subj_ID'], sunshine.tech[sunshine.tech$FLAG=='Fail', 'Subj_ID'])
#flagged.samples <- read.table('WashU_QC_flagged_failed_samples.tsv', header =T, stringsAsFactors = F, check.names = F) 
flagged.samples <- gsub('MAP_', 'ID_', flagged.samples)
omics.data <- lapply(omics.data, function(x) x[, !colnames(x) %in% flagged.samples])
sapply(omics.data, dim)

# filter features with no variance at all
for (i in 1:length(omics.data)) {
  omics.data[[i]] = omics.data[[i]][apply(omics.data[[i]], 1, var) > 0,]
}
sapply(omics.data, dim)

## select only AD and CO
co.ids <- paste0("ID_", clin.data$Subj_ID[clin.data$Status %in% c('Neuro_CO')])
ad.ids <- paste0("ID_", clin.data$Subj_ID[clin.data$Status %in% c('Neuro_AD')])
adad.ids <- paste0("ID_", clin.data$Subj_ID[clin.data$Status %in% c('ADAD_carrier', 'ADAD_Carrier')])
pd.dlb.ids <- paste0("ID_", clin.data$Subj_ID[clin.data$Status %in% c('Neuro_PD', 'Neuro_DLB')]) # 'Neuro_DLB', 'Neuro_PD'
prs.ids <- paste0("ID_", clin.data$Subj_ID[clin.data$Status %in% c('Neuro_Presymptomatic')])
trem2.ids <- paste0("ID_", clin.data$Subj_ID[!is.na(clin.data$TREM2_all_variants)])

omics.data <- lapply(omics.data, function(x) x[, colnames(x) %in% c(ad.ids)])
sapply(omics.data, dim)

## select top variant features for all omics data
keep.high.var.features <- function(omic, num.features=2000) {
  if (nrow(omic) < num.features) {
    return(omic)
  } else {
    feature.vars = apply(omic, 1, var)
    threshold = feature.vars[order(feature.vars, decreasing = T)][num.features]
    #return(omic[feature.vars >= threshold | rownames(omic) %in% shared.genes, ]) 
    return(omic[feature.vars >= threshold, ]) 
  }
}
### filter features 
omics.data = lapply(omics.data, keep.high.var.features)
## check the dimensions of all datasets after filtering 
sapply(omics.data, dim)


### plot rnaseq expression before integration 
exp.data <- log2(omics.data[["rnaseq"]]+1)
colBreaks = seq(-5, 5, 0.5)
#col.pan <- colorpanel(length(colBreaks)-1, "blue", "white", "red")
#exp.data[exp.data > 2.5] = 2.5
#exp.data[exp.data < -2.5] = -2.5
pdf('omics_integration/iCluster_output/Figures/rnaseq_exp_heatmap.pdf', width = 15, height = 15)
par(cex.main=1.8)
heatmap.2(as.matrix(exp.data), breaks = colBreaks, 
          col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(length(colBreaks)-1)),
          Rowv = T, Colv = T, labRow = "", labCol = "", scale="row", 
          #key.par=list(mar=c(0,5,20,10)), 
          trace="none", dendrogram="both", main= "RNA-Seq Expression",
          cexRow=1, cexCol=1.4, density.info="none", 
          key.title="", key.xlab="Expression (FPKM)",
          margin=c(10,10), lhei=c(2,20), lwid=c(2,5))
dev.off()

## plot protein expression before integration   
prt.exp <- log2(omics.data[["proteomic"]]+1)
colBreaks = seq(-2, 2, 0.1)
#col.pan <- colorpanel(length(colBreaks)-1, "blue", "white", "red")
pdf('omics_integration/iCluster_output/Figures/protein_exp_heatmap.pdf', width = 15, height = 15)
par(cex.main=1.8)
heatmap.2(as.matrix(prt.exp), breaks = colBreaks, 
          col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(length(colBreaks)-1)),
          Rowv = T, Colv = T, labRow = "", labCol = "", scale="row", 
          #key.par=list(mar=c(0,5,20,10)), 
          trace="none", dendrogram="both", main= "Protein Expression",
          cexRow=1, cexCol=1.4, density.info="none", 
          key.title="", key.xlab="Expression (FPKM)",
          margin=c(10,10), lhei=c(2,20), lwid=c(2,5))
dev.off()

## plot metabolomics expression before integration      
metab.readings <- omics.data[["metab"]]
#metab.readings <-t(log2(omics.data[[3]]))
colBreaks = seq(-5, 5, 0.5)
#colBreaks = seq(-10, 10, 1)
#col.pan <- colorpanel(length(colBreaks)-1, "blue", "white", "red")
pdf('omics_integration/iCluster_output/Figures/metab_readings_heatmap.pdf', width = 15, height = 15)
par(cex.main=1.8)
heatmap.2(as.matrix(metab.readings), breaks = colBreaks, 
          col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(length(colBreaks)-1)),
          Rowv = T, Colv = T, labRow = "", labCol = "", scale="row", 
          #key.par=list(mar=c(0,5,20,10)), 
          trace="none", dendrogram="both", main= "Metabolomics Readings",
          cexRow=1, cexCol=1.4, density.info="none", 
          key.title="", key.xlab="Expression (FPKM)",
          margin=c(10,10), lhei=c(2,20), lwid=c(2,5))
dev.off()

################################################################################
############## run/tune iClusterBayes     
################################################################################
set.seed(5152) 
MAX.NUM.CLUSTERS = 11
dev.ratios = c()
icluster.rets = list()
num.omics <- length(omics.data)
icluster.ret = tune.iClusterBayes(cpus=(MAX.NUM.CLUSTERS - 1), 
                                    t(log2(omics.data[["rnaseq"]]+1)), 
                                    t(10^omics.data[["proteomic"]]), 
                                    t(omics.data[["metab"]]), 
                                    #t(omics.data[["meth"]]), 
                                    K=1:(MAX.NUM.CLUSTERS - 1), 
                                    type=rep('gaussian', num.omics),
                                    n.burnin=12000, n.draw=22000,
                                    prior.gamma=rep(0.1, num.omics),
                                    pp.cutoff = 0.5,
                                    sdev=0.5,
                                    thin=1
                                  )$fit

# save the result object
saveRDS(icluster.ret, file="omics_integration/iCluster_output/icluster.res_seed_5151_shared_genes_only.rds")
################################################################################

################################################
###### extract dev.ratio's & BIC
############################################### 
dev.ratios = lapply(1:(MAX.NUM.CLUSTERS - 1), function(i) icluster.ret[[i]]$dev.ratio)
allBICs = lapply(1:(MAX.NUM.CLUSTERS - 1), function(i) icluster.ret[[i]]$BIC)

## extract best cluster 
optimal.solution = icluster.ret[[which.max(dev.ratios)]] 
best.clusters = optimal.solution$clusters
#save(optimal.solution, file="omics_integration/best.solution.rds")

## plot clusters vs dev.ratio 
dd1 = data.frame(k=1:(MAX.NUM.CLUSTERS - 1), dev.ratio= unlist(dev.ratios))
dd2 = data.frame(k=1:(MAX.NUM.CLUSTERS - 1), bic= unlist(allBICs))

p1 = ggplot(dd1, aes(x=k, y=dev.ratio)) + geom_line(color="orange3",  lwd=1) + 
  geom_point(color="orange3", size=3) + theme_bw() +
  labs(x="Number of clusters",  y="Deviance ratio") +
  ggtitle("Clusters vs. Deviance ratio") + 
  scale_x_continuous(breaks = 1:(MAX.NUM.CLUSTERS - 1) ) + 
  geom_vline(xintercept = which.max(dev.ratios), color="red", linetype="dashed") + 
  theme (axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),
         plot.title = element_text(size=14, hjust=0.5, face="bold"))

p2 = ggplot(dd2, aes(x=k, y=bic)) + geom_line(color="orange3",  lwd=1) + 
  geom_point(color="orange3", size=3) + theme_bw() +
  labs(x="Number of clusters",  y="BIC") +
  ggtitle("Clusters vs. Bayesian information criteria (BIC)") + 
  scale_x_continuous(breaks = 1:(MAX.NUM.CLUSTERS - 1) ) + 
  geom_vline(xintercept = which.min(allBICs), color="red", linetype="dashed") + 
  theme (axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),
         plot.title = element_text(size=14, hjust=0.5, face="bold")) 


pdf('omics_integration/Figures/k_vs_devRatio.pdf', width = 12, height = 5)
grid.arrange(p1, p2, ncol=2)
dev.off()


### plot the posterior probabity 
pdf('omics_integration/Figures/posterior_probability.pdf', width = 12, height = 10)
par(mfrow=c(3,1))
k=which.max(dev.ratios)
plot(icluster.ret[[k]]$beta.pp[[1]], xlab="Genes", ylab="Posterior probability", main="RNA-Seq expression")
plot(icluster.ret[[k]]$beta.pp[[2]], xlab="Genes", ylab="Posterior probability", main="Protein level")
plot(icluster.ret[[k]]$beta.pp[[3]], xlab="Genes", ylab="Posterior probability", main="Metabolomics")
#plot(icluster.ret[[k]]$beta.pp[[4]], xlab="Genes", ylab="Posterior probability", main="Methylation")

dev.off()

###########################################################################
################### feature selection ################### 
###########################################################################
features = alist()
features[["rnaseq"]] = colnames(t(omics.data[["rnaseq"]]))
features[["proteomic"]] = colnames(t(omics.data[["proteomic"]]))
features[["metab"]] = colnames(t(omics.data[["metab"]]))

sigfeatures=alist()
for(i in 1:num.omics){
  rowsum=apply(abs(optimal.solution$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]] = (features[[i]])[which(rowsum>upper)]
}
names(sigfeatures)=c("expression","proteomics", "metab")
#print a few examples of selected features
head(sigfeatures[[1]])
#write.table(sigfeatures$expression, file ="omics_integration/sigfeatures_expression.txt", quote = F, row.names = F, col.names = F )
#write.table(sigfeatures$proteomics, file ="omics_integration/sigfeatures_proteomics.txt", quote = F, row.names = F, col.names = F )

### plot the heatmap 
exp.image=scale(log2(t(omics.data[["rnaseq"]])+1))
exp.image[exp.image > 2.5] = 1.5
exp.image[exp.image < -2.5] = -1.5

prt.image = scale(t(omics.data[["proteomic"]]))
prt.image[prt.image > 2.5 ] = 2.5
prt.image[prt.image< -2.5 ] = --2.5

metab.image = scale(t(omics.data[["metab"]]))
metab.image[metab.image > 2.5 ] = 2.5
metab.image[metab.image< -2.5 ] = --2.5

#meth.image = scale(t(omics.data[["meth"]])+1)
#meth.image[meth.image > 2.5 ] = 2.5
#meth.image[meth.image< -2.5 ] = --2.5

col.scheme = alist()
col.scheme[[1]] = rev(colorRampPalette(brewer.pal(10, "RdBu"))(148))
col.scheme[[2]] = rev(colorRampPalette(brewer.pal(10, "RdBu"))(148))
col.scheme[[3]] = rev(colorRampPalette(brewer.pal(10, "RdBu"))(148))
#col.scheme[[4]] = rev(colorRampPalette(brewer.pal(10, "RdBu"))(148))

pdf('omics_integration/iCluster_output/Figures/iClusterPlus_Heatmap_high_nBurnin_nDraw_stdev_0.5_Sep21.pdf', width = 12, height = 8)
plotHMBayes(fit=optimal.solution, 
             datasets=list(exp.image, prt.image, metab.image), 
             type=rep("gaussian",num.omics),
             scale = rep(F,num.omics), 
             #col.scheme = col.scheme, 
             #sample.order <- c(c1.ss, c2.ss, c3.ss, c4.ss),
             threshold=c(0.5,0.5, 0.009),
             row.order=rep(T,num.omics),  
             sparse=rep(T,num.omics),
             cap=rep(T,num.omics)
           )
dev.off()

## write features contributed to the cluster significantly 
for(i in 1:num.omics){ 
  ff.res <- as.data.frame(optimal.solution$beta.pp[[i]])
  if (i==1) {
    ff.res$feature.name <- rownames(omics.data[["rnaseq"]])
    colnames(ff.res) <- c('pp', 'feature.name')
    write.table(ff.res$feature.name[ff.res$pp > 0.5], file='omics_integration/iCluster_output/rnaseq_sig_gene.tsv', quote = F, row.names = F, sep="\t", col.names = F)
  } else if (i==2){
    ff.res$feature.name <- rownames(omics.data[["proteomic"]])
    colnames(ff.res) <- c('pp', 'feature.name')
    write.table(ff.res$feature.name[ff.res$pp > 0.5], file='omics_integration/iCluster_output/protein_sig_gene.tsv', quote = F, row.names = F, sep="\t", col.names = F)
  } else if (i==3) {
    ff.res$feature.name <- rownames(omics.data[["metab"]])
    colnames(ff.res) <- c('pp', 'feature.name')
    write.table(ff.res$feature.name[ff.res$pp > 0.009], file='omics_integration/iCluster_output/metab_sig_features.tsv', quote = F, row.names = F, sep="\t", col.names = F)
  }
}

# plotHeatmap(fit=optimal.solution,
#             datasets=list(exp.image, prt.image), 
#             type=c("gaussian","gaussian"), col.scheme = col.scheme,
#             row.order=c(T,T), sparse=c(T,T),cap=c(F,F))


## extract cluster membership of the best class 
all.clusters <- matrix(data= NA, nrow= nrow(t(omics.data[["rnaseq"]])), ncol=MAX.NUM.CLUSTERS - 1)
for (k in 1:(MAX.NUM.CLUSTERS - 1)) {
  cc = icluster.ret[[k]]$clusters
  cc = as.matrix(cc)
  all.clusters[, k] = cc
}

rownames(all.clusters) = colnames(omics.data[["rnaseq"]])
colnames(all.clusters) = paste0('K', 1:(MAX.NUM.CLUSTERS - 1))
all.clusters = as.data.frame(all.clusters)
all.clusters$Subj_ID = rownames(all.clusters)
rownames(all.clusters) = NULL
all.clusters = all.clusters[, c('Subj_ID',  colnames(all.clusters)[colnames(all.clusters) !="Subj_ID"])]
## get best cluster membership
best.cl= unique(c(which.max(dev.ratios), which.min(allBICs)))
if (length(best.cl) > 1) {
  stop('Error: something is wrong!!')
}

best.cluster.memership <- all.clusters[, c('Subj_ID',  paste0('K', best.cl))]
colnames(best.cluster.memership) <- c('Subj_ID', 'best.cluster')
best.cluster.memership <- best.cluster.memership[order(best.cluster.memership$best.cluster), ]
### will be deleted 
table(best.cluster.memership$best.cluster[best.cluster.memership$Subj_ID %in% c4])
#write best result solution membership samples 
#write.table(best.cluster.memership, 'omics_integration/best.cluster.memership.tsv', sep="\t", quote = F, row.names = F)


###################################################################
################ check the phenotype association #################
###################################################################
bb <- read.table("omics_integration/iCluster_output/best.cluster.memership_v1_full_AD_only.tsv", header = T, sep="\t", stringsAsFactors = F)
## remove flagged samples 
bb <- bb[!bb$Subj_ID %in% as.integer(gsub("ID_","", flagged.samples)), ]

## calculate age differences 
aao.pval <- format(t.test(bb$AAO[bb$best.cluster== 4], bb$AAO[bb$best.cluster != 4])$p.value, scientific = T)
aao.diff <- abs(mean(bb$AAO[bb$best.cluster == 4]) - mean(bb$AAO[bb$best.cluster != 4], na.rm = T))
aod.pval <- format(t.test(bb$AOD[bb$best.cluster== 4], bb$AOD[bb$best.cluster != 4])$p.value, scientific = T)
aod.diff <- abs(mean(bb$AOD[bb$best.cluster == 4]) - mean(bb$AOD[bb$best.cluster != 4], na.rm = T))

#best.cluster.memership = read.table("omics_integration/iCluster_output/best.cluster.memership_v1.tsv", header = T, stringsAsFactors = F)
# best.cluster.memership$Subj_ID = gsub('ID_', '', best.cluster.memership$Subj_ID)
# 
# #################### merge with clinical #################
# bb = merge (best.cluster.memership, clin.data)
# bb$duration = bb$AOD - bb$AAO
# bb[!is.na(bb$TREM2_all_variants), 'TREM2_all_variants'] = "Mutated"
# bb[is.na(bb$TREM2_all_variants), 'TREM2_all_variants'] = "Not mutated"
# bb[!is.na(bb$PLD3_all_variants), 'PLD3_all_variants'] = "Mutated"
# bb[is.na(bb$PLD3_all_variants), 'PLD3_all_variants'] = "Not mutated"
# bb[!is.na(bb$MS4A4A_1), 'MS4A4A_1'] = "Mutated"
# bb[is.na(bb$MS4A4A_1), 'MS4A4A_1'] = "Not mutated"
# bb[!is.na(bb$PSEN1), 'PSEN1'] = "Mutated"
# bb[is.na(bb$PSEN1), 'PSEN1'] = "Not mutated"
# bb[!is.na(bb$PSEN2), 'PSEN2'] = "Mutated"
# bb[is.na(bb$PSEN2), 'PSEN2'] = "Not mutated"
# bb$TREM2.GRP = "mutated"
# bb[is.na(bb$TREM2), 'TREM2.GRP'] = "Not mutated"

## plot the scores
#bc = melt(best.cl.mem[,c('Subj_ID', 'best.cluster', 'BraakTau')], id.vars = 'Subj_ID')
#bb1= bb[bb$Status %in% c('Neuro_AD', 'Neuro_CO', 'Neuro_DLB', 'Neuro_PD'), ] # 'Neuro_AD_DLB', 'Neuro_AD_PD'

## compute pvalue 
pvals.df <- NULL
aod_df <- NULL
for (cc in c('1vs2', '1vs3', '1vs4', '2vs3', '2vs4', '3vs4', '123vs4')) {
  g1 <- as.numeric(unlist(strsplit(cc, "vs"))[1])
  g2 <- as.numeric(unlist(strsplit(cc, "vs"))[2])
  if (cc =='123vs4') {
    ff <- bb
    ff$grp <- 'C4'
    ff[ff$best.cluster %in% c(1,2,3), 'grp'] <- 'OT'
    glm.res <- glm(CDRe ~ grp + Sex + AOD, data = ff, family = 'gaussian')
    pval <- data.frame(comp=cc, PValue= signif(coef(summary(glm.res))[,4][["grpOT"]], digits = 3), y_pos = max(ff$CDRe, na.rm= T), best.cluster = NA)
  } else {
    ff <- bb[bb$best.cluster %in% c(g1, g2), ] 
    glm.res <- glm(CDRe ~ best.cluster + Sex + AOD, data = ff, family = 'gaussian')
    pval <- data.frame(comp=cc, PValue= signif(coef(summary(glm.res))[,4][["best.cluster"]], digits = 3), y_pos = max(ff$CDRe, na.rm= T), best.cluster = NA)
  }
  pvals.df <- rbind(pvals.df, pval)
}
#pvals.df$PValue <- as.numeric(format(pvals.df$PValue , scientific = T))
comp.list <- strsplit(as.character(pvals.df$comp[pvals.df$PValue < 0.05]), "vs")
#pvals.list <- paste0('p=',format(pvals.df$PValue, scientific = T))
pvals.df$PValue <-  paste0('p=',format(pvals.df$PValue, scientific = T))

##### for CDR 
tmp <- aggregate(Subj_ID ~ best.cluster + CDRe, bb, length)
tmp2 <- as.data.frame(table(bb$best.cluster))
colnames(tmp2) <- c('best.cluster', 'Num.Cases')
tmp <- merge(tmp, tmp2, by ="best.cluster")
tmp$pct <- tmp$Subj_ID/tmp$Num.Cases
mean.dist <- (sum(unique(tmp$Subj_ID[tmp$CDRe==3.0])))/sum(table(bb$best.cluster))
######

pheno.name <- 'CDRe'
p1 <- ggplot(bb, aes(x=as.factor(best.cluster), y= as.factor(CDRe), color=as.factor(best.cluster))) + 
             geom_point(size = 2,  shape=21, position=position_jitter(width=0.2, height=0.1)) + theme_bw() + 
            #geom_jitter(width = 0.6, size=1.5, aes(color = factor(best.cluster))) + theme_bw() +      
            #geom_boxplot(aes(color = factor(best.cluster)), outlier.shape=NA, outlier.size=0.5) +
            #geom_jitter(position=position_jitter(0.3), size=1.5, aes(color = factor(best.cluster))) + theme_bw() +
            labs(x='', y='CDRe Score') + ggtitle('') + scale_color_brewer(palette = "Dark2") +  
            theme(axis.text.x=element_text(size=18, vjust=0.5, color="black", angle=90),
                 axis.text.y=element_text(size=18, color="black"),
                 axis.title.x=element_text(size=19, face="bold"),
                 axis.title.y=element_text(size=18, face="bold"),
                 #panel.background=element_rect(color="black"),
                 plot.title = element_text(size = 14, hjust=0.5, color="black", face="bold"), legend.position="none") + 
  
            #scale_x_discrete(labels=c("1"="Cluster1", "2"="Cluster2", "3"="Cluster3", "4"="Cluster4")) +  
            #scale_y_continuous(breaks = c(0.5,1,2,3)) + 
      #scale_x_discrete(labels=c("1"=paste0("Cluster1\n(n=",length(bb$Subj_ID[bb$best.cluster==1]),')'), "2"= paste0("Cluster2\n(n=",length(bb$Subj_ID[bb$best.cluster==2]),')'),
      #                               "3"=paste0("Cluster3\n(n=",length(bb$Subj_ID[bb$best.cluster==3]),')'), "4"= paste0("Cluster4\n(n=",length(bb$Subj_ID[bb$best.cluster==4]),')'))) +
      geom_signif(data = pvals.df, aes(xmin = 1, xmax = 4, annotations = PValue[comp=="1vs4"], y_position = y_pos[comp=="2vs4"]+0.60), 
              textsize = 5.5, size=0.6, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") +
      geom_signif(data = pvals.df, aes(xmin = 2, xmax = 4, annotations = PValue[comp=="2vs4"], y_position = y_pos[comp=="2vs4"]+1.19), 
            textsize = 5.5, size=0.6, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") +
     geom_signif(data = pvals.df, aes(xmin = 3, xmax = 4, annotations = PValue[comp=="3vs4"], y_position = y_pos[comp=="3vs4"]+1.75), 
            textsize = 5.5, size=0.6, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") 
      ## for BraakTau
      # geom_signif(data = pvals.df, aes(xmin = 1, xmax = 3, annotations ='', y_position = pvals.df$y_pos[pvals.df$comp=="123vs4"]-0.2), 
      #         textsize = 5.5, size=0.6, step_increase=0.1, vjust = -0.2, manual = T, fontface="bold", lty="dashed") +
      # geom_signif(data = pvals.df, aes(xmin = 2, xmax = 4, annotations = paste0('p=',pvals.df$PValue[pvals.df$comp=="123vs4"]), y_position = pvals.df$y_pos[pvals.df$comp=="123vs4"]+0.1), 
      #         textsize = 5.5, size=0.6, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") 
p1 + stat_pvalue_manual(pvals.df, label = "PValue", y.position = 3.1)

pdf(paste0('omics_integration/iCluster_output/Figures/Final/phenotype_assoication_',pheno.name,'_Final_V2.pdf'), width = 4, height = 4, useDingbats = F)
p1
dev.off()


###############################################################
################ Plot metablomics profiles ####################
###############################################################
# read PC1 
bb <- read.table("omics_integration/iCluster_output/best.cluster.memership_v1_full.tsv", header = T, sep="\t", stringsAsFactors = F)
## remove flagged samples 
bb <- bb[!bb$Subj_ID %in% as.integer(gsub("ID_","", flagged.samples)), ]
bb[bb$Status=="Neuro_CO", 'best.cluster'] <- 0

#pc1 <- read.csv('omics_integration/data/17metab_PC1_AD_ADAD_TREM2_CO.csv', header =T, sep=",", check.names = F)
pc1 <- read.csv('omics_integration/data/Metabolomics/16metab_PC1_AD_ADAD_TREM2_CO.csv', header =T, sep=",", check.names = F)
pc1$Subj_ID <- gsub("MAP_", "", pc1$Subj_ID)
pc1$Subj_ID <- gsub("DIAN_", "", pc1$Subj_ID)
colnames(pc1) <- c('Subj_ID', 'Metab.PC1')
bb <- merge(bb, pc1, all.x= T)
## relace missing values by min controls 
bb[is.na(bb$Metab.PC1), 'Metab.PC1'] <- min(bb$Metab.PC1[bb$best.cluster==0], na.rm = T)

# ## add metab readings 
# mm = metab
# mm = as.data.frame(colMeans(mm))
# mm$ID = rownames(mm)
# rownames(mm) <- NULL
# colnames(mm) <- c('metab.readings', 'Subj_ID')
# mm$Subj_ID <- gsub("ID_", "", mm$Subj_ID)
# bb <- merge(bb, mm, sort =F)

#cmp.list = as.list(as.data.frame(combn(1: length(unique(bb$best.cluster)),2)))
# m <- ggplot(bb, aes(x=factor(best.cluster), y=metab.readings, group=factor(best.cluster))) + geom_boxplot(outlier.shape=NA, outlier.size=0.5) 
# m <- m + geom_jitter(position=position_jitter(0.3),  size=2) + theme_bw() + guides(color = FALSE) 
# m <- m + labs(x='', y='Readings') + ggtitle('Distribution of the means of metablomities by cluster')
# #p1 <- p1 + geom_text(label = ifelse(bb$best.cluster==2,bb$Subj_ID,""), size=4, position = position_jitter(seed = 0.5))
# m <- m + theme(axis.text.x=element_text(size=11, vjust=0.5, color="black"),
#                axis.text.y=element_text(size=12, color="black"), 
#                axis.title.y=element_text(size=12),
#                #panel.background=element_rect(color="black"),
#                plot.title = element_text(size = 14, hjust=0.5, color="black", face="bold"),
#                legend.position="right", 
#                panel.border = element_rect(linetype='solid', color='black'),
#                plot.margin=unit(c(1,1,1,5), 'mm'))
# #p1 <- p1 + scale_color_manual(values=c("#E69F00", "#56B4E9")) 
# m <- m + scale_x_discrete(labels=c("0"=paste0("Control\n(n=",length(bb$Subj_ID[bb$best.cluster==0]),')'), 
#                                    "1"=paste0("Cluster1\n(n=",length(bb$Subj_ID[bb$best.cluster==1]),')'), 
#                                    "2"= paste0("Cluster2\n(n=",length(bb$Subj_ID[bb$best.cluster==2]),')'),
#                                    "3"=paste0("Cluster3\n(n=",length(bb$Subj_ID[bb$best.cluster==3]),')'), 
#                                    "4"= paste0("Cluster4\n(n=",length(bb$Subj_ID[bb$best.cluster==4]),')')))
#m <- m + geom_signif(comparisons= cmp.list, step_increase=0.1, textsize = 4, map_signif_level=function(p)sprintf("p = %.2g", p), show.legend=F)
pvals.df <- NULL
aod_df <- NULL
for (cc in c('0vs1', '0vs2','0vs3', '0vs4','1vs2', '1vs3', '1vs4', '2vs3', '2vs4', '3vs4')) {
  g1 <- as.numeric(unlist(strsplit(cc, "vs"))[1])
  g2 <- as.numeric(unlist(strsplit(cc, "vs"))[2])
  ff <- bb[bb$best.cluster %in% c(g1, g2), ]
  glm.res <- glm(Metab.PC1 ~ best.cluster + Sex + AOD, data = ff, family = 'gaussian')
  pval <- data.frame(comp=cc, PValue= signif(coef(summary(glm.res))[,4][["best.cluster"]], digits = 3), y_pos = max(ff$Metab.PC1, na.rm= T)+1, best.cluster = NA)
  pvals.df <- rbind(pvals.df, pval)
}
comp.list <- strsplit(as.character(pvals.df$comp[pvals.df$PValue < 0.05]), "vs")
pvals.list <- paste0('p=',format(pvals.df$PValue[pvals.df$PValue < 0.05], scientific = T))
bb$best.cluster <- factor(bb$best.cluster, levels=c(0,3,1,2,4))
m2 <- ggplot(bb, aes(x=factor(best.cluster), y=Metab.PC1, group=factor(best.cluster))) + 
      geom_boxplot(aes(fill = factor(best.cluster)), outlier.shape=NA, outlier.size=0.5) + theme_bw() + 
      #geom_jitter(position=position_jitter(0.3), size=1.5, aes(fill = factor(best.cluster)))  
      labs(x='', y='PC1') + ggtitle('') +
      theme(axis.text.x=element_text(size=18, vjust=1, hjust=1, color="black", angle=45),
                 axis.text.y=element_text(size=18, color="black"), 
                 axis.title.y=element_text(size=24, face="bold"),
                 #panel.background=element_rect(color="black", size=1),
                 plot.title = element_text(size = 14, hjust=0.5, color="black", face="bold"),
                 legend.position="none", panel.border = element_rect(linetype='solid', color='black'),
                 plot.margin=unit(c(1,1,1,5), 'mm')) +
      scale_fill_manual(values=c("0"="gray60", "1"="#1B9E77", "2"="#D95F02", "3"="#7570B3", "4"="#E7298A")) +
      scale_x_discrete(labels=c("0"="Control", "1"="Knight-C1", "2"="Knight-C2", "3"="Knight-C3", "4"="Knight-C4")) + 
      # scale_x_discrete(labels=c("0"=paste0("Control\n(n=",length(bb$Subj_ID[bb$best.cluster==0]),')'), 
      #                              "1"=paste0("Cluster1\n(n=",length(bb$Subj_ID[bb$best.cluster==1]),')'), 
      #                              "2"= paste0("Cluster2\n(n=",length(bb$Subj_ID[bb$best.cluster==2]),')'),
      #                              "3"=paste0("Cluster3\n(n=",length(bb$Subj_ID[bb$best.cluster==3]),')'), 
      #                              "4"= paste0("Cluster4\n(n=",length(bb$Subj_ID[bb$best.cluster==4]),')'))) +
       # geom_signif(data=bb, aes(xmin = 2, xmax = 5, annotations =  'ordinal regression pval = 6.02e-03', y_position = max(bb$Metab.PC1)+0.5),   
       #        textsize = 4.5, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold", tip_length = 0)
       annotate(geom="text", x=3.4, y= 6.75, label="Ordinal regression (p=6.0e-03)", size=5.8, fontface="bold") + 
       geom_segment(aes(x = 2, y = 6.1, xend = 5, yend = 6.1), size = 1, arrow = arrow(length = unit(0.35, "cm")))

      #geom_signif(comparisons = comp.list, annotations = pvals.list, map_signif_level = TRUE, textsize=4, step_increase=0.1, fontface="bold")
m2

pdf('omics_integration/iCluster_output/Figures/Final/16metab_PC1_dist_per_cluster_ordinal_regression_v2.pdf', width = 6, height = 4.8, useDingbats = F)
m2
dev.off()

###############################################################
########################### Plot top genes ####################
###############################################################
## assign ids ti batches 
sunshine.ids <- sunshine.clin$Subj_ID[sunshine.clin$Status %in% c('Neuro_AD', 'Neuro_CO')]
mend.ids <- mend.clin$Subj_ID[mend.clin$Status %in% c('Neuro_AD', 'Neuro_CO')]
ov.ids <- intersect(mend.ids, sunshine.ids)
bb$batch = 'NA'
bb[bb$Subj_ID %in% mend.ids, 'batch'] <- 'M'
bb[bb$Subj_ID %in% sunshine.ids, 'batch'] <- 'S'
bb[bb$Subj_ID %in% ov.ids, 'batch'] <- 'OV'
zz = bb[, c('Subj_ID', 'batch')]
colnames(zz) = c('sample', 'batch')
zz$sample = paste0("ID_", zz$sample)
## prepare layout matrix 
mat = matrix(ncol=2, nrow=3)
mat[, 1] = c(1,3,5)
mat[, 2] = c(2,4,5)

all.genes <- list()
pdf('omics_integration/Figures/sig_genes_plots_v3.pdf', width=10, height = 12)
for (i in 1:length(unique(c(sigfeatures$expression, sigfeatures$proteomics)))) {
  gene = sigfeatures$expression[i]
  
  ## extract RNAseq expression
  g.exp = melt(rnaseq[rownames(rnaseq) == gene, gsub("ID_","", colnames(rnaseq)) %in% bb$Subj_ID])
  colnames(g.exp) = c('sample', 'rnaseq.exp')
  ## merge with batches 
  g.exp = merge(g.exp, zz, sort =F)
  
  g.prt = melt(prt.data[rownames(prt.data) == gene, gsub("ID_","", colnames(prt.data)) %in% bb$Subj_ID])
  colnames(g.prt) = c('sample', 'prt.exp')
  
  ## merge the two datasets
  g.exp = merge(g.exp, g.prt)
  g.exp$sample = gsub("ID_", "", g.exp$sample)
  
  ## add Status group
  g.exp$status.group = "CO"
  g.exp[g.exp$sample %in% bb$Subj_ID[bb$Status=="Neuro_AD"], 'status.group'] = "AD"
  g.exp$status.group = factor(g.exp$status.group, levels = c('CO','AD'))
  
  ## add clusters group 
  g.exp$group = 'CO'
  g.exp[g.exp$sample %in% bb$Subj_ID[bb$best.cluster==1 & bb$Status=="Neuro_AD"], 'group'] <- 'AD_C1'
  g.exp[g.exp$sample %in% bb$Subj_ID[bb$best.cluster==2 & bb$Status=="Neuro_AD"], 'group'] <- 'AD_C2'
  g.exp[g.exp$sample %in% bb$Subj_ID[bb$best.cluster==3 & bb$Status=="Neuro_AD"], 'group'] <- 'AD_C3'
  g.exp[g.exp$sample %in% bb$Subj_ID[bb$best.cluster==4 & bb$Status=="Neuro_AD"], 'group'] <- 'AD_C4'
  
  ## add ADAD, Prs_symptomatic, and TREM2 carrier 
  g.exp.adad.r <- melt(rnaseq[rownames(rnaseq) == gene, colnames(rnaseq) %in% adad.ids])
  g.exp.adad.p <- melt(prt.data[rownames(prt.data) == gene, colnames(prt.data) %in% adad.ids])
  colnames(g.exp.adad.r) <- c('sample', 'rnaseq.exp')
  colnames(g.exp.adad.p) <- c('sample', 'prt.exp')
  g.exp.adad <- merge(g.exp.adad.r, g.exp.adad.p)
  g.exp.adad$status.group <- 'ADAD'
  g.exp.adad$group <- 'ADAD'
  
  g.exp.prs.r <- melt(rnaseq[rownames(rnaseq) == gene, colnames(rnaseq) %in% prs.ids])
  g.exp.prs.p <- melt(prt.data[rownames(prt.data) == gene, colnames(prt.data) %in% prs.ids])
  colnames(g.exp.prs.r) <- c('sample', 'rnaseq.exp')
  colnames(g.exp.prs.p) <- c('sample', 'prt.exp')
  g.exp.prs <- merge(g.exp.prs.r, g.exp.prs.p)
  g.exp.prs$status.group <- 'PreSymp'
  g.exp.prs$group <- 'PreSymp'
  
  #g.exp.trem2.r <- melt(rnaseq[rownames(rnaseq) == gene, colnames(rnaseq) %in% trem2.ids])
  #g.exp.trem2.p <- melt(prt.data[rownames(prt.data) == gene, colnames(prt.data) %in% trem2.ids])
  #colnames(g.exp.trem2.r) <- c('sample', 'rnaseq.exp')
  #colnames(g.exp.trem2.p) <- c('sample', 'prt.exp')
  #g.exp.trem2 <- merge(g.exp.trem2.r, g.exp.trem2.p)
  #g.exp.trem2$status.group <- 'TREM2'
  #g.exp.trem2$group <- 'TREM2'
  
  ## merge with previous 
  g.exp$batch <- NULL
  g.exp <- rbind(g.exp, g.exp.adad)
  g.exp <- rbind(g.exp, g.exp.prs)
  #g.exp <- rbind(g.exp, g.exp.trem2)
  
  ## plot status
  p1 = ggplot(g.exp, aes(x=status.group, y=log2(rnaseq.exp+1))) + geom_boxplot(aes(fill=status.group), show.legend = F) +
    labs(x="", y="Log2(FPKM+1)") + ggtitle(paste0('RNAseq expression for ', gene)) + theme_bw() + 
    geom_jitter(position=position_jitter(0.3), size=1) +
    theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
          axis.text.x = element_text( vjust= 1, size=12, face="bold"), 
          axis.text.y = element_text(size=12, face="bold"),
          axis.title.y = element_text(size=12, face="bold"),
          legend.position = "none") +
    #scale_color_manual(name="", values= c('M'='gray40', 'S'='black', 'OV'='gray80'), labels=c('M'='Mend', 'S'='SUNSHINE', 'OV'='Overlap')) +
    geom_signif(comparisons= list(c('AD', 'CO'), c('ADAD', 'CO'), c('PreSymp', 'CO')), step_increase=0.1, textsize = 4, fontface="bold", map_signif_level=function(p)sprintf("p = %.2g", p), show.legend=F)
  
  ## plot clusters 
  g.exp$group = factor(g.exp$group, levels = c('CO', 'ADAD', 'PreSymp','AD_C1', 'AD_C2','AD_C3', 'AD_C4'))
  
  p2 = ggplot(g.exp, aes(x=group, y=log2(rnaseq.exp+1))) + geom_boxplot(aes(fill=group), show.legend=F) +  
    geom_jitter(position=position_jitter(0.3), size=1) +
    labs(x="", y="") + ggtitle(paste0('RNAseq expression for ', gene,' per cluster')) + theme_bw() +
    theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
          axis.text.x = element_text(vjust= 0.5, size=12, face="bold"), 
          axis.text.y = element_text(size=12, face="bold"),
          legend.position = "none") + 
    #scale_color_manual(name="", values= c('M'='gray40', 'S'='black', 'OV'='gray80'), labels=c('M'='Mend', 'S'='SUNSHINE', 'OV'='Overlap')) +
    geom_signif(comparisons= list(c('CO', 'AD_C1'), c('CO', 'AD_C2'), c('CO', 'AD_C3'), c('CO', 'AD_C4'), c('AD_C4', 'ADAD'), c('AD_C4', 'PreSymp')), 
                step_increase=0.1, textsize = 4, fontface="bold", map_signif_level=function(p)sprintf("p = %.2g", p), show.legend=F)
  
  ############### extract protein expression 
  ## plot status
  p3 = ggplot(g.exp[g.exp$prt.exp >0, ], aes(x=status.group, y=log2(prt.exp+1))) + geom_boxplot(aes(fill=status.group), show.legend = F) + 
    geom_jitter(position=position_jitter(0.3), size=1) +
    labs(x="", y="Log2(prt.exp+1)") + ggtitle(paste0('Protein expression for ', gene)) + theme_bw() +
    theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
          axis.text.x = element_text( vjust= 1, size=12, face="bold"), 
          axis.text.y = element_text(size=12, face="bold"),
          axis.title.y = element_text(size=12, face="bold"),
          legend.position = "none") +
    #scale_color_manual(name="", values= c('M'='gray40', 'S'='black', 'OV'='gray80'), labels=c('M'='Mend', 'S'='SUNSHINE', 'OV'='Overlap')) +
    geom_signif(comparisons= list(c('AD', 'CO'), c('ADAD', 'CO'), c('PreSymp', 'CO')), step_increase=0.1, textsize = 4, fontface="bold", map_signif_level=function(p)sprintf("p = %.2g", p), show.legend=F)
  
  ### plot per cluster 
  p4 = ggplot(g.exp[g.exp$prt.exp >0, ], aes(x=group, y=log2(prt.exp+1))) + geom_boxplot(aes(fill=group), show.legend = F) +  
    geom_jitter(position=position_jitter(0.3), size=1) +
    labs(x="", y="") + ggtitle(paste0('Protein expression for ', gene,' per cluster')) + theme_bw() +
    theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
          axis.text.x = element_text(vjust= 0.5, size=12, face="bold"), 
          axis.text.y = element_text(size=12, face="bold")) +
          #legend.position = "bottom", legend.text = element_text(size = 10)) +
    #scale_color_manual(name="", values= c('M'='gray40', 'S'='black', 'OV'='gray80'), 
    #                   labels=c('M'=paste0('Mend (',table(g.exp$batch)[['M']],')'), 'S'= paste0('SUNSHINE (',table(g.exp$batch)[['S']],')'),
    #                            'OV'=paste0('Overlap (', table(g.exp$batch)[['OV']],')')), guide= guide_legend(override.aes = list(size=5, alpha = 1))) +
    geom_signif(comparisons= list(c('CO', 'AD_C1'), c('CO', 'AD_C2'), c('CO', 'AD_C3'), c('CO', 'AD_C4'), c('AD_C4', 'ADAD'), c('AD_C4', 'PreSymp')),
                step_increase=0.1, textsize = 4, fontface="bold", map_signif_level=function(p)sprintf("p = %.2g", p), show.legend=F)
  
  ## plot the correlation between protein exp and RNAseq 
  cor.rr <- by(g.exp[,2:3], g.exp$group, function(x) {cor(x$rnaseq.exp, x$prt.exp, method="spearman")})
  cor.rr <- signif(cor.rr,digits = 2)
                      
  p5 <- ggplot(g.exp, aes(x=log2(prt.exp+1), y=log2(rnaseq.exp+1))) + geom_point(aes(color=group)) + 
    facet_wrap( ~ group, nrow = 1, scales = "free", labeller = labeller(group = 
                                                                                     c("Control" = paste0("CO (",cor.rr["Control"],")"),
                                                                                       "ADAD" = paste0("ADAD (",cor.rr["ADAD"],")"),
                                                                                       "PRS" = paste0("PRS (",cor.rr["PRS"],")"),
                                                                                       "TREM2" = paste0("TREM2 (",cor.rr["TREM2"],")"),
                                                                                       "AD_C1" = paste0("AD_C1 (",cor.rr["AD_C1"],")"),
                                                                                       "AD_C2" = paste0("AD_C2 (",cor.rr["AD_C2"],")"),
                                                                                       "AD_C3" = paste0("AD_C3 (",cor.rr["AD_C3"],")"),
                                                                                       "AD_C4" = paste0("AD_C4 (",cor.rr["AD_C4"],")")
                                                                                       )
                                                                                      )) +
        labs(x="Log2(protein.exp+1)", y="Log2(rnaseq.exp+1)") + ggtitle(paste0('Correlartion between RNASeq and protein expression for ', gene)) + theme_bw() +
        theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
              axis.text.x = element_text(vjust= 1, size=8, face="bold"), 
              axis.text.y = element_text(size=8, face="bold"),
              legend.position = "none", 
              strip.text = element_text(size=10, face="bold"))
  
  #grid.arrange(arrangeGrob(grobs=list(p1, p2, p3, p4, p5), ncol=2, layout_matrix=mat, heights = c(4,4,2), widths=c(2,4)))
  pdf('omics_integration/iCluster_output/Figures/SNCA_exp_profiles.pdf', width=12, height = 10)
  grid.arrange(arrangeGrob(grobs=list(p1, p2, p3, p4), ncol=2, widths=c(2,4)))
  dev.off()
  
}

#cor(g.exp$rnaseq.exp[g.exp$group=="Cluster1"], g.exp$prt.exp[g.exp$group=="Cluster1"], method="spearman")
#pdf('omics_integration/Figures/sig_genes_plots_v2.pdf', width=20, height = 100)
#grid.arrange(grobs=all.genes, ncol=10)
dev.off()
##############################################################


##############################################################
##################### Survival Association ###################
##############################################################
library(surviplot)
library(survival)
#library(survminer)
clin = bb[, c('Subj_ID','AOD')]
colnames(clin) = c('sample', 'time')
## add clusters group 
clin$group = 'Control'
clin[clin$sample %in% bb$Subj_ID[bb$best.cluster %in% c(1,2,3) & bb$Status=="Neuro_AD"], 'group'] <- 'Knight-C1-3'
clin[clin$sample %in% bb$Subj_ID[bb$best.cluster==4 & bb$Status=="Neuro_AD"], 'group'] <- 'Knight-C4'
max.time = max(clin$time, na.rm = T)

# #groups = c('Cluster1', 'Cluster2', 'Cluster3', 'Cluster4')
# groups = c('Cluster123', 'Cluster4')
# cox.res = NULL
# for (ii in 1:(length(groups)-1)){
#   for (jj in (ii+1):length(groups)){
#     cat(groups[ii], '_', groups[jj], '\n')
#     
#     thisy = clin[clin$group %in% c(groups[ii], groups[jj]),]
#     cox = summary(coxph(Surv(time) ~ group, data=thisy))
#     p = cox$sctest['pvalue']
#     hr = cox$conf.int[1, c(1, 3, 4)]
#     d = data.frame(comp=paste0(groups[ii], '_',groups[jj]), pval=p, HR_exp_coef=hr[1], HR_lower_0.95=hr[2], HR_upper_0.95=hr[3], stringsAsFactors = F, row.names = NULL )
#     cox.res = rbind(cox.res, d)
#   }
# }
# 
# ## extract cluster 4 
# cox.res_cl4 = cox.res[grepl("Cluster4", cox.res$comp), ]

cox = summary(coxph(Surv(time) ~ group, data=clin[clin$group !="Control", ]))
p <- cox$sctest['pvalue']
hr <- cox$conf.int[1, c(1, 3, 4)]

diff.months <- abs(mean(bb$duration[bb$best.cluster==4]) - mean(bb$duration[bb$best.cluster %in% c(1,2,3)], na.rm= T))
diff.months <- round(diff.months * 12)

## plot KM curve 
pdf('omics_integration/iCluster_output/Figures/Final/survival_AOD_two_groups.pdf', width=5, height = 5)
surviplot(Surv(time) ~ group, data=clin[clin$group !="Control", ], ylab='Survival Proportion',xlim=c(55,max.time), main ="", 
          xlab='Age at onset (years)', cex.main=1, mark.time=TRUE, col=c("skyblue3", '#E7298A'), lwd = 3, cex.axis=2, cex.lab= 1.5)
dev.off()


########################################################
######################### PCA ######################### 
########################################################
rnaseq.pca = as.data.frame(t(rnaseq[rownames(rnaseq) %in% sigfeatures$expression, ]))
rnaseq.pca$Subj_ID = gsub("ID_", "", rownames(rnaseq.pca))
rnaseq.pca =  merge(rnaseq.pca, bb[, c('Subj_ID', 'best.cluster')])
pca = prcomp(rnaseq.pca[, !colnames(rnaseq.pca) %in% c('Subj_ID', 'best.cluster')])
#pp = as.data.frame(pca$x)
plot(pca$x[,1], pca$x[,2])
rnaseq.pca$best.cluster = as.factor(rnaseq.pca$best.cluster)
autoplot(pca, data = rnaseq.pca, colour = 'best.cluster')


###
#metab.m = melt(metab, id.vars = c('Subj_ID', 'best.cluster'))
#p = ggplot(metab.m, aes(x=factor(best.cluster), y=value, group=factor(best.cluster))) + geom_boxplot(outlier.shape=NA, outlier.size=0.5) + geom_jitter()


### plot expression heatmap based on the clusters 
rnaseq.clusters = rnaseq
colnames(rnaseq.clusters) = gsub("ID_", "", colnames(rnaseq.clusters))
rnaseq.clusters = rnaseq.clusters[sigfeatures$expression, 
                                  c(as.character(bb$Subj_ID[bb$Status=="Neuro_CO"]),
                                    as.character(bb$Subj_ID[bb$best.cluster==1 & bb$Status=="Neuro_AD"]), 
                                    as.character(bb$Subj_ID[bb$best.cluster==2 & bb$Status=="Neuro_AD"]), 
                                    as.character(bb$Subj_ID[bb$best.cluster==3 & bb$Status=="Neuro_AD"]), 
                                    as.character(bb$Subj_ID[bb$best.cluster==4 & bb$Status=="Neuro_AD"])
                                  )]
exp.data2 <-log2(rnaseq.clusters+1)

colBreaks = seq(-5, 5, 1)
#col.pan <- colorpanel(length(colBreaks)-1, "blue", "white", "red")
col.pan <- colorRampPalette(c("red","white","darkgreen"))(length(colBreaks)-1)

pdf('omics_integration/Figures/xxxxx/exp_heatmap_with_clusters2.pdf', width = 40, height = 100)
heatmap.2(as.matrix(exp.data2), 
          dendrogram = "row", 
          Colv = FALSE, 
          Rowv = TRUE, cexRow=1,
          scale = "row", 
          col = col.pan, breaks= colBreaks,
          key = TRUE, keysize=0.75, 
          key.par=list(cex=2, mar=c(5,1,5,5)),
          #key.par=list(mar=c(0,5,25,10)), 
          density.info = "none", 
          key.title = NA, 
          key.xlab = "Expression",
          trace = "none",
          colsep = c(23, 137, 183, 236), 
          #ColSideColors = c(rep('gray80', 23), rep('#dd1c77', 114), rep('#2ca25f', 46), rep('#43a2ca', 53) ,rep('#d95f0e', 45) ),
          sepwidth=0.3, 
          sepcolor = 'blue',
          margin=c(25,20), 
          lhei=c(2,25), lwid=c(2,10))
#margins = c(7, 15))
dev.off()

######################################################################
################ prepare data for gsea
######################################################################
bb <- bb <- read.table("omics_integration/iCluster_output/best.cluster.memership_v1_full_AD_only.tsv", header = T, sep="\t", stringsAsFactors = F)
bb$group <- 'OT'
bb[!bb$best.cluster %in% c(1,2,3), 'group'] <- 'C4'
bb$Subj_ID <- paste0('ID_', bb$Subj_ID)

res <- read.table('omics_integration/iCluster_output/DE/C4vsC123/de_res_sig_C4vsC123.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')
res <- rnaseq[rnaseq$GeneName %in% res$GeneName, c('GeneName', bb$Subj_ID[bb$group=='C4'], bb$Subj_ID[bb$group=="OT"])]
res$DESCRIPTION = NA
res <- res[, c('GeneName', 'DESCRIPTION', bb$Subj_ID[bb$group=='C4'], bb$Subj_ID[bb$group=="OT"])]
colnames(res)[colnames(res)=="GeneName"] <- 'NAME'
write.table(res, file="omics_integration/data/exp_data_C4vsC123.gct", sep="\t", quote = F, row.names = F)

## make the phenotype file 
ph <- c(rep('C4', length(bb$Subj_ID[bb$group=='C4']) ), rep('OT', length(bb$Subj_ID[bb$group=="OT"])) )
ph <- paste(ph, collapse = " ")
write.table(ph, file="omics_integration/data/phenotype_class_C4vsC123.cls", sep=" ", quote = F, row.names = F)
#################################################################################

#######################################################################
## check results in metabolomics
#######################################################################
metab.m <- melt(metab, id.vars = "COMP.ID")
colnames(metab.m) <- c('COMP.ID', 'Subj_ID', 'reading')
metab.m$Subj_ID <- gsub("ID_", "", metab.m$Subj_ID)
metab.m <- merge( bb[, c('Subj_ID', 'best.cluster')], metab.m)

p <- ggplot(metab.m, aes(x=factor(best.cluster), y=reading, group=factor(best.cluster))) + geom_boxplot(outlier.shape=NA, outlier.size=0.5) 
p <- p + geom_jitter(position=position_jitter(0.3),  size=1, color=1) + theme_bw() + guides(color = FALSE) 
p <- p + labs(x='', y='Reading') + ggtitle('')
p

######################################################################
################ Perform DE analysis ################ 
######################################################################
dir.create("omics_integration/iCluster_output/DE")
## read cluster solution data 
bb = read.table("omics_integration/iCluster_output/best.cluster.memership_v1_full.tsv", header = T, sep="\t", stringsAsFactors = F)
#bb = read.table("omics_integration/iCluster_output/best.cluster.memership_v1_full_AD_only.tsv", header = T, sep="\t", stringsAsFactors = F) # #for supplementary table
## remove flagged samples 
bb <- bb[!bb$Subj_ID %in% as.integer(gsub("ID_","", flagged.samples)), ]

## add cell lines proportion information 
algs = c( "ssFrobenius", "meanProfile")
sunshine.path <- "/home/general/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/05.-Analyses/deconvolution_analysis"
mend.path <- "/home/general/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/05.-Analyses/deconvolution_analysis"

cell.prop.res <- NULL
for (alg in algs) {
  s.res <- read.table(paste0(sunshine.path, '/deconvolution_', alg,'_results.tsv'), header = T, check.names = F, stringsAsFactors = F)
  m.res <- read.table(paste0(mend.path, '/deconvolution_', alg,'_results.tsv'), header = T, check.names = F, stringsAsFactors = F)
  cell.prop.res <- rbind(cell.prop.res, s.res)
  cell.prop.res <- rbind(cell.prop.res, m.res)
}
colnames(cell.prop.res)[colnames(cell.prop.res)=="Subject"] <- "Sample_Name"
cell.prop.res$Sample_Name <- gsub("H_VY.", "H_VY-", cell.prop.res$Sample_Name) 

## merge count data
cell.prop.res$Subj_ID <- NA
rnase_counts <- merge(mend.cts, sunshine.cts)
for (i in 1:length(colnames(rnase_counts[,-1:-3]))) {
  col.name <- colnames(rnase_counts[,-1:-3])[i]
  subj.id <- sunshine.tech[sunshine.tech$Sample_Name==col.name, 'Subj_ID']
  if (length(subj.id) !=0) {
    colnames(rnase_counts)[colnames(rnase_counts)==col.name] <- subj.id
    cell.prop.res[cell.prop.res$Sample_Name==col.name, 'Subj_ID'] <- subj.id
  } else {
    subj.id <- mend.tech[mend.tech$Sample_ID==col.name, 'Subj_ID']
    colnames(rnase_counts)[colnames(rnase_counts)==col.name] <- subj.id
    cell.prop.res[cell.prop.res$Sample_Name==col.name, 'Subj_ID'] <- subj.id
  }
}

colnames(rnase_counts) <- gsub("MAP_", "", colnames(rnase_counts))
colnames(rnase_counts) <- gsub("MAC_", "", colnames(rnase_counts))
colnames(rnase_counts) <- gsub("DIAN_", "", colnames(rnase_counts))

cell.prop.res$Subj_ID <- gsub("MAP_", "", cell.prop.res$Subj_ID)
cell.prop.res$Subj_ID <- gsub("MAC_", "", cell.prop.res$Subj_ID)
cell.prop.res$Subj_ID <- gsub("DIAN_", "", cell.prop.res$Subj_ID)

## remove QC failed samples 
bb <- bb[!bb$Subj_ID %in% flagged.samples, ]
rnase_counts <- rnase_counts[, !colnames(rnase_counts) %in% flagged.samples]

## merge with cell proportions 
# get the max percentage for duplicated IDs
cell.prop.res <- cell.prop.res[cell.prop.res$Algorithm=='ssFrobenius', c('Subj_ID', 'Astrocyte', 'Neuron')]
xx = aggregate(cell.prop.res$Astrocyte, by = list(cell.prop.res$Subj_ID), max)
colnames(xx) <- c('Subj_ID', 'Astrocyte')
cell.prop.res <- merge(cell.prop.res, xx, by=c('Subj_ID', 'Astrocyte'))
bb <- merge(bb, cell.prop.res)

## write results for supplementary table 
res_for_supp <- bb[, c('Subj_ID', 'best.cluster', 'Algorithm', 'Astrocyte', 'Neuron', 'Microglia', 'Oligodendrocyte')]
colnames(res_for_supp) <- c('Subj_ID', 'Cluster', 'Algorithm', 'Astrocyte', 'Neuron', 'Microglia', 'Oligodendrocyte')
write.table(res_for_supp, file='omics_integration/Cell_proportions_KnightADRC.tsv', sep="\t", quote = F, row.names = F)

## make groups 
co.samples <-  as.character(bb[bb$Status=="Neuro_CO", "Subj_ID"])
c1.samples <-  as.character(bb[bb$Status=="Neuro_AD" & bb$best.cluster==1, "Subj_ID"]) 
c2.samples <-  as.character(bb[bb$Status=="Neuro_AD" & bb$best.cluster==2, "Subj_ID"]) 
c3.samples <-  as.character(bb[bb$Status=="Neuro_AD" & bb$best.cluster==3, "Subj_ID"]) 
c4.samples <-  as.character(bb[bb$Status=="Neuro_AD" & bb$best.cluster==4, "Subj_ID"]) 

#adad.samples <- gsub("ID_", "", adad.ids)
#preysym.samples <- gsub("ID_", "", prs.ids)

## for AD vs CO
# bb <- clin.data
# bb <- bb[!is.na(bb$Status) & bb$Status %in% c('Neuro_CO', 'Neuro_AD'), ]
# co.samples <-  as.character(bb[!is.na(bb$Status) & bb$Status=="Neuro_CO", "Subj_ID"])
# ad.samples <-  as.character(bb[!is.na(bb$Status) & bb$Status=="Neuro_AD", "Subj_ID"])
# rownames(bb) <- bb$Subj_ID
# bb <- bb[c(co.samples, ad.samples), ]

groups = factor(c( rep('CO', length(co.samples)), rep('C1', length(c1.samples)), rep('C2', length(c2.samples)),
                   rep('C3', length(c3.samples)), rep('C4', length(c4.samples)), 
                   rep('ADAD', length(adad.samples)), rep('PreSymp', length(preysym.samples))), 
                   levels = c("CO","ADAD", "PreSymp", "C1", "C2", "C3", "C4")) 
bb$group <- 'CO'
bb[bb$Subj_ID %in% c1.samples, 'group'] <- 'C1'
bb[bb$Subj_ID %in% c2.samples, 'group'] <- 'C2'
bb[bb$Subj_ID %in% c3.samples, 'group'] <- 'C3'
bb[bb$Subj_ID %in% c4.samples, 'group'] <- 'C4'

bb$group <- factor(bb$group, levels =c('CO', 'C1', 'C2', 'C3', 'C4'))

## define groups for comparing one cluster to the others 
bb$groupC1 <- 'CO'
bb[bb$Subj_ID %in% c1.samples, 'groupC1'] <- 'C1'
bb[bb$Subj_ID %in% c(c2.samples, c3.samples, c4.samples), 'groupC1'] <- 'C234'
bb$groupC1 <- factor(bb$groupC1, levels =c('CO', 'C1', 'C234'))

bb$groupC2 <- 'CO'
bb[bb$Subj_ID %in% c2.samples, 'groupC2'] <- 'C2'
bb[bb$Subj_ID %in% c(c1.samples, c3.samples, c4.samples), 'groupC2'] <- 'C134'
bb$groupC2 <- factor(bb$groupC2, levels =c('CO', 'C2', 'C134'))

bb$groupC3 <- 'CO'
bb[bb$Subj_ID %in% c3.samples, 'groupC3'] <- 'C3'
bb[bb$Subj_ID %in% c(c1.samples, c2.samples, c4.samples), 'groupC3'] <- 'C124'
bb$groupC3 <- factor(bb$groupC3, levels =c('CO', 'C3', 'C124'))

bb$groupC4 <- 'CO'
bb[bb$Subj_ID %in% c4.samples, 'groupC4'] <- 'C4'
bb[bb$Subj_ID %in% c(c1.samples, c2.samples, c3.samples), 'groupC4'] <- 'C123'
bb$groupC4 <- factor(bb$groupC4, levels =c('CO', 'C4', 'C123'))

## order data based on groups
rownames(rnase_counts) <- paste0(rnase_counts$GeneID,'_', rnase_counts$GeneName)
rnase_counts$GeneID <- NULL
rnase_counts$GeneName <- NULL
rnase_counts$GeneBiotype <- NULL
rnase_counts <- rnase_counts[, c(co.samples, c1.samples, c2.samples,c3.samples,c4.samples)]
#rnase_counts <- rnase_counts[, c(co.samples, ad.samples)] ## for AD vs CO
rownames(bb) = bb$Subj_ID
bb <- bb[c(co.samples, c1.samples, c2.samples,c3.samples,c4.samples), ]
## check samples consistency 
all(rownames(bb) == colnames(rnase_counts) )


## function to filter lowly expressed genes 
rnase_counts_norm <- cpm(rnase_counts)
selectGenes <- function(grp1, grp2) {
  grp1.t <- round(length(grp1) * 0.25)
  grp2.t <- round(length(grp2) * 0.25)
  keep <- (rowSums(rnase_counts_norm[,grp1] > 0.5) >= grp1.t) | (rowSums (rnase_counts_norm[,grp2]> 0.5) >= grp2.t)
  if (all(rownames(rnase_counts) != rownames(rnase_counts_norm))) {
    stop('rownames of the normalized counts are not the same as the raw counts ')
  } else {
    cmp.counts <- rnase_counts[keep, ]
  }
  return (cmp.counts)
}

## Function for fitting the DESeq2 model
DESeqModelFit <- function(count_data, model.group) {
   dds <- DESeqDataSetFromMatrix(countData = count_data, 
                                 colData = bb, 
                                 design = formula(paste("~ Sex + AOD + Astrocyte + Neuron +", model.group)) )
  
  # dds <- DESeqDataSetFromMatrix(countData = count_data, 
  #                               colData = bb, 
  #                               design = formula(~Sex + AOD + Status))
  
  ## Test for DE 
  dds <- DESeq(dds)
  resultsNames(dds)
  
  return (dds)
}

######################################################################
# check GWAS hits 
reported_genes <- as.data.frame(read_excel("omics_integration/data/Reported_Disease_Genes.xlsx", col_names = T, sheet =1, trim_ws = TRUE))
prioritized_genes <- as.data.frame(read_excel("omics_integration/data/Reported_Disease_Genes.xlsx", col_names = T, sheet =2, trim_ws = TRUE))

add_gwas_hits <- function (de.res) {

    de.res$GWAS.Num.Hits <- 0
    de.res$GWAS.Locus <- NA
    de.res$GWAS.FunctionalEvidence <- NA
    de.res$In.Pioritized.Genes.List <- 'N'
    for (i in 1:nrow(de.res)) {
      g <- de.res$GeneName[i]
      g.gwas <-  nrow(reported_genes[reported_genes$Gene==g & !is.na(reported_genes$Locus), ])
      in.pgl <-  nrow(prioritized_genes[prioritized_genes$Gene==g, ])
      if (in.pgl !=0) { de.res$In.Pioritized.Genes.List[i] = 'Y'}
      g.gwas.locus <- paste(reported_genes[reported_genes$Gene==g & !is.na(reported_genes$Locus), 'Locus'], collapse = "|")
      g.gwas.FE <- unique(reported_genes[reported_genes$Gene==g & !is.na(reported_genes$Locus), 'Functional Evidence'])
      g.gwas.FE <- as.character(na.exclude(g.gwas.FE))
      
      de.res$GWAS.Num.Hits[i] <- g.gwas
      if (length(g.gwas.locus) !=0) {
        de.res$GWAS.Locus[i] <- g.gwas.locus      
      } else {
        de.res$GWAS.Locus[i] <- NA
      }
      
      if (length(g.gwas.FE) !=0) {
        de.res$GWAS.FunctionalEvidence[i] <- g.gwas.FE     
      } else {
        de.res$GWAS.FunctionalEvidence[i] <- NA
      }
      
    }
    ## write results with 
    return (de.res) 
}

# ----------------------------------------------------------------------------------------------------------------------

## Function to perform pathways analysis using gage 
## prepare data 
#kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
kegg.sets.hs <- kegg.gsets(species="human")
kegg.sets.hs <- kegg.sets.hs$kg.sets
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

run_kegg_go <- function (res, cc, plot.top=0) {
  
  dir.create(paste0('omics_integration/iCluster_output/DE/', cc,'/KEGG_Pathways'))
  dir.create(paste0('omics_integration/iCluster_output/DE/', cc,'/GO_Pathways'))
  
  # patheay analysis using gage
  res$entrez <- mapIds(org.Hs.eg.db,
                           keys=res$GeneName,
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
  
  foldchanges <- res$log2FoldChange
  names(foldchanges) <- res$entrez
  
  # Get KEGG results
  kegg_res <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
  
  ## extract pathways per direction 
  kegg_res_up <- data.frame(id=rownames(kegg_res$greater), kegg_res$greater)
  kegg_res_dn <- data.frame(id=rownames(kegg_res$less), kegg_res$less)
  
  write.table(kegg_res_up, file=paste0('omics_integration/iCluster_output/DE/',cc, '/KEGG_Pathways/', cc, '_KEGG_UP.tsv'), sep="\t", quote = F, row.names = F)
  write.table(kegg_res_dn, file=paste0('omics_integration/iCluster_output/DE/',cc, '/KEGG_Pathways/', cc, '_KEGG_DN.tsv'), sep="\t", quote = F, row.names = F)
  
  if (plot.top !=0) {
    ## plot the top 5 IDs
    top_up_kegg_ids <- substr(rownames(kegg_res_up)[1:plot.top], start=1, stop=8)
    top_dn_kegg_ids <- substr(rownames(kegg_res_dn)[1:plot.top], start=1, stop=8)
    
    # Define plotting function for applying later
    plot_pathway <- function(pid) {
      pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
    }
    # plot multiple pathways (plots saved to disk and returns a throwaway list object)
    setwd(paste0('omics_integration/iCluster_output/DE/', cc, '/KEGG_Pathways'))
    tmp <- sapply(c(top_up_kegg_ids, top_dn_kegg_ids), function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))
    setwd('/home/eteleeb/projects')
  }
 
  ## Get GO results 
  #go.sets.hs <- go.sets.hs[go.subs.hs$BP]. ## for BP only 
  go_bp_res <- gage(foldchanges, gsets=go.bp.gs, same.dir=TRUE)
  go_cc_res <- gage(foldchanges, gsets=go.cc.gs, same.dir=TRUE)
  go_mf_res <- gage(foldchanges, gsets=go.mf.gs, same.dir=TRUE)
  
  #lapply(go_res, head)
  ## extract pathways per direction 
  go_bp_up <- data.frame(id=rownames(go_bp_res$greater), go_bp_res$greater)
  go_cc_up <- data.frame(id=rownames(go_cc_res$greater), go_cc_res$greater)
  go_mf_up <- data.frame(id=rownames(go_mf_res$greater), go_mf_res$greater)
  
  go_bp_dn <- data.frame(id=rownames(go_bp_res$less), go_bp_res$less)
  go_cc_dn <- data.frame(id=rownames(go_cc_res$less), go_cc_res$less)
  go_mf_dn <- data.frame(id=rownames(go_mf_res$less), go_mf_res$less)
  
  write.table(go_bp_up, file=paste0('omics_integration/iCluster_output/DE/',cc, '/GO_Pathways/', cc,'_GO_BP_UP.tsv'), sep="\t", quote = F, row.names = F)
  write.table(go_cc_up, file=paste0('omics_integration/iCluster_output/DE/',cc, '/GO_Pathways/', cc,'_GO_CC_UP.tsv'), sep="\t", quote = F, row.names = F)
  write.table(go_mf_up, file=paste0('omics_integration/iCluster_output/DE/',cc, '/GO_Pathways/', cc,'_GO_MF_UP.tsv'), sep="\t", quote = F, row.names = F)
  
  write.table(go_bp_dn, file=paste0('omics_integration/iCluster_output/DE/',cc, '/GO_Pathways/', cc,'_GO_BP_DN.tsv'), sep="\t", quote = F, row.names = F)
  write.table(go_cc_dn, file=paste0('omics_integration/iCluster_output/DE/',cc, '/GO_Pathways/', cc,'_GO_CC_DN.tsv'), sep="\t", quote = F, row.names = F)
  write.table(go_mf_dn, file=paste0('omics_integration/iCluster_output/DE/',cc, '/GO_Pathways/', cc,'_GO_MF_DN.tsv'), sep="\t", quote = F, row.names = F)

}

###############################################
## uenrichR
dbs <- listEnrichrDbs()
run_enrichr <- function (de, cc) {
  #gene.list <- unlist(strsplit(rownames(res_sig), "_"))
  #gene.list <- gene.list[!grepl('ENSG', gene.list)]
  en.up <- enrichr(de$GeneName[de$direction =="up"], databases = dbs$libraryName)  
  en.dn <- enrichr(de$GeneName[de$direction =="dn"], databases = dbs$libraryName)  
  save(en.up, file=paste0('omics_integration/iCluster_output/DE/', cc,'/', cc,'_enrichR_res_up.RData'))
  save(en.dn, file=paste0('omics_integration/iCluster_output/DE/', cc,'/', cc,'_enrichR_res_dn.RData'))
}
###############################################

# ## Checking the normalization
# pdf('omics_integration/iCluster_output/DE/raw_norm_density_counts.pdf', width = 15, height = 15)
# par(mfrow=c(2,2),cex.lab=0.7)
# boxplot(log2(counts(dds)+1),  col=as.numeric(groups), cex.axis=0.7, las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
# boxplot(log2(counts(dds, normalized=TRUE)+1),  col=as.numeric(groups), cex.axis=0.7, las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
# 
# plotDensity(log2(counts(dds)+1),  col=as.numeric(groups), xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
# plotDensity(log2(counts(dds, normalized=TRUE)+1), col=as.numeric(groups), xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid()) 
# dev.off()

## loop through the comparisons 
comps <- c('C4vsCO', 'C4vsC1', 'C4vsC2', 'C4vsC3', 'C1vsCO','C2vsCO','C3vsCO','C1vsC2','C1vsC3', 'C2vsC3', 'C1vsC234', 'C2vsC134', 'C3vsC124', 'C4vsC123')
for (cmp in comps) {
  cat('Running for:', cmp, '\n')
  dir.create(paste0('omics_integration/iCluster_output/DE/', cmp))
  
  ## extract groups 
  mygroup <- "group"
  if (cmp %in% c('C1vsC234', 'C2vsC134', 'C3vsC124', 'C4vsC123') ) {
    grp1.samples <- eval(parse(text = paste0(tolower(unlist(strsplit(cmp, "vs"))[1]), '.samples')))
    xx <- unlist(strsplit(tolower(unlist(strsplit(cmp, "vs"))[2]),""))
    xx <- paste0('c', xx[xx !="c"],'.samples')
    grp2.samples <- c(eval(parse(text =xx[1])), eval(parse(text =xx[2])), eval(parse(text =xx[3])))
    mygroup <- paste0("group", unlist(strsplit(cmp, "vs"))[1])
  } else {
    grp1.samples <- eval(parse(text = paste0(tolower(unlist(strsplit(cmp, "vs"))[1]), '.samples')))
    grp2.samples <-  eval(parse(text = paste0(tolower(unlist(strsplit(cmp, "vs"))[2]), '.samples')))   
  }
  
  ## remove lowly expressed genes 
  grp.counts <- selectGenes(grp1.samples, grp2.samples)
  ## fit DESeq model
  dds <- DESeqModelFit(grp.counts, mygroup)
  
  ## extract result table
  dds.res <- results(dds, alpha=0.05, contrast = c(mygroup, unlist(strsplit(cmp, "vs"))[1], unlist(strsplit(cmp, "vs"))[2]), tidy = F)
  #dds.res <- results(dds, alpha=0.05, contrast = c(mygroup, 'Neuro_AD', 'Neuro_CO'), tidy = F)
  summary(dds.res)
  
  ## add annotation 
  gg <- as.data.frame(str_split_fixed(rownames(dds.res), "_", 2))
  colnames(gg) <- c('GeneID', "GeneName")
  dds.res <- cbind(dds.res, gg)
  dds.res <- as.data.frame(merge(meta.cols , dds.res))
  
  ## add status 
  dds.res$direction='nc'
  dds.res[dds.res$log2FoldChange > 0 & !is.na(dds.res$padj) & dds.res$padj < 0.05, 'direction'] <- 'up'
  dds.res[dds.res$log2FoldChange < 0 & !is.na(dds.res$padj) & dds.res$padj < 0.05, 'direction'] <- 'dn'
  ## sort genes based on padj
  dds.res <- dds.res[order(dds.res$padj), ]
  
  ## extract the mean expression of each group 
  normalized_counts <- counts(dds, normalized=TRUE)
  normalized_counts = normalized_counts[paste0(dds.res$GeneID,'_', dds.res$GeneName), ]
  dds.res$grp1 <- rowMeans(normalized_counts[, grp1.samples], na.rm=TRUE)
  dds.res$grp2 <- rowMeans(normalized_counts[, grp2.samples], na.rm = TRUE)
  colnames(dds.res)[colnames(dds.res) =="grp1"] <- paste0(unlist(strsplit(cmp, "vs"))[1], '.ExpMean')
  colnames(dds.res)[colnames(dds.res) =="grp2"] <- paste0(unlist(strsplit(cmp, "vs"))[2], '.ExpMean')
  
  ## plot volcano plot
  top.10.genes <- dds.res$GeneName[1:10]
  pdf(paste0('omics_integration/iCluster_output/DE/',cmp, '/volcano_plot_',cmp,'2.pdf'), width = 6, height = 8)
  par(mfrow=c(1,1))
  with(dds.res, plot(log2FoldChange, -log10(pvalue), pch=20, main=paste0("Volcano plot - ", cmp), xlim=c(min(dds.res$log2FoldChange), max(dds.res$log2FoldChange))))
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(dds.res, direction=="up"), points(log2FoldChange, -log10(pvalue), pch=20, col="orange3"))
  with(subset(dds.res, direction=="dn") , points(log2FoldChange, -log10(pvalue), pch=20, col="skyblue3"))
  with(subset(dds.res,direction=="nc") , points(log2FoldChange, -log10(pvalue), pch=20, col="gray60"))
  with(dds.res, text(log2FoldChange, -log10(pvalue), labels=ifelse(GeneName %in% c(top.10.genes, 'CLU'), GeneName, ''), cex=0.7, offset =1, adj=c(0.5,0.01)))
  legend ("topright", c('Down', 'Up', 'NC'), col=c('skyblue3', 'orange3', 'gray60'), pch=c(20,20), pt.cex=2.5)
  dev.off()
  
  ## extract significant results  
  dds.res.sig <- dds.res[dds.res$direction !="nc", ]
  dds.res.sig <- dds.res.sig[order(dds.res.sig$padj), ]
  
  ## check with GWAS hits
  cat('Checking with GWAS hits for:', cmp, '\n')
  dds.res.sig <- add_gwas_hits (dds.res.sig)
  
  write.table(dds.res, file=paste0('omics_integration/iCluster_output/DE/',cmp,'/de_res_', cmp,'.tsv'), sep="\t", quote = F, row.names = F)
  write.table(dds.res.sig, file=paste0('omics_integration/iCluster_output/DE/',cmp, '/de_res_sig_',cmp,'.tsv'), sep="\t", quote = F, row.names = F)
  
  ## run kegg pathways 
  run_kegg_go(dds.res.sig, cmp, 0)
  run_enrichr(dds.res.sig, cmp)
             
}


#plotCounts(dds, "ENSG00000145335.16_SNCA", intgroup = "group", normalized = TRUE)
## ----------------------------------------------------------------------------------------------------------------------------------
## Generate venn diagrams for clusters vs CO and among clusters 
## ----------------------------------------------------------------------------------------------------------------------------------
# plotExpMean <- function (d, c4.uniq) {
#   C1vsCO.all <- read.table('omics_integration/iCluster_output/DE/C1vsCO/de_res_C1vsCO.tsv', header =T, stringsAsFactors = F, sep="\t")
#   C2vsCO.all <- read.table('omics_integration/iCluster_output/DE/C2vsCO/de_res_C2vsCO.tsv', header =T, stringsAsFactors = F, sep="\t")
#   C3vsCO.all <- read.table('omics_integration/iCluster_output/DE/C3vsCO/de_res_C3vsCO.tsv', header =T, stringsAsFactors = F, sep="\t")
#   C4vsCO.all <- read.table('omics_integration/iCluster_output/DE/C4vsCO/de_res_C4vsCO.tsv', header =T, stringsAsFactors = F, sep="\t")
#   
#   exp.mean <- NULL
#   #for (cc in c('C1', 'C2', 'C3', 'C4')) {
#   c1 <- C1vsCO.all[C1vsCO.all$GeneName %in% c4.uniq, c('GeneName', 'C1.ExpMean')] 
#   c2 <- C2vsCO.all[C2vsCO.all$GeneName %in% c4.uniq, c('GeneName', 'C2.ExpMean')] 
#   c3 <- C3vsCO.all[C3vsCO.all$GeneName %in% c4.uniq, c('GeneName', 'C3.ExpMean')] 
#   c4 <- C4vsCO.all[C4vsCO.all$GeneName %in% c4.uniq, c('GeneName', 'C4.ExpMean')] 
#   co <- C4vsCO.all[C4vsCO.all$GeneName %in% c4.uniq, c('GeneName', 'CO.ExpMean')] 
#   exp.mean <- merge(merge(merge(merge(co, c1) , c2),c3), c4)
#   
#   colnames(exp.mean) <- gsub('.ExpMean', '', colnames(exp.mean))
#   exp.mean <- reshape::melt(exp.mean)
#   colnames(exp.mean) <- c('GeneName', 'Cluster', 'Exp')
#   p <- ggplot(exp.mean, aes(x=Cluster, y=log2(Exp+1), group=GeneName)) + geom_line(color ='#440154FF', show.legend = T) + 
#           labs(x="Cluster", y="Log2(Expression)") + ggtitle('') + theme_bw() + 
#           theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
#                 axis.text.x = element_text( vjust= 1, size=12, face="bold", color="black"), 
#                 axis.text.y = element_text(size=12, face="bold"), 
#                 axis.title.x = element_text(size=14, face="bold"),
#                 axis.title.y = element_text(size=14, face="bold"),
#                 legend.text = element_text(size=12), 
#                 legend.position = "bottom") + scale_color_viridis_d(name="", direction = -1)
#   pdf(paste0("omics_integration/iCluster_output/Figures/Final/exp_mean_C4_genes_", d, ".pdf"), width = 9, height = 5)
#   print (p)
#   dev.off()
#   
# }

# transcriptome 
C4vsCO <- read.table('omics_integration/iCluster_output/DE/C4vsCO/de_res_sig_C4vsCO.tsv', header =T, stringsAsFactors = F, sep="\t")
C1vsCO <- read.table('omics_integration/iCluster_output/DE/C1vsCO/de_res_sig_C1vsCO.tsv', header =T, stringsAsFactors = F, sep="\t")
C2vsCO <- read.table('omics_integration/iCluster_output/DE/C2vsCO/de_res_sig_C2vsCO.tsv', header =T, stringsAsFactors = F, sep="\t")
C3vsCO <- read.table('omics_integration/iCluster_output/DE/C3vsCO/de_res_sig_C3vsCO.tsv', header =T, stringsAsFactors = F, sep="\t")

## extract 
C4vsC1 <- read.table('omics_integration/iCluster_output/DE/C4vsC1/de_res_sig_C4vsC1.tsv', header =T, stringsAsFactors = F, sep="\t")
C4vsC2 <- read.table('omics_integration/iCluster_output/DE/C4vsC2/de_res_sig_C4vsC2.tsv', header =T, stringsAsFactors = F, sep="\t")
C4vsC3 <- read.table('omics_integration/iCluster_output/DE/C4vsC3/de_res_sig_C4vsC3.tsv', header =T, stringsAsFactors = F, sep="\t")
C1vsC2 <- read.table('omics_integration/iCluster_output/DE/C1vsC2/de_res_sig_C1vsC2.tsv', header =T, stringsAsFactors = F, sep="\t")
C1vsC3 <- read.table('omics_integration/iCluster_output/DE/C1vsC3/de_res_sig_C1vsC3.tsv', header =T, stringsAsFactors = F, sep="\t")
C2vsC3 <- read.table('omics_integration/iCluster_output/DE/C2vsC3/de_res_sig_C2vsC3.tsv', header =T, stringsAsFactors = F, sep="\t")

# proteomics 
C1vsCO <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_1vsCO.csv', header =T, stringsAsFactors = F, sep=",")
C2vsCO <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_2vsCO.csv', header =T, stringsAsFactors = F, sep=",")
C3vsCO <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_3vsCO.csv', header =T, stringsAsFactors = F, sep=",")
C4vsCO <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_4vsCO.csv', header =T, stringsAsFactors = F, sep=",")

# metabolomics 
C1vsCO <- read.csv('omics_integration/data/Metabolomics/DE_results_metabolomics/03-metab_effect_pval_1vsCO.csv', header =T, stringsAsFactors = F, sep=",") 
C2vsCO <- read.csv('omics_integration/data/Metabolomics/DE_results_metabolomics/03-metab_effect_pval_2vsCO.csv', header =T, stringsAsFactors = F, sep=",")
C3vsCO <- read.csv('omics_integration/data/Metabolomics/DE_results_metabolomics/03-metab_effect_pval_3vsCO.csv', header =T, stringsAsFactors = F, sep=",")
C4vsCO <- read.csv('omics_integration/data/Metabolomics/DE_results_metabolomics/03-metab_effect_pval_4vsCO.csv', header =T, stringsAsFactors = F, sep=",")

for (d in c('up', 'dn')) {
  
    C1vsCO.res <- C1vsCO$GeneName[C1vsCO$direction == d]
    C2vsCO.res <- C2vsCO$GeneName[C2vsCO$direction == d]
    C3vsCO.res <- C3vsCO$GeneName[C3vsCO$direction == d]
    C4vsCO.res <- C4vsCO$GeneName[C4vsCO$direction == d]
    
  ## for proteomics 
  if (d == "up") {
    C1vsCO.res <- unique(C1vsCO$metab[C1vsCO$effect > 0 & C1vsCO$padj < 0.05])
    C2vsCO.res <- unique(C2vsCO$metab[C2vsCO$effect > 0 & C2vsCO$padj < 0.05])
    C3vsCO.res <- unique(C3vsCO$metab[C3vsCO$effect > 0 & C3vsCO$padj < 0.05])
    C4vsCO.res <- unique(C4vsCO$metab[C4vsCO$effect > 0 & C4vsCO$padj < 0.05])
  } else {
    C1vsCO.res <- unique(C1vsCO$metab[C1vsCO$effect < 0 & C1vsCO$padj < 0.05])
    C2vsCO.res <- unique(C2vsCO$metab[C2vsCO$effect < 0 & C2vsCO$padj < 0.05])
    C3vsCO.res <- unique(C3vsCO$metab[C3vsCO$effect < 0 & C3vsCO$padj < 0.05])
    C4vsCO.res <- unique(C4vsCO$metab[C4vsCO$effect < 0 & C4vsCO$padj < 0.05])
  }
  
  cat.names4 <- c(paste0("C1vsCO\n(",length(C1vsCO.res),")") , paste0("C2vsCO\n(",length(C2vsCO.res),")") ,
                 paste0("C3vsCO\n(",length(C3vsCO.res),")"), paste0("C4vsCO\n(",length(C4vsCO.res),")"))
  
  #c4.uniq <- C4vsCO.res[!C4vsCO.res %in% c(C1vsCO.res, C2vsCO.res, C3vsCO.res)]
  #c4.uniq <- C4vsCO$GeneName[C4vsCO$GeneName %in% c4.uniq][1:100]
  #plotExpMean(d, c4.uniq)
  myCol <- brewer.pal(4, "Dark2")
  #myCol <- viridis::viridis(4)
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  vv <- venn.diagram( 
    x = list(C1vsCO.res, C2vsCO.res, C3vsCO.res, C4vsCO.res),
    category.names = cat.names4, 
    filename = NULL, main.fontface = "bold", 
    height = 600, width = 600 , resolution = 300, compression = "lzw",
    lwd = 2, lty =1, main.cex = 1.2, 
    #col=c('#fde725ff', '#21908dff', "#440154ff"),
    #fill = c(alpha('#fde725ff',1), alpha("#440154ff",1)),  ## blue=#7AA6DCFF/#003C67FF, yellow = #EFC000FF
    fill= myCol,  
    cex = 5.5, cat.cex = 3,
    fontfamily = "sans", fontface ="bold", 
    cat.default.pos = "outer",
    scale = F, 
    #cat.pos = c(-0.1, -15),
    #cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans", cat.fontface = "bold", 
    #cat.col = c('#fde725ff', "#440154ff"),  
    cat.col = myCol
    #rotation = 1
  )
  
  pdf(paste0('omics_integration/iCluster_output/Figures/Final/VennDiagramFinal/ClustersVsControl_venn_diagramm_Metab_',d,'.pdf'), width=12, height=12)
  grid.draw(vv)
  dev.off()
  
  #######################################
  ## 6-circle venn diagram 
  ######################################
  ## for pairwise comparisons 
  C4vsC1.res <- C4vsC1$GeneName[C4vsC1$direction == d]
  C4vsC2.res <- C4vsC2$GeneName[C4vsC2$direction == d]
  C4vsC3.res <- C4vsC3$GeneName[C4vsC3$direction == d]
  C1vsC2.res <- C1vsC2$GeneName[C1vsC2$direction == d]
  C1vsC3.res <- C1vsC3$GeneName[C1vsC3$direction == d] 
  C2vsC3.res <- C2vsC3$GeneName[C2vsC3$direction == d]
  cat.names6 <- c(paste0("C4vsC1\n(",length(C4vsC1.res),")"),paste0("C4vsC2\n(",length(C4vsC2.res),")"), paste0("C4vsC3\n(",length(C4vsC3.res),")"), 
                 paste0("C1vsC2\n(",length(C1vsC2.res),")"),paste0("C1vsC3\n(",length(C1vsC3.res),")"), paste0("C2vsC3\n(",length(C2vsC3.res),")"))
  data.sets6 <- list(C4vsC1.res, C4vsC2.res, C4vsC3.res, C1vsC2.res, C1vsC3.res, C2vsC3.res )
  myCol6 <- viridis::viridis(length(data.sets))
  
  pdf(paste0('omics_integration/iCluster_output/Figures/Final/BetweenClusters_venn_diagramm_',d,'.pdf'), width = 8, height = 8, useDingbats = F)
  venn(data.sets6, ilab=TRUE, zcolor = myCol6, snames= cat.names6, ilcs=1.5, sncs = 2, box =F, opacity = 0.3, lwd =4,  ggplot = F)
  dev.off()
  
}

## proteomics 6 sets venn diagram 
C4vsC1 <- read.csv('omics_integration/data/DE_results_proteomics/prot_res_4vs1.csv', header =T, stringsAsFactors = F, sep=",")
C4vsC2 <- read.csv('omics_integration/data/DE_results_proteomics/prot_res_4vs2.csv', header =T, stringsAsFactors = F, sep=",")
C4vsC3 <- read.csv('omics_integration/data/DE_results_proteomics/prot_res_4vs3.csv', header =T, stringsAsFactors = F, sep=",")
C2vsC1 <- read.csv('omics_integration/data/DE_results_proteomics/prot_res_2vs1.csv', header =T, stringsAsFactors = F, sep=",")
C3vsC1 <- read.csv('omics_integration/data/DE_results_proteomics/prot_res_3vs1.csv', header =T, stringsAsFactors = F, sep=",")
C3vsC2 <- read.csv('omics_integration/data/DE_results_proteomics/prot_res_3vs2.csv', header =T, stringsAsFactors = F, sep=",")

for (d in c('up', 'dn')) { 
  if (d == "up") {
    C4vsC1.res <- C4vsC1$EntrezGeneSymbol[C4vsC1$effect > 0 & C4vsC1$padj < 0.05]
    C4vsC2.res <- C4vsC2$EntrezGeneSymbol[C4vsC2$effect > 0 & C4vsC2$padj < 0.05]
    C4vsC3.res <- C4vsC3$EntrezGeneSymbol[C4vsC3$effect > 0 & C4vsC3$padj < 0.05]
    C2vsC1.res <- C2vsC1$EntrezGeneSymbol[C2vsC1$effect > 0 & C2vsC1$padj < 0.05]
    C3vsC1.res <- C3vsC1$EntrezGeneSymbol[C3vsC1$effect > 0 & C3vsC1$padj < 0.05]
    C3vsC2.res <- C3vsC2$EntrezGeneSymbol[C3vsC2$effect > 0 & C3vsC2$padj < 0.05]
  } else {
    C4vsC1.res <- C4vsC1$EntrezGeneSymbol[C4vsC1$effect < 0 & C4vsC1$padj < 0.05]
    C4vsC2.res <- C4vsC2$EntrezGeneSymbol[C4vsC2$effect < 0 & C4vsC2$padj < 0.05]
    C4vsC3.res <- C4vsC3$EntrezGeneSymbol[C4vsC3$effect < 0 & C4vsC3$padj < 0.05]
    C2vsC1.res <- C2vsC1$EntrezGeneSymbol[C2vsC1$effect < 0 & C2vsC1$padj < 0.05]
    C3vsC1.res <- C3vsC1$EntrezGeneSymbol[C3vsC1$effect < 0 & C3vsC1$padj < 0.05]
    C3vsC2.res <- C3vsC2$EntrezGeneSymbol[C3vsC2$effect < 0 & C3vsC2$padj < 0.05]
  }
  
  cat.names6 <- c(paste0("C4vsC1\n(",length(C4vsC1.res),")"),paste0("C4vsC2\n(",length(C4vsC2.res),")"), paste0("C4vsC3\n(",length(C4vsC3.res),")"), 
                  paste0("C2vsC1\n(",length(C2vsC1.res),")"),paste0("C3vsC1\n(",length(C3vsC1.res),")"), paste0("C3vsC2\n(",length(C3vsC2.res),")"))
  data.sets6 <- list(C4vsC1.res, C4vsC2.res, C4vsC3.res, C1vsC2.res, C1vsC3.res, C2vsC3.res )
  myCol6 <- viridis::viridis(length(data.sets))
  pdf(paste0('omics_integration/iCluster_output/Figures/Final/BetweenClusters_venn_diagramm_Proteomics_',d,'.pdf'), width = 8, height = 8, useDingbats = F)
  venn(data.sets6, ilab=TRUE, zcolor = myCol6, snames= cat.names6, ilcs=1.5, sncs = 2, box =F, opacity = 0.3, lwd =4,  ggplot = F)
  dev.off()
}

## metabolomics 6 sets venn diagram 
C4vsC1 <- read.csv('omics_integration/data/DE_results_metabolomics/03-metab_effect_pval_4vs1.csv', header =T, stringsAsFactors = F, sep=",")
C4vsC2 <- read.csv('omics_integration/data/DE_results_metabolomics/03-metab_effect_pval_4vs2.csv', header =T, stringsAsFactors = F, sep=",")
C4vsC3 <- read.csv('omics_integration/data/DE_results_metabolomics/03-metab_effect_pval_4vs3.csv', header =T, stringsAsFactors = F, sep=",")
C2vsC1 <- read.csv('omics_integration/data/DE_results_metabolomics/03-metab_effect_pval_2vs1.csv', header =T, stringsAsFactors = F, sep=",")
C3vsC1 <- read.csv('omics_integration/data/DE_results_metabolomics/03-metab_effect_pval_3vs1.csv', header =T, stringsAsFactors = F, sep=",")
C3vsC2 <- read.csv('omics_integration/data/DE_results_metabolomics/03-metab_effect_pval_3vs2.csv', header =T, stringsAsFactors = F, sep=",")

for (d in c('up', 'dn')) { 
  if (d == "up") {
    C4vsC1.res <- C4vsC1$metab[C4vsC1$effect > 0 & C4vsC1$padj < 0.05]
    C4vsC2.res <- C4vsC2$metab[C4vsC2$effect > 0 & C4vsC2$padj < 0.05]
    C4vsC3.res <- C4vsC3$metab[C4vsC3$effect > 0 & C4vsC3$padj < 0.05]
    C2vsC1.res <- C2vsC1$metab[C2vsC1$effect > 0 & C2vsC1$padj < 0.05]
    C3vsC1.res <- C3vsC1$metab[C3vsC1$effect > 0 & C3vsC1$padj < 0.05]
    C3vsC2.res <- C3vsC2$metab[C3vsC2$effect > 0 & C3vsC2$padj < 0.05]
  } else {
    C4vsC1.res <- C4vsC1$metab[C4vsC1$effect < 0 & C4vsC1$padj < 0.05]
    C4vsC2.res <- C4vsC2$metab[C4vsC2$effect < 0 & C4vsC2$padj < 0.05]
    C4vsC3.res <- C4vsC3$metab[C4vsC3$effect < 0 & C4vsC3$padj < 0.05]
    C2vsC1.res <- C2vsC1$metab[C2vsC1$effect < 0 & C2vsC1$padj < 0.05]
    C3vsC1.res <- C3vsC1$metab[C3vsC1$effect < 0 & C3vsC1$padj < 0.05]
    C3vsC2.res <- C3vsC2$metab[C3vsC2$effect < 0 & C3vsC2$padj < 0.05]
  }
  
  cat.names6 <- c(paste0("C4vsC1\n(",length(C4vsC1.res),")"),paste0("C4vsC2\n(",length(C4vsC2.res),")"), paste0("C4vsC3\n(",length(C4vsC3.res),")"), 
                  paste0("C2vsC1\n(",length(C2vsC1.res),")"),paste0("C3vsC1\n(",length(C3vsC1.res),")"), paste0("C3vsC2\n(",length(C3vsC2.res),")"))
  data.sets6 <- list(C4vsC1.res, C4vsC2.res, C4vsC3.res, C1vsC2.res, C1vsC3.res, C2vsC3.res )
  myCol6 <- viridis::viridis(length(data.sets))
  pdf(paste0('omics_integration/iCluster_output/Figures/Final/BetweenClusters_venn_diagramm_Metab_',d,'.pdf'), width = 8, height = 8, useDingbats = F)
  venn(data.sets6, ilab=TRUE, zcolor = myCol6, snames= cat.names6, ilcs=1.5, sncs = 2, box =F, opacity = 0.3, lwd =4,  ggplot = F)
  dev.off()
}

################################################################
## plot the number of unique DE for each comparison 
################################################################
dd.all <- NULL
for (d in c('up', 'dn')) {
  C1vsCO.res <- C1vsCO$GeneName[C1vsCO$direction == d]
  C2vsCO.res <- C2vsCO$GeneName[C2vsCO$direction == d]
  C3vsCO.res <- C3vsCO$GeneName[C3vsCO$direction == d]
  C4vsCO.res <- C4vsCO$GeneName[C4vsCO$direction == d]
  dd <- data.frame(dir=d, 
                   C1 = length(unique(C1vsCO.res[!C1vsCO.res %in% c(C2vsCO.res, C3vsCO.res, C4vsCO.res)]))/length(unique(C1vsCO.res)),
                   C2 = length(unique(C2vsCO.res[!C2vsCO.res %in% c(C1vsCO.res, C3vsCO.res, C4vsCO.res)]))/length(unique(C2vsCO.res)),
                   C3 = length(unique(C3vsCO.res[!C3vsCO.res %in% c(C1vsCO.res, C2vsCO.res, C4vsCO.res)]))/length(unique(C3vsCO.res)),
                   C4 = length(unique(C4vsCO.res[!C4vsCO.res %in% c(C1vsCO.res, C2vsCO.res, C3vsCO.res)]))/length(unique(C4vsCO.res)), stringsAsFactors = F)
  dd.all <- rbind(dd.all, dd)
}
dd.all$dir[dd.all$dir=="up"] <- "Up"
dd.all$dir[dd.all$dir=="dn"] <- "Down"

p <- ggplot(melt(dd.all), aes(x=variable, y=value,fill=variable)) + 
      geom_bar(stat="identity", position = "stack", color = "white", show.legend = F, width = 0.7) + 
      labs(x="", y="Percentage of DEGs") + ggtitle('') + theme_bw() + 
      geom_text(label=ifelse(melt(dd.all)$value !=0, as.character(melt(dd.all)$dir),''), position = position_stack(vjust = 0.5), hjust=0.5, color="black", size=3.5) + 
      theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.text.x = element_text( vjust= 1, size=12, face="bold", color="black"), 
        axis.text.y = element_text(size=12, face="bold"), 
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.text = element_text(size=12), 
        legend.position = "bottom") + scale_fill_brewer(palette = "Dark2") 

pdf('omics_integration/iCluster_output/Figures/Final/ClustersVsControl_Unique_DEGs_Prot.pdf', width = 3, height = 4)
print (p)
dev.off()
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------


#################################################################################################
#################### Function to plot gene expression in cell lines #################### 
#################################################################################################
plot.celllines.exp <- function (gg) {
  cell.lines <- sunshine.tech[sunshine.tech$Tissue=="Cell Line", c('Sample_Name', 'Cell.type.Tissue.type')]
  gg <- 'CXCR3'
  gg.fpkm <- melt(sunshine.fpkm[sunshine.fpkm$GeneName==gg, !colnames(sunshine.fpkm) %in% c('GeneID', 'GeneBiotype')], id.vars = "GeneName")
  gg.tpm <- melt(sunshine.tpm[sunshine.tpm$GeneName==gg, !colnames(sunshine.tpm) %in% c('GeneID', 'GeneBiotype')], id.vars = "GeneName")
  colnames(gg.fpkm) <- c('GeneName', 'Sample_Name', 'fpkm')
  colnames(gg.tpm) <- c('GeneName', 'Sample_Name', 'tpm')
  
  ## merge with cell lines information  
  gg.exp <- merge(gg.fpkm, gg.tpm)
  gg.exp <- merge(gg.exp, cell.lines)
  
  g1 = ggplot(gg.exp, aes(x=Cell.type.Tissue.type, y=log2(fpkm+1))) + geom_boxplot(aes(fill=Cell.type.Tissue.type), show.legend=F) +  
    geom_jitter(position=position_jitter(0.3), size=1) +
    labs(x="", y="Log2(FPKM+1)") + ggtitle(paste0('Expression profiles for ', gg,' in cell lines - FPKM')) + theme_bw() +
    theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
          axis.text.x = element_text(vjust= 0.5, size=12, face="bold"), 
          axis.text.y = element_text(size=12, face="bold"),
          legend.position = "none") 
  
  g2 = ggplot(gg.exp, aes(x=Cell.type.Tissue.type, y=log2(tpm+1))) + geom_boxplot(aes(fill=Cell.type.Tissue.type), show.legend=F) +  
    geom_jitter(position=position_jitter(0.3), size=1) +
    labs(x="", y="Log2(TPM+1)") + ggtitle(paste0('Expression profiles for ', gg,' in cell lines - TPM')) + theme_bw() +
    theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
          axis.text.x = element_text(vjust= 0.5, size=12, face="bold"), 
          axis.text.y = element_text(size=12, face="bold"),
          legend.position = "none") 
  
  pdf(paste0('omics_integration/iCluster_output/Figures/', gg, '_exp_profiles_cell_lines.pdf'), width=8, height = 7)
  grid.arrange(g1, g2, nrow=2)
  dev.off()
}
#################################################################################################


#################################################################################################
#################### check association with PRS #################### 
#################################################################################################
bb <- read.table("omics_integration/iCluster_output/best.cluster.memership_v1_full_AD_only.tsv", header = T, sep="\t", stringsAsFactors = F)
bb <- bb[!bb$Subj_ID %in% as.integer(gsub("ID_","", flagged.samples)), ]
prs.data <- read.table('/home/data/WashU_Data/Metabolomics_Data/03.-Phenotype/201910_Metabolon_BrainTissue_CULEBRA//metabolomics_phenoID_PRS.csv', header =T, sep=',', stringsAsFactors = T, fill=TRUE)
colnames(prs.data)[colnames(prs.data)=="short_ID"] <- "Subj_ID"
data.cols <- c('pT_5eminus8', 'pT_1eminus5', 'pT_5eminus2', 'pT_5eminus1')

## PD
prs.data.pd <- data.frame(read_xlsx('omics_integration/data/PRS_for_AD_and_PD/Umber_METAL_PD_SE_PRS.xlsx', sheet =1), stringsAsFactors = F) ## PD
colnames(prs.data.pd)[colnames(prs.data.pd)=="IID"] <- 'GWAs_ID'
prs.data.pd <- merge(prs.data.pd, prs.data[, c('GWAs_ID', 'Subj_ID')])
prs.data <- prs.data.pd
data.cols <- c('pT_0.00000005', 'pT_0.00001', 'pT_0.05', 'pT_0.5')

## merge with clusters 
prs.data <- merge(prs.data, bb[, c('Subj_ID', 'best.cluster', 'Sex', 'AOD', 'AAO', 'Status', 'APOE')])
prs.data <- prs.data[, c('Subj_ID', 'best.cluster', 'Sex', 'AOD', 'AAO', 'Status', 'APOE','GWAs_ID', data.cols)]
prs.data[prs.data$Status=="Neuro_CO", 'best.cluster'] <- 0
prs.data$best.cluster <- as.factor(prs.data$best.cluster)

## add SNCA expression levels 
snca.exp.t = reshape2::melt(rnaseq[rownames(rnaseq) == 'SNCA', ])
snca.exp.p = reshape2::melt(prt.data[rownames(prt.data) == 'SNCA', ])
colnames(snca.exp.t) <- c('Subj_ID', 'SNCA.Exp.T')
colnames(snca.exp.p) <- c('Subj_ID', 'SNCA.Exp.P') 
snca.exp.t$Subj_ID <- gsub('ID_', '', snca.exp.t$Subj_ID)
snca.exp.p$Subj_ID <- gsub('ID_', '', snca.exp.p$Subj_ID)
prs.data <- merge(prs.data, snca.exp.t)
prs.data <- merge(prs.data, snca.exp.p)

## check PRS association with SNCA expression levels 
for (c in c(1:4)) {
  cluster.data <- prs.data[prs.data$best.cluster==c, ]
  #ff <- glm(SNCA.Exp.T ~ pT_5eminus8 + Sex + AOD, data = cluster.data, family = 'gaussian')
  ff <- glm(SNCA.Exp.P ~ pT_0.00000005 + Sex + AOD, data = cluster.data, family = 'gaussian')
  #ff <- glm.nb(SNCA.Exp.T ~ pT_5eminus8 + Sex + AOD, data = cluster.data)
  cat('----------------\nCluster:', c, '\n----------------')
  print(coef(summary(ff)))
}

ff <- glm(SNCA.Exp.T ~ pT_5eminus8, data = prs.data[prs.data$best.cluster != 0,], family = 'gaussian')
ff <- glm(SNCA.Exp.T ~ pT_0.00000005, data = prs.data[prs.data$best.cluster != 0,], family = 'gaussian')
coef(summary(ff))

## all together 
ff <- glm(SNCA.Exp.T ~ pT_5eminus8 + best.cluster + Sex + AOD, data = prs.data[prs.data$best.cluster != 0,], family = 'gaussian')
ff <- glm(SNCA.Exp.T ~ 0 + pT_0.00000005 + best.cluster + Sex + AOD, data = prs.data[prs.data$best.cluster != 0,], family = 'gaussian')
#ff <- glm.nb(SNCA.Exp.T ~ pT_5eminus8 + best.cluster + Sex + AOD, data = prs.data)
coef(summary(ff))

## run GLM for cluster vs control
for (c in c(1:4)) {
  cluster.data <- prs.data[prs.data$best.cluster %in% c(c,0), ]
  #ff <- glm(pT_5eminus8 ~ best.cluster + Sex + AOD, data = cluster.data, family = 'gaussian')
  ff <- glm(pT_0.00000005 ~ best.cluster + Sex + AOD, data = cluster.data, family = 'gaussian')
  cat('----------------\nCluster:', c, '\n----------------')
  print(coef(summary(ff)))
}

##
prs.data.m <- reshape2::melt(prs.data[, c('Subj_ID', 'best.cluster', data.cols)],id.vars = c("Subj_ID","best.cluster"), variable="threshold")
prs.data.m$threshold <- gsub('minus', '-', prs.data.m$threshold)

cmp.list <- as.list(as.data.frame(combn(1: length(unique(bb$best.cluster)),2)))
prs.p <- ggplot(prs.data.m, aes(x=factor(best.cluster), y=value, group=factor(best.cluster))) + geom_boxplot(outlier.shape=NA, outlier.size=0.5, aes(color=as.factor(best.cluster)), show.legend = F) 
prs.p <- prs.p + geom_jitter(position=position_jitter(0.3), size=1.8,  aes(color=as.factor(best.cluster)), show.legend = F) + theme_bw() + facet_wrap(~ threshold, scales = "free")
prs.p <- prs.p + labs(x='', y='PRS Score') + ggtitle('Distribution of polygenic risk score (PRS) across AD clusters')
prs.p <- prs.p + theme(axis.text.x=element_text(size=14, vjust=0.5, color="black", face="bold"),
                 axis.text.y=element_text(size=14, color="black", face="bold"), 
                 axis.title.y=element_text(size=20, face="bold", color="black"),
                 plot.title = element_text(size = 16, hjust=0.5, color="black", face="bold"),
                 strip.text = element_text(size=14, face="bold"), 
                 legend.position="right", panel.border = element_rect(linetype='solid', color='black')) 
prs.p <- prs.p + scale_color_manual(values=c("0"="gray10", "1"="#1B9E77", "2"="#D95F02", "3"="#7570B3", "4"="#E7298A")) 
prs.p <- prs.p + scale_x_discrete(labels=c("1"=paste0("Cluster1\n(n=",length(bb$Subj_ID[bb$best.cluster==1 & bb$Subj_ID %in% prs.data$Subj_ID]),')'), 
                                     "2"= paste0("Cluster2\n(n=",length(bb$Subj_ID[bb$best.cluster==2 & bb$Subj_ID %in% prs.data$Subj_ID]),')'),
                                     "3"=paste0("Cluster3\n(n=",length(bb$Subj_ID[bb$best.cluster==3 & bb$Subj_ID %in% prs.data$Subj_ID]),')'), 
                                     "4"= paste0("Cluster4\n(n=",length(bb$Subj_ID[bb$best.cluster==4 & bb$Subj_ID %in% prs.data$Subj_ID]),')')))
prs.p <- prs.p + geom_signif(comparisons= cmp.list, vjust = -0.2, fontface="bold", step_increase=0.1, textsize = 4.5, size=0.4, map_signif_level=function(p)sprintf("p = %.2g", p))

pdf('omics_integration/iCluster_output/Figures/Final/polygenic_risk_score.pdf', width=12, height = 10)
prs.p
dev.off()

## check the PRS association with PD 
prs.data.pd <- data.frame(read_xlsx('omics_integration/data/PRS_for_AD_and_PD/Umber_METAL_PD_SE_PRS.xlsx', sheet =1), stringsAsFactors = F) ## PD
colnames(prs.data.pd)[colnames(prs.data.pd)=="IID"] <- 'GWAs_ID'
prs.data.pd <- merge(prs.data.pd, prs.data[, c('GWAs_ID', 'Subj_ID')])

#################################################################################################
#################### cell proportion #################### 
#################################################################################################
best.cluster.memership = read.table("omics_integration/iCluster_output/best.cluster.memership_v1_full_AD_only.tsv", header = T, sep="\t", stringsAsFactors = F)
best.cluster.memership <- best.cluster.memership[!best.cluster.memership$Subj_ID %in% as.integer(gsub("ID_","", flagged.samples)), ]

sunshine.tech$Subj_ID <- gsub("MAP_", "", sunshine.tech$Subj_ID)
sunshine.tech$Subj_ID <- gsub("MAC_", "", sunshine.tech$Subj_ID)
sunshine.tech$Subj_ID <- gsub("DIAN_", "", sunshine.tech$Subj_ID)

mend.tech$Subj_ID <- gsub("MAP_", "", mend.tech$Subj_ID)
mend.tech$Subj_ID <- gsub("MAC_", "", mend.tech$Subj_ID)
mend.tech$Subj_ID <- gsub("DIAN_", "", mend.tech$Subj_ID)

sunshine.ids <- as.integer(unique(sunshine.clin$Subj_ID[sunshine.clin$Subj_ID %in% best.cluster.memership$Subj_ID]))
mend.ids <- as.integer(unique(mend.clin$Subj_ID[mend.clin$Subj_ID %in% best.cluster.memership$Subj_ID]))

dd <- NULL
for (i in 1:nrow(best.cluster.memership)) {
  subj.id <- best.cluster.memership$Subj_ID[i]
  
  if (subj.id %in% sunshine.tech$Subj_ID) { 
      batch <- 'S'
      saample_name <- sunshine.tech$Sample_Name[sunshine.tech$Subj_ID==subj.id]
  } else if (subj.id %in% mend.tech$Subj_ID) {
      batch <- 'M'
      saample_name <- mend.tech$Sample_ID[mend.tech$Subj_ID==subj.id]
  }
  
  d <- data.frame(Subj_ID=subj.id, Sample_Name=saample_name, Batch = batch)
  dd <- rbind(dd, d)
}

## merge with best cluste solution 
best.cluster.memership <- merge(best.cluster.memership, dd, sort =F)

## plot cell-type deconvolution
algs = c( "ssFrobenius", "meanProfile")
sunshine.path <- "/home/general/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/05.-Analyses/deconvolution_analysis"
mend.path <- "/home/general/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/05.-Analyses/deconvolution_analysis"

all.res <- NULL
for (alg in algs) {
  s.res <- read.table(paste0(sunshine.path, '/deconvolution_', alg,'_results.tsv'), header = T, check.names = F, stringsAsFactors = F)
  m.res <- read.table(paste0(mend.path, '/deconvolution_', alg,'_results.tsv'), header = T, check.names = F, stringsAsFactors = F)
  all.res <- rbind(all.res, s.res)
  all.res <- rbind(all.res, m.res)
}
colnames(all.res)[colnames(all.res)=="Subject"] <- "Sample_Name"
all.res$Sample_Name <- gsub("H_VY.", "H_VY-", all.res$Sample_Name) 

# melt 
all.res <- reshape2::melt(all.res)
colnames(all.res) <- c('Sample_Name', 'Algorithm', 'Cell_type', 'Proportion')

# merge with the best clusters 
all.res <- merge(all.res, best.cluster.memership)
all.res$Subj_ID <- NULL
all.res$Batch <- NULL

## plot per cell type 
#cmp.list = as.list(as.data.frame(combn(1: length(unique(all.res$best.cluster)),2)))
## compute pvalues using GLM
all.pvals <- NULL
for (a in c( "ssFrobenius", "meanProfile")) {
  for (cell in as.character(unique(all.res$Cell_type)) ) {
    for (cc in c('1vs2', '1vs3', '1vs4', '2vs3', '2vs4', '3vs4')) {
      g1 <- as.numeric(unlist(strsplit(cc, "vs"))[1])
      g2 <- as.numeric(unlist(strsplit(cc, "vs"))[2])
      ff <- all.res[all.res$best.cluster %in% c(g1, g2) & all.res$Algorithm==a & all.res$Cell_type==cell, ]
      glm.res <- glm(Proportion ~ best.cluster + Sex + AOD, data = ff, family = 'gaussian')
      coef(summary(glm.res))[,4]
      pval <- coef(summary(glm.res))[,4][["best.cluster"]]
      dd <- data.frame(comp =cc, Algorithm=a, Cell_type = cell, PValue = signif(pval, digits = 3), best.cluster =NA, y_pos = max(ff$Proportion))
      all.pvals <- rbind(all.pvals, dd)
    }
  }
}
all.pvals$PValue <- format(all.pvals$PValue, scientific = T, digits = 3)

for (alg in c( "ssFrobenius", "meanProfile")) {
  all.plots <- list()
  for (cell in as.character(unique(all.res$Cell_type)) ) {
    dd <- all.pvals[all.pvals$Algorithm==alg & all.pvals$Cell_type==cell, ]
    p <- ggplot(all.res[all.res$Algorithm==alg & all.res$Cell_type==cell, ], aes(x=factor(best.cluster), y=Proportion, group=factor(best.cluster))) + 
         geom_boxplot(aes(fill=factor(best.cluster))) +  theme_bw() + 
         #geom_jitter(position=position_jitter(0.3),  size=1.8, aes(color=factor(best.cluster))) + theme_bw() + 
         #facet_wrap(~ Cell_type, scales = "free") + 
        scale_fill_brewer(palette = "Dark2") + labs(x='', y='') + ggtitle(cell) + 
        theme(axis.text.x=element_text(size=18, vjust=1, color="black", angle=45, hjust=1),
            #strip.text = element_text(size = 14, face="bold"), 
            #strip.background =element_rect(fill="white"),
            axis.text.y=element_text(size=18, color="black"), 
            #axis.title.y=element_text(size=12, face="bold"),
            plot.title = element_text(size = 16, hjust=0.5, color="black", face="bold"),
            legend.position="none", panel.border = element_rect(linetype='solid', color='black') ) +
        #scale_x_discrete(labels=c("1"=paste0("Cluster1\n(n=",length(unique(best.cluster.memership$Subj_ID[best.cluster.memership$best.cluster==1])),')'), 
        #                        "2"= paste0("Cluster2\n(n=",length(unique(best.cluster.memership$Subj_ID[best.cluster.memership$best.cluster==2])),')'),
        #                        "3"=paste0("Cluster3\n(n=",length(unique(best.cluster.memership$Subj_ID[best.cluster.memership$best.cluster==3])),')'), 
        #                        "4"= paste0("Cluster4\n(n=",length(unique(best.cluster.memership$Subj_ID[best.cluster.memership$best.cluster==4])),')'))) +
        scale_x_discrete(labels=c("1"="Knight-C1", "2"= "Knight-C2", "3"="Knight-C3","4"= "Knight-C4")) +
        # geom_signif(comparisons = list(c("1", "2"), c("1", "3"), c("1","4"), c("2","3"), c("2","4"), c("3","4")), 
        #           annotations = c(paste0('p=',dd$PValue[dd$comp=="1vs2"]), paste0('p=',dd$PValue[dd$comp=="1vs3"]), paste0('p=', dd$PValue[dd$comp=="1vs4"]),
        #                           paste0('p=', dd$PValue[dd$comp=="2vs3"]), paste0('p=', dd$PValue[dd$comp=="2vs4"]), paste0('p=', dd$PValue[dd$comp=="3vs4"])
        #           ), map_signif_level = TRUE, textsize=4, step_increase=0.1,  fontface="bold", margin_top=0.1, size=0.4)
        geom_signif(comparisons = list(c("1","4"), c("2","4"), c("3","4")), 
                annotations = c(paste0('p=', dd$PValue[dd$comp=="1vs4"]), paste0('p=', dd$PValue[dd$comp=="2vs4"]), paste0('p=', dd$PValue[dd$comp=="3vs4"])
                ), map_signif_level = TRUE, textsize=6, step_increase=0.1, margin_top=0.1, size=0.4)
    
    all.plots[[cell]] <- p
    
  }
  
  ## print 
  myleft <- textGrob('Cell Proportion', gp=gpar(fontsize=18, fontface="bold"), rot=90)
  mybottom <- textGrob('(Cluster1=114, Cluster2=46, Cluster3=53, Cluster4=42)', gp=gpar(fontsize=18, fontface="bold"))
  pdf(paste0("omics_integration/iCluster_output/Figures/Final/deconvolution_celllines_per_cluster_",alg,"_v2.pdf"), width = 11, height = 6.5, useDingbats = F)
  grid.arrange(grobs=all.plots, nrow=1, ncol=4, left = myleft)
  #grid.arrange(grobs=all.plots[c("Neuron", "Astrocyte")], nrow=1, ncol=2, left = myleft, bottom=mybottom)
  dev.off()
}


######### mendelian
# all.res2 <- NULL
# for (alg in algs) {
#   deconv.res2 <- read.table(paste0(path2, '/deconvolution_', alg,'_results.tsv'), header = T, check.names = F, stringsAsFactors = F)
#   all.res2 <- rbind(all.res2, deconv.res2)
# }
# colnames(all.res2)[colnames(all.res2)=="Subject"] <- "Sample_Name"
# all.res2$Sample_Name <- gsub("H_VY.", "H_VY-", all.res2$Sample_Name) 
# # merge with the best clusters 
# all.res2 <- merge(all.res2, best.cluster.memership[best.cluster.memership$Batch=="M",])
# all.res2$Subj_ID <- NULL
# all.res2$Batch <- NULL
# # melt 
# all.res2 <- melt(all.res2, id.vars = c("Sample_Name", "best.cluster", "Algorithm"))
# colnames(all.res2) <- c('Sample_Name', 'best.cluster', 'Algorithm', 'Cell_type', 'Proportion')
# 
# pdf("omics_integration/iCluster_output/Figures/deconvolution_celllines_dist_MendelianVsSporadic.pdf", width = 12, height = 8)
# ggplot(all.res2 ,aes(x = Cell_type, y = Proportion)) + 
#   geom_jitter(size=1.5, shape=1, width=0.3, aes(color=factor(best.cluster))) + theme_bw(base_size = 12) + facet_wrap( ~ Algorithm) +
#   #geom_boxplot(alpha = 0.01, lwd=0.03, position=position_dodge()) +
#   labs(x='') + guides(colour = guide_legend(override.aes = list(size = 6, stroke = 2))) + 
#   theme(axis.text.x=element_text(size=14, angle =90, vjust=0.5, face="bold"), axis.title.x=element_text(size=18, face="bold"),
#         axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
#         plot.title = element_text(size = 20, hjust=0.5, face="bold"),
#         legend.position="top",
#         axis.ticks=element_blank(), panel.grid.minor = element_blank(),
#         legend.title=element_blank(), legend.background = element_rect(color = "black"),
#         panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 
# 
# dev.off()

genotype_data <- as.data.frame(fread('omics_integration/data/Multi-omics_samples_matrix.raw'))
coi <- colnames(genotype_data)[grepl('HET', colnames(genotype_data))]
genotype_data$total_het <- genotype_data[rowSums(genotype_data[m=, coi]),]


############################################################################
## check assoication with APOE 4
############################################################################
bb <- read.table("omics_integration/iCluster_output/best.cluster.memership_v1_full_AD_only.tsv", header = T, sep="\t", stringsAsFactors = F)
bb <- bb[!bb$Subj_ID %in% gsub('ID_','',flagged.samples), ]

d = bb[!is.na(bb$TREM2), c('Subj_ID', 'best.cluster', 'TREM2')]
d = d[order(d$best.cluster), ]
write.table(d, "omics_integration/iCluster_output/TREM2_subjects_dis.tsv", sep='\t', quote = F, row.names = F)

d2 = bb[!is.na(bb$APOE), c('Subj_ID', 'best.cluster', 'APOE')]
d2 = d2[order(d2$best.cluster), ]
d2 = aggregate (Subj_ID ~ APOE + best.cluster, data=d2, length)

pdf("omics_integration/iCluster_output/Figures/distribution_of_APOE_across_cluster.pdf", width = 6, height = 5)
ggplot(d2, aes(x=best.cluster, y=Subj_ID, fill=as.factor(APOE))) + theme_bw(base_size = 12) +
  geom_bar(stat="identity")  +
  labs(x='', y = 'Number of subjects', title = 'Distribution of APOE across cluster') +
  theme(axis.text.x=element_text(size=14, vjust=0.5, face="bold"),
        axis.text.y=element_text(size=14, face="bold"),
        axis.title.y=element_text(size=16, face="bold"),
        plot.title = element_text(size = 18, hjust=0.5, face="bold"),
        axis.ticks=element_blank(), panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.background = element_rect(color = "black"),
        panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) +
  scale_fill_brewer(palette = 'Dark2' )
dev.off()

colnames(d2) <- c('APOE', 'Cluster', 'Num.Subjects') 
d2 <- d2[, c('Cluster', 'APOE', 'Num.Subjects')]  
write.table(d2, "omics_integration/iCluster_output/distribution_of_APOE_per_cluster.tsv", sep='\t', quote = F, row.names = F)

##  boxplot for each cluster 
d3 <- bb[!is.na(bb$APOE), c('Subj_ID', 'best.cluster', 'APOE')]
d3$group <- 'APOE4-'
d3[d3$APOE %in% c(24,34,44), 'group'] <- 'APOE4+'
d3 <- aggregate (Subj_ID ~ best.cluster + group, data=d3, length)
d3 <- d3[order(d3$best.cluster), ]
colnames(d3) <- c('cluster', 'grp', 'num.subjects') 


# Creating a data frame
all.pvalues <- NULL
for (clu in c(1,2,3,4)) {
  conTable = data.frame(
    "Cluster" = c(d3$num.subjects[d3$cluster==clu & d3$grp=="APOE4+"], d3$num.subjects[d3$cluster==clu & d3$grp=="APOE4-"]),
    "Others" = c(sum(d3$num.subjects[d3$cluster !=clu & d3$grp=="APOE4+"]), sum(d3$num.subjects[d3$cluster !=clu & d3$grp=="APOE4-"]))
  )
  rownames(conTable) <- c('APOE4+', 'APOE4-')
  pval <- chisq.test(conTable)$p.value
  d <- data.frame(cluster=clu, Pvalue = pval)
  all.pvalues <- rbind(all.pvalues, d)
}

p <- ggplot(d3, aes(x=as.factor(cluster), y=factor(num.subjects), fill= as.factor(grp))) + 
      geom_bar(stat="identity", width=0.4, position = position_dodge(width=0.5)) +
     labs(x="", y="Number of subjects") + ggtitle('Distribution of APOE across clusters') + theme_bw() +
     theme(plot.title = element_text(hjust=0.5, size=18, face="bold"),
        axis.text.x = element_text( vjust= 1, size=14, face="bold", color="black"), 
        axis.text.y = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=16, face="bold"), 
        legend.position = "right", legend.text = element_text(size=14, face = "bold"), legend.title=element_blank()) + scale_fill_d3() +
     scale_x_discrete(labels=c("0"=paste0("Control\n(n=",length(bb$Subj_ID[bb$best.cluster==0]),')'), 
                                 "1"=paste0("Cluster1\n(n=",length(bb$Subj_ID[bb$best.cluster==1]),')'), 
                                 "2"= paste0("Cluster2\n(n=",length(bb$Subj_ID[bb$best.cluster==2]),')'),
                                 "3"=paste0("Cluster3\n(n=",length(bb$Subj_ID[bb$best.cluster==3]),')'), 
                                 "4"= paste0("Cluster4\n(n=",length(bb$Subj_ID[bb$best.cluster==4]),')')))  

pdf("omics_integration/iCluster_output/Figures/distribution_of_APOE_across_cluster_V2.pdf", width = 6, height = 5)
p
dev.off()

########################################################
### checksynaptic genes 
########################################################
data.name <- 'Synaptome'
if (data.name == 'Synaptome') {
  synaptic_genes <- as.data.frame(read_excel("omics_integration/data/Synaptome_updated 10_31_2017.xlsx", sheet =1))  
  dir.create(paste0('omics_integration/iCluster_output/Figures/Overlap_with_', data.name)) 
} else if (data.name == 'SynGO') {
  ## from SynGo
  synaptic_genes <- as.data.frame(read_excel("omics_integration/data/SynGO_Data/SynGO_bulk_download_release_20210225/syngo_annotations.xlsx", sheet =1))
  colnames(synaptic_genes)[colnames(synaptic_genes)=="hgnc_symbol"] <- "Symbol"
  dir.create(paste0('omics_integration/iCluster_output/Figures/Overlap_with_', data.name)) 
}

comps <- c('C1vsCO','C2vsCO','C3vsCO', 'C4vsCO', 'C1vsC234', 'C2vsC134', 'C3vsC124', 'C4vsC123')
for (cmp in comps) {
  cat('Running for:', cmp, '\n')
  #dir.create(paste0('omics_integration/iCluster_output/Figures/Overlap_with_SynGO/', cmp))
  
  #"Transcriptomcs"
  C4.res <- read.table(paste0('omics_integration/iCluster_output/DE/',cmp, '/de_res_sig_',cmp,'.tsv'), header =T, stringsAsFactors = F, sep="\t")
  all.gg <- read.table(paste0('omics_integration/iCluster_output/DE/',cmp, '/de_res_',cmp,'.tsv'), header =T, stringsAsFactors = F, sep="\t")[,2]
  
  # --------------------------------------------------------------
  # proteomics 
  C4.res  <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_4vsCO.csv', header =T, stringsAsFactors = F, sep=",")
  #C4.res <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_4vsall.csv', header =T, stringsAsFactors = F, sep=",") 
  all.gg <- unique(C4.res$EntrezGeneSymbol)
  C4.res <- C4.res[C4.res$padj < 0.05,]
  colnames(C4.res)[colnames(C4.res) =="EntrezGeneSymbol"] <- 'GeneName'
  prt.gg <- unlist(strsplit(C4.res$GeneName, " "))
  prt.gg <- unlist(strsplit(prt.gg, ","))
  myCatNames <- c(paste0("ADRC (",cmp,")\n(",length(prt.gg),")"), paste0(data.name,"\n(",length(unique(synaptic_genes$Symbol)),")"))
  myCompGroup <- list(prt.gg , synaptic_genes$Symbol)
  #mytitle <- 'Up-regulated Genes'
  kk <- prt.gg
  # --------------------------------------------------------------
  
  # for (d in c('up', 'dn')) {
  #   if (d=="up") {
  #     myCatNames <- c(paste0("ADRC (",cmp,")\n(",length(C4.res$GeneName[C4.res$direction=="up"]),")") , 
  #                     paste0("SynGO\n(",length(synaptic_genes$Symbol),")"))
  #     myCompGroup <- list(C4.res$GeneName[C4.res$direction=="up"], synaptic_genes$Symbol)
  #     mytitle <- 'Up-regulated Genes'
  #     kk <- C4.res$GeneName[C4.res$direction=="up"]
  #     tt <- 'Up-regulated'
  #   } else {
      # myCatNames <- c(paste0("ADRC (",cmp,")\n(",length(C4.res$GeneName[C4.res$direction=="dn"]),")") ,
      #                 paste0("SynGO\n(",length(synaptic_genes$Symbol),")"))
      # myCompGroup <- list(C4.res$GeneName[C4.res$direction=="dn"], synaptic_genes$Symbol)
      # mytitle <- 'Down-regulated Genes'
      # kk <- C4.res$GeneName[C4.res$direction=="dn"]
      # tt <- 'Down-regulated'
  #   }
    
     myCatNames <- c(paste0("ADRC (",cmp,")\n(",length(C4.res$GeneName),")"), paste0(data.name,"\n(",length(unique(synaptic_genes$Symbol)),")"))
     myCompGroup <- list(C4.res$GeneName, synaptic_genes$Symbol)
     #mytitle <- 'Up-regulated Genes'
     kk <- C4.res$GeneName
                  
    ## write results 
    shared.genes <- C4.res[C4.res$GeneName %in% intersect(kk, synaptic_genes$Symbol), ]
    shared.genes <- shared.genes[order(shared.genes$log2FoldChange), ]
    
    ## run hypergeometirc test 
    test.mat <- matrix(c( length(unique(all.gg)) - length(union(kk, unique(synaptic_genes$Symbol))), length(setdiff(kk, unique(synaptic_genes$Symbol))), 
                          length(setdiff(unique(synaptic_genes$Symbol), kk)), length(intersect(kk, unique(synaptic_genes$Symbol)))), nrow=2)
    fisher.paval <- format(fisher.test(test.mat, alternative = "two.sided")$p.value, scientific = T, digits = 3)
    if (as.numeric(fisher.paval) ==0) { fisher.paval ='p < 2.2e-16'} else {fisher.paval = paste0('p = ',fisher.paval) }
    
    ## generate venn diagram 
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    myCol <- brewer.pal(2, "Dark2")
    vv <- venn.diagram( 
      x = myCompGroup,
      category.names = c(myCatNames),
      filename = NULL, 
      main = paste0('\n(Fisher\'s pval = ',fisher.paval,')'), 
      main.cex = 1.2, main.col = "red", main.fontface = "bold", 
      height = 300, width = 300 , resolution = 300, compression = "lzw",
      #col=c('#fde725ff', '#21908dff', "#440154ff"),
      fill = pal_jco()(2),  
      cex = 2, cat.cex = 1.2,
      fontfamily = "sans", fontface ="bold", 
      cat.default.pos = "outer",
      scale = F, lwd = 2, 
      cat.pos = c(-0.1, -0.5),
      cat.dist = c(0.055, 0.055),
      cat.fontfamily = "sans", cat.fontface = "bold", 
      cat.col = pal_jco()(2),  
      #rotation = 1
    )
    
    ## write results 
    pdf(paste0('omics_integration/iCluster_output/Figures/Overlap_with_',data.name,'/',cmp,'/Overlap_with_',data.name,'_',cmp,'.pdf'), width=5, height = 5)
    grid.draw(vv)
    dev.off()
    
#  }
}




################################################
#### DE for C4 vs ADAD and PreSymp
################################################
bb <- read.table("omics_integration/iCluster_output/best.cluster.memership_v1_full.tsv", header = T, sep="\t", stringsAsFactors = F)
bb <- bb[, c("Subj_ID", "Sex","AOD","Status", "best.cluster")]
bb[bb$Status=="Neuro_CO", 'best.cluster'] <- 0
adad_presym_clin <- clin.data[clin.data$Status %in% c('ADAD_Carrier','Neuro_Presymptomatic'), c("Subj_ID", "Sex","AOD","Status")]
adad_presym_clin$best.cluster <- 'NA'
bb <- rbind(bb, adad_presym_clin )
bb <- bb[!bb$Subj_ID %in% as.integer(gsub("ID_","", flagged.samples)), ]

## read mend 2019 technical file for ADAD
mend.tech19 <- read.csv('/40/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/03.-Phenotype/2019_09_25_WashU_MendelianVsSporadics_technical.csv', header =T, sep=",", stringsAsFactors = F, check.names = F)

## add cell lines proportion information 
algs = c( "ssFrobenius", "meanProfile")
sunshine.path <- "/home/general/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/05.-Analyses/deconvolution_analysis"
mend.path <- "/home/general/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/05.-Analyses/deconvolution_analysis"

cell.prop.res <- NULL
for (alg in algs) {
  s.res <- read.table(paste0(sunshine.path, '/deconvolution_', alg,'_results.tsv'), header = T, check.names = F, stringsAsFactors = F)
  m.res <- read.table(paste0(mend.path, '/deconvolution_', alg,'_results.tsv'), header = T, check.names = F, stringsAsFactors = F)
  cell.prop.res <- rbind(cell.prop.res, s.res)
  cell.prop.res <- rbind(cell.prop.res, m.res)
}
colnames(cell.prop.res)[colnames(cell.prop.res)=="Subject"] <- "Sample_Name"
cell.prop.res$Sample_Name <- gsub("H_VY.", "H_VY-", cell.prop.res$Sample_Name) 

## merge count data
#cell.prop.res$Subj_ID <- NA
rnase_counts <- merge(mend.cts, sunshine.cts)
for (i in 1:length(colnames(rnase_counts[,-1:-3]))) {
  col.name <- colnames(rnase_counts[,-1:-3])[i]
  subj.id <- sunshine.tech[sunshine.tech$Sample_Name==col.name, 'Subj_ID']
  if (length(subj.id) !=0) {
    colnames(rnase_counts)[colnames(rnase_counts)==col.name] <- subj.id
    cell.prop.res[cell.prop.res$Sample_Name==col.name, 'Subj_ID'] <- subj.id
  } else {
    subj.id <- mend.tech19[mend.tech19$Sample_ID==col.name, 'Subj_ID']
    colnames(rnase_counts)[colnames(rnase_counts)==col.name] <- subj.id
    cell.prop.res[cell.prop.res$Sample_Name==col.name, 'Subj_ID'] <- subj.id
  }
}

colnames(rnase_counts) <- gsub("MAP_", "", colnames(rnase_counts))
colnames(rnase_counts) <- gsub("MAC_", "", colnames(rnase_counts))
colnames(rnase_counts) <- gsub("DIAN_", "", colnames(rnase_counts))
cell.prop.res$Subj_ID <- gsub("MAP_", "", cell.prop.res$Subj_ID)
cell.prop.res$Subj_ID <- gsub("MAC_", "", cell.prop.res$Subj_ID)
cell.prop.res$Subj_ID <- gsub("DIAN_", "", cell.prop.res$Subj_ID)

## remove QC failed samples 
rnase_counts <- rnase_counts[, !colnames(rnase_counts) %in% gsub("ID_","", flagged.samples)]

## merge with cell proportions 
# get the max percentage for duplicated IDs
cell.prop.res <- cell.prop.res[cell.prop.res$Algorithm=='ssFrobenius', c('Subj_ID', 'Astrocyte', 'Neuron')]
xx = aggregate(cell.prop.res$Astrocyte, by = list(cell.prop.res$Subj_ID), max)
colnames(xx) <- c('Subj_ID', 'Astrocyte')
cell.prop.res <- merge(cell.prop.res, xx, by=c('Subj_ID', 'Astrocyte'))
bb <- merge(bb, cell.prop.res)

## make groups 
co.samples <-  as.character(bb[bb$Status=="Neuro_CO", "Subj_ID"])
c1.samples <-  as.character(bb[bb$Status=="Neuro_AD" & bb$best.cluster==1, "Subj_ID"]) 
c2.samples <-  as.character(bb[bb$Status=="Neuro_AD" & bb$best.cluster==2, "Subj_ID"]) 
c3.samples <-  as.character(bb[bb$Status=="Neuro_AD" & bb$best.cluster==3, "Subj_ID"]) 
c4.samples <-  as.character(bb[bb$Status=="Neuro_AD" & bb$best.cluster==4, "Subj_ID"]) 
adad.samples <- as.character(bb[bb$Status=="ADAD_Carrier", "Subj_ID"]) 
preysym.samples <- as.character(bb[bb$Status=="Neuro_Presymptomatic", "Subj_ID"]) 

rownames(bb) <- bb$Subj_ID
bb <- bb[c(adad.samples, preysym.samples, c4.samples), ]
bb$Status[bb$Status=="Neuro_AD"] <- 'C4'
bb$Status[bb$Status=="ADAD_Carrier"] <- 'ADAD'
bb$Status[bb$Status=="Neuro_Presymptomatic"] <- 'PreSymp'

# bb$group <- 'CO'
# bb[bb$Subj_ID %in% c1.samples, 'group'] <- 'C1'
# bb[bb$Subj_ID %in% c2.samples, 'group'] <- 'C2'
# bb[bb$Subj_ID %in% c3.samples, 'group'] <- 'C3'
# bb[bb$Subj_ID %in% c4.samples, 'group'] <- 'C4'
# 
# bb$group <- factor(bb$group, levels =c('CO', 'C1', 'C2', 'C3', 'C4'))

## order data based on groups
rownames(rnase_counts) <- paste0(rnase_counts$GeneID,'_', rnase_counts$GeneName)
rnase_counts$GeneID <- NULL
rnase_counts$GeneName <- NULL
rnase_counts$GeneBiotype <- NULL
rnase_counts <- rnase_counts[, c(adad.samples, preysym.samples, c4.samples)]
## check if sample names match 
all(rownames(bb) == colnames(rnase_counts) )

## function to filter lowly expressed genes 
rnase_counts_norm <- cpm(rnase_counts)
selectGenes <- function(grp1, grp2) {
  grp1.t <- round(length(grp1) * 0.25)
  grp2.t <- round(length(grp2) * 0.25)
  keep <- (rowSums(rnase_counts_norm[,grp1] > 0.5) >= grp1.t) | (rowSums (rnase_counts_norm[,grp2]> 0.5) >= grp2.t)
  if (all(rownames(rnase_counts) != rownames(rnase_counts_norm))) {
    stop('rownames of the normalized counts are not the same as the raw counts ')
  } else {
    cmp.counts <- rnase_counts[keep, ]
  }
  return (cmp.counts)
}

## Function for fitting the DESeq2 model
DESeqModelFit <- function(count_data) {
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = bb,
                                design = formula(~ Sex + AOD + Status))
  ## Test for DE 
  dds <- DESeq(dds)
  resultsNames(dds)
  return (dds)
}

## loop through the comparisons 
comps <- c('C4vsADAD', 'C4vsPreSymp')
for (cmp in comps) {
  cat('Running for:', cmp, '\n')
  dir.create(paste0('omics_integration/iCluster_output/DE/', cmp))
  
  ## extract groups 
  grp1.samples <- eval(parse(text = paste0(tolower(unlist(strsplit(cmp, "vs"))[1]), '.samples')))
  grp2.samples <-  eval(parse(text = paste0(tolower(unlist(strsplit(cmp, "vs"))[2]), '.samples')))   
  
  ## remove lowly expressed genes 
  grp.counts <- selectGenes(grp1.samples, grp2.samples)
  ## fit DESeq model
  dds <- DESeqModelFit(grp.counts)
  
  ## extract result table
  dds.res <- results(dds, alpha=0.05, contrast = c('Status', unlist(strsplit(cmp, "vs"))[1], unlist(strsplit(cmp, "vs"))[2]), tidy = F)
  summary(dds.res)
  
  dds.res[grepl('ENSG00000145335', rownames(dds.res)), ]
  
  ## add annotation 
  gg <- as.data.frame(str_split_fixed(rownames(dds.res), "_", 2))
  colnames(gg) <- c('GeneID', "GeneName")
  dds.res <- cbind(dds.res, gg)
  dds.res <- as.data.frame(merge(meta.cols , dds.res))
  
  ## add status 
  dds.res$direction='nc'
  dds.res[dds.res$log2FoldChange > 0 & !is.na(dds.res$padj) & dds.res$padj < 0.05, 'direction'] <- 'up'
  dds.res[dds.res$log2FoldChange < 0 & !is.na(dds.res$padj) & dds.res$padj < 0.05, 'direction'] <- 'dn'
  ## sort genes based on padj
  dds.res <- dds.res[order(dds.res$padj), ]
  
  ## extract the mean expression of each group 
  normalized_counts <- counts(dds, normalized=TRUE)
  normalized_counts = normalized_counts[paste0(dds.res$GeneID,'_', dds.res$GeneName), ]
  dds.res$grp1 <- rowMeans(normalized_counts[, grp1.samples], na.rm=TRUE)
  dds.res$grp2 <- rowMeans(normalized_counts[, grp2.samples], na.rm = TRUE)
  colnames(dds.res)[colnames(dds.res) =="grp1"] <- paste0(unlist(strsplit(cmp, "vs"))[1], '.ExpMean')
  colnames(dds.res)[colnames(dds.res) =="grp2"] <- paste0(unlist(strsplit(cmp, "vs"))[2], '.ExpMean')
  
  ## plot volcano plot
  top.10.genes <- dds.res$GeneName[1:10]
  pdf(paste0('omics_integration/iCluster_output/DE/',cmp, '/volcano_plot_',cmp,'2.pdf'), width = 6, height = 8)
  par(mfrow=c(1,1))
  with(dds.res, plot(log2FoldChange, -log10(pvalue), pch=20, main=paste0("Volcano plot - ", cmp), xlim=c(min(dds.res$log2FoldChange), max(dds.res$log2FoldChange))))
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(dds.res, direction=="up"), points(log2FoldChange, -log10(pvalue), pch=20, col="orange3"))
  with(subset(dds.res, direction=="dn") , points(log2FoldChange, -log10(pvalue), pch=20, col="skyblue3"))
  with(subset(dds.res,direction=="nc") , points(log2FoldChange, -log10(pvalue), pch=20, col="gray60"))
  with(dds.res, text(log2FoldChange, -log10(pvalue), labels=ifelse(GeneName %in% c(top.10.genes, 'CLU'), GeneName, ''), cex=0.7, offset =1, adj=c(0.5,0.01)))
  legend ("topright", c('Down', 'Up', 'NC'), col=c('skyblue3', 'orange3', 'gray60'), pch=c(20,20), pt.cex=2.5)
  dev.off()
  
  ## extract significant results  
  dds.res.sig <- dds.res[dds.res$direction !="nc", ]
  dds.res.sig <- dds.res.sig[order(dds.res.sig$padj), ]
  
  ## check with GWAS hits
  cat('Checking with GWAS hits for:', cmp, '\n')
  dds.res.sig <- add_gwas_hits (dds.res.sig)
  
  write.table(dds.res, file=paste0('omics_integration/iCluster_output/DE/',cmp,'/de_res_', cmp,'.tsv'), sep="\t", quote = F, row.names = F)
  write.table(dds.res.sig, file=paste0('omics_integration/iCluster_output/DE/',cmp, '/de_res_sig_',cmp,'.tsv'), sep="\t", quote = F, row.names = F)
  
  ## run kegg pathways 
  run_kegg_go(dds.res.sig, cmp, 0)
  run_enrichr(dds.res.sig, cmp)
  
}

#########################################################################################
#### check the levels of NFL/SNAP25 levels 
#########################################################################################
csf.markers <- as.data.frame(read_excel("AD_Staging/data/csf_markers_june22.xlsx"))
colnames(csf.markers)[colnames(csf.markers)=="ID"] <- 'Subj_ID'
csf.markers <- merge(bb[,c('Subj_ID', 'best.cluster')], csf.markers)

cc <- 'comment|Comment|COMMENT|result|CSF_LP_DATE|Subj_ID|best.cluster|STUDYID|CSF_LP_DATE|RBARCODE|RESULTU'
col.to.remove <- colnames(csf.markers)[grepl(cc, colnames(csf.markers))]

all.plts <- list()
for (mm in colnames(csf.markers)[!colnames(csf.markers) %in% col.to.remove]) {
  var.data <- csf.markers[, c('best.cluster', mm)]
  colnames(var.data) <- c('best.cluster', 'Var')
  if (all(is.na(var.data$Var))) { next }
  p <- ggplot(var.data, aes(x=factor(best.cluster), y=Var, group=factor(best.cluster))) + 
        geom_boxplot(aes(fill = factor(best.cluster)), outlier.shape=NA, outlier.size=0.5) + theme_bw() + 
        geom_jitter(position=position_jitter(0.3), size=1.5, aes(fill = factor(best.cluster))) +  
        labs(x='', y=mm) + ggtitle(mm) +
        theme(axis.text.x=element_text(size=12, vjust=1, hjust=1, color="black", angle=45),
          axis.text.y=element_text(size=12, color="black"), 
          axis.title.y=element_text(size=14, face="bold"),
          #panel.background=element_rect(color="black", size=1),
          plot.title = element_text(size = 16, hjust=0.5, color="black", face="bold"),
          legend.position="none", panel.border = element_rect(linetype='solid', color='black')) +
    scale_fill_manual(values=c("0"="gray60", "1"="#1B9E77", "2"="#D95F02", "3"="#7570B3", "4"="#E7298A")) +
    scale_x_discrete(labels=c("0"="Control", "1"="Knight-C1", "2"="Knight-C2", "3"="Knight-C3", "4"="Knight-C4")) 
  all.plts[[mm]] <- p
}

pdf('/home/eteleeb/projects/omics_integration/iCluster_output/Figures/csf_markers.pdf', width = 20, height = 30)
marrangeGrob(grobs=all.plts, ncol=4, nrow=5, left = '')
dev.off()


