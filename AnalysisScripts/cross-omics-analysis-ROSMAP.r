setwd("/home/eteleeb/projects")
#library("iCluster")
library("iClusterPlus")
library(ggplot2)
library(gridExtra)
library(parallel)
library(reshape2)
library(data.table)
library(gplots)
library(ggsignif)
library(RColorBrewer)
library(lattice)
library(enrichR)
library(DESeq2)
library(edgeR)

#### run icluster on rosmap 
dir.create('omics_integration/Replication/ROSMAP')
#dir.create(paste0('omics_integration/rosmap/', r))

###################################################
## read expression data 
###################################################
data.path <- "/home/data/Public_Data/bulkRNASeq_Data/Hg38/201506_ROSMAP/dorsolateral_prefrontal_cortex/02.-ProcessedData/Hg38"
salmon.exp <- read.table(paste0(data.path,'/04.-Salmon/combined_exp_matricies/gene_quant_matrix_TPM.tsv'), header =T, check.names = F, stringsAsFactors = F)
meta.cols <- salmon.exp[, 1:3]
fpkm.exp <- read.table(paste0(data.path,'/03.-STAR/combined_count_matrix/gene_fpkm_matrix.tsv'), header =T, check.names = F, stringsAsFactors = F)
colnames(fpkm.exp)[colnames(fpkm.exp)=="gene"] <- "GeneID"
fpkm.exp <- merge(meta.cols, fpkm.exp)
sapply(list(salmon.exp, fpkm.exp), dim)

## remove 57 samples sequenced in batch 2
batch2.samples <- scan('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201506_ROSMAP/dorsolateral_prefrontal_cortex/batch2.samples_57.tsv', 'character')
fpkm.exp <- fpkm.exp[, !colnames(fpkm.exp) %in% batch2.samples]

## generate PCA plot 
exp.4.pca <- fpkm.exp
rownames(exp.4.pca) <- exp.4.pca$GeneID
exp.4.pca$GeneID <- NULL
exp.4.pca$GeneName <- NULL
exp.4.pca$GeneBiotype <- NULL
#pdf('omics_integration/Replication/ROSMAP/PCA+plot_fpkm_without_outliers.pdf', width =10, height = 10)
#plotPCA(as.matrix(exp.4.pca[, !colnames(exp.4.pca) %in% outliers ]), col = as.numeric(groups), labels = F, cex=1.5)
#legend ("topright" ,c("B0","B1", "B2", "B3", "B4", "B5", "B6", "B8"), col=1:8, pch=c(20,20), cex=1)
#dev.off()
outliers <- c('507_120515','629_120524', '367_120502', '500_120515', '380_120503')

## read clinical and technical data 
#rosmap.clin <- read.csv('/40/Public_Data/bulkRNASeq/201506_ROSMAP/Gene_Expression/03.-Phenotype/2020_12_01_ROSMAP_clinical.csv', sep=",", header =T, stringsAsFactors = F)
rosmap.clin <- unique(read.csv('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201506_ROSMAP/03.-Phenotype/2020_12_01_ROSMAP_clinical.csv', sep=",", header =T, stringsAsFactors = F))
rosmap.tech <- read.csv('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201506_ROSMAP/03.-Phenotype/2021_06_04_ROSMAP_technical.csv', sep=",", header =T, stringsAsFactors = F)
fastQC.failed <- rosmap.tech[rosmap.tech$FastQC_Num_Fails > 4, 'Sample.Name']
#rosmap.projid.rnaseq <- read.table('/40/Public_Data/bulkRNASeq/201506_ROSMAP/Gene_Expression/03.-Phenotype/rosmap_projid_rnaseq_bam_converter.txt', sep="\t", header =T, stringsAsFactors = F, check.names = F)

# remove outliers and samples > 4 failed FastQC 
# sample "150_120419_0_merged" has no clinical data associated with
rosmap.clin <- rosmap.clin[!rosmap.clin$rnaseq_id %in% c(outliers, fastQC.failed, '150_120419_0_merged', batch2.samples), ]
fpkm.exp <- fpkm.exp[, !colnames(fpkm.exp) %in% c(outliers, fastQC.failed, '150_120419_0_merged')]
rosmap.exp <- fpkm.exp

## extract AD samples and CO samples
samples.to.keep <- rosmap.clin$rnaseq_id[rosmap.clin$cogdx %in% c(1,4) & rosmap.clin$BraakTau !=9]
ad.samples <- rosmap.clin$rnaseq_id[rosmap.clin$Status=="AD" & rosmap.clin$rnaseq_id %in% samples.to.keep]
r.co.samples <- rosmap.clin$rnaseq_id[rosmap.clin$Status=="control" & rosmap.clin$rnaseq_id %in% samples.to.keep]
pd.samples <- rosmap.clin$rnaseq_id[rosmap.clin$Status=="PD" & rosmap.clin$rnaseq_id %in% samples.to.keep]
other.samples <- rosmap.clin$rnaseq_id[rosmap.clin$Status=="other" & rosmap.clin$rnaseq_id %in% samples.to.keep]

## extract shared samples 
#ad.samples <- intersect(ad.samples, shared.samples)
#r.co.samples <- intersect(r.co.samples, shared.samples)
sapply(list(AD=ad.samples, CO=r.co.samples), length)

## ## extract AD cases only for transcriptomics 
rosmap.ad <- fpkm.exp[, colnames(fpkm.exp) %in% c('GeneID','GeneName', ad.samples)]
mean.expr = rowMeans(rosmap.ad[, -c(1,2)], na.rm = T)
rosmap.ad = rosmap.ad[order(mean.expr, decreasing=T),]
rosmap.ad = rosmap.ad[!duplicated(rosmap.ad[["GeneName"]]),]
rownames(rosmap.ad) <- rosmap.ad$GeneName
rosmap.ad$GeneID <- NULL
rosmap.ad$GeneName <- NULL
dim(rosmap.ad)

## filter lowly expressed genes 
#keep <- rowSums(rosmap.ad > 0.5) >= floor(length(colnames(rosmap.ad))*0.20)
#rosmap.ad <- rosmap.ad[keep, ]
#dim(rosmap.ad)

## replace rnaseq_id with projid
for (col in colnames(rosmap.ad)) {
  sample.sub_id <- rosmap.clin$Subj_ID[rosmap.clin$rnaseq_id == col]
  colnames(rosmap.ad)[colnames(rosmap.ad) == col] <- sample.sub_id
}

##################################################
## read metabolomics data 
###################################################
#rosmap.metab <- readRDS('omics_integration/data/ROSMAP/03-model_data_cerad_recover2.rds')
#colnames(rosmap.metab)[colnames(rosmap.metab)=="projid"] <- "Subj_ID"
#dim(rosmap.metab)
#clin.cols <- colnames(rosmap.metab)[1:20]
#metab.cols <- colnames(rosmap.metab)[!colnames(rosmap.metab) %in% clin.cols]
#rosmap.metab <- rosmap.metab[, c('Subj_ID', metab.cols)]
#rownames(rosmap.metab) <- rosmap.metab$Subj_ID
#rosmap.metab$Subj_ID <- NULL
#rosmap.metab <- as.data.frame(t(rosmap.metab))
## check the number of missing values - exclude the ones with >=20%
#rosmap.metab$na_count <- apply(rosmap.metab, 1, function(x) sum(is.na(x)))
#rosmap.metab <- rosmap.metab[rosmap.metab$na_count < length(metab.cols) * 0.20, ]
#rosmap.metab$na_count <- NULL
## replace all NA with the lowest value
#for(i in 1:nrow(rosmap.metab)){
#  indx = which(is.na(rosmap.metab[i,]))
#  r.min <- rowMeans(rosmap.metab[i,], na.rm = T)
#  rosmap.metab[i, indx] <- r.min
#}
#dim(rosmap.metab)

## extract shares samples 
#shared.samples <- intersect(colnames(rosmap.ad), colnames(rosmap.metab))
#rosmap.ad <- rosmap.ad[, shared.samples]
#rosmap.metab <-rosmap.metab[, shared.samples]
#sapply(list(rosmap.ad, rosmap.metab), dim)
#all(colnames(rosmap.ad) == colnames(rosmap.metab))

## filter features based on discovery cluster 4 + cluster DE genes  
C4vsALL.res <- read.table('omics_integration/iCluster_output/DE/C4vsC123/de_res_sig_C4vsC123.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')
C2vsALL.res <- read.table('omics_integration/iCluster_output/DE/C2vsC134/de_res_sig_C2vsC134.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')

#C4vsALL.res.prt <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_4vsall.csv', header =T, stringsAsFactors = F, check.names = F, sep=',')
#C2vsALL.res.prt <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_2vsall.csv', header =T, stringsAsFactors = F, check.names = F, sep=',')

## extract only genes DE in C4vsALL/C2vsALL for transcriptomics 
gg <- c(C4vsALL.res$GeneName, C2vsALL.res$GeneName)
rosmap.ad <- rosmap.ad[rownames(rosmap.ad) %in% gg, ]
dim(rosmap.ad)
#sapply(list(rosmap.ad, rosmap.metab), dim)

# ## select top variant features for all omics data
# feature.vars = apply(rosmap.ad, 1, var)
# threshold = feature.vars[order(feature.vars, decreasing = T)] [floor(length(feature.vars)*0.05)]
# #threshold = feature.vars[order(feature.vars, decreasing = T)] [floor(length(feature.vars)*0.20)]
# rosmap.ad <- rosmap.ad[feature.vars >= threshold, ]
# 

################################################################################
############## run tune iclusterBayes     
################################################################################
seed.val <- 4025
set.seed(seed.val) ## 357 the first best solution 
MAX.NUM.CLUSTERS = 11
num.omics <- 1
icluster.res = tune.iClusterBayes(cpus=(MAX.NUM.CLUSTERS - 1), 
                                  t(log2(rosmap.ad+1)),
                                  K=1:(MAX.NUM.CLUSTERS - 1),
                                  type=rep('gaussian', num.omics),  # poisson
                                  #n.burnin=12000, n.draw=18000,
                                  prior.gamma=rep(0.5, num.omics),
                                  pp.cutoff = 0.5,
                                  sdev=0.05,
                                  thin=3
)$fit
# save the result object
saveRDS(icluster.res, file=paste0("omics_integration/Replication/ROSMAP/iCluster.res.", seed.val,".rds"))
################################################################################

######################################
## extract dev.ratio's & BIC 
#####################################
dev.ratios = lapply(1:(MAX.NUM.CLUSTERS - 1), function(i) icluster.res[[i]]$dev.ratio)
allBICs = lapply(1:(MAX.NUM.CLUSTERS - 1), function(i) icluster.res[[i]]$BIC)

## extract best cluster 
optimal.solution = icluster.res[[which.max(dev.ratios)]] 
best.clusters = optimal.solution$clusters

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

pdf('omics_integration/Replication/ROSMAP/Figures/k_vs_devRatio.pdf', width = 12, height = 5)
grid.arrange(p1, p2, ncol=2)
dev.off()

k=which.max(dev.ratios)
plot(icluster.res[[k]]$beta.pp[[1]], xlab="Genes", ylab="Posterior probability", main="RNA-Seq expression")
plot(icluster.res[[k]]$beta.pp[[2]], xlab="Genes", ylab="Posterior probability", main="Metabolomic readings")

exp.image=scale(log2(t(rosmap.ad+1)))
exp.image[exp.image > 2.5] = 2.5
exp.image[exp.image < -2.5] = -2.5

metab.image = scale(t(rosmap.metab))
#metab.image[metab.image > 2.5 ] = 2.5
#metab.image[metab.image< -2.5 ] = --2.5

#col.scheme = alist()
col.scheme = rev(colorRampPalette(brewer.pal(10, "RdBu"))(300))
#col.scheme[[2]] = rev(colorRampPalette(brewer.pal(10, "RdBu"))(300))

####################################
#### Plot solution heatmap
###################################
pdf('omics_integration/Replication/ROSMAP/iClusterPlus_Heatmap_RNAseq.pdf', width=10, height = 2.5)
plotHMBayes(fit=optimal.solution, 
            datasets= list(exp.image), 
            type=rep("gaussian",num.omics),
            scale = rep(F,num.omics), 
            #col.scheme = col.scheme, 
            threshold=c(0.5),
            row.order=rep(T,num.omics),  
            sparse=rep(T,num.omics),
            cap=rep(T,num.omics)
)
dev.off()

########################################################
## extract cluster membership of the best class 
########################################################
all.clusters <- matrix(data= NA, nrow= nrow(t(rosmap.ad)), ncol=MAX.NUM.CLUSTERS - 1)
for (k in 1:(MAX.NUM.CLUSTERS - 1)) {
  cc = icluster.res[[k]]$clusters
  cc = as.matrix(cc)
  all.clusters[, k] = cc
}

rownames(all.clusters) = colnames(rosmap.ad)
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
write.table(best.cluster.memership, file='omics_integration/Replication/ROSMAP/best.cluster.memership.tsv', sep="\t", quote = F, row.names = F)

##### plot phenotype association 
best.cluster.memership <- read.table('omics_integration/Replication/ROSMAP/best.cluster.memership.tsv', header =T, stringsAsFactors = F, check.names = F)
cc.data <- data.frame(Subj_ID= r.co.samples, best.cluster= 0)
best.cluster.memership <- rbind(best.cluster.memership, cc.data)
colnames(best.cluster.memership) <- c('rnaseq_id', 'best.cluster')
bb = merge (best.cluster.memership, rosmap.clin)
bb$AOD = round(as.numeric(bb$AOD))
bb$AAO = round(as.numeric(bb$AAO))

## compute pvalue 
pvals.df <- NULL
aod_df <- NULL
for (cc in c('1vs2')) {
  g1 <- as.numeric(unlist(strsplit(cc, "vs"))[1])
  g2 <- as.numeric(unlist(strsplit(cc, "vs"))[2])
  ff <- bb[bb$best.cluster %in% c(g1, g2), ]
  glm.res <- glm(nft ~ best.cluster + Sex + AOD, data = ff, family = 'gaussian')
  pval <- data.frame(comp=cc, PValue= signif(coef(summary(glm.res))[,4][["best.cluster"]], digits = 3), y_pos = max(ff$nft, na.rm= T), best.cluster = NA)
  pvals.df <- rbind(pvals.df, pval)
  
  ## for AOD 
  #AOD_pval <- t.test(bb$AOD[bb$best.cluster==g1], bb$AOD[bb$best.cluster==g2])$p.value
  #AAO_pval <- t.test(bb$AAO[bb$best.cluster==g1], bb$AAO[bb$best.cluster==g2])$p.value
  #aao_d <- data.frame(comp=cc, PValue= signif(AAO_pval, digits = 3), y_pos = max(bb$AAO, na.rm= T), best.cluster = NA)
  #aod_df <- rbind(aod_df, aod_d)
}

ph.name <- 'nft'
p <- ggplot(bb[bb$best.cluster !=0,], aes(x=factor(best.cluster), y=nft, group=factor(best.cluster))) + 
      geom_boxplot(aes(color = factor(best.cluster)), outlier.shape=NA, outlier.size=0.5) +
      geom_jitter(position=position_jitter(0.3), size=1.5, aes(color = factor(best.cluster))) + theme_bw() +
      labs(x='', y='Neurofibrillary tangle burden') + ggtitle('') + scale_color_brewer(palette = "Dark2") + 
      theme(axis.text.x=element_text(size=12, vjust=0.5, face="bold", color="black"),
            axis.text.y=element_text(size=12, color="black", face="bold"),
            #axis.title.x=element_text(size=19, face="bold"),
            axis.title.y=element_text(size=14, face="bold"),
            #panel.background=element_rect(color="black"),
            plot.title = element_text(size = 14, hjust=0.5, color="black", face="bold"),
            legend.position="none", panel.border = element_rect(linetype='solid', color='black'),
            plot.margin=unit(c(1,1,1,5), 'mm')) +
      scale_x_discrete(labels=c("1"=paste0("Cluster1\n(n=",length(bb$Subj_ID[bb$best.cluster==1]),')'), "2"= paste0("Cluster2\n(n=",length(bb$Subj_ID[bb$best.cluster==2]),')'))) +
      geom_signif(data =pval, aes(xmin = 1, xmax = 2, annotations =  paste0('p=',PValue), y_position = y_pos+0.5), textsize = 3.5, step_increase=0.1, vjust = -0.2, manual = T, fontface="bold", margin_top=0.5) 

p

pdf(paste0('omics_integration/Replication/ROSMAP/phenotype_assoication_', ph.name,'.pdf'), width = 3, height = 5, useDingbats = F)
p
dev.off()

###########################################################################
##### feature selection  
###########################################################################
features = alist()
features[["rnaseq"]] = colnames(t(rosmap.ad))
features[["proteomic"]] = colnames(t(prt.data))

sigfeatures=alist()
for(i in 1:num.omics){
  rowsum=apply(abs(optimal.solution$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]] = (features[[i]])[which(rowsum>upper)]
}
names(sigfeatures)=c("expression","proteomics")

## extract features based on pp.cutoff
sig.feat <- optimal.solution$beta.pp[[1]]
sig.feat <- (features[["rnaseq"]])[which(sig.feat > 0.5)]
write.table(sig.feat, file='omics_integration/Replication/ROSMAP/sig.features.pp.rnaseq.tsv', sep="\t", quote = F, row.names = F, col.names = F)

##############################################################
##### Survival Association 
##############################################################
library(surviplot)
library(survival)

clin = bb[, c('Subj_ID', 'AAO', 'Sex', 'best.cluster')]
clin = clin[clin$AAO < 90 & !is.na(clin$AAO), ]
#clin$AOD <- as.numeric(clin$AOD)
#clin$AAO <- as.numeric(clin$AAO)
colnames(clin) = c('sample', 'time', 'Sex', 'best.cluster')

## add clusters group 
clin$group = 'NA'
clin[clin$best.cluster ==0, 'group'] <- 'CO'
clin[clin$best.cluster ==1, 'group'] <- 'C1'
clin[clin$best.cluster ==2, 'group'] <- 'C2'

max.time = max(clin$time, na.rm = T)

#groups = c('Cluster1', 'Cluster2', 'Cluster3', 'Cluster4')
groups = unique(clin$best.cluster)
cox.res = NULL
for (ii in 1:(length(groups)-1)){
  for (jj in (ii+1):length(groups)){
    cat(groups[ii], '_', groups[jj], '\n')
    
    thisy = clin[clin$group %in% c(groups[ii], groups[jj]),]
    cox = summary(coxph(Surv(time) ~ group, data=thisy))
    p = cox$sctest['pvalue']
    hr = cox$conf.int[1, c(1, 3, 4)]
    d = data.frame(comp=paste0(groups[ii], '_',groups[jj]), pval=p, HR_exp_coef=hr[1], HR_lower_0.95=hr[2], HR_upper_0.95=hr[3], stringsAsFactors = F, row.names = NULL )
    cox.res = rbind(cox.res, d)
  }
}

## extract cluster 4 
cox.res_cl4 = cox.res[grepl("Cluster4", cox.res$comp), ]

## plot KM curve 
pdf('omics_integration/Replication/ROSMAP/survival_AAO_two_groups.pdf', width=5, height = 5)
surviplot(Surv(time) ~ group, data=clin[clin$group !='CO' & clin$time <= 90, ], ylab='Survival Proportion', xlim=c(25,max.time), 
          main ="Molecular Subtype and Overall Survival", xlab='Age at onset (years)', cex.main=1,
          mark.time=TRUE, col=c("skyblue3", 'orange2'), lwd = 2.5)

dev.off()

########################################################
####### plot genes 
########################################################
gene <- 'SNCA'
rosmap.g.exp <- as.data.frame(t(fpkm.exp[fpkm.exp$GeneName==gene, -c(1:3)]))
#rosmap.g.exp <- as.data.frame(t(rosmap.ad[rownames(rosmap.ad)==gene, ]))
rosmap.g.exp$rnaseq_id <- rownames(rosmap.g.exp)
rownames(rosmap.g.exp) <- NULL
colnames(rosmap.g.exp) <- c('exp', 'rnaseq_id')
rosmap.g.exp <- merge(rosmap.g.exp, bb[, c('rnaseq_id', 'best.cluster', 'Sex', 'AOD')], sort =F)
#rosmap.g.exp$best.cluster <- factor(rosmap.g.exp$best.cluster, levels =c('0', unique(rosmap.g.exp$best.cluster)[unique(rosmap.g.exp$best.cluster) !=0]))
rosmap.g.exp = rosmap.g.exp[order(rosmap.g.exp$best.cluster),]

# ## compute pvalue 
# pvals.df <- NULL
# for (cc in c('1vs0', '2vs0', '1vs2')) {
#   g1 <- as.numeric(unlist(strsplit(cc, "vs"))[1])
#   g2 <- as.numeric(unlist(strsplit(cc, "vs"))[2])
#   ff <- rosmap.g.exp[rosmap.g.exp$best.cluster %in% c(g1, g2), ]
#   glm.res <- glm(exp ~ best.cluster + Sex + AOD, data = ff, family = 'gaussian')
#   #coef(summary(glm.res))
#   pval <- data.frame(comp=cc, PValue= signif(coef(summary(glm.res))[,4][["best.cluster"]], digits = 3), y_pos = max(log2(rosmap.g.exp$exp+1), na.rm= T), best.cluster = NA)
#   pvals.df <- rbind(pvals.df, pval)
# }

## get pvalues from DE analysis 
Tran.C1vsCO <- read.table('omics_integration/Replication/ROSMAP//DE/C1vsCO/de_res_C1vsCO.tsv', header =T, stringsAsFactors = F, sep="\t")
Tran.C2vsCO <- read.table('omics_integration/Replication/ROSMAP/DE/C2vsCO/de_res_C2vsCO.tsv', header =T, stringsAsFactors = F, sep="\t")
Tran.C1vsC2 <- read.table('omics_integration/Replication/ROSMAP/DE/C1vsC2/de_res_C1vsC2.tsv', header =T, stringsAsFactors = F, sep="\t")
Tran.ADvsCO <- read.table('omics_integration/Replication/ROSMAP/DE/ADvsCO/de_res_ADvsCO.tsv', header =T, stringsAsFactors = F, sep="\t")

pvals <- data.frame(comp=c('ADvsCO', 'C1vsCO','C2vsCO','C1vsC2'), 
                    pval = c(Tran.ADvsCO[Tran.ADvsCO$GeneName==gene, 'FDR'], Tran.C1vsCO[Tran.C1vsCO$GeneName==gene, 'padj'], 
                             Tran.C2vsCO[Tran.C2vsCO$GeneName==gene, 'padj'], Tran.C1vsC2[Tran.C1vsC2$GeneName==gene, 'padj']))
pvals$pval <- format(pvals$pval, scientific = T, digits = 3)

rosmap.g.exp$status.group <- 'CO'
rosmap.g.exp[rosmap.g.exp$best.cluster %in% c(1,2), 'status.group'] <- 'AD'
rosmap.g.exp$status.group <- factor(rosmap.g.exp$status.group, levels = c('CO', 'AD'))
rosmap.g.exp$status.group <- factor(rosmap.g.exp$status.group, levels = unique(rosmap.g.exp$status.group))

rp1 = ggplot(rosmap.g.exp, aes(x=as.factor(status.group), y=log2(exp+1))) + 
        geom_boxplot(aes(fill=as.factor(status.group)), show.legend = F) +
        labs(x="", y="Log2(FPKM+1)") + ggtitle('') + theme_bw() + 
        #geom_jitter(aes(color=as.factor(status.group)), position=position_jitter(0.3), size=1.5) +
        theme(plot.title = element_text(hjust=0.5, size=28, face="bold"),
            axis.text.x = element_text( vjust= 1, size=22,  color="black"), 
            axis.text.y = element_text(size=28, color="black"),
            axis.title.y = element_text(size=28, color="black"),
            panel.background = element_rect(colour = "black", size=1),legend.position = "none") + scale_fill_jco() + 
        scale_x_discrete(labels=c("CO"=paste0("Control\n(n=",table(rosmap.g.exp$status.group)[["CO"]], ')'), 
                            "AD"=paste0("AD\n(n=",table(rosmap.g.exp$status.group)[["AD"]],')'))) +
       geom_signif(data = pvals, aes(xmin = 1, xmax = 2, annotations =  paste0('p=',pvals$pval[pvals$comp=="ADvsCO"]), y_position = log2(max(rosmap.g.exp$exp))+0.12), 
              textsize = 7, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, color="black") 

rp2 = ggplot(rosmap.g.exp, aes(x=as.factor(best.cluster), y=log2(exp+1))) + 
      geom_boxplot(aes(fill=as.factor(best.cluster)), show.legend = F) +
      labs(x="", y="") + ggtitle(paste0(gene, '- ROSMAP')) + theme_bw() + 
      #geom_jitter(aes(color=as.factor(best.cluster)), position=position_jitter(0.3), size=1.5) +
      theme(plot.title = element_text(hjust=0.5, size=20, face="bold"),
        axis.text.x = element_text( vjust= 1, size=18, color="black"), 
        #axis.text.x = element_blank(),  
        axis.text.y = element_text(size=28, color="black"),
        axis.title.y = element_text(size=20, face="bold"),
        panel.background = element_rect(colour = "black", size=1), legend.position = "none")  +
      scale_fill_manual(values=c("0"="gray80", "1"="#E7298A", "2"="#D95F02")) +
      scale_x_discrete(labels=c("0"=paste0("Control\n(n=",length(rosmap.g.exp$rnaseq_id[rosmap.g.exp$best.cluster==0]),')'), 
                            "1"=paste0("Cluster1\n(n=",length(rosmap.g.exp$rnaseq_id[rosmap.g.exp$best.cluster==1]),')'), 
                            "2"= paste0("Cluster2\n(n=",length(rosmap.g.exp$rnaseq_id[rosmap.g.exp$best.cluster==2]),')'))) +
      geom_signif(data = pvals, aes(xmin = 1, xmax = 1.97, annotations =  paste0('p=',pvals$pval[pvals$comp=="C1vsCO"]), y_position = log2(max(rosmap.g.exp$exp))+0.15), 
              textsize = 7, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, color="black") + 
      geom_signif(data = pvals, aes(xmin = 1, xmax = 3, annotations =  paste0('p=',pvals$pval[pvals$comp=="C2vsCO"]), y_position = log2(max(rosmap.g.exp$exp))+0.35), 
              textsize = 7, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, color="black") +
      geom_signif(data = pvals, aes(xmin = 2.03, xmax = 3, annotations =  paste0('p=',pvals$pval[pvals$comp=="C1vsC2"]), y_position = log2(max(rosmap.g.exp$exp))+0.15), 
              textsize = 7, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, color="black") 

pdf(paste0('omics_integration/Replication/ROSMAP/', gene, '_exp.pdf'), width=11, height = 5.5)
grid.arrange(arrangeGrob(grobs=list(rp1, rp2), ncol=2, widths=c(1.5,3)))
dev.off()

# dd <- data.frame(samples = bb$sampleIdentifier[bb$best.cluster==2], age=bb$AOD[clin$best.cluster==2])
# dd$age <- as.numeric(gsub('\\+', '', dd$age))
# ggplot(dd, aes(x=samples, y=age)) + geom_point(size=5) + geom_hline(yintercept = 80, color="red") +
#   labs(title="Distribution of AOD for C1", xlab="samples", y='Age of Death')
# 
# t.test(clin.data$AOD[clin.data$Subj_ID %in% gsub('ID_','', c4)], clin$time[clin$best.cluster==1])

#########################################################################################
#### Check cell proportions - Cell proportion files can be found in supplementary table 2
#########################################################################################
## plot cell-type deconvolution
algs = c( "ssFrobenius", "meanProfile")
cellProp.res <- NULL
for (alg in algs) {
  s.res <- read.table(paste0('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201506_ROSMAP/dorsolateral_prefrontal_cortex/05.-Analyses/deconvolution_analysis/deconvolution_', alg,'_results.tsv'), header = T, check.names = F, stringsAsFactors = F)
  cellProp.res <- rbind(cellProp.res, s.res)
}

# merge with the best clusters 
colnames(cellProp.res)[colnames(cellProp.res)=="Subject"] <- "rnaseq_id"
cellProp.res <- merge(cellProp.res, best.cluster.memership)
#cellProp.res$sampleIdentifier <- NULL
cell.prop <- cellProp.res ## one copy for the DE analysis
cellProp.res <- cellProp.res[cellProp.res$best.cluster !=0, ]

## -----
res_for_supp <- cellProp.res
res_for_supp <- res_for_supp[res_for_supp$Algorithm=="meanProfile", -2]
res_for_supp <- res_for_supp[, c('rnaseq_id','best.cluster','Astrocyte', 'Microglia', 'Neuron', 'Oligodendrocyte')]
write.table(res_for_supp, file='omics_integration/Replication/ROSMAP/Cell_proportions_ROSMAP.tsv', sep="\t", quote = F, row.names = F)
## -----

cellProp.res <- reshape2::melt(cellProp.res, id.vars = c("rnaseq_id", "best.cluster", "Algorithm"))
colnames(cellProp.res) <- c('rnaseq_id', 'best.cluster', 'Algorithm', 'Cell_type', 'Proportion')

## compute pvalues using GLM
#colnames(cellProp.res)[colnames(cellProp.res)=="Sample_Name"] <- 'Subj_ID'
cellProp.res <- merge(cellProp.res, bb)
all.pvals <- NULL
for (a in c( "ssFrobenius", "meanProfile")) {
  for (cell in as.character(unique(cellProp.res$Cell_type)) ) {
    ff <- cellProp.res[cellProp.res$best.cluster %in% c(1,2) & cellProp.res$Algorithm==a & cellProp.res$Cell_type==cell, ]
    glm.res <- glm(Proportion ~ best.cluster + Sex + AOD, data = ff, family = 'gaussian')
    coef(summary(glm.res))[,4]
    pval <- coef(summary(glm.res))[,4][["best.cluster"]]
    dd <- data.frame(Algorithm=a, Cell_type = cell, Pvalue = signif(pval, digits = 3), best.cluster =NA, y_pos = max(ff$Proportion))
    all.pvals <- rbind(all.pvals, dd)
  }
}
all.pvals$Pvalue <- format(all.pvals$Pvalue, scientific = T, digits = 3)

for (alg in c( "ssFrobenius", "meanProfile")) {
  all.plots <- list()
  for (cell in as.character(unique(all.res$Cell_type)) ) {
    dd <- all.pvals[all.pvals$Algorithm==alg & all.pvals$Cell_type==cell, ]
    p <- ggplot(cellProp.res[cellProp.res$Algorithm==alg & cellProp.res$Cell_type==cell, ], aes(x=factor(best.cluster), y=Proportion, group=factor(best.cluster))) + 
      geom_boxplot(outlier.shape=NA, outlier.size=0.5, aes(color=factor(best.cluster))) +
      geom_jitter(position=position_jitter(0.3),  size=1.8, aes(color=factor(best.cluster))) + theme_bw() + 
      #facet_wrap(~ Cell_type, scales = "free") + 
      scale_color_brewer(palette = "Dark2") + labs(x='', y='') + ggtitle(cell) + 
      theme(axis.text.x=element_text(size=10, vjust=1, color="black", face="bold", angle=45, hjust=1),
            axis.text.y=element_text(size=12, color="black", face="bold"), 
            plot.title = element_text(size = 14, hjust=0.5, color="black", face="bold"),
            legend.position="none", panel.border = element_rect(linetype='solid', color='black')) +
      scale_x_discrete(labels=c("1"=paste0("Cluster1\n(n=",length(unique(best.cluster.memership$rnaseq_id[best.cluster.memership$best.cluster==1])),')'), 
                                "2"= paste0("Cluster2\n(n=",length(unique(best.cluster.memership$rnaseq_id[best.cluster.memership$best.cluster==2])),')'))) +
      geom_signif(data = dd, aes(xmin = 1, xmax = 2, annotations =  paste0('p=', Pvalue), y_position = y_pos), textsize = 3, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold", margin_top=0.2) 
    
    all.plots[[cell]] <- p
    
  }
  
  ## print 
  myleft <- textGrob('Cell Proportion', gp=gpar(fontsize=16, fontface="bold"), rot=90)
  pdf(paste0("omics_integration/Replication/ROSMAP/deconvolution_celllines_per_cluster_",alg,".pdf"), width = 12, height = 6, useDingbats = F)
  grid.arrange(grobs=all.plots, nrow=1, ncol=4, left = myleft)
  dev.off()
}

# ## plot per cell type 
# #cmp.list = as.list(as.data.frame(combn(1: length(unique(cellProp.res$best.cluster)),2)))
# pdf("omics_integration/Replication//ROSMAP/deconvolution_celllines_per_cluster_all.pdf", width = 12, height = 8)
# for (alg in c( "ssFrobenius", "meanProfile")) {
#   dd <- all.pvals[all.pvals$Algorithm==alg, ]
#   p <- ggplot(cellProp.res[cellProp.res$Algorithm==alg, ], aes(x=factor(best.cluster), y=Proportion, group=factor(best.cluster))) + 
#         geom_boxplot(outlier.shape=NA, outlier.size=0.5, aes(color=factor(best.cluster))) +
#         geom_jitter(position=position_jitter(0.3),  size=1.8, aes(color=factor(best.cluster))) + theme_bw() + 
#         facet_wrap(. ~ Cell_type, scales = "free") + scale_color_brewer(palette = "Dark2") +
#         labs(x='', y='Cell Proportion') + ggtitle(paste0('Overall Cell Proportion Per Cluster - (Algorithm: ', alg,')')) + 
#         theme(axis.text.x=element_text(size=11, vjust=0.5, color="black"),
#               strip.text = element_text(size = 14, face="bold"), 
#               axis.text.y=element_text(size=11, color="black"), 
#               axis.title.y=element_text(size=14, face="bold"),
#               plot.title = element_text(size = 16, hjust=0.5, color="black", face="bold"),
#               legend.position="none", panel.border = element_rect(linetype='solid', color='black')) + 
#     #geom_signif(comparisons= cmp.list, step_increase=0.1, textsize = 3.5, map_signif_level=function(p)sprintf("p = %.2g", p), show.legend=F)
#     #geom_text(data = dd, mapping = aes(x = -Inf, y = -Inf , label = paste0('p=',Pvalue)), hjust=-0.50, vjust=-0.9)
#         geom_signif(data = dd, aes(xmin = 1, xmax = 2, annotations =  paste0('p=',signif(Pvalue, digits = 3)), y_position = y_pos), textsize = 3, step_increase=0.1, vjust = -0.2,manual = T) 
#   
#   print(p)
#   
# }
# dev.off()


################################################################
##### DE Analysis  
################################################################
dir.create('omics_integration/Replication/ROSMAP/DE')

## read raw counts 
#rosmap.cts <- read.table('/home/general/Public_Data/bulkRNASeq/201506_ROSMAP/dorsolateral_prefrontal_cortex/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies/gene_quant_matrix_NumReads.tsv', header =T, check.names = F, stringsAsFactors = F)
rosmap.cts <- read.table(paste0(data.path,'/03.-STAR/combined_count_matrix/gene_count_matrix.tsv'), header =T, check.names = F, stringsAsFactors = F)
colnames(rosmap.cts)[colnames(rosmap.cts)=="gene"] <- "GeneID"
rosmap.cts <- merge(meta.cols, rosmap.cts)
rownames(rosmap.cts) <- paste0(rosmap.cts$GeneID,'_', rosmap.cts$GeneName)
rosmap.cts$GeneID <- NULL
rosmap.cts$GeneName <- NULL
rosmap.cts$GeneBiotype <- NULL

## merge with cell proportions 
# get the max percentage for duplicated IDs
cell.prop.res <- cell.prop
cell.prop.res <- cell.prop.res[cell.prop.res$Algorithm=='ssFrobenius', c('rnaseq_id', 'Astrocyte', 'Neuron')]
# xx = aggregate(cell.prop.res$Astrocyte, by = list(cell.prop.res$rnaseq_id), max)
# colnames(xx) <- c('rnaseq_id', 'Astrocyte')
# cell.prop.res <- merge(cell.prop.res, xx, by=c('rnaseq_id', 'Astrocyte'))
bb.rosmap <- merge(bb, cell.prop.res)
rownames(bb.rosmap) <- bb.rosmap$rnaseq_id

## organize data 
c1.samples <- bb.rosmap$rnaseq_id[bb.rosmap$best.cluster==1]
c2.samples <- bb.rosmap$rnaseq_id[bb.rosmap$best.cluster==2]
rosmap.cts <- rosmap.cts[, c(r.co.samples, c1.samples, c2.samples)]
bb.rosmap <- bb.rosmap[c(r.co.samples, c1.samples, c2.samples), ]
## check samples consistency 
all(rownames(bb.rosmap) == colnames(rosmap.cts) )

## function to filter lowly expressed genes 
selectGenes <- function(grp1, grp2) {
  rosmap.cts_norm <- cpm(rosmap.cts)
  grp1.t <- round(length(grp1) * 0.25)
  grp2.t <- round(length(grp2) * 0.25)
  keep <- (rowSums (rosmap.cts_norm[,grp1]> 0.5) >= grp1.t) | (rowSums (rosmap.cts_norm[,grp2]> 0.5) >= grp2.t)
  #keep <- (rowSums (rosmap.cts_norm[,grp1]> 0.5) >= grp1.t) | (rowSums (rosmap.cts_norm[,grp2]> 0.5) >= grp2.t)
  if (all(rownames(rosmap.cts) != rownames(rosmap.cts_norm))) {
    stop('rownames of the normalized counts are not the same as the raw counts ')
  } else {
    cmp.counts <- rosmap.cts[keep, ]
  }
  return (cmp.counts)
}

## Function for fitting the DESeq2 model
DESeqModelFit <- function(count_data) {
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = bb.rosmap,
                                design = formula(~ Sex + AOD + Astrocyte + Neuron + best.cluster))
  ## Test for DE 
  dds <- DESeq(dds)
  print(resultsNames(dds))
  return (dds)
}

###############################################
## enrichR
dbs <- listEnrichrDbs()
run_enrichr <- function (de, cc) {
  #gene.list <- unlist(strsplit(rownames(res_sig), "_"))
  #gene.list <- gene.list[!grepl('ENSG', gene.list)]
  en.up <- enrichr(de$GeneName[de$direction =="up"], databases = dbs$libraryName)  
  en.dn <- enrichr(de$GeneName[de$direction =="dn"], databases = dbs$libraryName)  
  save(en.up, file=paste0('omics_integration/Replication/ROSMAP/DE/', cc,'/', cc,'_enrichR_res_up.RData'))
  save(en.dn, file=paste0('omics_integration/Replication/ROSMAP/DE/', cc,'/', cc,'_enrichR_res_dn.RData'))
}
###############################################

###############################################
## edgeR 
###############################################
## make groups 
bb.rosmap$best.cluster[bb.rosmap$best.cluster==0] <- "CO"
bb.rosmap$best.cluster[bb.rosmap$best.cluster==1] <- "C1"
bb.rosmap$best.cluster[bb.rosmap$best.cluster==2] <- "C2"

groups  <- factor(c(rep('CO', length(bb.rosmap$best.cluster[bb.rosmap$best.cluster=='CO'])), 
                    rep('C1', length(bb.rosmap$best.cluster[bb.rosmap$best.cluster=='C1'])),  
                    rep('C2', length(bb.rosmap$best.cluster[bb.rosmap$best.cluster=='C2']))), levels = c('CO', 'C1', 'C2'))

#groups <- factor(c(rep('control', length(bb.rosmap$Status[bb.rosmap$Status=="control"])), rep('AD', length(bb.rosmap$Status[bb.rosmap$Status=="AD"]))), levels = c('control', 'AD'))

age <- round(as.numeric(bb.rosmap$AOD))
sex <- bb.rosmap$Sex
astrocyte <- bb.rosmap$Astrocyte
neuron <- bb.rosmap$Neuron

run_edger <- function (count_data, grp1, grp2) {
  
  if (!all(sapply(list(groups, age, sex, astrocyte, neuron), length) == ncol(count_data))) { stop('covariates do not have the same length!!')} 
  
  d <- DGEList(counts=count_data, group=factor(groups))
  d <- calcNormFactors(d)
  
  #plotMDS(d[, c(c1.samples, c2.samples)], labels = NULL, pch=20, col=c(rep(1, length(c1.samples)), rep(2, length(c2.samples))))
  #legend("bottomleft", as.character(unique(d$samples$group[d$samples$group !="CO"])), col=1:3, pch=20)
  
  design.mat <- model.matrix(~ groups + sex + age + astrocyte + neuron)
  colnames(design.mat) <- c(levels(d$samples$group), 'sex', 'age', 'astrocyte', 'neuron')
  ## for AD vs CO
  #design.mat <- model.matrix(~ 0 + groups + sex + age)
  rownames(design.mat) <- colnames(d)
  
  d2 <- estimateGLMCommonDisp(d,design.mat)
  d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
  d2 <- estimateGLMTagwiseDisp(d2,design.mat)
  
  # test 
  myContrast <- makeContrasts(eval(parse(text=grp1)) - eval(parse(text=grp2)), levels=design.mat)
  #myContrast <- makeContrasts(ADvsCO = groupsAD - groupscontrol, levels=design.mat)
  fit <- glmFit(d2, design.mat)
  lrt <- glmLRT(fit, contrast=myContrast) ## 2 and 3
  print(summary(de.number <- decideTestsDGE(lrt)))
  de.res <- topTags(lrt, n=nrow(d))$table
  de.res <-de.res[order(de.res$FDR), ]   # order by FDR
  
  return (de.res)
}
###############################################

comps <- c('C1vsCO', 'C2vsCO', 'C1vsC2')
for (cmp in comps) {
  cat('Running for:', cmp, '\n')
  dir.create(paste0('omics_integration/Replication/ROSMAP/DE/', cmp))
  
  ## extract groups 
  #mygroup <- "best.cluster"
  grp1.samples <- bb.rosmap$rnaseq_id[bb.rosmap$best.cluster== unlist(strsplit(cmp, "vs"))[1] ]
  grp2.samples <- bb.rosmap$rnaseq_id[bb.rosmap$best.cluster== unlist(strsplit(cmp, "vs"))[2] ]
  
  ## remove lowly expressed genes 
  grp.counts <- selectGenes(grp1.samples, grp2.samples)
  dim(grp.counts)
  
  ## fit DESeq model
  dds <- DESeqModelFit(grp.counts)
  ## extract result table
  dds.res <- results(dds, alpha=0.05, contrast = c('best.cluster', unlist(strsplit(cmp, "vs"))[1], unlist(strsplit(cmp, "vs"))[2] ), tidy = F)
  summary(dds.res)
  
 ## run edgeR 
 # dds.res <- run_edger(grp.counts, unlist(strsplit(cmp, "vs"))[1], unlist(strsplit(cmp, "vs"))[2])
  
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
  normalized_counts <- cpm(grp.counts)
  normalized_counts = normalized_counts[paste0(dds.res$GeneID,'_', dds.res$GeneName), ]
  dds.res$grp1 <- rowMeans(normalized_counts[, grp1.samples], na.rm=TRUE)
  dds.res$grp2 <- rowMeans(normalized_counts[, grp2.samples], na.rm = TRUE)
  colnames(dds.res)[colnames(dds.res) =="grp1"] <- paste0('grp.', unlist(strsplit(cmp, "vs"))[1], '.ExpMean')
  colnames(dds.res)[colnames(dds.res) =="grp2"] <- paste0('grp.', unlist(strsplit(cmp, "vs"))[2], '.ExpMean')
  
  ## plot volcano plot
  top.10.genes <- dds.res$GeneName[1:10]
  pdf(paste0('omics_integration/Replication/ROSMAP/DE/',cmp, '/volcano_plot_',cmp,'.pdf'), width = 6, height = 8)
  par(mfrow=c(1,1))
  with(dds.res, plot(log2FoldChange, -log10(pvalue), pch=20, main=paste0("Volcano plot - ", cmp), xlim=c(min(dds.res$log2FoldChange), max(dds.res$log2FoldChange))))
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(dds.res, direction=="up"), points(log2FoldChange, -log10(pvalue), pch=20, col="orange3"))
  with(subset(dds.res, direction=="dn") , points(log2FoldChange, -log10(pvalue), pch=20, col="skyblue3"))
  with(subset(dds.res,direction=="nc") , points(log2FoldChange, -log10(pvalue), pch=20, col="gray60"))
  with(dds.res, text(log2FoldChange, -log10(pvalue), labels=ifelse(GeneName %in% c(top.10.genes, 'CLU', 'SNCA', 'GFAP', 'APOE', 'APP'), GeneName, ''), cex=0.7, offset =1, adj=c(0.5,0.01)))
  legend ("topright", c('Down', 'Up', 'NC'), col=c('skyblue3', 'orange3', 'gray60'), pch=c(20,20), pt.cex=2.5)
  dev.off()
  
  ## extract significant results  
  dds.res.sig <- dds.res[!is.na(dds.res$padj) & dds.res$padj < 0.05, ]
  dds.res.sig <- dds.res.sig[order(dds.res.sig$padj), ]
  
  write.table(dds.res, file=paste0('omics_integration/Replication/ROSMAP/DE/',cmp,'/de_res_', cmp,'.tsv'), sep="\t", quote = F, row.names = F)
  write.table(dds.res.sig, file=paste0('omics_integration/Replication/ROSMAP/DE/',cmp, '/de_res_sig_',cmp,'.tsv'), sep="\t", quote = F, row.names = F)
  
  ## run pathway analysis
  #run_kegg_go(dds.res.sig, cmp, 0)
  run_enrichr(dds.res.sig, cmp)
  
}

#####################################################
for (d in c('up', 'dn')) {
  
  Discovery <- read.table('omics_integration/iCluster_output/DE/C4vsC123/de_res_sig_C4vsC123.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')
  Discovery <- Discovery$GeneName[Discovery$direction==d]
  
  Replicate <- read.table('omics_integration/Replication/ROSMAP/DE/C1vsC2/de_res_sig_C1vsC2.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')
  Replicate <- Replicate$GeneName[Replicate$direction==d]
  
  ## run hypergeometirc test 
  all.genes <- union(read.table('omics_integration/iCluster_output/DE/C4vsC123/de_res_C4vsC123.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')[,2], 
                     read.table('omics_integration/Replication/ROSMAP/DE/C1vsC2/de_res_C1vsC2.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')[,2])
  #all.genes <- nrow(msbb.cts)
  test.mat <- matrix(c( length(all.genes) - length(union(Discovery, Replicate)), length(setdiff(Discovery, Replicate)), length(setdiff(Replicate, Discovery)), length(intersect(Discovery, Replicate))), nrow=2)
  fisher.pval <- fisher.test(test.mat, alternative = "two.sided")$p.value
  if (fisher.pval ==0) { fisher.pval ='P-value < 2.2e-16'} else {fisher.pval = paste0('P-value = ',signif(fisher.pval, digits = 4)) }
  
  ## generate venn diagram 
  fileName <- 'KnightADRC_C4vsC123_vs_ROSMAP_C1vsC2'
  myCatNames <- c(paste0("Knight ADRC\n(",length(Discovery),")") , paste0("ROSMAP\n(",length(Replicate),")"))
  myCompGroup <- list(Discovery, Replicate)
  
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  myCol <- brewer.pal(2, "Dark2")
  vv <- venn.diagram( 
    x = myCompGroup,
    category.names = myCatNames, 
    filename = NULL, main = fisher.pval, main.fontface = "bold", 
    height = 300, width = 300 , resolution = 300, compression = "lzw",
    lwd = 2, lty = 'blank', main.cex = 1.2, 
    #col=c('#fde725ff', '#21908dff', "#440154ff"),
    #fill = c(alpha('#fde725ff',1), alpha("#440154ff",1)),  ## blue=#7AA6DCFF/#003C67FF, yellow = #EFC000FF
    fill= pal_jco()(2),  
    cex = 1.5, cat.cex = 1,
    fontfamily = "sans", fontface ="bold", 
    cat.default.pos = "outer",
    scale = F, 
    cat.pos = c(-0.1, -0.1),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans", cat.fontface = "bold", 
    #cat.col = c('#fde725ff', "#440154ff"),  
    cat.col = pal_jco()(2)
    #rotation = 1
  )
  
  pdf(paste0('omics_integration/Replication/ROSMAP/DE/', fileName, '_venn_', d,'.pdf'), width=2.8, height=2.8)
  grid.draw(vv)
  dev.off()
  
}

#################################################################
##### run pathways analysis for GO
#################################################################
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
dbs <- listEnrichrDbs()
comps <- c('C1vsCO', 'C2vsCO', 'C1vsC2')
rm(en.up); rm(en.dn)
for (cc in comps) {
  
  de.res <- read.table(paste0('omics_integration/Replication/ROSMAP/DE/',cc, '/de_res_sig_',cc,'.tsv'), header =T, stringsAsFactors = F)
  
  # en.up <- enrichr(de.res$GeneName[de.res$direction =="up"], databases = dbs$libraryName)  
  # en.dn <- enrichr(de.res$GeneName[de.res$direction =="dn"], databases = dbs$libraryName)  
  # save(en.up, file=paste0('omics_integration/Replication/MSBB/', r, '/2Clusters/DE/', cc,'/', cc,'_enrichR_res_up.RData'))
  # save(en.dn, file=paste0('omics_integration/Replication/MSBB/',r,'/2Clusters/DE/', cc,'/', cc,'_enrichR_res_dn.RData'))
  
  #using gage 
  # patheay analysis using gage
  de.res$entrez <- mapIds(org.Hs.eg.db, keys=de.res$GeneName, column="ENTREZID",keytype="SYMBOL",multiVals="first")
  foldchanges <- de.res$log2FoldChange
  names(foldchanges) <- de.res$entrez
  go_bp_res <- gage(foldchanges, gsets=go.bp.gs, same.dir=TRUE)
  go_bp_up <- data.frame(id=rownames(go_bp_res$greater), go_bp_res$greater)
  go_bp_dn <- data.frame(id=rownames(go_bp_res$less), go_bp_res$less)
  write.table(go_bp_up, file=paste0('omics_integration/Replication/ROSMAP/DE/', cc,'/',cc,'_GO_BP_UP.tsv'), sep="\t", quote = F, row.names = F)
  write.table(go_bp_dn, file=paste0('omics_integration/Replication/ROSMAP/DE/', cc,'/',cc,'_GO_BP_DN.tsv'), sep="\t", quote = F, row.names = F)

 }

#################################################################
##### check homeostatic and actived mic genes 
#################################################################
## function to make stat table
make_stat_table <- function (g) {
  for (cmp in  c('C1vsCO', 'C2vsCO', 'C1vsC2')) {
    res <- read.table(paste0('omics_integration/Replication/ROSMAP/DE/',cmp,'/de_res_',cmp,'.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')
    res <- res[res$GeneName==g, c('log2FoldChange', 'pvalue', 'padj')]
    colnames(res) <- c('LogFC', 'PValue', 'FDR')
    res <- reshape2::melt(res)
    colnames(res) <- c('Data', cmp)
    res[,2] <- signif(res[,2], digits = 3)
    if (cmp == 'C1vsCO') {
      rr <- res 
    } else {
      rr <- merge(rr, res, sort =F )
    }
  }
  
  if (nrow(rr) !=0) {
    return(rr)
  }
  
}

mic.genes <- c('APOE', 'SPP1', 'TREM2', 'CXCR3', 'TMEM119', 'BIN1', 'MED12L', 'SELPLG', 'ABCA1', 'IL1B', 'CD68', 'RELB',
               'CSF1', 'ITGAM', 'P2RY12', 'TYROBP', 'AIF1', 'AXL', 'LAG3', 'SPI1', 'IRF8', 'CIITA', 'IFI16', 'ITGB2', 'TLR2', 'CX3CR1', 'IL1A')

homeostatic <- c('BIN1', 'CX3CR1', 'MED12L', 'SELPLG', 'P2RY12', 'P2RY13', 'TMEM119')
#bruno_genes_of_interest <- toupper(c('Spp1', 'Itgax', 'Axl', 'Lilrb4', 'Clec7a', 'Csf1','Apoe'))
bruno_homeostatic_genes <- toupper(c('P2ry12', 'Tmem119', 'Gpr34', 'Jun', 'Olfml3', 'Csf1r', 'Hexb', 'Mertk', 'Rhob', 'Cx3Cr1', 'Tgfbr1', 'Tgfb1'))
homeostatic <- unique(c(homeostatic, bruno_homeostatic_genes))
activated <- c('ABCA1', 'C5AR1', 'FCGR2B', 'GPNMB', 'CD68', 'CD83', 'JUN', 'LGMN', 'LPL', 'TNFAIP3')

## plot genes 
pdf('omics_integration/Replication/ROSMAP/mic.genes.all.pdf', width=6, height = 4)
for (gg in mic.genes) {
  
  gg.exp <- reshape2::melt(fpkm.exp[fpkm.exp$GeneName==gg, !colnames(fpkm.exp) %in% c('GeneID', 'GeneBiotype')], id.vars = "GeneName")
  colnames(gg.exp) <- c('GeneName', 'rnaseq_id', 'fpkm')
  
  # merge with the best clusters 
  gg.exp <- merge(gg.exp, bb)
  
  ## compute pvalue 
  pvals.df <- NULL
  for (cc in c('1vs0', '2vs0', '1vs2')) {
    g1 <- as.numeric(unlist(strsplit(cc, "vs"))[1])
    g2 <- as.numeric(unlist(strsplit(cc, "vs"))[2])
    ff <- gg.exp[gg.exp$best.cluster %in% c(g1, g2), ]
    glm.res <- glm(fpkm ~ best.cluster + Sex + AOD, data = ff, family = 'gaussian')
    pval <- data.frame(comp=cc, PValue= signif(coef(summary(glm.res))[,4][["best.cluster"]], digits = 3), y_pos = max(log2(ff$fpkm+1), na.rm= T), best.cluster = NA)
    pvals.df <- rbind(pvals.df, pval)
  }
  
  g1 = ggplot(gg.exp, aes(x=as.factor(best.cluster), y=log2(fpkm+1))) + geom_boxplot(aes(fill=as.factor(best.cluster)), show.legend=F) +  
      geom_jitter(position=position_jitter(0.3), size=1) +
      labs(x="", y="Log2(FPKM+1)") + ggtitle(paste0('RNA-Seq expression profiles for ', gg)) + theme_bw() +
      theme(plot.title = element_text(hjust=0.5, size=13, face="bold"),
          axis.text.x = element_text(vjust= 0.5, size=12, face="bold"), 
          axis.text.y = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"), 
          legend.position = "none") + scale_fill_brewer(palette = 'Dark2' ) +  
      geom_signif(comparisons = list(c("0", "1"), c("0", "2"), c("1","2")), 
              annotations = c(paste0('p=',pvals.df$PValue[pvals.df$comp=="1vs0"]), paste0('p=',pvals.df$PValue[pvals.df$comp=="2vs0"]), paste0('p=', pvals.df$PValue[pvals.df$comp=="1vs2"])), 
              map_signif_level = TRUE, textsize=3, step_increase=0.1)
  print(g1)
  
}
dev.off()

################################################################################
## prepare data for GSEA 
dir.create('omics_integration/data/4GSEA/ROSMAP')
res <- read.table('omics_integration/Replication/ROSMAP/DE/C1vsC2/de_res_sig_C1vsC2.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')
res <- fpkm.exp[fpkm.exp$GeneName %in% res$GeneName, c('GeneName', bb.rosmap$rnaseq_id[bb.rosmap$best.cluster=='C1'], bb.rosmap$rnaseq_id[bb.rosmap$best.cluster=='C2'])]
res$DESCRIPTION = NA
res <- res[, c('GeneName', 'DESCRIPTION', bb.rosmap$rnaseq_id[bb.rosmap$best.cluster=='C1'], bb.rosmap$rnaseq_id[bb.rosmap$best.cluster=='C2'])]
colnames(res)[colnames(res)=="GeneName"] <- 'NAME'
write.table(res, file="omics_integration/data/4GSEA/ROSMAP/exp_data_C1vsC2_rosmap.gct", sep="\t", quote = F, row.names = F)

## make the phenotype file 
ph <- c(rep('C1', length(bb.rosmap$rnaseq_id[bb.rosmap$best.cluster=='C1']) ), rep('C2', length(bb.rosmap$rnaseq_id[bb.rosmap$best.cluster=='C2'])) )
ph <- paste(ph, collapse = " ")
write.table(ph, file="omics_integration/data/4GSEA/ROSMAP/phenotype_class_C1vsC2_rosmap.cls", sep=" ", quote = F, row.names = F)
#################################################################################

## check with chip-seq data 
chip.seq <- as.data.frame(read_excel("omics_integration/data/ROSMAP/41593_2018_291_MOESM3_ESM.xls", sheet = 1))[-1,c(1,20)]
chip.seq <- chip.seq[as.integer(chip.seq$`Sample Information`) %in% bb$Subj_ID, ]
colnames(chip.seq) <- c('Subj_ID', 'Cross_Corr')
chip.seq$Subj_ID <- as.integer(chip.seq$Subj_ID)
chip.seq$Cross_Corr <- signif(as.numeric(chip.seq$Cross_Corr),digits = 3)
chip.seq <- merge(chip.seq, bb[, c('Subj_ID', 'best.cluster')])

ggplot(chip.seq, aes(x=factor(best.cluster), y=Cross_Corr, group=factor(best.cluster))) + geom_boxplot(outlier.shape=NA, outlier.size=0.5) +
      geom_jitter(position=position_jitter(0.3), size=2) + theme_bw() + # + guides(color = FALSE)  - aes(color=as.factor(cogdx))
      labs(x='', y='') + ggtitle(paste0('Cross Correlation')) +
      theme(axis.text.x=element_text(size=11, vjust=0.5, color="black"),
                 axis.text.y=element_text(size=12, color="black"), 
                 axis.title.y=element_text(size=12),
                 #panel.background=element_rect(color="black"),
                 plot.title = element_text(size = 14, hjust=0.5, color="black", face="bold"),
                 legend.position="right", 
                 panel.border = element_rect(linetype='solid', color='black'),
                 plot.margin=unit(c(1,1,1,5), 'mm'))+ scale_color_brewer(name="APOE", palette = "Dark2") +
      scale_y_continuous(breaks = seq(min(chip.seq$Cross_Corr), max(chip.seq$Cross_Corr), 0.05))


##############################################################
####### Run DE for Metabolomics  
##############################################################
dir.create('omics_integration/Replication/ROSMAP/Metab_DE')
dir.create('omics_integration/Replication/ROSMAP/Metab_DE/1vs0')
dir.create('omics_integration/Replication/ROSMAP/Metab_DE/2vs0')
dir.create('omics_integration/Replication/ROSMAP/Metab_DE/1vs2')
dir.create('omics_integration/Replication/ROSMAP/Metab_DE/2vs1')

## read cluster information
best.cluster.memership <- read.table('omics_integration/Replication/ROSMAP/best.cluster.memership.tsv', header =T, stringsAsFactors = F, check.names = F)
cc.data <- data.frame(Subj_ID= r.co.samples, best.cluster= 0)
best.cluster.memership <- rbind(best.cluster.memership, cc.data)
colnames(best.cluster.memership) <- c('rnaseq_id', 'best.cluster')
bb = merge (best.cluster.memership, rosmap.clin)
bb$AOD = round(as.numeric(bb$AOD))
bb$AAO = round(as.numeric(bb$AAO))

## read rosmap metabolomics 
rosmap.metab <- readRDS('omics_integration/data/ROSMAP/03-model_data_cerad_recover2.rds')
colnames(rosmap.metab)[colnames(rosmap.metab)=="projid"] <- "Subj_ID"
rosmap.metab <- merge(rosmap.metab, bb[, c('Subj_ID', 'best.cluster')])
dim(rosmap.metab)

clin.cols <- c(colnames(rosmap.metab)[1:20], 'best.cluster')
metab.cols <- colnames(rosmap.metab)[!colnames(rosmap.metab) %in% clin.cols]

## handle missing values 
rosmap.metab$na_count <- apply(rosmap.metab[, metab.cols], 1, function(x) sum(is.na(x)))
rosmap.metab <- rosmap.metab[rosmap.metab$na_count < (length(metab.cols) * 0.20), ]
rosmap.metab$na_count <- NULL
## replace all NA with the lowest value
for(mm in metab.cols){
  indx = which(is.na(rosmap.metab[, mm]))
  r.min <- mean(rosmap.metab[, mm], na.rm = T)
  rosmap.metab[indx, mm] <- r.min
}

## run GLM model
#reading ~ best.cluster + Age + Sex + PMI
for (cc in c('1vs0', '2vs0', '1vs2')) {
  g1 <- as.numeric(unlist(strsplit(cc, "vs"))[1])
  g2 <- as.numeric(unlist(strsplit(cc, "vs"))[2])
  ff <- rosmap.metab[rosmap.metab$best.cluster %in% c(g1, g2), ]
  
  ff$best.cluster[ff$best.cluster==g1] <- 1
  ff$best.cluster[ff$best.cluster==g2] <- 0
  
  ## loop through metabs
  mm.res.all <- NULL
  for (mm in metab.cols) {
    mm.reading <- ff[, c('best.cluster', 'age_death', 'msex', 'pmi', mm)]
    colnames(mm.reading) <- c('best.cluster', 'age_death', 'msex', 'pmi', 'reading')
    glm.res <- glm(reading ~ best.cluster + msex + age_death + pmi, data = mm.reading, family = 'gaussian')
    mm.res= as.data.frame(coef(summary(glm.res)))[rownames(as.data.frame(coef(summary(glm.res))))=="best.cluster",]
    colnames(mm.res) <- c('estimate', 'std.error','t.value', 'p.value')
    rownames(mm.res) <- mm
    mm.res <- signif(mm.res, digits = 3)
    mm.res <- data.frame(mm.res, stringsAsFactors = F)
    mm.res.all <- rbind(mm.res.all, mm.res)
  }
  ## write results 
  ## correct for false discovery 
  mm.res.all <- cbind(mm.res.all, padj = p.adjust(mm.res.all$p.value, method='BH'))
  write.table(mm.res.all, file=paste0('omics_integration/Replication/ROSMAP/Metab_DE/',cc, '/', cc, '.res.all.tsv'), sep="\t", quote = F, row.names = T)
  
}


########################################################
### check synaptic genes 
########################################################
data.name <- 'Synaptome'
if (data.name == 'Synaptome') {
  synaptic_genes <- as.data.frame(read_excel("omics_integration/data/Synaptome_updated 10_31_2017.xlsx", sheet =1))  
  dir.create(paste0('omics_integration/Replication/ROSMAP//Overlap_with_', data.name)) 
} else if (data.name == 'SynGO') {
  ## from SynGo
  synaptic_genes <- as.data.frame(read_excel("omics_integration/data/SynGO_Data/SynGO_bulk_download_release_20210225/syngo_annotations.xlsx", sheet =1))
  colnames(synaptic_genes)[colnames(synaptic_genes)=="hgnc_symbol"] <- "Symbol"
  dir.create(paste0('omics_integration/Replication/ROSMAP/Overlap_with_', data.name)) 
}

comps <- c('C1vsCO','C2vsCO','C1vsC2')
for (cmp in comps) {
  cat('Running for:', cmp, '\n')
  #dir.create(paste0('omics_integration/iCluster_output/Figures/Overlap_with_SynGO/', cmp))
  
  #"Transcriptomics"
  C4.res <- read.table(paste0('omics_integration/Replication/ROSMAP/DE/',cmp,'/de_res_sig_',cmp,'.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')
  all.gg <- read.table(paste0('omics_integration/Replication/ROSMAP/DE/',cmp,'/de_res_',cmp,'.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')[,2]
  
  # --------------------------------------------------------------
  # proteomics 
  # C4.res  <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_4vsCO.csv', header =T, stringsAsFactors = F, sep=",")
  # #C4.res <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_4vsall.csv', header =T, stringsAsFactors = F, sep=",") 
  # all.gg <- unique(C4.res$EntrezGeneSymbol)
  # C4.res <- C4.res[C4.res$padj < 0.05,]
  # colnames(C4.res)[colnames(C4.res) =="EntrezGeneSymbol"] <- 'GeneName'
  # prt.gg <- unlist(strsplit(C4.res$GeneName, " "))
  # prt.gg <- unlist(strsplit(prt.gg, ","))
  # myCatNames <- c(paste0("ADRC (",cmp,")\n(",length(prt.gg),")"), paste0(data.name,"\n(",length(unique(synaptic_genes$Symbol)),")"))
  # myCompGroup <- list(prt.gg , synaptic_genes$Symbol)
  # #mytitle <- 'Up-regulated Genes'
  # kk <- prt.gg
  # --------------------------------------------------------------
  
  myCatNames <- c(paste0("ROSMAP (",cmp,")\n(",length(C4.res$GeneName),")"), paste0(data.name,"\n(",length(unique(synaptic_genes$Symbol)),")"))
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
  pdf(paste0('omics_integration/Replication/ROSMAP/Overlap_with_',data.name,'/Overlap_with_',data.name,'_',cmp,'.pdf'), width=5, height = 5)
  grid.draw(vv)
  dev.off()
  
  #  }
}


