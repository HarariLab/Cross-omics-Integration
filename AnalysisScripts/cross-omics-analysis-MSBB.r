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
library(lattice)

#### run icluster on MSBB region BM36
r <- 'BM36'
dir.create('omics_integration/MSBB')
dir.create(paste0('omics_integration/MSBB/', r))

## read expression data 
msbb.tpm <- read.table(paste0('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201812_MSBB/',r, '/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies/gene_quant_matrix_TPM.tsv'), header =T, check.names = F, stringsAsFactors = F)
msbb.meta.cols <- msbb.tpm[, 1:3]
msbb.fpkm <- read.table(paste0('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201812_MSBB/',r, '/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_fpkm_matrix.tsv'), header =T, check.names = F, stringsAsFactors = F)
colnames(msbb.fpkm)[colnames(msbb.fpkm)=="gene"] <- "GeneID"
msbb <- merge(msbb.meta.cols, msbb.fpkm)

## read clinical and technical data 
msbb.clin <- read.csv('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201812_MSBB/03.-Phenotype/2020_03_25_MSBB_clinical.csv', sep=",", header =T, stringsAsFactors = F)
msbb.tech <- read.csv('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201812_MSBB/03.-Phenotype/2020_10_29_MSBB_technical.csv', sep=",", header =T, stringsAsFactors = F)
## remove FLAGed samples and extract region data only 
#msbb.tech <- msbb.tech[msbb.tech$Tissue==r & msbb.tech$FLAG !='FLAG' & !msbb.tech$Action %in% c('Exclude', 'Remap') & msbb.tech$PCA_outlier !="Yes", ]
msbb.tech <- msbb.tech[msbb.tech$Tissue==r & msbb.tech$FLAG !='FLAG' & !msbb.tech$Action %in% c('Exclude', 'Remap'), ]

### extract unique individuals 
subjects_to_keep <- c()
individualIDs <- unique(msbb.tech$individualID)
for (j in 1:length(individualIDs)) {
  individualID <- individualIDs[j]
  sample.id <- msbb.tech$sampleIdentifier[msbb.tech$individualID==individualID]
  s.exp <- msbb[, c('GeneID', 'GeneName', sample.id)]
  
  if (length(sample.id) > 1) {
    if (any(grepl('resequenced', sample.id))) {
      sample.id <- sample.id[grepl('resequenced', sample.id)]
      s.exp <- msbb[, c('GeneID', 'GeneName', sample.id)]
      subjects_to_keep <- c(subjects_to_keep, sample.id)
    } else {
      #s.exp$mean.exp <- rowMeans(s.exp[, -1:-2])
      #s.exp <- s.exp[, c('GeneID', 'GeneName', 'mean.exp')]   
      ss <- colSums(s.exp[,-1:-2])
      max.val <- order(ss, decreasing = T)
      rep.sample <- names(ss)[max.val[1]]
      s.exp <- s.exp[, c('GeneID', 'GeneName', rep.sample)] 
      subjects_to_keep <- c(subjects_to_keep, rep.sample)
    }
  } else {
    subjects_to_keep <- c(subjects_to_keep, sample.id)
  }
  
  ## make the MAPID as column name 
  colnames(s.exp) <- c('GeneID', 'GeneName', individualID)
  
  ## merge all 
  if (j ==1 ) {
    msbb.exp <- s.exp
  } else {
    msbb.exp <- merge(msbb.exp, s.exp, sort = F)
  }
}

## move data to msbb object 
msbb <- msbb.exp

## ------------------------------------------------------------------
## read proteomic data from BM36
prt.data <- read.table('omics_integration/data/MSBB_BM36/2a.MSBB_PHG-unregressed_log2(ratio)_9395x150.txt', header =T, check.names = F, stringsAsFactors = F)
prt.pheno <- read.csv('omics_integration/data/MSBB_BM36/0.MSBB_PHG-Traits.csv', sep=",", header =T, check.names = F, stringsAsFactors = F)
cols <- colnames(prt.data)
for (k in 1:length(cols)){
  col.id <- cols[k]
  ind.id <- prt.pheno$SampleID[prt.pheno$batch.channel==col.id]
  colnames(prt.data)[colnames(prt.data)==col.id] <- ind.id
}
sapply(list(msbb, prt.data), dim)

## check missing values 
prt.data$na_count <- apply(prt.data, 1, function(x) sum(is.na(x)))
prt.data <- prt.data[prt.data$na_count < length(colnames(prt.data)[colnames(prt.data) !='na_count']) * 0.20, ]
prt.data$na_count <- NULL
## replace all NA with the average value
for(i in 1:nrow(prt.data)){
  indx = which(is.na(prt.data[i,]))
  row.val <- rowMeans(prt.data[i,], na.rm = T)
  prt.data[i, indx] <- row.val
}
prt.exp <- prt.data
## ------------------------------------------------------------------

## extract shared samples 
shared.samples <- intersect(colnames(msbb), colnames(prt.data))

## extract AD samples and CO samples
ad.samples <- msbb.clin$Subj_ID[msbb.clin$Status=="AD"]
co.samples <- msbb.clin$Subj_ID[msbb.clin$Status=="control"]
## extract shared samples 
ad.samples <- intersect(ad.samples, shared.samples)
co.samples <- intersect(co.samples, shared.samples)
sapply(list(ad.samples, co.samples), length)
msbb.co.samples <- co.samples

## ## extract AD cases only for transcriptomics 
msbb.ad <- msbb[, colnames(msbb) %in% c('GeneID','GeneName', ad.samples)]
mean.expr = rowMeans(msbb.ad[, -c(1,2)], na.rm = T)
msbb.ad = msbb.ad[order(mean.expr, decreasing=T),]
msbb.ad = msbb.ad[!duplicated(msbb.ad[["GeneName"]]),]
rownames(msbb.ad) <- msbb.ad$GeneName
msbb.ad$GeneID <- NULL
msbb.ad$GeneName <- NULL

## filter lowly expressed genes 
#keep <- rowSums(msbb.ad > 0.5) >= floor(length(colnames(msbb.ad))*0.20)
#msbb.ad <- msbb.ad[keep, ]
#dim(msbb.ad)

## extract AD cases only for proteomics
prt.data <- prt.data[, ad.samples]
sapply(list(msbb.ad, prt.data), dim)

## use iCluster features 
#iCluster.feat <- read.table('omics_integration/V1_4Clusters/exp_top_genes_no_co.tsv', header =T, stringsAsFactors = F, check.names = F)
#msbb.ad <- msbb.ad[rownames(msbb.ad) %in% iCluster.feat$NAME, ]

################ using clsuer 4 DE genes 
## function to remove shared genes 
# comps <- c('C1', 'C2', 'C3', 'C4')
# remove_shared_genes <- function (cc) {
#   cc <- gsub('vsCO', '', cc)
#   rr <- NULL
#   for (cmp in comps[comps != cc]) {
#     res <- read.table(paste0('omics_integration/iCluster_output/DE/',cmp,'vsCO/de_res_sig_',cmp,'vsCO.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')
#     rr <- c(rr, res$GeneName)
#   }
#   return(rr)
# }

## read significant features from Kmight ADRC cohort 
C4vsALL.res <- read.table('omics_integration/iCluster_output/DE/C4vsC123/de_res_sig_C4vsC123.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')
C2vsALL.res <- read.table('omics_integration/iCluster_output/DE/C2vsC134/de_res_sig_C2vsC134.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')

C4vsALL.res.prt <- read.csv('omics_integration/data/DE_results_proteomics/prot_res_4vsall.csv', header =T, stringsAsFactors = F, check.names = F, sep=',')
C2vsALL.res.prt <- read.csv('omics_integration/data/DE_results_proteomics/prot_res_2vsall.csv', header =T, stringsAsFactors = F, check.names = F, sep=',')

#shared_genes_C4 <- remove_shared_genes('C4vsCO')
#shared_genes_C2 <- remove_shared_genes('C2vsCO')
#gg <- c(C4vsCO.res$GeneName[!C4vsCO.res$GeneName %in% shared_genes_C4],  C2vsCO.res$GeneName[!C2vsCO.res$GeneName %in% shared_genes_C2])
#gg <- c(C4vsCO.res$GeneName, C2vsCO.res$GeneName)
gg.tran <- c(C4vsALL.res$GeneName, C2vsALL.res$GeneName)
gg.prt <- c(C4vsALL.res.prt$EntrezGeneSymbol[C4vsALL.res.prt$padj < 0.05], C2vsALL.res.prt$EntrezGeneSymbol[C2vsALL.res.prt$padj < 0.05])

## extract only genes DE in C4vsALL/C2vsALL for transcriptomics 
msbb.ad <- msbb.ad[rownames(msbb.ad) %in% gg.tran, ]
## extract only genes DE in C4vsALL/C2vsALL for proteomics 
prt.data$target.gene <- gsub('\\|.*', '', rownames(prt.data))
prt.data <- prt.data[prt.data$target.gene %in% gg.prt, ]
prt.data$target.gene <- NULL
sapply(list(msbb.ad, prt.data), dim)

## select top variant features for all omics data
feature.vars = apply(msbb.ad, 1, var)
threshold = feature.vars[order(feature.vars, decreasing = T)] [floor(length(feature.vars)*0.20)]
msbb.ad <- msbb.ad[feature.vars >= threshold, ]

feature.vars = apply(prt.data, 1, var)
threshold = feature.vars[order(feature.vars, decreasing = T)] [floor(length(feature.vars)*0.20)]
prt.data <- prt.data[feature.vars >= threshold, ]
sapply(list(msbb.ad, prt.data), dim)

################################################################################
############## run/tune iClusterBayes     
################################################################################
set.seed(137) ## 137 for the current solution 
MAX.NUM.CLUSTERS = 11
num.omics <- 2
icluster.res = tune.iClusterBayes(cpus=(MAX.NUM.CLUSTERS - 1), 
                                  t(log2(msbb.ad+1)), t(prt.data),
                                  K=1:(MAX.NUM.CLUSTERS - 1),
                                  type=rep('gaussian', num.omics),  # poisson
                                  #n.burnin=12000, n.draw=18000,
                                  prior.gamma=rep(0.5, num.omics),
                                  pp.cutoff = 0.5,
                                  sdev=0.5,
                                  thin=3
)$fit
# save the result object
saveRDS(icluster.res, file="omics_integration/Replication/MSBB/BM36/2Clusters/iCluster.res.rds")
################################################################################

##########################################
## extract dev.ratio's & BIC 
#########################################
dev.ratios = lapply(1:(MAX.NUM.CLUSTERS - 1), function(i) icluster.res[[i]]$dev.ratio)
allBICs = lapply(1:(MAX.NUM.CLUSTERS - 1), function(i) icluster.res[[i]]$BIC)

## extract best cluster 
optimal.solution = icluster.res[[which.max(dev.ratios)]] 
best.clusters = optimal.solution$clusters

#dir.create('omics_integration/Replication/MSBB/BM36/2Clusters')
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

#pdf('omics_integration/Figures/k_vs_devRatio.pdf', width = 12, height = 5)
grid.arrange(p1, p2, ncol=2)
#dev.off()

k=which.max(dev.ratios)
plot(icluster.res[[k]]$beta.pp[[1]], xlab="Genes", ylab="Posterior probability", main="RNA-Seq expression")

exp.image=scale(log2(t(msbb.ad+1)))
#exp.image[exp.image > 2.5] = 2.5
#exp.image[exp.image < -2.5] = -2.5

prt.image = scale(t(prt.data))
#prt.image[prt.image > 2.5 ] = 2.5
#prt.image[prt.image< -2.5 ] = -2.5
col.scheme = alist()
col.scheme[[1]] = rev(colorRampPalette(brewer.pal(10, "RdBu"))(300))
col.scheme[[2]] = rev(colorRampPalette(brewer.pal(10, "RdBu"))(300))

pdf('omics_integration/Replication/MSBB/BM36/2Clusters/iClusterPlus_Heatmap_BM36.pdf', width=10, height = 5)
plotHMBayes(fit=optimal.solution, 
            datasets= list(exp.image, prt.image), 
            type=rep("gaussian",num.omics),
            scale = rep(F,num.omics), 
            #col.scheme = col.scheme, 
            threshold=c(0.5,0.2),
            row.order=rep(T,num.omics),  
            sparse=rep(T,num.omics),
            cap=rep(T,num.omics)
)
dev.off()

########################################################
## extract cluster membership of the best class 
########################################################
all.clusters <- matrix(data= NA, nrow= nrow(t(msbb.ad)), ncol=MAX.NUM.CLUSTERS - 1)
for (k in 1:(MAX.NUM.CLUSTERS - 1)) {
  cc = icluster.res[[k]]$clusters
  cc = as.matrix(cc)
  all.clusters[, k] = cc
}

rownames(all.clusters) = colnames(msbb.ad)
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
write.table(best.cluster.memership, file='omics_integration/Replication/MSBB/BM36/2Clusters/best.cluster.memership.tsv', sep="\t", quote = F, row.names = F)

#############################################
##### plot phenotype association 
############################################
best.cluster.memership <- read.table('omics_integration/Replication/MSBB/BM36/2Clusters/best.cluster.memership.tsv', sep="\t", header=T, check.names = F, stringsAsFactors = F)
cc.data <- data.frame(Subj_ID= co.samples, best.cluster= 0)
best.cluster.memership <- rbind(best.cluster.memership, cc.data)
bb = merge (best.cluster.memership, msbb.clin)
bb$AOD <- as.numeric(gsub('\\+', '', bb$AOD))
bb.de <- bb ## bb data frame for de 
#cmp.list = as.list(as.data.frame(combn(1: length(unique(bb$best.cluster)),2)))
write.table(bb, 'omics_integration/Replication/MSBB/BM36/2Clusters/best.cluster.memership.with.control.tsv', sep="\t", quote = F, row.names= F)
           
## compute pvalues using GLM
## for CDRe
ff <- bb[bb$best.cluster !=0, ]
ff$best.cluster[ff$best.cluster==1] <- 1
ff$best.cluster[ff$best.cluster==2] <- 0
glm.res <- glm(CDRe ~ best.cluster + Sex + AOD, data = ff, family = 'gaussian')
#coef(summary(glm.res))
pval <- data.frame(PValue= signif(coef(summary(glm.res))[,4][["best.cluster"]], digits = 3), y_pos = max(bb$CDRe, na.rm= T), best.cluster = NA)
pval$PValue <-format(pval$PValue, scientific = T)
## for AOD 
#pval <- t.test(bb$AOD[bb$best.cluster==1 & bb$AOD < 90], bb$AOD[bb$best.cluster==2 & bb$AOD < 90])$p.value
#pval <- data.frame(PValue= signif(pval, digits = 3), y_pos = max(bb$AOD), best.cluster = NA)

p <- ggplot(bb[bb$best.cluster !=0 & !is.na(bb$CDRe),], aes(x=as.factor(best.cluster), y=as.factor(CDRe), color=as.factor(best.cluster))) + 
        #geom_boxplot(aes(color = factor(best.cluster)), outlier.shape=NA, outlier.size=0.5) +
        #geom_jitter(position=position_jitter(0.3), size=1.5, aes(color = factor(best.cluster))) + theme_bw() +
        geom_point(size = 2,  shape=21, position=position_jitter(width=0.2, height=0.1)) + theme_bw() + 
        labs(x='', y='CDRe Score') + ggtitle('') + scale_color_brewer(palette = "Dark2") + 
        theme(axis.text.x=element_text(size=12, vjust=0.5, face="bold", color="black"),
                    axis.text.y=element_text(size=12, color="black", face="bold"),
                    #axis.title.x=element_text(size=19, face="bold"),
                    axis.title.y=element_text(size=14, face="bold"),
                    #panel.background=element_rect(color="black"),
                    plot.title = element_text(size = 14, hjust=0.5, color="black", face="bold"),
                    legend.position="none", panel.border = element_rect(linetype='solid', color='black'),
                    plot.margin=unit(c(1,1,1,5), 'mm')) +
         scale_x_discrete(labels=c("1"=paste0("Cluster1\n(n=",length(bb$Subj_ID[bb$best.cluster==1]),')'), "2"= paste0("Cluster2\n(n=",length(bb$Subj_ID[bb$best.cluster==2]),')'))) +
         geom_signif(data =pval, aes(xmin = 1, xmax = 2, annotations =  paste0('p=',PValue), y_position = y_pos+0.6), textsize = 3.5, step_increase=0.1, vjust = -0.2, manual = T, fontface="bold", margin_top=0.5) 


pdf('omics_integration/Replication//MSBB/BM36/2Clusters/phenotype_assoication_CDR_V2.pdf', width = 3, height = 5, useDingbats = F)
p
dev.off()

##########################################################
##### feature selection 
#########################################################
features = alist()
features[["rnaseq"]] = colnames(t(msbb.ad))
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
write.table(sig.feat, file='omics_integration/Replication/MSBB/BM36/2Clusters/sig.features.pp.rnaseq.tsv', sep="\t", quote = F, row.names = F, col.names = F)

##############################################################
##### Survival Association 
##############################################################
library(surviplot)
library(survival)

clin = bb[, c('Subj_ID', 'AOD', 'Sex', 'best.cluster')]
clin = clin[clin$AOD < 80, ]
colnames(clin) = c('sample', 'time', 'Sex', 'best.cluster')

## add clusters group 
clin$group = 'NA'
clin[clin$best.cluster ==0, 'group'] <- 'CO'
clin[clin$best.cluster ==1, 'group'] <- 'C1'
clin[!clin$best.cluster %in% c(0,1), 'group'] <- 'OT'

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
pdf(paste0('omics_integration/Replication/MSBB/', r, '/2Clusters/survival_AOD_two_groups_80.pdf'), width=5, height = 5)
surviplot(Surv(time) ~ group, data=clin[clin$group !='CO', ], ylab='Survival Proportion', xlim=c(25,max.time), 
          main ="Molecular Subtype and Overall Survival", xlab='Age of death (years)', cex.main=1,
          mark.time=TRUE, col=c("skyblue3", 'orange2'), lwd = 2.5)

dev.off()

########################################################
####### plot genes 
########################################################
# read cluster data 
best.cluster.memership <- read.table('omics_integration/Replication/MSBB/BM36/2Clusters/best.cluster.memership.tsv', sep="\t", header=T, check.names = F, stringsAsFactors = F)
cc.data <- data.frame(Subj_ID= co.samples, best.cluster= 0)
best.cluster.memership <- rbind(best.cluster.memership, cc.data)
bb = merge (best.cluster.memership, msbb.clin)
bb$AOD <- as.numeric(gsub('\\+', '', bb$AOD))

## extract expression for knight ADRC
gene <- 'NEAT1'
msbb.g.exp <- reshape2::melt(msbb[msbb$GeneName==gene, -c(1:2)], variable.name = 'Subj_ID', value.name = 'exp')
msbb.g.exp <- merge(msbb.g.exp, bb[, c('Subj_ID', 'best.cluster', 'Sex', 'AOD')], sort =F)
msbb.g.exp = msbb.g.exp[order(msbb.g.exp$best.cluster),]

## get pvalues from DE analysis 
Tran.C1vsCO <- read.table('omics_integration/Replication/MSBB/BM36/2Clusters/DE/C1vsCO/de_res_C1vsCO.tsv', header =T, stringsAsFactors = F, sep="\t")
Tran.C2vsCO <- read.table('omics_integration/Replication/MSBB/BM36/2Clusters/DE/C2vsCO/de_res_C2vsCO.tsv', header =T, stringsAsFactors = F, sep="\t")
Tran.C1vsC2 <- read.table('omics_integration/Replication/MSBB/BM36/2Clusters/DE/C1vsC2/de_res_C1vsC2.tsv', header =T, stringsAsFactors = F, sep="\t")
#Tran.ADvsCO <- read.table('omics_integration/Replication/MSBB/BM36/2Clusters/DE/iCluster_output/DE/ADvsCO/de_res_ADvsCO.tsv', header =T, stringsAsFactors = F, sep="\t")
Tran.ADvsCO <- data.frame(GeneName=gene, padj=wilcox.test(msbb.g.exp$exp[msbb.g.exp$best.cluster %in% c(1,2)], msbb.g.exp$exp[msbb.g.exp$best.cluster==0] )$p.value)

TMT.ADvsCO <- read.table('omics_integration/data/MSBB_BM36/DE_results_proteomics_TMT/01-MSBB_prot_effect_pval_ADvsCO.csv', sep=',', header =T, stringsAsFactors = F)
TMT.C1vsCO <- read.table('omics_integration/data/MSBB_BM36/DE_results_proteomics_TMT/01-MSBB_prot_effect_pval_C1vsCO.csv', sep=',', header =T, stringsAsFactors = F)
TMT.C2vsCO <- read.table('omics_integration/data/MSBB_BM36/DE_results_proteomics_TMT/01-MSBB_prot_effect_pval_C2vsCO.csv', sep=',', header =T, stringsAsFactors = F)
TMT.C1vsC2 <- read.table('omics_integration/data/MSBB_BM36/DE_results_proteomics_TMT/01-MSBB_prot_effect_pval_C1vsC2.csv', sep=',', header =T, stringsAsFactors = F)

#tmt.gene.name <- 'SNCA|P37840'
tmt.gene.name <- "GFAP|K7EPI4" #   "SNAP25|P60880" "GFAP|P14136"    
pvals <- data.frame(comp=c('ADvsCO', 'C1vsCO','C2vsCO','C1vsC2', 'TMT.ADvsCO', 'TMT.C1vsCO', 'TMT.C1vsC2', 'TMT.C2vsCO'), 
                    pval = c(Tran.ADvsCO[Tran.ADvsCO$GeneName==gene, 'padj'], Tran.C1vsCO[Tran.C1vsCO$GeneName==gene, 'padj'], 
                             Tran.C2vsCO[Tran.C2vsCO$GeneName==gene, 'padj'], Tran.C1vsC2[Tran.C1vsC2$GeneName==gene, 'padj'],  
                             TMT.ADvsCO[TMT.ADvsCO$analyte==tmt.gene.name, 'padj'] , TMT.C1vsCO[TMT.C1vsCO$analyte==tmt.gene.name, 'padj'], 
                             TMT.C1vsC2[TMT.C1vsC2$analyte==tmt.gene.name, 'padj'], TMT.C2vsCO[TMT.C2vsCO$analyte==tmt.gene.name, 'padj'] ))
pvals$pval <- format(pvals$pval, scientific = T, digits = 3)

msbb.g.exp$status.group <- 'CO'
msbb.g.exp[msbb.g.exp$best.cluster %in% c(1,2), 'status.group'] <- 'AD'
msbb.g.exp$status.group <- factor(msbb.g.exp$status.group, levels = c('CO', 'AD'))
msbb.g.exp$status.group <- factor(msbb.g.exp$status.group, levels = unique(msbb.g.exp$status.group))

# --------------------------------
## for other regions 
pvals <- NULL
for (cc in c('1vs0', '2vs0', '1vs2')) {
  g1 <- as.numeric(unlist(strsplit(cc, "vs"))[1])
  g2 <- as.numeric(unlist(strsplit(cc, "vs"))[2])
  ff <- msbb.g.exp[msbb.g.exp$best.cluster %in% c(g1, g2), ]
  glm.res <- glm(log2(exp+1) ~ best.cluster + Sex + AOD, data = ff, family = 'gaussian')
  pval <- data.frame(comp=cc, pval= signif(coef(summary(glm.res))[,4][["best.cluster"]], digits = 3), 
                     y_pos = max(log2(ff$exp+1), na.rm= T), best.cluster = NA, stringsAsFactors = F)
  pvals <- rbind(pvals, pval)
}
pvals[pvals$comp=="1vs0", 'comp'] <- 'C1vsCO'
pvals[pvals$comp=="2vs0", 'comp'] <- 'C2vsCO'
pvals[pvals$comp=="1vs2", 'comp'] <- 'C1vsC2'
pvals$pval <- format(pvals$pval, scientific = T, digits = 3)
# --------------------------------

rp1 = ggplot(msbb.g.exp, aes(x=as.factor(status.group), y=log2(exp+1))) + 
      geom_boxplot(notch = F, aes(fill=as.factor(status.group)), show.legend = F) +
      labs(x="", y="Log2(FPKM+1)") + ggtitle('') + theme_bw() + 
      #geom_jitter(aes(color=as.factor(status.group)), position=position_jitter(0.3), size=1.5) +
      theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.text.x = element_text( vjust= 1, size=22,  color="black"), 
        axis.text.y = element_text(size=28, color="black"),
        axis.title.y = element_text(size=16, face="bold"),
        panel.background = element_rect(colour = "black", size=1),legend.position = "none") + scale_fill_jco() + 
      scale_x_discrete(labels=c("CO"=paste0("Control\n(n=",table(msbb.g.exp$status.group)[["CO"]], ')'), 
                            "AD"=paste0("AD\n(n=",table(msbb.g.exp$status.group)[["AD"]],')'))) +
      geom_signif(data = pvals, aes(xmin = 1, xmax = 2, annotations =  paste0('p=',pvals$pval[pvals$comp=="ADvsCO"]), y_position = log2(max(msbb.g.exp$exp))+0.12), 
              textsize = 4.5, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") 

#scale_fill_manual(values=wes_palette("BottleRocket2"))
#scale_fill_viridis_d()
#scale_fill_brewer(palette = "Dark2")

rp44 = ggplot(msbb.g.exp, aes(x=as.factor(best.cluster), y=log2(exp+1))) + 
        geom_boxplot(aes(color=as.factor(best.cluster)), show.legend = F) +
        labs(x="", y="") + ggtitle(r) + theme_bw() + 
        #geom_jitter(aes(color=as.factor(best.cluster)), position=position_jitter(0.3), size=1.5) +
        theme(plot.title = element_text(hjust=0.5, size=20, face="bold"),
            axis.text.x = element_text( vjust= 1, size=16, face="bold", color="black"), 
            axis.text.y = element_text(size=16, face="bold"),
            axis.title.y = element_text(size=18, face="bold"),
            legend.position = "none")  +
        #scale_color_manual(values=c("0"="gray10", "1"="#E7298A", "2"="#00AFBB")) + #
        
        scale_x_discrete(labels=c("0"=paste0("Control\n(n=",length(msbb.g.exp$Subj_ID[msbb.g.exp$best.cluster==0]),')'), 
                            "1"=paste0("Cluster1\n(n=",length(msbb.g.exp$Subj_ID[msbb.g.exp$best.cluster==1]),')'), 
                            "2"= paste0("Cluster2\n(n=",length(msbb.g.exp$Subj_ID[msbb.g.exp$best.cluster==2]),')'))) +
        geom_signif(data = pvals, aes(xmin = 1, xmax = 1.97, annotations =  paste0('p=',pval[comp=="C1vsCO"]), y_position = log2(max(msbb.g.exp$exp))+0.15), 
              textsize = 4.5, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") + 
        geom_signif(data = pvals, aes(xmin = 1, xmax = 3, annotations =  paste0('p=',pval[comp=="C2vsCO"]), y_position = log2(max(msbb.g.exp$exp))+0.43), 
              textsize = 4.5, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") +
        geom_signif(data = pvals, aes(xmin = 2.03, xmax = 3, annotations =  paste0('p=',pval[comp=="C1vsC2"]), y_position = log2(max(msbb.g.exp$exp))+0.15), 
              textsize = 4.5, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") 
      

#prt.g.exp <- as.data.frame(t(prt.exp[rownames(prt.exp)==gene, ]))
#gene2 = "SNAP25"  # GFAP|P14136 GFAP|P14136-3  GFAP|K7EKD1 GFAP|K7EPI4
prt.g.exp <- as.data.frame(t(prt.exp[rownames(prt.exp)==tmt.gene.name, ]))
#prt.g.exp <- as.data.frame(t(prt.exp[grepl(gene, rownames(prt.exp)), ]))
prt.g.exp$Subj_ID <- rownames(prt.g.exp)
rownames(prt.g.exp) <- NULL
colnames(prt.g.exp) <- c('prt.exp', 'Subj_ID')
prt.g.exp <- merge(prt.g.exp, bb[, c('Subj_ID', 'best.cluster')], sort =F)
prt.g.exp = prt.g.exp[order(prt.g.exp$best.cluster),]

prt.g.exp$status.group <- 'CO'
prt.g.exp[prt.g.exp$best.cluster %in% c(1,2), 'status.group'] <- 'AD'
prt.g.exp$status.group <- factor(prt.g.exp$status.group, levels = c('CO', 'AD'))
prt.g.exp$status.group <- factor(prt.g.exp$status.group, levels = unique(prt.g.exp$status.group))

g.prt1 = ggplot(prt.g.exp, aes(x=as.factor(status.group), y=prt.exp)) + geom_boxplot(aes(color=as.factor(status.group)), show.legend = F) +
        labs(x="", y="Protein Expression") + ggtitle('') + theme_bw() + 
        geom_jitter(aes(color=as.factor(status.group)), position=position_jitter(0.3), size=1.5) +
        theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
              axis.text.x = element_text( vjust= 1, size=14, face="bold", color="black"), 
              axis.text.y = element_text(size=14, face="bold"),
              axis.title.y = element_text(size=16, face="bold"),
              legend.position = "none") + scale_color_jco() + 
        scale_x_discrete(labels=c("CO"=paste0("Control\n(n=",table(prt.g.exp$status.group)[["CO"]], ')'), 
                            "AD"=paste0("AD\n(n=",table(prt.g.exp$status.group)[["AD"]],')'))) +
        geom_signif(data = pvals, aes(xmin = 1, xmax = 2, annotations =  paste0('p=',pvals$pval[pvals$comp=="TMT.ADvsCO"]), y_position = max(prt.g.exp$prt.exp)+0.12), 
              textsize = 4.5, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") 

g.prt2 <- ggplot(prt.g.exp, aes(x=as.factor(best.cluster), y=prt.exp)) + geom_boxplot(aes(color=as.factor(best.cluster)), show.legend = F) +
           labs(x="", y="") + ggtitle('') + theme_bw() + 
           geom_jitter(aes(color=as.factor(best.cluster)), position=position_jitter(0.3), size=1.5) +
           theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
                axis.text.x = element_text( vjust= 1, size=14, face="bold", color="black"), 
                axis.text.y = element_text(size=14, face="bold"),
                axis.title.y = element_text(size=16, face="bold"), legend.position = "none") + 
          scale_color_manual(values=c("0"="gray10", "1"="#E7298A", "2"="#00AFBB")) + #
          scale_x_discrete(labels=c("0"=paste0("Control\n(n=",length(prt.g.exp$Subj_ID[prt.g.exp$best.cluster==0]),')'), 
                            "1"=paste0("Cluster1\n(n=",length(prt.g.exp$Subj_ID[prt.g.exp$best.cluster==1]),')'), 
                            "2"= paste0("Cluster2\n(n=",length(prt.g.exp$Subj_ID[prt.g.exp$best.cluster==2]),')'))) +
          geom_signif(data = pvals, aes(xmin = 1, xmax = 1.97, annotations =  paste0('p=',pval[comp=="TMT.C1vsCO"]), y_position = max(prt.g.exp$prt.exp)+0.08), 
              textsize = 4.5, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") + 
          geom_signif(data = pvals, aes(xmin = 1, xmax = 3, annotations =  paste0('p=',pval[comp=="TMT.C2vsCO"]), y_position = max(prt.g.exp$prt.exp)+0.34), 
              textsize = 4.5, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") +
          geom_signif(data = pvals, aes(xmin = 2.03, xmax = 3, annotations =  paste0('p=',pval[comp=="TMT.C1vsC2"]), y_position = max(prt.g.exp$prt.exp)+0.08), 
              textsize = 4.5, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") 

pdf(paste0('omics_integration/Replication/MSBB/', r,'/2Clusters/', gene, '_exp_BM10_22_44.pdf'), width=12, height = 7)
#grid.arrange(arrangeGrob(grobs=list(rp1, rp2, g.prt1, g.prt2), ncol=2, widths=c(1.5,3)))
grid.arrange(arrangeGrob(grobs=list(rp10,rp22,rp44), ncol=3, widths=c(4,4,4)))
dev.off()

# dd <- data.frame(samples = bb$sampleIdentifier[bb$best.cluster==2], age=bb$AOD[clin$best.cluster==2])
# dd$age <- as.numeric(gsub('\\+', '', dd$age))
# ggplot(dd, aes(x=samples, y=age)) + geom_point(size=5) + geom_hline(yintercept = 80, color="red") +
#   labs(title="Distribution of AOD for C1", xlab="samples", y='Age of Death')
# 
# t.test(clin.data$AOD[clin.data$Subj_ID %in% gsub('ID_','', c4)], clin$time[clin$best.cluster==1])

################################################################
#### Check cell proportions 
###############################################################
## plot cell-type deconvolution
algs = c( "ssFrobenius", "meanProfile")
cellProp.res <- NULL
for (alg in algs) {
  s.res <- read.table(paste0('/home/general/Public_Data/bulkRNASeq/201812_MSBB/BM36/05.-Analyses/deconvolution_analysis/deconvolution_', alg,'_results.tsv'), header = T, check.names = F, stringsAsFactors = F)
  cellProp.res <- rbind(cellProp.res, s.res)
}

## keep representative samples 
cellProp.res <- cellProp.res[cellProp.res$Subject %in% subjects_to_keep, ]

## merge with clinical data 
colnames(cellProp.res)[colnames(cellProp.res)=="Subject"] <- "sampleIdentifier"
cellProp.res <- merge(cellProp.res, msbb.tech[,c('sampleIdentifier','individualID')])

# merge with the best clusters 
colnames(cellProp.res)[colnames(cellProp.res)=="individualID"] <- "Subj_ID"
cellProp.res <- merge(cellProp.res, best.cluster.memership)
cellProp.res$sampleIdentifier <- NULL
cell.prop <- cellProp.res ## one copy for the DE analysis

## write for supplementary table 
res_for_supp <- cellProp.res
res_for_supp <- res_for_supp[res_for_supp$Algorithm=="ssFrobenius", -2]
res_for_supp <- res_for_supp[, c('Subj_ID','best.cluster','Astrocyte', 'Microglia', 'Neuron', 'Oligodendrocyte')]
write.table(res_for_supp, file='omics_integration/Replication/MSBB/BM36/2Clusters/Cell_proportions_MSBB_BM36.tsv', sep="\t", quote = F, row.names = F)
## ------------------

cellProp.res <- cellProp.res[cellProp.res$best.cluster !=0, ]
cellProp.res <- reshape2::melt(cellProp.res, id.vars = c("Subj_ID", "best.cluster", "Algorithm"))
colnames(cellProp.res) <- c('Subj_ID', 'best.cluster', 'Algorithm', 'Cell_type', 'Proportion')

## compute pvalues using GLM
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
      theme(axis.text.x=element_text(size=18, vjust=1, color="black", face="bold", angle=45, hjust=1),
            axis.text.y=element_text(size=18, color="black", face="bold"), 
            axis.title.y=element_text(size=18, color="black", face="bold"),
            plot.title = element_text(size = 22, hjust=0.5, color="black", face="bold"),
            legend.position="none", panel.border = element_rect(linetype='solid', color='black')) +
      #scale_x_discrete(labels=c("1"=paste0("Cluster1\n(n=",length(unique(best.cluster.memership$Subj_ID[best.cluster.memership$best.cluster==1])),')'), 
      #                          "2"= paste0("Cluster2\n(n=",length(unique(best.cluster.memership$Subj_ID[best.cluster.memership$best.cluster==2])),')'))) +
      scale_x_discrete(labels=c("1"="Cluster1", "2"= "Cluster2")) + 
      geom_signif(data = dd, aes(xmin = 1, xmax = 2, annotations =  paste0('p=', Pvalue), y_position = y_pos), textsize = 4, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold", margin_top=0.2) 
    
    all.plots[[cell]] <- p
    
  }
  
  ## print 
  myleft <- textGrob('Cell Proportion', gp=gpar(fontsize=22, fontface="bold"), rot=90)
  mybottom <- textGrob(paste0("(Cluster1=", length(bb$Subj_ID[bb$best.cluster==1]), ", Cluster2=", length(bb$Subj_ID[bb$best.cluster==2]),')'), gp=gpar(fontsize=22, fontface="bold"))
  pdf(paste0("omics_integration/Replication/MSBB/BM36/2Clusters/deconvolution_celllines_per_cluster_",alg,".pdf"), width = 12, height = 6, useDingbats = F)
  grid.arrange(grobs=all.plots, nrow=1, ncol=4, left = myleft, bottom=mybottom)
  dev.off()
}


## plot per cell type 
#cmp.list = as.list(as.data.frame(combn(1: length(unique(cellProp.res$best.cluster)),2)))
#pdf("omics_integration/Replication/MSBB/BM36/2Clusters/deconvolution_celllines_per_cluster_all.pdf", width = 12, height = 8)
# for (alg in c( "ssFrobenius", "meanProfile")) {
#   dd <- all.pvals[all.pvals$Algorithm==alg, ]
#   p <- ggplot(cellProp.res[cellProp.res$Algorithm==alg, ], aes(x=factor(best.cluster), y=Proportion, group=factor(best.cluster))) + 
#       geom_boxplot(outlier.shape=NA, outlier.size=0.5, aes(color=factor(best.cluster))) +
#       geom_jitter(position=position_jitter(0.3),  size=1.8, aes(color=factor(best.cluster))) + theme_bw() + 
#       facet_wrap(. ~ Cell_type, scales = "free") + scale_color_brewer(palette = "Dark2") +
#       labs(x='', y='Cell Proportion') + ggtitle(paste0('Overall Cell Proportion Per Cluster - (Algorithm: ', alg,')')) + 
#       theme(axis.text.x=element_text(size=11, vjust=0.5, color="black"),
#           strip.text = element_text(size = 14, face="bold"), 
#           axis.text.y=element_text(size=11, color="black"), 
#           axis.title.y=element_text(size=14, face="bold"),
#           plot.title = element_text(size = 16, hjust=0.5, color="black", face="bold"),
#           legend.position="none", panel.border = element_rect(linetype='solid', color='black')) + 
#      #geom_signif(comparisons= cmp.list, step_increase=0.1, textsize = 3.5, map_signif_level=function(p)sprintf("p = %.2g", p), show.legend=F)
#      #geom_text(data = dd, mapping = aes(x = -Inf, y = -Inf , label = paste0('p=',Pvalue)), hjust=-0.50, vjust=-0.9)
#      geom_signif(data = dd, aes(xmin = 1, xmax = 2, annotations =  paste0('p=', Pvalue), y_position = y_pos), textsize = 3, step_increase=0.1, vjust = -0.2,manual = T) 
# 
#   print(p)
#   
# }
# dev.off()


################################################################
##### DE Analysis  
################################################################
dir.create(paste0('omics_integration/Replication//MSBB/', r, '/2Clusters/DE'))

## read raw counts 
#msbb.cts <- read.table(paste0('/home/general/Public_Data/bulkRNASeq/201812_MSBB/',r, '/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies/gene_quant_matrix_NumReads.tsv'), header =T, check.names = F, stringsAsFactors = F)
msbb.cts <- read.table(paste0('/home/general/Public_Data/bulkRNASeq/201812_MSBB/', r, '/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_count_matrix.tsv'), header =T, check.names = F, stringsAsFactors = F)
colnames(msbb.cts)[colnames(msbb.cts)=="gene"] <- "GeneID"
msbb.cts <- merge(msbb.meta.cols, msbb.cts)

## replace msbb count header 
msbb.cts <- msbb.cts[, c('GeneID', 'GeneName', 'GeneBiotype', subjects_to_keep)]
for (k in 1:length(subjects_to_keep)) {
  s.id <- subjects_to_keep[k]
  ind.id <- msbb.tech$individualID[msbb.tech$sampleIdentifier==s.id]
  colnames(msbb.cts)[colnames(msbb.cts)==s.id] <- ind.id
}
rownames(msbb.cts) <- paste0(msbb.cts$GeneID,'_', msbb.cts$GeneName)
msbb.cts$GeneID <- NULL
msbb.cts$GeneName <- NULL
msbb.cts$GeneBiotype <- NULL

## merge with cell proportions 
# get the max percentage for duplicated IDs
cell.prop.res <- cell.prop
cell.prop.res <- cell.prop.res[cell.prop.res$Algorithm=='ssFrobenius', c('Subj_ID', 'Astrocyte', 'Neuron')]
xx = aggregate(cell.prop.res$Astrocyte, by = list(cell.prop.res$Subj_ID), max)
colnames(xx) <- c('Subj_ID', 'Astrocyte')
cell.prop.res <- merge(cell.prop.res, xx, by=c('Subj_ID', 'Astrocyte'))
bb <- merge(bb, cell.prop.res)


## organize data 
c1.samples <- bb$Subj_ID[bb$best.cluster==1]
c2.samples <- bb$Subj_ID[bb$best.cluster==2]
msbb.cts <- msbb.cts[, c(co.samples, c1.samples, c2.samples)]

## prepare data for ADvsCO
ad.samples <- msbb.clin$Subj_ID[msbb.clin$Status=="AD"]
co.samples <- msbb.clin$Subj_ID[msbb.clin$Status=="control"]
ad.samples <- ad.samples[ad.samples %in% colnames(msbb.cts)]
co.samples <- co.samples[co.samples %in% colnames(msbb.cts)]
msbb.cts <- msbb.cts[, c(co.samples, ad.samples)]
bb <- msbb.clin[grepl("AMPAD_", bb$Subj_ID), ]
rownames(bb) <- bb$Subj_ID
bb <- bb[c(co.samples, ad.samples), ] 
all(rownames(bb) == colnames(msbb.cts) )

bb <- bb[c(co.samples, c1.samples, c2.samples), ]
rownames(bb) <- bb$Subj_ID
## check samples consistency 
all(rownames(bb) == colnames(msbb.cts) )

## function to filter lowly expressed genes 
msbb.cts_norm <- cpm(msbb.cts)
selectGenes <- function(grp1, grp2) {
  grp1.t <- round(length(grp1) * 0.25)
  grp2.t <- round(length(grp2) * 0.25)
  keep <- (rowSums (msbb.cts_norm[,grp1]> 0.5) >= grp1.t) | (rowSums (msbb.cts_norm[,grp2]> 0.5) >= grp2.t)
  if (all(rownames(msbb.cts) != rownames(msbb.cts_norm))) {
    stop('rownames of the normalized counts are not the same as the raw counts ')
  } else {
    cmp.counts <- msbb.cts[keep, ]
  }
  return (cmp.counts)
}

## Function for fitting the DESeq2 model
#ph.data$AOD <- as.numeric(gsub("\\+", "", ph.data$AOD))
DESeqModelFit <- function(count_data, model.group) {
  dds <- DESeqDataSetFromMatrix(countData = count_data, 
                                colData = bb, 
                                design = formula(paste("~ Sex + AOD + Astrocyte + Neuron +", model.group)))
  ## for ADvsCO
  dds <- DESeqDataSetFromMatrix(countData = count_data, 
                                colData = bb, 
                                design = formula(~ Sex + AOD + Status))
  ## Test for DE 
  dds <- DESeq(dds)
  print(resultsNames(dds))
  return (dds)
}

## make groups 
bb$best.cluster[bb$best.cluster==1] <- "C1"
bb$best.cluster[bb$best.cluster==2] <- "C2"
bb$best.cluster[bb$best.cluster==0] <- "CO"

comps <- c('C1vsCO', 'C2vsCO', 'C1vsC2')
for (cmp in comps) {
  cat('Running for:', cmp, '\n')
  dir.create(paste0('omics_integration/Replication//MSBB/',r,'/2Clusters/DE/', cmp))
  
  ## extract groups 
  mygroup <- "best.cluster"
  grp1.samples <- bb$Subj_ID[bb$best.cluster== unlist(strsplit(cmp, "vs"))[1] ]
  grp2.samples <- bb$Subj_ID[bb$best.cluster== unlist(strsplit(cmp, "vs"))[2] ]
 
  ## remove lowly expressed genes 
  grp.counts <- selectGenes(grp1.samples, grp2.samples)
  grp.counts <- selectGenes(co.samples, ad.samples) ## for AD vs CO
  
  ## fit DESeq model
  dds <- DESeqModelFit(grp.counts, mygroup)
  #resultsNames(dds)
  
  ## extract result table
  dds.res <- results(dds, alpha=0.05, contrast = c(mygroup, unlist(strsplit(cmp, "vs"))[1], unlist(strsplit(cmp, "vs"))[2] ), tidy = F)
  summary(dds.res)
  
  ## add annotation 
  gg <- as.data.frame(str_split_fixed(rownames(dds.res), "_", 2))
  colnames(gg) <- c('GeneID', "GeneName")
  dds.res <- cbind(dds.res, gg)
  dds.res <- as.data.frame(merge(msbb.meta.cols , dds.res))
  
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
  colnames(dds.res)[colnames(dds.res) =="grp1"] <- paste0('grp.', unlist(strsplit(cmp, "vs"))[1], '.ExpMean')
  colnames(dds.res)[colnames(dds.res) =="grp2"] <- paste0('grp.', unlist(strsplit(cmp, "vs"))[2], '.ExpMean')
  
  ## plot volcano plot
  top.10.genes <- dds.res$GeneName[1:10]
  pdf(paste0('omics_integration/Replication/MSBB/', r,'/2Clusters/DE/',cmp, '/volcano_plot_',cmp,'2.pdf'), width = 6, height = 8)
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

  write.table(dds.res, file=paste0('omics_integration/Replication/MSBB/',r, '/2Clusters/DE/',cmp,'/de_res_', cmp,'.tsv'), sep="\t", quote = F, row.names = F)
  write.table(dds.res.sig, file=paste0('omics_integration/Replication/MSBB/', r,'/2Clusters/DE/',cmp, '/de_res_sig_',cmp,'.tsv'), sep="\t", quote = F, row.names = F)
  
}

#################################################################
##### run pathways analysis 
#################################################################
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
dbs <- listEnrichrDbs()
comps <- c('C1vsCO', 'C2vsCO', 'C1vsC2')
rm(en.up); rm(en.dn)
for (cc in comps) {
  
    de.res <- read.table(paste0('omics_integration/Replication/MSBB/', r,'/2Clusters/DE/',cc, '/de_res_sig_',cc,'.tsv'), header =T, stringsAsFactors = F)
    
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
    write.table(go_bp_up, file=paste0('omics_integration/Replication/MSBB/', r, '/2Clusters/DE/', cc,'/',cc,'_GO_BP_UP.tsv'), sep="\t", quote = F, row.names = F)
    write.table(go_bp_dn, file=paste0('omics_integration/Replication/MSBB/', r, '/2Clusters/DE/', cc,'/',cc,'_GO_BP_DN.tsv'), sep="\t", quote = F, row.names = F)
}


#####################################################
#### check overlap with Knight ADRC
#####################################################
# transcriptomics
Discovery <- read.table('omics_integration/iCluster_output/DE/C4vsC123/de_res_C4vsC123.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')
Replicate <- read.table('omics_integration/Replication/MSBB/BM36/2Clusters/DE/C1vsC2/de_res_C1vsC2.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')
all.genes <- union(read.table('omics_integration/iCluster_output/DE/C4vsC123/de_res_C4vsC123.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')[,2], 
                   read.table('omics_integration/Replication/MSBB/BM36/2Clusters/DE/C1vsC2/de_res_C1vsC2.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')[,2])
# proteomics 
Discovery <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_4vsall.csv', header =T, stringsAsFactors = F, check.names = F, sep=',')
Replicate <- read.csv('omics_integration/data/MSBB_BM36/DE_results_proteomics_TMT/01-MSBB_prot_effect_pval_C1vsC2.csv', header =T, stringsAsFactors = F, check.names = F, sep=',')
disc.gg <- unlist(strsplit(Discovery$EntrezGeneSymbol, " "))
disc.gg <- unlist(strsplit(disc.gg, ","))
all.genes <- union(unique(disc.gg), unique(gsub("\\|.*","", Replicate$analyte)))


for (d in c('up', 'dn')) {
  
  ## transcriptomics
  Discovery <- Discovery$GeneName[Discovery$direction==d & Discovery$padj < 0.05]
  Replicate <- Replicate$GeneName[Replicate$direction==d & Replicate$padj < 0.05]
  #write.table(intersect(Discovery, Replicate), file = paste0('omics_integration/Replication/MSBB/BM36/2Clusters/DE/knightADRC_BM36_gene_overlap_',d,'_pval.tsv'), quote = F, sep="\t", row.names = F)

  #proteomics 
  if (d=="up") { 
    Discovery <- Discovery$EntrezGeneSymbol[Discovery$effect > 0 & Discovery$padj < 0.05] 
    Replicate <- Replicate$analyte[Replicate$effect > 0 & Replicate$padj < 0.05]
  } else {
    Discovery <- Discovery$EntrezGeneSymbol[Discovery$effect < 0 & Discovery$padj < 0.05] 
    Replicate <- Replicate$analyte[Replicate$effect < 0 & Replicate$padj < 0.05]
  }
  Discovery <- unlist(strsplit(Discovery, " "))
  Discovery <- unique(unlist(strsplit(Discovery, ",")))
  Replicate <- unique(gsub("\\|.*","", Replicate))
  
  ## run hypergeometirc test 
  #all.genes <- nrow(msbb.cts)
  test.mat <- matrix(c( length(all.genes) - length(union(Discovery, Replicate)), length(sediff(Discovery, Replicate)), length(setdiff(Replicate, Discovery)), length(intersect(Discovery, Replicate))), nrow=2)
  fisher.pval <- fisher.test(test.mat, alternative = "two.sided")$p.value
  if (fisher.pval ==0) { fisher.pval ='Fisher\'s pval < 2.2e-16'} else {fisher.pval = paste0('Fisher\'s pval = ',signif(fisher.pval, digits = 4)) }
  
  ## generate venn diagram 
  fileName <- 'Disc_C4vsOT_vs_Rep_C1vsOT'
  myCatNames <- c(paste0("Knight ADRC\n(",length(Discovery),")") , paste0("MSBB (BM36)\n(",length(Replicate),")"))
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
    cat.pos = c(-0.1, -15),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans", cat.fontface = "bold", 
    #cat.col = c('#fde725ff', "#440154ff"),  
    cat.col = pal_jco()(2)
    #rotation = 1
  )
  
  pdf(paste0('omics_integration/Replication/MSBB/', r, '/2Clusters/DE/', fileName, '_venn_', d,'.pdf'), width=2.7, height=2.7)
  grid.draw(vv)
  dev.off()
  
}


#################################################################
##### check homeostatic and activated mic genes 
#################################################################
## function to make stat table
make_stat_table <- function (g) {
  for (cmp in  c('C1vsCO', 'C2vsCO', 'C1vsC2')) {
    res <- read.table(paste0('omics_integration/Replication/MSBB/BM36/2Clusters/DE/',cmp,'/de_res_',cmp,'.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')
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

## process prt.exp data
prt.exp$target.gene <- gsub('\\|.*', '', rownames(prt.exp))

## plot genes 
pdf('omics_integration/Replication/MSBB/BM36/2Clusters/mic.activated.genes.pdf', width=12, height = 6)
for (gg in activated) {
  
  gg.fpkm <- reshape2::melt(msbb.fpkm[msbb.fpkm$GeneName==gg, !colnames(msbb.fpkm) %in% c('GeneID', 'GeneBiotype')], id.vars = "GeneName")
  gg.tpm <- reshape2::melt(msbb.tpm[msbb.tpm$GeneName==gg, !colnames(msbb.tpm) %in% c('GeneID', 'GeneBiotype')], id.vars = "GeneName")
  colnames(gg.fpkm) <- c('GeneName', 'Sample_Name', 'fpkm')
  colnames(gg.tpm) <- c('GeneName', 'Sample_Name', 'tpm')
  
  ## merge with cell lines information  
  gg.exp <- merge(gg.fpkm, gg.tpm)
  
  ## merge with clinical data 
  gg.exp <- gg.exp[gg.exp$Sample_Name %in% subjects_to_keep, ]
  colnames(gg.exp)[colnames(gg.exp)=="Sample_Name"] <- "sampleIdentifier"
  gg.exp <- merge(gg.exp, msbb.tech[,c('sampleIdentifier','individualID')])
  # merge with the best clusters 
  colnames(gg.exp)[colnames(gg.exp)=="individualID"] <- "Subj_ID"
  gg.exp <- merge(gg.exp, best.cluster.memership)
  
  ## get proteomics data
  g.prt = reshape2::melt(prt.exp[prt.exp$target.gene == gg, colnames(prt.exp) !="target.gene"])
  if (nrow(g.prt) > 0) {
    colnames(g.prt) = c('Subj_ID', 'prt.exp')
    gg.exp = merge(gg.exp, g.prt, sort =F) 
  }
  
  ## make the DE info table
  mytable1 <- make_stat_table(gg)
  if (!is.null(mytable1)) {
    #mytable1 <- ggtexttable(mytable1, rows = NULL, theme = ttheme("mBlue", base_size = 6, padding = unit(c(2, 2),"mm")))  
    mytitle <- paste0('\nLogFC (C1vsCO, C2vsCO, C1vsC2): ', paste(mytable1[1, -1], collapse = '; '), 
                      '\nPval (C1vsCO, C2vsCO, C1vsC2): ',  paste(mytable1[2, -1], collapse = ';, '),
                      '\nFDR (C1vsCO, C2vsCO, C1vsC2): ',   paste(mytable1[3, -1], collapse = '; ')
    )
  }
  
  g1 = ggplot(gg.exp, aes(x=as.factor(best.cluster), y=log2(fpkm+1))) + geom_boxplot(aes(fill=as.factor(best.cluster)), show.legend=F) +  
        geom_jitter(position=position_jitter(0.3), size=1) +
        labs(x="", y="Log2(FPKM+1)") + ggtitle(paste0('RNA-Seq expression profiles for ', gg, mytitle)) + theme_bw() +
        theme(plot.title = element_text(hjust=0.5, size=13, face="bold"),
            axis.text.x = element_text(vjust= 0.5, size=12, face="bold"), 
            axis.text.y = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"), 
            legend.position = "none") + scale_fill_brewer(palette = 'Dark2' ) 
  # geom_signif(comparisons= list(c('CO', 'C1'), c('CO', 'C2'), c('CO', 'C3'), c('CO', 'C4')),  
  #         step_increase=0.1, textsize = 4, fontface="bold", map_signif_level=function(p)sprintf("p = %.2g", p), show.legend=F) 
  #if (!is.null(mytable1)) { g4 = g4 + annotation_custom(ggplotGrob(mytable1), xmin= 0.5, xmax=8.5, ymin=max(log2(g.exp$rnaseq.exp+1))-0.5, ymax=max(log2(g.exp$rnaseq.exp+1))+0.5) }
  
  if (nrow(g.prt) > 0) {
    g2 = ggplot(gg.exp, aes(x=as.factor(best.cluster), y=prt.exp)) + geom_boxplot(aes(fill=as.factor(best.cluster)), show.legend=F) +  
        geom_jitter(position=position_jitter(0.3), size=1) +
        labs(x="", y="protein.exp") + ggtitle(paste0('Proteomic expression profiles for ', gg)) + theme_bw() +
        theme(plot.title = element_text(hjust=0.5, size=13, face="bold"),
            axis.text.x = element_text(vjust= 0.5, size=12, face="bold"), 
            axis.text.y = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"), 
            legend.position = "none") + scale_fill_brewer(palette = 'Dark2' ) 
        #geom_signif(comparisons= list(c('CO', 'C1'), c('CO', 'C2'), c('CO', 'C3'), c('CO', 'C4')),  
        #          step_increase=0.1, textsize = 4, fontface="bold", map_signif_level=function(p)sprintf("p = %.2g", p), show.legend=F)
  } else {
    g2 = ggplot()+geom_point() + annotate("text", x = 0, y = 0, label = "Gene was not found in proteomics", size=6) + labs(x='', y='')
  }
  
  #grid.arrange(g4, g5, g1, g2, g3, ncol=2, nrow =3, heights = c(4,4,5))
  grid.arrange(g1, g2, ncol=2, nrow =1) #
  
}
dev.off()

#################################################################################
## prepare data for GSEA 
dir.create('omics_integration/data/4GSEA/BM36/')
res <- read.table('omics_integration/Replication/MSBB/BM36/2Clusters/DE/C1vsC2/de_res_sig_C1vsC2.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')
res <- msbb[msbb$GeneName %in% res$GeneName, c('GeneName', bb$Subj_ID[bb$best.cluster==1], bb$Subj_ID[bb$best.cluster==2])]
res$DESCRIPTION = NA
res <- res[, c('GeneName', 'DESCRIPTION', bb$Subj_ID[bb$best.cluster==1], bb$Subj_ID[bb$best.cluster==2])]
colnames(res)[colnames(res)=="GeneName"] <- 'NAME'
write.table(res, file="omics_integration/data/4GSEA/BM36//exp_data_C1vsC2.gct", sep="\t", quote = F, row.names = F)

## make the phenotype file 
ph <- c(rep('C1', length(bb$Subj_ID[bb$best.cluster==1]) ), rep('C2', length(bb$Subj_ID[bb$best.cluster==2])) )
ph <- paste(ph, collapse = " ")
write.table(ph, file="omics_integration/data/4GSEA/BM36/phenotype_class_C1vsC2.cls", sep=" ", quote = F, row.names = F)
#################################################################################

#################################################################################################
#### Single-cell integration 
#################################################################################################
setwd('/home/eteleeb/projects')
dir.create('omics_integration/Replication/MSBB/BM36/2Clusters/cell_specific_pathways')

folder_names <- list.dirs('/home/delaguilaj/DE_soma/', full.names = F, recursive = F)

## function to remove shared genes 
comps <- c('C1', 'C2')
remove_shared_genes <- function (cc) {
  cc <- gsub('vsCO', '', cc)
  rr <- NULL
  for (cmp in comps[comps != cc]) {
    res <- read.table(paste0('omics_integration/Replication/MSBB/BM36/2Clusters/DE/',cmp,'vsCO/de_res_sig_',cmp,'vsCO.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')
    rr <- c(rr, res$GeneName)
  }
  return(rr)
}


var <- 'estimate'
ll <- 'Estimate'
tt <- 'estimate'

## loop through the whole data 
for (d in c('all', 'up', 'dn')) {
    ff.plots <- list()
    tmp = NULL
    for (c in c('C1vsCO', 'C2vsCO')) {
      cat('   Running for: ', c, '\n')
      rr1 <- read.table(paste0('omics_integration/Replication/MSBB/BM36/2Clusters/DE/', c, '/de_res_sig_', c,'.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')
      #dir.create(paste0('omics_integration/iCluster_output/cell_specific_pathways/', c))
      
      #rr1 <- rr1[abs(rr1$log2FoldChange) > 1, ]
      
      ## remove shared genes 
      shared_genes <- remove_shared_genes(c)
      rr1 <- rr1[!rr1$GeneName %in% shared_genes, ]
      
      ## combine all sets 
      folder.data <- NULL
      for (clu in c(0,1,2,3,4,5)) {
        cluster.data <- NULL
        for (ff in folder_names) {
          sn.data <- read.table(list.files(paste0('/home/delaguilaj/DE_soma/', ff) , pattern=paste0('glmmTMB_celltype_', clu), full.names = T), sep="\t", header =F, stringsAsFactors = F)
          colnames(sn.data) <- c('GeneName', 'estimate', 'Std.Error', 'ZVvalue', 'PValue')
          sn.data <- sn.data[sn.data$PValue < 0.05, ]
          sn.data$cluster <- clu
          cluster.data <- rbind(cluster.data, sn.data)
        }
          ## extract top 10 
          upper_cutoff <- quantile(cluster.data$estimate, prob=0.90)
          cluster.data <- cluster.data[cluster.data$estimate > upper_cutoff, ]
       
          folder.data <- rbind(folder.data, cluster.data)
      }
      
      ## merge with omics data 
      if (d == 'all') {
        folder.data <- merge(folder.data, rr1)
      } else if (d=='up') {
        folder.data <- merge(folder.data, rr1[rr1$direction=="up",])
      } else if (d=='dn') {
        folder.data <- merge(folder.data, rr1[rr1$direction=="dn",])
      }
      
      num.genes <- length(unique(folder.data$GeneName))
      ## plot the results 
      folder.data[folder.data$cluster==0, 'cluster'] <- 'Oligo'
      folder.data[folder.data$cluster==1, 'cluster'] <- 'Neuron'
      folder.data[folder.data$cluster==2, 'cluster'] <- 'Microglia'
      folder.data[folder.data$cluster==3, 'cluster'] <- 'Astrocytes'
      folder.data[folder.data$cluster==4, 'cluster'] <- 'Opc'
      folder.data[folder.data$cluster==5, 'cluster'] <- 'Endo'
      folder.data$cluster <- factor(folder.data$cluster, levels =c("Astrocytes", "Microglia", "Neuron", "Oligo", "Opc", "Endo"))
      
      ## plot results 
      p <-  ggplot(folder.data, aes(x=factor(cluster), y= -log10(get(var)), group = cluster)) + theme_bw(base_size = 12)
      #p <-  ggplot(folder.data, aes(x=factor(cluster), y= get(var), group = cluster)) + theme_bw(base_size = 12) +
      p <- p + geom_boxplot(lwd =0.5, aes(color=factor(cluster)), show.legend = F) + 
        geom_jitter(size=1.5, alpha= 0.5, shape=1, width=0.3, aes(color=factor(cluster)), show.legend = F) +  
        labs(x='', y = '', title = paste0(c, ' (n=', num.genes,')')) + 
        theme(axis.text.x=element_text(size=14, angle =90, vjust=0.5, face="bold"), 
              axis.text.y=element_text(size=14, face="bold"), 
              axis.title.y=element_text(size=16, face="bold"),
              plot.title = element_text(size = 18, hjust=0.5, face="bold"),
              axis.ticks=element_blank(), panel.grid.minor = element_blank(),
              legend.title=element_blank(), legend.background = element_rect(color = "black"),
              panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 
      
      ff.plots[[c]] <- p
      folder.data$comp <- c
      tmp = rbind(tmp, folder.data[, c('GeneName', 'estimate', 'Std.Error', 'ZVvalue', 'PValue', 'cluster', 'comp')])
    } ## end comparisons 
    
    ## write results 
    write.table(tmp, paste0('omics_integration/Replication/MSBB/BM36/2Clusters/cell_specific_pathways/', d,'_top10_cluster_specific.tsv'), sep="\t", row.names = F, quote = F)
    
    ## plot all 
    mytitle = textGrob(paste0('Distribution of the ', tt,' for the top 10 DE cell type - (', d, ')\n'), gp=gpar(fontsize=20, fontface="bold"))
    myleft = textGrob(ll, gp=gpar(fontsize=20, fontface="bold"), rot = 90)
    pdf(paste0('omics_integration/Replication/MSBB/BM36/2Clusters/cell_specific_pathways/', d, '_top10_cluster_specific.pdf'), width = 13, height = 10)
    grid.arrange(arrangeGrob(grobs=ff.plots, ncol=2, top=mytitle, left= myleft))
    dev.off()
    
    
    ## plot cell type specific genes 
    tmp$index <- paste0('index_', 1:nrow(tmp))
    keep <- NULL
    for (c in c('C1vsCO', 'C2vsCO')) {
      
      for (j in 1:nrow(tmp[tmp$comp==c,])) {
        gg <- tmp$GeneName[tmp$comp==c][j]
        mm = tmp[tmp$GeneName==gg & tmp$comp==c, c('cluster', 'estimate', 'index')]
        mm = mm[order(mm$estimate, decreasing = T), ]
        
        if ( ((mm$estimate[1] - mm$estimate[2])/mm$estimate[1] > 0.5) | nrow(mm) ==1 ) {
          keep <- c(keep, mm$index[mm$estimate == mm$estimate[1]])
        } else {
          next
        }
        
      }
    }
    
    ## extract cell type specific genes 
    tmp <- tmp[tmp$index %in% keep, ]
    gene.counts <- aggregate(GeneName ~ cluster + comp, data = tmp, FUN = length)
    colnames(gene.counts) <- c('cluster', 'comp', 'num.genes')
    totals <- aggregate(. ~ comp, data = gene.counts, sum )
    gene.counts$pct <- 'NA'
    gene.counts[gene.counts$comp=='C1vsCO', 'pct'] <- (gene.counts$num.genes[gene.counts$comp=='C1vsCO'])/(totals$num.genes[totals$comp=='C1vsCO'])
    gene.counts[gene.counts$comp=='C2vsCO', 'pct'] <- (gene.counts$num.genes[gene.counts$comp=='C2vsCO'])/(totals$num.genes[totals$comp=='C2vsCO'])
    gene.counts$pct <- as.numeric(gene.counts$pct)
    gene.counts$comp <- gsub('vsCO', '', gene.counts$comp)
    
    
    p <-  ggplot(gene.counts, aes(x=factor(comp), y= pct)) + theme_bw(base_size = 12) + 
          geom_bar(stat="identity", aes(fill=factor(cluster))) + coord_polar(theta = "y" ) + 
          #geom_jitter(size=1.5, alpha= 0.5, shape=1, width=0.3, aes(color=factor(cluster)), show.legend = F) +
          labs(x='', y = 'Percentage of DE genes', title = paste0('Percentage of cell-type hits per cell-type (', d,')')) +
          theme(axis.text.x=element_text(size=14, vjust=0.5, face="bold"),
                axis.text.y=element_text(size=14, face="bold"),
                axis.title.y=element_text(size=16, face="bold"),
                plot.title = element_text(size = 18, hjust=0.5, face="bold"),
                axis.ticks=element_blank(), panel.grid.minor = element_blank(),
                legend.title=element_blank(), legend.background = element_rect(color = "black"),
                panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) +
          scale_fill_brewer(palette = 'Dark2' )
    
    pdf(paste0('omics_integration/Replication/MSBB/BM36/2Clusters/cell_specific_pathways/percentage_hits_per_cell_type_', d, '_top10.pdf'), width = 10, height = 7)
    p
    dev.off()
}
######################################################################################################################################


########################################################
### check synaptic genes 
########################################################
data.name <- 'SynGO'
if (data.name == 'Synaptome') {
  synaptic_genes <- as.data.frame(read_excel("omics_integration/data/Synaptome_updated 10_31_2017.xlsx", sheet =1))  
  dir.create(paste0('omics_integration/Replication/MSBB/BM36/2Clusters/Overlap_with_', data.name)) 
} else if (data.name == 'SynGO') {
  ## from SynGo
  synaptic_genes <- as.data.frame(read_excel("omics_integration/data/SynGO_Data/SynGO_bulk_download_release_20210225/syngo_annotations.xlsx", sheet =1))
  colnames(synaptic_genes)[colnames(synaptic_genes)=="hgnc_symbol"] <- "Symbol"
  dir.create(paste0('omics_integration/Replication/MSBB/BM36/2Clusters/Overlap_with_', data.name)) 
}

comps <- c('C1vsCO','C2vsCO','C1vsC2')
for (cmp in comps) {
  cat('Running for:', cmp, '\n')
  #dir.create(paste0('omics_integration/iCluster_output/Figures/Overlap_with_SynGO/', cmp))
  
  #"Transcriptomics"
  C4.res <- read.table(paste0('omics_integration/Replication/MSBB/BM36/2Clusters/DE/',cmp,'/de_res_sig_',cmp,'.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')
  all.gg <- read.table(paste0('omics_integration/Replication/MSBB/BM36/2Clusters/DE/',cmp,'/de_res_',cmp,'.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')[,2]
  
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
  
  myCatNames <- c(paste0("MSBB (",cmp,")\n(",length(C4.res$GeneName),")"), paste0(data.name,"\n(",length(unique(synaptic_genes$Symbol)),")"))
  myCompGroup <- list(C4.res$GeneName, synaptic_genes$Symbol)
  #mytitle <- 'Up-regulated Genes'
  kk <- C4.res$GeneName
  
  ## write results 
  shared.genes <- C4.res[C4.res$GeneName %in% intersect(kk, synaptic_genes$Symbol), ]
  shared.genes <- shared.genes[order(shared.genes$log2FoldChange), ]
  
  ## run hypergeometirc test 
  test.mat <- matrix(c( length(unique(all.gg)) - length(union(kk, unique(synaptic_genes$Symbol))), length(sediff(kk, unique(synaptic_genes$Symbol))), 
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
  pdf(paste0('omics_integration/Replication/MSBB/BM36/2Clusters/Overlap_with_',data.name,'/Overlap_with_',data.name,'_',cmp,'.pdf'), width=5, height = 5)
  grid.draw(vv)
  dev.off()
  
  #  }
}
