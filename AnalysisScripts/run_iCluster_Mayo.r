source("omics_integration/scripts/plotHMBayes.r")
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
library(readxl)

dir.create('omics_integration/Replication/Mayo/')
## read clinical data 
mayo.clin <- read.csv('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201606_MayoRNAseq/03.-Phenotype/2019_09_25_Mayo_RNAseq_clinical.csv', sep=",", header =T, stringsAsFactors = F)
mayo.clin <- mayo.clin[mayo.clin$Status %in% c('AD', 'Control'),]
mayo.tech <- read.csv('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201606_MayoRNAseq/03.-Phenotype/2020_10_20_Mayo_technical.csv', sep=",", header =T, stringsAsFactors = F)
## add IDs 
mayo.sub_ids <- as.data.frame(str_split_fixed(mayo.tech$Sample.Name, "_", 2))
mayo.tech <- cbind(mayo.tech, mayo.sub_ids[,1])
colnames(mayo.tech)[ncol(mayo.tech)] <- 'Subj_ID'
## filter AD and CO cases 
mayo.tech <- mayo.tech[mayo.tech$Subj_ID %in% mayo.clin$Subj_ID, ]
sapply(list(mayo.clin, mayo.tech), dim)

########################################################################################################################
mayo.region <- "temporal_cortex"
dir.create(paste0('omics_integration/Replication/Mayo/',mayo.region))
dir.create(paste0('omics_integration/Replication/Mayo/',mayo.region,'/Figures'))
dir.create(paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results'))

#### read expression & clinical files 
path <- paste0('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201606_MayoRNAseq/',mayo.region,'/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies')
## read count & exp data
mayo.cts <- as.data.frame(fread(paste0(path,'/gene_quant_matrix_NumReads.tsv'), header = 'auto', stringsAsFactors = F, check.names = F))
mayo.exp <- as.data.frame(fread(paste0(path, '/gene_quant_matrix_TPM.tsv'), header = 'auto', stringsAsFactors = F, check.names = F))
exp.meta.cols <- mayo.exp[, 1:3]

## read clinical data 
## extract region samples 
reg.clin <- mayo.tech[mayo.tech$Tissue == mayo.region, ]
reg.clin <- merge(mayo.clin, reg.clin[, c('Subj_ID', 'Sample.Name')])
reg.ad <- reg.clin$Sample.Name[reg.clin$Status == "AD"]
reg.co <- reg.clin$Sample.Name[reg.clin$Status=="Control"]
sapply(list(reg.ad, reg.co), length)

## keep region samples only 
reg.exp <- mayo.exp[, colnames(mayo.exp) %in% c('GeneID', 'GeneName', reg.ad)]
## make rownames 
mean.expr = rowMeans(reg.exp[, -c(1,2)], na.rm = T)
reg.exp = reg.exp[order(mean.expr, decreasing=T),]
reg.exp = reg.exp[!duplicated(reg.exp[["GeneName"]]),]
rownames(reg.exp) <- reg.exp$GeneName
reg.exp$GeneID <- NULL
reg.exp$GeneName <- NULL
dim(reg.exp)

## QC 
pdf(paste0('omics_integration/Replication/Mayo/',mayo.region,'/Figures/PCA_plot.pdf'), width =10, height = 10)
plotPCA(as.matrix(reg.exp), labels = F, cex=1.5, col = 'red', main=paste0('Mayo - ',mayo.region), cex.main = 3)
dev.off()
#acc.outliers <- c('SM-CSEEA.final')
## remove outlines from exp, count, and clincial data 
#reg.exp <- reg.exp[, !colnames(reg.exp) %in% acc.outliers]
#reg.clin <- reg.clin[!reg.clin$Sample.Name %in% acc.outliers, ]

## ------------------------------------------------------------------
## read proteomic data 
mayo.prt.data <- readRDS('omics_integration/data/Mayo/Mayo_Proteomics_TC_proteinoutput_corrected.rds')
mayo.prt.pheno <- readRDS('omics_integration/data/Mayo/Mayo_Proteomics_TC_clinical_data.rds')
## keep only shared samples 
shared.samples <- mayo.prt.pheno$data.cols[mayo.prt.pheno$individualID %in% intersect(reg.clin$Subj_ID, mayo.prt.pheno$individualID)]
mayo.prt.data <- mayo.prt.data[, shared.samples]
mayo.prt.pheno <- mayo.prt.pheno[shared.samples, ]
sapply(list(mayo.prt.data, mayo.prt.pheno), dim)

## remove possible outliers 
# p.outliers <- c('b1_141_04', 'b3_030_06', 'b3_023_04', 'b3_003_07')
# mayo.prt.data <- mayo.prt.data[, !colnames(mayo.prt.data) %in% p.outliers]
# mayo.prt.pheno <- mayo.prt.pheno[!rownames(mayo.prt.pheno) %in% p.outliers, ]
# sapply(list(mayo.prt.data, mayo.prt.pheno), dim)

## replace columns names 
cols <- colnames(mayo.prt.data)
for (k in 1:length(cols)){
  col.id <- cols[k]
  ind.id <- mayo.prt.pheno$individualID[mayo.prt.pheno$data.cols==col.id]
  sample.name <- mayo.tech$Sample.Name[mayo.tech$Subj_ID == ind.id & mayo.tech$Tissue=="temporal_cortex"]
  colnames(mayo.prt.data)[colnames(mayo.prt.data)==col.id] <- sample.name
}
# extract shared samples 
shared.samples <- intersect(colnames(reg.exp), colnames(mayo.prt.data))
mayo.prt.data <- mayo.prt.data[, shared.samples]
reg.exp <- reg.exp[, shared.samples]
sapply(list(reg.exp, mayo.prt.data), dim)
all(colnames(reg.exp) %in% colnames(mayo.prt.data))

# ## check missing values 
# mayo.prt.data$na_count <- apply(mayo.prt.data, 1, function(x) sum(is.na(x)))
# mayo.prt.data <- mayo.prt.data[mayo.prt.data$na_count < length(colnames(mayo.prt.data)[colnames(mayo.prt.data) !='na_count']) * 0.20, ]
# mayo.prt.data$na_count <- NULL
# ## replace all NA with the average value
# for(i in 1:nrow(mayo.prt.data)){
#   indx = which(is.na(mayo.prt.data[i,]))
#   row.val <- rowMeans(mayo.prt.data[i,], na.rm = T)
#   mayo.prt.data[i, indx] <- row.val
# }
# prt.exp <- mayo.prt.data
## ------------------------------------------------------------------

#

## filter features based on discovery cluster 4 + cluster DE genes  
C4vsALL.res <- read.table('omics_integration/iCluster_output/DE/C4vsC123/de_res_sig_C4vsC123.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')
C2vsALL.res <- read.table('omics_integration/iCluster_output/DE/C2vsC134/de_res_sig_C2vsC134.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')

C4vsALL.res.prt <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_4vsall.csv', header =T, stringsAsFactors = F, check.names = F, sep=',')
C2vsALL.res.prt <- read.csv('omics_integration/data/Proteomics/DE_results_proteomics/prot_res_2vsall.csv', header =T, stringsAsFactors = F, check.names = F, sep=',')
gg.prt <- c(C4vsALL.res.prt$EntrezGeneSymbol[C4vsALL.res.prt$padj < 0.05], 
            C2vsALL.res.prt$EntrezGeneSymbol[C2vsALL.res.prt$padj < 0.05])
gg.prt <- unlist(strsplit(gg.prt, " ","", ))
gg.prt <- unlist(strsplit(gg.prt, ",","", ))

## extract only genes DE in C4vsALL/C2vsALL for transcriptomics 
gg <- c(C4vsALL.res$GeneName, C2vsALL.res$GeneName)
reg.exp <- reg.exp[rownames(reg.exp) %in% gg, ]

## extract only shared proteins 
prt.genes <- read.table('omics_integration/data/Mayo/gene_names.txt', header =T, sep="\t", stringsAsFactors = F, check.names = F)
prt.genes$GeneNames <- gsub("_HUMAN", "", prt.genes$GeneNames)
shared.prt.geens <- unique(unlist(strsplit(prt.genes$GeneNames, "\\|", "")))
shared.prt.geens <- shared.prt.geens[!shared.prt.geens %in% c('>sp', '>tr')]
shared.prt.geens <- intersect(gg.prt, shared.prt.geens)
keep <- NULL
for (row in 1:nrow(prt.genes)) {
  g.list <- prt.genes$GeneNames[row]
  if (is.na(g.list)) { next }
  
  if (any(grepl(g.list, shared.prt.geens) == TRUE)) {
          keep <- c(keep, row)
  } 
}
prt.genes <- prt.genes[keep, ]
dim(prt.genes)
mayo.prt.data <- mayo.prt.data[rownames(mayo.prt.data) %in% prt.genes$MajorityProteinIDs, ]
## check both datasets
sapply(list(reg.exp, mayo.prt.data), dim)

# filter features with no variance at all
reg.exp <- reg.exp[apply(reg.exp, 1, var) > 0,]
mayo.prt.data <- mayo.prt.data[apply(mayo.prt.data, 1, var) > 0,]
sapply(list(reg.exp, mayo.prt.data), dim)

################################################################################
############## run tune iclusterBayes     
################################################################################
seed.val <- 1020
set.seed(seed.val) ## 357 the first best solution 
MAX.NUM.CLUSTERS = 11
num.omics <- 2
icluster.res = tune.iClusterBayes(cpus=(MAX.NUM.CLUSTERS - 1), 
                                  t(log2(reg.exp+1)), t(mayo.prt.data),
                                  K=1:(MAX.NUM.CLUSTERS - 1),
                                  type=rep('gaussian', num.omics),  # poisson
                                  #n.burnin=12000, n.draw=18000,
                                  prior.gamma=rep(0.5, num.omics),
                                  pp.cutoff = 0.5,
                                  sdev=0.05,
                                  thin=3
)$fit
# save the result object
saveRDS(icluster.res, file=paste0("omics_integration/Replication/Mayo/",mayo.region,"/iCluster.res.RNAseqProtein.rds"))
################################################################################

## extract dev.ratio's & BIC 
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

pdf(paste0('omics_integration/Replication/Mayo/',mayo.region,'/Figures/k_vs_devRatio.pdf'), width = 12, height = 5)
grid.arrange(p1, p2, ncol=2)
dev.off()

k=which.max(dev.ratios)
pdf(paste0('omics_integration/Replication/Mayo/',mayo.region,'/Figures/Posterior_probability_dist.pdf'), width = 8, height = 5)
par(mfrow=c(2,1))
plot(icluster.res[[k]]$beta.pp[[1]], xlab="Genes", ylab="Posterior probability", main="RNA-Seq expression")
abline(h = 0.5, v = 0, col = "red")
plot(icluster.res[[k]]$beta.pp[[2]], xlab="Genes", ylab="Posterior probability", main="Protein expression")
abline(h = 0.5, v = 0, col = "red")
dev.off()

reg.exp.image=scale(log2(t(reg.exp+1)))
reg.exp.image[reg.exp.image > 2.5] = 2.5
reg.exp.image[reg.exp.image < -2.5] = -2.5

mayo.prt.image=scale(t(mayo.prt.data))
mayo.prt.image[mayo.prt.image > 2.5] = 2.5
mayo.prt.image[mayo.prt.image < -2.5] = -2.5

#col.scheme = alist()
#col.scheme = rev(colorRampPalette(brewer.pal(10, "RdBu"))(300))
#col.scheme[[2]] = rev(colorRampPalette(brewer.pal(10, "RdBu"))(300))
pdf(paste0('omics_integration/Replication/Mayo/',mayo.region,'/Figures/iClusterPlus_Heatmap_RNAseqProtein.pdf'), width=10, height = 5)
MyplotHMBayes(fit=optimal.solution, 
            datasets= list(reg.exp.image, mayo.prt.image), 
            type=rep("gaussian",num.omics),
            scale = rep(F,num.omics), 
            #col.scheme = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
            threshold=c(0.5, 0.005),
            row.order=rep(T,num.omics),  
            sparse=rep(T,num.omics),
            cap=rep(T,num.omics),
            title = mayo.region,
            dd.types = c(paste0('(RNA-seq (pp=',0.5,')'), paste0('(Protein (pp=',0.005,')'))
)
dev.off()

########################################################
## Extract cluster membership of the best class 
########################################################
all.clusters <- matrix(data= NA, nrow= nrow(t(reg.exp)), ncol=MAX.NUM.CLUSTERS - 1)
for (k in 1:(MAX.NUM.CLUSTERS - 1)) {
  cc = icluster.res[[k]]$clusters
  cc = as.matrix(cc)
  all.clusters[, k] = cc
}
rownames(all.clusters) = colnames(reg.exp)
colnames(all.clusters) = paste0('K', 1:(MAX.NUM.CLUSTERS - 1))
all.clusters = as.data.frame(all.clusters)
all.clusters$Sample.Name = rownames(all.clusters)
rownames(all.clusters) = NULL
all.clusters = all.clusters[, c('Sample.Name',  colnames(all.clusters)[colnames(all.clusters) !="Sample.Name"])]
## get best cluster membership
best.cl= unique(c(which.max(dev.ratios), which.min(allBICs)))
if (length(best.cl) > 1) {
  stop('Error: something is wrong!!')
}
best.cluster.membership <- all.clusters[, c('Sample.Name',  paste0('K', best.cl))]
colnames(best.cluster.membership) <- c('Sample.Name', 'best.cluster')
best.cluster.membership <- best.cluster.membership[order(best.cluster.membership$best.cluster), ]
## add controls 
cc.data <- data.frame(Sample.Name= reg.co, best.cluster= 0)
best.cluster.membership <- rbind(best.cluster.membership, cc.data)
best.cluster.membership <- merge(reg.clin, best.cluster.membership)
#colnames(best.cluster.membership)[colnames(best.cluster.membership) == "projid"] <- 'Subj_ID'
write.table(best.cluster.membership, file=paste0('omics_integration/Replication/Mayo/',mayo.region,'/best.cluster.membership.tsv'), sep="\t", quote = F, row.names = F)

########################################################
##### plot phenotype association 
########################################################
bb <- read.table(paste0('omics_integration/Replication/Mayo/',mayo.region,'/best.cluster.membership.tsv'), header =T, sep="\t", stringsAsFactors = F, check.names = F)
bb$AOD <- as.numeric(gsub('_or_above','', bb$AOD))
## merge with other clinical data 
colnames(mayo.prt.pheno)[colnames(mayo.prt.pheno) == "individualID"] <- "Subj_ID"
bb <- merge(bb, mayo.prt.pheno, all.x = TRUE)
bb$ageDeath <- as.numeric(gsub('_or_above','', bb$ageDeath))

## compute pvalue 
ff <- bb[bb$best.cluster %in% c('1', '2'), ]

glm.res <- glm(Braak ~ best.cluster + Sex + AOD, data = ff, family = 'gaussian')
pval <- data.frame(PValue= signif(coef(summary(glm.res))[,4][["best.cluster"]], digits = 3), 
                   y_pos = max(ff$Braak, na.rm= T), best.cluster = NA)

## for AOD 
#pval <- t.test(bb$AOD[bb$best.cluster==1], bb$AOD[bb$best.cluster==2])$p.value
#pval <- data.frame(PValue= signif(pval, digits = 3), y_pos = max(bb$AOD), best.cluster = NA)

ph.name <- 'APOE4_Status'

## for APOE 
bb$APOE4_Status <- 'APOE4-'
bb[bb$APOE %in% unique(bb$APOE[grepl(4, bb$APOE)]), 'APOE4_Status'] <- 'APOE4+'

p <- ggplot(bb, aes(x=factor(best.cluster), y=as.factor(APOE4_Status), color=factor(best.cluster))) + 
      geom_point(size = 4, shape=21, position=position_jitter(width=0.2, height=0.1)) + theme_bw() + 
      #geom_jitter(position=position_jitter(0.3), size=1.5, aes(color = factor(best.cluster))) + theme_bw() +
      labs(x='', y=ph.name) + ggtitle(ph.name) + scale_color_brewer(palette = "Dark2") + 
      theme(axis.text.x=element_text(size=12, vjust=0.5, face="bold", color="black"),
        axis.text.y=element_text(size=12, color="black", face="bold"),
        #axis.title.x=element_text(size=19, face="bold"),
        axis.title.y=element_text(size=18, face="bold"),
        #panel.background=element_rect(color="black"),
        plot.title = element_text(size = 22, hjust=0.5, color="black", face="bold"),
        legend.position="none", panel.border = element_rect(linetype='solid', color='black'),
        plot.margin=unit(c(1,1,1,5), 'mm')) +
      scale_x_discrete(labels=c("0"=paste0("Control\n(n=",length(bb$Subj_ID[bb$best.cluster==0]),')'),
                            "1"=paste0("Cluster1\n(n=",length(bb$Subj_ID[bb$best.cluster==1]),')'), 
                            "2"= paste0("Cluster2\n(n=",length(bb$Subj_ID[bb$best.cluster==2]),')'))) +
      geom_signif(data = pval, aes(xmin = 2, xmax = 3, annotations = paste0('p=',PValue), y_position = y_pos), 
              textsize = 7, size=0.6, step_increase=0.1, vjust = -0.2,manual = TRUE, inherit.aes = FALSE) 

pdf(paste0('omics_integration/Replication/Mayo/',mayo.region,'/Figures/phen_assoc_', ph.name,'.pdf'), width = 6, height = 8, useDingbats = F)
p
dev.off()

################################################################
#### Check cell proportions 
###############################################################
bb <- read.table(paste0('omics_integration/Replication/Mayo/',mayo.region,'/best.cluster.membership.tsv'), header =T, sep="\t", stringsAsFactors = F, check.names = F)
bb$AOD <- as.numeric(gsub('_or_above','', bb$AOD))

algs = c( "ssNMF", "meanProfile")
cellProp.res <- NULL
for (alg in algs) {
  s.res <- read.table(paste0('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201606_MayoRNAseq/05.-Analyses/deconvolution/deconvolution_', alg,'.tsv'), header = T, check.names = F, stringsAsFactors = F)
  cellProp.res <- rbind(cellProp.res, s.res)
}
# merge with the best clusters 
colnames(cellProp.res)[colnames(cellProp.res)=="Sample"] <- "Sample.Name"
cellProp.res <- reshape2::melt(cellProp.res, id.vars = c("Sample.Name", "Algorithm"))
colnames(cellProp.res) <- c('Sample.Name', 'Algorithm', 'Cell_type', 'Proportion')
# merge with clusters 
cellProp.res <- merge(cellProp.res, bb)
cellProp.res <- cellProp.res[cellProp.res$best.cluster !=0, ]
dim(cellProp.res)

## compute pvalues using GLM
all.pvals <- NULL
for (a in c( "ssNMF", "meanProfile")) {
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

for (alg in c( "ssNMF", "meanProfile")) {
  all.plots <- list()
  for (cell in as.character(unique(cellProp.res$Cell_type)) ) {
    dd <- all.pvals[all.pvals$Algorithm==alg & all.pvals$Cell_type==cell, ]
    p <- ggplot(cellProp.res[cellProp.res$Algorithm==alg & cellProp.res$Cell_type==cell, ], 
                aes(x=factor(best.cluster), y=Proportion, group=factor(best.cluster))) + 
      geom_boxplot(aes(fill=factor(best.cluster))) +  theme_bw() +  
      scale_color_brewer(palette = "Dark2") + labs(x='', y='') + ggtitle(cell) + 
      theme(axis.text.x=element_text(size=18, vjust=1, color="black", hjust=1, angle= 45),
            axis.text.y=element_text(size=18, color="black"), 
            plot.title = element_text(size = 20, hjust=0.5, color="black"),
            legend.position="none", panel.border = element_rect(linetype='solid', color='black')) +
      scale_x_discrete(labels=c("1"=paste0("Cluster1\n(n=",length(unique(bb$Sample.Name[bb$best.cluster==1])),')'), 
                                "2"= paste0("Cluster2\n(n=",length(unique(bb$Sample.Name[bb$best.cluster==2])),')'))) +
      geom_signif(data = dd, aes(xmin = 1, xmax = 2, annotations =  paste0('p=', Pvalue), y_position = y_pos), 
                  textsize = 6, step_increase=0.1, vjust = -0.2,manual = T, fontcolor="bold", margin_top=0.2) 
    
    all.plots[[cell]] <- p
    
  }
  ## print 
  myleft <- textGrob('Cell Proportion', gp=gpar(fontsize=22), rot=90, vjust = 0.5)
  mytop <- textGrob(paste0('May - ',mayo.region,' (', alg,')\n'), gp=gpar(fontsize=22))
  pdf(paste0("omics_integration/Replication/Mayo/",mayo.region,"/Figures/deconv_",alg,".pdf"), width = 12, height = 8, useDingbats = F)
  grid.arrange(grobs=all.plots, nrow=1, ncol=4, left = myleft, top = mytop)
  dev.off()
}

##############################################################
##################### Survival Association ###################
##############################################################
library(surviplot)
library(survival)

clin = bb[, c('Subj_ID', 'AOD', 'Sex', 'best.cluster')]
#clin = clin[clin$AOD > 70, ]
colnames(clin) = c('sample', 'time', 'Sex', 'best.cluster')

## add clusters group 
clin$group = 'NA'
clin[clin$best.cluster ==0, 'group'] <- 'CO'
clin[clin$best.cluster ==1, 'group'] <- 'C1'
clin[clin$best.cluster ==2, 'group'] <- 'C2'

max.time = max(clin$time, na.rm = T)

#groups = c('Cluster1', 'Cluster2', 'Cluster3', 'Cluster4')
groups = unique(clin$group)

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


## plot KM curve 
pdf(paste0('omics_integration/Replication/Mayo/', mayo.region, '/Figures/survival_proportion.pdf'), width=5, height = 10)
par(mfrow=c(3,1))
surviplot(Surv(time) ~ group, data=clin[clin$group !='CO', ], ylab='Survival Proportion', xlim=c(50,max.time), 
          main ="Overall Survival (C2vsC1)", xlab='Age of death (years)', cex.main=1.5,
          mark.time=TRUE, col=c("skyblue3", 'orange2'), lwd = 3, cex.lab = 1.5, cex.axis = 1.5)
surviplot(Surv(time) ~ group, data=clin[clin$group !='C1', ], ylab='Survival Proportion', xlim=c(50,max.time), 
          main ="Overall Survival (C2vsCO)", xlab='Age of death (years)', cex.main=1.5,
          mark.time=TRUE, col=c("skyblue3", 'orange2'), lwd = 3, cex.lab = 1.5, cex.axis = 1.5)
surviplot(Surv(time) ~ group, data=clin[clin$group !='C2', ], ylab='Survival Proportion', xlim=c(50,max.time), 
          main ="Overall Survival (C1vsCO)", xlab='Age of death (years)', cex.main=1.5,
          mark.time=TRUE, col=c("skyblue3", 'orange2'), lwd = 3, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

################################################################
################## DE Analysis - RNA-Seq ####################### 
################################################################
## read raw counts 
reg.counts <- mayo.cts
meta.cols <- reg.counts[, 1:3]
rownames(reg.counts) <- paste0(reg.counts$GeneID,'_', reg.counts$GeneName)
reg.counts$GeneID <- NULL
reg.counts$GeneName <- NULL
reg.counts$GeneBiotype <- NULL
reg.counts <- round(reg.counts)

## cell proportion  
cell.res <- read.table('/home/data/Public_Data/bulkRNASeq_Data/Hg38/201606_MayoRNAseq/05.-Analyses/deconvolution/deconvolution_ssNMF.tsv', header = T, check.names = F, stringsAsFactors = F)
colnames(cell.res)[colnames(cell.res)=="Sample"] <- "Sample.Name"
bb.de <- merge(bb, cell.res[, c('Sample.Name', 'Astrocyte', 'Neuron', 'Oligodendrocyte')])
## organize data 
rownames(bb.de) <- bb.de$Sample.Name
c1.samples <- bb.de$Sample.Name[bb.de$best.cluster==1]
c2.samples <- bb.de$Sample.Name[bb.de$best.cluster==2]
reg.counts <- reg.counts[, c(reg.co, c1.samples, c2.samples)]
bb.de <- bb.de[c(reg.co, c1.samples, c2.samples), ]
all(rownames(bb.de) == colnames(reg.counts) )
sapply(list(bb.de, reg.counts), dim)

## function to filter lowly expressed genes 
selectGenes <- function(grp1, grp2) {
  norm.cts <- cpm(reg.counts)
  grp1.t <- round(length(grp1) * 0.50)
  grp2.t <- round(length(grp2) * 0.50)
  keep <- (rowSums (norm.cts[,grp1]> 0.5) >= grp1.t) | (rowSums (norm.cts[,grp2]> 0.5) >= grp2.t)
  #keep <- (rowSums (rosmap.cts_norm[,grp1]> 0.5) >= grp1.t) | (rowSums (rosmap.cts_norm[,grp2]> 0.5) >= grp2.t)
  if (all(rownames(reg.counts) != rownames(norm.cts))) {
    stop('rownames of the normalized counts are not the same as the raw counts ')
  } else {
    cmp.counts <- reg.counts[keep, ]
  }
  return (cmp.counts)
}

## Function for fitting the DESeq2 model
bb.de$AOD <- round(as.numeric(bb.de$AOD))
bb.de$Sex <- bb.de$Sex
bb.de$best.cluster <- as.factor(bb.de$best.cluster)
DESeqModelFit <- function(count_data) {
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = bb.de,
                                design = formula(~ Sex + AOD + Astrocyte + Neuron + best.cluster))
  ## Test for DE 
  dds <- DESeq(dds)
  print(resultsNames(dds))
  return (dds)
}

## enrichR
dbs <- listEnrichrDbs()
run_enrichr <- function (de, cc) {
  #gene.list <- unlist(strsplit(rownames(res_sig), "_"))
  #gene.list <- gene.list[!grepl('ENSG', gene.list)]
  en.up <- enrichr(de$GeneName[de$direction =="up"], databases = dbs$libraryName)  
  en.dn <- enrichr(de$GeneName[de$direction =="dn"], databases = dbs$libraryName)  
  save(en.up, file=paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/', cc,'/', cc,'_enrichR_res_up.RData'))
  save(en.dn, file=paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/', cc,'/', cc,'_enrichR_res_dn.RData'))
}

## run DE
comps <- c('1vs0', '2vs0', '2vs1')
for (cmp in comps) {
  cat('Running for:', cmp, '\n')
  dir.create(paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/', cmp))
  
  ## extract groups 
  #mygroup <- "best.cluster"
  grp1.samples <- bb.de$Sample.Name[bb.de$best.cluster== unlist(strsplit(cmp, "vs"))[1] ]
  grp2.samples <- bb.de$Sample.Name[bb.de$best.cluster== unlist(strsplit(cmp, "vs"))[2] ]
  
  ## remove lowly expressed genes 
  grp.counts <- selectGenes(grp1.samples, grp2.samples)
  dim(grp.counts)
  
  ## fit DESeq model
  dds <- DESeqModelFit(grp.counts)
  ## extract result table
  dds.res <- results(dds, alpha=0.05, contrast = c('best.cluster', unlist(strsplit(cmp, "vs"))[1], unlist(strsplit(cmp, "vs"))[2] ), tidy = F)
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
  normalized_counts <- cpm(grp.counts)
  normalized_counts = normalized_counts[paste0(dds.res$GeneID,'_', dds.res$GeneName), ]
  dds.res$grp1 <- rowMeans(normalized_counts[, grp1.samples], na.rm=TRUE)
  dds.res$grp2 <- rowMeans(normalized_counts[, grp2.samples], na.rm = TRUE)
  colnames(dds.res)[colnames(dds.res) =="grp1"] <- paste0('grp.', unlist(strsplit(cmp, "vs"))[1], '.ExpMean')
  colnames(dds.res)[colnames(dds.res) =="grp2"] <- paste0('grp.', unlist(strsplit(cmp, "vs"))[2], '.ExpMean')
  
  ## plot volcano plot
  top.10.genes <- dds.res$GeneName[1:10]
  pdf(paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/',cmp, '/volcano_plot_',cmp,'.pdf'), width = 6, height = 8)
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
  
  write.table(dds.res, file=paste0('omics_integration/Replication/Mayo/', mayo.region,'/DE Results/',cmp,'/de_res_', cmp,'.tsv'), sep="\t", quote = F, row.names = F)
  write.table(dds.res.sig, file=paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/',cmp, '/de_res_sig_',cmp,'.tsv'), sep="\t", quote = F, row.names = F)
  
  ## run pathway analysis
  #run_kegg_go(dds.res.sig, cmp, 0)
  run_enrichr(dds.res.sig, cmp)
  
}

#####################################################
####### check the overlap with discovery
#####################################################
cmp <- '2vs1'
for (d in c('up', 'dn')) {
  
  Discovery <- read.table('omics_integration/iCluster_output/DE/C4vsC123/de_res_sig_C4vsC123.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')
  Discovery <- Discovery$GeneName[Discovery$direction==d]
  
  Replicate <- read.table(paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/',cmp, '/de_res_sig_',cmp,'.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')
  Replicate <- Replicate$GeneName[Replicate$direction==d]
  
  ## run hypergeometirc test 
  all.genes <- union(read.table('omics_integration/iCluster_output/DE/C4vsC123/de_res_C4vsC123.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t')[,2], 
                     read.table(paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/',cmp, '/de_res_sig_',cmp,'.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')[,2])
  
  #all.genes <- nrow(msbb.cts)
  test.mat <- matrix(c( length(all.genes) - length(union(Discovery, Replicate)), length(setdiff(Discovery, Replicate)), length(setdiff(Replicate, Discovery)), length(intersect(Discovery, Replicate))), nrow=2)
  fisher.pval <- fisher.test(test.mat, alternative = "two.sided")$p.value
  if (fisher.pval ==0) { fisher.pval ='P-value < 2.2e-16'} else {fisher.pval = paste0('P-value = ',signif(fisher.pval, digits = 4)) }
  
  ## generate venn diagram 
  fileName <- paste0('Knight-C4_vs_Mayo-', mayo.region,'-C2')
  myCatNames <- c(paste0("Knight ADRC\n(",length(Discovery),")") , paste0("Mayo-TCX\n(",length(Replicate),")"))
  myCompGroup <- list(Discovery, Replicate)
  tt <- ifelse (d == "up", 'Up-regulated', 'Down-regulated') 
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  myCol <- brewer.pal(2, "Dark2")
  vv <- venn.diagram( 
    x = myCompGroup,
    category.names = myCatNames, 
    filename = NULL, main = paste0(tt,'\n(',fisher.pval,')\n'), main.fontface = "bold", 
    height = 300, width = 300 , resolution = 300, compression = "lzw",
    lwd = 2, lty = 'blank', main.cex = 0.8, 
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
  
  pdf(paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/',cmp, '/', fileName, '_venn_', d,'.pdf'), width=2.8, height=2.8)
  grid.draw(vv)
  dev.off()
  
}

########################################################
### overlap with synaptic genes 
########################################################
data.name <- 'SynGO'
if (data.name == 'Synaptome') {
  synaptic_genes <- as.data.frame(read_excel("omics_integration/data/Synaptome_updated 10_31_2017.xlsx", sheet =1))  
} else if (data.name == 'SynGO') {
  synaptic_genes <- as.data.frame(read_excel("omics_integration/data/SynGO_Data/SynGO_bulk_download_release_20210225/syngo_annotations.xlsx", sheet =1))
  colnames(synaptic_genes)[colnames(synaptic_genes)=="hgnc_symbol"] <- "Symbol"
}

comps <- c('1vs0','2vs0','2vs1')
for (cmp in comps) {
  cat('Running for:', cmp, '\n')
  #dir.create(paste0('omics_integration/iCluster_output/Figures/Overlap_with_SynGO/', cmp))
  #"Transcriptomics"
  clu.res <- read.table(paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/',cmp,'/de_res_sig_',cmp,'.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')
  clu.res <- clu.res[clu.res$GeneBiotype == "protein_coding",]
  all.gg <- read.table(paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/',cmp, '/de_res_',cmp,'.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')
  all.gg <- all.gg$GeneName[all.gg$GeneBiotype == "protein_coding"]
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
  
  myCatNames <- c(paste0("Mayo (",cmp,")\n(",length(clu.res$GeneName),")"), paste0(data.name,"\n(",length(unique(synaptic_genes$Symbol)),")"))
  myCompGroup <- list(clu.res$GeneName, synaptic_genes$Symbol)
  #mytitle <- 'Up-regulated Genes'
  kk <- clu.res$GeneName
  
  ## write results 
  shared.genes <- clu.res[clu.res$GeneName %in% intersect(kk, synaptic_genes$Symbol), ]
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
    main = paste0('\n(Fisher\'s ',fisher.paval,')'), 
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
  pdf(paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/',cmp, '/Overlap_with_',data.name, '.pdf'), width=5, height = 5)
  grid.draw(vv)
  dev.off()
}

########################################################
####### plot genes
########################################################
gene <- 'SNCA' # IGF1, NRXN3, and YWHAZ
## extract RNA-seq exp
g.exp <- as.data.frame(t(mayo.exp[mayo.exp$GeneName==gene, -c(1:3)]))
g.exp$Sample.Name <- rownames(g.exp)
rownames(g.exp) <- NULL
colnames(g.exp) <- c('exp', 'Sample.Name')
g.exp <- merge(g.exp, bb[, c('Sample.Name', 'best.cluster', 'Sex', 'AOD')], sort =F)
g.exp <- g.exp[order(g.exp$best.cluster),]

## extract protein exp
gene.line <- prt.genes$MajorityProteinIDs[prt.genes$GeneNames=="SNCA"]
p.exp <- reshape2::melt(mayo.prt.data[rownames(mayo.prt.data)==gene.line, ])
colnames(p.exp) <- c('Sample.Name', 'prt.exp')
p.exp <- merge(p.exp, bb[, c('Sample.Name', 'best.cluster', 'Sex', 'AOD')], sort =F)
p.exp <- p.exp[order(p.exp$best.cluster),]

## get pvalues from DE analysis 
Tran.C1vsCO <- read.table(paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/1vs0/de_res_1vs0.tsv'), header =T, stringsAsFactors = F, sep="\t")
Tran.C2vsCO <- read.table(paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/2vs0/de_res_2vs0.tsv'), header =T, stringsAsFactors = F, sep="\t")
Tran.C1vsC2 <- read.table(paste0('omics_integration/Replication/Mayo/',mayo.region,'/DE Results/1vs2/de_res_1vs2.tsv'), header =T, stringsAsFactors = F, sep="\t")
#Tran.ADvsCO <- read.table('omics_integration/Replication/ROSMAP/DE/ADvsCO/de_res_ADvsCO.tsv', header =T, stringsAsFactors = F, sep="\t")

pvals <- data.frame(comp=c('C1vsCO','C2vsCO','C1vsC2'), 
                    pval = c(Tran.C1vsCO[Tran.C1vsCO$GeneName==gene, 'padj'], 
                             Tran.C2vsCO[Tran.C2vsCO$GeneName==gene, 'padj'], 
                             Tran.C1vsC2[Tran.C1vsC2$GeneName==gene, 'padj']))
pvals$pval <- format(pvals$pval, scientific = T, digits = 3)

p <- ggplot(g.exp, aes(x=as.factor(best.cluster), y=log2(exp+1))) + 
      geom_boxplot(aes(fill=as.factor(best.cluster)), show.legend = F) +
      labs(x="", y="Log2(Exp)") + ggtitle(paste0(gene, ' - Mayo (',mayo.region,')')) + theme_bw() + 
      #geom_jitter(aes(color=as.factor(best.cluster)), position=position_jitter(0.3), size=1.5) +
      theme(plot.title = element_text(hjust=0.5, size=20, face="bold"),
          axis.text.x = element_text( vjust= 1, size=18, color="black"), 
          #axis.text.x = element_blank(),  
          axis.text.y = element_text(size=18, color="black"),
          axis.title.y = element_text(size=20, face="bold"),
          panel.background = element_rect(colour = "black", size=1), legend.position = "none")  +
      scale_fill_manual(values=c("0"="gray80", "1"="#D95F02", "2"="#E7298A")) +
      scale_x_discrete(labels=c("0"=paste0("Control\n(n=",length(g.exp$Sample.Name[g.exp$best.cluster==0]),')'), 
                            "1"=paste0(toupper(substring(mayo.region, 1,3)),"-C1\n(n=",length(g.exp$Sample.Name[g.exp$best.cluster==1]),')'), 
                            "2"= paste0(toupper(substring(mayo.region, 1,3)), "-C2\n(n=",length(g.exp$Sample.Name[g.exp$best.cluster==2]),')'))) +
      geom_signif(data = pvals, aes(xmin = 1, xmax = 1.97, annotations =  paste0('p=',pval[comp=="C1vsCO"]), y_position = log2(max(g.exp$exp))+0.15), 
              textsize = 7, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, color="black") + 
      geom_signif(data = pvals, aes(xmin = 1, xmax = 3, annotations =  paste0('p=',pval[comp=="C2vsCO"]), y_position = log2(max(g.exp$exp))+0.45), 
              textsize = 7, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, color="black") +
      geom_signif(data = pvals, aes(xmin = 2.03, xmax = 3, annotations =  paste0('p=',pval[comp=="C1vsC2"]), y_position = log2(max(g.exp$exp))+0.15), 
              textsize = 7, size=0.4, step_increase=0.1, vjust = -0.2,manual = T, color="black") 

pdf(paste0('omics_integration/Replication/Mayo/',mayo.region, '/Figures/', gene, '_exp.pdf'), width=6, height = 10)
p
dev.off()

ggplot(p.exp, aes(x=as.factor(best.cluster), y=prt.exp)) + 
  geom_boxplot(aes(fill=as.factor(best.cluster)), show.legend = F)
