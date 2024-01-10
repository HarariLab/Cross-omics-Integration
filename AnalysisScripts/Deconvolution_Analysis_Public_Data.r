setwd('/home/data')
load('/home/eteleeb/scripts/Deconvolution_ml_model.RData')
library(CellMix)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library('NMF')
library(readxl)
library(ggsignif)

cohort <- 'MSBB'
r <- 'cerebellum'

## create analysis directory 
#dir.create('Public_Data/bulkRNASeq/201812_MSBB/BM36/05.-Analyses/deconvolution_analysis')
if (cohort == 'rosmap') {
  dir.create('Public_Data/bulkRNASeq/201506_ROSMAP/dorsolateral_prefrontal_cortex/05.-Analyses/deconvolution_analysis')  
} else if (cohort == "MSBB") {
  dir.create(paste0('Public_Data/bulkRNASeq/201812_MSBB/',r,'/05.-Analyses/deconvolution_analysis'))
} else if (cohoer =="Mayo") {
  #dir.create(paste0('Public_Data/bulkRNASeq_Data/Hg38/201606_MayoRNAseq/05.-Analyses/deconvolution//201812_MSBB/',r,'/05.-Analyses/deconvolution_analysis'))
}

# load in deconvolution algorithm #
#source( 'deconvolute_alg.R')
deconvolute.alg = function ( expression, reference, forceExtractMarkers=F, 
                             selected.algs = c( "qprog", "cs-qprog", "DSA", "ssFrobenius", "meanProfile", "deconf"), 
                             types.reference = NA){
  oldw <- getOption("warn")
  options(warn=-1)
  decon.res = list()
  for( alg in selected.algs) {
    cat( "\n\nAlgorithm ", alg, "\n")
    tryCatch( {
      if ( forceExtractMarkers) {
        #work on this, 
        ml=extractMarkers(reference, types.reference, "HSD")
        res = ged(as.matrix(expression), ml, alg, verbose=T)
      } else {
        res = ged(as.matrix(expression), reference, alg, verbose=T)
      }
      decon.res [[alg]] = scoef(res)
      #decon.res [[alg]] = CellMix::coef(res)
    }, error = function( cond) {}) 
  }
  options(warn=oldw)
  return( decon.res)
}


# load in sample gene expression file #
if (cohort == 'rosmap') {
  exp.tpm <- read.table('Public_Data/bulkRNASeq//201506_ROSMAP/dorsolateral_prefrontal_cortex/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies/gene_quant_matrix_TPM.tsv', header =T, check.names = F, stringsAsFactors = F)
  exp.fpkm <- read.table('Public_Data/bulkRNASeq//201506_ROSMAP/dorsolateral_prefrontal_cortex/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_fpkm_matrix.tsv', header =T, check.names = F, stringsAsFactors = F)
} else if (cohort == 'MSBB'){
  exp.tpm <- read.table(paste0('Public_Data/bulkRNASeq/201812_MSBB/',r,'/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies/gene_quant_matrix_TPM.tsv'), header =T, check.names = F, stringsAsFactors = F)
  exp.fpkm <- read.table(paste0('Public_Data/bulkRNASeq/201812_MSBB/',r, '/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_fpkm_matrix.tsv'), header =T, check.names = F, stringsAsFactors = F)
} else if (cohort == "Mayo") {
  exp.tpm <- read.table(paste0('Public_Data/bulkRNASeq_Data/Hg38/201606_MayoRNAseq/',r,'/02.-ProcessedData/Hg38/04.-Salmon/combined_exp_matricies/gene_quant_matrix_TPM.tsv'), header =T, check.names = F, stringsAsFactors = F)
  exp.fpkm <- read.table(paste0('Public_Data/bulkRNASeq_Data/Hg38/201606_MayoRNAseq/',r,'/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_fpkm_matrix.tsv'), header =T, check.names = F, stringsAsFactors = F)
}

meta.cols <- exp.tpm[, 1:3]
#colnames(exp.fpkm)[colnames(exp.fpkm)=="gene"] <- "GeneID"
#exp.fpkm <- merge(meta.cols, exp.fpkm)
#sapply(list(exp.tpm, exp.fpkm), dim)
exp <- exp.tpm

# load in deconvolution learning model #
#load( file="Deconvolution_ml_model.RData")

# load in gene marker panel #
gene.panel = read.csv('/home/eteleeb/deconvolution_gene_panel_files/Brain_panel_genes_drop_PLP1_CCL2.csv', stringsAsFactors = F, row.names = 1)
gene.panel = gene.panel[, c('cellType', 'Human')]
gene.panel = subset( gene.panel, !is.na(gene.panel$Human) )
# read Celeste list of genes 
gene.panel2 = read.table('/home/eteleeb/deconvolution_gene_panel_files/human_microglia_markers.tsv', header =T, stringsAsFactors = F)
gene.panel = unique(rbind(gene.panel, gene.panel2))
colnames(gene.panel) <- c('cellType', 'gene')
## leave-one-out (e.g. GFAP)
#gene.panel <- gene.panel[!gene.panel$gene %in% c("GFAP", "AQP4"), ]

# extract gene marker expression from the samples  #
gene.panel.expr = exp[exp$GeneName %in% gene.panel$gene,]
rownames(gene.panel.expr) = gene.panel.expr$GeneName
gene.panel.expr$GeneID = NULL
gene.panel.expr$GeneName = NULL
gene.panel.expr$GeneBiotype = NULL
#gene.panel.expr$GeneAnnotationtLevel=NULL
#gene.panel.expr$GeneFullName = NULL 

# deconvolution using extracted gene marker expression #
deconv.gene.panel.expr = deconvolute.alg(gene.panel.expr, reference = ml.model)

# save deconvolution results in csv files #
plot.data = NULL
for (alg in c("qprog", "cs-qprog", "DSA", "ssFrobenius", "meanProfile", "deconf")){
  #for (alg in c("ssFrobenius", "meanProfile")){
  result = data.frame(deconv.gene.panel.expr[alg], check.names = F)
  rownames(result) = c("Astrocyte","Microglia","Neuron","Oligodendrocyte")
  result = data.frame(t(result))
  rownames(result) = gsub("cs.qprog","cs-qprog",rownames(result))
  #result$Algorithm = sapply(strsplit(rownames(result),"[.]"), `[`, 1)
  result$Algorithm = str_split_fixed(rownames(result), "[.]", 2)[,1] 
  #result$Subject = sapply(sapply(strsplit(rownames(result),"[.]"), `[`, -1), paste, collapse = ".")
  result$Subject = str_split_fixed(rownames(result), "[.]", 2)[,2] 
  result = result[,c("Subject","Algorithm","Astrocyte","Microglia","Neuron","Oligodendrocyte")]
  rownames(result) = NULL
  if (is.null(plot.data)){
    plot.data = result
  } else {
    plot.data = rbind(plot.data,result)
  }
  if (r == 'rosmap') {
    write.table(result, paste("Public_Data/bulkRNASeq/201506_ROSMAP/dorsolateral_prefrontal_cortex/05.-Analyses/deconvolution_analysis/deconvolution_",alg,"_results.tsv",sep =""), quote = F, row.names = F, sep="\t")
  } else {
    write.table(result, paste("Public_Data/bulkRNASeq/201812_MSBB/",r,"/05.-Analyses/deconvolution_analysis/deconvolution_",alg,"_results.tsv",sep =""), quote = F, row.names = F, sep="\t")
  }
}

# plot deconvolution results #
plot.data.m = reshape2::melt(plot.data,id.vars = c("Subject","Algorithm"))
colnames(plot.data.m)[colnames(plot.data.m) == "variable"] = "Celltype"
colnames(plot.data.m)[colnames(plot.data.m) == "value"] = "Proportion"

if (r == 'rosmap') {
  pdf("Public_Data/bulkRNASeq/201506_ROSMAP/dorsolateral_prefrontal_cortex/05.-Analyses/deconvolution_analysis/average_cell_proportion.pdf", width = 10, height = 6)
} else {
  pdf(paste0("Public_Data/bulkRNASeq/201812_MSBB/",r,"/05.-Analyses/deconvolution_analysis/average_cell_proportion.pdf"), width = 10, height = 6)
}
ggplot(plot.data.m ,aes(x = Celltype, y = Proportion, fill = Celltype)) + 
  geom_bar(stat="summary",fun.y="mean",position = "dodge",width = 0.7) + facet_wrap( ~ Algorithm) +
  geom_hline(yintercept=0) + scale_fill_manual(values = c("#F8766D","#7CAE00","#00BFC4","#C77CFF")) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width=0.7), width = 0.3, size=0.5) + 
  stat_summary(aes(label=abs(round(..y..,3))),  position = position_dodge(width=0.7), fun.y=mean, geom="text", size=4,vjust = -0.9) + 
  xlab("") + theme_bw() + ggtitle('MSBB - BM36') + 
  theme(legend.position="top",strip.text.x = element_text(colour = "black", face = "bold"), 
        plot.title = element_text(size = 16, hjust=0.5, face="bold"), 
        axis.text.x = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="black", fill=NA), 
        panel.border = element_rect(colour = "black")) + 
  guides(fill=guide_legend(title=""))

dev.off()

## plot per cell type
if (r == 'rosmap') {
  pdf("Public_Data/bulkRNASeq/201506_ROSMAP/dorsolateral_prefrontal_cortex/05.-Analyses/deconvolution_analysis/deconvolution_per_sample.pdf", width = 10, height = 8) 
} else {
  pdf(paste0("Public_Data/bulkRNASeq/201812_MSBB/",r,"/05.-Analyses/deconvolution_analysis/deconvolution_per_sample.pdf"), width = 10, height = 8)
}
ggplot(plot.data.m ,aes(x = Celltype, y = Proportion)) + 
  geom_jitter(size=2, width=0.3, aes(color=factor(Celltype))) + theme_bw(base_size = 12) + facet_wrap( ~ Algorithm) +
  #geom_boxplot(alpha = 0.01, lwd=0.03, position=position_dodge()) +
  labs(x='') + ggtitle('MSBB - BM36') + 
  theme(axis.text.x=element_text(size=14, angle =90, vjust=0.5, face="bold"), axis.title.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
        plot.title = element_text(size = 16, hjust=0.5, face="bold"),
        legend.position="none",
        axis.ticks=element_blank(), panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.background = element_rect(color = "black"),
        panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 
#  scale_color_manual(name="", values=c("p"="#df8640","n"="gray40"), labels=c("p"="Tumor","n"="Normal"),
dev.off()

## plot meanprofile & ssNMF
if (r == 'rosmap') {
  pdf("Public_Data/bulkRNASeq/201506_ROSMAP/dorsolateral_prefrontal_cortex/05.-Analyses/deconvolution_analysis/deconvolution_meanP_ssNMF.pdf", width = 10, height = 8) 
} else {
  pdf(paste0("Public_Data/bulkRNASeq/201812_MSBB/",r,"/05.-Analyses/deconvolution_analysis/deconvolution_meanP_ssNMF.pdf"), width = 10, height = 8)
}
ggplot(plot.data.m[plot.data.m$Algorithm %in% c('meanProfile', 'ssFrobenius'),] ,aes(x = Celltype, y = Proportion)) + 
  geom_jitter(size=2, width=0.3, aes(color=factor(Celltype))) + theme_bw(base_size = 12) + facet_wrap( ~ Algorithm, nrow=2) +
  #geom_text(label= ifelse (sunshine.all[sunshine.all$Algorithm %in% c('meanProfile', 'ssFrobenius'),'Cell_type']=="iPSC-Neurons", 
  #                         sunshine.all[sunshine.all$Algorithm %in% c('meanProfile', 'ssFrobenius'),'Subject'],''), size=2, angle=45, vjust=1) + 
  labs(x='', title='MSBB - BM36') + guides(colour = guide_legend(override.aes = list(size = 6, stroke = 2))) + 
  theme(axis.text.x=element_text(size=14, angle =90, vjust=0.5, face="bold"), axis.title.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
        plot.title = element_text(size = 16, hjust=0.5, face="bold"),
        legend.position="top",
        axis.ticks=element_blank(), panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.background = element_rect(color = "black"),
        panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 
dev.off()




