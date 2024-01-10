setwd("/dragon-home/general/Cruchaga_Data/bulkRNASeq")
library(CellMix)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)

#study = "201703_MendelianVsSporadics/pipeline-v1-hg38"
study = "201904_MGI_bulkRNAseq_SUNSHINE"
study="201703_MendelianVsSporadics"

# load in deconvolution algorithm #
#source( 'deconvolute_alg.R')
deconvolute.alg = function ( expression, reference, forceExtractMarkers=F, 
                             selected.algs = c( "qprog", "cs-qprog", "DSA", "ssFrobenius", "meanProfile", "deconf"), 
                             types.reference = NA){
  
  #crashes again: , "ssKL"
  
  
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
# rowname has to be gene name #
gene <- read.table(paste0(study, '/02.-ProcessedData/Hg38/03.-STAR/combined_count_matrix/gene_fpkm_matrix_2nd_read_strand.tsv'), header =T, check.names = F, stringsAsFactors = F)

# load in deconvolution learning model #
#load( file="Deconvolution_ml_model.RData")

# load in gene marker panel #
gene.panel = read.csv('/home/liz/RNAseq/Decon_paper/Synapse_data/website/Brain_panel_genes_drop_PLP1_CCL2.csv', stringsAsFactors = F, row.names = 1)
gene.panel = gene.panel[, c('cellType', 'Human')]
gene.panel = subset( gene.panel, !is.na(gene.panel$Human) )
# read Celeste list of genes 
gene.panel2 = read.table('201904_MGI_bulkRNAseq_SUNSHINE/05.-Analyses/deconvolution_analysis/human_microglia_markers.tsv', header =T, stringsAsFactors = F)
gene.panel = unique(rbind(gene.panel, gene.panel2))
colnames(gene.panel) <- c('cellType', 'gene')
## leave-one-out (e.g. GFAP)
gene.panel <- gene.panel[!gene.panel$gene %in% c("GFAP", "AQP4"), ]

# extract gene marker expression from the samples  #
gene.panel.expr = gene[gene$GeneName %in% gene.panel$gene,]
rownames(gene.panel.expr) = gene.panel.expr$GeneName
gene.panel.expr$GeneID = NULL
gene.panel.expr$GeneName = NULL

gene.panel.expr$GeneBiotype = NULL
gene.panel.expr$GeneAnnotationtLevel=NULL
gene.panel.expr$GeneFullName = NULL 

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
  
  write.table(result, paste(study,"/03.-Analyses/deconvolution_analysis/deconvolution_",alg,"_results.tsv",sep =""), quote = F, row.names = F, sep="\t")
  
}

# plot deconvolution results #
plot.data.m = melt(plot.data,id.vars = c("Subject","Algorithm"))
colnames(plot.data.m)[colnames(plot.data.m) == "variable"] = "Celltype"
colnames(plot.data.m)[colnames(plot.data.m) == "value"] = "Proportion"


# ggplot(plot.data.m, aes(x = Subject, y = Proportion, color = Celltype)) + geom_point() +
#   facet_wrap( ~ Algorithm) + theme_bw() + 
#   theme(legend.position="top",strip.text.x = element_text(colour = "black", face = "bold"), 
#         axis.ticks.x=element_blank(),axis.text.x=element_blank(),
#         panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#         strip.background = element_rect(colour="black", fill=NA), panel.border = element_rect(colour = "black")) + 
#   guides(fill=guide_legend(title="Cell Type"))

pdf(paste0(study,"/05.-Analyses/deconvolution_analysis/deconvolution_demo.pdf"), width = 12, height = 8)
    
ggplot(plot.data.m ,aes(x = Celltype, y = Proportion, fill = Celltype)) + 
  geom_bar(stat="summary",fun.y="mean",position = "dodge",width = 0.7) + facet_wrap( ~ Algorithm) +
  geom_hline(yintercept=0) + scale_fill_manual(values = c("#F8766D","#7CAE00","#00BFC4","#C77CFF")) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width=0.7), width = 0.3, size=0.5) + 
  stat_summary(aes(label=abs(round(..y..,3))),  position = position_dodge(width=0.7), fun.y=mean, geom="text", size=4,vjust = -0.9) + 
  xlab("") + theme_bw() + ggtitle('')
  theme(legend.position="top",strip.text.x = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(face = "bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="black", fill=NA), panel.border = element_rect(colour = "black")) + 
  guides(fill=guide_legend(title="Cell Type"))

dev.off()


############################################################################################################
#### compute the correlation between pipeline vs Zeran's for Mendelian 
zeran.res = read.csv('/40/Cruchaga_Data/bulkRNASeq/201703_MendelianVsSporadics/03.-Phenotype/archive/2019_03_19_WashU_MendelianVsSporadics_clinical_v3.csv' , header =T, sep=",", stringsAsFactors = F)
zeran.res[is.na(zeran.res$brain_region) & grepl('Pooled_RNA', zeran.res$ID_RNAseq), 'brain_region'] = "Pooled_RNA"

Celeste.Neuron = zeran.res[grepl('_Neuron_', zeran.res$ID_RNAseq), 'ID_RNAseq']
Celeste.Astro = zeran.res[grepl('_Astro_[0-9]', zeran.res$ID_RNAseq), 'ID_RNAseq']
Celeste.Astro_IGF = zeran.res[grepl('_Astro_IGF_', zeran.res$ID_RNAseq), 'ID_RNAseq']


## extract cell lines
cell.lines = unique(gene.panel$cellType)
cell.lines = cell.lines[cell.lines %in% colnames(zeran.res)]

### read pipeline results for the ssNMF + meanProfile 
algorithm =c('ssFrobenius', 'meanProfile')
pdf(paste0(study,"/05.-Analyses/deconvolution_analysis/deconvolution_corr.pdf"), width = 12, height = 8)

for (alg in algorithm) { 
  
  pipeline.res = read.table(paste0(study,"/05.-Analyses/deconvolution_analysis/deconvolution_",alg,"_results.tsv"), header =T, check.names = F, stringsAsFactors = F)
  pipeline.res$Subject = gsub("H_VY.", "H_VY-", pipeline.res$Subject)
  
  cl.plots = list()
  for (i in 1:length(cell.lines)) {
    cl = cell.lines[i]
    zz = zeran.res[, c('ID_RNAseq','brain_region', cl)]
    colnames(zz) = c('Subject', 'brain_region', 'Z' )
    pp = pipeline.res[, c('Subject', cl)]
    colnames(pp) = c('Subject', 'P')
    ## merge 
    cl.res = merge(zz, pp)
    cl.cor.p = signif(cor(cl.res$Z, cl.res$P, method ="pearson"), digits = 3)
    cl.cor.s = signif(cor(cl.res$Z, cl.res$P, method ="spearman"), digits = 3)
    
    ## plot the correlation 
    #cl.res$comp = 'N'
    #cl.res[cl.res$P >= cl.res$Z, 'comp'] = 'Y'
    ## add cell lines information 
    cl.res$region = 'Parietal'
    cl.res[cl.res$brain_region=="insular", 'region'] = "Insular"
    cl.res[cl.res$brain_region=="TBD/possibly_insular", 'region'] = "Possibly_insular"
    cl.res[cl.res$brain_region=="skin;_lower_back_dermal_fibroblasts" & cl.res$Subject %in% Celeste.Neuron, 'region'] = "Neuron"
    cl.res[cl.res$brain_region=="skin;_lower_back_dermal_fibroblasts" & cl.res$Subject %in% Celeste.Astro, 'region'] = "Astrocyte"
    cl.res[cl.res$brain_region=="skin;_lower_back_dermal_fibroblasts" & cl.res$Subject %in% Celeste.Astro_IGF, 'region'] = "Astrocyte_IGF"
    cl.res[cl.res$brain_region=="Pooled_RNA", 'region'] = "Pooled_RNA"
    
    p = ggplot(cl.res, aes(P, Z, color= factor(region))) + geom_point(shape=19, size=1.5) + theme_bw()
    p = p + labs(x='Pipepline', y = 'Zeran\'s')
    p = p + ggtitle(paste(cl, '\nPearson = ', cl.cor.p,', Spearman = ', cl.cor.s) )
    p = p + theme(axis.text.x=element_text(size=12, color="black"),
                  axis.text.y=element_text(size=12, color="black"), legend.title=element_blank(),
                  axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
                  plot.title = element_text(size = 12, hjust=0.5, color="black", face="bold"),
                  legend.position="right", legend.key.size = unit(5,'mm'), legend.text=element_text(size=10), 
                  panel.border = element_rect(linetype='solid', color='black'))
    p = p + scale_color_brewer(name="", palette = "Dark2", direction = -1 )
    #p = p + scale_color_manual(name ="PipelineProp > Zeran", values =c("N"="#F8766D", "Y"="#00BFC4"), guide =  guide_legend(override.aes = list(size = 5))) 
    
    cl.plots[[i]] <- p
    
  }
  
  mytitle = textGrob(paste0('Cell proportion for Mendelian brains using ', alg, ' algorithm\n'), gp=gpar(fontsize=20, fontface="bold"))
  grid.arrange(grobs=cl.plots, ncol=2, nrwo=2, top= mytitle)  
}

dev.off()


###### plot all samples as boxplot for Menedlian 
## add the region group 
plot.data.m$Subject = gsub("H_VY.", "H_VY-", plot.data.m$Subject)

plot.data.m$group = 'Brain'
plot.data.m[plot.data.m$Subject %in% Celeste.Neuron, 'group'] = 'Neuron'
plot.data.m[plot.data.m$Subject %in% Celeste.Astro, 'group'] = 'Astro'
plot.data.m[plot.data.m$Subject %in% Celeste.Astro_IGF, 'group'] = 'Astro_IGF'

pdf(paste0(study,"/05.-Analyses/deconvolution_analysis/deconvolution_per_sample.pdf"), width = 12, height = 8)

ggplot(plot.data.m ,aes(x = Celltype, y = Proportion)) + 
  geom_jitter(size=1.5, shape=1, width=0.3, aes(color=factor(group))) + theme_bw(base_size = 12) + facet_wrap( ~ Algorithm) +
  #geom_boxplot(alpha = 0.01, lwd=0.03, position=position_dodge()) +
  labs(x='') + 
  theme(axis.text.x=element_text(size=14, angle =90, vjust=0.5, face="bold"), axis.title.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
        plot.title = element_text(size = 20, hjust=0.5, face="bold"),
        legend.position="right",
        axis.ticks=element_blank(), panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.background = element_rect(color = "black"),
        panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 
#  scale_color_manual(name="", values=c("p"="#df8640","n"="gray40"), labels=c("p"="Tumor","n"="Normal"),
dev.off()

#########################################################################
############################# SUNSHINE ##################################
#########################################################################

##### plot all samples as boxplot for Sunshine 
clin = read.csv(paste0(study, '/03.-Phenotype/2019_09_25_WashU_SUNSHINE_clinical.csv') , header =T, sep=",", stringsAsFactors = F)
ss = read.csv(paste0(study, '/01.-RawData/fq1_fq2_fullsm.csv'), header =T, sep=",", stringsAsFactors = F)
ss = ss[grepl('HW', ss$fq1), ]
ss$Subj_ID = str_split_fixed(ss$fullsm, "\\^", 2)[,1]
ss$Barcode = str_split_fixed(ss$fullsm, "\\^", 2)[,2]
ss$Barcode = gsub('\\^SUNSHINE', '', ss$Barcode)
ss$Sample_ID = str_split_fixed(ss$fq1, "_R1", 2)[,1]

Piccio.ids = ss[grepl('Piccio', ss$Subj_ID), 'Sample_ID']
brain.ids = ss[grepl('MAP_', ss$Subj_ID), 'Sample_ID']
Celeste.ids = ss$Sample_ID[!ss$Sample_ID %in% c(Piccio.ids, brain.ids)]

## add groups 
plot.data.m$study = 'WASHU'
plot.data.m[plot.data.m$Subject %in% Celeste.ids, 'study'] = 'Celeste'
plot.data.m[plot.data.m$Subject %in% Piccio.ids, 'study'] = 'Piccio'

pdf(paste0(study,"/05.-Analyses/deconvolution_analysis/deconvolution_per_sample.pdf"), width = 12, height = 8)

ggplot(plot.data.m ,aes(x = Celltype, y = Proportion)) + 
  geom_jitter(size=1.5, shape=1, width=0.3, aes(color=factor(study))) + theme_bw(base_size = 12) + facet_wrap( ~ Algorithm) +
  #geom_boxplot(alpha = 0.01, lwd=0.03, position=position_dodge()) +
  labs(x='') + 
  theme(axis.text.x=element_text(size=14, angle =90, vjust=0.5, face="bold"), axis.title.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
        plot.title = element_text(size = 20, hjust=0.5, face="bold"),
        legend.position="right",
        axis.ticks=element_blank(), panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.background = element_rect(color = "black"),
        panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 
#  scale_color_manual(name="", values=c("p"="#df8640","n"="gray40"), labels=c("p"="Tumor","n"="Normal"),
dev.off()

### check cell lines (extract cell lines)
cell.data = read.table(paste0(study, '/03.-Phenotype/sunshine_celllines_info.tsv') , header =T, sep="\t", stringsAsFactors = F)
#cell.lines.info = data.frame(Subject = ss[grepl("F0510", ss$Subj_ID), 'Sample_ID'], Cell_type = "iPSC-Neurons", Barcode = NA, Subj_ID =cell.line.id, stringsAsFactors = F)

cell.lines.info = NULL
for (k in 1:nrow(cell.data)) {
  print (k)
  cell.line.id = cell.data$SUNSHINE_Individual_Name[k]
  #if (cell.line.id=="F0510") { next }
  cell.barcdoe = as.numeric(gsub("^.*_", "", cell.data$SUNSHINE_Sample_Name[k]))
  cell.type = cell.data$Cell.type[k]
  
  # extract sample name 
  cell.sample = ss[grepl(cell.barcdoe, ss$fullsm), 'Sample_ID']
  
  ### Fix swapped cell ines for SUNSHIN 
  if (study == "201904_MGI_bulkRNAseq_SUNSHINE") {
    if (cell.barcdoe == 1195060761) {
      cell.type = "iPSC-Neurons"
      cell.sample = ss[grepl(cell.barcdoe, ss$fullsm), 'Sample_ID']
      cell.barcdoe = 1195060476
    } else if (cell.barcdoe == 1195060476) {
      cell.type = "iMGL"
      cell.sample = ss[grepl(cell.barcdoe, ss$fullsm), 'Sample_ID']
      cell.barcdoe = 1195060761
    }
  }
 
  
  #cell.line.gt = cell.data$Genotype[k]
  
  
  if (length(cell.sample) !=0) {
    cell.lines.res = data.frame(Subject = cell.sample , Cell_type = cell.type, Barcode = cell.barcdoe, Subj_ID =cell.line.id, stringsAsFactors = F)  
  }
  
  cell.lines.info = rbind(cell.lines.info, cell.lines.res)
}

### merge with previous results 
sunshine.all = merge(plot.data.m, cell.lines.info, all.x = T, sort =F)
sunshine.all[is.na(sunshine.all$Cell_type), 'Cell_type'] = 'Brain'

### Fix swapped cell lines 

pdf(paste0(study,"/05.-Analyses/deconvolution_analysis/deconvolution_celllines_dist.pdf"), width = 12, height = 8)
ggplot(sunshine.all ,aes(x = Celltype, y = Proportion)) + 
  geom_jitter(size=1.5, shape=1, width=0.3, aes(color=factor(Cell_type))) + theme_bw(base_size = 12) + facet_wrap( ~ Algorithm) +
  geom_boxplot(alpha = 0.5, lwd=0.03, position=position_dodge()) +
  labs(x='') + guides(colour = guide_legend(override.aes = list(size = 6, stroke = 2))) + 
  theme(axis.text.x=element_text(size=14, angle =90, vjust=0.5, face="bold"), axis.title.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
        plot.title = element_text(size = 20, hjust=0.5, face="bold"),
        legend.position="top",
        axis.ticks=element_blank(), panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.background = element_rect(color = "black"),
        panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 

dev.off()

### plot meanProfile & ssFrobenius
pdf(paste0(study,"/05.-Analyses/deconvolution_analysis/Figures/deconvolution_celllines_meanP_ssNMF.pdf"), width = 12, height = 8)
ggplot(sunshine.all[sunshine.all$Algorithm %in% c('meanProfile', 'ssFrobenius'),] ,aes(x = Celltype, y = Proportion)) + 
  geom_jitter(size=1.5, shape=1, width=0.3, aes(color=factor(Cell_type))) + theme_bw(base_size = 12) + facet_wrap( ~ Algorithm, nrow=2) +
  geom_text(label= ifelse (sunshine.all[sunshine.all$Algorithm %in% c('meanProfile', 'ssFrobenius'),'Cell_type']=="iPSC-Neurons", 
                           sunshine.all[sunshine.all$Algorithm %in% c('meanProfile', 'ssFrobenius'),'Subject'],''), size=2, angle=45, vjust=1) + 
  labs(x='') + guides(colour = guide_legend(override.aes = list(size = 6, stroke = 2))) + 
  theme(axis.text.x=element_text(size=14, angle =90, vjust=0.5, face="bold"), axis.title.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
        plot.title = element_text(size = 20, hjust=0.5, face="bold"),
        legend.position="top",
        axis.ticks=element_blank(), panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.background = element_rect(color = "black"),
        panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 

dev.off()

### plot per cell 
pdf(paste0(study,"/05.-Analyses/deconvolution_analysis/Figures/deconvolution_celllines_per_cell_meanP_fixed.pdf"), width = 12, height = 8)
ggplot(sunshine.all[sunshine.all$Algorithm =='meanProfile',] ,aes(x = Subject, y = Proportion)) + 
  geom_jitter(size=1.5, shape=1, width=0.5, aes(color=factor(Cell_type))) + theme_bw(base_size = 12) + facet_wrap( ~ Celltype, nrow=2) +
  ggtitle('Cell Proportion using meanProfile Algorithm') +
  labs(x='Sample') + guides(colour = guide_legend(override.aes = list(size = 6, stroke = 2))) + 
  #geom_text(label=sunshine.all[sunshine.all$Algorithm =='meanProfile', 'Subject']) +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
        plot.title = element_text(size = 20, hjust=0.5, face="bold"),
        legend.position="top", strip.text = element_text(face="bold"), 
        axis.ticks=element_blank(), panel.grid.major = element_line(colour = "white", size = 0.3), 
        legend.title=element_blank(), legend.background = element_rect(color = "black"),
        panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 

dev.off()

pdf(paste0(study,"/05.-Analyses/deconvolution_analysis/Figures/deconvolution_celllines_per_cell_ssNMF_fixed.pdf"), width = 12, height = 8)
ggplot(sunshine.all[sunshine.all$Algorithm =='ssFrobenius',] ,aes(x = Subject, y = Proportion)) + 
  geom_jitter(size=1.5, shape=1, width=0.3, aes(color=factor(Cell_type))) + theme_bw(base_size = 12) + facet_wrap( ~ Celltype, nrow=2) +
  ggtitle('Cell Proportion using ssNMF Algorithm') +
  labs(x='Sample') + guides(colour = guide_legend(override.aes = list(size = 6, stroke = 2))) + 
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
        plot.title = element_text(size = 20, hjust=0.5, face="bold"),
        legend.position="top", panel.grid.major = element_line(colour = "white", size = 0.3), 
        axis.ticks=element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_text(face="bold"), 
        legend.title=element_blank(), legend.background = element_rect(color = "black"),
        panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 

dev.off()




##############################################################################################################
#### Deconvultion analysis for Celeste 
##############################################################################################################
setwd("/storage1/fs1/karchc/Active/fenix/40/Karch_Data/202009_MGI_bulkRNAseq_/ngi-pipeline")
study="keBYB3d4"
study="zdw3Pbdw"
dir.create(paste(study,"/03.-Analyses/deconvolution_analysis",sep =""))
dir.create(paste(study,"/03.-Analyses/deconvolution_analysis/Results",sep =""))
dir.create(paste(study,"/03.-Analyses/deconvolution_analysis/Figures",sep =""))

# load in sample gene expression file #
gene <- read.table(paste0(study, '/02.-ProcessedData/04.-Salmon/Hg38/combined_exp_matricies/gene_quant_matrix_NumReads.tsv'), header =T, check.names = F, stringsAsFactors = F)

# load in gene marker panel #
gene.panel = read.csv('/home/liz/RNAseq/Decon_paper/Synapse_data/website/Brain_panel_genes_drop_PLP1_CCL2.csv', stringsAsFactors = F, row.names = 1)
gene.panel = gene.panel[, c('cellType', 'Human')]
gene.panel = subset( gene.panel, !is.na(gene.panel$Human) )

# read Celeste list of genes 
gene.panel2 = read.table('/40/Cruchaga_Data/bulkRNASeq/201904_MGI_bulkRNAseq_SUNSHINE/05.-Analyses/deconvolution_analysis/human_microglia_markers.tsv', header =T, stringsAsFactors = F)
gene.panel = unique(rbind(gene.panel, gene.panel2))
colnames(gene.panel) <- c('cellType', 'gene')
## leave-one-out (e.g. GFAP)
#gene.panel <- gene.panel[!gene.panel$gene %in% c("GFAP", "AQP4"), ]

# extract gene marker expression from the samples  #
gene.panel.expr = gene[gene$GeneName %in% gene.panel$gene,]
rownames(gene.panel.expr) = gene.panel.expr$GeneName
gene.panel.expr$GeneID = NULL
gene.panel.expr$GeneName = NULL

gene.panel.expr$GeneBiotype = NULL
gene.panel.expr$GeneAnnotationtLevel=NULL
gene.panel.expr$GeneFullName = NULL 

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
  
  write.table(result, paste(study,"/03.-Analyses/deconvolution_analysis/Results/deconvolution_",alg,"_results.tsv",sep =""), quote = F, row.names = F, sep="\t")
  
}


# plot deconvolution results #
plot.data.m = melt(plot.data,id.vars = c("Subject","Algorithm"))
colnames(plot.data.m)[colnames(plot.data.m) == "variable"] = "Celltype"
colnames(plot.data.m)[colnames(plot.data.m) == "value"] = "Proportion"


# ggplot(plot.data.m, aes(x = Subject, y = Proportion, color = Celltype)) + geom_point() +
#   facet_wrap( ~ Algorithm) + theme_bw() + 
#   theme(legend.position="top",strip.text.x = element_text(colour = "black", face = "bold"), 
#         axis.ticks.x=element_blank(),axis.text.x=element_blank(),
#         panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#         strip.background = element_rect(colour="black", fill=NA), panel.border = element_rect(colour = "black")) + 
#   guides(fill=guide_legend(title="Cell Type"))

pdf(paste0(study,"/03.-Analyses/deconvolution_analysis/Figures/deconvolution_demo.pdf"), width = 12, height = 8)
ggplot(plot.data.m ,aes(x = Celltype, y = Proportion, fill = Celltype)) + 
  geom_bar(stat="summary",fun.y="mean",position = "dodge",width = 0.7, show.legend = F) + facet_wrap( ~ Algorithm) +
  geom_hline(yintercept=0) + scale_fill_manual(values = c("#F8766D","#7CAE00","#00BFC4","#C77CFF")) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width=0.7), width = 0.3, size=0.5) + 
  stat_summary(aes(label=abs(round(..y..,3))),  position = position_dodge(width=0.7), fun.y=mean, geom="text", size=4,vjust = -0.9) + 
  xlab("") + theme_bw() + ggtitle('')
theme(legend.position="top",strip.text.x = element_text(colour = "black", face = "bold"), 
      axis.text.x = element_text(face = "bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill=NA), panel.border = element_rect(colour = "black")) + 
  guides(fill=guide_legend(title="Cell Type"))
dev.off()

## plot per cell type  
pdf(paste0(study,"/03.-Analyses/deconvolution_analysis/Figures/deconvolution_per_sample.pdf"), width = 12, height = 8)
ggplot(plot.data.m ,aes(x = Celltype, y = Proportion)) + 
  geom_jitter(size=1.5, shape=1, width=0.3, aes(color=factor(study))) + theme_bw(base_size = 12) + facet_wrap( ~ Algorithm) +
  #geom_boxplot(alpha = 0.01, lwd=0.03, position=position_dodge()) +
  labs(x='') + 
  theme(axis.text.x=element_text(size=14, angle =90, vjust=0.5, face="bold"), axis.title.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
        plot.title = element_text(size = 20, hjust=0.5, face="bold"),
        legend.position="right",
        axis.ticks=element_blank(), panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.background = element_rect(color = "black"),
        panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 
#  scale_color_manual(name="", values=c("p"="#df8640","n"="gray40"), labels=c("p"="Tumor","n"="Normal"),
dev.off()

## plot meanprofile & ssNMF
pdf(paste0(study,"/03.-Analyses/deconvolution_analysis/Figures/deconvolution_meanP_ssNMF.pdf"), width = 12, height = 8)
ggplot(plot.data.m[plot.data.m$Algorithm %in% c('meanProfile', 'ssFrobenius'),] ,aes(x = Celltype, y = Proportion)) + 
  geom_jitter(size=1.5, shape=1, width=0.3, aes(color=factor(Celltype))) + theme_bw(base_size = 12) + facet_wrap( ~ Algorithm, nrow=2) +
  #geom_text(label= ifelse (sunshine.all[sunshine.all$Algorithm %in% c('meanProfile', 'ssFrobenius'),'Cell_type']=="iPSC-Neurons", 
  #                         sunshine.all[sunshine.all$Algorithm %in% c('meanProfile', 'ssFrobenius'),'Subject'],''), size=2, angle=45, vjust=1) + 
  labs(x='') + guides(colour = guide_legend(override.aes = list(size = 6, stroke = 2))) + 
  theme(axis.text.x=element_text(size=14, angle =90, vjust=0.5, face="bold"), axis.title.x=element_text(size=18, face="bold"),
        axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
        plot.title = element_text(size = 20, hjust=0.5, face="bold"),
        legend.position="top",
        axis.ticks=element_blank(), panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.background = element_rect(color = "black"),
        panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 

dev.off()

# ### plot per cell 
# pdf(paste0(study,"/05.-Analyses/deconvolution_analysis/Figures/deconvolution_celllines_per_cell_meanP_fixed.pdf"), width = 12, height = 8)
# ggplot(plot.data.m[plot.data.m$Algorithm =='meanProfile',] ,aes(x = Subject, y = Proportion)) + 
#   geom_jitter(size=1.5, shape=1, width=0.5, aes(color=factor(Cell_type))) + theme_bw(base_size = 12) + facet_wrap( ~ Celltype, nrow=2) +
#   ggtitle('Cell Proportion using meanProfile Algorithm') +
#   labs(x='Sample') + guides(colour = guide_legend(override.aes = list(size = 6, stroke = 2))) + 
#   #geom_text(label=sunshine.all[sunshine.all$Algorithm =='meanProfile', 'Subject']) +
#   theme(axis.text.x=element_blank(), axis.title.x=element_text(size=18, face="bold"),
#         axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
#         plot.title = element_text(size = 20, hjust=0.5, face="bold"),
#         legend.position="top", strip.text = element_text(face="bold"), 
#         axis.ticks=element_blank(), panel.grid.major = element_line(colour = "white", size = 0.3), 
#         legend.title=element_blank(), legend.background = element_rect(color = "black"),
#         panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 
# 
# dev.off()

# pdf(paste0(study,"/03.-Analyses/deconvolution_analysis/Figures/deconvolution_celllines_per_cell_ssNMF_fixed.pdf"), width = 12, height = 8)
# ggplot(sunshine.all[sunshine.all$Algorithm =='ssFrobenius',] ,aes(x = Subject, y = Proportion)) + 
#   geom_jitter(size=1.5, shape=1, width=0.3, aes(color=factor(Cell_type))) + theme_bw(base_size = 12) + facet_wrap( ~ Celltype, nrow=2) +
#   ggtitle('Cell Proportion using ssNMF Algorithm') +
#   labs(x='Sample') + guides(colour = guide_legend(override.aes = list(size = 6, stroke = 2))) + 
#   theme(axis.text.x=element_blank(), axis.title.x=element_text(size=18, face="bold"),
#         axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=18, face="bold"),
#         plot.title = element_text(size = 20, hjust=0.5, face="bold"),
#         legend.position="top", panel.grid.major = element_line(colour = "white", size = 0.3), 
#         axis.ticks=element_blank(), panel.grid.minor = element_blank(),
#         strip.text = element_text(face="bold"), 
#         legend.title=element_blank(), legend.background = element_rect(color = "black"),
#         panel.background=element_blank(),legend.text=element_text(size=16, face="bold")) 
# 
# dev.off()



