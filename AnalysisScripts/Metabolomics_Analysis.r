setwd("/home/eteleeb/projects")
library(ggplot2)
library(gridExtra)
library(reshape2)
library(data.table)
library(gplots)
library(ggsignif)
library(RColorBrewer)
library(VennDiagram)

###########################################################################
## extract shared metaboloties between Knight ADRC and ROSMAP  
###########################################################################
rosmap.clin <- unique(read.csv('/40/Public_Data/bulkRNASeq/201506_ROSMAP/Gene_Expression/03.-Phenotype/new/2020_12_01_ROSMAP_clinical.csv', sep=",", header =T, stringsAsFactors = F))

## read ROSMAP results for plotting and merge iwht clusters 
#rosmap.clin <- unique(read.csv('/40/Public_Data/bulkRNASeq/201506_ROSMAP/Gene_Expression/03.-Phenotype/new/2020_12_01_ROSMAP_clinical.csv', sep=",", header =T, stringsAsFactors = F))
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

# ## handle missing values 
# rosmap.metab$na_count <- apply(rosmap.metab[, metab.cols], 1, function(x) sum(is.na(x)))
# rosmap.metab <- rosmap.metab[rosmap.metab$na_count < (length(metab.cols) * 0.20), ]
# rosmap.metab$na_count <- NULL
# ## replace all NA with the lowest value
# for(mm in metab.cols){
#   indx = which(is.na(rosmap.metab[, mm]))
#   r.min <- mean(rosmap.metab[, mm], na.rm = T)
#   rosmap.metab[indx, mm] <- r.min
# }

## ---------------------------------------------------------------------------------------------------
## extract intersection between Knight ADRC and ROSMAP
## ---------------------------------------------------------------------------------------------------
## read might ADRC and rosmap data 
knight_adrc_meta <- readRDS('omics_integration/data/Metabolomics/04-metab_meta.rds')
rosmap.meta <- readRDS('omics_integration/data/ROSMAP/01-metab_meta_filtered_cerad_recover2.rds')

#data.type <- 'C4vsall_C1vsC2'  
data.type <- 'C4vsCO_C1vsCO'  
grp1 <- unlist(strsplit(data.type, '_'))[1]
grp2 <- unlist(strsplit(data.type, '_'))[2]

## read knight ADRC results
if (data.type =='C4vsall_C1vsC2') {
  knight_adrc <- read.csv('omics_integration/data/Metabolomics/DE_results_metabolomics/03-metab_effect_pval_4vsall.csv', header =T, stringsAsFactors = F, sep=",")
  rosmap <- read.table('omics_integration/Replication/ROSMAP/Metab_DE/1vs2/1vs2.res.all.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t', row.names = 1, quote = "")
  dt <- 'all'
} else if (data.type=='C4vsCO_C1vsCO') {
  knight_adrc <- read.csv('omics_integration/data/Metabolomics/DE_results_metabolomics/03-metab_effect_pval_4vsCO.csv', header =T, stringsAsFactors = F, sep=",")
  rosmap <- read.table('omics_integration/Replication/ROSMAP/Metab_DE/1vs0/1vs0.res.all.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t', row.names = 1, quote = "")
  dt <- 'CO'
}

colnames(knight_adrc)[colnames(knight_adrc)=="metab"] <- "BIOCHEMICAL"
knight_adrc <- merge(knight_adrc, knight_adrc_meta[, c('BIOCHEMICAL', 'CHEMICAL.ID')])

## read AD vs CO results 
KnightADRC_ADvsCO <- read.csv('omics_integration/data/Metabolomics/DE_results_metabolomics/03-metab_effect_pval_ADvsCO.csv', header =T, stringsAsFactors = F, sep=",")

## read ROSMAP results 
rosmap$SHORT_NAME <- rownames(rosmap)
rosmap <- merge(rosmap, rosmap.meta[, c('SHORT_NAME', 'CHEM_ID')])

## extract all features for fisher's test 
all.metab <- union(knight_adrc$BIOCHEMICAL, rosmap$SHORT_NAME)

## extract significant results 
knight_adrc <- knight_adrc[knight_adrc$pval < 0.05,]
KnightADRC_ADvsCO <- KnightADRC_ADvsCO[KnightADRC_ADvsCO$pval < 0.05, ]
rosmap <- rosmap[rosmap$p.value < 0.05, ]

## read ADADvsCO results 
adad.res <- read.csv('omics_integration/data/Metabolomics/ADADvsCO_KnightADRC_Matabolomics.csv', header =T, sep=",", stringsAsFactors = F)
metab.HMDB.ID <- adad.res[, c('Metabolite.Name', 'HMDB.ID')]
#pathways.res <- as.data.frame(read_excel('omics_integration/data/Metabolomics/Supplementary_Tables.xlsx', sheet = 7))
#colnames(pathways.res) = pathways.res[2, ]

metab.res.other <- function(d, dt) {
  all.de.res <- NULL
  for (c in 1:3) {
    clu.res <- read.csv(paste0('omics_integration/data/Metabolomics/DE_results_metabolomics//03-metab_effect_pval_',c,'vs',dt,'.csv'), header =T, stringsAsFactors = F, sep=",")
    if (d == "INC") {
      clu.res <- clu.res$metab[clu.res$pval < 0.05 & clu.res$effect > 0]
    } else {
      clu.res <- clu.res$metab[clu.res$pval < 0.05 & clu.res$effect < 0]
    }
    all.de.res <- c(all.de.res, clu.res)
  }
   return(unique(all.de.res))
}

get_other_pathways <- function(d) {
    all.pathways.res <- NULL
    for (c in 1:3) {
      #c.res <- read.csv(paste0('omics_integration/data/Metabolomics/DE_results_metabolomics/Pathways/pathway_results_',c, 'vsall_',d,'.csv'), header =T, stringsAsFactors = F, sep=",")
      c.res <- read.csv(paste0('omics_integration/data/Metabolomics/DE_results_metabolomics/Pathways/pathway_results_',c, 'vsCO_',d,'.csv'), header =T, stringsAsFactors = F, sep=",")
      all.pathways.res <- c(all.pathways.res, c.res$X[c.res$Raw.p < 0.05])
    }
    return (all.pathways.res)
}

metab_16 <- c('aspartate', 'gamma-glutamylthreonine', 'beta-citrylglutamate', 'glutamate',
              'N-acetylglutamate', 'ergothioneine', '3-hydroxy-2-ethylpropionate', '1,5-anhydroglucitol (1,5-AG)',
              '2-methylcitrate/homocitrate', 'glutarate (C5-DC)', 'CDP-choline', 'CDP-ethanolamine',
              'glycerophosphoinositol*', 'nicotinamide', 'alpha-tocopherol', 'retinol (Vitamin A)')

for (d in c('INC', 'DEC')) {
  
  if (d=="INC") {
    Discovery <- knight_adrc$CHEMICAL.ID[knight_adrc$effect > 0]
    Replicate <- rosmap$CHEM_ID[rosmap$estimate > 0]
    shared.metab <- intersect(knight_adrc$CHEMICAL.ID[knight_adrc$effect > 0], rosmap$CHEM_ID[rosmap$estimate > 0])
    shared.metab <- knight_adrc$BIOCHEMICAL[knight_adrc$CHEMICAL.ID %in% shared.metab]
    inADvsCO <- knight_adrc$CHEMICAL.ID[knight_adrc$BIOCHEMICAL %in% KnightADRC_ADvsCO$metab_name[KnightADRC_ADvsCO$effect > 0]]
    
    shared.metab <- shared.metab[!shared.metab %in% knight_adrc$BIOCHEMICAL[knight_adrc$CHEMICAL.ID %in% inADvsCO] ]
    myCatNames <- c(paste0("KnightADRC (",grp1,")\n(",length(knight_adrc$CHEMICAL.ID[knight_adrc$effect > 0 & !knight_adrc$CHEMICAL.ID %in% inADvsCO]),")") , paste0("ROSMAP (", grp2,")\n(",length(rosmap$CHEM_ID[rosmap$estimate > 0]),")"))
    myCompGroup <- list(knight_adrc$CHEMICAL.ID[knight_adrc$effect > 0 & !knight_adrc$CHEMICAL.ID %in% inADvsCO], rosmap$CHEM_ID[rosmap$estimate > 0])
    mytitle <- 'Increased metabolites'
    ADAvsCO_res <- adad.res$Metabolite.Name[adad.res$q.value < 0.05 & adad.res$Effect > 0]
    p.res <- read.csv(paste0('omics_integration/Metabolomics_Analysis/', data.type, '/MetaboAnalyst/pathway_results_INC.csv'), sep=',', header =T, check.names = F, stringsAsFactors = F, row.names = 1)
  } else {
    Discovery <- knight_adrc$CHEMICAL.ID[knight_adrc$effect < 0]
    Replicate <- rosmap$CHEM_ID[rosmap$estimate < 0]
    shared.metab <- intersect(knight_adrc$CHEMICAL.ID[knight_adrc$effect < 0], rosmap$CHEM_ID[rosmap$estimate < 0])
    shared.metab <- knight_adrc$BIOCHEMICAL[knight_adrc$CHEMICAL.ID %in% shared.metab]
    inADvsCO <- knight_adrc$CHEMICAL.ID[knight_adrc$BIOCHEMICAL %in% KnightADRC_ADvsCO$metab_name[KnightADRC_ADvsCO$effect < 0]]
    
    shared.metab <- shared.metab[!shared.metab %in% knight_adrc$BIOCHEMICAL[knight_adrc$CHEMICAL.ID %in% inADvsCO] ]
    myCatNames <- c(paste0("KnightADRC(",grp1,")\n(",length(knight_adrc$CHEMICAL.ID[knight_adrc$effect < 0 & !knight_adrc$CHEMICAL.ID %in% inADvsCO]),")") , paste0("ROSMAP(",grp2,")\n(",length(rosmap$CHEM_ID[rosmap$estimate < 0]),")"))
    myCompGroup <- list(knight_adrc$CHEMICAL.ID[knight_adrc$effect < 0 & !knight_adrc$CHEMICAL.ID %in% inADvsCO], rosmap$CHEM_ID[rosmap$estimate < 0])
    mytitle <- 'Decreased metabolites'
    ADAvsCO_res <- adad.res$Metabolite.Name[adad.res$q.value < 0.05 & adad.res$Effect < 0]
    p.res <- read.csv(paste0('omics_integration/Metabolomics_Analysis/', data.type, '/MetaboAnalyst/pathway_results_DEC.csv'), sep=',', header =T, check.names = F, stringsAsFactors = F, row.names = 1)
  }
  
  if (length(shared.metab) !=0) {
    
    ss1 <- knight_adrc[knight_adrc$BIOCHEMICAL %in% shared.metab, c('BIOCHEMICAL', 'CHEMICAL.ID', 'effect', 'pval', 'padj')]
    colnames(ss1) <- c('Metab.name', 'CHEMICAL.ID', 'KnightADRC_effect', 'KnightADRC_pvalue', 'KnightADRC_padj')
    ss2 <- rosmap[rosmap$SHORT_NAME %in% shared.metab, c('estimate', 'p.value', 'padj')] 
    colnames(ss2) <- c('ROSMAP_effect', 'ROSMAP_pvalue', 'ROSMAP_padj')
    ss <- cbind(ss1, ss2)
    hmdb.ids <- adad.res[adad.res$Metabolite.Name %in% ss$Metab.name, c('Metabolite.Name', 'HMDB.ID')]
    colnames(hmdb.ids) <- c('Metab.name', 'HMDB.ID')
    ss <- merge(ss, hmdb.ids)
    
    ## check cluster specificity 
    others.res <- metab.res.other(d, dt)
    c4vsall.specific <- shared.metab[!shared.metab %in% others.res]
    ss$C4.Specific <- 'N'
    ss[ss$Metab.name %in% c4vsall.specific, 'C4.Specific'] <- 'Y'
    
    ## more results 
    ss$Shared.with.ADAD <- 'N'
    ss[ss$Metab.name %in% ADAvsCO_res, 'Shared.with.ADAD'] <- 'Y'
    ss$Shared.with.16Metab <- 'N'
    ss[ss$Metab.name %in% metab_16, 'Shared.with.16Metab'] <- 'Y'
    
    ## write resutls 
    write.table(ss, file = paste0('omics_integration/Metabolomics_Analysis/',data.type,'/', data.type, '_shared.metab_', d,'.tsv'), quote = F, sep="\t", row.names = F) 
    
    ## run hypergeometirc test 
    test.mat <- matrix(c( length(all.metab) - length(union(Discovery, Replicate)), length(setdiff(Discovery, Replicate)), length(setdiff(Replicate, Discovery)), length(intersect(Discovery, Replicate))), nrow=2)
    fisher.paval <- format(fisher.test(test.mat, alternative = "two.sided")$p.value, scientific = T, digits = 3)
    if (as.numeric(fisher.paval) ==0) { fisher.paval ='p < 2.2e-16'} else {fisher.paval = paste0('p = ',fisher.paval) }
    
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    myCol <- brewer.pal(2, "Dark2")
    vv <- venn.diagram( 
      x = myCompGroup,
      category.names = myCatNames, 
      filename = NULL, main = fisher.pval, main.fontface = "bold", 
      height = 300, width = 300 , resolution = 300, compression = "lzw",
      lwd = 2, lty = 1, main.cex = 1.2, 
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
    
    myleft <- textGrob(paste(shared.metab, collapse = "\n"), gp=gpar(fontsize=3, fontface="bold"))
    myright <- NULL
    if (length(shared.metab) != length(ss$Metab.name[ss$C4.Specific=='Y'])) { myright <- textGrob(paste(ss$Metab.name[ss$C4.Specific=='N'], collapse = "\n"), gp=gpar(fontsize=3, fontface="bold")) }
    mytop <- textGrob(paste(mytitle, collapse = "\n"), gp=gpar(fontsize=14, fontface="bold"))
    pdf(paste0('omics_integration/Metabolomics_Analysis/',data.type,'/', data.type, '_overlapped_metabs_', d,'.pdf'), width=2.7, height=2.7)
    #grid.draw(vv)
    grid.arrange(gTree(children=vv), top=mytop, left=myleft, right=myright)
    #grid.arrange(gTree(children=vv), top=mytop)
    dev.off()
  }
  
  ## plot readings 
  pdf(paste0('omics_integration/Metabolomics_Analysis/', data.type,'/', data.type, '_reading_plots_for_shared_matbs_', d,'.pdf'), width=4.3, height=3.5)
  for (m in shared.metab) {
    r.reading <- rosmap.metab[, c('Subj_ID', 'best.cluster', m)]
    colnames(r.reading) <- c('Subj_ID', 'best.cluster', 'reading')
    p = ggplot(r.reading, aes(x=as.factor(best.cluster), y=reading)) + geom_boxplot(aes(color=as.factor(best.cluster)), show.legend = F) +
      labs(x="", y="Reading") + ggtitle(m) + theme_bw() + 
      #geom_jitter(aes(color=as.factor(best.cluster)), position=position_jitter(0.3), size=1.5) +
      theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
            axis.text.x = element_text( vjust= 1, size=12, face="bold", color="black"), 
            axis.text.y = element_text(size=12, face="bold"),
            axis.title.y = element_text(size=12, face="bold"),
            legend.position = "none") +
      scale_color_manual(values=c("0"="gray80", "1"="#1B9E77", "2"="#D95F02")) + 
      scale_x_discrete(labels=c("0"=paste0("Control\n(n=",length(r.reading$Subj_ID[r.reading$best.cluster==0]),')'), 
                                "1"=paste0("Cluster1\n(n=",length(r.reading$Subj_ID[r.reading$best.cluster==1]),')'), 
                                "2"= paste0("Cluster2\n(n=",length(r.reading$Subj_ID[r.reading$best.cluster==2]),')')))
    print(p)
  }
  dev.off()
  
  ## extract C4 specfic pathways 
  p.res <- p.res[p.res$`Raw p` < 0.05, ]
  other.pathways <- get_other_pathways (d)
  p.res <- p.res[!rownames(p.res) %in% other.pathways, ]
  p.res$pathway_id <- rownames(p.res)
  p.res <- p.res[,c('pathway_id', colnames(p.res)[colnames(p.res) !="pathway_id"])]
  write.table(p.res, file = paste0('omics_integration/Metabolomics_Analysis/',data.type,'/', data.type, '_C4_specific_pathways_', d,'.tsv'), quote = F, sep="\t", row.names = F) 
  
}

# ## generate HMDB_IDs for patyhways analyses
# colnames(metab.HMDB.ID) <- c('metab', 'HMDB.ID')
# dd.type <- 'CO'
# #dd.type <- 'all'
# for (c in c(1,2,3,4)) {
#   cluster.res <- read.csv(paste0('omics_integration/data/Metabolomics/DE_results_metabolomics/03-metab_effect_pval_',c,'vs',dd.type,'.csv'), header =T, stringsAsFactors = F, sep=",")
#   cluster.res <- merge(cluster.res, metab.HMDB.ID)
#   cluster.res.up <- cluster.res$HMDB.ID[cluster.res$effect > 0 & cluster.res$pval < 0.05]
#   cluster.res.up <-  cluster.res.up[cluster.res.up!=""]
#   cluster.res.dn <- cluster.res$HMDB.ID[cluster.res$effect < 0 & cluster.res$pval < 0.05]
#   cluster.res.dn <-  cluster.res.dn[cluster.res.dn!=""]
#   cat(sapply(list(cluster.res.up, cluster.res.dn), length),'\n')
#   write.table(cluster.res.up, file=paste0('omics_integration/data/Metabolomics/DE_results_metabolomics/HMDB_IDs/hmdb_id_',c,'vs',dd.type,'_up.txt'), sep="\t", quote = F, row.names = F, col.names = F)
#   write.table(cluster.res.dn, file=paste0('omics_integration/data/Metabolomics/DE_results_metabolomics/HMDB_IDs/hmdb_id_',c,'vs',dd.type,'_dn.txt'), sep="\t", quote = F, row.names = F, col.names = F)
# }

## extract shared metabolities between intersections 
dec.metab.CO <- read.table('omics_integration/Metabolomics_Analysis/C4vsCO_C1vsCO/C4vsCO_C1vsCO_shared.metab_DEC.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t', quote = "")
inc.metab.CO <- read.table('omics_integration/Metabolomics_Analysis/C4vsCO_C1vsCO/C4vsCO_C1vsCO_shared.metab_INC.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t', quote = "")

dec.metab.OT <- read.table('omics_integration/Metabolomics_Analysis/C4vsall_C1vsC2/C4vsall_C1vsC2_shared.metab_DEC.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t', quote = "")
inc.metab.OT <- read.table('omics_integration/Metabolomics_Analysis/C4vsall_C1vsC2/C4vsall_C1vsC2_shared.metab_INC.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t', quote = "")

C4vsCO_OV_INC <- intersect(inc.metab.CO$Metab.name, inc.metab.OT$Metab.name)
C4vsCO_OV_DEC <- intersect(dec.metab.CO$Metab.name, dec.metab.OT$Metab.name)


## extract pavalues 
mm <- 'glucose'
dec.metab.CO[dec.metab.CO$Metab.name==mm, c('Metab.name', 'KnightADRC_pvalue', 'ROSMAP_pvalue')]
dec.metab.OT[dec.metab.OT$Metab.name==mm, c('Metab.name', 'KnightADRC_pvalue', 'ROSMAP_pvalue')]

rbind(CO=inc.metab.CO[inc.metab.CO$Metab.name==mm, c('Metab.name', 'KnightADRC_effect', 'ROSMAP_effect', 'KnightADRC_pvalue', 'ROSMAP_pvalue')],
      OT=inc.metab.OT[inc.metab.OT$Metab.name==mm, c('Metab.name', 'KnightADRC_effect', 'ROSMAP_effect', 'KnightADRC_pvalue', 'ROSMAP_pvalue')])


##########################################################
##### Plot Supplementary Figures 
#########################################################
# Knight ADRC 
## read cluster solution
bb1 <- read.table("omics_integration/iCluster_output/best.cluster.memership_v1_full.tsv", header = T, sep="\t", stringsAsFactors = F)
## remove flagged samples 
bb1 <- bb1[!bb1$Subj_ID %in% as.integer(gsub("ID_","", flagged.samples)), ]
bb1[bb1$Status=="Neuro_CO", 'best.cluster'] <- 0

##### read metab data 
metab <- read.table('omics_integration/data/QCd_metabolites_for_RNAseq_brains.csv', header =T, sep=',', stringsAsFactors = F, check.names = F)
metab.data.cols <- colnames(metab)[14:ncol(metab)]
metab <- metab[, c('BIOCHEMICAL', 'COMP.ID', metab.data.cols)]
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
  r.min <- rowMeans(metab[i, -1], na.rm = T)
  metab[i, indx] <- r.min
}
colnames(metab)[colnames(metab)=="ID_BIOCHEMICAL"] <- "BIOCHEMICAL"

metab.name <- 'NAD+' 
rr <- reshape2::melt(metab[metab$BIOCHEMICAL==metab.name, gsub("ID_","", colnames(metab)) %in% bb1$Subj_ID], variable.name="Subj_ID", value.name="reading")
rr$Subj_ID <- gsub("ID_", "", rr$Subj_ID)
rr <- merge(bb1[, c('Subj_ID', 'best.cluster')], rr, all.x = T)

C4vsall.m <- read.csv('omics_integration/data/Metabolomics/DE_results_metabolomics//03-metab_effect_pval_4vsall.csv', header =T, stringsAsFactors = F, sep=",")
C4vsCO.m <- read.csv('omics_integration/data/Metabolomics/DE_results_metabolomics/03-metab_effect_pval_4vsCO.csv', header =T, stringsAsFactors = F, sep=",")

pvals <- data.frame(comp=c('C4vsCO','C4vsAll'), 
                    pval = c(C4vsCO.m[C4vsCO.m$metab==metab.name, 'pval'], 
                             C4vsall.m[C4vsall.m$metab==metab.name, 'pval']))
pvals$pval <- format(pvals$pval, scientific = T, digits = 2)

kp <- ggplot(rr, aes(x=as.factor(best.cluster), y=reading)) + geom_boxplot(aes(color=as.factor(best.cluster)), show.legend = F) +
      labs(x="", y="Reading") + ggtitle('Knight ADRC') + theme_bw() + 
      geom_jitter(aes(color=as.factor(best.cluster)), position=position_jitter(0.3), size=1.5) +
      theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
          axis.text.x = element_text( vjust= 1, size=14, face="bold", color="black"), 
          axis.text.y = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=16, face="bold"),
          panel.background = element_rect(colour = "black", size=1),legend.position = "none")  +
      scale_color_manual(values=c("0"="gray80", "1"="#1B9E77", "2"="#D95F02", "3"="#7570B3", "4"="#E7298A")) + 
      scale_x_discrete(labels=c("0"=paste0("Control\n(n=",length(rr$Subj_ID[rr$best.cluster==0]),')'), 
                            "1"=paste0("Cluster1\n(n=",length(rr$Subj_ID[rr$best.cluster==1]),')'), 
                            "2"= paste0("Cluster2\n(n=",length(rr$Subj_ID[rr$best.cluster==2]),')'),
                            "3"=paste0("Cluster3\n(n=",length(rr$Subj_ID[rr$best.cluster==3]),')'), 
                            "4"= paste0("Cluster4\n(n=",length(rr$Subj_ID[rr$best.cluster==4]),')'))) + 
      geom_signif(data = pvals, aes(xmin = 1, xmax = 5, annotations =  paste0('p=',pvals$pval[pvals$comp=="C4vsCO"]), y_position = max(rr$reading, na.rm=T)+0.1), 
              textsize = 5.5, size=0.6, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") + 
      geom_signif(data = pvals, aes(xmin = 1.6, xmax = 4.4, annotations ='', y_position = max(rr$reading, na.rm=T)), 
              textsize = 5.5, size=0.6, step_increase=0.1, vjust = -0.2, manual = T, fontface="bold", lty="dashed") +
      geom_signif(data = pvals, aes(xmin = 3, xmax = 5, annotations = paste0('p=',pvals$pval[pvals$comp=="C4vsAll"]), y_position = max(rr$reading, na.rm=T)+0.05), 
              textsize = 5.5, size=0.6, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") 


# ROSMAP 
## get pvalues from DE analysis 
r.c1.vs.all <- read.table('omics_integration/Replication/ROSMAP/Metab_DE/1vs2/1vs2.res.all.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t', row.names = 1, quote = "")
r.c1.vs.co <- read.table('omics_integration/Replication/ROSMAP/Metab_DE/1vs0/1vs0.res.all.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t', row.names = 1, quote = "")
r.c2.vs.co <- read.table('omics_integration/Replication/ROSMAP/Metab_DE/2vs0/2vs0.res.all.tsv', header =T, stringsAsFactors = F, check.names = F, sep='\t', row.names = 1, quote = "")

sig <- data.frame(comp=c('C1vsCO','C2vsCO','C1vsC2'),
                    pvalue = c(r.c1.vs.co[rownames(r.c1.vs.co) ==metab.name, 'p.value'],
                             r.c2.vs.co[rownames(r.c2.vs.co) ==metab.name, 'p.value'],
                             r.c1.vs.all[rownames(r.c1.vs.all) ==metab.name, 'p.value']))
sig$pvalue <- format(sig$pvalue, scientific = T, digits = 2)

r.reading <- rosmap.metab[, c('Subj_ID', 'best.cluster', metab.name)]
colnames(r.reading) <- c('Subj_ID', 'best.cluster', 'reading')

rp <- ggplot(r.reading, aes(x=as.factor(best.cluster), y=reading)) + geom_boxplot(aes(color=as.factor(best.cluster)), show.legend = F) +
          labs(x="", y="") + ggtitle('ROSMAP') + theme_bw() +
          geom_jitter(aes(color=as.factor(best.cluster)), position=position_jitter(0.3), size=1.5) +
          theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
                axis.text.x = element_text( vjust= 1, size=14, face="bold", color="black"),
                axis.text.y = element_text(size=14, face="bold"),
                axis.title.y = element_text(size=16, face="bold"), legend.position = "none")  +
          scale_color_manual(values=c("0"="gray80", "1"="#E7298A", "2"="#00AFBB")) +
          scale_x_discrete(labels=c("0"=paste0("Control\n(n=",length(r.reading$Subj_ID[r.reading$best.cluster==0]),')'),
                            "1"=paste0("Cluster1\n(n=",length(r.reading$Subj_ID[r.reading$best.cluster==1]),')'),
                            "2"= paste0("Cluster2\n(n=",length(r.reading$Subj_ID[r.reading$best.cluster==2]),')'))) +
          geom_signif(data = sig, aes(xmin = 1, xmax = 1.97, annotations =  paste0('p=',pvalue[comp=="C1vsCO"]), y_position = max(r.reading$reading, na.rm=T)+0.03),
              textsize = 5.5, size=0.6, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") +
          geom_signif(data = sig, aes(xmin = 1, xmax = 3, annotations =  paste0('p=',pvalue[comp=="C2vsCO"]), y_position =max(r.reading$reading, na.rm=T)+0.065),
              textsize = 5.5, size=0.6, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold") +
          geom_signif(data = sig, aes(xmin = 2.03, xmax = 3, annotations =  paste0('p=',pvalue[comp=="C1vsC2"]), y_position =max(r.reading$reading, na.rm=T)+0.03),
              textsize = 5.5, size=0.6, step_increase=0.1, vjust = -0.2,manual = T, fontface="bold")

pdf(paste0('omics_integration/Metabolomics_Analysis/Figures/', metab.name, '_reading_v2.pdf'), width=12, height = 7)
grid.arrange(arrangeGrob(grobs=list(kp, rp), ncol=2, widths=c(3,1.8)))
#kp
dev.off()
