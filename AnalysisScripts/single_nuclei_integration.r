setwd('/home/eteleeb/projects')
library(ggplot2)
library(ggpubr)
library(ggsci)

dir.create('omics_integration/iCluster_output/cell_specific_pathways')

folder_names <- list.dirs('/home/delaguilaj/DE_soma/', full.names = F, recursive = F)

## function to remove shared genes 
comps <- c('C1', 'C2', 'C3', 'C4')
remove_shared_genes <- function (cc) {
  cc <- gsub('vsCO', '', cc)
  rr <- NULL
  for (cmp in comps[comps != cc]) {
    res <- read.table(paste0('omics_integration/iCluster_output/DE/',cmp,'vsCO/de_res_sig_',cmp,'vsCO.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')
    rr <- c(rr, res$GeneName)
  }
  return(rr)
}

#for (dd.type in c('pval', 'est', 'zval')) {
#  cat(dd.type,'\n')
#  if (dd.type =="pval") {
#    var <- 'PValue'
#    ll <- '-Log10(P-values)'
#    tt <- 'pvalues'
#  } else if (dd.type == 'est') {
    dd.type <- 'est'
    var <- 'estimate'
    ll <- 'Estimate'
    tt <- 'estimate'
#  } else if (dd.type =='zval') {
#    var <- 'ZVvalue'
#    ll <- 'Z-Vvalue'
#    tt <- 'zvalues'
#  }
  
  ## loop through the whole data 
  for (d in c('all', 'up', 'dn')) {
  
    ff.plots <- list()
    tmp = NULL
    for (c in c('C1vsCO', 'C2vsCO', 'C3vsCO', 'C4vsCO')) {
      cat('   Running for: ', c, '\n')
      rr1 <- read.table(paste0('omics_integration/iCluster_output/DE/', c, '/de_res_sig_', c,'.tsv'), header =T, stringsAsFactors = F, check.names = F, sep='\t')
      #dir.create(paste0('omics_integration/iCluster_output/cell_specific_pathways/', c))
      
      #rr1 <- rr1[abs(rr1$log2FoldChange) > 1, ]
      
      ## remove shared genes 
      if (d=="dn") {
        shared_genes <- remove_shared_genes(c)
        rr1 <- rr1[!rr1$GeneName %in% shared_genes, ]    
      }
    
      
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
        
        #if (dd.type =="est") {
          upper_cutoff <- quantile(cluster.data$estimate, prob=0.90)
          cluster.data <- cluster.data[cluster.data$estimate > upper_cutoff, ]
        #} else if (dd.type =="zval") {
        #  upper_cutoff <- quantile(cluster.data$ZVvalue, prob=0.90)
        #  cluster.data <- cluster.data[cluster.data$ZVvalue > upper_cutoff, ]
        #}
        
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
      #if (dd.type == "pval") {
      #  p <-  ggplot(folder.data, aes(x=factor(cluster), y= -log10(get(var)), group = cluster)) + theme_bw(base_size = 12)
      #} else {
      p <-  ggplot(folder.data, aes(x=factor(cluster), y= get(var), group = cluster)) + theme_bw(base_size = 12)
      #}
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
    }
   
    ## write results 
    write.table(tmp, paste0('omics_integration/iCluster_output/cell_specific_pathways/', d, '_', dd.type, '_top10_cluster_specific.tsv'), sep="\t", row.names = F, quote = F)
    
    ## plot all 
    mytitle = textGrob(paste0('Distribution of the ', tt,' for the top 10 DE cell type - (', d, ')\n'), gp=gpar(fontsize=20, fontface="bold"))
    myleft = textGrob(ll, gp=gpar(fontsize=20, fontface="bold"), rot = 90)
    pdf(paste0('omics_integration/iCluster_output/cell_specific_pathways/', d, '_', dd.type, '_top10_cluster_specific.pdf'), width = 13, height = 10)
    grid.arrange(arrangeGrob(grobs=ff.plots, ncol=2, top=mytitle, left= myleft))
    dev.off()
  
  }
  
#}


######################################################################################################################################
tmp$index <- paste0('index_', 1:nrow(tmp))
keep <- NULL
for (c in c('C1vsCO', 'C2vsCO', 'C3vsCO', 'C4vsCO')) {
 
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

## extract keep genes 
tmp <- tmp[tmp$index %in% keep, ]
gene.counts <- aggregate(GeneName ~ cluster + comp, data = tmp, FUN = length)
colnames(gene.counts) <- c('cluster', 'comp', 'num.genes')
totals <- aggregate(. ~ comp, data = gene.counts, sum )
gene.counts$pct <- 'NA'
gene.counts[gene.counts$comp=='C1vsCO', 'pct'] <- (gene.counts$num.genes[gene.counts$comp=='C1vsCO'])/(totals$num.genes[totals$comp=='C1vsCO'])
gene.counts[gene.counts$comp=='C2vsCO', 'pct'] <- (gene.counts$num.genes[gene.counts$comp=='C2vsCO'])/(totals$num.genes[totals$comp=='C2vsCO'])
gene.counts[gene.counts$comp=='C3vsCO', 'pct'] <- (gene.counts$num.genes[gene.counts$comp=='C3vsCO'])/(totals$num.genes[totals$comp=='C3vsCO'])
gene.counts[gene.counts$comp=='C4vsCO', 'pct'] <- (gene.counts$num.genes[gene.counts$comp=='C4vsCO'])/(totals$num.genes[totals$comp=='C4vsCO'])
gene.counts$pct <- as.numeric(gene.counts$pct)
gene.counts$comp <- gsub('vsCO', '', gene.counts$comp)
gene.counts$pct <- gene.counts$pct*100

## for final figure 
p <- ggplot(gene.counts, aes(x=factor(comp), y= pct)) + theme_bw(base_size = 12) + 
        geom_bar(stat="identity", aes(fill=factor(cluster)), width=0.5, color="black") + #coord_polar(theta = "y") + 
        labs(x='', y = 'Percentage of DE genes', title ="Underexpressed hits") +
        theme(axis.text.x=element_text(size=18, vjust=1, hjust =1, color="black", angle=45),
            axis.text.y=element_text(size=18, color="black"),
            axis.title.x=element_text(size=19, color="black"),
            axis.title.y=element_text(size=19, color="black"),
            #panel.background=element_rect(color="black"),
            plot.title = element_text(size = 20, hjust=0.5, color="black", face="bold"), 
            legend.position="right", legend.title = element_text(size=14), legend.text=element_text(size=14)) +
        scale_fill_jco() + 
        scale_x_discrete(labels=c("C1"="Knight-C1", "C2"="Knight-C2", "C3"="Knight-C3", "C4"="Knight-C4"))# scale_fill_viridis_d(name="") + #scale_fill_brewer(palette = 'Dark2' )

pdf(paste0('omics_integration/iCluster_output/Figures/Final/cell_types/percentage_hits_per_cell_type_', d,'_v2.pdf'), width = 5, height = 5)
p
dev.off()
        
######################################################################################################################################




#####################################################################################################################################
## run patheays 
#####################################################################################################################################
#tmp <-  read.table(paste0('omics_integration/iCluster_output/cell_specific_pathways/', d, '_', dd.type, '_top10_cluster_specific.tsv'), header =T, stringsAsFactors = F)
int.dbs <- c('KEGG_2019_Human', 'GO_Biological_Process_2018', 'GO_Molecular_Function_2018', 'GO_Cellular_Component_2018')
gg= as.character(tmp[tmp$comp=="C4vsCO" & tmp$cluster=="Astrocytes", 'GeneName'])
ee = enrichr(gg, databases = int.dbs)  
#ee <- ee$GO_Biological_Process_2018
#ee$Overlap <- as.numeric(sub("/.*", "", ee$Overlap))
#ee <- ee[ee$P.value < 0.05 & ee$Overlap > 5, ]
all.terms <- list()
for (term in int.dbs) {
  #plotEnrich(ee, title ='Top 20 GO-BP results ordered by P-value')
  p <- plotEnrich(ee[[term]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", xlab = '', ylab = '', title = paste0('Top 20 for ', term))
  all.terms[[term]] <- p
}

mytop <- textGrob('Oligo-C2 only (Up-regulated genes only)', gp=gpar(fontsize=16, fontface="bold"))
pdf(paste0('omics_integration/iCluster_output/cell_specific_pathways/top_pathways_Oligo_C2_', d, '_', dd.type, '_top10.pdf'), width = 18, height = 10)
grid.arrange(grobs = all.terms, nrow=2, ncol=2, top=mytop)
dev.off()


## generate similarity matrix 
# go_ids <- ee$Term
# if (nrow(go_ids) > 0) {
#   go_ids <- gsub('.*\\(', '', go_ids)
#   go_ids <- gsub('\\)', '', go_ids)
#   mat = GO_similarity(go_ids, toupper(cat))
#   go_sim_p = simplifyGO(mat)
#   pdf(paste0('omics_integration/iCluster_output/pathways_analysis/',cc, '/pathways_plots/go_sim_', d, '_', c,'.pdf'), width = 8, height = 5)
#   simplifyGO(mat, word_cloud_grob_param = list(max_width = 50))
#   dev.off()
# }

#######################################################################
## check cell-type specific pathways 
#######################################################################
for (d in c('up', 'dn')) {
  tmp <-  read.table(paste0('omics_integration/iCluster_output/cell_specific_pathways/', d, '_', dd.type, '_top10_cluster_specific.tsv'), header =T, stringsAsFactors = F)
  
  for (c in c('C1vsCO', 'C2vsCO', 'C3vsCO', 'C4vsCO')) {
    dir.create(paste0('omics_integration/iCluster_output/cell_specific_pathways/', c))
    rm(en.up)
    rm(en.dn)
    load(paste0('omics_integration/iCluster_output/DE/', c, '/', c, '_enrichR_res_up.RData'))
    load(paste0('omics_integration/iCluster_output/DE/', c, '/', c, '_enrichR_res_dn.RData'))
    
    for (cell in unique(tmp$cluster)) {
      gg <- as.character(tmp[tmp$comp==c & tmp$cluster==cell, 'GeneName'])
      prepare_gsea_data(gg)
      
      if (length(gg) > 0) {
        dir.create(paste0('omics_integration/iCluster_output/cell_specific_pathways/', c,'/', cell))
        cat(paste0('Running for ', d, '_',c, '_', cell,'\n'))
        plot_top_pathways(d, c, cell, en.up, en.dn, gg)
      } else {
        cat(paste0('No genes were found for ',  d, '_',c, '_', cell))
        next 
      }
      
    }  ## end of cells 
    
  }
}

# c <- 'C1vsCO'
# d <- 'dn'
# cell <- 'Microglia' # Astrocytes\Neuron\Microglia
# dir.create(paste0('omics_integration/iCluster_output/cell_specific_pathways/', c))
# dir.create(paste0('omics_integration/iCluster_output/cell_specific_pathways/', c,'/', c))
# 
# rm(en.up)
# rm(en.dn)
# load(paste0('omics_integration/iCluster_output/DE/', c, '/', c, '_enrichR_res_up.RData'))
# load(paste0('omics_integration/iCluster_output/DE/', c, '/', c, '_enrichR_res_dn.RData'))
# 
# tmp <-  read.table(paste0('omics_integration/iCluster_output/cell_specific_pathways/', d, '_', dd.type, '_top10_cluster_specific.tsv'), header =T, stringsAsFactors = F)
# gg <- as.character(tmp[tmp$comp==c & tmp$cluster==cell, 'GeneName'])

plot_top_pathways <- function (d, c, cell, en.up, en.dn, gg) {
  all.terms <- NULL
  for (pathway in int.dbs) {
    
    if (d =="up") { en.res <- en.up[[pathway]] } else { en.res <- en.dn[[pathway]] }
    en.res <- en.res[en.res$P.value < 0.05, ]
    d.res <- NULL
    p.gene.dis.plots <- list()
    for (i in 1:nrow(en.res)) {
      p.genes <- unlist(strsplit(en.res[i, 'Genes'], ";"))
      dd <- data.frame(term = en.res$Term[i], pval = en.res$P.value[i], adj.pval = en.res$Adjusted.P.value[i], 
                       total.ov.genes = paste(p.genes, collapse = "|"), num.shared.genes =  length(intersect(p.genes, gg)), 
                       pct.shared.genes =  paste0(signif((length(intersect(p.genes, gg))/length(gg))*100, digits = 3)),    
                       shared.genes =  paste(intersect(p.genes, gg), collapse = "|"), stringsAsFactors = F)
      if (dd$num.shared.genes > 0) { 
        d.res <- rbind(d.res, dd)
        ## generate pathway genes dist plot 
        #p2 <- plot_p.genes.dis(en.res$Term[i], p.genes, length(intersect(p.genes, gg)), c )
        #p.gene.dis.plots[[en.res$Term[i]]] <- p2
      } else {
        next
      }
    }
    
    ## order pathways 
    if (!is.null(d.res)) {
      d.res <- d.res[order(d.res$num.shared.genes, decreasing = T), ]
      if (nrow(d.res) > 20) { d.res <- d.res[1:20, ]}
      
      ## plot the pathway genes dist 
      p.gene.dis.plots <- p.gene.dis.plots[!sapply(p.gene.dis.plots, is.null)]
      p.gene.dis.plots <- p.gene.dis.plots[names(p.gene.dis.plots) %in% d.res$term]
      mytop <- textGrob( paste0('Genes distribution of the top pathways for ', cell, ' genes (',c,'-', d,')\n', pathway), gp=gpar(fontsize=16, fontface="bold"))
      mybb <- textGrob( paste0('Percentage of genes'), gp=gpar(fontsize=16, fontface="bold"))
      nywidth <- ifelse(length(p.gene.dis.plots) <= 4, 15, length(p.gene.dis.plots))
      myheight <- ifelse(length(p.gene.dis.plots) <= 4, 8, length(p.gene.dis.plots)-5)  
      pdf(paste0('omics_integration/iCluster_output/cell_specific_pathways/',c,'/', cell, '/', cell, '_genes_dist_', pathway, '_', c, '_', d, '.pdf'), width = nywidth, height = myheight)
      grid.arrange(grobs = p.gene.dis.plots, ncol=ifelse( length(p.gene.dis.plots) >= 4, 4,length(p.gene.dis.plots)), top=mytop, bottom = mybb)
      dev.off()
      
      ## plot the top 20 
      d.res$term <- gsub('\\(.*', '', d.res$term)
      d.res$term <- with(d.res, factor(d.res$term, levels=rev(unique(d.res$term))))
      # p <- ggplot(d.res, aes(x=term, y=num.shared.genes)) + geom_bar(stat="identity", aes(fill=pval)) + coord_flip() + theme_bw() + 
      #   labs(x='', y = 'Percentage of genes', title = paste0('Top pathways for ', cell, ' genes (',c,'-', d,')\n', pathway)) +
      #   theme(axis.text.x=element_text(size=12, vjust=0.5, face="bold"),
      #         axis.text.y=element_text(size=12, face="bold"),
      #         axis.title.x=element_text(size=16, face="bold"),
      #         plot.title = element_text(size = 14, hjust=0.5, face="bold"),
      #         axis.ticks=element_blank(), panel.grid.minor = element_blank(),
      #         legend.title=element_text(size=16, face="bold"), legend.text=element_text(size=16, face="bold"),
      #         panel.background=element_blank()) + scale_fill_viridis_c()
      
      ## as points 
      p <- ggplot(d.res, aes(term, num.shared.genes)) + geom_point(aes(size=-log10(pval)), colour="blue", stroke=1.5) + coord_flip() +
                labs(x='', y='Number of genes', title='') + theme_bw(base_size=12) +
                theme(axis.text.x=element_text(size=16, vjust=0.5, face="bold", color="black"),
                    axis.text.y=element_text(size=18, face="bold", color="black", vjust=0.3),
                    axis.title.x=element_text(size=18, face="bold", color="black"),
                    plot.title = element_text(size = 14, hjust=0.5, face="bold"),
                    axis.ticks=element_blank(), panel.grid.minor = element_blank(),
                    legend.title=element_text(size=14, face="bold"), legend.text=element_text(size=14, face="bold"),
                    panel.background=element_blank()) + scale_size_continuous(name="-Log10(PValue)") +
                scale_y_continuous(breaks = c(seq(min(d.res$num.shared.genes), max(d.res$num.shared.genes), 1)), limits=c(min(d.res$num.shared.genes), max(d.res$num.shared.genes)) )
      
      ### for manuscript ###
      pdf(paste0('omics_integration/iCluster_output/Figures/Final/cell_types/Top_20_',cell,'_',pathway,'_',c,'_',d,'.pdf'), width = 12.3, height = 6)
      p
      dev.off()
      ###
      all.terms[[pathway]] <- p
    }
  }
  
  glist <- lapply(all.terms, ggplotGrob)
  multi.page <- ggpubr::ggarrange(plotlist = glist, nrow = 1, ncol = 1, widths = 10, heights = 10) 
  ggpubr::ggexport(multi.page, filename = paste0('omics_integration/iCluster_output/cell_specific_pathways/',c, '/', cell, '/', cell, '_genes_top_pathways_', c, '_', d, '_top10.pdf'), width = 10, height = 10)
  
}

###### plot gene distriubtion 
plot_p.genes.dis <- function (pathway, p.genes, ov, c) {
  ## chech the distriubtion of genes in other cells 
  p.genes.dis <- tmp[tmp$comp==c & tmp$GeneName %in% p.genes, ]
  p.genes.dis <- aggregate(GeneName ~ cluster, data = p.genes.dis, length)
  colnames(p.genes.dis) <- c('Cluster', 'Num.Genes')
  diff <-  length(p.genes) - sum(p.genes.dis[,2])
  if (diff > 0) {
    p.genes.dis <- rbind(p.genes.dis, data.frame(Cluster ="NonAvail",  Num.Genes = diff))
  }
  p.genes.dis$PCT <- (p.genes.dis$Num.Genes/length(p.genes)) * 100
  p.genes.dis <- p.genes.dis[order(p.genes.dis$PCT, decreasing = T), ]
  p.genes.dis$Cluster <- with(p.genes.dis, factor(p.genes.dis$Cluster, levels=rev(unique(p.genes.dis$Cluster))))
  
  p2 <- ggplot(p.genes.dis, aes(x=Cluster, y=PCT)) + 
    geom_bar(stat="identity", width = 0.8, fill=ifelse(p.genes.dis$Cluster=="NonAvail", "gray80", "#440154FF")) + 
    coord_flip() + theme_bw() + labs(x='', y = '', title = paste0(pathway,'\n(n=', length(p.genes), ')')) +
    geom_text(label=p.genes.dis$Num.Genes, hjust = 0.2, size=5) + 
    theme(axis.text.x=element_text(size=10, vjust=0.5, face="bold"),
          axis.text.y=element_text(size=10, face="bold"),
          plot.title = element_text(size=12, hjust=0.5, face="bold"), 
          axis.ticks=element_blank(), panel.grid.minor = element_blank(),
          panel.background=element_blank(), legend.position = 'none')  
    #scale_y_continuous(breaks = seq(0,max(p.genes.dis$PCT),2) ) 
  
  return (p2)
}









################################################################################
### generate data fir GSEA analysis 
################################################################################
bb <- read.table("omics_integration/iCluster_output/best.cluster.memership_v1_full.tsv", header = T, sep="\t", stringsAsFactors = F)
bb$group <- 'CO'
bb[bb$best.cluster == 1 & bb$Status =='Neuro_AD', 'group'] <- 'C1'
bb[bb$best.cluster == 2 & bb$Status =='Neuro_AD', 'group'] <- 'C2'
bb[bb$best.cluster == 3 & bb$Status =='Neuro_AD', 'group'] <- 'C3'
bb[bb$best.cluster == 4 & bb$Status =='Neuro_AD', 'group'] <- 'C4'

## read count & expression 
rnaseq <- read.table('omics_integration/data/brain_exp_matrix_FPKM_2nd_read_strand.tsv', header =T, stringsAsFactors = F, check.names = F)
mean.expr = rowMeans(rnaseq[, -c(1,2)], na.rm = T)
rnaseq = rnaseq[order(mean.expr, decreasing=T),]
rnaseq = rnaseq[!duplicated(rnaseq[["GeneName"]]),]
rownames(rnaseq) <- rnaseq$GeneName
rnaseq$GeneID <- NULL
rnaseq$GeneName <- NULL
colnames(rnaseq) <- gsub('ID_','', colnames(rnaseq))
rnaseq <- rnaseq[, colnames(rnaseq) %in% bb$Subj_ID]

tt <- as.data.frame(t(rnaseq))
tt$Subj_ID <- rownames(tt)
tt <- merge(bb[, c('Subj_ID', 'group')], tt)

prepare_gsea_data <- function (gg) {
  
  data.subset = rnaseq[gg,
    c(as.character(bb$Subj_ID[bb$Status=="Neuro_CO"]),
      as.character(bb$Subj_ID[bb$best.cluster==4 & bb$Status=="Neuro_AD"])
    )]
  
  dd = colnames(data.subset)
  data.subset$NAME = rownames(data.subset)
  rownames(data.subset) = NULL
  data.subset$DESCRIPTION = NA
  data.subset = data.subset[, c('NAME', 'DESCRIPTION', dd)]
  
  ph <- c(rep('CO', table(bb$group)[["CO"]]), rep('C4', table(bb$group)[["C4"]]))
  ph <- paste(ph, collapse = " ")
  
  write.table(data.subset, file=paste0('omics_integration/iCluster_output/cell_specific_pathways/',c, '/', cell, '/', cell, '_genes_for_GSEA_', c, '_', d, '.gct'), sep="\t", quote = F, row.names = F)
  write.table(ph, file=paste0('omics_integration/iCluster_output/cell_specific_pathways/',c, '/', cell, '/', cell, '_phenotype_class_', c, '_', d, '.cls'), quote = F, row.names = F)
  
  #ph = c(rep('Control', 23), rep('Cluster123', 114), rep('Cluster123', 46), rep('Cluster123', 53) ,rep('Cluster4', 45) )
  #ph = c(rep('Cluster123', 114), rep('Cluster123', 46), rep('Cluster123', 53) ,rep('Cluster4', 45) )

  
}

write.table(data.subset, file=paste0('omics_integration/iCluster_output/cell_specific_pathways/',c, '/', 'homeostatic_genes_for_GSEA_', c, '_', d, '.gct'), sep="\t", quote = F, row.names = F)




## combine all sets 
folder.data <- NULL
for (clu in c(0,1,2,3,4,5)) {
  cluster.data <- NULL
  for (ff in folder_names) {
    sn.data <- read.table(list.files(paste0('/home/delaguilaj/DE_soma/', ff) , pattern=paste0('glmmTMB_celltype_', clu), full.names = T), sep="\t", header =F, stringsAsFactors = F)
    colnames(sn.data) <- c('GeneName', 'estimate', 'Std.Error', 'ZVvalue', 'PValue')
    sn.data$cell.type <- clu
    cluster.data <- rbind(cluster.data, sn.data)
  }
  
  folder.data <- rbind(folder.data, cluster.data)
  
}

## plot the results 
folder.data[folder.data$cell.type==0, 'cell.type'] <- 'Oligo'
folder.data[folder.data$cell.type==1, 'cell.type'] <- 'Neuron'
folder.data[folder.data$cell.type==2, 'cell.type'] <- 'Microglia'
folder.data[folder.data$cell.type==3, 'cell.type'] <- 'Astrocytes'
folder.data[folder.data$cell.type==4, 'cell.type'] <- 'Opc'
folder.data[folder.data$cell.type==5, 'cell.type'] <- 'Endo'

write.table(folder.data, file='omics_integration/data/Jorge_DE_per_cell_types.tsv', sep="\t", quote = F, row.names = F)

