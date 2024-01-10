## Survival for gene
library(plyr)
library(dplyr)
library(survival)
library(surviplot)
library(ggplot2)
library(ggsignif)

setwd("/home/bnovotny/proteomics/proteomics_analysis/CSF")

## Read in protein data
features <- readRDS("data/00-csf_feature_data.rds")

# Read in protein expression data and reconfigure
log_readings <- readRDS("data/01-csf_ComBat_1254x713.rds") # Subjects in rows
log_readings <- data.frame(log_readings)
all(features$SOMAseqID == colnames(log_readings))
colnames(log_readings) <- features$SOMAseqID
log_readings$ExtIdentifier <- row.names(log_readings)

# Read in phenotype data and reconfigure
pheno <- readRDS("data/01-csf_pheno_no43.rds")
pheno$Status <- as.character(pheno$Status)
pheno$Status[pheno$trem2 == TRUE] <- "TREM2"
pheno$Status <- factor(pheno$Status, levels = c("CO", "AD", "ADAD", "TREM2", "FTD", "OT"))

# Put together relevant pheno data and expression data
pheno_small <- pheno %>% select(ExtIdentifier, MAPID, Status, LP.num, Age.at.LP)

cdr_data <- read.csv("data/Summary_Proteomics_Somalogic_final_V4_CYcleaned_sheet6_CDR.csv", stringsAsFactors = F)
# Make sure CDR measurements are ordered correctly
cdr_data <- cdr_data %>% group_by(ID) %>% arrange(age_at_CDR_SB_day, .by_group = TRUE)

cut_to_LP_date <- function(id) {
  df <- cdr_data[cdr_data$ID == id,] # Get all CDR observations for one individual
  if(!all(is.na(df$cdr)) & any(df$cdr >= 1)) { # If any CDRs exist and are greater than 1
    # df
    df <- filter(df, df$age_at_CDR_SB_year >= pheno$Age.at.LP[pheno$MAPID == id & pheno$LP.num == "p1"]) # Take all LPs greater than or equal to age at LP1
    if (!all(is.na(df$cdr)) & df$cdr[1] == 0) df#[1:min(which(df$cdr > 0)),] # If the first CDR is 0 keep
    else NULL
  } else if (!all(is.na(df$cdr)) & all(df$cdr == 0)) df # Also keep people who never converted
  else NULL # Otherwise delete
}


cdr_cut <- ldply(unique(pheno_small$MAPID), cut_to_LP_date)
cdr_cut_small <- cdr_cut %>% select("ID", "cdr", "age_at_CDR_SB_year")

# Put expression, cdr and pheno info together
pheno_small_cdr <- merge(pheno_small, cdr_cut_small, by.x = "MAPID", by.y = "ID") #############



run_survival <- function(gene) {
  gene_df <- data.frame(ExtIdentifier = row.names(log_readings), gene = log_readings[,gene])
  model_data_gene <- merge(pheno_small_cdr, gene_df, by = "ExtIdentifier")
  
  # model_data_gene <- model_data_gene %>% filter(LP.num == "p1")
  # Add time elapsed between LP and CDR
  model_data_gene$time_elapsed <- model_data_gene$age_at_CDR_SB_year - model_data_gene$Age.at.LP
  model_data_gene <- filter(model_data_gene, time_elapsed >= 0)
  
  # Add gene expression level
  model_data_gene$expression <- "> median"
  model_data_gene$expression[model_data_gene$gene <= median(log_readings[,gene])] <- "<= median"
  
  if(length(unique(model_data_gene$expression)) > 1) {
    # Each person can only have one expression level, so change all expression for
    # longitudinal subjects to their expression level at p1
    long_ids <- pheno_small$MAPID[pheno_small$LP.num == "p2"]
    long_id_p1_expr <- model_data_gene %>% filter(MAPID %in% long_ids, LP.num == "p1") %>% select("MAPID", "expression") %>% unique()
    row.names(long_id_p1_expr) <- long_id_p1_expr$MAPID
    
    
    for(i in 1:nrow(model_data_gene)) {
      id <- model_data_gene[i, "MAPID"]
      if (id %in% long_id_p1_expr$MAPID) {
        model_data_gene[i,"expression"] <- long_id_p1_expr[id,"expression"]
      }
    }
    
    # Survival for CDR by time since LP
    surv_summ <- summary(coxph(Surv(time_elapsed, cdr >= 0.5) ~ expression, data = model_data_gene))
    if (as.numeric(surv_summ$sctest["pvalue"]) < 0.05) {
      pdf(paste0("/home/eteleeb/projects/AD_Staging/Figures/", features$EntrezGeneSymbol[features$SOMAseqID ==gene], "_", gene, "_surv.pdf"), 
          width = 5, height = 5, useDingbats = F)
      surviplot(Surv(time_elapsed, cdr >= 0.5) ~ expression, data = model_data_gene, 
                main = paste0(features$EntrezGeneSymbol[features$SOMAseqID ==gene], "/", features$Target[features$SOMAseqID ==gene]),
                ylab='Survival Proportion',xlab='Time since LP', cex.main=1, mark.time=F, col=c("skyblue3", '#E7298A'), 
                lwd = 4, cex.axis=1.2, cex.lab= 1.5, hr.pos="left")
      #legend('bottomleft', legend=c('A', 'B'), col=c("skyblue3", '#E7298A'), pch=16, cex=1)
      dev.off()
    }
    return(data.frame(protein = gene, HR = surv_summ$coefficients[,"exp(coef)"], pval = surv_summ$sctest["pvalue"]))
  }
  else return(data.frame(protein = gene, HR = NA, pval = NA))
}

genes.of.int <- data.frame(gene.name =c('IGF1','NRXN3','YWHAZ'), gene.id = c('X2952.75', 'X5111.15', 'X5858.6'))
           
for (g in genes.of.int$gene.id) {
  run_survival(g)
}


### write results 
surv_all <- ldply(colnames(log_readings), run_survival)
surv_all$padj <- p.adjust(surv_all$pval, method = "BH")
surv_all <- surv_all[order(surv_all$padj),]

colnames(surv_all)[1] <-"SOMAseqID"
surv_all_joined <- inner_join(surv_all, features[,c(1,4,5,8)])

write.csv(surv_all_joined,"~/multi-omics/data/survival_allproteins_allcdr0_and_onecdr1_event0.5.csv", row.names = F)

#surviplot(Surv(time) ~ group, data=clin[clin$group !="Control", ], ylab='Survival Proportion',xlim=c(55,max.time), main ="", 
#          xlab='Age at onset (years)', cex.main=1, mark.time=TRUE, col=c("skyblue3", '#E7298A'), lwd = 3, cex.axis=2, cex.lab= 1.5)
