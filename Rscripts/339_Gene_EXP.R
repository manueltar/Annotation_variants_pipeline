
suppressMessages(library("plyr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("data.table", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("crayon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("farver", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("labeling", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("optparse", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("dplyr", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("backports", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("broom", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rstudioapi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cli", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tzdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tidyverse", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
suppressMessages(library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggeasy", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("sandwich", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("digest", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("desiR", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/")



opt = NULL

options(warn = 1)

Data_wrangling_GENE_EXP = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform tracking_variants ----
  
  tracking_variants = unlist(strsplit(opt$tracking_variants, split=","))
  
  cat("tracking_variants_\n")
  cat(sprintf(as.character(tracking_variants)))
  cat("\n")
  
  #### READ and transform relevant_not_relevant_weights ----
  
  relevant_not_relevant_weights = as.numeric(unlist(strsplit(opt$relevant_not_relevant_weights, split=",")))
  
  cat("relevant_not_relevant_weights_\n")
  cat(sprintf(as.character(relevant_not_relevant_weights)))
  cat("\n")
  cat(str(relevant_not_relevant_weights))
  cat("\n")
  
  #### READ and transform desiR_weights ----
  
  desiR_weights = as.numeric(unlist(strsplit(opt$desiR_weights, split=",")))
  
  cat("desiR_weights_\n")
  cat(sprintf(as.character(desiR_weights)))
  cat("\n")
  cat(str(desiR_weights))
  cat("\n")
  
  #### Read GENE_EXP_pre_ranked file ----
  
  GENE_EXP_pre_ranked<-as.data.frame(fread(file=opt$GENE_EXP_pre_ranked,sep="\t") , stringsAsFactors=F)
  
  cat("GENE_EXP_pre_ranked_0\n")
  cat(str(GENE_EXP_pre_ranked))
  cat("\n")
  cat(str(unique(GENE_EXP_pre_ranked$VAR)))
  cat("\n")
  cat(str(unique(GENE_EXP_pre_ranked$ensembl_gene_id)))
  cat("\n")
  cat(str(unique(GENE_EXP_pre_ranked$phenotype)))
  cat("\n")
  
  ### melt GENE EXP and Tag separately and GENE_EXP_pre_ranked.m later ----
  
  GENE_EXP_pre_ranked.m<-melt(GENE_EXP_pre_ranked, id.vars=c("VAR","ensembl_gene_id","HGNC","phenotype"), value.name="value", variable.name="variable")
  
  GENE_EXP_pre_ranked.m$variable<-as.character(GENE_EXP_pre_ranked.m$variable)
  
  cat("GENE_EXP_pre_ranked.m_0\n")
  cat(str(GENE_EXP_pre_ranked.m))
  cat("\n")
  cat(str(unique(GENE_EXP_pre_ranked.m$VAR)))
  cat("\n")
  cat(str(unique(GENE_EXP_pre_ranked.m$ensembl_gene_id)))
  cat("\n")
  cat(str(unique(GENE_EXP_pre_ranked.m$phenotype)))
  cat("\n")
  
  GENE_EXP_pre_ranked.m$Cell_Type<-GENE_EXP_pre_ranked.m$variable
  
  GENE_EXP_pre_ranked.m$Cell_Type<-gsub("mean_GENE_EXP_|Tag_","", GENE_EXP_pre_ranked.m$Cell_Type)
  
  
  GENE_EXP_pre_ranked.m$variable[grep("mean_GENE_EXP_", GENE_EXP_pre_ranked.m$variable)]<-"mean_GENE_EXP"
  GENE_EXP_pre_ranked.m$variable[grep("Tag_", GENE_EXP_pre_ranked.m$variable)]<-"Tag"
  
  
  cat("GENE_EXP_pre_ranked.m_1\n")
  cat(str(GENE_EXP_pre_ranked.m))
  cat("\n")
  cat(str(unique(GENE_EXP_pre_ranked.m$VAR)))
  cat("\n")
  cat(str(unique(GENE_EXP_pre_ranked.m$ensembl_gene_id)))
  cat("\n")
  cat(str(unique(GENE_EXP_pre_ranked.m$phenotype)))
  cat("\n")
  
  GENE_EXP_pre_ranked.m_wide<-as.data.frame(pivot_wider(GENE_EXP_pre_ranked.m,
                                                        id_cols=c("VAR","HGNC","ensembl_gene_id","phenotype","Cell_Type"),
                                                        names_from=variable,
                                                        values_from=value), stringsAsFactors=F)
  
  GENE_EXP_pre_ranked.m_wide$mean_GENE_EXP<-as.numeric(GENE_EXP_pre_ranked.m_wide$mean_GENE_EXP)
  
  
  cat("GENE_EXP_pre_ranked.m_wide_1\n")
  cat(str(GENE_EXP_pre_ranked.m_wide))
  cat("\n")
  cat(str(unique(GENE_EXP_pre_ranked.m_wide$VAR)))
  cat("\n")
  cat(str(unique(GENE_EXP_pre_ranked.m_wide$ensembl_gene_id)))
  cat("\n")
  cat(str(unique(GENE_EXP_pre_ranked.m_wide$phenotype)))
  cat("\n")
  
  
  check_1<-GENE_EXP_pre_ranked.m_wide[which(GENE_EXP_pre_ranked.m_wide$VAR%in%tracking_variants),]
  
  
  cat("check_1_0\n")
  cat(str(check_1))
  cat("\n")
  cat(sprintf(as.character(check_1$VAR)))
  cat("\n")
  cat(sprintf(as.character(check_1$nCells)))
  cat("\n")
 
 
  # #######################################################
  # quit(status = 1)
  
  
  #### Aggregate GENE EXP by VAR, ENSG and Tag, forget phenotype  ----
  
  GENE_EXP_pre_ranked.m_wide.dt<-data.table(GENE_EXP_pre_ranked.m_wide, key=c("VAR","HGNC","ensembl_gene_id","Tag"))
  

  Aggregation_table<-as.data.frame(GENE_EXP_pre_ranked.m_wide.dt[,.(Aggregate_GENE_EXP=sum(mean_GENE_EXP),
                                                          nCells=.N), by=key(GENE_EXP_pre_ranked.m_wide.dt)], stringsAsFactors=F)
  
  cat("Aggregation_table_0\n")
  cat(str(Aggregation_table))
  cat("\n")
  cat(str(unique(Aggregation_table$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Aggregation_table$Aggregate_GENE_EXP)))))
  cat("\n")
  cat(sprintf(as.character(summary(Aggregation_table$Aggregate_GENE_EXP))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Aggregation_table$nCells)))))
  cat("\n")
  cat(sprintf(as.character(summary(Aggregation_table$nCells))))
  cat("\n")
  
  
  
  #### Multiply by relevant/ not relevant weights ----
  
  Aggregation_table$Multiplier<-NA
  
  Aggregation_table$Multiplier[which(Aggregation_table$Tag == "Relevant")]<-relevant_not_relevant_weights[1]
  Aggregation_table$Multiplier[which(Aggregation_table$Tag == "Not_relevant")]<-relevant_not_relevant_weights[2]
  
  
 
  
  Aggregation_table$Aggregate_GENE_EXP_multiplied<-Aggregation_table$Aggregate_GENE_EXP*Aggregation_table$Multiplier
 
  
  cat("Aggregation_table_1\n")
  cat(str(Aggregation_table))
  cat("\n")
  cat(str(unique(Aggregation_table$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Aggregation_table$Aggregate_GENE_EXP_multiplied)))))
  cat("\n")
  cat(sprintf(as.character(summary(Aggregation_table$Aggregate_GENE_EXP_multiplied))))
  cat("\n")

  #### Divide by the number of relevant and not relevant cells ----
  
  Aggregation_table$Aggregate_GENE_EXP_normalised<-Aggregation_table$Aggregate_GENE_EXP_multiplied/Aggregation_table$nCells
  
  cat("Aggregation_table_2\n")
  cat(str(Aggregation_table))
  cat("\n")
  cat(str(unique(Aggregation_table$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Aggregation_table$Aggregate_GENE_EXP_normalised)))))
  cat("\n")
  cat(sprintf(as.character(summary(Aggregation_table$Aggregate_GENE_EXP_normalised))))
  cat("\n")
  
  check_2<-Aggregation_table[which(Aggregation_table$VAR%in%tracking_variants),]
  
  
  cat("check_2_0\n")
  cat(str(check_2))
  cat("\n")
  cat(sprintf(as.character(check_2$VAR)))
  cat("\n")
  cat(sprintf(as.character(check_2$nCells)))
  cat("\n")
 
  
  
  #### Finally add up per variant and gene  ----
  
  Aggregation_table.dt<-data.table(Aggregation_table, key=c("VAR","ensembl_gene_id","HGNC"))
  
  
  Aggregation_table_FINAL<-as.data.frame(Aggregation_table.dt[,.(Aggregate_GENE_EXP_FINAL=sum(Aggregate_GENE_EXP_normalised)), by=key(Aggregation_table.dt)], stringsAsFactors=F)
  
  cat("Aggregation_table_FINAL_0\n")
  cat(str(Aggregation_table_FINAL))
  cat("\n")
  cat(sprintf(as.character(names(summary(Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL)))))
  cat("\n")
  cat(sprintf(as.character(summary(Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL))))
  cat("\n")
  
  
  
  check_3<-Aggregation_table_FINAL[which(Aggregation_table_FINAL$VAR%in%tracking_variants),]
  
  
  cat("check_3_0\n")
  cat(str(check_3))
  cat("\n")
  cat(sprintf(as.character(check_3$VAR)))
  cat("\n")
  cat(sprintf(as.character(check_3$nCells)))
  cat("\n")

  
  # #### check ----
  
  setwd(out)
  

  write.table(check_1, file="check_1.tsv", sep="\t", quote = F, row.names = F)
  write.table(check_2, file="check_2.tsv", sep="\t", quote = F, row.names = F)
  write.table(check_3, file="check_3.tsv", sep="\t", quote = F, row.names = F)
  
  
  # ################################################################
  # quit(status = 1)
  
  
  
  ################################## desiR GENE_EXP ---------------------------------------------
  
  graph_path<-paste(out,'desiR_graphs','/',sep='')
  setwd(graph_path)
  
  pdf(file="desiR_Aggregate_GENE_EXP_FINAL_PRE.pdf")
  hist(Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL, breaks=50, col="grey", border="white", main="",
       xlab="Aggregate_GENE_EXP_FINAL")
  dev.off()
  
  
  
  
  # ######################################################################
  # quit(status = 1)
  
  Aggregate_GENE_EXP_FINAL_LOW<-desiR_weights[1]
  Aggregate_GENE_EXP_FINAL_HIGH<-desiR_weights[2]
  
  cat("desiR_weights_GENE_EXP\n")
  cat(sprintf(as.character(c(Aggregate_GENE_EXP_FINAL_LOW,Aggregate_GENE_EXP_FINAL_HIGH))))
  cat("\n")
 
    
  pdf(file="desIR_Aggregate_GENE_EXP_FINAL.pdf")
  hist(Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL, breaks=50, col="grey", border="white", main="",
       xlab="GENE_EXP-Seq in a Relevant Cell Type")
  des.line(Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL, "d.high", des.args=c(cut1=Aggregate_GENE_EXP_FINAL_LOW, cut2=Aggregate_GENE_EXP_FINAL_HIGH, scale=0.5))
  dev.off()
  
  Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL_component <- d.high(Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL, cut1=Aggregate_GENE_EXP_FINAL_LOW, cut2=Aggregate_GENE_EXP_FINAL_HIGH, scale=0.5)
  
  cat("Aggregation_table_FINAL_1\n")
  cat(str(Aggregation_table_FINAL))
  cat("\n")
  cat(sprintf(as.character(names(summary(Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL_component)))))
  cat("\n")
  cat(sprintf(as.character(summary(Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL_component))))
  cat("\n")
  
  #### Overall desirability ---- 
  
  
  Overall_Aggregate_GENE_EXP_FINAL<-desiR_weights[3]
  
  cat("desiR_weights_Overall\n")
  cat(sprintf(as.character(c(Overall_Aggregate_GENE_EXP_FINAL))))
  cat("\n")
  
  
  
  Aggregation_table_FINAL$Overall_weight <- d.overall(Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL_component, 
                                          weights=c(Overall_Aggregate_GENE_EXP_FINAL))
  
  setwd(graph_path)
  
  par(las=1)
  pdf(file= "desIR_GENE_EXP_overall.pdf")
  plot(rev(sort(Aggregation_table_FINAL$Overall_weight)), type="l", xlab="Rank", ylab="Overall Desirability")
  dev.off()
  
  cat("Aggregation_table_FINAL_2\n")
  cat(str(Aggregation_table_FINAL))
  cat("\n")
  cat(sprintf(as.character(names(summary(Aggregation_table_FINAL$Overall_weight)))))
  cat("\n")
  cat(sprintf(as.character(summary(Aggregation_table_FINAL$Overall_weight))))
  cat("\n")
  
  
 
  
  
  check<-Aggregation_table_FINAL[which(Aggregation_table_FINAL$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(check$VAR)))
  cat("\n")
  cat(sprintf(as.character(check$nCells)))
  cat("\n")
  cat(sprintf(as.character(check$Aggregate_GENE_EXP_FINAL)))
  cat("\n")

  
  # #### SAVE ----
  
  setwd(out)
  
  write.table(Aggregation_table_FINAL, file="GENE_EXP_GLOBAL_Ranked.tsv", sep="\t", quote = F, row.names = F)
  
  write.table(check, file="check_GENE_EXP_GLOBAL_Ranked.tsv", sep="\t", quote = F, row.names = F)
  
  # ################################################################
  # quit(status = 1)
  
}

Z_score_normalization_and_convergence = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform tracking_variants ----
  
  tracking_variants = unlist(strsplit(opt$tracking_variants, split=","))
  
  cat("tracking_variants_\n")
  cat(sprintf(as.character(tracking_variants)))
  cat("\n")
  
  
  #### Read GENE_EXP_ranked file ----
  
  setwd(out)
  
  GENE_EXP_ranked<-as.data.frame(fread(file="GENE_EXP_GLOBAL_Ranked.tsv",sep="\t") , stringsAsFactors=F)
  
  cat("GENE_EXP_ranked_0\n")
  cat(str(GENE_EXP_ranked))
  cat("\n")
  cat(str(unique(GENE_EXP_ranked$VAR)))
  cat("\n")
  
 
  #### calculate mean and sd -----
  
  GENE_EXP_ranked$mean_Overall_weight<-mean(GENE_EXP_ranked$Overall_weight, na.rm =T)
  GENE_EXP_ranked$sd_Overall_weight<-sd(GENE_EXP_ranked$Overall_weight, na.rm =T)
  
  
  cat("GENE_EXP_ranked_POST_GENE_EXP_pre_ranked.m\n")
  cat(str(GENE_EXP_ranked))
  cat("\n")
  cat(str(unique(GENE_EXP_ranked$VAR)))
  cat("\n")
  
  
  GENE_EXP_ranked$Overall_weight_Z_score<-(GENE_EXP_ranked$Overall_weight-GENE_EXP_ranked$mean_Overall_weight)/GENE_EXP_ranked$sd_Overall_weight
  
  
  GENE_EXP_ranked$Overall_weight_Z_score<-(GENE_EXP_ranked$Overall_weight-GENE_EXP_ranked$mean_Overall_weight)/GENE_EXP_ranked$sd_Overall_weight
  
  cat("GENE_EXP_ranked_POST_GENE_EXP_pre_ranked.m_Z_score\n")
  cat(str(GENE_EXP_ranked))
  cat("\n")
  cat(str(unique(GENE_EXP_ranked$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(GENE_EXP_ranked$Overall_weight_Z_score)))))
  cat("\n")
  cat(sprintf(as.character(summary(GENE_EXP_ranked$Overall_weight_Z_score))))
  cat("\n")
  
  #### select columns and change names to converge ----
  
  indx.int<-c(which(colnames(GENE_EXP_ranked) == "VAR"),which(colnames(GENE_EXP_ranked) == "ensembl_gene_id"),
              which(colnames(GENE_EXP_ranked) == "HGNC"),
              which(colnames(GENE_EXP_ranked) == "Overall_weight"),which(colnames(GENE_EXP_ranked) == "Overall_weight_Z_score"))
  
  GENE_EXP_ranked_subset<-unique(GENE_EXP_ranked[,indx.int])
  
  cat("GENE_EXP_ranked_subset_0\n")
  cat(str(GENE_EXP_ranked_subset))
  cat("\n")
  cat(str(unique(GENE_EXP_ranked_subset$VAR)))
  cat("\n")
  
 
  colnames(GENE_EXP_ranked_subset)[which(colnames(GENE_EXP_ranked_subset) == "Overall_weight")]<-"value"
  colnames(GENE_EXP_ranked_subset)[which(colnames(GENE_EXP_ranked_subset) == "Overall_weight_Z_score")]<-"value_Z_score"
  
  GENE_EXP_ranked_subset$variable<-paste('Rank_GENE_EXP',sep='')
  
  
  cat("GENE_EXP_ranked_subset_1\n")
  cat(str(GENE_EXP_ranked_subset))
  cat("\n")
  cat(str(unique(GENE_EXP_ranked_subset$VAR)))
  cat("\n")
  
  check<-GENE_EXP_ranked_subset[which(GENE_EXP_ranked_subset$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  # #### SAVE ----
  
  setwd(out)
  
  saveRDS(GENE_EXP_ranked_subset, file="Prepared_file_GENE_EXP.rds")
  
  write.table(check, file="check_GENE_EXP_GLOBAL_Ranked.tsv", sep="\t", quote = F, row.names = F)
  
  
  
}

printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--GENE_EXP_pre_ranked"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--desiR_weights"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--relevant_not_relevant_weights"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--tracking_variants"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
       make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
        make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  Data_wrangling_GENE_EXP(opt)
  Z_score_normalization_and_convergence(opt)
  
}


###########################################################################

system.time( main() )
