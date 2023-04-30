
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

Data_wrangling_ATAC = function(option_list)
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
  
  #### READ and transform Open_in_CT_threshold ----
  
  Open_in_CT_threshold = opt$Open_in_CT_threshold
  
  cat("Open_in_CT_threshold_\n")
  cat(str(Open_in_CT_threshold))
  cat("\n")
  
  #### READ and transform tracking_variants ----
  
  tracking_variants = unlist(strsplit(opt$tracking_variants, split=","))
  
  cat("tracking_variants_\n")
  cat(sprintf(as.character(tracking_variants)))
  cat("\n")
  
  #### READ and transform desiR_weights ----
  
  desiR_weights = as.numeric(unlist(strsplit(opt$desiR_weights, split=",")))
  
  cat("desiR_weights_\n")
  cat(sprintf(as.character(desiR_weights)))
  cat("\n")
  cat(str(desiR_weights))
  cat("\n")
  
  #### Read ATAC_pre_ranked file ----
  
  ATAC_pre_ranked<-as.data.frame(fread(file=opt$ATAC_pre_ranked,sep="\t") , stringsAsFactors=F)
  
  cat("ATAC_pre_ranked_0\n")
  cat(str(ATAC_pre_ranked))
  cat("\n")
  cat(str(unique(ATAC_pre_ranked$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(ATAC_pre_ranked$Lineage))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(ATAC_pre_ranked$Lineage)))))
  cat("\n")
  
  #### Converge lineage lymph ----
  
  ATAC_pre_ranked$Lineage[which(ATAC_pre_ranked$Lineage == 'lymph_lineage.CD4')]<-"lymph_lineage"
  ATAC_pre_ranked$Lineage[which(ATAC_pre_ranked$Lineage == 'lymph_lineage.CD8')]<-"lymph_lineage"
  ATAC_pre_ranked$Lineage[which(ATAC_pre_ranked$Lineage == 'lymph_lineage.B')]<-"lymph_lineage"
  ATAC_pre_ranked$Lineage[which(ATAC_pre_ranked$Lineage == 'lymph_lineage.NK')]<-"lymph_lineage"
  
  ATAC_pre_ranked_unified_lymph<-unique(ATAC_pre_ranked)
  
  cat("ATAC_pre_ranked_unified_lymph_0\n")
  cat(str(ATAC_pre_ranked_unified_lymph))
  cat("\n")
  cat(str(unique(ATAC_pre_ranked_unified_lymph$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(ATAC_pre_ranked_unified_lymph$Lineage))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(ATAC_pre_ranked_unified_lymph$Lineage)))))
  cat("\n")
  
  
  # quit(status = 1)
  
  #### Threshold per cell type -----
  
  ATAC_pre_ranked_unified_lymph_Thresholded<-ATAC_pre_ranked_unified_lymph[which(ATAC_pre_ranked_unified_lymph$value >= Open_in_CT_threshold),]
  
  cat("ATAC_pre_ranked_unified_lymph_Thresholded_0\n")
  cat(str(ATAC_pre_ranked_unified_lymph_Thresholded))
  cat("\n")
  cat(str(unique(ATAC_pre_ranked_unified_lymph_Thresholded$VAR)))
  cat("\n")
  
  
  #### how many CT open per lineage per variant ----
  
  ATAC_pre_ranked_unified_lymph_Thresholded.dt<-data.table(ATAC_pre_ranked_unified_lymph_Thresholded, key=c("VAR","Lineage"))
  

  Freq_table<-as.data.frame(ATAC_pre_ranked_unified_lymph_Thresholded.dt[,.(nCells=.N), by=key(ATAC_pre_ranked_unified_lymph_Thresholded.dt)], stringsAsFactors=F)
  
  cat("Freq_table_0\n")
  cat(str(Freq_table))
  cat("\n")
  cat(str(unique(Freq_table$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Freq_table$nCells)))))
  cat("\n")
  cat(sprintf(as.character(summary(Freq_table$nCells))))
  cat("\n")
  
  #### max signal per lineage per variant ----
  
  ATAC_pre_ranked_unified_lymph.dt<-data.table(ATAC_pre_ranked_unified_lymph, key=c("VAR","Lineage"))
  
  
  MAX_table<-as.data.frame(ATAC_pre_ranked_unified_lymph.dt[,.SD[which.max(value)], by=key(ATAC_pre_ranked_unified_lymph.dt)], stringsAsFactors=F)
  
  cat("MAX_table_0\n")
  cat(str(MAX_table))
  cat("\n")
  cat(str(unique(MAX_table$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(MAX_table$value)))))
  cat("\n")
  cat(sprintf(as.character(summary(MAX_table$value))))
  cat("\n")
  
  #### merge tables ----
  
  
  Merge_table<-merge(Freq_table,
                     MAX_table,
                     by=c("VAR","Lineage"),
                     all=T)
  
  Merge_table$nCells[is.na(Merge_table$nCells)]<-0
  
  cat("Merge_table_0\n")
  cat(str(Merge_table))
  cat("\n")
  cat(str(unique(Merge_table$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Merge_table$nCells)))))
  cat("\n")
  cat(sprintf(as.character(summary(Merge_table$nCells))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Merge_table$value)))))
  cat("\n")
  cat(sprintf(as.character(summary(Merge_table$value))))
  cat("\n")
  
  # ################################################
  # quit(status=1)
  
  ################################## desiR nCells ---------------------------------------------
  
  nCells_LOW<-desiR_weights[1]
  nCells_HIGH<-desiR_weights[2]
  
  cat("desiR_weights_nCells\n")
  cat(sprintf(as.character(c(nCells_LOW,nCells_HIGH))))
  cat("\n")
  
  graph_path<-paste(out,'desiR_graphs','/',sep='')
  setwd(graph_path)
  
  pdf(file="desIR_ATAC_nCells.pdf")
  hist(Merge_table$nCells, breaks=50, col="grey", border="white", main="",
       xlab="ATAC-Seq in a Relevant Cell Type")
  des.line(Merge_table$nCells, "d.high", des.args=c(cut1=nCells_LOW, cut2=nCells_HIGH, scale=0.5))
  dev.off()
  
  Merge_table$nCells_weight <- d.high(Merge_table$nCells, cut1=nCells_LOW, cut2=nCells_HIGH, scale=0.5)
  
  
  ################################## desiR value ---------------------------------------------
  
  ATAC_value_LOW<-desiR_weights[3]
  ATAC_value_HIGH<-desiR_weights[4]
  
  cat("desiR_weights_ATAC_value\n")
  cat(sprintf(as.character(c(ATAC_value_LOW,ATAC_value_HIGH))))
  cat("\n")
  
  
  setwd(graph_path)
  
  pdf(file="desIR_ATAC_value.pdf")
  hist(Merge_table$value, breaks=50, col="grey", border="white", main="",
       xlab="ATAC-Seq in a Relevant Cell Type")
  des.line(Merge_table$value, "d.high", des.args=c(cut1=ATAC_value_LOW, cut2=ATAC_value_HIGH, scale=0.5))
  dev.off()
  
  Merge_table$ATAC_value_weight <- d.high(Merge_table$value, cut1=ATAC_value_LOW, cut2=ATAC_value_HIGH, scale=0.5)
  
  
  #### Overall desirability ---- 
  
  Overall_nCells<-desiR_weights[5]
  Overall_ATAC_value<-desiR_weights[6]
  
  cat("desiR_weights_Overall\n")
  cat(sprintf(as.character(c(Overall_nCells,Overall_ATAC_value))))
  cat("\n")
  
  
  
  Merge_table$Overall_weight <- d.overall(Merge_table$nCells_weight, Merge_table$ATAC_value_weight, 
                                          weights=c(Overall_nCells,Overall_ATAC_value))
  
  setwd(graph_path)
  
  par(las=1)
  pdf(file= "desIR_ATAC_overall.pdf")
  plot(rev(sort(Merge_table$Overall_weight)), type="l", xlab="Rank", ylab="Overall Desirability")
  dev.off()
  
  
  check<-Merge_table[which(Merge_table$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(check$VAR)))
  cat("\n")
  cat(sprintf(as.character(check$nCells)))
  cat("\n")
  cat(sprintf(as.character(check$value)))
  cat("\n")

  
  # #### SAVE ----
  
  setwd(out)
  
  write.table(Merge_table, file="ATAC_GLOBAL_Ranked.tsv", sep="\t", quote = F, row.names = F)
  
  write.table(check, file="check_ATAC_GLOBAL_Ranked.tsv", sep="\t", quote = F, row.names = F)
  
  
  
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
  
  
  #### Read ATAC_ranked file ----
  
  setwd(out)
  
  ATAC_ranked<-as.data.frame(fread(file="ATAC_GLOBAL_Ranked.tsv",sep="\t") , stringsAsFactors=F)
  
  cat("ATAC_ranked_0\n")
  cat(str(ATAC_ranked))
  cat("\n")
  cat(str(unique(ATAC_ranked$VAR)))
  cat("\n")
  
  #### calculate mean and sd per lineage -----
  
  ATAC_ranked.dt<-data.table(ATAC_ranked, key=c("Lineage"))
  
  
  ATAC_ranked_Lineage_parameters<-as.data.frame(ATAC_ranked.dt[,.(mean_Overall_weight=mean(Overall_weight, na.rm =T),
                                                                  sd_Overall_weight=sd(Overall_weight, na.rm =T)),
                                                               by=key(ATAC_ranked.dt)], stringsAsFactors=F)
  
  cat("ATAC_ranked_Lineage_parameters_0\n")
  cat(str(ATAC_ranked_Lineage_parameters))
  cat("\n")
  cat(str(unique(ATAC_ranked_Lineage_parameters$VAR)))
  cat("\n")
  
  #### Merge and calculate Z-score -----
  
  ATAC_ranked<-merge(ATAC_ranked,
                     ATAC_ranked_Lineage_parameters,
                     by=c("Lineage"))
  
  cat("ATAC_ranked_POST_merge\n")
  cat(str(ATAC_ranked))
  cat("\n")
  cat(str(unique(ATAC_ranked$VAR)))
  cat("\n")
  
  
  ATAC_ranked$Overall_weight_Z_score<-(ATAC_ranked$Overall_weight-ATAC_ranked$mean_Overall_weight)/ATAC_ranked$sd_Overall_weight
  
  cat("ATAC_ranked_POST_merge_Z_score\n")
  cat(str(ATAC_ranked))
  cat("\n")
  cat(str(unique(ATAC_ranked$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ATAC_ranked$Overall_weight_Z_score)))))
  cat("\n")
  cat(sprintf(as.character(summary(ATAC_ranked$Overall_weight_Z_score))))
  cat("\n")
  
  #### select columns and change names to converge ----
  
  indx.int<-c(which(colnames(ATAC_ranked) == "VAR"),which(colnames(ATAC_ranked) == "Lineage"),which(colnames(ATAC_ranked) == "Overall_weight"),which(colnames(ATAC_ranked) == "Overall_weight_Z_score"))
  
  ATAC_ranked_subset<-unique(ATAC_ranked[,indx.int])
  
  cat("ATAC_ranked_subset_0\n")
  cat(str(ATAC_ranked_subset))
  cat("\n")
  cat(str(unique(ATAC_ranked_subset$VAR)))
  cat("\n")
  
  colnames(ATAC_ranked_subset)[which(colnames(ATAC_ranked_subset) == "Lineage")]<-"variable"
  colnames(ATAC_ranked_subset)[which(colnames(ATAC_ranked_subset) == "Overall_weight")]<-"value"
  colnames(ATAC_ranked_subset)[which(colnames(ATAC_ranked_subset) == "Overall_weight_Z_score")]<-"value_Z_score"
  
  ATAC_ranked_subset$variable<-paste('Rank_ATAC_',ATAC_ranked_subset$variable, sep='')
  
  
  cat("ATAC_ranked_subset_1\n")
  cat(str(ATAC_ranked_subset))
  cat("\n")
  cat(str(unique(ATAC_ranked_subset$VAR)))
  cat("\n")
  
  check<-ATAC_ranked_subset[which(ATAC_ranked_subset$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  # #### SAVE ----
  
  setwd(out)
  
  saveRDS(ATAC_ranked_subset, file="Prepared_file_ATAC.rds")
  
  write.table(check, file="check_ATAC_GLOBAL_Ranked.tsv", sep="\t", quote = F, row.names = F)
  
  
  
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
    make_option(c("--ATAC_pre_ranked"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--desiR_weights"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Open_in_CT_threshold"), type="numeric", default=NULL, 
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
  
  Data_wrangling_ATAC(opt)
  Z_score_normalization_and_convergence(opt)
  
}


###########################################################################

system.time( main() )
