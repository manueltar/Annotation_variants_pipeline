
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
  
  #### Read multi_ATAC_ranked file ----
  
  multi_ATAC_ranked<-as.data.frame(fread(file=opt$multi_ATAC_ranked,sep="\t") , stringsAsFactors=F)
  
  cat("multi_ATAC_ranked_0\n")
  cat(str(multi_ATAC_ranked))
  cat("\n")
  cat(str(unique(multi_ATAC_ranked$VAR)))
  cat("\n")
  
  #### nLineages open ----
  
  multi_ATAC_ranked.dt<-data.table(multi_ATAC_ranked, key="VAR")
  
  Freq_table<-as.data.frame(multi_ATAC_ranked.dt[,.(nLineages=.N), by=key(multi_ATAC_ranked.dt)], stringsAsFactors=F)
  
  cat("Freq_table_0\n")
  cat(str(Freq_table))
  cat("\n")
  cat(str(unique(Freq_table$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Freq_table$nLineages)))))
  cat("\n")
  cat(sprintf(as.character(summary(Freq_table$nLineages))))
  cat("\n")
  
  #### mean ATAC ----
  
  multi_ATAC_ranked.dt<-data.table(multi_ATAC_ranked, key="VAR")
  
  mean_table<-as.data.frame(multi_ATAC_ranked.dt[,.(mean_Rank_ATAC=mean(Overall_weight)), by=key(multi_ATAC_ranked.dt)], stringsAsFactors=F)
  
  cat("mean_table_0\n")
  cat(str(mean_table))
  cat("\n")
  cat(str(unique(mean_table$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(mean_table$mean_Rank_ATAC)))))
  cat("\n")
  cat(sprintf(as.character(summary(mean_table$mean_Rank_ATAC))))
  cat("\n")
  
  #### merge tables ----
  
  
  Merge_table<-merge(Freq_table,
                     mean_table,
                     by="VAR",
                     all=T)
  
  
  cat("Merge_table_0\n")
  cat(str(Merge_table))
  cat("\n")
  cat(str(unique(Merge_table$VAR)))
  cat("\n")
  
  
  # quit(status = 1)
  
  
  ################################## desiR nLineages ---------------------------------------------
  
  nLineages_LOW<-desiR_weights[1]
  nLineages_HIGH<-desiR_weights[2]
  
  cat("desiR_weights_nLineages\n")
  cat(sprintf(as.character(c(nLineages_LOW,nLineages_HIGH))))
  cat("\n")
  
  
  graph_path<-paste(out,'desiR_graphs','/',sep='')
  setwd(graph_path)
  
  pdf(file="desIR_multi_ATAC_nLineages.pdf")
  hist(Merge_table$nLineages, breaks=50, col="grey", border="white", main="",
       xlab="ATAC-Seq in a Relevant Cell Type")
  des.line(Merge_table$nLineages, "d.high", des.args=c(cut1=nLineages_LOW, cut2=nLineages_HIGH, scale=0.5))
  dev.off()
  
  Merge_table$nLineages_weight <- d.high(Merge_table$nLineages, cut1=nLineages_LOW, cut2=nLineages_HIGH, scale=0.5)
  
  
  ################################## desiR mean_Rank_ATAC ---------------------------------------------
  
  mean_Rank_ATAC_LOW<-desiR_weights[3]
  mean_Rank_ATAC_HIGH<-desiR_weights[4]
  
  cat("desiR_weights_mean_Rank_ATAC\n")
  cat(sprintf(as.character(c(mean_Rank_ATAC_LOW,mean_Rank_ATAC_HIGH))))
  cat("\n")
  
  
  setwd(graph_path)
  
  pdf(file="desIR_multi_ATAC_mean_Rank_ATAC.pdf")
  hist(Merge_table$mean_Rank_ATAC, breaks=50, col="grey", border="white", main="",
       xlab="ATAC-Seq in a Relevant Cell Type")
  des.line(Merge_table$mean_Rank_ATAC, "d.high", des.args=c(cut1=mean_Rank_ATAC_LOW, cut2=mean_Rank_ATAC_HIGH, scale=0.5))
  dev.off()
  
  Merge_table$mean_Rank_ATAC_weight <- d.high(Merge_table$mean_Rank_ATAC, cut1=mean_Rank_ATAC_LOW, cut2=mean_Rank_ATAC_HIGH, scale=0.5)
  
  
  #### Overall desirability ---- 
  
  Overall_nLineages<-desiR_weights[5]
  Overall_mean_Rank_ATAC<-desiR_weights[6]
  
  cat("desiR_weights_Overall\n")
  cat(sprintf(as.character(c(Overall_nLineages,Overall_mean_Rank_ATAC))))
  cat("\n")
  
  
  
  Merge_table$Overall_weight <- d.overall(Merge_table$nLineages_weight, Merge_table$mean_Rank_ATAC_weight, 
                                          weights=c(Overall_nLineages,Overall_mean_Rank_ATAC))
  
  setwd(graph_path)
  
  par(las=1)
  pdf(file= "desIR_overall_multi_ATAC.pdf")
  plot(rev(sort(Merge_table$Overall_weight)), type="l", xlab="Rank", ylab="Overall Desirability")
  dev.off()
  
  
  check<-Merge_table[which(Merge_table$VAR%in%tracking_variants),]

  
  # #### SAVE ----
  
  setwd(out)
  
  write.table(Merge_table, file="multi_ATAC_GLOBAL_Ranked.tsv", sep="\t", quote = F, row.names = F)
  
  write.table(check, file="check_multi_ATAC_GLOBAL_Ranked.tsv", sep="\t", quote = F, row.names = F)
  
  
  
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
  
  
  #### Read multi_ATAC_ranked file ----
  
  setwd(out)
  
  multi_ATAC_ranked<-as.data.frame(fread(file="multi_ATAC_GLOBAL_Ranked.tsv",sep="\t") , stringsAsFactors=F)
  
  cat("multi_ATAC_ranked_0\n")
  cat(str(multi_ATAC_ranked))
  cat("\n")
  cat(str(unique(multi_ATAC_ranked$VAR)))
  cat("\n")
  
  #### calculate mean and sd per lineage -----
 
  multi_ATAC_ranked$mean_Overall_weight<-mean(multi_ATAC_ranked$Overall_weight, na.rm =T)
  multi_ATAC_ranked$sd_Overall_weight<-sd(multi_ATAC_ranked$Overall_weight, na.rm =T)
  
  
  cat("multi_ATAC_ranked_POST_merge\n")
  cat(str(multi_ATAC_ranked))
  cat("\n")
  cat(str(unique(multi_ATAC_ranked$VAR)))
  cat("\n")
  
  
  multi_ATAC_ranked$Overall_weight_Z_score<-(multi_ATAC_ranked$Overall_weight-multi_ATAC_ranked$mean_Overall_weight)/multi_ATAC_ranked$sd_Overall_weight
  
  cat("multi_ATAC_ranked_POST_merge_Z_score\n")
  cat(str(multi_ATAC_ranked))
  cat("\n")
  cat(str(unique(multi_ATAC_ranked$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(multi_ATAC_ranked$Overall_weight_Z_score)))))
  cat("\n")
  cat(sprintf(as.character(summary(multi_ATAC_ranked$Overall_weight_Z_score))))
  cat("\n")
  
  #### select columns and change names to converge ----
  
  indx.int<-c(which(colnames(multi_ATAC_ranked) == "VAR"),which(colnames(multi_ATAC_ranked) == "Overall_weight"),which(colnames(multi_ATAC_ranked) == "Overall_weight_Z_score"))
  
  multi_ATAC_ranked_subset<-unique(multi_ATAC_ranked[,indx.int])
  
  cat("multi_ATAC_ranked_subset_0\n")
  cat(str(multi_ATAC_ranked_subset))
  cat("\n")
  cat(str(unique(multi_ATAC_ranked_subset$VAR)))
  cat("\n")
  
  
  colnames(multi_ATAC_ranked_subset)[which(colnames(multi_ATAC_ranked_subset) == "Overall_weight")]<-"value"
  colnames(multi_ATAC_ranked_subset)[which(colnames(multi_ATAC_ranked_subset) == "Overall_weight_Z_score")]<-"value_Z_score"
  
  multi_ATAC_ranked_subset$variable<-"multi_lineage_ATAC"
  
  
  cat("multi_ATAC_ranked_subset_1\n")
  cat(str(multi_ATAC_ranked_subset))
  cat("\n")
  cat(str(unique(multi_ATAC_ranked_subset$VAR)))
  cat("\n")
  
  check<-multi_ATAC_ranked_subset[which(multi_ATAC_ranked_subset$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  # #### SAVE ----
  
  setwd(out)
  
  saveRDS(multi_ATAC_ranked_subset, file="Prepared_file_multi_lineage_ATAC.rds")
  
  write.table(check, file="check_multi_ATAC_GLOBAL_Ranked.tsv", sep="\t", quote = F, row.names = F)
  
  
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
    make_option(c("--multi_ATAC_ranked"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--desiR_weights"), type="character", default=NULL, 
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
