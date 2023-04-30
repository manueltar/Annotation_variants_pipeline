
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

Data_wrangling_COGS = function(option_list)
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
  
  # #### READ and transform desiR_weights ----
  # 
  # desiR_weights = as.numeric(unlist(strsplit(opt$desiR_weights, split=",")))
  # 
  # cat("desiR_weights_\n")
  # cat(sprintf(as.character(desiR_weights)))
  # cat("\n")
  # cat(str(desiR_weights))
  # cat("\n")
  
  #### Read COGS file ----
  
  COGS<-as.data.frame(fread(file=opt$COGS,sep="\t") , stringsAsFactors=F)
  
  cat("COGS_0\n")
  cat(str(COGS))
  cat("\n")
  cat(str(unique(COGS$VAR)))
  cat("\n")
  
 
  
 
  #### Represent the score ----
  
  graph_path<-paste(out,'COGS_graphs','/',sep='')
  setwd(graph_path)
  
  pdf(file="COGS_MAX.pdf")
  hist(COGS$cogs, breaks=50, col="grey", border="white", main="",
       xlab="COGS Max value per gene")
  dev.off()
  
  #### NOT Maximize COGS per variant and gene ----
 
  graph_path<-paste(out,'COGS_graphs','/',sep='')
  setwd(graph_path)
  
  pdf(file="COGS_ALL.pdf")
  hist(COGS$cogs, breaks=50, col="grey", border="white", main="",
       xlab="COGS value per gene")
  dev.off()
 
  # ######################################################################################################################
  # quit(status = 1)                            
  # 
 
  # ################################## desiR COGS_value ---------------------------------------------
  # 
  # COGS_value_LOW<-desiR_weights[3]
  # COGS_value_HIGH<-desiR_weights[4]
  # 
  # cat("desiR_weights_COGS_value\n")
  # cat(sprintf(as.character(c(COGS_value_LOW,COGS_value_HIGH))))
  # cat("\n")
  # 
  # 
  # setwd(graph_path)
  # 
  # pdf(file="desIR_COGS_value.pdf")
  # hist(COGS$COGS_value, breaks=50, col="grey", border="white", main="",
  #      xlab="COGS-Seq in a Relevant Cell Type")
  # des.line(COGS$COGS_value, "d.high", des.args=c(cut1=COGS_value_LOW, cut2=COGS_value_HIGH, scale=0.5))
  # dev.off()
  # 
  # COGS$COGS_value_weight <- d.high(COGS$COGS_value, cut1=COGS_value_LOW, cut2=COGS_value_HIGH, scale=0.5)
  # 
  # 
  # #### Overall desirability ---- 
  # 
  # Overall_nCells<-desiR_weights[5]
  # Overall_COGS_value<-desiR_weights[6]
  # 
  # cat("desiR_weights_Overall\n")
  # cat(sprintf(as.character(c(Overall_nCells,Overall_COGS_value))))
  # cat("\n")
  # 
  # 
  # 
  # COGS$cogs <- d.overall(COGS$nCells_weight, COGS$COGS_value_weight, 
  #                                         weights=c(Overall_nCells,Overall_COGS_value))
  # 
  # setwd(graph_path)
  # 
  # par(las=1)
  # pdf(file= "desIR_COGS_overall.pdf")
  # plot(rev(sort(COGS$cogs)), type="l", xlab="Rank", ylab="Overall Desirability")
  # dev.off()
  
  
  check<-COGS[which(COGS$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(check$VAR)))
  cat("\n")
  cat(sprintf(as.character(check$nCells)))
  cat("\n")
  cat(sprintf(as.character(check$COGS_value)))
  cat("\n")

  # check<-COGS[which(COGS$VAR%in%tracking_variants),]
  # 
  # 
  # cat("check_0\n")
  # cat(str(check))
  # cat("\n")
  # cat(sprintf(as.character(check$VAR)))
  # cat("\n")
  # cat(sprintf(as.character(check$nCells)))
  # cat("\n")
  # cat(sprintf(as.character(check$COGS_value)))
  # cat("\n")
  
  # #### SAVE ----
  
  setwd(out)
  
  # write.table(COGS, file="COGS_GLOBAL_assessed.tsv", sep="\t", quote = F, row.names = F)
  
  # write.table(COGS, file="COGS_GLOBAL_assessed.tsv", sep="\t", quote = F, row.names = F)
  
  
  write.table(check, file="check_COGS_GLOBAL_assessed.tsv", sep="\t", quote = F, row.names = F)
  
  
  
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
  
  
  #### Read COGS file ----
  
  setwd(out)
  
  COGS<-as.data.frame(fread(file=opt$COGS,sep="\t") , stringsAsFactors=F)
  
  
  cat("COGS_0\n")
  cat(str(COGS))
  cat("\n")
  cat(str(unique(COGS$VAR)))
  cat("\n")
  
  #### calculate mean and sd per lineage -----
  
  COGS$mean_cogs<-mean(COGS$cogs, na.rm =T)
  COGS$sd_cogs<-sd(COGS$cogs, na.rm =T)
  
  
  cat("COGS_POST_merge\n")
  cat(str(COGS))
  cat("\n")
  cat(str(unique(COGS$VAR)))
  cat("\n")
  
  
  COGS$cogs_Z_score<-(COGS$cogs-COGS$mean_cogs)/COGS$sd_cogs
  
  cat("COGS_POST_merge_Z_score\n")
  cat(str(COGS))
  cat("\n")
  cat(str(unique(COGS$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(COGS$cogs_Z_score)))))
  cat("\n")
  cat(sprintf(as.character(summary(COGS$cogs_Z_score))))
  cat("\n")
  
  #### select columns and change names to converge ----
  
  indx.int<-c(which(colnames(COGS) == "VAR"),which(colnames(COGS) == "ensembl_gene_id"),which(colnames(COGS) == "HGNC"),which(colnames(COGS) == "phenotype"),which(colnames(COGS) == "cogs"),which(colnames(COGS) == "cogs_Z_score"))
  
  COGS_subset<-COGS[,indx.int]
  
  cat("COGS_subset_0\n")
  cat(str(COGS_subset))
  cat("\n")
  cat(str(unique(COGS_subset$VAR)))
  cat("\n")
  
  
  colnames(COGS_subset)[which(colnames(COGS_subset) == "cogs")]<-"value"
  colnames(COGS_subset)[which(colnames(COGS_subset) == "cogs_Z_score")]<-"value_Z_score"
  
  COGS_subset$variable<-"COGS"
  
  
  cat("COGS_subset_1\n")
  cat(str(COGS_subset))
  cat("\n")
  cat(str(unique(COGS_subset$VAR)))
  cat("\n")
  
  check<-COGS_subset[which(COGS_subset$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  # #### SAVE ----
  
  setwd(out)
  
  saveRDS(COGS_subset, file="Prepared_file_COGS.rds")
  
  write.table(check, file="check_COGS.tsv", sep="\t", quote = F, row.names = F)
  
  
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
    make_option(c("--COGS"), type="character", default=NULL, 
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
  
  Data_wrangling_COGS(opt)
  Z_score_normalization_and_convergence(opt)
  
}


###########################################################################

system.time( main() )
