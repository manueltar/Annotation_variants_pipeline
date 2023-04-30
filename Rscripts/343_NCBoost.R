
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
  
  
  #### Read NCBoost file ----
  
  setwd(out)
  
  NCBoost<-as.data.frame(fread(file=opt$NCBoost) , stringsAsFactors=F)
  
  cat("NCBoost_0\n")
  cat(str(NCBoost))
  cat("\n")
  cat(str(unique(NCBoost$VAR)))
  cat("\n")
  
  #### calculate mean and sd per lineage -----
  
  NCBoost$mean_NCBoost<-mean(NCBoost$NCBoost, na.rm =T)
  NCBoost$sd_NCBoost<-sd(NCBoost$NCBoost, na.rm =T)
  
  
  cat("NCBoost_POST_merge\n")
  cat(str(NCBoost))
  cat("\n")
  cat(str(unique(NCBoost$VAR)))
  cat("\n")
  
  
  NCBoost$NCBoost_Z_score<-(NCBoost$NCBoost-NCBoost$mean_NCBoost)/NCBoost$sd_NCBoost
  
  cat("NCBoost_POST_merge_Z_score\n")
  cat(str(NCBoost))
  cat("\n")
  cat(str(unique(NCBoost$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(NCBoost$NCBoost_Z_score)))))
  cat("\n")
  cat(sprintf(as.character(summary(NCBoost$NCBoost_Z_score))))
  cat("\n")
  
  #### select columns and change names to converge ----
  
  indx.int<-c(which(colnames(NCBoost) == "VAR"),which(colnames(NCBoost) == "NCBoost"),which(colnames(NCBoost) == "NCBoost_Z_score"))
  
  NCBoost_subset<-NCBoost[,indx.int]
  
  cat("NCBoost_subset_0\n")
  cat(str(NCBoost_subset))
  cat("\n")
  cat(str(unique(NCBoost_subset$VAR)))
  cat("\n")
  
  
  colnames(NCBoost_subset)[which(colnames(NCBoost_subset) == "NCBoost")]<-"value"
  colnames(NCBoost_subset)[which(colnames(NCBoost_subset) == "NCBoost_Z_score")]<-"value_Z_score"
  
  NCBoost_subset$variable<-"NCBoost"
  
  
  cat("NCBoost_subset_1\n")
  cat(str(NCBoost_subset))
  cat("\n")
  cat(str(unique(NCBoost_subset$VAR)))
  cat("\n")
  
  check<-NCBoost_subset[which(NCBoost_subset$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  # #### SAVE ----
  
  setwd(out)
  
  saveRDS(NCBoost_subset, file="Prepared_file_NCBoost.rds")
  
  write.table(check, file="check_NCBoost.tsv", sep="\t", quote = F, row.names = F)
  
  
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
    make_option(c("--NCBoost"), type="character", default=NULL, 
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
  
 
  Z_score_normalization_and_convergence(opt)
  
}


###########################################################################

system.time( main() )
