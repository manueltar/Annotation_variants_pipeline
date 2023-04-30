
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
  
  
  #### Read CADD file ----
  
  setwd(out)
  
  CADD<-as.data.frame(fread(file=opt$CADD) , stringsAsFactors=F)
  
  cat("CADD_0\n")
  cat(str(CADD))
  cat("\n")
  cat(str(unique(CADD$VAR)))
  cat("\n")
  
  #### calculate mean and sd per lineage -----
  
  CADD$mean_CADD<-mean(CADD$cadd_raw, na.rm =T)
  CADD$sd_CADD<-sd(CADD$cadd_raw, na.rm =T)
  
  
  cat("CADD_POST_merge\n")
  cat(str(CADD))
  cat("\n")
  cat(str(unique(CADD$VAR)))
  cat("\n")
  
  
  CADD$cadd_raw_Z_score<-(CADD$cadd_raw-CADD$mean_CADD)/CADD$sd_CADD
  
  cat("CADD_POST_merge_Z_score\n")
  cat(str(CADD))
  cat("\n")
  cat(str(unique(CADD$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(CADD$cadd_raw_Z_score)))))
  cat("\n")
  cat(sprintf(as.character(summary(CADD$cadd_raw_Z_score))))
  cat("\n")
  
  #### select columns and change names to converge ----
  
  indx.int<-c(which(colnames(CADD) == "VAR"),which(colnames(CADD) == "cadd_raw"),which(colnames(CADD) == "cadd_raw_Z_score"))
  
  CADD_subset<-CADD[,indx.int]
  
  cat("CADD_subset_0\n")
  cat(str(CADD_subset))
  cat("\n")
  cat(str(unique(CADD_subset$VAR)))
  cat("\n")
  
  
  colnames(CADD_subset)[which(colnames(CADD_subset) == "cadd_raw")]<-"value"
  colnames(CADD_subset)[which(colnames(CADD_subset) == "cadd_raw_Z_score")]<-"value_Z_score"
  
  CADD_subset$variable<-"CADD_raw"
  
  
  cat("CADD_subset_1\n")
  cat(str(CADD_subset))
  cat("\n")
  cat(str(unique(CADD_subset$VAR)))
  cat("\n")
  
  check<-CADD_subset[which(CADD_subset$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  # #### SAVE ----
  
  setwd(out)
  
  saveRDS(CADD_subset, file="Prepared_file_CADD.rds")
  
  write.table(check, file="check_CADD.tsv", sep="\t", quote = F, row.names = F)
  
  
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
    make_option(c("--CADD"), type="character", default=NULL, 
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
