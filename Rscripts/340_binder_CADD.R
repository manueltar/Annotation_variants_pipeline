
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
suppressMessages(library("R.methodsS3", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("R.oo", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("R.utils", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))


opt = NULL

options(warn = 1)

Data_wrangling_CADD = function(option_list)
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
  
  #### Read ALL_dB file ----
  
  ALL_dB<-as.data.frame(fread(file=opt$ALL_dB,sep="\t") , stringsAsFactors=F)
  
  cat("ALL_dB_0\n")
  cat(str(ALL_dB))
  cat("\n")
  
  
  indx.dep<-c(which(colnames(ALL_dB) == "maf_origin"))
  
  ALL_dB_subset<-unique(ALL_dB[,-indx.dep])
  
  # cat("ALL_dB_subset\n")
  # cat(str(ALL_dB_subset))
  # cat("\n")
  
  rm(ALL_dB)
  
  ALL_dB_subset$Allelic_Series_ID<-paste(ALL_dB_subset$phenotype,ALL_dB_subset$block_no,sep='_')
  
  cat("ALL_dB_subset_1\n")
  cat(str(ALL_dB_subset))
  cat("\n")
  cat(str(unique(ALL_dB_subset$VAR)))
  cat("\n")
  
 
  ### Read the CADD_result ----
  
 
  
  CADD_result = as.data.frame(fread(opt$CADD_result, sep=",", header = T, drop=1) ,stringsAsFactors=F)
  
  CADD_result$label<-gsub("\"","",CADD_result$label)
  colnames(CADD_result)[which(colnames(CADD_result) == "label")]<-"rs"
  
  CADD_result$allele<-gsub("\"","",CADD_result$allele)
  colnames(CADD_result)[which(colnames(CADD_result) == "allele")]<-"allele"
  
  colnames(CADD_result)[which(colnames(CADD_result) == "pos")]<-"pos37"
  
  CADD_result$chr<-paste('chr',CADD_result$chr, sep='')
  
  
  cat("CADD_result_0\n")
  cat(str(CADD_result))
  cat("\n")
  
  
  #### Merge ------
  
  Merge_table<-merge(ALL_dB_subset,
                     CADD_result,
                     by=c("chr","pos37","rs"))
  
  cat("Merge_table_0\n")
  cat(str(Merge_table))
  cat("\n")
  cat(str(unique(Merge_table$VAR)))
  cat("\n")
  
  indx.int<-c(which(colnames(Merge_table) == "VAR"),which(colnames(Merge_table) == "cadd_phred"),which(colnames(Merge_table) == "cadd_raw"))
  
  Merge_table_subset<-unique(Merge_table[,indx.int])
  
  cat("Merge_table_subset_0\n")
  cat(str(Merge_table_subset))
  cat("\n")
  cat(str(unique(Merge_table_subset$VAR)))
  cat("\n")
 
  check<-Merge_table_subset[which(Merge_table_subset$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  #### exclude NA's ----
  
  Merge_table_subset_NO_NA<-Merge_table_subset[!is.na(Merge_table_subset$cadd_raw),]
  
  
  cat("Merge_table_subset_NO_NA_0\n")
  cat(str(Merge_table_subset_NO_NA))
  cat("\n")
  cat(str(unique(Merge_table_subset_NO_NA$VAR)))
  cat("\n")
  
  # #### SAVE ----
  
  setwd(out)
  
  write.table(Merge_table_subset_NO_NA, file="CADD_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
  write.table(check, file="check_CADD_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
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
    make_option(c("--ALL_dB"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CADD_result"), type="character", default=NULL, 
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
  
  Data_wrangling_CADD(opt)
  
}


###########################################################################

system.time( main() )
