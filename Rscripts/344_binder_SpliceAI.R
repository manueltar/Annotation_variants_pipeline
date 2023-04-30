
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
suppressMessages(library("splitstackshape", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))


opt = NULL

options(warn = 1)

Data_wrangling_SpliceAI = function(option_list)
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
  
 
  #### Read GENE_table file ----
  
  GENE_table<-as.data.frame(fread(file=opt$GENE_table,sep="\t") , stringsAsFactors=F)
  
  cat("GENE_table\n")
  cat(str(GENE_table))
  cat("\n")
  
  colnames(GENE_table)<-c("chr","start","end","ensembl_gene_id","HGNC")
  
  
  indx.int<-c(which(colnames(GENE_table) == "ensembl_gene_id"),which(colnames(GENE_table) == "HGNC"))
  
  
  GENE_table_subset<-unique(GENE_table[,indx.int])
  
  
  cat("GENE_table_subset\n")
  cat(str(GENE_table_subset))
  cat("\n")
  
  
  
  
  
  ### Read the SpliceAI_result ----
  
 
  
  SpliceAI_result = as.data.frame(fread(opt$SpliceAI_result, sep="\t", header = F, skip=30, fill=TRUE) ,stringsAsFactors=F)
  colnames(SpliceAI_result)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  SpliceAI_result$VAR<-paste(paste('chr',SpliceAI_result$CHROM,sep=''),SpliceAI_result$POS,SpliceAI_result$REF,SpliceAI_result$ALT, sep="_")
  
  cat("SpliceAI_result_0\n")
  cat(str(SpliceAI_result))
  cat("\n")
  cat(str(unique(SpliceAI_result$VAR)))
  cat("\n")
  
  
  # quit(status = 1)
  
  row.with.results<-grep('SpliceAI',SpliceAI_result$INFO)
  
  SpliceAI_result_subset<-SpliceAI_result[row.with.results,]
  
  cat("SpliceAI_result_subset_1\n")
  cat(str(SpliceAI_result_subset))
  cat("\n")
  cat(str(unique(SpliceAI_result_subset$VAR)))
  cat("\n")
  
  
  
  
  #### Merge ------
  
  
  Merge_table<-merge(ALL_dB_subset,
                     SpliceAI_result_subset,
                     by=c("VAR"))
  
  cat("Merge_table_0\n")
  cat(str(Merge_table))
  cat("\n")
  cat(str(unique(Merge_table$VAR)))
  cat("\n")
  
  VARS_check<-c("chr1_111682307_G_C","chr1_11252716_T_C","chr1_11254006_AC_A")
  
  check<-Merge_table[which(Merge_table$VAR%in%VARS_check),]
  
  
  cat("check_1\n")
  cat(str(check))
  cat("\n")
  cat(str(unique(check$VAR)))
  cat("\n")
  
  
  Merge_table_separated_LONG<-unique(as.data.frame(cSplit(Merge_table, splitCols = "INFO",
                                        sep = ",", direction = "long", drop = F),stringsAsFactors=F))
  
  cat("Merge_table_separated_LONG_0\n")
  cat(str(Merge_table_separated_LONG))
  cat("\n")
  cat(str(unique(Merge_table_separated_LONG$VAR)))
  cat("\n")
  
  check<-Merge_table_separated_LONG[which(Merge_table_separated_LONG$VAR%in%VARS_check),]
  
  
  cat("check_1\n")
  cat(str(check))
  cat("\n")
  cat(str(unique(check$VAR)))
  cat("\n")
  
  
  Merge_table_separated_WIDE<-unique(as.data.frame(cSplit(Merge_table_separated_LONG, splitCols = "INFO",
                                                          sep = "|", direction = "wide", drop = F),stringsAsFactors=F))
  
  cat("Merge_table_separated_WIDE_0\n")
  cat(str(Merge_table_separated_WIDE))
  cat("\n")
  cat(str(unique(Merge_table_separated_WIDE$VAR)))
  cat("\n")
  
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_01")]<-"ALLELE"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_02")]<-"HGNC"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_03")]<-"SpliceAI_AG"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_04")]<-"SpliceAI_AL"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_05")]<-"SpliceAI_DG"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_06")]<-"SpliceAI_DL"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_07")]<-"SpliceAI_AG_Pos"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_08")]<-"SpliceAI_AL_Pos"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_09")]<-"SpliceAI_DG_Pos"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_10")]<-"SpliceAI_DL_Pos"


  Merge_table_separated_WIDE$ALLELE<-gsub("SpliceAI=","", Merge_table_separated_WIDE$ALLELE)
  
  cat("Merge_table_separated_WIDE_1\n")
  cat(str(Merge_table_separated_WIDE))
  cat("\n")
  cat(str(unique(Merge_table_separated_WIDE$VAR)))
  cat("\n")
  
  check<-Merge_table_separated_WIDE[which(Merge_table_separated_WIDE$VAR%in%VARS_check),]
  
  
  cat("check_2\n")
  cat(str(check))
  cat("\n")
  cat(str(unique(check$VAR)))
  cat("\n")
  
  indx.int<-c(which(colnames(Merge_table_separated_WIDE) == "VAR"),which(colnames(Merge_table_separated_WIDE) == "HGNC"),
              which(colnames(Merge_table_separated_WIDE) == "SpliceAI_AG"),which(colnames(Merge_table_separated_WIDE) == "SpliceAI_AL"),which(colnames(Merge_table_separated_WIDE) == "SpliceAI_DG"),which(colnames(Merge_table_separated_WIDE) == "SpliceAI_DL"),
              which(colnames(Merge_table_separated_WIDE) == "SpliceAI_AG_Pos"),which(colnames(Merge_table_separated_WIDE) == "SpliceAI_AL_Pos"),which(colnames(Merge_table_separated_WIDE) == "SpliceAI_DG_Pos"),which(colnames(Merge_table_separated_WIDE) == "SpliceAI_DL_Pos"))
  
  Merge_table_separated_WIDE_subset<-unique(Merge_table_separated_WIDE[,indx.int])
  
  cat("Merge_table_separated_WIDE_1\n")
  cat(str(Merge_table_separated_WIDE_subset))
  cat("\n")
  cat(str(unique(Merge_table_separated_WIDE_subset$VAR)))
  cat("\n")
  
  
  Merge_table_separated_WIDE_subset<-merge(Merge_table_separated_WIDE_subset,
                                           GENE_table_subset,
                                           by="HGNC",
                                           all.x=T)
  
  
  cat("Merge_table_separated_WIDE_2\n")
  cat(str(Merge_table_separated_WIDE_subset))
  cat("\n")
  cat(str(unique(Merge_table_separated_WIDE_subset$VAR)))
  cat("\n")
  
 
  check<-Merge_table_separated_WIDE_subset[which(Merge_table_separated_WIDE_subset$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  # #### SAVE ----
  
  setwd(out)
  
  write.table(Merge_table_separated_WIDE_subset, file="SpliceAI_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
  write.table(check, file="check_SpliceAI_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
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
    make_option(c("--GENE_table"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--SpliceAI_result"), type="character", default=NULL, 
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
                        --SpliceAI FILE.txt
                        --SpliceAI FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  Data_wrangling_SpliceAI(opt)
  
}


###########################################################################

system.time( main() )
