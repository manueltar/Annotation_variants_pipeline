
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
  
  
  #### Read SpliceAI file ----
  
 
  SpliceAI<-as.data.frame(fread(file=opt$SpliceAI) , stringsAsFactors=F)
  
  cat("SpliceAI_0\n")
  cat(str(SpliceAI))
  cat("\n")
  cat(str(unique(SpliceAI$VAR)))
  cat("\n")
  
  
  indx.keep<-c(which(colnames(SpliceAI) == "VAR"),which(colnames(SpliceAI) == "ensembl_gene_id"),which(colnames(SpliceAI) == "HGNC"))
  
  SpliceAI.m<-melt(SpliceAI, id.vars=colnames(SpliceAI)[indx.keep], variable.name="variable", value.name="value")
  
  SpliceAI.m$value<-as.numeric(SpliceAI.m$value)
  
  cat("SpliceAI.m_0\n")
  cat(str(SpliceAI.m))
  cat("\n")
  cat(str(unique(SpliceAI.m$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(SpliceAI.m$variable)))))
  cat("\n")
  cat(sprintf(as.character(summary(SpliceAI.m$variable))))
  cat("\n")
  
  
  indx.dep<-grep("_Pos",SpliceAI.m$variable)
  
  SpliceAI.m_subset<-droplevels(SpliceAI.m[-indx.dep,])
  
  cat("SpliceAI.m_subset_0\n")
  cat(str(SpliceAI.m_subset))
  cat("\n")
  cat(str(unique(SpliceAI.m_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(SpliceAI.m_subset$variable)))))
  cat("\n")
  cat(sprintf(as.character(summary(SpliceAI.m_subset$variable))))
  cat("\n")
  
  check_NA<-SpliceAI.m_subset[is.na(SpliceAI.m_subset$value),]
  
  cat("check_NA_0\n")
  cat(str(check_NA))
  cat("\n")
  cat(str(unique(check_NA$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(check_NA$variable)))))
  cat("\n")
  cat(sprintf(as.character(summary(check_NA$variable))))
  cat("\n")
  
  SpliceAI.m_subset_NO_NA<-SpliceAI.m_subset[!is.na(SpliceAI.m_subset$value),]
  
  cat("SpliceAI.m_subset_NO_NA_0\n")
  cat(str(SpliceAI.m_subset_NO_NA))
  cat("\n")
  cat(str(unique(SpliceAI.m_subset_NO_NA$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(SpliceAI.m_subset_NO_NA$variable)))))
  cat("\n")
  cat(sprintf(as.character(summary(SpliceAI.m_subset_NO_NA$variable))))
  cat("\n")
  
  
  #### calculate mean and sd per lineage -----
  
  SpliceAI.m_subset_NO_NA.dt<-data.table(SpliceAI.m_subset_NO_NA, key=c("variable"))
  
  
  SpliceAI.m_subset_NO_NA_SpliceAI_Type_parameters<-as.data.frame(SpliceAI.m_subset_NO_NA.dt[,.(mean_value=mean(value, na.rm =T),
                                                                  sd_value=sd(value, na.rm =T)),
                                                               by=key(SpliceAI.m_subset_NO_NA.dt)], stringsAsFactors=F)
  
  cat("SpliceAI.m_subset_NO_NA_SpliceAI_Type_parameters_0\n")
  cat(str(SpliceAI.m_subset_NO_NA_SpliceAI_Type_parameters))
  cat("\n")
  cat(str(unique(SpliceAI.m_subset_NO_NA_SpliceAI_Type_parameters$VAR)))
  cat("\n")
  
  #### Merge and calculate Z-score -----
  
  SpliceAI.m_subset_NO_NA<-merge(SpliceAI.m_subset_NO_NA,
                     SpliceAI.m_subset_NO_NA_SpliceAI_Type_parameters,
                     by=c("variable"))
  
  cat("SpliceAI.m_subset_NO_NA_POST_merge\n")
  cat(str(SpliceAI.m_subset_NO_NA))
  cat("\n")
  cat(str(unique(SpliceAI.m_subset_NO_NA$VAR)))
  cat("\n")
  
  
  SpliceAI.m_subset_NO_NA$value_Z_score<-(SpliceAI.m_subset_NO_NA$value-SpliceAI.m_subset_NO_NA$mean_value)/SpliceAI.m_subset_NO_NA$sd_value
  
  cat("SpliceAI.m_subset_NO_NA_POST_merge_Z_score\n")
  cat(str(SpliceAI.m_subset_NO_NA))
  cat("\n")
  cat(str(unique(SpliceAI.m_subset_NO_NA$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(SpliceAI.m_subset_NO_NA$value_Z_score)))))
  cat("\n")
  cat(sprintf(as.character(summary(SpliceAI.m_subset_NO_NA$value_Z_score))))
  cat("\n")
  
  
 
  # #### SAVE ----
  
  setwd(out)
  
  saveRDS(SpliceAI.m_subset_NO_NA, file="Prepared_file_SpliceAI.rds")
  
  
  
  
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
    make_option(c("--SpliceAI"), type="character", default=NULL, 
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
                        --SpliceAI FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
 
  Z_score_normalization_and_convergence(opt)
  
}


###########################################################################

system.time( main() )
