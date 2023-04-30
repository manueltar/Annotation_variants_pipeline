
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

Data_wrangling_GWAS = function(option_list)
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
  cat(str(unique(ALL_dB$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB$finemap_beta)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB$finemap_beta))))
  cat("\n")
  
  #### Obtain  Absolute_effect_size and change finemap_prob to PP ----
  
  ALL_dB$Absolute_effect_size<-abs(ALL_dB$finemap_beta)
  colnames(ALL_dB)[which(colnames(ALL_dB) == "finemap_prob")]<-"PP"
  colnames(ALL_dB)[which(colnames(ALL_dB) == "maf_origin")]<-"MAF"
  
  
  cat("ALL_dB_1\n")
  cat(str(ALL_dB))
  cat("\n")
  cat(str(unique(ALL_dB$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB$Absolute_effect_size)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB$Absolute_effect_size))))
  cat("\n")
  
  #### Max AF ----
  
  ALL_dB.dt<-data.table(ALL_dB, key="VAR")
  
  
  ALL_dB_MAX_AF<-as.data.frame(ALL_dB.dt[,.SD[which.max(MAF)], by=key(ALL_dB.dt)], stringsAsFactors=F)
  
  cat("ALL_dB_MAX_AF_0\n")
  cat(str(ALL_dB_MAX_AF))
  cat("\n")
  cat(str(unique(ALL_dB_MAX_AF$VAR)))
  cat("\n")
  
  
  indx.int<-c(which(colnames(ALL_dB_MAX_AF) == "VAR"),which(colnames(ALL_dB_MAX_AF) == "MAF"))
  
  ALL_dB_MAX_AF_subset<-unique(ALL_dB_MAX_AF[,indx.int])
  
  cat("ALL_dB_MAX_AF_subset_0\n")
  cat(str(ALL_dB_MAX_AF_subset))
  cat("\n")
  cat(str(unique(ALL_dB_MAX_AF_subset$VAR)))
  cat("\n")
  
  ALL_dB_MAX_AF_subset$mean_value<-mean(ALL_dB_MAX_AF_subset$MAF)
  ALL_dB_MAX_AF_subset$sd_value<-sd(ALL_dB_MAX_AF_subset$MAF)
  
  ALL_dB_MAX_AF_subset$MAF_Z_score<-(ALL_dB_MAX_AF_subset$MAF-ALL_dB_MAX_AF_subset$mean_value)/ALL_dB_MAX_AF_subset$sd_value
  
  
  cat("ALL_dB_MAX_AF_subset_1\n")
  cat(str(ALL_dB_MAX_AF_subset))
  cat("\n")
  cat(str(unique(ALL_dB_MAX_AF_subset$VAR)))
  cat("\n")
  
  indx.int<-c(which(colnames(ALL_dB_MAX_AF_subset) == "VAR"),which(colnames(ALL_dB_MAX_AF_subset) == "MAF"),which(colnames(ALL_dB_MAX_AF_subset) == "MAF_Z_score"))
  
  
  ALL_dB_MAX_AF_subset_double<-unique(ALL_dB_MAX_AF_subset[,indx.int])
  
  cat("ALL_dB_MAX_AF_subset_double_0\n")
  cat(str(ALL_dB_MAX_AF_subset_double))
  cat("\n")
  cat(str(unique(ALL_dB_MAX_AF_subset_double$VAR)))
  cat("\n")
  
  colnames(ALL_dB_MAX_AF_subset_double)[which(colnames(ALL_dB_MAX_AF_subset_double) == "MAF")]<-"value"
  colnames(ALL_dB_MAX_AF_subset_double)[which(colnames(ALL_dB_MAX_AF_subset_double) == "MAF_Z_score")]<-"value_Z_score"
  
  ALL_dB_MAX_AF_subset_double$variable<-"MAF"
  
  
  cat("ALL_dB_MAX_AF_subset_double_1\n")
  cat(str(ALL_dB_MAX_AF_subset_double))
  cat("\n")
  cat(str(unique(ALL_dB_MAX_AF_subset_double$VAR)))
  cat("\n")
  
 
  
  ############ Z score per phenotype Absolute_effect_size and PP-----------
  
  indx.int<-c(which(colnames(ALL_dB) == "VAR"),which(colnames(ALL_dB) == "phenotype"),which(colnames(ALL_dB) == "PP"),which(colnames(ALL_dB) == "Absolute_effect_size"))
  
  
  ALL_dB_subset<-unique(ALL_dB[,indx.int])
  
  cat("ALL_dB_subset_0\n")
  cat(str(ALL_dB_subset))
  cat("\n")
  cat(str(unique(ALL_dB_subset$VAR)))
  cat("\n")
  
  rm(ALL_dB)
  
 
  ALL_dB_subset.m<-melt(ALL_dB_subset, id.vars=c("VAR","phenotype"))
  
  
  cat("ALL_dB_subset.m_0\n")
  cat(str(ALL_dB_subset.m))
  cat("\n")
  cat(str(unique(ALL_dB_subset.m$VAR)))
  cat("\n")
  cat(sprintf(as.character(unique(ALL_dB_subset.m$variable))))
  cat("\n")
  
  #### calculate mean and sd per phenotype -----
  
  ALL_dB_subset.m.dt<-data.table(ALL_dB_subset.m, key=c("phenotype","variable"))
  
  
  ALL_dB_subset.m_GWAS_Type_parameters<-as.data.frame(ALL_dB_subset.m.dt[,.(mean_value=mean(value, na.rm =T),
                                                          sd_value=sd(value, na.rm =T)),
                                                       by=key(ALL_dB_subset.m.dt)], stringsAsFactors=F)
  
  cat("ALL_dB_subset.m_GWAS_Type_parameters_0\n")
  cat(str(ALL_dB_subset.m_GWAS_Type_parameters))
  cat("\n")
  
  
  #### Merge and calculate Z-score -----
  
  ALL_dB_subset.m<-merge(ALL_dB_subset.m,
                              ALL_dB_subset.m_GWAS_Type_parameters,
                by=c("phenotype","variable"))
  
  cat("ALL_dB_subset.m_POST_merge\n")
  cat(str(ALL_dB_subset.m))
  cat("\n")
  cat(str(unique(ALL_dB_subset.m$VAR)))
  cat("\n")
  
  
  ALL_dB_subset.m$value_Z_score<-(ALL_dB_subset.m$value-ALL_dB_subset.m$mean_value)/ALL_dB_subset.m$sd_value
   
  
  cat("ALL_dB_subset.m_POST_merge_Z_score\n")
  cat(str(ALL_dB_subset.m))
  cat("\n")
  cat(str(unique(ALL_dB_subset.m$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB_subset.m$value_Z_score)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB_subset.m$value_Z_score))))
  cat("\n")
 
  ##### subselect colums ----
  
  indx.int<-c(which(colnames(ALL_dB_subset.m) == "VAR"),which(colnames(ALL_dB_subset.m) == "phenotype"),which(colnames(ALL_dB_subset.m) == "variable"),which(colnames(ALL_dB_subset.m) == "value"),which(colnames(ALL_dB_subset.m) == "value_Z_score"))
  
  ALL_dB_subset.m_subset<-unique(ALL_dB_subset.m[,indx.int])
  
  cat("ALL_dB_subset.m_subset_0\n")
  cat(str(ALL_dB_subset.m_subset))
  cat("\n")
  cat(str(unique(ALL_dB_subset.m_subset$VAR)))
  cat("\n")
  cat(str(unique(ALL_dB_subset.m_subset$variable)))
  cat("\n")
  
  
  # #### Merge with AF and subselect columns ----
  
  
  
  cat("ALL_dB_MAX_AF_subset_double_REMIND\n")
  cat(str(ALL_dB_MAX_AF_subset_double))
  cat("\n")
  cat(str(unique(ALL_dB_MAX_AF_subset_double$VAR)))
  cat("\n")
  cat(str(unique(ALL_dB_MAX_AF_subset_double$variable)))
  cat("\n")
  
  
  # 
  # ALL_dB_subset.m_subset<-merge(ALL_dB_subset.m_subset,
  #                             ALL_dB_MAX_AF_subset,
  #                             by="VAR",
  #                             all.x=T)
  # 
  # 
  # cat("ALL_dB_subset.m_subset_1\n")
  # cat(str(ALL_dB_subset.m_subset))
  # cat("\n")
  # cat(str(unique(ALL_dB_subset.m_subset$VAR)))
  # cat("\n")
  # 
  # 
  # 
  # check<-ALL_dB_subset.m_subset[which(ALL_dB_subset.m_subset$VAR%in%tracking_variants),]
  # 
  # CONDITION_DEBUG <- 1
  # 
  # if(CONDITION_DEBUG == 1)
  # {
  #   cat("check_0\n")
  #   cat(str(check))
  #   cat("\n")
  #   cat(str(unique(check$VAR)))
  #   cat("\n")
  # }
  
  # #######################################################################
  # quit(status = 1)
  # #### SAVE ----
  
  setwd(out)
  
  write.table(ALL_dB_subset.m_subset, file="GWAS_GLOBAL_per_traits.tsv", sep="\t", quote = F, row.names = F)
  
  write.table(ALL_dB_MAX_AF_subset_double, file="MAF_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
  
  # write.table(check, file="check_GWAS_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
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
                        --GWAS FILE.txt
                        --GWAS FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  Data_wrangling_GWAS(opt)
  
}


###########################################################################

system.time( main() )
