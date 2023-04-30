
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
  
  #### READ and transform excluded_phenotypes ----
  
  excluded_phenotypes = unlist(strsplit(opt$excluded_phenotypes, split=","))
  
  cat("excluded_phenotypes_\n")
  cat(sprintf(as.character(excluded_phenotypes)))
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
  
  ### Exclude percentage phenotypes of white cells -----
  
  ALL_dB_subset_restricted<-unique(ALL_dB_subset[-which(ALL_dB_subset$phenotype%in%excluded_phenotypes),])
  
  
  cat("ALL_dB_subset_restricted_0\n")
  cat(str(ALL_dB_subset_restricted))
  cat("\n")
  cat(str(unique(ALL_dB_subset_restricted$VAR)))
  cat("\n")
  
  indx.int<-c(which(colnames(ALL_dB_subset_restricted) == "VAR"),which(colnames(ALL_dB_subset_restricted) == "phenotype"))
  
  ALL_dB_double_subset<-unique(ALL_dB_subset_restricted[,indx.int])
  
  cat("ALL_dB_double_subset_0\n")
  cat(str(ALL_dB_double_subset))
  cat("\n")
  cat(str(unique(ALL_dB_double_subset$VAR)))
  cat("\n")
  
  #### Read Trait_to_Lineage_table file ----
  
  Trait_to_Lineage_table<-as.data.frame(fread(file=opt$Trait_to_Lineage_table,sep="\t") , stringsAsFactors=F)
  
  cat("Trait_to_Lineage_table_0\n")
  cat(str(Trait_to_Lineage_table))
  cat("\n")
  
  #### Read Lineage_to_Cell_table file ----
  
  Lineage_to_Cell_table<-as.data.frame(fread(file=opt$Lineage_to_Cell_table,sep="\t") , stringsAsFactors=F)
  
  cat("Lineage_to_Cell_table_0\n")
  cat(str(Lineage_to_Cell_table))
  cat("\n")
 
  #### Read ATAC_INITIAL file ----
  
  ATAC_INITIAL<-as.data.frame(fread(file=opt$ATAC_INITIAL,sep=",") , stringsAsFactors=F)
  
  cat("ATAC_INITIAL_0\n")
  cat(str(ATAC_INITIAL))
  cat("\n")
  cat(str(unique(ATAC_INITIAL$VAR)))
  cat("\n")
  
 indx.dep<-which(colnames(ATAC_INITIAL) == "LineageDEF")
 
 ATAC_INITIAL_subset<-unique(ATAC_INITIAL[,-indx.dep])
 
 cat("ATAC_INITIAL_subset_0\n")
 cat(str(ATAC_INITIAL_subset))
 cat("\n")
 cat(str(unique(ATAC_INITIAL_subset$VAR)))
 cat("\n")
 
  
  ##### LOOP TO GET AGGREGATE measures-----
  
  phenotypes_array<-unique(ALL_dB_double_subset$phenotype)
  
  CONDITION_DEBUG<-0
  
  FINAL_df<-data.frame()
  
  for(i in 1:length(phenotypes_array))
  {
    phenotypes_array_sel<-phenotypes_array[i]
    
    cat("------------------------------------------->\t")
    cat(sprintf(as.character(phenotypes_array_sel)))
    cat("\n")
    
    if(phenotypes_array_sel == "mono")
    {
      CONDITION_DEBUG<-1

    }else{

      CONDITION_DEBUG<-0
    }
    
    ALL_dB_double_subset_sel<-ALL_dB_double_subset[which(ALL_dB_double_subset$phenotype == phenotypes_array_sel),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("ALL_dB_double_subset_sel_0\n")
      cat(str(ALL_dB_double_subset_sel))
      cat("\n")
      cat(str(unique(ALL_dB_double_subset_sel$VAR)))
      cat("\n")
      
    }
   
    
    Trait_to_Lineage_table_sel<-Trait_to_Lineage_table[which(Trait_to_Lineage_table$Trait == phenotypes_array_sel),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Trait_to_Lineage_table_sel_0\n")
      cat(str(Trait_to_Lineage_table_sel))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(unique(Trait_to_Lineage_table_sel$Factor4)))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(unique(Trait_to_Lineage_table_sel$Factor4))))))
      cat("\n")
      
    }
    
    Lineage_sel<-unique(Trait_to_Lineage_table_sel$Factor4)
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Lineage_sel_0\n")
      cat(str(Lineage_sel))
      cat("\n")
    }
    
    Lineage_to_Cell_table_sel<-Lineage_to_Cell_table[which(Lineage_to_Cell_table$Lineage%in%Trait_to_Lineage_table_sel$Factor4),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Lineage_to_Cell_table_sel_0\n")
      cat(str(Lineage_to_Cell_table_sel))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(unique(Lineage_to_Cell_table_sel$CellType)))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(unique(Lineage_to_Cell_table_sel$CellType))))))
      cat("\n")
      
    }
    
    
    ATAC_INITIAL_subset_sel<-ATAC_INITIAL_subset[which(ATAC_INITIAL_subset$VAR%in%ALL_dB_double_subset_sel$VAR),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("ATAC_INITIAL_subset_sel_0\n")
      cat(str(ATAC_INITIAL_subset_sel))
      cat("\n")
      cat(str(unique(ATAC_INITIAL_subset_sel$VAR)))
      cat("\n")
     
    
    }
    
    
    
    ATAC_INITIAL_subset_sel$Tag<-NA
    
    
    ATAC_INITIAL_subset_sel$Tag[which(ATAC_INITIAL_subset_sel$ATAC_Cell_Type%in%Lineage_to_Cell_table_sel$CellType)]<-"Relevant"
    ATAC_INITIAL_subset_sel$Tag[-which(ATAC_INITIAL_subset_sel$ATAC_Cell_Type%in%Lineage_to_Cell_table_sel$CellType)]<-"Not_relevant"
    
    if(CONDITION_DEBUG == 1)
    {
      cat("ATAC_INITIAL_subset_sel_1\n")
      cat(str(ATAC_INITIAL_subset_sel))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(ATAC_INITIAL_subset_sel$Tag))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ATAC_INITIAL_subset_sel$Tag)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(ATAC_INITIAL_subset_sel$ATAC_Cell_Type))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ATAC_INITIAL_subset_sel$ATAC_Cell_Type)))))
      cat("\n")
    }
    #### keep only cells associated to the lineage ----
    
    ATAC_INITIAL_subset_double_sel<-ATAC_INITIAL_subset_sel[which(ATAC_INITIAL_subset_sel$Tag == "Relevant"),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("ATAC_INITIAL_subset_double_sel_1\n")
      cat(str(ATAC_INITIAL_subset_double_sel))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(ATAC_INITIAL_subset_double_sel$Tag))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ATAC_INITIAL_subset_double_sel$Tag)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(ATAC_INITIAL_subset_double_sel$ATAC_Cell_Type))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ATAC_INITIAL_subset_double_sel$ATAC_Cell_Type)))))
      cat("\n")
    }
    
    
    indx.int<-c(which(colnames(ATAC_INITIAL_subset_double_sel) == "VAR"),
                which(colnames(ATAC_INITIAL_subset_double_sel) == "ATAC_Cell_Type"),
                which(colnames(ATAC_INITIAL_subset_double_sel) == "value"))
    
    
    ATAC_INITIAL_subset_FINAL_sel<-unique(ATAC_INITIAL_subset_double_sel[,indx.int])
    
    ATAC_INITIAL_subset_FINAL_sel$Lineage<-Lineage_sel
    
    if(CONDITION_DEBUG == 1)
    {
      cat("ATAC_INITIAL_subset_FINAL_sel_1\n")
      cat(str(ATAC_INITIAL_subset_FINAL_sel))
      cat("\n")
      
    }
    
    FINAL_df<-unique(rbind(FINAL_df,ATAC_INITIAL_subset_FINAL_sel))
    
    if(CONDITION_DEBUG == 1)
    {
      cat("FINAL_df_1\n")
      cat(str(FINAL_df))
      cat("\n")
      
    }
   
    # if(dim(ATAC_INITIAL_subset_double_sel)[1] >0) 
    # {
      
      # quit(status = 1)
    # }
   
   
  }# i in 1:length(phenotypes_array)
  
  
  
  
  CONDITION_DEBUG<-1
  
  if(CONDITION_DEBUG == 1)
  {
    cat("FINAL_df_0\n")
    cat(str(FINAL_df))
    cat("\n")
    cat(str(unique(FINAL_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df$value)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df$value))))
    cat("\n")
    
  }
  
  check_NA<-FINAL_df[is.na(FINAL_df$value),]
  
  if(CONDITION_DEBUG == 1)
  {
    cat("check_NA_0\n")
    cat(str(check_NA))
    cat("\n")
    cat(str(unique(check_NA$VAR)))
    cat("\n")
    cat(str(unique(check_NA$phenotype)))
    cat("\n")
  }
  
  FINAL_df_NO_NA<-FINAL_df[!is.na(FINAL_df$value),]
  
  if(CONDITION_DEBUG == 1)
  {
    cat("FINAL_df_NO_NA_0\n")
    cat(str(FINAL_df_NO_NA))
    cat("\n")
    cat(str(unique(FINAL_df_NO_NA$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df_NO_NA$value)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df_NO_NA$value))))
    cat("\n")
    
  }
  
  
  check<-FINAL_df_NO_NA[which(FINAL_df_NO_NA$VAR%in%tracking_variants),]
  
  if(CONDITION_DEBUG == 1)
  {
    cat("check_0\n")
    cat(str(check))
    cat("\n")
    cat(str(unique(check$VAR)))
    cat("\n")
  }
  
  
  
  
  
  
  
  
 ######################### SAVE -----
  
  setwd(out)
  
  write.table(FINAL_df_NO_NA, file="ATAC_GLOBAL_preranked.tsv", sep="\t", quote = F, row.names = F)
  
  write.table(check, file="check_ATAC_GLOBAL_preranked.tsv", sep="\t", quote = F, row.names = F)
  
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
    make_option(c("--ATAC_INITIAL"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Trait_to_Lineage_table"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Lineage_to_Cell_table"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--excluded_phenotypes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--tracking_variants"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--relevant_not_relevant_weights"), type="character", default=NULL, 
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
  
}


###########################################################################

system.time( main() )
