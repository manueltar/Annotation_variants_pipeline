
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

Data_wrangling_PCHiC = function(option_list)
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
  
  #### Read Trait_to_CT_table file ----
  
  Trait_to_CT_table<-as.data.frame(fread(file=opt$Trait_to_CT_table,sep="\t") , stringsAsFactors=F)
  
  cat("Trait_to_CT_table_0\n")
  cat(str(Trait_to_CT_table))
  cat("\n")
  
  
  #### convert cell types to labels in PCHiC file-----
  
  Trait_to_CT_table$Cell_Type<-NA
  
  Trait_to_CT_table$Cell_Type<-revalue(Trait_to_CT_table$RelevantCellType,
                                       c("MK"="MK",
                                         "Ery"="EB",
                                         "tCD4"="tCD4",
                                         "nCD4"="VB-CD4",
                                         "aCD4"="aCD4",
                                         "naCD4"="naCD4",
                                         "tCD8"="tCD8",
                                         "nCD8"="nCD8",
                                         "tB"="tB",
                                         "nB"="nB",
                                         "Mon"="VB-MONO",
                                         "Mac0"="MAC0",
                                         "Mac1"="MAC1",
                                         "Mac2"="MAC2",
                                         "Neu"="VB-NEUT",
                                         "EP"="EP",
                                         "FoeT"="FoeT"))
  
  cat("Trait_to_CT_table_1\n")
  str(Trait_to_CT_table)
  cat("\n")
  
  check_NA<-Trait_to_CT_table[is.na(Trait_to_CT_table$Cell_Type),]
  
  cat("check_NA_0\n")
  str(check_NA)
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(check_NA$RelevantCellType))))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(check_NA$RelevantCellType))))))
  cat("\n")
  
  #### Read PCHiC_INITIAL file ----
  
  PCHiC_INITIAL<-as.data.frame(fread(file=opt$PCHiC_INITIAL,sep=",") , stringsAsFactors=F)
  
  cat("PCHiC_INITIAL_0\n")
  cat(str(PCHiC_INITIAL))
  cat("\n")
  cat(str(unique(PCHiC_INITIAL$VAR)))
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
   
    
    Trait_to_CT_table_sel<-Trait_to_CT_table[which(Trait_to_CT_table$Trait == phenotypes_array_sel),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Trait_to_CT_table_sel_0\n")
      cat(str(Trait_to_CT_table_sel))
      cat("\n")
    }
    
    
    PCHiC_INITIAL_sel<-PCHiC_INITIAL[which(PCHiC_INITIAL$VAR%in%ALL_dB_double_subset_sel$VAR),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("PCHiC_INITIAL_sel_0\n")
      cat(str(PCHiC_INITIAL_sel))
      cat("\n")
      cat(str(unique(PCHiC_INITIAL_sel$VAR)))
      cat("\n")
    }
    
    PCHiC_INITIAL_sel$Tag<-NA
    
    
    PCHiC_INITIAL_sel$Tag[which(PCHiC_INITIAL_sel$CellType_DEF%in%Trait_to_CT_table_sel$Cell_Type)]<-"Relevant"
    PCHiC_INITIAL_sel$Tag[-which(PCHiC_INITIAL_sel$CellType_DEF%in%Trait_to_CT_table_sel$Cell_Type)]<-"Not_relevant"
    
    if(CONDITION_DEBUG == 1)
    {
      cat("PCHiC_INITIAL_sel_1\n")
      cat(str(PCHiC_INITIAL_sel))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(PCHiC_INITIAL_sel$Tag))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(PCHiC_INITIAL_sel$Tag)))))
      cat("\n")
    }
    
   
    
    #### Aggregate Chicago score values by VAR and Tag  ----
    
    PCHiC_INITIAL_sel.dt<-data.table(PCHiC_INITIAL_sel, key=c("VAR","Tag"))
    
    
    Aggregation_table<-as.data.frame(PCHiC_INITIAL_sel.dt[,.(Aggregate_PCHiC=sum(value),
                                                                  nCells=.N), by=key(PCHiC_INITIAL_sel.dt)], stringsAsFactors=F)
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Aggregation_table_0\n")
      cat(str(Aggregation_table))
      cat("\n")
      cat(str(unique(Aggregation_table$VAR)))
      cat("\n")
      cat(sprintf(as.character(names(summary(Aggregation_table$Aggregate_PCHiC)))))
      cat("\n")
      cat(sprintf(as.character(summary(Aggregation_table$Aggregate_PCHiC))))
      cat("\n")
      cat(sprintf(as.character(names(summary(Aggregation_table$nCells)))))
      cat("\n")
      cat(sprintf(as.character(summary(Aggregation_table$nCells))))
      cat("\n")
    }
    
    #### Multiply by relevant/ not relevant values ----
    
    Aggregation_table$Multiplier<-NA
    
    Aggregation_table$Multiplier[which(Aggregation_table$Tag == "Relevant")]<-relevant_not_relevant_weights[1]
    Aggregation_table$Multiplier[which(Aggregation_table$Tag == "Not_relevant")]<-relevant_not_relevant_weights[2]
    
    Aggregation_table$Aggregate_PCHiC_multiplied<-Aggregation_table$Aggregate_PCHiC*Aggregation_table$Multiplier
    
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Aggregation_table_1\n")
      cat(str(Aggregation_table))
      cat("\n")
      cat(str(unique(Aggregation_table$VAR)))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(Aggregation_table$Multiplier))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(Aggregation_table$Multiplier)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(Aggregation_table$Aggregate_PCHiC_multiplied)))))
      cat("\n")
      cat(sprintf(as.character(summary(Aggregation_table$Aggregate_PCHiC_multiplied))))
      cat("\n")
    }
    
    
    # 
    
    #### Divide by the number of relevant and not relevant cells ----
    
    Aggregation_table$Aggregate_PCHiC_normalised<-Aggregation_table$Aggregate_PCHiC_multiplied/Aggregation_table$nCells
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Aggregation_table_2\n")
      cat(str(Aggregation_table))
      cat("\n")
      cat(str(unique(Aggregation_table$VAR)))
      cat("\n")
      cat(sprintf(as.character(names(summary(Aggregation_table$Aggregate_PCHiC_normalised)))))
      cat("\n")
      cat(sprintf(as.character(summary(Aggregation_table$Aggregate_PCHiC_normalised))))
      cat("\n")
    }
    
    #### Finally add up per variant ----
    
    Aggregation_table.dt<-data.table(Aggregation_table, key=c("VAR"))
    
    
    Aggregation_table_FINAL<-as.data.frame(Aggregation_table.dt[,.(Aggregate_PCHiC_FINAL=sum(Aggregate_PCHiC_normalised)), by=key(Aggregation_table.dt)], stringsAsFactors=F)
    
    Aggregation_table_FINAL$phenotype<-phenotypes_array_sel  
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Aggregation_table_FINAL_0\n")
      cat(str(Aggregation_table_FINAL))
      cat("\n")
      cat(sprintf(as.character(names(summary(Aggregation_table_FINAL$Aggregate_PCHiC_FINAL)))))
      cat("\n")
      cat(sprintf(as.character(summary(Aggregation_table_FINAL$Aggregate_PCHiC_FINAL))))
      cat("\n")
    }
    
   
    
    FINAL_df<-unique(rbind(FINAL_df,Aggregation_table_FINAL))
    
    if(CONDITION_DEBUG == 1)
    {
      cat("FINAL_df_0\n")
      cat(str(FINAL_df))
      cat("\n")
    }
     # quit(status = 1)
    
    
  }# i in 1:length(phenotypes_array)
  
  
  
  
  CONDITION_DEBUG<-1
  
  if(CONDITION_DEBUG == 1)
  {
    cat("FINAL_df_0\n")
    cat(str(FINAL_df))
    cat("\n")
    cat(str(unique(FINAL_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df$Aggregate_PCHiC_FINAL)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df$Aggregate_PCHiC_FINAL))))
    cat("\n")
    
  }
  
  check_NA<-FINAL_df[is.na(FINAL_df$Aggregate_PCHiC_FINAL),]
  
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
  
  FINAL_df_NO_NA<-FINAL_df[!is.na(FINAL_df$Aggregate_PCHiC_FINAL),]
  
  if(CONDITION_DEBUG == 1)
  {
    cat("FINAL_df_NO_NA_0\n")
    cat(str(FINAL_df_NO_NA))
    cat("\n")
    cat(str(unique(FINAL_df_NO_NA$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df_NO_NA$Aggregate_PCHiC_FINAL)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df_NO_NA$Aggregate_PCHiC_FINAL))))
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
  
  write.table(FINAL_df_NO_NA, file="PCHiC_GLOBAL_preranked.tsv", sep="\t", quote = F, row.names = F)
  
  write.table(check, file="check_PCHiC_GLOBAL_preranked.tsv", sep="\t", quote = F, row.names = F)
  
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
    make_option(c("--PCHiC_INITIAL"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Trait_to_CT_table"), type="character", default=NULL, 
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
  
  Data_wrangling_PCHiC(opt)
  
}


###########################################################################

system.time( main() )
