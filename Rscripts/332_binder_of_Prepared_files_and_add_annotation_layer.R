
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

Data_wrangling = function(option_list)
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
  
  
  #### Read ATAC_data file ----
  
  ATAC_data<-readRDS(file=opt$ATAC_data)
  
  cat("ATAC_data_0\n")
  cat(str(ATAC_data))
  cat("\n")
  cat(str(unique(ATAC_data$VAR)))
  cat("\n")
  
  #### Read multi_ATAC_data file ----
  
  multi_ATAC_data<-readRDS(file=opt$multi_ATAC_data)
  
  cat("multi_ATAC_data_0\n")
  cat(str(multi_ATAC_data))
  cat("\n")
  cat(str(unique(multi_ATAC_data$VAR)))
  cat("\n")
  
  #### Read GWAS file ----
  
  MAF_GLOBAL<-as.data.frame(fread(file=opt$MAF_GLOBAL, sep="\t", header=T), stringsAsFactors=F)
  
  cat("MAF_GLOBAL_0\n")
  cat(str(MAF_GLOBAL))
  cat("\n")
  cat(str(unique(MAF_GLOBAL$VAR)))
  cat("\n")
  
  GWAS_GLOBAL_per_traits<-as.data.frame(fread(file=opt$GWAS_GLOBAL_per_traits, sep="\t", header=T), stringsAsFactors=F)
  
  cat("GWAS_GLOBAL_per_traits_0\n")
  cat(str(GWAS_GLOBAL_per_traits))
  cat("\n")
  cat(str(unique(GWAS_GLOBAL_per_traits$VAR)))
  cat("\n")
  
  ind.int<-c(which(colnames(GWAS_GLOBAL_per_traits) == "VAR"),which(colnames(GWAS_GLOBAL_per_traits) == "value"),which(colnames(GWAS_GLOBAL_per_traits) == "value_Z_score"),which(colnames(GWAS_GLOBAL_per_traits) == "variable"))
  
  GWAS_GLOBAL_per_traits_subset_1<-GWAS_GLOBAL_per_traits[,ind.int]
  
  cat("GWAS_GLOBAL_per_traits_subset_1_0\n")
  cat(str(GWAS_GLOBAL_per_traits_subset_1))
  cat("\n")
  cat(str(unique(GWAS_GLOBAL_per_traits_subset_1$VAR)))
  cat("\n")
  cat(sprintf(as.character(unique(GWAS_GLOBAL_per_traits_subset_1$variable))))
  cat("\n")
  
  #### Read Constraint_Z file ----
  
  Constraint_Z<-readRDS(file=opt$Constraint_Z)
  
  cat("Constraint_Z_0\n")
  cat(str(Constraint_Z))
  cat("\n")
  cat(str(unique(Constraint_Z$VAR)))
  cat("\n")
  
  
  #### Read NCBoost file ----
  
  NCBoost<-readRDS(file=opt$NCBoost)
  
  cat("NCBoost_0\n")
  cat(str(NCBoost))
  cat("\n")
  cat(str(unique(NCBoost$VAR)))
  cat("\n")
  
  ind.int<-c(which(colnames(NCBoost) == "VAR"),which(colnames(NCBoost) == "value"),which(colnames(NCBoost) == "value_Z_score"),which(colnames(NCBoost) == "variable"))
  
  NCBoost_subset<-NCBoost[,ind.int]
  
  cat("NCBoost_subset_0\n")
  cat(str(NCBoost_subset))
  cat("\n")
  cat(str(unique(NCBoost_subset$VAR)))
  cat("\n")
  
  #### Read SpliceAI file ----
  
  SpliceAI<-readRDS(file=opt$SpliceAI)
  
  cat("SpliceAI_0\n")
  cat(str(SpliceAI))
  cat("\n")
  cat(str(unique(SpliceAI$VAR)))
  cat("\n")
  
  ind.int<-c(which(colnames(SpliceAI) == "VAR"),which(colnames(SpliceAI) == "value"),which(colnames(SpliceAI) == "value_Z_score"),which(colnames(SpliceAI) == "variable"))
  
  SpliceAI_subset<-SpliceAI[,ind.int]
  
  cat("SpliceAI_subset_0\n")
  cat(str(SpliceAI_subset))
  cat("\n")
  cat(str(unique(SpliceAI_subset$VAR)))
  cat("\n")
  
  #### Read CADD file ----
  
  CADD<-readRDS(file=opt$CADD)
  
  cat("CADD_0\n")
  cat(str(CADD))
  cat("\n")
  cat(str(unique(CADD$VAR)))
  cat("\n")
  
  ind.int<-c(which(colnames(CADD) == "VAR"),which(colnames(CADD) == "value"),which(colnames(CADD) == "value_Z_score"),which(colnames(CADD) == "variable"))
  
  CADD_subset<-CADD[,ind.int]
  
  cat("CADD_subset_0\n")
  cat(str(CADD_subset))
  cat("\n")
  cat(str(unique(CADD_subset$VAR)))
  cat("\n")
  
  #### Read oe_LOF file ----
  
  oe_LOF<-readRDS(file=opt$oe_LOF)
  
  cat("oe_LOF_0\n")
  cat(str(oe_LOF))
  cat("\n")
  cat(str(unique(oe_LOF$VAR)))
  cat("\n")
  
  ind.int<-c(which(colnames(oe_LOF) == "VAR"),which(colnames(oe_LOF) == "value"),which(colnames(oe_LOF) == "value_Z_score"),which(colnames(oe_LOF) == "variable"))
  
  oe_LOF_subset<-oe_LOF[,ind.int]
  
  cat("oe_LOF_subset_0\n")
  cat(str(oe_LOF_subset))
  cat("\n")
  cat(str(unique(oe_LOF_subset$VAR)))
  cat("\n")
  
  #### Read COGS file ----
  
  COGS<-readRDS(file=opt$COGS)
  
  cat("COGS_0\n")
  cat(str(COGS))
  cat("\n")
  cat(str(unique(COGS$VAR)))
  cat("\n")
  
  ind.int<-c(which(colnames(COGS) == "VAR"),which(colnames(COGS) == "value"),which(colnames(COGS) == "value_Z_score"),which(colnames(COGS) == "variable"))
  
  COGS_subset<-COGS[,ind.int]
  
  cat("COGS_subset_0\n")
  cat(str(COGS_subset))
  cat("\n")
  cat(str(unique(COGS_subset$VAR)))
  cat("\n")
  
  #### Read GENE_EXP file ----
  
  GENE_EXP<-readRDS(file=opt$GENE_EXP)
  
  cat("GENE_EXP_0\n")
  cat(str(GENE_EXP))
  cat("\n")
  cat(str(unique(GENE_EXP$VAR)))
  cat("\n")
  
  ind.int<-c(which(colnames(GENE_EXP) == "VAR"),which(colnames(GENE_EXP) == "value"),which(colnames(GENE_EXP) == "value_Z_score"),which(colnames(GENE_EXP) == "variable"))
  
  GENE_EXP_subset<-GENE_EXP[,ind.int]
  
  cat("GENE_EXP_subset_0\n")
  cat(str(GENE_EXP_subset))
  cat("\n")
  cat(str(unique(GENE_EXP_subset$VAR)))
  cat("\n")
  
  #### Read PCHiC file ----
  
  PCHiC<-readRDS(file=opt$PCHiC)
  
  cat("PCHiC_0\n")
  cat(str(PCHiC))
  cat("\n")
  cat(str(unique(PCHiC$VAR)))
  cat("\n")
  
  #### Read chromstates file ----
  
  chromstates<-readRDS(file=opt$chromstates)
  
  cat("chromstates_0\n")
  cat(str(chromstates))
  cat("\n")
  cat(str(unique(chromstates$VAR)))
  cat("\n")
  
  #### Rbind ----
  
  
  df_1<-rbind(MAF_GLOBAL,GWAS_GLOBAL_per_traits_subset_1)
  
  cat("df_1_0\n")
  cat(str(df_1))
  cat("\n")
  
  check_NA_1<-df_1[is.na(df_1$value),]
  
  cat("check_NA_1_0\n")
  cat(str(check_NA_1))
  cat("\n")
  cat(sprintf(as.character(unique(check_NA_1$variable))))
  cat("\n")
  
  df_2<-rbind(CADD_subset,Constraint_Z,NCBoost_subset,SpliceAI_subset)
  
  cat("df_2_0\n")
  cat(str(df_2))
  cat("\n")
  
  
  check_NA_2<-df_2[is.na(df_2$value),]
  
  cat("check_NA_2_0\n")
  cat(str(check_NA_2))
  cat("\n")
  cat(sprintf(as.character(unique(check_NA_2$variable))))
  cat("\n")
  
  
  df_3<-rbind(ATAC_data,multi_ATAC_data,PCHiC,chromstates,COGS_subset,oe_LOF_subset,GENE_EXP_subset)
  
  cat("df_3_0\n")
  cat(str(df_3))
  cat("\n")
  
  check_NA_3<-df_3[is.na(df_3$value),]
  
  cat("check_NA_3_0\n")
  cat(str(check_NA_3))
  cat("\n")
  cat(sprintf(as.character(unique(check_NA_3$variable))))
  cat("\n")
  
  
  Master_file<-rbind(df_1,df_2,df_3)

  cat("Master_file_0\n")
  cat(str(Master_file))
  cat("\n")
  cat(str(unique(Master_file$VAR)))
  cat("\n")

  array_variables<-levels(as.factor(Master_file$variable))

  cat("array_variables_0\n")
  cat(str(array_variables))
  cat("\n")
  cat(sprintf(as.character(array_variables)))
  cat("\n")

  Master_file$variable<-factor(Master_file$variable,
                               levels=c("MAF","PP","Absolute_effect_size","CADD_raw",'Constraint_Z',"NCBoost","SpliceAI_DG","SpliceAI_DL","SpliceAI_AG","SpliceAI_AL",
                                        "Rank_ATAC_erythroid_lineage","Rank_ATAC_mega_lineage","Rank_ATAC_gran_mono_lineage","Rank_ATAC_lymph_lineage","multi_lineage_ATAC","Rank_PCHiC","Rank_chromstates",
                                        "COGS","oe_lof","Rank_GENE_EXP"),
                               ordered=T)


  cat("Master_file_1\n")
  cat(str(Master_file))
  cat("\n")
  cat(str(unique(Master_file$VAR)))
  cat("\n")
  cat(str(unique(Master_file$variable)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Master_file$variable)))))
  cat("\n")
  cat(sprintf(as.character(summary(Master_file$variable))))
  cat("\n")


  check_NA<-Master_file[is.na(Master_file$value),]

  cat("check_NA_0\n")
  cat(str(check_NA))
  cat("\n")
  cat(str(unique(check_NA$VAR)))
  cat("\n")
  cat(str(unique(check_NA$variable)))
  cat("\n")
  cat(sprintf(as.character(names(summary(check_NA$variable)))))
  cat("\n")
  cat(sprintf(as.character(summary(check_NA$variable))))
  cat("\n")
  
  
  # ##########################################################
  # quit(status = 1)
  
  ### Read VAR_Prioritization_dB----
  
  
  VAR_Prioritization_dB<-as.data.frame(readRDS(file=opt$VAR_Prioritization_dB) , stringsAsFactors=F)
  
  cat("VAR_Prioritization_dB_0\n")
  cat(str(VAR_Prioritization_dB))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB$VAR)))
  cat("\n")
  
  indx.int<-c(which(colnames(VAR_Prioritization_dB) == "VAR"),which(colnames(VAR_Prioritization_dB) == "Fig1_Annot_Category"))
  
  VAR_Prioritization_dB_subset<-unique(VAR_Prioritization_dB[,indx.int])
  
  cat("VAR_Prioritization_dB_subset_0\n")
  cat(str(VAR_Prioritization_dB_subset))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB_subset$VAR)))
  cat("\n")
  
 
  ### Read Table_S6----
  
  
  Table_S6<-as.data.frame(readRDS(file=opt$Table_S6) , stringsAsFactors=F)
  
  cat("Table_S6_0\n")
  cat(str(Table_S6))
  cat("\n")
  
  
  indx.int<-c(which(colnames(Table_S6) == "VAR"),which(colnames(Table_S6) == "Mechanistic_Class"),which(colnames(Table_S6) == "Manual_curation"),which(colnames(Table_S6) == "MPRA_CLASS"),which(colnames(Table_S6) == "genIE_CLASS"),which(colnames(Table_S6) == "Multi_Lineage"))
  
  Table_S6_subset<-unique(Table_S6[,indx.int])
  
  cat("Table_S6_subset_0\n")
  cat(str(Table_S6_subset))
  cat("\n")
  cat(str(unique(Table_S6_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$Mechanistic_Class))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$Manual_curation)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$Manual_curation))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$MPRA_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$genIE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$genIE_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$Multi_Lineage))))
  cat("\n")
  
  
  
  Table_S6_subset$M_and_M<-interaction(Table_S6_subset$Mechanistic_Class,Table_S6_subset$Manual_curation, sep="|",lex.order = T)#, na.omit=T)
  
  
  
  Table_S6_subset<-droplevels(Table_S6_subset)
  
  
  cat("Table_S6_subset_1\n")
  cat(str(Table_S6_subset))
  cat("\n")
  cat(str(unique(Table_S6_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$Mechanistic_Class))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$Manual_curation)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$Manual_curation))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$MPRA_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$genIE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$genIE_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$Multi_Lineage))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$M_and_M)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$M_and_M))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$MPRA_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$MPRA_CLASS_2))))
  cat("\n")
 
  
  #### Merge all labels ----
  
  Table_of_labels<-merge(VAR_Prioritization_dB_subset,
                         Table_S6_subset,
                         by="VAR",
                         all=T)
 
  cat("Table_of_labels_1\n")
  cat(str(Table_of_labels))
  cat("\n")
  cat(str(unique(Table_of_labels$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_of_labels$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_of_labels$Fig1_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_of_labels$MPRA_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_of_labels$MPRA_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_of_labels$genIE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_of_labels$genIE_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_of_labels$M_and_M)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_of_labels$M_and_M))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_of_labels$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_of_labels$Multi_Lineage))))
  cat("\n")
  
  ############## SAVE -------------------------
  
  setwd(out)
  
  saveRDS(Master_file,file="Master_file_scores.rds")
  
  saveRDS(Table_of_labels,file="Table_of_labels.rds")
  
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
    make_option(c("--Table_S6"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VAR_Prioritization_dB"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ATAC_data"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--multi_ATAC_data"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--COGS"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--oe_LOF"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CADD"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--NCBoost"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Constraint_Z"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--SpliceAI"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MAF_GLOBAL"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GWAS_GLOBAL_per_traits"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GENE_EXP"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--PCHiC"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--chromstates"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--tracking_variants"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--RMV_common"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--RMV_labels"), type="character", default=NULL, 
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
  
  Data_wrangling(opt)
  
}


###########################################################################

system.time( main() )
