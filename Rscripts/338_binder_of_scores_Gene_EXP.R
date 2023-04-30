
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

Data_wrangling_GENE_EXP = function(option_list)
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
  
  #### Read GENE_EXP ----
  
  GENE_EXP = read.table(opt$GENE_EXP, sep="\t", stringsAsFactors = F, header = T)
  
  cat("GENE_EXP_\n")
  str(GENE_EXP)
  cat("\n")
  
  #### Read TOME_correspondence ----
  
  TOME_correspondence = read.table(opt$TOME_correspondence, sep="\t", stringsAsFactors = F, header = T)
  
  cat("TOME_correspondence_\n")
  str(TOME_correspondence)
  cat("\n")
  
  #### Read VEP_CSQ file ----
  
  VEP_CSQ<-as.data.frame(fread(file=opt$VEP_CSQ,sep=",") , stringsAsFactors=F)
  
  cat("VEP_CSQ_0\n")
  cat(str(VEP_CSQ))
  cat("\n")
  cat(str(unique(VEP_CSQ$VAR)))
  cat("\n")
  
  indx.int<-c(which(colnames(VEP_CSQ) == "VAR"),
              which(colnames(VEP_CSQ) == "ensembl_gene_id"),
              which(colnames(VEP_CSQ) == "HGNC"))
  
  VEP_CSQ_subset<-unique(VEP_CSQ[,indx.int])
  
  cat("VEP_CSQ_subset_0\n")
  cat(str(VEP_CSQ_subset))
  cat("\n")
  cat(str(unique(VEP_CSQ_subset$VAR)))
  cat("\n")
  
  VEP_CSQ_subset_with_genes<-unique(VEP_CSQ_subset[which(VEP_CSQ_subset$ensembl_gene_id != "DUMMY"),])
  
  cat("VEP_CSQ_subset_with_genes_0\n")
  cat(str(VEP_CSQ_subset_with_genes))
  cat("\n")
  cat(str(unique(VEP_CSQ_subset_with_genes$VAR)))
  cat("\n")
  
  
  
  #### Read PCHiC file ----
  
  PCHiC<-as.data.frame(fread(file=opt$PCHiC,sep=",") , stringsAsFactors=F)
  
  cat("PCHiC_0\n")
  cat(str(PCHiC))
  cat("\n")
  cat(str(unique(PCHiC$VAR)))
  cat("\n")
  
  indx.int<-c(which(colnames(PCHiC) == "VAR"),
              which(colnames(PCHiC) == "ensembl_gene_id"),
              which(colnames(PCHiC) == "HGNC"))
  
  PCHiC_subset<-unique(PCHiC[,indx.int])
  
  cat("PCHiC_subset_0\n")
  cat(str(PCHiC_subset))
  cat("\n")
  cat(str(unique(PCHiC_subset$VAR)))
  cat("\n")
  
  ####### rbind ---------
  
  VEP_and_PCHiC<-unique(rbind(VEP_CSQ_subset_with_genes,PCHiC_subset))
  
  cat("VEP_and_PCHiC_0\n")
  cat(str(VEP_and_PCHiC))
  cat("\n")
  cat(str(unique(VEP_and_PCHiC$VAR)))
  cat("\n")
  
  ### LOOP and retrieve Gene EXP ----
  
  phenotypes_array<-unique(ALL_dB_subset_restricted$phenotype)
  
  cat("phenotypes_array_0\n")
  cat(str(phenotypes_array))
  cat("\n")
  
  CONDITION_DEBUG<-0
    
  Gather<-data.frame()
  
  for(i in 1:length(phenotypes_array))
  {
    phenotypes_array_sel<-phenotypes_array[i]
    
    cat("--------------------------->\t")
    cat(sprintf(as.character(phenotypes_array_sel)))
    cat("\n")
    
    
    ALL_dB_subset_restricted_sel<-ALL_dB_subset_restricted[which(ALL_dB_subset_restricted$phenotype%in%phenotypes_array_sel),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("ALL_dB_subset_restricted_sel_0\n")
      cat(str(ALL_dB_subset_restricted_sel))
      cat("\n")
      cat(str(unique(ALL_dB_subset_restricted_sel$VAR)))
      cat("\n")
    }
    
    #### Merge with genes with consequences ----
    
    ALL_dB_subset_restricted_sel<-merge(ALL_dB_subset_restricted_sel,
                                        VEP_and_PCHiC,
                                        by="VAR")
    
    
    if(CONDITION_DEBUG == 1)
    {
      cat("ALL_dB_subset_restricted_sel_1\n")
      cat(str(ALL_dB_subset_restricted_sel))
      cat("\n")
      cat(str(unique(ALL_dB_subset_restricted_sel$VAR)))
      cat("\n")
    }
    
    VEP_and_PCHiC_sel<-VEP_and_PCHiC[which(VEP_and_PCHiC$VAR%in%ALL_dB_subset_restricted_sel$VAR),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("VEP_and_PCHiC_sel_1\n")
      cat(str(VEP_and_PCHiC_sel))
      cat("\n")
      cat(str(unique(VEP_and_PCHiC_sel$VAR)))
      cat("\n")
      cat(str(unique(VEP_and_PCHiC_sel$ensembl_gene_id)))
      cat("\n")
    }
    
    
    TOME_correspondence_sel<-TOME_correspondence[which(TOME_correspondence$phenotype == phenotypes_array_sel),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("TOME_correspondence_sel_0\n")
      cat(str(TOME_correspondence_sel))
      cat("\n")
    }
    
    GENE_EXP_sel<-GENE_EXP[which(GENE_EXP$ensembl_gene_id%in%VEP_and_PCHiC_sel$ensembl_gene_id),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("GENE_EXP_sel_1\n")
      cat(str(GENE_EXP_sel))
      cat("\n")
      cat(str(unique(GENE_EXP_sel$ensembl_gene_id)))
      cat("\n")
    }
    
    
    indx.int<-c(which(colnames(GENE_EXP_sel) == "ensembl_gene_id"),grep("mean_",colnames(GENE_EXP_sel)))
    
    GENE_EXP_sel_subset<-unique(GENE_EXP_sel[,indx.int])
    
    if(CONDITION_DEBUG == 1)
    {
      cat("GENE_EXP_sel_subset_0\n")
      cat(str(GENE_EXP_sel_subset))
      cat("\n")
      cat(str(unique(GENE_EXP_sel_subset$ensembl_gene_id)))
      cat("\n")
    }
    
    GENE_EXP_sel_subset.m<-melt(GENE_EXP_sel_subset, id.vars="ensembl_gene_id", variable.name= "Cell_Type",value.name = "mean_GENE_EXP")
    
    if(CONDITION_DEBUG == 1)
    {
      cat("GENE_EXP_sel_subset.m_0\n")
      cat(str(GENE_EXP_sel_subset.m))
      cat("\n")
      cat(str(unique(GENE_EXP_sel_subset.m$ensembl_gene_id)))
      cat("\n")
    }
    
    GENE_EXP_sel_subset.m$Cell_Type<-gsub("mean_","",GENE_EXP_sel_subset.m$Cell_Type)
    GENE_EXP_sel_subset.m$Cell_Type<-gsub("\\.","-",GENE_EXP_sel_subset.m$Cell_Type)
    
    if(CONDITION_DEBUG == 1)
    {
      cat("GENE_EXP_Cell_Types_0\n")
      cat(sprintf(as.character(unique(GENE_EXP_sel_subset.m$Cell_Type))))
      cat("\n")
    }
    
    if(CONDITION_DEBUG == 1)
    {
      cat("TOME_Cell_Type_0\n")
      cat(sprintf(as.character(unique(TOME_correspondence$TOME_Cell_Type))))
      cat("\n")
      cat(sprintf(as.character(unique(TOME_correspondence_sel$TOME_Cell_Type))))
      cat("\n")
    }
    
    
    GENE_EXP_sel_subset.m$Tag<-NA
    
    GENE_EXP_sel_subset.m$Tag[which(GENE_EXP_sel_subset.m$Cell_Type%in%TOME_correspondence$TOME_Cell_Type)]<-"Not_relevant"
    GENE_EXP_sel_subset.m$Tag[which(GENE_EXP_sel_subset.m$Cell_Type%in%TOME_correspondence_sel$TOME_Cell_Type)]<-"Relevant"
    
    if(CONDITION_DEBUG == 1)
    {
      cat("GENE_EXP_sel_subset.m_Tag\n")
      cat(sprintf(as.character(unique(GENE_EXP_sel_subset.m$Cell_Type))))
      cat("\n")
      cat(sprintf(as.character(unique(GENE_EXP_sel_subset.m$Tag))))
      cat("\n")
      cat("NA_CT\n")
      cat(sprintf(as.character(unique(GENE_EXP_sel_subset.m$Cell_Type[is.na(GENE_EXP_sel_subset.m$Tag)]))))
      cat("\n")
    }
    
    GENE_EXP_sel_subset.m_NO_NA<-GENE_EXP_sel_subset.m[!is.na(GENE_EXP_sel_subset.m$Tag),]
    GENE_EXP_sel_subset.m_NO_NA$phenotype<-phenotypes_array_sel
    
    if(CONDITION_DEBUG == 1)
    {
      cat("GENE_EXP_sel_subset.m_NO_NA_0\n")
      cat(str(GENE_EXP_sel_subset.m_NO_NA))
      cat("\n")
      cat(str(unique(GENE_EXP_sel_subset.m_NO_NA$ensembl_gene_id)))
      cat("\n")
      cat(str(unique(GENE_EXP_sel_subset.m_NO_NA$phenotype)))
      cat("\n")
    }
    
    GENE_EXP_sel_subset.m_NO_NA<-merge(VEP_and_PCHiC_sel,
                                       GENE_EXP_sel_subset.m_NO_NA,
                                       by="ensembl_gene_id")
    
    if(CONDITION_DEBUG == 1)
    {
      cat("GENE_EXP_sel_subset.m_NO_NA_1\n")
      cat(str(GENE_EXP_sel_subset.m_NO_NA))
      cat("\n")
      cat(str(unique(GENE_EXP_sel_subset.m_NO_NA$ensembl_gene_id)))
      cat("\n")
      cat(str(unique(GENE_EXP_sel_subset.m_NO_NA$phenotype)))
      cat("\n")
    }
    
    GENE_EXP_sel_subset.m_NO_NA_wide<-as.data.frame(pivot_wider(GENE_EXP_sel_subset.m_NO_NA, id_cols=c("VAR","ensembl_gene_id","HGNC","phenotype"),
                                                  names_from=Cell_Type,
                                                  values_from=c("mean_GENE_EXP","Tag")), stringsAsFactors=F)
    
    if(CONDITION_DEBUG == 1)
    {
      cat("GENE_EXP_sel_subset.m_NO_NA_wide_0\n")
      cat(str(GENE_EXP_sel_subset.m_NO_NA_wide))
      cat("\n")
      cat(str(unique(GENE_EXP_sel_subset.m_NO_NA_wide$VAR)))
      cat("\n")
      cat(str(unique(GENE_EXP_sel_subset.m_NO_NA_wide$ensembl_gene_id)))
      cat("\n")
      cat(str(unique(GENE_EXP_sel_subset.m_NO_NA_wide$phenotype)))
      cat("\n")
    }
    
    Gather<-unique(rbind(GENE_EXP_sel_subset.m_NO_NA_wide,Gather))
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Gather_0\n")
      cat(str(Gather))
      cat("\n")
      cat(str(unique(Gather$VAR)))
      cat("\n")
      cat(str(unique(Gather$ensembl_gene_id)))
      cat("\n")
      cat(str(unique(Gather$phenotype)))
      cat("\n")
    }
    
    
    # ##########################################################
    # quit(status = 1)
    
    
    
  }# i in 1:length(phenotypes_array
  
  
  
  
  
 
  
  check<-Gather[which(Gather$VAR%in%tracking_variants),]
  
  if(CONDITION_DEBUG == 1)
  {
    cat("check_0\n")
    cat(str(check))
    cat("\n")
    cat(str(unique(check$VAR)))
    cat("\n")
  }
  
  # #### SAVE ----
  
  setwd(out)
  
  write.table(Gather, file="GENE_EXP_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
  write.table(check, file="check_GENE_EXP_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
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
    make_option(c("--GENE_EXP"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TOME_correspondence"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VEP_CSQ"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--PCHiC"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--excluded_phenotypes"), type="character", default=NULL, 
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
  
  Data_wrangling_GENE_EXP(opt)
  
}


###########################################################################

system.time( main() )
