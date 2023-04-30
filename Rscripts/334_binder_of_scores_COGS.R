
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
  
  ### Retrieve COGS ----
  
  
  master_path<-"/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/COGS/rCOGS_out/COGS_scores/"
  
  df<-data.frame()
  
  if(file.exists(master_path))
  {
    setwd(master_path)
    
    file_list <- list.files(path=master_path, include.dirs = FALSE)
    
    cat("file_list_0\n")
    cat(str(file_list))
    cat("\n")
    
    file_list_sel<-file_list[grep("_interest_cell_types_",file_list)]
    
    cat("file_list_sel_0\n")
    cat(str(file_list_sel))
    cat("\n")
    
    phenotypes<-gsub("_COGS_interest_cell_types_with_gene_names","", file_list_sel)
    
    cat("phenotypes_0\n")
    cat(str(phenotypes))
    cat("\n")
    cat(sprintf(as.character(phenotypes)))
    cat("\n")
    
    df<-as.data.frame(cbind(phenotypes,file_list_sel), stringsAsFactors=F)
    
    colnames(df)<-c("phenotype","file")
    
    cat("df_0\n")
    cat(str(df))
    cat("\n")
    
    ### Exclude percentage phenotypes of white cells -----
    
    df_restricted<-unique(df[-which(df$phenotype%in%excluded_phenotypes),])
    
    
    cat("df_restricted_0\n")
    cat(str(df_restricted))
    cat("\n")
   
    
    phenotypes_array<-unique(df_restricted$phenotype)
    
    cat("phenotypes_array_0\n")
    cat(str(phenotypes_array))
    cat("\n")
    
    
    Gather<-data.frame()
    
    for(i in 1:length(phenotypes_array))
    {
      phenotypes_array_sel<-phenotypes_array[i]
      
      cat("-------------------------------------------------------------------------------->\t")
      cat(sprintf(as.character(phenotypes_array_sel)))
      cat("\n")
      
      df_restricted_sel<-df_restricted[which(df_restricted$phenotype == phenotypes_array_sel),]
      
      cat("df_restricted_sel_0\n")
      cat(str(df_restricted_sel))
      cat("\n")
      
      file_selected_to_open<-unique(df_restricted_sel$file)
      
      cat("file_selected_to_open_0\n")
      cat(str(file_selected_to_open))
      cat("\n")
      
      if(file.exists(file_selected_to_open))
      {
        COGS<-as.data.frame(fread(file=file_selected_to_open,sep="\t") , stringsAsFactors=F)
        
        colnames(COGS)[which(colnames(COGS) == "ensg")]<-"ensembl_gene_id"
        
        cat("COGS\n")
        cat(str(COGS))
        cat("\n")
        
        indx.int<-c(which(colnames(COGS) == "ensembl_gene_id"),
                    which(colnames(COGS) == "cogs"))
        
        COGS_subset<-unique(COGS[,indx.int])
        
        cat("COGS_subset_0\n")
        cat(str(COGS_subset))
        cat("\n")
       
        ALL_dB_subset_restricted_sel<-ALL_dB_subset_restricted[which(ALL_dB_subset_restricted$phenotype == phenotypes_array_sel),]
        
        cat("ALL_dB_subset_restricted_sel_0\n")
        cat(str(ALL_dB_subset_restricted_sel))
        cat("\n")
        cat(str(unique(ALL_dB_subset_restricted_sel$VAR)))
        cat("\n")
       
        VEP_and_PCHiC_sel<-VEP_and_PCHiC[which(VEP_and_PCHiC$VAR%in%ALL_dB_subset_restricted_sel$VAR),]
        
        cat("VEP_and_PCHiC_sel_0\n")
        cat(str(VEP_and_PCHiC_sel))
        cat("\n")
        cat(str(unique(VEP_and_PCHiC_sel$VAR)))
        cat("\n")
        
        COGS_subset<-merge(VEP_and_PCHiC_sel,
                           COGS_subset,
                           by="ensembl_gene_id")
        
        COGS_subset$phenotype<-phenotypes_array_sel
        
        cat("COGS_subset_0\n")
        cat(str(COGS_subset))
        cat("\n")
        cat(str(unique(COGS_subset$VAR)))
        cat("\n")
        
        Gather<-rbind(COGS_subset,Gather)
        
        # cat("df_0\n")
        # cat(str(df))
        # cat("\n")
        # cat(str(unique(df$VAR)))
        # cat("\n")
        
        # #########################################################
        # quit(status = 1)
        
      }# file.exists(file_selected_to_open)
      
      
    }# i in 1:length(phenotypes_array)
    
  }# file.exists(master_path)
  
  
  
  CONDITION_DEBUG<-1
  
  if(CONDITION_DEBUG == 1)
  {
    cat("Gather_0\n")
    cat(str(Gather))
    cat("\n")
    cat(str(unique(Gather$VAR)))
    cat("\n")
  }
  
  check<-Gather[which(Gather$VAR%in%tracking_variants),]
  
  if(CONDITION_DEBUG == 1)
  {
    cat("check_0\n")
    cat(str(check))
    cat("\n")
    cat(str(unique(check$VAR)))
    cat("\n")
  }
  
  #### exclude NA's ----
  
  Gather_NO_NA<-Gather[!is.na(Gather$cogs),]
  
  
  cat("Gather_NO_NA_0\n")
  cat(str(Gather_NO_NA))
  cat("\n")
  cat(str(unique(Gather_NO_NA$VAR)))
  cat("\n")
  
  # #### SAVE ----
  
  setwd(out)
  
  write.table(Gather_NO_NA, file="COGS_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
  write.table(check, file="check_COGS_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
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
  
  Data_wrangling_ATAC(opt)
  
}


###########################################################################

system.time( main() )
