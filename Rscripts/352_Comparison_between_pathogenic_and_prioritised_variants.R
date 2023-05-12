
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
  
  #### READ and transform finemap_prob_Threshold ----
  
  finemap_prob_Threshold = opt$finemap_prob_Threshold
  
  cat("finemap_prob_Threshold_\n")
  cat(sprintf(as.character(finemap_prob_Threshold)))
  cat("\n")
 

  #### RMV DGKQ ----
  
  
  RMV_common = opt$RMV_common
  
  cat("RMV_common_\n")
  cat(sprintf(as.character(RMV_common)))
  cat("\n")
  
  #### RMV DGKQ ----
  
  
  RMV_labels = unlist(strsplit(opt$RMV_labels, split=","))
  
  cat("RMV_labels_\n")
  cat(sprintf(as.character(RMV_labels)))
  cat("\n")
  
 
  #### Read GWAS_GLOBAL_per_traits -----
  
  setwd(out)
  
  GWAS_GLOBAL_per_traits<-as.data.frame(fread(file=paste("GWAS_GLOBAL_per_traits_Thresholded",'_',finemap_prob_Threshold,'.tsv', sep=''), sep="\t", header=T), stringsAsFactors=F)
  
  cat("GWAS_GLOBAL_per_traits_0\n")
  cat(str(GWAS_GLOBAL_per_traits))
  cat("\n")
  cat(str(unique(GWAS_GLOBAL_per_traits$VAR)))
  cat("\n")
  cat(sprintf(as.character(unique(GWAS_GLOBAL_per_traits$variable))))
  cat("\n")
  
  
  
  
  ### Read List_of_pathogenic_variants----
  
  
  List_of_pathogenic_variants<-as.data.frame(fread(file=opt$List_of_pathogenic_variants, sep="\t", header=T) , stringsAsFactors=F)
  
  cat("List_of_pathogenic_variants_0\n")
  cat(str(List_of_pathogenic_variants))
  cat("\n")
  cat(str(unique(List_of_pathogenic_variants$rs)))
  cat("\n")
  
  
  ### Read VAR_Prioritization_dB----
  
  
  VAR_Prioritization_dB<-as.data.frame(readRDS(file=opt$VAR_Prioritization_dB) , stringsAsFactors=F)
  
  cat("VAR_Prioritization_dB_0\n")
  cat(str(VAR_Prioritization_dB))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_Prioritization_dB$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_Prioritization_dB$Fig1_Annot_Category))))
  cat("\n")
  
  
  check<-VAR_Prioritization_dB[which(VAR_Prioritization_dB$rs%in%List_of_pathogenic_variants$rs),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  cat(str(unique(check$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$Fig1_Annot_Category))))
  cat("\n")
  
  
  VAR_Prioritization_dB$category_INTRO<-NA
  
  VAR_Prioritization_dB$category_INTRO<-as.character(VAR_Prioritization_dB$Fig1_Annot_Category)
  VAR_Prioritization_dB$category_INTRO[which(VAR_Prioritization_dB$rs%in%List_of_pathogenic_variants$rs)]<-"Pathogenic_variant"
  
  VAR_Prioritization_dB$category_INTRO<-factor(VAR_Prioritization_dB$category_INTRO,
                                               levels=c("Common_variant","RV_C","RV_NC_lowPP","RV_NC_highPP_lowEffectSize","RV_NC_highPP_highEffectSize","Pathogenic_variant"),
                                               ordered=T)
  
  
  cat("VAR_Prioritization_dB_1\n")
  cat(str(VAR_Prioritization_dB))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(VAR_Prioritization_dB$category_INTRO))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(VAR_Prioritization_dB$category_INTRO)))))
  cat("\n")
  
  #### Merge ----
  
  
  GWAS_GLOBAL_per_traits<-merge(GWAS_GLOBAL_per_traits,
                                VAR_Prioritization_dB,
                                by="VAR")
  
  
  cat("GWAS_GLOBAL_per_traits_1\n")
  cat(str(GWAS_GLOBAL_per_traits))
  cat("\n")
  cat(str(unique(GWAS_GLOBAL_per_traits$VAR)))
  cat("\n")
  
  GWAS_GLOBAL_per_traits_subset<-droplevels(GWAS_GLOBAL_per_traits[which(GWAS_GLOBAL_per_traits$variable == 'Absolute_effect_size'),])
  
  cat("GWAS_GLOBAL_per_traits_subset_0\n")
  cat(str(GWAS_GLOBAL_per_traits_subset))
  cat("\n")
  cat(str(unique(GWAS_GLOBAL_per_traits_subset$VAR)))
  cat("\n")
  
  GWAS_GLOBAL_per_traits_subset.dt<-data.table(GWAS_GLOBAL_per_traits_subset, key=c("category_INTRO","variable"))
  
  Summary_table<-as.data.frame(GWAS_GLOBAL_per_traits_subset.dt[,.(n=.N,
                                                                  Min=round(as.numeric(summary(value_Z_score)[1]),3),
                                                                  Q1=round(as.numeric(summary(value_Z_score)[2]),3),
                                                                  M=round(as.numeric(summary(value_Z_score)[3]),3),
                                                                  Q3=round(as.numeric(summary(value_Z_score)[5]),3),
                                                                  Max=round(as.numeric(summary(value_Z_score)[6]),3)),by=key(GWAS_GLOBAL_per_traits_subset.dt)], stringsAsFactors=F)
  
  colnames(Summary_table)[which(colnames(Summary_table) == 'category_INTRO')]<-'category2'
  
  cat("Summary_table\n")
  cat(str(Summary_table))
  cat("\n")
  
  CONDITION_DEBUG<-1
  
  
  indx.category_INTRO<-which(colnames(GWAS_GLOBAL_per_traits_subset)=='category_INTRO')
  
  levels_category_vector<-levels(GWAS_GLOBAL_per_traits_subset[,indx.category_INTRO])
  
  if(CONDITION_DEBUG == 1)
  {
    cat("levels_category_vector\n")
    cat(str(levels_category_vector))
    cat("\n")
    
  }
  
  Pivotal_variety<-levels(droplevels(GWAS_GLOBAL_per_traits_subset[,indx.category_INTRO]))
  
  
  
  if(CONDITION_DEBUG == 1)
  {
    cat("Pivotal_variety\n")
    cat(str(Pivotal_variety))
    cat("\n")
    
  }
  
  
  List_variables<-list()
  
  if(length(Pivotal_variety) >1)
  {
    ### Pair wilcox test between levels of the category
    
    PW_category<-pairwise.wilcox.test(GWAS_GLOBAL_per_traits_subset$value_Z_score, GWAS_GLOBAL_per_traits_subset[,indx.category_INTRO],
                                      p.adjust.method = "BH")
    
    
    if(CONDITION_DEBUG == 1)
    {
      cat("PW_category\n")
      cat(str(PW_category))
      cat("\n")
    }
    
    PW_category_pvalue_df<-as.data.frame(PW_category$p.value, stringsAsFactors=F)
    
    list_cols<-list()
    
    
    
    for(PW_iteration in 1:dim(PW_category_pvalue_df)[2])
    {
      colnames_sel<-colnames(PW_category_pvalue_df)[PW_iteration]
      
      if(CONDITION_DEBUG == 1)
      {
        cat("----------------->colnames_sel\n")
        cat(sprintf(as.character(colnames_sel)))
        cat("\n")
      }
      
      list_rows<-list()
      
      
      for(PW_iteration_k in 1:dim(PW_category_pvalue_df)[1])
      {
        rownames_sel<-row.names(PW_category_pvalue_df)[PW_iteration_k]
        
        if(CONDITION_DEBUG == 1)
        {
          cat("--->rownames_sel\n")
          cat(sprintf(as.character(rownames_sel)))
          cat("\n")
        }
        
        PW_Wilcox_pvalue<-PW_category_pvalue_df[PW_iteration_k,PW_iteration]
        
        if(CONDITION_DEBUG == 1)
        {
          cat("PW_Wilcox_pvalue\n")
          cat(sprintf(as.character(PW_Wilcox_pvalue)))
          cat("\n")
        }
        
        
        log_pval_PW_Wilcox<-round(-1*log10(PW_Wilcox_pvalue),4)
        
        if(CONDITION_DEBUG == 1)
        {
          cat("log_pval_PW_Wilcox\n")
          cat(sprintf(as.character(log_pval_PW_Wilcox)))
          cat("\n")
        }
        
        
        vector_final_comparisons<-paste(sort(c(colnames_sel,rownames_sel)), collapse=";")
        
        
        a.dt<-as.data.frame(cbind(vector_final_comparisons,PW_Wilcox_pvalue,log_pval_PW_Wilcox), stringsAsFactors=F)
        
        colnames(a.dt)<-c("string_comp",'pval','MINUS_logpval')
        
        if(CONDITION_DEBUG == 1)
        {
          cat("a.dt\n")
          cat(str(a.dt))
          cat("\n")
        }
        
        list_rows[[PW_iteration_k]]<-a.dt
        
      }#PW_iteration_k
      
      df_col = as.data.frame(data.table::rbindlist(list_rows, fill=T), stringsAsFactors=F)
      
      if(CONDITION_DEBUG == 1)
      {
        cat("df_col\n")
        cat(str(df_col))
        cat("\n")
      }
      
      list_cols[[PW_iteration]]<-df_col
      
      
    }#PW_iteration
    
    PW_category = as.data.frame(data.table::rbindlist(list_cols, fill=T), stringsAsFactors=F)
    
    PW_category[,which(colnames(PW_category) == 'pval')]<-as.numeric(PW_category[,which(colnames(PW_category) == 'pval')])
    PW_category[,which(colnames(PW_category) == 'MINUS_logpval')]<-as.numeric(PW_category[,which(colnames(PW_category) == 'MINUS_logpval')])
    
    if(CONDITION_DEBUG == 1)
    {
      cat("PW_category\n")
      cat(str(PW_category))
      cat("\n")
    }
    
    
    PW_category_NO_NA<-PW_category[!is.na(PW_category[,which(colnames(PW_category) == 'MINUS_logpval')]),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("PW_category_NO_NA_0\n")
      cat(str(PW_category_NO_NA))
      cat("\n")
    }
    
    if(dim(PW_category_NO_NA)[1] > 0)
    {
      PW_category_NO_NA$category1<-gsub(";.+$","",PW_category_NO_NA$string_comp)
      PW_category_NO_NA$category2<-gsub("^[^;]+;","",PW_category_NO_NA$string_comp)
      
      if(CONDITION_DEBUG == 1)
      {
        cat("PW_category_NO_NA_1\n")
        cat(str(PW_category_NO_NA))
        cat("\n")
      }
      
      PW_category_NO_NA<-PW_category_NO_NA[,c(which(colnames(PW_category_NO_NA) == "category2"),
                                              which(colnames(PW_category_NO_NA) == "category1"),
                                              which(colnames(PW_category_NO_NA) == 'pval'),
                                              which(colnames(PW_category_NO_NA) == 'MINUS_logpval'))]
      
      if(CONDITION_DEBUG == 1)
      {
        cat("PW_category_NO_NA_2\n")
        cat(str(PW_category_NO_NA))
        cat("\n")
      }
      
      PW_category_NO_NA$category1<-factor(PW_category_NO_NA$category1,
                                          levels=levels_category_vector,
                                          ordered=T)
      
      PW_category_NO_NA$category2<-factor(PW_category_NO_NA$category2,
                                          levels=levels_category_vector,
                                          ordered=T)
      
      
      
      PW_category_NO_NA$variable<-'Absolute_effect_size'
      PW_category_NO_NA$Annotation<-'category_INTRO'
      
      if(CONDITION_DEBUG == 1)
      {
        cat("PW_category_NO_NA_3\n")
        cat(str(PW_category_NO_NA))
        cat("\n")
      }
      
      
      
      PW_category_NO_NA<-merge(PW_category_NO_NA,
                               Summary_table,
                               by=c("category2","variable"))
      
      colnames(Summary_table)[which(colnames(Summary_table) == "category2")]<-"category1"
      
      PW_category_NO_NA<-merge(PW_category_NO_NA,
                               Summary_table,
                               by=c("category1","variable"),
                               all.x=T)
      
      if(CONDITION_DEBUG == 1)
      {
        cat("PW_category_NO_NA_4\n")
        cat(str(PW_category_NO_NA))
        cat("\n")
      }
      
      
      PW_category_NO_NA<-PW_category_NO_NA[order(PW_category_NO_NA$category2,PW_category_NO_NA$category1),]
      
      colnames(PW_category_NO_NA)<-gsub("\\.x","_category2",colnames(PW_category_NO_NA))
      colnames(PW_category_NO_NA)<-gsub("\\.y","_category1",colnames(PW_category_NO_NA))
      
      if(CONDITION_DEBUG == 1)
      {
        cat("PW_category_NO_NA_5\n")
        cat(str(PW_category_NO_NA))
        cat("\n")
      }
      
      PW_category_NO_NA<-PW_category_NO_NA[,c("category2","category1","variable",
                                              "pval","MINUS_logpval",
                                              "n_category2","n_category1",
                                              "M_category2","M_category1",
                                              "Min_category2","Min_category1",
                                              "Q1_category2","Q1_category1",
                                              "Q3_category2","Q3_category1",
                                              "Max_category2","Max_category1","Annotation")]
      
      
      
      if(CONDITION_DEBUG == 1)
      {
        cat("PW_category_NO_NA_6\n")
        cat(str(PW_category_NO_NA))
        cat("\n")
      }
      
      ############## SAVE -------------------------
      
      setwd(out)
      
      write.table(PW_category_NO_NA, file=paste("Pathogenic_Comparison_Absolute_Effect_size_Z_score_with_stats",'_',finemap_prob_Threshold,'.tsv', sep=''), 
                  sep="\t",quote=F, row.names = F)
      
      
      
      
    }# dim(PW_category_NO_NA)[1] > 0
  }#length(Pivotal_variety) >1
 
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
    make_option(c("--Table_S4"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--finemap_prob_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--List_of_pathogenic_variants"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VAR_Prioritization_dB"), type="character", default=NULL, 
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
