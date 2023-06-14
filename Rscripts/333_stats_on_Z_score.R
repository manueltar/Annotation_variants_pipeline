
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

Stats_function = function(option_list)
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
  
  #### Read Master_file ----
  
  Master_file<-as.data.frame(readRDS(file=opt$Master_file) , stringsAsFactors=F)
  
  cat("Master_file_0\n")
  cat(str(Master_file))
  cat("\n")
  cat(str(unique(Master_file$VAR)))
  cat("\n")
  
  #### Read Table_of_labels ----
  
  Table_of_labels<-as.data.frame(readRDS(file=opt$Table_of_labels) , stringsAsFactors=F)
  
  cat("Table_of_labels_0\n")
  cat(str(Table_of_labels))
  cat("\n")
  cat(str(unique(Table_of_labels$VAR)))
  cat("\n")
 
  cat("Table_of_labels_0\n")
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
  
  ############## STATS LOOP -----------------
  
  
  ####   LOOP Through variables ----
  
  category_vector<-c("Fig1_Annot_Category","MPRA_CLASS_2","M_and_M","Multi_Lineage")
  
  
  List_DEF<-list()
  
  CONDITION_DEBUG<-1
  
  for(i in 1:length(category_vector))
  {
    category_sel<-category_vector[i]
    
    cat("-------------------------------->\t")
    cat(sprintf(as.character(category_sel)))
    cat("\n")
    
    indx.category_sel<-which(colnames(Table_of_labels) == category_sel)
    
    if(CONDITION_DEBUG == 1)
    {
      cat("indx.category_sel\n")
      cat(str(indx.category_sel))
      cat("\n")
      
    }
    
    ind.int<-c(which(colnames(Table_of_labels) == "VAR"),indx.category_sel)
    
    if(category_sel == "Multi_Lineage")
    {
      ind.int<-c(which(colnames(Table_of_labels) == "VAR"),which(colnames(Table_of_labels) == "M_and_M"),indx.category_sel)
      
    }#category_sel == "Multi_Lineage")
    
    Table_of_labels_subset<-unique(Table_of_labels[,ind.int])
    
    indx.category_sel_subset<-which(colnames(Table_of_labels_subset) == category_sel)
    
    if(CONDITION_DEBUG == 1)
    {
      cat("indx.category_sel_subset\n")
      cat(str(indx.category_sel_subset))
      cat("\n")
      
    }
    
    
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Table_of_labels_subset_0\n")
      cat(str(Table_of_labels_subset))
      cat("\n")
      cat(str(unique(Table_of_labels_subset$VAR)))
      cat("\n")
      cat(sprintf(as.character(names(summary(Table_of_labels_subset[,indx.category_sel_subset])))))
      cat("\n")
      cat(sprintf(as.character(summary(Table_of_labels_subset[,indx.category_sel_subset]))))
      cat("\n")
    }
    
   
    Table_of_labels_subset_restricted<-Table_of_labels_subset[!is.na(Table_of_labels_subset[,indx.category_sel_subset]),]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Table_of_labels_subset_restricted_0\n")
      cat(str(Table_of_labels_subset_restricted))
      cat("\n")
      cat(str(unique(Table_of_labels_subset_restricted$VAR)))
      cat("\n")
      cat(sprintf(as.character(names(summary(Table_of_labels_subset_restricted[,indx.category_sel_subset])))))
      cat("\n")
      cat(sprintf(as.character(summary(Table_of_labels_subset_restricted[,indx.category_sel_subset]))))
      cat("\n")
    }
    
    
    array_variable<-levels(Master_file$variable)
    
    if(CONDITION_DEBUG == 1)
    {
      cat("array_variable\n")
      cat(str(array_variable))
      cat("\n")
      
    }
    
    if(category_sel == "Multi_Lineage")
    {
      Figure_4_categories_of_M_and_M<-c("Transcriptional_Regulation|R_in_candidate","Alternative_Transcript_Usage|R_in_candidate","Transcriptional_Regulation_and_ATU|R_in_candidate")
      
      Table_of_labels_subset_restricted<-Table_of_labels_subset_restricted[which(Table_of_labels_subset_restricted$M_and_M%in%Figure_4_categories_of_M_and_M),]
      
      if(CONDITION_DEBUG == 1)
      {
        cat("Table_of_labels_subset_restricted_SPECIAL\n")
        cat(str(Table_of_labels_subset_restricted))
        cat("\n")
        cat(str(unique(Table_of_labels_subset_restricted$VAR)))
        cat("\n")
        cat(sprintf(as.character(names(summary(Table_of_labels_subset_restricted[,indx.category_sel_subset])))))
        cat("\n")
        cat(sprintf(as.character(summary(Table_of_labels_subset_restricted[,indx.category_sel_subset]))))
        cat("\n")
      }
      
      
     # quit(status = 1) 
    }else{
      
      if(category_sel == "M_and_M")
      {
        Figure_4_categories_of_M_and_M<-c('Transcriptional_Regulation|R_in_candidate','No_Regulation_Detected|Unknown_mechanism','No_Regulation_Detected|Different_RNAseq_readout_needed')
        
        Table_of_labels_subset_restricted<-Table_of_labels_subset_restricted[which(Table_of_labels_subset_restricted$M_and_M%in%Figure_4_categories_of_M_and_M),]
        
        if(CONDITION_DEBUG == 1)
        {
          cat("Table_of_labels_subset_restricted_SPECIAL\n")
          cat(str(Table_of_labels_subset_restricted))
          cat("\n")
          cat(str(unique(Table_of_labels_subset_restricted$VAR)))
          cat("\n")
          cat(sprintf(as.character(names(summary(Table_of_labels_subset_restricted[,indx.category_sel_subset])))))
          cat("\n")
          cat(sprintf(as.character(summary(Table_of_labels_subset_restricted[,indx.category_sel_subset]))))
          cat("\n")
        }
      }else{
        
        # Do nothing
        
      }#category_sel == "M_and_M"
    }
      
    
    
    List_variables<-list()
    
    for(k in 1:length(array_variable))
    {
      
      array_variable_sel<-array_variable[k]
      
      cat("-------------------------------->\t")
      cat(sprintf(as.character(array_variable_sel)))
      cat("\n")
      
      Master_file_sel<-Master_file[which(Master_file$variable == array_variable_sel),]
      
      if(CONDITION_DEBUG == 1)
      {
        cat("Master_file_sel_0\n")
        cat(str(Master_file_sel))
        cat("\n")
        cat(str(unique(Master_file_sel$VAR)))
        cat("\n")
        cat(sprintf(as.character(names(summary(Master_file_sel$value_Z_score)))))
        cat("\n")
        cat(sprintf(as.character(summary(Master_file_sel$value_Z_score))))
        cat("\n")
      }
      
      
      ### Add labels
      
      Master_file_sel<-merge(Master_file_sel,
                             Table_of_labels_subset_restricted,
                             by="VAR")
      
      if(CONDITION_DEBUG == 1)
      {
        cat("Master_file_sel_1\n")
        cat(str(Master_file_sel))
        cat("\n")
        cat(str(unique(Master_file_sel$VAR)))
        cat("\n")
        cat(sprintf(as.character(names(summary(Master_file_sel$value_Z_score)))))
        cat("\n")
        cat(sprintf(as.character(summary(Master_file_sel$value_Z_score))))
        cat("\n")
      }
      
      
      if(dim(Master_file_sel)[1] >0)
      {
        
       
        
        ### calculate distribution parameters by label
        
        indx.category_Master_sel<-which(colnames(Master_file_sel) == category_sel)
        
        if(CONDITION_DEBUG == 1)
        {
          cat("indx.category_Master_sel\n")
          cat(str(indx.category_Master_sel))
          cat("\n")
          
        }
        
        levels_category_vector<-levels(Master_file_sel[,indx.category_Master_sel])
        
        if(CONDITION_DEBUG == 1)
        {
          cat("levels_category_vector\n")
          cat(str(levels_category_vector))
          cat("\n")
          
        }
        
        Master_file_sel.dt<-data.table(Master_file_sel, key=colnames(Master_file_sel)[indx.category_Master_sel])
        
        # if(CONDITION_DEBUG == 1)
        # {
        #   cat("Master_file_sel.dt_1\n")
        #   cat(str(Master_file_sel.dt))
        #   cat("\n")
        # }
        
        
        Summary_table<-as.data.frame(Master_file_sel.dt[,.(n=.N,
                                                           Min=round(as.numeric(summary(value_Z_score)[1]),3),
                                                           Q1=round(as.numeric(summary(value_Z_score)[2]),3),
                                                           M=round(as.numeric(summary(value_Z_score)[3]),3),
                                                           Q3=round(as.numeric(summary(value_Z_score)[5]),3),
                                                           Max=round(as.numeric(summary(value_Z_score)[6]),3)),, by=key(Master_file_sel.dt)], stringsAsFactors=F)
        
        colnames(Summary_table)[which(colnames(Summary_table) == colnames(Master_file_sel)[indx.category_Master_sel])]<-"category2"
        Summary_table$variable<-array_variable_sel
        
        if(CONDITION_DEBUG == 1)
        {
          cat("Summary_table_1\n")
          cat(str(Summary_table))
          cat("\n")
        }
        
        
        
        Pivotal_variety<-levels(droplevels(Master_file_sel[,indx.category_Master_sel]))
        
        if(CONDITION_DEBUG == 1)
        {
          cat("Pivotal_variety\n")
          cat(str(Pivotal_variety))
          cat("\n")
          
        }
        
        CONDITION_DEBUG<-0
        
        if(length(Pivotal_variety) >1)
        {
          ### Pair wilcox test between levels of the category
          
          PW_category<-pairwise.wilcox.test(Master_file_sel$value, Master_file_sel[,indx.category_Master_sel],
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
            
            
            
            PW_category_NO_NA$variable<-array_variable_sel
            
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
            
            
            List_variables[[k]]<-PW_category_NO_NA
            
            # ##########################################################################
            # quit(status = 1)
            
          }# dim(PW_category_NO_NA)[1] > 0
        }#length(Pivotal_variety) >1
        
      }# dim(Master_file_sel)[1] >0
      
      CONDITION_DEBUG<-1
    }# k in 1:array_variable
    
   
    CONDITION_DEBUG<-1
    
    if(length(List_variables) >1 )
    {
      Wilcox_df_3 = unique(as.data.frame(data.table::rbindlist(List_variables, fill = T)))
      
      Wilcox_df_3$Annotation<-category_sel
      
      if(CONDITION_DEBUG == 1)
      {
        cat("Wilcox_df_3_0\n")
        cat(str(Wilcox_df_3))
        cat("\n")
      }
      
      Wilcox_df_3<-Wilcox_df_3[,c("category2","category1","variable",
                                  "pval","MINUS_logpval",
                                  "n_category2","n_category1",
                                  "M_category2","M_category1",
                                  "Min_category2","Min_category1",
                                  "Q1_category2","Q1_category1",
                                  "Q3_category2","Q3_category1",
                                  "Max_category2","Max_category1","Annotation")]
      
      if(category_sel == "Fig1_Annot_Category")
      {

        Wilcox_df_3<-Wilcox_df_3[which(Wilcox_df_3$category1 == "index_variants" |
                                                             Wilcox_df_3$category2 == "index_variants"),]

      }
      
      if(CONDITION_DEBUG == 1)
      {
        cat("Wilcox_df_3_1\n")
        cat(str(Wilcox_df_3))
        cat("\n")
      }
      List_DEF[[i]]<-Wilcox_df_3
      
    }# length(List_variables) >1
  }#i category_sel
  
  if(length(List_DEF) >1 )
  {
    Wilcox_df_4 = unique(as.data.frame(data.table::rbindlist(List_DEF, fill = T)))
    
       
    cat("Wilcox_df_4_0\n")
    cat(str(Wilcox_df_4))
    cat("\n")
    
    
    Wilcox_df_4$Related_figure<-NA
    
    
    Wilcox_df_4$Related_figure[which(Wilcox_df_4$Annotation == "Fig1_Annot_Category")]<-"Figure_1"
    Wilcox_df_4$Related_figure[which(Wilcox_df_4$Annotation == "MPRA_CLASS_2")]<-"Figure_S2"
    Wilcox_df_4$Related_figure[which(Wilcox_df_4$Annotation == "M_and_M")]<-"Figure_S4"
    Wilcox_df_4$Related_figure[which(Wilcox_df_4$Annotation == "Multi_Lineage")]<-"Figure_5"
    
    Wilcox_df_4$Related_figure<-factor(Wilcox_df_4$Related_figure,
                                       levels=c("Figure_1","Figure_S2","Figure_S4","Figure_5"),
                                       ordered=T)
    
    cat("Wilcox_df_4_1\n")
    cat(str(Wilcox_df_4))
    cat("\n")
    
    
   
    
    setwd(out)
    
    write.table(Wilcox_df_4[,-which(colnames(Wilcox_df_4) == "Annotation")],file="Table_S2.tsv", sep="\t", row.names=F,quote=F)
    
    saveRDS(Wilcox_df_4,file="STATS_value_Z_score.rds")
    
  }# length(List_DEF) >1
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
    make_option(c("--Table_of_labels"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Master_file"), type="character", default=NULL, 
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
  
  Stats_function(opt)
  
}


###########################################################################

system.time( main() )
