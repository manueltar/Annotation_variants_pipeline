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
suppressMessages(library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))


opt = NULL

Contigency_annotation_ALL = function(option_list)
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
  
  cat("OUT\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform tracking_variants ----
  
  tracking_variants = unlist(strsplit(opt$tracking_variants, split=","))
  
  cat("tracking_variants_\n")
  cat(sprintf(as.character(tracking_variants)))
  cat("\n")
  
  #### Table_S6_Manual_curation ----
  
  Table_S6<-readRDS(opt$Table_S6_Manual_curation)
  
  Table_S6<-droplevels(Table_S6[-which(Table_S6$Mechanistic_Class == 'No_RNA_Seq_HET_carriers'),])
  
  cat("Table_S6_0\n")
  cat(str(Table_S6))
  cat("\n")
  cat(str(unique(Table_S6$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6$MPRA_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6$genIE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6$genIE_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6$Mechanistic_Class))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6$Manual_curation)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6$Manual_curation))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6$Multi_Lineage))))
  cat("\n")
  
  levels_MPRA_CLASS<-levels(Table_S6$MPRA_CLASS)
  
  cat("levels_MPRA_CLASS\n")
  cat(str(levels_MPRA_CLASS))
  cat("\n")
  
  levels_genIE_CLASS<-levels(Table_S6$genIE_CLASS)
  
  cat("levels_genIE_CLASS\n")
  cat(str(levels_genIE_CLASS))
  cat("\n")
  
  levels_Mechanistic_Class<-levels(Table_S6$Mechanistic_Class)
  
  cat("levels_Mechanistic_Class\n")
  cat(str(levels_Mechanistic_Class))
  cat("\n")
  
  levels_Manual_curation<-levels(Table_S6$Manual_curation)
  
  cat("levels_Manual_curation\n")
  cat(str(levels_Manual_curation))
  cat("\n")
  
  levels_Multi_Lineage<-levels(Table_S6$Multi_Lineage)
  
  cat("levels_Multi_Lineage\n")
  cat(str(levels_Multi_Lineage))
  cat("\n")
  
 
  #### Read range_analysis ----
  
  
  range_analysis<-as.data.frame(fread(file=opt$range_analysis,sep="\t",header=T, stringsAsFactors = F))
  
  cat("range_analysis\n")
  str(range_analysis)
  cat("\n")
  
  #### Contingency_table ----
  
  Contingency_table<-Table_S6
  
  Contingency_table$regulatory_status<-NA
  
  
  Contingency_table$regulatory_status[which(Contingency_table$Mechanistic_Class != 'NRD')]<-'Regulated'
  Contingency_table$regulatory_status[which(Contingency_table$Mechanistic_Class == 'NRD')]<-'NRD'
  
  Contingency_table$regulatory_status<-factor(Contingency_table$regulatory_status,
                                              levels=c('NRD','Regulated'),
                                              ordered = T)
  
  Contingency_table$interaction_1<-interaction(Contingency_table$MPRA_CLASS,Contingency_table$regulatory_status, sep="|",lex.order = T)
  
  cat("Contingency_table_0\n")
  cat(str(Contingency_table))
  cat("\n")
  cat(str(unique(Contingency_table$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Contingency_table$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Contingency_table$MPRA_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Contingency_table$regulatory_status)))))
  cat("\n")
  cat(sprintf(as.character(summary(Contingency_table$regulatory_status))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Contingency_table$interaction_1)))))
  cat("\n")
  cat(sprintf(as.character(summary(Contingency_table$interaction_1))))
  cat("\n")
  
  Contingency_table_subset<-Contingency_table[,c(which(colnames(Contingency_table) == 'VAR'),
                                                 which(colnames(Contingency_table) == 'interaction_1'))]
  
  cat("Contingency_table_subset_0\n")
  cat(str(Contingency_table_subset))
  cat("\n")
  cat(str(unique(Contingency_table_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Contingency_table_subset$interaction_1)))))
  cat("\n")
  cat(sprintf(as.character(summary(Contingency_table_subset$interaction_1))))
  cat("\n")
  
  #### Merge with range analysis ----
  
  Contingency_table_subset<-merge(Contingency_table_subset,
                                  range_analysis,
                                  by='VAR')
  
  cat("Contingency_table_subset_1\n")
  cat(str(Contingency_table_subset))
  cat("\n")
  cat(str(unique(Contingency_table_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Contingency_table_subset$interaction_1)))))
  cat("\n")
  cat(sprintf(as.character(summary(Contingency_table_subset$interaction_1))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Contingency_table_subset$RNASeq_source))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Contingency_table_subset$RNASeq_source)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Contingency_table_subset$nHET)))))
  cat("\n")
  cat(sprintf(as.character(summary(Contingency_table_subset$nHET))))
  cat("\n")
  
  Contingency_table_subset_WB<-droplevels(unique(Contingency_table_subset[which(Contingency_table_subset$RNASeq_source == 'Whole blood'),]))
  
  cat("Contingency_table_subset_WB_0\n")
  cat(str(Contingency_table_subset_WB))
  cat("\n")
  cat(str(unique(Contingency_table_subset_WB$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Contingency_table_subset_WB$interaction_1)))))
  cat("\n")
  cat(sprintf(as.character(summary(Contingency_table_subset_WB$interaction_1))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Contingency_table_subset_WB$RNASeq_source))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Contingency_table_subset_WB$RNASeq_source)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Contingency_table_subset_WB$nHET)))))
  cat("\n")
  cat(sprintf(as.character(summary(Contingency_table_subset_WB$nHET))))
  cat("\n")
  
  
  
  range_analysis_table<-unique(Contingency_table_subset_WB[,c(which(colnames(Contingency_table_subset_WB) == 'VAR'),
                                                       which(colnames(Contingency_table_subset_WB) == 'nHET'))])
  
  cat("range_analysis_table_0\n")
  cat(str(range_analysis_table))
  cat("\n")
  cat(str(unique(range_analysis_table$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(range_analysis_table$nHET)))))
  cat("\n")
  cat(sprintf(as.character(summary(range_analysis_table$nHET))))
  cat("\n")
  
  check<-range_analysis_table[which(range_analysis_table$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  
  # #### SAVE ----
  
  setwd(out)
  
  write.table(range_analysis_table, file="range_analysis_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
  
  
  Contigency_annotation_ALL<-unique(Contingency_table_subset_WB[,c(which(colnames(Contingency_table_subset_WB) == 'VAR'),
                                                              which(colnames(Contingency_table_subset_WB) == 'interaction_1'))])
  
  colnames(Contigency_annotation_ALL)[which(colnames(Contigency_annotation_ALL) == 'interaction_1')]<-'Contigency_annotation_ALL'
  
  cat("Contigency_annotation_ALL_0\n")
  cat(str(Contigency_annotation_ALL))
  cat("\n")
  cat(str(unique(Contigency_annotation_ALL$VAR)))
  cat("\n")
  
  write.table(Contigency_annotation_ALL, file="Contigency_annotation_ALL.tsv", sep="\t", quote = F, row.names = F)

 
}

statistics_range_analysis = function(option_list)
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
  
  cat("OUT\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Read the contigency table and the n per variant of the previous analysis----
  
  setwd(out)
  
  filename<-"range_analysis_GLOBAL.tsv"
  
  range_analysis_table<-as.data.frame(fread(file=filename, sep="\t", header=T), stringsAsFactors=F)
  
  cat("range_analysis_table_0\n")
  cat(str(range_analysis_table))
  cat("\n")
  cat(str(unique(range_analysis_table$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(range_analysis_table$nHET)))))
  cat("\n")
  cat(sprintf(as.character(summary(range_analysis_table$nHET))))
  cat("\n")
  
  filename<-"Contigency_annotation_ALL.tsv"
  
  Contigency_annotation_ALL<-as.data.frame(fread(file=filename, sep="\t", header=T), stringsAsFactors=F)
  
  cat("Contigency_annotation_ALL_0\n")
  cat(str(Contigency_annotation_ALL))
  cat("\n")
  cat(str(unique(Contigency_annotation_ALL$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Contigency_annotation_ALL$Contigency_annotation_ALL))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Contigency_annotation_ALL$Contigency_annotation_ALL)))))
  cat("\n")
  
  colnames(Contigency_annotation_ALL)[which(colnames(Contigency_annotation_ALL) == 'Contigency_annotation_ALL')]<-'CLASS'
  
  Contigency_annotation_ALL$CLASS<-factor(Contigency_annotation_ALL$CLASS,
                                          levels=c('no MPRA hit|Regulated','MPRA hit|Regulated','no MPRA hit|NRD','MPRA hit|NRD'),
                                          ordered=T)
  
  cat("Contigency_annotation_ALL_1\n")
  cat(str(Contigency_annotation_ALL))
  cat("\n")
  cat(str(unique(Contigency_annotation_ALL$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Contigency_annotation_ALL$CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Contigency_annotation_ALL$CLASS))))
  cat("\n")
  
  
  #### merge variants ----
  
  Contigency_annotation_ALL<-merge(Contigency_annotation_ALL,
                                   range_analysis_table,
                                   by="VAR")
  
  cat("Contigency_annotation_ALL_2\n")
  cat(str(Contigency_annotation_ALL))
  cat("\n")
  cat(str(unique(Contigency_annotation_ALL$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Contigency_annotation_ALL$CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Contigency_annotation_ALL$CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Contigency_annotation_ALL$nHET)))))
  cat("\n")
  cat(sprintf(as.character(summary(Contigency_annotation_ALL$nHET))))
  cat("\n")
  
  #### Summary table ----
  
  
  Contigency_annotation_ALL.dt<-data.table(Contigency_annotation_ALL, key=c("CLASS"))
  
  Summary_table<-as.data.frame(Contigency_annotation_ALL.dt[,.(n=.N,
                                                                 Min=round(as.numeric(summary(nHET)[1]),3),
                                                                 Q1=round(as.numeric(summary(nHET)[2]),3),
                                                                 M=round(as.numeric(summary(nHET)[3]),3),
                                                                 Q3=round(as.numeric(summary(nHET)[5]),3),
                                                                 Max=round(as.numeric(summary(nHET)[6]),3)),, by=key(Contigency_annotation_ALL.dt)], stringsAsFactors=F)
  
  
  cat("Summary_table_0\n")
  cat(str(Summary_table))
  cat("\n")
  
  #### Pairwise comparisons wilcoxon ----
  
  CONDITION_DEBUG <- 1
  
  PW_comparison<-pairwise.wilcox.test(Contigency_annotation_ALL$nHET, Contigency_annotation_ALL$CLASS,
                                      p.adjust.method = "BH")
  
  
  if(CONDITION_DEBUG == 1)
  {
    cat("PW_comparison\n")
    cat(str(PW_comparison))
    cat("\n")
  }
  
  PW_comparison_pvalue_df<-as.data.frame(PW_comparison$p.value, stringsAsFactors=F)
  
  if(CONDITION_DEBUG == 1)
  {
    cat("PW_comparison_pvalue_df\n")
    cat(str(PW_comparison_pvalue_df))
    cat("\n")
  }
  
  list_cols<-list()
  
  for(PW_iteration in 1:dim(PW_comparison_pvalue_df)[2])
  {
    colnames_sel<-colnames(PW_comparison_pvalue_df)[PW_iteration]
    
    if(CONDITION_DEBUG == 1)
    {
      cat("----------------->colnames_sel\n")
      cat(sprintf(as.character(colnames_sel)))
      cat("\n")
    }
    
    list_rows<-list()
    
    
    for(PW_iteration_k in 1:dim(PW_comparison_pvalue_df)[1])
    {
      rownames_sel<-row.names(PW_comparison_pvalue_df)[PW_iteration_k]
      
      if(CONDITION_DEBUG == 1)
      {
        cat("--->rownames_sel\n")
        cat(sprintf(as.character(rownames_sel)))
        cat("\n")
      }
      
      PW_Wilcox_pvalue<-PW_comparison_pvalue_df[PW_iteration_k,PW_iteration]
      
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
      
      FLAG_NA<-sum(is.na(log_pval_PW_Wilcox))
      
      
      if(CONDITION_DEBUG == 1)
      {
        cat("FLAG_NA\n")
        cat(sprintf(as.character(FLAG_NA)))
        cat("\n")
      }
      
      if(FLAG_NA ==0)
      {
        vector_final_comparisons<-paste(sort(c(colnames_sel,rownames_sel)), collapse=";")
        
        
        a.dt<-as.data.frame(cbind(vector_final_comparisons,PW_Wilcox_pvalue,log_pval_PW_Wilcox), stringsAsFactors=F)
        
        colnames(a.dt)<-c("comparison",'pval','MINUS_logpval')
        
        if(CONDITION_DEBUG == 1)
        {
          cat("a.dt\n")
          cat(str(a.dt))
          cat("\n")
        }
        
        list_rows[[PW_iteration_k]]<-a.dt
        
      }#FLAG_NA ==0
      
      
    }#PW_iteration_k
    
    if(length(list_rows) >0)
    {
      
      df_col = as.data.frame(data.table::rbindlist(list_rows, fill=T), stringsAsFactors=F)
      
      if(CONDITION_DEBUG == 1)
      {
        cat("df_col\n")
        cat(str(df_col))
        cat("\n")
      }
      
      list_cols[[PW_iteration]]<-df_col
    }#length(list_rows) >0
    
  }#PW_iteration
  
  if(length(list_cols) >0)
  {
    PW_comparison = as.data.frame(data.table::rbindlist(list_cols, fill=T), stringsAsFactors=F)
    
    PW_comparison[,which(colnames(PW_comparison) == 'pval')]<-as.numeric(PW_comparison[,which(colnames(PW_comparison) == 'pval')])
    PW_comparison[,which(colnames(PW_comparison) == 'MINUS_logpval')]<-as.numeric(PW_comparison[,which(colnames(PW_comparison) == 'MINUS_logpval')])
    
    if(CONDITION_DEBUG == 1)
    {
      cat("PW_comparison\n")
      cat(str(PW_comparison))
      cat("\n")
    }
    
    
    PW_comparison_NO_NA<-PW_comparison[!is.na(PW_comparison[,which(colnames(PW_comparison) == 'MINUS_logpval')]),]
    
    PW_comparison_NO_NA$c1<-gsub(";.+$","",PW_comparison_NO_NA$comparison)
    PW_comparison_NO_NA$c2<-gsub("^[^;]+;","",PW_comparison_NO_NA$comparison)
    
    
    if(CONDITION_DEBUG == 1)
    {
      cat("PW_comparison_NO_NA_0\n")
      cat(str(PW_comparison_NO_NA))
      cat("\n")
    }
    
    Summary_table_c1<-Summary_table
    colnames(Summary_table_c1)[which(colnames(Summary_table_c1) == 'CLASS')]<-'c1'
    colnames(Summary_table_c1)[-1]<-paste('c1',colnames(Summary_table_c1)[-1], sep='_')
    
    Summary_table_c2<-Summary_table
    colnames(Summary_table_c2)[which(colnames(Summary_table_c2) == 'CLASS')]<-'c2'
    colnames(Summary_table_c2)[-1]<-paste('c2',colnames(Summary_table_c2)[-1], sep='_')
    
    
    PW_comparison_NO_NA<-merge(PW_comparison_NO_NA,
                               Summary_table_c2,
                               by='c2')
    
    PW_comparison_NO_NA<-merge(PW_comparison_NO_NA,
                               Summary_table_c1,
                               by='c1')
    
    ##### Save ----
    
    setwd(out)
    
    # write.table(Summary_table,file='test.tsv', sep="\t", quote = F, row.names = F)
    write.table(PW_comparison_NO_NA,file='BK_MPRA_vs_RNA_Range_analysis.tsv', sep="\t", quote = F, row.names = F)
    
    saveRDS(file='Contigency_annotation_ALL_with_range.rds',Contigency_annotation_ALL)
    
    
  }#length(list_cols) >0
  
  
}


violin_plots = function(option_list)
{
  library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
  
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
  
  cat("OUT\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Read the contigency table and the n per variant of the previous analysis----
  
  setwd(out)
  
  filename<-"Contigency_annotation_ALL_with_range.rds"
  
  Contigency_annotation_ALL<-readRDS(file=filename)
  
  cat("Contigency_annotation_ALL_0\n")
  cat(str(Contigency_annotation_ALL))
  cat("\n")
  cat(str(unique(Contigency_annotation_ALL$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Contigency_annotation_ALL$CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Contigency_annotation_ALL$CLASS))))
  cat("\n")
  
  
  #### sina plot graph ----
  
  
  summary_nHET<-summary(Contigency_annotation_ALL$nHET[!is.na(Contigency_annotation_ALL$nHET)])
  
  max_nHET<-max(summary_nHET)
  min_nHET<-min(summary_nHET)
  
  breaks.Rank<-seq(from=0,max_nHET+10, by=10)
  labels.Rank<-as.character(breaks.Rank)
  
  cat("labels.Rank\n")
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
  
  
  graph<-ggplot(data=Contigency_annotation_ALL,
                aes(x=CLASS, y=nHET, color=CLASS)) +
    geom_sina()+
    geom_boxplot(width = 0.2)+
    theme_classic()+
    theme(plot.title = element_text(size=11)) +
    theme(plot.title=element_text(size=11))+
    scale_x_discrete(name=NULL, drop=F)+
    theme(axis.text.x=element_text(angle=45,size=14,colour="black",vjust=1,hjust=1), 
          axis.title.y=element_text(size=16),
          legend.title=element_text(size=16),
          legend.text=element_text(size=12,colour="black"),
          axis.text.y=element_text(colour="black",size=12))+
    scale_y_continuous(name="Number of heterozygous carriers",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    labs(color=paste("MPRA|RNA regulation","classification", sep="\n"))+
    scale_color_manual(values = c('#32A852','#1877C9','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6'),
                       drop=F)+
    ggeasy::easy_center_title()
  
  
  cat("graph_DONE\n")
  
  setwd(out)
  
  svgname<-paste("sina_plot_","MPRA_vs_RNA_","number_of_HETs",".svg",sep='')
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= graph,
           device="svg",
           height=10, width=12)
  }
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
  traceback()
  options(show.error.locations = TRUE)
  
  option_list <- list(
    make_option(c("--range_analysis"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--tracking_variants"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Table_S6_Manual_curation"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
       make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
        make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out_2"), type="character", default=NULL, 
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
  
 
 # Contigency_annotation_ALL(opt)
 # statistics_range_analysis(opt)
 violin_plots(opt)
 
}




###########################################################################

system.time( main() )
