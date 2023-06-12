
suppressMessages(library("plyr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("data.table", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("crayon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
library("farver", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("labeling", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
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
suppressMessages(library("colorspace", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("digest", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
# library("backports", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
# library("Formula", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
# library("Hmisc", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")


library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
opt = NULL

options(warn = 1)

data_wrangling_graph_and_stats = function(option_list)
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
  
  ### Read initial selection----
  
 
  VAR_db<-as.data.frame(readRDS(file=opt$Initial_Selection) , stringsAsFactors=F)
  
  VAR_db<-unique(VAR_db[,-c(which(colnames(VAR_db) == "maf_origin"),
                            which(colnames(VAR_db) == "Fig0_Annot_Category"))])

  
  cat("VAR_db_0\n")
  cat(str(VAR_db))
  cat("\n")
  cat(str(unique(VAR_db$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_db$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_db$Fig1_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_db$Fig2_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_db$Fig2_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_db$Fig3_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_db$Fig3_Annot_Category))))
  cat("\n")
  
  Fig1<-droplevels(VAR_db[!is.na(VAR_db$Fig1_Annot_Category),])
  
  cat("Fig1_0\n")
  cat(str(Fig1))
  cat("\n")
  cat(str(unique(Fig1$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Fig1$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Fig1$Fig1_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Fig1$Fig2_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Fig1$Fig2_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Fig1$Fig3_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Fig1$Fig3_Annot_Category))))
  cat("\n")
  
  Fig2<-droplevels(VAR_db[!is.na(VAR_db$Fig2_Annot_Category),])
  
  cat("Fig2_0\n")
  cat(str(Fig2))
  cat("\n")
  cat(str(unique(Fig2$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Fig2$Fig2_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Fig2$Fig2_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Fig2$Fig2_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Fig2$Fig2_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Fig2$Fig3_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Fig2$Fig3_Annot_Category))))
  cat("\n")
  
  Fig3<-droplevels(VAR_db[!is.na(VAR_db$Fig3_Annot_Category),])
  
  cat("Fig3_0\n")
  cat(str(Fig3))
  cat("\n")
  cat(str(unique(Fig3$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Fig3$Fig3_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Fig3$Fig3_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Fig3$Fig2_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Fig3$Fig2_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Fig3$Fig3_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Fig3$Fig3_Annot_Category))))
  cat("\n")
  cat("genIE_vars\n")
  cat(sprintf(as.character(Fig3$VAR)))
  cat("\n")
  
  
  #### Create a new graph for a global representation -----
  
  Fig1_subset<-Fig1[,c(which(colnames(Fig1) == "VAR"),which(colnames(Fig1) == "rs"),which(colnames(Fig1) == "VEP_DEF_LABELS_wCSQ"),
                       which(colnames(Fig1) == "Fig1_Annot_Category"))]
  
  colnames(Fig1_subset)[which(colnames(Fig1_subset) == "Fig1_Annot_Category")]<-"variable"
  
  Fig1_subset$variable<-as.character(Fig1_subset$variable)
  
  Fig1_subset$variable[which(Fig1_subset$variable == "RV_NC_highPP_highEffectSize")]<-"index_variants"
  
  cat("Fig1_subset_0\n")
  cat(str(Fig1_subset))
  cat("\n")
  cat(str(unique(Fig1_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Fig1_subset$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Fig1_subset$variable)))))
  cat("\n")
  
  Fig2_subset<-Fig2[,c(which(colnames(Fig2) == "VAR"),which(colnames(Fig2) == "rs"),which(colnames(Fig2) == "VEP_DEF_LABELS_wCSQ"),
                       which(colnames(Fig2) == "Fig2_Annot_Category"))]
  
  colnames(Fig2_subset)[which(colnames(Fig2_subset) == "Fig2_Annot_Category")]<-"variable"
  
  Fig2_subset$variable<-as.character(Fig2_subset$variable)
  

  cat("Fig2_subset_0\n")
  cat(str(Fig2_subset))
  cat("\n")
  cat(str(unique(Fig2_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Fig2_subset$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Fig2_subset$variable)))))
  cat("\n")
  
  Fig2_subset$variable<-"MPRA_screened"
  
  cat("Fig2_subset_1\n")
  cat(str(Fig2_subset))
  cat("\n")
  cat(str(unique(Fig2_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Fig2_subset$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Fig2_subset$variable)))))
  cat("\n")
  
  
  Fig3_subset<-Fig3[,c(which(colnames(Fig3) == "VAR"),which(colnames(Fig3) == "rs"),which(colnames(Fig3) == "VEP_DEF_LABELS_wCSQ"),
                       which(colnames(Fig3) == "Fig3_Annot_Category"))]
  
  colnames(Fig3_subset)[which(colnames(Fig3_subset) == "Fig3_Annot_Category")]<-"variable"
  
  Fig3_subset$variable<-as.character(Fig3_subset$variable)
  
  
  cat("Fig3_subset_0\n")
  cat(str(Fig3_subset))
  cat("\n")
  cat(str(unique(Fig3_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Fig3_subset$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Fig3_subset$variable)))))
  cat("\n")
  
  Fig3_subset$variable<-"genIE_screened"
  
  cat("Fig3_subset_1\n")
  cat(str(Fig3_subset))
  cat("\n")
  cat(str(unique(Fig3_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Fig3_subset$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Fig3_subset$variable)))))
  cat("\n")
  
  rescue_SH2B3<-Fig1_subset[which(Fig1_subset$VAR == "chr12_111844956_C_T"),]
  rescue_SH2B3$variable<-"genIE_screened"
  
  cat("rescue_SH2B3_0\n")
  cat(str(rescue_SH2B3))
  cat("\n")
  cat(str(unique(rescue_SH2B3$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(rescue_SH2B3$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(rescue_SH2B3$variable)))))
  cat("\n")
  
  Fig3_subset<-rbind(Fig3_subset,rescue_SH2B3)
  
  cat("Fig3_subset_2\n")
  cat(str(Fig3_subset))
  cat("\n")
  cat(str(unique(Fig3_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Fig3_subset$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Fig3_subset$variable)))))
  cat("\n")
  
  #### LONG df ----
  
  LONG_df<-rbind(Fig1_subset,Fig2_subset,Fig3_subset)
  
  cat("LONG_df_0\n")
  cat(str(LONG_df))
  cat("\n")
  cat(str(unique(LONG_df$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(LONG_df$variable))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(LONG_df$variable)))))
  cat("\n")
  
  LONG_df$variable<-factor(LONG_df$variable,
                           levels=c("Common_variant","RV_C","RV_NC_lowPP","RV_NC_highPP_lowEffectSize","index_variants",
                                       "MPRA_screened","genIE_screened"),
                           ordered=T)
  
  cat("LONG_df_1\n")
  cat(str(LONG_df))
  cat("\n")
  cat(str(unique(LONG_df$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(LONG_df$variable)))))
  cat("\n")
  cat(sprintf(as.character(summary(LONG_df$variable))))
  cat("\n")
  
 
  
 
 
  
 # quit(status=1)

  #### path2 ----
  
  path2<-paste(out,'VEP_MSC','/', sep='')
  
  cat("path2\n")
  cat(sprintf(as.character(path2)))
  cat("\n")
  
  if (file.exists(path2)){
    
    
    
    
  } else {
    dir.create(file.path(path2))
    
  }
  

  
    
  #### Categories colors ----
  
  
  DEF_colors<-readRDS(file=opt$Categories_colors)
  
  cat("DEF_colors_0\n")
  cat(str(DEF_colors))
  cat("\n")
  
  # quit(status = 1)
  

  
  #### CSQ_colors colors ----
  
  
  df_CSQ_colors<-readRDS(file=opt$CSQ_colors)
  
  df_CSQ_colors$color[df_CSQ_colors$VEP_DEF_LABELS == "TFBS"]<-"red"
  
  
  cat("df_CSQ_colors_0\n")
  cat(str(df_CSQ_colors))
  cat("\n")
  
  
 
  
  
  #### graph_1 ----
  
  LONG_df.dt<-data.table(LONG_df,
                                    key=c("VEP_DEF_LABELS_wCSQ","variable"))
  
  Freq_table_sel<-as.data.frame(LONG_df.dt[, .N,
                                                      by=key(LONG_df.dt)]
                                ,stringsAsFactors=F)
  
  
  
  
  colnames(Freq_table_sel)[2]<-"variable"
  colnames(Freq_table_sel)[which(colnames(Freq_table_sel) == "N")]<-"Associations"
  
  
  cat("Freq_table_sel_0\n")
  cat(str(Freq_table_sel))
  cat("\n")
  
  
  Freq_table_sel.dt<-data.table(Freq_table_sel, key="variable")
  
  
  
  
  
  Freq_table_TOTAL<-as.data.frame(Freq_table_sel.dt[, .(Total_associations_per_category=sum(Associations)), 
                                                    by=key(Freq_table_sel.dt)], stringsAsFactors=F)
  
  
  
  cat("Freq_table_TOTAL_\n")
  cat(str(Freq_table_TOTAL))
  cat("\n")
  
  
  Plot_table<-as.data.frame(merge(Freq_table_sel,
                                  Freq_table_TOTAL,
                                  by="variable",
                                  all.x=T), stringsAsFactors=F)
  
  cat("Plot_table_0_line365\n")
  cat(str(Plot_table))
  cat("\n")
  
  levels_CSQ_sel<-levels(LONG_df$VEP_DEF_LABELS_wCSQ)
  
  cat("levels_CSQ_sel\n")
  cat(str(levels_CSQ_sel))
  cat("\n")
  
  levels_category<-levels(LONG_df$variable)
  
  cat("levels_category\n")
  cat(str(levels_category))
  cat("\n")
  
  indx.int_2<-which(colnames(Plot_table) == "variable")
  
  cat("indx.int_2:\n")
  cat(sprintf(as.character(indx.int_2)))
  cat("\n")
  
  indx.CSQ<-which(colnames(Plot_table) == "VEP_DEF_LABELS_wCSQ")
  
  
  Plot_table[,indx.int_2]<-factor(Plot_table[,indx.int_2],
                                  levels=rev(levels_category),
                                  ordered=T)
  
  Plot_table[,indx.CSQ]<-factor(Plot_table[,indx.CSQ],
                                levels=levels_CSQ_sel,
                                ordered=T)
  
  Plot_table$Perc<-round(100*(Plot_table$Associations/Plot_table$Total_associations_per_category),1)
  
  
  cat("Plot_table_1\n")
  cat(str(Plot_table))
  cat("\n")
  
  Table_K<-Plot_table
  
  colnames(Table_K)[indx.int_2]<-"LEVELS"
  colnames(Table_K)[indx.CSQ]<-"VEP_DEF_LABELS"
  colnames(Table_K)[which(colnames(Table_K) == "Associations")]<-"n"
  Table_K$Category<-"variable"
  
  cat("Table_K_1\n")
  cat(str(Table_K))
  cat("\n")
  

  
  
  df_CSQ_colors_sel<-df_CSQ_colors[which(df_CSQ_colors$VEP_DEF_LABELS%in%levels_CSQ_sel),]
  
  cat("df_CSQ_colors_sel_1\n")
  cat(str(df_CSQ_colors_sel))
  cat("\n")
  
  breaks.Rank<-seq(0,100,by=10)
  labels.Rank<-as.character(breaks.Rank)
  #   
  
  
  graph_1<-Plot_table %>%
    mutate(myaxis = paste0(Plot_table[,indx.int_2], "\n", "n=", Total_associations_per_category)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Plot_table[,indx.int_2]))) %>%
    ggplot(aes(x=Perc, y=myaxis, fill=Plot_table[,indx.CSQ])) +
    geom_bar(stat="identity",colour='black')+
    theme_bw()+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
          axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
          legend.title=element_text(size=16,color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_y_discrete(name=NULL, drop=F)+
    scale_x_continuous(name="Percentage",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
    scale_fill_manual(values=df_CSQ_colors_sel$color,drop=F)+
    theme(legend.position="right")+
    guides(fill=guide_legend(title=paste("VEP MSC", collapse="\n")))+
    ggeasy::easy_center_title()
  
  # graph_1<-graph_1 +coord_flip()
  
  
  setwd(path2)
  
  svgname<-paste("VEP_DEF_LABELS_wCSQ","_","variable",".svg",sep='')
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= graph_1,
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
  option_list <- list(
    make_option(c("--Initial_Selection"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CSQ_colors"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Categories_colors"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VEP_CSQ"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ALL_dB"), type="character", default=NULL, 
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
  
  data_wrangling_graph_and_stats(opt)
  
}


###########################################################################

system.time( main() )
