
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
library("cowplot",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")



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
 
  #### RMV DGKQ ----
  
  
  RMV_not_index = unlist(strsplit(opt$RMV_not_index, split=","))
  
  cat("RMV_not_index_\n")
  cat(sprintf(as.character(RMV_not_index)))
  cat("\n")
  
  #### RMV DGKQ ----
  
  
  RMV_labels = unlist(strsplit(opt$RMV_labels, split=","))
  
  cat("RMV_labels_\n")
  cat(sprintf(as.character(RMV_labels)))
  cat("\n")
  
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
  
 
  ### Read Table_S4----
  
  
  Table_S4<-as.data.frame(readRDS(file=opt$Table_S4) , stringsAsFactors=F)
  
  cat("Table_S4_0\n")
  cat(str(Table_S4))
  cat("\n")
  
  Table_S4<-Table_S4[-which(Table_S4$VAR%in%RMV_not_index),]
  
  cat("Table_S4_RMV_not_index\n")
  cat(str(Table_S4))
  cat("\n")
  
  Table_S4$genIE_CLASS_2<-revalue(Table_S4$genIE_CLASS,
                                         c("genIE_INACTIVE" = "no_genIE_hit",
                                           "genIE_ACTIVE" = "genIE_hit"))
  
  Table_S4$MPRA_CLASS_2<-revalue(Table_S4$MPRA_CLASS,
                                        c("NO_enhancer_activity" = "no_MPRA_hit",
                                          "AT_LEAST_1_TILE_with_enhancer_activity" = "no_MPRA_hit",
                                          "AT_LEAST_1_TILE_with_E_Plus_ASE_activity" = "MPRA_hit"))
  
  Table_S4$Mechanistic_Class<-revalue(Table_S4$Mechanistic_Class,
                                             c('Transcriptional_Regulation' = 'Transcriptional R.',
                                               'Alternative_Transcript_Usage' = 'ATU',
                                               'Transcriptional_Regulation_and_ATU'='Transcriptional R. + ATU',
                                               'No_Regulation_Detected' = 'NRD'))
  
  
  indx.int<-c(which(colnames(Table_S4) == "VAR"),which(colnames(Table_S4) == "Mechanistic_Class"),which(colnames(Table_S4) == "Manual_curation"),
              which(colnames(Table_S4) == "MPRA_CLASS_2"),which(colnames(Table_S4) == "genIE_CLASS_2"),which(colnames(Table_S4) == "Multi_Lineage"))
  
  Table_S4_subset<-unique(Table_S4[,indx.int])
  
  cat("Table_S4_subset_0\n")
  cat(str(Table_S4_subset))
  cat("\n")
  cat(str(unique(Table_S4_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$Mechanistic_Class))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$Manual_curation)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$Manual_curation))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$MPRA_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$genIE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$genIE_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$Multi_Lineage))))
  cat("\n")
  
  
  
  
  Table_S4_subset$MPRA_CLASS_2[which(Table_S4_subset$MPRA_CLASS_2 == RMV_labels[1])]<-NA
  
  Table_S4_subset$genIE_CLASS[which(Table_S4_subset$genIE_CLASS == RMV_labels[2])]<-NA
  
  Table_S4_subset$Mechanistic_Class[which(Table_S4_subset$Mechanistic_Class == RMV_labels[3])]<-NA
  
  Table_S4_subset$Manual_curation[which(Table_S4_subset$Manual_curation == RMV_labels[3])]<-NA
  
  Table_S4_subset$M_and_M<-interaction(Table_S4_subset$Mechanistic_Class,Table_S4_subset$Manual_curation, sep="|",lex.order = T)#, na.omit=T)
  
  
  
  Table_S4_subset<-droplevels(Table_S4_subset)
  
  
  cat("Table_S4_subset_1\n")
  cat(str(Table_S4_subset))
  cat("\n")
  cat(str(unique(Table_S4_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$Mechanistic_Class))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$Manual_curation)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$Manual_curation))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$MPRA_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$genIE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$genIE_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$Multi_Lineage))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$M_and_M)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$M_and_M))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$MPRA_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$MPRA_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_subset$genIE_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_subset$genIE_CLASS_2))))
  cat("\n")
 
  
  #### Merge all labels ----
  
  Table_of_labels<-merge(VAR_Prioritization_dB_subset,
                         Table_S4_subset,
                         by="VAR",
                         all=T)
  
  # Table_of_labels$chr<-gsub("_.+$","",Table_of_labels$VAR)
  # 
  # Table_of_labels$pos<-gsub("^[^_]+_","",Table_of_labels$VAR)
  # Table_of_labels$pos<-gsub("_.+$","",Table_of_labels$pos)
  # 
 
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
  
  
  ##### define sets to overlap -----
  
  index_set<-droplevels(Table_of_labels[which(Table_of_labels$Fig1_Annot_Category == 'RV_NC_highPP_highEffectSize'),])
  
  cat("index_set_0\n")
  cat(str(index_set))
  cat("\n")
  cat(str(unique(index_set$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(index_set$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(index_set$Fig1_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(index_set$MPRA_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(index_set$MPRA_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(index_set$genIE_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(index_set$genIE_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(index_set$M_and_M)))))
  cat("\n")
  cat(sprintf(as.character(summary(index_set$M_and_M))))
  cat("\n")
  cat(sprintf(as.character(names(summary(index_set$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(index_set$Multi_Lineage))))
  cat("\n")
  
 
  MPRA_set<-droplevels(Table_of_labels[!is.na(Table_of_labels$MPRA_CLASS_2),])
  
  cat("MPRA_set_0\n")
  cat(str(MPRA_set))
  cat("\n")
  cat(str(unique(MPRA_set$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_set$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_set$Fig1_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_set$MPRA_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_set$MPRA_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_set$genIE_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_set$genIE_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_set$M_and_M)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_set$M_and_M))))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_set$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_set$Multi_Lineage))))
  cat("\n")
  
  genIE_set<-droplevels(Table_of_labels[which(Table_of_labels$genIE_CLASS_2 != 'NOT_SCREENED_genIE'),])
  
  cat("genIE_set_0\n")
  cat(str(genIE_set))
  cat("\n")
  cat(str(unique(genIE_set$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(genIE_set$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(genIE_set$Fig1_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(genIE_set$genIE_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(genIE_set$genIE_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(genIE_set$genIE_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(genIE_set$genIE_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(genIE_set$M_and_M)))))
  cat("\n")
  cat(sprintf(as.character(summary(genIE_set$M_and_M))))
  cat("\n")
  cat(sprintf(as.character(names(summary(genIE_set$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(genIE_set$Multi_Lineage))))
  cat("\n")
 
  rna_seq_set<-droplevels(Table_of_labels[!is.na(Table_of_labels$Mechanistic_Class),])
  
  cat("rna_seq_set_0\n")
  cat(str(rna_seq_set))
  cat("\n")
  cat(str(unique(rna_seq_set$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(rna_seq_set$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(rna_seq_set$Fig1_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(rna_seq_set$MPRA_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(rna_seq_set$MPRA_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(rna_seq_set$genIE_CLASS_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(rna_seq_set$genIE_CLASS_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(rna_seq_set$M_and_M)))))
  cat("\n")
  cat(sprintf(as.character(summary(rna_seq_set$M_and_M))))
  cat("\n")
  cat(sprintf(as.character(names(summary(rna_seq_set$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(rna_seq_set$Multi_Lineage))))
  cat("\n")
  
  
  ###### overlaps 1 -----
  
  
  Ov_index_MPRA<-index_set[which(index_set$VAR%in%MPRA_set$VAR),]
  
  cat("Ov_index_MPRA_0\n")
  cat(str(Ov_index_MPRA))
  cat("\n")
  cat(str(unique(Ov_index_MPRA$VAR)))
  cat("\n")
  
  Ov_index_genIE<-index_set[which(index_set$VAR%in%genIE_set$VAR),]
  
  cat("Ov_index_genIE_0\n")
  cat(str(Ov_index_genIE))
  cat("\n")
  cat(str(unique(Ov_index_genIE$VAR)))
  cat("\n")
  
  Ov_index_rna_seq<-index_set[which(index_set$VAR%in%rna_seq_set$VAR),]
  
  cat("Ov_index_rna_seq_0\n")
  cat(str(Ov_index_rna_seq))
  cat("\n")
  cat(str(unique(Ov_index_rna_seq$VAR)))
  cat("\n")
 
  Ov_genIE_MPRA<-genIE_set[which(genIE_set$VAR%in%MPRA_set$VAR),]
  
  cat("Ov_genIE_MPRA_0\n")
  cat(str(Ov_genIE_MPRA))
  cat("\n")
  cat(str(unique(Ov_genIE_MPRA$VAR)))
  cat("\n")
  
  Ov_genIE_rna_seq<-genIE_set[which(genIE_set$VAR%in%rna_seq_set$VAR),]
  
  cat("Ov_genIE_rna_seq_0\n")
  cat(str(Ov_genIE_rna_seq))
  cat("\n")
  cat(str(unique(Ov_genIE_rna_seq$VAR)))
  cat("\n")
  
  Ov_MPRA_rna_seq<-MPRA_set[which(MPRA_set$VAR%in%rna_seq_set$VAR),]
  
  cat("Ov_MPRA_rna_seq_0\n")
  cat(str(Ov_MPRA_rna_seq))
  cat("\n")
  cat(str(unique(Ov_MPRA_rna_seq$VAR)))
  cat("\n")
  
 
  Minus_Ov_rna_seq_index<-rna_seq_set[-which(rna_seq_set$VAR%in%index_set$VAR),]
  
  cat("Minus_Ov_rna_seq_index_0\n")
  cat(str(Minus_Ov_rna_seq_index))
  cat("\n")
  cat(str(unique(Minus_Ov_rna_seq_index$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Minus_Ov_rna_seq_index$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Minus_Ov_rna_seq_index$Fig1_Annot_Category))))
  cat("\n")
  
  ################### Set def ----
  
  Ov_rna_seq_index<-rna_seq_set[which(rna_seq_set$VAR%in%index_set$VAR),]
  
  cat("Ov_rna_seq_index_0\n")
  cat(str(Ov_rna_seq_index))
  cat("\n")
  cat(str(unique(Ov_rna_seq_index$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Ov_rna_seq_index$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Ov_rna_seq_index$Fig1_Annot_Category))))
  cat("\n")
  
  Ov_rna_seq_index.m<-melt(Ov_rna_seq_index[,c(which(colnames(Ov_rna_seq_index) == 'VAR'),
                                               which(colnames(Ov_rna_seq_index) == 'Multi_Lineage'),
                                               which(colnames(Ov_rna_seq_index) == 'MPRA_CLASS_2'),
                                               which(colnames(Ov_rna_seq_index) == 'genIE_CLASS_2'),
                                               which(colnames(Ov_rna_seq_index) == 'Mechanistic_Class'),
                                               which(colnames(Ov_rna_seq_index) == 'Manual_curation'))], id.vars=c("VAR",'Multi_Lineage'), variable.name = "CLASS", value.name="VALUE")
  
  Ov_rna_seq_index.m<-Ov_rna_seq_index.m[order(Ov_rna_seq_index.m$VAR,Ov_rna_seq_index.m$CLASS),]
  
  cat("Ov_rna_seq_index.m_0\n")
  cat(str(Ov_rna_seq_index.m))
  cat("\n")
  
  
  Ov_rna_seq_index.m.dt<-data.table(Ov_rna_seq_index.m,key=c("VAR",'Multi_Lineage'))
  
  cat("Ov_rna_seq_index.m.dt\n")
  cat(str(Ov_rna_seq_index.m.dt))
  cat("\n")
  
  
  
  Ov_rna_seq_index.m_CLASS_summarised<-as.data.frame(Ov_rna_seq_index.m.dt[,.(CLASS_string=paste(VALUE,collapse="|")), by=key(Ov_rna_seq_index.m.dt)],stringsAsFactors=F)
  
  cat("Ov_rna_seq_index.m_CLASS_summarised\n")
  cat(str(Ov_rna_seq_index.m_CLASS_summarised))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Ov_rna_seq_index.m_CLASS_summarised$CLASS_string))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Ov_rna_seq_index.m_CLASS_summarised$CLASS_string)))))
  cat("\n")
  
  Ov_rna_seq_index.m_CLASS_summarised$Classif_DEF<-"NA"
  
  Ov_rna_seq_index.m_CLASS_summarised$Classif_DEF<-factor(Ov_rna_seq_index.m_CLASS_summarised$CLASS_string,
                                           levels = c('MPRA_hit|genIE_hit|NRD|Unknown_mechanism','MPRA_hit|genIE_hit|Transcriptional R. + ATU|R_in_candidate','MPRA_hit|NOT_SCREENED_genIE|ATU|Coding_proxy',
                                                      'MPRA_hit|NOT_SCREENED_genIE|ATU|R_in_candidate','MPRA_hit|NOT_SCREENED_genIE|ATU|R_in_non_candidate','MPRA_hit|NOT_SCREENED_genIE|NRD|Coding_proxy',
                                                      'MPRA_hit|NOT_SCREENED_genIE|NRD|Other_tissue_needed','MPRA_hit|NOT_SCREENED_genIE|NRD|Regulatory_proxy','MPRA_hit|NOT_SCREENED_genIE|NRD|Unknown_mechanism',
                                                      'MPRA_hit|NOT_SCREENED_genIE|Transcriptional R. + ATU|Other_tissue_needed','MPRA_hit|NOT_SCREENED_genIE|Transcriptional R.|Coding_proxy','MPRA_hit|NOT_SCREENED_genIE|Transcriptional R.|Other_tissue_needed',
                                                      'MPRA_hit|NOT_SCREENED_genIE|Transcriptional R.|R_in_candidate',
                                                      'no_MPRA_hit|genIE_hit|NRD|Unknown_mechanism','no_MPRA_hit|genIE_hit|Transcriptional R.|R_in_candidate','no_MPRA_hit|genIE_hit|Transcriptional R.|R_in_non_candidate',
                                                      'no_MPRA_hit|no_genIE_hit|NRD|Unknown_mechanism','no_MPRA_hit|NOT_SCREENED_genIE|ATU|R_in_candidate','no_MPRA_hit|NOT_SCREENED_genIE|ATU|R_in_non_candidate',
                                                      'no_MPRA_hit|NOT_SCREENED_genIE|NRD|Coding_proxy','no_MPRA_hit|NOT_SCREENED_genIE|NRD|Other_tissue_needed','no_MPRA_hit|NOT_SCREENED_genIE|NRD|Regulatory_proxy','no_MPRA_hit|NOT_SCREENED_genIE|NRD|Unknown_mechanism',
                                                      'no_MPRA_hit|NOT_SCREENED_genIE|Transcriptional R. + ATU|Other_tissue_needed','no_MPRA_hit|NOT_SCREENED_genIE|Transcriptional R. + ATU|R_in_candidate','no_MPRA_hit|NOT_SCREENED_genIE|Transcriptional R.|Coding_proxy','no_MPRA_hit|NOT_SCREENED_genIE|Transcriptional R.|Other_tissue_needed',
                                                      'no_MPRA_hit|NOT_SCREENED_genIE|Transcriptional R.|R_in_candidate','no_MPRA_hit|NOT_SCREENED_genIE|Transcriptional R.|R_in_non_candidate'),
                                           ordered=T)
  
  cat(sprintf(as.character(names(summary(Ov_rna_seq_index.m_CLASS_summarised$Classif_DEF)))))
  cat("\n")
  cat(sprintf(as.character(summary(Ov_rna_seq_index.m_CLASS_summarised$Classif_DEF))))
  cat("\n")
  
  Ov_rna_seq_index.m_CLASS_summarised_Classif_DEF_TOTAL.dt<-data.table(Ov_rna_seq_index.m_CLASS_summarised, key=c("CLASS_string"))
  
  
  Ov_rna_seq_index.m_CLASS_summarised_Classif_DEF_TOTAL<-as.data.frame(Ov_rna_seq_index.m_CLASS_summarised_Classif_DEF_TOTAL.dt[,.N,by=key(Ov_rna_seq_index.m_CLASS_summarised_Classif_DEF_TOTAL.dt)], stringsAsFactors=F)
  
  colnames(Ov_rna_seq_index.m_CLASS_summarised_Classif_DEF_TOTAL)[which(colnames(Ov_rna_seq_index.m_CLASS_summarised_Classif_DEF_TOTAL) == "N")]<-"instances"
  
  cat("Ov_rna_seq_index.m_CLASS_summarised_Classif_DEF_TOTAL_\n")
  cat(str(Ov_rna_seq_index.m_CLASS_summarised_Classif_DEF_TOTAL))
  cat("\n")
  
  step<-round(max(Ov_rna_seq_index.m_CLASS_summarised_Classif_DEF_TOTAL$instances)/10,0)
  
  cat("--------->\t")
  cat(sprintf(as.character(step)))
  cat("\n")
  
  if(step == 0){
    
    step=1
  }
  
  
  breaks.y<-seq(0,max(Ov_rna_seq_index.m_CLASS_summarised_Classif_DEF_TOTAL$instances)+step, by=step)
  labels.y<-as.character(breaks.y)
  
  
  cat(sprintf(as.character(breaks.y)))
  cat("\n")
  
  graph<-ggplot(data=Ov_rna_seq_index.m_CLASS_summarised_Classif_DEF_TOTAL,
                aes(x=CLASS_string,
                    y=instances)) +
    geom_bar(stat="identity",colour='black')+
    theme_bw()+
    theme_classic()+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
          axis.text.x=element_text(angle=90,vjust=1,hjust=1,size=8, color="black", family="sans"))+
    scale_y_continuous(name='Number of variants',breaks=breaks.y,labels=labels.y,
                       limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
    scale_x_discrete(name=NULL, drop=F)+
    ggeasy::easy_center_title()
  
  
  #### subgraph ----
  
  subgraph.dt<-data.table(Ov_rna_seq_index.m,key=c('VALUE'))
  
  cat("subgraph.dt\n")
  cat(str(subgraph.dt))
  cat("\n")
  
  
  subgraph_CLASS_summarised<-as.data.frame(subgraph.dt[,.(instances=.N), by=key(subgraph.dt)],stringsAsFactors=F)
  
  subgraph_CLASS_summarised$VALUE_2<-factor(subgraph_CLASS_summarised$VALUE,
                                          levels=rev(c('MPRA_hit','no_MPRA_hit',
                                                   'genIE_hit','no_genIE_hit','NOT_SCREENED_genIE',
                                                   'Transcriptional R.','Transcriptional R. + ATU','ATU','NRD',
                                                   'R_in_candidate','R_in_non_candidate','Regulatory_proxy','Coding_proxy','Other_tissue_needed','Unknown_mechanism')),
                                          ordered=T)
  
 
  cat("subgraph_CLASS_summarised\n")
  cat(str(subgraph_CLASS_summarised))
  cat("\n")
  cat(sprintf(as.character(names(summary(subgraph_CLASS_summarised$VALUE_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(subgraph_CLASS_summarised$VALUE_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(subgraph_CLASS_summarised$instances)))))
  cat("\n")
  cat(sprintf(as.character(summary(subgraph_CLASS_summarised$instances))))
  cat("\n")
  
  check<-subgraph_CLASS_summarised[is.na(subgraph_CLASS_summarised$VALUE_2),]
  
  cat("check\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$VALUE_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$VALUE_2))))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$instances)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$instances))))
  cat("\n")
  
  step<-round(max(subgraph_CLASS_summarised$instances)/5,0)
  
  cat("--------->\t")
  cat(sprintf(as.character(step)))
  cat("\n")
  
  if(step == 0){
    
    step=1
  }
  breaks.x<-rev(c(seq(0,max(subgraph_CLASS_summarised$instances)+step, by=step)))
  breaks.x<-rev(c(seq(0,50+1, by=1)))
  labels.x<-as.character(breaks.x)
  
  
  cat(sprintf(as.character(breaks.x)))
  cat("\n")
  
 
  
  
  subgraph<-ggplot(data=subgraph_CLASS_summarised[which(subgraph_CLASS_summarised$VALUE_2 != "NOT_SCREENED_genIE"),],
                   aes(y=VALUE_2,
                       x=instances,
                       fill=VALUE_2)) +
    geom_bar(stat="identity",colour='black')+
    theme_bw()+
    theme_classic()+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_text(angle=0,size=8, color="black", family="sans"),
          axis.text.x=element_text(angle=0,size=8, color="black", family="sans"))+
    scale_x_reverse(name="Number of variant",breaks=breaks.x,labels=labels.x,
                    limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
    scale_y_discrete(position="right",name=NULL, drop=T)+
    theme(legend.position="hidden")+
    ggeasy::easy_center_title()
  
  cat("subgraph DONE\n")
  
  
  
  
  
  graph_FINAL<-plot_grid(NULL,graph,subgraph,NULL,
                         nrow = 2,
                         ncol=2,
                         rel_heights = c(1, 0.25),
                         rel_widths=c(0.5,1))
  
  
 #### Main table susbstitutions -----
  
  
  indx.int<-c(which(colnames(Table_S4) == "rs"),which(colnames(Table_S4) == "VAR"),which(colnames(Table_S4) == "Mechanistic_Class"),which(colnames(Table_S4) == "Candidate_effector"),
              which(colnames(Table_S4) == "Multi_Lineage"),which(colnames(Table_S4) == "phenotype_DEF_string"),which(colnames(Table_S4) == "Lineage_string"))
  
 Main_table<-unique(Table_S4[which(Table_S4$Manual_curation == 'R_in_candidate'),indx.int])
 
 cat("Main_table_0\n")
 cat(str(Main_table))
 cat("\n")
 cat(str(unique(Main_table$VAR)))
 cat("\n")
 
 Main_table$VAR<-sub("_",":",Main_table$VAR)
 Main_table$VAR<-sub("_"," ",Main_table$VAR)
 Main_table$VAR<-sub("_",'>',Main_table$VAR)
 
 Main_table$Multi_Lineage<-gsub("_",' ',Main_table$Multi_Lineage)
 Main_table$Lineage_string<-gsub("_",' ',Main_table$Lineage_string)
 
 
 colnames(Main_table)[which(colnames(Main_table) == 'rs')]<-'Rsid'
 colnames(Main_table)[which(colnames(Main_table) == 'VAR')]<-'chr:pos ref > alt'
 colnames(Main_table)[which(colnames(Main_table) == 'Mechanistic_Class')]<-'Mechanism'
 colnames(Main_table)[which(colnames(Main_table) == 'Candidate_effector')]<-'Regulated candidate gene'
 colnames(Main_table)[which(colnames(Main_table) == 'Multi_Lineage')]<-'Lineage classification'
 colnames(Main_table)[which(colnames(Main_table) == 'phenotype_DEF_string')]<-'Phenotypes associated in GWAS'
 colnames(Main_table)[which(colnames(Main_table) == 'Lineage_string')]<-'Lineages associated'
 
 
 cat("Main_table_1\n")
 cat(str(Main_table))
 cat("\n")
 
  
  # ############## SAVE -------------------------
  # 
  setwd(out)
  
  
  pdf(paste('Breakdown_Screening_mechanism_and_manual_curation','.pdf',sep=''))
  print(graph)
  print(subgraph)
  dev.off()
  
  cat("Print tables\n")
  
  write.table(subgraph_CLASS_summarised,file='test_table.tsv', sep="\t", row.names = F,  quote=F)
  
  write.table(Ov_rna_seq_index.m_CLASS_summarised_Classif_DEF_TOTAL,file='test_table_2.tsv', sep="\t", row.names = F,  quote=F)
  
  Table_S4_for_printing<-Table_S4
  
  Table_S4_for_printing$Multi_Lineage<-gsub("_",' ',Table_S4_for_printing$Multi_Lineage)
  Table_S4_for_printing$Lineage_string<-gsub("_",' ',Table_S4_for_printing$Lineage_string)
  
  Table_S4_for_printing$Manual_curation<-gsub("_",' ',Table_S4_for_printing$Manual_curation)
  Table_S4_for_printing$MPRA_CLASS_2<-gsub("_",' ',Table_S4_for_printing$MPRA_CLASS_2)
  Table_S4_for_printing$genIE_CLASS_2<-gsub("_",' ',Table_S4_for_printing$genIE_CLASS_2)
  
  
  indx.dep<-c(which(colnames(Table_S4_for_printing) == 'MPRA_CLASS'),
              which(colnames(Table_S4_for_printing) == 'genIE_CLASS'))
  
  Table_S4_for_printing<-unique(Table_S4_for_printing[,-indx.dep])
  
  colnames(Table_S4_for_printing)[which(colnames(Table_S4_for_printing) == 'MPRA_CLASS_2')]<-'MPRA_CLASS'
  colnames(Table_S4_for_printing)[which(colnames(Table_S4_for_printing) == 'genIE_CLASS_2')]<-'genIE_CLASS'
  
  indx.reorder<-c(which(colnames(Table_S4_for_printing) == 'VAR'),
                  which(colnames(Table_S4_for_printing) == 'rs'),
                  which(colnames(Table_S4_for_printing) == 'Mechanistic_Class'),
                  which(colnames(Table_S4_for_printing) == 'Manual_curation'),
                  which(colnames(Table_S4_for_printing) == 'Candidate_effector'),
                  which(colnames(Table_S4_for_printing) == 'OpenTargets_QTL'),
                  which(colnames(Table_S4_for_printing) == 'Whole_blood_DE_HGNC_string'),
                  which(colnames(Table_S4_for_printing) == 'Whole_blood_DTU_HGNC_string'),
                  which(colnames(Table_S4_for_printing) == 'Monocyte_DE_HGNC_string'),
                  which(colnames(Table_S4_for_printing) == 'Tcell_DE_HGNC_string'),
                  which(colnames(Table_S4_for_printing) == 'Neutrophil_DE_HGNC_string'),
                  which(colnames(Table_S4_for_printing) == 'Neutrophil_DTU_HGNC_string'),
                  which(colnames(Table_S4_for_printing) == 'Tcell_DTU_HGNC_string'),
                  which(colnames(Table_S4_for_printing) == 'Monocyte_DTU_HGNC_string'),
                  which(colnames(Table_S4_for_printing) == 'Proxy_rsid_string'),
                  which(colnames(Table_S4_for_printing) == 'MPRA_CLASS'),
                  which(colnames(Table_S4_for_printing) == 'genIE_CLASS'),
                  which(colnames(Table_S4_for_printing) == 'Replication_OT_QTL'),
                  which(colnames(Table_S4_for_printing) == 'VAR_38'),
                  which(colnames(Table_S4_for_printing) == 'Multi_Lineage'),
                  which(colnames(Table_S4_for_printing) == 'Lineage_string'),
                  which(colnames(Table_S4_for_printing) == 'phenotype_DEF_string'),
                  which(colnames(Table_S4_for_printing) == 'maf_origin'),
                  which(colnames(Table_S4_for_printing) == 'VEP_DEF_LABELS_wCSQ'))
  
  Table_S4_for_printing<-Table_S4_for_printing[,indx.reorder]         
  
 
  
  cat("Table_S4_for_printing_0\n")
  cat(str(Table_S4_for_printing))
  cat("\n")
  cat(str(unique(Table_S4_for_printing$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_S4_for_printing$MPRA_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S4_for_printing$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_S4_for_printing$genIE_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S4_for_printing$genIE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_for_printing$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_for_printing$Mechanistic_Class))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_S4_for_printing$Manual_curation))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S4_for_printing$Manual_curation)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_S4_for_printing$Multi_Lineage))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S4_for_printing$Multi_Lineage)))))
  cat("\n")
  
  Table_S4_for_printing$MPRA_CLASS<-factor(Table_S4_for_printing$MPRA_CLASS,
                                           levels=c('MPRA hit','no MPRA hit','NOT SCREENED MPRA'),
                                           ordered=T)
  
  Table_S4_for_printing$genIE_CLASS<-factor(Table_S4_for_printing$genIE_CLASS,
                                           levels=c('genIE hit','no genIE hit','NOT SCREENED genIE'),
                                           ordered=T)
  
  Table_S4_for_printing$Manual_curation<-factor(Table_S4_for_printing$Manual_curation,
                                           levels=c('R in candidate','R in non candidate','Coding proxy','Regulatory proxy','Other tissue needed','Unknown mechanism','No RNA Seq HET carriers'),
                                           ordered=T)
  
  Table_S4_for_printing$Multi_Lineage<-factor(Table_S4_for_printing$Multi_Lineage,
                                                levels=c('Multi Lineage','Lineage restricted'),
                                                ordered=T)
  
  cat("Table_S4_for_printing_1\n")
  cat(str(Table_S4_for_printing))
  cat("\n")
  cat(str(unique(Table_S4_for_printing$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_for_printing$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_for_printing$MPRA_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_for_printing$genIE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_for_printing$genIE_CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_for_printing$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_for_printing$Mechanistic_Class))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_for_printing$Manual_curation)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_for_printing$Manual_curation))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S4_for_printing$Multi_Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S4_for_printing$Multi_Lineage))))
  cat("\n")
  
  setwd(out)
  
  write.table(Table_S4_for_printing,file='Table_S6_Manual_curation.tsv', sep="\t", row.names = F,  quote=F)
  
  write.table(Main_table,file='Main_table.tsv', sep="\t", row.names = F,  quote=F)
  
  saveRDS(Table_S4_for_printing, file='Table_S6_Manual_curation.rds')
  
 
  
}

Barplot_per_manual_category = function(option_list)
{
  
  library("RColorBrewer",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
  
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
  
  setwd(out)
  
  Table_S6<-readRDS(file='Table_S6_Manual_curation.rds')
  
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
  
  
  #### Freq table S6  ----
  
  Condition_DEBUG <- 1
  
  
  Table_S6.dt<-data.table(Table_S6, key=c("Mechanistic_Class","Manual_curation"))
  
  Freq_instances<-as.data.frame(Table_S6.dt[,.(instances=.N),by=key(Table_S6.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_instances_0\n")
    cat(str(Freq_instances))
    cat("\n")
    #quit(status = 1)
  }
  
  
  Table_S6.dt<-data.table(Table_S6, key=c("Mechanistic_Class"))
  
  Freq_TOTAL<-as.data.frame(Table_S6.dt[,.(TOTAL=.N),by=key(Table_S6.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_instances<-merge(Freq_instances,
                        Freq_TOTAL,
                        by=c("Mechanistic_Class"))
  
  # Freq_instances$Mechanistic_Class<-factor(Freq_instances$Mechanistic_Class,
  #                                          levels_Mechanistic_Class,
  #                                          ordered=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_instances_1\n")
    cat(str(Freq_instances))
    cat("\n")
    #quit(status = 1)
  }
  
  vector_colors_MandM<-brewer.pal(length(levels_Manual_curation), "Set1")
  
  vector_colors_MandM<-c(vector_colors_MandM[3],vector_colors_MandM[2],vector_colors_MandM[1],vector_colors_MandM[4:length(vector_colors_MandM)])
  
  A<-summary(Freq_TOTAL$TOTAL)
  
  breaks.Rank<-seq(0,max(A)+2, by=5)
  labels.Rank<-as.character(breaks.Rank)
  
  if(Condition_DEBUG == 1)
  {
    cat("labels.Rank\n")
    cat(str(labels.Rank))
    cat("\n")
    #quit(status = 1)
  }
  
  MandM_graph<-ggplot(data=Freq_instances,
                            aes(x=Mechanistic_Class, y=instances,
                                fill=Manual_curation)) +
    geom_bar(stat="identity",colour='black')+
    theme_classic()+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=14,hjust=1,vjust=1, color="black", family="sans"))+
    scale_y_continuous(name=paste("# variants",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    scale_x_discrete(name=NULL)+
    theme(legend.position="right",legend.title=element_blank(), legend.text = element_text(size=10))+
    scale_fill_manual(values=vector_colors_MandM)+
    ggeasy::easy_center_title()
  
  setwd(out)
  
  svglite(paste('MandM_graph','.svg',sep=''), width = 8, height = 8)
  print(MandM_graph)
  dev.off()
  
  cat("MandM_graph DONE\n")

  
  write.table(Freq_instances, file='Frequency_table_MandM_barplot.tsv', sep="\t", quote=F, row.names = F)
  
}

Barplot_per_manual_category_and_MPRA = function(option_list)
{
  
  library("RColorBrewer",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
  
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
  
  setwd(out)
  
  Table_S6<-readRDS(file='Table_S6_Manual_curation.rds')
  
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
  
  
  #### Freq table S6  ----
  
  Condition_DEBUG <- 1
  
  
  Table_S6.dt<-data.table(Table_S6, key=c("Mechanistic_Class","Manual_curation","MPRA_CLASS"))
  
  Freq_instances<-as.data.frame(Table_S6.dt[,.(instances=.N),by=key(Table_S6.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_instances_0\n")
    cat(str(Freq_instances))
    cat("\n")
    #quit(status = 1)
  }
  
  
  Table_S6.dt<-data.table(Table_S6, key=c("Mechanistic_Class","Manual_curation"))
  
  Freq_TOTAL<-as.data.frame(Table_S6.dt[,.(TOTAL=.N),by=key(Table_S6.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_instances<-merge(Freq_instances,
                        Freq_TOTAL,
                        by=c("Mechanistic_Class","Manual_curation"))
  
  Freq_instances$Perc<-100*(Freq_instances$instances/Freq_instances$TOTAL)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_instances_1\n")
    cat(str(Freq_instances))
    cat("\n")
  }
  

  vector_colors_MPRA<-brewer.pal(3, "Blues")
  vector_colors_MPRA<-c(vector_colors_MPRA[3],vector_colors_MPRA[1],vector_colors_MPRA[2])
  

  # A<-summary(Freq_TOTAL$TOTAL)
  
  breaks.Rank<-seq(0,100, by=10)
  labels.Rank<-as.character(breaks.Rank)
  
  if(Condition_DEBUG == 1)
  {
    cat("labels.Rank\n")
    cat(str(labels.Rank))
    cat("\n")
    #quit(status = 1)
  }
  
  MPRA_graph<-ggplot(data=Freq_instances,
                      aes(x=Manual_curation, y=Perc,
                          fill=MPRA_CLASS)) +
    geom_bar(stat="identity",colour='black')+
    theme_classic()+
    scale_fill_manual(values=vector_colors_MPRA)+
    ggeasy::easy_center_title()
  
  MPRA_graph<-MPRA_graph+
    facet_grid(. ~ Mechanistic_Class, drop=F)+
    theme_cowplot(font_size = 14)+
    theme( strip.background = element_blank(),
           strip.placement = "inside",
           strip.text = element_text(size=14),
           panel.spacing = unit(0.2, "lines"), 
           panel.background=element_rect(fill="white"),
           panel.border=element_rect(colour="black",size=1),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
    theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
          legend.key.height = unit(1.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=14))+ #change legend text font size
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=14,hjust=1,vjust=1, color="black", family="sans"))+
    scale_y_continuous(name=paste("# variants",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    theme(legend.position="right",legend.title=element_blank(), legend.text = element_text(size=10))+
    scale_x_discrete(name=NULL, drop=T)
  
  
  setwd(out)
  
  svglite(paste('MPRA_graph','.svg',sep=''), width = 8, height = 8)
  print(MPRA_graph)
  dev.off()
  
  cat("MPRA_graph DONE\n")
  
  
  write.table(Freq_instances, file='Frequency_table_MandM_MPRA_barplot.tsv', sep="\t", quote=F, row.names = F)
  
}

Stats_per_manual_category_and_MPRA = function(option_list)
{
  
  library("RColorBrewer",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
  
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
  
  setwd(out)
  
  Table_S6<-readRDS(file='Table_S6_Manual_curation.rds')
  
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
  
  
  #### STATS table S6  ----
  
  Condition_DEBUG <- 1
  
  
  Table_S6.dt<-data.table(Table_S6, key=c("Mechanistic_Class","Manual_curation","MPRA_CLASS"))
  
  Freq_instances<-as.data.frame(Table_S6.dt[,.(instances=.N),by=key(Table_S6.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_instances_0\n")
    cat(str(Freq_instances))
    cat("\n")
    #quit(status = 1)
  }
  
  
  Table_S6.dt<-data.table(Table_S6, key=c("Mechanistic_Class","Manual_curation"))
  
  Freq_TOTAL<-as.data.frame(Table_S6.dt[,.(TOTAL=.N),by=key(Table_S6.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_instances<-merge(Freq_instances,
                        Freq_TOTAL,
                        by=c("Mechanistic_Class","Manual_curation"))
  
  Freq_instances$Perc<-100*(Freq_instances$instances/Freq_instances$TOTAL)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_instances_1\n")
    cat(str(Freq_instances))
    cat("\n")
  }
  
  list_levels_Mechanistic_Class<-list()
  
 for(i in 1:length(levels_Mechanistic_Class))
 {
   levels_Mechanistic_Class_sel<-levels_Mechanistic_Class[i]
   
   cat("--------------------------------------------------------------------------------->\t")
   cat(sprintf(as.character(levels_Mechanistic_Class_sel)))
   cat("\n")
   
   Freq_instances_SELECTED_LEVEL<-Freq_instances[which(Freq_instances$Mechanistic_Class == levels_Mechanistic_Class_sel),]
   
   
   if(Condition_DEBUG == 1)
   {
     cat("Freq_instances_SELECTED_LEVEL_1\n")
     cat(str(Freq_instances_SELECTED_LEVEL))
     cat("\n")
   }
   
   if(dim(Freq_instances_SELECTED_LEVEL)[1] >0)
   {
     list_levels_Manual_curation<-list()
     
     for(k in 1:length(levels_Manual_curation))
     {
       levels_Manual_curation_sel<-levels_Manual_curation[k]
       
       cat("---------------------------------->\t")
       cat(sprintf(as.character(levels_Manual_curation_sel)))
       cat("\n")
       
       Freq_instances_SELECTED_LEVEL_MANUAL_SEL<-Freq_instances_SELECTED_LEVEL[which(Freq_instances_SELECTED_LEVEL$Manual_curation == levels_Manual_curation_sel),]
       
       
       if(Condition_DEBUG == 1)
       {
         cat("Freq_instances_SELECTED_LEVEL_MANUAL_SEL_1\n")
         cat(str(Freq_instances_SELECTED_LEVEL_MANUAL_SEL))
         cat("\n")
       }
       
       if(dim(Freq_instances_SELECTED_LEVEL_MANUAL_SEL)[1] >0)
       {
         
         MPRA_hit_BAIT<-0
         
         if(length(Freq_instances_SELECTED_LEVEL_MANUAL_SEL$MPRA_CLASS[which(Freq_instances_SELECTED_LEVEL_MANUAL_SEL$MPRA_CLASS == 'MPRA hit')]) >0)
         {
           
           MPRA_hit_BAIT<-as.numeric(Freq_instances_SELECTED_LEVEL_MANUAL_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MANUAL_SEL$MPRA_CLASS == 'MPRA hit')])
         }
         
         MPRA_no_hit_BAIT<-0
         
         if(length(Freq_instances_SELECTED_LEVEL_MANUAL_SEL$MPRA_CLASS[which(Freq_instances_SELECTED_LEVEL_MANUAL_SEL$MPRA_CLASS == 'no MPRA hit')]) >0)
         {
           
           MPRA_no_hit_BAIT<-as.numeric(Freq_instances_SELECTED_LEVEL_MANUAL_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MANUAL_SEL$MPRA_CLASS == 'no MPRA hit')])
         }
         
         if(Condition_DEBUG == 1)
         {
           cat("MPRA_hit_BAIT_1\n")
           cat(str(MPRA_hit_BAIT))
           cat("\n")
           
           cat("MPRA_no_hit_BAIT_1\n")
           cat(str(MPRA_no_hit_BAIT))
           cat("\n")
         }
         
         
         REST<-droplevels(Freq_instances[-which(Freq_instances$Mechanistic_Class == levels_Mechanistic_Class_sel &
                                       Freq_instances$Manual_curation ==  levels_Manual_curation_sel),])
         
         if(Condition_DEBUG == 1)
         {
           cat("REST_1\n")
           cat(str(REST))
           cat("\n")
         }
         
         
         REST_levels_Mechanistic_Class<-levels(REST$Mechanistic_Class)
         
         cat("REST_levels_Mechanistic_Class\n")
         cat(str(REST_levels_Mechanistic_Class))
         cat("\n")
         
         REST_levels_Manual_curation<-levels(REST$Manual_curation)
         
         cat("REST_levels_Manual_curation\n")
         cat(str(REST_levels_Manual_curation))
         cat("\n")
         
         list_REST_SELECTED_LEVEL<-list()
         
         for(l in 1:length(REST_levels_Mechanistic_Class))
         {
           REST_levels_Mechanistic_Class_sel<-REST_levels_Mechanistic_Class[l]
           
           cat("----------->\t")
           cat(sprintf(as.character(REST_levels_Mechanistic_Class_sel)))
           cat("\n")
           
           REST_SELECTED_LEVEL<-REST[which(REST$Mechanistic_Class == REST_levels_Mechanistic_Class_sel),]
           
           
           if(Condition_DEBUG == 1)
           {
             cat("REST_SELECTED_LEVEL_0\n")
             cat(str(REST_SELECTED_LEVEL))
             cat("\n")
           }
           
           if(dim(REST_SELECTED_LEVEL)[1] >0)
           {
             
             list_REST_levels_Manual_curation<-list()
             
             for(m in 1:length(REST_levels_Manual_curation))
             {
               REST_levels_Manual_curation_sel<-REST_levels_Manual_curation[m]
               
               cat("---->\t")
               cat(sprintf(as.character(REST_levels_Manual_curation_sel)))
               cat("\n")
               
               REST_SELECTED_LEVEL_MANUAL<-REST_SELECTED_LEVEL[which(REST_SELECTED_LEVEL$Manual_curation == REST_levels_Manual_curation_sel),]
               
               
               if(Condition_DEBUG == 1)
               {
                 cat("REST_SELECTED_LEVEL_MANUAL_0\n")
                 cat(str(REST_SELECTED_LEVEL_MANUAL))
                 cat("\n")
               }
               
               if(dim(REST_SELECTED_LEVEL_MANUAL)[1] >0)
               {
                 
                 MPRA_hit_PREY<-0
                 
                 if(length(REST_SELECTED_LEVEL_MANUAL$MPRA_CLASS[which(REST_SELECTED_LEVEL_MANUAL$MPRA_CLASS == 'MPRA hit')]) >0)
                 {
                   
                   MPRA_hit_PREY<-as.numeric(REST_SELECTED_LEVEL_MANUAL$instances[which(REST_SELECTED_LEVEL_MANUAL$MPRA_CLASS == 'MPRA hit')])
                 }
                 
                 MPRA_no_hit_PREY<-0
                 
                 if(length(REST_SELECTED_LEVEL_MANUAL$MPRA_CLASS[which(REST_SELECTED_LEVEL_MANUAL$MPRA_CLASS == 'no MPRA hit')]) >0)
                 {
                   
                   MPRA_no_hit_PREY<-as.numeric(REST_SELECTED_LEVEL_MANUAL$instances[which(REST_SELECTED_LEVEL_MANUAL$MPRA_CLASS == 'no MPRA hit')])
                 }
                 
                 if(Condition_DEBUG == 1)
                 {
                   cat("MPRA_hit_PREY_1\n")
                   cat(str(MPRA_hit_PREY))
                   cat("\n")
                   
                   cat("MPRA_no_hit_PREY_1\n")
                   cat(str(MPRA_no_hit_PREY))
                   cat("\n")
                 }
                 
                 A.df<-as.data.frame(rbind(cbind(MPRA_no_hit_BAIT,MPRA_hit_BAIT),
                                           cbind(MPRA_no_hit_PREY,MPRA_hit_PREY)))
                 
                 
                 colnames(A.df)<-c("MPRA_inactive","MPRA_active")
                 
                 tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                 
                 if(Condition_DEBUG == 1)
                 {
                   # cat("tab.chisq.test\n")
                   # cat(str(tab.chisq.test))
                   # cat("\n")
                 }
                 
                 pval<-as.numeric(tab.chisq.test$p.value)
                 log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                 
                 if(Condition_DEBUG == 1)
                 {
                   cat("log_pval\n")
                   cat(str(log_pval))
                   cat("\n")
                 }
                   
                 
                 comparison_df<-as.data.frame(cbind(levels_Mechanistic_Class_sel,levels_Manual_curation_sel,REST_levels_Mechanistic_Class_sel,REST_levels_Manual_curation_sel,
                                                    MPRA_no_hit_BAIT,MPRA_hit_BAIT,MPRA_no_hit_PREY,MPRA_hit_PREY,log_pval), stringsAsFactors=F)
                 colnames(comparison_df)<-c("MECH_c1","CURATION_c1","MECH_c2","CURATION_c2","MPRA_OTHER_LEV_c1","MPRA_E_Plus_ASE_c1","MPRA_OTHER_LEV_c2","MPRA_E_Plus_ASE_c2","minuslogpval")
                 
                 if(Condition_DEBUG == 1)
                 {
                   cat("comparison_df\n")
                   cat(str(comparison_df))
                   cat("\n")
                 }
                 
                 list_REST_levels_Manual_curation[[m]]<-comparison_df
                 # quit(status = 1)
                 
               }#dim(REST_SELECTED_LEVEL_MANUAL)[1] >0
             }#m in 1:length(REST_levels_Manual_curation)
             
             if(length(list_REST_levels_Manual_curation) >0)
             {
               df_REST_levels_Manual_curation <- as.data.frame(data.table::rbindlist(list_REST_levels_Manual_curation, fill=T), stringsAsFactors=F)
               
               if(Condition_DEBUG == 1)
               {
                 cat("df_REST_levels_Manual_curation\n")
                 cat(str(df_REST_levels_Manual_curation))
                 cat("\n")
               }
               
               list_REST_SELECTED_LEVEL[[l]]<-df_REST_levels_Manual_curation 
             }
             
           }#dim(REST_SELECTED_LEVEL)[1] >0
         }#l in 1:length(REST_levels_Mechanistic_Class)
         
        
         if(length(list_REST_SELECTED_LEVEL) >0)
         {
           df_REST_SELECTED_LEVEL <- as.data.frame(data.table::rbindlist(list_REST_SELECTED_LEVEL, fill=T), stringsAsFactors=F)
           
           if(Condition_DEBUG == 1)
           {
             cat("df_REST_SELECTED_LEVEL\n")
             cat(str(df_REST_SELECTED_LEVEL))
             cat("\n")
           }
           
           list_levels_Manual_curation[[k]]<-df_REST_SELECTED_LEVEL 
           
           
         }#length(list_REST_SELECTED_LEVEL) >0
         
       }# dim(Freq_instances_SELECTED_LEVEL_MANUAL_SEL)[1] >0
     }#k in 1:length(levels_Manual_curation)
     
     if(length(list_levels_Manual_curation) >0)
     {
       df_REST_SELECTED_LEVEL <- as.data.frame(data.table::rbindlist(list_levels_Manual_curation, fill=T), stringsAsFactors=F)
       
       if(Condition_DEBUG == 1)
       {
         cat("df_REST_SELECTED_LEVEL\n")
         cat(str(df_REST_SELECTED_LEVEL))
         cat("\n")
       }
       
       list_levels_Mechanistic_Class[[i]]<-df_REST_SELECTED_LEVEL 
       
       
     }#length(list_levels_Manual_curation) >0
     
   }#dim(Freq_instances_SELECTED_LEVEL)[1] >0
 }#i in 1:length(levels_Manual_curation)
  
  if(length(list_levels_Mechanistic_Class) >0)
  {
    df_FINAL <- as.data.frame(data.table::rbindlist(list_levels_Mechanistic_Class, fill=T), stringsAsFactors=F)
    
    if(Condition_DEBUG == 1)
    {
      cat("df_FINAL\n")
      cat(str(df_FINAL))
      cat("\n")
    }
    
    setwd(out)
    
    write.table(df_FINAL, file='manual_category_and_MPRA_stats.tsv', sep="\t", quote=F, row.names = F)
    
  }#length(list_levels_Mechanistic_Class) >0
  
 
  
  
}

Barplot_MPRA_per_manual_category = function(option_list)
{
  
  library("RColorBrewer",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
  
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
  
  setwd(out)
  
  Table_S6<-readRDS(file='Table_S6_Manual_curation.rds')
  
  Table_S6<-droplevels(Table_S6[-which(Table_S6$MPRA_CLASS == 'NOT SCREENED MPRA'),])
  
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
  
  levels_MPRA_CLASS<-rev(levels(Table_S6$MPRA_CLASS))
  
  Table_S6$MPRA_CLASS<-factor(Table_S6$MPRA_CLASS,
                              levels=levels_MPRA_CLASS,
                              ordered=T)
  
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
  
  
  #### Freq table S6  ----
  
  Condition_DEBUG <- 1
  
  
  Table_S6.dt<-data.table(Table_S6, key=c("MPRA_CLASS","Mechanistic_Class","Manual_curation"))
  
  Freq_instances<-as.data.frame(Table_S6.dt[,.(instances=.N),by=key(Table_S6.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_instances_0\n")
    cat(str(Freq_instances))
    cat("\n")
    #quit(status = 1)
  }
  
  
  Table_S6.dt<-data.table(Table_S6, key=c("MPRA_CLASS","Mechanistic_Class"))
  
  Freq_TOTAL<-as.data.frame(Table_S6.dt[,.(TOTAL=.N),by=key(Table_S6.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_instances<-merge(Freq_instances,
                        Freq_TOTAL,
                        by=c("MPRA_CLASS","Mechanistic_Class"))
  
  # Freq_instances$MPRA_CLASS<-factor(Freq_instances$MPRA_CLASS,
  #                                          levels_MPRA_CLASS,
  #                                          ordered=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_instances_1\n")
    cat(str(Freq_instances))
    cat("\n")
    #quit(status = 1)
  }
  
  vector_colors_MPRA_and_manual<-brewer.pal(length(levels_Manual_curation), "Set1")
  
  vector_colors_MPRA_and_manual<-c(vector_colors_MPRA_and_manual[3],vector_colors_MPRA_and_manual[2],vector_colors_MPRA_and_manual[1],vector_colors_MPRA_and_manual[4:length(vector_colors_MPRA_and_manual)])
  
  A<-summary(Freq_TOTAL$TOTAL)
  
  breaks.Rank<-seq(0,max(A)+1, by=1)
  labels.Rank<-as.character(breaks.Rank)
  
  if(Condition_DEBUG == 1)
  {
    cat("labels.Rank\n")
    cat(str(labels.Rank))
    cat("\n")
    #quit(status = 1)
  }
  
  MPRA_and_manual<-ggplot(data=Freq_instances,
                      aes(x=MPRA_CLASS, y=instances,
                          fill=Manual_curation)) +
    geom_bar(stat="identity",colour='black')+
    theme_classic()+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=14,hjust=1,vjust=1, color="black", family="sans"))+
    scale_y_continuous(name=paste("# variants",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    scale_x_discrete(name=NULL)+
    theme(legend.position="right",legend.title=element_blank(), legend.text = element_text(size=10))+
    scale_fill_manual(values=vector_colors_MPRA_and_manual)+
    ggeasy::easy_center_title()
  
  MPRA_and_manual<-MPRA_and_manual+
    facet_grid(. ~ Mechanistic_Class, drop=F)+
    theme_cowplot(font_size = 14)+
    theme( strip.background = element_blank(),
           strip.placement = "inside",
           strip.text = element_text(size=14),
           panel.spacing = unit(0.2, "lines"), 
           panel.background=element_rect(fill="white"),
           panel.border=element_rect(colour="black",size=1),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
    theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
          legend.key.height = unit(1.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=14))+ #change legend text font size
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=14,hjust=1,vjust=1, color="black", family="sans"))+
    scale_y_continuous(name=paste("# variants",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    theme(legend.position="right",legend.title=element_blank(), legend.text = element_text(size=10))+
    scale_x_discrete(name=NULL, drop=T)
  
  
  setwd(out)
  
  svglite(paste('MPRA_and_manual_graph','.svg',sep=''), width = 8, height = 8)
  print(MPRA_and_manual)
  dev.off()
  
  cat("MPRA_and_manual DONE\n")
  
  
  write.table(Freq_instances, file='Frequency_table_MPRA_and_manual_barplot.tsv', sep="\t", quote=F, row.names = F)
  
}

Stats_MPRA_per_manual_category = function(option_list)
{
  
  library("RColorBrewer",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
  
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
  
  setwd(out)
  
  Table_S6<-readRDS(file='Table_S6_Manual_curation.rds')
  
  Table_S6<-droplevels(Table_S6[-which(Table_S6$MPRA_CLASS == 'NOT SCREENED MPRA'),])
  
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
  
  
  
  levels_Multi_Lineage<-levels(Table_S6$Multi_Lineage)
  
  cat("levels_Multi_Lineage\n")
  cat(str(levels_Multi_Lineage))
  cat("\n")
  
  
  #### STATS table S6  ----
 
  Condition_DEBUG <- 1
  
  
  Table_S6.dt<-data.table(Table_S6, key=c("MPRA_CLASS","Mechanistic_Class","Manual_curation"))
  
  Freq_instances<-as.data.frame(Table_S6.dt[,.(instances=.N),by=key(Table_S6.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_instances_0\n")
    cat(str(Freq_instances))
    cat("\n")
    #quit(status = 1)
  }
  
  
  Table_S6.dt<-data.table(Table_S6, key=c("MPRA_CLASS","Mechanistic_Class"))
  
  Freq_TOTAL<-as.data.frame(Table_S6.dt[,.(TOTAL=.N),by=key(Table_S6.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_instances<-merge(Freq_instances,
                        Freq_TOTAL,
                        by=c("MPRA_CLASS","Mechanistic_Class"))
  
  # Freq_instances$MPRA_CLASS<-factor(Freq_instances$MPRA_CLASS,
  #                                          levels_MPRA_CLASS,
  #                                          ordered=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_instances_1\n")
    cat(str(Freq_instances))
    cat("\n")
    #quit(status = 1)
  }
  
  list_levels_Mechanistic_Class<-list()
  
  for(i in 1:length(levels_Mechanistic_Class))
  {
    levels_Mechanistic_Class_sel<-levels_Mechanistic_Class[i]
    
    cat("--------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(levels_Mechanistic_Class_sel)))
    cat("\n")
    
    Freq_instances_SELECTED_LEVEL<-Freq_instances[which(Freq_instances$Mechanistic_Class == levels_Mechanistic_Class_sel),]
    
    
    if(Condition_DEBUG == 1)
    {
      cat("Freq_instances_SELECTED_LEVEL_1\n")
      cat(str(Freq_instances_SELECTED_LEVEL))
      cat("\n")
    }
    
    if(dim(Freq_instances_SELECTED_LEVEL)[1] >0)
    {
      list_levels_MPRA_CLASS<-list()
      
      for(k in 1:length(levels_MPRA_CLASS))
      {
        levels_MPRA_CLASS_sel<-levels_MPRA_CLASS[k]
        
        cat("---------------------------------->\t")
        cat(sprintf(as.character(levels_MPRA_CLASS_sel)))
        cat("\n")
        
        Freq_instances_SELECTED_LEVEL_MPRA_SEL<-Freq_instances_SELECTED_LEVEL[which(Freq_instances_SELECTED_LEVEL$MPRA_CLASS == levels_MPRA_CLASS_sel),]
        
        
        if(Condition_DEBUG == 1)
        {
          cat("Freq_instances_SELECTED_LEVEL_MPRA_SEL_1\n")
          cat(str(Freq_instances_SELECTED_LEVEL_MPRA_SEL))
          cat("\n")
        }
        
        if(dim(Freq_instances_SELECTED_LEVEL_MPRA_SEL)[1] >0)
        {
          
          R_in_candidate_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'R in candidate')]) >0)
          {
            
            R_in_candidate_BAIT<-as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'R in candidate')])
          }
          
          R_not_in_candidate_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'R in candidate')]) >0)
          {
            
            R_not_in_candidate_BAIT<-sum(as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'R in candidate')]))
          }
          
          if(Condition_DEBUG == 1)
          {
            cat("R_in_candidate_BAIT_1\n")
            cat(str(R_in_candidate_BAIT))
            cat("\n")
            
            cat("R_not_in_candidate_BAIT_1\n")
            cat(str(R_not_in_candidate_BAIT))
            cat("\n")
          }
          
          R_in_non_candidate_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'R in non candidate')]) >0)
          {
            
            R_in_non_candidate_BAIT<-as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'R in non candidate')])
          }
          
          R_not_in_non_candidate_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'R in non candidate')]) >0)
          {
            
            R_not_in_non_candidate_BAIT<-sum(as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'R in non candidate')]))
          }
          
          if(Condition_DEBUG == 1)
          {
            cat("R_in_non_candidate_BAIT_1\n")
            cat(str(R_in_non_candidate_BAIT))
            cat("\n")
            
            cat("R_not_in_non_candidate_BAIT_1\n")
            cat(str(R_not_in_non_candidate_BAIT))
            cat("\n")
          }
          
          R_in_Proxy_coding_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'Coding proxy')]) >0)
          {
            
            R_in_Proxy_coding_BAIT<-as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'Coding proxy')])
          }
          
          R_not_in_Proxy_coding_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'Coding proxy')]) >0)
          {
            
            R_not_in_Proxy_coding_BAIT<-sum(as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'Coding proxy')]))
          }
          
          if(Condition_DEBUG == 1)
          {
            cat("R_in_Proxy_coding_BAIT_1\n")
            cat(str(R_in_Proxy_coding_BAIT))
            cat("\n")
            
            cat("R_not_in_Proxy_coding_BAIT_1\n")
            cat(str(R_not_in_Proxy_coding_BAIT))
            cat("\n")
          }
          
          R_in_Proxy_regulatory_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'Regulatory proxy')]) >0)
          {
            
            R_in_Proxy_regulatory_BAIT<-as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'Regulatory proxy')])
          }
          
          R_not_in_Proxy_regulatory_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'Regulatory proxy')]) >0)
          {
            
            R_not_in_Proxy_regulatory_BAIT<-sum(as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'Regulatory proxy')]))
          }
          
          if(Condition_DEBUG == 1)
          {
            cat("R_in_Proxy_regulatory_BAIT_1\n")
            cat(str(R_in_Proxy_regulatory_BAIT))
            cat("\n")
            
            cat("R_not_in_Proxy_regulatory_BAIT_1\n")
            cat(str(R_not_in_Proxy_regulatory_BAIT))
            cat("\n")
          }
          
          
          R_in_Unknown_mechanism_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'Unknown mechanism')]) >0)
          {
            
            R_in_Unknown_mechanism_BAIT<-as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'Unknown mechanism')])
          }
          
          R_not_in_Unknown_mechanism_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'Unknown mechanism')]) >0)
          {
            
            R_not_in_Unknown_mechanism_BAIT<-sum(as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'Unknown mechanism')]))
          }
          
          if(Condition_DEBUG == 1)
          {
            cat("R_in_Unknown_mechanism_BAIT_1\n")
            cat(str(R_in_Unknown_mechanism_BAIT))
            cat("\n")
            
            cat("R_not_in_Unknown_mechanism_BAIT_1\n")
            cat(str(R_not_in_Unknown_mechanism_BAIT))
            cat("\n")
          }
          
          R_in_Other_tissue_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'Other tissue needed')]) >0)
          {
            
            R_in_Other_tissue_BAIT<-as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'Other tissue needed')])
          }
          
          R_not_in_Other_tissue_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'Other tissue needed')]) >0)
          {
            
            R_not_in_Other_tissue_BAIT<-sum(as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'Other tissue needed')]))
          }
          
          if(Condition_DEBUG == 1)
          {
            cat("R_in_Other_tissue_BAIT_1\n")
            cat(str(R_in_Other_tissue_BAIT))
            cat("\n")
            
            cat("R_not_in_Other_tissue_BAIT_1\n")
            cat(str(R_not_in_Other_tissue_BAIT))
            cat("\n")
          }
          
          R_in_No_RNA_Seq_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'No RNA Seq HET carriers')]) >0)
          {
            
            R_in_No_RNA_Seq_BAIT<-as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation == 'No RNA Seq HET carriers')])
          }
          
          R_not_in_No_RNA_Seq_BAIT<-0
          
          if(length(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'No RNA Seq HET carriers')]) >0)
          {
            
            R_not_in_No_RNA_Seq_BAIT<-sum(as.numeric(Freq_instances_SELECTED_LEVEL_MPRA_SEL$instances[which(Freq_instances_SELECTED_LEVEL_MPRA_SEL$Manual_curation != 'No RNA Seq HET carriers')]))
          }
          
          if(Condition_DEBUG == 1)
          {
            cat("R_in_No_RNA_Seq_BAIT_1\n")
            cat(str(R_in_No_RNA_Seq_BAIT))
            cat("\n")
            
            cat("R_not_in_No_RNA_Seq_BAIT_1\n")
            cat(str(R_not_in_No_RNA_Seq_BAIT))
            cat("\n")
          }
          
          REST<-droplevels(Freq_instances[-which(Freq_instances$Mechanistic_Class == levels_Mechanistic_Class_sel &
                                                   Freq_instances$MPRA_CLASS ==  levels_MPRA_CLASS_sel),])
          
          if(Condition_DEBUG == 1)
          {
            cat("REST_1\n")
            cat(str(REST))
            cat("\n")
          }
          
          
         
          
          
          REST_levels_Mechanistic_Class<-levels(REST$Mechanistic_Class)
          
          cat("REST_levels_Mechanistic_Class\n")
          cat(str(REST_levels_Mechanistic_Class))
          cat("\n")
          
          REST_levels_MPRA_CLASS<-levels(REST$MPRA_CLASS)
          
          cat("REST_levels_MPRA_CLASS\n")
          cat(str(REST_levels_MPRA_CLASS))
          cat("\n")
          
          list_REST_SELECTED_LEVEL<-list()
          
          for(l in 1:length(REST_levels_Mechanistic_Class))
          {
            REST_levels_Mechanistic_Class_sel<-REST_levels_Mechanistic_Class[l]
            
            cat("----------->\t")
            cat(sprintf(as.character(REST_levels_Mechanistic_Class_sel)))
            cat("\n")
            
            REST_SELECTED_LEVEL<-REST[which(REST$Mechanistic_Class == REST_levels_Mechanistic_Class_sel),]
            
            
            if(Condition_DEBUG == 1)
            {
              cat("REST_SELECTED_LEVEL_0\n")
              cat(str(REST_SELECTED_LEVEL))
              cat("\n")
            }
            
            if(dim(REST_SELECTED_LEVEL)[1] >0)
            {
              
              list_REST_levels_MPRA_CLASS<-list()
              
              for(m in 1:length(REST_levels_MPRA_CLASS))
              {
                REST_levels_MPRA_CLASS_sel<-REST_levels_MPRA_CLASS[m]
                
                cat("---->\t")
                cat(sprintf(as.character(REST_levels_MPRA_CLASS_sel)))
                cat("\n")
                
                REST_SELECTED_LEVEL_MPRA<-REST_SELECTED_LEVEL[which(REST_SELECTED_LEVEL$MPRA_CLASS == REST_levels_MPRA_CLASS_sel),]
                
                
                if(Condition_DEBUG == 1)
                {
                  cat("REST_SELECTED_LEVEL_MPRA_0\n")
                  cat(str(REST_SELECTED_LEVEL_MPRA))
                  cat("\n")
                }
                
                Condition_DEBUG <- 0
                
                if(dim(REST_SELECTED_LEVEL_MPRA)[1] >0)
                {
                  
                  R_in_candidate_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'R in candidate')]) >0)
                  {
                    
                    R_in_candidate_PREY<-as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'R in candidate')])
                  }
                  
                  R_not_in_candidate_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'R in candidate')]) >0)
                  {
                    
                    R_not_in_candidate_PREY<-sum(as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'R in candidate')]))
                  }
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("R_in_candidate_PREY_1\n")
                    cat(str(R_in_candidate_PREY))
                    cat("\n")
                    
                    cat("R_not_in_candidate_PREY_1\n")
                    cat(str(R_not_in_candidate_PREY))
                    cat("\n")
                  }
                  
                  R_in_non_candidate_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'R in non candidate')]) >0)
                  {
                    
                    R_in_non_candidate_PREY<-as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'R in non candidate')])
                  }
                  
                  R_not_in_non_candidate_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'R in non candidate')]) >0)
                  {
                    
                    R_not_in_non_candidate_PREY<-sum(as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'R in non candidate')]))
                  }
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("R_in_non_candidate_PREY_1\n")
                    cat(str(R_in_non_candidate_PREY))
                    cat("\n")
                    
                    cat("R_not_in_non_candidate_PREY_1\n")
                    cat(str(R_not_in_non_candidate_PREY))
                    cat("\n")
                  }
                  
                  R_in_Unknown_mechanism_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'Unknown mechanism')]) >0)
                  {
                    
                    R_in_Unknown_mechanism_PREY<-as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'Unknown mechanism')])
                  }
                  
                  R_not_in_Unknown_mechanism_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'Unknown mechanism')]) >0)
                  {
                    
                    R_not_in_Unknown_mechanism_PREY<-sum(as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'Unknown mechanism')]))
                  }
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("R_in_Unknown_mechanism_PREY_1\n")
                    cat(str(R_in_Unknown_mechanism_PREY))
                    cat("\n")
                    
                    cat("R_not_in_Unknown_mechanism_PREY_1\n")
                    cat(str(R_not_in_Unknown_mechanism_PREY))
                    cat("\n")
                  }
                  
                  R_in_Proxy_coding_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'Coding proxy')]) >0)
                  {
                    
                    R_in_Proxy_coding_PREY<-as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'Coding proxy')])
                  }
                  
                  R_not_in_Proxy_coding_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'Coding proxy')]) >0)
                  {
                    
                    R_not_in_Proxy_coding_PREY<-sum(as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'Coding proxy')]))
                  }
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("R_in_Proxy_coding_PREY_1\n")
                    cat(str(R_in_Proxy_coding_PREY))
                    cat("\n")
                    
                    cat("R_not_in_Proxy_coding_PREY_1\n")
                    cat(str(R_not_in_Proxy_coding_PREY))
                    cat("\n")
                  }
                  
                  R_in_Proxy_regulatory_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'Regulatory proxy')]) >0)
                  {
                    
                    R_in_Proxy_regulatory_PREY<-as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'Regulatory proxy')])
                  }
                  
                  R_not_in_Proxy_regulatory_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'Regulatory proxy')]) >0)
                  {
                    
                    R_not_in_Proxy_regulatory_PREY<-sum(as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'Regulatory proxy')]))
                  }
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("R_in_Proxy_regulatory_PREY_1\n")
                    cat(str(R_in_Proxy_regulatory_PREY))
                    cat("\n")
                    
                    cat("R_not_in_Proxy_regulatory_PREY_1\n")
                    cat(str(R_not_in_Proxy_regulatory_PREY))
                    cat("\n")
                  }
                  
                  R_in_Other_tissue_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'Other tissue needed')]) >0)
                  {
                    
                    R_in_Other_tissue_PREY<-as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'Other tissue needed')])
                  }
                  
                  R_not_in_Other_tissue_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'Other tissue needed')]) >0)
                  {
                    
                    R_not_in_Other_tissue_PREY<-sum(as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'Other tissue needed')]))
                  }
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("R_in_Other_tissue_PREY_1\n")
                    cat(str(R_in_Other_tissue_PREY))
                    cat("\n")
                    
                    cat("R_not_in_Other_tissue_PREY_1\n")
                    cat(str(R_not_in_Other_tissue_PREY))
                    cat("\n")
                  }
                  
                
                  R_in_No_RNA_Seq_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'No RNA Seq HET carriers')]) >0)
                  {
                    
                    R_in_No_RNA_Seq_PREY<-as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation == 'No RNA Seq HET carriers')])
                  }
                  
                  R_not_in_No_RNA_Seq_PREY<-0
                  
                  if(length(REST_SELECTED_LEVEL_MPRA$Manual_curation[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'No RNA Seq HET carriers')]) >0)
                  {
                    
                    R_not_in_No_RNA_Seq_PREY<-sum(as.numeric(REST_SELECTED_LEVEL_MPRA$instances[which(REST_SELECTED_LEVEL_MPRA$Manual_curation != 'No RNA Seq HET carriers')]))
                  }
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("R_in_No_RNA_Seq_PREY_1\n")
                    cat(str(R_in_No_RNA_Seq_PREY))
                    cat("\n")
                    
                    cat("R_not_in_No_RNA_Seq_PREY_1\n")
                    cat(str(R_not_in_No_RNA_Seq_PREY))
                    cat("\n")
                  }
                  
                  ##### all comparisons ----
                  
                  ALL_COMP<-data.frame()
                  
                  ### R_in_candidate
                  
                  A.df<-as.data.frame(rbind(cbind(R_not_in_candidate_BAIT,R_in_candidate_BAIT),
                                            cbind(R_not_in_candidate_PREY,R_in_candidate_PREY)))
                  
                  
                  colnames(A.df)<-c("R_not_in_candidate","R_in_candidate")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  if(Condition_DEBUG == 1)
                  {
                    # cat("tab.chisq.test\n")
                    # cat(str(tab.chisq.test))
                    # cat("\n")
                  }
                  
                  pval<-as.numeric(tab.chisq.test$p.value)
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("log_pval\n")
                    cat(str(log_pval))
                    cat("\n")
                  }
                  
                  
                  comparison_df<-as.data.frame(cbind(levels_Mechanistic_Class_sel,levels_MPRA_CLASS_sel,REST_levels_Mechanistic_Class_sel,REST_levels_MPRA_CLASS_sel,
                                                     R_not_in_candidate_BAIT,R_in_candidate_BAIT,R_not_in_candidate_PREY,R_in_candidate_PREY,log_pval), stringsAsFactors=F)
                  colnames(comparison_df)<-c("MECH_c1","MPRA_c1","MECH_c2","MPRA_c2","MANUAL_OTHER_LEV_c1","MANUAL_SEL_LEV_c1","MANUAL_OTHER_LEV_c2","MANUAL_SEL_LEV_c2","minuslogpval")
                  
                  comparison_df$SEL_LEVEL<-"R_in_candidate"
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("comparison_df\n")
                    cat(str(comparison_df))
                    cat("\n")
                  }
                  
                  ALL_COMP<-rbind(comparison_df,ALL_COMP)
                  
                  ### R_in_non_candidate
                  
                 
                  A.df<-as.data.frame(rbind(cbind(R_not_in_non_candidate_BAIT,R_in_non_candidate_BAIT),
                                            cbind(R_not_in_non_candidate_PREY,R_in_non_candidate_PREY)))
                  
                  
                  colnames(A.df)<-c("R_not_in_non_candidate","R_in_non_candidate")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  if(Condition_DEBUG == 1)
                  {
                    # cat("tab.chisq.test\n")
                    # cat(str(tab.chisq.test))
                    # cat("\n")
                  }
                  
                  pval<-as.numeric(tab.chisq.test$p.value)
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("log_pval\n")
                    cat(str(log_pval))
                    cat("\n")
                  }
                  
                  
                  comparison_df<-as.data.frame(cbind(levels_Mechanistic_Class_sel,levels_MPRA_CLASS_sel,REST_levels_Mechanistic_Class_sel,REST_levels_MPRA_CLASS_sel,
                                                     R_not_in_non_candidate_BAIT,R_in_non_candidate_BAIT,R_not_in_non_candidate_PREY,R_in_non_candidate_PREY,log_pval), stringsAsFactors=F)
                  colnames(comparison_df)<-c("MECH_c1","MPRA_c1","MECH_c2","MPRA_c2","MANUAL_OTHER_LEV_c1","MANUAL_SEL_LEV_c1","MANUAL_OTHER_LEV_c2","MANUAL_SEL_LEV_c2","minuslogpval")
                  
                  comparison_df$SEL_LEVEL<-"R_in_non_candidate"
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("comparison_df\n")
                    cat(str(comparison_df))
                    cat("\n")
                  }
                  
                  ALL_COMP<-rbind(comparison_df,ALL_COMP)
                  
                
                  
                  ### R_in_Proxy_coding
                  
                  A.df<-as.data.frame(rbind(cbind(R_not_in_Proxy_coding_BAIT,R_in_Proxy_coding_BAIT),
                                            cbind(R_not_in_Proxy_coding_PREY,R_in_Proxy_coding_PREY)))
                  
                  
                  colnames(A.df)<-c("R_not_in_Proxy_coding","R_in_Proxy_coding")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  if(Condition_DEBUG == 1)
                  {
                    # cat("tab.chisq.test\n")
                    # cat(str(tab.chisq.test))
                    # cat("\n")
                  }
                  
                  pval<-as.numeric(tab.chisq.test$p.value)
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("log_pval\n")
                    cat(str(log_pval))
                    cat("\n")
                  }
                  
                  
                  comparison_df<-as.data.frame(cbind(levels_Mechanistic_Class_sel,levels_MPRA_CLASS_sel,REST_levels_Mechanistic_Class_sel,REST_levels_MPRA_CLASS_sel,
                                                     R_not_in_Proxy_coding_BAIT,R_in_Proxy_coding_BAIT,R_not_in_Proxy_coding_PREY,R_in_Proxy_coding_PREY,log_pval), stringsAsFactors=F)
                  colnames(comparison_df)<-c("MECH_c1","MPRA_c1","MECH_c2","MPRA_c2","MANUAL_OTHER_LEV_c1","MANUAL_SEL_LEV_c1","MANUAL_OTHER_LEV_c2","MANUAL_SEL_LEV_c2","minuslogpval")
                  
                  comparison_df$SEL_LEVEL<-"Proxy_coding"
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("comparison_df\n")
                    cat(str(comparison_df))
                    cat("\n")
                  }
                  
                  ALL_COMP<-rbind(comparison_df,ALL_COMP)
                  
                  ### Proxy_regulatory
                  
                  A.df<-as.data.frame(rbind(cbind(R_not_in_Proxy_regulatory_BAIT,R_in_Proxy_regulatory_BAIT),
                                            cbind(R_not_in_Proxy_regulatory_PREY,R_in_Proxy_regulatory_PREY)))
                  
                  
                  colnames(A.df)<-c("R_not_Proxy_regulatory","R_Proxy_regulatory")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  if(Condition_DEBUG == 1)
                  {
                    # cat("tab.chisq.test\n")
                    # cat(str(tab.chisq.test))
                    # cat("\n")
                  }
                  
                  pval<-as.numeric(tab.chisq.test$p.value)
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("log_pval\n")
                    cat(str(log_pval))
                    cat("\n")
                  }
                  
                  
                  comparison_df<-as.data.frame(cbind(levels_Mechanistic_Class_sel,levels_MPRA_CLASS_sel,REST_levels_Mechanistic_Class_sel,REST_levels_MPRA_CLASS_sel,
                                                     R_not_in_Proxy_regulatory_BAIT,R_in_Proxy_regulatory_BAIT,R_not_in_Proxy_regulatory_PREY,R_in_Proxy_regulatory_PREY,log_pval), stringsAsFactors=F)
                  colnames(comparison_df)<-c("MECH_c1","MPRA_c1","MECH_c2","MPRA_c2","MANUAL_OTHER_LEV_c1","MANUAL_SEL_LEV_c1","MANUAL_OTHER_LEV_c2","MANUAL_SEL_LEV_c2","minuslogpval")
                  
                  comparison_df$SEL_LEVEL<-"Proxy_regulatory"
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("comparison_df\n")
                    cat(str(comparison_df))
                    cat("\n")
                  }
                  
                  ALL_COMP<-rbind(comparison_df,ALL_COMP)
                  
                  ### Other_tissue
                  
                  A.df<-as.data.frame(rbind(cbind(R_not_in_Other_tissue_BAIT,R_in_Other_tissue_BAIT),
                                            cbind(R_not_in_Other_tissue_PREY,R_in_Other_tissue_PREY)))
                  
                  
                  colnames(A.df)<-c("R_not_Other_tissue","R_Other_tissue")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  if(Condition_DEBUG == 1)
                  {
                    # cat("tab.chisq.test\n")
                    # cat(str(tab.chisq.test))
                    # cat("\n")
                  }
                  
                  pval<-as.numeric(tab.chisq.test$p.value)
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("log_pval\n")
                    cat(str(log_pval))
                    cat("\n")
                  }
                  
                  
                  comparison_df<-as.data.frame(cbind(levels_Mechanistic_Class_sel,levels_MPRA_CLASS_sel,REST_levels_Mechanistic_Class_sel,REST_levels_MPRA_CLASS_sel,
                                                     R_not_in_Other_tissue_BAIT,R_in_Other_tissue_BAIT,R_not_in_Other_tissue_PREY,R_in_Other_tissue_PREY,log_pval), stringsAsFactors=F)
                  colnames(comparison_df)<-c("MECH_c1","MPRA_c1","MECH_c2","MPRA_c2","MANUAL_OTHER_LEV_c1","MANUAL_SEL_LEV_c1","MANUAL_OTHER_LEV_c2","MANUAL_SEL_LEV_c2","minuslogpval")
                  
                  comparison_df$SEL_LEVEL<-"Other_tissue_needed"
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("comparison_df\n")
                    cat(str(comparison_df))
                    cat("\n")
                  }
                  
                  ALL_COMP<-rbind(comparison_df,ALL_COMP)
                  
                  ### Unknown_mechanism
                  
                  A.df<-as.data.frame(rbind(cbind(R_not_in_Unknown_mechanism_BAIT,R_in_Unknown_mechanism_BAIT),
                                            cbind(R_not_in_Unknown_mechanism_PREY,R_in_Unknown_mechanism_PREY)))
                  
                  
                  colnames(A.df)<-c("R_not_Unknown_mechanism","R_Unknown_mechanism")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  if(Condition_DEBUG == 1)
                  {
                    # cat("tab.chisq.test\n")
                    # cat(str(tab.chisq.test))
                    # cat("\n")
                  }
                  
                  pval<-as.numeric(tab.chisq.test$p.value)
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("log_pval\n")
                    cat(str(log_pval))
                    cat("\n")
                  }
                  
                  
                  comparison_df<-as.data.frame(cbind(levels_Mechanistic_Class_sel,levels_MPRA_CLASS_sel,REST_levels_Mechanistic_Class_sel,REST_levels_MPRA_CLASS_sel,
                                                     R_not_in_Unknown_mechanism_BAIT,R_in_Unknown_mechanism_BAIT,R_not_in_Unknown_mechanism_PREY,R_in_Unknown_mechanism_PREY,log_pval), stringsAsFactors=F)
                  colnames(comparison_df)<-c("MECH_c1","MPRA_c1","MECH_c2","MPRA_c2","MANUAL_OTHER_LEV_c1","MANUAL_SEL_LEV_c1","MANUAL_OTHER_LEV_c2","MANUAL_SEL_LEV_c2","minuslogpval")
                  
                  comparison_df$SEL_LEVEL<-"Unknown_mechanism"
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("comparison_df\n")
                    cat(str(comparison_df))
                    cat("\n")
                  }
                  
                  ALL_COMP<-rbind(comparison_df,ALL_COMP)
                  
                  
                  ### No_RNA_Seq
                  
                  A.df<-as.data.frame(rbind(cbind(R_not_in_No_RNA_Seq_BAIT,R_in_No_RNA_Seq_BAIT),
                                            cbind(R_not_in_No_RNA_Seq_PREY,R_in_No_RNA_Seq_PREY)))
                  
                  
                  colnames(A.df)<-c("R_not_No_RNA_Seq","R_No_RNA_Seq")
                  
                  tab.chisq.test<-chisq.test(A.df,correct = TRUE)
                  
                  if(Condition_DEBUG == 1)
                  {
                    # cat("tab.chisq.test\n")
                    # cat(str(tab.chisq.test))
                    # cat("\n")
                  }
                  
                  pval<-as.numeric(tab.chisq.test$p.value)
                  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("log_pval\n")
                    cat(str(log_pval))
                    cat("\n")
                  }
                  
                  
                  comparison_df<-as.data.frame(cbind(levels_Mechanistic_Class_sel,levels_MPRA_CLASS_sel,REST_levels_Mechanistic_Class_sel,REST_levels_MPRA_CLASS_sel,
                                                     R_not_in_No_RNA_Seq_BAIT,R_in_No_RNA_Seq_BAIT,R_not_in_No_RNA_Seq_PREY,R_in_No_RNA_Seq_PREY,log_pval), stringsAsFactors=F)
                  colnames(comparison_df)<-c("MECH_c1","MPRA_c1","MECH_c2","MPRA_c2","MANUAL_OTHER_LEV_c1","MANUAL_SEL_LEV_c1","MANUAL_OTHER_LEV_c2","MANUAL_SEL_LEV_c2","minuslogpval")
                  
                  comparison_df$SEL_LEVEL<-"No_RNA_Seq_carriers"
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("comparison_df\n")
                    cat(str(comparison_df))
                    cat("\n")
                  }
                  
                  ALL_COMP<-rbind(comparison_df,ALL_COMP)
                  
                  ### FINAL 
                  
                  
                  
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("ALL_COMP\n")
                    cat(str(ALL_COMP))
                    cat("\n")
                  }
                  
                  list_REST_levels_MPRA_CLASS[[m]]<-ALL_COMP
                 
                  
                }#dim(REST_SELECTED_LEVEL_MPRA)[1] >0
                
                Condition_DEBUG <- 1
              }#m in 1:length(REST_levels_MPRA_CLASS)
              
              if(length(list_REST_levels_MPRA_CLASS) >0)
              {
                df_REST_levels_MPRA_CLASS <- as.data.frame(data.table::rbindlist(list_REST_levels_MPRA_CLASS, fill=T), stringsAsFactors=F)
                
                if(Condition_DEBUG == 1)
                {
                  cat("df_REST_levels_MPRA_CLASS\n")
                  cat(str(df_REST_levels_MPRA_CLASS))
                  cat("\n")
                }
                
                list_REST_SELECTED_LEVEL[[l]]<-df_REST_levels_MPRA_CLASS 
                
               
              }
              
            }#dim(REST_SELECTED_LEVEL)[1] >0
          }#l in 1:length(REST_levels_Mechanistic_Class)
          
          
          if(length(list_REST_SELECTED_LEVEL) >0)
          {
            df_REST_SELECTED_LEVEL <- as.data.frame(data.table::rbindlist(list_REST_SELECTED_LEVEL, fill=T), stringsAsFactors=F)
            
            if(Condition_DEBUG == 1)
            {
              cat("df_REST_SELECTED_LEVEL\n")
              cat(str(df_REST_SELECTED_LEVEL))
              cat("\n")
            }
            
            list_levels_MPRA_CLASS[[k]]<-df_REST_SELECTED_LEVEL 
            
            
          }#length(list_REST_SELECTED_LEVEL) >0
          
        }# dim(Freq_instances_SELECTED_LEVEL_MPRA_SEL)[1] >0
      }#k in 1:length(levels_MPRA_CLASS)
      
      if(length(list_levels_MPRA_CLASS) >0)
      {
        df_REST_SELECTED_LEVEL <- as.data.frame(data.table::rbindlist(list_levels_MPRA_CLASS, fill=T), stringsAsFactors=F)
        
        if(Condition_DEBUG == 1)
        {
          cat("df_REST_SELECTED_LEVEL\n")
          cat(str(df_REST_SELECTED_LEVEL))
          cat("\n")
        }
        
        list_levels_Mechanistic_Class[[i]]<-df_REST_SELECTED_LEVEL 
        
        
      }#length(list_levels_MPRA_CLASS) >0
      
    }#dim(Freq_instances_SELECTED_LEVEL)[1] >0
  }#i in 1:length(levels_MPRA_CLASS)
  
  if(length(list_levels_Mechanistic_Class) >0)
  {
    df_FINAL <- as.data.frame(data.table::rbindlist(list_levels_Mechanistic_Class, fill=T), stringsAsFactors=F)
    
    if(Condition_DEBUG == 1)
    {
      cat("df_FINAL\n")
      cat(str(df_FINAL))
      cat("\n")
    }
    
    setwd(out)
    
    write.table(df_FINAL, file='MPRA_per_manual_category_stats.tsv', sep="\t", quote=F, row.names = F)
    
  }#length(list_levels_Mechanistic_Class) >0
  
  
  
  
}

Barplot_genIE_per_manual_category = function(option_list)
{
  
  library("RColorBrewer",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
  
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
  
  setwd(out)
  
  Table_S6<-readRDS(file='Table_S6_Manual_curation.rds')
  
  Table_S6<-droplevels(Table_S6[-which(Table_S6$genIE_CLASS == 'NOT SCREENED genIE'),])
  
  cat("Table_S6_0\n")
  cat(str(Table_S6))
  cat("\n")
  cat(str(unique(Table_S6$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6$genIE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6$genIE_CLASS))))
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
  
  levels_genIE_CLASS<-rev(levels(Table_S6$genIE_CLASS))
  
  Table_S6$genIE_CLASS<-factor(Table_S6$genIE_CLASS,
                              levels=levels_genIE_CLASS,
                              ordered=T)
  
  cat("levels_genIE_CLASS\n")
  cat(str(levels_genIE_CLASS))
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
  
  
  #### Freq table S6  ----
  
  Condition_DEBUG <- 1
  
  
  Table_S6.dt<-data.table(Table_S6, key=c("genIE_CLASS","Mechanistic_Class","Manual_curation"))
  
  Freq_instances<-as.data.frame(Table_S6.dt[,.(instances=.N),by=key(Table_S6.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_instances_0\n")
    cat(str(Freq_instances))
    cat("\n")
    #quit(status = 1)
  }
  
  
  Table_S6.dt<-data.table(Table_S6, key=c("genIE_CLASS","Mechanistic_Class"))
  
  Freq_TOTAL<-as.data.frame(Table_S6.dt[,.(TOTAL=.N),by=key(Table_S6.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_instances<-merge(Freq_instances,
                        Freq_TOTAL,
                        by=c("genIE_CLASS","Mechanistic_Class"))
  
  # Freq_instances$genIE_CLASS<-factor(Freq_instances$genIE_CLASS,
  #                                          levels_genIE_CLASS,
  #                                          ordered=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_instances_1\n")
    cat(str(Freq_instances))
    cat("\n")
    #quit(status = 1)
  }
  
  vector_colors_genIE_and_manual<-brewer.pal(length(levels_Manual_curation), "Set1")
  
  vector_colors_genIE_and_manual<-c(vector_colors_genIE_and_manual[3],vector_colors_genIE_and_manual[2],vector_colors_genIE_and_manual[1],vector_colors_genIE_and_manual[4:length(vector_colors_genIE_and_manual)])
  
  A<-summary(Freq_TOTAL$TOTAL)
  
  breaks.Rank<-seq(0,max(A)+1, by=1)
  labels.Rank<-as.character(breaks.Rank)
  
  if(Condition_DEBUG == 1)
  {
    cat("labels.Rank\n")
    cat(str(labels.Rank))
    cat("\n")
    #quit(status = 1)
  }
  
  genIE_and_manual<-ggplot(data=Freq_instances,
                          aes(x=genIE_CLASS, y=instances,
                              fill=Manual_curation)) +
    geom_bar(stat="identity",colour='black')+
    theme_classic()+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=14,hjust=1,vjust=1, color="black", family="sans"))+
    scale_y_continuous(name=paste("# variants",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    scale_x_discrete(name=NULL)+
    theme(legend.position="right",legend.title=element_blank(), legend.text = element_text(size=10))+
    scale_fill_manual(values=vector_colors_genIE_and_manual)+
    ggeasy::easy_center_title()
  
  genIE_and_manual<-genIE_and_manual+
    facet_grid(. ~ Mechanistic_Class, drop=F)+
    theme_cowplot(font_size = 14)+
    theme( strip.background = element_blank(),
           strip.placement = "inside",
           strip.text = element_text(size=14),
           panel.spacing = unit(0.2, "lines"), 
           panel.background=element_rect(fill="white"),
           panel.border=element_rect(colour="black",size=1),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
    theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
          legend.key.height = unit(1.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=14))+ #change legend text font size
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=14,hjust=1,vjust=1, color="black", family="sans"))+
    scale_y_continuous(name=paste("# variants",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    theme(legend.position="right",legend.title=element_blank(), legend.text = element_text(size=10))+
    scale_x_discrete(name=NULL, drop=T)
  
  
  setwd(out)
  
  svglite(paste('genIE_and_manual_graph','.svg',sep=''), width = 8, height = 8)
  print(genIE_and_manual)
  dev.off()
  
  cat("genIE_and_manual DONE\n")
  
  
  write.table(Freq_instances, file='Frequency_table_genIE_and_manual_barplot.tsv', sep="\t", quote=F, row.names = F)
  
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
    make_option(c("--VAR_Prioritization_dB"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
   
    make_option(c("--RMV_not_index"), type="character", default=NULL, 
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
  
  # Data_wrangling(opt)
  # Barplot_per_manual_category(opt)
  # Barplot_per_manual_category_and_MPRA(opt)
  # Stats_per_manual_category_and_MPRA(opt)
  Barplot_MPRA_per_manual_category(opt)
  # Stats_MPRA_per_manual_category(opt)
  Barplot_genIE_per_manual_category(opt)
  
}


###########################################################################

system.time( main() )
