
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
suppressMessages(library("ggdendro", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("digest", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")


opt = NULL

options(warn=1)


Initial_Selection = function(option_list)
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
  
  #### READ and transform path_Patrick ----
  
  path_Patrick = opt$path_Patrick
  
  cat("path_Patrick_\n")
  cat(sprintf(as.character(path_Patrick)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  #### READ and transform maf_Threshold ----
  
  maf_Threshold = opt$maf_Threshold
  
  cat("maf_Threshold\n")
  cat(sprintf(as.character(maf_Threshold)))
  cat("\n")
  
  
  
  #### READ and transform Prob_Threshold ----
  
  Prob_Threshold = opt$Prob_Threshold
  
  cat("Prob_Threshold\n")
  cat(sprintf(as.character(Prob_Threshold)))
  cat("\n")
  
  #### Read genIE_input file ----
  
  paleo_file<-as.data.frame(fread(file=opt$paleo_file,sep=",") , stringsAsFactors=F)
  
  paleo_file$VAR<-paste(paste('chr',paleo_file$chr,sep=''),paleo_file$pos,paleo_file$ref,paleo_file$alt,sep="_")
  
  cat("paleo_file\n")
  cat(str(paleo_file))
  cat("\n")
  
  
    
    
    NEW_db_unique.rare<-paleo_file[which(paleo_file$maf<= 0.01),]
    
   
    
    HF1<-NEW_db_unique.rare[(NEW_db_unique.rare$info >=0.8 &
                               !is.na(NEW_db_unique.rare$info)),]
    
    
    
    HF2<-HF1[(HF1$prob >= 0.9),]
    
    
    
    
    HF2.dt<-data.table(HF2, key="VAR")
    
    
    phenotype.parameters_paleo <- as.data.frame(HF2.dt[,.(max_FOR_TRAIT = max(beta),
                                                    min_FOR_TRAIT = min(beta),
                                                    Q1=summary(beta)[2],
                                                    Q3=summary(beta)[5]), by = phenotype])
    
    
    
    
    HF2_parameters<-merge(HF2, phenotype.parameters_paleo,by="phenotype")
    
    
    
    
    
    HF2_parameters.Q3<-HF2_parameters[(HF2_parameters$beta >= HF2_parameters$Q3),]
    HF2_parameters.Q1<-HF2_parameters[(HF2_parameters$beta <= HF2_parameters$Q1),]
    
    
    HF3<-HF2_parameters[(HF2_parameters$beta >= HF2_parameters$Q3 |
                           HF2_parameters$beta <= HF2_parameters$Q1),]
    
    
   
    HF2_excluded<-HF2_parameters[-which(HF2_parameters$VAR%in%HF3$VAR),]
    
    cat("HF2_excluded\n")
    cat(str(HF2_excluded))
    cat("\n")
    cat(str(unique(HF2_excluded$VAR)))
    cat("\n")
    
    
    # ###################################################
    # quit(status = 1)
  
 
  #### Read genIE_input file ----
  
  genIE_input<-as.data.frame(fread(file=opt$genIE_input,sep=",") , stringsAsFactors=F)
  
  # cat("genIE_input\n")
  # cat(str(genIE_input))
  # cat("\n")
  
  genIE_input_subset<-unique(genIE_input[which(genIE_input$HGNC != ""),])
  
  # cat("genIE_input_subset_0\n")
  # cat(str(genIE_input_subset))
  # cat("\n")
  # cat(sprintf(as.character(genIE_input_subset$HGNC)))
  # cat("\n")
  
  
  HGNC_tested<-c("FOXP1","CUX1","TNRC6A","SH2B3","BRAP","UGCG","BID","NBN","C2CD5","EPB41")

  genIE_input_subset<-unique(genIE_input_subset[which(genIE_input_subset$HGNC%in%HGNC_tested),])
  
  # cat("genIE_input_subset_1\n")
  # cat(str(genIE_input_subset))
  # cat("\n")
  # cat(sprintf(as.character(genIE_input_subset$HGNC)))
  # cat("\n")
  
  VARS_genIE_screened<-unique(genIE_input_subset$VAR)
  
  # cat("VARS_genIE_screened_0\n")
  # cat(str(VARS_genIE_screened))
  # cat("\n")
  
  
  # 
  # genIE_input_subset<-genIE_input_subset[which(genIE_input_subset$HGNC),]
  
  

  
  #### Read MPRA_Tier_1 file ----
 
  MPRA_Tier_0 = readRDS(opt$MPRA_Tier_0)
  
  # cat("MPRA_Tier_0_\n")
  # cat(str(MPRA_Tier_0))
  # cat("\n")
  
  
   
  MPRA_Tier_1 = readRDS(opt$MPRA_Tier_1)
  
  # cat("MPRA_Tier_1_\n")
  # cat(str(MPRA_Tier_1))
  # cat("\n")
  
  
  
  VARS_MPRA_ACTIVE<-unique(MPRA_Tier_1$VAR)
  
  # cat("VARS_MPRA_ACTIVE_\n")
  # cat(str(VARS_MPRA_ACTIVE))
  # cat("\n")
  
  VARS_MPRA_INACTIVE<-unique(MPRA_Tier_0$VAR)
  
  # cat("VARS_MPRA_INACTIVE_0\n")
  # cat(str(VARS_MPRA_INACTIVE))
  # cat("\n")
  
  # VARS_MPRA_INACTIVE<-VARS_MPRA_INACTIVE[-which(VARS_MPRA_INACTIVE%in%VARS_MPRA_ACTIVE)]
  # 
  # 
  # cat("VARS_MPRA_INACTIVE_1\n")
  # cat(str(VARS_MPRA_INACTIVE))
  # cat("\n")
  
  
  VARS_MPRA_screened<-c(VARS_MPRA_ACTIVE,VARS_MPRA_INACTIVE)
  
  # cat("VARS_MPRA_screened_1\n")
  # cat(str(VARS_MPRA_screened))
  # cat("\n")
  
  
  MPRA_Tier_0_sel<-MPRA_Tier_0[which(MPRA_Tier_0$VAR%in%VARS_MPRA_INACTIVE),]
  
  MPRA_Tier_0_sel$K562<-"INACTIVE"
  MPRA_Tier_0_sel$CHRF<-"INACTIVE"
  MPRA_Tier_0_sel$THP1<-"INACTIVE"
  MPRA_Tier_0_sel$HL60<-"INACTIVE"
  
  MPRA_Tier_0_sel$CLASSIFICATION<-"INACTIVE"
  
  
  # cat("MPRA_Tier_0_sel_\n")
  # cat(str(MPRA_Tier_0_sel))
  # cat("\n")
  
  MPRA_Tier_1$CLASSIFICATION<-"enhancer_and_ASE"
  
  MPRA_Tier_1$K562[which(MPRA_Tier_1$K562 == 0)]<-"INACTIVE"
  MPRA_Tier_1$K562[which(MPRA_Tier_1$K562 == 1)]<-"enhancer_and_ASE"
  MPRA_Tier_1$CHRF[which(MPRA_Tier_1$CHRF == 0)]<-"INACTIVE"
  MPRA_Tier_1$CHRF[which(MPRA_Tier_1$CHRF == 1)]<-"enhancer_and_ASE"
  MPRA_Tier_1$HL60[which(MPRA_Tier_1$HL60 == 0)]<-"INACTIVE"
  MPRA_Tier_1$HL60[which(MPRA_Tier_1$HL60 == 1)]<-"enhancer_and_ASE"
  MPRA_Tier_1$THP1[which(MPRA_Tier_1$THP1 == 0)]<-"INACTIVE"
  MPRA_Tier_1$THP1[which(MPRA_Tier_1$THP1 == 1)]<-"enhancer_and_ASE"
  
  # cat("MPRA_Tier_1_\n")
  # cat(str(MPRA_Tier_1))
  # cat("\n")
  # 
  
  MPRA_merge<-rbind(MPRA_Tier_0_sel,
                    MPRA_Tier_1)
  
  # cat("MPRA_merge_0\n")
  # cat(str(MPRA_merge))
  # cat("\n")
  
   
  #### Read genIE results ----
  
  
  
  genIE_RESULTS<-readRDS(opt$genIE_RESULTS)
  
  genIE_RESULTS$VAR<-as.character(genIE_RESULTS$VAR)
  
  
  # cat("genIE_RESULTS\n")
  # str(genIE_RESULTS)
  # cat("\n")
  
  
  genIE_ACTIVE<-genIE_RESULTS[which(genIE_RESULTS$DEF_CLASS == "SIGNIFICANT"),]
  
  # cat("genIE_ACTIVE\n")
  # str(genIE_ACTIVE)
  # cat("\n")
  
  genIE_INACTIVE<-genIE_RESULTS[which(genIE_RESULTS$DEF_CLASS == "NON_SIGNIFICANT"),]
  
  # cat("genIE_INACTIVE\n")
  # str(genIE_INACTIVE)
  # cat("\n")
  
  genIE_EDITION_DROPOUT<-genIE_RESULTS[which(genIE_RESULTS$DEF_CLASS == "EDITION_DROPOUT"),]
  
  # cat("genIE_EDITION_DROPOUT\n")
  # str(genIE_EDITION_DROPOUT)
  # cat("\n")
  
  
  VARS_genIE_ACTIVE<-unique(genIE_ACTIVE$VAR)#[which(genIE_screened$KEY%in%genIE_Tier_1$KEY)])
  
  # cat("VARS_genIE_ACTIVE_\n")
  # cat(str(VARS_genIE_ACTIVE))
  # cat("\n")
  
  VARS_genIE_INACTIVE<-unique(genIE_INACTIVE$VAR)#[-which(genIE_screened$KEY%in%genIE_Tier_1$KEY)])
  
  # cat("VARS_genIE_INACTIVE_\n")
  # cat(str(VARS_genIE_INACTIVE))
  # cat("\n")
  
  
  VARS_genIE_EDITION_DROPOUT<-unique(genIE_EDITION_DROPOUT$VAR)#[-which(genIE_screened$KEY%in%genIE_Tier_1$KEY)])
  
  # cat("VARS_genIE_EDITION_DROPOUT_\n")
  # cat(str(VARS_genIE_EDITION_DROPOUT))
  # cat("\n")
  
  
 
  
   # quit(status = 1)

  
  
  #### ALL_dB_FINAL_subset ####
  
  setwd(path_Patrick)
  ALL_dB_FINAL_subset<-as.data.frame(fread(file="ALL_db_Plus_most_severe_CSQ.tsv",sep="\t") , stringsAsFactors=F)
  
  cat("ALL_dB_FINAL_subset_1\n")
  cat(str(ALL_dB_FINAL_subset))
  cat("\n")
  
  ALL_dB_for_Table_Nicole<-unique(ALL_dB_FINAL_subset[,c(which(colnames(ALL_dB_FINAL_subset) == "VAR"),
                                                  which(colnames(ALL_dB_FINAL_subset) == "rs"))])
  
  cat("ALL_dB_for_Table_Nicole_1\n")
  cat(str(ALL_dB_for_Table_Nicole))
  cat("\n")
  
 
  
  
  
  MPRA_merge<-merge(MPRA_merge,
                    ALL_dB_for_Table_Nicole,
                    by="VAR",
                    all.x=T)
  
  
  cat("MPRA_merge_1\n")
  cat(str(MPRA_merge))
  cat("\n")
  
  
  
  indx.dep<-c(which(colnames(MPRA_merge) == "KEY"),which(colnames(MPRA_merge) == "VEP_DEF_LABELS"))
  
  MPRA_merge<-MPRA_merge[,-indx.dep]
  
  colnames(MPRA_merge)[which(colnames(MPRA_merge) == "rs")]<-"rsid"
  
  indx.reorder<-c(which(colnames(MPRA_merge) == "rsid"),which(colnames(MPRA_merge) == "VAR"),which(colnames(MPRA_merge) == "CLASSIFICATION"),
                  which(colnames(MPRA_merge) == "K562"),which(colnames(MPRA_merge) == "CHRF"),which(colnames(MPRA_merge) == "HL60"),which(colnames(MPRA_merge) == "THP1"))
  
  MPRA_merge<-MPRA_merge[,indx.reorder]
  
  
  cat("MPRA_merge_2\n")
  cat(str(MPRA_merge))
  cat("\n")
  
  cat("genIE_ACTIVE\n")
  str(genIE_ACTIVE)
  cat("\n")
  cat(sprintf(as.character(unique(genIE_ACTIVE$HGNC))))
  cat("\n")
  
  #### Screened_variants ----
  
  SCREENED_VARIANTS<-unique(c(VARS_genIE_screened,VARS_MPRA_screened))
  
 
  cat("SCREENED_VARIANTS_2\n")
  cat(str(SCREENED_VARIANTS))
  cat("\n")
  
  

  
  
 
  
  
 
  #### Filter 1 SELECT RARE VARIANTS ----
  
  CODING_LABELS <-c("LOF","MISS","SYN","UTR5","UTR3")
  
  

  ALL_dB_filter_IS<-ALL_dB_FINAL_subset
  
  
  cat("ALL_dB_filter_IS\n")
  cat(str(ALL_dB_filter_IS))
  cat("\n")
  cat(str(unique(ALL_dB_filter_IS$VAR)))
  cat("\n")
  
 
  
  
  
  ALL_dB_filter_maf<-ALL_dB_filter_IS[which(ALL_dB_filter_IS$maf_origin <= maf_Threshold),] 
  
  cat("ALL_dB_filter_maf_\n")
  str(ALL_dB_filter_maf)
  cat("\n")
  cat(str(unique(ALL_dB_filter_maf$VAR)))
  cat("\n") #1476
  
 
 
  common_maf<-ALL_dB_filter_IS[which(ALL_dB_filter_IS$maf_origin > maf_Threshold),]

  # common_maf<-ALL_dB_filter_IS[which(ALL_dB_filter_IS$maf_origin > 0.05),]


  cat("common_maf_\n")
  str(common_maf)
  cat("\n")
  cat(str(unique(common_maf$VAR)))
  cat("\n")#177414
  
  
  
  SCREENED_VARIANTS_RV<-SCREENED_VARIANTS[which(SCREENED_VARIANTS%in%ALL_dB_filter_maf$VAR)]
  
  cat("SCREENED_VARIANTS_RV_\n")
  str(SCREENED_VARIANTS_RV)
  cat("\n")
  
  SCREENED_VARIANTS_OUT<-SCREENED_VARIANTS[-which(SCREENED_VARIANTS%in%SCREENED_VARIANTS_RV)]
  
  cat("SCREENED_VARIANTS_OUT_\n")
  str(SCREENED_VARIANTS_OUT)
  cat("\n")
  cat(sprintf(as.character(SCREENED_VARIANTS_OUT)))
  cat("\n")
  
  
  check_VAR<-"chr3_71355240_G_C" # FOXP1-IT1 variant 0.0274
  
  check_VAR<-"chr4_1008212_C_T" # DGQK 0.0421
  
  check_VAR_df<-ALL_dB_filter_IS[which(ALL_dB_filter_IS$VAR == check_VAR),]
  
  cat("check_VAR_df_\n")
  str(check_VAR_df)
  cat("\n")
  
  cat("check_VAR_\n")
  str(check_VAR)
  cat("\n")
  
  check_out_1<-ALL_dB_filter_maf[which(ALL_dB_filter_maf$VAR%in%SCREENED_VARIANTS_OUT),]

  cat("------------------------------>check_out_1_\n")
  str(check_out_1)
  cat("\n")
  
  #### Filter 2 SELECT NON_coding RARE VARIANTS ----
  
  CODING_LABELS <-c("LOF","MISS","SYN","UTR5","UTR3")
   
   
  ALL_dB_filter_NC<-ALL_dB_filter_maf[-which(ALL_dB_filter_maf$VEP_DEF_LABELS_wCSQ%in%CODING_LABELS),]
  
  cat("ALL_dB_filter_NC_\n")
  str(ALL_dB_filter_NC)
  cat("\n")
  cat(str(unique(ALL_dB_filter_NC$VAR)))
  cat("\n") #1476
  
  # LOW_FREQ_maf<-ALL_dB_filter_IS[which(ALL_dB_filter_IS$maf_origin > maf_Threshold &  ALL_dB_filter_IS$maf_origin <= 0.05),]
  # 
  # cat("LOW_FREQ_maf_\n")
  # str(LOW_FREQ_maf)
  # cat("\n")
  # cat(str(unique(LOW_FREQ_maf$VAR)))
  # cat("\n")#177414
  # 
  # 
  # LOW_FREQ_maf_CODING<-LOW_FREQ_maf[which(LOW_FREQ_maf$VEP_DEF_LABELS_wCSQ%in%CODING_LABELS),]
  # 
  # cat(str("LOW_FREQ_maf_CODING\n"))
  # str(LOW_FREQ_maf_CODING)
  # cat("\n")
  # cat(str(unique(LOW_FREQ_maf_CODING$VAR)))
  # cat("\n")
  # 
  # 
  # LOW_FREQ_maf_NON_CODING<-LOW_FREQ_maf[-which(LOW_FREQ_maf$VEP_DEF_LABELS_wCSQ%in%CODING_LABELS),]
  # 
  # cat(str("LOW_FREQ_maf_NON_CODING\n"))
  # str(LOW_FREQ_maf_NON_CODING)
  # cat("\n")
  # cat(str(unique(LOW_FREQ_maf_NON_CODING$VAR)))
  # cat("\n")
  
  
  
  
  
  
  
  SCREENED_VARIANTS_NC<-SCREENED_VARIANTS[which(SCREENED_VARIANTS%in%ALL_dB_filter_NC$VAR)]
  
  cat("SCREENED_VARIANTS_NC_\n")
  str(SCREENED_VARIANTS_NC)
  cat("\n")
  
  SCREENED_VARIANTS_OUT<-SCREENED_VARIANTS[-c(which(SCREENED_VARIANTS%in%SCREENED_VARIANTS_NC),
                                              which(SCREENED_VARIANTS%in%SCREENED_VARIANTS_RV))]
  
  cat("SCREENED_VARIANTS_OUT_\n")
  str(SCREENED_VARIANTS_OUT)
  cat("\n")
  cat(sprintf(as.character(SCREENED_VARIANTS_OUT)))
  cat("\n")
  
  
  check_VAR<-"chr3_71355240_G_C" # FOXP1-IT1 variant 0.0274
  
  check_VAR<-"chr4_1008212_C_T" # DGQK 0.0421
  
  check_VAR_df<-ALL_dB_filter_IS[which(ALL_dB_filter_IS$VAR == check_VAR),]
  
  cat("check_VAR_df_\n")
  str(check_VAR_df)
  cat("\n")
  
  cat("check_VAR_\n")
  str(check_VAR)
  cat("\n")
  
  check_out_2<-ALL_dB_filter_NC[which(ALL_dB_filter_NC$VAR%in%SCREENED_VARIANTS_OUT),]
  
  cat("------------------------------>check_out_2_\n")
  str(check_out_2)
  cat("\n")
  
  # 
  #### Filter 3  PP >= 0.9 Effect sizes NON_coding RARE VARIANTS ----
  
  ALL_dB_filter_3<-ALL_dB_filter_NC[which(ALL_dB_filter_NC$finemap_prob >= Prob_Threshold),]
  
  
  
  cat("ALL_dB_filter_3_\n")
  str(ALL_dB_filter_3)
  cat("\n")
  cat(str(unique(ALL_dB_filter_3$VAR)))
  cat("\n") #238
  
  
  SCREENED_VARIANTS_PP<-SCREENED_VARIANTS[which(SCREENED_VARIANTS%in%ALL_dB_filter_3$VAR)]
  
  cat("SCREENED_VARIANTS_PP_\n")
  str(SCREENED_VARIANTS_PP)
  cat("\n")
  
  
  SCREENED_VARIANTS_OUT<-SCREENED_VARIANTS[-c(which(SCREENED_VARIANTS%in%SCREENED_VARIANTS_RV),
                                              which(SCREENED_VARIANTS%in%SCREENED_VARIANTS_NC),
                                              which(SCREENED_VARIANTS%in%SCREENED_VARIANTS_PP))]
  
  cat("SCREENED_VARIANTS_OUT_\n")
  str(SCREENED_VARIANTS_OUT)
  cat("\n")
  cat(sprintf(as.character(SCREENED_VARIANTS_OUT)))
  cat("\n")
  
  
  check_out_3<-ALL_dB_filter_3[which(ALL_dB_filter_3$VAR%in%SCREENED_VARIANTS_OUT),]
  
  cat("------------------------------>check_out_3_\n")
  str(check_out_3)
  cat("\n")
  
  
  
  #### Filter 4 Effect sizes NON_coding RARE VARIANTS ----
  
  ALL_dB_filter_IS_FINEMAP_BETA<-ALL_dB_filter_IS[!is.na(ALL_dB_filter_IS$finemap_beta),]
  
  # cat("ALL_dB_filter_IS_FINEMAP_BETA\n")
  # str(ALL_dB_filter_IS_FINEMAP_BETA)
  # cat("\n")
  
  ALL_dB_filter_IS_FINEMAP_BETA.dt<-data.table(ALL_dB_filter_IS_FINEMAP_BETA, key=c("VAR","rs","maf_origin","VEP_DEF_LABELS_wCSQ"))
  
  
  phenotype.parameters <- as.data.frame(ALL_dB_filter_IS_FINEMAP_BETA.dt[,.(max_FOR_TRAIT = summary(finemap_beta)[6]+0.1,
                                                  min_FOR_TRAIT = summary(finemap_beta)[1]-0.1,
                                                  Q1=summary(finemap_beta)[2],
                                                  Q3=summary(finemap_beta)[5]), by = phenotype])
  
  # cat("phenotype.parameters_\n")
  # str(phenotype.parameters)
  # cat("\n")
  
  Phenotypes_array<-levels(as.factor(phenotype.parameters$phenotype))
  
  # cat("Phenotypes_array_0\n")
  # cat(str(Phenotypes_array))
  # cat("\n")
  
  
  list_phenotypes<-list()
  
  for(iteration_Phenotypes_array in 1:length(Phenotypes_array))
  {
    Phenotypes_array_sel<-Phenotypes_array[iteration_Phenotypes_array]
    
    # cat("------------------>\t")
    # cat(sprintf(as.character(Phenotypes_array_sel)))
    # cat("\n")
    
    ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel<-ALL_dB_filter_IS_FINEMAP_BETA[which(ALL_dB_filter_IS_FINEMAP_BETA$phenotype == Phenotypes_array_sel),]
    
    # cat("ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel_0\n")
    # cat(str(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel))
    # cat("\n")
    
    phenotype.parameters_phenotype_sel<-phenotype.parameters[which(phenotype.parameters$phenotype == Phenotypes_array_sel),]
    
    # cat("phenotype.parameters_phenotype_sel\n")
    # cat(str(phenotype.parameters_phenotype_sel))
    # cat("\n")
    # 
    
    
    
    A_points<-c(phenotype.parameters_phenotype_sel$min_FOR_TRAIT,
                phenotype.parameters_phenotype_sel$Q1,
                phenotype.parameters_phenotype_sel$Q3,
                phenotype.parameters_phenotype_sel$max_FOR_TRAIT)
    
   
    
    # cat("A_points_0\n")
    # cat(sprintf(as.character(A_points)))
    # cat("\n")
    
    duplicated_A_points_index<-which(duplicated(A_points))
    
    # cat("duplicated_A_points_index\n")
    # cat(str(duplicated_A_points_index))
    # cat("\n")
    
   
    
    if(length(duplicated_A_points_index) >0)
    {
      for(iteration_duplicated_A_points_index in 1:length(duplicated_A_points_index))
      {
        if(iteration_duplicated_A_points_index == 1)
        {
          sum_term=0.1
        }else{
          sum_term=sum_term+0.1
          
        }
        
        duplicated_A_points_index_sel<-duplicated_A_points_index[iteration_duplicated_A_points_index]
        
        # cat("----------------------------->\t")
        # cat(sprintf(as.character(duplicated_A_points_index_sel)))
        # cat("\t")
        
        A_points_duplicated<-A_points[duplicated_A_points_index_sel]
        
        # cat("----------------------------->\t")
        # cat(sprintf(as.character(A_points_duplicated)))
        # cat("\n")
        
        A_points_duplicated_updated<-A_points[duplicated_A_points_index_sel]+sum_term
        
        # cat("--->\t")
        # cat(sprintf(as.character(A_points_duplicated_updated)))
        # cat("\n")
        
        A_points[duplicated_A_points_index_sel]<-A_points_duplicated_updated
      }
      
    }
    
    # cat("A_points_1\n")
    # cat(sprintf(as.character(A_points)))
    # cat("\n")
    
    
   
    phenotype.parameters_paleo_sel<-phenotype.parameters_paleo[which(phenotype.parameters_paleo$phenotype == Phenotypes_array_sel),]
    
    # cat("phenotype.parameters_paleo_sel\n")
    # cat(str(phenotype.parameters_paleo_sel))
    # cat("\n")
    
    if(dim(phenotype.parameters_paleo_sel)[1] >0)
    {
      A_points<-c(A_points[1],
                  phenotype.parameters_paleo_sel$Q1,
                  phenotype.parameters_paleo_sel$Q3,
                  A_points[4])
      
      # cat("A_points_paleo\n")
      # cat(sprintf(as.character(A_points)))
      # cat("\n")
      
      ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$classif_breaks<-cut(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$finemap_beta,breaks = A_points,right = FALSE)
      
      
      # cat("ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel_1\n")
      # cat(str(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel))
      # cat("\n")
      
      ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$DEFINITIVE_CLASS<-"NA"
      
      ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$DEFINITIVE_CLASS[(as.numeric(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$classif_breaks) == 1)]<-"MIN_Q1"
      ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$DEFINITIVE_CLASS[(as.numeric(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$classif_breaks) == 2)]<-"Q1_Q3"
      ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$DEFINITIVE_CLASS[(as.numeric(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$classif_breaks) == 3)]<-"Q3_MAX"
      
      
      # cat("ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel_2\n")
      # cat(str(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel))
      # cat("\n")
      
      
      
      
      indx.int<-c(which(colnames(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel) == "VAR"),
                  which(colnames(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel) == "phenotype"),
                  which(colnames(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel) == "DEFINITIVE_CLASS"))
      
      # quit(status = 1)
      
    }else{
      
      ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$classif_breaks<-cut(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$finemap_beta,breaks = A_points,right = FALSE)
      
      
      # cat("ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel_1\n")
      # cat(str(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel))
      # cat("\n")
      
      ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$DEFINITIVE_CLASS<-"NA"
      
      ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$DEFINITIVE_CLASS[(as.numeric(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$classif_breaks) == 1)]<-"MIN_Q1"
      ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$DEFINITIVE_CLASS[(as.numeric(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$classif_breaks) == 2)]<-"Q1_Q3"
      ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$DEFINITIVE_CLASS[(as.numeric(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$classif_breaks) == 3)]<-"Q3_MAX"
      
      
      # cat("ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel_2\n")
      # cat(str(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel))
      # cat("\n")
      
      
      
      
      indx.int<-c(which(colnames(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel) == "VAR"),
                  which(colnames(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel) == "phenotype"),
                  which(colnames(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel) == "DEFINITIVE_CLASS"))
    }
    
    
    
    ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$PCHIC_classif<-cut(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$finemap_beta,breaks = A_points,right = FALSE)
    
    
    # cat("ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel_1\n")
    # cat(str(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel))
    # cat("\n")
    
    ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$DEFINITIVE_CLASS<-"NA"
    
    ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$DEFINITIVE_CLASS[(as.numeric(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$PCHIC_classif) == 1)]<-"MIN_Q1"
    ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$DEFINITIVE_CLASS[(as.numeric(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$PCHIC_classif) == 2)]<-"Q1_Q3"
    ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$DEFINITIVE_CLASS[(as.numeric(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel$PCHIC_classif) == 3)]<-"Q3_MAX"
   
    
    # cat("ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel_2\n")
    # cat(str(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel))
    # cat("\n")
    
    
   
    
    indx.int<-c(which(colnames(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel) == "VAR"),
                which(colnames(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel) == "phenotype"),
                which(colnames(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel) == "DEFINITIVE_CLASS"))
    
    #which(colnames(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel) == "lineage"),
    
    subset<-unique(ALL_dB_filter_IS_FINEMAP_BETA_phenotype_sel[,indx.int])
    
    # cat("subset_2\n")
    # cat(str(subset))
    # cat("\n")
    
    list_phenotypes[[iteration_Phenotypes_array]]<-subset
    
    # #####################################
    # quit(status = 1)
    
  }# iteration_Phenotypes_array
  
  
  Phenotypes_df = as.data.frame(data.table::rbindlist(list_phenotypes, fill=T), stringsAsFactors=F)
  
  
  cat("Phenotypes_df\n")
  cat(str(Phenotypes_df))
  cat("\n")
  
  ALL_dB_filter_IS_FINEMAP_BETA<-merge(ALL_dB_filter_IS_FINEMAP_BETA,
                                       Phenotypes_df,
                                       by=c("phenotype","VAR"))
  
  cat("ALL_dB_filter_IS_FINEMAP_BETA_0\n")
  str(ALL_dB_filter_IS_FINEMAP_BETA)
  cat("\n")
  cat(str(unique(ALL_dB_filter_IS_FINEMAP_BETA$VAR)))
  cat("\n") # 178890
  cat(sprintf(as.character(names(summary(as.factor(ALL_dB_filter_IS_FINEMAP_BETA$DEFINITIVE_CLASS))))))
  cat("\n") 
  cat(sprintf(as.character(summary(as.factor(ALL_dB_filter_IS_FINEMAP_BETA$DEFINITIVE_CLASS)))))
  cat("\n")
  # MIN_Q1 Q1_Q3 Q3_MAX
  # 16202 300696 50330
  
  
  
  
  
  
  
  ALL_dB_filter_3_parameters<-merge(ALL_dB_filter_3, 
                                     Phenotypes_df,
                                    by=c("phenotype","VAR"))
  
  cat("ALL_dB_filter_3_parameters_0\n")
  str(ALL_dB_filter_3_parameters)
  cat("\n")
  cat(str(unique(ALL_dB_filter_3_parameters$VAR)))
  cat("\n") #196
  cat(sprintf(as.character(names(summary(as.factor(ALL_dB_filter_3_parameters$DEFINITIVE_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(ALL_dB_filter_3_parameters$DEFINITIVE_CLASS)))))
  cat("\n")
  
  # MIN_Q1 Q1_Q3 Q3_MAX
  # 130 173 102
  
  
  
  
  ALL_dB_filter_3_parameters.Q3<-ALL_dB_filter_3_parameters[which(ALL_dB_filter_3_parameters$DEFINITIVE_CLASS == "Q3_MAX"),]
  ALL_dB_filter_3_parameters.Q1<-ALL_dB_filter_3_parameters[which(ALL_dB_filter_3_parameters$DEFINITIVE_CLASS == "MIN_Q1"),]
  
  
  ALL_dB_filter_4<-rbind(ALL_dB_filter_3_parameters.Q1,ALL_dB_filter_3_parameters.Q3)
  
  cat("ALL_dB_filter_4_0\n")
  str(ALL_dB_filter_4)
  cat("\n")
  cat(str(unique(ALL_dB_filter_4$VAR)))
  cat("\n") #123
  cat(sprintf(as.character(names(summary(as.factor(ALL_dB_filter_4$DEFINITIVE_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(ALL_dB_filter_4$DEFINITIVE_CLASS)))))
  cat("\n")
  
  # MIN_Q1 Q3_MAX
  # 130 102
  
  
 
  
  SCREENED_VARIANTS_Effect_size<-SCREENED_VARIANTS[which(SCREENED_VARIANTS%in%ALL_dB_filter_4$VAR)]
  
  cat("SCREENED_VARIANTS_Effect_size_\n")
  str(SCREENED_VARIANTS_Effect_size)
  cat("\n")#91
  
  SCREENED_VARIANTS_OUT<-SCREENED_VARIANTS[-c(which(SCREENED_VARIANTS%in%ALL_dB_filter_4$VAR))]
  
  cat("SCREENED_VARIANTS_OUT_\n")
  str(SCREENED_VARIANTS_OUT)
  cat("\n")
  cat(sprintf(as.character(SCREENED_VARIANTS_OUT)))
  cat("\n")
  
  
  
  
  check_VAR<-"chr12_111844956_C_T" # SH2B3

  check_VAR_df<-ALL_dB_filter_3_parameters[which(ALL_dB_filter_3_parameters$VAR == check_VAR),]

  cat("check_VAR_df_\n")
  str(check_VAR_df)
  cat("\n")

  phenotype.parameters_paleo_sel<-phenotype.parameters_paleo[which(phenotype.parameters_paleo$phenotype == check_VAR_df$phenotype),]
  
  cat("phenotype.parameters_paleo_sel_\n")
  str(phenotype.parameters_paleo_sel)
  cat("\n")
  
 
  
  
  #### Fig0_Annot_Category  ----
  
  ALL_dB_CSQ_IS_filtered<-unique(ALL_dB_filter_IS[,c(which(colnames(ALL_dB_filter_IS) == "VAR"),which(colnames(ALL_dB_filter_IS) == "rs"),
                                                             which(colnames(ALL_dB_filter_IS) == "maf_origin"),which(colnames(ALL_dB_filter_IS) == "VEP_DEF_LABELS_wCSQ"))])
  
  cat("ALL_dB_CSQ_IS_filtered_0:\n")
  cat(str(ALL_dB_CSQ_IS_filtered))
  cat("\n")
  
  
  
  

  
  ALL_dB_CSQ_IS_filtered_CODING<-ALL_dB_CSQ_IS_filtered[which(ALL_dB_CSQ_IS_filtered$VEP_DEF_LABELS_wCSQ%in%CODING_LABELS),]
  
  cat(str("ALL_dB_CSQ_IS_filtered_CODING\n"))
  str(ALL_dB_CSQ_IS_filtered_CODING)
  cat("\n")
  cat(str(unique(ALL_dB_CSQ_IS_filtered_CODING$VAR)))
  cat("\n")
  
  
  ALL_dB_CSQ_IS_filtered_NON_CODING<-ALL_dB_CSQ_IS_filtered[-which(ALL_dB_CSQ_IS_filtered$VEP_DEF_LABELS_wCSQ%in%CODING_LABELS),]
  
  cat(str("ALL_dB_CSQ_IS_filtered_NON_CODING\n"))
  str(ALL_dB_CSQ_IS_filtered_NON_CODING)
  cat("\n")
  cat(str(unique(ALL_dB_CSQ_IS_filtered_NON_CODING$VAR)))
  cat("\n")
  
  
  
  
  ALL_dB_CSQ_IS_filtered$Fig0_Annot_Category<-NA
  
  ALL_dB_CSQ_IS_filtered$Fig0_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%ALL_dB_CSQ_IS_filtered_CODING$VAR &
                                                     ALL_dB_CSQ_IS_filtered$VAR%in%ALL_dB_filter_maf$VAR)]<-"RV_CODING"
 
  ALL_dB_CSQ_IS_filtered$Fig0_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%ALL_dB_CSQ_IS_filtered_NON_CODING$VAR &
                                                     ALL_dB_CSQ_IS_filtered$VAR%in%ALL_dB_filter_maf$VAR)]<-"RV_NC"
  
  
  ALL_dB_CSQ_IS_filtered$Fig0_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%common_maf$VAR)]<-"Common_variant"
  
  

  
  ### LEVELS ORDERED
  
  
  ALL_dB_CSQ_IS_filtered$Fig0_Annot_Category<-factor(ALL_dB_CSQ_IS_filtered$Fig0_Annot_Category,
                                                  levels=c("Common_variant","RV_CODING","RV_NC"),
                                                  ordered=T)
  
  
  cat("-------------------------------------------->Fig0_Annot_Category\n")
  cat(sprintf(as.character(names(summary(ALL_dB_CSQ_IS_filtered$Fig0_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB_CSQ_IS_filtered$Fig0_Annot_Category))))
  cat("\n")
  
 
  
  #### Fig1_Annot_Category  ----
  
  ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category<-"Common_variant"
  
  ### RV 

  
  ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%ALL_dB_filter_maf$VAR)]<-"RV"
  
  
  cat("Classification_RV_vs_common\n")
  str(ALL_dB_CSQ_IS_filtered)
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category))))
  cat("\n")
  
  ## NC
  
  ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%ALL_dB_filter_maf$VAR)]<-"RV_C"
  
  
  ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%ALL_dB_filter_NC$VAR)]<-"RV_NC"
  
  cat("Classification_RV_C_vs_RV_NC\n")
  str(ALL_dB_CSQ_IS_filtered)
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category))))
  cat("\n")
  
  
  
  ## PP 
  
  ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%ALL_dB_filter_NC$VAR)]<-"RV_NC_lowPP"
  
  
  ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%ALL_dB_filter_3$VAR)]<-"RV_NC_highPP"
  
  cat("Classification_RV_NC_lowPPhighEffect_vs_RV_NC_highPPlowEffect\n")
  str(ALL_dB_CSQ_IS_filtered)
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category))))
  cat("\n")
  
  ## Effect_size
  
  ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%ALL_dB_filter_3$VAR)]<-"RV_NC_highPP_lowEffectSize"
  
  ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%ALL_dB_filter_4$VAR)]<-"index_variants"
  
  cat("Classification_RV_NC_lowPP_vs_RV_NC_highPP\n")
  str(ALL_dB_CSQ_IS_filtered)
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category))))
  cat("\n")
  
  ### Order levels
  
  ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category<-factor(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category,
                                                     levels=c("Common_variant","RV_C","RV_NC_lowPP","RV_NC_highPP_lowEffectSize","index_variants"),
                                                     ordered=T)
  
  cat("-------------------------------------------->Fig1_Annot_Category\n")
  cat("ALL_dB_CSQ_IS_filtered_MIDDLE_PASS\n")
  str(ALL_dB_CSQ_IS_filtered)
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category))))
  cat("\n")
  

  
  check<-sum(summary(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category))
  
  cat("check_MIDDLE_PASS\n")
  str(check)
  cat("\n")
  
  check2<-droplevels(ALL_dB_CSQ_IS_filtered[which(ALL_dB_CSQ_IS_filtered$VAR%in%SCREENED_VARIANTS),])
  
  cat("check2\n")
  str(check2)
  cat("\n")
  cat(sprintf(as.character(names(summary(check2$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(check2$Fig1_Annot_Category))))
  cat("\n")
  
  
  ### Fig2_Annot_Category ----
  
 # SCREENED_VARIANTS<-unique(c(VARS_genIE_screened,VARS_MPRA_screened))
  
  VARS_MPRA_screened_corrected<-VARS_MPRA_screened[which(VARS_MPRA_screened%in%ALL_dB_filter_4$VAR)]
  
  VARS_MPRA_INACTIVE_corrected<-VARS_MPRA_INACTIVE[which(VARS_MPRA_INACTIVE%in%VARS_MPRA_screened_corrected)]
  VARS_MPRA_ACTIVE_corrected<-VARS_MPRA_ACTIVE[which(VARS_MPRA_ACTIVE%in%VARS_MPRA_screened_corrected)]
  

  ALL_dB_CSQ_IS_filtered$Fig2_Annot_Category<-"NA"
  
  
  ALL_dB_CSQ_IS_filtered$Fig2_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%VARS_MPRA_INACTIVE_corrected)]<-"NON_ACTIVE"
  ALL_dB_CSQ_IS_filtered$Fig2_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%VARS_MPRA_ACTIVE_corrected)]<-"ACTIVE"
  
  indx.prioritised<-which(as.numeric(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category) == 5)
  
  cat("indx.prioritised\n")
  str(indx.prioritised)
  cat("\n")
  
  ALL_dB_CSQ_IS_filtered$Fig2_Annot_Category[-indx.prioritised]<-"NA"
  
  
  #### SWITCH ON/OFF LINE 
  
  ALL_dB_CSQ_IS_filtered$Fig2_Annot_Category[which(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category == "Common_variant")]<-"NA"
    

  ALL_dB_CSQ_IS_filtered$Fig2_Annot_Category<-factor(ALL_dB_CSQ_IS_filtered$Fig2_Annot_Category,
                                                     levels=c("NON_ACTIVE","ACTIVE"),
                                                     ordered=T)
  
  cat("-------------------------------------------->Fig2_Annot_Category\n")
  cat(sprintf(as.character(names(summary(ALL_dB_CSQ_IS_filtered$Fig2_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB_CSQ_IS_filtered$Fig2_Annot_Category))))
  cat("\n")
  
  # quit(status = 1)
  
  ### Fig3_Annot_Category----
  
  VARS_genIE_screened_corrected<-VARS_genIE_screened[which(VARS_genIE_screened%in%ALL_dB_filter_3_parameters$VAR)]
  
  VARS_genIE_INACTIVE_corrected<-VARS_genIE_INACTIVE[which(VARS_genIE_INACTIVE%in%VARS_genIE_screened_corrected)]
  VARS_genIE_ACTIVE_corrected<-VARS_genIE_ACTIVE[which(VARS_genIE_ACTIVE%in%VARS_genIE_screened_corrected)]
  VARS_genIE_EDITION_DROPOUT_corrected<-VARS_genIE_EDITION_DROPOUT[which(VARS_genIE_EDITION_DROPOUT%in%VARS_genIE_screened_corrected)]
  
  
  
  ALL_dB_CSQ_IS_filtered$Fig3_Annot_Category<-"NA"
  
  ALL_dB_CSQ_IS_filtered$Fig3_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%VARS_genIE_INACTIVE)]<-"NON_ACTIVE"
  ALL_dB_CSQ_IS_filtered$Fig3_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%VARS_genIE_ACTIVE)]<-"ACTIVE"
  ALL_dB_CSQ_IS_filtered$Fig3_Annot_Category[which(ALL_dB_CSQ_IS_filtered$VAR%in%VARS_genIE_EDITION_DROPOUT)]<-"EDITION_DROPOUT"
  
  indx.prioritised<-which(as.numeric(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category) == 5)
  
  cat("indx.prioritised\n")
  str(indx.prioritised)
  cat("\n")
  
  ALL_dB_CSQ_IS_filtered$Fig3_Annot_Category[-indx.prioritised]<-"NA"
  
  #### genIE Absent----
  
  genIE_RESULTS_ABSENT<-genIE_RESULTS[-which(genIE_RESULTS$VAR%in%ALL_dB_CSQ_IS_filtered$VAR),]
  
  cat("genIE_RESULTS_ABSENT\n")
  cat(str(genIE_RESULTS_ABSENT))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(genIE_RESULTS_ABSENT$HGNC))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(genIE_RESULTS_ABSENT$HGNC)))))
  cat("\n")
  
  
  
  #### SWITCH ON/OFF LINE

  # ALL_dB_CSQ_IS_filtered$Fig3_Annot_Category[which(ALL_dB_CSQ_IS_filtered$Fig1_Annot_Category == "Common_variant")]<-"NA"
  
  
  ALL_dB_CSQ_IS_filtered$Fig3_Annot_Category<-factor(ALL_dB_CSQ_IS_filtered$Fig3_Annot_Category,
                                                     levels=c("EDITION_DROPOUT","NON_ACTIVE","ACTIVE"),
                                                     ordered=T)
  
  cat("-------------------------------------------->Fig3_Annot_Category\n")
  cat(sprintf(as.character(names(summary(ALL_dB_CSQ_IS_filtered$Fig3_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB_CSQ_IS_filtered$Fig3_Annot_Category))))
  cat("\n")
  

  
  
  #### GET RID OF VEP_DEF_LABELS_wCSQ ----
  
  cat("ALL_dB_CSQ_IS_filtered\n")
  cat(str(ALL_dB_CSQ_IS_filtered))
  cat("\n")
  
  indx.dep<-which(colnames(ALL_dB_CSQ_IS_filtered) == "VEP_DEF_LABELS_wCSQ")
  
  cat("indx.dep\n")
  cat(str(indx.dep))
  cat("\n")
  
  Labelling_df_FINAL<-unique(ALL_dB_CSQ_IS_filtered[,-c(indx.dep)])
  
  cat("Labelling_df_FINAL\n")
  cat(str(Labelling_df_FINAL))
  cat("\n")
  
  
 
  
 
  
    #### Categories colors ----
  
  ###
  
  Fig0_Annot_Category_levels<-levels(Labelling_df_FINAL$Fig0_Annot_Category)
  colors_Fig0_Annot_Category_levels<-c('gray','#CC3333','#3399FF','#D45E85','gray','gray','gray','gray')
  
  
  df.Fig0_Annot_Category<-as.data.frame(cbind(Fig0_Annot_Category_levels,colors_Fig0_Annot_Category_levels[1:length(Fig0_Annot_Category_levels)]), stringsAsFactors=F)
  colnames(df.Fig0_Annot_Category)<-c("Category","colors")
  
  df.Fig0_Annot_Category$Annotation<-"Fig0_Annot_Category"
  
  
  cat("df.Fig0_Annot_Category_0\n")
  cat(str(df.Fig0_Annot_Category))
  cat("\n")
  
  ###
  
  Fig1_Annot_Category_levels<-levels(Labelling_df_FINAL$Fig1_Annot_Category)
  colors_Fig1_Annot_Category_levels<-c('#32A852','#1877C9','#553B68','#62D07F','#C9244B', '#D45E85','#6DB2EE')
  
  
  df.Fig1_Annot_Category<-as.data.frame(cbind(Fig1_Annot_Category_levels,colors_Fig1_Annot_Category_levels[1:length(Fig1_Annot_Category_levels)]), stringsAsFactors=F)
  colnames(df.Fig1_Annot_Category)<-c("Category","colors")
  
  df.Fig1_Annot_Category$Annotation<-"Fig1_Annot_Category"
  
  cat("df.Fig1_Annot_Category_0\n")
  cat(str(df.Fig1_Annot_Category))
  cat("\n")
  
  
  ###
  
  
  Fig2_Annot_Category_levels<-levels(Labelling_df_FINAL$Fig2_Annot_Category)
  colors_Fig2_Annot_Category_levels<-c('red','greenyellow','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6')
  
  
  df.Fig2_Annot_Category<-as.data.frame(cbind(Fig2_Annot_Category_levels,colors_Fig2_Annot_Category_levels[1:length(Fig2_Annot_Category_levels)]), stringsAsFactors=F)
  colnames(df.Fig2_Annot_Category)<-c("Category","colors")
  
  df.Fig2_Annot_Category$Annotation<-"Fig2_Annot_Category"
  
  
  cat("df.Fig2_Annot_Category_0\n")
  cat(str(df.Fig2_Annot_Category))
  cat("\n")
  
  ###
  
  
  Fig3_Annot_Category_levels<-levels(Labelling_df_FINAL$Fig3_Annot_Category)
  colors_Fig3_Annot_Category_levels<-c('gray','red','greenyellow','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6')
  
  
  df.Fig3_Annot_Category<-as.data.frame(cbind(Fig3_Annot_Category_levels,colors_Fig3_Annot_Category_levels[1:length(Fig3_Annot_Category_levels)]), stringsAsFactors=F)
  colnames(df.Fig3_Annot_Category)<-c("Category","colors")
  
  df.Fig3_Annot_Category$Annotation<-"Fig3_Annot_Category"
  
  cat("df.Fig3_Annot_Category_0\n")
  cat(str(df.Fig3_Annot_Category))
  cat("\n")
  
  ###
  
  
  
  DEF_colors<-rbind(df.Fig0_Annot_Category,df.Fig1_Annot_Category,df.Fig2_Annot_Category,df.Fig3_Annot_Category)
  
  cat("DEF_colors_0\n")
  cat(str(DEF_colors))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(DEF_colors$Annotation))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(DEF_colors$Annotation)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(DEF_colors$Category))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(DEF_colors$Category)))))
  cat("\n")
  
  
  #### SAVE ----
  
  cat("---------------------------------------------->Labelling_df_FINAL\n")
  str(Labelling_df_FINAL)
  cat("\n")
  cat(sprintf(as.character(names(summary(Labelling_df_FINAL$Fig0_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Labelling_df_FINAL$Fig0_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Labelling_df_FINAL$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Labelling_df_FINAL$Fig1_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Labelling_df_FINAL$Fig2_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Labelling_df_FINAL$Fig2_Annot_Category))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Labelling_df_FINAL$Fig3_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(Labelling_df_FINAL$Fig3_Annot_Category))))
  cat("\n")
  
  
  Table_S1<-Labelling_df_FINAL[,c(which(colnames(Labelling_df_FINAL) == "VAR"),which(colnames(Labelling_df_FINAL) == "rs"),
                                  which(colnames(Labelling_df_FINAL) == "maf_origin"),which(colnames(Labelling_df_FINAL) == "Fig1_Annot_Category"))]
  
  
  colnames(Table_S1)[which(colnames(Table_S1) == "Fig1_Annot_Category")]<-"Prioritisation_bins"
  
  setwd(out)
  
  filename2<-paste("Table_S1",".tsv",sep='')
  write.table(Table_S1,file=filename2,sep="\t",row.names = F,quote = F)
  
  filename2<-paste(type,"_Initial_Selection",".txt",sep='')
  write.table(Labelling_df_FINAL,file=filename2,sep="\t",row.names = F,quote = F)
  
  # quit(status = 1)
  
  filename2<-paste(type,"_Initial_Selection",".rds",sep='')
  saveRDS(Labelling_df_FINAL,file=filename2)
  
  filename2<-paste(type,"_Categories_colors",".rds",sep='')
  saveRDS(DEF_colors,file=filename2)
  
  ALL_dB_CSQ_IS_filtered$VEP_DEF_LABELS_wCSQ<-factor(ALL_dB_CSQ_IS_filtered$VEP_DEF_LABELS_wCSQ,
                                                          levels=c("LOF","MISS","SYN","UTR5","UTR3",
                                                                   "INTRON","INTERGENIC","UPSTREAM","DOWNSTREAM","REGULATORY",
                                                                   "TFBS","SPLICE","OTHER","NMD","NCT"),ordered = T)
  
  filename2<-paste(type,"_SCREENED_VARIANTS",".rds",sep='')
  saveRDS(SCREENED_VARIANTS,file=filename2)
  
  filename2<-paste(type,"_Initial_Selection_with_CSQ_labels",".rds",sep='')
  saveRDS(ALL_dB_CSQ_IS_filtered,file=filename2)
  
  filename2<-paste("Table_MPRA",".tsv",sep='')
  write.table(MPRA_merge,file=filename2,sep="\t",row.names = F,quote = F)
  
  filename2<-paste("Table_genIE",".tsv",sep='')
  write.table(genIE_ACTIVE,file=filename2,sep="\t",row.names = F,quote = F)
  
  
  ########## MAF parameters ----
  
  ALL_dB_per_category_values<-Labelling_df_FINAL
  
  cat("ALL_dB_per_category_values\n")
  str(ALL_dB_per_category_values)
  cat("\n")
  
  ALL_dB_per_category_values.dt<-data.table(ALL_dB_per_category_values, key=c("Fig0_Annot_Category"))
  
  
  MAF_parameters_Fig0 <- as.data.frame(ALL_dB_per_category_values.dt[,.(min_MAF = as.numeric(summary(maf_origin)[1]),
                                                                        Q1_MAF=as.numeric(summary(maf_origin)[2]),
                                                                        mean_MAF=as.numeric(summary(maf_origin)[3]),
                                                                        median_MAF=as.numeric(summary(maf_origin)[4]),
                                                                        Q3_MAF=as.numeric(summary(maf_origin)[5]),
                                                                        max_MAF = as.numeric(summary(maf_origin)[6])
  ),by=key(ALL_dB_per_category_values.dt)])
  
  MAF_parameters_Fig0$group<-"Fig0_Annot_Category"
  colnames(MAF_parameters_Fig0)[which(colnames(MAF_parameters_Fig0) == "Fig0_Annot_Category")]<-"levels"
  MAF_parameters_Fig0$levels<-as.character(MAF_parameters_Fig0$levels)
  
  cat("MAF_parameters_Fig0_\n")
  str(MAF_parameters_Fig0)
  cat("\n")
  
  ALL_dB_per_category_values.dt<-data.table(ALL_dB_per_category_values, key=c("Fig1_Annot_Category"))
  
  
  MAF_parameters_Fig1 <- as.data.frame(ALL_dB_per_category_values.dt[,.(min_MAF = as.numeric(summary(maf_origin)[1]),
                                                                         Q1_MAF=as.numeric(summary(maf_origin)[2]),
                                                                         mean_MAF=as.numeric(summary(maf_origin)[3]),
                                                                         median_MAF=as.numeric(summary(maf_origin)[4]),
                                                                         Q3_MAF=as.numeric(summary(maf_origin)[5]),
                                                                         max_MAF = as.numeric(summary(maf_origin)[6])
                                                                         ),by=key(ALL_dB_per_category_values.dt)])
  MAF_parameters_Fig1$group<-"Fig1_Annot_Category"
  colnames(MAF_parameters_Fig1)[which(colnames(MAF_parameters_Fig1) == "Fig1_Annot_Category")]<-"levels"
  MAF_parameters_Fig1$levels<-as.character(MAF_parameters_Fig1$levels)
  
  cat("MAF_parameters_Fig1_\n")
  str(MAF_parameters_Fig1)
  cat("\n")
  
  
  check_MAF_parameters_Fig1<-ALL_dB_per_category_values[which(ALL_dB_per_category_values$Fig1_Annot_Category == "Common_variant" &
                                                                ALL_dB_per_category_values$maf_origin < maf_Threshold  ),]
  cat("check_MAF_parameters_Fig1_\n")
  str(check_MAF_parameters_Fig1)
  cat("\n")
  
  
  
  ALL_dB_per_category_values.dt<-data.table(ALL_dB_per_category_values, key=c("Fig2_Annot_Category"))
  
  
  MAF_parameters_Fig2 <- as.data.frame(ALL_dB_per_category_values.dt[,.(min_MAF = as.numeric(summary(maf_origin)[1]),
                                                                        Q1_MAF=as.numeric(summary(maf_origin)[2]),
                                                                        mean_MAF=as.numeric(summary(maf_origin)[3]),
                                                                        median_MAF=as.numeric(summary(maf_origin)[4]),
                                                                        Q3_MAF=as.numeric(summary(maf_origin)[5]),
                                                                        max_MAF = as.numeric(summary(maf_origin)[6])
  ),by=key(ALL_dB_per_category_values.dt)])
  
  MAF_parameters_Fig2$group<-"Fig2_Annot_Category"
  colnames(MAF_parameters_Fig2)[which(colnames(MAF_parameters_Fig2) == "Fig2_Annot_Category")]<-"levels"
  MAF_parameters_Fig2$levels<-as.character(MAF_parameters_Fig2$levels)
  
  cat("MAF_parameters_Fig2_\n")
  str(MAF_parameters_Fig2)
  cat("\n")
  
  ALL_dB_per_category_values.dt<-data.table(ALL_dB_per_category_values, key=c("Fig3_Annot_Category"))
  
  
  MAF_parameters_Fig3 <- as.data.frame(ALL_dB_per_category_values.dt[,.(min_MAF = as.numeric(summary(maf_origin)[1]),
                                                                        Q1_MAF=as.numeric(summary(maf_origin)[2]),
                                                                        mean_MAF=as.numeric(summary(maf_origin)[3]),
                                                                        median_MAF=as.numeric(summary(maf_origin)[4]),
                                                                        Q3_MAF=as.numeric(summary(maf_origin)[5]),
                                                                        max_MAF = as.numeric(summary(maf_origin)[6])
  ),by=key(ALL_dB_per_category_values.dt)])
  
  MAF_parameters_Fig3$group<-"Fig3_Annot_Category"
  colnames(MAF_parameters_Fig3)[which(colnames(MAF_parameters_Fig3) == "Fig3_Annot_Category")]<-"levels"
  MAF_parameters_Fig3$levels<-as.character(MAF_parameters_Fig3$levels)
  
  cat("MAF_parameters_Fig3_\n")
  str(MAF_parameters_Fig3)
  cat("\n")
  
  
  DEF<-rbind(MAF_parameters_Fig0,MAF_parameters_Fig1,
             MAF_parameters_Fig2,MAF_parameters_Fig3)
  
  
  cat("DEF_\n")
  str(DEF)
  cat("\n")
  
  filename2<-paste(type,"_MAF_parameters",".txt",sep='')
  write.table(DEF,file=filename2,sep="\t",row.names = F,quote = F)
  
  # ####################################
  # quit(status = 1)

  
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
    make_option(c("--ALL_dB_Hannes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--paleo_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MPRA_Tier_0"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MPRA_Tier_1"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--genIE_RESULTS"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--genIE_Tier_1"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--genIE_input"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--maf_Threshold"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Info_Score_Threshold"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--finemap_prob_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Prob_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Initial_Selection"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
       make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--path_Patrick"), type="character", default=NULL, 
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
  
  Initial_Selection(opt)
 

}


###########################################################################

system.time( main() )
