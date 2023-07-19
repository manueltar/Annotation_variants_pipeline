
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
suppressMessages(library("cowplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))



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
  
  
  #### READ and transform subset_variables_1 ----
  
  subset_variables_1 = unlist(strsplit(opt$subset_variables_1, split=","))
  
  cat("subset_variables_1_\n")
  cat(sprintf(as.character(subset_variables_1)))
  cat("\n")
  
  df_subset_variables_1<-as.data.frame(cbind(rep("Variant_based_features", length(subset_variables_1)),subset_variables_1), stringsAsFactors=F)
  colnames(df_subset_variables_1)<-c("GROUP","variable")
  
  cat("df_subset_variables_1_0\n")
  cat(str(df_subset_variables_1))
  cat("\n")
  
  #### READ and transform subset_variables_2 ----
  
  subset_variables_2 = unlist(strsplit(opt$subset_variables_2, split=","))
  
  cat("subset_variables_2_\n")
  cat(sprintf(as.character(subset_variables_2)))
  cat("\n")
  
  df_subset_variables_2<-as.data.frame(cbind(rep("Gene_based_features", length(subset_variables_2)),subset_variables_2), stringsAsFactors=F)
  colnames(df_subset_variables_2)<-c("GROUP","variable")
  
  cat("df_subset_variables_2_0\n")
  cat(str(df_subset_variables_2))
  cat("\n")
  
  #### Merge ----
  
  df_subset_variables<-rbind(df_subset_variables_1,
                             df_subset_variables_2)
  
  cat("df_subset_variables_0\n")
  cat(str(df_subset_variables))
  cat("\n")
  cat(str(unique(df_subset_variables$GROUP)))
  cat("\n")
  cat(str(unique(df_subset_variables$variable)))
  cat("\n")
  
  #### Read Table_of_labels ----
  
  Table_of_labels<-droplevels(as.data.frame(readRDS(file=opt$Table_of_labels) , stringsAsFactors=F))
  
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
  cat(sprintf(as.character(names(summary(Table_of_labels$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_of_labels$MPRA_CLASS))))
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
  
  
  levels_MPRA_CLASS_selected<-levels(Table_of_labels$MPRA_CLASS)[c(1,2)]
  
  cat("levels_MPRA_CLASS_selected\n")
  cat(sprintf(as.character(levels_MPRA_CLASS_selected)))
  cat("\n")
  
  levels_M_and_M_selected<-levels(Table_of_labels$M_and_M)[c(1,13)]
  
  cat("levels_M_and_M_selected\n")
  cat(sprintf(as.character(levels_M_and_M_selected)))
  cat("\n")
  
  levels_Multi_Lineage_selected<-levels(Table_of_labels$Multi_Lineage)[c(1,2)]
  
  cat("levels_Multi_Lineage_selected\n")
  cat(sprintf(as.character(levels_Multi_Lineage_selected)))
  cat("\n")
  
  
  #### Categories colors ----
  
  DEF_colors<-readRDS(file=opt$Categories_colors)
  
  cat("DEF_colors_0\n")
  cat(str(DEF_colors))
  cat("\n")
  cat(sprintf(as.character(unique(DEF_colors$Annotation))))
  cat("\n")
  
  DEF_colors_sel_1<-DEF_colors[which(DEF_colors$Annotation == "Fig1_Annot_Category"),]
  
  cat("DEF_colors_sel_1_0\n")
  cat(str(DEF_colors_sel_1))
  cat("\n")
  cat(sprintf(as.character(unique(DEF_colors_sel_1$Annotation))))
  cat("\n")
  
  template_color<-DEF_colors[which(DEF_colors$Annotation == "Fig2_Annot_Category"),]
  
  cat("DEF_colors_sel_2_0\n")
  cat(str(template_color))
  cat("\n")
  cat(sprintf(as.character(unique(template_color$Annotation))))
  cat("\n")
  cat(sprintf(as.character(unique(template_color$colors))))
  cat("\n")
  
  DEF_colors_sel_2<-template_color
  
  DEF_colors_sel_2$Annotation[which(DEF_colors_sel_2$Annotation == "Fig2_Annot_Category")]<-"MPRA_CLASS"
  DEF_colors_sel_2$Category[which(DEF_colors_sel_2$Category == "NON_ACTIVE")]<-levels_MPRA_CLASS_selected[2]
  DEF_colors_sel_2$Category[which(DEF_colors_sel_2$Category == "ACTIVE")]<-levels_MPRA_CLASS_selected[1]
  
  DEF_colors_sel_2$colors[which(DEF_colors_sel_2$colors == "red")]<-"gray"
  
  cat("DEF_colors_sel_2_1\n")
  cat(str(DEF_colors_sel_2))
  cat("\n")
  cat(sprintf(as.character(unique(DEF_colors_sel_2$Annotation))))
  cat("\n")
  cat(sprintf(as.character(unique(DEF_colors_sel_2$colors))))
  cat("\n")
  
  
  
  DEF_colors_sel_3<-template_color
  
  DEF_colors_sel_3$Annotation[which(DEF_colors_sel_3$Annotation == "Fig2_Annot_Category")]<-"M_and_M"
  DEF_colors_sel_3$Category[which(DEF_colors_sel_3$Category == "NON_ACTIVE")]<-levels_M_and_M_selected[1]
  DEF_colors_sel_3$Category[which(DEF_colors_sel_3$Category == "ACTIVE")]<-levels_M_and_M_selected[2]
  
  DEF_colors_sel_3$colors[which(DEF_colors_sel_3$colors == "red")]<-"aquamarine4"
  DEF_colors_sel_3$colors[which(DEF_colors_sel_3$colors == "greenyellow")]<-"deeppink2"
  
  
  
  
  cat("DEF_colors_sel_3_1\n")
  cat(str(DEF_colors_sel_3))
  cat("\n")
  cat(sprintf(as.character(unique(DEF_colors_sel_3$Annotation))))
  cat("\n")
  cat(sprintf(as.character(unique(DEF_colors_sel_3$colors))))
  cat("\n")
  
  DEF_colors_sel_4<-template_color
  
  DEF_colors_sel_4$Annotation[which(DEF_colors_sel_4$Annotation == "Fig2_Annot_Category")]<-"Multi_Lineage"
  DEF_colors_sel_4$Category[which(DEF_colors_sel_4$Category == "NON_ACTIVE")]<-levels_Multi_Lineage_selected[2]
  DEF_colors_sel_4$Category[which(DEF_colors_sel_4$Category == "ACTIVE")]<-levels_Multi_Lineage_selected[1]
  
  # DEF_colors_sel_4$colors[which(DEF_colors_sel_4$colors == "red")]<-"gray"
  
  
  cat("DEF_colors_sel_4_1\n")
  cat(str(DEF_colors_sel_4))
  cat("\n")
  cat(sprintf(as.character(unique(DEF_colors_sel_4$Annotation))))
  cat("\n")
  cat(sprintf(as.character(unique(DEF_colors_sel_4$colors))))
  cat("\n")
  
  
  COLORS_DEF<-rbind(DEF_colors_sel_1,DEF_colors_sel_2,DEF_colors_sel_3,DEF_colors_sel_4)
  
  cat("COLORS_DEF_1\n")
  cat(str(COLORS_DEF))
  cat("\n")
  cat(sprintf(as.character(unique(COLORS_DEF$Annotation))))
  cat("\n")
  cat(sprintf(as.character(unique(COLORS_DEF$colors))))
  cat("\n")
  
  #### READ and transform tracking_variants ----
  
  tracking_variants = unlist(strsplit(opt$tracking_variants, split=","))
  
  cat("tracking_variants_\n")
  cat(sprintf(as.character(tracking_variants)))
  cat("\n")
  
  #### Read Master_file ----
  
  Master_file<-as.data.frame(readRDS(file=opt$Master_file) , stringsAsFactors=F)
  
  Master_file$variable<-factor(Master_file$variable,
                               levels=rev(levels(Master_file$variable)),
                                          ordered=T)
  
  cat("Master_file_0\n")
  cat(str(Master_file))
  cat("\n")
  cat(str(unique(Master_file$VAR)))
  cat("\n")
  
  #### Read Stats_file ----
  
  Stats_file<-as.data.frame(readRDS(file=opt$Stats_file) , stringsAsFactors=F)
  
  cat("Stats_file_0\n")
  cat(str(Stats_file))
  cat("\n")
  cat(sprintf(as.character(unique(Stats_file$Annotation))))
  cat("\n")
  
 
  #### Read Master_file ----
  
  Master_file<-as.data.frame(readRDS(file=opt$Master_file) , stringsAsFactors=F)
  
  cat("Master_file_0\n")
  cat(str(Master_file))
  cat("\n")
  cat(str(unique(Master_file$VAR)))
  cat("\n")
  
  ############## Graph LOOP -----------------
  
  
  ####   LOOP Through category_vector ----
  
  category_vector<-c("Fig1_Annot_Category","MPRA_CLASS","M_and_M","Multi_Lineage")
  
  
  List_DEF<-list()
  
  CONDITION_DEBUG<-1
  
  for(i in 1:length(category_vector))
  {
    category_sel<-category_vector[i]
    
    cat("----------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(category_sel)))
    cat("\n")
    
    Stats_file_single_sel<-droplevels(Stats_file[which(Stats_file$Annotation ==  category_sel &
                                                         Stats_file$variable%in%df_subset_variables$variable),])
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Stats_file_single_sel_0\n")
      cat(str(Stats_file_single_sel))
      cat("\n")
      cat(str(unique(Stats_file_single_sel$Annotation)))
      cat("\n")
    }
    
    
    categories_in_stats<-unique(as.character(c(Stats_file_single_sel$category2,Stats_file_single_sel$category1)))
    
    if(CONDITION_DEBUG == 1)
    {
      cat("categories_in_stats_0\n")
      cat(sprintf(as.character(categories_in_stats)))
      cat("\n")
    }
    
    ### Table_of_labels
    
    indx.category_sel<-which(colnames(Table_of_labels) == category_sel)

    if(CONDITION_DEBUG == 1)
    {
      cat("indx.category_sel\n")
      cat(str(indx.category_sel))
      cat("\n")

    }

    ind.int<-c(which(colnames(Table_of_labels) == "VAR"),indx.category_sel)

   

    Table_of_labels_subset<-unique(Table_of_labels[,ind.int])
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Table_of_labels_subset_0\n")
      cat(str(Table_of_labels_subset))
      cat("\n")
      
    }
    
    indx.category_sel_subset<-which(colnames(Table_of_labels_subset) == category_sel)
    
    if(CONDITION_DEBUG == 1)
    {
      cat("indx.category_sel_subset_0\n")
      cat(str(indx.category_sel_subset))
      cat("\n")
      
    }
    
    Table_of_labels_subset_restricted<-droplevels(Table_of_labels_subset[which(Table_of_labels_subset[,indx.category_sel_subset]%in%categories_in_stats),])
    
   
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
    
    indx.category_sel_subset<-which(colnames(Table_of_labels_subset_restricted) == category_sel)
    
    if(CONDITION_DEBUG == 1)
    {
      cat("indx.category_sel_subset_1\n")
      cat(str(indx.category_sel_subset))
      cat("\n")
      
    }
    
    if(category_sel == "Multi_Lineage")
    {
      Table_of_labels_subset_restricted[,indx.category_sel_subset]<-factor(Table_of_labels_subset_restricted[,indx.category_sel_subset],
                                                                           levels=rev(levels(Table_of_labels_subset_restricted[,indx.category_sel_subset])),
                                                                           ordered=T)
      
      

    }#category_sel == "Multi_Lineage")
    
    order_levels_category<-levels(Table_of_labels_subset_restricted[,indx.category_sel_subset])
    
    if(CONDITION_DEBUG == 1)
    {
      cat("order_levels_category_0\n")
      cat(str(order_levels_category))
      cat("\n")
    }
    
    

   #### DotPlot ----
    
   
    Stats_file_single_sel$variable<-factor(Stats_file_single_sel$variable,
                                           levels=rev(df_subset_variables$variable),
                                           ordered=T)
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Stats_file_single_sel_1\n")
      cat(str(Stats_file_single_sel))
      cat("\n")
      cat(sprintf(as.character(unique(Stats_file_single_sel$variable))))
      cat("\n")
    }
    
    if(category_sel == "Fig1_Annot_Category")
    {
      
      Stats_file_single_sel<-Stats_file_single_sel[which(Stats_file_single_sel$category1 == "index_variants" |
                                                           Stats_file_single_sel$category2 == "index_variants"),]
      
    }
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Stats_file_single_sel_2\n")
      cat(str(Stats_file_single_sel))
      cat("\n")
      cat(str(unique(Stats_file_single_sel$Annotation)))
      cat("\n")
    }
    
    indx.int<-c(which(colnames(Stats_file_single_sel) == "category2"),which(colnames(Stats_file_single_sel) == "category1"),which(colnames(Stats_file_single_sel) == "variable"),which(colnames(Stats_file_single_sel) == "MINUS_logpval"))
    
    Stats_file_single_sel_subset<-unique(Stats_file_single_sel[,indx.int])
    
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Stats_file_single_sel_subset_0\n")
      cat(str(Stats_file_single_sel_subset))
      cat("\n")
    }
    
    
    
    Stats_file_single_sel_subset$x_lab<-NA
    
    Stats_file_single_sel_subset$x_lab[which(Stats_file_single_sel_subset$category2 != order_levels_category[length(order_levels_category)])]<-as.character(Stats_file_single_sel_subset$category2[which(Stats_file_single_sel_subset$category2 != order_levels_category[length(order_levels_category)])])
    Stats_file_single_sel_subset$x_lab[which(Stats_file_single_sel_subset$category1 != order_levels_category[length(order_levels_category)])]<-as.character(Stats_file_single_sel_subset$category1[which(Stats_file_single_sel_subset$category1 != order_levels_category[length(order_levels_category)])])
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Stats_file_single_sel_subset_1\n")
      cat(str(Stats_file_single_sel_subset))
      cat("\n")
    }
    
    Stats_file_single_sel_subset$z_lab<-NA
    
    Stats_file_single_sel_subset$z_lab[which(Stats_file_single_sel_subset$category2 == order_levels_category[length(order_levels_category)])]<-as.character(Stats_file_single_sel_subset$category2[which(Stats_file_single_sel_subset$category2 == order_levels_category[length(order_levels_category)])])
    Stats_file_single_sel_subset$z_lab[which(Stats_file_single_sel_subset$category1 == order_levels_category[length(order_levels_category)])]<-as.character(Stats_file_single_sel_subset$category1[which(Stats_file_single_sel_subset$category1 == order_levels_category[length(order_levels_category)])])
    
    order_levels_category_RMV<-order_levels_category[-length(order_levels_category)]
    
    
    if(CONDITION_DEBUG == 1)
    {
      cat("order_levels_category_RMV_0\n")
      cat(str(order_levels_category_RMV))
      cat("\n")
    }
    
    Stats_file_single_sel_subset$x_lab<-factor(Stats_file_single_sel_subset$x_lab,
                                               levels=order_levels_category_RMV,
                                               ordered=T)
    
    if(CONDITION_DEBUG == 1)
    {
      cat("Stats_file_single_sel_subset_2\n")
      cat(str(Stats_file_single_sel_subset))
      cat("\n")
    }
    
    Stats_file_single_sel_subset<-Stats_file_single_sel_subset[order(Stats_file_single_sel_subset$variable,Stats_file_single_sel_subset$x_lab),]
    
    Stats_file_single_sel_subset$Significance<-NA
    
    
    Stats_file_single_sel_subset$Significance[which(Stats_file_single_sel_subset$MINUS_logpval >= 1.3)]<-"SIG"
    Stats_file_single_sel_subset$Significance[which(Stats_file_single_sel_subset$MINUS_logpval < 1.3)]<-"NO_SIG"
    
    
    
    # setwd(out)
    # 
    # write.table(Stats_file_single_sel_subset, file="test.tsv", sep="\t", quote=F,row.names=F)
    
    graph_path<-paste(out,'violin_plots_Z_score','/',category_sel,'/',sep='')
    
    if (file.exists(graph_path)){
      
    }else{
      
      dir.create(file.path(graph_path))
    }
    
    setwd(graph_path)
    
    
    
    dotplot<-ggplot(data=Stats_file_single_sel_subset,
                    aes(y=variable,
                        x=z_lab,
                        color=Significance)) +
      geom_point(aes(size=MINUS_logpval), stroke=1)+
      scale_size(range = c(0,10), name='-log10pval') +
      scale_color_manual(values=c("gray","black"),name='Significant') + 
      scale_y_discrete(name=NULL, drop=F)+
      scale_x_discrete(name=NULL, drop=F)+
      facet_wrap('x_lab', nrow=1, scales='free_x')+
      theme_bw()+
      theme_classic()+
      theme(axis.title.y=element_text(size=16, family="sans"),
            axis.text.y=element_text(angle=0,size=16, color="black", family="sans"),
            axis.text.x=element_blank())+
      ggeasy::easy_center_title()
    
  
    
    svgname<-paste("dotplot_value_Z_score",".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= dotplot,
             device="svg",
             height=10, width=12)
    }
    
    
    
    # quit(status = 1)
    
    
    #### Representation by subgroups of categories ----
    
    array_GROUP<-unique(df_subset_variables$GROUP)
    
    if(CONDITION_DEBUG == 1)
    {
      cat("array_GROUP\n")
      cat(str(array_GROUP))
      cat("\n")
      
    }
    
    
    list_array_GROUP<-list()
    
    for(k in 1:length(array_GROUP))
    {
      array_GROUP_sel<-array_GROUP[k]
      
      cat("------------>\t")
      cat(sprintf(as.character(array_GROUP_sel)))
      cat("\n")
      
      df_subset_variables_sel<-df_subset_variables[which(df_subset_variables$GROUP == array_GROUP_sel),]
      
      cat("df_subset_variables_sel\n")
      cat(str(df_subset_variables_sel))
      cat("\n")
      cat(sprintf(as.character(unique(df_subset_variables_sel$variable))))
      cat("\n")
      
      
      Master_file_sel<-droplevels(Master_file[which(Master_file$variable%in%df_subset_variables_sel$variable),])
      
      if(CONDITION_DEBUG == 1)
      {
        cat("Master_file_sel_0\n")
        cat(str(Master_file_sel))
        cat("\n")
        cat(str(unique(Master_file_sel$VAR)))
        cat("\n")
        cat(sprintf(as.character(unique(Master_file_sel$variable))))
        cat("\n")
        cat(sprintf(as.character(names(summary(Master_file_sel$value_Z_score)))))
        cat("\n")
        cat(sprintf(as.character(summary(Master_file_sel$value_Z_score))))
        cat("\n")
      }
      
      ### Add labels
      
      Master_file_sel<-droplevels(merge(Master_file_sel,
                             Table_of_labels_subset_restricted,
                             by="VAR"))
      
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
      
      Master_file_sel$variable<-factor(Master_file_sel$variable,
                                       levels=rev(df_subset_variables_sel$variable),
                                       ordered=T)
      
      if(CONDITION_DEBUG == 1)
      {
        cat("Master_file_sel_2\n")
        cat(str(Master_file_sel))
        cat("\n")
        cat(str(unique(Master_file_sel$VAR)))
        cat("\n")
        cat(sprintf(as.character(names(summary(Master_file_sel$variable)))))
        cat("\n")
        cat(sprintf(as.character(summary(Master_file_sel$variable))))
        cat("\n")
      }
      
      
      ind.cat<-c(which(colnames(Master_file_sel) == category_sel))
      
      if(CONDITION_DEBUG == 1)
      {
        cat("ind.cat\n")
        cat(str(ind.cat))
        cat("\n")
      }
    
      
      check_NA<-Master_file_sel[is.na(Master_file_sel$value),]
      
      if(CONDITION_DEBUG == 1)
      {
        cat("check_NA_0\n")
        cat(str(check_NA))
        cat("\n")
        cat(str(unique(check_NA$VAR)))
        cat("\n")
        cat(sprintf(as.character(names(summary(check_NA$value_Z_score)))))
        cat("\n")
        cat(sprintf(as.character(summary(check_NA$value_Z_score))))
        cat("\n")
      }
      
      
      
      ########################## Violin plot ------------------------
      
      
      COLORS_DEF_sel<-COLORS_DEF[which(COLORS_DEF$Annotation == category_sel),]
      
      if(CONDITION_DEBUG == 1)
      {
        cat("COLORS_DEF_sel_0\n")
        cat(str(COLORS_DEF_sel))
        cat("\n")
      }
      
      if(array_GROUP_sel == "Variant_based_features")
      {
        breaks.Rank<-unique(sort(c(0,seq(from= -3, to= 5,by=2))))
        labels.Rank<-as.character(breaks.Rank)
      }
      if(array_GROUP_sel == "Gene_based_features")
      {
        breaks.Rank<-unique(sort(c(0,seq(from= -1, to= 2.3,by=1))))
        labels.Rank<-as.character(breaks.Rank)
      }
      
      if(CONDITION_DEBUG == 1)
      {
        cat("labels.Rank_0\n")
        cat(str(labels.Rank))
        cat("\n")
      }
      
      if(CONDITION_DEBUG == 1)
      {
        cat("Master_file_sel_0\n")
        cat(str(Master_file_sel))
        cat("\n")
      }
      
      violin_plot <- ggplot(Master_file_sel, aes(value_Z_score, 
                                                variable, 
                                                fill = Master_file_sel[,ind.cat])) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        stat_summary(fun = median, fun.min = median, fun.max = median,
                     geom = "crossbar", width = 0.25)+
        scale_x_continuous(name="Z-Score normalised value",
                           breaks=breaks.Rank,labels=labels.Rank, 
                           limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
        scale_y_discrete(name=NULL, drop=F)+
        facet_grid(cols = vars(Master_file_sel[,ind.cat]), scales = "free") +
        theme_cowplot(font_size = 6) +
        theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              plot.title = element_text(hjust = 0.5),
              panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              strip.text.y.left = element_text(angle = 0)) +
        theme(axis.title.y=element_text(size=24, family="sans"),
              axis.title.x=element_text(size=18, family="sans"),
              axis.text.y=element_text(angle=0,size=4, color="black", family="sans"),
              axis.text.x=element_text(angle=0, size=14, color="black", family="sans"),
              legend.title=element_text(size=16,color="black", family="sans"),
              legend.text=element_text(size=12,color="black", family="sans"))+
        scale_fill_manual(values = COLORS_DEF_sel$colors,
                          drop=F)+
        theme(legend.position = "hidden")
      
      cat("violin_plot_0\n")
      
      if(array_GROUP_sel == "Variant_based_features")
      {
        cat("violin_plot_0_Variant_based_features\n")
        graph_DEF<-plot_grid(violin_plot,NULL,
                             nrow = 2,
                             ncol = 1,
                             rel_heights=c(8,1.5))
      }
      
      if(array_GROUP_sel == "Gene_based_features")
      {
        cat("violin_plot_0_Gene_based_features\n")
        
        violin_plot <- violin_plot +
          theme(legend.position = "none", panel.spacing = unit(0, "lines"),
                plot.title = element_text(hjust = 0.5),
                panel.background = element_rect(fill = NA, color = "black"),
                strip.background = element_blank(),
                strip.text = element_blank(),
                strip.text.y.left = element_blank()) +
          theme(axis.title.x=element_blank())
        
        graph_DEF<-plot_grid(NULL,violin_plot,
                             nrow = 2,
                             ncol = 1,
                             rel_heights=c(8,1.5))
      }
      
    
      
      setwd(graph_path)
      
      svgname<-paste("violin_plot_value_Z_score","_",array_GROUP_sel,".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= graph_DEF,
               device="svg",
               height=10, width=12)
      }
      
      cat("END\n")
      
    }#k in 1:length(array_GROUP)
  }#i category_sel

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
    make_option(c("--Stats_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--subset_variables_1"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--subset_variables_2"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Categories_colors"), type="character", default=NULL, 
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
