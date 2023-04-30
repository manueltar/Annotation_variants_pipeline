
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

Data_wrangling_chromstates = function(option_list)
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
  
 
  
  #### READ and transform desiR_weights ----
  
  desiR_weights = as.numeric(unlist(strsplit(opt$desiR_weights, split=",")))
  
  cat("desiR_weights_\n")
  cat(sprintf(as.character(desiR_weights)))
  cat("\n")
  cat(str(desiR_weights))
  cat("\n")
  
 
  
  #### Read chromstates_pre_ranked file ----
  
  chromstates_pre_ranked<-as.data.frame(fread(file=opt$chromstates_pre_ranked,sep="\t") , stringsAsFactors=F)
  
  cat("chromstates_pre_ranked_0\n")
  cat(str(chromstates_pre_ranked))
  cat("\n")
  cat(str(unique(chromstates_pre_ranked$VAR)))
  cat("\n")
  
  
  #### Aggregate Chicago weights by VAR and Tag  ----
  
  CONDITION_DEBUG<-1
  
  chromstates_pre_ranked.dt<-data.table(chromstates_pre_ranked, key=c("VAR"))
  
  
  Aggregation_table<-as.data.frame(chromstates_pre_ranked.dt[,.(Total_Aggregate_chromstates=sum(Aggregate_chromstates_FINAL),
                                                                 nPhenotypes=.N), by=key(chromstates_pre_ranked.dt)], stringsAsFactors=F)
  
  
  Aggregation_table$normalised_Total_Aggregate_chromstates<-Aggregation_table$Total_Aggregate_chromstates/Aggregation_table$nPhenotypes
  
  if(CONDITION_DEBUG == 1)
  {
    cat("Aggregation_table_0\n")
    cat(str(Aggregation_table))
    cat("\n")
    cat(str(unique(Aggregation_table$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(Aggregation_table$Total_Aggregate_chromstates)))))
    cat("\n")
    cat(sprintf(as.character(summary(Aggregation_table$Total_Aggregate_chromstates))))
    cat("\n")
    cat(sprintf(as.character(names(summary(Aggregation_table$nPhenotypes)))))
    cat("\n")
    cat(sprintf(as.character(summary(Aggregation_table$nPhenotypes))))
    cat("\n")
    
    cat(sprintf(as.character(names(summary(Aggregation_table$normalised_Total_Aggregate_chromstates)))))
    cat("\n")
    cat(sprintf(as.character(summary(Aggregation_table$normalised_Total_Aggregate_chromstates))))
    cat("\n")
  }
  
  
  check<-Aggregation_table[which(Aggregation_table$VAR%in%tracking_variants),]
  
  if(CONDITION_DEBUG == 1)
  {
    cat("check_0\n")
    cat(str(check))
    cat("\n")
    cat(str(unique(check$VAR)))
    cat("\n")
  }
  
  
  
  
  
  
  
  
  ######################### SAVE -----
  
  setwd(out)
  

  write.table(check, file="check_chromstates_GLOBAL_preranked.tsv", sep="\t", quote = F, row.names = F)
  
  
  
  
  ################################## desiR chromstates ---------------------------------------------
  
  graph_path<-paste(out,'desiR_graphs','/',sep='')
  setwd(graph_path)
  
  pdf(file="desiR_normalised_Total_Aggregate_chromstates_FINAL_PRE.pdf")
  hist(Aggregation_table$normalised_Total_Aggregate_chromstates, breaks=50, col="grey", border="white", main="",
       xlab="normalised_Total_Aggregate_chromstates_FINAL")
  dev.off()
  
  
  normalised_Total_Aggregate_chromstates_FINAL_LOW<-desiR_weights[1]
  normalised_Total_Aggregate_chromstates_FINAL_HIGH<-desiR_weights[2]
  
  cat("desiR_weights_chromstates\n")
  cat(sprintf(as.character(c(normalised_Total_Aggregate_chromstates_FINAL_LOW,normalised_Total_Aggregate_chromstates_FINAL_HIGH))))
  cat("\n")
 
    
  pdf(file="desIR_normalised_Total_Aggregate_chromstates_FINAL.pdf")
  hist(Aggregation_table$normalised_Total_Aggregate_chromstates, breaks=50, col="grey", border="white", main="",
       xlab="chromstates-Seq in a Relevant Cell Type")
  des.line(Aggregation_table$normalised_Total_Aggregate_chromstates, "d.high", des.args=c(cut1=normalised_Total_Aggregate_chromstates_FINAL_LOW, cut2=normalised_Total_Aggregate_chromstates_FINAL_HIGH, scale=0.5))
  dev.off()
  
  
 
  Aggregation_table$normalised_Total_Aggregate_chromstates_component <- d.high(Aggregation_table$normalised_Total_Aggregate_chromstates, cut1=normalised_Total_Aggregate_chromstates_FINAL_LOW, cut2=normalised_Total_Aggregate_chromstates_FINAL_HIGH, scale=0.5)
  
  cat("Aggregation_table_1\n")
  cat(str(Aggregation_table))
  cat("\n")
  cat(sprintf(as.character(names(summary(Aggregation_table$normalised_Total_Aggregate_chromstates_component)))))
  cat("\n")
  cat(sprintf(as.character(summary(Aggregation_table$normalised_Total_Aggregate_chromstates_component))))
  cat("\n")
  
  #### Overall desirability ---- 
  
  
  Overall_normalised_Total_Aggregate_chromstates_FINAL<-desiR_weights[3]
  
  cat("desiR_weights_Overall\n")
  cat(sprintf(as.character(c(Overall_normalised_Total_Aggregate_chromstates_FINAL))))
  cat("\n")
  
  
  
  Aggregation_table$Overall_weight <- d.overall(Aggregation_table$normalised_Total_Aggregate_chromstates_component, 
                                          weights=c(Overall_normalised_Total_Aggregate_chromstates_FINAL))
  
  setwd(graph_path)
  
  par(las=1)
  pdf(file= "desIR_chromstates_overall.pdf")
  plot(rev(sort(Aggregation_table$Overall_weight)), type="l", xlab="Rank", ylab="Overall Desirability")
  dev.off()
  
  cat("Aggregation_table_2\n")
  cat(str(Aggregation_table))
  cat("\n")
  cat(sprintf(as.character(names(summary(Aggregation_table$Overall_weight)))))
  cat("\n")
  cat(sprintf(as.character(summary(Aggregation_table$Overall_weight))))
  cat("\n")
  
  
 
  
  
  check<-Aggregation_table[which(Aggregation_table$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(check$VAR)))
  cat("\n")
  cat(sprintf(as.character(check$nPhenotypes)))
  cat("\n")
  cat(sprintf(as.character(check$normalised_Total_Aggregate_chromstates_FINAL)))
  cat("\n")

  
  # #### SAVE ----
  
  setwd(out)
  
  write.table(Aggregation_table, file="chromstates_GLOBAL_Ranked.tsv", sep="\t", quote = F, row.names = F)
  
  write.table(check, file="check_chromstates_GLOBAL_Ranked.tsv", sep="\t", quote = F, row.names = F)
  
  # ################################################################
  # quit(status = 1)
  
}

Z_score_normalization_and_convergence = function(option_list)
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
  
  
  #### Read chromstates_ranked file ----
  
  setwd(out)
  
  chromstates_ranked<-as.data.frame(fread(file="chromstates_GLOBAL_Ranked.tsv",sep="\t") , stringsAsFactors=F)
  
  cat("chromstates_ranked_0\n")
  cat(str(chromstates_ranked))
  cat("\n")
  cat(str(unique(chromstates_ranked$VAR)))
  cat("\n")
  
 
  #### calculate mean and sd -----
  
  chromstates_ranked$mean_Overall_weight<-mean(chromstates_ranked$Overall_weight, na.rm =T)
  chromstates_ranked$sd_Overall_weight<-sd(chromstates_ranked$Overall_weight, na.rm =T)
  
  
  cat("chromstates_ranked_POST_merge\n")
  cat(str(chromstates_ranked))
  cat("\n")
  cat(str(unique(chromstates_ranked$VAR)))
  cat("\n")
  
  
  chromstates_ranked$Overall_weight_Z_score<-(chromstates_ranked$Overall_weight-chromstates_ranked$mean_Overall_weight)/chromstates_ranked$sd_Overall_weight
  
  
  chromstates_ranked$Overall_weight_Z_score<-(chromstates_ranked$Overall_weight-chromstates_ranked$mean_Overall_weight)/chromstates_ranked$sd_Overall_weight
  
  cat("chromstates_ranked_POST_merge_Z_score\n")
  cat(str(chromstates_ranked))
  cat("\n")
  cat(str(unique(chromstates_ranked$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(chromstates_ranked$Overall_weight_Z_score)))))
  cat("\n")
  cat(sprintf(as.character(summary(chromstates_ranked$Overall_weight_Z_score))))
  cat("\n")
  
  #### select columns and change names to converge ----
  
  indx.int<-c(which(colnames(chromstates_ranked) == "VAR"),which(colnames(chromstates_ranked) == "Overall_weight"),which(colnames(chromstates_ranked) == "Overall_weight_Z_score"))
  
  chromstates_ranked_subset<-unique(chromstates_ranked[,indx.int])
  
  cat("chromstates_ranked_subset_0\n")
  cat(str(chromstates_ranked_subset))
  cat("\n")
  cat(str(unique(chromstates_ranked_subset$VAR)))
  cat("\n")
  
 
  colnames(chromstates_ranked_subset)[which(colnames(chromstates_ranked_subset) == "Overall_weight")]<-"value"
  colnames(chromstates_ranked_subset)[which(colnames(chromstates_ranked_subset) == "Overall_weight_Z_score")]<-"value_Z_score"
  
  chromstates_ranked_subset$variable<-paste('Rank_chromstates',sep='')
  
  
  cat("chromstates_ranked_subset_1\n")
  cat(str(chromstates_ranked_subset))
  cat("\n")
  cat(str(unique(chromstates_ranked_subset$VAR)))
  cat("\n")
  
  check<-chromstates_ranked_subset[which(chromstates_ranked_subset$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  # #### SAVE ----
  
  setwd(out)
  
  saveRDS(chromstates_ranked_subset, file="Prepared_file_chromstates.rds")
  
  write.table(check, file="check_chromstates_GLOBAL_Ranked.tsv", sep="\t", quote = F, row.names = F)
  
  
  
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
    make_option(c("--chromstates_pre_ranked"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--desiR_weights"), type="character", default=NULL, 
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
  
  Data_wrangling_chromstates(opt)
  Z_score_normalization_and_convergence(opt)
  
}


###########################################################################

system.time( main() )
