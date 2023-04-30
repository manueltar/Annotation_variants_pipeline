
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

Data_wrangling_GENE_PLOEUF = function(option_list)
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
  
  
  #### Read GENE_PLOEUF ----
  
  GENE_PLOEUF = read.table(opt$GENE_PLOEUF, sep="\t", stringsAsFactors = F, header = T)
  
  colnames(GENE_PLOEUF)[which(colnames(GENE_PLOEUF) == "transcript")]<-"transcript_id"
  
  cat("GENE_PLOEUF_0\n")
  str(GENE_PLOEUF)
  cat("\n")
  cat(str(unique(GENE_PLOEUF$gene))) #19658
  cat("\n")
  cat(str(unique(GENE_PLOEUF$transcript_id)))#19704
  cat("\n")
  cat(sprintf(as.character(names(summary(GENE_PLOEUF$oe_lof)))))
  cat("\n")
  cat(sprintf(as.character(summary(GENE_PLOEUF$oe_lof))))
  cat("\n")
  
  GENE_PLOEUF_NO_NA<-GENE_PLOEUF[!is.na(GENE_PLOEUF$oe_lof),]
  
  cat("GENE_PLOEUF_NO_NA_0\n")
  str(GENE_PLOEUF_NO_NA)
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA$gene))) #19155
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA$transcript_id))) #19197
  cat("\n")
  cat(sprintf(as.character(names(summary(GENE_PLOEUF_NO_NA$oe_lof)))))
  cat("\n")
  cat(sprintf(as.character(summary(GENE_PLOEUF_NO_NA$oe_lof))))
  cat("\n")
  
  
  #### Read TRANSCRIPTS_table file ----
  
  TRANSCRIPTS_table<-as.data.frame(fread(file=opt$TRANSCRIPTS_table,sep="\t") , stringsAsFactors=F)
  
  cat("TRANSCRIPTS_table\n")
  cat(str(TRANSCRIPTS_table))
  cat("\n")
  
  colnames(TRANSCRIPTS_table)<-c("ensembl_gene_id","HGNC","transcript_id","transcript_name")
  
  
  indx.int<-c(which(colnames(TRANSCRIPTS_table) == "ensembl_gene_id"),which(colnames(TRANSCRIPTS_table) == "HGNC"),which(colnames(TRANSCRIPTS_table) == "transcript_id"))
  
  
  TRANSCRIPTS_table_subset<-unique(TRANSCRIPTS_table[,indx.int])
  
  
  cat("TRANSCRIPTS_table_subset\n")
  cat(str(TRANSCRIPTS_table_subset))
  cat("\n")
  
  #### merge pLOEUF and gene ID ----
  
  GENE_PLOEUF_NO_NA<-merge(GENE_PLOEUF_NO_NA,
                     TRANSCRIPTS_table_subset,
                     by="transcript_id")
  
  cat("GENE_PLOEUF_NO_NA_1\n")
  str(GENE_PLOEUF_NO_NA)
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA$gene))) #18223
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA$transcript_id))) #18223
  cat("\n")
  cat(sprintf(as.character(names(summary(GENE_PLOEUF_NO_NA$oe_lof)))))
  cat("\n")
  cat(sprintf(as.character(summary(GENE_PLOEUF_NO_NA$oe_lof))))
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
  cat(str(unique(VEP_and_PCHiC$ensembl_gene_id)))
  cat("\n")
 
  #### Merge with GENE_PLOEUF_NO_NA ----
  
  
  GENE_PLOEUF_NO_NA<-merge(GENE_PLOEUF_NO_NA,
                     VEP_and_PCHiC,
                     by=c("ensembl_gene_id","HGNC"))
  
  cat("GENE_PLOEUF_NO_NA_2\n")
  str(GENE_PLOEUF_NO_NA)
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA$VAR)))
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA$gene)))
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA$transcript_id)))
  cat("\n")
  
  
  indx.int<-c(which(colnames(GENE_PLOEUF_NO_NA) == "VAR"),which(colnames(GENE_PLOEUF_NO_NA) == "ensembl_gene_id"),which(colnames(GENE_PLOEUF_NO_NA) == "HGNC"),which(colnames(GENE_PLOEUF_NO_NA) == "oe_lof"))
  
  GENE_PLOEUF_NO_NA_subset<-unique(GENE_PLOEUF_NO_NA[,indx.int])
  
  cat("GENE_PLOEUF_NO_NA_subset_0\n")
  str(GENE_PLOEUF_NO_NA_subset)
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA_subset$VAR)))
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA_subset$ensembl_gene_id)))
  cat("\n")
  cat(sprintf(as.character(names(summary(GENE_PLOEUF_NO_NA_subset$oe_lof)))))
  cat("\n")
  cat(sprintf(as.character(summary(GENE_PLOEUF_NO_NA_subset$oe_lof))))
  cat("\n")
 
  
  #### melt and Z.score ----
  
  GENE_PLOEUF_NO_NA_subset.m<-melt(GENE_PLOEUF_NO_NA_subset, id.vars=c("VAR","HGNC","ensembl_gene_id"))
  
  GENE_PLOEUF_NO_NA_subset.m$variable<-as.character(GENE_PLOEUF_NO_NA_subset.m$variable)
  
  GENE_PLOEUF_NO_NA_subset.m$variable<-'oe_lof'
  
  cat("GENE_PLOEUF_NO_NA_subset.m_0\n")
  str(GENE_PLOEUF_NO_NA_subset.m)
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA_subset.m$VAR)))
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA_subset.m$ensembl_gene_id)))
  cat("\n")
  cat(sprintf(as.character(names(summary(GENE_PLOEUF_NO_NA_subset.m$value)))))
  cat("\n")
  cat(sprintf(as.character(summary(GENE_PLOEUF_NO_NA_subset.m$value))))
  cat("\n")
  
  
  #### Zscore ----
  
  mean_LOEUF<-mean(GENE_PLOEUF_NO_NA_subset.m$value)
  sd_LOEUF<-sd(GENE_PLOEUF_NO_NA_subset.m$value)
  
  cat("mean_LOEUF\n")
  str(mean_LOEUF)
  cat("\n")
  str(sd_LOEUF)
  cat("\n")
  
  GENE_PLOEUF_NO_NA_subset.m$value_Z_score<-(GENE_PLOEUF_NO_NA_subset.m$value-mean_LOEUF)/sd_LOEUF
  
  cat("GENE_PLOEUF_NO_NA_subset.m_1\n")
  str(GENE_PLOEUF_NO_NA_subset.m)
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA_subset.m$VAR)))
  cat("\n")
  cat(str(unique(GENE_PLOEUF_NO_NA_subset.m$ensembl_gene_id)))
  cat("\n")
  cat(sprintf(as.character(names(summary(GENE_PLOEUF_NO_NA_subset.m$value)))))
  cat("\n")
  cat(sprintf(as.character(summary(GENE_PLOEUF_NO_NA_subset.m$value))))
  cat("\n")
  cat(sprintf(as.character(names(summary(GENE_PLOEUF_NO_NA_subset.m$value_Z_score)))))
  cat("\n")
  cat(sprintf(as.character(summary(GENE_PLOEUF_NO_NA_subset.m$value_Z_score))))
  cat("\n")
  
  check<-GENE_PLOEUF_NO_NA_subset.m[which(GENE_PLOEUF_NO_NA_subset.m$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  # #################################################
  # quit(status = 1)
  # 
 
  # #### SAVE ----
  
  setwd(out)
  
  write.table(GENE_PLOEUF_NO_NA_subset.m, file="GENE_PLOEUF_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
  saveRDS(GENE_PLOEUF_NO_NA_subset.m, file="Prepared_file_GENE_PLOEUF.rds")
  
  
  write.table(check, file="check_GENE_PLOEUF_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
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
    make_option(c("--GENE_PLOEUF"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TRANSCRIPTS_table"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VEP_CSQ"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--PCHiC"), type="character", default=NULL, 
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
  
  Data_wrangling_GENE_PLOEUF(opt)
  
}


###########################################################################

system.time( main() )
