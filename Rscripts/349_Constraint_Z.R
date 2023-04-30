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
suppressMessages(library("BiocGenerics", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("S4Vectors", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("IRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Biobase", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicFeatures", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("OrganismDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GO.db", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Homo.sapiens", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("gwascat", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rtracklayer", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("liftOver",lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))




opt = NULL

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
  
  cat("OUT\n")
  cat(sprintf(as.character(out)))
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
  cat(str(unique(ALL_dB$VAR)))
  cat("\n")
  
  
  indx.dep<-c(which(colnames(ALL_dB) == "maf_origin"))
  
  ALL_dB_subset<-unique(ALL_dB[,-indx.dep])
  
  cat("ALL_dB_subset\n")
  cat(str(ALL_dB_subset))
  cat("\n")
  cat(str(unique(ALL_dB_subset$VAR)))
  cat("\n")
  
  rm(ALL_dB)
  
  VAR_df<-unique(data.frame(chr=ALL_dB_subset$chr,
                            pos37=ALL_dB_subset$pos37,
                            ref=ALL_dB_subset$ref,
                            alt=ALL_dB_subset$alt,
                            VAR=ALL_dB_subset$VAR,
                            stringsAsFactors = F))
  
  Condition_DEBUG <- 1
  
  if(Condition_DEBUG == 1)
  {
    cat("VAR_df_\n")
    str(VAR_df)
    cat("\n")
    cat(str(unique(VAR_df$VAR)))
    cat("\n")
  }
  
  gr_VARS <-unique(GRanges(
    seqnames = as.character(VAR_df$chr),
    name2=rep("VAR",length(VAR_df$VAR)),
    ranges=IRanges(
      start=as.numeric(VAR_df$pos37),
      end=as.numeric(VAR_df$pos37),
      names = VAR_df$VAR)))
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("gr_VARS_\n")
    str(gr_VARS)
    cat("\n")
    cat(str(unique(gr_VARS$VAR))) #179324
    cat("\n")
  }
  
  
 
  
  #### LiftOver variants ----
  
 
  ch = import.chain("/nfs/team151/software/manuel_R_ext_data_4_1/hg19ToHg38.over.chain")
  
  seqlevelsStyle(gr_VARS) = "UCSC"  # necessary
  gr_VARS38 = liftOver(gr_VARS, ch)
  gr_VARS38 = unlist(gr_VARS38)
  genome(gr_VARS38) = "hg38"
  
  if(length(gr_VARS38) >0)
  {
    
    chr_38<-as.character(seqnames(gr_VARS38))
    names_38<-as.character(names(gr_VARS38))
    
    ref_VAR38<-gsub("^chr[^_]+_[0-9]+_","",names_38)
    ref_VAR38<-gsub("_.+$","",ref_VAR38)
    
    
    # cat("ref_VAR38\n")
    # cat(sprintf(as.character(ref_VAR38)))
    # cat("\n")
    
    alt_VAR38<-gsub("^chr[^_]+_[0-9]+_[^_]+_","",names_38)
    # alt_VAR38<-gsub("_.+$","",alt_VAR38)
    
    
    # cat("alt_VAR38\n")
    # cat(sprintf(as.character(alt_VAR38)))
    # cat("\n")
    
    
    
    
    VAR_38_df<-data.frame(chr=as.character(seqnames(gr_VARS38)),
                          pos_38=start(gr_VARS38),
                          ref=ref_VAR38,
                          alt=alt_VAR38,
                          VAR=names(gr_VARS38),
                          stringsAsFactors = F)
    
    VAR_38_df$VAR_38<-paste(VAR_38_df$chr,VAR_38_df$pos_38,VAR_38_df$ref,VAR_38_df$alt,sep='_')
    
    if(Condition_DEBUG == 1)
    {
      cat("VAR_38_df_1\n")
      str(VAR_38_df)
      cat("\n")
      cat("VAR_df_\n")
      str(VAR_df)
      cat("\n")
      cat(str(unique(VAR_df$VAR)))
      cat("\n")
    }
    
    
    VAR_ALL_dB_subset_df<-unique(merge(VAR_df,
                                                           VAR_38_df,
                                                           by=c("chr","ref","alt","VAR"),
                                                           all=T))
    
    VAR_ALL_dB_subset_df$VAR_38[is.na(VAR_ALL_dB_subset_df$VAR_38)]<-"ABSENT"
    
    if(Condition_DEBUG == 1)
    {
      cat("VAR_ALL_dB_subset_df_2\n")
      str(VAR_ALL_dB_subset_df)
      cat("\n")
    }
    
    check.ABSENT<-VAR_ALL_dB_subset_df[which(VAR_ALL_dB_subset_df$VAR_38 == "ABSENT"),]
    #
    if(Condition_DEBUG == 1)
    {
      cat("check.ABSENT\n")
      str(check.ABSENT) # 151 absent variants
      cat("\n")
    }
   
    
  }else{
    
    
    stop("NO_LIFT_OVER\n")
    
  }# length(gr_VARS38) >0
  
  
  #### Subset for liftedOver variants ----
  
  
  VAR_ALL_dB_LiftedOver<-VAR_ALL_dB_subset_df[which(VAR_ALL_dB_subset_df$VAR_38 != "ABSENT"),]
  
  if(Condition_DEBUG == 1)
  {
    cat("VAR_ALL_dB_LiftedOver_1\n")
    str(VAR_ALL_dB_LiftedOver)
    cat("\n")
    cat(str(unique(VAR_ALL_dB_LiftedOver$VAR)))
    cat("\n")
  }
  
  gr_VARS_38 <-unique(GRanges(
    seqnames = as.character(VAR_ALL_dB_LiftedOver$chr),
    name2=rep("VAR",length(VAR_ALL_dB_LiftedOver$VAR_38)),
    ranges=IRanges(
      start=as.numeric(VAR_ALL_dB_LiftedOver$pos_38),
      end=as.numeric(VAR_ALL_dB_LiftedOver$pos_38),
      names = VAR_ALL_dB_LiftedOver$VAR_38)))
  
  
  
  # if(Condition_DEBUG == 1)
  # {
  #   cat("gr_VARS_38_0\n")
  #   str(gr_VARS_38)
  #   cat("\n")
  #  
  # }
  
  
  #### Read Constraint_Z ----
  
  
  Constraint_Z<-as.data.frame(fread(file=opt$Constraint_Z,sep="\t",header=T, stringsAsFactors = F))
  
  cat("Constraint_Z\n")
  str(Constraint_Z)
  cat("\n")
  
  gr_Constraint_Z <-unique(GRanges(
    seqnames = as.character(Constraint_Z$chrom),
    name2=rep("1Kb_window",length(Constraint_Z$element_id)),
    ranges=IRanges(
      start=as.numeric(Constraint_Z$start),
      end=as.numeric(Constraint_Z$end),
      names = Constraint_Z$element_id)))
  
  
  
  # if(Condition_DEBUG == 1)
  # {
  #   cat("gr_Constraint_Z_0\n")
  #   str(gr_Constraint_Z)
  #   cat("\n")
  # 
  # }
  # 
  ################## find regions that overlap -------------------
  
  m <- findOverlaps(gr_Constraint_Z,gr_VARS_38)
  
  if(Condition_DEBUG == 1)
  {
    cat("m\n")
    str(m)
    cat("\n")
   
  }
  
  # Add gc.name to subject GRanges (i.e. left join)
  mcols(gr_VARS_38)$gc.name <- names(gr_VARS_38)
  mcols(gr_Constraint_Z)$gc.name <- names(gr_Constraint_Z)
  mcols(gr_VARS_38)[subjectHits(m), "gc.name"] = mcols(gr_Constraint_Z)[queryHits(m), "gc.name"]
  
  # if(Condition_DEBUG == 1)
  # {
  #   cat("gr_VARS_38_1\n")
  #   str(gr_VARS_38)
  #   cat("\n")
  #   
  # }
  # 
  # if(Condition_DEBUG == 1)
  # {
  #   cat("gr_Constraint_Z_1\n")
  #   str(gr_Constraint_Z)
  #   cat("\n")
  #   
  # }
  
  
  df2 <- data.frame(Feature=gr_VARS_38$gc.name,
                    OV_chr=seqnames(gr_VARS_38),
                    OV_start=as.integer(start(gr_VARS_38)-1),
                    OV_end=as.integer(end(gr_VARS_38)),
                    VAR_38=names(gr_VARS_38), stringsAsFactors = F)
  
  cat("df2_PRE\n")
  cat(str(df2))
  cat("\n")
  
  indx.overlaps<-grep("^[^-]+-[^-]+-.+",df2$Feature)
  
  cat("indx.overlaps\n")
  cat(str(indx.overlaps))
  cat("\n")
  
  Selected_overlaps<-df2[indx.overlaps,]
  
  cat("Selected_overlaps_0\n")
  cat(str(Selected_overlaps))
  cat("\n")
  cat(str(unique(Selected_overlaps$names)))
  cat("\n") #113687 variants with overlaps 
  
  Selected_not_overlaps<-df2[-indx.overlaps,]
  
  cat("Selected_not_overlaps\n")
  cat(str(Selected_not_overlaps))
  cat("\n")
  cat(str(unique(Selected_not_overlaps$names)))
  cat("\n") #65486 variants with no overlap
  
  # setwd(out)
  # write.table(Selected_overlaps,file="test.tsv", sep="\t", row.names = F, quote=F)
  # 
  # write.table(Selected_not_overlaps,file="NOT_OVERLAP_test.tsv", sep="\t", row.names = F, quote=F)
  
  #########################################################   Constraint_Z is not genome wide #######################################################################################################
  #########################################################   Constraint_Z is not genome wide #######################################################################################################
  #########################################################   Constraint_Z is not genome wide #######################################################################################################
  
  
  #Z$ head -15145 constraint_z_genome_1kb.qc.download.txt|tail -4
  # chr1    26304000        26305000        chr1-26304000-26305000  2320.0  165.25799501684213      140.0   0.8471602235385466      1.964796761051548
  # chr1    26305000        26306000        chr1-26305000-26306000  1950.0  122.1003267797037       110.0   0.9008984898006422      1.0950619424372166
  # chr1    26314000        26315000        chr1-26314000-26315000  2993.0  204.155544089216        177.0   0.866986007113532       1.900543961162472
  # chr1    26317000        26318000        chr1-26317000-26318000  2077.0  140.6101222139205       118.0   0.8391998964375972      1.9067537841989481
  
  # Constraint_Z is incomplete here it jumps from 26306000 to 26314000 and misses the variant chr1_26308275_A_AAAAAG
  
  # Constraint_Z$ head -132056 constraint_z_genome_1kb.qc.download.txt|tail -4
  # chr1    205287000       205288000       chr1-205287000-205288000        2995.0  197.63908473689685      173.0   0.8753329344259151      1.752621577118162
  # chr1    205292000       205293000       chr1-205292000-205293000        2796.0  205.145711458397        196.0   0.9554184613786004      0.6385373040024369
  
  # Constraint_Z is incomplete here it jumps from 205288000 to 205292000 and misses the variant chr1_205281091_A_AAAAC
  
  
  
  colnames(Selected_overlaps)[which(colnames(Selected_overlaps) == "Feature")]<-"element_id"
  
  
  cat("Selected_overlaps_1\n")
  cat(str(Selected_overlaps))
  cat("\n")
  cat(str(unique(Selected_overlaps$VAR_38)))
  cat("\n") #113687 variants with overlaps 
 
  
  Constraint_Z_hits<-merge(Constraint_Z,
                           Selected_overlaps,
                           by="element_id")
  
  cat("Constraint_Z_hits_0\n")
  cat(str(Constraint_Z_hits))
  cat("\n")
  cat(str(unique(Constraint_Z_hits$VAR_38)))
  cat("\n") # 113687
  cat(sprintf(as.character(names(summary(Constraint_Z_hits$z)))))
  cat("\n")
  cat(sprintf(as.character(summary(Constraint_Z_hits$z))))
  cat("\n")
  
  
  Constraint_Z_hits<-merge(Constraint_Z_hits,
                           VAR_ALL_dB_LiftedOver,
                           by="VAR_38")
  
  
  cat("Constraint_Z_hits_0\n")
  cat(str(Constraint_Z_hits))
  cat("\n")
  cat(str(unique(Constraint_Z_hits$VAR)))
  cat("\n") #
  cat(str(unique(Constraint_Z_hits$VAR_38)))
  cat("\n") #
  cat(sprintf(as.character(names(summary(Constraint_Z_hits$z)))))
  cat("\n")
  cat(sprintf(as.character(summary(Constraint_Z_hits$z))))
  cat("\n")
  
  check<-Constraint_Z_hits[which(Constraint_Z_hits$VAR%in%tracking_variants),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
 
  
  # #### SAVE ----
  
  setwd(out)
  
  write.table(Constraint_Z_hits, file="Constraint_Z_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
  write.table(check, file="check_Constraint_Z_GLOBAL.tsv", sep="\t", quote = F, row.names = F)
  
  ### create the rds version for stats ----
  
  indx.int<-c(which(colnames(Constraint_Z_hits) == "VAR"),which(colnames(Constraint_Z_hits) == "z"))
  
  Constraint_Z_hits_subset<-unique(Constraint_Z_hits[,indx.int])
  
  colnames(Constraint_Z_hits_subset)[which(colnames(Constraint_Z_hits_subset) == 'z')]<-"value"
  Constraint_Z_hits_subset$value_Z_score<-Constraint_Z_hits_subset$value
  
  Constraint_Z_hits_subset$variable<-'Constraint_Z'
  
  cat("Constraint_Z_hits_subset_0\n")
  cat(str(Constraint_Z_hits_subset))
  cat("\n")
  
  setwd(out)
  
  saveRDS(Constraint_Z_hits_subset,file="Prepared_file_Constraint_Z.rds")

 
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
    make_option(c("--Constraint_Z"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ALL_dB"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--tracking_variants"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--AS_PreC"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GENE_EXP"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Rank_file"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GL_COGS_CSQ_EXP"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GL_TWAS_COLOCS"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GL_OpenTargets_score"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Partition_START"), type="numeric", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Partition_END"), type="numeric", default=NULL,
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
  
 
 Data_wrangling(opt)
  
}




###########################################################################

system.time( main() )
