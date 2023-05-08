# Annotation_variants_pipeline


#  Master script command line:

$ bash 313_new_calculator_of_scores.sh \<Path_to_output_files\> \<Memory\> \<Processors\> \<queue\>

e.g. $ bash 313_new_calculator_of_scores.sh Annotation_results/ 4000 1 normal

# R library requirements

"AnnotationDbi","backports","Biobase","BiocGenerics","broom","cli","cowplot","crayon","data.table","desiR","digest","dplyr","farver","GenomeInfoDb","GenomicFeatures","GenomicRanges","ggeasy","ggforce","ggplot2","GO.db","gwascat","Homo.sapiens","IRanges","labeling","liftOver","optparse","OrganismDbi","org.Hs.eg.db","plyr","reshape2","R.methodsS3","R.oo","rstudioapi","rtracklayer","R.utils","S4Vectors","sandwich","splitstackshape","svglite","tidyverse","TxDb.Hsapiens.UCSC.hg19.knownGene","tzdb","withr"


# Subscript: 346_binder_GWAS_and_Z_score.R 

################### Input files: ALL_dB.tsv (annotation of GWAS blood traits)

$ head /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv

VAR     variant_id      rs      chr     pos37   ref     alt     maf_origin      block_no        phenotype       ncondind        credset_size    finemap_beta    finemap_se      finemap_z       finemap_prob    finemap_log10bf info    is_cond_ind     chromstates
chr1_1689164_C_T        4154    rs2235536       chr1    1689164 C       T       0.425205        1       ret     1       107     0.0131952       0.00215928      6.11093 0.00330496      1.44972 0.996973        0       E2C1:E2C3:E2C5:E2C6:E2C7:E2C8:E2C9:E2C10:E2C11:E2C12:E9C14:E2C15:E2C16:E2C17:E9C18:E2C19:E2C20:E2C21:E2C22:E2C23:E2C24:E2C25:E2C26:E2C27:E2C28:E2C29:E2C30:E2C31:E2C32:E2C33:E2C34:E2C35:E2C36:E2C37:E2C38:E2C39

################### Code Main points:

  colnames(ALL_dB)[which(colnames(ALL_dB) == "finemap_prob")]<-"PP" # change name to PP (Posterior Probability)
  colnames(ALL_dB)[which(colnames(ALL_dB) == "maf_origin")]<-"MAF" # change name to MAF (Minimum Allele Frequency)

  ALL_dB$Absolute_effect_size<-abs(ALL_dB$finemap_beta) # calculate the absolute effect size from the effect size

  ALL_dB.dt<-data.table(ALL_dB, key="VAR")
  ALL_dB_MAX_AF<-as.data.frame(ALL_dB.dt[,.SD[which.max(MAF)], by=key(ALL_dB.dt)], stringsAsFactors=F) # Calculate the maximum MAF per variant (some variants have two values slightly different)

  ALL_dB_MAX_AF_subset$mean_value<-mean(ALL_dB_MAX_AF_subset$MAF)
  ALL_dB_MAX_AF_subset$sd_value<-sd(ALL_dB_MAX_AF_subset$MAF)
  ALL_dB_MAX_AF_subset$MAF_Z_score<-(ALL_dB_MAX_AF_subset$MAF-ALL_dB_MAX_AF_subset$mean_value)/ALL_dB_MAX_AF_subset$sd_value # Z score normalisation of the MAF values

  ALL_dB_subset.m<-melt(ALL_dB_subset, id.vars=c("VAR","phenotype"))
  ALL_dB_subset.m.dt<-data.table(ALL_dB_subset.m, key=c("phenotype","variable"))
  ALL_dB_subset.m_GWAS_Type_parameters<-as.data.frame(ALL_dB_subset.m.dt[,.(mean_value=mean(value, na.rm =T),
  sd_value=sd(value, na.rm =T)),
  by=key(ALL_dB_subset.m.dt)], stringsAsFactors=F)

  ALL_dB_subset.m<-merge(ALL_dB_subset.m,
  ALL_dB_subset.m_GWAS_Type_parameters,
  by=c("phenotype","variable"))
  ALL_dB_subset.m$value_Z_score<-(ALL_dB_subset.m$value-ALL_dB_subset.m$mean_value)/ALL_dB_subset.m$sd_value # Z score per phenotype Absolute_effect_size and PP 

  write.table(ALL_dB_subset.m_subset, file="GWAS_GLOBAL_per_traits.tsv", sep="\t", quote = F, row.names = F)
  write.table(ALL_dB_MAX_AF_subset_double, file="MAF_GLOBAL.tsv", sep="\t", quote = F, row.names = F) # Save the two files separately

################### Output files:

$ head GWAS_GLOBAL_per_traits.tsv
VAR     phenotype       variable        value   value_Z_score
chr2_181894386_CT_C     baso    Absolute_effect_size    0.0147284       -0.409226791052507

$ head MAF_GLOBAL.tsv
VAR     value   value_Z_score   variable
chr10_101222300_G_A     0.104047        -1.20875664040154       MAF

# Subscript: 340_binder_CADD.R 

################### Input files: ALL_dB.tsv (annotation of GWAS blood traits)
################### Input files: Vuckovic_VEP_scores.csv.gz (annotation of CADD results)

$ zcat /lustre/scratch123/hgi/teams/soranzo/projects/ALL_dB/csv_files/Vuckovic_VEP_scores.csv.gz|head -2
"variant_id","label","chr","pos","allele","cadd_phred","cadd_raw","conservation"
1,"rs35254660","1",1442413,"T",0.791,-0.185485,""

################### Code Main points:

     Merge_table<-merge(ALL_dB_subset,
			 CADD_result,
			 by=c("chr","pos37","rs")) # Merge CADD result with GWAS results and keep the overlap
      Merge_table_subset_NO_NA<-Merge_table_subset[!is.na(Merge_table_subset$cadd_raw),] # exclude NA's 

################### Output files:

$ head CADD_GLOBAL.tsv
    VAR     cadd_phred      cadd_raw
    chr10_101222300_G_A     0.191   -0.393199
    chr10_101237412_T_C     3.693   0.066067

# Subscript: 341_CADD.R

################### Input files: CADD_GLOBAL.tsv
    VAR     cadd_phred      cadd_raw
    chr10_101222300_G_A     0.191   -0.393199
    chr10_101237412_T_C     3.693   0.066067


################### Code Main points:

 CADD$mean_CADD<-mean(CADD$cadd_raw, na.rm =T)
  CADD$sd_CADD<-sd(CADD$cadd_raw, na.rm =T)
CADD$cadd_raw_Z_score<-(CADD$cadd_raw-CADD$mean_CADD)/CADD$sd_CADD # Z-score normalization across all the variants

################### Output files:

$ Prepared_file_CADD.rds

VAR     value   value_Z_score   variable
chr1_202129205_G_A      1.933218        4.49399301821548        CADD_raw

# Subscript: 342_binder_NCBoost.R

################### Input files: ALL_dB.tsv (annotation of GWAS blood traits)
################### Input files: haemvar_manuel_200903_NCBOOST.tsv (annotation of NCBoost results)

$ head /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NCBoost/haemvar_manuel_200903_NCBOOST.tsv
chr     pos     ref     alt     annovar_annotation      closest_gene_name       gene_type       pLI     familyMemberCount       ncRVIS  ncGERP  RVIS_percentile slr_dnds        GDI     gene_age        UTR3    UTR5    downstream      intergenic      intronic        upstream        GC      CpG     priPhCons       mamPhCons    verPhCons       priPhyloP       mamPhyloP       verPhyloP       GerpN   GerpS   bStatistic      TajimasD_YRI_pvalue     TajimasD_CEU_pvalue     TajimasD_CHB_pvalue     FuLisD_YRI_pvalue       FuLisD_CEU_pvalue       FuLisD_CHB_pvalue       FuLisF_YRI_pvalue       FuLisF_CEU_pvalue       FuLisF_CHB_pvalue   meanDaf1000G     meanHet1000G    meanMAF1000G    meanMAFGnomAD   meanMAF_AFRGnomAD       meanMAF_AMRGnomAD       meanMAF_ASJGnomAD       meanMAF_EASGnomAD       meanMAF_FINGnomAD       meanMAF_NFEGnomAD       meanMAF_OTHGnomAD       CDTS    NCBoost
1       1689164 C       T       intronic        NADK    protein_coding  0.003878536     0       -0.3783625      -0.7206092      45.35858        0.06153 176.4511        0       0       0       0       0       1       0       0.6     0.04    0.028   0.008   0.026   0.457   1.003   1.317   1.86    1.86    662     1.023368     0.1377089       0.08570405      1.024593        0.5245587       0.8809586       1.149016        0.3383352       0.4969653       0.1478  0.01082257      0.00159 0.000709        0.00159 0.000294        0.000432        0.000645        0.000452        0.000302        0.000347        -0.126957       0.06201531
1

################### Code Main points:

  NCBoost_result$VAR<-paste(paste('chr',NCBoost_result$chr, sep=''),NCBoost_result$pos,NCBoost_result$ref,NCBoost_result$alt,sep="_")
 Merge_table<-merge(ALL_dB_subset,
                     NCBoost_result,
                     by=c("VAR")) # Merge NCBoost result with GWAS results and keep the overlap

Merge_table_subset_NO_NA<-Merge_table_subset[!is.na(Merge_table_subset$NCBoost),] # exclude NA's


################### Output files:

$ head NCBoost_GLOBAL.tsv
VAR     closest_gene_name       NCBoost
chr10_101222300_G_A     GOT1    0.015575

# Subscript: 343_NCBoost.R

################### Input files: NCBoost_GLOBAL.tsv
VAR     closest_gene_name       NCBoost
chr10_101222300_G_A     GOT1    0.015575

################### Code Main points:

  NCBoost$mean_NCBoost<-mean(NCBoost$NCBoost, na.rm =T)
  NCBoost$sd_NCBoost<-sd(NCBoost$NCBoost, na.rm =T)
  NCBoost$NCBoost_Z_score<-(NCBoost$NCBoost-NCBoost$mean_NCBoost)/NCBoost$sd_NCBoost # Z-score normalization across all the variants

 indx.int<-c(which(colnames(NCBoost) == "VAR"),which(colnames(NCBoost) == "NCBoost"),which(colnames(NCBoost) == "NCBoost_Z_score"))

  NCBoost_subset<-NCBoost[,indx.int] # Deselect closest_gene_name without doing a unique

################### Output files:

$ Prepared_file_NCBoost.rds

VAR     value   value_Z_score   variable
chr1_202129205_G_A      0.6034472       7.45634719560261        NCBoost
chr12_111844956_C_T     0.7309633       9.18665236911807        NCBoost

# Subscript: 349_Constraint_Z.R

################### Input files: ALL_dB.tsv (annotation of GWAS blood traits)
################### Input files: constraint_z_genome_1kb.qc.download.txt (1 Kb qc ConstraintZ windows nc_constraint_gnomad_v31)

chrom   start   end     element_id      possible        expected        observed        oe      z
chr1    783000  784000  chr1-783000-784000      2443.0  159.09212499779832      142.0   0.8925646068400002      1.3551011656250318

################### Code Main points:

VAR_df<-unique(data.frame(chr=ALL_dB_subset$chr,
                            pos37=ALL_dB_subset$pos37,
                            ref=ALL_dB_subset$ref,
                            alt=ALL_dB_subset$alt,
                            VAR=ALL_dB_subset$VAR,
                            stringsAsFactors = F))
ch = import.chain("/nfs/team151/software/manuel_R_ext_data_4_1/hg19ToHg38.over.chain")
gr_VARS38 = liftOver(gr_VARS, ch) # LiftOver the GWAS variants to GRCh38
check.ABSENT<-VAR_ALL_dB_subset_df[which(VAR_ALL_dB_subset_df$VAR_38 == "ABSENT"),]# 151 absent variants, 151 variants didn't lift over


m <- findOverlaps(gr_Constraint_Z,gr_VARS_38)
  mcols(gr_VARS_38)$gc.name <- names(gr_VARS_38)
  mcols(gr_Constraint_Z)$gc.name <- names(gr_Constraint_Z)
  mcols(gr_VARS_38)[subjectHits(m), "gc.name"] = mcols(gr_Constraint_Z)[queryHits(m), "gc.name"]
 df2 <- data.frame(Feature=gr_VARS_38$gc.name,
                    OV_chr=seqnames(gr_VARS_38),
                    OV_start=as.integer(start(gr_VARS_38)-1),
                    OV_end=as.integer(end(gr_VARS_38)),
                    VAR_38=names(gr_VARS_38), stringsAsFactors = F) # Use GRanges to find overlap between the Constraint Z windows and the variants in GRCh38

indx.overlaps<-grep("^[^-]+-[^-]+-.+",df2$Feature)
Selected_overlaps<-df2[indx.overlaps,] # Select the overlaps between variants and windows. 113687 variants with overlaps and 65486 variants with no overlap
 Constraint_Z_hits<-merge(Constraint_Z_hits,
                           VAR_ALL_dB_LiftedOver,
                           by="VAR_38")

 indx.int<-c(which(colnames(Constraint_Z_hits) == "VAR"),which(colnames(Constraint_Z_hits) == "z"))

  Constraint_Z_hits_subset<-unique(Constraint_Z_hits[,indx.int])

  colnames(Constraint_Z_hits_subset)[which(colnames(Constraint_Z_hits_subset) == 'z')]<-"value"
  Constraint_Z_hits_subset$value_Z_score<-Constraint_Z_hits_subset$value # Because Constraint Z is already a Z score I don't Z score normalise it


################### Output files 1:

$ head Constraint_Z_GLOBAL.tsv
VAR_38  element_id      chrom   start   end     possible        expected        observed        oe      z       OV_chr  OV_start        OV_end  chr     ref     alt     VAR     pos37   pos_38
chr10_100037190_G_A     chr10-100037000-100038000       chr10   100037000       100038000       2383    164.466955961191        142     0.863395319565026       1.75188222761269        chr10   100037189       100037190       chr10   G       A       chr10_101796947_G_A     101796947       100037190

################### Output files 2:

Prepared_file_Constraint_Z.rds

  VAR     value value_Z_score     variable
chr10_101796947_G_A 1.7518822     1.7518822 Constraint_Z
chr10_101797412_A_G 1.7518822     1.7518822 Constraint_Z

# Subscript: 344_binder_SpliceAI.R

################### Input files: ALL_dB.tsv (annotation of GWAS blood traits)
################### Input files: haemvar_manuel_200903_spliceAI_output_fixed.vcf (annotation of SpliceAI results, SpliceAIv1.3.1)

##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
1       11249132        .       G       C       .       .       SpliceAI=C|MTOR|0.00|0.00|0.00|0.00|12|41|-21|41
1       11252716        .       T       C       .       .       SpliceAI=C|MTOR|0.00|0.00|0.00|0.00|46|-15|-26|9,C|ANGPTL7|0.00|0.00|0.00|0.00|-6|13|21|-38

################### Input files: Homo_sapiens.GRCh37.87_GENES_table.txt (ensembl_gene_id equivalence release 87 ENSEMBL)

$ head /nfs/users/nfs_m/mt19/RareVar_Dragana/Homo_sapiens.GRCh37.87_GENES_table.txt
chr1    11869   14412   ENSG00000223972 DDX11L1

################### Code Main points:

SpliceAI_result$VAR<-paste(paste('chr',SpliceAI_result$CHROM,sep=''),SpliceAI_result$POS,SpliceAI_result$REF,SpliceAI_result$ALT, sep="_")
row.with.results<-grep('SpliceAI',SpliceAI_result$INFO)
SpliceAI_result_subset<-SpliceAI_result[row.with.results,] # Add VAR field to SpliceAI results and select the rows in the vcf with SpliceAI predictions

Merge_table<-merge(ALL_dB_subset,
                     SpliceAI_result_subset,
                     by=c("VAR")) # Merge with the variants dataset

Merge_table_separated_LONG<-unique(as.data.frame(cSplit(Merge_table, splitCols = "INFO",
                                        sep = ",", direction = "long", drop = F),stringsAsFactors=F)) # Split long by ',' to separate in different rows with the same variant Splice AI predictions in two different genes
Merge_table_separated_WIDE<-unique(as.data.frame(cSplit(Merge_table_separated_LONG, splitCols = "INFO",
                                                          sep = "|", direction = "wide", drop = F),stringsAsFactors=F)) # Split wide by '|' to separate in different columns all the fields of an Splice AI prediction in a gene

  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_01")]<-"ALLELE"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_02")]<-"HGNC"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_03")]<-"SpliceAI_AG"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_04")]<-"SpliceAI_AL"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_05")]<-"SpliceAI_DG"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_06")]<-"SpliceAI_DL"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_07")]<-"SpliceAI_AG_Pos"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_08")]<-"SpliceAI_AL_Pos"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_09")]<-"SpliceAI_DG_Pos"
  colnames(Merge_table_separated_WIDE)[which(colnames(Merge_table_separated_WIDE) == "INFO_10")]<-"SpliceAI_DL_Pos" # Name the columns with the appropriate field names

################### Output files:

$ head SpliceAI_GLOBAL.tsv
HGNC    VAR     SpliceAI_AG     SpliceAI_AL     SpliceAI_DG     SpliceAI_DL     SpliceAI_AG_Pos SpliceAI_AL_Pos SpliceAI_DG_Pos SpliceAI_DL_Pos ensembl_gene_id
A3GALT2 chr1_33777524_T_C       0.00    0.06    0.00    0.00    -1      -24     -3      -24     ENSG00000184389

# Subscript: 345_SpliceAI.R

################### Input files: SpliceAI_GLOBAL.tsv
HGNC    VAR     SpliceAI_AG     SpliceAI_AL     SpliceAI_DG     SpliceAI_DL     SpliceAI_AG_Pos SpliceAI_AL_Pos SpliceAI_DG_Pos SpliceAI_DL_Pos ensembl_gene_id
A3GALT2 chr1_33777524_T_C       0.00    0.06    0.00    0.00    -1      -24     -3      -24     ENSG00000184389


################### Code Main points:

indx.keep<-c(which(colnames(SpliceAI) == "VAR"),which(colnames(SpliceAI) == "ensembl_gene_id"),which(colnames(SpliceAI) == "HGNC"))
SpliceAI.m<-melt(SpliceAI, id.vars=colnames(SpliceAI)[indx.keep], variable.name="variable", value.name="value")
SpliceAI.m$value<-as.numeric(SpliceAI.m$value)
indx.dep<-grep("_Pos",SpliceAI.m$variable)
SpliceAI.m_subset<-droplevels(SpliceAI.m[-indx.dep,]) # Melt SpliceAI data frame and select the variables SpliceAI_AG SpliceAI_AL SpliceAI_DG SpliceAI_DL

SpliceAI.m_subset_NO_NA<-SpliceAI.m_subset[!is.na(SpliceAI.m_subset$value),] # Exclude NA's

SpliceAI.m_subset_NO_NA.dt<-data.table(SpliceAI.m_subset_NO_NA, key=c("variable"))
SpliceAI.m_subset_NO_NA_SpliceAI_Type_parameters<-as.data.frame(SpliceAI.m_subset_NO_NA.dt[,.(mean_value=mean(value, na.rm =T),
                                                                  sd_value=sd(value, na.rm =T)),by=key(SpliceAI.m_subset_NO_NA.dt)], stringsAsFactors=F)
SpliceAI.m_subset_NO_NA<-merge(SpliceAI.m_subset_NO_NA,
                     SpliceAI.m_subset_NO_NA_SpliceAI_Type_parameters,
                     by=c("variable"))							       
SpliceAI.m_subset_NO_NA$value_Z_score<-(SpliceAI.m_subset_NO_NA$value-SpliceAI.m_subset_NO_NA$mean_value)/SpliceAI.m_subset_NO_NA$sd_value # Z score normalisation per variable


################### Output files:

$ Prepared_file_SpliceAI.rds
variable               VAR ensembl_gene_id    HGNC value  mean_value
SpliceAI_AG chr1_33777524_T_C ENSG00000184389 A3GALT2  0.00 0.001449955 
SpliceAI_AL chr1_33777524_T_C ENSG00000184389 A3GALT2  0.06 0.001397634
SpliceAI_DG chr1_33777524_T_C ENSG00000184389 A3GALT2  0.00 0.001639763
SpliceAI_DL chr1_33777524_T_C ENSG00000184389 A3GALT2  0.00 0.001467607
sd_value value_Z_score
0.01418904   -0.10218832
0.01487098    3.94072031
0.01925533   -0.08515894
0.01753616   -0.08369034

# 329_binder_of_scores_unranked_chromstates_v2.R

################### Input files: ALL_dB.tsv (annotation of GWAS blood traits)
################### Input files: chrmstates_graphs.csv
state_ordered,CellType_DEF,E_state,VAR
E1: Transcription Low signal H3K36me3,VB-Class-switched-B,E1,chr21_36280376_A_G
E1: Transcription Low signal H3K36me3,VB-Eff-mem-CD8,E1,chr21_36280376_A_G
E1: Transcription Low signal H3K36me3,TON-PC,E1,chr21_36280376_A_G
################### Input files: Weights_chromatin_states.tsv
Criteria# 1     Chromatin states        Weight
Active enhancer Active Enhancer High Signal H3K4me1 & H3K27Ac   1
Active promoter Active TSS High Signal H3K4me3 & H3K27Ac        1
Enhancer        Enhancer High Signal H3K4me1    0.1
Promoter        Active TSS High Signal H3K4me3 & H3K4me1        0.1
Repression states       Repressed Polycomb TSS High Signal H3K27me3 & H3K4me3 & H3K4me1 0.05
Repression states       Repressed Polycomb High signal H3K27me3 0.05
Repression states       Repressed Polycomb Low signal H3K27me3  0.05
Repression states       Heterochromatin High Signal H3K9me3     0.05
Transcription states    Transcription High signal H3K36me3      0.01
################### Input files: CellType_Trait_table_generation.txt
RelevantCellType        phenotype       label
CD34-negative.CD41-positive.CD42-positive.megakaryocyte.cell_cord.blood plt     C3
CD34-negative.CD41-positive.CD42-positive.megakaryocyte.cell_cord.blood mpv     C3
CD34-negative.CD41-positive.CD42-positive.megakaryocyte.cell_cord.blood pdw     C3
################### Input parameters: excluded_phenotypes=$(echo "wbc,eo_p,mono_p,neut_p,lymph_p,baso_p")
################### Input parameters: relevant_not_relevant_weights=$(echo "1,0.01")

################### Code Main points:

ALL_dB_subset_restricted<-unique(ALL_dB_subset[-which(ALL_dB_subset$phenotype%in%excluded_phenotypes),]) # Exclude phenotypes of white blood cells that are correlated
for(i in 1:length(phenotypes_array)) # Loop per phenotype
ALL_dB_double_subset_sel<-ALL_dB_double_subset[which(ALL_dB_double_subset$phenotype == phenotypes_array_sel),] # Select all variants associated to that phenotype
Trait_to_CT_table_sel<-Trait_to_CT_table[which(Trait_to_CT_table$phenotype == phenotypes_array_sel),] # Select all cell types relevant for that phenotype
chromstates_INITIAL_sel<-chromstates_INITIAL[which(chromstates_INITIAL$VAR%in%ALL_dB_double_subset_sel$VAR),] # Select all variants with chromatin state predictions associated to that phenotype

chromstates_INITIAL_sel$Tag[which(chromstates_INITIAL_sel$CellType_DEF%in%Trait_to_CT_table_sel$Cell_Type)]<-"Relevant"
chromstates_INITIAL_sel$Tag[-which(chromstates_INITIAL_sel$CellType_DEF%in%Trait_to_CT_table_sel$Cell_Type)]<-"Not_relevant" # Label every cell type prediction as Relevant or Not relevant

chromstates_INITIAL_sel<-merge(chromstates_INITIAL_sel,
                                   matrix_weighted_regulatory_states,
                                   by="state") # Merge the predictions with the table of weights per chromatin state

chromstates_INITIAL_sel.dt<-data.table(chromstates_INITIAL_sel, key=c("VAR","Tag"))
Aggregation_table<-as.data.frame(chromstates_INITIAL_sel.dt[,.(Aggregate_chromstates=sum(Weight),
                                                                  nCells=.N), by=key(chromstates_INITIAL_sel.dt)], stringsAsFactors=F) # Aggregate weights per variant and Tag ('Relevant', 'Not_relevant')

Aggregation_table$Multiplier[which(Aggregation_table$Tag == "Relevant")]<-relevant_not_relevant_weights[1]
Aggregation_table$Multiplier[which(Aggregation_table$Tag == "Not_relevant")]<-relevant_not_relevant_weights[2]
Aggregation_table$Aggregate_chromstates_multiplied<-Aggregation_table$Aggregate_chromstates*Aggregation_table$Multiplier # Multiply the aggregate weight of the chromatin states per variant and Tag by the weight of relevant (1) vs Not relevant cell type (0.01)

Aggregation_table$Aggregate_chromstates_normalised<-Aggregation_table$Aggregate_chromstates_multiplied/Aggregation_table$nCells # Normalise by the number of relvant and not relevant cells per variant for that phenotype

Aggregation_table.dt<-data.table(Aggregation_table, key=c("VAR"))
Aggregation_table_FINAL<-as.data.frame(Aggregation_table.dt[,.(Aggregate_chromstates_FINAL=sum(Aggregate_chromstates_normalised)), by=key(Aggregation_table.dt)], stringsAsFactors=F)
Aggregation_table_FINAL$phenotype<-phenotypes_array_sel # Finally aggregate the Relevant and No_relevant values to one value per variant and add the phenotype field

################### Output files: chromstates_GLOBAL_preranked.tsv
VAR     Aggregate_chromstates_FINAL     phenotype
chr10_104698523_G_A     0.0100148648648649      ret
chr10_104701039_T_TTATAA        0.0100124324324324      ret

# 337_Rank_chromstates_v2.R

################### Input files: chromstates_GLOBAL_preranked.tsv
VAR     Aggregate_chromstates_FINAL     phenotype
chr10_104698523_G_A     0.0100148648648649      ret
chr10_104701039_T_TTATAA        0.0100124324324324      ret
################### Input parameters: desiR_weights=$(echo "0.001,1.2,1")

################### Code main points First Function:

chromstates_pre_ranked.dt<-data.table(chromstates_pre_ranked, key=c("VAR"))
Aggregation_table<-as.data.frame(chromstates_pre_ranked.dt[,.(Total_Aggregate_chromstates=sum(Aggregate_chromstates_FINAL),
                                                                 nPhenotypes=.N), by=key(chromstates_pre_ranked.dt)], stringsAsFactors=F)
Aggregation_table$normalised_Total_Aggregate_chromstates<-Aggregation_table$Total_Aggregate_chromstates/Aggregation_table$nPhenotypes # Aggregate the Aggregate_chromstates_FINAL per variant (adding the number in different phenotypes) and normalise by the number of phenotypes per variant

Aggregation_table$normalised_Total_Aggregate_chromstates_component <- d.high(Aggregation_table$normalised_Total_Aggregate_chromstates, cut1=normalised_Total_Aggregate_chromstates_FINAL_LOW, cut2=normalised_Total_Aggregate_chromstates_FINAL_HIGH, scale=0.5)
Aggregation_table$Overall_weight <- d.overall(Aggregation_table$normalised_Total_Aggregate_chromstates_component,
                                          weights=c(Overall_normalised_Total_Aggregate_chromstates_FINAL)) # To get a continues value between 1 and 0 apply the d.high function of the desiR package. Below 0.001 everything is flattened to 0 and above 1.2 everything is flattened to 1.

################### Intermediate Output files: chromstates_GLOBAL_Ranked.tsv
VAR     Total_Aggregate_chromstates     nPhenotypes     normalised_Total_Aggregate_chromstates  normalised_Total_Aggregate_chromstates_component        Overall_weight
chr10_101244819_G_A     0.0004475       1       0.0004475       0       0

################### Code main points Second Function:

chromstates_ranked$mean_Overall_weight<-mean(chromstates_ranked$Overall_weight, na.rm =T)
chromstates_ranked$sd_Overall_weight<-sd(chromstates_ranked$Overall_weight, na.rm =T)
chromstates_ranked$Overall_weight_Z_score<-(chromstates_ranked$Overall_weight-chromstates_ranked$mean_Overall_weight)/chromstates_ranked$sd_Overall_weight # Z-score normalisation for all the variants

################### Output files: Prepared_file_chromstates.rds
VAR     value   value_Z_score   variable
chr12_111844956_C_T     0.915374733518353       4.58495209415341        Rank_chromstates
chr16_86016328_C_T      0.378500238084892       1.54630294787155        Rank_chromstates

# 329_binder_of_scores_unranked_but_thresholded_PCHiC_v2.R

################### Input files: ALL_dB.tsv (annotation of GWAS blood traits)
################### Input files: $ head PCHIC_ChicagoScore_graphs.csv. Chicago Score values for PCHiC filtered to aminimum of 5 Chicago Score per PIR that overlaps the variant.
VAR,HGNC,ensembl_gene_id,value,CellType_DEF
chr21_36789420_C_G,AP000330.8,ENSG00000234380,5.509692830975,aCD4
chr21_36789420_C_G,AP000330.8,ENSG00000234380,5.59914383391615,tCD4


################### Input files: PCHIC_part_II_Javierre_corresp_generation.txt
Trait   RelevantCellType
mono    Mac0
wbc     Mac0
mono_p  Mac0

################### Input parameters: excluded_phenotypes=$(echo "wbc,eo_p,mono_p,neut_p,lymph_p,baso_p")
################### Input parameters: relevant_not_relevant_weights=$(echo "1,0.01")

################### Code Main points:

ALL_dB_subset_restricted<-unique(ALL_dB_subset[-which(ALL_dB_subset$phenotype%in%excluded_phenotypes),]) # Exclude phenotypes of white blood cells that are correlated
for(i in 1:length(phenotypes_array)) # Loop per phenotype
ALL_dB_double_subset_sel<-ALL_dB_double_subset[which(ALL_dB_double_subset$phenotype == phenotypes_array_sel),] # Select all variants associated to that phenotype
Trait_to_CT_table_sel<-Trait_to_CT_table[which(Trait_to_CT_table$phenotype == phenotypes_array_sel),] # Select all cell types relevant for that phenotype

PCHiC_INITIAL_sel<-PCHiC_INITIAL[which(PCHiC_INITIAL$VAR%in%ALL_dB_double_subset_sel$VAR),] # Select the variants with PCHiC values associated to that phenotype

PCHiC_INITIAL_sel$Tag[which(PCHiC_INITIAL_sel$CellType_DEF%in%Trait_to_CT_table_sel$Cell_Type)]<-"Relevant"
PCHiC_INITIAL_sel$Tag[-which(PCHiC_INITIAL_sel$CellType_DEF%in%Trait_to_CT_table_sel$Cell_Type)]<-"Not_relevant" # Label every cell type prediction as Relevant or Not releva

PCHiC_INITIAL_sel.dt<-data.table(PCHiC_INITIAL_sel, key=c("VAR","Tag"))
Aggregation_table<-as.data.frame(PCHiC_INITIAL_sel.dt[,.(Aggregate_PCHiC=sum(value),
                                                                  nCells=.N), by=key(PCHiC_INITIAL_sel.dt)], stringsAsFactors=F) # Aggregate the Chicago Score value of all the interactions per variant and Tag

Aggregation_table$Multiplier[which(Aggregation_table$Tag == "Relevant")]<-relevant_not_relevant_weights[1]
Aggregation_table$Multiplier[which(Aggregation_table$Tag == "Not_relevant")]<-relevant_not_relevant_weights[2]
Aggregation_table$Aggregate_PCHiC_multiplied<-Aggregation_table$Aggregate_PCHiC*Aggregation_table$Multiplier # Multiply the aggregate Chicago scores per variant and Tag by the weight of relevant (1) vs Not relevant cell type (0.01

Aggregation_table$Aggregate_PCHiC_normalised<-Aggregation_table$Aggregate_PCHiC_multiplied/Aggregation_table$nCells # Normalise by the number of relvant and not relevant cells per variant for that phenotype

Aggregation_table_FINAL<-as.data.frame(Aggregation_table.dt[,.(Aggregate_PCHiC_FINAL=sum(Aggregate_PCHiC_normalised)), by=key(Aggregation_table.dt)], stringsAsFactors=F)
Aggregation_table_FINAL$phenotype<-phenotypes_array_sel # Finally aggregate the Relevant and No_relevant values to one value per variant and add the phenotype field

FINAL_df_NO_NA<-FINAL_df[!is.na(FINAL_df$Aggregate_PCHiC_FINAL),] # Exclude NA's

################### Output files: PCHiC_GLOBAL_preranked.tsv
VAR     Aggregate_PCHiC_FINAL   phenotype
chr10_104698523_G_A     0.0607758101027744      ret
chr10_104701039_T_TTATAA        0.0607758101027744      ret

# 336_Rank_PCHiC_v2.R

################### Input files: PCHiC_GLOBAL_preranked.tsv
VAR     Aggregate_PCHiC_FINAL   phenotype
chr10_104698523_G_A     0.0607758101027744      ret
chr10_104701039_T_TTATAA        0.0607758101027744      ret

################### Input parameters: desiR_weights=$(echo "0.1,10,1")

################### Code main points First Function:

PCHiC_pre_ranked.dt<-data.table(PCHiC_pre_ranked, key=c("VAR"))
Aggregation_table<-as.data.frame(PCHiC_pre_ranked.dt[,.(Total_Aggregate_PCHiC=sum(Aggregate_PCHiC_FINAL),
                                                                nPhenotypes=.N), by=key(PCHiC_pre_ranked.dt)], stringsAsFactors=F)
Aggregation_table$normalised_Total_Aggregate_PCHiC<-Aggregation_table$Total_Aggregate_PCHiC/Aggregation_table$nPhenotypes # Aggregate the Aggregate_PCHiC_FINAL per variant (adding the number in different phenotypes) and normalise by the number of phenotypes per variant

Aggregation_table$normalised_Total_Aggregate_PCHiC_component <- d.high(Aggregation_table$normalised_Total_Aggregate_PCHiC, cut1=normalised_Total_Aggregate_PCHiC_FINAL_LOW, cut2=normalised_Total_Aggregate_PCHiC_FINAL_HIGH, scale=0.5)
Aggregation_table$Overall_weight <- d.overall(Aggregation_table$normalised_Total_Aggregate_PCHiC_component,
                                          weights=c(Overall_normalised_Total_Aggregate_PCHiC_FINAL)) # To get a continues value between 1 and 0 apply the d.high function of the desiR package. Below 0.1 everything is flattened to 0 and above 10 everything is flattened to 1.

################### Intermediate Output files: PCHiC_GLOBAL_Ranked.tsv
VAR     Total_Aggregate_PCHiC   nPhenotypes     normalised_Total_Aggregate_PCHiC        normalised_Total_Aggregate_PCHiC_component      Overall_weight
chr10_101277639_TA_T    0.12268012554431        1       0.12268012554431        0.0478635745860307      0.0478635745860307

################### Code main points Second Function:

PCHiC_ranked$mean_Overall_weight<-mean(PCHiC_ranked$Overall_weight, na.rm =T)
PCHiC_ranked$sd_Overall_weight<-sd(PCHiC_ranked$Overall_weight, na.rm =T)
PCHiC_ranked$Overall_weight_Z_score<-(PCHiC_ranked$Overall_weight-PCHiC_ranked$mean_Overall_weight)/PCHiC_ranked$sd_Overall_weight # Z-score normalisation for all the variants

################### Output files: Prepared_file_PCHiC.rds
VAR     value   value_Z_score   variable
chr16_86016328_C_T      0.724737400667775       0.710225805378219       Rank_PCHiC
chr17_38764524_T_A      0.71661332492833        0.691028536045906       Rank_PCHiC



# 329_binder_of_scores_unranked_ATAC_v2.R

################### Input files: ALL_dB.tsv (annotation of GWAS blood traits)
################### Input files: ATAC_Seq_Ranked_graphs.csv. ATAC-seq values per cell type normalised by cell type.
VAR,BP_CELL_LABELS,LineageDEF,ATAC_Cell_Type,variable,value
chr21_36789420_C_G,CMP,gran_mono_lineage,CMP,ATAC_value,0.235284614543108

################### Input files: ATAC_scaled_trait_table_generation.txt (association lineage to blood traits)
Trait   Factor4
rbc     erythroid_lineage
mcv     erythroid_lineage

################### Input files: ATAC_scaled_Lineage_hierarchy_generation.txt (association between cell types and lineages)
CellType        Lineage Factor5
HSC     mega_lineage    HSC
MPP     mega_lineage    MPP

################### Input parameters: excluded_phenotypes=$(echo "wbc,eo_p,mono_p,neut_p,lymph_p,baso_p")
################### Input parameters: relevant_not_relevant_weights=$(echo "1,0.01")

################### Code Main points:

ALL_dB_subset_restricted<-unique(ALL_dB_subset[-which(ALL_dB_subset$phenotype%in%excluded_phenotypes),]) # Exclude phenotypes of white blood cells that are correlated
for(i in 1:length(phenotypes_array)) # Loop per phenotype
ALL_dB_double_subset_sel<-ALL_dB_double_subset[which(ALL_dB_double_subset$phenotype == phenotypes_array_sel),] # Select all variants associated to that phenotype

Trait_to_Lineage_table_sel<-Trait_to_Lineage_table[which(Trait_to_Lineage_table$Trait == phenotypes_array_sel),] 
Lineage_sel<-unique(Trait_to_Lineage_table_sel$Factor4) # Select the lineage associated with that phenotype
Lineage_to_Cell_table_sel<-Lineage_to_Cell_table[which(Lineage_to_Cell_table$Lineage%in%Trait_to_Lineage_table_sel$Factor4),] # Select all the cell types associated to that lineage

ATAC_INITIAL_subset_sel<-ATAC_INITIAL_subset[which(ATAC_INITIAL_subset$VAR%in%ALL_dB_double_subset_sel$VAR),] # Select all the variants associated to that phenotype

ATAC_INITIAL_subset_sel$Tag[which(ATAC_INITIAL_subset_sel$ATAC_Cell_Type%in%Lineage_to_Cell_table_sel$CellType)]<-"Relevant"
ATAC_INITIAL_subset_sel$Tag[-which(ATAC_INITIAL_subset_sel$ATAC_Cell_Type%in%Lineage_to_Cell_table_sel$CellType)]<-"Not_relevant" # Label every cell type prediction as Relevant or Not relevant

ATAC_INITIAL_subset_double_sel<-ATAC_INITIAL_subset_sel[which(ATAC_INITIAL_subset_sel$Tag == "Relevant"),] # Important,differently from the other rankings, we keep only the ATAC values of the relevant cell types for the lineage

indx.int<-c(which(colnames(ATAC_INITIAL_subset_double_sel) == "VAR"),
                which(colnames(ATAC_INITIAL_subset_double_sel) == "ATAC_Cell_Type"),
                which(colnames(ATAC_INITIAL_subset_double_sel) == "value"))
ATAC_INITIAL_subset_FINAL_sel<-unique(ATAC_INITIAL_subset_double_sel[,indx.int])
ATAC_INITIAL_subset_FINAL_sel$Lineage<-Lineage_sel
FINAL_df<-unique(rbind(FINAL_df,ATAC_INITIAL_subset_FINAL_sel)) # Important, we add the lineage field not the phnotype field and make a unique function. This means that one variant associated to two traits of the same lineage will be collapsed to just one set of unique values of ATAC seq of all the cell types relevant for that lineage

FINAL_df_NO_NA<-FINAL_df[!is.na(FINAL_df$value),] # Exclude NA's

################### Output files: ATAC_GLOBAL_preranked.tsv

chr1_202129205_G_A      CMP     0.305548373494111       erythroid_lineage
chr1_202129205_G_A      HSC     0.280952138034825       erythroid_lineage
chr1_202129205_G_A      MPP     0.293918507628686       erythroid_lineage
chr1_202129205_G_A      Ery     0.28599250535214        erythroid_lineage
chr1_202129205_G_A      MEP     0.329119986078948       erythroid_lineage
chr1_202129205_G_A      CMP     0.305548373494111       mega_lineage
chr1_202129205_G_A      HSC     0.280952138034825       mega_lineage
chr1_202129205_G_A      MPP     0.293918507628686       mega_lineage
chr1_202129205_G_A      MEP     0.329119986078948       mega_lineage
chr1_202129205_G_A      Mega    0.394347811595599       mega_lineage


# 330_Rank_ATAC_v2.R

################### Input files: ATAC_GLOBAL_preranked.tsv

chr1_202129205_G_A      CMP     0.305548373494111       erythroid_lineage
chr1_202129205_G_A      HSC     0.280952138034825       erythroid_lineage
chr1_202129205_G_A      MPP     0.293918507628686       erythroid_lineage
chr1_202129205_G_A      Ery     0.28599250535214        erythroid_lineage
chr1_202129205_G_A      MEP     0.329119986078948       erythroid_lineage
chr1_202129205_G_A      CMP     0.305548373494111       mega_lineage
chr1_202129205_G_A      HSC     0.280952138034825       mega_lineage
chr1_202129205_G_A      MPP     0.293918507628686       mega_lineage
chr1_202129205_G_A      MEP     0.329119986078948       mega_lineage
chr1_202129205_G_A      Mega    0.394347811595599       mega_lineage


################### Input parameters: Open_in_CT_threshold=$(echo "0.1")
################### Input parameters: desiR_weights=$(echo "0.5,5,0.05,0.8,1,3")

################### Code main points First Function:

ATAC_pre_ranked$Lineage[which(ATAC_pre_ranked$Lineage == 'lymph_lineage.CD4')]<-"lymph_lineage"
ATAC_pre_ranked$Lineage[which(ATAC_pre_ranked$Lineage == 'lymph_lineage.CD8')]<-"lymph_lineage"
ATAC_pre_ranked$Lineage[which(ATAC_pre_ranked$Lineage == 'lymph_lineage.B')]<-"lymph_lineage"
ATAC_pre_ranked$Lineage[which(ATAC_pre_ranked$Lineage == 'lymph_lineage.NK')]<-"lymph_lineage" # Convert all lymph sublineages to lymph_lineage

ATAC_pre_ranked_unified_lymph_Thresholded<-ATAC_pre_ranked_unified_lymph[which(ATAC_pre_ranked_unified_lymph$value >= Open_in_CT_threshold),]
ATAC_pre_ranked_unified_lymph_Thresholded.dt<-data.table(ATAC_pre_ranked_unified_lymph_Thresholded, key=c("VAR","Lineage"))
Freq_table<-as.data.frame(ATAC_pre_ranked_unified_lymph_Thresholded.dt[,.(nCells=.N), by=key(ATAC_pre_ranked_unified_lymph_Thresholded.dt)], stringsAsFactors=F) # Use a Threshold of 0.1 of scaled signal to filter all the cells per variant and lineage above that value and count how many they are


ATAC_pre_ranked_unified_lymph.dt<-data.table(ATAC_pre_ranked_unified_lymph, key=c("VAR","Lineage"))
MAX_table<-as.data.frame(ATAC_pre_ranked_unified_lymph.dt[,.SD[which.max(value)], by=key(ATAC_pre_ranked_unified_lymph.dt)], stringsAsFactors=F) # Separately get the maximum ATAC seq signal per variant and phenotype

Merge_table<-merge(Freq_table,
                     MAX_table,
                     by=c("VAR","Lineage"),
                    all=T)
Merge_table$nCells[is.na(Merge_table$nCells)]<-0 # Merge both tables and make 0 the values for variants that don;t have any cell type above the threshold for a given lineage

Merge_table$nCells_weight <- d.high(Merge_table$nCells, cut1=nCells_LOW, cut2=nCells_HIGH, scale=0.5) # Calculate a continuous value for the number of cells per lineage above the threshold. Values above 5 are flattened to 1
Merge_table$ATAC_value_weight <- d.high(Merge_table$value, cut1=ATAC_value_LOW, cut2=ATAC_value_HIGH, scale=0.5) # Calculate a continuous value for the maximum ATAC signal per variant and lineage. Values below 0.05 are flattened to 0 and above 0.8 are flattened to 1.

Merge_table$Overall_weight <- d.overall(Merge_table$nCells_weight, Merge_table$ATAC_value_weight,
                                          weights=c(Overall_nCells,Overall_ATAC_value)) # With the previous two calculate a composite measure giving a relative weight of 1 to the Overall_nCells and of 3 to the Overall_ATAC_value

################### Intermediate Output files: ATAC_GLOBAL_Ranked.tsv
VAR     Lineage nCells  ATAC_Cell_Type  value   nCells_weight   ATAC_value_weight       Overall_weight
chr10_101245230_T_C     mega_lineage    1       HSC     0.100621412475289       0.333333333333333       0.259798415379537       0.276501224964101
chr10_101248979_A_G     erythroid_lineage       5       HSC     0.536134835602244       1       0.805096131404811       0.849935178410447
chr10_101248979_A_G     mega_lineage    4       HSC     0.536134835602244       0.881917103688197       0.805096131404811       0.823650082327669


################### Code main points Second Function:

ATAC_ranked.dt<-data.table(ATAC_ranked, key=c("Lineage"))
ATAC_ranked_Lineage_parameters<-as.data.frame(ATAC_ranked.dt[,.(mean_Overall_weight=mean(Overall_weight, na.rm =T),
                                                                  sd_Overall_weight=sd(Overall_weight, na.rm =T)),
                                                               by=key(ATAC_ranked.dt)], stringsAsFactors=F)
ATAC_ranked<-merge(ATAC_ranked,
                     ATAC_ranked_Lineage_parameters,
                     by=c("Lineage"))
ATAC_ranked$Overall_weight_Z_score<-(ATAC_ranked$Overall_weight-ATAC_ranked$mean_Overall_weight)/ATAC_ranked$sd_Overall_weight # Z-score normalisation per Lineage


colnames(ATAC_ranked_subset)[which(colnames(ATAC_ranked_subset) == "Lineage")]<-"variable"
colnames(ATAC_ranked_subset)[which(colnames(ATAC_ranked_subset) == "Overall_weight")]<-"value"
colnames(ATAC_ranked_subset)[which(colnames(ATAC_ranked_subset) == "Overall_weight_Z_score")]<-"value_Z_score"
ATAC_ranked_subset$variable<-paste('Rank_ATAC_',ATAC_ranked_subset$variable, sep='') # Subset columns and create the variable Rank_ATAC_ per lineage


################### Output files: Prepared_file_ATAC.rds

VAR     variable        value   value_Z_score
chr1_202129205_G_A      Rank_ATAC_erythroid_lineage     0.690277355735851       0.795936299393506
chr1_202129205_G_A      Rank_ATAC_mega_lineage  0.74683734086938        0.999126381544563

# 331_Rank_multi_lineage_ATAC.R

################### Input files: ATAC_GLOBAL_Ranked.tsv
VAR     Lineage nCells  ATAC_Cell_Type  value   nCells_weight   ATAC_value_weight       Overall_weight
chr10_101245230_T_C     mega_lineage    1       HSC     0.100621412475289       0.333333333333333       0.259798415379537       0.276501224964101
chr10_101248979_A_G     erythroid_lineage       5       HSC     0.536134835602244       1       0.805096131404811       0.849935178410447
chr10_101248979_A_G     mega_lineage    4       HSC     0.536134835602244       0.881917103688197       0.805096131404811       0.823650082327669
################### Input parameters: desiR_weights=$(echo "0.9,2.5,0.1,0.7,3,1")

################### Code main points First Function:

multi_ATAC_ranked.dt<-data.table(multi_ATAC_ranked, key="VAR")
Freq_table<-as.data.frame(multi_ATAC_ranked.dt[,.(nLineages=.N), by=key(multi_ATAC_ranked.dt)], stringsAsFactors=F) # Get a table of how many lineages are associated to every variant with ATAC seq signal

multi_ATAC_ranked.dt<-data.table(multi_ATAC_ranked, key="VAR")
mean_table<-as.data.frame(multi_ATAC_ranked.dt[,.(mean_Rank_ATAC=mean(Overall_weight)), by=key(multi_ATAC_ranked.dt)], stringsAsFactors=F) # Get a table of what's the mean Overall_weight as calculated in the script 330_Rank_ATAC_v2.R. It is a mean of the ATAC ranking across lineages for the same variant.

Merge_table<-merge(Freq_table,
                     mean_table,
                     by="VAR",
                     all=T)# Merge both tables

Merge_table$nLineages_weight <- d.high(Merge_table$nLineages, cut1=nLineages_LOW, cut2=nLineages_HIGH, scale=0.5) # Give a continuous value between 0 to 1 to the number of lineages per variant. Values above 2 are flattened to 1.

Merge_table$mean_Rank_ATAC_weight <- d.high(Merge_table$mean_Rank_ATAC, cut1=mean_Rank_ATAC_LOW, cut2=mean_Rank_ATAC_HIGH, scale=0.5) # Give a continuous value between 0 to 1 to the mean_Rank_ATAC across lineages for the same variant. Values below 0.1 are flattened to 0 and above 0.7 are flattened to 1.

Merge_table$Overall_weight <- d.overall(Merge_table$nLineages_weight, Merge_table$mean_Rank_ATAC_weight,
                                          weights=c(Overall_nLineages,Overall_mean_Rank_ATAC))  # With the previous two calculate a composite measure giving a relative weight of 3 to the nLineages_weight and of 1 to the Overall_mean_Rank_ATAC

################### Intermediate Output files: multi_ATAC_GLOBAL_Ranked.tsv

VAR     nLineages       mean_Rank_ATAC  nLineages_weight        mean_Rank_ATAC_weight   Overall_weight
chr10_101245230_T_C     1       0.276501224964101       0.25    0.542373218617496       0.303409773098511
chr10_101248979_A_G     2       0.836792630369058       0.82915619758885        1       0.868914937906556



################### Code main points Second Function:

multi_ATAC_ranked$mean_Overall_weight<-mean(multi_ATAC_ranked$Overall_weight, na.rm =T)
multi_ATAC_ranked$sd_Overall_weight<-sd(multi_ATAC_ranked$Overall_weight, na.rm =T)
multi_ATAC_ranked$Overall_weight_Z_score<-(multi_ATAC_ranked$Overall_weight-multi_ATAC_ranked$mean_Overall_weight)/multi_ATAC_ranked$sd_Overall_weight # Z-score normalisation

indx.int<-c(which(colnames(multi_ATAC_ranked) == "VAR"),which(colnames(multi_ATAC_ranked) == "Overall_weight"),which(colnames(multi_ATAC_ranked) == "Overall_weight_Z_score"))
multi_ATAC_ranked_subset<-unique(multi_ATAC_ranked[,indx.int])
colnames(multi_ATAC_ranked_subset)[which(colnames(multi_ATAC_ranked_subset) == "Overall_weight")]<-"value"
colnames(multi_ATAC_ranked_subset)[which(colnames(multi_ATAC_ranked_subset) == "Overall_weight_Z_score")]<-"value_Z_score"
multi_ATAC_ranked_subset$variable<-"multi_lineage_ATAC" # Subset columns and create the variable multi_lineage_ATAC

################### Output files: Prepared_file_multi_lineage_ATAC.rds

VAR     value   value_Z_score   variable
chr1_202129205_G_A      0.868914937906556       2.33058889131214        multi_lineage_ATAC
chr12_111844956_C_T     0.353553390593274       -0.114069886812755      multi_lineage_ATAC


# 338_binder_of_scores_Gene_EXP.R

################### Input files: ALL_dB.tsv (annotation of GWAS blood traits)
################### Input files: GENE_EXP_Harmonization.tsv Mean gene expression values for all the BluePrint cell types. See Table Sn for the actual samples included.

d-B     sd_VB-Class-switched-B  cv_VB-Class-switched-B  Rank_GENE_VB-Class-switched-B   Rank_discretized_VB-Class-switched-B    mean_VB-mem-B   sd_VB-mem-B     cv_VB-mem-B     Rank_GENE_VB-mem-B      Rank_discretized_VB-mem-B       mean_CB-CD4     sd_CB-CD4       cv_CB-CD4       Rank_GENE_CB-CD4        Rank_discretized_CB-CD4      mean_VB-CD4     sd_VB-CD4       cv_VB-CD4       Rank_GENE_VB-CD4        Rank_discretized_VB-CD4 mean_VB-Mem-CD4 sd_VB-Mem-CD4   cv_VB-Mem-CD4   Rank_GENE_VB-Mem-CD4    Rank_discretized_VB-Mem-CD4     mean_VB-Eff-mem-CD4     sd_VB-Eff-mem-CD4       cv_VB-Eff-mem-CD4       Rank_GENE_VB-Eff-mem-CD4    Rank_discretized_VB-Eff-mem-CD4  mean_VB-T-Reg   sd_VB-T-Reg     cv_VB-T-Reg     Rank_GENE_VB-T-Reg      Rank_discretized_VB-T-Reg       mean_CB-CD8     sd_CB-CD8       cv_CB-CD8       Rank_GENE_CB-CD8        Rank_discretized_CB-CD8 mean_VB-Mem-CD8 sd_VB-Mem-CD8   cv_VB-Mem-CD8   Rank_GENE_VB-Mem-CD8    Rank_discretized_VB-Mem-CD8  mean_VB-Eff-mem-CD8     sd_VB-Eff-mem-CD8       cv_VB-Eff-mem-CD8       Rank_GENE_VB-Eff-mem-CD8        Rank_discretized_VB-Eff-mem-CD8 mean_CB-NK      sd_CB-NK        cv_CB-NK        Rank_GENE_CB-NK Rank_discretized_CB-NK  mean_VB-NK      sd_VB-NK        cv_VB-NK        Rank_GENE_VB-NK Rank_discretized_VB-NK       mean_CB-MONO    sd_CB-MONO      cv_CB-MONO      Rank_GENE_CB-MONO       Rank_discretized_CB-MONO        mean_VB-MONO    sd_VB-MONO      cv_VB-MONO      Rank_GENE_VB-MONO       Rank_discretized_VB-MONO        mean_CB-MAC     sd_CB-MAC       cv_CB-MAC       Rank_GENE_CB-MAC        Rank_discretized_CB-MAC      mean_VB-MAC     sd_VB-MAC       cv_VB-MAC       Rank_GENE_VB-MAC        Rank_discretized_VB-MAC mean_VB-Inf-MAC sd_VB-Inf-MAC   cv_VB-Inf-MAC   Rank_GENE_VB-Inf-MAC    Rank_discretized_VB-Inf-MAC     mean_VB-Alt-MAC sd_VB-Alt-MAC   cv_VB-Alt-MAC   Rank_GENE_VB-Alt-MAC    Rank_discretized_VB-Alt-MAC     mean_CB-Alt-MAC      sd_CB-Alt-MAC   cv_CB-Alt-MAC   Rank_GENE_CB-Alt-MAC    Rank_discretized_CB-Alt-MAC     mean_CB-NEUT    sd_CB-NEUT      cv_CB-NEUT      Rank_GENE_CB-NEUT       Rank_discretized_CB-NEUT        mean_VB-NEUT    sd_VB-NEUT      cv_VB-NEUT      Rank_GENE_VB-NEUT       Rank_discretized_VB-NEUT        mean_HSC     sd_HSC  cv_HSC  Rank_GENE_HSC   Rank_discretized_HSC    mean_HMP        sd_HMP  cv_HMP  Rank_GENE_HMP   Rank_discretized_HMP    mean_MEP        sd_MEP  cv_MEP  Rank_GENE_MEP   Rank_discretized_MEP    mean_CMP        sd_CMP  cv_CMP  Rank_GENE_CMP   Rank_discretized_CMP    mean_GMP        sd_GMP  cv_GMP  Rank_GENE_GMP        Rank_discretized_GMP    mean_CLP        sd_CLP  cv_CLP  Rank_GENE_CLP   Rank_discretized_CLP    mean_CB-DC      sd_CB-DC        cv_CB-DC        Rank_GENE_CB-DC Rank_discretized_CB-DC  mean_CB-Prolif-Endo     sd_CB-Prolif-Endo       cv_CB-Prolif-Endo       Rank_GENE_CB-Prolif-Endo        Rank_discretized_CB-Prolif-Endo      mean_CB-Resting-Endo    sd_CB-Resting-Endo      cv_CB-Resting-Endo      Rank_GENE_CB-Resting-Endo       Rank_discretized_CB-Resting-Endo        mean_VB-stim-0d-MONO    sd_VB-stim-0d-MONO      cv_VB-stim-0d-MONO      Rank_GENE_VB-stim-0d-MONO       Rank_discretized_VB-stim-0d-MONO        mean_VB-stim-none-MONO       sd_VB-stim-none-MONO    cv_VB-stim-none-MONO    Rank_GENE_VB-stim-none-MONO     Rank_discretized_VB-stim-none-MONO      mean_VB-stim-6d-untreated-MAC   sd_VB-stim-6d-untreated-MAC     cv_VB-stim-6d-untreated-MAC     Rank_GENE_VB-stim-6d-untreated-MAC      Rank_discretized_VB-stim-6d-untreated-MAC    mean_VB-stim-6d-Glucan-MAC      sd_VB-stim-6d-Glucan-MAC        cv_VB-stim-6d-Glucan-MAC        Rank_GENE_VB-stim-6d-Glucan-MAC Rank_discretized_VB-stim-6d-Glucan-MAC  mean_VB-stim-6d-LPS-MAC sd_VB-stim-6d-LPS-MAC   cv_VB-stim-6d-LPS-MAC   Rank_GENE_VB-stim-6d-LPS-MAC    Rank_discretized_VB-stim-6d-LPS-MAC
ENSG00000000003 TSPAN6  1.90653670949975        0.642239497896798       33.6861857784691        0.336711933472643       [1.64,3.04)     6.19232744667757        3.34852137865528        54.0753280166393        0.581400938121943       [3.36,8.61)     0       0       0       0       0       0       0       0       0   00       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0   00       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0   00       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0   00       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       11.956365       0       0       0.617106886477477       [8.75,21.1)     17.902283       0       0       0.669598572825032       [7.8,18)        0       0       0       0       0   00       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0


################### Input files: Correspondence_phenotype_TOME.txt
phenotype       TOME_Cell_Type  Lineage
rbc     EB      erythroid_lineage
mcv     EB      erythroid_lineage
hct     EB      erythroid_lineage
mch     EB      erythroid_lineage

################### Input files: PCHIC_ChicagoScore_graphs.csv
VAR,HGNC,ensembl_gene_id,value,CellType_DEF
chr21_36789420_C_G,AP000330.8,ENSG00000234380,5.509692830975,aCD4
chr21_36789420_C_G,AP000330.8,ENSG00000234380,5.59914383391615,tCD4

################### Input files: VEP_consequence_graphs.csv VEP consequences per gene per variant

transcript_id,ensembl_gene_id,HGNC,VEP_DEF_LABELS,VAR
ENST00000300305,ENSG00000159216,RUNX1,INTRON,chr21_36280376_A_G
ENST00000416754,ENSG00000159216,RUNX1,INTRON,chr21_36280376_A_G

################### Input parameters: excluded_phenotypes=$(echo "wbc,eo_p,mono_p,neut_p,lymph_p,baso_p")


################### Code Main points:

ALL_dB_subset_restricted<-unique(ALL_dB_subset[-which(ALL_dB_subset$phenotype%in%excluded_phenotypes),]) # Exclude phenotypes of white blood cells that are correlated

VEP_CSQ_subset_with_genes<-unique(VEP_CSQ_subset[which(VEP_CSQ_subset$ensembl_gene_id != "DUMMY"),]) # Exclude VEP consequences that are not linked to genes

VEP_and_PCHiC<-unique(rbind(VEP_CSQ_subset_with_genes,PCHiC_subset)) # Get a data frame with the genes linked to a variant by VEP consequences and PCHiC interactions

for(i in 1:length(phenotypes_array)) # Loop per phenotype
ALL_dB_double_subset_sel<-ALL_dB_double_subset[which(ALL_dB_double_subset$phenotype == phenotypes_array_sel),] # Select all variants associated to that phenotype


VEP_and_PCHiC_sel<-VEP_and_PCHiC[which(VEP_and_PCHiC$VAR%in%ALL_dB_subset_restricted_sel$VAR),] # Select all variants associated to that phenotype with VEP and/or PCHiC linked genes

GENE_EXP_sel<-GENE_EXP[which(GENE_EXP$ensembl_gene_id%in%VEP_and_PCHiC_sel$ensembl_gene_id),] 
indx.int<-c(which(colnames(GENE_EXP_sel) == "ensembl_gene_id"),grep("mean_",colnames(GENE_EXP_sel)))
GENE_EXP_sel_subset<-unique(GENE_EXP_sel[,indx.int]) # Select gene expression data for those genes, only the mean expression parameter per cell type

GENE_EXP_sel_subset.m$Cell_Type<-gsub("mean_","",GENE_EXP_sel_subset.m$Cell_Type)
GENE_EXP_sel_subset.m$Cell_Type<-gsub("\\.","-",GENE_EXP_sel_subset.m$Cell_Type) # Obtain the Cell Type field

TOME_correspondence_sel<-TOME_correspondence[which(TOME_correspondence$phenotype == phenotypes_array_sel),] # Obtain the cell types relevant for the phenotype

GENE_EXP_sel_subset.m$Tag[which(GENE_EXP_sel_subset.m$Cell_Type%in%TOME_correspondence$TOME_Cell_Type)]<-"Not_relevant"
GENE_EXP_sel_subset.m$Tag[which(GENE_EXP_sel_subset.m$Cell_Type%in%TOME_correspondence_sel$TOME_Cell_Type)]<-"Relevant" # Label the expression data in the cell type as Relevant or Not_relevant

GENE_EXP_sel_subset.m_NO_NA$phenotype<-phenotypes_array_sel # Add the phenotype field

GENE_EXP_sel_subset.m_NO_NA_wide<-as.data.frame(pivot_wider(GENE_EXP_sel_subset.m_NO_NA, id_cols=c("VAR","ensembl_gene_id","HGNC","phenotype"),
                                                  names_from=Cell_Type,
                                                  values_from=c("mean_GENE_EXP","Tag")), stringsAsFactors=F) # convert to the wide format to make the intermediate output smaller

################### Output files: GENE_EXP_GLOBAL.tsv

VAR     ensembl_gene_id HGNC    phenotype       mean_GENE_EXP_MK        mean_GENE_EXP_VB-CD4    mean_GENE_EXP_CB-Alt-MAC        mean_GENE_EXP_CMP       mean_GENE_EXP_VB-Mem-CD4        mean_GENE_EXP_HMP       mean_GENE_EXP_VB-B      mean_GENE_EXP_CLP       mean_GENE_EXP_EB        mean_GENE_EXP_VB-MAC    mean_GENE_EXP_VB-NK     mean_GENE_EXP_CB-NEUT   mean_GENE_EXP_GMP       mean_GENE_EXP_HSC       mean_GENE_EXP_CB-CD4    mean_GENE_EXP_VB-Inf-MAC        mean_GENE_EXP_VB-Alt-MAC        mean_GENE_EXP_VB-Eff-mem-CD8    mean_GENE_EXP_VB-Eff-mem-CD4    mean_GENE_EXP_CB-MAC    mean_GENE_EXP_CB-B      mean_GENE_EXP_MEP       mean_GENE_EXP_VB-NEUT   mean_GENE_EXP_VB-MONO   mean_GENE_EXP_CB-MONO   mean_GENE_EXP_VB-Class-switched-B       mean_GENE_EXP_VB-Mem-CD8        mean_GENE_EXP_VB-T-Reg  mean_GENE_EXP_VB-mem-B  mean_GENE_EXP_CB-NK     mean_GENE_EXP_CB-CD8    Tag_MK  Tag_VB-CD4      Tag_CB-Alt-MAC  Tag_CMP Tag_VB-Mem-CD4  Tag_HMP Tag_VB-B
        Tag_CLP Tag_EB  Tag_VB-MAC      Tag_VB-NK       Tag_CB-NEUT     Tag_GMP Tag_HSC Tag_CB-CD4      Tag_VB-Inf-MAC  Tag_VB-Alt-MAC  Tag_VB-Eff-mem-CD8      Tag_VB-Eff-mem-CD4      Tag_CB-MAC      Tag_CB-B        Tag_MEP Tag_VB-NEUT     Tag_VB-MONO     Tag_CB-MONO     Tag_VB-Class-switched-B Tag_VB-Mem-CD8  Tag_VB-T-Reg    Tag_VB-mem-B    Tag_CB-NK       Tag_CB-CD8
chr3_50332697_G_A       ENSG00000001617 SEMA3F  baso    0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       2.33386433333333        0       0       0       0       0       0
       Not_relevant    Not_relevant    Not_relevant    Relevant        Not_relevant    Relevant        Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Relevant        Relevant        Relevant        Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Relevant        Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant

# 339_Gene_EXP.R

################### Input files: GENE_EXP_GLOBAL.tsv

VAR     ensembl_gene_id HGNC    phenotype       mean_GENE_EXP_MK        mean_GENE_EXP_VB-CD4    mean_GENE_EXP_CB-Alt-MAC        mean_GENE_EXP_CMP       mean_GENE_EXP_VB-Mem-CD4        mean_GENE_EXP_HMP       mean_GENE_EXP_VB-B      mean_GENE_EXP_CLP       mean_GENE_EXP_EB        mean_GENE_EXP_VB-MAC    mean_GENE_EXP_VB-NK     mean_GENE_EXP_CB-NEUT   mean_GENE_EXP_GMP       mean_GENE_EXP_HSC       mean_GENE_EXP_CB-CD4    mean_GENE_EXP_VB-Inf-MAC        mean_GENE_EXP_VB-Alt-MAC        mean_GENE_EXP_VB-Eff-mem-CD8    mean_GENE_EXP_VB-Eff-mem-CD4    mean_GENE_EXP_CB-MAC    mean_GENE_EXP_CB-B      mean_GENE_EXP_MEP       mean_GENE_EXP_VB-NEUT   mean_GENE_EXP_VB-MONO   mean_GENE_EXP_CB-MONO   mean_GENE_EXP_VB-Class-switched-B       mean_GENE_EXP_VB-Mem-CD8        mean_GENE_EXP_VB-T-Reg  mean_GENE_EXP_VB-mem-B  mean_GENE_EXP_CB-NK     mean_GENE_EXP_CB-CD8    Tag_MK  Tag_VB-CD4      Tag_CB-Alt-MAC  Tag_CMP Tag_VB-Mem-CD4  Tag_HMP Tag_VB-B
        Tag_CLP Tag_EB  Tag_VB-MAC      Tag_VB-NK       Tag_CB-NEUT     Tag_GMP Tag_HSC Tag_CB-CD4      Tag_VB-Inf-MAC  Tag_VB-Alt-MAC  Tag_VB-Eff-mem-CD8      Tag_VB-Eff-mem-CD4      Tag_CB-MAC      Tag_CB-B        Tag_MEP Tag_VB-NEUT     Tag_VB-MONO     Tag_CB-MONO     Tag_VB-Class-switched-B Tag_VB-Mem-CD8  Tag_VB-T-Reg    Tag_VB-mem-B    Tag_CB-NK       Tag_CB-CD8
chr3_50332697_G_A       ENSG00000001617 SEMA3F  baso    0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       2.33386433333333        0       0       0       0       0       0
       Not_relevant    Not_relevant    Not_relevant    Relevant        Not_relevant    Relevant        Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Relevant        Relevant        Relevant        Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Relevant        Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant    Not_relevant


################### Input parameters: relevant_not_relevant_weights=$(echo "1,0.01")
################### Input parameters: desiR_weights=$(echo "0.1,20,1")

################### Code main points First Function:

GENE_EXP_pre_ranked.m<-melt(GENE_EXP_pre_ranked, id.vars=c("VAR","ensembl_gene_id","HGNC","phenotype"), value.name="value", variable.name="variable") # Go from the wide matrix to the long matrix per VAR, gene and phenotype

GENE_EXP_pre_ranked.m$variable<-as.character(GENE_EXP_pre_ranked.m$variable)
GENE_EXP_pre_ranked.m$Cell_Type<-GENE_EXP_pre_ranked.m$variable
GENE_EXP_pre_ranked.m$Cell_Type<-gsub("mean_GENE_EXP_|Tag_","", GENE_EXP_pre_ranked.m$Cell_Type) # Obtain the Cell_Type field from the variable field

GENE_EXP_pre_ranked.m$variable[grep("mean_GENE_EXP_", GENE_EXP_pre_ranked.m$variable)]<-"mean_GENE_EXP"
GENE_EXP_pre_ranked.m$variable[grep("Tag_", GENE_EXP_pre_ranked.m$variable)]<-"Tag" # Obtain the variable field with mean and Tag

GENE_EXP_pre_ranked.m_wide<-as.data.frame(pivot_wider(GENE_EXP_pre_ranked.m,
                                                        id_cols=c("VAR","HGNC","ensembl_gene_id","phenotype","Cell_Type"),
                                                        names_from=variable,
                                                        values_from=value), stringsAsFactors=F) # Pivot to wider to have the mean_Gene_EXP and Tag in different columns per variant, gene, phenotype and cell type

GENE_EXP_pre_ranked.m_wide.dt<-data.table(GENE_EXP_pre_ranked.m_wide, key=c("VAR","HGNC","ensembl_gene_id","Tag"))
Aggregation_table<-as.data.frame(GENE_EXP_pre_ranked.m_wide.dt[,.(Aggregate_GENE_EXP=sum(mean_GENE_EXP),
                                                          nCells=.N), by=key(GENE_EXP_pre_ranked.m_wide.dt)], stringsAsFactors=F) # Aggregate mean gene expression per variant, gene and Tag. Don't use phenotypes as they will be proportional to the number of cell types in the normalisation

Aggregation_table$Multiplier[which(Aggregation_table$Tag == "Relevant")]<-relevant_not_relevant_weights[1]
Aggregation_table$Multiplier[which(Aggregation_table$Tag == "Not_relevant")]<-relevant_not_relevant_weights[2]
Aggregation_table$Aggregate_GENE_EXP_multiplied<-Aggregation_table$Aggregate_GENE_EXP*Aggregation_table$Multiplier # Introduce the Relevant / Not Relevant weights 1/0.01

Aggregation_table$Aggregate_GENE_EXP_normalised<-Aggregation_table$Aggregate_GENE_EXP_multiplied/Aggregation_table$nCells # Normalise by the number of cell types relevant and not relevant
Aggregation_table_FINAL<-as.data.frame(Aggregation_table.dt[,.(Aggregate_GENE_EXP_FINAL=sum(Aggregate_GENE_EXP_normalised)), by=key(Aggregation_table.dt)], stringsAsFactors=F) # Get a final number per variant and gene

Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL_component <- d.high(Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL, cut1=Aggregate_GENE_EXP_FINAL_LOW, cut2=Aggregate_GENE_EXP_FINAL_HIGH, scale=0.5) 
Aggregation_table_FINAL$Overall_weight <- d.overall(Aggregation_table_FINAL$Aggregate_GENE_EXP_FINAL_component,
                                          weights=c(Overall_Aggregate_GENE_EXP_FINAL)) # Obtain a continuous value from 0 to 1 for the Aggregate_GENE_EXP_FINAL. Values below 0.1 are flattened to 0 and above 20 are flattened to 1

################### Intermediate Output files: GENE_EXP_GLOBAL_Ranked.tsv
VAR     ensembl_gene_id HGNC    Aggregate_GENE_EXP_FINAL        Aggregate_GENE_EXP_FINAL_component      Overall_weight
chr10_101277639_TA_T    ENSG00000120053 GOT1    6.2849340093879 0.55749502597722        0.55749502597722


################### Code main points Second Function:

GENE_EXP_ranked$mean_Overall_weight<-mean(GENE_EXP_ranked$Overall_weight, na.rm =T)
GENE_EXP_ranked$sd_Overall_weight<-sd(GENE_EXP_ranked$Overall_weight, na.rm =T)
GENE_EXP_ranked$Overall_weight_Z_score<-(GENE_EXP_ranked$Overall_weight-GENE_EXP_ranked$mean_Overall_weight)/GENE_EXP_ranked$sd_Overall_weight # Z spcre normalisation of the Overall_weight across all genes and variants

################### Output files: Prepared_file_GENE_EXP.rds
VAR     ensembl_gene_id HGNC    value   value_Z_score   variable
chr12_111844956_C_T     ENSG00000111252 SH2B3   1       2.22500841360019        Rank_GENE_EXP
chr12_111844956_C_T     ENSG00000257595 RP3-473L9.4     0.187394957124891       -0.354927879046733      Rank_GENE_EXP

# 334_binder_of_scores_COGS.R

################### Input files: ALL_dB.tsv (annotation of GWAS blood traits)
################### Input files: COGS scores calculated in 7 cell types of PCHiC interactions (Ery,MK,Mono,Neut,tCD4,tCD8 and tB)
ensg    name    cogs
ENSG00000116285 ERRFI1  1.09892749422524e-06
ENSG00000116288 PARK7   1.40998324127395e-14
ENSG00000049249 TNFRSF9 0
ENSG00000131686 CA6     0

################### Input files: PCHIC_ChicagoScore_graphs.csv
VAR,HGNC,ensembl_gene_id,value,CellType_DEF
chr21_36789420_C_G,AP000330.8,ENSG00000234380,5.509692830975,aCD4
chr21_36789420_C_G,AP000330.8,ENSG00000234380,5.59914383391615,tCD4

################### Input files: VEP_consequence_graphs.csv VEP consequences per gene per variant

transcript_id,ensembl_gene_id,HGNC,VEP_DEF_LABELS,VAR
ENST00000300305,ENSG00000159216,RUNX1,INTRON,chr21_36280376_A_G
ENST00000416754,ENSG00000159216,RUNX1,INTRON,chr21_36280376_A_G

################### Input parameters: excluded_phenotypes=$(echo "wbc,eo_p,mono_p,neut_p,lymph_p,baso_p")

################### Code Main points:

ALL_dB_subset_restricted<-unique(ALL_dB_subset[-which(ALL_dB_subset$phenotype%in%excluded_phenotypes),]) # Exclude phenotypes of white blood cells that are correlated

VEP_CSQ_subset_with_genes<-unique(VEP_CSQ_subset[which(VEP_CSQ_subset$ensembl_gene_id != "DUMMY"),]) # Exclude VEP consequences that are not linked to genes

VEP_and_PCHiC<-unique(rbind(VEP_CSQ_subset_with_genes,PCHiC_subset)) # Get a data frame with the genes linked to a variant by VEP consequences and PCHiC interactions


master_path<-"/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/COGS/rCOGS_out/COGS_scores/"
file_list <- list.files(path=master_path, include.dirs = FALSE)
file_list_sel<-file_list[grep("_interest_cell_types_",file_list)]
phenotypes<-gsub("_COGS_interest_cell_types_with_gene_names","", file_list_sel)
df<-as.data.frame(cbind(phenotypes,file_list_sel), stringsAsFactors=F)
df_restricted<-unique(df[-which(df$phenotype%in%excluded_phenotypes),]) # Obtain a data frame with the phenotypes and the filename for the COGS files calculated for the interest cell types (Ery,MK,Mono,Neut,tCD4,tCD8 and tB)

phenotypes_array<-unique(df_restricted$phenotype)
for(i in 1:length(phenotypes_array)) # Loop per phenotype

df_restricted_sel<-df_restricted[which(df_restricted$phenotype == phenotypes_array_sel),]
file_selected_to_open<-unique(df_restricted_sel$file)
COGS<-as.data.frame(fread(file=file_selected_to_open,sep="\t") , stringsAsFactors=F) # Open the COGs interest cell types (Ery,MK,Mono,Neut,tCD4,tCD8 and tB) file for that phenotype

ALL_dB_subset_restricted_sel<-ALL_dB_subset_restricted[which(ALL_dB_subset_restricted$phenotype == phenotypes_array_sel),] # Select all variants associated to that phenotype

VEP_and_PCHiC_sel<-VEP_and_PCHiC[which(VEP_and_PCHiC$VAR%in%ALL_dB_subset_restricted_sel$VAR),] # Select all the genes connected to those variants by VEP consequences or PCHiC

COGS_subset<-merge(VEP_and_PCHiC_sel,
                           COGS_subset,
                           by="ensembl_gene_id") # Merge COGs values with the genes; we only keep the COGs values of the genes linked to variants by VEP and/or PCHiC

COGS_subset$phenotype<-phenotypes_array_sel # Add the phenotype field
Gather_NO_NA<-Gather[!is.na(Gather$cogs),] # Exclude NAs

################### Output files: COGS_GLOBAL.tsv

ensembl_gene_id VAR     HGNC    cogs    phenotype
ENSG00000001167 chr6_15254908_G_A       NFYA    6.26313588159011e-06    ret_p
ENSG00000002549 chr4_18018944_C_T       LAP3    0.377508796649921       ret_p
ENSG00000002549 chr4_18032493_G_A       LAP3    0.377508796649921       ret_p
ENSG00000002549 chr4_18019572_C_T       LAP3    0.377508796649921       ret_p

# 335_COGS.R

################### Input files: COGS_GLOBAL.tsv

ensembl_gene_id VAR     HGNC    cogs    phenotype
ENSG00000001167 chr6_15254908_G_A       NFYA    6.26313588159011e-06    ret_p
ENSG00000002549 chr4_18018944_C_T       LAP3    0.377508796649921       ret_p
ENSG00000002549 chr4_18032493_G_A       LAP3    0.377508796649921       ret_p
ENSG00000002549 chr4_18019572_C_T       LAP3    0.377508796649921       ret_p

################### Code main points Second Function:


COGS$mean_cogs<-mean(COGS$cogs, na.rm =T)
COGS$sd_cogs<-sd(COGS$cogs, na.rm =T)
COGS$cogs_Z_score<-(COGS$cogs-COGS$mean_cogs)/COGS$sd_cogs # Z-score normalisation


################### Output files: Prepared_file_COGS.rds
VAR     ensembl_gene_id HGNC    phenotype       value   value_Z_score   variable
chr18_60920854_C_T      ENSG00000171791 BCL2    rbc     0.997050456310743       1.97861233188362        COGS
chr3_184091102_T_G      ENSG00000090534 THPO    plt     0       -0.842945974185377      COGS
chr1_202129205_G_A      ENSG00000143851 PTPN7   plt     6.995914381136e-06      -0.842926176410645      COGS
chr3_184091102_T_G      ENSG00000145191 EIF2B5  plt     0       -0.842945974185377      COGS

# 350_binder_of_scores_pLOEUF.R

################### Input files: gnomAD LOF_oe data gnomad.v2.1.1.lof_metrics.by_gene.txt

gene    transcript      obs_mis exp_mis oe_mis  mu_mis  possible_mis    obs_mis_pphen   exp_mis_pphen   oe_mis_pphen    possible_mis_pphen      obs_syn exp_syn oe_syn  mu_syn  possible_syn    obs_lof mu_lof  possible_lof    exp_lof pLI     pNull   pRec    oe_lof  oe_syn_lower    oe_syn_upper    oe_mis_lower    oe_mis_upper oe_lof_lower    oe_lof_upper    constraint_flag syn_z   mis_z   lof_z   oe_lof_upper_rank       oe_lof_upper_bin        oe_lof_upper_bin_6      n_sites classic_caf     max_af  no_lofs obs_het_lof     obs_hom_lof     defined p       exp_hom_lof     classic_caf_afr classic_caf_amr classic_caf_asj classic_caf_eas      classic_caf_fin classic_caf_nfe classic_caf_oth classic_caf_sas p_afr   p_amr   p_asj   p_eas   p_fin   p_nfe   p_oth   p_sas   transcript_type gene_id transcript_level        cds_length      num_coding_exons        gene_type       gene_length     exac_pLI        exac_obs_lof    exac_exp_lof    exac_oe_lof brain_expression chromosome      start_position  end_position
MED13   ENST00000397786 871     1.1178e+03      7.7921e-01      5.5598e-05      14195   314     5.2975e+02      5.9273e-01      6708    422     3.8753e+02      1.0890e+00      1.9097e-05      4248    0       4.9203e-06      1257    9.8429e+01      1.0000e+00      8.9436e-40      1.8383e-16      0.0000e+00      1.0050e+00   1.1800e+00      7.3600e-01      8.2400e-01      0.0000e+00      3.0000e-02              -1.3765e+00     2.6232e+00      9.1935e+00      0       0       0       2       1.2058e-05      8.0492e-06      124782  3       0       124785  1.2021e-05      1.8031e-05      0.0000e+00      0.0000e+00      0.0000e+00  0.0000e+00       9.2812e-05      8.8571e-06      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      9.2760e-05      8.8276e-06      0.0000e+00      0.0000e+00      protein_coding  ENSG00000108510 2       6522    30      protein_coding  122678  1.0000e+00      0       6.4393e+01   0.0000e+00      NA      17      60019966        60142643

################### Input files: PCHIC_ChicagoScore_graphs.csv
VAR,HGNC,ensembl_gene_id,value,CellType_DEF
chr21_36789420_C_G,AP000330.8,ENSG00000234380,5.509692830975,aCD4
chr21_36789420_C_G,AP000330.8,ENSG00000234380,5.59914383391615,tCD4

################### Input files: VEP_consequence_graphs.csv VEP consequences per gene per variant

transcript_id,ensembl_gene_id,HGNC,VEP_DEF_LABELS,VAR
ENST00000300305,ENSG00000159216,RUNX1,INTRON,chr21_36280376_A_G
ENST00000416754,ENSG00000159216,RUNX1,INTRON,chr21_36280376_A_G

################### Code Main points:

GENE_PLOEUF_NO_NA<-GENE_PLOEUF[!is.na(GENE_PLOEUF$oe_lof),] # Exclude NAs in oe_lof

VEP_CSQ_subset_with_genes<-unique(VEP_CSQ_subset[which(VEP_CSQ_subset$ensembl_gene_id != "DUMMY"),]) # Exclude VEP consequences that are not linked to genes

VEP_and_PCHiC<-unique(rbind(VEP_CSQ_subset_with_genes,PCHiC_subset)) # Get a data frame with the genes linked to a variant by VEP consequences and PCHiC interactions

GENE_PLOEUF_NO_NA<-merge(GENE_PLOEUF_NO_NA,
                     VEP_and_PCHiC,
                     by=c("ensembl_gene_id","HGNC")) # Keep the oe_lof values for the genes that can be linked to the variants by VEP or PCHiC

GENE_PLOEUF_NO_NA_subset.m<-melt(GENE_PLOEUF_NO_NA_subset, id.vars=c("VAR","HGNC","ensembl_gene_id")) # melt by "VAR","HGNC","ensembl_gene_id"

mean_LOEUF<-mean(GENE_PLOEUF_NO_NA_subset.m$value)
sd_LOEUF<-sd(GENE_PLOEUF_NO_NA_subset.m$value)
GENE_PLOEUF_NO_NA_subset.m$value_Z_score<-(GENE_PLOEUF_NO_NA_subset.m$value-mean_LOEUF)/sd_LOEUF # Z-score normalisation of the oe_lof values

################### Output files: Prepared_file_GENE_PLOEUF.rds

VAR     HGNC    ensembl_gene_id variable        value   value_Z_score
chr3_128322617_G_A      SEC61A1 ENSG00000058262 oe_lof  0.11405 -0.889375188828683
chr3_184091102_T_G      THPO    ENSG00000090534 oe_lof  0.27137 -0.487842804843122

# 332_binder_of_Prepared_files_and_add_annotation_layer.R

################### Input files from the pipeline:

ATAC_data=$(echo "Prepared_file_ATAC.rds")
multi_ATAC_data=$(echo "Prepared_file_multi_lineage_ATAC.rds")
COGS=$(echo "Prepared_file_COGS.rds")
oe_LOF=$(echo "Prepared_file_GENE_PLOEUF.rds")
GENE_EXP=$(echo "Prepared_file_GENE_EXP.rds")
PCHiC=$(echo "Prepared_file_PCHiC.rds")
chromstates=$(echo "Prepared_file_chromstates.rds")
CADD=$(echo "Prepared_file_CADD.rds")
Constraint_Z=$(echo "Prepared_file_Constraint_Z.rds")
NCBoost=$(echo "Prepared_file_NCBoost.rds")
SpliceAI=$(echo "Prepared_file_SpliceAI.rds")
GWAS_GLOBAL_per_traits=$(echo "GWAS_GLOBAL_per_traits.tsv")
MAF_GLOBAL=$(echo "MAF_GLOBAL.tsv")

################### Input files ER_Labelling_Initial_Selection.rds Variant prioritisation bins

'data.frame':   178890 obs. of  7 variables:
 $ VAR                : chr  "chr10_101222300_G_A" "chr10_101237412_T_C" "chr10_101244819_G_A" "chr10_101245230_T_C" ...
 $ rs                 : chr  "rs117746067" "rs117156560" "rs2494662" "rs75541431" ...
 $ maf_origin         : num  0.104 0.104 0.272 0.107 0.275 ...
 $ Fig0_Annot_Category: Ord.factor w/ 3 levels "Common_variant"<..: 1 1 1 1 1 1 1 1 1 1 ...
 $ Fig1_Annot_Category: Ord.factor w/ 5 levels "Common_variant"<..: 1 1 1 1 1 1 1 1 1 1 ...
 $ Fig2_Annot_Category: Ord.factor w/ 2 levels "NON_ACTIVE"<"ACTIVE": NA NA NA NA NA NA NA NA NA NA ...
 $ Fig3_Annot_Category: Ord.factor w/ 3 levels "EDITION_DROPOUT"<..: NA NA NA NA NA NA NA NA NA NA ...

################### Input files Supp_Table_4_CURATED_Plus_phenotypes.rds Table with the MPRA, curation and Lineage restricted vs multilineage labels

> str(df)
'data.frame':   100 obs. of  22 variables:
 $ VAR                        : chr  "chr10_9032555_C_G" "chr1_101744646_C_G" "chr1_158613314_G_A" "chr1_158638260_G_A" ...
 $ rs                         : chr  "rs72782676" "rs77046277" "rs200830867" "rs140446749" ...
 $ Mechanistic_Class          : Ord.factor w/ 5 levels "Transcriptional_Regulation"<..: 5 4 4 4 1 1 5 1 4 4 ...
 $ Manual_curation            : Ord.factor w/ 7 levels "R_in_candidate"<..: 7 3 3 4 1 1 7 1 3 3 ...
 $ Candidate_effector         : chr  "ABSENT" "S1PR1" "SPTA1" "SPTA1" ...
 $ OpenTargets_QTL            : chr  "NOT FOUND" "NOT FOUND" "NOT FOUND" "NOT FOUND" ...
 $ Observations               : chr  "" "" "" "rs138055271__0.74__SPTA1(UPSTREAM,NCT,DOWNSTREAM,MISS)|rs148912436__1__SPTA1(MISS)" ...
 $ Whole_blood_DE_HGNC_string : chr  NA NA NA NA ...
 $ Whole_blood_DTU_HGNC_string: chr  NA NA NA NA ...
 $ Monocyte_DE_HGNC_string    : chr  NA NA NA NA ...
 $ Tcell_DE_HGNC_string       : chr  NA NA NA NA ...
 $ Neutrophil_DE_HGNC_string  : chr  NA NA NA NA ...
 $ Proxy_rsid_string          : chr  NA "rs77372977__0.86__PPIAP7(DOWNSTREAM)" "rs369186418__0.69__NA(INTERGENIC)|rs555280966__0.67__NA(REGULATORY);EI24P2(DOWNSTREAM)|rs183408858__0.85__NA(IN"| __truncated__ "rs140554024__0.61__OR10R2(DOWNSTREAM);RP11-144L1.4(NCT)|rs138055271__0.74__SPTA1(UPSTREAM,NCT,DOWNSTREAM,MISS)|"| __truncated__ ...
 $ MPRA_CLASS                 : Ord.factor w/ 4 levels "NO_enhancer_activity"<..: 2 1 2 3 3 3 2 4 2 3 ...
 $ genIE_CLASS                : Ord.factor w/ 3 levels "genIE_INACTIVE"<..: 3 3 3 3 3 3 3 2 1 2 ...
 $ Replication_OT_QTL         : chr  "No_RNA_Seq_HET_carriers" NA NA NA ...
 $ VAR_38                     : chr  "chr10_8990592_C_G" "chr1_101279090_C_G" "chr1_158643524_G_A" "chr1_158668470_G_A" ...
 $ Multi_Lineage              : Ord.factor w/ 2 levels "Lineage_restricted"<..: 1 1 1 1 1 2 1 1 1 1 ...
 $ Lineage_string             : chr  "granulocyte_monocyte_lineage" "lymphocyte_lineage" "erythroid_lineage" "erythroid_lineage" ...
 $ phenotype_DEF_string       : chr  "eo;eo_p" "mono_p;neut_p;lymph;lymph_p;wbc" "mchc;ret;ret_p;hlr;hlr_p;mrv;mscv" "mcv;mch;mchc;rdw_cv;ret;ret_p;hlr;hlr_p;mrv;mscv" ...
 $ maf_origin                 : num  0.00698 0.00151 0.00909 0.01156 0.0099 ...
 $ VEP_DEF_LABELS_wCSQ        : Ord.factor w/ 15 levels "LOF"<"MISS"<"SYN"<..: 7 7 6 6 6 6 7 6 6 6

################### Input parameters: RMV_common=$(echo "chr4_1008212_C_T") Removes the variant chr4_1008212_C_T that went above 3% MAF
################### Input parameters: RMV_labels=$(echo 'NOT_SCREENED_MPRA,NOT_SCREENED_genIE,No_RNA_Seq_HET_carriers') removes these labels in the appropriate figures to simplify the representation


################### Code Main points:

ind.int<-c(which(colnames(NCBoost) == "VAR"),which(colnames(NCBoost) == "value"),which(colnames(NCBoost) == "value_Z_score"),which(colnames(NCBoost) == "variable"))
NCBoost_subset<-NCBoost[,ind.int]
ind.int<-c(which(colnames(SpliceAI) == "VAR"),which(colnames(SpliceAI) == "value"),which(colnames(SpliceAI) == "value_Z_score"),which(colnames(SpliceAI) == "variable"))
SpliceAI_subset<-SpliceAI[,ind.int] 
ind.int<-c(which(colnames(oe_LOF) == "VAR"),which(colnames(oe_LOF) == "value"),which(colnames(oe_LOF) == "value_Z_score"),which(colnames(oe_LOF) == "variable"))
oe_LOF_subset<-oe_LOF[,ind.int] 
ind.int<-c(which(colnames(GENE_EXP) == "VAR"),which(colnames(GENE_EXP) == "value"),which(colnames(GENE_EXP) == "value_Z_score"),which(colnames(GENE_EXP) == "variable"))
GENE_EXP_subset<-GENE_EXP[,ind.int] #For NCBoost, SpliceAI, oe_lof and GENE_EXP do not select the associated gene

ind.int<-c(which(colnames(COGS) == "VAR"),which(colnames(COGS) == "value"),which(colnames(COGS) == "value_Z_score"),which(colnames(COGS) == "variable"))
COGS_subset<-COGS[,ind.int] #For COGS do not select the associated gene and phenotype

ind.int<-c(which(colnames(GWAS_GLOBAL_per_traits) == "VAR"),which(colnames(GWAS_GLOBAL_per_traits) == "value"),which(colnames(GWAS_GLOBAL_per_traits) == "value_Z_score"),which(colnames(GWAS_GLOBAL_per_traits) == "variable"))
GWAS_GLOBAL_per_traits_subset_1<-GWAS_GLOBAL_per_traits[,ind.int] # For GWAS do not select the associated phenotype

Master_file<-rbind(df_1,df_2,df_3)
Master_file$variable<-factor(Master_file$variable,
                               levels=c("MAF","PP","Absolute_effect_size","CADD_raw",'Constraint_Z',"NCBoost","SpliceAI_DG","SpliceAI_DL","SpliceAI_AG","SpliceAI_AL",
                                        "Rank_ATAC_erythroid_lineage","Rank_ATAC_mega_lineage","Rank_ATAC_gran_mono_lineage","Rank_ATAC_lymph_lineage","multi_lineage_ATAC","Rank_PCHiC","Rank_chromstates",
                                        "COGS","oe_lof","Rank_GENE_EXP"),
                               ordered=T) # bind all the features and order them

Table_S4<-Table_S4[-which(Table_S4$VAR == RMV_common),]
Table_S4_subset$MPRA_CLASS_2[which(Table_S4_subset$MPRA_CLASS_2 == RMV_labels[1])]<-NA
Table_S4_subset$genIE_CLASS[which(Table_S4_subset$genIE_CLASS == RMV_labels[2])]<-NA
Table_S4_subset$Mechanistic_Class[which(Table_S4_subset$Mechanistic_Class == RMV_labels[3])]<-NA
Table_S4_subset$Manual_curation[which(Table_S4_subset$Manual_curation == RMV_labels[3])]<-NA
Table_S4_subset$M_and_M<-interaction(Table_S4_subset$Mechanistic_Class,Table_S4_subset$Manual_curation, sep="|",lex.order = T)
Table_S4_subset$MPRA_CLASS_2<-revalue(Table_S4_subset$MPRA_CLASS,
                                        c("NO_enhancer_activity" = "no_hit",
                                          "AT_LEAST_1_TILE_with_enhancer_activity" = "no_hit",
                                          "AT_LEAST_1_TILE_with_E_Plus_ASE_activity" = "hit")) # Process the Table_S4

Table_of_labels<-merge(VAR_Prioritization_dB_subset,
                         Table_S4_subset,
                         by="VAR",
                         all=T) # Merge all the variant prioritisation tables of labels

################### Output files: Master_file_scores.rds

'data.frame':   3632292 obs. of  4 variables:
 $ VAR          : chr  "chr10_101222300_G_A" "chr10_101237412_T_C" "chr10_101244819_G_A" "chr10_101245230_T_C" ...
 $ value        : num  0.104 0.104 0.272 0.107 0.275 ...
 $ value_Z_score: num  -1.20876 -1.20612 0.00997 -1.18414 0.03356 ...
 $ variable     : Ord.factor w/ 20 levels "MAF"<"PP"<"Absolute_effect_size"<..: 1 1 1 1 1 1 1 1 1 1 ...

################### Output files: Table_of_labels.rds

'data.frame':   178890 obs. of  9 variables:
 $ VAR                : chr  "chr10_101222300_G_A" "chr10_101237412_T_C" "chr10_101244819_G_A" "chr10_101245230_T_C" ...
 $ Fig1_Annot_Category: Ord.factor w/ 5 levels "Common_variant"<..: 1 1 1 1 1 1 1 1 1 1 ...
 $ Mechanistic_Class  : Ord.factor w/ 4 levels "Transcriptional_Regulation"<..: NA NA NA NA NA NA NA NA NA NA ...
 $ Manual_curation    : Ord.factor w/ 6 levels "R_in_candidate"<..: NA NA NA NA NA NA NA NA NA NA ...
 $ MPRA_CLASS         : Ord.factor w/ 4 levels "NO_enhancer_activity"<..: NA NA NA NA NA NA NA NA NA NA ...
 $ genIE_CLASS        : Ord.factor w/ 2 levels "genIE_INACTIVE"<..: NA NA NA NA NA NA NA NA NA NA ...
 $ Multi_Lineage      : Ord.factor w/ 2 levels "Lineage_restricted"<..: NA NA NA NA NA NA NA NA NA NA ...
 $ MPRA_CLASS_2       : Ord.factor w/ 2 levels "no_hit"<"hit": NA NA NA NA NA NA NA NA NA NA ...
 $ M_and_M            : Factor w/ 14 levels "Transcriptional_Regulation|R_in_candidate",..: NA NA NA NA NA NA NA NA NA NA ...

# 333_stats_on_Z_score.R

################### Input files: Master_file_scores.rds

'data.frame':   3632292 obs. of  4 variables:
 $ VAR          : chr  "chr10_101222300_G_A" "chr10_101237412_T_C" "chr10_101244819_G_A" "chr10_101245230_T_C" ...
 $ value        : num  0.104 0.104 0.272 0.107 0.275 ...
 $ value_Z_score: num  -1.20876 -1.20612 0.00997 -1.18414 0.03356 ...
 $ variable     : Ord.factor w/ 20 levels "MAF"<"PP"<"Absolute_effect_size"<..: 1 1 1 1 1 1 1 1 1 1 ...

################### Input files: Table_of_labels.rds

'data.frame':   178890 obs. of  9 variables:
 $ VAR                : chr  "chr10_101222300_G_A" "chr10_101237412_T_C" "chr10_101244819_G_A" "chr10_101245230_T_C" ...
 $ Fig1_Annot_Category: Ord.factor w/ 5 levels "Common_variant"<..: 1 1 1 1 1 1 1 1 1 1 ...
 $ Mechanistic_Class  : Ord.factor w/ 4 levels "Transcriptional_Regulation"<..: NA NA NA NA NA NA NA NA NA NA ...
 $ Manual_curation    : Ord.factor w/ 6 levels "R_in_candidate"<..: NA NA NA NA NA NA NA NA NA NA ...
 $ MPRA_CLASS         : Ord.factor w/ 4 levels "NO_enhancer_activity"<..: NA NA NA NA NA NA NA NA NA NA ...
 $ genIE_CLASS        : Ord.factor w/ 2 levels "genIE_INACTIVE"<..: NA NA NA NA NA NA NA NA NA NA ...
 $ Multi_Lineage      : Ord.factor w/ 2 levels "Lineage_restricted"<..: NA NA NA NA NA NA NA NA NA NA ...
 $ MPRA_CLASS_2       : Ord.factor w/ 2 levels "no_hit"<"hit": NA NA NA NA NA NA NA NA NA NA ...
 $ M_and_M            : Factor w/ 14 levels "Transcriptional_Regulation|R_in_candidate",..: NA NA NA NA NA NA NA NA NA NA ...

################### Code Main points:

category_vector<-c("Fig1_Annot_Category","MPRA_CLASS_2","M_and_M","Multi_Lineage")
for(i in 1:length(category_vector)) # LOOP 1 per category
for(k in 1:length(array_variable)) # LOOP 2 per feature

Summary_table<-as.data.frame(Master_file_sel.dt[,.(n=.N,
                                                           Min=round(as.numeric(summary(value_Z_score)[1]),3),
                                                           Q1=round(as.numeric(summary(value_Z_score)[2]),3),
                                                           M=round(as.numeric(summary(value_Z_score)[3]),3),
                                                           Q3=round(as.numeric(summary(value_Z_score)[5]),3),
                                                           Max=round(as.numeric(summary(value_Z_score)[6]),3)),, by=key(Master_file_sel.dt)], stringsAsFactors=F)
colnames(Summary_table)[which(colnames(Summary_table) == colnames(Master_file_sel)[indx.category_Master_sel])]<-"category2"
Summary_table$variable<-array_variable_sel # calculate key parameters of the distribution of the feature

PW_category<-pairwise.wilcox.test(Master_file_sel$value, Master_file_sel[,indx.category_Master_sel],
                                            p.adjust.method = "BH") # Pairwise Wilcoxon test for all the comparisons

################### Output files: STATS_value_Z_score.rds

'data.frame':   137 obs. of  19 variables:
 $ category2     : Factor w/ 23 levels "Transcriptional_Regulation|R_in_candidate",..: 17 18 19 19 17 18 19 19 17 18 ...
 $ category1     : Factor w/ 23 levels "Transcriptional_Regulation|R_in_candidate",..: 19 19 15 16 19 19 15 16 19 19 ...
 $ variable      : chr  "MAF" "MAF" "MAF" "MAF" ...
 $ pval          : num  2.31e-28 6.43e-12 1.19e-81 4.25e-10 7.31e-223 ...
 $ MINUS_logpval : num  27.64 11.19 80.92 9.37 222.14 ...
 $ n_category2   : int  5052 73 123 123 7426 153 394 394 7410 153 ...
 $ n_category1   : int  123 123 173077 565 394 394 357974 1281 394 394 ...
 $ M_category2   : num  -1.836 -1.848 -1.916 -1.916 -0.255 ...
 $ M_category1   : num  -1.916 -1.916 0.108 -1.869 7.112 ...
 $ Min_category2 : num  -1.96 -1.95 -1.96 -1.96 -0.37 ...
 $ Min_category1 : num  -1.964 -1.964 -1.964 -1.963 -0.336 ...
 $ Q1_category2  : num  -1.884 -1.892 -1.947 -1.947 -0.309 ...
 $ Q1_category1  : num  -1.95 -1.95 -0.71 -1.93 6.42 ...
 $ Q3_category2  : num  -1.78 -1.81 -1.88 -1.88 0.05 ...
 $ Q3_category1  : num  -1.882 -1.882 0.879 -1.811 7.725 ...
 $ Max_category2 : num  -1.75 -1.75 -1.75 -1.75 7.79 ...
 $ Max_category1 : num  -1.75 -1.75 1.67 -1.75 8.89 ...
 $ Annotation    : chr  "Fig1_Annot_Category" "Fig1_Annot_Category" "Fig1_Annot_Category" "Fig1_Annot_Category" ...
 $ Related_figure: Ord.factor w/ 4 levels "Figure_1"<"Figure_S2"<..: 1 1 1 1 1 1 1 1 1 1 ...

# 348_Z_score_violin_plots.R

################### Input files from pipeline

Master_file=$(echo "Master_file_scores.rds")
Table_of_labels=$(echo "Table_of_labels.rds")
Stats_file=$(echo "STATS_value_Z_score.rds")
Categories_colors=$(echo "ER_Labelling_Categories_colors.rds")

################### Input parameters, which variants to represent in two separate graphs one for the variant based features and one for the gene based features

subset_variables_1=$(echo "PP,Absolute_effect_size,CADD_raw,Constraint_Z,NCBoost,SpliceAI_DG,SpliceAI_DL,SpliceAI_AG,SpliceAI_AL,Rank_ATAC_erythroid_lineage,Rank_ATAC_mega_lineage,Rank_ATAC_gran_mono_lineage,Rank_ATAC_lymph_lineage,multi_lineage_ATAC,Rank_PCHiC,Rank_chromstates")
subset_variables_2=$(echo "COGS,oe_lof,Rank_GENE_EXP")

################### Code Main points:

################### Dotplot for the minus logpval values and violin plots to show the distributions