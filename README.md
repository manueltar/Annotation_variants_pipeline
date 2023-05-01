# Annotation_variants_pipeline


#  Master script command line:

$ bash 313_new_calculator_of_scores.sh Path_to_output_files Memory Processors queue

e.g. $ bash 313_new_calculator_of_scores.sh Annotation_results/ 4000 1 normal


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

























