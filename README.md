# Annotation_variants_pipeline


#  Master script command line:

$ bash 313_new_calculator_of_scores.sh Path_to_output_files Memory Processors queue

e.g. $ bash 313_new_calculator_of_scores.sh Annotation_results/ 4000 1 normal


# Subscript: 346_binder_GWAS_and_Z_score.R 

################### Input files: ALL_dB.tsv (annotation of GWAS blood traits)

$ head /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv

VAR     variant_id      rs      chr     pos37   ref     alt     maf_origin      block_no        phenotype       ncondind        credset_size    finemap_beta    finemap_se      finemap_z       finemap_prob    finemap_log10bf info    is_cond_ind     chromstates
chr1_1689164_C_T        4154    rs2235536       chr1    1689164 C       T       0.425205        1       ret     1       107     0.0131952       0.00215928      6.11093 0.00330496      1.44972 0.996973        0       E2C1:E2C3:E2C5:E2C6:E2C7:E2C8:E2C9:E2C10:E2C11:E2C12:E9C14:E2C15:E2C16:E2C17:E9C18:E2C19:E2C20:E2C21:E2C22:E2C23:E2C24:E2C25:E2C26:E2C27:E2C28:E2C29:E2C30:E2C31:E2C32:E2C33:E2C34:E2C35:E2C36:E2C37:E2C38:E2C39

################### Pseudocode:

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

################### Pseudocode:

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


################### Pseudocode:

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

################### Pseudocode:

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

################### Pseudocode:

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

################### Pseudocode:

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
################### Input files: haemvar_manuel_200903_NCBOOST.tsv (annotation of NCBoost results)

$ head /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NCBoost/haemvar_manuel_200903_NCBOOST.tsv
chr     pos     ref     alt     annovar_annotation      closest_gene_name       gene_type       pLI     familyMemberCount       ncRVIS  ncGERP  RVIS_percentile slr_dnds        GDI     gene_age        UTR3    UTR5    downstream      intergenic      intronic        upstream        GC      CpG     priPhCons       mamPhCons    verPhCons       priPhyloP       mamPhyloP       verPhyloP       GerpN   GerpS   bStatistic      TajimasD_YRI_pvalue     TajimasD_CEU_pvalue     TajimasD_CHB_pvalue     FuLisD_YRI_pvalue       FuLisD_CEU_pvalue       FuLisD_CHB_pvalue       FuLisF_YRI_pvalue       FuLisF_CEU_pvalue       FuLisF_CHB_pvalue   meanDaf1000G     meanHet1000G    meanMAF1000G    meanMAFGnomAD   meanMAF_AFRGnomAD       meanMAF_AMRGnomAD       meanMAF_ASJGnomAD       meanMAF_EASGnomAD       meanMAF_FINGnomAD       meanMAF_NFEGnomAD       meanMAF_OTHGnomAD       CDTS    NCBoost
1       1689164 C       T       intronic        NADK    protein_coding  0.003878536     0       -0.3783625      -0.7206092      45.35858        0.06153 176.4511        0       0       0       0       0       1       0       0.6     0.04    0.028   0.008   0.026   0.457   1.003   1.317   1.86    1.86    662     1.023368     0.1377089       0.08570405      1.024593        0.5245587       0.8809586       1.149016        0.3383352       0.4969653       0.1478  0.01082257      0.00159 0.000709        0.00159 0.000294        0.000432        0.000645        0.000452        0.000302        0.000347        -0.126957       0.06201531
1

################### Pseudocode:

  NCBoost_result$VAR<-paste(paste('chr',NCBoost_result$chr, sep=''),NCBoost_result$pos,NCBoost_result$ref,NCBoost_result$alt,sep="_")
 Merge_table<-merge(ALL_dB_subset,
                     NCBoost_result,
                     by=c("VAR")) # Merge NCBoost result with GWAS results and keep the overlap

Merge_table_subset_NO_NA<-Merge_table_subset[!is.na(Merge_table_subset$NCBoost),] # exclude NA's


################### Output files:

$ head NCBoost_GLOBAL.tsv
VAR     closest_gene_name       NCBoost
chr10_101222300_G_A     GOT1    0.015575

