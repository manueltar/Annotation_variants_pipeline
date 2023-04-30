# Annotation_variants_pipeline


#  Master script command line:

$ bash 313_new_calculator_of_scores.sh <Path_to_output_files> <Memory allocation (Mb)> <Processors> <queue>

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

$ head GWAS_GLOBAL_per_traits.tsv
VAR     phenotype       variable        value   value_Z_score
chr2_181894386_CT_C     baso    Absolute_effect_size    0.0147284       -0.409226791052507

$ head MAF_GLOBAL.tsv
VAR     value   value_Z_score   variable
chr10_101222300_G_A     0.104047        -1.20875664040154       MAF

# Subscript: 340_binder_CADD.R 



