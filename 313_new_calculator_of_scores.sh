#!/bin/bash>

MASTER_ROUTE=$1
output_dir=$MASTER_ROUTE
mem=$2
pc=$3
queue=$4



output="/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/313_CREATION.sh"

touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output


#### Rscript

Rscript=/software/R-4.1.0/bin/Rscript




echo "###################################################### MASTER_LABELLING ########### MASTER_LABELLING ########### MASTER_LABELLING #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_ER_Labelling=/nfs/users/nfs_m/mt19/Scripts/R/235_ER_Labelling_v9.R



MPRA_Tier_1=$(echo "E_Plus_ASE_ACTIVE_AT_LEAST_1")
genIE_Tier_1=$(echo "Tier_del")

batch=$(echo "MPRA_""$MPRA_Tier_1""_genIE_""$genIE_Tier_1")
finemap_prob_Threshold=$(echo '0.1')

echo "BATCH:$batch"
echo "finemap_prob_Threshold:$finemap_prob_Threshold"

maf_Threshold=$(echo '0.03')
Info_Score_Threshold=$(echo '0.8')
Prob_Threshold=$(echo '0.9')


ALL_dB_Hannes=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/mt19/ALL_db.tsv')
ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
path_Patrick=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/')




type=$(echo "ER_Labelling")
Initial_Selection=$(echo "$output_dir""ER_Labelling_Initial_Selection.rds")

outfile_ER_Labelling=$(echo "$output_dir""outfile_""$type""_""$batch""_""$finemap_prob_Threshold"".out")
touch $outfile_ER_Labelling
echo -n "" > $outfile_ER_Labelling
name_ER_Labelling=$(echo "$type""_""$batch""_""$finemap_prob_Threshold""_job")


MPRA_Tier_0_name=$MPRA_Tier_1
MPRA_Tier_0_name=$(echo $MPRA_Tier_0_name|sed -e 's/_AT_LEAST_.*$/_0/g')

echo "$MPRA_Tier_0_name"

MPRA_Tier_0=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/GLOBAL_ANALYSIS_Feb_2022_FC01/""$MPRA_Tier_0_name"".rds")
MPRA_Tier_1=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/GLOBAL_ANALYSIS_Feb_2022_FC01/""$MPRA_Tier_1"".rds")
genIE_input=$(echo '/nfs/users/nfs_m/mt19/analysis_genIE/input.csv')
genIE_RESULTS=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/GLOBAL_ANALYSIS/genIE_Tiers_DEF.rds")


paleo_file=$(echo "/nfs/users/nfs_m/mt19/RareVar_2019/varbase190212_full.csv")

step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"



echo "bsub -G team151 -o $outfile_ER_Labelling -M $step_mem  -J $name_ER_Labelling -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_ER_Labelling \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--ALL_dB_Hannes $ALL_dB_Hannes \\" >> $output
echo "--paleo_file $paleo_file \\" >> $output
echo "--path_Patrick $path_Patrick \\" >> $output
echo "--MPRA_Tier_0 $MPRA_Tier_0 \\" >> $output
echo "--MPRA_Tier_1 $MPRA_Tier_1 \\" >> $output
echo "--genIE_RESULTS $genIE_RESULTS \\" >> $output
echo "--genIE_input $genIE_input \\" >> $output
echo "--maf_Threshold $maf_Threshold \\" >> $output
echo "--Info_Score_Threshold $Info_Score_Threshold \\" >> $output
echo "--Prob_Threshold $Prob_Threshold \\" >> $output
echo "--Initial_Selection $Initial_Selection \\" >> $output
echo "--finemap_prob_Threshold $finemap_prob_Threshold \\" >> $output
echo "--type $type --out $output_dir\"" >> $output




echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "####################################################### MSC GRAPH #############################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output



Rscript_VEP_CSQ=/nfs/users/nfs_m/mt19/Scripts/R/235_ER_VEP_CSQ_v5.R

VEP_CSQ=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/VEP_consequence_graphs.csv')
ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
CSQ_colors=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/df_CSQ_colors.rds")


Initial_Selection_WITH_CSQ_labels=$(echo "$output_dir""ER_Labelling_Initial_Selection_with_CSQ_labels.rds")
Categories_colors=$(echo "$output_dir""ER_Labelling_Categories_colors.rds")

type=$(echo "VEP_CSQ")

outfile_ER_VEP_CSQ=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_ER_VEP_CSQ
echo -n "" > $outfile_ER_VEP_CSQ
name_ER_VEP_CSQ=$(echo "$type""_job")


echo "bsub -G team151 -o $outfile_ER_VEP_CSQ -M $mem -w\"done($name_ER_Labelling)\" -J $name_ER_VEP_CSQ -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_VEP_CSQ \\" >> $output
echo "--Initial_Selection $Initial_Selection_WITH_CSQ_labels \\" >> $output
echo "--Categories_colors $Categories_colors \\" >> $output
echo "--CSQ_colors $CSQ_colors \\" >> $output
echo "--VEP_CSQ $VEP_CSQ \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--type $type --out $output_dir\"" >> $output






echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "####################################################### VARIANT BASED FEATURES #############################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output


echo "###################################################### binder_of_scores_GWAS  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_binder_of_scores_GWAS=/nfs/users/nfs_m/mt19/Scripts/R/346_binder_GWAS_and_Z_score.R


ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
type=$(echo "binder_of_scores_GWAS")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")

outfile_binder_of_scores_GWAS=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_binder_of_scores_GWAS
echo -n "" > $outfile_binder_of_scores_GWAS
name_binder_of_scores_GWAS=$(echo "$type""_job")

 
step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"

 

echo "bsub -G team151 -o $outfile_binder_of_scores_GWAS -M $step_mem  -J $name_binder_of_scores_GWAS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_binder_of_scores_GWAS \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output


echo "###################################################### binder_of_scores_CADD  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_binder_of_scores_CADD=/nfs/users/nfs_m/mt19/Scripts/R/340_binder_CADD.R


ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
type=$(echo "binder_of_scores_CADD")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
CADD_result=$(echo "/lustre/scratch123/hgi/teams/soranzo/projects/ALL_dB/csv_files/Vuckovic_VEP_scores.csv.gz")


outfile_binder_of_scores_CADD=$(echo "$output_dir""outfile""_""$type""_""$finemap_prob_Threshold"".out")
touch $outfile_binder_of_scores_CADD
echo -n "" > $outfile_binder_of_scores_CADD
name_binder_of_scores_CADD=$(echo "$type""_""$finemap_prob_Threshold""_job")


step_mem=$(expr $mem \* 2)
step_pc=$(expr $pc \* 2)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"

 

echo "bsub -G team151 -o $outfile_binder_of_scores_CADD -M $step_mem  -J $name_binder_of_scores_CADD -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_binder_of_scores_CADD \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--CADD_result $CADD_result \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output



echo "###################################################### CADD  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_CADD=/nfs/users/nfs_m/mt19/Scripts/R/341_CADD.R

CADD=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/CADD_GLOBAL.tsv')
type=$(echo "CADD")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")


outfile_CADD=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_CADD
echo -n "" > $outfile_CADD
name_CADD=$(echo "$type""_job")


step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"



echo "bsub -G team151 -o $outfile_CADD -M $step_mem -w\"done($name_binder_of_scores_CADD)\" -J $name_CADD -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CADD -M $step_mem  -J $name_CADD -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CADD \\" >> $output
echo "--CADD $CADD \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output



echo "###################################################### binder_of_scores_NCBoost  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_binder_of_scores_NCBoost=/nfs/users/nfs_m/mt19/Scripts/R/342_binder_NCBoost.R


ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
type=$(echo "binder_of_scores_NCBoost")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
NCBoost_result=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NCBoost/haemvar_manuel_200903_NCBOOST.tsv")


outfile_binder_of_scores_NCBoost=$(echo "$output_dir""outfile""_""$type""_""$finemap_prob_Threshold"".out")
touch $outfile_binder_of_scores_NCBoost
echo -n "" > $outfile_binder_of_scores_NCBoost
name_binder_of_scores_NCBoost=$(echo "$type""_""$finemap_prob_Threshold""_job")


step_mem=$(expr $mem \* 2)
step_pc=$(expr $pc \* 2)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"

 

echo "bsub -G team151 -o $outfile_binder_of_scores_NCBoost -M $step_mem  -J $name_binder_of_scores_NCBoost -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_binder_of_scores_NCBoost \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--NCBoost_result $NCBoost_result \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output

echo "###################################################### NCBoost  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_NCBoost=/nfs/users/nfs_m/mt19/Scripts/R/343_NCBoost.R

NCBoost=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/NCBoost_GLOBAL.tsv')
type=$(echo "NCBoost")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")


outfile_NCBoost=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_NCBoost
echo -n "" > $outfile_NCBoost
name_NCBoost=$(echo "$type""_job")


step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"



echo "bsub -G team151 -o $outfile_NCBoost -M $step_mem -w\"done($name_binder_of_scores_NCBoost)\" -J $name_NCBoost -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_NCBoost -M $step_mem  -J $name_NCBoost -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_NCBoost \\" >> $output
echo "--NCBoost $NCBoost \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output

echo "###################################################### binder_of_scores_constraint_Z  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_binder_of_scores_constraint_Z=/nfs/users/nfs_m/mt19/Scripts/R/349_Constraint_Z.R


ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
type=$(echo "binder_of_scores_constraint_Z")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
Constraint_Z=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Constraint_Z/constraint_z_genome_1kb.qc.download.txt")


outfile_binder_of_scores_constraint_Z=$(echo "$output_dir""outfile""_""$type""_""$finemap_prob_Threshold"".out")
touch $outfile_binder_of_scores_constraint_Z
echo -n "" > $outfile_binder_of_scores_constraint_Z
name_binder_of_scores_constraint_Z=$(echo "$type""_""$finemap_prob_Threshold""_job")


step_mem=$(expr $mem \* 2)
step_pc=$(expr $pc \* 2)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"

 

echo "bsub -G team151 -o $outfile_binder_of_scores_constraint_Z -M $step_mem  -J $name_binder_of_scores_constraint_Z -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_binder_of_scores_constraint_Z \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--Constraint_Z $Constraint_Z \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output




echo "###################################################### binder_of_scores_SpliceAI  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_binder_of_scores_SpliceAI=/nfs/users/nfs_m/mt19/Scripts/R/344_binder_SpliceAI.R


ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
type=$(echo "binder_of_scores_SpliceAI")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
SpliceAI_result=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/spliceAI/haemvar_manuel_200903_spliceAI_output_fixed.vcf")
GENE_table=$(echo "/nfs/users/nfs_m/mt19/RareVar_Dragana/Homo_sapiens.GRCh37.87_GENES_table.txt")

outfile_binder_of_scores_SpliceAI=$(echo "$output_dir""outfile""_""$type""_""$finemap_prob_Threshold"".out")
touch $outfile_binder_of_scores_SpliceAI
echo -n "" > $outfile_binder_of_scores_SpliceAI
name_binder_of_scores_SpliceAI=$(echo "$type""_""$finemap_prob_Threshold""_job")


step_mem=$(expr $mem \* 2)
step_pc=$(expr $pc \* 2)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"

 

echo "bsub -G team151 -o $outfile_binder_of_scores_SpliceAI -M $step_mem  -J $name_binder_of_scores_SpliceAI -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_binder_of_scores_SpliceAI \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--GENE_table $GENE_table \\" >> $output
echo "--SpliceAI_result $SpliceAI_result \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output



echo "###################################################### SpliceAI  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_SpliceAI=/nfs/users/nfs_m/mt19/Scripts/R/345_SpliceAI.R

SpliceAI=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/SpliceAI_GLOBAL.tsv')
type=$(echo "SpliceAI")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")


outfile_SpliceAI=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_SpliceAI
echo -n "" > $outfile_SpliceAI
name_SpliceAI=$(echo "$type""_job")


step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"



echo "bsub -G team151 -o $outfile_SpliceAI -M $step_mem -w\"done($name_binder_of_scores_SpliceAI)\" -J $name_SpliceAI -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_SpliceAI -M $step_mem  -J $name_SpliceAI -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_SpliceAI \\" >> $output
echo "--SpliceAI $SpliceAI \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output




echo "###################################################### binder_of_scores_unranked_chromstates  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_binder_of_scores_unranked_chromstates=/nfs/users/nfs_m/mt19/Scripts/R/329_binder_of_scores_unranked_chromstates_v2.R


ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
type=$(echo "binder_of_scores_unranked_chromstates")
excluded_phenotypes=$(echo "wbc,eo_p,mono_p,neut_p,lymph_p,baso_p")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
relevant_not_relevant_weights=$(echo "1,0.01")
matrix_weighted_regulatory_states=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Weights_chromatin_states.tsv")
chromstates_INITIAL=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/chrmstates_graphs.csv")
Trait_to_CT_table=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/ALL_dB/CellType_Trait_table_generation.txt")



outfile_binder_of_scores_unranked_chromstates=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_binder_of_scores_unranked_chromstates
echo -n "" > $outfile_binder_of_scores_unranked_chromstates
name_binder_of_scores_unranked_chromstates=$(echo "$type""_job")


step_mem=$(expr $mem \* 2)
step_pc=$(expr $pc \* 2)
step_queue=$(echo "normal")

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"
echo "$queue""->""$step_queue"

 

echo "bsub -G team151 -o $outfile_binder_of_scores_unranked_chromstates -M $step_mem  -J $name_binder_of_scores_unranked_chromstates -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $step_queue -- \\" >> $output
echo "\"$Rscript $Rscript_binder_of_scores_unranked_chromstates \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--matrix_weighted_regulatory_states $matrix_weighted_regulatory_states \\" >> $output
echo "--relevant_not_relevant_weights $relevant_not_relevant_weights \\" >> $output
echo "--chromstates_INITIAL $chromstates_INITIAL \\" >> $output
echo "--Trait_to_CT_table $Trait_to_CT_table \\" >> $output
echo "--excluded_phenotypes $excluded_phenotypes \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output




echo "###################################################### chromstates Rank  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_chromstates=/nfs/users/nfs_m/mt19/Scripts/R/337_Rank_chromstates_v2.R

chromstates_pre_ranked=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/chromstates_GLOBAL_preranked.tsv')
type=$(echo "chromstates")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
desiR_weights=$(echo "0.001,1.2,1")



outfile_chromstates=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_chromstates
echo -n "" > $outfile_chromstates
name_chromstates=$(echo "$type""_job")


step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"


echo "bsub -G team151 -o $outfile_chromstates -M $step_mem -w\"done($name_binder_of_scores_unranked_chromstates)\" -J $name_chromstates -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_chromstates -M $step_mem  -J $name_chromstates -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_chromstates \\" >> $output
echo "--chromstates_pre_ranked $chromstates_pre_ranked \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--desiR_weights $desiR_weights \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output


echo "###################################################### binder_of_scores_unranked_PCHiC  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_binder_of_scores_unranked_PCHiC=/nfs/users/nfs_m/mt19/Scripts/R/329_binder_of_scores_unranked_but_thresholded_PCHiC_v2.R


ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
type=$(echo "binder_of_scores_unranked_PCHiC")
excluded_phenotypes=$(echo "wbc,eo_p,mono_p,neut_p,lymph_p,baso_p")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
relevant_not_relevant_weights=$(echo "1,0.01")
PCHiC_INITIAL=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/PCHIC_ChicagoScore_graphs.csv")
Trait_to_CT_table=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/ALL_dB/PCHIC_part_II_Javierre_corresp_generation.txt")



outfile_binder_of_scores_unranked_PCHiC=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_binder_of_scores_unranked_PCHiC
echo -n "" > $outfile_binder_of_scores_unranked_PCHiC
name_binder_of_scores_unranked_PCHiC=$(echo "$type""_job")


step_mem=$(expr $mem \* 2)
step_pc=$(expr $pc \* 2)
step_queue=$(echo "normal")

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"
echo "$queue""->""$step_queue"

 

echo "bsub -G team151 -o $outfile_binder_of_scores_unranked_PCHiC -M $step_mem  -J $name_binder_of_scores_unranked_PCHiC -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $step_queue -- \\" >> $output
echo "\"$Rscript $Rscript_binder_of_scores_unranked_PCHiC \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--relevant_not_relevant_weights $relevant_not_relevant_weights \\" >> $output
echo "--PCHiC_INITIAL $PCHiC_INITIAL \\" >> $output
echo "--Trait_to_CT_table $Trait_to_CT_table \\" >> $output
echo "--excluded_phenotypes $excluded_phenotypes \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output




echo "###################################################### PCHiC Rank  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_PCHiC=/nfs/users/nfs_m/mt19/Scripts/R/336_Rank_PCHiC_v2.R

PCHiC_pre_ranked=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/PCHiC_GLOBAL_preranked.tsv')
type=$(echo "PCHiC")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
desiR_weights=$(echo "0.1,10,1")



outfile_PCHiC=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_PCHiC
echo -n "" > $outfile_PCHiC
name_PCHiC=$(echo "$type""_job")


step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"


echo "bsub -G team151 -o $outfile_PCHiC -M $step_mem -w\"done($name_binder_of_scores_unranked_PCHiC)\" -J $name_PCHiC -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_PCHiC -M $step_mem  -J $name_PCHiC -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_PCHiC \\" >> $output
echo "--PCHiC_pre_ranked $PCHiC_pre_ranked \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--desiR_weights $desiR_weights \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output


echo "###################################################### binder_of_scores_unranked_ATAC  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_binder_of_scores_unranked_ATAC=/nfs/users/nfs_m/mt19/Scripts/R/329_binder_of_scores_unranked_ATAC_v2.R


ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
type=$(echo "binder_of_scores_unranked_ATAC")
excluded_phenotypes=$(echo "wbc,eo_p,mono_p,neut_p,lymph_p,baso_p")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
relevant_not_relevant_weights=$(echo "1,0.01")
ATAC_INITIAL=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/ATAC_Seq_Ranked_graphs.csv")
Trait_to_Lineage_table=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/ALL_dB/ATAC_scaled_trait_table_generation.txt")
Lineage_to_Cell_table=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/ALL_dB/ATAC_scaled_Lineage_hierarchy_generation.txt")



outfile_binder_of_scores_unranked_ATAC=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_binder_of_scores_unranked_ATAC
echo -n "" > $outfile_binder_of_scores_unranked_ATAC
name_binder_of_scores_unranked_ATAC=$(echo "$type""_job")


step_mem=$(expr $mem \* 2)
step_pc=$(expr $pc \* 2)
step_queue=$(echo "normal")

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"
echo "$queue""->""$step_queue"

 

echo "bsub -G team151 -o $outfile_binder_of_scores_unranked_ATAC -M $step_mem  -J $name_binder_of_scores_unranked_ATAC -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $step_queue -- \\" >> $output
echo "\"$Rscript $Rscript_binder_of_scores_unranked_ATAC \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--relevant_not_relevant_weights $relevant_not_relevant_weights \\" >> $output
echo "--ATAC_INITIAL $ATAC_INITIAL \\" >> $output
echo "--Trait_to_Lineage_table $Trait_to_Lineage_table \\" >> $output
echo "--Lineage_to_Cell_table $Lineage_to_Cell_table \\" >> $output
echo "--excluded_phenotypes $excluded_phenotypes \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output


 
echo "###################################################### ATAC_Rank  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_ATAC_Rank=/nfs/users/nfs_m/mt19/Scripts/R/330_Rank_ATAC_v2.R

ATAC_pre_ranked=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/ATAC_GLOBAL_preranked.tsv')
type=$(echo "ATAC_Rank")
Open_in_CT_threshold=$(echo "0.1")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
desiR_weights=$(echo "0.5,5,0.05,0.8,1,3")


outfile_ATAC_Rank=$(echo "$output_dir""outfile""_""$type""_""$Open_in_CT_threshold"".out")
touch $outfile_ATAC_Rank
echo -n "" > $outfile_ATAC_Rank
name_ATAC_Rank=$(echo "$type""_""$Open_in_CT_threshold""_job")


step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"


echo "bsub -G team151 -o $outfile_ATAC_Rank -M $step_mem -w\"done($name_binder_of_scores_unranked_ATAC)\"  -J $name_ATAC_Rank -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_ATAC_Rank -M $step_mem  -J $name_ATAC_Rank -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_ATAC_Rank \\" >> $output
echo "--ATAC_pre_ranked $ATAC_pre_ranked \\" >> $output
echo "--Open_in_CT_threshold $Open_in_CT_threshold \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--desiR_weights $desiR_weights \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output


echo "###################################################### multi_lineage_ATAC  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_multi_lineage_ATAC=/nfs/users/nfs_m/mt19/Scripts/R/331_Rank_multi_lineage_ATAC.R

ATAC_ranked=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/ATAC_GLOBAL_Ranked.tsv')
type=$(echo "multi_lineage_ATAC")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
desiR_weights=$(echo "0.9,2.5,0.1,0.7,3,1")


outfile_multi_lineage_ATAC=$(echo "$output_dir""outfile""_""$type""_""$Open_in_CT_threshold"".out")
touch $outfile_multi_lineage_ATAC
echo -n "" > $outfile_multi_lineage_ATAC
name_multi_lineage_ATAC=$(echo "$type""_""$Open_in_CT_threshold""_job")


step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"


echo "bsub -G team151 -o $outfile_multi_lineage_ATAC -M $step_mem -w\"done($name_ATAC_Rank)\" -J $name_multi_lineage_ATAC -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_multi_lineage_ATAC -M $step_mem  -J $name_multi_lineage_ATAC -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_multi_lineage_ATAC \\" >> $output
echo "--ATAC_ranked $ATAC_ranked \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--desiR_weights $desiR_weights \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output



echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "####################################################### GENE BASED FEATURES #############################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output


echo "###################################################### binder_of_scores_GENE_EXP  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_binder_of_scores_GENE_EXP=/nfs/users/nfs_m/mt19/Scripts/R/338_binder_of_scores_Gene_EXP.R


ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
type=$(echo "binder_of_scores_GENE_EXP")
excluded_phenotypes=$(echo "wbc,eo_p,mono_p,neut_p,lymph_p,baso_p")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
VEP_CSQ=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/VEP_consequence_graphs.csv')
PCHiC=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/PCHIC_ChicagoScore_graphs.csv')
GENE_EXP=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/BP_Iso_Reanalysis/OLD/GENE_EXP_Harmonization.tsv")
TOME_correspondence=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/ALL_dB/Correspondence_phenotype_TOME.txt")





outfile_binder_of_scores_GENE_EXP=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_binder_of_scores_GENE_EXP
echo -n "" > $outfile_binder_of_scores_GENE_EXP
name_binder_of_scores_GENE_EXP=$(echo "$type""_job")


step_mem=$(expr $mem \* 4)
step_pc=$(expr $pc \* 4)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"

 

echo "bsub -G team151 -o $outfile_binder_of_scores_GENE_EXP -M $step_mem  -J $name_binder_of_scores_GENE_EXP -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_binder_of_scores_GENE_EXP \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--VEP_CSQ $VEP_CSQ \\" >> $output
echo "--PCHiC $PCHiC \\" >> $output
echo "--GENE_EXP $GENE_EXP \\" >> $output
echo "--TOME_correspondence $TOME_correspondence \\" >> $output
echo "--excluded_phenotypes $excluded_phenotypes \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output

echo "###################################################### GENE_EXP_Rank  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_GENE_EXP_Rank=/nfs/users/nfs_m/mt19/Scripts/R/339_Gene_EXP.R

GENE_EXP_pre_ranked=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/GENE_EXP_GLOBAL.tsv')
type=$(echo "GENE_EXP_Rank")
relevant_not_relevant_weights=$(echo "1,0.01")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
desiR_weights=$(echo "0.1,20,1")


outfile_GENE_EXP_Rank=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_GENE_EXP_Rank
echo -n "" > $outfile_GENE_EXP_Rank
name_GENE_EXP_Rank=$(echo "$type""_job")


step_mem=$(expr $mem \* 3)
step_pc=$(expr $pc \* 3)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"


echo "bsub -G team151 -o $outfile_GENE_EXP_Rank -M $step_mem -w\"done($name_binder_of_scores_GENE_EXP)\"  -J $name_GENE_EXP_Rank -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_GENE_EXP_Rank -M $step_mem  -J $name_GENE_EXP_Rank -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_GENE_EXP_Rank \\" >> $output
echo "--GENE_EXP_pre_ranked $GENE_EXP_pre_ranked \\" >> $output
echo "--relevant_not_relevant_weights $relevant_not_relevant_weights \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--desiR_weights $desiR_weights \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output

echo "###################################################### binder_of_scores_GENE_PLOEUF  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_binder_of_scores_GENE_PLOEUF=/nfs/users/nfs_m/mt19/Scripts/R/350_binder_of_scores_pLOEUF.R


type=$(echo "binder_of_scores_GENE_PLOEUF")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
VEP_CSQ=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/VEP_consequence_graphs.csv')
PCHiC=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/PCHIC_ChicagoScore_graphs.csv')
GENE_PLOEUF=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/pLOEUF/gnomad.v2.1.1.lof_metrics.by_gene.txt")
TRANSCRIPTS_table=$(echo "/nfs/users/nfs_m/mt19/RareVar_Dragana/Homo_sapiens.GRCh37.87_Transcripts_table.txt")



outfile_binder_of_scores_GENE_PLOEUF=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_binder_of_scores_GENE_PLOEUF
echo -n "" > $outfile_binder_of_scores_GENE_PLOEUF
name_binder_of_scores_GENE_PLOEUF=$(echo "$type""_job")


step_mem=$(expr $mem \* 2)
step_pc=$(expr $pc \* 2)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"

 

echo "bsub -G team151 -o $outfile_binder_of_scores_GENE_PLOEUF -M $step_mem  -J $name_binder_of_scores_GENE_PLOEUF -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_binder_of_scores_GENE_PLOEUF \\" >> $output
echo "--VEP_CSQ $VEP_CSQ \\" >> $output
echo "--PCHiC $PCHiC \\" >> $output
echo "--GENE_PLOEUF $GENE_PLOEUF \\" >> $output
echo "--TRANSCRIPTS_table $TRANSCRIPTS_table \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output



echo "###################################################### binder_of_scores_COGS  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_binder_of_scores_COGS=/nfs/users/nfs_m/mt19/Scripts/R/334_binder_of_scores_COGS.R


ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
type=$(echo "binder_of_scores_COGS")
excluded_phenotypes=$(echo "wbc,eo_p,mono_p,neut_p,lymph_p,baso_p")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
VEP_CSQ=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/VEP_consequence_graphs.csv')
PCHiC=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/PCHIC_ChicagoScore_graphs.csv')


outfile_binder_of_scores_COGS=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_binder_of_scores_COGS
echo -n "" > $outfile_binder_of_scores_COGS
name_binder_of_scores_COGS=$(echo "$type""_job")


step_mem=$(expr $mem \* 2)
step_pc=$(expr $pc \* 2)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"

 

echo "bsub -G team151 -o $outfile_binder_of_scores_COGS -M $step_mem  -J $name_binder_of_scores_COGS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_binder_of_scores_COGS \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--VEP_CSQ $VEP_CSQ \\" >> $output
echo "--PCHiC $PCHiC \\" >> $output
echo "--excluded_phenotypes $excluded_phenotypes \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output

echo "###################################################### COGS  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_COGS=/nfs/users/nfs_m/mt19/Scripts/R/335_COGS_v2.R

COGS=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/COGS_GLOBAL.tsv')
type=$(echo "COGS")
COGS_Threshold=$(echo '0')
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")
#desiR_weights=$(echo "0.5,10,0.05,0.8,1,3")


outfile_COGS=$(echo "$output_dir""outfile""_""$type""_""$Open_in_CT_threshold"".out")
touch $outfile_COGS
echo -n "" > $outfile_COGS
name_COGS=$(echo "$type""_""$Open_in_CT_threshold""_job")


step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"



echo "bsub -G team151 -o $outfile_COGS -M $step_mem -w\"done($name_binder_of_scores_COGS)\" -J $name_COGS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_COGS -M $step_mem  -J $name_COGS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_COGS \\" >> $output
echo "--COGS $COGS \\" >> $output
echo "--COGS_Threshold $COGS_Threshold \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output



echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "####################################################### Comparison with pathogenic variants from Vuckovic et al #########################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output


echo "###################################################### binder_of_abs_finemap_for_Pathogenic_variant_comparison  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output



finemap_prob_Threshold_string=$(echo '0;''0.5;''0.9;')
#echo $finemap_prob_Threshold_string
declare -a c=($(echo "$finemap_prob_Threshold_string" | tr ";" '\n'))

arraylength=${#c[@]}



echo "array_length: $arraylength"

 
for (( i=0; i<${arraylength}; i++ ));
do
    finemap_prob_Threshold=${c[$i]}
    echo "index: $i, finemap_prob_Threshold: $finemap_prob_Threshold"

    output_dir_Pathogenic=$(echo "$output_dir""Pathogenic_variant_comparison""_""$finemap_prob_Threshold""/")

    echo "$output_dir_Pathogenic"
    
    rm -rf $output_dir_Pathogenic
    mkdir -p $output_dir_Pathogenic



    Rscript_binder_of_abs_finemap_for_Pathogenic_variant_comparison=/nfs/users/nfs_m/mt19/Scripts/R/353_New_binder_for_Absolute_effect_size_for_Pathogenic_variants_comparison.R


    ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
    type=$(echo "binder_of_abs_finemap_for_Pathogenic_variant_comparison")



    outfile_binder_of_abs_finemap_for_Pathogenic_variant_comparison=$(echo "$output_dir_Pathogenic""outfile""_""$type""_""$finemap_prob_Threshold"".out")
    touch $outfile_binder_of_abs_finemap_for_Pathogenic_variant_comparison
    echo -n "" > $outfile_binder_of_abs_finemap_for_Pathogenic_variant_comparison
    name_binder_of_abs_finemap_for_Pathogenic_variant_comparison=$(echo "$type""_""$finemap_prob_Threshold""_job")

    
    step_mem=$(expr $mem \* 1)
    step_pc=$(expr $pc \* 1)

    echo "$mem""->""$step_mem"
    echo "$pc""->""$step_pc"

    

    echo "bsub -G team151 -o $outfile_binder_of_abs_finemap_for_Pathogenic_variant_comparison -M $step_mem  -J $name_binder_of_abs_finemap_for_Pathogenic_variant_comparison -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
    echo "\"$Rscript $Rscript_binder_of_abs_finemap_for_Pathogenic_variant_comparison \\" >> $output
    echo "--ALL_dB $ALL_dB \\" >> $output
    echo "--finemap_prob_Threshold $finemap_prob_Threshold \\" >> $output
    echo "--type $type --out $output_dir_Pathogenic\"" >> $output




    Rscript_compare_pathogenic_vs_prioritised_variants=/nfs/users/nfs_m/mt19/Scripts/R/352_Comparison_between_pathogenic_and_prioritised_variants.R

    List_of_pathogenic_variants=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Pathogenic_variants.txt")
    VAR_Prioritization_dB=$(echo "$output_dir""ER_Labelling_Initial_Selection.rds")
    Table_S4=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/Fig4_pannels/Supp_Table_4_CURATED_Plus_phenotypes.rds")
    RMV_common=$(echo "chr4_1008212_C_T")
    RMV_labels=$(echo 'NOT_SCREENED_MPRA,NOT_SCREENED_genIE,No_RNA_Seq_HET_carriers')


    type=$(echo "compare_pathogenic_vs_prioritised_variants")




    outfile_compare_pathogenic_vs_prioritised_variants=$(echo "$output_dir_Pathogenic""outfile""_""$type""_""$finemap_prob_Threshold"".out")
    touch $outfile_compare_pathogenic_vs_prioritised_variants
    echo -n "" > $outfile_compare_pathogenic_vs_prioritised_variants
    name_compare_pathogenic_vs_prioritised_variants=$(echo "$type""_""$finemap_prob_Threshold""_job")


    step_mem=$(expr $mem \* 1)
    step_pc=$(expr $pc \* 1)

    echo "$mem""->""$step_mem"
    echo "$pc""->""$step_pc"
    

    echo "bsub -G team151 -o $outfile_compare_pathogenic_vs_prioritised_variants -M $step_mem -w\"done($name_ER_Labelling) && done($name_binder_of_abs_finemap_for_Pathogenic_variant_comparison)\" -J $name_compare_pathogenic_vs_prioritised_variants -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
    #echo "bsub -G team151 -o $outfile_compare_pathogenic_vs_prioritised_variants -M $step_mem  -J $name_compare_pathogenic_vs_prioritised_variants -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
    echo "\"$Rscript $Rscript_compare_pathogenic_vs_prioritised_variants \\" >> $output
    echo "--List_of_pathogenic_variants $List_of_pathogenic_variants \\" >> $output
    echo "--finemap_prob_Threshold $finemap_prob_Threshold \\" >> $output
    echo "--VAR_Prioritization_dB $VAR_Prioritization_dB \\" >> $output
    echo "--Table_S4 $Table_S4 \\" >> $output
    echo "--RMV_common $RMV_common \\" >> $output
    echo "--RMV_labels $RMV_labels \\" >> $output
    echo "--type $type --out $output_dir_Pathogenic\"" >> $output


done






echo  "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "####################################################### DATA WRANGLING, STATS & GRAPHs  #################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

echo "###################################################### Prepared_files_and_add_annotation_layer  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_Prepared_files_and_add_annotation_layer=/nfs/users/nfs_m/mt19/Scripts/R/332_binder_of_Prepared_files_and_add_annotation_layer.R

ATAC_data=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Prepared_file_ATAC.rds")
multi_ATAC_data=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Prepared_file_multi_lineage_ATAC.rds")
COGS=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Prepared_file_COGS.rds")
oe_LOF=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Prepared_file_GENE_PLOEUF.rds")
GENE_EXP=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Prepared_file_GENE_EXP.rds")
PCHiC=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Prepared_file_PCHiC.rds")
chromstates=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Prepared_file_chromstates.rds")
CADD=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Prepared_file_CADD.rds")
Constraint_Z=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Prepared_file_Constraint_Z.rds")
NCBoost=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Prepared_file_NCBoost.rds")
SpliceAI=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Prepared_file_SpliceAI.rds")
GWAS_GLOBAL_per_traits=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/GWAS_GLOBAL_per_traits.tsv")
MAF_GLOBAL=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/MAF_GLOBAL.tsv")


VAR_Prioritization_dB=$(echo "$output_dir""ER_Labelling_Initial_Selection.rds")
Table_S6=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/Table_S6_Manual_curation.rds")


type=$(echo "Prepared_files_and_add_annotation_layer")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")



outfile_Prepared_files_and_add_annotation_layer=$(echo "$output_dir""outfile""_""$type""_""$Open_in_CT_threshold"".out")
touch $outfile_Prepared_files_and_add_annotation_layer
echo -n "" > $outfile_Prepared_files_and_add_annotation_layer
name_Prepared_files_and_add_annotation_layer=$(echo "$type""_""$Open_in_CT_threshold""_job")


step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"
 

echo "bsub -G team151 -o $outfile_Prepared_files_and_add_annotation_layer -M $step_mem -w\"done($name_ER_Labelling) && done($name_binder_of_scores_GWAS) && done($name_binder_of_scores_GENE_PLOEUF) && done($name_binder_of_scores_constraint_Z) && done($name_ATAC_Rank) && done($name_multi_lineage_ATAC) && done($name_SpliceAI) && done($name_CADD) && done($name_NCBoost) && done($name_COGS) && done($name_PCHiC) && done($name_chromstates) && done($name_GENE_EXP_Rank)\" -J $name_Prepared_files_and_add_annotation_layer -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_Prepared_files_and_add_annotation_layer -M $step_mem -w\"done($name_COGS)\" -J $name_Prepared_files_and_add_annotation_layer -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_Prepared_files_and_add_annotation_layer -M $step_mem  -J $name_Prepared_files_and_add_annotation_layer -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_Prepared_files_and_add_annotation_layer \\" >> $output
echo "--ATAC_data $ATAC_data \\" >> $output
echo "--multi_ATAC_data $multi_ATAC_data \\" >> $output
echo "--Constraint_Z $Constraint_Z \\" >> $output
echo "--NCBoost $NCBoost \\" >> $output
echo "--SpliceAI $SpliceAI \\" >> $output
echo "--GWAS_GLOBAL_per_traits $GWAS_GLOBAL_per_traits \\" >> $output
echo "--MAF_GLOBAL $MAF_GLOBAL \\" >> $output
echo "--CADD $CADD \\" >> $output
echo "--COGS $COGS \\" >> $output
echo "--oe_LOF $oe_LOF \\" >> $output
echo "--GENE_EXP $GENE_EXP \\" >> $output
echo "--PCHiC $PCHiC \\" >> $output
echo "--chromstates $chromstates \\" >> $output
echo "--VAR_Prioritization_dB $VAR_Prioritization_dB \\" >> $output
echo "--Table_S6 $Table_S6 \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output



echo "###################################################### Stats_Z_score  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_Stats_Z_score=/nfs/users/nfs_m/mt19/Scripts/R/333_stats_on_Z_score.R

Master_file=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Master_file_scores.rds")
Table_of_labels=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Table_of_labels.rds")



type=$(echo "Stats_Z_score")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")


 
outfile_Stats_Z_score=$(echo "$output_dir""outfile""_""$type""_""$Open_in_CT_threshold"".out")
touch $outfile_Stats_Z_score
echo -n "" > $outfile_Stats_Z_score
name_Stats_Z_score=$(echo "$type""_""$Open_in_CT_threshold""_job")


step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"


echo "bsub -G team151 -o $outfile_Stats_Z_score -M $step_mem -w\"done($name_Prepared_files_and_add_annotation_layer)\" -J $name_Stats_Z_score -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_Stats_Z_score -M $step_mem  -J $name_Stats_Z_score -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_Stats_Z_score \\" >> $output
echo "--Master_file $Master_file \\" >> $output
echo "--Table_of_labels $Table_of_labels \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output




echo "###################################################### Z_score_violin_plots  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_Stats_Z_score=/nfs/users/nfs_m/mt19/Scripts/R/348_Z_score_violin_plots_v2.R

Master_file=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Master_file_scores.rds")
Table_of_labels=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/Table_of_labels.rds")
Stats_file=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/STATS_value_Z_score.rds")
Categories_colors=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/new_desiR_scores/ER_Labelling_Categories_colors.rds")
subset_variables_1=$(echo "PP,Absolute_effect_size,CADD_raw,Constraint_Z,NCBoost,SpliceAI_DG,SpliceAI_DL,SpliceAI_AG,SpliceAI_AL,Rank_ATAC_erythroid_lineage,Rank_ATAC_mega_lineage,Rank_ATAC_gran_mono_lineage,Rank_ATAC_lymph_lineage,multi_lineage_ATAC,Rank_PCHiC,Rank_chromstates")
subset_variables_2=$(echo "COGS,oe_lof,Rank_GENE_EXP")

type=$(echo "Z_score_violin_plots")
tracking_variants=$(echo "chr1_202129205_G_A,chr12_111844956_C_T,chr18_60880701_T_C,chr18_60920854_C_T,chr3_128317978_C_T,chr3_128322617_G_A,chr3_184091102_T_G,chr17_38764524_T_A,chr3_71355240_G_C,chr16_86016328_C_T")


 
outfile_Z_score_violin_plots=$(echo "$output_dir""outfile""_""$type""_""$Open_in_CT_threshold"".out")
touch $outfile_Z_score_violin_plots
echo -n "" > $outfile_Z_score_violin_plots
name_Z_score_violin_plots=$(echo "$type""_""$Open_in_CT_threshold""_job")


step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"


echo "bsub -G team151 -o $outfile_Z_score_violin_plots -M $step_mem -w\"done($name_Stats_Z_score)\" -J $name_Z_score_violin_plots -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_Z_score_violin_plots -M $step_mem  -J $name_Z_score_violin_plots -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_Stats_Z_score \\" >> $output
echo "--Master_file $Master_file \\" >> $output
echo "--Table_of_labels $Table_of_labels \\" >> $output
echo "--subset_variables_1 $subset_variables_1 \\" >> $output
echo "--subset_variables_2 $subset_variables_2 \\" >> $output
echo "--Stats_file $Stats_file \\" >> $output
echo "--Categories_colors $Categories_colors \\" >> $output
echo "--tracking_variants $tracking_variants \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output


bash $output
 
