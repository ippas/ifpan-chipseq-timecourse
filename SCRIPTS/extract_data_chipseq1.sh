#!/bin/bash

ACCESS_PATH_1=~/ChIP-seq/EXTRACT
ACCESS_PATH_2=~/ChIP-seq/DATA
ACCESS_PATH_3=~/ifpan-chipseq-timecourse/DATA
ACCESS_PATH_4=~/ifpan-chipseq-timecourse/SCRIPTS


############################
# extract range promotores #
############################
echo "extract range promotores"

$ACCESS_PATH_4/./bigwig_genomic_bucket500_extract_normalize_to_tsv.sh $ACCESS_PATH_3/promotores_peaks_info.tsv > $ACCESS_PATH_2/promotores_peaks_value.tsv


##########################################################
# prepare file BED with the same peaks from three sample #
##########################################################

#GSM2421929_ENCFF835HHK_peaks_GRCh38.bed.gz
#GSM2421930_ENCFF044MLR_peaks_GRCh38.bed.gz
#GSM2421931_ENCFF597LEE_peaks_GRCh38.bed.gz


FILE_1=$ACCESS_PATH_3/peaks_part_of_common_three_file_NR3C1_time60.bed
FILE_2=$ACCESS_PATH_3/tmp_significant_random.bed
FILE_3=$ACCESS_PATH_2/tmp_file1_file2.bed
FILE_4=$ACCESS_PATH_3/enhancer_peaks_info.tsv


bedtools intersect -a $ACCESS_PATH_1/GSM2421929_ENCFF835HHK_peaks_GRCh38.bed.gz -b $ACCESS_PATH_1/GSM2421930_ENCFF044MLR_peaks_GRCh38.bed.gz  > $FILE_3
bedtools intersect -a $FILE_3 -b $ACCESS_PATH_1/GSM2421931_ENCFF597LEE_peaks_GRCh38.bed.gz > $FILE_1


cat $ACCESS_PATH_3/enhancer_info.tsv |
   tail +2 |
   awk '{print "chr"$3"\t"$4"\t"$5"\t"$0}' |
   sed 's/\t-[0-9]*/\t0/g' > $FILE_2


bedtools intersect -a $FILE_2 -b $FILE_1 -wo |
   cut -f4-12 |
   awk 'OFS="\t" {print $1, $2, $3,$8, $9, $6}' |
   sed '1 i\ensemblid\tgene.name\tchromosome\tstart_peak\tend_peak\tgene.regulation' > $FILE_4

rm $FILE_1
rm $FILE_2
rm $FILE_3


####################################
# extract amplitudes for enhancers #
####################################

echo "extract amplitudes for enhancers"

$ACCESS_PATH_4/./bigwig_genomic_amplitude_extract_normalize_to_tsv.sh $FILE_4 > $ACCESS_PATH_2/enhancer_amplitude_value.tsv


###########################
# extract range enhancers #
###########################

echo "extract range enhancers"

$ACCESS_PATH_4/./bigwig_genomic_range_extract_normalize_to_tsv.sh $FILE_4 > $ACCESS_PATH_2/enhancer_peaks_value.tsv



rm $ACCESS_PATH_3/$FILE_4

