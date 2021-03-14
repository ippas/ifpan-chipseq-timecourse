
#!/bin/bash
#file bed
#GSM2421929_ENCFF835HHK_peaks_GRCh38.bed.gz 
#GSM2421930_ENCFF044MLR_peaks_GRCh38.bed.gz 
#GSM2421931_ENCFF597LEE_peaks_GRCh38.bed.gz

ACCESS_PATH_1=~/ChIP-seq/EXTRACT
ACCESS_PATH_2=~/ChIP-seq/DATA
ACCESS_PATH_3=~/ifpan-chipseq-timecourse/DATA

FILE_1=$ACCESS_PATH_3/peaks_part_of_common_three_file_NR3C1_time60.bed
FILE_2=$ACCESS_PATH_3/tmp_range_all_genes.bed
FILE_3=$ACCESS_PATH_2/tmp_file1_file2.bed


bedtools intersect -a $ACCESS_PATH_1/GSM2421929_ENCFF835HHK_peaks_GRCh38.bed.gz -b $ACCESS_PATH_1/GSM2421930_ENCFF044MLR_peaks_GRCh38.bed.gz  > $FILE_3
bedtools intersect -a $FILE_3 -b $ACCESS_PATH_1/GSM2421931_ENCFF597LEE_peaks_GRCh38.bed.gz > $FILE_1


cat $ACCESS_PATH_3/range_all_genes.bed |
   tail +2 |
   awk '{print "chr"$3"\t"$4"\t"$5"\t"$0}' > $FILE_2


bedtools intersect -a $FILE_2 -b $FILE_1 -wo |
   cut -f4-12 |
   awk 'OFS="\t" {print $1, $2, $3,$8, $9, $6}' |
   sed '1 i\ensemblid\tgene.name\tchromosome\tstart_peak\tend_peak\tgene.regulation' > $ACCESS_PATH_3/peaks_all_genes.tsv

echo "start extract data from chipseq"

~/ifpan-chipseq-timecourse/SCRIPTS/./bigwig_genomic_amplitude_extract_normalize_to_tsv.NR3C1-EP300.sh ~/ifpan-chipseq-timecourse/DATA/peaks_all_genes.tsv > ~/ChIP-seq/DATA/peaks_all_genes_ep300_nr3c1_amplitude.tsv


#delecting unnecessary files
rm $FILE_1
rm $FILE_2
rm $FILE_3
rm $ACCESS_PATH_3/peaks_all_genes.tsv


