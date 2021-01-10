#!/bin/bash

#Skript make normalization bigwig file and extract the reading form chip-seq and also bucket values

# imput file should contein in columns:
#  1. ensemblid
#  2. gene name
#  3. chromosome
#  4. start - position from which start extract from BedGraph
#  5. end - positon from which end extract from BedGraph
#  6. gene.regulation - information about gene regulation

# output file contein columns:
#  1. gene name
#  2. chromosome
#  3. start - position from which start extract from BedGraph
#  4. end - positon from which end extract from BedGraph
#  5. gene.regulation - information about gene regulation
#  6. TF - transcript factor with which the cells were treated
#  7. time - time in which cells were consolidate
#  8. file - name file from which the readings chip-seq were extract
#  9-48. values the readings chip-seq which are bucket after 500, last bucket after 501

function get_coverage {
    CHROMOSOME=$1
    START=$2
    END=$3
    NORMALIZE_FILE=$4
    $ACCESS_PATH_3/./bigWigToBedGraph -chrom=$CHROMOSOME -start=$START -end=$END $NORMALIZE_FILE stdout |
       awk -v OFS="\t" '{for (i = $2; i < $3; i++) print $1"\t"i"\t"(i+1)"\t"$4"\t"}' | cut -f4 | awk -v RS= -v OFS='\t' '{$1 = $1} 1' |
       sed 's/\t/\n/500;P;D' |
       tr '\n' ':' |
       sed 's/:/\t/40' |
       sed 's/:/\n/g' |
       awk '{for(i=1;i<=500;i++){NUM=NUM?NUM+$i:$i};$(501+1)=NUM;NUM=""} 1' |
       awk '{print $NF}' |
       tr "\n" "\t" |
       awk '{ for(i=1; i<NF; ++i) printf $i/500 FS; print $NF/501}' | 
       sed 's/ /\t/g'
       #awk '{ for(i=1; i<NF; ++i) printf $i/500"\t"FS; print """\t"$NF/501}' |
       #sed 's/\t \t/\t/g' |
       #sed 's/\t /\t/g'
       #awk -v OFS="\t" '{for (i = $2; i < $3; i++) print $1"\t"i"\t"(i+1)"\t"$4"\t"}' | cut -f4 | awk -v RS= -v OFS='\t' '{$1 = $1} 1'
}

function get_coverages_for_file {
    GENES_INFO_FILE_NAME=$1
    NORMALIZE_FILE=$2
    cat $GENES_INFO_FILE_NAME | tail -n +2 | while read GENE_INFO
    do
       CHROMOSOME="chr"$(echo $GENE_INFO | awk '{print $3}')
       START=$(echo $GENE_INFO | awk '{print $4}')
       END=$(echo $GENE_INFO | awk '{print $5}')
       SYMBOL=$(echo $GENE_INFO | awk '{print $2}')
       GENE_REGULATION=$(echo $GENE_INFO | awk '{print $6}')
       PREFIX=$SYMBOL"\t"$CHROMOSOME"\t"$START"\t"$END"\t"$GENE_REGULATION"\t"$TF_NAME"\t"$TIME"\t"$FILE_NAME"\t"
       get_coverage $CHROMOSOME $START $END $NORMALIZE_FILE | awk -v prefix="$PREFIX" '{print prefix $0}'
    done
}
INPUT_FILE=$1
ACCESS_PATH_1=~/ChIP-seq/EXTRACT
ACCESS_PATH_2=~/ifpan-chipseq-timecourse/DATA
ACCESS_PATH_3=~/ifpan-chipseq-timecourse/SCRIPTS
ACCESS_PATH_4=~/ChIP-seq/DOWNLOAD

cat $ACCESS_PATH_2/chipseq-file-info.tsv | grep -v Control | grep -P 'NR3C1|EP300|H3K27ac|H3K4me1' |  while read LINE #grep -P 'NR3C1|EP300|H3K27ac|H3K4me1' it' to remove
do
  TF_NAME=$(echo $LINE | cut -d $' ' -f 1)
  TIME=$(echo $LINE | cut -d $' ' -f 2)
  GSE=$(echo $LINE | cut -d $' ' -f 4)
  TAR=$ACCESS_PATH_4/$GSE"_RAW.tar"
  REPLICATES=$(tar -tf $TAR | cut -d '_' -f 1 | uniq)
  for REPLICATE in $REPLICATES
  do
    BIGWIG_FILE=$(ls $ACCESS_PATH_1/$REPLICATE*Wig)
    FILE_NAME=$(echo $BIGWIG_FILE | awk -F "/" '{print $6}')
    #bigwig file normalization
    #convert bigwig to bedgraph
    $ACCESS_PATH_3/./bigWigToBedGraph $BIGWIG_FILE /dev/fd/1 > $ACCESS_PATH_1/tmp.unnormalize.$REPLICATE.bgd
    #normalization bedgraph
    python $ACCESS_PATH_3/normalize_bedgraph.py --to-mean-signal 1.0 $ACCESS_PATH_1/tmp.unnormalize.$REPLICATE.bgd > $ACCESS_PATH_1/tmp.normalize.$REPLICATE.bdg
    #convert bedgraph to bigwig
    $ACCESS_PATH_3/./bedGraphToBigWig $ACCESS_PATH_1/tmp.normalize.$REPLICATE.bdg $ACCESS_PATH_2/hg38.chrom.sizes $ACCESS_PATH_1/tmp.normalize.$REPLICATE.bw
    NORMALIZE_FILE=$ACCESS_PATH_1/tmp.normalize.$REPLICATE.bw

    get_coverages_for_file $INPUT_FILE $NORMALIZE_FILE
    rm $ACCESS_PATH_1/tmp.unnormalize.$REPLICATE.bgd $ACCESS_PATH_1/tmp.normalize.$REPLICATE.bdg $ACCESS_PATH_1/tmp.normalize.$REPLICATE.bw
    done
done

