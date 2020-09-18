#!/bin/bash

function get_coverage {
    CHROMOSOME=$1
    START=$2
    END=$3
    NORMALIZED_FILE=$4
    ./bigWigToBedGraph -chrom=$CHROMOSOME -start=$START -end=$END $ACCESS_PATH_1/tmp.normalize.$REPLICATE.bw stdout | awk -v OFS="\t" '{for (i = $2; i < $3; i++) print $1"\t"i"\t"(i+1)"\t"$4"\t"}' | cut -f4 | awk -v RS= -v OFS='\t' '{$1 = $1} 1'
}

function get_coverages_for_file {
    GENES_INFO_FILE_NAME=$1
    cat $GENES_INFO_FILE_NAME | tail -n +2 | while read GENE_INFO
    do
       CHROMOSOME="chr"$(echo $GENE_INFO | awk '{print $3}')
       START=$(echo $GENE_INFO | awk '{print $4}')
       END=$(echo $GENE_INFO | awk '{print $5}')
       SYMBOL=$(echo $GENE_INFO | awk '{print $2}')
       PREFIX=$SYMBOL"\t"$CHROMOSOME"\t"$START"\t"$END"\t"$TF_NAME"\t"$TIME"\t"$FILE_NAME"\t"
       get_coverage $CHROMOSOME $START $END | awk -v prefix="$PREFIX" '{print prefix $0}'
    done
}
INPUT_FILE=$1
ACCESS_PATH_1=~/ChIP-seq
ACCESS_PATH_2=~/ifpan-chipseq-timecourse/DATA
ACCESS_PATH_3=~/ifpan-chipseq-timecourse/SKRIPTS

cat $ACCESS_PATH_2/chipseq-file-info.txt | grep -v Control | while read LINE
do
  TF_NAME=$(echo $LINE | cut -d $' ' -f 1)
  TIME=$(echo $LINE | cut -d $' ' -f 2)
  GSE=$(echo $LINE | cut -d $' ' -f 4)
  TAR=$ACCESS_PATH_1/$GSE"_RAW.tar"
  REPLICATES=$(tar -tf $TAR | cut -d '_' -f 1 | uniq)
  for REPLICATE in $REPLICATES
  do
    BIGWIG_FILE=$(ls $ACCESS_PATH_1/$REPLICATE*Wig)
    FILE_NAME=$(echo $BIGWIG_FILE | awk -F "/" '{print $5}')
    #bigwig file normalization
    #convert bigwig to bedgraph
    $ACCESS_PATH_3/./bigWigToBedGraph $BIGWIG_FILE /dev/fd/1 > $ACCESS_PATH_1/tmp.unnormalize.$REPLICATE.bgd
    #normalization bedgraph
    python $ACCESS_PATH_3/normalize_bedgraph.py --to-mean-signal 1.0 $ACCESS_PATH_1/tmp.unnormalize.$REPLICATE.bgd > $ACCESS_PATH_1/tmp.normalize.$REPLICATE.bdg
    #convert bedgraph to bigwig
    $ACCESS_PATH_3/./bedGraphToBigWig $ACCESS_PATH_1/tmp.normalize.$REPLICATE.bdg $ACCESS_PATH_2/hg38.chrom.sizes $ACCESS_PATH_1/tmp.normalize.$REPLICATE.bw

    get_coverages_for_file $INPUT_FILE

    rm $ACCESS_PATH_1/tmp.unnormalize.$REPLICATE.bgd $ACCESS_PATH_1/tmp.normalize.$REPLICATE.bdg $ACCESS_PATH_1/tmp.normalize.$REPLICATE.bw
    done
done

