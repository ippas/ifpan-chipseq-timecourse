#nano coverage.txt
#NUMBER_OF_FILES=$(ls costam | wc -l)
#NUMBER_OF_GENES=$(cat important_gene.txt | tail -n +2 | wc -l )
#nazwa pliku do którego zapisać wynik dla signification_gene: big.table.normalize.txt 
#nazwa pliku do którego zapisać wynik dla randome_gene: random_gene_normalize.txt
function get_coverage {
    CHROMOSOME=$1
    START=$2
    END=$3
    #echo "end: " $END
    FILE=$4
    #echo "file: "$FILE
    ./bigWigToBedGraph -chrom=$CHROMOSOME -start=$START -end=$END normalize.bw stdout | awk -v OFS="\t" '{for (i = $2; i < $3; i++) print $1"\t"i"\t"(i+1)"\t"$4"\t"}' | cut -f4 | awk -v RS= -v OFS='\t' '{$1 = $1} 1'
}

function get_coverages_for_file {
    FILE=$1
    #for signification_gene file is: signification_gene.txt     important_gene.txt
    #for random_gene file is: random_gene.txt
    GENES_INFO_FILE_NAME=random_gene.txt
#important_gene.txt
    cat $GENES_INFO_FILE_NAME | tail -n +2 | while read GENE_INFO
    do
       START=$(echo $GENE_INFO | awk '{print $3}')
       END=$(echo $GENE_INFO | awk '{print $4}')
       CHROMOSOME="chr"$(echo $GENE_INFO | awk '{print $2}')
       SYMBOL=$(echo $GENE_INFO | awk '{print $5}')
       #echo $END $FILE
       PREFIX=$TF_NAME"\t"$TIME"\t"$SYMBOL"\t"
       get_coverage $CHROMOSOME $START $END $FILE $FILE2 | awk -v prefix="$PREFIX" '{print prefix $0}'
    done
}


cat chipseq-file-info.txt | grep -v Control | while read LINE
do
  TF_NAME=$(echo $LINE | cut -d $' ' -f 1)
  TIME=$(echo $LINE | cut -d $' ' -f 2)
  GSE=$(echo $LINE | cut -d $' ' -f 4)
  TAR=$GSE"_RAW.tar"
  REPLICATES=$(tar -tf $TAR | cut -d '_' -f 1 | uniq)
  #echo $TF_NAME $TIME $GSE
  for REPLICATE in $REPLICATES
  do
    #echo $(ls $REPLICATE*Wig)
    BIGWIG_FILE=$(ls $REPLICATE*Wig)
    #normalizacja big wig
    #echo "zamiana biWig to bedGraph"
    ./bigWigToBedGraph $BIGWIG_FILE /dev/fd/1 > tmp.unnormalize.bdg
    #echo "normalizacja bedGraph"
    python normalize_bedgraph.py --to-mean-signal 1.0 tmp.unnormalize.bdg > tmp.normalize.bdg
    #echo "zamiana bedGraph to bigWig"
    #echo "zamiana bedGraph to bigWig"
    ./bedGraphToBigWig tmp.normalize.bdg hg38.chrom.sizes normalize.bw
    #echo "save to object"
    #BIGWIG_FILE="normalize.bw"
    #echo "file_normalize: " $BIGWIG_FILE_normalise
    #jeśli źle do tego miejsca usunąć
    get_coverages_for_file $BIGWIG_FILE
    rm tmp.unnormalize.bdg
    rm tmp.normalize.bdg
    rm normalize.bw
  done
done


exit 0

for ((j=1; $j <= 35; j++)); do
for ((i=2; $i <= 736; i++)); do
echo "i" $i
echo "j" $j
start=$(awk '{print $3}' important_gene.txt |sed "${i}q;d") 
end=$(awk '{print $4}' important_gene.txt | sed "${i}q;d")
chromosome=$(awk '{print "chr"$2}' important_gene.txt | sed "${i}q;d")
file=$(grep 'NR3C1' chipseq-file-info.txt | sort -k2 -n | awk '{print $4"_RAW.tar"}' | xargs -i bash -c "tar -tf {} | grep bigWig" | sed "${j}q;d")
awk '{print $5}' important_gene.txt | sed "${i}q;d" | awk '{print $1"_"'$j'}'> gene.txt 
echo $start
echo $file
#chromosome=$(cut -f10 important_gene.txt | sed $i,$i\p)
./bigWigToBedGraph -chrom=$chromosome -start=$start -end=$end $file stdout | 
awk -v OFS="\t" '{for (i = $2; i < $3; i++) print $1"\t"i"\t"(i+1)"\t"$4"\t"}' | cut -f4 | awk -v RS= -v OFS='\t' '{$1 = $1} 1' > line1.txt
#paste coverage.txt line.txt
paste gene.txt line1.txt >> coverage_NR3C1_numbersample.txt
rm gene.txt
rm line1.txt
done
done
#GSM2421792_ENCFF778FPM_signal_of_unique_reads_GRCh38.bigWig stdout
