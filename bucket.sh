#!/bin/bash
NUMBER_BUCKET=500
#dla signification_gene FILE_INPUT=big.table.normalize.txt 
#dla signification_gene FILE_OUTPUT=bigtablebucket_normalize.tsv
#dla randowmowych genów FILE_INPUT=random_gene_normalize.txt
#dla randomowych genów FILE_OUTPUT=random_gene_normalize_bucket.tsv
FILE_INPUT=random_gene_normalize.txt
FILE_OUTPUT=random_gene_normalize_bucket.tsv
echo $FILE_INPUT
echo $FILE_OUTPUT
#$FILE_OUTPUT
echo "Tworzenie piku: bigtablebucket_normalize.tsv"


exit 0 
#wycięcie pierwszych kolumn i utworzenie i zapisanie do pliku
cat $FILE_INPUT | cut -f1-3 > $FILE_OUTPUT
for ((j = 4; $j <= 19504; j=(($j+$NUMBER_BUCKET)))); do
   echo j = $j
   #wybieranie 500 komórek
   z=$(($j+499))
   echo z = $z
   #wybranie 501 komórek dla ostatniego obrotu pętli
   t=$(($j+500))
   if [ $j -lt 19504 ]; then
      #sumowanie 500 komorek z każdego wiersza i obliczanie średniej
      awk '{for(i='$j';i<='$z';i++){NUM=NUM?NUM+$i:$i};$(500+1)=NUM;NUM=""} 1' $FILE_INTPUT | awk '{print $501}' | awk '{print $1/500}' > mean1_head.txt
      paste $FILE_OUTPUT mean1_head.txt > tmp_file.txt
      #zmiana nazwy pliku
      mv tmp_file.txt $FILE_OUTPUT
   else
      echo t = $t 
      #sumowanie 501 komórek z każdego wiersza i obliczanie średniej dla ostatniego obrotu
      awk '{for(i='$j';i<='$t';i++){NUM=NUM?NUM+$i:$i};$(501+1)=NUM;NUM=""} 1' $FILE_INTPUT | awk '{print $502}' | awk '{print $1/501}' > mean1_head.txt
      paste $FILE_OUTPUT mean1_head.txt > tmp_file.txt
      #zmina nazwy pliku
      mv tmp_file.txt $FILE_OUTPUT
fi
done
rm mean1_head.txt
