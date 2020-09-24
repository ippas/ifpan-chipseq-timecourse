# ifpan-chipseq-timecourse
###### Mateusz Zięba
---

### RNA-seq

Dane do RNA-seq dla dexamethazonu (12 plików) pobrano z [https://www.ncbi.nlm.nih.gov/gds/?term=tim+reddy+dexamethasone+rna-seq](https://www.ncbi.nlm.nih.gov/gds/?term=tim+reddy+dexamethasone+rna-seq) 

Pobrano plik tekstowy powyższej strony i zapisano jako info-RNA-seq-to-download.txt i zapisano w ~/ifpan-chipseq-timecourse/DATA

Przy pomocy komendy przygotowano plik (mRNA_seq-file-info.tsv) z informacjami potrzebnymi do pobrania plików.

```bash
cat ~/ifpan-chipseq-timecourse/DATA/info-RNA-seq-to-download.txt | 
   tail +2 | 
   sed 'N;N;N;N;N;N;N;s/\n/ /g' | 
   grep 'mRNA-seq'  | 
   awk '{print $19"\t"$21*60"\t"$56"suppl/"$59"_RAW.tar\t"$59}' | 
   sort -n -k2 > ~/ifpan-chipseq-timecourse/DATA/mRNA_seq-file-info.tsv
```

Przy pomocy polecenia pobrano pliki RNA-seq dla deksametazonu

```bash
cat ~/ifpan-chipseq-timecourse/DATA/mRNA_seq-file-info.tsv | 
   cut -f3 | 
   xargs -i bash -c 'wget {} -P ~/dexamethasone/DOWNLOAD/'
```

Rozpakowano do folderu EXTRACT, komendą:

```bash
ls ~/dexamethasone/DOWNLOAD/*tar | 
   xargs -i bash -c 'tar -C ~/dexamethasone/EXTRACT -xvf {}'
```

Tworzenie pliku sample.info.txt, który zawiera: "samplied", "time", "replicate" i "file". Wykonano przy użyciu komendy:


```bash
ls ~/dexamethasone/EXTRACT/*.tsv.gz | 
   xargs -i bash -c ' zgrep -H "rep" {} | 
   head -1' | 
   cut -d "/" -f6,36 | 
   sed  's/:/\t/; s/\//\t/' | 
   awk '{print $6"\t"$1}' | 
   sed 's/_/\t/' | 
   cut -c2- | 
   awk '{print "t"$1"_"$2"\t"$1"\t"$2"\t"$3}' | 
   sed 's/05/0\.5/g2' | 
   awk '{print $1"\t"$2*60"\t"$3"\t"$4}' | 
   sort -k2 -n |
   sed '1 i\samplied\ttime\treplicate\tfile ' > ~/ifpan-chipseq-timecourse/DATA/sample.info.tsv
```

Przypomocy pliku sample.info.tsv i plików *gene_quantifications_GRCh38.tsv.gz przygotowano plik raw_expression_matrix_dexamethasone.tsv dla RNAseq. Wykonano przy użyciu komendy:

```bash
awk 'FNR==NR { a[FNR""] = $0; next } { print a[FNR""]"\t" $0 }' <(cat ~/ifpan-chipseq-timecourse/DATA//sample.info.tsv | 
   head -2 | 
   tail +2 | 
   cut -f4 | 
   xargs -i bash -c 'zcat ~/dexamethasone/EXTRACT/{}' | 
      tail +3 | 
      cut -f1,6 )  <(cat ~/ifpan-chipseq-timecourse/DATA/sample.info.tsv | 
   cut -f4 | 
   tail +2 | 
   xargs -i bash -c 'zcat ~/dexamethasone/EXTRACT/{} | 
      tail +3 | 
      cut -f7' | 
   awk -v row=$(cat ~/ifpan-chipseq-timecourse/DATA/sample.info.tsv | 
      head -2 | 
      tail -1 | 
      cut -f4 | 
      xargs -i bash -c 'zcat ~/dexamethasone/EXTRACT/{}' | 
         wc -l | xargs -i bash -c 'echo $(({}-2))') '{A[(NR-1)%row]=A[(NR-1)%row]$0"\t ";next}END{for(i in A)print A[i]}') | 
   awk -i inplace -v first=$(cat ~/ifpan-chipseq-timecourse/DATA/sample.info.tsv | 
      cut -f1 | 
      tail +2 | 
      sed '1 i\Geneid\nLength' | tr "\n" ":" ) 'BEGINFILE{print first}{print}' |  
   sed 's/:/\t/g' | 
   sed s'/ $//' | 
   sed s'/\t$//' > ~/ifpan-chipseq-timecourse/DATA/raw_expression_matrix_dexamethasone.tsv
```

Z esembla ściągnięto plik zwierający: 
-Gene stable ID
-Gene stable ID version
-Gene name
plik z nazwy mart.export.txt zmieniono na ID_ID.version_gene.tsv i zapisano ~/ifpan-chipseq-timecourse/DATA/

Z esembla ściągnięto plik zawierający:
- Gene.stable.ID
- Chromosome.scaffold.name
- Gene.start..bp.
- Gene.end..bp.
- Gene.name
- Strand

Zmieniono nazwę pliku z mart.export.txt na gene_chromosome_start_end_strand.tsv i zapisano w ~/ifpan-chipseq-timecourse/DATA/

Uruchomić skrypt z R: skript_R_clean.R (od 1-141 lini) skrypt wczytuje  pliki raw_macierz.txt (zapisuje do raw.data),  sample.info.tsv i ID_ID.version_gene.tsv. Wykonuje anove na raw.data, i przy FDR_THRESHOLD=0.001, zostaje wybranych 640 genów (dla dwóch nie została przypisana nazwa, została odrzucone i zostało 638).  Skrypt tworzy heatmap dla RNA-seq dla wybranych transkryptów (z dwoma klastrami), 

![Heatmap pokazująca zmieany transkryptów w czasie](PLOTS/heatmap_significant_genes.jpeg)

oraz wykres liniowy pokazujący jak zmienia się zawartość transkryptów dla obu klastrów w czasie.

![Kiku](PLOTS/lineplot_up_down_regulation_significant_genes.jpeg)

### Chip-seq

Dane dla Chip-seq ściągnięto z:
[https://www.ncbi.nlm.nih.gov/gds/?term=tim+reddy+dexamethasone+chip-seq+nr3c1](https://www.ncbi.nlm.nih.gov/gds/?term=tim+reddy+dexamethasone+chip-seq+nr3c1),
Pobrano plik tekstowy gds_result.txt którego nazwę zmieniono na info-chip-seq-to-download.txt i zapisano w ~/ifpan-chipseq-timecourse/DATA/, przypomocy tego pliku przygotowano plik zawierający nazwę czynnika transkrypcyjnego, czas w którym nastąpiło utrwalenie komórek, adress ftp do strony z której trzeba pobrać plik z danymi. Plik wykonano przy użyciu komendy:

```bash
cat ~/ifpan-chipseq-timecourse/DATA/info-chip-seq-to-download.txt | 
   sed -e 's/ /\t/g' | 
   grep -P -B 1 -A 6 "ChIP-seq\ton" | 
   grep -oP 'GSE[0-9]*|[0-9\.]*.hours | supplied\).*ChIP-seq' | 
   xargs -n5 -d'\n' | 
   sed 's/(GR)\t//' | 
   awk '{print $2 "\t" $4*60 "\t""ftp://ftp.ncbi.nlm.nih.gov/geo/series/"$6"nnn/"$7"/suppl/"$7"_RAW.tar""\t"$7}' > ~/ifpan-chipseq-timecourse/DATA/chipseq-file-info.tsv
```

Przy pomocy pliku chipseq-file-info.tsv pobrano pliki przy użyciu komendy:

```bash
cat ~/ifpan-chipseq-timecourse/DATA/chipseq-file-info.tsv | 
   cut -f3 | 
   xargs -i bash -c 'wget {} -P ~/ChIP-seq/DOWNLOAD'
```

Z esembla ściągnąć plik zawierający 
- Gene.stable.ID
- Gene.stable.ID.version
- Transcript.length
- gene.name

Nazwę pliku zmienić na transcript_length.tsv
Uruchomić fragment skryptu  skript_R_clean.R (od 126-139), skrypt wczytuje plik gene_chromosome_start_end_strand.txt, i tworzy plik signification_gene.txt
W skrypcie bigwig_genomic_range_extract_normalize_totsv.sh, do GENES_INFO_FILE_NAME przypisać plik signification_gene.txt, wynik działania skryptu zapisać big.table.normalize.tsv. Następnie w skrypcie bucket.sh, do do FILE_INPUT przypisać big.table.normalize.txt, a do FILE_OUTPUT bigtablebucket_normalize.tsv. 

```bash
~/ifpan-chipseq-timecourse/SKRIPTS/./bigwig_genomic_range_extract_normalize_to_tsv_bucket.sh ~/ifpan-chipseq-timecourse/DATA/significant_genes_ensemblid_genename_chromosome_start_end.tsv > ~/ChIP-seq/DATA/significant_genes_chip-seq_gene_chromosome_start_end_TF_time_file.tsv

~/ifpan-chipseq-timecourse/SCRIPTS/./bigwig_genomic_range_extract_normalize_to_tsv_bucket.sh ~/ifpan-chipseq-timecourse/DATA/significant_and_random_genes_ensemblid_genename_chromosome_start_end.tsv > ~/ChIP-seq/DATA/significant_and_random_genes_chip-seq_gene_chromosome_start_end_TF_time_file.tsv


```

Skrypty do normalizacji ściągnięto z: 
[https://github.com/porchard/normalize_bedgraph](https://github.com/porchard/normalize_bedgraph)
oraz skrypty: bigWigToBedGraph, bedGraphToBigWig.

Uruchomić skrypt z R (linie 141-173) w którym zostaje obliczone RPKM (średnia ilość transkryptu) dla signification_gene.txt zgodnie ze wzorem znajdującym sie na stronie: [https://www.biostars.org/p/273537/](https://www.biostars.org/p/273537/), oraz wygenerowany zostaje wykres w skali logarytmiczej,  ukazujący średnią długość transkryptów.

![Heatmap, pokazuje jak zmienia się ilość transkryptu w czasie](PLOTS/logmean_transcriptlength_signification_gene.jpeg)

Następnie urochomienie kodu od 176-196 skript_R_clean.R w którym zostaje obliczone RPKM dla wszystkich transkryptów, oraz zostaje wylosowana próba losowa 1000 genów (których RPKM znajduję sie pomiędzy 2 a 8192 - ponieważ w takiej wiekości mieściły się transkrypty dla signification_gene, odczytano z wykresu)

Uruchomienie skript_R_clean.R  (linie 199 -209) przygotuje stworzy plik random_gene.txt zawierający:
- numer chromosomu
- miejsce 10000 przed miejscem startu transkrypcji
- miejsce 10000 po miejscu startu transkrypcji
- nazwę genu

W skrypcie bigwig_genomic_range_extract_normalize_totsv.sh do GENES_INFO_FILE_NAME przypisać plik random_gene.txt, wynik działania skryptu zapisać do random_gene_normalize.txt.  
Następnie w skrypcie bucket.sh, do do FILE_INPUT przypisać random_gene_normalize.txt, a do FILE_OUTPUT random_gene_normalize_bucket.tsv. 
Uruchomić skript_R_clean.R (linie 211-253) skrypt wczytuje random_gene_normalize_bucket.tsv, tworzy tabelę dla znaczących i losowych genów i toworzy wykres.

![Kiku](PLOTS/signification_random_gene_plot.jpeg)

![Kiku](PLOTS/signification_random_gene_normalize.jpeg)

![Kiku](PLOTS/signification_random_gene_relative_change.jpeg)


### Peak dla enhancerów
przygotowanie do peaków
```bash
cat ~/ifpan-chipseq-timecourse/DATA/significant_and_random_genes_ensemblid_genename_chromosome_start_end.tsv | tail +2 | awk 'BEGIN {OFS ="\t"}{print "chr"$3,$4,$5, $0}' | grep -v -F 'chrMT'
cat ~/ifpan-chipseq-timecourse/DATA/chipseq-file-info.tsv | grep 'NR3C1.60.f' | awk '{print $4"_RAW.tar"}' | xargs -i bash -c 'tar -tf {}' | grep 'bed'
 bedtools intersect -a tmp_significant_random.bed -b $three_file -wb -wa | sort
 bedtools intersect -a tmp_significant_random.bed -b $three_file -wo | grep '3.chr' 

three_file=$(cat ~/ifpan-chipseq-timecourse/DATA/chipseq-file-info.tsv | grep 'NR3C1.60.f' | awk '{print $4"_RAW.tar"}' | xargs -i bash -c 'tar -tf {}' | grep 'bed' | tr "\n" " ")


bedtools intersect -a tmp_significant_random.bed -b $three_file -wo  | tr "\n" ":" | grep -iPo 'chr[0-9\t]*ENSG[0-9]*\t[0-9a-z.]*[0-9\t]*[a-z-]*\t1\tchr[0-9a-z\t_.]*:chr[0-9\t]*ENSG[0-9]*\t[0-9a-z.]*[0-9\t]*[a-z-]*\t2\tchr[0-9a-z\t_.]*:chr[0-9\t]*ENSG[0-9]*\t[0-9a-z.]*[0-9\t]*[a-z-]*\t3\tchr[0-9a-z\t_.]*:' | sed 's/:/\n/g' | sed -r '/^\s*$/d'     

bedtools intersect -a tmp_significant_random.bed -b $three_file -wo | tr "\n" ":" |  grep -ioP 'chr[0-9x-z\t]*ENSG[0-9]*\t[0-9a-z.-]*\t[0-9x][0-9 \t]*\t[a-z0-]*\t1\tchr[0-9a-z_.\t]*:chr[0-9x-z\t]*ENSG[0-9]*\t[0-9a-z.-]*\t[0-9x][0-9 \t]*\t[a-z0-]*\t2\tchr[0-9a-z_.\t]*:chr[0-9x-z\t]*ENSG[0-9]*\t[0-9a-z.-]*\t[0-9x][0-9 \t]*\t[a-z0-]*\t3\tchr[0-9a-z_.\t]*:'  | sed 's/:/\n/g' | sed -r '/^\s*$/d' 

```



Wycięto peaki dla enhancerów dla genów upregulowanych, dla NR3C1 z pierwszej godziny przy pomocy polecenia:
```console
bedtools intersect -a GSM2421929_ENCFF835HHK_peaks_GRCh38.bed.gz -b up_regulated_gene.BED -wb | cut -f1,2,3,14 | awk '{print $1"\t"$2-1000"\t"$3+1000"\t"$4}'  > upregulatedgne_peaks_NR3C1_time60_file1.tsv
```
dla trzech plików:
```console
grep 'NR3C1.60.f' chipseq-file-info.txt | sort -k2 -n | awk '{print $4"_RAW.tar"}' | xargs -i bash -c "tar -tf {} | grep bed.gz" |  xargs -i bash -c 'bedtools intersect -a {} -b up_regulated_gene.BED'
```
Następnie uruchomić skrypt:
```console
./peak_NR3C1_normalize.sh > peaks_NR3C1_normalize_file1_time60.tsv
```
Wybieranie randomowych peaków (154) przy pomocy:
```console
zcat GSM2421929_ENCFF835HHK_peaks_GRCh38.bed.gz | shuf -n 154 | cut -f1,2,3 | awk '{print $1"\t"$2-1000"\t"$3+1000"\tRANDOM"}' > random_peaks_NR3C1_time60_file1.tsv

```

Uruchomić skrypt wyciągający peaki dla randomowych:
```console
./peak_NR3C1_normalize.sh > random_peaks_NR3C1_time60_file1_normalize.tsv 
```
Pliki: 
- peaks_NR3C1_normalize_file1_time60.tsv
- random_peaks_NR3C1_time60_file1_normalize.tsv  
wczytano do R i zrobiono wykresy: 

![Kiku](PLOTS/upregulatedgene_peak_NR3C1_file1_time60.jpeg)

![Kiku](PLOTS/random_peak_NR3C1_file1_time60.jpeg)

zapisanie peaków wspólnych, z trzech plików dla NR3C1 po jednej godzinie 

```console
bedtools intersect -a up_regulated_gene.BED -b GSM2421929_ENCFF835HHK_peaks_GRCh38.bed.gz GSM2421930_ENCFF044MLR_peaks_GRCh38.bed.gz GSM2421931_ENCFF597LEE_peaks_GRCh38.bed.gz -u > upregulated_gene_allfiles_peaks_NR3C1.tsv

```
wyciągnięcie peaków, i przygotowanie do puszczenia w skrypcie
```console
bedtools intersect -a upregulated_gene_allfiles_peaks_NR3C1.tsv -b GSM2421929_ENCFF835HHK_peaks_GRCh38.bed.gz -wb | awk '{print $1"\t"$6-1000"\t"$7+1000"\t"$4}' > upregulatedgene_peak_NR3C1_time60_allfile_intoskript.tsv
```


Uruchomienie skryptu:

```console
./peak_NR3C1_normalize.sh > upregulated_peaks_NR3C1_allfile_time60_allTF.tsv

```
grep 'NR3C1.60.f' chipseq-file-info.txt | awk '{print $4"_RAW.tar"}' | xargs -i bash -c "tar -tf {} | grep bed.gz" | echo $(cat) | xargs -i bash -c 'bedtools intersect -a down_regulated_gene.BED -b {} -wb -wa' | head

tworzenie pliku części wspólnej plików: 
- GSM2421929_ENCFF835HHK_peaks_GRCh38.bed.gz 
- GSM2421930_ENCFF044MLR_peaks_GRCh38.bed.gz 
- GSM2421931_ENCFF597LEE_peaks_GRCh38.bed.gz
oraz wyciągnięcie i zapisanie peaków down i up regulowanych

```console
bedtools intersect -a  GSM2421929_ENCFF835HHK_peaks_GRCh38.bed.gz -b GSM2421930_ENCFF044MLR_peaks_GRCh38.bed.gz > tmp_peak_file1-2.tsv
bedtools intersect -a GSM2421931_ENCFF597LEE_peaks_GRCh38.bed.gz -b tmp_peak_file1-2.tsv > peaks_allfile_NR3C1_time60.bed
bedtools intersect -a up_regulated_gene.BED -b peaks_allfile_NR3C1_time60.bed -wb | awk '{print $1"\t"$6-1000"\t"$7+1000"\t"$4}' > upregulated_peak_NR3C1_time60_allfile.tsv
down_regulated_gene.BED -b peaks_allfile_NR3C1_time60.bed -wb | awk '{print $1"\t"$6-1000"\t"$7+1000"\t"$4}' > downregulated_peak_NR3C1_time60_allfile.tsv

```
puszczenie skryptu
```console
./peak_NR3C1_normalize.sh > upregulated_peaks_NR3C1_allfile_time60_allTF_normalize.tsv
```
sprawdzić i poprawić co nie potrzebne

komenda wyciągająca największą wartość z wiersza, test
```bash
for((i=1;i<=100;i+=1)); do cat upregulated_peaks_NR3C1_allfile_time60_allTF_normalize.tsv | cut -f8-3000 | head -100 | head -n$i | tail -1 | sed 's/\t/\n/g' | sort -n | tail -1; done 
```
komenda łącząca dwa wyniki
```bash
awk 'FNR==NR { a[FNR""] = $0; next } { print a[FNR""]"\t" $0 }' <(cat upregulated_peaks_NR3C1_allfile_time60_allTF_normalize.tsv | cut -f1-7 | head -10) <(for((i=1;i<=10;i+=1)); do cat upregulated_peaks_NR3C1_allfile_time60_allTF_normalize.tsv | cut -f8-3000 | head -100 | head -n$i | tail -1 | sed 's/\t/\n/g' | sort | tail -1; done) 
```
```bash
for((i=1;i<=$(awk 'END {print NR}' important_gene.txt); i+=1));do echo $i;done
cat upregulated_peaks_NR3C1_allfile_time60_allTF_normalize.tsv | head -1 | cut -f1-$(cat upregulated_peaks_NR3C1_allfile_time60_allTF_normalize.tsv | awk '{print NF}' | sort -n | tail -1)
```

```bash
awk 'FNR==NR { a[FNR""] = $0; next } { print a[FNR""]"\t" $0 }' <(cat upregulated_peaks_NR3C1_allfile_time60_allTF_normalize.tsv | cut -f1-7) <(for((i=1;i<=395968;i+=1)); do cat upregulated_peaks_NR3C1_allfile_time60_allTF_normalize.tsv | cut -f8-3560 | head -n$i | tail -1 | sed 's/\t/\n/g' | sort -n | tail -1; done) > upregulated_peaks_NR3C1_allTF_amplitude_time60.tsv
```



Wyciągnięcie 1000 randomowych peaków
```bash
cat peaks_allfile_NR3C1_time60.bed | shuf -n 1000 | cut -f1,2,3 | awk '{print $1"\t"$2-1000"\t"$3+1000"\tRANDOM"}' > random1000_peak.tsv
```
uruchomienie skryptu
```bash
./peak_amplitude_NR3C1_normalize.sh random1000_peak.tsv > random1000_peaks_NR3C1_allTF_amplitude_time60.tsv
```

przygotowanie peaków dla randomowych genów
```bash
bedtools intersect -a random_gene_into_peak.tsv -b peaks_allfile_NR3C1_time60.bed -wb | awk '{print $1"\t"$6-1000"\t"$7+1000"\t"$4}' > random_gen_peak_NR3C1_time60_allfile.tsv
```

uruchomienie skryptu wyciągającego amplitudę dla peaku
```bash
```


![Kiku](PLOTS/barplot_upregulated_amplitude_peak_allTF_NR3C1_time60.jpeg)
![Kiku](PLOTS/barplot_downregulated_amplitude_peak_allTF_NR3C1_time60.jpeg)
![Kiku](PLOTS/barplot_random1000_amplitude_peak_allTF_NR3C1_time60.jpeg)
![Kiku](PLOTS/barplot_random_gen_amplitude_peak_allTF_NR3C1_time60.jpeg)

![Kiku](PLOTS/boxplot_upregulated_amplitude_peak_allTF_NR3C1_time60.jpeg)
![Kiku](PLOTS/boxplot_downregulated_amplitude_peak_allTF_NR3C1_time60.jpeg)
![Kiku](PLOTS/boxplot_random1000_amplitude_peak_allTF_NR3C1_time60.jpeg)
![Kiku](PLOTS/boxplot_random_gen_amplitude_peak_allTF_NR3C1_time60.jpeg)

![Kiku](PLOTS/boxplot_freey_upregulated_amplitude_peak_allTF_NR3C1_time60.jpeg)
![Kiku](PLOTS/boxplot_freey_downregulated_amplitude_peak_allTF_NR3C1_time60.jpeg)
![Kiku](PLOTS/boxplot_freey_random1000_amplitude_peak_allTF_NR3C1_time60.jpeg)
![Kiku](PLOTS/boxplot_freey_random_gen_amplitude_peak_allTF_NR3C1_time60.jpeg)
