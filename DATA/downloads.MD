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
-Gene name) 
plik z nazwy mart.export.txt zmieniono na ID_ID.version_gene.tsv i zapisano ~/ifpan-chipseq-timecourse/DATA/

Z esembla ściągnięto plik zawierający:
- Gene.stable.ID
- Chromosome.scaffold.name
- Gene.start..bp.
- Gene.end..bp.
- Gene.name
- Strand

Zmieniono nazwę pliku z mart.export.txt na gene_chromosome_start_end_strand.tsv i zapisano w ~/ifpan-chipseq-timecourse/DATA/
