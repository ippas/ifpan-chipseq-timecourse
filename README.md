# ifpan-chipseq-timecourse
# ifpan-chipseq-timecourse
# header
## subheader
-------

Dane do RNA-seq dla dexamethazonu (12 plików) pobrano z [https://www.ncbi.nlm.nih.gov/gds/?term=tim+reddy+dexamethasone+rna-seq](https://www.ncbi.nlm.nih.gov/gds/?term=tim+reddy+dexamethasone+rna-seq)
na podstawie plików przygotowano plik raw_macierz.txt i sample.info.txt
Z esembla ściągnięto plik zwierający: 
-Gene stable ID
-Gene stable ID version
-Gene name
plik z nazwy mart.export.txt zmieniono na ID_ID.version_gene.txtx
Z esembla ściągnięto plik zawierający:
-Gene.stable.ID
-Chromosome.scaffold.name
-Gene.start..bp.
-Gene.end..bp.
-Gene.name
-Strand
Zmieniono nazwę pliku z mart.export.txt na gene_chromosome_start_end_strand.txt
Uruchomić skrypt z R: skript_R_clean.R (od 1-123 lini) skrypt wczytuje  pliki raw_macierz.txt (zapisuje do raw.data),  sample.info.txt(zapisuje do samples) i ID_ID.version_gene.txt (zapisuje do ID_ID.version_gene). Wykonuje anove na raw.data, i przy FDR_THRESHOLD=0.001, zostaje wybranych 737 genów (dla dwóch nie została przypisana nazwa, została odrzucone i zostało 735).  Skrypt tworzy heatmap dla RNA-seq dla wybranych transkryptów (z dwoma klastrami), oraz wykres liniowy pokazujący jak zmienia się zawartość transkryptów dla obu klastrów w czasie.


Chip-seq
Dane dla Chip-seq ściągnięto z:
[https://www.ncbi.nlm.nih.gov/gds/?term=tim+reddy+dexamethasone+chip-seq+nr3c1](https://www.ncbi.nlm.nih.gov/gds/?term=tim+reddy+dexamethasone+chip-seq+nr3c1)
dane dla poszczególnych plików znajdują się w pliku: chipseq-file-info.txt (przygotowano go przy pomocy pliku gds.results.txt ściągniętego z powyższej strony i komendy: cat gds_result.txt |  sed -e 's/ /\t/g' | grep -P -B 1 -A 6 "ChIP-seq\ton" | grep -oP 'GSE[0-9]*|[0-9\.]*.hours|supplied\).*ChIP-seq' | xargs -n5 -d'\n' | sed 's/(GR)\t//' | awk '{print $2 "\t" $4*60 "\t""ftp://ftp.ncbi.nlm.nih.gov/geo/series/"$6"nnn/"$7"/suppl/"$7"_RAW.tar""\t"$7}' ) i pobrano pliki przy użyciu komendy: UZUPEŁNIĆ KOMENDĘ

Uruchomić fragment skryptu  skript_R_clean.R (od 126-139), skrypt wczytuje plik
