# ifpan-chipseq-timecourse
###### Mateusz Zieba
---

### RNA-seq

RNA-seq data from the ENCODE project were downloaded with [this link](https://www.ncbi.nlm.nih.gov/gds/?term=tim+reddy+dexamethasone+rna-seq) saved as `info-RNA-seq-to-download.txt` in `DATA`. See [this file](DATA/downloads.MD) for details on how the data was downloaded.



Uruchomić skrypt z R: RNA-seq.R (od 1-141 lini) skrypt wczytuje  pliki raw_expression_matrix_dexamethasone.tsv,  sample.info.tsv i ID_ID.version_gene.tsv. Wykonuje jednoczynnikową ANOVE dla każdego transkryptu, i przy FDR_THRESHOLD=0.001, zostaje wybranych 640 genów (dla dwóch nie została przypisana nazwa, została odrzucone i zostało 638).  Skrypt tworzy heatmap dla RNA-seq dla wybranych transkryptów (z dwoma klastrami), 

![Heatmap pokazująca zmieany transkryptów w czasie](PLOTS/heatmap_significant_genes.jpeg)

oraz wykres liniowy pokazujący jak zmienia się zawartość transkryptów dla obu klastrów w czasie.

![Kiku](PLOTS/lineplot_up_down_regulation_significant_genes.jpeg)


Z esembla ściągnięto plik zawierający:
- Gene.stable.ID
- Chromosome.scaffold.name 
- Gene.start..bp.
- Gene.end..bp.
- Gene.name
- Strand

Uruchomienie RNA-seq.R (od 162-212) wczyta plik ~/ifpan-chipseq-timecourse/DATA/gene_chromosome_start_end_strand.tsv i wygeneruje histogram ze środnią ilością transkryptu dal genów down-regulated i up-regulated. RPKM został obliczony przy pomocy wzoru podanego na stronie: [https://www.biostars.org/p/273537/](https://www.biostars.org/p/273537/)

![Kuku](PLOTS/logmean_transcriptlength_signification_gene.jpeg)

![Kiku](PLOTS/histogram_significant_gene_logmean_transcriptlength.jpeg)

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

Z esembla ściągnięto plik zawierający 
- Gene.stable.ID
- Gene.stable.ID.version
- Transcript.length
- gene.name

Nazwę pliku zmieniono na transcript_length.tsv i zapisano w ~/ifpan-chipseq-timecourse/DATA/
Uruchomiono fragment skryptu RNA-seq.R (od 235-260), skrypt wczytuje plik gene_chromosome_start_end_strand.txt, i tworzy tabelę z significant i random genes tablę, którą zapisuje do pliku ~/ifpan-chipseq-timecourse/DATA/significant_and_random_genes_ensemblid_genename_chromosome_start_end.tsv

Przy pomocy pliku significant_and_random_genes_ensemblid_genename_chromosome_start_end.tsv zostają wyciągnięte zakresy lokalizacji dla genów, wartości zostają znormalizowane.
 
```bash
~/ifpan-chipseq-timecourse/SCRIPTS/./bigwig_genomic_range_extract_normalize_to_tsv_bucket.sh ~/ifpan-chipseq-timecourse/DATA/significant_and_random_genes_ensemblid_genename_chromosome_start_end.tsv > ~/ChIP-seq/DATA/significant_random_genes_chip-seq_normalized_bucket_gene_chromosome_start_end_TF_time_file.tsv
```

Skrypty do normalizacji ściągnięto z: 
[https://github.com/porchard/normalize_bedgraph](https://github.com/porchard/normalize_bedgraph)
Skrypty: bigWigToBedGraph, bedGraphToBigWig, bigWigSummary, zciągnięto ze strony: [http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/)

Uruchomienie skryptu RNA-seq.R (od 262-307) wczyta plik /ChIP-seq/DATA/significant_random_genes_chip-seq_normalized_bucket_gene_chromosome_start_end_TF_time_file.tsv i wygeneruje wykresy:

![Kiku](PLOTS/lineplot_significant_random_genes_normalized_bucket.jpeg)

![Kiku](PLOTS/lineplot_significant_random_genes_normalized_bucket_relative_changes.jpeg)


### Peak dla enhancerów
Stworzono plik ze wpólnymi pekami z trzech próbek dla NR3C1 z time60 przy pomocy skryptu creation_file_significant_random_gene_peaks.sh , który tworzy plik tymczasowy ze wspólnymi pekami dl trzech próbek, tworzy tymczasowy plik bed z sinificant i random genes, a na koniec wyciąga peaki dla significant i random genes i zapisuje do pliku: ~/ifpan-chipseq-timecourse/DATA/significant_and_random_genes_ensemblid_genename_chromosome_start-peak_end-peak_regulation.tsv. Pliki niepotrzene usuwa.

```bash
~/ifpan-chipseq-timecourse/SCRIPTS/./creation_file_significant_random_gene_peaks.sh 
```
Następnie przy pomocy polecenia wyciągnięto amplitudy dla peaków i zapisano do pliku:

```bash
~/ifpan-chipseq-timecourse/SCRIPTS/./bigwig_genomic_amplitude_extract_normalize_to_tsv.sh ~/ifpan-chipseq-timecourse/DATA/significant_random_genes_ensemblid_genename_chromosome_start-peak_end-peak_regulation.tsv > ~/ChIP-seq/DATA/significant_random_genes_chip-seq_normalized_gene_chromosome_start-peak_end-peak_TF_time_file_amplitude.tsv
```

Przy pomocy RNA-seq.R (od 309-375) i pliku  ~/ChIP-seq/DATA/significant_random_genes_chip-seq_normalized_gene_chromosome_start-peak_end-peak_TF_time_file_amplitude.tsv przygotowano wykresy:

![Kiku](PLOTS/boxplot_significant_random_genes_strongest_peak.jpeg)

![Kiku](PLOTS/boxplot_significant_random_genes_mean_peaks.jpeg)

![Kiku](PLOTS/barplot_significant_random_genes_strongest_peak.jpeg)

