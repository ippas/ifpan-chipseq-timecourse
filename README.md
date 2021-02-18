# ifpan-chipseq-timecourse
###### Mateusz Zieba
---

  

### RNA-seq

RNA-seq data from the ENCODE project were downloaded with [this link](https://www.ncbi.nlm.nih.gov/gds/?term=tim+reddy+dexamethasone+rna-seq) saved as `info-RNA-seq-to-download.txt` in `DATA`. See [this file](DATA/downloads.MD) for details on how the data was downloaded.

[analysis_transcriptome.R](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/analysis_transcriptome.R) contains RNAseq analysis code:

1. Downloaded data are loaded (including data downloaded from genebank)
2. Raw read counts are normalised (normalize quantiles)
3. Differentially Expressed Genes (DEGs) are selected based on one-way ANOVA (n = 745 genes with FDR < 0.0000001 treshold)
4. DEGs are divided into two clusters - upregulated (n = 370) and downregulated (n = 375)
5. A list of random genes with similar expression levels is chosen (n = 745)
6. Following graphs are prepared: line graph, heatmap, boxplot of weighted timepoint of max change, RPKM histogram (calculated according to [this info](https://www.biostars.org/p/273537/), transcript length histogram (accordint to: [transcript.length.tsv](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/DATA/transcript_length.tsv)
7. Additional tables for further analyses are prepared: [promotores_peaks_info.tsv]()https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/DATA/promotores_peaks_info.tsv and [enhancer_info.tsv](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/DATA/enhancer_info.tsv) do znalezienia enhancerów - these ables contain information about gene identifiers, location and cluster.

![Kiku](PLOTS/heatmap_expression_significant.svg)

![Kiku](PLOTS/lineplot_change_expression.svg)

![Kiku](PLOTS/histogram_RPKM.svg)

![Kiku](PLOTS/histogram_transcript_lenght.svg)

### ChIP-seq

Dane dla Chip-seq [ściągnięto z](https://www.ncbi.nlm.nih.gov/gds/?term=tim+reddy+dexamethasone+chip-seq+nr3c1),
Pobrano plik tekstowy gds_result.txt którego nazwę zmieniono na [info-chip-seq-to-download.txt](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/DATA/info-chip-seq-to-download.txt). Przypomocy tego pliku stworzono plik zawierający: 
- nazwę czynnika transkrypcyjnego, 
- czas w którym nastąpiło utrwalenie komórek, 
- adress ftp do pobrania danych. 
Plik wykonano przy użyciu komendy:

```bash
cat ~/ifpan-chipseq-timecourse/DATA/info-chip-seq-to-download.txt | 
   sed -e 's/ /\t/g' | 
   grep -P -B 1 -A 6 "ChIP-seq\ton" | 
   grep -oP 'GSE[0-9]*|[0-9\.]*.hours | supplied\).*ChIP-seq' | 
   xargs -n5 -d'\n' | 
   sed 's/(GR)\t//' | 
   awk '{print $2 "\t" $4*60 "\t""ftp://ftp.ncbi.nlm.nih.gov/geo/series/"$6"nnn/"$7"/suppl/"$7"_RAW.tar""\t"$7}' > ~/ifpan-chipseq-timecourse/DATA/chipseq-file-info.tsv
```

Używając [info-chip-seq-to-download.txt](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/DATA/info-chip-seq-to-download.txt) pobrano pliki korzystając z komendy: 

```bash
cat ~/ifpan-chipseq-timecourse/DATA/chipseq-file-info.tsv | 
   cut -f3 | 
   xargs -i bash -c 'wget {} -P ~/ChIP-seq/DOWNLOAD'
```

2. Uruchomić skrypt [extract_data_chipseq1.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/extract_data_chipseq1.sh) komendą:
 
 ```bash
 ~/ifpan-chipseq-timecourse/SCRIPTS/extract_data_chipseq1.sh
 ```
> Skrypt wykonuje:
 - przy pomocy [promotores_peaks_info.tsv](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/DATA/promotores_peaks_info.tsv) i [bigwig_genomic_bucket500_extract_normalize_to_tsv.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigwig_genomic_bucket500_extract_normalize_to_tsv.sh) wyciąga informacje dotyczące przyłączania się TF do DNA w zakresie +/-10000 od TSS i zapisuje do ~/ChIP-seq/DATApromotores_peaks_value.tsv
 - z trzech próbek tworzy część wspólną peaków w genomie dla NR3C1 dla time = 60
 - przy pomocy wyznaczonej części wspólnej peaków, wyciąga enhancery znajdujące się w zakresach znajdujących się w [enhancer_info.tsv](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/DATA/enhancer_info.tsv) (+/-100000 od TSS)
 - dla wybranych enhancerów korzystając z [bigwig_genomic_amplitude_extract_normalize_to_tsv.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigwig_genomic_amplitude_extract_normalize_to_tsv.sh) wyciąga maksymalną wartość przyłączania się TF do peaku i zapisuje dane do: ~/ChIP-seq/DATA/enhancer_amplitude_value.tsv 
 - dla wybranych enhancerów korzystając z [bigwig_genomic_range_extract_normalize_to_tsv.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigwig_genomic_range_extract_normalize_to_tsv.sh) wyciąga wartości przyłączeń się TF na całym zakresie peaku i zapisuje dane do: ~/ChIP-seq/DATA/enhancer_peaks_value.tsv

3. Uruchomić skrypt [visualization_promoters.R](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/visualization_promoters.R), który wykonuje wykresy:
- wczytuje plik ~/ChIP-seq/DATA/promotores_peaks_value.tsv, który jest wynikiem skryptu [extract_data_chipseq1.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/extract_data_chipseq1.sh)
- wykres liniowy predstawiający średnią z grupy przyłączenia się TF do promotorów w pozycjach =/- 10000 od TSS
- wykres liniowy predstawiający średnią z grupy przyłączenia się TF do promotorów w pozycjach =/- 10000 od TSS, w odniesieniu do genów z grupy random
- wykres liniowy predstawiający średnią z grupy przyłączenia się TF (wybrane cztery TF: EP300, H3K27ac, H3K4me1, NR3C1) do promotorów w pozycjach =/- 10000 od TSS,
- wykres liniowy predstawiający średnią z grupy przyłączenia się TF (wybrane cztery TF: EP300, H3K27ac, H3K4me1, NR3C1) do promotorów w pozycjach =/- 10000 od TSS, w odniesieniu do genów z grupy random


![Kiku](PLOTS/lineplot_promotores_fourTF.svg)

![Kiku](PLOTS/lineplot_promotores_relative_changes_fourTF.svg)

4. Uruchomić skrypt [visualization_amplitudes.R](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/visualization_amplitudes.R), który:
- wczytuje plik ~/ChIP-seq/DATA/enhancer_amplitude_value.tsv, będący wynikiem skryptu [extract_data_chipseq1.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/extract_data_chipseq1.sh)
- wykres zmian amplitud przyłączania TF w czasie dla najsilniejszych peaków w genie
- wykres zmian amplitud przyłączania TF (wybrane cztery TF: EP300, H3K27ac, H3K4me1, NR3C1) w czasie dla najsilniejszych peaków w genie
- wykres przedstawiający Mean Weighted Time (MWT) dla TF i gene regulation
- wykres przedstawiający Mean Weighted Time (MWT) dla EP300, H3K27ac, H3K4me1, NR3C1 i genów up-regulated


![Kiku](PLOTS/boxplot_enhancer_amplitude_fourTF.svg)

![Kiku](PLOTS/histogram_distance_enhancer_TSS.svg)

![Kiku](PLOTS/boxplot_MWT_MCTP.svg)

![Kiku](PLOTS/boxlineplot_expression_NR3C1_EP300.svg)


 
Uruchomić [extract_data_chipseq2.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/extract_data_chipseq2.sh) komendą:

```bash
 ~/ifpan-chipseq-timecourse/SCRIPTS/extract_data_chipseq1.sh
 ```
 
 
 Skrypt [[korzysta z]](#1)

Uruchomić [visualization_enhancer_range.R](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/visualization_enhancer_range.R)

![Kiku](PLOTS/lineplot_enhancer_fourTF.svg)

![Kiku](PLOTS/lineplot_enhancer_relative_changes_fourTF.svg)

![Kiku](PLOTS/heatmap_enhancer.png)





### Analiza EP300 


#### Delta EP300
Uruchomić skrypt [delta_ep300.R](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/delta_ep300.R), którego kod wykonuje:
- wczytuje 
- dokonuje wyboru peaków do analizy poprzez obliczenie różnicy pomiędzy minimalną wartością przyłączania się EP300 do DNA, a maksymalną wartością, przypisuje enhancery do najbliższego genu, a gdy dwa peaki zostaną przypisane do tego samego genu to wybiera ten którego start znajduje się bliżej TSS i wybiera top 100 peaków/genów, jednak w dalszej analizie jeden gen zostje usunięty z powodu, że 
- wykonuje heatmape ekspressji genów

![Kiku](PLOTS/heatmap_expression_delta_ep300.svg)

![Kiku](PLOTS/lineplot_change_expression_delta_ep300.svg)

![Kiku](PLOTS/boxplot_amplitude_delta_ep300.svg)

![Kiku](PLOTS/boxplot_MCTP_MWT_delta_ep300.svg)

Uruchomić skrypt [extract_data_chipseq2_ep300.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/extract_data_chipseq2_ep300.sh), który wykonuje:
- przy pomocy [bigwig_genomic_range_extract_normalize_to_tsv_bucket10.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigwig_genomic_range_extract_normalize_to_tsv_bucket10.sh) i pliku [enhancer_bigrange_top_ep300_info.tsv](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/DATA/enhancer_bigrange_top_ep300_info.tsv) wyciąga dane z przyłączania TF w enhancerach, a wynik zapisuje do pliku ~/ChIP-seq/DATA/enhancer_bigrange_top_ep300.tsv
- rzy pomocy [bigwig_genomic_range_extract_normalize_to_tsv_bucket10.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigwig_genomic_range_extract_normalize_to_tsv_bucket10.sh) i pliku [enhancer_bigrange_delta_ep300_info.tsv](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/DATA/enhancer_bigrange_delta_ep300_info.tsv) wyciąga dane z przyłączania TF w enhancerach, a wynik zapisuje do pliku ~/ChIP-seq/DATA/enhancer_bigrange_delta_ep300.tsv


Uruchomić skrypt [analysis_enhancer_range_ep300.R](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/analysis_enhancer_range_ep300.R), który wykonuje:
- wczytuje dane dla genów top ep300 ~/ChIP-seq/DATA/enhancer_bigrange_top_ep300.tsv i dla genów delta ep300 ~/ChIP-seq/DATA/enhancer_bigrange_delta_ep300.tsv, oba pliki zawierają wartości przyłączeń TF +/- 10000 od pozycji amplitudy dla peaków NR3C1
-wykonuje uśrednione wykresy liniowy zmian przyłączania TF w zakresie +/- 2000 od pozycji amplitudy peaków dla NR3C1 dla czterech TF EP300, H3K27ac, H3K4me1, NR3C1, zarówno dla genów top ep300 jak i delta ep300
- tworzy heatmapy zmian przyłączania TF w zakresie +/- 2000 od pozycji amplitudy peaków dla NR3C1, dla dancyh top ep300 jak i delta ep300
- Wykresy enhancerów dla genów top ep300

![Kiku](PLOTS/lineplot_enhancer_range_top_ep300.jpeg)

![Kiku](PLOTS/heatmap_enhancer_top_ep300.jpeg)

- Wykresy enhancerów dla genów delta ep300

![Kiku](PLOTS/lineplot_enhancer_range_delta_ep300.jpeg)

![Kiku](PLOTS/heatmap_enhancer_delta_ep300.jpeg)



Do wyciągania danych dla Chip-seq z plików bigWig przygotowano skrypty:
- <a id="1">[1]</a>[bigwig_genomic_bucket500_extract_normalize_to_tsv.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigwig_genomic_bucket500_extract_normalize_to_tsv.sh)
    - wyciąga dane z plików bigwig o zakresie 20001, wynik zbiera po 500 pozycji i liczy średnią w celu zmniejszenia ilości danych i wyeliminowania szumów, w przypadku ostaniego baketowania liczy średnią z 501 pozycji
    - wyciąganie zakresów przy pomocy [bigWigToBedGraph](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigWigToBedGraph) [pobrane z](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/)
    - wyciąga dane dla wszystkich czynników transkrypcyjnych znajdujących się w pliku [chipseq-file-info.tsv](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/DATA/chipseq-file-info.tsv)
- [bigwig_genomic_amplitude_extract_normalize_to_tsv.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigwig_genomic_amplitude_extract_normalize_to_tsv.sh)
    - wyciąga amplitudy z podanych zakresów w plikach bigwig
    - wyciąga amplitudy przy pomocy [bigWigSummary](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigWigSummary) [pobrane z](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/)
    - wyciąga dane dla wszystkich czynników transkrypcyjnych znajdujących się w pliku [chipseq-file-info.tsv](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/DATA/chipseq-file-info.tsv)
- [bigwig_genomic_range_extract_normalize_to_tsv.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigwig_genomic_range_extract_normalize_to_tsv.sh)
    - wyciąga cały podany zakres z plików bigwig
    - wyciąganie zakresów przy pomocy [bigWigToBedGraph](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigWigToBedGraph) [pobrane z](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/)
    - wyciąga dane dla czterech czynników transkrypcyjnych: EP300, H3K27ac, H3K4me1, NR3C1
- [bigwig_genomic_range_extract_normalize_to_tsv_bucket10.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigwig_genomic_range_extract_normalize_to_tsv_bucket10.sh)
    - wyciąga dane z plików bigwig o zakresie 20001, wynik zbiera po 10 i liczy średnią w celu zmniejszenia ilości danych i wyeliminowania szumów, w przypadku ostatniego baketowania liczy średnią z 11 pozycji
    - wyciąganie zakresów przy pomocy [bigWigToBedGraph](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigWigToBedGraph) [pobrane z](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/)
    - - wyciąga dane dla czterech czynników transkrypcyjnych: EP300, H3K27ac, H3K4me1, NR3C1
- [bigwig_genomic_amplitude_extract_normalize_to_tsv.NR3C1-EP300.sh](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigwig_genomic_amplitude_extract_normalize_to_tsv.NR3C1-EP300.sh)
    - wyciąga amplitudy z podanych zakresów w plikach bigwig
    - wyciąga amplitudy przy pomocy [bigWigSummary](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigWigSummary) [pobrane z](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/)
    - wyciąga dane dla dwóch czynników transkrypcyjnych: NR3C1 i EP300
- wszystkie powyższe skrypty wykonują normalizację danych korzystając ze skryptów:
    - [bigWigToBedGraph](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bigWigToBedGraph) [pobrane z](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/)
    - [normalize_bedgraph.py](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/normalize_bedgraph.py) [pobrane z](https://github.com/porchard/normalize_bedgraph)
    - [bedGraphToBigWig](https://github.com/ippas/ifpan-chipseq-timecourse/blob/master/SCRIPTS/bedGraphToBigWig) [pobrane z](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/)

-----------------------------------------------Starsza wersja--------------------------------------------------------

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
~/ChIP-seq/DATA/enhancer_bigrange_delta_ep300.tsv

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


Przy pomocy komendy ściągnięto plik gtf:
```bash
wget  ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz
```
Następnie wybrano geny kodujące białka:
```bash
zcat ~/ChIP-seq/DATA/gencode.v36.annotation.gtf | 
    grep "gene_type \"protein_coding\"" | 
    awk '{print $10"\t"$16}' | 
    sort -u | 
    grep -v "[1-3];" | 
    sed 's/"//g' | 
    sed 's/;//g' | 
    sed 's/\./\t/' | 
    cut -f1,3 > ~/ifpan-chipseq-timecourse/DATA/Homo_sapiens.GRCh38.95.protein_coding.gtf
```

#####
Przy pomocy komendy wyciągnięto amplitudy enhancerów dla wszystkich genów:
```bash
 ~/ifpan-chipseq-timecourse/SCRIPTS/./bigwig_genomic_amplitude_extract_normalize_to_tsv.NR3C1-EP300.sh ~/ifpan-chipseq-timecourse/DATA/peaks_all_genes_without_promoters.tsv > ~/ifpan-chipseq-timecourse/DATA/peaks_all_genes_without_promoters_amplitude.tsv
```
~/ifpan-chipseq-timecourse/SCRIPTS/./search_peaks_all_genes.sh


~/ifpan-chipseq-timecourse/SCRIPTS/./search_peaks_all_genes.sh



---


