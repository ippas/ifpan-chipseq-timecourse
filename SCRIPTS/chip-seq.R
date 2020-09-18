#wczytanie pliku pobranego z esmbl, z nazwami genów...
gene_chromosome_start_end_strand <- read.delim("/home/mateusz/ChIP-seq/gene_chromosome_start_end_strand.txt", header = TRUE, stringsAsFactors = FALSE)
#przygotowanie data_frame dla signification_gene
signification_gene <- annotations[order(results$pvalue)[1:number_signification_genes],]
signification_gene$gene.start.bp <- gene_chromosome_start_end_strand$Gene.start..bp.[match(signification_gene$ensemblid,gene_chromosome_start_end_strand$Gene.stable.ID)]
signification_gene$gene.end.bp <- gene_chromosome_start_end_strand$Gene.end..bp.[match(signification_gene$ensemblid,gene_chromosome_start_end_strand$Gene.stable.ID)]
signification_gene$strand <- gene_chromosome_start_end_strand$Strand[match(signification_gene$ensemblid,gene_chromosome_start_end_strand$Gene.stable.ID)]
signification_gene$chromosome <- gene_chromosome_start_end_strand$Chromosome.scaffold.name[match(signification_gene$ensemblid,gene_chromosome_start_end_strand$Gene.stable.ID)]
#usunięcie wierwszy z brakującymi nazwami genów
signification_gene <- na.omit(signification_gene)
#przygotowanie blików BAM do wyciągania peaków, up_regulated_gene i down_regulated_gene
regulated_gene <- signification_gene
regulated_gene$branch <- cutrees$branches[match(regulated_gene$genename, cutrees$gene_name)]
up_regulated_gene <- regulated_gene %>% 
  mutate(pos=gene.start.bp * (strand == 1) + gene.end.bp * (strand == -1)) %>% 
  mutate(start = pos - 100000, end = pos + 100001) %>% 
  filter(branch == 1) %>% 
  select(chromosome, start, end, genename)
up_regulated_gene$chromosome <- sub("^", "chr", up_regulated_gene$chromosome)
colnames(up_regulated_gene) <- NULL
down_regulated_gene <- regulated_gene %>% 
  mutate(pos=gene.start.bp * (strand == 1) + gene.end.bp * (strand == -1)) %>% 
  mutate(start = pos - 100000, end = pos + 100001) %>% 
  filter(branch == 2) %>% 
  select(chromosome, start, end, genename)
down_regulated_gene$chromosome <- sub("^", "chr", down_regulated_gene$chromosome)
colnames(down_regulated_gene) <- NULL
fwrite(up_regulated_gene, "/home/mateusz/ChIP-seq/up_regulated_gene.BED",sep="\t")
fwrite(down_regulated_gene, "/home/mateusz/ChIP-seq/down_regulated_gene.BED",sep="\t")
signification_gene <- signification_gene %>% mutate(pos=gene.start.bp * (strand == 1) + gene.end.bp * (strand == -1)) %>% mutate(start = pos - 10000, end = pos + 10001) %>% select(chromosome, start, end, genename, ensemblid) #%>% head
#zapisanie signification_gene do pliku, wyszukanie, normalizacja w bashu z wykorzystaniem pliku signification_gene.txt
fwrite(signification_gene,"/home/mateusz/ChIP-seq/signification_gene.txt", sep="\t",col.names = TRUE, row.names = TRUE)

###średnia ilość transkryptu
#wczytanie bliku z długościami transkryptu pobranymi z esembl
transcript_length <- as.data.frame(read.delim("/home/mateusz/ChIP-seq/transcript_length.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE))
#średnia długość transkryptu, obliczanie
transcript_length <- transcript_length %>% 
  group_by(.$Gene.stable.ID) %>% 
  summarize(Transcript.length..including.UTRs.and.CDS. = median(Transcript.length..including.UTRs.and.CDS.))
#tworzenie tabeli dla istotnych genów
RPKM_signification_gene <- data[order(results$pvalue)[1:number_signification_genes],]
rownames(RPKM_signification_gene) <- rownames(RPKM_signification_gene) %>% 
  as.data.frame() %>% 
  separate(., ., c("start", "end"), sep ="[.]") %>% .[, 1]
RPKM_signification_gene <- as.data.frame(RPKM_signification_gene)
#dodanie do tabeli kolumny z długością transkryptu i usunięcie wierszy dla których nie ma długości
RPKM_signification_gene$length_transcript <- transcript_length$Transcript.length..including.UTRs.and.CDS.[match(row.names(RPKM_signification_gene), transcript_length$`.$Gene.stable.ID`)] 
RPKM_signification_gene <- na.omit(RPKM_signification_gene)
#utworzenie obiektów z sumą kolumn, nazwami wierszy i dla długości transkryptów
tmp_colums_sum <- colSums(RPKM_signification_gene[,1:46], na.rm = FALSE, dims =1) %>% as.matrix()
tmp_rownames <- rownames(RPKM_signification_gene)
tmp_transcript_length <- RPKM_signification_gene[,47] %>% as.matrix()
#oblicznia
RPKM_signification_gene <- {RPKM_signification_gene[,1:46] / rep(tmp_colums_sum, each = nrow(RPKM_signification_gene))}
RPKM_signification_gene <- RPKM_signification_gene * 10^3 * 10^6
RPKM_signification_gene <- RPKM_signification_gene / tmp_transcript_length
RPKM_signification_gene <- as.data.frame(RPKM_signification_gene)
#tworzenie kolumny ze średnią dla kaźdego wiersza
RPKM_signification_gene <- mutate(RPKM_signification_gene, mean = rowSums(RPKM_signification_gene) %>% as.matrix() /46)
rownames(RPKM_signification_gene) <- tmp_rownames
RPKM_signification_gene$genename <- signification_gene$genename[match(row.names(RPKM_signification_gene), signification_gene$ensemblid)]
rm(tmp_rownames, tmp_colums_sum, tmp_transcript_length)
RPKM_signification_gene$number_branches <- cutrees$branches[match(RPKM_signification_gene$genename, cutrees$gene_name)]
#wykres dla RPKM
RPKM_signification_gene[, c("mean","genename","number_branches")] %>% mutate(logmean = log2(mean)) %>% group_by(number_branches) %>% {ggplot(.,aes(x = logmean)) + geom_histogram()} 

########Średnia ilość transkryptu dla wszystkich genów
RPKM_all_transcript <- data
RPKM_all_transcript <- as.data.frame(RPKM_all_transcript) %>% 
  mutate(., ensemblid = {rownames(data) %>%
      as.data.frame() %>% separate(., ., c("start", "end"), sep ="[.]") %>% 
      .[, 1]}) %>% mutate(genename = annotations$genename) %>% 
  na.omit() 
tmp_esimbl <- RPKM_all_transcript[, 47]
tmp_namegenes <- RPKM_all_transcript[, 48]
RPKM_all_transcript <- RPKM_all_transcript[,1:46] / rep({RPKM_all_transcript[, 1:46] %>% colSums(., na.rm = FALSE, dims = 1)}, each =nrow(RPKM_all_transcript)) * 10^6 * 10^3
RPKM_all_transcript <- mutate(RPKM_all_transcript, genename = tmp_namegenes)
RPKM_all_transcript <- mutate(RPKM_all_transcript, esimbl = tmp_esimbl)
row.names(RPKM_all_transcript) <- tmp_esimbl
RPKM_all_transcript$length_transcript <- transcript_length$Transcript.length..including.UTRs.and.CDS.[match(row.names(RPKM_all_transcript), transcript_length$`.$Gene.stable.ID`)]
RPKM_all_transcript <- RPKM_all_transcript[,1:46] / RPKM_all_transcript[,49]
RPKM_all_transcript <- mutate(RPKM_all_transcript, mean = rowSums(RPKM_all_transcript) %>% as.matrix() /46)
RPKM_all_transcript <- mutate(RPKM_all_transcript, genename = tmp_namegenes)
#wyciąganie losowej próby
random_gene <- RPKM_all_transcript[{RPKM_all_transcript$mean > 2 & RPKM_all_transcript$mean < 8192},] %>%
  .[sample(nrow(.), 1000), 48] %>% 
  as.data.frame() 
colnames(random_gene) <- c("genename")

#przygotowanie tabeli
random_gene$gene.start.bp <- gene_chromosome_start_end_strand$Gene.start..bp.[match(random_gene$genename ,gene_chromosome_start_end_strand$Gene.name)]
random_gene$gene.end.bp <- gene_chromosome_start_end_strand$Gene.end..bp.[match(random_gene$genename,gene_chromosome_start_end_strand$Gene.name)]
random_gene$strand <- gene_chromosome_start_end_strand$Strand[match(random_gene$genename, gene_chromosome_start_end_strand$Gene.name)]
random_gene$chromosome <- gene_chromosome_start_end_strand$Chromosome.scaffold.name[match(random_gene$genename ,gene_chromosome_start_end_strand$Gene.name)]
random_gene <- na.omit(random_gene)
#objekt dla genóœ randomowych do peakóœ
random_gene_into_peak <- random_gene %>% mutate(pos=gene.start.bp * (strand == 1) + gene.end.bp * (strand == -1)) %>% mutate(start = pos - 100000, end = pos + 100001) %>% select(chromosome, start, end, genename)  %>% 
  mutate(start = if_else(start < 0, 0, start)) 
random_gene_into_peak$chromosome <- sub("^", "chr", random_gene_into_peak$chromosome)
colnames(random_gene_into_peak) <- NULL
fwrite(random_gene_into_peak,"/home/mateusz/ChIP-seq/random_gene_into_peak.tsv", sep="\t",col.names = TRUE, row.names = FALSE) 

random_gene <- random_gene %>% mutate(pos=gene.start.bp * (strand == 1) + gene.end.bp * (strand == -1)) %>% mutate(start = pos - 10000, end = pos + 10001) %>% select(chromosome, start, end, genename)

#zapisanie do pliku
tmp_random_gene <- random_gene
fwrite(random_gene,"/home/mateusz/ChIP-seq/random_gene.txt", sep="\t",col.names = TRUE, row.names = TRUE)
#dodanie długości transkryptów
rm(tmp_namegenes, tmp_esimbl)
rm(transcript_length)
#dla jednego genu po normalizacji
bigtable_bucket_normalize <- as.data.frame(read.delim("/home/mateusz/ChIP-seq/bigtablebucket_normalize.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE))
colnames(bigtable_bucket_normalize) <- c("factor_name", "time", "gene_name", 1:40)
bigtable_bucket_normalize  %>% 
  gather(., "bucket", "value", -c("factor_name", "time", "gene_name")) %>%
  mutate(number_branches = cutrees$branches[match(.$gene_name, cutrees$gene_name)]) %>% 
  group_by(bucket, number_branches, time, factor_name) %>%
  summarize(value = mean(value)) %>%  
  ggplot(., aes(x = as.numeric(bucket)*500, y = value, color = factor(number_branches))) + 
  geom_line(size = 0.5) +
  facet_grid(factor_name ~ time, scales = "free_y") + 
  theme(legend.position = "bottom") +
  ggtitle("big_table_normalize")

bigtable_bucket_normalize <- bigtable_bucket_normalize %>%
  gather(., "bucket", "value", -c("factor_name", "time", "gene_name")) %>%
  mutate(number_branches = cutrees$branches[match(.$gene_name, cutrees$gene_name)])

random_gene_normalize_bucket <- as.data.frame(read.delim("/home/mateusz/ChIP-seq/random_gene_normalize_bucket.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE))
colnames(random_gene_normalize_bucket) <- c("factor_name", "time", "gene_name", 1:40)
random_gene_normalize_bucket %>%
  gather(., "bucket", "value", -c("factor_name", "time", "gene_name")) %>%
  group_by(bucket, time, factor_name) %>%
  summarize(value = mean(value)) %>%  
  ggplot(., aes(x = as.numeric(bucket)*500, y = value)) + 
  geom_line(size = 0.5) +
  facet_grid(factor_name ~ time, scales = "free_y") + 
  theme(legend.position = "bottom") +
  ggtitle("random_gene_normalize_bucket")

random_gene_normalize_bucket <- random_gene_normalize_bucket %>%
  gather(., "bucket", "value", -c("factor_name", "time", "gene_name")) %>% 
  mutate(number_branches = rep(3, 21089600))

rbind(bigtable_bucket_normalize, random_gene_normalize_bucket) %>%
  group_by(bucket, number_branches, time, factor_name) %>%
  summarize(value = mean(value)) %>%  
  ggplot(., aes(x = as.numeric(bucket)*500, y = value, color = factor(number_branches))) + 
  geom_line(size = 0.5) +
  facet_grid(factor_name ~ time, scales = "free_y") + 
  theme(legend.position = "bottom") +
  ggtitle("importantant_and_random_gene")
  
rm(tmp_esimbl, tmp_namegenes)

tmp_data_EP300_peak1 <- read.delim("/home/mateusz/ChIP-seq/tmp_data_normalize_peak1_EP300.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(tmp_data_EP300_peak1) <- c("factor_name", "time", 1:2805)
tmp_data_EP300_peak1 %>% gather(., "position", "value", -c("factor_name", "time")) %>% group_by(position, time) %>% 
  summarize(mean = mean(value)) %>% as.data.frame %>% mutate(time = as.factor(time), position = as.numeric(position)) %>% 
  ggplot(., aes(x = position, y = mean)) + geom_line(size = 0.5) + facet_grid(. ~ time, scales = "free_y") + 
  ggtitle("pierwszy_peak_dla_EP300")

tmp_normalize_data_EP300_peak1 <- read.delim("/home/mateusz/ChIP-seq/tmp_data_normalize_peak1_EP300.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(tmp_normalize_data_EP300_peak1) <- c("factor_name", "time", 1:2805)

tmp_normalize_data_EP300_peak1 %>% 
  gather(., "position", "value", -c("factor_name", "time")) %>% 
  group_by(position, time) %>% 
  summarize(mean = mean(value)) %>% 
  as.data.frame %>% 
  mutate(time = as.factor(time), position = as.numeric(position)) %>% 
  ggplot(., aes(x = position, y = mean)) + geom_line(size = 0.5) + facet_grid(. ~ time, scales = "free_y") + 
  ggtitle("pierwszy_peak_dla_EP300_znormalizowany")





###peaks
peaks_NR3C1_normalize_file1_time60 <- as.data.frame(read.delim("/home/mateusz/ChIP-seq/peaks_NR3C1_normalize_file1_time60.tsv", header = FALSE, sep = "\t", fill = TRUE, col.names = c ("TF", "time", "gene_name", "file", 1:3707), stringsAsFactors = FALSE))
peaks_NR3C1_amplitude<- peaks_NR3C1_normalize_file1_time60 %>% 
  mutate(amplitude = {rowMaxs({peaks_NR3C1_normalize_file1_time60[,5:3711] %>% 
  as.matrix()}, na.rm = T)}) %>% 
  select(TF, time, gene_name, file, amplitude)

peaks_NR3C1_amplitude %>% group_by(time) %>% ggplot(., aes(s = time, y = amplitude)) + geom_boxplot()  + facet_grid(. ~ time)
peaks_NR3C1_amplitude %>% 
  group_by(time, gene_name, file) %>% 
  summarise(max = max(amplitude)) %>% 
  ggplot(., aes(s = time, y = max)) + 
  geom_boxplot()  + 
  facet_grid(. ~ time) + 
  ggtitle("upregulated_gene_peak_NR3C1_file1")
#random peaks NR3C1 .rm 
random_peaks_NR3C1_time60_file1_normalize <- as.data.frame(read.delim("/home/mateusz/ChIP-seq/random_peaks_NR3C1_time60_file1_normalize.tsv", header = FALSE, sep = "\t", fill = TRUE, col.names = c ("TF", "time", "gene_name", "file", 1:3096), stringsAsFactors = FALSE))
random_peaks_NR3C1_time60_file1_amplitude <- random_peaks_NR3C1_time60_file1_normalize %>%
  mutate(amplitude = {rowMaxs({random_peaks_NR3C1_time60_file1_normalize[5:3100] %>%
    as.matrix()}, na.rm = T)}) %>% 
  select(TF, time, gene_name, file, amplitude)

random_peaks_NR3C1_time60_file1_amplitude %>% group_by(time) %>% 
  ggplot(., aes(s = time, y = amplitude)) + 
  geom_boxplot()  + 
  facet_grid(. ~ time) +
  ggtitle("random_peak_NR3C1_file1")

###peaki allTF NR3C1 allfile time60
upregulated_peaks_NR3C1_allTF_amplitude_time60 <- as.data.frame(read.delim("/home/mateusz/ChIP-seq/upregulated_peaks_NR3C1_allTF_amplitude_time60.tsv", header = FALSE, sep = "\t", fill = TRUE, col.names = c ("chromosome", "start_peak", "end_peak", "gene_name", "time", "TF", "file", "amplitude"), stringsAsFactors = FALSE))


upregulated_peaks_NR3C1_allTF_amplitude_time60 %>% 
  group_by(start_peak, gene_name, time, TF) %>%
  summarise(mean = mean(amplitude)) %>% 
  #filter(mean >= 25) %>%
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_boxplot(aes(group = cut_width(as.factor(time), width = 0.2))) +
  facet_wrap(TF ~ ., ncol = 4) +
  ggtitle("uporegulated_peak_NR3C1_time60")


upregulated_peaks_NR3C1_allTF_amplitude_time60 %>% 
  #filter(TF == "NR3C1") %>%
  group_by(start_peak, gene_name, time, TF) %>%
  summarise(mean_amplitude = mean(amplitude)) %>%
  group_by(time, TF) %>%
  summarize(mean = mean(mean_amplitude), sd = sd(mean_amplitude)) %>%
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(TF ~ ., ncol =4) +
  ggtitle("uporegulated_peak_NR3C1_allTF_time60")



upregulated_peaks_NR3C1_allTF_amplitude_time60 %>%
  #filter(gene_name == "SCNN1A", TF == "RAD21") %>%
  group_by(start_peak, gene_name, time, TF) %>%
  summarise(mean_alltime_amplitude = mean(amplitude)) %>%
  group_by(time, TF) %>%
  summarize(mean = mean(mean_amplitude), sd = sd(mean_amplitude)) %>%
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(TF ~ ., ncol =4) +
  ggtitle("uporegulated_peak_NR3C1_allTF_time60")

  
#wczytanie downregulowanych peaków
downregulated_peaks_NR3C1_allTF_amplitude_time60 <- as.data.frame(read.delim("/home/mateusz/ChIP-seq/downregulated_peaks_NR3C1_allfile_time60_allTF_normalize.tsv", 
                                                                             header = FALSE, 
                                                                             sep = "\t", 
                                                                             fill = TRUE, 
                                                                             col.names = c ("chromosome", "start_peak", "end_peak", "gene_name", "time", "TF", "file", "amplitude"), 
                                                                             stringsAsFactors = FALSE))

#wczytanie random1000 peaków
random1000_peaks_NR3C1_allTF_amplitude_time60 <- as.data.frame(read.delim("/home/mateusz/ChIP-seq/random1000_peaks_NR3C1_allTF_amplitude_time60.tsv", 
                                                                          header = FALSE, 
                                                                          sep = "\t", 
                                                                          fill = TRUE, 
                                                                          col.names = c ("chromosome", "start_peak", "end_peak", "gene_name", "time", "TF", "file", "amplitude"), 
                                                                          stringsAsFactors = FALSE))



###wykresy dla największych peaków
#upregulowane
barplot_upregulated_amplitude_peak_allTF_NR3C1_time60 <- merge(x = upregulated_peaks_NR3C1_allTF_amplitude_time60, y = {upregulated_peaks_NR3C1_allTF_amplitude_time60 %>% 
      #filter(gene_name == "SCNN1A", TF == "RAD21") %>%
      group_by(start_peak, gene_name, TF) %>%
      summarise(tmp_mean_alltime_amplitude = mean(amplitude))}, by = c("start_peak", "gene_name", "TF"), all = TRUE) %>% 
  group_by(gene_name, time, TF, tmp_mean_alltime_amplitude) %>%
  summarize(mean_amplitude_gene = mean(amplitude)) %>%
  filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
  #summarize(mean = mean(mean_amplitude_gene), sd = sd(mean_amplitude_gene)) %>%
  group_by(time, TF) %>%
  summarize(mean = mean(mean_amplitude_gene), sd = sd(mean_amplitude_gene)) %>%
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(TF ~ ., ncol =4) +
  ggtitle("barplot_upregulated_amplitude_peak_allTF_NR3C1_time60")

boxplot_upregulated_amplitude_peak_allTF_NR3C1_time60 <- merge(x = upregulated_peaks_NR3C1_allTF_amplitude_time60, y = {upregulated_peaks_NR3C1_allTF_amplitude_time60 %>% 
    group_by(start_peak, gene_name, TF) %>%
    summarise(tmp_mean_alltime_amplitude = mean(amplitude))}, by = c("start_peak", "gene_name", "TF"), all = TRUE) %>% 
  group_by(gene_name, time, TF, tmp_mean_alltime_amplitude) %>%
  filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
  group_by(start_peak, gene_name, time, TF) %>%
  summarise(mean = mean(amplitude)) %>% 
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_boxplot(aes(group = cut_width(as.factor(time), width = 0.2))) +
  facet_wrap(TF ~ ., ncol = 4) +
  ggtitle("boxplot_upregulated_amplitude_peak_allTF_NR3C1_time60")

boxplot_freey_upregulated_amplitude_peak_allTF_NR3C1_time60 <- merge(x = upregulated_peaks_NR3C1_allTF_amplitude_time60, y = {upregulated_peaks_NR3C1_allTF_amplitude_time60 %>% 
    group_by(start_peak, gene_name, TF) %>%
    summarise(tmp_mean_alltime_amplitude = mean(amplitude))}, by = c("start_peak", "gene_name", "TF"), all = TRUE) %>% 
  group_by(gene_name, time, TF, tmp_mean_alltime_amplitude) %>%
  filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
  group_by(start_peak, gene_name, time, TF) %>%
  summarise(mean = mean(amplitude)) %>% 
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_boxplot(aes(group = cut_width(as.factor(time), width = 0.2))) +
  facet_wrap(TF ~ ., ncol = 4, scales = "free_y") +
  ggtitle("boxplot_free-y_upregulated_amplitude_peak_allTF_NR3C1_time60")


###downregulated peak

boxplot_downregulated_amplitude_peak_allTF_NR3C1_time60 <- merge(x = downregulated_peaks_NR3C1_allTF_amplitude_time60, y = {downregulated_peaks_NR3C1_allTF_amplitude_time60 %>% 
    group_by(start_peak, gene_name, TF) %>%
    summarise(tmp_mean_alltime_amplitude = mean(amplitude))}, by = c("start_peak", "gene_name", "TF"), all = TRUE) %>% 
  group_by(gene_name, time, TF, tmp_mean_alltime_amplitude) %>%
  filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
  group_by(start_peak, gene_name, time, TF) %>%
  summarise(mean = mean(amplitude)) %>% 
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_boxplot(aes(group = cut_width(as.factor(time), width = 0.2))) +
  facet_wrap(TF ~ ., ncol = 4) +
  ggtitle("boxplot_downregulated_amplitude_peak_allTF_NR3C1_time60")

boxplot_freey_downregulated_amplitude_peak_allTF_NR3C1_time60 <- merge(x = downregulated_peaks_NR3C1_allTF_amplitude_time60, y = {downregulated_peaks_NR3C1_allTF_amplitude_time60 %>% 
    group_by(start_peak, gene_name, TF) %>%
    summarise(tmp_mean_alltime_amplitude = mean(amplitude))}, by = c("start_peak", "gene_name", "TF"), all = TRUE) %>% 
  group_by(gene_name, time, TF, tmp_mean_alltime_amplitude) %>%
  filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
  group_by(start_peak, gene_name, time, TF) %>%
  summarise(mean = mean(amplitude)) %>% 
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_boxplot(aes(group = cut_width(as.factor(time), width = 0.2))) +
  facet_wrap(TF ~ ., ncol = 4, scales = "free_y") +
  ggtitle("boxplot_free-y_downregulated_amplitude_peak_allTF_NR3C1_time60")

barplot_downregulated_amplitude_peak_allTF_NR3C1_time60 <- merge(x = downregulated_peaks_NR3C1_allTF_amplitude_time60, y = {downregulated_peaks_NR3C1_allTF_amplitude_time60 %>% 
    group_by(start_peak, gene_name, TF) %>%
    summarise(tmp_mean_alltime_amplitude = mean(amplitude))}, by = c("start_peak", "gene_name", "TF"), all = TRUE) %>% 
  group_by(gene_name, time, TF, tmp_mean_alltime_amplitude) %>%
  summarize(mean_amplitude_gene = mean(amplitude)) %>%
  filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
  group_by(time, TF) %>%
  summarize(mean = mean(mean_amplitude_gene), sd = sd(mean_amplitude_gene)) %>%
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(TF ~ ., ncol =4) +
  ggtitle("barplot_downregulated_amplitude_peak_allTF_NR3C1_time60")

#random1000 wykresy
# boxplot_random1000_amplitude_peak_allTF_NR3C1_time60 <- merge(x = random1000_peaks_NR3C1_allTF_amplitude_time60, y = {random1000_peaks_NR3C1_allTF_amplitude_time60 %>% 
#     group_by(start_peak, gene_name, TF) %>%
#     summarise(tmp_mean_alltime_amplitude = mean(amplitude))}, by = c("start_peak", "gene_name", "TF"), all = TRUE) %>% 
#   group_by(gene_name, time, TF, tmp_mean_alltime_amplitude) %>%
#   filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
#   group_by(start_peak, gene_name, time, TF) %>%
#   summarise(mean = mean(amplitude)) %>% 
#   ggplot(., aes(x = as.factor(time), y = mean)) +
#   geom_boxplot(aes(group = cut_width(as.factor(time), width = 0.2))) +
#   facet_wrap(TF ~ ., ncol = 4) +
#   ggtitle("boxplot_random1000_amplitude_peak_allTF_NR3C1_time60")

# boxplot_free-y_random1000_amplitude_peak_allTF_NR3C1_time60 <- merge(x = downregulated_peaks_NR3C1_allTF_amplitude_time60, y = {downregulated_peaks_NR3C1_allTF_amplitude_time60 %>% 
#     group_by(start_peak, gene_name, TF) %>%
#     summarise(tmp_mean_alltime_amplitude = mean(amplitude))}, by = c("start_peak", "gene_name", "TF"), all = TRUE) %>% 
#   group_by(gene_name, time, TF, tmp_mean_alltime_amplitude) %>%
#   summarize(mean_amplitude_gene = mean(amplitude)) %>%
#   filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
#   group_by(time, TF) %>%
#   summarize(mean = mean(mean_amplitude_gene), sd = sd(mean_amplitude_gene)) %>%
#   ggplot(., aes(x = as.factor(time), y = mean)) +
#   geom_bar(position=position_dodge(), stat="identity") +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
#   facet_wrap(TF ~ ., ncol =4, scales = "free_y") +
#   ggtitle("boxplot_free-y_random1000_amplitude_peak_allTF_NR3C1_time60")

#randomowe1000 peaków
boxplot_random1000_amplitude_peak_allTF_NR3C1_time60 <- random1000_peaks_NR3C1_allTF_amplitude_time60 %>% 
  group_by(start_peak, TF, time) %>%
  summarise(mean_amplitude = mean(amplitude)) %>%
  ggplot(., aes(x = as.factor(time), y = mean_amplitude)) +
  geom_boxplot(aes(group = cut_width(as.factor(time), width = 0.2))) +
  facet_wrap(TF ~ ., ncol = 4) +
  ggtitle("boxplot_random1000_amplitude_peak_allTF_NR3C1_time60")

boxplot_freey_random1000_amplitude_peak_allTF_NR3C1_time60 <- random1000_peaks_NR3C1_allTF_amplitude_time60 %>% 
  group_by(start_peak, TF, time) %>%
  summarise(mean_amplitude = mean(amplitude)) %>%
  ggplot(., aes(x = as.factor(time), y = mean_amplitude)) +
  geom_boxplot(aes(group = cut_width(as.factor(time), width = 0.2))) +
  facet_wrap(TF ~ ., ncol = 4, scales = "free_y") +
  ggtitle("boxplot_freey_random1000_amplitude_peak_allTF_NR3C1_time60")


barplot_random1000_amplitude_peak_allTF_NR3C1_time60 <- random1000_peaks_NR3C1_allTF_amplitude_time60 %>% 
  group_by(start_peak, TF, time) %>%
  summarise(mean_amplitude = mean(amplitude)) %>%
  group_by(time, TF) %>%
  summarize(mean = mean(mean_amplitude), sd = sd(mean_amplitude)) %>%
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(TF ~ ., ncol =4) +
  ggtitle("barplot_random1000_amplitude_peak_allTF_NR3C1_time60")
#zapisanie wykresów do plików
tmp_name_file <- c("barplot_upregulated_amplitude_peak_allTF_NR3C1_time60",
                   "boxplot_upregulated_amplitude_peak_allTF_NR3C1_time60", 
  "boxplot_freey_upregulated_amplitude_peak_allTF_NR3C1_time60", 
  "boxplot_downregulated_amplitude_peak_allTF_NR3C1_time60",
  "boxplot_freey_downregulated_amplitude_peak_allTF_NR3C1_time60", 
  "barplot_downregulated_amplitude_peak_allTF_NR3C1_time60",
  "boxplot_random1000_amplitude_peak_allTF_NR3C1_time60", 
  "boxplot_freey_random1000_amplitude_peak_allTF_NR3C1_time60",
  "barplot_random1000_amplitude_peak_allTF_NR3C1_time60") 
  

jpeg("/home/mateusz/ifpan-chipseq-timecourse/PLOTS/bboxplot_freey_random1000_amplitude_peak_allTF_NR3C1_time60.jpeg", width = 1400, height = 802)
print(boxplot_freey_random1000_amplitude_peak_allTF_NR3C1_time60)
dev.off()



for(tmp_file in tmp_name_file) { 
  jpeg(print(paste("/home/mateusz/ifpan-chipseq-timecourse/PLOTS/", tmp_file, ".jpeg", sep = "")),width = 1400, height = 802)
  print(tmp_file %>% get())
  dev.off()
  #print(print(tmp_file))
  tmp_file %>% as.name() %>% rm()
  }
# jpeg("/home/mateusz/ifpan-chipseq-timecourse/PLOTS/barplot_random1000_amplitude_peak_allTF_NR3C1_time60.jpeg",width = 1400, height = 802)
# print(barplot_random1000_amplitude_peak_allTF_NR3C1_time60)
# dev.off()
# rm(tmp_name_file)

random_gen_peaks_NR3C1_allTF_amplitude_time60 <- as.data.frame(read.delim("/home/mateusz/ChIP-seq/random_gen_peak_NR3C1_allTF_amplitude_time60.tsv", 
                                                                          header = FALSE, 
                                                                          sep = "\t", 
                                                                          fill = TRUE, 
                                                                          col.names = c ("chromosome", "start_peak", "end_peak", "gene_name", "time", "TF", "file", "amplitude"), 
                                                                          stringsAsFactors = FALSE))

barplot_random_gen_amplitude_peak_allTF_NR3C1_time60 <- merge(x = random_gen_peaks_NR3C1_allTF_amplitude_time60, y = {random_gen_peaks_NR3C1_allTF_amplitude_time60 %>% 
    #filter(gene_name == "SCNN1A", TF == "RAD21") %>%
    group_by(start_peak, gene_name, TF) %>%
    summarise(tmp_mean_alltime_amplitude = mean(amplitude))}, by = c("start_peak", "gene_name", "TF"), all = TRUE) %>% 
  group_by(gene_name, time, TF, tmp_mean_alltime_amplitude) %>%
  summarize(mean_amplitude_gene = mean(amplitude)) %>%
  filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
  #summarize(mean = mean(mean_amplitude_gene), sd = sd(mean_amplitude_gene)) %>%
  group_by(time, TF) %>%
  summarize(mean = mean(mean_amplitude_gene), sd = sd(mean_amplitude_gene)) %>%
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  facet_wrap(TF ~ ., ncol =4) +
  ggtitle("barplot_random_gen_amplitude_peak_allTF_NR3C1_time60")

boxplot_random_gen_amplitude_peak_allTF_NR3C1_time60 <- merge(x = random_gen_peaks_NR3C1_allTF_amplitude_time60, y = {random_gen_peaks_NR3C1_allTF_amplitude_time60 %>% 
    group_by(start_peak, gene_name, TF) %>%
    summarise(tmp_mean_alltime_amplitude = mean(amplitude))}, by = c("start_peak", "gene_name", "TF"), all = TRUE) %>% 
  group_by(gene_name, time, TF, tmp_mean_alltime_amplitude) %>%
  filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
  group_by(start_peak, gene_name, time, TF) %>%
  summarise(mean = mean(amplitude)) %>% 
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_boxplot(aes(group = cut_width(as.factor(time),o width = 0.2))) +
  facet_wrap(TF ~ ., ncol = 4) +
  ggtitle("boxplot_random_gen_amplitude_peak_allTF_NR3C1_time60")

boxplot_freey_random_gen_amplitude_peak_allTF_NR3C1_time60 <- merge(x = random_gen_peaks_NR3C1_allTF_amplitude_time60, y = {random_gen_peaks_NR3C1_allTF_amplitude_time60 %>% 
    group_by(start_peak, gene_name, TF) %>%
    summarise(tmp_mean_alltime_amplitude = mean(amplitude))}, by = c("start_peak", "gene_name", "TF"), all = TRUE) %>% 
  group_by(gene_name, time, TF, tmp_mean_alltime_amplitude) %>%
  filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
  group_by(start_peak, gene_name, time, TF) %>%
  summarise(mean = mean(amplitude)) %>% 
  ggplot(., aes(x = as.factor(time), y = mean)) +
  geom_boxplot(aes(group = cut_width(as.factor(time), width = 0.2))) +
  facet_wrap(TF ~ ., ncol = 4, scales = "free_y") +
  ggtitle("boxplot_freey_random_gen_peaks_NR3C1_allTF_amplitude_time60")

tmp_name_plot <- c("barplot_random_gen_amplitude_peak_allTF_NR3C1_time60",
                   "boxplot_random_gen_amplitude_peak_allTF_NR3C1_time60",
                   "boxplot_freey_random_gen_amplitude_peak_allTF_NR3C1_time60")

for(tmp_plot in tmp_name_plot) { 
  jpeg(print(paste("/home/mateusz/ifpan-chipseq-timecourse/PLOTS/", tmp_plot, ".jpeg", sep = "")),width = 1400, height = 802)
  print(tmp_plot %>% get())
  dev.off()
  #print(print(tmp_file))
  #rm(as.name(tmp_plot))
  #tmp_plot %>% as.name() %>% rm()
}
