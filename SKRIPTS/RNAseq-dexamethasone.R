require(dplyr)
require(magrittr)
require(reshape2)
require(tidyr)
require(gplots)
require(dendextend)
require(data.table)
install.packages("matrixStats")
require(matrixStats)
FDR_THRESHOLD=0.001

raw.data <- read.delim("raw_data/raw_tar/extrakt/raw_macierz.txt", header = TRUE, stringsAsFactors = FALSE)
rownames(raw.data) <- raw.data[, 1]
samples <- read.delim("raw_data/raw_tar/extrakt/sample.info.txt", header = TRUE, stringsAsFactors = FALSE)
data <- raw.data[, 3:48] %>% as.matrix()
ID_ID.version_gene <- read.delim("raw_data/raw_tar/extrakt/ID_ID.version_gene.txt", header = TRUE, stringsAsFactors = FALSE)

pvalue <- anova(aov(data[1,] ~ samples$time))[[1,5]]
pvalue_matrix <- matrix(pvalue)
 
i <- 2
##do przerobienia z applay
while (i <= 60483) {
  pvalue <- anova(aov(data[i,] ~ samples$time))[[1,5]]
  pvalue_matrix <- rbind(pvalue_matrix, pvalue)
  i <- i+1
}
rm(i, pvalue)
colnames(pvalue_matrix) <- c("pvalue")
pvalue_matrix[is.na(pvalue_matrix)] <- 1

annotations <- raw.data[, -(3:48)]
annotations["ensemblid"] <- annotations %>% rownames %>% strsplit("[.]") %>% unlist %>%  magrittr::extract(seq(1, length(.), 2))
annotations$genename <- ID_ID.version_gene$Gene.name[match(annotations$ensemblid, ID_ID.version_gene$Gene.stable.ID)]  

results <- raw.data[, -(3:48)] 
results["ensemblid"] <- results %>% rownames %>% strsplit("[.]") %>% unlist %>%  magrittr::extract(seq(1, length(.), 2))
results <- cbind(results, pvalue_matrix)
results$genename <- ID_ID.version_gene$Gene.name[match(results$ensemblid, ID_ID.version_gene$Gene.stable.ID)]  
rm(pvalue_matrix)
#błąd pierwszego rodzaju
sum(results$fdr < FDR_THRESHOLD, na.rm = T)

number_signification_genes <- 737
#plik z danymi dla heatmap
tmp_data_for_plot <- data[order(results$pvalue)[1:number_signification_genes],] %>%
  apply(1, scale) %>% 
  t %>%
  apply(1, function(x, threshold){x[x > threshold] <- threshold; x[x < -threshold] <- -threshold; x}, threshold = 2.0) %>%
  t %>%
  {rownames(.) <- results$genename[order(results$pvalue)[1:737]]; .}  %>%
  {colnames(.) <- colnames(data); .}

number_clusters <- 2
tmp_cols_branches <- c("darkred", "forestgreen", "yellow")
tmp_dend1<- as.dist(1-cor(t(tmp_data_for_plot))) %>% hclust %>% as.dendrogram
tmp_dend1 <- color_branches(tmp_dend1, k = number_clusters)
col_labels <- get_leaves_branches_col(tmp_dend1)
col_labels <- col_labels[order(order.dendrogram(tmp_dend1))]

#heatmap
heatmap.2(tmp_data_for_plot,
          trace="none",
          margins =c(5,7),
          Colv = FALSE,
          col=bluered(20),
          dendrogram="row",      
          Rowv = tmp_dend1,  
          key.xlab = "Concentration (index)",
          distfun = function(x) as.dist(1-cor(t(x))),
          RowSideColors = col_labels, # to add nice colored strips        
          colRow = col_labels # to add nice colored labels - only for qplots 2.17.0 and higher
)
tmp_dend <- tmp_data_for_plot %>% t %>% {1 - cor(.)} %>% as.dist %>% 
  hclust %>% 
  as.dendrogram #%>%  
  #set(., "labels_to_character") %>% 
  #color_branches(k=number_clusters)

#geny z informacją do której gałęzi należą
cutrees <- cutree(tmp_dend, k = number_clusters) %>% as.matrix(.) %>%
{x <- .; names = names(x); names[is.na(names)] <- "noname"; names(x) <- names; x} %>%
  as.data.frame()
cutrees <- cbind(rownames(cutrees), cutrees) 
colnames(cutrees) <- NULL
colnames(cutrees) <- c("gene_name", "branches")
rownames(cutrees) <- NULL
cutrees <- cutrees %>% as.data.frame() %>% na.omit()

#Wykres liniowy
tmp_data_for_line_plot <- melt(tmp_data_for_plot) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  separate(., Var2, c("start", "end"), sep="[_]") %>% 
  .[, c("Var1", "start", "value")] %>% 
  group_by(Var1, time = factor(start, levels=c("t00","t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12"))) %>% 
  summarize(mean = mean(value)) %>% 
  arrange(time)
tmp_data_for_line_plot$number_branches <- cutrees$branches[match(cutrees$gene_name, tmp_data_for_line_plot$Var1)]
ggplot(tmp_data_for_line_plot, aes(x = time, y = mean, color = factor(number_branches))) + 
  geom_line(aes(group = Var1), alpha = 0.05) + 
  geom_smooth(aes(group = number_branches), se = FALSE, size =2) + 
  theme(legend.position="bottom")


#Tworzenie Barplot dla wybranego genu
data.frame(exprs = data[which(results$genename == "NYAP1"),], 
           time = factor(samples$time, 
                         levels = c("0", "30", "60", "120", "180", "240", "300", "360", "420", "480", "600", "720"))) %>% 
  group_by(time) %>% 
  summarise(mean = mean(exprs), sd = sd(exprs)) %>%
  {ggplot(., aes(x=time, y=mean)) + 
      geom_bar(position=position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))}
#usuwanie zbędnych ziennych
rm(tmp_cols_branches)
rm(tmp_data_for_line_plot)
rm(tmp_data_for_plot)
rm(number_clusters)
rm(ID_ID.version_gene)
rm(tmp_dend1)
rm(tmp_dend)


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
random_gene <- random_gene %>% mutate(pos=gene.start.bp * (strand == 1) + gene.end.bp * (strand == -1)) %>% mutate(start = pos - 10000, end = pos + 10001) %>% select(chromosome, start, end, genename)
#zapisanie do pliku
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
tmp_buckets <- rbind(bigtable_bucket_normalize, random_gene_normalize_bucket)
#rbind(bigtable_bucket_normalize, random_gene_normalize_bucket) %>%
#before normalize
tmp_buckets %>% 
  group_by(bucket, number_branches, time, factor_name) %>%
  summarize(value = mean(value)) %>%
  ggplot(., aes(x = as.numeric(bucket)*500, y = value, color = factor(number_branches))) + 
  geom_line(size = 0.5) +
  facet_grid(factor_name ~ time, scales = "free_y") + 
  theme(legend.position = "bottom") +
  ggtitle("signification_and_random_gene")

#plot whit normalize in ggplot
  tmp_buckets %>% 
  group_by(bucket, number_branches, time, factor_name) %>%
  summarize(value = mean(value)) %>%
  mutate(control = ifelse(number_branches == 3, 1, 0)) %>% 
  group_by(time, factor_name)  %>% 
  mutate(max = mean(value * control)) %>% 
  as.data.frame() %>%
  mutate(value = value / max) %>%
  ggplot(., aes(x = as.numeric(bucket)*500, y = value, color = factor(number_branches))) + 
  geom_line(size = 0.5) +
  facet_grid(factor_name ~ time, scales = "free_y") + 
  theme(legend.position = "bottom") +
  ggtitle("signification_and_random_gene_afternormalizeggplot")

#plot for change relative
tmp_buckets %>% 
  group_by(bucket, number_branches, time, factor_name) %>%
  summarize(value = mean(value)) %>%
  mutate(control = ifelse(number_branches == 3, 1, 0)) %>% 
  group_by(time, factor_name, bucket)  %>% 
  mutate(max = max(value * control)) %>% 
  as.data.frame() %>%
  mutate(value = value / max) %>%
  ggplot(., aes(x = as.numeric(bucket)*500, y = value, color = factor(number_branches))) + 
  geom_line(size = 0.5) +
  facet_grid(factor_name ~ time, scales = "free_y") + 
  theme(legend.position = "bottom") +
  ggtitle("signification_and_random_gene_relative_change")

rm(tmp_buckets)
rm(tmp_esimbl, tmp_namegenes)