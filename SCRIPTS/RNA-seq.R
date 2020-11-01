require(dplyr)
require(magrittr)
require(reshape2)
require(tidyr)
require(gplots)
require(dendextend)
require(data.table)
install.packages("matrixStats")
require(matrixStats)
install.packages("rebus")
require(rebus)
FDR_THRESHOLD=0.001
install.packages("BiocManager")
require(BiocManager)
BiocManager::install("preprocessCore")

raw.data <- read.delim("~/ifpan-chipseq-timecourse/DATA/raw_expression_matrix_dexamethasone.tsv",  
                       header = TRUE, 
                       stringsAsFactors = FALSE) %>% 
  set_rownames(.$Geneid)

samples <- read.delim("~/ifpan-chipseq-timecourse/DATA/sample.info.tsv", 
                      header = TRUE, 
                      stringsAsFactors = FALSE)

data <- raw.data[, 3:48] %>% 
  as.matrix()

raw.data[, 3:48] %>% 
  as.matrix() %>% 
  normalize.quantiles() %>% 
  set_rownames(., rownames(data)) %>% 
  set_colnames(colnames(data)) %>%
  as.data.frame() %>% 
  rownames_to_column(var="Geneid") %>% 
  gather(., "time.sample", "value.RNAseq", -c(Geneid)) %>%
  mutate(log.value.RNAseq = log2(value.RNAseq + 1)) %>%
  filter(log.value.RNAseq > 1) %>%
  mutate(number.sample = .$time.sample %>%
           str_split(., "_", n =2, simplify = TRUE) %>%
           .[,2]) %>%
  mutate(time = .$time.sample %>%
           str_split(., "_", n =2, simplify = TRUE) %>%
           .[,1]) %>%
           {ggplot(., aes(x=number.sample, y=log.value.RNAseq)) +
               geom_boxplot() + 
               facet_wrap(time ~., ncol = 4)}


raw.data[, 3:48] %>% 
  as.matrix() %>% 
  normalize.quantiles() %>% 
  set_rownames(., rownames(data)) %>% 
  set_colnames(colnames(data)) -> data


ID_ID.version_gene <- read.delim("~/ifpan-chipseq-timecourse/DATA/ID_ID.version_gene.tsv", 
                                 header = TRUE, 
                                 stringsAsFactors = FALSE)

tmp_pvalue_dataframe <- data %>% 
  apply(., 1, function(i){anova(aov(i~samples$time))[[1,5]]}) %>% 
  as.matrix() %>% 
  set_colnames("pvalue") %>%
  replace_na(., 1) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Geneid")

results <- raw.data[, -(3:48)] %>%
  mutate(Gene.stable.ID = .$Geneid %>% 
           str_split(., fixed("."), n = 2, simplify = TRUE) %>% 
           .[, 1]) %>%
  left_join(., ID_ID.version_gene, by = "Gene.stable.ID") %>%
  rename(ensemblid = Gene.stable.ID,
         gene.name = Gene.name) %>%
  select(-Gene.stable.ID.version) %>%
  left_join(., tmp_pvalue_dataframe, by = "Geneid") %>%
  select("Geneid", "ensemblid", "Length", "gene.name", "pvalue")
  
rm(tmp_pvalue_dataframe)
#results %>% mutate(fdr = p.adjust(results$pvalue, method = "fdr")) %>% filter(fdr < 0.001) %>% nrow()
number_signification_genes <- results %>% 
  mutate(fdr = p.adjust(results$pvalue, method = "fdr")) %>% 
  filter(fdr < 0.00000001) %>% 
  nrow()

results %>% 
  mutate(fdr = p.adjust(results$pvalue, method = "fdr")) %>% 
  filter(fdr < 0.00000001)  -> results.filtered

tmp_data_for_heatmap <- data[order(results$pvalue)[1:number_signification_genes],] %>%
  apply(1, scale) %>% 
  t %>%
  apply(1, function(x, threshold){x[x > threshold] <- threshold; x[x < -threshold] <- -threshold; x}, threshold = 2.0) %>%
  t %>%
  {rownames(.) <- results$gene.name[order(results$pvalue)[1:number_signification_genes]]; .}  %>%
  {colnames(.) <- colnames(data); .}

number_clusters <- 2

tmp_dend <- as.dist(1-cor(t(tmp_data_for_heatmap))) %>% 
  hclust %>% 
  as.dendrogram %>% 
  color_branches(., k = number_clusters)

tmp_col_labels <- get_leaves_branches_col(tmp_dend) %>% 
  .[order(order.dendrogram(tmp_dend))]


#Create heatmap and save into file
jpeg("~/ifpan-chipseq-timecourse/PLOTS/heatmap_significant_genes.jpeg", 
     width = 1400, 
     height = 802)

heatmap.2(tmp_data_for_heatmap,
          trace="none",
          margins =c(5,7),
          Colv = FALSE,
          col=bluered(20),
          dendrogram="row",      
          Rowv = tmp_dend,  
          key.xlab = "Concentration (index)",
          distfun = function(x) as.dist(1-cor(t(x))),
          RowSideColors = tmp_col_labels, 
          colRow = tmp_col_labels 
) 

dev.off()


gene_regulation <- tmp_data_for_plot %>% 
  t %>% 
  {1 - cor(.)} %>% 
  as.dist %>% 
  hclust %>% 
  as.dendrogram %>%
  cutree(., k = number_clusters) %>% 
  as.table() %>% 
  as.data.frame() %>%  
  set_colnames(c("gene.name", "number_regulation")) %>% 
  mutate(gene.regulation = recode(number_regulation, 
                                  `1` = "up-regulated", 
                                  `2` = "down-regulated")) %>%
  select(gene.name, gene.regulation)

#Create lineplot for up and down regulation significant genes and save into file
#Create lineplot for up and down regulation significant genes and save into file

#but first lets prepare random genes the same way

data %>% data.frame() %>%
  filter(rownames(data) %in% random$Geneid) %>%
  apply(1, scale) %>% 
  t %>% 
  apply(1, function(x, threshold){x[x > threshold] <- threshold; x[x < -threshold] <- -threshold; x}, threshold = 2.0) %>%
  t %>%
  {rownames(.) <- results$gene.name[match(rownames(data)[rownames(data) %in% random$Geneid], results$Geneid)]; .}  %>%
  {colnames(.) <- colnames(data); .} -> random.prepared 


random.prepared %>% melt() %>% 
  na.omit() %>% 
  as.data.frame() %>%
  mutate(Var2 = str_remove_all(Var2, pattern = "_rep" %R% one_or_more(DGT))) %>%
  mutate(Var2 = str_remove(Var2, pattern = START %R% "t")) %>% 
  mutate(Var2 = 60 * as.numeric(recode(Var2, 
                                       "00" = "0", 
                                       "05" = "0.5")))  %>%
  set_colnames(c("gene.name", "time", "value")) %>% 
  mutate(., gene.regulation = c(rep("random", nrow(.)))) %>%
  group_by(gene.name, time, gene.regulation) %>% 
  summarize(mean = mean(value)) %>%
  arrange(time) %>% data.frame() -> random.prepared



jpeg("~/ifpan-chipseq-timecourse/PLOTS/lineplot_up_down_regulation_significant_genes.jpeg", 
     width = 1400, 
     height = 802)

tmp_data_for_heatmap %>% 
  melt() %>% 
  na.omit() %>% 
  as.data.frame() %>%
  mutate(Var2 = str_remove_all(Var2, pattern = "_rep" %R% one_or_more(DGT))) %>%
  mutate(Var2 = str_remove(Var2, pattern = START %R% "t")) %>% 
  mutate(Var2 = 60 * as.numeric(recode(Var2, 
                                       "00" = "0", 
                                       "05" = "0.5")))  %>%
  set_colnames(c("gene.name", "time", "value")) %>%
  left_join(., gene_regulation, by = "gene.name") %>% 
  group_by(gene.name, time, gene.regulation) %>% 
  summarize(mean = mean(value)) %>%
  arrange(time) %>% na.omit() %>% 
  ggplot(., 
         aes(x = as.numeric(time), y = mean, color = gene.regulation)) +
  geom_line(data = random.prepared, aes(group = gene.name), alpha = 0.01, color="grey") +
  geom_smooth(data = random.prepared, aes(group = gene.regulation), se = FALSE, size = 1, color="grey") +
  geom_line(aes(group = gene.name), alpha = 0.01) +
  geom_smooth(aes(group = gene.regulation), se = FALSE, size = 1) +
  theme(legend.position = "bottom")

# tmp_data_for_heatmap %>% 
#   melt() %>% 
#   na.omit() %>% 
#   as.data.frame() %>% #head %>%
#   mutate(Var2 = str_remove_all(Var2, pattern = "_rep" %R% one_or_more(DGT))) %>%
#   mutate(Var2 = str_remove(Var2, pattern = START %R% "t")) %>% 
#   mutate(Var2 = 60 * as.numeric(recode(Var2, 
#                        "00" = "0", 
#                        "05" = "0.5")))  %>%
#   set_colnames(c("gene.name", "time", "value")) %>%
#   left_join(., gene_regulation, by = "gene.name") %>%
#   group_by(gene.name, time, gene.regulation) %>% 
#   summarize(mean = mean(value)) %>%
#   arrange(time) %>%
#   ggplot(., 
#        aes(x = as.factor(time), y = mean, color = gene.regulation)) +
#   geom_line(aes(group = gene.name), alpha = 0.05) +
#   geom_smooth(aes(group = gene.regulation), se = FALSE, size = 2) +
#   theme(legend.position = "bottom")

dev.off()
# 
# #Create Barplot for a single gene
# data.frame(exprs = data[which(results$genename == "NYAP1"),], 
#            time = factor(samples$time, 
#                          levels = c("0", "30", "60", "120", "180", "240", "300", "360", "420", "480", "600", "720"))) %>% 
#   group_by(time) %>% 
#   summarise(mean = mean(exprs), sd = sd(exprs)) %>%
#   ggplot(., aes(x=time, y=mean)) + 
#       geom_bar(position=position_dodge(), stat="identity") +
#       geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
#                     width=.2,
#                     position=position_dodge(.9))
# 
# dev.off()

# remove tmp variable
rm(tmp_data_for_heatmap,
   tmp_dend,
   tmp_col_labels)

gene_chromosome_start_end_strand <- read.delim("~/ifpan-chipseq-timecourse/DATA/gene_chromosome_start_end_strand.tsv", 
                                               header = TRUE, 
                                               stringsAsFactors = FALSE)

# results[order(results$pvalue)[1:number_signification_genes],] %>% 
#   left_join(., gene_chromosome_start_end_strand, by = c("ensemblid" = "Gene.stable.ID")) %>% 
#   mutate(pos=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% 
#   mutate(start = pos - 10000, end = pos + 10001) %>% 
#   select("ensemblid", "gene.name", "Chromosome.scaffold.name", "start", "end") %>%
#   set_colnames(c("ensemblid", "gene.name", "chromosome", "start", "end")) %>% 
#   left_join(., gene_regulation, by = "gene.name") %>% 
#   na.omit() %>% 
#   fwrite("~/ifpan-chipseq-timecourse/DATA/significant_genes_ensemblid_genename_chromosome_start_end.tsv", 
#          sep="\t", 
#          col.names = TRUE, 
#          row.names = FALSE)
  

tmp_transcript_length <- as.data.frame(read.delim("~/ifpan-chipseq-timecourse/DATA/transcript_length.tsv", 
                                                  header = TRUE, 
                                                  sep = "\t", 
                                                  stringsAsFactors = FALSE)) %>% 
  group_by(.$Gene.stable.ID) %>% 
  summarize(median = median(Transcript.length..including.UTRs.and.CDS.)) %>% 
  set_colnames(c("ensemblid", "median"))

jpeg("~/ifpan-chipseq-timecourse/PLOTS/histogram_significant_gene_logmean_transcriptlength.jpeg", 
     width = 1400, 
     height = 802)

data[order(results$pvalue)[1:number_signification_genes],] %>% 
  set_rownames(rownames(data[order(results$pvalue)[1:number_signification_genes],]) %>% 
                 str_split(., fixed("."), n = 2, simplify = TRUE) %>%
                 .[, 1]) %>% 
  as.data.frame() %>% 
  mutate(ensemblid = rownames(.)) %>%
  left_join(., tmp_transcript_length, by = "ensemblid") %>%
  na.omit()  %>%
  {.[,1:46] <- ((.[,1:46]/rep(colSums(.[,1:46]), each = nrow(.)))*10^3 * 10^6)/.$median; .} %>% 
  mutate(., mean = rowSums(.[,1:46])/46)  %>% 
  left_join(.,results[, c("ensemblid", "gene.name")], by = "ensemblid" ) %>%
  column_to_rownames(., var = "ensemblid") %>%  
  left_join(., gene_regulation, by = "gene.name") %>% 
  mutate(logmean = log2(mean))  %>% 
  {ggplot(.,aes(x = logmean)) + 
      geom_histogram() +
      facet_grid(~gene.regulation) +
      ggtitle("histogram_significant_gene_logmean_transcriptlength")
  } 

dev.off()

#random
data %>% 
  as.data.frame() %>% 
  mutate(., ensemblid = {row.names(.) %>% 
           #as.data.frame() %>% 
           str_split(., fixed("."), n = 2, simplify = TRUE) %>% 
           .[, 1]}) %>%
  left_join(., results[, c("ensemblid", "gene.name")], by = "ensemblid") %>%
  left_join(., tmp_transcript_length, by = "ensemblid") %>% 
  na.omit()  %>% 
  {.[,1:46] <- ((.[,1:46]/rep(colSums(.[,1:46]), each = nrow(.)))*10^3 * 10^6)/.$median; .} %>% 
  mutate(., mean = rowSums(.[,1:46])/46) %>% 
  .[{.$mean > 2 & .$mean < 8192},] %>%
  .[sample(nrow(.), 1000), ] %>% 
  left_join(gene_chromosome_start_end_strand, ., by = c(c("Gene.stable.ID"="ensemblid"),c("Gene.name"="gene.name"))) %>% 
  na.omit() %>% 
  mutate(logmean = log2(mean))  %>% 
  {ggplot(.,aes(x = logmean)) + 
      geom_histogram()}
      

#create file with significant and random genes
read.table("~/ifpan-chipseq-timecourse/DATA/random_genes_geneid_ensemblid_length_gene.name_pvalue_mean.expression.tsv", 
           header = TRUE,
           row.names = 1,
           stringsAsFactors = FALSE) %>%
  select(ensemblid, gene.name) %>%
  left_join(., gene_chromosome_start_end_strand, by = c("ensemblid" = "Gene.stable.ID")) %>%
  mutate(pos=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% 
  mutate(start = pos - 10000, end = pos + 10001) %>%
  select("ensemblid", "gene.name", "Chromosome.scaffold.name", "start", "end") %>%
  set_colnames(c("ensemblid", "gene.name", "chromosome", "start", "end")) %>%
  mutate(gene.regulation = "random") %>% na.omit() %>% 
  rbind(., {results[order(results$pvalue)[1:number_signification_genes],] %>% 
               left_join(., gene_chromosome_start_end_strand, by = c("ensemblid" = "Gene.stable.ID")) %>% 
               mutate(pos=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% 
               mutate(start = pos - 10000, end = pos + 10001) %>% 
               select("ensemblid", "gene.name", "Chromosome.scaffold.name", "start", "end") %>%
               set_colnames(c("ensemblid", "gene.name", "chromosome", "start", "end")) %>% 
               left_join(., gene_regulation, by = "gene.name") %>% 
               na.omit()}) %>% 
  mutate(start = if_else(start < 0, 0, start)) %>% 
  fwrite("~/ifpan-chipseq-timecourse/DATA/significant_and_random_genes_ensemblid_genename_chromosome_start_end.tsv", 
         sep="\t", 
         col.names = TRUE, 
         row.names = FALSE)


jpeg("~/ifpan-chipseq-timecourse/PLOTS/lineplot_significant_random_genes_normalized_bucket.jpeg", 
     width = 1400, 
     height = 802)

read.table("~/ChIP-seq/DATA/significant_random_genes_chip-seq_normalized_bucket_gene_chromosome_start_end_TF_time_file.tsv", 
           header = FALSE, 
           sep = "\t", 
           stringsAsFactors = FALSE) %>% 
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", 1:40)) %>% 
  gather(., "bucket.range", "value", -c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file")) %>% 
  group_by(bucket.range, time, TF, gene.regulation) %>%
  summarize(mean.value = mean(value)) %>%
  ggplot(., aes(x = as.numeric(bucket.range)*500, y = mean.value, color = as.factor(gene.regulation))) + 
     geom_line(size = 0.5) + 
     facet_grid(TF~time, scales = "free_y") +
     theme(legend.position = "bottom") +
     ggtitle("lineplot_significant_random_genes_normalized_bucket")

dev.off()


jpeg("~/ifpan-chipseq-timecourse/PLOTS/lineplot_significant_random_genes_normalized_bucket_relative_changes.jpeg", 
     width = 1400, 
     height = 802)

read.table("~/ChIP-seq/DATA/tmp_significant_random_genes_chip-seq_normalized_bucket_gene_chromosome_start_end_TF_time_file.tsv",
           header = FALSE,
           sep = "\t",
           stringsAsFactors = FALSE) %>%
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", 1:40)) %>%
  gather(., "bucket.range", "value", -c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file")) %>%
  group_by(bucket.range, time, TF, gene.regulation) %>%
  summarize(mean.value = mean(value)) %>%
  mutate(number.regulation=c("down-regulated"=1, "random"=3, "up-regulated"=2)) %>%
  mutate(control = ifelse(number.regulation == 3, 1, 0)) %>%
  group_by(time, TF, bucket.range)  %>% 
  mutate(max = max(mean.value * control)) %>% 
  mutate(relative.value = mean.value / max) %>%
  ggplot(., aes(x = as.numeric(bucket.range)*500, y = relative.value, color = as.factor(gene.regulation))) + 
    geom_line(size = 0.5) + 
    facet_grid(TF~time, scales = "free_y") +
    theme(axis.text.x = element_text(size=6),
          legend.position = "bottom") +
    ggtitle("lineplot_significant_random_genes_normalized_bucket_relative_changes")

dev.off()


### Enhancer peak ########################################################

tmp_significant_random_genes_peak_normalized_amplitude <- read.table("~/ChIP-seq/DATA/significant_random_genes_chip-seq_normalized_gene_chromosome_start-peak_end-peak_TF_time_file_amplitude.tsv", 
           header = FALSE, 
           sep = "\t", 
           stringsAsFactors = FALSE) %>% 
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", "amplitude"))
  
jpeg("~/ifpan-chipseq-timecourse/PLOTS/boxplot_significant_random_genes_strongest_peak.jpeg", 
     width = 1400, 
     height = 802)

#making boxplot for strongest peak for each gene
merge(x = tmp_significant_random_genes_peak_normalized_amplitude,
      y = {tmp_significant_random_genes_peak_normalized_amplitude %>% 
          group_by(gene.name, start.range, TF, gene.regulation) %>%
          summarise(tmp_mean_alltime_amplitude = mean (amplitude))}) %>%
  group_by(gene.name, TF, gene.regulation) %>%
  filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.max.peak = mean(amplitude)) %>% 
  ggplot(., aes(x = as.factor(time), y = log(mean.max.peak), color = gene.regulation)) +
  geom_boxplot(position = position_dodge(), outlier.size = 0) +  
  facet_wrap(TF ~ ., ncol = 4, scales = "free_y" ) +
  theme(legend.position = "bottom") 
  ggtitle("boxplot_significant_random_genes_strongest_peak")
  
dev.off()

jpeg("~/ifpan-chipseq-timecourse/PLOTS/boxplot_significant_random_genes_mean_peaks.jpeg", 
       width = 1400, 
       height = 802)

#making boxplot for all peaks for each gene
tmp_significant_random_genes_peak_normalized_amplitude %>%
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.max.peak = mean(amplitude)) %>% 
  ggplot(., aes(x = as.factor(time), y = log(mean.max.peak), color = gene.regulation)) +
  #geom_boxplot(aes(group = cut_width(as.factor(time),  width = 0.02)), position=position_dodge(), stat="identity") +
  geom_boxplot(position = position_dodge(), outlier.size = 0) +  
  facet_wrap(TF ~ ., ncol = 4, scales = "free_y" ) +
  theme(legend.position = "bottom") +
  ggtitle("boxplot_significant_random_genes_mean_peaks")

dev.off()

jpeg("~/ifpan-chipseq-timecourse/PLOTS/barplot_significant_random_genes_strongest_peak.jpeg", 
     width = 1400, 
     height = 802)

#making barplot for strongest peak for each gene
merge(x = tmp_significant_random_genes_peak_normalized_amplitude,
      y = {tmp_significant_random_genes_peak_normalized_amplitude %>% 
          group_by(gene.name, start.range, TF, gene.regulation) %>%
          summarise(tmp_mean_alltime_amplitude = mean (amplitude))}) %>%
  group_by(gene.name, TF, gene.regulation) %>%
  filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
  #group_by(gene.name, start.range, time, TF, gene.regulation) %>%
  group_by(time, TF, gene.regulation) %>% 
  summarize(mean = mean(amplitude), sd = sd(amplitude)) %>% #filter(TF=="BCL3", gene.regulation == "random")
  #summarize(mean = mean(mean_amplitude_gene), sd = sd(mean_amplitude_gene)) %>%
  ggplot(., aes(x = as.factor(time), y = mean, fill = gene.regulation)) +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
    facet_wrap(TF ~ ., ncol = 4, scales = "free_y") +
    theme(legend.position = "bottom") +
    ggtitle("barplot_significant_random_genes_strongest_peak")

dev.off()


#Two way ANOVA for the strogest peaks
merge(x = tmp_significant_random_genes_peak_normalized_amplitude,
      y = {tmp_significant_random_genes_peak_normalized_amplitude %>% 
          group_by(gene.name, start.range, TF, gene.regulation) %>%
          summarise(tmp_mean_alltime_amplitude = mean (amplitude))}) %>%
  group_by(gene.name, start.range, time, TF, gene.regulation) %>%
  mutate(mean_three_sample_peak = mean(amplitude)) %>% 
  select(gene.name, start.range, time, TF, gene.regulation, mean_three_sample_peak) %>%
  mutate(log.mean.max.peak = log(mean_three_sample_peak + 1))  -> tmp_maxpeak_data_to_ANOVA
  

lapply(setNames({tmp_significant_random_genes_peak_normalized_amplitude %>% select(TF) %>% .[order(.$TF),] %>% unique()}, 
                {tmp_significant_random_genes_peak_normalized_amplitude %>% select(TF) %>% .[order(.$TF),] %>% unique()}),  
       function(x) {tmp_maxpeak_data_to_ANOVA %>% 
           filter(TF == x)%>% 
           aov(log.mean.max.peak~time*gene.regulation, data = .) %>% 
           summary() %>% 
           .[[1]] %>% 
           .[1:3,5]}) %>% 
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  set_colnames(c("TF", "anova.time", "anova.gene.regulation", "anova.interaction")) %>%
  mutate(fdr.time = p.adjust(.$anova.time, method = "fdr"),
         fdr.gene.regulation = p.adjust(.$anova.gene.regulation, method = "fdr"),
         fdr.interaction = p.adjust(.$anova.interaction, method = "fdr")) %>%
  fwrite("~/ifpan-chipseq-timecourse/DATA/max-peak_TF_ANOVA-time-gene.regulation-interaction_fdr-time-gene.regulation-interaction.tsv", 
         sep="\t", 
         col.names = TRUE, 
         row.names = FALSE)

rm(tmp_maxpeak_data_to_ANOVA)

#Two way ANOVA for all peaks
tmp_significant_random_genes_peak_normalized_amplitude %>%
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.peak = mean(amplitude)) %>%
  select(gene.name, start.range, time, TF, gene.regulation, mean.peak) %>%
  mutate(log.mean.peak = log(mean.peak + 1)) -> tmp_allpeak_data_to_ANOVA

lapply(setNames({tmp_significant_random_genes_peak_normalized_amplitude %>% select(TF) %>% .[order(.$TF),] %>% unique()}, 
                {tmp_significant_random_genes_peak_normalized_amplitude %>% select(TF) %>% .[order(.$TF),] %>% unique()}),  
       function(x) {tmp_allpeak_data_to_ANOVA %>% 
           filter(TF == x)%>% 
           aov(log.mean.peak~time*gene.regulation, data = .) %>% 
           summary() %>% 
           .[[1]] %>% 
           .[1:3,5]}) %>% 
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  set_colnames(c("TF", "anova.time", "anova.gene.regulation", "anova.interaction")) %>%
  mutate(fdr.time = p.adjust(.$anova.time, method = "fdr"),
         fdr.gene.regulation = p.adjust(.$anova.gene.regulation, method = "fdr"),
         fdr.interaction = p.adjust(.$anova.interaction, method = "fdr")) %>%
  fwrite("~/ifpan-chipseq-timecourse/DATA/all-peak_TF_ANOVA-time-gene.regulation-interaction_fdr-time-gene.regulation-interaction.tsv", 
         sep="\t", 
         col.names = TRUE, 
         row.names = FALSE)

rm(tmp_allpeak_data_to_ANOVA)

#calculation mean weighted time
#all peak
tmp_significant_random_genes_peak_normalized_amplitude %>%
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.max.peak = mean(amplitude)) %>% 
  spread(., key = "time", value = "mean.max.peak") %>%
  as.data.frame() %>% mutate(., mean.weighted.time = rowSums(t(t(.[6:16])*{.[6:16] %>% colnames() %>% as.numeric()/60}))/rowSums(.[,6:16])) %>% 
  select(gene.name, start.range, TF, gene.regulation, mean.weighted.time) %>% 
  group_by(gene.regulation, TF) %>% 
  arrange(gene.regulation, TF, desc(mean.weighted.time)) %>% 
  filter(mean.weighted.time <= quantile(mean.weighted.time, 0.7)) %>% 
  #as.data.frame() %>% mutate(., mean.weighted.time = rowSums(t(t(.[5:16])*c(1:12)))/rowSums(.[,5:16])) %>%  
  #as.data.frame() %>% mutate(., mean.weighted.time = rowSums(t(t(.[5:16])*c(1:12)))) %>%  mutate(., tmp = rowSums(.[,5:16]))
  {ggplot(., aes(x = gene.regulation, y = mean.weighted.time, color = gene.regulation)) +
      geom_boxplot() +
      theme(axis.title.x = element_blank()) +
      facet_wrap(TF ~ ., ncol =8) +
      ggtitle("all_peaks")
  }

# #strogest peak
# merge(x = tmp_significant_random_genes_peak_normalized_amplitude,
#       y = {tmp_significant_random_genes_peak_normalized_amplitude %>% 
#           group_by(gene.name, start.range, TF, gene.regulation) %>%
#           summarise(tmp_mean_alltime_amplitude = mean (amplitude))}) %>%
#   group_by(gene.name, TF, gene.regulation) %>%
#   filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
#   group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
#   summarise(mean.max.peak = mean(amplitude))  %>%
#   spread(., key = "time", value = "mean.max.peak") %>%
#   as.data.frame() %>% 
#   mutate(., mean.weighted.time = rowSums(t(t(.[5:16])*c(1:12)))/rowSums(.[,5:16])) %>%  
#   {ggplot(., aes(x = gene.regulation, y = mean.weighted.time, color = gene.regulation)) +
#       geom_boxplot() +
#       theme(axis.title.x = element_blank()) +
#       facet_wrap(TF ~ ., ncol =8) +
#       ggtitle("strongest_peak")
#   }
  
merge(x = tmp_significant_random_genes_peak_normalized_amplitude,
      y = {tmp_significant_random_genes_peak_normalized_amplitude %>% 
          group_by(gene.name, start.range, TF, gene.regulation) %>%
          summarise(tmp_mean_alltime_amplitude = mean (amplitude))}) %>%
  group_by(gene.name, TF, gene.regulation) %>%
  filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.max.peak = mean(amplitude))  %>%
  spread(., key = "time", value = "mean.max.peak") %>%
  as.data.frame() %>%
  mutate(., mean.weighted.time = rowSums(t(t(.[6:16])*{.[6:16] %>% colnames() %>% as.numeric()/60}))/rowSums(.[,6:16])) %>%  
  select(gene.name, start.range, TF, gene.regulation, mean.weighted.time) %>% 
  group_by(gene.regulation, TF) %>% 
  arrange(gene.regulation, TF, desc(mean.weighted.time)) %>% 
  filter(mean.weighted.time <= quantile(mean.weighted.time, 0.7)) %>% str
  #mutate(., mean.weighted.time = rowSums(t(t(.[5:16])*c(1:12)))/rowSums(.[,5:16])) %>%  
  {ggplot(., aes(x = gene.regulation, y = log(mean.weighted.time), color = gene.regulation)) +
      geom_boxplot() +
      theme(axis.title.x = element_blank()) +
      facet_wrap(TF ~ ., ncol =8) +
      ggtitle("strongest_peak")
  }


#checking if normalization is performed on RNA-seq
tmp_data_for_heatmap %>%
  as.data.frame() %>%
  mutate(gene.name = rownames(.)) %>%
  gather(., "time.sample", "value.RNAseq", -c(gene.name)) %>% 
  mutate(number.sample = .$time.sample %>%
           str_split(., "_", n =2, simplify = TRUE) %>%
           .[,2]) %>%
  mutate(time = .$time.sample %>%
           str_split(., "_", n =2, simplify = TRUE) %>%
           .[,1]) %>%
  {ggplot(., aes(x=number.sample, y=value.RNAseq)) +
      geom_boxplot() + 
      facet_wrap(time ~., ncol = 4)}



#fragment do usnięcia, tylko żeby mieć blisko wcześniej jest w kodzie
tmp_data_for_heatmap <- data[order(results$pvalue)[1:number_signification_genes],] %>%
  apply(1, scale) %>% 
  t %>%
  apply(1, function(x, threshold){x[x > threshold] <- threshold; x[x < -threshold] <- -threshold; x}, threshold = 2.0) %>%
  t %>%
  {rownames(.) <- results$gene.name[order(results$pvalue)[1:number_signification_genes]]; .}  %>%
  {colnames(.) <- colnames(data); .}

data[order(results$pvalue)[1:number_signification_genes],] %>%
  as.data.frame() %>%
  mutate(gene.name = rownames(.)) %>%
  gather(., "time.sample", "value.RNAseq", -c(gene.name)) %>% 
  mutate(number.sample = .$time.sample %>%
           str_split(., "_", n =2, simplify = TRUE) %>%
           .[,2]) %>%
  mutate(time = .$time.sample %>%
           str_split(., "_", n =2, simplify = TRUE) %>%
           .[,1]) %>%
           {ggplot(., aes(x=number.sample, y=log2(value.RNAseq + 1))) +
               geom_boxplot() + 
               facet_wrap(time ~., ncol = 4)}

raw.data %>% 
  select(-c(Length)) %>% 
  as.data.frame() %>% 
  gather(., "time.sample", "value.RNAseq", -c(Geneid)) %>% 
  mutate(number.sample = .$time.sample %>%
           str_split(., "_", n =2, simplify = TRUE) %>%
           .[,2]) %>%
  mutate(time = .$time.sample %>%
           str_split(., "_", n =2, simplify = TRUE) %>%
           .[,1]) %>% 
           {ggplot(., aes(x=number.sample, y=log10(value.RNAseq + 1))) +
               geom_boxplot() + 
               facet_wrap(time ~., ncol = 4)}

###quantile.normalization tmp

        
  
  
