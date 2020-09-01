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


raw.data <- read.delim("~/ifpan-chipseq-timecourse/DATA/raw_expression_matrix_dexamethasone.tsv",  
                       header = TRUE, 
                       stringsAsFactors = FALSE) %>% 
  set_rownames(.$Genedit)

samples <- read.delim("~/ifpan-chipseq-timecourse/DATA/sample.info.tsv", 
                      header = TRUE, 
                      stringsAsFactors = FALSE)

data <- raw.data[, 3:48] %>% 
  as.matrix()

ID_ID.version_gene <- read.delim("~/ifpan-chipseq-timecourse/DATA/ID_ID.version_gene.txt", 
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
#results %>% mutate(fdr = p.adjust(results$pvalue, method = "fdr"))
#sum(results$fdr < 0.0016, na.rm = T)
number_signification_genes <- 737

tmp_data_for_heatmap <- data[order(results$pvalue)[1:number_signification_genes],] %>%
  apply(1, scale) %>% 
  t %>%
  apply(1, function(x, threshold){x[x > threshold] <- threshold; x[x < -threshold] <- -threshold; x}, threshold = 2.0) %>%
  t %>%
  {rownames(.) <- results$genename[order(results$pvalue)[1:number_signification_genes]]; .}  %>%
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
          Rowv = tmp_dend1,  
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
jpeg("~/ifpan-chipseq-timecourse/PLOTS/lineplot_up_down_regulation_significant_genes.jpeg", 
     width = 1400, 
     height = 802)

tmp_data_for_heatmap %>% 
  melt() %>% 
  na.omit() %>% 
  as.data.frame() %>% #head %>%
  mutate(Var2 = str_remove_all(Var2, pattern = "_rep" %R% one_or_more(DGT))) %>%
  mutate(Var2 = str_remove(Var2, pattern = START %R% "t")) %>% 
  mutate(Var2 = 60 * as.numeric(recode(Var2, 
                       "00" = "0", 
                       "05" = "0.5")))  %>%
  set_colnames(c("gene.name", "time", "value")) %>%
  group_by(gene.name, time) %>% 
  summarize(mean = mean(value)) %>%
  arrange(time) %>%
  left_join(., gene_regulation, by = "gene.name") %>%
  ggplot(., 
       aes(x = as.factor(time), y = mean, color = gene.regulation)) +
  geom_line(aes(group = gene.name), alpha = 0.05) +
  geom_smooth(aes(group = gene.regulation), se = FALSE, size = 2) +
  theme(legend.position = "bottom")

#Create Barplot for a single gene
data.frame(exprs = data[which(results$genename == "NYAP1"),], 
           time = factor(samples$time, 
                         levels = c("0", "30", "60", "120", "180", "240", "300", "360", "420", "480", "600", "720"))) %>% 
  group_by(time) %>% 
  summarise(mean = mean(exprs), sd = sd(exprs)) %>%
  ggplot(., aes(x=time, y=mean)) + 
      geom_bar(position=position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                    width=.2,
                    position=position_dodge(.9))

dev.off()

# remove tmp variable
rm(tmp_data_for_heatmap,
   tmp_dend,
   tmp_col_labels)

gene_chromosome_start_end_strand <- read.delim("~/ifpan-chipseq-timecourse/DATA/gene_chromosome_start_end_strand.tsv", 
                                               header = TRUE, 
                                               stringsAsFactors = FALSE)

results[order(results$pvalue)[1:number_signification_genes],] %>%
  left_join(., gene_chromosome_start_end_strand, by = c("ensemblid" = "Gene.stable.ID")) %>% 
  mutate(pos=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% 
  mutate(start = pos - 10000, end = pos + 10001) %>% 
  select("ensemblid", "gene.name", "Chromosome.scaffold.name", "start", "end") %>%
  set_colnames(c("ensemblid", "gene.name", "chromosome", "start", "end")) %>% 
  na.omit() %>% 
  fwrite("~/ifpan-chipseq-timecourse/DATA/significant_genes_ensemblid_genename_chromosome_start_end.tsv", 
         sep="\t", 
         col.names = TRUE, 
         row.names = FALSE)
  

tmp_transcript_length <- as.data.frame(read.delim("~/ifpan-chipseq-timecourse/DATA/transcript_length.tsv", 
                                                  header = TRUE, 
                                                  sep = "\t", 
                                                  stringsAsFactors = FALSE)) %>% 
  group_by(.$Gene.stable.ID) %>% 
  summarize(median = median(Transcript.length..including.UTRs.and.CDS.)) %>% 
  set_colnames(c("ensemblid", "median"))


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
  mutate(logmean = log2(mean))  %>% 
  {ggplot(.,aes(x = logmean)) + 
      geom_histogram()} 


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
  column_to_rownames(., var = "ensemblid") %>% 
  .[{.$mean > 2 & .$mean < 8192},] %>%
  .[sample(nrow(.), 1000), 47] %>%
  as.data.frame() %>% 
  set_colnames("Gene.name") %>% 
  left_join(., gene_chromosome_start_end_strand, by = "Gene.name") %>% 
  mutate(pos=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% 
  mutate(start = pos - 10000, end = pos + 10001) %>%  
  select("Gene.stable.ID", "Gene.name", "Chromosome.scaffold.name", "start", "end") %>% 
  set_colnames(c("ensemblid", "gene.name", "chromosome", "start", "end")) %>% head

  


