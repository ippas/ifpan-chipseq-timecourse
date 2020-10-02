require(dplyr)
require(magrittr)
require(reshape2)
require(tidyr)
require(gplots)
require(ggplot2)
require(dendextend)
require(data.table)
require(matrixStats)
require(rebus)
require(tibble)
require(stringr)

raw.data <- read.delim(#"~/ifpan-chipseq-timecourse/DATA/raw_expression_matrix_dexamethasone.tsv", 
                       "DATA/raw_expression_matrix_dexamethasone.tsv",
                       header = TRUE, 
                       stringsAsFactors = FALSE) %>% 
  set_rownames(.$Geneid)

samples <- read.delim(#"~/ifpan-chipseq-timecourse/DATA/sample.info.tsv",
                      "DATA/sample.info.tsv", 
                      header = TRUE, 
                      stringsAsFactors = FALSE)

data <- raw.data[, 3:48] %>% 
  as.matrix()

ID_ID.version_gene <- read.delim(#"~/ifpan-chipseq-timecourse/DATA/ID_ID.version_gene.tsv",
                                 "DATA/ID_ID.version_gene.tsv", 
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

#set FDR treshold
results %>% 
  mutate(fdr = p.adjust(results$pvalue, method = "fdr")) %>% 
  filter(fdr < 0.001) -> results.filtered


#add 640 random genes to plotting and further analyses


#part 1: see the distribution of expression and generate histogram:

raw.data[match(results.filtered$Geneid, rownames(raw.data)),3:48] %>%
  rowMeans() -> results.filtered$mean.expression

raw.data[match(results$Geneid, rownames(raw.data)),3:48] %>%
  rowMeans() -> results$mean.expression


hist(log2(as.numeric(results.filtered$mean.expression)+1), breaks=10) -> expression.pattern
sum(expression.pattern$counts)

#generate randomly expressed genes with the same distribution of mean expression:

 
random <- data.frame()

for (i in seq_along(expression.pattern$mids)) {
  print(i)
  start <- expression.pattern$mids[i] - 1
  stop <- expression.pattern$mids[i] + 1 
  results %>% 
    mutate(fdr = p.adjust(results$pvalue, method = "fdr")) %>%
    filter(fdr > 0.1) %>% 
    filter((log2(as.numeric(mean.expression)+1) > start) &
             (log2(as.numeric(mean.expression)+1) <= stop)) %>% 
    sample_n(expression.pattern$counts[i]) -> temp.random
  print(nrow(temp.random))
  random <- rbind(random, temp.random)
}
  
write.table(random, 
          # "~/ifpan-chipseq-timecourse/DATA/random_genes_geneid_ensemblid_length_gene.name_pvalue_mean.expression.tsv", 
             "DATA/random_genes_geneid_ensemblid_length_gene.name_pvalue_mean.expression.tsv", 
             row.names = FALSE)

hist(log2(as.numeric(random$mean.expression)+1), breaks=10)

# at this point we have : results.filtered - this has our chosen genes and random - with 640 random genes


to.plot <- data[order(results$pvalue)[1:640],] %>%
  apply(1, scale) %>% 
  t %>%
  apply(1, function(x, threshold){x[x > threshold] <- threshold; x[x < -threshold] <- -threshold; x}, threshold = 2.0) %>%
  t %>%
  {rownames(.) <- results$gene.name[order(results$pvalue)[1:640]]; .}  %>%
  {colnames(.) <- colnames(data); .}

number_clusters <- 2

tmp_dend <- as.dist(1-cor(t(to.plot))) %>% 
  hclust %>% 
  as.dendrogram %>% 
  color_branches(., k = number_clusters)

tmp_col_labels <- get_leaves_branches_col(tmp_dend) %>% 
  .[order(order.dendrogram(tmp_dend))]


#Create heatmap

heatmap.2(to.plot,
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



gene_regulation <- to.plot %>% 
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





to.plot %>% 
  melt() %>% 
  na.omit() %>% 
  as.data.frame() %>% #head %>%
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
  geom_line(data = random.prepared, aes(group = gene.name), alpha = 0.02, color="grey") +
  geom_smooth(data = random.prepared, aes(group = gene.regulation), se = FALSE, size = 1, color="grey") +
  geom_line(aes(group = gene.name), alpha = 0.02) +
  geom_smooth(aes(group = gene.regulation), se = FALSE, size = 1) +
  theme(legend.position = "bottom")

