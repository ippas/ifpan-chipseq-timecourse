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
install.packages('svglite')
require(svglite)



raw.data <- read.delim("~/ifpan-chipseq-timecourse/DATA/raw_expression_matrix_dexamethasone.tsv",  
                       header = TRUE, 
                       stringsAsFactors = FALSE) %>% 
  set_rownames(.$Geneid)

samples <- read.delim("~/ifpan-chipseq-timecourse/DATA/sample.info.tsv", 
                      header = TRUE, 
                      stringsAsFactors = FALSE)

data <- raw.data[, 3:48] %>% 
#   as.matrix()
# 
# raw.data[, 3:48] %>% 
#   as.matrix() %>% 
#   normalize.quantiles() %>% 
#   set_rownames(., rownames(data)) %>% 
#   set_colnames(colnames(data)) %>%
#   as.data.frame() %>% 
#   rownames_to_column(var="Geneid") %>% 
#   gather(., "time.sample", "value.RNAseq", -c(Geneid)) %>%
#   mutate(log.value.RNAseq = log2(value.RNAseq + 1)) %>%
#   filter(log.value.RNAseq > 1) %>%
#   mutate(number.sample = .$time.sample %>%
#            str_split(., "_", n =2, simplify = TRUE) %>%
#            .[,2]) %>%
#   mutate(time = .$time.sample %>%
#            str_split(., "_", n =2, simplify = TRUE) %>%
#            .[,1]) %>%
#            {ggplot(., aes(x=number.sample, y=log.value.RNAseq)) +
#                geom_boxplot() + 
#                facet_wrap(time ~., ncol = 4)}

#quantile normalization
raw.data[, 3:48] %>% 
  as.matrix() %>% 
  normalize.quantiles() %>%
  set_rownames(., rownames(data)) %>% 
  set_colnames(colnames(data)) -> data


genes.names <- read.delim("~/ifpan-chipseq-timecourse/DATA/ID_ID.version_gene.tsv", 
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
  left_join(., genes.names, by = "Gene.stable.ID") %>%
  rename(ensemblid = Gene.stable.ID,
         gene.name = Gene.name) %>%
  select(-Gene.stable.ID.version) %>%
  left_join(., tmp_pvalue_dataframe, by = "Geneid") %>%
  select("Geneid", "ensemblid", "Length", "gene.name", "pvalue")

rm(tmp_pvalue_dataframe)


#choose number genes
number_signification_genes <- results %>% 
  mutate(fdr = p.adjust(results$pvalue, method = "fdr")) %>% 
  filter(fdr < 0.0000001) %>% 
  nrow()

#test before corelation data to heat map
results %>% 
  mutate(fdr = p.adjust(results$pvalue, method = "fdr")) %>% 
  filter(fdr < 0.0000001) -> results.filtered

##########################################################
# selecting significant genes and dividing into clusters #
##########################################################

to.plot <- data[order(results$pvalue)[1:number_signification_genes],] %>%
  as.matrix() %>% 
  #normalize.quantiles() %>% 
  set_rownames(., rownames(raw.data[match(results.filtered$Geneid, rownames(raw.data)),3:48])) %>% 
  set_colnames(colnames(raw.data[match(results.filtered$Geneid, rownames(raw.data)),3:48])) %>%
  
  apply(1, scale) %>% 
  t %>%
  apply(1, function(x, threshold){x[x > threshold] <- threshold; x[x < -threshold] <- -threshold; x}, threshold = 2.0) %>%
  t %>%
  {rownames(.) <- results$gene.name[order(results$pvalue)[1:number_signification_genes]]; .}  %>%
  {colnames(.) <- colnames(data); .}

number_clusters <- 2

tmp_dend <- as.dist(1-cor(t(to.plot), method = "spearman")) %>% 
  hclust %>% 
  as.dendrogram %>% 
  color_branches(., k = number_clusters, col=c("firebrick", "dodgerblue"))

#Create heatmap changes in expression significant genes
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/heatmap_expression_significant.svg", 
        width = 10,
        height = 8)

heatmap.2(to.plot,
          trace="none",
          margins =c(5,7),
          Colv = FALSE,
          col=bluered(20),
          dendrogram="row",      
          Rowv = tmp_dend,  
          key.xlab = "Concentration (index)",
          distfun = function(x) as.dist(1-cor(t(x)), method = "spearman"),
          RowSideColors =  (get_leaves_branches_col(tmp_dend) %>% 
                              .[order(order.dendrogram(tmp_dend))]), 
          colRow = (get_leaves_branches_col(tmp_dend) %>% 
                      .[order(order.dendrogram(tmp_dend))])
) 

dev.off()

gene_regulation <- to.plot %>% 
  t %>% 
  {1 - cor(.)} %>% 
  as.dist %>% 
  hclust %>% 
  as.dendrogram %>%
  cutree(., k = number_clusters) %>% 
  as.table() %>% 
  as.data.frame() %>%  
  set_colnames(c("gene.name", "number_regulation")) %>% mutate(gene.regulation=number_regulation) %>% 
  mutate(gene.regulation = recode(number_regulation,
                                  `1` = "up-regulated",
                                  `2` = "down-regulated")) %>%
  select(gene.name, gene.regulation)

##########################################################

#############################
# selection of random genes #
#############################

#generate randomly expressed genes with the same distribution of mean expression:
#but first lets prepare random genes the same way

#part 1: see the distribution of expression and generate histogram:

raw.data[match(results.filtered$Geneid, rownames(raw.data)),3:48] %>%
  rowMeans() -> results.filtered$mean.expression

raw.data[match(results$Geneid, rownames(raw.data)),3:48] %>%
  rowMeans() -> results$mean.expression


hist(log2(as.numeric(results.filtered$mean.expression)+1), breaks=10) -> expression.pattern
sum(expression.pattern$counts)

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
  rm(temp.random)
}

#################################################################
###########################################
# reconstruction random genes # to remove #
###########################################
read.delim("~/ifpan-chipseq-timecourse/DATA/enhancer_info.tsv", 
           header = TRUE, 
           stringsAsFactors = FALSE) %>% 
  filter(gene.regulation == "random") %>%
  select(gene.name) %>%
  left_join(., results, by = "gene.name") %>%
  select(-mean.expression) -> random
#################################################################



hist(log2(as.numeric(random$mean.expression)+1), breaks=10)

dev.off()

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

##########################################################


###############################################################################
# Create line plot relative changes in expression significant and random genes#
###############################################################################
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/lineplot_change_expression.svg", 
        width = 10,
        height = 8)

to.plot %>%  
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
  arrange(time) %>% 
  na.omit() %>% 
  ungroup() %>%
  mutate(gene.regulation=as.character(gene.regulation)) %>% 
  rbind(., random.prepared) %>% 
  group_by(gene.name, gene.regulation, time) %>% 
  arrange(gene.name, time) %>% 
  group_by(gene.name, gene.regulation) %>%
  mutate(mean = mean - mean[1]) %>% 
  ungroup -> tmp_expression_normalize

tmp_expression_normalize %>%
  {ggplot(., 
          aes(x = as.numeric(time), y = mean, color = as.factor(gene.regulation))) +
      geom_line(aes(group = gene.name), alpha = 0.01) +
      geom_smooth(aes(group = gene.regulation), se = FALSE, size = 1) +
      theme(
            legend.position = "bottom") +
      scale_color_manual(values = c("random" = "gray",
                                    "up-regulated" = "firebrick", 
                                    "down-regulated" = "dodgerblue")) +
      theme(legend.position = "bottom") +
      scale_x_continuous(limits=c(0, 720), breaks = c(0, 120, 240, 360, 480, 600, 720)) +
      labs(x = "Time [min]",
           y = "Relative changes in gene expression")}

dev.off()

tmp_expression_normalize %>%
  fwrite("~/ifpan-chipseq-timecourse/DATA/expression_normalize_significant_random.tsv", 
         sep="\t", 
         col.names = TRUE, 
         row.names = FALSE)


# remove tmp variable
rm(to.plot,
   tmp_dend,
   tmp_expression_normalize)

#########################
# Max change time point #
# #########################
# column_differences <- function(z){
#   m <- c(z[1], z[2])
#   for(i in c(4:14)){
#     y <- abs(as.numeric(z[i]) - as.numeric(z[i-1]))
#     m <-c(m, y)
#   }
#   return(m)
# }
# 
# 
# data[order(results$pvalue)[1:number_signification_genes],] %>%
#   as.matrix() %>% 
#   set_rownames(., rownames(raw.data[match(results.filtered$Geneid, rownames(raw.data)),3:48])) %>% 
#   set_colnames(colnames(raw.data[match(results.filtered$Geneid, rownames(raw.data)),3:48])) %>% 
#   {rownames(.) <- results$gene.name[order(results$pvalue)[1:number_signification_genes]]; .}  %>%
#   as.data.frame() %>% 
#   rownames_to_column(., var = "gene.name") %>% 
#   gather(., key = "samplied", value = "value_expression", -gene.name) %>%
#   left_join(., {samples %>% select(samplied, time)}) %>% 
#   select(-samplied) %>% 
#   group_by(gene.name, time) %>%
#   summarise(mean_value_expression = as.numeric(mean(value_expression))) %>% 
#   left_join(., gene_regulation, by = "gene.name") %>% 
#   spread(., key = "time", value = "mean_value_expression") %>% 
#   as.data.frame() %>% 
#   apply(., 1, column_differences) %>% 
#   t %>% 
#   set_colnames(c("gene.name", "gene.regulation", "15","45", "90", "150", "210", "270", "330", "390", "450", "540", "660")) %>% 
#   as.data.frame() %>%
#   gather(., key = "time", value = "value", -c(gene.name, gene.regulation)) %>%
#   mutate(value = as.numeric(value)) %>% 
#   spread(., key = "time", value = "value") %>% 
#   mutate(sum_row = apply(., 1, function(x){sum(as.numeric(x[3:13]))})) %>% 
#   gather(., key = "time", value = "value", -c(gene.name, gene.regulation, sum_row)) %>% 
#   mutate(time = as.numeric(time)) %>% 
#   mutate(time_value = time *  value) %>%
#   group_by(gene.name, gene.regulation, sum_row) %>%
#   summarise(sum_time_value = sum(time_value)) %>%
#   mutate(max_change_time_point = sum_time_value/sum_row) -> max.change.time.point.significant

# jpeg("~/ifpan-chipseq-timecourse/PLOTS/boxplot_MCTP.jpeg", 
#      width = 1400, 
#      height = 802)
# svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/boxplot_MCTP.jpeg", 
#         width = 10,
#         height = 8)
# 
# max.change.time.point.significant %>%
# {ggplot(., aes(x = gene.regulation, y = max_change_time_point)) +
#     geom_boxplot() +
#     labs(x = "Gene regulation",
#          y = "Time [min]") +
#     ggtitle("Max change time point")}
# 
# dev.off()


###########################################
# prepare data to calculate RPKM and FPKM #
###########################################
gene_chromosome_start_end_strand <- read.delim("~/ifpan-chipseq-timecourse/DATA/gene_chromosome_start_end_strand.tsv", 
                                               header = TRUE, 
                                               stringsAsFactors = FALSE)


tmp_transcript_length <- as.data.frame(read.delim("~/ifpan-chipseq-timecourse/DATA/transcript_length.tsv", 
                                                  header = TRUE, 
                                                  sep = "\t", 
                                                  stringsAsFactors = FALSE)) %>% 
  group_by(.$Gene.stable.ID) %>% 
  summarize(median = median(Transcript.length..including.UTRs.and.CDS.)) %>% 
  set_colnames(c("ensemblid", "median"))


svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/histogram_RPKM.svg", 
        width = 10,
        height = 8)

# Historgram RPKM (some number ensemblid do not have match)
data[order(results$pvalue)[1:number_signification_genes],] %>% 
  rbind(., data[random$Geneid,]) %>%
  as.data.frame() %>%
  mutate(ensemblid = {rownames(.) %>% str_split(., fixed("."), n=2, simplify = TRUE) %>% .[,1]}) %>%
  left_join(., tmp_transcript_length, by = "ensemblid") %>% 
  na.omit()  %>% 
  {.[,1:46] <- ((.[,1:46]/rep(colSums(.[,1:46]), each = nrow(.)))*10^3 * 10^6)/.$median; .} %>% 
  mutate(., mean = rowSums(.[,1:46])/46)  %>% 
  left_join(.,results[, c("ensemblid", "gene.name")], by = "ensemblid" ) %>%
  select(ensemblid, median, mean, gene.name) %>%
  left_join(., rbind(gene_regulation, random  %>% select(gene.name) %>% mutate(gene.regulation = "random")), 
            by = "gene.name") %>% 
  mutate(log2mean = log2(mean + 1)) %>% 
  {ggplot(.,aes(x = log2mean)) + 
      geom_density(aes(y=..density.., color = gene.regulation)) +
      geom_histogram(aes(y=..density.., fill = gene.regulation), position="identity", alpha=0.4, bins = 30) +
      scale_color_manual(values = c("up-regulated" = "firebrick",
                                    "random" = "gray20",
                                    "down-regulated" = "dodgerblue")) +
      scale_fill_manual(values = c("up-regulated" = "firebrick",
                                   "random" = "gray",
                                   "down-regulated" = "dodgerblue")) +
      ggtitle("RPKM")
  } 

dev.off()


svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/histogram_transcript_lenght.svg", 
        width = 10,
        height = 8)

#historgram for FPKM
data[order(results$pvalue)[1:number_signification_genes],] %>% 
  rbind(., data[random$Geneid,]) %>% 
  as.data.frame() %>%
  mutate(ensemblid = {rownames(.) %>% str_split(., fixed("."), n=2, simplify = TRUE) %>% .[,1]}) %>%
  left_join(., tmp_transcript_length, by = "ensemblid") %>% 
  na.omit()  %>% 
  left_join(.,results[, c("ensemblid", "gene.name")], by = "ensemblid" ) %>%
  select(ensemblid, median, gene.name) %>% 
  left_join(., rbind(gene_regulation, random  %>% select(gene.name) %>% mutate(gene.regulation = "random")), 
            by = "gene.name") %>%  
  mutate(log.median = log(median)) %>%
  filter(median < 90000) %>%
  {ggplot(.,aes(x = median)) + 
      geom_histogram(aes(y=..density.., fill = gene.regulation), position="identity", alpha=0.4) +
      scale_color_manual(values = c("up-regulated" = "firebrick",
                                    "random" = "gray20",
                                    "down-regulated" = "dodgerblue")) +
      scale_fill_manual(values = c("up-regulated" = "firebrick",
                                   "random" = "gray",
                                   "down-regulated" = "dodgerblue")) +
      
      ggtitle("Transcript length")
  } 

dev.off()
###########################################################

###########################################################
# save to file information to extract range for promoters #
###########################################################

random %>%
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
  fwrite("~/ifpan-chipseq-timecourse/DATA/promotores_peaks_info.tsv", 
         sep="\t", 
         col.names = TRUE, 
         row.names = FALSE)

###########################################################
# save to file information to finding enhancers for genes #
###########################################################

random %>% 
  select(ensemblid, gene.name) %>%
  left_join(., gene_chromosome_start_end_strand, by = c("ensemblid" = "Gene.stable.ID")) %>%
  mutate(pos=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% 
  mutate(start = pos - 100000, end = pos + 100001) %>%
  select("ensemblid", "gene.name", "Chromosome.scaffold.name", "start", "end") %>%
  set_colnames(c("ensemblid", "gene.name", "chromosome", "start", "end")) %>% 
  mutate(gene.regulation = "random") %>% na.omit() %>% 
  rbind(., {results[order(results$pvalue)[1:number_signification_genes],] %>% 
      left_join(., gene_chromosome_start_end_strand, by = c("ensemblid" = "Gene.stable.ID")) %>% 
      mutate(pos=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% 
      mutate(start = pos - 100000, end = pos + 100001) %>% 
      select("ensemblid", "gene.name", "Chromosome.scaffold.name", "start", "end") %>%
      set_colnames(c("ensemblid", "gene.name", "chromosome", "start", "end")) %>% 
      left_join(., gene_regulation, by = "gene.name") %>% 
      na.omit()}) %>% 
  mutate(start = if_else(start < 0, 0, start)) %>% 
  fwrite("~/ifpan-chipseq-timecourse/DATA/enhancer_info.tsv", 
         sep="\t", 
         col.names = TRUE, 
         row.names = FALSE)