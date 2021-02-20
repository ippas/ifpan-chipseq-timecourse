####################################poprawić kod, usunąć promotory po wczytaniu pliku
#  choose peaks without promoters  #
####################################
# peaks.enhancer.all.genes <- read.table('~/ifpan-chipseq-timecourse/DATA/peaks_all_genes.tsv', 
#                                        header = TRUE, 
#                                        sep = '\t')
# 
# peaks.enhancer.all.genes %>% 
#   left_join(., gene_chromosome_start_end_strand, by = c("gene.name" = "Gene.name")) %>% 
#   select(-c(Chromosome.scaffold.name, Gene.stable.ID)) %>% 
#   mutate(pos.TSS=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% 
#   mutate(start.range.promoter = pos.TSS - 2000, end.range.promoter = pos.TSS + 2001) %>% 
#   filter(end_peak < start.range.promoter | start_peak > end.range.promoter) %>% 
#   .[, 1:5] %>%
#   mutate(gene.regulation = "NA") %>%
#   fwrite("~/ifpan-chipseq-timecourse/DATA/peaks_all_genes_without_promoters.tsv",
#          sep="\t", 
#          col.names = TRUE, 
#          row.names = FALSE)
 

# peaks_all_genes_without_promoters_amplitude <- read.table("~/ChIP-seq/DATA/peaks_all_genes_without_promoters_amplitude.tsv", #zmiana pliku do wczytywania danych
#                                                           header = FALSE, 
#                                                           sep = "\t", 
#                                                           stringsAsFactors = FALSE) %>% 
#   set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", "amplitude"))

peaks_all_genes_ep300_nr3c1_amplitude <- read.table('~/ifpan-chipseq-timecourse/DATA/peaks_all_genes_ep300_nr3c1_amplitude.tsv',
                                                    header = TRUE, 
                                                    sep = "\t") %>%
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", "amplitude")) %>% 
  remove_peak_promoter() 

#####################################################################
#  choose and plot for too 100 peaks ep300 max differences min-max  #
#####################################################################

peaks_all_genes_ep300_nr3c1_amplitude %>% 
  select(-c(file)) %>% 
  group_by(gene.name, start.range, end.range, TF, time, gene.regulation) %>% 
  summarise(amplitude_mean = mean(amplitude)) %>% 
  spread(., key = "time", value = "amplitude_mean") %>%
  as.data.frame() %>% 
  mutate(max_value = apply(., 1, function(x){as.numeric(x[7:17]) %>% max()})) %>%
  mutate(min_value = apply(., 1, function(x){as.numeric(x[7:17]) %>% min()})) %>% 
  mutate(differ_min_max = max_value - min_value) %>% 
  filter(TF == "EP300") %>% 
  select(-gene.name) %>% 
  .[order(.$differ_min_max, decreasing = TRUE),] %>% #sort by differ_min_max
  unique() %>% 
  #head(150) %>% # choose top 150 peaks
  left_join(., {peaks_all_genes_without_promoters_amplitude %>% 
      select(-c(file)) %>% 
      group_by(gene.name, start.range, end.range, TF, time, gene.regulation) %>% 
      summarise(amplitude_mean = mean(amplitude)) %>% 
      spread(., key = "time", value = "amplitude_mean") %>%
      filter(TF == "EP300") %>% .[,1:3]}, by = c("start.range", "end.range")) %>%
  select(gene.name, everything()) %>% 
  left_join(., {gene_chromosome_start_end_strand %>% filter(Gene.stable.ID  %in% {results %>% .[,2]})}, 
            by = c("gene.name" = "Gene.name")) %>% 
  mutate(pos.TSS=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% #marking position TSS for genes
  select(-c("Chromosome.scaffold.name", "Gene.start..bp.", "Gene.end..bp.", "Strand")) %>% # remove unnecessary columns
  mutate(diff.TSS.start.peak = abs(pos.TSS-start.range)) %>% # calculate disctance between start enhancer and position TSS
  group_by(start.range) %>% 
  filter(diff.TSS.start.peak == min(diff.TSS.start.peak)) %>% # choose gene which are the closest TSS
  as.data.frame() %>% 
  ungroup() %>% 
  group_by(gene.name) %>% 
  filter(diff.TSS.start.peak == min(diff.TSS.start.peak)) %>% # choose enhancer closer to the TSS
  as.data.frame() %>% 
  head(100) %>% # choose top 100 peaks
  .[, 1:17] %>% 
  gather(., key = "time", value = "amplitude", 
         -c(gene.name, 
            start.range, 
            end.range, TF, 
            gene.regulation)) -> delta_ep300

###########################
# prepare data to heatmap #
###########################
data %>%
  as.data.frame() %>%
  rownames_to_column(., var = "Geneid") %>%
  left_join(., results[,c(1,4)], by = "Geneid") %>%
  filter(gene.name %in% {delta_ep300 %>% select(gene.name) %>% unique %>% .[, 1]}) %>%
  column_to_rownames(., var = "gene.name") %>%
  select(-Geneid) %>%
  #filter(gene.name == "ALDH1A1")
  as.matrix() %>%
  apply(1, scale) %>%
  t %>%
  apply(1, function(x, threshold){x[x > threshold] <- threshold; x[x < -threshold] <- -threshold; x}, threshold = 2.0) %>%
  t %>%
  {colnames(.) <- colnames(data); .} %>%
  as.data.frame() %>% na.omit() %>%
  as.matrix() -> to.heatmap
  
# data[order(results$pvalue),] %>%
#   as.matrix() %>% 
#   #set_rownames(., rownames(raw.data[match(results$Geneid, rownames(raw.data)),3:48])) %>% head
#   #set_colnames(colnames(raw.data[match(results$Geneid, rownames(raw.data)),3:48])) %>%
#   apply(1, scale) %>% 
#   t %>%
#   apply(1, function(x, threshold){x[x > threshold] <- threshold; x[x < -threshold] <- -threshold; x}, threshold = 2.0) %>%
#   t %>%
#   {rownames(.) <- results$Geneid[order(results$pvalue)]; .}  %>%
#   {colnames(.) <- colnames(data); .} %>%
#   as.data.frame() %>% 
#   rownames_to_column() %>%
#   filter(rowname %in% {delta_ep300 %>%
#       select(gene.name) %>% unique %>%
#       .[, 1]}) %>%
#   column_to_rownames(., var = "rowname") %>%
#   as.matrix() %>%
#   na.omit -> to.heatmap

number_clusters <- 2

# heatmap delta ep300 #
tmp_dend <- as.dist(1-cor(t(to.heatmap), method = "spearman")) %>% 
  hclust %>% 
  as.dendrogram %>% 
  color_branches(., k = number_clusters, col = c("skyblue4", "lightpink4"))


svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/heatmap_expression_delta_ep300.svg", 
        width = 10,
        height = 8)

heatmap.2(to.heatmap,
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

rm(tmp_dend)

##########################
# divide on two clusters #
##########################
gene_regulation_ep300 <- to.heatmap %>% 
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
                                  `1` = "up-delta_ep300",
                                  `2` = "down-delta_ep300")) %>%
  select(gene.name, gene.regulation) 

###############################
# line plot change expression #
###############################
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/lineplot_change_expression_delta_ep300.svg", 
        width = 10,
        height = 8)

to.heatmap %>%  
  melt() %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  mutate(Var2 = str_remove_all(Var2, pattern = "_rep" %R% one_or_more(DGT))) %>%
  mutate(Var2 = str_remove(Var2, pattern = START %R% "t")) %>% 
  mutate(Var2 = 60 * as.numeric(recode(Var2, 
                                       "00" = "0", 
                                       "05" = "0.5")))  %>%
  set_colnames(c("gene.name", "time", "value")) %>%
  left_join(., gene_regulation_ep300, by = "gene.name") %>% 
  group_by(gene.name, time, gene.regulation) %>% 
  summarize(mean = mean(value)) %>% 
  arrange(time) %>% 
  na.omit() %>% 
  ungroup() %>%
  mutate(gene.regulation=as.character(gene.regulation)) %>% 
  group_by(gene.name, gene.regulation, time) %>% 
  arrange(gene.name, time) %>% 
  group_by(gene.name, gene.regulation) %>%
  mutate(mean = mean - mean[1]) %>% 
  ungroup %>%
  {ggplot(., 
          aes(x = as.numeric(time), y = mean, color = as.factor(gene.regulation))) +
      geom_line(aes(group = gene.name), alpha = 0.1) +
      stat_smooth(aes(group = gene.regulation), se = FALSE, size = 1) +
      scale_color_manual(values = c("up-delta_ep300" = "lightpink4", 
                                    "down-delta_ep300" = "skyblue4")) +
      theme(legend.position = "bottom") +
      scale_x_continuous(limits=c(0, 720), breaks = c(0, 120, 240, 360, 480, 600, 720)) +
      labs(x = "Time [min]",
           y = "Relative changes in gene expression") +
      ggtitle("delta ep300")}

dev.off()



#####################################################################
# Boxplot, which show change amplitude for enhancer during the time #
#####################################################################
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/boxplot_amplitude_delta_ep300.svg", 
        width = 10,
        height = 8)


peaks_all_genes_without_promoters_amplitude %>%
  select(-c(file)) %>% 
  group_by(gene.name, start.range, end.range, TF, time, gene.regulation) %>% 
  summarise(amplitude_mean = mean(amplitude)) %>%
  rename(amplitude = amplitude_mean) %>%
  ungroup %>%
  mutate(id = paste(gene.name, start.range, end.range, sep = '.')) %>% 
  filter(id %in% {delta_ep300 %>% #filter by delta_ep300
      mutate(id = paste(gene.name, start.range, end.range, sep = '.')) %>%
      .[, 8]}) %>% 
  select(-id) %>% 
  group_by(gene.name, start.range, time, TF) %>% 
  summarise(mean.max.peak = mean(amplitude)) %>% 
  left_join(., gene_regulation_ep300, by = "gene.name") %>%
  na.omit %>%
  {ggplot(., aes(x = as.factor(time), y = log(mean.max.peak), color = gene.regulation)) +
      geom_boxplot(position = position_dodge(), outlier.size = 0) +  
      facet_wrap(TF ~ ., ncol = 2, scales = "free_y" ) +
      theme(legend.position = "bottom") +
      scale_color_manual(values = c("up-delta_ep300" = "lightpink4", 
                                    "down-delta_ep300" = "skyblue4")) +
      labs(x = "Time [min]",
           y = "Log amplitude values") +
      ggtitle("Delta ep300")}

dev.off()


###########################################
# Preapare data with MCTP and MWT to plot #
###########################################
data[order(results$pvalue),] %>%
  as.matrix() %>% 
  set_rownames(., rownames(raw.data[match(results$Geneid, rownames(raw.data)),3:48])) %>% 
  set_colnames(colnames(raw.data[match(results$Geneid, rownames(raw.data)),3:48])) %>% 
  {rownames(.) <- results$gene.name[order(results$pvalue)]; .} %>%
  as.data.frame() %>%
  rownames_to_column(., var = "gene.name") %>% 
  filter(gene.name %in% {delta_ep300 %>% .[,1]}) %>%
  gather(., key = "samplied", value = "value_expression", -gene.name) %>%
  left_join(., {samples %>% select(samplied, time)}) %>% 
  select(-samplied) %>% 
  group_by(gene.name, time) %>%
  summarise(mean_value_expression = as.numeric(mean(value_expression))) %>% 
  left_join(., gene_regulation_ep300, by = "gene.name") %>% 
  spread(., key = "time", value = "mean_value_expression") %>% 
  as.data.frame() %>%
  apply(., 1, column_differences) %>% 
  t %>% 
  set_colnames(c("gene.name", "gene.regulation", "15","45", "90", "150", "210", "270", "330", "390", "450", "540", "660")) %>% 
  as.data.frame() %>%
  gather(., key = "time", value = "value", -c(gene.name, gene.regulation)) %>%
  mutate(value = as.numeric(value)) %>% 
  spread(., key = "time", value = "value") %>% 
  mutate(sum_row = apply(., 1, function(x){sum(as.numeric(x[3:13]))})) %>% 
  gather(., key = "time", value = "value", -c(gene.name, gene.regulation, sum_row)) %>% 
  mutate(time = as.numeric(time)) %>% 
  mutate(time_value = time *  value) %>% 
  group_by(gene.name, gene.regulation, sum_row) %>%
  summarise(sum_time_value = sum(time_value)) %>%
  mutate(max_change_time_point = sum_time_value/sum_row) %>% 
  na.omit %>% 
  as.data.frame() %>%
  filter(gene.regulation == "up-delta_ep300") %>% 
  mutate(type.data = "expression") %>%
  rename(value = max_change_time_point) %>%
  select(gene.name, type.data ,gene.regulation, value) -> MCTP_delta_ep300

 
peaks_all_genes_without_promoters_amplitude %>%
  select(-c(file)) %>% 
  group_by(gene.name, start.range, end.range, TF, time, gene.regulation) %>% 
  summarise(amplitude_mean = mean(amplitude)) %>%
  rename(amplitude = amplitude_mean) %>%
  ungroup %>%
  mutate(id = paste(gene.name, start.range, end.range, sep = '.')) %>% 
  filter(id %in% {delta_ep300 %>% #filter by delta_ep300
      mutate(id = paste(gene.name, start.range, end.range, sep = '.')) %>%
      .[, 8]}) %>% 
  select(-id) %>% 
  select(-gene.regulation) %>%
  left_join(., gene_regulation_ep300, by = "gene.name") %>% 
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.max.peak = mean(amplitude)) %>% 
  spread(., key = "time", value = "mean.max.peak") %>%
  as.data.frame() %>% 
  mutate(., mean.weighted.time = rowSums(t(t(.[5:16])*{.[5:16] %>% colnames() %>% as.numeric()/60}))/rowSums(.[,5:16])) %>% 
  select(gene.name, start.range, TF, gene.regulation, mean.weighted.time) %>% 
  filter(gene.regulation == "up-delta_ep300") %>% 
  mutate(mean.weighted.time = mean.weighted.time*60) %>% 
  select(-start.range) %>%
  rename(type.data = TF, value = mean.weighted.time) -> MWT_delta_ep300
  
#########################
# Plot for MCTP and MWT #
#########################
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/boxplot_MCTP_MWT_delta_ep300.svg", 
        width = 10,
        height = 8)

rbind(MCTP_delta_ep300, MWT_delta_ep300) %>% 
{ggplot(., aes(x = type.data, y = value, color = type.data)) + 
    geom_boxplot() + 
    theme(axis.title.x = element_blank())}

dev.off()

rm(MWT_delta_ep300, MCTP_delta_ep300)

###########################################################
# Function prepare dataframe for extract data to enhancer #
###########################################################
file_to_search_enhancer <- function(data_frame){
  data_frame %>% 
    select(-c(gene.regulation, amplitude, TF, time)) %>% #remove diff_min_max, code to check!!!
    unique() %>% 
    left_join(., {peaks_all_genes_without_promoters_amplitude %>% 
        select(gene.name, chromosome) %>% unique()}, by = "gene.name") %>% 
    na.omit() %>%
    mutate(chromosome = str_replace(.$chromosome, "chr", "")) %>% 
    left_join(., gene_regulation_ep300, by = "gene.name") %>%
    na.omit() %>%
    mutate(ensemblid = "NA") %>%#add select, change number of column 
    select(ensemblid, gene.name, chromosome, start.range, end.range, gene.regulation) 
}


###########################################################
# save to file information to finding enhancers for genes #
###########################################################
file_to_search_enhancer(delta_ep300) %>% 
  fwrite("~/ifpan-chipseq-timecourse/DATA/enhancer_delta_ep300_info.tsv", 
         sep="\t", 
         col.names = TRUE, 
         row.names = FALSE)

rm(delta_ep300, 
   to.heatmap, 
   number_clusters, 
   tmp_dend, 
   gene_regulation_ep300)
