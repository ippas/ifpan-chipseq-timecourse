####################################poprawić kod, usunąć promotory po wczytaniu pliku
#  choose peaks without promoters  #
####################################
peaks_all_genes_ep300_nr3c1_amplitude <- read.table('~/ChIP-seq/DATA/peaks_all_genes_ep300_nr3c1_amplitude.tsv',
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
  left_join(., {peaks_all_genes_ep300_nr3c1_amplitude %>% 
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
  rownames_to_column(., var = "ensemblid") %>%
  left_join(., results[,c(2,4)], by = "ensemblid") %>%
  filter(gene.name %in% {delta_ep300 %>% select(gene.name) %>% unique %>% .[, 1]}) %>%
  column_to_rownames(., var = "gene.name") %>%
  select(-ensemblid) %>%
  as.matrix() %>%
  apply(1, scale) %>%
  t %>%
  apply(1, function(x, threshold){x[x > threshold] <- threshold; x[x < -threshold] <- -threshold; x}, threshold = 2.0) %>%
  t %>%
  {colnames(.) <- colnames(data); .} %>%
  as.data.frame() %>% na.omit() %>%
  as.matrix() -> to.heatmap
  
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



# #####################################################################
# # Boxplot, which show change amplitude for enhancer during the time #
# #####################################################################
# svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/boxplot_amplitude_delta_ep300.svg", 
#         width = 10,
#         height = 8)
# 
# 
# peaks_all_genes_ep300_nr3c1_amplitude %>%
#   select(-c(file)) %>% 
#   group_by(gene.name, start.range, end.range, TF, time, gene.regulation) %>% 
#   summarise(amplitude_mean = mean(amplitude)) %>%
#   rename(amplitude = amplitude_mean) %>%
#   ungroup %>%
#   mutate(id = paste(gene.name, start.range, end.range, sep = '.')) %>% 
#   filter(id %in% {delta_ep300 %>% #filter by delta_ep300
#       mutate(id = paste(gene.name, start.range, end.range, sep = '.')) %>%
#       .[, 8]}) %>% 
#   select(-id) %>% 
#   group_by(gene.name, start.range, time, TF) %>% 
#   summarise(mean.max.peak = mean(amplitude)) %>% 
#   left_join(., gene_regulation_ep300, by = "gene.name") %>%
#   na.omit %>%
#   {ggplot(., aes(x = as.factor(time), y = log(mean.max.peak), color = gene.regulation)) +
#       geom_boxplot(position = position_dodge(), outlier.size = 0) +  
#       facet_wrap(TF ~ ., ncol = 2, scales = "free_y" ) +
#       theme(legend.position = "bottom") +
#       scale_color_manual(values = c("up-delta_ep300" = "lightpink4", 
#                                     "down-delta_ep300" = "skyblue4")) +
#       labs(x = "Time [min]",
#            y = "Log amplitude values") +
#       ggtitle("Delta ep300")}
# 
# dev.off()


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

 
peaks_all_genes_ep300_nr3c1_amplitude %>%
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
  mutate(type.data = factor(.$type.data, levels = c("NR3C1", "expression", "EP300"))) %>%
  {ggplot(., aes(x = type.data, y = value, color = type.data)) + 
    geom_boxplot() + 
    theme(axis.title.x = element_blank())}

dev.off()

#####################################
# Statistic: ANOVA, pairwise.t.test #
#####################################
sink("~/ifpan-chipseq-timecourse/DATA/MWT_MCTP_delta_ep300_ANOVA.txt")
rbind(MCTP_delta_ep300, MWT_delta_ep300) %>%
  aov(value ~ type.data + Error(gene.name), data = .) %>%
  summary() %>% print()
sink()  

# calculate quantiles: Q1, Q2, Q3, and max and min value

rbind(MCTP_delta_ep300, MWT_delta_ep300) %>% 
  rename(group = type.data) -> tmp_MCTP_MWT_delta_ep300
tapply(tmp_MCTP_MWT_delta_ep300$value, tmp_MCTP_MWT_delta_ep300$group, summary) %>%
  lapply(., function(x) { x %>% 
      t %>% 
      t %>% 
      as.data.frame() %>% 
      select(-Var2) %>% 
      set_colnames(c("type_summary", "value"))}) %>%
  melt() %>% 
  select(-variable) %>% 
  rename(type_data = L1) %>%
  select(type_data, type_summary, value) %>% 
  spread(., key = "type_summary", value = "value") %>%
  fwrite("~/ifpan-chipseq-timecourse/DATA/MWT_MCTP_delta_ep300_basic_summary.tsv", 
         sep="\t", 
         col.names = TRUE, 
         row.names = FALSE)

rm(MWT_delta_ep300, MCTP_delta_ep300, tmp_MCTP_MWT_delta_ep300)


#########################################################
# timecourse graph, expression and combine EP300, NR3C1 #
#########################################################

gene_regulation_ep300 %>%
  left_join(., {genes.names %>% filter(Gene.stable.ID %in% {results %>% select(ensemblid) %>% unique %>% .[,1]})}, by = c("gene.name" = "Gene.name")) %>%
  select(-Gene.stable.ID.version) %>%
  rename(gene.id = Gene.stable.ID) -> gene_regulation_ep300


up_delta_ep300 <- data[((rownames(data) %>% str_split(., fixed("."), n = 2, simplify = TRUE) %>% .[, 1]) %in% gene_regulation_ep300$gene.id[which(gene_regulation_ep300$gene.regulation == "up-delta_ep300")]),] %>%
  apply(1, fc) %>% t  
colnames(up_delta_ep300) <- colnames(raw.data[3:48])


##################
# for expression #
##################

# ## PART 2 ##
# ## prepare data to plot ###
up_delta_ep300 %>% 
  t() %>%
  melt() %>% 
  mutate(time = Var1 %>% str_split(., "_", simplify = TRUE) %>% .[ ,1]) %>%
  mutate(id = rep(rownames(up_delta_ep300), each = ncol(up_delta_ep300))) %>%
  mutate(time = replace(time, time=="t00", 0) %>%
           replace(time=="t05", 0.5) %>%
           replace(time=="t1", 1) %>%
           replace(time=="t10", 10) %>%
           replace(time=="t12", 12) %>%
           replace(time=="t2", 2) %>%
           replace(time=="t3", 3) %>%
           replace(time=="t4", 4) %>%
           replace(time=="t5", 5) %>%
           replace(time=="t6", 6) %>%
           replace(time=="t7", 7) %>%
           replace(time=="t8", 8) %>% as.numeric) %>%
  aggregate(. ~ id * time, data = ., FUN = . %>% median) %>% 
  mutate(ensemblid = id %>% str_split(., fixed("."), simplify = TRUE) %>% .[ ,1]) %>% 
  select(ensemblid, time, value) %>%
  left_join(., {gene_chromosome_start_end_strand[, c(1, 5)] %>% filter(Gene.stable.ID  %in% {results %>% .[,2]})}, 
            by = c("ensemblid" = "Gene.stable.ID")) %>%
  select(time, value) -> spline.exp.delta.ep300


## PART 3 ##
## prepare boxplot ###

with(spline.exp.delta.ep300, boxplot(value ~ time))
colnames(spline.exp.delta.ep300) <- c("time", "z")
spline.exp.delta.ep300$group <- c(rep('expression', nrow(spline.exp.delta.ep300)))


#remove.outliers <- function(x) {!x %in% boxplot.stats(x)$out}
spline.exp.delta.ep300 <- spline.exp.delta.ep300 %>%
  group_by(time) %>%
  filter(remove.outliers(z))

## PART 4 ##
## prepare lines to plot - binding at the time for NR3C1 and EP300 and derivative for expression

degree = 3
mod.exp.delta.ep300 <- with(spline.exp.delta.ep300, lm(z ~ bs(time, degree = degree)))
pdat.exp.delta.ep300 <- with(spline.exp.delta.ep300, data.frame(time = seq(min(time), max(time), length = 100)))
pdat.exp.delta.ep300 <- transform(pdat.exp.delta.ep300, yhat = predict(mod.exp.delta.ep300, newdata = pdat.exp.delta.ep300))
pdat.exp.delta.ep300 <- pdat.exp.delta.ep300 %>% mutate(deriv = 100 * (lead(yhat) - yhat))
pdat.exp.delta.ep300$group <- c(rep('expression changes', nrow(pdat.exp.delta.ep300)))


## PART 5 ##
## plot boxplot and lines on top ##

svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/boxlineplot_derivative_expression_delta_ep300.svg",
        width = 10,
        height = 8)

ggplot(spline.exp.delta.ep300, aes(x = time, y = z)) + 
  geom_boxplot(data = spline.exp.delta.ep300, aes(group = time), outlier.shape = NA, alpha=0.4, color = 'grey20') +
  ylim(-1,5) +
  geom_line(inheret.aes = FALSE, data = pdat.exp.delta.ep300, aes(x = time, y = yhat), alpha=0.4, color = 'grey20') +
  geom_line(inheret.aes = FALSE, data = pdat.exp.delta.ep300, aes(x = time, y = deriv), size=1.5)

dev.off()

#######################
# for ep300 and nr3c1 #
#######################
spline.nr3c1.delta.ep300 <- peaks_all_genes_ep300_nr3c1_amplitude %>% 
  mutate(id = paste(gene.name, start.range, end.range, sep = "_")) %>% 
  filter(id %in% {delta_ep300 %>% mutate(id = paste(gene.name, start.range, end.range, sep = "_")) %>% select(id) %>% unique() %>% .[,1]}) %>%
  select(-c(id, gene.regulation)) %>%
  left_join(., gene_regulation_ep300, by = "gene.name") %>% 
  filter(gene.regulation == "up-delta_ep300") %>% 
  filter(TF == "NR3C1") %>% 
  mutate(id = paste(chromosome, start.range, end.range, sep = ".")) %>%
  select(id, time, amplitude) %>%
  mutate(time = time / 60) %>%
  group_by(id) %>%
  data.frame %>%
  aggregate(. ~ id * time, data = ., FUN = . %>% median) %>%
  group_by(id) %>% 
  mutate(z = scale(amplitude)) %>%
  ungroup() %>%
  select(time, z)


spline.ep300.delta.ep300 <- peaks_all_genes_ep300_nr3c1_amplitude %>% 
  mutate(id = paste(gene.name, start.range, end.range, sep = "_")) %>% 
  filter(id %in% {delta_ep300 %>% mutate(id = paste(gene.name, start.range, end.range, sep = "_")) %>% select(id) %>% unique() %>% .[,1]}) %>%
  select(-c(id, gene.regulation)) %>%
  left_join(., gene_regulation_ep300, by = "gene.name") %>% 
  filter(gene.regulation == "up-delta_ep300") %>% 
  filter(TF == "EP300") %>%
  mutate(id = paste(chromosome, start.range, end.range, sep = ".")) %>%
  select(id, time, amplitude) %>%
  mutate(time = time / 60) %>%
  group_by(id) %>%
  data.frame %>%
  aggregate(. ~ id * time, data = ., FUN = . %>% median) %>%
  group_by(id) %>% 
  mutate(z = scale(amplitude)) %>%
  ungroup() %>%
  select(time, z)

## PART 3 ##
## prepare boxplot ###
spline.ep300.delta.ep300$group <- c(rep('EP300', nrow(spline.ep300.delta.ep300)))
spline.nr3c1.delta.ep300$group <- c(rep('NR3C1', nrow(spline.nr3c1.delta.ep300)))

boxplot.data.delta.ep300 <- rbind(spline.nr3c1.delta.ep300, spline.ep300.delta.ep300) %>% mutate(id = paste(time, group, sep = "_"))

## PART 4 ##
## prepare lines to plot - binding at the time for NR3C1 and EP300 and derivative for expression

mod.gr.delta.ep300 <- with(spline.nr3c1.delta.ep300, lm(z ~ bs(time, degree = degree)))
pdat.gr.delta.ep300 <- with(spline.nr3c1.delta.ep300, data.frame(time = seq(min(time), max(time), length = 100)))
pdat.gr.delta.ep300 <- transform(pdat.gr.delta.ep300, z = scale(predict(mod.gr.delta.ep300, newdata = pdat.gr.delta.ep300)))

mod.ep.delta.ep300 <- with(spline.ep300.delta.ep300, lm(z ~ bs(time, degree = degree)))
pdat.ep.delta.ep300 <- with(spline.ep300.delta.ep300, data.frame(time = seq(min(time), max(time), length = 100)))
pdat.ep.delta.ep300 <- transform(pdat.ep.delta.ep300, z = scale(predict(mod.ep.delta.ep300, newdata = pdat.ep.delta.ep300)))

pdat.ep.delta.ep300$group <- c(rep('EP300', nrow(pdat.ep.delta.ep300)))
pdat.gr.delta.ep300$group <- c(rep('NR3C1', nrow(pdat.gr.delta.ep300)))

lineplot.data.delta.ep300 <- rbind(pdat.ep.delta.ep300, pdat.gr.delta.ep300)


svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/boxlineplot_changes_NR3C1_EP300_delta_ep300.svg",
        width = 10,
        height = 8)

ggplot() + 
  geom_boxplot(data = boxplot.data.delta.ep300, aes(x=as.numeric(time), y=z, fill=group, group = id), outlier.shape = NA, alpha=0.4, color = 'grey') +
  geom_line(data = lineplot.data.delta.ep300, aes(x = as.numeric(time), y = z, color = group), size=1.5)

dev.off()


rm(up_delta_ep300,
   degree, 
   mod.exp.delta.ep300,
   pdat.exp.delta.ep300,
   spline.nr3c1.delta.ep300,
   spline.ep300.delta.ep300,
   boxplot.data.delta.ep300,
   mod.gr.delta.ep300,
   pdat.gr.delta.ep300,
   mod.ep.delta.ep300,
   pdat.ep.delta.ep300,
   lineplot.data.delta.ep300)




###########################################################
# Function prepare dataframe for extract data to enhancer #
###########################################################
file_to_search_enhancer <- function(data_frame){
  data_frame %>% 
    select(-c(gene.regulation, amplitude, TF, time)) %>% #remove diff_min_max, code to check!!!
    unique() %>% 
    left_join(., {peaks_all_genes_ep300_nr3c1_amplitude %>% 
        select(gene.name, chromosome) %>% unique()}, by = "gene.name") %>% 
    na.omit() %>%
    mutate(chromosome = str_replace(.$chromosome, "chr", "")) %>% 
    left_join(., gene_regulation_ep300, by = "gene.name") %>%
    na.omit() %>%
    mutate(ensemblid = "NA") %>% #add select, change number of column 
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
