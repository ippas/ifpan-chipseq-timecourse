###############################################
# making plot with peaks promoters for all TF #
###############################################
# jpeg("~/ifpan-chipseq-timecourse/PLOTS/lineplot_promotores.jpeg", 
# #jpeg("~/dexamethasone/lineplot_promotores.jpeg", 
#      width = 1400, 
#      height = 802)
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/lineplot_promotores.svg", 
        width = 10,
        height = 8)

read.table("~/ChIP-seq/DATA/promotores_peaks_value.tsv", 
           header = FALSE, 
           sep = "\t", 
           stringsAsFactors = FALSE) %>% 
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", 1:40)) %>% 
  gather(., "bucket.range", "value", -c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file")) %>% 
  group_by(bucket.range, time, TF, gene.regulation) %>%
  summarize(mean.value = mean(value)) %>% 
  {ggplot(., aes(x = as.numeric(bucket.range)*500, y = mean.value, color = as.factor(gene.regulation))) + 
      geom_line(size = 0.5) + 
      facet_grid(TF~time, scales = "free_y") +
      labs(fill = "Gene regulation") + 
      theme(axis.text.x = element_text(angle=45, hjust = 1, size = 14),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            strip.text.x = element_text(size = 16),
            strip.text.y = element_text(size = 10),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16),
            legend.position = "bottom") +
      #labs(fill = "Gene regulation") + 
      scale_color_manual(values = c("up-regulated" = "firebrick", 
                                    "random" = "gray", 
                                    "down-regulated" = "dodgerblue")) +
      scale_x_continuous(limits=c(0, 20001), 
                         breaks = c(1, 10001, 20001), 
                         labels= c("-10000", "0", "10000")) +
      
      ggtitle("Peaks for promoters")}

dev.off()

##################################################################
# making plot with peaks promoters for all TF - relative changes #
##################################################################
# jpeg("~/ifpan-chipseq-timecourse/PLOTS/lineplot_promotores_relative_changes.jpeg", 
#      width = 1400, 
#      height = 802)
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/lineplot_promotores_relative_changes.svg", 
        width = 10,
        height = 8)

read.table("~/ChIP-seq/DATA/promotores_peaks_value.tsv",
           header = FALSE,
           sep = "\t",
           stringsAsFactors = FALSE) %>%
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", 1:40)) %>%
  gather(., "bucket.range", "value", -c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file")) %>%
  group_by(bucket.range, time, TF, gene.regulation) %>%
  summarize(mean.value = mean(value)) %>%
  mutate(number.regulation=c("down"=1, "random"=3, "up"=2)) %>%
  mutate(control = ifelse(number.regulation == 3, 1, 0)) %>%
  mutate(max = max(mean.value * control)) %>% 
  mutate(relative.value = mean.value / max) %>%
  {ggplot(., aes(x = as.numeric(bucket.range)*500, y = relative.value, color = as.factor(gene.regulation))) + 
      geom_line(size = 0.5) + 
      facet_grid(TF~time, scales = "free_y") +
      theme(axis.text.x = element_text(angle=45, hjust = 1, size = 14),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            strip.text.x = element_text(size = 16),
            strip.text.y = element_text(size = 10),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16),
            legend.position = "bottom") +
      scale_color_manual(values = c("up-regulated" = "firebrick",
                                    "random" = "gray",
                                    "down-regulated" = "dodgerblue")) +
      scale_x_continuous(limits=c(0, 20001), 
                         breaks = c(1, 10001, 20001), 
                         labels= c("-10000", "0", "10000")) +
      ggtitle("Relative peak changes for promoters")}

dev.off()

##################
# Choose four TF #
##################
filtered_TF <- tmp_significant_random_genes_peak_normalized_amplitude %>% 
  select(TF) %>% 
  unique() %>% 
  .$TF %>% 
  sort %>% 
  .[5:14] %>% 
  .[-2] %>% 
  .[-7:-8] %>%
  .[-4:-6]

################################################
# making plot with peaks promoters for four TF #
################################################
# jpeg("~/ifpan-chipseq-timecourse/PLOTS/lineplot_promotores_fourTF.jpeg", 
#      width = 1400, 
#      height = 802)
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/lineplot_promotores_fourTF.svg", 
        width = 10,
        height = 8)

read.table("~/ChIP-seq/DATA/promotores_peaks_value.tsv", 
                                 header = FALSE, 
                                 sep = "\t", 
                                 stringsAsFactors = FALSE) %>% 
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", 1:40)) %>% 
  gather(., "bucket.range", "value", -c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file")) %>% 
  filter(TF %in% filtered_TF) %>% 
  group_by(bucket.range, time, TF, gene.regulation) %>%
  summarize(mean.value = mean(value)) %>% 
  {ggplot(., aes(x = as.numeric(bucket.range)*500, y = mean.value, color = as.factor(gene.regulation))) + 
      geom_line(size = 0.5) + 
      facet_grid(TF~time, scales = "free_y") +
      theme(axis.text.x = element_text(angle=45, hjust = 1, size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            strip.text.x = element_text(size = 16),
            strip.text.y = element_text(size = 16),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16),
            legend.position = "bottom") +
      scale_color_manual(values = c("up-regulated" = "firebrick", 
                                    "random" = "gray", 
                                    "down-regulated" = "dodgerblue")) +
      scale_x_continuous(limits=c(0, 20001), 
                         breaks = c(1, 10001, 20001), 
                         labels= c("-10000", "1", "10000")) +
      ggtitle("Peaks for promoters")}

dev.off()

###################################################################
# making plot with peaks promoters for four TF - relative changes #
###################################################################
# jpeg("~/ifpan-chipseq-timecourse/PLOTS/lineplot_promotores_relative_changes_fourTF.jpeg", 
#      width = 1400, 
#      height = 802)
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/lineplot_promotores_relative_changes_fourTF.svg", 
        width = 10,
        height = 8)

read.table("~/ChIP-seq/DATA/promotores_peaks_value.tsv",
                                 header = FALSE,
                                 sep = "\t",
                                 stringsAsFactors = FALSE) %>%
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", 1:40)) %>%
  gather(., "bucket.range", "value", -c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file")) %>%
  filter(TF %in% filtered_TF) %>% 
  group_by(bucket.range, time, TF, gene.regulation) %>%
  summarize(mean.value = mean(value)) %>%
  mutate(number.regulation=c("down"=1, "random"=3, "up"=2)) %>%
  mutate(control = ifelse(number.regulation == 3, 1, 0)) %>%
  mutate(max = max(mean.value * control)) %>% 
  mutate(relative.value = mean.value / max) %>%
  {ggplot(., aes(x = as.numeric(bucket.range)*500, y = relative.value, color = as.factor(gene.regulation))) + 
      geom_line(size = 0.5) + 
      facet_grid(TF~time, scales = "free_y") +
      theme(axis.text.x = element_text(angle=45, hjust = 1, size = 14),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            strip.text.x = element_text(size = 16),
            strip.text.y = element_text(size = 16),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16),
            legend.position = "bottom") +
      scale_color_manual(values = c("up-regulated" = "firebrick",
                                    "random" = "gray",
                                    "down-regulated" = "dodgerblue")) +
      scale_x_continuous(limits=c(0, 20001), 
                         breaks = c(1, 10001, 20001), 
                         labels= c("-10000", "1", "10000")) +
      ggtitle("Relative peak changes for promoters")}

dev.off()