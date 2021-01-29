enhancer_amplitude <- read.table("~/ChIP-seq/DATA/enhancer_amplitude_value.tsv", 
                                                                     header = FALSE, 
                                                                     sep = "\t", 
                                                                     stringsAsFactors = FALSE) %>% 
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", "amplitude")) 


##################################################
# function that selects the maximum peak per gen #
##################################################
choose_strongest_peak <- function(data_frame){
  merge(x = data_frame,
        y = {data_frame %>%
            group_by(gene.name, chromosome, start.range, end.range, TF, gene.regulation) %>%
            summarise(tmp_mean_alltime_amplitude = mean (amplitude))}) %>%
    group_by(gene.name, TF, gene.regulation) %>% 
    filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>%
    ungroup() %>%
    select(-tmp_mean_alltime_amplitude)
}


############################################################
# making boxplot amplitude changes for the strongest peaks #
############################################################
jpeg("~/ifpan-chipseq-timecourse/PLOTS/boxplot_enhancer_amplitude.jpeg",
#jpeg("~/dexamethasone/boxplot_enhancer_amplitude.jpeg", 
     width = 1400, 
     height = 802)

#making boxplot for strongest peak for each gene
choose_strongest_peak(enhancer_amplitude) %>%
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.max.peak = mean(amplitude)) %>% 
  {ggplot(., aes(x = as.factor(time), y = log(mean.max.peak), color = gene.regulation)) +
      geom_boxplot(position = position_dodge(), outlier.size = 0) +  
      facet_wrap(TF ~ ., ncol = 4, scales = "free_y" ) +
      theme(legend.position = "bottom") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
           axis.text.y = element_text(size = 16),
           axis.title.x = element_text(size = 22),
           axis.title.y = element_text(size = 22),
           strip.text.x = element_text(size = 16),
           legend.title = element_text(size = 20),
           legend.text = element_text(size = 18),
           legend.position = "bottom") +
      scale_color_manual(values = c("up-regulated" = "firebrick",
                                    "random" = "gray",
                                    "down-regulated" = "dodgerblue")) +
      labs(x = "Time [min]",
           y = "Logarithmic mean for the amplitude") +
      ggtitle("Strongest peak for gene")}

dev.off()


#########################################################################
# making boxplot amplitude changes for the strongest peaks, for four TF #
#########################################################################
jpeg("~/ifpan-chipseq-timecourse/PLOTS/boxplot_enhancer_amplitude_fourTF.jpeg",
#jpeg("~/dexamethasone/boxplot_enhancer_amplitude_fourTF.jpeg", 
     width = 1400, 
     height = 802)

choose_strongest_peak(enhancer_amplitude) %>%
  filter(TF %in% filtered_TF) %>% 
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.max.peak = mean(amplitude)) %>% 
  {ggplot(., aes(x = as.factor(time), y = log(mean.max.peak), color = gene.regulation)) +
      geom_boxplot(position = position_dodge(), outlier.size = 0) +  
      facet_wrap(TF ~ ., ncol = 4, scales = "free_y" ) +
      theme(legend.position = "bottom") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.x = element_text(size = 22),
            axis.title.y = element_text(size = 22),
            strip.text.x = element_text(size = 16),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 18),
            legend.position = "bottom") +
      scale_color_manual(values = c("up-regulated" = "firebrick",
                                    "random" = "gray",
                                    "down-regulated" = "dodgerblue")) +
      labs(x = "Time [min]",
           y = "Logarithmic mean for the amplitude") +
      ggtitle("Strongest peak for gene")}

dev.off()

#############################################################
# making boxplot for mean weighted time for strongest peaks #
#############################################################
jpeg("~/ifpan-chipseq-timecourse/PLOTS/boxplot_enhancer_MWT.jpeg",
#jpeg("~/dexamethasone/boxplot_enhancer_MWT.jpeg", 
     width = 1400, 
     height = 802)

choose_strongest_peak(enhancer_amplitude) %>%
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.max.peak = mean(amplitude)) %>% 
  spread(., key = "time", value = "mean.max.peak") %>%
  as.data.frame() %>% 
  mutate(., mean.weighted.time = rowSums(t(t(.[5:16])*{.[5:16] %>% colnames() %>% as.numeric()/60}))/rowSums(.[,5:16])) %>% #calculate with time = 0
  select(gene.name, start.range, TF, gene.regulation, mean.weighted.time) %>% 
  {ggplot(., aes(x = gene.regulation, y = mean.weighted.time, color = gene.regulation)) +
      geom_boxplot() +
      theme(axis.title.x = element_blank()) +
      facet_wrap(.~TF, scales = "free_x", ncol =8) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(size = 16),
            axis.title.x = element_text(size = 22),
            axis.title.y = element_text(size = 22),
            strip.text.x = element_text(size = 16),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 18)) +
      scale_color_manual(values = c("up-regulated" = "firebrick",
                                    "random" = "gray",
                                    "down-regulated" = "dodgerblue")) +
      labs(y = "Mean weighted time",
           x = "Gene regulation") +
      ggtitle("Mean weighted time for strongest peak")}

dev.off()

#######################################################
# Mean weighted time for four TF, for strongest peaks #
#######################################################
jpeg("~/ifpan-chipseq-timecourse/PLOTS/boxplot_enhancer_MWT_fourTF.jpeg",
#jpeg("~/dexamethasone/boxplot_enhancer_MWT_fourTF.jpeg", 
     width = 1400, 
     height = 802)

choose_strongest_peak(enhancer_amplitude) %>% 
  filter(TF %in% filtered_TF) %>% 
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.max.peak = mean(amplitude)) %>% 
  spread(., key = "time", value = "mean.max.peak") %>%
  as.data.frame() %>% 
  mutate(., mean.weighted.time = rowSums(t(t(.[5:16])*{.[5:16] %>% colnames() %>% as.numeric()/60}))/rowSums(.[,5:16])) %>% #calculate with time = 0
  select(gene.name, start.range, TF, gene.regulation, mean.weighted.time) %>% 
  filter(gene.regulation == "up-regulated") %>% 
  {ggplot(., aes(x = TF, y = mean.weighted.time, color = TF)) +
      geom_boxplot() +
      theme(axis.title.x = element_blank()) +
      theme(axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.x = element_text(size = 22),
            axis.title.y = element_text(size = 22),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 18)) +
      labs(y = "Mean weighted time [hour]",
           x = "Transcription factor") +
      ggtitle("Mean weighted time for strongest peak")}

dev.off()
