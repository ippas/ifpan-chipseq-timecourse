#################################################################################################
#####################################
# analysis enhancer for delta ep300 #
#####################################


###############################################
# read file with data for ehnacer delta ep300 #
###############################################
tmp.enhancer.bigrange.delta.ep300 <- read.table("~/ChIP-seq/DATA/enhancer_bigrange_delta_ep300.tsv",
                                              header = FALSE,
                                              sep = "\t",
                                              stringsAsFactors = FALSE) %>%
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", 1:2000)) %>%
  gather(., "bucket.range", "value", -c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file"))



###################################
# lineplot for enhancer top ep300 #
###################################
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/lineplot_enhancer_range_delta_ep300.svg", 
        width = 10,
        height = 8)

tmp.enhancer.bigrange.delta.ep300 %>% 
  group_by(bucket.range, time, TF, gene.regulation) %>%
  summarize(mean.value = mean(value)) %>% 
  ungroup() %>%  
  mutate(bucket.range = as.numeric(bucket.range)) %>%
  filter(bucket.range >= 800, bucket.range <= 1200)  %>% 
  {ggplot(., aes(x = bucket.range, y = mean.value, color = gene.regulation)) + 
      geom_line(size = 0.5) + 
      facet_grid(TF~time) +
      theme(legend.position = "bottom") +
      scale_color_manual(values = c("down-delta_ep300" = "skyblue4", 
                                    "up-delta_ep300" = "lightpink4")) +
      theme(axis.text.x = element_text(angle=45, hjust = 1),
            legend.position = "bottom") +
      scale_x_continuous(limits=c(800, 1200),breaks = c(800, 1001, 1200), labels = c("-2000","0", "2000")) +
      ggtitle("Peaks for enhancer delta ep300 without promoters")}

dev.off()

##################################
# heatmap for enhancer top ep300 #
##################################
for(list.vector in list(c("up-delta_ep300", "lightpink4", "Enhancers up-delta_ep300", "tmp_heatmap_up_delta_ep300"), 
                        c("down-delta_ep300", "skyblue4", "Enhancers down-delta_ep300", "tmp_heatmap_down_delta_ep300"))){
  tmp.enhancer.bigrange.delta.ep300 %>% 
    mutate(bucket.range = as.numeric(bucket.range)) %>%
    filter(gene.regulation == list.vector[1]) %>%
    #group_by(time, TF) %>%
    mutate(scale.value = scale(value)) %>% 
    mutate(scale.value = {scale.value[scale.value > threshold] <- threshold; scale.value[scale.value < -threshold] <- -threshold; scale.value}) %>%
    ungroup() %>%
    mutate(value = scale.value) %>%
    filter(bucket.range >= 800, bucket.range <= 1200) -> tmp.heatmap
  
  assign(list.vector[4], {tmp.heatmap %>%
  {ggplot(., aes(x=bucket.range, y=gene.name, fill=value)) +
      geom_tile(aes(x=bucket.range, y=reorder(gene.name, value), fill=value)) +
      #scale_fill_gradient(low="white", high="red") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
      facet_grid(TF~time) +
      theme(axis.text.y = element_text(size = 3, colour = list.vector[2]),
            axis.text.x = element_text(angle=45, hjust = 1)) +
      scale_x_continuous(limits=c(800, 1200),breaks = c(800, 1001, 1200), labels = c("-2000","0", "2000")) +
      ggtitle(list.vector[3])
  }})
}

png("~/ifpan-chipseq-timecourse/PLOTS/heatmap_enhancer_delta_ep300.png", 
     width = 1400, 
     height = 802)

grid.arrange(tmp_heatmap_up_delta_ep300, tmp_heatmap_down_delta_ep300, ncol=2)
dev.off()

rm(tmp_heatmap_up_delta_ep300, tmp_heatmap_down_delta_ep300, tmp.heatmap, tmp.enhancer.bigrange.delta.ep300)

#################################################################################################