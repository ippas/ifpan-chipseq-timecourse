#############################################################################
# loading the data for enhacer with birange (+/-10000 from the center peak) #
#############################################################################
enhancer_bigrange <- read.table("~/ChIP-seq/DATA/enhancer_bigrange_value.tsv",
                                    header = FALSE,
                                    sep = "\t",
                                    stringsAsFactors = FALSE) %>%
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", 1:2000)) %>%
  left_join(., {gene_chromosome_start_end_strand %>% filter(Gene.stable.ID  %in% {results %>% .[,2]})}, 
            by = c("gene.name" = "Gene.name")) %>% 
  mutate(pos.TSS=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% #marking position TSS for genes
  select(-c("Chromosome.scaffold.name", "Gene.start..bp.", "Gene.end..bp.", "Strand")) %>% # remove unnecessary columns
  mutate(diff.TSS.start.peak = abs(pos.TSS-start.range)) %>% # calculate disctance between start enhancer and position TSS
  group_by(start.range) %>% 
  filter(diff.TSS.start.peak == min(diff.TSS.start.peak)) %>% # choose gene which are the closest TSS
  as.data.frame() %>% 
  ungroup() %>%
  select(-c("Gene.stable.ID", "pos.TSS", "diff.TSS.start.peak")) %>%
  gather(., "bucket.range", "value", -c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file"))


###############################################
# making plot with peaks enhancer for four TF #
###############################################
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/lineplot_enhancer_fourTF.svg", 
        width = 10,
        height = 8)

enhancer_bigrange %>% 
  group_by(bucket.range, time, TF, gene.regulation) %>%
  summarize(mean.value = mean(value)) %>% 
  mutate(gene.regulation = factor(gene.regulation, c("up-regulated", "random", "down-regulated"))) %>% 
  ungroup() %>%  
  mutate(bucket.range = as.numeric(bucket.range)) %>%
  filter(bucket.range >= 800, bucket.range <= 1200)  %>% 
  {ggplot(., aes(x = bucket.range, y = mean.value, color = gene.regulation)) + 
      geom_line(size = 0.5) + 
      facet_grid(TF~time) +
      theme(legend.position = "bottom") +
      scale_color_manual(values = c("up-regulated" = "firebrick", 
                                    "random" = "gray", 
                                    "down-regulated" = "dodgerblue")) +
      theme(axis.text.x = element_text(angle=45, hjust = 1),
            legend.position = "bottom") +
      scale_x_continuous(limits=c(800, 1200),breaks = c(800, 1001, 120), labels = c("-2000","0", "2000")) +
      ggtitle("Peaks for enhancer without promoters")}

dev.off()


##################################################################
# making plot with peaks enhancer for four TF - relative changes #
##################################################################
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/lineplot_enhancer_relative_changes_fourTF.svg", 
        width = 10,
        height = 8)

enhancer_bigrange %>%
  group_by(bucket.range, time, TF, gene.regulation) %>%
  summarize(mean.value = mean(value)) %>%
  mutate(number.regulation=c("down-regulated"=1, "random"=3, "up-regulated"=2)) %>%
  mutate(control = ifelse(number.regulation == 3, 1, 0)) %>%
  mutate(max = max(mean.value * control)) %>% 
  mutate(relative.value = mean.value / max) %>%
  ungroup() %>%  
  mutate(bucket.range = as.numeric(bucket.range)) %>%
  filter(bucket.range >= 800, bucket.range <= 1200)  %>% 
  {ggplot(., aes(x = bucket.range, y = relative.value, color = as.factor(gene.regulation))) + 
      geom_line(size = 0.5) + 
      facet_grid(TF~time) +
      theme(legend.position = "bottom") +
      scale_color_manual(values = c("up-regulated" = "firebrick", 
                                    "random" = "gray", 
                                    "down-regulated" = "dodgerblue")) +
      theme(axis.text.x = element_text(angle=45, hjust = 1),
            legend.position = "bottom") +
      scale_x_continuous(limits=c(800, 1200),breaks = c(800, 1001, 1200), labels = c("-2000","0", "2000")) +
      ggtitle("Relative peak changes for enhancer without promoters")}

dev.off()

############################################################
# heatmapa bindding TF to DNA +/-2000 from the center peak #
############################################################
install.packages('tidytext')
library(tidytext)
library(reshape2)

for(list.vector in list(c("up-regulated", "firebrick", "Enhancers up-regulated", "up_regulated_plot"), 
                        c("down-regulated", "dodgerblue", "Enhancers down-regulated", "down_regulated_plot"), 
                        c("random", "gray", "Enhancers random", "random_plot"))){
  threshold = 5.0
  enhancer_bigrange %>% 
    mutate(bucket.range = as.numeric(bucket.range)) %>%
    filter(gene.regulation == list.vector[1]) %>%
    #filter(TF == "H3K27ac") %>%
    #group_by(time, TF) %>%
    mutate(scale.value = scale(value)) %>% 
    na.omit() %>%
    mutate(scale.value = {scale.value[scale.value > threshold] <- threshold; scale.value[scale.value < -threshold] <- -threshold; scale.value}) %>%
    ungroup() %>%
    mutate(value = scale.value) %>%
    filter(bucket.range >= 800, bucket.range <= 1200) -> tmp.heatmap
  
  assign(list.vector[4], {tmp.heatmap %>%
  {ggplot(., aes(x=bucket.range, y=gene.name, fill=value)) +
      geom_tile(aes(x=bucket.range, y=reorder(gene.name, value), fill=value)) +
      #scale_fill_gradient2(midpoint = 0, low ="blue", mid = "white", high = "red")+
      #scale_fill_gradient(low="white", high="red") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
      facet_grid(TF~time) +
      theme(axis.text.y = element_text(size = 2, colour = list.vector[2]),
            axis.text.x = element_text(angle=45, hjust = 1)) +
      scale_x_continuous(limits=c(800, 1200),breaks = c(800, 1001, 1200), labels = c("-2000","0", "2000")) +
      ggtitle(list.vector[3])
  }})
}
# 
# jpeg("~/ifpan-chipseq-timecourse/PLOTS/heatmap_enhancer_top_ep300.jpeg", 
#      width = 1400, 
#      height = 802)

png("~/ifpan-chipseq-timecourse/PLOTS/heatmap_enhancer.png", 
    width = 1400, 
    height = 802)

grid.arrange(up_regulated_plot, down_regulated_plot, random_plot, ncol=2)

dev.off()

rm(up_regulated_plot, down_regulated_plot, random_plot)

# spradzić zmienić nazwę i zmodyfikować skrypty w bashu
#########################################
# Changes in amplitude for GR and EP300 #
#########################################

gtf2 <- read.table('~/ifpan-chipseq-timecourse/DATA/Homo_sapiens.GRCh38.95.protein_coding.gtf', 
                   header = FALSE, sep = '\t', col.names = c("ensemblid", "gene.name"))


results %>%
  left_join(., {gene_chromosome_start_end_strand %>% filter(Gene.stable.ID  %in% {results %>% .[,2]})}, 
            by = c("ensemblid" = "Gene.stable.ID")) %>% 
  mutate(pos=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% 
  mutate(start = pos - 100000, end = pos + 100001) %>% 
  select("ensemblid", "gene.name", "Chromosome.scaffold.name", "start", "end") %>% 
  na.omit() %>%  
  left_join(., gtf2, by="gene.name") %>%
  na.omit() %>% 
  select(-ensemblid.y) %>%
  set_colnames(c("ensemblid", "gene.name", "chromosome", "start", "end")) %>% 
  mutate(gene.regulation = "NA") %>%
  mutate(start = ifelse(start < 0, 0, start)) %>% 
  fwrite('~/ifpan-chipseq-timecourse/DATA/range_all_genes.bed',
         sep="\t",
         col.names = TRUE,
         row.names = FALSE)

rm(gtf2, enhancer_bigrange)



