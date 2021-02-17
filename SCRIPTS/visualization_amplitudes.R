
##################################################
# function that selects the maximum peak per gen #
##################################################
choose_strongest_peak <- function(data_frame){
  data_frame %>%
    mutate(id = paste(gene.name, start.range, end.range, sep = '.')) %>% 
    filter(id %in% {
      merge(x = data_frame,
            y = {data_frame %>% 
                group_by(gene.name, start.range, TF, gene.regulation) %>%
                summarise(tmp_mean_alltime_amplitude = mean (amplitude))}) %>% filter(TF == "NR3C1") %>% 
        group_by(gene.name, TF, gene.regulation) %>%  
        filter(tmp_mean_alltime_amplitude == max(tmp_mean_alltime_amplitude)) %>% 
        ungroup() %>% mutate(id = paste(gene.name, start.range, end.range, sep = '.')) %>% 
        as.data.frame() %>%.[,11]}) %>%
    select(-id)
}


remove_peak_promoter <- function(data_frame){
  data_frame %>% 
    left_join(., {gene_chromosome_start_end_strand %>% filter(Gene.stable.ID %in% {results %>% .[,2]})}, by = c("gene.name" = "Gene.name")) %>%
    select(-Chromosome.scaffold.name) %>%
    mutate(pos.TSS=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>%
    mutate(start.range.promoter = pos.TSS - 2000, end.range.promoter = pos.TSS + 2001) %>% 
    filter(end.range < start.range.promoter | start.range > end.range.promoter) %>% 
    select(c(gene.name, chromosome, start.range, end.range, gene.regulation, TF, time, file, amplitude))
}

# zminić, nazwę tak aby móc korzystać z genów etc
enhancer_amplitude <- read.table("~/ChIP-seq/DATA/enhancer_amplitude_value.tsv", 
                                 header = FALSE, 
                                 sep = "\t", 
                                 stringsAsFactors = FALSE) %>% 
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", "amplitude")) %>% 
  remove_peak_promoter() %>%
  choose_strongest_peak()


tmp_enhancer_amplitude <- enhancer_amplitude %>% mutate(gene.name = "NA") %>% unique 

# ############################################################
# # making boxplot amplitude changes for the strongest peaks #
# ############################################################
# svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/boxplot_enhancer_amplitude.svg", 
#         width = 10,
#         height = 8)
# 
# #making boxplot for strongest peak for each gene
# tmp_enhancer_amplitude %>%
#   group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
#   summarise(mean.max.peak = mean(amplitude)) %>% ungroup() %>%  
#   {ggplot(., aes(x = as.factor(time), y = log(mean.max.peak), color = gene.regulation)) +
#       geom_boxplot(position = position_dodge(), outlier.size = 0) +  
#       facet_wrap(TF ~ ., ncol = 4, scales = "free_y" ) +
#       theme(legend.position = "bottom") +
#       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
#            axis.text.y = element_text(size = 16),
#            axis.title.x = element_text(size = 22),
#            axis.title.y = element_text(size = 22),
#            strip.text.x = element_text(size = 16),
#            legend.title = element_text(size = 20),
#            legend.text = element_text(size = 18),
#            legend.position = "bottom") +
#       scale_color_manual(values = c("up-regulated" = "firebrick",
#                                     "random" = "gray",
#                                     "down-regulated" = "dodgerblue")) +
#       labs(x = "Time [min]",
#            y = "Logarithmic mean for the amplitude") +
#       ggtitle("Strongest peak for gene")}
# 
# dev.off()


#########################################################################
# making boxplot amplitude changes for the strongest peaks, for four TF #
#########################################################################
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/boxplot_enhancer_amplitude_fourTF.svg", 
        width = 10,
        height = 8)


tmp_enhancer_amplitude %>%
  filter(TF %in% filtered_TF) %>% 
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.max.peak = mean(amplitude)) %>% 
  {ggplot(., aes(x = as.factor(time), y = log(mean.max.peak), color = gene.regulation)) +
      geom_boxplot(position = position_dodge(), outlier.size = 0) +  
      facet_wrap(TF ~ ., ncol = 4, scales = "free_y") +
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

# #############################################################
# # making boxplot for mean weighted time for strongest peaks #
# #############################################################
# svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/boxplot_enhancer_MWT.svg", 
#         width = 10,
#         height = 8)
# 
# tmp_enhancer_amplitude %>%
#   group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
#   summarise(mean.max.peak = mean(amplitude)) %>% 
#   spread(., key = "time", value = "mean.max.peak") %>%
#   as.data.frame() %>% 
#   mutate(., mean.weighted.time = rowSums(t(t(.[5:16])*{.[5:16] %>% colnames() %>% as.numeric()/60}))/rowSums(.[,5:16])) %>% #calculate with time = 0
#   select(gene.name, start.range, TF, gene.regulation, mean.weighted.time) %>% 
#   {ggplot(., aes(x = gene.regulation, y = mean.weighted.time, color = gene.regulation)) +
#       geom_boxplot() +
#       theme(axis.title.x = element_blank()) +
#       facet_wrap(.~TF, scales = "free_x", ncol =8) +
#       theme(axis.text.x = element_blank(),
#             axis.text.y = element_text(size = 16),
#             axis.title.x = element_text(size = 22),
#             axis.title.y = element_text(size = 22),
#             strip.text.x = element_text(size = 16),
#             legend.title = element_text(size = 20),
#             legend.text = element_text(size = 18)) +
#       scale_color_manual(values = c("up-regulated" = "firebrick",
#                                     "random" = "gray",
#                                     "down-regulated" = "dodgerblue")) +
#       labs(y = "Mean weighted time",
#            x = "Gene regulation") +
#       ggtitle("Mean weighted time for strongest peak")}
# 
# dev.off()

#######################################################
# Mean weighted time for four TF, for strongest peaks #
#######################################################
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/boxplot_enhancer_MWT_fourTF.svg", 
        width = 10,
        height = 8)


tmp_enhancer_amplitude %>% 
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

rm(tmp_enhancer_amplitude)

#########################################################
# timecourse graph, expression and combine EP300, NR3C1 #
#########################################################

gene.regulation <- gene_regulation %>%
  left_join(., genes.names, by = c("gene.name" = "Gene.name")) %>%
  select(-Gene.stable.ID.version) %>%
  rename(gene.id = Gene.stable.ID)


upregulated <- data[((rownames(data) %>% str_split(., fixed("."), n = 2, simplify = TRUE) %>% .[, 1]) %in% gene.regulation$gene.id[which(gene.regulation$gene.regulation == "up-regulated")]),] %>%
  apply(1, scale) %>% t
colnames(upregulated) <- colnames(raw.data[3:48])

genes.names$Gene.stable.ID.version[match(gene.regulation$gene.name, genes.names$Gene.name)] %>%
  str_split(., fixed("."), n = 2, simplify = TRUE) %>% .[, 1]



# tmp_significant_random_genes_peak_normalized_amplitude <- read.table("/home/mateusz/ChIP-seq/DATA/significant_random_genes_chip-seq_normalized_gene_chromosome_start-peak_end-peak_TF_time_file_amplitude.tsv",
#                                                                      header = FALSE,
#                                                                      sep = "\t",
#                                                                      stringsAsFactors = FALSE) %>%
#   set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", "amplitude")) %>%
#   #replace gene.regulation from 4 cluster to 2 - this is to remove in the future
#   mutate(gene.regulation = {.$gene.regulation %>%
#       str_split(., "-", n = 2, simplify = TRUE) %>%
#       .[,1]}) %>%
#   left_join(., {gene_chromosome_start_end_strand %>% filter(Gene.stable.ID %in% {results %>% .[,2]})}, by = c("gene.name" = "Gene.name")) %>%
#   select(-Chromosome.scaffold.name) %>%
#   mutate(pos.TSS=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>%
#   mutate(start.range.promoter = pos.TSS - 2000, end.range.promoter = pos.TSS + 2001) %>%
#   filter(end.range < start.range.promoter | start.range > end.range.promoter) %>%
#   select(c(gene.name, chromosome, start.range, end.range, gene.regulation, TF, time, file, amplitude)) %>% 
#   choose_strongest_peak() %>% 
#   ungroup #%>%
#   #select(-c( tmp_mean_alltime_amplitude))

## PART 2 ##
## prepare data to plot ###

spline.data <- upregulated %>%
  t() %>%
  melt() %>% 
  mutate(time = Var1 %>% str_split(., "_", simplify = TRUE) %>% .[ ,1]) %>%
  mutate(id = rep(rownames(upregulated), each = 46)) %>%
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
  select(time, value)


spline.nr3c1 <- enhancer_amplitude %>%
  filter(gene.regulation == "up-regulated") %>%
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

spline.ep300 <- enhancer_amplitude %>%
  filter(gene.regulation == "up-regulated") %>%
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

colnames(spline.data) <- c("time", "z")

spline.data$group <- c(rep('expression', nrow(spline.data)))
spline.ep300$group <- c(rep('EP300', nrow(spline.ep300)))
spline.nr3c1$group <- c(rep('NR3C1', nrow(spline.nr3c1)))

boxplot.data <- rbind(spline.data, spline.nr3c1, spline.ep300)

## PART 4 ##
## prepare lines to plot - binding at the time for NR3C1 and EP300 and derivative for expression


mod.exp <- with(spline.data, lm(z ~ bs(time, degree = 3)))
pdat.exp <- with(spline.data, data.frame(time = seq(min(time), max(time), length = 100)))
pdat.exp <- transform(pdat.exp, yhat = predict(mod.exp, newdata = pdat.exp))
pdat.exp <- pdat.exp %>% mutate(deriv = scale(lead(yhat) - yhat))

pdat.exp %>% select(time, deriv) -> pdat.exp

colnames(pdat.exp) <-c("time", 'z')

degree = 3

mod.gr <- with(spline.nr3c1, lm(z ~ bs(time, degree = degree)))
pdat.gr <- with(spline.nr3c1, data.frame(time = seq(min(time), max(time), length = 100)))
pdat.gr <- transform(pdat.gr, z = scale(predict(mod.gr, newdata = pdat.gr)))

mod.ep <- with(spline.ep300, lm(z ~ bs(time, degree = degree)))
pdat.ep <- with(spline.ep300, data.frame(time = seq(min(time), max(time), length = 100)))
pdat.ep <- transform(pdat.ep, z = scale(predict(mod.ep, newdata = pdat.ep)))

pdat.exp$group <- c(rep('expression changes', nrow(pdat.exp)))
pdat.ep$group <- c(rep('EP300', nrow(pdat.ep)))
pdat.gr$group <- c(rep('NR3C1', nrow(pdat.gr)))

lineplot.data <- rbind(pdat.exp, pdat.ep, pdat.gr)

## PART 5 ##
## plot boxplot and lines on top ##

svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/boxlineplot_expression_NR3C1_EP300.svg", 
        width = 10,
        height = 8)

ggplot() + 
  geom_boxplot(data = boxplot.data, aes(x=as.factor(time), y=z, fill=group), outlier.shape = NA, alpha=0.4, color = 'grey') +
  geom_line(data = lineplot.data, aes(x = time, y = z, color = group), size=1.5)

dev.off()

rm(gene.regulation, 
   upregulated, 
   spline.data, 
   spline.nr3c1, 
   spline.ep300, 
   boxplot.data,
   mod.exp, pdat.exp,
   mod.gr, pdat.gr,
   mod.ep, pdat.ep,
   degree,
   lineplot.data)

   

##############################
# distance enhancer from TSS #
##############################

read.table("~/ChIP-seq/DATA/enhancer_peaks_value.tsv",
           header = FALSE,
           sep = "\t",
           col.names = c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", seq_len((number.column-8)*2)),
           stringsAsFactors = FALSE,
           fill=TRUE) %>% 
  #remove_peak_promoter() %>%
  mutate(id = paste(gene.name, start.range, end.range, sep = '.')) %>%
  filter(id %in% {enhancer_amplitude  %>% mutate(id = paste(gene.name, start.range, end.range, sep = '.')) %>% .[,10]}) %>%
  select(-id) %>%
  mutate(index.max.value = apply(., 1, function(x){as.numeric(x[9:(number.column*2-8)], na.rm=T) %>% which.max})) %>%
  select(gene.name, chromosome, start.range, end.range, TF, gene.regulation,  time, index.max.value) %>%
  #to remove filter because in the future one TF
  filter(TF == "NR3C1", time != 0) %>%
  group_by(gene.name, chromosome, start.range, end.range, TF, gene.regulation) %>% 
  summarise(index.max.value.mean = {mean(index.max.value) %>% round}) %>% 
  rename(index.max.value = index.max.value.mean) %>% 
  ungroup() -> tmp.data.enhancer.range 
  
random %>% 
  select(ensemblid, gene.name) %>%
  left_join(., gene_chromosome_start_end_strand, by = c("ensemblid" = "Gene.stable.ID")) %>%
  mutate(pos=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% 
  select(ensemblid, gene.name, Chromosome.scaffold.name , pos) %>%
  mutate(gene.regulation = "random") %>%
  rbind(., results[order(results$pvalue)[1:number_signification_genes],] %>% 
          left_join(., gene_chromosome_start_end_strand, by = c("ensemblid" = "Gene.stable.ID")) %>% 
          mutate(pos=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>%
          select(ensemblid, gene.name, Chromosome.scaffold.name , pos) %>%
          left_join(., gene_regulation, by = "gene.name")) %>% na.omit() %>% 
  set_colnames(c("ensemblid", "gene.name", "chromosome", "pos_TSS", "gene.regulation")) %>% 
  left_join(., {tmp.data.enhancer.range %>%
      unique() %>%
      mutate(position.amplitude = start.range + index.max.value)  %>% 
      select(gene.name, position.amplitude) %>% 
      unique()}, by = "gene.name") %>% 
  na.omit() -> enhancer.distance.from.TSS


svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/histogram_distance_enhancer_TSS.svg", 
        width = 10,
        height = 8)
enhancer.distance.from.TSS %>%  
  .[,3:6] %>% 
  unique() %>% 
  mutate(distance_to_TSS = pos_TSS - position.amplitude) %>%
  {ggplot(.,aes(x = distance_to_TSS)) + 
      geom_histogram(aes(fill = gene.regulation), position="identity", alpha=0.4, bins = 50) +
      scale_color_manual(values = c("up-regulated" = "firebrick",
                                    "random" = "gray20",
                                    "down-regulated" = "dodgerblue")) +
      scale_fill_manual(values = c("up-regulated" = "firebrick",
                                   "random" = "gray",
                                   "down-regulated" = "dodgerblue")) +
      scale_x_continuous(labels = comma) +
      ggtitle("Enhance distance from TSS")
  } 

dev.off()

rm(enhancer.distance.from.TSS, tmp.data.enhancer.range)
###################################################################################


###################################################################################
# do sprawdzenia co jest potrzebne
####################################################
# centering peaks for top ep300 without promotores #
####################################################

#aligned to the highest value

results.list <- lapply(readLines('/home/mateusz/ChIP-seq/DATA/enhancer_peaks_value.tsv'), strsplit, split = '\t')
lapply(results.list, lengths) -> results.lenghts

number.column <- results.list %>% lapply(lengths) %>% unlist() %>% max() # result = number.column - so length of table: (number.column x 2) + 1 = number.column x 2 + 1 so the middle value is number.column - this is where we put max

# adjust_numeric_columns <- function(x) {
#   
#   new.line <- integer(number.column * 2) # make a vector of zeros
#   x %>% unlist() -> single.row # unlist the line
#   single.row[c(9:length(single.row))] -> peak #get just the numeric valuse
#   peak %>% as.numeric() %>% which.max() -> amp.index #find where the max value of the peak is
#   new.line[c((number.column - amp.index + 1): (number.column + (length(peak) - amp.index)))] <- peak #put the peak into the new.line with the max set to index 1144
#   full.line <- c(single.row[c(1:8)], new.line)
#   return(full.line)
# }
# 
# results.list %>% lapply(adjust_numeric_columns) %>% data.frame(row.names = NULL) -> adjusted.peaks.table
# 
# names(adjusted.peaks.table) <- NULL
# 
# adjusted.peaks.table %>% t() %>% data.frame() -> adjusted.peaks.table
# 
# adjusted.peaks.table %>%
#   set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", "position.amplitude", 1:(number.column*2-1))) %>%  
#   gather(., "range", "value", -c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file")) -> gather.adjusted.peaks.table 
# 
# 
# #line plot for enhancer peak without promoters
# gather.adjusted.peaks.table %>% 
#   mutate(range = as.numeric(range), value = as.numeric(value)) %>%
#   group_by(gene.name, gene.regulation, TF, time, range) %>%
#   summarise(value = mean(value)) %>%
#   ungroup %>%
#   group_by(gene.name, gene.regulation, TF, time, range) %>%
#   summarize(value = mean(value)) %>%
#   group_by(gene.regulation, TF, time, range) %>%
#   summarize(value = mean(value)) -> tmp.to.plot
# tmp.to.plot %>% 
#   ungroup() %>%
#   mutate(time = as.numeric(as.character(time))) %>% 
#   {ggplot(., aes(x = range, y = value, color = gene.regulation)) +
#       geom_line() +
#       # scale_color_manual(values = c("up" = "firebrick",
#       #                               "random" = "gray",
#       #                               "down" = "dodgerblue")) +
#       facet_grid(TF~time)}
# 
# rm(gather.adjusted.peaks.table, gather.adjusted.peaks.table)


########################################################################
# prepare +/-10000 range for amplitude enhancer NR3C1 and save to file #
########################################################################

read.table("~/ChIP-seq/DATA/enhancer_peaks_value.tsv",
           header = FALSE,
           sep = "\t",
           col.names = c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", seq_len((number.column-8)*2)),
           stringsAsFactors = FALSE,
           fill=TRUE) %>% 
  mutate(id = paste(gene.name, start.range, end.range, sep = '.')) %>%
  filter(id %in% {enhancer_amplitude  %>% mutate(id = paste(gene.name, start.range, end.range, sep = '.')) %>% .[,10]}) %>%
  select(-id) %>%
  mutate(index.max.value = apply(., 1, function(x){as.numeric(x[9:(number.column*2-8)], na.rm=T) %>% which.max})) %>%
  select(gene.name, chromosome, start.range, end.range, TF, gene.regulation,  time, index.max.value) %>%
  #to remove filter because in the future one TF
  filter(TF == "NR3C1", time != 0) %>%
  group_by(gene.name, chromosome, start.range, end.range, TF, gene.regulation) %>% 
  summarise(index.max.value.mean = {mean(index.max.value) %>% round}) %>% 
  rename(index.max.value = index.max.value.mean) %>% 
  ungroup() %>%
  mutate(position.amplitude = start.range + index.max.value) %>% 
  mutate(start = position.amplitude - 10000, end = position.amplitude + 10001, ensemblid = "NA") %>%
  select(ensemblid, gene.name, chromosome, start, end, gene.regulation) %>% 
  mutate(chromosome = str_replace(.$chromosome, "chr", "")) %>%
  fwrite("~/ifpan-chipseq-timecourse/DATA/enhancer_bigrange_info.tsv", 
         sep="\t", 
         col.names = TRUE, 
         row.names = FALSE)



#########################################################################
# boxplot for Max Change Time Point (MCTP) and Mean Weighted Time (MWT) #
#########################################################################  

enhancer_amplitude %>% 
  filter(TF %in% filtered_TF) %>% 
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.max.peak = mean(amplitude)) %>% 
  spread(., key = "time", value = "mean.max.peak") %>%
  as.data.frame() %>% 
  mutate(., mean.weighted.time = rowSums(t(t(.[5:16])*{.[5:16] %>% colnames() %>% as.numeric()/60}))/rowSums(.[,5:16])) %>% #calculate with time = 0
  select(gene.name, start.range, TF, gene.regulation, mean.weighted.time) %>% 
  filter(gene.regulation == "up-regulated") %>%
  left_join(., {gene_chromosome_start_end_strand %>% filter(Gene.stable.ID  %in% {results %>% .[,2]})}, 
            by = c("gene.name" = "Gene.name")) %>% 
  mutate(pos=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% 
  select(-c("Chromosome.scaffold.name", "Gene.start..bp.", "Gene.end..bp.", "Strand")) %>%
  rename(pos.TSS = pos) %>% 
  mutate(diff.TSS.start.peak = abs(pos.TSS-start.range)) %>% 
  group_by(start.range) %>%
  filter(diff.TSS.start.peak == min(diff.TSS.start.peak)) %>% 
  as.data.frame() %>% 
  select(gene.name, start.range, TF, gene.regulation, mean.weighted.time) %>%
  filter(gene.regulation == "up-regulated") %>% 
  mutate(mean.weighted.time = mean.weighted.time*60) %>%
  select(-start.range) %>%
  rename(value = mean.weighted.time, group = TF) %>%  
  rbind(.,
        {max.change.time.point.significant %>% 
            filter(gene.regulation == "up-regulated") %>%
            mutate(group = "expression") %>%
            rename(value = max_change_time_point) %>% 
            select(c(gene.name, group, gene.regulation, value)) %>% 
            as.data.frame() %>%
            filter(gene.name %in% {tmp_filtered_enhancer %>% .[,1]})}) -> tmp_combine_MWT_MCTP


#####################
# plot MWT and MCTP #
#####################
svglite(file = "~/ifpan-chipseq-timecourse/PLOTS/boxplot_MWT_MCTP.svg", 
        width = 10,
        height = 8)

tmp_combine_MWT_MCTP %>% 
  mutate(group = factor(.$group, levels = c("NR3C1", "H3K4me1", "H3K27ac", "expression", "EP300"))) %>%
            {ggplot(., aes(x = group, y = value, color = group)) +
                geom_boxplot() +
                labs(x = "Gene regulation",
                     y = "Time [min]")}
dev.off()


#####################################
# Statistic: ANOVA, pairwise.t.test #
#####################################
tmp_combine_MWT_MCTP %>%
  aov(value ~ group + Error(gene.name), data = .) %>%
  summary()


pairwise.t.test(tmp_combine_MWT_MCTP$value, tmp_combine_MWT_MCTP$group, p.adjust.method = "none") %>% 
  .[3] %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "group1") %>% 
  gather(., "group2", "p.value", -group1) %>%
  mutate(group2 = str_replace(.$group2, "p.value.", "")) %>%
  filter(group1 == "expression" | group1 == "NR3C1" | group2 == "expression" | group2 == "NR3C1") %>% 
  na.omit() %>%
  mutate(p.value = p.adjust(.$p.value, method = "bonferroni"))


########################################
# heatmap with label GR, NOGR and NUGR #
########################################

# gene GR
enhancer_amplitude %>% 
  filter(TF %in% filtered_TF) %>% 
  group_by(gene.name, start.range, time, TF, gene.regulation) %>% 
  summarise(mean.max.peak = mean(amplitude)) %>% 
  spread(., key = "time", value = "mean.max.peak") %>%
  as.data.frame() %>% 
  mutate(., mean.weighted.time = rowSums(t(t(.[5:16])*{.[5:16] %>% colnames() %>% as.numeric()/60}))/rowSums(.[,5:16])) %>% #calculate with time = 0
  select(gene.name, start.range, TF, gene.regulation, mean.weighted.time) %>% 
  filter(gene.regulation != "random") %>%
  left_join(., {gene_chromosome_start_end_strand %>% filter(Gene.stable.ID  %in% {results %>% .[,2]})}, 
            by = c("gene.name" = "Gene.name")) %>% 
  mutate(pos=Gene.start..bp. * (Strand == 1) + Gene.end..bp. * (Strand == -1)) %>% 
  select(-c("Chromosome.scaffold.name", "Gene.start..bp.", "Gene.end..bp.", "Strand")) %>%
  rename(pos.TSS = pos) %>% 
  mutate(diff.TSS.start.peak = abs(pos.TSS-start.range)) %>% 
  group_by(start.range) %>%
  filter(diff.TSS.start.peak == min(diff.TSS.start.peak)) %>% 
  as.data.frame() %>% 
  select(gene.name) %>% 
  unique() %>%
  mutate(type.GR = "GR") -> gene.GR

# gene NOGR
enhancer_amplitude %>% 
  filter(gene.regulation != "random") %>%
  anti_join(., gene.GR, by = "gene.name") %>%
  select(gene.name) %>%
  unique() %>% 
  mutate(type.GR = "NUGR") -> gene.NUGR

# gene NUGR
max.change.time.point.significant %>% 
  ungroup() %>%
  anti_join(., enhancer_amplitude, by = "gene.name") %>%
  select(gene.name) %>%
  unique() %>% 
  mutate(type.GR = "NOGR") -> gene.NOGR

rbind(gene.GR, gene.NOGR, gene.NUGR) -> tmp.type.GR


col1 <- c("blue","red",  "green")

cutree(tmp_dend, k=2) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "gene.name") %>% 
  left_join(., tmp.type.GR, by = "gene.name") -> tmp_to_color

tmp_dend <- as.dist(1-cor(t(to.plot), method = "spearman")) %>% 
  hclust %>% 
  as.dendrogram %>% 
  color_branches(., k = number_clusters, col=c("firebrick", "dodgerblue"))


heatmap.2(to.plot,
          trace="none",
          margins =c(5,7),
          Colv = FALSE,
          col=bluered(20),
          dendrogram="row",      
          Rowv = tmp_dend,  
          key.xlab = "Concentration (index)",
          distfun = function(x) as.dist(1-cor(t(x)), method = "spearman"),
          RowSideColors = col1[as.factor(tmp_to_color$type.GR)],
          colRow = (get_leaves_branches_col(tmp_dend) %>% 
                       .[order(order.dendrogram(tmp_dend))])
) 

legend(y=1.1, x=.25, xpd=TRUE,      
       legend = unique(tmp_to_color$type.GR),
       col = unique(col1[as.factor(tmp_to_color$type.GR)]), 
       lty= 1,             
       lwd = 5,           
       cex=.7
)

rm(gene.GR, gene.NOGR, gene.NUGR, tmp_to_color, col1, tmp_dend, tmp_combine_MWT_MCTP)
rm(tmp_enhancer_amplitude, tmp.data.enhancer.range)
