###################################################
# Prepare +/-10000 for enhancer without promoters #
###################################################



####################################################
# centering peaks for top ep300 without promotores #
####################################################

#aligned to the highest value

results.list <- lapply(readLines('/home/mateusz/ChIP-seq/DATA/enhancer_peaks_top_ep300.tsv'), strsplit, split = '\t')
lapply(results.list, lengths) -> results.lenghts

number.column <- results.list %>% lapply(lengths) %>% unlist() %>% max() # result = number.column - so length of table: (number.column x 2) + 1 = number.column x 2 + 1 so the middle value is number.column - this is where we put max

adjust_numeric_columns <- function(x) {
  
  new.line <- integer(number.column * 2) # make a vector of zeros
  x %>% unlist() -> single.row # unlist the line
  single.row[c(9:length(single.row))] -> peak #get just the numeric valuse
  peak %>% as.numeric() %>% which.max() -> amp.index #find where the max value of the peak is
  new.line[c((number.column - amp.index + 1): (number.column + (length(peak) - amp.index)))] <- peak #put the peak into the new.line with the max set to index 1144
  full.line <- c(single.row[c(1:8)], new.line)
  return(full.line)
}

results.list %>% lapply(adjust_numeric_columns) %>% data.frame(row.names = NULL) -> adjusted.peaks.table

names(adjusted.peaks.table) <- NULL

adjusted.peaks.table %>% t() %>% data.frame() -> adjusted.peaks.table

adjusted.peaks.table %>%
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", "position.amplitude", 1:(number.column*2-1))) %>%  
gather(., "range", "value", -c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file")) -> gather.adjusted.peaks.table 


#line plot for enhancer peak without promoters
gather.adjusted.peaks.table %>% 
  mutate(range = as.numeric(range), value = as.numeric(value)) %>%
  group_by(gene.name, gene.regulation, TF, time, range) %>%
  summarise(value = mean(value)) %>%
  ungroup %>%
  group_by(gene.name, gene.regulation, TF, time, range) %>%
  summarize(value = mean(value)) %>%
  group_by(gene.regulation, TF, time, range) %>%
  summarize(value = mean(value)) -> tmp.to.plot
tmp.to.plot %>% 
  ungroup() %>%
  mutate(time = as.numeric(as.character(time))) %>% 
{ggplot(., aes(x = range, y = value, color = gene.regulation)) +
    geom_line() +
    # scale_color_manual(values = c("up" = "firebrick",
    #                               "random" = "gray",
    #                               "down" = "dodgerblue")) +
    facet_grid(TF~time)}


#########################################################################################################


# gather.adjusted.peaks.table %>% 
#   mutate(range = as.numeric(range), value = as.numeric(value)) %>%
#   group_by(gene.name, gene.regulation, TF, time, range) %>%
#   summarise(value = mean(value)) %>%
#   ungroup %>%
#   #filter(TF=="NR3C1") %>% 
#   mutate(time = as.numeric(as.character(time))) %>% 
#   {ggplot(., aes(x = range, y = value)) +
#       
#       geom_tile(aes(x=range, y=reorder(gene.name, value), fill=value)) +
#       # scale_color_manual(values = c("up" = "firebrick",
#       #                               "random" = "gray",
#       #                               "down" = "dodgerblue")) +
#       facet_grid(TF~time)}
#   

#prepare +/-10000 range for amplitude enhancer NR3C1 and save to file
read.table("~/ChIP-seq/DATA/enhancer_peaks_top_ep300.tsv",
           header = FALSE,
           sep = "\t",
           col.names = c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", seq_len((number.column-8)*2)),
           stringsAsFactors = FALSE,
           fill=TRUE) %>%
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
fwrite("~/ifpan-chipseq-timecourse/DATA/enhancer_bigrange_top_ep300_info.tsv", 
       sep="\t", 
       col.names = TRUE, 
       row.names = FALSE)



######################################################
# centering peaks for delta ep300 without promotores #
######################################################

#aligned to the highest value

results.list <- lapply(readLines('/home/mateusz/ChIP-seq/DATA/enhancer_peaks_delta_ep300.tsv'), strsplit, split = '\t')
lapply(results.list, lengths) -> results.lenghts

number.column <- results.list %>% lapply(lengths) %>% unlist() %>% max() # result = number.column - so length of table: (number.column x 2) + 1 = number.column x 2 + 1 so the middle value is number.column - this is where we put max

adjust_numeric_columns <- function(x) {
  
  new.line <- integer(number.column * 2) # make a vector of zeros
  x %>% unlist() -> single.row # unlist the line
  single.row[c(9:length(single.row))] -> peak #get just the numeric valuse
  peak %>% as.numeric() %>% which.max() -> amp.index #find where the max value of the peak is
  new.line[c((number.column - amp.index + 1): (number.column + (length(peak) - amp.index)))] <- peak #put the peak into the new.line with the max set to index 1144
  full.line <- c(single.row[c(1:8)], new.line)
  return(full.line)
}

results.list %>% lapply(adjust_numeric_columns) %>% data.frame(row.names = NULL) -> adjusted.peaks.table

names(adjusted.peaks.table) <- NULL

adjusted.peaks.table %>% t() %>% data.frame() -> adjusted.peaks.table

adjusted.peaks.table %>%
  set_colnames(c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", "position.amplitude", 1:(number.column*2-1))) %>%  
  gather(., "range", "value", -c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file")) -> gather.adjusted.peaks.table 


#line plot for enhancer peak without promoters
gather.adjusted.peaks.table %>% 
  mutate(range = as.numeric(range), value = as.numeric(value)) %>%
  group_by(gene.name, gene.regulation, TF, time, range) %>%
  summarise(value = mean(value)) %>%
  ungroup %>%
  group_by(gene.name, gene.regulation, TF, time, range) %>%
  summarize(value = mean(value)) %>%
  group_by(gene.regulation, TF, time, range) %>%
  summarize(value = mean(value)) -> tmp.to.plot
tmp.to.plot %>% 
  ungroup() %>%
  mutate(time = as.numeric(as.character(time))) %>% 
  {ggplot(., aes(x = range, y = value, color = gene.regulation)) +
      geom_line() +
      # scale_color_manual(values = c("up" = "firebrick",
      #                               "random" = "gray",
      #                               "down" = "dodgerblue")) +
      facet_grid(TF~time)}


#########################################################################################################



#prepare +/-10000 range for amplitude enhancer NR3C1 and save to file
read.table("~/ChIP-seq/DATA/enhancer_peaks_delta_ep300.tsv",
           header = FALSE,
           sep = "\t",
           col.names = c("gene.name", "chromosome", "start.range", "end.range", "gene.regulation", "TF", "time", "file", seq_len((number.column-8)*2)),
           stringsAsFactors = FALSE,
           fill=TRUE) %>%
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
  fwrite("~/ifpan-chipseq-timecourse/DATA/enhancer_bigrange_delta_ep300_info.tsv", 
         sep="\t", 
         col.names = TRUE, 
         row.names = FALSE)

