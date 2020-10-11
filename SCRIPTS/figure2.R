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

## note, running this script requires variables created by script: figure1.R


# compare transcript lengths:


tmp_transcript_length <- as.data.frame(read.delim("DATA/transcript_length.tsv", 
                                                  header = TRUE, 
                                                  sep = "\t", 
                                                  stringsAsFactors = FALSE)) %>% 
  group_by(.$Gene.stable.ID) %>% 
  summarize(median = median(Transcript.length..including.UTRs.and.CDS.)) %>% 
  set_colnames(c("ensemblid", "median"))


results.filtered %>%
  left_join(., tmp_transcript_length, by = "ensemblid") %>%
  left_join(., gene_regulation, by = "gene.name") -> results.filtered


random %>%
  left_join(., tmp_transcript_length, by = "ensemblid") %>%  mutate(gene.regulation = "random") -> random


#create histogram of gene expression:

expression.hist <- rbind(results.filtered, sample_n(random, 320))
                 
ggplot(expression.hist, aes(x = log2(as.numeric(mean.expression)+1), fill = gene.regulation)) +     
  geom_histogram(alpha = 0.25, bins = 20, position = "identity")


#histogram of median transcript lenght

ggplot(expression.hist, aes(x = median, fill = gene.regulation)) +     
  geom_histogram(alpha = 0.25, bins = 50, position = "identity") + scale_x_continuous(trans='log2')


