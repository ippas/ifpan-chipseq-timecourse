#!/bin/bash

#bigwig_genomic_range_extract_normalize_to_tsv_bucket10.sh
#~/ifpan-chipseq-timecourse/DATA/enhancer_bigrange_delta_ep300_info.tsv


~/ifpan-chipseq-timecourse/SCRIPTS/./bigwig_genomic_range_extract_normalize_to_tsv_bucket10.sh ~/ifpan-chipseq-timecourse/DATA/enhancer_bigrange_top_ep300_info.tsv > ~/ChIP-seq/DATA/enhancer_bigrange_top_ep300.tsv


~/ifpan-chipseq-timecourse/SCRIPTS/./bigwig_genomic_range_extract_normalize_to_tsv_bucket10.sh ~/ifpan-chipseq-timecourse/DATA/enhancer_bigrange_delta_ep300_info.tsv > ~/ChIP-seq/DATA/enhancer_bigrange_delta_ep300.tsv

