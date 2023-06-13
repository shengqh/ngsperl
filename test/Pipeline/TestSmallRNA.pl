#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Test::More tests => 1;
use CQS::StringUtils;

{
  my $result_files = ["/gpfs52/data/stein_lab/mjo_sRNA_data/20230613_10065_AP_smRNAseq_human_truseq/host_genome/deseq2_miRNA_TotalReads/result/Jo_vs_Healthy_min5_pvalue0.05_DESeq2_volcanoEnhanced.png",
    "/gpfs52/data/stein_lab/mjo_sRNA_data/20230613_10065_AP_smRNAseq_human_truseq/host_genome/deseq2_miRNA_TotalReads/result/nonJo_vs_Healthy_min5_pvalue0.05_DESeq2_volcanoEnhanced.png"];
  
  my $comparison = "Jo_vs_Healthy";
  my $pattern = "[\/]${comparison}_.+_volcanoEnhanced.png";

  my $filtered = filter_array( $result_files, $pattern, 1 );

  is_deeply([$result_files->[0]], $filtered);
}


1;
