#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use Test::Simple;

my $data = read_cufflinks_fpkm("../../data/genes.fpkm_tracking");

ok (scalar(keys %{$data}) eq 12);
ok ($data->{ENSG00000186092} eq "0");
ok ($data->{ENSG00000236679} eq "558215");
ok ($data->{"CUFF.8"} eq "3.14195e+06");
