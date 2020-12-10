#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Test::More tests => 3;
use Pipeline::PipelineUtils;

my $chromosomes = readChromosomeFromDictFile(dirname(__FILE__) . "/../../data/GRCm38.p6.genome.dict");
is (scalar(@$chromosomes), 22);
is (join(",", @$chromosomes), 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX,chrY,chrM');

my $config = {};
addFilesFromSraRunTable($config, dirname(__FILE__) . "/../../data/SraRunTable.txt");

is_deeply($config->{files}, {
                       'HGADFN166_ABE10day_1' => [
                                                   'SRR11657427'
                                                 ],
                       'Untreated_Parent167_1' => [
                                                    'SRR11657417'
                                                  ],
                       'HGADFN166_ABE20day_3' => [
                                                   'SRR11657408'
                                                 ],
                       'HGADFN188_ABE10day_3' => [
                                                   'SRR11657405'
                                                 ],
                       'HGADFN188_ABE20day_1' => [
                                                   'SRR11657404'
                                                 ],
                       'HGADFN188_ABE20day_2' => [
                                                   'SRR11657425'
                                                 ],
                       'Untreated_HGADFN188_3' => [
                                                    'SRR11657418'
                                                  ],
                       'Untreated_HGADFN166_3' => [
                                                    'SRR11657421'
                                                  ],
                       'Untreated_Parent167_2' => [
                                                    'SRR11657416'
                                                  ],
                       'HGADFN166_ABE10day_2' => [
                                                   'SRR11657426'
                                                 ],
                       'HGADFN166_ABE20day_2' => [
                                                   'SRR11657409'
                                                 ],
                       'HGADFN188_ABE10day_2' => [
                                                   'SRR11657406'
                                                 ],
                       'Untreated_HGADFN188_1' => [
                                                    'SRR11657420'
                                                  ],
                       'Untreated_HGADFN166_1' => [
                                                    'SRR11657423'
                                                  ],
                       'HGADFN188_ABE20day_3' => [
                                                   'SRR11657424'
                                                 ],
                       'HGADFN188_ABE10day_1' => [
                                                   'SRR11657407'
                                                 ],
                       'Untreated_HGADFN166_2' => [
                                                    'SRR11657422'
                                                  ],
                       'Untreated_HGADFN188_2' => [
                                                    'SRR11657419'
                                                  ],
                       'HGADFN166_ABE20day_1' => [
                                                   'SRR11657410'
                                                 ],
                       'HGADFN166_ABE10day_3' => [
                                                   'SRR11657415'
                                                 ],
                       'Untreated_Parent167_3' => [
                                                    'SRR11657414'
                                                  ]
                     }
        );

1;
