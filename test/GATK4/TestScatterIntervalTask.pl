#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use Data::Dumper;

my $config = {
  "test1" => {
    class       => "GATK4::ScatterIntervalsTask",
    interval_list_file => "/scratch/jbrown_lab/references/genome_AAV1_AAV2/genome_AAV1_AAV2.intervals_list"
  },
  "test2" => {
    class       => "GATK4::HaplotypeCallerScatter",
    interval_list_file => "/scratch/jbrown_lab/references/genome_AAV1_AAV2/genome_AAV1_AAV2.intervals_list"
  },
};

my $myclass    = instantiate( $config->{test2}{class} );
#print("_docker_prefix=" . $myclass->{_docker_prefix});

#my $myclass2    = instantiate( $config->{test2}{class} );
#print("_docker_prefix=" . $myclass2->{_docker_prefix});


#my $map = $myclass->get_interval_file_map($config, "test");

#print(Dumper($map));
$myclass->get_docker_value();

1;
