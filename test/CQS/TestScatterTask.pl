#!/usr/bin/perl
package CQS::TestScatterTask;

use strict;
use warnings;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::StringUtils;
use CQS::ScatterTask;
use Data::Dumper;
use Test::More tests => 3;


our @ISA = qw(CQS::ScatterTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_st";
  bless $self, $class;
  return $self;
}

sub get_sample_names {
  return (["G1", "G2"]);
}

sub get_scatter_names {
  return (["chr1", "chr2", "chr3"]);
}

sub get_result_files {
  my ( $self, $config, $section, $result_dir, $sample_name, $scatter_name, $key_name ) = @_;
  return [$result_dir . "/" . $key_name . ".csv", $result_dir . "/" . $key_name . ".png"];
}

my $test = CQS::TestScatterTask->new();

my $config = {
  general => {
    task_name => "scatter",
  },
  "test" => {
    class => "CQS::TestScatterTask",
    target_dir => "/test",
    source     => {
      "G1" => [ "S1_1", "S1_2" ],
      "G2" => [ "S2_1", "S2_2" ],
    },
    pbs => {
      email => "a\@b.c",
    }
  },
};

#test result
my $section = "test";
my $source = $config->{$section}{source};
my $scatter_names = $test->get_scatter_names($config, $section);
my $key_sample_map = {};
for my $sample_name (sort keys %$source){
  for my $scatter_name (@$scatter_names){
    my $key = $sample_name . "." . $scatter_name;
    $key_sample_map->{$key} = [$sample_name];
  }
}

#test result
my $expect_result = {};
for my $key (sort keys %$key_sample_map){
  $expect_result->{$key} = ["/test/result/" . $key . ".csv", "/test/result/" . $key . ".png"];
}
my $actual_result = $test->result( $config, "test" );
is_deeply( $actual_result, $expect_result );

#test pbs
my $expect_pbs = {};
for my $key (sort keys %$key_sample_map){
  $expect_pbs->{$key} = "/test/pbs/" . $key . "_st.pbs";
}
my $actual_pbs = $test->get_pbs_files( $config, "test" );
is_deeply( $actual_pbs, $expect_pbs );

#test pbs source
my $expect_source = {};
for my $key (sort keys %$key_sample_map){
  $expect_source->{"/test/pbs/" . $key . "_st.pbs"} = $key_sample_map->{$key};
}
my $actual_source = $test->get_pbs_source( $config, "test" );
is_deeply( $actual_source, $expect_source );

1;