#!/usr/bin/perl
package CQS::TestGatherTask;

use strict;
use warnings;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::StringUtils;
use CQS::GatherTask;
use Data::Dumper;
use Test::More tests => 4;

our @ISA = qw(CQS::GatherTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_gt";
  bless $self, $class;
  return $self;
}

sub get_gather_map {
  return {
    "G1" => ["G1_chr1", "G1_chr2"],
    "G2" => ["G2_chr1", "G2_chr2"],
  };
}

sub get_result_files {
  my ( $self, $config, $section, $result_dir, $gather_name ) = @_;
  return [$result_dir . "/" . $gather_name . ".csv", $result_dir . "/" . $gather_name . ".png"];
}

my $test = CQS::TestGatherTask->new();

my $config = {
  general => {
    task_name => "gather",
  },
  "test" => {
    class => "CQS::TestGatherTask",
    target_dir => "/test",
    source     => {
      "G1_chr1" => [ "G1_chr1.csv" ],
      "G1_chr2" => [ "G1_chr2.csv" ],
      "G2_chr1" => [ "G2_chr1.csv" ],
      "G2_chr2" => [ "G2_chr2.csv" ],
    },
    pbs => {
      email => "a\@b.c",
    }
  },
};

#test result
my $section = "test";
my $gather_map = $test->get_gather_map($config, $section);
my $expect_result = {};
for my $key (sort keys %$gather_map){
  $expect_result->{$key} = ["/test/result/" . $key . ".csv", "/test/result/" . $key . ".png"];
}
my $actual_result = $test->result( $config, "test" );
is_deeply( $actual_result, $expect_result );

#test pbs
my $expect_pbs = {};
for my $key (sort keys %$gather_map){
  $expect_pbs->{$key} = "/test/pbs/" . $key . "_gt.pbs";
}
my $actual_pbs = $test->get_pbs_files( $config, "test" );
is_deeply( $actual_pbs, $expect_pbs );

#test pbs source
my $expect_source = {};
for my $key (sort keys %$gather_map){
  $expect_source->{"/test/pbs/" . $key . "_gt.pbs"} = $gather_map->{$key};
}
my $actual_source = $test->get_pbs_source( $config, "test" );
is_deeply( $actual_source, $expect_source );

#test result_pbs
my $expect_result_pbs_map = {};
for my $key (sort keys %$gather_map){
  $expect_result_pbs_map->{$key} = "/test/pbs/" . $key . "_gt.pbs";
}

my $result_pbs_map = $test->get_result_pbs($config, "test");
is_deeply( $result_pbs_map, $expect_result_pbs_map );

1;