#!/usr/bin/perl
package CQS::TestTask;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::Task;
use Test::More tests => 3;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_t";
  bless $self, $class;
  return $self;
}

sub get_result_files {
  my ( $self, $config, $section, $result_dir, $sample_name ) = @_;
  return [$result_dir . "/" . $sample_name . ".csv", $result_dir . "/" . $sample_name . ".png"];
}

my $config = {
  general => {
    task_name => "task",
  },
  "test" => {
    class => "CQS::Task",
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

my $test = CQS::TestTask->new();

#test result
my $section = "test";
my $source = $config->{$section}{source};
my $key_sample_map = {};
for my $sample_name (sort keys %$source){
  $key_sample_map->{$sample_name} = [$sample_name];
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
  $expect_pbs->{$key} = "/test/pbs/" . $key . "_t.pbs";
}
my $actual_pbs = $test->get_pbs_files( $config, "test" );
is_deeply( $actual_pbs, $expect_pbs );

#test pbs source
my $expect_source = {};
for my $key (sort keys %$key_sample_map){
  $expect_source->{"/test/pbs/" . $key . "_t.pbs"} = $key_sample_map->{$key};
}
my $actual_source = $test->get_pbs_source( $config, "test" );
is_deeply( $actual_source, $expect_source );

1