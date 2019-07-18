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
use Test::More tests => 4;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_tt";
  bless $self, $class;
  return $self;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my $result = {
    S3 => ["/temp/S3.tsv", "/temp/S3.csv"],
    S2 => ["/temp/S2.tsv", "/temp/S2.csv"],
    S1 => ["/temp/S1.tsv", "/temp/S1.csv"],
  };
  return($result);
}

my $test = CQS::TestTask->new();

is($test->get_final_file({}, "", "/temp"), "S3.csv");
is($test->get_final_file({}, "", "/temp/"), "S3.csv");
is($test->get_final_file({}, "", "/temp", "S1"), "S1.csv");
is($test->get_final_file({}, "", "/temp/", "S1"), "S1.csv");

1;


1