#!/usr/bin/perl
package CQS::TestGroupTask;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::GroupTask;
use Test::More tests => 1;
use Data::Dumper;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_tt";
  bless $self, $class;
  return $self;
}

sub get_pbs_files {
  my ( $self, $config, $section, $pattern ) = @_;
  my $result = {
    G2 => "/temp/G2.pbs",
    G1 => "/temp/G1.pbs",
  };
  return($result);
}

my $test = CQS::TestGroupTask->new();

my $config = {
  "test" => {
    groups => {
      "G1" => ["S1_1", "S1_2"],
      "G2" => ["S2_1", "S2_2"],
    }
  }
};

my $source = $test->get_pbs_source($config, "test");

#print(Dumper(%$source));

is_deeply($source, 
	{"/temp/G1.pbs" => ["S1_1", "S1_2"],
     "/temp/G2.pbs" => ["S2_1", "S2_2"]});

1;


1
