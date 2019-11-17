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
use Chipseq::MACS;
use Test::More tests => 1;
use Data::Dumper;

my $test = Chipseq::MACS->new();

my $config = {
	"general" => {
		task_name => "mytest",
	},
  "test" => {
  	target_dir => "/temp",
    groups => {
      "G1" => ["S1_1"],
      "G2" => ["S2_1"],
    },
    controls => {
      "G1" => ["S1_2"],
      "G2" => ["S2_2"],
    },
    
    pbs => {
    	email => "a\@c.d",
    },
  }
};

my $source = $test->get_pbs_source($config, "test");

print(Dumper(%$source));

is_deeply($source, 
	{"/temp/pbs/G1_macs.pbs" => ["S1_1", "S1_2"],
     "/temp/pbs/G2_macs.pbs" => ["S2_1", "S2_2"]});

1;


1
