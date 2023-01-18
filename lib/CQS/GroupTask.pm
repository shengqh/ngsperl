#!/usr/bin/perl
package CQS::GroupTask;

use strict;
use warnings;
use CQS::Task;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use Data::Dumper;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_pbskey} = "groups";
  $self->{_depend_all} = 0;
  $self->{_group_keys} = ["groups"],
  bless $self, $class;
  return $self;
}

sub get_pbs_key {
  my ($self, $config, $section) = @_;
  if ( has_raw_files( $config, $section, "groups" ) ) {
    return "groups";
  }else{
    return "source";
  }
}

sub get_pbs_source {
  my ( $self, $config, $section ) = @_;

  my $pbsFiles = $self->get_pbs_files( $config, $section );
  my $result = {};
  
  my $group_keys = $self->{"_group_keys"};
  for my $group_key (@$group_keys){
    if (has_raw_files($config, $section, $group_key)) {
  	  my $groups = get_raw_files($config, $section, $group_key);
	    for my $resKey ( keys %$pbsFiles ) {
	    	my $samples = $groups->{$resKey};
        if(!is_array($samples)){
          $samples = [$samples];
        }
	    	if(defined $result->{ $pbsFiles->{$resKey}}){
	    		my $oldSamples = $result->{ $pbsFiles->{$resKey}};
     		  push(@$oldSamples, @$samples);
	    	}else{
	    	  $result->{ $pbsFiles->{$resKey}} = $samples;
	    	}
	    }
    }
	}

	return $result;
}

1;
