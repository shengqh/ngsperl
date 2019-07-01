#!/usr/bin/perl
package CQS::ClusterTorqueSingapore;

use strict;
use warnings;
use File::Basename;
use CQS::ClusterScript;

our @ISA = qw(CQS::ClusterScript);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  bless $self, $class;
  return $self;
}

sub get_cluster_desc {
  my $walltime = "24";
  my $mem      = "40G";
  my $nodes    = "10";
  my $project = "21070175";

  my ($self, $pbsParamHashRef, $config) = @_;
  my $generalOptions = $self->get_general_options($config);
  my $email = $generalOptions->{email};
  
  if ( defined $pbsParamHashRef ) {
    my %hash = %{$pbsParamHashRef};
    foreach ( keys %hash ) {
      if ( $_ eq "walltime" ) {
        $walltime = $hash{$_};
      }
      elsif ( $_ eq "email" ) {
        $email = $hash{$_};
      }
      elsif ( $_ eq "mem" ) {
        $mem = $hash{$_};
      }
      elsif ( $_ eq "nodes" ) {
        $nodes = $hash{$_};
      }elsif ($_ eq "project"){
        $project = $hash{$_};
      }
    }
  }
  
  if ($nodes =~ /ppn/){
    $nodes =~ s/1:ppn=//g;
  }
  
  my $newmem = $mem;
  $newmem =~ s/\D//g;
  my $newnodes = $newmem / 4;
  if($newnodes > $nodes){
    $nodes = $newnodes;
  }

  die "Assign email address in hash (\"email\" => \"youremail\") and pass hash as parameter to get_cluster_desc" if ( $email eq "" );

  my $pbs_desc = <<PBS;
#!/bin/bash
#PBS -q normal
#PBS -l select=1:ncpus=$nodes:mem=$mem
#PBS -l walltime=${walltime}:00:00
#PBS -P $project
PBS

  return ($pbs_desc);
}

sub get_log_description {
  my ($self, $log_file) = @_;
  
  my ($file,$dir,$ext) = fileparse($log_file, qr/\.[^.]*/);
  my $result = <<PBS;
#PBS -o $log_file
#PBS -j oe
#PBS -N $file
PBS
  
  return ($result);
}

sub get_submit_command {
  return ("qsub");
}

1;
