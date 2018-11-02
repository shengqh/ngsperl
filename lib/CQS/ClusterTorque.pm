#!/usr/bin/perl
package CQS::ClusterTorque;

use strict;
use warnings;
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
  my $walltime = "48";
  my $mem      = "15000mb";
  my $nodes    = "1";

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
      }
    }
  }

  die "Assign email address in hash (\"email\" => \"youremail\") and pass hash as parameter to get_cluster_desc" if ( $email eq "" );

  my $pbs_desc = <<PBS;
#!/bin/bash
#Beginning of PBS bash script
#PBS -M $email
#Status/Progress Emails to be sent
#PBS -m bae
#Email generated at b)eginning, a)bort, and e)nd of jobs
#PBS -l nodes=$nodes
#Processors needed
#PBS -l mem=$mem
#Total job memory required (specify how many megabytes)
#PBS -l walltime=${walltime}:00:00
#You must specify Wall Clock time (hh:mm:ss) [Maximum allowed 30 days = 720:00:00]
#PBS -q all
PBS

  return ($pbs_desc);
}

sub get_log_description {
  my ($self, $pbs_file) = @_;
  
  my $result = <<PBS;
#PBS -o $pbs_file
#PBS -j oe
PBS
  
  return ($result);
}

sub get_submit_command {
  return ("qsub");
}

1;
