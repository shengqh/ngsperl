#!/usr/bin/perl
package CQS::ClusterSLURM;

use strict;
use warnings;
use CQS::ClusterScript;

our @ISA = qw(CQS::ClusterScript);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "ClusterSLURM";
  bless $self, $class;
  return $self;
}

sub get_cluster_desc {
  my $walltime = "48";
  my $email    = "";
  my $mem      = "15000M";
  my $nodes    = "1";
  my $ntasks   = "";

  my ($pbsParamHashRef, $pbsfile) = @_;
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
      elsif ( $_ eq "ntasks" ) {
        $ntasks = $hash{$_};
      }
    }
  }
  
  die "Assign email address in hash (\"email\" => \"youremail\") and pass hash as parameter to get_cluster_desc" if ( $email eq "" );

  if ( $ntasks eq "" ) {
  	$ntasks = $nodes;
  	$nodes = "1";
  }
  
  $mem =~ s/mb/M/g;

  my $pbsDesc = <<SBATCH;
#!/bin/bash
#SBATCH --mail-user=$email
#SBATCH --mail-type=ALL
#SBATCH --nodes=$nodes
#SBATCH --ntasks=$ntasks
#SBATCH --time=$walltime
#SBATCH --mem=$mem
#SBATCH -o $pbsfile
SBATCH

  return ($pbsDesc);
}


sub get_log_desc {
  my ($pbsfile) = @_;
  
  my $result = <<SBATCH;
#SBATCH -o $pbsfile
SBATCH
  
  return ($result);
}

sub get_submit_command {
  return ("sbatch");
}

1;
