#!/usr/bin/perl
package CQS::ClusterSLURM;

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
  my $email    = "";
  my $mem      = "15000M";
  my $nodes    = "1";
  my $ntasks   = "";

  my ( $self, $pbsParamHashRef ) = @_;
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
    my $pos = index( $nodes, ":" );
    if ( $pos >= 0 ) {
      $ntasks = substr( $nodes, $pos + 1 );
      $nodes = substr( $nodes, 0, $pos );
      $pos = index( $ntasks, "ppn=" );
      if ( $pos >= 0 ) {
        $ntasks = substr( $ntasks, $pos + 4 );
      }
    }
    else {
      $nodes  = "1";
      $ntasks = $nodes;
    }
  }

  if ( index( $walltime, ':' ) == -1 ) {
    $walltime = $walltime . ":00:00";
  }

  $mem =~ s/mb/M/g;
  $mem =~ s/gb/G/g;

  my $pbs_desc = <<SBATCH;
#!/bin/bash
#SBATCH --mail-user=$email
#SBATCH --mail-type=ALL
#SBATCH --nodes=$nodes
#SBATCH --ntasks=$ntasks
#SBATCH --time=$walltime
#SBATCH --mem=$mem
SBATCH

  return ($pbs_desc);
}

sub get_log_description {
  my ( $self, $pbs_file ) = @_;

  my $result = <<SBATCH;
#SBATCH -o $pbs_file
SBATCH

  return ($result);
}

sub get_submit_command {
  return ("sbatch");
}

1;
