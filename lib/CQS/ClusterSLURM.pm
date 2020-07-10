#!/usr/bin/perl
package CQS::ClusterSLURM;

use strict;
use warnings;
use CQS::ClusterScript;

our @ISA = qw(CQS::ClusterScript);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = __PACKAGE__;
  bless $self, $class;
  return $self;
}

sub get_cluster_desc {
  my $walltime = "48";
  my $mem      = "15000M";
  my $nodes    = "1";
  my $ntasks   = "";

  my ( $self, $pbsParamHashRef, $config ) = @_;
  my $generalOptions = $self->get_general_options($config);
  my $email          = $generalOptions->{email};
  my $emailType      = $generalOptions->{emailType};
  my $constraint     = $generalOptions->{constraint};
  my $account        = $generalOptions->{account};

  if ( defined $pbsParamHashRef ) {
    my %hash = %{$pbsParamHashRef};
    foreach ( keys %hash ) {
      if ( $_ eq "walltime" ) {
        $walltime = $hash{$_};
      }
      elsif ( $_ eq "email" ) {
        $email = $hash{$_};
      }
      elsif ( $_ eq "emailType" ) {
        if ( defined $hash{$_} ) {
          $emailType = $hash{$_};
        }
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
#SBATCH --mail-type=$emailType
#SBATCH --nodes=$nodes
#SBATCH --ntasks=$ntasks
#SBATCH --time=$walltime
#SBATCH --mem=$mem
SBATCH
  if ( defined $constraint ) {
    $pbs_desc = $pbs_desc . "#SBATCH --constraint=$constraint\n";
  }
  if ( defined $account ) {
    $pbs_desc = $pbs_desc . "#SBATCH --account=$account\n";
  }

  return ($pbs_desc);
}

sub get_log_description {
  my ( $self, $pbs_file ) = @_;

  my $result = <<SBATCH;
#SBATCH -o $pbs_file

smemwatch -k 99 -d 50 \$\$ \& 
SBATCH

  return ($result);
}

sub get_submit_command {
  return ("sbatch");
}

1;
