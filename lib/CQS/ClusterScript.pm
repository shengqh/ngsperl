#!/usr/bin/perl
package CQS::ClusterScript;

use strict;
use warnings;
use CQS::ConfigUtils;

sub new {
  my ($class) = @_;
  my $self = { _name => __PACKAGE__ };
  bless $self, $class;
  return $self;
}

sub name {
  my ($self) = @_;
  return $self->{_name};
}

sub get_cluster_thread {
  my $nodes    = "1";
  my $ntasks   = "";

  my ( $self, $pbsParamHashRef ) = @_;
  if ( defined $pbsParamHashRef ) {
    my %hash = %{$pbsParamHashRef};
    foreach ( keys %hash ) {
      if ( $_ eq "nodes" ) {
        $nodes = $hash{$_};
      }
      elsif ( $_ eq "ntasks" ) {
        $ntasks = $hash{$_};
      }
    }
  }

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
  
  return ($ntasks);
}

sub get_cluster_memory {
  my $result    = "10G";

  my ( $self, $pbsParamHashRef ) = @_;
  if ( defined $pbsParamHashRef ) {
    my %hash = %{$pbsParamHashRef};
    foreach ( keys %hash ) {
      if ( $_ eq "mem" ) {
        $result = $hash{$_};
      }
    }
  }
  
  $result =~ s/mb/M/g;
  $result =~ s/gb/G/g;
  
  return ($result);
}

sub get_cluster_desc {
}

sub get_log_desc {
}

sub get_submit_command {
}

1;
