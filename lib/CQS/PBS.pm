#!/usr/bin/perl
package CQS::PBS;

use strict;
use warnings;
use File::Basename;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(init_dir get_pbs_desc get_pbs_thread)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub init_dir {
  my ( $rootDir, $create ) = @_;

  if ( !defined $create ) {
    $create = 1;
  }

  #defined several folders
  my $pbsDir    = "$rootDir/pbs";
  my $resultDir = "$rootDir/result";
  my $logDir    = "$rootDir/log";

  if ($create) {
    create_directory_or_die($rootDir);
    create_directory_or_die($pbsDir);
    create_directory_or_die($resultDir);
    create_directory_or_die($logDir);
  }

  return ( $logDir, $pbsDir, $resultDir );
}

sub get_pbs_desc {
  my $walltime = "48";
  my $email    = "";
  my $mem      = "15000mb";
  my $nodes    = "1";

  my ($pbsParamHashRef) = @_;
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

  die "Assign email address in hash (\"email\" => \"youremail\") and pass hash as parameter to get_pbs_desc" if ( $email eq "" );

  my $pbsDesc = <<PBS;
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

  return ($pbsDesc);
}

sub get_pbs_thread {
  my $nodes = "1";

  my ($pbsParamHashRef) = @_;
  if ( defined $pbsParamHashRef ) {
    my %hash = %{$pbsParamHashRef};
    foreach ( keys %hash ) {
      if ( $_ eq "nodes" ) {
        $nodes = $hash{$_};
      }
    }
  }

  if ( $nodes =~ /1\:ppn=(\d+)/ ) {
    $nodes = $1;
  }

  return ($nodes);
}
1;

