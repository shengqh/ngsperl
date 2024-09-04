#!/usr/bin/perl
package CQS::ClassFactory;

use strict;
use warnings;
use Data::Dumper;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(instantiate getTaskClass performTask performTaskByPattern performConfig performTrace getResult)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.02';

sub instantiate {
  my $class = shift;
  my $location;
  if ( $class =~ m'::' ) {
    $location = $class;
    $location =~ s/::/\//g;
  }
  else {
    $location = "CQS/" . $class;
    $class    = "CQS::" . $class;
  }
  $location = $location . ".pm";

  #print "require $location \n";
  require $location;
  return $class->new(@_);
}

sub getTaskClass {
  my ( $config, $section ) = @_;
  if ( !defined $config->{$section} ) {
    die "No section $section defined.";
  }
  my $classname = $config->{$section}{class} or die "No class defined in section $section.";
  my $obj = instantiate($classname);
  $obj->{_config} = $config;
  $obj->{_section} = $section;
  return ($obj);
}

sub performTask {
  my ( $config, $section, $runImmediately ) = @_;
  my $obj = getTaskClass($config, $section);
  $obj->perform( $config, $section );
  if(defined $runImmediately && $runImmediately){
    my $pbs = $obj->get_pbs_files($config, $section);
    for my $pbskey (keys %$pbs){
      my $pbsfile = $pbs->{$pbskey};
      `bash $pbsfile 1 `;
    }
  }
}

sub performTaskByPattern {
  my ( $config, $section_pattern, $runImmediately ) = @_;
  for my $section (sort keys %$config){
    if ($section =~ /$section_pattern/){
      if (ref $config->{$section} eq ref {}){
        if (defined $config->{$section}{class}){
          print("Performing $section\n");
          performTask($config, $section, $runImmediately);
        }
      }
    }
  }
}

sub getResult {
  my ( $config, $section, $pattern ) = @_;
  my $obj = getTaskClass($config, $section);
  my $result = $obj->result( $config, $section, $pattern );
  return ($result);
}

sub performConfig {
  my ( $config, $pattern, $force ) = @_;

  my @sections = sort keys %$config;
  foreach my $section ( sort keys %$config ) {
    my $cursection = $config->{$section};
    if (ref $cursection ne ref {}){
      next;
    }
    if ( !defined $pattern || $section =~ m/$pattern/ ) {
      my $classname = $config->{$section}{class};
      if ( defined $classname ) {
        if ( $classname =~ /SequenceTask/ ) {
          next;
        }

        my $perform;
        if ( defined $force && $force ) {
          $perform = 1;
        }
        else {
          $perform = $config->{$section}{perform};
        }

        if ( !defined $perform || $perform ) {
          print "Performing $section by $classname \n";
          my $obj = instantiate($classname);
          $obj->{_config} = $config;
          $obj->{_section} = $section;
          $obj->perform( $config, $section );
        }
      }
    }
  }

  foreach my $section ( sort keys %$config ) {
    my $cursection = $config->{$section};
    if (ref $cursection ne ref {}){
      next;
    }
    if ( !defined $pattern || $section =~ m/$pattern/ ) {
      my $classname = $config->{$section}{class};
      if ( defined $classname ) {
        if ( $classname =~ /SequenceTask/ ) {

          my $perform;
          if ( defined $force && $force ) {
            $perform = 1;
          }
          else {
            $perform = $config->{$section}{perform};
          }

          if ( defined $perform && $perform ) {
            print "Performing $section by $classname \n";
            my $obj = instantiate($classname);
            $obj->{_config} = $config;
            $obj->{_section} = $section;
            $obj->perform( $config, $section );
          }
        }
      }
    }
  }
}

sub performTrace {
  my ( $config, $pattern, $force ) = @_;

  foreach my $section ( sort keys %{$config} ) {
    my $cursection = $config->{$section};
    foreach my $key ( keys %{$cursection} ) {
      if ( $key =~ /_ref$/ ) {
        my $refSectionName = $cursection->{$key};
        if ( is_array($refSectionName) ) {
          my @parts = @{$refSectionName};
          $refSectionName = $parts[0];
        }
        print "$section <- $refSectionName \n";
      }
    }
  }
}

1;
