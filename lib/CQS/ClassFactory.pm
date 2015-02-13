#!/usr/bin/perl
package CQS::ClassFactory;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(instantiate performTask performConfig performTrace getPbsFiles)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.02';

sub instantiate {
  my $class = shift;
  my $location;
  if ( $class =~ m'::' ) {
    $location = $class;
    $location =~ s/::/\//;
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

sub performTask {
  my ( $config, $section ) = @_;
  if ( !defined $config->{$section} ) {
    die "No section $section defined.";
  }
  my $classname = $config->{$section}{class} or die "No class defined in section $section.";
  my $obj = instantiate($classname);
  $obj->perform( $config, $section );
}

sub getPbsFiles {
  my ( $config, $section ) = @_;
  if ( !defined $config->{$section} ) {
    die "No section $section defined.";
  }
  my $classname = $config->{$section}{class} or die "No class defined in section $section.";
  my $obj = instantiate($classname);
  $obj->pbsfiles( $config, $section );
}

sub performConfig {
  my ( $config, $pattern, $force ) = @_;
  foreach my $section ( keys %{$config} ) {
    if ( !defined $pattern || $section =~ m/$pattern/ ) {
      my $classname = $config->{$section}{class};
      if ( defined $classname ) {
        if ( $classname eq "CQS::SequenceTask" || $classname eq "SequenceTask" ) {
          next;
        }

        my $perform;
        if ( defined $force && $force ) {
          $perform = 1;
        }
        else {
          $perform = $config->{$section}{perform};
        }

        if ( defined $perform && $perform ) {
          print "Preforming $section by $classname \n";
          my $obj = instantiate($classname);
          $obj->perform( $config, $section );
        }
      }
    }
  }

  foreach my $section ( keys %{$config} ) {
    if ( !defined $pattern || $section =~ m/$pattern/ ) {
      my $classname = $config->{$section}{class};
      if ( defined $classname ) {
        if ( $classname eq "CQS::SequenceTask" || $classname eq "SequenceTask" ) {

          my $perform;
          if ( defined $force && $force ) {
            $perform = 1;
          }
          else {
            $perform = $config->{$section}{perform};
          }

          if ( defined $perform && $perform ) {
            print "Preforming $section by $classname \n";
            my $obj = instantiate($classname);
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
        if ( ref($refSectionName) eq 'ARRAY' ) {
          my @parts = @{$refSectionName};
          $refSectionName = $parts[0];
        }
        print "$section <- $refSectionName \n";
      }
    }
  }
}

1;
