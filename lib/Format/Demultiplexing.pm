#!/usr/bin/perl
package Format::Demultiplexing;

use strict;
use warnings;
use File::Basename;
use Data::Table;
use String::Util qw(trim);
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_dem";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my %mapFiles = %{ get_raw_files( $config, $section, "maps" ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $fastqfile    = $sample_files[0];
    my $summaryfile  = $sample_name . "_indecies.tsv";

    my @maps    = @{ $mapFiles{$sample_name} };
    my $mapfile = $maps[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );
    print $pbs "cqstools fastq_demultiplex $option -m $mapfile -i $fastqfile -o . -s $summaryfile";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->name() . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %mapFiles = %{ get_raw_files( $config, $section, "maps" ) };

  my $result = {};
  for my $sample_name ( keys %mapFiles ) {
    my @maps    = @{ $mapFiles{$sample_name} };
    my $mapfile = $maps[0];

    my $table = Data::Table::fromTSV( $mapfile, 0 );
    if ( $table->nofCol() == 3 ) {
      foreach my $i ( 0 .. $table->lastRow ) {
        my $name     = trim( $table->elm( $i, 2 ) );
        my $filename = trim( $table->elm( $i, 1 ) );
        $result->{$name} = [ $result_dir . "/" . $filename ];
      }
    }
    else {
      foreach my $i ( 0 .. $table->lastRow ) {
        my $filename = trim( $table->elm( $i, 1 ) );
        my $name = change_extension_gzipped($filename);
        $result->{$name} = [ $result_dir . "/" . $filename ];
      }
    }
  }
  return $result;
}

1;
