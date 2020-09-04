#!/usr/bin/perl
package Proteomics::Summary::AccumulateBuildSummary;

use strict;
use warnings;
use File::Basename;
use File::Slurp;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;
use Proteomics::Summary::AbstractBuildSummary;

our @ISA = qw(Proteomics::Summary::AbstractBuildSummary);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_abs";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $proteomicstools = get_param_file( $config->{$section}{proteomicstools}, "proteomicstools", 1, not $self->using_docker() );

  my $bin_size = get_option( $config, $section, "bin_size", 10 );
  my $bins     = $config->{$section}{"bins"};
  my @bins     = ();
  if ( defined $bins ) {
    @bins = @{$bins};
  }

  my ( $datasets, $lines, $dataset ) = $self->get_datasets( $config, $section );
  my %datasets = %{$datasets};
  my @lines    = @{$lines};
  my @dataset  = @{$dataset};

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  my @sample_names     = sort keys %datasets;
  my $sample_nameCount = scalar(@sample_names);
  my $numlen           = length($sample_nameCount);

  my $index      = 0;
  my $sampleSize = 0;

  while ( $sampleSize <= $sample_nameCount ) {
    if ( $index < scalar(@bins) ) {
      $sampleSize = $bins[$index];
      $index++;
    }
    else {
      $sampleSize += $bin_size;
    }
    if ( $sampleSize > $sample_nameCount ) {
      $sampleSize = $sample_nameCount;
    }
    my $currentTaskName = sprintf( "${task_name}_%0${numlen}d", $sampleSize );

    my $currentParamFile = $result_dir . "/" . $currentTaskName . ".param";
    open( my $out1, ">$currentParamFile" ) or die $!;

    for ( my $index = 0 ; $index < scalar(@lines) ; $index++ ) {
      print $out1 $lines[$index] . "\n";
      if ( $lines[$index] =~ "<Datasets>" ) {
        last;
      }
    }
    for ( my $index = 0 ; $index < $sampleSize ; $index++ ) {
      my $sample_name = $sample_names[$index];
      print $out1 "    <Dataset>\n";
      print $out1 "      <Name>$sample_name</Name>\n";
      foreach my $dsline (@dataset) {
        print $out1 $dsline . "\n";
      }
      print $out1 "      <PathNames>\n";
      my @sample_files = @{ $datasets{$sample_name} };
      for my $sampleFile (@sample_files) {
        print $out1 "        <PathName>$sampleFile</PathName>\n";
      }
      print $out1 "      </PathNames>\n";
      print $out1 "    </Dataset>\n";
    }

    for ( my $index = 0 ; $index < scalar(@lines) ; $index++ ) {
      if ( $lines[$index] =~ "</Datasets>" ) {
        for ( my $nextindex = $index ; $nextindex < scalar(@lines) ; $nextindex++ ) {
          print $out1 $lines[$nextindex] . "\n";
        }
        last;
      }
    }

    close $out1;

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $currentTaskName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $currentTaskName );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $result_file = $result_dir . "/" . $currentTaskName . ".noredundant";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $result_file );
    print $pbs "mono $proteomicstools buildsummary -i $currentParamFile";
    $self->close_pbs($pbs);

    print "$pbs_file created \n";

    if ( $sampleSize == scalar(@sample_names) ) {
      last;
    }
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all ", $self->{_name}, " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $bin_size = get_option( $config, $section, "bin_size", 10 );
  my $bins     = $config->{$section}{"bins"};
  my @bins     = ();
  if ( defined $bins ) {
    @bins = @{$bins};
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result       = {};
  my @result_files = ();

  my $sample_nameCount = scalar( keys %raw_files );
  my $numlen           = length($sample_nameCount);

  my $index      = 0;
  my $sampleSize = 0;

  while ( $sampleSize <= $sample_nameCount ) {
    if ( $index < scalar(@bins) ) {
      $sampleSize = $bins[$index];
      $index++;
    }
    else {
      $sampleSize += $bin_size;
    }
    if ( $sampleSize > $sample_nameCount ) {
      $sampleSize = $sample_nameCount;
    }
    my $currentTaskName = sprintf( "${task_name}_%0${numlen}d", $sampleSize );

    my $currentNoredundantFile = $result_dir . "/" . $currentTaskName . ".noredundant";
    push( @result_files, $currentNoredundantFile );

    if ( $sampleSize == $sample_nameCount ) {
      last;
    }

    $sampleSize += $bin_size;
    if ( $sampleSize > $sample_nameCount ) {
      $sampleSize = $sample_nameCount;
    }
  }
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}
1;
