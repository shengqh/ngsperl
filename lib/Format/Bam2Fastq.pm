#!/usr/bin/perl
package Format::Bam2Fastq;

use strict;
use warnings;
use File::Basename;
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
  $self->{_suffix} = "_b2q";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $ispaired = get_option( $config, $section, "ispaired", 0 );
  my $unzipped = get_option( $config, $section, "unzipped", 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $final_file = $ispaired ? $sample_name . "_1.fastq" : $sample_name . ".fastq";
    if ( !$unzipped ) {
      $final_file = $final_file . ".gz";
    }

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

    print $pbs "
if [ ! -s $final_file ]; then
  bam2fastq -o ${sample_name}#.fastq $bam_file
";
    if ( !$unzipped ) {
      if ($ispaired) {
        print $pbs "  gzip ${sample_name}_1.fastq
  gzip ${sample_name}_2.fastq ";
      }
      else {
        print $pbs "  gzip ${sample_name}.fastq";
      }
    }
    print $pbs "
fi
";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all Bam2Fastq tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $ispaired = get_option( $config, $section, "ispaired", 0 );
  my $unzipped = get_option( $config, $section, "unzipped", 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {

    my @result_files = ();
    if ($ispaired) {
      if ($unzipped) {
        push( @result_files, $result_dir . "/" . $sample_name . "_1.fastq" );
        push( @result_files, $result_dir . "/" . $sample_name . "_2.fastq" );
      }
      else {
        push( @result_files, $result_dir . "/" . $sample_name . "_1.fastq.gz" );
        push( @result_files, $result_dir . "/" . $sample_name . "_2.fastq.gz" );
      }
    }
    else {
      if ($unzipped) {
        push( @result_files, $result_dir . "/" . $sample_name . ".fastq" );
      }
      else {
        push( @result_files, $result_dir . "/" . $sample_name . ".fastq.gz" );
      }
    }

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
