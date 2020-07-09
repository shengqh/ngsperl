#!/usr/bin/perl
package GATK4::SamToFastq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use GATK4::GATK4Task;

our @ISA = qw(GATK4::GATK4Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_sf";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

  my $java_option = $self->get_java_option( $config, $section, $memory );

  #parameter files
#  $self->get_docker_value(1);

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $interleave_fastq = get_option( $config, $section, "interleave_fastq", 1 );
#  $interleave_fastq = ( $interleave_fastq ) ? "true" : "false";

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile   = $sample_files[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $unmapped_fastq_file =$sample_name . ".fastq.gz";
    my $unmapped_fastq_file1 =$sample_name . ".1.fastq.gz";
    my $unmapped_fastq_file2 =$sample_name . ".2.fastq.gz";
    my $final_file = ( $interleave_fastq ) ? $unmapped_fastq_file : $unmapped_fastq_file1;

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir,$final_file, $init_command );

    my $pbsCmd="";
    if ($interleave_fastq) { #paired together
        $pbsCmd="gatk --java-options \"$java_option\" SamToFastq -I $sampleFile -F ${unmapped_fastq_file} -INTERLEAVE true"
    } else { #paired seprately
        $pbsCmd="gatk --java-options \"$java_option\" SamToFastq -I $sampleFile -F ${unmapped_fastq_file1} -F2 ${unmapped_fastq_file2} -INTERLEAVE false"
    }
    print $pbs "  

cd $result_dir

${pbsCmd}
";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $interleave_fastq = get_option( $config, $section, "interleave_fastq", 1 );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    if ($interleave_fastq) { #paired together
        my $unmapped_fastq_file =$sample_name . ".fastq.gz";
        push( @result_files, "${result_dir}/${unmapped_fastq_file}" );
    } else { #paired seprately
        my $unmapped_fastq_file1 =$sample_name . ".1.fastq.gz";
        my $unmapped_fastq_file2 =$sample_name . ".2.fastq.gz";
        push( @result_files, "${result_dir}/${unmapped_fastq_file1}" );
        push( @result_files, "${result_dir}/${unmapped_fastq_file2}" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
