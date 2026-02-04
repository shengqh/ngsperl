#!/usr/bin/perl
package GATK4::CollectReadCounts;

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
  $self->{_suffix} = "_crc";
  bless $self, $class;
  return $self;
}

#https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/cnv_common_tasks.wdl
sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $java_option = $self->get_java_option($config, $section, $memory);

  #parameter files
  $self->get_docker_value(1);

  my $preprocessed_intervals = parse_param_file( $config, $section, "preprocessed_intervals", 1 );

  my $ref_fasta_dict = get_param_file( $config->{$section}{ref_fasta_dict}, "ref_fasta_dict", 1 );
  my $ref_fasta      = get_param_file( $config->{$section}{ref_fasta},      "ref_fasta",      1 );

  my $count_format = get_option( $config, $section, "count_format", "HDF5" );

  #make PBS
  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $final_file = ( $count_format eq "HDF5" ) ? $sample_name . ".count.hdf5" : $sample_name . ".count.tsv";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command );

    my $other_exts = $sample_files[0] =~ /.cram$/ ? [".crai"] : [".bai"];

    my $localized_files = [];
    @sample_files = @{$self->localize_files_in_tmp_folder($pbs, \@sample_files, $localized_files, $other_exts)};
    my $sampleFile    = $sample_files[0];

    my $cur_option;
    if(is_string($option)) {
      $cur_option = $option;
    }else{
      $cur_option = $option->{$sample_name};
    }
    print $pbs "  

gatk --java-options \"$java_option\" CollectReadCounts $cur_option \\
  -L ${preprocessed_intervals} \\
  --input $sampleFile \\
  --reference $ref_fasta \\
  --format $count_format \\
  --interval-merging-rule OVERLAPPING_ONLY \\
  --output $final_file

status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.failed
  rm -f $final_file
else
  touch $sample_name.succeed
fi
            
rm -rf .cache .conda .config .theano
";

    $self->clean_temp_files($pbs, $localized_files);

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all GATK CNV tasks. \n ";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $count_format = get_option( $config, $section, "count_format", "HDF5" );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $final_file = ( $count_format eq "HDF5" ) ? $sample_name . ".count.hdf5" : $sample_name . ".count.tsv";
    my @result_files = ();
    push( @result_files, "${result_dir}/${final_file}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }

  return $result;
}

1;
