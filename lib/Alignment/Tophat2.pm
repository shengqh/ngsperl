#!/usr/bin/perl
package Alignment::Tophat2;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_th2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  $option = $option . " --keep-fasta-order --no-coverage-search";

  my $sort_by_query = get_option( $config, $section, "sort_by_query", 0 );
  my $rename_bam    = get_option( $config, $section, "rename_bam",    0 );

  my $bowtie2_index = get_option( $config, $section, "bowtie2_index" );
  my %fqFiles = %{ get_raw_files( $config, $section ) };

  my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "${section}::transcript_gtf", 0 );
  my $transcript_gtf_index = $config->{$section}{transcript_gtf_index};
  if ( !defined $transcript_gtf ) {
    $transcript_gtf = get_param_file( $config->{general}{transcript_gtf}, "general::transcript_gtf", 0 );
  }
  if ( !defined $transcript_gtf_index ) {
    $transcript_gtf_index = $config->{general}{transcript_gtf_index};
  }

  my $has_gtf_file   = file_exists($transcript_gtf);
  my $has_index_file = transcript_gtf_index_exists($transcript_gtf_index);

  if ( $has_gtf_file && !defined $transcript_gtf_index ) {
    die "transcript_gtf was defined but transcript_gtf_index was not defined, you should defined transcript_gtf_index to cache the parsing result.";
  }

  if ( defined $transcript_gtf_index && !$has_index_file ) {
    if ($has_gtf_file) {
      print "transcript_gtf_index $transcript_gtf_index defined but not exists, you may run the script once to cache the index.\n";
    }
    else {
      die "transcript_gtf_index $transcript_gtf_index defined but not exists, and transcript_gtf is not defined.\n";
    }
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %fqFiles ) {
    my @sample_files = @{ $fqFiles{$sample_name} };
    my $samples = join( " ", @sample_files );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir     = create_directory_or_die( $result_dir . "/$sample_name" );
    my $rgline      = "--rg-id $sample_name --rg-sample $sample_name --rg-library $sample_name";
    my $tophat2file = "accepted_hits.bam";

    my $gtfstr = "";
    if ($has_gtf_file) {
      $gtfstr = "-G $transcript_gtf --transcriptome-index=$transcript_gtf_index";
    }
    elsif ($has_index_file) {
      $gtfstr = "--transcriptome-index=$transcript_gtf_index";
    }

    my $final_file     = $rename_bam    ? "${sample_name}.bam"                                                                       : "accepted_hits.bam";
    my $sort_cmd       = $sort_by_query ? "samtools sort -n accepted_hits.bam ${sample_name}.sortedname"                  : "";
    my $rename_bam_cmd = $rename_bam    ? "mv accepted_hits.bam ${sample_name}.bam\nmv accepted_hits.bam.bai ${sample_name}.bam.bai" : "";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    print $pbs "
tophat2 $option $rgline $gtfstr -o . $bowtie2_index $samples

if [ -s $tophat2file ]; then
  samtools index $tophat2file
  $sort_cmd
  $rename_bam_cmd
  samtools flagstat $final_file > ${final_file}.stat
fi 
   
";

    $self->close_pbs( $pbs, $pbs_file );

    print $sh "\$MYCMD ./$pbs_name \n";
  }
  print $sh "exit 0\n";
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $sort_by_query = get_option_value( $config->{$section}{sort_by_query}, 0 );
  my $rename_bam    = get_option_value( $config->{$section}{rename_bam},    0 );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    if ($rename_bam) {
      push( @result_files, "${result_dir}/${sample_name}/${sample_name}.bam" );
    }
    else {
      push( @result_files, "${result_dir}/${sample_name}/accepted_hits.bam" );
    }

    if ($sort_by_query) {
      push( @result_files, "${result_dir}/${sample_name}/${sample_name}.sortedname.bam" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
