#!/usr/bin/perl
package Count::FeatureCounts;

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
  $self->{_suffix} = "_fc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $is_paired_end = get_is_paired_end_option( $config, $section);
  if ($is_paired_end) {
    $option = $option . " -p --countReadPairs";
  }

  if ( $option !~ /-F SAF/ ) {
    if ( $option !~ /-g/ ) {
      $option = $option . " -g gene_id";
    }

    if ( $option !~ /-t/ ) {
      $option = $option . " -t exon";
    }
  }

  my $gffFile = parse_param_file( $config, $section, "gff_file", 1 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @bam_files  = @{ $raw_files{$sample_name} };
    my $bam_file   = join( " ", @bam_files );
    my $final_file = "${sample_name}.count";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";
    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    print $pbs "featureCounts $option -T $thread -a $gffFile -o $final_file $bam_file 

featureCounts -v 2>\&1 | grep featureCounts | cut -d ' ' -f2 | awk '{print \"featureCounts,\"\$1}' > ${sample_name}.count.featureCounts.version
";

    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all mirna_count tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    push( @result_files, "${result_dir}/${sample_name}.count" );
    push( @result_files, "${result_dir}/${sample_name}.count.summary" );
    push( @result_files, "${result_dir}/${sample_name}.count.featureCounts.version" );
    
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
