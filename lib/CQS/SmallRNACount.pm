#!/usr/bin/perl
package CQS::SmallRNACount;

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
  $self->{_suffix} = "_sc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqstools = get_cqstools( $config, $section, 1 );
  my $coordinate_file = get_param_file( $config->{$section}{coordinate_file}, "coordinate_file", 1 );
  my $fastaFile       = get_param_file( $config->{$section}{fasta_file},      "fasta_file",      0 );

  if ( defined $fastaFile ) {
    $option = $option . " -f $fastaFile";
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my %seqcount_files = ();
  if ( has_raw_files( $config, $section, "seqcount" ) ) {
    %seqcount_files = %{ get_raw_files( $config, $section, "seqcount" ) };
  }

  my %fastqFiles = ();
  if ( has_raw_files( $config, $section, "fastq_files" ) ) {
    %fastqFiles = %{ get_raw_files( $config, $section, "fastq_files" ) };
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @bam_files  = @{ $raw_files{$sample_name} };
    my $bam_file   = $bam_files[0];
    my $final_file = $sample_name . ".count";
    my $final_xml_file = $sample_name . ".count.mapped.xml";

    my $seqcountFile = "";
    if ( defined $seqcount_files{$sample_name} ) {
      my @seqcounts = @{ $seqcount_files{$sample_name} };
      my $seqcount  = $seqcounts[0];
      $seqcountFile = " -c $seqcount";
    }

    my $fastqFile = "";
    if ( defined $fastqFiles{$sample_name} ) {
      my @files = @{ $fastqFiles{$sample_name} };
      my $file  = $files[0];
      $fastqFile = " -q $file";
    }

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_xml_file );
    print $pbs "mono $cqstools smallrna_count $option -i $bam_file -g $coordinate_file $seqcountFile $fastqFile -o $final_file";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $unmapped_fastq = $option =~ /--unmapped_fastq/;

  my %seqcount_files = ();
  if ( defined $config->{$section}{"seqcount"} || defined $config->{$section}{"seqcount_ref"} ) {
    %seqcount_files = %{ get_raw_files( $config, $section, "seqcount" ) };
  }

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $cur_dir = $result_dir . "/$sample_name";

    my @result_files     = ();
    my $countFile        = "${cur_dir}/${sample_name}.count";
    my $tRNAPositionFile = "${cur_dir}/${sample_name}.tRNA.position";
    push( @result_files, $countFile );
    push( @result_files, $tRNAPositionFile );
    push( @result_files, "${countFile}.mapped.xml" );
    push( @result_files, "${cur_dir}/${sample_name}.info" );

    if ($unmapped_fastq) {
      my $unmapped = change_extension( $countFile, ".unmapped.fastq.gz" );
      push( @result_files, $unmapped );
      if ( defined $seqcount_files{$sample_name} ) {
        push( @result_files, $unmapped . ".dupcount" );
      }
    }

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
