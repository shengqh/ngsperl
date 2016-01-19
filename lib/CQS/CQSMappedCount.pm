#!/usr/bin/perl
package CQS::CQSMappedCount;

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
  $self->{_suffix} = "_ct";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqstools = get_cqstools( $config, $section, 1 );
  my $samtools  = get_param_file( $config->{$section}{samtools},   "samtools",   1 );
  my $gffFile   = get_param_file( $config->{$section}{gff_file},   "gff_file",   1 );
  my $fastaFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 0 );

  my $min_overlap = $config->{$section}{min_overlap};
  if ( defined $min_overlap ) {
    $option = $option . " --min_overlap $min_overlap";
  }

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
    my $countFile = $sample_name . ".count";

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
    my $log     = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);
    
    my $name = $self->{_name};

    open( my $out, ">$pbs_file" ) or die $!;
    print $out "$pbs_desc
$log_desc

$path_file

cd $cur_dir

if [ -s $countFile ]; then
  echo job has already been done. if you want to do again, delete $countFile and submit job again.
  exit 0
fi
  
echo ${name}_start=`date` 

mono-sgen $cqstools mapped_count $option --samtools $samtools -i $bam_file -g $gffFile -o $countFile $seqcountFile $fastqFile

echo ${name}_end=`date`

exit 0 
";

    close $out;

    print "$pbs_file created \n";
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all trna_count tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $fasta_format = $config->{$section}{fasta_format};
  if ( !defined $fasta_format ) {
    $fasta_format = 0;
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $cur_dir = $result_dir . "/$sample_name";

    my @bam_files = @{ $raw_files{$sample_name} };

    my @result_files = ();
    my $countFile   = "${cur_dir}/${sample_name}.count";
    push( @result_files, $countFile );
    push( @result_files, "${countFile}.mapped.xml" );
    push( @result_files, "${cur_dir}/${sample_name}.info" );

    my $unmapped;
    if ($fasta_format) {
      $unmapped = change_extension( $countFile, ".unmapped.fasta.gz" );
    }
    else {
      $unmapped = change_extension( $countFile, ".unmapped.fastq.gz" );
    }
    push( @result_files, $unmapped );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
