#!/usr/bin/perl
package ATACseq::CleanBam;

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
  $self->{_suffix} = "_cb";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1 );
  my $remove_chromosome = get_option( $config, $section, "remove_chromosome" );
  my $keep_chromosome   = get_option( $config, $section, "keep_chromosome", "" );
  my $minimum_maq       = get_option( $config, $section, "minimum_maq", 10 );
  my $isSortedByName    = get_option( $config, $section, "is_sorted_by_name" );

  if ( $keep_chromosome !~ /^\s*$/ ) {
    $keep_chromosome = "| grep $keep_chromosome";
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile   = $sample_files[0];

    my $fixFile     = $sample_name . ".fix.bam";
    my $fixSortFile = $sample_name . ".fix.sorted.bam";
    my $redupFile   = $sample_name . ".fix.rmdup.bam";
    my $finalFile  = $sample_name . ".fix.rmdup.noChr" . $remove_chromosome . ".bam";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $finalFile );
    my $rmlist = "";
    if ( !$isSortedByName ) {
      my $sorted = $sample_name . ".sortByName.bam";
      print $pbs "
if [ ! -s $sorted ]; then
  echo SortByName=`date` 
  samtools sort -n -o $sorted $sampleFile 
fi
";
      $sampleFile = $sorted;
      $rmlist     = $sorted;
    }

    print $pbs "
if [[ -s $sampleFile && ! -s $fixFile ]]; then
  echo FixMate=`date` 
  samtools fixmate $sampleFile  > $fixFile 
fi

if [[ -s fixFile && ! -s $fixSortFile ]]; then
  echo SortByCoordinate=`date` 
  samtools sort -@ $thread -m $memory -o $fixSortFile $fixFile  
fi

if [[ -s $fixSortFile && ! -s $redupFile ]]; then
  echo RemoveDuplicate=`date` 
  java $option -jar $picard_jar MarkDuplicates I=$fixSortFile O=$redupFile ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT M=${redupFile}.metrics
fi

if [[ -s $redupFile && ! -s $finalFile ]]; then
  echo RemoveChromosome=`date` 
  samtools idxstats $fixSortFile | cut -f 1 | grep -v $remove_chromosome $keep_chromosome | xargs samtools view -b -q $minimum_maq $redupFile > $finalFile 
fi

if [[ -s $finalFile && ! -s ${finalFile}.bai ]]; then
  echo BamIndex=`date` 
  samtools index $finalFile
  samtools flagstat $finalFile > ${finalFile}.stat
  rm $rmlist $fixFile $fixSortFile $redupFile
fi
";

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $remove_chromosome = get_option( $config, $section, "remove_chromosome" );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $finalFile  = $sample_name . ".fix.rmdup.noChr" . $remove_chromosome . ".bam";

    my @result_files = ();
    push( @result_files, "${result_dir}/${finalFile}" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
