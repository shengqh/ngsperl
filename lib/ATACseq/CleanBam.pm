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
use Alignment::AlignmentUtils;

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $sort_memory = $thread == 1 ? $memory : "4G";

  my $remove_chromosome = get_option( $config, $section, "remove_chromosome", "M" );
  my $keep_chromosome   = get_option( $config, $section, "keep_chromosome",   "" );
  my $minimum_maq       = get_option( $config, $section, "minimum_maq",       30 );
  my $isSortedByCoordinate = get_option( $config, $section, "is_sorted_by_coordinate" );
  my $blacklistfile = get_param_file( $config->{$section}{"blacklist_file"}, "blacklist_file", 0 );
  my $maxInsertSize = get_option( $config, $section, "maximum_insert_size", 0 );
  my $is_paired_end = get_is_paired_end_option( $config, $section, 1 );

  my $mark_duplicates = hasMarkDuplicate( $config->{$section}, 1 );
  my $picard_jar      = "";
  if ($mark_duplicates) {
    $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1, not $self->using_docker() );
  }

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

    my $redupFile = undef;
    my $finalFile = undef;
    if ($mark_duplicates) {
      $redupFile = $sample_name . ".rmdup.bam";
      $finalFile = $sample_name . ".rmdup.noChr" . $remove_chromosome . ".bam";
    }
    else {
      $finalFile = $sample_name . ".noChr" . $remove_chromosome . ".bam";
    }
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $finalFile . ".bai" );
    my $rmlist = "";
    if ( !$isSortedByCoordinate ) {
      my $sorted = $sample_name . ".sortedByCoordinate.bam";
      print $pbs "
if [ ! -s $sorted ]; then
  echo SortByCoordinate=`date` 
  samtools sort -m $sort_memory -o $sorted $sampleFile 
fi
";
      $sampleFile = $sorted;
      $rmlist     = $sorted;
    }

    my $filterFile;
    if ($is_paired_end) {
      $filterFile = "${sample_name}.filter.bam";
    }
    else {
      $filterFile = $finalFile;
    }

    my $filterOption;
    if ( defined $blacklistfile ) {
      $filterOption = "| samtools view -b -U $filterFile -o ${sample_name}.discard.bam -L $blacklistfile";
      $rmlist       = $rmlist . " ${sample_name}.discard.bam";
    }
    else {
      $filterOption = "> $filterFile";
    }

    my $input = $sampleFile;
    if ($mark_duplicates) {
      print $pbs "
if [[ -s $sampleFile && ! -s ${redupFile}.bai ]]; then
  echo RemoveDuplicate=`date` 
  java -jar $picard_jar MarkDuplicates I=$sampleFile O=$redupFile ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT M=${redupFile}.metrics
  samtools index $redupFile
fi
";

      $rmlist = $rmlist . " $redupFile ${redupFile}.bai";

      $input = $redupFile;
    }

    if ( defined $maxInsertSize && $maxInsertSize > 0 ) {
      my $insertFile = $sample_name . ".insertsize.bam";
      print $pbs "
if [[ -s $input && ! -s ${insertFile}.bai ]]; then
  echo FilterInsertSize=`date`
  samtools view -h $input | awk -F '\\t' 'function abs(v) {return v < 0 ? -v : v} {if(substr(\$1,1,1) == \"\@\") {print \$0} else { if(abs(\$9) < $maxInsertSize) print \$0}}' | samtools view -b > $insertFile
  samtools index $insertFile
fi
";
      $input  = $insertFile;
      $rmlist = $rmlist . " $insertFile  ${insertFile}.bai";
    }

    print $pbs "
if [[ -s $input && ! -s ${filterFile}.bai ]]; then 
  echo FilterBam=`date` 
  samtools idxstats $input | cut -f 1 | grep -v $remove_chromosome $keep_chromosome | xargs samtools view $option -b -q $minimum_maq $input $filterOption 
  samtools index $filterFile 
fi
";
    if ($is_paired_end) {
      print $pbs "
if [[ -s ${filterFile}.bai && ! -s ${finalFile}.bai ]]; then 
  echo RemoveUnpaired=`date` 
  samtools sort -n $filterFile | samtools fixmate -O bam - -| samtools view $option -b | samtools sort -T $sample_name -o $finalFile 
  samtools index $finalFile 
  samtools idxstats $finalFile > ${finalFile}.chromosome.count
fi
";
      $rmlist = $rmlist . " $filterFile ${filterFile}.bai";
    }

    print $pbs "
if [[ -s ${finalFile}.bai && ! -s ${finalFile}.stat ]]; then
  echo Stat=`date`
  samtools flagstat $finalFile > ${finalFile}.stat
fi
";

    print $pbs "
if [ -s ${finalFile}.stat ]; then 
  rm $rmlist  
fi
";

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks . \n ";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $remove_chromosome = get_option( $config, $section, "remove_chromosome" );
  my $mark_duplicates = hasMarkDuplicate( $config->{$section} );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $finalFile = $mark_duplicates ? $sample_name . ".rmdup.noChr" . $remove_chromosome . ".bam" : $sample_name . ".noChr" . $remove_chromosome . ".bam";

    my @result_files = ();
    push( @result_files, "${result_dir}/${finalFile}" );
    push( @result_files, "${result_dir}/${finalFile}.stat" );
    push( @result_files, "${result_dir}/${finalFile}.chromosome.count" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
