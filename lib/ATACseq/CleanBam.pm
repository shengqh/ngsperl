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
  $self->{_use_tmp_folder} = 1;
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

  my $remove_duplicates = get_option( $config, $section, "remove_duplicates", 1 );
  my $mark_duplicates = get_option( $config, $section, "mark_duplicates", 1 );
  my $picard_jar      = "";
  if ($remove_duplicates || $mark_duplicates) {
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

    my $redupFile = undef;
    my $finalFile = undef;
    my $remove_duplicates_option = "";
    if($remove_duplicates){
      $remove_duplicates_option = "REMOVE_DUPLICATES=true";
      $redupFile = $sample_name . ".rmdup.bam";
      $finalFile = $sample_name . ".rmdup.noChr" . $remove_chromosome . ".bam";
    } elsif($mark_duplicates) {
      $remove_duplicates_option = "REMOVE_DUPLICATES=false";
      $redupFile = $sample_name . ".markeddup.bam";
      $finalFile = $sample_name . ".markeddup.noChr" . $remove_chromosome . ".bam";
    } else {
      $finalFile = $sample_name . ".noChr" . $remove_chromosome . ".bam";
    }
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $finalFile . ".bai" );

    print $sh "
if [[ ! -s $result_dir/${finalFile}.bai ]]; then
  \$MYCMD ./$pbs_name 
fi
";

    my $localized_files = [];
    @sample_files = @{$self->localize_files_in_tmp_folder($pbs, \@sample_files, $localized_files, [".bai"])};
    my $sampleFile   = $sample_files[0];

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
      $filterOption = "| samtools view -b -U ${filterFile}.tmp.bam -o ${sample_name}.discard.bam -L $blacklistfile";
      $rmlist       = $rmlist . " ${sample_name}.discard.bam";
    }
    else {
      $filterOption = "> ${filterFile}.tmp.bam";
    }

    my $input = $sampleFile;
    if ($remove_duplicates || $mark_duplicates) {
      print $pbs "
if [[ -s $sampleFile && ! -s ${redupFile}.bai ]]; then
  rm -f ${redupFile}.failed ${redupFile}.succeed

  echo MarkDuplicates=`date` 
  java -jar $picard_jar MarkDuplicates I=${sampleFile} O=${redupFile}.tmp.bam $remove_duplicates_option ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT M=${redupFile}.metrics

  status=\$?
  if [[ \$status -ne 0 ]]; then
    echo MarkDuplicates failed with exit code \$status
    echo \$status > ${redupFile}.failed
    rm -f ${redupFile}.tmp.bam*
    exit \$status
  else
    mv ${redupFile}.tmp.bam ${redupFile}.failed
    touch ${redupFile}.succeed
  fi

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
  rm -f ${insertFile}.failed ${insertFile}.succeed

  echo FilterInsertSize=`date`
  samtools view -h $input | awk -F '\\t' 'function abs(v) {return v < 0 ? -v : v} {if(substr(\$1,1,1) == \"\@\") {print \$0} else { if(abs(\$9) < $maxInsertSize) print \$0}}' | samtools view -b > $insertFile.tmp.bam

  status=\$?
  if [[ \$status -ne 0 ]]; then
    echo FilterInsertSize failed with exit code \$status
    echo \$status > ${insertFile}.failed
    rm -f $insertFile.tmp.bam
    exit \$status
  else
    mv $insertFile.tmp.bam $insertFile
    touch ${insertFile}.succeed
  fi

  samtools index $insertFile
fi
";
      $input  = $insertFile;
      $rmlist = $rmlist . " $insertFile  ${insertFile}.bai";
    }

    print $pbs "
if [[ -s $input && ! -s ${filterFile}.bai ]]; then 
  rm -f ${filterFile}.failed ${filterFile}.succeed

  echo FilterBam=`date` 
  samtools idxstats $input | cut -f 1 | grep -v $remove_chromosome $keep_chromosome | xargs samtools view $option -b -q $minimum_maq $input $filterOption 

  status=\$?
  if [[ \$status -ne 0 ]]; then
    echo FilterBam failed with exit code \$status
    echo \$status > ${filterFile}.failed
    rm -f ${filterFile}.tmp.bam*
    exit \$status
  else
    mv ${filterFile}.tmp.bam $filterFile
    touch ${filterFile}.succeed
  fi

  samtools index $filterFile 
fi
";
    if ($is_paired_end) {
      print $pbs "
if [[ -s ${filterFile}.bai && ! -s ${finalFile}.bai ]]; then 
  rm -f ${finalFile}.failed ${finalFile}.succeed

  echo RemoveUnpaired=`date` 
  samtools sort -n $filterFile | samtools fixmate -O bam - -| samtools view $option -b | samtools sort -T $sample_name -o ${finalFile}.tmp.bam 
  
  status=\$?
  if [[ \$status -ne 0 ]]; then
    echo RemoveUnpaired failed with exit code \$status
    echo \$status > ${finalFile}.failed
    rm -f ${finalFile}.tmp.bam*
    exit \$status
  else
    mv ${finalFile}.tmp.bam ${finalFile}
    touch ${finalFile}.succeed
  fi

  samtools index $finalFile 
fi
";
      $rmlist = $rmlist . " $filterFile ${filterFile}.bai";
    }

    print $pbs "
if [[ -s ${finalFile}.bai && ! -s ${finalFile}.stat ]]; then
  echo Stat=`date`
  samtools idxstats $finalFile > ${finalFile}.chromosome.count
  samtools flagstat $finalFile > ${finalFile}.stat
fi
";

    print $pbs "
if [ -s ${finalFile}.stat ]; then 
  rm -f $rmlist  
fi
";

    $self->clean_temp_files($pbs, $localized_files);

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
  my $remove_duplicates = get_option( $config, $section, "remove_duplicates", 1 );
  my $mark_duplicates = get_option( $config, $section, "mark_duplicates", 1 );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $finalFile = undef;
    if($remove_duplicates){
      $finalFile = $sample_name . ".rmdup.noChr" . $remove_chromosome . ".bam";
    } elsif($mark_duplicates) {
      $finalFile = $sample_name . ".markeddup.noChr" . $remove_chromosome . ".bam";
    } else {
      $finalFile = $sample_name . ".noChr" . $remove_chromosome . ".bam";
    }

    my @result_files = ();
    push( @result_files, "${result_dir}/${finalFile}" );
    push( @result_files, "${result_dir}/${finalFile}.stat" );
    push( @result_files, "${result_dir}/${finalFile}.chromosome.count" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
