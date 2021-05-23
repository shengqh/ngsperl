#!/usr/bin/perl
package Alignment::BWA;

use strict;
use warnings;
use File::Basename;
use POSIX;
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
  $self->{_suffix} = "_bwa";
  $self->{_use_tmp_folder} = 1;
  $self->init_docker_prefix(__PACKAGE__);
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my ($sort_memory, $isMB) = getMemoryPerThread($memory, $thread);
  if ($isMB) {
    $sort_memory = $sort_memory . "M";
  }else{
    $sort_memory = $sort_memory . "G";
  }

  my $selfname = $self->{_name};

  my $use_sambamba = get_option( $config, $section, "use_sambamba",                0 );
  my $cleansam                = get_option( $config, $section, "cleansam",                0 );
  my $chromosome_grep_pattern = get_option( $config, $section, "chromosome_grep_pattern", "" );
  my $sortByCoordinate        = get_option( $config, $section, "sort_by_coordinate",      1 );
  my $bwa_only        = get_option( $config, $section, "bwa_only",      0 );
  my $alignmentOnly        = get_option( $config, $section, "alignmentOnly",      0 );
  my $mark_duplicates         = hasMarkDuplicate( $config->{$section} );
  my $rg_name_regex        = get_option( $config, $section, "rg_name_regex",      "" );
  my $rg_id_regex        = get_option( $config, $section, "rg_id_regex",      "" );

  my $sort_by_thread        = get_option( $config, $section, "sort_by_thread",     0 );
  my $samtools_sort_thread = $sort_by_thread ? "-@ $thread":"";
  my $sambamba_sort_thread = $sort_by_thread ? "-t $thread":"";

  my $index_by_thread        = get_option( $config, $section, "index_by_thread",     0 );
  my $samtools_index_thread = $index_by_thread ? "-@ $thread":"";
  my $sambamba_index_thread = $index_by_thread ? "-t $thread":"";

  $option = $option . " -M";

  if ( !( $option =~ /\s-t\s/ ) ) {
    if ( $thread > 1 ) {
      $option = $option . " -t " . $thread;
    }
  }

  my $bwa_index = $config->{$section}{bwa_index};
  if ( !defined $bwa_index ) {
    $bwa_index = $config->{$section}{fasta_file} or die "define ${section}::bwa_index first";
  }

  my $py_script = dirname(__FILE__) . "/bamStat.py";
  if ( !-e $py_script ) {
    die "File not found : " . $py_script;
  }

  my $picard_jar;
  
  if($cleansam || $mark_duplicates) {
    $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1, not $self->using_docker() );
  } 

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $unsorted_bam_file = $sample_name . ".unsorted.bam";
    my $bam_stat = $sample_name . ".bamstat";
    my $clean_bam_file    = $sample_name . ".unsorted.clean.bam";

    my $sorted_bam_file    = $sample_name . ($sortByCoordinate? ".sortedByCoord.bam":".sortedByQuery.bam");
    my $sorted_bam_file_index    = "${sorted_bam_file}.bai";

    my $rmdup_bam_file    = $sample_name . ($sortByCoordinate? ".sortedByCoord.rmdup.bam":".sortedByQuery.rmdup.bam");
    my $rmdup_bam_file_index    = change_extension($rmdup_bam_file, ".bai");

    my $tag               = get_bam_tag($sample_name);

    my $rg_sample_name = $sample_name;
    if ($rg_name_regex ne ""){
      if ($sample_name =~ /$rg_name_regex/ ) {
        $rg_sample_name = $1;
      }
    }
    my $rg_sample_id = "1";
    if ($rg_id_regex ne ""){
      if ($sample_name =~ /$rg_id_regex/ ) {
        $rg_sample_id = $1;
      }
    }

    my $rg;
    if ($alignmentOnly) { #only alignment, no sort or other works. For UMI pipeline
      $rg ="";
    } else {
      $rg = "-R \"\@RG\\tID:${rg_sample_id}\\tPU:${rg_sample_name}\\tLB:${rg_sample_name}\\tSM:${rg_sample_name}\\tPL:ILLUMINA\"";
    }

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $log_desc = $cluster->get_log_description($log);

    my $final_file =  $mark_duplicates ? $rmdup_bam_file : $sorted_bam_file;
    my $check_file = $sortByCoordinate? ($mark_duplicates?$rmdup_bam_file_index:$sorted_bam_file_index) :$final_file;

    my $sample_file_0 = $sample_files[0];
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $check_file, $init_command, 0, $sample_file_0 );

    print $sh "
if [[ ! -s $result_dir/$final_file ]]; then
  \$MYCMD ./$pbs_name 
fi
";

    my $localized_files = [];
    @sample_files = @{$self->localize_files_in_tmp_folder($pbs, \@sample_files, $localized_files)};
    $sample_file_0 = $sample_files[0];
    my $sample_files_str = ( scalar(@sample_files) == 2 ) ? "\"" . $sample_file_0 . "\" \"" . $sample_files[1] . "\"" : "\"" . $sample_file_0 . "\"";

    print $pbs "
if [[ ! -s $unsorted_bam_file ]]; then
  echo bwa_mem=`date`
  bwa mem $option $rg $bwa_index $sample_files_str | samtools view -bS -o $unsorted_bam_file
  bwa 2>\&1 | grep Version | cut -d ' ' -f2 | cut -d '-' -f1 | awk '{print \"bwa,v\"\$1}' > ${sample_name}.bwa.version
fi
";

    if ($bwa_only) {
      $self->close_pbs( $pbs, $pbs_file );
      next;
    }

    print $pbs "
if [[ -s $unsorted_bam_file ]]; then
  echo bamStat=`date` 
  python3 $py_script -i $unsorted_bam_file -o $bam_stat
fi
";

    my $rmlist = "";

    if ($alignmentOnly) { #only alignment, no sort or other works. For UMI pipeline
      print $pbs "
if [ -s $unsorted_bam_file ]; then
  echo flagstat=`date` 
  samtools flagstat $unsorted_bam_file > ${unsorted_bam_file}.stat 
fi
";
      if ($rmlist ne "") {
        print $pbs "
if [ -s $unsorted_bam_file ]; then
  rm $rmlist
fi
";
      }

      $self->clean_temp_files($pbs, $localized_files);

      $self->close_pbs( $pbs, $pbs_file );
      next;
    }

    if ($cleansam) {
      print $pbs "
if [[ (-s $unsorted_bam_file) && ((1 -eq \$1) || (! -s $clean_bam_file)) ]]; then
  echo CleanSam=`date`
  java -jar $picard_jar CleanSam VALIDATION_STRINGENCY=SILENT I=$unsorted_bam_file O=$clean_bam_file
fi
";
      $rmlist            = $rmlist . " " . $unsorted_bam_file;
      $unsorted_bam_file = $clean_bam_file;
    }

    if ($sortByCoordinate) {
      my $chromosome_grep_command = getChromosomeFilterCommand( $sorted_bam_file, $chromosome_grep_pattern );
      if($use_sambamba){
        print $pbs "    
if [[ (-s $unsorted_bam_file) && ((1 -eq \$1) || (! -s $sorted_bam_file)) ]]; then
  echo sort_bam=`date`
  sambamba sort -m $sort_memory $sambamba_sort_thread --tmpdir tmp_${sample_name} -o $sorted_bam_file $unsorted_bam_file
fi

if [[ (-s $sorted_bam_file) && ((1 -eq \$1) || (! -s ${sorted_bam_file}.bai)) ]]; then
  echo index_bam=`date`
  sambamba index $sambamba_index_thread $sorted_bam_file 
  samtools idxstats $sorted_bam_file > ${sorted_bam_file}.chromosome.count
fi

$chromosome_grep_command
";
      }else{
        print $pbs "    
if [[ (-s $unsorted_bam_file) && ((1 -eq \$1) || (! -s $sorted_bam_file)) ]]; then
  echo sort_bam=`date`
  samtools sort -m $sort_memory $samtools_sort_thread -T /tmp/bwa_${sample_name} -o $sorted_bam_file $unsorted_bam_file
fi

if [[ (-s $sorted_bam_file) && ((1 -eq \$1) || (! -s ${sorted_bam_file}.bai)) ]]; then
  echo index_bam=`date`
  samtools index $samtools_index_thread $sorted_bam_file 
  samtools idxstats $sorted_bam_file > ${sorted_bam_file}.chromosome.count
fi

$chromosome_grep_command

";
      }
    }
    else {
      if($use_sambamba){
        print $pbs "    
if [[ (-s $unsorted_bam_file) && ((1 -eq \$1) || (! -s $sorted_bam_file)) ]]; then
  echo sort_bam=`date`
  sambamba sort -m $sort_memory $sambamba_sort_thread --sort-picard -o $sorted_bam_file $unsorted_bam_file
fi
";
      }else{
        print $pbs "    
if [[ (-s $unsorted_bam_file) && ((1 -eq \$1) || (! -s $sorted_bam_file)) ]]; then
  echo sort_bam=`date`
  samtools sort -m $sort_memory $samtools_sort_thread -T /tmp/bwa_${sample_name} -n -o $sorted_bam_file $unsorted_bam_file
fi
";
      }
    }
    $rmlist = $rmlist . " " . $unsorted_bam_file;

    if ($mark_duplicates) {
      my $ASSUME_SORT_ORDER = $sortByCoordinate ? "coordinate":"queryname";
      print $pbs "
if [[ (-s $sorted_bam_file) && ((1 -eq \$1) || (! -s $rmdup_bam_file)) ]]; then
  echo MarkDuplicate=`date` 
  java -jar $picard_jar MarkDuplicates \\
    INPUT=$sorted_bam_file \\
    OUTPUT=$rmdup_bam_file \\
    METRICS_FILE=${rmdup_bam_file}.metrics \\
    VALIDATION_STRINGENCY=SILENT \\
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \\
    ASSUME_SORT_ORDER=\"$ASSUME_SORT_ORDER\" \\
    CLEAR_DT=false \\
    ADD_PG_TAG_TO_READS=false \\
    CREATE_INDEX=true \\
    REMOVE_DUPLICATES=true
fi

samtools idxstats $rmdup_bam_file > ${rmdup_bam_file}.chromosome.count

";
      $rmlist = $rmlist . " " . $sorted_bam_file;
    }

    print $pbs "
if [[ (-s $final_file) && (-s $check_file) ]]; then
  echo flagstat=`date` 
  samtools flagstat $final_file > ${final_file}.stat 
fi
";

  if ($rmlist ne "") {
    print $pbs "
if [[ -s $check_file && -s $bam_stat ]]; then
  rm $rmlist tmp_${sample_name}
fi
";
  }

    $self->clean_temp_files($pbs, $localized_files);

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all BWA tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $sortByCoordinate = get_option( $config, $section, "sort_by_coordinate", 1 );
  my $mark_duplicates  = hasMarkDuplicate( $config->{$section} );
  my %raw_files        = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $final_file;
    my $alignmentOnly        = get_option( $config, $section, "alignmentOnly",      0 );
    if ($alignmentOnly) { #only alignment, no sort or other works. For UMI pipeline
      $final_file=$sample_name.".unsorted.bam";
    } else { #all other bwa tasks, with sort or mark_duplicates
      my $sorted_bam_file    = $sample_name . ($sortByCoordinate? ".sortedByCoord.bam":".sortedByQuery.bam");
      my $rmdup_bam_file    = $sample_name . ($sortByCoordinate? ".sortedByCoord.rmdup.bam":".sortedByQuery.rmdup.bam");
      $final_file =  $mark_duplicates ? $rmdup_bam_file : $sorted_bam_file;
    }
    $final_file = "${result_dir}/$final_file";
    my @result_files = ();
    push( @result_files, $final_file );
    push( @result_files, $final_file . ".stat" );
    push( @result_files, $final_file . ".chromosome.count" );
    push( @result_files, "${result_dir}/${sample_name}.bamstat" );
    push( @result_files, "${result_dir}/${sample_name}.bwa.version" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
