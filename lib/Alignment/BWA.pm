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
  my $chromosome_grep_pattern = get_option( $config, $section, "chromosome_grep_pattern", "" );
  my $sortByCoordinate        = get_option( $config, $section, "sort_by_coordinate",      1 );
  my $bwa_only        = get_option( $config, $section, "bwa_only",      0 );
  my $alignmentOnly        = get_option( $config, $section, "alignmentOnly",      0 );
  my $rg_name_regex        = get_option( $config, $section, "rg_name_regex",      "" );
  my $rg_id_regex        = get_option( $config, $section, "rg_id_regex",      "" );

  my $output_unmapped_fastq = get_option( $config, $section, "output_unmapped_fastq", 0 );

  my $sort_by_thread        = get_option( $config, $section, "sort_by_thread",     0 );
  my $samtools_sort_thread = $sort_by_thread ? "-@ $thread":"";
  my $sambamba_sort_thread = $sort_by_thread ? "-t $thread":"";

  my $index_by_thread        = get_option( $config, $section, "index_by_thread",     0 );
  my $samtools_index_thread = $index_by_thread ? "-@ $thread":"";
  my $sambamba_index_thread = $index_by_thread ? "-t $thread":"";

  my $do_bam_stat = get_option( $config, $section, "do_bam_stat",     1 );

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

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $unsorted_bam_file = $sample_name . ".unsorted.bam";
    my $bam_stat = $sample_name . ".bamstat";

    my $sorted_bam_file    = $sample_name . ($sortByCoordinate? ".sortedByCoord.bam":".sortedByQuery.bam");
    my $sorted_bam_file_index    = "${sorted_bam_file}.bai";

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

    my $final_file =  $sorted_bam_file;
    my $check_file = $sortByCoordinate? $sorted_bam_file_index :$final_file;

    my $sample_file_0 = $sample_files[0];
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $check_file, $init_command, 0, $sample_file_0 );

    print $sh "
if [[ ! -s $result_dir/$final_file ]]; then
  \$MYCMD ./$pbs_name 
fi
";

    print $pbs "
rm -f ${unsorted_bam_file}.failed
";

    my $localized_files = [];
    @sample_files = @{$self->localize_files_in_tmp_folder($pbs, \@sample_files, $localized_files)};
    $sample_file_0 = $sample_files[0];
    my $sample_files_str = ( scalar(@sample_files) == 2 ) ? "\"" . $sample_file_0 . "\" \"" . $sample_files[1] . "\"" : "\"" . $sample_file_0 . "\"";

    print $pbs "

if [[ ! -e ${unsorted_bam_file}.succeed ]]; then
  echo bwa_mem=`date`
  bwa mem $option $rg $bwa_index $sample_files_str | samtools view -bS -o $unsorted_bam_file

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch ${unsorted_bam_file}.succeed
    rm -f ${unsorted_bam_file}.failed
  else
    rm -f $unsorted_bam_file ${unsorted_bam_file}.succeed
    touch ${unsorted_bam_file}.failed
  fi

  bwa 2>\&1 | grep Version | cut -d ' ' -f2 | cut -d '-' -f1 | awk '{print \"bwa,v\"\$1}' > ${sample_name}.bwa.version
fi
";

    if($output_unmapped_fastq){
      #http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/1040-2/
      print $pbs "

if [[ -e ${unsorted_bam_file}.succeed ]]; then
  echo bwa_unmapped=`date`
  samtools view -u  -f 4 -F 264 $unsorted_bam_file  > ${sample_name}.tmps1.bam
  samtools view -u -f 8 -F 260 $unsorted_bam_file  > ${sample_name}.tmps2.bam
  samtools view -u -f 12 -F 256 $unsorted_bam_file > ${sample_name}.tmps3.bam

  samtools merge -u - ${sample_name}.tmps[123].bam | samtools sort -n - -o ${sample_name}.unmapped.bam
  status=\$?
  if [[ \$status -eq 0 ]]; then
    rm ${sample_name}.tmps[123].bam
    samtools fastq -1 ${sample_name}.unmapped.1.fq.gz -2 ${sample_name}.unmapped.2.fq.gz ${sample_name}.unmapped.bam 
    status=\$?
    if [[ \$status -eq 0 ]]; then
      rm ${sample_name}.unmapped.bam
    fi
  fi
fi

";
    }

    if ($bwa_only) {
      $self->close_pbs( $pbs, $pbs_file );
      next;
    }

    my $bam_stat_condition = "";
    if($do_bam_stat){
      print $pbs "
if [[ -e ${unsorted_bam_file}.succeed && ((1 -eq \$1) || (! -s $bam_stat)) ]]; then
  echo bamStat=`date` 
  python3 $py_script -i $unsorted_bam_file -o $bam_stat
fi
";
      $bam_stat_condition = "&& -s $bam_stat";
    }

    my $rmlist = "";

    if ($alignmentOnly) { #only alignment, no sort or other works. For UMI pipeline
      print $pbs "
if [[ -e ${unsorted_bam_file}.succeed ]]; then
  echo flagstat=`date` 
  samtools flagstat $unsorted_bam_file > ${unsorted_bam_file}.stat 
fi
";
      if ($rmlist ne "") {
        print $pbs "
if [[ -e ${unsorted_bam_file}.succeed $bam_stat_condition ]]; then
  rm $rmlist
fi
";
      }

      $self->clean_temp_files($pbs, $localized_files);

      $self->close_pbs( $pbs, $pbs_file );
      next;
    }

    if ($sortByCoordinate) {
      my $chromosome_grep_command = getChromosomeFilterCommand( $sorted_bam_file, $chromosome_grep_pattern );
      if($use_sambamba){
        print $pbs "    
if [[ (-e ${unsorted_bam_file}.succeed) && ((1 -eq \$1) || (! -s $sorted_bam_file)) ]]; then
  echo sort_bam=`date`
  sambamba sort -m $sort_memory $sambamba_sort_thread --tmpdir bwa_${sample_name} --sort-picard -o $sorted_bam_file $unsorted_bam_file

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch ${sorted_bam_file}.succeed
    rm -rf bwa_${sample_name} ${sorted_bam_file}.failed
  else
    rm -rf bwa_${sample_name} $sorted_bam_file ${sorted_bam_file}.succeed
    touch ${sorted_bam_file}.failed
  fi
fi

if [[ (-e ${sorted_bam_file}.succeed) && ((1 -eq \$1) || (! -s ${sorted_bam_file}.bai)) ]]; then
  echo index_bam=`date`
  sambamba index $sambamba_index_thread $sorted_bam_file 
  samtools idxstats $sorted_bam_file > ${sorted_bam_file}.chromosome.count
fi

$chromosome_grep_command
";
      }else{
        print $pbs "    
if [[ (-e ${unsorted_bam_file}.succeed) && ((1 -eq \$1) || (! -s $sorted_bam_file)) ]]; then
  echo sort_bam=`date`
  samtools sort -m $sort_memory $samtools_sort_thread -T bwa_${sample_name} -o $sorted_bam_file $unsorted_bam_file

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch ${sorted_bam_file}.succeed
    rm -rf bwa_${sample_name} ${sorted_bam_file}.failed
  else
    rm -rf bwa_${sample_name} $sorted_bam_file ${sorted_bam_file}.succeed
    touch ${sorted_bam_file}.failed
  fi
fi

if [[ (-e ${sorted_bam_file}.succeed) && ((1 -eq \$1) || (! -s ${sorted_bam_file}.bai)) ]]; then
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
if [[ (-e ${unsorted_bam_file}.succeed) && ((1 -eq \$1) || (! -s $sorted_bam_file)) ]]; then
  echo sort_bam=`date`
  sambamba sort -m $sort_memory $sambamba_sort_thread --tmpdir bwa_${sample_name} --sort-picard -o $sorted_bam_file $unsorted_bam_file && touch ${sorted_bam_file}.succeed
  rm -rf bwa_${sample_name}
  if [[ ! -e ${sorted_bam_file}.succeed ]]; then
    rm $sorted_bam_file
  fi
fi
";
      }else{
        print $pbs "    
if [[ (-e ${unsorted_bam_file}.succeed) && ((1 -eq \$1) || (! -s $sorted_bam_file)) ]]; then
  echo sort_bam=`date`
  samtools sort -m $sort_memory $samtools_sort_thread -T bwa_${sample_name} -n -o $sorted_bam_file $unsorted_bam_file && touch ${sorted_bam_file}.succeed
  rm -rf bwa_${sample_name}
  if [[ ! -e ${sorted_bam_file}.succeed ]]; then
    rm -f $sorted_bam_file
  fi
fi
";
      }
    }
    $rmlist = $rmlist . " " . $unsorted_bam_file;

    print $pbs "
if [[ -e ${sorted_bam_file}.succeed && -s $check_file ]]; then
  echo flagstat=`date` 
  samtools flagstat $final_file > ${final_file}.stat 
fi

rm -rf bwa_${sample_name} 
";

  if ($rmlist ne "") {
    print $pbs "
if [[ -e ${sorted_bam_file}.succeed $bam_stat_condition ]]; then
  rm -f $rmlist
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
  my %raw_files        = %{ get_raw_files( $config, $section ) };
  my $output_unmapped_fastq = get_option( $config, $section, "output_unmapped_fastq", 0 );

  my $do_bam_stat = get_option( $config, $section, "do_bam_stat",     1 );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $final_file;
    my $alignmentOnly        = get_option( $config, $section, "alignmentOnly",      0 );
    if ($alignmentOnly) { #only alignment, no sort or other works. For UMI pipeline
      $final_file=$sample_name.".unsorted.bam";
    } else { #all other bwa tasks, with sort or mark_duplicates
      my $sorted_bam_file    = $sample_name . ($sortByCoordinate? ".sortedByCoord.bam":".sortedByQuery.bam");
      $final_file = $sorted_bam_file;
    }
    $final_file = "${result_dir}/$final_file";
    my @result_files = ();
    push( @result_files, $final_file );
    push( @result_files, $final_file . ".stat" );
    if($sortByCoordinate){
      push( @result_files, $final_file . ".chromosome.count" );
    }
    if($output_unmapped_fastq){
      push( @result_files, "${result_dir}/${sample_name}.unmapped.1.fq.gz" );
      push( @result_files, "${result_dir}/${sample_name}.unmapped.2.fq.gz" );
    }
    if($do_bam_stat){
      push( @result_files, "${result_dir}/${sample_name}.bamstat" );
    }
    push( @result_files, "${result_dir}/${sample_name}.bwa.version" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
