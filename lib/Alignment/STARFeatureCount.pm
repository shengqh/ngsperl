#!/usr/bin/perl
package Alignment::STARFeatureCount;

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
use List::Util qw[min max];

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_sf";
  $self->{_use_tmp_folder} = 1;
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  #  if ( $option !~ /outSAMprimaryFlag/ ) {
  #    $option = $option . "  --outSAMprimaryFlag AllBestScore";
  #  }
  #

  my $py_script = dirname(__FILE__) . "/bamStat.py";
  if ( !-e $py_script ) {
    die "File not found : " . $py_script;
  }

  my ($sort_memory, $isMB) = getMemoryPerThread($memory, $thread);
  if ($isMB) {
    $sort_memory = $sort_memory . "M";
  }else{
    $sort_memory = $sort_memory . "G";
  }

  my $chromosome_grep_pattern = get_option( $config, $section, "chromosome_grep_pattern", "" );
  my $star                    = get_option( $config, $section, "star_location",           "STAR" );

  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder", 1 );
  my $output_sort_by_coordinate = getSortByCoordinate( $config, $section, 1 );
  my $delete_star_featureCount_bam  = get_option( $config, $section, "delete_star_featureCount_bam", 0 );

  if ($delete_star_featureCount_bam) {
    $output_sort_by_coordinate = 0;
  }

  my $output_unsorted = get_option( $config, $section, "output_unsorted", 0 );
  if ( !$output_sort_by_coordinate && !$output_unsorted ) {
    $output_unsorted = 1;
  }
  my $output_format = "--outSAMtype BAM";

  #always output unsorted
  $output_format = $output_format . " Unsorted";

  my $star_index;
  if ( defined $config->{$section}{"genome_dir"} ) {
    $star_index = parse_param_file( $config, $section, "genome_dir", 1 );
  }
  else {
    $star_index = parse_param_file( $config, $section, "star_index", 1 );
  }

  my $featureCountOption = get_option( $config, $section, "featureCount_option", "" );
  my $ispaired = get_is_paired_end_option( $config, $section );
  if ($ispaired) {
    $featureCountOption = $featureCountOption . " -p --countReadPairs";
  }

  if ( $featureCountOption !~ /-g/ ) {
    $featureCountOption = $featureCountOption . " -g gene_id";
  }

  if ( $featureCountOption !~ /-t/ ) {
    $featureCountOption = $featureCountOption . " -t exon";
  }

  my $gffFile = parse_param_file( $config, $section, "gff_file", 1 );

  my %fqFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %fqFiles ) {
    my @sample_files = @{ $fqFiles{$sample_name} };
    my $sample_file_1 = $sample_files[0];

    my $uncompress = ( $sample_file_1 =~ /.gz$/ ) ? " --readFilesCommand zcat" : "";

    my $samples = join( " ", @sample_files );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );
    my $rgline = "ID:$sample_name SM:$sample_name LB:$sample_name PL:ILLUMINA PU:ILLUMINA";

    my $unsorted = $sample_name . "_Aligned.out.bam";
    my $bam_stat = $sample_name . ".bamstat";

    my $final_bam = $output_sort_by_coordinate ? $sample_name . "_Aligned.sortedByCoord.out.bam" : $unsorted;

    my $final_file = "${sample_name}.count";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    my $chromosome_grep_command = $output_sort_by_coordinate ? getChromosomeFilterCommand( $final_bam, $chromosome_grep_pattern ) : "";

    print $pbs "

if [[ ! -s $final_bam && -s $sample_file_1 ]]; then
  echo performing star ...
  $star $option --outSAMattrRGline $rgline --runThreadN $thread --genomeDir $star_index --readFilesIn $samples $uncompress --outFileNamePrefix ${sample_name}_ $output_format
  $star --version | awk '{print \"STAR,v\"\$1}' > ${final_file}.star.version
  rm -rf ${sample_name}__STARgenome ${sample_name}__STARpass1 ${sample_name}_Log.progress.out
fi
";

    print $pbs "
if [ -s $unsorted ]; then
  echo bamStat=`date` 
  python3 $py_script -i $unsorted -o $bam_stat
fi
";


  if ($output_sort_by_coordinate) {
    print $pbs "
if [ ! -s ${final_bam} ]; then
  samtools sort -m $sort_memory -T ${sample_name} -t $thread -o $final_bam $unsorted
  samtools index $final_bam
  samtools idxstats $final_bam > ${final_bam}.chromosome.count
fi  
";
  }

print $pbs "

if [ -s $final_bam ]; then
  $chromosome_grep_command
  
  if [ ! -s ${sample_name}.splicing.bed ]; then
    awk {'if(\$4==\"2\") print \"\"\$1\"\\t\"\$2-\$9-1\"\\t\"\$3+\$9\"\\tJUNC000\"NR\"\\t\"\$7+\$8\"\\t-\\t\"\$2-\$9-1\"\\t\"\$3+\$9\"\\t255,0,0\\t2\\t\"\$9\",\"\$9\"\\t\",\"0,\"\$3-\$2+\$9+1; else if(\$4==\"1\") print \"\"\$1\"\\t\"\$2-\$9-1\"\\t\"\$3+\$9\"\\tJUNC000\"NR\"\\t\"\$7+\$8\"\\t+\\t\"\$2-\$9-1\"\\t\"\$3+\$9\"\\t0,0,255\\t2\\t\"\$9\",\"\$9\"\\t\",\"0,\"\$3-\$2+\$9+1'} ${sample_name}_SJ.out.tab > ${sample_name}.splicing.bed
  fi
fi

if [ -s $unsorted ]; then
  echo performing featureCounts ...
  featureCounts $featureCountOption -T $thread -a $gffFile -o $final_file $unsorted
  featureCounts -v 2>\&1 | grep featureCounts | cut -d ' ' -f2 | awk '{print \"featureCounts,\"\$1}' > ${final_file}.featureCounts.version
fi 
";

    if ( !$output_unsorted ) {
      print $pbs "
if [ -s $final_file && -s $bam_stat ]; then
  rm $unsorted 
fi
";
    }

    if ($delete_star_featureCount_bam){
      print $pbs "
if [ -s $final_file ]; then
  rm ${sample_name}_Aligned.out.bam ${sample_name}_Aligned.out.bam.bai
  rm ${sample_name}_SJ.out.tab ${sample_name}.splicing.bed
fi
";
    }

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
  my $output_sort_by_coordinate = getSortByCoordinate( $config, $section );
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder", 1 );
  my $delete_star_featureCount_bam  = get_option( $config, $section, "delete_star_featureCount_bam", 0 );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $cur_dir = $output_to_same_folder ? $result_dir : $result_dir . "/$sample_name";

    my @result_files              = ();
    my $output_sort_by_coordinate = get_option( $config, $section, "output_sort_by_coordinate", 1 );
    my $output_unsorted           = get_option( $config, $section, "output_unsorted", 0 );
    if ( !$output_sort_by_coordinate && !$output_unsorted ) {
      $output_unsorted = 1;
    }

    push( @result_files, "$cur_dir/${sample_name}.count.summary" );
    push( @result_files, "$cur_dir/${sample_name}_Log.final.out" );
    push( @result_files, "$cur_dir/${sample_name}.bamstat" );
    if (!$delete_star_featureCount_bam) {
      if ($output_sort_by_coordinate) {
        push( @result_files, "$cur_dir/${sample_name}_Aligned.sortedByCoord.out.bam" );
        push( @result_files, "$cur_dir/${sample_name}_Aligned.sortedByCoord.out.bam.chromosome.count" );
      }
      if ($output_unsorted) {
        push( @result_files, "$cur_dir/${sample_name}_Aligned.out.bam" );
      }
    }
    push( @result_files, "$cur_dir/${sample_name}.count" );
    push( @result_files, "$cur_dir/${sample_name}.count.star.version" );
    push( @result_files, "$cur_dir/${sample_name}.count.featureCounts.version" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
