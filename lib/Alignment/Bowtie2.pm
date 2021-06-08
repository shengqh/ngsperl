#!/usr/bin/perl
package Alignment::Bowtie2;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use Alignment::AbstractBowtie;
use Alignment::AlignmentUtils;

our @ISA = qw(Alignment::AbstractBowtie);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_bt2";
  $self->{_use_tmp_folder} = 1;
  $self->{log_result} = 1;
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $bowtie2_index = $config->{$section}{bowtie2_index} or die "define ${section}::bowtie2_index first";
  my $chromosome_grep_pattern = get_option( $config, $section, "chromosome_grep_pattern", "" );
  my $outputToSameFolder      = get_option( $config, $section, "output_to_same_folder",   1 );

  my $mark_duplicates = hasMarkDuplicate( $config->{$section} );
  my $picard_jar      = "";
  if ($mark_duplicates) {
    $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1, not $self->using_docker() );
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct), "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sam_file     = $sample_name . ".sam";
    my $bam_file     = $sample_name . ".bam";
    my $log_file     = $sample_name . ".log";
    my $sort_log_file     = $sample_name . ".sort.log";

    my $pbs_name = $self->pbs_name($sample_name);
    my $pbs_file = $pbs_dir . "/$pbs_name";
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir  = $outputToSameFolder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    my $final_file = $mark_duplicates ? $sample_name . ".rmdup.bam" : $bam_file;
    my $chromosome_grep_command = getChromosomeFilterCommand( $bam_file, $chromosome_grep_pattern );


    print $sh "
if [[ ! -s $result_dir/${final_file}.bai ]]; then
  \$MYCMD ./$pbs_name 
fi
";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    my $localized_files = [];
    @sample_files = @{$self->localize_files_in_tmp_folder($pbs, \@sample_files, $localized_files)};

    my $indent = "";
    my $tag    = "--sam-RG ID:$sample_name --sam-RG LB:$sample_name --sam-RG SM:$sample_name --sam-RG PL:ILLUMINA";

    my $input = "";
    if ( scalar(@sample_files) == 2 ) {    #pairend data
      $input = "-1 " . $sample_files[0] . " -2 " . $sample_files[1];
    }
    else {
      my $fastqs = join( ',', @sample_files );
      $input = "-U " . $fastqs;
    }
    my $bowtie2_aln_command = "bowtie2 -p $thread $option -x $bowtie2_index $input $tag -S $sam_file 2> $log_file";

    my $index_command = get_index_command( $bam_file, $indent );
    my $stat_command = get_stat_command( $bam_file, $indent );

    print $pbs "

if [ ! -s $sam_file ]; then
  $bowtie2_aln_command
  bowtie2 --version | grep bowtie2 | grep version | cut -d ' ' -f3 | awk '{print \"bowtie2,v\"\$1}' > ${final_file}.version
fi

if [ -s $log_file ]; then
  isSucceed=\$(cat $log_file | grep -c \"overall alignment rate\")
  if [ \$isSucceed ]; then
    echo alignment succeed, sorting sam to bam ...
    samtools view -Shu -F 256 $sam_file | samtools sort -o $bam_file -T $sample_name - 
    if [[ -s $bam_file ]]; then
      echo index bam ...
      samtools index $bam_file 
      $chromosome_grep_command
      samtools idxstats $bam_file > ${bam_file}.chromosome.count
      rm $sam_file
    fi
  fi
fi
";

    if ($mark_duplicates) {
      print $pbs "
if [ -s $bam_file ]; then
  echo RemoveDuplicate=`date` 
  java -jar $picard_jar MarkDuplicates I=$bam_file O=$final_file ASSUME_SORTED=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT M=${final_file}.metrics
  if [ -s $final_file ]; then
    rm $bam_file ${bam_file}.bai ${bam_file}.chromosome.count
    samtools index $final_file 
    samtools idxstats $final_file > ${final_file}.chromosome.count
  fi
fi
";
    } else { #prepare to do flagstat if not mark_duplicates
      $final_file=$bam_file;
    }

    print $pbs "
if [[ (-s $final_file) ]]; then
  echo flagstat=`date` 
  samtools flagstat $final_file > ${final_file}.stat 
fi
";


    $self->clean_temp_files($pbs, $localized_files);

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

1;
