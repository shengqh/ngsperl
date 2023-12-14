#!/usr/bin/perl
package Alignment::Abismal;

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
  $self->{_suffix} = "_abismal";
  $self->{_docker_prefix} = "dnmtools_";
  $self->{_docker_shell} = "sh";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $selfname = $self->{_name};

  my $dnmtools_command = get_option($config, $section, "dnmtools_command", "dnmtools");

  print("dnmtools_command = " . $dnmtools_command);

  my $preseq_command = get_option($config, $section, "preseq_command", "preseq");

  my $abismal_index = get_option($config, $section, "abismal_index");
  my $chr_fasta = get_option_file($config, $section, "chr_fasta");
  my $addqual_pythonFile = get_option_file($config, $section, "addqual_pythonFile");
  my $interval_list = get_option_file($config, $section, "interval_list");

  my $picard = $config->{$section}{picard};

  my %raw_files = %{ get_raw_files( $config, $section ) };
  
  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sample_files_str = ( scalar(@sample_files) == 2 ) ? "-1 " . $sample_files[0] . " -2 " . $sample_files[1]: "-1 " . $sample_files[0];

    my $raw_bam = $sample_name . ".raw.bam";
    my $sorted_bam = $sample_name . ".sorted.bam";
    my $result_file_addqual     = $sample_name . ".sorted.addqual.sam";
    my $result_file_addqual_bam     = $sample_name . ".sorted.addqual.bam";
    my $result_file_addqual_bai     = $sample_name . ".sorted.addqual.bam.bai";

    my $uniq_bam = $sample_name . ".uniq.bam";

    my $map_stat                 = $sample_name . ".mapstats";

    my $hs_metrics               = $sample_name . "_hs_metrics.txt";
    my $per_target_coverage      = $sample_name . "_hs_metrics_per_interval.txt";
    my $rrbs_metrics             = $sample_name . ".rrbs_summary_metrics";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $rmlist = "";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, "${sample_name}.rrbs_summary_metrics" );

		foreach my $sampleFile (@sample_files) {
			print $pbs "
if [ ! -s $sampleFile ]; then
  echo missing input file $sampleFile!
  exit 1;
fi
";
    }

    $sample_files_str = ( scalar(@sample_files) == 2 ) ? " " . $sample_files[0] . " " . $sample_files[1]: " " . $sample_files[0];
	
    print $pbs "
if [[ ! -s $raw_bam ]]; then
  echo abismal=`date`
  $dnmtools_command abismal -B -v -t $thread -i $abismal_index -s $map_stat $sample_files_str -o $raw_bam
fi

if [[ -s $raw_bam && ! -s ${sample_name}.uniq.bam.dupstats ]]; then
  if [[ ! -s ${sample_name}.formatted.sorted.bam ]]; then
    #dnmtools format requires unsorted bam file
    #-F: This option forces the format command to process paired-end reads even if it is unable to detect mates. 
    echo dnmtools format=`date`
    $dnmtools_command format -t $thread -B -F -v -f abismal -stdout $raw_bam | samtools sort -@ $thread -o ${sample_name}.formatted.sorted.bam
  fi

  echo dnmtools uniq=`date`
  $dnmtools_command uniq -t $thread -B -v -S ${sample_name}.uniq.bam.dupstats ${sample_name}.formatted.sorted.bam ${sample_name}.uniq.bam
  samtools index -@ $thread ${sample_name}.uniq.bam

  if [[ -s ${sample_name}.uniq.bam.dupstats ]]; then
    rm ${sample_name}.formatted.sorted.bam
  fi
fi

if [[ -s $raw_bam && ! -s ${sample_name}.rrbs_summary_metrics ]]; then
  if [[ ! -s $result_file_addqual_bai ]]; then
    if [[ ! -s $sorted_bam ]]; then
      echo sort=`date`
      samtools sort -@ $thread -O BAM -o $sorted_bam $raw_bam
    fi

    echo add_qual =`date`
    python3 $addqual_pythonFile $sorted_bam $result_file_addqual_bam
    samtools index $result_file_addqual_bam $result_file_addqual_bai
  fi

  if [[ -s $result_file_addqual_bai ]]; then
   if [[ ! -s ${sample_name}.c_curve ]]; then
      echo preseq c_curve=`date`
      $preseq_command c_curve -B $result_file_addqual_bam > ${sample_name}.c_curve
    fi

    if [[ ! -s ${sample_name}.lc_extrap ]]; then
      echo preseq lc_extrap=`date`
      $preseq_command lc_extrap -B -D $result_file_addqual_bam > ${sample_name}.lc_extrap
    fi

    if [[ ! -s $hs_metrics ]]; then
      echo picard CollectHsMetrics =`date`
      java -jar $picard CollectHsMetrics \\
        I=$result_file_addqual_bam \\
        O=$hs_metrics \\
        R=$chr_fasta \\
        BAIT_INTERVALS=$interval_list \\
        TARGET_INTERVALS=$interval_list
    fi

    if [[ ! -s $rrbs_metrics ]]; then
      echo picard CollectRrbsMetrics =`date`
      java -jar $picard CollectRrbsMetrics \\
        I=$result_file_addqual_bam \\
        M=$sample_name \\
        R=$chr_fasta

      gzip ${sample_name}.rrbs_detail_metrics
    fi

    if [[ -s ${sample_name}.rrbs_summary_metrics ]]; then
      rm -f $result_file_addqual_bam $result_file_addqual_bai
    fi
  fi
fi

if [[ -s ${sample_name}.rrbs_summary_metrics && -s ${sample_name}.uniq.bam.dupstats ]]; then
  rm -f $raw_bam
fi

";

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all abismal tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    push( @result_files, "${result_dir}/${sample_name}.uniq.bam" );
    push( @result_files, "${result_dir}/${sample_name}.uniq.bam.dupstats" );
    push( @result_files, "${result_dir}/${sample_name}.mapstats" );
    push( @result_files, "${result_dir}/${sample_name}_hs_metrics.txt" );
    push( @result_files, "${result_dir}/${sample_name}.rrbs_summary_metrics" );
    
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
