#!/usr/bin/perl
package Alignment::AbismalFilterIntervals;

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

  my $picard = get_option_file($config, $section, "picard");
  my $gatk = get_option_file($config, $section, "gatk");

  my %raw_files = %{ get_raw_files( $config, $section ) };
  
  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sample_files_str = ( scalar(@sample_files) == 2 ) ? "-1 " . $sample_files[0] . " -2 " . $sample_files[1]: "-1 " . $sample_files[0];

    my $raw_bam = $sample_name . ".raw.bam";

    my $map_stat                 = $sample_name . ".mapstats";

    my $hs_metrics               = $sample_name . "_hs_metrics.txt";
    my $rrbs_metrics             = $sample_name . ".rrbs_summary_metrics";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $rmlist = "";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $rrbs_metrics );

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
# Mapping with abismal    
if [[ ! -s $raw_bam || ! -e ${sample_name}.abismal.succeed ]]; then
  rm -f ${sample_name}.abismal.failed ${sample_name}.abismal.succeed $raw_bam

  echo abismal=`date`
  $dnmtools_command abismal -B -v -t $thread -i $abismal_index -s $map_stat $sample_files_str -o $raw_bam
  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch ${sample_name}.abismal.succeed
  else
    touch ${sample_name}.abismal.failed
    rm -f $raw_bam
    exit \$status
  fi
fi

# Performed CollectHsMetrics on raw data level.
if [[ ! -s $hs_metrics ]]; then
  if [[ ! -s $sample_name.sorted.bam || ! -e ${sample_name}.sort.succeed ]]; then
    rm -f $sample_name.sorted.bam ${sample_name}.sort.failed ${sample_name}.sort.succeed

    echo sort=`date`
    samtools sort --threads $thread -O BAM -o $sample_name.sorted.bam $raw_bam

    status=\$?
    if [[ \$status -eq 0 ]]; then
      touch ${sample_name}.sort.succeed
    else
      touch ${sample_name}.sort.failed
      rm -f $sample_name.sorted.bam
      exit \$status
    fi
  fi

  # Missing the qual data would cause error. We will have to add the qual data first.
  if [[ ! -s $sample_name.sorted.addqual.bam.bai || ! -e ${sample_name}.sorted.addqual.succeed ]]; then
    rm -f $sample_name.sorted.addqual.bam $sample_name.sorted.addqual.bam.bai ${sample_name}.sorted.addqual.succeed ${sample_name}.sorted.addqual.failed

    echo add_qual =`date`
    python3 $addqual_pythonFile $sample_name.sorted.bam $sample_name.sorted.addqual.bam

    status=\$?
    if [[ \$status -eq 0 ]]; then
      touch ${sample_name}.sorted.addqual.succeed
    else
      touch ${sample_name}.sorted.addqual.failed
      rm -f $sample_name.sorted.addqual.bam
      exit \$status
    fi    

    echo samtools index=`date`
    samtools index --threads $thread $sample_name.sorted.addqual.bam
  fi

  rm -f ${sample_name}.CollectHsMetrics.failed ${sample_name}.CollectHsMetrics.succeed

  echo picard CollectHsMetrics =`date`
  java -jar $picard CollectHsMetrics \\
    I=$sample_name.sorted.addqual.bam \\
    O=$hs_metrics \\
    R=$chr_fasta \\
    BAIT_INTERVALS=$interval_list \\
    TARGET_INTERVALS=$interval_list

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch ${sample_name}.CollectHsMetrics.succeed
    #rm -f $sample_name.sorted.addqual.bam $sample_name.sorted.addqual.bam.bai $sample_name.sorted.bam $sample_name.sorted.bam.bai
  else
    touch ${sample_name}.CollectHsMetrics.failed
    exit \$status
  fi
fi

# Prepare final bam file for methylation calling
if [[ ! -s ${sample_name}.intervals.uniq.addqual.bam.bai || ! -e ${sample_name}.intervals.uniq.addqual.succeed ]]; then
  if [[ ! -s ${sample_name}.formatted.sorted.bam.bai || ! -s ${sample_name}.formatted.sorted.succeed ]]; then
    rm -f ${sample_name}.formatted.sorted.bam ${sample_name}.formatted.sorted.bam.bai ${sample_name}.formatted.sorted.bam.flagstat ${sample_name}.formatted.sorted.succeed ${sample_name}.formatted.sorted.failed

    #dnmtools format requires unsorted bam file
    #-F: This option forces the format command to process paired-end reads even if it is unable to detect mates. 
    echo dnmtools format=`date`
    $dnmtools_command format -t $thread -B -F -v -f abismal -stdout $raw_bam | samtools sort --threads $thread -o ${sample_name}.formatted.sorted.bam
    
    status=\$?
    if [[ \$status -eq 0 ]]; then
      touch ${sample_name}.formatted.sorted.succeed
    else
      touch ${sample_name}.formatted.sorted.failed
      rm -f ${sample_name}.formatted.sorted.bam
      exit \$status
    fi

    echo samtools index=`date`
    samtools index --threads $thread ${sample_name}.formatted.sorted.bam

    echo samtools flagstat=`date`
    samtools flagstat ${sample_name}.formatted.sorted.bam > ${sample_name}.formatted.sorted.bam.flagstat
  fi

  # Keep the reads only in the target regions
  if [[ ! -s ${sample_name}.intervals.bam.bai || ! -e ${sample_name}.intervals.succeed ]]; then
    rm -f ${sample_name}.intervals.bam ${sample_name}.intervals.bam.bai ${sample_name}.intervals.bam.flagstat ${sample_name}.intervals.succeed ${sample_name}.intervals.failed

    # Disable WellformedReadFilter to keep all reads, otherwise, it would remove all reads.
    echo gatk filter intervals PrintReads=`date`
    java -jar $gatk PrintReads \\
      --disable-read-filter WellformedReadFilter \\
      --input ${sample_name}.formatted.sorted.bam \\
      --output ${sample_name}.intervals.bam \\
      --intervals $interval_list \\
      --reference $chr_fasta

    status=\$?
    if [[ \$status -eq 0 ]]; then
      touch ${sample_name}.intervals.succeed
    else
      touch ${sample_name}.intervals.failed
      rm -f ${sample_name}.intervals.bam
      exit \$status
    fi

    echo samtools flagstat=`date`
    samtools flagstat ${sample_name}.intervals.bam > ${sample_name}.intervals.bam.flagstat
  fi

  # Remove duplicates
  if [[ ! -s ${sample_name}.intervals.uniq.bam || ! -e ${sample_name}.intervals.uniq.succeed ]]; then
    rm -f ${sample_name}.intervals.uniq.bam ${sample_name}.intervals.uniq.failed ${sample_name}.intervals.uniq.succeed

    echo dnmtools uniq=`date`
    $dnmtools_command uniq -t $thread -B -v -S ${sample_name}.intervals.uniq.bam.dupstats ${sample_name}.intervals.bam ${sample_name}.intervals.uniq.bam

    status=\$?
    if [[ \$status -eq 0 ]]; then
      touch ${sample_name}.intervals.uniq.succeed
    else
      touch ${sample_name}.intervals.uniq.failed
      rm -f ${sample_name}.intervals.uniq.bam
      exit \$status
    fi
  fi  

  # add qual to uniq bam file in case it would be used in downstream analysis

  rm -f ${sample_name}.intervals.uniq.addqual.bam ${sample_name}.intervals.uniq.addqual.bam.bai ${sample_name}.intervals.uniq.addqual.failed ${sample_name}.intervals.uniq.addqual.succeed

  echo add_qual =`date`
  python3 $addqual_pythonFile ${sample_name}.intervals.uniq.bam ${sample_name}.intervals.uniq.addqual.bam

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch ${sample_name}.intervals.uniq.addqual.succeed
  else
    touch ${sample_name}.intervals.uniq.addqual.failed
    rm -f ${sample_name}.intervals.uniq.addqual.bam
    exit \$status
  fi

  echo samtools index=`date`
  samtools index --threads $thread ${sample_name}.intervals.uniq.addqual.bam

  echo samtools flagstat=`date`
  samtools flagstat ${sample_name}.intervals.uniq.addqual.bam > ${sample_name}.intervals.uniq.addqual.bam.flagstat

  if [[ -s ${sample_name}.intervals.uniq.addqual.bam.bai && -e ${sample_name}.intervals.uniq.addqual.succeed ]]; then
    #rm -f $raw_bam ${sample_name}.formatted.sorted.bam ${sample_name}.formatted.sorted.bam.bai ${sample_name}.intervals.bam ${sample_name}.intervals.bam.bai ${sample_name}.intervals.uniq.bam ${sample_name}.intervals.uniq.bam.bai
  fi
fi

# For rrbs summary metrics, we will perform on the filtered uniq bam file since we want to know 
# the QC result of the final usable data.
if [[ ! -s ${sample_name}.rrbs_summary_metrics ]]; then
  if [[ ! -s ${sample_name}.c_curve ]]; then
    echo preseq c_curve=`date`
    $preseq_command c_curve -B ${sample_name}.intervals.uniq.addqual.bam > ${sample_name}.c_curve
  fi

  if [[ ! -s ${sample_name}.lc_extrap ]]; then
    echo preseq lc_extrap=`date`
    $preseq_command lc_extrap -B -D ${sample_name}.intervals.uniq.addqual.bam > ${sample_name}.lc_extrap
  fi

  if [[ ! -s $rrbs_metrics ]]; then
    echo picard CollectRrbsMetrics =`date`
    java -jar $picard CollectRrbsMetrics \\
      I=${sample_name}.intervals.uniq.addqual.bam \\
      M=$sample_name \\
      R=$chr_fasta

    gzip ${sample_name}.rrbs_detail_metrics
  fi
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
    push( @result_files, "${result_dir}/${sample_name}.intervals.uniq.addqual.bam" );
    push( @result_files, "${result_dir}/${sample_name}.intervals.uniq.addqual.bam.dupstats" );
    push( @result_files, "${result_dir}/${sample_name}.mapstats" );
    push( @result_files, "${result_dir}/${sample_name}_hs_metrics.txt" );
    push( @result_files, "${result_dir}/${sample_name}.rrbs_summary_metrics" );
    
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
