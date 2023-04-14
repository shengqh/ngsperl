#!/usr/bin/perl
package Alignment::abismal;

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
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $selfname = $self->{_name};

  my $abismal_index = $config->{$section}{abismal_index};

  my $addqual_perlFile = $config->{$section}{addqual_perlFile};

  my $picard = $config->{$section}{picard};
  my $chrDir=$config->{$section}{chr_dir};
  my $interval_list=$config->{$section}{interval_list};
  if ( !defined $chrDir ) {
    die "define ${section}::chr_dir first";
  }
  if ( !defined $abismal_index ) {
    die "define ${section}::abismal_index first";
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };
  
  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sample_files_str = ( scalar(@sample_files) == 2 ) ? "-1 " . $sample_files[0] . " -2 " . $sample_files[1]: "-1 " . $sample_files[0];

    my $result_file_raw          = $sample_name . ".raw.sam";
    my $result_file              = $sample_name . ".sorted.sam";
    my $result_file_bam          = $sample_name . ".sorted.bam";
    my $result_file_addqual     = $sample_name . ".sorted.addqual.sam";
    my $result_file_addqual_bam     = $sample_name . ".sorted.addqual.bam";
    my $result_file_addqual_bai     = $sample_name . ".sorted.addqual.bam.bai";
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
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $result_file );

#			print $pbs "
#if [ ! -s $result_file ]; then
#  echo zcat=`date`
#";
	if ($sample_files[0]=~/\.gz$/) { #.gz fle, need to zcat
		my @sample_files_unzip=();
		foreach my $sampleFile (@sample_files) {
			my $sampleFileUnzip = basename(change_extension( $sampleFile, "" ));
			push @sample_files_unzip,$sampleFileUnzip;
			print $pbs "
  if [ ! -s $sampleFileUnzip ]; then
      echo zcat=`date`
      zcat $sampleFile > $sampleFileUnzip
  fi
";
			$rmlist=$rmlist. " $sampleFileUnzip"
		}
#					print $pbs "
#fi
#";
		$sample_files_str = ( scalar(@sample_files_unzip) == 2 ) ? " " . $sample_files_unzip[0] . " " . $sample_files_unzip[1]: " " . $sample_files_unzip[0];
	} else { #NOT .gz fle
	  $sample_files_str = ( scalar(@sample_files) == 2 ) ? " " . $sample_files[0] . " " . $sample_files[1]: " " . $sample_files[0];
	}
	
    print $pbs "
if [[ ! -s $result_file_raw && ! -s $result_file_bam ]]; then
  echo abismal=`date`
  dnmtools abismal -t $thread -i $abismal_index $sample_files_str -s $map_stat -o $result_file_raw
  samtools sort -@ $thread -O BAM -o $result_file_bam $result_file_raw
fi
";
    $rmlist=$rmlist ." $result_file_raw";
    if ($rmlist ne "") {
      print $pbs "
if [[ -s $result_file_bam ]]; then
  rm $rmlist
fi
";
    }

    print $pbs "
if [[ ( -s $result_file_bam ) && ( ! -s $result_file_addqual_bam) ]]; then
  echo add_qual =`date`
  samtools view -h -O SAM $result_file_bam | perl $addqual_perlFile - $result_file_addqual
  samtools view -b -o $result_file_addqual_bam $result_file_addqual
  samtools index $result_file_addqual_bam $result_file_addqual_bai
fi
";
    $rmlist=" $result_file_addqual";

    if ($rmlist ne "") {
    	    print $pbs "
if [[ -s $result_file_addqual_bam ]]; then
  rm $rmlist 
fi
";
    }

#preseq
    print $pbs "
if [ ! -s ${sample_name}.c_curve ]; then
  echo preseq c_curve=`date`
  preseq c_curve -B $result_file_addqual_bam > ${sample_name}.c_curve
fi

if [ ! -s ${sample_name}.lc_extrap ]; then
  echo preseq lc_extrap=`date`
  preseq lc_extrap -B -D $result_file_addqual_bam > ${sample_name}.lc_extrap
fi

";

    print $pbs "
if [[ ( -s $result_file_addqual_bam ) && ( ! -s $hs_metrics ) ]]; then
  echo picard CollectHsMetrics =`date`
  java -jar $picard CollectHsMetrics \\
    I=$result_file_addqual_bam \\
    O=$hs_metrics \\
    R=$chrDir \\
    #PER_TARGET_COVERAGE=$per_target_coverage \\
    BAIT_INTERVALS=$interval_list \\
    TARGET_INTERVALS=$interval_list
fi
";

    print $pbs "
if [[ ( -s $result_file_addqual_bam ) && ( ! -s $rrbs_metrics ) ]]; then
  echo picard CollectRrbsMetrics =`date`
  java -jar $picard CollectRrbsMetrics \\
    I=$result_file_addqual_bam \\
    M=$sample_name \\
    R=$chrDir
fi
";

    print $pbs "
if [[ ( -s $hs_metrics ) && ( -s $rrbs_metrics ) ]]; then
  rm $result_file_addqual_bam $result_file_addqual_bai
fi
    ";

    print $pbs "
if [[ ( -s $hs_metrics ) && ( -s $rrbs_metrics ) ]]; then
  echo samtools final result =`date`
  samtools view -h -O SAM -o $result_file $result_file_bam
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
    my $sam_file     = "${result_dir}/${sample_name}.sorted.sam";
    my @result_files = ();
    push( @result_files, $sam_file );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
