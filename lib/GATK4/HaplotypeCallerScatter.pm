#!/usr/bin/perl
package GATK4::HaplotypeCallerScatter;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use GATK4::IntervalsScatterTask;

our @ISA = qw(GATK4::IntervalsScatterTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_hcs";
  $self->{_docker_prefix}   = "gatk4_";
  bless $self, $class;

  return $self;
}

#The results from whole file and by chromosome may be different at BaseQRankSum, MQRankSum or ReadPosRankSum
#< 1     1291126 .       G       A,<NON_REF>     6.79    .       BaseQRankSum=0.358;ClippingRankSum=0.358;DP=6;MLEAC=1,0;MLEAF=0.500,0.00;MQ=60.00;MQRankSum=0.358;ReadPosRankSum=-1.231 GT:AD:DP:GQ:PL:SB     0/1:3,2,0:5:34:34,0,59,43,65,109:2,1,1,1
#> 1     1291126 .       G       A,<NON_REF>     6.79    .       BaseQRankSum=1.231;ClippingRankSum=0.358;DP=6;MLEAC=1,0;MLEAF=0.500,0.00;MQ=60.00;MQRankSum=0.358;ReadPosRankSum=-0.358 GT:AD:DP:GQ:PL:SB     0/1:3,2,0:5:34:34,0,59,43,65,109:2,1,1,1
# based on https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS&p=2
sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $scatter_map = get_interval_file_map($config, $section);

  my $gvcf =get_option( $config, $section, "gvcf", 1 );
  if($gvcf){
    $option = $option . " -ERC GVCF";
  }

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $blacklist_intervals = get_param_file( $config->{$section}{blacklist_file}, "blacklist_file", 0 );
  my $blacklist_intervals_option = $blacklist_intervals ? "-XL " . $blacklist_intervals : "";

  my $extension = get_option( $config, $section, "extension", ".g.vcf" );

  my $java_option = $self->get_java_option($config, $section, $memory);

  $self->get_docker_value(0);

  my $interval_padding   = get_option( $config, $section, "interval_padding", 0 );
  my $restrict_intervals="";
  if ($interval_padding!=0) {
    $restrict_intervals="-ip $interval_padding";
  }
  
  my %bam_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %bam_files ) {
    my @sample_files = @{ $bam_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    for my $scatter_name (sort keys %$scatter_map) {
      my $interval_file = $scatter_map->{$scatter_name};
      my $prefix = get_key_name($sample_name, $scatter_name);
      my $snvOut = $prefix . $extension;
      my $snvTmp = $prefix . ".tmp". $extension;
      my $snvTmpIndex;

      #if the program throw exception, the idx file will not be generated.
      my $snvOutIndex;
      if ($extension =~ ".gz\$") {
        $snvOutIndex = $snvOut . ".tbi";
        $snvTmpIndex = $snvTmp . ".tbi";
      } else {
        $snvOutIndex = $snvOut . ".idx";
        $snvTmpIndex = $snvTmp . ".idx";
      }

      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $prefix );
      my $pbs_name = basename($pbs_file);
      my $log      = $self->get_log_filename( $log_dir, $prefix );

      print $sh "if [[ ! -s $cur_dir/$snvOutIndex ]]; then 
  \$MYCMD ./$pbs_name 
fi
";

      my $log_desc = $cluster->get_log_description($log);
      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $snvOutIndex );
      print $pbs "
gatk --java-options \"$java_option\" \\
  HaplotypeCaller $option \\
  -L $interval_file $blacklist_intervals_option \\
  --native-pair-hmm-threads $thread \\
  -R $faFile \\
  -I $bam_file \\
  -O $snvTmp

status=\$?
if [[ \$status -eq 0 ]]; then
  touch ${snvOut}.succeed
  rm -f ${snvOut}.failed
  mv $snvTmp $snvOut
  mv $snvTmpIndex $snvOutIndex
else
  touch ${snvOut}.failed
  rm -f ${snvOut}.succeed
  rm -f $snvTmp
  rm -f $snvTmpIndex
fi
";
      
      $self->close_pbs( $pbs, $pbs_file );
    }
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all GATK4::HaplotypeCallerScatter tasks. \n ";
}

sub get_result_files {
  my ( $self, $config, $section, $result_dir, $sample_name, $scatter_name, $key_name ) = @_;
  my $extension = get_option( $config, $section, "extension", ".g.vcf" );
  my $final_file = "${result_dir}/${sample_name}/${key_name}${extension}";
  return [$final_file];
}

1;
