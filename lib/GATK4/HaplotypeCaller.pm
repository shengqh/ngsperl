#!/usr/bin/perl
package GATK4::HaplotypeCaller;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use GATK4::GATK4Task;

our @ISA = qw(GATK4::GATK4Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_hc";
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

  my $by_chromosome = get_option( $config, $section, "by_chromosome", 0 );
  my $gvcf =get_option( $config, $section, "gvcf", 1 );
  if($gvcf){
    $option = $option . " -ERC GVCF";
  }

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $blacklist_intervals = get_param_file( $config->{$section}{blacklist_file}, "blacklist_file", 0 );
  my $blacklist_intervals_option = $blacklist_intervals ? "-XL " . $blacklist_intervals : "";

  my $extension = get_option( $config, $section, "extension", ".g.vcf.gz" );

  my $java_option = $self->get_java_option($config, $section, $memory);

  $self->get_docker_value(1);

  my $bedFile            = get_param_file( $config->{$section}{bed_file}, "bed_file", 0 );
  my $interval_padding   = get_option( $config, $section, "interval_padding", 0 );
  my $restrict_intervals="";
  if (defined $bedFile and $bedFile ne "") {
  	if (defined $interval_padding and $interval_padding!=0) {
  		$restrict_intervals="-L $bedFile -ip $interval_padding";
  	} else {
  		$restrict_intervals="-L $bedFile";
  	}
  }
  
  my %bam_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %bam_files ) {
    my @sample_files = @{ $bam_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    my $snvOut    = $sample_name . $extension;

    #if the program throw exception, the idx file will not be generated.
    my $snvOutIndex;
    if ($extension =~ ".gz\$") {
      $snvOutIndex = $snvOut . ".tbi";
    } else {
      $snvOutIndex = $snvOut . ".idx";
    }

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );


    print $sh "if [[ ! -s $result_dir/$snvOutIndex ]]; then
  \$MYCMD ./$pbs_name 
fi
    
";

    my $log_desc = $cluster->get_log_description($log);
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $snvOutIndex );

    my $localized_files = [];
    @sample_files = @{$self->localize_files_in_tmp_folder($pbs, \@sample_files, $localized_files, [".bai"])};
    $bam_file    = $sample_files[0];

    print $pbs "
gatk --java-options \"$java_option\" \\
  HaplotypeCaller $option $restrict_intervals $blacklist_intervals_option \\
  --native-pair-hmm-threads $thread \\
  -R $faFile \\
  -I $bam_file \\
  -O $snvOut

status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.hc.failed
  rm -f $snvOut $snvOutIndex
else
  touch $sample_name.hc.succeed
fi

";

    $self->clean_temp_files($pbs, $localized_files);

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all GATK4::HaplotypeCaller tasks. \n ";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  my $result = {};

  my $extension = get_option( $config, $section, "extension", ".g.vcf.gz" );

  my %bam_files = %{ get_raw_files( $config, $section ) };
  for my $sample_name ( sort keys %bam_files ) {
    my $snvOut       = $sample_name . $extension;
    my @result_files = ();
    push( @result_files, "${result_dir}/${snvOut}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
