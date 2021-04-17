#!/usr/bin/perl
package ParClip::Bowtie1PARalyzer;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_bp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  #bowtie option
  my $bowtie1_index = $config->{$section}{bowtie1_index} or die "define ${section}::bowtie1_index first";
  my $mappedonlyoption = "-F 4";

  #paralyzer option
  my $genome2bit = get_param_file( $config->{$section}{genome2bit}, "genome2bit", 1 );
  my $mirna_db   = get_param_file( $config->{$section}{mirna_db},   "mirna_db",   1 );
  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $sam_file = $sample_name . ".sam";
    my $bam_file = $sample_name . ".bam";
    my $alignlog = $sample_name . ".log";

    #prepare paralyzer configuration file
    my $iniFile = "${sample_name}.ini";
    open( INI, ">${result_dir}/${iniFile}" ) or die "Cannot create ${result_dir}/${iniFile}";
    print INI "
BANDWIDTH=3
CONVERSION=T>C
MINIMUM_READ_COUNT_PER_GROUP=5
MINIMUM_READ_COUNT_PER_CLUSTER=5
MINIMUM_READ_COUNT_FOR_KDE=5
MINIMUM_CLUSTER_SIZE=10
MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER=2
MINIMUM_CONVERSION_COUNT_FOR_CLUSTER=1
MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION=5
MINIMUM_READ_LENGTH=16
MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES=5

ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL=5

SAM_FILE=${result_dir}/$sam_file
GENOME_2BIT_FILE=$genome2bit
FIND_MIRNA_SEEDMATCHES=$mirna_db

OUTPUT_DISTRIBUTIONS_FILE=${sample_name}.distribution.csv
OUTPUT_GROUPS_FILE=${sample_name}.groups.csv
OUTPUT_CLUSTERS_FILE=${sample_name}.cluster.csv
OUTPUT_MIRNA_TARGETS_FILE=${sample_name}.target.csv
    
";
    close(INI);

    my $m_option = ( $option =~ /\-m/ ) ? "--max ${bam_file}.max.txt" : "";

    my $indent = "";
    my $tag    = "--sam-RG ID:$sample_name --sam-RG LB:$sample_name --sam-RG SM:$sample_name --sam-RG PL:ILLUMINA --sam-RG PU:$sample_name";

    my $cmd_file_exists = check_file_exists_command( \@sample_files, "" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $final_file = "${sample_name}.target.csv";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    my $collapserFasta = $sample_name . ".collapser.fasta";
    my $fastqs         = join( ' ', @sample_files );
    my $cat            = ( $sample_files[0] =~ /.gz/ ) ? "zcat" : "cat";
    print $pbs "
$cmd_file_exists
$cat $fastqs | fastx_collapser -o $collapserFasta 

if [[ -s $collapserFasta ]]; then
  bowtie $option $m_option -S $tag -f $bowtie1_index $collapserFasta | samtools view -hSF4 > $sam_file
fi

if [[ -s $sam_file ]]; then
  PARalyzer $memory $iniFile
fi

if [[ -s $final_file && -s $sam_file ]]; then
  samtools sort -T ${sample_name}_tmp -o $bam_file $sam_file
  if [ -s $bam_file ]; then
    rm $collapserFasta $sam_file 
    samtools index $bam_file 
  fi
fi
";

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all ParClip::PARalyzer tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    push( @result_files, "${result_dir}/${sample_name}/${sample_name}.bam" );
    push( @result_files, "${result_dir}/${sample_name}/${sample_name}.bam.log" );
    push( @result_files, "${result_dir}/${sample_name}/${sample_name}.cluster.csv" );
    push( @result_files, "${result_dir}/${sample_name}/${sample_name}.distribution.csv" );
    push( @result_files, "${result_dir}/${sample_name}/${sample_name}.groups.csv" );
    push( @result_files, "${result_dir}/${sample_name}/${sample_name}.target.csv" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
