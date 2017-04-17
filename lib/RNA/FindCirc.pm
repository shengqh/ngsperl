#!/usr/bin/perl
package RNA::FindCirc;

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
  $self->{_suffix} = "_fc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $FindCircDir = get_option( $config, $section, "FindCircDir", dirname(__FILE__) );
  my $bowtie2_index = $config->{$section}{bowtie2_index} or die "define ${section}::bowtie2_index first";
  my $genomeFa = $config->{$section}{genomeFile} or die "define ${section}::genomeFile first";

  my $python_script1 = $FindCircDir . "/unmapped2anchors.py";
  if ( !-e $python_script1 ) {
    die "File not found : " . $python_script1;
  }
  my $python_script2 = $FindCircDir . "/find_circ.py";
  if ( !-e $python_script2 ) {
    die "File not found : " . $python_script2;
  }

  my %raw_files     = %{ get_raw_files( $config, $section ) };

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile   = $sample_files[0];
    my $final_file   = "${sample_name}_splice_sites.bed";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
    print $pbs "
cd  $result_dir

echo find_circ step1
bowtie2 -p16 --very-sensitive --score-min=C,-15,0 --mm -x $bowtie2_index -q -U $sampleFile >${sample_name}.sam 2> ${sample_name}.bowtie2.log 
samtools view -hbuS ${sample_name}.sam >${sample_name}.bam
samtools sort ${sample_name}.bam >${sample_name}.sorted.bam


echo find_circ step2
# get the unmapped and pipe through unmapped2anchors.py
samtools view -hf 4 ${sample_name}.sorted.bam | samtools view -Sb >${sample_name}.sorted.unmapped.bam
$python_script1 ${sample_name}.sorted.unmapped.bam | gzip > ${sample_name}_anchors.fastq.gz

      
echo find_circ step3
bowtie2 -p 16 --score-min=C,-15,0 --reorder --mm \
        -q -U ${sample_name}_anchors.fastq.gz -x $bowtie2_index |\
            $python_script2 \
                --genome=$genomeFa \
                --prefix=${sample_name}_ \
                --name=${sample_name}_sample \
                --stats=${sample_name}_stats.txt \
                --reads=${sample_name}_spliced_reads.fa \
                    > $final_file
";
    $self->close_pbs( $pbs, $pbs_file );

  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $splice_file     = "${result_dir}/${sample_name}_splice_sites.bed";
    my @result_files = ();
    push( @result_files, $splice_file );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
