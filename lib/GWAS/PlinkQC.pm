#!/usr/bin/perl
package GWAS::PlinkQC;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_pq";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $rawFiles = get_raw_files( $config, $section );
  
  my $thres_geno = get_option( $config, $section, "thres_geno", 0.01 );
  my $thres_maf = get_option( $config, $section, "thres_maf", 0.01 );
  my $thres_hwe = get_option( $config, $section, "thres_hwe", 0.0001 );
  my $thres_mind = get_option( $config, $section, "thres_mind", 0.03 );

  for my $sampleName ( sort keys %$rawFiles ) {
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sampleName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sampleName );
    my $log_desc = $cluster->get_log_description($log);

    my $files = $rawFiles->{$sampleName};
    my $file = $files->[0];
    my $filePrefix = $file;
    $filePrefix =~ s{\.[^.]+$}{};

    my $qc1_maf = $sampleName . "_1_maf";
    my $qc2_geno = $sampleName . "_2_geno";
    my $qc3_hwe = $sampleName . "_3_hwe";
    my $qc4_mind = $sampleName . "_4_mind";
    my $qc5_resetid = $sampleName . "_5_resetid";
    my $output_file = $sampleName . "_clean";
    
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, "${output_file}.bed" );
    print $pbs "
plink2 --bfile $filePrefix --maf $thres_maf --make-bed --out $qc1_maf
plink2 --bfile $qc1_maf --geno $thres_geno --make-bed --out $qc2_geno
plink2 --bfile $qc2_geno --hwe $thres_maf --make-bed --out $qc3_hwe
plink2 --bfile $qc3_hwe --mind $thres_maf --make-bed --out $qc4_mind
plink2 --bfile $qc4_mind --set-all-var-ids \@:#\\\$r_\\\$a --make-bed --out $qc5_resetid 
plink2 --bfile $qc5_resetid --rm-dup force-first --make-bed --out $output_file

wc -l *.bim *.fam > ${sampleName}_qc_count.txt
rm ${qc1_maf}* ${qc2_geno}*  ${qc3_hwe}*  ${qc4_mind}*  ${qc5_resetid}* 

";
    $self->close_pbs( $pbs, $pbs_file );
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $rawFiles = get_raw_files( $config, $section );

  my $result = {};
  for my $sampleName ( sort keys %$rawFiles ) {
    my @result_files = ("$result_dir/${sampleName}_clean.bed");
    $result->{$sampleName} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
