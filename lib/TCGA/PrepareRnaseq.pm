#!/usr/bin/perl
package TCGA::PrepareRnaseq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::Task;
use Pipeline::PipelineUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_pr";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $cancerNames = get_option( $config, $section, "source" );
  my $geneLengthFile = get_param_file($config->{$section}{"gene_length_file"},"gene_length_file",1);

  my $script = dirname(__FILE__) . "/PrepareRnaseq.r";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  my $expectFiles = $self->result( $config, $section );
  for my $cancerName (@$cancerNames) {
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $cancerName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $cancerName );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $cur_dir = create_directory_or_die( $result_dir . "/" . $cancerName );

    my $final_file = $expectFiles->{$cancerName}->[0];

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    print $pbs "R --vanilla -f $script --args $cancerName $cur_dir \"$geneLengthFile\" 
rm -rf GDCdata
";
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section, 0 );

  my $cancerNames = get_option( $config, $section, "source" );

  my $result = {};

  for my $cancerName (@$cancerNames) {
    my @result_files = ();
    my $curdir       = $result_dir . "/" . $cancerName;
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.fpkm.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.fpkm.Solid_Tissue_Normal.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.fpkm.Primary_Solid_Tumor.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.fpkm.Solid_Tissue_Normal_paired.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.fpkm.Primary_Solid_Tumor_paired.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.tpm.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.tpm.Solid_Tissue_Normal.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.tpm.Primary_Solid_Tumor.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.tpm.Solid_Tissue_Normal_paired.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.tpm.Primary_Solid_Tumor_paired.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.count.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.count.Solid_Tissue_Normal.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.count.Primary_Solid_Tumor.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.count.Solid_Tissue_Normal_paired.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.count.Primary_Solid_Tumor_paired.tsv" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.gene.pos" );
    push( @result_files, "${curdir}/${cancerName}.rnaseq2.count.rdata" );
    $result->{$cancerName} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $cancerNames = get_option( $config, $section, "source" );

  my $result = {};
  for my $cancerName (@$cancerNames) {
    $result->{$cancerName} = $self->get_pbs_filename( $pbs_dir, $cancerName );
  }

  return $result;
}

1;
