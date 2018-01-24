#!/usr/bin/perl
package Annotation::WebGestaltR;

use strict;
use warnings;
use File::Basename;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::Task;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_wr";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $raw_files = get_raw_files( $config, $section );
  my $organism         = get_option( $config, $section, "organism" );
  my $interestGeneType = get_option( $config, $section, "interestGeneType", "genesymbol" );
  my $referenceSet     = get_option( $config, $section, "referenceSet", "genome" );

  my $script = dirname(__FILE__) . "/WebGestaltR.r";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . " \n";
  
  my $expect_result = $self->result($config, $section);

  for my $sample_name ( sort keys %$raw_files ) {
    my $cur_dir = scalar( keys %$raw_files ) == 1 ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";
    my $final_file = $expect_result->{$sample_name}[0];
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    my $inputFile = $raw_files->{$sample_name}->[0];

    print $pbs " 
R --vanilla -f $script --args $organism $sample_name $inputFile . $interestGeneType $referenceSet
rm */*.tar.gz
";
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );
  my $raw_files = get_raw_files( $config, $section );

  my $result       = {};
  for my $sample_name ( sort keys %$raw_files ) {
    my $cur_dir = scalar( keys %$raw_files ) == 1 ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );
    my $finalFile = "Project_" . $sample_name . "_geneontology_Biological_Process/Report_" . $sample_name . "_geneontology_Biological_Process.html";
    my @result_files = ();
    push( @result_files, "$cur_dir/Project_${sample_name}_geneontology_Biological_Process/enrichment_results_${sample_name}_geneontology_Biological_Process.txt" );
    push( @result_files, "$cur_dir/Project_${sample_name}_geneontology_Cellular_Component/enrichment_results_${sample_name}_geneontology_Cellular_Component.txt" );
    push( @result_files, "$cur_dir/Project_${sample_name}_geneontology_Molecular_Function/enrichment_results_${sample_name}_geneontology_Molecular_Function.txt" );
    push( @result_files, "$cur_dir/Project_${sample_name}_pathway_KEGG/enrichment_results_${sample_name}_pathway_KEGG.txt" );
    push( @result_files, "$cur_dir/Project_${sample_name}_geneontology_Biological_Process" );
    push( @result_files, "$cur_dir/Project_${sample_name}_geneontology_Cellular_Component" );
    push( @result_files, "$cur_dir/Project_${sample_name}_geneontology_Molecular_Function" );
    push( @result_files, "$cur_dir/Project_${sample_name}_pathway_KEGG" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
