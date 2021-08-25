#!/usr/bin/perl
package Annotation::WebGestaltR;

use strict;
use warnings;
use File::Basename;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $comparisons = get_raw_files( $config, $section, undef, undef, 1 );
  my @comparison_names = keys %{$comparisons};

  #writeParameterSampleFile( $config, $section, $result_dir, 1 );

  my $organism         = get_option( $config, $section, "organism" );
  my $interestGeneType = get_option( $config, $section, "interestGeneType", "genesymbol" );
  my $referenceSet     = get_option( $config, $section, "referenceSet", "genome" );

  my $script1 = dirname(__FILE__) . "/WebGestaltReportFunctions.r";
  my $script2 = dirname(__FILE__) . "/WebGestaltR.r";
  if ( (!-e $script1) | (!-e $script2) ) {
     die "File not found : " . $script1." or ".$script2;
  }

  my $script = $result_dir . "/${task_name}.WebGestaltR.r";
  open( my $rf, ">$script" ) or die "Cannot create $script";
  open( my $rt, "<$script1" ) or die $!;
  while ( my $row = <$rt> ) {
      chomp($row);
      $row =~ s/\r//g;
      print $rf "$row\n";
  }
  close($rt);
  open( $rt, "<$script2" ) or die $!;
  while ( my $row = <$rt> ) {
      chomp($row);
      $row =~ s/\r//g;
      print $rf "$row\n";
  }
  close($rt);
  close($rf);

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);
  my $expect_result = $self->result($config, $section, ".txt\$");

  my $pbs                = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );
  for my $sample_name ( sort keys %$comparisons ) {
    my $cur_dir = scalar( keys %$comparisons ) == 1 ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    my $final_file = $expect_result->{$sample_name}[-1];

    my $inputFile = $comparisons->{$sample_name}->[0];

    print $pbs " 
if [[ ! -s $final_file || ! -d $final_file ]]; then
  cd $cur_dir 
  if [[ -f $inputFile ]]; then
    if [[ -s $inputFile ]]; then
      R --vanilla -f $script --args $organism $sample_name $inputFile . $interestGeneType $referenceSet
      rm */*/*.zip
    else
      echo \"Empty gene file\" > ${sample_name}.empty
    fi 
  else
    echo \"Gene file not exist: $inputFile\" .
  fi
fi

";
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 1 );
  my $comparisons = get_raw_files( $config, $section, undef, undef, 1 );

  my $result       = {};
  for my $sample_name ( sort keys %$comparisons ) {
    my $cur_dir = scalar( keys %$comparisons ) == 1 ? $result_dir : $result_dir . "/$sample_name" ;
    my @result_files = ();
    push( @result_files, "$cur_dir/Project_${sample_name}_geneontology_Biological_Process/enrichment_results_${sample_name}_geneontology_Biological_Process.txt" );
    push( @result_files, "$cur_dir/Project_${sample_name}_geneontology_Cellular_Component/enrichment_results_${sample_name}_geneontology_Cellular_Component.txt" );
    push( @result_files, "$cur_dir/Project_${sample_name}_geneontology_Molecular_Function/enrichment_results_${sample_name}_geneontology_Molecular_Function.txt" );
    push( @result_files, "$cur_dir/Project_${sample_name}_pathway_KEGG/enrichment_results_${sample_name}_pathway_KEGG.txt" );
    push( @result_files, "$cur_dir/Project_${sample_name}_geneontology_Biological_Process" );
    push( @result_files, "$cur_dir/Project_${sample_name}_geneontology_Cellular_Component" );
    push( @result_files, "$cur_dir/Project_${sample_name}_geneontology_Molecular_Function" );
    push( @result_files, "$cur_dir/Project_${sample_name}_pathway_KEGG" );
    push( @result_files, "$cur_dir/WebGestaltR.version" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern, 1 );
  }
  return $result;
}

1;
