#!/usr/bin/perl
package CQS::BuildReport;

use strict;
use warnings;
use File::Basename;
use File::Copy;
use Data::Dumper;
use String::Util ':all';
use CQS::UniqueTask;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_br";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $final_pbs      = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $final_log      = $self->get_log_filename( $log_dir, $task_name );
  my $final_log_desp = $cluster->get_log_description($final_log);

  my $rtemplate = dirname(__FILE__) . "/" . $config->{$section}{report_rmd_file};
  my $rfile     = $result_dir . "/${task_name}.Rmd";
  copy( $rtemplate, $rfile ) or die "Copy failed: $!";

  if ( defined $config->{$section}{additional_rmd_files} ) {
    my @additional_rtemplates = split( ';', $config->{$section}{additional_rmd_files} );
    for my $additional (@additional_rtemplates) {
      my $additional_rtempalte = dirname(__FILE__) . "/" . trim($additional);
      my $additional_rfile     = $result_dir . "/" . basename($additional);
      copy( $additional_rtempalte, $additional_rfile ) or die "Copy failed: $!";
    }
  }

  my $raw_file_list = get_raw_file_list($config, $section, "parameterSampleFile1");
  my $raw_file_names = $config->{$section}{parameterSampleFile1_names};
  
  die "File lists (". scalar(@$raw_file_list) . ") is not equals to file names (" . scalar(@$raw_file_names) . ")" if scalar(@$raw_file_list) != scalar(@$raw_file_names);
  
  open( my $list, ">$result_dir/fileList1.txt" ) or die "Cannot create $result_dir/fileList1.txt";
  while (my ($index, $element) = each(@$raw_file_list)) {
    print $list "$element\t" . $raw_file_names->[$index] . "\n";
  }
  close($list);
  
  writeParameterSampleFile($config, $section, $result_dir, 2);

  my $copy_file_list = get_raw_file_list($config, $section, "parameterSampleFile3");
  my $report_folder = create_directory_or_die("$result_dir/$task_name");
  my $final_file     = "$task_name/${task_name}.html";
  my $final = $self->open_pbs( $final_pbs, $pbs_desc, $final_log_desp, $path_file, $result_dir );
  
  for my $copy_file (@$copy_file_list){
    my $to_folder = $task_name;
    if ($copy_file =~ /_WebGestalt/){
      create_directory_or_die("$report_folder/Functional_enrichment/");
      create_directory_or_die("$report_folder/Functional_enrichment/webGestalt");
      $to_folder = "$task_name/Functional_enrichment/webGestalt";
    }elsif($copy_file =~ /_GSEA/){
      create_directory_or_die("$report_folder/Functional_enrichment/");
      create_directory_or_die("$report_folder/Functional_enrichment/gsea");
      $to_folder = ("$task_name/Functional_enrichment/gsea");
    }
    print $final "cp -r -u $copy_file $to_folder \n";
  }
  
  print $final "
R --slave -e \"library(knitr);rmarkdown::render('" . "${task_name}.Rmd');\"
";

  $self->close_pbs( $final, $final_pbs );
}

1;
