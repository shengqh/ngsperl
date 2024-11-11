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
  $self->{_docker_prefix} = "report_";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $final_pbs = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $final_log = $self->get_log_filename( $log_dir, $task_name );
  my $final_log_desp = $cluster->get_log_description($final_log);

  my $rtemplate = dirname(__FILE__) . "/" . $config->{$section}{report_rmd_file};
  my $rfile     = $result_dir . "/${task_name}.Rmd";

  if ( defined $config->{$section}{function_r_files} ) {
    my @function_r_files = split( ';', $config->{$section}{function_r_files} );

    open( my $rmd, ">$rfile" )     or die "Cannot create $rfile";
    open( my $rin, "<$rtemplate" ) or die "Cannot open $rtemplate";
    my $function_done = 0;
    while (<$rin>) {
      my $rinline = $_;
      $rinline = s/\r|\n//g;
      if ( $rinline =~ /^```\{/ ) {
        if ( !$function_done ) {
          print $rmd '```{r, include=FALSE} \n';
          for my $func_file (@function_r_files) {
            $func_file = dirname(__FILE__) . "/" . trim($func_file);
            open( my $fin, "<$func_file" ) or die "Cannot open $func_file";
            while (<$fin>) {
              my $finline = $_;
              $finline = s/\r|\n//g;
              print $rmd $finline . "\n";
            }
          }
          print $rmd '``` \n';
          $function_done = 1;
        }
      }
    }
  }
  else {
    copy( $rtemplate, $rfile ) or die "Copy failed: $! : " . $rtemplate;
  }

  my $additional_rmd_files = $config->{$section}{additional_rmd_files};
  if ( defined $additional_rmd_files ) {
    my @additional_rtemplates = split( ';', $additional_rmd_files );
    for my $additional (@additional_rtemplates) {
      my $additional_rtempalte = dirname(__FILE__) . "/" . trim($additional);
      my $additional_rfile     = $result_dir . "/" . basename($additional);
      copy( $additional_rtempalte, $additional_rfile ) or die "Copy failed for $additional_rtempalte: $!";
    }
  }
  my $raw_file_list = get_raw_file_list( $config, $section, "parameterSampleFile1", 1 );
  my $raw_file_names = $config->{$section}{parameterSampleFile1_names};

  if(defined($raw_file_names)){
    if ( scalar(@$raw_file_list) != scalar(@$raw_file_names) ) {
      print( "\nRaw file list = \n" . join( "\n", @$raw_file_list ) );
      print( "\nRaw file names = \n" . join( "\n", @$raw_file_names ) );
      die "File lists (" . scalar(@$raw_file_list) . ") is not equals to file names (" . scalar(@$raw_file_names) . ")";
    }

    open( my $list, ">$result_dir/fileList1.txt" ) or die "Cannot create $result_dir/fileList1.txt";
    while ( my ( $index, $element ) = each(@$raw_file_list) ) {
      print $list "$element\t" . $raw_file_names->[$index] . "\n";
    }
    close($list);
  }else{
    writeParameterSampleFile( $config, $section, $result_dir, 1 );
  }

  writeParameterSampleFile( $config, $section, $result_dir, 2 );
  writeParameterSampleFile( $config, $section, $result_dir, 4, 1 );
  writeParameterSampleFile( $config, $section, $result_dir, 5, 1 );
  writeParameterSampleFile( $config, $section, $result_dir, 6, 1 );
  writeParameterSampleFile( $config, $section, $result_dir, 7, 1 );
  writeParameterSampleFile( $config, $section, $result_dir, 8, 1 );
  writeParameterSampleFile( $config, $section, $result_dir, 9, 1 );

  my $final_file = "${task_name}.html";
  my $final = $self->open_pbs( $final_pbs, $pbs_desc, $final_log_desp, $path_file, $result_dir );

  if(has_raw_files($config, $section, "parameterSampleFile3")){
    my $copy_file_list = get_raw_file_list( $config, $section, "parameterSampleFile3" );
    if ( scalar(@$copy_file_list) > 0 ) {
      my $report_folder = create_directory_or_die("$result_dir/$task_name");

      for my $copy_file ( sort @$copy_file_list ) {
        my $to_folder = $task_name;
        #print($copy_file . "\n");
        if ( $copy_file =~ /.txt.html/ ) {
          $to_folder = $task_name;
          #print("  " . $to_folder . "\n");
        } elsif ( $copy_file =~ /_WebGestalt/ ) {
          create_directory_or_die("$report_folder/webGestalt");
          $to_folder = "$task_name/webGestalt";
        } elsif ( $copy_file =~ /_GSEA/ & $copy_file !~ /_GSEA.rnk$/ ) {
          create_directory_or_die("$report_folder/gsea");
          $to_folder = ("$task_name/gsea");
        } elsif ( $copy_file =~ /_homer/ ) {
          my $treatment=basename(dirname($copy_file));
          create_directory_or_die("$report_folder/homer");
          create_directory_or_die("$report_folder/homer/$treatment");
          $to_folder = ("$task_name/homer/$treatment");
        }
        print $final "cp -r -u $copy_file $to_folder \n";
      }
    }
  }
  print $final "
R -e \"library(knitr);rmarkdown::render('${task_name}.Rmd');\"

rm -rf .cache

";

  $self->close_pbs( $final, $final_pbs );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result       = {};
  my @result_files = ();
  push( @result_files, "${result_dir}/${task_name}.html" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
