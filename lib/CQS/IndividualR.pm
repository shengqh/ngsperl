#!/usr/bin/perl
package CQS::IndividualR;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::ProgramWrapperOneToOne;
use CQS::StringUtils;
use File::Spec;
use File::Copy;

our @ISA = qw(CQS::ProgramWrapperOneToOne);

my $directory;

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "";
  $self->{_output_to_same_folder} = 0;
  bless $self, $class;
  return $self;
}

sub getRcode {
  my ( $self, $config, $section ) = @_;
  return get_option( $config, $section, "rCode", "" );
}

sub prepare_data {
  my ( $self, $config, $section ) = @_;
  return;
}

sub get_init_pbs {
  my ( $self, $config, $section ) = @_;
  return "";
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );
  my $task_suffix = $self->{_task_suffix};

  $self->prepare_data($config, $section);
  
  my $init_command = get_option($config, $section, "init_command", "");
  my $rCode = $self->getRcode( $config, $section );
  my $output_file     = get_option( $config, $section, "output_file",     "" );
  my $output_to_result_directory = get_option( $config, $section, "output_to_result_directory", 0 );

  my $run_rmd_independent = get_option( $config, $section, "run_rmd_independent", 0 );
  my $rmd_ext = get_option( $config, $section, "rmd_ext", ".html" );
  my $out_report_at_root_folder = get_option( $config, $section, "out_report_at_root_folder", 1 );

  my $use_vanilla = get_option( $config, $section, "use_vanilla", 1 );

  my $copy_template = get_option( $config, $section, "copy_template", 1 );

  my $removeEmpty = get_option( $config, $section, "remove_empty_parameter", 0 );

  my $raw_files = get_raw_files($config, $section, "parameterSampleFile1");

  my @sample_names = ( sort keys %$raw_files );
  my $has_multi_samples = scalar(@sample_names) > 1;

  my $shfile;
  my $sh;
  if ($has_multi_samples) {
    $shfile = $self->get_task_filename( $pbs_dir, $task_name );
    open( $sh, ">$shfile" ) or die "Cannot create $shfile";
    print $sh get_run_command($sh_direct) . "\n";
  }

  for my $sample_name (@sample_names){
    my $cur_dir = create_directory_or_die($result_dir . "/" . $sample_name);

    my $paramSampleFiles = {};
    for my $myidx (1..10) {
      my $cur_sample_name = $myidx == 1 ? $sample_name : undef;
      my $paramSampleFile = writeParameterSampleFile( $config, $section, $cur_dir, $myidx, $removeEmpty, $cur_sample_name );
      if (($paramSampleFile ne "") or ($myidx <= 3)){
        $paramSampleFiles->{"parSampleFile" . $myidx} = $paramSampleFile;
      }
    }

    my $rfile = $cur_dir . "/${sample_name}${task_suffix}.r";
    open( my $rf, ">$rfile" ) or die "Cannot create $rfile";
    print $rf "rm(list=ls()) \n";
    print $rf "sample_name='$sample_name'\n";

    my $all_results =  $self->result( $config, $section );
    my $cur_results = $all_results->{$sample_name};
    #print(Dumper($all_results));
    my $final_file = $cur_results->[0];

    my $output_file_r = $sample_name . $output_file;

    my $rParameter = "outFile='$output_file_r'\n";
    for my $file_key (sort keys %$paramSampleFiles) {
      my $cur_file = $paramSampleFiles->{$file_key};
      $rParameter = $rParameter . $file_key . "='$cur_file'\n";
    }

    for my $index (1..10){
      my $key = "parameterFile" . $index;
      my $parameterFile = parse_param_file( $config, $section, $key, 0 );
      if ( !defined($parameterFile) && $index < 4 ) {
        $parameterFile = "";
      }
      if ( defined($parameterFile) ) {
        $rParameter = $rParameter . "parFile$index='$parameterFile'\n";
      }
    }
    if ($output_to_result_directory) {
      $rParameter = $rParameter . "outputDirectory='.'\n";
    }

    if ( $rParameter ne "" ) {
      print $rf $rParameter;
    }
    if ( defined($rCode) ) {
      print $rf $rCode . "\n";
    }
    
    print $rf "\nsetwd('$cur_dir')\n";
    
    my $setting_line = "### Parameter setting end ###";
    print $rf "\n$setting_line\n\n";

    my $rtemplates = get_option( $config, $section, "rtemplate" );
    my @rtemplates = split( ",|;", $rtemplates );

    my $rnum = scalar(@rtemplates);

    my $output_sessioninfo = 1;
    my $ignore_lines={};
    for my $index (0 .. $#rtemplates){
      my $rtemplate = $rtemplates[$index];
      my $is_absolute = File::Spec->file_name_is_absolute($rtemplate);
      if ( !$is_absolute ) {
        $rtemplate = dirname(__FILE__) . "/$rtemplate";
      }
      if ( !( -e $rtemplate ) ) {
        die("rtemplate $rtemplate defined but not exists!");
      }

      if($copy_template) {
        if ($index < $rnum-1){
          my $remote_r = $cur_dir . "/" . basename($rtemplate);
          copy($rtemplate, $remote_r);
          if($rtemplate =~ '.R$' || $rtemplate =~ '.r$' ){
            my $line = 'source("' . basename($rtemplate) . '")';
            $ignore_lines->{$line} = 1;
            $line =~ s/"/'/g;
          $ignore_lines->{$line} = 0;
          }
          next;
        }else{
          for my $line (keys %$ignore_lines){
            if($ignore_lines->{$line}){
              print $rf "$line\n";
            }
          }
        }
      }
      
      my @valid_lines = ();
      
      open( my $rt, "<$rtemplate" ) or die $!;
      while ( my $row = <$rt> ) {
        chomp($row);
        $row =~ s/\r//g;
        if($row eq $setting_line){
          @valid_lines = ();
          next;
        }

        if($copy_template && ($row =~/^source/)){
          if(defined $ignore_lines->{$row}){
            next;
          }
        }

        push(@valid_lines, $row);
      }
      close($rt);

      my $first = 1;
      for my $row (@valid_lines){
        if ($row eq ""){
          if ($first){
            next
          }
        }else{
          $first = 0;
        }
        print $rf "$row\n";

        if ($row =~ /sessionInfo/){
          $output_sessioninfo = 0;
        }
      }
    }

    if ($output_sessioninfo){
      print $rf "writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')\n";
    }

    close($rf);

    my $rmd_file = undef;
    my $rReportTemplates = get_option( $config, $section, "rReportTemplate", "" );
    if ( $rReportTemplates ne "" ) {
      my @rReportTemplates = split( ",|;", $rReportTemplates );
      for my $rtemplate (@rReportTemplates){
        my $is_absolute = File::Spec->file_name_is_absolute($rtemplate);
        if ( !$is_absolute ) {
          $rtemplate = dirname(__FILE__) . "/$rtemplate";
        }
        if ( !( -e $rtemplate ) ) {
          die("rmd template $rtemplate defined but not exists!");
        }

        my $remote_r = $cur_dir . "/" . basename($rtemplate);
        copy($rtemplate, $remote_r);

        if(!defined $rmd_file){
          $rmd_file = basename($rtemplate);
        }
      }
    }

    my $additional_rmd_files = get_option( $config, $section, "additional_rmd_files","" );
    if ( $additional_rmd_files ne "" ) {
      my @additional_rmd_files = split( ",|;", $additional_rmd_files );
      foreach my $additional_rmd_files (@additional_rmd_files) {
        my $is_absolute = File::Spec->file_name_is_absolute($additional_rmd_files);
        if ( !$is_absolute ) {
          $additional_rmd_files = dirname(__FILE__) . "/$additional_rmd_files";
        }
        if ( !( -e $additional_rmd_files ) ) {
          die("additional_rmd_files $additional_rmd_files defined but not exists!");
        }
        my $additional_rmd_files_inResult = $cur_dir . "/" . basename($additional_rmd_files);
        open( my $rf, ">$additional_rmd_files_inResult" ) or die "Cannot create $additional_rmd_files_inResult";
        open( my $rt, "<$additional_rmd_files" ) or die $!;
        while ( my $row = <$rt> ) {
          chomp($row);
          $row =~ s/\r//g;
          print $rf "$row\n";
        }
        close($rt);
        close($rf);
      }
    }

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    if ( $has_multi_samples ) {
      print $sh "if [[ ! -s $final_file ]]; then
  \$MYCMD ./$pbs_name 
fi
";
    }

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file, $init_command, $self->can_result_be_empty_file($config, $section) );

    my $rlibs = get_option_include_general($config, $section, "R_LIBS", "");
    if($rlibs ne ""){
      print $pbs "export R_LIBS=$rlibs \n\n";
    }

    print $pbs "$init_command\n";
    print $pbs $self->get_init_pbs($config, $section) . "\n";
    
    my $rscript = get_option_include_general($config, $section, "Rscript", "Rscript");

    my $vanilla_option = $use_vanilla ? "--vanilla ":"";
    if ( defined($option) and $option ne "" ) {
      print $pbs "$rscript $vanilla_option " . basename($rfile) . " $option";
    }
    else {
      print $pbs "$rscript $vanilla_option " . basename($rfile);
    }

    if($run_rmd_independent){
      if(defined $rmd_file){
        my $folder = $out_report_at_root_folder ? "../" : "";
        my $rmd_command = "$rscript $vanilla_option -e \"library('rmarkdown');rmarkdown::render('$rmd_file',output_file='$folder${sample_name}${rmd_ext}')\"";
        print $pbs "

  if [[ -s $final_file ]]; then
    $rmd_command
  fi
  ";
        my $report_sh = "$pbs_dir/${sample_name}_report.sh";
        open( my $rs, ">$report_sh" ) or die $!;
        print $rs "
  cd $cur_dir

  $rmd_command
  ";
        close($rs);
      }
    }
    $self->close_pbs( $pbs, $pbs_file );
  }

  if ( $has_multi_samples ) {
    close $sh;
    if ( is_linux() ) {
      chmod 0755, $shfile;
    }

    print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
  }
}

1;
