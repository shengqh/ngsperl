#!/usr/bin/perl
package CQS::IndividualRmd;

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
  my $post_command = get_option($config, $section, "post_command", "");
  
  my $output_ext = get_option( $config, $section, "output_ext", ".html" );

  my $use_vanilla = get_option( $config, $section, "use_vanilla", 1 );

  my $removeEmpty = get_option( $config, $section, "remove_empty_parameter", 0 );

  my $source_key = $self->get_pbs_key($config, $section);
  my $raw_files = get_raw_files($config, $section, $source_key);

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

    my $final_file = $self->get_absolute_final_file($config, $section, $sample_name);

    my $rmd_file = undef;
    my $rReportTemplates = get_option( $config, $section, "rReportTemplate");
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

    print $pbs $self->get_init_pbs($config, $section) . "\n";
    
    my $rscript = get_option_include_general($config, $section, "Rscript", "Rscript");

    my $vanilla_option = $use_vanilla ? "--vanilla " : "";

    print $pbs "

$rscript $vanilla_option -e \"library('rmarkdown');rmarkdown::render('$rmd_file',output_file='${sample_name}${output_ext}')\"

";

    print $pbs "\n\n$post_command\n\n";
    print $pbs "\n\nrm -rf .local .cache .java\n";

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
