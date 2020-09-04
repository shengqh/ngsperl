#!/usr/bin/perl
package Chipseq::Encode;

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
use Data::Dumper;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ec";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $chip_wdl = get_option_file($config, $section, "chip.wdl");
  my $chip_genome_tsv = get_option_file($config, $section, "chip.genome_tsv");
  my $singularity_file = get_option_file($config, $section, "singularity.file");
  #my $ispairend = get_is_paired_end_option( $config, $section );

  my $raw_files = get_grouped_raw_files($config, $section, "treatments");
  my $control_files = get_grouped_raw_files( $config, $section, "inputs" );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %$raw_files ) {
    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    my $jsonfile = "$cur_dir/${sample_name}.json";
    open(my $fjson, ">:utf8", $jsonfile) or die $!;
    print $fjson "{
    \"chip.pipeline_type\" : \"histone\",
    \"chip.genome_tsv\" : \"$chip_genome_tsv\",
";

    my $sample_files = $raw_files->{$sample_name};
    if (scalar(@$sample_files) == 2) {
      my $sample_file_1 = $sample_files->[0];
      my $sample_file_2 = $sample_files->[1];
      print $fjson "    \"chip.paired_end\" : true,
    \"chip.fastqs_rep1_R1\" : [\"$sample_file_1\" ],
    \"chip.fastqs_rep1_R2\" : [\"$sample_file_2\" ],
";
    }else{
      my $sample_file_1 = $sample_files->[0];
      print $fjson "    \"chip.paired_end\" : false,
    \"chip.fastqs_rep1_R1\" : [\"$sample_file_1\" ],
";
    }

    my $control_files = $control_files->{$sample_name};
    if (scalar(@$control_files) == 2) {
      my $control_file_1 = $control_files->[0];
      my $control_file_2 = $control_files->[1];
      print $fjson "    \"chip.ctl_paired_end\" : true,
    \"chip.ctl_fastqs_rep1_R1\" : [\"$control_file_1\" ],
    \"chip.ctl_fastqs_rep1_R2\" : [\"$control_file_2\" ],
";
    }else{
      my $control_file_1 = $control_files->[0];
      print $fjson "    \"chip.ctl_paired_end\" : false,
    \"chip.ctl_fastqs_rep1_R1\" : [\"$control_file_1\" ],
";
    }

    print $fjson "    \"chip.always_use_pooled_ctl\" : false,
    \"chip.title\" : \"$sample_name\",
    \"chip.description\" : \"$sample_name\"
}
";

    close($fjson);

    my $final_file = "${sample_name}_peaks.txt";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    print $pbs "
caper run $chip_wdl $option -i $jsonfile --singularity $singularity_file

";

    $self->close_pbs( $pbs, $pbs_file );

    print $sh "\$MYCMD ./$pbs_name \n";
  }

  print $sh "exit 0\n";
  close $sh;
  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $raw_files = get_grouped_raw_files($config, $section, "treatments");

  my $result = {};
  for my $sample_name ( sort keys %$raw_files ) {
    my $cur_dir      = $result_dir . "/$sample_name/chip";
    my @result_files = ();
    push( @result_files, $cur_dir );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
