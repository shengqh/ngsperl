#!/usr/bin/perl
package Alignment::Shrimp2;

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
  $self->{_suffix} = "_srp2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $shrimp2_index = $config->{$section}{shrimp2_index} or die "define ${section}::shrimp2_index first";
  die "shrimp2_index ${shrimp2_index}.genome not exist" if ( !-e "${shrimp2_index}.genome" );

  my $is_mirna = $config->{$section}{is_mirna} or die "define ${section}::is_mirna first";
  my $mirna = "-M mirna" if $is_mirna or "";

  my $output_bam = $config->{$section}{output_bam} or die "define ${section}::output_bam first";

  my %raw_files = %{ get_raw_files( $config, $section, "source", ".fastq\$" ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile   = $sample_files[0];

    my $shrimpFile = $sample_name . ".shrimp";
    my $sam_file   = $sample_name . ".sam";
    my $bam_file   = $sample_name . ".bam";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    my $log_desc = $cluster->get_log_description($log);

    my $final_file = $output_bam ? $bam_file : $shrimpFile;
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    if ($output_bam) {
      print $pbs "gmapper -L $shrimp2_index $sampleFile $mirna $option --extra-sam-fields > $sam_file

if [ -s $sam_file ]; then
  samtools view -S -b $sam_file | samtools sort - $sample_name
  samtools index $bam_file 
  samtools flagstat $bam_file > ${bam_file}.stat 
fi

echo finished=`date`
";
    }
    else {
      print $pbs "gmapper -L $shrimp2_index $sampleFile $mirna $option --pretty >$shrimpFile
      
echo finished=`date` 
";
    }

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $output_bam = $config->{$section}{output_bam} or die "define ${section}::output_bam first";

  my %raw_files = %{ get_raw_files( $config, $section, "source", ".fastq\$" ) };

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my $cur_dir      = $result_dir . "/$sample_name/";
    my @result_files = ();
    if ($output_bam) {
      push( @result_files, $cur_dir . $sample_name . ".bam" );
    }
    else {
      push( @result_files, $cur_dir . $sample_name . ".shrimp" );
    }

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
