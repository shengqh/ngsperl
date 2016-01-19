#!/usr/bin/perl
package Variants::GlmvcAnnotation;

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
  $self->{_suffix} = "_ga";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $glmvcfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1 );

  my $rnaediting_db = get_directory( $config, $section, "rnaediting_db", 0 );
  if ( defined $rnaediting_db ) {
    $option = $option . " --rnaediting_db $rnaediting_db ";
  }

  my $annovar_buildver = $config->{$section}{annovar_buildver};
  if ( defined $annovar_buildver ) {
    $option = $option . " --annovar_buildver $annovar_buildver ";

    my $annovar_protocol = $config->{$section}{annovar_protocol};
    if ( defined $annovar_protocol ) {
      $option = $option . " --annovar_protocol $annovar_protocol ";
    }

    my $annovar_operation = $config->{$section}{annovar_operation};
    if ( defined $annovar_operation ) {
      $option = $option . " --annovar_operation $annovar_operation ";
    }

    my $annovar_db = $config->{$section}{annovar_db};
    if ( defined $annovar_db ) {
      $option = $option . " --annovar_db $annovar_db ";
    }
  }

  my $distance_exon_gtf = get_param_file( $config->{$section}{distance_exon_gtf}, "distance_exon_gtf", 0 );
  if ( defined $distance_exon_gtf ) {
    $option = $option . " --distance_exon_gtf $distance_exon_gtf ";
  }

  my $anno = defined $rnaediting_db || defined $annovar_buildver || defined $annovar_buildver;

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";
  print $sh "cd $pbs_dir \n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile  = $sample_files[0];
    my $cur_dir      = create_directory_or_die( $result_dir . "/$sample_name" );

    my $pbs_file  = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name  = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);
    my $final    = "${sample_name}.annotation.tsv";

    print $sh "\$MYCMD ./$pbs_name \n";

    open( my $out, ">$pbs_file" ) or die $!;
    print $out "$pbs_desc
$log_desc

$path_file 

echo Glmvc=`date` 

cd $cur_dir

if [ -s $final ]; then
  echo job has already been done. if you want to do again, delete ${cur_dir}/${final} and submit job again.
  exit 0;
fi      
      
mono $glmvcfile annotation $option -i $sampleFile -o ${final}

echo finished=`date`
";

    close $out;

    print "$pbs_file created \n";

  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $raw_files = get_raw_files( $config, $section );

  my $result = {};
  for my $sample_name ( keys %{$raw_files} ) {
    my @result_files = ();
    my $cur_dir      = $result_dir . "/$sample_name";
    push( @result_files, "$cur_dir/${sample_name}.annotation.tsv" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}
1;
