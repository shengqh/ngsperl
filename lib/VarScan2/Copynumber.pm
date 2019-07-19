#!/usr/bin/perl
package VarScan2::Copynumber;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_cp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $rtemplate = dirname(__FILE__) . "/Copynumber.r";
  if ( !-e $rtemplate ) {
    die "File not found : " . $rtemplate;
  }

  my $varscan2_jar = get_param_file( $config->{$section}{VarScan2_jar}, "VarScan2_jar", 1, $self->using_docker() );
  my $faFile       = get_param_file( $config->{$section}{fasta_file},   "fasta_file",   1 );

  my $mpileup_options = get_option( $config, $section, "mpileup_options", "-q 1" );
  my $call_options    = get_option( $config, $section, "call_options",    "--amp-threshold 0.3 --del-threshold 0.3" );

  my %group_sample_map = %{ get_group_sample_map( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";
  print $sh "cd $pbs_dir\n";

  my $java_option = get_option( $config, $section, "java_option", "" );

  for my $group_name ( sort keys %group_sample_map ) {
    my @sample_files = @{ $group_sample_map{$group_name} };
    my $sampleCount  = scalar(@sample_files);

    if ( $sampleCount != 2 ) {
      die "SampleFile should be normal,tumor paired.";
    }

    my $cur_dir = create_directory_or_die( $result_dir . "/$group_name" );

    my $normal = $sample_files[0][1];
    my $tumor  = $sample_files[1][1];

    my $cpRawFile  = "${group_name}.copynumber";
    my $cpCallFile = "${group_name}.call";
    my $cpSegFile  = "${group_name}.segment";

    my $rfile = $cur_dir . "/${group_name}.r";
    open( RF, ">$rfile" ) or die "Cannot create $rfile";
    print RF "setwd(\"$cur_dir\")  
inputfile<-\"$cpCallFile\"
outputfile<-\"$cpSegFile\"
";
    open RT, "<$rtemplate" or die $!;
    while (<RT>) {
      print RF $_;
    }
    close(RT);
    close(RF);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

    print $pbs "
if [ ! -s $cpRawFile ]; then
  if [ ! -s ${normal}.bai ]; then
    samtools index ${normal}
  fi

  if [ ! -s ${tumor}.bai ]; then
    samtools index ${tumor}
  fi

  samtools mpileup $mpileup_options -f $faFile $normal $tumor | awk 'NF==9 && \$4!=0' | java $java_option -jar $varscan2_jar copynumber - $group_name --mpileup 1 $option  
fi

if [[ -s $cpRawFile && ! -s $cpCallFile ]]; then
  java $java_option -jar $varscan2_jar copyCaller $cpRawFile $call_options --output-file $cpCallFile 
fi

if [[ -s $cpCallFile && ! -s $cpSegFile ]]; then
  R --vanilla -f $rfile
fi

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $group_name ( keys %{$groups} ) {
    my @result_files = ();
    my $cur_dir      = $result_dir . "/$group_name";
    push( @result_files, "$cur_dir/${group_name}.copynumber" );
    push( @result_files, "$cur_dir/${group_name}.call" );
    push( @result_files, "$cur_dir/${group_name}.segment" );
    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
