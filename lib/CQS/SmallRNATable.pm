#!/usr/bin/perl
package CQS::SmallRNATable;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_srt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );
  my $is_tRH = get_option( $config, $section, "is_tRH", 0 );

  my $python_script = dirname(__FILE__) . "/smallRNATable.py";
  if ( !-e $python_script ) {
    die "File not found : " . $python_script;
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $noCategory = $option =~ /noCategory/;

  if ( defined $config->{$section}{groups} || defined $config->{$section}{groups_ref} ) {
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $group_name ( sort keys %{$groups} ) {
      my $filelist   = $self->get_file( $pbs_dir,    "${task_name}_${group_name}", ".filelist", 0 );
      my $outputfile = $self->get_file( $result_dir, "${task_name}_${group_name}", ".count",    0 );
      my $outputname = basename($outputfile);

      my $ntaFile     = basename( $self->get_file( $result_dir, "${task_name}_${group_name}", ".miRNA.NTA.count",      0 ) );
      my $ntaBaseFile = basename( $self->get_file( $result_dir, "${task_name}_${group_name}", ".miRNA.NTA.base.count", 0 ) );

      my @samples = @{ $groups->{$group_name} };
      open( my $fl, ">$filelist" ) or die "Cannot create $filelist";
      for my $sample_name ( sort @samples ) {
        my @count_files = @{ $raw_files{$sample_name} };
        my $countFile   = $count_files[0];
        print $fl $sample_name, "\t", $countFile, "\n";
      }
      close($fl);

      my $pythonCode = $noCategory ? "" : ($is_tRH? "": "python3 $python_script \"$ntaFile\" \"$ntaBaseFile\"");
      print $pbs "
if [ ! -s $outputname ]; then
  cqstools smallrna_table $option -o $outputname -l $filelist
  $pythonCode
fi
";
    }
    $self->close_pbs( $pbs, $pbs_file );
  }
  else {
    my $filelist   = $self->get_file( $pbs_dir,    ${task_name}, ".filelist", 0 );
    my $outputfile = $self->get_file( $result_dir, ${task_name}, ".count",    0 );
    my $outputname = basename($outputfile);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $outputname );

    open( my $fl, ">$filelist" ) or die "Cannot create $filelist";
    for my $sample_name ( sort keys %raw_files ) {
      my @count_files = @{ $raw_files{$sample_name} };
      my $countFile   = $count_files[0];
      print $fl $sample_name, "\t", $countFile, "\n";
    }
    close($fl);

    my $ntaFile     = basename( $self->get_file( $result_dir, ${task_name}, ".miRNA.NTA.count",      0 ) );
    my $ntaBaseFile = basename( $self->get_file( $result_dir, ${task_name}, ".miRNA.NTA.base.count", 0 ) );

    my $pythonCode = $noCategory ? "" : ($is_tRH?"": "python3 $python_script \"$ntaFile\" \"$ntaBaseFile\"");
    print $pbs "
cqstools smallrna_table $option -o $outputname -l $filelist
$pythonCode
";
    $self->close_pbs( $pbs, $pbs_file );
  }
}

sub addOutput {
  my ( $self, $result_files, $result_dir, $pbs_dir, $key, $option, $is_tRH ) = @_;
  my $noCategory = $option =~ /noCategory/;
  push( @$result_files, $self->get_file( $result_dir, $key, ".count", 0 ) );
  push( @$result_files, $self->get_file( $result_dir, $key, ".read.count", 0 ) );

  if ( !$noCategory ) {
    if ( !$is_tRH ) {
      push( @$result_files, $self->get_file( $result_dir, $key, ".miRNA.count",            0 ) );
      push( @$result_files, $self->get_file( $result_dir, $key, ".miRNA.count.position",   0 ) );
      push( @$result_files, $self->get_file( $result_dir, $key, ".miRNA.read.count",       0 ) );
      push( @$result_files, $self->get_file( $result_dir, $key, ".miRNA.isomiR.count",     0 ) );
      push( @$result_files, $self->get_file( $result_dir, $key, ".miRNA.isomiR_NTA.count", 0 ) );
      push( @$result_files, $self->get_file( $result_dir, $key, ".miRNA.NTA.count",        0 ) );
      push( @$result_files, $self->get_file( $result_dir, $key, ".miRNA.NTA.base.count",   0 ) );
    }

    push( @$result_files, $self->get_file( $result_dir, $key, ".tRNA.count",                    0 ) );
    push( @$result_files, $self->get_file( $result_dir, $key, ".tRNA.count.position",           0 ) );
    push( @$result_files, $self->get_file( $result_dir, $key, ".tRNA.read.count",               0 ) );
    push( @$result_files, $self->get_file( $result_dir, $key, ".tRNA.aminoacid.count",          0 ) );
    push( @$result_files, $self->get_file( $result_dir, $key, ".tRNA.aminoacid.count.position", 0 ) );
    push( @$result_files, $self->get_file( $result_dir, $key, ".tRNA.count.startpos", 0 ) );

    if ( !$is_tRH ) {

      if ( $option =~ /exportYRNA/ ) {
        push( @$result_files, $self->get_file( $result_dir, $key, ".yRNA.count",          0 ) );
        push( @$result_files, $self->get_file( $result_dir, $key, ".yRNA.count.position", 0 ) );
        push( @$result_files, $self->get_file( $result_dir, $key, ".yRNA.read.count",     0 ) );
      }
      if ( $option =~ /exportSnRNA/ ) {
        push( @$result_files, $self->get_file( $result_dir, $key, ".snRNA.count",          0 ) );
        push( @$result_files, $self->get_file( $result_dir, $key, ".snRNA.count.position", 0 ) );
        push( @$result_files, $self->get_file( $result_dir, $key, ".snRNA.read.count",     0 ) );
      }
      if ( $option =~ /exportSnoRNA/ ) {
        push( @$result_files, $self->get_file( $result_dir, $key, ".snoRNA.count",          0 ) );
        push( @$result_files, $self->get_file( $result_dir, $key, ".snoRNA.count.position", 0 ) );
        push( @$result_files, $self->get_file( $result_dir, $key, ".snoRNA.read.count",     0 ) );
      }
      push( @$result_files, $self->get_file( $result_dir, $key, ".rRNA.count",       0 ) );
      push( @$result_files, $self->get_file( $result_dir, $key, ".rRNA.read.count",  0 ) );
      push( @$result_files, $self->get_file( $result_dir, $key, ".other.count",      0 ) );
      push( @$result_files, $self->get_file( $result_dir, $key, ".other.read.count", 0 ) );
    }
  }
  push( @$result_files, $self->get_file( $pbs_dir, $key, ".filelist", 0 ) );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $is_tRH = get_option( $config, $section, "is_tRH", 0 );

  my $result = {};

  my @result_files = ();
  if ( defined $config->{$section}{groups} || defined $config->{$section}{groups_ref} ) {
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $group_name ( sort keys %{$groups} ) {
      $self->addOutput( \@result_files, $result_dir, $pbs_dir, "${task_name}_${group_name}", $option, $is_tRH );
    }
  }
  else {
    $self->addOutput( \@result_files, $result_dir, $pbs_dir, $task_name, $option, $is_tRH );
  }
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
