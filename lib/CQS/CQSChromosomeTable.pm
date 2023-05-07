#!/usr/bin/perl
package CQS::CQSChromosomeTable;

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
  $self->{_suffix} = "_ct";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $force_species_file = get_option($config, $section, "force_species_file", 0);

  if ( defined $config->{$section}{groups} || defined $config->{$section}{groups_ref} ) {
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $group_name ( sort keys %{$groups} ) {
      my @samples = @{ $groups->{$group_name} };

      my $filelist   = $self->get_file( $pbs_dir,    "${task_name}_${group_name}", ".filelist", 0 );
      my $outputfile = $self->get_file( $result_dir, "${task_name}_${group_name}", ".count",    0 );
      my $outputname = basename($outputfile);

      my $specics_file = basename($self->get_file( $result_dir, "${task_name}_${group_name}", ".Species.count", 0 ));

      open( FL, ">$filelist" ) or die "Cannot create $filelist";
      for my $sample_name ( sort @samples ) {
        my @count_files = @{ $raw_files{$sample_name} };
        my $countFile   = $count_files[0];
        print FL $sample_name, "\t", $countFile, "\n";
      }
      close(FL);

      print $pbs "
cqstools chromosome_table $option -o $outputname -l $filelist

";
      if($force_species_file){
        print $pbs "
  if [[ ! -s $specics_file && -s $outputname ]]; then
    cp $outputname $specics_file
  fi
  ";
      }
    }

    $self->close_pbs( $pbs, $pbs_file );
  }
  else {
    my $filelist   = $self->get_file( $pbs_dir,    ${task_name}, ".filelist", 0 );
    my $outputfile = $self->get_file( $result_dir, ${task_name}, ".count",    0 );
    my $outputname = basename($outputfile);

    my $specics_file = basename($self->get_file( $result_dir, $task_name, ".Species.count", 0 ));

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $outputname );

    open( FL, ">$filelist" ) or die "Cannot create $filelist";
    for my $sample_name ( sort keys %raw_files ) {
      my @count_files = @{ $raw_files{$sample_name} };
      my $countFile   = $count_files[0];
      print FL $sample_name, "\t", $countFile, "\n";
    }
    close(FL);

    print $pbs "
cqstools chromosome_table $option -o $outputname -l $filelist

";

    if($force_species_file){
      print $pbs "
if [[ ! -s $specics_file && -s $outputname ]]; then
  cp $outputname $specics_file
fi
";
    }
    $self->close_pbs( $pbs, $pbs_file );
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result         = {};
  my $has_category   = $option =~ /--categoryMapFile/;
  my $has_read_count = $option =~ /--outputReadTable/;

  my $has_contig = $option =~ /--outputReadContigTable/;

  my @result_files = ();
  if ( defined $config->{$section}{groups} || defined $config->{$section}{groups_ref} ) {
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $group_name ( sort keys %{$groups} ) {
      push( @result_files, $self->get_file( $result_dir, "${task_name}_${group_name}", ".count",     0 ) );
      push( @result_files, $self->get_file( $pbs_dir, "${task_name}_${group_name}", ".filelist", 0 ) );
      push( @result_files, $self->get_file( $result_dir, "${task_name}_${group_name}", ".count.xml", 0 ) );
      if ($has_category) {
        push( @result_files, $self->get_file( $result_dir, "${task_name}_${group_name}", ".Species.count", 0 ) );
      }
      if ($has_read_count) {
        push( @result_files, $self->get_file( $result_dir, "${task_name}_${group_name}", ".read.count", 0 ) );
      }
      if ($has_contig) {
        push( @result_files, $self->get_file( $result_dir, "${task_name}_${group_name}", ".contig.count",               0 ) );
        push( @result_files, $self->get_file( $result_dir, "${task_name}_${group_name}", ".contig.count.details",       0 ) );
        push( @result_files, $self->get_file( $result_dir, "${task_name}_${group_name}", ".contig.count.details.depth", 0 ) );
      }
    }
  }
  else {
    push( @result_files, $self->get_file( $result_dir, $task_name, ".count", 0 ) );
    push( @result_files, $self->get_file( $pbs_dir, $task_name, ".filelist", 0 ) );
    push( @result_files, $self->get_file( $result_dir, $task_name, ".count.xml", 0 ) );
    if ($has_category) {
      push( @result_files, $self->get_file( $result_dir, $task_name, ".Species.count", 0 ) );
    }
    if ($has_read_count) {
      push( @result_files, $self->get_file( $result_dir, $task_name, ".read.count", 0 ) );
    }
    if ($has_contig) {
      push( @result_files, $self->get_file( $result_dir, $task_name, ".contig.count",               0 ) );
      push( @result_files, $self->get_file( $result_dir, $task_name, ".contig.count.details",       0 ) );
      push( @result_files, $self->get_file( $result_dir, $task_name, ".contig.count.details.depth", 0 ) );
    }
  }

  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
