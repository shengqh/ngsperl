#!/usr/bin/perl
package CQS::Task;

use strict;
use warnings;
use CQS::ConfigUtils;

sub new {
  my ($class) = @_;
  my $self = { _name => undef, _suffix => "", _task_prefix => "", _task_suffix => "", _pbskey => "source" };
  print __PACKAGE__;
  bless $self, $class;
  return $self;
}

sub name {
  my ($self) = @_;
  return $self->{_name};
}

sub perform {
}

sub result {
  my $result = {};
  return $result;
}

sub require {
  my $result = [];
  return $result;
}

sub getname {
  my ( $self, $name, $extension, $hassuffix ) = @_;
  if ( !defined $extension ) {
    $extension = "";
  }
  if ( !defined $hassuffix ) {
    $hassuffix = 1;
  }

  if ($hassuffix) {
    return $self->{_task_prefix} . $name . $self->{_task_suffix} . $self->{_suffix} . $extension;
  }
  else {
    return $self->{_task_prefix} . $name . $self->{_task_suffix} . $extension;
  }
}

sub getfile {
  my ( $self, $dir, $name, $extension, $hassuffix ) = @_;
  if ( !defined $extension ) {
    $extension = "";
  }
  if ( !defined $hassuffix ) {
    $hassuffix = 1;
  }

  return $dir . "/" . $self->getname( $name, $extension, $hassuffix );
}

sub pbs_name {
  my ( $self, $sample_name ) = @_;
  return $self->getname( $sample_name, ".pbs" );
}

sub pbs_file {
  my ( $self, $dir, $sample_name ) = @_;
  return $self->getfile( $dir, $sample_name, ".pbs" );
}

sub log_file {
  my ( $self, $dir, $sample_name ) = @_;
  return $self->getfile( $dir, $sample_name, ".log" );
}

sub task_file {
  my ( $self, $dir, $task_name ) = @_;
  return $self->getfile( $dir, $task_name, ".sh" );
}

sub pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_description, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  #print  "task_name = " . $task_name . "\n";

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};
  if ( $self->{_pbskey} eq "" ) {
    $result->{$task_name} = $self->pbs_file( $pbs_dir, $task_name );
  }
  else {
    my %fqFiles = %{ get_raw_files( $config, $section, $self->{_pbskey} ) };

    for my $sample_name ( sort keys %fqFiles ) {
      $result->{$sample_name} = $self->pbs_file( $pbs_dir, $sample_name );
    }
  }

  return $result;
}

sub write_pbs_start {
  my ( $self, $pbs, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $module_name ) = @_;

  print $pbs "$pbs_desc
$log_desc

$path_file

cd $result_dir

if [ -s $final_file ]; then
  echo job has already been done. if you want to do again, delete $final_file and submit job again.
  exit 0
fi

echo ${module_name}_start=`date`
 
";
}

sub write_pbs_end {
  my ( $self, $pbs, $module_name ) = @_;

  print $pbs "
echo ${module_name}_end=`date`

exit 1
 
";
}

1;
