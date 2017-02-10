#!/usr/bin/perl
package CQS::Task;

use strict;
use warnings;
use CQS::ConfigUtils;
use Data::Dumper;

sub new {
  my ($class) = @_;
  my $self = { _name => __PACKAGE__, _suffix => "", _task_prefix => "", _task_suffix => "", _pbskey => "source" };
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

sub get_clear_map {
  my $self   = shift;
  my $result = $self->result(@_);
  for my $key ( keys %$result ) {
    my $values = $result->{$key};
    my @newvalues = grep { !/\/pbs\// } @$values;
    if ( scalar(@newvalues) > 0 ) {
      $result->{$key} = \@newvalues;
    }
    else {
      $result->{$key} = undef;
    }
  }
  return $result;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  #print  "task_name = " . $task_name . "\n";

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};
  if ( $self->{_pbskey} eq "" ) {
    $result->{$task_name} = $self->get_pbs_filename( $pbs_dir, $task_name );
  }
  else {
    my %fqFiles = %{ get_raw_files( $config, $section, $self->{_pbskey} ) };

    for my $sample_name ( sort keys %fqFiles ) {
      $result->{$sample_name} = $self->get_pbs_filename( $pbs_dir, $sample_name );
    }
  }

  return $result;
}

sub html {
  my ( $self, $config, $section ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );
  my $result = $self->result( $config, $section );

  my $htmlfile = $target_dir . "/" . $section . ".html";
  open( my $html, ">$htmlfile" ) or die "Cannot create $htmlfile";
  print $html "<HTML>
<HEAD>
<TITLE>$task_name : $section</TITLE>
</HEAD>
<BODY>
<h4>$task_name : $section</h4>
<table>
";

  for my $expect_name ( sort keys %{$result} ) {
    my $expect_files = $result->{$expect_name};
    print "<tr><td>$expect_name</td><td>\n";
    for my $expect_file ( @{$expect_files} ) {
      my $commonLen = 0;
      ( $target_dir ^ $expect_file ) =~ /^(\0*)/;
      $commonLen = $+[0];
      substr $expect_file, 0, $commonLen, "";
      print "$expect_file\n";
    }
    print "</td></tr>\n";
  }

  print $html "</table>  
</BODY>
</HTML>
";

  return ($htmlfile);
}

sub require {
  my $result = [];
  return $result;
}

sub get_name {
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

sub get_file {
  my ( $self, $dir, $name, $extension, $hassuffix ) = @_;
  if ( !defined $extension ) {
    $extension = "";
  }
  if ( !defined $hassuffix ) {
    $hassuffix = 1;
  }

  return $dir . "/" . $self->get_name( $name, $extension, $hassuffix );
}

sub pbs_name {
  my ( $self, $sample_name ) = @_;
  return $self->get_name( $sample_name, ".pbs" );
}

sub get_pbs_filename {
  my ( $self, $dir, $sample_name ) = @_;
  return $self->get_file( $dir, $sample_name, ".pbs" );
}

sub get_log_filename {
  my ( $self, $dir, $sample_name ) = @_;
  return $self->get_file( $dir, $sample_name, ".log" );
}

sub get_task_filename {
  my ( $self, $dir, $task_name ) = @_;
  return $self->get_file( $dir, $task_name, ".sh" );
}

sub open_pbs {
  my ( $self, $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command, $final_file_can_empty ) = @_;

  if ( !defined $init_command ) {
    $init_command = "";
  }

  my $module_name = $self->{_name};

  open( my $pbs, ">$pbs_file" ) or die $!;

  print $pbs "$pbs_desc
$log_desc

$path_file
$init_command

cd $result_dir

";
  if ( defined $final_file ) {
    my $delete_file = ( $final_file =~ /^\// ) ? $final_file : "${result_dir}/${final_file}";

    my $checkFile = $final_file_can_empty ? "-e" : "-s";
    print $pbs "
if [[ ( $checkFile $final_file ) || ( -d $final_file ) ]]; then
  echo job has already been done. if you want to do again, delete $delete_file and submit job again.
  exit 0
fi
";
  }

  print $pbs "
echo ${module_name}_start=`date`
echo working in $result_dir ...
 
";
  return $pbs;
}

sub close_pbs {
  my ( $self, $pbs, $pbs_file ) = @_;

  my $module_name = $self->{_name};

  print $pbs "

echo ${module_name}_end=`date`

exit 0
 
";

  close $pbs;

  print "$pbs_file created. \n";
}

1;
