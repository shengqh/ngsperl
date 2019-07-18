#!/usr/bin/perl
package CQS::Task;

use strict;
use warnings;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use Data::Dumper;

sub new {
  my ($class) = @_;
  my $self = { _name => __PACKAGE__, _suffix => "", _task_prefix => "", _task_suffix => "", _pbskey => "source", _docker_prefix => "" };
  bless $self, $class;
  return $self;
}

sub name {
  my ($self) = @_;
  return $self->{_name};
}

sub perform {
}

sub get_final_file {
  my ( $self, $config, $section, $result_dir, $sample ) = @_;
  my $expect = $self->result($config, $section);
  
  if(not defined $sample){
    my @samples = sort keys %$expect;
    @samples = sort {$b cmp $a} @samples;
    $sample = $samples[0];
  }
  
  my $final_files_ref = $expect->{$sample};
  my @final_files = @$final_files_ref;
  my $result = $final_files[-1];
  
  if (rindex($result, $result_dir) == 0){
    $result = substr($result, length($result_dir));
    my $firstChar = substr($result, 0, 1);
    if (($firstChar eq '/') or ($firstChar eq '\\')){
      $result = substr($result, 1);
    }
  }
  return($result);
}

sub result {
  my $result = {};
  return $result;
}

sub acceptSample {
  my ( $self, $config, $section, $sampleName ) = @_;
  return (1);
}

sub can_result_be_empty_file {
  my ( $self, $config, $section ) = @_;
  my $curSection = get_config_section( $config, $section );
  if ( defined $curSection->{can_result_be_empty_file} ) {
    return $curSection->{can_result_be_empty_file};
  }
  return 0;
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

sub get_pbs_key {
  my ( $self, $config, $section ) = @_;
  return $self->{_pbskey};
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  #print  "task_name = " . $task_name . "\n";

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};

  my $pbsKey = $self->get_pbs_key( $config, $section );
  if ( $pbsKey eq "" ) {
    $result->{$task_name} = $self->get_pbs_filename( $pbs_dir, $task_name );
  }
  else {
    my %fqFiles = %{ get_raw_files( $config, $section, $pbsKey ) };

    for my $sample_name ( sort keys %fqFiles ) {
      if ( $self->acceptSample( $config, $section, $sample_name ) ) {
        $result->{$sample_name} = $self->get_pbs_filename( $pbs_dir, $sample_name );
      }
    }
  }

  return $result;
}

#get pbs source map which indicates which sample name the pbs file comes from
#for impute2 which generate multiple pbs files and multiple result files from 1 sample name,
#the multiple pbs should be mapped to same sample name
sub get_pbs_source {
  my ( $self, $config, $section ) = @_;

  my $pbsFiles = $self->get_pbs_files( $config, $section );
  my $result = {};
  for my $resKey ( keys %$pbsFiles ) {
    $result->{ $pbsFiles->{$resKey} } = [$resKey];
  }
  return $result;
}

#get result pbs map which indicates which pbs the result name related.
#for bed file split which has 1 pbs and multiple result files, the multiple result names
#should be mapped to same pbs
sub get_result_dependent_pbs {
  my ( $self, $config, $section ) = @_;

  return $self->get_pbs_files( $config, $section );
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

sub do_get_docker_value {
  my ($self, $keyName) = @_;
  my $result = undef;

  if ( defined $self->{_config} ) {
    if ( ( defined $self->{_section} ) and ( defined $self->{_config}{ $self->{_section} }{$keyName} ) ) {
      $result = $self->{_config}{ $self->{_section} }{$keyName};
      if ( defined $result ) {
        return ($result);
      }
    }

    if ( ( defined $self->{_config} ) and ( defined $self->{_config}{general} ) and ( defined $self->{_config}{general}{$keyName} ) ) {
      $result = $self->{_config}{general}{$keyName};
      if ( defined $result ) {
        return ($result);
      }
    }
  }
  
  return ($result);
}

sub get_docker_value {
  my ($self, $required) = @_;
  my $command = undef;
  my $init = undef;

  my $commandKey = $self->{_docker_prefix} . "docker_command";
  my $initKey = $self->{_docker_prefix} . "docker_init";

  $command = $self->do_get_docker_value($commandKey);
  if (defined $command){
    $init = $self->do_get_docker_value($initKey);
    return ($command, $init);
  }
  
  my $baseCommandKey = "docker_command";
  my $baseInitKey = "docker_init";

  $command = $self->do_get_docker_value($baseCommandKey);
  if (defined $command){
    $init = $self->do_get_docker_value($baseInitKey);
    return ($command, $init);
  }
  
  if(defined $required and $required){
    die "Define $commandKey for task " . $self->{_name};
  }
  
  return (undef, undef);
}

sub open_pbs {
  my ( $self, $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command, $final_file_can_empty, $input_file ) = @_;

  if ( !defined $init_command ) {
    $init_command = "";
  }

  my $module_name = $self->{_name};

  open( my $pbs, ">$pbs_file" ) or die $!;

  print $pbs "$pbs_desc
$log_desc

$path_file
$init_command

cd \"$result_dir\"

";
  if ( defined $final_file ) {
    my $delete_file = ( $final_file =~ /^\// ) ? $final_file : "${result_dir}/${final_file}";

    my $checkFile = $final_file_can_empty ? "-e" : "-s";
    print $pbs "
if [[ !(1 -eq \$1) ]]; then
  if [[ ( $checkFile $final_file ) || ( -d $final_file ) ]]; then
    echo job has already been done. if you want to do again, delete $delete_file and submit job again.
    exit 0
  fi
fi
";
  }

  if ( defined $input_file ) {
    print $pbs "
if [[ ! -s $input_file ]]; then
  echo input file not exist: $input_file
  exit 1
fi
";
  }

  print $pbs "
echo ${module_name}_start=`date`
echo working in $result_dir ...
 
";

  my ($docker_command, $docker_init) = $self->get_docker_value();
  my $is_sequenceTask = ( $module_name =~ /SequenceTask/ );
  if ( ( defined $docker_command ) and ( not $is_sequenceTask ) ) {
    if(not defined $docker_init){
      $docker_init = "";
    }
    my $sh_file = $pbs_file . ".sh";

    print $pbs "
export R_LIBS=
export PYTHONPATH=
export JAVA_HOME=
 
$docker_command bash $sh_file 

echo ${module_name}_end=`date`

exit 0

";
    close $pbs;
    open( $pbs, ">$sh_file" ) or die $!;
    print $pbs "
export HOME=$result_dir

cd $result_dir

$docker_init
";

  }

  return $pbs;
}

sub close_pbs {
  my ( $self, $pbs, $pbs_file ) = @_;

  my $module_name = $self->{_name};

  my $docker_command = $self->get_docker_value();

  if ( not defined $docker_command ) {
    print $pbs "

echo ${module_name}_end=`date`

exit 0
 
";
  }

  close $pbs;

  print "$pbs_file created. \n";
}

sub get_java_option {
  my ( $self, $config, $section, $memory ) = @_;
  my $result = $config->{$section}{java_option};
  if ( !defined $result || $result eq "" ) {
    $result = "-Xmx${memory}";
  }
  return ($result);
}

1;
