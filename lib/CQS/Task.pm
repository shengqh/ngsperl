#!/usr/bin/perl
package CQS::Task;

use strict;
use warnings;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::StringUtils;
use CQS::NGSCommon;
use File::Basename;
use Data::Dumper;
use List::MoreUtils qw(uniq);

sub new {
  my ($class) = @_;

  my $self = {
    _name          => __PACKAGE__,
    _suffix        => "",
    _task_prefix   => "",
    _task_suffix   => "",
    _pbskey        => "source",
    _docker_prefix => "",
    _can_use_docker => 1,
    _forbid_tmp_folder  => 0,
    _use_tmp_folder => 0,
    _localize_to_local_folder => 0,
    _final_file_in_last => 1,
  };
  bless $self, $class;
  return $self;
}

sub init_docker_prefix {
  my ($self, $package) = @_;
  my $docker_prefix = $package;
  $docker_prefix =~ s/.+://g;
  $docker_prefix = $docker_prefix . "_";
  $self->{_docker_prefix} = $docker_prefix;
}

sub name {
  my ($self) = @_;
  return $self->{_name};
}

sub perform {
}


sub get_absolute_final_file {
  my ( $self, $config, $section, $sample ) = @_;
  my $expect = $self->result( $config, $section, "(?<!version)\$", 1 );

  my @samples = sort keys %$expect;
  @samples = sort { $b cmp $a } @samples;

  if ( not defined $sample ) {
    $sample  = $samples[0];
  }

  if (not defined $expect->{$sample}){
    $sample  = $samples[0];
  }

  my $final_files_ref = $expect->{$sample};
  my @final_files     = @$final_files_ref;

  my @no_filelists = grep(!/.filelist$/, @final_files);
  if(scalar(@no_filelists) == 0){
    die "Cannot find final file of $sample in section $section";
  }

  my $result          = $self->{_final_file_in_last} ? $no_filelists[-1] :  $no_filelists[0];

  return ($result);
}

sub get_final_file {
  my ( $self, $config, $section, $result_dir, $sample ) = @_;

  my $result          = $self->get_absolute_final_file($config, $section, $sample);

  if ( rindex( $result, $result_dir ) == 0 ) {
    $result = substr( $result, length($result_dir) );
    my $firstChar = substr( $result, 0, 1 );
    if ( ( $firstChar eq '/' ) or ( $firstChar eq '\\' ) ) {
      $result = substr( $result, 1 );
    }
  }

  return ($result);
}

sub get_result_files {
  my ( $self, $config, $section, $result_dir, $sample_name ) = @_;

  die "Override get_result_files of " . $self->{_name} . " first.";
}

sub result {
  my ( $self, $config, $section, $pattern, $removeEmpty ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section );

  my $result = {};

  my $raw_files = get_raw_files($config, $section);
  for my $sample_name (sort keys %$raw_files){
    if ( $self->acceptSample( $config, $section, $sample_name ) ) {
      my $result_files = $self->get_result_files( $config, $section, $result_dir, $sample_name );
      my $filtered_files = filter_array( $result_files, $pattern, 1 );
      $result->{$sample_name} = $filtered_files;
    }
  }

  return $result;
}

sub acceptSample {
  my ( $self, $config, $section, $sampleName ) = @_;
  return (1);
}

sub can_result_be_empty_file {
  my ( $self, $config, $section, $filename ) = @_;
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
    my $values    = $result->{$key};
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

sub init_tmp_folder {
  my ( $self, $pbs, $result_dir ) = @_;
  if( $self->{"_use_tmp_folder"}){
    print $pbs "
res_dir='$result_dir'
tmp_dir=\$(mktemp -d -t ci-\$(date +\%Y-\%m-\%d-\%H-\%M-\%S)-XXXXXXXXXX)

tmp_cleaner()
{
rm -rf \${tmp_dir}
exit -1
}
trap 'tmp_cleaner' TERM

echo using tmp_dir=\$tmp_dir
cd \$tmp_dir

";
    return(1);
  }else{
    return(0);
  }
}

sub clean_tmp_folder {
  my ( $self, $pbs ) = @_;
  if( $self->{"_use_tmp_folder"}){
    print $pbs "
if [[ -d \$tmp_dir && \$tmp_dir != '/' ]]; then
  echo copy result from \$tmp_dir to \$res_dir
  #if the pbs was generated again during task is running, copy may be unpredictable. 
  #make sure to change to tmp_dir before copy result

  cd \$tmp_dir
  cd \$tmp_dir
  cd \$tmp_dir
  cd \$tmp_dir
  cd \$tmp_dir
  cd \$tmp_dir
  cd \$tmp_dir
  cd \$tmp_dir
  cd \$tmp_dir
  cd \$tmp_dir

  cp -p -r * \$res_dir
  cd \$res_dir
  echo delete tmp folder \$tmp_dir
  rm -rf \$tmp_dir
  echo move file and clean tmp folder done.
fi
";
  }
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

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
  my $result   = {};
  for my $resKey ( keys %$pbsFiles ) {
    $result->{ $pbsFiles->{$resKey} } = [$resKey];
  }
  return $result;
}

#get result pbs map which indicates which pbs the result name related.
#for bed file split which has 1 pbs and multiple result files, the multiple result names
#should be mapped to same pbs
sub get_result_pbs {
  my ( $self, $config, $section ) = @_;

  return $self->get_pbs_files( $config, $section );
}

#get the task dependent pbs list map
sub get_dependent_pbs_map {
  my ( $self, $config, $section ) = @_;

  my $result = {};

  my $task_section = $config->{$section};
  for my $key ( keys %$task_section ) {
    my $mapname = $key;
    if ( $mapname !~ /_ref$/ ) {
      next;
    }

    $mapname =~ s/_config_ref//g;
    $mapname =~ s/_ref//g;
    my $refpbsmap = get_ref_section_pbs( $config, $section, $mapname );
    for my $refkey ( keys %$refpbsmap ) {
      my $refpbs = $refpbsmap->{$refkey};
      my $curpbs = $result->{$refkey};
      if ( !defined $curpbs ) {
        $curpbs = {};
      }
      for my $eachrefpbs (@$refpbs) {
        $curpbs->{$eachrefpbs} = 1;
      }

      $result->{$refkey} = $curpbs;
    }
  }
  return ($result);
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
  my ( $self, $keyName ) = @_;
  my $result = undef;

  if ( defined $self->{_config} ) {
    if (  ( defined $self->{_section} )
      and ( defined $self->{_config}{ $self->{_section} }{$keyName} ) )
    {
      $result = $self->{_config}{ $self->{_section} }{$keyName};
      if ( defined $result ) {
        return ($result);
      }
    }

    if (  ( defined $self->{_config} )
      and ( defined $self->{_config}{general} )
      and ( defined $self->{_config}{general}{$keyName} ) )
    {
      $result = $self->{_config}{general}{$keyName};
      if ( defined $result ) {
        return ($result);
      }
    }
  }

  return ($result);
}

sub can_use_docker(){
  my ($self) = @_;
  return($self->{_can_use_docker});
}

sub using_docker {
  my ($self) = @_;

  if (not $self->{_can_use_docker}){
    return(0);
  }

  my ( $docker_command, $docker_init ) = $self->get_docker_value();
  my $is_sequenceTask = ( $self->{_name} =~ /SequenceTask/ );
  return ( ( defined $docker_command ) and ( not $is_sequenceTask ) );
}

sub get_docker_value {
  my ( $self, $required ) = @_;
  my $command = undef;
  my $init    = undef;

  if (not $self->{_can_use_docker}){
    return ( $command, $init );
  }

  if ( defined $self->{_config} ) {
    if ( defined $self->{_section} ) {
      if ( $self->{_config}{ $self->{_section} }{"no_docker"} ) {
        return ( $command, $init );
      }

      if (defined $self->{_config}{ $self->{_section} }{docker_prefix}) {
        $self->{_docker_prefix} = $self->{_config}{ $self->{_section} }{docker_prefix};
      }
    }
  }

  my $commandKey = $self->{_docker_prefix} . "docker_command";
  my $initKey    = $self->{_docker_prefix} . "docker_init";

  $command = $self->do_get_docker_value($commandKey);
  if ( defined $command ) {
    $init = $self->do_get_docker_value($initKey);
    return ( $command, $init );
  }

  my $baseCommandKey = "docker_command";
  my $baseInitKey    = "docker_init";

  $command = $self->do_get_docker_value($baseCommandKey);
  if ( defined $command ) {
    $init = $self->do_get_docker_value($baseInitKey);
    return ( $command, $init );
  }

  if ( defined $required and $required ) {
    die "Define $commandKey for task " . $self->{_name};
  }

  return ( undef, undef );
}

sub localize_files {
  my ($self, $pbs, $sample_files, $localized_files, $other_exts) = @_;
  if($self->{_use_tmp_folder} || $self->{_localize_to_local_folder} ){
    my $target_dir = $self->{_use_tmp_folder} ? '$res_dir/':'';
    if($sample_files->[0] =~ /.bam$/){
      if(defined $other_exts){
        push(@$other_exts, ".bai");
      }else{
        $other_exts = [".bai"];
      }
      my @unique_words = uniq @$other_exts;
      $other_exts = \@unique_words;
    }

    my $result = [];
    print $pbs "
echo localize start at `date`
";
    for my $old_file (@$sample_files){
      my $new_file = basename($old_file);

      push(@$result, $new_file);
      push(@$localized_files, $new_file);

      my $all_local_file = join(" ", @$localized_files);

      print $pbs "
echo $old_file      
if [[ ! -s $old_file ]]; then
  echo file not exists: $old_file
  touch ${target_dir}${new_file}.not.exist
  rm $all_local_file
  exit 1
fi

for i in {1..5}; do 
  echo iteration \$i ...
  cp -fL $old_file $new_file
  diff $old_file $new_file
  status=\$?
  if [[ \$status -eq 0 ]]; then
    break
  fi
done

if [[ \$status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch ${target_dir}${new_file}.copy.failed
  rm $all_local_file
  exit 1
fi
";

      if(defined $other_exts){
        for my $other_ext (@$other_exts){
          my $old_ext_file = $old_file . $other_ext;
          my $new_ext_file = $new_file . $other_ext;
          my $before_local_file = join(" ", @$localized_files);
          push(@$localized_files, $new_ext_file);
          my $after_local_file = join(" ", @$localized_files);
          print $pbs "
if [[ ! -s $old_ext_file ]]; then
  echo file not exists: $old_ext_file
  touch ${target_dir}${new_ext_file}.not.exist
  rm $before_local_file
  exit 1
fi

for i in {1..5}; do 
  echo iteration \$i ...
  cp -fL $old_ext_file $new_ext_file
  diff $old_ext_file $new_ext_file
  status=\$?
  if [[ \$status -eq 0 ]]; then
    break
  fi
done

if [[ \$status -ne 0 ]]; then
  echo file copied is not idential to original file, quit.
  touch ${target_dir}${new_ext_file}.copy.failed
  rm $all_local_file
  exit 1
fi
";
        }
      }
    }
      print $pbs "
ls *
echo localize end at `date`

";
    return($result);
  }else{
    return($sample_files);
  }
}

sub localize_files_in_tmp_folder {
  return( localize_files(@_)); 
}

sub clean_temp_files {
  my ($self, $pbs, $temp_files) = @_;
  if(scalar(@$temp_files) > 0){
    my $rmstr = join(" ",  @$temp_files);
    print $pbs "
rm $rmstr
";
  }
}

sub open_pbs {
  my ( $self, $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command, $final_file_can_empty, $input_file, $sh_command ) = @_;

  if ( !defined $init_command ) {
    $init_command = "";
  }

  if(!defined $sh_command){
    $sh_command = "bash";
  }

  my $module_name = $self->{_name};

  open( my $pbs, ">$pbs_file" ) or die $!;

  print $pbs "${pbs_desc}${log_desc}

$path_file
$init_command

set -o pipefail

cd '$result_dir'

";
  if ( defined $final_file ) {
    my $checkFile = $final_file_can_empty ? "-e" : "-s";
    if (is_array($final_file)){
      my @final_files = @$final_file;
      if(scalar(@final_files) > 1){
        my $final_files_1 = $final_files[0];
        my $final_files_2 = $final_files[1];

        my $delete_file = ( $final_files_1 =~ /^\// ) ? $final_files_1 : "${result_dir}/${final_files_1}";
        print $pbs "
if [[ !(1 -eq \$1) ]]; then
  if [[ ( $checkFile $final_files_1 && $checkFile $final_files_2 ) || ( -d $final_files_1 && -d $final_files_2 ) ]]; then
    echo job has already been done. if you want to do again, delete $delete_file and submit job again.
    exit 0
  fi
fi
";
      }else{
        $final_file = $final_files[0];
      }
    }
    
    if(!is_array($final_file)){
      my $delete_file = ( $final_file =~ /^\// ) ? $final_file : "${result_dir}/${final_file}";

      print $pbs "
if [[ !(1 -eq \$1) ]]; then
  if [[ ( $checkFile $final_file ) || ( -d $final_file ) ]]; then
    echo job has already been done. if you want to do again, delete $delete_file and submit job again.
    exit 0
  fi
fi
";
    }
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

  my ( $docker_command, $docker_init ) = $self->get_docker_value();
  my $is_sequenceTask = ( $module_name =~ /SequenceTask/ );
  if (  ( defined $docker_command )
    and ( ( not $is_sequenceTask ) or ( $pbs_file =~ /report/ ) ) )
  {
    if ( not defined $docker_init ) {
      $docker_init = "";
    }

    my $sing = "singularity exec -e";
    if (substr($docker_command, 0, length($sing)) eq $sing) {
      my $additional = "";
      if($docker_command !~ / -H /) {
        $additional = " -H $result_dir ";
      }
      if ($docker_command !~ / -B /) {
        $additional = $additional . " -B `pwd` -B /home ";
      }
      $docker_command = $sing . $additional . substr($docker_command, length($sing));
    }

    my $sh_file = $pbs_file . ".sh";
    my $sh_base_file = basename($sh_file);

    print $pbs "
export R_LIBS=
export PYTHONPATH=
export JAVA_HOME=
 
$docker_init

$docker_command $sh_command $sh_file 

exitcode=\$?

echo ${module_name}_end=`date`

exit \$exitcode

";
    close $pbs;
    open( $pbs, ">$sh_file" ) or die $!;

    if(!$self->init_tmp_folder($pbs, $result_dir)){
      print $pbs "
cd '$result_dir'
";
    }
    print $pbs "
set -o pipefail

$docker_init
";

  }else{
    $self->init_tmp_folder($pbs, $result_dir);
  }

  return $pbs;
}

sub close_pbs {
  my ( $self, $pbs, $pbs_file ) = @_;

  $self->clean_tmp_folder($pbs);

  my $module_name = $self->{_name};

  my ( $docker_command, $docker_init ) = $self->get_docker_value();

  if ( not defined $docker_command ) {
    #print "not defined docker_command";
    print $pbs "

exitcode=\$?

echo ${module_name}_end=`date`

exit \$exitcode
 
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

sub init_parameter {
  my ( $self, $config, $section, $create_directory ) = @_;

  $self->{_docker_prefix} = get_option( $config, $section, "docker_prefix", $self->{_docker_prefix} );
  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );
  $self->{_localize_to_local_folder} = $config->{general}{localize_to_local_folder} || get_option( $config, $section, "localize_to_local_folder", $self->{_localize_to_local_folder});
  if($self->{_localize_to_local_folder} ){
    $self->{_use_tmp_folder} = 0;  
  }elsif($self->{_forbid_tmp_folder}){
    $self->{_use_tmp_folder} = 0;
  }elsif(defined $config->{$section}{"use_tmp_folder"}){
    $self->{_use_tmp_folder} = $config->{$section}{"use_tmp_folder"};
  }elsif(defined $config->{general}{use_tmp_folder}){
    $self->{_use_tmp_folder} = $config->{general}{"use_tmp_folder"};  
  }elsif(should_use_tmp_folder($config->{$section}{target_dir})){
    $self->{_use_tmp_folder} = get_option( $config, $section, "use_tmp_folder", $self->{_use_tmp_folder} );
  }else{
    $self->{_use_tmp_folder} = 0;  
  }

  #print("target_dir=" . $config->{$section}{target_dir} . "\n");
  #print("_use_tmp_folder=" . $self->{_use_tmp_folder} . "\n");

  if ($self->{_task_suffix} ne ""){
    $self->{_suffix} = "";
  }

  return (get_parameter( $config, $section, $create_directory ));
}

1;
