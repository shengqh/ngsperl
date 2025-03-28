#!/usr/bin/perl
package CQS::ConfigUtils;

use strict;
use warnings;
use POSIX;
use Carp qw<longmess>;
use File::Basename;
use File::Copy;
use CQS::FileUtils;
use CQS::PBS;
use CQS::ClassFactory;
use CQS::StringUtils;
use CQS::CQSDebug;
use Data::Dumper;
use Hash::Merge qw( merge );
use Tie::IxHash;
use String::Util qw(trim);
use List::MoreUtils qw(uniq);
use Text::CSV qw( csv );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(
      get_ignore_sample_map
      get_all_sample_names
      print_trace
      is_string
      is_array
      is_not_array
      is_hash
      is_not_hash
      copy_files
      getValueDefaultUndef
      getValue
      get_config_section
      has_config_section
      has_option
      get_option
      get_option_include_general
      get_option_file
      get_java
      get_task_name
      get_cluster
      get_parameter
      get_param_file
      get_directory
      get_result_file
      get_first_result_file
      parse_param_file
      has_raw_files
      get_raw_files_and_keys
      do_get_unsorted_raw_files
      get_raw_files
      get_raw_files_keys
      get_raw_files_attributes
      get_raw_files2
      get_run_command
      get_option_value
      get_pair_groups
      get_pair_groups_names
      get_group_sample_map
      get_group_samplefile_map
      do_get_group_samplefile_map
      get_group_samplefile_map_key
      get_grouped_raw_files
      save_parameter_sample_file
      saveConfig 
      writeFileList
      initDefaultValue
      get_pure_pairs
      writeParameterSampleFile
      get_raw_file_list
      fix_task_name
      get_parameter_file
      get_parameter_sample_files
      is_paired_end
      get_is_paired_end_option
      get_ref_section_pbs
      get_rawfiles_option
      get_parameter_options
      get_pair_group_sample_map
      get_version_files
      get_output_ext
      get_output_ext_list
      get_program_param
      get_program
      get_interation_sample_subsample_map
      get_interation_subsample_sample_map
      get_groups
      get_unique_groups
      get_correlation_groups_by_pattern
      get_covariances
      getGroupPickResult
      getMemoryPerThread
      option_contains_arg
      left_pad
      get_key_name
      get_interval_file_map
      read_table
      merge_hash_left_precedent
      merge_hash_right_precedent
      get_groups_by_pattern_dic      
      get_groups_by_pattern
      get_groups_by_file
      get_covariances_by_pattern
      create_covariance_file_by_pattern
      create_covariance_file_by_file_pattern
      write_HTO_sample_file
      get_parameter_file_option
      get_hash_level2
      get_expanded_genes
      parse_curated_genes
      get_joined_files
      get_joined_names
      process_parameter_sample_file
      get_task_dep_pbs_map
      add_bind
      get_final_file_by_task_name
      get_groups_from_covariance_file
      )
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub is_string {
  my $data = shift;
  return(ref $data eq ref "");
}

sub is_array {
  my $data = shift;
  return(ref $data eq ref []);
}

sub is_not_array {
  return not is_array(@_);
}

sub is_hash {
  my $data = shift;
  return(ref $data eq ref {});
}

sub is_not_hash {
  return not is_hash(@_);
}

sub getValueDefaultUndef {
  my ( $def, $name ) = @_;
  if ( defined $def->{$name} ) {
    return $def->{$name};
  }

  if(is_array($name)){
    for my $n (@$name){
      #print("check $n\n");
      if ( defined $def->{$n} ) {
        #print("find " . $def->{$n}. "\n");
        return $def->{$n};
      }
    }
  }
  return undef;
}

sub getValue {
  my ( $def, $name, $defaultValue ) = @_;
  if ( defined $def->{$name} ) {
    return $def->{$name};
  }

  if(is_array($name)){
    for my $n (@$name){
      #print("check $n\n");
      if ( defined $def->{$n} ) {
        #print("find " . $def->{$n}. "\n");
        return $def->{$n};
      }
    }
  }

  if ( defined $defaultValue ) {
    return $defaultValue;
  }

  print Dumper(longmess());
  die "Define $name in user definition first.";
}

####
# Get section definition from config definition. section can be multiple layers, such as "task1::task2" which indicates $config->{task1}->{task2}
####
sub get_config_section {
  my ( $config, $section ) = @_;

  if(!defined $section){
    print Dumper(longmess());
    die "section not defined";
  }
  my @sections = split( '::', $section );
  if ( is_hash($section) ){
    my $mess = longmess();
    print Dumper( $mess );    
  }
  #print(Dumper($section));
  my $result   = $config;
  for my $curSection (@sections) {
    $result = $result->{$curSection};
    if ( !defined $result ) {
      print Dumper(longmess());
      die "Cannot find section $section!";
    }
  }
  return ($result);
}

sub print_trace {
  print STDERR "Stack Trace:\n";
  my $i = 1;
  while ( (my @call_details = (caller($i++))) ){
    print STDERR $call_details[1].":".$call_details[2]." in function ".$call_details[3]."\n";
  }    
}

sub has_config_section {
  my ( $config, $section ) = @_;

  if(!defined $section){
    print Dumper(longmess());
    die "section not defined";
  }

  my @sections = split( '::', $section );
  my $result   = $config;
  for my $curSection (@sections) {
    $result = $result->{$curSection};
    return 0 if ( !defined $result );
  }
  return 1;
}

sub has_option {
  my ( $config, $section, $key ) = @_;
  if (!defined $section){
    print Dumper(longmess());
    die "section not defined";
  }

  my $curSection = get_config_section( $config, $section );
  my $result     = $curSection->{$key};
  return ( defined $result );
}

sub get_option {
  my ( $config, $section, $key, $default ) = @_;

  my $curSection = get_config_section( $config, $section );

  my $result = $curSection->{$key};
  if ( !defined $result ) {
    if ( !defined $default ) {
      print Dumper(longmess());
      die "Define ${section}::${key} first!";
    }
    else {
      $result = $default;
    }
  }

  return ($result);
}

#get option from task section and general section
sub get_option_include_general {
  my ( $config, $section, $key, $default ) = @_;

  my $curSection = get_config_section( $config, $section );

  my $result = $curSection->{$key};
  if ( !defined $result ) {
    $result = $config->{general}{$key};
    if (!defined $result) {
      if ( !defined $default ) {
        print Dumper(longmess());
        die "Define ${section}::${key} first!";
      }
      else {
        $result = $default;
      }
    }
  }

  return ($result);
}

sub get_option_file {
  my ( $config, $section, $key, $default ) = @_;

  my $curSection = get_config_section( $config, $section );

  my $result = $curSection->{$key};
  if ( !defined $result ) {
    if (defined $config->{general}){
      $result = $config->{general}{$key};
    }
  }

  if ( !defined $result ) {
    print Dumper(longmess());
    die "Define ${section}::${key} first!";
  }

  if ( !is_debug() && !-e $result ) {
    print Dumper(longmess());
    die "$key $result defined but not exists!";
  }
  return ($result);
}

sub get_cluster {
  my ( $config, $section ) = @_;

  my $curSection = get_config_section( $config, $section );

  my $cluster_name;
  if ( defined $curSection->{cluster} ) {
    $cluster_name = get_option_value( $curSection->{cluster}, "slurm" );
  }
  else {
    $cluster_name = get_option_value( $config->{general}{cluster}, "slurm" );
  }

  my $cluster;
  if ( $cluster_name eq "torque" ) {
    $cluster = instantiate("CQS::ClusterTorque");
  }
  elsif ( $cluster_name eq "torqueSingapore" ) {
    $cluster = instantiate("CQS::ClusterTorqueSingapore");
  }
  else {
    $cluster = instantiate("CQS::ClusterSLURM");
  }

  return ($cluster);
}

sub get_value_in_section_or_general {
  my ( $config, $section, $name, $defaultvalue ) = @_;

  my $curSection = get_config_section( $config, $section );

  my $result;
  if ( defined $curSection->{$name} ) {
    $result = get_option_value( $curSection->{$name}, $defaultvalue );
  }
  else {
    $result = get_option_value( $config->{general}{$name}, $defaultvalue );
  }

  return ($result);
}

sub get_java {
  my ( $config, $section ) = @_;
  return ( get_value_in_section_or_general( $config, $section, "java", "java" ) );
}

sub get_task_name {
  my ( $config, $section, $create_directory ) = @_;
  my $task_name = get_option( $config, $section, "task_name", "" );
  if ( $task_name eq "" ) {
    $task_name = get_option( $config, "general", "task_name" );
  }
  $task_name =~ s/\s/_/g;
  return ($task_name);
}

sub get_parameter {
  my ( $config, $section, $create_directory ) = @_;

  my $curSection = get_config_section( $config, $section );

  $create_directory = 1 if !defined($create_directory);

  my $task_name = get_task_name( $config, $section );

  my $cluster = get_cluster(@_);

  my $path_file = get_param_file( $curSection->{path_file}, "path_file", 0 );
  if ( !defined $path_file ) {
    $path_file = get_param_file( $config->{general}{path_file}, "path_file", 0 );
  }
  if ( defined $path_file && -e $path_file ) {
    $path_file = "source $path_file";
  }
  else {
    $path_file = "";
  }

  my $refPbs     = get_option( $config, $section, "pbs" );
  my $target_dir = get_option( $config, $section, "target_dir" );
  $target_dir =~ s|//|/|g;
  $target_dir =~ s|/$||g;
  my ( $log_dir, $pbs_dir, $result_dir ) = init_dir( $target_dir, $create_directory, $config, $section );

  my $pbs_desc = $cluster->get_cluster_desc( $refPbs, $config );

  my $option    = get_option( $config, $section, "option",    "" );
  my $sh_direct = get_option( $config, $section, "sh_direct", 0 );

  my $init_command = get_option( $config, $section, "init_command", "" );

  if ($sh_direct) {
    $sh_direct = "sh";
  }
  else {
    $sh_direct = $cluster->get_submit_command();
  }

  my $thread = $cluster->get_cluster_thread($refPbs);
  my $memory = $cluster->get_cluster_memory($refPbs);
  
  if ($option =~ /__MEMORY__/) {
    $option =~ s/__MEMORY__/$memory/g;
  }
  
  if ($option =~ /__THREAD__/) {
    $option =~ s/__THREAD__/$thread/g;
  }

  return ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command );
}

#get parameter which indicates a file. If required, not defined or not exists, die. If defined but not exists, die.
#returned file either undefined or exists.
sub get_param_file {
  my ( $file, $name, $required, $checkExists ) = @_;

  if ( not defined $checkExists ) {
    $checkExists = 1;
  }

  my $result = $file;

  if ($required) {
    if ( !defined $file ) {
      print Dumper(longmess());
      die "$name was not defined!";
    }

    if ( !is_debug() and ( $checkExists and !-e $file ) ) {
      print Dumper(longmess());
      die "$name $file defined but not exists!";
    }
  }
  else {
    if ( defined($file) ) {
      if ( $file eq "" ) {
        undef($result);
      }
      elsif ( !is_debug() and ( $checkExists and !-e $file ) ) {
        print Dumper(longmess());
        die "$name $file defined but not exists!";
      }
    }
  }
  return ($result);
}

#get parameter which indicates a file. If required, not defined or not exists, die. If defined but not exists, die.
#returned file either undefined or exists.
sub get_directory {
  my ( $config, $section, $name, $required ) = @_;

  my $curSection = get_config_section( $config, $section );

  if (!defined $name){
    print Dumper(longmess());
    die "parameter name must be defined!";
  }

  my $result = $curSection->{$name};

  if ($required) {
    if ( !defined $result ) {
      print Dumper(longmess());
      die "$name was not defined!";
    }

    if ( !is_debug() && !-e $result ) {
      print Dumper(longmess());
      die "$name $result defined but not exists!";
    }
  }
  else {
    if ( defined($result) ) {
      if ( $result eq "" ) {
        undef($result);
      }
      elsif ( !is_debug() && !-e $result ) {
        print Dumper(longmess());
        die "$name $result defined but not exists!";
      }
    }
  }
  return ($result);
}

sub get_result_file {
  my ( $config, $refSectionName, $pattern ) = @_;
  my $refSection = get_config_section( $config, $refSectionName );
  if ( defined $refSection->{class} ) {
    my $myclass = instantiate( $refSection->{class} );
    my $result  = $myclass->result( $config, $refSectionName, $pattern, 1 );
    return ($result);
  }

  return (undef);
}

sub get_first_result_file {
  my ( $config, $refSectionName, $pattern ) = @_;
  my $refSection = get_config_section( $config, $refSectionName );

  if (is_string($refSection)) {
    if ((defined $pattern) and ($pattern ne "")){
      if ($refSection =~ /$pattern/){
        return($refSection);
      }else{
        print Dumper(longmess());
        die "$refSection doesn't match pattern $pattern";
      }
    }else{
      return($refSection);
    }
  }

  if (is_array($refSection)) {
    if ( scalar(@$refSection) > 0 ) {
      if ((defined $pattern) and ($pattern ne "")){
        for my $file (@$refSection){
          if ($file =~ /$pattern/){
            return($file);
          } 
        }
        print Dumper(longmess());
        die "no file in $refSectionName matches pattern $pattern";
      }else{
        return($refSection->[0]);
      }
    }else{
      print Dumper(longmess());
      die "no file in $refSectionName";
    }
  }

  if (is_hash($refSection)) {
    my $result;
    if ( defined $refSection->{class} ) {
      my $myclass = instantiate( $refSection->{class} );
      $result  = $myclass->result( $config, $refSectionName, $pattern, 1 );
    }else{
      $result = $refSection;
    }
    foreach my $k ( sort keys %{$result} ) {
      my @files = @{ $result->{$k} };
      if ( scalar(@files) > 0 ) {
        return $files[0];
      }
    }
    print Dumper(longmess());
    die "section $refSectionName return nothing for pattern $pattern !";
  }

  print Dumper(longmess());
  die "section $refSectionName return nothing for pattern $pattern !";
}

sub parse_param_file {
  my ( $config, $section, $key, $required ) = @_;

  my $curSection = get_config_section( $config, $section );

  if(!defined $key){
    print Dumper(longmess());
    die "parameter key must be defined!";
  }

  if ( defined $curSection->{$key} ) {
    return $curSection->{$key};
  }

  my $key_ref = $key . "_ref";
  if ( defined $curSection->{$key_ref} ) {
    my $refSectionName = $curSection->{$key_ref};
    my $pattern;
    if ( is_array($refSectionName) ) {
      my @parts = @{$refSectionName};
      if ( scalar(@parts) == 2 ) {
        $pattern        = $parts[1];
        $refSectionName = $parts[0];
      }
      else {
        $refSectionName = $parts[0];
      }
    }

    my $result = get_first_result_file( $config, $refSectionName, $pattern );

    if ( defined $result ) {
      return ($result);
    }
  }

  if ($required) {
    print Dumper(longmess());
    die "define ${section}::${key} first.";
  }

  return undef;
}

sub has_raw_files {
  my ( $config, $section, $mapname ) = @_;

  my $curSection = get_config_section( $config, $section );

  if ( !defined $mapname ) {
    $mapname = "source";
  }

  my $mapname_ref        = $mapname . "_ref";
  my $mapname_config_ref = $mapname . "_config_ref";

  return ( defined $curSection->{$mapname} ) || ( defined $curSection->{$mapname_ref} ) || ( defined $curSection->{$mapname_config_ref} );
}

sub get_refmap {
  my ( $config, $section, $mapname, $pattern ) = @_;

  my $curSection         = get_config_section( $config, $section );
  my $mapname_ref        = $mapname . "_ref";
  my $mapname_config_ref = $mapname . "_config_ref";

  my $result = {};
  if ( defined $curSection->{$mapname_ref} ) {

    #in same config
    my $targetSection = $curSection->{$mapname_ref};

    if ( is_hash($targetSection) ) {
      return ( $result, 0 );
    }

    if ( is_array($targetSection) ) {
      my @parts      = @{$targetSection};
      my $partlength = scalar(@parts);
      for ( my $index = 0 ; $index < $partlength ; ) {
        if ( !has_config_section( $config, $parts[$index] ) ) {
          print Dumper(longmess());
          die "undefined section " . Dumper($parts[$index]) . " in $section for $mapname ";
        }
        get_config_section( $config, $parts[$index] );

        if ( $index == ( $partlength - 1 ) || has_config_section( $config, $parts[ $index + 1 ] ) ) {
          $result->{$index} = { config => $config, section => $parts[$index], pattern => $pattern };
          $index++;
        }
        else {
          $result->{$index} = { config => $config, section => $parts[$index], pattern => $parts[ $index + 1 ] };
          $index += 2;
        }
      }
    }
    else {
      if ( !has_config_section( $config, $targetSection ) ) {
        print Dumper(longmess());
        die "undefined section $targetSection";
      }
      $result->{1} = { config => $config, section => $targetSection, pattern => $pattern };
    }
  }
  else {

    #in another config, has to be array
    my $refSectionName = $curSection->{$mapname_config_ref};
    if ( is_not_array($refSectionName) ) {
      print Dumper(longmess());
      die "$mapname_config_ref has to be defined as ARRAY with [config, section, pattern]";
    }
    my @parts      = @{$refSectionName};
    my $partlength = scalar(@parts);
    for ( my $index = 0 ; $index < $partlength - 1 ; ) {
      my $targetConfig  = $parts[$index];
      my $targetSection = $parts[ $index + 1 ];

      if ( is_not_hash($targetConfig) ) {
        print Dumper(longmess());
        die
"$mapname_config_ref has to be defined as ARRAY with [config1, section1, pattern1,config2, section2, pattern2] or [config1, section1,config2, section2] format. config should be hash and section should be string";
      }

      if ( !has_config_section( $targetConfig, $targetSection ) ) {
        print Dumper(longmess());
        die "undefined section $targetSection in $mapname_config_ref of $section";
      }

      if ( $index == ( $partlength - 2 ) || is_hash( $parts[ $index + 2 ] ) ) {
        $result->{$index} = { config => $targetConfig, section => $targetSection, pattern => $pattern };
        $index += 2;
      }
      else {
        $result->{$index} = { config => $targetConfig, section => $targetSection, pattern => $parts[ $index + 2 ] };
        $index += 3;
      }
    }
  }
  return ( $result, 0 );
}

sub get_raw_file_list {
  my ( $config, $section, $mapname, $mustBeOne ) = @_;

  my $curSection = get_config_section( $config, $section );

  if ( !defined $mustBeOne ) {
    $mustBeOne = 0;
  }

  if ( !defined $mapname ) {
    $mapname = "source";
  }
  my $mapname_ref        = $mapname . "_ref";
  my $mapname_config_ref = $mapname . "_config_ref";

  if ( defined $curSection->{$mapname} ) {
    return $curSection->{$mapname};
  }

  if ( defined $curSection->{$mapname_ref} || defined $curSection->{$mapname_config_ref} ) {
    my ( $refmap, $returnNow ) = get_refmap( $config, $section, $mapname, undef );
    if ($returnNow) {
      return $refmap;
    }

    my $result     = [];
    my @sortedKeys = sort { $a <=> $b } keys %$refmap;
    for my $index (@sortedKeys) {
      my $values       = $refmap->{$index};
      my $targetConfig = $values->{config};
      my $section      = $values->{section};
      if(not defined $section){
        print Dumper(longmess());
        die("ERROR: no section defined for index " . $index);
      }
      my $pattern      = $values->{pattern};

      my $targetSection = get_config_section( $targetConfig, $section );

      my %myres = ();
      if ( defined $targetSection->{class} ) {
        my $myclass = instantiate( $targetSection->{class} );
        %myres = %{ $myclass->result( $targetConfig, $section, $pattern, 1 ) };
      }
      else {
        my ( $res, $issource ) = do_get_unsorted_raw_files( $targetConfig, $section, 1, undef, $pattern, 1 );
        %myres = %{$res};
      }

      my $bFound    = 0;
      my @curResult = ();
      for my $myvalues ( values %myres ) {
        if(is_not_array($myvalues)){
          print Dumper(longmess());
          die "Return value should be array.";
        }
        
        if ( scalar(@$myvalues) > 0 ) {
          push( @curResult, @$myvalues );
          $bFound = 1;
        }
      }
      if ( not $bFound ) {
        print Dumper(longmess());
        die "Cannot find file for " . $values->{section} . " and pattern " . $values->{pattern};
      }

      if ( $mustBeOne && scalar(@curResult) > 1 ) {
        print Dumper(longmess());
        die "Only one result allowed but multiple found " . $values->{section} . " and pattern " . $values->{pattern} . "\n" . join( "\n  ", @curResult );
      }

      push( @$result, @curResult );
    }

    return $result;
  }

  print Dumper(longmess());
  die "define $mapname or $mapname_ref or $mapname_config_ref for $section";
}

sub get_ignore_sample_map {
  my ($config, $section) = @_;

  my $result = {};
  if($config->{ignore_samples}){
    if(!get_option($config, $section, "use_all_samples", 0)){
      my $ignore_samples = $config->{ignore_samples};
      if(defined $ignore_samples){
        #print(Dumper($ignore_samples));
        if(is_array($ignore_samples)){
          my %ignore_map = map { $_ => 1 } @$ignore_samples;
          $result = \%ignore_map;
        }
      }
    }
  }
  
  return($result)
}

sub do_get_unsorted_raw_files_no_ignored {
  my ( $config, $section, $returnself, $mapname, $pattern, $removeEmpty ) = @_;

  if(!defined $removeEmpty){
    if(!defined $mapname || $mapname eq "source" || $mapname eq "parameterSampleFile1"){
      $removeEmpty = get_option($config, $section, "remove_empty_source", 0);
    }
  }

  my $curSection = get_config_section( $config, $section );

  if ( !defined $mapname ) {
    $mapname = "source";
  }
  my $mapname_ref        = $mapname . "_ref";
  my $mapname_config_ref = $mapname . "_config_ref";

  if ( defined $curSection->{$mapname} ) {
    my $result = $curSection->{$mapname};
    if (is_hash($result)) {
      return ( $curSection->{$mapname}, 1 );
    }else{
      warn("Did you forget to use _ref for $mapname of $section?");
      return ( $curSection->{$mapname}, 1 );
    }
  }

  if ( defined $curSection->{$mapname_ref} || defined $curSection->{$mapname_config_ref} ) {
    my ( $refmap, $returnNow ) = get_refmap( $config, $section, $mapname, $pattern );
    if ($returnNow) {
      return ( $refmap, 1 );
    }

    my %result     = ();
    my @sortedKeys = sort { $a <=> $b } keys %$refmap;
    for my $index (@sortedKeys) {
      my $values       = $refmap->{$index};
      my $targetConfig = $values->{config};
      my $section      = $values->{section};
      my $pattern      = $values->{pattern};

      my $targetSection = get_config_section( $targetConfig, $section );

      my %myres = ();
      if ( defined $targetSection->{class} ) {
        my $myclass = instantiate( $targetSection->{class} );
        %myres = %{ $myclass->result( $targetConfig, $section, $pattern, $removeEmpty ) };
      }
      else {
        my ( $res, $issource ) = do_get_unsorted_raw_files_no_ignored( $targetConfig, $section, 1, undef, $pattern, $removeEmpty );
        %myres = %{$res};
      }

      my $refcount = keys %myres;
      for my $mykey ( keys %myres ) {
        my $myvalues = $myres{$mykey};
        if ( is_string($myvalues) ) {
          $myvalues = [$myvalues];
        }

        if ( is_array($myvalues) && ( scalar( @{$myvalues} ) > 0 ) ) {
          if ( exists $result{$mykey} ) {
            my $oldvalues = $result{$mykey};
            if ( is_array($oldvalues) ) {
              my @merged = ( @{$oldvalues}, @{$myvalues} );

              #print "merged ARRAY ", Dumper(\@merged);
              $result{$mykey} = \@merged;
            }
            else {
              print Dumper(longmess());
              die "The source of $section->$mapname should be all HASH or all ARRAY";
            }
          }
          else {
            $result{$mykey} = $myvalues;
          }
        }

        if ( is_hash($myvalues) && ( scalar( keys %{$myvalues} ) > 0 ) ) {
          if ( exists $result{$mykey} ) {
            my $oldvalues = $result{$mykey};
            if ( is_hash($oldvalues) ) {
              $result{$mykey} = merge_hash_right_precedent( $oldvalues, $myvalues );
            }
            else {
              print Dumper(longmess());
              die "The source of $section->$mapname should be all HASH or all ARRAY";
            }
          }
          else {
            $result{$mykey} = $myvalues;
          }
        }
      }

      #print "--------------- $section, $mapname, $index ----------------\n";
      #print Dumper(%result);
    }

    my $final = \%result;
    return ( $final, 0 );
  }

  if ($returnself) {
    if ( defined $pattern ) {
      my $result = {};
      for my $key ( sort keys %$curSection ) {
        my $values = $curSection->{$key};
        $result->{$key} = filter_array( $values, $pattern );
      }
      return ( $result, 0 );
    }
    else {
      return ( $curSection, 0 );
    }
  }
  else {
    print Dumper(longmess());
    die "define $mapname or $mapname_ref or $mapname_config_ref for $section";
  }
}

sub do_get_unsorted_raw_files {
  my ( $config, $section, $returnself, $mapname, $pattern, $removeEmpty ) = @_;

  my ($res, $is_source) = do_get_unsorted_raw_files_no_ignored($config, $section, $returnself, $mapname, $pattern, $removeEmpty);

  if (is_hash($res)){
    my $ignored = get_ignore_sample_map($config, $section);
    if (%$ignored){
      my $newres = {};
      for my $key (keys %$res){
        if(!exists($ignored->{$key})){
          $newres->{$key} = $res->{$key};
        }
      }
      return($newres, $is_source);
    }
  }

  return($res, $is_source);
}

sub get_ref_section_pbs {
  my ( $config, $section, $mapname ) = @_;

  my $result = {};

  my $curSection = get_config_section( $config, $section );

  my $mapname_ref        = $mapname . "_ref";
  my $mapname_config_ref = $mapname . "_config_ref";

  if ( defined $curSection->{$mapname_ref} || defined $curSection->{$mapname_config_ref} ) {
    my ( $refmap, $unused ) = get_refmap( $config, $section, $mapname );

    my @sortedKeys = sort { $a <=> $b } keys %$refmap;
    for my $index (@sortedKeys) {
      my $values       = $refmap->{$index};
      my $targetConfig = $values->{config};
      my $section      = $values->{section};

      my $targetSection = get_config_section( $targetConfig, $section );
      if ( is_not_hash($targetSection) ) {
        next;
      }

      if ( defined $targetSection->{class} ) {
        my $myclass = instantiate( $targetSection->{class} );

        #print ($targetSection->{class} . " " . $section. "\n");
        my $result_pbs_map = $myclass->get_result_pbs( $config, $section );
        for my $resultKey ( keys %$result_pbs_map ) {
          my $newpbs  = $result_pbs_map->{$resultKey};
          my $pbslist = $result->{$resultKey};
          if ( defined $pbslist ) {
            push @$pbslist, $newpbs;
          }
          else {
            $pbslist = [$newpbs];
          }
          $result->{$resultKey} = $pbslist;
        }
      }
      else {
        next;
      }
    }
  }

  return ($result);
}

sub do_get_raw_files_keys {
  my $resultUnsorted = shift;
  my @keys           = grep { $_ !~ /^\./ } keys %$resultUnsorted;
  my @result;
  if ( exists $resultUnsorted->{".order"} ) {
    my $orders = $resultUnsorted->{".order"};
    @result = @$orders;
    print Dumper(longmess());
    die "number of key defined in .order not equals to actual keys for @_" if ( scalar(@result) != scalar(@keys) );
  }
  else {
    @result = sort @keys;
  }
  return \@result;
}

sub get_raw_files_keys {
  my ( $resultUnsorted, $issource ) = do_get_unsorted_raw_files(@_);
  return do_get_raw_files_keys($resultUnsorted);
}

sub get_sorted_raw_files {
  my $resultUnsorted = shift;
  if (!is_hash($resultUnsorted)){
    print("resultUnsorted = " . Dumper($resultUnsorted));
    print_trace();
  }
  my @keys           = grep { $_ !~ /^\./ } keys %$resultUnsorted;
  my %result;
  tie %result, 'Tie::IxHash';
  my @orderedKeys;
  if ( exists $resultUnsorted->{".order"} ) {
    my $orders = $resultUnsorted->{".order"};
    @orderedKeys = @$orders;
    if ( scalar(@orderedKeys) != scalar(@keys) ){
      print Dumper(longmess());
      die "number of key defined in .order not equals to actual keys for @_" ;
    }
  }
  else {
    @orderedKeys = sort @keys;
  }
  for my $key (@orderedKeys) {
    $result{$key} = $resultUnsorted->{$key};
  }
  return \%result;
}

sub do_get_raw_files {
  my ( $resultUnsorted, $issource ) = do_get_unsorted_raw_files(@_);
  my $result = get_sorted_raw_files($resultUnsorted);
  return ( $result, $issource );
}

sub get_raw_files_attributes {
  my %result;

  my ( $resultUnsorted, $issource ) = do_get_unsorted_raw_files(@_);

  my @orderedKeys = keys %$resultUnsorted;

  for my $key (@orderedKeys) {
    if ( $key =~ /^\./ ) {
      $result{$key} = $resultUnsorted->{$key};
    }
  }
  return ( \%result );
}

sub get_raw_files_and_keys {
  my ( $resultUnsorted, $issource ) = do_get_unsorted_raw_files(@_);
  my $result      = get_sorted_raw_files($resultUnsorted);
  my $orderedKeys = do_get_raw_files_keys($resultUnsorted);
  return ( $result, $orderedKeys );
}

sub get_raw_files {
  my ( $config, $section, $mapname, $pattern, $removeEmpty ) = @_;
  my ( $result, $issource ) = do_get_raw_files( $config, $section, 0, $mapname, $pattern, $removeEmpty );
  return $result;
}

#return raw files and if the raw files are extracted from source directly
sub get_raw_files2 {
  my ( $config, $section, $mapname, $pattern ) = @_;
  return do_get_raw_files( $config, $section, 0, $mapname, $pattern );
}

sub get_run_command {
  my $sh_direct = shift;
  return ("MYCMD=\"$sh_direct\" \n");
}

sub get_run_command_old {
  my $sh_direct = shift;
  if ($sh_direct) {
    return ("MYCMD=\"bash\" \n");
  }
  else {
    return ("type -P qsub &>/dev/null && MYCMD=\"qsub\" || MYCMD=\"bash\" \n");
  }
}

sub get_option_value {
  my ( $value, $defaultValue ) = @_;
  if ( !defined $value ) {
    return ($defaultValue);
  }
  else {
    return ($value);
  }
}

sub get_pair_groups {
  my ( $pairs, $pair_name ) = @_;
  my $group_names;
  my $ispaired       = 0;
  my $tmpgroup_names = $pairs->{$pair_name};
  if ( is_hash($tmpgroup_names) ) {
    $group_names = $tmpgroup_names->{"groups"};
    $ispaired    = $tmpgroup_names->{"paired"};
  }
  else {
    $group_names = $tmpgroup_names;
  }
  if ( !defined $ispaired ) {
    $ispaired = 0;
  }
  return ( $ispaired, $group_names );
}

sub get_pair_groups_names {
  my ( $pairs, $pair_name ) = @_;
  my $group_names;
  my $pairedNames;
  my $tmpgroup_names = $pairs->{$pair_name};
  if ( is_hash($tmpgroup_names) ) {
    $group_names = $tmpgroup_names->{"groups"};
    $pairedNames = $tmpgroup_names->{"paired"};
  }
  else {
    $group_names = $tmpgroup_names;
  }
  return ( $pairedNames, $group_names );
}

##################################################
#Get pure pair-group information for UniqueR tasks
##################################################
sub get_pure_pairs {
  my ($pairs) = @_;

  if ( !defined $pairs ) {
    return (undef);
  }

  my $result = get_sorted_raw_files($pairs);

  for my $pair_name ( keys %$result ) {
    my $tmpgroup_names = $result->{$pair_name};
    if ( is_hash($tmpgroup_names) ) {
      my $group_names = $tmpgroup_names->{"groups"};
      $result->{$pair_name} = $group_names;
    }
  }

  return ($result);
}

#Return
#{
#  groupName1 => [
#    [sample_name1_1, sampleFile1_1_1, sampleFile1_1_2],
#    [sample_name1_2, sampleFile1_2_1, sampleFile1_2_2],
#  ],
#  groupName2 => [
#    [sample_name2_1, sampleFile2_1_1, sampleFile2_1_2],
#    [sample_name2_2, sampleFile2_2_1, sampleFile2_2_2],
#  ],
#}
sub get_group_sample_map {
  my ( $config, $section, $samplePattern ) = @_;

  my $raw_files = get_raw_files( $config, $section, "source", $samplePattern );
  my $groups    = get_raw_files( $config, $section, "groups" );
  my %group_sample_map = ();
  for my $group_name ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$group_name} };
    my @gfiles  = ();
    foreach my $sample_name (@samples) {
      my $sample_files = $raw_files->{$sample_name};
      if(!defined $sample_files){
        print Dumper(longmess());
        die "Cannot find file of $sample_name for group $group_name";
      }
      my @bam_files = @$sample_files;
      my @sambam    = ( $sample_name, @bam_files );
      push( @gfiles, \@sambam );
    }
    $group_sample_map{$group_name} = \@gfiles;
  }

  return \%group_sample_map;
}

#Return
#{
#  groupName1 => [sampleFile1_1_1, sampleFile1_1_2, sampleFile1_2_1, sampleFile1_2_2],
#  groupName2 => [sampleFile2_1_1, sampleFile2_1_2, sampleFile2_2_1, sampleFile2_2_2],
#}
sub get_group_samplefile_map {
  my ( $config, $section, $sample_pattern ) = @_;
  return get_group_samplefile_map_key( $config, $section, $sample_pattern, "groups" );
}

sub do_get_group_samplefile_map {
  my ($groups, $raw_files) = @_;
  my %group_sample_map = ();
  for my $group_name ( sort keys %{$groups} ) {
    my $gfiles        = [];
    $group_sample_map{$group_name} = $gfiles;
    my $group_samples = $groups->{$group_name};
    if ( ref $group_samples eq ref "" ) {
      if(! defined $raw_files->{$group_samples}){
        push( @$gfiles, $group_samples );
        continue
      }else{
        $group_samples = [$group_samples];
      }
    }

    my @samples = @{$group_samples};
    foreach my $sample_name (@samples) {
      if ( not defined $raw_files->{$sample_name} ) {
        print Dumper(longmess());
        die "Cannot find $sample_name of group $group_name, check your sample names";
      }
      my @bam_files = @{ $raw_files->{$sample_name} };
      foreach my $bam_file (@bam_files) {
        push( @$gfiles, $bam_file );
      }
    }
  }
  return \%group_sample_map;
}

sub get_group_samplefile_map_key {
  my ( $config, $section, $sample_pattern, $group_key ) = @_;

  my $raw_files = get_raw_files( $config, $section, "source", $sample_pattern );
  my $groups    = get_raw_files( $config, $section, $group_key );

  return(do_get_group_samplefile_map($groups, $raw_files));
}

sub get_pair_group_sample_map {
  my ( $pairs, $groups ) = @_;
  my $pure_pairs = get_pure_pairs($pairs);
  my $result     = {};
  for my $pair_name ( keys %$pure_pairs ) {
    my $pair_map = {};
    $result->{$pair_name} = $pair_map;
    my $groups_names = $pure_pairs->{$pair_name};
    for my $group_name (@$groups_names) {
      my $samples = $groups->{$group_name};
      $pair_map->{$group_name} = $samples;
    }
  }

  return ($result);
}

sub get_grouped_raw_files {
  my ( $config, $section, $group_key ) = @_;
  my $raw_files;
  if ( has_raw_files( $config, $section, $group_key ) ) {
    $raw_files = get_group_samplefile_map_key( $config, $section, "", $group_key );
  }
  else {
    $raw_files = get_raw_files( $config, $section );
  }
  return $raw_files;
}

sub save_parameter_sample_file {
  my ( $config, $section, $key, $outputFile ) = @_;

  my $curSection = get_config_section( $config, $section );

  if ( has_raw_files( $config, $section, $key ) ) {
    my %temp = %{ get_raw_files( $config, $section, $key ) };
    my @orderedSampleNames;
    my $sampleFileOrder = $curSection->{ $key . "Order" };
    if ( defined $sampleFileOrder ) {
      @orderedSampleNames = @{$sampleFileOrder};
    }
    else {
      @orderedSampleNames = sort keys %temp;
    }

    my $fileOnly = get_option( $config, $section, $key . "_fileonly", 0 );
    my $fileFirst = get_option( $config, $section, $key . "_fileFirst", 1 );
    my $join_delimiter = get_option( $config, $section, $key . "_join_delimiter", "" );

    open( my $list, '>', $outputFile ) or die "Cannot create $outputFile";
    foreach my $sample_name (@orderedSampleNames) {
      my $subSampleFiles = $temp{$sample_name};
      my $refstr         = ref($subSampleFiles);
      if ( $refstr eq 'HASH' ) {
        foreach my $groupName ( sort keys %$subSampleFiles ) {
          my $groupSampleNames = $subSampleFiles->{$groupName};
          for my $groupSampleName (@$groupSampleNames) {
            if ($fileOnly){
              print $list $groupSampleName . "\n";
            }else{
              if($fileFirst){
                print $list "${groupSampleName}\t${groupName}\t${sample_name}\n";
              }else{
                print $list "${sample_name}\t${groupName}\t${groupSampleName}\n";
              }
            }
          }
        }
      }
      elsif ( $refstr eq 'ARRAY' ) {
        if($join_delimiter ne ""){
          my $sample_str = join($join_delimiter, @$subSampleFiles);
          if ($fileOnly){
            print $list $sample_str . "\n";
          }else{
            if($fileFirst){
              print $list "$sample_str\t$sample_name\n";
            }else{
              print $list "$sample_name\t$sample_str\n";
            }
          }
        }else{
          foreach my $subSampleFile (@$subSampleFiles) {
            if ($fileOnly){
              print $list $subSampleFile . "\n";
            }else{
              if($fileFirst){
                print $list "$subSampleFile\t$sample_name\n";
              }else{
                print $list "$sample_name\t$subSampleFile\n";
              }
            }
          }
        }
      }
      else {
        if(!defined $subSampleFiles){
          $subSampleFiles = "";
        }
        if ($fileOnly){
          print $list $subSampleFiles . "\n";
        }else{
          if($fileFirst){
            print $list "$subSampleFiles\t$sample_name\n";
          }else{
            print $list "$sample_name\t$subSampleFiles\n";
          }
        }
      }
    }
    close($list);
    return ($outputFile);
  }
  else {
    return ("");
  }
}

sub copy_files {
  my ($files_str, $result_dir) = @_;

  my @files = split( ",|;", $files_str );

  my $rnum = scalar(@files);

  for my $index (0 .. $#files){
    my $cur_file = $files[$index];
    my $is_absolute = File::Spec->file_name_is_absolute($cur_file);
    if ( !$is_absolute ) {
      $cur_file = dirname(__FILE__) . "/$cur_file";
    }
    if ( !( -e $cur_file ) ) {
      die("file $cur_file defined but not exists!");
    }

    if ($index < $rnum){
      my $remote_file = $result_dir . "/" . basename($cur_file);
      copy($cur_file, $remote_file);
    }
  }
}

sub saveConfig {
  my ( $def, $config ) = @_;

  my $targetFile = $def->{target_dir} . "/" . basename($0);
  print( $0 . " => " . $targetFile . "\n" );
  copy( $0, $targetFile );

  my $def_file;
  if ( $def->{target_dir} =~ /\/$/ ) {
    $def_file = $def->{target_dir} . $def->{task_name} . '.def';
  }
  else {
    $def_file = $def->{target_dir} . '/' . $def->{task_name} . '.def';
  }

  open( my $sh1, ">$def_file" ) or die "Cannot create $def_file";
  print $sh1 Dumper($def);
  close $sh1;
  print "Saved user definition file to " . $def_file . "\n";

  my $config_file;
  if ( $def->{target_dir} =~ /\/$/ ) {
    $config_file = $def->{target_dir} . $def->{task_name} . '.config';
  }
  else {
    $config_file = $def->{target_dir} . '/' . $def->{task_name} . '.config';
  }

  open( my $sh2, ">$config_file" ) or die "Cannot create $config_file";
  print $sh2 Dumper($config);
  close $sh2;
  print "Saved configuration file to " . $config_file . "\n";
}

sub writeFileList {
  my ( $fileName, $fileMap, $exportAllFiles, $fileFirst, $fileOnly ) = @_;
  if (not defined $exportAllFiles) {
    $exportAllFiles = 0;
  }

  if (not defined $fileFirst) {
    $fileFirst = 0;
  }

  if (not defined $fileOnly) {
    $fileOnly = 0;
  }

  open( my $fl, ">$fileName" ) or die "Cannot create $fileName";
  for my $sample_name ( sort keys %$fileMap ) {
    my $files = $fileMap->{$sample_name};
    if( is_array($files) ){
      if ( not $exportAllFiles ) {
        $files = [$files->[0]]
      }  
    }else{
      $files = [$files]
    }

    for my $eachFile (@$files) {
      if ($fileOnly) {
        print $fl $eachFile, "\n";
      }else {
        if ($fileFirst) {
          print $fl $eachFile, "\t", $sample_name, "\n";
        }else{
          print $fl $sample_name, "\t", $eachFile, "\n";
        }
      }
    }
  }
  close($fl);
}

sub initDefaultValue {
  my ( $def, $name, $defaultValue ) = @_;
  if ( !defined $def->{$name} ) {
    if ( !defined $defaultValue ) {
      print Dumper(longmess());
      die "defaultValue cannot be undefined for $name.";
    }
    $def->{$name} = $defaultValue;
  }
  return $def;
}

sub writeParameterSampleFile {
  my ( $config, $section, $resultDir, $index, $removeEmpty, $cur_sample_name ) = @_;
  my $result      = "";
  my $task_suffix = get_option( $config, $section, "suffix", "" );
  my $key         = "parameterSampleFile" . $index;
  if ( has_raw_files( $config, $section, $key ) ) {
    my $temp = get_raw_files( $config, $section, $key, undef, $removeEmpty );
    my @orderedSampleNames;
    my $keyOrder                 = $key . "Order";
    my $parameterSampleFileOrder = $config->{$section}{$keyOrder};
    if ( defined $parameterSampleFileOrder ) {
      @orderedSampleNames = @{$parameterSampleFileOrder};
    }
    else {
      @orderedSampleNames = sort keys %$temp;
    }

    my @outputNames              = ();
    my $keyNames                 = $key . "Names";
    my $parameterSampleFileNames = $config->{$section}{$keyNames};
    if ( defined $parameterSampleFileNames ) {
      @outputNames = @$parameterSampleFileNames;
    }

    my $delimiter = get_option( $config, $section, $key . "_delimiter", "\t" );
    my $name_first = get_option( $config, $section, $key . "_name_first", 0 );
    $result = "fileList${index}${task_suffix}.txt";
    open( my $list, ">$resultDir/$result" ) or die "Cannot create $result";
    my $nameIndex = -1;
    foreach my $sample_name (@orderedSampleNames) {
      if(defined $cur_sample_name){
        if($cur_sample_name ne $sample_name){
          next;
        }
      }

      my $subSampleFiles = $temp->{$sample_name};
      if(not defined $subSampleFiles) {
        if($name_first){
          print $list "$sample_name${delimiter}\n";
        }else{
          print $list "${delimiter}$sample_name\n";
        }
        next;
      }

      my $refstr         = ref($subSampleFiles);
      if ( $refstr eq 'HASH' ) {
        foreach my $groupName ( sort keys %$subSampleFiles ) {
          my $groupSampleNames = $subSampleFiles->{$groupName};
          if (is_string($groupSampleNames)){
            if($name_first){
              print $list "${sample_name}${delimiter}$groupSampleNames${delimiter}${groupName}\n";
            }else{
              print $list "$groupSampleNames${delimiter}${groupName}${delimiter}${sample_name}\n";
            }
          }else{
            for my $groupSampleName (@$groupSampleNames) {
              if($name_first){
                print $list "${sample_name}${delimiter}${groupSampleName}${delimiter}${groupName}\n";
              }else{
                print $list "${groupSampleName}${delimiter}${groupName}${delimiter}${sample_name}\n";
              }
            }
          }
        }
      }
      elsif ( $refstr eq 'ARRAY' ) {
        foreach my $subSampleFile (@$subSampleFiles) {
          my $curSampleName = $sample_name;
          if ( scalar(@outputNames) > 0 ) {
            $nameIndex     = $nameIndex + 1;
            $curSampleName = $outputNames[$nameIndex];
          }
          if($name_first){
            print $list $subSampleFile . "$curSampleName${delimiter}\n";
          }else{
            print $list $subSampleFile . "${delimiter}$curSampleName\n";
          }
        }
      }
      else {
        if(!defined $subSampleFiles) {
          print $list "${delimiter}$sample_name\n";
        }else{
          print $list $subSampleFiles . "${delimiter}$sample_name\n";
        }
      }
    }
    close($list);
  }
  return $result;
}

sub fix_task_name {
  my $def      = shift;
  my $taskName = $def->{task_name};
  $taskName =~ s/\s/_/g;
  $def->{task_name} = $taskName;
}

sub get_parameter_file {
  my ( $config, $section, $key ) = @_;
  my $result    = parse_param_file( $config, $section, $key, 0 );
  my $resultArg = get_option( $config, $section, $key . "_arg", "" );
  if ( !defined($result) ) {
    $result = "";
  }
  else {
    $result = "\"" . $result . "\"";
  }
  return ( $result, $resultArg );
}

sub get_parameter_sample_files {
  my ( $config, $section, $key ) = @_;
  my $result = {};
  if ( has_raw_files( $config, $section, $key ) ) {
    $result = get_raw_files( $config, $section, $key );
  }
  my $resultArg           = get_option( $config, $section, $key . "_arg",            "" );
  my $resultJoinDelimiter = get_option( $config, $section, $key . "_join_delimiter", "," );
  my $resultNameArg = get_option( $config, $section, $key . "_name_arg", "" );
  my $resultNameJoinDelimiter = get_option( $config, $section, $key . "_name_join_delimiter", "," );

  #print($key . " delimiter=" . $resultJoinDelimiter . "\n");

  return ( $result, $resultArg, $resultJoinDelimiter, $resultNameArg, $resultNameJoinDelimiter );
}

sub is_paired_end {
  my ($def) = @_;
  my $result = $def->{is_paired_end};
  if ( not defined $result ) {
    $result = $def->{is_pairend};
  }
  if ( not defined $result ) {
    $result = $def->{is_paired};
  }
  if ( not defined $result ) {
    $result = $def->{ispairend};
  }
  if ( not defined $result ) {
    $result = $def->{pairend};
  }
  return ($result);
}

sub get_is_paired_end_option {
  my ( $config, $section, $default ) = @_;

  my $curSection = get_config_section( $config, $section );

  my $result = is_paired_end($curSection);
  if ( not defined $result and defined $default ) {
    $result = $default;
  }

  return ($result);
}

sub get_rawfiles_option {
  my ( $raw_files, $option_string, $intern ) = @_;
  my $result = "";
  if ( not defined $intern ) {
    $intern = " ";
  }
  for my $sample_name ( sort keys %$raw_files ) {
    my @sample_files = @{ $raw_files->{$sample_name} };
    my $sampleFile   = $sample_files[0];
    $result = $result . $intern . $option_string . " " . $sampleFile;
  }
  return $result;
}

sub get_parameter_options {
  my ( $config, $section, $prefix, $parameters, $defaults ) = @_;

  my $curSection = get_config_section( $config, $section );
  my $result     = "";
  foreach my $i ( 0 .. ( scalar(@$parameters) - 1 ) ) {
    my $parameter = $parameters->[$i];
    my $value     = $curSection->{$parameter};
    if ( not defined $value and defined $defaults ) {
      $value = $defaults->[$i];
    }
    if ( defined $value ) {
      $result = $result . " \\\n  " . $prefix . $parameter . " " . $value;
    }
  }

  return ($result);
}

sub get_version_files {
  my ($config) = @_;

  my $result = {};
  for my $task_section ( keys %$config ) {
    if(!is_hash($config->{$task_section})) {
      next;
    }
    my $classname = $config->{$task_section}{class};
    if ( !defined $classname ) {
      next;
    }
    my $myclass = instantiate($classname);

    my $expect_file_map = $myclass->result( $config, $task_section, "version\$" );
    if (%$expect_file_map) {
      my $versionFiles = [];
      for my $key ( keys %$expect_file_map ) {
        my $expectFiles = $expect_file_map->{$key};
        push @$versionFiles, @$expectFiles;
      }
      if ( scalar(@$versionFiles) > 0 ) {
        $result->{$task_section} = $versionFiles;
      }
    }
  }
  return ($result);
}

sub ext_to_list {
  my $ext = shift;
  my $result = [];

  my @output_other_exts = split( "[,;]+", $ext );
  for my $other_ext (@output_other_exts) {
    my $trim_ext = trim($other_ext);
    if ($trim_ext ne "") {
      push( @$result, $trim_ext );
    }
  }

  return($result);
}

sub get_output_ext {
  my ( $config, $section, $defaultValue ) = @_;
  my $ext_list = get_output_ext_list($config, $section);
  if (($ext_list->[0] ne "") or (not defined $defaultValue)) {
    return $ext_list->[0];
  }else{
    return $defaultValue;
  }
}

sub get_output_ext_list {
  my ( $config, $section ) = @_;

  my $result = [];
  for my $key ("output_file_ext", "output_ext", "output_other_ext") {
    my $other_result = ext_to_list(get_option( $config, $section, $key, "" ));
    if ($other_result ne "") {
      push (@$result, @$other_result);
    }
  }

  if (scalar(@$result) == 0){
    push(@$result, "");
  }

  return ($result);
}

sub get_program_param {
  my ( $files, $arg, $joinDelimiter, $sample_name, $result_dir, $index ) = @_;
  
  my $result = "";
  my $pfiles = $files->{$sample_name};
  if(defined $pfiles){
    if (ref $pfiles eq ref {}){
      my $configfile = $result_dir . "/fileList_${index}_${sample_name}.txt";
      open( my $list, ">$configfile" ) or die "Cannot create $configfile";
      foreach my $sub_name (sort keys %$pfiles) {
        my $subSampleFiles = $pfiles->{$sub_name};
        my $refstr         = ref($subSampleFiles);
        if ( $refstr eq 'HASH' ) {
          foreach my $groupName ( sort keys %$subSampleFiles ) {
            my $groupSampleNames = $subSampleFiles->{$groupName};
            for my $groupSampleName (@$groupSampleNames) {
              print $list "${groupSampleName}\t${groupName}\t${sub_name}\n";
            }
          }
        }
        elsif ( $refstr eq 'ARRAY' ) {
          foreach my $subSampleFile (@$subSampleFiles) {
            my $curSampleName = $sample_name;
            print $list $subSampleFile . "\t$curSampleName\n";
          }
        }
        else {
          print $list $subSampleFiles . "\t$sub_name\n";
        }
      }
      close($list);
      $result = trim($arg . " " . $configfile);
    }else{
      my $pfile;
      if(is_array($pfiles)){
        $pfile  = join($joinDelimiter, @$pfiles );
      }else{
        $pfile = $pfiles;
      }
      $result = trim($arg . " " . $pfile);
    }
  }
  
  return ($result);
}

sub get_program {
  my ($config, $section) = @_;  

  my $result     = get_option( $config, $section, "program" );

  if ( get_option( $config, $section, "check_program", 1 ) ) {
    if ( !File::Spec->file_name_is_absolute($result) ) {
      $result = dirname(__FILE__) . "/$result";
    }
    if ( !( -e $result ) ) {
      print Dumper(longmess());
      die("program $result defined but not exists!");
    }
  }
  
  return ($result);
}

sub get_interation_sample_subsample_map {
  my $source_files = shift;
  my $sample_name_map      = {};
  for my $individual_sample_name ( sort keys %$source_files ) {
    my $sample_name = $individual_sample_name;
    $sample_name =~ s/_ITER_.*//g;
    if ( not defined $sample_name_map->{$sample_name} ) {
      $sample_name_map->{$sample_name} = [$individual_sample_name];
    }
    else {
      my $names = $sample_name_map->{$sample_name};
      push( @$names, $individual_sample_name );
    }
  }
  return ($sample_name_map);
}

sub get_interation_subsample_sample_map {
  my $source_files = shift;
  my $result      = {};
  for my $individual_sample_name ( sort keys %$source_files ) {
    my $sample_name = $individual_sample_name;
    $sample_name =~ s/_ITER_.*//g;
    $result->{$individual_sample_name} = $sample_name;
  }
  return ($result);
}

sub get_groups {
  my ($files, $pattern) = @_;
  my $groups = {};
  for my $name (keys %$files){
    if ( $name =~ $pattern ) {
      my $group = $1;

      my $names = $groups->{$group};
      if (defined $names) {
        push @$names, ($name);
      }else{
        $groups->{$group} = [$name];
      }
    }else{
      print Dumper(longmess());
      die "Cannot find pattern $pattern in $name.";
    }
  }

  print("  groups => {\n");
  for my $group (sort keys %$groups){
    my $names = $groups->{$group};
    my @sorted_names = sort @$names;
    print('    "' . $group . '" => ["' . join('","', @sorted_names) . '"],' . "\n");
  }
  print("  },\n");
}

sub get_covariances {
  my ($pairs, $groups, $covname, $covpattern) = @_;

  print("  \"pairs\" => {\n");
  for my $pairName (sort keys %$pairs){
    print("    \"$pairName\" => {\n");
    my $groupNames = $pairs->{$pairName}{groups};
    print('      groups => ["'. join('", "', @$groupNames) . '"], ' . "\n");
    my $covs = [];
    for my $groupName (@$groupNames){
      my $samples = $groups->{$groupName};
      for my $sampleName (@$samples){
        $sampleName =~ /$covpattern/;
        my $cov = $1;
        push(@$covs, $cov);
      }
    }
    print('      ' . $covname . ' => ["'. join('", "', @$covs) . '"] '. "\n");
    print("    }, \n");
  }
  print("  },\n");
}

sub getGroupPickResult {
  my ($config, $source_ref, $sample_index_in_group, $pattern ) = @_;

  my $temp_section = "temp_section";
  $config->{$temp_section} = {     
    "class" => "CQS::GroupPickTask",
    "source_ref" => [$source_ref],
    "groups_ref" => ["groups"],
    "sample_index_in_group" => $sample_index_in_group, 
  };
  my $myclass = instantiate("CQS::GroupPickTask");
  my $result = $myclass->result($config, $temp_section, $pattern);
  delete $config->{$temp_section};
  return ($result);
}

sub getMemoryPerThread {
  my ($memory_in_gb, $thread) = @_;
  my $result = $memory_in_gb;
  $result =~ /(\d+)(\S+)/;
  my $memNum = $1;
  $result = floor($memNum / $thread);
  my $isMB = 0;
  if ($result < 1) {
    $result = floor($result * 1024);
    $isMB = 1;
  }
  return($result, $isMB);
}

sub option_contains_arg {
  my ($option, $arg) = @_;

  if ($arg eq "") {
    return(0);
  }
  
  if (index($option, " " . $arg . " ") != -1) {
    return(1);
  }

  if (substr($option, 0, length($arg)) eq $arg){
    return(1);
  }

  if (substr($option, -length($arg)) eq $arg){
    return(1);
  }

  return(0);
}

sub left_pad {
  my ($iter, $max_length) = @_;
  return("0" x ($max_length -length("$iter")) . $iter);
}

sub get_key_name {
  my ($sample_name, $scatter_name) = @_;
  return ($sample_name . "." . $scatter_name);
}

sub get_interval_file_map {
  my ($config, $section) = @_;
  my $interval_list_file = get_option_file($config, $section, "interval_list_file");
  my $files = [];
  open(FH, '<', $interval_list_file) or die $!;
  while(<FH>){
    my $file = $_;
    chomp($file);
    push(@$files, $file);
  }

  #print(Dumper($files));

  my $result = {};
  tie %$result, 'Tie::IxHash';

  my $iter_end = scalar(@$files) - 1;
  for my $iter (0..$iter_end){
    my $key = left_pad($iter, length("$iter_end"));
    $result->{$key} = $files->[$iter]
  }
  return $result;
}

sub read_table {
  my ($filename, $name_index) = @_;
  if(not defined $name_index){
    $name_index = 0;
  }
  
  my $result = {};
  tie %$result, 'Tie::IxHash';

  open(my $fh, '<', $filename) or die $!;
  my $header = <$fh>;
  chomp($header);
  my @headers = split('\t', $header);
  map { s/^\s+|\s+$//g; } @headers;
  #print(Dumper(@headers));
  
  my $names = {};
  tie %$names, 'Tie::IxHash';
  for my $idx (0..(scalar(@headers)-1)){
    if ($idx == $name_index){
      next;
    }
    $names->{$headers[$idx]} = 1;
  }
  
  while(<$fh>){
    my $line = $_;
    #print($line);
    my @parts = split('\t', $line);
    #print(Dumper(@parts));
    my $dic = {};
    tie %$dic, 'Tie::IxHash';
    for my $idx (0..(scalar(@headers)-1)){
      if ($idx == $name_index){
        next;
      }
      
      my $column = $headers[$idx];
      my $value = $parts[$idx];
      chomp($value);
      $dic->{$column} = $value;
    }
    my $name = $parts[$name_index];
    $result->{$name} = $dic;
  }
  
  #print(Dumper($result));
  return($result, $names);
}

sub merge_hash_left_precedent {
  my ($a, $b) = @_;
  my $merge_c = Hash::Merge->new('LEFT_PRECEDENT');
  return $merge_c->merge($a, $b);
}

sub merge_hash_right_precedent {
  my ($a, $b) = @_;
  my $merge_c = Hash::Merge->new('RIGHT_PRECEDENT');
  return $merge_c->merge($a, $b);
}

sub get_all_sample_names {
  my $def = shift;

  my $files = $def->{files};

  my $sample_names = {};
  for my $sample_name (keys %$files){
    $sample_names->{$sample_name} = $sample_name;
  }

  if (getValue($def, "perform_split_hto_samples", 0)){
    my $HTO_samples = getValue($def, "HTO_samples");
    for my $hto (keys %$HTO_samples){
      my $hto_map = $HTO_samples->{$hto};
      my @hto_sub_samples = (values %$hto_map);
      delete $sample_names->{$hto};
      for my $hto_sub_sample (@hto_sub_samples) {
        $sample_names->{$hto_sub_sample} = $hto_sub_sample;
      }
    }
  }

  if(getValue($def, "pool_sample", 0)){
    my $sample_pool = getValue($def, "pool_sample_groups");
    my $sample_to_pool = {};
    for my $pool (keys %$sample_pool){
      my $samples = $sample_pool->{$pool};
      for my $s (@$samples){
        $sample_names->{$s} = $pool;
      }
    }
  }

  my @samplenames = sort(uniq(values %$sample_names));
  return(\@samplenames);
}

sub get_groups_by_pattern_dic {
  my ($def, $gpattern_dic) = @_;
  
  if(!defined $gpattern_dic){
    $gpattern_dic = $def->{groups_pattern};
  }

  my $samplenames = get_all_sample_names($def);

  my $groups = {};
  for my $groupname (sort keys %$gpattern_dic){
    my $gpattern = $gpattern_dic->{$groupname};
    for my $samplename (@$samplenames) {
      if($samplename =~ /$gpattern/){
        if (not defined $groups->{$groupname}){
          $groups->{$groupname} = [$samplename];
        }
        else{
          my $samples = $groups->{$groupname};
          push (@$samples, $samplename);
        }
      }
    }
  }
  return ($groups);
}

sub get_groups_by_pattern_value {
  my ($def, $gpattern) = @_;
  
  if(!defined $gpattern){
    $gpattern = $def->{groups_pattern};
  }

  my $samplenames = get_all_sample_names($def);
  
  #print($gpattern);
  #print(Dumper($files));
  my $groups = {};
  for my $samplename (@$samplenames) {
    my $groupname = $samplename;
    if($samplename =~ /$gpattern/){
      $groupname = $1;
      if(defined $2){
        $groupname = $groupname . $2;
      }
      if(defined $3){
        $groupname = $groupname . $3;
      }
      #print($groupname . " : " . $samplename . "\n");
      if (not defined $groups->{$groupname}){
        $groups->{$groupname} = [$samplename];
      }
      else{
        my $samples = $groups->{$groupname};
        push (@$samples, $samplename);
      }
    }
  }
  return ($groups);
}

sub get_groups_by_pattern_array {
  my ($def, $gpatterns) = @_;
  
  if(!defined $gpatterns){
    $gpatterns = $def->{groups_pattern};
  }

  my $groups = {};
  for my $gpattern (@$gpatterns){
    my $subgroups = get_groups_by_pattern_value($def, $gpattern);
    $groups = merge_hash_right_precedent($groups, $subgroups);
  }
  return($groups);
}

sub get_groups_by_pattern {
  my ($def, $gpattern) = @_;
  if(not defined $gpattern){
    $gpattern = $def->{groups_pattern};
  }
  if (is_hash($gpattern)){
    return(get_groups_by_pattern_dic($def, $gpattern));
  }elsif (is_array($gpattern)){
    return(get_groups_by_pattern_array($def, $gpattern));
  }else{
    return(get_groups_by_pattern_value($def, $gpattern));
  }
}

sub get_groups_by_file {
  my ($group_file) = @_;

  my $result = {};
  tie %$result, 'Tie::IxHash';

  open(my $fh, '<', $group_file) or die $!;
  
  while(<$fh>){
    my $line = $_;
    chomp($line);
    #print($line);
    my @parts = split('\t', $line);
    #print(Dumper(@parts));
    my $value = $parts[0];
    my $key = $parts[1];

    if (not defined $result->{$key}){
      $result->{$key} = [$value];
    }else{
      my $values = $result->{$key};
      push(@$values, $value);
      $result->{$key} = $values;
    }
  }
  
  return($result);
}

sub get_unique_groups {
  my $groups = shift;
  my $file_dic = {};
  for my $g (sort keys %$groups){
    my $samples = $groups->{$g};
    for my $s (@$samples){
      if(not defined $file_dic->{$s}){
        $file_dic->{$s} = [$g];
      }else{
        my $cur_groups= $file_dic->{$s};
        push @$cur_groups, $g;
        $file_dic->{$s} = $cur_groups;
      }
    }
  }
  my $result = {};
  for my $s (sort keys %$file_dic){
    my $groups = $file_dic->{$s};
    my $groups_str = join(":", sort @$groups);
    if(not defined $result->{$groups_str}){
      $result->{$groups_str} = [$s];
    }else{
      my $cur_files= $result->{$groups_str};
      push @$cur_files, $s;
      $result->{$groups_str} = $cur_files;
    }
  }

  return($result);
}

sub get_correlation_groups_by_pattern {
  my ($def) = @_;

  if(defined $def->{correlation_groups}){
    return($def->{correlation_groups});
  }

  my $gpattern_dic = $def->{correlation_groups_pattern_dic};
  if(defined $gpattern_dic){
    if(not is_hash($gpattern_dic)){
      print Dumper(longmess());
      die "correlation_groups_pattern has to be hash";
    }

    my $result = {};
    for my $gpattern_key (sort keys %$gpattern_dic){
      my $gpattern = $gpattern_dic->{$gpattern_key};
      $result->{$gpattern_key} = get_groups_by_pattern($def, $gpattern);
    }
    return($result);
  }else{
    if(defined $def->{unique_groups}){
      return({all => $def->{unique_groups}});
    }else{
      return({all => get_unique_groups($def->{groups})});
    }
  }
}

sub get_covariances_by_pattern {
  my ($def) = @_;

  my $files = $def->{files};
  if(getValue($def, "perform_split_hto_samples", 0)){
    $files = {};
    my $HTO_samples = getValue($def, "HTO_samples");
    for my $exp (keys %$HTO_samples){
      my $samples = $HTO_samples->{$exp};
      for my $sample (values %$samples){
        $files->{$sample} = [];
      }
    }
  }
  my $covariance_patterns = $def->{covariance_patterns};
  my $covariances = [sort keys %$covariance_patterns];
  my $samplenames = [sort keys %$files];
  my $cov_map = {};
  for my $covariance (@$covariances) {
    my $cov_pattern_def = $covariance_patterns->{$covariance};

    my $cov_pattern;
    my $cov_prefix;
    if (ref $cov_pattern_def eq 'HASH'){
      $cov_pattern = getValue($cov_pattern_def, "pattern");
      $cov_prefix = getValue($cov_pattern_def, "prefix", "");
    }
    else{
      $cov_pattern = $cov_pattern_def;
      $cov_prefix = "";
    }

    $cov_map->{$covariance} = {};
    for my $samplename (@$samplenames) {
      my $cov_value = $samplename;
      if ($samplename =~ /$cov_pattern/){
        $cov_value = $1;
      }
      $cov_map->{$covariance}{$samplename} = $cov_prefix . $cov_value;
    }
  }
  return ($cov_map, $covariances, $samplenames);
}


sub write_covariance_file {
  my ($def, $cov_map, $covariances, $samplenames) = @_;

  my $target_dir = getValue($def, "target_dir");
  my $cov_file = $target_dir . "/covariance.txt";
  open( my $cov, ">$cov_file" ) or die "Cannot create $cov_file";
  print $cov "Sample";
  for my $covariance (@$covariances) {
    print $cov "\t" . $covariance;
  }
  print $cov "\n";
  for my $samplename (@$samplenames) {
    print $cov "$samplename";
    for my $covariance (@$covariances) {
      print $cov "\t" . $cov_map->{$covariance}{$samplename};
    }
    print $cov "\n";
  }
  close($cov);

  return($cov_file);
}

sub create_covariance_file_by_pattern {
  my ($def) = @_;

  my ($cov_map, $covariances, $samplenames) = get_covariances_by_pattern($def);

  my $cov_file = write_covariance_file($def, $cov_map, $covariances, $samplenames);

  return($cov_file);
}

sub get_covariances_by_file_pattern {
  my ($def) = @_;

  my $files = $def->{files};
  my $covariance_patterns = $def->{covariance_file_patterns};
  my $covariances = [sort keys %$covariance_patterns];
  my $samplenames = [sort keys %$files];
  my $cov_map = {};
  for my $covariance (@$covariances) {
    my $cov_pattern_def = $covariance_patterns->{$covariance};

    my $cov_pattern;
    my $cov_prefix;
    if (ref $cov_pattern_def eq 'HASH'){
      $cov_pattern = getValue($cov_pattern_def, "pattern");
      $cov_prefix = getValue($cov_pattern_def, "prefix", "");
    }
    else{
      $cov_pattern = $cov_pattern_def;
      $cov_prefix = "";
    }

    $cov_map->{$covariance} = {};
    for my $samplename (@$samplenames) {
      my $cov_value = $samplename;
      my $file_path = $files->{$samplename}->[0];
      if ($file_path =~ /$cov_pattern/){
        $cov_value = $1;
      }
      $cov_map->{$covariance}{$samplename} = $cov_prefix . $cov_value;
    }
  }
  return ($cov_map, $covariances, $samplenames);
}

sub create_covariance_file_by_file_pattern {
  my ($def) = @_;

  my ($cov_map, $covariances, $samplenames) = get_covariances_by_file_pattern($def);

  my $cov_file = write_covariance_file($def, $cov_map, $covariances, $samplenames);

  return($cov_file);
}

sub write_HTO_sample_file {
  my ($def) = @_;

  my $HTO_samples = $def->{HTO_samples};
  my $target_dir = getValue($def, "target_dir");
  my $cov_file = $target_dir . "/hto_samples.txt";
  open( my $cov, ">$cov_file" ) or die "Cannot create $cov_file";
  print $cov "File\tTagname\tSample\n";
  for my $fname (sort keys %$HTO_samples) {
    my $htos = $HTO_samples->{$fname};
    for my $hto (sort keys %$htos) {
      my $sample = $htos->{$hto};
      print $cov "${fname}\t${hto}\t${sample}\n";
    }
  }
  close($cov);

  return($cov_file);
}

sub get_parameter_file_option {
  my ($config, $section, $option) = @_;
  my $result = $option;
  for my $index (1..10){
    my $key = "parameterFile" . $index;
    my $parameterFile = parse_param_file( $config, $section, $key, 0 );

    if (defined($parameterFile)){
      my $param_key = "__${key}__";
      if ($result =~ /$param_key/){
        $result =~ s/$param_key/$parameterFile/g;
      }else{
        my $parameterFileArg = get_option($config, $section, "${key}_arg", "");
        $result = $result . " " . $parameterFileArg . " " . $parameterFile;
      }
    }
  }
  return($result);
}

sub get_hash_level2 {
  my ($hash, $level2_key) = @_;
  my $result = {};
  for my $level1_key (keys %$hash){
    my $value = $hash->{$level1_key};
    $result->{$level1_key} = $value->{$level2_key};
  }
  return($result);
}

sub get_expanded_genes {
  my $curated_gene_def = shift;
  my $result = {};

  if (is_array($curated_gene_def)) {
    $result->{"interesting_genes"} = {
      clusters => [],
      genes => $curated_gene_def
    };
    return($result);
  }

  if(is_hash($curated_gene_def)) {
    for my $key (keys %$curated_gene_def) {
      my $value = $curated_gene_def->{$key};
      if (is_array($value)) {
        $result->{$key} = {
          clusters => [],
          genes => $value
        };
      }else{
        $result->{$key} = $value;
      }
    }

    return($result);
  }

  print Dumper(longmess());
  die "curated_gene_def should be either array or hash";
}

sub parse_curated_genes {
  my $curated_gene_def = shift;
  my $result = get_expanded_genes($curated_gene_def);

  my $clusters = get_hash_level2($result, "clusters");
  my $genes = get_hash_level2($result, "genes");

  return($result, $clusters, $genes)
}

sub get_joined_files {
  my ( $files, $join_delimiter ) = @_;
  my $pfiles                  = [];
  for my $individual_sample_name (sort keys %$files) {
    my $p_invividual_files = $files->{$individual_sample_name};
    my $p_invividual_file  = $p_invividual_files->[0];
    push( @$pfiles, $p_invividual_file );
  }
  my $result = join( $join_delimiter, @$pfiles );
  return ($result);
}

sub get_joined_names {
  my ( $files, $join_delimiter ) = @_;
  my @pfiles = sort keys %$files;
  my $result = join( $join_delimiter, @pfiles );
  return ($result);
}

sub process_parameter_sample_file {
  my ($config, $section, $result_dir, $task_name, $task_suffix, $cur_option, $source_key, $index) = @_;
  my ( $fileMap, $fileArg, $fileJoinDelimiter, $nameArg, $nameJoinDelimiter) = get_parameter_sample_files( $config, $section, $source_key );

  if($index == 1){
    if ($cur_option =~ /__SAMPLE_NAMES__/){
      my $sample_names = get_joined_names($fileMap, $nameJoinDelimiter);
      $cur_option =~ s/__SAMPLE_NAMES__/$sample_names/g;
    }
  }

  my $input = "";
  if (defined $config->{$section}{$source_key . "_type"} && ($config->{$section}{$source_key . "_type"} eq "array")){
    $input = get_joined_files($fileMap, $fileJoinDelimiter);

    if((defined $nameArg) && ($nameArg ne "")){
      my $sample_names = get_joined_names($fileMap, $nameJoinDelimiter);
      if ($config->{$section}{$source_key . "_name_has_comma"}){
        $input = $input . " " . $nameArg . " \"" . $sample_names . "\"";
      }else{
        $input = $input . " " . $nameArg . " " . $sample_names;
      }
    }
  }else{
    my $no_prefix = get_option( $config, $section, "no_prefix", 0 );

    my $filelist_name;
    if($no_prefix){
      $filelist_name = "fileList${index}.list";
    }else{
      $filelist_name = "${task_name}_${task_suffix}_fileList${index}.list";
    }
    my $list_file = save_parameter_sample_file( $config, $section, $source_key, "${result_dir}/${filelist_name}" );

    if($list_file ne ""){
      $input = basename($list_file);
    }
  }

  my $file_key = $index == 1 ? "__FILE__" : "__FILE${index}__";
  if ($cur_option =~ /$file_key/){
    $cur_option =~ s/$file_key/$input/g;
  } elsif (option_contains_arg($cur_option, $fileArg)) {
  } else{
    $cur_option = $cur_option . " " . $fileArg . " " . $input;
  }
  return($cur_option, $input);
}

sub get_task_dep_pbs_map {
  my ($config, $task_section_name) = @_;
  my $task_section = $config->{$task_section_name};
  
  if ( not defined $task_section->{class} ) {
    next
  }
  
  my $myclass = instantiate( $task_section->{class} );
  my $pbs_sample_map = $myclass->get_pbs_source( $config, $task_section_name );

  my $allSampleNameMap = {};
  for my $pbs (keys %$pbs_sample_map){
    my $sample_names = $pbs_sample_map->{$pbs};
    for my $sample_name (@$sample_names){
      $allSampleNameMap->{$sample_name} = 1;
    }
  }
  my @allSampleNames = (sort keys %$allSampleNameMap);

  my $taskdeppbsmap = {};
  for my $key ( keys %$task_section ) {
    if ($key eq "gather_name_ref"){
      next;
    }

    my $mapname = $key;
    if ( $mapname =~ /_ref$/ ) {
      $mapname =~ s/_config_ref//g;
      $mapname =~ s/_ref//g;
      my $refpbsmap = get_ref_section_pbs( $config, $task_section_name, $mapname );
      my @refNames = (sort keys %$refpbsmap);

      my $keyEquals = 1;
      if (scalar(@allSampleNames) != scalar(@refNames)){
        $keyEquals = 0;
      }else{
        for my $idx (0..(scalar(@allSampleNames)-1)){
          if ($allSampleNames[$idx] ne $refNames[$idx]) {
            $keyEquals = 0;
            last;
          }
        }
      }

      for my $pbs (keys %$pbs_sample_map){
        my $curpbs = $taskdeppbsmap->{$pbs};
        if ( !defined $curpbs ) {
          $curpbs = {};
        }
        
        if ($keyEquals) {
          my $sample_names = $pbs_sample_map->{$pbs};
          for my $sample_name (@$sample_names) {
            my $ref_pbs_list = $refpbsmap->{$sample_name};
            if (defined $ref_pbs_list){
              for my $ref_pbs (@$ref_pbs_list){
                $curpbs->{$ref_pbs} = 1;
              }
            }
          }
        }else{
          for my $sample_name (sort keys %$refpbsmap){
            my $ref_pbs_list = $refpbsmap->{$sample_name};
            for my $ref_pbs (@$ref_pbs_list){
              $curpbs->{$ref_pbs} = 1;
            }
          }
        }

        $taskdeppbsmap->{$pbs} = $curpbs;
      }
    }
  }

  return($taskdeppbsmap);
}

sub add_bind{
  my ($config, $bind) = @_;

  for my $section_key (keys %$config){
    my $section = $config->{$section_key};
    if (is_hash($section)){
      for my $key (keys %$section){
        if ($key =~ /docker_command$/){
          my $value = $section->{$key};
          $value =~ s/exec/exec -B $bind/;
          $section->{$key} = $value;
          #print($value . "\n");
        }
      }
    }
  }

  return($config);
}

sub get_final_file_by_task_name {
  my ($all_results, $task_name) = @_;
  my $result_files = $all_results->{$task_name};
  my $final_file;
  if((not defined $result_files) || (scalar(@$result_files) == 0)){
    #get last sample_name result file
    for my $cur_key (sort keys %$all_results){
      if($cur_key eq $task_name){
        next;
      }

      my $cur_files = $all_results->{$cur_key};
      my @no_filelists = grep(!/.filelist$/, @$cur_files);
      $final_file = $no_filelists[-1];
      #print($final_file . "\n");
    }
  }else{
    #get last task_name result file
    my @no_filelists = grep(!/.filelist$/, @$result_files);
    $final_file = $no_filelists[-1];
  }
  return($final_file);
}

sub get_groups_from_covariance_file {
  my ($covariance_file, $sample_column, $group_column) = @_;
  my $result = {};
  my $aoh;
  if($covariance_file =~ /.csv$/){
    $aoh = csv(in => $covariance_file, headers => "auto"); 
  }else{
    $aoh = csv(in => $covariance_file, sep_char => "\t", headers => "auto"); 
  }
  for my $dic (@$aoh){
    my $cur_sample = $dic->{$sample_column};
    my $cur_group = $dic->{$group_column};
    #print($cur_sample . " => " . $cur_group . "\n");
    if (!defined $result->{$cur_group}){
      $result->{$cur_group} = [$cur_sample];
    }else{
      my $old_samples = $result->{$cur_group};
      push(@$old_samples, $cur_sample);
    }
  }
  return($result);
}

1;
