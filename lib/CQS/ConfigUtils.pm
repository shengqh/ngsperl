#!/usr/bin/perl
package CQS::ConfigUtils;

use strict;
use warnings;
use File::Basename;
use CQS::FileUtils;
use CQS::PBS;
use CQS::ClassFactory;
use CQS::StringUtils;
use CQS::CQSDebug;
use Data::Dumper;
use Hash::Merge qw( merge );
use Tie::IxHash;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(get_config_section has_config_section get_option get_java get_cluster get_parameter get_param_file get_directory parse_param_file has_raw_files get_raw_files_and_keys get_raw_files get_raw_files_keys get_raw_files_attributes get_raw_files2 get_run_command get_option_value get_pair_groups
      get_pair_groups_names get_cqstools get_group_sample_map get_group_samplefile_map get_group_samplefile_map_key save_parameter_sample_file saveConfig writeFileList initDefaultValue get_pure_pairs writeParameterSampleFile)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

####
# Get section definition from config definition. section can be multiple layers, such as "task1::task2" which indicates $config->{task1}->{task2}
####
sub get_config_section {
  my ( $config, $section ) = @_;
  my @sections = split( '::', $section );
  my $result = $config;
  for my $curSection (@sections) {
    $result = $result->{$curSection};
    die "Cannot find section $section!" if ( !defined $result );
  }
  return ($result);
}

sub has_config_section {
  my ( $config, $section ) = @_;
  my @sections = split( '::', $section );
  my $result = $config;
  for my $curSection (@sections) {
    $result = $result->{$curSection};
    return 0 if ( !defined $result );
  }
  return 1;
}

sub get_option {
  my ( $config, $section, $key, $default ) = @_;

  my $curSection = get_config_section( $config, $section );

  my $result = $curSection->{$key};
  if ( !defined $result ) {
    if ( !defined $default ) {
      die "Define ${section}::${key} first!";
    }
    else {
      $result = $default;
    }
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

sub get_parameter {
  my ( $config, $section, $create_directory ) = @_;

  my $curSection = get_config_section( $config, $section );

  $create_directory = 1 if !defined($create_directory);

  my $task_name = get_option( $config, $section, "task_name", "" );
  if ( $task_name eq "" ) {
    $task_name = get_option( $config, "general", "task_name" );
  }

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
  my ( $log_dir, $pbs_dir, $result_dir ) = init_dir( $target_dir, $create_directory );
  my ($pbs_desc) = $cluster->get_cluster_desc($refPbs);

  my $option    = get_option( $config, $section, "option",    "" );
  my $sh_direct = get_option( $config, $section, "sh_direct", 0 );

  my $init_command = get_option( $config, $section, "init_command", "" );

  if ($sh_direct) {
    $sh_direct = "bash";
  }
  else {
    $sh_direct = $cluster->get_submit_command();
  }

  my $thread = $cluster->get_cluster_thread($refPbs);
  my $memory = $cluster->get_cluster_memory($refPbs);

  return ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command );
}

#get parameter which indicates a file. If required, not defined or not exists, die. If defined but not exists, die.
#returned file either undefined or exists.
sub get_param_file {
  my ( $file, $name, $required ) = @_;

  my $result = $file;

  if ($required) {
    if ( !defined $file ) {
      die "$name was not defined!";
    }

    if ( !is_debug() && !-e $file ) {
      die "$name $file defined but not exists!";
    }
  }
  else {
    if ( defined($file) ) {
      if ( $file eq "" ) {
        undef($result);
      }
      elsif ( !is_debug() && !-e $file ) {
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

  die "parameter name must be defined!" if !defined $name;

  my $result = $curSection->{$name};

  if ($required) {
    if ( !defined $result ) {
      die "$name was not defined!";
    }

    if ( !is_debug() && !-e $result ) {
      die "$name $result defined but not exists!";
    }
  }
  else {
    if ( defined($result) ) {
      if ( $result eq "" ) {
        undef($result);
      }
      elsif ( !is_debug() && !-e $result ) {
        die "$name $result defined but not exists!";
      }
    }
  }
  return ($result);
}

sub get_cqstools {
  my ( $config, $section, $required ) = @_;
  my $curSection = get_config_section( $config, $section );

  my $cqstools = get_param_file( $curSection->{cqs_tools}, "cqs_tools", 0 );
  if ( !defined $cqstools ) {
    $cqstools = get_param_file( $curSection->{cqstools}, "cqstools", $required );
  }
  return ($cqstools);
}

sub parse_param_file {
  my ( $config, $section, $key, $required ) = @_;

  my $curSection = get_config_section( $config, $section );

  die "parameter key must be defined!" if !defined $key;

  if ( defined $curSection->{$key} ) {
    return $curSection->{$key};
  }

  my $key_ref = $key . "_ref";
  if ( defined $curSection->{$key_ref} ) {
    my $refSectionName = $curSection->{$key_ref};
    my $pattern;
    if ( ref($refSectionName) eq 'ARRAY' ) {
      my @parts = @{$refSectionName};
      if ( scalar(@parts) == 2 ) {
        $pattern        = $parts[1];
        $refSectionName = $parts[0];
      }
      else {
        $refSectionName = $parts[0];
      }
    }

    my $refSection = get_config_section( $config, $refSectionName );
    if ( defined $refSection->{class} ) {
      my $myclass = instantiate( $refSection->{class} );
      my $result = $myclass->result( $config, $refSectionName, $pattern, 1 );
      foreach my $k ( sort keys %{$result} ) {
        my @files = @{ $result->{$k} };
        if ( scalar(@files) > 0 ) {
          return $files[0];
        }
      }
      die "section $refSectionName return nothing!";
    }
  }

  if ($required) {
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

sub do_get_unsorted_raw_files {
  my ( $config, $section, $returnself, $mapname, $pattern ) = @_;

  my $curSection = get_config_section( $config, $section );

  if ( !defined $mapname ) {
    $mapname = "source";
  }
  my $mapname_ref        = $mapname . "_ref";
  my $mapname_config_ref = $mapname . "_config_ref";

  if ( defined $curSection->{$mapname} ) {
    return ( $curSection->{$mapname}, 1 );
  }

  if ( defined $curSection->{$mapname_ref} || defined $curSection->{$mapname_config_ref} ) {
    my $refmap = {};
    if ( defined $curSection->{$mapname_ref} ) {

      #in same config
      my $targetSection = $curSection->{$mapname_ref};

      if ( ref($targetSection) eq 'HASH' ) {
        return ( $targetSection, 1 );
      }

      if ( ref($targetSection) eq 'ARRAY' ) {
        my @parts      = @{$targetSection};
        my $partlength = scalar(@parts);
        for ( my $index = 0 ; $index < $partlength ; ) {
          if ( !has_config_section( $config, $parts[$index] ) ) {
            die "undefined section $parts[$index]";
          }
          get_config_section( $config, $parts[$index] );

          if ( $index == ( $partlength - 1 ) || has_config_section( $config, $parts[ $index + 1 ] ) ) {
            $refmap->{$index} = { config => $config, section => $parts[$index], pattern => $pattern };
            $index++;
          }
          else {
            $refmap->{$index} = { config => $config, section => $parts[$index], pattern => $parts[ $index + 1 ] };
            $index += 2;
          }
        }
      }
      else {
        if ( !has_config_section( $config, $targetSection ) ) {
          die "undefined section $targetSection";
        }
        $refmap->{1} = { config => $config, section => $targetSection, pattern => $pattern };
      }
    }
    else {

      #in another config, has to be array
      my $refSectionName = $curSection->{$mapname_config_ref};
      if ( !( ref($refSectionName) eq 'ARRAY' ) ) {
        die "$mapname_config_ref has to be defined as ARRAY with [config, section, pattern]";
      }
      my @parts      = @{$refSectionName};
      my $partlength = scalar(@parts);
      for ( my $index = 0 ; $index < $partlength - 1 ; ) {
        my $targetConfig  = $parts[$index];
        my $targetSection = $parts[ $index + 1 ];

        if ( !( ref($targetConfig) eq 'HASH' ) ) {
          die
"$mapname_config_ref has to be defined as ARRAY with [config1, section1, pattern1,config2, section2, pattern2] or [config1, section1,config2, section2] format. config should be hash and section should be string";
        }

        if ( !has_config_section( $targetConfig, $targetSection ) ) {
          die "undefined section $targetSection in $mapname_config_ref of $section";
        }

        if ( $index == ( $partlength - 2 ) || ref( $parts[ $index + 2 ] ) eq 'HASH' ) {
          $refmap->{$index} = { config => $targetConfig, section => $targetSection, pattern => $pattern };
          $index += 2;
        }
        else {
          $refmap->{$index} = { config => $targetConfig, section => $targetSection, pattern => $parts[ $index + 2 ] };
          $index += 3;
        }
      }
    }

    #print Dumper($refmap);

    my %result = ();
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
        %myres = %{ $myclass->result( $targetConfig, $section, $pattern ) };
      }
      else {
        my ( $res, $issource ) = do_get_unsorted_raw_files( $targetConfig, $section, 1, undef, $pattern );
        %myres = %{$res};
      }

      my $refcount = keys %myres;
      for my $mykey ( keys %myres ) {
        my $myvalues = $myres{$mykey};
        if ( ref($myvalues) eq '' ) {
          $myvalues = [$myvalues];
        }

        if ( ( ref($myvalues) eq 'ARRAY' ) && ( scalar( @{$myvalues} ) > 0 ) ) {
          if ( exists $result{$mykey} ) {
            my $oldvalues = $result{$mykey};
            if ( ref($oldvalues) eq 'ARRAY' ) {
              my @merged = ( @{$oldvalues}, @{$myvalues} );

              #print "merged ARRAY ", Dumper(\@merged);
              $result{$mykey} = \@merged;
            }
            else {
              die "The source of $section->$mapname should be all HASH or all ARRAY";
            }
          }
          else {
            $result{$mykey} = $myvalues;
          }
        }

        if ( ( ref($myvalues) eq 'HASH' ) && ( scalar( keys %{$myvalues} ) > 0 ) ) {
          if ( exists $result{$mykey} ) {
            my $oldvalues = $result{$mykey};
            if ( ref($oldvalues) eq 'HASH' ) {
              $result{$mykey} = merge( $oldvalues, $myvalues );
            }
            else {
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
    die "define $mapname or $mapname_ref or $mapname_config_ref for $section";
  }
}

sub do_get_raw_files_keys {
  my $resultUnsorted = shift;
  my @keys = grep { $_ !~ /^\./ } keys %$resultUnsorted;
  my @result;
  if ( exists $resultUnsorted->{".order"} ) {
    my $orders = $resultUnsorted->{".order"};
    @result = @$orders;
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
  my @keys = grep { $_ !~ /^\./ } keys %$resultUnsorted;
  my %result;
  tie %result, 'Tie::IxHash';
  my @orderedKeys;
  if ( exists $resultUnsorted->{".order"} ) {
    my $orders = $resultUnsorted->{".order"};
    @orderedKeys = @$orders;
    die "number of key defined in .order not equals to actual keys for @_" if ( scalar(@orderedKeys) != scalar(@keys) );
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
  my ( $config, $section, $mapname, $pattern ) = @_;
  my ( $result, $issource ) = do_get_raw_files( $config, $section, 0, $mapname, $pattern );
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
  if ( ref($tmpgroup_names) eq 'HASH' ) {
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
  if ( ref($tmpgroup_names) eq 'HASH' ) {
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
    if ( ref($tmpgroup_names) eq 'HASH' ) {
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
  my $groups = get_raw_files( $config, $section, "groups" );
  my %group_sample_map = ();
  for my $group_name ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$group_name} };
    my @gfiles  = ();
    foreach my $sample_name (@samples) {
      my @bam_files = @{ $raw_files->{$sample_name} };
      my @sambam = ( $sample_name, @bam_files );
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

sub get_group_samplefile_map_key {
  my ( $config, $section, $sample_pattern, $group_key ) = @_;

  my $raw_files = get_raw_files( $config, $section, "source", $sample_pattern );
  my $groups = get_raw_files( $config, $section, $group_key );
  my %group_sample_map = ();
  for my $group_name ( sort keys %{$groups} ) {
    my @gfiles        = ();
    my $group_samples = $groups->{$group_name};
    if ( ref $group_samples eq ref "" ) {
      push( @gfiles, $group_samples );
    }
    else {
      my @samples = @{$group_samples};
      foreach my $sample_name (@samples) {
        my @bam_files = @{ $raw_files->{$sample_name} };
        foreach my $bam_file (@bam_files) {
          push( @gfiles, $bam_file );
        }
      }
    }
    $group_sample_map{$group_name} = \@gfiles;
  }
  return \%group_sample_map;
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
    open( my $list, '>', $outputFile ) or die "Cannot create $outputFile";
    foreach my $sample_name (@orderedSampleNames) {
      foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
        print $list $subSampleFile . "\t$sample_name\n";
      }
    }
    close($list);
    return ($outputFile);
  }
  else {
    return ("");
  }
}

sub saveConfig {
  my ( $def, $config ) = @_;

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
  my ( $fileName, $fileMap, $exportAllFiles ) = @_;
  open( my $fl, ">$fileName" ) or die "Cannot create $fileName";
  for my $sample_name ( sort keys %$fileMap ) {
    my @files = @{ $fileMap->{$sample_name} };
    if ( defined $exportAllFiles and $exportAllFiles ) {
      for my $eachFile (@files) {
        print $fl $sample_name, "\t", $eachFile, "\n";
      }
    }
    else {
      print $fl $sample_name, "\t", $files[0], "\n";
    }
  }
  close($fl);
}

sub initDefaultValue {
  my ( $def, $name, $defaultValue ) = @_;
  if ( !defined $def->{$name} ) {
    if ( !defined $defaultValue ) {
      die "defaultValue cannot be undefined for $name.";
    }
    $def->{$name} = $defaultValue;
  }
  return $def;
}

sub writeParameterSampleFile {
  my ( $config, $section, $resultDir, $index ) = @_;
  my $result      = "";
  my $task_suffix = get_option( $config, $section, "suffix", "" );
  my $key         = "parameterSampleFile" . $index;
  if ( has_raw_files( $config, $section, $key ) ) {
    my $temp = get_raw_files( $config, $section, $key );
    my @orderedSampleNames;
    my $keyOrder                 = $key . "Order";
    my $parameterSampleFileOrder = $config->{$section}{$keyOrder};
    if ( defined $parameterSampleFileOrder ) {
      @orderedSampleNames = @{$parameterSampleFileOrder};
    }
    else {
      @orderedSampleNames = keys %$temp;
    }
    $result = "fileList${index}${task_suffix}.txt";
    open( my $list, ">$resultDir/$result" ) or die "Cannot create $result";
    foreach my $sample_name (@orderedSampleNames) {
      foreach my $subSampleFile ( @{ ${$temp}{$sample_name} } ) {
        print $list $subSampleFile . "\t$sample_name\n";
      }
    }
    close($list);
  }
  return $result;
}

1;
