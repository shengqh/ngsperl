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

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(get_option get_java get_cluster get_parameter get_param_file parse_param_file get_raw_files get_raw_files2 get_run_command get_option_value get_pair_groups get_pair_groups_names get_cqstools)]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub get_option {
  my ( $config, $section, $key, $default ) = @_;

  die "no section $section found!" if !defined $config->{$section};

  my $result = $config->{$section}{$key};
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

  my $cluster_name;
  if ( defined $config->{$section}{cluster} ) {
    $cluster_name = get_option_value( $config->{$section}{cluster}, "slurm" );
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

  my $result;
  if ( defined $config->{$section}{$name} ) {
    $result = get_option_value( $config->{$section}{$name}, $defaultvalue );
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
  my ( $config, $section ) = @_;

  die "no section $section found!" if !defined $config->{$section};

  my $task_name = get_option( $config, "general", "task_name" );

  my $cluster = get_cluster(@_);

  my $path_file = get_param_file( $config->{$section}{path_file}, "path_file", 0 );
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
  my ( $logDir, $pbsDir, $resultDir ) = init_dir($target_dir);
  my ($pbsDesc) = $cluster->get_cluster_desc($refPbs);

  my $option    = get_option( $config, $section, "option",    "" );
  my $sh_direct = get_option( $config, $section, "sh_direct", 0 );

  if ($sh_direct) {
    $sh_direct = "bash";
  }
  else {
    $sh_direct = $cluster->get_submit_command();
  }

  my $thread = $cluster->get_cluster_thread($refPbs);

  return ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread );
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

sub get_cqstools {
  my ( $config, $section, $required ) = @_;
  my $cqsFile = get_param_file( $config->{$section}{cqs_tools}, "cqs_tools", 0 );
  if ( !defined $cqsFile ) {
    $cqsFile = get_param_file( $config->{$section}{cqstools}, "cqstools", $required );
  }
  return ($cqsFile);
}

sub parse_param_file {
  my ( $config, $section, $key, $required ) = @_;

  die "section $section was not defined!" if !defined $config->{$section};
  die "parameter key must be defined!" if !defined $key;

  if ( defined $config->{$section}{$key} ) {
    return $config->{$section}{$key};
  }

  my $key_ref = $key . "_ref";
  if ( defined $config->{$section}{$key_ref} ) {
    my $refSectionName = $config->{$section}{$key_ref};
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
    die "section $refSectionName was not defined!" if !defined $config->{$refSectionName};
    if ( defined $config->{$refSectionName}{class} ) {
      my $myclass = instantiate( $config->{$refSectionName}{class} );
      my $result = $myclass->result( $config, $refSectionName, $pattern );
      foreach my $k ( sort keys %{$result} ) {
        my @files = @{ $result->{$k} };
        return $files[0];
      }
    }
  }

  if ($required) {
    die "define ${section}::${key} first.";
  }

  return undef;
}

sub do_get_raw_files {
  my ( $config, $section, $returnself, $mapname, $pattern ) = @_;

  die "section $section was not defined!" if !defined $config->{$section};

  if ( !defined $mapname ) {
    $mapname = "source";
  }
  my $mapname_ref        = $mapname . "_ref";
  my $mapname_config_ref = $mapname . "_config_ref";

  if ( defined $config->{$section}{$mapname} ) {
    return ( $config->{$section}{$mapname}, 1 );
  }

  if ( defined $config->{$section}{$mapname_ref} || defined $config->{$section}{$mapname_config_ref} ) {
    my $refmap = {};
    if ( defined $config->{$section}{$mapname_ref} ) {

      #in same config
      my $targetSection = $config->{$section}{$mapname_ref};

      if ( ref($targetSection) eq 'HASH' ) {
        return ( $targetSection, 1 );
      }

      if ( ref($targetSection) eq 'ARRAY' ) {
        my @parts      = @{$targetSection};
        my $partlength = scalar(@parts);
        for ( my $index = 0 ; $index < $partlength ; ) {
          if ( !defined $config->{ $parts[$index] } ) {
            die "undefined section $parts[$index]";
          }

          if ( $index == ( $partlength - 1 ) || defined $config->{ $parts[ $index + 1 ] } ) {
            $refmap->{ $parts[$index] } = { config => $config, pattern => $pattern };
            $index++;
          }
          else {
            $refmap->{ $parts[$index] } = { config => $config, pattern => $parts[ $index + 1 ] };
            $index += 2;
          }
        }
      }
      else {
        if ( !defined $config->{$targetSection} ) {
          die "undefined section $targetSection";
        }
        $refmap->{$targetSection} = { config => $config, pattern => $pattern };
      }
    }
    else {

      #in another config, has to be array
      my $refSectionName = $config->{$section}{$mapname_config_ref};
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

        if ( !defined $targetConfig->{$targetSection} ) {
          die "undefined section $targetSection in $mapname_config_ref of $section";
        }

        if ( $index == ( $partlength - 2 ) || ref( $parts[ $index + 2 ] ) eq 'HASH' ) {
          $refmap->{$targetSection} = { config => $targetConfig, pattern => $pattern };
          $index += 2;
        }
        else {
          $refmap->{$targetSection} = { config => $targetConfig, pattern => $parts[ $index + 2 ] };
          $index += 3;
        }
      }
    }

    my %result = ();
    for my $refsec ( keys %{$refmap} ) {
      my $values       = $refmap->{$refsec};
      my $targetConfig = $values->{config};
      my $pattern      = $values->{pattern};

      my %myres = ();
      if ( defined $targetConfig->{$refsec}{class} ) {
        my $myclass = instantiate( $targetConfig->{$refsec}{class} );
        %myres = %{ $myclass->result( $targetConfig, $refsec, $pattern ) };
      }
      else {
        my ( $res, $issource ) = do_get_raw_files( $targetConfig, $refsec, 1 );
        %myres = %{$res};
      }

      my $refcount = keys %myres;
      if ( $refcount > 0 ) {
        @result{ keys %myres } = values %myres;
      }
    }

    my $final = \%result;
    return ( $final, 0 );
  }

  if ($returnself) {
    return ( $config->{$section}, 0 );
  }
  else {
    die "define $mapname or $mapname_ref or $mapname_config_ref for $section";
  }
}

sub get_raw_files {
  my ( $config, $section, $mapname, $pattern ) = @_;
  my ( $result, $issource ) = do_get_raw_files( $config, $section, 0, $mapname, $pattern );
  return ($result);
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
  my ( $pairs, $pairName ) = @_;
  my $groupNames;
  my $ispaired      = 0;
  my $tmpGroupNames = $pairs->{$pairName};
  if ( ref($tmpGroupNames) eq 'HASH' ) {
    $groupNames = $tmpGroupNames->{"groups"};
    $ispaired   = $tmpGroupNames->{"paired"};
  }
  else {
    $groupNames = $tmpGroupNames;
  }
  if ( !defined $ispaired ) {
    $ispaired = 0;
  }
  return ( $ispaired, $groupNames );
}

sub get_pair_groups_names {
  my ( $pairs, $pairName ) = @_;
  my $groupNames;
  my $pairedNames;
  my $tmpGroupNames = $pairs->{$pairName};
  if ( ref($tmpGroupNames) eq 'HASH' ) {
    $groupNames  = $tmpGroupNames->{"groups"};
    $pairedNames = $tmpGroupNames->{"paired"};
  }
  else {
    $groupNames = $tmpGroupNames;
  }
  return ( $pairedNames, $groupNames );
}

1;

