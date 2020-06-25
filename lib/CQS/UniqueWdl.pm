#!/usr/bin/perl
package CQS::UniqueWdl;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;
use JSON;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_uwdl";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $cromwell_config_file = get_option_file($config, $section, "cromwell_config_file");
  my $cromwell_jar = get_option_file($config, $section, "cromwell_jar");
 
  my $wdl_file = get_option_file( $config, $section, "wdl_file");
  my $input_option_file = get_option_file( $config, $section, "input_option_file" );

  my $input_json_file = get_option_file( $config, $section, "input_json_file" );

  #softlink singularity_image_files to result folder
  my $singularity_image_files = get_raw_files( $config, $section, "singularity_image_files" ); 
  for my $image_name ( sort keys %$singularity_image_files ) {
    my $simgSoftlinkCommand="cp -P ".$singularity_image_files->{$image_name}[0]." ".$result_dir."/".$image_name;
    print($simgSoftlinkCommand."\n");
  #  print $singularity_image_files->{$image_name}[0]."\n";
    system($simgSoftlinkCommand);
    #`ls -la`;
  }

  my $raw_files = get_raw_files( $config, $section );
  
  my $input_parameters = get_option($config, $section, "input_parameters");
  die "input_parameters should be hash" if(ref($input_parameters) ne 'HASH');

  my $input_list = get_option($config, $section, "input_list", {});
  die "$input_list should be hash" if(ref($input_list) ne 'HASH');

  my $input_single = get_option($config, $section, "input_single", {});
  die "$input_single should be hash" if(ref($input_single) ne 'HASH');

  my $input_array = get_option($config, $section, "input_array", {});
  die "$input_array should be hash" if(ref($input_array) ne 'HASH');

  my $replace_dics = {};
  my $replace_values = {};

  for my $input_key (keys %$input_parameters){
    my $input_value = $input_parameters->{$input_key};
    if (ref($input_value) ne "ARRAY"){
      $replace_values->{$input_key} = $input_value;
      next;
    }
    
    if (substr($input_key, -4) eq "_ref") {
      my $input_name = $input_key;
      $input_name =~ s/_config_ref$//g;
      $input_name =~ s/_ref$//g;
      $config->{$section}{$input_key} = $input_parameters->{$input_key};
      $replace_dics->{$input_name} = get_raw_files( $config, $section, $input_name );
      delete $config->{$section}{$input_key};
    }else{
      $replace_values->{$input_key} = $input_parameters->{$input_key};
    }
  }

  for my $input_key (keys %$input_array){
    my $input_name = $input_key;
    my $input_value = $input_array->{$input_key};
    if (ref($input_value) eq "ARRAY"){
      die "$input_key should include _ref suffix" if (substr($input_key, -4) ne "_ref");
      $input_name =~ s/_config_ref$//g;
      $input_name =~ s/_ref$//g;
      $config->{$section}{$input_key} = $input_array->{$input_key};
      $input_value = get_raw_files( $config, $section, $input_name );
      delete $config->{$section}{$input_key};
    }

    die "$input_key should point to hash" if (ref($input_value) ne "HASH");
    
    my $res_array = [];
    for my $sample_name (keys %$input_value){
      my $sample_values = $input_value->{$sample_name};
      die "value for $sample_name in input_list should be array" if (ref($sample_values) ne "ARRAY");
      push @$res_array, @$sample_values;
    }
    
    $replace_dics->{$input_name} = $res_array;
  }

  for my $input_key (keys %$input_list){
    my $input_name = $input_key;
    my $input_value = $input_list->{$input_key};
    if (ref($input_value) eq "ARRAY"){
      die "$input_key should include _ref suffix" if (substr($input_key, -4) ne "_ref");
      $input_name =~ s/_config_ref$//g;
      $input_name =~ s/_ref$//g;
      $config->{$section}{$input_key} = $input_parameters->{$input_key};
      $input_value = get_raw_files( $config, $section, $input_name );
      delete $config->{$section}{$input_key};
    }

    die "$input_key should point to hash" if (ref($input_value) ne "HASH");
    
    my $list_file = $self->get_file( $result_dir, $task_name, "." . $input_key . ".list" );
    open( my $list, ">$list_file" ) or die "Cannot create $list_file";
    for my $sample_name (keys %$input_value){
      my $sample_values = $input_value->{$sample_name};
      die "value for $sample_name in input_list should be array" if (ref($sample_values) ne "ARRAY");
      
      for my $sample_value (@$sample_values) {
        print $list $sample_value . "\n";
      }
    }
    close($list);
    
    $replace_dics->{$input_name} = $list_file;
  }
  
  for my $input_key (keys %$input_single){
    my $input_value = $input_single->{$input_key};
    
    if ((substr($input_key, -4) ne "_ref") and (ref($input_value) eq ref(""))){
      $replace_values->{$input_key} = $input_value;
      next;
    }
    die "value of $input_key in input_single should be ARRAY" if (ref($input_value) ne "ARRAY");
    
    if (substr($input_key, -4) eq "_ref") {
      my $input_name = $input_key;
      $input_name =~ s/_config_ref$//g;
      $input_name =~ s/_ref$//g;
      $config->{$section}{$input_key} = $input_parameters->{$input_key};
      my $singles = get_raw_files( $config, $section, $input_name );
      delete $config->{$section}{$input_key};
      my $single_file = $singles->{$task_name}[0];
      $replace_values->{$input_key} = $single_file;
    }else{
      $replace_values->{$input_key} = $input_single->{$input_key}[0];
    }
  }
  
  my $json_dic = read_json($input_json_file);
  for my $input_key (keys %$replace_dics){
    if (!defined $json_dic->{$input_key}){
      die "Cannot find " . $input_key . " in json file";
    }
  }

  my $json = JSON->new;
  $json = $json->pretty([1]);

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);
  
  for my $input_key (keys %$replace_dics){
    $json_dic->{$input_key} = $replace_dics->{$input_key};
  }
  
  my $task_input_file = $self->get_file( $result_dir, $task_name, ".inputs.json" );
  open my $fh, ">", $task_input_file;
  print $fh $json->encode($json_dic);
  close $fh;
  
  my $input_file = basename($task_input_file);
  
  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

  print $pbs "
java -Dconfig.file=$cromwell_config_file -jar $cromwell_jar run $wdl_file --inputs $input_file --options $input_option_file
  
";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $cur_dir = $result_dir . "/cromwell_finalOutputs";
  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $output_exts = get_output_ext_list( $config, $section );
  
  my $result = {};
  my @result_files = ();
  for my $output_ext (@$output_exts) {
    if ( $output_ext ne "" ) {
      my $result_file = $task_name . $output_ext;
      push( @result_files, "${cur_dir}/$result_file" );
    }
  }
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}


1;
