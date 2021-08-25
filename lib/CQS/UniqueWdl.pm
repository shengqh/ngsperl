#!/usr/bin/perl
package CQS::UniqueWdl;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::StringUtils;
use JSON;
use List::MoreUtils qw(uniq);

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_uwdl";
  $self->{_can_use_docker} = 0;
  $self->{_use_tmp_folder} = 0;
  bless $self, $class;
  return $self;
}

sub can_use_docker(){
  my ($self) = @_;
  return(0);
}

sub get_input_value {
  my ($config, $section, $input_hash, $input_key) = @_;
  my $input_name = $input_key;
  my $result = $input_hash->{$input_name};
  
  my $sample_names = $config->{$section}{sample_names};
  my $is_key_ref = substr($input_key, -4) eq "_ref";
  
  if ($is_key_ref){
    die "Key-ref $input_key should point to array!" if ( is_not_array($result) );
    $input_name =~ s/_config_ref$//g;
    $input_name =~ s/_ref$//g;
    $config->{$section}{$input_key} = $result;
    $result = get_raw_files( $config, $section, $input_name );
    
    if (defined $sample_names){
      my $oldresult = $result;
      $result = {};
      for my $sample_name (@$sample_names){
        $result->{$sample_name} = $oldresult->{$sample_name};
      }
    }
    
    delete $config->{$section}{$input_key};
  }

  return ($input_name, $result, $is_key_ref);
}

sub prepare_wdl_parameters {
  my ($self, $config, $section, $task_name, $result_dir, $replace_dics, $replace_values) = @_;
  
  my $input_parameters = get_option($config, $section, "input_parameters");
  die "input_parameters should be hash in section $section!" if(is_not_hash($input_parameters));

  for my $input_key (keys %$input_parameters){
    my ($input_name, $input_value, $is_key_ref) = get_input_value($config, $section, $input_parameters, $input_key);
    
    if ($is_key_ref) {
      $replace_values->{$input_name} = $input_value->{$task_name};
    }else{
      $replace_values->{$input_name} = $input_value;
    }
  }
}

sub prepare_wdl_array {
  my ($self, $config, $section, $task_name, $result_dir, $replace_dics, $replace_values) = @_;
  
  my $input_array = get_option($config, $section, "input_array", {});
  die "$input_array should be hash" if(is_not_hash($input_array));

  for my $input_key (keys %$input_array){
    my ($input_name, $input_value, $is_key_ref) = get_input_value($config, $section, $input_array, $input_key);

    if ((not $is_key_ref) and (is_not_hash($input_value))) {
      die "Key $input_key in input_array of $section should point to HASH!";
    }
    
    my $res_array = [];
    for my $sample_name (sort keys %$input_value){
      my $sample_values = $input_value->{$sample_name};
      die "value for $sample_name in input_list should be array" if ( is_not_array($sample_values) );
      push @$res_array, @$sample_values;
    }
    
    $res_array=[uniq(@$res_array)];
    $replace_values->{$input_name} = $res_array;
  }
}

sub prepare_wdl_list {
  my ($self, $config, $section, $task_name, $result_dir, $replace_dics, $replace_values) = @_;
  
  my $input_list = get_option($config, $section, "input_list", {});
  die "$input_list should be hash" if(is_not_hash($input_list));

  for my $input_key (keys %$input_list){
    my ($input_name, $input_value, $is_key_ref) = get_input_value($config, $section, $input_list, $input_key);

    if ((not $is_key_ref) and (is_not_hash($input_value))) {
      die "Key $input_key in input_list of $section should point to HASH!";
    }
    
    my $list_file = $self->get_file( $result_dir, $task_name, "." . $input_key . ".list" );
    open( my $list, ">$list_file" ) or die "Cannot create $list_file";
    for my $sample_name (keys %$input_value){
      my $sample_values = $input_value->{$sample_name};
      die "value for $sample_name in input_list should be array" if (is_not_array($sample_values));
      
      for my $sample_value (@$sample_values) {
        print $list $sample_value . "\n";
      }
    }
    close($list);
    
    $replace_values->{$input_name} = $list_file;
  }
}

sub prepare_wdl_single {
  my ($self, $config, $section, $task_name, $result_dir, $replace_dics, $replace_values) = @_;
  
  my $input_single = get_option($config, $section, "input_single", {});
  die "$input_single should be hash" if(is_not_hash($input_single));

  for my $input_key (keys %$input_single){
    my ($input_name, $input_value, $is_key_ref) = get_input_value($config, $section, $input_single, $input_key);

    if ((not $is_key_ref) and (is_not_hash($input_value))) {
      die "Key $input_key in input_array of $section should point to HASH!";
    }
    
    if (not defined $input_value->{$task_name}){
      die "Value for key $input_key in $section should be at task level.";
    }
    
    $replace_values->{$input_name} = $input_value->{$task_name}[0];
  }
}

sub prepare_wdl_values {
  my ($self, $config, $section, $task_name, $result_dir) = @_;

  my $replace_dics = {};
  my $replace_values = {};
  
  $self->prepare_wdl_parameters($config, $section, $task_name, $result_dir, $replace_dics, $replace_values);
  
  $self->prepare_wdl_single($config, $section, $task_name, $result_dir, $replace_dics, $replace_values);

  $self->prepare_wdl_list($config, $section, $task_name, $result_dir, $replace_dics, $replace_values);
    
  $self->prepare_wdl_array($config, $section, $task_name, $result_dir, $replace_dics, $replace_values);
  
  return($replace_dics, $replace_values);
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $cromwell_config_file = get_option_file($config, $section, "cromwell_config_file");
  my $cromwell_jar = get_option_file($config, $section, "cromwell_jar");
 
  my $wdl_file = get_option_file( $config, $section, "wdl_file");
  my $input_option_file = get_option_file( $config, $section, "input_option_file" );

  my $input_json_file = get_option_file( $config, $section, "input_json_file" );

  #softlink singularity_image_files to result folder
  my $singularity_image_files = get_raw_files( $config, $section, "singularity_image_files" ); 
  for my $image_name ( sort keys %$singularity_image_files ) {
    my $source_image=$singularity_image_files->{$image_name};
    if( is_array($source_image) ) {
      $source_image=${$source_image}[0];
    }
    my $target_image = $result_dir."/".$image_name;
    if (! -e $target_image) {
      my $simgSoftlinkCommand;
      if (-l $singularity_image_files->{$image_name}) { #softlink, copy
          $simgSoftlinkCommand="cp -P ".$source_image." ".$target_image;
      } else { #file, make softlink
          $simgSoftlinkCommand="ln -s ".$source_image." ".$target_image;
      }
      print($simgSoftlinkCommand."\n");
      system($simgSoftlinkCommand);
    }
  }

  my ($replace_dics, $replace_values) = $self->prepare_wdl_values($config, $section, $task_name, $result_dir);
  
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

  for my $input_key (keys %$replace_values){
    $json_dic->{$input_key} = $replace_values->{$input_key};
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $cur_dir = $result_dir . "/cromwell_finalOutputs";
#  my %raw_files = %{ get_raw_files( $config, $section ) };
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
