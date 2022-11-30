#!/usr/bin/perl
package CQS::Wdl;

use strict;
use warnings::register;
use File::Basename;
use CQS::PBS;
use CQS::NGSCommon;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;
use JSON;
use Data::Dumper;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_wdl";
  $self->{_can_use_docker} = 0;
  #$self->{_forbid_tmp_folder}  = 1;
  bless $self, $class;
  return $self;
}

sub can_use_docker(){
  my ($self) = @_;
  return(0);
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $cromwell_config_file = get_option_file($config, $section, "cromwell_config_file");
  my $cromwell_jar = get_option_file($config, $section, "cromwell_jar");
 
  my $wdl_file = get_option_file( $config, $section, "wdl_file");
  my $input_option_file = get_option_file( $config, $section, "input_option_file" );

  my $input_json_file = get_option_file( $config, $section, "input_json_file" );

  my $sample_name_regex = get_option( $config, $section, "sample_name_regex", "" );

  my $check_output_file_pattern = get_option( $config, $section, "check_output_file_pattern", "" );

  my $use_caper = get_option( $config, $section, "use_caper", 0 );

  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder", 1);

  #softlink singularity_image_files to result folder
  my $singularity_image_files = get_raw_files( $config, $section, "singularity_image_files" ); 
  my $singularity_option = "";
  for my $image_name ( sort keys %$singularity_image_files ) {
    my $source_image=$singularity_image_files->{$image_name};
    if( is_array($source_image) ) {
      $source_image=${$source_image}[0];
    }

    $singularity_option="--singularity $source_image";

    # print $image_name."\n";
    # print $source_image."\n";
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

  my $raw_files = get_raw_files( $config, $section );
  
  my $input_parameters = get_option($config, $section, "input_parameters");
  die "input_parameters should be hash" if is_not_hash($input_parameters);

  #write input to list file
  my $input_list = get_option($config, $section, "input_list", {});
  die "$input_list should be hash" if is_not_hash($input_list);

  my $input_parameters_is_vector = get_option($config, $section, "input_parameters_is_vector", {});
  die "$input_parameters_is_vector should be hash" if is_not_hash($input_parameters_is_vector);

  #replace value with single value/file
  my $input_single = get_option($config, $section, "input_single", {});
  die "$input_single should be hash" if is_not_hash($input_single);

  my $replace_dics = {};
  my $replace_values = {};

  for my $input_key (keys %$input_parameters){
    my $input_value = $input_parameters->{$input_key};
    if ( is_not_array($input_value) ){
      $replace_values->{$input_key} = $input_value;
      next;
    }
    
    if (substr($input_key, -4) eq "_ref") {
      my $input_name = $input_key;
      $input_name =~ s/_config_ref$//g;
      $input_name =~ s/_ref$//g;
      $config->{$section}{$input_key} = $input_parameters->{$input_key};
      $replace_dics->{$input_name} = get_raw_files( $config, $section, $input_name );

      #print($input_name);
      #print(Dumper($replace_dics->{$input_name}));

      delete $config->{$section}{$input_key};
    }else{
      $replace_values->{$input_key} = $input_parameters->{$input_key};
    }
  }

  for my $input_key (keys %$input_list){
    my $input_name = $input_key;
    my $input_value = $input_list->{$input_key};
    if ( is_array($input_value) ){
      die "$input_key should include _ref suffix" if (substr($input_key, -4) ne "_ref");
      $input_name =~ s/_config_ref$//g;
      $input_name =~ s/_ref$//g;
      $config->{$section}{$input_key} = $input_list->{$input_key};
      $input_value = get_raw_files( $config, $section, $input_name );
      delete $config->{$section}{$input_key};
    }

    die "$input_key should point to hash" if is_not_hash($input_value);
    
    my $cur_dic = {};
    for my $sample_name (keys %$input_value){
      my $sample_values = $input_value->{$sample_name};
      die "value for $sample_name in input_list should be array" if is_not_array($sample_values);
      
      my $list_file = $self->get_file( $result_dir, $sample_name, "." . $input_key . ".list" );
      open( my $list, ">$list_file" ) or die "Cannot create $list_file";
      for my $sample_value (@$sample_values) {
        print $list $sample_value . "\n";
      }
      close($list);
      
      $cur_dic->{$sample_name} = [$list_file];
    }
    
    $replace_dics->{$input_name} = $cur_dic;
  }

  for my $input_key (keys %$input_single){
    my $input_value = $input_single->{$input_key};
    
    if ((substr($input_key, -4) ne "_ref") and is_string($input_value)){
      $replace_values->{$input_key} = $input_value;
      next;
    }
    die "value of $input_key in input_single should be ARRAY" if is_not_array($input_value);
    
    if (substr($input_key, -4) eq "_ref") {
      my $input_name = $input_key;
      $input_name =~ s/_config_ref$//g;
      $input_name =~ s/_ref$//g;
      $config->{$section}{$input_key} = $input_single->{$input_key};
      my $singles = get_raw_files( $config, $section, $input_name );
      delete $config->{$section}{$input_key};
      my $single_file = $singles->{$task_name}[0];
      $replace_values->{$input_name} = $single_file;
    }else{
      $replace_values->{$input_key} = $input_single->{$input_key}[0];
    }
  }
  
  my $json_dic = read_json($input_json_file);
  for my $input_key (keys %$replace_dics){
    if (!defined $json_dic->{$input_key}){
  #    die "Cannot find " . $input_key . " in json file";
      warnings::warn("Cannot find " . $input_key . " in json file");
    }
  }
  
  my @json_keys_toSampleNames=();
  for my $replace_key (keys %$replace_values){
    $json_dic->{$replace_key} = $replace_values->{$replace_key};
    if ($replace_values->{$replace_key} eq "SAMPLE_NAME") {
      push @json_keys_toSampleNames,$replace_key;
    }
  }

  my $json = JSON->new;
  $json = $json->pretty([1]);

  my $expected = $self->result($config, $section);

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %$raw_files ) {
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );
    my $cromwell_finalOutputs = $use_caper ? 0 : get_option($config, $section, "cromwell_finalOutputs", 1);
    my $final_dir = $cromwell_finalOutputs ? $cur_dir . "/cromwell_finalOutputs" : $cur_dir;

    my $expect_file = $expected->{$sample_name}[0];

    print $sh "if [[ ! -s $expect_file ]]; then
  \$MYCMD ./$pbs_name
fi
";

    my $log_desc = $cluster->get_log_description($log);
    
    # for my $json_key (keys %$json_dic){
    #   if ($json_dic->{$json_key} =~ /SAMPLE_NAME/){
    #     $json_dic->{$json_key} =~ s/SAMPLE_NAME/$sample_name/g;
    #   }
    # }
    if (scalar(@json_keys_toSampleNames) != 0) {
      my $cur_sample_name = $sample_name;
      if ($sample_name_regex ne ""){
        $cur_sample_name = $1 if ($sample_name =~ /$sample_name_regex/);
      }
      foreach my $json_key_toSampleNames (@json_keys_toSampleNames) {
        $json_dic->{$json_key_toSampleNames} = $cur_sample_name;
      }
    }
    
    for my $input_key (keys %$replace_dics){
      my $input_values = $replace_dics->{$input_key}{$sample_name};
      if (is_array($json_dic->{$input_key})){
        $json_dic->{$input_key} = $input_values;
      }else{
        if($input_parameters_is_vector->{$input_key}){
          $json_dic->{$input_key} = $input_values;
        }else{
          $json_dic->{$input_key} = $input_values->[0];
        }
      }
    }
    
    #print("Before deletion: " . Dumper($json_dic));
    my @keys = keys %$json_dic;
    for my $key (@keys){
      if ($json_dic->{$key} eq ""){
        delete($json_dic->{$key});
        next;
      }

      if(is_array($json_dic->{$key})){
        my $jarray = $json_dic->{$key};
        if (scalar(@$jarray) == 0){
          delete($json_dic->{$key});
          next;
        }

        #print(@$jarray);
      }
    }
    #print("After deletion: " . Dumper($json_dic));

    my $sample_input_file = $self->get_file( $cur_dir, $sample_name, ".inputs.json" );
    open my $fh, ">", $sample_input_file;
    print $fh $json->encode($json_dic);
    close $fh;
    
    my $input_file = basename($sample_input_file);
    
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $expect_file );

    if( $self->{"_use_tmp_folder"}){
      print $pbs "
cp -P \$res_dir/*.simg .
cp $sample_input_file $input_file
";
    }

    if($check_output_file_pattern ne ""){
      print $pbs "
file_count=\$(find . -name $check_output_file_pattern | wc -l) 
if [[ \$file_count -gt 0 ]]; then
  echo \"Warning: $check_output_file_pattern found \$file_count times in $cur_dir !\"
  exit 0
fi
"      
    }

    print $pbs "
if [[ -e $final_dir/${sample_name}.failed ]]; then
  rm $final_dir/${sample_name}.failed
fi
";
    if($use_caper){
      print $pbs "
caper run $wdl_file $option -i $input_file $singularity_option -m $cur_dir/metadata.json
    
";
    }else{
      print $pbs "
java -Dconfig.file=$cromwell_config_file \\
  -jar $cromwell_jar \\
  run $wdl_file $option \\
  --inputs $input_file \\
  --options $input_option_file
";
    }

    print $pbs "
status=\$?
if [[ \$status -ne 0 ]]; then
  touch $final_dir/${sample_name}.failed
fi
";

    if( $self->{"_use_tmp_folder"}){
      print $pbs "
rm *.simg
rm $input_file
";
      if(not $use_caper){
      print $pbs "
rm -rf cromwell-executions
";
      }
    }
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $use_caper = get_option( $config, $section, "use_caper", 0 );
  my $output_to_same_folder = get_option($config, $section, "output_to_same_folder", 1);
  my $cromwell_finalOutputs = $use_caper ? 0 : get_option($config, $section, "cromwell_finalOutputs", 1);
  my $use_filename_in_result = get_option($config, $section, "use_filename_in_result", 0);

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $output_exts = get_output_ext_list( $config, $section );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $cur_dir = $output_to_same_folder ? $result_dir : $result_dir . "/$sample_name";
    my $final_dir = $cromwell_finalOutputs ? $cur_dir . "/cromwell_finalOutputs" : $cur_dir;

    my $sample_file = $raw_files{$sample_name}[0];
    my $sample_prefix = $use_filename_in_result ? change_extension(basename($sample_file), "") : $sample_name;
    my @result_files = ();
    for my $output_ext (@$output_exts) {
      if ( $output_ext ne "" ) {
        if($cromwell_finalOutputs){
          push( @result_files, "${final_dir}/${sample_prefix}${output_ext}" );
        }else{
          push( @result_files, "${final_dir}/$output_ext" );
        }
      }
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}


1;
