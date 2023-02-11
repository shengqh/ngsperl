#!/usr/bin/perl
package CQS::WdlSimple;

#this module doesn't handle input json generation, input should be json file

use strict;
use warnings::register;
use File::Basename;
use CQS::PBS;
use CQS::NGSCommon;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Wdl;
use CQS::StringUtils;
use JSON;
use Data::Dumper;

our @ISA = qw(CQS::Wdl);

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

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $cromwell_config_file = get_option_file($config, $section, "cromwell_config_file");
  my $cromwell_jar = get_option_file($config, $section, "cromwell_jar");
 
  my $wdl_file = get_option_file( $config, $section, "wdl_file");

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
    
    my $sample_input_file = $raw_files->{$sample_name};
    if(is_array($sample_input_file)){
      $sample_input_file = $sample_input_file->[0];
    }

    my $input_file = basename($sample_input_file);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $expect_file );

    print $pbs "
cp $sample_input_file $input_file
";

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
      my $input_option_file = get_option_file( $config, $section, "input_option_file" );

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

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

1;
