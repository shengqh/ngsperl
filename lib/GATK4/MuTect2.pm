#!/usr/bin/perl
package GATK4::MuTect2;

use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::Fasta;
use CQS::GroupTask;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_mt2";
  $self->{_docker_prefix} = "gatk4_";
  $self->{_use_tmp_folder} = 1;
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $suffix = ".vcf.gz";

  my $germline_resource = get_param_file( $config->{$section}{germline_resource}, "germline_resource", 0 );
  if(defined $germline_resource and ($germline_resource ne "")){
    $option = $option . " --germline-resource $germline_resource ";
  }

  my $panel_of_normals = parse_param_file( $config, $section, "panel_of_normals", 0 );
  #print("panel_of_normals=" . $panel_of_normals . "\n");
  if(defined $panel_of_normals and ($panel_of_normals ne "")){
    $option = $option . " --panel-of-normals $panel_of_normals ";
  }

  #target region
  my $intervals_option = "";
  my $intervals = get_param_file( $config->{$section}{intervals}, "intervals", 0 );
  if ( defined $intervals and ($intervals ne "") ) {
    $intervals_option = "--interval-set-rule INTERSECTION -L $intervals";
    $option = $option . " -L $intervals";
  }

  my $m2_extra_filtering_args = get_option($config, $section, "m2_extra_filtering_args", "");

  my $variants_for_contamination = get_param_file( $config->{$section}{variants_for_contamination}, "variants_for_contamination", 0 );

  my $java_option = get_option( $config, $section, "java_option", "-Xms" . lc($memory) );

  #Sample and Group
  my %raw_files = %{ get_raw_files( $config, $section ) };
  my %group_sample_map;
  if ( has_raw_files( $config, $section, "groups" ) ) {
    %group_sample_map = %{ get_group_sample_map( $config, $section ) };
  }
  else {
    %group_sample_map = %raw_files;
  }
  #print(Dumper(%group_sample_map));

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";
  print $sh "cd $pbs_dir\n";

  for my $group_name ( sort keys %group_sample_map ) {
    my $cur_dir  = create_directory_or_die( $result_dir . "/$group_name" );

    my @sample_files = @{ $group_sample_map{$group_name} };
    my $sampleCount  = scalar(@sample_files);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    my $vcf            = "${group_name}.unfiltered${suffix}";
    my $filtered_vcf   = "${group_name}.filtered${suffix}";
    my $passvcf        = "${group_name}.pass${suffix}";

    print $sh "if [[ ! -s $result_dir/$passvcf ]]; then
  \$MYCMD ./$pbs_name 
fi
    
";

    my $log_desc = $cluster->get_log_description($log);
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $passvcf );

    my $localized_files = [];

    my $normal = undef;
    my $tumor = undef;
    my $sample_param = undef;
    if ( $sampleCount == 1 ) {
      my $tumor_ref = $sample_files[0];

      my $tumor_name;
      if(is_string($tumor_ref)){
        $tumor_name = $group_name;
        $tumor = $tumor_ref;
      }else{
        $tumor_name = $tumor_ref->[0];
        $tumor = $tumor_ref->[1];
      }

      my $local_files = $self->localize_files_in_tmp_folder($pbs, [$tumor], $localized_files, [".bai"]);
      $tumor = $local_files->[0];

      $sample_param = "-I $tumor -tumor $tumor_name";
    }
    elsif ( $sampleCount != 2 ) {
      die "SampleFile should be tumor only or normal,tumor paired.";
    }
    else {
      my $normal_name = $sample_files[0][0];
      $normal      = $sample_files[0][1];
      my $tumor_name = $sample_files[1][0];
      $tumor       = $sample_files[1][1];

      my $local_files = $self->localize_files_in_tmp_folder($pbs, [$normal, $tumor], $localized_files, [".bai"]);
      $normal = $local_files->[0]; 
      $tumor = $local_files->[1];

      $sample_param = "-I $tumor -tumor $tumor_name -I $normal -normal $normal_name";
    }

    print $pbs "$init_command 

mkdir tmp_${group_name}    
";

    if ( $sampleCount == 2 ) {
      print $pbs "
if [ ! -s ${normal}.bai ]; then
  samtools index ${normal}
fi

";
    }
    print $pbs "
if [ ! -s ${tumor}.bai ]; then
  samtools index ${tumor}
fi

echo Mutect2 ...
gatk --java-options \"-Djava.io.tmpdir=`pwd`/tmp_${group_name} $java_option\" Mutect2 $option \\
  -R $faFile \\
  $sample_param \\
  --f1r2-tar-gz f1r2.tar.gz \\
  -O $vcf

m2_exit_code=\$?
if [[ \$m2_exit_code -eq 0 ]]; then
  echo LearnReadOrientationModel ...
  gatk --java-options \"-Djava.io.tmpdir=`pwd`/tmp_${group_name} $java_option\" LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz
";

  my $contamination_table_option="";
  if (defined $variants_for_contamination){
    print $pbs "
  echo GetPileupSummaries tumor ...
  gatk --java-options \"-Djava.io.tmpdir=`pwd`/tmp_${group_name} $java_option\" GetPileupSummaries \\
    -R $faFile $intervals_option \\
    -I $tumor \\
    -V $variants_for_contamination \\
    -L $variants_for_contamination \\
    -O ${group_name}.tumor_pileups.table 
";
    my $normal_option = "";
    if (defined $normal){
      print $pbs "
  echo GetPileupSummaries normal ...
  gatk --java-options \"-Djava.io.tmpdir=`pwd`/tmp_${group_name} $java_option\" GetPileupSummaries \\
    -R $faFile $intervals_option \\
    -I $normal \\
    -V $variants_for_contamination \\
    -L $variants_for_contamination \\
    -O ${group_name}.normal_pileups.table 
";
      $normal_option = "-matched ${group_name}.normal_pileups.table";
    }

    print $pbs "
  echo CalculateContamination ...
  gatk --java-options \"-Djava.io.tmpdir=`pwd`/tmp_${group_name} $java_option\" CalculateContamination \\
    -I ${group_name}.tumor_pileups.table \\
    -O ${group_name}.contamination.table \\
    --tumor-segmentation ${group_name}.segments.table $normal_option
";

    $contamination_table_option = "--contamination-table ${group_name}.contamination.table --tumor-segmentation ${group_name}.segments.table";
  }

  print $pbs "
  echo FilterMutectCalls ...
  gatk --java-options \"-Djava.io.tmpdir=`pwd`/tmp_${group_name} $java_option\" FilterMutectCalls \\
    -V ${vcf} \\
    -R $faFile \\
    -O $filtered_vcf $contamination_table_option $m2_extra_filtering_args \\
    --stats ${vcf}.stats \\
    --ob-priors read-orientation-model.tar.gz \\
    --filtering-stats ${filtered_vcf}.stats

  filter_exit_code=\$?
  if [[ \$filter_exit_code -ne 0 ]]; then
    rm ${filtered_vcf}.*
    touch ${filtered_vcf}.failed
  else
    echo SelectVariants ...
    gatk --java-options \"-Djava.io.tmpdir=`pwd`/tmp_${group_name} $java_option\" SelectVariants \\
      --exclude-filtered \\
      -V $filtered_vcf \\
      -O $passvcf

    select_exit_code=\$?
    if [[ \$select_exit_code -ne 0 ]]; then
      rm ${passvcf}.*
      touch ${passvcf}.failed
    fi
  fi
else
  rm ${vcf}.*
  touch ${vcf}.failed
fi

rm -rf tmp_${group_name} f1r2.tar.gz read-orientation-model.tar.gz

";

    $self->clean_temp_files($pbs, $localized_files);

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $suffix = ".vcf.gz";

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my %group_sample_map;
  if ( has_raw_files( $config, $section, "groups" ) ) {
    %group_sample_map = %{ get_group_sample_map( $config, $section ) };
  }
  else {
    %group_sample_map = %raw_files;
  }

  my $result = {};
  for my $group_name ( sort keys %group_sample_map ) {
    my @result_files = ();
    push( @result_files, "$result_dir/${group_name}/${group_name}.pass${suffix}" );
    push( @result_files, "$result_dir/${group_name}/${group_name}.filtered${suffix}" );
    push( @result_files, "$result_dir/${group_name}/${group_name}.unfiltered${suffix}" );
    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
