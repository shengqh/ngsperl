#!/usr/bin/perl
package GATK4::MuTect2indel;

use strict;
use warnings;
use File::Basename;
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
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $ERC_mode = get_option($config, $section, "ERC_mode", "NONE");
  my $suffix = ($ERC_mode eq "GCVF") ? ".g.vcf.gz" : ".vcf.gz";

  my $germline_resource = get_param_file( $config->{$section}{germline_resource}, "germline_resource", 0 );
  if(defined $germline_resource and ($germline_resource ne "")){
    $option = $option . " --germline-resource $germline_resource ";
  }

  my $panel_of_normals = get_param_file( $config->{$section}{panel_of_normals}, "panel_of_normals", 0 );
  if(defined $panel_of_normals and ($panel_of_normals ne "")){
    $option = $option . " --panel-of-normals $panel_of_normals ";
  }

  #target region
  my $intervals = get_param_file( $config->{$section}{intervals}, "intervals", 0 );
  my $interval_padding = get_option( $config, $section, "interval_padding", 0 );
  my $restrict_intervals = "";
  if ( defined $intervals and ($intervals ne "") ) {
    $option = $option . " -L $intervals";
    if ( defined $interval_padding and $interval_padding != 0 ) {
      $option = $option . " -ip $interval_padding";
    }
  }

  my $java_option = get_option( $config, $section, "java_option", "-Xms" . lc($memory) );

  #Sample and Group
  my %raw_files = %{ get_raw_files( $config, $section ) };
  my %group_sample_map;
  if ( has_raw_files( $config, $section, "groups" ) ) {

    #%group_sample_map = %{ get_raw_files( $config, $section, "groups" ) };
    %group_sample_map = %{ get_group_sample_map( $config, $section ) };
  }
  else {
    %group_sample_map = %raw_files;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";
  print $sh "cd $pbs_dir\n";

  for my $group_name ( sort keys %group_sample_map ) {
    my @sample_files = @{ $group_sample_map{$group_name} };
    my $sampleCount  = scalar(@sample_files);

    my $normal;
    my $tumor;
    my $sample_param;
    if ( $sampleCount == 1 ) {
      $normal      = "";
      $tumor       = $sample_files[0];
      $sample_param = "-I $tumor";
    }
    elsif ( $sampleCount != 2 ) {
      die "SampleFile should be tumor only or normal,tumor paired.";
    }
    else {
      my $normal_name = $sample_files[0][0];
      $normal      = $sample_files[0][1];
      $tumor       = $sample_files[1][1];
      $sample_param = "-I $tumor -I $normal -normal $normal_name";
    }

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $tmp_vcf        = "${group_name}.somatic.tmp${suffix}";
    my $vcf            = "${group_name}.somatic${suffix}";
    my $passvcf        = "${group_name}.somatic.pass${suffix}";
    my $snpPass        = "${group_name}.somatic.snp.pass${suffix}";
    my $indelPass      = "${group_name}.somatic.indel.pass${suffix}";


    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $indelPass );
    print $pbs "$init_command \n";

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

if [ ! -s $vcf ]; then
  echo calling variation ...
  gatk Mutect2 --java-options \"$java_option\" $option \\
    -R $faFile \\
    $sample_param \\
    -ERC $ERC_mode \\
    -O $tmp_vcf

  if [[ -s ${tmp_vcf}.tbi ]]; then
    mv $tmp_vcf $vcf
    mv ${tmp_vcf}.tbi ${vcf}.tbi
  fi
fi 

if [[ -s $vcf && ! -s $passvcf ]]; then
  echo filtering pass ...
  gatk SelectVariants --java-options \"$java_option\" \\
    --exclude-filtered \\
    -V $vcf \\
    -O $passvcf
fi

if [[ -s $passvcf && ! -s $indelPass ]]; then
  echo filtering snp ...
  gatk SelectVariants --java-options \"$java_option\" \\
    -select-type SNP \\
    -V $passvcf \\
    -O $snpPass

  echo filtering indel ...
  gatk SelectVariants --java-options \"$java_option\" \\
    -select-type INDEL \\
    -V $passvcf \\
    -O $indelPass
fi

";

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

  my $ERC_mode = get_option($config, $section, "ERC_mode", "None");
  my $suffix = ($ERC_mode eq "GCVF") ? ".g.vcf.gz" : ".vcf.gz";

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $group_name ( keys %{$groups} ) {
    my @result_files = ();
    push( @result_files, "$result_dir/${group_name}.somatic.pass.indel${suffix}" );
    push( @result_files, "$result_dir/${group_name}.somatic.pass.snp${suffix}" );
    push( @result_files, "$result_dir/${group_name}.somatic.pass${suffix}" );
    push( @result_files, "$result_dir/${group_name}.somatic${suffix}" );
    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
