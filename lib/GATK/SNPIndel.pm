#!/usr/bin/perl
package GATK::SNPIndel;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_snv";
  bless $self, $class;
  return $self;
}

sub getGroupSampleMap {
  my ( $config, $section ) = @_;

  my $raw_files = get_raw_files( $config, $section );
  my %group_sample_map = ();
  if ( defined $config->{$section}{groups} || defined $config->{$section}{groups_ref} ) {
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $group_name ( sort keys %{$groups} ) {
      my @samples = @{ $groups->{$group_name} };
      my @gfiles  = ();
      my $index   = 0;
      foreach my $sample_name (@samples) {
        my @bam_files = @{ $raw_files->{$sample_name} };
        push( @gfiles, $bam_files[0] );
      }
      $group_sample_map{$group_name} = \@gfiles;
    }
  }
  else {
    %group_sample_map = %{$raw_files};
  }

  return ( \%group_sample_map );
}

#RNASeq: based on https://software.broadinstitute.org/gatk/guide/article?id=3891
sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $faFile   = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1, not $self->using_docker() );

  my $call_option = get_option( $config, $section, "is_rna" ) ? "-stand_emit_conf 20 -stand_call_conf 20" : "--genotyping_mode DISCOVERY  -stand_emit_conf 10 -stand_call_conf 30";
  my $snp_filter =
    get_option( $config, $section, "is_rna" )
    ? "-window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\""
    : "--filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" -filterName \"snp_filter\"";
  my $indel_filter =
    ( get_option( $config, $section, "is_rna" ) ? "-window 35 -cluster 3" : "" ) . " --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" -filterName \"indel_filter\"";

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option ) {
    $java_option = "";
  }

  my %bam_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  #print $sh "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $sample_name ( sort keys %bam_files ) {
    my @sample_files = @{ $bam_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    my $snvOut  = $sample_name . "_snv.vcf";
    my $snvStat = $sample_name . "_snv.stat";

    my $snpOut       = $sample_name . "_snp.vcf";
    my $snpStat      = $sample_name . "_snp.stat";
    my $snpFilterOut = $sample_name . "_snp_filtered.vcf";
    my $snpPass      = $sample_name . "_snp_filtered.pass.vcf";

    my $indelOut       = $sample_name . "_indel.vcf";
    my $indelStat      = $sample_name . "_indel.stat";
    my $indelFilterOut = $sample_name . "_indel_filtered.vcf";
    my $indelPass      = $sample_name . "_indel_filtered.pass.vcf";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

    print $pbs "

if [[ -s $snpOut && ! -s $indelOut ]]; then
  mv $snpOut $snvOut
  if [ -s $snpStat ]; then
    mv $snpStat $snvStat
  fi 
fi

if [ ! -s $snvOut ]; then
  java $java_option -jar $gatk_jar -T HaplotypeCaller -R $faFile -I $bam_file $option -dontUseSoftClippedBases $call_option -nct $thread --out $snvOut
fi

if [ -s $snvOut ]; then
  java $java_option -Xmx${memory} -jar $gatk_jar -T SelectVariants -R $faFile -V $snvOut -selectType SNP -o $snpOut 
  java $java_option -Xmx${memory} -jar $gatk_jar -T VariantFiltration -R $faFile -V $snpOut $snp_filter -o $snpFilterOut 
  cat $snpFilterOut | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $snpPass

  java $java_option -Xmx${memory} -jar $gatk_jar -T SelectVariants -R $faFile -V $snvOut -selectType INDEL -o $indelOut 
  java $java_option -Xmx${memory} -jar $gatk_jar -T VariantFiltration -R $faFile -V $indelOut $indel_filter -o $indelFilterOut 
  cat $indelFilterOut | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $indelPass
fi 

";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all GATK SnpInDel tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  my $result = {};

  my %bam_files = %{ get_raw_files( $config, $section ) };
  for my $sample_name ( sort keys %bam_files ) {
    my $cur_dir      = $result_dir . "/$sample_name";
    my $snpPass      = $sample_name . "_snp_filtered.pass.vcf";
    my $indelPass    = $sample_name . "_indel_filtered.pass.vcf";
    my @result_files = ();
    push( @result_files, "${cur_dir}/${snpPass}" );
    push( @result_files, "${cur_dir}/${indelPass}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result = {};
  if ( $self->{_pbskey} eq "" ) {
    $result->{$task_name} = $self->get_pbs_filename( $pbs_dir, $task_name );
  }
  else {
    my %group_sample_map = %{ getGroupSampleMap( $config, $section ) };

    for my $sample_name ( sort keys %group_sample_map ) {
      my @result_files = ();
      push( @result_files, $self->get_pbs_filename( $pbs_dir, $sample_name ) );
      $result->{$sample_name} = \@result_files;
    }
  }

  return $result;
}

1;
