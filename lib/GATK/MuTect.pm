#!/usr/bin/perl
package GATK::MuTect;

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
  $self->{_suffix} = "_mt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

  my $muTect_jar = get_param_file( $config->{$section}{muTect_jar}, "muTect_jar", 1, not $self->using_docker() );
  my $faFile     = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $dbsnpfile = get_param_file( $config->{$section}{dbsnp_file}, "dbsnp_file", 0 );
  my $dbsnpOption = defined $dbsnpfile ? "--dbsnp $dbsnpfile" : "";

  my $cosmicfile = get_param_file( $config->{$section}{cosmic_file}, "cosmic_file", 0 );
  my $cosmicOption = defined $cosmicfile ? "--cosmic $cosmicfile" : "";

  my $java = get_java( $config, $section );

  my $java_option = get_option( $config, $section, "java_option", "" );

  my %group_sample_map = %{ get_group_sample_map( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";
  print $sh "cd $pbs_dir\n";

  for my $group_name ( sort keys %group_sample_map ) {
    my @sample_files = @{ $group_sample_map{$group_name} };
    my $sampleCount  = scalar(@sample_files);
    my $cur_dir      = create_directory_or_die( $result_dir . "/$group_name" );

    if ( $sampleCount != 2 ) {
      die "SampleFile should be normal,tumor paired.";
    }

    my $normal = $sample_files[0][1];
    my $tumor  = $sample_files[1][1];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $out_file = "${group_name}.somatic.out";
    my $vcf      = "${group_name}.somatic.vcf";
    my $passvcf  = "${group_name}.somatic.pass.vcf";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $passvcf, $init_command );
    print $pbs "
if [ ! -s ${normal}.bai ]; then
  samtools index ${normal}
fi

if [ ! -s ${tumor}.bai ]; then
  samtools index ${tumor}
fi

if [ ! -s $vcf ]; then
  $java $java_option -jar $muTect_jar --analysis_type MuTect $option --reference_sequence $faFile $cosmicOption $dbsnpOption --input_file:normal $normal --input_file:tumor $tumor -o $out_file --vcf $vcf --enable_extended_output
fi 

if [[ -s $vcf && ! -s $passvcf ]]; then
  if [ -s $out_file ]; then
    awk '{if (\$1 ~ /^##INFO/) {exit;} else {print;}}' $vcf > $passvcf
    echo \"##INFO=<ID=LOD,Number=1,Type=Float,Description=\\\"Log of (likelihood tumor event is real / likelihood event is sequencing error )\\\">\" >> $passvcf
    awk 'BEGIN{info=0}{if (\$1 ~ /^##INFO/) {info=1;} if(info) {print;}}' $vcf | grep \"^#\" >> $passvcf
    grep -v \"^##\" $out_file | awk -vCOLM=t_lod_fstar 'NR==1 { for (i = 1; i <= NF; i++) { if (\$i == COLM) { cidx = i; } } } NR > 1 {print \";LOD=\"\$cidx}' > ${out_file}.lod
    grep -v \"^#\" $vcf | cut -f1-8 > ${vcf}.first
    paste -d \"\" ${vcf}.first ${out_file}.lod > ${vcf}.first_lod
    grep -v \"^#\" $vcf | cut -f9- > ${vcf}.second
    paste ${vcf}.first_lod ${vcf}.second | grep -v \"^#CHROM\" | grep -v REJECT >> $passvcf
    rm ${out_file}.lod ${vcf}.first ${vcf}.first_lod ${vcf}.second
    grep -v \"^#\" $passvcf | cut -f1 | uniq -c | awk '{print \$2\"\\t\"\$1}' > ${passvcf}.chromosome
    
  else
    grep -v REJECT $vcf > $passvcf
  fi
fi

";

    $self->close_pbs( $pbs, $pbs_file );

    print "$pbs_file created \n";
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $group_name ( keys %{$groups} ) {
    my @result_files = ();
    my $cur_dir      = $result_dir . "/$group_name";
    push( @result_files, "$cur_dir/${group_name}.somatic.pass.vcf" );
    push( @result_files, "$cur_dir/${group_name}.somatic.vcf" );
    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
