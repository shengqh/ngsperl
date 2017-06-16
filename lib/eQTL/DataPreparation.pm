#!/usr/bin/perl
package eQTL::DataPreparation;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_dp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $rnaseq_files = get_raw_files( $config, $section );
  my $plink_prefix_names = get_raw_files( $config, $section, "plink_prefix_names" );

  my $find_common_sample_script = dirname(__FILE__) . "/findCommonSample.py";
  if ( !-e $find_common_sample_script ) {
    die "File not found : " . $find_common_sample_script;
  }

  for my $sample_name ( sort keys %$rnaseq_files ) {
    my $rnaseq_file  = $rnaseq_files->{$sample_name}->[0];
    my $plink_prefix = $plink_prefix_names->{$sample_name}->[0];
    my $fam_file     = $plink_prefix . ".fam";

    my $common_fam_file = $sample_name . "_common.fam";
    my $common_bed_file = $sample_name . "_common.bed";

    my $cur_dir  = create_directory_or_die( $result_dir . "/" . $sample_name );
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

    print $pbs "
if [ ! -s $common_fam_file ]; then
  echo findCommonSample=`date`
  python $find_common_sample_script -f $fam_file -r $rnaseq_file -o $common_fam_file 
fi

if [ ! -s $common_bed_file ]; then
  echo extractCommonSample=`date`
  plink -bfile $plink_prefix --keep $common_fam_file --out ${sample_name}_common --make-bed 
fi

if [ ! -s related/${sample_name}_common_related.genome.genome ]; then
  echo finding related samples ...
  mkdir related
  cd related
  pyGenClean_find_related_samples  --bfile ../${sample_name}_common --out ${sample_name}_common_related
  cd ..
fi

if [ ! -s ${sample_name}_common_filtered.bed ]; then
  echo filtering SNP ...
  plink --bfile ${sample_name}_common --out ${sample_name}_common_filtered --mind 0.05 --geno 0.05 --maf 0.05 --hwe 0.001 --make-bed
fi

if [ ! -s ${sample_name}_common_filtered.traw ]; then
  echo prepareing data for MatrixQTL ...
  plink --bfile ${sample_name}_common_filtered --out ${sample_name}_common_filtered --recode A-transpose tab
fi
";

    #if [ ! -s ${sample_name}_common_filtered.tsv ]; then
    #  echo convering ${name}_brca_common_filtered.traw to ${name}_brca_common_filtered.tsv
    #  python /home/shengq1/program/projects/guoyan/20160826_guoyan_tcga_multi_omics/formatSNPData.py -i ${name}_brca_common_filtered.traw -o ${name}_brca_common_filtered.tsv
    #fi

    #";
    $self->close_pbs( $pbs, $pbs_file );
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();

    my $final_file = $sample_name . ".bed";
    push( @result_files, "${result_dir}/${final_file}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
