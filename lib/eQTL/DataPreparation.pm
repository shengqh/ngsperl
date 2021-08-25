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

sub acceptSample {
  my ( $self, $config, $section, $sample_name ) = @_;
  my $sample_set = get_option( $config, $section, "sample_set", [] );
  if ( ( scalar(@$sample_set) > 0 ) & !grep( /^$sample_name$/, @$sample_set ) ) {
    return (0);
  }
  return (1);
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );
  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $prefix = get_option( $config, $section, "prefix", "" );

  my $fpkm_files = get_raw_files( $config, $section );
  my $fpkm_common_files = get_raw_files( $config, $section, "rnaseq_common_samples" );

  my $plink_bed_files    = get_raw_files( $config, $section, "plink_bed_files" );
  my $plink_common_files = get_raw_files( $config, $section, "plink_common_samples" );

  my $genepos_files = get_raw_files( $config, $section, "gene_pos" );
  my $remove_temp_files = get_option( $config, $section, "remove_temp_files", 1 );
  my $pattern           = get_option( $config, $section, "pattern",           1 );

  my $format_snp_data_script = dirname(__FILE__) . "/formatSNPData.py";
  if ( !-e $format_snp_data_script ) {
    die "File not found : " . $format_snp_data_script;
  }

  my $prepare_gene_script = dirname(__FILE__) . "/prepareGeneExpression.r";
  if ( !-e $prepare_gene_script ) {
    die "File not found : " . $prepare_gene_script;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %$fpkm_files ) {
    if ( !$self->acceptSample( $config, $section, $sample_name ) ) {
      next;
    }

    my $fpkm_file          = $fpkm_files->{$sample_name}->[0];
    my $common_rnaseq_file = $fpkm_common_files->{$sample_name}->[0];
    my $genepos_file       = $genepos_files->{$sample_name}->[0];

    my $plink_prefix = $plink_bed_files->{$sample_name}->[0];
    $plink_prefix =~ s/\.[^.]*$//;
    my $common_snp_file = $plink_common_files->{$sample_name}->[0];

    my $file_name = $prefix . $sample_name;

    my $common_prefix   = $file_name . "_common";
    my $common_bed_file = $common_prefix . ".bed";

    my $cur_dir  = create_directory_or_die( $result_dir . "/" . $sample_name );
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

    print $pbs "
if [[ -s ${file_name}_snp.pos && -s ${file_name}_snp.genotype && -s ${file_name}_gene.pos && -s ${file_name}_gene.expression ]]; then
  echo $sample_name has been prepared. Delete $cur_dir/${file_name}_snp.pos to regenerate the result.
  exit 1;
fi

if [ ! -s temp_dir ]; then
  mkdir temp_dir
fi

cd temp_dir    

if [ ! -s $common_bed_file ]; then
  echo extractCommonSample=`date`
  plink -bfile $plink_prefix --keep $common_snp_file --out $common_prefix --make-bed --indiv-sort f $common_snp_file
fi

if [ ! -s ${common_prefix}_filtered.bed ]; then
  echo filtering SNP ...
  plink --bfile $common_prefix --out ${common_prefix}_filtered --mind 0.05 --geno 0.05 --maf 0.05 --hwe 0.001 --make-bed
fi

if [ ! -s ${common_prefix}_filtered.traw ]; then
  echo prepareing data for MatrixQTL ...
  plink --bfile ${common_prefix}_filtered --out ${common_prefix}_filtered --recode A-transpose tab
fi

if [ ! -s ${common_prefix}_filtered.tsv ]; then
  echo convering ${common_prefix}_filtered.traw to ${common_prefix}_filtered.tsv
  python3 $format_snp_data_script -i ${common_prefix}_filtered.traw -o ${common_prefix}_filtered.tsv
fi

cd ..

if [ ! -s ${file_name}_snp.pos ]; then
  echo converting ${common_prefix}_filtered.tsv ${file_name}_snp.pos
  awk '{ if (NR==1){print \"snp\\tchr\\tpos\";}else{print \$2 \"\\t\" \$1 \"\\t\" \$4;} }' temp_dir/${common_prefix}_filtered.tsv > ${file_name}_snp.pos
  cp temp_dir/${common_prefix}_filtered.bim ${file_name}_snp.bim

  echo converting ${common_prefix}_filtered.tsv ${file_name}_snp.genotype
  cut -f2,7- temp_dir/${common_prefix}_filtered.tsv | awk -v OFS='\\t' '{ if (\$1 == \"SNP\") \$1 = \"id\"; print}' > ${file_name}_snp.genotype
fi

if [ ! -s ${file_name}_gene.bed.gz ]; then
  echo prepareGeneExpression ${file_name}_gene.bed.gz
  R --vanilla -f $prepare_gene_script --args $fpkm_file $common_rnaseq_file $genepos_file ${file_name}_gene.bed
  bgzip ${file_name}_gene.bed && tabix -p bed ${file_name}_gene.bed.gz
fi

if [[ -s ${file_name}_gene.bed.gz && ! -s ${file_name}_gene.pos ]]; then
  echo preparing ${file_name}_gene.expression
  zcat ${file_name}_gene.bed.gz | cut -f4- > ${file_name}_gene.expression
  
  echo preparing ${file_name}_gene.pos
  zcat ${file_name}_gene.bed.gz | awk '{ if (NR==1){print \"geneid\\tchr\\ts1\\ts2\";}else{print \$4 \"\\t\" \$1 \"\\t\" \$2 \"\\t\" \$3;} }' > ${file_name}_gene.pos
fi

";

    if ($remove_temp_files) {
      print $pbs "if [[ -s ${file_name}_snp.pos && -s ${file_name}_snp.genotype && -s ${file_name}_gene.pos && -s ${file_name}_gene.expression ]]; then
  rm -rf temp_dir
fi
";
    }
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all fastx_trimmer tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $prefix = get_option( $config, $section, "prefix", "" );
  my $sample_set = get_option( $config, $section, "sample_set", [] );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    if ( !$self->acceptSample( $config, $section, $sample_name ) ) {
      next;
    }

    my @result_files = ();
    my $cur_dir      = $result_dir . "/" . $sample_name;

    my $final_name = $prefix . $sample_name;
    push( @result_files, "${cur_dir}/${final_name}_snp.genotype" );
    push( @result_files, "${cur_dir}/${final_name}_snp.pos" );
    push( @result_files, "${cur_dir}/${final_name}_snp.bim" );
    push( @result_files, "${cur_dir}/${final_name}_gene.expression" );
    push( @result_files, "${cur_dir}/${final_name}_gene.pos" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
