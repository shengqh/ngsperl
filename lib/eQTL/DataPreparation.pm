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

  my $plink_prefix_names = get_raw_files( $config, $section, "plink_prefix_names" );
  my $prefix = get_option( $config, $section, "prefix", "" );

  my $fpkm_files = get_raw_files( $config, $section );
  my $genepos_files = get_raw_files( $config, $section, "gene_pos" );
  my $remove_temp_files = get_option( $config, $section, "remove_temp_files", 1 );

  my $find_common_sample_script = dirname(__FILE__) . "/findCommonSample.py";
  if ( !-e $find_common_sample_script ) {
    die "File not found : " . $find_common_sample_script;
  }

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
    my $fpkm_file    = $fpkm_files->{$sample_name}->[0];
    my $genepos_file = $genepos_files->{$sample_name}->[0];

    my $plink_prefix = $plink_prefix_names->{$sample_name}->[0];
    my $fam_file     = $plink_prefix . ".fam";

    my $file_name = $prefix . $sample_name;

    my $common_fam_file = $file_name . "_common.fam";
    my $common_bed_file = $file_name . "_common.bed";

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

if [ ! -s $common_fam_file ]; then
  echo findCommonSample=`date`
  python $find_common_sample_script -f $fam_file -r $fpkm_file -o $common_fam_file 
fi

if [ ! -s $common_bed_file ]; then
  echo extractCommonSample=`date`
  plink -bfile $plink_prefix --keep $common_fam_file --out ${file_name}_common --make-bed 
fi

if [ ! -s related/${file_name}_common_related.genome.genome ]; then
  echo finding related samples ...
  mkdir related
  cd related
  pyGenClean_find_related_samples  --bfile ../${file_name}_common --out ${file_name}_common_related
  cd ..
fi

if [ ! -s ${file_name}_common_filtered.bed ]; then
  echo filtering SNP ...
  plink --bfile ${file_name}_common --out ${file_name}_common_filtered --mind 0.05 --geno 0.05 --maf 0.05 --hwe 0.001 --make-bed
fi

if [ ! -s ${file_name}_common_filtered.traw ]; then
  echo prepareing data for MatrixQTL ...
  plink --bfile ${file_name}_common_filtered --out ${file_name}_common_filtered --recode A-transpose tab
fi

if [ ! -s ${file_name}_common_filtered.tsv ]; then
  echo convering ${file_name}_brca_common_filtered.traw to ${file_name}_common_filtered.tsv
  python $format_snp_data_script -i ${file_name}_common_filtered.traw -o ${file_name}_common_filtered.tsv
fi

cd ..

if [ ! -s ${file_name}_snp.pos ]; then
  echo converting ${file_name}_common_filtered.tsv ${file_name}_snp.pos
  awk '{ if (NR==1){print \"snp\\tchr\\tpos\";}else{print \$2 \"\\t\" \$1 \"\\t\" \$4;} }' temp_dir/${file_name}_common_filtered.tsv > ${file_name}_snp.pos

  echo converting ${file_name}_common_filtered.tsv ${file_name}_snp.genotype
  cut -f2,7- temp_dir/${file_name}_common_filtered.tsv | awk -v OFS='\\t' '{ if (\$1 == \"SNP\") \$1 = \"id\"; print}' > ${file_name}_snp.genotype
fi

if [ ! -s ${file_name}_gene.bed.gz ]; then
  echo prepareGeneExpression ${file_name}_gene.bed.gz
  R --vanilla -f $prepare_gene_script --args $fpkm_file temp_dir/$common_fam_file $genepos_file ${file_name}_gene.bed
  bgzip ${file_name}_gene.bed && tabix -p bed ${file_name}_gene.bed.gz
fi

if [[ -s ${file_name}_gene.bed.gz && ! -s ${file_name}_gene.pos ]]; then
  echo preparing ${file_name}_gene.expression
  zcat ${file_name}_gene.bed.gz | cut -f4- > ${file_name}_gene.expression
  
  echo preparing ${file_name}_gene.pos
  zcat ${file_name}_gene.bed.gz | awk '{ if (NR==1){print \"geneid\\tchr\\ts1\\ts2\";}else{print \$4 \"\\t\" \$1 \"\\t\" \$2 \"\\t\" \$3;} }' > ${file_name}_gene.pos
fi

";

if($remove_temp_files){
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $prefix = get_option( $config, $section, "prefix", "" );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    my $cur_dir  = $result_dir . "/" . $sample_name;

    my $final_name = $prefix . $sample_name;
    push( @result_files, "${cur_dir}/${final_name}_snp.genotype" );
    push( @result_files, "${cur_dir}/${final_name}_snp.pos" );
    push( @result_files, "${cur_dir}/${final_name}_gene.expression" );
    push( @result_files, "${cur_dir}/${final_name}_gene.pos" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
