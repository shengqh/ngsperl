#!/usr/bin/perl
package Pipeline::TargetWGS;

use strict;
use warnings;
use File::Basename;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Pipeline::WdlPipeline;
use Pipeline::WGS;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performTargetWGS) ] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  getValue($def, "target_genes");

  #check files
  getValue($def, "ref_fasta");
  getValue($def, "dbsnp");
  getValue($def, "known_indels_sites_VCFs");

  #start from FASTQ
  initDefaultValue( $def, "perform_preprocessing", 1 );

  #let's start from bam file
  #initDefaultValue( $def, "perform_preprocessing", 0 );
  #initDefaultValue( $def, "bam_file_section", "files" );

  initDefaultValue( $def, "gatk_prefix", "gatk4_" );

  initDefaultValue( $def, "perform_replace_read_group", 0);
  initDefaultValue( $def, "replace_read_group_by_gatk4", 0);
  
  initDefaultValue( $def, "perform_mark_duplicates", 1);

  initDefaultValue( $def, "perform_gvcf_to_genotype", 1);

  initDefaultValue( $def, "perform_filter_and_merge", 1);

  initDefaultValue( $def, "callvariants_vqsr_mode", 1 );
  
  initDefaultValue( $def, "max_thread", 8 );
  initDefaultValue( $def, "subdir",     0 );

  initDefaultValue( $def, "sra_to_fastq", 0 );

  initDefaultValue( $def, "aligner", "bwa" );
  if ( $def->{aligner} eq "bwa" ) {
    initDefaultValue( $def, "bwa_option", "-K 100000000 -v 3" );
  }

  initDefaultValue( $def, "mark_duplicates_use_tmp_folder", 0 );

  if(defined $def->{"ref_fasta_dict"} && (! defined $def->{"chromosome_names"})){
    my $dictFile = getValue($def, "ref_fasta_dict");
    my $primary_chromosome_only = getValue($def, "primary_chromosome_only", 1);
    my $chromosomes = readChromosomeFromDictFile($dictFile, $primary_chromosome_only);
    initDefaultValue( $def, "chromosome_names" , join(",", @$chromosomes));
  }

  initDefaultValue( $def, "filter_variants_by_allele_frequency",            1 );
  initDefaultValue( $def, "filter_variants_by_allele_frequency_percentage", 0.9 );
  initDefaultValue( $def, "filter_variants_by_allele_frequency_maf",        0.3 );
  initDefaultValue( $def, "filter_variants_fq_equal_1",                     0 );

  initDefaultValue( $def, "perform_cnv_gatk4_cohort", 1 );
  initDefaultValue( $def, "gatk4_cnv_by_scatter",     1 );
  initDefaultValue( $def, "gatk4_cnv_scatter_count",  100 );
  initDefaultValue( $def, "perform_cnv_gatk4_somatic", 0 );

  initDefaultValue( $def, "recalibration_annotation_values", ["QD", 
  "MQRankSum",
  "ReadPosRankSum",
  "FS",
  "MQ",
  "SOR",
  "DP"]);

  initDefaultValue( $def, "snp_recalibration_tranche_values", ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ]);
  initDefaultValue( $def, "snp_recalibration_annotation_values", ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"]);
  initDefaultValue( $def, "indel_recalibration_tranche_values", ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"]);
  initDefaultValue( $def, "indel_recalibration_annotation_values",  ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"]);
  initDefaultValue( $def, "snp_filter_level", 99.7 );
  initDefaultValue( $def, "indel_filter_level", 99.7 );
  initDefaultValue( $def, "SNP_VQSR_downsampleFactor", 10 );

  initDefaultValue( $def, "perform_muTect",       0 );
  initDefaultValue( $def, "perform_muTect2indel", 0 );

  initDefaultValue( $def, "perform_annovar",      1 );
  initDefaultValue( $def, "maximum_freq_values",  "0.05,0.01" );

  initDefaultValue( $def, "perform_cnv",          1 );
  initDefaultValue( $def, "perform_vep",          0 );
  initDefaultValue( $def, "perform_IBS",          0 );

  if ( $def->{perform_muTect} || $def->{perform_muTect2indel} ) {
    if ( defined $def->{mills} ) {
      initDefaultValue( $def, "indel_realignment", 1 );
      initDefaultValue( $def, "indel_vcf_files",   $def->{mills} );
    }
    else {
      initDefaultValue( $def, "indel_realignment", 0 );
    }
  }

  initDefaultValue( $def, "remove_duplicate", 1 );
  initDefaultValue( $def, "perform_multiqc", 0 );

  my $default_onco_options = {
    "picture_width" => "0",
    "picture_height" => "0",
    "sampleNamePattern" => ".",
    "BACKGROUND_color" => "lightgray",
    "MISSENSE_color" => "darkgreen",
    "MISSENSE_height" => 0.75,
    "TRUNC_color" => "red",
    "TRUNC_height" => 0.5, 
    "DUP_color" => "brown", 
    "DUP_height" => 0.25,
    "DEL_color" => "blue",
    "DEL_height" => 0.25,
    "MANUAL_order" => 0,
  };

  if (defined $def->{onco_options}) {
    $def->{onco_options} = merge_hash_left_precedent($def->{onco_options}, $default_onco_options);
  }else{
    $def->{onco_options} = $default_onco_options;
  }

  my $perform_split_fastq = getValue($def, "perform_split_fastq", 0);
  if($perform_split_fastq) {#[ 0, "by_dynamic", "by_file", "by_scatter"],
    if ($perform_split_fastq eq "by_dynamic"){ 
      initDefaultValue($def, "split_fastq_min_file_size_gb", 10); #only the file with file size larger than this number would be splitted
      initDefaultValue($def, "split_fastq_trunk_file_size_gb", 5); #the splitted file will be smaller than this file size
    }elsif($perform_split_fastq eq "by_scatter"){ 
      initDefaultValue($def, "aligner_scatter_count", 10); #split data into equal number of small files
    }elsif($perform_split_fastq eq "by_file"){ 
      #nothing need to be set
    }else{
      die 'wrong perform_split_fastq value, it should be one of [ 0, "by_dynamic", "by_file", "by_scatter"]';
    }
  }

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [@$individual, @$summary];
  $individual = $tasks;
  $summary = $tasks;

  my $ref_fasta = getValue($def, "ref_fasta");

  my $gene_bed = do_add_gene_locus($config, $def, $tasks, $target_dir, "target_gene_locus", $def->{target_genes});

  my $gatk_prefix = getValue($def, "gatk_prefix");
  my $gatk_index_snv = "SNV_index";
  my $callvariants_vqsr_mode = getValue($def, "callvariants_vqsr_mode");

  my $splitFastq = "bwa_00_splitFastq";
  add_split_fastq_dynamic($config, $def, $tasks, $target_dir, $splitFastq, $source_ref);

  my $bwa_option = getValue( $def, "bwa_option" );
  my $bwa_index = getValue( $def, "bwa_fasta" );
  my $rg_name_regex = "(.+)_ITER_";

  my $bwa_01_alignment = "bwa_01_alignment";
  add_BWA_WGS($config, $def, $tasks, $target_dir, $bwa_01_alignment, $splitFastq, $rg_name_regex);

  my $bwa_02_markduplicates = "bwa_02_markduplicates";

  my $mem = getValue($def, "mark_duplicates_mem", "40gb");
  #replace gb with empty and minus 2
  my $java_memory_size = $mem;
  $java_memory_size =~ s/gb//g;
  $java_memory_size = $java_memory_size - 2;

  $config->{$bwa_02_markduplicates} = {
    class                 => "CQS::ProgramWrapperManyToOne",
    perform               => 1,
    target_dir            => "$target_dir/$bwa_02_markduplicates",
    option                => "
echo MarkDuplicates=`date`
gatk --java-options \"-Dsamjdk.compression_level=2 -Xms${java_memory_size}g\" \\
  MarkDuplicates \\
  --INPUT __INPUT__ \\
  --REFERENCE_SEQUENCE $ref_fasta \\
  --METRICS_FILE __NAME__.duplicates_metrics.txt \\
  --VALIDATION_STRINGENCY SILENT \\
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
  --ASSUME_SORT_ORDER coordinate \\
  --CLEAR_DT false \\
  --ADD_PG_TAG_TO_READS false \\
  --OUTPUT tmp.__OUTPUT__

status=\$?
if [[ \$status -ne 0 ]]; then
  touch __OUTPUT__.failed
  rm -f __OUTPUT__.succeed tmp.__OUTPUT__
else
  rm -f __OUTPUT__.failed
  touch __OUTPUT__.succeed
  mv tmp.__OUTPUT__ __OUTPUT__

  echo samtools_index=`date`
  samtools index __OUTPUT__
fi
",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $bwa_01_alignment,
    source_arg            => "",
    source_join_delimiter => " \\\n  --INPUT ",
    output_to_same_folder => 1,
    output_arg            => "--OUTPUT",
    output_file_prefix    => ".duplicates_marked.cram",
    output_file_ext       => ".duplicates_marked.cram",
    use_tmp_folder => getValue($def, "mark_duplicates_use_tmp_folder", 0),
    sh_direct             => 0,
    docker_prefix        => "gatk4_",
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => getValue($def, "mark_duplicates_walltime", "48"),
      "mem"       => getValue($def, "mark_duplicates_mem", "40gb")
    },
  };

  push(@$tasks, $bwa_02_markduplicates);

  my $bam_task = $bwa_02_markduplicates;

  my $target_gene_cram = "target_gene_cram";
  $config->{$target_gene_cram} = {
    class => "CQS::ProgramWrapperOneToOne",
    perform => 1,
    target_dir => $target_dir . "/$target_gene_cram",
    interpretor => "",
    program => "",
    check_program => 0,
    option => "
echo target_gene_cram=`date`
samtools view --cram \\
  --reference $ref_fasta \\
  --target-file __parameterFile1__ \\
  --output __OUTPUT__ \\
  __FILE__ 

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __OUTPUT__.succeed
  rm -f __OUTPUT__.failed

  echo samtools_index=`date`
  samtools index __OUTPUT__
else
  touch __OUTPUT__.failed
  rm -f __OUTPUT__.succeed __OUTPUT__
fi
",
    source_arg => "",
    source_ref => $bam_task,
    source_join_delimiter => " ",
    parameterFile1_ref => $gene_bed,
    output_arg => "--output",
    output_file_prefix => ".genes.cram",
    output_file_ext => ".genes.cram",
    no_output => 1,
    sh_direct => 0,
    pbs => {
      "nodes" => "1:ppn=1",
      "walltime" => "10",
      "mem" => "10gb"
    },
  };

  $bam_task = ["target_gene_cram"];
  push @$tasks, $target_gene_cram;

  my $dbsnp = getValue($def, "dbsnp");
  my $known_indels_sites_VCFs = getValue($def, "known_indels_sites_VCFs");
  my $all_know_vcfs = [ $dbsnp, @$known_indels_sites_VCFs ];
  my $known_indels_sites_VCFs_option = "--known-sites " . join(" \\\n  --known-sites ", @$all_know_vcfs);

  my $gatk4_BaseRecalibrator = $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_BaseRecalibrator";
  $config->{$gatk4_BaseRecalibrator} = {
    class => "CQS::ProgramWrapperOneToOne",
    perform => 1,
    target_dir => $target_dir . "/$gatk4_BaseRecalibrator",
    interpretor => "",
    program => "",
    check_program => 0,
    docker_prefix => "gatk4_",
    option => "
gatk --java-options \"-Xmx20G\" \\
  BaseRecalibrator \\
  -R $ref_fasta \\
  -I __FILE__ \\
  --use-original-qualities \\
  -L __parameterFile1__ \\
  $known_indels_sites_VCFs_option \\
  -O tmp.__OUTPUT__ \\

status=\$?
if [[ \$status -eq 0 ]]; then
  rm -f __OUTPUT__.failed
  touch __OUTPUT__.succeed
  mv tmp.__OUTPUT__ __OUTPUT__
else
  touch __OUTPUT__.failed
  rm -f __OUTPUT__.succeed tmp.__OUTPUT__
fi
",
    source_arg => "",
    source_ref => $bam_task,
    source_join_delimiter => " ",
    parameterFile1_ref => $gene_bed,
    output_arg => "",
    output_file_prefix => ".recal_data.csv",
    output_file_ext => ".recal_data.csv",
    no_output => 1,
    sh_direct => 0,
    pbs => {
      "nodes" => "1:ppn=1",
      "walltime" => "10",
      "mem" => "20gb"
    },
  };   
  push(@$tasks, $gatk4_BaseRecalibrator); 

  my $gatk4_02_ApplyBQSR = $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_ApplyBQSR";
  $config->{$gatk4_02_ApplyBQSR} = {
    class => "CQS::ProgramWrapperOneToOne",
    perform => 1,
    target_dir => $target_dir . "/$gatk4_02_ApplyBQSR",
    interpretor => "",
    program => "",
    check_program => 0,
    docker_prefix => "gatk4_",
    option => "
gatk --java-options \"-Xmx20G\" \\
  ApplyBQSR \\
  --add-output-sam-program-record \\
  --create-output-bam-index false \\
  -R $ref_fasta \\
  -I __FILE__ \\
  --use-original-qualities \\
  -bqsr __FILE2__ \\
  -L __parameterFile1__ \\
  -O tmp.__OUTPUT__

status=\$?
if [[ \$status -ne 0 ]]; then
  touch __OUTPUT__.failed
  rm -f __OUTPUT__.succeed tmp.__OUTPUT__
else
  touch __OUTPUT__.succeed
  rm -f __OUTPUT__.failed
  mv tmp.__OUTPUT__ __OUTPUT__

  echo samtools_index=`date`
  samtools index __OUTPUT__
fi
",
    source_arg => "",
    source_ref => $bam_task,
    source_join_delimiter => " ",
    parameterSampleFile2_ref => [ $gatk4_BaseRecalibrator, ".recal_data.csv" ],
    parameterFile1_ref => $gene_bed,
    output_arg => "",
    output_file_prefix => ".recalibrated.cram",
    output_file_ext => ".recalibrated.cram",
    no_output => 1,
    sh_direct => 0,
    pbs => {
      "nodes" => "1:ppn=1",
      "walltime" => "10",
      "mem" => "20gb"
    },
  };   
  push(@$tasks, $gatk4_02_ApplyBQSR); 

  my $blacklist_file = getValue($def, "blacklist_file");
  my $gatk4_03_HaplotypeCaller = $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_HaplotypeCaller";
  $config->{$gatk4_03_HaplotypeCaller} = {
    class => "CQS::ProgramWrapperOneToOne",
    perform => 1,
    target_dir => $target_dir . "/$gatk4_03_HaplotypeCaller",
    interpretor => "",
    program => "",
    check_program => 0,
    docker_prefix => "gatk4_",
    option => "
gatk --java-options \"-Xmx40G\" \\
  HaplotypeCaller \\
  -contamination 0 \\
  -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \\
  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \\
  --soft-clip-low-quality-ends true \\
  --dont-use-soft-clipped-bases true \\
  -ERC GVCF \\
  -L __parameterFile1__ \\
  -XL $blacklist_file \\
  --native-pair-hmm-threads 8 \\
  -R $ref_fasta \\
  -I __FILE__ \\
  -O tmp.__OUTPUT__

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __OUTPUT__.succeed
  rm -f __OUTPUT__.failed
  mv tmp.__OUTPUT__ __OUTPUT__
  mv tmp.__OUTPUT__.tbi __OUTPUT__.tbi
else
  touch __OUTPUT__.failed
  rm -f __OUTPUT__.succeed tmp.__OUTPUT__ tmp.__OUTPUT__.tbi
fi
",
    source_arg => "",
    source_ref => $gatk4_02_ApplyBQSR,
    source_join_delimiter => " ",
    parameterFile1_ref => $gene_bed,
    output_arg => "",
    output_file_prefix => ".g.vcf.gz",
    output_file_ext => ".g.vcf.gz",
    no_output => 1,
    sh_direct => 0,
    pbs => {
      "nodes" => "1:ppn=1",
      "walltime" => "10",
      "mem" => "40gb"
    },
  };   
  push(@$tasks, $gatk4_03_HaplotypeCaller); 

  my $gatk4_04_GenomicsDBImport = $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_GenomicsDBImport";
  $config->{$gatk4_04_GenomicsDBImport} = {
    class => "CQS::ProgramWrapper",
    perform => 1,
    target_dir => $target_dir . "/$gatk4_04_GenomicsDBImport",
    interpretor => "",
    program => "",
    check_program => 0,
    docker_prefix => "gatk4_",
    option => "
rm -rf __NAME__

gatk --java-options \"-Xmx8g -Xms8g\"  \\
  GenomicsDBImport \\
  --genomicsdb-workspace-path __NAME__ \\
  --batch-size 50 \\
  -L __parameterFile1__ \\
  --sample-name-map __FILE__ \\
  --reader-threads 5 \\
  --merge-input-intervals \\
  --consolidate

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.succeed
  rm -f __NAME__.failed
else
  touch __NAME__.failed
  rm -rf __NAME__.succeed __NAME__
fi
",
    source_ref => $gatk4_03_HaplotypeCaller,
    source_fileFirst => 0,
    parameterFile1_ref => $gene_bed,
    output_arg => "",
    output_file_ext => " ",
    no_output => 1,
    sh_direct => 0,
    pbs => {
      "nodes" => "1:ppn=6",
      "walltime" => "10",
      "mem" => "40gb"
    },
  };   
  push(@$tasks, $gatk4_04_GenomicsDBImport); 

  my $gatk4_05_GenotypeGVCFs = $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_GenotypeGVCFs";
  $config->{$gatk4_05_GenotypeGVCFs} = {
    class => "CQS::ProgramWrapperOneToOne",
    perform => 1,
    target_dir => $target_dir . "/$gatk4_05_GenotypeGVCFs",
    interpretor => "",
    program => "",
    check_program => 0,
    docker_prefix => "gatk4_",
    option => "
gatk --java-options \"-Xms8g -Xmx40g\" \\
  GenotypeGVCFs \\
  -R $ref_fasta \\
  -D $dbsnp \\
  -G StandardAnnotation -G AS_StandardAnnotation \\
  --only-output-calls-starting-in-intervals \\
  -V gendb://__FILE__ \\
  -L __parameterFile1__ \\
  --merge-input-intervals \\
  -O tmp.__OUTPUT__

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __OUTPUT__.succeed
  rm -f __OUTPUT__.failed
  mv tmp.__OUTPUT__ __OUTPUT__
  mv tmp.__OUTPUT__.tbi __OUTPUT__.tbi
else
  touch __OUTPUT__.failed
  rm -f __OUTPUT__.succeed tmp.__OUTPUT__ tmp.__OUTPUT__.tbi
fi
",
    source_ref => $gatk4_04_GenomicsDBImport,
    parameterFile1_ref => $gene_bed,
    output_arg => "",
    output_file_prefix => ".g.vcf.gz",
    output_file_ext => ".g.vcf.gz",
    no_output => 1,
    sh_direct => 0,
    pbs => {
      "nodes" => "1:ppn=1",
      "walltime" => "10",
      "mem" => "40gb"
    },
  };   
  push(@$tasks, $gatk4_05_GenotypeGVCFs); 

  my $gatk4_06_VariantFilterHard = $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_VariantFilterHard";
  $config->{$gatk4_06_VariantFilterHard} = {
    class => "CQS::ProgramWrapperOneToOne",
    perform => 1,
    target_dir => $target_dir . "/$gatk4_06_VariantFilterHard",
    interpretor => "",
    program => "",
    check_program => 0,
    docker_prefix => "gatk4_",
    option => "
echo SelectVariants_SNP=`date`
rm -f __NAME__.snps.vcf.gz.failed __NAME__.snps.vcf.gz.succeed
gatk SelectVariants \\
  -V __FILE__ \\
  -select-type SNP \\
  -O tmp.__NAME__.snps.vcf.gz

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.snps.vcf.gz.succeed
  rm -f __NAME__.snps.vcf.gz.failed

  mv tmp.__NAME__.snps.vcf.gz __NAME__.snps.vcf.gz
  mv tmp.__NAME__.snps.vcf.gz.tbi __NAME__.snps.vcf.gz.tbi
else
  touch __NAME__.snps.vcf.gz.failed
  rm -f __NAME__.snps.vcf.gz.succeed tmp.__NAME__.snps.vcf.gz tmp.__NAME__.snps.vcf.gz.tbi
fi

echo VariantFiltration_SNP=`date`
rm -f __NAME__.snps_filtered.succeed __NAME__.snps_filtered.failed
if [[ -e __NAME__.snps.vcf.gz.succeed ]]; then
  gatk VariantFiltration \\
    -V __NAME__.snps.vcf.gz \\
    -filter \"QD < 2.0\" --filter-name \"QD2\" \\
    -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \\
    -filter \"SOR > 3.0\" --filter-name \"SOR3\" \\
    -filter \"FS > 60.0\" --filter-name \"FS60\" \\
    -filter \"MQ < 40.0\" --filter-name \"MQ40\" \\
    -filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \\
    -filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\"  \\
    -O tmp.__NAME__.snps_filtered.vcf.gz

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch __NAME__.snps_filtered.succeed
    rm -f __NAME__.snps_filtered.failed __NAME__.snps.vcf.gz __NAME__.snps.vcf.gz.tbi __NAME__.snps.vcf.gz.succeed

    mv tmp.__NAME__.snps_filtered.vcf.gz __NAME__.snps_filtered.vcf.gz
    mv tmp.__NAME__.snps_filtered.vcf.gz.tbi __NAME__.snps_filtered.vcf.gz.tbi
  else
    touch __NAME__.snps_filtered.failed
    rm -f __NAME__.snps_filtered.succeed tmp.__NAME__.snps_filtered.vcf.gz tmp.__NAME__.snps_filtered.vcf.gz.tbi
  fi
fi

echo SelectVariants_INDEL=`date`
rm -f __NAME__.indel.vcf.gz.failed __NAME__.indel.vcf.gz.succeed
gatk SelectVariants \\
  -V __FILE__ \\
  -select-type INDEL \\
  -select-type MIXED \\
  -O tmp.__NAME__.indel.vcf.gz

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.indel.vcf.gz.succeed
  rm -f __NAME__.indel.vcf.gz.failed

  mv tmp.__NAME__.indel.vcf.gz __NAME__.indel.vcf.gz
  mv tmp.__NAME__.indel.vcf.gz.tbi __NAME__.indel.vcf.gz.tbi
else
  touch __NAME__.indel.vcf.gz.failed
  rm -f __NAME__.indel.vcf.gz.succeed tmp.__NAME__.indel.vcf.gz tmp.__NAME__.indel.vcf.gz.tbi
fi

echo VariantFiltration_INDEL=`date`
rm -f __NAME__.indels_filtered.succeed __NAME__.indels_filtered.failed
if [[ -e __NAME__.indel.vcf.gz.succeed ]]; then
  gatk VariantFiltration \\
    -V __NAME__.indel.vcf.gz \\
    -filter \"QD < 2.0\" --filter-name \"QD2\" \\
    -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \\
    -filter \"FS > 200.0\" --filter-name \"FS200\" \\
    -filter \"ReadPosRankSum < -20.0\" --filter-name \"ReadPosRankSum-20\"   \\
    -O tmp.__NAME__.indels_filtered.vcf.gz

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch __NAME__.indels_filtered.succeed 
    rm -f __NAME__.indels_filtered.failed __NAME__.indel.vcf.gz __NAME__.indel.vcf.gz.tbi __NAME__.indel.vcf.gz.succeed

    mv tmp.__NAME__.indels_filtered.vcf.gz __NAME__.indels_filtered.vcf.gz
    mv tmp.__NAME__.indels_filtered.vcf.gz.tbi __NAME__.indels_filtered.vcf.gz.tbi
  else
    touch __NAME__.indels_filtered.failed
    rm -f __NAME__.indels_filtered.succeed tmp.__NAME__.indels_filtered.vcf.gz tmp.__NAME__.indels_filtered.vcf.gz.tbi
  fi
fi

echo MergeVcfs=`date`
rm -f __NAME__.filtered.vcf.gz.failed __NAME__.filtered.vcf.gz.succeed
if [[ -e __NAME__.snps_filtered.succeed && -e __NAME__.indels_filtered.succeed ]]; then
  gatk MergeVcfs \\
    -I __NAME__.snps_filtered.vcf.gz \\
    -I __NAME__.indels_filtered.vcf.gz \\
    -O tmp.__NAME__.filtered.vcf.gz

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch __NAME__.filtered.vcf.gz.succeed
    rm -f __NAME__.filtered.vcf.gz.failed __NAME__.snps_filtered.vcf.gz __NAME__.snps_filtered.vcf.gz.tbi __NAME__.indels_filtered.vcf.gz __NAME__.indels_filtered.vcf.gz.tbi __NAME__.snps_filtered.succeed __NAME__.indels_filtered.succeed

    mv tmp.__NAME__.filtered.vcf.gz __NAME__.filtered.vcf.gz
    mv tmp.__NAME__.filtered.vcf.gz.tbi __NAME__.filtered.vcf.gz.tbi
  else
    touch __NAME__.filtered.vcf.gz.failed
    rm -f __NAME__.filtered.vcf.gz.succeed tmp.__NAME__.filtered.vcf.gz tmp.__NAME__.filtered.vcf.gz.tbi
  fi
fi

echo SelectVariants_PASS=`date`
rm -f __OUTPUT__.failed __OUTPUT__.succeed
if [[ -e __NAME__.filtered.vcf.gz.succeed ]]; then
  gatk SelectVariants \\
    --exclude-filtered \\
    -V __NAME__.filtered.vcf.gz \\
    -O tmp.__OUTPUT__

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch __OUTPUT__.succeed
    rm -f __OUTPUT__.failed __NAME__.filtered.vcf.gz __NAME__.filtered.vcf.gz.tbi __NAME__.filtered.vcf.gz.succeed

    mv tmp.__OUTPUT__ __OUTPUT__
    mv tmp.__OUTPUT__.tbi __OUTPUT__.tbi
  else
    touch __OUTPUT__.failed
    rm -f __OUTPUT__.succeed tmp.__OUTPUT__ tmp.__OUTPUT__.tbi
  fi
fi
",
    source_ref => $gatk4_05_GenotypeGVCFs,
    output_arg => "",
    output_file_prefix => ".hardfilter.pass.vcf.gz",
    output_file_ext => ".hardfilter.pass.vcf.gz",
    no_output => 1,
    sh_direct => 0,
    pbs => {
      "nodes" => "1:ppn=1",
      "walltime" => "10",
      "mem" => "20gb"
    },
  };   
  push(@$tasks, $gatk4_06_VariantFilterHard); 

  my $gatk4_07_LeftTrim = $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_LeftTrim";
  my $lefttrim_script = dirname(__FILE__) . "/../GATK4/fixLeftTrimDeletion.py";
  $config->{$gatk4_07_LeftTrim} = {
    class => "CQS::ProgramWrapperOneToOne",
    perform => 1,
    target_dir => $target_dir . "/$gatk4_07_LeftTrim",
    interpretor => "",
    program => "",
    check_program => 0,
    option => "
echo bcftools_split=`date`
bcftools norm -m- \\
  -o tmp.__NAME__.split.vcf \\
  __FILE__

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.split.vcf.succeed
  rm -f __NAME__.split.vcf.failed

  mv tmp.__NAME__.split.vcf __NAME__.split.vcf
else
  touch __NAME__.split.vcf.failed
  rm -f __NAME__.split.vcf.succeed tmp.__NAME__.split.vcf
fi

rm -f __NAME__.norm.vcf.failed __NAME__.norm.vcf.succeed
if [[ -e __NAME__.split.vcf.succeed ]]; then
  echo bcftools_norm=`date`
  bcftools norm -f $ref_fasta \\
    -o tmp.__NAME__.norm.vcf \\
    __NAME__.split.vcf

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch __NAME__.norm.vcf.succeed
    rm -f __NAME__.norm.vcf.failed __NAME__.split.vcf __NAME__.split.vcf.succeed

    mv tmp.__NAME__.norm.vcf __NAME__.norm.vcf
  else
    touch __NAME__.norm.vcf.failed
    rm -f __NAME__.norm.vcf.succeed
  fi
fi

if [[ -e __NAME__.norm.vcf.succeed ]]; then
  echo noSpanDeletion=`date`
  python3 $lefttrim_script \\
    -i __NAME__.norm.vcf \\
    -o tmp.__NAME__.hardfilter.pass.norm.nospan.vcf

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch __OUTPUT__.succeed
    rm -f __OUTPUT__.failed __NAME__.norm.vcf __NAME__.norm.vcf.succeed

    mv tmp.__NAME__.hardfilter.pass.norm.nospan.vcf __NAME__.hardfilter.pass.norm.nospan.vcf
    bgzip __NAME__.hardfilter.pass.norm.nospan.vcf
    tabix -p vcf __NAME__.hardfilter.pass.norm.nospan.vcf.gz
  else
    touch __NAME__.norm.vcf.failed
    rm -f __NAME__.norm.vcf.succeed tmp.__NAME__.hardfilter.pass.norm.nospan.vcf
  fi
fi
",
    source_ref => $gatk4_06_VariantFilterHard,
    output_arg => "",
    output_file_prefix => ".hardfilter.pass.norm.nospan.vcf.gz",
    output_file_ext => ".hardfilter.pass.norm.nospan.vcf.gz",
    no_output => 1,
    sh_direct => 0,
    pbs => {
      "nodes" => "1:ppn=1",
      "walltime" => "10",
      "mem" => "20gb"
    },
  };   
  push(@$tasks, $gatk4_07_LeftTrim); 

#  For target genes, we cannot do VQSR mode. it would failed due to very limited number of known variants in target genes for modelling
#   my $gatk4_06_HardFilterAndMakeSitesOnlyVcf = "gatk4_06_HardFilterAndMakeSitesOnlyVcf";
#   $config->{$gatk4_06_HardFilterAndMakeSitesOnlyVcf} = {
#     class => "CQS::ProgramWrapperOneToOne",
#     perform => 1,
#     target_dir => $target_dir . "/$gatk4_06_HardFilterAndMakeSitesOnlyVcf",
#     interpretor => "",
#     program => "",
#     check_program => 0,
#     docker_prefix => "gatk4_",
#     option => "
# echo VariantFiltration ...
# gatk --java-options \"-Xms10g -Xmx18g\" \\
#   VariantFiltration \\
#   --filter-expression \"ExcessHet > 54.69\" \\
#   --filter-name ExcessHet \\
#   -V __FILE__ \\
#   -O tmp.__OUTPUT__

# status=\$?
# if [[ \$status -eq 0 ]]; then
#   touch __OUTPUT__.succeed
#   rm -f __OUTPUT__.failed

#   mv tmp.__OUTPUT__ __OUTPUT__
#   mv tmp.__OUTPUT__.tbi __OUTPUT__.tbi

#   echo MakeSitesOnlyVcf ...
#   gatk --java-options \"-Xms10g -Xmx18g\" \\
#     MakeSitesOnlyVcf \\
#     -I __OUTPUT__ \\
#     -O tmp.__NAME__.variant_filtered.sites_only.vcf.gz

#   status=\$?
#   if [[ \$status -eq 0 ]]; then
#     touch __NAME__.variant_filtered.sites_only.vcf.gz.succeed
#     rm -f __NAME__.variant_filtered.sites_only.vcf.gz.failed
#     mv tmp.__NAME__.variant_filtered.sites_only.vcf.gz __NAME__.variant_filtered.sites_only.vcf.gz
#     mv tmp.__NAME__.variant_filtered.sites_only.vcf.gz.tbi __NAME__.variant_filtered.sites_only.vcf.gz.tbi
#   else
#     touch __NAME__.variant_filtered.sites_only.vcf.gz.failed
#     rm -f __NAME__.variant_filtered.sites_only.vcf.gz.succeed tmp.__NAME__.variant_filtered.sites_only.vcf.gz tmp.__NAME__.variant_filtered
# .sites_only.vcf.gz.tbi
#   fi
# else
#   touch __OUTPUT__.failed
#   rm -f __OUTPUT__.succeed tmp.__OUTPUT__ tmp.__OUTPUT__.tbi
# fi

# ",
#     source_ref => $gatk4_05_GenotypeGVCFs,
#     output_arg => "",
#     output_file_prefix => ".variant_filtered.vcf.gz",
#     output_file_ext => ".variant_filtered.sites_only.vcf.gz",
#     output_other_ext => ".variant_filtered.vcf.gz",
#     no_output => 1,
#     sh_direct => 0,
#     pbs => {
#       "nodes" => "1:ppn=1",
#       "walltime" => "10",
#       "mem" => "20gb"
#     },
#   };   
#   push(@$tasks, $gatk4_06_HardFilterAndMakeSitesOnlyVcf); 

#   my $mills_resource_vcf = getValue($def, "mills");
#   my $dbsnp_resource_vcf = getValue($def, "dbsnp");
#   my $axiomPoly_resource_vcf = getValue($def, "axiomPoly");

#   my $gatk4_07_IndelsVariantRecalibrator = "gatk4_07_IndelsVariantRecalibrator";
#   $config->{$gatk4_07_IndelsVariantRecalibrator} = {
#     class => "CQS::ProgramWrapperOneToOne",
#     perform => 1,
#     target_dir => $target_dir . "/$gatk4_07_IndelsVariantRecalibrator",
#     interpretor => "",
#     program => "",
#     check_program => 0,
#     docker_prefix => "gatk4_",
#     option => "
# gatk --java-options \"-Xms24g -Xms25g\" \\
#   VariantRecalibrator \\
#   -AS \\
#   -V __FILE__ \\
#   --tranches-file __NAME__.indels.tranches \\
#   --trust-all-polymorphic \\
#   -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \\
#   -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \\
#   -mode INDEL \\
#   --max-gaussians 4 \\
#   --resource:mills,known=false,training=true,truth=true,prior=12 $mills_resource_vcf \\
#   --resource:axiomPoly,known=false,training=true,truth=false,prior=10 $axiomPoly_resource_vcf \\
#   --resource:dbsnp,known=true,training=false,truth=false,prior=2 $dbsnp_resource_vcf \\
#   -O tmp.__OUTPUT__

# status=\$?
# if [[ \$status -eq 0 ]]; then
#   touch __OUTPUT__.succeed
#   rm -f __OUTPUT__.failed
#   mv tmp.__OUTPUT__ __OUTPUT__
#   mv tmp.__OUTPUT__.tbi __OUTPUT__.tbi
# else
#   touch __OUTPUT__.failed
#   rm -f __OUTPUT__.succeed tmp.__OUTPUT__ tmp.__OUTPUT__.tbi __NAME__.indels.tranches
# fi
# ",
#     source_ref => [ $gatk4_06_HardFilterAndMakeSitesOnlyVcf, ".variant_filtered.sites_only.vcf.gz"],
#     output_arg => "",
#     output_file_prefix => ".indels.recal.vcf.gz",
#     output_file_ext => ".indels.recal.vcf.gz",
#     output_other_ext => ".indels.tranches",
#     no_output => 1,
#     sh_direct => 0,
#     pbs => {
#       "nodes" => "1:ppn=1",
#       "walltime" => "10",
#       "mem" => "30gb"
#     },
#   };   
#   push(@$tasks, $gatk4_07_IndelsVariantRecalibrator); 

#   my $hapmap_resource_vcf = getValue($def, "hapmap");
#   my $omni_resource_vcf = getValue($def, "omni");
#   my $one_thousand_genomes_resource_vcf = getValue($def, "g1000");

#   my $gatk4_08_SNPsVariantRecalibratorClassic = "gatk4_08_SNPsVariantRecalibratorClassic";
#   $config->{$gatk4_08_SNPsVariantRecalibratorClassic} = {
#     class => "CQS::ProgramWrapperOneToOne",
#     perform => 1,
#     target_dir => $target_dir . "/$gatk4_08_SNPsVariantRecalibratorClassic",
#     interpretor => "",
#     program => "",
#     check_program => 0,
#     docker_prefix => "gatk4_",
#     option => "
# gatk --java-options \"-Xmx24g -Xms24g\" \\
#   VariantRecalibrator \\
#   -AS \\
#   -V __FILE__ \\
#   --tranches-file __NAME__.snps.tranches \\
#   --trust-all-polymorphic \\
#   -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \\
#   -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \\
#   -mode SNP \\
#   --max-gaussians 6 \\
#   --resource:hapmap,known=false,training=true,truth=true,prior=15 $hapmap_resource_vcf \\
#   --resource:omni,known=false,training=true,truth=true,prior=12 $omni_resource_vcf \\
#   --resource:1000G,known=false,training=true,truth=false,prior=10 $one_thousand_genomes_resource_vcf \\
#   --resource:dbsnp,known=true,training=false,truth=false,prior=7 $dbsnp_resource_vcf \\
#   -O tmp.__OUTPUT__

# status=\$?
# if [[ \$status -eq 0 ]]; then
#   touch __OUTPUT__.succeed
#   rm -f __OUTPUT__.failed
#   mv tmp.__OUTPUT__ __OUTPUT__
#   mv tmp.__OUTPUT__.tbi __OUTPUT__.tbi
# else
#   touch __OUTPUT__.failed
#   rm -f __OUTPUT__.succeed tmp.__OUTPUT__ tmp.__OUTPUT__.tbi __NAME__.snps.tranches
# fi
# ",
#     source_ref => [ $gatk4_06_HardFilterAndMakeSitesOnlyVcf, ".variant_filtered.sites_only.vcf.gz"],
#     output_arg => "",
#     output_file_prefix => ".snps.recal.vcf.gz",
#     output_file_ext => ".snps.recal.vcf.gz",
#     output_other_ext => ".snps.tranches",
#     no_output => 1,
#     sh_direct => 0,
#     pbs => {
#       "nodes" => "1:ppn=1",
#       "walltime" => "10",
#       "mem" => "30gb"
#     },
#   };   
#   push(@$tasks, $gatk4_08_SNPsVariantRecalibratorClassic); 

#   my $gatk4_09_ApplyRecalibration = "gatk4_09_ApplyRecalibration";
#   $config->{$gatk4_09_ApplyRecalibration} = {
#     class => "CQS::ProgramWrapperOneToOne",
#     perform => 1,
#     target_dir => $target_dir . "/$gatk4_09_ApplyRecalibration",
#     interpretor => "",
#     program => "",
#     check_program => 0,
#     docker_prefix => "gatk4_",
#     option => "
# gatk --java-options \"-Xms5000m -Xmx6500m\" \\
#   ApplyVQSR \\
#   -AS \\
#   -V __FILE__ \\
#   --recal-file __parameterFile1__ \\
#   --use-allele-specific-annotations \\
#   --tranches-file __parameterFile2__ \\
#   --truth-sensitivity-filter-level 99.7 \\
#   --create-output-variant-index true \\
#   -mode INDEL \\
#   -O __NAME__.indel.recalibrated.vcf.gz

# status=\$?
# if [[ \$status -eq 0 ]]; then
#   touch __NAME__.indel.recalibrated.vcf.gz.succeed
#   rm -f __NAME__.indel.recalibrated.vcf.gz.failed

#   gatk --java-options \"-Xms5000m -Xmx6500m\" \\
#     ApplyVQSR \\
#     -V __NAME__.indel.recalibrated.vcf.gz \\
#     --recal-file __parameterFile3__ \\
#     --use-allele-specific-annotations \\
#     --tranches-file __parameterFile4__ \\
#     --truth-sensitivity-filter-level 99.7 \\
#     --create-output-variant-index true \\
#     -mode SNP \\
#     -O tmp.__OUTPUT__

#   status=\$?
#   if [[ \$status -eq 0 ]]; then
#     touch __OUTPUT__.succeed
#     rm -f __OUTPUT__.failed __NAME__.indel.recalibrated.vcf.gz __NAME__.indel.recalibrated.vcf.gz.tbi __NAME__.indel.recalibrated.vcf.gz.succeed
#     mv tmp.__OUTPUT__ __OUTPUT__
#     mv tmp.__OUTPUT__.tbi __OUTPUT__.tbi
#   else
#     touch __OUTPUT__.failed
#     rm -f __OUTPUT__.succeed tmp.__OUTPUT__ tmp.__OUTPUT__.tbi
#   fi
# else
#   touch __NAME__.indel.recalibrated.vcf.gz.failed
#   rm -f __NAME__.indel.recalibrated.vcf.gz.succeed __NAME__.indel.recalibrated.vcf.gz __NAME__.indel.recalibrated.vcf.gz.tbi
# fi
# ",
#     source_ref => [ $gatk4_06_HardFilterAndMakeSitesOnlyVcf, ".variant_filtered.vcf.gz"],
#     parameterFile1_ref => [ $gatk4_07_IndelsVariantRecalibrator, ".indels.recal.vcf.gz" ],
#     parameterFile2_ref => [ $gatk4_07_IndelsVariantRecalibrator, ".indels.tranches" ],
#     parameterFile3_ref => [ $gatk4_08_SNPsVariantRecalibratorClassic, ".snps.recal.vcf.gz" ],
#     parameterFile4_ref => [ $gatk4_08_SNPsVariantRecalibratorClassic, ".snps.tranches" ],
#     output_arg => "",
#     output_file_prefix => ".filtered.vcf.gz",
#     output_file_ext => ".filtered.vcf.gz",
#     no_output => 1,
#     sh_direct => 0,
#     pbs => {
#       "nodes" => "1:ppn=1",
#       "walltime" => "10",
#       "mem" => "30gb"
#     },
#   };   
#   push(@$tasks, $gatk4_09_ApplyRecalibration); 

#   my $ref_fasta_dict = getValue($def, "ref_fasta_dict");
#   my $gatk4_10_CollectVariantCallingMetrics = "gatk4_10_CollectVariantCallingMetrics";
#   $config->{$gatk4_10_CollectVariantCallingMetrics} = {
#     class => "CQS::ProgramWrapperOneToOne",
#     perform => 1,
#     target_dir => $target_dir . "/$gatk4_10_CollectVariantCallingMetrics",
#     interpretor => "",
#     program => "",
#     check_program => 0,
#     docker_prefix => "gatk4_",
#     option => "
# gatk --java-options \"-Xms6000m -Xmx7000m\" \\
#   CollectVariantCallingMetrics \\
#   --INPUT __FILE__ \\
#   --DBSNP $dbsnp_resource_vcf \\
#   --SEQUENCE_DICTIONARY $ref_fasta_dict \\
#   --THREAD_COUNT 8 \\
#   --TARGET_INTERVALS __parameterFile1__ \\
#   --OUTPUT tmp.__NAME__

# status=\$?
# if [[ \$status -eq 0 ]]; then
#   touch __NAME__.succeed
#   rm -f __NAME__.failed
#   mv tmp.__NAME__.variant_calling_summary_metrics __NAME__.variant_calling_summary_metrics
#   mv tmp.__NAME__.variant_calling_detail_metrics __NAME__.variant_calling_detail_metrics
# else
#   touch __NAME__.failed
#   rm -f __NAME__.succeed tmp.__NAME__.variant_calling_detail_metrics tmp.__NAME__.variant_calling_summary_metrics
# fi
# ",
#     source_ref => $gatk4_09_ApplyRecalibration,
#     parameterFile1_ref => $gene_bed,
#     output_arg => "",
#     output_file_ext => ".variant_calling_summary_metrics",
#     output_other_ext => ".variant_calling_detail_metrics",
#     no_output => 1,
#     sh_direct => 0,
#     pbs => {
#       "nodes" => "1:ppn=8",
#       "walltime" => "10",
#       "mem" => "8gb"
#     },
#   };   
#   push(@$tasks, $gatk4_10_CollectVariantCallingMetrics); 

  my $annovar_name = addAnnovar($config, $def, $summary, $target_dir, $gatk4_07_LeftTrim, '.gz$', $gatk_prefix, $def, $gatk_index_snv );

  if ( $def->{annovar_param} =~ /exac/ || $def->{annovar_param} =~ /1000g/ || $def->{annovar_param} =~ /gnomad/ ) {
    my $annovar_filter_name = addAnnovarFilter( $config, $def, $summary, $target_dir, $annovar_name, $gatk_prefix, $def, $gatk_index_snv);
    my $mafreport = addAnnovarMafReport($config, $def, $summary, $target_dir, $annovar_filter_name, $gatk_prefix, $def, $gatk_index_snv);
  }

  $config->{sequencetask} = {
    class => "CQS::SequenceTaskSlurmSlim",
    perform => 1,
    target_dir            => "${target_dir}/sequencetask",
    option                => "",
    source                => {
      processing => $tasks,
    },
    sh_direct             => 1,
    pbs                   => {
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  };

  return($config);
};

sub performTargetWGS {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

1;
