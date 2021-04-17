#!/usr/bin/perl
package TCGA::RNAseq;

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
  $self->{_suffix} = "_tcgar";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $picardbin   = get_option( $config, $section, "picard_dir" );
  my $bedtoolsbin = get_option( $config, $section, "bedtools_dir" );
  my $tcgabin     = get_option( $config, $section, "tcga_bin_dir" );

  my $ubuoption = "-cp ${tcgabin}/ubu-1.0.jar:${tcgabin}/jopt-simple-4.6.jar:${tcgabin}/commons-collections-3.2.1.jar:${picardbin}/sam-1.90.jar edu.unc.bioinf.ubu.Ubu";

  my $samtools     = "${tcgabin}/samtools-0.1.19/samtools";
  my $mapsplicebin = "${tcgabin}/MapSplice_multithreads_12_07/bin";

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    if ( scalar(@sample_files) != 2 ) {
      die "only pair-end data allowed, error sample: " . $sample_name;
    }
    my $sample1 = $sample_files[0];
    my $sample2 = $sample_files[1];

    my $sample1name = $sample1;
    my $sample2name = $sample2;
    my $sample1command;
    my $sample2command;
    if ( $sample1 =~ /\.gz$/ ) {
      $sample1name    = "working/data_1.fastq";
      $sample1command = "gunzip -d -c $sample1 > $sample1name
java -Xmx512M $ubuoption fastq-format --phred33to64 --strip --suffix /1 -in $sample1name --out working/prep_1.fastq
rm $sample1name";
    }
    else {
      $sample1command = "java -Xmx512M $ubuoption fastq-format --phred33to64 --strip --suffix /1 -in $sample1name --out working/prep_1.fastq";
    }
    if ( $sample2 =~ /\.gz$/ ) {
      $sample2name    = "working/data_2.fastq";
      $sample2command = "gunzip -d -c $sample2 > $sample2name
java -Xmx512M $ubuoption fastq-format --phred33to64 --strip --suffix /1 -in $sample2name --out working/prep_2.fastq
rm $sample2name";
    }
    else {
      $sample2command = "java -Xmx512M $ubuoption fastq-format --phred33to64 --strip --suffix /1 -in $sample2name --out working/prep_2.fastq";
    }

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

    print $pbs "
if [ ! -d working ]; then
  mkdir working
fi 

echo 1. Format fastq 1 for Mapsplice
$sample1command

echo 2. Format fastq 2 for Mapsplice
$sample2command

echo 3. Mapsplice
python3 $mapsplicebin/mapsplice_multi_thread.py --fusion --all-chromosomes-files ${tcgabin}/hg19_M_rCRS/hg19_M_rCRS.fa --pairend -X 8 -Q fq --chromosome-files-dir ${tcgabin}/hg19_M_rCRS/chromosomes --Bowtieidx ${tcgabin}/hg19_M_rCRS/ebwt/humanchridx_M_rCRS -1 working/prep_1.fastq -2 working/prep_2.fastq -o . 

echo 4. Add read groups
java -Xmx2G -jar $picardbin/AddOrReplaceReadGroups.jar INPUT=alignments.bam OUTPUT=working/rg_alignments.bam RGSM=$sample_name RGID=$sample_name RGLB=TruSeq RGPL=illumina RGPU=$sample_name VALIDATION_STRINGENCY=SILENT TMP_DIR=./add_rg_tag_tmp

echo 5. Convert back to phred33
java -Xmx512M $ubuoption sam-convert --phred64to33 --in working/rg_alignments.bam -out working/phred33_alignments.bam

echo 6. Sort by coordinate
$samtools sort working/phred33_alignments.bam ${sample_name}

echo 7. Flagstat 
$samtools flagstat ${sample_name}.bam > ${sample_name}.bam.flagstat 
 
echo 8. Index 
$samtools index ${sample_name}.bam

echo 9. Sort by chromosome, then read id
perl ${tcgabin}/sort_bam_by_reference_and_name.pl --input ${sample_name}.bam -output working/sorted_by_chr_read.bam --temp-dir . -samtools $samtools

echo 10. Translate to transcriptome coords
java -Xms3G -Xmx3G $ubuoption sam-xlate --bed ${tcgabin}/unc_hg19.bed -in working/sorted_by_chr_read.bam --out working/transcriptome_alignments.bam -order ${tcgabin}/rsem_ref/hg19_M_rCRS_ref.transcripts.fa --xgtags --reverse

echo 11. Filter indels, large inserts, zero mapping quality from transcriptome bam
java -Xmx512M $ubuoption sam-filter --in working/transcriptome_alignments.bam -out working/transcriptome_alignments_filtered.bam --strip-indels --max-insert 10000 --mapq 1

echo 12. RSEM
$tcgabin/rsem_ref/rsem/rsem-calculate-expression --paired-end --bam --estimate-rspd -p 8 working/transcriptome_alignments_filtered.bam $tcgabin/rsem_ref/hg19_M_rCRS_ref ${sample_name}.rsem

echo 13. Strip trailing tabs from rsem.isoforms.results
perl ${tcgabin}/strip_trailing_tabs.pl --input ${sample_name}.rsem.isoforms.results --temp working/orig.isoforms.results

echo 14. Prune isoforms from gene quant file 
mv ${sample_name}.rsem.genes.results orig.genes.results; sed /^uc0/d orig.genes.results > ${sample_name}.rsem.genes.results

echo 15. Normalize gene quant
perl ${tcgabin}/quartile_norm.pl -c 2 -q 75 -t 1000 -o ${sample_name}.rsem.genes.normalized_results ${sample_name}.rsem.genes.results 

echo 16. Normalize isoform quant
perl ${tcgabin}/quartile_norm.pl -c 2 -q 75 -t 300 -o ${sample_name}.rsem.isoforms.normalized_results ${sample_name}.rsem.isoforms.results

echo 17. Junction counts
java -Xmx512M $ubuoption sam-junc --junctions ${sample_name}.splice_junctions.txt --in ${sample_name}.bam --out ${sample_name}.junction_quantification.txt

echo 18. Exon counts
$bedtoolsbin/coverageBed -split -abam ${sample_name}.bam -b ${tcgabin}/composite_exons.bed | perl ${tcgabin}/normalizeBedToolsExonQuant.pl ${tcgabin}/composite_exons.bed

echo 19. Cleanup large intermediate output
#rm alignments.bam working/phred33_alignments.bam working/rg_alignments.bam working/sorted_by_chr_read.bam working/transcriptome_alignments.bam working/transcriptome_alignments_filtered.bam working/prep_1.fastq working/prep_2.fastq

";
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all fastx_trimmer tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $extension = get_option( $config, $section, "extension" );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my @result_files = ();

    for my $sampleFile (@sample_files) {
      my $trimFile = get_trim_file( $sampleFile, $extension );
      push( @result_files, "${result_dir}/${trimFile}" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
