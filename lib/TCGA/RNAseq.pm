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
  $self->{_name} = "TCGA::RNAseq";
  $self->{_suffix} = "_tcgar";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  
  my $picardbin = get_option($config, $section, "picard_dir");
  my $bedtoolsbin = get_option($config, $section, "bedtools_dir");
  my $tcgabin = get_option($config, $section, "tcga_bin_dir");
  my $ubujar = "${tcgabin}/ubu-1.0.jar";
  my $ubuoption = "-cp $ubujar edu.unc.bioinf.ubu.Ubu";
    
  my $mapslice_version = get_option($config, $section, "mapslice_version", "2");
  my $mapsplicebin;
  my $samtools;
  if($mapslice_version eq "1" || $mapslice_version == 1){
    $mapsplicebin = "${tcgabin}/MapSplice_multi_threads_2.0.1.9/bin";
     $samtools = "${tcgabin}/MapSplice_multi_threads_2.0.1.9/samtools-0.1.9/samtools";
  }else{
    $mapsplicebin = "${tcgabin}/MapSplice_multithreads_12_07/bin";
    $samtools = "${tcgabin}/MapSplice_multithreads_12_07/samtools-0.1.9/samtools";
  }
    
  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    
    if(scalar(@sampleFiles) != 2){
      die "only pair-end data allowed, error sample: " . $sampleName;
    }
    my $sample1 = $sampleFiles[0];
    my $sample2 = $sampleFiles[1];

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );
    
    my $curDir  = create_directory_or_die( $resultDir . "/$sampleName" );
    
    print SH "\$MYCMD ./$pbsName \n";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

if [ ! -d working ]; then
  mkdir working
fi 

#1. Format fastq 1 for Mapsplice
java -Xmx512M $ubuoption fastq-format --phred33to64 --strip --suffix /1 -in $sample1 --out working/prep_1.fastq > working/mapsplice_prep1.log

#2. Format fastq 2 for Mapsplice
java -Xmx512M $ubuoption fastq-format --phred33to64 --strip --suffix /1 -in $sample2 --out working/prep_2.fastq > working/mapsplice_prep2.log

#3. Mapsplice
python $mapsplicebin/mapsplice_multi_thread.py --fusion --all-chromosomes-files ${tcgabin}/hg19_M_rCRS/hg19_M_rCRS.fa --pairend -X 8 -Q fq --chromosome-files-dir ${tcgabin}/hg19_M_rCRS/chromosomes --Bowtieidx ${tcgabin}/hg19_M_rCRS/ebwt/humanchridx_M_rCRS -1 working/prep_1.fastq -2 working/prep_2.fastq -o $sampleName 2> working/mapsplice.log

#4. Add read groups
java -Xmx2G -jar $picardbin/AddOrReplaceReadGroups.jar INPUT=alignments.bam OUTPUT=working/rg_alignments.bam RGSM=$sampleName RGID=$sampleName RGLB=TruSeq RGPL=illumina RGPU=$sampleName VALIDATION_STRINGENCY=SILENT TMP_DIR=./add_rg_tag_tmp > working/add_rg_tag.log 2> working/add_rg_tag.log

#5. Convert back to phred33
java -Xmx512M $ubuoption sam-convert --phred64to33 --in working/rg_alignments.bam -out working/phred33_alignments.bam > working/sam_convert.log 2> working/sam_convert.log

#6. Sort by coordinate
$samtools sort working/phred33_alignments.bam ${sampleName}.bam

#7. Flagstat 
$samtools flagstat ${sampleName}.bam > ${sampleName}.bam.flagstat 
 
#8. Index 
$samtools index ${sampleName}.bam

#9. Sort by chromosome, then read id
perl ${tcgabin}/sort_bam_by_reference_and_name.pl --input ${sampleName}.bam -output working/sorted_by_chr_read.bam --temp-dir . -samtools $samtools > working/sorted_by_chr_read.log 2>working/sorted_by_chr_read.log

#10. Translate to transcriptome coords
java -Xms3G -Xmx3G $ubuoption sam-xlate --bed ${tcgabin}/unc_hg19.bed -in working/sorted_by_chr_read.bam --out working/transcriptome_alignments.bam -order ${tcgabin}/rsem_ref/hg19_M_rCRS_ref.transcripts.fa --xgtags --reverse >working/genome_to_transcriptome.log 2> working/genome_to_transcriptome.log

#11. Filter indels, large inserts, zero mapping quality from transcriptome bam
java -Xmx512M $ubuoption sam-filter --in working/transcriptome_alignments.bam -out working/transcriptome_alignments_filtered.bam --strip-indels --max-insert 10000 --mapq 1 > working/sam_filter.log 2> working/sam_filter.log

#12. RSEM
$tcgabin/rsem-1.1.13/rsem-calculate-expression --gcr-output-file --paired-end --bam --estimate-rspd -p 8 working/transcriptome_alignments_filtered.bam $tcgabin/rsem_ref/hg19_M_rCRS_ref ${sampleName}.rsem > working/rsem.log 2> working/rsem.log

#13. Strip trailing tabs from rsem.isoforms.results
perl ${tcgabin}/strip_trailing_tabs.pl --input ${sampleName}.rsem.isoforms.results --temp working/orig.isoforms.results > working/trim_isoform_tabs.log 2>working/trim_isoform_tabs.log 

#14. Prune isoforms from gene quant file 
mv ${sampleName}.rsem.genes.results orig.genes.results; sed /^uc0/d orig.genes.results > ${sampleName}.rsem.genes.results

#15. Normalize gene quant
perl ${tcgabin}/quartile_norm.pl -c 2 -q 75 -t 1000 -o ${sampleName}.rsem.genes.normalized_results ${sampleName}.rsem.genes.results 

#16. Normalize isoform quant
perl ${tcgabin}/quartile_norm.pl -c 2 -q 75 -t 300 -o ${sampleName}.rsem.isoforms.normalized_results ${sampleName}.rsem.isoforms.results

#17. Junction counts
java -Xmx512M $ubuoption sam-junc --junctions ${sampleName}.splice_junctions.txt --in ${sampleName}.bam --out ${sampleName}.junction_quantification.txt >${sampleName}.junction_quantification.log 2> ${sampleName}.junction_quantification.log

#18. Exon counts
$bedtoolsbin/coverageBed -split -abam ${sampleName}.bam -b ${tcgabin}/composite_exons.bed | perl ${tcgabin}/normalizeBedToolsExonQuant.pl ${tcgabin}/composite_exons.bed> ${sampleName}.bt.exon_quantification.txt 2> ${sampleName}.bt_exon_quantification.log

#19. Cleanup large intermediate output
rm alignments.bam working/phred33_alignments.bam working/rg_alignments.bam working/sorted_by_chr_read.bam working/transcriptome_alignments.bam working/transcriptome_alignments_filtered.bam working/prep_1.fastq working/prep_2.fastq > working/cleanup.log

echo finished=`date`

exit 0 
";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all fastx_trimmer tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $extension = get_option( $config, $section, "extension" );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my @resultFiles = ();

    for my $sampleFile (@sampleFiles) {
      my $trimFile = get_trim_file( $sampleFile, $extension );
      push( @resultFiles, "${resultDir}/${trimFile}" );
    }
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
