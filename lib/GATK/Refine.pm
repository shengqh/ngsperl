#!/usr/bin/perl
package GATK::Refine;

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
	$self->{_suffix} = "_rf";
	bless $self, $class;
	return $self;
}

sub perform {
	my ( $self, $config, $section ) = @_;

	my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

	my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
	my $vcfFiles           = $config->{$section}{indel_vcf_files} or die "Define indel_vcf_files in section $section first.";
	my @vcfFiles           = @{$vcfFiles};
	my $sitesVcfFiles      = $config->{$section}{known_vcf_files} or die "Define known_vcf_files in section $section first.";
	my @sitesVcfFiles      = @{$sitesVcfFiles};
	my $gatk_jar           = get_param_file( $config->{$section}{gatk_jar}, "gatk_jar", 1 );
	my $picard_jar         = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1 );
	my $fixMisencodedQuals = get_option( $config, $section, "fixMisencodedQuals", 0 ) ? "-fixMisencodedQuals" : "";
	my $baq                = get_option( $config, $section, "samtools_baq_calibration", 0 );
	my $slim               = get_option( $config, $section, "samtools_slim", 0 );
	my $remove_duplicate   = get_option( $config, $section, "remove_duplicate", 1 );

	my $bedFile = get_param_file( $config->{$section}{bed_file}, "bed_file", 0 );
	my $interval_padding = get_option( $config, $section, "interval_padding", 0 );
	my $restrict_intervals = "";
	if ( defined $bedFile and $bedFile ne "" ) {
		if ( defined $interval_padding and $interval_padding != 0 ) {
			$restrict_intervals = "-L $bedFile -ip $interval_padding";
		}
		else {
			$restrict_intervals = "-L $bedFile";
		}
	}

	my $indel_vcf = "";
	foreach my $vcf (@vcfFiles) {
		$indel_vcf = $indel_vcf . " -known $vcf";
	}

	my $knownsitesvcf = "";
	my %vcfs          = {};
	foreach my $vcf ( @vcfFiles, @sitesVcfFiles ) {
		if ( $vcfs{$vcf} ) {
			next;
		}
		$knownsitesvcf = $knownsitesvcf . " -knownSites $vcf";
		$vcfs{$vcf} = 1;
	}

	my %raw_files = %{ get_raw_files( $config, $section ) };

	my $sorted = get_option( $config, $section, "sorted", 0 );

	my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
	open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
	print $sh get_run_command($sh_direct) . "\n";

	for my $sample_name ( sort keys %raw_files ) {
		my @sample_files = @{ $raw_files{$sample_name} };
		my $sampleFile   = $sample_files[0];

		my $inputFile     = $sampleFile;
		my $presortedFile = "";
		my $sortCmd       = "";
		if ( !$sorted ) {
			my $presortedPrefix = $sample_name . ".sorted";
			$presortedFile = $presortedPrefix . ".bam";
			$sortCmd       = "samtools sort -@ $thread -m 4G $sampleFile $presortedPrefix";
			$inputFile     = $presortedFile;
		}

		my $rmdupFile       = $inputFile;
		my $rmdupFileIndex  = "";
		my $rmdupFilesToDel = "";
		my $rmdupResultName = "";
		if ($remove_duplicate) {
			$rmdupFile       = $sample_name . ".rmdup.bam";
			$rmdupFileIndex  = change_extension( $rmdupFile, ".bai" );
			$rmdupFilesToDel = "$rmdupFile $rmdupFileIndex ${rmdupFile}.metrics";
			$rmdupResultName = ".rmdup";
		}
		else {

		}

		my $intervalFile = $sample_name . "$rmdupResultName.intervals";

		my $realignedFile = $sample_name . "$rmdupResultName.realigned.bam";
		my $realignedFileIndex = change_extension( $realignedFile, ".bai" );

		my $grpFile = $realignedFile . ".grp";

		my $recalFile = $sample_name . "$rmdupResultName.realigned.recal.bam";
		my $recalFileIndex = change_extension( $recalFile, ".bai" );

		my $final_file = $recalFile;

		my $baqcmd   = "";
		my $slimCmd  = "";
		my $slimFile = "";
		my $rmlist   = "";

		if ($slim) {
			$slimFile   = $sample_name . "$rmdupResultName.realigned.recal.slim.bam";
			$final_file = $slimFile;
			$slimCmd    = "
if [[ -s $recalFile && ! -s $slimFile ]]; then
  echo slim=`date` 
  samtools view -h $recalFile | sed 's/\\tBD\:Z\:[^\\t]*//' | sed 's/\\tPG\:Z\:[^\\t]*//' | sed 's/\\tBI\:Z\:[^\\t]*//' | samtools view -S -b > $slimFile
  samtools index $slimFile
fi  
";
		}

		if ($baq) {
			if ($slim) {
				$final_file = $sample_name . "$rmdupResultName.realigned.recal.slim.baq.bam";
				$baqcmd     = "
if [[ -s $slimFile && ! -s $final_file ]]; then
  echo baq=`date` 
  samtools calmd -Abr $slimFile $faFile > $final_file
  samtools index $final_file
fi      
";
				$rmlist = "$slimFile ${slimFile}.bai";
			}
			else {
				$final_file = $sample_name . "$rmdupResultName.realigned.recal.baq.bam";
				$baqcmd     = "
if [[ -s $recalFile && ! -s $final_file ]]; then
  echo baq=`date` 
  samtools calmd -Abr $recalFile $faFile > $final_file
  samtools index $final_file
fi      
";
			}
		}
		my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
		my $pbs_name = basename($pbs_file);
		my $log      = $self->get_log_filename( $log_dir, $sample_name );

		my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

		print $sh "\$MYCMD ./$pbs_name \n";

		my $log_desc = $cluster->get_log_description($log);

		my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

		if ($remove_duplicate) {
			print $pbs "
if [ ! -s $rmdupFile ]; then
  echo MarkDuplicates=`date` 
  $sortCmd
  java $option -jar $picard_jar MarkDuplicates I=$inputFile O=$rmdupFile ASSUME_SORTED=true REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${rmdupFile}.metrics
fi
";
		}

		print $pbs "
if [[ -s $rmdupFile && ! -s $intervalFile ]]; then
  echo RealignerTargetCreator=`date` 
  java $option -jar $gatk_jar -T RealignerTargetCreator -nt $thread $fixMisencodedQuals -I $rmdupFile -R $faFile $indel_vcf -o $intervalFile $restrict_intervals
fi

if [[ -s $intervalFile && ! -s $realignedFile ]]; then
  echo IndelRealigner=`date` 
  #InDel parameter referenced: http://www.broadinstitute.org/gatk/guide/tagged?tag=local%20realignment
  java $option -Djava.io.tmpdir=tmpdir -jar $gatk_jar -T IndelRealigner $fixMisencodedQuals -I $rmdupFile -R $faFile -targetIntervals $intervalFile $indel_vcf --consensusDeterminationModel USE_READS -LOD 0.4 -o $realignedFile 
fi

if [[ -s $realignedFile && ! -s $grpFile ]]; then
  echo BaseRecalibrator=`date` 
  java $option -jar $gatk_jar -T BaseRecalibrator -nct $thread -rf BadCigar -R $faFile -I $realignedFile $knownsitesvcf -o $grpFile $restrict_intervals
fi

if [[ -s $grpFile && ! -s $recalFile ]]; then
  echo PrintReads=`date`
  java $option -jar $gatk_jar -T PrintReads -nct $thread -rf BadCigar -R $faFile -I $realignedFile -BQSR $grpFile -o $recalFile 
fi

$slimCmd

$baqcmd

if [[ -s $final_file && ! -s ${final_file}.stat ]]; then
  echo flagstat=`date` 
  samtools flagstat $final_file > ${final_file}.stat
  rm $presortedFile $rmdupFilesToDel $intervalFile $realignedFile $realignedFileIndex $grpFile $recalFile $recalFileIndex $rmlist
fi
";
		$self->close_pbs( $pbs, $pbs_file );
	}
	close $sh;

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	print "!!!shell file $shfile created, you can run this shell file to submit all GATK refine tasks.\n";
}
sub result {
	my ( $self, $config, $section, $pattern ) = @_;

	my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

	my %raw_files = %{ get_raw_files( $config, $section ) };
	my $baq              = get_option( $config, $section, "samtools_baq_calibration", 0 );
	my $slim               = get_option( $config, $section, "samtools_slim", 0 );
	my $remove_duplicate = get_option( $config, $section, "remove_duplicate",         1 );

	my $result          = {};
	my $rmdupResultName = "";
	if ($remove_duplicate) {
		$rmdupResultName = ".rmdup";
	}
	for my $sample_name ( keys %raw_files ) {
		my $final_file ="";
		if ($slim) {
			if ($baq) {
				$final_file="$sample_name.$rmdupResultName.realigned.recal.slim.baq.bam";
			} else {
				$final_file="$sample_name.$rmdupResultName.realigned.recal.slim.bam";
			}
		} else {
			if ($baq) {
                $final_file="$sample_name.$rmdupResultName.realigned.recal.baq.bam";
            } else {
                $final_file="$sample_name.$rmdupResultName.realigned.recal.bam";
            }
		}
		my @result_files = ();
		push( @result_files, "${result_dir}/${sample_name}/${final_file}" );
		$result->{$sample_name} = filter_array( \@result_files, $pattern );
	}
	return $result;
}

1;
