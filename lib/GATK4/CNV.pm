#!/usr/bin/perl
package GATK4::CNV;

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
  $self->{_suffix} = "_cnv";
  bless $self, $class;
  return $self;
}

#Based on http://gatkforums.broadinstitute.org/gatk/discussion/6791/description-and-examples-of-the-steps-in-the-cnv-case-and-cnv-pon-creation-workflows#
sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

  #parameter files
  my $gatkJar   = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1 );
  my $targetFile   = get_param_file( $config->{$section}{bed_file},   "bed_file",   1 );
  my $gslLibraryFile   = get_param_file( $config->{$section}{gslLibraryFile},   "gslLibraryFile",   1 );
  my $hdfViewFolder   = get_param_file( $config->{$section}{hdfViewFolder},   "hdfViewFolder",   1 );
  
  my $genoemFile   = get_param_file( $config->{$section}{fasta_file},   "fasta_file",   1 );
  my $PONFile   = get_param_file( $config->{$section}{PanelOfNormal},   "PanelOfNormal",   1 );
  
  #make PBS  
  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile   = $sample_files[0];

    my $final_file      = "${sample_name}.coverage.cr.tsv.segfile.tsv.called.tsv";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command );

    my $inputFile    = $sampleFile;
    my $resultPrefix = $sample_name;
  
    #step1: Collect proportional coverage
    my $step1OutputFile=$resultPrefix.".coverage";
    print $pbs "
if [[ -s $inputFile && ! -s $step1OutputFile ]]; then
  echo Step1=`date` 
  java $option -jar $gatkJar CalculateTargetCoverage -I $inputFile -O $step1OutputFile --targets $targetFile -R $genoemFile -transform PCOV --targetInformationColumns FULL -groupBy SAMPLE -keepdups
fi
";
    #step2: Create coverage profile
    my $step2OutputFileTN=$step1OutputFile.".cr.tsv";
    my $step2OutputFileFNO=$step1OutputFile.".fnt.tsv";
    my $step2OutputFileBHO=$step1OutputFile.".beta_hats.tsv";
    my $step2OutputFilePTN=$step1OutputFile.".pre_tangent_normalization_cr.tsv";
        
print $pbs "
if [[ -s $step1OutputFile && ! -s $step2OutputFileTN ]]; then
  echo Step2=`date` 
  export LD_PRELOAD=/usr/local/gsl/latest/x86_64/gcc46/nonet/lib/libgslcblas.so
  java $option -Djava.library.path=$hdfViewFolder -jar $gatkJar NormalizeSomaticReadCounts -I $step1OutputFile -T $targetFile -PON $PONFile -TN $step2OutputFileTN -FNO $step2OutputFileFNO -BHO $step2OutputFileBHO -PTN $step2OutputFilePTN
fi
";

	#Step 3. Segment coverage profile
	my $step3OutputFile=$step2OutputFileTN.".segfile.tsv";
    print $pbs "
if [[ -s $step2OutputFileTN && ! -s $step3OutputFile ]]; then
  echo Step1=`date` 
  java $option -jar $gatkJar PerformSegmentation -TN $step2OutputFileTN -O $step3OutputFile --log2Input
fi
";	
	
	#Step 4. Plot coverage profile
  my $step4OutputFile=$step3OutputFile.".segfile.tsv";
  my $step4ImagePre=$sample_name.".CNV";
    print $pbs "
if [[ -s $step3OutputFile && ! -s $step4OutputFile ]]; then
  echo Step1=`date` 
  java $option -jar $gatkJar PlotSegmentedCopyRatio -TN $step2OutputFileTN -PTN $step2OutputFilePTN -S $step3OutputFile -pre $step4ImagePre -O . --log2Input
fi
";  

	#Step 5. Call segments
	 my $step5OutputFile=$step3OutputFile.".called.tsv";
    print $pbs "
if [[ -s $step3OutputFile && ! -s $step5OutputFile ]]; then
  echo Step1=`date` 
  java $option -jar $gatkJar CallSegments -TN $step2OutputFileTN -S $step3OutputFile -O $step5OutputFile
fi
";
	
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all GATK CNV tasks. \n ";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $final_file      = "${sample_name}.coverage.cr.tsv.segfile.tsv.called.tsv";
    my @result_files = ();
    push( @result_files, "${result_dir}/${final_file}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
