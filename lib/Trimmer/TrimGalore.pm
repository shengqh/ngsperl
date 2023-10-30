#!/usr/bin/perl
package Trimmer::TrimGalore;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;
use HTML::Parser;

our @ISA = qw(CQS::Task);

sub new {
	my ($class) = @_;
	my $self = $class->SUPER::new();
	$self->{_name}   = __PACKAGE__;
	$self->{_suffix} = "_cut";
	bless $self, $class;
	return $self;
}

sub parsePairedSamples {
	my ($samples)    = @_;
	my @sample_files = @{$samples};
	my @result       = ();
	for my $sample (@sample_files) {
		if ( $sample =~ /,/ ) {
			my @files = split( ',', $sample );
			for my $file (@files) {
				push( @result, $file );
			}
		}
		else {
			push( @result, $sample );
		}
	}

	return @result;
}

sub remove_fq_file_extension {
	my ($fileName) = @_;
	$fileName=basename($fileName);
	if ( $fileName =~ /\.gz$/ ) {
		$fileName =~ s/\.gz$//g;
	}
	if ( $fileName =~ /\.fastq$/ ) {
		$fileName =~ s/\.fastq$//g;
	}
	if ( $fileName =~ /\.fq$/ ) {
		$fileName =~ s/\.fq$//g;
	}
	return ($fileName);
}

sub get_final_files {
	my ( $self, $ispairend, $sampleFilesRef, $extension, $fastqextension ) = @_;

	if ($ispairend) {
		my $sampleFileName1 = remove_fq_file_extension( ${$sampleFilesRef}[0] );
		my $sampleFileName2 = remove_fq_file_extension( ${$sampleFilesRef}[1] );
		return (
			$sampleFileName1 . $extension . "_1" . $fastqextension . ".gz",
			$sampleFileName2 . $extension . "_2" . $fastqextension . ".gz"
		);
	}
	else {
		my $sampleFileName1 = remove_fq_file_extension( ${$sampleFilesRef}[0] );
		my $finalName       = $sampleFileName1 . $extension . $fastqextension;
		my $final_file      = "${finalName}.gz";
		my $finalShortFile  = $finalName . ".short.gz";
		my $finalLongFile   = $finalName . ".long.gz";
		return ( $final_file, $finalShortFile, $finalLongFile );
	}
}

sub get_extension {
	my ( $self, $config, $section ) = @_;

	my $extension = $config->{$section}{extension}
	  or die "define ${section}::extension first";
	if ( $extension =~ /\.gz$/ ) {
		$extension =~ s/\.gz$//g;
	}

	my $fastqextension = ".fastq";
	if ( $extension =~ /\.fastq$/ ) {
		$extension =~ s/\.fastq$//g;
	}

	if ( $extension =~ /\.fq$/ ) {
		$fastqextension = ".fq";
		$extension =~ s/\.fq$//g;
	}

	return ( $extension, $fastqextension );
}

sub perform {
	my ( $self, $config, $section ) = @_;

	my (
		$task_name, $path_file,  $pbs_desc, $target_dir, $log_dir,
		$pbs_dir,   $result_dir, $option,   $sh_direct,  $cluster
	) = $self->init_parameter( $config, $section );

	my $ispairend = get_is_paired_end_option( $config, $section, 0 );

	my ( $extension, $fastqextension ) = $self->get_extension( $config, $section );
  my $do_fastqc = get_option($config, $section, "do_fastqc", 0);
  my $fastqc_option = $do_fastqc ? "--fastqc --fastqc_args \"--outdir .\"" : "";

	my %raw_files = %{ get_raw_files( $config, $section ) };

	my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
	open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
	print $sh get_run_command($sh_direct);

	for my $sample_name ( sort keys %raw_files ) {

		my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
		my $pbs_name = basename($pbs_file);
		my $log      = $self->get_log_filename( $log_dir, $sample_name );

		print $sh "\$MYCMD ./$pbs_name \n";

		my $log_desc = $cluster->get_log_description($log);

		my @sample_files = @{ $raw_files{$sample_name} };

		my @final_files =
		  $self->get_final_files( $ispairend, \@sample_files, $extension,
			$fastqextension );
		my $pbs =
		  $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file,
			$result_dir, $final_files[0] );

		if ($ispairend) {
			die "should be pair-end data but not!" if ( scalar(@sample_files) != 2 );

			#pair-end data
			my $read1file = $sample_files[0];
			my $read2file = $sample_files[1];

			print $pbs "
trim_galore $option $fastqc_option --paired --retain_unpaired --output_dir . $read1file $read2file
";
		}
		else {    #single end
			my ( $final_file, $finalShortFile, $finalLongFile ) =
			  $self->get_final_files( $ispairend, $sample_name, $extension,
				$fastqextension );

			if ( scalar(@sample_files) == 1 ) {
				print $pbs "
trim_galore $option $fastqc_option --retain_unpaired --output_dir . $sample_files[0]
";
			}
			else {
				die "should be single-end data but not!";
			}
		}
		$self->close_pbs( $pbs, $pbs_file );
	}
	close $sh;

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	print
"!!!shell file $shfile created, you can run this shell file to submit all "
	  . $self->{_name}
	  . " tasks.\n";

	#`qsub $pbs_file`;
}

sub result {
	my ( $self, $config, $section, $pattern ) = @_;

	my (
		$task_name, $path_file,  $pbs_desc, $target_dir, $log_dir,
		$pbs_dir,   $result_dir, $option,   $sh_direct
	) = $self->init_parameter( $config, $section, 0 );

	my $ispairend = get_option( $config, $section, "pairend", 0 );

	my ( $extension, $fastqextension ) =
	  $self->get_extension( $config, $section );

	my %raw_files = %{ get_raw_files( $config, $section ) };

	my $result = {};
	for my $sample_name ( keys %raw_files ) {
		my @sample_files = @{ $raw_files{$sample_name} };

		my @result_files = ();
		if ($ispairend) {
			my $sampleFileName1 = remove_fq_file_extension( $sample_files[0] );
			my $sampleFileName2 = remove_fq_file_extension( $sample_files[1] );

			my $resultFile1 =
			  $sampleFileName1 . $extension . "_1" . $fastqextension . ".gz";
			push( @result_files,
				"${result_dir}/${resultFile1}" );
			my $resultFile2 =
			  $sampleFileName2 . $extension . "_2" . $fastqextension . ".gz";
			push( @result_files,
				"${result_dir}/${resultFile2}" );
		}
		else {
			my $sampleFileName1 =
			  remove_fq_file_extension( $sample_files[0] );
			my $resultFile1 =
			  $sampleFileName1 . $extension . $fastqextension . ".gz";
			push( @result_files,
				"${result_dir}/${resultFile1}" );
		}

		$result->{$sample_name} = filter_array( \@result_files, $pattern );
	}
	return $result;
}

1;
