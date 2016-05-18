#!/usr/bin/perl
package CQS::RunProgram;

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
use File::Spec;

our @ISA = qw(CQS::Task);

sub new {
	my ($class) = @_;
	my $self = $class->SUPER::new();
	$self->{_name}   = __PACKAGE__;
	$self->{_suffix} = "_perl";
	bless $self, $class;
	return $self;
}

sub perform {
	my ( $self, $config, $section ) = @_;

	my (
		$task_name, $path_file,  $pbs_desc, $target_dir, $log_dir,
		$pbs_dir,   $result_dir, $option,   $sh_direct,  $cluster
	) = get_parameter( $config, $section );

	my $runProgram = get_option( $config, $section, "runProgram" );
	my $is_absolute = File::Spec->file_name_is_absolute($runProgram);
	if ( !$is_absolute ) {
		$runProgram = dirname(__FILE__) . "/$runProgram";
	}
	if ( !( -e $runProgram ) ) {
		die("runProgram $runProgram defined but not exists!");
	}

	my $output_file_label =
	  get_option( $config, $section, "output_file_label", "" );
	my $output_file = get_option( $config, $section, "output_file",     "" );
	my $output_ext  = get_option( $config, $section, "output_file_ext", "" );

	my %raw_files = %{ get_raw_files( $config, $section, "source1" ) };

	my %parameterFiles2 = ();
	if ( has_raw_files( $config, $section, "source2" ) ) {
		%parameterFiles2 = %{ get_raw_files( $config, $section, "source2" ) };
	}
	my %parameterFiles3 = ();
	if ( has_raw_files( $config, $section, "source3" ) ) {
		%parameterFiles3 = %{ get_raw_files( $config, $section, "source3" ) };
	}
	my %parameterFiles4 = ();
	if ( has_raw_files( $config, $section, "source4" ) ) {
		%parameterFiles4 = %{ get_raw_files( $config, $section, "source4" ) };
	}
	my %parameterFiles5 = ();
	if ( has_raw_files( $config, $section, "source5" ) ) {
		%parameterFiles5 = %{ get_raw_files( $config, $section, "source5" ) };
	}

	my $shfile = $pbs_dir . "/${task_name}_perl.sh";
	open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
	if ($sh_direct) {
		print $sh "export MYCMD=\"bash\" \n";
	}
	else {
		print $sh
"type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
	}

	for my $sample_name ( sort keys %raw_files ) {
		my $parameterFile1 = "";
		my @sample_files   = @{ $raw_files{$sample_name} };
		my $parameterFileLabel =
		  get_option( $config, $section, "source1Label", "" );
		my $parameterFileLabelEach =
		  get_option( $config, $section, "source1LabelEach", "" );
		my $parameterFileRegex =
		  get_option( $config, $section, "source1Regex", "" );
		foreach my $fileEach (@sample_files) {
			if ( $parameterFileRegex ne "" ) {
				if ( $fileEach =~ /$parameterFileRegex/ ) {
					$fileEach = $1;
				}
			}
			if ( $parameterFile1 eq "" ) {
				$parameterFile1 = $parameterFileLabelEach . $fileEach;
			}
			else {
				$parameterFile1 =
				  $parameterFile1 . "," . $parameterFileLabelEach . $fileEach;
			}
		}
		$parameterFile1 = $parameterFileLabel . $parameterFile1;

		my $parameterFile2 = "";
		if ( defined $parameterFiles2{$sample_name} ) {
			my $parameterFileLabel =
			  get_option( $config, $section, "source2Label", "" );
			my $parameterFileLabelEach =
			  get_option( $config, $section, "source2LabelEach", "" );
			my $parameterFileRegex =
			  get_option( $config, $section, "source2Regex", "" );
			my @files = @{ $parameterFiles2{$sample_name} };
			foreach my $fileEach (@files) {
				if ( $parameterFileRegex ne "" ) {
					if ( $fileEach =~ /$parameterFileRegex/ ) {
						$fileEach = $1;
					}
				}
				if ( $parameterFile2 eq "" ) {
					$parameterFile2 = $parameterFileLabelEach . $fileEach;
				}
				else {
					$parameterFile2 =
					    $parameterFile2 . ","
					  . $parameterFileLabelEach
					  . $fileEach;
				}
			}
			$parameterFile2 = $parameterFileLabel . $parameterFile2;
		}

		my $parameterFile3 = "";
		if ( defined $parameterFiles3{$sample_name} ) {
			my $parameterFileLabel =
			  get_option( $config, $section, "source3Label", "" );
			my $parameterFileLabelEach =
			  get_option( $config, $section, "source3LabelEach", "" );
			my $parameterFileRegex =
			  get_option( $config, $section, "source3Regex", "" );
			my @files = @{ $parameterFiles3{$sample_name} };
			foreach my $fileEach (@files) {
				if ( $parameterFileRegex ne "" ) {
					if ( $fileEach =~ /$parameterFileRegex/ ) {
						$fileEach = $1;
					}
				}
				if ( $parameterFile3 eq "" ) {
					$parameterFile3 = $parameterFileLabelEach . $fileEach;
				}
				else {
					$parameterFile3 =
					    $parameterFile3 . ","
					  . $parameterFileLabelEach
					  . $fileEach;
				}
			}
			$parameterFile3 = $parameterFileLabel . $parameterFile3;
		}

		my $pbs_file   = $self->get_pbs_filename( $pbs_dir, $sample_name );
		my $pbs_name   = basename($pbs_file);
		my $log        = $self->get_log_filename( $log_dir, $sample_name );
		my $final_file = $sample_name . $output_file . $output_ext;

		my $log_desc = $cluster->get_log_description($log);
		my $pbs      = $self->open_pbs(
			$pbs_file,  $pbs_desc,   $log_desc,
			$path_file, $result_dir, $final_file
		);
		print $pbs "
$runProgram $option ${output_file_label}$final_file $parameterFile1 $parameterFile2 $parameterFile3
";

		$self->close_pbs( $pbs, $pbs_file );

		print $sh "\$MYCMD ./$pbs_name \n";
	}

	close $sh;

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	print
"!!!shell file $shfile created, you can run this shell file to submit all Perl tasks.\n";
}

sub result {
	my ( $self, $config, $section, $pattern ) = @_;

	my (
		$task_name, $path_file,  $pbs_desc, $target_dir, $log_dir,
		$pbs_dir,   $result_dir, $option,   $sh_direct
	) = get_parameter( $config, $section, 0 );

	my %raw_files = %{ get_raw_files( $config, $section ) };
	my $output_ext = get_option( $config, $section, "output_ext", "" );
	my $output_other_ext =
	  get_option( $config, $section, "output_other_ext", "" );
	my @output_other_exts;
	if ( $output_other_ext ne "" ) {
		@output_other_exts = split( ",", $output_other_ext );
	}

	my $result = {};
	for my $sample_name ( keys %raw_files ) {

		#    my @sample_files = @{ $raw_files{$sample_name} };
		#    my @result_files = ();
		#    for my $sampleFile (@sample_files) {
		#      my $name = basename($sampleFile);
		#		push( @result_files, "${result_dir}/${name}${output_ext}" );
		#    }
		my @result_files = ();
		push( @result_files, "${result_dir}/${sample_name}${output_ext}" );
		if ( $output_other_ext ne "" ) {
			foreach my $output_other_ext_each (@output_other_exts) {
				push( @result_files,
					"${result_dir}/${sample_name}${output_other_ext_each}" );
			}
		}

		$result->{$sample_name} = filter_array( \@result_files, $pattern );
	}
	return $result;
}

1;
