#!/usr/bin/perl
package CQS::UniqueR;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::UniqueTask;
use CQS::StringUtils;
use File::Spec;

our @ISA = qw(CQS::UniqueTask);

my $directory;

sub new {
	my ($class) = @_;
	my $self = $class->SUPER::new();
	$self->{_name}   = __PACKAGE__;
	$self->{_suffix} = "_uniqueR";
	bless $self, $class;
	return $self;
}

sub perform {
	my ( $self, $config, $section ) = @_;

	my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

	$self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
	my $task_suffix=get_option( $config, $section, "suffix", "" );
	$self->{_task_suffix} = $task_suffix;

	my $rCode           = get_option( $config, $section, "rCode",           "" );
	my $output_file     = get_option( $config, $section, "output_file",     "" );
	my $output_file_ext = get_option( $config, $section, "output_file_ext", "" );

	my $parametersample_files1 = "";
	if ( has_raw_files( $config, $section, "parameterSampleFile1" ) ) {
	  
		my %temp = %{ get_raw_files( $config, $section, "parameterSampleFile1" ) };
    print Dumper(%temp);
		open( LIST, ">$result_dir/fileList1${task_suffix}.txt" ) or die "Cannot create fileList1.txt";
		foreach my $sample_name ( keys %temp ) {
			foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
				print LIST $subSampleFile . "\t$sample_name\n";
			}
		}
		$parametersample_files1 = "fileList1${task_suffix}.txt";
		close(LIST);
	}
	my $parametersample_files2 = "";
	if ( has_raw_files( $config, $section, "parameterSampleFile2" ) ) {
		my %temp = %{ get_raw_files( $config, $section, "parameterSampleFile2" ) };
		open( LIST, ">$result_dir/fileList2${task_suffix}.txt" ) or die "Cannot create fileList2.txt";
		foreach my $sample_name ( keys %temp ) {
			foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
				print LIST $subSampleFile . "\t$sample_name\n";
			}
		}
		$parametersample_files2 = "fileList2${task_suffix}.txt";
		close(LIST);
	}
	my $parametersample_files3 = "";
	if ( has_raw_files( $config, $section, "parameterSampleFile3" ) ) {
		my %temp = %{ get_raw_files( $config, $section, "parameterSampleFile3" ) };
		open( LIST, ">$result_dir/fileList3${task_suffix}.txt" ) or die "Cannot create fileList3.txt";
		foreach my $sample_name ( keys %temp ) {
			foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
				print LIST $subSampleFile . "\t$sample_name\n";
			}
		}
		$parametersample_files3 = "fileList3${task_suffix}.txt";
		close(LIST);
	}

	my $parameterFile1 = parse_param_file( $config, $section, "parameterFile1", 0 );
	my $parameterFile2 = parse_param_file( $config, $section, "parameterFile2", 0 );
	my $parameterFile3 = parse_param_file( $config, $section, "parameterFile3", 0 );
	if ( !defined($parameterFile1) ) {
		$parameterFile1 = "";
	}
	if ( !defined($parameterFile2) ) {
		$parameterFile2 = "";
	}
	if ( !defined($parameterFile3) ) {
		$parameterFile3 = "";
	}

	my $rfile = $result_dir . "/${task_name}${task_suffix}.r";
	open( my $rf, ">$rfile" ) or die "Cannot create $rfile";

	my $final_file    = "";
	my $output_file_r = "";
	if ( $output_file eq "parameterSampleFile1" or $output_file eq "parameterSampleFile2" or $output_file eq "parameterSampleFile3" ) {
		if ( has_raw_files( $config, $section, $output_file ) ) {
			my %temp = %{ get_raw_files( $config, $section, $output_file ) };
			foreach my $sample_name ( keys %temp ) {
				foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
					$final_file = "${subSampleFile}${output_file_ext}";
				}
			}
		}
		else {
			die "output_file defined as " . $output_file . ", but " . $output_file . " not in configure\n";
		}
	}
	else {
		$output_file_r = $task_name . $output_file;
		$final_file    = "${task_name}${output_file}${output_file_ext}";
	}

	my $rParameter = "outFile='$output_file_r'\n";
	if ( defined($parametersample_files1) ) {
		$rParameter = $rParameter . "parSampleFile1='$parametersample_files1'\n";
	}
	if ( defined($parametersample_files2) ) {
		$rParameter = $rParameter . "parSampleFile2='$parametersample_files2'\n";
	}
	if ( defined($parametersample_files3) ) {
		$rParameter = $rParameter . "parSampleFile3='$parametersample_files3'\n";
	}
	if ( defined($parameterFile1) ) {
		$rParameter = $rParameter . "parFile1='$parameterFile1'\n";
	}
	if ( defined($parameterFile2) ) {
		$rParameter = $rParameter . "parFile2='$parameterFile2'\n";
	}
	if ( defined($parameterFile3) ) {
		$rParameter = $rParameter . "parFile3='$parameterFile3'\n";
	}

	if ( $rParameter ne "" ) {
		print $rf $rParameter;
	}
	if ( defined($rCode) ) {
		print $rf $rCode . "\n";
	}

	my $rtemplates = get_option_value( $config->{$section}{rtemplate}, "rtemplate", 1 );
	my @rtemplates = split( ",|;", $rtemplates );
	foreach my $rtemplate (@rtemplates) {
		my $is_absolute = File::Spec->file_name_is_absolute($rtemplate);
		if ( !$is_absolute ) {
			$rtemplate = dirname(__FILE__) . "/$rtemplate";
		}
		if ( !( -e $rtemplate ) ) {
			die("rtemplate $rtemplate defined but not exists!");
		}
		open( my $rt, "<$rtemplate" ) or die $!;
		while (<$rt>) {
			print $rf $_;
		}
		close($rt);
	}
	close($rf);

	my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
	my $pbs_name = basename($pbs_file);
	my $log      = $self->get_log_filename( $log_dir, $task_name );

	my $log_desc = $cluster->get_log_description($log);

	my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
	if ( defined($option) and $option ne "" ) {
		print $pbs "R --vanilla --slave -f $rfile --args $option";

	 #    "R --vanilla --slave -f $rfile --args $task_name$output_file $option $parametersample_files1 $parametersample_files2 $parametersample_files3 $parameterFile1 $parameterFile2 $parameterFile3";
	}
	else {
		print $pbs "R --vanilla --slave -f $rfile";
	}
	$self->close_pbs( $pbs, $pbs_file );
}

sub result {
	my ( $self, $config, $section, $pattern ) = @_;

	my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

	my $output_file     = get_option( $config, $section, "output_file",     "" );
	my $output_file_ext = get_option( $config, $section, "output_file_ext", "" );
	my $result          = {};
	my @result_files    = ();

	if ( $output_file eq "parameterSampleFile1" or $output_file eq "parameterSampleFile2" or $output_file eq "parameterSampleFile3" ) {
		if ( has_raw_files( $config, $section, $output_file ) ) {
			my %temp = %{ get_raw_files( $config, $section, $output_file ) };
			foreach my $sample_name ( keys %temp ) {
				foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
					push( @result_files, "${subSampleFile}${output_file_ext}" );
				}
			}
		}
		else {
			die "output_file defined as " . $output_file . ", but " . $output_file . " not in configure\n";
		}
	}
	else {
		push( @result_files, "${result_dir}/${task_name}${output_file}${output_file_ext}" );
	}

	$result->{$task_name} = filter_array( \@result_files, $pattern );
	return $result;
}

1;
