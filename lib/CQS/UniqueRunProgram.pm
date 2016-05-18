#!/usr/bin/perl
package CQS::UniqueRunProgram;

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
	$self->{_suffix} = "_uniqueRunPro";
	bless $self, $class;
	return $self;
}

sub perform {
	my ( $self, $config, $section ) = @_;

	my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

	$self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
	my $task_suffix=get_option( $config, $section, "suffix", "" );
	$self->{_task_suffix} = $task_suffix;

    my $output_file_label     = get_option( $config, $section, "output_file_label",     "" );
	my $output_file     = get_option( $config, $section, "output_file",     "" );
	my $output_file_ext = get_option( $config, $section, "output_file_ext", "" );

	my $parametersample_files1 = "";
	if ( has_raw_files( $config, $section, "parameterSampleFile1" ) ) {
		my %temp = %{ get_raw_files( $config, $section, "parameterSampleFile1" ) };
		my @orderedSampleNames;
		my $parameterSampleFileOrder=$config->{$section}{parameterSampleFile1Order};
		if (defined $parameterSampleFileOrder) {
			@orderedSampleNames=@{$parameterSampleFileOrder};
		} else {
			@orderedSampleNames=sort keys %temp;
		}
		my $parameterSampleFileLabel = get_option( $config, $section, "parameterSampleFile1Label", "" );
		my $parameterSampleFileLabelEach = get_option( $config, $section, "parameterSampleFile1LabelEach", "" );
		my $parameterSampleFileSepEach = get_option( $config, $section, "parameterSampleFile1SepEach", " " );
		my $parameterSampleFileRegex = get_option( $config, $section, "parameterSampleFile1Regex", "" );
		foreach my $sample_name ( @orderedSampleNames ) {
			foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
				if ($parameterSampleFileRegex ne "") {
					if ($subSampleFile=~/$parameterSampleFileRegex/) {
						$subSampleFile=$1;
					}
				}
				if ( $parametersample_files1 eq "" ) {
					$parametersample_files1 = $parameterSampleFileLabelEach.$subSampleFile;
				}
				else {
					$parametersample_files1 =
					    $parametersample_files1 . $parameterSampleFileSepEach
					  . $parameterSampleFileLabelEach.$subSampleFile;
				}
			}
		}
		if ($parameterSampleFileLabel ne "") {
			$parametersample_files1=$parameterSampleFileLabel.$parametersample_files1
		}
	}
	
    my $parametersample_files2 = "";
    if ( has_raw_files( $config, $section, "parameterSampleFile2" ) ) {
        my %temp = %{ get_raw_files( $config, $section, "parameterSampleFile2" ) };
        my @orderedSampleNames;
        my $parameterSampleFileOrder=$config->{$section}{parameterSampleFile2Order};
        if (defined $parameterSampleFileOrder) {
            @orderedSampleNames=@{$parameterSampleFileOrder};
        } else {
            @orderedSampleNames=sort keys %temp;
        }
        my $parameterSampleFileLabel = get_option( $config, $section, "parameterSampleFile2Label", "" );
        my $parameterSampleFileLabelEach = get_option( $config, $section, "parameterSampleFile2LabelEach", "" );
        my $parameterSampleFileSepEach = get_option( $config, $section, "parameterSampleFile2SepEach", " " );
        my $parameterSampleFileRegex = get_option( $config, $section, "parameterSampleFile2Regex", "" );
        foreach my $sample_name ( @orderedSampleNames ) {
            foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
                if ($parameterSampleFileRegex ne "") {
                    print($subSampleFile."\n");
                    if ($subSampleFile=~/$parameterSampleFileRegex/) {
                        $subSampleFile=$1;
                    }
                }
				if ( $parametersample_files2 eq "" ) {
					$parametersample_files2 = $parameterSampleFileLabelEach.$subSampleFile;
				}
				else {
					$parametersample_files2 =
					    $parametersample_files2 . $parameterSampleFileSepEach
					  . $parameterSampleFileLabelEach.$subSampleFile;
				}            }
        }
        if ($parameterSampleFileLabel ne "") {
            $parametersample_files2=$parameterSampleFileLabel.$parametersample_files2
        }
    }
    
    my $parametersample_files3 = "";
    if ( has_raw_files( $config, $section, "parameterSampleFile3" ) ) {
        my %temp = %{ get_raw_files( $config, $section, "parameterSampleFile3" ) };
        my @orderedSampleNames;
        my $parameterSampleFileOrder=$config->{$section}{parameterSampleFile3Order};
        if (defined $parameterSampleFileOrder) {
            @orderedSampleNames=@{$parameterSampleFileOrder};
        } else {
            @orderedSampleNames=sort keys %temp;
        }
        my $parameterSampleFileLabel = get_option( $config, $section, "parameterSampleFile3Label", "" );
        my $parameterSampleFileLabelEach = get_option( $config, $section, "parameterSampleFile3LabelEach", "" );
        my $parameterSampleFileSepEach = get_option( $config, $section, "parameterSampleFile3SepEach", " " );
        my $parameterSampleFileRegex = get_option( $config, $section, "parameterSampleFile3Regex", "" );
        foreach my $sample_name ( @orderedSampleNames ) {
            foreach my $subSampleFile ( @{ $temp{$sample_name} } ) {
                if ($parameterSampleFileRegex ne "") {
                    print($subSampleFile."\n");
                    if ($subSampleFile=~/$parameterSampleFileRegex/) {
                        $subSampleFile=$1;
                    }
                }
				if ( $parametersample_files3 eq "" ) {
					$parametersample_files3 = $parameterSampleFileLabelEach.$subSampleFile;
				}
				else {
					$parametersample_files3 =
					    $parametersample_files3 . $parameterSampleFileSepEach
					  . $parameterSampleFileLabelEach.$subSampleFile;
				}            }
        }
        if ($parameterSampleFileLabel ne "") {
            $parametersample_files3=$parameterSampleFileLabel.$parametersample_files3
        }
    }
    
    my $parametersample_files=$parametersample_files1.' '.$parametersample_files2.' '.$parametersample_files3;

    my $parameterFile1Label = get_option( $config, $section, "parameterFile1Label", "" );
	my $parameterFile1 = parse_param_file( $config, $section, "parameterFile1", 0 );
    my $parameterFile2Label = get_option( $config, $section, "parameterFile2Label", "" );
	my $parameterFile2 = parse_param_file( $config, $section, "parameterFile2", 0 );
	my $parameterFile3 = parse_param_file( $config, $section, "parameterFile3", 0 );
	if ( defined($parameterFile1) ) {
		if ($parameterFile1Label ne "") {
			$parameterFile1=$parameterFile1Label.$parameterFile1
		}
	} else {
		$parameterFile1 = "";
	}
	if ( defined($parameterFile2) ) {
		if ($parameterFile2Label ne "") {
			$parameterFile2=$parameterFile2Label.$parameterFile2
		}
	} else {
		$parameterFile2 = "";
	}
	if ( !defined($parameterFile3) ) {
		$parameterFile3 = "";
	}
    my $parameterFiles=$parameterFile1.' '.$parameterFile2.' '.$parameterFile3;

	my $final_file    = "${task_name}${output_file}${output_file_ext}";

	my $runProgram = get_option( $config, $section, "runProgram");
	my $is_absolute = File::Spec->file_name_is_absolute($runProgram);
	if ( !$is_absolute ) {
			$runProgram = dirname(__FILE__) . "/$runProgram";
	}
	if ( !( -e $runProgram ) ) {
			die("runProgram $runProgram defined but not exists!");
	}

	my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
	my $pbs_name = basename($pbs_file);
	my $log      = $self->get_log_filename( $log_dir, $task_name );

	my $log_desc = $cluster->get_log_description($log);

	my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
	print $pbs "$runProgram $option $parametersample_files $parameterFiles ${output_file_label}$final_file";
	$self->close_pbs( $pbs, $pbs_file );
}

sub result {
	my ( $self, $config, $section, $pattern ) = @_;

	my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

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
