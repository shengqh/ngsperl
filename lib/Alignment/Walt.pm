#!/usr/bin/perl
package Alignment::Walt;

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
  $self->{_suffix} = "_walt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $selfname = $self->{_name};

  my $walt_index = $config->{$section}{walt_index};
  if ( !defined $walt_index ) {
    die "define ${section}::walt_index first";
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sample_files_str = ( scalar(@sample_files) == 2 ) ? "-1 " . $sample_files[0] . " -2 " . $sample_files[1]: "-1 " . $sample_files[0];

    my $result_file_unsorted          = $sample_name . ".unsorted.mr";
	my $result_file          = $sample_name . ".mr";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $rmlist = "";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $result_file );

			print $pbs "
if [ ! -s $result_file ]; then
  echo zcat=`date`
";
	if ($sample_files[0]=~/\.gz$/) { #.gz fle, need to zcat
		my @sample_files_unzip=();
		foreach my $sampleFile (@sample_files) {
			my $sampleFileUnzip = basename(change_extension( $sampleFile, "" ));
			push @sample_files_unzip,$sampleFileUnzip;
			print $pbs "
  if [ ! -s $sampleFileUnzip ]; then
      zcat $sampleFile > $sampleFileUnzip
  fi
";
			$rmlist=$rmlist. " $sampleFileUnzip"
		}
					print $pbs "
fi
";
		$sample_files_str = ( scalar(@sample_files_unzip) == 2 ) ? "-1 " . $sample_files_unzip[0] . " -2 " . $sample_files_unzip[1]: "-1 " . $sample_files_unzip[0];
	}
	
    print $pbs "
if [ ! -s $result_file_unsorted && ! -s $result_file ]; then
  echo walt=`date`
  walt -i $walt_index $sample_files_str -o $result_file_unsorted
fi
";
    print $pbs "
if [ ! -s $result_file ]; then
  echo walt=`date`
  sort -k1,1 -k2,2g -k3,3g -k6,6 $result_file_unsorted -o $result_file;
fi
";
    $rmlist=$rmlist ." $result_file_unsorted";

    if ($rmlist ne "") {
    	    print $pbs "
if [ -s $result_file_unsorted && -s $result_file ]; then
  rm $rmlist 
fi
";
    }
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all Walt tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $bam_file     = "${result_dir}/${sample_name}.mr";
    my @result_files = ();
    push( @result_files, $bam_file );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
