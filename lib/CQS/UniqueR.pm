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
use File::Spec;

our @ISA = qw(CQS::UniqueTask);

my $directory;

sub new {
	my ($class) = @_;
	my $self = $class->SUPER::new();
	$self->{_name}   = "UniqueR";
	$self->{_suffix} = "_uniqueR";
	bless $self, $class;
	return $self;
}

sub perform {
	my ( $self, $config, $section ) = @_;

	my (
		$task_name, $path_file, $pbsDesc, $target_dir, $logDir,
		$pbsDir,    $resultDir, $option,  $sh_direct,  $cluster
	) = get_parameter( $config, $section );

	my $rtemplate =
	  get_option_value( $config->{$section}{rtemplate}, "rtemplate", 1 );
	my $is_absolute = File::Spec->file_name_is_absolute($rtemplate);
	if ( !$is_absolute ) {
		$rtemplate = dirname(__FILE__) . "/$rtemplate";
	}
	if ( !( -e $rtemplate ) ) {
		die("rtemplate $rtemplate defined but not exists!");
	}
	my $rCode = get_option( $config, $section, "rCode","" );
	my $output_file = get_option( $config, $section, "output_file", 0 );
	
	my $rfile = $resultDir . "/${task_name}.r";
	open( RF, ">$rfile" ) or die "Cannot create $rfile";
	open RT, "<$rtemplate" or die $!;
	
	if (defined($rCode)) {
		print RF $rCode."\n";
	}
	while(<RT>) {
		print RF $_;
	}

	my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
	my $pbsName = basename($pbsFile);
	my $log     = $self->logfile( $logDir, $task_name );

	my $log_desc = $cluster->get_log_desc($log);

	open( OUT, ">$pbsFile" ) or die $!;
	print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

R --vanilla --slave -f $rfile $output_file $option

echo finished=`date`
";
	close(OUT);

	print "!!!shell file $pbsFile created.\n";
}

sub result {
	my ( $self, $config, $section, $pattern ) = @_;

	my (
		$task_name, $path_file, $pbsDesc, $target_dir, $logDir,
		$pbsDir,    $resultDir, $option,  $sh_direct
	) = get_parameter( $config, $section );

	my $output_file = get_option( $config, $section, "output_file", 0 );
	my $result      = {};
	my @resultFiles = ();

	push( @resultFiles, "${resultDir}/${output_file}" );

	$result->{$task_name} = filter_array( \@resultFiles, $pattern );
	return $result;
}

1;
