#!/usr/bin/perl
package GATK::MuTect2;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::Fasta;
use CQS::GroupTask;

our @ISA = qw(CQS::GroupTask);

sub new {
	my ($class) = @_;
	my $self = $class->SUPER::new();
	$self->{_name}   = __PACKAGE__;
	$self->{_suffix} = "_mt";
	bless $self, $class;
	return $self;
}

sub perform {
	my ( $self, $config, $section ) = @_;

	my (
		$task_name, $path_file, $pbs_desc,   $target_dir,
		$log_dir,   $pbs_dir,   $result_dir, $option,
		$sh_direct, $cluster,   $thread
	) = get_parameter( $config, $section );

	my $gatk_jar =
	  get_param_file( $config->{$section}{gatk_jar}, "gatk_jar", 1, not $self->using_docker() );
	my $faFile =
	  get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
	my $dbsnpfile =
	  get_param_file( $config->{$section}{dbsnp_file}, "dbsnp_file", 1 );

	#mouse genome has no corresponding cosmic database
	my $cosmic_param = "";
	my $cosmicfile =
	  get_param_file( $config->{$section}{cosmic_file}, "cosmic_file", 0 );
	if ( defined $cosmicfile ) {
		$cosmic_param = "--cosmic $cosmicfile";
	}

	#normal panel vcf
	my $normalPanel_param = "";
	my $normalPanelfile =
	  get_param_file( $config->{$section}{normal_panel_file}, "normal_panel_file", 0 );
	if ( defined $normalPanelfile ) {
		$normalPanel_param = "--normal_panel $normalPanelfile";
	}
	
	#target region
	my $bedFile =
	  get_param_file( $config->{$section}{bed_file}, "bed_file", 0 );
	my $interval_padding =
	  get_option( $config, $section, "interval_padding", 0 );
	my $restrict_intervals = "";
	if ( defined $bedFile and $bedFile ne "" ) {
		if ( defined $interval_padding and $interval_padding != 0 ) {
			$restrict_intervals = "-L $bedFile -ip $interval_padding";
		}
		else {
			$restrict_intervals = "-L $bedFile";
		}
	}

	my $java = get_java( $config, $section );

	my $java_option = get_option( $config, $section, "java_option", "" );

	#Sample and Group
	my %raw_files = %{ get_raw_files( $config, $section ) };
	my %group_sample_map;
	if ( has_raw_files( $config, $section, "groups" ) ) {

		#%group_sample_map = %{ get_raw_files( $config, $section, "groups" ) };
		%group_sample_map = %{ get_group_sample_map( $config, $section ) };
	}
	else {
		%group_sample_map = %raw_files;
	}

	my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
	open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
	print $sh get_run_command($sh_direct) . "\n";
	print $sh "cd $pbs_dir\n";

	for my $group_name ( sort keys %group_sample_map ) {
		my @sample_files = @{ $group_sample_map{$group_name} };
		print join( "\n", @sample_files );
		my $sampleCount = scalar(@sample_files);
		my $cur_dir = create_directory_or_die( $result_dir . "/$group_name" );

		my $normal;
		my $tumor;
		my $sample_parm;
		if ( $sampleCount == 1 ) {
			$normal      = "";
			$tumor       = $sample_files[0];
			$sample_parm = "-I:tumor $tumor";
			if ($normalPanel_param eq "") { #Only one sample, no normal panel, then is is a normal only sample, need to add --artifact_detection_mode   
				$sample_parm=" --artifact_detection_mode ".$sample_parm;
			}
		}
		elsif ( $sampleCount != 2 ) {
			die "SampleFile should be tumor only or normal,tumor paired.";
		}
		else {
			$normal      = $sample_files[0][1];
			$tumor       = $sample_files[1][1];
			$sample_parm = "-I:tumor $tumor -I:normal $normal";
		}
		print "\n" . "Normal: " . $normal . "\nTumor: " . $tumor;

		my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
		my $pbs_name = basename($pbs_file);
		my $log      = $self->get_log_filename( $log_dir, $group_name );

		print $sh "\$MYCMD ./$pbs_name \n";

		#		my $out_file = "${group_name}.somatic.out";
		my $vcf = "${group_name}.somatic.vcf";
		my $vcfPASS="${group_name}.somatic.PASS.vcf";

		#		my $passvcf  = "${group_name}.somatic.pass.vcf";

		my $log_desc = $cluster->get_log_description($log);

		my $pbs =
		  $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file,
			$cur_dir );
		if ( $sampleCount == 2 ) {
			print $pbs "
if [ ! -s ${normal}.bai ]; then
  samtools index ${normal}
fi

";
		}
		print $pbs "
if [ ! -s ${tumor}.bai ]; then
  samtools index ${tumor}
fi

if [ ! -s $vcf ]; then
  $java $java_option -jar $gatk_jar $option -T MuTect2 -R $faFile $cosmic_param $normalPanel_param --dbsnp $dbsnpfile $sample_parm $restrict_intervals -o $vcf
  cat $vcf | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $vcfPASS
fi 
";

		$self->close_pbs( $pbs, $pbs_file );
		print "$pbs_file created \n";
	}
	close $sh;

	if ( is_linux() ) {
		chmod 0755, $shfile;
	}

	print
"!!!shell file $shfile created, you can run this shell file to submit all "
	  . $self->{_name}
	  . " tasks.\n";
}

sub result {
	my ( $self, $config, $section, $pattern ) = @_;

	my (
		$task_name, $path_file,  $pbs_desc, $target_dir, $log_dir,
		$pbs_dir,   $result_dir, $option,   $sh_direct
	) = get_parameter( $config, $section, 0 );

	my $groups = get_raw_files( $config, $section, "groups" );

	my $result = {};
	for my $group_name ( keys %{$groups} ) {
		my @result_files = ();
		my $cur_dir      = $result_dir . "/$group_name";
		push( @result_files, "$cur_dir/${group_name}.somatic.PASS.vcf" );
		push( @result_files, "$cur_dir/${group_name}.somatic.vcf" );
		$result->{$group_name} = filter_array( \@result_files, $pattern );
	}
	return $result;
}

1;
