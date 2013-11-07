#!/usr/bin/perl
package VarScan2::Mpileup2snp;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::AbstractSomaticMutation;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::AbstractSomaticMutation);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "VarScan2::Mpileup2snp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $varscan2_jar = get_param_file( $config->{$section}{VarScan2_jar},  "VarScan2_jar",  1 );
  my $faFile     = get_param_file( $config->{$section}{fasta_file},  "fasta_file",  1 );

  my $mpileup_options = get_option($config, $section, "mpileup_options", "");
  my $java_option = get_option($config, $section, "java_option", "");

  my %group_sample_map = %{$self->get_group_sample_map ($config, $section)};

  my $shfile = $pbsDir . "/${task_name}_vs2.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $groupName ( sort keys %group_sample_map ) {
    my @sampleFiles = @{ $group_sample_map{$groupName} };
  
    my $curDir      = create_directory_or_die( $resultDir . "/$groupName" );

    my $normal = $sampleFiles[0][0];
    my $snpvcf = "${groupName}.snp.vcf";
    my $pbsName = "${groupName}_vs2.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${groupName}_vs2.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file 

echo varscan2=`date` 

cd $curDir

if [ ! -s ${normal}.bai ]; then
  samtools index ${normal}
fi

if [ ! -s $snpvcf ]; then
  samtools mpileup $mpileup_options -f $faFile $normal | java $java_option -jar $varscan2_jar mpileup2snp $option --output-vcf $snpvcf
fi

echo finished=`date` \n";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all VarScan2 tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $groupName ( keys %{$groups} ) {
    my @resultFiles = ();
    my $curDir      = $resultDir . "/$groupName";
    my $snpvcf = "${groupName}.snp.vcf";
    push( @resultFiles, "$curDir/${snpvcf}.Somatic.hc" );
    $result->{$groupName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
