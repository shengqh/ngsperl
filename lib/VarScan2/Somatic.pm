#!/usr/bin/perl
package VarScan2::Somatic;

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
  $self->{_name}   = "VarScan2::Somatic";
  $self->{_suffix} = "_vs2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $varscan2_jar = get_param_file( $config->{$section}{VarScan2_jar}, "VarScan2_jar", 1 );
  my $faFile       = get_param_file( $config->{$section}{fasta_file},   "fasta_file",   1 );

  my $mpileup_options = get_option( $config, $section, "mpileup_options", "" );

  my %group_sample_map = %{ $self->get_group_sample_map( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  my $java_option = get_option( $config, $section, "java_option", "" );

  for my $groupName ( sort keys %group_sample_map ) {
    my @sampleFiles = @{ $group_sample_map{$groupName} };
    my $sampleCount = scalar(@sampleFiles);

    if ( $sampleCount != 2 ) {
      die "SampleFile should be normal,tumor paired.";
    }

    my $curDir = create_directory_or_die( $resultDir . "/$groupName" );

    my $normal = $sampleFiles[0][1];
    my $tumor  = $sampleFiles[1][1];

    my $normalfile = $sampleFiles[0][0];
    my $tumorfile  = $sampleFiles[1][0];

    my $snpvcf   = "${groupName}.snp.vcf";
    my $indelvcf = "${groupName}.indel.vcf";

    my $normal_pileup = "samtools mpileup -q 1 -f $faFile $normal";
    my $tumor_pileup  = "samtools mpileup -q 1 -f $faFile $tumor";

    my $pbsFile = $self->pbsfile( $pbsDir, $groupName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

    print SH "\$MYCMD ./$pbsName \n";
    print SH "cd $pbsDir\n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file 

echo varscan2=`date` 

cd $curDir

if [ ! -s $snpvcf ]; then
  if [ ! -s ${normal}.bai ]; then
    samtools index ${normal}
  fi

  if [ ! -s ${tumor}.bai ]; then
    samtools index ${tumor}
  fi

  java $java_option -jar $varscan2_jar somatic <($normal_pileup) <($tumor_pileup) $groupName $option --output-vcf 
fi

java $java_option -jar $varscan2_jar processSomatic $snpvcf
grep -v \"^#\" $snpvcf | cut -f1 | uniq -c | awk '{print \$2\"\t\"\$1}' > ${snpvcf}.chromosome

java $java_option -jar $varscan2_jar processSomatic $indelvcf
grep -v \"^#\" $indelvcf | cut -f1 | uniq -c | awk '{print \$2\"\t\"\$1}' > ${indelvcf}.chromosome

echo finished=`date`
";
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
    push( @resultFiles, "$curDir/${groupName}.snp.Somatic.hc.vcf" );
    push( @resultFiles, "$curDir/${groupName}.indel.Somatic.hc.vcf" );
    $result->{$groupName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
