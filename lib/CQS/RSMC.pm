#!/usr/bin/perl
package CQS::RSMC;

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
  $self->{_name} = "RSMC";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $rsmcfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1 );
  my $source_type = $config->{$section}{source_type} or die "source_type is not defined in $section";

  my $rnaediting_db = $config->{$section}{rnaediting_db};
  if ( defined $rnaediting_db ) {
    $option = $option . " --rnaediting --rnaediting_db $rnaediting_db ";
  }

  my $rawFiles = get_raw_files( $config, $section );

  my %group_sample_map = ();

  my $fafile           = "";
  my $mpileupParameter = "";
  my $isbam            = lc($source_type) eq "bam";
  if ($isbam) {
    $fafile = get_param_file( $config->{$section}{fasta_file}, "fasta_file (for mpileup)", 1 );
    $mpileupParameter = $config->{$section}{mpileup_option};
    if ( defined $mpileupParameter ) {
      if ( $mpileupParameter eq "" ) {
        undef($$mpileupParameter);
      }
    }

    my $groups = get_raw_files( $config, $section, "groups" );
    for my $groupName ( sort keys %{$groups} ) {
      my @samples = @{ $groups->{$groupName} };
      my @gfiles  = ();
      my $index   = 0;
      foreach my $sampleName (@samples) {
        my @bamFiles = @{ $rawFiles->{$sampleName} };
        push( @gfiles, $bamFiles[0] );
      }
      $group_sample_map{$groupName} = \@gfiles;
    }
  }
  else {
    %group_sample_map = %{$rawFiles};
  }

  my $shfile = $pbsDir . "/${task_name}_rs.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $groupName ( sort keys %group_sample_map ) {
    my @sampleFiles = @{ $group_sample_map{$groupName} };
    my $sampleCount = scalar(@sampleFiles);
    my $curDir      = create_directory_or_die( $resultDir . "/$groupName" );

    my $pbsName = "${groupName}_rs.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";
    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${groupName}_rs.log";
    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file 

echo RSMC=`date` 

cd $curDir
";

    if ($isbam) {
      if ( $sampleCount != 2 ) {
        die "SampleFile should be normal,tumor paired.";
      }

      my $normal = $sampleFiles[0];
      my $tumor  = $sampleFiles[1];

      my $cmd;
      if ( defined $mpileupParameter ) {
        $cmd = "samtools mpileup -f $fafile $mpileupParameter $normal $tumor | mono-sgen $rsmcfile all -t console $option -o $curDir";
      }
      else {
        $cmd = "mono-sgen $rsmcfile all -t bam -f $fafile $option -b $normal,$tumor -o $curDir";
      }

      print OUT "
if [ ! -s ${normal}.bai ]; then
  samtools index ${normal}
fi

if [ ! -s ${tumor}.bai ]; then
  samtools index ${tumor}
fi

$cmd

echo finished=`date` \n";
    }
    else {
      print OUT "mono-sgen $rsmcfile all -t mpileup -m $sampleFiles[0] $option";
    }

    close OUT;

    print "$pbsFile created \n";

  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all RSMC tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $groupName ( keys %{$groups} ) {
    my @resultFiles = ();
    my $curDir      = $resultDir . "/$groupName";
    push( @resultFiles, "$curDir/${groupName}.somatic.pass.vcf" );
    $result->{$groupName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
