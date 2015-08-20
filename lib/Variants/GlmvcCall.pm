#!/usr/bin/perl
package Variants::GlmvcCall;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "GlmvcCall";
  $self->{_suffix} = "_gc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $glmvcfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1 );
  my $source_type = $config->{$section}{source_type} or die "source_type is not defined in $section";

  my $rnaediting_db = get_directory( $config, $section, "rnaediting_db", 0 );
  if ( defined $rnaediting_db ) {
    $option = $option . " --rnaediting_db $rnaediting_db ";
  }

  my $annovar_buildver = $config->{$section}{annovar_buildver};
  if ( defined $annovar_buildver ) {
    $option = $option . " --annovar_buildver $annovar_buildver ";
  }

  my $distance_exon_gtf = get_param_file( $config->{$section}{distance_exon_gtf}, "distance_exon_gtf", 0 );
  if ( defined $annovar_buildver ) {
    $option = $option . " --distance_exon_gtf $distance_exon_gtf ";
  }

  my $anno = defined $rnaediting_db || defined $annovar_buildver || defined $annovar_buildver;

  my $rawFiles = get_raw_files( $config, $section );

  my %group_sample_map = ();

  my $fafile           = "";
  my $mpileupParameter = "";
  my $isbam            = lc($source_type) eq "bam";
  if ($isbam) {
    $fafile = get_param_file( $config->{$section}{fasta_file}, "fasta_file (for mpileup)", 1 );
    $mpileupParameter = $config->{$section}{mpileup_option};
    
    print "$mpileupParameter \n";
    if ( defined $mpileupParameter ) {
      if ( ! ($mpileupParameter eq "") ) {
        $mpileupParameter = "--mpileup \"" . $mpileupParameter . "\"";
      }
    }
    print "$mpileupParameter \n";

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

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";
  print SH "cd $pbsDir \n";

  for my $groupName ( sort keys %group_sample_map ) {
    my @sampleFiles = @{ $group_sample_map{$groupName} };
    my $sampleCount = scalar(@sampleFiles);
    my $curDir      = create_directory_or_die( $resultDir . "/$groupName" );

    my $pbsFile = $self->pbsfile( $pbsDir, $groupName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file 

echo Glmvc=`date` 

cd $curDir
";

    my $finalvcf = "${groupName}.vcf";

    if ($isbam) {
      if ( $sampleCount != 2 ) {
        die "SampleFile should be normal,tumor paired.";
      }

      my $normal = $sampleFiles[0];
      my $tumor  = $sampleFiles[1];
      my $final  = $anno ? "${groupName}.annotation.tsv" : "${groupName}.tsv";

      my $cmd = "mono-sgen $glmvcfile call -c $thread -t bam -f $fafile $option $mpileupParameter --normal $normal --tumor $tumor -o ${curDir}/${groupName}";

      print OUT "
if [ -s $final ]; then
  echo job has already been done. if you want to do again, delete ${curDir}/${final} and submit job again.
  exit 0;
fi      
      
if [ ! -s ${normal}.bai ]; then
  samtools index ${normal}
fi

if [ ! -s ${tumor}.bai ]; then
  samtools index ${tumor}
fi

$cmd

";
    }
    else {
      print OUT "mono-sgen $glmvcfile call -t mpileup -m $sampleFiles[0] $option -o ${curDir}/${groupName} \n";
    }

    print OUT "grep -v \"^#\" $finalvcf | cut -f1 | uniq -c | awk '{print \$2\"\\t\"\$1}' > ${finalvcf}.chromosome 
    
echo finished=`date`
";

    close OUT;

    print "$pbsFile created \n";

  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $rnaediting_db     = get_directory( $config, $section, "rnaediting_db", 0 );
  my $annovar_buildver  = $config->{$section}{annovar_buildver};
  my $distance_exon_gtf = get_param_file( $config->{$section}{distance_exon_gtf}, "distance_exon_gtf", 0 );
  my $anno              = defined $rnaediting_db || defined $annovar_buildver || defined $annovar_buildver;

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $groupName ( keys %{$groups} ) {
    my @resultFiles = ();
    my $curDir      = $resultDir . "/$groupName";
    push( @resultFiles, "$curDir/${groupName}.vcf" );
    push( @resultFiles, "$curDir/${groupName}.tsv" );
    if ($anno) {
      push( @resultFiles, "$curDir/${groupName}.annotation.tsv" );
    }
    $result->{$groupName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

sub pbsfiles {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $pairs = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $pairName ( sort keys %{$pairs} ) {
    my $pbsName = $self->pbsname($pairName);
    my $pbsFile = $pbsDir . "/$pbsName";
    $result->{$pairName} = $pbsFile;
  }
  return $result;
}

1;
