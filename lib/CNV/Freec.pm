#!/usr/bin/perl
package CNV::Freec;

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
  $self->{_name}   = "CNV::Freec";
  $self->{_suffix} = "_fc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my $snpFile = parse_param_file( $config, $section, "SNPfile", 0 );
  if ( defined $snpFile ) {
    $self->performPileup( $config, $section );
  }
  else {
    $self->performBAM( $config, $section );
  }
}

sub performBAM {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $chrLenFile             = $config->{$section}{chrLenFile}             or die "define ${section}::chrLenFile first";
  my $chrFiles               = $config->{$section}{chrFiles}               or die "define ${section}::chrFiles first, it is the path of directory contains chromosome fasta files";
  my $ploidy                 = $config->{$section}{ploidy}                 or die "define ${section}::ploidy first, such as 2";
  my $coefficientOfVariation = $config->{$section}{coefficientOfVariation} or die "define ${section}::coefficientOfVariation first, such as 0.05";
  my $inputFormat            = $config->{$section}{inputFormat}            or die "define ${section}::inputFormat first";

  my $mateOrientation = $config->{$section}{mateOrientation};
  if ( !defined $mateOrientation ) {
    die "define ${section}::mateOrientation first";    #0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)
  }

  my $bedfile = parse_param_file( $config, $section, "bedfile", 0 );

  my $raw_files = get_raw_files( $config, $section );
  my $groups;
  if ( defined $config->{$section}{"groups"} ) {
    $groups = get_raw_files( $config, $section, "groups" );
  }
  else {
    $groups = {};
    for my $sample_name ( sort keys %{$raw_files} ) {
      $groups->{$sample_name} = [$sample_name];
    }
  }

  my $shfile = $pbs_dir . "/${task_name}.sh";
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $group_name ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$group_name} };
    my $bam_file1;
    my $bam_file2;
    if ( scalar(@samples) > 1 ) {
      $bam_file2 = $raw_files->{ $samples[0] }[0];    #control
      $bam_file1 = $raw_files->{ $samples[1] }[0];    #sample
    }
    else {
      $bam_file1 = $raw_files->{ $samples[0] }[0];    #sample
    }

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    print $sh "\$MYCMD ./$pbs_name \n";
    my $cur_dir    = create_directory_or_die( $result_dir . "/$group_name" );
    my $configName = "${group_name}.conf";
    my $configFile = ${cur_dir} . "/$configName";

    my $log_desc   = $cluster->get_log_description($log);
    my $final_file = "${task_name}.call";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    print $pbs "   
freec -conf $configName 
";

    $self->close_pbs( $pbs, $pbs_file );

    open( my $con, ">$configFile" ) or die $!;

    print $con "[general] 
chrLenFile=$chrLenFile 
ploidy=$ploidy 
coefficientOfVariation=$coefficientOfVariation 
chrFiles=$chrFiles
BedGraphOutput=TRUE

[sample]
mateFile=$bam_file1 
inputFormat=$inputFormat 
mateOrientation=$mateOrientation 

";
    if ( defined $bam_file2 ) {
      print $con "[control] 
mateFile=$bam_file2 
inputFormat=$inputFormat 
mateOrientation=$mateOrientation 

";
    }

    if ( defined $bedfile ) {
      print $con "[target]
captureRegions = $bedfile

";
    }

    close $con;
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all freec tasks.\n";
}

sub performPileup {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $chrLenFile             = $config->{$section}{chrLenFile}             or die "define ${section}::chrLenFile first";
  my $chrFiles               = $config->{$section}{chrFiles}               or die "define ${section}::chrFiles first, it is the path of directory contains chromosome fasta files";
  my $ploidy                 = $config->{$section}{ploidy}                 or die "define ${section}::ploidy first, such as 2";
  my $coefficientOfVariation = $config->{$section}{coefficientOfVariation} or die "define ${section}::coefficientOfVariation first, such as 0.05";
  my $inputFormat            = $config->{$section}{inputFormat}            or die "define ${section}::inputFormat first";

  my $mateOrientation = $config->{$section}{mateOrientation}
    or die "define ${section}::mateOrientation first";    #0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)

  my $snpFile = parse_param_file( $config, $section, "SNPfile", 1 );
  my $bedfile = parse_param_file( $config, $section, "bedfile", 0 );

  my $fastaFile = parse_param_file( $config, $section, "fasta_file", 1 );

  my $raw_files = get_raw_files( $config, $section );
  my $groups = get_raw_files( $config, $section, "groups" );

  my $shfile = $pbs_dir . "/${task_name}.sh";
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $group_name ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$group_name} };
    my $bam_file1;
    my $bam_file2;
    my $pileup1;
    my $pileup2;
    if ( scalar(@samples) > 1 ) {
      $bam_file2 = $raw_files->{ $samples[0] }[0];    #control
      $pileup2   = $samples[0] . ".pileup.gz";
      $bam_file1 = $raw_files->{ $samples[1] }[0];    #sample
      $pileup1   = $samples[1] . ".pileup.gz";
    }
    else {
      $bam_file1 = $raw_files->{ $samples[0] }[0];    #sample
      $pileup1   = $samples[0] . ".pileup.gz";
    }

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    print $sh "\$MYCMD ./$pbs_name \n";
    my $cur_dir    = create_directory_or_die( $result_dir . "/$group_name" );
    my $configName = "${group_name}.conf";
    my $configFile = ${cur_dir} . "/$configName";

    my $log_desc   = $cluster->get_log_description($log);
    my $final_file = "${task_name}.call";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    print $pbs "
if [ ! -s $pileup1 ]; then
  samtools mpileup -A -f $fastaFile $bam_file1 | gzip > $pileup1
fi
";

    if ( defined $pileup2 ) {
      print $pbs "
if [ ! -s $pileup2 ]; then
  samtools mpileup -A -f $fastaFile $bam_file2 | gzip > $pileup2
fi
";
    }

    print $pbs "
freec -conf $configName 
";

    $self->close_pbs( $pbs, $pbs_file );

    open( my $con, ">$configFile" ) or die $!;

    print $con "[general] 
chrLenFile=$chrLenFile 
ploidy=$ploidy 
coefficientOfVariation=$coefficientOfVariation 
chrFiles=$chrFiles
BedGraphOutput=TRUE

[sample]
mateFile=$pileup1 
inputFormat=pileup 
mateOrientation=$mateOrientation 

";
    if ( defined $bam_file2 ) {
      print $con "[control] 
mateFile=$pileup2 
inputFormat=pileup 
mateOrientation=$mateOrientation 

";
    }

    if ( defined $bedfile ) {
      print $con "[target]
captureRegions = $bedfile

";
    }

    print $con "[BAF]
SNPfile = $snpFile
";
    close($con);
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all freec tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = { $task_name => [ $result_dir . "/${task_name}.call" ] };

  return $result;
}

1;
