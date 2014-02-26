#!/usr/bin/perl
package CQS::DESeq2;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

my $directory;

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "DESeq2";
  $self->{_suffix} = "_de2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $pairs = get_raw_files( $config, $section );
  my $totalPair = scalar( keys %{$pairs} );
  if ( 0 == $totalPair ) {
    die "No pair defined!";
  }

  my $groups = get_raw_files( $config, $section, "groups" );

  my $countfile = parse_param_file( $config, $section, "countfile", 1 );

  my $rtemplate = dirname(__FILE__) . "/DESeq2.r";
  if ( !-e $rtemplate ) {
    die "File not found : " . $rtemplate;
  }

  my %tpgroups = ();
  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    $tpgroups{$groupName} = "\"" . join( "\",\"", @samples ) . "\"";
  }

  my $rfile = $resultDir . "/${task_name}.r";
  open( RF, ">$rfile" ) or die "Cannot create $rfile";
  open RT, "<$rtemplate" or die $!;

  print RF "
setwd(\"$resultDir\")  
  
data<-read.table(\"$countfile\",row.names=1, header=T, check.names=F)

pairs=list(
";
  my $first = 0;
  for my $pairName ( sort keys %{$pairs} ) {
    $first++;
    my ($ispaired, $gNames) = get_pair_groups( $pairs, $pairName );
    my @groupNames = @{$gNames};
    if ( scalar(@groupNames) != 2 ) {
      die "Comparison in pair $pairName should contains and only contains two groups!";
    }

    my $g1 = $groupNames[0];
    my $g2 = $groupNames[1];
    my $s1 = $tpgroups{$g1};
    my $s2 = $tpgroups{$g2};
    if ($ispaired) {
      print RF "  \"$pairName\" = list(\"$g1\" = c($s1), \"$g2\" = c($s2), \"paired\" = TRUE)";
    }
    else {
      print RF "  \"$pairName\" = list(\"$g1\" = c($s1), \"$g2\" = c($s2))";
    }

    if ( $first != $totalPair ) {
      print RF ", \n";
    }
    else {
      print RF " \n";
    }
  }
  print RF ") \n\n";

  while (<RT>) {
    if ( $_ =~ '^#' ) {
      next;
    }
    last;
  }
  while (<RT>) {
    print RF $_;
  }
  close(RT);
  close(RF);

  my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName = basename($pbsFile);
  my $log     = $self->logfile( $logDir, $task_name );

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

R --vanilla -f $rfile
";
  close(OUT);

  print "!!!shell file $pbsFile created.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $pairs = get_raw_files( $config, $section );

  my $result = {};
  for my $pairName ( sort keys %{$pairs} ) {
    my @resultFiles = ();
    push( @resultFiles, $resultDir . "/${pairName}.csv" );
    push( @resultFiles, $resultDir . "/${pairName}.png" );
    $result->{$pairName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
