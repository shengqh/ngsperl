#!/usr/bin/perl
package CQS::DESeq2;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;

our @ISA = qw(CQS::Task);

my $directory;

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "DESeq2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $pairs = get_raw_files( $config, $section );

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
  print RF "
setwd(\"$resultDir\")  
  
data<-read.table(\"$countfile\",row.names=1, header=T, check.names=F)

pairs=list(
";
  my $first = 0;
  my $totalGroup = scalar(keys %{$pairs});
  for my $pairName ( sort keys %{$pairs} ) {
    $first ++;
    my @groupNames = @{ $pairs->{$pairName} };

    if ( scalar(@groupNames) != 2 ) {
      die "Comparison in pair $pairName should contains and only contains two groups!";
    }

    my $g1 = $groupNames[0];
    my $g2 = $groupNames[1];
    my $s1 = $tpgroups{$g1};
    my $s2 = $tpgroups{$g2};
    print RF "  \"$pairName\" = list(\"$g1\" = c($s1), \"$g2\" = c($s2))";
    if ( $first != $totalGroup) {
      print RF ", \n";
    }
    else{
      print RF " \n";
    }
  }
  print RF ") \n\n";

  open RT, "<$rtemplate" or die $!;
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

  my $shfile = $pbsDir . "/${task_name}_de2.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "cd $resultDir
R --vanilla -f $rfile
";
  close(SH);

  print "!!!shell file $shfile created, you can run this shell file to do calculation.\n";
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
