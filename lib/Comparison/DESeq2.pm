#!/usr/bin/perl
package Comparison::DESeq2;

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

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $comparisons = get_raw_files( $config, $section );
  my $totalPair = scalar( keys %{$comparisons} );
  if ( 0 == $totalPair ) {
    die "No pair defined!";
  }

  my $groups = get_raw_files( $config, $section, "groups" );

  my $countfile = parse_param_file( $config, $section, "countfile", 1 );

  my $rtemplate = dirname(__FILE__) . "/DESeq2.r";
  if ( !-e $rtemplate ) {
    die "File not found : " . $rtemplate;
  }

  my $showLabelInPCA = get_option( $config, $section, "showLabelInPCA", 1 );
  my $showDEGeneCluster = get_option( $config, $section, "showLabelInPCA", 0 );
  my $pvalue = get_option( $config, $section, "pvalue", 0.05 );
  my $foldChange = get_option( $config, $section, "foldChange", 2.0 );

  my %tpgroups = ();
  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    $tpgroups{$groupName} = "\"" . join( "\",\"", @samples ) . "\"";
  }

  my $rfile = $resultDir . "/${task_name}.r";
  open( RF, ">$rfile" ) or die "Cannot create $rfile";
  open RT, "<$rtemplate" or die $!;

  my $readfunc;
  if ( $countfile =~ /csv$/ ) {
    $readfunc = "read.csv";
  }
  else {
    $readfunc = "read.table";
  }
  print RF "
setwd(\"$resultDir\")  
  
data<-${readfunc}(\"$countfile\",row.names=1, header=T, check.names=F)

showLabelInPCA<-${showLabelInPCA}
showDEGeneCluster<-${showDEGeneCluster}
pvalue<-${pvalue}

comparisons=list(";
  my $first = 0;
  for my $comparisonName ( sort keys %{$comparisons} ) {
    $first++;

    my $gNames = $comparisons->{$comparisonName};
    my @groupNames;
    
    my $paired = 0;
    my $groupIds = {};
    if ( ref $gNames eq ref {} ) {
      @groupNames = @{$gNames->{groups}};
      $paired = defined $gNames->{paired};
      if($paired){
        $groupIds = $gNames->{paired};
      }
    }
    else {
      @groupNames = @{$gNames};
    }
    
    print(Dumper(@groupNames));
    
    if ( scalar(@groupNames) != 2 ) {
      die "Comparison of $comparisonName should contains and only contains two groups!";
    }
    
    my $g1 = $groupNames[0];
    my $g2 = $groupNames[1];
    my @s1 = @{ $groups->{$g1} };
    my @s2 = @{ $groups->{$g2} };

    my $filename = "${comparisonName}.design";
    if ( $first != 1 ) {
      print RF ",";
    }
    print RF "
  \"${comparisonName}\" = c(\"$filename\", \"$g1\", \"$g2\")";

    my $cdfile = $resultDir . "/$filename";
    open( CD, ">$cdfile" ) or die "Cannot create $cdfile";
    if ($paired) {
      print CD "Sample\tPaired\tCondition\n";
      for my $i ( 0 .. $#s1 ) {
        my $sname = $s1[$i];
        my $id    = $groupIds->[$i];
        print CD "${sname}\t${id}\t${g1}\n";
      }
      for my $i ( 0 .. $#s2 ) {
        my $sname = $s2[$i];
        my $id    = $groupIds->[$i];
        print CD "${sname}\t${id}\t${g2}\n";
      }
    }
    else {
      print CD "Sample\tCondition\n";
      for my $i ( 0 .. $#s1 ) {
        my $sname = $s1[$i];
        print CD "${sname}\t${g1}\n";
      }
      for my $i ( 0 .. $#s2 ) {
        my $sname = $s2[$i];
        print CD "${sname}\t${g2}\n";
      }
    }
    close(CD);
  }
  print RF "
) \n\n";

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

  my $log_desc = $cluster->get_log_desc($log);

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

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

  my $comparisons = get_raw_files( $config, $section );

  my $result = {};
  for my $comparisonName ( sort keys %{$comparisons} ) {
    my @resultFiles = ();
    push( @resultFiles, $resultDir . "/${comparisonName}.csv" );
    push( @resultFiles, $resultDir . "/${comparisonName}.png" );
    $result->{$comparisonName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
