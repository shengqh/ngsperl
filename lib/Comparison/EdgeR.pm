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
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_er";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $comparisons = get_raw_files( $config, $section );
  my $totalPair = scalar( keys %{$comparisons} );
  if ( 0 == $totalPair ) {
    die "No pair defined!";
  }

  my $groups = get_raw_files( $config, $section, "groups" );

  my $countfile = parse_param_file( $config, $section, "countfile", 1 );

  my $rtemplate = dirname(__FILE__) . "/EdgeR.r";
  if ( !-e $rtemplate ) {
    die "File not found : " . $rtemplate;
  }

  my $showLabelInPCA = get_option( $config, $section, "show_label_PCA", 1 );
  my $showDEGeneCluster = get_option( $config, $section, "show_DE_gene_cluster", 0 );
  my $pvalue = get_option( $config, $section, "pvalue", 0.05 );
  my $foldChange = get_option( $config, $section, "fold_change", 2.0 );
  my $minMedianInGroup = get_option( $config, $section, "min_median_read", 0 );

  my %tpgroups = ();
  for my $group_name ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$group_name} };
    $tpgroups{$group_name} = "\"" . join( "\",\"", @samples ) . "\"";
  }

  my $rfile = $result_dir . "/${task_name}.r";
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
setwd(\"$result_dir\")  
  
data<-${readfunc}(\"$countfile\",row.names=1, header=T, check.names=F)

showLabelInPCA<-$showLabelInPCA
showDEGeneCluster<-$showDEGeneCluster
pvalue<-$pvalue
foldChange<-$foldChange
minMedianInGroup<-$minMedianInGroup

comparisons=list(";
  my $first = 0;
  for my $comparisonName ( sort keys %{$comparisons} ) {
    $first++;

    my $gNames = $comparisons->{$comparisonName};
    my @group_names;
    
    my $paired = 0;
    my $groupIds = {};
    if ( ref $gNames eq ref {} ) {
      @group_names = @{$gNames->{groups}};
      $paired = defined $gNames->{paired};
      if($paired){
        $groupIds = $gNames->{paired};
      }
    }
    else {
      @group_names = @{$gNames};
    }
    
    print(Dumper(@group_names));
    
    if ( scalar(@group_names) != 2 ) {
      die "Comparison of $comparisonName should contains and only contains two groups!";
    }
    
    my $g1 = $group_names[0];
    my $g2 = $group_names[1];
    my @s1 = @{ $groups->{$g1} };
    my @s2 = @{ $groups->{$g2} };

    my $filename = "${comparisonName}.design";
    if ( $first != 1 ) {
      print RF ",";
    }
    print RF "
  \"${comparisonName}\" = c(\"$filename\", \"$g1\", \"$g2\")";

    my $cdfile = $result_dir . "/$filename";
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

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log     = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  open( my $out, ">$pbs_file" ) or die $!;
  print $out "$pbs_desc
$log_desc

$path_file

cd $result_dir

R --vanilla -f $rfile
";
  close $out;

  print "!!!shell file $pbs_file created.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $comparisons = get_raw_files( $config, $section );

  my $result = {};
  for my $comparisonName ( sort keys %{$comparisons} ) {
    my @result_files = ();
    push( @result_files, $result_dir . "/${comparisonName}.csv" );
    push( @result_files, $result_dir . "/${comparisonName}.png" );
    $result->{$comparisonName} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
