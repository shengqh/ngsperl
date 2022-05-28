#!/usr/bin/perl
package SRA::SRAUtils;

use strict;
use warnings;
use File::Basename;
use String::Util qw(trim);
use CQS::ConfigUtils;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(
  GsmToSrr
  SrrToUrl
  getGsmSrrMap
  getSraFiles)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub GsmToSrr {
  my ( $gsm ) = @_;

  #my $cmd = "grep $gsm $sraTable | grep -e \"^SRR\" | cut -f1";
  my $cmd = "esearch -db sra -query $gsm |efetch -format docsum |xtract -pattern DocumentSummary -element Run\@acc";
  print $cmd . "\n";
  my $res = ` $cmd `;
  $res = trim($res);
  $res =~ s/\n/ /g;
  return $res;
}

sub SrrToUrl {
  my ( $srr ) = @_;

  #my $cmd = "grep $gsm $sraTable | grep -e \"^SRR\" | cut -f1";
  my $cmd = "curl \"https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$srr&result=read_run&fields=fastq_ftp&format=tsv&download=true\"";
  print $cmd . "\n";
  my $res = ` $cmd `;
  $res = trim($res);
  $res =~ s/\n/ /g;
  $res =~ s/.+\sftp.sra.ebi.ac.uk/ftp.sra.ebi.ac.uk/g;
  my @result = split(';', $res);
  return \@result;
}

sub getGsmSrrMap {
  my $listFile = shift;
  open my $list_handle, "<$listFile";
  my $first_line = <$list_handle>;
  close $list_handle;
  my @columns = split( '\t', $first_line );
  my $srrIndex = first { $columns[$_] eq 'Run_s' } 0 .. $#columns;
  my $gsmIndex = first { $columns[$_] eq 'Sample_Name_s' } 0 .. $#columns;
  my $taMap = readDictionaryByIndex( $listFile, $gsmIndex, $srrIndex, 1 );
  return $taMap;
}

sub getSraFiles {
  my ( $config, $section ) = @_;
  my $result;
  if ( defined $config->{$section}{"list_file"} ) {
    my $taMap = getGsmSrrMap( $config->{$section}{"list_file"} );
    for my $gsm ( keys %$taMap ) {
      $result->{$gsm} = [ $taMap->{$gsm} ];
    }
  }
  else {
    if ( defined $config->{$section}{"source"} ) {
      my $res = $config->{$section}{"source"};
      if ( is_array($res) ) {
        for my $gsm (@$res) {
          $result->{$gsm} = [$gsm];
        }
      }
      else {
        $result = $res;
      }
    }
    else {
      my $fileSection = $config->{$section}{"source_ref"}[0];
      my $files       = $config->{$fileSection};
      if ( is_array($files) ) {
        for my $gsm (@$files) {
          $result->{$gsm} = [$gsm];
        }
      }
      else {
        $result = get_raw_files( $config, $section );
      }
    }
  }
  return $result;
}

1;
