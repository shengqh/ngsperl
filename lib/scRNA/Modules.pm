#!/usr/bin/perl
package scRNA::Modules;

use strict;
use warnings;
require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(addEnclone addClonotypeMerge)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub addEnclone {
  my ( $config, $def, $tasks, $taskName, $parentDir, $sourceRef ) = @_;

  $config->{$taskName} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$parentDir/$taskName",
    init_command          => '
dn=`dirname __FILE__`
',
    option                => "TCR=\${dn} POUT=__NAME__.csv > __NAME__.log",
    interpretor           => "",
    check_program         => 0,
    program               => "enclone",
    source_ref            => $sourceRef,
    source_arg            => "TCR=",
    source_join_delimiter => " ",
    output_to_same_folder => 1,
    output_to_folder      => 1,
    output_arg            => ">",
    output_file_ext       => ".csv",
    output_other_ext      => ".csv",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "4",
      "mem"       => "10gb"
    },
  };

  push(@$tasks, $taskName);
}

sub addClonotypeMerge {
  my ( $config, $def, $tasks, $target_dir, $taskname, $source_ref ) = @_;

  $config->{$taskname} = {
    class                    => "CQS::ProgramWrapper",
    perform                  => 1,
    target_dir               => "${target_dir}/$taskname",
    option                   => "",
    interpretor              => "python",
    program                  => "../scRNA/clonotype_merge.py",
    source_arg               => "-i",
    source_ref               => $source_ref,
    output_arg               => "-o",
    output_file              => "all_contig_annotations.json",
    output_no_name           => 1,
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push @$tasks, $taskname;
}

1;
