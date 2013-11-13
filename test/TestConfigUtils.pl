#!/usr/bin/perl
use strict;
use warnings;

use CQS::ConfigUtils;
use Test::Simple;
use Test::Deep;

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => "Task" },
  fastqfiles => {
    "P2177-01" => [ "/data/cqs/shengq1/2177/2177-WE-1_1_sequence.txt", "/data/cqs/shengq1/2177/2177-WE-1_2_sequence.txt" ],
    "P2177-02" => [ "/data/cqs/shengq1/2177/2177-WE-2_1_sequence.txt", "/data/cqs/shengq1/2177/2177-WE-2_2_sequence.txt" ],
  },
  fastqfiles2 => {
    "P2177-03" => [ "/data/cqs/shengq1/2177/2177-WE-3_1_sequence.txt", "/data/cqs/shengq1/2177/2177-WE-3_2_sequence.txt" ],
    "P2177-04" => [ "/data/cqs/shengq1/2177/2177-WE-4_1_sequence.txt", "/data/cqs/shengq1/2177/2177-WE-4_2_sequence.txt" ],
  },
  tophat2 => {
    class      => "Tophat2",
    target_dir => "tophat2",
    source_ref => "fastqfiles",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  tophat2_2 => {
    class      => "Tophat2",
    target_dir => "tophat2_2",
    source_ref => "fastqfiles2",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  cufflink => {
    class      => "Cufflinks",
    target_dir => "cufflink",
    source_ref => [ "tophat2", "tophat2_2" ],
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  cuffmerge => {
    class      => "Cuffmerge",
    target_dir => "cuffmerge",
    source_ref => "cufflink",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  cufflink_pattern => {
    class      => "Cufflinks",
    target_dir => "cufflink",
    source_ref => [ "tophat2", "tsv\$", "tophat2_2" ],
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

my $files = get_raw_files( $config, "tophat2" );
ok( eq_deeply( $files, $config->{fastqfiles} ) );

$files = get_raw_files( $config, "tophat2_2" );
ok( eq_deeply( $files, $config->{fastqfiles2} ) );

$files = get_raw_files( $config, "cufflink" );
ok(
  eq_deeply(
    $files,
    {
      "P2177-01" => ['tophat2/result/P2177-01/accepted_hits.bam'],
      "P2177-02" => ['tophat2/result/P2177-02/accepted_hits.bam'],
      "P2177-03" => ['tophat2_2/result/P2177-03/accepted_hits.bam'],
      "P2177-04" => ['tophat2_2/result/P2177-04/accepted_hits.bam']
    }
  )
);

$files = get_raw_files( $config, "cuffmerge" );
ok(
  eq_deeply(
    $files, $files,
    {
      "P2177-01" => ['cufflink/result/P2177-01/transcripts.gtf'],
      "P2177-02" => ['cufflink/result/P2177-02/transcripts.gtf'],
      "P2177-03" => ['cufflink/result/P2177-03/transcripts.gtf'],
      "P2177-04" => ['cufflink/result/P2177-04/transcripts.gtf'],
    }
  )
);

$files = get_raw_files( $config, "cufflink_pattern" );
ok(
  eq_deeply(
    $files,
    {
      "P2177-01" => [],
      "P2177-02" => [],
      "P2177-03" => ['tophat2_2/result/P2177-03/accepted_hits.bam'],
      "P2177-04" => ['tophat2_2/result/P2177-04/accepted_hits.bam']
    }
  )
);

