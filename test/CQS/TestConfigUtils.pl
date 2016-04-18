#!/usr/bin/perl
package CQS::TestConfigUtils;

use strict;
use warnings;

use CQS::ConfigUtils;
use Test::More;
use Test::Deep;

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config1 = {
  general    => { task_name => "Task" },
  macs2peaks => {
    "EC_H3K27AC_CON" => [
      "/scratch/cqs/shengq1/chipseq/20151208_gse53998/macs2callpeak/result/EC_H3K27AC_CON/EC_H3K27AC_CON_treat_pileup.bdg",
      "/scratch/cqs/shengq1/chipseq/20151208_gse53998/macs2callpeak/result/EC_H3K27AC_CON/EC_H3K27AC_CON_control_lambda.bdg"
    ],
    "d17_static_iPSctrl2_H3K27ac" => [
      "/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/macs2callpeak/result/d17_static_iPSctrl2_H3K27ac/d17_static_iPSctrl2_H3K27ac_treat_pileup.bdg",
      "/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/macs2callpeak/result/d17_static_iPSctrl2_H3K27ac/d17_static_iPSctrl2_H3K27ac_control_lambda.bdg"
    ],
  },
  macs2bdgdiff => {
    class      => "Chipseq::MACS2Bdgdiff",
    perform    => 1,
    target_dir => "macs2bdgdiff",
    option     => "",
    source_ref => "macs2peaks",
    groups_ref => "groups",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },  
};

my $files1 = get_raw_files( $config1, "macs2bdgdiff", "source", "_treat_pileup.bdg");
ok(
  eq_deeply(
    $files1,
    {
      "EC_H3K27AC_CON" => ['/scratch/cqs/shengq1/chipseq/20151208_gse53998/macs2callpeak/result/EC_H3K27AC_CON/EC_H3K27AC_CON_treat_pileup.bdg'],
      "d17_static_iPSctrl2_H3K27ac" => ['/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/macs2callpeak/result/d17_static_iPSctrl2_H3K27ac/d17_static_iPSctrl2_H3K27ac_treat_pileup.bdg'],
    }
  )
);

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
    class      => "Alignment::Tophat2",
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
    class      => "Alignment::Tophat2",
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
    class      => "Cufflinks::Cufflinks",
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
    class      => "Cufflinks::Cuffmerge",
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
    class      => "Cufflinks::Cufflinks",
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
      "P2177-03" => ['tophat2_2/result/P2177-03/accepted_hits.bam'],
      "P2177-04" => ['tophat2_2/result/P2177-04/accepted_hits.bam']
    }
  )
);

my $pairs = {
  "HiSeq_vs_MiSeq" => {
    groups => [ "MiSeq", "HiSeq" ],
    paired => 1
  }
};

my ( $ispaired, $gNames ) = get_pair_groups( $pairs, "HiSeq_vs_MiSeq" );
ok($ispaired);

ok( eq_deeply( $gNames, [ "MiSeq", "HiSeq" ] ) );

my $config2 = {
  general    => { task_name => "Task" },
  cufflink => {
    class      => "Cufflinks",
    target_dir => "cufflink",
    source_config_ref => [$config, "tophat2", $config, "tophat2_2" ],
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

$files = get_raw_files( $config2, "cufflink" );
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

done_testing();

