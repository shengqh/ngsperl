#!/usr/bin/perl
use strict;
use warnings;

use CQS::ConfigUtils;
use Test::Simple;

my $config = {
	fastqfiles => {
		"P2177-01" => [ "/data/cqs/shengq1/2177/2177-WE-1_1_sequence.txt", "/data/cqs/shengq1/2177/2177-WE-1_2_sequence.txt" ],
		"P2177-02" => [ "/data/cqs/shengq1/2177/2177-WE-2_1_sequence.txt", "/data/cqs/shengq1/2177/2177-WE-2_2_sequence.txt" ],
	},
	fastqc => {
		target_dir => "fastqc",
		source     => {
			"P2177-01" => [ "/data/cqs/shengq1/2177/2177-WE-1_1_sequence.txt", "/data/cqs/shengq1/2177/2177-WE-1_2_sequence.txt" ],
			"P2177-02" => [ "/data/cqs/shengq1/2177/2177-WE-2_1_sequence.txt", "/data/cqs/shengq1/2177/2177-WE-2_2_sequence.txt" ],
		},
	},
	tophat2 => {
		target_dir => "tophat2",
		source_ref => "fastqfiles",
	},
	cufflink => {
		target_dir => "cufflink",
		source_ref => "tophat2",
	},
	cuffmerge => {
		target_dir => "cuffmerge",
		source_ref => "cufflink",
	},
};

my ( $files, $issource ) = get_raw_files( $config, "fastqc" );
ok( $issource eq 1 );

( $files, $issource ) = get_raw_files( $config, "tophat2" );
ok( $issource eq 0 );
ok ($files eq $config->{fastqfiles});

( $files, $issource ) = get_raw_files( $config, "cufflink" );
ok( $issource eq 0 );
ok ($files eq $config->{fastqfiles});

( $files, $issource ) = get_raw_files( $config, "cuffmerge" );
ok( $issource eq 0 );
ok ($files eq $config->{fastqfiles});
