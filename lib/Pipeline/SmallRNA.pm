#!/usr/bin/perl
package Pipeline::SmallRNA;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performSmallRNA performSmallRNATask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub getConfig{
	my ($def) = @_;

	create_directory_or_die( $def->{target_dir} );

	my $config = {
		general => { "task_name" => $def->{task_name}, },
		files   => $def->{files},
		fastq_remove_N => {
			class      => "CQS::FastqTrimmer",
			perform    => 1,
			target_dir => $def->{target_dir} . "/fastq_remove_N",
			option     => "-n -z",
			extension  => "_trim.fastq.gz",
			source_ref => "files",
			cqstools   => $def->{cqstools},
			sh_direct  => 1,
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "2",
				"mem"      => "10gb"
			},
		},
		fastqc_pre_trim => {
			class      => "QC::FastQC",
			perform    => 1,
			target_dir => $def->{target_dir} . "/fastqc_pre_trim",
			option     => "",
			source_ref => "fastq_remove_N",
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "2",
				"mem"      => "10gb"
			},
		},
		cutadapt => {
			class      => "Cutadapt",
			perform    => 1,
			target_dir => $def->{target_dir} . "/cutadapt",
			option     => "-O 10 -m " . $def->{min_read_length},
			source_ref => "fastq_remove_N",
			adaptor    => "TGGAATTCTCGGGTGCCAAGG",
			extension  => "_clipped.fastq",
			sh_direct  => 1,
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "24",
				"mem"      => "20gb"
			},
		},
		fastqc_post_trim => {
			class      => "QC::FastQC",
			perform    => 1,
			target_dir => $def->{target_dir} . "/fastqc_post_trim",
			option     => "",
			source_ref => [ "cutadapt", ".fastq.gz" ],
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "2",
				"mem"      => "10gb"
			},
		},
		fastq_len => {
			class      => "FastqLen",
			perform    => 1,
			target_dir => $def->{target_dir} . "/fastq_len",
			option     => "",
			source_ref => "cutadapt",
			cqstools   => $def->{cqstools},
			sh_direct  => 1,
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "24",
				"mem"      => "20gb"
			},
		},
		identical => {
			class      => "FastqIdentical",
			perform    => 1,
			target_dir => $def->{target_dir} . "/identical",
			option     => "",
			source_ref => [ "cutadapt", ".fastq.gz" ],
			cqstools   => $def->{cqstools},
			extension  => "_clipped_identical.fastq.gz",
			sh_direct  => 1,
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "24",
				"mem"      => "20gb"
			},
		},
		identical_NTA => {
			class        => "CQS::FastqMirna",
			perform      => 1,
			target_dir   => $def->{target_dir} . "/identical_NTA",
			option       => "-l " . $def->{min_read_length},
			source_ref   => [ "identical", ".fastq.gz\$" ],
			seqcount_ref => [ "identical", ".dupcount\$" ],
			cqstools     => $def->{cqstools},
			extension    => "_clipped_identical_NTA.fastq.gz",
			sh_direct    => 1,
			pbs          => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "24",
				"mem"      => "20gb"
			},
		},

		#not identical, for IGV
		bowtie1_genome_1mm_notidentical => {
			class         => "Bowtie1",
			perform       => 1,
			target_dir    => $def->{target_dir} . "/bowtie1_genome_1mm_notidentical",
			option        => $def->{bowtie1_option_1mm},
			source_ref    => [ "cutadapt", ".fastq.gz\$" ],
			bowtie1_index => $def->{bowtie1_index},
			samonly       => 0,
			sh_direct     => 0,
			pbs           => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=" . $def->{max_thread},
				"walltime" => "72",
				"mem"      => "40gb"
			},
		},

		#1 mismatch search, NTA
		bowtie1_genome_1mm_NTA => {
			class         => "Bowtie1",
			perform       => 1,
			target_dir    => $def->{target_dir} . "/bowtie1_genome_1mm_NTA",
			option        => $def->{bowtie1_option_1mm},
			source_ref    => [ "identical_NTA", ".fastq.gz\$" ],
			bowtie1_index => $def->{bowtie1_index},
			samonly       => 0,
			sh_direct     => 1,
			pbs           => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=" . $def->{max_thread},
				"walltime" => "72",
				"mem"      => "40gb"
			},
		},
		bowtie1_genome_1mm_NTA_mirna_count => {
			class           => "MirnaCount",
			perform         => 1,
			target_dir      => $def->{target_dir} . "/bowtie1_genome_1mm_NTA_mirna_count",
			option          => $def->{mirnacount_option},
			source_ref      => "bowtie1_1mm_NTA",
			fastq_files_ref => "identical_NTA",
			seqcount_ref    => [ "identical_NTA", ".dupcount\$" ],
			cqs_tools       => $def->{cqstools},
			gff_file        => $def->{mirna_coordinate},
			fasta_file      => $def->{mirna_fasta},
			samtools        => $def->{samtools},
			sh_direct       => 1,
			pbs             => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "72",
				"mem"      => "40gb"
			},
		},
		bowtie1_genome_1mm_NTA_mirna_table => {
			class      => "CQS::CQSMirnaNTATable",
			perform    => 1,
			target_dir => $def->{target_dir} . "/bowtie1_genome_1mm_NTA_mirna_table",
			option     => "",
			source_ref => [ "bowtie1_genome_1mm_NTA_mirna_count", ".mapped.xml" ],
			cqs_tools  => $def->{cqstools},
			prefix     => "miRNA_1mm_NTA_",
			sh_direct  => 1,
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "10",
				"mem"      => "10gb"
			},
		},

		#1 mismatch search
		bowtie1_genome_1mm => {
			class         => "Bowtie1",
			perform       => 1,
			target_dir    => $def->{target_dir} . "/bowtie1_genome_1mm",
			option        => $def->{bowtie1_option_1mm},
			source_ref    => [ "identical", ".fastq.gz\$" ],
			bowtie1_index => $def->{bowtie1_index},
			samonly       => 0,
			sh_direct     => 1,
			pbs           => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=" . $def->{max_thread},
				"walltime" => "72",
				"mem"      => "40gb"
			},
		},
		bowtie1_genome_1mm_miRNA_overlap => {
			class           => "CQSMappedCount",
			perform         => 1,
			target_dir      => $def->{target_dir} . "/bowtie1_genome_1mm_miRNA_overlap",
			option          => $def->{mirna_overlap_count_option},
			source_ref      => "bowtie1_genome_1mm",
			fastq_files_ref => "identical",
			seqcount_ref    => [ "identical", ".dupcount\$" ],
			cqs_tools       => $def->{cqstools},
			gff_file        => $def->{mirna_coordinate},
			fasta_file      => $def->{mirna_fasta},
			samtools        => $def->{samtools},
			sh_direct       => 1,
			pbs             => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "72",
				"mem"      => "20gb"
			},
		},
		bowtie1_genome_1mm_miRNA_overlap_position => {
			class      => "CQSMappedPosition",
			perform    => 1,
			target_dir => $def->{target_dir} . "/bowtie1_genome_1mm_miRNA_overlap_position",
			option     => "-o " . $def->{task_name} . "_miRNA.position",
			source_ref => "bowtie1_genome_1mm_miRNA_overlap",
			cqs_tools  => $def->{cqstools},
			sh_direct  => 1,
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "10",
				"mem"      => "10gb"
			},
		},
		bowtie1_genome_1mm_tRNA_count => {
			class           => "CQSMappedCount",
			perform         => 1,
			target_dir      => $def->{target_dir} . "/bowtie1_genome_1mm_tRNA_count",
			option          => $def->{trnacount_option},
			source_ref      => "bowtie1_genome_1mm",
			fastq_files_ref => "identical",
			seqcount_ref    => [ "identical", ".dupcount\$" ],
			cqs_tools       => $def->{cqstools},
			gff_file        => $def->{trna_coordinate},
			fasta_file      => $def->{trna_fasta},
			samtools        => $def->{samtools},
			sh_direct       => 1,
			pbs             => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "72",
				"mem"      => "20gb"
			},
		},
		bowtie1_genome_1mm_tRNA_table => {
			class      => "CQSMappedTable",
			perform    => 1,
			target_dir => $def->{target_dir} . "/bowtie1_genome_1mm_tRNA_table",
			option     => "",
			source_ref => [ "bowtie1_genome_1mm_tRNA_count", ".xml" ],
			cqs_tools  => $def->{cqstools},
			prefix     => "tRNA_1mm_",
			sh_direct  => 1,
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "10",
				"mem"      => "10gb"
			},
		},
		bowtie1_genome_1mm_tRNA_position => {
			class      => "CQSMappedPosition",
			perform    => 1,
			target_dir => $def->{target_dir} . "/bowtie1_genome_1mm_tRNA_position",
			option     => "-o " . $def->{task_name} . "_tRNA.position",
			source_ref => "bowtie1_genome_1mm_tRNA_count",
			cqs_tools  => $def->{cqstools},
			sh_direct  => 1,
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "10",
				"mem"      => "10gb"
			},
		},
		bowtie1_genome_1mm_smallRNA_count => {
			class           => "CQSMappedCount",
			perform         => 1,
			target_dir      => $def->{target_dir} . "/bowtie1_genome_1mm_smallRNA_count",
			option          => $def->{trnacount_option},
			source_ref      => "bowtie1_genome_1mm",
			fastq_files_ref => "identical",
			seqcount_ref    => [ "identical", ".dupcount\$" ],
			cqs_tools       => $def->{cqstools},
			gff_file        => $def->{smallrna_coordinate},
			samtools        => $def->{samtools},
			sh_direct       => 1,
			pbs             => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "72",
				"mem"      => "20gb"
			},
		},
		bowtie1_genome_1mm_smallRNA_table => {
			class      => "CQSMappedTable",
			perform    => 1,
			target_dir => $def->{target_dir} . "/bowtie1_genome_1mm_smallRNA_table",
			option     => "",
			source_ref => [ "bowtie1_genome_1mm_smallRNA_count", ".xml" ],
			cqs_tools  => $def->{cqstools},
			prefix     => "smallRNA_1mm_",
			sh_direct  => 1,
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "10",
				"mem"      => "10gb"
			},
		},
		bowtie1_genome_1mm_NTA_smallRNA_category => {
			class           => "CQSSmallRNACategory",
			perform         => 1,
			target_dir      => $def->{target_dir} . "/bowtie1_genome_1mm_NTA_smallRNA_category",
			option          => "",
			source_ref      => [ "bowtie1_genome_1mm_smallRNA_count", ".mapped.xml\$" ],
			mirna_count_ref => [ "bowtie1_genome_1mm_NTA_mirna_count", ".mapped.xml\$" ],
			cqs_tools       => $def->{cqstools},
			sh_direct       => 1,
			pbs             => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "72",
				"mem"      => "40gb"
			},
		},

		#2 perfect match search to mirbase only
		bowtie1_genome_pmnames => {
			class      => "Samtools::PerfectMappedReadNames",
			perform    => 1,
			target_dir => $def->{target_dir} . "/bowtie1_genome_pmnames",
			option     => "",
			source_ref => "bowtie1_genome_cutadapt_topN_1mm",
			sh_direct  => 1,
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=" . $def->{max_thread},
				"walltime" => "72",
				"mem"      => "40gb"
			},
		},
		bowtie1_miRbase_pm => {
			class         => "Alignment::Bowtie1",
			perform       => 1,
			target_dir    => $def->{target_dir} . "/bowtie1_miRbase_pm",
			option        => $def->{bowtie1_option_pm},
			source_ref    => [ "identical", ".fastq.gz\$" ],
			bowtie1_index => $def->{bowtie1_miRBase_index},
			samonly       => 0,
			sh_direct     => 1,
			pbs           => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=" . $def->{max_thread},
				"walltime" => "72",
				"mem"      => "40gb"
			},
		},
		bowtie1_miRbase_pm_count => {
			class                   => "CQS::CQSChromosomeCount",
			perform                 => 1,
			target_dir              => $def->{target_dir} . "/bowtie1_miRbase_pm_count",
			option                  => "",
			source_ref              => "bowtie1_miRbase_pm",
			seqcount_ref            => [ "identical", ".dupcount\$" ],
			perfect_mapped_name_ref => "bowtie1_genome_pmnames",
			cqs_tools               => $def->{cqstools},
			samtools                => $def->{samtools},
			sh_direct               => 1,
			pbs                     => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "72",
				"mem"      => "40gb"
			},
		},
		bowtie1_miRbase_pm_table => {
			class      => "CQS::CQSChromosomeTable",
			perform    => 1,
			target_dir => $def->{target_dir} . "/bowtie1_miRbase_pm_table",
			option     => "",
			source_ref => [ "bowtie1_miRbase_pm_count", ".xml" ],
			cqs_tools  => $def->{cqstools},
			prefix     => "miRBase_pm_",
			sh_direct  => 1,
			pbs        => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=1",
				"walltime" => "10",
				"mem"      => "10gb"
			},
		},
		sequencetask => {
			class      => "CQS::SequenceTask",
			perform    => 1,
			target_dir => $def->{target_dir} . "/sequencetask",
			option     => "",
			source     => {
				individual => [

					#data preparation
					"fastq_remove_N",   "fastqc_pre_trim", "cutadapt", "fastqc_post_trim", "fastq_len",
					"identical", "identical_NTA",

					#NTA data analysis
					"bowtie1_genome_1mm_NTA", "bowtie1_genome_1mm_NTA_miRNA_count",

					#non-NTA count
					"bowtie1_genome_1mm", "bowtie1_genome_1mm_miRNA_overlap",
					"bowtie1_genome_1mm_tRNA_count",                   "bowtie1_genome_1mm_smallRNA_count",

					#miRBase analysis
					"bowtie1_genome_pmnames",
					"bowtie1_miRbase_pm", "bowtie1_miRbase_pm_count",

					#for IGV
					"bowtie1_genome_1mm_notidentical",
				],
				summary => [

					#NTA table
					"bowtie1_genome_1mm_NTA_miRNA_table",

					#non-NTA table and graph
					"bowtie1_genome_1mm_NTA_miRNA_position", 
					"bowtie1_genome_1mm_NTA_tRNA_table",
					"bowtie1_genome_1mm_NTA_smallRNA_table",         
					"bowtie1_genome_1mm_NTA_smallRNA_category",
					"bowtie1_genome_1mm_tRNA_position", 

					#miRBase
					"bowtie1_miRbase_pm_table"
				],
			},
			sh_direct => 0,
			pbs       => {
				"email"    => $def->{email},
				"nodes"    => "1:ppn=" . $def->{max_thread},
				"walltime" => "72",
				"mem"      => "40gb"
			},
		},
	};

	return ($config);	
}

sub performSmallRNA {
	my ($def) = @_;
	
	my $config = getConfig($def);
	
	performConfig($config);

	1;
};

sub performSmallRNATask {
	my ($def, $task) = @_;
	
	my $config = getConfig($def);
	
	performTask($config, $task);

	1;
};
