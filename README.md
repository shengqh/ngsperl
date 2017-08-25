NGSPERL : A semi-automated framework for large scale next generation sequencing data analysis
==========
* [Introduction](#introduction)
* [Citation](#citation)
* [Download and install](#download)
* [Framework](#framework)
* [Quick Start](#example)
* [Pipelines](#pipelines)
* [Modules](#module)

<a name="introduction"/>

# Introduction 

High-throughput sequencing technologies have been widely used in the research field, especially in cancer biology. With the huge amounts of sequencing data being generated, data analysis has become the bottle-neck of the research procedure. A lot of tools have been developed for different data types and different data analysis purposes while new tools are still being published every month. A software framework which not only supports large scale data analysis using the existing pipeline on cluster but can also easily replace/extend old modules in the pipeline will help solve the problem. We have designed and implemented NGSPERL, a semi-automated module-based framework, for high-throughput sequencing data analysis. Three major analysis pipelines with multiple tasks have been developed for RNA sequencing, exome sequencing, and small RNA sequencing data. The pipelines cover the tasks from raw data pre-processing, quality control, mapping, and comparison to report. Each task in the pipelines was developed as module. The module uses the output from the previous task in the pipeline as the input parameter to generate the corresponding portable batch system (PBS) scripts with other user-defined parameters. The module with such a trace-back design can be easily plugged into or unplugged from existing pipelines. The PBS scripts generated at each task can be submitted to cluster or run directly based on user choice. Multiple tasks can also be combined together as a single task to simplify the data analysis. Such a flexible framework will significantly accelerate the speed of large scale sequencing data analysis.

<a name="citation"/>

# Citation

Sheng Q, Zhao S, Guo M, Shyr Y: NGSPERL: a semi-automated framework for large scale next generation sequencing data analysis. International Journal of Computational Biology and Drug Design 2015, 8(3):203-211.

<a name="download"/>

# Download and install 

You can download NGSPERL package from [github](https://github.com/shengqh/ngsperl/).

`git clone https://github.com/shengqh/ngsperl/`

Assume you download the NGSPERL package to "/home/user/ngsperl", add "/home/user/ngsperl/lib" into your your perl library path.

`export PERL5LIB=/home/user/ngsperl/lib:$PERL5LIB`

NGSPERL just provide the interface to generate PBS scripts but not running the actual scripts. So, the tools used in the pipeline should be installed individually before the scripts being executed, for example, bwa for genome alignment.

<a name="framework"/>

# Framework

Our object oriented module-based framework includes three parts: modules, configurations, and a module parser. Each task in the pipeline will be implemented as a module. A configuration will be used for each specific research project which includes multiple tasks with user-defined parameters. A module parser will be used to parse the configuration to generate PBS scripts. Corresponding Linux shell scripts will also be generated for submitting the PBS scripts to the Linux cluster for lengthy tasks or running the PBS scripts directly for shorter tasks based on user choice.  

In order to allow the module integration, each module must implement three functions: result, perform, and get_pbs_files. 

The function result will provide the expected result files as if the task were executed.
 
The function perform is used to generate the corresponding PBS scripts based on user-defined parameters. Once two tasks are joined together in configuration, the downstream task will take the expected result files of the upstream task as the input parameters. The results from multiple tasks can also be used as input files in downstream tasks which makes the framework more flexible.
The function get_pbs_files will return the file names of the corresponding PBS scripts that function perform generated. With function get_pbs_files, the user can merge multiple shorter tasks into a single task through module "SequenceTask", then submit the PBS script of this single task to the Linux cluster. This feature is very useful for a pipeline with many quick tasks, such as smallRNA sequencing data analysis. 

Once the required tools are implemented as modules, for each real project, a Perl configuration structure will be defined to join those tools together. Then a module parser will parse this configuration and generate individual PBS scripts for each task. A shell script will also be generated to help the user submit multiple PBS scripts to the cluster or execute those PBS scripts sequentially. 
  
<a name="example"/>

# Quick start

The following code indicates a configuration of the simplest differentially expressed gene comparison. A tophat2 task is followed by a cuffdiff task. After executing the configuration script, the PBS scripts for both tophat2 and cuffdiff tasks will be generated under the user defined directory. Then the user can submit the generated tophat2 PBS scripts to the cluster first and submit cuffdiff PBS scripts after the tophat2 tasks are finished and the tophat2 results are validated.

	#!/usr/bin/perl
	use strict;
	use warnings;
	use CQS::ClassFactory;
	use CQS::FileUtils;
	my $target_dir           = create_directory_or_die("pipeline2");
	my $fasta_file           = "bowtie2_index/mm10.fa";
	my $bowtie2_index        = "bowtie2_index/mm10";
	my $transcript_gtf       = "Mus_musculus.GRCm38.68_chr1-22-X-Y-M.gtf";
	my $transcript_gtf_index = "Mus_musculus.GRCm38.68";
	my $email                = "quanhu.sheng\@vanderbilt.edu";
	my $task                 = "pipeline";
	my $config               = {
  		general => { task_name => $task },
  		files   => {
		  "S1" => ["rawdata/s1_sequence.txt"],
		  "S2" => ["rawdata/s2_sequence.txt"],
	    	  "S3" => ["rawdata/s3_sequence.txt"],
	    	  "S4" => ["rawdata/s4_sequence.txt"],
	    	  "S5" => ["rawdata/s5_sequence.txt"],
	    	  "S6" => ["rawdata/s6_sequence.txt"],
  		},
  		groups => {
	    	  "G1" => [ "S1", "S2", "S3" ],
		  "G2" => [ "S4", "S5", "S6" ],
  		},
	  	pairs   => { "G2_vs_G1" => [ "G1", "G2" ] },
  		tophat2 => {
		    class                => "Alignment::Tophat2",
    		    perform              => 1,
		    target_dir           => "${target_dir}/tophat2",
		    option               => "--segment-length 25 -r 0 -p 6",
		    source_ref           => "files",
		    bowtie2_index        => $bowtie2_index,
		    transcript_gtf       => $transcript_gtf,
		    transcript_gtf_index => $transcript_gtf_index,
		    rename_bam           => 1,
		    sh_direct            => 1,
		    pbs                  => {
    			"email"    => $email,
    	  		"nodes"    => "1:ppn=6",
    	  		"walltime" => "72",
    	  		"mem"      => "30gb"
    		    },
  		},
  		cuffdiff => {
    		class          => "Cufflinks::Cuffdiff",
    		perform        => 1,
    		target_dir     => "${target_dir}/cuffdiff",
    		option         => "-p 8 -u -N",
    		transcript_gtf => $transcript_gtf,
    		source_ref     => "tophat2",
    		groups_ref     => "groups",
    		pairs_ref      => "pairs",
    		fasta_file     => $fasta_file,
    		sh_direct      => 1,
    		pbs            => {
    	  		"email"    => $email,
			"nodes"    => "1:ppn=8",
    	  		"walltime" => "720",
    	  		"mem"      => "40gb"
    		},
  	    },
	};
	performConfig($config);

<a name="pipelines"/>

# Pipelines

We initialized a few pipelines:

|Pipeline|Example|Software|Description|
|---|---|---|---|
|RNASeq|[simple](https://github.com/shengqh/ngsperl/blob/master/examples/pipeline_rnaseq_simple.pl)|||
||[advanced](https://github.com/shengqh/ngsperl/blob/master/examples/pipeline_rnaseq_advanced.pl)|||
|ExomeSeq|[simple](https://github.com/shengqh/ngsperl/blob/master/examples/pipeline_exomeseq_simple.pl)|||
||[advanced](https://github.com/shengqh/ngsperl/blob/master/examples/pipeline_exomeseq_advanced.pl)|||
|ChIPSeq|[simple](https://github.com/shengqh/ngsperl/blob/master/examples/pipeline_chipseq_simple.pl)|||
||[advanced](https://github.com/shengqh/ngsperl/blob/master/examples/pipeline_chipseq_advanced.pl)|||

<a name="module"/>

# Implemented Modules

|Category|Module|Software|Description|
|---|---|---|---|
|Preprocessing|Format::Demultiplex|cqstools fastq_demultiplex|Demultiplexing the fastq file based on barcodes|
||CQS::FastqTrimmer|cqstools fastq_trimmer|Trimming 'N' from both 5' and 3' terminal of reads|
||Trimmer::Cutadapt|cutadapt|Removing adapter sequences from reads|
|Mapping|Alignment::BWA|bwa|Bwa genome mapping algorithm|
||Alignment::Bowtie1|bowtie1|Bowtie1 genome mapping algorithm|
||Alignment::Bowtie2|bowtie2|Bowtie2 genome mapping algorithm|
||Alignment::Tophat2|tophat2|Tophat2 RNAseq data assembler|
|Refinement|GATK::Refine|gatk|Realignment, base calibration and removing duplication|
|QC|QC::FastQC|fastqc|Quality control of fastq file|
||QC::RNASeQC|RNASeQC|Quality control of bam file|
|Count|Count::HTSeqCount|HTSeq-count|Counting gene reads|
||Count::DexseqCount|DEXSeq|Count exon reads|
|Summarize|CQS::CQSDatatable|cqstools data_table|Build count table from multiple counting result|
|Comparison|Comparison::DESeq2|DESeq2|Differential expression comparison of count data|
||Cufflinks::Cufflinks|cufflinks|Transcript assembly, differential expression, and differential regulation for RNA-Seq|
||Cufflinks::Cuffmerge|cuffmerge||
||Cufflinks::Cuffdiff|cuffdiff||
|Variants|GATK::MuTect|mutect|Somatic mutation caller|
||GATK::SNPIndel|gatk|SNP, indel caller|
||VarScan2::Mpileup2snp|Varscan2|SNP caller|
||VarScan2::Somatic|Varscan2|Somatic mutation caller|
|Annotation|Annotation::Annovar|annovar|Annotating SNP, indel and somatic mutation|
