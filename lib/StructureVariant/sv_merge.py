"""
this script is for merging the SV call results across samples into a single vcf file
would genotype the bam file of each sample
"""

import os

import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('-pw', help="""working path""")
ps.add_argument('-out', help="""the cohort name for output""")
ps.add_argument('-gm_build', help="""gm build, used for choose fasta file and the exclude bed file""",
choices=['hg19', 'hg38'])
args = ps.parse_args()

# load logging
import logging
import logging.config
import yaml
home = os.path.expanduser('~')
fl_log_conf = f'{home}/jb/config/logging_setting.yaml'
log_prefix = 'sv_call_pipeline'
log_prefix = log_prefix + '_' if log_prefix else ''
with open(fl_log_conf) as f:
    cfg = yaml.safe_load(f.read())
    cfg['handlers']['file']['filename'] = log_prefix + cfg['handlers']['file']['filename']
    cfg['handlers']['error']['filename'] = log_prefix + cfg['handlers']['error']['filename']
    logging.config.dictConfig(cfg)
logger = logging.getLogger('main')

logger.debug(f'input args={args}')

pw = args.pw
cohort = args.out
genome_build = args.gm_build
dock_file = '/scratch/cqs/chenh19/dock/centos.sif'
dock_prefix = f'singularity exec dock_file '

# make dir
os.system(f'mkdir -p {pw}/sv_genotyped/ 2>/dev/null')
os.system(f'mkdir -p {pw}/sv_genotyped_merge/ 2>/dev/null')
os.system(f'mkdir -p {pw}/sv_annotate/ 2>/dev/null')

# the the reference file
if genome_build == 'hg19':
    genome = '/scratch/h_vangard_1/chenh19/ref/hg19/Homo_sapiens_assembly19.fasta'
    gff3 = '/scratch/h_vangard_1/chenh19/ref/hg19/Homo_sapiens.GRCh37.87.gff3.gz'
elif genome_build == 'hg38':
    genome = '/scratch/h_vangard_1/chenh19/ref/hg38/Homo_sapiens_assembly38.fasta'
    gff3 = '/scratch/h_vangard_1/chenh19/ref/hg38/Homo_sapiens.GRCh38.99.gff3.gz'

# merge all the SV sites
# output = merged.sites.vcf.gz
os.system(f"{dock_prefix} smoove merge --name merged -f {genome} --outdir {pw}/call_sv {pw}/call_sv/result/*.genotyped.vcf.gz")

# genotype the sampleat the called SV sites
# parallel run the tasks
thread = 5
os.system(f"""ls {pw}/bam_merge/result/*.merge.bam | parallel -j  {thread} 'bamname={{}}; lb=${{bamname%.merge.bam}}; {dock_prefix} smoove genotype -d -x -p 1 --name $lb-joint --outdir {pw}/sv_genotyped/ --fasta {genome} --vcf merged.sites.vcf.gz {{}}""")


# merge all the samples together
os.system(f'{dock_prefix} smoove paste --name {cohort} --outdir {pw}/sv_genotyped_merge/ sv_genotyped/*.vcf.gz')

# annotate the final vcf with exons, UTRs
os.system(f'{dock_prefix} smoove annotate --gff {gff3} {pw}/sv_genotyped_merge/{cohort}.smoove.square.vcf.gz| bgzip -c > {cohort}.smoove.square.anno.vcf.gz')
