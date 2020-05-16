"""
input = final bam file of a single sample file
output = $pw/call_sv/result/$sample-smoove.genotyped.vcf.gz
purpose = get the SV sites
https://github.com/brentp/smoove/blob/master/README.md
"""
import os
import re
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('-pw', help="""working path""")
ps.add_argument('-bam', help="""bam file for SV call""")
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

step = 'call_sv'
genome_build = args.gm_build
pw = args.pw
bam_in = args.bam
dock_file = '/scratch/cqs/chenh19/dock/centos.sif'
dock_prefix = f'singularity exec {dock_file} '
lb = re.sub(".(gz|fastq|txt|sam|bam|fastq.gz|fq.gz|tsv)$", '', bam_in.rsplit('/', 1)[-1])

# make dir
os.system(f'mkdir -p {pw}/{step}/result 2>/dev/null')

# the the reference file
if genome_build == 'hg19':
    genome = '/scratch/h_vangard_1/chenh19/ref/hg19/Homo_sapiens_assembly19.fasta'
    excludebed = '/scratch/h_vangard_1/chenh19/ref/svcall/ceph18.b37.lumpy.exclude.2014-01-15.bed'
elif genome_build == 'hg38':
    genome = '/scratch/h_vangard_1/chenh19/ref/hg38/Homo_sapiens_assembly38.fasta'
    excludebed = '/scratch/h_vangard_1/chenh19/ref/svcall/exclude.cnvnator_100bp.GRCh38.20170403.bed'


os.system(f'{dock_prefix} smoove call --outdir {pw}/{step}/result --exclude {excludebed} --name {lb} --fasta {genome} -p 1 --excludechroms "hs37d5,~:,~^GL,~decoy,~alt" --genotype {bam_in}')
