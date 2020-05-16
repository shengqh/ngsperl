"""
purpose= recalibrate the bam file based on the known SV/variant
result = '{pw}/{step}/result/{lb}.recal.markdup.bam'
"""
import os
import re
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('-pw', help="""working path""")
ps.add_argument('-bam', help="""the input bam file""")
ps.add_argument('-gm_fa', help="""genome fasta file""")
ps.add_argument('-mills_indel', help="""mills indel file for recalibration""")
ps.add_argument('-dbsnp', help="""dbsnp file for recalibration""")
ps.add_argument('-known_indel', help="""optional, known_indel file for recalibration""", nargs='?')
ps.add_argument('-pw_tmp', help="""optional, the path for the temp file, default is under system temp folder""", nargs='?', default=None)
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

step = 'bam_refine'
pw = args.pw
tmp_root = args.pw_tmp or '/temp'
bam_in = args.bam
genome = args.gm_fa  # gm, dbsnp, mills_indel, known_indel
mills = args.mills_indel
dbsnp = args.dbsnp
known_indel = args.known_indel

dock_file = '/scratch/cqs/chenh19/dock/centos.sif'
dock_prefix = f'singularity exec {dock_file} '
lb = re.sub(".(gz|fastq|txt|sam|bam|fastq.gz|fq.gz|tsv)$", '', bam_in.rsplit('/', 1)[-1])

# mkdir
pwtmp = f'{tmp_root}/gatk'
os.system(f'mkdir -p {pwtmp} 2>/dev/null')
os.system(f'mkdir -p {pw}/{step}/result 2>/dev/null')


# the known indel file is optional
# bacause if the build is hg19, this file is not available
known_sites_indel = f'-knownSites {known_indel}' if known_indel else ''
known_indel = f'-known {known_indel}' if known_indel else ''

# filename
f_recal_indel_interval = f'{pw}/{step}/result/{lb}.indel.intervals'
f_recal_snp_table = f'{pw}/{step}/result/{lb}.BQSR.table'
bam_recal_indel = f'{pwtmp}/{lb}.recal.indel.bam'
bam_recal_snp = f'{pwtmp}/{lb}.recal.dbsnp.bam'
bam_markdup = f'{pw}/{step}/result/{lb}.recal.markdup.bam'

# build pbs.sh
shell_file = f'{pw}/{step}/tmp_{lb}_{step}.sh'
with open(shell_file, 'w') as out:
    print(f"""#!/usr/bin/env bash
# realign indel
echo "gatk RealignerTargetCreator"
date
if [ ! -s {f_recal_indel_interval} ];then
gatk3 -Xmx2g -T RealignerTargetCreator \\
    -R {genome} \\
    -I {bam_in} \\
    -known {mills} {known_indel} \\
    -o {f_recal_indel_interval}
else
echo {lb}.indel.intervals already exist {f_recal_indel_interval}
fi

echo -e "\\n\\ngatk IndelRealigner"
date
if [ ! -s {bam_recal_indel} ];then
gatk3 -Xmx4g -T IndelRealigner \\
    -R {genome} \\
    -I {bam_in} \\
    -targetIntervals {f_recal_indel_interval} \\
    -known {mills} {known_indel} \\
    -o {bam_recal_indel}
else
echo file already exist {bam_recal_indel}
fi

# base recal BQSR  base quality scores recalibrate
echo -e "\\n\\ngatk BaseRecalibrator"
date

if [ ! -s {f_recal_snp_table} ];then
gatk3 -Xmx4g -T BaseRecalibrator -R {genome} -I {bam_recal_indel} -knownSites {dbsnp} -knownSites {mills} {known_sites_indel} -o  {f_recal_snp_table}
else
    echo file already exist {f_recal_snp_table}
fi


echo -e "\\n\\ngatk PrintReads"
date
if [ ! -s {bam_recal_snp} ];then
    gatk3 -Xmx2g -T PrintReads -R {genome} -I {bam_recal_indel} --BQSR {f_recal_snp_table} -o {bam_recal_snp}
else
    echo file already exist {bam_recal_snp}
fi

# MarkDuplicates
echo "picard MarkDuplicates"
date
if [ ! -s {bam_markdup} ];then
picard -Xmx2g MarkDuplicates I={bam_recal_snp} O={bam_markdup} M={pw}/{step}/result/{lb}.markdup.metrics.txt
else
echo file already exist {bam_markdup}
fi

# cleaning
# rm {bam_recal_indel}* 2>/dev/null
# rm {bam_recal_snp}* 2>/dev/null


    """, file=out)

# run this step
os.system(f'{dock_prefix} bash {shell_file}')
