"""
merge the split bam file
"""
import os
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('-pw', help="""working path""")
ps.add_argument('-out', help="""outname prefix for the result bam file""")
ps.add_argument('-bamlist', help="""bam file chunks, sep by space""", nargs='+')
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

step = 'bam_merge'
pw = args.pw
bam_list = args.bamlist
bam_list_str = ' '.join(bam_list)
dock_file = '/scratch/cqs/chenh19/dock/centos.sif'
dock_prefix = f'singularity exec dock_file '
lb = args.out
result = f'{pw}/{step}/result/{lb}.merge.bam'


logger.info('merge bam')
os.system(f"{dock_prefix} sambamba merge {result}  {bam_list_str}")
logger.info('merge bam completed')

logger.info('build index')
os.system(f'{dock_prefix} sambamba index {result}')
logger.info('build index completed')
