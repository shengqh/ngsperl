#!/usr/bin/env python3
"""
get subregional bam
input = bam file / folder  + region(string or bed file)
"""

import sys
import os
import re
# import glob
import argparse as arg
ps = arg.ArgumentParser(description=__doc__)

ps.add_argument(
    '-region',
    help='region to substract, format=chr1:123-456:regionlb (1-indexed, label at the end could be absent),or bed file(0-indexed)')
ps.add_argument('-bam', help='bam file or folder, default=all bam file under this folder', nargs='?', default='.')
ps.add_argument('-lbcol', help='optinal, valid if input is bed file, 1-indexed. column number for the region label',
                nargs='?', type=int, default=-999)
ps.add_argument('-pad', help='extend left and right', nargs='?', type=int, default=-999)
ps.add_argument(
    '-pw', '-out', dest='pw',
    help='the path for the subregion result, if not exist, would create one, defaultname = "subregioinal_bam"',
    nargs='?', default="subregional_bam")
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


region = args.region
chrmap = {1: '1', 2: '2', 3: '3', 4: '4', 5: '5', 6: '6', 7: '7', 8: '8', 9: '9', 10: '10', 11: '11', 12: '12',
          13: '13', 14: '14', 15: '15', 16: '16', 17: '17', 18: '18', 19: '19', 20: '20', 21: '21', 22: '22', 23: 'X',
          24: 'Y', 25: '25', 26: '26', '1': '1', '2': '2', '3': '3', '4': '4', '5': '5', '6': '6', '7': '7', '8': '8',
          '9': '9', '10': '10', '11': '11', '12': '12', '13': '13', '14': '14', '15': '15', '16': '16', '17': '17',
          '18': '18', '19': '19', '20': '20', '21': '21', '22': '22', '23': 'X', '24': 'Y', '25': '25', '26': '26',
          'chr1': '1', 'chr2': '2', 'chr3': '3', 'chr4': '4', 'chr5': '5', 'chr6': '6', 'chr7': '7', 'chr8': '8',
          'chr9': '9', 'chr10': '10', 'chr11': '11', 'chr12': '12', 'chr13': '13', 'chr14': '14', 'chr15': '15',
          'chr16': '16', 'chr17': '17', 'chr18': '18', 'chr19': '19', 'chr20': '20', 'chr21': '21', 'chr22': '22',
          'chr23': 'X', 'chr24': 'Y', 'chr25': '25', 'chr26': '26', 'chr01': '1', 'chr02': '2', 'chr03': '3',
          'chr04': '4', 'chr05': '5', 'chr06': '6', 'chr07': '7', 'chr08': '8', 'chr09': '9', '01': '1', '02': '2',
          '03': '3', '04': '4', '05': '5', '06': '6', '07': '7', '08': '8', '09': '9', 'x': 'X', 'y': 'Y', 'xy': '25',
          'mt': '26', '0x': 'X', '0y': 'Y', 'X': 'X', 'Y': 'Y', 'XY': '25', 'MT': '26', '0X': 'X', '0Y': 'Y', 'chrx':
          'X', 'chry': 'Y', 'chrxy': '25', 'chrmt': '26', 'chr0x': 'X', 'chr0y': 'Y', 'chrX': 'X', 'chrY': 'Y',
          'chrXY': '25', 'chrMT': '26', 'chr0X': 'X', 'chr0Y': 'Y', 'chrM': '25', 'mito': '25', 'M': '25', 'm': '25',
          'chrm': '25', 'Mito': '25'}

# get bam list
bam = args.bam

if os.path.isfile(bam):
    bamlist = [bam]
elif os.path.isdir(bam):
    bam = os.path.abspath(bam)
    bamlist = os.popen(f'find {bam} -iname "*.bam"').read().split('\n')
    bamlist = [_.strip() for _ in bamlist if _.strip()]
else:
    print('no bam found: ', bam)
    sys.exit()

if len(bamlist) < 1:
    print('no bam found under folder:', bam)
    sys.exit()

pw = args.pw
if pw[0] == '/' and pw.split('/') == 2:
    print('you cannot export the result to the root!')
    sys.exit()

os.system(f'mkdir {pw} 2>/dev/null')


if os.path.isfile(region):
    region_list_raw = [re.split(r'\s+', _.strip()) for _ in open(region)]
    region_list = []
    for i in region_list_raw:
        if args.lbcol != -999:
            try:
                lb = i[args.lbcol - 1]
            except:
                print(f'wrong column number for region label (-lbcol={args.lbcol}), first row is {i}')
                sys.exit()
        else:
            lb = ''
        try:
            i[1] = int(i[1])
            i[2] = int(i[2])
            if args.pad != -999:
                i[1] -= args.pad
                i[2] += args.pad
            i[0] = 'chr' + chrmap[i[0]]
            lb = lb or f'{i[0]}_{i[1]}_{i[2]}'
            lb = f'region_{lb}' if re.match(r'^\d+$', lb) else lb
            region_list.append([i[0], i[1], i[2], lb])
        except:
            print('wrong chr/coordinate in bed file:', i)
            sys.exit()

else:
    m = re.split('[:_,-]', region)
    if len(m) < 3:
        print('bad regioin format, should be like  chr1:123-456')
        print(m)
        sys.exit()
    else:
        try:
            m[0] = 'chr' + chrmap[m[0]]
            m[1] = int(m[1]) - 1
            m[2] = int(m[2]) - 1
            if args.pad != -999:
                m[1] -= args.pad
                m[2] += args.pad

        except:
            print('wrong chr/coordinate in bed file:', m)
            sys.exit()
        try:
            lb = m[3]
        except:
            lb = f'{m[0]}_{m[1]}_{m[2]}'
        lb = f'region_{lb}' if re.match(r'^\d+$', lb) else lb
        region_list = [m[:3] + [lb]]


for ibam in bamlist:
    n = 0
    bam_sn = ibam.rsplit('/', 1)[-1].rsplit('.', 1)[0]
    for irg in region_list:
        n += 1
        lb_region = irg[3]
        bam_out = f'{pw}/{lb_region}_{bam_sn}.bam'
        os.system(f'samtools view -bS {ibam} "{irg[0]}:{irg[1]}-{irg[2]}" -o {bam_out}')
        os.system(f'sambamba index {bam_out} 2>/dev/null')
        print(f'{ibam}\t{irg}')
