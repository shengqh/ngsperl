"""
dependency = curl,  parallel, fastq-dump
purpose: ge the corresponding fastq file for the SRR
input = SRRxxxx
output = [{pw}/sra2fastq/result/SRRxxx/SRRxxx_1.fastq.gz, {pw}/sra2fastq/result/SRRxxx/SRRxxx_2.fastq.gz]
"""
import os
import sys
from subprocess import Popen, PIPE
import re

import argparse as arg
from argparse import RawTextHelpFormatter
from common_utils import checkFileExists, runCmd, initializeLogger

def check_fastq_integrity(file_local, srr, ftp_folder):
    """
    return 1 if the file is complete
    return 0 if the filesize is not the same as EBI
    """

    size_local = os.popen(f"ls -l {file_local} | awk '{{print $5}}'").read().strip()
    size_ebi = os.popen(f"""curl {ftp_folder}|awk '{{print $5}}'""").read().strip()

    if size_local != size_ebi:
        logger.warning(f'{ifl}: incomplete fastq file. EBI size={size_ebi}, local size={size_local}')
        return 0
    else:
        return 1

def sra_ebi_query(logger, sra_acc):
    """
    check if the fastq file exist on the EBI server
    """
    part1 = sra_acc[:6]
    part2 = f'00{sra_acc[-1]}'
    ftp_folder = f'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{part1}/{part2}/{sra_acc}/'
    flist, stderr = Popen(['curl', '-l', ftp_folder], stderr=PIPE, stdout=PIPE, encoding='utf-8').communicate()
    flist = [_.strip() for _ in flist.split('\n') if _.strip()]
    # e.g  ['SRR8670693_1.fastq.gz', 'SRR8670693_2.fastq.gz']

    if stderr.find('Access failed') > -1:
        logger.warning(f'{sra_acc}: no fastq found on EBI server: would run fastq-dump locally')
        return 1
    elif len(flist) not in [1, 2]:
        logger.warning(f'{sra_acc}: wrong fastq count on EBI n = {len(flist)}, would run fastq-dump locally')
        return 1
    else:
        return flist

def download(logger, srr_raw, is_single_end, output_folder): 
    paired = not is_single_end

    # validate the SRR
    m = re.match(r'.*([SED]RR\d+)', srr_raw)
    if not m:
        logger.fatal(f'Input {srr_raw} is not a valid srr ID, please have a check.')
        sys.exit(1)
    srr = m.group(1)

    part1 = srr[:6]
    part2 = f'00{srr[-1]}'
    ftp_folder = f'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{part1}/{part2}/{srr}/'
    logger.debug(f'srr={srr}, ftp folder={ftp_folder}')

    # check availablity of the sra fastq on EBI
    res = sra_ebi_query(logger, srr)
    res_new = []
    logger.debug(f'EBI fastq query result={res}')

    if res == 1:   # fail to get fastq file on EBI server
        logger.info('now running fastq dump')
        split = '--split-e' if paired else ''
        runCmd(logger, f'fastq-dump --gzip --helicos {split} --origfmt -O {output_folder} {srr}')
    else:
        # first check if the file already exist
        for ifl in res:
            ifl_local = f'{output_folder}/{ifl}'
            if os.path.exists(ifl_local) and check_fastq_integrity(ifl_local, srr, ftp_folder):
                logger.info('file already exist and complete: {ifl_local}')
                # the file already exist, so don't download it
            else:
                res_new.append(ifl)

        # download the file
        url = ' '.join([f'{ftp_folder}{ifl}' for ifl in res_new])
        # logger.debug(url)
        # the actual download command
        runCmd(logger, f"cd {output_folder}; parallel 'wget -c {{}}' ::: {url} ")

        # verify the downloaded file
        fail = 0
        for ifl in res:
            valid = check_fastq_integrity(ifl, srr, ftp_folder)
            if valid:
                logger.debug(f'{ifl} is complete')
            else:
                fail = 1

        return fail

def main():
    DEBUG = False
    NOT_DEBUG = not DEBUG

    ps = arg.ArgumentParser(description="Download FASTQ files based on SRR ID", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ps.add_argument("-i", "--input", help="Input SRR ID", required=NOT_DEBUG)
    ps.add_argument("-s", "--is_single_end", action='store_true', help="Is the sample single-end data?")
    ps.add_argument("-o", "--output", help="Output file prefix", required=NOT_DEBUG)
    args = ps.parse_args()

    logger = initializeLogger(args.output + ".log", "download_fastq")
    logger.debug(f'input args={args}')

    download(logger, args.input, args.is_single_end, args.output)

if __name__ == "__main__":
    main()
