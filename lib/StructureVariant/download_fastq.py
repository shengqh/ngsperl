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
import argparse
from common_utils import checkFileExists, runCmd, initializeLogger

def check_fastq_integrity(logger, file_local, srr, ftp_folder):
    """
    return 1 if the file is complete
    return 0 if the filesize is not the same as EBI
    """
    if not os.path.isfile(file_local):
        return 0

    size_local = os.popen(f"ls -l {file_local} | awk '{{print $5}}'").read().strip()
    size_ebi = os.popen(f"""curl {ftp_folder}|awk '{{print $5}}'""").read().strip()

    if size_local != size_ebi:
        logger.warning(f'{ifl}: incomplete fastq file. EBI size={size_ebi}, local size={size_local}')
        return 0
    else:
        return 1

def sra_ebi_query(logger, sra_acc, is_single_end):
    """
    check if the fastq file exist on the EBI server
    """
    part1 = sra_acc[:6]
    part2 = f'00{sra_acc[-1]}'
    ftp_folder = f'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{part1}/{part2}/{sra_acc}/'
    flist, stderr = Popen(['curl', '-l', ftp_folder], stderr=PIPE, stdout=PIPE, encoding='utf-8').communicate()

    flist = [file_name.strip() for file_name in flist.split('\n') if file_name.strip()]
    logger.debug(f"files in ftp server: {flist}")

    # e.g  ['SRR8670693_1.fastq.gz', 'SRR8670693_2.fastq.gz']

    if stderr.find('Access failed') > -1:
        logger.warning(f'{sra_acc}: no fastq found on EBI server: would run fastq-dump locally')
        return 1
    else:
        file_expected = [f"{sra_acc}.fastq.gz"] if is_single_end else [f"{sra_acc}_1.fastq.gz", f"{sra_acc}_2.fastq.gz"]
        if any([fe for fe in file_expected if fe not in flist]):
            logger.debug(f"expect file {file_expected} not found in ftp server: {flist}")
            return 1
        else:
            return file_expected

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
    res = sra_ebi_query(logger, srr, is_single_end)
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
            if os.path.isfile(ifl_local) and check_fastq_integrity(logger, ifl_local, srr, ftp_folder):
                logger.info(f'file already exist and complete: {ifl_local}')
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
        for ifl in res_new:
            ifl_local = f'{output_folder}/{ifl}'
            valid = check_fastq_integrity(logger, ifl_local, srr, ftp_folder)
            if valid:
                logger.debug(f'{ifl} is complete')
            else:
                fail = 1

        return fail

def main():
    DEBUG = True
    NOT_DEBUG = not DEBUG

    ps = argparse.ArgumentParser(description="Download FASTQ files based on SRR ID", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ps.add_argument("-i", "--input", help="Input SRR ID", required=NOT_DEBUG)
    ps.add_argument("-s", "--is_single_end", action='store_true', help="Is the sample single-end data?")
    ps.add_argument("-o", "--output", help="Output folder", required=NOT_DEBUG)
    args = ps.parse_args()

    if DEBUG:
        #args.input = "SRR4253621" #single end
        #args.is_single_end = False
        args.input = "SRR4032155"
        args.is_single_end = False
        args.output = "/home/shengq2/test"

    logger = initializeLogger(os.path.join(args.output, args.input + ".log"), "download_fastq")
    logger.debug(f'input args={args}')

    download(logger, args.input, args.is_single_end, args.output)

if __name__ == "__main__":
    main()
