import argparse
import json
from pathlib import Path

ENZYMES = ("HindIII", "DpnII", "MboI", "none")

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input fastq files",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output json file",
    )
    parser.add_argument(
        "-c",
        "--config_file",
        required=True,
        help="Input config file",
    )
    return parser

def get_input_json(fastqs_str, config_file):
    options = {}
    with open(config_file, "rt") as fin:
      for line in fin:
        parts = line.rstrip().split('\t')
        options[parts[1]] = parts[0]

    fastqs_list = fastqs_str.split(',')
    fastqs = []
    i = 0
    while(i < len(fastqs_list)):
      f1 = fastqs_list[i]
      f2 = fastqs_list[i+1]
      i+=2
      fastqs.append([{"read_1":f1, "read_2":f2}])
    
    input_json = {
        "hic.fastq": fastqs,
        "hic.assembly_name": options["hic.assembly_name"],
        "hic.chrsz": options["hic.chrsz"],
        "hic.reference_index": options["hic.reference_index"],
        "hic.align_num_cpus": options["hic.align_num_cpus"]
    }

    restriction_enzymes = options["hic.restriction_enzymes"].split(',')
    input_json["hic.restriction_enzymes"] = restriction_enzymes

    if options["hic.restriction_enzymes"] != "none":
      input_json["hic.restriction_sites"] = options["hic.restriction_sites"]

    return input_json

def write_json_to_file(data, outfile):
    Path(outfile).write_text(json.dumps(data, indent=2, sort_keys=True))

def main():
    parser = get_parser()
    args = parser.parse_args()

    input_json = get_input_json(
        fastqs_str=args.input,
        config_file=args.config_file
    )

    write_json_to_file(input_json, args.output)

if __name__ == "__main__":
    main()
