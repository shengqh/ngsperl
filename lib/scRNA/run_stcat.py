import argparse
import logging
import os
import sys
import STCAT
import scanpy as sc


def build_parser():
  parser = argparse.ArgumentParser(
    description="Run STCAT on an h5ad file and write per-cell predictions.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument("-i", "--input", required=True, help="Input h5ad file")
  parser.add_argument(
    "--output_prefix",
    required=True,
    help="Output prefix. The script writes <output_prefix>.meta.csv and <output_prefix>.log",
  )
  parser.add_argument(
    "--reduction",
    required=True,
    help="Reduction to use for neighbor search",
  )
  return parser


def initialize_logger(logfile):
  logger = logging.getLogger("STCAT")
  logger.setLevel(logging.INFO)
  logger.handlers = []
  logger.propagate = False

  formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)-8s - %(message)s")

  stream_handler = logging.StreamHandler()
  stream_handler.setLevel(logging.INFO)
  stream_handler.setFormatter(formatter)
  logger.addHandler(stream_handler)

  file_handler = logging.FileHandler(logfile, "w")
  file_handler.setLevel(logging.INFO)
  file_handler.setFormatter(formatter)
  logger.addHandler(file_handler)

  return logger


def check_file(filename, parser, label):
  if not os.path.isfile(filename):
    print(f"error: {label} file does not exist: {filename}", file=sys.stderr)
    parser.print_help()
    sys.exit(1)


def get_output_paths(output_prefix):
  if output_prefix.endswith(".meta.csv"):
    base_prefix = output_prefix[:-9]
  elif output_prefix.endswith(".csv"):
    base_prefix = output_prefix[:-4]
  else:
    base_prefix = output_prefix

  output_file = base_prefix + ".meta.csv"
  log_file = base_prefix + ".log"
  return output_file, log_file


def stcat_predict(args, reduction, output_file, logger):
  logger.info("reading input %s", args.input)
  adata = sc.read_h5ad(args.input)

  if "connectivities" not in adata.obsp:
    if reduction != "pca":
      adata.obsm.keys()

      use_rep = "X_" + reduction
      if use_rep not in adata.obsm:
        logger.error("reduction %s not found in input data", reduction)
        sys.exit(1)
      
      logger.info("computing neighbors using %s ...", use_rep)
      sc.pp.neighbors(
        adata,
        use_rep=use_rep,   # key step
        n_neighbors=15,
        metric="euclidean"
      )

      if "connectivities" not in adata.obsp:
        logger.error("neighbors not computed correctly, 'connectivities' not found in adata.obsp")
        sys.exit(1)

  logger.info("run STCAT ...")
  results = STCAT.STCAT(adata) # also a AnnData

  # get meta data columns
  print(results.obs.columns)

  # 'Prediction', 'Cluster', 'Uncertainty score'

  results.obs['Prediction'].value_counts()
  
  output_dir = os.path.dirname(output_file)
  if output_dir:
    os.makedirs(output_dir, exist_ok=True)

  logger.info("writing %s", output_file)
  results.obs[['Prediction', 'Cluster', 'Uncertainty score']].to_csv(output_file, index=True, index_label="cell")
  logger.info("done")


def main():
  parser = build_parser()

  if len(sys.argv) == 1:
    parser.print_help()
    return 1

  args = parser.parse_args()
  check_file(args.input, parser, "input")

  output_file, log_file = get_output_paths(args.output_prefix)
  logger = initialize_logger(log_file)

  try:
    stcat_predict(args, args.reduction, output_file, logger)
  except Exception:
    logger.exception("STCAT annotation failed")
    return 1

  return 0


if __name__ == "__main__":
  sys.exit(main())
