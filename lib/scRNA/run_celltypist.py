import argparse
import logging
import os
import sys
import celltypist
import scanpy as sc


def build_parser():
  parser = argparse.ArgumentParser(
    description="Run CellTypist on an h5ad file and write per-cell predictions.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument("-i", "--input", required=True, help="Input h5ad file")
  parser.add_argument("--model_file", required=True, help="CellTypist model file (.pkl)")
  parser.add_argument(
    "--output_prefix",
    required=True,
    help="Output prefix. The script writes <output_prefix>.meta.csv and <output_prefix>.log",
  )
  parser.add_argument(
    "--majority_voting",
    action="store_true",
    help="Enable CellTypist majority voting refinement in addition to raw predictions",
  )
  return parser


def initialize_logger(logfile):
  logger = logging.getLogger("celltypist")
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


def celltypist_predict(args, output_file, logger):
  logger.info("reading input %s", args.input)
  adata = sc.read_h5ad(args.input)

  # 1. Store raw counts in .layers if you haven't already (good practice)
  adata.layers["counts"] = adata.X.copy()

  # 2. Normalize each cell to 10,000 counts
  sc.pp.normalize_total(adata, target_sum=1e4)

  # 3. Logarithmize the data (log(x+1))
  sc.pp.log1p(adata)

  logger.info(f"run celltypist using model {args.model_file} with majority_voting={args.majority_voting} ...")
  predictions = celltypist.annotate(adata, model = args.model_file, majority_voting = args.majority_voting)
  
  prediction_df = predictions.predicted_labels.copy()
  if getattr(prediction_df, "ndim", 1) == 1:
    prediction_df = prediction_df.to_frame(name="predicted_labels")

  prediction_df.index.name = "cell"

  output_dir = os.path.dirname(output_file)
  if output_dir:
    os.makedirs(output_dir, exist_ok=True)

  logger.info("writing %s", output_file)
  prediction_df.to_csv(output_file, index=True, index_label="cell")
  logger.info("done")


def main():
  parser = build_parser()

  if len(sys.argv) == 1:
    parser.print_help()
    return 1

  args = parser.parse_args()
  check_file(args.input, parser, "input")
  check_file(args.model_file, parser, "model")

  output_file, log_file = get_output_paths(args.output_prefix)
  logger = initialize_logger(log_file)

  try:
    celltypist_predict(args, output_file, logger)
  except Exception:
    logger.exception("CellTypist annotation failed")
    return 1

  return 0


if __name__ == "__main__":
  sys.exit(main())
