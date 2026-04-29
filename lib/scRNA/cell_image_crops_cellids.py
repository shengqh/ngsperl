#!/usr/bin/env python3

#conda create -n spatial_311 python=3.11 pyproj geopandas numpy pandas tifffile rasterio shapely ipykernel

import argparse
import logging
import os
import subprocess
import geopandas as gpd
import numpy as np
import pandas as pd
import tifffile
from PIL import Image

np.random.seed(20260427)


def render_rmarkdown_report(output_df, output_prefix, template_rmd, rscript="Rscript"):
    output_csv = f"{output_prefix}.image_report.csv"
    output_html = f"{output_prefix}.image_report.html"
    output_dir = os.path.dirname(os.path.abspath(output_html))

    output_df.to_csv(output_csv, index=False)

    render_code = "\n".join([
        "args <- commandArgs(trailingOnly = TRUE)",
        "template <- normalizePath(args[1], mustWork = TRUE)",
        "report_csv <- normalizePath(args[2], mustWork = TRUE)",
        "output_html <- normalizePath(args[3], winslash = '/', mustWork = FALSE)",
        "output_dir <- normalizePath(dirname(output_html), mustWork = TRUE)",
        "rmarkdown::render(",
        "  input = template,",
        "  output_file = basename(output_html),",
        "  output_dir = output_dir,",
        "  params = list(report_csv = report_csv),",
        "  knit_root_dir = dirname(report_csv),",
        "  envir = new.env(parent = globalenv()),",
        "  quiet = TRUE",
        ")",
    ])

    command = [
        rscript,
        "--vanilla",
        "-e",
        render_code,
        os.path.abspath(template_rmd),
        os.path.abspath(output_csv),
        os.path.abspath(output_html),
    ]

    result = subprocess.run(
        command,
        cwd=output_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        error_message = "\n".join(part for part in [result.stdout.strip(), result.stderr.strip()] if part)
        raise RuntimeError(f"Failed to render R Markdown report with {rscript}:\n{error_message}")

    return output_csv, output_html


def extract_crop(gdf, cell_id, padding, full_image, save_file_path=""):
    idx = gdf.index[gdf.cell_id == cell_id][0]
    minx, miny, maxx, maxy = gdf.geometry.bounds.values[idx]

    # Convert coordinates to pixel indices (adding padding)
    left, top = int(minx) - padding, int(miny) - padding
    right, bottom = int(maxx) + padding, int(maxy) + padding

    # Crop the numpy array: image[y_range, x_range]
    crop = full_image[max(0, top):bottom, max(0, left):right]
    crop_pil = Image.fromarray(crop).convert("RGB")

    if save_file_path:
        crop_pil.save(save_file_path)

    return crop_pil

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract cell crops from a VisiumHD image and segmentation geojson."
    )
    parser.add_argument(
        "--cell_geojson", 
        required=True, 
        help="Path to cell_segmentations.geojson")
    parser.add_argument(
        "--nucleus_geojson",
        help="Path to nucleus_segmentations.geojson (optional)",
    )
    parser.add_argument(
        "--dhsr_tiff", 
        required=True, 
        help="Path to the source TIFF image")
    parser.add_argument(
        "--cellid_csv",
        required=True,
        help="CSV file with cell_id in the first column and group/cell type in the second column",
    )
    parser.add_argument(
        "--output_prefix",
        required=True,
        help="Prefix for output filenames",
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=0,
        help="Padding value in pixels (default: 0)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=50,
        help="Maximum number of cell_ids to sample per group (default: 50, set to 0 for no limit)",
    )
    parser.add_argument(
        "--engine",
        default="pyogrio",
        help="GeoJSON reader engine (default: pyogrio)",
    )
    return parser.parse_args()


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    logger = logging.getLogger("cell_image_crops")

    args = parse_args()

    # make sure all files exists
    if not os.path.exists(args.cell_geojson):
        raise Exception('File not found: ' + args.cell_geojson)

    if args.nucleus_geojson and not os.path.exists(args.nucleus_geojson):
        raise Exception('File not found: ' + args.nucleus_geojson)

    if not os.path.exists(args.dhsr_tiff):
        raise Exception('File not found: ' + args.dhsr_tiff)

    if not os.path.exists(args.cellid_csv):
      raise Exception('File not found: ' + args.cellid_csv)

    logger.info("Reading cell geojson data...")
    cell_gdf = gpd.read_file(args.cell_geojson, engine=args.engine)

    nucleus_gdf = None
    if args.nucleus_geojson:
        logger.info("Reading nucleus geojson data...")
        nucleus_gdf = gpd.read_file(args.nucleus_geojson, engine=args.engine)
    
    logger.info("Reading image data...")
    full_image = tifffile.imread(args.dhsr_tiff)

    logger.info("Reading cell id information...")
    cell_df = pd.read_csv(args.cellid_csv)
    if cell_df.shape[1] < 2:
        raise ValueError("cellid_csv must contain at least two columns: cell_id and seurat_cell_type")

    cell_df = cell_df.iloc[:, :2].copy()
    cell_df.columns = ["cell_id", "cell_group"]

    if cell_df.empty:
        logger.info("No cell ids found in %s", args.cellid_csv)
        return

    selected_groups = []
    for group_name, group_df in cell_df.groupby("cell_group", sort=False):
        sample_size = len(group_df)
        if args.limit and args.limit > 0:
            sample_size = min(sample_size, args.limit)

        logger.info("Selecting %d of %d cell_ids for group %s", sample_size, len(group_df), group_name)
        if sample_size < len(group_df):
            selected_groups.append(group_df.sample(n=sample_size))
        else:
            selected_groups.append(group_df)

    selected_df = pd.concat(selected_groups, ignore_index=True)
    output_rows = []

    for _, row in selected_df.iterrows():
        cellgroup = row["cell_group"]
        cellid = row["cell_id"]

        cell_image_name = f"{args.output_prefix}_{cellid}.cell.padding_{args.padding}.png"
        cell_image = extract_crop(cell_gdf, cellid, args.padding, full_image, save_file_path=cell_image_name)

        output_row = {
            "cell_group": cellgroup,
            "cell_id": cellid,
            "cell_image": cell_image_name,
            "cell_width": cell_image.size[0],
            "cell_height": cell_image.size[1],
        }

        if nucleus_gdf is not None:
            nucleus_image_name = f"{args.output_prefix}_{cellid}.nucleus.padding_{args.padding}.png"
            nucleus_image = extract_crop(nucleus_gdf, cellid, args.padding, full_image, save_file_path=nucleus_image_name)
            output_row.update({
                "nucleus_image": nucleus_image_name,
                "nucleus_width": nucleus_image.size[0],
                "nucleus_height": nucleus_image.size[1],
            })

        output_rows.append(output_row)

    output_df = pd.DataFrame(output_rows)

    template_rmd = f"{os.path.abspath(__file__)}.rmd"
    if not os.path.exists(template_rmd):
        raise FileNotFoundError(f"R Markdown template not found: {template_rmd}")

    report_csv, output_html = render_rmarkdown_report(output_df, args.output_prefix, template_rmd)
    logger.info("Saved report data to %s", report_csv)
    logger.info("Saved image report to %s", output_html)

if __name__ == "__main__":
    main()
