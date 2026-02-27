#!/usr/bin/env python3

#conda create -n spatial_311 python=3.11 pyproj geopandas numpy pandas tifffile rasterio shapely ipykernel

import argparse
import logging
import os
import geopandas as gpd
import numpy as np
import pandas as pd
import tifffile
from PIL import Image
from rasterio import features
from shapely import affinity


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


def extract_masked_crop(gdf, cell_id, padding, full_image, save_file_path=""):
    idx = gdf.index[gdf.cell_id == cell_id][0]

    # Get the specific geometry (not just the bounds)
    geom = gdf.geometry.iloc[idx]

    minx, miny, maxx, maxy = gdf.geometry.bounds.values[idx]

    # Convert coordinates to pixel indices (adding padding)
    left, top = int(minx) - padding, int(miny) - padding
    right, bottom = int(maxx) + padding, int(maxy) + padding

    crop = full_image[max(0, top):bottom, max(0, left):right]

    # Create a mask for the crop
    # Shift geometry coordinates to match the local crop coordinates
    shifted_geom = affinity.translate(geom, xoff=-left, yoff=-top)

    # Create a binary mask (1 inside the shape, 0 outside)
    mask = features.rasterize(
        [shifted_geom],
        out_shape=crop.shape[:2],
        all_touched=True,
        fill=0,
        default_value=1,
    )

    # Apply the mask to the array
    if len(crop.shape) == 3:
        masked_array = crop * mask[:, :, np.newaxis]
    else:
        masked_array = crop * mask

    crop_pil = Image.fromarray(masked_array).convert("RGB")

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
        "--nuclear_geojson",
        required=True,
        help="Path to nucleus_segmentations.geojson",
    )
    parser.add_argument(
        "--dhsr_tiff", 
        required=True, 
        help="Path to the source TIFF image")
    parser.add_argument(
        "--celltype_csv",
        required=True,
        help="CSV with cell_id and cell_type columns",
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
        default=0,
        help="Limit number of rows to process (0 = all)",
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

    if not os.path.exists(args.nuclear_geojson):
      raise Exception('File not found: ' + args.nuclear_geojson)

    if not os.path.exists(args.dhsr_tiff):
      raise Exception('File not found: ' + args.dhsr_tiff)

    if not os.path.exists(args.celltype_csv):
      raise Exception('File not found: ' + args.celltype_csv)

    logger.info("Reading cell geojson data...")
    cell_gdf = gpd.read_file(args.cell_geojson, engine=args.engine)

    logger.info("Reading nuclear geojson data...")
    nuclear_gdf = gpd.read_file(args.nuclear_geojson, engine=args.engine)
    
    logger.info("Reading image data...")
    full_image = tifffile.imread(args.dhsr_tiff)

    logger.info("Reading cell type information...")
    celltype_df = pd.read_csv(args.celltype_csv)

    total_rows = len(celltype_df)
    if args.limit and args.limit > 0:
        total_rows = min(total_rows, args.limit)

    for rowidx in range(total_rows):
        cellid = celltype_df.cell_id.values[rowidx]
        celltype = celltype_df.cell_type.values[rowidx].replace(" ", "_")

        cell_out_name = f"{args.output_prefix}{celltype}_{cellid}.cell.padding_{args.padding}.png"
        extract_crop(cell_gdf, cellid, args.padding, full_image, save_file_path=cell_out_name)

        nuclear_out_name = f"{args.output_prefix}{celltype}_{cellid}.nuclear.padding_{args.padding}.png"
        extract_crop(nuclear_gdf, cellid, args.padding, full_image, save_file_path=nuclear_out_name)

        cell_out_name = f"{args.output_prefix}{celltype}_{cellid}.cell.masked.padding_{args.padding}.png"
        extract_masked_crop(cell_gdf, cellid, args.padding, full_image, save_file_path=cell_out_name)

        nuclear_out_name = f"{args.output_prefix}{celltype}_{cellid}.nuclear.masked.padding_{args.padding}.png"
        extract_masked_crop(nuclear_gdf, cellid, args.padding, full_image, save_file_path=nuclear_out_name)

if __name__ == "__main__":
    main()
