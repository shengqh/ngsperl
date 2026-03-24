#!/usr/bin/env python3
"""
Extract cell images from a VisiumHD spatial dataset.

Usage:
  python extract_image.py \
    --geojson_path /path/to/cell_segmentations.geojson \
    --image_path /path/to/image.tiff \
    --cell_file /path/to/cells.txt \
    --output_dir /path/to/output \
    [--padding 0]
"""

import argparse
import os
import re
import geopandas as gpd
import numpy as np
import pandas as pd
from PIL import Image
import tifffile
from rasterio import features
from shapely import affinity


def extract_crop(gdf, cell_id, padding, full_image, save_file_path=""):
    idx = gdf.index[gdf.cell_id == cell_id][0]
    minx, miny, maxx, maxy = gdf.geometry.bounds.values[idx]
    print(f"Cell ID: {cell_id}, Bounds: ({minx}, {maxx}, {miny}, {maxy})")

    left, top = int(minx) - padding, int(miny) - padding
    right, bottom = int(maxx) + padding, int(maxy) + padding

    crop = full_image[max(0, top):bottom, max(0, left):right]
    crop_pil = Image.fromarray(crop).convert('RGB')

    if save_file_path:
        crop_pil.save(save_file_path)

    return crop_pil


def extract_masked_crop(gdf, cell_id, padding, full_image, save_file_path=""):
    idx = gdf.index[gdf.cell_id == cell_id][0]
    geom = gdf.geometry.iloc[idx]

    minx, miny, maxx, maxy = gdf.geometry.bounds.values[idx]

    left, top = int(minx) - padding, int(miny) - padding
    right, bottom = int(maxx) + padding, int(maxy) + padding

    crop = full_image[max(0, top):bottom, max(0, left):right]

    shifted_geom = affinity.translate(gdf.geometry.iloc[idx], xoff=-left, yoff=-top)

    mask = features.rasterize(
        [shifted_geom],
        out_shape=crop.shape[:2],
        all_touched=True,
        fill=0,
        default_value=1,
    )

    if len(crop.shape) == 3:
        masked_array = crop * mask[:, :, np.newaxis]
    else:
        masked_array = crop * mask

    crop_pil = Image.fromarray(masked_array).convert('RGB')

    if save_file_path:
        crop_pil.save(save_file_path)

    return crop_pil


def main():
    parser = argparse.ArgumentParser(description="Extract cell images from VisiumHD spatial data.")
    parser.add_argument("--geojson_path", required=True, help="Path to cell_segmentations.geojson")
    parser.add_argument("--image_path", required=True, help="Path to the TIFF image file")
    parser.add_argument("--cell_file", required=True, help="Path to text file listing cell IDs (one per line, no header)")
    parser.add_argument("--output_dir", required=True, help="Directory to write output PNG images")
    parser.add_argument("--padding", type=int, default=0, help="Pixel padding around each cell bounding box (default: 0)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Reading cell list...")
    cells = pd.read_csv(args.cell_file, header=None, sep="\t")
    cells.columns = ["cell"]
    cells['cell_id'] = cells['cell'].apply(
        lambda x: pd.to_numeric(re.sub(r"cellid_0+", "", x).replace("-1", ""))
    )
    print(f"  {len(cells)} cells loaded.")

    print("Reading geojson data...")
    gdf = gpd.read_file(args.geojson_path)

    print("Reading image data...")
    full_image = tifffile.imread(args.image_path)

    print("Extracting crops...")
    for rowidx in range(len(cells)):
        cellid = cells.cell_id.values[rowidx]
        out_path = os.path.join(args.output_dir, f"{cellid}.padding_{args.padding}.png")
        extract_crop(gdf, cellid, args.padding, full_image, save_file_path=out_path)

    print(f"Done. {len(cells)} images saved to {args.output_dir}")


if __name__ == "__main__":
    main()
