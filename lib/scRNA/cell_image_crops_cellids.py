#!/usr/bin/env python3

#conda create -n spatial_311 python=3.11 pyproj geopandas numpy pandas tifffile rasterio shapely ipykernel

import argparse
import html
import logging
import os
import geopandas as gpd
import numpy as np
import pandas as pd
import tifffile
from PIL import Image

np.random.seed(20260427)


def format_cell_id(cell_id):
    if pd.isna(cell_id):
        return ""

    try:
        return str(int(cell_id))
    except (TypeError, ValueError):
        return str(cell_id)


def write_html_report(output_df, output_html):
    report_dir = os.path.dirname(os.path.abspath(output_html))
    image_specs = [
        ("Cell Images", "cell_image", "cell_width", "cell_height"),
    ]
    if "nucleus_image" in output_df.columns:
        image_specs.append(("Nucleus Images", "nucleus_image", "nucleus_width", "nucleus_height"))

    max_dimension = 0
    for _, _, width_col, height_col in image_specs:
        max_dimension = max(max_dimension, int(output_df[width_col].max()), int(output_df[height_col].max()))

    scale_factor = 1.0
    if max_dimension > 180:
        scale_factor = 180.0 / max_dimension

    html_lines = [
        "<!DOCTYPE html>",
        "<html lang=\"en\">",
        "<head>",
        "  <meta charset=\"utf-8\">",
        "  <title>Cell Image Report</title>",
        "  <style>",
        "    body { font-family: Arial, sans-serif; margin: 24px; color: #1f2933; }",
        "    h1, h2, h3 { margin-bottom: 8px; }",
        "    p { margin-top: 0; }",
        "    .group-section { margin-bottom: 32px; }",
        "    .grid-wrapper { overflow-x: auto; padding-bottom: 8px; }",
        "    .image-grid { display: grid; grid-template-columns: repeat(10, max-content); gap: 12px; align-items: start; }",
        "    .image-card { margin: 0; text-align: center; }",
        "    .image-card img { display: block; border: 1px solid #d2d6dc; background: #ffffff; }",
        "    .image-card figcaption { margin-top: 4px; font-size: 12px; white-space: nowrap; }",
        "  </style>",
        "</head>",
        "<body>",
        "  <h1>Cell Image Report</h1>",
        f"  <p>All images are displayed using the same scale factor: {scale_factor:.3f}x.</p>",
    ]

    for cell_group, group_df in output_df.groupby("cell_group", sort=True):
        html_lines.append(f"  <section class=\"group-section\"><h2>{html.escape(str(cell_group))}</h2>")
        html_lines.append(f"  <p>{len(group_df)} cells</p>")

        for section_title, image_col, width_col, height_col in image_specs:
            html_lines.append(f"  <h3>{html.escape(section_title)}</h3>")
            html_lines.append("  <div class=\"grid-wrapper\"><div class=\"image-grid\">")

            for _, row in group_df.iterrows():
                rel_path = os.path.relpath(row[image_col], report_dir)
                display_width = max(1, int(round(row[width_col] * scale_factor)))
                display_height = max(1, int(round(row[height_col] * scale_factor)))
                formatted_cell_id = format_cell_id(row["cell_id"])
                cell_id_label = html.escape(formatted_cell_id)
                alt_text = html.escape(f"{cell_group} {formatted_cell_id} {section_title}")
                src_path = html.escape(rel_path)
                html_lines.append(
                    "    <figure class=\"image-card\">"
                    f"<img src=\"{src_path}\" alt=\"{alt_text}\" width=\"{display_width}\" height=\"{display_height}\" loading=\"lazy\">"
                    f"<figcaption>{cell_id_label}</figcaption></figure>"
                )

            html_lines.append("  </div></div>")

        html_lines.append("  </section>")

    html_lines.extend(["</body>", "</html>"])

    with open(output_html, "w", encoding="utf-8") as out_handle:
        out_handle.write("\n".join(html_lines))

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
    manifest_columns = ["cell_group", "cell_id", "cell_image"]
    if "nucleus_image" in output_df.columns:
        manifest_columns.append("nucleus_image")
    manifest_df = output_df[manifest_columns]
    output_csv = f"{args.output_prefix}.image_names.csv"
    manifest_df.to_csv(output_csv, index=False)
    logger.info("Saved image manifest to %s", output_csv)

    output_html = f"{args.output_prefix}.image_report.html"
    write_html_report(output_df, output_html)
    logger.info("Saved image report to %s", output_html)

if __name__ == "__main__":
    main()
