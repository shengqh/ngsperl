#!/usr/bin/env python3
"""Extract per-sample nucleus distances and plot combined distributions."""

import argparse
import csv
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree


def nearest_nucleus_distances(gdf, id_column):
    """Return each nucleus's nearest neighbor, radius, and distance/radius ratio."""
    valid = gdf.loc[
        gdf.geometry.notna() & ~gdf.geometry.is_empty & gdf[id_column].notna(),
        [id_column, "geometry"],
    ].copy()
    areas = np.array([geometry.area for geometry in valid.geometry])
    valid = valid.loc[areas > 0]
    areas = areas[areas > 0]
    if len(valid) < 2:
        raise ValueError("At least two nuclei with IDs and positive-area geometries are required.")

    cell_ids = valid[id_column].astype(str).to_numpy()
    if len(set(cell_ids)) != len(cell_ids):
        raise ValueError(f"The '{id_column}' field contains duplicate cell IDs.")

    centroids = [geometry.centroid for geometry in valid.geometry]
    coordinates = np.array([(centroid.x, centroid.y) for centroid in centroids])
    _, neighbor_indices = cKDTree(coordinates).query(coordinates, k=2, workers=-1)

    row_indices = np.arange(len(cell_ids))
    nearest_indices = np.where(
        neighbor_indices[:, 0] == row_indices,
        neighbor_indices[:, 1],
        neighbor_indices[:, 0],
    )
    distances = np.linalg.norm(coordinates - coordinates[nearest_indices], axis=1)
    radii = np.sqrt(areas / np.pi)

    return [
        (
            cell_ids[index],
            cell_ids[neighbor_indices],
            float(distances[index]),
            float(radii[index]),
            float(distances[index] / radii[index]),
        )
        for index, neighbor_indices in enumerate(nearest_indices)
    ]


def mutual_nearest_neighbors(records):
    """Keep directed records whose nearest-neighbor relationship is reciprocal."""
    nearest_by_cell = {record[0]: record[1] for record in records}
    return [
        record
        for record in records
        if nearest_by_cell.get(record[1]) == record[0]
    ]


def remove_large_distance_outliers(records, outlier_factor):
    """Apply Tukey's upper fence to nearest-neighbor distances."""
    distances = np.array([record[2] for record in records])
    lower_quartile, upper_quartile = np.quantile(distances, [0.25, 0.75])
    cutoff = upper_quartile + outlier_factor * (upper_quartile - lower_quartile)
    return [record for record in records if record[2] <= cutoff], float(cutoff)


def read_file_map(file_map):
    """Read whitespace-delimited filepath/name rows without a header."""
    entries = []
    with file_map.open() as map_file:
        for line_number, line in enumerate(map_file, start=1):
            line = line.strip()
            if not line:
                continue
            fields = line.rsplit(maxsplit=1)
            if len(fields) != 2:
                raise ValueError(f"Invalid file-map row {line_number}: expected filepath and name")
            entries.append((Path(fields[0]), fields[1]))
    if not entries:
        raise ValueError(f"No input rows found in file map: {file_map}")
    return entries


def write_records(records, output_path):
    """Write extracted nucleus records as a CSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(
            ["cellid1", "cellid2", "distance", "nucleus_radius", "distance_to_radius_ratio"]
        )
        writer.writerows(records)


def plot_combined_histograms(sample_values, output_path, title, xlabel, color):
    """Plot one histogram per sample in a single height-scaled column."""
    sample_names = list(sample_values)
    all_values = [value for values in sample_values.values() for value in values]
    if all_values:
        spread = float(np.ptp(all_values))
        scale = max(float(np.max(np.abs(all_values))), 1.0)
        bins = 1 if spread <= 1e-9 * scale else np.histogram_bin_edges(all_values, bins="auto")
    else:
        bins = 1

    figure_height = max(3.2, 2.6 * len(sample_names) + 0.8)
    figure, axes = plt.subplots(
        len(sample_names),
        1,
        figsize=(10, figure_height),
        sharex=True,
        squeeze=False,
    )
    for axis, sample_name in zip(axes[:, 0], sample_names):
        values = sample_values[sample_name]
        if values:
            axis.hist(values, bins=bins, color=color, edgecolor="white")
        else:
            axis.text(0.5, 0.5, "No mutual nearest-neighbor pairs", ha="center", va="center")
            axis.set_yticks([])
        axis.set_ylabel("Count")
        axis.set_title(sample_name, loc="left", fontsize=10)
    axes[-1, 0].set_xlabel(xlabel)
    figure.suptitle(title, y=0.995)
    figure.tight_layout(rect=(0, 0, 1, 0.99))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(str(output_path), dpi=200)
    plt.close(figure)


def extract_sample(geojson_path, sample_name, output_dir, outlier_factor):
    """Calculate and save standard and mutual-neighbor records for one sample."""
    if not geojson_path.is_file():
        raise FileNotFoundError(f"GeoJSON for {sample_name} does not exist: {geojson_path}")
    gdf = gpd.read_file(geojson_path)
    id_column = next((name for name in ("cell_id", "cellid", "id") if name in gdf.columns), None)
    if id_column is None:
        raise ValueError(f"No cell ID property (cell_id, cellid, or id) in {geojson_path}")

    records, cutoff = remove_large_distance_outliers(
        nearest_nucleus_distances(gdf, id_column), outlier_factor
    )
    if not records:
        raise ValueError(f"Outlier filtering removed all nucleus distances for {sample_name}")
    mutual_records = mutual_nearest_neighbors(records)
    safe_name = Path(sample_name).name
    write_records(records, output_dir / f"{safe_name}_nearest_nucleus_distances.csv")
    write_records(mutual_records, output_dir / f"{safe_name}_mutual_nearest_neighbors.csv")
    return records, mutual_records, cutoff


def main():
    parser = argparse.ArgumentParser(
        description="Calculate per-sample nucleus distances and plot combined histograms."
    )
    parser.add_argument(
        "--file-map",
        type=Path,
        required=True,
        help="Headerless filepath/name map",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        required=True,
        help="Prefix for combined figure filenames; per-sample CSVs use sample-name prefixes in its directory",
    )
    parser.add_argument(
        "--outlier-factor",
        type=float,
        default=1.5,
        help="IQR multiplier for the per-sample upper distance fence (default: 1.5)",
    )
    args = parser.parse_args()
    if args.outlier_factor < 0:
        parser.error("--outlier-factor must be non-negative")
    output_prefix = args.output_prefix
    output_dir = output_prefix.parent
    entries = read_file_map(args.file_map)

    all_distances = {}
    all_ratios = {}
    mutual_distances = {}
    mutual_ratios = {}
    for geojson_path, sample_name in entries:
        records, mutual_records, cutoff = extract_sample(
            geojson_path, sample_name, output_dir, args.outlier_factor
        )
        all_distances[sample_name] = [record[2] for record in records]
        all_ratios[sample_name] = [record[4] for record in records]
        mutual_distances[sample_name] = [record[2] for record in mutual_records]
        mutual_ratios[sample_name] = [record[4] for record in mutual_records]
        print(
            f"{sample_name}: extracted {len(records)} nuclei and "
            f"{len(mutual_records)} mutual-neighbor rows; distance cutoff {cutoff:.6g}"
        )

    plot_combined_histograms(
        all_distances,
        output_prefix.with_name(f"{output_prefix.name}_combined_nearest_nucleus_distances.png"),
        "Nearest-nucleus distances by sample",
        "Nearest-neighbor centroid distance (GeoJSON coordinate units)",
        "#397a72",
    )
    plot_combined_histograms(
        all_ratios,
        output_prefix.with_name(f"{output_prefix.name}_combined_distance_radius_ratios.png"),
        "Nearest-nucleus distance-to-radius ratios by sample",
        "Nearest-neighbor distance / nucleus equivalent radius",
        "#b16a45",
    )
    plot_combined_histograms(
        mutual_distances,
        output_prefix.with_name(
            f"{output_prefix.name}_combined_mutual_nearest_neighbor_distances.png"
        ),
        "Mutual nearest-neighbor distances by sample",
        "Nearest-neighbor centroid distance (GeoJSON coordinate units)",
        "#397a72",
    )
    plot_combined_histograms(
        mutual_ratios,
        output_prefix.with_name(
            f"{output_prefix.name}_combined_mutual_nearest_neighbor_distance_radius_ratios.png"
        ),
        "Mutual nearest-neighbor distance-to-radius ratios by sample",
        "Nearest-neighbor distance / nucleus equivalent radius",
        "#b16a45",
    )
    print(f"Saved per-sample CSVs to {output_dir}")
    print(f"Saved four combined figures with prefix {output_prefix}")


if __name__ == "__main__":
    main()