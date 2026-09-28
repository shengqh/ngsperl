#!/usr/bin/env python3
"""Measure gene expression from Visium HD bins covered by nuclei."""

import argparse
import csv
import logging
import math
from pathlib import Path

import geopandas as gpd
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


LOGGER = logging.getLogger("nucleus_cell_gene_expression_ratio")
BIN_RESOLUTION = "square_002um"
CELL_BARCODE_PATTERN = r"^cellid_(?P<cell_id>\d+)-1$"
OUTPUT_COLUMNS = [
    "sample_name",
    "rank",
    "cell_id",
    "cell_area",
    "nucleus_area",
    "nucleus_to_cell_area_ratio",
    "total_gene_expression",
    "nucleus_bin_gene_expression",
    "nucleus_expression_fraction",
    "nucleus_expression_percent",
    "bins_assigned_to_cell",
    "bins_covered_by_nucleus",
]


def read_cell_ids_and_areas(geojson_path):
    """Read unique cell IDs and polygon areas in GeoJSON coordinate units squared."""
    LOGGER.info("Reading cell IDs and polygon areas from %s", geojson_path)
    gdf = gpd.read_file(geojson_path)
    if "cell_id" not in gdf.columns:
        raise ValueError(f"Missing cell_id property in {geojson_path}")
    valid = gdf["cell_id"].notna() & gdf.geometry.notna() & ~gdf.geometry.is_empty
    gdf = gdf.loc[valid]
    ids = gdf["cell_id"].astype("int64").to_numpy()
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate cell_id values in {geojson_path}")
    areas = np.array([geometry.area for geometry in gdf.geometry], dtype=np.float64)
    if np.any(areas <= 0):
        invalid_id = int(ids[np.flatnonzero(areas <= 0)[0]])
        raise ValueError(f"Non-positive polygon area for cell_id {invalid_id} in {geojson_path}")
    ordering = np.argsort(ids)
    ids = ids[ordering]
    areas = areas[ordering]
    LOGGER.info("Read %d unique cell IDs and areas from %s", len(ids), geojson_path)
    return ids, areas


def read_bin_counts(matrix_path):
    """Sum all gene counts per sparse HDF5 barcode and return expressed bins."""
    LOGGER.info("Loading sparse matrix data and column pointers from %s", matrix_path)
    with h5py.File(matrix_path, "r") as h5_file:
        matrix = h5_file["matrix"]
        data = matrix["data"][:]
        indptr = matrix["indptr"][:]
        LOGGER.info(
            "Loaded %d stored gene counts across %d bin barcodes; summing counts per bin",
            len(data),
            len(indptr) - 1,
        )
        prefix_sum = np.empty(len(data) + 1, dtype=np.int64)
        prefix_sum[0] = 0
        np.cumsum(data, dtype=np.int64, out=prefix_sum[1:])
        total_counts = prefix_sum[indptr[1:]] - prefix_sum[indptr[:-1]]
        nonzero_columns = np.flatnonzero(total_counts > 0)
        LOGGER.info(
            "Found %d expressed bins; reading all barcodes sequentially to avoid slow indexed HDF5 reads",
            len(nonzero_columns),
        )
        all_barcodes = matrix["barcodes"][:]
        expressed_barcodes = all_barcodes[nonzero_columns]
        del all_barcodes

    LOGGER.info("Parsing coordinates from %d expressed-bin barcodes", len(expressed_barcodes))
    barcode_bytes = expressed_barcodes.view(np.uint8).reshape(
        len(expressed_barcodes), expressed_barcodes.dtype.itemsize
    )
    prefix = np.frombuffer(b"s_002um_", dtype=np.uint8)
    valid_format = (
        np.all(barcode_bytes[:, : len(prefix)] == prefix, axis=1)
        & (barcode_bytes[:, 13] == ord("_"))
        & (barcode_bytes[:, 19] == ord("-"))
        & (barcode_bytes[:, 20] == ord("1"))
    )
    row_digits = barcode_bytes[:, 8:13].astype(np.int32) - ord("0")
    column_digits = barcode_bytes[:, 14:19].astype(np.int32) - ord("0")
    valid_digits = (
        np.all((row_digits >= 0) & (row_digits <= 9), axis=1)
        & np.all((column_digits >= 0) & (column_digits <= 9), axis=1)
    )
    valid = valid_format & valid_digits
    if not np.all(valid):
        invalid_index = int(np.flatnonzero(~valid)[0])
        invalid_barcode = expressed_barcodes[invalid_index].decode(errors="replace")
        raise ValueError(f"Unexpected {BIN_RESOLUTION} barcode: {invalid_barcode}")

    digit_weights = np.array([10000, 1000, 100, 10, 1], dtype=np.int32)
    coordinates = np.empty((len(expressed_barcodes), 2), dtype=np.int32)
    coordinates[:, 0] = row_digits @ digit_weights
    coordinates[:, 1] = column_digits @ digit_weights

    counts = total_counts[nonzero_columns]
    del total_counts, prefix_sum, data, indptr, expressed_barcodes, barcode_bytes
    LOGGER.info("Finished reading counts and parsing coordinates for %d expressed bins", len(counts))
    return coordinates, counts


def accumulate_counts_by_cell(
    coordinates,
    counts,
    mapping_path,
    cell_ids,
    batch_size=10000,
):
    """Aggregate expressed-bin counts using Space Ranger's segmentation mapping."""
    mapping_file = pq.ParquetFile(mapping_path)
    row_count = mapping_file.metadata.num_rows
    grid_width = int(round(np.sqrt(row_count)))
    if grid_width * grid_width != row_count:
        raise ValueError(
            f"Expected a square {BIN_RESOLUTION} mapping grid, found {row_count} rows"
        )
    if np.any(coordinates < 0) or np.any(coordinates >= grid_width):
        raise ValueError("A bin barcode has row/column coordinates outside the mapping grid")

    linear_indices = coordinates[:, 0].astype(np.int64) * grid_width + coordinates[:, 1]
    sort_order = np.argsort(linear_indices)
    linear_indices = linear_indices[sort_order]
    counts = counts[sort_order]
    coordinates = coordinates[sort_order]

    group_lengths = np.array(
        [mapping_file.metadata.row_group(i).num_rows for i in range(mapping_file.metadata.num_row_groups)],
        dtype=np.int64,
    )
    group_offsets = np.concatenate(([0], np.cumsum(group_lengths)))
    group_starts = np.searchsorted(linear_indices, group_offsets[:-1], side="left")
    group_ends = np.searchsorted(linear_indices, group_offsets[1:], side="left")
    used_groups = np.flatnonzero(group_ends > group_starts)
    LOGGER.info(
        "Mapping %d expressed bins through %d of %d parquet row groups",
        len(linear_indices),
        len(used_groups),
        mapping_file.metadata.num_row_groups,
    )

    cell_totals = np.zeros(len(cell_ids), dtype=np.int64)
    nucleus_totals = np.zeros(len(cell_ids), dtype=np.int64)
    cell_bin_counts = np.zeros(len(cell_ids), dtype=np.int64)
    nucleus_bin_counts = np.zeros(len(cell_ids), dtype=np.int64)
    cell_id_to_index = {int(cell_id): index for index, cell_id in enumerate(cell_ids)}

    mapped_bin_count = 0
    mapped_expression = 0
    nucleus_bin_count = 0
    nucleus_expression = 0
    cell_barcode_ids = np.full(len(counts), -1, dtype=np.int64)
    in_nucleus_flags = np.zeros(len(counts), dtype=bool)

    progress_interval = max(1, len(used_groups) // 20)
    for group_number, group_index in enumerate(used_groups, start=1):
        selected = np.arange(group_starts[group_index], group_ends[group_index])
        row_start = int(group_offsets[group_index])
        local_rows = linear_indices[selected] - row_start
        table = mapping_file.read_row_group(
            int(group_index),
            columns=["cell_id", "in_nucleus"],
        )
        selected_rows = pa.array(local_rows, type=pa.int64())
        mapped_barcodes = table.column("cell_id").take(selected_rows)
        mapped_nuclei = table.column("in_nucleus").take(selected_rows)

        extracted_ids = pc.struct_field(
            pc.extract_regex(mapped_barcodes, CELL_BARCODE_PATTERN), "cell_id"
        )
        parsed_ids = pc.cast(extracted_ids, pa.int64(), safe=True)
        malformed = pc.and_(pc.is_valid(mapped_barcodes), pc.is_null(parsed_ids))
        if pc.any(malformed).as_py():
            raise ValueError("Unexpected segmented cell barcode in barcode_mappings.parquet")

        parsed_ids = pc.fill_null(parsed_ids, -1).to_numpy(zero_copy_only=False)
        covered_by_nucleus = pc.fill_null(mapped_nuclei, False).to_numpy(
            zero_copy_only=False
        )
        valid_ids = parsed_ids >= 0
        if np.any(valid_ids & ~np.isin(parsed_ids, cell_ids)):
            unknown_id = int(parsed_ids[valid_ids & ~np.isin(parsed_ids, cell_ids)][0])
            raise ValueError(
                f"Mapped bin refers to cell_id {unknown_id}, absent from cell GeoJSON"
            )
        cell_barcode_ids[selected] = parsed_ids
        in_nucleus_flags[selected] = covered_by_nucleus

        valid = cell_barcode_ids[selected] >= 0
        if np.any(valid):
            selected_bins = selected[valid]
            cell_indices = np.fromiter(
                (cell_id_to_index[cell_id] for cell_id in cell_barcode_ids[selected_bins]),
                dtype=np.int64,
                count=len(selected_bins),
            )
            bin_counts = counts[selected_bins]
            np.add.at(cell_totals, cell_indices, bin_counts)
            np.add.at(cell_bin_counts, cell_indices, 1)
            mapped_bin_count += len(selected_bins)
            mapped_expression += int(bin_counts.sum())

            covered = in_nucleus_flags[selected_bins]
            if np.any(covered):
                nucleus_indices = cell_indices[covered]
                nucleus_counts_for_bins = bin_counts[covered]
                np.add.at(nucleus_totals, nucleus_indices, nucleus_counts_for_bins)
                np.add.at(nucleus_bin_counts, nucleus_indices, 1)
                nucleus_bin_count += int(np.count_nonzero(covered))
                nucleus_expression += int(nucleus_counts_for_bins.sum())

        if group_number % progress_interval == 0 or group_number == len(used_groups):
            LOGGER.info(
                "Mapped %d / %d row groups; %d / %d expressed bins assigned to cells",
                group_number,
                len(used_groups),
                mapped_bin_count,
                len(linear_indices),
            )

    return (
        cell_totals,
        nucleus_totals,
        cell_bin_counts,
        nucleus_bin_counts,
        mapped_bin_count,
        mapped_expression,
        nucleus_bin_count,
        nucleus_expression,
    )


def write_results(output_prefix, sample_results):
    """Write all cells and per-sample expression summary rows."""
    LOGGER.info("Writing results for all cells across %d samples", len(sample_results))
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    cells_path = output_prefix.with_name(f"{output_prefix.name}.all_cells.csv")
    summary_path = output_prefix.with_name(f"{output_prefix.name}.sample_summary.csv")
    summary_fields = [
        "sample_name",
        "cells_in_geojson",
        "cells_with_expression",
        "all_mapped_gene_expression",
        "all_mapped_nucleus_bin_gene_expression",
        "all_mapped_nucleus_expression_percent",
        "mapped_expressed_bins",
        "nucleus_covered_expressed_bins",
    ]

    with cells_path.open("w", newline="") as cells_file, summary_path.open(
        "w", newline=""
    ) as summary_file:
        cell_writer = csv.DictWriter(cells_file, fieldnames=OUTPUT_COLUMNS)
        summary_writer = csv.DictWriter(summary_file, fieldnames=summary_fields)
        cell_writer.writeheader()
        summary_writer.writeheader()

        for result in sample_results:
            sample_name = result["sample_name"]
            cell_ids = result["cell_ids"]
            cell_areas = result["cell_areas"]
            nucleus_areas = result["nucleus_areas"]
            totals = result["cell_totals"]
            nucleus_totals = result["nucleus_totals"]
            bin_counts = result["cell_bin_counts"]
            nucleus_bin_counts = result["nucleus_bin_counts"]
            ordering = np.lexsort((cell_ids, -totals))
            percentages = np.divide(
                nucleus_totals,
                totals,
                out=np.zeros(len(totals), dtype=np.float64),
                where=totals > 0,
            ) * 100
            result["percentages"] = percentages

            for rank, index in enumerate(ordering, start=1):
                cell_writer.writerow(
                    {
                        "sample_name": sample_name,
                        "rank": rank,
                        "cell_id": int(cell_ids[index]),
                        "cell_area": float(cell_areas[index]),
                        "nucleus_area": float(nucleus_areas[index]),
                        "nucleus_to_cell_area_ratio": float(
                            nucleus_areas[index] / cell_areas[index]
                        ),
                        "total_gene_expression": int(totals[index]),
                        "nucleus_bin_gene_expression": int(nucleus_totals[index]),
                        "nucleus_expression_fraction": percentages[index] / 100,
                        "nucleus_expression_percent": percentages[index],
                        "bins_assigned_to_cell": int(bin_counts[index]),
                        "bins_covered_by_nucleus": int(nucleus_bin_counts[index]),
                    }
                )

            total_expression = int(totals.sum())
            nucleus_expression = int(nucleus_totals.sum())
            summary_writer.writerow(
                {
                    "sample_name": sample_name,
                    "cells_in_geojson": len(cell_ids),
                    "cells_with_expression": int(np.count_nonzero(totals)),
                    "all_mapped_gene_expression": total_expression,
                    "all_mapped_nucleus_bin_gene_expression": nucleus_expression,
                    "all_mapped_nucleus_expression_percent": (
                        100 * nucleus_expression / total_expression if total_expression else 0.0
                    ),
                    "mapped_expressed_bins": result["mapped_bin_count"],
                    "nucleus_covered_expressed_bins": result["nucleus_bin_count"],
                }
            )

    LOGGER.info("Wrote all-cell results to %s", cells_path)
    LOGGER.info("Wrote per-sample summaries to %s", summary_path)
    return cells_path, summary_path


def plot_sample_violins(output_prefix, sample_results):
    """Plot each sample's per-cell nucleus-expression percentage distribution."""
    plot_path = output_prefix.with_name(
        f"{output_prefix.name}.nucleus_expression_percent_violin.png"
    )
    values = [result["percentages"] for result in sample_results]
    labels = [result["sample_name"] for result in sample_results]
    figure_width = max(8, 1.2 * len(labels) + 4)
    figure, axis = plt.subplots(figsize=(figure_width, 6))
    violin = axis.violinplot(values, showmeans=False, showmedians=True, showextrema=True)
    for body in violin["bodies"]:
        body.set_facecolor("#397a72")
        body.set_edgecolor("#24564f")
        body.set_alpha(0.75)
    axis.set_xticks(np.arange(1, len(labels) + 1), labels, rotation=30, ha="right")
    axis.set_ylabel("Nucleus-covered expression (% of cell's assigned bin expression)")
    axis.set_title("Per-cell nucleus-covered expression by sample")
    axis.set_ylim(0, 100)
    axis.grid(axis="y", alpha=0.25)
    figure.tight_layout()
    figure.savefig(str(plot_path), dpi=200)
    plt.close(figure)
    LOGGER.info("Wrote cross-sample percentage violin plot to %s", plot_path)
    return plot_path


def plot_sample_area_ratio_violins(output_prefix, sample_results):
    """Plot each sample's per-cell nucleus-to-cell area-ratio distribution."""
    plot_path = output_prefix.with_name(
        f"{output_prefix.name}.nucleus_to_cell_area_ratio_violin.png"
    )
    values = [
        result["nucleus_areas"] / result["cell_areas"]
        for result in sample_results
    ]
    labels = [result["sample_name"] for result in sample_results]
    figure_width = max(8, 1.2 * len(labels) + 4)
    figure, axis = plt.subplots(figsize=(figure_width, 6))
    violin = axis.violinplot(values, showmeans=False, showmedians=True, showextrema=True)
    for body in violin["bodies"]:
        body.set_facecolor("#b16a45")
        body.set_edgecolor("#78452e")
        body.set_alpha(0.75)
    axis.set_xticks(np.arange(1, len(labels) + 1), labels, rotation=30, ha="right")
    axis.set_ylabel("Nucleus area / cell area (GeoJSON coordinate units)")
    axis.set_title("Per-cell nucleus-to-cell area ratio by sample")
    axis.grid(axis="y", alpha=0.25)
    figure.tight_layout()
    figure.savefig(str(plot_path), dpi=200)
    plt.close(figure)
    LOGGER.info("Wrote cross-sample area-ratio violin plot to %s", plot_path)
    return plot_path


def plot_sample_ratio_hotspots(output_prefix, sample_results):
    """Plot log-scaled cell-density hotspots by sample."""
    plot_path = output_prefix.with_name(
        f"{output_prefix.name}.area_vs_expression_ratio_hotspots.png"
    )
    sample_count = len(sample_results)
    column_count = min(2, sample_count)
    row_count = math.ceil(sample_count / column_count)
    figure, axes = plt.subplots(
        row_count,
        column_count,
        figsize=(6.5 * max(row_count, column_count),) * 2,
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    hexbin = None
    for axis, result in zip(axes.flat, sample_results):
        cell_areas = result["cell_areas"]
        area_ratios = result["nucleus_areas"] / cell_areas
        expression_percentages = np.divide(
            result["nucleus_totals"],
            result["cell_totals"],
            out=np.zeros(len(cell_areas), dtype=np.float64),
            where=result["cell_totals"] > 0,
        ) * 100
        hexbin = axis.hexbin(
            area_ratios,
            expression_percentages,
            gridsize=55,
            mincnt=1,
            bins="log",
            cmap="magma",
            linewidths=0,
        )
        axis.set_title(result["sample_name"])
        axis.set_box_aspect(1.0)
        axis.set_ylim(0, 100)
        axis.grid(alpha=0.12)

    for axis in axes[-1, :]:
        axis.set_xlabel("Nucleus area / cell area")
    for axis in axes[:, 0]:
        axis.set_ylabel("Nucleus-covered gene expression (%)")
    for axis in axes.flat[sample_count:]:
        axis.set_visible(False)
    figure.suptitle("Cell density: area ratio versus gene-expression ratio by sample")
    figure.subplots_adjust(left=0.08, right=0.88, bottom=0.09, top=0.91, wspace=0.12, hspace=0.2)
    color_axis = figure.add_axes([0.91, 0.17, 0.018, 0.68])
    figure.colorbar(hexbin, cax=color_axis, label="Cells per hexagon (log scale)")
    figure.savefig(str(plot_path), dpi=200)
    plt.close(figure)
    LOGGER.info("Wrote per-sample ratio hotspot figure to %s", plot_path)
    return plot_path


def read_file_map(file_map):
    """Read headerless nucleus-GeoJSON path and sample-name rows."""
    entries = []
    with file_map.open() as map_file:
        for line_number, line in enumerate(map_file, start=1):
            line = line.strip()
            if not line:
                continue
            fields = line.rsplit(maxsplit=1)
            if len(fields) != 2:
                raise ValueError(f"Invalid file-map row {line_number}: expected path and sample name")
            nucleus_path = Path(fields[0])
            sample_name = fields[1]
            sample_dir = nucleus_path.parent.parent
            cell_path = nucleus_path.with_name("cell_segmentations.geojson")
            matrix_path = sample_dir / "binned_outputs" / BIN_RESOLUTION / "filtered_feature_bc_matrix.h5"
            mapping_path = sample_dir / "barcode_mappings.parquet"
            entries.append((sample_name, sample_dir, cell_path, nucleus_path, matrix_path, mapping_path))
    if not entries:
        raise ValueError(f"No sample rows found in file map: {file_map}")
    names = [entry[0] for entry in entries]
    if len(names) != len(set(names)):
        raise ValueError("Sample names in file map must be unique")
    return entries


def process_sample(sample_name, sample_dir, cell_geojson, nucleus_geojson, matrix_path, mapping_path):
    """Read and aggregate one sample's expressed bins for all segmented cells."""
    for required_path in (cell_geojson, nucleus_geojson, matrix_path, mapping_path):
        if not required_path.is_file():
            raise FileNotFoundError(f"Required input for {sample_name} does not exist: {required_path}")
    LOGGER.info("[%s] Validated inputs in %s", sample_name, sample_dir)
    cell_ids, cell_areas = read_cell_ids_and_areas(cell_geojson)
    nucleus_ids, nucleus_areas = read_cell_ids_and_areas(nucleus_geojson)
    if not np.array_equal(cell_ids, nucleus_ids):
        raise ValueError(f"Cell and nucleus GeoJSONs have different cell_id sets for {sample_name}")
    LOGGER.info("[%s] Cell and nucleus GeoJSONs contain matching IDs", sample_name)

    LOGGER.info("[%s] Reading sparse 2 µm binned gene-expression counts", sample_name)
    coordinates, counts = read_bin_counts(matrix_path)
    LOGGER.info("[%s] Found %d expressed bins", sample_name, len(counts))
    LOGGER.info("[%s] Mapping expressed bins to cells and nuclei", sample_name)
    (
        cell_totals,
        nucleus_totals,
        cell_bin_counts,
        nucleus_bin_counts,
        mapped_bin_count,
        mapped_expression,
        nucleus_bin_count,
        nucleus_expression,
    ) = accumulate_counts_by_cell(coordinates, counts, mapping_path, cell_ids)
    LOGGER.info(
        "[%s] Mapped %d bins (%d UMIs); %d bins (%d UMIs) nucleus-covered",
        sample_name,
        mapped_bin_count,
        mapped_expression,
        nucleus_bin_count,
        nucleus_expression,
    )
    return {
        "sample_name": sample_name,
        "cell_ids": cell_ids,
        "cell_areas": cell_areas,
        "nucleus_areas": nucleus_areas,
        "cell_totals": cell_totals,
        "nucleus_totals": nucleus_totals,
        "cell_bin_counts": cell_bin_counts,
        "nucleus_bin_counts": nucleus_bin_counts,
        "mapped_bin_count": mapped_bin_count,
        "nucleus_bin_count": nucleus_bin_count,
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate 2 µm bin expression for all segmented cells and compare "
            "nucleus-covered expression percentages across samples."
        )
    )
    parser.add_argument(
        "--file-map",
        type=Path,
        required=True,
        help="Headerless rows of nucleus_segmentations.geojson path and sample name",
    )
    parser.add_argument("--output-prefix", type=Path, required=True, help="Prefix for output files")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    file_map_entries = read_file_map(args.file_map)
    sample_results = [
        process_sample(*entry)
        for entry in file_map_entries
    ]
    cells_path, summary_path = write_results(args.output_prefix, sample_results)
    violin_path = plot_sample_violins(args.output_prefix, sample_results)
    area_ratio_violin_path = plot_sample_area_ratio_violins(args.output_prefix, sample_results)
    hotspot_path = plot_sample_ratio_hotspots(args.output_prefix, sample_results)
    LOGGER.info("All-cell results: %s", cells_path)
    LOGGER.info("Per-sample summary: %s", summary_path)
    LOGGER.info("Cross-sample violin plot: %s", violin_path)
    LOGGER.info("Cross-sample area-ratio violin plot: %s", area_ratio_violin_path)
    LOGGER.info("Per-sample area/expression hotspot figure: %s", hotspot_path)


if __name__ == "__main__":
    main()