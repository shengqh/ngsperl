import os

# Environment Setup for Docker Permissions - MUST BE BEFORE TORCH IMPORTS
script_dir = os.getcwd()
torch_dir = os.path.join(script_dir, '.torch')
os.environ['TORCH_HOME'] = torch_dir
os.environ['TORCHINDUCTOR_CACHE_DIR'] = torch_dir

import json
import torch
import argparse
import torchvision.models as models
import torchvision.transforms as transforms
from torch.utils.data import Dataset, DataLoader
import geopandas as gpd
import pandas as pd
import numpy as np
from PIL import Image
import tifffile
from tqdm import tqdm

class VisiumHDCellDataset(Dataset):
    def __init__(self, gdf, full_image, transform=None, padding=2):
        self.gdf = gdf
        self.image = full_image
        self.transform = transform
        self.padding = padding
        # Pre-calculate bounds and IDs for speed
        self.bounds = gdf.geometry.bounds.values
        self.cell_ids = gdf['cell_id'].values

    def __len__(self):
        return len(self.gdf)

    def __getitem__(self, idx):
        minx, miny, maxx, maxy = self.bounds[idx]
        left, top = int(minx) - self.padding, int(miny) - self.padding
        right, bottom = int(maxx) + self.padding, int(maxy) + self.padding

        # Crop from the numpy array (y, x)
        crop = self.image[max(0, top):bottom, max(0, left):right]
        crop_img = Image.fromarray(crop).convert('RGB')

        if self.transform:
            crop_img = self.transform(crop_img)

        return crop_img, self.cell_ids[idx]

def calculate_morphology(gdf, prefix, um_per_px):
    """Calculates area, perimeter, and circularity for a GeoDataFrame."""
    if gdf.crs is None or gdf.crs.is_geographic:
        gdf.crs = "EPSG:3857"

    area_px = gdf.geometry.area
    perimeter_px = gdf.geometry.length

    results = pd.DataFrame({
        'cell_id': gdf['cell_id'].astype(str), # Ensure ID is string for merging
        f'{prefix}_area_um2': area_px * (um_per_px**2),
        f'{prefix}_perimeter_um': perimeter_px * um_per_px,
        f'{prefix}_circularity': (4 * np.pi * area_px) / (perimeter_px**2)
    })
    return results

def main():
    parser = argparse.ArgumentParser(description="Visium HD Cell & Nucleus Feature Extractor")

    # Paths
    parser.add_argument("--cell_geojson", type=str, required=True, help="Path to cell_segmentations.geojson")
    parser.add_argument("--nucleus_geojson", type=str, required=True, help="Path to nucleus_segmentations.geojson")
    parser.add_argument("--image", type=str, required=True, help="Path to tissue TIFF image")
    parser.add_argument("--scales", type=str, required=True, help="Path to scalefactors_json.json")
    parser.add_argument("--output", type=str, default="combined_features.parquet", help="Output path")

    # Settings
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--num_workers", type=int, default=4)
    parser.add_argument("--padding", type=int, default=2)

    args = parser.parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # --- 1. Load SCALES ---
    with open(args.scales) as f:
        scales = json.load(f)
    um_per_px = scales['microns_per_pixel']

    # --- 2. MORPHOLOGY CALCULATION ---
    print("Calculating morphology for cells and nuclei...")
    cell_gdf = gpd.read_file(args.cell_geojson)
    nuc_gdf = gpd.read_file(args.nucleus_geojson)

    cell_morph = calculate_morphology(cell_gdf, "cell", um_per_px)
    nuc_morph = calculate_morphology(nuc_gdf, "nucleus", um_per_px)

    # Merge morphology dataframes on cell_id (both now strings)
    combined_morph = pd.merge(cell_morph, nuc_morph, on='cell_id', how='left')

    # --- 3. IMAGE FEATURE EXTRACTION ---
    print(f"Loading Image: {args.image}")
    full_image = tifffile.imread(args.image)

    # Model definition - using V2 weights for better accuracy
    weights = models.ResNet50_Weights.IMAGENET1K_V2
    preprocess = weights.transforms()
    model = models.resnet50(weights=weights)
    model = torch.nn.Sequential(*(list(model.children())[:-1]))
    model = model.to(device).eval()

    dataset = VisiumHDCellDataset(cell_gdf, full_image, transform=preprocess, padding=args.padding)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False,
                        num_workers=args.num_workers, pin_memory=True)

    all_features, all_ids = [], []
    print(f"Extracting ResNet features for {len(cell_gdf)} cells...")
    with torch.no_grad():
        for batch_imgs, batch_ids in tqdm(loader):
            batch_imgs = batch_imgs.to(device)
            features = model(batch_imgs).view(batch_imgs.size(0), -1)
            all_features.append(features.cpu().numpy())

            # Convert IDs to strings immediately to avoid merge conflicts
            if torch.is_tensor(batch_ids):
                all_ids.extend([str(i) for i in batch_ids.tolist()])
            else:
                all_ids.extend([str(i.item()) if torch.is_tensor(i) else str(i) for i in batch_ids])

    # --- 4. FINAL INTEGRATION ---
    print("Merging all features...")
    df_resnet = pd.DataFrame(np.vstack(all_features).astype(np.float32))
    df_resnet.insert(0, 'cell_id', all_ids)

    # Merge ResNet features with combined morphology (both keys are now strings)
    final_df = pd.merge(df_resnet, combined_morph, on='cell_id', how='left')

    # Reorder columns: Metadata first, then deep features
    meta_cols = ['cell_id', 'cell_area_um2', 'cell_perimeter_um', 'cell_circularity',
                 'nucleus_area_um2', 'nucleus_perimeter_um', 'nucleus_circularity']
    resnet_cols = [c for c in final_df.columns if c not in meta_cols]
    final_df = final_df[meta_cols + resnet_cols]

    print(f"Saving results to: {args.output}")
    final_df.to_parquet(args.output, index=False, engine='pyarrow')

    print(f"Saving first ten cells to: {args.output}.csv for quick inspection")
    final_df.head(10).to_csv(f"{args.output}.csv", index=False)

    print("Process Complete.")

if __name__ == "__main__":
    main()


