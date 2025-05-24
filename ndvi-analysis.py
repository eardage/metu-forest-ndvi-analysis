import rasterio as rs
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
from rasterio.mask import mask
from shapely.geometry import MultiPolygon # Keep if MultiPolygon is explicitly used/needed

# --- 1. Data Loading and Basic Operations Functions ---

def load_ndvi_data(file_path):

    with rs.open(file_path) as src:
        ndvi_array = src.read(1).astype(float) # Convert to float type is important
        profile = src.profile # Get meta information
        print(f"Loaded {file_path}: shape={ndvi_array.shape}, min={np.nanmin(ndvi_array):.4f}, max={np.nanmax(ndvi_array):.4f}")
        return ndvi_array, src.transform, src.crs, profile

def calculate_ndvi_difference(ndvi_may, ndvi_april):

    if ndvi_may.shape != ndvi_april.shape:
        raise ValueError("NDVI array shapes are not equal. Cannot calculate difference.")
    return ndvi_may - ndvi_april

def print_ndvi_difference_stats(ndvi_diff_array, title="NDVI Difference"):

    # Get statistics excluding NaN values
    valid_values = ndvi_diff_array[~np.isnan(ndvi_diff_array)]
    print(f"\n--- {title} Statistics ---")
    print(f"Min value (excluding NaN): {np.nanmin(valid_values):.4f}")
    print(f"Max value (excluding NaN): {np.nanmax(valid_values):.4f}")
    print(f"Mean value (excluding NaN): {np.nanmean(valid_values):.4f}")
    print(f"Median value (excluding NaN): {np.nanmedian(valid_values):.4f}")
    print(f"Standard deviation (excluding NaN): {np.nanstd(valid_values):.4f}")

def visualize_ndvi_difference(ndvi_diff_array, title="NDVI Difference Map", vmin=-1, vmax=1):

    plt.figure(figsize=(10, 8))
    plt.imshow(ndvi_diff_array, cmap='RdYlGn', vmin=vmin, vmax=vmax)
    plt.colorbar(label='NDVI Difference')
    plt.title(title)
    plt.xlabel('X (Pixel)')
    plt.ylabel('Y (Pixel)')
    plt.show()

# --- 2. Masking Functions ---

def load_and_validate_geojson(geojson_path):

    forest_gdf = gpd.read_file(geojson_path)
    geometries = [geom for geom in forest_gdf.geometry if geom.geom_type in ['Polygon', 'MultiPolygon']]
    if not geometries:
        raise ValueError(f"No suitable (Polygon or MultiPolygon) geometry found in GeoJSON file: {geojson_path}")
    return geometries

def mask_raster_with_geojson(ndvi_array, transform, crs, profile, geometries_to_mask):

    # Update profile: single band (count=1), float type, correct CRS
    updated_profile = profile.copy()
    updated_profile.update(dtype=ndvi_array.dtype, count=1, crs=crs)

    with rs.MemoryFile() as memfile:
        with memfile.open(**updated_profile) as temp_dataset:
            temp_dataset.write(ndvi_array, 1) # Write NDVI difference to temporary dataset

            # Perform the masking operation
            masked_data, out_transform = mask(temp_dataset, geometries_to_mask, crop=True, nodata=np.nan)

            if masked_data.ndim == 3: # After mask, it might be (band, height, width)
                return masked_data[0] # Take the first band
            return masked_data

def plot_histogram(data_array, title="Histogram of NDVI Difference", xlabel="NDVI Difference", bins=50, range_vals=(-1, 1)):

    valid_values = data_array[~np.isnan(data_array)]
    plt.figure(figsize=(10, 6))
    plt.hist(valid_values, bins=bins, range=range_vals, color='skyblue', edgecolor='black')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('Number of Pixels')
    plt.grid(axis='y', alpha=0.75)
    plt.axvline(x=0, color='red', linestyle='dashed', linewidth=1)
    plt.show()

# --- 3. Classification and Area Calculation Functions ---

def classify_ndvi_difference(ndvi_diff_array, threshold_loss, threshold_growth):

    classified_map = np.full(ndvi_diff_array.shape, np.nan, dtype=float)
    classified_map[ndvi_diff_array > threshold_growth] = 1 # Growth
    classified_map[ndvi_diff_array < threshold_loss] = -1 # Loss
    # For no significant change: non-NaN values that are not in other classes
    no_change_mask = (ndvi_diff_array >= threshold_loss) & \
                     (ndvi_diff_array <= threshold_growth) & \
                     (~np.isnan(ndvi_diff_array))
    classified_map[no_change_mask] = 0 # No significant change

    print(f"\nChange Classification:")
    print(f"Loss (-1) Threshold: < {threshold_loss}")
    print(f"Growth (+1) Threshold: > {threshold_growth}")
    print(f"No Change (0) Threshold: {threshold_loss} to {threshold_growth}")
    return classified_map

def calculate_and_print_areas(classified_map, pixel_area_sqm=100):

    total_pixels_in_forest = np.sum(~np.isnan(classified_map)) # Count of non-NaN pixels
    loss_pixels = np.sum(classified_map == -1)
    growth_pixels = np.sum(classified_map == 1)
    no_change_pixels = np.sum(classified_map == 0) # Include 0 class

    # Areas in square meters
    area_loss_sqm = loss_pixels * pixel_area_sqm
    area_growth_sqm = growth_pixels * pixel_area_sqm
    area_no_change_sqm = no_change_pixels * pixel_area_sqm
    total_forest_area_sqm = total_pixels_in_forest * pixel_area_sqm

    # Convert to hectares (1 ha = 10,000 m^2)
    area_loss_ha = area_loss_sqm / 10000
    area_growth_ha = area_growth_sqm / 10000
    area_no_change_ha = area_no_change_sqm / 10000
    total_forest_area_ha = total_forest_area_sqm / 10000

    print(f"\n--- Area Calculations ---")
    print(f"Total ODTU Forest Area (Masked): {total_forest_area_ha:.2f} hectares ({total_pixels_in_forest} pixels)")
    print(f"Area of Loss: {area_loss_ha:.2f} hectares ({loss_pixels} pixels)")
    print(f"Area of Growth: {area_growth_ha:.2f} hectares ({growth_pixels} pixels)")
    print(f"Area of No Change: {area_no_change_ha:.2f} hectares ({no_change_pixels} pixels)")

    return {
        "total_ha": total_forest_area_ha,
        "loss_ha": area_loss_ha,
        "growth_ha": area_growth_ha,
        "no_change_ha": area_no_change_ha
    }

def visualize_classified_map(classified_map, title="Classified Change Map",
                             threshold_loss=-0.05, threshold_growth=0.05):

    plt.figure(figsize=(10, 8))
    # RdYlGn (Red-Yellow-Green) palette is suitable for -1, 0, 1 values.
    # Red (-1) -> Yellow (0) -> Green (1)
    cmap = plt.cm.RdYlGn
    # Set vmin and vmax according to classification values
    norm = plt.Normalize(vmin=-1, vmax=1)
    
    plt.imshow(classified_map, cmap=cmap, norm=norm)
    cbar = plt.colorbar(ticks=[-1, 0, 1], label='NDVI Change Class')
    cbar.ax.set_yticklabels(['Loss (-1)', 'No Change (0)', 'Growth (1)'])

    plt.title(f'{title} (Threshold: {threshold_loss} / {threshold_growth})')
    plt.xlabel('X (Pixels)')
    plt.ylabel('Y (Pixels)')
    plt.show()

# --- Main Workflow ---

def main():

    # File paths
    file_path_may = "/Users/egeardaozturk/METU_Forest_NDVI_Analysis/data/2025-05-15-00_00_2025-05-15-23_59_Sentinel-2_L2A_NDVI.tiff"
    file_path_april = "/Users/egeardaozturk/METU_Forest_NDVI_Analysis/data/2025-04-15-00:00_2025-04-15-23:59_Sentinel-2_L2A_NDVI.tiff"
    forest_boundary_geojson_path = "/Users/egeardaozturk/METU_Forest_NDVI_Analysis/METU_Forest_boundary.geojson"

    # --- 1. Load NDVI Data and Calculate Difference ---
    try:
        print("Loading NDVI data...")
        ndvi_may_array, may_transform, may_crs, may_profile = load_ndvi_data(file_path_may)
        ndvi_april_array, april_transform, april_crs, april_profile = load_ndvi_data(file_path_april)

        print("Calculating NDVI difference...")
        ndvi_difference = calculate_ndvi_difference(ndvi_may_array, ndvi_april_array)
        print_ndvi_difference_stats(ndvi_difference, "Overall NDVI Difference (Before Masking)")
        visualize_ndvi_difference(ndvi_difference, "NDVI Difference (May - April) - Before Masking")

    except Exception as e:
        print(f"Error during NDVI loading or difference calculation: {e}")
        return # Stop program if error occurs

    # --- 2. Load Forest Boundary and Mask NDVI Difference ---
    try:
        print("\nLoading forest boundary GeoJSON...")
        geometries_to_mask = load_and_validate_geojson(forest_boundary_geojson_path)

        print("Masking NDVI difference map according to forest boundary...")
        # We can use may_profile for masking, as their dimensions are the same.
        masked_ndvi_difference = mask_raster_with_geojson(
            ndvi_difference, may_transform, may_crs, may_profile, geometries_to_mask
        )
        print_ndvi_difference_stats(masked_ndvi_difference, "Masked NDVI Difference (ODTU Forest)")
        visualize_ndvi_difference(masked_ndvi_difference, "ODTU Forest NDVI Difference Map (Masked)", vmin=-0.4, vmax=0.4) # Narrower range

    except Exception as e:
        print(f"Error during masking operation: {e}")
        return # Stop program if error occurs

    # --- 3. Plot Histogram ---
    print("\nPlotting histogram of masked NDVI difference...")
    plot_histogram(masked_ndvi_difference, title='ODTU Forest Masked NDVI Difference Histogram (May - April)',
                   xlabel='NDVI Difference Value (May - April)', bins=50, range_vals=(-1, 1))

    # --- 4. Classify Changes and Calculate Area ---
    # These thresholds are determined based on the 2025 April-May histogram.
    threshold_growth = 0.05
    threshold_loss = -0.05

    print("\nClassifying NDVI difference...")
    classified_map = classify_ndvi_difference(masked_ndvi_difference, threshold_loss, threshold_growth)
    
    print("Calculating areas...")
    calculated_areas = calculate_and_print_areas(classified_map, pixel_area_sqm=10*10)

    # --- 5. Visualize Classified Map ---
    print("\nVisualizing classified map...")
    visualize_classified_map(classified_map, title='ODTU Forest Classified Change Map (May - April)',
                             threshold_loss=threshold_loss, threshold_growth=threshold_growth)

if __name__ == "__main__":
    main()