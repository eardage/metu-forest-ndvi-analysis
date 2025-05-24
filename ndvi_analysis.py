import rasterio as rs
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
from rasterio.mask import mask 
from shapely.geometry import MultiPolygon 


file_path_may = "/Users/egeardaozturk/forest_detection/2025-05-15-00_00_2025-05-15-23_59_Sentinel-2_L2A_NDVI.tiff"

file_path_april = "/Users/egeardaozturk/forest_detection/2025-04-15-00:00_2025-04-15-23:59_Sentinel-2_L2A_NDVI.tiff"

try:

    with rs.open(file_path_may) as src_may:

        ndvi_may_array = src_may.read(1)

        # Array'in boyutunu kontrol et
        print(f"May NDVI array shape: {ndvi_may_array.shape}")
        print(f"May NDVI array data type: {ndvi_may_array.dtype}")


        print("May NDVI array's first 5x5 pixels:")
        print(ndvi_may_array[:5, :5])
    
    with rs.open(file_path_april) as src_april:
        ndvi_april_array = src_april.read(1)
        print(f"April NDVI array shape: {ndvi_april_array.shape}")
        print(f"April NDVI array data type: {ndvi_april_array.dtype}")


        print("April NDVI array's first 5x5 pixels::")
        print(ndvi_april_array[:5, :5])

        if ndvi_may_array.shape != ndvi_april_array.shape:
            print("Shapes are not equal, which is not possible to calculate the difference.")
        else:
            print("Shapes are equal, we can calculate the ndvi difference!")
except Exception as e:
    print(f"Error reading the file: {e}")

ndvi_difference = ndvi_may_array - ndvi_april_array
print(f"NDVI Difference Data Type: {ndvi_difference.dtype}")
print("NDVI Difference first 5x5 pixels")
print(ndvi_difference[:5, :5])

print(f"\nNDVI Difference min value: {np.nanmin(ndvi_difference)}")
print(f"NDVI Difference max value: {np.nanmax(ndvi_difference)}")

print(f"NDVI Difference mean value: {np.nanmean(ndvi_difference)}")
print(f"NDVI Difference median value: {np.nanmedian(ndvi_difference)}")
print(f"NDVI Difference standard deviation value: {np.nanstd(ndvi_difference)}")

plt.imshow(ndvi_difference, cmap='RdYlGn', vmin=-1, vmax=1)
plt.colorbar(label='NDVI Difference')
plt.title('NDVI Difference between May and April')
plt.show()



# GeoJSON dosyasının adı ve yolu
forest_boundary_geojson_path = "/Users/egeardaozturk/forest_detection/ODTÜ_Ormanı_Polgon.geojson"


try:
    # GeoJSON dosyasını GeoDataFrame olarak oku
    forest_gdf = gpd.read_file(forest_boundary_geojson_path)


    geometries_to_mask = [geom for geom in forest_gdf.geometry if geom.geom_type in ['Polygon', 'MultiPolygon']]

    if not geometries_to_mask:
        raise ValueError("GeoJSON dosyasında maskeleme için uygun (Polygon veya MultiPolygon) geometri bulunamadı.")


    with rs.open(file_path_may) as src_may:
        masked_may_ndvi, _ = mask(src_may, geometries_to_mask, crop=True, nodata=np.nan)
        masked_may_ndvi = masked_may_ndvi[0] if masked_may_ndvi.ndim == 3 else masked_may_ndvi

    with rs.open(file_path_april) as src_april:
        masked_april_ndvi, _ = mask(src_april, geometries_to_mask, crop=True, nodata=np.nan)
        masked_april_ndvi = masked_april_ndvi[0] if masked_april_ndvi.ndim == 3 else masked_april_ndvi

    # SONRA farkı alması gerekirdi
    # BURADA BOYUT UYUMLULUK SORUNLARI OLABİLİR EĞER CROP=TRUE İKİSİNİ FARKLI KIRPARSA
        masked_ndvi_difference = masked_may_ndvi - masked_april_ndvi
        


        print(f"Maskelenmiş NDVI Fark array'inin boyutu: {masked_ndvi_difference.shape}")
        print(f"Maskelenmiş NDVI Fark array'inin veri tipi: {masked_ndvi_difference.dtype}")
        print(f"Maskelenmiş NDVI Fark array'inde NaN (NoData) değeri var mı? {np.isnan(masked_ndvi_difference).any()}")

        # Statistics of masked NDVI difference (Without NaN values)
        print(f"\nMaskelenmiş NDVI Fark min value: {np.nanmin(masked_ndvi_difference)}")
        print(f"Maskelenmiş NDVI Fark max value: {np.nanmax(masked_ndvi_difference)}")
        print(f"Maskelenmiş NDVI Fark mean value: {np.nanmean(masked_ndvi_difference)}")
        print(f"Maskelenmiş NDVI Fark median value: {np.nanmedian(masked_ndvi_difference)}")
        print(f"Maskelenmiş NDVI Fark standard deviation value: {np.nanstd(masked_ndvi_difference)}")

        # Visualize the masked NDVI difference
        plt.figure(figsize=(10, 8))

        plt.imshow(masked_ndvi_difference, cmap='RdYlGn', vmin=-0.5, vmax=0.5) # veya -0.5, 0.5 dene
        plt.colorbar(label='NDVI Fark (Mayıs - Nisan)')
        plt.title('ODTÜ Ormanı NDVI Fark Haritası (Maskelenmiş)')
        plt.xlabel('X (Piksel)')
        plt.ylabel('Y (Piksel)')
        plt.show()

except Exception as e:
    print(f"Maskeleme işlemi sırasında bir hata oluştu: {e}")
    print("GeoJSON dosya yolunu ve içeriğini kontrol et.")
    print("Orijinal NDVI dosyalarının doğru bir şekilde açıldığından emin ol.")
 # İndirdiğin GeoJSON dosyasının tam yolunu buraya yaz



 #Histogram

valid_ndvi_diff_values = masked_ndvi_difference[~np.isnan(masked_ndvi_difference)]


plt.figure(figsize=(10, 6))
plt.hist(valid_ndvi_diff_values, bins=50, range=(-1, 1), color='skyblue', edgecolor='black')
plt.title('Histogram of NDVI Difference')
plt.xlabel('NDVI Difference (May - April)')
plt.ylabel('Number of Pixels')
plt.grid(axis='y', alpha=0.75)
plt.axvline(x=0, color='red', linestyle='dashed', linewidth=1)
plt.text(0.05, plt.ylim()[1]*0.9, 'No Difference', color='red', ha='center', va='center', transform=plt.gca().transAxes)
plt.show()


print("Max value of Masked NDVI Difference: ", np.max(masked_ndvi_difference))
print("Min value of Masked NDVI Difference: ", np.min(masked_ndvi_difference))
print("Mean value of Masked NDVI Difference: ", np.mean(masked_ndvi_difference))
print("Median value of Masked NDVI Difference: ", np.median(masked_ndvi_difference))
print("Standard deviation of Masked NDVI Difference: ", np.std(masked_ndvi_difference))

threshold_growth = 0.05
threshold_loss = -0.05

classified_map = np.full(masked_ndvi_difference.shape, np.nan, dtype=float)

#Significant growth
classified_map[masked_ndvi_difference > threshold_growth] = 1

#Significant loss
classified_map[masked_ndvi_difference < threshold_loss] = -1

#No significant change  
classified_map[np.logical_and(masked_ndvi_difference >= threshold_loss, masked_ndvi_difference <= threshold_growth)] = 0

print("Classified ChangeMap:")
print(f"Loss (-1) threshold: {threshold_loss}")
print(f"Growth (+1) threshold: {threshold_growth}")
print(f"No change (0) threshold: {threshold_loss} to {threshold_growth}")

#Area calculation

pixel_area_sqm = 10*10 # 10 meter x 10 meter = 10 m^2 in Sentinel-2 resolution

#Total pixel count

total_pixels_in_forest = np.sum(~np.isnan(masked_ndvi_difference))
loss_pixels = np.sum(classified_map == -1)
growth_pixels = np.sum(classified_map == 1)
no_change_pixels = total_pixels_in_forest - (loss_pixels + growth_pixels)

#Areas in square meters
loss_area_sqm = loss_pixels * pixel_area_sqm
growth_area_sqm = growth_pixels * pixel_area_sqm
no_change_area_sqm = no_change_pixels * pixel_area_sqm

#Total area
total_area_sqm = loss_area_sqm + growth_area_sqm + no_change_area_sqm

#Convert to hectares (1 ha = 10,000 m^2 1 km^2 = 1,000,000 m^2)
loss_area_ha = loss_area_sqm / 10000
growth_area_ha = growth_area_sqm / 10000
no_change_area_ha = no_change_area_sqm / 10000
total_area_ha = total_area_sqm / 10000

#Convert to square kilometers
loss_area_km2 = loss_area_sqm / 1000000
growth_area_km2 = growth_area_sqm / 1000000
no_change_area_km2 = no_change_area_sqm / 1000000
total_area_km2 = total_area_sqm / 1000000

#Print results 
print(f"Total area of the forest: {total_area_ha:.2f} hectares, pixels: {total_pixels_in_forest}")
print(f"Area of loss: {loss_area_ha:.2f} hectares, pixels: {loss_pixels}")
print(f"Area of growth: {growth_area_ha:.2f} hectares, pixels: {growth_pixels}")
print(f"Area of no change: {no_change_area_ha:.2f} hectares, pixels: {no_change_pixels}")

#Create a color map for the classified map
plt.figure(figsize=(10, 8))
plt.imshow(classified_map, cmap='RdYlGn', vmin=-1, vmax=1)
plt.colorbar(label='NDVI Difference Classification (-1: Loss, 0: No Change, 1: Growth)', ticks=[-1, 0, 1])
plt.title(f'Classified Change Map (May - April) - Threshold: {threshold_loss} to {threshold_growth}')
plt.xlabel('X (Pixels)')
plt.ylabel('Y (Pixels)')
plt.show()





















