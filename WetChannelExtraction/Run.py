
import Pre_Processing
import Wet_Channel_Extraction
import os
## Inputs

# 1. the path to the bare-earth DEM
ground_DEM = 'C:\\Data\\Codes\\Comb_Code\\Input\\gro_DEM.tif'

# 2. the path to the vegetation  DEM (created from all LiDAR points)
veg_DEM = 'C:\\Data\\Codes\\Comb_Code\\Input\\veg_DEM.tif'

# 3. the path to the intensity  raster
int_raster = 'C:\\Data\\Codes\\Comb_Code\\Input\\gro_INT.tif'

# 4. valley network data
flowdir_raster = "C:\\Data\\Codes\\Comb_Code\\Output\\driange_network\\maps\\Modified_Fdir_8.tif"
valley_raster = "C:\\Data\\Codes\\Comb_Code\\Output\\driange_network\\maps\\valley_network.tif"

# 5. output folder
output_path = 'C:\\Data\\Codes\\Comb_Code\\Output'

# 6. resolution data
unit = 'm'
pix_per_m = 0.5
ele_per_meter = 1

# 7. paramters 
area_thresh = 1 #m^2
diff_thresh = 2 # m for vegitation
min_slope = 0 # m/m
max_slope = 2 # m/m
min_pos_slope = 0.75 # m/m
max_gap = 15 # m


os.mkdir(output_path + '\\fixed_input')
Pre_Processing.Input_Make(flowdir_raster, valley_raster, \
                                  ground_DEM , veg_DEM , int_raster, output_path + '\\fixed_input')
os.mkdir(output_path + '\\wet_output')
Wet_Channel_Extraction.Wet_Channel_Extraction(output_path + '\\fixed_input' , unit ,area_thresh, diff_thresh, min_slope, max_slope, \
                       min_pos_slope, max_gap,pix_per_m , ele_per_meter , output_path + '\\wet_output')
