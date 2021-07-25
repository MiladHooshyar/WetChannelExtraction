import Pre_Processing
import Wet_Channel_Extraction
import os

# 1. the path to the bare-earth DEM
ground_DEM = '<enter the path to ground elevation dem>'

# 2. the path to the vegetation  DEM (created from all LiDAR points)
veg_DEM = '<enter the path to top of vegetation elevation dem>'

# 3. the path to the intensity  raster
int_raster = '<enter the path to top of LiDAR intensity raster>'

# 4. valley network data
flowdir_raster = '<enter the path to flow direction raster>'
valley_raster = '<enter the path to valley network skeleton raster>'

# 5. output folder
output_path = '<enter the output folder>'

# 6. resolution data
unit = 'm'
pix_per_m = 0.5
ele_per_meter = 1

# 7. parameters
area_thresh = 1  # m^2
diff_thresh = 2  # m for vegetation
min_slope = 0  # m/m
max_slope = 2  # m/m
min_pos_slope = 0.75  # m/m
max_gap = 15  # m

if __name__ == "__main__":
    os.mkdir(output_path + '/fixed_input')
    Pre_Processing.Input_Make(flowdir_raster, valley_raster,
                              ground_DEM, veg_DEM, int_raster,
                              output_path + '/fixed_input')
    os.mkdir(output_path + '/wet_output')
    Wet_Channel_Extraction.Wet_Channel_Extraction(output_path + '/fixed_input', unit,
                                                  area_thresh, diff_thresh, min_slope,
                                                  max_slope, min_pos_slope, max_gap,
                                                  pix_per_m, ele_per_meter,
                                                  output_path + '/wet_output')
