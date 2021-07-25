import arcpy
import Wet_Extraction_Fun
import os


def Wet_Channel_Extraction(input_path, unit, area_thresh, diff_thresh, min_slope, max_slope, min_pos_slope, max_gap,
                           pix_per_m, ele_per_meter, output_path):
    diff_thresh = diff_thresh / ele_per_meter

    os.makedirs(output_path + '/Temp')
    arcpy.env.scratchWorkspace = output_path + '/Temp'
    arcpy.env.workspace = output_path

    arcpy.env.extent = input_path + '/Fdir_8.tif'
    arcpy.env.outputCoordinateSystem = input_path + '/Fdir_8.tif'

    print('4. finding vegetation extent')

    Wet_Extraction_Fun.veg_bound_extraction(input_path + '/veg_DEM.tif', input_path + '/gro_DEM.tif',
                                            output_path + '/No_Veg.tif', diff_thresh)
    Wet_Extraction_Fun.cut_raster(input_path + '/gro_INT.tif', output_path + '/No_veg.tif',
                                  output_path + '/No_Neg_INT.tif')

    print('5. making intensity PDF')
    I_w, I_d = Wet_Extraction_Fun.GMM_fitting(output_path + '/No_Neg_INT.tif')
    edge_diff = (I_w - I_d) / 2

    # With Vegetation Cut
    print('6. With Vegetation Cut')
    Wet_Extraction_Fun.veg_fix(output_path + '/No_veg.tif', output_path + '/veg.tif', area_thresh, pix_per_m)
    Wet_Extraction_Fun.edge_find(output_path + '/No_Neg_INT.tif', output_path + '/veg.tif',
                                 output_path + '/fill_edge.tif', edge_diff, I_d)
    Wet_Extraction_Fun.make_wet_bound(output_path + '/No_Neg_INT.tif', output_path + '/fill_edge.tif',
                                      output_path + '/wet.tif', area_thresh, I_w, pix_per_m)
    Wet_Extraction_Fun.width_cal(output_path + '/wet.tif', pix_per_m, output_path + '/width.tif',
                                 output_path + '/wet_pos_width.tif')
    Wet_Extraction_Fun.depth_cal(output_path + '/wet_pos_width.tif', output_path + '/width.tif',
                                 input_path + '/gro_DEM.tif',
                                 pix_per_m, min_slope, max_slope, min_pos_slope, output_path + '/wet_pos_slope.tif',
                                 output_path + '/water_ele.tif')
    # Without Vegetation Cut
    print('7. Without Vegetation Cut')
    Wet_Extraction_Fun.veg_bound_extraction(input_path + '/veg_DEM.tif', input_path + '/gro_DEM.tif',
                                            output_path + '/veg_All.tif', 500)
    Wet_Extraction_Fun.edge_find(input_path + '/gro_INT.tif', output_path + '/veg_All.tif',
                                 output_path + '/fill_edge_All.tif', edge_diff, I_d)
    Wet_Extraction_Fun.make_wet_bound(input_path + '/gro_INT.tif', output_path + '/fill_edge_All.tif',
                                      output_path + '/wet_All.tif', area_thresh, I_w, pix_per_m)
    Wet_Extraction_Fun.width_cal(output_path + '/wet_All.tif', pix_per_m, output_path + '/width_All.tif',
                                 output_path + '/wet_pos_width_All.tif')
    Wet_Extraction_Fun.depth_cal(output_path + '/wet_pos_width_All.tif', output_path + '/width_All.tif',
                                 input_path + '/gro_DEM.tif',
                                 pix_per_m, min_slope, max_slope, min_pos_slope, output_path + '/wet_pos_slope_All.tif',
                                 output_path + '/water_ele_All.tif')

    print('8. Combine')
    Wet_Extraction_Fun.combine_wet(output_path + '/wet_pos_slope_All.tif', output_path + '/wet_pos_slope.tif', \
                                   output_path + '/width_All.tif', output_path + '/width.tif',
                                   output_path + '/water_ele_All.tif', output_path + '/water_ele.tif',
                                   output_path + '/comb_wet.tif', output_path + '/comb_width.tif',
                                   output_path + '/comb_water_ele.tif')

    Wet_Extraction_Fun.connect_wet_segments(input_path + '/Valley_Network.tif', output_path + '/comb_water_ele.tif', \
                                            output_path + '/comb_wet.tif', input_path + '/Fdir_8.tif',
                                            input_path + '/gro_DEM.tif'
                                            , output_path + '/comb_width.tif', output_path + '/order.tif',
                                            output_path + '/wet_connected.tif'
                                            , max_gap, pix_per_m)
