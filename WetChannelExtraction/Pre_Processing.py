import arcpy
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")

def extraction(bound_file , in_raster ,  out_raster):
    arcpy.sa.ExtractByMask(in_raster, bound_file).save(out_raster)

def bound_make(in_raster ,  bound_file):
    temp_raster = Con(Raster(in_raster) >= 0 , 0)
    arcpy.RasterToPolygon_conversion(temp_raster, bound_file)

def Input_Make(flowdir_raster, valley_raster, ground_DEM , veg_DEM , int_raster, output_file_path):
    
    arcpy.env.extent = flowdir_raster
    arcpy.env.outputCoordinateSystem  = flowdir_raster

    print '1. making bound shape file'
    bound_shapefile = output_file_path + '\\bound.shp'
    bound_make(ground_DEM , bound_shapefile)

    print '2. Copying Valley Data'
    extraction(bound_shapefile , flowdir_raster ,  output_file_path + '\\Fdir_8.tif')
    ## extraction(bound_shapefile , flowacc_raster ,  output_file_path + '\\Facc_8.tif')
    extraction(bound_shapefile , valley_raster  ,  output_file_path + '\\Valley_Network.tif')

    print '3. Copying Intesity Data'
    extraction(bound_shapefile  , veg_DEM ,  output_file_path + '\\veg_DEM.tif')
    extraction(bound_shapefile  ,  ground_DEM ,  output_file_path + '\\gro_DEM.tif')
    extraction(bound_shapefile  , int_raster ,  output_file_path + '\\gro_INT.tif')
