## WetChannelExtraction
This is a Python code to extract wet segments using high resolution Digital Elevation Models (DEMs) and Intensity raster.

## A. Installation

The code itself does not require installation; however, there are some packages that are required to be installed before running the code.

1. “arcpy” library: this library is imbedded in ArcGIS package. If you install ArcGIS in your computer, the “arcpy” library and Python 2.7 will be installed automatically. It is recommended to use 64-bit python which is available on new versions of ArcGIS. 64-bit background geoprocessing is also available for older versions at http://resources.arcgis.com/en/help/main/10.1/index.html#/Background_Geoprocessing_64_bit/002100000040000000/.

2. •	“numpy”, “scipy”, “Scikit-learn”, and “matplotlib” libraries: 64-bit version of these libraries can be found at (http://www.lfd.uci.edu/~gohlke/pythonlibs/). Make sure you download the libraries for Python 2.7. The whl files can be installed using pip as illustrated in this video (https://www.youtube.com/watch?v=zPMr0lEMqpo).


## B. Code structure

1. “Run.py” is the file for setting the parameters and running the code.

2. “Wet_Channel_Extraction.py” is the main code which calls the functions from “Wet_Extraction_Fun.py” to delineate wet segments.

3. “Pre_Processing.py” contains the functions to fix the extent of the input rasters.


## C. Inputs

Before running the code, there are some parameters which should be set in “Run.py” including.

1. The path to the bare-earth DEM (ground_DEM).

2. The path to the vegetation DEM which is created from all LiDAR points (veg_DEM)

3. The path to the intensity  raster (int_raster)

4. The path to the flow direction grid (flowdir_raster) and valley network (valley_raster). These two rasters can be extracted using the the code available at https://github.com/MiladHooshyar/DrainageNetworkExtraction.git and are saved as “Modified_Fdir_8.tif” and “valley_network.tif”, respectively.

5. The unit of the DEM (unit). It is ‘m’ for meter and ‘ft’ for feet.

6. The size of the DEM in meters (pix_per_m).

7. The elevation per meter (ele_per_meter). It is 1 if the unit is meter, 0.3045 if the unit is feet.

8. The minimum area of wet segments (area_thresh).

9. The elevation different between the top of dense vegetation and land surface (diff_thresh).

10. Slope thresholds for differentiation erroneous wet segments (min_slope, max_slope , and min_pos_slope).

11. The maximum distance between two erroneously disconnected segment (max_gap).


After setting the parameters, one can execute “Run.py” to extract the wet segments. The output files will be saved in a folder called “wet_output” in the specified Output folder path. “comb_wet.tif” and “wet_connected.tif” are original and connected wet segments, respectively.

## D. Publications

1. Hooshyar, M., S. Kim, D. Wang, and S. C. Medeiros (2015), Wet channel network extraction by integrating LiDAR intensity and elevation data, Water Resour. Res., 51, 10029–10046, doi:10.1002/2015WR018021.
