import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np
import scipy
from scipy import ndimage, sparse
import math
import os
import sklearn
from sklearn import mixture
import matplotlib.pyplot as pl


arcpy.CheckOutExtension("Spatial")

def veg_bound_extraction(veg_DEM , gro_DEM , no_veg_raster , diff_thresh):
    DIF_RASTER = Raster(veg_DEM) - Raster(gro_DEM)
    NO_VEG = Con(DIF_RASTER < diff_thresh , 1 , 0)
    NO_VEG.save(no_veg_raster)

def cut_raster(in_raster, bound_raster , out_raster):
    Con(Raster(bound_raster)  == 1 , Raster(in_raster)).save(out_raster)

def GMM_fitting(int_raster):
    data = arcpy.RasterToNumPyArray(int_raster , nodata_to_value=0)
    data.reshape((1 , data.shape[0] *  data.shape[1]))
    data = data[data > 0]
    data = data[data < 250]

    # Maximum Likelihood
    s_data = shrink_data(data , 5000)
    fitting_result = prob_fit(s_data)

    x = np.arange(np.max(data))
    p_all = np.zeros_like(x)
    p = np.zeros((3 , x.shape[0]))
    for i in range (0 , 3):
        n_dist = scipy.stats.norm(fitting_result[1 , i], fitting_result[2 , i])
        
        p[i , :] = n_dist.pdf(x) * fitting_result[0 , i]
        p_all = p_all + p[i , :]
        
    flage_dry = 0
    flage_wet = 0
    for i in range (0 , x.shape[0]):
        if p[2 , i] <= p[1 , i] and flage_wet == 0:
            flage_wet = 1
            print 'wet trasition = ' , x[i]
            wet_tras = x[i]
        if p[1 , i] <= p[0 , i] and flage_dry == 0:
            flage_dry = 1
            print 'dry trasition = ' , x[i]
            dry_tras = x[i]
    return wet_tras , dry_tras

def prob_fit(int_Array):
    gmix = mixture.GMM(n_components = 3, covariance_type='diag' , n_init  = 20 , n_iter = 15000 , tol = 0.00000001)
    gmix.fit(int_Array)
    fitting_result = np.zeros((3 , 3))

    for i in range(0 , 3):
        fitting_result[0 , i] = gmix.weights_[np.argmax(gmix.means_)]
        fitting_result[1 , i] = gmix.means_[np.argmax(gmix.means_)]
        fitting_result[2 , i] = gmix.covars_[np.argmax(gmix.means_)] ** 0.5
        gmix.means_[np.argmax(gmix.means_)] = 0
    
    return fitting_result

def shrink_data(data , new_size):
    bins = range(0 , np.max(data) , 1)
    hist , xx = np.histogram(data , bins)

    hist = (np.round(hist.astype(np.float32) / data.shape[0] * new_size , 0)).astype(np.uint16)

    pl.plot(xx[0:xx.size - 1] , hist.astype(np.float32) / np.sum(hist))
    
    s_data = []
    for i in range (0 , len(bins) - 1):
        for j in range (0 , hist[i]):
            s_data.append(xx[i])
            
    s_data = np.asarray(s_data).reshape((len(s_data),1))

    return s_data

def veg_fix(veg_raster , out_veg_raster , area_thresh , pix_per_m):
    corner = arcpy.Point(arcpy.Describe(veg_raster).Extent.XMin,arcpy.Describe(veg_raster).Extent.YMin)
    dx = arcpy.Describe(veg_raster).meanCellWidth

    veg = arcpy.RasterToNumPyArray(veg_raster  , nodata_to_value=0)
    bin_veg , _ = delet_small_area(veg, area_thresh , pix_per_m)
    arcpy.NumPyArrayToRaster(bin_veg , corner, dx, dx , value_to_nodata=0).save(out_veg_raster)

def edge_find(int_raster , veg_raster , edge_raster ,  edge_diff , dry_thresh):

    corner = arcpy.Point(arcpy.Describe(int_raster).Extent.XMin,arcpy.Describe(int_raster).Extent.YMin)
    dx = arcpy.Describe(int_raster).meanCellWidth

    veg_Array = arcpy.RasterToNumPyArray(veg_raster , nodata_to_value= 0).astype(np.uint8)
    veg_Array = scipy.ndimage.morphology.binary_erosion(veg_Array , structure = np.ones((3 , 3))).astype(veg_Array.dtype)

    org_inte = arcpy.RasterToNumPyArray(int_raster  , nodata_to_value=0)
   
    # Find edge
    for i in range (0 , 2):
        # find the max of a moving 3*3 window
        max_raster = FocalStatistics(int_raster , NbrRectangle( 3 , 3 , "CELL") , "MAXIMUM" , "DATA")
        max_inte = arcpy.RasterToNumPyArray(max_raster  , nodata_to_value=0)
        max_raster = None

        inte = arcpy.RasterToNumPyArray(int_raster  , nodata_to_value=0)
        # find max - int >= thresh
        if i == 0 :
            edge = np.zeros_like(max_inte).astype(np.uint8)
        edge = np.where((max_inte - inte) >= edge_diff , 1 , edge)
        edge = np.where(veg_Array == 0 , 0 , edge).astype(edge.dtype)
        edge = np.where(org_inte >= dry_thresh , 0 , edge).astype(edge.dtype)
        # replace the edges with max
        inte = np.where(edge == 1 , max_inte , inte)
        int_raster = arcpy.NumPyArrayToRaster(inte , corner,dx,dx , value_to_nodata=0)

    arcpy.NumPyArrayToRaster(edge , corner,dx,dx , value_to_nodata=0).save(edge_raster)

def make_wet_bound(int_raster , edge_raster , out_wet_raster,  area_thresh, int_thresh , pix_per_m):
    corner = arcpy.Point(arcpy.Describe(int_raster).Extent.XMin,arcpy.Describe(int_raster).Extent.YMin)
    dx = arcpy.Describe(int_raster).meanCellWidth

    inte = arcpy.RasterToNumPyArray(int_raster  , nodata_to_value=0)
    edge = arcpy.RasterToNumPyArray(edge_raster  , nodata_to_value=0).astype(np.uint8)
    
    wet = np.where((inte <= int_thresh) & (inte > 0) , 1, 0)

    wet = np.where((wet == 1) | (edge == 1) , 1, 0)
    bin_wet , area_wet = delet_small_area(wet, area_thresh , pix_per_m)
    arcpy.NumPyArrayToRaster(bin_wet , corner, dx, dx , value_to_nodata=0).save(out_wet_raster)
    

def width_cal(wet_raster, pix_per_m , width_raster , wet_pos_width_raster):

    corner = arcpy.Point(arcpy.Describe(wet_raster).Extent.XMin,arcpy.Describe(wet_raster).Extent.YMin)
    dx = arcpy.Describe(wet_raster).meanCellWidth

    org_wet = arcpy.RasterToNumPyArray(wet_raster  , nodata_to_value=0)
    wet = ndimage.morphology.binary_fill_holes(org_wet).astype(np.uint8)

    er_wet_1 = ndimage.binary_erosion(wet).astype(np.uint8)
    bound_1 = np.where((wet == 1) & (er_wet_1 == 0), 1, 0).astype(np.uint8)

    Lab_wet, num_label = ndimage.label(wet , structure = np.ones((3 , 3)))

    labels = np.arange(1, num_label + 1)
    area_label = ndimage.labeled_comprehension(wet, Lab_wet, labels, np.sum, float, 0) *  pix_per_m ** 2
    per_label = ndimage.labeled_comprehension(bound_1, Lab_wet, labels, np.sum, float, 0) * pix_per_m

    len_label = np.zeros_like(area_label)
    wid_label = np.zeros_like(area_label)
    for i in range (0 , len(area_label)):
        if (-1 * per_label[i] / 2) ** 2 - 4 * 1 *  area_label[i] > 0:
            temp_root = np.roots([1 , -1 * per_label[i] / 2 , area_label[i]])
            len_label[i] = np.max(temp_root)
            wid_label[i] = np.min(temp_root)
        else:
            len_label[i] = 0
            wid_label[i] = 0
            
    width = np.zeros_like(wet).astype(np.float32)
    width = wid_label[Lab_wet - 1].astype(np.float32)
    width = np.where(Lab_wet == 0 , 0 , width)
    arcpy.NumPyArrayToRaster(width , corner,dx,dx , value_to_nodata=0).save(width_raster)
    wet = np.where((org_wet == 1) & (width > 0) , 1 , 0).astype(np.uint8)
    arcpy.NumPyArrayToRaster(wet , corner,dx,dx , value_to_nodata=0).save(wet_pos_width_raster)
    del wet , width

def depth_cal(wet_raster , width_raster , ele_raster, pix_per_m , min_slope , max_slope , min_pos_slope , wet_pos_depth_raster , water_ele_raster ):

    corner = arcpy.Point(arcpy.Describe(wet_raster).Extent.XMin,arcpy.Describe(wet_raster).Extent.YMin)
    dx = arcpy.Describe(wet_raster).meanCellWidth

    org_wet = arcpy.RasterToNumPyArray(wet_raster  , nodata_to_value=0)
    wet = ndimage.morphology.binary_fill_holes(org_wet).astype(np.uint8)
    ele = arcpy.RasterToNumPyArray(ele_raster  , nodata_to_value=0)

    commen_border = find_sep_border(org_wet , 7 , corner , dx)
    wet = ndimage.binary_dilation(org_wet , iterations = 2).astype(np.uint8)
    wet = np.where(commen_border == 1 , 0 , wet)
    commen_border = find_sep_border(wet , 7 , corner , dx)
    wet = ndimage.morphology.binary_fill_holes(wet).astype(np.uint8)
    wet = np.where(commen_border == 1 , 0 , wet)
    wet = np.where(org_wet == 1 , 1 , wet)

    er_wet_1 = ndimage.binary_erosion(wet).astype(np.uint8)
    bound_1 = np.where((wet == 1) & (er_wet_1 == 0), 1, 0).astype(np.uint8)
    ele_1 = np.where(bound_1 == 1 , ele , 0)
    temp_ele_1 = arcpy.NumPyArrayToRaster(ele_1 , corner,dx,dx , value_to_nodata=0)
    ele_1_thick = FocalStatistics(temp_ele_1 , NbrRectangle( 3 , 3 , "CELL") , "MEAN" , "DATA")

    del ele_1

    er_wet_2 = ndimage.binary_erosion(er_wet_1).astype(np.uint8)
    bound_2 = np.where((er_wet_1 == 1) & (er_wet_2 == 0), 1, 0).astype(np.uint8)
    ele_2 = np.where(bound_2 == 1 , ele , 0)
        
    ele_1 = arcpy.RasterToNumPyArray(ele_1_thick  , nodata_to_value=0)
    ele_1 = np.where(bound_2 == 0 , 0 , ele_1)

    slope = (ele_1 - ele_2) / pix_per_m
    slope_sign = np.where(slope > 0 , 1 , 0)

    Lab_wet, num_label = ndimage.label(wet , structure = np.ones((3 , 3)))
    labels = np.arange(1, num_label + 1)
    slope_label = ndimage.labeled_comprehension(slope , Lab_wet, labels, np.sum, float, 0)
    slope_sign_label = ndimage.labeled_comprehension(slope_sign , Lab_wet, labels, np.sum, float, 0)
    count_2 = ndimage.labeled_comprehension(bound_2 , Lab_wet, labels, np.sum, float, 0)
    water_ele_label = ndimage.labeled_comprehension(ele , Lab_wet, labels, np.min, float, 0)

    slope_seg = np.zeros_like(org_wet).astype(np.float32)
    slope_seg = slope_label[Lab_wet - 1].astype(np.float32) / count_2[Lab_wet - 1].astype(np.float32)
    slope_seg = np.where((Lab_wet == 0) | (org_wet == 0), 0 , slope_seg)
   
    new_wet = np.where((slope_seg > min_slope) & (slope_seg < max_slope), 1 , 0).astype(np.uint8)

    width = arcpy.RasterToNumPyArray(width_raster  , nodata_to_value=0)
    depth = slope_seg * width / 2

    del slope_seg , slope_label, width

    slope_sign_seg = np.zeros_like(org_wet).astype(np.float32)
    slope_sign_seg = slope_sign_label[Lab_wet - 1].astype(np.float32) / count_2[Lab_wet - 1].astype(np.float32)
    slope_sign_seg = np.where((Lab_wet == 0) | (org_wet == 0), 0 , slope_sign_seg)

    new_wet = np.where((slope_sign_seg <= min_pos_slope), 0 , new_wet).astype(np.uint8)
    del slope_sign_seg , slope_sign_label , count_2
    
    arcpy.NumPyArrayToRaster(new_wet , corner,dx,dx , value_to_nodata=0).save(wet_pos_depth_raster)

    water_ele = water_ele_label[Lab_wet - 1].astype(np.float32) + depth
    water_ele = np.where((Lab_wet == 0) | (new_wet == 0), 0 , water_ele)
    arcpy.NumPyArrayToRaster(water_ele , corner,dx,dx , value_to_nodata=0).save(water_ele_raster)

    del new_wet , water_ele , depth


def combine_wet(wet_raster_1 , wet_raster_2 , width_raster_1 , width_raster_2 , water_ele_raster_1 , \
                water_ele_raster_2 , wet_raster_out , width_raster_out, water_ele_raster_out):
    corner = arcpy.Point(arcpy.Describe(wet_raster_1).Extent.XMin,arcpy.Describe(wet_raster_1).Extent.YMin)
    dx = arcpy.Describe(wet_raster_1).meanCellWidth

    wet_1 = arcpy.RasterToNumPyArray(wet_raster_1  , nodata_to_value=0).astype(np.uint8)
    wet_2 = arcpy.RasterToNumPyArray(wet_raster_2  , nodata_to_value=0).astype(np.uint8)
    wet_out = np.where((wet_1 == 1) | (wet_2 == 1) , 1 , 0)
    arcpy.NumPyArrayToRaster(wet_out , corner,dx,dx , value_to_nodata = 0).save(wet_raster_out)
    del wet_out
    
    width_1 = arcpy.RasterToNumPyArray(width_raster_1  , nodata_to_value=0)
    width_2 = arcpy.RasterToNumPyArray(width_raster_2  , nodata_to_value=0)
    width_out = np.where((wet_1 == 1) , width_1 , width_2)
    arcpy.NumPyArrayToRaster(width_out , corner,dx,dx , value_to_nodata = 0).save(width_raster_out)
    del width_1 , width_2 , width_out 
    
    water_ele_1 = arcpy.RasterToNumPyArray(water_ele_raster_1  , nodata_to_value=0)
    water_ele_2 = arcpy.RasterToNumPyArray(water_ele_raster_2  , nodata_to_value=0)
    water_ele_out = np.where((wet_1 == 1) , water_ele_1 , water_ele_2)
    arcpy.NumPyArrayToRaster(water_ele_out , corner,dx,dx , value_to_nodata = 0).save(water_ele_raster_out)
    del water_ele_1 , water_ele_2 , water_ele_out

def connect_wet_segments(valley_raster , water_ele_raster , wet_raster , flowdir_raster , ele_raster , \
                         width_raster , order_raster, connected_wet_raster, max_gap, pix_per_m):


    corner = arcpy.Point(arcpy.Describe(valley_raster).Extent.XMin,arcpy.Describe(valley_raster).Extent.YMin)
    dx = arcpy.Describe(valley_raster).meanCellWidth

    valley = arcpy.RasterToNumPyArray(valley_raster  , nodata_to_value=0).astype(np.uint8)
    dep_ele = arcpy.RasterToNumPyArray(water_ele_raster  , nodata_to_value=0)
    wet = arcpy.RasterToNumPyArray(wet_raster  , nodata_to_value=0).astype(np.uint8)
    ele = arcpy.RasterToNumPyArray(ele_raster  , nodata_to_value=0)
    width = arcpy.RasterToNumPyArray(width_raster  , nodata_to_value=0)
    width = np.where(wet == 1 , width , 0)
    dep_ele = np.where(wet == 1 , dep_ele , 0)

    Lab_wet, num_label_wet = ndimage.label(wet , structure = np.ones((3 , 3)))


    # find dry and wet valley
    tran = disconnect_valley_segments(valley_raster , flowdir_raster , 3 , order_raster)

    dry_valley = np.where((valley == 1) & (wet == 0) , 1 , 0).astype(np.uint8)
    wet_valley = np.where((valley == 1) & (wet == 1) , 1 , 0).astype(np.uint8)

    temp_dry_valley = ndimage.binary_dilation(dry_valley , iterations = 1 , structure = np.ones((3 , 3))).astype(np.uint8)
    dry_valley = np.where(((wet_valley == 1) & (temp_dry_valley == 1)) | (dry_valley == 1), 1 , 0).astype(np.uint8)
    dry_valley = np.where(tran == 1, 0 , dry_valley).astype(np.uint8)
    del temp_dry_valley 

    Lab_dry, num_label_dry = ndimage.label(dry_valley , structure = np.ones((3 , 3)))
    labels = np.arange(1, num_label_dry + 1)

    # 1. connection of dry segments
    num_conn_label = ndimage.labeled_comprehension(Lab_wet, Lab_dry, labels, count_unique , int, 0)
    num_conn = np.zeros_like(wet)
    num_conn = num_conn_label[Lab_dry - 1]
    num_conn = np.where(Lab_dry == 0 , 0 , num_conn)
    conn_dry_valley = np.where(num_conn >= 3 , 1 , 0)
    del num_conn, num_conn_label

    # 2. elevation of dry segments
    Lab_dry, num_label_dry = ndimage.label(conn_dry_valley , structure = np.ones((3 , 3)))
    labels = np.arange(1, num_label_dry + 1)
    ave_ele_label = ndimage.labeled_comprehension(dep_ele, Lab_dry, labels, np.max , float, 0)
    max_ele_label = ndimage.labeled_comprehension(ele, Lab_dry, labels, np.max , float, 0)
    diff_ele = np.zeros_like(ele)
    diff_ele = ave_ele_label[Lab_dry - 1] - max_ele_label[Lab_dry - 1]
    diff_ele = np.where(Lab_dry == 0 , 0 , diff_ele)
    conn_dry_valley = np.where(diff_ele < 0 , 0 , conn_dry_valley)

    del diff_ele, max_ele_label , ave_ele_label


    # 3. length of dry segment
    Lab_dry, num_label_dry = ndimage.label(conn_dry_valley , structure = np.ones((3 , 3)))
    labels = np.arange(1, num_label_dry + 1)
    area_gap_label = ndimage.labeled_comprehension(conn_dry_valley, Lab_dry, labels, np.sum, float, 0) *  pix_per_m ** 2
    area_dry = np.zeros_like(ele)
    area_dry = area_gap_label[Lab_dry - 1]
    area_dry = np.where((Lab_dry == 0) , 0 , area_dry)
    conn_dry_valley = np.where(area_dry > max_gap , 0 , conn_dry_valley)
    del area_dry , ele
    
    # find width of added segments
    Lab_dry, num_label_dry = ndimage.label(conn_dry_valley , structure = np.ones((3 , 3)))
    labels = np.arange(1, num_label_dry + 1)
    width_gap_label = ndimage.labeled_comprehension(width, Lab_dry, labels, np.sum , float, 0)
    num_conn_label = ndimage.labeled_comprehension(Lab_wet, Lab_dry, labels, count_unique , int, 0)
    width_conn = np.zeros_like(width)
    width_conn = (width_gap_label[Lab_dry - 1] /  num_conn_label[Lab_dry - 1]).astype(np.float32)
    width_conn = np.where(Lab_dry == 0 , 0 , width_conn)

    add_wet = add_width(width_conn)
    wet = np.where((add_wet == 1) | (wet == 1) , 1 , 0)

    arcpy.NumPyArrayToRaster(wet , corner,dx,dx).save(connected_wet_raster)


def disconnect_valley_segments(valley , flowdir , size , order):
    stream_order(valley , flowdir , order)
    orde = arcpy.RasterToNumPyArray(order  , nodata_to_value=0).astype(np.uint8)
    temp_raster = FocalStatistics(order , NbrRectangle( size , size , "CELL") , "VARIETY" , "DATA")
    tran_1 = arcpy.RasterToNumPyArray(temp_raster, nodata_to_value=0)
    temp_raster = None
    temp_raster = FocalStatistics(order , NbrRectangle( size , size , "CELL") , "MAXIMUM" , "DATA")
    tran_2 = arcpy.RasterToNumPyArray( temp_raster, nodata_to_value=0)
    temp_raster = None
    tran = np.where((tran_1 >= 2) & (tran_2 != orde) , 1 , 0)
    del tran_1 , tran_2
    return tran

def count_unique(a):
    c = np.unique(a).shape[0]
    return int(c)

def add_width(a):
    b = np.where(a > 0 , 1 , 0).astype(np.uint8)
    number_row = b.shape[0]
    number_col = b.shape[1]
    cx = sparse.coo_matrix(a)
    for i , j , v in zip(cx.row, cx.col, cx.data):
        R = int(round(v /2,1))
        for ii in range (-R , R):
            for jj in range (-R , R):
                new_i = ii + i 
                new_j = jj + j
                if bond_check(new_i , 0 , number_row) == 1 and bond_check(new_j , 0 , number_col):
                    b[new_i , new_j] = 1
    return b

def find_sep_border(a , size , corner , dx):
    Lab, _ = ndimage.label(a , structure = np.ones((3 , 3)))
    comm_raster = arcpy.NumPyArrayToRaster(Lab , corner , dx , dx , value_to_nodata=0)
    del Lab
    commen_border = arcpy.RasterToNumPyArray(FocalStatistics(comm_raster , NbrRectangle( size , size , "CELL") , "VARIETY" , "DATA") , nodata_to_value=0)
    commen_border = np.where(commen_border >= 2 , 1 , 0).astype(np.uint8)
    return commen_border

def bond_check(i , min_i , max_i):
    flage = 0
    if i >= min_i and i < max_i:
        flage = 1
    return flage

def delet_small_area(A, area_thresh , pix_per_m):
    Lab_Array, num_label = ndimage.label(A , structure = np.ones((3 , 3)))
    labels = np.arange(1, num_label + 1)
    area_label = ndimage.labeled_comprehension(A, Lab_Array, labels, np.sum, int, 0) * pix_per_m ** 2
    area = np.zeros_like(A).astype(np.float32)
    for index , value in np.ndenumerate(Lab_Array):
        if value > 0:
            area[index] = area_label[int(value) - 1]
    bin_A = np.where((area > area_thresh) & (A > 0) , 1, 0).astype(np.uint8)
    return bin_A , area

def stream_order(st_raster , dir_raster , out_raster):    
    StreamOrder(st_raster , dir_raster).save(out_raster)
