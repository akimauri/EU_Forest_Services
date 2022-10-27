import rasterio
import csv
from osgeo import gdal_array,gdal,ogr,osr,gdalconst
from numpy import where
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from numpy import linspace,meshgrid,nan,where,isnan,ma,array,log,zeros,mean
import rasterio
from matplotlib.patches import Path, PathPatch
import numpy as np
from random import sample
from scipy.spatial.distance import mahalanobis
import os



spp = os.listdir('./')
spp.remove('preprocessing_rasters.py')
spp.sort()

rcps = ['rcp45_','rcp85_']


for rcp in rcps:
	suf = ['_cur.tif','_'+rcp+'fut1_bin.tif',
	'_'+rcp+'fut2_bin.tif',
	'_'+rcp+'fut3_bin.tif',
	'_'+rcp+'fut1_disp.tif',
	'_'+rcp+'fut2_disp.tif',
	'_'+rcp+'fut3_disp.tif']
	referencefile = './'+spp[0]+'/rasters/'+spp[0]+suf[0]#Path to reference file
	reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
	referenceProj = reference.GetProjection()
	referenceTrans = reference.GetGeoTransform()
	x = reference.RasterXSize
	y = reference.RasterYSize
	ref_raster = rasterio.open('./'+spp[0]+'/rasters/'+spp[0]+suf[0])
	meta = ref_raster.meta
	meta.update(nodata=-9999)
	meta.update(dtype='int16')
	rast = rcp+'future_disp_full.tif'
	for sp in spp:
		inputfile = './'+sp+'/rasters/'+sp+'_'+rast
		input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
		inputProj = input.GetProjection()
		inputTrans = input.GetGeoTransform()
		outputfile = './'+sp+'/rasters/'+sp+'_ok_'+rast#Path to output file
		driver= gdal.GetDriverByName('GTiff')
		bandreference = input.GetRasterBand(1)
		output = driver.Create(outputfile,x,y,1,bandreference.DataType)
		output.SetGeoTransform(referenceTrans)
		output.SetProjection(referenceProj)
		gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_NearestNeighbour)
		del output
		r0 = rasterio.open('./'+sp+'/rasters/'+sp+'_cur.tif').read(1)
		r0[where(ref_raster<1)]=0
		r0[where(r0<1)]=0
		#fut1
		r1 = rasterio.open('./'+sp+'/rasters/'+sp+suf[1]).read(1) #fut1
		r1[where(r1<1)]=0
		mod = rasterio.open('./'+sp+'/rasters/'+sp+'_ok_'+rast).read(1)
		#mod[where(mod==-32768)]=0
		mod[where(mod>=30000)]=0
		mod[where(mod==-301)]=1 #decolonizzato step 3
		mod[where(mod==-201)]=1 #decolonizzato step 2
		mod[where((mod>1)*(mod<201))]=1
		mod[where(mod==-101)]=0 #decolonizzato step 1
		mod[where((mod>=201)*(mod<301))]=0
		mod[where((mod>=301)*(mod<401))]=0
		mod[where((r0>0)*(r1>0))] = 1
		out=rasterio.open('./'+sp+'/rasters/'+sp+'_ok_'+rcp+'fut1_disp.tif', 'w', **meta)
		out.write(mod.astype(rasterio.int16),1)
		out.close()
		mod1 = mod.copy()
		#fut2
		r1 = rasterio.open('./'+sp+'/rasters/'+sp+suf[2]).read(1) #fut1
		r1[where(r1<1)]=0
		mod = rasterio.open('./'+sp+'/rasters/'+sp+'_ok_'+rast).read(1)
		mod[where(mod==-32768)]=0
		mod[where(mod>=30000)]=0
		mod[where(mod==-301)]=1 #decolonizzato step 3
		mod[where((mod>1)*(mod<201))]=1
		mod[where(mod==-101)]=0 #decolonizzato step 1
		mod[where((mod>=201)*(mod<301))]=1
		mod[where(mod==-201)]=0 #decolonizzato step 2
		mod[where((mod>=301)*(mod<401))]=0
		mod[where((mod1>0)*(r1>0))] = 1
		out=rasterio.open('./'+sp+'/rasters/'+sp+'_ok_'+rcp+'fut2_disp.tif', 'w', **meta)
		out.write(mod.astype(rasterio.int16),1)
		out.close()
		mod2 = mod.copy()
		#fut3
		r1 = rasterio.open('./'+sp+'/rasters/'+sp+suf[3]).read(1) #fut1
		r1[where(r1<1)]=0
		mod = rasterio.open('./'+sp+'/rasters/'+sp+'_ok_'+rast).read(1)
		mod[where(mod==-32768)]=0
		mod[where(mod>=30000)]=0
		mod[where((mod>1)*(mod<201))]=1
		mod[where(mod==-101)]=0 #decolonizzato step 1
		mod[where((mod>=201)*(mod<301))]=1
		mod[where(mod==-201)]=0 #decolonizzato step 2
		mod[where((mod>=301)*(mod<401))]=1
		mod[where(mod==-301)]=0 #decolonizzato step 3
		mod[where((mod2>0)*(r1>0))] = 1
		out=rasterio.open('./'+sp+'/rasters/'+sp+'_ok_'+rcp+'fut3_disp.tif', 'w', **meta)
		out.write(mod.astype(rasterio.int16),1)
		out.close()
		mod3 = mod.copy()
		print sp



