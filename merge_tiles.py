
"""
Created on Wed Dec  1 15:53:38 2021

@author: stronagi
"""
import rasterio
import csv
from numpy import where,array,zeros,mean
import numpy as np
from random import sample
from scipy.spatial.distance import mahalanobis
import os
from itertools import combinations
from numpy import argsort
from time import sleep
from random import randrange
import sys


rcp = 'rcp85_'

res_folder = './'+rcp+'res_rasters_multi_merged/'
if not os.path.exists(res_folder):
	os.makedirs(res_folder)



ref_folder = './'+rcp+'res_rasters_multi_0'
fff = os.listdir(ref_folder)
ref_raster = rasterio.open(ref_folder+'/'+fff[0])
meta = ref_raster.meta
ref_raster = ref_raster.read(1)

for ff in fff:
	em = zeros(ref_raster.shape)
	for tile in range(24):
		f = 	rasterio.open('./'+rcp+'res_rasters_multi_'+str(tile)+'/'+ff).read(1)
		f[where(f==-9999)] = 0
		em+=f
	out=rasterio.open(res_folder+ff, 'w', **meta)
	out.write(em.astype(rasterio.float32),1)
	out.close()
	print (ff)


out = open('rcp85_results_multi.csv','w')
for tile in range(24):
	f = open('rcp85_results_multi_'+str(tile)+'.csv','r')
	if tile!=0:
		head = next(f)
	for i in f:
		out.write(i)


out.close()










