import rasterio
import csv
from osgeo import gdal_array,gdal,ogr,osr,gdalconst
from numpy import where
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from numpy import linspace,meshgrid,nan,where,isnan,ma,array,log,zeros,mean
from colormaps import plasma,viridis
from matplotlib.patches import Path, PathPatch
import numpy as np
from random import sample
from scipy.spatial.distance import mahalanobis
import os
import sys

rcp = sys.argv[1]+'_'


def av_dist(sp1,sp_pool,v):
	d = []
	for sp in sp_pool:
		d.append(mahalanobis(sp1, sp, v))
	return mean(d)


def func_div(sp_pool,v):
	d = []
	for sp1 in range(len(sp_pool)):
		for sp2 in range(len(sp_pool)):
			if sp1<sp2:
				d.append(mahalanobis(sp_pool[sp1], sp_pool[sp2], v))
	return mean(d)


def write_raster(l,all_locs,outfile,ref = './species/Abies_alba/rasters/Abies_alba_cur.tif'):
	ref_raster = rasterio.open(ref)
	meta = ref_raster.meta
	meta.update(nodata=-9999)
	meta.update(dtype='float32')
	em = zeros(ref_raster.shape)
	em-=9999
	for i in range(len(all_locs)):
		lat,lon = all_locs[i]
		em[lat][lon] = l[i]
	out=rasterio.open(outfile+'.tif', 'w', **meta)
	out.write(em.astype(rasterio.float32),1)
	out.close()
	os.system('gdalwarp '+outfile+'.tif '+outfile+'_wgs84.tif -t_srs "+proj=longlat +ellps=WGS84"')
	os.system('gdal_translate '+outfile+'_wgs84.tif '+outfile+'_wgs84_clipped.tif -projwin -13 71 45 34')


def plot_map(input_raster, output_png, title):
	ds = gdal.Open(input_raster)
	gt = ds.GetGeoTransform()
	proj = ds.GetProjection()
	data = ds.ReadAsArray()
	data = np.flipud(data)
	xres = gt[1]
	yres = gt[5]
	xmin = gt[0] + xres * 0.5
	xmax = gt[0] + (xres * ds.RasterXSize) - xres * 0.5
	ymin = gt[3] + (yres * ds.RasterYSize) + yres * 0.5
	ymax = gt[3] - yres * 0.5
	map_edges = array([[xmin,ymin],[xmax,ymin],[xmax,ymax],[xmin,ymax]])
	data[where(data==-9999)] = nan
	plt.rcParams['hatch.linewidth'] = 0.001
	m = Basemap(resolution='l', projection='cyl',
	            llcrnrlat=ymin, llcrnrlon=xmin,
				urcrnrlat=ymax, urcrnrlon=xmax)
	m.drawcoastlines(linewidth=0.1)
	polys = [p.boundary for p in m.landpolygons]
	polys = [map_edges]+polys[:]
	codes = [[Path.MOVETO] + [Path.LINETO for p in p[1:]] for p in polys]
	polys_lin = [v for p in polys for v in p]
	codes_lin = [c for cs in codes for c in cs]
	path = Path(polys_lin, codes_lin)
	patch = PathPatch(path,facecolor='white', lw=0)
	m.fillcontinents(color='darkgrey', lake_color='darkgrey',zorder=0)  # zorder=0 to paint over continents
	x = linspace(xmin, xmax, data.shape[1])
	y = linspace(ymin, ymax, data.shape[0])
	xx,yy = meshgrid(x, y)
	data_m = ma.masked_where(isnan(data),data)
	vmax = data_m.max()
	vmin = data_m.min()
	im = m.pcolormesh(xx, yy, data_m, cmap='plasma',
					  vmin=vmin,vmax=vmax)
	plt.gca().add_patch(patch)
	plt.colorbar(orientation = 'horizontal',pad=0.01)
	plt.title(title)
	plt.savefig(output_png, dpi=300,bbox_inches='tight')
	clear = [plt.clf() for i in range(10000)]


traits =  [i for i in csv.reader(open('traits.csv','r'))]
func =  [i for i in csv.reader(open('services_new.csv','r'))]
n_func = len(func[0])-1


spp = sorted(list(set([i[0] for i in func[1:]])))
#if 'Olea_europaea' in spp:
#	spp.remove('Olea_europaea')


sp_dict = dict([])
for i in spp:
	t,f = [],[]
	for j in traits:
		if i==j[0]:
			t = map(float,j[1:])
	for j in func:
		if i==j[0]:
			f = map(float,j[1:])
	sp_dict[i] = f+t



if not os.path.exists('./'+rcp+'res_rasters_single'):
	os.makedirs('./'+rcp+'res_rasters_single')


if not os.path.exists('./'+rcp+'res_maps_single'):
	os.makedirs('./'+rcp+'res_maps_single')



#####create lists
suf = ['_cur.tif',
'_'+rcp+'fut1_bin.tif',
'_'+rcp+'fut2_bin.tif',
'_'+rcp+'fut3_bin.tif',
'_ok_'+rcp+'fut1_disp.tif',
'_ok_'+rcp+'fut2_disp.tif',
'_ok_'+rcp+'fut3_disp.tif']


ref_raster = rasterio.open('mask_eu.tif')
ref_raster = ref_raster.read(1)
lat,lon = where(ref_raster>=0)
all_locs = zip(lat,lon)
#len(all_locs)


all_lists = []
for su in suf:
	l = [[] for i in range(len(all_locs))]
	for sp in spp:
		mod = rasterio.open('./species/'+sp+'/rasters/'+sp+su).read(1)
		for i in range(len(all_locs)):
			lat,lon = all_locs[i]
			if mod[lat][lon]>0:
				l[i].append(sp)
	all_lists.append(l)
	print (su)


cur_list,fut1_list,fut2_list,fut3_list,fut1_disp_list,fut2_disp_list,fut3_disp_list = all_lists
no_manag = [fut1_disp_list,fut2_disp_list,fut3_disp_list]

fut_scen = []
fut_c = 0
fl_n = 0
for fut_list in [fut1_list,fut2_list,fut3_list]:
	list_loss = no_manag[fl_n]
	fl_n+=1
	list_rand = []
	for i in range(len(cur_list)):
		list_loss[i] = list(set(list_loss[i])&set(fut_list[i]))
		div0 = len(cur_list[i])
		div1 = len(fut_list[i])
		ext = set(cur_list[i])-set(fut_list[i])
		aliv = set(list_loss[i])#doublecheck, should be ok
		pot = [sp for sp in list(set(fut_list[i])-(aliv|set(cur_list[i])))]
		n = len(ext)
		if n>=len(pot):
			list_rand.append(fut_list[i])
		else:
			list_rand.append(list(aliv)+sample(pot,n))
	###replace maximizing services
	list_best = []
	for serv in range(n_func):
		list_best_ = []
		for i in range(len(cur_list)):
			div0 = len(cur_list[i])
			div1 = len(fut_list[i])
			ext = set(cur_list[i])-set(fut_list[i])
			aliv = set(list_loss[i]) #doublecheck, should be ok
			pot = [[sp_dict[sp][serv],sp] for sp in list(set(fut_list[i])-(aliv|set(cur_list[i])))]
			sy,sn = [j[1] for j in pot if j[0]==1],[j[1] for j in pot if j[0]==0]
			n = len(ext)
			if n>=len(pot):
				list_best_.append(list(aliv|(set(fut_list[i])&set(cur_list[i]))))
			else:
				n1,n2 = len(sy),len(sn)
				if n1>=n:
					list_best_.append(list(aliv)+sample(sy,n))
				else:
					n2 = n-n1
					list_best_.append(list(aliv)+sample(sy,n1)+sample(sn,n2))
		list_best.append(list_best_)
	fut_scen.append([list_loss,list_rand]+list_best)
	fut_c+=1
	print (fut_c)




scen_names = ['fut1','fut2','fut3']
l_names = ['none','rand']+['opt_' + func[0][i+1].replace(' ','_') for i in range(n_func)] #none is also no-dispersal


res_0 = []
for i in range(len(all_locs)):
	sp_pool = array([sp_dict[sp][:n_func] for sp in cur_list[i]])
	if cur_list[i]!=[]:
		res_0.append([len(sp_pool)]+list(array(sp_pool).sum(0)))
	else:
		res_0.append([0.0 for j in range(n_func+1)])


res = []
for l_name in l_names:
	res+=[['cur',l_name]+i for i in res_0]



for scen in range(len(fut_scen)):
	for l in range(len(fut_scen[scen])):
		for i in range(len(all_locs)):
			sp_pool = [sp_dict[sp] for sp in fut_scen[scen][l][i]]
			row=[scen_names[scen],l_names[l]]
			if sp_pool!=[]:
				row += [len(sp_pool)]+[sum([j[k] for j in sp_pool]) for k in range(n_func)]
			else:
				row += [0.0 for j in range(n_func+1)]
			res.append(row)
		print (l_names[l])


names = ['scenario','management','div']+[func[0][i+1] for i in range(n_func)]

out = open(rcp+'results_single.csv','w')
out.write(','.join(names)+'\n')

for i in range(len(res)):
	out.write(','.join(map(str,res[i]))+'\n')


out.close()



#######make_maps

measure_names = ['div']+[func[0][i+1].replace(' ','_') for i in range(n_func)]
for scen_name in scen_names:
	for manag in l_names:
		for measure in range(len(measure_names)):
			to_map = [i[measure+2] for i in res if i[0]==scen_name and i[1]==manag]
			rast_name = './'+rcp+'res_rasters_single/'+scen_name+'_'+manag+'_'+measure_names[measure]
			map_name = './'+rcp+'res_maps_single/'+scen_name+'_'+manag+'_'+measure_names[measure]+'.png'
			map_title = scen_name+'_'+manag+'_'+measure_names[measure]
			write_raster(to_map,all_locs,rast_name,ref = './species/Abies_alba/rasters/Abies_alba_cur.tif')
			plot_map(rast_name+'_wgs84_clipped.tif', map_name, map_title)



for measure in range(len(measure_names)):
	to_map = [i[measure+2] for i in res if i[0]=='cur' and i[1]=='none']
	rast_name = './'+rcp+'res_rasters_single/current_'+measure_names[measure]
	map_name = './'+rcp+'res_maps_single/current_'+measure_names[measure]+'.png'
	map_title = 'current_'+measure_names[measure]
	write_raster(to_map,all_locs,rast_name,ref = './species/Abies_alba/rasters/Abies_alba_cur.tif')
	plot_map(rast_name+'_wgs84_clipped.tif', map_name, map_title)



###delete unused rasters
fff = os.listdir('./'+rcp+'res_rasters_single/')
for f in fff:
	if '_wgs84_clipped.tif' not in f:
		os.remove('./'+rcp+'res_rasters_single/'+f)

