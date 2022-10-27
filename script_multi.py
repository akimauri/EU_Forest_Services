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

rcp = sys.argv[1]+'_'
sleep(randrange(120))

def func_div(sp_pool_func,ext_func,v):
	d = []
	for sp1 in ext_func:
		d_ = []
		for sp2 in sp_pool_func:
			d_.append(mahalanobis(sp1, sp2, v))
		d.append(min(d_))
	return min(d)





def choice_opt(spp_,n,cooc_sc,sp_pool_func,ext_func,v):
	combs = combinations(range(len(spp_)),n)
	comb_res = []
	for co in combs:
		z = zeros(len(spp_[0][1]))
		for j in co:
			z+=spp_[j][1]
		sc0 = z.sum() #tot services
		sc1 = sum(z>0) #!= from sum; total number/diversity of services
		sc2 = sum([cooc_sc[j][1] for j in co])
		pot_func = sp_pool_func[list(co)]
		sc3 = func_div(pot_func,ext_func,v)
		comb_res.append([co,sc0,sc1,sc2,sc3])
	ranks = zeros(len(comb_res))
	for i in range(3):
		ranks+=argsort([j[i+1] for j in comb_res])[::-1]
	ranks+=argsort([j[4] for j in comb_res])
	best = sample(set(where(ranks == min(ranks))[0]),1)[0]
	return ([spp_[i][0] for i in comb_res[best][0]],len(comb_res),comb_res[best][1:])



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




fff = [i for i in os.listdir('./') if rcp+'res_rasters_multi_' in i]
tile = 0
if fff!=[]:
	tile = max([int(i.split('_')[-1]) for i in fff])+1


res_folder = './'+rcp+'res_rasters_multi_'+str(tile)
if not os.path.exists(res_folder):
	os.makedirs(res_folder)


ref_raster = rasterio.open('mask_eu.tif')
ref_raster = ref_raster.read(1)
lat,lon = where(ref_raster>=0)
all_locs = list(zip(lat,lon))[tile*2500:(tile+1)*2500]


traits =  [i for i in csv.reader(open('traits.csv','r'))]
func =  [i for i in csv.reader(open('services_new.csv','r'))]
n_func = len(func[0])-1


spp = sorted(list(set([i[0] for i in func[1:]])))
#out = open('tree_species_list.csv','w')
#for i in enumerate(spp):
#	out.write(','.join(map(str,i))+'\n')
#
#
#out.close()
#if 'Olea_europaea' in spp:
#	spp.remove('Olea_europaea')


sp_dict = dict([])
for i in spp:
	t,f = [],[]
	for j in traits:
		if i==j[0]:
			t = list(map(float,j[1:]))
	for j in func:
		if i==j[0]:
			f = list(map(float,j[1:]))
	sp_dict[i] = f+t



######################COOCCURRENCE
cooc = []
for i in range(len(spp)):
	sp0 = spp[i]
	r0 = rasterio.open('./species/'+sp0+'/rasters/'+sp0+'_cur.tif').read(1)
	r0[where(ref_raster<1)]=0
	r0[where(r0<0)]=0
	for j in range(len(spp)):
		if i<j:
			sp1 = spp[j]
			r1 = rasterio.open('./species/'+sp1+'/rasters/'+sp1+'_cur.tif').read(1)
			r1[where(ref_raster<1)]=0
			r1[where(r1<0)]=0
			ov = (r0*r1).sum()
			co = ov/float(min(r0.sum(),r1.sum()))
			cooc.append([tuple(sorted([sp0,sp1])),co])


cooc_dict = dict(cooc)



#####create lists
suf = ['_cur.tif',
'_'+rcp+'fut1_bin.tif',
'_'+rcp+'fut2_bin.tif',
'_'+rcp+'fut3_bin.tif',
'_ok_'+rcp+'fut1_disp.tif',
'_ok_'+rcp+'fut2_disp.tif',
'_ok_'+rcp+'fut3_disp.tif']


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


#for mahalanobis
data = array([i[n_func:] for i in sp_dict.values()])
v = np.linalg.inv(np.cov(data.T))
#####



###
fut_scen = []
fut_c = 0
att = 0
fl_n = 0
for fut_list in [fut1_list,fut2_list,fut3_list]:
	list_loss = no_manag[fl_n]
	fl_n+=1
	list_best = []
	for i in range(len(cur_list)):
		ext = set(cur_list[i])-set(fut_list[i])
		n = len(ext)
		if n == 0:
			list_best.append([list(set(list_loss[i])|(set(fut_list[i])&set(cur_list[i]))),['na','na','na','na',0]])
		else:
			aliv = set(list_loss[i])&set(fut_list[i])#check
			pot = [[sp,array(sp_dict[sp][:n_func])] for sp in list(set(fut_list[i])-(aliv|set(cur_list[i])))] #
			if n>=len(pot):
				list_best.append([fut_list[i],['na','na','na','na',0]])
			else:
				sp_pool_func = array([sp_dict[sp[0]][n_func:] for sp in pot])
				ext_func  = [sp_dict[sp][n_func:] for sp in ext]
				cooc_sc = []
				for sp0 in pot:
					cooc_sc_ = []
					for sp1 in cur_list[i]:
						cooc_sc_.append(cooc_dict[tuple(sorted([sp0[0],sp1]))])
					cooc_sc.append([sp0[0],mean(cooc_sc_)])
				ch_opt,sc,scores = choice_opt(pot,n,cooc_sc,sp_pool_func,ext_func,v)
				list_best.append([list(aliv)+ch_opt,scores+[len(ch_opt)]])
				att+=sc
	fut_scen.append([list_best])
	fut_c+=1
	print (fut_c)






scen_names = ['fut1','fut2','fut3']
l_names = ['opt']


res_0 = []
for i in range(len(all_locs)):
	sp_pool = array([sp_dict[sp][:n_func] for sp in cur_list[i]])
	sp_ids = '_'.join(map(str,[spp.index(j) for j in cur_list[i]]))
	if cur_list[i]!=[]:
		res_0.append([len(sp_pool)]+list(array(sp_pool).sum(0))+[sp_ids,'na','na','na','na',0])
	else:
		res_0.append(list(zeros(1+n_func))+['','na','na','na','na',0])


res = []
for l_name in l_names:
	res+=[['cur',l_name]+i for i in res_0]

for scen in range(len(fut_scen)):
	for l in range(len(fut_scen[scen])):
		for i in range(len(all_locs)):
			sp_pool = array([sp_dict[sp][:n_func] for sp in fut_scen[scen][l][i][0]])
			sp_ids = '_'.join(map(str,[spp.index(j) for j in fut_scen[scen][l][i][0]]))
			row=[scen_names[scen],l_names[l]]
			if fut_scen[scen][l][i][0]!=[]:
				row += [len(sp_pool)]+list(sp_pool.sum(0))+[sp_ids]+fut_scen[scen][l][i][1]
			else:
				row += [0.0 for j in range(n_func+1)]+[sp_ids]+fut_scen[scen][l][i][1]
			res.append(row)
		print (l_names[l])


names = ['scenario','management','div']+[func[0][i+1] for i in range(n_func)]+['local_species_pool','sc0','sc1','sc2','sc3','replaced']

out = open(rcp+'results_multi_'+str(tile)+'.csv','w')
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
			rast_name = res_folder+'/'+scen_name+'_'+manag+'_'+measure_names[measure]
			write_raster(to_map,all_locs,rast_name,ref = './species/Abies_alba/rasters/Abies_alba_cur.tif')



###delete unused rasters
fff = os.listdir(res_folder)
for f in fff:
	if '_wgs84_clipped.tif' not in f:
		os.remove(res_folder+'/'+f)



