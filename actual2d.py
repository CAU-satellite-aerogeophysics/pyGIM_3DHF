"""
Georg HÃ¼ttner - September 2024

This script performs two dimensional thermal modelling from provided lithospheric
interfaces. Based on this a surface geothermal heat flow can be calculated.
"""

import pygimli as pg
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib import ticker
import cartopy.crs as ccrs
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter)
import datetime
import pickle
import gepy

#%%___General parameters______________________________

area_coords = [1300000,1400000,-275000,-375000]
# profile_coord=([1300000,1400000],[-275000,-375000])

# profile_coord=([1300000,1400000],[-320000,-300000]) # LV prof 1
# profile_coord=([1300000,1400000],[-340000,-320000]) # LV prof 2
profile_coord=([1325000,1365000],[-275000,-375000]) # LV prof 3

# area_coords = [1355000,1435000,-855000,-985000] # Dome C

# profile_coord=([1405000,1435000],[-985000,-855000]) # DC prof 1
# profile_coord=([1360000,1415000],[-985000,-855000]) # DC prof 2
# profile_coord=([1355000,1420000],[-920000,-985000]) # DC prof 3
# profile_coord=([1355000,1435000],[-905000,-970000]) # DC prof 3
# 
# area_coords = [1800000,1900000,-650000,-750000] # ASB

# profile_coord=([1900000,1800000],[-650000,-750000]) # ASB prof 1
# profile_coord=([1800000,1900000],[-750000,-650000]) # ASB prof 1
# profile_coord=([1800000,1900000],[-700000,-700000]) # ASB prof 2

#_____Layers__________________________________________

layers = 4
do_sediments = True
do_hp = True

#_____Data Input______________________________________

data_input = []

data_input.append(np.load("data_example/topography.npz"))
data_input.append(np.load("data_example/Sediments.npz"))
data_input.append(np.load("data_example/Moho.npz"))
data_input.append(np.load("data_example/LAB.npz"))
data_input.append(np.load("data_example/HP.npz"))

data_resolution = [1000,9250,10000,10000,12500]

#%%___Thermal parameters______________________________

k0 = 1.5
k1 = 2.7
k2 = 3.5
hp_0 = 0.000001
# If you have input HP data: do not use hp_0 in any function.
# If you do not have input HP data: set a hp_0 value and give it in 

tc = [k0,k1,k2]

#%%___cut data to size & interpolate along profile____

data_list = []

for i,layer in enumerate(data_input):
    data = gepy.cut_and_interpolate_2d(layer, data_resolution, i, layers,area_coords,
                                        do_sediments=do_sediments, do_hp=do_hp, hp_0=hp_0)
    data_list.append(data)

if do_sediments:
    data_list[1,1] = data_list[0,1]-data_list[1,1]

"""
#%% plot that shit

vmin = np.nanmin(data_calc[0][2])
vmax = np.nanmax(data_calc[0][2])

plt.figure(dpi=300)
ax = plt.axes(projection=ccrs.AzimuthalEquidistant(central_latitude=-90))
norm = Normalize(vmin=vmin, vmax=vmax)
# sc = ax.pcolormesh(data_calc[0][0],data_calc[0][1],data_calc[0][2],cmap='jet',vmin=vmin,vmax=vmax)
ax.set_extent([250000,2750000,-2200000,-38000],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, rotate_labels = False)
gl.top_labels = False
gl.right_labels = False
gl.ylabel_style = {'rotation': 45, 'rotation_mode': 'anchor'}
gl.xformatter = LongitudeFormatter()
gl.yformatter = LatitudeFormatter()
pr = ax.plot(x_int_s,y_int_s,'m')
ax.contour(data_calc[1][0],data_calc[1][1],data_calc[1][2],levels=[0,1500,3000,4500],colors='k',linewidths=0.5)
# cbar = plt.colorbar(cm.ScalarMappable(norm=norm,cmap='jet'), ax=ax,fraction=0.038, pad=0.04)
# cbar.set_label("Sediment depth in m")
# cbar.formatter.set_powerlimits((-2, 2))
# cbar.formatter.set_useMathText(True)
# cbar.locator = ticker.MaxNLocator(nbins=7)
# cbar.update_ticks()
ax.coastlines(resolution='10m',linewidth=0.5)
plt.show()

"""

#%%___build the world_________________________________

mesh,line_list = gepy.build_world_2d(data_list, profile_coord, area=None,
                                      layers=layers, do_sediments=do_sediments)





#%%___apply correct markers to cells and nodes________

mesh,force = gepy.assign_markers_2d(data_list, mesh,
                         layers=layers, do_sediments=do_sediments, do_hp=do_hp, hp_0=hp_0, tc=tc)
    
# #%% quick world plot to check

# fig, ax = plt.subplots(figsize=(9,9),dpi=300)
# pg.show(mesh,showMesh=True,markers=True,ax=ax)
# # pg.show(mesh,data=T,showMesh=True, label='Temperature in °C', cMap="inferno",nCols=265,ax=ax)
# # # pg.show(mesh,data=gradient[:,1],showMesh=True, ax=ax,fitView=False,vmin=-0.04,vmax=-0.03)
# # # ax.set_aspect('equal','box')
# # ax.set_xlim([0,20000])
# # ax.set_ylim([-1000,-250])
# ax.set_xlabel('x in m')
# ax.set_ylabel('depth in m')

#%% run calc and get shf

T = gepy.calc_temp(mesh, "2D", force, layers,do_sediments = True,tc = tc)

# pg.show(mesh,data=T,showMesh=True, label='Temperature in °C', cMap="inferno",nCols=256)

ghf_2D = gepy.calc_ghf(T, mesh, data_list, "2D", tc)
       


#%% temperature depth anomaly plot
'''
# fig,ax = plt.subplots(figsize=(9,9),dpi=300)
# pg.show(mesh,showMesh=True,markers=True,ax=ax,clipBoundaryMarkers=True)
# ax.set_xlim([40000,60000])
# ax.set_ylim([-2000,2000])
# ax.set_xlabel('x in m')
# ax.set_ylabel('depth in m')

nodes_pos = [mesh.nodes()[i].pos() for i in range(len(mesh.nodes()))]

depth_bins = np.arange(np.min(-profile_lab),np.max(profile_topo)+100,100)
depth_mean = np.zeros(len(depth_bins))
T_arr = np.array(T)
n_pos = np.array(nodes_pos)
for i,d in enumerate(depth_bins[:-1]):
    sel = (n_pos[:,1] >= depth_bins[i]) & (n_pos[:,1] < depth_bins[i+1])
    depth_mean[i] = T_arr[sel].mean()

pred_T = np.interp(n_pos[:,1],depth_bins,depth_mean)
T_anom = T-pred_T

for i, node in enumerate(mesh.nodes()):
    x,y,z = node.pos()
    
    idx_t_0,idx_t_1 = np.argsort(np.abs(dist_prof_t-x))[0],np.argsort(np.abs(dist_prof_t-x))[1]
    
    m_t = (profile_topo[idx_t_0]-profile_topo[idx_t_1])/(dist_prof_t[idx_t_0]-dist_prof_t[idx_t_1])
    y_t = m_t * (x-dist_prof_t[idx_t_1]) + profile_topo[idx_t_1] 
    
    if y >= y_t:
        T_anom[i] = 0


# gradient_full = pg.solver.grad(mesh,T)


# gradient_nodes = pg.solver.grad(mesh, T, nodes_pos)
# gradient_nodes_norm = np.linalg.norm(gradient_nodes,axis=1)

fig,ax = plt.subplots(figsize=(9,4),dpi=300)
pg.show(mesh,data=T_anom,cMin=-25,cMax=25,showMesh=False,showBoundary=True,label='$\Delta$T in °C',orientation='vertical', cMap="RdYlBu_r",nCols=50,ax=ax,shading='flat')#,cMin=0,cMax=10
# pg.show(mesh,data=gradient_nodes_norm,showMesh=True,showBoundary=True,label='idk', cMap="plasma",nCols=50,ax=ax,shading='gouraud')#,cMin=0,cMax=10
# pg.show(mesh,linewidth=.5,data=np.linalg.norm(gradient_full,axis=1),cMin=0.01,cMax=0.04,showMesh=True,showBoundary=True,label='Tempterature gradient in K/m', cMap="plasma",nCols=50,ax=ax)#,cMin=0,cMax=10
# pg.show(mesh,data=gradient_full[:,1],showMesh=True,showBoundary=True,label='idk', cMap="plasma",nCols=50,ax=ax,shading='gouraud')#,cMin=0,cMax=10
# pg.show(mesh,data=T,showMesh=True,showBoundary=True,label='Temperature in °C',cMin=0,cMax=10, cMap="plasma",nCols=50,ax=ax)#,cMin=0,cMax=10
# ax.set_xlim([40000,65000])
ax.set_ylim([-12500,1000])
ax.set_xlabel('Profile length in m')
ax.set_ylabel('Depth in m')
ax.set_aspect('auto')

# from pygimli.viewer.mpl import drawStreams
# fig,ax = plt.subplots(figsize=(9,9),dpi=300)
# pg.show(mesh,data=T,showMesh=True,showBoundary=True,label='Temperature in °C',cMin=0,cMax=50, cMap="plasma",nCols=50,ax=ax,shading='gouraud')#,cMin=0,cMax=10
# drawStreams(ax,mesh,gradient_full*100000000,color='green',quiver=True)
# ax.set_xlim([40000,50000])
# ax.set_ylim([-500,500])
# ax.set_xlabel('x in m')
# ax.set_ylabel('depth in m')
# ax.set_aspect('auto')
'''

#%% Plot data and compare to 1D result

# plt.figure(dpi=300)
# plt.plot(dist_prof_t,gradientTop[:,1])   
# plt.plot(dist_prof_t,shf_2D)      
# plt.plot(dist_prof_A,profile_A)   

profile_moho_t = profile_moho_for_topo - profile_sed_for_topo + profile_topo
profile_lab_t = profile_lab_for_topo - profile_moho_for_topo

A = profile_A_for_topo
# A = 0.000001

shf_1D = ( 1315 - 0 + (A*profile_moho_t*profile_lab_t)/k2 + (1/2*A*profile_moho_t**2)/k1 )/( profile_lab_t/k2 + profile_moho_t/k1 + profile_sed_for_topo/k0 )*1000

plt.figure(figsize=(8,3),dpi=300)
plt.plot(dist_prof_t,shf_2D,label="2D")
# plt.plot(dist_prof_t,shf_2D_interp,label="2D interp")
plt.plot(dist_prof_t,shf_1D,label="1D")  
plt.xlabel('x in m')
plt.ylabel('GHF in mW/m$^2$')
plt.legend(loc='upper right')

# plt.figure(figsize=(8,3),dpi=300)
# plt.plot(dist_prof_t,shf_2D,label="2D")
# plt.plot(collect_cells_unique[:,1],collect_cells_unique[:,3]*-1.5*1000,label="2D interp")
# plt.plot(dist_prof_t,shf_1D,label="1D")  
# plt.xlabel('x in m')
# plt.ylabel('GHF in mW/m$^2$')
# plt.legend()

# plt.figure(dpi=300) 
# plt.plot(dist_prof_t,shf_2D-shf_1D)  

# plt.figure(dpi=300) 
# plt.plot(dist_prof_t,profile_topo,label='topo')  
# plt.plot(dist_prof_t,profile_topo-profile_sed_for_topo,label='sed depth')
# plt.xlabel('x in m')
# plt.ylabel('height in m')
# plt.legend() 


# interface (change acordingly):

# fig, (ax1, ax2) = plt.subplots(2, 1,sharex=True,dpi=300)

# ax1.plot(dist_prof_t,shf_2D,label="2D")
# ax1.plot(dist_prof_t,shf_1D,label="1D")  
# ax1.set_ylabel('GHF in mW/m$^2$')
# ax1.legend(loc='lower right')

# # ax2.plot(dist_prof_t,profile_topo/1000,color='slategrey')
# # ax2.plot(dist_prof_t,(profile_topo-profile_sed_for_topo)/1000,color='gold')
# # ax2.plot(dist_prof_m,-profile_moho/1000,color='darkorange')
# # ax2.plot(dist_prof_l,-profile_lab/1000,color='orangered')
# ax2.plot(dist_prof_t,profile_A_for_topo,color='darkorchid')
# ax2.set_xlabel('x in m')
# # ax2.set_ylabel('height in km')
# # ax2.set_ylabel('thickness in km')
# # ax2.set_ylabel('depth in km')
# ax2.set_ylabel('HP in W/m$^3$')


# all interfaces:

# fig, (ax1, ax2, ax3) = plt.subplots(3, 1,sharex=True,dpi=300)

# ax1.plot(dist_prof_t,(profile_topo-profile_sed_for_topo)/1000,color='yellow',label='Sediment depth')
# ax1.plot(dist_prof_t,profile_topo/1000,label='Topography')
# # ax1.set_xlabel('x in m')
# ax1.set_ylabel('height in km')
# # ax1.colorbar(label='mW/m$^2$')

# ax2.plot(dist_prof_m,-profile_moho/1000,color='orange',label='Moho')
# # ax2.set_xlabel('x in m')
# ax2.set_ylabel('depth in km')

# ax3.plot(dist_prof_l,-profile_lab/1000,color='red',label='LAB')
# ax3.set_xlabel('x in m')
# ax3.set_ylabel('depth in km')

# fig.legend()

#%% save
print(notes)

# still need to add seperate msp
save_dict = {'shf_2D': shf_2D,
             'topcoord': topcoord,
             'T': T_list,
             'shf_1D': shf_1D,
             'profile_topo': profile_topo,
             'dist_prof_t': dist_prof_t,
             'x_int_t': x_int_t,
             'y_int_t': y_int_t,
             'profile_sed': profile_sed_for_topo,
             'profile_A': profile_A,
             'dist_prof_A': dist_prof_A,
             'x_int_A': x_int_A,
             'y_int_A': y_int_A,
             'profile_moho': profile_moho,
             'dist_prof_m': dist_prof_m,
             'x_int_m': x_int_m,
             'y_int_m': y_int_m,
             'profile_lab': profile_lab,
             'dist_prod_l': dist_prof_l,
             'x_int_l': x_int_l,
             'y_int_l': y_int_l,}

now = datetime.datetime.now()
mesh.exportVTK('outputVTK/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes))
world.exportPLC('outputVTK/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes))
with open('output/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes)+'.pkl','wb+') as file:
    pickle.dump(save_dict,file)

