#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# First attempt to switch to pygmt from an old geophgraphic/tectonic figure. 2023.03.21

fstas = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/results_compiled/compiled_results_standard.mat'
finventory = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/DatabaseOBSPY/AttenBody_201910OBSPY/misc/inventory_june_all.xml'
region = [-88, -68, 26, 46] # Define region of interest 
# ffigpath = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/DatabaseOBSPY/AttenBody_201910OBSPY/figures/geog_map_pygmt.pdf';
ffigpath = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/map_view/figure_out/map_fig.pdf'
# Some colors
tcol = '0/0/255'
pcol = '190/0/0'
gcol = '0/150/0'
hacol = '255/0/150'
sgrcolor = '255/255/255'


import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np
import pandas as pd
import obspy
import matplotlib.patheffects as PathEffects
from matplotlib import patches 
import sys
# import inout
import scipy.io as sio
import pygmt
import shapefile

# from myinterp.gaus import gaus_smooth_irreg
# from mycolors.customColorMap import interpolateColormap

#%% Define some geographic features/shapes
app_bord = np.array([    [-73.9819443,  40.4970924], # appalachian border
       [-74.619146 ,  40.7306085],
       [-75.520073 ,  40.9467137],
       [-76.1132732,  41.1455697],
       [-76.7944389,  41.2695495],
       [-77.1789631,  41.2365112],
       [-77.6294266,  41.0959121],
       [-77.9920049,  40.8719878],
       [-78.2117582,  40.7056279],
       [-78.4383114,  40.3967643],
       [-79.1416076,  39.5548831],
       [-79.8888598,  38.4965935],
       [-80.2331436,  38.1259146],
       [-80.7605221,  37.7272803],
       [-81.2438831,  37.5445773],
       [-81.9140786,  37.3439591],
       [-82.4853569,  37.1603165],
       [-83.2434644,  36.8268747],
       [-84.3284171,  36.1733569],
       [-84.7700317,  35.880149 ],
       [-85.1106908,  35.4427709],
       [-85.5944904,  34.7777158],
       [-86.0889287,  34.1799976],
       [-87.330563 ,  32.9349287]])

app_bord2 = np.array([[-73.9819443,  40.4970924],
                      [-74.1098586,  40.996484 ],
       [-73.8461226,  42.4396742],
       [-73.3406284,  43.8503745],
       [-73.0768923,  44.8247083],
       [-72.4395302,  45.6140374]])

gre_bord = np.array([[-82.8790832,  43.9928145],
       [-83.5164453,  41.8040781],
       [-83.6483134,  40.8969058],
       [-84.0878735,  39.6733704],
       [-84.3516096,  38.9764925],
       [-85.714246 ,  36.6155276],
       [-87.3406184,  34.9940038],
       [-88.1098487,  34.3978449]])

#%% Get lons and lats of ENAM OBSs
inventory = obspy.read_inventory(finventory)
inventory = inventory.select(network='YO')[0]
lon_obs = []
lat_obs = []
for iobs in range(len(inventory)): 
    sta = inventory[iobs]
    lat_obs.append(sta.latitude)
    lon_obs.append(sta.longitude)

#%% Get lons at lats for MCMC stations
stamat = sio.loadmat(fstas)
stas = stamat['mdls']
lon_sta = stas[0]['lon'][0].astype('float64').ravel()
lat_sta = stas[0]['lat'][0].astype('float64').ravel()


#%% Load sample grid (3 arc-minutes global relief) in target area
grid = pygmt.datasets.load_earth_relief(resolution="01m", region=region)
#%%
#  Plot grid
fig = pygmt.Figure()
fig.basemap(region=region, projection="M12c", frame=["f"]) # frame=["f", "+tclipped grid"]
fig.grdimage(grid=grid, cmap="oleron")
fig.coast(shorelines=True, frame=True,
    borders=["1/thick,black", "2/thin,black", "3/thin,black"],
    )

# Scatter stations
fig.plot(x=lon_sta, y=lat_sta, style="t0.15c", fill="black", pen="black")
fig.plot(x=lon_obs, y=lat_obs, style="t0.25c", fill="yellow", pen="black")

# Plot gravity and magnetics
for ipdat, (pdat, color) in enumerate(zip(['grav_PGA', 'mag_ECMA', 'mag_BSMA', 'mag_BMA'], [gcol, pcol, pcol, pcol])): 
    pfile = 'potential_field_contours/'+pdat+'.mat'
    print(pfile)
    pdatx, pdaty = sio.loadmat(pfile)['theline'][0][0]
    pdatx = pdatx.ravel()
    pdaty = pdaty.ravel()
    fig.plot(x=pdatx, y=pdaty, close=True, pen="1.5p,"+color)#, 

# Plot South Georgia rift
for i in [1]: 
    sgr_poly = pd.read_csv('SGR_polygon/SGR_poly_pt' + str(i) + '.txt', header = None,  delim_whitespace=True)
    x, y = (sgr_poly[0].values, sgr_poly[1].values)
    fig.plot(x=x, y=y, close=False, pen="1.5p,"+sgrcolor)
    fig.plot(x=x[-1:0:-1], y=y[-1:0:-1], close=False, pen="1.5p,"+sgrcolor) # Glitch! pygmt only seems to plot a few points at a time. So plot the reversed lines also to make sure everything plots...

# Plot Grenville and Appalachian fronts
for ibord, bord in enumerate([app_bord, app_bord2]): # gre_bord
    x = bord[:,0]
    y = bord[:,1] 
    fig.plot(x=x, y=y, pen="1.5p,"+sgrcolor)

# Plot Harrisonburg anomaly circle
fig.plot(x=[-80], y=[38.3], size=[200], style='E-', pen="1.5p,"+hacol)

#%% Whitmeyer and Karlstrom stuff
def shp_to_xy(shape):
    xall = np.zeros(len(shape))
    yall = np.zeros(len(shape))
    for ipt in range(len(shape)):
        x, y = shape[ipt]
        xall[ipt] = x
        yall[ipt] = y
    return xall, yall

def plot_many_shapes(shapes, color=sgrcolor):
    for ishape, shape in enumerate(shapes):
        x, y = shp_to_xy(shape.points)
        fig.plot(x=x, y=y, pen="1.5p," + color)

# Mid-continent rift
shapes = shapefile.Reader('whitmeyer_karlstrom/layers/MCR.shx').shapes()
plot_many_shapes(shapes, hacol)

# Reelfoot rift
shapes = shapefile.Reader('whitmeyer_karlstrom/layers/Reelfoot.shx').shapes()
plot_many_shapes(shapes, hacol)

# Grenville front
shape = shapefile.Reader('whitmeyer_karlstrom/layers/grv_frt.shx').shapes()[0].points
x, y = shp_to_xy(shape)
fig.plot(x=x, y=y, pen="1.5p,"+sgrcolor)

shapes = shapefile.Reader('whitmeyer_karlstrom/layers/something_province.shx').shapes()#[0].points # This is misnamed, it is thrusts. But I can't change file names or I get a glitch.
plot_many_shapes(shapes, "155/155/155")


#%% Label features
fig.text(x=-83.3, y=40  , text='Grenville front' , font="12p,white", angle=70)
fig.text(x=-83, y=36  , text='Appalachian front', font="12p,white", angle=40)
fig.text(x=-84.5, y=32.8, text='SGR', font="12p,white")
fig.text(x=-71, y=40.5, text='PGA', font="12p,"+gcol)
fig.text(x=-78, y=33.5, text='BMA', font="12p,"+pcol)
fig.text(x=-77, y=30.8, text='ECMA', font="12p,"+pcol)
fig.text(x=-75, y=30, text='BSMA', font="12p,"+pcol)
fig.text(x=-80, y=38.3, text='HA', font="12p,"+hacol)
fig.text(x=-79.0298, y=35.8908, text='US.CEH', font="8p,255/255/255")

# Add inset
with fig.inset(position="jBR+w3.0c+o0.1/0.5c", margin=0): 
    fig.coast(
        region="g",
        projection="G-77/36/?",
        land="gray",
        water="white",
        borders=[1, 2],
        shorelines="1/thin",
        frame="afg",
    )

fig.show()
fig.savefig(ffigpath)