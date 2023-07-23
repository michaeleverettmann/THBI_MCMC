#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import obspy
import scipy.io as sio
import pygmt
import shapefile
import matplotlib.path as mplPath

# Define region to plot
region = [-88, -70, 26, 46]

# Paths
fstas = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/results_compiled/compiled_results_standard.mat' # sta data
finventory = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/DatabaseOBSPY/AttenBody_201910OBSPY/misc/inventory_june_all.xml' # inventory file, if wanting extra stations.
ffigpath = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/map_view/figure_out/map_fig.pdf' # output figure.
fsect_locs = '../many_stas/xsect_positions.mat'# cross-section locations.

fvs_averages = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/figures/plot_region_averages/vs_percentiles.eps'
faverages_locations = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/figures/plot_region_averages/region_locations.eps'
f_averages_borders = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/figures/plot_region_averages/choose_split_regions/tect_regions.mat'

# Colors
color_front    = '255/255/255' # Grenville Front, Appalachian front.
color_thrust   = '125/125/125'
color_mag      = '190/000/000'
color_grav     = '000/150/000'
color_ha       = '250/020/020' # '255/000/010' # Harrisonburg Anomaly
color_rift     = '255/191/0'

# The regional velocity profile colors
reg_keys = ['app', 'cra', 'gre']
poly_colors = ['black', 'blue', 'yellow']
clr_reg = ['14/207/18', '189/114/229', '201/124/24'] # Taken straight from my Matlab script. Can be very slightly different in map than velocity profiles, to enhance contrast.

# Any other plot parameters.
linewidth = "2.0"; # Base linewidth.

# Modify for matlab here.
# color_front    = [255, 255, 255]./255; % Grenville Front, Appalachian front.
# color_thrust   = [125, 125, 125]./255;
# color_mag      = [190, 000, 000]./255;
# color_grav     = [000, 150, 000]./255;
# color_ha       = [255, 000, 010]./255; % Harrisonburg Anomaly
# color_rift     = [255, 111, 000]./255;



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

app_bord2 = np.array([[-73.9819443,  40.4970924], # Appalachian border, second part.
       [-74.1098586,  40.996484 ],
       [-73.8461226,  42.4396742],
       [-73.3406284,  43.8503745],
       [-73.0768923,  44.8247083],
       [-72.4395302,  45.6140374]])

gre_bord = np.array([[-82.8790832,  43.9928145], # Grenville obrder.
       [-83.5164453,  41.8040781],
       [-83.6483134,  40.8969058],
       [-84.0878735,  39.6733704],
       [-84.3516096,  38.9764925],
       [-85.714246 ,  36.6155276],
       [-87.3406184,  34.9940038],
       [-88.1098487,  34.3978449]])

#%% load lons and lats of ENAM OBSs
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

#%% Load topography grid in target area
grid = pygmt.datasets.load_earth_relief(resolution="01m", region=region)
#%%

# Initialize figure
fig = pygmt.Figure()

# Plot grid
fig.basemap(region=region, projection="M12c", frame=["f"])
fig.grdimage(grid=grid, cmap="oleron")
fig.coast(shorelines=True, frame=True,     lakes="lightblue",
    borders=["1/thick,black", "2/thin,black", "3/thin,black"],
    )

# Scatter stations

# Load polygons for different regions where we plotted results.

av_locs = sio.loadmat(f_averages_borders)

which_poly = np.zeros((len(lon_sta)))-1
for ipolystr, polystr in enumerate(reg_keys):
    poly = av_locs[polystr]

    p = mplPath.Path(poly) # https://stackoverflow.com/questions/31542843/inpolygon-examples-of-matplotlib-path-path-contains-points-method
    in_this_poly = p.contains_points(np.array([lon_sta, lat_sta]).T)
    which_poly[in_this_poly] = ipolystr

    fig.plot(x=lon_sta[in_this_poly], y=lat_sta[in_this_poly],
             style="t0.28c", fill=clr_reg[ipolystr], pen="black")
fig.plot(x=lon_sta[which_poly==-1], y=lat_sta[which_poly==-1],
         style="t0.28c", fill='160/160/160', pen="black")

# fig.plot(x=lon_sta, y=lat_sta, style="t0.15c", fill="black", pen="black")
# fig.plot(x=lon_obs, y=lat_obs, style="t0.25c", fill="yellow", pen="black")

# # Plot gravity and magnetics
# for ipdat, (pdat, color) in enumerate(zip(['grav_PGA', 'mag_ECMA', 'mag_BSMA', 'mag_BMA'], [color_grav, color_mag, color_mag, color_mag])):
#     pfile = 'potential_field_contours/'+pdat+'.mat'
#     print(pfile)
#     pdatx, pdaty = sio.loadmat(pfile)['theline'][0][0]
#     pdatx = pdatx.ravel()
#     pdaty = pdaty.ravel()
#     fig.plot(x=pdatx, y=pdaty, close=True, pen="1.5p,"+color),

# Plot South Georgia rift
# for i in [1]:
#     sgr_poly = pd.read_csv('SGR_polygon/SGR_poly_pt' + str(i) + '.txt', header = None,  delim_whitespace=True)
#     x, y = (sgr_poly[0].values, sgr_poly[1].values)
#     fig.plot(x=x, y=y, close=False, pen=linewidth+"p,"+color_rift)
#     fig.plot(x=x[-1:0:-1], y=y[-1:0:-1], close=False, pen=linewidth+"p,"+color_rift) # Glitch! pygmt only seems to plot a few points at a time. So plot the reversed lines also to make sure everything plots...

# Plot Grenville and Appalachian fronts
for ibord, bord in enumerate([app_bord, app_bord2]): # gre_bord
    x = bord[:,0]
    y = bord[:,1] 
    fig.plot(x=x, y=y, pen=linewidth+"p,"+color_front)

# Plot Harrisonburg anomaly circle
# fig.plot(x=[-80], y=[38.3], size=[200], style='E-', pen=linewidth+"p,"+color_ha)

#%% Whitmeyer and Karlstrom stuff
def shp_to_xy(shape):
    xall = np.zeros(len(shape))
    yall = np.zeros(len(shape))
    for ipt in range(len(shape)):
        x, y = shape[ipt]
        xall[ipt] = x
        yall[ipt] = y
    return xall, yall

def plot_many_shapes(shapes, color=color_rift):
    for ishape, shape in enumerate(shapes):
        x, y = shp_to_xy(shape.points)
        fig.plot(x=x, y=y, pen=linewidth+"p," + color)

# # Mid-continent rift
# shapes = shapefile.Reader('whitmeyer_karlstrom/layers/MCR.shx').shapes()
# plot_many_shapes(shapes, color_rift)
#
# # Reelfoot rift
shapes = shapefile.Reader('whitmeyer_karlstrom/layers/Reelfoot.shx').shapes()
shapes = [shapes[1]] # Just plot RT, not the other one
plot_many_shapes(shapes, color_rift)

# Grenville front
shape = shapefile.Reader('whitmeyer_karlstrom/layers/grv_frt.shx').shapes()[0].points
x, y = shp_to_xy(shape)
fig.plot(x=x, y=y, pen=linewidth+"p,"+color_front)

# shapes = shapefile.Reader('whitmeyer_karlstrom/layers/something_province.shx').shapes()#[0].points # This is misnamed, it is thrusts. But I can't change file names or I get a glitch.
# plot_many_shapes(shapes, color_thrust)


#%% Label features
fig.text(x=-83.2,  y=40,    text='GF',  font="14p,"+color_front, angle=80)
fig.text(x=-83,    y=36.5,  text='AF',  font="14p,"+color_front, angle=35)
# fig.text(x=-84.1,  y=31.6,  text='SGR', font="12p,"+color_rift ,         )
fig.text(x=-83.5,  y=37.75, text='RT' , font="14p,"+color_rift ,         )
# fig.text(x=-87.25, y=36.75, text='RR' , font="12p,"+color_rift ,         )
fig.text(x=-84.9,    y=35,  text='SAA',  font="14p,"+color_ha   ,         )
fig.text(x=-80,    y=38.3,  text='CAA',  font="14p,"+color_ha   ,         )
fig.text(x=-72.4,  y=43.4,  text='NAA',  font="14p,"+color_ha   ,         )
# fig.text(x=-85.25, y=42.5,  text='MCR', font="12p,"+color_rift ,         )
# fig.text(x=-71, y=40.5, text='PGA', font="12p,"+color_grav)
# fig.text(x=-78, y=33.5, text='BMA', font="12p,"+color_mag)
# fig.text(x=-77, y=30.8, text='ECMA', font="12p,"+color_mag)
# fig.text(x=-75, y=30, text='BSMA', font="12p,"+color_mag)
# fig.text(x=-79.0298, y=35.8908, text='US.CEH', font="8p,255/255/255")

#%% Add inset
# with fig.inset(position="jBR+w3.0c+o0.1/0.5c", margin=0): # Good if not putting velocity average figure here
with fig.inset(position="jBL+w2.5c+o0.2/0.2c", margin=0.1): # Move to left side of figure if showing velocity average figure.

    # State and coast lines.
    fig.coast(
        region="g",
        projection="G-77/36/?",
        land="gray",
        water="white",
        borders=[1, 2],
        shorelines="1/thin",
        frame="afg",
    )

    # Plot a rectangle around our area.
    rectangle = [[region[0], region[2], region[1], region[3]]]
    fig.plot(data=rectangle, style="r+s", pen=linewidth+"p,blue")

fig.image(fvs_averages, position="jBR+w3.8c+o0.15/0.15c", box="+pblack+g255/255/255")
# fig.image(fvs_averages, position="jBR+w4.0c+o0.15/0.15c")
# fig.image(faverages_locations, position="jBR+w2c+o0.008/7.13c", box="+pblack+g0/0/0+c0/0")


#%%
fig.show()

#%%
fig.savefig(ffigpath)