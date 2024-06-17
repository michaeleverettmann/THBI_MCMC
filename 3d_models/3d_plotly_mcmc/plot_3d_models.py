# BRB 2023/04/17 Modified from Brunsvik 2021 GGG
# I deleted many parts, but they can be easily gotten back.
#%%
import numpy as np 
import scipy.io as sio
from scipy.interpolate import griddata
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import plotly.graph_objects as go
import plotly.express as px
import shapefile
import pickle
from copy import deepcopy
import scipy.interpolate as interp
import pygmt
from scipy.spatial import ConvexHull, convex_hull_plot_2d


modelPath = '../surface_colated_b1_V7.mat'
stationPath = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/results_compiled/compiled_results_standard.mat'
ymin = -1300; ymax = 1900;
xmin = -1400; xmax = 1100;
zmin = -350; zmax = -50;
scount = 4
vmin_plt = 4.45
vmax_plt = 4.75
prtb = 3.4
opacity = 0.4 # Lower is more see through
dist_sta_plt = 100 # Make a contour showing where a region is further than this from any station

#%% Load station info
msta = sio.loadmat(stationPath)
slat = msta['mdls']['lat'][0][0].flatten().astype('float64')
slon = msta['mdls']['lon'][0][0].flatten().astype('float64')

#%% Load MCMC output model
mat = sio.loadmat(modelPath)

lat3d      =  mat['lat3d'     ].astype('float64')
latgrid    =  mat['latgrid'   ].astype('float64')
lon3d      =  mat['lon3d'     ].astype('float64')
longrid    =  mat['longrid'   ].astype('float64')
mgrid3d    =  mat['mgrid3d'   ].astype('float64')
xgrid      =  mat['xgrid'     ].astype('float64')
ygrid      =  mat['ygrid'     ].astype('float64')
depths     = -mat['depths'    ].astype('float64')
z3d        = -mat['z3d'       ].astype('float64')
zmoh_surf  =  mat['zmoh_surf' ].astype('float64')
zsed_surf  =  mat['zsed_surf' ].astype('float64')

#%%
mshape = mgrid3d.shape # nx x ny x nz

# These are the things we will plot. Get them shaped nicely here.
xmdl = np.zeros(mshape) + np.nan
ymdl = np.zeros(mshape) + np.nan
zmdl = np.zeros(mshape) + np.nan
vmdl = np.zeros(mshape) + np.nan

# Repeat x and y for each depth.
for iz in range(mshape[2]):
    xmdl[:,:,iz] = xgrid
    ymdl[:,:,iz] = ygrid

zmdl = z3d
vmdl = mgrid3d

# Have to flatten for plotly
xmdl = xmdl.flatten()
ymdl = ymdl.flatten()
zmdl = zmdl.flatten()
vmdl = vmdl.flatten()

# Units. brb2023/04/17 Temporary, convert roughly to km
lat_to_km = 6371*2*np.pi/360 # about 111, did the math for... fun?
dykm = (latgrid.max() - latgrid.min()) * lat_to_km
dymap = (ygrid.max() - ygrid.min())
ymdl = ymdl * dykm / dymap
xmdl = xmdl * dykm / dymap

#%% Make object to convert from lon lat to our x y system
xterp = xgrid.flatten() * dykm / dymap# Our original xgrid. Flatten it and use it for interpolation. Use dykm/dymap to convert to km.
yterp = ygrid.flatten() * dykm / dymap
lonterp = longrid.flatten() # Our original non-linear longitude grid. Flatten and use it for interpolation.
latterp = latgrid.flatten()
lontobj = interp.CloughTocher2DInterpolator(np.array([xterp, yterp]).T, lonterp) # Interpolates from original x and y to lon
lattobj = interp.CloughTocher2DInterpolator(np.array([xterp, yterp]).T, latterp)
xtobj = interp.CloughTocher2DInterpolator(np.array([lonterp, latterp]).T, xterp) # Interpolates from original lon and lat to x
ytobj = interp.CloughTocher2DInterpolator(np.array([lonterp, latterp]).T, yterp) # Is piecewise cubic smooth

# Only keep within desired bounds.
keep = (xmdl < xmax) * (xmin < xmdl) * (ymdl < ymax) * (ymin < ymdl) * (zmdl < zmax) * (zmin < zmdl)
xmdl = xmdl[keep]; ymdl = ymdl[keep]; zmdl = zmdl[keep]; vmdl = vmdl[keep]

#%% Convert station positions to new coordinate system
xsta = xtobj(slon, slat)
ysta = ytobj(slon, slat)

min_sta_dist = 100 # km
sta_dist = np.zeros(xmdl.shape)
for ista, (_xmdl, _ymdl) in enumerate(zip(xmdl, ymdl)):
    sta_dist[ista] = np.sqrt( ((xsta-_xmdl)**2 + (ysta-_ymdl)**2).min() )

tofar = sta_dist > dist_sta_plt
closesta = ~tofar

#%%
distobj = interp.CloughTocher2DInterpolator(np.array([xmdl, ymdl]).T, sta_dist)
xdst = np.linspace(xmdl.min(), xmdl.max(), 250)
ydst = np.linspace(ymdl.min(), ymdl.max(), 250)
xdst, ydst = np.meshgrid(xdst, ydst)
dist_grid = distobj(xdst, ydst)
cntr = plt.contour(xdst, ydst, dist_grid, levels = [dist_sta_plt,9999999]) # A contour at dist_sta_plt. To tell matplotlib we aren't specifyin g anumber of levels, we need to give a list of increasing values. do 999999 since nothing is that far.
p1 = cntr.collections[0].get_paths()

#%%
dataDist = []
for icnt in range(len(p1)):
    dist_cntr = p1[icnt].vertices
    x_dist_cntr = dist_cntr[:,0]
    y_dist_cntr = dist_cntr[:,1]

    dataDist += [go.Scatter3d(
                x = x_dist_cntr,
                y = y_dist_cntr,
                z = (zmax-5) * np.ones(x_dist_cntr.shape),
                showlegend = False,
                mode = 'lines',
                name = '',
                line=dict(
                    color='orange',
                    width=6,
                )
            )]

#%% Remove mean velocity for each layer. Need to do this after removing non-pertinent edges of model
z_all = np.sort(np.unique(zmdl))
dvmdl = vmdl.copy()
for iz, z in enumerate(z_all):
    zb = zmdl == z # z boolean

    # Remove the mean
    # When calculating mean, only utilize nodes close to stations.
    dvmdl[zb] = (dvmdl[zb] - np.mean(dvmdl[zb*closesta])) / np.mean(dvmdl[zb*closesta]) * 100

#%% Get state boundaries information.
thisfile = 'data/other_states/cb_2020_us_state_20m.shp'
shape = shapefile.Reader(thisfile) # Get us government file state border info.
shape = shape.shapes()

latstate = [] # extract lon and lat points from PITA shapefile.
lonstate = []
xstate = []
ystate = []
for ishape, _shape in enumerate(shape):
    points = _shape.points
    parts = _shape.parts
    parts = np.append(parts, len(points))
    parttodo = len(parts)-1 if len(parts) > 2 else 1
    for ipart in range(parttodo): # Have to loop through disconnected parts of states.
        subpoints = points[parts[ipart]:parts[ipart+1]] # Should be points[parts[ipart]:parts[ipart+1]+1] but this gives glitchy lines... bad dataset?
        lats = []
        lons = []
        for ipt, pt in enumerate(subpoints):
            lats.append(pt[1])
            lons.append(pt[0])
        lats = np.array(lats)
        lons = np.array(lons)
        _x = xtobj(lons, lats) # Interpolate to new coordinate system.
        _y = ytobj(lons, lats)
        keep = np.ones(_x.shape, dtype = 'bool')
        keep[np.isnan(_x) * np.isnan(_y)] = False
        _x = _x[keep]
        _y = _y[keep]
        lats = lats[keep]
        lons = lons[keep]
        latstate.append(lats)
        lonstate.append(lons)
        xstate.append( _x )
        ystate.append( _y )

dataStates = []
istate = 0
lineDownsample = 1
for istate in range(len(xstate)):
    _xstate = xstate[istate][::lineDownsample]
    _ystate = ystate[istate][::lineDownsample]
    _zstate = zmdl.max() * np.ones(_ystate.shape)
    dataStates.append(
        go.Scatter3d(
            x = _xstate,
            y = _ystate,
            z = _zstate,
            showlegend = False,
            mode = 'lines',
            name = '',
            line=dict(
                color='black',
                width=2,
            )
        )
    )

#%% Function to update layout to be consistent. 
def update_layout(fig):
    xrng = [xmdl.min(), xmdl.max()] # Lazy use of globals
    yrng = [ymdl.min(), ymdl.max()]
    zrng = [zmdl.min(), zmdl.max()] # TODO extend vertical limit.
    scaler = 1/np.diff(xrng)[0]
    fig.update_layout(
        template = 'plotly_white',
        scene = dict(
                        xaxis_title='East (km)',
                        yaxis_title='North (km)',
                        zaxis_title='Depth (km)',
                        xaxis = dict(nticks=4, range=xrng,),
                        yaxis = dict(nticks=4, range=yrng,),
                        zaxis = dict(nticks=4, range=zrng,),
                        aspectmode='manual',
                        aspectratio={'x':np.diff(xrng)[0]*scaler,
                                    'y':np.diff(yrng)[0]*scaler,
                                    'z':np.diff(zrng)[0]*scaler},
        ),
    )
    return fig

cbar = dict(x=-.02, title = '%dVs')

dataVel = [
        go.Isosurface(
        value = vmdl, isomin = vmin_plt, isomax = vmax_plt, colorscale = 'rdbu', reversescale = False, cmin = vmin_plt, cmax = vmax_plt,
        colorbar = cbar,
        x = xmdl, y = ymdl, z = zmdl,
        opacity=opacity, surface_count=2,
        caps=dict(x_show=False, y_show=False, z_show=False),
                                    )]

d_v_step = (vmax_plt-vmin_plt)/3; # If there are four contours, add 1/3 diff of max and min for step size. [0/3, 1/3, 2/3, 3/3]*(max-min)+min is all contour values.
dataVel += [
        go.Isosurface(
        value = vmdl, isomin = vmin_plt+d_v_step, isomax = vmax_plt-d_v_step, colorscale = 'rdbu', reversescale = False, cmin = vmin_plt, cmax = vmax_plt,
        colorbar = cbar,
        x = xmdl, y = ymdl, z = zmdl,
        opacity=.15, surface_count=2,
        caps=dict(x_show=False, y_show=False, z_show=False),
                                    )]

d_prtb_step = (2*prtb)/3; # If there are four contours, add 1/3 diff of max and min for step size. [0/3, 1/3, 2/3, 3/3]*(max-min)+min is all contour values.
# # Fine tuned version.
dataVelR = [ # Relative
        go.Isosurface(
        value = dvmdl, isomin = -prtb, isomax = prtb, colorscale = 'rdbu', reversescale = False, cmin = -prtb, cmax = prtb,
        colorbar = cbar,
        x = xmdl, y = ymdl, z = zmdl,
        opacity=opacity, surface_count=2,
        caps=dict(x_show=False, y_show=False, z_show=False),
                                    )]
dataVelR +=[ # Relative
        go.Isosurface(
        value = dvmdl, isomin = -prtb+d_prtb_step, isomax = prtb-d_prtb_step, colorscale = 'rdbu', reversescale = False, cmin = -prtb, cmax = prtb,
        colorbar = cbar,
        x = xmdl, y = ymdl, z = zmdl,
        opacity=.15, surface_count=2,
        caps=dict(x_show=False, y_show=False, z_show=False),
                                    )]

# Temporary
fig = update_layout( go.Figure(data=dataVel+dataStates+dataDist) )
layout = fig.layout
layout['scene']['camera']['eye'] = {'x':0, 'y':-1, 'z':.65}
fig.update({'layout':layout})
fig.show()
fig.write_html('figs/3d_mcmc_absolute.html')

fig = update_layout( go.Figure(data=dataVelR+dataStates+dataDist) )
layout = fig.layout
layout['scene']['camera']['eye'] = {'x':0, 'y':-1, 'z':.65}
fig.update({'layout':layout})
fig.show()
fig.write_html('figs/3d_mcmc_relative.html')
