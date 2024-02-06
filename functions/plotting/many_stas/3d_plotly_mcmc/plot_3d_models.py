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
# import matlab.engine
# import matlab
import shapefile
import pickle
from copy import deepcopy
# from mysave import inout
import scipy.interpolate as interp
import pygmt
from scipy.spatial import ConvexHull, convex_hull_plot_2d


modelPath = '../surface_colated_b1_V7.mat'
stationPath = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/results_compiled/compiled_results_standard.mat'
# v_style = 'relative'
# ymin = -611; ymax = 1300;
# xmin = -1000; xmax = 250;
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
#
# # Make a grid of evenly spaces latitude and longitude coordinates for plotly isosurface
# nlon = 200; nlat = 201;
# lonnew = np.linspace(lonterp.min(), lonterp.max(), nlon)
# latnew = np.linspace(latterp.min(), latterp.max(), nlat)
# lonnew, latnew = np.meshgrid(lonnew, latnew, indexing = 'ij')
#
# xnew = xtobj(lonnew, latnew) # Interpolate from regular lonlat grid to non-linear x and y coordinates.
# ynew = ytobj(lonnew, latnew)
#
# # nz = z3d.shape[2]
# zunique = np.sort(np.unique(z3d))
# nz = len(zunique)
# vmdl = np.zeros((nlon, nlat, nz))+np.nan # Interpolated m values, on regular lat lon grid.
#
# for iz in range(nz):
#     mtobj = interp.RegularGridInterpolator((xgrid[:,0], ygrid[0,:]), mgrid3d[:,:,iz], bounds_error=False,
#                                            method = 'linear') # Object to interpolate given x and y coordinates to m
#     vmdl[:,:,iz] = mtobj(np.array([xnew.flatten(), ynew.flatten()]).T).reshape(vmdl[:,:,iz].shape) # Interpolate from new, nonlinear x and y coordinates based on old linear x and y coordinates. Cooresponds to linear lat lon grid.
#     # plt.scatter(xnew, ynew, c = mgrid3dterp[:,:,iz])
#     # plt.show()
#
# # Get lon, lat, z in same shape and gridding as vmdl
# lonmdl = np.zeros((nlon, nlat, nz))
# latmdl = np.zeros((nlon, nlat, nz))
# zmdl   = np.zeros((nlon, nlat, nz))
# for iz in range(nz):
#     lonmdl[:,:,iz] = lonnew
#     latmdl[:,:,iz] = latnew
#     zmdl  [:,:,iz] = zunique[iz]
#
# # Flatten grid, for plotly.
# lonmdl = lonmdl.flatten()
# latmdl = latmdl.flatten()
# zmdl   = zmdl  .flatten()
# vmdl   = vmdl  .flatten()
#
# keep = (xmdl < xmax) * (xmin < xmdl) * (ymdl < ymax) * (ymin < ymdl) * (zmdl < zmax) * (zmin < zmdl)
# xmdl = xmdl[keep]; ymdl = ymdl[keep]; zmdl = zmdl[keep]; vmdl = vmdl[keep]

# Only keep within desired bounds.
keep = (xmdl < xmax) * (xmin < xmdl) * (ymdl < ymax) * (ymin < ymdl) * (zmdl < zmax) * (zmin < zmdl)
xmdl = xmdl[keep]; ymdl = ymdl[keep]; zmdl = zmdl[keep]; vmdl = vmdl[keep]

#%% Convert station positions to new coordinate system
xsta = xtobj(slon, slat)
ysta = ytobj(slon, slat)

#%%
# min_sta_dist = 100 # km
# sta_dist = np.zeros(xmdl.shape)
# for ista, (_xsta, _ysta) in enumerate(zip(xsta, ysta)):
#     sta_dist[ista] = np.sqrt( ((_xsta-xmdl)**2 + (_ysta-ymdl)**2).min() )

min_sta_dist = 100 # km
sta_dist = np.zeros(xmdl.shape)
for ista, (_xmdl, _ymdl) in enumerate(zip(xmdl, ymdl)):
    sta_dist[ista] = np.sqrt( ((xsta-_xmdl)**2 + (ysta-_ymdl)**2).min() )

tofar = sta_dist > dist_sta_plt
closesta = ~tofar

# #%% Convex hull so we can plot a polygon showing where we have stations.
# hull = ConvexHull(  np.array([xmdl[closesta], ymdl[closesta]]).T  )
# # xhull = hull.points[:,0]
# # yhull = hull.points[:,1]
# xhull = xmdl[closesta][hull.vertices]
# yhull = ymdl[closesta][hull.vertices]
#
# dataHull = [go.Scatter3d(
#             x = xhull,
#             y = yhull,
#             z = (zmax-5) * np.ones(xhull.shape),
#             showlegend = False,
#             mode = 'lines',
#             name = '',
#             line=dict(
#                 color='black',
#                 width=4,
#             )
#         )]

# #%% Plot distance to stations
# dataDist = [go.Surface(x=xmdl, y = ymdl, z = sta_dist,
#             showscale=False,
#             opacity = 0,
#             contours = {"z": {"show": True, "start": 150, "end": 150, "size": 5,
#                               "color":"black"}},
#             # cmin = -30, cmax = 2,
#             name = 'Distance', # visible = True,
#              )
#             ]

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

#%% Set to constant values outside regions of resolution
# for iz, z in enumerate(z_all):
#     zb = zmdl == z # z boolean

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


# # Convert to current reference frame. Very time consuming, something about calling the matlab functions within project_to_xy.
# states = {'statex':[], 'statey':[]}
# for iline, (latplt, lonplt) in enumerate(zip(lats, lons)):
#     print('State number {:3.0f}'.format(iline))
#     x, y = project_to_xy(par, lonplt, latplt)
#     states['statex'].append(x)
#     states['statey'].append(y)



dataStates = []
istate = 0
lineDownsample = 1
for istate in range(len(xstate)):
    # figa = px.line_3d(x = states['statex'][istate], y = states['statey'][istate], z = np.zeros(states['statey'][istate].shape) )
    _xstate = xstate[istate][::lineDownsample]
    _ystate = ystate[istate][::lineDownsample]
    # _zstate = np.zeros(ystate[istate].shape)[::lineDownsample] + maxdepth # Using topography
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

#%%
# plt.figure()
# for ista in range(len(lons)): 
#     # plt.plot(states['statex'][ista], states['statey'][ista]) 
#     plt.plot(lons[ista], lats[ista])
# plt.show()

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
                        #  aspectratio = {}
        ),
        # showlegend = True,
    )
    return fig
# #%% Moho
# if False: # Convert matlab file into something more manageable. A bit slow, so save it and reload it unless the file has changed. 
#     mohof = eng.load('../../get_crustal_values/crust_depth_3d.mat')['dat_for_colton']; # Moho depth file. 
#     shape = np.array(mohof['lon']).shape 
#     mohof['lon'] = np.array(mohof['lon']).ravel() 
#     mohof['lat'] = np.array(mohof['lat']).ravel() 
#     mohof['h'  ] = np.array(mohof['h'  ]).ravel() 
#     mohox, mohoy = project_to_xy(par, np.array(mohof['lon']), np.array(mohof['lat']))
#     mohox = mohox     .reshape(shape) 
#     mohoy = mohoy     .reshape(shape) 
#     mohoz = mohof['h'].reshape(shape)
#     inout.saveDict('data/moho.pkl', {'mohox':mohox, 'mohoy':mohoy, 'mohoz':mohoz})
# else: 
#     mohof = inout.loadDict('data/moho.pkl')
#     mohox = mohof['mohox']
#     mohoy = mohof['mohoy']
#     mohoz = mohof['mohoz']

#     downsampMoho = 4
#     mohox = mohox[::downsampMoho,::downsampMoho]
#     mohoy = mohoy[::downsampMoho,::downsampMoho]
#     mohoz = mohoz[::downsampMoho,::downsampMoho]

# visibleMoho = False 
# dataMoho = [go.Surface(x=mohox, y = mohoy, z = -mohoz,
#             showscale=False,
#             opacity = .1,
#             contours = {"z": {"show": True, "start": -50, "end": 0, "size": 5, 
#                               "color":"black"}},
#             cmin = -30, cmax = 2, 
#             name = 'Moho', visible = visibleMoho,
#              )
#             ]
# print("setting moho visibility to false for now.")
# fig = update_layout( go.Figure(data=dataMoho) )
# # fig = add_buttons(fig)
# fig.show()


# #%% Topography
# from myearthplot import getTopo
# if False: # Re load and plot topo data 
#     topox, topoy, topoz = getTopo.getTopo(
#         min_lon = -92, max_lon = -67, min_lat = 25, max_lat = 47.5,
#     )
#     shape = topox.shape 
#     topox, topoy = project_to_xy(par, topox.ravel(), topoy.ravel())
#     topox = topox.reshape(shape) 
#     topoy = topoy.reshape(shape)
#     inout.saveDict('data/topo.pkl', dict(topox=topox, topoy=topoy, topoz=topoz))
# else: 
#     downsampTopo = 20
#     topof = inout.loadDict('data/topo.pkl')
#     topox = topof['topox'][::downsampTopo,::downsampTopo]
#     topoy = topof['topoy'][::downsampTopo,::downsampTopo] 
#     topoz = topof['topoz'][::downsampTopo,::downsampTopo]/1000

# dataTopo = [go.Surface(x=topox, y = topoy, z = -topoz,
#             showscale=False,
#             opacity = .1,
#             contours = {"z": {"show": True, "start": -2, "end": 7, "size": .2, 
#                               "color":"black"}},
#             cmin = -30, cmax = 2, 
#             name = 'Topo', visible = True,
#              )
#             ]
# fig = update_layout( go.Figure(data=dataTopo) )
# # fig = add_buttons(fig)
# fig.show()
# #%% Button stuff
# def add_buttons(fig, labels = ['Moho', 'V', 'A'], indices = [[1], [2,3], [4,5]], visStart = [0, 1, 1]): 
#     buttons = []
#     for ilabel, label in enumerate(labels): 
#         # Whether first arg is to make visible or invisible depends on whether the thing starts off as visible or not. 
#         buttons += [
#             dict(method='restyle',
#                         label=label,
#                         visible=True,
#                         args=[{'visible': not visStart[ilabel]}, indices[ilabel] ],
#                         args2 = [{'visible':visStart[ilabel]}, indices[ilabel] ],
#                         )
#         ]

#     fig.update_layout(
#         updatemenus=[
#             dict(
#                 buttons = buttons, 
#                 type = 'buttons',
#             )
#         ]
#     )

#     return fig 
# %%

#%% Plot of everything

# dataVelo1 = [
#         go.Isosurface(
#         value = vv, isomin = -1.2, isomax = 1.2, colorscale = 'jet', reversescale = True, cmin = -5, cmax = 5, 
#         colorbar = dict(x=-.02, title = '%dVs'),
#         x = mx, y = my, z = -mz, 
#         opacity=0.2, surface_count=2, 
#         caps=dict(x_show=False, y_show=False, z_show=False),
#                                     )]
# dataVelo2 = [
#         go.Isosurface(
#         value = vv, isomin = -2.2, isomax = 2.2, colorscale = 'jet', reversescale = True, cmin = -5, cmax = 5, 
#         showscale=False,
#         x = mx, y = my, z = -mz, 
#         opacity=1, surface_count=2, 
#         caps=dict(x_show=False, y_show=False, z_show=False),
#                                     )]
# visibleAnis = False
# dataAnis1 = [go.Isosurface(
#         value = aa, isomin = -.5, isomax = .5, colorscale = 'rdbu', cmin = -1, cmax = 1, name = 'A1',visible = visibleAnis,
#         colorbar = dict(x=1.02, title = '%Anis'),
#         x = mx, y = my, z = -mz, 
#         opacity=.2, surface_count=2, # showlegend = True, # Note that using show legend puts something in the legend. If it's in the legend, you can click it to turn it off. 
#         caps=dict(x_show=False, y_show=False, z_show=False),
#                                     )] 
# dataAnis2 = [go.Isosurface(
#         value = aa, isomin = -1, isomax = 1, colorscale = 'rdbu', cmin = -1, cmax = 1, name = 'A2',visible = visibleAnis,
#         showscale=False,
#         x = mx, y = my, z = -mz, 
#         opacity=1, surface_count=2, 
#         caps=dict(x_show=False, y_show=False, z_show=False),
#                                     )]

#%%
# TODO change from jet to turbo

# xgrid3d = np.ones(z3d.shape) * xgrid.reshape(xgrid.shape[0], xgrid.shape[1], 1)
# ygrid3d = np.ones(z3d.shape) * ygrid.reshape(ygrid.shape[0], ygrid.shape[1], 1)

# dataVel = [
#         go.Isosurface(
#         value = vmdl, isomin = 4, isomax = 5, colorscale = 'jet', reversescale = True, cmin = 4, cmax = 5,
#         colorbar = dict(x=-.02, title = 'Vs (km/s)'),
#         x = xmdl, y = ymdl, z = zmdl,
#         opacity=0.5, surface_count=8,
#         caps=dict(x_show=False, y_show=False, z_show=False),
#                                     )]


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

# dataVelR = [ # Relative
#         go.Isosurface(
#         value = dvmdl, isomin = -prtb, isomax = prtb, colorscale = 'rdbu', reversescale = False, cmin = -prtb, cmax = prtb,
#         colorbar = dict(x=-.02, title = 'Vs (km/s)'),
#         x = xmdl, y = ymdl, z = zmdl,
#         opacity=opacity, surface_count=scount,
#         caps=dict(x_show=False, y_show=False, z_show=False),
#                                     )]

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

# fig = update_layout( go.Figure(data=dataMoho+dataTopo+dataVelo1+dataVelo2+dataAnis1+dataAnis2+dataStates) )
# fig = update_layout( go.Figure(data=dataVel) )
# fig = update_layout( go.Figure(data=dataVel+dataVelR) )
# fig = add_buttons(fig, labels = ['Moho', 'Topo','Vel', 'Anis'], indices = [[0], [1], [2,3], [4,5]], visStart = [visibleMoho, True, True, visibleAnis])
# layout = fig.layout
# layout['scene']['camera']['eye'] = {'x':0, 'y':-1, 'z':.65}
# fig.update({'layout':layout})
#
# fig.show()
# fig.write_html('figs/3d_mcmc.html')

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
