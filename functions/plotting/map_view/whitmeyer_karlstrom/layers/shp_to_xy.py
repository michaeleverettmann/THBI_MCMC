#%% Convert shape files to xy files. A shapefile with 10 features will produce 10 files.
#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np
import shapefile

#%% Whitmeyer and Karlstrom stuff
def shp_to_xy(shape):
    xall = np.zeros(len(shape))
    yall = np.zeros(len(shape))
    for ipt in range(len(shape)):
        x, y = shape[ipt]
        xall[ipt] = x
        yall[ipt] = y
    return xall, yall

def xy_many_shapes(shapes, savename):
    for ishape, shape in enumerate(shapes):
        x, y = shp_to_xy(shape.points)
        # np.savetxt('%s_{3.0f}.txt'.format(savename, ishape), np.array([x.T,y.T]).T )
        # np.savetxt('%s_{3.0f}.txt'.format(savename, ishape), np.array([x.T,y.T]).T )
        np.savetxt(
            f'{savename}_{ishape}.txt',
            np.array([x.T, y.T]).T
        )

# Mid-continent rift
shapes = shapefile.Reader('MCR.shx').shapes()
xy_many_shapes(shapes, 'MCR')

# # Reelfoot rift
shapes = shapefile.Reader('Reelfoot.shx').shapes()
xy_many_shapes(shapes, 'Reelfoot')

# Grenville front
shapes = shapefile.Reader('grv_frt.shx').shapes()
xy_many_shapes(shapes, 'grv_frt')

# ?
shapes = shapefile.Reader('something_province.shx').shapes() # This is misnamed, it is thrusts. But I can't change file names or I get a glitch.
xy_many_shapes(shapes, 'something_province')
