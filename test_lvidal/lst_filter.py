import gc
import glob
import s3fs
import netCDF4
import numpy as np
import pandas as pd
import matplotlib as mpl
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from datetime import datetime, timedelta, timezone

import matplotlib.colors as mcolors
import matplotlib.cm as cm

import scipy.ndimage

import cartopy
import cartopy.crs as ccrs
from pyproj import Proj

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

import matplotlib.dates as mdates
from datetime import datetime as dt

from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

#------------------------------------------------------------------------------#

# Defino los limites de paises (descarga por unica vez)
provincias = cartopy.feature.NaturalEarthFeature(category='cultural',
                                                 name='admin_1_states_provinces_lines',
                                                 scale='10m',
                                                 facecolor='none')

paises = cartopy.feature.NaturalEarthFeature(category='cultural',
                                             name='admin_0_countries',
                                             scale='10m',
                                             facecolor='none')

#------------------------------------------------------------------------------#

# PALETA DE COLORES PRODUCTO LST (Autor: DSS/SMN)
lst_colors = [ [255/255,255/255,255/255], # +12 to +10°C
               [ 74/255,  8/255,143/255], # +10 to  +8°C
               [ 51/255,  4/255,180/255], #  +8 to  +6°C
               [ 45/255, 46/255,255/255], #  +6 to  +4°C
               [ 38/255,105/255,247/255], #  +4 to  +2°C
               [  2/255,128/255,250/255], #  +2 to   0°C
               [ 82/255,176/255,250/255], #   0 to  -2°C
               [132/255,215/255,249/255], #  -2 to  -4°C
               [160/255,227/255,196/255], #  -4 to  -6°C
               [182/255,235/255,199/255], #  -6 to  -8°C
               [207/255,236/255,204/255], #  -8 to -10°C
               [221/255,242/255,209/255], # -10 to -12°C
               [233/255,242/255,209/255], # -12 to -14°C
               [218/255,220/255,217/255]  # -14 to -16°C
             ]
cmap_lst = mcolors.ListedColormap(lst_colors)
bounds = np.array([-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16])
norm_topo = mcolors.BoundaryNorm(boundaries=bounds, ncolors=14)

#------------------------------------------------------------------------------#

import time
t = time.time()

# ROI (Argentina)
x1_lon = 520
x2_lon = 785 #765
y1_lat = 650 #735 
y2_lat = 1025

# FECHA Y HORA DE INTERES PARA ANALIZAR
# NOTA: este producto se disponibiliza uno por hora
fecha = '202311221000' # Formato: yyyymmddHHMM (La hora en UTC)

fecha0 = datetime.strptime(fecha, '%Y%m%d%H%M')
print('Fecha y hora de interés: ' + fecha0.strftime('%Y/%m/%d %H:%M UTC'))
print('======================================================================')

fs = s3fs.S3FileSystem(anon=True)
files = fs.ls('noaa-goes16/ABI-L2-LSTF/'+fecha0.strftime('%Y')+'/'+fecha0.strftime('%j')+'/'+fecha0.strftime('%H')+'/')
fname=files[0]
print(fname)

print('======================================================================')

with fs.open(fname) as f:

  with netCDF4.Dataset(fname.split('/')[-1], memory=f.read()) as goes:

    # print(goes)

    H = goes.variables['goes_imager_projection'].getncattr('perspective_point_height')
    lon_0 = goes.variables['goes_imager_projection'].getncattr('longitude_of_projection_origin')
    sat_sweep = goes.variables['goes_imager_projection'].getncattr('sweep_angle_axis')
    x = goes.variables['x'][x1_lon:x2_lon] * H
    y = goes.variables['y'][y1_lat:y2_lat] * H
    xv, yv = np.meshgrid(np.array(x), np.array(y))
    # Doc: https://proj.org/operations/projections/geos.html
    geo = Proj(proj='geos', h=H, lon_0=lon_0, sweep=sat_sweep)
    lons_lst, lats_lst = geo(xv, yv, inverse=True)
    print(lons_lst)
    print(lats_lst)

    LST = goes.variables['LST'][y1_lat:y2_lat,x1_lon:x2_lon][::1,::1]
    DQF = goes.variables['DQF'][y1_lat:y2_lat,x1_lon:x2_lon][::1,::1]
    PQI = goes.variables['PQI'][y1_lat:y2_lat,x1_lon:x2_lon][::1,::1]

    add_seconds = int(goes.variables['time_bounds'][0])
    date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)

elapsed = time.time() - t
print('Tiempo de lectura: '+str(np.round(elapsed,1))+' sec')



########################

# Definimos proyeccion
projection = ccrs.PlateCarree()
#extent = [-76.5, -52, -56.5, -20.]
extent = [-76.5, -38., -56.5, -9.]

fig=plt.figure(figsize=[11,10])

ax = fig.add_subplot(111, projection=projection)
ax.set_extent(extents=extent, crs=projection)

img=plt.pcolormesh(lons_lst,
                   lats_lst,
                   LST-273,
                   vmin=-16,
                   vmax=12,
                   cmap=cmap_lst,
                   transform=ccrs.PlateCarree())
plt.colorbar(img, pad=0.01, aspect=42, shrink=0.5, ticks=np.arange(-16,12.1,2))

ax.add_feature(provincias, linewidth=0.25, edgecolor='k', facecolor='None')
ax.add_feature(paises, linewidth=1.0, edgecolor='k', facecolor='None')

# Agregamos lineas y marcas lat/lon
gl = ax.gridlines(crs=ccrs.PlateCarree(),
                  color='gray',
                  alpha=1.0,
                  linestyle='--',
                  linewidth=0.1,
                 xlocs=np.arange(-180, 180, 5),
                 ylocs=np.arange(-90, 90, 5),
                  draw_labels=True)
gl.top_labels = False
gl.right_labels = False

plt.title('GOES-16 LST', fontsize=14, loc='left')
plt.title(date.strftime('%Y-%m-%d %H:%M UTC'), fontsize=14, loc='right')

plt.show()

###############################################################
# Paleta de 4 colores para indicar la calidad del dato LST
dqf_colors = [
               [ 51/255,204/255,  0/255], # Alta
               [255/255,153/255,  0/255], # Media
               [255/255,  0/255,  0/255], # Baja
               [218/255,220/255,217/255]  # No Data
             ]
cmap_dqf = mcolors.ListedColormap(dqf_colors)

# Definimos proyeccion
projection = ccrs.PlateCarree()
#extent = [-76.5, -52, -56.5, -20.]
extent = [-76.5, -38., -56.5, -9.]
fig=plt.figure(figsize=[11,10])

ax = fig.add_subplot(111, projection=projection)
ax.set_extent(extents=extent, crs=projection)

img=plt.pcolormesh(lons_lst,
                   lats_lst,
                   DQF,
                   vmin=0,
                   vmax=4,
                   cmap=cmap_dqf,
                   transform=ccrs.PlateCarree())
cbar = plt.colorbar(img, pad=0.01, aspect=18, shrink=0.5)
cbar.ax.get_yaxis().set_ticks([])
for j, lab in enumerate(['$ALTA$','$MEDIA$','$BAJA$','$ND$']):
  cbar.ax.text(.5, j + 0.5, lab, ha='center', va='center', color='k', rotation=270)

ax.add_feature(provincias, linewidth=0.25, edgecolor='k', facecolor='None')
ax.add_feature(paises, linewidth=1.0, edgecolor='k', facecolor='None')

# Agregamos lineas y marcas lat/lon
gl = ax.gridlines(crs=ccrs.PlateCarree(),
                  color='gray',
                  alpha=1.0,
                  linestyle='--',
                  linewidth=0.1,
                 xlocs=np.arange(-180, 180, 5),
                 ylocs=np.arange(-90, 90, 5),
                  draw_labels=True)
gl.top_labels = False
gl.right_labels = False

plt.title('GOES-16 LST Calidad [DQF]', fontsize=14, loc='left')
plt.title(date.strftime('%Y-%m-%d %H:%M UTC'), fontsize=14, loc='right')

plt.show()


#==============================================================================#

# Lo primero que hago es crea una mascara que tendra 0 y 1, del mismo tamaño
# que la matriz DQF
mask_lst = np.copy(DQF)

# Asigno 1 donde DQF==0 (calidad buena de la estimacion LST)
mask_lst[DQF<0.5]=1

# Asigno 0 donde la calidad de la estimacion de LST es media o baja.
mask_lst[DQF>=0.5]=0

#==============================================================================#

# Definimos proyeccion
projection = ccrs.PlateCarree()
#extent = [-76.5, -52, -56.5, -20.]
extent = [-76.5, -38., -56.5, -9.]
fig=plt.figure(figsize=[11,10])

ax = fig.add_subplot(111, projection=projection)
ax.set_extent(extents=extent, crs=projection)

img=plt.pcolormesh(lons_lst,
                   lats_lst,
                   mask_lst,
                   vmin=0,
                   vmax=1,
                   transform=ccrs.PlateCarree())
plt.colorbar(img, pad=0.01, aspect=40, shrink=0.5)

ax.add_feature(provincias, linewidth=0.25, edgecolor='k', facecolor='None')
ax.add_feature(paises, linewidth=1.0, edgecolor='k', facecolor='None')

# Agregamos lineas y marcas lat/lon
gl = ax.gridlines(crs=ccrs.PlateCarree(),
                  color='gray',
                  alpha=1.0,
                  linestyle='--',
                  linewidth=0.1,
                 xlocs=np.arange(-180, 180, 5),
                 ylocs=np.arange(-90, 90, 5),
                  draw_labels=True)
gl.top_labels = False
gl.right_labels = False

plt.title('GOES-16 LST mask', fontsize=14, loc='left')
plt.title(date.strftime('%Y-%m-%d %H:%M UTC'), fontsize=14, loc='right')

plt.show()



##################


# DATO LST FILTRADO
LST_bueno = LST * mask_lst

# Definimos proyeccion
projection = ccrs.PlateCarree()
#extent = [-76.5, -52, -56.5, -20.]
extent = [-76.5, -38., -56.5, -9.]
fig=plt.figure(figsize=[11,10])

ax = fig.add_subplot(111, projection=projection)
ax.set_extent(extents=extent, crs=projection)

img=plt.pcolormesh(lons_lst,
                   lats_lst,
                   LST_bueno-273,
                   vmin=-16,
                   vmax=12,
                   cmap=cmap_lst,
                   transform=ccrs.PlateCarree())
plt.colorbar(img, pad=0.01, aspect=42, shrink=0.5, ticks=np.arange(-16,12.1,2))

ax.add_feature(provincias, linewidth=0.25, edgecolor='k', facecolor='None')
ax.add_feature(paises, linewidth=1.0, edgecolor='k', facecolor='None')

# Agregamos lineas y marcas lat/lon
gl = ax.gridlines(crs=ccrs.PlateCarree(),
                  color='gray',
                  alpha=1.0,
                  linestyle='--',
                  linewidth=0.1,
                 xlocs=np.arange(-180, 180, 5),
                 ylocs=np.arange(-90, 90, 5),
                  draw_labels=True)
gl.top_labels = False
gl.right_labels = False

plt.title('GOES-16 LST filtrado', fontsize=14, loc='left')
plt.title(date.strftime('%Y-%m-%d %H:%M UTC'), fontsize=14, loc='right')

plt.show()