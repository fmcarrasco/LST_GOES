import time
import numpy as np
import numpy.ma as ma
import pandas as pd
import datetime as dt

from pyproj import Transformer, transform

import matplotlib.pyplot as plt

import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

from cartopy.io.shapereader import natural_earth


import sys
sys.path.append('./lib')
from lst_class_horario_consultoria import lst_horario

def get_nearest_index(lons, lats, loni, lati):
        # Transformar de EPSG:4326 a ECEF (EPSG:4978)
        # GtiFF viene en EPSG:4326 y pasamos a coord cartesianas EPSG:4978
        transformer = Transformer.from_crs(4326, 4978)
        X, Y = transformer.transform(lons, lats)
        x, y = transformer.transform(loni,lati)
        distance = ((x-X)**2 + (y-Y)**2)**0.5
        pos_index = np.unravel_index(distance.argmin(), distance.shape)

        return pos_index

def subset_goes(lons, lats, lon_a, lat_a):
        lon1, lon2 = lon_a
        lat1, lat2 = lat_a
        x1, y1 = get_nearest_index(lons, lats, lon1, lat1)
        x2, y2 = get_nearest_index(lons, lats, lon2, lat1)
        x3, y3 = get_nearest_index(lons, lats, lon2, lat2)
        x4, y4 = get_nearest_index(lons, lats, lon1, lat2)
        x = [x1, x2, x3, x4]
        y = [y1, y2, y3, y4]
        return [min(x), max(x), min(y), max(y)]

start = time.time()

# Datos para la figura
fname1 = '../SMN_SQPE/Mapas_semanales/IGN/departamento/departamento.shp'
shape_feature1 = ShapelyFeature(Reader(fname1).geometries(),
                                ccrs.PlateCarree(), facecolor='none')
provincias = cartopy.feature.NaturalEarthFeature(category='cultural',
                                                 name='admin_1_states_provinces_lines',
                                                 scale='10m',
                                                 facecolor='none')
paises = cartopy.feature.NaturalEarthFeature(category='cultural',
                                             name='admin_0_countries',
                                             scale='10m',
                                             facecolor='none')
lst_colors = [ [ 74/255,  8/255,143/255], # +10 to  +8°C
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
               ]
lst_colors = ["#4B088A","#3104B4","#2E2EFE","#2E64FE",
                      "#0080FF","#58ACFA","#81DAF5","#A3E2C3",
                      "#BAE9C6","#CDEECA","#DEF1D0","#EAF2D6"]
cmap_lst = mcolors.ListedColormap(lst_colors)
cmap_lst.set_bad(color='white')
cmap_lst.set_extremes(under='#4B088A', over='#DADCDA')
bounds = np.array([-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14])
norm_topo = mcolors.BoundaryNorm(boundaries=bounds, ncolors=12)
##################
# Definimos proyeccion
projection = ccrs.PlateCarree()
extent = [-63.94, -56.54, -41.43, -33.07]




mes = 4
nmes = 'Abril'

lon_a = [-61.5, -58.5]
lat_a = [-38.1, -36.1]

if mes in [10,11,12]:
    years = np.arange(2019,2024)
elif mes in [1,2,3,4]:
    years = np.arange(2019,2025)

for year in years:
    print('Trabajando el año', year)
    ini = str(year)+'0401'
    fin = str(year)+'0501'
    fechas = pd.date_range(start=ini, end=fin, freq='H')
    fechasd = pd.date_range(start=ini, end=fin)
    #print((fechasd[0:-1]))
    ntime = len(fechasd[0:-1])
    ny = 290
    nx = 245
    print('./salidas/LST_' + nmes + '_' + str(year) + '.npy')
    LST = np.load('./salidas/LST_' + nmes + '_' + str(year) + '.npy')
    DQF = np.load('./salidas/DQF_' + nmes + '_' + str(year) + '.npy')
    lats = np.load('./salidas/lats.npy')
    lons = np.load('./salidas/lons.npy')
    LST[DQF>0] = np.nan # Filtramos aquellos datos con problemas
    # Area de extraccion de datos
    x1, x2, y1, y2 = subset_goes(lons, lats, lon_a, lat_a)
    latsub = lats[x1:x2,y1:y2]
    lonsub = lons[x1:x2,y1:y2]
    dato = LST[:, x1:x2,y1:y2]
    #
    LSTmin = np.empty((ntime, ny, nx))
    LSTmin[:] = np.nan
    conteo = np.empty((ntime, ny, nx))
    conteo[:] = np.nan
    for dia, it in enumerate(np.arange(0,len(fechas[0:-1]),24)):
          conteo[dia,:,:] = np.count_nonzero(~np.isnan(LST[it:it+13,:,:]), axis=0)
          LSTmin[dia,:,:] = np.nanmin(LST[it:it+13,:,:], axis=0)
    
    nfile0 = './salidas/LSTmin_' + nmes + '_' + str(year) + '.npy'
    nfile1 = './salidas/COUNTmin_' + nmes + '_' + str(year) + '.npy'
    np.save(nfile0, LSTmin)
    np.save(nfile1, conteo)
         

    #LSTmin_sub = LSTmin[:,x1:x2,y1:y2]
    #conteo_sub = conteo[:,x1:x2,y1:y2]
    #print(LSTmin_sub[:,5,9])
    #print(conteo_sub[:,5,9])

'''
fig=plt.figure(figsize=[11,10])
ax = fig.add_subplot(111, projection=projection)
ax.set_extent(extents=extent, crs=projection)
img = ax.pcolormesh(lonsub, latsub, LSTmin[20,x1:x2,y1:y2], vmin=-14, vmax=10, cmap=cmap_lst,
                    transform=ccrs.PlateCarree())
#img = ax.pcolormesh(lons, lats, LST[500,:,:], vmin=-14, vmax=10, cmap=cmap_lst,
#                    transform=ccrs.PlateCarree())
ax.plot([lon_a[0], lon_a[1], lon_a[1], lon_a[0], lon_a[0]], [lat_a[0], lat_a[0], lat_a[1], lat_a[1], lat_a[0]],
        color='black', linewidth=1, marker='.',transform=ccrs.Geodetic(), #remove this line to get straight lines
        )
plt.colorbar(img, pad=0.01, aspect=42, shrink=0.5, ticks=np.arange(-14,10.1,2), extend='both')
#ax.set_title('Temperatura de superficie GOES-16 23/11/2023 05 UTC (2 AM Hora Local)', loc='left')
ax.add_feature(provincias, linewidth=0.25, edgecolor='k', facecolor='None')
ax.add_feature(paises, linewidth=1.0, edgecolor='k', facecolor='None')
ax.add_feature(shape_feature1, linestyle='-', linewidth=0.3)
plt.show()
#plt.savefig('./salidas/Figura2.jpg', dpi=150, bbox_inches='tight')
#plt.close(fig)
''' 

    



end = time.time()

secs = np.round(end - start, 2)
print('Se demoro ', secs, ' segundos.')