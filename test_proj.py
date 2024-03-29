import numpy as np
import numpy.ma as ma
import pandas as pd
import datetime as dt

import sys
sys.path.append('./lib')
from lst_class_horario_consultoria import lst_horario
from lst_class_diario import lst_class_diario
from figure_functions import mapa_base, colores_tmin, legend_temp

from pyproj import Transformer, transform

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

def get_nearest_index(ds, loni, lati):
        # Transformar de EPSG:4326 a ECEF (EPSG:4978)
        # GtiFF viene en EPSG:4326 y pasamos a coord cartesianas EPSG:4978
        transformer = Transformer.from_crs(4326, 4978)
        X, Y = transformer.transform(ds.lons, ds.lats)
        x, y = transformer.transform(loni,lati)
        distance = ((x-X)**2 + (y-Y)**2)**0.5
        pos_index = np.unravel_index(distance.argmin(), distance.shape)

        return pos_index

def subset_goes(ds, lon_a, lat_a):
        lon1, lon2 = lon_a
        lat1, lat2 = lat_a
        x1, y1 = get_nearest_index(ds, lon1, lat1)
        x2, y2 = get_nearest_index(ds, lon2, lat1)
        x3, y3 = get_nearest_index(ds, lon2, lat2)
        x4, y4 = get_nearest_index(ds, lon1, lat2)
        x = [x1, x2, x3, x4]
        y = [y1, y2, y3, y4]
        return [min(x), max(x), min(y), max(y)]





#Objeto para leer IMERG

fecha = '201705280500' # Formato: yyyymmddHHMM (La hora en UTC)
print(fecha)
a = lst_horario(fecha, './salidas/')
print(a.LST.shape)
exit()
# Punto de extracción test
lati = -37.46
loni = -59.09
# Area de extraccion de datos
lon_a = [-61.5, -58.5]
lat_a = [-38.1, -36.1]
x0, y0 = get_nearest_index(a, loni, lati)
x1, x2, y1, y2 = subset_goes(a, lon_a, lat_a)
lats = a.lats[x1:x2,y1:y2]
lons = a.lons[x1:x2,y1:y2]
dato = a.LST_filtrado[x1:x2,y1:y2]

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

##########################################################################
############## Figura 1 ##################################################
fig=plt.figure(figsize=[11,10])

ax = fig.add_subplot(111, projection=projection)
ax.set_extent(extents=extent, crs=projection)

img = ax.pcolormesh(a.lons, a.lats, a.LST_filtrado, vmin=-14, vmax=10, cmap=cmap_lst,
                           transform=ccrs.PlateCarree())
ax.plot([lon_a[0], lon_a[1], lon_a[1], lon_a[0], lon_a[0]], [lat_a[0], lat_a[0], lat_a[1], lat_a[1], lat_a[0]],
         color='black', linewidth=1, marker='.',
         transform=ccrs.Geodetic(), #remove this line to get straight lines
         )

plt.colorbar(img, pad=0.01, aspect=42, shrink=0.5, ticks=np.arange(-14,10.1,2), extend='both')
ax.set_title('Temperatura de superficie GOES-16 23/11/2023 05 UTC (2 AM Hora Local)', loc='left')
ax.add_feature(provincias, linewidth=0.25, edgecolor='k', facecolor='None')
ax.add_feature(paises, linewidth=1.0, edgecolor='k', facecolor='None')
ax.add_feature(shape_feature1, linestyle='-', linewidth=0.3)
plt.savefig('./salidas/Figura1.jpg', dpi=150, bbox_inches='tight')
plt.close(fig)

##########################################################################
############## Figura 2 ##################################################
fig=plt.figure(figsize=[11,10])

ax = fig.add_subplot(111, projection=projection)
ax.set_extent(extents=extent, crs=projection)

img = ax.pcolormesh(lons, lats, dato, vmin=-14, vmax=10, cmap=cmap_lst,
                           transform=ccrs.PlateCarree())
ax.plot([lon_a[0], lon_a[1], lon_a[1], lon_a[0], lon_a[0]], [lat_a[0], lat_a[0], lat_a[1], lat_a[1], lat_a[0]],
         color='black', linewidth=1, marker='.',
         transform=ccrs.Geodetic(), #remove this line to get straight lines
         )

plt.colorbar(img, pad=0.01, aspect=42, shrink=0.5, ticks=np.arange(-14,10.1,2), extend='both')
ax.set_title('Temperatura de superficie GOES-16 23/11/2023 05 UTC (2 AM Hora Local)', loc='left')
ax.add_feature(provincias, linewidth=0.25, edgecolor='k', facecolor='None')
ax.add_feature(paises, linewidth=1.0, edgecolor='k', facecolor='None')
ax.add_feature(shape_feature1, linestyle='-', linewidth=0.3)
plt.savefig('./salidas/Figura2.jpg', dpi=150, bbox_inches='tight')
plt.close(fig)

##########################################################################
############## Figura 3 ##################################################
# Azul, Benito Juarez, Olavarria, Tandil
lat_est = [-36.833333, -37.716667, -36.883333, -37.233333]
lon_est = [-59.883333, -59.783333, -60.216667, -59.25]

lat_recorte = lats.ravel()
lon_recorte = lons.ravel()

fig=plt.figure(figsize=[11,10])

ax = fig.add_subplot(111, projection=projection)
ax.set_extent(extents=extent, crs=projection)


ax.scatter(lon_est, lat_est, 13, marker='s', color='r', label='Estaciones', zorder=2)

ax.scatter(lon_recorte, lat_recorte, 0.5, color='b', transform=ccrs.Geodetic(), label=u'lat/lon retícula', zorder=0)


ax.plot([lon_a[0], lon_a[1], lon_a[1], lon_a[0], lon_a[0]], [lat_a[0], lat_a[0], lat_a[1], lat_a[1], lat_a[0]],
         color='black', linewidth=1, marker='.',
         transform=ccrs.Geodetic(), #remove this line to get straight lines
         label=u'Área de interés', zorder=1)
ax.legend()
ax.set_title('Area de recorte, estaciones (rojo) y puntos de reticula (azul)', loc='left')
ax.add_feature(provincias, linewidth=0.25, edgecolor='k', facecolor='None')
ax.add_feature(paises, linewidth=1.0, edgecolor='k', facecolor='None')
ax.add_feature(shape_feature1, linestyle='-', linewidth=0.3)
plt.savefig('./salidas/Figura3.jpg', dpi=150, bbox_inches='tight')
plt.close(fig)


##########################################################################
############## Figura 4 ##################################################
# Azul, Benito Juarez, Olavarria, Tandil
lat_est = [-36.833333, -37.716667, -36.883333, -37.233333]
lon_est = [-59.883333, -59.783333, -60.216667, -59.25]

lat_recorte = lats.ravel()
lon_recorte = lons.ravel()

lon_pix = [lons[5,9], lons[8,6], lons[11,11], lons[11,21], lons[8,16], lons[5,19]]
lat_pix = [lats[5,9], lats[8,6], lats[11,11], lats[11,21], lats[8,16], lats[5,19]]
numeros = [str(i+1) for i in np.arange(len(lon_pix))]

fig=plt.figure(figsize=[11,10])

ax = fig.add_subplot(111, projection=projection)
ax.set_extent(extents=extent, crs=projection)


ax.scatter(lon_est, lat_est, 13, marker='s', color='r', label='Estaciones', zorder=2)

ax.scatter(lon_pix, lat_pix, 0.5, color='b', transform=ccrs.Geodetic(), label=u'lat/lon retícula', zorder=0)
ax.scatter(lon_pix, lat_pix, 12, marker='s', color='y', zorder=-1)

for i, txt in enumerate(numeros):
    if i == 3:
          ax.annotate(txt, (lon_pix[i] + 0.09, lat_pix[i]))
    else:
          ax.annotate(txt, (lon_pix[i] - 0.16, lat_pix[i]))


ax.plot([lon_a[0], lon_a[1], lon_a[1], lon_a[0], lon_a[0]], [lat_a[0], lat_a[0], lat_a[1], lat_a[1], lat_a[0]],
         color='black', linewidth=1, marker='.',
         transform=ccrs.Geodetic(), #remove this line to get straight lines
         label=u'Área de interés', zorder=1)
ax.legend()
ax.set_title('Area de recorte, estaciones (rojo) y puntos de reticula seleccionados (amarillo)', loc='left')
ax.add_feature(provincias, linewidth=0.25, edgecolor='k', facecolor='None')
ax.add_feature(paises, linewidth=1.0, edgecolor='k', facecolor='None')
ax.add_feature(shape_feature1, linestyle='-', linewidth=0.3)
plt.savefig('./salidas/Figura4.jpg', dpi=150, bbox_inches='tight')
plt.close(fig)