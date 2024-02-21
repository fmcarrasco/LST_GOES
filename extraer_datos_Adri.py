import time
import numpy as np
import numpy.ma as ma
import pandas as pd
import datetime as dt

from pyproj import Transformer, transform

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


mes = 10
nmes = 'Octubre'

lon_a = [-61.5, -58.5]
lat_a = [-38.1, -36.1]

if mes in [10,11,12]:
    years = np.arange(2019,2024)
elif mes in [1,2,3,4]:
    years = np.arange(2019,2025)

pares = [(5,9), (8,6), (11,11), (11,21), (8,16), (5,19)]
nomb = ['LST pixel 1', 'conteo pixel 1', 'LST pixel 2', 'conteo pixel 2', 'LST pixel 3', 'conteo pixel 3',
        'LST pixel 4', 'conteo pixel 4', 'LST pixel 5', 'conteo pixel 5', 'LST pixel 6', 'conteo pixel 6'  ]


lista_lst = []
lista_cnt = []
lista_fec = []
for year in years:
    print('Trabajando el año', year)
    ini = str(year)+'1001'
    fin = str(year)+'1101'
    fechasd = pd.date_range(start=ini, end=fin)
    print((fechasd[0:-1]))
    LSTmin = np.load('./salidas/LSTmin_' + nmes + '_' + str(year) + '.npy')
    conteo = np.load('./salidas/COUNTmin_' + nmes + '_' + str(year) + '.npy')
    
    lista_lst.append(LSTmin)
    lista_cnt.append(conteo)
    lista_fec.append(fechasd[0:-1])
LSTmin = np.concatenate(lista_lst, axis=0)
conteo = np.concatenate(lista_cnt, axis=0)
indice = lista_fec[0].union(lista_fec[1]).union(lista_fec[2]).union(lista_fec[3]).union(lista_fec[4])

#lista_fec[0].union(lista_fec[1:])
#print(lista_fec[0])
lats = np.load('./salidas/lats.npy')
lons = np.load('./salidas/lons.npy')
# Area de extraccion de datos
x1, x2, y1, y2 = subset_goes(lons, lats, lon_a, lat_a)
latsub = lats[x1:x2,y1:y2]
lonsub = lons[x1:x2,y1:y2]
LSTmin_sub = LSTmin[:, x1:x2,y1:y2]
conteo_sub = conteo[:, x1:x2,y1:y2]
# preparamos el Dataframe para excel
test = np.empty((len(years)*31, 12))
test[:] = np.nan
for i, par in enumerate(pares):
    test[:,i*2] = LSTmin_sub[:,par[0], par[1] ]
    test[:,i*2+1] = conteo_sub[:,par[0], par[1] ]

df = pd.DataFrame(index=indice, data=test, columns= nomb)
df.to_excel('./salidas/LST_pixels.xlsx')


end = time.time()

secs = np.round(end - start, 2)
print('Se demoro ', secs, ' segundos.')