import datetime as dt
import pandas as pd
import numpy as np
import os
import sys
import sys
sys.path.append('./lib')
from lst_class import lst_class
from figure_functions import mapa_temperatura_minima_superficial

#
f1 = sys.argv[1]
fecha_i = dt.datetime.strptime(f1, '%Y%m%d')
f2 = sys.argv[2]
fecha_f = dt.datetime.strptime(f2, '%Y%m%d')
fechas = pd.date_range(fecha_i, fecha_f)
nt = len(fechas)
#
carpeta = './salidas/'
os.makedirs(carpeta, exist_ok=True)
#

im_data = lst_class(f1, f2 )

if im_data.SinDatos:
    print('No hay Archivos TIFF para el periodo correspondiente!')
    exit()

x = im_data.lon
y = im_data.lat
LSTmin = np.nanmin(im_data.lst_min, axis=0)
l_lat = [-60., -18.]
l_lon = [-79., -47.]

nfig = carpeta + 'LSTmin_' + fecha_i.strftime('%Y%m%d') + '_' + fecha_f.strftime('%Y%m%d') + '.jpg'

mapa_temperatura_minima_superficial(x, y, LSTmin, l_lat, l_lon, nfig)