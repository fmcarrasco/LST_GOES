import numpy as np
import numpy.ma as ma
import pandas as pd
import datetime as dt

import sys
sys.path.append('./lib')
from lst_class_horario import lst_class_horario
from lst_class_diario import lst_class_diario
from figure_functions import mapa_base, colores_tmin, legend_temp

import matplotlib.pyplot as plt
#import cartopy.crs as ccrs

#Objeto para leer IMERG
'''
fecha = '202311231000' # Formato: yyyymmddHHMM (La hora en UTC)
a = lst_class_horario(fecha, './salidas/')
print(a.archivo)
print(a.carpeta)
print(a.H)
a.save_map_lst(opt=1)
a.save_map_lst(opt=2)
a.save_map_DQF()

fecha = '20230810'
b = lst_class_diario(fecha, './salidas/')
print(b.LSTmin)
b.mapa_helada_ORA()
b.mapa_helada_diario()


for hora in range(11,12):
    print(hora)
    fecha = '20220622'+str(hora).zfill(2)+'00'
    a = lst_class_horario(fecha, './salidas/')
    a.save_gtiff(file_name=fecha+'.tiff')
'''

fechad = '20231012'
b = lst_class_diario(fechad, './salidas/')
#b.save_gtiff(file_name='LSTmin_' + fechad + '.tiff')
b.mapa_helada_ORA()
b.mapa_helada_diario()


