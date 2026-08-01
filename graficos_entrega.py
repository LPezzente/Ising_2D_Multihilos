'''
Created on 14/07/2018

@author: leandro
'''

from scipy.special._ufuncs import ellipe
from scipy.special.basic import ellipk

import matplotlib.pyplot as plt
import numpy as np


def plotearCurvas( titulo, leyendaEjeX, leyendaEjeY, minEjeX, maxEjeX , vectorCurvaX, vectorCurvaY, configuracionPlot ):
    fig, ax = plt.subplots()
    ax.cla()
    ax.set_title( titulo )
    ax.set_xlim( minEjeX, maxEjeX )
    ax.set_xlabel( leyendaEjeX )
    ax.set_ylabel( leyendaEjeY )
    ax.minorticks_off()
    ax.grid( True )
    ax.plot( vectorCurvaX , vectorCurvaY , configuracionPlot )
    plt.savefig( leyendaEjeY, fmt = "png" )
    fig.canvas.draw()

def plotear_mag( N, data, w, h, img_dpi ):

    # Datos del plot
    vectorCurvaX = data[:, 0]
    vectorCurvaY = data[:, 1]
    vectorCurvaY = np.divide( vectorCurvaY, N );

    # Parametros del pLot
    minEjeX = min( vectorCurvaX )
    maxEjeX = max( vectorCurvaX )

    # Configuracion del plot
    configuracionPlot = "o-"

    # Texto del grafico
    titulo = "Magnetizacion vs T"
    leyendaEjeX = 'Temperatura'
    leyendaEjeY = 'Magnetizacion'

    onsager = np.power ( np.ones( len( vectorCurvaX ) ) - np.power( np.sinh( np.divide( 2.0 , vectorCurvaX ) ) , -4.0 ) , 1.0 / 8.0 )
    for index in range( len( onsager ) ):
        if np.isnan( onsager[index] ):
            onsager[index] = 0
        else:
            pass

    campo_medio = np.zeros( len( vectorCurvaX ) )
    for index in range( len( campo_medio ) ):
        m = 0.5
        T_iter = vectorCurvaX[index]
        for iter in range( 1000 ):
            m = np.tanh( 4 * m / T_iter )
        campo_medio[index] = m

    fig, ax = plt.subplots( dpi = img_dpi )
    DPI = fig.get_dpi()
    fig.set_size_inches( h / float( DPI ), w / float( DPI ) )
    ax.cla()
    ax.set_title( titulo )
    ax.set_xlim( minEjeX, maxEjeX )
    ax.set_xlabel( leyendaEjeX )
    ax.set_ylabel( leyendaEjeY )
    ax.minorticks_off()
    ax.grid( True )

    ax.set_ylim( [0, 1.1] )
    ax.plot( vectorCurvaX , vectorCurvaY , 'D', color = "lightblue" , label = 'Metropolis' )
    ax.plot( vectorCurvaX , onsager , '-', color = "red" , label = 'Onsager' )
    ax.plot( vectorCurvaX , campo_medio , '--', color = "blue" , label = 'Campo Medio' )
    ax.legend( loc = 'upper right' )

    plt.savefig( leyendaEjeY, fmt = "png" )
    fig.canvas.draw()

def plotear_energia( N, data, w, h, img_dpi ):

    # Datos del plot
    vectorCurvaX = data[:, 0]
    vectorCurvaY = data[:, 1]
    vectorCurvaY = np.divide( vectorCurvaY, N );

    # Parametros del pLot
    minEjeX = min( vectorCurvaX )
    maxEjeX = max( vectorCurvaX )

    # Configuracion del plot
    configuracionPlot = "o-"

    # Texto del grafico
    titulo = "Energia vs T"
    leyendaEjeX = 'Temperatura'
    leyendaEjeY = 'Energia'

    onsager = np.zeros( len( vectorCurvaX ) )
    for index in range( len( onsager ) ):
        T = vectorCurvaX[index]
        alpha = np.power( np.tanh( 2.0 / T ) , -1 )
        k0 = ( 2.0 * np.sinh( 2.0 / T ) ) / np.power( ( np.cosh( 2.0 / T ) ) , 2 )
        k1 = 2.0 * np.power( np.tanh( 2.0 / T ) , 2 ) - 1.0
        onsager[index] = -alpha * ( 1 + ( 2.0 / np.pi ) * k1 * ellipk( k0 ) )

    mag_cm = np.zeros( len( vectorCurvaX ) )
    for index in range( len( mag_cm ) ):
       m = 0.5
       T_iter = vectorCurvaX[index]
       for iter in range( 1000 ):
           m = np.tanh( 4 * m / T_iter )
       mag_cm[index] = m

    campo_medio = np.zeros( len( vectorCurvaX ) )
    for index in range( len( campo_medio ) ):
        T = vectorCurvaX[index]
        m = mag_cm[index]
        campo_medio[index] = -2.0 * np.power( m , 2 )

    fig, ax = plt.subplots( dpi = img_dpi )
    DPI = fig.get_dpi()
    fig.set_size_inches( h / float( DPI ), w / float( DPI ) )
    ax.cla()
    ax.set_title( titulo )
    ax.set_xlim( minEjeX, maxEjeX )
    ax.set_xlabel( leyendaEjeX )
    ax.set_ylabel( leyendaEjeY )
    ax.minorticks_off()
    ax.grid( True )

    ax.plot( vectorCurvaX , vectorCurvaY , 'D', color = "lightblue" , label = 'Metropolis' )
    ax.plot( vectorCurvaX , onsager , '-', color = "red" , label = 'Onsager' )
    ax.plot( vectorCurvaX , campo_medio , '--', color = "blue" , label = 'Campo Medio' )
    ax.legend( loc = 'upper left' )

    plt.savefig( leyendaEjeY, fmt = "png" )
    fig.canvas.draw()

def plotear_cv( N, data, w, h, img_dpi ):

    # Datos del plot
    vectorCurvaX = data[:, 0]
    vectorCurvaY = data[:, 1]
    vectorCurvaY = np.divide( vectorCurvaY, N );

    # Parametros del pLot
    minEjeX = min( vectorCurvaX )
    maxEjeX = max( vectorCurvaX )

    # Configuracion del plot
    configuracionPlot = "o-"

    # Texto del grafico
    titulo = "Capacidad Calorifica vs T"
    leyendaEjeX = 'Temperatura'
    leyendaEjeY = 'Capacidad Calorifica'

    onsager = np.zeros( len( vectorCurvaX ) )
    for index in range( len( onsager ) ):
        T = vectorCurvaX[index]
        alpha = ( 2.0 / np.pi ) * np.power( T, -2 ) * np.power( np.tanh( 2.0 / T ) , -2 )
        k0 = ( 2.0 * np.sinh( 2.0 / T ) ) / np.power( ( np.cosh( 2.0 / T ) ) , 2 )
        k1 = 2.0 * np.power( np.tanh( 2.0 / T ) , 2 ) - 1.0
        onsager[index] = alpha * ( 2.0 * ellipk( k0 ) - 2.0 * ellipe( k0 ) - ( 1 - k1 ) * ( ( np.pi / 2.0 ) + k1 * ellipk( k0 ) ) )

    mag_cm = np.zeros( len( vectorCurvaX ) )
    for index in range( len( mag_cm ) ):
       m = 0.5
       T_iter = vectorCurvaX[index]
       for iter in range( 1000 ):
           m = np.tanh( 4 * m / T_iter )
       mag_cm[index] = m

    campo_medio = np.zeros( len( vectorCurvaX ) )
    for index in range( len( campo_medio ) ):
        T = vectorCurvaX[index]
        m = mag_cm[index]
        R = 4.0 * np.power( T , -2 ) * np.power( np.cosh( 4.0 * m / T ), -2 )
        dm_dT = np.divide( -m * R, 1 - T * R )
        campo_medio[index] = -4.0 * m * dm_dT

    fig, ax = plt.subplots( dpi = img_dpi )
    DPI = fig.get_dpi()
    fig.set_size_inches( h / float( DPI ), w / float( DPI ) )
    ax.cla()
    ax.set_title( titulo )
    ax.set_ylim( [0, max( vectorCurvaY[:int( len( vectorCurvaY ) * 0.8 )] ) * 1.5 ] )
    ax.set_xlim( minEjeX, maxEjeX )
    ax.set_xlabel( leyendaEjeX )
    ax.set_ylabel( leyendaEjeY )
    ax.minorticks_off()
    ax.grid( True )

    ax.plot( vectorCurvaX , vectorCurvaY , 'D', color = "lightblue" , label = 'Metropolis' )
    ax.plot( vectorCurvaX , onsager , '-', color = "red" , label = 'Onsager' )
    ax.plot( vectorCurvaX , campo_medio , '--', color = "blue" , label = 'Campo Medio' )
    ax.legend( loc = 'upper right' )

    plt.savefig( leyendaEjeY, fmt = "png" )
    fig.canvas.draw()

def plotear_chi( N, data, w, h, img_dpi ):

    # Datos del plot
    vectorCurvaX = data[:, 0]
    vectorCurvaY = data[:, 1]
    vectorCurvaY = np.divide( vectorCurvaY, N );

    # Parametros del pLot
    minEjeX = min( vectorCurvaX )
    maxEjeX = max( vectorCurvaX )

    # Configuracion del plot
    configuracionPlot = "o-"

    # Texto del grafico
    titulo = "Susceptibilidad Magnetica vs T"
    leyendaEjeX = 'Temperatura'
    leyendaEjeY = 'Susceptibilidad Magnetica'

    onsager_m = np.power ( np.abs( np.ones( len( vectorCurvaX ) ) - np.power( np.sinh( np.divide( 2.0 , vectorCurvaX ) ) , -4.0 ) ) , 1.0 / 8.0 )

    onsager = np.zeros( len( vectorCurvaX ) )
    for index in range( len( onsager ) ):
        T = vectorCurvaX[index]
        m_0 = np.real( onsager_m[index] )
        onsager[index] = np.power( np.power( T, 2 ) * np.power( m_0, 7 ) * np.power( np.sinh( 2.0 / T ), 5 ) , -1 ) * np.cosh( 2.0 / T )

    mag_cm = np.zeros( len( vectorCurvaX ) )
    for index in range( len( mag_cm ) ):
       m = 0.5
       T_iter = vectorCurvaX[index]
       for iter in range( 1000 ):
           m = np.tanh( 4 * m / T_iter )
       mag_cm[index] = m

    campo_medio = np.zeros( len( vectorCurvaX ) )
    for index in range( len( campo_medio ) ):
        T = vectorCurvaX[index]
        m = mag_cm[index]
        R = np.power( T , -1 ) * np.power( np.cosh( 4.0 * m / T ), -2 )
        campo_medio[index] = np.divide( R, 1 - 4.0 * R )

    fig, ax = plt.subplots( dpi = img_dpi )
    DPI = fig.get_dpi()
    fig.set_size_inches( h / float( DPI ), w / float( DPI ) )
    ax.cla()
    ax.set_title( titulo )
    ax.set_ylim( [0, max( vectorCurvaY[:int( len( vectorCurvaY ) * 0.8 )] ) * 1.5 ] )
    ax.set_xlim( minEjeX, maxEjeX )
    ax.set_xlabel( leyendaEjeX )
    ax.set_ylabel( leyendaEjeY )
    ax.minorticks_off()
    ax.grid( True )

    ax.plot( vectorCurvaX , vectorCurvaY , 'D', color = "lightblue" , label = 'Metropolis' )
    ax.plot( vectorCurvaX , onsager , '-', color = "red" , label = 'Onsager' )
    ax.plot( vectorCurvaX , campo_medio , '--', color = "blue" , label = 'Campo Medio' )
    ax.legend( loc = 'upper right' )

    plt.savefig( leyendaEjeY, fmt = "png" )
    fig.canvas.draw()

def plotear_C_r( N, dT, data, w, h, img_dpi ):

    # Datos del plot
    vectorCurvaX = np.arange( L )
    matriz_C_r = data[:, 1: ]
    vectorT = data[:, 0]
    num_puntos = len( vectorT )
    num_plotss = int( ( np.max( vectorT ) / dT ) )

    # Parametros del pLot
    minEjeX = min( vectorCurvaX )
    maxEjeX = max( vectorCurvaX )

    # Configuracion del plot
    configuracionPlot = "o-"

    # Texto del grafico
    titulo = "Correlacion vs r en funcion de T"
    leyendaEjeX = "Radio r"
    leyendaEjeY = "Correlacion"

    corr_T = np.zeros( [ L  , num_plotss ] )
    corr_labels = []
    n = 0
    ind_T = 0
    for T in vectorT :
        if( T % dT ) == 0 :
            fila = matriz_C_r[ind_T, : ]
            columna = np.transpose( fila )
            corr_T[:, n] = columna
            corr_labels.append( 'Metropolis - {0} K'.format( T ) )
            n += 1
        else:
            pass
        ind_T += 1

    mag_cm = np.zeros( len( vectorT ) )
    for index in range( len( mag_cm ) ):
       m = 0.5
       T_iter = vectorT[index]
       for iter in range( 1000 ):
           m = np.tanh( 4 * m / T_iter )
       mag_cm[index] = m

    corr_cm_T = np.zeros( [L  , num_plotss ] )
    n = 0
    ind_T = 0
    for T in vectorT :
        if( T % dT ) == 0 :
            for r in vectorCurvaX :
                corr_cm_T[r, n] = np.power( np.tanh( 4.0 * mag_cm[ind_T] / T ) , r )
            corr_labels.append( 'Campo Medio - {0} K'.format( T ) )
            n += 1
        else:
            pass
        ind_T += 1

    fig, ax = plt.subplots( dpi = img_dpi )
    DPI = fig.get_dpi()
    fig.set_size_inches( h / float( DPI ), w / float( DPI ) )
    ax.cla()
    ax.set_title( titulo )
    ax.set_xlim( minEjeX, maxEjeX )
    ax.set_ylim( -0.8 , 1.1 )
    ax.set_xlabel( leyendaEjeX )
    ax.set_ylabel( leyendaEjeY )
    ax.minorticks_off()
    ax.grid( True )

    ax.plot( vectorCurvaX , corr_T , 'd-' )
    ax.plot( vectorCurvaX , corr_cm_T , '--' )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend( corr_labels, ncol = 2, loc = 'lower left', fontsize = 11 )

    plt.savefig( leyendaEjeY, fmt = "png" )
    fig.canvas.draw()

def plotear_long_corr( N, data, w, h, img_dpi ):

    # Datos del plot
    vectorCurvaX = data[:, 0]
    vectorCurvaY = data[:, 1]
    num_puntos = len( vectorCurvaX )

    # Parametros del pLot
    minEjeX = min( vectorCurvaX )
    maxEjeX = max( vectorCurvaX )

    # Configuracion del plot
    configuracionPlot = "o-"

    # Texto del grafico
    titulo = "Longitud de Correlacion vs T"
    leyendaEjeX = 'Temperatura'
    leyendaEjeY = 'Longitud de Correlacion'


    fig, ax = plt.subplots( dpi = img_dpi )
    DPI = fig.get_dpi()
    fig.set_size_inches( h / float( DPI ), w / float( DPI ) )
    ax.cla()
    ax.set_title( titulo )
    ax.set_xlim( minEjeX, maxEjeX )
    ax.set_xlabel( leyendaEjeX )
    ax.set_ylabel( leyendaEjeY )
    ax.minorticks_off()
    ax.grid( True )

    ax.plot( vectorCurvaX , vectorCurvaY , 'D', color = "lightblue" , label = 'Metropolis' )
    ax.legend( loc = 'upper right' )

    plt.savefig( leyendaEjeY , fmt = "png" )
    fig.canvas.draw()

if __name__ == '__main__':
    L = 16
    N = np.float( L * L )

    dT = 1.0

    # Leo los datos del disco
    magnet = np.loadtxt( 'magnetizacion.txt' )
    var_magnet = np.loadtxt( 'varianza_magnetizacion.txt' )
    energia = np.loadtxt( 'energia.txt' )
    var_energia = np.loadtxt( 'varianza_energia.txt' )
    cv = np.loadtxt( 'capacidad_calorifica.txt' )
    chi = np.loadtxt( 'susceptibilidad_magnetica.txt' )
    long_corr = np.loadtxt( 'longitud_de_correlacion.txt' )
    C_r = np.loadtxt( 'matriz_de_correlacion.txt' )

    print( "Carga finalizada" )

    w = 768.0
    h = 1024.0
    img_dpi = 80


    # Ploteo las curvas
    plotear_mag( N, magnet, w, h, img_dpi )
    plotear_energia( N, energia, w, h, img_dpi )
    plotear_cv( N , cv, w, h, img_dpi )
    plotear_chi( N , chi, w, h, img_dpi )
    plotear_C_r( L, dT, C_r, w, h, img_dpi )
    plotear_long_corr( N, long_corr, w, h, img_dpi )

    # plt.show( block = True )


    pass
