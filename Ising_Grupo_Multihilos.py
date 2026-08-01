
# -*- coding: utf-8 -*-
# script general para hacer una corrida a un set de parámetros,
# beta, tamaño de la red,

from copy import deepcopy
from multiprocessing.pool import ThreadPool
import ctypes
import hashlib
import multiprocessing
import os
import time

from mpl_toolkits.mplot3d.axes3d import Axes3D
from scipy.integrate.quadrature import simps

import matplotlib.pyplot as plt
import numpy as np

############################################################################
# Setup del sistema para habilitar la afinidad al core - solo valido en Linux - desactivar en windows
# os.system( "taskset -p 0xff %d" % os.getpid() )  # esto solo anda para Linux en Windows con Python 3.3 descomentar la linea con sched_setaffinity en ising-worker


# atplotlib.use('tkAgg')
plt.ion()

#############################################################
# estas son las funciones que tienen que escribir ustedes
# calcMagnet ya esta lista

def calcMagnet( S ):

    M = np.sum( S )

    return  M

def calcEnergia( S ):
    # Calcula la Energia de una configuracion dada
    energia = 0
    N = len( S )
    for i in range( len( S ) ):
        for j in range( len( S ) ):
            s = S[i, j]
            nv = S[( i + 1 ) % N, j] + S[i, ( j + 1 ) % N]
            energia += -( nv * s )
    return  energia

#############################################################
# Implementacion del algoritmo de Metropolis - Monte Carlo

def ising2Dpaso( S, beta ):
    N = len( S );
    curr_proc_pid = multiprocessing.current_process().ident
    curr_seed = int( time.time() * 100 * curr_proc_pid )
    seed = np.frombuffer( hashlib.sha256( str( curr_seed ) ).digest() , dtype = 'int32' )
    generator = np.random.RandomState( seed )

    i = generator.randint( 0, N )
    j = generator.randint( 0, N )
    s = S[i, j]
    s_ini = s
    nvEnergia = S[( i + 1 ) % N, j] + S[i, ( j + 1 ) % N] + S[( i - 1 ) % N, j] + S[i, ( j - 1 ) % N]
    energiaInicial = -s_ini * nvEnergia  # Hay un -1 delante del J
    dE = -2 * energiaInicial  # Converge al estado de Minima Energia, Todos los Espines Apuntan en la misma direccion a baja T

    if dE <= 0:
        s *= -1

    else:
        variableAleatoria = generator.uniform( 0.0, 1.0 )  # random()
        exponencial = np.exp( -dE * beta )
        azar = variableAleatoria < exponencial

        if azar :
            s *= -1
    S[i, j] = s
    energiaFinal = -s * nvEnergia  # J = -1
    dM = s - s_ini
    dE = energiaFinal - energiaInicial

    return dE , dM


# Este metodo inicializa:
# + La matriz de los espines S
# + el delta de energia dE,
# + el dellta de magnetizacion dM,
# + el numero de pasos de pretermalizacion npre,
# + la energia inicial energia.
# + la magnetizacion inicial magnet,
# + el numero de pasos totales npasos


# Este metodo inicializa:
# + La matriz de los espines S
# + el delta de energia dE,
# + el dellta de magnetizacion dM,
# + el numero de pasos de pretermalizacion npre,
# + la energia inicial energia.
# + la magnetizacion inicial magnet,
# + el numero de pasos totales npasos

def inicializarEstado( L, T, beta, N , esCorridaInicial ) :

    # me genero arrays vacios, a ser llenados con los valores de energia y magnetizacion
    energia = np.zeros( npasos )
    magnet = np.zeros( npasos )

    S = None

    # propongo un estado inicial al azar
    # S es una matriz de 1 y -1 indicando las dos proyecciones de
    # espin
    # Solo se ejecuta si es la primera iteracion
    if esCorridaInicial:
        # Genera el seed para la matriz random muntiplicando el tiempo por el valor del pid y obteniendo un hash a partir de eso
        # Es necesario hacer esto porque sino los multiples procesos generan todos el mismo estado inicial en vez de generar
        # estados iniciales diferentes
        curr_proc_pid = multiprocessing.current_process().ident
        curr_seed = int( time.time() * 100 * curr_proc_pid )
        seed = np.frombuffer( hashlib.sha256( str( curr_seed ) ).digest() , dtype = 'int32' )
        generator = np.random.RandomState( seed )

        S = 2 * ( generator.rand( L, L ) > 0.5 ) - 1;

        # pretermalizo
        # ising2Dpaso hace un nuevo elemento de la cadena de Markov
        # la tienen que escribir Uds...
        for n in range( npre ):
            dE, dM = ising2Dpaso( S, beta )

    return S, npre, npasos, energia, magnet

# Parte del script original, setea el pyplot para visualizar la matriz de espines
def inicializarVisualizacion( matriz_S ):

    # muestro el estado inicial
    fig_s, ax_s = plt.subplots()
    ax_s.imshow( matriz_S, interpolation = 'none' )

    # Preparamos la figura para graficar
    # la energia y magnetizacion
    fig_em, ax_list = plt.subplots( 2, 1 )

    return fig_s, ax_s, fig_em, ax_list

# Parte del script original, setea el pyplot para visualizar la curva de magnetizacion y de energia
def visualizarDatos ( fig_s, ax_s, fig_em, ax_list, npasos, n, beta, matriz_S, magnet, energia ):

    ax_e, ax_m = ax_list

    # grafico el estado de la red, actualizandolo en cada iteracion
    ax_s.clear()
    ax_s.set_title( "n=%i beta=%.2f mag=%.2f energia=%.2f" % ( n, beta, magnet[n], energia[n] ) )
    ax_s.imshow( matriz_S, interpolation = 'none' )
    fig_s.canvas.draw()
    plt.pause( 0.01 )

    # graficamos la energia y magnetizacion hasta el paso actual
    ax_e.cla()
    ax_m.cla()

    ax_e.set_xlim( 0, npasos )
    ax_e.set_ylabel( 'Energia' )
    ax_m.set_xlim( 0, npasos )
    ax_m.set_ylabel( 'Magnetizacion' )

    ax_e.grid( True )
    ax_m.grid( True )

    ax_e.plot( energia[:n] )
    ax_m.plot( magnet[:n] )
    fig_em.canvas.draw()
    plt.pause( 0.01 )

    pass

#########################################################
# La funcion que realiza la simulacion propiamente dicha
def calcularIsing( L, T, npre, npasos, SInicial = None ):
    beta = 1 / T
    N = L * L

    # Determina si es la primera iteracion, correspondiente a T = +infinito
    # Esto se hace determinando si la matriz Inicial pasada como parametro es None o no
    esCorridaInicial = not ( type( SInicial ) == np.ndarray )

    # Funcion de inicializacion de los parametros, solo efecuta la pretermalizacion si esCorridaInicial = True
    S, npre, npasos, energia, magnet = inicializarEstado( L, T, beta, N, esCorridaInicial )

    # Si no es la primera iteracion, la matriz a partir de la cual hay que iterar es la de la iteracion anterior
    if not esCorridaInicial:
        S = SInicial

    matricesSpin = [np.array( S )]

    energia[0] = calcEnergia( S )
    magnet[0] = calcMagnet( S )

    # Para cada iteracion entre n = 0 y n = npasos-1 aplico el algoritmo y registro los datos
    for n in range( 0, npasos - 1, 1 ):
        actualizarValores( S, n, beta, energia, magnet )
        copyS = np.array( S )
        matricesSpin.append( copyS )

        # print( "iteracion numero : %i" % n )

    vectorResultados = list( [energia, magnet] )

    return matricesSpin, vectorResultados, beta

# Parte del script original, es para ver la evolucion de la matriz de espines en la simulacion
def visualizarIsing( pasoMuestreo, matricesSpin, vectorResultados, beta, npre, npasos ):

    _S = matricesSpin[0]
    fig_s, ax_s, fig_em, ax_list = inicializarVisualizacion( _S )
    fig_s.show()

    energia = vectorResultados[0]
    magnet = vectorResultados[1]

    for n in range( 1, len( matricesSpin ) - 1, 1 ) :
        if n % pasoMuestreo == 0:
            _S = matricesSpin[n]
            visualizarDatos ( fig_s, ax_s, fig_em, ax_list, npasos, n, beta, _S, magnet, energia )

    pass

# Funcion generica  para plotear las diferentes curvas
# esta pensada para obtener graficos preliminares y ver que todo este bien
# para generar los graficos finales usar graficos_entrega.py
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

# La matriz de Correlacion en funcion de la temperatura se visualiza diferente
def plotearCorrelacion( titulo, leyendaEjeX, leyendaEjeY, leyendaEjeZ, minEjeX, maxEjeX , vectorCurvaX , vectorCurvaY , matrizZ ):
    fig, ax = plt.subplots()
    ax.cla()
    ax.set_title( titulo )
    ax.set_ylabel( leyendaEjeY )
    ax.set_xlabel( leyendaEjeX )
    ax.minorticks_on
    ax.grid( True )
    ax.imshow( matrizZ, interpolation = 'none' )
    plt.locator_params( axis = 'y', nbins = len( vectorCurvaY ) + 1 )
    plt.locator_params( axis = 'x', nbins = len( vectorCurvaX ) + 1 )
    fig.canvas.draw()

    ax.set_xticklabels( vectorCurvaX )
    ax.set_yticklabels( vectorCurvaY )
    ax.tick_params( labelsize = 8 )
    fig.canvas.draw()
    plt.savefig( leyendaEjeZ, fmt = "png" )

# Actualiza los valores de energia y magnetizacion en cada iteracion de la simulacion
def actualizarValores( S, n, beta, energia, magnet ):

    dE, dM = ising2Dpaso( S, beta )
    energia[n + 1] = energia[n] + np.float( dE )
    magnet[n + 1] = magnet[n] + np.float( dM )

    pass


def calcular_long_corr( L, iteracionesT, matriz_C_r ):
    # Datos de la matriz de Correlacion

    long_corr = np.zeros( iteracionesT )
    for indice_T in range( iteracionesT ) :
        indices = len( matriz_C_r[indice_T, :] )
        rango = 3 * indices / 4  # Tomo submuestras entre L/4 y L
        temp_l = np.zeros( rango )
        rmax = indices - 1

        # Metodo de la integral usando micro-muestreo
        for r in range( 1, indices ):

            # Calculo el factor de escala de la submuestra
            c_lambda = np.float( r ) / np.float( L )
            expo = np.divide( 1.0, c_lambda )

            # Valore de la funcion de correlacion en r= 0 y r = rmax
            C_0 = matriz_C_r[indice_T, 0]
            C_r = matriz_C_r[indice_T, r]

            # prefactor que multiplica a la integral
            delta = np.float( r ) / ( np.float( L ) * ( C_0 - C_r ) )

            # Sumas superior en inferior de la integral
            # Nota importante, como los datos de la correlacion tienen fluctuaciones por debajo de cero
            # hay que calcular la raiz de C(r) como numero complejo y despues tomar la parte real
            S1 = np.real( np.sum( ( matriz_C_r[indice_T, 0:( rmax - 1 ) ] + 0j ) ** expo ) )
            S2 = np.real( np.sum( ( matriz_C_r[indice_T, 1:rmax ] + 0j ) ** expo ) )

            # Almaceno el valor de la longitud de correlacion para la submuestra
            temp_l[r - 1 - ( rango )] = np.divide ( S1 + S2 , ( 2.0 * delta ) )

            if not( np.isfinite( temp_l[r - ( rango )] ) ):
                temp_l[r - ( rango )] = 0.0
            else:
                pass

        # Calculo la longitud de correlacion a la temperatura dada como promedio de todas las submuestras
        idx_tmp = len( temp_l ) - 1
        res = np.divide( np.sum( temp_l[0: idx_tmp - 1] ) + np.sum( temp_l[1: idx_tmp] ) , 2 * ( idx_tmp + 1 ) )

        long_corr[indice_T] = res

    return long_corr


if __name__ == '__main__':

    # Aca defino los parámetros y corro la cadena de markov
    # Lado de la red,
    L = 16

#    Parte del scipt original, no son necesarios a menos que se quiera visualizar la simulacion
#    lo cual no es recomendable hacer siendo que ahora el script se ejecuta en paralelo
#    T = 1.25
#    pasoMuestreo = 1000

    # defino la cantidad de iteraciones de cada etapa
    npre = 150000
    npasos = 30000

    # Asumo que estoy partiendo de Temperatura +infinito
    SInicial = None

    # # Aca comienza la rutina de visualizacion de los datos almacenados en matricesSpin y vectorResultados
    # print( "\nIniciando visualizacion\n" )
    # visualizarIsing( pasoMuestreo, matricesSpin, vectorResultados, beta, npre,npasos )

    # Defino los parametros de la simulacion
    Tmin = 0.01
    Tmax = 5.00
    Tpaso = 0.01

    # Numero de corridas a ejecutar entre Tmax y Tmin
    numeroCorridas = 30
    # Numero de iteraciones de termalizacion a promediar
    valoresPromediadios = 20000
    inicioPromedio = npasos - valoresPromediadios - 1
    iteracionesT = int( ( ( Tmax - Tmin ) / Tpaso ) + 1 )
    # Numero de procesos por corrida a ejecutar, se define un numero maximo para no aniquilar la maquina
    maxBatch = 10
    numProc = min( maxBatch, numeroCorridas )

    # Defino que fila de espines voy a samplear para calcular la correlaciom , por defecto, tomo la fila a la mitad de la muestra
    # Tomo que r0 = 0 y rmax = L
    filaEspines = int( L / 2 )  # toma como fila de espines para calcular la correlacion la posicion a la mitad de la matriz

    ejeT = np.linspace( Tmax, Tmin, iteracionesT )

    # Creo Las array de Memoria Compartida
    #  Aca van a ir a parar todos los datos de todas las simulaciones que se corren en paralelo
    # Es decir, cada vez que se llega a una termalizacion a un T dada, cada proceso en paralelo
    # guarda los datos en una zona comun de memoria
    shared_L = multiprocessing.Value( ctypes.c_int , L )
    shared_iteracionesT = multiprocessing.Value( ctypes.c_int , iteracionesT )
    shared_valoresPromediadios = multiprocessing.Value( ctypes.c_int, valoresPromediadios )
    shared_npre = multiprocessing.Value( ctypes.c_int , npre )
    shared_npasos = multiprocessing.Value( ctypes.c_int , npasos )
    shared_filaEspines = multiprocessing.Value( ctypes.c_int , filaEspines )

    shared_valoresT_base = multiprocessing.RawArray( ctypes.c_longdouble, iteracionesT * numeroCorridas )
    shared_valoresT = np.frombuffer( shared_valoresT_base, dtype = np.dtype( np.longdouble ) )

    shared_valoresMagnet_base = multiprocessing.RawArray( ctypes.c_longdouble, iteracionesT * numeroCorridas )
    shared_valoresMagnet = np.frombuffer( shared_valoresMagnet_base, dtype = np.dtype( np.longdouble ) )

    shared_valoresEnergia_base = multiprocessing.RawArray( ctypes.c_longdouble, iteracionesT * numeroCorridas )
    shared_valoresEnergia = np.frombuffer( shared_valoresEnergia_base, dtype = np.dtype( np.longdouble ) )

    shared_arrayEspinesCorr_base = multiprocessing.RawArray( ctypes.c_float, 2 * iteracionesT * numeroCorridas * L )
    shared_arrayEspinesCorr = np.frombuffer( shared_arrayEspinesCorr_base, dtype = np.dtype( np.float ) )

    shared_arrayParesCorr_base = multiprocessing.RawArray( ctypes.c_float, 2 * iteracionesT * numeroCorridas * L )
    shared_arrayParesCorr = np.frombuffer( shared_arrayParesCorr_base, dtype = np.dtype( np.float ) )

    # Un "worker" es una seccion de codigo que se ejecuta en paralelo
    # En particular esta es la ejecuta toda la corrida entre Tmax y Tmin
    def ising_worker( corrida, def_params = [shared_filaEspines, shared_npre, shared_npasos, shared_valoresPromediadios, shared_iteracionesT, shared_L, shared_valoresT, shared_valoresMagnet, shared_valoresEnergia, shared_arrayEspinesCorr, shared_arrayParesCorr] ):

        # #Setup de la afinidad del proceso y ajuste de la propidad de este - solo valido para Linux - desactivar en Windows
#        curr_proc = multiprocessing.current_process()
#        pid = curr_proc.pid
#        cpu_count = multiprocessing.cpu_count()
#        affinity = corrida % cpu_count
#        os.system( "taskset -p {0} {1}".format( cpu_count , pid ) )


        # Setup de la afinidad del proceso - solo valido en Windows con Python 3.3 o superior - desactivar en Linux
#        import win32api, win32process, win32con
#        pid = win32api.GetCurrentProcessId()
#        handle = win32api.OpenProcess( win32con.PROCESS_ALL_ACCESS, True, pid )
#        win32process.SetPriorityClass( handle, win32process.HIGH_PRIORITY_CLASS )
#        cores_list = os.sched_getaffinity(0) #Esto solo anda en Python 3.3 o superior
#        os.sched_setaffinity({ corrida % len(cores_list) } , cores_list  ) #Esto solo anda en Python 3.3 o superior

        # # Setea el nivel de prioridad del proceso a "Alto" , para matar la maquina poner os.nice(-20)
        # # En Linux se necesitan nivel de administrador su o sudo
        os. nice( -10 )

        # #Valores iniciales de la parametros
        SInicial = None
        indiceT = 0

        # #defino los workers para los ser ejecutados en los hilos
        # # Estos workers se ejecutan de manera asincronica, es decir, se manda a calcular y paso a la siguiente iteracion
        # # como estoy usando memoria compartida no hay riesgo de pisar los datos de otro proceso
        # # por lo que no espero a termine el calculo antes de pasar a la siguiente iteracion
        def temperatura_worker( shared_valoresT, indice_actual, T ):
            shared_valoresT[indice_actual] = T

        def magnet_worker( vectorResultados, shared_valoresMagnet, shared_valoresPromediadios ):
            magnet = deepcopy( vectorResultados[1] )
            # shared_valoresMagnet[indice_actual] = np.sum( magnet[inicioPromedio :-1] ) / shared_valoresPromediadios.value  # obtengo el promedio de la magnetizacion
            shared_valoresMagnet[indice_actual] = magnet[inicioPromedio :].mean()

        def energia_worker( vectorResultados, shared_valoresEnergia, shared_valoresPromediadios ):
            energia = deepcopy( vectorResultados[0] )
            # shared_valoresEnergia[indice_actual] = np.sum( energia[inicioPromedio :-1] ) / shared_valoresPromediadios.value  # obtengo el promedio de la energia
            shared_valoresEnergia[indice_actual] = energia[inicioPromedio :].mean()

        def corrEspines_worker( SInicial, shared_arrayEspinesCorr, shared_L, indice_actual , shared_filaEspines ):
            # Ahora obtengo la fila de espines a partir de la cual voy a calcular la corellacion para el espin (1,L/2)
            fila = deepcopy( SInicial[shared_filaEspines.value, : ] )
            for indice_espin in range( shared_L.value ):
                shared_arrayEspinesCorr[indice_actual * shared_L.value + indice_espin] = fila[indice_espin]

        def paresEspines_worker( SInicial, shared_arrayEspinesCorr, shared_L, indice_actual , shared_filaEspines ):
            s0 = deepcopy( SInicial[shared_filaEspines.value, 0] )
            for indice_par_espin in range( shared_L.value ):
                shared_arrayParesCorr[indice_actual * shared_L.value + indice_par_espin] = s0 * deepcopy( SInicial[shared_filaEspines.value, indice_par_espin] )

        # #Setup del pool de hilos
        threadpool = ThreadPool( 20 )

        for T in ejeT :
            # Este el indice de la memoria compartida donde voy a poner los datos
            indice_actual = corrida * shared_iteracionesT.value + indiceT
            print( "Iteracion para %i de T : %f" % ( corrida + 1, T ) )
            matricesSpin, vectorResultados, beta = calcularIsing( shared_L.value, T , shared_npre.value , shared_npasos.value, SInicial )  # calculo la matriz de Ising
            SInicial = deepcopy( matricesSpin[-1] )


            # Defino los hilos y los hago correr
            threadpool.apply_async( temperatura_worker, args = [ shared_valoresT, indice_actual, T ] )
            threadpool.apply_async( magnet_worker, args = [  vectorResultados, shared_valoresMagnet, shared_valoresPromediadios ] )
            threadpool.apply_async( energia_worker, args = [ vectorResultados, shared_valoresEnergia, shared_valoresPromediadios ] )
            threadpool.apply_async( corrEspines_worker, args = [ SInicial, shared_arrayEspinesCorr, shared_L, indice_actual , shared_filaEspines ] )
            threadpool.apply_async( paresEspines_worker, args = [ SInicial, shared_arrayEspinesCorr, shared_L, indice_actual , shared_filaEspines ] )

            indiceT += 1

        # Cierro los el pool de hilos
        threadpool.close()
        threadpool.join()

    # Para cada corrida se calculan la magnetizacion y la energia promedio para una temperatura dada entre Tmax y Tmin
    # El calculo se hace de la temperatura mas alta a la mas baja, tomando como matriz de espines inicial a una temperatura dada
    # la ultima matriz espines termalizada de la temperatura anterior ( T_anterior > T )
    # aca va el el proceso con multihilos

    pool = multiprocessing.Pool( processes = numProc )
    pool.map( ising_worker, range( numeroCorridas ) )

    print( "calculo finalizado" )

    # # Expreso el arreglo con las tiras de espines como tira de espines para un T dado en funcion de las corridas
    # # Rompemos el vector en matrices M[Temperatura,Parametro]
    matrizMagnet = deepcopy( shared_valoresMagnet.reshape( ( numeroCorridas, iteracionesT ) ) )
    matrizEnergia = deepcopy ( shared_valoresEnergia.reshape( ( numeroCorridas, iteracionesT ) ) )

    valoresT = deepcopy ( shared_valoresT.reshape( ( 1, iteracionesT * numeroCorridas ) ) )
    matrizEspinesCorr = deepcopy( shared_arrayEspinesCorr.reshape( ( numeroCorridas, iteracionesT * L ) ) )
    paresEspinesCorr = deepcopy( shared_arrayParesCorr.reshape( ( numeroCorridas, iteracionesT * L ) ) )

    # Ahora calculo el promedio sobre todas las corridas
    # #matrizEspinesCorr = np.divide( np.sum ( matrizEspinesCorr, axis = 0 ) , numeroCorridas )
    matrizEspinesCorr = np.average( matrizEspinesCorr, axis = 0 )

    # idem para los pares de espines
    # #paresEspinesCorr = np.divide( np.sum( paresEspinesCorr, axis = 0 ) , numeroCorridas )
    paresEspinesCorr = np.average( paresEspinesCorr, axis = 0 )

    # Y separo la matriz en una mariz de filas de espines en funcion de la temperatura
    matrizEspinesCorr = matrizEspinesCorr.reshape( iteracionesT, L )

    # idem para la matriz de pares de espines
    paresEspinesCorr = paresEspinesCorr.reshape( iteracionesT, L )


    # Ahora calculo la matriz de C(r) en funcion de T para r entro 0 y L
    matrizC_r = np.zeros( ( iteracionesT, L ) )
    for indice_T in range( iteracionesT ):
        s0 = matrizEspinesCorr[indice_T, 0]  # Para la temperatura dada, tomo el valor medio del primer spin de la fila
        for indice_r in range( L ) :
            sr = matrizEspinesCorr[indice_T, indice_r]  # tomo valor medio del espin en la posision r
            s0r = paresEspinesCorr[indice_T, indice_r]  # tomo el valor medio del par s0*sr
            matrizC_r[ indice_T, indice_r ] = s0r - s0 * sr  # C(R) = <s_0>*<s_r> - <s_0*s_r>

    # Ahora calculo la longitud de correlacion por metodo de la integral con micro-muestreo
    longCorr = np.zeros( iteracionesT )
    longCorr = calcular_long_corr( L, iteracionesT, matrizC_r )

    ejeTemperatura = np.array( valoresT )

    ###########################################################################################################
    # De todas las Corridas, calculamos el valor medio y el desvio standard de la magnetizacion y de la Energia

    # #valorMedioMagnet = np.sum( np.abs( matrizMagnet ) , axis = 0 ) / numeroCorridas  # Calcula <M> a partir de la matriz
    # #valorCuadradoMagnet = np.sum( np.power( matrizMagnet, 2 ) , axis = 0 ) / numeroCorridas  # Calcula <M^2> a partir de la matriz
    # #sigmaCuadradoMagnet = valorCuadradoMagnet - np.power( valorMedioMagnet, 2 )  # Calcula sigma cuadrado de M

    valorMedioMagnet = np.average( np.abs( matrizMagnet ) , axis = 0 )
    sigmaCuadradoMagnet = np.var( np.abs( matrizMagnet ) , axis = 0 )

    # #valorMedioEnergia = np.sum( matrizEnergia  , axis = 0 ) / numeroCorridas  # Calcula <E> a partir de la matriz
    # #valorCuadradoEnergia = np.sum( np.power( matrizEnergia, 2 ) , axis = 0 ) / numeroCorridas  # Calcula <E^2> a partir de la matriz
    # #sigmaCuadradoEnergia = valorCuadradoEnergia - np.power( valorMedioEnergia, 2 )  # Calcula sigma de E

    valorMedioEnergia = np.average( matrizEnergia  , axis = 0 )
    sigmaCuadradoEnergia = np.var( matrizEnergia  , axis = 0 )

    # Calculamos la Capacidad Calorifica como la derivada de la densidad de energia en funcion de la temperatura
    curvaCv = np.divide( sigmaCuadradoEnergia , np.power( np.linspace( Tmax, Tmin, iteracionesT ) , 2 ) )
    curvaChi = np.divide( sigmaCuadradoMagnet , np.linspace( Tmax, Tmin, iteracionesT ) )

    # plt.close()

    ejeXNormal = np.linspace( Tmax, Tmin, iteracionesT )
    ejePosiciones = range( 0, L + 1 )
    list.reverse( ejePosiciones )
    for idx in range( len( ejePosiciones ) ) : ejePosiciones[idx] = str( ejePosiciones[idx] )
    ejeTemperatura = np.linspace( Tmax, Tmin, iteracionesT + 1 )
    for idx in range( len( ejeTemperatura ) ) : ejeTemperatura[idx] = str( ejeTemperatura[idx] )[0:3]

    # Grabo los datos a disco
    np.savetxt( 'magnetizacion.txt', np.column_stack( ( ejeXNormal, valorMedioMagnet ) ) , delimiter = '\t' )
    np.savetxt( 'varianza_magnetizacion.txt', np.column_stack( ( ejeXNormal, sigmaCuadradoMagnet ) ) , delimiter = '\t' )
    np.savetxt( 'energia.txt', np.column_stack( ( ejeXNormal, valorMedioEnergia ) ) , delimiter = '\t' )
    np.savetxt( 'varianza_energia.txt', np.column_stack( ( ejeXNormal, sigmaCuadradoEnergia ) ) , delimiter = '\t' )
    np.savetxt( 'capacidad_calorifica.txt', np.column_stack( ( ejeXNormal, curvaCv ) ) , delimiter = '\t' )
    np.savetxt( 'susceptibilidad_magnetica.txt', np.column_stack( ( ejeXNormal, curvaChi ) ) , delimiter = '\t' )
    np.savetxt( 'longitud_de_correlacion.txt', np.column_stack( ( ejeXNormal, longCorr ) ) , delimiter = '\t' )
    np.savetxt( 'matriz_de_correlacion.txt', np.column_stack( ( ejeXNormal, matrizC_r ) ) , delimiter = '\t' )

    # Ploteo las curvas
    plotearCurvas( "Magnetizacion vs T", 'Temperatura', 'Magnetizacion', Tmin, Tmax , ejeXNormal  , valorMedioMagnet, "o-" )
    plotearCurvas( "Varianza de la Magnetizacion vs T", 'Temperatura', 'Varianza de la Magnetizacion', Tmin, Tmax , ejeXNormal , sigmaCuadradoMagnet, "o-" )
    plotearCurvas( "Energia vs T", 'Temperatura', 'Energia', Tmin, Tmax , ejeXNormal , valorMedioEnergia, "o-" )
    plotearCurvas( "Varianza de la Energia vs T", 'Temperatura', 'Varianza de la Energia', Tmin, Tmax , ejeXNormal ,
sigmaCuadradoEnergia, "o-" )
    plotearCurvas( "Capacidad Calorifica vs T", 'Temperatura', 'Capacidad Calorifica', Tmin, Tmax, ejeXNormal, curvaCv, "o-" )
    plotearCurvas( "Susceptibilidad Magnetica vs T", 'Temperatura', 'Susceptibilidad Magnetica', Tmin, Tmax, ejeXNormal, curvaChi, "o-" )
    plotearCurvas( "Longitud de Correlacion vs T", 'Temperatura', 'Longitud de Correlacion', Tmin, Tmax, ejeXNormal, longCorr, "o-" )
    plotearCorrelacion( "Matriz de Correlacion vs T", "Temperatura", "Sitio de Espin", "Correlacion", Tmin, Tmax , ejeTemperatura, ejePosiciones , np.transpose( matrizC_r ) )

    plt.show( block = True )

pass
