# Ising_2D_Multihilos

Simulaciones del modelo de Ising en 2D usando Monte Carlo (Metropolis) con opciones secuenciales y multihilo/multiproceso.

Este repositorio contiene implementaciones educativas y de investigación para estudiar la magnetización, energía, capacidad calorífica, susceptibilidad y longitud de correlación en redes 2D de espines.

Principales componentes

- Ising2D0bis.py — versión didáctica y simplificada. Incluye plantillas para implementar calcEnergia e ising2Dpaso y visualizar la evolución en tiempo real.
- Ising_Grupo_Multihilos.py — implementación más completa y orientada a corridas en paralelo (multiprocessing + ThreadPool). Calcula promedios sobre múltiples corridas, guarda resultados y genera datos para graficar.
- graficos_entrega.py — utilidades para leer los ficheros de salida y generar gráficos comparando resultados con soluciones analíticas (Onsager, campo medio, etc.).
- integrador.py — script independiente que integra un sistema 2D de ODEs y dibuja retratos de fase (ejemplo adicional, no parte del Ising directamente).
- LICENSE — licencia del proyecto.

Requisitos

- Python 3.6+ (recomendado Python 3.8+)
- numpy
- scipy
- matplotlib

Instalación rápida (virtualenv recomendado):

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install numpy scipy matplotlib
```

Cómo ejecutar

1) Versión didáctica (visualización interactiva)

```bash
python3 Ising2D0bis.py
```

- Edita las funciones `calcEnergia(S)` e `ising2Dpaso(S, beta)` dentro de Ising2D0bis.py para completar el algoritmo de Metropolis si quieres probar implementaciones propias.
- Este script muestra la red y las curvas de energía/magnetización en tiempo real (usa matplotlib en modo interactivo).

2) Versión para corridas (multihilos/multiproceso, orientada a cálculos intensivos)

```bash
python3 Ising_Grupo_Multihilos.py
```

- Archivo principal para realizar barridos de temperatura, ejecutar varias corridas en paralelo y guardar resultados.
- Parámetros principales a editar en el script: L (tamaño de la red), npre (pasos de pretermalización), npasos (pasos de muestreo), Tmin/Tmax/Tpaso y numeroCorridas.
- El script crea archivos de salida en el directorio de trabajo:
  - magnetizacion.txt
  - varianza_magnetizacion.txt
  - energia.txt
  - varianza_energia.txt
  - capacidad_calorifica.txt
  - susceptibilidad_magnetica.txt
  - longitud_de_correlacion.txt
  - matriz_de_correlacion.txt

3) Generar gráficos a partir de los ficheros de salida

```bash
python3 graficos_entrega.py
```

- `graficos_entrega.py` lee los archivos generados por `Ising_Grupo_Multihilos.py` y guarda gráficos en PNG en el directorio actual.

Notas y recomendaciones

- Afinidad/priority: El script multihilo contiene llamadas a `os.nice()` y líneas para configurar afinidad al CPU; esas secciónes están pensadas para Linux. Si ejecutas en Windows, comenta o adapta las llamadas (hay comentarios en el código).
- Memoria y CPU: Las corridas grandes (L grande, muchos npasos y numeroCorridas elevado) consumen mucha memoria/CPU. Ajusta `numProc` y `maxBatch` en `Ising_Grupo_Multihilos.py` para no saturar la máquina.
- Semillas aleatorias: Para evitar que procesos hijos generen la misma semilla se usa un hash del tiempo y PID; puedes cambiar el mecanismo si quieres reproducibilidad exacta (usa np.random.seed con valores fijos y pasa semillas reproducibles a cada proceso).
- Compatibilidad: Código probado en entornos Linux con Python 3.x. Algunas partes usan funciones y módulos que requieren SciPy y Matplotlib.

Desarrollo

- Si quieres implementar mejoras:
  - Añade un `requirements.txt` o `pyproject.toml` para fijar dependencias.
  - Extrae funciones comunes a un módulo `ising/` para facilitar tests unitarios.
  - Añade argparse a los scripts para controlar parámetros desde la línea de comandos en lugar de editar el código.

Licencia

Este repositorio incluye un archivo LICENSE. Consulta LICENSE para los términos exactos.

Contacto

Autor: Leandro (usuario GitHub: @LPezzente)

Si querés que actualice este README con más detalles (ejemplos de output, benchmarks, instrucciones para ejecutar en clúster/SLURM o un requirements.txt), decime qué preferís y lo agrego.