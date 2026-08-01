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