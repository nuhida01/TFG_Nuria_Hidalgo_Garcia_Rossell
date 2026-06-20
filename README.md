# TFG
Código asociado al Trabajo de Fin de Grado de Nuria Hidalgo García Rossell:  
**“Descomposición en Valores Singulares. Métodos de cálculo y aplicaciones”**.

El repositorio contiene:
- Implementaciones en MATLAB de algoritmos para el cálculo de la SVD.
- Una aplicación de la SVD al análisis de firmas génicas en cáncer de mama mediante PCA y PCA disperso.

Todo el código está implementado en **MATLAB**.


## Contenido del repositorio

### Capítulo 2 – Algoritmos de cálculo de la SVD

**Funciones principales**

- `QRrayleigh.m`  
  Implementación del algoritmo QR con traslación de Rayleigh.

- `SVD1.m`  
  Primer algoritmo de Golub–Kahan para el cálculo de valores singulares de una matriz.

- `SVD2.m`  
  Segundo algoritmo de Golub–Kahan para el cálculo de valores singulares de una matriz.

**Funciones auxiliares**

- `bidiag.m`  
  Reducción de una matriz a forma bidiagonal.

- `hessenberg.m`  
  Reducción de una matriz a forma de Hessenberg.

- `QRhess.m`  
  Implementación del algoritmo QR para matrices de Hessenberg.


### Capítulo 3 – Aplicación de la SVD al análisis de firmas génicas en cáncer de mama mediante PCA y PCA disperso

**Funciones**

- `ClassicPCA.m`  
  Cálculo de las componentes principales de una matriz de datos (PCA clásico).

- `ThreSPCA.m`  
  Implementación del algoritmo ThreSPCA para el PCA disperso de una matriz.

**Scripts**

- `TCGA.m`  
  Script principal que ejecuta el estudio completo sobre el conjunto de datos de expresión génica de cáncer de mama. Nota: el conjunto de datos de expresión génica no se incuyen en el repositorio.



## Referencia

> Nuria Hidalgo García Rossell, “Descomposición en Valores Singulares. Métodos de cálculo y aplicaciones”, Trabajo de Fin de Grado, Universidad Complutense de Madrid, Facultad de Ciencias Matemáticas, 2026
