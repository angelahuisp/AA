# Análisis de SNPs mediante Random Forest y selección de variables

Este repositorio contiene el código y los datos empleados en un caso de estudio de selección de SNPs utilizando Random Forest y validación cruzada de importancia de variables.

## Contenido del repositorio

- `Trabajo_AA.R` : Script principal donde se realiza el análisis completo mediante Random Forest.
- `correlacion.R` : Script donde se realiza el análisis de regresión entre los pesos de los SNPs y su importancia.
- `Dataset_SNPs.xlsx` : Conjunto de datos con los genotipos de los SNPs utilizados en el análisis.
- `importancia_CVPVI.txt` : Resultados de importancia de los SNPs calculados mediante el método CVPVI.
- `SNPs_seleccionados.txt` : Listado de los SNPs seleccionados como significativos tras el análisis NTA.
- `pesos.xlsx` : Archivo con los pesos asignados a cada SNP.

## Descripción general

El objetivo de este trabajo es aplicar un modelo de Random Forest sobre un conjunto de datos simulado de SNPs, identificar las variables (SNPs) más relevantes mediante un procedimiento de validación cruzada de importancia de variables (CVPVI), y evaluar su significancia estadística mediante el algoritmo NTA.

Posteriormente, se realizó un análisis de regresión para estudiar la relación entre los pesos de los SNPs y su importancia estimada.
