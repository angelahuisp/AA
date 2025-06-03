#Script para realizar el análisis de regresión entre los pesos de los SNPs y sus respectivas importancias

set.seed(123)

library(readxl)

setwd("C:/Users/angel/OneDrive/Escritorio/AA")

importancia <- read.delim("importancia_CVPVI.txt") #cargamos los valores de importancia
seleccionados <- read.delim("SNPs_seleccionados.txt") #cargamos los SNPS seleccionados como importantes
pesos <-read_excel("pesos.xlsx") #cargamos los pesos de los SNP
n_snps <- nrow(pesos)

color <- c(0) #asignamos un color a cada SNP para la representación
for (i in 1:n_snps) {
  if (importancia$`CV_PerVarImp`[i] <= 0) {
    color[i] <- "red" #rojo si su importancia es negativa
  } else {
    if (importancia$`X`[i] %in% seleccionados$SNPS) {
      color[i] <- "#00a000" #verde si se ha seleccionado
    } else {
      color[i] <- 'black' #negro en el caso contrario
    }
  }
}

pesosVI <- data.frame(
  pesos = abs(pesos), #lo que nos importa es el valor absoluto
  VI = c(importancia$`CV_PerVarImp`),
  SNP = importancia$X) #construimos un data.frame con los datos para la regresión lineal

lmPesos = lm(VI~pesos, data = pesosVI) #realizamos la regresión lineal
summary(lmPesos)

plot(pesosVI$pesos, pesosVI$VI, col=color, pch=16, main="Relación pesos e importancia de variables",
     ylab="Importancia variables", xlab="pesos (valor absoluto)") #gráfica peso/importancia
abline(lmPesos, lwd=4, col="#0099ff") #añadimos la recta de la regresión
legend(x="bottomright", pch=16, legend=c("VI <= 0", "VI > 0 & p > 0.05", "p < 0.05"), 
       col=c("red", "black","#00a000"), cex=0.9) #añadimos la leyenda
