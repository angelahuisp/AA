
### Funciones del paquete randomForest (https://github.com/cran/rfUtilities/blob/master/R/rf.crossValidation.R)
#Se han incluido en el script, debido a que no están ya disponibles dentro del paquete randomForest

rf.crossValidation <- function(x, p=0.10, n=99, seed=NULL, normalize = FALSE, 
                               bootstrap = FALSE, p.threshold = 0.60, 
                               trace = FALSE) {
  if (!inherits(x, c("randomForest", "ranger"))) 
    stop("x is not randomForest class object")
  if(length(grep("~", x$call[[2]])) > 0)
    stop("This package does not support a formula interface, please use x, y arguments")
  if(inherits(x, "randomForest")) {
    if (x$type == "unsupervised") 
      stop("Unsupervised classification not supported")
    mtype <- tolower(x$type)
  } else if(inherits(x, "ranger")) {
    mtype <- tolower(x$treetype)
  }
  if(!is.null(seed)) { set.seed(seed) }
  if(bootstrap) 
    cat("Bootstrap sampling is being applied,", paste("p",p,sep="="), 
        "argument is ignored", "\n")
  
  # formating call and pulling data
  a <- as.list(x$call)[-1] 
  xdata = eval(a[[which(names(a) == "x")]]) 
  ydata = eval(a[[which(names(a) == "y")]]) 
  dat <- data.frame(y=ydata, xdata) 
  
  #*********************************	
  if(mtype == "regression") {
    #*********************************
    cat("running:", mtype, "cross-validation", "with", n, "iterations", "\n")
    # Validation statistics (RMSE, MBE, MAE)
    # Root Mean Square Error (RMSE) 
    rmse <- function(y, x, norm = FALSE){
      if( length(y[is.na(y)]) > 0) stop("NA values present in y data")
      if( length(x[is.na(x)]) > 0) stop("NA values present in x data")
      e <- sqrt(mean((y - x)^2)) 
      if( norm ) e <- e / diff(range(y, na.rm = TRUE))
      return( e )		   
    }          
    # Mean Bias Error (MBE) 
    mbe <- function(y, x, norm = FALSE){
      if( length(y[is.na(y)]) > 0) stop("NA values present in y data")
      if( length(x[is.na(x)]) > 0) stop("NA values present in x data")
      e <- mean(x - y)
      # e <- mean(x - y) / mean(y) * 100 
      if( norm ) e <- e / diff(range(y, na.rm = TRUE))
      return( e )		   
    }     
    # Mean Absolute Error (MAE) 
    mae <- function(y, x, norm = FALSE){
      if( length(y[is.na(y)]) > 0) stop("NA values present in y data")
      if( length(x[is.na(x)]) > 0) stop("NA values present in x data") 
      e <- mean(abs(y - x))
      if( norm ) e <- e / diff(range(y))
      return( e )		   
    }
    # Kolmogorov-Smirnov Test (1=D, 2=p.value)
    ks <- function(y, x, s = c(1,2)) {
      stats::ks.test(x, stats::ecdf(y))[s[1]]
    }
    # Define validation vectors
    y.rmse <- rep(NA, n)  
    y.mae <- rep(NA, n)
    y.mbe <- rep(NA, n)
    model.varExp <- rep(NA, n)
    model.mse <- rep(NA, n)
    ks.p <- rep(NA, n)
    ks.d <- rep(NA, n)
    if(bootstrap) boot.sample.size <- rep(NA, n)	
    sample.size = round( (length(ydata) * p), digits=0)
    #**************************
    # cross-validation for loop
    #**************************
    for(i in 1:n) {
      if(trace) cat("running iteration:", i, "of", n, "\n")
      # Draw random sample		
      if(!bootstrap) {
        sidx <- sample(1:nrow(dat), sample.size)   
        dat.sub <- dat[-sidx,]                   
        dat.cv <- dat[sidx,]		            
      } else {	
        dat.sub <- dat[sample(1:nrow(dat), replace=TRUE),]              
        dat.cv <- dat[which(!rownames(dat) %in% rownames(dat.sub)),]	   
      }		
      a[["y"]] <- dat.sub[,"y"]
      a[["x"]] <- dat.sub[,2:ncol(dat.sub)]
      if(inherits(x, "ranger")) {    
        rf.fit <- do.call(ranger::ranger, a)
        model.mse[i] <- rf.fit$prediction.error
        model.varExp[i] <- rf.fit$r.squared
        prd <- stats::predict(rf.fit,data = dat.cv[,2:ncol(dat.cv)])$predictions 			
        
      } else if(inherits(x,"randomForest")) {
        rf.fit <- do.call(randomForest::randomForest, a)
        model.mse[i] <- rf.fit$mse[length(rf.fit$mse)]
        model.varExp[i] <- round(100*rf.fit$rsq[length(rf.fit$rsq)], digits=2) 
        prd <- stats::predict(rf.fit, newdata = dat.cv[,2:ncol(dat.cv)]) 			
      }
      y.rmse[i] <- rmse(dat.cv[,"y"], prd, norm=normalize) 
      y.mbe[i] <- mbe(dat.cv[,"y"], prd, norm=normalize) 
      y.mae[i] <- mae(dat.cv[,"y"], prd, norm=normalize) 
      ks.p[i] <- as.numeric(ks(dat.cv[,"y"], prd, s=2))  
      ks.d[i] <- as.numeric(ks(dat.cv[,"y"], prd, s=1))
      if(bootstrap) boot.sample.size[i] <- nrow(dat.cv) 
    }		   
    if(inherits(x, "ranger")) {
      fit.var.exp = x$r.squared
      fit.mse = x$prediction.error
    } else if(inherits(x,"randomForest")) {
      fit.var.exp = round(100*x$rsq[length(x$rsq)], digits=2) 
      fit.mse = stats::median(x$mse)  
    }
    r.cv <- list(fit.var.exp=fit.var.exp, 
                 fit.mse=fit.mse, y.rmse = y.rmse, 
                 y.mbe = y.mbe, y.mae = y.mae, D = ks.d, 
                 p.val = ks.p, model.mse = model.mse, 
                 model.varExp = model.varExp )	
    class(r.cv) <- c("rf.cv", "regression")
    
    #*********************************	
  } else if(any(mtype %in% c("classification","probability estimation"))) {
    #*********************************		
    cat("running:", mtype, "cross-validation", "with", n, "iterations", "\n")
    if(inherits(x, "ranger")) { 
      if("probability estimation" == mtype) {  
        ppred <- ifelse(x$predictions[,2] < p.threshold, 0, 1)
      } else {
        ppred <- x$predictions 
      }
      cm <- accuracy(ydata, ppred) 
    } else if(inherits(x,"randomForest")) {
      cm <- accuracy(ydata, x$predicted)
    }
    mdl.oob <- cm$PCC	
    mdl.kappa <- cm$kappa	
    mdl.ua <- cm$users.accuracy 
    mdl.pa <- cm$producers.accuracy 
    model.error = data.frame(model.users.accuracy = mdl.ua, 
                             model.producers.accuracy=mdl.pa, 
                             model.oob=mdl.oob, model.kappa=mdl.kappa) 
    if(bootstrap) boot.sample.size <- vector()	
    classes <- as.vector(levels( ydata ))	
    sample.size = round( (length(ydata) * p) / length(x$classes), digits=0) 
    cv.ua <- as.data.frame(array(0, dim=c(0,length(classes))))
    cv.pa <- as.data.frame(array(0, dim=c(0,length(classes))))
    mdl.ua <- as.data.frame(array(0, dim=c(0,length(classes)))) 
    mdl.pa <- as.data.frame(array(0, dim=c(0,length(classes))))
    mdl.oob <- as.data.frame(array(0, dim=c(0,2+length(classes))))
    cv.oob <- as.data.frame(array(0, dim=c(0,2+length(classes))))	
    
    # Create class-level sample sizes
    nclass <- length(unique(ydata)) 
    p.class = p / nclass
    sample.sizes <- vector()
    for(s in 1:nclass) { 
      sample.sizes[s] <- round(length(ydata[ydata == unique(ydata)[s]]) * 
                                 p.class, digits=0)    
    }
    cv.oob <- c(0,0,0)
    #**************************
    # cross-validation for loop
    #**************************	  
    for(i in 1:n) {
      if(trace) cat("running iteration:", i, "\n")
      # Draw random sample		
      if(!bootstrap) {
        sidx <- list()				
        for(s in 1:nclass) {
          sidx[[s]] <- sample(which(ydata %in% unique(ydata)[s]), sample.sizes[s])
        }
        sidx <- unlist(sidx)  
        tx <- xdata[sidx,]
        ty <- ydata[sidx]
        mx <- xdata[-sidx,]
        my <- ydata[-sidx]
      } else {		
        tx <- dat[sample(1:nrow(dat), replace=TRUE),]
        ty <- tx$y 
        tx <- tx[,2:ncol(tx)]
        mx <- dat[-which(!rownames(dat) %in% rownames(tx)),]
        my <- mx$y
        mx <- mx[,2:ncol(mx)] 
      }
      
      dat.sub <- data.frame(y=my, mx)
      dat.cv <- data.frame(y=ty, tx) 
      a[["y"]] <- dat.sub[,"y"]
      a[["x"]] <- dat.sub[,2:ncol(dat.sub)]
      if(inherits(x, "ranger")) {    
        rf.fit <- do.call(ranger::ranger, a)
        prd <- stats::predict(rf.fit,data = dat.cv[,2:ncol(dat.cv)])$predictions
        if("probability estimation" == mtype) {	
          prd <- ifelse(prd[,2] < p.threshold, 0, 1)
        } 
        cv.acc <- accuracy(prd,  ty)
        mdl.oob <- rf.fit$prediction.error
      } else if(inherits(x,"randomForest")) {
        rf.fit <- do.call(randomForest::randomForest, a)
        prd <- stats::predict(rf.fit, newdata = dat.cv[,2:ncol(dat.cv)])
        cv.acc <- accuracy(prd,  ty)
        mdl.oob <- c(apply(rf.fit$err.rate, MARGIN = 2, stats::median))[1] 
      }
      if(bootstrap) boot.sample.size <- append(boot.sample.size, length(my))
      cv.oob <- rbind(cv.oob, c(mdl.oob, cv.acc$PCC, cv.acc$kappa))	
      cv.ua <- rbind(cv.ua, cv.acc$users.accuracy) 
      cv.pa <- rbind(cv.pa, cv.acc$producers.accuracy) 
      #if(!length(classes) == length(unique(my))) {
      #### add code to check for presence of classes and set 
      ####  to NA if absent classes occur   
      #}
    }  
    names(cv.ua) <- c(classes) 
    names(cv.pa)  <- c(classes)
    cv.oob <- as.data.frame(cv.oob)	  
    names(cv.oob) <- c("Model.PCC", "CV.PCC", "CV.kappa")
    cv.oob <- cv.oob[-1,]
    rownames(cv.oob) <- 1:nrow(cv.oob)		  
    r.cv <- list(cv.users.accuracy = cv.ua, cv.producers.accuracy=cv.pa, 
                 cv.oob = cv.oob, model.error = model.error)
    if(bootstrap) r.cv[["boot.sample.size"]] <- boot.sample.size
    class(r.cv) <- c("rf.cv",  "classification")
  }   
  return( r.cv )
}  

accuracy <- function (x, y) {
  if(inherits(x, c("table", "matrix"))) {
    t.xy <- x 
  } else {
    t.xy <- table(x, y)
  }
  if(inherits(t.xy, "matrix")) {
    if(is.null(rownames(t.xy))) {	
      rownames(t.xy) <- c("0","1")
      colnames(t.xy) <- c("0","1")
    }
  }
  tc <- match(colnames(t.xy), rownames(t.xy))
  mtc <- matrix(ncol = ncol(t.xy), nrow = length(tc[tc == "NA"]), 0)
  nrn <- colnames(t.xy)[is.na(tc) == TRUE]
  rownames(mtc) <- nrn
  t1 <- rbind(t.xy, mtc)
  tr <- match(rownames(t1), colnames(t1))
  mtr <- matrix(nrow = nrow(t1), ncol = length(tr[tr == "NA"]), 0)
  ncn <- rownames(t1)[is.na(tr) == TRUE]
  colnames(mtr) <- ncn
  t2 <- cbind(t1, mtr)
  sr <- sort(rownames(t2))
  mr <- match(sr, rownames(t2))
  t3 <- t(t2[mr, ])
  sc <- sort(rownames(t3))
  mc <- match(sc, rownames(t3))
  t4 <- t(t3[mc, ])
  agree <- diag(t4)
  prod1 <- apply(t4, 1, sum)
  prod2 <- agree / prod1
  user1 <- apply(t4, 2, sum)
  user2 <- agree / user1
  N <- sum(t4)
  k1 <- sum(agree)
  k2 <- sum(prod1 * user1)
  khat <- abs(((N * k1) - k2) / (N^2 - k2))	
  if( dim(t.xy)[1] == 2 ) {
    # TP(1)  FN(3)  
    # FP(2)  TN(4)  	
    n = sum(t.xy)  # N	
    TP <- t.xy[1]  # True Positives, Power 
    FP <- t.xy[2]  # False Positives, Type-I error 
    FN <- t.xy[3]  # False Negatives, Type-II error 
    TN <- t.xy[4]  # True Negatives
    # prevalence <- TP / n
    precision <- TP / (TP + FP)  	  
    tpr <- TP / (TP + FN)  # true positive rate (aka, sensitivity, recall)
    tnr <- TN / (TN + FP)  # true negative rate (aka, specificity, selectivity)	 
    fpr <- FP / (FP + TN) 
    fnr <- FN / (FN + TP)  # Beta       
    type1.error <- 1 - tnr                     
    type2.error <- 1 - tpr                     
    plr <- tpr / (1 - tnr)                     
    nlr <- (1 - tpr) / tnr 
    auc <- (tpr - fpr + 1) / 2
    gini <- 2 * auc - 1  
    f.score <- 2 * (precision * tpr) / (precision + tpr)
    true.skill <- ( (t.xy[1] * t.xy[4]) - (t.xy[3] * t.xy[2]) ) / 
      ( (t.xy[1] + t.xy[2]) * (t.xy[3] + t.xy[4]) )
    gain <- precision / ( (t.xy[1] + t.xy[4]) / n )
    mcc <- (TP * TN - FP * FN) / sqrt(prod(c((TP + FP),(TP + FN),(TN + FP),(TN + FN))))  
    confusion <- matrix(c(paste0("True positive(", TP, ")"),
                          paste0("False positive(", FP, ")"),
                          paste0("False negative(", FN, ")"),
                          paste0("True negative(", TN, ")")),
                        nrow=2, byrow=TRUE)
    rownames(confusion) <- rownames(t.xy)
    colnames(confusion) <- colnames(t.xy)
    acc <- list( PCC = (sum(diag(t4))/sum(t4)) * 100,
                 auc = auc, 	
                 users.accuracy = round(user2 * 100, 1),  
                 producers.accuracy = round(prod2 * 100, 1),
                 kappa = round(khat, 4),
                 true.skill = true.skill, 
                 sensitivity = tpr,
                 specificity = tnr,
                 plr = plr,
                 nlr = nlr,
                 typeI.error = type1.error,
                 typeII.error = type2.error,
                 gini = gini,
                 f.score = f.score,
                 gain = gain,
                 matthews = mcc,
                 confusion = confusion )	
  } else {
    acc <- list( PCC = (sum(diag(t4))/sum(t4)) * 100, 
                 users.accuracy = round(user2 * 100, 1),  
                 producers.accuracy = round(prod2 * 100, 1),
                 kappa = round(khat, 4),
                 confusion = t.xy )
  }			 
  class(acc) <- c("accuracy", "list") 			   
  return( acc )
}

# =========================================================================================================
# Simulación y análisis de datos genómicos mediante Random Forest (RF) y selección de variables significativas

# Este script realiza la simulación de un conjunto de datos genómicos, ajusta modelos de RF con validación cruzada (CV),
# evalúa la mportancia de las variables y selecciona SNPs significativos usando CV y pruebas estadísticas

#Se pretende evaluar las métricas de los modelos RF. Un modelo se entrenará con todo el conjunto de variables dispnibles, 
#mientras que el otro modelo únicamente utilizará las variables significativas resultantes del test estadístico. 
# ============================ Carga de paquetes requeridos ===============================================

library(randomForest)
library(randomForestSRC)
library(ranger)
library(caret)
library(tibble)
library(openxlsx)
library(writexl)
library(readxl)
library(vita)
library(dplyr)

# ============================ Simulación de datos genómicos ==============================================

setwd("C:/Users/angel/OneDrive/Escritorio/AA")

data <-read_excel("genotipo_dataset_2.xlsx", col_types = "text")
data[1] <- NULL  # Eliminar primera columna (nombres de fila)
data$fenotipo <- as.factor(data$fenotipo) #convertir a factor 


niveles<- c("0", "1", "2") #Los niveles que deben de tener los SNPs: 0- homocigoto para alelo 1; 1- heterocigoto; 2-homocigoto alelo 2
cols_a_factor <- setdiff(names(data), "fenotipo")
data[cols_a_factor] <- lapply(data[cols_a_factor], function(col) {
  col <- as.character(col)
  col[!col %in% niveles] <- NA  # poner NA si hay valores fuera de 0,1,2
  factor(col, levels = niveles)
})
sum(is.na(data))


set.seed(123)
train_indices <- sample(1:nrow(data), size = 0.7 * nrow(data)) #Se separa los datos en 2 conjuntos, train (70%) y test (30%)
data_train <- data[train_indices, ]
data_test  <- data[-train_indices, ]
# ============================ Búsqueda del mejor valor de mtry ===========================================
#mtry = nº de variables que se seleccionan aleatoriamente en cada split de cada árbol

set.seed(123)
par(mar=c(4, 4, 2, 1))
params <- tuneRF(x = data_train[, -which(names(data_train) == "fenotipo")],  
                 y = data_train$fenotipo,
                 mtryStart  = 2,  # Empieza probando con mtry = 2.
                 stepFactor = 1.5, #Cada vez multiplica el mtry anterior por 1.5
                 ntreeTry   = 5000, #Nº de árboles de cada modelo que prueba
                 improve    = 0.01,  #sigue probando si la mejora en el error OOB es al menos del 1%.
                 plot=TRUE,#Dibuja el gráfico de error OOB vs mtry
                 trace=TRUE)


cat(params)

# Selección del mejor mtry.
best_mtry <- params[which.min(params[, "OOBError"]), "mtry"]
cat(best_mtry)



# ============================ Entrenamiento del modelo Random Forest =====================================
# Ajuste del balance de las clases, ya que un desbalance de clases puede provocar que el modelo 
#tienda a acertar mejor la clase mayoritaria que la minoritaria. 

# Preparar x e y por separado
x_train <- data_train[, setdiff(names(data_train), "fenotipo")]
y_train <- data_train$fenotipo
table(y_train) #ver proporción de clases en los datos de train
min_class_size <- min(table(y_train)) #se guarda la clase minoritaria
min_class_size
# Crear el modelo de Random Forest con interfaz x, y
rf.mdl <- randomForest(
  x = x_train,
  y = y_train,
  strata = y_train, 
  sampsize = c(min_class_size,min_class_size),#Toma el mismo número de muestras por clase en cada árbol
  mtry = best_mtry,
  ntree = 5000, #RF con 5000 árboles
  seed=123) #Número total de árboles a entrenar


# Validación cruzada del modelo 
rf_cv <- rf.crossValidation(
  x = rf.mdl,
  p=0.2, #Proporción de datos usada como conjunto de validación en cada iteración de la cross-validation.
  n = 5, #número repeticiones de validación cruzada
  seed = 123,
  trace = TRUE, #Muestra progreso
)


#Validación en conjunto test 

# Predicción y matriz de confusión
predicciones <- predict(object = rf.mdl, newdata = data_test)
cm <- caret::confusionMatrix(predicciones, data_test$fenotipo)

# ============================ Selección de variables significativas (CVPVI) ==============================
# Novel Test approach
# Prueba sobre la importancia de las variables, para evaluar qué variables tienen significancia en el modelo y cuáles no.

X <- within(data_train, fenotipo <- NULL )
X = data.frame(X) #generar dataframe para lo que quieres predecir y otra para lo demas
Y = (data_train$fenotipo)

cv_vi = CVPVI(X,Y,k = 2, mtry = best_mtry, ntree =5000, seed=123)
#cv_vi$cv_varim #si se quiere observar la importancia calculada para cada snp

cv_p = NTA(cv_vi$cv_varim)
cv_p_summary <- summary(cv_p)  
# Resumen de los resultados con pvalue < 0.05
cmat <- cv_p_summary$cmat #matriz 'cmat' que contiene los valores de importancia y p-value
cv_p_df <- as.data.frame(cmat)
colnames(cv_p_df) <- c("CVPVI_Importancia", "p-value")
significant_vars <- cv_p_df[cv_p_df$`p-value` < 0.05, ] # Filtrar las variables cuyo p-value sea menor que 0.05
significant_var_names <- rownames(significant_vars) # Obtener los nombres de las variables significativas
significant_var_names



# Subconjunto de variables significativas
#Se extrae de los datos únicamente las variables significativas, se ha establecido un umbral de valor p < 0.05
data_train_NTA <- data_train[, c(significant_var_names, "fenotipo")]
data_test_NTA <- data_test[, c(significant_var_names, "fenotipo")]


# ============================ Nuevo Random Forest con SNPs significativos ================================

set.seed(123)
data_train_NTA$fenotipo <- as.factor(data_train_NTA$fenotipo)
par(mar=c(4, 4, 2, 1))
params_NTA <- tuneRF(x = data_train_NTA[, -which(names(data_train_NTA) == "fenotipo")],  
                     y = data_train_NTA$fenotipo,
                     mtryStart  = 2,  # Reducir el número inicial de variables
                     stepFactor = 1.5,
                     ntreeTry   = 10000,
                     improve    = 0.01,  
                     plot=TRUE)

cat(params)

# Extraer el mtry con menor error
best_mtry_NTA <- params[which.min(params[, "OOBError"]), "mtry"]
cat(best_mtry_NTA)

# Preparar x e y por separado
x_train_snps <- data_train_NTA[, setdiff(names(data_train_NTA), "fenotipo")]
y_train_snps <- data_train_NTA$fenotipo



rf.mdl_NTA <- randomForest(
  x = x_train_snps,
  y = y_train_snps,
  strata = y_train_snps, 
  sampsize = c(min_class_size,min_class_size),
  mtry = best_mtry_NTA,
  ntree = 5000
)

# Ejecución la validación cruzada 
rf_cv_NTA <- rf.crossValidation(
  x = rf.mdl_NTA,
  p = 0.20,
  n = 5,
  seed = 123,
  trace = TRUE
)




# ============================ Evaluación final y visualización de matrices de confusión ==================
# Predicciones en el conjunto de prueba reducido

predicciones_NTA <- predict(object = rf.mdl_NTA, data_test_NTA)
cm_NTA <- caret::confusionMatrix(predicciones_NTA, data_test_NTA$fenotipo)


#Métricas para modelo completo (todas las variables)
print(cm$table) #matriz de confusión 
cat("\nExactitud :", round(cm$overall["Accuracy"], 4), "\n")
cat("Sensibilidad :", round(cm$byClass["Sensitivity"], 4), "\n")
cat("Especificidad :", round(cm$byClass["Specificity"], 4), "\n")


#Métricas para modelo con variables significativas

print(cm_NTA$table)
cat("\nExactitud :", round(cm_NTA$overall["Accuracy"], 4), "\n")
cat("Sensibilidad :", round(cm_NTA$byClass["Sensitivity"], 4), "\n")
cat("Especificidad :", round(cm_NTA$byClass["Specificity"], 4), "\n")
#Matrices de Confusión

library(ggplot2)
library(reshape2)
plot_conf_matrix <- function(cm, title) {
  cm_table <- as.data.frame(cm$table)
  colnames(cm_table) <- c("Referencia", "Predicción", "Frecuencia")
  
  p <- ggplot(data = cm_table, aes(x = Predicción, y = Referencia, fill = Frecuencia)) +  # <-- intercambiados
    geom_tile() +
    geom_text(aes(label = Frecuencia), color = "white", size = 6) +
    scale_fill_gradient(low = "lightblue", high = "blue") +
    theme_minimal() +
    labs(title = title, x = "Clase Predicha", y = "Clase Real")
  
  print(p)
}


# Matrices de confusión
plot_conf_matrix(cm, "Matriz de Confusión - Modelo Completo")
plot_conf_matrix(cm_NTA, "Matriz de Confusión - Modelo con SNPs signifcativos")


#Archivo con las métricas y los
sink("Resultados.txt")  
cat("Datos de entrenamiento y de validación\n")
cat("Train:", nrow(data_train), "individuos")
cat("Distribución en Train:\n")
print(table(data_train$fenotipo))
cat("Test:", nrow(data_test), "individuos")
cat("Distribución en Test:\n")
print(table(data_test$fenotipo))
cat("Resultados RF modelo completo\n")
print(cm$table)

# Extraer métricas específicas
cat("\nExactitud :", round(cm$overall["Accuracy"], 4), "\n")
cat("Sensibilidad :", round(cm$byClass["Sensitivity"], 4), "\n")
cat("Especificidad :", round(cm$byClass["Specificity"], 4), "\n")

cat("\nResultados Validación cruzada modelo completo:\n")

cat("Accuracy en cada iteración:", rf_cv$cv.oob$CV.PCC)
cat("\nAccuracy medio: ", mean(rf_cv$cv.oob$CV.PCC))
cat("\nDesviación estándar: ", sd(rf_cv$cv.oob$CV.PCC))

cat("\n\nVariables significativas:")
print(cv_p_summary)
cat("\n\nResultados RF modelo completo\n")

print(cm_NTA$table)

# Extraer métricas específicas
cat("\nExactitud :", round(cm_NTA$overall["Accuracy"], 4), "\n")
cat("Sensibilidad :", round(cm_NTA$byClass["Sensitivity"], 4), "\n")
cat("Especificidad :", round(cm_NTA$byClass["Specificity"], 4), "\n")

cat("\nResultados Validación cruzada modelo con variables significativas\n:")
cat("Accuracy en cada iteración:", rf_cv_NTA$cv.oob$CV.PCC)
cat("\nAccuracy medio: ", mean(rf_cv_NTA$cv.oob$CV.PCC))
cat("\nDesviación estándar: ", sd(rf_cv_NTA$cv.oob$CV.PCC))
sink()















