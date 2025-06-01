
### Funciones Random Forest (incluir la referencia del Github donde se han obtenido estas funciones)


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

# --------------------- Load packages -----------------------------------------------------------------------------------------------

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

#CREACIÓN DE LOS DATOS

set.seed(123)

n_individuos <- 1000
n_snps <- 80 #Número de  variables 

# Simulación de genotipos
genotype_matrix <- replicate(n_snps,sample(0:2, n_individuos, replace = TRUE, prob = c(0.33, 0.34, 0.33))) #Probabilidad de que el genotipo sea 0,1 o 2 (0- homo alelo 1; 1- hetero; 2-homo alelo 2)

genotipo <- as.data.frame(genotype_matrix)
colnames(genotipo) <- paste0("SNP", 1:n_snps)
rownames(genotipo) <- paste0("Ind", 1:n_individuos)

# Calcular fenotipo según la mayoría por fila
fenotipo <- apply(genotipo, 1, function(row) {
  n0 <- sum(row == 0)
  n1 <- sum(row == 1)    
  n2 <- sum(row == 2)
  
  if (n0 > n1 & n0 > n2) {
    return(0) #fenotipo = 0 si la mayoría de los SNPs son 0
  } else if (n2 > n0 & n2 > n1) {
    return(1) #fenotipo = 1 si la mayoría son 2
  } else {
    return(sample(0:1, 1))  # fenotipo = 1 si caso empate o mayoría de 1
  }
})


# Añadir columnas "fenotipo" e "individuo"
genotipo$fenotipo <- fenotipo
genotipo$individuo <- rownames(genotipo)

# Reordenar columnas: individuo al principio
genotipo <- genotipo[, c("individuo", paste0("SNP", 1:n_snps), "fenotipo")]
setwd("//wsl.localhost/Ubuntu/home/angela/example/TRABAJO")
genotipo
write_xlsx(genotipo,"genotipo_dataset.xlsx")

genotipo[1] <- NULL  # Eliminar primera columna (nombres de fila)
genotipo <- within(genotipo, {
  fenotipo <- as.factor(fenotipo)
})




niveles<- c("0", "1", "2") #Los niveles que deben de tener los SNPs: 0- homo alelo 1; 1- hetero; 2-homo alelo 2
excluir <- c("fenotipo") #establecemos las variables como factores, excepto fenotipo porque solo tiene dos niveles
cols_a_modificar <- setdiff(names(genotipo), excluir)

# Aplicar la transformación solo a esas columnas
genotipo[cols_a_modificar] <- lapply(
  genotipo[cols_a_modificar],
  function(col) {
    col <- as.character(col)  # <- conversión crucial
    col[!col %in% niveles] <- NA
    factor(col, levels = niveles)
  }
)
str(genotipo, list.len=ncol(genotipo))  


#SEPARACIÓN DEL CONJUNTO EN 70% TRAIN Y 30% TEST
set.seed(123)  # Para que la división sea reproducible

# Total de individuos
n_individuos <- nrow(genotipo)

# Crear índices aleatorios para el 70% de entrenamiento
train_indices <- sample(1:n_individuos, size = 0.7 * n_individuos)

# Dividir el dataset
data_train_70 <- genotipo[train_indices, ]
data_test_30 <- genotipo[-train_indices, ]

# Verifica las dimensiones
cat("Train:", nrow(data_train_70), "individuos\n")
cat("Test:", nrow(data_test_30), "individuos\n")


#BÚSQUEDA DEL MEJOR MTRY
set.seed(123)
data_train_70$fenotipo <- as.factor(data_train_70$fenotipo)
par(mar=c(4, 4, 2, 1))
params <- tuneRF(x = data_train_70[, -which(names(data_train_70) == "fenotipo")],  
                 y = data_train_70$fenotipo,
                 mtryStart  = 2,  # Reducir el número inicial de variables
                 stepFactor = 1.5,
                 ntreeTry   = 10000,
                 improve    = 0.01,  
                 plot=TRUE,
                 trace=TRUE, 
                 importance=TRUE)


print(params)

# Extraer el mtry con menor error
best_mtry <- params[which.min(params[, "OOBError"]), "mtry"]
print(best_mtry)


# Preparar x e y por separado
x_train <- data_train_70[, setdiff(names(data_train_70), "fenotipo")]
y_train <- data_train_70[, "fenotipo"]

min_class_size <- min(table(y_train))
min_class_size
# Crear el modelo de Random Forest con interfaz x, y
rf.mdl <- randomForest(
  x = x_train,
  y = y_train,
  strata = y_train, 
  sampsize = c(min_class_size,min_class_size),
  mtry = best_mtry,
  keep.inbag = TRUE,
  importance = "permutation",
  ntree = 10000, 
  prox = TRUE,
  localImp = TRUE,
  norm.votes = TRUE,
  importanceSD = TRUE,
  confusion = TRUE
)

save(rf.mdl, file = "rf.mdl")

# Ahora puedes ejecutar la validación cruzada
rf_cv <- rf.crossValidation(
  x = rf.mdl,
  p = 0.20,
  n = 5,
  seed = 123,
  normalize = FALSE,
  bootstrap = FALSE,
  trace = TRUE
)


sink("resultados_rf_cv_1.txt")
# Accuracy general (Porce_NTAje de Clasificación Correcta)
cat("=== ACCURACY GENERAL (PCC) ===\n")
print(rf_cv$cv.oob$CV.PCC)  # Cross-validated Percent Correct Classification
cat("\nMedia PCC: ", mean(rf_cv$cv.oob$CV.PCC), "\n")
cat("SD PCC: ", sd(rf_cv$cv.oob$CV.PCC), "\n\n")

# Kappa (medida de precisión que tiene en cue_NTA la precisión esperada por azar)
cat("=== ESTADÍSTICO KAPPA ===\n")
print(rf_cv$cv.oob$CV.kappa)
cat("\nMedia Kappa: ", mean(rf_cv$cv.oob$CV.kappa), "\n")
cat("SD Kappa: ", sd(rf_cv$cv.oob$CV.kappa), "\n\n")

# Precisión de usuarios por clase (Users Accuracy)
cat("=== PRECISIÓN DE USUARIOS POR CLASE ===\n")
print(rf_cv$cv.users.accuracy)
cat("\nMedia de Users Accuracy por clase:\n")
print(apply(rf_cv$cv.users.accuracy, 2, mean))
cat("\nSD de Users Accuracy por clase:\n")
print(apply(rf_cv$cv.users.accuracy, 2, sd))
cat("\n")

# Precisión de productores por clase (Producers Accuracy)
cat("=== PRECISIÓN DE PRODUCTORES POR CLASE ===\n")
print(rf_cv$cv.producers.accuracy)
cat("\nMedia de Producers Accuracy por clase:\n")
print(apply(rf_cv$cv.producers.accuracy, 2, mean))
cat("\nSD de Producers Accuracy por clase:\n")
print(apply(rf_cv$cv.producers.accuracy, 2, sd))
cat("\n")

# Tabla completa de resultados OOB
cat("=== TABLA COMPLETA DE RESULTADOS OOB ===\n")
print(rf_cv$cv.oob)
cat("\n")

# Cerrar la conexión
sink()





#Extrae la importancia de las variables basada en MeanDecreaseAccuracy y MeanGini
importancia_pred <- rf.mdl$importance %>% enframe(name = row.names, value = "MeanDecreaseAccuracy")
importancia <- as.data.frame(rf.mdl$importance) %>%
  rownames_to_column(var = "Variable") 
write_xlsx(importancia, "importancia.xlsx")

# Mostrar importancia de variables
svg("importancia.svg")
varImpPlot(rf.mdl, sort = TRUE, main = "Importancia de variables")
dev.off()



#-----Validación 
predicciones <- predict(object = rf.mdl, newdata = data_test_30)

cm <- caret::confusionMatrix(predicciones, data_test_30$fenotipo)


#-CVPVI
X <- within(data_train_70, fenotipo <- NULL )
X = data.frame(X) #generar dataframe para lo que quieres predecir y otra para lo demas
Y = (data_train_70$fenotipo)

set.seed(123)
cv_vi = CVPVI(X,Y,k = 10, mtry = best_mtry, ntree =10000)
save(cv_vi, file = "cv_vi")
summary(cv_vi)
cv_vi$cv_varim

# Novel Test approach, 
#prueba sobre la importancia de las variables, para evaluar qué variables tienen significancia en el modelo y cuáles no.

cv_p = NTA(cv_vi$cv_varim)
sink("pvalue.txt")
print(summary(cv_p,pless = 0.05))
sink()
cv_p_summary <- summary(cv_p, pless = 0.05)  # Resumen de los resultados
cmat <- cv_p_summary$cmat #matriz 'cmat' que contiene los valores de importancia y p-value
cv_p_df <- as.data.frame(cmat)
colnames(cv_p_df) <- c("CV-PerVarImp", "p-value")
ncol(cv_p_df)
significant_vars <- cv_p_df[cv_p_df$`p-value` < 0.05, ] # Filtrar las variables cuyo p-value sea menor que 0.05

significant_var_names <- rownames(significant_vars) # Obtener los nombres de las variables significativas
significant_var_names



#OBTENERMOS LAS VARIABLES SIGNIFICATIVAS


data_train_snps <- data_train_70[, c(significant_var_names, "fenotipo")]
data_test_NTA <- data_test_30[, c(significant_var_names, "fenotipo")]


set.seed(123)
#REALIZAR RANDOM FOREST AHORA PARA LOS SIGNIFICATIVOS
data_train_snps$fenotipo <- as.factor(data_train_snps$fenotipo)
par(mar=c(4, 4, 2, 1))
params_NTA <- tuneRF(x = data_train_snps[, -which(names(data_train_snps) == "fenotipo")],  
                      y = data_train_snps$fenotipo,
                      mtryStart  = 2,  # Reducir el número inicial de variables
                      stepFactor = 1.5,
                      ntreeTry   = 10000,
                      improve    = 0.01,  
                      plot=TRUE,
                      trace=TRUE, 
                      importance=TRUE)

#load("rf.ml_NTA")
print(params)

# Extraer el mtry con menor error
best_mtry_NTA <- params[which.min(params[, "OOBError"]), "mtry"]
print(best_mtry_NTA)



# Preparar x e y por separado
x_train_snps <- data_train_snps[, setdiff(names(data_train_snps), "fenotipo")]
y_train_snps <- data_train_snps[, "fenotipo"]


summary(data_train_snps$fenotipo)

min_class_size <- min(table(y_train_snps))
min_class_size


rf.mdl_NTA <- randomForest(
  x = x_train_snps,
  y = y_train_snps,
  strata = y_train_snps, 
  sampsize = c(min_class_size,min_class_size),
  mtry = best_mtry_NTA,
  keep.inbag = TRUE,
  importance = "permutation",
  ntree = 10000, 
  prox = TRUE,
  localImp = TRUE,
  norm.votes = TRUE,
  importanceSD = TRUE,
  confusion = TRUE
)

# Ahora puedes ejecutar la validación cruzada
rf_cv_NTA <- rf.crossValidation(
  x = rf.mdl_NTA,
  p = 0.20,
  n = 5,
  seed = 123,
  normalize = FALSE,
  bootstrap = FALSE,
  trace = TRUE
)




#Extrae la importancia de las variables basada en MeanDecreaseAccuracy.
importancia_pred_NTA <- rf.mdl_NTA$importance %>% enframe(name = row.names, value = "MeanDecreaseAccuracy")
importancia_NTA <- as.data.frame(rf.mdl_NTA$importance) %>%
  rownames_to_column(var = "Variable") 
write_xlsx(importancia_NTA, "importancia_pvalue.xlsx")

svg("importancia_pvalue.svg")
varImpPlot(rf.mdl_NTA, sort = TRUE, main = "Importancia de variables")
dev.off()
save(rf.mdl, file = "rf.mld")
save(rf.mdl_NTA, file = "rf.mld_NTA")

save(rf.mdl_NTA, file = "rf.mdl_NTA")
#load("rf_cv.mdl_NTA")

save(rf_cv_NTA, file = "rf_cv.mdl_NTA")

# Abrir conexión para guardar resultados en un archivo txt
sink("resultados_rf_cv_2.txt")

# Accuracy general (Porce_NTAje de Clasificación Correcta)
cat("=== ACCURACY GENERAL (PCC) ===\n")
print(rf_cv_NTA$cv.oob$CV.PCC)  # Cross-validated Percent Correct Classification
cat("\nMedia PCC: ", mean(rf_cv_NTA$cv.oob$CV.PCC), "\n")
cat("SD PCC: ", sd(rf_cv_NTA$cv.oob$CV.PCC), "\n\n")

# Kappa (medida de precisión que tiene en cue_NTA la precisión esperada por azar)
cat("=== ESTADÍSTICO KAPPA ===\n")
print(rf_cv_NTA$cv.oob$CV.kappa)
cat("\nMedia Kappa: ", mean(rf_cv_NTA$cv.oob$CV.kappa), "\n")
cat("SD Kappa: ", sd(rf_cv_NTA$cv.oob$CV.kappa), "\n\n")

# Precisión de usuarios por clase (Users Accuracy)
cat("=== PRECISIÓN DE USUARIOS POR CLASE ===\n")
print(rf_cv_NTA$cv.users.accuracy)
cat("\nMedia de Users Accuracy por clase:\n")
print(apply(rf_cv_NTA$cv.users.accuracy, 2, mean))
cat("\nSD de Users Accuracy por clase:\n")
print(apply(rf_cv_NTA$cv.users.accuracy, 2, sd))
cat("\n")

# Precisión de productores por clase (Producers Accuracy)
cat("=== PRECISIÓN DE PRODUCTORES POR CLASE ===\n")
print(rf_cv_NTA$cv.producers.accuracy)
cat("\nMedia de Producers Accuracy por clase:\n")
print(apply(rf_cv_NTA$cv.producers.accuracy, 2, mean))
cat("\nSD de Producers Accuracy por clase:\n")
print(apply(rf_cv_NTA$cv.producers.accuracy, 2, sd))
cat("\n")

# Tabla completa de resultados OOB
cat("=== TABLA COMPLETA DE RESULTADOS OOB ===\n")
print(rf_cv_NTA$cv.oob)
cat("\n")

# Cerrar la conexión
sink()




###### VALIDATION ###
#--------------------------------COMPROBACIÓN DE AMBOS CONJUNTOS DE DATOS----------------------- 
#Se hacen las predicciones
predicciones_NTA <- predict(object = rf.mdl_NTA, data_test_NTA)

cm_NTA <- caret::confusionMatrix(predicciones_NTA, data_test_NTA$fenotipo)
print(cm)
print(cm_NTA)

library(ggplot2)
library(reshape2)

# Función para graficar la matriz de confusión y guardarla como SVG
plot_and_save_conf_matrix <- function(cm, title, filename) {
  cm_table <- as.data.frame(cm$table)  # Convertir a data.frame
  colnames(cm_table) <- c("Referencia", "Predicción", "Frecuencia")
  
  p <- ggplot(data = cm_table, aes(x = Referencia, y = Predicción, fill = Frecuencia)) +
    geom_tile() +
    geom_text(aes(label = Frecuencia), color = "white", size = 6) +
    scale_fill_gradient(low = "lightblue", high = "blue") +
    theme_minimal() +
    labs(title = title, x = "Clase Real", y = "Clase Predicha")
  
  # Guardar en formato .svg
  ggsave(filename = filename, plot = p, width = 6, height = 6, device = "svg")
}

# Guardar ambas matrices de confusión
plot_and_save_conf_matrix(cm, "Matriz de Confusión - Modelo Completo", "conf_matrix_completo.svg")
plot_and_save_conf_matrix(cm_NTA, "Matriz de Confusión - Modelo SNPs", "conf_matrix_NTA.svg")

#Guardar exactitud
accuracy_completo <- cm$overall["Accuracy"]
accuracy_NTA <- cm_NTA$overall["Accuracy"]

sink("confusion_matrix_results.txt")  
cat("### Matriz de Confusión - Modelo Completo ###\n")
print(cm)
cat("\nExactitud del modelo completo:", round(accuracy_completo, 4), "\n\n")

cat("### Matriz de Confusión - Modelo SNPs ###\n")
print(cm_NTA)
cat("\nExactitud del modelo SNPs:", round(accuracy_NTA, 4), "\n")
sink()  # Cerrar el archiv


















