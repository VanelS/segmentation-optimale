library(Rcpp)
library(dtw)
library(timeSeries)
library(TSclust)
library(caret)

sourceCpp("./fonctions_de_distances.cpp")


ZnormalisationDeTS <- function(t){
  return (t - mean(t)) / sd(t) 
}


#Entrée : 
# - t : une série temporelle
#Sortie : 
# - série temporelle sans les valeurs maquantes
SupprimerNA <- function(t){
  a <- data.frame(t(t))
  b <- na.omit(a)
  c <- t(b)
  c
}


# Prédire la classe d'une série temporelle à partir de celle de ses k plus proche voisins
#Entrée : 
# - ind.TS : l'indice de la série temporelle dont il faut prédire la classe
# - training_d : le training_set
# - clId_training_d : classe du training set
# - matrice.distance : matrice de distance entre les TS du test_set (en ligne) et celles du training_set (en colonne)
#Sortie : 
# - La classe prédite
PredirClassPAA<-function(ind.TS,training_d,k, clId_training_d, matrice.distance)
{
  distances <- matrix(nrow = dim(training_d)[1], ncol = 1)
  distances <- matrice.distance[ind.TS, ]
  
  s <- sort(as.vector(distances), index.return=TRUE)
  tb <- as.vector(clId_training_d[s$ix[1:(k)]])
  
  cl <- ClasseMajoritaire(tb, clId_training_d)
  
  return(cl)
}


#Calcul la somme du carré des Erreurs
#Entrée : 
# - Un vecteur de reels
#Sortie : 
# Somme du carré des erreurs
SSE <- function(v){
  v_numeric <- as.numeric(v)
  taille    <- length(v_numeric)
  moy       <- mean(v_numeric)
  som       <- 0
  
  for(i in 1:taille){
    som <- som + (moy - v_numeric[i])^2
  }
  return(som)
}


# Calculer la somme des carrées des erreurs d'un segment
#Entree : 
# - v : un vecteur de points
# - nbPoints : le nombre de points d'un segment
# - ind_debut : indice de début du segment
#Sortie : 
#somme du carré des erreurs d'un segment
SSE_segment <- function(v, nbPoints, ind_debut){
  if((ind_debut + nbPoints) > length(v)){
    print("Erreur dans la fonction SSE_segment")
    print("La longueur du segment est supérieure à celle de la série temporelle")
    return(-1)
  }
  else{
    
  }
  return(SSE(v[ind_debut:(ind_debut + nbPoints)]))
}


#Calcul la somme des carrés des erreurs tous les segments d'une série temporelle
#Entrées : 
# - v : une série temporelle
# - nbPoints : la longueur d'un segment
#Sortie : 
#Somme des erreurs carrée d'un segment
somme_SSE <- function( v, nbPoints){
  n         <- length(v)
  ind_debut <- 1
  aux_se    <- 0
  
  while((ind_debut + nbPoints) <= n){
    aux_se <- aux_se + SSE_segment(v, nbPoints, ind_debut)
    ind_debut <- ind_debut + nbPoints
  }
  return(aux_se)
}


#Calcul de la longueur optimale de chaque segment
# - long_min : la longueur minimale à considérer pour un segment
# - long_max : la longueur maximale à considérer pour un segment
# - v : un vecteur de points
longueur_optimale <- function(long_min, long_max, v){
  len_v <- length(v)
  n     <- long_max - long_min + 1
  j     <- 1
  
  if((long_min < 1) || (long_max > len_v)){
    print("Longueur minimale ou maximale du segment non adapt?e")
    print(paste("long_min", as.character(long_min), sep = " = "))
    print(paste("long_max", as.character(long_max), sep = " = "))
    return(-1)
  }
  
  x <- matrix(nrow = n, ncol = 1)
  y <- matrix(nrow = n, ncol = 1)
  z <- matrix(nrow = n, ncol = 1)
  
  for(i in long_min:long_max){
    x[j, 1] <- i
    z[j, 1] <- -1 * (1/len_v) * somme_SSE(v, i)
    j <- j + 1
  }  
  
  ind_max <- indice_maximun(z)
  ind_max <- ind_max + 1
  
  return(floor(len_v/x[ind_max, 1]))
  
}


#Calcul de la matrice de distance entre les TS de longueure optimale
#Entrees: 
#training_d : jeu de données d'apprentissage
#test_d : jeu de données de test
#code_distance : entier pour le choix de la distance
#          si code_distance = 1 alors PaaDist
#          si code_distance = 2 alors dtw
#Sortie:
#une matrice de distance avec dim(test_d)[1] lignes et dim(training_d)[1] colonnes
MatriceDeDistance_w_opt <- function(training_d, test_d, code_distance){
  
  #initialisation des variables
  distances          <- matrix(nrow= dim(test_d)[1],ncol=dim(training_d)[1])
  w_op_training      <- matrix(nrow= dim(training_d)[1],ncol=1)
  w_op_test          <- matrix(nrow= dim(test_d)[1],ncol=1)
  
  # W_op pour le training set
  for (i in 1:(dim(training_d)[1])){
    
    a2               <- SupprimerNA(training_d[i,])
    b2               <- as.numeric(unlist(a2))
    t2               <- ZnormalisationDeTS(b2)# a voir
    
    n                <- length(b2)
    lg_min           <- 1
    lg_max           <- floor(n/ceiling((n * 0.05) + 1))
    w_op_training[i] <- longueur_optimale(lg_min, lg_max, t2) # a voir znormalized
    
  }
  
  # W_op pour le test set
  for (i in 1:(dim(test_d)[1])){
    
    a1 <- SupprimerNA(test_d[i,])
    b1 <- as.numeric(unlist(a1))
    t1 <- ZnormalisationDeTS(b1)# a voir
    
    n            <- length(b1)
    lg_min       <-  1
    lg_max       <- floor(n/ceiling((n * 0.05) + 1))
    w_op_test[i] <- longueur_optimale(lg_min, lg_max, t1) # a voir znormalized
    
  }
  
  plot(w_op_test)
  
  #------- par curiosité -------------
  print("monyenne w_op_training = ")
  print(mean(w_op_training))
  
  print("Premier quartiles training = ")
  print(quantile(w_op_training))
  
  print("monyenne w_op_test = ")
  print(mean(w_op_test))
  
  print("Premier quartiles test= ")
  print(quantile(w_op_test))
  
  
  #Calcul de la distance entre les s?ries temporelles
  for (i in 1:dim(test_d)[1]){
    a1 <- SupprimerNA(test_d[i,])
    b1 <- as.numeric(unlist(a1))
    t1 <- ZnormalisationDeTS(b1)# a voir
    t1.paa    <- as.vector(PAA(t1, w_op_test[i]))
    
    for (j in 1:dim(training_d)[1]){
      a2     <- SupprimerNA(training_d[j,])
      b2     <- as.numeric(unlist(a2))
      t2     <- ZnormalisationDeTS(b2)# a voir
      t2.paa <- as.vector(PAA(t2, w_op_training[j]))
      
      if(code_distance == 1){
        distances[i, j] <- PaaDist(t1.paa, t2.paa, n, w_op_test[i]) # on fait l'hypoth?se que toutes les TS ont la m?me longueur d'ou n
      }else{
        distances[i, j] <- dynamictw(t1.paa, t2.paa)
      }
      
    }
  }
  as.matrix(distances)
}


#Tester les effets du choix d'une longueure optimale sur la classification
#Entrées : 
# training_set : le jeux de données d'apprentissage
# test_set : le jeux de données de test
# code_distance : si code_distance  1 alors PaaDist sinon alors DTW
#Sorties : 
#rien en sortie
test_w_op <- function(training_set, test_set, k, code_distance){
  
  #initialisation des variables
  pClassId         <- c()
  classId_training <- training_set[,1]
  training_data    <- training_set[,-(1)]
  classId_test     <- test_set[,1]
  test_data        <- test_set[,-(1)]
  
  
  heure1 <- Sys.time()
  
  matrice.distance <- MatriceDeDistance_w_opt(training_data, test_data, code_distance)
  
  for (i in 1:dim(test_data)[1]){
    pClassId[i] <- PredirClassPAA(i,training_data,k, classId_training, matrice.distance)
  }
  
  print(table(classId_test, pClassId))
  # accuracy
  acc <- (sum(classId_test==pClassId)) / nrow(test_data)
  
  heure2<-Sys.time()
  
  t <- difftime(heure2, heure1, units = "mins")
  
  return(list(accuracy = acc, temps = t))
  
}

#test_w_op(Computers_TRAIN, Computers_TEST, 1, 1)