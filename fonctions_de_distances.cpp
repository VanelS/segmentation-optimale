#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
using namespace Rcpp;

/*
* Normaliser les sériest temporelles
*   - prend en paramètre une série temporelle
*   - retourne la série temporelle normalisée 
*/
// [[Rcpp::export]]
NumericVector ZnormalisationDeTS(NumericVector t){
  return (t - mean(t)) / sd(t);
}

/*
 * Détermine les classes distinctes contenues dans le tableau
 *  - classes_des instances : est un vecteur de classes pour les instances
 *  - Retourne les classes distinctes
 */
NumericVector tableau_classe_distincte(NumericVector classe_des_instances){
  
  int j = 0;
  int nb_instances = classe_des_instances.length();
  NumericVector tableau_classes_distinctes(nb_instances);
  NumericVector copie_classe_des_instances(nb_instances);
  
  for(int i = 0; i < nb_instances; i++){
    copie_classe_des_instances[i] = classe_des_instances[i];
  }
  
  
  for(int i = 0; i < (nb_instances - 1); i++){
    
    if(copie_classe_des_instances[i] != -1){// Chaque fois qu'on r?cup?re une classe, pour assurer l'unicit?, on replace tout ses doublons par -1
      
      tableau_classes_distinctes[j] = copie_classe_des_instances[i];
      
      for(int k = (i + 1); k < nb_instances; k++){ // Mettre tout les doublons ? -1
        if(copie_classe_des_instances[i] == copie_classe_des_instances[k]){
          copie_classe_des_instances[k] = -1;
        }else{}
      }
      
      j++;
      
    }else{}
    
  }
  if(copie_classe_des_instances[(nb_instances - 1)] != -1){
    tableau_classes_distinctes[j] = copie_classe_des_instances[(nb_instances - 1)];
    j++;
  }
  
  
  NumericVector out(j);
  int l;
  
  for(l = 0; l < j; l++){
    out[l] = tableau_classes_distinctes[l];
  }
  
  return out;
}



/*
* Calcul l'effectif de chaque classe
*  - classes : ensemble de classes distinctes
*  - classe de chaqu'une des instances
*  - retourne : l'effectif de chaque classes
*/
NumericVector effectif_chaque_classe(NumericVector classes, NumericVector classes_des_instances){
  
  int nb_classe = classes.length();
  int nb_instances = classes_des_instances.length();
  NumericVector effectif_par_classe(nb_classe);
  
  
  for(int i = 0; i < nb_instances; i++){
    for(int j = 0; j < nb_classe; j++){
      if(classes[j] == classes_des_instances[i]){
        effectif_par_classe[j] = effectif_par_classe[j] + 1;
      }else{}
    }
  }
  
  return effectif_par_classe;
}


/*
* Recherche l'indice du maximun
* - V : un vecteur de reelle
* - Retourne l'indice du maximun
*/
// [[Rcpp::export]]
int indice_maximun(NumericVector v){
  
  int n = v.length();
  int ind_max = 0;
  int max = v[0];
  
  for(int i = 1; i < n; i++){
    if(max < v[i]){
      max = v[i];
      ind_max = i;
    }else{}
  }
  
  return ind_max;
}


/*
* Recherche la classe majoritaire
* - classes_voisins : classes des k plus proches voisins
* - classes_instances : classes des instances de la base de données
*/
// [[Rcpp::export]]
int ClasseMajoritaire(NumericVector classes_voisins, NumericVector classes_instances){
  NumericVector classes_distinctes = tableau_classe_distincte(classes_instances);
  NumericVector e = effectif_chaque_classe(classes_distinctes, classes_voisins);
  int ind_max = indice_maximun(e);
  
  return classes_distinctes[ind_max];
  
}



/*
* Calcul la distance euclidienne entre deux vecteurs de reelles
* - x et y : des vecteurs de reelles
*/
double distance(double x, double y) {
  return sqrt(pow((x - y), 2));
}


/*
* Calcul le minimun entre deux valeurs
*  - x et y deux valeurs
*/
// [[Rcpp::export]]
double min(double x, double y){
  return (x<y)?x:y;
}


/*
* Calcul l'alignement ()Dynamic Time Warping) entre deux séries temporelles
* - t1 et t2 : deux séries temporelles
* - retourne la distance entre t1 et t2
*/
// [[Rcpp::export]]
double dynamictw(const std::vector<double>& t1, const std::vector<double>& t2) {
  int m = t1.size();
  int n = t2.size();
  
  //variable that content the distance
  double d;
  
  // Create the cost matrix
  double** cost = new double*[m];
  for(int i = 0; i < m; ++i)
    cost[i] = new double[n];
  
  cost[0][0] = distance(t1[0], t2[0]);
  
  // calculate first row
  for(int i = 1; i < m; i++){
    cost[i][0] = cost[i-1][0] + distance(t1[i], t2[0]);
  }
  
  // calculate first column
  for(int j = 1; j < n; j++)
    cost[0][j] = cost[0][j-1] + distance(t1[0], t2[j]);
  // fill matrix
  for(int i = 1; i < m; i++){
    for(int j = 1; j < n; j++){
      cost[i][j] = min(cost[i-1][j], min(cost[i][j-1], cost[i-1][j-1])) + distance(t1[i],t2[j]);
    }
  }
  
  d = cost[m-1][n-1];
  
  // Delate the cost matrix
  for(int i = 0; i < m; ++i)
    delete[] cost[i] ;
  delete[] cost;
  
  
  return d;
}
