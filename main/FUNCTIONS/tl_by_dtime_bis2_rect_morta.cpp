#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
DataFrame rcpp_tl_by_dtime2_rect_morta(double sst, double dtime, double morta) {
  double ref = 2.;
  double a=20.19*exp(0.041*sst);
  double b=3.258;
 
  int i=1;
  NumericVector dtimevec(1);
  NumericVector tlbisvec(1); //est ce que ça créer un vecteur de longueur non definie qui va s'incrementer au fur et a mesure?
  tlbisvec[0]=ref;
  dtimevec[0]=0;
 
  while (ref < 5.5 ) {
        // while (ref < 5.5 && ref1<5.5) {
        tlbisvec.push_back(ref + a*pow(ref,(-b))*(1+morta)*dtime);  //cf feuille demonstration
        // tlbisvec.push_back(ref + a*pow(ref,(-b))*dtime + morta*dtime);
        // testvec.push_back(ref + a*pow(ref+ morta,(-b))*dtime); //push_back pour ajouter une case au vecteur ou dataframe ou...
        //push_back pour ajouter une case au vecteur ou dataframe ou...
        dtimevec.push_back(i*dtime);
        
                // tlbisvec[i]+=ref + a*pow(ref,(-b));   
        ref=tlbisvec[i];
        // ref=testvec[i];
                // dtimevec[i]+=i*dtime;
        
        i++;
      } 
  dtimevec=round(dtimevec,3);
  tlbisvec[i-1]=5.5; // force last value to be tl=5.5
  // Rcout << "The value of new : " << tlbisvec << "\n";
  // Rcout << "The value of old : " << testvec << "\n";
  // Rcout << "The value of a : " << a << "\n";   Rcout to see all values & number of values taken by variable
  // Rcout << "The value of b : " << b << "\n";
  // Rcout << "The value of tlbisvec : " << tlbisvec << "\n";
  // Rcout << "The value of i : " << i << "\n";
  
  NumericVector TLinf(dtimevec.size());
  
  TLinf[0]=2;
  
  for(int i = 1; i < tlbisvec.size(); ++i) {
    TLinf[i]=tlbisvec[i-1];
  }
  
  DataFrame df = DataFrame::create( Named("tl") = tlbisvec,
                                    // Named("old_tl") = testvec,// simple assign
                                    Named("time") = dtimevec,
                                    Named("TLinf") = TLinf,
                                    Named("TLsup") = clone(tlbisvec));
  // using clone()
  
  // attention si tu fais des simple assign, la colonne de ton DataFrame
  // est un pointeur vers (ici) sst. Donc si tu modifie sst après ça, le 
  // DataFrame est modifié aussi.
  // Si tu utilises clone c'est une copie.
  // comme ici tu en as certainement besoin juste pour créer le dataFrame à la fin 
  // il ne devrait pas y avoir de conflit et je conseille d'utiliser le 'simple assign'
  // afin d'être plus économe sur la mémoire.
  return df;
  
}

