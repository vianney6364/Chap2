#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
double te_rcpp(String eco_type,double ssti) {
  double a;
  double b;
  double c;
  double d;
  double e;
  double te_value;
  if (eco_type == "polar"){b=0; a=0;c=0;d=0;e=0;}
  if (eco_type == "temperate"){b=(0.142); a=(-0.013);c=0;d=0;e=0;}
  if (eco_type == "tropical"){b=(-0.352); a=(0.015);c=0;d=0;e=0;}
  if (eco_type == "upwelling"){b=(0.167); a=(-0.032);c=0;d=0;e=0;}
  if (eco_type == "moyen"){b=0;a=0;c=0.05592;d=0.00325;e=0.003502;}
  te_value=exp((-2.162-c) + (-0.025+d)*ssti + b + a*ssti)*(1.038013+e);
  return te_value;
}

// double te_rcpp(String eco_type,double ssti) {
//   double a;
//   double b;
//   double c;
//   double te_value;
//   if (eco_type == "polar"){b=0; a=0;c=0;}
//   if (eco_type == "temperate"){b=(0.142); a=(-0.013);c=0;}
//   if (eco_type == "tropical"){b=(-0.352); a=(0.015);c=0;}
//   if (eco_type == "upwelling"){b=(0.167); a=(-0.032);c=0;}
//   if (eco_type == "moyen"){b=0; a=0;c=0.003502;}
//   te_value=exp((-2.162) + (-0.025)*ssti + b + a*ssti)*(1.038013+c);
//   return te_value;
// }
// [[Rcpp::export]]
DataFrame compute_new_time_rcpp_rect_morta_a_fixed (DataFrame actualize_t1,String eco_type, double sst_kin,
                                double nppi, double ssti, double time_stepi, double dtime,double mhwi, double morta) {
  
  NumericVector tl = actualize_t1["tl"];
  NumericVector time = actualize_t1["time"];
  NumericVector TLinf = actualize_t1["TLinf"];
  NumericVector TLsup = actualize_t1["TLsup"];
  NumericVector prod = actualize_t1["production"];
  int nligne = actualize_t1.nrows();
  // On créé des raccourcis (pointeurs) pour chacune des autres colonnes
  
  NumericVector npp(nligne); 
  NumericVector sst(nligne); 
  // as a parameter in kinetic fixed -> needed sst ref
  NumericVector ref_sst(nligne); 
  NumericVector te(nligne); 
  NumericVector fishing(nligne); 
  NumericVector mortality(nligne); 
  NumericVector heat_tho(nligne); 
  LogicalVector MHW(nligne); 
  NumericVector kinetic(nligne); 
  NumericVector kinetic_heat(nligne);
  NumericVector delta_tl_dyn(nligne); 
  NumericVector cumsum(nligne); 
  NumericVector flow_m_per_tl(nligne); 
  NumericVector flow_tl_2to7(nligne); 
  NumericVector biomass(nligne); 
  NumericVector time_step(nligne); 
  
  
  // Ci-dessus on créé un vecteur du nom de my_prod et de longueur ET_T2.nrows()
  // pour obtenir le nombre de ligne de ET_T2 on écrit ET_T2.nrows() plutôt que nrows(ET_T2)
  // on pourrait aussi stocker directement le nombre de ligne de ET_T2 dans une variable:
  // double matrix_length = ET_T2.nrows();
  // pour éviter de répéter cette opération
  // "double" est un type classique du C++ qui correspond à un réel.
  // Rcout << "je suis pas encore dans la boucle";
     
// Rcout << "je suis dans la boucle mhwi==true";
    for(int i = 0; i < nligne; ++i) {
      npp[i]=nppi;
      sst[i]=ssti;
      te[i]=te_rcpp(eco_type,ssti);
      // # fishing= case_when(tl<=3~0,TRUE~0.2),
      mortality[i]=morta;
      MHW[i]=mhwi;
      fishing[i]=0;
      // #'*without fishing*
      // # kinetic=(20.19*(tl^(-3.26))*exp(.041*sst))*dtime,
      // #'*with fishing*
      // #'to fix kinetic replace sst by sst_ref
      kinetic[i]=(20.19*(pow(tl[i],(-3.26)))*exp(.041*sst_kin)+ fishing[i]);
      heat_tho[i]=kinetic[i]*mortality[i];
      kinetic_heat[i]=kinetic[i]+heat_tho[i];
        
      delta_tl_dyn[i]=TLsup[i]-TLinf[i];
      flow_m_per_tl[i]=prod[i]/delta_tl_dyn[i];
      flow_tl_2to7[i]=(delta_tl_dyn[i]*(-log(te[i]*(1-fishing[i]/kinetic[i])*(1-heat_tho[i]/kinetic[i])))*flow_m_per_tl[i])/
        (1-pow(te[i]*(1-fishing[i]/kinetic[i])*(1-heat_tho[i]/kinetic[i]),delta_tl_dyn[i]));
      // flow_tl_2to7[i]=(delta_tl_dyn[i]*(-log(te[i])+fishing[i]/kinetic[i]+mortality[i]/kinetic[i])*flow_m_per_tl[i])/
      //   (1-exp(-(-log(te[i])+(fishing[i]/kinetic[i]+mortality[i]/kinetic[i]))*delta_tl_dyn[i]));
      biomass[i]=prod[i]/kinetic_heat[i];
      time_step[i]=time_stepi;
      
      // Rcout << "je suis à tl=" << i << "\n";
      
      // ***** NOTE IMPORTANTE *****
      // les indices commence à 0 en C++ et à 1 en R
      // donc en R x[1] c'est équivalent en C++ à x[0]
      
    }
    
  
  // cumsum des delta_tl_dyn i=1 a i...
  cumsum[0]=0;  
  for(int i = 1; i < nligne; ++i) {
    cumsum[i]=(cumsum[i-1]+delta_tl_dyn[i]); 
  }
  
  prod[0]=nppi*te[0];
  biomass[0]=prod[0]/kinetic_heat[0];
  flow_m_per_tl[0]=prod[0]/delta_tl_dyn[1];
  flow_tl_2to7[0]=(flow_m_per_tl[0]*delta_tl_dyn[1]*(-log(te[0]*(1-fishing[1]/kinetic[0])*(1-heat_tho[1]/kinetic[0]))))/
    (1-pow(te[0]*(1-fishing[1]/kinetic[0])*(1-heat_tho[1]/kinetic[0]),delta_tl_dyn[1]));
  
  // Rcout << "The value of b : " << biomass << "\n";
  
  actualize_t1.push_back(npp, "npp");
  actualize_t1.push_back(sst, "sst");
  actualize_t1.push_back(te, "te");
  actualize_t1.push_back(fishing, "fishing");
  actualize_t1.push_back(morta, "mortality");
  actualize_t1.push_back(heat_tho, "heat_tho");
  actualize_t1.push_back(MHW, "MHW");
  actualize_t1.push_back(kinetic, "kinetic");
  actualize_t1.push_back(kinetic_heat, "kinetic_heat");
  actualize_t1.push_back(delta_tl_dyn, "delta_tl_dyn");
  actualize_t1.push_back(cumsum, "cumsum");
  actualize_t1.push_back(flow_m_per_tl, "flow_m_per_tl");
  actualize_t1.push_back(flow_tl_2to7, "flow_tl_2to7");
  actualize_t1.push_back(biomass, "biomass");
  actualize_t1.push_back(time_step, "time_step");
  
  
  // ici on ajoute le vecteur my_prod comme nouvelle (et dernière) colonne de ET_T2,
  //  grace à la fonction push_back. On peut spécifier le nom, ici "production"
  return actualize_t1;
}



