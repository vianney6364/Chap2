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

// [[Rcpp::export]]
//Model morta biologique et morta50=cat4


// model by cat
double morta_bio_rcpp(String eco_type, double ssti, double clim, double anomaly, double ratio)  {
  double b;
  double d;
  double e;
  double c;
  double mhw_cat;
  double affected_value;
  mhw_cat=(ssti-clim)/anomaly; b=-3.09+0.474*clim-0.0224*pow(clim,2);c=0; d=1;e=3.09-0.0938*clim; //stable value from ~/ownCloud/Thèse/Script/~/ownCloud/Thèse/Script/SCRIPT EcoTroph/script ecotroph dynamique/Virtual dynamique EcoTroph GOMPERTZ morta/morta uniform food web/data/creation mean_anomaly by clim temperature by biome_RData.R
  // and filter(df,eco_type=="polar")$percentage[1]/100
  affected_value= c+(d-c)*exp(-exp(b*ratio*(mhw_cat-e/ratio)));
  return affected_value;
}
