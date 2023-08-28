#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
DataFrame phi_pass_ct1_rcpp_rect_morta(DataFrame ET) {
  
  NumericVector flow_tl_2to7 = ET["flow_tl_2to7"];
  NumericVector te = ET["te"];
  NumericVector fishing = ET["fishing"];
  NumericVector morta = ET["mortality"];
  NumericVector heat_tho = ET["heat_tho"];
  NumericVector kinetic = ET["kinetic"];
  NumericVector kinetic_heat = ET["kinetic_heat"];
  NumericVector delta_tl_dyn = ET["delta_tl_dyn"];
  int nligne = ET.nrows();
  
  
  NumericVector lag_flow_tl_2to7(nligne); 
  NumericVector flow_actualize(nligne); 
  NumericVector flow_m_actualize(nligne); 
  NumericVector prod(nligne); 
  
  lag_flow_tl_2to7[0]=R_NaN; 
  for(int i = 1; i < te.size(); ++i) {
    lag_flow_tl_2to7[i]=flow_tl_2to7[i-1];
  }

  for(int j = 0; j < te.size(); ++j) {
    flow_actualize[j]=lag_flow_tl_2to7[j]*pow(te[j]*(1-fishing[j]/kinetic[j])*(1-heat_tho[j]/kinetic[j]),delta_tl_dyn[j]);
    flow_m_actualize[j]=flow_actualize[j]*(1-pow(te[j]*(1-fishing[j]/kinetic[j])*(1-heat_tho[j]/kinetic[j]),delta_tl_dyn[j]))/
      (delta_tl_dyn[j]*(-log(te[j]*(1-fishing[j]/kinetic[j])*(1-heat_tho[j]/kinetic[j]))));
    prod[j]=flow_m_actualize[j]*delta_tl_dyn[j];
    // -log(te*(1-mortality[1]/kinetic[1])*(1-fishing[1]/kinetic[1]))
    // flow_actualize[j]=lag_flow_tl_2to7[j]*exp(-(-log(te[j])+fishing[j]/kinetic[j]+morta[j]/kinetic[j])*delta_tl_dyn[j]);
    // flow_m_actualize[j]=flow_actualize[j]*(1-exp(-(-log(te[j])+fishing[j]/kinetic[j]+morta[j]/kinetic[j])*delta_tl_dyn[j]))/
    //   (delta_tl_dyn[j]*(-log(te[j])+fishing[j]/kinetic[j]+morta[j]/kinetic[j]));
    // prod[j]=flow_m_actualize[j]*delta_tl_dyn[j];
  }

  ET.push_back(lag_flow_tl_2to7, "lag");
  ET.push_back(flow_actualize, "flow_actualize");
  ET.push_back(flow_m_actualize, "flow_m_actualize");
  ET.push_back(prod, "production");

  return ET;
  
}

