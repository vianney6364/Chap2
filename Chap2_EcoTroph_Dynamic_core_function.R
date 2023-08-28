########################################################################################
#creation EcoTroph_dynamic reference state with no equal length of trophic classes (CT)#
#with a same time step for passing to one CT at next one                               #
########################################################################################
# rm(list=ls()) 
#'* 1. Study framework*
#' no fishing & so no differentiation between accessible biomass or not
#' NPP fixed at 10000000
#' SST fixe at 15°C
#' Kinetic of biomass flows depending on TL an SST. (Gascuel 2008 DOI: https://doi.org/10.1016/j.ecolmodel.2008.05.012)
# i.e K(tl,SST)=20.19*(tl^(-3.26))*exp(.041*sst)
# transfert efficiency depending only to time (SST evolution)(based on du Pontavice et al., 2019 DOI: 10.1111/gcb.14944)
# i.e TE= exp((-2.162) + (-0.025)*sst + b + a*sst)*1.038013 with a and b depending of ecosystem type


#Load necessary functions : 
#load("DATA/example_data.RData")

library(Rcpp)
library(tidyverse)
library(data.table)

#' #'*Matrice evolution envi data*
source("FUNCTIONS/transfer_efficiency.R") # calulculate TE HTL
source("FUNCTIONS/dyn_reference state_add_rect_morta.R") 
# source("FUNCTIONS/dyn_reference state_add_rect_morta_by_tl.R")
# source("FUNCTIONS/reference tl_by_dtime.R")
# source("FUNCTIONS/dyn_dispatch_prod_loop.R")
# source("FUNCTIONS/dyn_dispatch_phi_per_ct.R") 
sourceCpp("FUNCTIONS/dyn_dispatch_prod_rect_morta.cpp")
sourceCpp("FUNCTIONS/tl_by_dtime_bis2_rect_morta.cpp")
sourceCpp("FUNCTIONS/compute_new_time_rcpp_rect_morta_a_fixed.cpp")
sourceCpp("FUNCTIONS/phi_pass_ct1_rcpp_rect_morta.cpp")
sourceCpp("FUNCTIONS/add_morta_biologic_coef.cpp")


#' #'*to test function execute le 4 next row*
# nb_time_step=366;dtime=1/24;eco_type="tropical";data_envi=data_envi_tropical;ratio=0.8
# nb_time_step=184;dtime=1/24;eco_type="moyen";data_envi=data_envi_tropical;sc_morta=no_morta;temperature=22
# choose ecotype between "temperate", "tropical","polar","upwelling" or "moyen" cf te_function

EcoTroph_dyn_iterative_v5_morta_biologic<-function(data_envi,nb_time_step,dtime,eco_type,clim,ratio,anomaly){ #with eco_type<-c("polar","temperate","tropical","upwelling") 
  #'*1*
  #'#'*create virtual environment data_frame*
  # data_envi<-filter(clim_ref2,clim==temperature) %>% slice(rep(1:n(), each = nb_time_step)) %>%
  # select(sst=clim,category1,category2,category3,category4) %>%
  # mutate(sst2=c(rep(temperature,20),sst[1]+category1[1],rep(temperature,40),sst[1]+category2[1],rep(temperature,40),
  #               sst[1]+category3[1],rep(temperature,40),sst[1]+category4[1],rep(temperature,40)),
  #        npp=rep(15500000,nb_time_step),
  #        time_step=seq(0,nb_time_step-1,1),
  #        mhw=c(rep(F,20),T,rep(F,40),T,rep(F,40),T,rep(F,40),T,rep(F,40))) %>%
  # select(-c(sst,category1,category2,category3,category4)) %>%
  # mutate(sst=round(sst2,1)) %>% select(-sst2) %>%
  # as.data.frame()
  #'*Reference state*
  npp_ref<-data_envi$npp[1]
  sst_ref<-data_envi$sst[1]
  mhw_ref<-data_envi$mhw[1]
  sst_year1<-data_envi$sst_year[1]
  mhw_day1<-data_envi$mhw_day[1]
  # eco_type=data_envi$eco_type[1]
  # ratio=data_envi$ratio[1]
  
  #'* 2. Creation of unequal trophic classes length for same time step to reach*
  #'*    the upper trophic classe*
  if (mhw_ref==F){
    
    morta_ref<-0  
  } else {
    morta_ref<-morta_bio_rcpp(eco_type=eco_type,sst=sst_ref,clim=clim,anomaly = anomaly,ratio=ratio)*mhw_day1
    
  }
  #Compute unequal trophic classes division
  CT_dataframe_ref_v4_rect_morta<-rcpp_tl_by_dtime2_rect_morta(sst=sst_year1,dtime=dtime,morta=morta_ref)
  # # CT_dataframe_ref2<-tl_by_dtime2(sst_ref,dtime)
  # CT_dataframe_ref_v4_morta<-CT_dataframe_ref_v4_morta %>% mutate(TLinf=lag(tl),TLsup=tl,time=round(time,3))
  # CT_dataframe_ref_v4_morta$TLinf[1]<-2.000 
  # head(CT_dataframe_ref1)
  
  # tail(CT_dataframe_ref1)
  # dim(CT_dataframe_ref1)
  #'*3. Compute production and biomass for each trophic class*
  dyn_eco_ref_fishing_v4_rect_morta<-ecotroph_core_dyn_ref_rect_morta(CT_dataframe_ref_v4_rect_morta,eco_type,dtime,npp_ref,sst_ref,sst_year1,mhw_ref,morta=morta_ref)
  # head(dyn_eco_ref_fishing)
  # dim(dyn_eco_ref_fishing)
  #------------------------------------------------------
  #STEP T to T+1 
  #----------------------------------------
  ET_ref_rect_morta<-dyn_eco_ref_fishing_v4_rect_morta%>%mutate(time_step=0,mhw_day=mhw_day1) 
  # ET<-head(ET_ref)
  # ET<-ET_ref %>% select(-production)
  # result_rect_morta<-tibble()
  result_rect_morta<-list()   # work on list and no daframe to gain a lot of time!!!
  # result_rect_morta<-rbind(ET_ref_rect_morta)
  result_rect_morta[[paste('time_step',0)]]<-ET_ref_rect_morta
  sst<-data_envi$sst[-1]
  npp<-data_envi$npp[-1]
  mhw<-data_envi$mhw[-1]
  sst_year<-data_envi$sst_year[-1]
  mhw_day<-data_envi$mhw_day[-1]
  
  
  
  #'*computing next time_step*
    for (i in 1:(nb_time_step-1)){
      
      if (mhw[i]==F){
        
        mortality<-0
        mortality2<-mortality
      } else {
        mortality<-morta_bio_rcpp(eco_type=eco_type,sst=sst[i],clim=clim,anomaly = anomaly,ratio = ratio)*mhw_day[i]
        mortality2<-mortality
      
        ssts=sst[i]

        ET_ref_rect_morta<-ET_ref_rect_morta %>% mutate(sst=ssts,
                                                        te=te_funct(eco_type,sst),
                                                        mortality=mortality2,
                                                        mhw=mhw[i],
                                                        fishing=0,
                                                        #'*with mortality*
                                                        kinetic=(20.19*(tl^(-3.26))*exp(.041*sst_year[i])+fishing),
                                                        heat_tho=(kinetic*mortality),
                                                        kinetic_heat=kinetic+heat_tho)
      }
      #2.2 creating offset of phi and compute new production value
      phi_pass_ct1_v4_rect_morta<-phi_pass_ct1_rcpp_rect_morta(ET=select(ET_ref_rect_morta,-production))
      # phi_pass_ct1_v4_morta<-ET%>%mutate(flow_actualize=lag(flow_tl_2to7)*exp(-(-log(te)+fishing/kinetic)*delta_tl_dyn),
      #                           flow_m_actualize=flow_actualize*(1-exp(-(-log(te)+fishing/kinetic)*delta_tl_dyn))/
      #                             (delta_tl_dyn*(-log(te)+fishing/kinetic)),
      #                           production=flow_m_actualize*delta_tl_dyn)%>%
      #   as.data.frame()
      # print(head(phi_pass_ct1_v4_rect_morta))
      
      #'* Compute limit of Trophic classes time t+1 : *
      #'*need to pick new time sst and npp!!!!!!!!!!*
      #Compute unequal trophic classes division at time t+1
      # dtime<-1/24 # possible to change value to see influence on biomass flow through the food web
      CT_dataframe_t1_v4_rect_morta<-rcpp_tl_by_dtime2_rect_morta(sst=sst_year[i],dtime=dtime,
                                                                  morta=mortality)
      # print(head(CT_dataframe_t1_v4_rect_morta))
      # CT_dataframe_t2<-tl_by_dtime2(sst[i],dtime)
      # CT_dataframe_t1<-CT_dataframe_t1 %>% mutate(TLinf=lag(tl),TLsup=tl,time=round(time,3))
      # CT_dataframe_t1$TLinf[1]<-2.000 
      # tail(CT_dataframe_t1_v4_morta)
      # dim(CT_dataframe_t1)
      
      #2.5 actualisation of phi by CT
      actualize_t1_rect_morta<-dyn_dispatch_prod_rect_morta_rcpp(ET_T1=phi_pass_ct1_v4_rect_morta,ET_T2=CT_dataframe_t1_v4_rect_morta)
      # print(head(actualize_t1_rect_morta))
      # dim(actualize_t1)
      
      #2.6 compute t+1 te,k,P and biomass 
      test_t1_rect_morta<-compute_new_time_rcpp_rect_morta_a_fixed(actualize_t1=actualize_t1_rect_morta,eco_type=eco_type,sst_kin=sst_year[i],
                                                                   nppi=npp[i],ssti=sst[i],time_stepi=i,dtime=dtime,mhwi=mhw[i],morta=mortality) %>% mutate(mhw_day=mhw_day[i])
      # print(head(test_t1_rect_morta))
      # result_rect_morta<-rbind(result_rect_morta,test_t1_rect_morta)
      # result_morta[[paste('time_step',i)]]<-test_t1 #add each new time_temp to the list
      result_rect_morta[[paste('time_step',i)]]<-test_t1_rect_morta
      
      ET_ref_rect_morta<-test_t1_rect_morta
    }

  result_rect_morta<-rbindlist(result_rect_morta,use.names = TRUE) #finally,transform list into dataframe data.table::rbindlist()
  return(result_rect_morta)
}



#' 
#' #' # to test function
#' output_ecotroph_dyn_iterative_temperate2<-EcoTroph_dyn_iterative("entre le num d'iteration que tu veux i.e 1500",1/24,"temperate",data_envi)
#'profvis::profvis(output_ecotroph_dyn_iterative_temperate2<-EcoTroph_dyn_iterative(10000,1/24,"temperate",data_envi))
