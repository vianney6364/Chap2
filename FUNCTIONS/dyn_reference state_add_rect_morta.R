# ---------------------------------------------------------------------
# ECoTroph Core function 
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# Tl= Trophic Level
# TE: Transfer efficiency = TE HTL: higher trophic levels transfer efficiency (du pontavice et al., 2019 DOI: https://doi.org/10.1111/gcb.14944)
# NPP: Net Primary Production
# SST: Sea Surface Temperature
# ---------------------------------------------------------------------
# load required library
# install.packages("tidyverse")
# library(tidyverse)
# library(dplyr)
# ---------------------------------------------------------------------
# load several fonction 
# source("FUNCTIONS/flow_kinetic.R") # calculate the flow kinetic
# source("FUNCTIONS/transfer_efficiency.R") # calulculate TE HTL
# source("FUNCTIONS/dyn_production by trophic class geometry sequence.R") #process to calculate consumer production between TL=2 and TL=7
# geometric sequence to speed up the process to calculate consumer production
# ---------------------------------------------------------------------
# 1.ref dataframe first time step (one single year)
# 2.sst=15 (in °C)
# 3.npp (in mole of Carbon .m-2.s-1)set to 10000000
# eco_type : 4 choices of ecosystem types: "polar","temperate","tropical","upwelling" # the sensivty of te_htl to temperature is varying among ecosystem types (du pontavice et al., 2019 DOI: https://doi.org/10.1111/gcb.14944)
# surface_area: Surface of the ecosystem in m²
#----------------
#' #'*Enter parameter for ecotroph_core function
# data_envi_rect=CT_dataframe_ref_v4_rect_morta
# eco_type
# dtimenpp_ref,sst_ref,mhw_ref,morta=morta_ref
#' 
#-------------------------

ecotroph_core_dyn_ref_rect_morta<-function(data_envi_rect,eco_type,dtime,npp_ref,sst_ref,sst_year,mhw_ref,morta) { 
  # if (mhw_ref==F){
  data_envi_rect<-data_envi_rect%>%
    mutate(npp=npp_ref,sst=sst_ref,
           te=te_funct(eco_type,sst),
           mortality=morta,
           MHW=mhw_ref,
           fishing=0,
           #'*without mortality*
           # kinetic=(20.19*(tl^(-3.26))*exp(.041*sst))*dtime,
           #'*with mortality*
           kinetic=(20.19*(tl^(-3.26))*exp(.041*sst_year)+fishing)*dtime,
           heat_tho=(kinetic*mortality),
           # heat_tho=mortality,
           kinetic_heat=kinetic+heat_tho*dtime,
           delta_tl_dyn=TLsup-TLinf,
           cumsum=cumsum(delta_tl_dyn))
  #' } else {
  #'   data_envi_rect<-data_envi_rect%>%
  #'     mutate(npp=npp_ref,sst=sst_ref,te=te_funct(eco_type,sst),
  #'            # LT50=4+filter(event_biome,biome==eco_type)$clim,
  #'            # mortality=0+(1-0)*exp(-(exp(-1*(sst-LT50)))),
  #'            mortality=morta,
  #'            MHW=T,
  #'            fishing=0,
  #'            #'*without mortality*
  #'            # kinetic=(20.19*(tl^(-3.26))*exp(.041*sst))*dtime,
  #'            #'*with mortality*
  #'            kinetic=(20.19*(tl^(-3.26))*exp(.041*sst)+fishing)*dtime,
  #'            heat_tho=(kinetic*mortality),
  #'            kinetic_heat=kinetic+heat_tho,
  #'            delta_tl_dyn=TLsup-TLinf,
  #'            cumsum=cumsum(delta_tl_dyn))
  #'     }
  ##---------------------------------------------------------------------
  # 1st step: Calculate production from TL=2 to TL=7
  # ---------------------------------------------------------------------
  # Production at TL = 2: because by convention prod2/NPP = TE LTL
  data_envi_rect$prod2<- data_envi_rect$npp* data_envi_rect$te
  # mean biomass flow between 2 and tl sup
  data_envi_rect$flow_m2<-data_envi_rect$prod2/(data_envi_rect$delta_tl_dyn[2])
  #'*biomass flow at CT=1 without add mortality and fishing*
  # data_envi_rect$flow2<-data_envi_rect$flow_m2[1]*(data_envi_rect$delta_tl_dyn[2]*(-log(data_envi_rect$te)))/(1-data_envi_rect$te^(data_envi_rect$delta_tl_dyn[2]))
  #'*biomass flow at CT=1 with add mortality and fishing*
  data_envi_rect<-data_envi_rect%>%
    mutate(flow2=flow_m2[1]*delta_tl_dyn[2]*(-log(te*(1-heat_tho[1]/kinetic[1])*(1-fishing[1]/kinetic[1])))/
             (1-(te*(1-heat_tho[1]/kinetic[1])*(1-fishing[1]/kinetic[1]))**delta_tl_dyn[2]))
  # # production between TL=2 and TL=7 and mean biomass flow between TL=2 and TL=3
  #'*without add mortality and fishing*
  # data_envi_rect<-data_envi_rect%>%
  #   mutate(flow_tl_2to7=flow2*(te^cumsum),
  #          flow_m_per_tl=flow_tl_2to7*(1-(te^delta_tl_dyn))/(delta_tl_dyn*(-log(te))),
  #          production=flow_m_per_tl*delta_tl_dyn)%>%
  #   # select(-flow_tl_2to7)%>%
  #   as.data.frame()
  #'*with add mortality and fishing*
  data_envi_rect<-data_envi_rect%>%
    mutate(flow_tl_2to7=flow2*(te*(1-heat_tho[1]/kinetic[1])*(1-fishing[1]/kinetic[1]))**cumsum,
           flow_m_per_tl=flow_tl_2to7*(1-(te*(1-heat_tho[1]/kinetic[1])*(1-fishing[1]/kinetic[1]))**delta_tl_dyn)/
             (delta_tl_dyn*(-log(te*(1-heat_tho[1]/kinetic[1])*(1-fishing[1]/kinetic[1])))),
           production=flow_m_per_tl*delta_tl_dyn)%>%
    # mutate(flow_tl_2to7=flow2*exp(-(-log(te)+(heat_tho/kinetic)+fishing/kinetic)*cumsum),
    #        flow_m_per_tl=flow_tl_2to7*(1-exp(-(-log(te)+(heat_tho/kinetic)+fishing/kinetic)*delta_tl_dyn))/
    #          (delta_tl_dyn*(-log(te)+(heat_tho/kinetic)+fishing/kinetic)),
    #        production=flow_m_per_tl*delta_tl_dyn)%>%
    # select(-flow_tl_2to7)%>%
    as.data.frame()
  #  production
  data_envi_rect$production[1]<-data_envi_rect$prod2[1]
  data_envi_rect$flow_m_per_tl[1]<-data_envi_rect$flow_m2[1]
  
  # approximation : dernier TL ~environ egal au precedent : au niveau 7 l'erreur devient assez negligeable cf Paul Gatti work
  #' data_envi[length(data_envi$tl),"production"]<-data_envi[length(data_envi$tl)-1,"production"]
  #' #'[Check that with Didier, & try to understand why issue for last line of production ]
  #' # Biomass : Production(tl)/Kinetic(tl)
  data_envi_rect$biomass<-(data_envi_rect$production/data_envi_rect$kinetic_heat)
  result<-data_envi_rect%>%dplyr::select(-prod2,-flow2,-flow_m2)
  return(result)
}

#------------------------
##offset function a pour but de faire un décalage de phi entre classe trophic
#mais toujours avec le même découpage : 1ere étape pour passer à t+1
#to test : 
# data_envi<-dyn_data_ref_mortality
#' offset<-function(data_envi) { 
#'   # data_envi$npp<-npp  #new time npp
#'   #data_envi$flow_m_per_tl[1]<-(data_envi$npp[1]* data_envi$te[1])/(data_envi$delta_tl_dyn[2]) #new phi2 [(new time npp)*(new time te)]/delta.tl_dyn[2] (always equal to 1)
#'   # data_envi$flow_tl_2to7[2:(length(data_envi$flow_tl_2to7))]<-data_envi$flow_tl_2to7[-length(data_envi$flow_tl_2to7)]*exp(-(-log(data_envi$te[-length(data_envi$te)])+(data_envi$mortality[-length(data_envi$mortality)]/data_envi$kinetic[-length(data_envi$kinetic)]))*data_envi$delta_tl_dyn[-length(data_envi$delta_tl_dyn)]) #réalisation du décallage de phi
#'   data_envi<-data_envi%>%mutate(flow_actualize=lag(flow_tl_2to7,1),
#'                        flow_m_actualize=flow_actualize*(1-exp(-(-log(te)+mortality/kinetic))*delta_tl_dyn)/
#'                   (delta_tl_dyn*(-log(te)+mortality/kinetic)),
#'                   production=flow_m_actualize*delta_tl_dyn)%>%
#'     select(TLinf,TLsup,tl,production)%>%as.data.frame()
#'        
#'   # data_envi$flow_tl_2to7[2:(length(data_envi$flow_tl_2to7))]<-data_envi$flow_tl_2to7[-length(data_envi$flow_tl_2to7)]
#'   #'*passer en dplyr et Pipe*
#'   # data_envi$flow_m_per_tl[2:(length(data_envi$flow_m_per_tl))]<-data_envi$flow_tl_2to7[2:(length(data_envi$flow_tl_2to7))]*(1-exp(-(-log(data_envi$te[2:(length(data_envi$te))])+(data_envi$mortality[2:(length(data_envi$mortality))]/data_envi$kinetic[2:(length(data_envi$kinetic))]))*data_envi$delta_tl_dyn[2:(length(data_envi$delta_tl_dyn))]))/(data_envi$delta_tl_dyn[2:(length(data_envi$delta_tl_dyn))]*(-log(data_envi$te[2:(length(data_envi$te))])+(data_envi$mortality[2:(length(data_envi$mortality))]/data_envi$kinetic[2:(length(data_envi$kinetic))])))
#'   # data_envi$production[2:(length(data_envi$production))]<-data_envi$flow_m_per_tl[2:(length(data_envi$flow_m_per_tl))]*data_envi$delta_tl_dyn[2:(length(data_envi$flow_m_per_tl))]
#'   #réalisation du décallage de phi
#'   # data_envi$flow_m_per_tl[2:(length(data_envi$flow_m_per_tl))]<-data_envi$flow_m_per_tl[-length(data_envi$flow_m_per_tl)] #réalisation du décallage de phi
#'   result<-data_envi%>%select(TLinf,TLsup,tl,production)%>%as.data.frame()
#'   return(result)
#' }

