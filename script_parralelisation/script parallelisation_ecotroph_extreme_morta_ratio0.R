# ---------------------------------------------------------------------
rm(list=ls())
#-------------------------------------------------------------------------------
library(RPostgreSQL)
library(tidyverse)
library(parallel)
# ---------------------------------------------------------------------
source('~/CHAP2_/ET_dyn_iterative_v5_morta_biologic_a_fixe.R')
# --------------------------------------------------------------------
# linking mt database
library(rpostgis)
drv <- dbDriver("PostgreSQL")
liaisons_db_via<- dbConnect(drv, dbname="vianney_db",host="halieut.agrocampus-ouest.fr",
                            port=5432,user="vianney",password="vianneystage!")

#--------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
run_et<-function(id) {
  # NPP in mol m-2 s-1
  # SST in 
  # ---------------------------------------------------------------------
  # environmental data: see  "~/ownCloud/Thèse/CHAP2_/test_prelimilaire_envi_data_preparation.R")
  # data_envi<-filter(data,cell_no==id)
  # ---------------------------------------------------------------------
  #'*enter the wished ratio*
  ratio<-0
  #'*retrieve the biome associated to each cell  
  cell_info<-(dbGetQuery(liaisons_db_via,paste("select distinct cell_no,trunc_clim as clim,threshold,abs(threshold-trunc_clim) as anomaly, eco_type
               from daily_sst_from_avhrr.fix_yearly_baseline_threshold_clim_90th  inner join geo.biogeo3 using (cell_no)
               where cell_no=",id,"",sep="")))
  
  
  eco_type<-cell_info$eco_type
  clim<-cell_info$clim
  anomaly<-cell_info$anomaly 
  
  
  dummy_year<-dbGetQuery(liaisons_db_via,paste("
  select cell_no,year,year_month,quinzaine,sst as sst,npp_interpolated as npp from observed_forcing.fausse_annee_data_envi_quinzaine_2007
  where cell_no =",id,"",sep=""))
  
  dummy_year<-rbind(dummy_year,dummy_year,dummy_year,dummy_year,dummy_year,dummy_year,dummy_year,dummy_year,dummy_year,dummy_year,dummy_year,dummy_year)%>%
    mutate(mhw=F,mhw_day=0,time_step=seq(0,287)) %>% 
    group_by(cell_no,year) %>% mutate(sst_year=mean(sst)) #'*add mutate(mhw=F) for ratio=0*
  # head(dummy_year)
  
  data_envi_alpha0<-dbGetQuery(liaisons_db_via,paste("
  select cell_no,year,year_month,quinzaine,sst as sst,npp_interpolated as npp,mhw,percent_mhw_day as mhw_day from observed_forcing.long_term_envi_data_quinzaine
  where cell_no =",id,"",sep="")) %>% replace_na(list(mhw = FALSE)) %>% mutate(time_step=seq(288,863),mhw_day=mhw_day/16) %>% 
    group_by(cell_no,year) %>% mutate(sst_year=mean(sst)) #'*add mutate(mhw=F) for ratio=0*
  # head(data_envi_alpha0)
  
  data_envi_alpha0<-rbind(dummy_year,data_envi_alpha0)
  # ---------------------------------------------------------------------
  cell_id_alpha<- data_envi_alpha0 %>% distinct(cell_no,year,time_step)
  #-------------------------------------------------------------------------
  #'* EcoTroph-dyn*
  output_alpha<-EcoTroph_dyn_iterative_v5_morta_biologic(data_envi=data_envi_alpha0%>% mutate(mhw=F),nb_time_step = 864,dtime=1/24,clim=clim,
                                                         eco_type=eco_type,ratio=ratio,anomaly=anomaly) %>%
    inner_join(cell_id_alpha) %>% relocate(cell_no,year,time_step) %>% as.data.frame() %>% select(-c(sst,npp))
  names(output_alpha)<-tolower(names(output_alpha))
  
  #-------------------------------------------------------------------------
  output_tl_alpha<-filter(output_alpha,tl>=2,tl<=5.5)%>%
    mutate(tl=round(tl,1)) %>% 
    group_by(cell_no,time_step,tl) %>%
    # mutate(tl_moyen=tlinf+(tlsup-tlinf)/2,biomass_tl_moyen=biomass*tl_moyen) %>%
    # group_by(cell_no,time_step,quinzaine,tl) %>%
    # summarise(tl_moyen=sum(biomass_tl_moyen)/sum(biomass)
    summarise(tot_b=sum(biomass),
              tot_p=sum(production),
              te=mean(te),
              kin=mean(kinetic),
              phi=mean(flow_m_per_tl),
              kin_heat=mean(kinetic_heat),
              mortality=mean(mortality),
              heat_tho=mean(heat_tho),
              mhw=mhw[1])
  # mhw_day=mhw_day[1])
  dbWriteTable(liaisons_db_via,c("ecotroph_output","ecotroph_extreme_morta_ratio0"),output_tl_alpha, row.names=FALSE, append=TRUE) # append = T pour ajouter des lignes
  
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
id_cell_no <- as.data.frame(dbGetQuery(liaisons_db_via,"select distinct cell_no from observed_forcing.long_term_envi_data_quinzaine")) 
dim(id_cell_no)

system.time({
  cl <- makeCluster(60) 
  clusterExport(cl, c('compute_new_time_rcpp_rect_morta_a_fixed','dyn_dispatch_prod_rect_morta_rcpp','ecotroph_core_dyn_ref_rect_morta',
                      'morta_bio_rcpp','te_funct','phi_pass_ct1_rcpp_rect_morta','rcpp_tl_by_dtime2_rect_morta','EcoTroph_dyn_iterative_v5_morta_biologic'))
  
  clusterEvalQ(cl, {  # ici la boucle et pour faire charger d'ans l'environnement des 20 R ce qui est nécessaire pour faire tourner EcoTroph
    library(Rcpp)
    library(RPostgreSQL)
    library(data.table)
    library(tidyverse)
    # linking mt database
    drv <- dbDriver("PostgreSQL")
    liaisons_db_via<- dbConnect(drv, dbname="vianney_db",host="halieut.agrocampus-ouest.fr",
                                port=5432,user="vianney",password="vianneystage!")
    source('ET_dyn_iterative_v5_morta_biologic_a_fixe.R')
    sourceCpp("rcpp function add morta/dyn_dispatch_prod_rect_morta.cpp")
    sourceCpp("rcpp function add morta/tl_by_dtime_bis2_rect_morta.cpp")
    sourceCpp("rcpp function add morta/compute_new_time_rcpp_rect_morta_a_fixed.cpp")
    sourceCpp("rcpp function add morta/phi_pass_ct1_rcpp_rect_morta.cpp")
    sourceCpp("rcpp function add morta/add_morta_biologic_old.cpp")
    
  })
  values <-id_cell_no[,1] # pour récupérer l'ensemble des cellules dans un vecteur
  result <- parLapply(cl,values,run_et)   # paralel execution
  clusterEvalQ(cl, {
    dbDisconnect(liaisons_db_via)
  })
  stopCluster(cl)
})


