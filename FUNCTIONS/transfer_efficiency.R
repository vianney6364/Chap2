# ---------------------------------------------------------------------
# function to calculate Transfer Efficiency up to 3 based on du Pontavice et al., 2019
# ---------------------------------------------------------------------
# TE is varying according to SST and ecosystem type
te_funct <- function(eco_type,sst) {
  if (eco_type == "polar"){b=0; a=0;c=0;d=0;e=0}
  if (eco_type == "temperate"){b=(0.142); a=(-0.013);c=0;d=0;e=0}
  if (eco_type == "tropical"){b=(-0.352); a=(0.015);c=0;d=0;e=0}
  if (eco_type == "upwelling"){b=(0.167); a=(-0.032);c=0;d=0;e=0}
  if (eco_type == "moyen"){b=0;a=0;c=0.05592;d=0.00325;e=0.003502}
  return(exp((-2.162-c) + (-0.025+d)*sst + b + a*sst)*(1.038013+e))
}

# te_funct <- function(eco_type,sst) {
#   if (eco_type == "polar"){b=0; a=0;c=0}
#   if (eco_type == "temperate"){b=(0.142); a=(-0.013);c=0}
#   if (eco_type == "tropical"){b=(-0.352); a=(0.015);c=0}
#   if (eco_type == "upwelling"){b=(0.167); a=(-0.032);c=0}
#   if (eco_type == "moyen"){b=0; a=0;c=0.003502}
#   
#   return(exp((-2.162) + (-0.025)*sst + b + a*sst)*(1.038013+c))
# }

# sst<-seq(-2,33,0.001)
# 
# df<-data.frame(sst)
# df_polar<-filter(df,sst<8) %>% mutate(te=te_funct("polar",sst))
# df_temperate<-filter(df,sst>=7&sst<22) %>% mutate(te=te_funct("temperate",sst))
# df_upwelling<-filter(df,sst>=8&sst<25) %>% mutate(te=te_funct("upwelling",sst))
# df_tropical<-filter(df,sst>=22&sst<33) %>% mutate(te=te_funct("tropical",sst))
# df_moyen<-df %>% mutate(te=te_funct("moyen",sst))
# 
# ggplot()+
#   geom_line(data=df_polar,aes(x = sst,y=te),color="blue")+
#   geom_line(data=df_temperate,aes(x = sst,y=te),color="goldenrod2")+
#   geom_line(data=df_upwelling,aes(x = sst,y=te),color="green2")+
#   geom_line(data=df_tropical,aes(x = sst,y=te),color="red2")+
#   geom_line(data=df_moyen,aes(x = sst,y=te),color="black")
  
  

