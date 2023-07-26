library(Rcpp)
library(RcppEigen)
library(RobustGaSP)
library(L1pack)
library(gnorm)

apolar_vicsek=T 

source('functions/function_cell.R')

file_number=441 ## please choose a number 
delete_first_2_hours=T
unsorted_frameID=F
unsorted_particleID=F
##
file_directory=paste0('real_data/',file_number,'-velocity_pos_record.csv')
velocity_loc_frame_record=read.csv(file_directory,header=T)

attach(velocity_loc_frame_record)
head(velocity_loc_frame_record)



dim(velocity_loc_frame_record)

mean(abs(velocity_loc_frame_record$vx))
mean(abs(velocity_loc_frame_record$vy))

##1.1 one particle cannot have two locations, so delete duplication 
for(i in unique(velocity_loc_frame_record$frameID)){
  index_selected=which(velocity_loc_frame_record$frameID==i)
  p_here=cbind(velocity_loc_frame_record$px[index_selected],velocity_loc_frame_record$py[index_selected])  ##
  v_here=cbind(velocity_loc_frame_record$vx[index_selected],velocity_loc_frame_record$vy[index_selected])  ##
  
  delete_index=delete_duplicate_pos(p_here,v_here)
  if(length(delete_index)>0){
    velocity_loc_frame_record=velocity_loc_frame_record[-index_selected[delete_index],]
  }
  
}

unique(velocity_loc_frame_record$frame)
# 1.2 the first few has tracking problems for 441
if(delete_first_2_hours){
  index_row_bad_tracking=which(velocity_loc_frame_record$frameID<6)
  ###delete the number in this crack
  velocity_loc_frame_record=velocity_loc_frame_record[-index_row_bad_tracking,]
  velocity_loc_frame_record$frameID=velocity_loc_frame_record$frameID-6 ###change this back to the original frame
}
if(file_number==161){
  # #the below is only for 161, delete the crack
  index_row_crack=which(velocity_loc_frame_record$py>1180&velocity_loc_frame_record&px<1000)
  ##delete the number in this crack
  velocity_loc_frame_record=velocity_loc_frame_record[-index_row_crack,]
  delta_t=1 ##is it?
  index_plot=which(velocity_loc_frame_record[,5]==3) ##some problem when py is close to 1000
  plot(velocity_loc_frame_record[index_plot,3],velocity_loc_frame_record[index_plot,4],
       cex=.2,pch=20,xlab=expression(p[x]),ylab=expression(p[y]))
  arrows(x0=velocity_loc_frame_record[index_plot,3],y0=velocity_loc_frame_record[index_plot,4],
         x1=velocity_loc_frame_record[index_plot,3]+velocity_loc_frame_record[index_plot,1]*delta_t*10000,
         y1=velocity_loc_frame_record[index_plot,4]+velocity_loc_frame_record[index_plot,2]*delta_t*10000,
         length=0.03
  )
  
}
###sort frameID if not
if(unsorted_frameID){
  velocity_loc_frame_record_unsorted_frameID=velocity_loc_frame_record
  frameID_unique= unique(velocity_loc_frame_record$frameID)
  count=0
  for(i in frameID_unique){
    index_frameID_here=which(velocity_loc_frame_record_unsorted_frameID$frameID==i)
    velocity_loc_frame_record[count+(1:length(index_frameID_here)),]=  velocity_loc_frame_record_unsorted_frameID[index_frameID_here,]
    count=count+length(index_frameID_here)
  }
  rm(velocity_loc_frame_record_unsorted_frameID)
  
}
###sort particleID if not
if(unsorted_particleID){
  frameID_unique= unique(velocity_loc_frame_record$frameID)
  for(i in frameID_unique){
    index_frameID_here=which(velocity_loc_frame_record$frameID==i)
    sort_particleID_here_all=sort(velocity_loc_frame_record$particleID[index_frameID_here],index.return=T)
    velocity_loc_frame_record[index_frameID_here,]=velocity_loc_frame_record[index_frameID_here[sort_particleID_here_all$ix],]
  }
}

index_plot=which(velocity_loc_frame_record$frame==20) ##some problem when py is close to 1000

plot(velocity_loc_frame_record[index_plot,]$px,velocity_loc_frame_record[index_plot,]$py,
     cex=.2,pch=20,xlab=expression(p[x]),ylab=expression(p[y]))
arrows(x0=velocity_loc_frame_record[index_plot,]$px,y0=velocity_loc_frame_record[index_plot,]$py,
       x1=velocity_loc_frame_record[index_plot,]$px+velocity_loc_frame_record[index_plot,]$vx*10,
       y1=velocity_loc_frame_record[index_plot,]$py+velocity_loc_frame_record[index_plot,]$vy*10,
       length=0.03
)


# 
S=1
D=2
T_time=max(velocity_loc_frame_record$frameID)  # time steps


##2. put  position and velocity into lists by time; form pairs
pos_v_theta_all_list=get_pos_v_theta_all_vec(velocity_loc_frame_record,T_time)
v_all=pos_v_theta_all_list$v_all
v_all_end=pos_v_theta_all_list$v_all_end   ##end 
theta_all=pos_v_theta_all_list$theta_all ###alignment vectors 

pos_all_vec=pos_v_theta_all_list$pos_all_vec ##position vectors
n_record=pos_v_theta_all_list$n_record   ###number of cell at each frame

vx_end_all=v_all_end[seq(1,length(v_all_end),2)]
vy_end_all=v_all_end[seq(2,length(v_all_end),2)]

vx_all=v_all[seq(1,length(v_all),2)]
vy_all=v_all[seq(2,length(v_all),2)]

theta_end_all=atan2(vy_end_all,vx_end_all)

pos_x_all=pos_all_vec[seq(1,length(pos_all_vec),2)]
pos_y_all=pos_all_vec[seq(2,length(pos_all_vec),2)]

sigma_2_x=rep(NA,T_time)
sigma_2_y=rep(NA,T_time)
vel_x_end_list=as.list(1:T_time)
vel_y_end_list=as.list(1:T_time)
vel_x_list=as.list(1:T_time)
vel_y_list=as.list(1:T_time)

pos_x_list=as.list(1:T_time)
pos_y_list=as.list(1:T_time)

theta_list=as.list(1:T_time)
theta_end_list=as.list(1:T_time)

count_sum=0
theta_0_pi_list=as.list(1:T_time)
theta_end_0_pi_list=as.list(1:T_time)

##new look at this l1 loss 
mean_abs_loss_x=rep(NA,T_time)
mean_abs_loss_y=rep(NA,T_time)

for(t in 1:T_time){
  sigma_2_x[t]=var(vx_all[count_sum+1:n_record[t]])
  sigma_2_y[t]=var(vy_all[count_sum+1:n_record[t]])
  
  ##
  mean_abs_loss_x[t]=mean(abs(vx_all[count_sum+1:n_record[t]]-mean(vx_all[count_sum+1:n_record[t]]) ))
  mean_abs_loss_y[t]=mean(abs(vy_all[count_sum+1:n_record[t]]-mean(vy_all[count_sum+1:n_record[t]]) ))
  
  vel_x_list[[t]]=(vx_all[count_sum+1:n_record[t]])
  vel_y_list[[t]]=(vy_all[count_sum+1:n_record[t]])
  
  vel_x_end_list[[t]]=(vx_end_all[count_sum+1:n_record[t]])
  vel_y_end_list[[t]]=(vy_end_all[count_sum+1:n_record[t]])
  
  
  pos_x_list[[t]]=(pos_x_all[count_sum+1:n_record[t]])
  pos_y_list[[t]]=(pos_y_all[count_sum+1:n_record[t]])
  
  theta_list[[t]]=(theta_all[count_sum+1:n_record[t]])
  theta_0_pi_list[[t]]= theta_list[[t]]
  theta_0_pi_list[[t]][which(theta_list[[t]]<0)]= theta_0_pi_list[[t]][which(theta_list[[t]]<0)]+pi
  
  theta_end_list[[t]]=(theta_end_all[count_sum+1:n_record[t]])
  theta_end_0_pi_list[[t]]= theta_end_list[[t]]
  theta_end_0_pi_list[[t]][which(theta_end_list[[t]]<0)]= theta_end_0_pi_list[[t]][which(theta_end_list[[t]]<0)]+pi
  
  
  count_sum=count_sum+n_record[t]
  
}
#
order_param=rep(NA,T_time)

for(t in 1:T_time){
  order_param[t]=sum(cos(2*theta_0_pi_list[[t]]))/n_record[t]
}

##save cell counts
csv_directory=paste0('csv_variables/cell_count_real_',file_number,'.csv')
write.csv(n_record,file=csv_directory,row.names=F)

##save experimental obs 
csv_directory=paste0('csv_variables/order_param_real_',file_number,'.csv')
write.csv(order_param,file=csv_directory,row.names=F)

csv_directory=paste0('csv_variables/sigma_2_x_real_',file_number,'.csv')
write.csv(sigma_2_x,file=csv_directory,row.names=F)

csv_directory=paste0('csv_variables/sigma_2_y_real_',file_number,'.csv')
write.csv(sigma_2_y,file=csv_directory,row.names=F)

csv_directory=paste0('csv_variables/mean_abs_loss_x_real_',file_number,'.csv')
write.csv(mean_abs_loss_x,file=csv_directory,row.names=F)

csv_directory=paste0('csv_variables/mean_abs_loss_y_real_',file_number,'.csv')
write.csv(mean_abs_loss_y,file=csv_directory,row.names=F)


##trend of sd
plot(sqrt(sigma_2_x))
plot(sqrt(sigma_2_y))

px_min=min(velocity_loc_frame_record[,3])
px_max=max(velocity_loc_frame_record[,3])
py_min=min(velocity_loc_frame_record[,4])
py_max=max(velocity_loc_frame_record[,4])


nx=30
ny=20

grid_boundary_info=get_boundary_grid(px_min=px_min,px_max=px_max,
                                     py_min=py_min,py_max=py_max,nx=nx,ny=ny)

cut_r_max=100


T_sim=T_time
h=1200
##no duplication now
p0=cbind(pos_x_list[[1]],pos_y_list[[1]])  ##
v0=cbind(vel_x_list[[1]],vel_y_list[[1]])
n_sim=dim(p0)[1]


#####36 simulations
#####only 1 repetitation  
###save orientational order 
apolar_viscek_vec = c(T,F)
fixed_weight_vec=c(F,T)
residual_type_vec=c('Gaussian','Laplace','GGD') 
cut_r_vec = c(25, 50, 75)
#sigma_0=.5
#settings for all experiment
sim_set = expand.grid(apolar_viscek=apolar_viscek_vec,fixed_weight=fixed_weight_vec,residual_type=residual_type_vec,cut_r=cut_r_vec)
#n_repeat=1 ##only repetition first 

order_param_rec = matrix(NA,dim(sim_set)[1],T_sim,
                     dimnames = list(paste("apolar_viscek=",sim_set[,1],",fixed_weight=",sim_set[,2],",residual_type=",
                                           sim_set[,3],",cut_r=",  sim_set[,4],sep=""),
                                     paste("T_sim",1:T_sim,sep = "")))

sigma_2_x_rec =  matrix(NA,dim(sim_set)[1],T_sim,
                        dimnames = list(paste("apolar_viscek=",sim_set[,1],",fixed_weight=",sim_set[,2],",residual_type=",
                                              sim_set[,3],",cut_r=",  sim_set[,4],sep=""),
                                        paste("T_sim",1:T_sim,sep = "")))

sigma_2_y_rec = matrix(NA,dim(sim_set)[1],T_sim,
                       dimnames = list(paste("apolar_viscek=",sim_set[,1],",fixed_weight=",sim_set[,2],",residual_type=",
                                             sim_set[,3],",cut_r=",  sim_set[,4],sep=""),
                                       paste("T_sim",1:T_sim,sep = "")))

mean_abs_loss_x_rec =  matrix(NA,dim(sim_set)[1],T_sim,
                        dimnames = list(paste("apolar_viscek=",sim_set[,1],",fixed_weight=",sim_set[,2],",residual_type=",
                                              sim_set[,3],",cut_r=",  sim_set[,4],sep=""),
                                        paste("T_sim",1:T_sim,sep = "")))

mean_abs_loss_y_rec = matrix(NA,dim(sim_set)[1],T_sim,
                       dimnames = list(paste("apolar_viscek=",sim_set[,1],",fixed_weight=",sim_set[,2],",residual_type=",
                                             sim_set[,3],",cut_r=",  sim_set[,4],sep=""),
                                       paste("T_sim",1:T_sim,sep = "")))


#fixe_time_forecast=T
#fixed_time_frame=20
for(sim in 1:dim(sim_set)[1]){
  
  print(sim)
  apolar_vicsek = sim_set[sim,1]
  fixed_weight = sim_set[sim,2]
  residual_type = sim_set[sim,3]
  cut_r = sim_set[sim,4]

  set.seed((sim))
  
  ##1. find maximum neighbor
  max_neighbors_list=find_max_neighbors_fast_grid(pos_x_list,pos_y_list, vel_x_list,vel_y_list,n_record,T_time,
                                                  grid_boundary_info,cut_r_max=cut_r_max,
                                                  apolar_vicsek=apolar_vicsek)
  
  v_max_neighbor_x_record=max_neighbors_list$v_neighbor_x_record
  v_max_neighbor_y_record=max_neighbors_list$v_neighbor_y_record
  
  d_pos_max_vec=max_neighbors_list$d_pos_vec
  num_neighbors_max_vec=max_neighbors_list$num_neighbors_vec
  
  ##2. param estimation
  param_est=get_param_est_with_max_neighbor(vel_x_list,vel_y_list,vel_x_end_list,vel_y_end_list,n_record,T_time,
                                            v_max_neighbor_x_record,v_max_neighbor_y_record,
                                            d_pos_max_vec,num_neighbors_max_vec,cut_r=cut_r,
                                            fixed_weight=fixed_weight,residual_type=residual_type)
  
  sigma_2_0_x_hat=param_est$sigma_2_0_x_hat
  sigma_2_0_y_hat=param_est$sigma_2_0_y_hat
  wx_hat=param_est$wx_hat
  wy_hat=param_est$wy_hat
  
  alpha_x_hat=NULL
  alpha_y_hat=NULL
  beta_x_hat=NULL
  beta_y_hat=NULL
  if(residual_type=='GGD'){
    alpha_x_hat=param_est$alpha_x_hat
    alpha_y_hat=param_est$alpha_y_hat
    beta_x_hat=param_est$beta_x_hat
    beta_y_hat=param_est$beta_y_hat
  }
  

  ##3. simulation 
  m_simulation=simulate_cells_anisotropic_with_neighbor_fast_grid(cut_r,p0,v0,sigma_2_0_x_hat,sigma_2_0_y_hat,
                                                                  T_sim,h,wx_hat,wy_hat,grid_boundary_info,
                                                                  apolar_vicsek=apolar_vicsek,
                                                                  alpha_x_hat=alpha_x_hat,alpha_y_hat=alpha_y_hat,beta_x_hat=beta_x_hat,beta_y_hat=beta_y_hat,
                                                                  residual_type=residual_type)
  order_param_rec[sim,]=m_simulation$record_ensemble_order

  sigma_2_x_sim=rep(NA,T_time)
  sigma_2_y_sim=rep(NA,T_time)
  
  mean_abs_loss_x_sim=rep(NA,T_time)
  mean_abs_loss_y_sim=rep(NA,T_time)
  
  for(i_time in 1:T_time){
    sigma_2_x_sim[i_time]=var(m_simulation$vel_x_sim_record[,i_time])
    sigma_2_y_sim[i_time]=var(m_simulation$vel_y_sim_record[,i_time])
    
    mean_abs_loss_x_sim[i_time]=mean(abs(m_simulation$vel_x_sim_record[,i_time]-mean(m_simulation$vel_x_sim_record[,i_time]) ))
    mean_abs_loss_y_sim[i_time]=mean(abs(m_simulation$vel_y_sim_record[,i_time]-mean(m_simulation$vel_y_sim_record[,i_time]) ))
    
  }
  
  sigma_2_x_rec[sim,]=sigma_2_x_sim
  sigma_2_y_rec[sim,]=sigma_2_y_sim
  
  mean_abs_loss_x_rec[sim,]=mean_abs_loss_x_sim
  mean_abs_loss_y_rec[sim,]=mean_abs_loss_y_sim
  
  
}


write.csv(sim_set,file="csv_variables/simul_setting.csv",row.names = F)

csv_directory=paste0('csv_variables/order_param_rec_',file_number,'.csv')
write.csv(order_param_rec,file=csv_directory)

csv_directory=paste0('csv_variables/mean_abs_loss_x_rec_',file_number,'.csv')
write.csv(mean_abs_loss_x_rec,file=csv_directory)

csv_directory=paste0('csv_variables/mean_abs_loss_y_rec_',file_number,'.csv')
write.csv(mean_abs_loss_y_rec,file=csv_directory)

# csv_directory=paste0('csv_variables/sigma_2_x_rec_',file_number,'.csv')
# write.csv(sigma_2_x_rec,file=csv_directory)
# 
# csv_directory=paste0('csv_variables/sigma_2_y_rec_',file_number,'.csv')
# write.csv(sigma_2_y_rec,file=csv_directory)

