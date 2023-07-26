library(Rcpp)
library(RcppEigen)
library(RobustGaSP)
library(gnorm)
library(L1pack)

source('functions/function_cell.R')

######################################
############ load dataset ############
######################################
file_number=441 ## please choose a number, 
delete_first_2_hours=T  
unsorted_frameID=F
unsorted_particleID=F
interval_est=F
run_permutation_F_test=F
run_test_neighbor=F
run_fit_marginal_velocity=F
run_fit_conditional_velocity_simulation=F

##load dataset
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

T_time=max(velocity_loc_frame_record$frameID)  # time steps

S=1
D=2


##2. put position and velocity into lists by time; form pairs
pos_v_theta_all_list=get_pos_v_theta_all_vec(velocity_loc_frame_record,T_time)
v_all=pos_v_theta_all_list$v_all
v_all_end=pos_v_theta_all_list$v_all_end   ##end 
theta_all=pos_v_theta_all_list$theta_all ###alignment vectors 

pos_all_vec=pos_v_theta_all_list$pos_all_vec ##position vectors
n_record=pos_v_theta_all_list$n_record   ###number of cell at each frame

#separate position and velocity for two directions
vx_end_all=v_all_end[seq(1,length(v_all_end),2)]
vy_end_all=v_all_end[seq(2,length(v_all_end),2)]

vx_all=v_all[seq(1,length(v_all),2)]
vy_all=v_all[seq(2,length(v_all),2)]

pos_x_all=pos_all_vec[seq(1,length(pos_all_vec),2)]
pos_y_all=pos_all_vec[seq(2,length(pos_all_vec),2)]

theta_end_all=atan2(vy_end_all,vx_end_all)

#put all variables into list
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

index_frame=60
data_time=data.frame(vel_x=vel_x_list[[index_frame]],vel_y=vel_y_list[[index_frame]],
                     vel_x_end=vel_x_end_list[[index_frame]],vel_y_end=vel_y_end_list[[index_frame]],
                     pos_x=pos_x_list[[index_frame]],pos_y=pos_y_list[[index_frame]],
                     theta=theta_list[[index_frame]],theta_0_pi=theta_0_pi_list[[index_frame]],
                     theta_end=theta_end_list[[index_frame]],theta_end_0_pi=theta_end_0_pi_list[[index_frame]])
write.csv(data_time,file=paste0('csv_variables/data_at_time_',index_frame,'_exp_',file_number,'.csv'),row.names = F)

## use Guassian, Laplace, or GGD to fit the marginal velocity at time frame=60
if(run_fit_marginal_velocity){
  ##one should let others to choose 
  
  index_frame=60
  par(mfrow=c(1,2))
  
  #x-diraction
  density_x=density(vel_x_end_list[[index_frame]])
  plot(density_x$x,log10(density_x$y),ylim=c(0,log10(110)),type='l')
  mean_x_norm=mean(vel_x_end_list[[index_frame]])
  sd_x_norm=sd(vel_x_end_list[[index_frame]])
  ##best fit normal 
  density_norm_x=dnorm(density_x$x,mean=mean_x_norm,sd=sd_x_norm)
  lines(density_x$x,log10(density_norm_x),lty=4,col='orange')
  ##best fit laplace
  m_laplace_vx=lad(vel_x_end_list[[index_frame]]~1)
  density_laplace_x=dlaplace(density_x$x, location = m_laplace_vx$coefficients, 
                             scale = m_laplace_vx$scale, log = FALSE)
  lines(density_x$x,log10(density_laplace_x),col='red')
  ##best fit GGD
  param=c(1.9,.5)
  m_GGD_x=optim(param,fn=GGD_log_likelihood_profile,gr=GGD_log_likelihood_profile_grad,
                vel_mean=rep(1,length(vel_x_end_list[[index_frame]])),vel=vel_x_end_list[[index_frame]],method = "L-BFGS-B",
                lower = c(.1,0),upper = c(2,1),control=list(fnscale=-1))
  n_t=length(vel_x_end_list[[index_frame]])
  beta_x_ggd=m_GGD_x$par[1]
  w_x_ggd=m_GGD_x$par[2]
  alpha_x_ggd=(sum(abs(vel_x_end_list[[index_frame]]-w_x_ggd*rep(1,length(vel_x_end_list[[index_frame]])))^beta_x_ggd)*beta_x_ggd/n_t)^(1/beta_x_ggd)
  density_ggd_x=dgnorm(density_x$x, mu = w_x_ggd, 
                       alpha = alpha_x_ggd,beta=beta_x_ggd, log = FALSE)
  lines(density_x$x,log10(density_ggd_x),col='blue')
  
  
  #y-diraction
  density_y=density(vel_y_end_list[[index_frame]])
  plot(density_y$x,log10(density_y$y),ylim=c(0,log10(300)),type='l')
  mean_y_here=mean(vel_y_end_list[[index_frame]])
  sd_y_here=sd(vel_y_end_list[[index_frame]])
  ##best fit Gaussian
  density_norm_y=dnorm(density_y$x,mean=mean_y_here,sd=sd_y_here)
  lines(density_y$x,log10(density_norm_y),lty=4,col='orange')
  ##best fit Laplace
  m_laplace_vy=lad(vel_y_end_list[[index_frame]]~1)
  density_laplace_y=dlaplace(density_y$x, location = m_laplace_vy$coefficients, 
                             scale = m_laplace_vy$scale, log = FALSE)
  lines(density_y$x,log10(density_laplace_y),col='red')
  ##best fit ggd
  m_GGD_y=optim(param,fn=GGD_log_likelihood_profile,gr=GGD_log_likelihood_profile_grad,
                vel_mean=rep(1,length(vel_y_end_list[[index_frame]])),vel=vel_y_end_list[[index_frame]],method = "L-BFGS-B",
                lower = c(.1,0),upper = c(2,1),control=list(fnscale=-1))
  n_t=length(vel_y_end_list[[index_frame]])
  beta_y_ggd=m_GGD_y$par[1]
  w_y_ggd=m_GGD_y$par[2]
  alpha_y_ggd=(sum(abs(vel_y_end_list[[index_frame]]-w_y_ggd*rep(1,length(vel_y_end_list[[index_frame]])))^beta_y_ggd)*beta_y_ggd/n_t)^(1/beta_y_ggd)
  density_ggd_y=dgnorm(density_y$x, mu = w_y_ggd, 
                       alpha = alpha_y_ggd,beta=beta_y_ggd, log = FALSE)
  lines(density_y$x,log10(density_ggd_y),col='blue')
  velocity_density_mat=cbind(density_x$x,density_x$y,density_norm_x,density_laplace_x,density_ggd_x,
                             density_y$x,density_y$y,density_norm_y,density_laplace_y,density_ggd_y)
  colnames(velocity_density_mat)=c('input_x','density_x','density_norm_x','density_laplace_x','density_ggd_x',
                                   'input_y','density_y','density_norm_y','density_laplace_y','density_ggd_y')
  
  write.csv(velocity_density_mat,file=paste0('csv_variables/marginal_velocity_density_time_',index_frame,'_experiment_',file_number,'.csv'),row.names = F)
  
}

##f test for equality of var
if(run_permutation_F_test){
  
  F_test_statistics_p_value_record=matrix(NA,T_time,2)
  B=10^4
  F_test_sample_B_record=matrix(NA,T_time,B)
  
  record_time_F_test=system.time(
    for(i_t in 1:T_time){
      print(i_t)
      res_permutation_F_test=permutation_F_test(vx=vel_x_list[[i_t]], vy=vel_y_list[[i_t]], permutation=T,B=B,record_permutation=T)
      F_test_statistics_p_value_record[i_t,]=c(res_permutation_F_test$statistics,res_permutation_F_test$p_value)
      F_test_sample_B_record[i_t,]=res_permutation_F_test$F_statistics_B
      print( F_test_statistics_p_value_record[i_t,])
    }
  )
  colnames(F_test_statistics_p_value_record)=c('F_test_statistics','p_value')
  write.csv(F_test_statistics_p_value_record,file=paste0('csv_variables/F_test_statistics_p_value_record_experiment',file_number,'.csv'),row.names = F)
  write.csv(F_test_sample_B_record,file=paste0('csv_variables/F_test_sample_B_record_experiment',file_number,'.csv'),row.names = F)
  
}

##trend of sd
plot(sqrt(sigma_2_x))
plot(sqrt(sigma_2_y))

##normality test 
##x
record_normality_test_x=matrix(NA,T_time,2)
for(i in 1:T_time){
  test=shapiro.test(vel_x_list[[i]])
  record_normality_test_x[i,1]=test$statistic
  record_normality_test_x[i,2]=test$p.value
}
colnames(record_normality_test_x)=c('W','p-value')
write.csv(record_normality_test_x,file=paste0('csv_variables/record_normality_test_x_experiment_',file_number,'.csv'),row.names = F)
##y 
record_normality_test_y=matrix(NA,T_time,2)
for(i in 1:T_time){
  test=shapiro.test(vel_y_list[[i]])
  record_normality_test_y[i,1]=test$statistic
  record_normality_test_y[i,2]=test$p.value
  
}
colnames(record_normality_test_y)=c('W','p-value')
write.csv(record_normality_test_y,file=paste0('csv_variables/record_normality_test_y_experiment_',file_number,'.csv'),row.names = F)


#####this is fast code by forming grid
px_min=min(velocity_loc_frame_record[,3])
px_max=max(velocity_loc_frame_record[,3])
py_min=min(velocity_loc_frame_record[,4])
py_max=max(velocity_loc_frame_record[,4])

# number of boxes on each direction
nx=30
ny=20
# get the information for each boxes
grid_boundary_info=get_boundary_grid(px_min=px_min,px_max=px_max,
                                     py_min=py_min,py_max=py_max,nx=nx,ny=ny)

#set a upper bound for interaction radius r
cut_r_max=100

###test neighbor for comparing correlation
if(run_test_neighbor){
  #apolar vicsek
  apolar_vicsek=T 
  max_neighbors_list_apolar_vicsek=find_max_neighbors_fast_grid(pos_x_list,pos_y_list, vel_x_list,vel_y_list,n_record,T_time,
                                                                grid_boundary_info,cut_r_max=cut_r_max,
                                                                apolar_vicsek=apolar_vicsek)
  
  v_max_neighbor_x_record=max_neighbors_list_apolar_vicsek$v_neighbor_x_record
  v_max_neighbor_y_record=max_neighbors_list_apolar_vicsek$v_neighbor_y_record
  
  d_pos_max_vec=max_neighbors_list_apolar_vicsek$d_pos_vec
  num_neighbors_max_vec=max_neighbors_list_apolar_vicsek$num_neighbors_vec
  
  cut_r=75
  residual_type=c('Laplace') #c('Gaussian','Laplace','GGD')
  param_est=get_param_est_with_max_neighbor(vel_x_list,vel_y_list,vel_x_end_list,vel_y_end_list,n_record,T_time,
                                            v_max_neighbor_x_record,v_max_neighbor_y_record,
                                            d_pos_max_vec,num_neighbors_max_vec,cut_r=cut_r,
                                            fixed_weight=F,residual_type=residual_type)
  
  #polar vicsek
  max_neighbors_list_polar_vicsek=find_max_neighbors_fast_grid(pos_x_list,pos_y_list, vel_x_list,vel_y_list,n_record,T_time,
                                                               grid_boundary_info,cut_r_max=cut_r_max,
                                                               apolar_vicsek=F)
  
  
  v_max_neighbor_x_record_polar=max_neighbors_list_polar_vicsek$v_neighbor_x_record
  v_max_neighbor_y_record_polar=max_neighbors_list_polar_vicsek$v_neighbor_y_record
  
  d_pos_max_vec_polar=max_neighbors_list_polar_vicsek$d_pos_vec
  num_neighbors_max_vec_polar=max_neighbors_list_polar_vicsek$num_neighbors_vec
  
  
  param_est_polar=get_param_est_with_max_neighbor(vel_x_list,vel_y_list,vel_x_end_list,vel_y_end_list,n_record,T_time,
                                                  v_max_neighbor_x_record_polar,v_max_neighbor_y_record_polar,
                                                  d_pos_max_vec_polar,num_neighbors_max_vec_polar,cut_r=cut_r,
                                                  fixed_weight=F,residual_type=residual_type)
  
  
  #record correlation
  correlation_record=matrix(NA,T_time,12)
  colnames(correlation_record)=c('apolar_vx_est','apolar_vx_lower_95','apolar_vx_upper_95',
                                 'apolar_vy_est','apolar_vy_lower_95','apolar_vy_upper_95',
                                 'polar_vx_est','polar_vx_lower_95','polar_vx_upper_95',
                                 'polar_vy_est','polar_vy_lower_95','polar_vy_upper_95')
  for(i_t in 1:T_time){
    #print(i_t)
    ##the correlation does not depend on w
    correlation_record[i_t,1]=cor(vel_x_end_list[[i_t]],param_est$wx_hat[i_t]*param_est$v_mean_neighbor_x_list[[i_t]])
    F_r=atanh(correlation_record[i_t,1])
    lower_F_r=F_r+qnorm(0.025)/sqrt(n_record[i_t]-3)
    upper_F_r=F_r+qnorm(0.975)/sqrt(n_record[i_t]-3)
    correlation_record[i_t,2]=tanh(lower_F_r)
    correlation_record[i_t,3]=tanh(upper_F_r)
    
    correlation_record[i_t,4]=cor(vel_y_end_list[[i_t]],param_est$wy_hat[i_t]*param_est$v_mean_neighbor_y_list[[i_t]])
    F_r=atanh(correlation_record[i_t,4])
    lower_F_r=F_r+qnorm(0.025)/sqrt(n_record[i_t]-3)
    upper_F_r=F_r+qnorm(0.975)/sqrt(n_record[i_t]-3)
    correlation_record[i_t,5]=tanh(lower_F_r)
    correlation_record[i_t,6]=tanh(upper_F_r)
    
    ###polar
    correlation_record[i_t,7]=cor(vel_x_end_list[[i_t]],param_est_polar$wx_hat[i_t]*param_est_polar$v_mean_neighbor_x_list[[i_t]])
    F_r=atanh(correlation_record[i_t,7])
    lower_F_r=F_r+qnorm(0.025)/sqrt(n_record[i_t]-3)
    upper_F_r=F_r+qnorm(0.975)/sqrt(n_record[i_t]-3)
    correlation_record[i_t,8]=tanh(lower_F_r)
    correlation_record[i_t,9]=tanh(upper_F_r)
    
    correlation_record[i_t,10]=cor(vel_y_end_list[[i_t]],param_est_polar$wy_hat[i_t]*param_est_polar$v_mean_neighbor_y_list[[i_t]])
    F_r=atanh(correlation_record[i_t,10])
    lower_F_r=F_r+qnorm(0.025)/sqrt(n_record[i_t]-3)
    upper_F_r=F_r+qnorm(0.975)/sqrt(n_record[i_t]-3)
    correlation_record[i_t,11]=tanh(lower_F_r)
    correlation_record[i_t,12]=tanh(upper_F_r)
    
  }
  
  write.csv(correlation_record,file=paste0('csv_variables/correlation_record_experiment_',file_number,'.csv'),row.names = F)
  
  index_frame=60
  data_time$apolar_mean_x=param_est$v_mean_neighbor_x_list[[index_frame]]
  data_time$apolar_mean_y=param_est$v_mean_neighbor_y_list[[index_frame]]
  data_time$polar_mean_x=param_est_polar$v_mean_neighbor_x_list[[index_frame]]
  data_time$polar_mean_y=param_est_polar$v_mean_neighbor_y_list[[index_frame]]
  write.csv(data_time,file=paste0('csv_variables/data_at_time_',index_frame,'_exp_',file_number,'.csv'),row.names = F)
  
}


## 95% interval of weights using Laplace fit
if(interval_est){
  #record weight
  weight_record=matrix(NA,T_time,12)
  colnames(weight_record)=c('apolar_wx_est','apolar_wx_lower_95','apolar_wx_upper_95',
                            'apolar_wy_est','apolar_wy_lower_95','apolar_wy_upper_95',
                            'polar_wx_est','polar_wx_lower_95','polar_wx_upper_95',
                            'polar_wy_est','polar_wy_lower_95','polar_wy_upper_95')
  
  apolar_vicsek=T 
  max_neighbors_list_apolar_vicsek=find_max_neighbors_fast_grid(pos_x_list,pos_y_list, vel_x_list,vel_y_list,n_record,T_time,
                                                                grid_boundary_info,cut_r_max=cut_r_max,
                                                                apolar_vicsek=apolar_vicsek)
  
  v_max_neighbor_x_record=max_neighbors_list_apolar_vicsek$v_neighbor_x_record
  v_max_neighbor_y_record=max_neighbors_list_apolar_vicsek$v_neighbor_y_record
  
  d_pos_max_vec=max_neighbors_list_apolar_vicsek$d_pos_vec
  num_neighbors_max_vec=max_neighbors_list_apolar_vicsek$num_neighbors_vec
  
  
  #4.get parameter est
  cut_r=75
  residual_type=c('Laplace') #c('Gaussian','Laplace','GGD')
  
  param_interval_est=get_param_est_with_max_neighbor(vel_x_list,vel_y_list,vel_x_end_list,vel_y_end_list,n_record,T_time,
                                                     v_max_neighbor_x_record,v_max_neighbor_y_record,
                                                     d_pos_max_vec,num_neighbors_max_vec,cut_r=cut_r,
                                                     fixed_weight=F,residual_type=residual_type,
                                                     interval = T,B=50)
  weight_record[,1]=param_interval_est$wx_hat
  weight_record[,2]=param_interval_est$w_x_LB95
  weight_record[,3]=param_interval_est$w_x_UB95
  weight_record[,4]=param_interval_est$wy_hat
  weight_record[,5]=param_interval_est$w_y_LB95
  weight_record[,6]=param_interval_est$w_y_UB95
  
  #polar
  max_neighbors_list_polar_vicsek=find_max_neighbors_fast_grid(pos_x_list,pos_y_list, vel_x_list,vel_y_list,n_record,T_time,
                                                               grid_boundary_info,cut_r_max=cut_r_max,
                                                               apolar_vicsek=F)
  
  
  v_max_neighbor_x_record_polar=max_neighbors_list_polar_vicsek$v_neighbor_x_record
  v_max_neighbor_y_record_polar=max_neighbors_list_polar_vicsek$v_neighbor_y_record
  
  d_pos_max_vec_polar=max_neighbors_list_polar_vicsek$d_pos_vec
  num_neighbors_max_vec_polar=max_neighbors_list_polar_vicsek$num_neighbors_vec
  
  
  param_interval_est_polar=get_param_est_with_max_neighbor(vel_x_list,vel_y_list,vel_x_end_list,vel_y_end_list,n_record,T_time,
                                                           v_max_neighbor_x_record_polar,v_max_neighbor_y_record_polar,
                                                           d_pos_max_vec_polar,num_neighbors_max_vec_polar,cut_r=cut_r,
                                                           fixed_weight=F,residual_type=residual_type,
                                                           interval = T,B=50)
  weight_record[,7]=param_interval_est_polar$wx_hat
  weight_record[,8]=param_interval_est_polar$w_x_LB95
  weight_record[,9]=param_interval_est_polar$w_x_UB95
  weight_record[,10]=param_interval_est_polar$wy_hat
  weight_record[,11]=param_interval_est_polar$w_y_LB95
  weight_record[,12]=param_interval_est_polar$w_y_UB95
  
  csv_directory=paste0('csv_variables/weight_interval_est_mat_',residual_type,'_',file_number,'.csv')
  write.csv(weight_record,file=csv_directory,row.names = F)
  
}


## distributions of random fluctuations between the observations and simulated models
if(run_fit_conditional_velocity_simulation){
  residual_type_vec=c('Gaussian','Laplace','GGD') 
  index_frame=60
  par(mfrow=c(1,2))
  
  apolar_vicsek=T 
  max_neighbors_list_apolar_vicsek=find_max_neighbors_fast_grid(pos_x_list,pos_y_list, vel_x_list,vel_y_list,n_record,T_time,
                                                                grid_boundary_info,cut_r_max=cut_r_max,
                                                                apolar_vicsek=apolar_vicsek)
  
  v_max_neighbor_x_record=max_neighbors_list_apolar_vicsek$v_neighbor_x_record
  v_max_neighbor_y_record=max_neighbors_list_apolar_vicsek$v_neighbor_y_record
  
  d_pos_max_vec=max_neighbors_list_apolar_vicsek$d_pos_vec
  num_neighbors_max_vec=max_neighbors_list_apolar_vicsek$num_neighbors_vec
  
  cut_r=75
  
  T_sim=T_time
  h=1200
  ##no duplication now
  p0=cbind(pos_x_list[[1]],pos_y_list[[1]])  ##
  v0=cbind(vel_x_list[[1]],vel_y_list[[1]])
  n_sim=dim(p0)[1]
  
  for(i_type in 1:length(residual_type_vec)){
    residual_type=residual_type_vec[i_type]
    set.seed(2)
    param_est=get_param_est_with_max_neighbor(vel_x_list,vel_y_list,vel_x_end_list,vel_y_end_list,n_record,T_time,
                                              v_max_neighbor_x_record,v_max_neighbor_y_record,
                                              d_pos_max_vec,num_neighbors_max_vec,cut_r=cut_r,
                                              fixed_weight=F,residual_type=residual_type)
    
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
    
    m_simulation=simulate_cells_anisotropic_with_neighbor_fast_grid(cut_r,p0,v0,sigma_2_0_x_hat,sigma_2_0_y_hat,
                                                                    T_sim,h,wx_hat,wy_hat,grid_boundary_info,
                                                                    apolar_vicsek=apolar_vicsek,
                                                                    alpha_x_hat=alpha_x_hat,alpha_y_hat=alpha_y_hat,beta_x_hat=beta_x_hat,beta_y_hat=beta_y_hat,
                                                                    residual_type=residual_type)
    
    residual_x=vel_x_end_list[[index_frame]]-param_est$wx_hat[index_frame]*param_est$v_mean_neighbor_x_list[[index_frame]]
    density_x=density(residual_x)
    plot(density_x$x,log10(density_x$y),ylim=c(0,log10(150)),type='l',main=residual_type)
    
    residual_x_simulated=m_simulation$vel_x_sim_record[,index_frame+1]-param_est$wx_hat[index_frame]*m_simulation$v_mean_neighbor_x_sim_record[,index_frame]
    density_x_simulated=density(residual_x_simulated)
    lines(density_x_simulated$x,log10(density_x_simulated$y),col='blue')
    
    residual_y=vel_y_end_list[[index_frame]]-param_est$wx_hat[index_frame]*param_est$v_mean_neighbor_y_list[[index_frame]]
    density_y=density(residual_y)
    plot(density_y$x,log10(density_y$y),ylim=c(0,log10(300)),type='l',main=residual_type)
    
    residual_y_simulated=m_simulation$vel_y_sim_record[,index_frame+1]-param_est$wx_hat[index_frame]*m_simulation$v_mean_neighbor_y_sim_record[,index_frame]
    density_y_simulated=density(residual_y_simulated)
    lines(density_y_simulated$x,log10(density_y_simulated$y),col='blue')
    
    velocity_conditional_density_mat=cbind(density_x$x,density_x$y,density_x_simulated$x,density_x_simulated$y,
                                           density_y$x,density_y$y,density_y_simulated$x,density_y_simulated$y)
    
    colnames(velocity_conditional_density_mat)=c('input_x','density_x','input_x_simulated','density_x_simulated',
                                                 'input_y','density_y','input_y_simulated','density_y_simulated')
    
    write.csv(velocity_conditional_density_mat,file=paste0('csv_variables/conditional_velocity_density_time_',index_frame,'_experiment_',file_number,'_residual_',residual_type,'.csv'),row.names = F)
    
  }
  
}


