### phase diagram of wx_hat/wy_hat and sigma_0_x_hat/sigma_0_y_hat
library(Rcpp)
library(RcppEigen)
library(RobustGaSP)
source('functions/function_cell.R')
library(L1pack)

file_number=441 ## please choose a number, 
delete_first_2_hours=T
unsorted_frameID=F
unsorted_particleID=F
check_convergence=T
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
       x1=velocity_loc_frame_record[index_plot,]$px+velocity_loc_frame_record[index_plot,]$vx*3600,
       y1=velocity_loc_frame_record[index_plot,]$py+velocity_loc_frame_record[index_plot,]$vy*3600,
       length=0.03
)


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

for(t in 1:T_time){
  sigma_2_x[t]=var(vx_all[count_sum+1:n_record[t]])
  sigma_2_y[t]=var(vy_all[count_sum+1:n_record[t]])
  
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

###this is original order
order_param=rep(NA,T_time)

for(t in 1:T_time){
  order_param[t]=sum(cos(2*theta_0_pi_list[[t]]))/n_record[t]
}

plot(order_param)
##trend of sd
plot(sqrt(sigma_2_x))
plot(sqrt(sigma_2_y))

plot(sqrt(sigma_2_x)/(sqrt(sigma_2_y)))


px_min=min(velocity_loc_frame_record[,3])
px_max=max(velocity_loc_frame_record[,3])
py_min=min(velocity_loc_frame_record[,4])
py_max=max(velocity_loc_frame_record[,4])


nx=30
ny=20

grid_boundary_info=get_boundary_grid(px_min=px_min,px_max=px_max,
                                     py_min=py_min,py_max=py_max,nx=nx,ny=ny)

cut_r_max=100
apolar_vicsek=T

max_neighbors_list_apolar_vicsek=find_max_neighbors_fast_grid(pos_x_list,pos_y_list, vel_x_list,vel_y_list,n_record,T_time,
                                                              grid_boundary_info,cut_r_max=cut_r_max,
                                                              apolar_vicsek=apolar_vicsek)

v_max_neighbor_x_record=max_neighbors_list_apolar_vicsek$v_neighbor_x_record
v_max_neighbor_y_record=max_neighbors_list_apolar_vicsek$v_neighbor_y_record

d_pos_max_vec=max_neighbors_list_apolar_vicsek$d_pos_vec
num_neighbors_max_vec=max_neighbors_list_apolar_vicsek$num_neighbors_vec


#4.get parameter est
##for Laplace, use L1pack for now for LAD, later we will implement IRLS
cut_r=75 ###this is fixed for now 
residual_type=c('Laplace') #c('Gaussian','Laplace','GGD')
param_est=get_param_est_with_max_neighbor(vel_x_list,vel_y_list,vel_x_end_list,vel_y_end_list,n_record,T_time,
                                          v_max_neighbor_x_record,v_max_neighbor_y_record,
                                          d_pos_max_vec,num_neighbors_max_vec,cut_r=cut_r,
                                          fixed_weight=F,residual_type=residual_type)

sigma_2_0_x_hat=param_est$sigma_2_0_x_hat
sigma_2_0_y_hat=param_est$sigma_2_0_y_hat
plot(sqrt(sigma_2_0_x_hat/sigma_2_0_y_hat))
min(sqrt(sigma_2_0_x_hat/sigma_2_0_y_hat))
max(sqrt(sigma_2_0_x_hat/sigma_2_0_y_hat))

wx_hat=param_est$wx_hat
wy_hat=param_est$wy_hat
plot(wx_hat)
plot(wy_hat)

plot(wx_hat/wy_hat)
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


###phase diagram
h=1200
##no duplication now
p0=cbind(pos_x_list[[1]],pos_y_list[[1]])  ##
v0=cbind(vel_x_list[[1]],vel_y_list[[1]])
n_sim=dim(p0)[1]
T_sim=70
T_sim_burn_in=20
n1=50
n2=50
order_diagram_record=matrix(NA,n1,n2)
apolar_vicsek=T
#fix sigma_2_0_x, vary the ratio between sigma_2_0_x and sigma_2_0_y
sigma_2_0_x_base=2*10^{-5}
ratio_x_y_sigma_0=seq(1.2,2,0.8/(n1-1)) ##ratio for sigma_0_x and sigma_0_y
#fix wx, vary the ratio between wx and wy
wx_base=0.7
ratio_x_y_w=seq(0.9,1.5,0.6/(n2-1)) ##ratio for sigma_0_x and sigma_0_y

if(check_convergence){
  set.seed(0)
  sample_size_check=10
  index_sample=sample(1:length(order_diagram_record),sample_size_check)
  
  i_sample=ceiling(index_sample/n2)
  j_sample=index_sample-(i_sample-1)*n2
  
  param_sample_record=matrix(NA,sample_size_check,4)
  colnames(param_sample_record)=c('wx_sim','wy_sim','sigma_0_x_sim','sigma_0_y_sim')
  order_sample_record=matrix(NA,sample_size_check,T_sim)
  
  for(k in 1:sample_size_check){
    print(k)
    set.seed(k)
    wx_sim=rep(wx_base,T_sim)
    wy_sim=rep(wx_base/ratio_x_y_w[i_sample[k]],T_sim) 
    
    sigma_2_0_x_sim=rep(sigma_2_0_x_base,T_sim)  ##tau in the paper 
    sigma_2_0_y_sim=rep(sigma_2_0_x_base/ratio_x_y_sigma_0[j_sample[k]]^2,T_sim) 
    param_sample_record[k,]=c(wx_sim[1],wy_sim[1],sqrt(sigma_2_0_x_sim[1]),sqrt(sigma_2_0_y_sim[1]))
    m_simulation=simulate_cells_anisotropic_with_neighbor_fast_grid(cut_r,p0,v0,sigma_2_0_x_sim,sigma_2_0_y_sim,
                                                                    T_sim,h,wx_sim,wy_sim,grid_boundary_info,
                                                                    apolar_vicsek=apolar_vicsek,
                                                                    alpha_x_hat=alpha_x_sim,alpha_y_hat=alpha_y_sim,beta_x_hat=beta_x_sim,beta_y_hat=beta_y_sim,
                                                                    residual_type=residual_type)
    order_sample_record[k,]=(m_simulation$record_ensemble_order)
    
  }
  
  matplot(t(order_sample_record),type='l')
}


if(residual_type=='Laplace'){
  for(i_sim in 1:n1){
    print(i_sim)
    for(j_sim in 1:n2){
      print(j_sim)
      
      wx_sim=rep(wx_base,T_sim)
      wy_sim=rep(wx_base/ratio_x_y_w[i_sim],T_sim) 
      
      sigma_2_0_x_sim=rep(sigma_2_0_x_base,T_sim)  ##tau in the paper 
      sigma_2_0_y_sim=rep(sigma_2_0_x_base/ratio_x_y_sigma_0[j_sim]^2,T_sim) 
      
      alpha_x_sim=alpha_y_sim=beta_x_sim=beta_y_sim=NULL
      
      ##p0 and v0 is from original, but they may not affect too much 
      
      m_simulation=simulate_cells_anisotropic_with_neighbor_fast_grid(cut_r,p0,v0,sigma_2_0_x_sim,sigma_2_0_y_sim,
                                                                      T_sim,h,wx_sim,wy_sim,grid_boundary_info,
                                                                      apolar_vicsek=apolar_vicsek,
                                                                      alpha_x_hat=alpha_x_sim,alpha_y_hat=alpha_y_sim,beta_x_hat=beta_x_sim,beta_y_hat=beta_y_sim,
                                                                      residual_type=residual_type)
      order_diagram_record[i_sim,j_sim]=mean(m_simulation$record_ensemble_order[(T_sim_burn_in+1):T_sim])
      print(order_diagram_record[i_sim,j_sim])
    }
  }
  
  w_ratio_data=wx_hat/wy_hat
  sigma_0_ratio_data=sqrt(sigma_2_0_x_hat/sigma_2_0_y_hat)
  library(plot3D)
  image2D(t(order_diagram_record),x=ratio_x_y_sigma_0,y=ratio_x_y_w,xlab='sigma_0_x/sigma_0_y',ylab='w_x/w_y')
  lines(sigma_0_ratio_data,w_ratio_data)
  w_sigma_0_order_diagram_data=cbind(sigma_0_ratio_data,w_ratio_data,order_param)
  w_sigma_0_diagram=cbind(ratio_x_y_sigma_0,ratio_x_y_w)
  
  colnames(w_sigma_0_diagram)=c('sigma_0_x/sigma_0_y','w_x/w_y') ##sigma_0 is tau in the paper
  write.csv(order_diagram_record,file=paste0('csv_variables/phase_diagram_sim_experiment_',file_number,'.csv'),row.names = F) 
  write.csv(w_sigma_0_diagram,file=paste0('csv_variables/w_sigma_0_diagram_sim_experiment_',file_number,'.csv'),row.names = F) 
  write.csv(w_sigma_0_order_diagram_data,file=paste0('csv_variables/w_sigma_0_order_diagram_data_experiment_',file_number,'.csv'),row.names = F) 
}


