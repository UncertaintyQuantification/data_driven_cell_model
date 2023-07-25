#obtain position, velocity, velocity angle (theta) from the dataset
#only consider the particles that appear in two consecutive time frames
get_pos_v_theta_all_vec<-function(velocity_loc_frame_record,T_time){
  
  start_frame=0
  N_tilde=0
  pos_all_vec= NULL ##position at previous frame
  v_all=NULL  ##position at current frame
  v_all_end=NULL  ##position at current frame
  theta_all=NULL ##velocity anagle
  n_record=rep(NA,T_time) ##record particles at each frame (they need to exist between two frames)
  count_frame=0
  
  for(i_frame in start_frame:(start_frame+T_time-1)){
    count_frame=count_frame+1
    index_start=which(velocity_loc_frame_record$frameID==i_frame)
    index_end=which(velocity_loc_frame_record$frameID==(i_frame+1))
    
    shared_particleID=intersect( velocity_loc_frame_record[index_start,]$particleID,
                                 velocity_loc_frame_record[index_end,]$particleID)
    n_record[count_frame]=length(shared_particleID)
    
    index_select_frame_start=NULL
    count=1
    for(i in 1:length(index_start)){ ##
      if(count<=length(shared_particleID)){
        if(velocity_loc_frame_record[index_start[i],]$particleID==shared_particleID[count]){
          index_select_frame_start=c(index_select_frame_start,i)
          count=count+1
        }
      }
      
    }
    
    index_select_frame_end=NULL
    #index_select_end=NULL
    count=1
    for(i in 1:length(index_end)){ ##
      if(count<=length(shared_particleID)){
        if(velocity_loc_frame_record[index_end[i],]$particleID==shared_particleID[count]){
          index_select_frame_end=c(index_select_frame_end,i)
          count=count+1
        }
      }
    }
    
    
    
    N_tilde=N_tilde+length(index_select_frame_end)*2 ###velocity obs
    ##velocity at the start frame 
    v_all=c(v_all,as.vector(t(velocity_loc_frame_record[index_start[index_select_frame_start],1:2]))) ##actually this velocity is the one that is between these two
    ##velocity at the end frame 
    v_all_end=c(v_all_end,as.vector(t(velocity_loc_frame_record[index_end[index_select_frame_end],1:2]))) ##actually this velocity is the one that is between these two
    
    ##velocity angle at the start frame
    theta_here=atan2(velocity_loc_frame_record[index_start[index_select_frame_start],2],
                     velocity_loc_frame_record[index_start[index_select_frame_start],1])
    
    theta_all=c(theta_all,theta_here)
    
    
    ##position at the start frame 
    pos_all_vec=c(pos_all_vec,as.vector(t(velocity_loc_frame_record[index_start[index_select_frame_start],3:4])))
  }
  
  ans_list=list()
  ans_list$v_all=v_all
  ans_list$v_all_end=v_all_end
  ans_list$theta_all=theta_all
  ans_list$pos_all_vec=pos_all_vec
  ans_list$n_record=n_record
  
  
  return(ans_list)
  
}


#return the indices of duplicated positions
delete_duplicate_pos<-function(pos,v){
  ##pos n times 2
  ##v n times 2
  n_ori=dim(pos)[1]
  
  ##I need they to be large and outside unique shouldn't be a unique number
  pos_x_unique=rep(-9^9,n_ori)
  pos_y_unique=rep(-9^9,n_ori) 
  pos_x_unique[1:length(unique(pos[,1]))]=unique(pos[,1])
  pos_y_unique[1:length(unique(pos[,2]))]=unique(pos[,2])
  delete_index_x=NULL
  delete_index_y=NULL
  
  count_x=1
  count_y=1
  
  for(i in 1:n_ori){
    if(pos[i,1]==pos_x_unique[count_x]){
      count_x=count_x+1
    }else{
      delete_index_x=c(delete_index_x,i)
    }
    if(pos[i,2]==pos_y_unique[count_y]){
      count_y=count_y+1
    }else{
      delete_index_y=c(delete_index_y,i)
    }
    
  }
  

  delete_index=intersect(delete_index_x,delete_index_y)
  return(delete_index)
}

#return the information of each grid/box
get_boundary_grid<-function(px_min,px_max,py_min,py_max,nx=30,ny=20){
  grid_boundary_mat=matrix(NA,nx*ny,4)
  colnames(grid_boundary_mat)=c('pos_x_min','pos_x_max','pos_y_min','pos_y_max')
  
  ##get slightly larger grid 
  len_x_ori=px_max-px_min
  len_y_ori=py_max-py_min
  delta_x=len_x_ori/nx*0.1
  delta_y=len_y_ori/ny*0.1
  
  px_seq=seq(px_min-delta_x,px_max+delta_x,(px_max-px_min+2*delta_x)/(nx))
  py_seq=seq(py_min-delta_y,py_max+delta_y,(py_max-py_min+2*delta_y)/(ny))
  grid_boundary_mat[,1]=rep(px_seq[1:nx],ny)
  grid_boundary_mat[,2]=rep(px_seq[2:(nx+1)],ny)
  grid_boundary_mat[,3]=as.numeric(t(matrix(py_seq[1:ny],ny,nx)))
  grid_boundary_mat[,4]=as.numeric(t(matrix(py_seq[2:(ny+1)],ny,nx)))

  my_grid=list()
  
  
  Lx_min = min(grid_boundary_mat[,1:2]);
  Lx_max = max(grid_boundary_mat[,1:2]);
  Ly_min = min(grid_boundary_mat[,3:4]);
  Ly_max = max(grid_boundary_mat[,3:4]);
  
  len_x=  (max(grid_boundary_mat[,1:2])-min(grid_boundary_mat[,1:2]))/nx
  len_y=  (max(grid_boundary_mat[,3:4])-min(grid_boundary_mat[,3:4]))/ny
  
  
  grid_boundary_info=list()
  grid_boundary_info$grid_boundary_mat=grid_boundary_mat
  grid_boundary_info$grid_info=as.matrix(c(Lx_min,Lx_max,Ly_min,Ly_max,nx,ny,len_x,len_y),8,1)
  rownames(  grid_boundary_info$grid_info)=c('Lx_min','Lx_max','Ly_min','Ly_max','nx','ny','len_x','len_y')
  
  return(grid_boundary_info)
  
}


initiate_grid<-function(grid_boundary_info,pos_x,pos_y,vel_x,vel_y,neighbor_index_list=NULL){
  
  #m_grid stores the velocities and positions of the cells 
  #in i-th grid and in the neighboring 9 grids of i-th grid
  m_grid=as.list(rep(NA,nx*ny))
  n_t=length(pos_x)
  
  for(i in 1:(nx*ny)){
    m_grid[[i]]=list()
    m_grid[[i]]$particle_pos=NA
    m_grid[[i]]$neighbor_pos=NA
    
    m_grid[[i]]$particle_vel=NA
    m_grid[[i]]$neighbor_vel=NA
  }

  
  Lx_min=grid_boundary_info$grid_info[1]
  Lx_max=grid_boundary_info$grid_info[2]
  Ly_min=grid_boundary_info$grid_info[3]
  Ly_max=grid_boundary_info$grid_info[4]
  nx=grid_boundary_info$grid_info[5]
  ny=grid_boundary_info$grid_info[6]
  len_x=grid_boundary_info$grid_info[7]
  len_y=grid_boundary_info$grid_info[8]

  for(i in 1:n_t){
    i_x=ceiling((pos_x[i]-Lx_min)/len_x)
    i_y=ceiling((pos_y[i]-Ly_min)/len_y)
    
    index_grid=(i_y-1)*nx+i_x
    
    ##update position and velocity for each cell in the corresponding grid
    if(is.na(m_grid[[index_grid]]$particle_pos)[1]){
      m_grid[[index_grid]]$particle_pos=(as.vector(c(pos_x[i],pos_y[i])))
      m_grid[[index_grid]]$particle_vel=(as.vector(c(vel_x[i],vel_y[i])))
    }else{
      m_grid[[index_grid]]$particle_pos=cbind((m_grid[[index_grid]]$particle_pos),c(pos_x[i],pos_y[i]))
      m_grid[[index_grid]]$particle_vel=cbind((m_grid[[index_grid]]$particle_vel),c(vel_x[i],vel_y[i]))
    }

  }
  
  ##neighbor_index_list records the indices of neighboring 9 grids
  ##from x first
  if(is.null(neighbor_index_list)){
    neighbor_index_list=as.list(1:(nx*ny))
    ##exterior
    for(i in 1:(nx*ny)){
      i_x=(i%%nx)
      if(i_x==0){
        i_x=nx
      }
      i_y=ceiling(i/nx)
      
      #contain itself
      #neighbor_index_list[[i]]=i
      
      if((i_x-1)>0&(i_y-1)>0){
        neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-2)*nx+(i_x-1))
      }
      
      if((i_y-1)>0){
        neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-2)*nx+(i_x))
      }
      
      if((i_x+1)<=nx&(i_y-1)>0){
        neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-2)*nx+(i_x+1))
      }
      
      if((i_x-1)>0){
        neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-1)*nx+(i_x-1))
      }
      
      if((i_x+1)<=nx){
        neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-1)*nx+(i_x+1))
      }
      
      if((i_x-1)>0&(i_y+1)<=ny){
        neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y)*nx+(i_x-1))
      }
      
      if((i_y+1)<=ny){
        neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y)*nx+(i_x))
      }
      if((i_x+1)<=nx&(i_y+1)<=ny){
        neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y)*nx+(i_x+1))
      }
    }
  }
  #update positions and velocities for neighboring cells in the corresponding grid
  for(i in 1:(nx*ny)){
    #print(i)
    num_neighbor=length( neighbor_index_list[[i]])
    for(j in 1:num_neighbor){
      if(!is.na(m_grid[[neighbor_index_list[[i]][j] ]]$particle_pos)[1]){
        if(is.na(m_grid[[i]]$neighbor_pos)[1]){
          m_grid[[i]]$neighbor_pos=as.matrix(m_grid[[neighbor_index_list[[i]][j] ]]$particle_pos)
          m_grid[[i]]$neighbor_vel=as.matrix(m_grid[[neighbor_index_list[[i]][j] ]]$particle_vel)
        }else{
          m_grid[[i]]$neighbor_pos=cbind(m_grid[[i]]$neighbor_pos,m_grid[[neighbor_index_list[[i]][j] ]]$particle_pos)
          m_grid[[i]]$neighbor_vel=cbind(m_grid[[i]]$neighbor_vel,m_grid[[neighbor_index_list[[i]][j] ]]$particle_vel)
          
        }
      }
    }
  }
  ans_list=list()
  ans_list$m_grid=m_grid
  ans_list$neighbor_index_list=neighbor_index_list

  return(ans_list)

}


#define max neighbor within a radius for vicsek model 
#fast 
##change pos_all_vec and v_all to list will make life easier
find_max_neighbors_fast_grid<-function(pos_x_list,pos_y_list, vel_x_list,vel_y_list,n_record,T_time,
                                        grid_boundary_info,cut_r_max=100,
                                       apolar_vicsek=F){

   Lx_min=grid_boundary_info$grid_info[1]
   Lx_max=grid_boundary_info$grid_info[2]
   Ly_min=grid_boundary_info$grid_info[3]
   Ly_max=grid_boundary_info$grid_info[4]
   nx=grid_boundary_info$grid_info[5]
   ny=grid_boundary_info$grid_info[6]
   len_x=grid_boundary_info$grid_info[7]
   len_y=grid_boundary_info$grid_info[8]
   
   count=0
  
   num_neighbors_vec=rep(NA,sum(n_record))
   v_neighbor_x_record=rep(NA,sum(n_record)*15)  ###this is to create something large
   v_neighbor_y_record=rep(NA,sum(n_record)*15)
   d_pos_vec=rep(NA,sum(n_record)*15)

   neighbor_index_list=NULL
  for(t in 1:T_time){
    if(t>1){
      index_start=sum(n_record[1:(t-1)])
    }else{
      index_start=0
    }
    
    initiate_grid_list=initiate_grid(grid_boundary_info=grid_boundary_info,pos_x=pos_x_list[[t]],pos_y=pos_y_list[[t]],
                                     vel_x=vel_x_list[[t]],vel_y=vel_y_list[[t]],neighbor_index_list=neighbor_index_list)
    m_grid_here=initiate_grid_list$m_grid
    neighbor_index_list=initiate_grid_list$neighbor_index_list
    
    pos_x_t=pos_x_list[[t]]
    pos_y_t=pos_y_list[[t]]
    vel_x_t=vel_x_list[[t]]
    vel_y_t=vel_y_list[[t]]
    
    for(i in 1:n_record[t]){
      input_pos_i=as.vector(c(pos_x_t[i],pos_y_t[i]))
      input_vel_i=as.vector(c(vel_x_t[i],vel_y_t[i]))

      i_x=ceiling((input_pos_i[1]-Lx_min)/len_x)
      i_y=ceiling((input_pos_i[2]-Ly_min)/len_y)

      index_grid=(i_y-1)*nx+i_x
      d_vec_here_all=input_pos_i-as.matrix(m_grid_here[[index_grid]]$neighbor_pos)

      d_here=sqrt(colSums(d_vec_here_all^2))
      if(apolar_vicsek==F){
        index_neighbor=which(d_here<cut_r_max)
      }else{
        index_neighbor=which(d_here<cut_r_max)
        index_same_v_direction=which(colSums(m_grid_here[[index_grid]]$neighbor_vel*input_vel_i)>=0)
        index_neighbor=intersect(index_neighbor,index_same_v_direction)
      }

      num_neighbors_vec[index_start+i]=length(index_neighbor)

      v_neighbor_x_record[count+(1:length(index_neighbor))]= m_grid_here[[index_grid]]$neighbor_vel[1,index_neighbor]
      v_neighbor_y_record[count+(1:length(index_neighbor))]=m_grid_here[[index_grid]]$neighbor_vel[2,index_neighbor]
      d_pos_vec[count+(1:length(index_neighbor))]=d_here[index_neighbor]

      count=count+length(index_neighbor)

    }

  }


  v_neighbor_x_record=v_neighbor_x_record[1:count]
  v_neighbor_y_record=v_neighbor_y_record[1:count]
  d_pos_vec=d_pos_vec[1:count]

  ans_list=list()
  ans_list$v_neighbor_x_record=v_neighbor_x_record
  ans_list$v_neighbor_y_record=v_neighbor_y_record
  ans_list$d_pos_vec=d_pos_vec
  ans_list$num_neighbors_vec=num_neighbors_vec
  ans_list$apolar_vicsek=apolar_vicsek
  
  return(ans_list)
}



##GGD
#log-likelihood
GGD_log_likelihood_profile=function(param,vel_mean,vel){
  beta_t=param[1]
  w_t=param[2]
  n_t=length(vel)
  log_L=n_t*log(beta_t/(2*gamma(1/beta_t)))-
    n_t/beta_t*log(beta_t/n_t*sum(abs(vel-w_t*vel_mean)^beta_t))-n_t/beta_t
  return(log_L)
}
#gradient
GGD_log_likelihood_profile_grad=function(param,vel_mean,vel){
  beta_t=param[1]
  w_t=param[2]
  n_t=length(vel)
  
  e_abs=abs(vel-w_t*vel_mean)
  
  zero_ind=which(e_abs==0)
  if(length(zero_ind)>0){
    vel_mean=vel_mean[-zero_ind]
    vel=vel[-zero_ind]
    n_t=n_t-length(zero_ind)
    e_abs=abs(vel-w_t*vel_mean)
  }
  
  e_beta_sum=sum(e_abs^beta_t)
  e_beta_log_sum=sum(e_abs^beta_t*log(e_abs))
  
  d_l_d_beta_t=n_t/beta_t+n_t/(beta_t^2)*digamma(1/beta_t)+
    n_t/(beta_t^2)*log(beta_t/n_t*e_beta_sum)-n_t/beta_t*e_beta_log_sum/e_beta_sum
  
  d_l_d_w_t=n_t*sum(e_abs^(beta_t-1)*sign(vel-w_t*vel_mean)*vel_mean)/e_beta_sum
  
  return(c(d_l_d_beta_t,d_l_d_w_t))
}


###
GGD_log_likelihood_profile_fixed_weight=function(param,vel_mean,vel){
  beta_t=param

  n_t=length(vel)
  log_L=n_t*log(beta_t/(2*gamma(1/beta_t)))-
    n_t/beta_t*log(beta_t/n_t*sum(abs(vel-1*vel_mean)^beta_t))-n_t/beta_t
  return(log_L)
}

###
get_param_est_with_max_neighbor<-function(vel_x_list,vel_y_list,vel_x_end_list,vel_y_end_list,n_record,T_time,
                                          v_max_neighbor_x_record,v_max_neighbor_y_record,
                                          d_pos_max_vec,num_neighbors_max_vec,cut_r,
                                          fixed_weight=F,residual_type=residual_type,
                                          interval=F,B=50){
  
  ans=list()
  ans$wx_hat=ans$wy_hat=rep(1,T_time)
  ans$sigma_2_0_x_hat=ans$sigma_2_0_y_hat=rep(NA,T_time)
  v_mean_neighbor_x_list=as.list(1:T_time)###
  v_mean_neighbor_y_list=as.list(1:T_time)######
  
  count=0
  
  if(residual_type=="GGD"){
    ans$alpha_x_hat=rep(NA,T_time)
    ans$alpha_y_hat=rep(NA,T_time)
    ans$beta_x_hat=rep(NA,T_time)
    ans$beta_y_hat=rep(NA,T_time)
    ans$x_iteration=rep(NA,T_time)
    ans$y_iteration=rep(NA,T_time)
  }
  if(interval){
    if(residual_type=="GGD"){
      ans$alpha_x_LB95=rep(NA,T_time)
      ans$alpha_x_UB95=rep(NA,T_time)
      ans$beta_x_LB95=rep(NA,T_time)
      ans$beta_x_UB95=rep(NA,T_time)
      ans$w_x_LB95=rep(NA,T_time)
      ans$w_x_UB95=rep(NA,T_time)

      ans$alpha_y_LB95=rep(NA,T_time)
      ans$alpha_y_UB95=rep(NA,T_time)
      ans$beta_y_LB95=rep(NA,T_time)
      ans$beta_y_UB95=rep(NA,T_time)
      ans$w_y_LB95=rep(NA,T_time)
      ans$w_y_UB95=rep(NA,T_time)
    }else if(residual_type=="Laplace"){
      ans$w_x_LB95=rep(NA,T_time)
      ans$w_x_UB95=rep(NA,T_time)
      ans$phi_x_LB95=rep(NA,T_time)
      ans$phi_x_UB95=rep(NA,T_time)

      ans$w_y_LB95=rep(NA,T_time)
      ans$w_y_UB95=rep(NA,T_time)
      ans$phi_y_LB95=rep(NA,T_time)
      ans$phi_y_UB95=rep(NA,T_time)
    }else if(residual_type=="Gaussian"){
      ans$w_x_LB95=rep(NA,T_time)
      ans$w_x_UB95=rep(NA,T_time)
      ans$sigma_2_x_LB95=rep(NA,T_time)
      ans$sigma_2_x_UB95=rep(NA,T_time)

      ans$w_y_LB95=rep(NA,T_time)
      ans$w_y_UB95=rep(NA,T_time)
      ans$sigma_2_y_LB95=rep(NA,T_time)
      ans$sigma_2_y_UB95=rep(NA,T_time)
    }
    
  }
  
  for(t in 1:T_time){
    if(t>1){
      index_start=sum(n_record[1:(t-1)])
    }else{
      index_start=0
    }
    
    v_x_neighbor=rep(NA,n_record[t])
    v_y_neighbor=rep(NA,n_record[t])
    
    #find the neighboring particle with cut_r
    for(i in 1:n_record[t]){
      
      #already selected neighbor
      index_selected=which(d_pos_max_vec[count+(1:num_neighbors_max_vec[index_start+i])]<cut_r)
      
      
      v_mean_x=mean(v_max_neighbor_x_record[count+index_selected])
      v_mean_y=mean(v_max_neighbor_y_record[count+index_selected])
      
      v_x_neighbor[i]=v_mean_x
      v_y_neighbor[i]=v_mean_y
      
      count=count+(num_neighbors_max_vec[index_start+i])
      
    }
    
    #estimate parameters
    if(residual_type=='Gaussian'|residual_type=='Nonparametric'){
      if(fixed_weight==F){
        
        ans$wx_hat[t]=sum(v_x_neighbor*vel_x_end_list[[t]])/sum(v_x_neighbor*v_x_neighbor)
        ans$wy_hat[t]=sum(v_y_neighbor*vel_y_end_list[[t]])/sum(v_y_neighbor*v_y_neighbor)
        
        ans$sigma_2_0_x_hat[t]=sum((vel_x_end_list[[t]]-ans$wx_hat[t]*v_x_neighbor)^2)/(n_record[t]-1)
        ans$sigma_2_0_y_hat[t]=sum((vel_y_end_list[[t]]-ans$wy_hat[t]*v_y_neighbor)^2)/(n_record[t]-1)
        
        if(interval & residual_type=='Gaussian'){
          #Gaussian: closed form of 95% confidence interval for estimated parameters
          sd_x=sqrt(ans$sigma_2_0_x_hat[t]/sum(v_x_neighbor^2))
          ans$w_x_LB95[t]=ans$wx_hat[t]+qt(0.025,n_record[t]-1)*sd_x
          ans$w_x_UB95[t]=ans$wx_hat[t]+qt(0.975,n_record[t]-1)*sd_x
          ans$sigma_2_x_LB95[t]=(n_record[t]-1)*ans$sigma_2_0_x_hat[t]/qchisq(0.975,n_record[t]-1)
          ans$sigma_2_x_UB95[t]=(n_record[t]-1)*ans$sigma_2_0_x_hat[t]/qchisq(0.025,n_record[t]-1)
          
          
          sd_y=sqrt(ans$sigma_2_0_y_hat[t]/sum(v_y_neighbor^2))
          ans$w_y_LB95[t]=ans$wy_hat[t]+qt(0.025,n_record[t]-1)*sd_y
          ans$w_y_UB95[t]=ans$wy_hat[t]+qt(0.975,n_record[t]-1)*sd_y
          ans$sigma_2_y_LB95[t]=(n_record[t]-1)*ans$sigma_2_0_y_hat[t]/qchisq(0.975,n_record[t]-1)
          ans$sigma_2_y_UB95[t]=(n_record[t]-1)*ans$sigma_2_0_y_hat[t]/qchisq(0.025,n_record[t]-1)
        }
      }else{
        ans$sigma_2_0_x_hat[t]=sum((vel_x_end_list[[t]]- v_x_neighbor)^2)/(n_record[t]-1)
        ans$sigma_2_0_y_hat[t]=sum((vel_y_end_list[[t]]-v_y_neighbor)^2)/(n_record[t]-1)
        
      }
    }else if(residual_type=='Laplace'){
      if(fixed_weight==F){
        m_laplace_x=lad(vel_x_end_list[[t]]~v_x_neighbor-1)
        m_laplace_y=lad(vel_y_end_list[[t]]~v_y_neighbor-1)
        
        ans$wx_hat[t]=m_laplace_x$coefficients
        ans$wy_hat[t]=m_laplace_y$coefficients
        ans$sigma_2_0_x_hat[t]=m_laplace_x$scale^2 ##this is the variance equal to below 
        ans$sigma_2_0_y_hat[t]=m_laplace_y$scale^2 ##this is the variance equal to below 

        if(interval){
          #Laplace: bootstrapped 95% confidence interval for estimated parameters
          #x direction
          res_x=vel_x_end_list[[t]]-ans$wx_hat[t]*v_x_neighbor
          
          w_x_rec=rep(NA,B)
          phi_x_rec=rep(NA,B)
          
          for(b in 1:B){
            vel_new=ans$wx_hat[t]*v_x_neighbor+sample(res_x,replace = T)
            m_x_new=lad(vel_new~v_x_neighbor-1)
            
            w_x_rec[b]=m_x_new$coefficients
            phi_x_rec[b]=m_x_new$scale
          }
          ans$w_x_LB95[t]=quantile(w_x_rec,probs = 0.025)
          ans$w_x_UB95[t]=quantile(w_x_rec,probs = 0.975)
          ans$phi_x_LB95[t]=quantile(phi_x_rec,probs = 0.025)
          ans$phi_x_UB95[t]=quantile(phi_x_rec,probs = 0.975)

          #y direction
          res_y=vel_y_end_list[[t]]-ans$wy_hat[t]*v_y_neighbor
          
          w_y_rec=rep(NA,B)
          phi_y_rec=rep(NA,B)
          
          for(b in 1:B){
            vel_new=ans$wy_hat[t]*v_y_neighbor+sample(res_y,replace = T)
            m_y_new=lad(vel_new~v_y_neighbor-1)
            
            w_y_rec[b]=m_y_new$coefficients
            phi_y_rec[b]=m_y_new$scale
          }
          ans$w_y_LB95[t]=quantile(w_y_rec,probs = 0.025)
          ans$w_y_UB95[t]=quantile(w_y_rec,probs = 0.975)
          ans$phi_y_LB95[t]=quantile(phi_y_rec,probs = 0.025)
          ans$phi_y_UB95[t]=quantile(phi_y_rec,probs = 0.975)
        }
        
      }else{
        ans$sigma_2_0_x_hat[t]=(sqrt(2)*sum(abs(vel_x_end_list[[t]]-v_x_neighbor))/n_record[t])^2
        ans$sigma_2_0_y_hat[t]=(sqrt(2)*sum(abs(vel_y_end_list[[t]]-v_y_neighbor))/n_record[t])^2
        
      }
      
    }else if (residual_type=='GGD'){
      if(fixed_weight==F){
        param=c(1,.5)
        m_GGD_x=optim(param,fn=GGD_log_likelihood_profile,gr=GGD_log_likelihood_profile_grad,
                      vel_mean=v_x_neighbor,vel=vel_x_end_list[[t]],method = "L-BFGS-B",
                      lower = c(.1,0),upper = c(2,1),control=list(fnscale=-1))
        m_GGD_y=optim(param,fn=GGD_log_likelihood_profile,gr=GGD_log_likelihood_profile_grad,
                      vel_mean=v_y_neighbor,vel=vel_y_end_list[[t]],method = "L-BFGS-B",
                      lower = c(.1,0),upper = c(2,1),control=list(fnscale=-1))
        
        n_t=length(v_x_neighbor)
        
        beta_x=m_GGD_x$par[1]
        w_x=m_GGD_x$par[2]
        alpha_x=(sum(abs(vel_x_end_list[[t]]-w_x*v_x_neighbor)^beta_x)*beta_x/n_t)^(1/beta_x)
        ans$wx_hat[t]=w_x
        ans$sigma_2_0_x_hat[t]=alpha_x^2*gamma(3/beta_x)/gamma(1/beta_x)
        ans$alpha_x_hat[t]=alpha_x
        ans$beta_x_hat[t]=beta_x
        ans$x_iteration[t]=m_GGD_x$counts[1]
        
        beta_y=m_GGD_y$par[1]
        w_y=m_GGD_y$par[2]
        alpha_y=(sum(abs(vel_y_end_list[[t]]-w_y*v_y_neighbor)^beta_y)*beta_y/n_t)^(1/beta_y)
        ans$wy_hat[t]=w_y
        ans$sigma_2_0_y_hat[t]=alpha_y^2*gamma(3/beta_y)/gamma(1/beta_y)
        ans$alpha_y_hat[t]=alpha_y
        ans$beta_y_hat[t]=beta_y
        ans$y_iteration[t]=m_GGD_y$counts[1]
        
        if(interval){
          #GGD: bootstrapped 95% confidence interval for estimated parameters
          #x direction
          res_x=vel_x_end_list[[t]]-w_x*v_x_neighbor
          
          alpha_x_rec=rep(NA,B)
          beta_x_rec=rep(NA,B)
          w_x_rec=rep(NA,B)
          
          for(b in 1:B){
            vel_new=w_x*v_x_neighbor+sample(res_x,replace = T)
            m_x_new=optim(param,fn=GGD_log_likelihood_profile,gr=GGD_log_likelihood_profile_grad,
                          vel_mean=v_x_neighbor,vel=vel_new,method = "L-BFGS-B",
                          lower = c(.1,0),upper = c(2,1),control=list(fnscale=-1))
            
            beta_x_new=m_x_new$par[1]
            w_x_new=m_x_new$par[2]
            alpha_x_new=(sum(abs(vel_new-w_x_new*v_x_neighbor)^beta_x_new)*beta_x_new/n_t)^(1/beta_x_new)
            
            alpha_x_rec[b]=alpha_x_new
            beta_x_rec[b]=beta_x_new
            w_x_rec[b]=w_x_new
          }
          ans$alpha_x_LB95[t]=quantile(alpha_x_rec,probs = 0.025)
          ans$alpha_x_UB95[t]=quantile(alpha_x_rec,probs = 0.975)
          ans$beta_x_LB95[t]=quantile(beta_x_rec,probs = 0.025)
          ans$beta_x_UB95[t]=quantile(beta_x_rec,probs = 0.975)
          ans$w_x_LB95[t]=quantile(w_x_rec,probs = 0.025)
          ans$w_x_UB95[t]=quantile(w_x_rec,probs = 0.975)

          #y direction
          res_y=vel_y_end_list[[t]]-w_y*v_y_neighbor
          
          alpha_y_rec=rep(NA,B)
          beta_y_rec=rep(NA,B)
          w_y_rec=rep(NA,B)
          
          for(b in 1:B){
            vel_new=w_y*v_y_neighbor+sample(res_y,replace = T)
            m_y_new=optim(param,fn=GGD_log_likelihood_profile,gr=GGD_log_likelihood_profile_grad,
                          vel_mean=v_y_neighbor,vel=vel_new,method = "L-BFGS-B",
                          lower = c(.1,0),upper = c(2,1),control=list(fnscale=-1))
            
            beta_y_new=m_y_new$par[1]
            w_y_new=m_y_new$par[2]
            alpha_y_new=(sum(abs(vel_new-w_y_new*v_y_neighbor)^beta_y_new)*beta_y_new/n_t)^(1/beta_y_new)
            
            alpha_y_rec[b]=alpha_y_new
            beta_y_rec[b]=beta_y_new
            w_y_rec[b]=w_y_new
          }
          ans$alpha_y_LB95[t]=quantile(alpha_y_rec,probs = 0.025)
          ans$alpha_y_UB95[t]=quantile(alpha_y_rec,probs = 0.975)
          ans$beta_y_LB95[t]=quantile(beta_y_rec,probs = 0.025)
          ans$beta_y_UB95[t]=quantile(beta_y_rec,probs = 0.975)
          ans$w_y_LB95[t]=quantile(w_y_rec,probs = 0.025)
          ans$w_y_UB95[t]=quantile(w_y_rec,probs = 0.975)
        }
      }else{##fixed weight
        m_GGD_x=optimize(f=GGD_log_likelihood_profile_fixed_weight,vel_mean=v_x_neighbor,vel=vel_x_end_list[[t]],lower=0.1,upper=2,maximum=T)
        m_GGD_y=optimize(f=GGD_log_likelihood_profile_fixed_weight,vel_mean=v_y_neighbor,vel=vel_y_end_list[[t]],lower=0.1,upper=2,maximum=T)
        
        n_t=length(v_x_neighbor)
        
        beta_x=m_GGD_x$maximum
        w_x=1
        alpha_x=(sum(abs(vel_x_end_list[[t]]-w_x*v_x_neighbor)^beta_x)*beta_x/n_t)^(1/beta_x)
        ans$sigma_2_0_x_hat[t]=alpha_x^2*gamma(3/beta_x)/gamma(1/beta_x)
        ans$alpha_x_hat[t]=alpha_x
        ans$beta_x_hat[t]=beta_x

        
        beta_y=m_GGD_y$maximum
        w_y=1
        alpha_y=(sum(abs(vel_y_end_list[[t]]-w_y*v_y_neighbor)^beta_y)*beta_y/n_t)^(1/beta_y)
        ans$sigma_2_0_y_hat[t]=alpha_y^2*gamma(3/beta_y)/gamma(1/beta_y)
        ans$alpha_y_hat[t]=alpha_y
        ans$beta_y_hat[t]=beta_y

      }
    }
    
    v_mean_neighbor_x_list[[t]]=v_x_neighbor
    v_mean_neighbor_y_list[[t]]=v_y_neighbor
    
    
  }
  
  ans$v_mean_neighbor_x_list=v_mean_neighbor_x_list
  ans$v_mean_neighbor_y_list=v_mean_neighbor_y_list
  
  return(ans)
  
  
}



simulate_cells_anisotropic_with_neighbor_fast_grid<-function(cut_r,p0,v0,sigma_2_0_x_hat,sigma_2_0_y_hat,
                                                        T_sim,h,wx_hat,wy_hat,grid_boundary_info,
                                                        apolar_vicsek=T,density_x_list=NULL,density_y_list=NULL,
                                                        alpha_x_hat=NULL,alpha_y_hat=NULL,beta_x_hat=NULL,beta_y_hat=NULL,
                                                        residual_type=residual_type){

  sigma_0_x_hat=sqrt(sigma_2_0_x_hat)
  sigma_0_y_hat=sqrt(sigma_2_0_y_hat)
  
  
  n_sim=dim(p0)[1] ##fix particle
  pos_x_sim_record=matrix(NA,dim(p0)[1],T_sim)
  vel_x_sim_record=matrix(NA,dim(v0)[1],T_sim)
  pos_y_sim_record=matrix(NA,dim(p0)[1],T_sim)
  vel_y_sim_record=matrix(NA,dim(v0)[1],T_sim)
  
  pos_x_sim_record[,1]=p0[,1]
  pos_y_sim_record[,1]=p0[,2]
  vel_x_sim_record[,1]=v0[,1]
  vel_y_sim_record[,1]=v0[,2]
  
  v_mean_neighbor_x_sim_record=matrix(NA,dim(p0)[1],T_sim)
  v_mean_neighbor_y_sim_record=matrix(NA,dim(p0)[1],T_sim)
  
  theta_all_sim=matrix(NA,n_sim,(T_sim))
  theta_all_sim[,1]=atan2(  vel_y_sim_record[,1],
                            vel_x_sim_record[,1])
  
  theta_all_0_pi_sim=matrix(NA,n_sim,(T_sim))
  theta_all_0_pi_sim[,1]=  theta_all_sim[,1]
  theta_all_0_pi_sim[which( theta_all_0_pi_sim[,1]<0),1]= theta_all_0_pi_sim[which( theta_all_0_pi_sim[,1]<0),1]+pi
  record_ensemble_order=rep(NA,(T_sim))
  record_ensemble_order[1]=sum(cos(2*theta_all_0_pi_sim[,1]))/length( theta_all_0_pi_sim[,1])
  
  Lx_min=grid_boundary_info$grid_info[1]
  Lx_max=grid_boundary_info$grid_info[2]
  Ly_min=grid_boundary_info$grid_info[3]
  Ly_max=grid_boundary_info$grid_info[4]
  nx=grid_boundary_info$grid_info[5]
  ny=grid_boundary_info$grid_info[6]
  len_x=grid_boundary_info$grid_info[7]
  len_y=grid_boundary_info$grid_info[8]
  
  #m_grid_here=m_grid
  neighbor_index_list=NULL
  for(t in 1:(T_sim-1)){
    #print(t)
    
    initiate_grid_list=initiate_grid(grid_boundary_info=grid_boundary_info,pos_x=pos_x_sim_record[,t],pos_y=pos_y_sim_record[,t],
                                     vel_x=vel_x_sim_record[,t],vel_y=vel_y_sim_record[,t],neighbor_index_list=neighbor_index_list)
    m_grid_here=initiate_grid_list$m_grid
    neighbor_index_list=initiate_grid_list$neighbor_index_list
    
    
    input_here=t(cbind(pos_x_sim_record[,t],pos_y_sim_record[,t]))
    input_v_here=t(cbind(vel_x_sim_record[,t],vel_y_sim_record[,t]))
    
    if(residual_type=='Nonparametric'){
      c_prop_x=3
      c_x_sample=max(density_x_list[[t]][,2]/dnorm(density_x_list[[t]][,1],mean=0,sd=(c_prop_x*sigma_0_x_hat[t]) ))
      c_prop_y=3
      c_y_sample=max(density_y_list[[t]][,2]/dnorm(density_y_list[[t]][,1],mean=0,sd=(c_prop_y*sigma_0_y_hat[t]) ))
    }
    
    
    for(i in 1:n_sim){
      #set.seed((t-1)*n_sim+i)
      #generage noise
      z_xy=rep(NA,2)
      if(residual_type=='Gaussian'){
        z_xy[1]=sigma_0_x_hat[t]*rnorm(1) 
        z_xy[2]=sigma_0_y_hat[t]*rnorm(1) 
      }else if(residual_type=='Laplace'){
        z_xy[1]=rmLaplace(n = 1,Scatter=matrix(sigma_0_x_hat[t]^2,1,1)) #may consider rgnorm in gnorm
        z_xy[2]=rmLaplace(n = 1,Scatter=matrix(sigma_0_y_hat[t]^2,1,1)) #may consider rgnorm in gnorm
      }else if(residual_type=='GGD'){
        z_xy[1]=rgnorm(n=1, mu = 0, alpha = alpha_x_hat[t], beta = beta_x_hat[t])
        z_xy[2]=rgnorm(n=1, mu = 0, alpha = alpha_y_hat[t], beta = beta_y_hat[t])
      }else if(residual_type=='Nonparameteric'){
        accept=F
        while(accept==F){
          z_xy[1]=c_prop_x*sigma_0_x_hat[t]*rnorm(1) 
          u=runif(1)
          abs_diff=abs(z_xy[1]-density_x_list[[t]][,1])
          index_selected=which(abs_diff==min(abs(abs_diff)))
          true_density=density_x_list[[t]][index_selected,2]
          proposal_density=dnorm(density_x_list[[t]][index_selected,1],mean=0,sd=c_prop_x*sigma_0_x_hat[t])
          if(u<=true_density/(c_x_sample*proposal_density)){
            accept=T
          }
        }
        accept=F
        while(accept==F){
          z_xy[2]=c_prop_y*sigma_0_y_hat[t]*rnorm(1) 
          u=runif(1)
          abs_diff=abs(z_xy[2]-density_y_list[[t]][,1])
          index_selected=which(abs_diff==min(abs(abs_diff)))
          true_density=density_y_list[[t]][index_selected,2]
          proposal_density=dnorm(density_y_list[[t]][index_selected,1],mean=0,sd=c_prop_y*sigma_0_y_hat[t])
          if(u<=true_density/(c_y_sample*proposal_density)){
            accept=T
          }
        }
      }
 
      input_i=c(pos_x_sim_record[i,t],pos_y_sim_record[i,t])
      
      
      i_x=ceiling((pos_x_sim_record[i,t]-Lx_min)/len_x)
      i_y=ceiling((pos_y_sim_record[i,t]-Ly_min)/len_y)
      
      index_grid=(i_y-1)*nx+i_x
      d_vec_here_all=input_i-as.matrix(m_grid_here[[index_grid]]$neighbor_pos)
      
      d_here=sqrt(colSums(d_vec_here_all^2))
      if(apolar_vicsek==F){
        index_neighbor=which(d_here<cut_r)
      }else{
        index_neighbor=which(d_here<cut_r)
        index_same_v_direction=which(colSums(m_grid_here[[index_grid]]$neighbor_vel*input_v_here[,i])>=0) ##Jan 11
        index_neighbor=intersect(index_neighbor,index_same_v_direction)
        
      }
  
      v_mean_neighbor_x_sim_record[i,t]=mean(m_grid_here[[index_grid]]$neighbor_vel[1,index_neighbor])
      v_mean_neighbor_y_sim_record[i,t]=mean(m_grid_here[[index_grid]]$neighbor_vel[2,index_neighbor])
      
      vel_x_sim_record[i,t+1]=wx_hat[t]*v_mean_neighbor_x_sim_record[i,t]+z_xy[1]
      pos_x_sim_record[i,t+1]=pos_x_sim_record[i,t]+vel_x_sim_record[i,t]*h
      
      vel_y_sim_record[i,t+1]=wy_hat[t]*v_mean_neighbor_y_sim_record[i,t]+z_xy[2]
      pos_y_sim_record[i,t+1]=pos_y_sim_record[i,t]+vel_y_sim_record[i,t]*h
    
    
      ###boundary 
      if(pos_x_sim_record[i,t+1]>=Lx_max | pos_x_sim_record[i,t+1]<=Lx_min ){
        vel_x_sim_record[i,t+1] = -vel_x_sim_record[i,t+1]
        ##
        pos_x_sim_record[i,t+1]=pos_x_sim_record[i,t]#+vel_x_sim_record[i,t]*h  ##stay there don't go out
        
      }
      
      if( pos_y_sim_record[i,t+1]>=Ly_max | pos_y_sim_record[i,t+1]<=Ly_min){
        vel_y_sim_record[i,t+1] = -vel_y_sim_record[i,t+1]
        ##
        pos_y_sim_record[i,t+1]=pos_y_sim_record[i,t] #+vel_y_sim_record[i,t]*h
        
      }
      theta_here=atan2(vel_y_sim_record[i,t+1],
                       vel_x_sim_record[i,t+1])
      
      theta_all_sim[i,t+1]=c(theta_here)
      
      
    }
    
    
    theta_all_0_pi_sim[,t+1]=  theta_all_sim[,t+1]
    theta_all_0_pi_sim[which( theta_all_0_pi_sim[,t+1]<0),t+1]= theta_all_0_pi_sim[which( theta_all_0_pi_sim[,t+1]<0),t+1]+pi
    record_ensemble_order[t+1]=sum(cos(2*theta_all_0_pi_sim[,t+1]))/length( theta_all_0_pi_sim[,t+1])

    
  }
  ans_list=list()
  ans_list$pos_x_sim_record=pos_x_sim_record
  ans_list$pos_y_sim_record=pos_y_sim_record  
  ans_list$v_mean_neighbor_x_sim_record=v_mean_neighbor_x_sim_record  
  
  ans_list$v_mean_neighbor_y_sim_record=v_mean_neighbor_y_sim_record  
  
  ans_list$vel_x_sim_record=vel_x_sim_record
  ans_list$vel_y_sim_record=vel_y_sim_record  
  
  ans_list$theta_all_sim=theta_all_sim
  ans_list$theta_all_0_pi_sim=theta_all_0_pi_sim
  
  ans_list$record_ensemble_order=record_ensemble_order
  
  return(ans_list)
  
}

permutation_F_test<-function(vx, vy, permutation=T,B=10^3,record_permutation=F){
  nx=length(vx)
  ny=length(vy)
  res=list()
  var_test=var.test(vx,vy, alt="two.sided")
  F_statistics=var_test$statistic
  res$statistics=F_statistics
  
  if(!permutation){ ##no permutation, just return two sided p-value
    res$p_value=var_test$p.value
  }else{
    F_statistics_B=rep(NA, B)
    v_xy=c(vx,vy)
    n_xy=length(v_xy)
    
    for(i_B in 1:B){
      sample_index=sample(1:length(v_xy),replace=F)
      
      vx_B=v_xy[sample_index[1:(n_xy/2)]]
      vy_B=v_xy[sample_index[(n_xy/2+1):n_xy]]
      
      F_statistics_B[i_B]=var(vx_B)/var(vy_B)
    }

    lower_n=length(which(F_statistics>F_statistics_B))
    upper_n=length(which(F_statistics<F_statistics_B))
    res$p_value=min(lower_n,upper_n)*2/B
    if(record_permutation){
      res$F_statistics_B=F_statistics_B
    }
  }
  return(res)
  
}


