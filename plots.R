library(ggplot2)
library(scales)
# theme_set(theme_linedraw())
# theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

########## load dataset
file_number=441
index_frame=60

#from real experiment
n_record=as.vector(as.matrix(read.csv(paste0('csv_variables/cell_count_real_',file_number,'.csv'))))
order_param=as.vector(as.matrix(read.csv(paste0('csv_variables/order_param_real_',file_number,'.csv'))))
sigma_2_x=as.vector(as.matrix(read.csv(paste0('csv_variables/sigma_2_x_real_',file_number,'.csv'))))
sigma_2_y=as.vector(as.matrix(read.csv(paste0('csv_variables/sigma_2_y_real_',file_number,'.csv'))))
mean_abs_loss_x=as.vector(as.matrix(read.csv(paste0('csv_variables/mean_abs_loss_x_real_',file_number,'.csv'))))
mean_abs_loss_y=as.vector(as.matrix(read.csv(paste0('csv_variables/mean_abs_loss_y_real_',file_number,'.csv'))))

#data of time=60 (hour = 2)
data_time=read.csv(paste0('csv_variables/data_at_time_',index_frame,'_exp_',file_number,'.csv'))

#marginal
velocity_density_mat=as.matrix(read.csv(paste0('csv_variables/marginal_velocity_density_time_',index_frame,'_experiment_',file_number,'.csv')))

#F test
F_test_statistics_p_value_record=as.matrix(read.csv(paste0('csv_variables/F_test_statistics_p_value_record_experiment',file_number,'.csv')))
F_test_sample_B_record=as.matrix(read.csv(paste0('csv_variables/F_test_sample_B_record_experiment',file_number,'.csv')))

#normality test
record_normality_test_x=as.matrix(read.csv(paste0('csv_variables/record_normality_test_x_experiment_',file_number,'.csv')))
record_normality_test_y=as.matrix(read.csv(paste0('csv_variables/record_normality_test_y_experiment_',file_number,'.csv')))

#correlation
correlation_record=as.matrix(read.csv(paste0('csv_variables/correlation_record_experiment_',file_number,'.csv')))

#weight
weight_record=as.matrix(read.csv(paste0('csv_variables/weight_interval_est_mat_Laplace_',file_number,'.csv')))

#from simulation
sim_set=as.matrix(read.csv("csv_variables/simul_setting.csv"))
order_param_rec=as.matrix(read.csv(paste0('csv_variables/order_param_rec_',file_number,'.csv'),row.names = 1))
mean_abs_loss_x_rec=as.matrix(read.csv(paste0('csv_variables/mean_abs_loss_x_rec_',file_number,'.csv'),row.names = 1))
mean_abs_loss_y_rec=as.matrix(read.csv(paste0('csv_variables/mean_abs_loss_y_rec_',file_number,'.csv'),row.names = 1))

#residual
velocity_conditional_density_Gaussian=as.matrix(read.csv(paste0('csv_variables/conditional_velocity_density_time_',index_frame,'_experiment_',file_number,'_residual_Gaussian.csv')))
velocity_conditional_density_Laplace=as.matrix(read.csv(paste0('csv_variables/conditional_velocity_density_time_',index_frame,'_experiment_',file_number,'_residual_Laplace.csv')))
velocity_conditional_density_GGD=as.matrix(read.csv(paste0('csv_variables/conditional_velocity_density_time_',index_frame,'_experiment_',file_number,'_residual_GGD.csv')))

#phase diagram
order_diagram_record=as.matrix(read.csv(paste0('csv_variables/phase_diagram_sim_experiment_',file_number,'.csv')))
w_sigma_0_diagram=as.matrix(read.csv(paste0('csv_variables/w_sigma_0_diagram_sim_experiment_',file_number,'.csv')))
w_sigma_0_order_diagram_data=as.matrix(read.csv(paste0('csv_variables/w_sigma_0_order_diagram_data_experiment_',file_number,'.csv')))


T_time=length(order_param)
# data_plot=data.frame(Time=1:T_time/3,n_record=n_record,
#                      order_param=order_param,
#                      sigma_2_x=sigma_2_x,sigma_2_y=sigma_2_y)


#### start to plot

####Figure 2
#(a)
par(mgp=c(2,1,0), mar=c(3,4,2,4.5)+.1)
barplot(c(0,n_record-2600),col='#edf1df',border='#edf1df',yaxt = "n",xaxt = "n",
        xlab="Time (hr)",ylab="",
        ylim=c(0,400),width=1/3,space=0)#ylim=c(2600,3000),
yticks <- seq(0L,400L,100L)
for (i in 1:length(yticks)){
  axis(4, at=yticks[i], labels=bquote(.(yticks[i]+2600)),cex.lab=1.2,las=1,col="#d9531a",col.axis="#d9531a")
}
xticks <- seq(0,40,10)
for (i in 1:length(xticks)){
  axis(1, at=xticks[i], labels=bquote(.(xticks[i])),cex.lab=1.2,las=1)
}
axis(1, at=c(-10,40),labels=c(),las=0)
axis(3, at=c(-10,40),labels=c(),las=0)
axis(2, at=c(-50,500),labels=c(),las=0,col="#1a96e3",col.axis="#1a96e3")
axis(4, at=c(-50,500),labels=c(),las=0,col="#d9531a",col.axis="#d9531a")
ylim_count=c(0,400)
ylim_sigma=c(0,7.2)*(10^(-3))
b = diff(ylim_count)/diff(ylim_sigma)
a = ylim_count[1] - b*ylim_sigma[1]
points(1:T_time/3,a+sqrt(sigma_2_x)*b,pch=21,bg='#eec164',lwd=.5,cex=1.5)
points(1:T_time/3,a+sqrt(sigma_2_y)*b,pch=24,bg='#46ada7',lwd=.5,cex=1.5)
yticks <- a+seq(0L,6L,2L)*10^(-3)*b
for (i in 1:length(yticks)){
  axis(2, at=yticks[i], labels=bquote(.((yticks[i]-a)/b)),cex.lab=1.2,las=1,col="#1a96e3",col.axis="#1a96e3")
}
text(par("usr")[2]*1.2,mean(par("usr")[3:4])*1.1, "Cell count", srt = -90, xpd = TRUE, pos = 4,col="#d9531a")
text(par("usr")[1]*7,mean(par("usr")[3:4])*.76, expression(paste(sigma[x],', ',sigma[y]," (",mu,m,"/",s,")")), srt = 90, xpd = TRUE, pos = 4,col="#1a96e3")
#(b) 
ylim_hist=c(0,200)
ylim_prob=c(0,150)
b = diff(ylim_hist)/diff(ylim_prob)
a = ylim_hist[1] - b*ylim_prob[1]
hist(data_time$vel_x,breaks = 100,ylim=c(0,200),col="#f4daa5",border="#f4daa5",main="",xlab=expression(paste(v[x]," (",mu,m,"/",s,")")))
lines(velocity_density_mat[,'input_x'],a+velocity_density_mat[,'density_norm_x']*b,lty=2)
lines(velocity_density_mat[,'input_x'],a+velocity_density_mat[,'density_laplace_x']*b,col="#d9531a")
axis(1, at=c(-1,1),labels=c(),las=0)
axis(3, at=c(-1,1),labels=c(),las=0)
axis(2, at=c(-50,500),labels=c(),las=0)
axis(4, at=c(-50,500),labels=c(),las=0,col="#d9531a",col.axis="#d9531a")
yticks <- a+seq(0L,150L,50L)*b
for (i in 1:length(yticks)){
  axis(4, at=yticks[i], labels=bquote(.((yticks[i]-a)/b)),cex.lab=1.2,las=1,col="#d9531a",col.axis="#d9531a")
}
text(par("usr")[2]*1.3,mean(par("usr")[3:4])*1.4, "Probability density", srt = -90, xpd = TRUE, pos = 4,col="#d9531a")
#(c)
ylim_hist=c(0,320)
ylim_prob=c(0,250)
b = diff(ylim_hist)/diff(ylim_prob)
a = ylim_hist[1] - b*ylim_prob[1]
hist(data_time$vel_y,breaks = 120,ylim=c(0,320),col="#95cecb",border="#95cecb",main="",xlab=expression(paste(v[y]," (",mu,m,"/",s,")")))
lines(velocity_density_mat[,'input_y'],a+velocity_density_mat[,'density_norm_y']*b,lty=2)
lines(velocity_density_mat[,'input_y'],a+velocity_density_mat[,'density_laplace_y']*b,col="#d9531a")
axis(1, at=c(-1,1),labels=c(),las=0)
axis(3, at=c(-1,1),labels=c(),las=0)
axis(2, at=c(-50,500),labels=c(),las=0)
axis(4, at=c(-50,500),labels=c(),las=0,col="#d9531a",col.axis="#d9531a")
yticks <- a+seq(0L,250L,50L)*b
for (i in 1:length(yticks)){
  axis(4, at=yticks[i], labels=bquote(.((yticks[i]-a)/b)),cex.lab=1.2,las=1,col="#d9531a",col.axis="#d9531a")
}
text(par("usr")[2]*1.3,mean(par("usr")[3:4])*1.4, "Probability density", srt = -90, xpd = TRUE, pos = 4,col="#d9531a")
#(d)
plot(1:T_time/3,order_param,pch=22,bg='#1a96e3',cex=1.5,lwd=.5,
     xlab='Time (hr)',ylab='Order Parameter S',xlim=c(0,40))
#(e) 
ggplot(data_time,aes(theta)) +
  geom_histogram(binwidth = pi/10,color='#1a96e3',fill='#a6cffe') + 
  coord_polar(start = pi/2,direction=-1)+
  scale_x_continuous(breaks = rev(seq(pi,-pi+.00001,-pi/6)),minor_breaks=NULL,
                     labels = expression(paste(-5,pi,'/',6),paste(-2,pi,'/',3),paste(-pi,'/',2),paste(-pi,'/',3),paste(-pi,'/',6),0,
                                         paste(pi,'/',6),paste(pi,'/',3),paste(pi,'/',2),paste(2,pi,'/',3),paste(5,pi,'/',6),paste(pi)))+
  theme_linedraw()


####Figure 3
plot_time=60
hist(F_test_sample_B_record[plot_time,],breaks=40,ylim=c(0,800),ylab="Count",xlab="F-test statistics",
     col="#cad5a1",border="#cad5a1",main="",
     xlim=range(F_test_sample_B_record[plot_time,],F_test_statistics_p_value_record[plot_time,1]))
axis(1, at=c(-10,5),labels=c(),las=0)
axis(3, at=c(-10,5),labels=c(),las=0)
axis(2, at=c(-50,1000),labels=c(),las=0)
axis(4, at=c(-50,1000),labels=c(),las=0)
points(F_test_statistics_p_value_record[plot_time,1],y=(par("usr")[3]),col='#d9531a',pch=17,cex=1.5,xpd = TRUE)
#abline(v=F_test_statistics_p_value_record[plot_time,1],lty=2,col='#d9531a')  



####Figure 4
#par(mgp=c(3,1,0), mar=c(3,4,2,1)+.1)
plot(1:T_time/3,record_normality_test_x[,2],log="y",ylim=10^c(-50,-10),pch=21,bg='#eec164',
     xlab="Time (hr)",ylab='p-value',yaxt = "n",cex=1.5,lwd=.5)
points(1:T_time/3,record_normality_test_y[,2],pch=24,bg='#46ada7',cex=1.5,lwd=.5)
yticks <- seq(-50,-10,10)
for (i in 1:length(yticks)){
  axis(2, at=10^yticks[i], labels=bquote(10^.(yticks[i])),cex.lab=1.2,las=1)
}
###Figure 4 inset
qq_y=qqnorm((data_time$vel_y-mean(data_time$vel_y))/sd(data_time$vel_y),plot.it = F)
qqnorm((data_time$vel_x-mean(data_time$vel_x))/sd(data_time$vel_x),pch=21,bg='#eec164',cex=.8,lwd=.5,xlab="Theoretical Quantiles")
points(qq_y$x,qq_y$y,pch=24,bg='#46ada7',cex=.8,lwd=.5)
abline(0,1)



####Figure 5
#(b)
plot(data_time$vel_x_end,data_time$polar_mean_x,pch=20,cex=.3,col="#eec164",xlim=c(-0.025,0.025),ylim=c(-0.025,0.025),
     xlab=expression(paste(bar(v)[i][","][x]," (",t-1,")")),ylab=expression(paste(v[i][","][x]," (",t,")")))
abline(0,1,lty=2)
w_x_polar=sum(data_time$polar_mean_x*data_time$vel_x_end)/sum(data_time$polar_mean_x*data_time$polar_mean_x)
abline(0,w_x_polar,col="#d9531a",lwd=1.5)
#(c)
plot(data_time$vel_y_end,data_time$polar_mean_y,pch=20,cex=.3,col="#46ada7",xlim=c(-0.02,0.02),ylim=c(-0.02,0.02),
     xlab=expression(paste(bar(v)[i][","][y]," (",t-1,")")),ylab=expression(paste(v[i][","][y]," (",t,")")))
abline(0,1,lty=2)
w_y_polar=sum(data_time$polar_mean_y*data_time$vel_y_end)/sum(data_time$polar_mean_y*data_time$polar_mean_y)
abline(0,w_y_polar,col="#d9531a",lwd=1.5)
#(e)
plot(data_time$vel_x_end,data_time$apolar_mean_x,pch=20,cex=.3,col="#eec164",xlim=c(-0.025,0.025),ylim=c(-0.025,0.025),
     xlab=expression(paste(bar(v)[i][","][x]," (",t-1,")")),ylab=expression(paste(v[i][","][x]," (",t,")")))
abline(0,1,lty=2)
w_x_apolar=sum(data_time$apolar_mean_x*data_time$vel_x_end)/sum(data_time$apolar_mean_x*data_time$apolar_mean_x)
abline(0,w_x_apolar,col="#d9531a",lwd=1.5)
#(f)
plot(data_time$vel_y_end,data_time$apolar_mean_y,pch=20,cex=.3,col="#46ada7",xlim=c(-0.02,0.02),ylim=c(-0.02,0.02),
     xlab=expression(paste(bar(v)[i][","][y]," (",t-1,")")),ylab=expression(paste(v[i][","][y]," (",t,")")))
abline(0,1,lty=2)
w_y_apolar=sum(data_time$apolar_mean_y*data_time$vel_y_end)/sum(data_time$apolar_mean_y*data_time$apolar_mean_y)
abline(0,w_y_apolar,col="#d9531a",lwd=1.5)
#(g)
plot(1:T_time/3,correlation_record[,'apolar_vx_est'],type="n",
     ylab="Correlation",ylim=c(.2,.7),xlab="Time (hour)")
polygon(c(1:T_time/3,rev(1:T_time/3)),c(correlation_record[,'apolar_vx_lower_95'],rev(correlation_record[,'apolar_vx_upper_95'])),
        col=alpha('#faedd2',.6),border = F)
lines(1:T_time/3,correlation_record[,'apolar_vx_est'],col='#eec164',lwd=2)
polygon(c(1:T_time/3,rev(1:T_time/3)),c(correlation_record[,'apolar_vy_lower_95'],rev(correlation_record[,'apolar_vy_upper_95'])),
        col=alpha('#cbe7e4',.6),border = F)
lines(1:T_time/3,correlation_record[,'apolar_vy_est'],col='#46ada7',lwd=2)
polygon(c(1:T_time/3,rev(1:T_time/3)),c(correlation_record[,'polar_vx_lower_95'],rev(correlation_record[,'polar_vx_upper_95'])),
        col=alpha('#f3f4da',.6),border = F)
lines(1:T_time/3,correlation_record[,'polar_vx_est'],col='#d5db7e',lwd=2)
polygon(c(1:T_time/3,rev(1:T_time/3)),c(correlation_record[,'polar_vy_lower_95'],rev(correlation_record[,'polar_vy_upper_95'])),
        col=alpha('#d3efe4',.6),border = F)
lines(1:T_time/3,correlation_record[,'polar_vy_est'],col='#65c9a6',lwd=2)
#(h)
plot(1:T_time/3,weight_record[,'apolar_wx_est'],type="n",ylim=c(0.4,0.9),
     ylab=expression(paste(hat(w)[y],", ",hat(w)[x])),xlab="Time (hour)")
polygon(c(1:T_time/3,rev(1:T_time/3)),c(weight_record[,'apolar_wx_lower_95'],rev(weight_record[,'apolar_wx_upper_95'])),
        col=alpha('#faedd2',.6),border = F)
lines(1:T_time/3,weight_record[,'apolar_wx_est'],col='#eec164',lwd=2)
polygon(c(1:T_time/3,rev(1:T_time/3)),c(weight_record[,'apolar_wy_lower_95'],rev(weight_record[,'apolar_wy_upper_95'])),
        col=alpha('#cbe7e4',.6),border = F)
lines(1:T_time/3,weight_record[,'apolar_wy_est'],col='#46ada7',lwd=2)
polygon(c(1:T_time/3,rev(1:T_time/3)),c(weight_record[,'polar_wx_lower_95'],rev(weight_record[,'polar_wx_upper_95'])),
        col=alpha('#f3f4da',.6),border = F)
lines(1:T_time/3,weight_record[,'polar_wx_est'],col='#d5db7e',lwd=2)
polygon(c(1:T_time/3,rev(1:T_time/3)),c(weight_record[,'polar_wy_lower_95'],rev(weight_record[,'polar_wy_upper_95'])),
        col=alpha('#d3efe4',.6),border = F)
lines(1:T_time/3,weight_record[,'polar_wy_est'],col='#65c9a6',lwd=2)



####Figure 6
cut_r=75
##Gaussian plot
Gaussian_ind=which(sim_set[,"cut_r"]==cut_r & sim_set[,"residual_type"]=="Gaussian")
#order parameter
plot(1:T_time/3,order_param,main='Gaussian',pch=0,col='#1a96e3',cex=0.7,lwd=2,
     xlab="Time (hr)", ylab="Order Parameter S",ylim=c(0.1,0.5))
for(i in Gaussian_ind){
  if(sim_set[i,"apolar_viscek"]==T & sim_set[i,"fixed_weight"]==F){
    lines(1:T_time/3,order_param_rec[i,],col='#003fbd')
  }
  if(sim_set[i,"apolar_viscek"]==T & sim_set[i,"fixed_weight"]==T){
    lines(1:T_time/3,order_param_rec[i,],col='#869fdf',lty=6)
  }
  if(sim_set[i,"apolar_viscek"]==F & sim_set[i,"fixed_weight"]==F){
    lines(1:T_time/3,order_param_rec[i,],col='#869fdf',lty=2)
  }
  if(sim_set[i,"apolar_viscek"]==F & sim_set[i,"fixed_weight"]==T){
    lines(1:T_time/3,order_param_rec[i,],col='#869fdf',lty=3)
  }
}
#absolute deviation
plot(1:T_time/3,mean_abs_loss_x,main='Gaussian',pch=1,col="#eec164",lwd=2,cex=.7,ylim=c(0,8)*10^(-3),
     xlab="Time (hr)",ylab=expression(paste(tilde(sigma)[x],",",tilde(sigma)[y]," (",mu,m,"/",s,")")))
points(1:T_time/3,mean_abs_loss_y,pch=2,col="#46ada7",lwd=2,cex=.7)
for(i in Gaussian_ind){
  if(sim_set[i,"apolar_viscek"]==T & sim_set[i,"fixed_weight"]==F){
    lines(1:T_time/3,mean_abs_loss_x_rec[i,],col="#c84601")
    lines(1:T_time/3,mean_abs_loss_y_rec[i,],col="#006400")
  }
  if(sim_set[i,"apolar_viscek"]==T & sim_set[i,"fixed_weight"]==T){
    lines(1:T_time/3,mean_abs_loss_x_rec[i,],col="#e6a48b",lty=6)
    lines(1:T_time/3,mean_abs_loss_y_rec[i,],col="#8db187",lty=6)
  }
  if(sim_set[i,"apolar_viscek"]==F & sim_set[i,"fixed_weight"]==F){
    lines(1:T_time/3,mean_abs_loss_x_rec[i,],col="#e6a48b",lty=2)
    lines(1:T_time/3,mean_abs_loss_y_rec[i,],col="#8db187",lty=2)
  }
  if(sim_set[i,"apolar_viscek"]==F & sim_set[i,"fixed_weight"]==T){
    lines(1:T_time/3,mean_abs_loss_x_rec[i,],col="#e6a48b",lty=3)
    lines(1:T_time/3,mean_abs_loss_y_rec[i,],col="#8db187",lty=3)
  }
}

##Laplace plot
Laplace_ind=which(sim_set[,"cut_r"]==cut_r & sim_set[,"residual_type"]=="Laplace")
#order parameter
plot(1:T_time/3,order_param,main='Laplace',pch=0,col='#1a96e3',cex=0.7,lwd=2,
     xlab="Time (hr)", ylab="Order Parameter S",ylim=c(0.1,0.5))
for(i in Laplace_ind){
  if(sim_set[i,"apolar_viscek"]==T & sim_set[i,"fixed_weight"]==F){
    lines(1:T_time/3,order_param_rec[i,],col='#003fbd')
  }
  if(sim_set[i,"apolar_viscek"]==T & sim_set[i,"fixed_weight"]==T){
    lines(1:T_time/3,order_param_rec[i,],col='#869fdf',lty=6)
  }
  if(sim_set[i,"apolar_viscek"]==F & sim_set[i,"fixed_weight"]==F){
    lines(1:T_time/3,order_param_rec[i,],col='#869fdf',lty=2)
  }
  if(sim_set[i,"apolar_viscek"]==F & sim_set[i,"fixed_weight"]==T){
    lines(1:T_time/3,order_param_rec[i,],col='#869fdf',lty=3)
  }
}
#absolute deviation
plot(1:T_time/3,mean_abs_loss_x,main='Laplace',pch=1,col="#eec164",lwd=2,cex=.7,ylim=c(0,8)*10^(-3),
     xlab="Time (hr)",ylab=expression(paste(tilde(sigma)[x],",",tilde(sigma)[y]," (",mu,m,"/",s,")")))
points(1:T_time/3,mean_abs_loss_y,pch=2,col="#46ada7",lwd=2,cex=.7)
for(i in Laplace_ind){
  if(sim_set[i,"apolar_viscek"]==T & sim_set[i,"fixed_weight"]==F){
    lines(1:T_time/3,mean_abs_loss_x_rec[i,],col="#c84601")
    lines(1:T_time/3,mean_abs_loss_y_rec[i,],col="#006400")
  }
  if(sim_set[i,"apolar_viscek"]==T & sim_set[i,"fixed_weight"]==T){
    lines(1:T_time/3,mean_abs_loss_x_rec[i,],col="#e6a48b",lty=6)
    lines(1:T_time/3,mean_abs_loss_y_rec[i,],col="#8db187",lty=6)
  }
  if(sim_set[i,"apolar_viscek"]==F & sim_set[i,"fixed_weight"]==F){
    lines(1:T_time/3,mean_abs_loss_x_rec[i,],col="#e6a48b",lty=2)
    lines(1:T_time/3,mean_abs_loss_y_rec[i,],col="#8db187",lty=2)
  }
  if(sim_set[i,"apolar_viscek"]==F & sim_set[i,"fixed_weight"]==T){
    lines(1:T_time/3,mean_abs_loss_x_rec[i,],col="#e6a48b",lty=3)
    lines(1:T_time/3,mean_abs_loss_y_rec[i,],col="#8db187",lty=3)
  }
}

##GGD plot
GGD_ind=which(sim_set[,"cut_r"]==cut_r & sim_set[,"residual_type"]=="GGD")
#order parameter
plot(1:T_time/3,order_param,main='GGD',pch=0,col='#1a96e3',cex=0.7,lwd=2,
     xlab="Time (hr)", ylab="Order Parameter S",ylim=c(0.1,0.5))
for(i in GGD_ind){
  if(sim_set[i,"apolar_viscek"]==T & sim_set[i,"fixed_weight"]==F){
    lines(1:T_time/3,order_param_rec[i,],col='#003fbd')
  }
  if(sim_set[i,"apolar_viscek"]==T & sim_set[i,"fixed_weight"]==T){
    lines(1:T_time/3,order_param_rec[i,],col='#869fdf',lty=6)
  }
  if(sim_set[i,"apolar_viscek"]==F & sim_set[i,"fixed_weight"]==F){
    lines(1:T_time/3,order_param_rec[i,],col='#869fdf',lty=2)
  }
  if(sim_set[i,"apolar_viscek"]==F & sim_set[i,"fixed_weight"]==T){
    lines(1:T_time/3,order_param_rec[i,],col='#869fdf',lty=3)
  }
}
#absolute deviation
plot(1:T_time/3,mean_abs_loss_x,main='GGD',pch=1,col="#eec164",lwd=2,cex=.7,ylim=c(0,8)*10^(-3),
     xlab="Time (hr)",ylab=expression(paste(tilde(sigma)[x],",",tilde(sigma)[y]," (",mu,m,"/",s,")")))
points(1:T_time/3,mean_abs_loss_y,pch=2,col="#46ada7",lwd=2,cex=.7)
for(i in GGD_ind){
  if(sim_set[i,"apolar_viscek"]==T & sim_set[i,"fixed_weight"]==F){
    lines(1:T_time/3,mean_abs_loss_x_rec[i,],col="#c84601")
    lines(1:T_time/3,mean_abs_loss_y_rec[i,],col="#006400")
  }
  if(sim_set[i,"apolar_viscek"]==T & sim_set[i,"fixed_weight"]==T){
    lines(1:T_time/3,mean_abs_loss_x_rec[i,],col="#e6a48b",lty=6)
    lines(1:T_time/3,mean_abs_loss_y_rec[i,],col="#8db187",lty=6)
  }
  if(sim_set[i,"apolar_viscek"]==F & sim_set[i,"fixed_weight"]==F){
    lines(1:T_time/3,mean_abs_loss_x_rec[i,],col="#e6a48b",lty=2)
    lines(1:T_time/3,mean_abs_loss_y_rec[i,],col="#8db187",lty=2)
  }
  if(sim_set[i,"apolar_viscek"]==F & sim_set[i,"fixed_weight"]==T){
    lines(1:T_time/3,mean_abs_loss_x_rec[i,],col="#e6a48b",lty=3)
    lines(1:T_time/3,mean_abs_loss_y_rec[i,],col="#8db187",lty=3)
  }
}



####Figure 7
##x-direction
ind_plot=seq(1,512,6)
#(a)Gaussian
plot(velocity_conditional_density_Gaussian[ind_plot,'input_x'],velocity_conditional_density_Gaussian[ind_plot,'density_x'],log="y",
     ylim=c(1,500),pch=1,col="#eec164",cex=.9,yaxt = "n",xaxt = "n",xlim=c(-.025,.025),lwd=1.8,
     xlab=expression(paste(v[x]," (",mu,m,"/",s,")")),ylab="PDF")
yticks <- seq(0L,2L,2L)
for (i in 1:length(yticks)){
  axis(2, at=10^yticks[i], labels=bquote(10^.(yticks[i])),cex.lab=1.2,las=1)
}
xticks <- seq(-0.02,0.02,0.02)
for (i in 1:length(xticks)){
  axis(1, at=xticks[i], labels=bquote(.(xticks[i])),cex.lab=1.2,las=1)
}
lines(velocity_conditional_density_Gaussian[ind_plot,'input_x_simulated'],velocity_conditional_density_Gaussian[ind_plot,'density_x_simulated'],col='#2c070c',lwd=1.5)
#(b)Laplace
plot(velocity_conditional_density_Laplace[ind_plot,'input_x'],velocity_conditional_density_Laplace[ind_plot,'density_x'],log="y",
     ylim=c(1,500),pch=1,col="#eec164",cex=.9,yaxt = "n",xaxt = "n",xlim=c(-.025,.025),lwd=1.8,
     xlab=expression(paste(v[x]," (",mu,m,"/",s,")")),ylab="PDF")
yticks <- seq(0L,2L,2L)
for (i in 1:length(yticks)){
  axis(2, at=10^yticks[i], labels=bquote(10^.(yticks[i])),cex.lab=1.2,las=1)
}
xticks <- seq(-0.02,0.02,0.02)
for (i in 1:length(xticks)){
  axis(1, at=xticks[i], labels=bquote(.(xticks[i])),cex.lab=1.2,las=1)
}
lines(velocity_conditional_density_Laplace[ind_plot,'input_x_simulated'],velocity_conditional_density_Laplace[ind_plot,'density_x_simulated'],col='#b44719',lwd=1.5)
#(c)GGD
plot(velocity_conditional_density_GGD[ind_plot,'input_x'],velocity_conditional_density_GGD[ind_plot,'density_x'],log="y",
     ylim=c(1,500),pch=1,col="#eec164",cex=.9,yaxt = "n",xaxt = "n",xlim=c(-.025,.025),lwd=1.8,
     xlab=expression(paste(v[x]," (",mu,m,"/",s,")")),ylab="PDF")
yticks <- seq(0L,2L,2L)
for (i in 1:length(yticks)){
  axis(2, at=10^yticks[i], labels=bquote(10^.(yticks[i])),cex.lab=1.2,las=1)
}
xticks <- seq(-0.02,0.02,0.02)
for (i in 1:length(xticks)){
  axis(1, at=xticks[i], labels=bquote(.(xticks[i])),cex.lab=1.2,las=1)
}
lines(velocity_conditional_density_GGD[ind_plot,'input_x_simulated'],velocity_conditional_density_GGD[ind_plot,'density_x_simulated'],col='#0000ff',lwd=1.5)

##y-direction
#(d)Gaussian
plot(velocity_conditional_density_Gaussian[ind_plot,'input_y'],velocity_conditional_density_Gaussian[ind_plot,'density_y'],log="y",
     ylim=c(1,500),pch=2,col="#46ada7",cex=.9,yaxt = "n",xaxt = "n",xlim=c(-.025,.025),lwd=1.8,
     xlab=expression(paste(v[x]," (",mu,m,"/",s,")")),ylab="PDF")
yticks <- seq(0L,2L,2L)
for (i in 1:length(yticks)){
  axis(2, at=10^yticks[i], labels=bquote(10^.(yticks[i])),cex.lab=1.2,las=1)
}
xticks <- seq(-0.02,0.02,0.02)
for (i in 1:length(xticks)){
  axis(1, at=xticks[i], labels=bquote(.(xticks[i])),cex.lab=1.2,las=1)
}
lines(velocity_conditional_density_Gaussian[ind_plot,'input_y_simulated'],velocity_conditional_density_Gaussian[ind_plot,'density_y_simulated'],col='#2c070c',lwd=1.5)
#(e)Laplace
plot(velocity_conditional_density_Laplace[ind_plot,'input_y'],velocity_conditional_density_Laplace[ind_plot,'density_y'],log="y",
     ylim=c(1,500),pch=2,col="#46ada7",cex=.9,yaxt = "n",xaxt = "n",xlim=c(-.025,.025),lwd=1.8,
     xlab=expression(paste(v[x]," (",mu,m,"/",s,")")),ylab="PDF")
yticks <- seq(0L,2L,2L)
for (i in 1:length(yticks)){
  axis(2, at=10^yticks[i], labels=bquote(10^.(yticks[i])),cex.lab=1.2,las=1)
}
xticks <- seq(-0.02,0.02,0.02)
for (i in 1:length(xticks)){
  axis(1, at=xticks[i], labels=bquote(.(xticks[i])),cex.lab=1.2,las=1)
}
lines(velocity_conditional_density_Laplace[ind_plot,'input_y_simulated'],velocity_conditional_density_Laplace[ind_plot,'density_y_simulated'],col='#b44719',lwd=1.5)
#(f)GGD
plot(velocity_conditional_density_GGD[ind_plot,'input_y'],velocity_conditional_density_GGD[ind_plot,'density_y'],log="y",
     ylim=c(1,500),pch=2,col="#46ada7",cex=.9,yaxt = "n",xaxt = "n",xlim=c(-.025,.025),lwd=1.8,
     xlab=expression(paste(v[x]," (",mu,m,"/",s,")")),ylab="PDF")
yticks <- seq(0L,2L,2L)
for (i in 1:length(yticks)){
  axis(2, at=10^yticks[i], labels=bquote(10^.(yticks[i])),cex.lab=1.2,las=1)
}
xticks <- seq(-0.02,0.02,0.02)
for (i in 1:length(xticks)){
  axis(1, at=xticks[i], labels=bquote(.(xticks[i])),cex.lab=1.2,las=1)
}
lines(velocity_conditional_density_GGD[ind_plot,'input_y_simulated'],velocity_conditional_density_GGD[ind_plot,'density_y_simulated'],col='#0000ff',lwd=1.5)



####Figure 8
image2D(t(order_diagram_record),x=w_sigma_0_diagram[,1],y=w_sigma_0_diagram[,2], col = hcl.colors(5000, hcl.pals()[33]),#33
        xlab=expression(paste(tau[x],' / ',tau[y])),ylab=expression(paste(w[x],' / ',w[y])))
lines(w_sigma_0_order_diagram_data[,1],w_sigma_0_order_diagram_data[,2],type="b",lty=2)
points(w_sigma_0_order_diagram_data[1,1],w_sigma_0_order_diagram_data[1,2],pch=4,col='red',lwd=3,cex=1.5)
points(w_sigma_0_order_diagram_data[T_time,1],w_sigma_0_order_diagram_data[T_time,2],pch=4,col='blue',lwd=3,cex=1.5)



