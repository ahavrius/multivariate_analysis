library(readxl)
require(CCA)
library(cluster)


canon_comp<-function(res.cl, data_tmp=dataset, which = c(1,3), dot.size=0.2){
  pal<-rainbow(max(res.cl)+1)
  k_tmp<-length(levels(as.factor(res.cl)))
  n<-nrow(data_tmp)
  C_matrix<-matrix(data=as.numeric(rep(res.cl,k_tmp)==rep(1:k_tmp, each=n)), ncol = k_tmp, nrow = n)
  cc_res<-rcc(data_tmp, C_matrix, 0.1, 0.1)
  return(cc_res$scores$xscores)
}

plot_2d<-function(x, y, color, dot.size=0.2){
  pal<-rainbow(15)
  plot(x, y, col=pal[color], cex=dot.size)
}

scaling<-function(data_temp, which.dist){
  dist_<-dist(data_temp, which.dist)
  fit<-cmdscale(dist_, eig = TRUE)
  return(fit$points)  
}

scaling_eig<-function(data_temp, which.dist){
  dist_<-dist(data_temp, which.dist)
  fit<-cmdscale(dist_, eig = TRUE)
  return(fit$points * fit$eig)  
}

dataset <- read.table('/home/asta/Documents/Multivariables/mult3.txt')
real_data<- read_excel("data_by_counties.xls")
num_data<-real_data[,3:10]
num_data<-scale(num_data)

km.res <- kmeans(dataset, 12, nstart = 25)
pam.res <- pam(dataset, 12)

dist.type<-c("euclidean", "maximum", "manhattan")
for (dist in dist.type){
  xy_scale<-scaling(dataset, dist)
  #xy_scale_eig<-scaling_eig(dataset, dist)
  xy_canon = canon_comp(km.res$cluster, data_tmp = xy_scale)
  
  plot_2d(xy_scale[,1], xy_scale[,2], km.res$cluster)
  #plot_2d(xy_scale_eig[,1], xy_scale_eig[,2], km.res$cluster)
  plot_2d(xy_canon[,1], xy_canon[,2], km.res$cluster)
}

for (dist in dist.type){
  xy_scale<-scaling(dataset, dist)
  #xy_scale_eig<-scaling_eig(dataset, dist)
  xy_canon<-canon_comp(pam.res$cluster, data_tmp = xy_scale)
  
  plot_2d(xy_scale[,1], xy_scale[,2], pam.res$clustering)
  #plot_2d(xy_scale_eig[,1], xy_scale_eig[,2], pam.res$clustering)
  plot_2d(xy_canon[,1], xy_canon[,2], pam.res$clustering)
}


km.res <- kmeans(num_data, 6, nstart = 25)
pam.res <- pam(num_data, 7)

for (dist in dist.type){
  xy_scale<-scaling(num_data, dist)
  xy_canon = canon_comp(km.res$cluster, data_tmp = xy_scale)
  plot_2d(xy_scale[,1], xy_scale[,2], km.res$cluster, dot.size = 1)
  plot_2d(xy_canon[,1], xy_canon[,2], km.res$cluster, dot.size = 1)
  
  #xy_scale_eig<-scaling_eig(num_data, dist)
  #xy_canon_eig = canon_comp(km.res$cluster, data_tmp = xy_scale_eig)
  #plot_2d(xy_scale_eig[,1], xy_scale_eig[,2], km.res$cluster, dot.size=1)
  #plot_2d(xy_canon_eig[,1], xy_canon_eig[,2], km.res$cluster, dot.size = 1)
}

for (dist in dist.type){
  xy_scale<-scaling(num_data, dist)
  xy_canon<-canon_comp(pam.res$cluster, data_tmp = xy_scale)
  plot_2d(xy_scale[,1], xy_scale[,2], pam.res$clustering, dot.size = 1)
  plot_2d(xy_canon[,1], xy_canon[,2], pam.res$clustering, dot.size = 1)
  
  #xy_scale_eig<-scaling_eig(num_data, dist)
  #xy_canon_eig<-canon_comp(pam.res$cluster, data_tmp = xy_scale_eig)
  #plot_2d(xy_scale_eig[,1], xy_scale_eig[,2], pam.res$clustering, dot.size=1)
  #plot_2d(xy_canon_eig[,1], xy_canon_eig[,2], pam.res$clustering, dot.size = 1)
}
