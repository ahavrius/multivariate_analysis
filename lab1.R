library(factoextra)
require(CCA)
library(cluster)
library(fossil)
library(MASS)


dataset <- read.table('/home/asta/Documents/Multivariables/mult3.txt')

calculate_ratio_ss<-function(res, data_tmp=dataset){
  within_ss<-res$tot.withinss
  between_ss<-km.res$betweenss
  ratio_ss = between_ss/within_ss
  cat("the ratio of ss is", ratio_ss)
  pal<-c(1:9)
  plot(data_tmp[,4], data_tmp[,5],cex=0.2,col=pal[res$cluster])
}

plot_by_method<-function(method, data_tmp=dataset){
  fviz_nbclust(data_tmp, method, method = "wss",k.max = 20)
  fviz_nbclust(data_tmp, method, method = "silhouette",k.max = 20)
}

prin_comp_plot<-function(res, which_comp, data_tmp=dataset){
  pal<-c("black","red","blue","green","magenta","chocolate", "darkblue","darkred","aquamarine","grey")
  plot(princomp(data_tmp)$scores[,which_comp],col=pal[res$cluster],cex=0.2)
}

canon_comp_plot=function(res_cluster, data_tmp=dataset){
  n<-nrow(data_tmp)
  k_tmp<-length(levels(as.factor(res_cluster)))
  C_matrix<-matrix(data=as.numeric(rep(res_cluster,k_tmp)==rep(1:k_tmp, each=n)), ncol = k_tmp, nrow = n)
  cc_res<-rcc(data_tmp, C_matrix, 0.1, 0.1)
  plot(cc_res$scores$xscores[,2:3],col=pal[res_cluster],cex=0.2)
  
}


#for kmeans
plot_by_method(kmeans)
k_kmean_best = c(3, 7, 13)
for (k in k_kmean_best){
  km.res <- kmeans(dataset, k, nstart = 25)
  calculate_ratio_ss(km.res, dataset)
  prin_comp_plot(km.res, 2:3)
}


pam.res<-pam(dataset, 13)

plot_by_method(pam)



k_pam_best = c(3, 7, 13)
for (k in k_pam_best){
  pam.res <- pam(dataset, k)
  calculate_ratio_ss(pam.res, dataset)
  #prin_comp_plot(pam.res$clustering, 2:3)
}



index<-rand.index(pam.res$clustering, km.res$cluster)
adj.index<-adj.rand.index(pam.res$clustering, km.res$cluster)

table(pam.res$clustering, km.res$cluster)
