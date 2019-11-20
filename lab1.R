library(factoextra)
require(CCA)
library(cluster)
library(fossil)
library(MASS)
library(readxl)
library(viridisLite)

dataset <- read.table('/home/asta/Documents/Multivariables/mult3.txt')

calculate_ratio_ss<-function(res, data_tmp=dataset){
  within_ss<-res$tot.withinss
  between_ss<-km.res$betweenss
  ratio_ss = between_ss/within_ss
  cat("the ratio of ss is", ratio_ss, "\n")
}

plot_scatter<-function(res.cl, data_tmp=dataset, which=c(3, 5), dot.size=0.2){
  pal<-rainbow(max(res.cl)+1)
  slice_data = data.frame(data_tmp[,which[1]], data_tmp[,which[2]])
  plot(slice_data, col=pal[res.cl], cex=dot.size,
       xlab=paste("X ", which[1]),
       ylab=paste("X ", which[2]),
       main="Scatter plot")
}

optimal_by_wss<-function(method, data_tmp=dataset){
  fviz_nbclust(data_tmp, method, method = "wss",k.max = 25)
}

optimal_by_silhouette<-function(method, data_tmp=dataset){
  fviz_nbclust(data_tmp, method, method = "silhouette",k.max = 25)
}

prin_comp_plot<-function(res.cl, which_comp, data_tmp=dataset, dot.size=0.2){
  pal<-rainbow(max(res.cl)+1)
  plot(princomp(data_tmp)$scores[,which_comp],col=pal[res.cl],cex=dot.size,
       main="Scatter plot of principal components")
}

canon_comp_plot=function(res.cl, data_tmp=dataset, which = c(1,2), dot.size=0.2){
  pal<-rainbow(max(res.cl)+1)
  k_tmp<-length(levels(as.factor(res.cl)))
  n<-nrow(data_tmp)
  C_matrix<-matrix(data=as.numeric(rep(res.cl,k_tmp)==rep(1:k_tmp, each=n)), ncol = k_tmp, nrow = n)
  cc_res<-rcc(data_tmp, C_matrix, 0.1, 0.1)
  plot(cc_res$scores$xscores[,which],col=pal[res.cl],cex=dot.size,
       xlab=paste("Comp.", which[1]),
       ylab=paste("Comp.", which[2]),
       main="Scatter plot of canonical components")
}


#for kmeans
optimal_by_wss(kmeans)
optimal_by_silhouette(kmeans)
k_kmean_best = c(3, 7, 13)
for (k in k_kmean_best){
  km.res <- kmeans(dataset, k, nstart = 25)
  calculate_ratio_ss(km.res, dataset)
  plot_scatter(km.res$cluster, dataset)
  prin_comp_plot(km.res$cluster, 1:2)
  prin_comp_plot(km.res$cluster, 2:3)
  canon_comp_plot(km.res$cluster)
}

optimal_by_wss(pam)
optimal_by_silhouette(pam)
k_pam_best = c(3, 12, 13)
for (k in k_pam_best){
  pam.res <- pam(dataset, k)
  plot_scatter(pam.res$clustering)
  prin_comp_plot(pam.res$clustering, 1:2)
  prin_comp_plot(pam.res$clustering, 2:3)
  canon_comp_plot(pam.res$clustering)
}


index<-rand.index(pam.res$clustering, km.res$cluster)
adj.index<-adj.rand.index(pam.res$clustering, km.res$cluster)
cat("rand index =",index, " adj index =", adj.index, "\n")
table(pam.res$clustering, km.res$cluster)


real_data<- read_excel("data_by_counties.xls")
num_data<-real_data[,3:10]
num_data<-scale(num_data)

optimal_by_wss(kmeans, num_data)
optimal_by_silhouette(kmeans, num_data)
k_means_opt<-c(3, 6, 10)
for (k in k_means_opt){
  km.res <- kmeans(num_data, k, nstart = 25)
  calculate_ratio_ss(km.res, num_data)
  plot_scatter(km.res$cluster, num_data, c(3, 8), dot.size=1)
  prin_comp_plot(km.res$cluster, 1:2, num_data, dot.size = 1)
  #prin_comp_plot(km.res$cluster, 2:3, num_data, dot.size = 1)
  canon_comp_plot(km.res$cluster, num_data, dot.size = 1)
}

optimal_by_wss(pam, num_data)
optimal_by_silhouette(pam, num_data)
k_pam_opt = c(2, 7, 12)
for (k in k_pam_opt){
  pam.res <- pam(num_data, k)
  plot_scatter(pam.res$clustering, num_data, dot.size = 1)
  prin_comp_plot(pam.res$clustering, 1:2, num_data, dot.size = 1)
  #prin_comp_plot(pam.res$clustering, 2:3, num_data, dot.size = 1)
  canon_comp_plot(pam.res$clustering, num_data, dot.size = 1)
}

km.res <- kmeans(num_data, 6, nstart = 25)
pam.res <- pam(num_data, 7)

index<-rand.index(pam.res$clustering, km.res$cluster)
adj.index<-adj.rand.index(pam.res$clustering, km.res$cluster)
cat("rand index =",index, " adj index =", adj.index, "\n")
table(pam.res$clustering, km.res$cluster)

real_data$kmeans = km.res$cluster
real_data$pam = pam.res$clustering
