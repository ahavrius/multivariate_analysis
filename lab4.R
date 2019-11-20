library(readxl)
library(rgl)
library(kernlab)

hier_clust<-function(data_tmp, dist.type="euclidean", aglo.method="single", label="none"){
  dist.data<-dist(data_tmp, method=dist.type)
  hc<-hclust(dist.data, method = aglo.method)
  plot(as.dendrogram(hc), leaflab=label)
  return(hc)
}

prin.comp<-function(data_tmp){
  prin.res<-princomp(data_tmp)
  plot(prin.res)
  summary(prin.res)
  plot(prin.res$scores[, c(1, 2)])
  plot3d(prin.res$scores[, c(1:3)])
  return(prin.res)  
}

plot_3d<-function(xyz, color, dot.size=0.2, title=''){
  pal<- sample(rainbow(15))
  plot3d(xyz, col=pal[color], cex=dot.size, main=title)
}

make_all_job<-function(data_tmp, k_clust=2, data_plot=data_tmp){
  hc<-hier_clust(data_tmp)
  hc.res<-cutree(hc, k=k_clust)
  #plot_3d(data_plot, hc.res, title='hier tree')
  spec.res<-specc(data_tmp, centers=k_clust)
  plot_3d(data_plot, spec.res, title='spect res')
  return(spec.res)
}

dataset <- read.table('/home/asta/Documents/Multivariables/F3c.txt')
real_data<- read_excel("data_by_counties.xls")
num_data<-real_data[,3:10]
num_data<-scale(num_data)


pairs(dataset)
prin.res.data1<-prin.comp(dataset)

make_all_job(dataset, data_plot=prin.res.data1$scores[, c(1:3)])
make_all_job(prin.res.data1$scores[, c(1:3)])


pairs(num_data)
prin.res.data2<-prin.comp(num_data)

make_all_job(num_data, data_plot=prin.res.data2$scores[, c(1:3)], k=5)
make_all_job(prin.res.data2$scores[, c(1:3)], k=5)
