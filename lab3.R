library(dendextend)
library(readxl)

cut_plot_subtree<-function(hc, h){
  dendro<-as.dendrogram(hc)
  upper_part<-cut(dendro, h=h)$upper
  lower_part<-cut(dendro, h=h)$lower
  
  num_iteam<-length(lower_part)
  label<-letters[1:num_iteam]
  labels(upper_part)<-label
  
  plot(upper_part, main = "upper subtree of cut")
  for (i in 1:num_iteam){
    labels(lower_part[[i]]) <-NULL
    plot(lower_part[[i]], main = label[i])
  }
    
}


coph_dist<-function(data_tmp, dist.type="euclidean", aglo.method="average"){
  dist.data<-dist(data_tmp, method=dist.type)
  hc<-hclust(dist.data, method = aglo.method)
  coph<-cophenetic(hc)
  cor.coph<-cor(coph, dist.data)
  cat("distance", dist.type, "aglo method", aglo.method, "cophenetic cor", cor.coph, "\n")
  return(hc)
}


dataset <- read.table('/home/asta/Documents/Multivariables/mult3.txt')
real_data<- read_excel("data_by_counties.xls")
num_data<-real_data[,3:10]
num_data<-scale(num_data)

dist.types<-c("euclidean", "maximum", "manhattan")
aglo.methods<-c("single", "average", "complete")

for (aglo.method in aglo.methods){
  for (dist.type in dist.types){
    hc<-coph_dist(dataset, dist.type, aglo.method)
    plot(hc, labels = FALSE)
  }
}

cut_plot_subtree(hc, h=50)


for (aglo.method in aglo.methods){
  for (dist.type in dist.types){
    hc<-coph_dist(num_data, dist.type, aglo.method)
    plot(hc, labels = FALSE)
  }
}

#hc<-coph_dist(num_data)
plot(hc, labels = FALSE)
cut_plot_subtree(hc, h=7)
clustering<-cutree(hc, h=7)
real_data$tree<-clustering

