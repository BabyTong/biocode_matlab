##############################################################
#函数随机游走
##############################################################
##RandomWalk function
library(igraph)

#子程序 

RandomWalk2igraph<-function(igraphM,VertexWeight,EdgeWeight=TRUE,gamma=0.7){
  if(EdgeWeight==TRUE){
    adjM<-get.adjacency(igraphM,attr="weight") # convert igraph object to a weight matrix
  }
  if(EdgeWeight==FALSE){
    adjM<-get.adjacency(igraphM) # convert igraph object to a conventional matrix
  }
  res<-rw(adjM,VertexWeight,gamma)
  #        print(Sys.time());
  return(drop(res))
}
rw<-function(W,p0,gamma) {
  ## perform a random walk;  
  p0<-t(p0)
  p0 <- p0/sum(p0)
  PT <- p0
  k <- 0
  delta <- 1
  Ng <- dim(W)[2]
  for (i in 1:Ng) {
    sumr<-sum(W[i,])
    if(sumr==0){
      W[i,] <-numeric(length=length(W[i,]))
    }
    if(sumr>0){
      W[i,] <- W[i,]/sum(W[i,])
    }
  }
  W<-as.matrix(W)
  W <- t(W)
  
  while(delta>1e-6) {
    PT1 <- (1-gamma)*W
    PT2 <- PT1 %*% t(PT)
    PT3 <- (gamma*p0)
    PT4 <- t(PT2) + PT3
    delta <- sum(abs(PT4 - PT))
    PT <- PT4
    k <- k + 1
  }
  PT<-t(PT)
  rownames(PT)<-NULL
  return(PT)
}

#主程序
#文件net包含三列：前两列代表你的网络，第三列是边的权重
#seed_id是你筛选的种子节点
#score是随机游走运行结束后所有节点的最终的打分

library(igraph)
edge1<-as.character(net[,1])
edge2<-as.character(net[,2])
weight<-as.numeric(net[,3])
relationship<-data.frame(edge1,edge2,weight)
point_0<-data.frame(name=unique(union(edge1,edge2)),size=0)
rownames(point_0)<-point_0[,1]
g1<-graph.data.frame(relationship,directed = FALSE,vertices = point_0)
adjweight<-get.adjacency(g1,sparse = T,attr = 'weight')
adjnotweight<-get.adjacency(g1,sparse = T)

# adjweight 有权邻接矩阵
# adjnotweight 无权邻接矩阵

point<-point_0
seed_id<-seed_id
point[seed_id,]$size<-1
score<-rw(W =adjweight,p0 = point[,2],gamma = 0.85)

#score是随机游走运行结束后所有节点的最终的打分