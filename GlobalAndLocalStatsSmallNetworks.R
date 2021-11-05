library(tidyverse)
library(igraph)
library("DirectedClustering")
library("tnet")

remove(list = ls())

########### Funciones glopbales #####

#Function to create directed weighted networks. Need a file with 3 columns"
#Source, Target and Weight
networkCreationDW<-function(x){
  el1 <-as.data.frame(x) #Verify that it's a data.frame
  #  el1<-subset(el1,source_epiunit_id!=destination_epiunit_id) #Avoid loops, self-connections
  names(el1)<-c("i","j","w")
  el1<-subset(el1,i!=j) #Avoid loops, self-connections
  #  el1<-subset(el1,Source!=Target) #Avoid loops, self-connections
  ell3<-as.matrix(el1)
  ell3[,1]=as.character(ell3[,1]) #Because the vertex IDs in this dataset are numbers, we make sure igraph knows these should be treated as characters. Otherwise, it'll create problems (see page on data import)
  ell3[,2]=as.character(ell3[,2])
  g=graph.edgelist(ell3[,1:2],directed = TRUE) #Create the directed network
  E(g)$weight=as.numeric(ell3[,3]) #weights using the number of cattle made
  g<-simplify(g) #Verify that there are not multiple connection or remaining loops
}
### Funciones para medidas globales redes pequenas
netsRealStatsSmall<-function(x){
  g<-networkCreationDW(x)
  pi<-transitivity(g,type="weighted",weights = E(g)$weight)
  is.na(pi)<-sapply(pi,is.infinite)
  #  A<-get.adjacency(g,sparse = F,attr="weight")
  #  coefAgrup<-ClustBCG(A,"directed")
  Distances<-distance_w(x,directed = T)
  reci<-data.frame(Nodes=vcount(g),Edges=ecount(g),
                   Clustering=mean(pi,na.rm=T),
                   #                   ClustCoef=coefAgrup$GlobaltotalCC,
                   GSCC=max(components(g,mode="strong")$csize),
                   GWCC=max(components(g,mode="weak")$csize),
                   Lambda1=eigen_centrality(g,directed = T,weights = E(g)$weight)$value,
                   Reciprocity=reciprocity(g),
                   SPL=mean(Distances,na.rm=T),
                   Diameter=max(Distances,na.rm=T),
                   Assort=assortativity(g, types1 = graph.strength(g), directed = T))
  reci
}
#Takes an edge list with weights (from, to, weight)
randomNetStatisticsSmall <- function(x){
  #create a random network using algorithm explained above and conver to igraph object
  prueba1<-rg_reshuffling_w(x,option = "links",directed = T) 
  esto<-rg_reshuffling_w(prueba1,option = "weights",directed = T) 
  #calculate network measures
  g<-networkCreationDW(esto)
  pi<-transitivity(g,type="weighted",weights = E(g)$weight)
  is.na(pi)<-sapply(pi,is.infinite)
  #  A<-get.adjacency(g,sparse = F,attr="weight")
  #  coefAgrup<-ClustBCG(A,"directed")
  Distances<-distance_w(esto,directed = T)
  reci<-data.frame(Nodes=vcount(g),Edges=ecount(g),
                   Clustering=mean(pi,na.rm=T),
                   GSCC=max(components(g,mode="strong")$csize),
                   GWCC=max(components(g,mode="weak")$csize),
                   #                  ClustCoef=coefAgrup$GlobaltotalCC,
                   Lambda1=eigen_centrality(g,directed = T,weights = E(g)$weight)$value,
                   Reciprocity=reciprocity(g),
                   SPL=mean(Distances,na.rm=T),
                   Diameter=max(Distances,na.rm=T),
                   Assort=assortativity(g, types1 = graph.strength(g), directed = T))
  reci
}

nodeStatsDW<-function(x){
  eig<-eigen_centrality(x,directed = T,weights = E(x)$weight)$vector
  inDegree<-degree(x,mode = "in")
  outDegree<-degree(x,mode = "out")
  inStrength<-strength(x,mode="in",weights = E(x)$weight)
  outStrength<-strength(x,mode="out",weights = E(x)$weight)
  betweenness<-betweenness(x, directed=T, weights=E(x)$weight,nobigint = TRUE, normalized = FALSE)
  inCoreness<-coreness(x,mode="in")
  outCoreness<-coreness(x,mode="out")
  centrality <- data.frame(Node = as.numeric(rownames(as.data.frame(inDegree))),
                           OutDegree    = outDegree,
                           InDegree    = inDegree,
                           OutStrength = outStrength,
                           InStrength   = inStrength,
                           OutCoreness = outCoreness,
                           InCoreness = inCoreness,
                           Betweenness  = betweenness,
                           EigenCent = eig
  )
  centrality
}

######################################
#I get the nales of all edgelists I will use
archivos0<-Sys.glob("~/Projects/FMD-Surveillance/Communities/REDESNuevo/RedesSemester/*.txt")
StatsLayersDistrictsSemester<-NULL #output data frame for global stats
for(i in archivos0){
  file0<-read.table(i,sep=" ",header = F) #for each name, I read the file
  que<-netsRealStatsSmall(file0) #Calculate the set of global measures (see function above)
  queRan<-NULL # data frame for the random equivalent set
  for (j in 1:20) { #In this case I am using 20 random equivalente networks
    #create a random equivalente network and calculate all measure as in the real
    #network (look above)
    prueba<-randomNetStatisticsSmall(file0) 
    queRan<-rbind(queRan,prueba)
  }
  #Averages, aggregation and names of each network, according to archivos0
  means<-sapply(queRan,mean) %>% round(.,4)
  sds<-sapply(queRan,sd) %>% round(.,4)
  #This changes depending on the path of the files (line 94)
  esto0<-as.vector(eval(i) %>% str_split(.,"/") %>% unlist(.))
  #This changes depending on the path of the files (line 94)
  esto1<-gsub('.{4}$', '', esto0[9]) %>% str_split("ter") %>% unlist()
  esto2<-as.numeric(esto1[2])
  #Put everything together
  esta<-data.frame(Measure=rownames(t(que)),
                   Real=t(que),RandMean=means,
                   RandSd=sds,Semester=esto2)
  rownames(StatsLayersDistrictsSemester) <- NULL
  #Add it to a data frame. Data frame with stats from real network (line 98)
  #And global stats from set of random equivalent (lines 100-115)
  StatsLayersDistrictsSemester<-rbind(StatsLayersDistrictsSemester,esta)
}

localStatsLayersDistrictsSemester<-NULL #Outout data frame for local stats (each node of each network)
for(i in archivos0){
  file0<-read.table(i,sep=" ",header = F)   #read the file
  red<-networkCreationDW(file0) #Create the network
  temp0<-nodeStatsDW(red) #Calculate local stats (look above)
  #Write the names (might change according to path of files)
  esto0<-as.vector(eval(i) %>% str_split(.,"/") %>% unlist(.))
  esto1<-gsub('.{4}$', '', esto0[9]) %>% str_split("mester") %>% unlist()
  esto2<-esto1[2]
  temp0 <- temp0 %>% mutate(Semester=esto2)
  #Aggregate in a data frame
  localStatsLayersDistrictsSemester<-rbind(localStatsLayersDistrictsSemester,temp0)
}
