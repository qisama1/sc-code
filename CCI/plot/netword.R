library(tidyverse)
library(RColorBrewer)
library(scales)
library(igraph)

pvalues=read.table("./test/pvalues.txt",header = T,sep = "t",stringsAsFactors = F)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
colnames(statdf)=c("number")

statdf$indexb=str_replace(rownames(statdf),"^.*\.","")
statdf$indexa=str_replace(rownames(statdf),"\..*$","")
rankname=sort(unique(statdf$indexa)) 

A=c()
B=c()
C=c()
remaining=rankname
for (i in rankname[-6]) {
  remaining=setdiff(remaining,i)
  for (j in remaining) {
    count=statdf[statdf$indexa == i & statdf$indexb == j,"number"]+
      statdf[statdf$indexb == i & statdf$indexa == j,"number"]
    A=append(A,i)
    B=append(B,j)
    C=append(C,count)
  }
}

statdf2=data.frame(indexa=A,indexb=B,number=C)
statdf2=statdf2 %>% rbind(statdf[statdf$indexa==statdf$indexb,c("indexa","indexb","number")])
statdf2=statdf2[statdf2$number > 0,] #过滤掉值为0的观测

#设置节点和连线的颜色
color1=c("#8DD3C7", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD")
names(color1)=rankname
color2=colorRampPalette(brewer.pal(9, "Reds")[3:7])(20) #将颜色分成多少份，取决于互作关系数目的最大值
names(color2)=1:20 #每一份颜色用对应的数字命名

#做网络图
##下面的四行代码相对固定
net <- graph_from_data_frame(statdf2[,c("indexa","indexb","number")])
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = order(membership(group)))

E(net)$width <- E(net)$number / 2 #将数值映射到连线的宽度，有时还需要微调，这里除以2就是这个目的
E(net)$color <- color2[as.character(ifelse(E(net)$number > 20,20,E(net)$number))] #用前面设置好的颜色赋给连线，颜色深浅对应数值大小
E(net)$label = E(net)$number #连线的标注
E(net)$label.color <- "black" #连线标注的颜色
V(net)$label.color <- "black" #节点标注的颜色
V(net)$color <- color1[names(V(net))] #节点的填充颜色，前面已经设置了；V(net)返回节点信息

#调整节点位置的线条角度
##如果没有这两行代码，节点位置的圆圈是向右的
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]

#pdf("interaction.num.3.pdf",width = 6,height = 6)
plot(net,
     edge.arrow.size = 0, #连线不带箭头
     edge.curved = 0, #连线不弯曲
     vertex.frame.color = "black", #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = 30) #节点大小
#dev.off()