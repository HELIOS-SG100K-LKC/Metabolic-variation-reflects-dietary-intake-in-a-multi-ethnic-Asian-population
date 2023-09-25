# To visualise partial correlation networks between foods, beverages and dietary biomarker panels using fruchterman-reingold layout
library(tidyverse)
library(ggraph)
library(readxl)

# Read files of example dataset
data_network=read.csv("./data/data_network.csv",row.names=1)

graph=as_tbl_graph(data_network,directed=FALSE)
graph %>% activate(nodes) %>% as.data.frame() %>% nrow()

# Create edge group to distinguish between positive and negative coefficients
data_network=data_network %>%
  mutate(r_dir=ifelse(r>=0,"pos","neg"))
data_network$r_dir=factor(data_network$r_dir, levels=c("pos","neg"))

# Visualise network
set.seed(188)
z=c("food1","food1","food1","food2","food3","food3","food3","food2","food2","food3")
s=c("triangle","circle","circle","triangle","triangle","circle","circle","circle","circle","circle")
ggraph(graph,layout='fr')+
  geom_edge_link(aes(width=data_network$r,colour=factor(data_network$r_dir)),alpha=0.3) + 
  scale_edge_width(range=c(0.1,4)) +
  scale_edge_colour_manual(values=c("pos"="#D79B9B","neg"="royalblue4"),name="Direction")+
  geom_node_point(size=6,aes(colour=factor(z),shape=s),alpha=0.8)+
  scale_colour_manual(values=c("food1"="dodgerblue2","food2"="deepskyblue2","food3"="darkblue"),
                      name="Node type")+
  theme_graph()

