#From Monte Neate-Clegg (2023-10-30) - Neate-Clegg et al. 2023 Current Biology Fig. 2:

#This is it, based on a dataframe, "phylo", containing the mean UAI of each family and reading in 500 sample trees ("alltrees"). Because I played around with colouring and scaling the bars, there are a few options there. Important thing is that the vectors match the phylogenetic order of the tree. Lemme know if you have questions

#CY: think this requires phytools and viridis

# Create a consensus tree
tree<-ls.consensus(alltrees,method='symmetric.difference')

spp<-tree$tip.label # Extract the tree tip labels (in order), use to reorder dataframe
tree$tip.label<-as.character(phylo[spp,]$Family) # Replace tip labels with family names

UAIs<-phylo[spp,]$UAI # Extract reordered UAIs
names(UAIs)<-tree$tip.label # Name elements to match tip labels
cols<-round(UAIs*100)+1 # Create integer colors based on UAI
heights<-phylo[spp,]$Bars # Extract height of bars: log1p(species richness)
names(heights)<-tree$tip.label # Name elements
cols2<-round(heights*100)-68 # Alt colors based on species richness

tiff(filename='Fig. 2.tiff',width=2000,height=2000,units="px",res=300)
par(mar=c(0,0,0,0))
plotTree.wBars(tree,x=UAIs,type='fan',fsize=.4,
               width=5,tip.labels=T,scale=20,
               col=viridis(max(cols2),option='C')[cols2],
               border=F)
ramp<-viridis(475,option='C')[round(log1p(1:227)*100)-68] # Define palette w species richness
gradientLegend(1:227,side=2,pos = c(-200,140,-185,200), # Color ramp
               color=ramp,coords=T,
               length=0.15,depth=0.025,n.seg=c(50,100,150),inside=T,cex=.5)
text(-190,200,'# Species',adj=0.5,pos=3,offset=.5,cex=.8)

arrows(x0=c(2.5,140),x1=c(2.5,120),y0=c(-210,200),y1=c(-180,180),
       length=0.15,lwd=4)
arrows(x0=c(2.5,140),x1=c(2.5,120),y0=c(-210,200),y1=c(-180,180),
       length=0.15,lwd=2,col='white')
text(140,200,'Highest UAI = 3.08',adj=0.5,pos=3,offset=.5,cex=.8)
text(5,-210,'Lowest UAI = 0.005',adj=0.5,pos=4,offset=.5,cex=.8)

dev.off()