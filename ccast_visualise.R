## This file has implemented the methods to run a ccast tree clustering

ccast_biaxialplot <-
function(optccastreeoutput,ylabel,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ##pdf(file="Biaxial tree node plots.pdf",width=10, height=10)
    ###tiff("Biaxial tree node plots%03d.tiff", height = 10, width = 10,units = 'cm',pointsize = 6,compression = "lzw", type="tifflzw", res=300)
    tiff(paste0(outputDir, "/Biaxial_tree_node_plots%03d.tiff"), height = 10, width = 10,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    cex <- 1.6
    tree33<-optccastreeoutput$tree
    dD33=optccastreeoutput$puredata[[length(optccastreeoutput$puredata)]]
    inner <- party:::nodes(tree33, c(1:max(unique(where(tree33))))[-unique(where(tree33))])
    ##layout(matrix(1:length(inner), ncol = round(sqrt(length(inner)),0)))
    ###par(mfrow=c(round(s/2)+1,2))
    out <- sapply(inner, function(i) {
        splitstat <- i$psplit$splitstatistic
        w1=which(i$weights==1)
        x <- dD33[[i$psplit$variableName]][w1]
        y<-dD33[[ylabel]][w1]
        
        if (bandwidth.nrd(x)>0&bandwidth.nrd(y)>0) {
        xydens <- kde2d(x,y,n=25,lims = c(range(x), range(y)))
        }else {
            print("Bandwidths=0")
           print(c(bandwidth.nrd(x), bandwidth.nrd(y)))
           xydens <- kde2d(x,y,h=0.5,n=25,lims = c(range(x), range(y)))
        }
        plot(x,y, main = paste("Node", i$nodeID,sep=" "),
        xlab = i$psplit$variableName, ylab = ylabel, ylim = c(0, max(y)+5),
        cex.axis = cex, cex.lab = cex, cex.main = cex,cex=0.2, pch=20,col="dark blue")
        abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
        #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
        contour(xydens,nlevels=10,add=T,col="dark orange")
    })
    dev.off();
    
}

#### biaxial plot for intermediate cels ####


##### plot biaxial tsne plots ####
ccast_tsne_plot1 <-function(x,y,outputDir, file="Biaxial T-SNE plots.tiff") {
    tiff(paste0(outputDir, file), height = 10, width = 10,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    cex <- 1.8
        xydens <- kde2d(x,y,n=100)
        plot(x,y, main ="2D T-SNE",
        xlab = "t-SNE1", ylab = "t-SNE2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = cex, cex.lab = cex, cex.main = cex,cex=0.2, pch=20,col="dark blue")
        ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
        #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
        contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    dev.off();

}
#####
##### plot biaxial pca plots ####
ccast_pca_plot1 <-function(x,y,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ##pdf(file="Biaxial T-SNE plots.pdf",width=10, height=10)
    tiff(paste0(outputDir, "Biaxial PCA plots.tiff"), height = 10, width = 10,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    cex <- 1.8
    xydens <- kde2d(x,y,n=100)
    plot(x,y, main ="2D PCA",
    xlab = "PC1", ylab = "PC2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = cex, cex.lab = cex, cex.main = cex,cex=0.2, pch=20,col="dark blue")
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    dev.off();
    
}
#####

##### plot biaxial tsne plots for all files ####
ccast_tsne_plot2 <-function(x,y,pop,poplabel,antibody,dat, outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ##pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir,"Biaxial T-SNE plots by subpopulations1.tiff"), height = 10, width = 10,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    cex_site<-c(0.2,0.2,0.8,0.6)[factor(pop)]
     pch_site<-c(20,20,19,18)[factor(pop)]
    plot(x,y, main ="2D T-SNE",
    xlab = "V1", ylab = "V2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8))
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*10)/100)]
        
        cex_site<-c(0.2,0.2,0.8,0.6)[factor(pop)]
        pch_site<-c(20,20,19,18)[factor(pop)]
        plot(x,y, main =antibody[r],
        xlab = "V1", ylab = "V2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
    }
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    dev.off();
    
}
##### plot biaxial tsne plots for all files ####
ccast_tsne_plot3 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste(outputDir,"/Biaxial T-SNE plots by subpopulations1.tiff"), height = 12, width = 12,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    cols[nn-1]="#FAA9A5"
    cols_t1<-cols[pop]
    cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D T-SNE",
    xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff("Biaxial T-SNE plots by subpopulations2.tiff", height = 12, width = 20,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,4))
    ss=c(1,3,4,2,2,5,6,7)
    for ( s in 1:(nn+1)) {
        s1=ss[s]
        pop1=ifelse(pop==s1,s1,0)
        nn1=length(unique(pop))
        cols<-brewer.pal(n=nn1,name="Set1")
        ##cols[nn-1]="#F9E429"
        cols[nn-1]="#FAA9A5"
        cols_t1<-rep(cols[s1],length(pop))
        cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
        plot(x,y, main = poplabel[s1],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t2)
        legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))

    }
     dev.off();
    
    
    tiff(paste0(outputDir, "Biaxial T-SNE plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)

    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
      
    }
    dev.off()
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")

}
#### END t-SNE  plot biaxial tsne plots for all files with no tissue #####
###### BEGIN Major  + Transition #####
ccast_tsne_plot_transition1 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir,"/Biaxial T-SNE plots by subpopulations1.tiff"), height = 12, width = 12,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    cols[nn]="#FAA9A5"
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D T-SNE",
    xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations2.tiff"), height = 12, width = 20,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,4))
    ss=c(1,5,6,2,7,3,4)
    for ( s in 1:(nn)) {
        s1=ss[s]
        pop1=ifelse(pop==s1,s1,0)
        nn1=length(unique(pop))
        cols<-brewer.pal(n=nn1,name="Set1")
        ##cols[nn-1]="#F9E429"
        cols[nn]="#FAA9A5"
        cols_t1<-rep(cols[s1],length(pop))
        cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
        plot(x,y, main = poplabel[s1],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t2)
        legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))
        
    }
    dev.off();
    
    
    tiff(paste0(outputDir, "Biaxial T-SNE plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    
}
###### END Major  + Transition #####

######## BEGIN ccast_pca_plot_transition1 ###############
ccast_pca_plot_transition1 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir, "/Biaxial PCA plots by subpopulations1.tiff"), height = 12, width = 12,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    cols[nn]="#FAA9A5"
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D PCA",
    xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations2.tiff"), height = 12, width = 20,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,4))
    ss=c(1,5,6,2,7,3,4)
    for ( s in 1:(nn)) {
        s1=ss[s]
        pop1=ifelse(pop==s1,s1,0)
        nn1=length(unique(pop))
        cols<-brewer.pal(n=nn1,name="Set1")
        ##cols[nn-1]="#F9E429"
        cols[nn]="#FAA9A5"
        cols_t1<-rep(cols[s1],length(pop))
        cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
        plot(x,y, main = poplabel[s1],
        xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t2)
        legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))
        
    }
    dev.off();
    
    
    tiff(paste0(outputDir, "/Biaxial PCA plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    
}
##### END ccast_pca_plot_transition1 #####
###### BEGIN TRANSITiON 2 ######################
ccast_tsne_plot_transition2 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations1.tiff"), height = 15, width = 15,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    cols[c(10,11,12)]=c("#F9E429","#FAA9A5","#F9E550")
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D T-SNE",
    xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations2.tiff"), height = 12, width = 20,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,6))
    ss=c(1,8,9,2,10,3,4,11,5,12,6,7)
    for ( s in 1:(nn)) {
        s1=ss[s]
        pop1=ifelse(pop==s1,s1,0)
        nn1=length(unique(pop))
        cols<-brewer.pal(n=nn1,name="Set1")
        ##cols[nn-1]="#F9E429"
        cols[c(10,11,12)]=c("#F9E429","#FAA9A5","#F9E550")
        cols_t1<-rep(cols[s1],length(pop))
        cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20,20,20,20,20)[factor(pop)]
        plot(x,y, main = poplabel[s1],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t2)
        if (s ==1){
        legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
        }
        
    }
    dev.off();
    
    
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20,20,20,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    
}
###### END Transition2
##### BEGIN ccast_pca_plot_transition2 #####################
ccast_pca_plot_transition2 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir, "/Biaxial PCA plots by subpopulations1.tiff"), height = 15, width = 15,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    cols[c(10,11,12)]=c("#F9E429","#FAA9A5","#F9E550")
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D PCA",
    xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff(paste0(outputDir, "/Biaxial PCA plots by subpopulations2.tiff"), height = 12, width = 20,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,6))
    ss=c(1,8,9,2,10,3,4,11,5,12,6,7)
    for ( s in 1:(nn)) {
        s1=ss[s]
        pop1=ifelse(pop==s1,s1,0)
        nn1=length(unique(pop))
        cols<-brewer.pal(n=nn1,name="Set1")
        ##cols[nn-1]="#F9E429"
        cols[c(10,11,12)]=c("#F9E429","#FAA9A5","#F9E550")
        cols_t1<-rep(cols[s1],length(pop))
        cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20,20,20,20,20)[factor(pop)]
        plot(x,y, main = poplabel[s1],
        xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t2)
        if(s==1){
        legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
        }
    }
    dev.off();
    
    
    tiff(paste0(outputDir,"/Biaxial PCA plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20,20,20,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    
}
######### END ccast_pca_plot_transition2 ########

##### BEGIN ccast_tsne_plot_transition3 #####################
ccast_tsne_plot_transition3 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir, "/Biaxial t-SNE plots by subpopulations1.tiff"), height = 15, width = 15,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    cols[6]="#FAA9A5"
    cols[c(10,11,12,13,14,15,16,17,18,19)]=c("#F8C429","#8AAADE","#C1B842","#FBC6A5","#9BBCDE","#000000","#CCFFD9","#B3FFFF","#FF5959","#FF9366")
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-rep(0.4,length(cols))[factor(pop)]
    pch_site<-rep(20,length(cols))[factor(pop)]
    plot(x,y, main ="2D T-SNE",
    xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=rep(20,length(poplabel)),bty="n",ncol=1,pt.cex=rep(0.8,length(poplabel)))
    dev.off();
    
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations2.tiff"), height = 12, width = 16,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,2))
    ss=list()
    ss[[1]]=c(1,8,9,2)
    ss[[2]]=c(2,10,3,4)
    ss[[3]]=c(4,11,12,5)
    ss[[4]]=c(5,13,6,7)
    subsetid=list()
    ###ss=c(1,8,9,2,10,3,4,11,5,12,6,7)
    for ( s in 1:(length(ss))) {
        
        pop1=ifelse(pop%in%c(ss[[s]]),pop[which(pop%in%c(ss[[s]]))],0)
        nn1=length(unique(pop))
        cols<-brewer.pal(n=nn1,name="Set1")
        cols<-brewer.pal(n=nn,name="Set1")
        cols[6]="#FAA9A5"
        cols[c(10,11,12,13,14,15,16,17,18,19)]=c("#F8C429","#8AAADE","#C1B842","#FBC6A5","#9BBCDE","#000000","#CCFFD9","#B3FFFF","#FF5959","#FF9366")

        cols_t1<-cols[pop]
        cols_t2=ifelse(pop1==0,"white",cols_t1)
        subsetid[[s]]=pop1
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-rep(0.4,length(cols))[factor(pop)]
        pch_site<-rep(20,length(cols))[factor(pop)]
        plot(x,y, main = paste(poplabel[ss[[s]]],collapse=""),
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.2, cex=cex_site, pch=pch_site,col=cols_t2)
        legend("topright",legend=poplabel[ss[[s]]],col=cols[ss[[s]]],pch=rep(20,length(poplabel))[1:length(ss[[s]])],bty="n",ncol=1,pt.cex=rep(0.8,length(poplabel))[1:length(ss[[s]])])
        
        
    }
    dev.off();
    
    
    
    for ( s in 1:(length(ss))) {
        tiff(paste(outputDir,"Biaxial T-SNE plots by subpopulations",s,"%03d.tiff",sep=""), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
        for ( r in 1:length(antibody)) {
            nrcolors = 100
            half = 1 + nrcolors/2
            colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
            rbPal = colorRampPalette(colpal)(nrcolors)
            z<-dat[,antibody[r]]
            z1<-(z - min(z))/diff(range(z))*100 + 1
            zcolor <- rbPal[z1]
            z1color2=ifelse(subsetid[[s]]==0,"white",zcolor)
            z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
            
            #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
            #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
            #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
            cex_site<-rep(0.4,length(cols))[factor(pop)]
            pch_site<-rep(20,length(cols))[factor(pop)]
            p=plot(x,y, main =antibody[r],
            xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
            cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=z1color2)
            ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
            legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        }
        dev.off()
    }
    
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    
}
######### END ccast_pca_plot_transition3 ########


##### BEGIN ccast_pca_plot_transition3 #####################
ccast_pca_plot_transition3 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir, "/Biaxial PCA plots by subpopulations1.tiff"), height = 15, width = 15,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    cols[c(10,11,12)]=c("#F8C429","#8AAADE","#C1B842")
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D PCA",
    xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff(paste0(outputDir, "/Biaxial PCA plots by subpopulations2.tiff"), height = 12, width = 16,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,2))
    ss=list()
    ss[[1]]=c(1,8,9,2)
    ss[[2]]=c(2,10,3,4)
    ss[[3]]=c(4,11,5)
    ss[[4]]=c(5,12,6,7)
    subsetid=list()
    ###ss=c(1,8,9,2,10,3,4,11,5,12,6,7)
    for ( s in 1:(length(ss))) {
        
        pop1=ifelse(pop%in%c(ss[[s]]),pop[which(pop%in%c(ss[[s]]))],0)
        nn1=length(unique(pop))
        cols<-brewer.pal(n=nn1,name="Set1")
        ##cols[nn-1]="#F9E429"
        cols[c(10,11,12)]=c("#F8C429","#8AAADE","#C1B842")
        cols_t1<-cols[pop]
        cols_t2=ifelse(pop1==0,"white",cols_t1)
        subsetid[[s]]=pop1
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20,20,20,20,20)[factor(pop)]
        plot(x,y, main = paste(poplabel[ss[[s]]],collapse=""),
        xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.2, cex=cex_site, pch=pch_site,col=cols_t2)
        legend("topright",legend=poplabel[ss[[s]]],col=cols[ss[[s]]],pch=c(20,20,20,20,20,20,20,20,20,20,20,20)[1:length(ss[[s]])],bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8)[1:length(ss[[s]])])
        
        
    }
    dev.off();
    
    
    
    for ( s in 1:(length(ss))) {
        tiff(paste(outputDir, "Biaxial PCA plots by subpopulations",s,"%03d.tiff",sep=""), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
        for ( r in 1:length(antibody)) {
            nrcolors = 100
            half = 1 + nrcolors/2
            colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
            rbPal = colorRampPalette(colpal)(nrcolors)
            z<-dat[,antibody[r]]
            z1<-(z - min(z))/diff(range(z))*100 + 1
            zcolor <- rbPal[z1]
            z1color2=ifelse(subsetid[[s]]==0,"white",zcolor)
            z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
            
            #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
            #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
            #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
            cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
            pch_site<-c(20,20,20,20,20,20,20,20,20,20,20,20)[factor(pop)]
            p=plot(x,y, main =antibody[r],
            xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
            cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=z1color2)
            ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
            legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        }
         dev.off()
    }
   
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    
}
######### END ccast_pca_plot_transition3 ########



######## Major subpopulations t-SNE plots ####
ccast_pca_plot_major1 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations1.tiff"), height = 12, width = 12,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    #cols[nn]="#FAA9A5"
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-c(0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D PCA",
    xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff("Biaxial PCA plots by subpopulations2.tiff", height = 12, width = 20,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,2))
    ss=c(1,2,3,4)
    for ( s in 1:(nn)) {
        s1=ss[s]
        pop1=ifelse(pop==s1,s1,0)
        nn1=length(unique(pop))
        cols<-brewer.pal(n=nn1,name="Set1")
        ##cols[nn-1]="#F9E429"
        #cols[nn]="#FAA9A5"
        cols_t1<-rep(cols[s1],length(pop))
        cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20)[factor(pop)]
        plot(x,y, main = poplabel[s1],
        xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t2)
        legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8))
        
    }
    dev.off();
    
    
    tiff(paste0(outputDir, "/Biaxial PCA plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    
}
##### END ccast_pca_plot_major1
##############################################################
ccast_tsne_plot_major1 <-function(x,y,pop,poplabel,antibody,dat, outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir,"Biaxial T-SNE plots by subpopulations1.tiff"), height = 12, width = 12,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    #cols[nn]="#FAA9A5"
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-c(0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D T-SNE",
    xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff("Biaxial T-SNE plots by subpopulations2.tiff", height = 12, width = 20,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,2))
    ss=c(1,2,3,4)
    for ( s in 1:(nn)) {
        s1=ss[s]
        pop1=ifelse(pop==s1,s1,0)
        nn1=length(unique(pop))
        cols<-brewer.pal(n=nn1,name="Set1")
        ##cols[nn-1]="#F9E429"
        #cols[nn]="#FAA9A5"
        cols_t1<-rep(cols[s1],length(pop))
        cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20)[factor(pop)]
        plot(x,y, main = poplabel[s1],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t2)
        legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8))
        
    }
    dev.off();
    
    
    tiff(plot0(outputDir, "/Biaxial T-SNE plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    
}
ccast_tsne_plot_major2 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations1.tiff"), height = 12, width = 12,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    cols[nn]="#FAA9A5"
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D T-SNE",
    xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations2.tiff"), height = 12, width = 20,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,4))
    
    ss=c(1,2,3,4,4,5,6,7)

    for ( s in 1:length(ss)) {
        s1=ss[s]
        pop1=ifelse(pop==s1,s1,0)
        nn1=length(unique(pop))
        cols<-brewer.pal(n=nn1,name="Set1")
        ##cols[nn-1]="#F9E429"
        cols[nn1]="#FAA9A5"
        cols_t1<-rep(cols[s1],length(pop))
        cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
        plot(x,y, main = poplabel[s1],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t2)
        legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))
        
    }
    dev.off();
    
    
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    
}

######## END major subpopulations t-SNE plots ####

######## BEGIN major subpopulations t-SNE with tissue plots ####
ccast_tsne_plot_major22 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir, "Biaxial T-SNE plots by subpopulations1.tiff"), height = 12, width = 12,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    cols[nn-2]="#FAA9A5"
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D T-SNE",
    xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations2.tiff"), height = 12, width = 20,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,4))
    
    ss=c(1,2,3,4,5,6,7,8)
    nn1=length(unique(pop))
    cols<-brewer.pal(n=nn1,name="Set1")
    cols[nn1-2]="#FAA9A5"
    for ( s in 1:length(ss)) {
        s1=ss[s]
        pop1=ifelse(pop==s1,s1,0)
        cols_t1<-rep(cols[s1],length(pop))
        cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20)[factor(pop)]
        plot(x,y, main = poplabel[s1],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t2)
        legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
        
    }
    dev.off();
    
    
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    
}
##### END of TSNE major with tissue
######## BEGIN ccast_tsne_plot_major222 with tissue ####

ccast_tsne_plot_major222 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations1.tiff"), height = 12, width = 12,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    cols[nn-3]="#FAA9A5"
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D T-SNE",
    xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff(paste0(outputDir, "Biaxial T-SNE plots by subpopulations2.tiff"), height = 15, width = 15,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    #par(mfrow=c(2,4))
    
    ss=c(1,2,3,4,5,6,7,8,9)
    nn1=length(unique(pop))
    cols<-brewer.pal(n=nn1,name="Set1")
    cols[nn1-3]="#FAA9A5"
    
    
    for ( s in 1:length(ss)) {
        
        if (s %in% c(8,9)) {
            s1=ss[s]
            # pop1=ifelse(pop==s1,s1,0)
            
            ##cols[nn-1]="#F9E429"
            
            #cols_t1<-rep(cols[s1],length(pop))
            #cols_t2=ifelse(pop1==s1,cols_t1,"white")
            ##xydens <- kde2d(x,y,n=100)
            #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
            #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
            #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
            cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
            pch_site<-c(20,20,20,20,20,20,20,20,20)[factor(pop)]
            #plot(x[which(pop==s1)],y[which(pop==s1)], main = "Cluster boundaries",
            #xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
            #cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols[s1])
            points(x[which(pop==s1)],y[which(pop==s1)], main = "Cluster boundaries",
            xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
            cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols[s1])
            #X=cbind(x[which(pop==s1)],y[which(pop==s1)])
            #hpts <- chull(X)
            #hpts <- c(hpts, hpts[1])
            #lines(X[hpts, ],lwd=2, col=cols[s1])
            Plot_ConvexHull(xcoord = x[which(pop==s1)], ycoord = y[which(pop==s1)], lcolor=cols[s1],linew=2)
        } else {
            s1=ss[s]
            # pop1=ifelse(pop==s1,s1,0)
            ##cols[nn-1]="#F9E429"
            #cols_t1<-rep(cols[s1],length(pop))
            #cols_t2=ifelse(pop1==s1,cols_t1,"white")
            ##xydens <- kde2d(x,y,n=100)
            #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
            #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
            #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
            cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
            pch_site<-c(20,20,20,20,20,20,20,20,20)[factor(pop)]
            if (s==1) {
                plot(x[which(pop==s1)],y[which(pop==s1)], main = "Cluster boundaries",
                xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
                cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols[s1])
                Plot_ConvexHull(xcoord = x[which(pop==s1)], ycoord = y[which(pop==s1)], lcolor=cols[s1],linew=2)
            }
            points(x[which(pop==s1)],y[which(pop==s1)], main = "Cluster boundaries",
            xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
            cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col="white")
            #X=cbind(x[which(pop==s1)],y[which(pop==s1)])
            #hpts <- chull(X)
            #hpts <- c(hpts, hpts[1])
            #lines(X[hpts, ],lwd=2, col=cols[s1])
            Plot_ConvexHull(xcoord = x[which(pop==s1)], ycoord = y[which(pop==s1)], lcolor=cols[s1],linew=2)
        }
        
    }
    
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    
    tiff(paste0(outputDir, "/Biaxial T-SNE plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
    
}

######## END ccast_tsne_plot_major222 with tissue ####

########################

ccast_pca_plot_major2 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir, "/Biaxial PCA plots by subpopulations1.tiff"), height = 12, width = 12,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    cols[nn-2]="#FAA9A5"
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D PCA",
    xlab = "PCA 1", ylab = "PCA 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff(paste0(outputDir,"Biaxial PCA plots by subpopulations2.tiff"), height = 12, width = 20,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,4))
    
    ss=c(1,2,3,4,4,5,6,7)
    
    for ( s in 1:length(ss)) {
        s1=ss[s]
        pop1=ifelse(pop==s1,s1,0)
        nn1=length(unique(pop))
        cols<-brewer.pal(n=nn1,name="Set1")
        ##cols[nn-1]="#F9E429"
        cols[nn1-2]="#FAA9A5"
        cols_t1<-rep(cols[s1],length(pop))
        cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
        plot(x,y, main = poplabel[s1],
        xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t2)
        legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))
        
    }
    dev.off();
    
    
    tiff(paste0(outputDir,"/Biaxial PCA plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]

        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
}
######## END ccast_pca_plot_major2####
######## BEGIN ccast_pca_plot_major2 with tissue ####

ccast_pca_plot_major22 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir,"Biaxial PCA plots by subpopulations1.tiff"), height = 12, width = 12,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    cols[nn-2]="#FAA9A5"
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
    cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D PCA",
    xlab = "PCA 1", ylab = "PCA 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff(paste0(outputDir,"/Biaxial PCA plots by subpopulations2.tiff"), height = 12, width = 20,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,4))
    
    ss=c(1,2,3,4,5,6,7,8)
    nn1=length(unique(pop))
    cols<-brewer.pal(n=nn1,name="Set1")
    cols[nn1-2]="#FAA9A5"
    for ( s in 1:length(ss)) {
        s1=ss[s]
        pop1=ifelse(pop==s1,s1,0)
        
        ##cols[nn-1]="#F9E429"
        
        cols_t1<-rep(cols[s1],length(pop))
        cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20)[factor(pop)]
        plot(x,y, main = poplabel[s1],
        xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t2)
        legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
        
    }
    dev.off();
    
    
    tiff(paste0(outputDir,"/Biaxial PCA plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
    ###legend("topright",legend=paste(rep(c("Small","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",ncol=2,cex=0.7,pt.cex=0.7)
    ####abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
    #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
    ##contour(xydens,nlevels=10,add=TRUE,col="dark orange")
    
}
######## END ccast_pca_plot_major2 with tissue ####
######## BEGIN ccast_pca_plot_major222 with tissue ####
Plot_ConvexHull<-function(xcoord, ycoord, lcolor,linew){
    hpts <- chull(x = xcoord, y = ycoord)
    hpts <- c(hpts, hpts[1])
    lines(xcoord[hpts], ycoord[hpts], col = lcolor,lwd=linew)
}

ccast_pca_plot_major222 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir, "/Biaxial PCA plots by subpopulations1.tiff"), height = 12, width = 12,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    ###cols[nn-1]="#F9E429"
    cols[nn-3]="#FAA9A5"
    cols_t1<-cols[pop]
    cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D PCA",
    xlab = "PCA 1", ylab = "PCA 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff(paste0(outputDir, "Biaxial PCA plots by subpopulations2.tiff"), height = 15, width = 15,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    #par(mfrow=c(2,4))
    
    ss=c(1,2,3,4,5,6,7,8,9)
    nn1=length(unique(pop))
    cols<-brewer.pal(n=nn1,name="Set1")
    cols[nn1-3]="#FAA9A5"
    
    
    for ( s in 1:length(ss)) {
        
        if (s %in% c(8,9)) {
            s1=ss[s]
            # pop1=ifelse(pop==s1,s1,0)
            
            ##cols[nn-1]="#F9E429"
            
            #cols_t1<-rep(cols[s1],length(pop))
            #cols_t2=ifelse(pop1==s1,cols_t1,"white")
            ##xydens <- kde2d(x,y,n=100)
            #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
            #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
            #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
            cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
            pch_site<-c(20,20,20,20,20,20,20,20,20)[factor(pop)]
            #plot(x[which(pop==s1)],y[which(pop==s1)], main = "Cluster boundaries",
            #xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
            #cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols[s1])
            points(x[which(pop==s1)],y[which(pop==s1)], main = "Cluster boundaries",
            xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
            cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols[s1])
            #X=cbind(x[which(pop==s1)],y[which(pop==s1)])
            #hpts <- chull(X)
            #hpts <- c(hpts, hpts[1])
            #lines(X[hpts, ],lwd=2, col=cols[s1])
            Plot_ConvexHull(xcoord = x[which(pop==s1)], ycoord = y[which(pop==s1)], lcolor=cols[s1],linew=2)
        } else {
        s1=ss[s]
        # pop1=ifelse(pop==s1,s1,0)
        ##cols[nn-1]="#F9E429"
        #cols_t1<-rep(cols[s1],length(pop))
        #cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20,20)[factor(pop)]
        if (s==1) {
        plot(x[which(pop==s1)],y[which(pop==s1)], main = "Cluster boundaries",
         xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols[s1])
        Plot_ConvexHull(xcoord = x[which(pop==s1)], ycoord = y[which(pop==s1)], lcolor=cols[s1],linew=2)
        }
            points(x[which(pop==s1)],y[which(pop==s1)], main = "Cluster boundaries",
            xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
            cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col="white")
            #X=cbind(x[which(pop==s1)],y[which(pop==s1)])
            #hpts <- chull(X)
            #hpts <- c(hpts, hpts[1])
            #lines(X[hpts, ],lwd=2, col=cols[s1])
            Plot_ConvexHull(xcoord = x[which(pop==s1)], ycoord = y[which(pop==s1)], lcolor=cols[s1],linew=2)
        }
        
        }
    
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    
    tiff(paste0(outputDir, "Biaxial PCA plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        #cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2)[factor(pop)]
        cex_site<-c(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "PC 1", ylab = "PC 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
    
}
######## END ccast_pca_plot_major222 with tissue ####
##### t-SNE PLOT with TISSUE #######
ccast_tsne_plot33 <-function(x,y,pop,poplabel,antibody,dat,outputDir) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ###pdf(file="Biaxial T-SNE plots by subpopulations.pdf",width=10, height=10)
    tiff(paste0(outputDir,"Biaxial T-SNE plots by subpopulations1.tiff"), height = 12, width = 12,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    #cex <- 1.8
    require(RColorBrewer)
    ###display.brewer.all()
    nn=length(unique(pop))
    #print(nn)
    cols<-brewer.pal(n=nn,name="Set1")
    #cols[nn-2]="#F9E429"
    cols[nn-2]="#FAA9A5"
    cols_t1<-cols[pop]
    ##xydens <- kde2d(x,y,n=100)
    #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
    #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
    cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2,0.2)[factor(pop)]
    pch_site<-c(20,20,20,20,20,20,20,20)[factor(pop)]
    plot(x,y, main ="2D T-SNE",
    xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t1)
    legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
    dev.off();
    
    tiff(paste0(outputDir,"Biaxial T-SNE plots by subpopulations2.tiff"), height = 12, width = 20,units = 'cm',pointsize = 5,compression = "lzw", type="cairo", res=300)
    par(mfrow=c(2,4))
    ss=c(1,3,4,2,5,6,7,8)
    for ( s in 1:nn) {
        s1=ss[s]
        pop1=ifelse(pop==s1,s1,0)
        nn1=length(unique(pop))
        cols<-brewer.pal(n=nn1,name="Set1")
        ##cols[nn-2]="#F9E429"
        cols[nn-2]="#FAA9A5"
        cols_t1<-rep(cols[s1],length(pop))
        cols_t2=ifelse(pop1==s1,cols_t1,"white")
        ##xydens <- kde2d(x,y,n=100)
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2,0.2)[factor(pop)]
        pch_site<-c(20,20,20,20,20,20,20,20)[factor(pop)]
        plot(x,y, main = poplabel[s1],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=cols_t2)
        legend("topright",legend=poplabel,col=cols,pch=c(20,20,20,20,20,20,20,20),bty="n",ncol=1,pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8))
        
    }
    dev.off();
    
    
    tiff(paste0(outputDir, "Biaxial T-SNE plots by subpopulations3%03d.tiff"), height = 15, width = 15,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    
    for ( r in 1:length(antibody)) {
        nrcolors = 100
        half = 1 + nrcolors/2
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
        rbPal = colorRampPalette(colpal)(nrcolors)
        z<-dat[,antibody[r]]
        z1<-(z - min(z))/diff(range(z))*100 + 1
        zcolor <- rbPal[z1]
        z1color<-rbPal[quantile(z1,  probs = c(1,c(1:10)*9.5,100)/100)]
        
        #cex_site<-c(0.2,0.2,0.6,0.4,0.6,0.4,0.2)[factor(pop)]
        #pch_site<-c(20,20,19,18,19,18,20)[factor(pop)]
        cex_site<-c(0.2,0.2,0.2,0.4,0.2,0.4,0.2,0.2)[factor(pop)]
        pch_site<-c(20,20,19,18,19,18,20,20)[factor(pop)]
        p=plot(x,y, main =antibody[r],
        xlab = "t-SNE 1", ylab = "t-SNE 2", ylim = c(min(y)-2, max(y)+5),xlim = c(min(x)-2, max(x)+5),
        cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8, cex=cex_site, pch=pch_site,col=zcolor)
        ##legend("topright",legend=poplabel,col=cols,pch=c(20,20,19,18),bty="n",ncol=1,pt.cex=c(0.2,0.2,0.8,0.6))
        legend("topright",title="Percentile",legend=c(1,c(1:10)*10),col =z1color,pch=20)
        
    }
    dev.off()
    
}


ccast_tsne_plot2b <-function(x,y,test,file1="Projection T-SNE plots2.tiff",sample="Tissue 2",window=NULL) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ##pdf(file="Biaxial T-SNE plots.pdf",width=10, height=10)
    if (is.null(window)){
        window = c(min(x)-2, max(x)+5, min(y)-2, max(y)+5)
    }
    tiff(file1, height = 10, width = 10,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    cex <- 1.8
    xydens <- kde2d(x,y,n=100)
    plot(x,y, main = sample,
    xlab = "t-SNE1", ylab = "t-SNE2", ylim = window[3:4],xlim = window[1:2],
    cex.axis = cex, cex.lab = cex, cex.main = cex,cex=0.2, pch=20,col="grey")

    require(plot3D)
   
    contour2D(z=xydens$z,x=xydens$x,y=xydens$y, add=TRUE, lwd = 2, colkey = TRUE)
    P <- voronoi.polygons(test)
    plot(test,add=TRUE,lwd=2)
    dev.off();
    
}


##### KEEP
ccast_biaxialplot2 <- function(optccastreeoutput,DD,ylabel,outputDir) {
    tiff(paste0(outputDir,"/Biaxial_tree_node_plots%03d.tiff"), height = 10, width = 10,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    cex <- 1.6
    tree33<-optccastreeoutput[[3]]
    dD33 <- DD
    inner <- party:::nodes(tree33, c(1:max(unique(where(tree33))))[-unique(where(tree33))])
    out <- sapply(inner, function(i) {
       ###splitstat <- i$psplit$splitstatistic
        w1=which(i$weights==1)
        x <- dD33[[i$psplit$variableName]][w1]
        y<-dD33[[ylabel]][w1]
        if (bandwidth.nrd(x)>0&bandwidth.nrd(y)>0) {
            xydens <- kde2d(x,y,n=25,lims = c(range(x), range(y)))
        }else {
            print("Bandwidths=0")
            print(c(bandwidth.nrd(x), bandwidth.nrd(y)))
            xydens <- kde2d(x,y,h=0.5,n=25,lims = c(range(x), range(y)))
        }

        plot(x,y, main = paste("Node", i$nodeID,sep=" "),
        xlab = i$psplit$variableName, ylab = ylabel, ylim = c(0, max(y)+5),
        cex.axis = cex, cex.lab = cex, cex.main = cex,cex=0.2, pch=20,col="dark blue")
        abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
        contour(xydens,nlevels=10,add=T,col="dark orange")
    })
    dev.off()
}

#### Plotting bar plots ###
bar_plot_ccast <- function(M,M2,file,outputDir) {
  pdf(paste(outputDir, file,sep="/"),width=15, height=10)
  u1<-ifelse(min(M)<0,min(M)-1,min(M))
  u2<-max(M)+1
  for (i in 1:dim(M)[1]) {
      mp<-barplot(M[i,]+u1,names.arg=colnames(M),main=rownames(M)[i],ylim=c(u1-2,(u2+1)),yaxt="n",yaxs="r",family= "sans",space=c(1,1),xlab="",ylab="",cex.lab=2)
      
    axis(1,mp, labels = FALSE, tick = T)
    axis(2,M[i,]+u1,labels = FALSE, tick = T)
    arrows(mp,M[i,]+u1,mp,(M[i,]+u1)+ M2[i,],length=0.1, angle = 90)
    arrows(mp,M[i,]+u1,mp,(M[i,]+u1)- M2[i,], length=0.1,angle = 90)
    mtext(min(M[i,]),side = 2, at =min(M[i,])+u1 ,family= "sans", line = 1)
    
  }
  dev.off()
}

#### Plotting heatmaps ###

ccast_heatmaplot <-function(optccastreeoutput,outputDir) {
    
    ##pdf(file="Heatmaps Homogeneous cells.pdf",width=10, height=10)
    tiff(paste0(outputDir, "/Heatmaps_Homogeneous_cells%03d.tiff"), height = 10, width = 10,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
    tree3<-optccastreeoutput$tree
    tt1=party:::nodes(tree3, unique(where(tree3)))
    dD3=optccastreeoutput$puredata[[length(optccastreeoutput$puredata)]]
    nrcolors = 100
    half = 1 + nrcolors/2
    colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
    colorpalette = colorRampPalette(colpal)(nrcolors)
    ##layout(matrix(1:length(tt1), ncol = round(sqrt(length(tt1)),0)))
    MM <- NULL
    MM2 <- NULL
    rname<-NULL
    rowids<-list()
    ###par(mfrow=c(4,2)))
    homodata=list()
    for( j in 1:length(tt1)) {
      mx=which.max(tt1[[j]]$prediction)
      w1=which(tt1[[j]]$weights==1)
      w2=w1[which(dD3$groups[w1]==mx)]
      print(table(dD3$groups[w2]))
      minv=min(dD3[,-dim(dD3)[2]],na.rm=TRUE)
      maxv=max(dD3[,-dim(dD3)[2]],na.rm=TRUE)
      key_color1=seq(minv,maxv,length.out=dim(dD3[,-dim(dD3)[2]])[2])
      rowids[[j]]=as.numeric(rownames(dD3[w2,-dim(dD3)[2]]))
      BB=as.matrix(rbind(dD3[w2,-dim(dD3)[2]],key_color1))
      rownames(BB)=rep(paste("Cell",j),dim(BB)[1])
      mm <- round(c(apply(BB[-dim(BB)[1],], 2, mean)),3)
      MM <- rbind(MM,mm)
      mm2 <- round(c(apply(BB[-dim(BB)[1],], 2, sd)),3)
      MM2 <- rbind(MM2,mm2)
      image(x=1:dim(BB)[2],y=1:dim(BB)[1],z=as.matrix(t(BB)),col = colorpalette,xlab = paste("Celltype",mx),xaxt="n",yaxt="n",ylab = "Cells",family="sans",cex.lab=1.4,main=paste("Node",tt1[[j]]$ nodeID))
      mtext(rownames(BB),at=c(1:dim(BB)[1]),side=4,las=2,line=1, family= "sans", cex=0.4)
      mtext(colnames(BB),at=c(1:dim(BB)[2]),side=1,las=2,line=1, family= "sans", cex=0.4)
      homodata[[j]]=BB[-dim(BB)[1],]
      rownames(homodata[[j]])=paste(rowids[[j]])
      rname=c(rname,paste("Node",tt1[[j]]$ nodeID))
    }
    
    rownames(MM)=rname
    rownames(MM2)=rname
    save(homodata,file=paste0(outputDir, "/combinedatalist.rdata"))
    save(rname,file=paste0(outputDir, "/NodeIDs.rdata"))
    save(rowids,file=paste0(outputDir, "/rowids.rdata"))
    
    dev.off();
    
    return(list(ccastmeans=MM,ccastsds=MM2))
    
  }

ccast_heatmaplot2 <-function(optccastreeoutput, DD, s.dat,outputDir) {
    tiff(paste0(outputDir, "/Heatmaps_Homogeneous_cells%03d.tiff"), height = 10, width = 10,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
    tree3 <- optccastreeoutput[[3]]
    tt1 <- optccastreeoutput[[1]]
    print("TT1")
    print(tt1)
    
    #return(tt1)
    
    dD3 <- DD
    
    nrcolors = 100
    half = 1 + nrcolors/2
    colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
    colorpalette = colorRampPalette(colpal)(nrcolors)
    
    MM <- MM2 <- rname <- NULL
    rowids <- homodata <- list()
    nd <- unique(tt1)
    for( j in 1:length(nd)) {
      ##w1 <- which(tt1 == nd[j])
      #w1 <- nd[j]
      w1 <- j
      minv=min(dD3[w1,-c((dim(dD3)[2]-1),dim(dD3)[2])],na.rm=TRUE)
      maxv=max(dD3[w1,-c((dim(dD3)[2]-1),dim(dD3)[2])],na.rm=TRUE)
      key_color1=seq(minv,maxv,length.out=dim(dD3[,-c((dim(dD3)[2]-1),dim(dD3)[2])])[2])
      
      rowids[[j]] <- as.numeric(rownames(dD3[w1,-c((dim(dD3)[2]-1),dim(dD3)[2])]))
      BB <- as.matrix(rbind(dD3[w1,-c((dim(dD3)[2]-1),dim(dD3)[2])],key_color1))
      rownames(BB) <- rep(paste("Cell",j),dim(BB)[1])
      print("BB")
      print(dim(BB))
      print(BB)
      mm <- round(c(apply(BB[-dim(BB)[1],], 2, mean)),3)
      MM <- rbind(MM,mm)
      mm2 <- round(c(apply(BB[-dim(BB)[1],], 2, sd)),3)
      MM2 <- rbind(MM2,mm2)
      image(x=1:dim(BB)[2],y=1:dim(BB)[1],z=as.matrix(t(BB)),col = colorpalette,xlab = paste("Subpopulationnode",nd[j]),xaxt="n",yaxt="n",ylab = "Cells",family="sans",cex.lab=1.4,main=paste("Node",nd[j]))
      mtext(rownames(BB),at=c(1:dim(BB)[1]),side=4,las=2,line=1, family= "sans", cex=0.4)
      mtext(colnames(BB),at=c(1:dim(BB)[2]),side=1,las=2,line=1, family= "sans", cex=0.4)
      homodata[[j]]=BB[-dim(BB)[1],]
      rownames(homodata[[j]])=paste(rowids[[j]])
      rname=c(rname,paste("Node",nd[j]))
    }
    rownames(MM)=rname
    rownames(MM2)=rname
    if(s.dat == TRUE) {
      save(homodata,file=paste0(outputDir, "combinedatalist.rdata"))
      save(rname,file=paste0(outputDir, "/NodeIDs.rdata"))
      save(rowids,file=paste0(outputDir, "/rowids.rdata"))
    }
    dev.off()
   return(list(ccastmeans=MM,ccastsds=MM2))
}
