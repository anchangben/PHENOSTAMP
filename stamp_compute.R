## This file has implemented the methods to run a ccast tree clustering
source("./stamp_visualise.R")

load_dataset <- function(filename, patient, colid=NULL, coln=NULL, rown=NULL, transformlogic=FALSE, asinhp=1,
                         subanalysis=TRUE, subsamplesize=NULL, colid_proteins=NULL, coln_proteins=NULL){
    # This function subsaples the files present in filename with extension fcs. It uses fcsimport3 function for the downsampling and aggregation.

    if(is.matrix(file)) { global_dataset <- filename } 
    else {
      if (!(grepl('fcs', filename))) setwd(filename)
      Data_all <- fcsimport3(file=filename,transformlogic=transformlogic,asinhp=asinhp,
                          colid=colid,coln=coln,rown=rown,patient=patient,subanalysis=subanalysis,
                          subsamplesize=subsamplesize)
      subsampled_dataset <- Data_all[[1]]
      #Dall is the subsampled columns and lines
      fcsfile <- Data_all[[3]]
      #fcsfile is a descriptive file
      global_dataset <- Data_all[[2]]
      #all data is the subsampled lines
      rownames(global_dataset) <- paste(1:dim(global_dataset)[1])
    }
    subsampled_dataset_processed <- asinh(global_dataset[,colid]/asinhp)
    colnames(subsampled_dataset_processed) <- coln
    rownames(subsampled_dataset_processed) <- paste(1:dim(subsampled_dataset_processed)[1])
    subsampled_dataset <- subsampled_dataset_processed

    global_proteins_dataset <- global_dataset[, colid_proteins]
    colnames(global_proteins_dataset) = coln_proteins
    rownames(global_proteins_dataset) <- paste(1:dim(global_proteins_dataset)[1])

    return(list(subsampled_dataset = subsampled_dataset, global_dataset = global_dataset, fcsfile=fcsfile, global_proteins_dataset = global_proteins_dataset))

}

ccast_tree_onepass <- function(dataset,groups, outputDir, max_depth = 15, filter=TRUE, s.plot=TRUE, figname='initial_tree.pdf', min_bucket=100){
    # This function creates an "unpruned" CCAST tree form a dataset and labels.
    # With filter, make tree deeper until each groups is majoritary in at least one node.
    # without , just fit CCAST using max_depth
    # Apply a filter on the size of each leaf with min_bucket
    # Plot the resulting tree using s.plot=True

  li_groups_maj <- 1
  ## list of the groups that are majoritary in at least a leaf of the tree
  num_groups = length(groups)
  if(filter){
    for (s in 1:max_depth) {
        if(length(li_groups_maj) < num_groups) {
        ##creates deeper trees until each prediction is majoritary in at least a group
        tree <- party::ctree(formula=as.factor(groups) ~ ., data = dataset,controls = ctree_control(maxdepth = s, minbucket = min(min_bucket,dim(dataset)[1]/2)))
        tree_nodes = party:::nodes(tree, unique(where(tree)))
        num_leaves=length(tree)
        li_groups_maj=unique(unlist(lapply(tree_nodes, function(x) which.max(x$prediction))))
        }
        else {
        break
        }
    }
  }
  else{
    s=max_depth + 1
    tree <- ctree(formula=as.factor(groups) ~ ., data = dataset,controls = ctree_control(maxdepth = max_depth, minbucket = min(min_bucket,dim(dataset)[1]/2)))
    tree_nodes = party:::nodes(tree, unique(where(tree)))
  }

  if(s.plot == TRUE) {
    pdf(file=paste(outputDir, figname, sep='/'),width=15, height=10)
    plot(tree,xlab="Celltypes",main="Initial_CCAST_tree")
    dev.off();
  }
  return(list(tree_nodes=tree_nodes,depth=s-1, tree=tree))
}

ccast_prune_tree <- function(tree_nodes, dataset, tree_depth, outputDir, filename="final_tree", s.save=FALSE, s.plot=FALSE, maxiter=10, min_bucket=100){
    # Prune a "raw" CCAST tree
    # The pruning is performed by removing each cell that is not in the dominant group of each leaf, and fitting the tree again until every leaf
    # if pure.

    # dataset must include the prediction as a last column
    last=dim(dataset)[2]
    if(colnames(dataset)[last]!="groups") colnames(dataset)[last]<-"groups"
    groups = dataset[, last]
    prop_maj_clust=unique(unlist(lapply(tree_nodes, function(x) max(x$prediction))))
    # This list holds the unique percentages ratio between the majoritarian pop in each cluster and its cluster size
    # At the end of the pruning, it will be 1 for each cluster 
    iter=1
    successive_pruned_datasets=list()
    successive_pruned_datasets[[iter]]=dataset
    minleaf=1:max(dataset[,last])
    if (length(prop_maj_clust)==1) {
        
        tree <- party::ctree(formula=as.factor(groups) ~ ., 
                             data = dataset,controls = ctree_control(maxdepth = tree_depth,
                                                                     mincriterion = 0.98,
                                                                     minbucket = min(min_bucket,dim(dataset)[1]/2)))

        range_indexes=c(1:dim(dataset)[1])
        aliencells <- secondlargestalienpop <- topalieansubpops <- largestalienpop <- secondlargestsubpop <- NULL
        largestsubpop <- cleanmainsubpops <- c(1:dim(dataset)[1])

        if(s.save == TRUE) {
          save(tree,file=paste(outputDir, "finaltree.rdata", sep='/'))
        }
        if(s.plot == TRUE) {
          pdf(file=paste(outputDir,"finaltree.pdf", sep='/'),width=15, height=10)
          plot(tree,xlab="Celltypes",main="Final CCAST tree")
          dev.off();
        }
        return(list(tree=tree, puredata=list(dataset)))
    } else {
        # If every cluster is not clean, we purge until it is the case
        li_largest_wrong_clu <- li_largest_correct_clu <- li_second_largest_correct_clu <- li_second_largest_wrong_clu <- list()
        while(((sum(prop_maj_clust==1)!=length(prop_maj_clust)) & (iter<=maxiter))) {
            #Do until all clusters are pure or more than 5 iterations

            dom_clu <- vec_dom_clu_idx <- NULL
            li_not_dom_clu <- li_dom_clu_idx <- list()
            for( j in 1:length(tree_nodes)) {
                name_dom_clu=which.max(tree_nodes[[j]]$prediction)
                ##name_dom_clu is the dominant prediction in the cluster j
                idx_in_clu=which(tree_nodes[[j]]$weights==1)
                ##idx_in_clu is the list of indices that are in cluster j
                temp_dom_clu_idx=idx_in_clu[which(successive_pruned_datasets[[iter]][[last]][idx_in_clu]==name_dom_clu)]
                ##temp_dom_clu_idx is the list of indices that are in cluster j and affected to the dominant cluster
                temp_dom_clu=names(table(successive_pruned_datasets[[iter]][[last]][temp_dom_clu_idx]))
                #temp_dom_clu is the name of the dominant cluster 
                dom_clu=c(dom_clu,temp_dom_clu)

                vec_dom_clu_idx=c(vec_dom_clu_idx,temp_dom_clu_idx)
                #vec_dom_clu_idx concatenates the temp_dom_clu_idx from the different clusters
                li_not_dom_clu[[j]]=idx_in_clu[-which(successive_pruned_datasets[[iter]][[last]][idx_in_clu]==name_dom_clu)]
                ##li_not_dom_clu is the list of indices in cluster j not in the main prediction
                li_dom_clu_idx[[j]]=temp_dom_clu_idx
                ##li_dom_clu_idx stores the temp_dom_clu_idx lists
            }
            largest_wrong_clu=which.max(lapply(li_not_dom_clu,length))
            ## largest_wrong_clu will be the cluster in which the most samples are wrongly affected
            largest_wrong_clu_idx=li_not_dom_clu[[largest_wrong_clu]]
            li_largest_wrong_clu[[iter]]=largest_wrong_clu_idx
            ## largest_wrong_clu_idx will be the list of indexes badly affected in the cluster in which the most are badly affected
            second_largest_wrong_clu=which.max(lapply(li_not_dom_clu[-largest_wrong_clu],length))
            ##Second cluster where th most are rongly affected
            second_largest_wrong_clu_idx=li_not_dom_clu[-largest_wrong_clu][[second_largest_wrong_clu]]
            ##list of the second most wrongly affected 
            li_second_largest_wrong_clu[[iter]]=second_largest_wrong_clu_idx

            largest_correct_clu=which.max(lapply(li_dom_clu_idx,length))
            largest_correct_clu_idx=li_dom_clu_idx[[largest_correct_clu]]
            ##largest_correct_clu is the name of the cluster in which the most are well affected
            li_largest_correct_clu[[iter]]=largest_correct_clu_idx
            second_largest_correct_clu=which.max(lapply(li_dom_clu_idx[-largest_correct_clu],length))
            ##second_largest_correct_clu will be the name of the second cluster where the most are well affected 
            second_largest_correct_clu_idx=li_dom_clu_idx[-largest_correct_clu][[second_largest_correct_clu]]
            ##second_largest_correct_clu_idx is the list of the indexes in the second cluster where the most are well affected
            li_second_largest_correct_clu[[iter]]=second_largest_correct_clu_idx

            if(all(minleaf %in% dom_clu)) {
                ##There is at least one cluster in which each of the prediction is dominant
                iter=iter+1
                successive_pruned_datasets[[iter]]=as.data.frame(as.matrix(successive_pruned_datasets[[iter-1]])[vec_dom_clu_idx,,drop=FALSE])
                
                tree <- party::ctree(as.factor(groups) ~ ., data = successive_pruned_datasets[[iter]],controls = ctree_control(maxdepth = tree_depth,mincriterion = 0.98,minbucket = min(min_bucket,dim(successive_pruned_datasets[[iter]])[1]/2)))
                tree_nodes=party:::nodes(tree, unique(where(tree)))
                prop_maj_clust=unique(unlist(lapply(tree_nodes, function(x) max(x$prediction))))
                
            } else {
                tree <- party::ctree(as.factor(groups) ~ ., data = successive_pruned_datasets[[iter]],controls = ctree_control(maxdepth = tree_depth,mincriterion = 0.98,minbucket = min(min_bucket,dim(successive_pruned_datasets[[iter]])[1]/2)))
                tree_nodes=party:::nodes(tree, unique(where(tree)))
                prop_maj_clust=unique(unlist(lapply(tree_nodes, function(x) max(x$prediction))))
                break
            }
        }
    
        if(length(li_largest_wrong_clu)>0) {
          ## There are alien cells in the clusters
                aliencells=c(1:dim(dataset)[1])[-unique(vec_dom_clu_idx)]
                iter_most_alien=which.max(lapply(li_largest_wrong_clu,length))
                ##Iteration at which there are the most alien cells
                iter_second_most_alien=which.max(lapply(li_second_largest_wrong_clu,length))
                ##Iteration at which the second sample set that are wrongly affected is the largest
                largestalienpop=unlist(li_largest_wrong_clu[[iter_most_alien]])
                ##List of the indexes wonrly affected in the iteration ww4
                secondlargestalienpop=unlist(li_second_largest_wrong_clu[[iter_second_most_alien]])
                topalieansubpops=unique(unlist(li_not_dom_clu))
            } else {
                print("No alien cells")
                largestalienpop <- secondlargestalienpop <- topalieansubpops <- NULL
            }
            
        if(length(li_largest_correct_clu)>0) {
                iter_most_correct=which.max(lapply(li_largest_correct_clu,length))
                iter_second_most_correct=which.max(lapply(li_second_largest_correct_clu,length))
                largestsubpop=unlist(li_largest_correct_clu[[iter_most_correct]])
                secondlargestsubpop=unlist(li_second_largest_correct_clu[[iter_second_most_correct]])
                cleanmainsubpops=unique(unlist(li_dom_clu_idx))
            } else {
                print("No target cells")
                cleanmainsubpops <- secondlargestsubpop <- largestsubpop <- NULL
            }
           
        if(s.save == TRUE) {save(tree,file=paste(outputDir, paste(filename, ".rdata"), sep='/'))}
        if(s.plot == TRUE) {
          pdf(file=paste(outputDir,paste(filename, ".pdf", sep=""), sep='/'),width=15, height=10)
          plot(tree,xlab="Celltypes",main="Final CCAST tree")
          dev.off()
        }
    }
    return(list(tree=tree,fprediction=prop_maj_clust,puredata=successive_pruned_datasets,targetcells=vec_dom_clu_idx,
                aliencells=aliencells,largestalienpop=largestalienpop,secondlargestalienpop=secondlargestalienpop,
                topalieansubpops=topalieansubpops,largestsubpop=largestsubpop,cleanmainsubpops=cleanmainsubpops,
                secondlargestsubpop=secondlargestsubpop))
}

ccast_tree_twopass <- function(dataset,outputDir,ylabel, k=NULL, groups=NULL, s.plot=FALSE,s.save=FALSE,
                              param=NULL,fcsdes=NULL,maxiter_pruning=10,
                              subanalysis=TRUE,fcsfile=NULL,
                              transition=TRUE, max_depth = 15, min_bucket=100){
    # Performs a hierarchical clustering on the dataset if no labels are provided.
    # Then, perform successively a CCAST tree and a pruning.
    # subanalysis computes metrics on thetree and is not necessary for the purpose of fitting

    if(is.data.frame(dataset) == FALSE && is.matrix(dataset) == FALSE){
        Data_all = dataset$subsampled_dataset
    }
    else{Data_all = as.data.frame(dataset)}
    #running a hierarchical clustering to get clusters in the data
    #hier_cluster is an object output of the method
    if(is.null(groups)){
        hier_cluster <- hclust.vector(Data_all, method='ward')
        groups <- cutree(hier_cluster, k=k)
        #groups is a list of the groups to which each sample is affected
        ### Create a dataframe of data and predicted cluster assignments ######
    }
    df_data_all <- as.data.frame(cbind(Data_all,groups))
    ##Determining CCAST tree and prunning level L of CCAST tree #########
    num_groups <- max(groups)

    if (transition == TRUE) {
      init_tree<-ccast_tree_onepass(df_data_all,groups,outputDir, max_depth=max_depth, min_bucket=min_bucket)
      #Returns the shallowest tree for which each prediction is dominant in at least a cluster 
      #ccast tree has the nodes as first arg, the level as second arg, the tree as third
      ##Optimizing CCAST tree #########
      first_tree_nodes <- init_tree$tree_nodes
      first_tree_depth <- init_tree$depth

      pruned_tree <- ccast_prune_tree(tree_nodes=first_tree_nodes, dataset=df_data_all, tree_depth=first_tree_depth, outputDir=outputDir,
                                      s.save=s.save, s.plot=s.plot, min_bucket=min_bucket, maxiter=maxiter_pruning)
      #optimizes the tree to have pure populations
      targetcellsid <- pruned_tree$largestsubpop
      aliencellsid <- pruned_tree$largestalienpop
      finalD <- pruned_tree$puredata[[length(pruned_tree$puredata)]]
      ###############################
      if(subanalysis == TRUE){
      if((is.null(param)) & (is.null(fcsdes))) {
        param <- flowCore:::parameters(dataset$fcsfile)
        cDD=dataset$global_dataset[as.numeric(rownames(dataset$subsampled_dataset)),]
        out_frame <- flowFrame(cDD, param, description = description(dataset$fcsfile))
        finalD <- pruned_tree$puredata[[length(pruned_tree$puredata)]]
        ##finalD is the dataset with only points affected to pure clusters
        finalfcsD<-dataset$global_dataset[as.numeric(rownames(finalD)),]
        out_frame2 <- flowFrame(finalfcsD, param, description = description(dataset$fcsfile))
        targetD<-dataset$global_dataset[as.numeric(targetcellsid),]
        alienD<-dataset$global_dataset[as.numeric(aliencellsid),]
        out_frame3 <- flowFrame(targetD, param, description = description(dataset$fcsfile))
        out_frame4 <- flowFrame(alienD, param, description = description(dataset$fcsfile))
      } else {
        cDD <- dataset$global_dataset
        out_frame <- flowFrame(cDD, param, description = fcsdes)
        finalD <- pruned_tree$puredata[[length(pruned_tree$puredata)]]
        finalfcsD <- dataset$global_dataset[as.numeric(rownames(finalD)),]
        out_frame2 <- flowFrame(finalfcsD, param, description = fcsdes)
        targetD <- dataset$global_dataset[as.numeric(targetcellsid),]
        alienD <- dataset$global_dataset[as.numeric(aliencellsid),]
        out_frame3 <- flowFrame(targetD, param, description = fcsdes)
        out_frame4 <- flowFrame(alienD, param, description = description(dataset$fcsfile))
      }
    }
    if(s.save == TRUE) {
      obj_to_save = dataset$global_dataset
      save(obj_to_save, file=paste(outputDir, "Allexprsdata.rda", sep='/'))
      obj_to_save = dataset$subsampled_dataset
      save(obj_to_save,file=paste(outputDir, "cleandata.rda", sep='/'))
      save(groups,file=paste(outputDir, "groups.rda", sep='/'))
      write.FCS(out_frame, paste(outputDir,"origfiltereddata.fcs", sep='/'))
      write.FCS(out_frame2, paste(outputDir,"finalfiltereddata.fcs", sep='/'))
      write.FCS(out_frame3, paste(outputDir,"majorfiltereddata.fcs", sep='/'))
      write.FCS(out_frame3, paste(outputDir,"alienfiltereddata.fcs", sep='/'))
    }
    } else {
      init_tree<-ccast_tree_onepass(df_data_all,groups,outputDir,filter=FALSE, max_depth=max_depth, min_bucket=min_bucket)
      targetcellsid=1:dim(df_data_all)[1]

      finalD <- df_data_all
      finalfcsD <- dataset$global_dataset
      targetD <- dataset$global_dataset
      alienD <- dataset$global_dataset
    }
    if(s.plot == TRUE) {
      ############## Biaxial plot of tree #########
      ccast_biaxialplot(optccastreeoutput=pruned_tree,ylabel,outputDir)
      #### Heatmaps and barplots for Homogenous subpopulations ###########
      heatmapdata <- ccast_heatmaplot(optccastreeoutput=pruned_tree,outputDir)

      #### Barplots for Homogenous subpopulations ###########
      MM <- heatmapdata[[1]]
      MM2 <- heatmapdata[[2]]

      bar_plot_ccast(MM, MM2, "Barplot_for_homogenous_cells.pdf",outputDir)
      ccast_biaxialplot2(optccastreeoutput=init_tree,df_data_all,ylabel,outputDir)      
    }

    if (subanalysis == FALSE){
        return(list(initialtree=init_tree$tree,finaltree=pruned_tree$tree, groups=groups, treeheight=first_tree_depth, Finaldata=finalD))
    }else{
    
    return(list(initialtree=init_tree$tree,finaltree=pruned_tree$tree, groups=groups,
                Initialdata=dataset$subsampled_dataset,Finaldata=finalD,allfinalexprdata=finalfcsD,
                allorigdata=dataset$global_dataset,treeheight=first_tree_depth,majordata=targetD,
                majoraliendata=alienD))
  }
}

### Train neural net function ####
#### neural net code from neuralnet package ####
trainmodel<-function(input1,input2,trainprop,hidden,seed,maxit, outputDir="./output") {
    
        resp=as.data.frame(input1)
        colnames(resp)=c("Y1","Y2")
        
        rand.vars=as.data.frame(input2)
        data<-data.frame(resp,rand.vars)
        # Train-test random splitting for linear model
        #set.seed(123)
        index <- sample(1:nrow(data),round(trainprop*nrow(data)))
        train <- as.data.frame(data[index,])
        test <- as.data.frame(data[-index,])
        
        maxs <- apply(data, 2, max)
        mins <- apply(data, 2, min)
        scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))
        
        # Train-test split
        train2 <- scaled[index,]
        test2 <- scaled[-index,]
        n <- names(train2)
        
        
        require(nnet) ### Fit single-hidden-layer neural network, possibly with skip-layer connections.
        set.seed(2)
        res2=train2[,1:2]
        rand.vars2=train2[,3:dim(train2)[2]]
        mod1<-nnet(rand.vars2,res2,data=train2,size=hidden,linout=T,maxit = maxit)
        
        test3=test2[,3:dim(test2)[2]]
        nnet.predict=predict(mod1, test3)
        tiff(paste(outputDir, "Neural network prediction plots.tiff", sep="/"), height = 15, width = 15,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
        layout(matrix(c(1:4), nrow=2, byrow=TRUE),c(rep(1,4)))
        plot(test2[,1], nnet.predict[,1],
        main="Predictions vs actual transformed",
        xlab="Actual")
        plot(test2[,2], nnet.predict[,2],
        main="Predictions vs actual transformed",
        xlab="Actual")
        
        orig.pred1 <- nnet.predict[,1]*(max(data$Y1)-min(data$Y1))+min(data$Y1)
        orig.test1 <- (test2[,1])*(max(data$Y1)-min(data$Y1))+min(data$Y1)
        orig.pred2 <- nnet.predict[,2]*(max(data$Y2)-min(data$Y2))+min(data$Y2)
        orig.test2 <- (test2[,2])*(max(data$Y2)-min(data$Y2))+min(data$Y2)
        par(mfrow=c(1,2))
        plot(orig.pred1,orig.test1,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
        abline(0,1,lwd=2)
        legend('bottomright',legend='NN',pch=18,col='red', bty='n')
        
        plot(orig.pred2,orig.test2,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
        abline(0,1,lwd=2)
        legend('bottomright',legend='NN',pch=18,col='red', bty='n')
        
        dev.off();
        MSE.nn1 <- sum((orig.test1 - orig.pred1)^2)/nrow(test2)
        MSE.nn2 <- sum((orig.test2 - orig.pred2)^2)/nrow(test2)
        
    
    return(list(nnmodel=mod1,MSEerrors=c(MSE.nn1,MSE.nn2),dataout=data))
}


##### Main function: Project each marker expression on tsne-map ####
## tsnedat: tsne output data
## newdat: EMT data with 6 markers
## newdat2: Data with all markers
## antibody: antibody names of data columns
## nn: neural network output from neuralnet R package
## vor: voronoi output with polygon boundaries from tripack package
## file1: character string for sample name
maptsne<-function (tsnedat,asinhp,newdat,newdat2,antibody,nn,vor,outputDir='./output', file1="Sample.tiff",sample="Sample", window=NULL) {
    Y1=tsnedat[,1]
    Y2=tsnedat[,2]
    if (!(is.null(asinhp))) {
        Dall=asinh(newdat/asinhp)
    }else{
        Dall=newdat
    }
    
    data=as.data.frame(Dall)
    maxs <- apply(data, 2, max)
    mins <- apply(data, 2, min)
    
    scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))
    
    nnet.predict=predict(nn, scaled)
    pr.nn1 <- nnet.predict[,1]
    pr.nn2 <- nnet.predict[,2]
    
    pr.nn11 <- pr.nn1*(max(Y1)-min(Y1))+min(Y1)
    pr.nn22 <- pr.nn2*(max(Y2)-min(Y2))+min(Y2)
    
    if (is.null(window)){
        window = c(min(Y1)-5, max(Y1)+5, min(Y2)-5, max(Y2)+5)
    }
    dir.create(paste(getwd(),outputDir,sample,sep="/"))

    
    sample_pred=cbind(Y1=pr.nn11,Y2=pr.nn22)
    colnames(sample_pred)<-c("Y1","Y2")
    save(sample_pred,file=paste(outputDir, "/", sample,".rda",sep=""))
    
    curr_dir = getwd()
    setwd(paste(getwd(),outputDir,sample,sep="/"))
    ccast_tsne_plot2b(sample_pred[,1],sample_pred[,2],vor,file1=file1,sample=sample, window=window)
    
    x=sample_pred[,1]
    y=sample_pred[,2]
    if (!(is.null(asinhp))) {
        dat=asinh(newdat2/asinhp)
    }else{
        dat=newdat2
    }
    colnames(dat)=antibody
    
    
    if (!is.null(antibody)){
        for ( r in 1:length(antibody)) {
            
            tiff(paste(antibody[r],".tiff",sep=""), height = 15, width = 15,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
            
            ##if (!(r%in% c(46,47,48))){
            z <- dat[,antibody[r]]
            colvar <- asinh(z)
            scatter2D(x,y,colvar=colvar,pch=20,cex=0.4,main=antibody[r],colkey = FALSE, xlim=window[1:2], ylim=window[3:4])
            
            ###plot(vor,add=TRUE,lwd=2)
            plot(vor,add=TRUE, lwd=2,wl='tess',wp='n')
            ## }
            dev.off();
        }

    } 

    setwd(curr_dir)
    
}

maptsne2<-function (tsnedat,asinhp,newdat,newdat2,antibody,nn,vor,outputDir='./output', file1="Sample.tiff",sample="Sample", window=NULL) {
    Y1=tsnedat[,1]
    Y2=tsnedat[,2]
    if (!(is.null(asinhp))) {
        Dall=asinh(newdat/asinhp)
    }else{
        Dall=newdat
    }
    
    data=as.data.frame(Dall)
    maxs <- apply(data, 2, max)
    mins <- apply(data, 2, min)
    
    scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))
    
    nnet.predict=predict(nn, scaled)
    pr.nn1 <- nnet.predict[,1]
    pr.nn2 <- nnet.predict[,2]
    
    pr.nn11 <- pr.nn1*(max(Y1)-min(Y1))+min(Y1)
    pr.nn22 <- pr.nn2*(max(Y2)-min(Y2))+min(Y2)
    
    if (is.null(window)){
        window = c(min(Y1)-5, max(Y1)+5, min(Y2)-5, max(Y2)+5)
    }
    dir.create(paste(getwd(),outputDir,sample,sep="/"))

    
    sample_pred=cbind(Y1=pr.nn11,Y2=pr.nn22)
    colnames(sample_pred)<-c("Y1","Y2")
    save(sample_pred,file=paste(outputDir, "/", sample,".rda",sep=""))
    
    curr_dir = getwd()
    setwd(paste(getwd(),outputDir,sample,sep="/"))
    ccast_tsne_plot2b(sample_pred[,1],sample_pred[,2],vor,file1=file1,sample=sample, window=window)
    
    x=sample_pred[,1]
    y=sample_pred[,2]
    if (!(is.null(asinhp))) {
        dat=asinh(newdat2/asinhp)
    }else{
        dat=newdat2
    }
    colnames(dat)=antibody
    
    
    if (!is.null(antibody)){
        for ( r in 1:length(antibody)) {
            
            tiff(paste(antibody[r],".tiff",sep=""), height = 15, width = 15,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
            
            ##if (!(r%in% c(46,47,48))){
            z <- dat[,antibody[r]]
            colvar <- asinh(z)
            scatter2D(x,y,colvar=colvar,pch=20,cex=0.4,main=antibody[r],colkey = FALSE, xlim=window[1:2], ylim=window[3:4])
            
            plot(vor,add=TRUE,lwd=2)
            ##plot.deldir(vor,add=TRUE, lwd=2,wl='tess',wp='n')
            ## }
            dev.off();
        }

    }

    setwd(curr_dir)
    
}
                                                    
voronoi_mapping <- function(groups, ds, dbins, points_to_add=list()){
    targetclust = unique(groups)
    seeds = list()
for ( i in 1:length(targetclust)) {
    id1=which(groups==targetclust[i])
    xysum <- condense(bin(ds[id1,1], dbins[i] ), bin(ds[id1,2],dbins[i]))
    id2=which.max(xysum[,3])
    cc1=xysum[id2,1:2]
    seeds=rbind(seeds,cc1)
}
x1=seeds[,1]
y1=seeds[,2]

if(length(points_to_add) > 0){
    x1 = c(x1, points_to_add[[1]])
    y1 = c(y1, points_to_add[[2]])
}

seeds = as.matrix(seeds)

z = deldir(seeds[, 1], seeds[, 2], rw=c(min(ds[, 1]), max(ds[, 1]), min(ds[, 2]), max(ds[, 2])))
return(z)
}

ccast_tree_analyze_and_clean <- function(tree,DD, s.save, outputDir, filename="Heatmaps_Homogeneous_cells%03d.tiff") {
    #'Creates heatmap and computes cluster averages from a tree and a dataset
    
    #'Parameters:
    #'tree : a fitted clustering tree
    #'DD : a dataset
    #'s.dat : boolean, True to save
    #'outputDir : output directory to save the figures
    #'filename : the generic filename to save the heatmaps
    
  #DD is the "pure" data, without the cluster affectation
  tiff(paste(outputDir, filename, sep="/"), height = 10, width = 10,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
  
  tt1=party:::nodes(tree, unique(where(tree)))
  #tt1 is the leaves of the tree
  dD3=DD
  
  nrcolors = 100
  half = 1 + nrcolors/2
  colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
  colorpalette = colorRampPalette(colpal)(nrcolors)
  
  MM <- MM2 <- rname <- NULL
  clusterid=rep(0,dim(dD3)[1])
  # list full of zeros, of length of the dataset
  rowids<-list()
  homodata=list()
  for( j in 1:length(tt1)) {
    #Iterate over each leaf
    w1 <- which(tt1[[j]]$weights==1)
    ##index of the samples belonging to this leaf
    
    clusterid[w1]=j
    #fill the cluster id with the current id
    minv=min(dD3[w1,],na.rm=TRUE)
    #Min value of the markers for the cells in this groups 
    maxv=max(dD3[w1,],na.rm=TRUE)
    #Min value of the markers for the cells in this groups
    key_color1=seq(minv,maxv,length.out=dim(dD3)[2])
    
    rowids[[j]]=as.numeric(rownames(dD3[w1,]))
    BB=as.matrix(rbind(dD3[w1,],key_color1))
    rownames(BB)=rep(paste("Cell",j),dim(BB)[1])
    mm <- round(c(apply(BB[-dim(BB)[1],], 2, mean)),3)
    #mean of each cluster
    MM <- rbind(MM,mm)
    mm2 <- round(c(apply(BB[-dim(BB)[1],], 2, sd)),3)
    #std of each cluster
    MM2 <- rbind(MM2,mm2)
    image(x=1:dim(BB)[2],y=1:dim(BB)[1],z=as.matrix(t(BB)),col = colorpalette,xlab = paste("Subpopulationnode",tt1[[j]]$nodeID),xaxt="n",yaxt="n",ylab = "Cells",family="sans",cex.lab=1.4,main=paste("Node",tt1[[j]]$nodeID))
    mtext(rownames(BB),at=c(1:dim(BB)[1]),side=4,las=2,line=1, family= "sans", cex=0.4)
    mtext(colnames(BB),at=c(1:dim(BB)[2]),side=1,las=2,line=1, family= "sans", cex=0.4)
    homodata[[j]] <- BB[-dim(BB)[1],]
    rownames(homodata[[j]]) <- paste(rowids[[j]])
    rname <- c(rname,paste("Node",tt1[[j]]$nodeID))
  }
  rownames(MM)=rname
  rownames(MM2)=rname
  if(s.save == TRUE) {
    save(homodata,file=paste0(outputDir, "/combinedatalist.rdata"))
    save(rname,file=paste0(outputDir, "/NodeIDs.rdata"))
    save(rowids,file=paste0(outputDir, "/rowids.rdata"))
    }
  dev.off()
  return(list(ccastmeans=MM,ccastsds=MM2,
              rids=rowids,noclusters=length(tt1),
              newclusters=clusterid))
}

optimalbinsize2 <- function(pop,timeid,ds,maxbinsize,targetclust=NULL) {
    l1=length(unique(pop))
    if (l1==length(targetclust) || is.null(targetclust)) {
        # l3 = 1:l1
        l3 = unique(pop)
    } else {
        l3=targetclust
    }
    l2=length(unique(timeid))
    bs=maxbinsize
    binmatrix=matrix(0,nrow=length(l3),ncol=bs)
    colnames(binmatrix)=paste("bin",1:dim(binmatrix)[2])
    rownames(binmatrix)=paste("State",1:dim(binmatrix)[1])
    cc=NULL
    for ( i in l3) {

           for ( k in c(1:bs) ){
                id1=which(pop==i)
                id2=timeid[id1]
                optstate=as.numeric(names(table(id2)))[which.max(table(id2))]
                id3=id1[grep(optstate,id2)]
                xysum <- condense(bin(ds[id3,1], k ), bin(ds[id3,2],k))
                
                ii=match(i,l3)
                binmatrix[ii,k]= sort(xysum[,3],decreasing=T)[1]-sort(xysum[,3],decreasing=T)[2]
                id2=which.max(xysum[,3])
                cc1=xysum[id2,1:2]
                cc=rbind(cc,cc1)
        }
    }
    binmatrix1=ifelse(binmatrix==0,1,binmatrix)
    folddiff=list()
    bestbin=NULL
    for ( i in l3)  {
        fc=NULL
        ii=match(i,l3)
        for (j in 1:(dim(binmatrix1)[2]-1)) {
            fc=c(fc,binmatrix1[ii,j+1]/binmatrix1[ii,j])
        }
        folddiff[[ii]]=fc
        bestbin=c(bestbin,c(which.max(folddiff[[ii]])+1))
    }
    return(bestbin)
}

######################################### HELPER FUNCTIONS ############################################################

SPADE.addDensity.downsample<- function (infilename, cols = colid, 
                                        arcsinh_cofactor = asinhp, 
                                        kernel_mult = 5, apprx_mult = 1.5, 
                                        med_samples = 2000, comp = T,
                                        exclude_pctile = 0.01, 
                                        target_pctile = 0.05, 
                                        desired_samples = 20000) {

    in_fcs <- read.FCS(infilename, transformation = FALSE)
    in_data1 <- exprs(in_fcs)
    rownames(in_data1) <- paste(1:dim(in_data1)[1])
    in_data <- in_data1[,cols]
    if(desired_samples > dim(in_data1)[1]) {
        out_downsample2=in_data1
        params <- flowCore:::parameters(in_fcs)
        desc <- description(in_fcs)
        pd <- pData(params)
        
        firstout <- paste0("results/", basename(infilename),".downsample.fcs")
        out_frame <- flowFrame(out_downsample2, params, description = desc)
        
    } else {
    ########Estimate density#############
    density <- spade:::SPADE.density(asinh(in_data1/arcsinh_cofactor),
    kernel_mult = kernel_mult, apprx_mult = apprx_mult, med_samples = med_samples)
    if (max(density) == 0)
    warning(paste(infilename, "has degenerate densities, possibly due to many identical observations",
    sep = " "))
    in_data <- cbind(in_data, density = density)
    ########## downsample#######
    d_idx <- match("density", colnames(in_data))
    if (is.na(d_idx)) {
        stop("No density parameter in FCS file")
    }
    boundary <- quantile(in_data[, d_idx], c(exclude_pctile,
    target_pctile), names = FALSE)
    out_data <- subset(in_data, in_data[, d_idx] > boundary[1])
    density <- out_data[, d_idx]
    if (is.null(desired_samples)) {
        boundary <- boundary[2]
        out_data <- subset(out_data, boundary/density > runif(nrow(out_data)))
    }
    else if (desired_samples < nrow(out_data)) {
        density_s <- sort(density)
        cdf <- rev(cumsum(1/rev(density_s)))
        boundary <- desired_samples/cdf[1]
        if (boundary > density_s[1]) {
            targets <- (desired_samples - 1:length(density_s))/cdf
            boundary <- targets[which.min(targets - density_s >
            0)]
        }
        
        out_data <- subset(out_data, boundary/density > runif(length(density)))
    }
    if(dim(out_data)[1] == 0) {
        out_downsample <- in_data1[sample(1:dim(in_data1)[1],desired_samples,replace = FALSE),]
    } else {
    out_downsample <- out_data[,-dim(out_data)[2]]
    }
    out_downsample2 <- in_data1[as.numeric(rownames(out_downsample)),]
    params <- flowCore:::parameters(in_fcs)
    desc <- description(in_fcs)
    pd <- pData(params)
        
    dir.create(paste0(dirname(getwd()),"/downsamples/"))
    firstout <- paste0(dirname(getwd()),"/downsamples/", basename(infilename), ".downsample.fcs")
    print(basename(infilename))
    out_frame <- flowFrame(out_downsample2, params, description = desc)
    write.FCS(out_frame, firstout)
    }
    return(out_downsample2)
}
##### import function for multiple fcs files#####
fcsimport3 <- function(file,transformlogic=FALSE,asinhp=asinhp, colid,coln=NULL,
                       rown=NULL,patient,subanalysis=FALSE,subsamplesize=subsamplesize) {
    set.seed(123)
    if (grepl('fcs', file)) {
        Dall <- dat1 <- ctype <- NULL
        if(subanalysis==TRUE) {
            d0= SPADE.addDensity.downsample(file,cols = colid, 
                                            arcsinh_cofactor = asinhp,
                                            desired_samples=subsamplesize)
            dat=d0
            d2=dat[,colid]
            if(!is.null(coln)) colnames(d2)=coln
            if(!is.null(rown))  {rownames(d2) = rown
            } else {
                rownames(d2) =rep(patient,dim(d2)[1])
            }
            Dall <- D <- d2
            dat1=dat
            ctype <- rep(1,dim(d2)[1])
        } else {
            d1=read.FCS(file,transformation=transformlogic)
            dat=exprs(d1)
            d2=dat[,colid]
            if(!is.null(coln)) colnames(d2)=coln
            if(!is.null(rown))  {rownames(d2) = rown
            } else { rownames(d2) =rep(file,dim(d2)[1]) }
            Dall <- D <- d2
            dat1=dat
            ctype=rep(1,dim(d2)[1])
            
        }
        
        d1=read.FCS(file,transformation=transformlogic)
        in_data <- dat1
        params <- flowCore::parameters(d1)
        desc <- description(d1)
        
        pd <- pData(params)
        typevec <- ctype
        
        firstout <- "pooled_treatment.fcs"
        
        channel_number <- ncol(in_data) + 1
        channel_id <- paste0("$P", channel_number)
        channel_name <- "treatment"
        channel_range <- length(file)
        plist <- matrix(c(channel_name, channel_name, channel_range,
        0, channel_range - 1))
        rownames(plist) <- c("name", "desc", "range", "minRange",
        "maxRange")
        colnames(plist) <- c(channel_id)
        pd <- rbind(pd, t(plist))
        pData(params) <- pd
        out_data <- cbind(in_data, treatment = typevec)
        out_frame <- flowFrame(out_data, params, description = desc)
        keyval <- list()
        keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
        keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
        keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
        keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
        keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
        keyword(out_frame) <- keyval

        write.FCS(out_frame,firstout)

    }else{
      ff <- list.files(file)
      patientff <- grep(patient,ff)
      filelist <- paste(file,ff[patientff],sep="/")
      
      d1 <- d2 <- list()
      Dall <- dat1 <- ctype <- NULL

    
    for ( i in 1:length(filelist)) {
        if(subanalysis==TRUE) {
            d0= SPADE.addDensity.downsample(filelist[i],cols = colid, arcsinh_cofactor = asinhp,desired_samples=subsamplesize)
            dat=d0
            d2=dat[,colid]
            if(!is.null(coln)) colnames(d2)=coln
            if(!is.null(rown))  {rownames(d2) = rown
            } else {
                rownames(d2) =rep(patient,dim(d2)[1])
            }
            D=d2
            Dall=rbind(Dall,D)
            dat1=rbind(dat1,dat)
            ctype=c(ctype,rep(i,dim(d2)[1]))
            
        } else {
            d1=read.FCS(filelist[i],transformation=transformlogic)
            dat=exprs(d1)
            d2=dat[,colid]
            if(!is.null(coln)) colnames(d2)=coln
            if(!is.null(rown))  {rownames(d2) = rown
            } else { rownames(d2) =rep(file,dim(d2)[1]) }
            D=d2
            Dall=rbind(Dall,D)
            dat1=rbind(dat1,dat)
            ctype=c(ctype,rep(i,dim(d2)[1]))
            
        }
    }
    
    
    d1=read.FCS(filelist[i],transformation=transformlogic)
    in_data <- dat1
    params <- flowCore::parameters(d1)
    desc <- description(d1)
    
    pd <- pData(params)
    typevec <- ctype
    
    firstout <- "pooled_treatment.fcs"
    
    channel_number <- ncol(in_data) + 1
    channel_id <- paste("$P", channel_number, sep = "")
    channel_name <- "treatment"
    channel_range <- length(filelist)
    plist <- matrix(c(channel_name, channel_name, channel_range,
    0, channel_range - 1))
    rownames(plist) <- c("name", "desc", "range", "minRange",
    "maxRange")
    colnames(plist) <- c(channel_id)
    pd <- rbind(pd, t(plist))
    pData(params) <- pd
    out_data <- cbind(in_data, treatment = typevec)
    out_frame <- flowFrame(out_data, params, description = desc)
    keyval <- list()
    keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
    keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
    keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
    keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
    keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
    keyword(out_frame) <- keyval
    
    #write.FCS(out_frame,firstout)
}
    return(list(biomarkerdata=Dall,allexprsdata=out_data,fcsfile=out_frame))
    
}
