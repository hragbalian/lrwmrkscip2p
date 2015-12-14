







p2pWrap<-function(
			SeqDataPath,			# Path to spss sequence data, should be cleaned if there's any cleaning that it needs, first column should contain ids
			DFA = FALSE,			# Should we build a classification algorithm based
			CapSeqLength = 20,		# Should we cap the length of sequences, if so, enter numeric
			SingleState= FALSE,		# Should sequences be defined by single state transitions, or allow the same event type to repeat multiple in a row
			HowManyClusters<-5:10,	# A range or scalar to return cluster solutions
			costMatrix=NULL			# the square indel cost matrix for each type state
		)
	{
	
	# Load external libraries
	message("Loading in built-in function")
	library(TraMineR)
	library(foreign)
	library(cluster)
	library(WriteXLS)
	library(MASS)
	library(klaR)
	
	# Load custom functions
	# Extract the classification coefficients  
    	coeffLDA<-function(LDAmodel,DV)
    		{
      
      ModelData<-model.matrix(LDAmodel)[,-1]
      
      priorsGet<-LDAmodel$prior
      
      gr <- length(table(DV))   ## groups might be factors or numeric
      v <- dim(ModelData)[2] ## variables
      m <- LDAmodel$means ## group means
      
      w <- array(NA, dim = c(v, v, gr))
      
      for(i in 1:gr){
        tmp <- scale(subset(ModelData, DV == unique(DV)[i]), scale = FALSE)
        w[,,i] <- t(tmp) %*% tmp
      }
      
      W <- w[,,1]
      for(i in 2:gr)
        W <- W + w[,,i]
      
      V <- W/(nrow(ModelData) - gr)
      iV <- solve(V)
      
      class.funs <- matrix(NA, nrow = v + 1, ncol = gr)
      colnames(class.funs) <- paste(names(table(DV)))
      rownames(class.funs) <- c("constant", colnames(ModelData))
      
      
      for(i in 1:gr) {
        class.funs[1, i] <- -.5 * (t(m[i,]) %*% iV %*% (m[i,])) + log(priorsGet)[i] # Intercept
        class.funs[2:(v+1) ,i] <- iV %*% (m[i,])
      }
      return(class.funs)
    }
  	
	# Crate Dummy variables out of the factor
		CreateDummy<-function(CategoricalVector,Varname,CreateBaseline=1) 
			{
		UniqueCategories<-unique(CategoricalVector)
		Store<-matrix(0,ncol=length(UniqueCategories),nrow=length(CategoricalVector))
		for (i in 1:length(CategoricalVector)) {
		  Store[i,which(UniqueCategories%in%CategoricalVector[i])]<-1
		}
	
		  #Process labels
		  UniqueCategories<-as.numeric(gsub("-","_", UniqueCategories))
		  Store<-Store[,order(UniqueCategories)]
		  UniqueCategories<-UniqueCategories[order(UniqueCategories)]
		  colnames(Store)<-paste(Varname,UniqueCategories,sep="_")
		  if (CreateBaseline!=0) {
			tempColnames<-colnames(Store)[-CreateBaseline]
			Store<-as.matrix(Store[,-CreateBaseline])
			colnames(Store)<-tempColnames
		  }
	  
		return(Store)
	
	  }
	
	
	# Load data
	message("Loading in data")
	suppressWarnings(Data<-read.spss(SeqDataPath, use.value.labels = F, to.data.frame = T))
	
		# Remove rows with only missing data
		if (length(which(apply(Data[,-c(1,2)],1,function(x) all(is.na(x)))))>0) Data<-Data[-which(apply(Data[,-c(1,2)],1,function(x) all(is.na(x)))),]
	
		# If there's been a cap placed on sequence length, remove sequences greater than that cap
		if (!is.null(CapSeqLength)) {
			Data<-Data[which(apply(Data[,-1],1,function(x) table(!is.na(x))[2]<=CapSeqLength)),]
			Data<-Data[,c(1:(CapSeqLength+1))]
			}
	
	
	# Start building sequences
		if (!is.null(CapSeqLength)) SeqData<-seqdef(Data,var=colnames(Data)[2:(CapSeqLength+1)])
		if (is.null(CapSeqLength)) SeqData<-seqdef(Data[,-1])
	
	# Check to see if single state is turned on, and if so, revise sequence data
		if (SingleState) SeqData<-seqdss(SeqData)

	# Check to see if there are any single-step sequences.  If so, remove them from the sequences
		Any1<-any(seqlength(SeqData)%in%1)	
		if (Any1) {
			message("Removing single-step sequences")
			Which1<-which(seqlength(SeqData)%in%1)	
			# Truncate sequences to those that are multi-step and those that are single step
			SeqDataSingle<-SeqData[Which1,] # can return this later
			SeqData<-SeqData[-(Which1),]
			}

	###############################	
	## Start of cluster analysis ##
	###############################
		
		if (is.null(costMatrix)) Cost<-seqsubm(SeqData, method = "CONSTANT", with.miss = TRUE)
		
		message("Calculating distances between paths")
		Distances<-seqdist(SeqData,method="OM",indel=1,norm=T,sm=Cost,with.missing=T)

		message("Clustering paths")
		clusterward <- agnes(Distances, diss = TRUE, method = "ward")

		message(paste("Creating", HowManyClusters[1], "to",HowManyClusters[length(HowManyClusters)], "cluster groups"))
		for (i in 1:length(HowManyClusters)) {
			eval(parse(text=paste("cluster",HowManyClusters[i]," <- cutree(clusterward, k=", HowManyClusters[i],")",sep="")))
			}
		
		# If DFA is true, generate an algorithm for each of the clustering solutions based on length, first event
		
		if (DFA) {
			
			# Retrieves indvidual-level measurements on the sequences
			Length<-seqlength(SeqData)
			FirstEvent<-CreateDummy(SeqData[,1],"FirstEvent")
			LastEvent<-matrix(,ncol=1,nrow=length(Length))
				for(i in 1:length(Length)) LastEvent[i]<-SeqData[i,Length[i]]
				LastEvent<-CreateDummy(LastEvent,"LastEvent")
			CountOfEachEvent<-seqistatd(SeqData)
			
			clusterModelAccuracy<-list()
			clusterModelCoefs<-list()
			for (hmc in 1:length(HowManyClusters)) {
				
				message(paste("Evaluating DFA for cluster",HowManyClusters[hmc],sep=""))
				
				# Check to see if any given sequence has a "common" subsequence of any given type
				# The check will need to be returned
				storeSubSeqMatch<-list()
				for (cl in 1:HowManyClusters[hmc]) {
					eval(parse(text=paste("baseObject<-seqefsub(seqecreate(SeqData[cluster",HowManyClusters[hmc],"%in%",cl,",]),minSupport=50)",sep="")))
					baseObject$seqe<-seqefsub(seqecreate(SeqData),minSupport=50)$seqe
		
					subseqCountCheck<-seqeapplysub(baseObject,method="presence")
					storeSubSeqMatch[[paste("cluster",cl,"of",HowManyClusters[hmc],sep="")]]<-apply(subseqCountCheck,1,sum)
					}
		
				SubSeqPresenceType<-Reduce("cbind",storeSubSeqMatch)
				colnames(SubSeqPresenceType)<-names(storeSubSeqMatch)
				
				DFAData<-cbind(FirstEvent,LastEvent,CountOfEachEvent,SubSeqPresenceType)
				
				eval(parse(text=paste("StepDFA<-greedy.wilks(cluster",HowManyClusters[hmc],"~DFAData)",sep="")))
				KeepFromStep<-gsub("DFAData","",as.character(StepDFA[[1]]$vars))
				
				# build model with gathered information
				eval(parse(text=paste("AccTable<-table(cluster",HowManyClusters[hmc],",predict(lda(cluster",HowManyClusters[hmc],"~DFAData[,KeepFromStep]))$class)",sep="")))
				clusterModelAccuracy[[paste("cluster",HowManyClusters[hmc],sep="")]]<-sum(diag(AccTable))/sum(AccTable)
			
				# store coefficients
				eval(parse(text=paste("tempCoefs<-coeffLDA(testModel,cluster",HowManyClusters[hmc],")",sep="")))
				rownames(tempCoefs)<-gsub("[]]","",gsub("[[, ]","",gsub("KeepFromStep","",gsub("DFAData","",rownames(tempCoefs)))))
				clusterModelCoefs[[paste("cluster",HowManyClusters[hmc],sep="")]]<-tempCoefs
			
				}
			}
		
	#################################################
	## Report the various profiles on the clusters ##	
	#################################################	
		
		ClusterNames<-paste("cluster",HowManyClusters,sep="")
		
		
		
				
	}



# Example
#SeqDataPath = "/users/hragbalian/desktop/goodyear/Goodyear w derotation.sav"
#Data<-Data[,-2]



	
	
	


	
#	IDsWithCodes<-cbind(Data[-(Which1),1],clusterK)
#		colnames(IDsWithCodes)<-c("Serial","Cluster")
#		IDsWithCodes<-as.data.frame(IDsWithCodes)
	
#	SingleStepIDsWithCodes<-cbind(Data[Which1,1],dssSeqDataSingle)
#	colnames(SingleStepIDsWithCodes)<-"Serial"
#	SingleStepIDsWithCodes<-as.data.frame(SingleStepIDsWithCodes)
	
#	WriteXLS(c("SeqType1df","SeqType2df","SeqType3df","SeqType4df",
#			"ProtoType1","ProtoType2","ProtoType3","ProtoType4",
#			"TransType1","TransType2","TransType3","TransType4",
#			"IDsWithCodes","SingleStepIDsWithCodes"),
#		paste(Directory,"Sequence Typing.xlsx",sep=""),row.names=T)




# Follow ups.
#NewData<-read.spss("/users/hragbalian/desktop/goodyear/Goodyear Follow Up Sequence Data.sav", use.value.labels = F, to.data.frame = T)
#NewData<-NewData[-which(apply(NewData[,-c(1,2)],1,function(x) all(is.na(x)))),]
#NewData<-NewData[which(apply(NewData[,-c(1,71)],1,function(x) table(!is.na(x))[2]<=20)),]
#NewData<-NewData[,1:21]

#SeqNewData<-seqdef(NewData,var=colnames(NewData)[2:21])
#dssSeqNewData<-seqdss(SeqNewData)
#dssSeqNewData<-dssSeqNewData[,1:19]




	# Add proto sequences from original analysis 
#	dssSeqNewData<-rbind(ProtoType1[1,],ProtoType2[1,],ProtoType3[1,],ProtoType4[1,],dssSeqNewData)

#	Distances<-seqdist(dssSeqNewData,method="OM",indel=1,norm=T,sm=Cost,with.missing=T)


#	ClusterAssign<-apply(Distances[1:4,],2,which.min)
#	ClusterAssign<-ClusterAssign[5:length(ClusterAssign)]
	
#	Return<-cbind(NewData$Respondent_Serial,ClusterAssign)
	
#	write.csv(Return,"/users/hragbalian/desktop/goodyear/Goodyear follow up - cluster assignments.csv")
	


# Second follow up. 
#NewData2<-read.spss("/users/hragbalian/desktop/goodyear/Goodyear P2P round2 follow up.sav", use.value.labels = F, to.data.frame = T)
#NewData2<-NewData2[-which(apply(NewData2[,-c(1,2)],1,function(x) all(is.na(x)))),]
#NewData2<-NewData2[which(apply(NewData2[,-c(1,71)],1,function(x) table(!is.na(x))[2]<=20)),]
#NewData2<-NewData2[,1:21]

#SeqNewData2<-seqdef(NewData2,var=colnames(NewData2)[2:21])
#dssSeqNewData2<-seqdss(SeqNewData2)
#dssSeqNewData2<-dssSeqNewData2[,1:19]

# Add proto sequences from original analysis 
#	dssSeqNewData2<-rbind(ProtoType1[1,],ProtoType2[1,],ProtoType3[1,],ProtoType4[1,],dssSeqNewData2)

#	Distances<-seqdist(dssSeqNewData2,method="OM",indel=1,norm=T,sm=Cost,with.missing=T)


#	ClusterAssign<-apply(Distances[1:4,],2,which.min)
#	ClusterAssign<-ClusterAssign[5:length(ClusterAssign)]
	
#	Return<-cbind(NewData2$Respondent_Serial,ClusterAssign)
	
	

#	write.csv(Return,"/users/hragbalian/desktop/goodyear/Goodyear second follow up - cluster assignments.csv")



    