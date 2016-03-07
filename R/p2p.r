

tpLabels=c("a","b","c","d","e","f","g","h","i","j","k","l","m","n")
SeqDataPath<-"/users/hragbalian/desktop/goodyear/Goodyear w derotation.sav"
ExcelOutPath<-"/users/hragbalian/desktop/goodyear/new test OUT.xlsx"


p2p_wrap(SeqDataPath=SeqDataPath,
	DFA=T,
	tpLabels=tpLabels,
	ExcelOutPath=ExcelOutPath)
	

#' Path to purchase wrapper
#'
#' 


#' @export

p2p_wrap<-function(SeqDataPath,				# Path to spss sequence data, should be cleaned if there's any cleaning that it needs, first column should contain ids
			DFA = FALSE,				# Should we build a classification algorithm based
			CapSeqLength = 20,			# Should we cap the length of sequences, if so, enter numeric
			SingleState= FALSE,			# Should sequences be defined by single state transitions, or allow the same event type to repeat multiple in a row
			HowManyClusters=5:10,		# A range or scalar to return cluster solutions
			costMatrix=NULL,			# the square indel cost matrix for each type state
			tpLabels=NULL,				# A character vector of the touchpoint labels
			OtherProfileBinariesPath=NULL, # A path to an SPSS file that contains profile variable binaries, first column should be IDs.
			ExcelOutPath=NULL,
			ConvertOutToSPSS=TRUE		# Convert output to an SPSS conformable means table?
		)
	{
	
	# Load external libraries
	message("Loading in built-in function")
	library(TraMineR)
	library(foreign)
	library(cluster)
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
	
	maxtpCode=max(Data[,-1],na.rm=T)
	
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
				eval(parse(text=paste("testModel<-lda(cluster",HowManyClusters[hmc],"~DFAData[,KeepFromStep])",sep="")))
				
				# store DFA versin of cluster
				eval(parse(text=paste("DFA_cluster",HowManyClusters[hmc],"<-predict(testModel)$class",sep="")))
				
				# retrieve cross table, predicted vs observed
				eval(parse(text=paste("AccTable<-table(cluster",HowManyClusters[hmc],",DFA_cluster",HowManyClusters[hmc],")",sep="")))
				
				# retrieve % accuracy
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
		
		if (!DFA) ClusterNames<-paste("cluster",HowManyClusters,sep="")
		if (DFA) ClusterNames<-paste("DFA_cluster",HowManyClusters,sep="")
		if (ConvertOutToSPSS) SPSSOut<-list()
		
		BindedTablesStore<-list()
		
		for (cN in ClusterNames) {
			
			currSol<-eval(parse(text=cN))
			maxCurrSol<-max(as.numeric(as.character(currSol)))
			
			FirstList<-list()
			LastList<-list()
			TotalList<-list()
			TransitionList<-list()
			
			BaseSizes<-list()
			for (rNg in 1:(maxCurrSol+1)) {
				
				if (rNg<=maxCurrSol) currClusSeqs<-SeqData[which(currSol==rNg),]
				if (rNg==(maxCurrSol+1)) currClusSeqs<-SeqData
				
				BaseSizes[[rNg]]<-dim(currClusSeqs)[1]
				
				currClusSeqsLength<-seqlength(currClusSeqs)

				# First and Last event
				FirstEventTable<-round(table(factor(currClusSeqs[,1],levels=1:maxtpCode))/dim(currClusSeqs)[1],digits=5)
				LastEventTable<-matrix(,ncol=1,nrow=length(currClusSeqsLength))
					for (i in 1:length(currClusSeqsLength)) LastEventTable[i]<-currClusSeqs[i,currClusSeqsLength[i]]
					LastEventTable<-table(factor(LastEventTable,levels=1:maxtpCode))/dim(currClusSeqs)[1]
				
				# Total frequency
				TotalPercent<-seqstatf(currClusSeqs)$Percent/100
				
				# Transition rates
				TransitionRate<-matrix(t(seqtrate(currClusSeqs)))
				rownames(TransitionRate)<-matrix(t(seqetm(currClusSeqs)))
				
				# Store
				FirstList[[rNg]]<-FirstEventTable
				LastList[[rNg]]<-LastEventTable
				TotalList[[rNg]]<-TotalPercent
				TransitionList[[rNg]]<-TransitionRate
				
				}
				
				# Store
				FirstStore<-as.data.frame(Reduce("cbind",FirstList))
					rownames(FirstStore)<-paste("FEV",1:length(tpLabels)," ",tpLabels,sep="")
					colnames(FirstStore)<-c(paste("JourneyType_",1:maxCurrSol,sep=""),"Total")
				LastStore<-as.data.frame(Reduce("cbind",LastList))
					rownames(LastStore)<-paste("LEV",1:length(tpLabels)," ",tpLabels,sep="")
					colnames(LastStore)<-c(paste("JourneyType_",1:maxCurrSol,sep=""),"Total")
				TotalStore<-as.data.frame(Reduce("cbind",TotalList))
					rownames(TotalStore)<-paste("TEV",1:length(tpLabels)," ",tpLabels,sep="")
					colnames(TotalStore)<-c(paste("JourneyType_",1:maxCurrSol,sep=""),"Total")
				TransitionStore<-as.data.frame(Reduce("cbind",TransitionList))
					colnames(TransitionStore)<-c(paste("JourneyType_",1:maxCurrSol,sep=""),"Total")
				
				# process labels
				TransSplit<-strsplit(rownames(TransitionStore),">","")
				for (i in 1:length(TransSplit)) TransSplit[[i]]<-paste(tpLabels[as.numeric(TransSplit[[i]])],collapse=">")
				rownames(TransitionStore)<-paste("TRA",1:length(unlist(TransSplit))," ",unlist(TransSplit),sep="")
			
				# Bind the tables together		
				BindedTables<-rbind(FirstStore,LastStore,TotalStore,TransitionStore)
				BindedTablesStore[[cN]]<-BindedTables
			

			# If convert to SPSS out is specified, convert the means table to SPSS conformable
		
				if (ConvertOutToSPSS) {
					Out<-list()
					for (j in 1:dim(BindedTables)[2]) Out[[j]]<-cbind(BindedTables[,j],BaseSizes[[j]],0)
					Out<-Reduce("cbind",Out)
					rownames(Out)<-rownames(BindedTables)
					
					BlankRow<-rep("",(maxCurrSol+1)*3)
					SecondRow<-BlankRow
						SecondRow[1]<-cN
					ThirdRow<-BlankRow
						ThirdRow[seq(1,(maxCurrSol+1)*3,by=3)]<-c(1:maxCurrSol,"Total")
			
					FourthRow<-rep(c("Mean","N","Std. Deviation"),maxCurrSol+1)
		
					StoreMeans<-rbind(BlankRow,SecondRow,ThirdRow,FourthRow,Out)
		
					rownames(StoreMeans)[1:4]<-c("Report","a","b","c")
					colnames(StoreMeans)<-BlankRow
					
					SPSSOut[[cN]]<-StoreMeans
				}
		
		}
		
		
		if (!is.null(ExcelOutPath))	{
			library(WriteXLS)
			for (spj in 1:length(SPSSOut)) eval(parse(text=paste("SPSS",ClusterNames[spj],"<-as.data.frame(SPSSOut[[spj]])",sep="")))
			WriteXLS(paste("SPSS",ClusterNames,sep=""),ExcelOutPath,row.names=T,col.names=F)
			}
			
		Return<-list()
			Return[["Main Table"]]<-BindedTablesStore
			if (ConvertOutToSPSS) Return[["SPSS Main Table"]]<-SPSSOut
			return(Return)		
	}

