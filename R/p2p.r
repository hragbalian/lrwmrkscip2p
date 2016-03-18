

SeqDataPath<-"/users/hragbalian/desktop/coty p2p/R151081_RAW Journey Data.sav"
ExcelOutPath<-"/users/hragbalian/desktop/coty p2p/new test OUT.xlsx"
OtherProfileBinariesPath="/users/hragbalian/desktop/coty p2p/Profile Variables.sav"



OUTPUT<-p2p_wrap(SeqDataPath=SeqDataPath,
	DFA=TRUE,
	CapSeqLength = 20,
	SingleState=FALSE,
	HowManyClusters=3:7,
	costMatrix=NULL,
	ConvertOutToSPSS=TRUE,
	ExcelOutPath=ExcelOutPath,
	OtherProfileBinariesPath=OtherProfileBinariesPath,
	Group=list("DV_Category_PipeIn",2),
	skewThreshold=.3)
	


	MeansTable<-TransitionStore

	assigner<-function(MeansTable,RangeThresh) {
				rownames(MeansTable)<-lapply(lapply(strsplit(rownames(MeansTable)," ",""),function(x) x[-1]),paste, collapse=" ")
				store<-list(); for (go in 1:(dim(MeansTable)[2]-1)) store[[go]]<-""
				Range<-apply(MeansTable[,1:(dim(MeansTable)[2]-1)],1,function(x) diff(range(x)))	
				Candidates<-which(Range>RangeThresh)
				tempMeansTable<-MeansTable[Candidates,1:(dim(MeansTable)[2]-1)]
				MeansTableMaxes<-apply(tempMeansTable,1,max)
				for (mtm in 1:dim(tempMeansTable)[1]) {
					whichGroup<-which(tempMeansTable[mtm,]%in%MeansTableMaxes[mtm])
					for (wg in 1:length(whichGroup)) {
						if (store[[whichGroup[wg]]][1]!="") store[[whichGroup[wg]]]<-c(store[[whichGroup[wg]]],rownames(tempMeansTable)[mtm])
						if (store[[whichGroup[wg]]][1]=="") store[[whichGroup[wg]]]<-rownames(tempMeansTable)[mtm]
						}
					}
				return(store)
			}



	assigner<-function(MeansTable,Thresh) {
		rownames(MeansTable)<-lapply(lapply(strsplit(rownames(MeansTable)," ",""),function(x) x[-1]),paste, collapse=" ")
		store<-list(); storeMeans<-list()
		
		tempMeansTable<-MeansTable[,1:(dim(MeansTable)[2]-1)]

		for (mtm in 1:dim(tempMeansTable)[2]) {
			store[[mtm]]<-rownames(MeansTable[which(tempMeansTable[,mtm]>Thresh),])	
			storeMeans[[mtm]]<-MeansTable[which(tempMeansTable[,mtm]>Thresh),mtm]
			}
		return(list(Arcs=store,ArcSizes=storeMeans))
	}
			
#' Path to purchase wrapper
#'
#' @param SeqDataPath A Character vector to an SPSS data file, Path to spss sequence 
#' data, should be cleaned if there's any cleaning that it needs, first column should contain ids.
#' 
#' @param DFA A logical. Default is True.  Should the cluster be based off of a DFA algorithm. If 
#' False, then cluster is derived from the ward hierarchical clustering algorithm. The DFA
#' algorithm is based off 4 characteristics of sequences: First event, Last event, 
#' frequency of each event type, and the presence/frequency of the top 10 common subsequences
#' for each cluster type. If DFA is True, additional information is released in the return
#' statement to reproduce the clusters. 
#'
#' @param CapSeqLength A numeric. Identifies how many events to cap sequences by, if any. 
#'
#' @param SingleState A logical. Should sequences be defined by single state transitions, or 
#' allow the same event type to repeat multiple in a row.
#'
#' @param HowManyClusters A numeric vector or scalar. Iterates through cluster solution sizes. 
#' default is 5:10.
#' 
#' @param costMatrix An nxn cost matrix, overriding the default costs. 
#' 
#' @param tpLabels

#' @export

p2p_wrap<-function(SeqDataPath,		
			DFA = TRUE,				
			CapSeqLength = 20,
			SingleState= FALSE,			
			HowManyClusters=5:10,		
			costMatrix=NULL,
			OtherProfileBinariesPath=NULL, # A path to an SPSS file that contains profile variable binaries, first column should be IDs.
			ExcelOutPath=NULL,
			ConvertOutToSPSS=TRUE,		# Convert output to an SPSS conformable means table?
			Group=NULL,				# Either Null, or a two slot list, one with variable name, one with variable value to select on, e.g. list("DV_Category_PipeIn",2)
			skewThreshold=.3
		)
	{
	
	# Load external libraries
	message("Loading in built-in function")
	library(TraMineR)
	library(foreign)
	library(cluster)
	library(MASS)
	library(klaR)
	library(igraph)
	library(networkD3)
	
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
	
	# retrieve event type labels
	tpLabels<-names(attr(Data[,2],"value.labels"))[order(as.numeric(attr(Data[,2],"value.labels")))]
	
	
	if (!is.null(Group)) { 
		# Validate the Group
		if (!is.character(Group[[1]]) || !is.numeric(Group[[2]])) stop("Check your Group object: first slot should be character, second slot numeric")
		Data<-Data[Data[,Group[[1]]]==Group[[2]],]
		}
	
	
		# Remove rows with only missing data
		if (length(which(apply(Data[,-c(1,2)],1,function(x) all(is.na(x)))))>0) Data<-Data[-which(apply(Data[,-c(1,2)],1,function(x) all(is.na(x)))),]
	
		# If there's been a cap placed on sequence length, remove sequences greater than that cap
		if (!is.null(CapSeqLength)) {
			Data<-Data[which(apply(Data[,-1],1,function(x) table(!is.na(x))[2]<=CapSeqLength)),]
			Data<-Data[,c(1:(CapSeqLength+1))]
			}
	
		# Record the ids of everyone who remains
		IDs<-Data[,1]
	
		# Retrieve the max code for the touchpoints
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
			IDs<-IDs[-(Which1)]
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
				if (length(table(LastEvent))==1) for(i in 1:length(Length)) LastEvent[i]<-SeqData[i,Length[i]-1]
				LastEvent<-CreateDummy(LastEvent,"LastEvent")
			CountOfEachEvent<-seqistatd(SeqData)
				colnames(CountOfEachEvent)<-paste("CountOfEventType",colnames(CountOfEachEvent),sep="")
			
			# Check for any constants and remove them
			anyConstant<-any(apply(CountOfEachEvent,2,var)==0)
			if (anyConstant) CountOfEachEvent<-CountOfEachEvent[,-which(apply(CountOfEachEvent,2,var)==0)]
			
			# Remove highly correlated variables at > .6
				#DFAData<-cbind(FirstEvent,LastEvent,CountOfEachEvent,SubSeqPresenceType)
				DFAData<-cbind(FirstEvent,LastEvent,CountOfEachEvent)
				corDFAData<-cor(DFAData)
				diag(corDFAData)<-0
				cordumDFAData<-matrix(0,dim(corDFAData)[1],dim(corDFAData)[2])
					colnames(cordumDFAData)<-colnames(corDFAData)
					rownames(cordumDFAData)<-colnames(corDFAData)
				cordumDFAData[abs(corDFAData)>.6]<-1
					corGraph<-graph_from_adjacency_matrix(cordumDFAData)
					corGraph<-as_edgelist(corGraph)

				# randomly remove one of the problematic pairs
				DFAData<-DFAData[,-which(colnames(DFAData)%in%apply(corGraph,1,sample,size=1))]
			
		
			clusterModelAccuracy<-list()
			clusterModelCoefs<-list()
			clusterSubSeqs<-list()
			for (hmc in 1:length(HowManyClusters)) {
				
				message(paste("Evaluating DFA for cluster",HowManyClusters[hmc],sep=""))
				
				# Check to see if any given sequence has a "common" subsequence of any given type
				# The check will need to be returned
				#storeSubSeqMatch<-list()
				#storeClustSubSeq<-list()
				#SeqEData<-seqecreate(SeqData) # Initialize event sequence data
				
				#for (cl in 1:HowManyClusters[hmc]) {
				#	eval(parse(text=paste("baseObject<-seqefsub(seqecreate(SeqData[cluster",HowManyClusters[hmc],"%in%",cl,",]),minSupport=2)",sep="")))
				#	if (length(baseObject$subseq)>10) {
				#		baseObject$subseq<-baseObject$subseq[1:10] # keep top 10 subseqs for the cluster
				#		baseObject$data<-baseObject$data[1:10,]
				#		}
					
				#	subseqCountCheck<-seqeapplysub(baseObject,method="presence")
					
				#	match(rownames(subseqCountCheck),unlist(lapply(unlist(SeqEData),as.character)))
					
				#	storeClustSubSeq[[paste("SubSeqFromCluster",cl,sep="")]]<-as.character(baseObject$subseq)
				#	storeSubSeqMatch[[paste("SubSeqCountFromCluster",cl,sep="")]]<-apply(subseqCountCheck,1,sum)
				#	}
				
				#storeClustSubSeq<-lapply(storeClustSubSeq , function(x) c(x,rep(NA,10-length(x))))
				#storeClustSubSeq<-Reduce("rbind",lapply(storeClustSubSeq,as.character))
				#	rownames(storeClustSubSeq)<-paste("TopSubSeqForClust",1:HowManyClusters[hmc],sep="")
				#	clusterSubSeqs[[paste("cluster",HowManyClusters[hmc],sep="")]]<-storeClustSubSeq
					
				#SubSeqPresenceType<-Reduce("cbind",storeSubSeqMatch)
				#colnames(SubSeqPresenceType)<-names(storeSubSeqMatch)		
								
				eval(parse(text=paste("StepDFA<-greedy.wilks(cluster",HowManyClusters[hmc],"~DFAData)",sep="")))
				KeepFromStep<-gsub("DFAData","",as.character(StepDFA[[1]]$vars))
				
				# build model with gathered information
				eval(parse(text=paste("testModel<-lda(cluster",HowManyClusters[hmc],"~DFAData[,KeepFromStep])",sep="")))
				
				# store DFA version of cluster
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
		
		# Import and process OtherProfileBinariesPath, if included
		if (!is.null(OtherProfileBinariesPath)) {
				suppressWarnings(OthProfileData<-read.spss(OtherProfileBinariesPath, use.value.labels = F, to.data.frame = T))
				OPLabels<-as.character(attr(OthProfileData,"variable.labels"))
				
				FilteredProfileData<-OthProfileData[match(IDs,OthProfileData[,1]),]
				FilteredProfileData<-FilteredProfileData[-1]
				OPVarNames<-colnames(FilteredProfileData)
				OPLabels<-OPLabels[-1]
				}
				
		if (!DFA) ClusterNames<-paste("cluster",HowManyClusters,sep="")
		if (DFA) ClusterNames<-paste("DFA_cluster",HowManyClusters,sep="")
		if (ConvertOutToSPSS) SPSSOut<-list()
		
		BindedTablesStore<-list()
		
		storeProto<-list()
		storeClusterSolutions<-list()
		masterStoreGraphs<-list()
		masterStoreGraphTables<-list()
		storeBaseSizes<-list()
		
		for (cN in ClusterNames) {
			
			## Compute means for profiling
			
			currSol<-eval(parse(text=cN))
			maxCurrSol<-max(as.numeric(as.character(currSol)))
			storeClusterSolutions[[cN]]<-cbind(IDs,currSol)
			
			## Identify PROTOTYPICAL sequences
			
			DistanceFromCenter<-disscenter(Distances,currSol)
			#DistanceFromCenter[which(currSol==1)]
			
			tempStoreProto<-list()
			for (rNg in 1:maxCurrSol) tempStoreProto[[rNg]]<-unique(head(SeqData[currSol%in%rNg,][order(DistanceFromCenter[currSol%in%rNg]),],20))
			storeProto[[paste("cluster",cN,sep="")]]<-tempStoreProto
			
			
			# ... continue with profiling	
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
					if (var(LastEventTable)==0) for (i in 1:length(currClusSeqsLength)) LastEventTable[i]<-currClusSeqs[i,currClusSeqsLength[i]-1]
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
				for (i in 1:length(TransSplit)) {
					pasterwhich<-as.numeric(TransSplit[[i]])
					if (length(pasterwhich)>1) TransSplit[[i]]<-paste(tpLabels[pasterwhich],collapse=">")
					if (length(pasterwhich)==1) TransSplit[[i]]<-paste(tpLabels[pasterwhich],tpLabels[pasterwhich],sep=">")
					}
				rownames(TransitionStore)<-paste("TRA",1:length(unlist(TransSplit))," ",unlist(TransSplit),sep="")
			
				# Bind the tables together		
				BindedTables<-rbind(FirstStore,LastStore,TotalStore,TransitionStore)
				BindedTablesStore[[cN]]<-BindedTables
			
				# Piece sequence map together
				FirstAssigned<-assigner(FirstStore,.1)
				Assigned<-assigner(TransitionStore,skewThreshold)
				
				storeGraphs<-list()
				storeGraphTables<-list()
				
				for (group in 1:maxCurrSol) {

					Converted<-balialab::convert_edgelist(Reduce("rbind",strsplit(Assigned$Arcs[[group]], ">","")))
				
					Links<-as.data.frame(cbind(Converted$Numeric-1,1))
						colnames(Links)<-c("from","to","value")
					Nodes<-as.data.frame(cbind(Converted$Character,1,1))
						colnames(Nodes)<-c("name","group","size")
						Nodes$group<-factor(Nodes$group,levels=c(1:3))
				
					#assign start
					Nodes$group[which(Nodes$name%in%FirstAssigned$Arcs[[group]])]<-3
				
					#assign size
					Nodes$size<-TotalStore[match(Nodes$name,tpLabels),group]*100
				
					#assign arc sizes
					Links$Value<-Assigned$ArcSizes[[group]]*100
					
					#consolidate tables
					LinksNodes<-list(Links=Links,Nodes=Nodes)
					storeGraphTables[[group]]<-LinksNodes
					
					storeGraphs[[group]]<-forceNetwork(Links,Nodes,Source = "from",
						Target = "to", Value = "value", NodeID = "name",
						Nodesize = "size",
							 Group = "group",opacity = 1, charge=-1000,opacityNoHover=1,fontSize=10,zoom=TRUE)	
					
				}
					
				masterStoreGraphs[[cN]]<-storeGraphs
				masterStoreGraphTables[[cN]]<-storeGraphTables
				storeBaseSizes[[cN]]<-cbind(1:maxCurrSol,unlist(BaseSizes)[-length(unlist(BaseSizes))])
				
			###########################################################################################
			## If there are additional profile variables to append, add those and obtain their means ##
			###########################################################################################
			
			if (!is.null(OtherProfileBinariesPath)) {	
				#temp data
				#OthProfileData<-as.data.frame(cbind(Data[,1],matrix(sample(c(0,1),dim(Data)[1]*dim(Data)[2],replace=T,prob=c(.8,.2)),nrow=dim(Data)[1],ncol=dim(Data)[2])))
	
				storeOtherProfileMeans<-list()
				for (op in 1:dim(FilteredProfileData)[2]){
					eval(parse(text=paste("storeOtherProfileMeans[[OPVarNames[op]]]<-tapply(FilteredProfileData[,op],",cN,",mean,na.rm=T)",sep="")))
					}
				
				#Reduce
				storeOtherProfileMeans<-Reduce("rbind",storeOtherProfileMeans)
				storeOtherProfileMeans<-cbind(storeOtherProfileMeans,apply(FilteredProfileData,2,mean,na.rm=T))
					rownames(storeOtherProfileMeans)<-paste(OPVarNames," ",OPLabels,sep="")
					colnames(storeOtherProfileMeans)<-c(paste("JourneyType_",1:maxCurrSol,sep=""),"Total")
					
				# append the means
					BindedTables<-rbind(BindedTables,storeOtherProfileMeans)
					BindedTablesStore[[cN]]<-rbind(BindedTablesStore[[cN]],storeOtherProfileMeans)
				}
			
			
				######################################################################################
				## If convert to SPSS out is specified, convert the means table to SPSS conformable ##
				######################################################################################
					
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
	
		
		} # close main cluster loop
		
		
		############
		## Return ##
		############
				
		if (!is.null(ExcelOutPath))	{
			library(WriteXLS)
			for (spj in 1:length(SPSSOut)) eval(parse(text=paste("SPSS",ClusterNames[spj],"<-as.data.frame(SPSSOut[[spj]])",sep="")))
			WriteXLS(paste("SPSS",ClusterNames,sep=""),ExcelOutPath,row.names=T,col.names=F)
			}
			
		Return<-list()
			Return[["MainTable"]]<-BindedTablesStore
			Return[["Prototypes"]]<-storeProto
			Return[["TPLabels"]]<-tpLabels
			Return[["Solutions"]]<-storeClusterSolutions
			Return[["Graphs"]]<-masterStoreGraphs
			Return[["GraphTables"]]<-masterStoreGraphTables
			Return[["BaseSizes"]]<-storeBaseSizes
			if (ConvertOutToSPSS) Return[["SPSSMainTable"]]<-SPSSOut
			if (DFA) {
				Return[["DFACoefs"]]<-clusterModelCoefs
				Return[["DFAModelAccuracy"]]<-clusterModelAccuracy
				Return[["SubSeqsForDFA"]]<-clusterSubSeqs
				}
				
		return(Return)		
	}



#' Survival based drivers for sequence data


SPSSPath<-"/users/hragbalian/desktop/coty p2p/basic.sav" # should have both
TPVarNames<-c(paste("Toucpoint_00",1:9,sep=""),paste("Toucpoint_0",10:50,sep=""))
TimeVarNames<-c(paste("SliderVal_00",1:9,sep=""),paste("SliderVal_0",10:50,sep=""))


Duration<-abs(floor(rnorm(57,mean=5,sd=2)))



sequence_driver<-function(SPSSPath,
	TPVarNames,
	TimeVarNames,
	DurationVarName,
	DecayTime=35,
	DecayRate=.95,
	) 
	{

	# decay function
	decay_func<-function(timelength,rate) {
		vect<-1 
		for (i in 2:timelength) vect<-c(vect,vect[length(vect)]*rate)
		return(vect)
		}

	# Data read
	suppressWarnings(Data<-read.spss(SPSSPath, use.value.labels = F, to.data.frame = T))
	Duration<-Data[,DurationVarName]
	
	# Grab decay
	Decay<-decay_func(DecayTime,DecayRate)
	
	MaxCode<-max(Data[,TPVarNames],na.rm=T)
	MaxTime<-max(Data[,TimeVarNames],na.rm=T)
	
	storeRspMats<-list()
	
	########
	## Start respondent loop
	for (rsp in 1:dim(Data)[1]) {
		
		rspTP<-Data[rsp,TPVarNames]
			rspTP<-rspTP[!is.na(rspTP)]
		rspTime<-Data[rsp,TimeVarNames]
			rspTime<-rspTime[!is.na(rspTime)]
		
		# Adjust all the times if the max time is not MaxTime
		if (max(rspTime)!=MaxTime) {
			adjustment<-MaxTime/max(rspTime)
			rspTime<-floor(rspTime*adjustment)
			if (max(rspTime)!=MaxTime) rspTime[length(rspTime)]<-MaxTime
			}
		
		# Make the time slider into a proportion
		rspTime<-rspTime/100

		# Initialize time decaying object for each event type
		rspMatrix<-matrix(0,nrow=Duration[rsp],ncol=MaxCode)
	
		# Map
		Map<-Duration[rsp]*rspTime
		
		# Map floored - to be able to enter into matrix, add 1
		MapFloored<-floor(Map)	
		if (min(MapFloored)==0) MapFloored[MapFloored==0]<-1
		
		
		# Store events
		for (p in 1:length(MapFloored)) rspMatrix[MapFloored[p],rspTP[p]]<-rspMatrix[MapFloored[p],rspTP[p]]+1
		
		# Apply decay to event matrix
		for (j in 1:dim(rspMatrix)[2]) {
			# check to see if there are any events of this type
				AnyNot0<-any(rspMatrix[,j]!=0)
				if (AnyNot0) {
					WhichNot0<-which(rspMatrix[,j]!=0)
					if (!all(WhichNot0%in%dim(rspMatrix)[1])) {
						storeEventDecay<-list()
						for (wn0 in 1:length(WhichNot0)) {
							initMat<-matrix(0,nrow=dim(rspMatrix)[1],ncol=1)
							initMat[WhichNot0[wn0]]<-rspMatrix[WhichNot0[wn0],j]
							if (WhichNot0[wn0]!=length(initMat)) {	
								decayloc<-2
								for (lp in 	(WhichNot0[wn0]+1):length(initMat)) {
									initMat[lp]<-initMat[lp-1]*Decay[decayloc]
									decayloc<-decayloc+1
									}
								}
							#store initMat
							storeEventDecay[[wn0]]<-initMat	
							}
							
						rspMatrix[,j]<-apply(Reduce("cbind",storeEventDecay),1,sum)
					} # close if whichnot doesn't equal the length of time
				} 
			} # close for j (the decay)
		
		# Add the duration variable to the matrix
		rspMatrix<-cbind(rspMatrix,1:Duration[rsp])
		
		#Store the respondent matrix
		storeRspMats[[rsp]]<-rspMatrix
		
		} # close the rsp loop (looping across respondents)
	
	### Reduce the data to a single data.frame
		survivalData<-as.data.frame(Reduce("rbind",storeRspMats))
		colnames(survivalData)<-c(paste("JourneyType_",1:MaxCode,sep=""),"Duration")
	
		survivalVars<-colnames(survivalData)
	
	### Estimate the survival model
		eval(parse(text=paste("Model<-glm(",survivalVars[length(survivalVars)-1],"~.,survivalData,family='binomial')",sep="")))
		
		
	}


#timeData$JourneyType_28[timeData$JourneyType_28>0 & timeData$JourneyType_28<1]<-1

