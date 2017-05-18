
#' Function to remove some specific journey points from the data, and format remaining data correctly

remove_journey_points<-function(SeqDataPath, 	
	WhichEventCodesRemove=NULL,
	TerminalEventCode,
	ReplaceTerminalWithCode=NULL) {
	
	library(foreign)

	if (is.character(SeqDataPath)) Data<-read.spss(SeqDataPath,use.value.labels = F, to.data.frame = T)
	if (is.data.frame(SeqDataPath)) Data<-SeqDataPath

	whichDP<-apply(Data[,2:51],1,function(x) which(x%in%WhichEventCodesRemove))
	DataDP<-Data
	for (i in 1:dim(Data)[1]) DataDP[i,whichDP[[i]]+1]<-NA

	whichValid<-apply(DataDP,1,function(x) which(!is.na(x)))

	DataDPFinal<-matrix(,ncol=51,nrow=dim(Data)[1])
	for (i in 1:dim(Data)[1]) {
		DataDPFinal[i,1:length(whichValid[[i]])]<-as.matrix(DataDP[i,whichValid[[i]]])
		DataDPFinal[i,which(DataDPFinal[i,]%in%TerminalEventCode)]<-ReplaceTerminalWithCode
		}
	colnames(DataDPFinal)<-colnames(Data)
	return(as.data.frame(DataDPFinal))

}



#' Apply a set of coefficients from the wrapper to a new dataset
#' 

apply_dfacoefs<-function(
	SeqDataPath,
	priorP2PwrapObj,		# Object from p2p_wrap
	whichCluster=NULL
	)
	{		
	
	if (is.null(whichCluster)) stop("Specify from which cluster you want to extract the DFA coefficients.")
	
	# Retrieve stuff from the priorP2PwrapObj
	
		CapSeqLength<-priorP2PwrapObj$CapSeqLength
		SeqMinLength<-priorP2PwrapObj$SeqMinLength
		SingleState<-priorP2PwrapObj$SingleState	
		DFACoefs<-eval(parse(text=paste("gap_m4_c8_15_ssF$DFACoefs$cluster",whichCluster,sep="")))

	#######################
	## Data manipulation ##
	#######################
	 	
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

	# Non sequence group
	NonSequenceGroup<-list()
		
	# Load data
	message("Loading in data")
	suppressWarnings(Data<-read.spss(SeqDataPath, use.value.labels = F, to.data.frame = T))
	
	# retrieve event type labels
	if (is.null(tpLabels)) tpLabels<-names(attr(Data[,2],"value.labels"))[order(as.numeric(attr(Data[,2],"value.labels")))]
	
	
		# Remove rows with only missing data
		if (length(which(apply(Data[,-c(1,2)],1,function(x) all(is.na(x)))))>0) Data<-Data[-which(apply(Data[,-c(1,2)],1,function(x) all(is.na(x)))),]
	
		# If there's been a cap placed on sequence length, remove sequences greater than that cap
		if (!is.null(CapSeqLength)) {
			NonSequenceGroup[[length(NonSequenceGroup)+1]]<-Data[which(apply(Data[,-1],1,function(x) table(!is.na(x))[2]<SeqMinLength)),1]
			Data<-Data[which(apply(Data[,-1],1,function(x) table(!is.na(x))[2]<=CapSeqLength & table(!is.na(x))[2]>=SeqMinLength)),]
			Data<-Data[,c(1:(CapSeqLength+1))]
			}
	
		# Record the ids of everyone who remains
		IDs<-Data[,1]
	
		# Retrieve the max code for the touchpoints
		maxtpCode=max(Data[,-1],na.rm=T)
		TerminalEventLabel<-tpLabels[maxtpCode]
		
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
			NonSequenceGroup[[length(NonSequenceGroup)+1]]<-Data[Which1,1]
			# Truncate sequences to those that are multi-step and those that are single step
			SeqDataSingle<-SeqData[Which1,] # can return this later
			SeqData<-SeqData[-(Which1),]
			IDs<-IDs[-(Which1)]
			}
	
	
	########################
	## Apply Coefficients ##
	########################
	
			# Retrieves indvidual-level measurements on the sequences
			Length<-seqlength(SeqData)
			FirstEvent<-CreateDummy(SeqData[,1],"FirstEvent")
			LastEvent<-matrix(,ncol=1,nrow=length(Length))
				for(i in 1:length(Length)) LastEvent[i]<-SeqData[i,Length[i]]
				if (length(table(LastEvent))==1) for(i in 1:length(Length)) LastEvent[i]<-SeqData[i,Length[i]-1]
				LastEvent<-CreateDummy(LastEvent,"LastEvent")
			CountOfEachEvent<-seqistatd(SeqData)
				colnames(CountOfEachEvent)<-paste("CountOfEventType",colnames(CountOfEachEvent),sep="")
				
			# Remove highly correlated variables at > .6
				DFAData<-cbind(FirstEvent,LastEvent,CountOfEachEvent)
			
		
			# Apply	
			WhichCols<-match(rownames(DFACoefs),colnames(DFAData))
			AlgoData<-DFAData[,WhichCols[2:length(WhichCols)]]	
			
			apply_coef<-function(tData,Coefs) Coefs[1]+sum(tData*Coefs[2:length(Coefs)])
				
			storeScores<-list()
			for (hmc in 1:whichCluster) storeScores[[hmc]]<-apply(AlgoData,1,apply_coef,DFACoefs[,hmc])
			storeScores<-Reduce("cbind",storeScores)
			
			Assigned<-rbind(cbind(IDs,max.col(storeScores)),cbind(unlist(NonSequenceGroup),99))
			colnames(Assigned)<-c("IDs","currSol")
			return(Assigned)
			
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
			SeqMinLength = 3,
			tpLabels = NULL,
			SingleState= FALSE,			
			HowManyClusters=5:10,
			SeriesTitle=NULL,		
			costMatrix=NULL,
			OtherProfileBinariesPath=NULL, # A path to an SPSS file that contains profile variable binaries, first column should be IDs.
			ExcelOutPath=NULL,
			ConvertOutToSPSS=TRUE,		# Convert output to an SPSS conformable means table?
			#Group=NULL,				# Either Null, or a two slot list, one with variable name, one with variable value to select on, e.g. list("DV_Category_PipeIn",2). Creates sequences within the group. 
			skewThreshold=.3,			# Threshold for skews to qualify as a connection in the journey map
			CustomCluster=NULL,			# Nx2 dataframe, first column containing IDs, second column the labels. 
			ReportBtwWithinPlot=TRUE,
			CustomClusterIsJourneyGroup=TRUE # is the custom cluster the journey group (allows us to generate full tables with other subgroups if FALSE)
			
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
	
	# Assign Sequence moments to a given cluster.
	assigner<-function(MeansTable,Thresh) 
		{
		rownames(MeansTable)<-lapply(lapply(strsplit(rownames(MeansTable)," ",""),function(x) x[-1]),paste, collapse=" ")
		store<-list(); storeMeans<-list()
		
		tempMeansTable<-MeansTable[,1:(dim(MeansTable)[2]-1)]

		for (mtm in 1:dim(tempMeansTable)[2]) {
			store[[mtm]]<-rownames(MeansTable[which(tempMeansTable[,mtm]>Thresh),])	
			storeMeans[[mtm]]<-MeansTable[which(tempMeansTable[,mtm]>Thresh),mtm]
			}
		return(list(Arcs=store,ArcSizes=storeMeans))
	}
	
	# Function to create a simple baliaviz visualization from a converted edgelist output
	gen_baliaviz_fromconvertedgelist<-function(theConvertedEdgeList) 
		{
				
				LINKS<-as.data.frame(theConvertedEdgeList$Numeric)
					colnames(LINKS)<-c("source","target","value")
				
				NODES<-as.data.frame(theConvertedEdgeList$Character)
					colnames(NODES)<-c("name","group","size")
					NODES$group<-as.numeric(as.character(NODES$group))
					NODES$size<-as.numeric(as.character(NODES$size))
				
				LINKS$source<-as.numeric(as.character(LINKS$source))
				LINKS$target<-as.numeric(as.character(LINKS$target))
				LINKS$value<-as.numeric(as.character(LINKS$value))
				NODES$group<-as.numeric(as.character(NODES$group))
				NODES$size<-as.numeric(as.character(NODES$size))
					
				if (!is.null(theConvertedEdgeList$Title)) return(baliaviz::forceNetwork(Nodes = NODES, Links = LINKS,title=theConvertedEdgeList$Title))
				if (is.null(theConvertedEdgeList$Title)) return(baliaviz::forceNetwork(Nodes = NODES, Links = LINKS))

				}
	
	# Non sequence group
	NonSequenceGroup<-list()
	NoJourneyGroup<-list()
		
	# Load data
	message("Loading in data")
	if (is.character(SeqDataPath)) 
		{
		suppressWarnings(Data<-read.spss(SeqDataPath, use.value.labels = F, to.data.frame = T))
		# retrieve event type labels
		if (is.null(tpLabels)) tpLabels<-names(attr(Data[,2],"value.labels"))[order(as.numeric(attr(Data[,2],"value.labels")))]
		}
	
	if (is.data.frame(SeqDataPath)) Data<-SeqDataPath
	
		# Remove rows with only missing data
		if (length(which(apply(Data[,-c(1,2)],1,function(x) all(is.na(x)))))>0 && CustomClusterIsJourneyGroup) {
			NoJourneyGroup[[1]]<-Data[which(apply(Data[,-c(1,2)],1,function(x) all(is.na(x)))),1]
			Data<-Data[-which(apply(Data[,-c(1,2)],1,function(x) all(is.na(x)))),]
			}
	
		# If there's been a cap placed on sequence length, remove sequences greater than that cap
		if (!is.null(SeqMinLength)) {
			NonSequenceGroup[[length(NonSequenceGroup)+1]]<-Data[which(apply(Data[,-1],1,function(x) table(!is.na(x))[2]<SeqMinLength)),1]
			Data<-Data[which(apply(Data[,-1],1,function(x) table(!is.na(x))[2]<=CapSeqLength & table(!is.na(x))[2]>=SeqMinLength)),]
			Data<-Data[,c(1:(CapSeqLength+1))]
			}
	
		# If customcluster has been specified and there are 98s or 99s as part of that solution, remove those codes
		#if (!is.null(CustomCluster) && any(CustomCluster[,2]%in%c(98,99))) {
		#	NonSequenceGroup[[length(NonSequenceGroup)+1]]<-Data[which(CustomCluster[,2]%in%c(98,99)),1]
		#	Data<-Data[-which(CustomCluster[,2]%in%c(98,99)),]
		#	CustomCluster<-CustomCluster[-which(CustomCluster[,2]%in%c(98,99)),]
		#	}
	
		# Record the ids of everyone who remains
		IDs<-Data[,1]
	
		# Retrieve the max code for the touchpoints
		maxtpCode=max(Data[,-1],na.rm=T)
		TerminalEventLabel<-tpLabels[maxtpCode]
		
	# Start building sequences
		if (!is.null(CapSeqLength)) SeqData<-seqdef(Data,var=colnames(Data)[2:(CapSeqLength+1)])
		if (is.null(CapSeqLength)) SeqData<-seqdef(Data[,-1])
	
	# Check to see if single state is turned on, and if so, revise sequence data
		if (SingleState) SeqData<-seqdss(SeqData)

	# Check to see if there are any single-step sequences.  If so, remove them from the sequences
		Any1<-any(seqlength(SeqData)%in%1)	
		if (Any1 && CustomClusterIsJourneyGroup) {
			message("Removing single-step sequences")
			Which1<-which(seqlength(SeqData)%in%1)	
			NonSequenceGroup[[length(NonSequenceGroup)+1]]<-Data[Which1,1]
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

	if (is.null(CustomCluster)) {
	
		message("Clustering paths")
		clusterward <- agnes(Distances, diss = TRUE, method = "ward")

		################################################
		### Create within vs. between distance plot. ###
		################################################
		if (ReportBtwWithinPlot) {
		
			message("Generating Between vs. Within Cluster Distances Plot")
			storeClusterDistances<-list()
		
				# Create plot metrics
				DistWithin<-matrix(,ncol=1,nrow=length(1:50))
				DisBtw<-matrix(,ncol=1,nrow=length(1:50))
				for (i in 2:50) {
					curCut<-cutree(clusterward, k=i)
			
					# Avg distance within
					tempStore<-matrix(,ncol=1,nrow=i)
					for (kk in 1:i) {
						tempDist<-Distances[which(curCut==kk),which(curCut==kk)]
						tempStore[kk]<-mean(tempDist[lower.tri(tempDist)])
						}
						DistWithin[i]<-mean(tempStore)
			
					# Distances between 
					tempGroupDistances<-matrix(0,i,i)
					for (kk1 in 1:i) {
						for (kk2 in 1:i) {
							if (kk1!=kk2) tempGroupDistances[kk1,kk2]<-mean(Distances[which(curCut==kk1),which(curCut==kk2)])
							}
						}
						# Store the cluster distance matrix if its specified in HowManyClusters
						if (any(i==HowManyClusters)) storeClusterDistances[[paste("cluster",i,"_dists",sep="")]]<-tempGroupDistances
						DisBtw[i]<-mean(tempGroupDistances[lower.tri(tempGroupDistances)])
					}
		
				# Look at rate of change of btw vs within cluster			
				DiffBtwWith<-diff(DisBtw-DistWithin)
				storeZ<-matrix(,ncol=1,nrow=length(DiffBtwWith))
				for (ll in 2:length(DiffBtwWith)) {
					storeZ[ll]<-(DiffBtwWith[ll]-mean(DiffBtwWith[1:ll],na.rm=T))/sd(DiffBtwWith[1:ll],na.rm=T)
					}
			
				library(ggplot2)
				plotData<-as.data.frame(cbind(DiffBtwWith,1:length(DiffBtwWith)))
					colnames(plotData)<-c("Distances","Clusters")
				plotter<-ggplot(plotData)
				microClusterPlot<-plotter+geom_line(aes(x=Clusters,y=Distances))+
					stat_smooth(aes(x=Clusters,y=Distances),se=F)+
					ggtitle("Difference Between vs. Within Cluster Distances")
			
			}	
			
		################################################
		################################################
		
		# Use HowManyClusters to control the search space

			if (is.null(HowManyClusters)) stop("Please specify the number(s) of cluster(s) you want to create")
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
				DFAData<-DFAData[,-which(colnames(DFAData)%in%apply(corGraph,1,function(x) x[1]))]
			
		
			clusterModelAccuracy<-list()
			clusterModelCoefs<-list()
			clusterSubSeqs<-list()
			storeClusterDistances<-list()
			for (hmc in 1:length(HowManyClusters)) {
				
				message(paste("Evaluating DFA for cluster",HowManyClusters[hmc],sep=""))
					
								
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
			
				# restore the cluster-level distances using DFA'ed groupings
				
				tempGroupDistances<-matrix(0,HowManyClusters[hmc],HowManyClusters[hmc])
				for (kk1 in 1:HowManyClusters[hmc]) {
					for (kk2 in 1:HowManyClusters[hmc]) {
						if (kk1!=kk2) eval(parse(text=paste("tempGroupDistances[kk1,kk2]<-mean(Distances[which(DFA_cluster",HowManyClusters[hmc],"==kk1),which(DFA_cluster",HowManyClusters[hmc],"==kk2)])",sep="")))
						}
					}
					# Store the cluster distance matrix if its specified in HowManyClusters
					storeClusterDistances[[paste("DFA_cluster",HowManyClusters[hmc],"_dists",sep="")]]<-tempGroupDistances
			
				}
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
				
		if (!DFA && is.null(CustomCluster)) ClusterNames<-paste("cluster",HowManyClusters,sep="")
		if (DFA && is.null(CustomCluster)) ClusterNames<-paste("DFA_cluster",HowManyClusters,sep="")
		if (!is.null(CustomCluster)) {
			ClusterNames<-"CustomCluster"
			CustomCluster<-CustomCluster[match(IDs,CustomCluster[,1]),2]
			}
		
		if (ConvertOutToSPSS) SPSSOut<-list()
		
		BindedTablesStore<-list()
		
		storeProto<-list()
		storeClusterSolutions<-list()
		masterStoreGraphs<-list()
		masterStoreGraphTables<-list()
		storeBaseSizes<-list()
		storeCompositeProtos<-list()
		storeCompositeProtosBaliaviz<-list()
		storeCompositeProtoDistances<-list()
		storeSimplyCompositeProtoDistances<-list()
		
		
		storeStandCPSDists<-list()
		storeStandClassAssign<-list()
		storeCPSClassAssignProp<-list()
		
		for (cN in ClusterNames) {
			
			## Compute means for profiling
			
			currSol<-eval(parse(text=cN))
			if (class(currSol)=="integer") currSol<-as.factor(currSol)
			maxCurrSol<-max(as.numeric(as.character(currSol)))
			if (length(unlist(NonSequenceGroup))==0) storeClusterSolutions[[cN]]<-cbind(IDs,currSol)
			if (length(unlist(NonSequenceGroup))>0 & length(unlist(NoJourneyGroup))==0) storeClusterSolutions[[cN]]<-rbind(cbind(IDs,currSol),cbind(unlist(NonSequenceGroup),99))
			if (length(unlist(NonSequenceGroup))==0 & length(unlist(NoJourneyGroup))>0) storeClusterSolutions[[cN]]<-rbind(cbind(IDs,currSol),cbind(unlist(NoJourneyGroup),98))
			if (length(unlist(NonSequenceGroup))>0 & length(unlist(NoJourneyGroup))>0) storeClusterSolutions[[cN]]<-rbind(cbind(IDs,currSol),cbind(unlist(NoJourneyGroup),98),cbind(unlist(NonSequenceGroup),99))
			
			## Identify PROTOTYPICAL sequences
			#if (is.null(CustomCluster)) {
				DistanceFromCenter<-disscenter(Distances,currSol)
				#DistanceFromCenter[which(currSol==1)]
			
				tempStoreProto<-list()
				for (rNg in 1:maxCurrSol) {
					#currProto<-unique(head(SeqData[currSol%in%rNg,][order(DistanceFromCenter[currSol%in%rNg]),],20))
					currProto<-SeqData[currSol%in%rNg,][order(DistanceFromCenter[currSol%in%rNg]),]
					currProto<-apply(currProto,2,as.character)
					storeCurrProto<-list()
					if (!is.null(dim(currProto))) for (prt in 1:dim(currProto)[1]) storeCurrProto[[prt]] <-suppressWarnings(paste(tpLabels[as.numeric(currProto[prt,])][!is.na(tpLabels[as.numeric(currProto[prt,])])],collapse=">"))
					if (is.null(dim(currProto))) storeCurrProto[[1]] <-suppressWarnings(paste(tpLabels[as.numeric(currProto)][!is.na(tpLabels[as.numeric(currProto)])],collapse=">"))
					tempStoreProto[[rNg]]<-unlist(storeCurrProto)
					}
				
				storeProto[[paste("cluster",cN,sep="")]]<-tempStoreProto
				#}
			
			# Map each prototype as a (temporally) directed graph - composite prototypes
			AvgProtoSeqLength<-unlist(lapply(lapply(lapply(tempStoreProto,function(x) strsplit(x,">","")),function(x) lapply(x,length)),function(x) floor(mean(unlist(x)))))
			ProtoLengths<-lapply(lapply(lapply(lapply(tempStoreProto,function(x) strsplit(x,">","")),function(x) lapply(x,length)),unlist),table)
			max2<-unlist(lapply(ProtoLengths,function(x) as.numeric(names(x)[x%in%max(x[-which(x%in%max(x))])])))
			longest<-unlist(ProtoLengths,function(x) as.numeric(names(x)[length(x)]))
			CompositeProto<-list()
			for (clu in 1:maxCurrSol) {
				pieceCompProto<-list()
				protoSplit<-lapply(tempStoreProto[[clu]],function(x) strsplit(x,">",""))
				#for (apsl in 1:AvgProtoSeqLength[clu]) {
				for (apsl in 1:max2[clu]) {
				#for (apsl in 1:longest[clu]) {
					tempTab<-table(unlist(lapply(protoSplit,function(x)x[[1]][apsl])))
					if (any(names(tempTab)%in%TerminalEventLabel)) tempTab<-tempTab[-which(names(tempTab)%in%TerminalEventLabel)]
					namesTempTab<-names(tempTab)
					maxTempTab<-max(tempTab)
					tempNamesStore<-namesTempTab[which(tempTab%in%maxTempTab)]
					#if (any(tempNamesStore%in%TerminalEventLabel)) tempNamesStore<-tempNamesStore[-which(tempNamesStore%in%TerminalEventLabel)]
					pieceCompProto[[apsl]]<-tempNamesStore
					}
					
				# Put the purchase event as the last event, and remove it from any other slot
				#pieceCompProto[[length(pieceCompProto)+1]]<-TerminalEventLabel
				if (any(unlist(lapply(pieceCompProto,length))==0)) pieceCompProto<-pieceCompProto[-which(unlist(lapply(pieceCompProto,length))==0)]
								
				
				#clean it out
				#for (cc in 1:(length(pieceCompProto)-1)) {
				#	if (any(pieceCompProto[[cc]]%in%pieceCompProto[[cc+1]])) {
				#		finder<-which(pieceCompProto[[cc+1]]%in%pieceCompProto[[cc]])
				#		pieceCompProto[[cc+1]]<-pieceCompProto[[cc+1]][-finder]
				#		}
				#	}
				
				#if (any(unlist(lapply(pieceCompProto,length))==0)) pieceCompProto<-pieceCompProto[-which(unlist(lapply(pieceCompProto,length))==0)]
			
				#piece edge list together
				adder<-1
				pieceCompProtoEdge<-list()
				simplyCompProtoEdge<-list()
				for (cc in 1:(length(pieceCompProto)-1)) {
					for (pcp in 1:length(pieceCompProto[[cc]])) {
							pieceCompProtoEdge[[adder]]<-paste(pieceCompProto[[cc]][pcp],pieceCompProto[[cc+1]],sep=">")
							simplyCompProtoEdge[[adder]]<-paste(paste("(",1:length(pieceCompProto[[cc]][pcp]),") ",pieceCompProto[[cc]],sep=""),paste(" (",1:length(pieceCompProto[[cc+1]]),") ",pieceCompProto[[cc+1]],sep=""),sep=" >>> ",collapse="")
							adder<-adder+1
							}
						}
					
				CompositeProto[[clu]]<-Reduce("rbind",strsplit(unlist(pieceCompProtoEdge),">",""))
				
				# Add in the terminal event
				if (length(CompositeProto[[clu]])>2) CompositeProto[[clu]]<-rbind(CompositeProto[[clu]],cbind(unique(matrix(CompositeProto[[clu]])),TerminalEventLabel))
				if (length(CompositeProto[[clu]])==2) CompositeProto[[clu]]<-rbind(CompositeProto[[clu]],cbind(unique(CompositeProto[[clu]]),TerminalEventLabel))
				
				# Remove duplicated edges
				CompositeProto[[clu]]<-Reduce("rbind",strsplit(unique(apply(CompositeProto[[clu]],1,paste,collapse=">")),">",""))
				
				if (length(CompositeProto[[clu]])==2) CompositeProto[[clu]]<-matrix(CompositeProto[[clu]],ncol=2,byrow=T)
				}
			
			# Convert the character edgelist for mapping
			storeCompositeProtos[[cN]]<-lapply(lapply(CompositeProto,balialab::convert_edgelist),function(x) list(Numeric=x$Numeric-1,Character=x$Character))
			

			
						
			# Identify distances between micro-clusters
			tempGroupDistances<-matrix(0,maxCurrSol,maxCurrSol)
			for (kk1 in 1:maxCurrSol) {
				for (kk2 in 1:maxCurrSol) {
					if (kk1!=kk2) tempGroupDistances[kk1,kk2]<-mean(Distances[which(currSol==kk1),which(currSol==kk2)])
					}
				}
			
			meanDist<-mean(tempGroupDistances[lower.tri(tempGroupDistances)])
			sdDist<-sd(tempGroupDistances[lower.tri(tempGroupDistances)])
			standDist<-round((tempGroupDistances-meanDist)/sdDist,digits=2)
			diag(standDist)<-""
			standDist[upper.tri(standDist)]<-""
			standDist<-as.data.frame(standDist)
			
			CPSNames<-paste("cps",1:maxCurrSol,sep="")
			colnames(standDist)<-CPSNames
			rownames(standDist)<-CPSNames
			
			
			# store it
			storeStandCPSDists[[cN]]<-standDist
			
			
			# Explore the sequence class membership	(range of classes is automatically selected, around a 1/5 of the number of sequences
			if (is.null(CustomCluster)) {
				RangeClassesExplore<-ceiling(maxCurrSol/5):(ceiling(maxCurrSol/5)+ceiling(maxCurrSol/5))
				if (any(RangeClassesExplore==1)) RangeClassesExplore<-RangeClassesExplore[-1] # remove 1 class to explore, cuz can't do it. 
				tempStoreClassAssign<-matrix(,ncol=length(RangeClassesExplore),nrow=length(currSol))
				tempTempStoreCPSClassAssignProp<-list()
				for (rce in RangeClassesExplore) {
					currClassMembership<-cutree(agnes(as.dist(tempGroupDistances), diss = TRUE, method = "ward"),k=rce)
					tempStoreCPSClassAssignProp<-list()
					for (tg in 1:rce) {
						tempStoreClassAssign[which(currSol%in%which(currClassMembership%in%tg)),which(RangeClassesExplore%in%rce)]<-tg
						tempTab<-table(currSol[currSol%in%which(currClassMembership%in%tg)])
						tempStoreCPSClassAssignProp[[tg]]<-round(tempTab/sum(tempTab),digits=2)
						}
					colnames(tempStoreClassAssign)<-paste("numClasses",RangeClassesExplore,sep="")		
					toStore<-Reduce("rbind",tempStoreCPSClassAssignProp)
					rownames(toStore)<-paste("class",1:rce,sep="")
					colnames(toStore)<-CPSNames
					tempTempStoreCPSClassAssignProp[[rce]]<-toStore
					}
			
				# store it
				storeCPSClassAssignProp[[cN]]<-tempTempStoreCPSClassAssignProp
				storeStandClassAssign[[cN]]<-cbind(IDs,tempStoreClassAssign)
			}
				
			
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
				TotalPercent<-table(factor(matrix(as.matrix(currClusSeqs)),levels=1:maxtpCode))/sum(table(factor(matrix(as.matrix(currClusSeqs)),levels=1:maxtpCode)))
				
				# Transition rates
				TransitionRate<-matrix(t(seqtrate(currClusSeqs)))
				rownames(TransitionRate)<-matrix(t(seqetm(currClusSeqs)))
				
				# Store
				FirstList[[rNg]]<-FirstEventTable
				LastList[[rNg]]<-LastEventTable
				TotalList[[rNg]]<-TotalPercent
				TransitionList[[rNg]]<-TransitionRate
				
					##################################################################################################################
					## Go back to composite prototype tables to fill in those tables with profiling information about the sequences ##
					##################################################################################################################
					if (rNg<=maxCurrSol) {
					
						tempNet<-storeCompositeProtos[[cN]][[rNg]]$Numeric+1
						tempNodes<-storeCompositeProtos[[cN]][[rNg]]$Character
						tempNet[,1]<-match(tempNodes,tpLabels)[tempNet[,1]]
						tempNet[,2]<-match(tempNodes,tpLabels)[tempNet[,2]]
						tempNet<-apply(tempNet,1,paste,collapse=">")
						tempLinkSizes<-TransitionRate[match(tempNet,rownames(TransitionRate))]
						if (any(is.na(tempLinkSizes))) tempLinkSizes[is.na(tempLinkSizes)]<-.5
						tempNodeSizes<-TotalPercent[match(tempNodes,tpLabels)]
						
						# Remove an edge if the transition probability is 0
						if (any(tempLinkSizes%in%0)) {
							whichtoremove<-which(tempLinkSizes%in%0)
							tempLinkSizes<-tempLinkSizes[-whichtoremove]
							storeCompositeProtos[[cN]][[rNg]]$Numeric<-storeCompositeProtos[[cN]][[rNg]]$Numeric[-whichtoremove,]
							}
					
						storeCompositeProtos[[cN]][[rNg]]$Numeric<- cbind(storeCompositeProtos[[cN]][[rNg]]$Numeric,2-(2-10)*((tempLinkSizes-min(tempLinkSizes,na.rm=T))/(max(tempLinkSizes,na.rm=T)-min(tempLinkSizes,na.rm=T))))
						storeCompositeProtos[[cN]][[rNg]]$Character<-cbind(tempNodes,c(5,rep(1,length(tempNodes)-2),3),2-(2-10)*((tempNodeSizes-min(tempNodeSizes))/(max(tempNodeSizes)-min(tempNodeSizes))))
						if (!is.null(SeriesTitle)) storeCompositeProtos[[cN]][[rNg]]$Title<-paste(SeriesTitle,"_",cN,"_gr",rNg,sep="")
					}
				
				}
					
					# Create the baliaviz visual
				#	storeCompositeProtosBaliaviz[[cN]]<-lapply(storeCompositeProtos[[cN]],gen_baliaviz_fromconvertedgelist)	
					
					
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
			if (is.null(CustomCluster)) Return[["Prototypes"]]<-storeProto
			Return[["TPLabels"]]<-tpLabels
			Return[["Solutions"]]<-storeClusterSolutions
			Return[["BaseSizes"]]<-storeBaseSizes
			if (ConvertOutToSPSS) Return[["SPSSMainTable"]]<-SPSSOut
			if (DFA) {
				Return[["DFACoefs"]]<-clusterModelCoefs
				Return[["DFAModelAccuracy"]]<-clusterModelAccuracy
				Return[["SubSeqsForDFA"]]<-clusterSubSeqs
				}
			if (ReportBtwWithinPlot) Return[["microClusterPlot"]]<-microClusterPlot
			Return[["CompositeProtos"]]<-storeCompositeProtos
			#Return[["CompositeProtosBaliaviz"]]<-storeCompositeProtosBaliaviz
			Return[["CPSDistances"]]<-storeStandCPSDists
			if (is.null(CustomCluster)) Return[["ClassAssigned"]]<-storeStandClassAssign
			if (is.null(CustomCluster)) Return[["ClassAssignProps"]]<-storeCPSClassAssignProp
			Return[["CapSeqLength"]]<-CapSeqLength
			Return[["SeqMinLength"]]<-SeqMinLength
			Return[["SingleState"]]<-SingleState
		return(Return)		
	}



#' Survival based drivers for sequence data


sequence_driver<-function(SPSSPath,
	TPVarNames,
	TimeVarNames,
	DurationVarName,
	#DecayTime=35,
	DecayRate=.95,
	tpLabels=NULL,
	TermEventNumCode=NULL,
	JourneySubgroups=NULL,				# this should be a specific Solutions slot from the p2p_wrap
	ExcelOutPath=NULL
	) 
	{

	if (is.null(TermEventNumCode)) stop("Please specify the numeric code for the terminal event")
	
	# decay function
	decay_func<-function(timelength,rate) {
		vect<-1 
		for (i in 2:timelength) vect<-c(vect,vect[length(vect)]*rate)
		return(vect)
		}
		
	  
	library(foreign)
	
	# Data read
	suppressWarnings(Data<-read.spss(SPSSPath, use.value.labels = F, to.data.frame = T))

	if (is.null(tpLabels)) tpLabels<-names(attr(Data[,TPVarNames][,1],"value.labels"))[order(as.numeric(attr(Data[,TPVarNames][,1],"value.labels")))]

	# Remove cases that are classified as having short journey groups
	if (!is.null(JourneySubgroups)) {
		GroupCodesAtData<-JourneySubgroups[match(Data[,1],JourneySubgroups[,1]),2]
		if (any(JourneySubgroups[,2]%in%99)) JourneySubgroups<-JourneySubgroups[-which(JourneySubgroups[,2]%in%99),]
		Data<-Data[-which(GroupCodesAtData%in%99 | is.na(GroupCodesAtData)),]
		}
	
	
	Duration<-Data[,DurationVarName]
	
	# Grab decay
	Decay<-decay_func(max(Duration),DecayRate)
	
	MaxCode<-max(Data[,TPVarNames],na.rm=T)
	MaxTime<-max(Data[,TimeVarNames],na.rm=T)
	
	storeRspMats<-list()
	storeSerials<-list()
	
	message("Expanding survival data with decay")
	###########################
	## Start respondent loop ##
	for (rsp in 1:dim(Data)[1]) {
		
		storeSerials[[rsp]]<-Data[rsp,1]
		
		rspTP<-Data[rsp,TPVarNames]
			rspTP<-rspTP[!is.na(rspTP)]
		rspTime<-Data[rsp,TimeVarNames]
			rspTime<-rspTime[!is.na(rspTime)]
		
		Serials<-Data[rsp,1]
		
		if (all(rspTime==0)) rspTime[length(rspTime)]<-MaxTime
		
		# Adjust all the times if the max time is not MaxTime
		if (max(rspTime)!=MaxTime) {
			adjustment<-MaxTime/max(rspTime)
			rspTime<-floor(rspTime*adjustment)
			if (max(rspTime)!=MaxTime) rspTime[length(rspTime)]<-MaxTime
			}
		
		# Make the time slider into a proportion
		rspTime<-rspTime/100

		# Initialize time decaying object for each event type
		rspMatrix<-matrix(0,nrow=Duration[rsp]+1,ncol=MaxCode)
	
		# Map
		Map<-Duration[rsp]*rspTime
		
		# Map floored - to be able to enter into matrix, add 1
		MapFloored<-floor(Map)	
		if (min(MapFloored)==0) MapFloored[MapFloored==0]<-1
		
		
		# Store events
		for (p in 1:length(MapFloored)) {
			if (!(rspTP[p]%in%TermEventNumCode)) rspMatrix[MapFloored[p],rspTP[p]]<-rspMatrix[MapFloored[p],rspTP[p]]+1
			if (rspTP[p]%in%TermEventNumCode) rspMatrix[dim(rspMatrix)[1],rspTP[p]]<-1
			}
	
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
		rspMatrix<-cbind(rspMatrix,0:Duration[rsp])
		
		#Store the respondent matrix
		storeRspMats[[rsp]]<-rspMatrix
		
		} # close the rsp loop (looping across respondents)
	
	### Reduce the data to a single data.frame
		StoreResults<-list()
		
		message("Running total level survival impact")
		survivalData<-as.data.frame(Reduce("rbind",storeRspMats))
		whichLowVar<-which(apply(survivalData,2,var,na.rm=T)<.009)
		if (length(whichLowVar)>0) survivalData[,whichLowVar]<-0
		
		colnames(survivalData)<-c(paste("JourneyType_",1:MaxCode,sep=""),"Duration")
		survivalVars<-colnames(survivalData)
	
		VarsShouldModel<-colnames(survivalData)
			VarsShouldModel<-VarsShouldModel[-which(VarsShouldModel%in%c("Duration",paste("JourneyType_",TermEventNumCode,sep="")))]
		
		eval(parse(text=paste("Model<-glm(",survivalVars[length(survivalVars)-1],"~.,survivalData,family='binomial')",sep="")))
		tempResults<-coef(summary(Model))[,c(1,4)]
			tempResults<-tempResults[match(VarsShouldModel,rownames(tempResults)),]
			tempResults<-cbind(tempResults[,1],abs(tempResults[,1]),tempResults[,2],round(100*abs(tempResults[,1])/mean(abs(tempResults[,1]),na.rm=T),digits=0))

		colnames(tempResults)<-paste(paste("Gr_Total",sep=""),c("Coef","absCoef","sig","Indexed"),sep="")
		rownames(tempResults)<-tpLabels[-TermEventNumCode]
		StoreResults[["TotalImpact"]]<-tempResults
		
		if (!is.null(JourneySubgroups)) {
			message("Running journey/sub group level impact")
			JourneySubgroups<-as.data.frame(JourneySubgroups)	
			HowManyGroups<-max(JourneySubgroups$currSol)
			GroupCodesAtData<-JourneySubgroups[match(Data[,1],JourneySubgroups[,1]),2]
			for (gr in 1:HowManyGroups) {
				message(paste("Group ", gr, " of ",HowManyGroups, sep=""))
				gr_filter<-which(GroupCodesAtData%in%gr)
				survivalData<-as.data.frame(Reduce("rbind",storeRspMats[gr_filter]))
				colnames(survivalData)<-c(paste("JourneyType_",1:MaxCode,sep=""),"Duration")
				whichLowVar<-which(apply(survivalData,2,var,na.rm=T)<.009)
				if (length(whichLowVar)>0) survivalData[,whichLowVar]<-0
		
				
				survivalVars<-colnames(survivalData)
				eval(parse(text=paste("Model<-glm(",survivalVars[length(survivalVars)-1],"~.,survivalData,family='binomial')",sep="")))
				tempResults<-coef(summary(Model))[,c(1,4)]
					tempResults<-tempResults[match(VarsShouldModel,rownames(tempResults)),]
					tempResults<-cbind(tempResults[,1],abs(tempResults[,1]),tempResults[,2],round(100*abs(tempResults[,1])/mean(abs(tempResults[,1]),na.rm=T),digits=0))
				colnames(tempResults)<-paste(paste("Gr_",gr,sep=""),c("Coef","absCoef","sig","Indexed"),sep="")
				rownames(tempResults)<-tpLabels[-TermEventNumCode]
				StoreResults[[paste("JGImpact",gr,sep="")]]<-tempResults

				}
			}
			
		Return<-as.data.frame(Reduce("cbind",StoreResults))
		#colnames(Return)<-matrix(matrix(c(paste(colnames(Return)[1],names(StoreResults),sep="_"),paste(colnames(Return)[2],names(StoreResults),sep="_")),nrow=2,ncol=length(StoreResults),byrow=T))
	
	if (!is.null(ExcelOutPath)) {
		library(WriteXLS)
		WriteXLS("Return",ExcelOutPath,row.names=T,col.names=T)
		}
		
	return(Return)
	
	}



