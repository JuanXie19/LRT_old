#' @rdname getClonotypeLineages
#'
#' @description This function infers the clonotype lineages structure from \code{CombinedDataSet}.
#'
#' @param data A class of \code{CombinedDataSet} object
#' @param start.clus (optional) character, indicates the starting cluster(s)
#'   from which lineages will be drawn.
#' @param end.clus (optional) character, indicates which cluster(s) will be
#'   forced to be leaf nodes in the graph.
#' @param dist.method (optional) character, specifies the method for calculating
#'   distances between clusters. Default is \code{"simple"}, see
#'   \code{\link[TrajectoryUtils]{createClusterMST}} for details.
#' @param use.median logical, whether to use the median (instead of mean) when
#'   calculating cluster centroid coordinates. Default is \code{'TRUE'}.
#' 
#' @details Given a \code{CombinedDataSet} object, this function constructs the minimum spanning tree(s) on cells sharing the same clonotype.
#'  These cells may belong to different clusters defined by scRNA-seq features, therefore a graph can be constructed with node being the center of a cluster,
#'  and weight of an edge being the distance between the nodes. MST is built based on such a graph.
#'
#' @details Once the MST is known, lineages are identified in
#'   any tree with at least two clusters. For a given tree, if there is an
#'   annotated starting cluster, every possible path out of a starting cluster
#'   and ending in a leaf that isn't another starting cluster will be returned.
#'   If no starting cluster is annotated, one will be chosen by a heuristic
#'   method, but this is not recommended.
#'
#' @import Seurat
#' @import dplyr
#' @import slingshot
#' @import dtwclust
#' @import cluster
#' @import TrajectoryUtils
#' @import progress
#' @import SingleCellExperiment
#'
#' @examples
#' TCR <-read.csv('/PATH/TO/YOUR/scTCR-seqData/',header=T)
#' load('Mice.sub.rda')
#' Combined <- getCombinedDataSet(TCR,Mice.sub)
#' Trajectory <- getClonotypeLineages(Combined,start.clus = NULL, end.clus = NULL, dist.method = 'simple', use.median = TRUE)
#'
#' @return A class of \code{TrajectoryDataSet} object
#'
#' @export 


setMethod(f = 'getClonotypeLineages',
	signature = signature('CombinedDataSet'),
	definition = function (CombinedDataSet,start.clus = NULL, end.clus = NULL, dist.method = 'simple', use.median = TRUE,...){
	
	CLONE <- unique(CombinedDataSet@seurat$cdr3)
	df <- list()
	NAMES <-vector()
	
	## set up the progress bar

	pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Inferring lineage.Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = length(CLONE),
                       complete = "=",   # Completion bar character
                       incomplete = "-", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE,    # If TRUE, clears the bar when finish
                       width = 100)      # Width of the progress bar	
	
	## median as center of cluster, and filter clusters with less than 3 cells
	for (i in 1:length(CLONE)){
		INDEX <- which(CombinedDataSet@seurat$cdr3==CLONE[i])   # 
		temp <- CombinedDataSet@seurat[,INDEX]   
		# only consider clusters with at least 3 cells
		COUNT <- as.matrix(table(temp$clusters))
		DELET <- names(which(table(temp$clusters)<3))
		INDEX1 <- which(temp$clusters%in%DELET)
		if(length(INDEX1) ==0){
			temp2 <- temp}
		else if(length(INDEX1)> 0 &length(INDEX1)< dim(temp)[2]){
			temp2 <- temp[,-INDEX1]
		}else{temp2 <- data.frame()}
 
		if(dim(temp2)[2]>3 & length(unique(temp2$clusters))>1){
			sce <- Seurat::as.SingleCellExperiment(temp2,assay='RNA')
			CLUSTERS <- SingleCellExperiment::colData(sce)$clusters
			sce <- slingshot::getLineages(sce, reducedDim = 'UMAP',use.median = use.median, 
					clusterLabels = CLUSTERS,
                     start.clus = start.clus, dist.method= dist.method)
			df[[i]] <- slingMST(sce,as.df=TRUE)
			df[[i]]$clone <-CLONE[i]
			
		}#else if(dim(temp)[2]>1){
			#print('too few cells')}
		#else{print('too few cells')} 
		
		# add the progress bar
		pb$tick()
		
	}
	
	
	df[sapply(df,is.null)] <- NULL
	for (j in 1:length(df)){
		NAMES[j] <- unique(df[[j]]$clone)
	}
	names(df) <-NAMES
	
	t <- lapply(df,function(x) x%>%group_split(Lineage)) # some MST may contain multiple lineages, here split
	t <- unlist(t,recursive=FALSE)
	t <- lapply(t,function(x) column_to_rownames(x, var = "Cluster"))	
	temp <-lapply(t, function(x) x[,1:2])

	
	Params <-list(start.clus = start.clus, end.clus = end.clus, dist.method = dist.method, use.median = use.median)
	
		
	TrajectoryData <- new('TrajectoryDataSet',lineages = temp, clonotype = NAMES,trajectoryParams = Params)
	return(TrajectoryData)
})


