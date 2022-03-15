#' get clustering input
#'
#' This function prepare the input for the \code{clustLineages} function
#' @param object A \code{TrajectoryDataSet} object, which contains a list of clonotype lineages
#'
#' @import dplyr
#' @import tidyverse
#'
#' @export
#'
#' @return A object of class \code{ClusterInput} with clonotype information as well as UMAP coordinates. 
#' 
#' @examples
#' TCR <-read.csv("/PATH/TO/YOUR/scTCR-seqData/",header=T)
#' load('Mice.sub.rda')
#' Combined <- getCombinedDataSet(TCR,Mice.sub)
#' Trajectory <- getClonotypeLineages(Combined,start.clus = NULL, end.clus = NULL, dist.method = 'simple', use.median = TRUE)
#' ClusterInput <-getClusteringInput(Trajectory)
#'
setMethod(f = 'getClusteringInput',
		signature = signature('TrajectoryDataSet'),
		definition = function(TrajectoryDataSet){
		
			df <- TrajectoryDataSet@lineages
			names(df) <- TrajectoryDataSet@clonotype
			
			t <- lapply(df,function(x) x%>%group_split(Lineage)) # some MST may contain multiple lineages, here split
			t <- unlist(t,recursive=FALSE)
			t <- lapply(t,function(x) column_to_rownames(x, var = "Cluster"))	
			
			temp <-lapply(t, function(x) x[,1:2])
			
			### try to automatic reverse the order of some dataframe
			for (i in 1:length(temp)){
				TT <- temp[[i]]
				START <- TT[1,1]       # UMAP_1 for the first cluster
				END <- TT[nrow(TT),1]  # UMAP_1 for the last cluster
				if(START <END){
					temp[[i]] <- as.data.frame(apply(temp[[i]],2,rev))
				}
			}
			
			clusterInput <- new('ClusterInput',lineages = temp)
			return(clusterInput)
			
		}
	)