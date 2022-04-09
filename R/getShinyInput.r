#' prepare Shiny input
#'
#' This function prepares the input for the two shiny apps
#' @param object A \code{CombinedDataSet} object 
#'
#' @details This function takes the \code{CombinedDataSet}object as the input and extract the UMAP coordinates for each clonotype, which will be used as input of the two shiny apps.
#'
#' @import Seurat
#' @import dplyr
#' @import tidyverse
#'
#' @export
#'
#' @return A data frame with cell clonotype and UMAP coordinates . 
#' 
#' @examples
#' TCR <-read.csv("/PATH/TO/YOUR/scTCR-seqData/",header=T)
#' load('Mice.sub.rda')
#' Combined <- getCombinedDataSet(TCR,Mice.sub)
#' Trajectory <- getClonotyoeLineages(Combined,start.clus = NULL, end.clus = NULL, dist.method = 'simple', use.median = TRUE)
#' clustInput <- getClusteringInput(Trajectory)
#' shinyInput <- shinyInput(clustInput)
#' write.csv(shinyInput,file='/path/to/folder/shiny.csv')
#' 
#'
setMethod(f = 'getShinyInput',
		signature = signature('CombinedDataSet'),
		definition = function(CombinedDataSet){
			UMAP1 <- CombinedDataSet@seurat@reductions$umap@cell.embeddings[,1]
			UMAP2 <- CombinedDataSet@seurat@reductions$umap@cell.embeddings[,2]
			cdr3 <- CombinedDataSet@seurat$cdr3
			cell.type <- CombinedDataSet@seurat$clusters
			G.df <- data.frame(UMAP1,UMAP2,cell.type,cdr3)
			return(G.df)
		}
	)