#' @title Initialize an object of class \code{CombinedDataSet}
#' @name getCombinedDataSet
#' @docType methods
#'
#' @description Construct a \code{CombinedDataSet} object with Seurat clustering results and clonotype information . Additional helper methods 
#'  for manipulating \code{CombinedDatSet} objects are also described below.
#'
#' @param scTCRseq The preprocessed scTCR-seq data. See details
#' @param seurat An integrated Seurat object. The object must contain clustering results and saved as 'clusters', and must have cell barcodes saved as 'cell.id'
#'
#' @details The input scTCR-seq needs to be a dataframe, and alpha chain and beta chain need to be paired. If you have multiple sets of scTCR-seq from more than one sample, the data needs to be combined into one.
#' @details The seurat object must contain clustering results and saved as 'clusters', contain cell barcodes saved as 'cell.id' and also need to have UMAP coordinates for cells.
#'
#' @import Seurat
#' @import dplyr
#'
#' @export
#'

setGeneric(
	name = 'getCombinedDataSet',
	signature = c('TCR','seurat'),
	def = function(TCR,seurat,...){
		standardGeneric('getCombinedDataSet')
		}
)


#' @title Infer clonotype lineage structure from \code{CombinedDatSet}
#' @name getClonotypeLineages
#'
#' @param CombinedDataSet The \code{CombinedDataSet} object obtained from \code{getCombinedDataSet}
#' @param ... Additional arguments to specify how lineages are constructed from
#'   clusters.
#' @export
setGeneric(
    name = "getClonotypeLineages",
	signature = 'CombinedDataSet',
    def = function(CombinedDataSet,...) {
        standardGeneric("getClonotypeLineages")
    }
)

#' @title prepare the input for \code{clustLineages}
#' @name getClusteringInput
#'
#' @param TrajectoryDataSet The \code{TrajectoryDataSet} object obtained from \code{getClonotypeLineages}
#' 
#' @export

setGeneric(
    name = "getClusteringInput",
	signature = 'TrajectoryDataSet',
    def = function(TrajectoryDataSet,...) {
        standardGeneric("getClusteringInput")
    }
)

