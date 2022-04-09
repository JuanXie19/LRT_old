# accesor methods
#' @describeIn CombinedDataSet returns the matrix representing the reduced
#'  dimensional dataset,i.e., UMAP coordinates
#' @param combined a \code{CombinedDataSet} object
#' @import Seurat
#' @export
setMethod(
	f = 'reducedDim',
	signature = 'CombinedDataSet',
	definition = function(CombinedDataSet) CombinedDataSet@seurat@reductions$umap@cell.embeddings

)

#' @rdname seuratClusterLabels
#' @import Seurat
#' @export
setMethod(
	f ='seuratClusterLabels',
	signature = 'CombinedDataSet',
	definition = function(CombinedDataSet) CombinedDataSet@seurat$clusters
)

#' @rdname getTCRs
#' @import Seurat
#' @export
setMethod(
	f = 'getTCRs',
	signature = 'CombinedDataSet',
	definition = function(CombinedDataSet) CombinedDataSet@seurat$cdr3
)

#' @rdname lrtLineages
#' @import slingshot
#' @import TrajectoryUtils
#' @export
setMethod(
	f = 'lrtLineages',
	signature = 'TrajectoryDataSet',
	definition = function(TrajectoryDataSet) TrajectoryDataSet@lineages
)

#' @rdname lrtParams
#' @import slingshot
#' @import TrajectoryUtils
#' @export
setMethod(
	f = 'lrtParams',
	signature = 'TrajectoryDataSet',
	definition = function(TrajectoryDataSet) TrajectoryDataSet@trajectoryParams
)