
#' 
#' generate a list from CombinedDataSet object
#' @param Combined  a \code{CombinedDataSet}object.
#' @param group.by The column header to group the new list. NULL will return clusters.
#' @export
list.input.return <- function(Combined,group.by) {
  if (inherits(x=Combined, what ="CombinedDataSet") | inherits(x=Combined, what ="SummarizedExperiment")) {
	if(is.null(group.by)){
      group.by <- "clusters"
    }

    df <- expression2List(Combined,group.by)
  } 
  return(df)
}


#' Allows users to take the meta data in Seurat and place it into a list 
#' that will work with all the functions
#'
#' Allows users to perform more fundamental measures of clonotype analysis 
#' using the meta data from the Seurat object. For Seurat objects the 
#' active identity is automatically added as "cluster". Remaining grouping 
#' parameters or Seurat objects must appear in the meta data.
#'
#' @examples
#' Combined <- getCombinedDataSet(TCR, Mice.sub2)
#' newList <- expression2List(Combined, group.by = "seurat_clusters")
#' 
#' @param Combined  a \code{CombinedDataSet}object.
#' @param group.by The column header to group the new list. NULL will return clusters.
#' @importFrom stringr str_sort
#' @export
#' @return list derived from the meta data of seurat object with 
#' elements divided by the group.by parameter
expression2List <- function(Combined,group.by) {
    if (!inherits(x=Combined, what ="CombinedDataSet")) {
            stop("Use a CombinedDataSet object to convert into a list")
    }
	cloneCall <- "cdr3"

    meta <- grabMeta(Combined,group.by)
    if(is.null(group.by)){
      group.by <- "clusters"
    }
    unique <- stringr::str_sort(as.character(unique(meta[,group.by])), numeric = TRUE)
    df <- NULL
    for (i in seq_along(unique)) {
        subset <- subset(meta, meta[,group.by] == unique[i])
        #subset <- subset(subset, !is.na(cloneType))
        df[[i]] <- subset
    }
    names(df) <- unique
    return(df)
}


#' This is to grab the meta data from a CombinedDataSet object
#' @param Combined  a \code{CombinedDataSet}object.
#' @param group.by The column header to group the new list. NULL will return clusters.
#' @importFrom SingleCellExperiment colData 
#' @import dplyr
#' @export
grabMeta <- function(Combined,group.by) {
    if (inherits(x=Combined, what ="CombinedDataSet")) {
	
	cloneCall <- "cdr3"

	
    meta <- data.frame(Combined@seurat[[]], slot(Combined@seurat, "active.ident"))
    if ("clusters" %in% colnames(meta)) {
          colnames(meta)[length(meta)] <- "cluster.active.ident"
    } else {
        colnames(meta)[length(meta)] <- "clusters"
     }
	
    return(meta)
	}
}

quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}


#' calculate the Morisita index for Overlap Analysis
morisitaIndex <- function(df, length, coef_matrix) {

	cloneCall <- 'cdr3'
    for (i in seq_along(length)){
        df.i <- df[[i]]
        df.i <- data.frame(table(df.i[,cloneCall]))
        colnames(df.i) <- c(cloneCall, 'Count')
        df.i[,2] <- as.numeric(df.i[,2])
        for (j in seq_along(length)){
            if (i >= j){ next }
            else { df.j <- df[[j]]
            df.j <- data.frame(table(df.j[,cloneCall]))
            colnames(df.j) <- c(cloneCall, 'Count')
            df.j[,2] <- as.numeric(df.j[,2])
            merged <- merge(df.i, df.j, by = cloneCall, all = TRUE)
            merged[is.na(merged)] <- 0
            X <- sum(merged[,2])
            Y <- sum(merged[,3])
            sum.df.i <- sum(df.i[,2]^2)
            sum.df.j <- sum(df.j[,2]^2)
            
            num <- 2 * sum(merged[, 2] * merged[, 3])
            den <- ((sum.df.i / (X^2) + sum.df.j / (Y^2)) * X * Y)
                
            coef.i.j <- num/den
            coef_matrix[i,j] <- coef.i.j
            }
        }
    }
    return(coef_matrix)
}


#' Calculate the Jaccard Similarity Index for Overlap Analysis
jaccardIndex <- function(df, length, coef_matrix) {
  	cloneCall <- 'cdr3'
 
  for (i in seq_along(length)){
    df.i <- df[[i]]
    df.i <- df.i[,c("barcode",cloneCall)]
    df.i_unique <- df.i[!duplicated(df.i[,cloneCall]),]
    for (j in seq_along(length)){
      if (i >= j){ next }
      else { 
        df.j <- df[[j]]
        df.j <- df.j[,c("barcode",cloneCall)]
        df.j_unique <- df.j[!duplicated(df.j[,cloneCall]),]
        overlap <- length(intersect(df.i_unique[,cloneCall], 
                                    df.j_unique[,cloneCall]))
        coef_matrix[i,j] <- 
          overlap/(sum(length(df.i_unique[,cloneCall]), 
                                  length(df.j_unique[,cloneCall]))-overlap)
      } 
    }
  }
  return(coef_matrix)
}


rawIndex <- function(df, length, coef_matrix) {
 	cloneCall <- 'cdr3'
 for (i in seq_along(length)){
    df.i <- df[[i]]
    df.i <- df.i[,c("barcode",cloneCall)]
    df.i_unique <- df.i[!duplicated(df.i[,cloneCall]),]
    for (j in seq_along(length)){
      if (i >= j){ next }
      else { 
        df.j <- df[[j]]
        df.j <- df.j[,c("barcode",cloneCall)]
        df.j_unique <- df.j[!duplicated(df.j[,cloneCall]),]
        overlap <- length(intersect(df.i_unique[,cloneCall], 
                                    df.j_unique[,cloneCall]))
        coef_matrix[i,j] <- overlap
      } 
    }
  }
  return(coef_matrix)
}




#' Calculate the Overlap Coefficient for Overlap Analysis
#' @author Nick Bormann, Nick Borcherding
overlapIndex <- function(df, length, coef_matrix) {
    	cloneCall <- 'cdr3'

	for (i in seq_along(length)){
        df.i <- df[[i]]
        df.i <- df.i[,c("barcode",cloneCall)]
        df.i_unique <- df.i[!duplicated(df.i[,cloneCall]),]
        for (j in seq_along(length)){
            if (i >= j){ next }
            else { df.j <- df[[j]]
            df.j <- df.j[,c("barcode",cloneCall)]
            df.j_unique <- df.j[!duplicated(df.j[,cloneCall]),]
            overlap <- length(intersect(df.i_unique[,cloneCall], 
                                        df.j_unique[,cloneCall]))
            coef_matrix[i,j] <- 
                overlap/min(length(df.i_unique[,cloneCall]), 
                length(df.j_unique[,cloneCall])) } } }
    return(coef_matrix)
}

theCall <- function(x) {
    if (x %in% c("CTnt", "CTgene", "CTaa", "CTstrict")) {
      x <- x
    }else if (x == "gene" | x == "genes") {
        x <- "CTgene"
    } else if(x == "nt" | x == "nucleotide") {
        x <- "CTnt"
    } else if (x == "aa" | x == "amino") {
        x <- "CTaa"
    } else if (x == "gene+nt" | x == "strict") {
        x <- "CTstrict"
    }
    return(x)
}

checkBlanks <- function(df) {
	cloneCall <- "cdr3"

    for (i in seq_along(df)) {
        if (length(df[[i]][,cloneCall]) == length(which(is.na(df[[i]][,cloneCall]))) | 
            length(which(!is.na(df[[i]][,cloneCall]))) == 0 | 
            nrow(df[[i]]) == 0) {
            df[[i]] <- NULL
        } else {
            next()
        }
    }
    return(df)
}

short.check <- function(df) {
	cloneCall <- "cdr3"

  min <- c()
  for (x in seq_along(df)) {
    min.tmp <- length(which(!is.na(unique(df[[x]][,cloneCall]))))
    min <- c(min.tmp, min)
  }
  min <- min(min)
  return(min)
}

#' Calculating diversity using Vegan R package
#' @importFrom vegan diversity estimateR
diversityCall <- function(data) {
    w <- diversity(data[,"Freq"], index = "shannon")
    x <- diversity(data[,"Freq"], index = "invsimpson")
    y <- estimateR(data[,"Freq"])[2] #Chao
    z <- estimateR(data[,"Freq"])[4] #ACE
    z2 <- diversity(data[,"Freq"], index = "shannon")/log(length(data[,"Freq"]))
    out <- c(w,x,y,z, z2)
    return(out)
}

checkList <- function(df) {
    df <- if(is(df)[1] != "list") list(df) else df
    return(df)
}

