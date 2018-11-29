#' getAllGOlists
#'
#' Returned as list of genes for each GO term
#'
#' @param species Two-letter species code (corresponding to the 'org.Xx.eg' package)
#' @param translate.cat.ids Logical; whether to translate GO category IDs to terms
#' @param gene.id.type The type of gene identifier to return; 
#' either 'symbol' (default) or 'ensembl'.
#' @param genes The subset of genes to use.
#' @param minSize The minimum set size (default 10).
#' @param maxSize The maximum set size (default 500).
#' @param returnGeneSets Logical; whether to return an object of class 
#' `GeneSetCollection` rather than a normal list (default FALSE).
#' @param categories The GO categories to use (default all)
#' @param ontologies The GO tree to use; any combination of "BP", "CC", and "MF" 
#' (default all)
#'
#' @return A named list of genes, or a `GeneSetCollection` if `returnGeneSets=T`
#' @export
getAllGOlists <- function( species = "Mm",
                           translate.cat.ids = TRUE,
                           gene.id.type = "symbol",
                           genes = NULL,
                           minSize = 10,
                           maxSize = 500,
                           returnGeneSets = FALSE,
                           categories = NULL,
                           ontologies = c("BP", "CC", "MF")){
  library(AnnotationDbi)
  library(GO.db)
  gene.id.type <- match.arg(tolower(gene.id.type), c("symbol","ensembl"))
  db <- paste0('org.',species,'.eg')
  library(package=paste0(db,'.db'), character.only = T)
  ontologies <- match.arg(ontologies, c("BP","MF","CC"), several.ok=T)
  allcats <- unique(unlist(lapply(paste0("GO",ontologies,"CHILDREN"), FUN=function(x){
    unique(unlist(as.list(get(x))))
  })))
  allcats <- unique(allcats[!is.na(allcats)])
  if(!is.null(categories)) allcats <- intersect(categories, allcats)
  x <- mget(allcats, get(paste0(db,"GO2ALLEGS")), ifnotfound=NA)
  x <- x[which(sapply(x,FUN=function(x){ !all(is.na(x)) }))]
  db2id <- paste0(db, switch(gene.id.type,
                             symbol="SYMBOL",
                             ensembl="ENSEMBL",
                             stop("unknown gene id type")))
  x <- lapply(x, db=db2id, FUN=function(x,db){ unique(as.character(unlist(mget(as.character(x), get(db))))) })
  if(!is.null(genes)){
    x <- lapply(x, y=genes, FUN=intersect)
  }
  x <- x[which(sapply(x, min=minSize, max=maxSize, FUN=function(x, min, max){ length(x) > min & length(x) < max }))]
  if(translate.cat.ids) names(x) <- as.character(Term(names(x)))
  if(returnGeneSets){
    library(GSEABase)
    x <- GeneSetCollection(lapply(1:length(x), FUN=function(i){ 
      GeneSet(x[[i]], setName=names(x)[i])
    }))
  }
  x
}



#' ORA
#'
#' A simple hypergeometric over-enrichment analysis
#'
#' @param set1 Character vector; the set to test
#' @param sets Named list containing the sets for which to test enrichment
#' @param universe A character vector containing the universe/background, 
#' or a numeric or length 1 containing the size of the universe.
#' @param topN The number of top enrichments to return (default 10)
#'
#' @return A data.frame
#' @export
ORA <- function(set1, sets, universe, topN=10){
  set1 <- unique(as.character(set1))
  sets <- lapply(sets, FUN=as.character)
  sets <- lapply(sets, FUN=unique)
  if (class(universe) == "character"){
    set1 <- intersect(set1, universe)
    sets <- lapply(sets, y=universe, FUN=intersect)
    universe <- length(unique(universe))
  }
  res <- t(sapply(sets, set1=set1, u=universe, FUN=function(set2, set1, u){
    ov <- sum(set1 %in% set2)
    p <- phyper(max(0, ov - 1), length(set1), u - length(set1), 
                length(set2), lower.tail=F)
    expected <- length(set1)*length(set2)/u
    c(n=ov, enrichment=round(ov/expected,2), pval=p)
  }))
  res[order(as.numeric(res[,"pval"]))[1:min(topN,nrow(res))],]
}


.getEnsemblFromConcat <- function(genes){
  sapply(as.character(genes),FUN=function(x){ strsplit(x,".",fixed=T)[[1]][[1]] })
}
.annotateMarkers <- function(markers, go=NULL, topN=10, ...){
  markers <- split(.getEnsemblFromConcat(markers$gene), markers$cluster)
  if(is.null(go)) go <- getAllGOlists(..., genes = unique(unlist(markers)))
  lapply(markers, sets=go, universe=unlist(markers), topN=topN, FUN=ORA)
}


#' annotateAll
#'
#' @param markerslist A list of markers table (one table per resolution), 
#' as produced by `Seurat::FindAllMarkers`
#' @param go An optional named list of genesets. If NULL, will be fetched 
#' using `getAllGOlists`.
#' @param topN The top number of terms to report (default 10)
#' @param ... Passed to `getAllGOlists` if `go=NULL`
#'
#' @return A nested list of data.frames.
#' @export
annotateAll <- function(markerslist, go=NULL, topN=10, ...){
  gall <- as.character(unlist(lapply(markerslist, FUN=function(x) x$gene)))
  gall <- unique(.getEnsemblFromConcat(gall))
  if(is.null(go)) go <- getAllGOlists(..., genes=gall)
  lapply(markerslist, go=go, topN=topN, FUN=.annotateMarkers)
}



#' prepSeuratForClustView
#'
#' @param seurat An object of class `Seurat`
#' @param markerslist A nested list of markers' tables (as produced by 
#' `Seurat::FindAllMarkers()`) for each method and resolution.
#' @param res_prefixes The list of clustering prefixes to use from 
#' `seurat@meta.data` (default all detected)
#' @param species Two-letter species code (corresponding to the 
#' 'org.Xx.eg' package), default 'Mm'.
#' @param ontologies The GO tree to use; any combination of 
#' "BP", "CC", and "MF" (default `c("BP", "CC")`)
#'
#' @return The `seurat` object, with cluster annotation in the slot
#' `seurat@misc$clusterAnnotation`
#' @export
prepSeuratForClustView <- function( seurat,
                                    markerslist=NULL,
                                    res_prefixes=NULL,
                                    species="Mm",
                                    ontologies=c("BP", "CC")
){
  # we extract the different clusterings from the seurat object
  cn <- colnames(seurat@meta.data)
  cn <- cn[grep("res\\.[0-9\\.]+$",cn)]
  cn2 <- t(sapply(cn, FUN=function(x){
    x <- strsplit(x,"res.",fixed=T)[[1]]
    c(paste(c(x[1:(length(x)-1)],""),collapse="res."),x[length(x)])
  }))
  if(is.null(res_prefixes)){
    res_prefixes <- unique(cn2[,1])
    message( paste("Using prefixes:",
                   paste(res_prefixes, collapse=", ")
    )
    )
  }else{
    if(!all(res_prefixes %in% cn2[,1])){
      stop(paste("The following prefix(es) could not be found:", 
                 paste(setdiff(res_prefixes, cn2[,1]),collapse=", ")))
    }
    cn2 <- cn2[which(cn2[,1] %in% res_prefixes),,drop=F]
  }
  # overview of the different clustering and resolutions
  print(aggregate(cn2[,2], by=list(prefix=cn2[,1]), FUN=length))
  
  if(length(res_prefixes)==1){
    if(!(res_prefixes %in% names(markerslist))){
      m <- list()
      m[[res_prefixes]] <- markerslist
      markerslist <- m
    }
  }
  
  # go enrichment analysis of the markers
  allg <- .getEnsemblFromConcat(as.character(unique(unlist(markerslist))))
  go <- list()
  for(o in ontologies){
    go[[o]] <- getAllGOlists(species=species, gene.id.type='ensembl', ontologies=o, genes=allg)
  }
  gores <- lapply(markerslist, go=go, FUN=function(m, go){
    lapply(go, m=m, FUN=function(go, m){
      annotateAll(m, go=go)
    })
  })
  
  seurat@misc$clusterAnnotation <- list()
  for(i in names(markerslist)){
    seurat@misc$clusterAnnotation[[i]] <- list(markers=markerslist[[i]], go=gores[[i]])
  }
  
  seurat
}
