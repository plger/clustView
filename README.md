# clustView
shiny viewer for browsing and annotating clusters

# Requirements

Required R packages for the app:
```{r}
c( "Seurat", "shiny","shinydashboard","shinycssloaders",
   "DT", "data.table", "ggplot2", "cowplot" )
```

Package `clustree` is also required if `enableClustree=TRUE`.

I tried to make it backward-compatible with Seurat v2, but right now it's only tested with V3.

Required packages for preparing the cluster annotation: `AnnotationDbi`, `GO.db`, and the appropriate `org.Xx.eg.db` for your species.

# Usage

## Prepare cluster annotation:

Assuming a seurat object named `se`, first get all markers for each cluster/resolution. For example, in v3 (assuming the default prefix after integration):

```{r}
# get all computed resolutions:
cn <- names(seurat@meta.data)
resolutions <- as.numeric(gsub("integrated_snn_res.","",cn[grep("^integrated_snn_res",cn)],fixed=T))

# get markers for each clustering:
markers <- list()
for(r in resolutions){
  markers[[as.character(r)]] <- FindAllMarkers( 
	  	object = se,
	  	only.pos = FALSE, 
	    min.pct = 0.25, 
	    resolution = r,
	    logfc.threshold = 0.5,
	    test.use = "wilcox",
	    max.cells.per.ident = 300 )
}

source("path/to/plg.GO.R")
se <- prepSeuratForClustView( se,
                              markerslist=markers,
                              species="Hs",
                              ontologies=c("BP", "CC") )
```

## Launching the app:

Assuming a seurat object named `se`:

```{r}
source("path/to/clustView.ui.R")
source("path/to/clustView.server.R")
clustView(se)
```

You can use `clustView(se, enableClustree=F)` to disable the clustering tree (saves loading time and does not require the `clustree` package.) 