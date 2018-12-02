# clustView
A shiny app for browsing and annotating scRNAseq clusters.

**This package is still experimental!**

<br/><br/>

### Next on the to-do list:
- generalization to `SingleCellExperiment` objects
- functionality to merge clusters
- importing cell type predictions based on references

<br/><br/>

## Requirements

Required R packages for the app:
```{r}
BiocManager::install(c( "Seurat", "shiny","shinydashboard","shinycssloaders",
   "DT", "data.table", "ggplot2", "cowplot", "AnnotationDbi", "GO.db" ) )
```

Package `clustree` is also required:
```{r}
devtools::install_github("lazappi/clustree", dependencies = TRUE)
```

In addition, preparing the cluster annotation will require the appropriate `org.Xx.eg.db` for your species, e.g. `org.Hs.eg.db`.

I tried to make it backward-compatible with Seurat v2, but right now it's only tested with V3. To install Seurat V3:
```{r}
devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')
```

<br/><br/>


## Usage

### Preparing cluster annotation:

Assuming a seurat object named `se`, first get all markers for each cluster/resolution. For example, in v3 (assuming the default prefix after integration):

```{r}
library(clustView)

# get all computed resolutions:
cn <- names(seurat@meta.data)
resolutions <- as.numeric(gsub("integrated_snn_res.","",cn[grep("^integrated_snn_res",cn)],fixed=T))

# get markers for each cluster/resolution:
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

se <- prepSeuratForClustView( se,
                              markerslist=markers,
                              species="Hs",
                              ontologies=c("BP", "CC") )
```

### Launching the app:

Assuming a seurat object named `se`:

```{r}
clustView(se)
```
